"""grid_refit.py - centre a plot grid on the rows that were actually drilled.

Generic across trials. Nothing about plot size, row count, row spacing or CRS is
hard-coded: everything is measured from the orthomosaic and the shapefile.

WHY
    A plot grid digitised by hand or generated from a design sits a few
    centimetres off the real sowing. Downstream that shows up as crops whose
    rows are not centred, plots that appear to lose a row, and - if row snapping
    is enabled - crops that follow the rows instead of the plot and so no longer
    correspond to the plot they are named after. Correcting the grid once fixes
    all of it, and makes row snapping unnecessary.

THE ONE THING THIS TOOL REFUSES TO DO
    If a trial is drilled CONTINUOUSLY - row pitch equal to plot width divided by
    the row count - then a plot boundary leaves no signature in the imagery at
    all, and "centring on the rows" is meaningless: every position is equally
    valid modulo one row. The tool measures this and REFUSES to refit, rather
    than moving polygons on the strength of detection noise. Muresk
    NUE/RootPhenotyping is such a trial (216 mm measured vs 1500/7 = 214.3);
    OZ Barley York is not (212 mm within a plot vs 264 mm to the neighbour).

USAGE
    python grid_refit.py --shp GRID.shp --ortho ORTHO.tif [--rows N]
                         [--out REFIT.shp] [--report-only]

OUTPUT
    <out>.shp                 the corrected grid, same attributes and CRS
    <out>_refit_audit.csv     per plot: rows found, gaps, shift applied, flags
    a diagnosis printed to stdout, including the feasibility test
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys

import numpy as np

# set once in main(): pitch hint (mm) and ortho pixel size (m).
# The sampler uses them only to choose a decimation factor.
_PITCH_FOR_DECIM = [200.0]
_PX_M = [0.001]


# ----------------------------------------------------------------------------
# signal helpers
# ----------------------------------------------------------------------------
def smooth(v, k):
    if k < 3:
        return v
    return np.convolve(v, np.ones(k, np.float32) / k, mode="same")


def find_peaks(p, min_dist, prom):
    """local maxima, greedy by height, with a minimum separation and prominence"""
    n = len(p)
    cand = [i for i in range(1, n - 1) if p[i] >= p[i - 1] and p[i] > p[i + 1]]
    cand.sort(key=lambda i: -p[i])
    keep = []
    for i in cand:
        if any(abs(i - j) < min_dist for j in keep):
            continue
        lo = max(0, i - int(min_dist))
        hi = min(n, i + int(min_dist) + 1)
        left = p[lo:i].min() if i > lo else 0.0
        right = p[i + 1:hi].min() if hi > i + 1 else 0.0
        if p[i] - max(left, right) < prom:
            continue
        keep.append(i)
    return sorted(keep)


# ----------------------------------------------------------------------------
# geometry
# ----------------------------------------------------------------------------
def plot_axes(ring_prj):
    """Return (A, across_unit, along_unit, width_m, length_m) in projected units.

    A is the corner from which the SHORT edge leaves; rows run along the long
    axis, so 'across' is the direction a row-to-row step takes.
    """
    P = np.asarray(ring_prj, dtype=float)[:4]
    d = [float(np.hypot(*(P[(i + 1) % 4] - P[i]))) for i in range(4)]
    a0 = 0 if (d[0] + d[2]) <= (d[1] + d[3]) else 1
    A = P[a0 % 4]
    B = P[(a0 + 1) % 4]
    D = P[(a0 + 3) % 4]
    w = float(np.hypot(*(B - A)))
    l = float(np.hypot(*(D - A)))
    u = (B - A) / (w or 1.0)
    v = (D - A) / (l or 1.0)
    return A, u, v, w, l


def sample_profile(src, ring_geo, ring_prj, extend, ns, nl, nodata_alpha=True):
    """Mean ExG across the rows, sampled in the plot's own axes and extended
    beyond both across-row edges. Returns (profile, mm_per_step, width_m)."""
    import rasterio

    A, u, v, w, l = plot_axes(ring_prj)
    G = np.asarray(ring_geo, dtype=float)[:4]
    d = [float(np.hypot(*(np.asarray(ring_prj)[:4][(i + 1) % 4]
                          - np.asarray(ring_prj)[:4][i]))) for i in range(4)]
    a0 = 0 if (d[0] + d[2]) <= (d[1] + d[3]) else 1
    Ag, Bg, Cg, Dg = (G[(a0 + k) % 4] for k in range(4))

    s = np.linspace(-extend, 1.0 + extend, ns)
    t = np.linspace(0.06, 0.94, nl)
    S, T = np.meshgrid(s, t, indexing="ij")
    x = (Ag[0] + S * (Bg[0] - Ag[0]) + T * (Dg[0] - Ag[0])
         + S * T * (Cg[0] - Bg[0] - Dg[0] + Ag[0]))
    y = (Ag[1] + S * (Bg[1] - Ag[1]) + T * (Dg[1] - Ag[1])
         + S * T * (Cg[1] - Bg[1] - Dg[1] + Ag[1]))
    rows, cols = rasterio.transform.rowcol(src.transform, x.ravel(), y.ravel())
    rows = np.clip(np.asarray(rows), 0, src.height - 1)
    cols = np.clip(np.asarray(cols), 0, src.width - 1)
    r0, r1 = int(rows.min()), int(rows.max()) + 1
    c0, c1 = int(cols.min()), int(cols.max()) + 1
    win = rasterio.windows.Window(c0, r0, c1 - c0, r1 - r0)
    # Read through an overview. A native-resolution ortho can be tens of
    # gigapixels (the OZ Barley one is 300k x 310k at 0.37 mm), so a full-res
    # window would be ~230 MB per plot and hundreds of GB over a trial. Row
    # detection only needs the pitch resolved to a few tens of pixels, so
    # decimate to about 30 px per pitch and let GDAL serve it from an overview.
    pitch_px = max(1.0, _PITCH_FOR_DECIM[0] / (_PX_M[0] * 1000.0))
    f = max(1, int(pitch_px / 30.0))
    oh = max(1, int(math.ceil((r1 - r0) / f)))
    ow = max(1, int(math.ceil((c1 - c0) / f)))
    a = src.read(window=win, out_shape=(src.count, oh, ow))
    rr = np.clip(((rows - r0) / f).astype(int), 0, oh - 1)
    cc = np.clip(((cols - c0) / f).astype(int), 0, ow - 1)
    R = a[0][rr, cc].reshape(S.shape).astype(np.float32)
    Gc = a[1][rr, cc].reshape(S.shape).astype(np.float32)
    B = a[2][rr, cc].reshape(S.shape).astype(np.float32)
    if nodata_alpha and a.shape[0] >= 4:
        valid = a[3][rr, cc].reshape(S.shape) > 0
    else:
        valid = np.maximum.reduce([R, Gc, B]) > 4
    if valid.mean() < 0.5:
        return None, None, w
    exg = 2 * Gc - R - B
    wsum = valid.sum(axis=1)
    prof = np.where(wsum > 0, (np.where(valid, exg, 0)).sum(axis=1)
                    / np.maximum(wsum, 1e-6), 0.0)
    mm_per_step = (w * 1000.0) * (1.0 + 2 * extend) / ns
    return prof, mm_per_step, w


# ----------------------------------------------------------------------------
# per-plot measurement
# ----------------------------------------------------------------------------
def measure(src, pid, ring_geo, ring_prj, pitch_hint, extend, ns, nl):
    out = {"plot_id": pid}
    prof, mmstep, w = sample_profile(src, ring_geo, ring_prj, extend, ns, nl)
    out["plot_w_mm"] = round(w * 1000.0, 1)
    if prof is None:
        out["flag"] = "nodata"
        return out
    md = max(3.0, pitch_hint / mmstep * 0.55)
    prof = smooth(prof, max(3, int(pitch_hint / mmstep * 0.30) | 1))
    p = (prof - prof.min()) / (np.ptp(prof) + 1e-9)
    pk = find_peaks(p, md, 0.09)
    if len(pk) < 3:
        out["flag"] = "rows_not_found"
        out["n_all"] = len(pk)
        return out
    # position of each detected row, mm from the plot's low across-row edge
    s = np.linspace(-extend, 1.0 + extend, ns)
    pos = np.array([s[i] * w * 1000.0 for i in pk])
    out["n_all"] = len(pos)
    out["rows_mm"] = ";".join(f"{v:.0f}" for v in pos)
    g = np.diff(pos)
    g = g[(g > pitch_hint * 0.55) & (g < pitch_hint * 1.8)]
    out["gap_med_mm"] = round(float(np.median(g)), 1) if g.size else ""
    inside = pos[(pos >= 0) & (pos <= w * 1000.0)]
    out["n_in"] = len(inside)
    out["_pos"] = pos
    return out


def pick_block(pos, n, pitch, centre_mm):
    """Best contiguous run of n rows: regular gaps first, then nearest the centre."""
    best, score = None, None
    for k in range(0, len(pos) - n + 1):
        wnd = pos[k:k + n]
        gaps = np.diff(wnd)
        reg = float(np.mean(np.abs(gaps - pitch)))
        cen = abs((wnd[0] + wnd[-1]) / 2.0 - centre_mm)
        sc = 3.0 * reg + cen
        if score is None or sc < score:
            best, score = k, sc
    if best is None:
        return None, None
    wnd = pos[best:best + n]
    return wnd, float(np.mean(np.abs(np.diff(wnd) - pitch)))


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shp", required=True, help="plot grid shapefile")
    ap.add_argument("--ortho", required=True, help="orthomosaic covering the trial")
    ap.add_argument("--rows", type=int, default=0,
                    help="rows per plot; 0 = infer from the data")
    ap.add_argument("--out", default="", help="output shapefile (default <shp>_refit.shp)")
    ap.add_argument("--extend", type=float, default=0.45,
                    help="fraction of plot width sampled beyond each across-row edge")
    ap.add_argument("--ns", type=int, default=320, help="samples across the rows")
    ap.add_argument("--nl", type=int, default=240, help="samples along the rows")
    ap.add_argument("--max-shift-frac", type=float, default=0.5,
                    help="reject a shift larger than this fraction of the row pitch")
    ap.add_argument("--report-only", action="store_true",
                    help="diagnose and measure, write no shapefile")
    ap.add_argument("--force", action="store_true",
                    help="refit even if the continuity test says it is meaningless")
    a = ap.parse_args(argv)

    import rasterio
    import shapefile
    from pyproj import CRS, Transformer

    out_shp = a.out or a.shp.replace(".shp", "_refit.shp")
    audit = os.path.splitext(out_shp)[0] + "_refit_audit.csv"

    scrs = CRS.from_wkt(open(a.shp.replace(".shp", ".prj")).read())
    with rasterio.open(a.ortho) as src:
        ocrs = CRS.from_user_input(src.crs)
        px = abs(src.transform.a)
        print(f"grid  : {scrs.name} (epsg {scrs.to_epsg()})")
        print(f"ortho : {ocrs.name} (epsg {ocrs.to_epsg()})  {src.width}x{src.height}px  "
              f"pixel {px:g} {'deg' if ocrs.is_geographic else 'units'}")
        to_o = Transformer.from_crs(scrs, ocrs, always_xy=True)

        r = shapefile.Reader(a.shp)
        recs = r.records()
        shps = r.shapes()
        fields = r.fields[1:]
        n_plots = len(recs)
        idf = next((f[0] for f in fields if f[0].lower() in
                    ("plot_id", "plotid", "id", "plot")), fields[0][0])
        print(f"plots : {n_plots}   id field {idf!r}")

        rings_prj, rings_geo, ids = [], [], []
        for i in range(n_plots):
            pts = np.array(shps[i].points[:5], dtype=float)
            if len(pts) < 4:
                pts = np.array((list(shps[i].points) + list(shps[i].points))[:5], float)
            gx, gy = to_o.transform(pts[:, 0], pts[:, 1])
            rings_prj.append(pts)
            rings_geo.append(np.column_stack([gx, gy]))
            ids.append(recs[i][idf])

        # plot geometry, and a first pitch guess from width / rows
        W = np.array([plot_axes(rp)[3] for rp in rings_prj])
        L = np.array([plot_axes(rp)[4] for rp in rings_prj])
        print(f"plot  : {np.median(W):.3f} x {np.median(L):.3f} m "
              f"(width sd {W.std()*1000:.0f} mm, length sd {L.std()*1000:.0f} mm)")

        rows_guess = a.rows if a.rows else 0
        pitch_hint = (np.median(W) * 1000.0 / rows_guess) if rows_guess else 200.0
        _PITCH_FOR_DECIM[0] = pitch_hint
        _PX_M[0] = px if not ocrs.is_geographic else px * 111320.0

        # pass 1: measure with the hint
        print(f"\npass 1: measuring rows (pitch hint {pitch_hint:.0f} mm) ...")
        res = []
        for i in range(n_plots):
            res.append(measure(src, ids[i], rings_geo[i], rings_prj[i],
                               pitch_hint, a.extend, a.ns, a.nl))
            if (i + 1) % 100 == 0:
                print(f"   {i+1}/{n_plots}")

        gaps = np.array([r_["gap_med_mm"] for r_ in res
                         if isinstance(r_.get("gap_med_mm"), float)])
        if gaps.size < 10:
            print("could not measure a row pitch - aborting")
            return 2
        pitch = float(np.median(gaps))
        _PITCH_FOR_DECIM[0] = pitch
        # rows per plot: mode of the count inside a plot
        nin = np.array([r_.get("n_in", 0) for r_ in res if r_.get("n_in")])
        rows_n = a.rows or int(np.bincount(nin).argmax())
        print(f"\nmeasured row pitch      : {pitch:.1f} mm  (sd {gaps.std():.1f})")
        print(f"rows per plot            : {rows_n}"
              + ("" if a.rows else "  (inferred as the commonest count inside a plot)"))

        # ---- the feasibility test -------------------------------------------
        continuous_pitch = float(np.median(W)) * 1000.0 / rows_n
        delta = pitch - continuous_pitch
        span = (rows_n - 1) * pitch
        gap_to_neighbour = float(np.median(W)) * 1000.0 - span
        print(f"\nCONTINUITY TEST")
        print(f"  measured pitch                    {pitch:8.1f} mm")
        print(f"  pitch if drilled continuously     {continuous_pitch:8.1f} mm  "
              f"(plot width / rows)")
        print(f"  difference                        {delta:+8.1f} mm")
        print(f"  {rows_n} rows span                       {span:8.1f} mm")
        print(f"  gap from outer row to neighbour   {gap_to_neighbour:8.1f} mm")
        # The discriminator is the OUTER-ROW-TO-NEIGHBOUR gap against the pitch,
        # not the pitch against width/rows. If drilling is continuous the plot
        # boundary sits half a pitch past the outer row on each side, so that gap
        # equals one pitch. A real alley makes it larger. Comparing pitch with
        # width/rows is the same statement but far less sensitive - on OZ Barley it
        # came out 11.4 mm apart on a 212 mm pitch and wrongly read as continuous,
        # while the gap ratio is an unambiguous 1.27.
        gap_ratio = gap_to_neighbour / pitch if pitch else 1.0
        print(f"  gap / pitch ratio                 {gap_ratio:8.3f}"
              f"   (1.00 = continuous, >1 = a real gap)")
        feasible = gap_ratio >= 1.12
        if feasible:
            print(f"  -> plot boundaries ARE distinguishable "
                  f"({gap_to_neighbour:.0f} mm gap vs {pitch:.0f} mm pitch): refit is meaningful")
        else:
            print(f"  -> rows are drilled CONTINUOUSLY across plots: a plot boundary has")
            print(f"     no signature in the imagery, so centring on rows is meaningless.")
            if not a.force:
                print(f"     REFUSING to refit. The existing grid is the only definition of")
                print(f"     a plot edge; use it as-is with row snapping OFF. (--force overrides.)")

        # ---- per-plot shift --------------------------------------------------
        keys = ["plot_id", "plot_w_mm", "n_all", "n_in", "gap_med_mm", "block_span_mm",
                "block_gap_err_mm", "offset_mm", "global_mm", "resid_mm", "shift_mm",
                "action", "flag", "rows_mm"]
        do_refit = feasible or a.force

        # First measure every plot, then split the correction in two:
        #   global  - the median offset, a systematic placement error of the whole
        #             grid, applied to EVERY plot including ones whose rows could
        #             not be measured;
        #   residual- what is left per plot, applied only where the detection is
        #             clean and the residual is small.
        # Without this split a systematic bias inflates every per-plot residual and
        # trips the safety guard, so the plots that most need correcting are the
        # ones that get skipped.
        pre = []
        for r_ in res:
            pos = r_.get("_pos")
            if pos is None or len(pos) < rows_n:
                continue
            centre = r_["plot_w_mm"] / 2.0
            wnd, gerr = pick_block(pos, rows_n, pitch, centre)
            if wnd is None:
                continue
            r_["_blk"] = (wnd, gerr)
            if gerr <= 0.18 * pitch:
                pre.append(float((wnd[0] + wnd[-1]) / 2.0 - centre))
        global_mm = float(np.median(pre)) if len(pre) >= 20 else 0.0
        print("")
        print(f"systematic grid offset (median over {len(pre)} clean plots):"
              f" {global_mm:+.1f} mm  -> applied to all plots")
        shifts = {}
        limit = a.max_shift_frac * pitch
        for r_ in res:
            pos = r_.pop("_pos", None)
            blk = r_.pop("_blk", None)
            r_.setdefault("flag", "")
            r_["global_mm"] = round(global_mm, 1)
            total = global_mm
            if blk is None:
                r_["action"] = ("global only - rows not measured" if do_refit
                                else "measured only (continuity test)")
                r_["shift_mm"] = round(total, 1) if do_refit else 0.0
                if do_refit and total:
                    shifts[r_["plot_id"]] = total / 1000.0
                continue
            wnd, gerr = blk
            centre = r_["plot_w_mm"] / 2.0
            off = float((wnd[0] + wnd[-1]) / 2.0 - centre)
            resid = off - global_mm
            r_["block_span_mm"] = round(float(wnd[-1] - wnd[0]), 1)
            r_["block_gap_err_mm"] = round(gerr, 1)
            r_["offset_mm"] = round(off, 1)
            r_["resid_mm"] = round(resid, 1)
            if gerr > 0.18 * pitch:
                r_["flag"] = "irregular"
                r_["action"] = "global only - irregular gaps"
            elif abs(resid) > limit:
                r_["flag"] = "large_residual"
                r_["action"] = f"global only - residual > {limit:.0f} mm"
            else:
                total = off
                r_["action"] = "global + residual"
            if not do_refit:
                r_["action"] = "measured only (continuity test)"
                r_["shift_mm"] = 0.0
                continue
            r_["shift_mm"] = round(total, 1)
            if total:
                shifts[r_["plot_id"]] = total / 1000.0

    with open(audit, "w", newline="") as fh:
        w_ = csv.DictWriter(fh, fieldnames=keys)
        w_.writeheader()
        for r_ in res:
            w_.writerow({k: r_.get(k, "") for k in keys})
    print(f"\naudit -> {audit}")

    off = np.array([r_["offset_mm"] for r_ in res if isinstance(r_.get("offset_mm"), float)])
    if off.size:
        print(f"\nBEFORE refit, block centre vs plot centre:")
        print(f"  mean {off.mean():+.1f}  sd {off.std():.1f}  |off| med {np.median(np.abs(off)):.1f} mm"
              f"  max {np.abs(off).max():.1f}")
    acts = {}
    for r_ in res:
        acts[r_.get("action", "?")] = acts.get(r_.get("action", "?"), 0) + 1
    for k in sorted(acts):
        print(f"  {k:34s} {acts[k]:4d}")

    if not do_refit or a.report_only:
        print("")
        print("no shapefile written"
              + (" (--report-only)" if a.report_only else " (continuity test refused the refit)"))
        # Exit code is the contract the GUI button switches on, so it must
        # depend on the VERDICT, not on how the tool was invoked:
        #   0 = a refit was written, or would be on a feasible trial
        #   3 = refused because the trial is drilled continuously
        # Returning 0 for --report-only regardless of feasibility conflated
        # the two and would show "Grid corrected" for a refused trial.
        return 0 if feasible else 3

    # ---- write the corrected grid ------------------------------------------
    import shapefile
    w_out = shapefile.Writer(out_shp.replace(".shp", ""), shapeType=shapefile.POLYGON)
    for f in fields:
        w_out.field(*f)
    moved = 0
    for i in range(n_plots):
        pid = ids[i]
        pts = np.array(shps[i].points, dtype=float)
        sh = shifts.get(pid)
        if sh:
            _, u, _, _, _ = plot_axes(rings_prj[i])
            pts = pts + u * sh
            moved += 1
        w_out.poly([[tuple(p) for p in pts]])
        w_out.record(*[recs[i][f[0]] for f in fields])
    w_out.close()
    with open(out_shp.replace(".shp", ".prj"), "w") as fh:
        fh.write(open(a.shp.replace(".shp", ".prj")).read())
    print(f"\nrefit grid -> {out_shp}   ({moved} of {n_plots} polygons moved)")
    sm = np.array([v * 1000.0 for v in shifts.values()])
    if sm.size:
        print(f"  shift applied: mean {sm.mean():+.1f}  sd {sm.std():.1f}  "
              f"|shift| med {np.median(np.abs(sm)):.1f}  max {np.abs(sm).max():.1f} mm")
    return 0


if __name__ == "__main__":
    sys.exit(main())
