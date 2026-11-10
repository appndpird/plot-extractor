"""crop_rows.py - crop-row lines inside the plots of a plot grid.

Two ways to place the rows, as in PhenoFuse 4.0.7 (phenoapp/core/grid_rows.py):

  count    a KNOWN number of rows per plot, evenly spaced by the seeder: spacing =
           plot width / rows (or a given spacing, rows centred), row k lies (k - 0.5)
           spacings from the plot's side edge and runs the full plot length. No image.
  detect   the rows found in the orthomosaic: each plot is resampled into its own frame
           (along x across the rows) at a few mm per pixel and row_detect.detect_rows
           (the PhenoFuse row method used by the tiles tab) finds the pitch, the tilt
           and every row.

Rows are saved as a polyline shapefile next to the grid, <grid>_rows.shp, one line per
row: Plot_ID, row_idx (1 = the row nearest the plot's first side), v_off_m (signed
offset of the row from the plot centre, across, metres), length_m, angle_deg (tilt
against the plot's long side), source (count / detect / edit) and contrast (detect).

Pure numpy + pyshp (+ rasterio for detection); no shapely / geopandas.
"""
from __future__ import annotations

import csv
import math
import os

import numpy as np

ROW_FIELDS = ("Plot_ID", "row_idx", "v_off_m", "length_m", "angle_deg", "source", "contrast")


# ---------------------------------------------------------------------------
# coordinates: the grid's CRS <-> local metres
# ---------------------------------------------------------------------------
class Metric:
    """Local metric frame for a set of rings: identity for a projected CRS, an
    equirectangular scaling around the grid centre for lon/lat degrees."""

    def __init__(self, rings):
        xs = [p[0] for r in rings for p in r]
        ys = [p[1] for r in rings for p in r]
        self.geographic = max(abs(min(xs)), abs(max(xs)), abs(min(ys)), abs(max(ys))) <= 360
        if self.geographic:
            lat = sum(ys) / len(ys)
            self.mx, self.my = 111320.0 * math.cos(math.radians(lat)), 110540.0
        else:
            self.mx = self.my = 1.0

    def to_m(self, p):
        return np.array([p[0] * self.mx, p[1] * self.my], float)

    def from_m(self, q):
        return (float(q[0]) / self.mx, float(q[1]) / self.my)


def plot_frame(ring_m, along="long"):
    """(centre, a_u, a_v, half_length, half_width) of a polygon in metres.

    The axes are those of the polygon's minimum-area rectangle. a_u runs ALONG the
    crop rows: the rectangle's long side (along='long', the usual plot), or its short
    side (along='short'). a_v = a_u turned +90 degrees (across the rows)."""
    P = np.asarray(ring_m, float)
    best = None
    for i in range(len(P)):
        e = P[(i + 1) % len(P)] - P[i]
        n = math.hypot(*e)
        if n < 1e-9:
            continue
        u = e / n
        v = np.array([-u[1], u[0]])
        a, b = P @ u, P @ v
        area = (a.max() - a.min()) * (b.max() - b.min())
        if best is None or area < best[0] - 1e-12:
            best = (area, u, v, a, b)
    _, u, v, a, b = best
    la, lb = a.max() - a.min(), b.max() - b.min()
    c = u * (a.max() + a.min()) / 2 + v * (b.max() + b.min()) / 2
    if (la >= lb) == (along == "long"):
        a_u, hl, hw = u, la / 2, lb / 2
    else:
        a_u, hl, hw = v, lb / 2, la / 2
    a_v = np.array([-a_u[1], a_u[0]])
    return c, a_u, a_v, hl, hw


def grid_frames(rings, along="long"):
    """plot_frame for every ring, with a_u pointing the same way in every plot (the
    first plot's direction), so row 1 is on the same side of every plot."""
    M = Metric(rings)
    out, ref = [], None
    for r in rings:
        c, a_u, a_v, hl, hw = plot_frame([M.to_m(p) for p in r], along)
        if ref is None:
            ref = a_u
        if float(a_u @ ref) < 0:
            a_u, a_v = -a_u, -a_v
        out.append((c, a_u, a_v, hl, hw))
    return M, out


def _line(M, c, a_u, a_v, hl, v, angle_deg=0.0):
    """World end points of a row v metres across from the centre, hl either side
    along, tilted by angle_deg."""
    t = math.tan(math.radians(angle_deg))
    p0 = c + a_v * (v - t * hl) - a_u * hl
    p1 = c + a_v * (v + t * hl) + a_u * hl
    return [M.from_m(p0), M.from_m(p1)]


# ---------------------------------------------------------------------------
# rows from a count (fixed spacing)
# ---------------------------------------------------------------------------
def rows_from_count(rings, ids, n_rows, pitch_m=None, along="long"):
    """Evenly spaced rows. n_rows: one int, or one per plot (0 = no rows there).
    pitch_m: a known spacing (rows centred) instead of width / n_rows.
    Returns (rows, per-plot pitch list)."""
    M, fr = grid_frames(rings, along)
    counts = [int(n_rows)] * len(rings) if np.isscalar(n_rows) else [int(x or 0) for x in n_rows]
    rows, pitches = [], []
    for pid, (c, a_u, a_v, hl, hw), n in zip(ids, fr, counts):
        if n <= 0:
            pitches.append(float("nan"))
            continue
        p = float(pitch_m) if pitch_m else 2.0 * hw / n
        pitches.append(p)
        for k in range(n):
            v = (k - (n - 1) / 2.0) * p
            rows.append(dict(Plot_ID=str(pid), row_idx=k + 1, v_off_m=round(v, 4),
                             length_m=round(2 * hl, 3), angle_deg=0.0, source="count",
                             contrast=float("nan"), line=_line(M, c, a_u, a_v, hl, v)))
    return rows, pitches


# ---------------------------------------------------------------------------
# rows detected in the orthomosaic
# ---------------------------------------------------------------------------
def _plot_window(ds, to_ds, M, c, a_u, a_v, hl, hw, pad, step):
    """The orthomosaic resampled into the plot's frame: 3 x n_across x n_along uint8,
    along = image x, across = image y (top row = -hw - pad)."""
    import rasterio
    n_l = max(8, int(round(2 * hl / step)))
    n_v = max(8, int(round(2 * (hw + pad) / step)))
    ls = -hl + (np.arange(n_l) + 0.5) * (2 * hl / n_l)
    vs = -hw - pad + (np.arange(n_v) + 0.5) * step
    L, V = np.meshgrid(ls, vs)                         # n_v x n_l
    X = c[0] + a_u[0] * L + a_v[0] * V
    Y = c[1] + a_u[1] * L + a_v[1] * V
    wx, wy = X.ravel() / M.mx, Y.ravel() / M.my        # grid CRS
    if to_ds is not None:
        wx, wy = to_ds(wx, wy)
    inv = ~ds.transform
    col, row = inv * (np.asarray(wx), np.asarray(wy))
    c0, c1 = int(np.floor(col.min())) - 1, int(np.ceil(col.max())) + 2
    r0, r1 = int(np.floor(row.min())) - 1, int(np.ceil(row.max())) + 2
    c0, r0 = max(c0, 0), max(r0, 0)
    c1, r1 = min(c1, ds.width), min(r1, ds.height)
    if c1 - c0 < 4 or r1 - r0 < 4:
        return None
    # read decimated so one ortho read pixel ~ one sample (GDAL uses overviews):
    # ortho pixels between two neighbouring samples
    d = math.hypot(col[1] - col[0], row[1] - row[0])
    f = max(1.0, 0.75 * d)
    oh, ow = max(1, int((r1 - r0) / f)), max(1, int((c1 - c0) / f))
    win = rasterio.windows.Window(c0, r0, c1 - c0, r1 - r0)
    nb = [1, 2, 3] if ds.count >= 3 else [1, 1, 1]
    a = ds.read(nb, window=win, out_shape=(3, oh, ow))
    alpha = ds.read(4, window=win, out_shape=(oh, ow)) if ds.count >= 4 else None
    rr = np.clip(((row - r0) * oh / (r1 - r0)).astype(int), 0, oh - 1)
    cc = np.clip(((col - c0) * ow / (c1 - c0)).astype(int), 0, ow - 1)
    out = a[:, rr, cc].reshape(3, n_v, n_l)
    if a.dtype != np.uint8:
        hi = float(np.percentile(out, 99.5)) or 1.0
        out = np.clip(out.astype(np.float32) / hi * 255, 0, 255)
    out = out.astype(np.uint8)
    if alpha is not None:
        out[:, alpha[rr, cc].reshape(n_v, n_l) == 0] = 0
    return out


def detect_rows_in_grid(rings, ids, ortho_path, grid_crs=None, rows_expected=0,
                        pitch_m=None, along="long", step_mm=5.0, pad_frac=0.15,
                        log=print, progress=None, cancel=None):
    """Find the rows of every plot in the orthomosaic. Returns (rows, report):
    rows as rows_from_count, report = one dict per plot (Plot_ID, ok, n_rows, pitch_m,
    angle_deg, reason)."""
    import rasterio
    from row_detect import detect_rows
    M, fr = grid_frames(rings, along)
    step = step_mm / 1000.0
    rows, report = [], []
    with rasterio.open(ortho_path) as ds:
        to_ds = None
        if ds.crs is not None:
            from rasterio.crs import CRS
            src = None
            if grid_crs:
                for fn in (CRS.from_wkt, CRS.from_user_input):
                    try:
                        src = fn(grid_crs)
                        break
                    except Exception:
                        pass
            if src is None and M.geographic:
                src = CRS.from_epsg(4326)
            if src is not None and src != ds.crs:
                from rasterio.warp import transform as _tf
                to_ds = (lambda x, y, _s=src, _d=ds.crs: _tf(_s, _d, list(x), list(y)))
        try:
            if ds.profile.get("blockysize", 512) == 1 or (ds.block_shapes and ds.block_shapes[0][0] == 1):
                log("[!] the orthomosaic is striped (one row per block): reading plots is slow. "
                    "'Make tiled copy' on the Plot Extraction tab makes a fast copy.")
        except Exception:
            pass
        for i, (pid, (c, a_u, a_v, hl, hw)) in enumerate(zip(ids, fr)):
            if cancel is not None and cancel():
                break
            rep = dict(Plot_ID=str(pid), ok=False, n_rows=0, pitch_m="", angle_deg="", reason="")
            pad = pad_frac * 2 * hw
            try:
                arr = _plot_window(ds, to_ds, M, c, a_u, a_v, hl, hw, pad, step)
            except Exception as e:
                arr, rep["reason"] = None, f"read failed: {e}"
            if arr is None:
                rep["reason"] = rep["reason"] or "plot outside the orthomosaic"
                report.append(rep)
                continue
            det = detect_rows(arr, across_axis=0, mm_px=step_mm,
                              pitch_mm=(pitch_m * 1000.0) if pitch_m else None,
                              rows_expected=int(rows_expected or 0))
            if not det.get("ok"):
                rep["reason"] = det.get("reason", "rows not found")
                report.append(rep)
                continue
            ang = float(det.get("angle_deg", 0.0))
            pr = det["plot_rows"]
            for k, j in enumerate(pr):
                v = -hw - pad + det["rows"][j] * step
                # along_centre is the image centre along: the reference line = plot centre
                rows.append(dict(Plot_ID=str(pid), row_idx=k + 1, v_off_m=round(v, 4),
                                 length_m=round(2 * hl, 3), angle_deg=round(ang, 3),
                                 source="detect", contrast=round(float(det["contrast"][j]), 4),
                                 line=_line(M, c, a_u, a_v, hl, v, ang)))
            rep.update(ok=True, n_rows=len(pr), pitch_m=round(det["pitch_px"] * step, 4),
                       angle_deg=round(ang, 2),
                       reason="" if all(det["present"][j] for j in pr) else "weak row(s)")
            report.append(rep)
            if progress:
                progress(i + 1, len(rings))
    return rows, report


# ---------------------------------------------------------------------------
# editing helpers
# ---------------------------------------------------------------------------
def relayout_plot(rows_of_plot, M, frame):
    """Recompute the line of every row of one plot from its v_off_m / angle_deg."""
    c, a_u, a_v, hl, hw = frame
    for r in rows_of_plot:
        r["length_m"] = round(2 * hl, 3)
        r["line"] = _line(M, c, a_u, a_v, hl, r["v_off_m"], r.get("angle_deg") or 0.0)


def renumber(rows):
    """row_idx 1..n by v_off_m inside every plot."""
    by = {}
    for r in rows:
        by.setdefault(r["Plot_ID"], []).append(r)
    for lst in by.values():
        lst.sort(key=lambda r: r["v_off_m"])
        for k, r in enumerate(lst):
            r["row_idx"] = k + 1


def even_spacing(rows_of_plot):
    """Space the rows of one plot evenly between its outer rows."""
    lst = sorted(rows_of_plot, key=lambda r: r["v_off_m"])
    if len(lst) < 3:
        return
    a, b = lst[0]["v_off_m"], lst[-1]["v_off_m"]
    for k, r in enumerate(lst):
        r["v_off_m"] = round(a + (b - a) * k / (len(lst) - 1), 4)
        if r["source"] != "count":
            r["source"] = "edit"


# ---------------------------------------------------------------------------
# files
# ---------------------------------------------------------------------------
def rows_path_for(grid_path):
    return os.path.splitext(grid_path)[0] + "_rows.shp"


def write_rows(path, rows, prj_from=None, report=None):
    """<path> polyline shapefile (+ .prj copied from prj_from's .prj) and
    <path>.csv with the rows per plot."""
    import shapefile
    w = shapefile.Writer(path, shapeType=shapefile.POLYLINE)
    w.field("Plot_ID", "C", 40)
    w.field("row_idx", "N", 4, 0)
    w.field("v_off_m", "N", 12, 4)
    w.field("length_m", "N", 10, 3)
    w.field("angle_deg", "N", 8, 3)
    w.field("source", "C", 10)
    w.field("contrast", "N", 10, 4)
    try:
        for r in rows:
            w.line([[list(p) for p in r["line"]]])
            ct = r.get("contrast")
            w.record(str(r["Plot_ID"]), int(r["row_idx"]), float(r["v_off_m"]),
                     float(r["length_m"]), float(r.get("angle_deg") or 0.0), r.get("source", ""),
                     None if ct is None or (isinstance(ct, float) and math.isnan(ct)) else float(ct))
    finally:
        w.close()
    if prj_from:
        src = os.path.splitext(prj_from)[0] + ".prj"
        if os.path.exists(src):
            with open(src) as f, open(os.path.splitext(path)[0] + ".prj", "w") as g:
                g.write(f.read())
    by = {}
    for r in rows:
        by.setdefault(r["Plot_ID"], []).append(r)
    rep = {d["Plot_ID"]: d for d in (report or [])}
    with open(os.path.splitext(path)[0] + ".csv", "w", newline="") as f:
        cw = csv.writer(f)
        cw.writerow(["Plot_ID", "n_rows", "spacing_m", "rows_v_off_m", "source", "note"])
        for pid in sorted(set(by) | set(rep), key=lambda s: (len(s), s)):
            lst = sorted(by.get(pid, []), key=lambda r: r["v_off_m"])
            vs = [r["v_off_m"] for r in lst]
            sp = float(np.median(np.diff(vs))) if len(vs) > 1 else ""
            cw.writerow([pid, len(lst), round(sp, 4) if sp != "" else "",
                         ";".join(f"{v:.3f}" for v in vs),
                         "/".join(sorted({r["source"] for r in lst})),
                         rep.get(pid, {}).get("reason", "")])
    return path


def read_rows(path):
    import shapefile
    rdr = shapefile.Reader(path)
    names = [f[0] for f in rdr.fields[1:]]
    rows = []
    for sr in rdr.iterShapeRecords():
        d = dict(zip(names, sr.record))
        pts = [tuple(p[:2]) for p in sr.shape.points]
        if len(pts) < 2:
            continue
        ct = d.get("contrast")
        rows.append(dict(Plot_ID=str(d.get("Plot_ID")).strip(), row_idx=int(d.get("row_idx") or 0),
                         v_off_m=float(d.get("v_off_m") or 0.0),
                         length_m=float(d.get("length_m") or 0.0),
                         angle_deg=float(d.get("angle_deg") or 0.0),
                         source=str(d.get("source") or "edit").strip(),
                         contrast=float(ct) if ct not in (None, "") else float("nan"),
                         line=[pts[0], pts[-1]]))
    return rows


def replace_count_rows(old_rows_path, new_grid_path, rings, ids, along="long"):
    """After the plots were moved / resized: rows that came from a count are placed
    again in the edited plots (spacing = new width / rows); detected or hand-edited
    rows stay where they are on the ground. Writes <new grid>_rows.shp."""
    old = read_rows(old_rows_path)
    if not old:
        return None
    counts = {}
    for r in old:
        counts[r["Plot_ID"]] = counts.get(r["Plot_ID"], 0) + 1
    if all(r["source"] == "count" for r in old):
        rows, _ = rows_from_count(rings, ids, [counts.get(str(i), 0) for i in ids], along=along)
    else:
        rows = old
    return write_rows(rows_path_for(new_grid_path), rows, prj_from=new_grid_path)
