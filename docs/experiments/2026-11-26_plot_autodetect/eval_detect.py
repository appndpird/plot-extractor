"""Score plot_detect against hand-checked grids.

usage: python eval_detect.py <case> [noquad] [counts]
"""
import json
import math
import os
import sys
import time

import numpy as np
import shapefile

APP = r"E:\Software\PlotExtractor-2026-11-26"
sys.path.insert(0, APP)
import plot_detect  # noqa: E402
from crop_rows import Metric, plot_frame  # noqa: E402

CASES = {
    "york2025": dict(
        ortho=r"F:\Ibrahim's Workspace2\2025_Calibrationtrial_I_SIF_DPIRD\Wheat_CalibrationTrial_York_20250625_Orthomosaic_tiled.tif",
        ref=r"F:\Ibrahim's Workspace2\2025_Calibrationtrial_I_SIF_DPIRD\PlotShapefile\edited_grid_edited_edited_edited_edited.shp"),
    "york2026": dict(
        ortho=r"F:\Ibrahim's Workspace2\Plot extraction processing\2026_CalibrationTrial_I_SIF_DPIRD\orthomosaic and point cloud\20260812_CalibrationTrial_SIF_York_Gobi_run00_YCD_Row1_RGB_Orthomosaic_10mm.tif",
        ref=r"F:\Ibrahim's Workspace2\Plot extraction processing\2026_CalibrationTrial_I_SIF_DPIRD\plot_gridrgb1_corrected.shp"),
    "northam": dict(
        ortho=r"I:\APPN2026\2026_NUEWheat_I_DPIRD\2026Northam_F\HIRES\20260703\run00\T1_proc\Orthomosaic\20260703_NUEWheat_Northam_I_DPIRD_Orthomosaic_7850.tif",
        ref=r"F:\Ibrahim's Workspace2\Plot extraction processing\2026_NUEWheat_I_DPIRD\PlotExtractor_2026-10-25_Northam_20260703\NUE_Northam_2026_fitted_grid.shp"),
}


def rings_of(path):
    r = shapefile.Reader(path)
    return [s.points[:4] for s in r.shapes()], [str(x[0]) for x in r.records()]


def ref_quad(rings, jitter=0.25, seed=1):
    """trial corners BL BR TR TL from a reference grid, jittered like hand clicks"""
    M = Metric(rings)
    P = np.array([M.to_m(p) for r in rings for p in r])
    fr = [plot_frame([M.to_m(p) for p in r]) for r in rings]
    au = fr[0][1]
    if au[1] < 0:
        au = -au
    av = np.array([au[1], -au[0]])          # across = along turned -90
    a, c = P @ au, P @ av
    rng = np.random.default_rng(seed)
    q = []
    for aa, cc in ((a.min(), c.min()), (a.max() * 0 + a.min(), c.max()), (a.max(), c.max()), (a.max(), c.min())):
        pt = aa * au + cc * av + rng.normal(0, jitter, 2)
        q.append(M.from_m(pt))
    # order BL, BR, TR, TL
    return [q[0], q[1], q[2], q[3]]


def score(det_rings, ref_rings):
    M = Metric(ref_rings)
    dc = [M.to_m(np.mean(r, 0)) for r in det_rings]
    dc = np.array(dc)
    out = []
    used = set()
    for r in ref_rings:
        c, au, av, hl, hw = plot_frame([M.to_m(p) for p in r])
        d = np.hypot(*(dc - c).T)
        k = int(np.argmin(d))
        if d[k] > hw:
            out.append(None)
            continue
        used.add(k)
        c2, au2, av2, hl2, hw2 = plot_frame([M.to_m(p) for p in det_rings[k]])
        off = c2 - c
        out.append(dict(across=float(off @ av), along=float(off @ au), dw=2 * (hw2 - hw), dl=2 * (hl2 - hl),
                        iou=iou([M.to_m(p) for p in r], [M.to_m(p) for p in det_rings[k]])))
    return out, len(det_rings) - len(used)


def iou(a, b):
    from shapely.geometry import Polygon
    A, B = Polygon(a), Polygon(b)
    return A.intersection(B).area / max(A.union(B).area, 1e-9)


def main():
    case = sys.argv[1]
    noquad = "noquad" in sys.argv
    counts = "counts" in sys.argv
    C = CASES[case]
    ref_rings, ref_ids = rings_of(C["ref"])
    quad = None if noquad else ref_quad(ref_rings)
    t = time.time()
    res = plot_detect.detect_plots(C["ortho"], quad,
                                   n_ranges=int(sys.argv[sys.argv.index("counts") + 1]) if counts else 0,
                                   n_cols=int(sys.argv[sys.argv.index("counts") + 2]) if counts else 0)
    det = [p["ring"] for p in res["plots"]]
    sc, extra = score(det, ref_rings)
    ok = [s for s in sc if s]
    A = np.array([[s["across"], s["along"], s["dw"], s["dl"], s["iou"]] for s in ok])
    print(json.dumps(res["info"]))
    print("%s quad=%s: ref %d, detected %d, matched %d, extra %d, %.1f s" % (
        case, not noquad, len(ref_rings), len(det), len(ok), extra, time.time() - t))
    if len(A):
        for k, nm in enumerate(["across", "along", "d_width", "d_length"]):
            print("  %-9s median %+6.1f cm  |abs| median %5.1f  p90 %5.1f  max %5.1f" % (
                nm, 100 * np.median(A[:, k]), 100 * np.median(np.abs(A[:, k])),
                100 * np.percentile(np.abs(A[:, k]), 90), 100 * np.abs(A[:, k]).max()))
        print("  IoU median %.3f  min %.3f" % (np.median(A[:, 4]), A[:, 4].min()))
    worst = sorted([(abs(x["along"]) + abs(x["across"]), rid, x) for rid, x in zip(ref_ids, sc) if x], reverse=True)[:10]
    print("  worst:", "; ".join("%s a%+.0f c%+.0f" % (rid, 100 * x["along"], 100 * x["across"]) for _, rid, x in worst))
    print("  unmatched ref:", [rid for rid, x in zip(ref_ids, sc) if not x][:20])
    fl = [p for p in res["plots"] if p["conf"] < plot_detect.FLAG_BELOW]
    print("  flagged %d:" % len(fl), "; ".join("%s %s" % (p["props"]["Plot_ID"], ",".join(p["flags"])) for p in fl[:12]))
    od = os.path.join(os.path.dirname(__file__), "out")
    os.makedirs(od, exist_ok=True)
    from grid_designer import write_grid_shapefile
    tag = case + ("_noquad" if noquad else "") + ("_counts" if counts else "")
    write_grid_shapefile(os.path.join(od, tag + ".shp"), res["plots"], res["crs"])
    plot_detect.write_report(os.path.join(od, tag + "_detect.csv"), res)


if __name__ == "__main__":
    main()
