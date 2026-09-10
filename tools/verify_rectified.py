"""Coverage verification for rectified crops: % of unfilled (pure black)
pixels per plot + georeferencing sanity. Writes rectified_qc.csv per day."""
import csv
import glob
import math
import os
import sys

import numpy as np
import rasterio

ROOT = (r"F:\Ibrahim's Workspace2\Plant Emergency Dataset UWA"
        r"\2025_NUE+RootPhenotyping_I_DPIRD-20260902\Rectified")

for day in (sys.argv[1:] or ["Day1", "Day2"]):
    folder = os.path.join(ROOT, day, "RawCrops")
    rows = []
    for f in sorted(glob.glob(os.path.join(folder, "*.tif"))):
        pid = os.path.basename(f).replace("plot_", "").replace(".tif", "")
        with rasterio.open(f) as ds:
            # full-res exact black count via decimation is lossy; read at /4
            # with nearest sampling - black regions are contiguous, so the
            # estimate is accurate to well under 0.1%
            a = ds.read(out_shape=(3, max(1, ds.height // 4),
                                   max(1, ds.width // 4)))
            black = float((a.max(axis=0) == 0).mean() * 100.0)
            t = ds.transform
            gsd_mm = (math.hypot(t.a, t.d) + math.hypot(t.b, t.e)) / 2 * 1000
            tag = ds.tags().get("PE_AFFINE_RESID", "")
            has_geo = ds.crs is not None and not t.is_identity
            w_m = math.hypot(t.a * ds.width, t.d * ds.width)
            h_m = math.hypot(t.b * ds.height, t.e * ds.height)
        rows.append(dict(
            plot_id=pid, width_px=ds.width, height_px=ds.height,
            gsd_mm_per_px=round(gsd_mm, 3),
            size_m=f"{w_m:.2f}x{h_m:.2f}",
            black_pct=round(black, 2),
            covered_pct=round(100 - black, 2),
            georeferenced=("yes" if has_geo else "NO"),
            resid=tag,
            flag=("INCOMPLETE >1% missing" if black > 1.0 else
                  "check" if black > 0.2 else "ok"),
        ))
    out = os.path.join(ROOT, day, "rectified_qc.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    b = np.array([r["black_pct"] for r in rows])
    geo = sum(1 for r in rows if r["georeferenced"] == "yes")
    print(f"\n=== {day}: {len(rows)} rectified crops ===")
    print(f"georeferenced: {geo}/{len(rows)}")
    print(f"missing-area %: p50={np.percentile(b,50):.2f} "
          f"p90={np.percentile(b,90):.2f} max={b.max():.2f}")
    print(f"fully covered (<=0.2% black): "
          f"{sum(1 for r in rows if r['black_pct'] <= 0.2)}")
    bad = [r for r in rows if r["black_pct"] > 1.0]
    print(f"INCOMPLETE (>1% missing): {len(bad)}")
    for r in sorted(bad, key=lambda r: -r["black_pct"])[:15]:
        print(f"  plot {r['plot_id']}: {r['black_pct']:.1f}% missing "
              f"({r['size_m']} m, {r['gsd_mm_per_px']} mm/px)")
    print(f"wrote {out}")
