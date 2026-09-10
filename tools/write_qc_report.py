"""Write georef_qc.csv next to manifest.csv: per-plot geotag accuracy after the
2026-09-01 GCP->affine fix, from the PE_AFFINE_RESID tags."""
import csv
import glob
import math
import os

import numpy as np
import rasterio

BASE = r"F:\Ibrahim's Workspace2\Plant Emergency Dataset UWA\2025_NUE+RootPhenotyping_I_DPIRD-20263108"
OUT = os.path.join(BASE, "georef_qc.csv")

rows = []
for f in sorted(glob.glob(os.path.join(BASE, "RawCrops", "*.tif"))):
    with rasterio.open(f) as ds:
        tag = ds.tags().get("PE_AFFINE_RESID", "")
        t = ds.transform
        gsd_mm = (math.hypot(t.a, t.d) + math.hypot(t.b, t.e)) / 2 * 1000
        crs = str(ds.crs)
        w, h = ds.width, ds.height
    mean = float(tag.split("mean=")[1].split("mm")[0]) if tag else None
    mx = float(tag.split("max=")[1].split("mm")[0]) if tag else None
    n = int(tag.rsplit("n=", 1)[1]) if tag else 0
    rows.append(dict(
        plot_id=os.path.basename(f).replace("plot_", "").replace(".tif", ""),
        file=os.path.basename(f), crs=crs, width_px=w, height_px=h,
        gsd_mm_per_px=round(gsd_mm, 4), n_gcps=n,
        resid_mean_mm=mean, resid_max_mm=mx,
        flag=("OBLIQUE_FRAME resid>50mm" if mx and mx > 50 else
              "check" if mx and mx > 25 else "ok"),
    ))

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

mx = np.array([r["resid_max_mm"] for r in rows if r["resid_max_mm"] is not None])
print(f"wrote {OUT} ({len(rows)} rows)")
print(f"max-residual mm: p50={np.percentile(mx,50):.1f} p90={np.percentile(mx,90):.1f} "
      f"max={mx.max():.1f}; flagged>50mm: {sum(1 for r in rows if 'OBLIQUE' in r['flag'])}")
