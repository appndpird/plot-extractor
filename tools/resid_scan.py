import glob
import os

import numpy as np
import rasterio

BASE = r"F:\Ibrahim's Workspace2\Plant Emergency Dataset UWA\2025_NUE+RootPhenotyping_I_DPIRD-20263108\RawCrops"

rows = []
for f in sorted(glob.glob(os.path.join(BASE, "*.tif"))):
    with rasterio.open(f) as ds:
        tag = ds.tags().get("PE_AFFINE_RESID", "")
    if not tag:
        continue
    mean = float(tag.split("mean=")[1].split("mm")[0])
    mx = float(tag.split("max=")[1].split("mm")[0])
    rows.append((os.path.basename(f), mean, mx))

mx = np.array([r[2] for r in rows])
print(f"{len(rows)} files with residual tags")
print(f"max-residual percentiles (mm): p50={np.percentile(mx,50):.1f} "
      f"p90={np.percentile(mx,90):.1f} p99={np.percentile(mx,99):.1f} "
      f"max={mx.max():.1f}")
for thr in (10, 25, 50, 100):
    print(f"files with max residual > {thr} mm: {(mx > thr).sum()}")
print("\nWorst 10:")
for name, mean, m in sorted(rows, key=lambda r: -r[2])[:10]:
    print(f"  {name}: mean={mean:.1f} mm max={m:.1f} mm")
