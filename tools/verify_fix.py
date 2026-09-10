import sys

import numpy as np
import rasterio

orig, fixed = sys.argv[1], sys.argv[2]

with rasterio.open(fixed) as ds:
    print("crs:", ds.crs)
    print("transform:", tuple(round(v, 6) for v in ds.transform)[:6])
    print("gcps left:", len(ds.gcps[0]))
    print("bounds:", ds.bounds)
    print("resid tag:", ds.tags().get("PE_AFFINE_RESID"))
    print("res (m/px):", ds.res)
    a = ds.read()

with rasterio.open(orig) as ds:
    b = ds.read()

print("pixels identical:", bool(np.array_equal(a, b)))

# distance sanity: crop width in metres via the affine vs mm/px * px
import math
with rasterio.open(fixed) as ds:
    t = ds.transform
    w_m = math.hypot(t.a * ds.width, t.d * ds.width)
    h_m = math.hypot(t.b * ds.height, t.e * ds.height)
    print(f"crop ground size: {w_m:.3f} m x {h_m:.3f} m  "
          f"(area {w_m*h_m:.3f} m2)")
