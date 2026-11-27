"""draw reference (red) and detected (cyan) plots on the orthomosaic, trial-wide"""
import os
import sys

import numpy as np
import rasterio
from PIL import Image, ImageDraw
from rasterio.enums import Resampling

from eval_detect import CASES, rings_of

case, tag = sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else sys.argv[1]
zoom = sys.argv[3] if len(sys.argv) > 3 else None     # "plotid" to zoom to a ref plot
C = CASES[case]
ref, ref_ids = rings_of(C["ref"])
det, det_ids = rings_of(os.path.join(os.path.dirname(__file__), "out", tag + ".shp"))
allp = np.array([p for r in ref + det for p in r])
x0, y0 = allp.min(0)
x1, y1 = allp.max(0)
if zoom:
    r = np.array(ref[ref_ids.index(zoom)])
    cx, cy = r.mean(0)
    sx, sy = (r.max(0) - r.min(0)) * 1.2
    s = max(sx, sy) * float(sys.argv[4]) if len(sys.argv) > 4 else max(sx, sy)
    x0, x1, y0, y1 = cx - s, cx + s, cy - s, cy + s
pw, ph = x1 - x0, y1 - y0
x0, x1, y0, y1 = x0 - 0.03 * pw, x1 + 0.03 * pw, y0 - 0.03 * ph, y1 + 0.03 * ph
with rasterio.open(C["ortho"]) as ds:
    win = rasterio.windows.from_bounds(x0, y0, x1, y1, ds.transform)
    W = 1800
    H = int(W * win.height / win.width)
    a = ds.read([1, 2, 3], window=win, out_shape=(3, H, W), resampling=Resampling.average, boundless=True)
im = Image.fromarray(np.transpose(a, (1, 2, 0)))
d = ImageDraw.Draw(im)


def px(p):
    return ((p[0] - x0) / (x1 - x0) * W, (y1 - p[1]) / (y1 - y0) * H)


for r in ref:
    d.polygon([px(p) for p in r], outline=(255, 40, 40), width=2)
for r, i in zip(det, det_ids):
    d.polygon([px(p) for p in r], outline=(0, 230, 255), width=2)
    c = px(np.mean(r, 0))
    d.text(c, i, fill=(255, 255, 0))
out = os.path.join(os.path.dirname(__file__), "diag", tag + ("_z" + zoom if zoom else "") + "_overlay.jpg")
im.save(out, quality=85)
print(out)
