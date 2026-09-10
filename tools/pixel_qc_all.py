"""Bit-exact pixel identity check: every RawCrops TIF vs its source JPEG.

Finds each crop's offset in the source frame by template match, then compares
the full crop. Writes pixel_qc.csv incrementally next to manifest.csv.
Workers keep a one-deep cache of the last decoded JPEG (files are dispatched
sorted by source so plots sharing a frame reuse the decode).
"""
import csv
import glob
import os
import sys
import time
from multiprocessing import Pool

BASE = r"F:\Ibrahim's Workspace2\Plant Emergency Dataset UWA\2025_NUE+RootPhenotyping_I_DPIRD-20263108"
OUT = os.path.join(BASE, "pixel_qc.csv")

_cache = {"path": None, "rgb": None, "gray": None}


def _load_source(path):
    import cv2
    import numpy as np
    from PIL import Image
    if _cache["path"] != path:
        Image.MAX_IMAGE_PIXELS = None
        rgb = np.asarray(Image.open(path))
        _cache.update(path=path, rgb=rgb,
                      gray=cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY))
    return _cache["rgb"], _cache["gray"]


def check(tif):
    import cv2
    import numpy as np
    import rasterio
    t0 = time.time()
    row = dict(file=os.path.basename(tif), status="", source="", x0="", y0="",
               diff_pixels="", max_abs_diff="", secs="")
    try:
        with rasterio.open(tif) as ds:
            src = ds.tags().get("SOURCE_IMAGE", "")
            crop = ds.read().transpose(1, 2, 0)
        row["source"] = os.path.basename(src)
        if not src or not os.path.exists(src):
            row["status"] = "SOURCE_MISSING"
            return row
        full, gray = _load_source(src)
        h, w = crop.shape[:2]
        cy, cx = h // 2, w // 2
        tpl = cv2.cvtColor(crop[cy - 64:cy + 64, cx - 64:cx + 64],
                           cv2.COLOR_RGB2GRAY)
        r = cv2.matchTemplate(gray, tpl, cv2.TM_SQDIFF)
        _, _, mn, _ = cv2.minMaxLoc(r)
        x0, y0 = mn[0] - (cx - 64), mn[1] - (cy - 64)
        row["x0"], row["y0"] = x0, y0
        if x0 < 0 or y0 < 0 or x0 + w > full.shape[1] or y0 + h > full.shape[0]:
            row["status"] = "OFFSET_OUT_OF_BOUNDS"
            return row
        region = full[y0:y0 + h, x0:x0 + w]
        if np.array_equal(region, crop):
            row["status"], row["diff_pixels"], row["max_abs_diff"] = "IDENTICAL", 0, 0
        else:
            nz = mx = 0
            for i in range(0, h, 512):          # chunked stats, bounded memory
                d = region[i:i + 512].astype(np.int16) - crop[i:i + 512]
                nz += int(np.count_nonzero(d.any(axis=2)))
                mx = max(mx, int(np.abs(d).max()))
            row["status"] = "DIFFERS"
            row["diff_pixels"], row["max_abs_diff"] = nz, mx
    except Exception as e:
        row["status"] = f"ERROR {type(e).__name__}: {e}"
    finally:
        row["secs"] = f"{time.time() - t0:.1f}"
    return row


def source_of(tif):
    import rasterio
    with rasterio.open(tif) as ds:
        return ds.tags().get("SOURCE_IMAGE", "")


if __name__ == "__main__":
    files = sorted(glob.glob(os.path.join(BASE, "RawCrops", "*.tif")))
    print(f"{len(files)} files; reading source mapping...", flush=True)
    files.sort(key=lambda f: (source_of(f), f))
    fields = ["file", "status", "source", "x0", "y0",
              "diff_pixels", "max_abs_diff", "secs"]
    done = 0
    t0 = time.time()
    with open(OUT, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=fields)
        wr.writeheader()
        with Pool(4) as pool:
            for row in pool.imap(check, files, chunksize=4):
                wr.writerow(row)
                fh.flush()
                done += 1
                if row["status"] != "IDENTICAL":
                    print(f"[!] {row['file']}: {row['status']} "
                          f"diff={row['diff_pixels']}", flush=True)
                if done % 20 == 0:
                    el = time.time() - t0
                    print(f"{done}/{len(files)}  {el/60:.1f} min elapsed, "
                          f"eta {el/done*(len(files)-done)/60:.0f} min", flush=True)
    print(f"\nDone in {(time.time()-t0)/60:.1f} min. Results: {OUT}", flush=True)
