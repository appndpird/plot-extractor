"""Quality + geotag accuracy check for RawCrops plot TIFs.

1. Pixel identity: crop TIF vs source JPEG (find offset by template match,
   then compare the full overlap exactly).
2. Affine fit to the embedded GCPs (EPSG:7850): RMS/max residual in px and mm.
"""
import glob
import os
import sys

import numpy as np
import rasterio

BASE = r"F:\Ibrahim's Workspace2\Plant Emergency Dataset UWA\2025_NUE+RootPhenotyping_I_DPIRD-20263108"


def affine_fit(gcps):
    # least-squares affine (col,row) -> (E,N)
    A = np.array([[g.col, g.row, 1.0] for g in gcps])
    E = np.array([g.x for g in gcps])
    N = np.array([g.y for g in gcps])
    ce, *_ = np.linalg.lstsq(A, E, rcond=None)
    cn, *_ = np.linalg.lstsq(A, N, rcond=None)
    pe = A @ ce - E
    pn = A @ cn - N
    res = np.hypot(pe, pn)              # metres in 7850
    # pixel size from the affine
    px = np.hypot(ce[0], cn[0])
    py = np.hypot(ce[1], cn[1])
    return ce, cn, res, (px + py) / 2.0


def geo_report(paths):
    print("=== GCP affine-fit residuals ===")
    allres = []
    for p in paths:
        with rasterio.open(p) as ds:
            gcps, crs = ds.gcps
        if len(gcps) < 4:
            print(f"{os.path.basename(p)}: only {len(gcps)} GCPs - skipped")
            continue
        ce, cn, res, gsd = affine_fit(gcps)
        allres.append(res)
        print(f"{os.path.basename(p)}: {len(gcps)} GCPs crs={crs} "
              f"gsd={gsd*1000:.3f} mm/px  resid rms={res.std()+0:.4f} "
              f"mean={res.mean()*1000:.1f} mm max={res.max()*1000:.1f} mm "
              f"(={res.max()/gsd:.1f} px)")
    if allres:
        r = np.concatenate(allres)
        print(f"\nOverall: mean residual {r.mean()*1000:.1f} mm, "
              f"max {r.max()*1000:.1f} mm over {len(r)} GCPs / {len(allres)} files")


def pixel_check(tif_path):
    print("\n=== Pixel identity: crop vs source JPEG ===")
    import cv2
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    with rasterio.open(tif_path) as ds:
        src = ds.tags().get("SOURCE_IMAGE")
        crop = ds.read().transpose(1, 2, 0)
    print("crop:", crop.shape, "source:", src)
    if not src or not os.path.exists(src):
        print("source JPEG not accessible - skipping pixel check")
        return
    im = Image.open(src)
    print("JPEG size:", im.size, "mode:", im.mode,
          "subsampling:", getattr(im, "layer", None))
    try:
        from PIL import JpegImagePlugin
        print("JPEG chroma subsampling code:", JpegImagePlugin.get_sampling(im))
    except Exception as e:
        print("subsampling probe failed:", e)
    full = np.asarray(im)
    # locate crop offset with a template from the crop centre
    cy, cx = crop.shape[0] // 2, crop.shape[1] // 2
    tpl = crop[cy - 64:cy + 64, cx - 64:cx + 64]
    r = cv2.matchTemplate(cv2.cvtColor(full, cv2.COLOR_RGB2GRAY),
                          cv2.cvtColor(tpl, cv2.COLOR_RGB2GRAY),
                          cv2.TM_SQDIFF)
    _, _, mn, _ = cv2.minMaxLoc(r)
    y0 = mn[1] - (cy - 64)
    x0 = mn[0] - (cx - 64)
    print(f"crop offset in source frame: x0={x0}, y0={y0}, min sqdiff={r.min():.1f}")
    h, w = crop.shape[:2]
    if y0 < 0 or x0 < 0 or y0 + h > full.shape[0] or x0 + w > full.shape[1]:
        print("offset out of bounds - crop extends past source?")
        return
    region = full[y0:y0 + h, x0:x0 + w]
    diff = (region.astype(np.int16) - crop.astype(np.int16))
    nz = np.count_nonzero(diff.any(axis=2))
    print(f"differing pixels: {nz} / {h*w} ({100.0*nz/(h*w):.4f}%)  "
          f"max abs diff: {np.abs(diff).max()}")
    if nz == 0:
        print("=> crop is a BIT-EXACT copy of the source JPEG pixels (lossless)")


if __name__ == "__main__":
    files = sorted(glob.glob(os.path.join(BASE, "RawCrops", "plot_*.tif")))
    sample = files[:: max(1, len(files) // 8)][:8]
    geo_report(sample)
    pixel_check(files[0])
