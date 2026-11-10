"""Experiment 2026-11-10, step 7: the two extraction routes of tab 3 give the same pixels, and
the raw-sensor output of tab 3b is exact.

  a) tab 2: two IIQ frames converted to 16-bit TIFF (Neutral), full size and 0.5 scale
  b) tab 3a (indirect) from the full-size TIFFs, 16-bit, bounding box -> crops of the first
     plots == the tab 3b crops (N16_geotiff, same Neutral development)?
  c) tab 3a from the 0.5-scale TIFFs -> cameras scaled x0.5, crop = 3b crop at half size?
  d) tab 3b, Pixels = raw sensor values, bounding box GeoTIFF, all 36 plots; QA (tab 4c
     check) reads every box again from the IIQ -> IDENTICAL?
"""
import csv, os, shutil, sys, threading, time
import numpy as np
import rasterio

APP = r"E:\Software\PlotExtractor-2026-11-10"
sys.path.insert(0, APP)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import iiq, qa
import plot_extractor as pe
from exp1_run_formats import BASE, ROOT
from exp5_raw_sensor import phase_shift

AN = os.path.join(ROOT, "analysis")
W = os.path.join(ROOT, "indirect")
FR = ["P0019790", "P0019784"]
L = []


def log_to(path):
    f = open(path, "w", encoding="utf-8")
    return lambda m: (f.write(time.strftime("%H:%M:%S ") + str(m) + "\n"), f.flush())


# a) conversion
src = os.path.join(W, "iiq_copy")
os.makedirs(src, exist_ok=True)
for fr in FR:
    d = os.path.join(src, fr + ".IIQ")
    if not os.path.exists(d):
        shutil.copy2(os.path.join(BASE["rawdir"], fr + ".IIQ"), d)
for tag, scale in (("tif16_full", 1.0), ("tif16_half", 0.5)):
    out = os.path.join(W, tag)
    if not os.path.isdir(out) or len([f for f in os.listdir(out) if f.endswith(".tif")]) < 2:
        iiq.convert_folder(dict(iiq_dir=src, recursive=False, out_dir=out, fmt="tiff", bits=16,
                                lzw=True, scale=scale, quality=95,
                                settings=os.path.join(APP, r"iiq_tools\presets\neutral.json"),
                                sdk=BASE["iiq_sdk"]), log_to(os.path.join(W, tag + "_log.txt")),
                           lambda a, b: None, threading.Event())
    with rasterio.open(os.path.join(out, FR[0] + ".tif")) as d:
        L.append(f"a) tab 2 conversion {tag}: {d.width}x{d.height} {d.dtypes[0]}")

# b, c) tab 3a from the converted TIFFs
n16 = {m["plot_id"]: m for m in csv.DictReader(open(os.path.join(ROOT, "N16_geotiff", "manifest.csv")))}
for tag in ("tif16_full", "tif16_half"):
    out = os.path.join(W, "3a_" + tag)
    if not os.path.exists(os.path.join(out, "manifest.csv")):
        p = dict(BASE, outdir=out, rawdir=os.path.join(W, tag), raw_source="image",
                 format="TIFF (lossless)", bits="16", compress="deflate", limit="8")
        pe.run_extraction(p, log_to(os.path.join(W, "3a_" + tag + "_log.txt")),
                          lambda a, b: None, threading.Event())
    man = list(csv.DictReader(open(os.path.join(out, "manifest.csv"))))
    same, total, notes = 0, 0, []
    for m in man:
        if m["raw_image"] not in FR or m["raw_image"] != n16[m["plot_id"]]["raw_image"]:
            continue
        with rasterio.open(os.path.join(out, "RawCrops", m["raw_file"])) as d:
            a = d.read()
        with rasterio.open(os.path.join(ROOT, "N16_geotiff", "RawCrops", f"plot_{m['plot_id']}.tif")) as d:
            b = d.read()
        total += 1
        if tag == "tif16_full":
            same += int(np.array_equal(a, b))
            notes.append(f"{m['plot_id']} {a.dtype} {a.shape[2]}x{a.shape[1]} vs {b.shape[2]}x{b.shape[1]}")
        else:
            # half-size: compare with the 3b crop averaged 2x2 (green, log), correlation
            h, w = min(a.shape[1], b.shape[1] // 2), min(a.shape[2], b.shape[2] // 2)
            bb = b[1, :2 * h, :2 * w].reshape(h, 2, w, 2).mean(axis=(1, 3))
            aa = np.log1p(a[1, :h, :w].astype(float))
            cy, cx, k = h // 2, w // 2, min(256, h // 2, w // 2)
            dy, dx, _pk = phase_shift(aa[cy - k:cy + k, cx - k:cx + k], np.log1p(bb)[cy - k:cy + k, cx - k:cx + k])
            off = 2 * max(abs(dx), abs(dy))            # in full-size pixels
            same += int(off <= 1.0)
            notes.append(f"{m['plot_id']} {a.shape[2]}x{a.shape[1]} vs 3b {b.shape[2]}x{b.shape[1]} offset {off:.2f} px")
    L.append(f"{'b' if tag == 'tif16_full' else 'c'}) tab 3a from {tag}: "
             f"{same}/{total} plots {'pixel-identical to the 3b crops' if tag == 'tif16_full' else 'line up with the 3b crops (offset <= 1 full-size px)'}"
             f"; {'; '.join(notes)}")
    L.append("   log: " + " | ".join(l.strip()[9:] for l in open(os.path.join(W, "3a_" + tag + "_log.txt"))
                                       if "scaled" in l or "16-bit" in l))

# d) tab 3b raw sensor values + QA
out = os.path.join(ROOT, "RAW_sensor_geotiff")
if not os.path.exists(os.path.join(out, "manifest.csv")):
    t = time.time()
    p = dict(BASE, outdir=out, format="GeoTIFF (geotagged, lossless)", bits="16", compress="deflate",
             iiq_pixels="raw_sensor", raw_mode="bbox")
    iiq.run_iiq_extraction(p, log_to(os.path.join(ROOT, "RAW_sensor_run_log.txt")), lambda a, b: None,
                           threading.Event())
    L.append(f"d) tab 3b raw sensor values: {round(time.time() - t)} s")
res = qa.check_crops(out, log=lambda m: None)
from collections import Counter
L.append(f"d) QA re-read from the IIQs: {dict(Counter(r['result'] for r in res))}")
with rasterio.open(os.path.join(out, "RawCrops", "plot_6006.tif")) as d:
    L.append(f"d) plot_6006.tif: {d.count} band {d.dtypes[0]}, CRS {d.crs}, tags "
             f"{ {k: v for k, v in d.tags().items() if k.startswith('PE_CFA') or k == 'PIXEL_SOURCE'} }")
sz = sum(os.path.getsize(os.path.join(out, "RawCrops", f)) for f in os.listdir(os.path.join(out, "RawCrops"))) / 1e6
L.append(f"d) 36 raw-sensor GeoTIFFs (Deflate): {sz:.0f} MB")
txt = "\n".join(L)
open(os.path.join(AN, "INDIRECT_RAW_SUMMARY.txt"), "w", encoding="utf-8").write(txt)
print(txt)
