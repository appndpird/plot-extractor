"""Experiment 2026-11-10, step 1: direct IIQ plot extraction of York 2025 (36 plots) in every
output format PlotExtractor offers, with the same frames and boxes in every run.

Runs (bounding box, margin 0.10 m - no mask, so every pixel of the box is kept):
  Neutral preset (no noise reduction, no sharpening: the SDK development is bit-repeatable, so any
  difference between runs is the file format):
    N16_geotiff  N16_tiff  N16_tiff_none  N16_png  N8_geotiff  N8_tiff  N8_png  N8_jpeg
  APPN HR preset (the app's default development):
    HR16_geotiff  HR8_geotiff
Only reads the project and the IIQ files. Python 3.11 (Metashape wheel):
  python exp1_run_formats.py [run names...]
"""
import json, os, sys, threading, time
APP = r"E:\Software\PlotExtractor-2026-11-10"
sys.path.insert(0, APP)
import iiq

ROOT = r"F:\Ibrahim's Workspace2\Plot extraction processing\2025_Calibrationtrial_I_SIF_DPIRD\Experiment_2026-11-10_IIQ_formats"
RUN00 = r"G:\APPN 2025\2025_Calibrationtrial_I_SIF_DPIRD\2025-York_F\PhaseOne\2025-06-25\run00"
BASE = dict(psx=os.path.join(RUN00, r"T1_Proc\MetaShapeProcessing.psx"), chunk="",
            plots=r"F:\Ibrahim's Workspace2\2025_Calibrationtrial_I_SIF_DPIRD\PlotShapefile\edited_grid_edited_edited_edited_edited.shp",
            id_field="Plot_ID", rawdir=os.path.join(RUN00, "T0_Raw"),
            margin="0.10", n_per_plot="1", limit="0", do_raw=True, do_ortho=False,
            crs_override="auto", raw_mode="bbox", to_srgb=False, workers="4", auto_fit=False,
            frame_choice="whole_first", raw_source="iiq", gcp_epsg="7850",
            iiq_sdk=os.path.join(APP, r"iiq_tools\sdk"), quality="95")
FMT = {"geotiff": "GeoTIFF (geotagged, lossless)", "tiff": "TIFF (lossless)",
       "png": "PNG (lossless)", "jpeg": "JPEG"}
RUNS = {}
for preset, tag in (("neutral", "N"), ("appn_hr", "HR")):
    for bits in (16, 8):
        for f in ("geotiff", "tiff", "png", "jpeg"):
            if bits == 16 and f == "jpeg":
                continue
            RUNS[f"{tag}{bits}_{f}"] = dict(preset=preset, bits=bits, fmt=f, compress="lzw")
RUNS["N16_tiff_none"] = dict(preset="neutral", bits=16, fmt="tiff", compress="none")
DEFAULT = ["N16_geotiff", "N16_tiff", "N16_tiff_none", "N16_png", "N8_geotiff", "N8_tiff",
           "N8_png", "N8_jpeg", "HR16_geotiff", "HR8_geotiff"]

if __name__ == "__main__":
    names = sys.argv[1:] or DEFAULT
    os.makedirs(ROOT, exist_ok=True)
    times = {}
    tf = os.path.join(ROOT, "run_times.json")
    if os.path.exists(tf):
        times = json.load(open(tf))
    for name in names:
        c = RUNS[name]
        out = os.path.join(ROOT, name)
        if os.path.exists(os.path.join(out, "manifest.csv")):
            print("skip (done)", name); continue
        p = dict(BASE, outdir=out, format=FMT[c["fmt"]], bits=str(c["bits"]), compress=c["compress"],
                 iiq_settings=os.path.join(APP, r"iiq_tools\presets", c["preset"] + ".json"))
        os.makedirs(out, exist_ok=True)
        logf = open(os.path.join(out, "run_log.txt"), "w", encoding="utf-8")
        def log(m, f=logf):
            f.write(time.strftime("%H:%M:%S ") + str(m) + "\n"); f.flush()
        t = time.time()
        print("RUN", name, flush=True)
        iiq.run_iiq_extraction(p, log, lambda a, b: None, threading.Event())
        times[name] = round(time.time() - t, 1)
        logf.close()
        json.dump(times, open(tf, "w"), indent=1)
        print("  done", name, times[name], "s", flush=True)
