"""Rectified extraction of Day1 + Day2: each plot deskewed to an upright
rectangle, gap-filled from neighbour frames with min-error seams. Same psx /
chunks / shapefiles / CRS as the native recipe; margin 0."""
import os
import sys
import threading
import time
import traceback

sys.path.insert(0, r"E:\Software\PlotExtractor")
import plot_extractor as pe

ROOT = (r"F:\Ibrahim's Workspace2\Plant Emergency Dataset UWA"
        r"\2025_NUE+RootPhenotyping_I_DPIRD-20260902\Rectified")
T1 = (r"G:/APPN 2025/2025_NUE+RootPhenotyping_I_DPIRD/2025-Muresk_F/PhaseOne"
      r"/2025-06-09/run00/T1_Proc")
SHP = T1 + "/PlotShape/PlotShapes_EPSG7850"

COMMON = dict(
    source="metashape",
    psx=T1 + "/20250609_NUE_RootPhenotyping_I_DPIRD_Muresk_F_PhaseOne.psx",
    id_field="Plot_ID", margin=0.0, n_per_plot=1, limit=0, quality=100,
    do_raw=True, raw_mode="rectify", engine="metashape", do_ortho=False,
    to_srgb=False, workers=32, auto_fit=False, rows_per_plot=0,
    fit_mode="width", crs_override="auto", gcp_epsg="7850",
    format="GeoTIFF (geotagged, lossless)",
    smart_seam=True, seam_map=False, mesh_step=96,
)

RUNS = [
    ("Day1", "RootPheno_Day1",
     SHP + "/NUE_RT_Muersk_Day1_final_plot_shapes_7850.shp", T1 + "/JPEG/Day1"),
    ("Day2", "RootPheno_Day2",
     SHP + "/NUE_RT_Muersk_Day2_final_plot_shapes_7850.shp", T1 + "/JPEG/Day2"),
]
if len(sys.argv) > 1:
    RUNS = [r for r in RUNS if r[0] in sys.argv[1:]]

failures = []
for day, chunk, shp, rawdir in RUNS:
    outdir = os.path.join(ROOT, day)
    os.makedirs(outdir, exist_ok=True)
    p = dict(COMMON, chunk=chunk, plots=shp, rawdir=rawdir, outdir=outdir)
    t0 = time.time()
    print(f"\n{'='*70}\n=== RECTIFIED {day}: chunk={chunk}  ->  {outdir}\n{'='*70}",
          flush=True)

    def log(msg, _d=day):
        print(f"[{_d}] {msg}", flush=True)

    def progress(done, total, _d=day, _s={"last": 0}):
        if time.time() - _s["last"] > 60:
            _s["last"] = time.time()
            print(f"[{_d}] progress {done}/{total}", flush=True)

    try:
        pe.run_extraction(p, log, progress, threading.Event())
        print(f"[{day}] finished in {(time.time()-t0)/60:.1f} min", flush=True)
    except Exception:
        failures.append(day)
        print(f"[{day}] FAILED after {(time.time()-t0)/60:.1f} min:", flush=True)
        traceback.print_exc()

print("\nALL DONE." if not failures else f"\nFAILED: {failures}", flush=True)
sys.exit(1 if failures else 0)
