"""Experiment 2026-11-10, step 2: the original APPN-UQ PhaseOne_Pipeline plot extraction
(Pipeline/03_backprojection.py + IIQ_tools/HR_clip.exe, unmodified) on York 2025, 36 plots.

Their functions are imported from the repo and called in the order of their main(). Only four
things differ from running `python 03_backprojection.py --config config.yaml`, all to keep the
York data untouched:
  * the 509 GB Metashape project is opened READ-ONLY (their main() opens it read-write; it never saves);
  * the photo paths (D:\BIPUL-WORKSPACE\..., gone) are re-pointed in memory to T1_Proc\JPEG;
  * the output folder is ROOT\APPN_pipeline\ instead of <project folder>\<name>_<chunk>_BackProjector;
  * HR_clip.exe writes next to the IIQ files, so it is given a folder holding COPIES of the IIQs
    their CoordinateFile picked (not T0_Raw).
scale_factor = 1.0: this project was aligned on the full-size T1_Proc JPEGs (12768 x 9564), not on
their 0.3-scale TIFFs, so their u/v are already full-size pixels.
Their shapefile needs a lowercase 'plot_id' field: a copy of the grid with plot_id = Plot_ID is written.
"""
import importlib.util, os, shutil, sys, time
import shapefile
import Metashape

REPO = r"F:\Ibrahim's Workspace2\PhaseOne_Pipeline-main\PhaseOne_Pipeline-main\Pipeline"
ROOT = r"F:\Ibrahim's Workspace2\Plot extraction processing\2025_Calibrationtrial_I_SIF_DPIRD\Experiment_2026-11-10_IIQ_formats"
RUN00 = r"G:\APPN 2025\2025_Calibrationtrial_I_SIF_DPIRD\2025-York_F\PhaseOne\2025-06-25\run00"
PSX = os.path.join(RUN00, r"T1_Proc\MetaShapeProcessing.psx")
IIQ = os.path.join(RUN00, "T0_Raw")
GRID = r"F:\Ibrahim's Workspace2\2025_Calibrationtrial_I_SIF_DPIRD\PlotShapefile\edited_grid_edited_edited_edited_edited.shp"

OUT = os.path.join(ROOT, "APPN_pipeline")
BP = os.path.join(OUT, "BackProjector")
IMG = os.path.join(OUT, "Images")          # copies of the chosen IIQs; HR_clip writes here
os.makedirs(BP, exist_ok=True)
os.makedirs(IMG, exist_ok=True)

spec = importlib.util.spec_from_file_location("appn_bp", os.path.join(REPO, "03_backprojection.py"))
bp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bp)

# grid copy with their lowercase plot_id field
src = shapefile.Reader(GRID)
names = [f[0] for f in src.fields[1:]]
shp = os.path.join(OUT, "york2025_plots_plot_id.shp")
w = shapefile.Writer(shp, shapeType=src.shapeType)
w.field("plot_id", "N", 10, 0)
for sr in src.shapeRecords():
    w.poly([sr.shape.points])
    w.record(int(sr.record["Plot_ID"]))
w.close()
shutil.copy(os.path.splitext(GRID)[0] + ".prj", os.path.splitext(shp)[0] + ".prj")

t0 = time.time()
doc = Metashape.Document()
doc.open(PSX, read_only=True)
chunk = doc.chunk
T = chunk.transform.matrix
print("chunk", chunk.label, len(chunk.cameras), "cameras", flush=True)
# the project's photo paths point to the D: drive it was processed on: re-point them (in memory
# only - the project is read-only and never saved) to the same JPEGs in T1_Proc\JPEG
JPG = os.path.join(RUN00, r"T1_Proc\JPEG")
moved = 0
for cam in chunk.cameras:
    if cam.photo is not None and not os.path.exists(cam.photo.path):
        q = os.path.join(JPG, os.path.basename(cam.photo.path))
        if os.path.exists(q):
            cam.photo.path = q
            moved += 1
print("photo paths re-pointed:", moved, flush=True)
ids = bp.import_plot_markers(chunk, shp, BP)
xyz = bp.interpolate_marker_elevation(chunk, T)
bp.backproject_markers(chunk, T, xyz, BP)
merged = bp.merge_plot_corners(BP, ids)
best = bp.select_best_camera_per_plot(chunk, BP, merged, ids)
coord = bp.write_coordinate_file(best, BP)
bp.generate_marker_check_and_clips(coord, BP)
t1 = time.time()
print(f"back-projection {t1 - t0:.0f} s; {len(best)} of {len(ids)} plots have a camera with all 4 corners",
      flush=True)

import pandas as pd
cams = sorted(set(pd.read_csv(coord)["Camera_ID"]))
for c in cams:
    d = os.path.join(IMG, c + ".IIQ")
    if not os.path.exists(d):
        shutil.copy2(os.path.join(IIQ, c + ".IIQ"), d)
# HR_clip reads config_hr.txt from the IIQ folder; their step 01 puts the template there
for f in ("config_hr.txt", "config_lr.txt"):
    shutil.copy2(os.path.join(REPO, "IIQ_tools", "template_configs", f), os.path.join(IMG, f))
t2 = time.time()
bp.run_hr_clip(os.path.join(REPO, "IIQ_tools"), IMG, coord, "1.0")
t3 = time.time()
print(f"copied {len(cams)} IIQs in {t2 - t1:.0f} s; HR_clip {t3 - t2:.0f} s", flush=True)
open(os.path.join(OUT, "times.txt"), "w").write(
    f"backprojection_s={t1 - t0:.1f}\ncopy_s={t2 - t1:.1f}\nhr_clip_s={t3 - t2:.1f}\nplots={len(ids)}\n"
    f"plots_with_camera={len(best)}\ncameras={len(cams)}\n")
