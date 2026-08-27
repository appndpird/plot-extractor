# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for Plot Extractor — one-folder Windows distribution.
import os
from PyInstaller.utils.hooks import collect_all, collect_submodules

APP = SPECPATH          # build whichever copy of the project this spec sits in

datas, binaries, hidden = [], [], []
# Big native packages: collect the module, its DLLs and data files.
for pkg in ("Metashape", "rasterio", "pyproj", "shapely", "PIL"):
    try:
        d, b, h = collect_all(pkg)
        datas += d; binaries += b; hidden += h
    except Exception as e:
        print("collect_all skip", pkg, e)

hidden += collect_submodules("rasterio")
hidden += ["PIL.ImageTk", "PIL.ImageCms", "PIL._tkinter_finder",
           "rasterio.sample", "rasterio._features", "rasterio.vrt",
           "rasterio.control", "shapely.geometry", "shapefile",
           # the grid-designer tab is imported lazily inside main(), so name it
           # explicitly rather than relying on bytecode scanning
           "grid_designer",
           # the "Correct grid to rows" button imports this lazily
           "grid_refit",
           # the "Training tiles" tab imports this lazily
           "tiler"]

# app assets (logos + icon) -> _internal/assets, found via sys._MEIPASS
datas += [(os.path.join(APP, "assets"), "assets")]

a = Analysis(
    [os.path.join(APP, "plot_extractor.py")],
    pathex=[APP],
    binaries=binaries,
    datas=datas,
    hiddenimports=hidden,
    hookspath=[],
    runtime_hooks=[],
    excludes=["matplotlib", "scipy", "pandas", "PyQt5", "PySide2", "PySide6",
              "IPython", "notebook", "pytest", "sphinx"],
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, [],
    exclude_binaries=True,
    name="PlotExtractor",
    console=True,                # shows progress/errors; GUI opens alongside
    icon=os.path.join(APP, "assets", "app.ico"),
    upx=False,                   # UPX corrupts GDAL/PROJ DLLs
)
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name="PlotExtractor")
