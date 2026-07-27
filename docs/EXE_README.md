# Plot Extractor — Standalone Executable

A one-folder Windows build. **No Python install needed** on the target PC.

## How to deploy
1. Copy the whole **`PlotExtractor\`** folder (from `dist\PlotExtractor\`) to the target PC —
   anywhere (Desktop, a network share, a USB drive). Keep the folder intact.
2. Double-click **`PlotExtractor.exe`**.

That's it. The app remembers settings in a `config.json` written next to the exe.

## What runs where

| Mode | Needs | Runs on… |
|---|---|---|
| **Orthomosaic GeoTIFF** (Source = ortho) | nothing extra | **any Windows 10/11 x64** |
| **Metashape project** (Source = .psx) | a valid **Metashape licence** reachable at runtime | machines that can reach your Agisoft floating-licence server (or have a node-locked licence) |

> The Metashape engine is bundled, but **licensing is not** — Agisoft requires a live licence to run.
> On a PC with no licence, use **Orthomosaic GeoTIFF** mode (crops plots straight from an ortho
> raster; no Metashape, no licence). Everything else — auto-fit, seam-match, ortho fallback — works
> the same.

## Licence server
By default the app points at `C:\Program Files\Agisoft\Licensing\licenses`. To use a different
floating server on another PC, set the environment variable **`agisoft_LICENSE`** to the server
(`port@host`) or licence folder before launching.

## Verify a build
Open a terminal in the folder and run:
```
PlotExtractor.exe --selftest
```
It imports every engine and reports the Metashape licence status, then exits — no GUI. `SELFTEST PASS`
means the bundle is healthy.

## Notes
- First launch may take a few seconds (unpacking native libraries).
- A small console window opens alongside the GUI — it streams progress and any errors. Safe to ignore.
- The folder is large (~1.5–2 GB) because it bundles Metashape, GDAL/rasterio and PROJ.
- To rebuild: `python -m PyInstaller PlotExtractor.spec --noconfirm` from the app folder.
