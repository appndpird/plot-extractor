# Plot Extractor

Desktop tool that extracts **one raw image per plot** (and optionally a deskewed
**orthomosaic crop** per plot) from an Agisoft Metashape project, using a
plot-boundary shapefile or GeoJSON. Crops are snapped to the planted rows,
radiometrically seam-matched, and each file is exactly one plot — built for
plant-emergence phenotyping.

*DPIRD × APPN · NUE + Root Phenotyping trials (Muresk / York).*

## Documentation & downloads

- **User guide, presentation & demo script:** see [`docs/`](docs/)
  (`PlotExtractor_User_Guide.html`, `PlotExtractor_Presentation.html` / `.pptx`,
  `DEMO_SCRIPT.md`).
- **Standalone Windows executable** (no Python needed) and **demo videos** are
  attached to the latest [**Release**](../../releases) — see
  [`docs/EXE_README.md`](docs/EXE_README.md) for deployment.

> The exe bundles the Metashape engine but not a licence: **Orthomosaic-GeoTIFF
> mode** runs anywhere; **Metashape-project mode** needs a live Agisoft licence.

## How to run

Double-click **`PlotExtractor.bat`** (or run `python plot_extractor.py` with a
Python that has the required packages). Build the exe with
`python -m PyInstaller PlotExtractor.spec --noconfirm`.

## Inputs

| Field | Required | Notes |
|---|---|---|
| Source | — | **Metashape project** (raw + ortho crops, needs a licence) or **Orthomosaic GeoTIFF (no Metashape)** — crops plots straight out of an ortho raster with no Metashape and no licence (ortho crops only) |
| Metashape project (.psx) | if Source = Metashape | Must be aligned + georeferenced; DEM needed for best accuracy, orthomosaic needed for ortho crops |
| Orthomosaic GeoTIFF | if Source = Ortho GeoTIFF | Any georeferenced ortho `.tif`; plot rings are auto-reprojected into the raster CRS |
| Plot boundaries (.shp / .geojson) | yes | **Any CRS is accepted** — read from the `.prj` / GeoJSON `crs` and converted automatically |
| Plot ID attribute | auto | Dropdown fills from the file's attributes; auto-picks `Plot_ID`-like fields, else numbers plots 1..N |
| Plot CRS | auto | Only set this (an EPSG code, e.g. `28350`) if the file has no CRS info |
| Raw images folder | optional | If set, images here are **preferred** over the .psx's own photo paths (use it to point at your original/best-quality JPEGs); also serves as fallback when project paths are stale. `manifest.csv` records the exact source file used per plot |
| Output folder | yes | Gets `RawCrops\`, `OrthoCrops\` and `manifest.csv` |
| Margin (m) | 0.2 | Extra area captured around each plot |
| Images per plot | 1 | Extra views saved as `plot_<id>_2.png`, … |
| Limit | 0 | Set e.g. 5 for a quick test run of the first plots |
| Format | GeoTIFF | **GeoTIFF (geotagged, lossless)**: raw crops carry WGS84 GCPs + the source image's EXIF (GDAL-readable via `gdalinfo`); ortho crops carry a true geotransform with everything outside the plot polygon set to nodata. PNG / TIFF are lossless without geotags; JPEG is smaller but lossy |
| Convert raw crops to sRGB | off | Converts each raw crop from its embedded ICC profile (e.g. **Adobe RGB**) to **sRGB**. Tick this when the raw source is Adobe RGB — e.g. Capture One / `T1_Proc` JPEGs — so downstream tools that assume sRGB read the greens correctly. Leave off if the source is already sRGB (e.g. iX-Capture TIFs) |
| CPU cores | 0 (auto) | Number of cores used for the raw image work (decode / crop / convert / save), which runs in parallel. `0` = auto (min of 8 and the machine's cores). Raise it for faster runs; each active core holds a full decoded frame (~370 MB) in memory |

## How it works

For each plot polygon: vertices are reprojected to the project CRS, dropped onto the
DEM, and back-projected into every aligned camera. The camera with the most plot
corners in-frame (then the least overshoot) wins, and the plot is cropped from that
raw photo. Ortho crops are exported at native resolution and deskewed so each file is
exactly one plot.

**Raw crop style** (dropdown next to "Raw image crop per plot"):

- **Native (no resample, no black)** (default) — a native-pixel crop centred on the
  plot, padded symmetrically so the plot sits in the middle with equal margin on every
  side and **no black borders and no resampling** (bit-identical source pixels).
- **Auto-straighten tilted** — keeps upright plots bit-identical to the raw image, and
  only rectifies plots that are tilted >3° in their frame, then trims the black.
- **Bounding box** — the whole tilted plot plus surrounding area (over-covers).
- **Masked to plot** — bounding box with everything outside the polygon set to black;
  100% bit-identical plot pixels, but tilted plots stay tilted-in-black.
- **Rectified rectangle** — every plot warped to an upright rectangle (resampled).

**Performance:** the Metashape geometry is resolved first (one pass), then all raw
crops are decoded, cropped, colour-converted and saved **in parallel across CPU cores**
(set by *CPU cores*). If a GPU is present it is enabled for Metashape's own compute
(e.g. ortho raster export); the bulk of a plot-extraction run is image decode/compress,
which is what the multi-core image stage accelerates.

GeoTIFF output geotags raw crops (GCPs, recomputed for rectified/straightened grids)
and ortho crops (a true geotransform that places them correctly in QGIS). Tick both
"Raw image crop" and "Orthomosaic crop" to generate both per plot.

**Edit plots…** opens an interactive editor: grow/shrink plots by metres (one or all,
+/- keys), drag corners to reshape (single selection), **drag the blue handle to
resize every selected plot at once** (select 2+ plots — a dashed blue box with corner
handles appears; dragging a handle out/in scales all selected plots about their own
centres, so they keep their positions and only change size), move one or many plots by
dragging, load an ortho GeoTIFF as a basemap to work over the imagery, then
"Save & use" writes an edited shapefile the extractor picks up.

**`manifest.csv` columns:** `plot_id`, `plot_area_m2` (true polygon area), `margin_m`,
`raw_image`, `raw_file`, `crop_w`, `crop_h`, `corners_in_frame`,
`raw_covered_pct` / `raw_covered_m2` (how much of the plot the chosen raw frame
actually images), `ortho_file`, `ortho_covered_pct` / `ortho_covered_m2`,
`source_path`, `note`. Coverage below 100% means part of the plot fell outside the
frame (raw) or the orthomosaic had holes (ortho).

## Requirements (already installed on this PC, Python 3.11)

```
pip install <Metashape wheel from agisoft.com/downloads> pillow pyshp
```

A valid Metashape license is required — the app uses the floating license server
configured in `C:\Program Files\Agisoft\Licensing\licenses`.

## Troubleshooting

- **"No valid Metashape license"** — license server unreachable (VPN/network).
- **"plots are ~X km away from the project centre"** — wrong CRS or wrong .psx.
- **"raw file missing"** — set the *Raw images folder* to where the photos live
  (e.g. the project's `JPEG` folder).
