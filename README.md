# Plot Extractor

Small desktop tool that extracts **one raw image per plot** (and optionally a
deskewed **orthomosaic crop** per plot) from an Agisoft Metashape project, using a
plot-boundary shapefile or GeoJSON.

## How to run

Double-click **`PlotExtractor.bat`** (or run `python plot_extractor.py` with a
Python that has the required packages).

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
DEM, and back-projected into every aligned camera. Each camera is scored by the
**exact share of the plot it actually images** (the projected plot quad clipped to the
frame rectangle), and the plot is cropped from the best frame. Ortho crops are
exported at native resolution and deskewed so each file is exactly one plot.

**Raw crop style** (dropdown next to "Raw image crop per plot"):

- **Native (no resample, no black)** (default) — a native-pixel crop centred on the
  plot, padded symmetrically so the plot sits in the middle with equal margin on every
  side and **no black borders and no resampling** (bit-identical source pixels).
- **Auto-straighten tilted** — keeps upright plots bit-identical to the raw image, and
  only rectifies plots that are tilted >3° in their frame, then trims the black.
- **Bounding box** — the whole tilted plot plus surrounding area (over-covers).
- **Masked to plot** — bounding box with everything outside the polygon set to black;
  100% bit-identical plot pixels, but tilted plots stay tilted-in-black.
- **Rectified rectangle** — every plot warped to an upright rectangle, stitching the
  fewest frames needed to cover the whole plot (see *Accurate rectification* below).
- **Best single frame (no stitching)** — one image per plot and nothing else: the frame
  that covers the largest part of the plot, rectified upright. Whatever that frame does
  not reach stays black, and `best_cover_pct` / `best_cover_m2` in the manifest record
  exactly how much of the plot it holds. Use this when a stitch is unacceptable at any
  price. `candidates.csv` lists **every** frame that sees each plot with its coverage,
  so you can see what the alternatives were.

### Accurate rectification, and why stitches used to duplicate plants

A PhaseOne frame at 12 m covers about 4.2 m of ground, so a 4.35 m plot **cannot** fit
in one frame — a rectified plot is normally built from 2–3 frames.

Rectification maps each output pixel back to a source pixel through a grid of ground
points projected with the **full camera model** (perspective *and* lens distortion),
replayed as a piecewise mesh. Every frame of a plot is sampled on the **same ground
grid**, so the frames land on each other.

This matters a lot. The earlier version fitted a single bilinear quad through only the
plot's four projected corners. Measured on OZ Barley York that mapping is wrong by
20–95 px in the middle of a plot, and wrong *differently* for each frame, so the two
frames of a stitch disagreed by up to **37 mm on the ground** — a plant sitting on the
boundary was drawn twice. With the mesh (default 96 px cells) the two frames agree to
**0–3 px (under 1 mm)**. Coarser cells are not enough: 256 px cells still leave up to
3 mm, because the seed ridges are ~1–3 cm high and ~20 cm apart, so the mesh has to be
finer than the ridge spacing to follow the relief. 48 px buys almost nothing over 96.

Frames are also chosen to **minimise the number of seams**: the primary frame is the
one holding the largest share of the plot, and a greedy set cover then adds only the
frames still needed to reach 100%. Measured over the whole OZ Barley York trial
(696 plots): **184 plots (26%) come out of a single frame** with no seam at all, 447
need two and 65 need three. Across the stitched plots the primary frame still supplies
a **median 93%** of the pixels (min 80%, max 100%) — the rest is a narrow band along
the seam.

**Seamless stitching** (checkbox, on by default) then adds two things:

- the hand-over boundary is placed by a **minimum-error seam**: dynamic programming
  routes it along the path where the two frames look most alike, i.e. across bare soil
  instead of through a seedling. The seam may only retreat a short distance into the
  overlap (`seam_slack_px`, default 300 px ≈ 10 cm — several times a seedling) and is
  biased back towards the primary frame's own edge, so it steps around a plant without
  handing a large area to a second frame. Letting it settle mid-overlap looks tidier but
  throws away the best frame's pixels over half the plot;
- the residual brightness/colour step is removed by a **feathered offset** — measured
  along the seam column by column, smoothed hard, and faded out over 600 px.

Every output pixel still comes from **exactly one frame**. Nothing is averaged or
blended, so no plant can become a ghost; only a smooth low-frequency correction is
applied, which is why the texture stays native raw. `frames_used` and `seam_y_pct` in
the manifest say how many frames a plot needed and where the seam sits; **save seam
map (QC)** writes `plot_<id>_seam.png`, a per-pixel map of which frame each pixel came
from, so any seam can be checked directly.

Runs are deterministic — the same inputs give bit-identical crops.

Residual red/cyan fringing on standing stubble is real 3-D parallax (objects above the
DEM surface), not a warp error — it cannot be removed by rectification, which is
precisely why the seam is routed around such features.

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

## Generate plot grid (tab)

The second tab builds the plot shapefile itself, from an **orthomosaic only** — no
Metashape project and no licence needed.

1. **Load orthomosaic…** (Ctrl+O) — any georeferenced GeoTIFF. It is drawn as the
   basemap; only the visible window is rendered, so zooming a 45000 x 69000 px ortho is
   safe. Wheel zooms to the cursor, Space/middle-drag/arrows pan, **F** fits the view.
2. Click the **4 corners of the trial** in the order shown in the dropdown
   (`BL → BR → TR → TL` or `TL → TR → BR → BL`). Right-click or Ctrl+Z undoes the last
   corner, Esc clears them.
3. Set the layout: **Ranges**, **Rows per range**, **Plot width across rows (m)**,
   **Plot length along range (m)**, **Gap between rows (m)**, **Gap between ranges (m)**.
   *Stretch grid to fit the 4 clicked corners* (default on) stretches the whole grid to
   fill the quadrilateral you clicked — best for rotated or slightly irregular trials,
   with the sizes used only as plot-to-alley proportions. Switch it off to use the exact
   metre dimensions, placed from the bottom-left corner.
4. **Preview grid** (G) draws it and reports the plot count and the measured mean plot
   width × length in metres, so you can sanity-check before saving.
5. **Save shapefile** (Ctrl+S) writes the `.shp` plus a `.prj` carrying the ortho's CRS
   (a `.geojson` extension is also accepted). **Save & use** (S) does the same and drops
   the path straight into *Plot boundaries* on the Extract tab.

Attributes written per plot: `Plot_ID` (= range × 1000 + row), `B/R`, `Range`, `Row`,
`Bank` — the same scheme the extractor and the LiDAR tool expect. Plot sizes are in real
metres even when the ortho is in a geographic (degrees) CRS.

Fine-tune the result on the Extract tab with **Edit plots…** (move/resize/shift, width
and length tools) before running.

## Output files

- `RawCrops\plot_<id>.<ext>` — one raw-image crop per plot.
- `OrthoCrops\plot_<id>.<ext>` — one deskewed orthomosaic crop per plot.
- `RawCrops\plot_<id>_seam.png` — optional per-pixel "which frame" QC map.
- `manifest.csv` — one row per plot.
- `candidates.csv` — one row per (plot, candidate frame): `rank`, `camera`,
  `cover_pct`, `cover_m2`, `corners_in_frame`, `aniso` (1.0 = nadir), `mm_per_px`,
  `used`. This is the record of **all** the images available for a plot and which one
  was chosen.

**`manifest.csv` columns:** `plot_id`, `plot_area_m2` (true polygon area), `margin_m`,
`raw_image`, `raw_file`, `crop_w`, `crop_h`, `corners_in_frame`,
`raw_covered_pct` / `raw_covered_m2` (how much of the plot the delivered crop actually
holds), `cand_frames` (how many raw frames see this plot), `best_frame` /
`best_cover_pct` / `best_cover_m2` (the single best frame and how much of the plot it
covers on its own — the "best raw image per plot" record), `frames_used` (1 = seamless,
2+ = stitched), `seam_y_pct` (where the seam sits; an `x` prefix means a vertical seam),
`fit_rows` / `fit_width_m` / `fit_length_m` / `fit_area_m2`, `ortho_file`,
`ortho_covered_pct` / `ortho_covered_m2`, `source_path`, `note`. Coverage below 100%
means part of the plot fell outside every usable frame (raw) or the orthomosaic had
holes (ortho).

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
