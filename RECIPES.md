# PlotExtractor recipes: raw native vs rectified plot crops

Two proven recipes for per-plot extraction from a Metashape project.
They share every setting except the crop style. Worked example values are
from the 2025 NUE+RootPhenotyping (Muresk) campaign.

## Shared settings (both recipes)

| Setting | Value | Why |
|---|---|---|
| Chunk | select **explicitly** | Never "(active)" in a multi-chunk psx — a wrong active chunk silently produces garbage (the log warns with `TINY CROP` when it happens). |
| Plot boundaries | shapefile of the **same trial area** as the chunk, one polygon per plot, `Plot_ID` field | Different flights can be different physical trials — verify before blaming the extractor. |
| Plot CRS | `auto` (read from .prj) | |
| Output GCP CRS | your projected EPSG, e.g. `7850` (GDA2020/MGA50) | Crops get a real affine geotransform + CRS — they open georeferenced in QGIS and distances measured on them are true ground distances. Fit quality is stored per file in the `PE_AFFINE_RESID` tag (expect ≤ ~2 cm). |
| Format | GeoTIFF (geotagged, lossless) | |
| Margin | `0.0` for the bare plot (`0.2` m if context wanted) | |
| Convert to sRGB | **off** for exact raw colours | The source ICC profile (e.g. Phase One Adobe RGB 1998) is embedded in the crops automatically, so colours display exactly as the source JPEGs. |
| Images per plot 1 · Limit 0 · Quality 100 · Auto-fit off · Ortho crops off · Engine Metashape API | | |

## Recipe 1 — RAW NATIVE (`Native (no resample, no black)`)

One crop per plot, cut from the single best near-nadir frame.
**Pixels are bit-exact copies of the source image** — zero resampling, zero
recompression. The frame ranking prefers near-nadir, full-resolution frames
(oblique or distant frames are demoted), so the crop holds one plot, not a
sweep of neighbours.

- Use as: the reference dataset for anything that needs untouched raw samples.
- Shape: bbox of the tilted plot — small corner triangles of surroundings are
  included; coverage is whatever the best single frame holds (often 77–100%).
- Verify: crop sizes ~10–12.5k px at the flight's native mm/px; `gdalinfo`
  shows your EPSG + geotransform; no `TINY CROP` warnings in the log.

## Recipe 2 — RECTIFIED (`Rectified rectangle`)

Each plot deskewed to an upright rectangle of its true ground dimensions at
native resolution, warped through the full camera model (lens distortion
included). What the best frame misses is gap-filled from neighbour frames
along minimum-error seams with feathered colour matching
(`smart_seam: true`, `mesh_step: 96`).

- Use as: uniform, directly comparable plot rasters — row analysis, counting
  grids, cross-plot mosaics. Whole plot in frame by design.
- Trade-off: pixels are **resampled** (bicubic) — keep the native crops as
  the bit-exact reference.
- Verify: unfilled area renders as pure black — measure black-pixel % per
  crop (`tools/verify_rectified.py` writes `rectified_qc.csv`); expect ≥ 99%
  covered away from data problems.

## Picking between them

| Need | Recipe |
|---|---|
| Untouched raw pixels (radiometry, training data provenance) | Raw native |
| Whole plot guaranteed in frame, upright, same size every day | Rectified |
| Accurate ground measurement in QGIS | Either (both are geotagged) |
| Absolute maximum sharpness | Raw native (no resampling) |

Batch/headless runs: `tools/rerun_both_days.py` (native) and
`tools/run_rectified.py` (rectified) show how to call `run_extraction()`
directly with these parameter dicts — a few minutes per 392-plot day on GPU.
