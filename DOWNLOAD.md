# Get PlotExtractor (ready-to-run Windows app)

**Download:** https://github.com/appndpird/plot-extractor/releases/latest
-> `PlotExtractor-<date>-win64.zip` (~124 MB)

Unzip anywhere and run `PlotExtractor\PlotExtractor.exe`. No Python
installation needed - everything is bundled. A Metashape Professional
license must be reachable for `.psx` extraction; **ortho-only mode**
(cropping plots straight from an orthomosaic GeoTIFF) needs no license.

The zip includes `README.txt` (quick start + the settings that matter)
and the full `PlotExtractor_User_Guide.html`.

> Why is the exe not committed on this branch? GitHub rejects files over
> 100 MB (the bundled `Metashape.pyd` is 109 MB), so binaries are
> distributed via Releases - this file is the pointer.

## Quick-start essentials
1. Select the **chunk explicitly** - never leave "(active)" in a
   multi-chunk project.
2. Plot shapefile must be for the **same trial area** as the chunk
   (a "TINY CROP" warning in the log means they don't match).
3. Set **Output GCP CRS** to your projected EPSG (e.g. 7850) - crops
   open georeferenced in QGIS with true ground distances.
4. Style **Native** = bit-exact source pixels; **Rectified rectangle** =
   upright deskewed plot, gap-filled from neighbouring frames.
5. Leave **Convert to sRGB off** for exact raw colours - the source ICC
   profile is embedded automatically.
