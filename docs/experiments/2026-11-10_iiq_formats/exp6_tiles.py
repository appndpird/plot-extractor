"""Experiment 2026-11-10, step 6: tiles (tab 4) cut from the plot images of each format, and how
the tiles compare with the plot images they came from.

Sources: N16_geotiff (16-bit), N8_tiff (8-bit lossless), N8_jpeg (8-bit JPEG q95).
Tiles of each source (bits 'same', LZW TIFF):
  pixel tiles  1024 px, overlap 0, whole plot band     - native pixels, never resampled
  ground 1 m   1 m x 1 m, one output scale per run     - plots off that scale are resampled (Lanczos)
  row tiles    4b: rows 'middle 3', 1 m along x row spacing
Checks per tile set: tiles, plots, bit depth kept, tiles that are exact copies of the plot image
(max diff 0), resampled tiles and their scale factor.
Writes ROOT\Tiles\<source>\<set>\ and ROOT\analysis\TILES_SUMMARY.txt.
"""
import csv, os, sys, time, collections
import numpy as np
import rasterio, rasterio.windows

APP = r"E:\Software\PlotExtractor-2026-11-10"
sys.path.insert(0, APP)
import tiler

ROOT = r"F:\Ibrahim's Workspace2\Plot extraction processing\2025_Calibrationtrial_I_SIF_DPIRD\Experiment_2026-11-10_IIQ_formats"
AN = os.path.join(ROOT, "analysis")
log = lambda m: None


def check(src, ix):
    rows = list(csv.DictReader(open(ix)))
    dt, exact, bad, res, fac = collections.Counter(), 0, 0, 0, []
    for r in rows:
        with rasterio.open(os.path.join(os.path.dirname(ix), r["rel_path"])) as t_:
            a = t_.read()
            dt[t_.dtypes[0]] += 1
        if r.get("resampled", "False") == "True":
            res += 1
            try:
                fac.append(float(r["out_w"]) / float(r["w"]) if r.get("out_w") else np.nan)
            except (KeyError, ValueError):
                pass
            continue
        with rasterio.open(os.path.join(src, r["src"])) as s:
            x, y, w, h = (int(r[k]) for k in ("x0", "y0", "w", "h"))
            b = s.read(window=rasterio.windows.Window(x, y, w, h))[:3]
        if r.get("turned_90") == "True":
            b = np.rot90(b, 1, axes=(1, 2))
        exact += 1
        bad += int(not np.array_equal(a, b))
    fac = [f for f in fac if np.isfinite(f)]
    return dict(tiles=len(rows), plots=len({r["plot_stem"] for r in rows}), dtypes=dict(dt),
                exact_copies=exact - bad, mismatched=bad, resampled=res,
                scale=(f"{min(fac):.3f}-{max(fac):.3f}" if fac else ""))


def check_rows(src, out):
    ix = os.path.join(out, "row_tiles_index.csv")
    if not os.path.exists(ix):
        return {}
    rows = list(csv.DictReader(open(ix)))
    dt = collections.Counter()
    for r in rows[:200]:
        with rasterio.open(os.path.join(out, r["rel_path"])) as t_:
            dt[t_.dtypes[0]] += 1
    res = sum(1 for r in rows if r.get("resampled") == "True")
    return dict(tiles=len(rows), dtypes_first200=dict(dt), resampled=res)


L = []
for name in ("N16_geotiff", "N8_tiff", "N8_jpeg"):
    src = os.path.join(ROOT, name, "RawCrops")
    for setn, kw in (("pixel_1024", dict(tile=1024, overlap=0, tile_m=0.0)),
                     ("ground_1m", dict(tile_m=1.0)),):
        out = os.path.join(ROOT, "Tiles", name, setn)
        t = time.time()
        msgs = []
        ix = os.path.join(out, "tile_index.csv")
        if not os.path.exists(ix):
            n, p, ix = tiler.make_tiles(src, out, rows_per_plot=5, fmt="tif", min_veg_pct=0,
                                        keep_empty_frac=1.0, write_qa=False, bits="same",
                                        rows_only=True, log=msgs.append, **kw)
        if not ix:
            L.append(f"{name:12s} {setn:11s} no tiles: "
                     + (" / ".join(m for m in msgs if "[!]" in m) or "no TIFF plot images to tile"))
            print(L[-1], flush=True)
            continue
        c = check(src, ix)
        L.append(f"{name:12s} {setn:11s} {c['tiles']:4d} tiles / {c['plots']} plots  {c['dtypes']}  "
                 f"exact copies {c['exact_copies']}, mismatched {c['mismatched']}, resampled "
                 f"{c['resampled']} {('(x' + c['scale'] + ')') if c['scale'] else ''}  {round(time.time() - t)} s")
        print(L[-1], flush=True)
    out = os.path.join(ROOT, "Tiles", name, "row_tiles_middle3")
    t = time.time()
    if name != "N16_geotiff":
        continue                     # row tiles need GeoTIFF crops (ground scale)
    try:
        if not os.path.exists(os.path.join(out, "row_tiles_index.csv")):
            tiler.make_row_tiles(src, out, rows="middle 3", rows_per_plot=5, save_rows=False,
                                 row_tiles=True, tile_len_m=1.0, fmt="tif", bits="same",
                                 write_qa=False, log=log)
        c = check_rows(src, out)
        L.append(f"{name:12s} row tiles   {c.get('tiles')} tiles, {c.get('dtypes_first200')}, resampled "
                 f"{c.get('resampled')}  {round(time.time() - t)} s")
    except Exception as e:
        L.append(f"{name:12s} row tiles   failed: {e}")
    print(L[-1], flush=True)
txt = "\n".join(L)
open(os.path.join(AN, "TILES_SUMMARY.txt"), "w", encoding="utf-8").write(txt)
print(txt)
