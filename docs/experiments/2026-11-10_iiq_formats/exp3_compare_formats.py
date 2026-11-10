"""Experiment 2026-11-10, step 3: what each output format keeps of the IIQ development.

For every plot of the runs of step 1 (same frame, same box in every run):
  A  16-bit containers   N16_geotiff vs N16_tiff vs N16_tiff_none vs N16_png       identical?
  B  8-bit containers    N8_geotiff vs N8_tiff vs N8_png                          identical?
  C  16 -> 8 bit         N8 vs round(N16 / 257); tonal levels used in the plot
  D  JPEG                N8_jpeg vs N8_tiff: PSNR, mean / max error, pixels changed
  E  development         HR16 vs N16 (APPN HR noise reduction + sharpening vs none)
  F  fidelity            6 plots developed again from the IIQ (Neutral, 16-bit) = saved file?
  G  georeference        GeoTIFF affine fit residuals (PE_AFFINE_RESID), CRS, pixel size
  H  compression         one 16-bit and one 8-bit crop written with LZW / Deflate / ZSTD /
                         none (horizontal predictor) - sizes, read back identical?
Writes formats_per_plot.csv and FORMATS_SUMMARY.txt into ROOT\analysis\.
"""
import csv, json, os, re, sys, time
import numpy as np
import rasterio
from PIL import Image

APP = r"E:\Software\PlotExtractor-2026-11-10"
sys.path.insert(0, APP)
ROOT = r"F:\Ibrahim's Workspace2\Plot extraction processing\2025_Calibrationtrial_I_SIF_DPIRD\Experiment_2026-11-10_IIQ_formats"
OUT = os.path.join(ROOT, "analysis")
os.makedirs(OUT, exist_ok=True)


def rd(path):
    if path.lower().endswith(".jpg"):
        return np.asarray(Image.open(path).convert("RGB")).transpose(2, 0, 1)
    with rasterio.open(path) as d:
        return d.read()


def crop_file(run, pid):
    d = os.path.join(ROOT, run, "RawCrops")
    for e in (".tif", ".png", ".jpg"):
        p = os.path.join(d, f"plot_{pid}{e}")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(run, pid)


man = list(csv.DictReader(open(os.path.join(ROOT, "N16_geotiff", "manifest.csv"))))
pids = [m["plot_id"] for m in man]
rows = []
t0 = time.time()
for m in man:
    pid = m["plot_id"]
    r = dict(plot=pid, frame=m["raw_image"], w=m["crop_w"], h=m["crop_h"])
    n16 = rd(crop_file("N16_geotiff", pid))
    r["A_16bit_tiff_same"] = bool(np.array_equal(n16, rd(crop_file("N16_tiff", pid))))
    r["A_16bit_tiff_none_same"] = bool(np.array_equal(n16, rd(crop_file("N16_tiff_none", pid))))
    r["A_16bit_png_same"] = bool(np.array_equal(n16, rd(crop_file("N16_png", pid))))
    n8 = rd(crop_file("N8_geotiff", pid))
    r["B_8bit_tiff_same"] = bool(np.array_equal(n8, rd(crop_file("N8_tiff", pid))))
    r["B_8bit_png_same"] = bool(np.array_equal(n8, rd(crop_file("N8_png", pid))))
    # C: 16 -> 8 bit
    q = np.round(n16.astype(np.float64) / 257.0).astype(np.int32)
    d = np.abs(n8.astype(np.int32) - q)
    r["C_8_vs_16div257_maxdiff"] = int(d.max())
    r["C_8_vs_16div257_pct_diff"] = round(100.0 * float((d > 0).mean()), 4)
    err = n16.astype(np.float64) - n8.astype(np.float64) * 257.0       # 16-bit units
    r["C_quant_rms_16bit_units"] = round(float(np.sqrt((err ** 2).mean())), 1)
    for k, c in enumerate("RGB"):
        r[f"C_levels16_{c}"] = int(np.unique(n16[k]).size)
        r[f"C_levels8_{c}"] = int(np.unique(n8[k]).size)
    # D: JPEG vs lossless 8-bit
    jp = np.asarray(Image.open(crop_file("N8_jpeg", pid)).convert("RGB")).transpose(2, 0, 1)
    dj = jp.astype(np.float64) - n8.astype(np.float64)
    mse = float((dj ** 2).mean())
    r["D_jpeg_psnr_db"] = round(10 * np.log10(255.0 ** 2 / mse), 2) if mse > 0 else 99.0
    r["D_jpeg_mean_abs"] = round(float(np.abs(dj).mean()), 3)
    r["D_jpeg_max_abs"] = int(np.abs(dj).max())
    r["D_jpeg_pct_changed"] = round(100.0 * float((np.abs(dj).max(axis=0) > 0).mean()), 2)
    # E: APPN HR vs neutral development (8-bit units)
    hr16 = rd(crop_file("HR16_geotiff", pid))
    de = (hr16.astype(np.float64) - n16.astype(np.float64)) / 257.0
    r["E_hr_vs_neutral_mean_abs_8bit"] = round(float(np.abs(de).mean()), 3)
    r["E_hr_vs_neutral_p99_abs_8bit"] = round(float(np.percentile(np.abs(de), 99)), 2)
    g16 = n16[1].astype(np.float64) / 257.0
    gh = hr16[1].astype(np.float64) / 257.0
    hp = lambda a: a[1:-1, 1:-1] * 4 - a[:-2, 1:-1] - a[2:, 1:-1] - a[1:-1, :-2] - a[1:-1, 2:]
    r["E_neutral_highpass_sd"] = round(float(hp(g16).std()), 3)
    r["E_hr_highpass_sd"] = round(float(hp(gh).std()), 3)
    hr8 = rd(crop_file("HR8_geotiff", pid))
    dh = np.abs(hr8.astype(np.int32) - np.round(hr16 / 257.0).astype(np.int32))
    r["E_hr8_vs_hr16div257_maxdiff"] = int(dh.max())
    r["E_hr8_vs_hr16div257_pct_diff"] = round(100.0 * float((dh > 0).mean()), 4)
    # G: georeference
    with rasterio.open(crop_file("N16_geotiff", pid)) as ds:
        tags = ds.tags()
        r["G_crs"] = ds.crs.to_string() if ds.crs else ""
        r["G_px_mm"] = round(float(np.hypot(ds.transform.a, ds.transform.d)) * 1000, 3)
        mm = re.findall(r"mean=([\d.]+)mm max=([\d.]+)mm n=(\d+)", tags.get("PE_AFFINE_RESID", ""))
        if mm:
            r["G_affine_resid_mean_mm"], r["G_affine_resid_max_mm"] = float(mm[0][0]), float(mm[0][1])
            r["G_gcps"] = int(mm[0][2])
        r["G_dtype"] = ds.dtypes[0]
        r["G_compress"] = ds.compression.name if ds.compression else "none"
    rows.append(r)
    print(pid, "ok", round(time.time() - t0), "s", flush=True)

with open(os.path.join(OUT, "formats_per_plot.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

# F: fresh development of 6 plots (Neutral, 16-bit) vs the saved files
import iiq
dev = iiq.Developer(os.path.join(APP, r"iiq_tools\presets\neutral.json"),
                    os.path.join(OUT, "_work"), sdk_hint=os.path.join(APP, r"iiq_tools\sdk"))
fres = []
for m in man[:6]:
    x0, y0, x1, y1 = (int(m[k]) for k in ("src_x0", "src_y0", "src_x1", "src_y1"))
    a = dev.develop(m["source_path"], (x0, y0, x1 - x0, y1 - y0), bits=16).transpose(2, 0, 1)
    b = rd(crop_file("N16_geotiff", m["plot_id"]))
    fres.append((m["plot_id"], bool(np.array_equal(a, b)), int(np.abs(a.astype(np.int32) - b).max())))
dev.cleanup()

# H: compression of the same pixels
comp = []
for name, run in (("16-bit", "N16_tiff_none"), ("8-bit", "N8_tiff")):
    src = crop_file(run, pids[0])
    a = rd(src)
    for c, kw in (("none", {}), ("lzw+pred2", dict(compress="lzw", predictor=2)),
                  ("lzw", dict(compress="lzw")),
                  ("deflate+pred2", dict(compress="deflate", predictor=2, zlevel=6)),
                  ("zstd+pred2", dict(compress="zstd", predictor=2, zstd_level=9))):
        p = os.path.join(OUT, f"_comp_{name}_{c}.tif")
        t = time.time()
        with rasterio.open(p, "w", driver="GTiff", width=a.shape[2], height=a.shape[1], count=3,
                           dtype=a.dtype.name, photometric="RGB", tiled=True, blockxsize=512,
                           blockysize=512, **kw) as d:
            d.write(a)
        tw = time.time() - t
        same = bool(np.array_equal(rd(p), a))
        comp.append((name, c, os.path.getsize(p) / 1e6, round(tw, 2), same))
        os.remove(p)

sizes = {}
for run in ("N16_geotiff", "N16_tiff", "N16_tiff_none", "N16_png", "N8_geotiff", "N8_tiff",
            "N8_png", "N8_jpeg", "HR16_geotiff", "HR8_geotiff"):
    d = os.path.join(ROOT, run, "RawCrops")
    sizes[run] = sum(os.path.getsize(os.path.join(d, f)) for f in os.listdir(d)) / 1e6
times = json.load(open(os.path.join(ROOT, "run_times.json")))

def allt(k):
    return sum(1 for r in rows if r[k])

def med(k):
    return float(np.median([r[k] for r in rows if r.get(k) not in (None, "")]))

def rng(k):
    v = [r[k] for r in rows if r.get(k) not in (None, "")]
    return f"{min(v)}-{max(v)}"

n = len(rows)
L = []
L.append(f"IIQ formats experiment - York 2025-06-25, {n} plots, direct IIQ extraction (bounding box, margin 0.10 m)")
L.append("")
L.append("A  16-bit containers (Neutral development)")
L.append(f"   GeoTIFF = TIFF (LZW): {allt('A_16bit_tiff_same')}/{n} plots pixel-identical")
L.append(f"   GeoTIFF = TIFF (uncompressed): {allt('A_16bit_tiff_none_same')}/{n}")
L.append(f"   GeoTIFF = PNG: {allt('A_16bit_png_same')}/{n}")
L.append("B  8-bit containers")
L.append(f"   GeoTIFF = TIFF: {allt('B_8bit_tiff_same')}/{n};  GeoTIFF = PNG: {allt('B_8bit_png_same')}/{n}")
L.append("C  16 -> 8 bit")
L.append(f"   8-bit file = round(16-bit / 257): max diff {max(r['C_8_vs_16div257_maxdiff'] for r in rows)}, "
         f"pixels differing (median) {med('C_8_vs_16div257_pct_diff'):.4f} %")
L.append(f"   tonal levels used in a plot, green: 16-bit {rng('C_levels16_G')} vs 8-bit {rng('C_levels8_G')}")
L.append(f"   rounding error of 8-bit: RMS {med('C_quant_rms_16bit_units'):.0f} of 65535 (= 1/sqrt(12) of one 8-bit step)")
L.append("D  JPEG (quality 95) vs lossless 8-bit")
L.append(f"   PSNR median {med('D_jpeg_psnr_db'):.1f} dB (range {rng('D_jpeg_psnr_db')}), mean abs error "
         f"{med('D_jpeg_mean_abs'):.2f}/255, max {max(r['D_jpeg_max_abs'] for r in rows)}/255, "
         f"pixels changed {med('D_jpeg_pct_changed'):.1f} %")
L.append("E  development settings (APPN HR vs Neutral, same box)")
L.append(f"   mean abs difference {med('E_hr_vs_neutral_mean_abs_8bit'):.2f}/255, 99th pct {med('E_hr_vs_neutral_p99_abs_8bit'):.1f}/255")
L.append(f"   high-pass SD (green, 8-bit units): neutral {med('E_neutral_highpass_sd'):.2f} vs APPN HR {med('E_hr_highpass_sd'):.2f}")
L.append(f"   HR 8-bit = round(HR 16-bit / 257): max diff {max(r['E_hr8_vs_hr16div257_maxdiff'] for r in rows)}, "
         f"pixels differing {med('E_hr8_vs_hr16div257_pct_diff'):.4f} % (noise reduction not bit-repeatable)")
L.append("F  fidelity: the box developed again from the IIQ (Neutral, 16-bit) vs the saved GeoTIFF")
for pid, same, mx in fres:
    L.append(f"   plot {pid}: {'identical' if same else 'DIFFERENT'} (max diff {mx})")
L.append("G  georeference of the GeoTIFFs")
L.append(f"   CRS {rows[0].get('G_crs')}, pixel {rng('G_px_mm')} mm, affine fit through {rows[0].get('G_gcps')} GCPs: "
         f"residual mean {med('G_affine_resid_mean_mm'):.2f} mm (median of plots), max {max(r.get('G_affine_resid_max_mm', 0) for r in rows):.2f} mm")
L.append("H  compression of the same pixels (plot %s, 512 px tiles)" % pids[0])
for name, c, mb, tw, same in comp:
    L.append(f"   {name:6s} {c:14s} {mb:8.1f} MB  write {tw:5.2f} s  read back identical: {same}")
L.append("")
L.append("Files of the 36 plots (MB) and run times (s, incl. opening the project)")
for k, v in sizes.items():
    L.append(f"   {k:14s} {v:8.0f} MB   {times.get(k, '')} s")
txt = "\n".join(L)
open(os.path.join(OUT, "FORMATS_SUMMARY.txt"), "w", encoding="utf-8").write(txt)
print(txt)
