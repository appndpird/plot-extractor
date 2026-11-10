"""Experiment 2026-11-10, step 5: the plot box as RAW SENSOR values (no development at all).

The IIQ holds one 16-bit value per photosite under a Bayer colour filter (R G / G B). Every
developed image - TIFF, GeoTIFF, PNG or JPEG, 8 or 16 bit - is computed from these values
(black level, demosaicing, white balance, colour matrix, tone curve, noise reduction,
sharpening). The only output that keeps the recorded data unchanged is the sensor values
themselves. Here LibRaw (rawpy) reads them for the same box PlotExtractor developed:
  - the box of every plot as a 1-band 16-bit TIFF (Deflate, lossless), read back identical?
  - value range / levels used (how many bits the sensor data really has)
  - the box lines up with the developed crop: 2x2 Bayer cells averaged to RGB vs the Neutral
    16-bit crop averaged 2x2, correlation and shift (phase correlation, in pixels)
Writes ROOT\RawSensor\plot_<id>_cfa.tif and ROOT\analysis\RAW_SENSOR_SUMMARY.txt.
"""
import csv, os, time
import numpy as np
import rasterio
import rawpy

ROOT = r"F:\Ibrahim's Workspace2\Plot extraction processing\2025_Calibrationtrial_I_SIF_DPIRD\Experiment_2026-11-10_IIQ_formats"
OUT = os.path.join(ROOT, "RawSensor")
AN = os.path.join(ROOT, "analysis")
os.makedirs(OUT, exist_ok=True)
man = list(csv.DictReader(open(os.path.join(ROOT, "N16_geotiff", "manifest.csv"))))


def phase_shift(a, b):
    """(dy, dx) shift of b against a, sub-pixel by parabola fit, plus the peak height."""
    a = (a - a.mean()) / (a.std() + 1e-9)
    b = (b - b.mean()) / (b.std() + 1e-9)
    F = np.fft.fft2(a) * np.conj(np.fft.fft2(b))
    r = np.fft.ifft2(F / (np.abs(F) + 1e-9)).real
    iy, ix = np.unravel_index(np.argmax(r), r.shape)
    def sub(m, i, n):
        l, c, h = m[(i - 1) % n], m[i], m[(i + 1) % n]
        d = (l - 2 * c + h)
        return i + (0.5 * (l - h) / d if d else 0.0)
    dy = sub(r[:, ix], iy, r.shape[0]); dx = sub(r[iy, :], ix, r.shape[1])
    dy = dy - r.shape[0] if dy > r.shape[0] / 2 else dy
    dx = dx - r.shape[1] if dx > r.shape[1] / 2 else dx
    return dy, dx, float(r.max())


if __name__ == "__main__":
    rows, cache = [], {}
    t0 = time.time()
    for m in man:
        pid, src = m["plot_id"], m["source_path"]
        x0, y0, x1, y1 = (int(m[k]) for k in ("src_x0", "src_y0", "src_x1", "src_y1"))
        raw = rawpy.imread(src)
        vis = raw.raw_image_visible
        pat = "".join("RGBG"[int(c)] for c in raw.raw_colors_visible[y0 % 2:y0 % 2 + 2, x0 % 2:x0 % 2 + 2].ravel())
        cfa = np.ascontiguousarray(vis[y0:y1, x0:x1])
        p = os.path.join(OUT, f"plot_{pid}_cfa.tif")
        with rasterio.open(p, "w", driver="GTiff", width=cfa.shape[1], height=cfa.shape[0], count=1,
                           dtype="uint16", compress="deflate", predictor=2, tiled=True,
                           blockxsize=512, blockysize=512) as d:
            d.write(cfa, 1)
            d.update_tags(PE_CFA_PATTERN=pat, PE_SOURCE_IIQ=os.path.basename(src),
                          PE_BOX=f"{x0},{y0},{x1},{y1}", PE_BLACK=str(raw.black_level_per_channel),
                          PE_WHITE=str(raw.white_level))
        with rasterio.open(p) as d:
            back = d.read(1)
        same = bool(np.array_equal(back, cfa))
        # 2x2 cells -> RGB (pattern of this box), vs the developed crop averaged 2x2
        hh, ww = (cfa.shape[0] // 2) * 2, (cfa.shape[1] // 2) * 2
        c = cfa[:hh, :ww].astype(np.float64)
        cells = {pat[0]: c[0::2, 0::2], pat[1]: c[0::2, 1::2], pat[3]: c[1::2, 1::2]}
        g2 = (c[0::2, 1::2] + c[1::2, 0::2]) / 2 if pat[1] == "G" else (c[0::2, 0::2] + c[1::2, 1::2]) / 2
        with rasterio.open(os.path.join(ROOT, "N16_geotiff", "RawCrops", f"plot_{pid}.tif")) as d:
            dev = d.read().astype(np.float64)[:, :hh, :ww]
        dg = (dev[1, 0::2, 0::2] + dev[1, 0::2, 1::2] + dev[1, 1::2, 0::2] + dev[1, 1::2, 1::2]) / 4
        # tone curves differ (linear sensor vs developed), so compare on log values
        a, b = np.log1p(g2), np.log1p(dg)
        cy, cx = a.shape[0] // 2, a.shape[1] // 2
        s = 512
        dy, dx, pk = phase_shift(a[cy - s:cy + s, cx - s:cx + s], b[cy - s:cy + s, cx - s:cx + s])
        corr = float(np.corrcoef(a.ravel(), b.ravel())[0, 1])
        rows.append(dict(plot=pid, iiq=os.path.basename(src), w=cfa.shape[1], h=cfa.shape[0], cfa_pattern=pat,
                         min=int(cfa.min()), max=int(cfa.max()), levels=int(np.unique(cfa).size),
                         size_mb=round(os.path.getsize(p) / 1e6, 1), read_back_identical=same,
                         corr_log_green_vs_developed=round(corr, 4),
                         shift_x_px=round(2 * dx, 2), shift_y_px=round(2 * dy, 2)))
        print(pid, rows[-1], round(time.time() - t0), "s", flush=True)
        raw.close()

    with open(os.path.join(AN, "raw_sensor_per_plot.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    n = len(rows)
    mx = max(r["max"] for r in rows)
    L = [f"Raw sensor values (LibRaw {rawpy.libraw_version}) of the plot boxes, York 2025, {n} plots",
         f"  written as 1-band 16-bit TIFF (Deflate): read back identical {sum(r['read_back_identical'] for r in rows)}/{n}",
         f"  values {min(r['min'] for r in rows)}-{mx} (16-bit container, ~{np.log2(mx + 1):.1f} bits used); "
         f"levels used per plot {min(r['levels'] for r in rows)}-{max(r['levels'] for r in rows)}",
         f"  Bayer pattern of the boxes: {sorted({r['cfa_pattern'] for r in rows})}",
         f"  line-up with the developed Neutral crop (2x2 cells, log green): correlation "
         f"{min(r['corr_log_green_vs_developed'] for r in rows):.3f}-{max(r['corr_log_green_vs_developed'] for r in rows):.3f}, "
         f"shift x {min(r['shift_x_px'] for r in rows)}..{max(r['shift_x_px'] for r in rows)} px, "
         f"y {min(r['shift_y_px'] for r in rows)}..{max(r['shift_y_px'] for r in rows)} px",
         f"  size: {sum(r['size_mb'] for r in rows):.0f} MB for {n} plots "
         f"(one value per pixel: a third of a 16-bit RGB crop before compression)"]
    txt = "\n".join(L)
    open(os.path.join(AN, "RAW_SENSOR_SUMMARY.txt"), "w", encoding="utf-8").write(txt)
    print(txt)
