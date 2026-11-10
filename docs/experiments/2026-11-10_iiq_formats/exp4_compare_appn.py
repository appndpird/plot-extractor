"""Experiment 2026-11-10, step 4: the APPN-UQ PhaseOne_Pipeline clips (step 2) against
PlotExtractor's direct IIQ crops (step 1) of the same plots.

Per plot: camera chosen by each tool, their box vs ours, and - where both used the same photo -
their HR_clip pixels against PlotExtractor's crop developed with the same APPN HR settings
(HR8_geotiff) and with the Neutral settings (N16_geotiff / 257), over their box.
Writes ROOT\analysis\appn_vs_plotextractor.csv and APPN_SUMMARY.txt.
"""
import csv, math, os
import numpy as np
import rasterio

ROOT = r"F:\Ibrahim's Workspace2\Plot extraction processing\2025_Calibrationtrial_I_SIF_DPIRD\Experiment_2026-11-10_IIQ_formats"
APPN = os.path.join(ROOT, "APPN_pipeline")
OUT = os.path.join(ROOT, "analysis")
os.makedirs(OUT, exist_ok=True)


def rd(p):
    with rasterio.open(p) as d:
        info = dict(dtype=d.dtypes[0], georef=(not d.transform.is_identity) or bool(d.gcps[0]),
                    compression=d.compression.name if d.compression else "none")
        return d.read(), info


import sys
sys.path.insert(0, r"E:\Software\PlotExtractor-2026-11-10")
import iiq
CONFIG_HR = os.path.join(r"F:\Ibrahim's Workspace2\PhaseOne_Pipeline-main\PhaseOne_Pipeline-main\Pipeline",
                         r"IIQ_tools\template_configs\config_hr.txt")
dev = iiq.Developer(CONFIG_HR, os.path.join(OUT, "_work_appn"),
                    sdk_hint=r"E:\Software\PlotExtractor-2026-11-10\iiq_tools\sdk")
coords = list(csv.DictReader(open(os.path.join(APPN, r"BackProjector\CoordinateFile\CoordinateFile.csv"))))
ours = {m["plot_id"]: m for m in csv.DictReader(open(os.path.join(ROOT, "HR8_geotiff", "manifest.csv")))}
rows = []
for c in coords:
    pid, cam = c["POI"], c["Camera_ID"]
    u = np.array([float(c[f"u{i}"]) for i in range(1, 5)])
    v = np.array([float(c[f"v{i}"]) for i in range(1, 5)])
    # HR_clip box at scale 1.0 (IIQHRclip: floor of the minimum, ceil of the extent)
    x0, y0 = math.floor(u.min()), math.floor(v.min())
    clip_p = os.path.join(APPN, "Images", "output_hr", f"{cam}_{pid}.tiff")
    th, tds = rd(clip_p)
    h, w = th.shape[1:]
    quad = np.stack([u, v], 1)
    area = 0.5 * abs(np.dot(quad[:, 0], np.roll(quad[:, 1], 1)) - np.dot(quad[:, 1], np.roll(quad[:, 0], 1)))
    o = ours.get(pid)
    r = dict(plot=pid, appn_camera=cam, pe_camera=o["raw_image"] if o else "",
             same_camera=bool(o and o["raw_image"] == cam),
             appn_box=f"{x0},{y0},{w},{h}", appn_bits=tds["dtype"], appn_georef=tds["georef"], appn_compression=tds["compression"],
             appn_plot_share_pct=round(100 * area / (w * h), 1),
             appn_size_mb=round(os.path.getsize(clip_p) / 1e6, 1))
    if o and r["same_camera"]:
        ox0, oy0 = int(o["src_x0"]), int(o["src_y0"])
        ox1, oy1 = int(o["src_x1"]), int(o["src_y1"])
        r["pe_box"] = f"{ox0},{oy0},{ox1 - ox0},{oy1 - oy0}"
        # overlap of their box with ours, in each image
        ix0, iy0 = max(x0, ox0), max(y0, oy0)
        ix1, iy1 = min(x0 + w, ox1), min(y0 + h, oy1)
        r["their_box_inside_ours"] = (ix0, iy0, ix1, iy1) == (x0, y0, x0 + w, y0 + h)
        a = th[:, iy0 - y0:iy1 - y0, ix0 - x0:ix1 - x0].astype(np.int32)
        hr8, _ = rd(os.path.join(ROOT, "HR8_geotiff", "RawCrops", f"plot_{pid}.tif"))
        b = hr8[:, iy0 - oy0:iy1 - oy0, ix0 - ox0:ix1 - ox0].astype(np.int32)
        d = np.abs(a - b)
        r["vs_pe_appnhr_identical_pct"] = round(100 * float((d.max(axis=0) == 0).mean()), 2)
        r["vs_pe_appnhr_max_diff"] = int(d.max())
        r["vs_pe_appnhr_mean_abs"] = round(float(d.mean()), 4)
        n16, _ = rd(os.path.join(ROOT, "N16_geotiff", "RawCrops", f"plot_{pid}.tif"))
        nb = np.round(n16[:, iy0 - oy0:iy1 - oy0, ix0 - ox0:ix1 - ox0] / 257.0).astype(np.int32)
        dn = np.abs(a - nb)
        r["vs_pe_neutral_mean_abs"] = round(float(dn.mean()), 2)
        r["vs_pe_neutral_p99"] = float(np.percentile(dn, 99))
        # like for like: THEIR box developed by PlotExtractor's developer with THEIR config_hr.txt
        dv = dev.develop(os.path.join(APPN, "Images", cam + ".IIQ"), (x0, y0, w, h),
                         bits=8).transpose(2, 0, 1).astype(np.int32)
        dd = np.abs(th.astype(np.int32) - dv)
        r["same_box_same_config_identical_pct"] = round(100 * float((dd.max(axis=0) == 0).mean()), 2)
        r["same_box_same_config_max_diff"] = int(dd.max())
        r["same_box_same_config_mean_abs"] = round(float(dd.mean()), 4)
        # clipping: pixels with a channel at 255, their clip vs PlotExtractor Neutral 8-bit
        r["appn_clipped_pct"] = round(100 * float((th.max(axis=0) == 255).mean()), 2)
        n8, _ = rd(os.path.join(ROOT, "N8_geotiff", "RawCrops", f"plot_{pid}.tif"))
        n8b = n8[:, iy0 - oy0:iy1 - oy0, ix0 - ox0:ix1 - ox0]
        r["pe_neutral_clipped_pct"] = round(100 * float((n8b.max(axis=0) == 255).mean()), 2)
        r["appn_mean_rgb"] = "/".join(f"{v:.0f}" for v in th.reshape(3, -1).mean(1))
        r["pe_neutral_mean_rgb"] = "/".join(f"{v:.0f}" for v in n8b.reshape(3, -1).mean(1))
    rows.append(r)
dev.cleanup()

keys = []
for r in rows:
    keys += [k for k in r if k not in keys]
with open(os.path.join(OUT, "appn_vs_plotextractor.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys)
    w.writeheader()
    w.writerows(rows)

same = [r for r in rows if r["same_camera"]]
done = {r["plot"] for r in rows}
missing = sorted(set(ours) - done, key=int)
L = [f"APPN-UQ PhaseOne_Pipeline (03_backprojection.py + HR_clip.exe) on York 2025-06-25, {len(ours)} plots",
     f"  plots with a clip: {len(rows)} / {len(ours)}; dropped (no photo holds all 4 corners): {', '.join(missing)}",
     f"  PlotExtractor (direct IIQ): {len(ours)} / {len(ours)} "
     f"(the dropped plots: whole plot in a photo for "
     f"{sum(1 for p in missing if float(ours[p]['raw_covered_pct']) >= 99.95)}, "
     f"otherwise the photo holding the largest part: "
     f"{', '.join(ours[p]['raw_covered_pct'] + ' %' for p in missing if float(ours[p]['raw_covered_pct']) < 99.95)})",
     f"  same photo chosen: {len(same)} / {len(rows)}",
     f"  their clip inside PlotExtractor's box (same photo): {sum(1 for r in same if r.get('their_box_inside_ours'))} / {len(same)}",
     f"  their output: {rows[0]['appn_bits']} TIFF, compression {rows[0]['appn_compression']}, georeferenced: "
     f"{any(r['appn_georef'] for r in rows)}, plot share of the clip {min(r['appn_plot_share_pct'] for r in rows)}-"
     f"{max(r['appn_plot_share_pct'] for r in rows)} % (no mask, no margin)"]
if same:
    L.append(f"  THEIR box developed by PlotExtractor's developer with THEIR config_hr.txt: identical pixels "
             f"{np.median([r['same_box_same_config_identical_pct'] for r in same]):.2f} % (median), max diff "
             f"{max(r['same_box_same_config_max_diff'] for r in same)}/255, mean abs "
             f"{np.median([r['same_box_same_config_mean_abs'] for r in same]):.4f}/255")
    L.append(f"  their config_hr.txt (repo template now): exposure +0.8 EV, shadow / highlight recovery 0.5, "
             f"colour NR 0.5, no luminance NR, no sharpening -> mean RGB {same[0]['appn_mean_rgb']} vs Neutral "
             f"{same[0]['pe_neutral_mean_rgb']} (plot {same[0]['plot']}); clipped pixels (a channel at 255) median "
             f"{np.median([r['appn_clipped_pct'] for r in same]):.2f} % vs Neutral "
             f"{np.median([r['pe_neutral_clipped_pct'] for r in same]):.2f} %")
    L.append(f"  pixels vs PlotExtractor's 'APPN HR_clip' preset (their earlier config: NR 0.8, sharpening 3.1, no exposure change): "
             f"identical pixels "
             f"{np.median([r['vs_pe_appnhr_identical_pct'] for r in same]):.2f} % (median), max diff "
             f"{max(r['vs_pe_appnhr_max_diff'] for r in same)}/255, mean abs "
             f"{np.median([r['vs_pe_appnhr_mean_abs'] for r in same]):.4f}/255")
    L.append(f"  pixels vs PlotExtractor Neutral (no noise reduction / sharpening): mean abs "
             f"{np.median([r['vs_pe_neutral_mean_abs'] for r in same]):.2f}/255, 99th pct "
             f"{np.median([r['vs_pe_neutral_p99'] for r in same]):.0f}/255")
txt = "\n".join(L)
open(os.path.join(OUT, "APPN_SUMMARY.txt"), "w", encoding="utf-8").write(txt)
print(txt)
