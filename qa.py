"""qa.py - quality checks and a viewer for an extraction run and its tiles.

Checks (run_checks), written to <tiles folder>/QA/ (or <run>/QA/ without tiles):

  crops     raw crops (any style except 'rectify') are compared pixel by pixel with
            the box of their source frame recorded in manifest.csv (src_x0..src_y1),
            decoded with the same decoder (Pillow / libjpeg). Every pixel inside the
            plot must be identical; the mask only sets the outside to 0. Rectified
            crops are re-rendered on a ground grid by design, so they are listed with
            their frames, not compared.
  tiles     every tile in tile_index.csv must be the exact pixels of its window in the
            crop (unless the tile was resampled to a fixed ground scale), and in
            middle-rows mode it must lie inside the usable plot length.
  rows      rows_report.csv: plots where rows were not found, where the number of rows
            differs from the run's usual number, or where the row spacing is more than
            8 % away from the run's median.
  row imgs  every saved row image (Rows/rows_index.csv) must equal its window in the
            crop wherever it is not masked.

The tab (build_qa_tab) lists the plots with their flags and shows, per plot, the row
detection with the middle band and the tiles, the tiles themselves, the saved rows, and
the QA report.
"""
from __future__ import annotations

import csv
import os
from collections import Counter

import numpy as np

QA_DIR = "QA"


def _read_csv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _write_csv(rows, path):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys or ["none"])
        w.writeheader()
        w.writerows(rows)


def _crops_dir(run_dir):
    d = os.path.join(run_dir, "RawCrops")
    return d if os.path.isdir(d) else run_dir


def find_run_for_tiles(tiles_dir):
    """The extraction folder a tiles folder was cut from (from tile_index.csv)."""
    idx = _read_csv(os.path.join(tiles_dir, "tile_index.csv"))
    par = os.path.dirname(os.path.abspath(tiles_dir))
    if os.path.isdir(os.path.join(par, "RawCrops")):
        return par
    return par if idx else ""


# ---------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------
def _iiq_developer(run_dir):
    """Developer with the settings and SDK an IIQ extraction run used, and whether
    those settings are neutral (no noise reduction -> bit-repeatable)."""
    import json
    import iiq
    info = {}
    try:
        with open(os.path.join(run_dir, "IIQ_run.json")) as f:
            info = json.load(f)
    except Exception:
        pass
    settings = os.path.join(run_dir, "IIQ_settings_used.json")
    dev = iiq.Developer(settings, os.path.join(run_dir, "QA", "_iiq_work"),
                        info.get("sdk_root", ""), parallel=1)
    try:
        with open(settings) as f:
            cfg = json.load(f)
        neutral = not cfg.get("noiseReductionEnabled") and not cfg.get("sharpeningEnabled")
    except Exception:
        neutral = False
    return dev, neutral


def check_crops(run_dir, log=print, cancel=None, progress=None):
    import rasterio
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    man = _read_csv(os.path.join(run_dir, "manifest.csv"))
    crops = _crops_dir(run_dir)
    res = []
    cache = {}
    for i, r in enumerate(man):
        if cancel is not None and cancel.is_set():
            break
        fn = r.get("raw_file", "")
        p = os.path.join(crops, fn)
        rec = dict(plot=r.get("plot_id", ""), file=fn, style=r.get("raw_style", ""),
                   source=os.path.basename(r.get("source_path", "")),
                   frames_used=r.get("frames_used", ""),
                   covered_pct=r.get("raw_covered_pct", ""))
        if not fn or not os.path.exists(p):
            rec.update(result="MISSING", note="crop file not found")
        elif r.get("raw_style", "") == "rectify":
            rec.update(result="RESAMPLED", note="rectified: re-rendered on a ground grid "
                                                "(not a pixel copy by design)")
        elif not all(r.get(k, "") != "" for k in ("src_x0", "src_y0", "src_x1", "src_y1")):
            rec.update(result="NO SOURCE BOX", note="manifest has no src_x0..src_y1")
        elif not os.path.exists(r.get("source_path", "")):
            rec.update(result="NO SOURCE", note="source frame not reachable")
        elif r.get("source_path", "").lower().endswith(".iiq"):
            # developed from an IIQ: develop the same box again with the run's settings
            # and SDK and compare
            try:
                if "iiq" not in cache:
                    cache["iiq"] = _iiq_developer(run_dir)
                dev, neutral = cache["iiq"]
                x0, y0, x1, y1 = (int(float(r[k])) for k in ("src_x0", "src_y0", "src_x1", "src_y1"))
                with rasterio.open(p) as s:
                    crop = s.read([1, 2, 3]).transpose(1, 2, 0)
                bits = 16 if crop.dtype == np.uint16 else 8
                ref = dev.develop(r["source_path"], (x0, y0, x1 - x0, y1 - y0), bits=bits)
                top = 65535 if bits == 16 else 255
                tol = 2 * (257 if bits == 16 else 1)      # 2/255 at either depth
                rec["bits"] = bits
                if ref.shape != crop.shape:
                    rec.update(result="FAIL", note=f"size {crop.shape} vs IIQ box {ref.shape}")
                else:
                    inside = crop.max(axis=2) > 0
                    diff = np.abs(crop.astype(np.int32) - ref.astype(np.int32)).max(axis=2)
                    n_in = int(inside.sum())
                    n_same = int((diff[inside] == 0).sum())
                    mx = int(diff[inside].max()) if n_in else 0
                    if n_same == n_in:
                        res_, note = "IDENTICAL", "same as a fresh IIQ development"
                    elif mx <= tol and not neutral:
                        res_, note = "MATCH", (f"within {mx}/{top} of a fresh development on "
                                               f"{100.0 * (n_in - n_same) / n_in:.1f}% of pixels "
                                               "(the SDK's noise reduction is not "
                                               "bit-repeatable)")
                    else:
                        res_, note = "FAIL", f"{n_in - n_same} pixels differ (max {mx})"
                    rec.update(inside_px=n_in, identical_px=n_same, max_diff=mx,
                               result=res_, note=note)
            except Exception as e:
                rec.update(result="ERROR", note=f"IIQ re-development: {e}")
        else:
            try:
                sp = r["source_path"]
                if sp not in cache:
                    cache.clear()
                    with Image.open(sp) as im:
                        cache[sp] = np.asarray(im.convert("RGB"))
                src = cache[sp]
                x0, y0, x1, y1 = (int(float(r[k])) for k in ("src_x0", "src_y0", "src_x1", "src_y1"))
                ref = src[y0:y1, x0:x1]
                with rasterio.open(p) as s:
                    crop = s.read([1, 2, 3]).transpose(1, 2, 0)
                if ref.shape != crop.shape:
                    rec.update(result="FAIL", note=f"size {crop.shape} vs source box {ref.shape}")
                else:
                    inside = crop.max(axis=2) > 0
                    diff = np.abs(crop.astype(np.int16) - ref.astype(np.int16)).max(axis=2)
                    n_in = int(inside.sum())
                    n_same = int((diff[inside] == 0).sum())
                    rec.update(inside_px=n_in, identical_px=n_same,
                               max_diff=int(diff[inside].max()) if n_in else 0,
                               result="IDENTICAL" if n_same == n_in else "FAIL",
                               note="" if n_same == n_in else f"{n_in - n_same} pixels differ")
            except Exception as e:
                rec.update(result="ERROR", note=str(e))
        res.append(rec)
        if progress:
            progress(i + 1, len(man))
    return res


def check_tiles(tiles_dir, crops_dir, log=print, cancel=None, progress=None):
    import rasterio
    import tiler
    from rasterio.windows import Window
    from PIL import Image
    idx = _read_csv(os.path.join(tiles_dir, "tile_index.csv"))
    rep = {r["plot_stem"]: r for r in _read_csv(os.path.join(tiles_dir, "rows_report.csv"))}
    res = []
    open_src = {}
    for i, t in enumerate(idx):
        if cancel is not None and cancel.is_set():
            break
        rec = dict(tile=t["tile"], plot=t["plot_stem"], x0=t["x0"], y0=t["y0"],
                   w=t["w"], h=t["h"], veg_pct=t.get("veg_pct", ""))
        try:
            tp = os.path.join(tiles_dir, t.get("rel_path") or t["tile"])
            sp = os.path.join(crops_dir, t["src"])
            x0, y0, w, h = int(t["x0"]), int(t["y0"]), int(t["w"]), int(t["h"])
            tile = tiler.read_rgb(tp)
            if t["src"] not in open_src:
                for s_ in open_src.values():
                    s_.close()
                open_src.clear()
                open_src[t["src"]] = rasterio.open(sp)
            s = open_src[t["src"]]
            if tile.shape[:2] != (h, w):
                rec.update(pixels="RESAMPLED", note=f"tile {tile.shape[1]}x{tile.shape[0]} from "
                                                    f"a {w}x{h} window (fixed ground scale)")
            else:
                ref = s.read([1, 2, 3], window=Window(x0, y0, w, h)).transpose(1, 2, 0)
                conv = ref.dtype != tile.dtype
                if conv:            # tile written at another bit depth on purpose
                    ref = tiler.convert_bits(ref, 16 if tile.dtype == np.uint16 else 8)
                same = bool(np.array_equal(ref, tile))
                rec.update(pixels=("EXACT CONVERSION" if conv else "IDENTICAL") if same
                           else "FAIL",
                           note="" if same else "tile differs from its window in the crop")
                if conv and same:
                    rec["note"] = (f"{8 * np.dtype(s.dtypes[0]).itemsize}-bit window converted "
                                   f"exactly to {8 * tile.dtype.itemsize}-bit")
            rec["georef"] = _tile_georef_check(tp, s, x0, y0, w, h)
            rr = rep.get(t["plot_stem"])
            if rr and rr.get("usable_lo", "") != "":
                ax = int(t.get("across_axis") or 0)
                a0 = x0 if ax == 0 else y0
                a1 = a0 + (w if ax == 0 else h)
                ul, uh = float(rr["usable_lo"]), float(rr["usable_hi"])
                inside = a0 >= ul - 1 and a1 <= uh + 1
                rec.update(along_from=a0, along_to=a1, usable_lo=int(ul), usable_hi=int(uh),
                           in_usable_length="yes" if inside else "NO")
        except Exception as e:
            rec.update(pixels="ERROR", note=str(e))
        res.append(rec)
        if progress:
            progress(i + 1, len(idx))
    for s_ in open_src.values():
        s_.close()
    return res


def _tile_georef_check(tile_path, s, x0, y0, w, h):
    """'' (tile has no georeference), 'OK' or 'WRONG ...': the tile's transform must be
    the crop's transform moved to the tile window (scaled if resampled)."""
    if not tile_path.lower().endswith((".tif", ".tiff")):
        return ""
    import warnings
    import rasterio
    from rasterio.transform import Affine
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with rasterio.open(tile_path) as t:
            if t.crs is None and not t.gcps[0]:
                return ""
            if t.transform.is_identity:
                return "GCPs"
            if s.crs is None or s.transform.is_identity:
                return "tile georeferenced, crop not"
            if t.crs != s.crs:
                return f"WRONG CRS {t.crs} vs {s.crs}"
            exp = s.transform * Affine.translation(x0, y0) * Affine.scale(
                w / float(t.width), h / float(t.height))
            gx, gy = exp * (t.width / 2.0, t.height / 2.0)
            hx, hy = t.transform * (t.width / 2.0, t.height / 2.0)
            off = float(np.hypot(gx - hx, gy - hy))
            return "OK" if off < 1e-6 + abs(exp.a) * 0.01 else f"WRONG by {off:.4f} map units"


def check_rows(tiles_dir):
    rep = _read_csv(os.path.join(tiles_dir, "rows_report.csv"))
    ok = [r for r in rep if str(r.get("ok")) == "True"]
    nrows = Counter(int(r["rows_found"]) for r in ok if r.get("rows_found", "") != "")
    usual = nrows.most_common(1)[0][0] if nrows else 0
    pitches = [float(r["pitch_mm"] or r["pitch_px"]) for r in ok
               if (r.get("pitch_mm") or r.get("pitch_px"))]
    med = float(np.median(pitches)) if pitches else 0.0
    res = []
    for r in rep:
        flags = []
        if str(r.get("ok")) != "True":
            flags.append("rows not found: " + (r.get("reason") or "?"))
        else:
            if usual and int(r.get("rows_found") or 0) != usual:
                flags.append(f"{r['rows_found']} rows (usual {usual})")
            pv = r.get("pitch_mm") or r.get("pitch_px")
            if med and pv and abs(float(pv) / med - 1) > 0.08:
                flags.append(f"row spacing {pv} vs median {med:.1f}")
            if r.get("reason"):
                flags.append(r["reason"])
            if str(r.get("tiles_drawn", "")) == "0":
                flags.append("no tiles")
        res.append(dict(plot=r["plot_stem"], rows_found=r.get("rows_found", ""),
                        pitch_mm=r.get("pitch_mm", ""), angle_deg=r.get("angle_deg", ""),
                        middle_rows=r.get("middle_rows", ""), tiles=r.get("tiles_drawn", ""),
                        veg_mode=r.get("veg_mode", ""),
                        flag="; ".join(flags)))
    return res, usual, med


def check_row_images(tiles_dir, crops_dir, cancel=None, progress=None):
    import rasterio
    import tiler
    from rasterio.windows import Window
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    rdir = os.path.join(tiles_dir, "Rows")
    idx = _read_csv(os.path.join(rdir, "rows_index.csv"))
    res = []
    for i, r in enumerate(idx):
        if cancel is not None and cancel.is_set():
            break
        rec = dict(file=r["file"], plot=r["plot_stem"], row=r["row"], middle=r.get("middle", ""))
        try:
            src = os.path.join(crops_dir, r["plot_stem"] + ".tif")
            if not os.path.exists(src):
                src = os.path.join(crops_dir, r["plot_stem"] + ".tiff")
            a = tiler.read_rgb(os.path.join(rdir, r.get("rel_path") or r["file"]))
            with rasterio.open(src) as s:
                ref = s.read([1, 2, 3], window=Window(int(r["x0"]), int(r["y0"]),
                                                      int(r["w"]), int(r["h"]))).transpose(1, 2, 0)
            if ref.dtype != a.dtype:
                ref = tiler.convert_bits(ref, 16 if a.dtype == np.uint16 else 8)
            kept = a.max(axis=2) > 0
            same = bool(ref.shape == a.shape and np.array_equal(ref[kept], a[kept]))
            rec.update(kept_px=int(kept.sum()), pixels="IDENTICAL" if same else "FAIL")
        except Exception as e:
            rec.update(pixels="ERROR", note=str(e))
        res.append(rec)
        if progress:
            progress(i + 1, len(idx))
    return res


def run_checks(run_dir, tiles_dir="", log=print, progress=None, cancel=None):
    """Run every check that applies; write QA_REPORT.txt + CSVs. Returns report path."""
    out = os.path.join(tiles_dir or run_dir, QA_DIR)
    os.makedirs(out, exist_ok=True)
    lines = ["PlotExtractor QA report", f"extraction: {run_dir or '-'}",
             f"tiles:      {tiles_dir or '-'}", ""]
    crops_dir = _crops_dir(run_dir) if run_dir else ""
    flags = {}
    if run_dir and os.path.exists(os.path.join(run_dir, "manifest.csv")):
        log("QA: checking the crops against the source frames ...")
        cr = check_crops(run_dir, log, cancel, progress)
        _write_csv(cr, os.path.join(out, "qa_crops.csv"))
        c = Counter(r["result"] for r in cr)
        lines.append("CROPS (" + str(len(cr)) + "): " + ", ".join(f"{k} {v}" for k, v in c.items()))
        iiq_src = any(r.get("source", "").lower().endswith(".iiq") for r in cr)
        if c.get("IDENTICAL"):
            px = sum(r.get("inside_px", 0) for r in cr if r["result"] == "IDENTICAL")
            lines.append(f"  {c['IDENTICAL']} raw crop(s) identical to "
                         + ("a fresh development of their IIQ box" if iiq_src
                            else "the source frames")
                         + f" ({px:,} pixels compared, every one equal).")
        if c.get("MATCH"):
            lines.append(f"  {c['MATCH']} crop(s) within 2/255 of a fresh IIQ development: "
                         "the settings use noise reduction, which the Phase One SDK does "
                         "not repeat bit for bit (neutral settings do).")
        if c.get("RESAMPLED"):
            lines.append("  rectified crops are re-rendered on a ground grid (resampled by "
                         "design); frames used: " + ", ".join(
                             f"{k} frame(s) x{v}" for k, v in
                             sorted(Counter(r.get('frames_used', '') for r in cr
                                            if r['result'] == 'RESAMPLED').items())))
        for r in cr:
            if r["result"] not in ("IDENTICAL", "RESAMPLED", "MATCH"):
                flags.setdefault("plot_" + str(r["plot"]), []).append("crop: " + r["result"])
        lines.append("")
    if tiles_dir and os.path.exists(os.path.join(tiles_dir, "tile_index.csv")):
        log("QA: checking the tiles against the crops ...")
        tr = check_tiles(tiles_dir, crops_dir, log, cancel, progress)
        _write_csv(tr, os.path.join(out, "qa_tiles.csv"))
        c = Counter(r.get("pixels", "") for r in tr)
        lines.append(f"TILES ({len(tr)}): " + ", ".join(f"{k} {v}" for k, v in c.items()))
        outside = [r for r in tr if r.get("in_usable_length") == "NO"]
        if any("in_usable_length" in r for r in tr):
            lines.append(f"  outside the usable plot length: {len(outside)}")
        per = Counter(r["plot"] for r in tr)
        if per:
            lines.append(f"  tiles per plot: min {min(per.values())}, max {max(per.values())}, "
                         f"{len(per)} plot(s)")
        for r in tr:
            if r.get("pixels") in ("FAIL", "ERROR") or r.get("in_usable_length") == "NO":
                flags.setdefault(r["plot"], []).append(f"tile {r['tile']}: {r.get('pixels')}"
                                                       + (" outside usable length"
                                                          if r.get("in_usable_length") == "NO" else ""))
        lines.append("")
        rr, usual, med = check_rows(tiles_dir)
        if rr:
            _write_csv(rr, os.path.join(out, "qa_rows.csv"))
            nf = [r for r in rr if r["flag"]]
            lines.append(f"ROWS ({len(rr)} plots): usual {usual} rows per plot, median row "
                         f"spacing {med:.1f} {'mm' if rr and rr[0]['pitch_mm'] else 'px'}; "
                         f"{len(rr) - len(nf)} plot(s) without a flag, {len(nf)} flagged")
            for r in nf:
                flags.setdefault(r["plot"], []).append("rows: " + r["flag"])
            lines.append("")
        if os.path.exists(os.path.join(tiles_dir, "Rows", "rows_index.csv")):
            log("QA: checking the saved row images ...")
            ri = check_row_images(tiles_dir, crops_dir, cancel, progress)
            _write_csv(ri, os.path.join(out, "qa_row_images.csv"))
            c = Counter(r.get("pixels", "") for r in ri)
            lines.append(f"ROW IMAGES ({len(ri)}): " + ", ".join(f"{k} {v}" for k, v in c.items())
                         + "  (pixels inside each row equal the crop)")
            for r in ri:
                if r.get("pixels") != "IDENTICAL":
                    flags.setdefault(r["plot"], []).append(f"row image {r['file']}: {r.get('pixels')}")
            lines.append("")
    lines.append(f"FLAGGED PLOTS: {len(flags)}")
    for k in sorted(flags):
        lines.append(f"  {k}: " + "; ".join(flags[k]))
    path = os.path.join(out, "QA_REPORT.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    _write_csv([dict(plot=k, flags="; ".join(v)) for k, v in sorted(flags.items())],
               os.path.join(out, "qa_flags.csv"))
    for ln in lines:
        log(ln)
    return path


# ---------------------------------------------------------------------------
# pictures
# ---------------------------------------------------------------------------
def rows_picture(crop_path, rows_expected=0, n_middle=3, max_side=1600):
    """Row detection of one crop drawn on a reduced copy (used when the tiles run
    wrote no _rows_qa picture). Returns (PIL image, detection dict)."""
    import rasterio
    import row_detect as rd
    import tiler
    with rasterio.open(crop_path) as s:
        W, H = s.width, s.height
        dec = max(1, int(round(min(W, H) / 1200.0)))
        small = tiler.to8(s.read([1, 2, 3], out_shape=(3, max(1, H // dec), max(1, W // dec))))
        mm = tiler._plot_mm_px(s)
    det = rd.detect_rows(small, mm_px=(mm * W / small.shape[2]) if mm else None,
                         rows_expected=rows_expected)
    band = rd.middle_band(det, n_middle) if n_middle else None
    return rd.draw_rows(np.transpose(small, (1, 2, 0)), det, band=band, max_side=max_side), det


def _thumb(path, box):
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    try:
        with Image.open(path) as im:
            f = max(1, int(max(im.size[0] / box[0], im.size[1] / box[1])))
            if f > 1:
                im = im.reduce(f)
            im = im.convert("RGB")
            im.thumbnail(box)
            return im
    except Exception:                  # 16-bit RGB (PIL has no such mode): via GDAL
        import tiler
        im = Image.fromarray(tiler.to8(tiler.read_rgb(path)))
        im.thumbnail(box)
        return im


def montage(paths, cell=(360, 360), cols=None, labels=None, gap=8):
    """Thumbnails of several images in a grid (tiles, row images)."""
    from PIL import Image, ImageDraw
    ims = [_thumb(p, cell) for p in paths]
    if not ims:
        return None
    n = len(ims)
    if cols is None:
        cols = n if n <= 8 else int(np.ceil(np.sqrt(n)))
    rows = int(np.ceil(n / cols))
    cw = max(i.size[0] for i in ims); ch = max(i.size[1] for i in ims) + 16
    out = Image.new("RGB", (cols * (cw + gap) + gap, rows * (ch + gap) + gap), (40, 40, 40))
    d = ImageDraw.Draw(out)
    for k, im in enumerate(ims):
        x = gap + (k % cols) * (cw + gap); y = gap + (k // cols) * (ch + gap)
        out.paste(im, (x, y + 16))
        d.text((x + 2, y + 2), (labels[k] if labels else os.path.basename(paths[k])),
               fill=(255, 255, 255))
    return out


# ---------------------------------------------------------------------------
# tab
# ---------------------------------------------------------------------------
def build_qa_tab(parent, defaults=None):
    """'QA' tab: open an extraction (and its tiles), browse plots, run the checks."""
    import queue
    import threading
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    from PIL import ImageTk

    d = defaults or {}
    frm = ttk.Frame(parent, padding=10)
    frm.columnconfigure(1, weight=1)
    q = queue.Queue()
    cancel_ev = threading.Event()
    st = {"busy": False, "img": None, "pil": None, "plots": [], "flags": {},
          "rep": {}, "tiles": {}, "rows": {}}
    V = {"run": tk.StringVar(value=d.get("run", "")),
         "tiles": tk.StringVar(value=d.get("tiles", "")),
         "view": tk.StringVar(value="rows"),
         "nmid": tk.StringVar(value="3"), "nrows": tk.StringVar(value="0")}

    r = 0
    ttk.Label(frm, text="Check an extraction and its tiles",
              font=("", 10, "bold")).grid(row=r, column=0, columnspan=3, sticky="w")
    r += 1
    for label, var, hint in (("Extraction folder *", V["run"],
                              "the folder holding RawCrops\\ and manifest.csv"),
                             ("Tiles folder", V["tiles"],
                              "the Training-tiles output (tile_index.csv, rows_report.csv, "
                              "_rows_qa\\, Rows\\); blank = crops only")):
        ttk.Label(frm, text=label).grid(row=r, column=0, sticky="w", pady=2)
        ttk.Entry(frm, textvariable=var).grid(row=r, column=1, sticky="ew", padx=6)

        def browse(v=var):
            p = filedialog.askdirectory()
            if p:
                v.set(p)
                if v is V["tiles"] and not V["run"].get().strip():
                    V["run"].set(find_run_for_tiles(p))
                load()
        ttk.Button(frm, text="Browse...", command=browse).grid(row=r, column=2)
        r += 1
        ttk.Label(frm, text=hint, foreground="grey").grid(row=r, column=1, sticky="w", padx=6)
        r += 1

    bar = ttk.Frame(frm)
    bar.grid(row=r, column=0, columnspan=3, sticky="w", pady=4)
    r += 1
    ttk.Button(bar, text="Load plots", command=lambda: load()).pack(side="left")
    run_btn = ttk.Button(bar, text="Run QA checks")
    run_btn.pack(side="left", padx=8)
    cancel_btn = ttk.Button(bar, text="Cancel", state="disabled",
                            command=lambda: cancel_ev.set())
    cancel_btn.pack(side="left")
    prog = ttk.Progressbar(bar, mode="determinate", length=220)
    prog.pack(side="left", padx=8)
    ttk.Label(bar, text="show:").pack(side="left", padx=(12, 4))
    for txt, val in (("rows + tiles", "rows"), ("tiles", "tiles"), ("saved rows", "rowimgs"),
                     ("crop", "crop"), ("QA report", "report")):
        ttk.Radiobutton(bar, text=txt, value=val, variable=V["view"],
                        command=lambda: show()).pack(side="left", padx=2)
    ttk.Label(bar, text="  middle rows:").pack(side="left")
    ttk.Spinbox(bar, from_=0, to=12, width=3, textvariable=V["nmid"]).pack(side="left")
    ttk.Label(bar, text=" sown rows:").pack(side="left")
    ttk.Spinbox(bar, from_=0, to=24, width=3, textvariable=V["nrows"]).pack(side="left")

    pw = ttk.PanedWindow(frm, orient="horizontal")
    pw.grid(row=r, column=0, columnspan=3, sticky="nsew", pady=4)
    frm.rowconfigure(r, weight=1)
    r += 1
    left = ttk.Frame(pw)
    lb = tk.Listbox(left, width=26, font=("Consolas", 9), exportselection=False)
    lsb = ttk.Scrollbar(left, command=lb.yview)
    lb.configure(yscrollcommand=lsb.set)
    lb.pack(side="left", fill="both", expand=True)
    lsb.pack(side="left", fill="y")
    pw.add(left, weight=0)
    right = ttk.Frame(pw)
    canvas = tk.Canvas(right, background="#202020", highlightthickness=0)
    canvas.pack(fill="both", expand=True)
    info = tk.Text(right, height=5, wrap="word", font=("Consolas", 9), state="disabled")
    info.pack(fill="x")
    pw.add(right, weight=1)

    logbox = tk.Text(frm, height=6, wrap="none", font=("Consolas", 9), state="disabled")
    logbox.grid(row=r, column=0, columnspan=3, sticky="ew")

    def put(m):
        q.put(("log", str(m)))

    def set_text(w, txt):
        w.configure(state="normal"); w.delete("1.0", "end"); w.insert("end", txt)
        w.configure(state="disabled")

    def crops_dir():
        return _crops_dir(V["run"].get().strip())

    def load():
        run = V["run"].get().strip(); td = V["tiles"].get().strip()
        if td and not run:
            run = find_run_for_tiles(td); V["run"].set(run)
        cd = _crops_dir(run) if run else ""
        plots = set()
        if cd and os.path.isdir(cd):
            plots |= {os.path.splitext(f)[0] for f in os.listdir(cd)
                      if f.lower().endswith((".tif", ".tiff"))}
        st["rep"] = {x["plot_stem"]: x for x in _read_csv(os.path.join(td, "rows_report.csv"))} if td else {}
        st["tiles"] = {}
        for t in (_read_csv(os.path.join(td, "tile_index.csv")) if td else []):
            st["tiles"].setdefault(t["plot_stem"], []).append(t)
        st["rows"] = {}
        for x in (_read_csv(os.path.join(td, "Rows", "rows_index.csv")) if td else []):
            st["rows"].setdefault(x["plot_stem"], []).append(x)
        plots |= set(st["rep"]) | set(st["tiles"])
        qa = os.path.join(td or run, QA_DIR, "qa_flags.csv") if (td or run) else ""
        st["flags"] = {x["plot"]: x["flags"] for x in _read_csv(qa)}
        st["plots"] = sorted(plots)
        lb.delete(0, "end")
        for p in st["plots"]:
            rr = st["rep"].get(p, {})
            mark = "!" if p in st["flags"] else (" " if not st["flags"] else "ok")
            lb.insert("end", f"{mark:>2} {p:<14} {rr.get('rows_found', ''):>2}r "
                             f"{len(st['tiles'].get(p, [])):>2}t")
            if p in st["flags"]:
                lb.itemconfigure("end", foreground="#c00000")
        put(f"{len(st['plots'])} plot(s) loaded"
            + (f", {len(st['flags'])} flagged by the last QA run" if st["flags"] else ""))
        if st["plots"]:
            lb.selection_clear(0, "end"); lb.selection_set(0); show()

    def current():
        sel = lb.curselection()
        return st["plots"][sel[0]] if sel and st["plots"] else None

    def draw(pil):
        st["pil"] = pil
        canvas.delete("all")
        if pil is None:
            return
        cw, ch = max(50, canvas.winfo_width()), max(50, canvas.winfo_height())
        im = pil.copy()
        im.thumbnail((cw, ch))
        st["img"] = ImageTk.PhotoImage(im)
        canvas.create_image(cw // 2, ch // 2, image=st["img"])

    def show(*_):
        p = current()
        td = V["tiles"].get().strip()
        view = V["view"].get()
        if view == "report":
            rp = os.path.join(td or V["run"].get().strip(), QA_DIR, "QA_REPORT.txt")
            txt = open(rp, encoding="utf-8").read() if os.path.exists(rp) else \
                "No QA report yet - press 'Run QA checks'."
            draw(None)
            canvas.create_text(10, 10, anchor="nw", text=txt, fill="white",
                               font=("Consolas", 9))
            return
        if not p:
            return
        rr = st["rep"].get(p, {})
        lines = [f"{p}   " + ("  ".join(f"{k}={rr[k]}" for k in
                                         ("rows_found", "pitch_mm", "angle_deg", "middle_rows",
                                          "usable_len_m", "tiles_drawn", "veg_mode")
                                         if rr.get(k, "") != ""))]
        if p in st["flags"]:
            lines.append("QA FLAG: " + st["flags"][p])
        frm.configure(cursor="watch"); frm.update_idletasks()
        try:
            pil = None
            crop = os.path.join(crops_dir(), p + ".tif")
            if view == "rows":
                qa = os.path.join(td, "_rows_qa", p + "_rows.jpg") if td else ""
                if qa and os.path.exists(qa):
                    from PIL import Image
                    pil = Image.open(qa).convert("RGB")
                    lines.append("green = rows (thick = the plot's rows, red = missing), yellow "
                                 "= middle band, grey = plot ends, orange = usable length, "
                                 "cyan = tiles written, blue = other candidates")
                elif os.path.exists(crop):
                    pil, det = rows_picture(crop, int(V["nrows"].get() or 0),
                                            int(V["nmid"].get() or 0))
                    lines.append(f"detected now: {len(det.get('plot_rows', []))} rows, "
                                 + (f"angle {det['angle_deg']:.2f} deg, " if det.get("ok") else "")
                                 + (det.get("reason") or "") +
                                 "  (green = rows, yellow = middle band)")
            elif view == "tiles":
                ts = st["tiles"].get(p, [])
                pil = montage([os.path.join(td, t.get("rel_path") or t["tile"]) for t in ts],
                              labels=[f"{k + 1}: {t['tile']}  veg {t.get('veg_pct', '')}%"
                                      for k, t in enumerate(ts)]) if ts else None
                if not ts:
                    lines.append("no tiles for this plot")
            elif view == "rowimgs":
                rs = st["rows"].get(p, [])
                if rs:
                    vertical = int(rs[0]["h"]) > int(rs[0]["w"])
                    cell = (220, 1400) if vertical else (1400, 160)
                    pil = montage([os.path.join(td, "Rows", x.get("rel_path") or x["file"]) for x in rs],
                                  cell=cell, cols=len(rs) if vertical else 1,
                                  labels=[f"row {x['row']}" + (" (middle)" if x.get("middle") == "True" else "")
                                          for x in rs])
                else:
                    lines.append("no saved rows (tick 'also save each row' on the Training "
                                 "tiles tab)")
            elif view == "crop" and os.path.exists(crop):
                pil = _thumb(crop, (1600, 1600))
            draw(pil)
        except Exception as e:
            lines.append("[!] " + str(e)); draw(None)
        finally:
            frm.configure(cursor="")
        set_text(info, "\n".join(lines))

    def on_run():
        run = V["run"].get().strip(); td = V["tiles"].get().strip()
        if not run and not td:
            messagebox.showwarning("QA", "Choose the extraction folder.")
            return
        if st["busy"]:
            return
        st["busy"] = True; cancel_ev.clear()
        run_btn.configure(state="disabled"); cancel_btn.configure(state="normal")

        def work():
            try:
                run_checks(run, td, log=put, progress=lambda a, b: q.put(("prog", (a, b))),
                           cancel=cancel_ev)
            except Exception as e:
                put("[!] " + str(e))
            q.put(("done", None))
        threading.Thread(target=work, daemon=True).start()

    run_btn.configure(command=on_run)
    lb.bind("<<ListboxSelect>>", show)
    canvas.bind("<Configure>", lambda e: draw(st["pil"]) if st["pil"] is not None else None)

    def poll():
        try:
            while True:
                kind, payload = q.get_nowait()
                if kind == "log":
                    logbox.configure(state="normal"); logbox.insert("end", payload + "\n")
                    logbox.see("end"); logbox.configure(state="disabled")
                elif kind == "prog":
                    prog["maximum"], prog["value"] = payload[1], payload[0]
                elif kind == "done":
                    st["busy"] = False
                    run_btn.configure(state="normal"); cancel_btn.configure(state="disabled")
                    load()
        except queue.Empty:
            pass
        frm.after(150, poll)

    poll()
    return frm
