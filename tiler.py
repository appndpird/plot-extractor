"""tiler.py - cut per-plot crops into training tiles, and put them back together.

Built for the labelling loop: cut tiles -> label them (SAM3 or by hand) -> stitch
the results back into one image per plot to see what you have.

Design notes
------------
* Tile names keep the convention already in use in the group's scripts,
  ``<plot stem>_tile_c<col>_r<row>.<ext>``, so tiles drop straight into
  sam3_batch_tiles.ipynb, which just reads a flat folder of images.

* Every run writes ``tile_index.csv``. That is what makes the tiles reversible:
  it records, per tile, which plot it came from and where it sat. Without it a
  detection in a tile cannot be traced back to a plot, and overlapping tiles
  double-count.

* ``rows_only`` restricts tiling to the band that actually holds the seed rows.
  A plot crop includes margin beyond the outermost rows; tiles cut there are
  bare soil and waste labelling effort. The band is measured per plot from the
  vegetation profile, not assumed.

* Tiles are written LOSSLESS by default (PNG). They are training data and are
  labelled pixel-by-pixel; JPEG ringing on a leaf edge becomes a wrong label.
"""
from __future__ import annotations

import csv
import math
import os
import random

import numpy as np


# ---------------------------------------------------------------------------
# row-band detection
# ---------------------------------------------------------------------------
def _smooth(v, k):
    if k < 3:
        return v
    return np.convolve(v, np.ones(k, np.float32) / k, mode="same")


def _peaks(p, min_dist, prom):
    n = len(p)
    cand = [i for i in range(1, n - 1) if p[i] >= p[i - 1] and p[i] > p[i + 1]]
    cand.sort(key=lambda i: -p[i])
    keep = []
    for i in cand:
        if any(abs(i - j) < min_dist for j in keep):
            continue
        lo = max(0, i - int(min_dist)); hi = min(n, i + int(min_dist) + 1)
        left = p[lo:i].min() if i > lo else 0.0
        right = p[i + 1:hi].min() if hi > i + 1 else 0.0
        if p[i] - max(left, right) < prom:
            continue
        keep.append(i)
    return sorted(keep)


def row_band(arr, across_axis, pitch_px, margin_frac=0.5):
    """Return (lo, hi) pixel bounds of the seed-row band along `across_axis`.

    Falls back to the full extent when rows cannot be found - a plot with too
    little emergence still gets tiled, it just is not trimmed.
    """
    R, G, B = arr[0].astype(np.float32), arr[1].astype(np.float32), arr[2].astype(np.float32)
    exg = 2 * G - R - B
    valid = arr.max(axis=0) > 4
    axis = 1 if across_axis == 0 else 0
    cnt = valid.sum(axis=axis)
    num = np.where(valid, exg, 0).sum(axis=axis)
    prof = np.where(cnt > 0, num / np.maximum(cnt, 1e-6), 0.0)
    n = len(prof)
    if n < 10 or pitch_px < 4:
        return 0, n
    prof = _smooth(prof, max(3, int(pitch_px * 0.30) | 1))
    rng = float(np.ptp(prof))
    if rng <= 1e-6:
        return 0, n
    p = (prof - prof.min()) / rng
    pk = _peaks(p, pitch_px * 0.55, 0.09)
    if len(pk) < 2:
        return 0, n
    pad = pitch_px * margin_frac
    return max(0, int(pk[0] - pad)), min(n, int(pk[-1] + pad) + 1)


# ---------------------------------------------------------------------------
# tiling
# ---------------------------------------------------------------------------
def _starts(a, b, tile, stride):
    """Tile origins covering [a, b), with the LAST tile flush against b.

    A plain range(a, b - tile + 1, stride) stops short whenever the span is not
    a whole number of strides, so the far edge of every plot goes untiled - on a
    2418 px band with 1024 px tiles and 768 stride that silently drops the last
    626 px. The final tile is therefore pinned to b - tile, overlapping its
    neighbour a little more than the nominal overlap.
    """
    if b - a <= tile:
        return [a]
    xs = list(range(a, b - tile + 1, stride))
    if not xs:
        xs = [a]
    if xs[-1] + tile < b:
        xs.append(b - tile)
    return xs


def to8(a):
    """uint8 view of image data for analysis (row detection, vegetation, previews):
    uint16 is mapped to 0-255 by v/257 rounded (the exact inverse of v*257)."""
    a = np.asarray(a)
    if a.dtype == np.uint8:
        return a
    return ((a.astype(np.uint32) + 128) // 257).clip(0, 255).astype(np.uint8)


def convert_bits(a, bits):
    """Pixel values at the asked bit depth. bits 'same' (or the array's own depth)
    returns the array unchanged - the exact copy. 16 -> 8 is v/257 rounded,
    8 -> 16 is v*257 (0 -> 0, 255 -> 65535; reversible)."""
    b = str(bits).strip().lower()
    if b == "8" and a.dtype != np.uint8:
        return to8(a)
    if b == "16" and a.dtype == np.uint8:
        return a.astype(np.uint16) * 257
    return a


def read_rgb(path):
    """HxWx3 array of an image file with its own bit depth (8 or 16): TIFF/PNG via
    GDAL, anything else via PIL."""
    if path.lower().endswith((".tif", ".tiff", ".png")):
        import warnings
        import rasterio
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with rasterio.open(path) as s:
                return s.read([1, 2, 3]).transpose(1, 2, 0)
    from PIL import Image
    with Image.open(path) as im:
        return np.asarray(im.convert("RGB"))


def write_tile(dest, a, ext, compress="lzw", crs=None, transform=None, gcps=None,
               tags=None):
    """Write a 3 x h x w tile, values unchanged (8 or 16 bit). TIFF and PNG are
    written with GDAL - TIFF as GeoTIFF when crs + transform (or gcps) are given,
    LZW with horizontal predictor or uncompressed; JPEG (8-bit only) with PIL."""
    if ext == ".jpg":
        from PIL import Image
        if a.dtype != np.uint8:
            raise ValueError("JPEG is 8-bit only")
        Image.fromarray(np.transpose(a, (1, 2, 0))).save(dest, quality=95, subsampling=0)
        return
    import warnings
    import rasterio
    prof = dict(height=a.shape[1], width=a.shape[2], count=3, dtype=a.dtype.name)
    if ext == ".png":
        prof["driver"] = "PNG"
    else:
        prof.update(driver="GTiff", photometric="RGB")
        if compress != "none":
            prof.update(compress="lzw", predictor=2)
        if crs is not None and transform is not None:
            prof.update(crs=crs, transform=transform)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")            # NotGeoreferencedWarning
        with rasterio.open(dest, "w", **prof) as dst:
            dst.write(a)
            if ext != ".png" and gcps and transform is None:
                dst.gcps = gcps
            if tags:
                dst.update_tags(**{k: str(v) for k, v in tags.items() if v != ""})


def _tile_georef(s, x, y, w, h, out_w, out_h):
    """(crs, transform, gcps) of the tile window (x, y, w, h) of the open crop s,
    written at out_w x out_h px. Affine crops: shifted (and scaled) transform;
    GCP-only crops: the GCPs moved into the tile. (None, None, None) if the crop
    has no georeference."""
    from rasterio.control import GroundControlPoint
    from rasterio.transform import Affine
    sx, sy = w / float(out_w), h / float(out_h)
    if s.crs is not None and not s.transform.is_identity:
        return s.crs, s.transform * Affine.translation(x, y) * Affine.scale(sx, sy), None
    g, gcrs = s.gcps
    if g and gcrs is not None:
        moved = [GroundControlPoint(row=(p.row - y) / sy, col=(p.col - x) / sx,
                                    x=p.x, y=p.y, z=p.z) for p in g]
        return gcrs, None, (moved, gcrs)
    return None, None, None


def _source_links(src_dir):
    """Per crop file: the source photo / IIQ and the crop's box in it, from the
    extraction's manifest.csv (next to RawCrops) and IIQ_run.json."""
    import json
    run = os.path.dirname(os.path.abspath(src_dir))
    links = {}
    for mf in (os.path.join(run, "manifest.csv"), os.path.join(src_dir, "manifest.csv")):
        if os.path.exists(mf):
            try:
                with open(mf, newline="", encoding="utf-8") as fh:
                    for rr in csv.DictReader(fh):
                        if rr.get("raw_file"):
                            links[rr["raw_file"]] = rr
            except Exception:
                pass
            break
    iiq = {}
    try:
        with open(os.path.join(run, "IIQ_run.json")) as fh:
            iiq = json.load(fh)
    except Exception:
        pass
    settings = os.path.join(run, "IIQ_settings_used.json")
    return links, iiq, (settings if os.path.exists(settings) else "")


def _link_fields(man, crop_tags, x, y, w, h, iiq_settings):
    """Index columns / GeoTIFF tags tying a tile to its source photo or IIQ. The
    pixel box in the photo exists for crops that are a 1:1 box of one frame
    (masked / bounding box / native); rectified or stitched crops were resampled."""
    src = (man or {}).get("source_path") or crop_tags.get("SOURCE_IMAGE", "")
    pix = (man or {}).get("pixel_source") or ("iiq" if src.lower().endswith(".iiq") else
                                             ("image" if src else ""))
    style = (man or {}).get("raw_style", "")
    box = ""
    one_frame = str((man or {}).get("frames_used") or "1").strip() in ("", "1")
    if man and style in ("mask", "bbox", "native") and one_frame and man.get("src_x0", "") != "":
        try:
            box = (int(float(man["src_x0"])) + x, int(float(man["src_y0"])) + y, w, h)
        except ValueError:
            box = ""
    return dict(source_image=src, pixel_source=pix,
                src_photo_x0=box[0] if box else "", src_photo_y0=box[1] if box else "",
                src_photo_w=box[2] if box else "", src_photo_h=box[3] if box else "",
                iiq_settings=(iiq_settings if pix == "iiq" else ""),
                link_note=("" if box else
                           ("crop was resampled (rectified/stitched): no 1:1 pixel box"
                            if man else "no manifest.csv: pixel box unknown")))


def _plot_mm_px(src):
    """Ground sampling of an open raster, mm per pixel, from its GCPs.

    Returns None when the file carries no usable georeferencing.
    """
    try:
        if src.crs is not None and not src.transform.is_identity:
            # affine-georeferenced crop (the normal GeoTIFF output): pixel size
            t = src.transform
            px = math.sqrt(abs(t.a * t.e - t.b * t.d))
            if src.crs.is_geographic:
                lat = math.radians((t * (src.width / 2.0, src.height / 2.0))[1])
                px *= math.sqrt(111320.0 * math.cos(lat) * 110540.0)
            else:
                px *= src.crs.linear_units_factor[1] if src.crs.linear_units_factor else 1.0
            return px * 1000.0 if px > 0 else None
        g, crs = src.gcps
        if len(g) < 2 or crs is None:
            return None
        geographic = crs.is_geographic
        pts = [(q.col, q.row, q.x, q.y) for q in g]
        sc = []
        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                dpx = math.hypot(pts[i][0] - pts[j][0], pts[i][1] - pts[j][1])
                if dpx < 50:
                    continue
                dx, dy = pts[i][2] - pts[j][2], pts[i][3] - pts[j][3]
                if geographic:
                    lat = math.radians((pts[i][3] + pts[j][3]) / 2.0)
                    dx *= 111320.0 * math.cos(lat)
                    dy *= 110540.0
                sc.append(math.hypot(dx, dy) * 1000.0 / dpx)
        return float(np.median(sc)) if len(sc) >= 1 else None
    except Exception:
        return None


def _usable_length(small, across_axis, det, band, mm_small, ext_margin_mm, end_margin_mm):
    """Along-plot span (reduced-copy pixels) where middle-row tiles may sit.

    The plot ends are bare soil (or the black mask / frame edge of a raw crop), and
    tiles there teach a model that empty ground is part of a plot. So:
      data ends     first/last along position where the band is fully imaged
      plot ends     data ends moved in by the extraction margin (manifest margin_m)
      planted ends  where the vegetation share in the band reaches 35 % of its
                    median over the central half of the plot (only used when the
                    planted length is at least half the plot: sparse emergence
                    must not shrink the plot to a few clumps)
      usable        the tighter of plot/planted ends, moved in by end_margin_mm.
    Returns dict(u0, u1, d0, d1, p0, p1, q0, q1, trim) in reduced-copy pixels;
    q0, q1 = the plot ends before any planted-length trim."""
    a = small[:3].astype(np.float32)
    if across_axis == 1:                        # make axis 1 = across, axis 2 = along
        a = a.transpose(0, 2, 1)
    R, G, B = a[0], a[1], a[2]
    valid = a.max(axis=0) > 4
    with np.errstate(all="ignore"):
        exg = (2 * G - R - B) / np.maximum(R + G + B, 1e-6)
    thr = det.get("veg_threshold", np.nan)
    veg = valid & (exg > (thr if np.isfinite(thr) else 0.05))
    na, nl = valid.shape
    t = math.tan(math.radians(det["angle_deg"]))
    lc = det["along_centre"]
    u = np.arange(na, dtype=np.float32)[:, None] - t * (np.arange(nl, dtype=np.float32)[None, :] - lc)
    inb = (u >= band[0]) & (u < band[1])
    bcnt = inb.sum(axis=0)
    vcnt = (valid & inb).sum(axis=0)
    data = (bcnt > 0) & (vcnt >= 0.9 * np.maximum(bcnt, 1))
    idx = np.nonzero(data)[0]
    if not len(idx):
        return dict(u0=0, u1=nl, d0=0, d1=nl, p0=0, p1=nl, q0=0, q1=nl, trim="none")
    d0, d1 = int(idx[0]), int(idx[-1]) + 1
    px_per_mm = (1.0 / mm_small) if mm_small else None
    em = int(round(ext_margin_mm * px_per_mm)) if (px_per_mm and ext_margin_mm) else 0
    p0, p1 = d0 + em, d1 - em
    if p1 - p0 < 8:
        p0, p1 = d0, d1
    q0, q1 = p0, p1
    trim = "plot ends"
    gfrac = np.where(data, (veg & inb).sum(axis=0) / np.maximum(vcnt, 1), 0.0)
    k = int(round(300.0 * px_per_mm)) if px_per_mm else max(5, nl // 25)
    g = _smooth(gfrac.astype(np.float32), max(3, k | 1))
    q = (p1 - p0) // 4
    ref = float(np.median(g[p0 + q:p1 - q])) if p1 - p0 > 8 else 0.0
    if ref > 0:
        on = np.nonzero((g >= 0.35 * ref) & data)[0]
        on = on[(on >= p0) & (on < p1)]
        if len(on) and (on[-1] + 1 - on[0]) >= 0.5 * (p1 - p0):
            if on[0] > p0 or on[-1] + 1 < p1:
                trim = "planted ends"
            p0, p1 = max(p0, int(on[0])), min(p1, int(on[-1]) + 1)
    if px_per_mm:
        e = int(round(end_margin_mm * px_per_mm))
    else:                                       # no ground scale: 8 % of the plot
        e = int(round(0.08 * (p1 - p0)))
    return dict(u0=p0 + e, u1=p1 - e, d0=d0, d1=d1, p0=p0, p1=p1, q0=q0, q1=q1, trim=trim)


def _pick(cand, n, how, rnd):
    """n of the candidate tiles: 'random' sample or 'spread' (evenly spaced along
    the plot). n = 0 or n >= len(cand) keeps them all."""
    if not n or n >= len(cand):
        return list(cand)
    if how == "spread":
        idx = sorted(set(int(round(i)) for i in np.linspace(0, len(cand) - 1, n)))
        return [cand[i] for i in idx]
    return sorted(rnd.sample(cand, n))


def _small_copy(s, W, H):
    """8-bit reduced copy (shorter side ~1200 px) used for row detection and QA."""
    dec = max(1, int(round(min(W, H) / 1200.0)))
    return to8(s.read([1, 2, 3], out_shape=(3, max(1, H // dec), max(1, W // dec))))


def _detect_plot(s, W, H, across_axis, src_mm, row_pitch_mm, rows_expected, small=None,
                 det=None):
    """Rows of one plot on a reduced copy: (small, det, f, fl); f / fl = across /
    along scale from the reduced copy to the crop. det: an earlier detection of the
    same copy (it is then only re-read, not re-detected)."""
    import row_detect as rd
    if small is None:
        small = _small_copy(s, W, H)
    sx = W / float(small.shape[2]); sy = H / float(small.shape[1])
    if det is None:
        det = rd.detect_rows(small, across_axis=across_axis,
                             mm_px=(src_mm * sx) if src_mm else None,
                             pitch_mm=row_pitch_mm or None, rows_expected=rows_expected)
    f = sy if across_axis == 0 else sx
    fl = sx if across_axis == 0 else sy
    return small, det, f, fl


def _middle_row_cands(s, W, H, across_axis, n_middle, rows_expected, src_mm,
                      row_pitch_mm, band_tile, win_px, overlap, qa_path=None,
                      ext_margin_mm=0.0, end_margin_mm=500.0, strict=False,
                      win_across=None, det_in=None):
    """Tile candidates covering the n_middle central crop rows of one plot.

    Rows are found with row_detect (the PhenoFuse RGB row detection) on a reduced
    copy of the crop; the band of the middle rows follows the rows' tilt, so it also
    works on unrectified raw crops. band_tile: square tiles whose side is the band
    width, laid end to end along the plot (native pixels, no overlap); otherwise
    tiles win_px along x win_across across (default square) centred on the band.
    Along the plot the tiles stay inside the usable length (_usable_length: away
    from the bare plot ends by end_margin_mm) and are centred in it. strict: the
    tile size is fixed (ground tiles); a plot narrower than one tile gets none
    instead of a smaller tile. det_in: (small, det, f, fl) from _detect_plot.

    Returns (cand, info): cand = [(col, row, x, y, w, h)], info = per-plot report;
    info['_qa'] is set on every path when qa_path is given (skipped plots too)."""
    import row_detect as rd
    small, det, f, fl = det_in or _detect_plot(s, W, H, across_axis, src_mm,
                                               row_pitch_mm, rows_expected)
    info = dict(ok=det.get("ok", False), reason=det.get("reason", ""),
                rows_found=len(det.get("plot_rows", [])),
                rows_present=det.get("n_rows_present", 0),
                pitch_px=round(det["pitch_px"] * f, 1) if det.get("ok") else "",
                pitch_mm=(round(det["pitch_px"] * f * src_mm, 1)
                          if (det.get("ok") and src_mm) else ""),
                angle_deg=round(det.get("angle_deg", 0.0), 2), band_lo="", band_hi="",
                veg_mode=det.get("veg_mode", ""),
                middle_rows="")
    qa = dict(path=qa_path, small=small, det=det, band=None, us=None, fl=fl)
    info["_qa"] = qa if qa_path else None
    if not det.get("ok"):
        info["reason"] = info["reason"] or "rows not found"
        return [], info
    info["_det"] = (det, f, fl)
    band = rd.middle_band(det, n_middle)
    if band is None:
        info["reason"] = (f"only {len(det['plot_rows'])} row(s) found, "
                          f"{n_middle} middle rows asked")
        return [], info
    qa["band"] = band
    lo, hi, blk = band[0] * f, band[1] * f, band[2]
    info["_det"] = (det, f, fl, blk)
    k0 = det["plot_rows"].index(blk[0])
    info.update(band_lo=int(lo), band_hi=int(hi),
                middle_rows=f"{k0 + 1}-{k0 + len(blk)} of {len(det['plot_rows'])}")
    if src_mm:
        info["band_mm"] = round((hi - lo) * src_mm, 1)
    across_len = H if across_axis == 0 else W
    along_len = W if across_axis == 0 else H
    if band_tile:
        side = across = int(round(hi - lo))
    else:
        side = int(win_px)                             # along
        across = int(win_across or win_px)
    if src_mm:
        info["tile_across_mm"] = round(across * src_mm, 1)
    info["rows_covered"] = bool(hi - lo <= across + 1)
    if strict and (across > across_len or side > along_len):
        info["reason"] = (f"plot is {across_len} px across, one tile is {across} px"
                          if across > across_len else
                          f"plot is {along_len} px long, one tile is {side} px")
        return [], info
    if not strict:
        side = max(8, min(side, across_len))
        across = side if (band_tile or not win_across) else max(8, min(across, across_len))
    stride = side if band_tile else max(1, side - max(0, overlap))
    # the rows may be tilted in a raw (unrectified) crop: the band centre moves
    # across by tan(angle) per pixel along (angle measured in the reduced copy)
    t = math.tan(math.radians(det["angle_deg"]))
    lc_full = det["along_centre"] * fl
    centre0 = (lo + hi) / 2.0
    # along the plot: only the usable length, tiles centred in it
    us = _usable_length(small, across_axis, det, band,
                        (src_mm * fl) if src_mm else None, ext_margin_mm, end_margin_mm)
    qa["us"] = us
    u0, u1 = us["u0"] * fl, us["u1"] * fl
    span = u1 - u0
    if band_tile:
        n_fit = int(span // side)
        if n_fit == 0 and span >= 0.6 * side:   # one tile, centred, ends just beyond
            n_fit = 1
        first = u0 + (span - n_fit * side) / 2.0
        starts = [int(round(first + i * side)) for i in range(n_fit)]
    elif strict:                                # ground tiles: end to end, centred
        starts = []
        if span >= side:
            nt = int((span - side) // stride) + 1
            first = u0 + (span - (side + (nt - 1) * stride)) / 2.0
            starts = [int(round(first + i * stride)) for i in range(nt)]
    else:
        starts = _starts(int(round(u0)), int(round(u1)), side, stride) if span >= side else []
    mm = src_mm or 0.0
    info.update(plot_len_m=round((us["p1"] - us["p0"]) * fl * mm / 1000.0, 2) if mm else "",
                usable_lo=int(u0), usable_hi=int(u1),
                usable_len_m=round(max(0.0, span) * mm / 1000.0, 2) if mm else "",
                end_trim=us["trim"], tiles_fit=len(starts))
    cand = []
    for c_i, l0 in enumerate(starts):
        l0 = max(0, min(l0, along_len - side))
        mid_along = l0 + side / 2.0
        a_c = centre0 + t * (mid_along - lc_full) * (f / fl)
        a0 = int(round(a_c - across / 2.0))
        a0 = max(0, min(a0, across_len - across))
        x, y, w, h = ((l0, a0, side, across) if across_axis == 0
                      else (a0, l0, across, side))
        cand.append((c_i, 0, x, y, w, h))
    if not cand:
        info["reason"] = (f"usable length {max(0.0, span) * mm / 1000.0:.2f} m is shorter "
                          "than one tile" if mm else "usable length shorter than one tile")
    return cand, info


def _qa_font(size):
    from PIL import ImageFont
    for name in ("arialbd.ttf", "arial.ttf", "DejaVuSans-Bold.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except Exception:
            pass
    return ImageFont.load_default()


def _write_middle_qa(qa, W, H, across_axis, cand_all, cand_kept, note=""):
    """QA picture of one plot: rows (green), middle band (yellow), plot ends (grey),
    usable length limits (orange), candidate tiles (thin blue) and the tiles written
    (thick cyan, numbered). A plot without tiles still gets its picture, with the
    reason written on it. Also used for the plain grid (no rows: det None)."""
    import row_detect as rd
    from PIL import Image, ImageDraw
    small, det, band, us, fl = qa["small"], qa["det"], qa["band"], qa["us"], qa["fl"]
    rgb = np.transpose(small, (1, 2, 0))
    im = (rd.draw_rows(rgb, det, band=band) if det is not None
          else Image.fromarray(np.ascontiguousarray(rgb)))
    if det is None and max(im.size) > 1600:
        k = 1600.0 / max(im.size)
        im = im.resize((max(1, int(im.size[0] * k)), max(1, int(im.size[1] * k))))
    sc = im.size[0] / float(W)
    dr = ImageDraw.Draw(im)

    def across_line(v, col, width):
        v = v * fl * sc
        if across_axis == 0:
            dr.line([(v, 0), (v, im.size[1])], fill=col, width=width)
        else:
            dr.line([(0, v), (im.size[0], v)], fill=col, width=width)
    if us:
        across_line(us["p0"], (150, 150, 150), 2)
        across_line(us["p1"], (150, 150, 150), 2)
        across_line(us["u0"], (255, 140, 0), 3)
        across_line(us["u1"], (255, 140, 0), 3)
    if qa.get("lohi"):                                   # plain grid: the tiled band
        lo, hi = qa["lohi"]
        for v in (lo, hi):
            v = v * sc
            if across_axis == 0:
                dr.line([(0, v), (im.size[0], v)], fill=(255, 210, 0), width=2)
            else:
                dr.line([(v, 0), (v, im.size[1])], fill=(255, 210, 0), width=2)
    kept = set((c[2], c[3]) for c in cand_kept)
    for (_c, _r, x, y, w, h) in cand_all:
        if (x, y) not in kept:
            dr.rectangle([x * sc, y * sc, (x + w) * sc, (y + h) * sc],
                         outline=(0, 120, 255), width=1)
    fnt = _qa_font(max(12, im.size[1] // 30))
    for k, (_c, _r, x, y, w, h) in enumerate(cand_kept, 1):
        dr.rectangle([x * sc, y * sc, (x + w) * sc, (y + h) * sc],
                     outline=(0, 230, 255), width=4)
        dr.text((x * sc + 6, y * sc + 4), str(k), fill=(0, 230, 255), font=fnt)
    if note:
        fs = _qa_font(max(14, im.size[1] // 22))
        dr.rectangle([0, 0, im.size[0], fs.size + 12 if hasattr(fs, "size") else 30],
                     fill=(0, 0, 0))
        dr.text((8, 4), note, fill=(255, 80, 80), font=fs)
    os.makedirs(os.path.dirname(qa["path"]), exist_ok=True)
    im.save(qa["path"], quality=88)


def _save_rows(s, W, H, across_axis, det_job, rows_dir, stem, ext, bits="same",
               compress="lzw", select=None, horizontal=False, names=None):
    """Write each of the plot's crop rows as its own image: <stem>_row<k><ext>, k = 1..n
    across the plot (from image top / left). The image is the axis-aligned box around
    the row's band (row centre +- half a pitch, following the rows' tilt), full plot
    length; pixels outside the band are set to 0, everything inside is the crop's own
    pixels (no resampling - the row strip of a raw crop is original camera pixels).
    select: row numbers k to write (None = all). horizontal: turn the image 90 deg
    when the rows run up-down in the crop, so every row image runs left-right
    (lossless). names: k -> file name (default <stem>_row<k><ext>).
    Returns one index record per row."""
    import rasterio.windows
    from PIL import Image
    det, f, fl = det_job[:3]
    blk = set(det_job[3]) if len(det_job) > 3 else set()
    across_len = H if across_axis == 0 else W
    along_len = W if across_axis == 0 else H
    P = det["pitch_px"] * f
    t = math.tan(math.radians(det["angle_deg"])) * (f / fl)
    lc = det["along_centre"] * fl
    os.makedirs(rows_dir, exist_ok=True)
    out = []
    al = np.arange(along_len, dtype=np.float64)
    for k, ri in enumerate(det["plot_rows"], 1):
        if select is not None and k not in select:
            continue
        c = det["rows"][ri] * f
        lo_l = c - P / 2.0 + t * (al - lc)               # band edges along the plot
        hi_l = c + P / 2.0 + t * (al - lc)
        a0 = int(max(0, math.floor(lo_l.min())))
        a1 = int(min(across_len, math.ceil(hi_l.max())))
        if a1 - a0 < 4:
            continue
        if across_axis == 0:
            win = rasterio.windows.Window(0, a0, W, a1 - a0)
        else:
            win = rasterio.windows.Window(a0, 0, a1 - a0, H)
        a = s.read([1, 2, 3], window=win)              # 3 x h x w
        acr = np.arange(a0, a1, dtype=np.float64)
        if across_axis == 0:                            # across = image y
            keep = (acr[:, None] >= lo_l[None, :]) & (acr[:, None] < hi_l[None, :])
        else:                                           # across = image x
            keep = (acr[None, :] >= lo_l[:, None]) & (acr[None, :] < hi_l[:, None])
        a = convert_bits(np.where(keep[None], a, 0).astype(a.dtype),
                         8 if ext == ".jpg" else bits)
        turned = bool(horizontal and across_axis == 1)
        if turned:
            a = np.ascontiguousarray(np.rot90(a, 1, axes=(1, 2)))
        name = (names or {}).get(k) or f"{stem}_row{k}{ext}"
        write_tile(os.path.join(rows_dir, name), a, ext, compress)
        x0, y0 = (0, a0) if across_axis == 0 else (a0, 0)
        out.append(dict(file=name, plot_stem=stem, row=k, of_rows=len(det["plot_rows"]),
                        middle=ri in blk, present=bool(det["present"][ri]),
                        contrast=round(float(det["contrast"][ri]), 3),
                        centre_px=round(c, 1), width_px=round(P, 1),
                        angle_deg=round(det["angle_deg"], 2),
                        x0=x0, y0=y0, w=(a.shape[1] if turned else a.shape[2]),
                        h=(a.shape[2] if turned else a.shape[1]), across_axis=across_axis,
                        turned_90=turned))
    return out


def make_tiles(src_dir, out_dir, tile=1024, overlap=0, rows_only=True,
               row_pitch_mm=0.0, plot_width_mm=0.0, fmt="png",
               min_veg_pct=0.0, keep_empty_frac=0.15, max_tiles_per_plot=0,
               seed=7, target_mm_px=0.0, per_plot_dirs=False,
               middle_rows=0, rows_per_plot=0, band_tile=True, write_qa=True,
               end_margin_m=0.5, tiles_pick="random", save_rows=False,
               bits="same", compress="lzw", georef=True, link_source=True,
               tile_m=0.0, overlap_m=0.0, cover_rows=False, row_buffer_m=0.05,
               tile_along_m=0.0,
               log=print, progress=None, cancel=None):
    """Cut every crop in `src_dir` into tiles.

    tile_m / overlap_m    > 0: GROUND tiles - every tile covers tile_m x tile_m metres
                          and has the same pixel size on every plot: the window is
                          tile_m / (plot mm/px) source pixels, written at the output
                          scale target_mm_px (0 = the median mm/px of the plots), i.e.
                          tile_m / target px. Plots at that scale are copied exactly,
                          others resampled (Lanczos). Only whole tiles are written.
                          tile / overlap / band_tile are then derived / ignored.
    cover_rows            middle-rows mode: the tile's size ACROSS the rows is set so
                          it holds the middle_rows rows: middle_rows x the median row
                          spacing of all plots (measured first) + row_buffer_m on each
                          side - one size for the whole run. Along the rows the tile
                          is tile_along_m (0 = the same, square). Non-square tiles of
                          plots whose rows run up-down are turned 90 deg, so every
                          tile has the same width x height. rows_report.csv says per
                          plot whether its band fits in the tile (rows_covered).
    tile / overlap        pixels. stride = tile - overlap.
    rows_only             restrict to the measured seed-row band.
    row_pitch_mm          row spacing; only used to find the band. 0 = guess from
                          plot_width_mm and a nominal 5 rows.
    min_veg_pct           drop tiles below this vegetation fraction, EXCEPT a
                          random `keep_empty_frac` of them, kept as hard
                          negatives - bare soil and stubble is what a seedling
                          detector gets wrong, so it has to see some.
    max_tiles_per_plot    tiles per plot, every mode. 0 = all. Otherwise picked by
                          tiles_pick; on the plain grid 'random' back-fills tiles
                          dropped as nodata / low vegetation.
    per_plot_dirs         write each plot's tiles into its own subfolder
                          out_dir/<plot stem>/ instead of one flat folder.
    middle_rows           > 0: tile only the band holding this many CENTRAL crop
                          rows of each plot (rows detected per plot, as in
                          PhenoFuse); rows_only is ignored. With
                          max_tiles_per_plot, that many tiles are drawn at random
                          from the band.
    rows_per_plot         sown rows per plot, if known (0 = the run of rows found).
    band_tile             middle-rows mode: square tiles as wide as the band
                          (native pixels, end to end); False = `tile` px tiles.
    write_qa              middle-rows mode: out_dir/_rows_qa/<plot>_rows.jpg with
                          the detected rows, the band and the tiles, plus
                          rows_report.csv.
    end_margin_m          middle-rows mode: keep tiles this far (metres) inside the
                          plot ends (after trimming to the planted length); the
                          tiles are centred in what is left.
    tiles_pick            with max_tiles_per_plot: 'random' or 'spread' (evenly
                          spaced along the plot).
    save_rows             middle-rows mode: also write every detected crop row of
                          each plot as its own image, out_dir/Rows/<plot>/
                          <plot>_row<k>.<ext> (original crop pixels, outside the
                          row masked to 0), indexed in Rows/rows_index.csv.
    bits                  'same' (default) = the crop's own values, an exact copy;
                          8 or 16 converts (16->8 v/257 rounded, 8->16 v*257).
    compress              TIFF compression: 'lzw' (lossless, smaller) or 'none'.
                          Both hold identical pixels.
    georef                TIFF tiles are GeoTIFFs: the crop's CRS and transform
                          shifted to the tile window (or its GCPs moved into it).
    link_source           tile_index.csv columns and GeoTIFF tags tying each tile
                          to its source photo / IIQ and its pixel box there.

    Returns (n_tiles, n_plots, index_path).
    """
    import rasterio
    from PIL import Image

    make_tiles.last_run = {}
    ext = {"png": ".png", "jpeg": ".jpg", "jpg": ".jpg", "tif": ".tif"}.get(fmt.lower(), ".png")
    bits = str(bits).strip().lower() or "same"
    compress = "none" if str(compress).lower().startswith(("none", "unc")) else "lzw"
    if ext == ".jpg" and bits == "16":
        raise ValueError("JPEG is 8-bit only: choose TIFF or PNG for 16-bit tiles.")
    if ext != ".tif" and georef:
        log("[i] georeference is written into TIFF tiles only (PNG/JPEG carry none); "
            "tile_index.csv still maps every tile to its plot")
    links, _iiq_run, iiq_settings = (_source_links(src_dir) if link_source
                                     else ({}, {}, ""))
    os.makedirs(out_dir, exist_ok=True)
    files = sorted(f for f in os.listdir(src_dir)
                   if f.lower().endswith((".tif", ".tiff")))
    if not files:
        log(f"[!] no GeoTIFF crops in {src_dir}")
        return 0, 0, ""
    cover_rows = bool(cover_rows and middle_rows > 0)
    uniform = tile_m > 0 or cover_rows
    out_along = out_across = 0                 # uniform: output tile px along / across
    det_cache = {}
    if uniform:
        # ground tiles: one output scale for the whole run, so every tile has the
        # same pixel size and the same ground size
        gsd = {}
        for fn in files:
            try:
                with rasterio.open(os.path.join(src_dir, fn)) as s0:
                    m = _plot_mm_px(s0)
                    if m is None and plot_width_mm > 0:
                        m = plot_width_mm / float(min(s0.width, s0.height))
                    if m:
                        gsd[fn] = m
            except Exception:
                pass
        if not gsd:
            log("[!] the plot images carry no ground scale (no georeference) and no "
                "plot width was given: tiles in metres cannot be cut")
            return 0, 0, ""
        miss = [f for f in files if f not in gsd]
        if miss:
            log(f"[!] {len(miss)} plot image(s) have no ground scale and are skipped: "
                + ", ".join(miss[:5]) + (" ..." if len(miss) > 5 else ""))
            files = [f for f in files if f in gsd]
        v = np.array(list(gsd.values()))
        if target_mm_px <= 0:                  # the plot nearest the median: its
            target_mm_px = float(v[np.argmin(abs(v - np.median(v)))])   # tiles are exact
        along_m = across_m = tile_m
        if cover_rows:
            # rows of every plot first: one tile size across, from the median spacing
            log(f"Measuring the rows of {len(files)} plot(s) to size the tiles ...")
            pitches = []
            for fi, fn in enumerate(files):
                if cancel is not None and cancel.is_set():
                    log("Tiling cancelled.")
                    return 0, 0, ""
                try:
                    with rasterio.open(os.path.join(src_dir, fn)) as s0:
                        ax0 = 0 if s0.height < s0.width else 1
                        _sm, det0, f0, fl0 = _detect_plot(s0, s0.width, s0.height, ax0,
                                                          gsd[fn], row_pitch_mm,
                                                          rows_per_plot)
                    det_cache[fn] = (det0, f0, fl0)
                    if det0.get("ok"):
                        pitches.append(det0["pitch_px"] * f0 * gsd[fn])
                except Exception as e:
                    log(f"[!] {fn}: {e}")
                if progress:
                    progress(fi + 1, 2 * len(files))
            if pitches:
                pitch_mm = float(np.median(pitches))
                log(f"Row spacing {min(pitches):.0f}-{max(pitches):.0f} mm on "
                    f"{len(pitches)} plot(s), median {pitch_mm:.0f} mm.")
            elif row_pitch_mm > 0:
                pitch_mm = row_pitch_mm
                log(f"[!] rows not found on any plot: using the given spacing {pitch_mm:g} mm")
            else:
                log("[!] rows not found on any plot and no row spacing given - "
                    "nothing written")
                return 0, 0, ""
            across_m = round(middle_rows * pitch_mm / 1000.0 + 2 * max(0.0, row_buffer_m), 3)
            along_m = tile_along_m if tile_along_m > 0 else across_m
            tile_m = along_m
        out_along = max(8, int(round(along_m * 1000.0 / target_mm_px)))
        out_across = max(8, int(round(across_m * 1000.0 / target_mm_px)))
        tile = out_along
        overlap = max(0, int(round(overlap_m * 1000.0 / target_mm_px)))
        if overlap >= tile:
            raise ValueError("Overlap must be smaller than the tile size.")
        band_tile = False
        if cover_rows:
            log(f"Tiles cover the {middle_rows} middle rows: across {middle_rows} x "
                f"{pitch_mm:.0f} mm + 2 x {row_buffer_m * 1000:.0f} mm buffer = "
                f"{across_m:.3f} m, along {along_m:.3f} m"
                + (" (square)" if along_m == across_m else "") + ".")
        log(f"Ground tiles: {along_m:g} m along x {across_m:g} m across, overlap "
            f"{overlap_m:g} m. Plot images "
            f"{v.min():.3f}-{v.max():.3f} mm/px (median {np.median(v):.3f}); output "
            f"{target_mm_px:.3f} mm/px -> every tile {out_along} x {out_across} px"
            + (" (rows left-right; tiles of plots with up-down rows are turned)"
               if out_along != out_across else "") + ".")
    stride = max(1, tile - max(0, overlap))
    rnd = random.Random(seed)
    rows_out = []
    row_reports = []
    row_files = []
    # extraction margin per crop (manifest.csv next to RawCrops): the plot's own
    # ends lie that far inside the crop ends
    ext_margin = {}
    for mf in (os.path.join(os.path.dirname(os.path.abspath(src_dir)), "manifest.csv"),
               os.path.join(src_dir, "manifest.csv")):
        if os.path.exists(mf):
            try:
                with open(mf, newline="", encoding="utf-8") as fh:
                    for rr in csv.DictReader(fh):
                        try:
                            ext_margin[rr.get("raw_file", "")] = float(rr.get("margin_m") or 0)
                        except ValueError:
                            pass
            except Exception:
                pass
            break
    if middle_rows > 0:
        log(f"Tiling {len(files)} plot(s): the {middle_rows} middle crop rows of each "
            f"plot (rows detected per plot"
            + (f", {rows_per_plot} sown rows" if rows_per_plot else "") + "), "
            + ("square tiles spanning the band, native pixels" if band_tile
               else (f"{tile_m:g} m x {out_across * target_mm_px / 1000.0:.3g} m tiles "
                     "centred on the band") if uniform
               else f"{tile}px tiles on the band")
            + f", tiles kept {end_margin_m:.2f} m inside the planted plot ends"
            + (f", {max_tiles_per_plot} tile(s) per plot ({tiles_pick})" if max_tiles_per_plot
               else ", every tile that fits")
            + f", {ext[1:].upper()} output.")
        _out_desc(log, ext, bits, compress, georef, link_source)
    else:
        log(f"Tiling {len(files)} plot(s): "
            + (f"{tile_m:g} m tiles, {overlap_m:g} m overlap " if uniform else
               f"{tile}px tiles, {overlap}px overlap ")
            + f"(stride {stride} px), {'row band only' if rows_only else 'whole crop'}"
            + (f", {max_tiles_per_plot} tile(s) per plot ({tiles_pick})" if max_tiles_per_plot
               else ", every tile")
            + f", {ext[1:].upper()} output.")
        _out_desc(log, ext, bits, compress, georef, link_source)

    det_cache_used = bool(det_cache)
    for fi, fn in enumerate(files):
        if cancel is not None and cancel.is_set():
            log("Tiling cancelled.")
            break
        stem = os.path.splitext(fn)[0]
        path = os.path.join(src_dir, fn)
        try:
            with rasterio.open(path) as s:
                W, H = s.width, s.height
                crop_tags = s.tags() if link_source else {}
                # rows run along the long axis; 'across' is the short one
                across_axis = 0 if H < W else 1
                # Ground-referenced tiling. Ground sampling varies plot to plot
                # (0.274-0.470 mm/px on OZ Barley), so a fixed pixel tile covers a
                # different amount of ground on every plot - 2.2 rows on one, 3.8 on
                # another. Reading a window sized in GROUND units and resampling it
                # to `tile` px makes every tile the same real size, which is what a
                # "cover ~3 rows" requirement actually asks for.
                src_mm = _plot_mm_px(s) if (target_mm_px > 0 or middle_rows > 0) else None
                if src_mm is None and middle_rows > 0 and plot_width_mm > 0:
                    src_mm = plot_width_mm / float(H if across_axis == 0 else W)
                if uniform:
                    src_mm = gsd[fn]
                win_px = tile
                win_across = None
                if target_mm_px > 0 and src_mm and src_mm > 0:
                    win_px = max(8, int(round(tile * target_mm_px / src_mm)))
                    if uniform and out_across != out_along:
                        win_across = max(8, int(round(out_across * target_mm_px / src_mm)))
                if row_pitch_mm > 0 and plot_width_mm > 0:
                    across_px = H if across_axis == 0 else W
                    pitch_px = row_pitch_mm / (plot_width_mm / across_px)
                else:
                    pitch_px = (H if across_axis == 0 else W) / 5.5
                lo, hi = 0, (H if across_axis == 0 else W)
                qa_job, qa_note, cand_all = None, "", []
                if middle_rows > 0:
                    qa = (os.path.join(out_dir, "_rows_qa", stem + "_rows.jpg")
                          if write_qa else None)
                    det_in = None
                    if fn in det_cache:
                        det_in = _detect_plot(s, W, H, across_axis, src_mm, row_pitch_mm,
                                              rows_per_plot, det=det_cache.pop(fn)[0])
                    cand, info = _middle_row_cands(
                        s, W, H, across_axis, middle_rows, rows_per_plot, src_mm,
                        row_pitch_mm, band_tile, win_px,
                        int(round(overlap * win_px / max(1, tile))), qa,
                        ext_margin_mm=1000.0 * ext_margin.get(fn, 0.0),
                        end_margin_mm=1000.0 * max(0.0, end_margin_m), strict=uniform,
                        win_across=win_across, det_in=det_in)
                    qa_job = info.pop("_qa", None)
                    det_job = info.pop("_det", None)
                    if save_rows and det_job:
                        recs = _save_rows(s, W, H, across_axis, det_job,
                                          os.path.join(out_dir, "Rows", stem), stem, ext,
                                          bits, compress)
                        for rc in recs:
                            rc["rel_path"] = stem + "/" + rc["file"]
                        row_files.extend(recs)
                        info["rows_saved"] = len(recs)
                    info = dict(plot_stem=stem, src=fn, tiles_in_band=len(cand), **info)
                    cand_all = cand
                    if not cand:
                        qa_note = f"no tiles: {info['reason']}"
                        log(f"[!] {fn}: no middle-row tiles ({info['reason']}) - skipped")
                    else:
                        cand = _pick(cand, max_tiles_per_plot, tiles_pick, rnd)
                        if max_tiles_per_plot and len(cand) < max_tiles_per_plot:
                            log(f"   {fn}: only {len(cand)} tile(s) fit, "
                                f"{max_tiles_per_plot} asked")
                        if info.get("rows_covered") is False:
                            qa_note = (f"band {info.get('band_mm', '?')} mm is wider than "
                                       f"the tile ({info.get('tile_across_mm', '?')} mm)")
                    info["tiles_drawn"] = len(cand)
                    row_reports.append(info)
                    if cand:
                        lo = min(c[3] if across_axis == 0 else c[2] for c in cand)
                        hi = max(c[3] + c[5] if across_axis == 0 else c[2] + c[4]
                                 for c in cand)
                else:
                    small = None
                    if rows_only:
                        dec = max(1, int(round(min(W, H) / 900)))
                        small = to8(s.read([1, 2, 3], out_shape=(3, max(1, H // dec), max(1, W // dec))))
                        l2, h2 = row_band(small, across_axis, pitch_px / dec)
                        lo, hi = l2 * dec, min(hi, h2 * dec)
                        if hi - lo < win_px:    # band thinner than a tile: centre it
                            c = (lo + hi) // 2
                            lo = max(0, c - win_px // 2)
                            hi = min((H if across_axis == 0 else W), lo + win_px)
                    if write_qa:
                        if small is None:
                            small = _small_copy(s, W, H)
                        qa_job = dict(path=os.path.join(out_dir, "_tiles_qa",
                                                        stem + "_tiles.jpg"),
                                      small=small, det=None, band=None, us=None, fl=1.0,
                                      lohi=(lo, hi))

                limit = 0
                if middle_rows <= 0:
                    y0b, y1b = (lo, hi) if across_axis == 0 else (0, H)
                    x0b, x1b = (0, W) if across_axis == 0 else (lo, hi)
                    # grid is laid out in SOURCE pixels; win_px == tile unless the
                    # tiles are ground-referenced, in which case it is the source
                    # window that resamples down to `tile`.
                    src_stride = max(1, int(round(stride * win_px / tile)))
                    cand = []
                    for r_i, y in enumerate(_starts(y0b, y1b, win_px, src_stride)):
                        for c_i, x in enumerate(_starts(x0b, x1b, win_px, src_stride)):
                            cand.append((c_i, r_i, x, y, win_px, win_px))
                    if not cand and not uniform:
                        cand = [(0, 0, x0b, y0b, win_px, win_px)]
                    if uniform:                 # whole tiles only
                        cand = [c for c in cand if c[2] + win_px <= W and c[3] + win_px <= H]
                        if not cand:
                            qa_note = f"no tiles: {W}x{H} px is smaller than one tile"
                            log(f"[!] {fn}: {W}x{H} px is smaller than one "
                                f"{tile_m:g} m tile ({win_px} px) - skipped")
                    cand_all = list(cand)
                    if max_tiles_per_plot and tiles_pick != "spread":
                        # random: walk a shuffled list until N tiles are kept, so
                        # tiles dropped as nodata / low vegetation are back-filled
                        rnd.shuffle(cand)
                        limit = max_tiles_per_plot
                    else:
                        cand = _pick(cand, max_tiles_per_plot, tiles_pick, rnd)

                kept = 0
                written = []
                for (c_i, r_i, x, y, sw, sh) in cand:
                    if middle_rows <= 0 and limit and kept >= limit:
                        break
                    # output size: ground tiles are out_along x out_across px on
                    # every plot; a band-spanning middle-rows tile keeps its native
                    # size unless a ground scale is fixed; else `tile`
                    if uniform:
                        fw, fh = ((out_along, out_across) if across_axis == 0
                                  else (out_across, out_along))
                    elif middle_rows > 0 and band_tile:
                        fw = fh = (max(8, int(round(sw * src_mm / target_mm_px)))
                                   if (target_mm_px > 0 and src_mm) else sw)
                    else:
                        fw = fh = tile
                    w = min(sw, W - x); h = min(sh, H - y)
                    if w < sw // 2 or h < sh // 2:
                        continue
                    if uniform and (w < sw or h < sh):
                        continue
                    a = s.read(window=rasterio.windows.Window(x, y, w, h))
                    if a.shape[0] < 3:
                        continue
                    rgb = a[:3]
                    v8 = to8(rgb)                     # analysis only; rgb is written
                    valid = v8.max(axis=0) > 4
                    if valid.mean() < 0.5:            # mostly nodata
                        continue
                    R, G, B = (v8[i].astype(np.float32) for i in range(3))
                    exg = 2 * G - R - B
                    veg = float((exg[valid] > 18).mean() * 100) if valid.any() else 0.0
                    if min_veg_pct > 0 and veg < min_veg_pct:
                        if rnd.random() > keep_empty_frac:
                            continue                   # dropped, but some kept
                    name = f"{stem}_tile_c{c_i}_r{r_i}{ext}"
                    if per_plot_dirs:
                        sub = os.path.join(out_dir, stem)
                        os.makedirs(sub, exist_ok=True)
                        dest = os.path.join(sub, name)
                        rel = stem + "/" + name
                    else:
                        dest = os.path.join(out_dir, name)
                        rel = name
                    out = convert_bits(rgb, 8 if ext == ".jpg" else bits)
                    # proportional, so a clipped window keeps the same scale; each
                    # channel resampled in float (works for 8 and 16 bit)
                    ow = max(1, int(round(w * fw / sw)))
                    oh = max(1, int(round(h * fh / sh)))
                    resampled = (ow, oh) != (w, h)
                    if resampled:
                        out = _resize(out, ow, oh)
                    crs = tf = gc = None
                    if georef and ext == ".tif":
                        crs, tf, gc = _tile_georef(s, x, y, w, h, ow, oh)
                    # non-square ground tiles all run left-right (rows along x)
                    turned = bool(uniform and across_axis == 1 and out_along != out_across)
                    if turned:
                        out = np.ascontiguousarray(np.rot90(out, 1, axes=(1, 2)))
                        if tf is not None or gc:
                            crs, tf, gc = _turn_georef(crs, tf, gc, ow)
                    lk = (_link_fields(links.get(fn), crop_tags, x, y, w, h, iiq_settings)
                          if link_source else {})
                    tags = None
                    if link_source and ext == ".tif":
                        tags = dict(PLOT_IMAGE=fn, PLOT_PIXEL_BOX=f"{x},{y},{w},{h}",
                                    SOURCE_IMAGE=lk.get("source_image", ""),
                                    PIXEL_SOURCE=lk.get("pixel_source", ""),
                                    IIQ_SETTINGS=os.path.basename(lk.get("iiq_settings", "")),
                                    SOURCE_PIXEL_BOX=(
                                        f"{lk['src_photo_x0']},{lk['src_photo_y0']},"
                                        f"{lk['src_photo_w']},{lk['src_photo_h']}"
                                        if lk.get("src_photo_x0", "") != "" else ""),
                                    BIT_DEPTH=f"{8 * out.dtype.itemsize}"
                                              + ("" if out.dtype == rgb.dtype else
                                                 f" (converted from {8 * rgb.dtype.itemsize})"))
                    if tags is not None and turned:
                        tags["TURNED_90"] = "True"
                    write_tile(dest, out, ext, compress, crs, tf, gc, tags)
                    written.append((c_i, r_i, x, y, w, h))
                    rows_out.append(dict(tile=name, rel_path=rel, plot_stem=stem, col=c_i, row=r_i,
                                         x0=x, y0=y, w=w, h=h,
                                         out_px=(out.shape[2] if out.shape[2] == out.shape[1]
                                                 else f"{out.shape[2]}x{out.shape[1]}"),
                                         out_w=out.shape[2], out_h=out.shape[1],
                                         turned_90=turned,
                                         src_mm_px=round(src_mm, 4) if src_mm else "",
                                         tile_mm_px=round(target_mm_px, 4) if (src_mm and target_mm_px) else "",
                                         tile_m=tile_m if uniform else "",
                                         tile_across_m=(round(out_across * target_mm_px / 1000.0, 4)
                                                        if uniform else ""),
                                         resampled=resampled,
                                         ground_mm=round(w * (src_mm or 0), 1) if src_mm else "",
                                         band_lo=lo, band_hi=hi,
                                         across_axis=across_axis,
                                         src_w=W, src_h=H, veg_pct=round(veg, 2),
                                         src=fn,
                                         bits=8 * out.dtype.itemsize,
                                         src_bits=8 * rgb.dtype.itemsize,
                                         georef=("affine" if tf is not None else
                                                 "gcps" if gc else ""),
                                         **lk))
                    kept += 1
                if qa_job:
                    try:
                        if cand_all and not written and not qa_note:
                            qa_note = "no tiles kept (nodata / low vegetation)"
                        _write_middle_qa(qa_job, W, H, across_axis, cand_all, written,
                                         note=qa_note)
                    except Exception as e:
                        log(f"[!] {fn}: QA picture: {e}")
        except Exception as e:
            log(f"[!] {fn}: {e}")
        if progress:
            progress((len(files) + fi + 1) if det_cache_used else (fi + 1),
                     (2 * len(files)) if det_cache_used else len(files))
        if (fi + 1) % 25 == 0:
            log(f"   {fi+1}/{len(files)} plots, {len(rows_out)} tiles")

    idx = os.path.join(out_dir, "tile_index.csv")
    keys = ["tile", "rel_path", "plot_stem", "col", "row", "x0", "y0", "w", "h", "out_px",
            "out_w", "out_h", "turned_90",
            "src_mm_px", "tile_mm_px", "ground_mm", "tile_m", "tile_across_m", "resampled",
            "band_lo", "band_hi", "across_axis", "src_w", "src_h", "veg_pct", "src",
            "bits", "src_bits", "georef"]
    if link_source:
        keys += ["source_image", "pixel_source", "src_photo_x0", "src_photo_y0",
                 "src_photo_w", "src_photo_h", "iiq_settings", "link_note"]
    with open(idx, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        wr.writeheader()
        for rr in rows_out:
            wr.writerow(rr)
    if row_files:
        ri = os.path.join(out_dir, "Rows", "rows_index.csv")
        with open(ri, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(row_files[0].keys()))
            wr.writeheader()
            wr.writerows(row_files)
        log(f"Saved {len(row_files)} row images from "
            f"{len({r['plot_stem'] for r in row_files})} plot(s) -> {os.path.dirname(ri)}")
    if row_reports:
        rep = os.path.join(out_dir, "rows_report.csv")
        rkeys = ["plot_stem", "src", "ok", "reason", "rows_found", "rows_present",
                 "pitch_px", "pitch_mm", "angle_deg", "middle_rows", "band_lo",
                 "band_hi", "band_mm", "tile_across_mm", "rows_covered",
                 "plot_len_m", "end_trim", "usable_lo", "usable_hi",
                 "usable_len_m", "tiles_fit", "tiles_in_band", "tiles_drawn",
                 "rows_saved", "veg_mode"]
        with open(rep, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=rkeys, extrasaction="ignore")
            wr.writeheader()
            for rr in row_reports:
                wr.writerow(rr)
        n_ok = sum(1 for rr in row_reports if rr["tiles_in_band"])
        log(f"Middle rows found on {n_ok}/{len(row_reports)} plot(s); per-plot rows: "
            f"{rep}" + (f"  (QA pictures of all {len(row_reports)} plots in _rows_qa\\)"
                        if write_qa else ""))
        unc = [rr["plot_stem"] for rr in row_reports
               if rr.get("tiles_drawn") and rr.get("rows_covered") is False]
        if unc:
            log(f"[!] on {len(unc)} plot(s) the {middle_rows} rows are wider than the tile "
                f"(wider row spacing than the median): " + ", ".join(unc[:6])
                + (" ..." if len(unc) > 6 else "")
                + " - raise the buffer to cover them (rows_report.csv: band_mm, "
                  "tile_across_mm)")
        elif cover_rows:
            log(f"Rows covered: the {middle_rows} middle rows fit inside the tile on every "
                "plot with tiles.")
    elif write_qa:
        log(f"QA pictures: {os.path.join(out_dir, '_tiles_qa')}")
    plots = len({r["plot_stem"] for r in rows_out})
    log(f"Wrote {len(rows_out)} tiles from {plots} plot(s) -> {out_dir}")
    log(f"Index: {idx}  (needed to map a tile back to its plot)")
    if rows_out:
        v = np.array([r["veg_pct"] for r in rows_out])
        log(f"Vegetation per tile: median {np.median(v):.2f}%  "
            f"empty (<0.5%) {100*(v<0.5).mean():.1f}%")
        if uniform:
            sz = {f"{r['out_w']}x{r['out_h']}" for r in rows_out}
            log(f"Tile size check: {len(rows_out)} tiles, pixel sizes {sorted(sz)} "
                f"({'all the same' if len(sz) == 1 else 'NOT uniform'}); "
                f"{sum(1 for r in rows_out if not r['resampled'])} exact copies, "
                f"{sum(1 for r in rows_out if r['resampled'])} resampled to "
                f"{target_mm_px:.3f} mm/px")
    # what the run actually used, for the recipe
    make_tiles.last_run = (dict(tile_px=[out_along, out_across],
                                tile_along_m=round(out_along * target_mm_px / 1000.0, 4),
                                tile_across_m=round(out_across * target_mm_px / 1000.0, 4),
                                out_mm_px=round(target_mm_px, 4))
                           if uniform else {})
    return len(rows_out), plots, idx


# ---------------------------------------------------------------------------
# whole rows and fixed-size tiles along each row (tab 6b)
# ---------------------------------------------------------------------------
def parse_rows(text):
    """Which crop rows to use: 'all' -> ('all', None); 'middle 3' -> ('middle', 3);
    '1,2,5-7' -> ('list', [1, 2, 5, 6, 7]). Rows are numbered 1..n across the plot
    from the image top (rows running left-right) or left (rows running up-down)."""
    t = str(text or "").strip().lower()
    if t in ("", "all", "*"):
        return "all", None
    if t.startswith("middle"):
        try:
            return "middle", max(1, int(t.split()[-1]))
        except ValueError:
            raise ValueError(f"rows '{text}': write e.g. 'middle 3'")
    ks = set()
    for part in t.replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            ks.update(range(int(a), int(b) + 1))
        else:
            ks.add(int(part))
    if not ks or min(ks) < 1:
        raise ValueError(f"rows '{text}': use 'all', 'middle N' or numbers like 1,2,3")
    return "list", sorted(ks)


def _resize(a, ow, oh):
    """Lanczos resample of a 3 x h x w array to ow x oh, per channel in float, same
    dtype (8 or 16 bit)."""
    from PIL import Image
    top = 255 if a.dtype == np.uint8 else 65535
    return np.stack([np.asarray(Image.fromarray(c.astype(np.float32), "F")
                                .resize((ow, oh), Image.LANCZOS))
                     for c in a]).round().clip(0, top).astype(a.dtype)


def _turn_georef(crs, tf, gc, ow):
    """Georeference of a tile (ow px wide before turning) after np.rot90(.., 1):
    turned pixel (x, y) is unturned pixel (ow - y, x)."""
    from rasterio.control import GroundControlPoint
    from rasterio.transform import Affine
    if tf is not None:
        tf = tf * Affine(0, -1, ow, 1, 0, 0)
    if gc:
        g, gcrs = gc
        gc = ([GroundControlPoint(row=ow - p.col, col=p.row, x=p.x, y=p.y, z=p.z)
               for p in g], gcrs)
    return crs, tf, gc


ROW_TILE_KEYS = ["tile", "rel_path", "plot_stem", "row", "of_rows", "tile_k", "x0", "y0",
                 "w", "h", "out_w", "out_h", "turned_90", "resampled", "src_mm_px",
                 "out_mm_px", "tile_len_m", "row_width_m", "row_centre_px", "angle_deg",
                 "row_lo_px", "row_hi_px", "veg_pct", "src", "bits", "src_bits", "georef"]
REPORT_KEYS = ["plot_stem", "src", "ok", "reason", "rows_found", "rows_used", "pitch_mm",
               "angle_deg", "src_mm_px", "row_len_m", "rows_saved", "row_tiles", "veg_mode"]


def _write_csv(path, recs, keys):
    with open(path, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        wr.writeheader()
        wr.writerows(recs)


def _write_rows_qa(path, small, det, W, boxes, note=""):
    """QA picture: detected rows (row_detect colours) and every row tile written
    (cyan), labelled with its row number; note = why a plot got nothing."""
    import row_detect as rd
    from PIL import ImageDraw
    os.makedirs(os.path.dirname(path), exist_ok=True)
    im = rd.draw_rows(np.transpose(small, (1, 2, 0)), det)
    sc = im.size[0] / float(W)
    dr = ImageDraw.Draw(im)
    fnt = _qa_font(max(12, im.size[1] // 30))
    for (x, y, w, h, k) in boxes:
        dr.rectangle([x * sc, y * sc, (x + w) * sc, (y + h) * sc], outline=(0, 230, 255), width=2)
    for (x, y, w, h, k) in {b[4]: b for b in boxes}.values():
        dr.text((x * sc + 4, y * sc + 2), f"row {k}", fill=(255, 255, 0), font=fnt)
    if note:
        fs = _qa_font(max(14, im.size[1] // 22))
        dr.rectangle([0, 0, im.size[0], getattr(fs, "size", 18) + 12], fill=(0, 0, 0))
        dr.text((8, 4), note, fill=(255, 80, 80), font=fs)
    im.save(path, quality=88)


def make_row_tiles(src_dir, out_dir, rows="all", rows_per_plot=0, row_pitch_mm=0.0,
                   plot_width_mm=0.0, save_rows=True, row_tiles=True, tile_len_m=1.0,
                   row_width_m=0.0, overlap_m=0.0, out_mm_px=0.0, end_margin_m=0.0,
                   horizontal=True, max_tiles_per_row=0, tiles_pick="spread", seed=7,
                   fmt="tif", bits="same", compress="lzw", georef=True, link_source=True,
                   write_qa=True, log=print, progress=None, cancel=None):
    """Whole crop rows, and tiles of the same ground size along every row.

    rows            'all', 'middle N' or a list such as '1,2,3' (parse_rows).
    rows_per_plot   sown rows per plot when known (0 = the run of rows found).
    save_rows       out/<plot>/rows/<plot>_row<k>.<ext>: the whole row, full plot
                    length, crop pixels inside the row's band (centre +- half a
                    pitch, following the tilt), 0 outside.
    row_tiles       out/<plot>/row_tiles/row<k>/<plot>_row<k>_t<j>.<ext>: tiles
                    tile_len_m long along the row and row_width_m across, centred on
                    the row line. The row's length is the plot's (crop ends minus the
                    extraction margin, minus end_margin_m) - no vegetation threshold or
                    planted-length trim. Tiles are laid end to end (overlap_m) and
                    centred in that length; only whole tiles.
    row_width_m     0 = the median row spacing measured over all plots.
    out_mm_px       output scale; 0 = the median mm/px of the plots. Every row tile
                    is then exactly tile_len_m / out_mm_px x row_width_m / out_mm_px
                    pixels on every row of every plot.
    horizontal      turn row images and tiles 90 deg (lossless) when the rows run
                    up-down in the crop, so the row always runs left-right.
    max_tiles_per_row  0 = every tile that fits; else that many per row
                    ('spread' evenly or 'random').
    Writes rows_index.csv, row_tiles_index.csv, rows_report.csv and _rows_qa/.
    Returns dict(rows=..., tiles=..., plots=...)."""
    import rasterio
    import rasterio.windows
    import row_detect as rd

    ext = {"png": ".png", "jpeg": ".jpg", "jpg": ".jpg", "tif": ".tif"}.get(fmt.lower(), ".tif")
    bits = str(bits).strip().lower() or "same"
    compress = "none" if str(compress).lower().startswith(("none", "unc")) else "lzw"
    if ext == ".jpg" and bits == "16":
        raise ValueError("JPEG is 8-bit only: choose TIFF or PNG for 16-bit output.")
    if not (save_rows or row_tiles):
        raise ValueError("Choose whole rows, row tiles or both.")
    mode, sel = parse_rows(rows)
    links, _iiq_run, iiq_settings = (_source_links(src_dir) if link_source
                                     else ({}, {}, ""))
    os.makedirs(out_dir, exist_ok=True)
    files = sorted(f for f in os.listdir(src_dir) if f.lower().endswith((".tif", ".tiff")))
    if not files:
        log(f"[!] no GeoTIFF crops in {src_dir}")
        return dict(rows=0, tiles=0, plots=0)
    # extraction margin per crop: the plot's own ends lie that far inside the crop ends
    ext_margin = {}
    for mf in (os.path.join(os.path.dirname(os.path.abspath(src_dir)), "manifest.csv"),
               os.path.join(src_dir, "manifest.csv")):
        if os.path.exists(mf):
            try:
                with open(mf, newline="", encoding="utf-8") as fh:
                    for rr in csv.DictReader(fh):
                        try:
                            ext_margin[rr.get("raw_file", "")] = float(rr.get("margin_m") or 0)
                        except ValueError:
                            pass
            except Exception:
                pass
            break
    log(f"Rows of {len(files)} plot(s): rows '{rows or 'all'}'"
        + (f", {rows_per_plot} sown rows per plot" if rows_per_plot else "")
        + (", whole row images" if save_rows else "")
        + (f", {tile_len_m:g} m tiles along each row" if row_tiles else "") + ".")
    _out_desc(log, ext, bits, compress, georef and row_tiles, link_source and row_tiles)

    # ---- pass 1: rows of every plot, ground scale, row length ----------------------
    plan, report = {}, []
    for fi, fn in enumerate(files):
        if cancel is not None and cancel.is_set():
            log("Cancelled.")
            return dict(rows=0, tiles=0, plots=0)
        stem = os.path.splitext(fn)[0]
        info = dict(plot_stem=stem, src=fn, ok=False, reason="")
        report.append(info)
        try:
            with rasterio.open(os.path.join(src_dir, fn)) as s:
                W, H = s.width, s.height
                across_axis = 0 if H < W else 1
                src_mm = _plot_mm_px(s)
                if src_mm is None and plot_width_mm > 0:
                    src_mm = plot_width_mm / float(H if across_axis == 0 else W)
                dec = max(1, int(round(min(W, H) / 1200.0)))
                small = to8(s.read([1, 2, 3], out_shape=(3, max(1, H // dec), max(1, W // dec))))
        except Exception as e:
            info["reason"] = str(e)
            log(f"[!] {fn}: {e}")
            continue
        sx = W / float(small.shape[2]); sy = H / float(small.shape[1])
        f = sy if across_axis == 0 else sx                 # across scale small -> full
        fl = sx if across_axis == 0 else sy                # along scale
        det = rd.detect_rows(small, across_axis=across_axis,
                             mm_px=(src_mm * sx) if src_mm else None,
                             pitch_mm=row_pitch_mm or None, rows_expected=rows_per_plot)
        info.update(rows_found=len(det.get("plot_rows", [])),
                    angle_deg=round(det.get("angle_deg", 0.0), 2),
                    veg_mode=det.get("veg_mode", ""),
                    src_mm_px=round(src_mm, 4) if src_mm else "")

        def skipped_qa():                       # skipped plots get their picture too
            if write_qa:
                try:
                    _write_rows_qa(os.path.join(out_dir, "_rows_qa", stem + "_rows.jpg"),
                                   small, det, W, [], note=f"no rows written: {info['reason']}")
                except Exception as e:
                    log(f"[!] {fn}: QA picture: {e}")
        if not det.get("ok"):
            info["reason"] = det.get("reason", "rows not found")
            log(f"[!] {fn}: rows not found ({info['reason']}) - skipped")
            skipped_qa()
            continue
        if not src_mm:
            info["reason"] = "no ground scale (not georeferenced; give the plot width)"
            log(f"[!] {fn}: {info['reason']} - skipped")
            skipped_qa()
            continue
        pr = det["plot_rows"]
        n = len(pr)
        if mode == "all":
            ks = list(range(1, n + 1))
        elif mode == "middle":
            band = rd.middle_band(det, sel)
            ks = [pr.index(b) + 1 for b in band[2]] if band else []
        else:
            ks = [k for k in sel if k <= n]
            if len(ks) < len(sel):
                info["reason"] = f"rows {[k for k in sel if k > n]} not in this plot ({n} rows)"
        if not ks:
            info["reason"] = info["reason"] or f"{n} row(s) found, rows '{rows}' asked"
            log(f"[!] {fn}: {info['reason']} - skipped")
            skipped_qa()
            continue
        P = det["pitch_px"]
        ends = {}
        for k in ks:
            c = det["rows"][pr[k - 1]]
            us = _usable_length(small, across_axis, det, (c - P / 2.0, c + P / 2.0),
                                src_mm * fl, 1000.0 * ext_margin.get(fn, 0.0), 0.0)
            e = 1000.0 * max(0.0, end_margin_m) / src_mm
            ends[k] = (us["q0"] * fl + e, us["q1"] * fl - e)
        info.update(ok=True, rows_used=" ".join(map(str, ks)),
                    pitch_mm=round(P * f * src_mm, 1),
                    row_len_m=round(float(np.median([b - a for a, b in ends.values()]))
                                    * src_mm / 1000.0, 2))
        plan[fn] = dict(W=W, H=H, across_axis=across_axis, src_mm=src_mm, det=det, f=f,
                        fl=fl, ks=ks, ends=ends, info=info,
                        small=small if write_qa else None)
        if progress:
            progress(fi + 1, 2 * len(files))
    if not plan:
        log("[!] no plot with rows - nothing written")
        _write_csv(os.path.join(out_dir, "rows_report.csv"), report, REPORT_KEYS)
        return dict(rows=0, tiles=0, plots=0)

    gsd = np.array([p["src_mm"] for p in plan.values()])
    pitches = np.array([p["info"]["pitch_mm"] for p in plan.values()])
    out_mm = (out_mm_px if out_mm_px > 0 else          # plot nearest the median scale
              float(gsd[np.argmin(abs(gsd - np.median(gsd)))]))
    row_w = row_width_m if row_width_m > 0 else round(float(np.median(pitches)) / 1000.0, 3)
    out_len = max(8, int(round(tile_len_m * 1000.0 / out_mm)))
    out_wid = max(4, int(round(row_w * 1000.0 / out_mm)))
    log(f"{len(plan)}/{len(files)} plot(s) with rows; plot images {gsd.min():.3f}-"
        f"{gsd.max():.3f} mm/px (median {np.median(gsd):.3f}); row spacing "
        f"{pitches.min():.0f}-{pitches.max():.0f} mm (median {np.median(pitches):.0f}).")
    if row_tiles:
        log(f"Row tiles: {tile_len_m:g} m along x {row_w * 1000:.0f} mm across"
            + (" (median row spacing)" if row_width_m <= 0 else "")
            + f", output {out_mm:.3f} mm/px -> every tile {out_len} x {out_wid} px"
            + (" (rows run left-right)" if horizontal else "") + ".")

    # ---- pass 2: write -------------------------------------------------------------
    rnd = random.Random(seed)
    row_recs, tile_recs = [], []
    for fi, fn in enumerate(files):
        if cancel is not None and cancel.is_set():
            log("Cancelled.")
            break
        p = plan.get(fn)
        if p is None:
            continue
        stem = os.path.splitext(fn)[0]
        det, f, fl, ks = p["det"], p["f"], p["fl"], p["ks"]
        W, H, ax, src_mm = p["W"], p["H"], p["across_axis"], p["src_mm"]
        info = p["info"]
        pdir = os.path.join(out_dir, stem)
        qa_boxes = []
        try:
            with rasterio.open(os.path.join(src_dir, fn)) as s:
                if save_rows:
                    recs = _save_rows(s, W, H, ax, (det, f, fl), os.path.join(pdir, "rows"),
                                      stem, ext, bits, compress, select=set(ks),
                                      horizontal=horizontal)
                    for rc in recs:
                        rc["rel_path"] = f"{stem}/rows/{rc['file']}"
                    row_recs.extend(recs)
                    info["rows_saved"] = len(recs)
                if row_tiles:
                    crop_tags = s.tags() if link_source else {}
                    across_len = H if ax == 0 else W
                    along_len = W if ax == 0 else H
                    Li = int(round(tile_len_m * 1000.0 / src_mm))
                    Wi = int(round(row_w * 1000.0 / src_mm))
                    stride = max(1, Li - int(round(overlap_m * 1000.0 / src_mm)))
                    t = math.tan(math.radians(det["angle_deg"])) * (f / fl)
                    lc = det["along_centre"] * fl
                    turned = bool(horizontal and ax == 1)
                    n_written = 0
                    for k in ks:
                        c0 = det["rows"][det["plot_rows"][k - 1]] * f
                        L0, L1 = p["ends"][k]
                        span = L1 - L0
                        if span < Li or Wi > across_len:
                            continue
                        nt = int((span - Li) // stride) + 1
                        first = L0 + (span - (Li + (nt - 1) * stride)) / 2.0
                        starts = [(j, int(round(first + j * stride))) for j in range(nt)]
                        starts = _pick(starts, max_tiles_per_row, tiles_pick, rnd)
                        rdir = os.path.join(pdir, "row_tiles", f"row{k}")
                        for j, l0 in starts:
                            l0 = max(0, min(l0, along_len - Li))
                            ac = c0 + t * (l0 + Li / 2.0 - lc)
                            a0 = int(round(ac - Wi / 2.0))
                            a0c = max(0, min(a0, across_len - Wi))
                            if abs(a0c - a0) > 0.25 * Wi:
                                continue                # the row runs off the crop here
                            a0 = a0c
                            x, y, w, h = (l0, a0, Li, Wi) if ax == 0 else (a0, l0, Wi, Li)
                            a = s.read(window=rasterio.windows.Window(x, y, w, h))
                            if a.shape[0] < 3:
                                continue
                            rgb = a[:3]
                            v8 = to8(rgb)               # analysis only; rgb is written
                            valid = v8.max(axis=0) > 4
                            R, G, B = (v8[i].astype(np.float32) for i in range(3))
                            veg = (float(((2 * G - R - B)[valid] > 18).mean() * 100)
                                   if valid.any() else 0.0)
                            out = convert_bits(rgb, 8 if ext == ".jpg" else bits)
                            ow, oh = (out_len, out_wid) if ax == 0 else (out_wid, out_len)
                            resampled = (w, h) != (ow, oh)
                            if resampled:
                                out = _resize(out, ow, oh)
                            crs = tf = gc = None
                            if georef and ext == ".tif":
                                crs, tf, gc = _tile_georef(s, x, y, w, h, ow, oh)
                            if turned:
                                out = np.ascontiguousarray(np.rot90(out, 1, axes=(1, 2)))
                                if tf is not None or gc:
                                    crs, tf, gc = _turn_georef(crs, tf, gc, ow)
                            lk = (_link_fields(links.get(fn), crop_tags, x, y, w, h,
                                               iiq_settings) if link_source else {})
                            tags = None
                            if link_source and ext == ".tif":
                                tags = dict(PLOT_IMAGE=fn, PLOT_PIXEL_BOX=f"{x},{y},{w},{h}",
                                            ROW=f"{k} of {len(det['plot_rows'])}",
                                            TURNED_90=str(turned),
                                            SOURCE_IMAGE=lk.get("source_image", ""),
                                            PIXEL_SOURCE=lk.get("pixel_source", ""),
                                            SOURCE_PIXEL_BOX=(
                                                f"{lk['src_photo_x0']},{lk['src_photo_y0']},"
                                                f"{lk['src_photo_w']},{lk['src_photo_h']}"
                                                if lk.get("src_photo_x0", "") != "" else ""))
                            os.makedirs(rdir, exist_ok=True)
                            name = f"{stem}_row{k}_t{j + 1:02d}{ext}"
                            write_tile(os.path.join(rdir, name), out, ext, compress,
                                       crs, tf, gc, tags)
                            qa_boxes.append((x, y, w, h, k))
                            n_written += 1
                            tile_recs.append(dict(
                                tile=name, rel_path=f"{stem}/row_tiles/row{k}/{name}",
                                plot_stem=stem, row=k, of_rows=len(det["plot_rows"]),
                                tile_k=j + 1, x0=x, y0=y, w=w, h=h,
                                out_w=out.shape[2], out_h=out.shape[1], turned_90=turned,
                                resampled=resampled, src_mm_px=round(src_mm, 4),
                                out_mm_px=round(out_mm, 4), tile_len_m=tile_len_m,
                                row_width_m=row_w, row_centre_px=round(ac, 1),
                                angle_deg=round(det["angle_deg"], 2),
                                row_lo_px=int(round(L0)), row_hi_px=int(round(L1)),
                                veg_pct=round(veg, 2), src=fn,
                                bits=8 * out.dtype.itemsize, src_bits=8 * rgb.dtype.itemsize,
                                georef=("affine" if tf is not None else "gcps" if gc else ""),
                                **lk))
                    info["row_tiles"] = n_written
        except Exception as e:
            info["reason"] = str(e)
            log(f"[!] {fn}: {e}")
        if write_qa and p.get("small") is not None:
            try:
                note = ""
                if row_tiles and not qa_boxes:
                    note = ("no row tiles: the row is shorter than one tile or the "
                            "tile is wider than the plot")
                _write_rows_qa(os.path.join(out_dir, "_rows_qa", stem + "_rows.jpg"),
                               p["small"], det, W, qa_boxes, note=note)
            except Exception as e:
                log(f"[!] {fn}: QA picture: {e}")
            p["small"] = None
        if progress:
            progress(len(files) + fi + 1, 2 * len(files))
        if (fi + 1) % 25 == 0:
            log(f"   {fi + 1}/{len(files)} plots, {len(row_recs)} rows, {len(tile_recs)} tiles")

    if row_recs:
        _write_csv(os.path.join(out_dir, "rows_index.csv"), row_recs, list(row_recs[0].keys()))
        log(f"Saved {len(row_recs)} whole-row images -> <out>\\<plot>\\rows\\ "
            "(index rows_index.csv)")
    if tile_recs:
        keys = ROW_TILE_KEYS + (["source_image", "pixel_source", "src_photo_x0",
                                 "src_photo_y0", "src_photo_w", "src_photo_h",
                                 "iiq_settings", "link_note"] if link_source else [])
        _write_csv(os.path.join(out_dir, "row_tiles_index.csv"), tile_recs, keys)
        sz = {(r["out_w"], r["out_h"]) for r in tile_recs}
        log(f"Saved {len(tile_recs)} row tiles -> <out>\\<plot>\\row_tiles\\row<k>\\ "
            "(index row_tiles_index.csv)")
        log(f"Tile size check: pixel sizes {sorted(sz)} "
            f"({'all the same' if len(sz) == 1 else 'NOT uniform'}); "
            f"{sum(1 for r in tile_recs if not r['resampled'])} exact copies, "
            f"{sum(1 for r in tile_recs if r['resampled'])} resampled to {out_mm:.3f} mm/px")
    _write_csv(os.path.join(out_dir, "rows_report.csv"), report, REPORT_KEYS)
    plots = len({r["plot_stem"] for r in row_recs} | {r["plot_stem"] for r in tile_recs})
    log(f"Done: {plots} plot(s) -> {out_dir}; per-plot rows in rows_report.csv"
        + ("  (pictures in _rows_qa\\)" if write_qa else ""))
    return dict(rows=len(row_recs), tiles=len(tile_recs), plots=plots,
                tile_px=[out_len, out_wid] if row_tiles else None,
                out_mm_px=round(out_mm, 4), row_width_m=row_w)


# ---------------------------------------------------------------------------
# stitching back
# ---------------------------------------------------------------------------
BITS_CHOICES = ["same as plot image (exact copy)", "8", "16"]
COMPRESS_CHOICES = ["LZW (lossless, smaller)", "None (uncompressed)"]


def bits_key(label):
    s = str(label).strip()
    return s if s in ("8", "16") else "same"


def _out_desc(log, ext, bits, compress, georef, link_source):
    if ext == ".jpg":
        log("Tiles: JPEG quality 95 (lossy, 8-bit) - not an exact copy")
        return
    log("Tiles: " + ext[1:].upper()
        + (", LZW" if (ext == ".tif" and compress == "lzw") else
           ", uncompressed" if ext == ".tif" else "")
        + (", the plot image's own values (exact copy)" if bits == "same" else
           f", converted to {bits}-bit")
        + (", georeferenced" if (georef and ext == ".tif") else "")
        + (", linked to the source photo / IIQ" if link_source else ""))


def _recipe_buttons(frame, V, section, log):
    """'Load recipe...' button: fills the tab's fields from a PlotExtractor_recipe.json."""
    from tkinter import ttk, filedialog
    import recipe as rcp

    def load():
        f = filedialog.askopenfilename(
            title="Load recipe", filetypes=[("PlotExtractor recipe", "*.json"), ("All", "*.*")])
        if not f:
            return
        vals = rcp.gui_section(f, section)
        if not vals:
            log(f"[!] {os.path.basename(f)} has no '{section}' section")
            return
        n = rcp.apply_gui(V, vals)
        log(f"Recipe loaded: {n} setting(s) from {f}")
    b = ttk.Button(frame, text="Load recipe...", command=load)
    b.pack(side="left", padx=(12, 0))
    return b


def _gui_values(V):
    return {k: v.get() for k, v in V.items()}


def build_tiles_tab(parent, defaults=None):
    """Tab 6: two sub-tabs - 6a plot tiles (fixed ground size) and 6b rows and row
    tiles."""
    from tkinter import ttk
    nb = ttk.Notebook(parent)
    a = ttk.Frame(nb)
    b = ttk.Frame(nb)
    nb.add(a, text="6a  Plot tiles (fixed ground size)")
    nb.add(b, text="6b  Rows and row tiles")
    _build_plot_tiles_tab(a, defaults).pack(fill="both", expand=True)
    _build_row_tiles_tab(b, defaults).pack(fill="both", expand=True)
    return nb


def _build_plot_tiles_tab(parent, defaults=None):
    """Tab 6a: cut fixed-size ground tiles, and put labelled tiles back together.

    Self-contained - its own log and progress, so it does not depend on the
    extract tab being visible. Long jobs run on a worker thread.
    """
    import queue
    import threading
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox

    d = defaults or {}
    frm = ttk.Frame(parent, padding=10)
    frm.columnconfigure(1, weight=1)
    q = queue.Queue()
    cancel_ev = threading.Event()
    busy = {"on": False}

    V = {
        "src": tk.StringVar(value=d.get("src", "")),
        "out": tk.StringVar(value=d.get("out", "")),
        "tile_m": tk.StringVar(value=str(d.get("tile_m", 1.0))),
        "overlap_m": tk.StringVar(value=str(d.get("overlap_m", 0.0))),
        "rows_only": tk.BooleanVar(value=bool(d.get("rows_only", True))),
        "pitch": tk.StringVar(value=str(d.get("pitch", ""))),
        "plotw": tk.StringVar(value=str(d.get("plotw", ""))),
        "minveg": tk.StringVar(value=str(d.get("minveg", 0.5))),
        "mmpx": tk.StringVar(value=str(d.get("mmpx", ""))),
        "subdirs": tk.BooleanVar(value=bool(d.get("subdirs", False))),
        "keepempty": tk.StringVar(value=str(d.get("keepempty", 15))),
        "ntiles": tk.StringVar(value=str(d.get("ntiles", 0))),
        "fmt": tk.StringVar(value=d.get("fmt", "tif")),
        "seed": tk.StringVar(value=str(d.get("seed", 7))),
        "mid_on": tk.BooleanVar(value=bool(d.get("mid_on", False))),
        "mid_n": tk.StringVar(value=str(d.get("mid_n", 3))),
        "mid_rows": tk.StringVar(value=str(d.get("mid_rows", 0))),
        "mid_qa": tk.BooleanVar(value=bool(d.get("mid_qa", True))),
        "mid_cover": tk.BooleanVar(value=bool(d.get("mid_cover", True))),
        "mid_buffer": tk.StringVar(value=str(d.get("mid_buffer", 0.05))),
        "mid_along": tk.StringVar(value=str(d.get("mid_along", ""))),
        "mid_pick": tk.StringVar(value=d.get("mid_pick", "random")),
        "mid_end": tk.StringVar(value=str(d.get("mid_end", 0.5))),
        "bits": tk.StringVar(value=str(d.get("bits", BITS_CHOICES[0]))),
        "compress": tk.StringVar(value=d.get("compress", COMPRESS_CHOICES[0])),
        "georef": tk.BooleanVar(value=bool(d.get("georef", True))),
        "link": tk.BooleanVar(value=bool(d.get("link", True))),
        "stitch_in": tk.StringVar(value=""),
        "stitch_out": tk.StringVar(value=""),
        "grid": tk.BooleanVar(value=True),
    }

    r = 0

    def path_row(label, var, kind, hint=""):
        nonlocal r
        ttk.Label(frm, text=label).grid(row=r, column=0, sticky="w", pady=3)
        e = ttk.Entry(frm, textvariable=var)
        e.grid(row=r, column=1, sticky="ew", padx=6)

        def browse():
            v = (filedialog.askdirectory() if kind == "dir"
                 else filedialog.askopenfilename())
            if v:
                var.set(v)
        ttk.Button(frm, text="Browse...", command=browse).grid(row=r, column=2)
        if hint:
            r += 1
            ttk.Label(frm, text=hint, foreground="grey").grid(
                row=r, column=1, sticky="w", padx=6)
        r += 1

    ttk.Label(frm, text="Cut plot crops into tiles of a fixed ground size",
              font=("", 10, "bold")).grid(row=r, column=0, columnspan=3, sticky="w")
    r += 1
    path_row("Plot crops folder *", V["src"], "dir",
             "the RawCrops folder of an extraction (GeoTIFF/TIFF; raw or rectified)")
    path_row("Tiles output folder *", V["out"], "dir",
             "tiles are written flat here, plus tile_index.csv")

    o1 = ttk.Frame(frm)
    o1.grid(row=r, column=0, columnspan=3, sticky="w", pady=4)
    r += 1
    for text, var, w in (("Tile size (m):", V["tile_m"], 6),
                         ("Overlap (m):", V["overlap_m"], 6),
                         ("Seed:", V["seed"], 4)):
        ttk.Label(o1, text=text).pack(side="left", padx=(0, 4))
        ttk.Entry(o1, textvariable=var, width=w).pack(side="left", padx=(0, 12))
    ttk.Label(o1, text="Format:").pack(side="left", padx=(0, 4))
    ttk.Combobox(o1, textvariable=V["fmt"], width=6, state="readonly",
                 values=["tif", "png", "jpeg"]).pack(side="left")

    o1b = ttk.Frame(frm)
    o1b.grid(row=r, column=0, columnspan=3, sticky="w", pady=4)
    r += 1
    ttk.Label(o1b, text="Bits:").pack(side="left", padx=(0, 4))
    ttk.Combobox(o1b, textvariable=V["bits"], width=28, state="readonly",
                 values=BITS_CHOICES).pack(side="left", padx=(0, 12))
    ttk.Label(o1b, text="TIFF compression:").pack(side="left", padx=(0, 4))
    ttk.Combobox(o1b, textvariable=V["compress"], width=22, state="readonly",
                 values=COMPRESS_CHOICES).pack(side="left", padx=(0, 12))
    o1c = ttk.Frame(frm)
    o1c.grid(row=r, column=0, columnspan=3, sticky="w", pady=(0, 4))
    r += 1
    ttk.Checkbutton(o1c, text="Georeference tiles (GeoTIFF, from the plot image)",
                    variable=V["georef"]).pack(side="left", padx=(0, 12))
    ttk.Checkbutton(o1c, text="Link tiles to the source photo / IIQ (tags + index)",
                    variable=V["link"]).pack(side="left")
    ttk.Label(frm, text="tif + 'same as plot image' = every tile is an exact copy of its "
                        "window in the plot image (8 or 16 bit); LZW and uncompressed "
                        "hold identical pixels, LZW files are smaller. JPEG is lossy, 8-bit.",
              foreground="grey").grid(row=r, column=0, columnspan=3, sticky="w")
    r += 1

    o1n = ttk.Frame(frm)
    o1n.grid(row=r, column=0, columnspan=3, sticky="w", pady=4)
    r += 1
    ttk.Label(o1n, text="Tiles per plot:", font=("", 9, "bold")).pack(
        side="left", padx=(0, 4))
    ttk.Spinbox(o1n, from_=0, to=500, textvariable=V["ntiles"], width=5).pack(
        side="left", padx=(0, 4))
    ttk.Label(o1n, text="(0 = every tile that fits)").pack(side="left", padx=(0, 12))
    ttk.Label(o1n, text="choose them:").pack(side="left", padx=(0, 4))
    ttk.Combobox(o1n, textvariable=V["mid_pick"], width=8, state="readonly",
                 values=["random", "spread"]).pack(side="left", padx=(0, 6))
    ttk.Label(o1n, text="(spread = evenly spaced along the plot)",
              foreground="grey").pack(side="left")

    o2 = ttk.Frame(frm)
    o2.grid(row=r, column=0, columnspan=3, sticky="w", pady=4)
    r += 1
    ttk.Checkbutton(o2, text="Tile only the seed-row band",
                    variable=V["rows_only"]).pack(side="left", padx=(0, 12))
    ttk.Label(o2, text="row spacing (mm):").pack(side="left", padx=(0, 4))
    ttk.Entry(o2, textvariable=V["pitch"], width=6).pack(side="left", padx=(0, 12))
    ttk.Label(o2, text="plot width (mm):").pack(side="left", padx=(0, 4))
    ttk.Entry(o2, textvariable=V["plotw"], width=7).pack(side="left")
    ttk.Label(frm, text="the band is measured per plot; the two numbers only "
                        "set the search scale (blank = guess)",
              foreground="grey").grid(row=r, column=0, columnspan=3, sticky="w")
    r += 1

    o2m = ttk.Frame(frm)
    o2m.grid(row=r, column=0, columnspan=3, sticky="w", pady=4)
    r += 1
    ttk.Checkbutton(o2m, text="Tile only the middle", variable=V["mid_on"]).pack(
        side="left", padx=(0, 4))
    ttk.Spinbox(o2m, from_=1, to=12, textvariable=V["mid_n"], width=3).pack(
        side="left", padx=(0, 4))
    ttk.Label(o2m, text="rows of each plot;  sown rows per plot (0 = detect):").pack(
        side="left", padx=(0, 4))
    ttk.Spinbox(o2m, from_=0, to=24, textvariable=V["mid_rows"], width=3).pack(
        side="left", padx=(0, 12))
    o2c = ttk.Frame(frm)
    o2c.grid(row=r, column=0, columnspan=3, sticky="w", pady=(0, 4))
    r += 1
    ttk.Label(o2c, text="      ").pack(side="left")
    ttk.Checkbutton(o2c, text="Size the tiles to cover these rows:  N x row spacing + buffer",
                    variable=V["mid_cover"]).pack(side="left", padx=(0, 4))
    ttk.Entry(o2c, textvariable=V["mid_buffer"], width=5).pack(side="left", padx=(0, 4))
    ttk.Label(o2c, text="m each side;  length along the rows").pack(side="left", padx=(0, 4))
    ttk.Entry(o2c, textvariable=V["mid_along"], width=5).pack(side="left", padx=(0, 4))
    ttk.Label(o2c, text="m (blank = square)", foreground="grey").pack(side="left")
    o2n = ttk.Frame(frm)
    o2n.grid(row=r, column=0, columnspan=3, sticky="w", pady=(0, 4))
    r += 1
    ttk.Label(o2n, text="      keep tiles").pack(side="left", padx=(0, 4))
    ttk.Entry(o2n, textvariable=V["mid_end"], width=5).pack(side="left", padx=(0, 4))
    ttk.Label(o2n, text="m from the plot ends").pack(side="left", padx=(0, 16))
    ttk.Label(o2n, text="(whole rows and row tiles: tab 6b)",
              foreground="grey").pack(side="left")
    ttk.Label(frm, text="rows are detected in every plot (PhenoFuse row detection, "
                        "tilt-aware, works on raw and rectified crops); the band "
                        "is centred on the N central rows. Along the plot, tiles "
                        "stay inside the planted length, the set distance away "
                        "from the bare plot ends, and are centred there; 'tiles "
                        "per plot' then takes that many (random, or 'spread' = "
                        "evenly spaced). 'Cover these rows' measures the row spacing "
                        "of every plot first and makes ONE tile size for the run "
                        "(e.g. 3 x 215 mm + 2 x 0.05 m = 0.745 m), so every tile holds "
                        "the N rows and all tiles have the same pixel size (Tile size "
                        "is then not used). Give the sown rows per plot when known - "
                        "it decides which rows are central when the grid clips an "
                        "outer row. Set a margin at extraction so outer rows are "
                        "whole.",
              foreground="grey", wraplength=900, justify="left").grid(
        row=r, column=0, columnspan=3, sticky="w")
    r += 1

    o2b = ttk.Frame(frm)
    o2b.grid(row=r, column=0, columnspan=3, sticky="w", pady=4)
    r += 1
    ttk.Label(o2b, text="Output scale (mm/px):").pack(side="left", padx=(0, 4))
    ttk.Entry(o2b, textvariable=V["mmpx"], width=7).pack(side="left")
    ttk.Label(frm, text="every tile covers the same ground (tile size x tile size m) AND "
                        "has the same pixel size on every plot: tile m / mm/px. Blank = "
                        "the median mm/px of the plot images; plots at that scale are "
                        "exact copies, the others are resampled (Lanczos). Only whole "
                        "tiles are written.",
              foreground="grey").grid(row=r, column=0, columnspan=3, sticky="w")
    r += 1

    o3 = ttk.Frame(frm)
    o3.grid(row=r, column=0, columnspan=3, sticky="w", pady=4)
    r += 1
    ttk.Label(o3, text="Drop tiles below vegetation %:").pack(side="left", padx=(0, 4))
    ttk.Entry(o3, textvariable=V["minveg"], width=5).pack(side="left", padx=(0, 12))
    ttk.Label(o3, text="but keep this % of them as hard negatives:").pack(side="left", padx=(0, 4))
    ttk.Entry(o3, textvariable=V["keepempty"], width=5).pack(side="left")

    o4 = ttk.Frame(frm)
    o4.grid(row=r, column=0, columnspan=3, sticky="w", pady=4)
    r += 1
    ttk.Checkbutton(o4, text="One folder per plot (out/<plot>/...) instead of "
                            "one flat folder",
                    variable=V["subdirs"]).pack(side="left", padx=(0, 16))
    ttk.Checkbutton(o4, text="QA pictures, one per plot (_rows_qa\\ with middle rows, "
                            "else _tiles_qa\\)",
                    variable=V["mid_qa"]).pack(side="left")

    btn1 = ttk.Frame(frm)
    btn1.grid(row=r, column=0, columnspan=3, sticky="w", pady=6)
    r += 1
    make_btn = ttk.Button(btn1, text="Make tiles")
    make_btn.pack(side="left")
    cancel_btn = ttk.Button(btn1, text="Cancel", state="disabled")
    cancel_btn.pack(side="left", padx=8)
    prog = ttk.Progressbar(btn1, mode="determinate", length=260)
    prog.pack(side="left", padx=8)
    _recipe_buttons(btn1, V, "plot_tiles", lambda m: put(m))
    ttk.Label(btn1, text="each run writes its settings to PlotExtractor_recipe.json "
                         "in the output folder", foreground="grey").pack(side="left", padx=8)

    ttk.Separator(frm, orient="horizontal").grid(
        row=r, column=0, columnspan=3, sticky="ew", pady=8)
    r += 1
    ttk.Label(frm, text="Rebuild one image per plot",
              font=("", 10, "bold")).grid(row=r, column=0, columnspan=3, sticky="w")
    r += 1
    path_row("Tiles or labelled tiles *", V["stitch_in"], "dir",
             "the tiles folder, or a folder of SAM3 overlays with the same names")
    path_row("Stitched output folder *", V["stitch_out"], "dir")
    b2 = ttk.Frame(frm)
    b2.grid(row=r, column=0, columnspan=3, sticky="w", pady=4)
    r += 1
    stitch_btn = ttk.Button(b2, text="Combine tiles into plot images")
    stitch_btn.pack(side="left")
    ttk.Checkbutton(b2, text="draw tile boundaries", variable=V["grid"]).pack(
        side="left", padx=10)

    box = tk.Text(frm, height=12, wrap="none", state="disabled",
                  font=("Consolas", 9))
    box.grid(row=r, column=0, columnspan=3, sticky="nsew")
    frm.rowconfigure(r, weight=1)
    sb = ttk.Scrollbar(frm, command=box.yview)
    sb.grid(row=r, column=3, sticky="ns")
    box.configure(yscrollcommand=sb.set)

    def put(m):
        q.put(("log", str(m)))

    def _f(var, default=0.0):
        try:
            return float(str(var.get()).strip() or default)
        except (TypeError, ValueError):
            return default

    def start(fn):
        if busy["on"]:
            return
        busy["on"] = True
        cancel_ev.clear()
        make_btn.configure(state="disabled")
        stitch_btn.configure(state="disabled")
        cancel_btn.configure(state="normal")
        prog["value"] = 0

        def wrap():
            try:
                fn()
            except Exception as e:
                put("[!] " + str(e))
            q.put(("done", None))
        threading.Thread(target=wrap, daemon=True).start()

    def on_make():
        src, out = V["src"].get().strip(), V["out"].get().strip()
        if not (src and os.path.isdir(src)):
            messagebox.showwarning("Make tiles", "Choose the plot crops folder.")
            return
        if not out:
            messagebox.showwarning("Make tiles", "Choose an output folder.")
            return
        if os.path.abspath(src) == os.path.abspath(out):
            messagebox.showwarning("Make tiles",
                                   "The output folder must not be the crops folder.")
            return
        tm = _f(V["tile_m"], 1.0); om = _f(V["overlap_m"], 0.0)
        cover = bool(V["mid_on"].get() and V["mid_cover"].get())
        if tm <= 0 and not cover:
            messagebox.showwarning("Make tiles", "Tile size must be more than 0 m.")
            return
        if om >= tm and not cover:
            messagebox.showwarning("Make tiles",
                                   "Overlap must be smaller than the tile size.")
            return
        kw = dict(
            tile_m=tm, overlap_m=om, cover_rows=cover,
            row_buffer_m=max(0.0, _f(V["mid_buffer"], 0.05)),
            tile_along_m=max(0.0, _f(V["mid_along"], 0.0)),
            rows_only=bool(V["rows_only"].get()),
            row_pitch_mm=_f(V["pitch"], 0.0), plot_width_mm=_f(V["plotw"], 0.0),
            fmt=V["fmt"].get(), min_veg_pct=_f(V["minveg"], 0.0),
            keep_empty_frac=max(0.0, min(1.0, _f(V["keepempty"], 15.0) / 100.0)),
            max_tiles_per_plot=max(0, int(_f(V["ntiles"], 0))),
            seed=int(_f(V["seed"], 7)), target_mm_px=_f(V["mmpx"], 0.0),
            per_plot_dirs=bool(V["subdirs"].get()),
            middle_rows=(int(_f(V["mid_n"], 3)) if V["mid_on"].get() else 0),
            rows_per_plot=int(_f(V["mid_rows"], 0)),
            band_tile=False, write_qa=bool(V["mid_qa"].get()),
            end_margin_m=_f(V["mid_end"], 0.5), tiles_pick=V["mid_pick"].get() or "random",
            save_rows=False,
            bits=bits_key(V["bits"].get()),
            compress="none" if V["compress"].get().lower().startswith("none") else "lzw",
            georef=bool(V["georef"].get()), link_source=bool(V["link"].get()))
        gui = _gui_values(V)

        def job():
            import recipe as rcp
            prm = dict(src_dir=src, out_dir=out, **kw)
            rcp.write_recipe(out, "plot_tiles", gui=gui, params=prm, src_dir=src)
            n, plots, idx = make_tiles(src, out, log=put,
                                       progress=lambda a, b: q.put(("prog", (a, b))),
                                       cancel=cancel_ev, **kw)
            rcp.write_recipe(out, "plot_tiles", gui=gui, params=prm,
                             results=dict(tiles=n, plots=plots, index=idx,
                                          cancelled=cancel_ev.is_set(),
                                          **getattr(make_tiles, "last_run", {})),
                             src_dir=src, log=put)
        start(job)

    def on_stitch():
        src = V["stitch_in"].get().strip() or V["out"].get().strip()
        out = V["stitch_out"].get().strip()
        if not (src and os.path.isdir(src)):
            messagebox.showwarning("Combine tiles", "Choose the tiles folder.")
            return
        if not out:
            messagebox.showwarning("Combine tiles", "Choose an output folder.")
            return
        idx = os.path.join(src, "tile_index.csv")
        if not os.path.exists(idx):
            idx = os.path.join(V["out"].get().strip() or src, "tile_index.csv")
        start(lambda: stitch_tiles(src, out, index_csv=idx, log=put,
                                   progress=lambda a, b: q.put(("prog", (a, b))),
                                   draw_grid=bool(V["grid"].get()),
                                   cancel=cancel_ev))

    make_btn.configure(command=on_make)
    stitch_btn.configure(command=on_stitch)
    cancel_btn.configure(command=lambda: (cancel_ev.set(), put("cancelling...")))

    def poll():
        try:
            while True:
                kind, payload = q.get_nowait()
                if kind == "log":
                    box.configure(state="normal")
                    box.insert("end", payload + "\n")
                    box.see("end")
                    box.configure(state="disabled")
                elif kind == "prog":
                    a, b = payload
                    prog["maximum"] = b
                    prog["value"] = a
                elif kind == "done":
                    busy["on"] = False
                    make_btn.configure(state="normal")
                    stitch_btn.configure(state="normal")
                    cancel_btn.configure(state="disabled")
        except queue.Empty:
            pass
        frm.after(150, poll)

    poll()
    frm._pe_vars = V
    return frm


def _build_row_tiles_tab(parent, defaults=None):
    """Tab 6b: whole crop rows (all, middle N or chosen rows) and tiles of one ground
    size along each row, written per plot:
        <out>/<plot>/rows/<plot>_row<k>.<ext>
        <out>/<plot>/row_tiles/row<k>/<plot>_row<k>_t<j>.<ext>"""
    import queue
    import threading
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox

    d = defaults or {}
    frm = ttk.Frame(parent, padding=10)
    frm.columnconfigure(1, weight=1)
    q = queue.Queue()
    cancel_ev = threading.Event()
    busy = {"on": False}
    out0 = d.get("out", "")
    V = {
        "src": tk.StringVar(value=d.get("src", "")),
        "out": tk.StringVar(value=(os.path.join(os.path.dirname(out0), "Rows")
                                   if out0 else "")),
        "rows": tk.StringVar(value="all"),
        "rows_per_plot": tk.StringVar(value="0"),
        "pitch": tk.StringVar(value=str(d.get("pitch", ""))),
        "plotw": tk.StringVar(value=str(d.get("plotw", ""))),
        "save_rows": tk.BooleanVar(value=True),
        "horizontal": tk.BooleanVar(value=True),
        "row_tiles": tk.BooleanVar(value=True),
        "tile_len": tk.StringVar(value="1.0"),
        "row_w": tk.StringVar(value=""),
        "overlap_m": tk.StringVar(value="0"),
        "end_m": tk.StringVar(value="0"),
        "ntiles": tk.StringVar(value="0"),
        "pick": tk.StringVar(value="spread"),
        "mmpx": tk.StringVar(value=""),
        "fmt": tk.StringVar(value="tif"),
        "bits": tk.StringVar(value=BITS_CHOICES[0]),
        "compress": tk.StringVar(value=COMPRESS_CHOICES[0]),
        "georef": tk.BooleanVar(value=True),
        "link": tk.BooleanVar(value=True),
        "qa": tk.BooleanVar(value=True),
        "seed": tk.StringVar(value="7"),
    }
    r = 0

    def path_row(label, var, hint=""):
        nonlocal r
        ttk.Label(frm, text=label).grid(row=r, column=0, sticky="w", pady=3)
        ttk.Entry(frm, textvariable=var).grid(row=r, column=1, sticky="ew", padx=6)

        def browse():
            v = filedialog.askdirectory()
            if v:
                var.set(v)
        ttk.Button(frm, text="Browse...", command=browse).grid(row=r, column=2)
        if hint:
            r += 1
            ttk.Label(frm, text=hint, foreground="grey").grid(
                row=r, column=1, sticky="w", padx=6)
        r += 1

    def line(pady=4):
        nonlocal r
        f = ttk.Frame(frm)
        f.grid(row=r, column=0, columnspan=3, sticky="w", pady=pady)
        r += 1
        return f

    def note(text):
        nonlocal r
        ttk.Label(frm, text=text, foreground="grey", wraplength=900, justify="left").grid(
            row=r, column=0, columnspan=3, sticky="w")
        r += 1

    def field(f, text, var, w=6, after=""):
        ttk.Label(f, text=text).pack(side="left", padx=(0, 4))
        ttk.Entry(f, textvariable=var, width=w).pack(side="left", padx=(0, 4))
        if after:
            ttk.Label(f, text=after).pack(side="left", padx=(0, 12))
        else:
            ttk.Frame(f, width=8).pack(side="left")

    ttk.Label(frm, text="Extract whole rows, and tiles of a fixed ground size along each row",
              font=("", 10, "bold")).grid(row=r, column=0, columnspan=3, sticky="w")
    r += 1
    path_row("Plot crops folder *", V["src"],
             "the RawCrops folder of an extraction (GeoTIFF/TIFF; raw or rectified)")
    path_row("Output folder *", V["out"],
             "one folder per plot: <plot>\\rows\\ and <plot>\\row_tiles\\row<k>\\")

    f = line()
    ttk.Label(f, text="Rows:", font=("", 9, "bold")).pack(side="left", padx=(0, 4))
    ttk.Combobox(f, textvariable=V["rows"], width=12,
                 values=["all", "middle 3", "middle 1", "1", "1,2", "1,2,3", "2,3,4"]).pack(
        side="left", padx=(0, 6))
    ttk.Label(f, text="('all', 'middle N', or row numbers such as 1,2,3 or 2-4)",
              foreground="grey").pack(side="left", padx=(0, 12))
    ttk.Label(f, text="sown rows per plot (0 = detect):").pack(side="left", padx=(0, 4))
    ttk.Spinbox(f, from_=0, to=24, textvariable=V["rows_per_plot"], width=3).pack(side="left")
    f = line()
    field(f, "row spacing (mm):", V["pitch"], 6)
    field(f, "plot width (mm):", V["plotw"], 7)
    note("Rows are found in every plot (PhenoFuse row detection, tilt-aware, raw and "
         "rectified crops) and numbered 1..n across the plot from the image top (or left "
         "when the rows run up-down). Give the sown rows per plot when known: it decides "
         "which rows count when the grid clips an outer row. Spacing and width only set "
         "the search scale (blank = guess).")

    f = line()
    ttk.Checkbutton(f, text="Save each chosen row as a whole image (<plot>\\rows\\)",
                    variable=V["save_rows"]).pack(side="left", padx=(0, 16))
    ttk.Checkbutton(f, text="Turn so the row runs left-right (lossless 90 deg)",
                    variable=V["horizontal"]).pack(side="left")
    note("A whole row is the full plot length; the row's own band (centre +- half the "
         "spacing, following the tilt) holds the crop's pixels unchanged, outside it is 0.")

    f = line()
    ttk.Checkbutton(f, text="Make row tiles:", variable=V["row_tiles"]).pack(
        side="left", padx=(0, 8))
    field(f, "tile length along the row (m):", V["tile_len"], 5)
    field(f, "tile width across (m):", V["row_w"], 6, "(blank = median row spacing)")
    f = line()
    ttk.Label(f, text="      ").pack(side="left")
    field(f, "overlap (m):", V["overlap_m"], 5)
    field(f, "keep tiles", V["end_m"], 5, "m from the plot ends")
    ttk.Label(f, text="tiles per row:").pack(side="left", padx=(0, 4))
    ttk.Spinbox(f, from_=0, to=100, textvariable=V["ntiles"], width=4).pack(
        side="left", padx=(0, 4))
    ttk.Label(f, text="(0 = all)").pack(side="left", padx=(0, 8))
    ttk.Combobox(f, textvariable=V["pick"], width=8, state="readonly",
                 values=["spread", "random"]).pack(side="left")
    f = line()
    ttk.Label(f, text="      ").pack(side="left")
    field(f, "Output scale (mm/px):", V["mmpx"], 7, "(blank = median of the plots)")
    note("Each row tile is centred on its row line and has the same ground size and the "
         "same pixel size on every row of every plot (length / mm/px x width / mm/px). "
         "The row length is the plot's: crop ends minus the extraction margin (and the "
         "'keep' distance) - no vegetation threshold or planted-length trim. Tiles are "
         "laid end to end and centred along the row; only whole tiles are written.")

    f = line()
    ttk.Label(f, text="Format:").pack(side="left", padx=(0, 4))
    ttk.Combobox(f, textvariable=V["fmt"], width=6, state="readonly",
                 values=["tif", "png", "jpeg"]).pack(side="left", padx=(0, 12))
    ttk.Label(f, text="Bits:").pack(side="left", padx=(0, 4))
    ttk.Combobox(f, textvariable=V["bits"], width=28, state="readonly",
                 values=BITS_CHOICES).pack(side="left", padx=(0, 12))
    ttk.Label(f, text="TIFF compression:").pack(side="left", padx=(0, 4))
    ttk.Combobox(f, textvariable=V["compress"], width=22, state="readonly",
                 values=COMPRESS_CHOICES).pack(side="left")
    f = line()
    ttk.Checkbutton(f, text="Georeference row tiles", variable=V["georef"]).pack(
        side="left", padx=(0, 12))
    ttk.Checkbutton(f, text="Link to the source photo / IIQ", variable=V["link"]).pack(
        side="left", padx=(0, 12))
    ttk.Checkbutton(f, text="QA pictures (_rows_qa\\)", variable=V["qa"]).pack(
        side="left", padx=(0, 12))
    field(f, "Seed:", V["seed"], 4)

    f = line(6)
    run_btn = ttk.Button(f, text="Extract rows")
    run_btn.pack(side="left")
    cancel_btn = ttk.Button(f, text="Cancel", state="disabled")
    cancel_btn.pack(side="left", padx=8)
    prog = ttk.Progressbar(f, mode="determinate", length=260)
    prog.pack(side="left", padx=8)

    box = tk.Text(frm, height=12, wrap="none", state="disabled", font=("Consolas", 9))

    def put(m):
        q.put(("log", str(m)))
    _recipe_buttons(f, V, "rows_and_row_tiles", put)
    ttk.Label(f, text="each run writes its settings to PlotExtractor_recipe.json",
              foreground="grey").pack(side="left", padx=8)

    box.grid(row=r, column=0, columnspan=3, sticky="nsew")
    frm.rowconfigure(r, weight=1)
    sb = ttk.Scrollbar(frm, command=box.yview)
    sb.grid(row=r, column=3, sticky="ns")
    box.configure(yscrollcommand=sb.set)

    def _f(var, default=0.0):
        try:
            return float(str(var.get()).strip() or default)
        except (TypeError, ValueError):
            return default

    def on_run():
        if busy["on"]:
            return
        src, out = V["src"].get().strip(), V["out"].get().strip()
        if not (src and os.path.isdir(src)):
            messagebox.showwarning("Extract rows", "Choose the plot crops folder.")
            return
        if not out:
            messagebox.showwarning("Extract rows", "Choose an output folder.")
            return
        if os.path.abspath(src) == os.path.abspath(out):
            messagebox.showwarning("Extract rows",
                                   "The output folder must not be the crops folder.")
            return
        try:
            parse_rows(V["rows"].get())
        except ValueError as e:
            messagebox.showwarning("Extract rows", str(e))
            return
        if not (V["save_rows"].get() or V["row_tiles"].get()):
            messagebox.showwarning("Extract rows", "Choose whole rows, row tiles or both.")
            return
        tl, om = _f(V["tile_len"], 1.0), _f(V["overlap_m"], 0.0)
        if V["row_tiles"].get() and (tl <= 0 or om >= tl):
            messagebox.showwarning("Extract rows",
                                   "Tile length must be > 0 and longer than the overlap.")
            return
        kw = dict(rows=V["rows"].get().strip() or "all",
                  rows_per_plot=int(_f(V["rows_per_plot"], 0)),
                  row_pitch_mm=_f(V["pitch"], 0.0), plot_width_mm=_f(V["plotw"], 0.0),
                  save_rows=bool(V["save_rows"].get()), row_tiles=bool(V["row_tiles"].get()),
                  tile_len_m=tl, row_width_m=_f(V["row_w"], 0.0), overlap_m=om,
                  out_mm_px=_f(V["mmpx"], 0.0), end_margin_m=_f(V["end_m"], 0.0),
                  horizontal=bool(V["horizontal"].get()),
                  max_tiles_per_row=max(0, int(_f(V["ntiles"], 0))),
                  tiles_pick=V["pick"].get() or "spread", seed=int(_f(V["seed"], 7)),
                  fmt=V["fmt"].get(), bits=bits_key(V["bits"].get()),
                  compress="none" if V["compress"].get().lower().startswith("none") else "lzw",
                  georef=bool(V["georef"].get()), link_source=bool(V["link"].get()),
                  write_qa=bool(V["qa"].get()))
        gui = _gui_values(V)
        busy["on"] = True
        cancel_ev.clear()
        run_btn.configure(state="disabled")
        cancel_btn.configure(state="normal")
        prog["value"] = 0

        def job():
            import recipe as rcp
            prm = dict(src_dir=src, out_dir=out, **kw)
            try:
                rcp.write_recipe(out, "rows_and_row_tiles", gui=gui, params=prm, src_dir=src)
                res = make_row_tiles(src, out, log=put,
                                     progress=lambda a, b: q.put(("prog", (a, b))),
                                     cancel=cancel_ev, **kw)
                res["cancelled"] = cancel_ev.is_set()
                rcp.write_recipe(out, "rows_and_row_tiles", gui=gui, params=prm,
                                 results=res, src_dir=src, log=put)
            except Exception as e:
                put("[!] " + str(e))
            q.put(("done", None))
        threading.Thread(target=job, daemon=True).start()

    run_btn.configure(command=on_run)
    cancel_btn.configure(command=lambda: (cancel_ev.set(), put("cancelling...")))

    def poll():
        try:
            while True:
                kind, payload = q.get_nowait()
                if kind == "log":
                    box.configure(state="normal")
                    box.insert("end", payload + "\n")
                    box.see("end")
                    box.configure(state="disabled")
                elif kind == "prog":
                    a, b = payload
                    prog["maximum"] = b
                    prog["value"] = a
                elif kind == "done":
                    busy["on"] = False
                    run_btn.configure(state="normal")
                    cancel_btn.configure(state="disabled")
        except queue.Empty:
            pass
        frm.after(150, poll)

    poll()
    frm._pe_vars = V
    frm._pe_run = on_run
    return frm


def stitch_tiles(tile_dir, out_dir, index_csv="", log=print, progress=None,
                 draw_grid=False, cancel=None):
    """Rebuild one image per plot from a folder of tiles.

    Works on the tiles themselves or on anything produced from them that kept the
    file names - SAM3 overlays, for instance - so you can look at a whole plot's
    labels instead of 60 separate tiles.

    Overlapping tiles simply overwrite; this is for looking at, not measuring.
    """
    from PIL import Image, ImageDraw

    index_csv = index_csv or os.path.join(tile_dir, "tile_index.csv")
    if not os.path.exists(index_csv):
        log(f"[!] no tile_index.csv - looked in {index_csv}")
        return 0
    rows = list(csv.DictReader(open(index_csv)))
    os.makedirs(out_dir, exist_ok=True)
    by_plot = {}
    for r in rows:
        by_plot.setdefault(r["plot_stem"], []).append(r)
    log(f"Rebuilding {len(by_plot)} plot(s) from {len(rows)} tile records.")
    made = 0
    for pi, (stem, rs) in enumerate(sorted(by_plot.items())):
        if cancel is not None and cancel.is_set():
            log("Stitch cancelled.")
            break
        try:
            ax = int(rs[0]["across_axis"])
            lo, hi = int(rs[0]["band_lo"]), int(rs[0]["band_hi"])
            W, H = int(rs[0]["src_w"]), int(rs[0]["src_h"])
            y0b, y1b = (lo, hi) if ax == 0 else (0, H)
            x0b, x1b = (0, W) if ax == 0 else (lo, hi)
            canvas = Image.new("RGB", (x1b - x0b, y1b - y0b), (26, 26, 26))
            drew = 0
            for r in rs:
                # accept any extension: the labelling step may have changed it
                base = os.path.splitext(r["tile"])[0]
                # tiles may be flat, or one folder per plot; rel_path records
                # which, and a labelling step may have changed the extension
                homes = [tile_dir, os.path.join(tile_dir, stem)]
                rel = r.get("rel_path") or r["tile"]
                cand = [os.path.join(tile_dir, rel.replace("/", os.sep))]
                cand += [os.path.join(h, base + e) for h in homes
                         for e in (".png", ".jpg", ".jpeg", ".tif")]
                src = next((c for c in cand if os.path.exists(c)), None)
                if not src:
                    continue
                tim = Image.fromarray(to8(read_rgb(src)))
                if str(r.get("turned_90", "")).lower() == "true":
                    tim = tim.transpose(Image.ROTATE_270)      # undo np.rot90(.., 1)
                # ground-referenced tiles were resampled off the source grid;
                # put them back at source scale so the plot reassembles square
                wh = (int(r["w"]), int(r["h"]))
                if tim.size != wh:
                    tim = tim.resize(wh, Image.LANCZOS)
                canvas.paste(tim, (int(r["x0"]) - x0b, int(r["y0"]) - y0b))
                drew += 1
            if not drew:
                continue
            if draw_grid:
                d = ImageDraw.Draw(canvas)
                for r in rs:
                    x, y = int(r["x0"]) - x0b, int(r["y0"]) - y0b
                    d.rectangle([x, y, x + int(r["w"]) - 1, y + int(r["h"]) - 1],
                                outline=(255, 210, 0), width=3)
            canvas.save(os.path.join(out_dir, f"{stem}_stitched.jpg"), quality=90)
            made += 1
        except Exception as e:
            log(f"[!] {stem}: {e}")
        if progress:
            progress(pi + 1, len(by_plot))
    log(f"Rebuilt {made} plot image(s) -> {out_dir}")
    return made
