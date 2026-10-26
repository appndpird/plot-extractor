"""row_detect.py - find the crop rows in one plot image (raw or rectified crop).

Ported from PhenoFuse 4.0 (phenoapp/core/grid_rows.py: bin_rgb_windows, _row_signal,
row_pitch_from_profile, rows_from_profile) so both tools find the same rows. PhenoFuse
works on a 2 cm raster in each plot's ground frame; here the same steps run in the
pixels of a single plot crop:

  vegetation    excess green on chromatic coordinates, ExG = (2G - R - B) / (R + G + B);
                vegetation = ExG above one Otsu threshold for the plot. The share of
                vegetation pixels per across-position is the row profile, and the
                continuous ExG is added as a second ridge signal (seedlings are a few
                pixels wide; the share alone is noisy at emergence).
  pitch         peak of the summed power spectra of short segments along the plot,
                inside the seeder range (17-32 cm when the ground sampling is known),
                or the given spacing; refined on the de-tilted whole-plot profile.
  tilt          raw crops are camera frames, so the rows are not always parallel to an
                image axis. The row phase is measured in each segment and its drift
                along the plot is fitted as a straight line (robust weighted LS);
                the slope is the tilt (kept within +-max_tilt_deg).
  comb          n teeth at that pitch fill the planted span; only the phase is fitted
                (ridge at the teeth minus signal at the mid-gaps). Every tooth gets a
                contrast; contrast > 0 = a real row (non-positive = missing / bare).
  plot rows     the run of present teeth (or the best block of rows_expected teeth);
                middle_band() returns the band holding the N central rows.

Pure numpy (the PyInstaller build excludes scipy).
"""
from __future__ import annotations

import math

import numpy as np

ROW_PITCH_RANGE_MM = (170.0, 320.0)   # cereal / canola plot seeders (PhenoFuse ROW_PITCH_RANGE)


# ---------------------------------------------------------------------------
# helpers (numpy stand-ins for scipy.ndimage.uniform_filter1d and Otsu)
# ---------------------------------------------------------------------------
def _smooth(p, k):
    p = np.asarray(p, float)
    fill = np.nanmean(p) if np.isfinite(p).any() else 0.0
    p = np.nan_to_num(p, nan=fill)
    k = max(int(round(k)), 1)
    if k < 2 or len(p) < 3:
        return p
    pad = k // 2
    q = np.pad(p, (pad, k - 1 - pad), mode="edge")
    return np.convolve(q, np.ones(k) / k, mode="valid")


def _otsu(x, bins=256):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < 10:
        return float(np.nanmedian(x)) if x.size else 0.0
    lo, hi = np.percentile(x, 0.5), np.percentile(x, 99.5)
    if hi <= lo:
        return float(lo)
    hist, edges = np.histogram(x, bins=bins, range=(lo, hi))
    mids = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * mids) / np.maximum(w0, 1)
    m1 = (np.sum(hist * mids) - np.cumsum(hist * mids)) / np.maximum(w1, 1)
    return float(mids[np.argmax(w0 * w1 * (m0 - m1) ** 2)])


def _row_signal(profile, step):
    """detrended, unit-variance ridge signal of an across profile (rows = peaks).
    step = profile bin size in mm (PhenoFuse: 4 cm smoothing, 0.5 m detrend)."""
    p = _smooth(profile, 40.0 / step)
    det = p - _smooth(p, 500.0 / step)
    sd = np.nanstd(det)
    return det / sd if sd > 0 else det


def row_pitch_from_profile(signal, step, lo, hi):
    """(pitch, power share) of the strongest periodicity of a ridge signal with a
    period between lo and hi (same units as step)."""
    s = np.nan_to_num(np.asarray(signal, float))
    s = s - s.mean()
    if len(s) < 20:
        return np.nan, 0.0
    n = max(4096, 1 << int(math.ceil(math.log2(len(s) * 4))))
    F = np.abs(np.fft.rfft(s * np.hanning(len(s)), n)) ** 2
    fq = np.fft.rfftfreq(n, d=step)
    band = (fq >= 1.0 / hi) & (fq <= 1.0 / lo)
    if not band.any() or F[band].sum() <= 0:
        return np.nan, 0.0
    k = int(np.argmax(np.where(band, F, 0)))
    return float(1.0 / fq[k]), float(F[k] / F[band].sum())


def rows_from_profile(v, sig, lo, hi, pitch, step, edge_gap):
    """Comb of rows at `pitch` filling [lo, hi] (PhenoFuse rows_from_profile, with the
    pitch already chosen). v = bin centres, sig = ridge signal (rows = peaks)."""
    n_fit = int(round((hi - lo - pitch / 2) / pitch)) + 1
    if n_fit < 1 or not np.isfinite(pitch) or pitch <= 0:
        return None
    half_w = max(step, 0.12 * pitch)

    def val(pos):
        m = np.abs(v - pos) < half_w
        return float(np.nanmean(sig[m])) if m.any() else np.nan

    centre = (lo + hi) / 2
    half = (n_fit - 1) * pitch / 2
    best = (None, -np.inf)
    for off in np.arange(-pitch / 2, pitch / 2, max(step / 2, pitch / 200)):
        teeth = centre + off + np.arange(n_fit) * pitch - half
        if teeth[0] < lo - edge_gap or teeth[-1] > hi + edge_gap:
            continue
        s = np.nansum([val(t) for t in teeth]) - np.nansum([val(t + pitch / 2) for t in teeth[:-1]])
        if s > best[1]:
            best = (teeth, s)
    teeth = best[0]
    if teeth is None:
        return None
    contrast = np.array([val(t) - np.nanmean([val(t - pitch / 2), val(t + pitch / 2)]) for t in teeth])
    present = np.isfinite(contrast) & (contrast > 0)
    return dict(teeth=teeth, contrast=contrast, present=present)


# ---------------------------------------------------------------------------
# main entry
# ---------------------------------------------------------------------------
GREEN_FLOOR = 0.05      # chromatic ExG of clearly green plant pixels (soil is ~ -0.1..0)


def detect_rows(arr, across_axis=None, mm_px=None, pitch_mm=None, rows_expected=0,
                max_tilt_deg=8.0, work_px=1200):
    """Find the crop rows of one plot crop (see _detect_core for the arguments).

    At emergence the plants are a fraction of a percent of the pixels and the Otsu
    split of ExG falls inside the SOIL colours: on ridged seedbeds (NUE Northam
    2026-07-03) the lighter ridges then read as 'vegetation' and the comb lands half a
    pitch off, between the plant rows. So when the Otsu threshold is below
    GREEN_FLOOR the rows are also found from the clearly green pixels only, and the
    result whose rows sit on more green (green_contrast) is kept. Where the plants
    dominate the split, both agree (York 2025-06-25 is unchanged)."""
    kw = dict(across_axis=across_axis, mm_px=mm_px, pitch_mm=pitch_mm,
              rows_expected=rows_expected, max_tilt_deg=max_tilt_deg, work_px=work_px)
    d0 = _detect_core(arr, **kw)
    if not (np.isfinite(d0.get("veg_threshold", np.nan)) and d0["veg_threshold"] < GREEN_FLOOR):
        return d0
    d1 = _detect_core(arr, green_floor=GREEN_FLOOR, **kw)
    if d1.get("ok") and (not d0.get("ok") or
                         d1.get("green_contrast", -1) > d0.get("green_contrast", -1)):
        return d1
    return d0


def _detect_core(arr, across_axis=None, mm_px=None, pitch_mm=None, rows_expected=0,
                 max_tilt_deg=8.0, work_px=1200, green_floor=None):
    """Find the crop rows of one plot crop.

    arr            3 x H x W uint8 (rasterio band order). Black (all bands <= 4) = no data
                   (the outside of a masked crop).
    across_axis    0 = rows run along x (across = image y), 1 = rows run along y
                   (across = image x). None: the short image side is across.
    mm_px          ground sampling of `arr` (mm per pixel), if known.
    pitch_mm       seeder row spacing, if known (then only the phase is fitted).
    rows_expected  sown rows per plot (0 = take the run of rows that are present).
    green_floor    None: vegetation = ExG above the plot's Otsu threshold, ridge signal =
                   share + mean ExG. A number: vegetation = ExG above it and the ridge
                   signal uses only the green excess (soil colour ignored).

    Returns a dict in FULL-resolution pixels of `arr`:
      ok, reason, across_axis, angle_deg (rotation of the rows from the image axis),
      pitch_px, rows (tooth positions across, measured on the centre line along),
      contrast, present, plot_rows (indices of the plot's rows into `rows`),
      along_centre (pixel position along of the reference line), veg_threshold,
      pitch_share (how periodic the profile is, 0-1).
    """
    C, H, W = arr.shape
    if across_axis is None:
        across_axis = 0 if H < W else 1
    out = dict(ok=False, reason="", across_axis=across_axis, angle_deg=0.0, pitch_px=np.nan,
               rows=[], contrast=[], present=[], plot_rows=[], along_centre=None,
               veg_threshold=np.nan, pitch_share=0.0,
               veg_mode="green" if green_floor is not None else "otsu")
    # work at <= work_px on the across side; one work pixel = dec source pixels
    dec = max(1, int(round(min(H, W) / float(work_px))))
    a = arr[:3, ::dec, ::dec].astype(np.float32)
    R, G, B = a[0], a[1], a[2]
    tot = R + G + B
    valid = a.max(axis=0) > 4
    if valid.sum() < 500:
        out["reason"] = "almost no image data"
        return out
    with np.errstate(all="ignore"):
        exg = np.where(valid & (tot > 0), (2 * G - R - B) / np.maximum(tot, 1e-6), np.nan)
    otsu = _otsu(exg[valid])
    thr = otsu if green_floor is None else float(green_floor)
    out["veg_threshold"] = otsu
    veg = valid & (exg > thr)
    green = valid & (exg > GREEN_FLOOR)
    h, w = valid.shape
    # pixel coordinates: across (a) and along (l), in work pixels
    yy, xx = np.nonzero(valid)
    if across_axis == 0:
        acr, alo = yy.astype(float), xx.astype(float)
        n_across, n_along = h, w
    else:
        acr, alo = xx.astype(float), yy.astype(float)
        n_across, n_along = w, h
    ev = exg[yy, xx]
    if green_floor is not None:                          # soil colour must not count
        ev = np.clip(ev - green_floor, 0.0, None)
    vv = veg[yy, xx].astype(float)
    gv = green[yy, xx].astype(float)
    lc = n_along / 2.0

    mm_work = (mm_px * dec) if mm_px else None          # mm per work pixel
    if pitch_mm and mm_work:
        p_lo = p_hi = pitch_mm / mm_work
    elif mm_work:
        p_lo, p_hi = ROW_PITCH_RANGE_MM[0] / mm_work, ROW_PITCH_RANGE_MM[1] / mm_work
    else:                                                # no scale: 3-14 rows across the crop
        p_lo, p_hi = n_across / 14.0, n_across / 3.0
    step = 1.0                                           # profile bin = 1 work pixel
    # mm-based smoothing in _row_signal needs mm per bin; without a scale assume the
    # pitch range midpoint is ~240 mm
    mm_bin = mm_work if mm_work else 240.0 / ((p_lo + p_hi) / 2.0)

    def profiles(ang_deg, A, L, E, V):
        t = math.tan(math.radians(ang_deg))
        u = A - t * (L - lc)                             # across position, rotated
        nb = n_across
        idx = np.floor(u / step).astype(np.int64)
        ok = (idx >= 0) & (idx < nb)
        cnt = np.bincount(idx[ok], minlength=nb).astype(float)
        sv = np.bincount(idx[ok], weights=V[ok], minlength=nb)
        se = np.bincount(idx[ok], weights=np.nan_to_num(E[ok]), minlength=nb)
        with np.errstate(all="ignore"):
            share = np.where(cnt > 3, sv / cnt, np.nan)
            mexg = np.where(cnt > 3, se / cnt, np.nan)
        return share, mexg, cnt

    # --- pitch and tilt from short segments along the plot -------------------
    # A tilt smears the whole-plot profile, but over a short segment the rows are
    # still sharp. Pitch = peak of the segments' summed power spectra (the phase
    # does not matter there, so the tilt cannot bias it); tilt = drift of the row
    # phase from segment to segment, fitted as a straight line (the same idea as
    # PhenoFuse's strip-long row lines: rows chained across banks, u = a + b v).
    x = (np.arange(n_across) + 0.5) * step
    K = int(np.clip(n_along / 250.0, 6, 24))
    seg = np.minimum((alo / n_along * K).astype(int), K - 1)
    segs = []
    for k in range(K):
        m = seg == k
        if m.sum() < 200:
            continue
        sh_k, me_k, _ = profiles(0.0, acr[m], alo[m], ev[m], vv[m])
        okk = np.isfinite(sh_k)
        if okk.sum() < 20:
            continue
        s_k = np.nan_to_num(_row_signal(sh_k, mm_bin) + _row_signal(me_k, mm_bin)) * okk
        segs.append((alo[m].mean(), s_k))
    if len(segs) < 3:
        out["reason"] = "too little image data along the plot"
        return out
    if p_lo == p_hi:
        P, share_pw = float(p_lo), 1.0
    else:
        nfft = 1 << 16
        fq = np.fft.rfftfreq(nfft, d=step)
        band = (fq >= 1.0 / p_hi) & (fq <= 1.0 / p_lo)
        hann = np.hanning(n_across)
        Fs = sum(np.abs(np.fft.rfft(s_k * hann, nfft)) ** 2 for _, s_k in segs)
        if not band.any() or Fs[band].sum() <= 0:
            out["reason"] = "no row periodicity found"
            return out
        kk = int(np.argmax(np.where(band, Fs, 0)))
        P, share_pw = float(1.0 / fq[kk]), float(Fs[kk] / Fs[band].sum())

    def seg_phase(P_):
        e = np.exp(-2j * np.pi * x / P_)
        ls, pos, amp = [], [], []
        for l_k, s_k in segs:
            z = (s_k * e).sum()
            ls.append(l_k); amp.append(abs(z) / len(s_k))
            pos.append((-np.angle(z) / (2 * np.pi) * P_) % P_)      # row position mod P
        return np.array(ls), np.array(pos), np.array(amp)

    ang = 0.0
    if max_tilt_deg > 0:
        ls, pos, amp = seg_phase(P)
        pu = pos.copy()
        for i in range(1, len(pu)):                # unwrap along the plot
            d = ((pos[i] - pu[i - 1] + P / 2) % P) - P / 2
            pu[i] = pu[i - 1] + d
        wts = amp.copy()
        sol = np.array([np.median(pu), 0.0])
        for _it in range(4):                       # weighted LS, outliers dropped
            sw = np.sqrt(np.maximum(wts, 0))
            if (sw > 0).sum() < 3:
                break
            A = np.c_[np.ones(len(ls)), ls - lc]
            sol, *_ = np.linalg.lstsq(A * sw[:, None], pu * sw, rcond=None)
            r_ = pu - A @ sol
            tol = max(P / 8.0, 2.5 * float(np.median(np.abs(r_[wts > 0]))))
            wts = amp * (np.abs(r_) < tol)
        ang = math.degrees(math.atan(float(sol[1])))
        if abs(ang) > max_tilt_deg:
            ang = 0.0
    # fine pitch on the de-tilted whole-plot profile (the FFT bin is coarse)
    share, mexg, cnt = profiles(ang, acr, alo, ev, vv)
    if p_lo != p_hi:
        sig0 = np.nan_to_num(_row_signal(share, mm_bin) + _row_signal(mexg, mm_bin)) * np.isfinite(share)
        best_p = (P, -1.0)
        for P_ in P * np.linspace(0.94, 1.06, 241):
            pw = abs((sig0 * np.exp(-2j * np.pi * x / P_)).sum())
            if pw > best_p[1]:
                best_p = (float(P_), pw)
        P = best_p[0]
    sig = _row_signal(share, mm_bin) + _row_signal(mexg, mm_bin)
    v = (np.arange(n_across) + 0.5) * step
    okb = np.isfinite(share)
    lo, hi = v[okb][0], v[okb][-1]
    res = rows_from_profile(v, sig, lo, hi, P, step, edge_gap=0.25 * P)
    if res is None:
        out["reason"] = "comb fit failed"
        return out
    teeth, contrast, present = res["teeth"], res["contrast"], res["present"]
    n = len(teeth)
    # a tooth the comb placed beyond the image data is no row: its "contrast" is an
    # edge effect of the detrended profile, and it can be the largest of all
    in_data = (teeth >= lo) & (teeth <= hi)
    present = present & in_data
    # the plot's own rows
    if rows_expected and rows_expected <= int(in_data.sum()):
        best_k, best_s = None, -np.inf
        for k in range(0, n - rows_expected + 1):
            if not in_data[k:k + rows_expected].all():
                continue
            s = np.nansum(contrast[k:k + rows_expected])
            if s > best_s:
                best_k, best_s = k, s
        plot_rows = list(range(best_k, best_k + rows_expected)) if best_k is not None else []
    else:
        # a tooth within 0.3 pitch of the image-data edge is a row the plot
        # boundary cuts through (usually shared with the neighbour plot): it is
        # kept in `rows` but is not one of this plot's rows
        inner = (teeth > lo + 0.3 * P) & (teeth < hi - 0.3 * P)
        idx = np.nonzero(present & inner)[0]
        plot_rows = list(range(int(idx[0]), int(idx[-1]) + 1)) if len(idx) else []
    if not plot_rows:
        out["reason"] = "no row with positive contrast"
        return out
    # how well the rows sit on clearly green pixels (used to pick the otsu / green
    # result): green share at the plot's teeth minus at the mid-gaps
    gshare, _, _ = profiles(ang, acr, alo, ev, gv)
    gs = np.nan_to_num(_smooth(gshare, max(1.0, 0.12 * P)))

    def gval(pos):
        i = int(round(pos / step - 0.5))
        return float(gs[i]) if 0 <= i < len(gs) else 0.0
    gc = [gval(teeth[k]) - 0.5 * (gval(teeth[k] - P / 2) + gval(teeth[k] + P / 2))
          for k in plot_rows]
    f = float(dec)
    out.update(ok=True, green_contrast=float(np.mean(gc)) if gc else -1.0,
               angle_deg=float(ang), pitch_px=float(P * f),
               rows=[float(t * f) for t in teeth], contrast=[float(c) for c in contrast],
               present=[bool(p) for p in present], plot_rows=plot_rows,
               along_centre=float(lc * f), pitch_share=float(share_pw),
               span_centre=float((lo + hi) / 2.0 * f),
               n_rows_present=int(np.sum(present[plot_rows[0]:plot_rows[-1] + 1])))
    return out


def middle_band(det, n_middle):
    """(lo, hi, row indices) across-pixel band holding the n_middle central rows of the
    plot, measured on the reference line (det['along_centre']); half a pitch is added
    outside the first and last chosen row, so a tile spanning the band holds exactly
    those rows. None when fewer rows were found."""
    if not det.get("ok"):
        return None
    pr = det["plot_rows"]
    if n_middle <= 0 or len(pr) < n_middle:
        return None
    rows = det["rows"]
    # choose the block whose centre is closest to the centre of the plot's rows; with an
    # odd/even mismatch (e.g. 3 of 6) the block nearer the plot centre wins
    centre = (rows[pr[0]] + rows[pr[-1]]) / 2.0
    span_c = det.get("span_centre", centre)
    best = None
    for k in range(0, len(pr) - n_middle + 1):
        blk = pr[k:k + n_middle]
        c = (rows[blk[0]] + rows[blk[-1]]) / 2.0
        key = (round(abs(c - centre), 3), abs(c - span_c))   # tie -> nearer the crop centre
        if best is None or key < best[0]:
            best = (key, blk)
    blk = best[1]
    P = det["pitch_px"]
    return rows[blk[0]] - P / 2.0, rows[blk[-1]] + P / 2.0, blk


def across_at(det, pos_across, pos_along):
    """Across coordinate of a line parallel to the rows that passes the reference line
    at `pos_across`, evaluated at `pos_along` (handles the tilt)."""
    t = math.tan(math.radians(det.get("angle_deg", 0.0)))
    return pos_across + t * (pos_along - det["along_centre"])


def draw_rows(rgb_hwc, det, band=None, max_side=1600):
    """QA picture: the plot with every detected row (green = present, red = missing),
    the plot's rows thicker, and the middle band shaded. Returns a PIL image."""
    from PIL import Image, ImageDraw
    im = Image.fromarray(rgb_hwc)
    W, H = im.size
    s = min(1.0, max_side / float(max(W, H)))
    if s < 1.0:
        im = im.resize((max(1, int(W * s)), max(1, int(H * s))), Image.BILINEAR)
    if not det.get("ok"):
        return im
    ov = Image.new("RGBA", im.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(ov)
    ax = det["across_axis"]
    L = W if ax == 0 else H

    def line(pos, col, width):
        a0 = across_at(det, pos, 0) * s
        a1 = across_at(det, pos, L) * s
        if ax == 0:
            d.line([(0, a0), (L * s, a1)], fill=col, width=width)
        else:
            d.line([(a0, 0), (a1, L * s)], fill=col, width=width)

    if band is not None:
        lo, hi = band[0], band[1]
        a = [across_at(det, lo, 0) * s, across_at(det, hi, 0) * s,
             across_at(det, hi, L) * s, across_at(det, lo, L) * s]
        if ax == 0:
            poly = [(0, a[0]), (0, a[1]), (L * s, a[2]), (L * s, a[3])]
        else:
            poly = [(a[0], 0), (a[1], 0), (a[2], L * s), (a[3], L * s)]
        d.polygon(poly, fill=(255, 210, 0, 60), outline=(255, 210, 0, 255))
    pr = set(det["plot_rows"])
    for i, (r, p) in enumerate(zip(det["rows"], det["present"])):
        col = (0, 255, 0, 230) if p else (255, 40, 40, 230)
        line(r, col, 3 if i in pr else 1)
    return Image.alpha_composite(im.convert("RGBA"), ov).convert("RGB")
