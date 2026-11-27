"""plot_detect.py - find the plots of a field trial in an RGB orthomosaic.

A trial is a regular layout: ranges stacked along the crop rows, separated by bare
alleys, and in every range plots side by side across the rows, separated by a
wider gap than the row spacing. The plots are found in three passes, all on one
decimated read of the orthomosaic (pure numpy; the PyInstaller build has no scipy):

  1 direction   vegetation = excess green on chromatic coordinates above one Otsu
                threshold (never below GREEN_FLOOR, so seedlings on bare soil work).
                The crop rows give the direction: the angle whose across-row
                profiles have the most power at seeder spacings (15-36 cm), summed
                over short segments so slightly bent rows and plots drilled with a
                different row phase still add up. Coarse 1 deg search, then 0.1 and
                0.01 deg. With 4 clicked corners only +-4 deg around them is searched.
  2 layout      the vegetation is resampled into the trial frame (a = along the rows,
                c = across). Along a, the share of vegetation drops in the alleys;
                across c, after smoothing over one row spacing, it drops in the gaps
                between plots. Crop = above half the plateau (half-height edges: for
                rows this is half a row spacing outside the outer rows, the same plot
                edge as 'width = rows x spacing'). Runs that are a whole number of
                pitches long are split, slivers are dropped.
  3 each plot   every plot's own sides and ends are measured the same way in a window
                around it (the alleys are not always straight, plots in a column can
                drift a few cm). All plots get the median size; each plot keeps its
                own measured centre, checked against the straight line through its
                column. Each plot gets a confidence 0-1 and the reason when it is low:
                unclear side / end, little crop (empty or failed plot), or out of line
                with its neighbours.

Works best a few weeks after emergence up to just before the canopy closes; it needs
bare gaps between plots across the rows (a trial drilled continuously across has
none: then give 'Plots per range' or the plot width and the plots are laid out at
that pitch inside the detected crop).

Public API
----------
detect_plots(ortho_path, quad=None, n_ranges=0, n_cols=0, plot_w=0, plot_l=0,
             log=print, progress=None, cancel=None) -> dict (see the function)
write_report(path, result)                     -> per-plot confidence CSV
"""
from __future__ import annotations

import csv
import math
import time

import numpy as np

from row_detect import GREEN_FLOOR, _otsu, _smooth

ROW_BAND_M = (0.15, 0.36)       # seeder row spacings looked for
STEP_M = 0.02                    # working sample (2 cm)
MAX_SAMPLES = 24e6               # cap on the decimated read (memory / speed)
FLAG_BELOW = 0.5                 # confidence below this = "check this plot"


class DetectError(Exception):
    """A layout could not be found; the message says what to give instead."""


# ---------------------------------------------------------------------------
# reading
# ---------------------------------------------------------------------------
def _metres_per_unit(crs, lat):
    try:
        geo = bool(crs is not None and crs.is_geographic)
    except Exception:
        geo = False
    if not geo:
        return 1.0, 1.0
    return 111320.0 * math.cos(math.radians(lat)), 110540.0


def read_field(path, quad=None, step_m=STEP_M, max_samples=MAX_SAMPLES):
    """One decimated read of the orthomosaic (inside the quad's bounding box when
    given). Returns dict(veg bool HxW, valid bool HxW, transform, mx, my, step, crs)."""
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.windows import from_bounds, Window
    with rasterio.open(path) as ds:
        b = ds.bounds
        if quad:
            xs, ys = [p[0] for p in quad], [p[1] for p in quad]
            px, py = 0.05 * (max(xs) - min(xs)), 0.05 * (max(ys) - min(ys))
            x0, x1 = max(min(xs) - px, b.left), min(max(xs) + px, b.right)
            y0, y1 = max(min(ys) - py, b.bottom), min(max(ys) + py, b.top)
            if not (x1 > x0 and y1 > y0):
                raise DetectError("the 4 corners are outside the orthomosaic")
        else:
            x0, y0, x1, y1 = b.left, b.bottom, b.right, b.top
        mx, my = _metres_per_unit(ds.crs, (y0 + y1) / 2.0)
        px_m = min(abs(ds.res[0]) * mx, abs(ds.res[1]) * my)
        area = (x1 - x0) * mx * (y1 - y0) * my
        step = max(step_m, px_m, math.sqrt(area / max_samples))
        win = from_bounds(x0, y0, x1, y1, ds.transform)
        c0, r0 = max(0, int(math.floor(win.col_off))), max(0, int(math.floor(win.row_off)))
        c1 = min(ds.width, int(math.ceil(win.col_off + win.width)))
        r1 = min(ds.height, int(math.ceil(win.row_off + win.height)))
        win = Window(c0, r0, c1 - c0, r1 - r0)
        ow = max(8, int(round((c1 - c0) * abs(ds.res[0]) * mx / step)))
        oh = max(8, int(round((r1 - r0) * abs(ds.res[1]) * my / step)))
        bands = [1, 2, 3] if ds.count >= 3 else [1, 1, 1]
        rgb = ds.read(bands, window=win, out_shape=(3, oh, ow),
                      resampling=Resampling.average).astype(np.float32)
        try:
            msk = ds.read_masks(1, window=win, out_shape=(oh, ow)) > 0
        except Exception:
            msk = np.ones((oh, ow), bool)
        tr = ds.window_transform(win) * rasterio.Affine.scale((c1 - c0) / ow, (r1 - r0) / oh)
        crs = ds.crs
    s = rgb.sum(0)
    valid = msk & (s > 12)
    exg = np.where(valid, (2 * rgb[1] - rgb[0] - rgb[2]) / np.maximum(s, 1.0), np.nan)
    ev = exg[valid]
    thr = _otsu(ev) if ev.size else GREEN_FLOOR
    thr = max(thr, GREEN_FLOOR)
    veg = valid & (exg > thr)
    # graded greenness 0..1: 0 at the median pixel (soil while the crop covers less
    # than half; inside the crop later - the gaps then stay 0 either way), 1 at the
    # greenest. A hard threshold keeps only the brightest seedlings and breaks plots.
    # The zero is the LOCAL median (2 m blocks): soil colour drifts across a field
    # (moisture, sun angle, mosaic seams) by more than seedlings add.
    med = float(np.median(ev)) if ev.size else 0.0
    hi = max(float(np.percentile(ev, 99.5)) if ev.size else 1.0, med + 0.08)
    bg = _block_median(exg, max(4, int(round(2.0 / step))))
    bg = np.where(np.isfinite(bg), bg, med)
    score = np.clip((np.nan_to_num(exg, nan=med) - bg) / (hi - med), 0, 1).astype(np.float32)
    score[~valid] = np.nan
    lum = (s / 3.0).astype(np.float32)
    lum[~valid] = np.nan
    return dict(veg=veg, score=score, lum=lum, valid=valid, transform=tr, mx=mx, my=my, step=step,
                crs=crs, exg_threshold=float(thr),
                cover=float(veg.sum() / max(valid.sum(), 1)))


def _block_median(a, b):
    """median of a (nan-aware) in b x b blocks, bilinearly expanded to a's shape"""
    H, W = a.shape
    nh, nw = int(math.ceil(H / b)), int(math.ceil(W / b))
    pad = np.full((nh * b, nw * b), np.nan, np.float32)
    pad[:H, :W] = a
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        m = np.nanmedian(pad.reshape(nh, b, nw, b).transpose(0, 2, 1, 3).reshape(nh, nw, b * b),
                         axis=2)
    if not np.isfinite(m).any():
        return np.full(a.shape, np.nan, np.float32)
    m = np.where(np.isfinite(m), m, np.nanmedian(m))
    # bilinear expansion: block centres at (i + 0.5) * b
    yi = np.clip((np.arange(H) + 0.5) / b - 0.5, 0, nh - 1)
    xi = np.clip((np.arange(W) + 0.5) / b - 0.5, 0, nw - 1)
    y0, x0 = np.floor(yi).astype(int), np.floor(xi).astype(int)
    y1, x1 = np.minimum(y0 + 1, nh - 1), np.minimum(x0 + 1, nw - 1)
    fy, fx = (yi - y0)[:, None], (xi - x0)[None, :]
    top = m[y0][:, x0] * (1 - fx) + m[y0][:, x1] * fx
    bot = m[y1][:, x0] * (1 - fx) + m[y1][:, x1] * fx
    return (top * (1 - fy) + bot * fy).astype(np.float32)


def _pixel_xy_m(F, rows, cols):
    """metric coordinates of raster pixel centres"""
    tr = F["transform"]
    x = tr.c + (cols + 0.5) * tr.a + (rows + 0.5) * tr.b
    y = tr.f + (cols + 0.5) * tr.d + (rows + 0.5) * tr.e
    return x * F["mx"], y * F["my"]


# ---------------------------------------------------------------------------
# 1 direction of the rows
# ---------------------------------------------------------------------------
def _row_power(xy, wt, theta, seg, step, band):
    """summed band power of the across-row profiles of segments seg long"""
    ca, sa = math.cos(theta), math.sin(theta)
    u = xy[:, 0] * ca + xy[:, 1] * sa
    v = -xy[:, 0] * sa + xy[:, 1] * ca
    iu = ((u - u.min()) / seg).astype(np.int64)
    iv = ((v - v.min()) / step).astype(np.int64)
    nu, nv = int(iu.max()) + 1, int(iv.max()) + 1
    H = np.bincount(iu * nv + iv, weights=wt, minlength=nu * nv).reshape(nu, nv).astype(np.float32)
    H -= H.mean(axis=1, keepdims=True)
    n = 1 << int(math.ceil(math.log2(max(nv, 16))))
    P = np.abs(np.fft.rfft(H, n, axis=1)) ** 2
    f = np.fft.rfftfreq(n, d=step)
    sel = (f >= 1.0 / band[1]) & (f <= 1.0 / band[0])
    return float(P[:, sel].sum())


def find_direction(F, around=None, log=print):
    """angle (radians, metric x-y frame) of the crop rows, and how clear it was."""
    sc = F["score"]
    rr, cc = np.nonzero(sc > 0.05)
    if rr.size < 500:
        raise DetectError("almost no green crop in the orthomosaic (%.2f %% of the pixels)"
                          % (100 * F["cover"]))
    rng = np.random.default_rng(0)
    if rr.size > 1_500_000:
        k = rng.choice(rr.size, 1_500_000, replace=False)
        rr, cc = rr[k], cc[k]
    wt = sc[rr, cc].astype(np.float64)
    x, y = _pixel_xy_m(F, rr, cc)
    xy = np.c_[x - x.mean(), y - y.mean()]
    step = max(F["step"], 0.01)
    if around is None:
        stages = [(np.radians(np.arange(0.0, 180.0, 1.0)), 1.5),
                  (None, 0.1, 1.5, 5.0), (None, 0.01, 0.15, 10.0)]
    else:
        stages = [(np.radians(np.arange(-4.0, 4.01, 0.25)) + around, 2.0),
                  (None, 0.05, 0.4, 5.0), (None, 0.01, 0.08, 10.0)]
    best, clarity = None, 0.0
    for st in stages:
        if st[0] is not None:
            angs, seg = st
        else:
            _, inc, half, seg = st
            angs = best + np.radians(np.arange(-half, half + inc / 2, inc))
        pw = np.array([_row_power(xy, wt, t, seg, step, ROW_BAND_M) for t in angs])
        best = float(angs[int(np.argmax(pw))])
        if st[0] is not None:
            clarity = float(pw.max() / max(np.median(pw), 1e-9))
    return best % math.pi, clarity


# ---------------------------------------------------------------------------
# 2 trial frame + layout
# ---------------------------------------------------------------------------
def _resample(F, O, ah, ch, a_rng, c_rng, step, key="score"):
    """F[key] (greenness, or 'lum' brightness) in the trial frame: float32
    (na x nc), nan = no data"""
    tr = F["transform"]
    inv = ~tr
    a0, a1 = a_rng
    c0, c1 = c_rng
    na, nc = int((a1 - a0) / step) + 1, int((c1 - c0) / step) + 1
    out = np.full((na, nc), np.nan, np.float32)
    C = c0 + (np.arange(nc) + 0.5) * step
    H, W = F[key].shape
    for i0 in range(0, na, 512):
        A = a0 + (np.arange(i0, min(na, i0 + 512)) + 0.5) * step
        AA, CC = np.meshgrid(A, C, indexing="ij")
        X = (O[0] + AA * ah[0] + CC * ch[0]) / F["mx"]
        Y = (O[1] + AA * ah[1] + CC * ch[1]) / F["my"]
        col, row = inv * (X, Y)
        col, row = np.floor(col).astype(np.int64), np.floor(row).astype(np.int64)
        ok = (col >= 0) & (col < W) & (row >= 0) & (row < H)
        blk = np.full(AA.shape, np.nan, np.float32)
        blk[ok] = F[key][row[ok], col[ok]]
        out[i0:i0 + blk.shape[0]] = blk
    return out


def _share(G, axis):
    """mean greenness per line of the trial raster (nan where mostly no data)"""
    valid = np.isfinite(G)
    n = valid.sum(axis=axis).astype(float)
    s = np.nansum(G, axis=axis, dtype=np.float64) if G.size else np.zeros(0)
    return np.where(n > 0.2 * G.shape[axis], s / np.maximum(n, 1), np.nan)


def _runs(mask):
    """(start, end) index pairs of True runs, end exclusive"""
    m = np.r_[False, np.asarray(mask, bool), False]
    d = np.diff(m.astype(np.int8))
    return list(zip(np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]))


def _envelope(p, win):
    """local plateau: running 90th percentile over win bins (block-decimated)"""
    n = len(p)
    d = max(1, int(win) // 40)
    nb = int(math.ceil(n / d))
    q = np.full(nb * d, np.nan)
    q[:n] = p
    with np.errstate(all="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            b = np.nanmean(q.reshape(nb, d), axis=1)
            k = max(3, int(win) // d) | 1
            pad = np.r_[np.full(k // 2, np.nan), b, np.full(k // 2, np.nan)]
            w = np.lib.stride_tricks.sliding_window_view(pad, k)
            env = np.nanpercentile(w, 90, axis=1)
    return np.repeat(env, d)[:n]


def _fill_nan(v):
    v = np.asarray(v, float)
    ok = np.isfinite(v)
    if ok.all() or not ok.any():
        return np.nan_to_num(v)
    return np.interp(np.arange(len(v)), np.nonzero(ok)[0], v[ok])


def _normalise(p, step, env_m=4.0):
    """profile / local plateau (running 90th percentile): 0 = bare, ~1 = crop"""
    env = _envelope(p, env_m / step)
    ok = np.isfinite(p) & np.isfinite(env) & (env > 0)
    q = np.full(len(p), np.nan)
    q[ok] = np.clip(p[ok] / env[ok], 0, 1.5)
    return q


def _line_phase(Sn, nph, loc_frac=0.08):
    """phase (in samples, unwrapped) of every profile: the best straight line
    through all profiles (offset + slope per profile) on the score table Sn
    (profiles x nph), then each profile within +-loc_frac of the line."""
    npf = Sn.shape[0]
    if npf < 3:
        return np.array([float(np.argmax(Sn[j])) for j in range(npf)])
    best = (-np.inf, 0, 0)
    jj = np.arange(npf)
    for sl in range(-(nph // 4), nph // 4 + 1):
        sh = (np.arange(nph)[None, :] + sl * jj[:, None]) % nph
        tot = Sn[jj[:, None], sh].sum(axis=0)
        t = int(np.argmax(tot))
        if tot[t] > best[0]:
            best = (tot[t], t, sl)
    _, t0, sl = best
    line = (t0 + sl * np.arange(npf)).astype(float)
    loc = max(1, int(loc_frac * nph))
    out = np.empty(npf)
    for j in range(npf):
        c = int(round(line[j]))
        cand = [(c + d) % nph for d in range(-loc, loc + 1)]
        k = max(cand, key=lambda q: Sn[j, q])
        d = ((k - c + nph // 2) % nph) - nph // 2
        out[j] = line[j] + d
    return out


def _boxcar_comb(qs, step, lo_m, hi_m):
    """_boxcar_fit at the autocorrelation pitch and at 1/2 and 1/3 of it (a profile
    that varies range by range - cloud shadow, two treatments - can make twice the
    pitch the stronger peak)."""
    P, acs = _pitch_ac(qs, step, lo_m, hi_m)
    if not np.isfinite(P):
        return None
    fits = []
    for d in (1, 2, 3):
        if P / d < lo_m:
            break
        r = _boxcar_fit(qs, step, P / d, acs)
        if r is not None:
            fits.append(r)
    if not fits:
        return None
    # the score is a mean per box, so a comb on every other range scores as well as
    # one on all of them: take the shortest pitch that scores near the best
    top = max(r["score"] for r in fits)
    return [r for r in fits if r["score"] >= 0.75 * top][-1]


def _boxcar_fit(qs, step, P, acs):
    """The ranges as a repeating box: plots of length L every P along each column
    (one P and L for the trial, a phase per column on a straight line through the
    columns, each plot then nudged on its own). The score is mean signal inside the
    box minus outside it, so both crop ends count and whatever grows in the middle of
    an alley (a drilled strip, weeds on the wheel track) does not matter.
    Returns None or dict(pitch, size, ph, score, cells=[(k, j, start, end, mean,
    contrast)])."""
    n = len(qs[0])
    Pp = P / step
    Cs = [np.r_[0.0, np.cumsum(np.nan_to_num(q))] for q in qs]
    Cn = [np.r_[0.0, np.cumsum(np.isfinite(q))] for q in qs]

    def wmean(j, a, b):
        a = np.clip(a, 0, n)
        b = np.clip(b, 0, n)
        cnt = Cn[j][b] - Cn[j][a]
        with np.errstate(all="ignore"):
            m = (Cs[j][b] - Cs[j][a]) / cnt
        return np.where(cnt >= 0.6 * np.maximum(b - a, 1), m, np.nan)

    nph = max(8, int(round(Pp)))
    K = int(n / Pp) + 3
    st = (np.arange(nph)[:, None] + np.arange(-1, K)[None, :] * Pp)   # nph x K
    import warnings

    def table(Ls):
        # inside the box minus two narrow strips just outside its ends
        gs = max(2.0, min(0.3 / step, (Pp - Ls) / 2))
        S = np.full((len(qs), nph), np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            for j in range(len(qs)):
                a = st.astype(int)
                b = (st + Ls).astype(int)
                mi = wmean(j, a, b)
                mo = np.nanmean([wmean(j, (st - gs).astype(int), a),
                                 wmean(j, b, (st + Ls + gs).astype(int))], axis=0)
                S[j] = np.nanmean(mi - mo, axis=1)
        return np.nan_to_num(S, nan=-1.0)

    # The shortest box that scores within 10 % of the best: cutting into the crop
    # costs a lot, swallowing clutter at the plot end (a drilled strip in the alley,
    # a headland turn) costs little, so the longest box is not the plot.
    res = []
    for frac in np.arange(0.4, 0.975, 0.01):
        S = table(frac * Pp)
        res.append((frac, float(S.max(axis=1).sum()), S))
    top = max(r[1] for r in res)
    frac, tot, S = next(r for r in res if r[1] >= 0.9 * top)
    noise = float(np.nanmedian([np.nanstd(np.diff(q)) for q in qs])) or 1e-9
    score = tot / len(qs) / noise
    Ls = frac * Pp
    ph = _line_phase(S, nph) * step                       # start of plot k = 0, metres
    L = Ls * step
    g = max(2 * step, min(0.3, (P - L) / 2))
    cells = []
    for j, q in enumerate(qs):
        for k in range(-1, K):
            s0 = ph[j] + k * P
            if s0 + L < 0 or s0 > n * step:
                continue
            bestc = None
            for dlt in np.arange(-0.1 * P, 0.1 * P + step / 2, step):
                a, b = s0 + dlt, s0 + dlt + L
                ia, ib = int(a / step), int(b / step)
                mi = wmean(j, np.array([ia]), np.array([ib]))[0]
                with warnings.catch_warnings():          # both flanks outside the data
                    warnings.simplefilter("ignore", RuntimeWarning)
                    fl = np.nanmean([wmean(j, np.array([int((a - g) / step)]), np.array([ia]))[0],
                                     wmean(j, np.array([ib]), np.array([int((b + g) / step)]))[0]])
                if not (np.isfinite(mi) and np.isfinite(fl)):
                    continue
                sc_ = mi - fl
                if bestc is None or sc_ > bestc[0]:
                    bestc = (sc_, a, b, mi, fl)
            if bestc is None:
                continue
            sc_, a, b, mi, fl = bestc
            cells.append((k, j, a, b, mi, float(np.clip(sc_ / mi, 0, 1)) if mi > 0 else 0.0))
    return dict(pitch=P, size=L, ph=ph, cells=cells, ac=acs, score=score)


def _comb(qs, step, lo_m, hi_m, log=None):
    """The gaps between plots as a comb: one pitch for all profiles, one phase per
    profile (unwrapped, so neighbouring profiles may drift - staggered alleys), then
    every gap refined on its own (+-0.1 pitch; +-0.2 for a single profile).

    A gap is a NARROW feature repeating at the pitch. It is usually a dip (bare soil
    between plots) but can be a peak: at emergence a weedy or better-watered strip
    between plots can be greener than the seedlings (NUE Northam 2026-07-03). The
    filter is narrow mean minus pitch-wide mean; the sign with the stronger comb wins.

    qs: profiles (nan = no data), smoothed so the crop rows are gone.
    Returns None (no periodicity) or dict(pitch, sign, ph, ks, gaps (K x n), sf
    (feature strength per gap), strength (median), gw (dip width or nan), f)."""
    P, acs = _pitch_ac(qs, step, lo_m, hi_m)
    if not np.isfinite(P):
        return None
    g = min(max(0.12 * P, 3 * step), 0.6)
    kg, kP = max(1, int(round(g / step))), max(3, int(round(P / step)))
    n = len(qs[0])
    x = (np.arange(n) + 0.5) * step
    fs = []
    for q in qs:
        v = _fill_nan(q)
        f = _smooth(v, kg) - _smooth(v, kP)
        f[~np.isfinite(q)] = np.nan
        fs.append(f)
    nph = max(8, int(round(P / step)))
    K = int(n / (P / step)) + 2
    T = np.arange(nph)[:, None] + np.arange(K)[None, :] * (P / step)
    T = T.astype(int)
    inr = T < n
    S = np.full((len(fs), nph), np.nan)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for j, f in enumerate(fs):
            vals = np.where(inr, f[np.minimum(T, n - 1)], np.nan)
            S[j] = np.nanmean(vals, axis=1)
        pos = np.nansum(np.nanmax(S, axis=1))
        neg = np.nansum(-np.nanmin(S, axis=1))
    sign = 1.0 if pos > neg else -1.0
    # Alleys run in straight lines across the trial, even when staggered: fit one
    # line (offset + slope per profile) through all profiles at once, maximising the
    # summed comb score, then let each profile move a little off it. Phases picked
    # profile by profile lock half a pitch off where the crop is thin.
    Sn = np.nan_to_num(sign * S, nan=0.0)
    npf = len(fs)
    if npf >= 3:
        best = (-np.inf, 0, 0)
        jj = np.arange(npf)
        for sl in range(-(nph // 4), nph // 4 + 1):
            sh = (np.arange(nph)[None, :] + sl * jj[:, None]) % nph
            tot = Sn[jj[:, None], sh].sum(axis=0)
            t = int(np.argmax(tot))
            if tot[t] > best[0]:
                best = (tot[t], t, sl)
        _, t0, sl = best
        line = (t0 + sl * np.arange(npf)).astype(float)       # in phase samples
        loc = max(1, int(0.08 * nph))
        ph = np.empty(npf)
        for j in range(npf):
            c = int(round(line[j]))
            cand = [(c + d) % nph for d in range(-loc, loc + 1)]
            k = max(cand, key=lambda q: Sn[j, q])
            d = ((k - c + nph // 2) % nph) - nph // 2
            ph[j] = (line[j] + d) * (P / nph)
    else:
        ph = np.array([np.argmax(Sn[j]) * (P / nph) for j in range(npf)], float)
        ph[~np.isfinite(S).any(axis=1)] = np.nan
        ok = np.isfinite(ph)
        if not ok.any():
            return None
        ang = np.interp(np.arange(npf), np.nonzero(ok)[0], ph[ok] / P * 2 * np.pi)
        ph = np.unwrap(ang) / (2 * np.pi) * P
    kmin = int(math.floor(-ph.max() / P)) - 1
    kmax = int(math.ceil((n * step - ph.min()) / P)) + 1
    ks = np.arange(kmin, kmax + 1)
    gaps = np.full((len(ks), len(fs)), np.nan)
    sf = np.full((len(ks), len(fs)), np.nan)
    widths = []
    qsm = [_smooth(_fill_nan(q), max(1, kg // 2)) for q in qs]
    for j, f in enumerate(fs):
        for ii, k in enumerate(ks):
            t = ph[j] + k * P
            win = 0.2 * P if npf < 3 else 0.1 * P
            i0, i1 = max(int((t - win) / step), 0), min(int((t + win) / step), n)
            if i1 - i0 < 3 or np.isfinite(f[i0:i1]).sum() < 0.6 * (i1 - i0):
                continue
            w = sign * f[i0:i1]
            m = int(np.nanargmax(w))
            d = 0.0
            if 0 < m < len(w) - 1 and np.isfinite(w[m - 1]) and np.isfinite(w[m + 1]):
                den = w[m - 1] - 2 * w[m] + w[m + 1]
                d = 0.5 * (w[m - 1] - w[m + 1]) / den if den < 0 else 0.0
            gaps[ii, j] = x[i0 + m] + float(np.clip(d, -0.5, 0.5)) * step
            sf[ii, j] = w[m]
            if sign < 0:
                # a dip: its edges at half depth; the centre is midway between them,
                # which is well defined however wide and flat the alley is
                q = qsm[j]
                c = i0 + m
                lo_i, hi_i = max(0, c - kP // 2), min(n, c + kP // 2)
                plat = float(np.percentile(q[lo_i:hi_i], 75))
                half = (plat + q[c]) / 2
                a_, b_ = c, c
                while a_ > lo_i and q[a_] < half:
                    a_ -= 1
                while b_ < hi_i - 1 and q[b_] < half:
                    b_ += 1
                if plat > q[c] and a_ > lo_i and b_ < hi_i - 1:
                    widths.append((b_ - a_) * step)
                    gaps[ii, j] = (x[a_] + x[b_]) / 2
    strength = float(np.nanmedian(sf)) if np.isfinite(sf).any() else 0.0
    gw = float(np.clip(np.median(widths), 2 * step, 0.45 * P)) if widths else np.nan
    return dict(pitch=P, sign=sign, ph=ph, ks=ks, gaps=gaps, sf=sf, strength=strength,
                gw=gw, f=fs, ac=acs)


def _cells(cb, qs, step, size):
    """plots between neighbouring gaps of every profile: list of (k, j, centre,
    mean of the profile inside, strength of the weaker gap)"""
    P, gaps, sf, ks = cb["pitch"], cb["gaps"], cb["sf"], cb["ks"]
    n = len(qs[0])
    out = []
    for j, q in enumerate(qs):
        for ii in range(len(ks) - 1):
            a_, b_ = gaps[ii, j], gaps[ii + 1, j]
            sa_, sb_ = sf[ii, j], sf[ii + 1, j]
            if np.isfinite(a_) and not np.isfinite(b_):
                b_, sb_ = a_ + P, 0.0
            elif np.isfinite(b_) and not np.isfinite(a_):
                a_, sa_ = b_ - P, 0.0
            elif not np.isfinite(a_):
                continue
            if not (0.6 * P < b_ - a_ < 1.4 * P):
                continue
            c = (a_ + b_) / 2
            i0, i1 = max(0, int((c - size / 2) / step)), min(n, int((c + size / 2) / step))
            seg = q[i0:i1]
            if i1 - i0 < 3 or np.isfinite(seg).sum() < 0.7 * size / step:
                continue
            out.append((int(ks[ii]), j, c, float(np.nanmean(seg)),
                        float(min(np.nan_to_num(sa_), np.nan_to_num(sb_)))))
    return out


def _trim(present, n_want, label, log):
    """contiguous index block holding the plots; with a count, drop the ends with the
    fewest plots first (outer columns are often buffers)"""
    ks = sorted(present)
    k0, k1 = ks[0], ks[-1]
    if n_want and (k1 - k0 + 1) > n_want:
        while (k1 - k0 + 1) > n_want:
            if present.get(k0, 0) <= present.get(k1, 0):
                k0 += 1
            else:
                k1 -= 1
        log("  %s: %d found, kept %d as asked (the outer ones are taken as buffers - "
            "check them)" % (label, ks[-1] - ks[0] + 1, n_want))
    elif n_want and (k1 - k0 + 1) < n_want:
        log("  [!] %s: only %d found, %d asked for" % (label, k1 - k0 + 1, n_want))
    return k0, k1


def _col_profiles(G, runs_c, step, row_pitch=None):
    """along-row signal of every column, each divided by its local plateau
    (0 = bare, ~1 = crop). With visible rows the signal is the ROW AMPLITUDE: the
    greenness across the column correlated with the rows at the trial's spacing and
    phase, summed over 20 cm along. Only rows parallel to the plot's add up, so
    weeds, stubble and headland turns drilled at another angle (often in the alleys)
    drop out. Without rows it is the mean greenness (smoothed 10 cm)."""
    out = []
    for cs, ce in runs_c:
        w = 0.1 * (ce - cs)
        ic0, ic1 = max(0, int((cs + w) / step)), int(math.ceil((ce - w) / step))
        ic1 = max(ic1, ic0 + 1)
        B = G[:, ic0:ic1].astype(np.float64)
        valid = np.isfinite(B).sum(axis=1) > 0.2 * B.shape[1]
        if row_pitch:
            cpos = (np.arange(ic0, ic1) + 0.5) * step
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                B = np.nan_to_num(B - np.nanmean(B, axis=1, keepdims=True))
            Z = B @ np.exp(-2j * np.pi * cpos / row_pitch) / max(B.shape[1], 1)
            k = max(1, int(round(0.2 / step)))
            p = np.abs(np.convolve(Z, np.ones(k) / k, "same"))
        else:
            prof = _share(G[:, ic0:ic1], 1)
            p = _smooth(prof, 0.10 / step)
        p[~valid] = np.nan
        env = _envelope(p, 4.0 / step)
        ok = np.isfinite(p) & np.isfinite(env) & (env > 0)
        q = np.full(len(p), np.nan)
        q[ok] = np.clip(p[ok] / env[ok], 0, 1.5)
        out.append(q)
    return out


def _dark_profiles(GL, strips, step):
    """soil darkness along every column: (local bright level - brightness) / it,
    smoothed 10 cm and normalised like the other profiles (bare alleys ~0)"""
    out = []
    for cs, ce in strips:
        ic0, ic1 = max(0, int(cs / step)), int(math.ceil(ce / step))
        b = _share(GL[:, ic0:max(ic1, ic0 + 1)], 1)
        b = _smooth(b, 0.10 / step)
        hi = _envelope(b, 6.0 / step)
        d = np.where(np.isfinite(b) & np.isfinite(hi) & (hi > 0), (hi - b) / hi, np.nan)
        out.append(_normalise(d, step))
    return out


def _pitch_ac(qs, step, lo_m=1.5, hi_m=40.0):
    """pitch: peak of the summed autocorrelation of the profiles in [lo_m, hi_m],
    preferring the shortest strong harmonic. Returns (pitch, mean peak height)."""
    n = max(len(q) for q in qs)
    nf = 1 << int(math.ceil(math.log2(2 * n)))
    ac = np.zeros(nf // 2)
    for q in qs:
        x = np.nan_to_num(q - np.nanmean(q)) if np.isfinite(q).any() else np.zeros(len(q))
        f = np.fft.rfft(x, nf)
        a = np.fft.irfft(f * np.conj(f), nf)[: nf // 2]
        if a[0] > 0:
            ac += a / a[0]
    # unbiased: fewer overlapping samples at long lags
    lag = np.arange(len(ac))
    ac = ac / np.maximum(1 - lag / n, 0.2)
    i0, i1 = int(lo_m / step), min(len(ac) - 1, int(min(hi_m, 0.6 * n * step) / step))
    if i1 <= i0 + 2:
        return np.nan, 0.0
    seg = _smooth(ac[i0:i1], max(1, 0.1 / step))
    pk = [i for i in range(1, len(seg) - 1) if seg[i] >= seg[i - 1] and seg[i] > seg[i + 1]]
    if not pk:
        return np.nan, 0.0
    best = max(pk, key=lambda i: seg[i])
    top = seg[best]
    for i in sorted(pk):                     # shortest peak nearly as strong
        if seg[i] >= 0.75 * top:
            best = i
            break
    return (i0 + best) * step, float(seg[best] / max(len(qs), 1))


def _row_pitch(profile, step):
    from row_detect import row_pitch_from_profile
    p = np.nan_to_num(profile, nan=float(np.nanmean(profile)) if np.isfinite(profile).any() else 0)
    det = p - _smooth(p, 0.5 / step)
    pitch, share = row_pitch_from_profile(det, step, *ROW_BAND_M)
    return pitch, share


# ---------------------------------------------------------------------------
# 3 each plot
# ---------------------------------------------------------------------------
def _edge(p, x, target, search, level, outward):
    """half-height crossing of profile p (positions x) nearest to target within
    +-search. outward = -1 for a low edge (crop on the right), +1 for a high edge.
    Returns position or None."""
    m = np.isfinite(p)
    best = None
    for i in range(len(p) - 1):
        if not (m[i] and m[i + 1]):
            continue
        a, b = p[i] - level, p[i + 1] - level
        if outward < 0 and not (a < 0 <= b):
            continue
        if outward > 0 and not (a >= 0 > b):
            continue
        t = a / (a - b) if a != b else 0.5
        pos = x[i] + t * (x[i + 1] - x[i])
        if abs(pos - target) <= search and (best is None or abs(pos - target) < abs(best - target)):
            best = pos
    return best


def _measure(G, a_lo, a_hi, c_lo, c_hi, step, axis, smooth_m, plateau_band, target,
             search, gap):
    """sides (axis=1, across) or ends (axis=0, along) of one plot.
    G window rows = a, cols = c. Returns (lo, hi, contrast_lo, contrast_hi, cover)."""
    na, nc = G.shape
    ia0, ia1 = max(0, int(a_lo / step)), min(na, int(math.ceil(a_hi / step)))
    ic0, ic1 = max(0, int(c_lo / step)), min(nc, int(math.ceil(c_hi / step)))
    if ia1 - ia0 < 4 or ic1 - ic0 < 4:
        return None, None, 0.0, 0.0, np.nan
    sub = G[ia0:ia1, ic0:ic1]
    prof = _share(sub, 0 if axis == 1 else 1)
    x0 = (ic0 if axis == 1 else ia0) * step
    x = x0 + (np.arange(len(prof)) + 0.5) * step
    ps = _smooth(prof, smooth_m / step)
    ps[~np.isfinite(prof)] = np.nan
    lo_t, hi_t = target
    inside = (x > lo_t + plateau_band * (hi_t - lo_t)) & (x < hi_t - plateau_band * (hi_t - lo_t))
    if not np.isfinite(ps[inside]).any():
        return None, None, 0.0, 0.0, np.nan
    plateau = float(np.nanmedian(ps[inside]))
    cover = float(np.nanmean(prof[(x > lo_t) & (x < hi_t)]))
    if plateau <= 0.02:
        return None, None, 0.0, 0.0, cover
    level = 0.5 * plateau
    lo = _edge(ps, x, lo_t, search, level, -1)
    hi = _edge(ps, x, hi_t, search, level, +1)

    def contrast(e, side):
        if e is None:
            return 0.0
        w = (x > e - gap) & (x < e) if side < 0 else (x > e) & (x < e + gap)
        if not np.isfinite(ps[w]).any():
            return 0.0
        return float(np.clip((plateau - np.nanmin(ps[w])) / plateau, 0, 1))
    return lo, hi, contrast(lo, -1), contrast(hi, +1), cover


def _ends_1d(q, step, target, search, gap):
    """plot ends on a column profile (already smoothed and normalised): half-height
    crossings nearest the expected ends. Returns (lo, hi, contrast_lo, contrast_hi)."""
    n = len(q)
    x = (np.arange(n) + 0.5) * step
    lo_t, hi_t = target
    inside = (x > lo_t + 0.15 * (hi_t - lo_t)) & (x < hi_t - 0.15 * (hi_t - lo_t))
    if not np.isfinite(q[inside]).any():
        return None, None, 0.0, 0.0
    plateau = float(np.nanmedian(q[inside]))
    if plateau <= 0.05:
        return None, None, 0.0, 0.0
    level = 0.5 * plateau
    i0 = max(0, int((lo_t - search - gap) / step))
    i1 = min(n, int((hi_t + search + gap) / step) + 1)
    lo = _edge(q[i0:i1], x[i0:i1], lo_t, search, level, -1)
    hi = _edge(q[i0:i1], x[i0:i1], hi_t, search, level, +1)

    def contrast(e, side):
        if e is None:
            return 0.0
        w = (x > e - gap) & (x < e) if side < 0 else (x > e) & (x < e + gap)
        if not np.isfinite(q[w]).any():
            return 0.0
        return float(np.clip((plateau - np.nanmin(q[w])) / plateau, 0, 1))
    return lo, hi, contrast(lo, -1), contrast(hi, +1)


def _crop_width(pc, step, col_c, Pc):
    """median half-height width of the crop of every column on the trial-wide
    across profile (rows smoothed away); nan when unclear"""
    x = (np.arange(len(pc)) + 0.5) * step
    ws = []
    for c in col_c:
        m = (x > c - 0.3 * Pc) & (x < c + 0.3 * Pc)
        if not np.isfinite(pc[m]).any():
            continue
        lvl = 0.5 * float(np.nanmedian(pc[m]))
        lo = _edge(pc, x, c - 0.4 * Pc, 0.3 * Pc, lvl, -1)
        hi = _edge(pc, x, c + 0.4 * Pc, 0.3 * Pc, lvl, +1)
        if lo is not None and hi is not None and hi > lo:
            ws.append(hi - lo)
    return float(np.median(ws)) if len(ws) >= max(2, len(col_c) // 3) else np.nan


def _rows_between_gaps(sc, step, col_c, Pc, row_pitch):
    """crop rows per plot: the row peaks of the trial-wide across profile (every
    range adds up, so even seedlings show) inside each column, kept half a row
    spacing clear of the gap centres (a green strip between plots is not a row).
    Returns (rows, columns that agree) or (0, 0)."""
    p = _smooth(_fill_nan(sc), max(1, int(round(0.04 / step))))
    p[~np.isfinite(sc)] = np.nan
    x = (np.arange(len(p)) + 0.5) * step
    sep = 0.6 * row_pitch / step
    ns = []
    for c in col_c:
        m = np.nonzero((x > c - Pc / 2 + 0.5 * row_pitch) & (x < c + Pc / 2 - 0.5 * row_pitch))[0]
        if m.size < 5 or not np.isfinite(p[m]).all():
            continue
        v = p[m]
        lo, hi = float(np.percentile(v, 10)), float(v.max())
        if not hi > lo:
            continue
        pk = [i for i in range(1, len(v) - 1) if v[i] >= v[i - 1] and v[i] > v[i + 1]
              and v[i] > lo + 0.35 * (hi - lo)]
        pk.sort(key=lambda i: -v[i])
        kept = []
        for i in pk:
            if all(abs(i - k) >= sep for k in kept):
                kept.append(i)
        if kept:
            ns.append(len(kept))
    if len(ns) < max(2, len(col_c) // 3):
        return 0, 0
    vals, cnt = np.unique(ns, return_counts=True)
    return int(vals[np.argmax(cnt)]), int(cnt.max())


def _count_rows(ortho_path, F, O, ah, ch, a0, c0, cc, ac, col_c, EA, Pc, row_pitch,
                L, log, n_sample=9):
    """crop rows per plot from row_detect on a few plots spread over the trial; each
    window is one plot pitch less one row wide, centred on the plot, so it holds the
    plot's rows and none of its neighbours'. Returns (rows, spacing) or (0, nan)."""
    import crop_rows
    nr, ncol = cc.shape
    cand = [(i, j) for i in range(nr) for j in range(ncol)
            if np.isfinite(cc[i, j]) and np.isfinite(ac[i, j])]
    if not cand:
        return 0, np.nan
    pick = [cand[int(k)] for k in np.linspace(0, len(cand) - 1, min(n_sample, len(cand)))]
    half_w = max(Pc - row_pitch, 2 * row_pitch) / 2
    rings, ids = [], []
    for i, j in pick:
        a_, c_ = ac[i, j], cc[i, j]
        ring = []
        for d_a, d_c in ((-0.3 * L, -half_w), (-0.3 * L, half_w), (0.3 * L, half_w),
                         (0.3 * L, -half_w)):
            Pt = O + (a0 + a_ + d_a) * ah + (c0 + c_ + d_c) * ch
            ring.append((float(Pt[0] / F["mx"]), float(Pt[1] / F["my"])))
        rings.append(ring)
        ids.append("%d_%d" % (i, j))
    try:
        crs = F["crs"].to_wkt() if F["crs"] is not None else None
        _, rep = crop_rows.detect_rows_in_grid(rings, ids, ortho_path, grid_crs=crs,
                                               pitch_m=row_pitch, pad_frac=0.0,
                                               log=lambda *a: None)
    except Exception as e:
        log("  [!] counting the rows failed: %s" % e)
        return 0, np.nan
    ns = [r["n_rows"] for r in rep if r.get("ok") and r.get("n_rows")]
    ps = [float(r["pitch_m"]) for r in rep if r.get("ok") and r.get("pitch_m")]
    if len(ns) < max(2, len(pick) // 2):
        return 0, np.nan
    vals, cnt = np.unique(ns, return_counts=True)
    n = int(vals[np.argmax(cnt)])
    p = float(np.median(ps)) if ps else row_pitch
    log("  %d crop rows per plot every %.1f cm (%d of %d sample plots agree)"
        % (n, 100 * p, int(cnt.max()), len(ns)))
    return n, p


def _robust_line(t, y, w):
    """y ~ k t + d, iteratively reweighted (Huber-like); returns predict(t)"""
    t, y, w = map(lambda q: np.asarray(q, float), (t, y, w))
    ok = np.isfinite(y) & (w > 0)
    if ok.sum() == 0:
        return lambda q: np.full(np.shape(q), np.nan)
    if ok.sum() < 3:
        m = float(np.average(y[ok], weights=w[ok]))
        return lambda q: np.full(np.shape(q), m)
    ww = w.copy()
    k, d = 0.0, float(np.median(y[ok]))
    for _ in range(6):
        A = np.c_[t[ok], np.ones(ok.sum())] * np.sqrt(ww[ok])[:, None]
        sol = np.linalg.lstsq(A, y[ok] * np.sqrt(ww[ok]), rcond=None)[0]
        k, d = float(sol[0]), float(sol[1])
        r = np.abs(y - (k * t + d))
        s = max(1.4826 * float(np.median(r[ok])), 0.01)
        ww = w * np.where(r < 2 * s, 1.0, 2 * s / np.maximum(r, 1e-9))
    return lambda q: k * np.asarray(q, float) + d


# ---------------------------------------------------------------------------
# main entry
# ---------------------------------------------------------------------------
def detect_plots(ortho_path, quad=None, n_ranges=0, n_cols=0, plot_w=0.0, plot_l=0.0,
                 log=print, progress=None, cancel=None):
    """Find the plots of one trial.

    ortho_path  RGB orthomosaic (any CRS; lon/lat is handled in local metres)
    quad        optional 4 trial corners BL, BR, TR, TL in the ortho's CRS: limits the
                search to the trial and fixes the numbering (range 1 at BL-BR, plot 1
                of each range at BL-TL, as in 1a). Without it the whole orthomosaic is
                used, ranges are numbered from the south (or west) end and plots from
                the left looking along the rows.
    n_ranges, n_cols  optional counts (0 = find them). With more plots found than
                asked for, the middle block is kept (outer ones = buffers).
    plot_w, plot_l    optional plot size in metres (0 = median of the measured plots)

    Returns dict(plots=[{ring, props, conf, flags, ...}], quad=[BL, BR, TR, TL] world,
    info={...}) - rings in the orthomosaic's CRS, in range-major order.
    """
    t0 = time.time()

    def prog(f, msg):
        if progress:
            progress(f, msg)
        if cancel is not None and cancel():
            raise DetectError("cancelled")

    prog(0.02, "reading the orthomosaic")
    F = read_field(ortho_path, quad)
    log("Auto-detect: orthomosaic read at %.1f cm, %.1f %% green (ExG > %.3f)"
        % (100 * F["step"], 100 * F["cover"], F["exg_threshold"]))
    if F["cover"] < 0.02 and F["step"] > 0.012:
        F = read_field(ortho_path, quad, step_m=max(0.01, F["step"] / 2))
        log("  little green - re-read at %.1f cm, %.1f %% green" % (100 * F["step"], 100 * F["cover"]))
    mx, my = F["mx"], F["my"]

    # ---- 1 direction --------------------------------------------------------
    prog(0.1, "finding the row direction")
    around = None
    Q = None
    if quad:
        Q = [np.array([p[0] * mx, p[1] * my]) for p in quad]          # BL BR TR TL
        up = (Q[3] - Q[0]) + (Q[2] - Q[1])
        around = math.atan2(up[1], up[0])
    theta, clarity = find_direction(F, around, log)
    ah = np.array([math.cos(theta), math.sin(theta)])                # along the rows
    if Q is not None:
        if float(ah @ up) < 0:
            ah = -ah
    elif ah[1] < -1e-6 or (abs(ah[1]) <= 1e-6 and ah[0] < 0):
        ah = -ah                                                      # point north / east
    ch = np.array([ah[1], -ah[0]])                                    # across: a turned -90
    if Q is not None and float(ch @ ((Q[1] - Q[0]) + (Q[2] - Q[3]))) < 0:
        ch = -ch
    log("  rows run at %.2f deg from east (row signal %.1fx background)"
        % (math.degrees(math.atan2(ah[1], ah[0])), clarity))

    # ---- 2 trial frame ------------------------------------------------------
    if Q is not None:
        pts = np.array(Q)
        # clicks a little inside the outer plots still keep them; a buffer column
        # beyond the corners does not come in (it would be a sliver, then dropped)
        pad = min(0.6, 0.02 * max(np.ptp(pts @ ah), np.ptp(pts @ ch)))
    else:
        H, W = F["veg"].shape
        rr = np.array([0, 0, H, H]) - 0.5
        cc = np.array([0, W, 0, W]) - 0.5
        x, y = _pixel_xy_m(F, rr, cc)
        pts = np.c_[x, y]
        pad = 0.0
    O = np.zeros(2)
    a_rng = (float((pts @ ah).min() - pad), float((pts @ ah).max() + pad))
    c_rng = (float((pts @ ch).min() - pad), float((pts @ ch).max() + pad))
    step = F["step"]
    prog(0.2, "resampling into the trial frame")
    G = _resample(F, O, ah, ch, a_rng, c_rng, step)
    GL = _resample(F, O, ah, ch, a_rng, c_rng, step, "lum")
    if Q is not None:
        # outside the clicked quad = no data (other trials, headlands)
        qa = np.array([q @ ah for q in Q]) - a_rng[0]
        qc = np.array([q @ ch for q in Q]) - c_rng[0]
        A = (np.arange(G.shape[0]) + 0.5) * step
        C = (np.arange(G.shape[1]) + 0.5) * step
        AA, CC = np.meshgrid(A, C, indexing="ij")
        # point-in-convex-quad test (sign consistent with the quad's turn)
        s = None
        inside = np.ones(G.shape, bool)
        for k in range(4):
            p0 = (qa[k], qc[k])
            p1 = (qa[(k + 1) % 4], qc[(k + 1) % 4])
            cr = (p1[0] - p0[0]) * (CC - p0[1]) - (p1[1] - p0[1]) * (AA - p0[0])
            if s is None:
                ctr = (qa.mean(), qc.mean())
                s = np.sign((p1[0] - p0[0]) * (ctr[1] - p0[1]) - (p1[1] - p0[1]) * (ctr[0] - p0[0]))
            # grow the quad by the pad so clicks just inside the crop still keep it
            inside &= (cr * s) >= -pad * math.hypot(p1[0] - p0[0], p1[1] - p0[1])
        G[~inside] = np.nan
        GL[~inside] = np.nan
        del AA, CC

    # ---- 2a across the rows: the columns ---------------------------------------
    # a column of plots is one seeder pass, straight through all the ranges, so one
    # profile over the whole trial (rows smoothed away) holds every gap between plots
    prog(0.3, "finding the plots across the rows")
    sc = _share(G, 0)
    row_pitch, row_share = _row_pitch(sc, step)
    rows_seen = bool(np.isfinite(row_pitch))
    if not rows_seen:
        row_pitch = 0.25
    pc = _smooth(sc, row_pitch / step)
    pc[~np.isfinite(sc)] = np.nan
    qc = _normalise(pc, step)
    cbc = _comb([qc], step, max(0.6, 3 * row_pitch), 8.0)
    continuous = cbc is None or cbc["ac"] < 0.12
    if not continuous:
        Pc = cbc["pitch"]
        W0 = Pc - cbc["gw"] if np.isfinite(cbc["gw"]) else 0.75 * Pc
        cells = _cells(cbc, [qc], step, 0.6 * Pc)
        lvl = np.median([c[3] for c in cells]) if cells else 0
        cells = [c for c in cells if c[3] > 0.4 * lvl]
        if len(cells) < 2:
            continuous = True
    if continuous:
        ok = np.nonzero(np.nan_to_num(qc) > 0.35)[0]
        if not ok.size:
            raise DetectError("no crop found across the rows")
        lo, hi = ok[0] * step, (ok[-1] + 1) * step
        if n_cols:
            n_c = n_cols
        elif plot_w:
            n_c = max(1, int(round((hi - lo) / plot_w)))
        else:
            raise DetectError(
                "no gaps between the plots across the rows (the trial looks drilled "
                "continuously). Enter 'Plots per range' (or the plot width), tick 'Use "
                "my counts and sizes' and run Auto-detect again.")
        Pc = (hi - lo) / n_c
        col_c = lo + (np.arange(n_c) + 0.5) * Pc
        col_s = np.full(n_c, np.nan)
        W0 = plot_w or Pc
        sign_c = 0.0
        log("  no gaps between plots across the rows: %d plots laid out evenly" % n_c)
    else:
        present = {}
        for k, _, c, _, _ in cells:
            present[k] = present.get(k, 0) + 1
        k0, k1 = _trim(present, n_cols, "plots per range", log)
        byk = dict((k, (c, st)) for k, _, c, _, st in cells)
        col_c = np.array([byk[k][0] if k in byk else np.nan for k in range(k0, k1 + 1)])
        # an empty column inside the trial: at the pitch from its neighbours
        idx = np.arange(len(col_c))
        okc = np.isfinite(col_c)
        col_c = np.interp(idx, idx[okc], col_c[okc] - idx[okc] * Pc) + idx * Pc
        sign_c = cbc["sign"]
        log("  plots every %.2f m across the rows (gaps show as %s)"
            % (Pc, "bare strips" if sign_c < 0 else "green strips"))
    ncol = len(col_c)

    # ---- 2b along the rows: the ranges, column by column -------------------------
    # The alleys are often staggered (each pass lifts the seeder at a slightly
    # different place), so every column has its own comb phase. With rows visible
    # the along signal is the row amplitude: only rows parallel to the plot's count,
    # so stubble, weeds and headland turns drilled at another angle drop out.
    prog(0.4, "finding the ranges in each column")
    strips = [(c - 0.35 * W0, c + 0.35 * W0) for c in col_c]
    use_rows = rows_seen and not continuous
    # Three candidate signals; the one whose plot boxes stand out most wins:
    #   row amplitude  rows parallel to the plot (stubble, weeds, headland turns out)
    #   greenness      crop vs bare alley
    #   soil darkness  bare alleys are lighter, wheel-tracked soil: at emergence the
    #                  drill rows often run straight through them, so only this shows
    cands = []
    if use_rows:
        cands.append(("row amplitude", _col_profiles(G, strips, step, row_pitch)))
    cands.append(("greenness", _col_profiles(G, strips, step, None)))
    cands.append(("soil darkness", _dark_profiles(GL, strips, step)))
    cba, qs, sig_a = None, None, ""
    for name, q_ in cands:
        r = _boxcar_comb(q_, step, 1.5, 40.0)
        if r is not None and (cba is None or r["score"] > cba["score"]):
            cba, qs, sig_a = r, q_, name
    if qs is None:
        qs = cands[0][1]
    if cba is None:
        # a single range: the crop run of every column
        La = []
        acen = np.full((1, ncol), np.nan)
        for j, q in enumerate(qs):
            r = _runs(np.nan_to_num(q) > 0.35)
            if r:
                s_, e_ = max(r, key=lambda t: t[1] - t[0])
                acen[0, j] = (s_ + e_) / 2 * step
                La.append((e_ - s_) * step)
        if not La:
            raise DetectError("no crop found along the rows")
        Pa, L0, nr = float(np.median(La)) + 0.5, float(np.median(La)), 1
        EA = np.where(np.isfinite(acen), acen, np.nanmedian(acen))
        a_con = np.full((1, ncol), 0.6)
    else:
        Pa, L0 = cba["pitch"], cba["size"]
        cells_a = cba["cells"]
        lvl = np.median([c[4] for c in cells_a]) if cells_a else 0
        cells_a = [c for c in cells_a if c[4] > 0.4 * lvl]
        if not cells_a:
            raise DetectError("no ranges found along the rows")
        present = {}
        for c in cells_a:
            present[c[0]] = present.get(c[0], 0) + 1
        k0, k1 = _trim(present, n_ranges, "ranges", log)
        nr = k1 - k0 + 1
        acen = np.full((nr, ncol), np.nan)
        a_con = np.zeros((nr, ncol))
        for k, j, s_, e_, _, con in cells_a:
            if k0 <= k <= k1:
                acen[k - k0, j] = (s_ + e_) / 2
                a_con[k - k0, j] = con
        EA = cba["ph"][None, :] + np.arange(k0, k1 + 1)[:, None] * Pa + L0 / 2
        log("  ranges every %.2f m, plots %.2f m long (from %s)" % (Pa, L0, sig_a))
    log("  %d ranges x %d plots per range" % (nr, ncol))

    # ---- 3 each plot ---------------------------------------------------------------
    prog(0.5, "measuring the plots")
    # along: the box of every plot (both crop ends at once)
    L = float(plot_l) if plot_l else L0
    ac = acen.copy()
    s_a = np.clip(np.nan_to_num(a_con) / 0.4, 0, 1)

    # across: every plot's own gaps, in the stretch of its range
    cc = np.full((nr, ncol), np.nan)
    s_c = np.full((nr, ncol), 0.6 if continuous else 0.0)
    if not continuous:
        g = min(max(0.12 * Pc, 3 * step), 0.6)
        kg, kP = max(1, int(round(g / step))), max(3, int(round(Pc / step)))
        x = (np.arange(G.shape[1]) + 0.5) * step
        for i in range(nr):
            prog(0.5 + 0.35 * i / max(nr, 1), "measuring range %d of %d" % (i + 1, nr))
            am = float(np.nanmedian(np.where(np.isfinite(ac[i]), ac[i], EA[i])))
            ia0 = max(0, int((am - 0.35 * L) / step))
            ia1 = min(G.shape[0], int((am + 0.35 * L) / step))
            if ia1 - ia0 < 4:
                continue
            p = _share(G[ia0:ia1], 0)
            p = _smooth(p, row_pitch / step)
            q = _normalise(p, step)
            v = _fill_nan(q)
            f = sign_c * (_smooth(v, kg) - _smooth(v, kP))
            f[~np.isfinite(q)] = np.nan
            ref = abs(cbc["strength"]) or 1e-9
            for j in range(ncol):
                feats = []
                for side in (-1, 1):
                    t = col_c[j] + side * Pc / 2
                    i0 = max(0, int((t - 0.1 * Pc) / step))
                    i1 = min(len(f), int((t + 0.1 * Pc) / step))
                    if i1 - i0 < 3 or not np.isfinite(f[i0:i1]).any():
                        feats.append((None, 0.0))
                        continue
                    m = int(np.nanargmax(f[i0:i1]))
                    feats.append((x[i0 + m], float(f[i0 + m]) / ref))
                (l_, sl), (r_, sr) = feats
                if l_ is not None and r_ is not None and min(sl, sr) > 0.3:
                    cc[i, j] = (l_ + r_) / 2
                elif l_ is not None and sl > 0.3:
                    cc[i, j] = l_ + Pc / 2
                elif r_ is not None and sr > 0.3:
                    cc[i, j] = r_ - Pc / 2
                s_c[i, j] = float(np.clip(min(sl, sr) / 0.5, 0, 1))
    else:
        cc[:] = col_c[None, :]

    # width: rows x spacing when the rows are seen (the 1a/1c plot definition)
    n_rows_plot = 0
    if plot_w:
        W = float(plot_w)
    else:
        W = W0
        We = _crop_width(pc, step, col_c, Pc) if sign_c < 0 else np.nan
        if use_rows:
            # width = rows x spacing: the plot edge half a spacing outside the outer
            # rows, as in 1a / 1c
            n_rows_plot, agree = _rows_between_gaps(sc, step, col_c, Pc, row_pitch)
            if n_rows_plot:
                W = n_rows_plot * row_pitch
                log("  %d crop rows per plot every %.1f cm (%d of %d columns agree)"
                    % (n_rows_plot, 100 * row_pitch, agree, ncol))
        if not n_rows_plot and np.isfinite(We):
            W = We
        if use_rows and not n_rows_plot:
            n_rows_plot, rp2 = _count_rows(ortho_path, F, O, ah, ch, a_rng[0], c_rng[0],
                                           cc, ac, col_c, EA, Pc, row_pitch, L, log)
            if n_rows_plot:
                W = n_rows_plot * rp2
                row_pitch = rp2
    prog(0.88, "checking the plots against their neighbours")

    # neighbour lines: a column is one seeder pass (straight through the ranges); the
    # plot ends of a range lie on a near-straight line across the columns
    gap_c = max(Pc - W, 0.0)
    gap_a = max(Pa - L, 0.0)
    tol_c = max(0.05, 0.35 * gap_c, 0.3 * row_pitch)
    tol_a = max(0.10, 0.35 * gap_a)
    cc_pred = np.full((nr, ncol), np.nan)
    ac_pred = np.full((nr, ncol), np.nan)
    # across, the neighbours are the other plots of the same range: its offset plus a
    # gentle tilt (robust line across the columns through the clear plots)
    jj = np.arange(ncol, dtype=float)
    for i in range(nr):
        base = col_c
        strong = np.isfinite(cc[i]) & (s_c[i] >= 0.6)
        if strong.sum() >= 4:
            cc_pred[i] = base + _robust_line(jj, cc[i] - base, np.clip(s_c[i], 0.1, 1))(jj)
        elif strong.any():
            cc_pred[i] = base + float(np.median((cc[i] - base)[strong]))
        else:
            cc_pred[i] = base
    for i in range(nr):
        c_t = np.where(np.isfinite(cc[i]), cc[i], col_c)
        ac_pred[i, :] = _robust_line(c_t, ac[i], np.clip(s_a[i], 0.1, 1))(c_t)
    cc_pred = np.where(np.isfinite(cc_pred), cc_pred, col_c[None, :])
    ac_pred = np.where(np.isfinite(ac_pred), ac_pred, EA)

    # crop cover of every plot (empty / failed plots)
    cover = np.full((nr, ncol), np.nan)
    for i in range(nr):
        for j in range(ncol):
            a_ = ac[i, j] if np.isfinite(ac[i, j]) else ac_pred[i, j]
            c_ = cc[i, j] if np.isfinite(cc[i, j]) else cc_pred[i, j]
            ia0, ia1 = max(0, int((a_ - 0.4 * L) / step)), int((a_ + 0.4 * L) / step)
            ic0, ic1 = max(0, int((c_ - 0.4 * W) / step)), int((c_ + 0.4 * W) / step)
            blk = G[ia0:ia1, ic0:ic1]
            if blk.size and np.isfinite(blk).mean() > 0.5:
                cover[i, j] = float(np.nanmean(blk))
    med_cover = float(np.nanmedian(cover)) if np.isfinite(cover).any() else 0.0
    # Without corners the layout spreads over everything in the orthomosaic: peel
    # off outer ranges / columns that are mostly empty (an empty plot INSIDE the
    # trial is kept and flagged).
    i_lo, i_hi, j_lo, j_hi = 0, nr - 1, 0, ncol - 1
    if quad is None and med_cover > 0:
        empty = ~(np.nan_to_num(cover) > 0.5 * med_cover)
        changed = True
        while changed and i_hi > i_lo and j_hi > j_lo:
            changed = False
            for side in ("i_lo", "i_hi", "j_lo", "j_hi"):
                blk = {"i_lo": empty[i_lo, j_lo:j_hi + 1], "i_hi": empty[i_hi, j_lo:j_hi + 1],
                       "j_lo": empty[i_lo:i_hi + 1, j_lo], "j_hi": empty[i_lo:i_hi + 1, j_hi]}[side]
                if blk.mean() > 0.6:
                    if side == "i_lo":
                        i_lo += 1
                    elif side == "i_hi":
                        i_hi -= 1
                    elif side == "j_lo":
                        j_lo += 1
                    else:
                        j_hi -= 1
                    changed = True
                    break
        if (i_lo, i_hi, j_lo, j_hi) != (0, nr - 1, 0, ncol - 1):
            log("  outer ranges / columns without crop left out: %d x %d plots kept"
                % (i_hi - i_lo + 1, j_hi - j_lo + 1))

    plots = []
    a0, c0 = a_rng[0], c_rng[0]
    for i in range(i_lo, i_hi + 1):
        for j in range(j_lo, j_hi + 1):
            flags = []
            dc = cc[i, j] - cc_pred[i, j]
            da = ac[i, j] - ac_pred[i, j]
            cov = cover[i, j]
            sv = float(np.clip((cov / med_cover) / 0.5, 0, 1)) if med_cover > 0 and np.isfinite(cov) else 0.0
            slc = float(np.clip(1 - abs(dc) / (2 * tol_c), 0, 1)) if np.isfinite(dc) else 0.4
            sla = float(np.clip(1 - abs(da) / (2 * tol_a), 0, 1)) if np.isfinite(da) else 0.4
            sc_, sa_ = float(s_c[i, j]), float(s_a[i, j])
            # use the measurement when it is clear or in line, else the neighbours' line
            c_use = cc[i, j] if (np.isfinite(dc) and (abs(dc) <= tol_c or
                                                      (sc_ >= 0.8 and abs(dc) <= 2 * tol_c))) else cc_pred[i, j]
            a_use = ac[i, j] if (np.isfinite(da) and (sa_ >= 0.5 or abs(da) <= tol_a)) else ac_pred[i, j]
            if sv < 0.5:
                flags.append("little crop (empty or failed plot?)")
            if not continuous and sc_ < 0.5:
                flags.append("side unclear" + ("" if np.isfinite(dc) else " - placed in line"))
            if sa_ < 0.5:
                flags.append("end unclear" + ("" if np.isfinite(da) else " - placed in line"))
            if np.isfinite(dc) and abs(dc) > tol_c:
                flags.append("%.0f cm out of line across" % (100 * abs(dc)))
            if np.isfinite(da) and abs(da) > tol_a:
                flags.append("%.0f cm out of line along" % (100 * abs(da)))
            conf = min(sc_, sa_, sv, slc, sla)
            ring = []
            for d_a, d_c in ((-L / 2, -W / 2), (-L / 2, W / 2), (L / 2, W / 2), (L / 2, -W / 2)):
                Pt = O + (a0 + a_use + d_a) * ah + (c0 + c_use + d_c) * ch
                ring.append((float(Pt[0] / mx), float(Pt[1] / my)))
            plots.append(dict(
                ring=ring,
                props={'Plot_ID': (i - i_lo + 1) * 1000 + (j - j_lo + 1),
                       'B/R': "Rg%dRw%d" % (i - i_lo + 1, j - j_lo + 1),
                       'Range': i - i_lo + 1, 'Row': j - j_lo + 1, 'Bank': i - i_lo + 1},
                conf=round(float(conf), 3), flags=flags,
                cover=round(float(cov), 4) if np.isfinite(cov) else None,
                shift_across_m=round(float(dc), 3) if np.isfinite(dc) else None,
                shift_along_m=round(float(da), 3) if np.isfinite(da) else None,
                centre_tf=(float(a_use), float(c_use))))
    nr, ncol = i_hi - i_lo + 1, j_hi - j_lo + 1
    ac = np.array([p["centre_tf"][0] for p in plots]).reshape(nr, ncol)
    cc = np.array([p["centre_tf"][1] for p in plots]).reshape(nr, ncol)
    # the trial quad (BL, BR, TR, TL) around the plots, world CRS
    a_lo_t = a0 + float(np.nanmin(ac)) - L / 2
    a_hi_t = a0 + float(np.nanmax(ac)) + L / 2
    c_lo_t = c0 + float(np.nanmin(cc)) - W / 2
    c_hi_t = c0 + float(np.nanmax(cc)) + W / 2
    quad_out = []
    for a_, c_ in ((a_lo_t, c_lo_t), (a_lo_t, c_hi_t), (a_hi_t, c_hi_t), (a_hi_t, c_lo_t)):
        P = O + a_ * ah + c_ * ch
        quad_out.append((float(P[0] / mx), float(P[1] / my)))
    n_flag = sum(1 for p in plots if p["conf"] < FLAG_BELOW)
    info = dict(n_ranges=nr, n_cols=ncol, plot_w=round(W, 3), plot_l=round(L, 3),
                gap_across=round(gap_c, 3), gap_along=round(gap_a, 3),
                range_pitch=round(Pa, 3), plot_pitch=round(Pc, 3),
                row_pitch=round(float(row_pitch), 4) if rows_seen else 0.0,
                crop_rows=int(n_rows_plot),
                angle_deg=round(math.degrees(math.atan2(ah[1], ah[0])), 3),
                direction_clarity=round(clarity, 2), continuous=continuous,
                step_m=round(step, 4), green_share=round(F["cover"], 4),
                n_flagged=n_flag, seconds=round(time.time() - t0, 1))
    prog(1.0, "done")
    log("  %d plots, %.3f m wide x %.3f m long; %d to check (confidence < %.1f); %.1f s"
        % (len(plots), W, L, n_flag, FLAG_BELOW, info["seconds"]))
    return dict(plots=plots, quad=quad_out, info=info, crs=F["crs"])


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def report_path_for(grid_path):
    import os
    return os.path.splitext(grid_path)[0] + "_detect.csv"


def write_report(path, result):
    """per-plot confidence CSV next to the grid (read by 1b 'Next to check')"""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Plot_ID", "Range", "Row", "confidence", "check", "reasons", "crop_cover",
                    "shift_across_m", "shift_along_m"])
        for p in result["plots"]:
            pr = p["props"]
            w.writerow([pr["Plot_ID"], pr["Range"], pr["Row"], p["conf"],
                        "yes" if p["conf"] < FLAG_BELOW else "",
                        "; ".join(p["flags"]), p["cover"], p["shift_across_m"],
                        p["shift_along_m"]])
    return path


def read_report(path):
    """{Plot_ID (str): (confidence, reasons)} from a _detect.csv, or {}"""
    import os
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            try:
                out[str(r["Plot_ID"])] = (float(r["confidence"]), r.get("reasons", ""))
            except (KeyError, ValueError):
                pass
    return out


if __name__ == "__main__":
    import argparse
    import json
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ortho")
    ap.add_argument("--quad", help="BLx,BLy;BRx,BRy;TRx,TRy;TLx,TLy in the ortho CRS")
    ap.add_argument("--ranges", type=int, default=0)
    ap.add_argument("--cols", type=int, default=0)
    ap.add_argument("--out", help="grid .shp to write")
    a = ap.parse_args()
    q = [tuple(map(float, s.split(","))) for s in a.quad.split(";")] if a.quad else None
    res = detect_plots(a.ortho, q, a.ranges, a.cols)
    print(json.dumps(res["info"], indent=1))
    if a.out:
        from grid_designer import write_grid_shapefile
        write_grid_shapefile(a.out, res["plots"], res["crs"])
        write_report(report_path_for(a.out), res)
