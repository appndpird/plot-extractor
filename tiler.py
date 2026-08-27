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


def _plot_mm_px(src):
    """Ground sampling of an open raster, mm per pixel, from its GCPs.

    Returns None when the file carries no usable georeferencing.
    """
    try:
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


def make_tiles(src_dir, out_dir, tile=1024, overlap=0, rows_only=True,
               row_pitch_mm=0.0, plot_width_mm=0.0, fmt="png",
               min_veg_pct=0.0, keep_empty_frac=0.15, max_tiles_per_plot=0,
               seed=7, target_mm_px=0.0, per_plot_dirs=False,
               log=print, progress=None, cancel=None):
    """Cut every crop in `src_dir` into tiles.

    tile / overlap        pixels. stride = tile - overlap.
    rows_only             restrict to the measured seed-row band.
    row_pitch_mm          row spacing; only used to find the band. 0 = guess from
                          plot_width_mm and a nominal 5 rows.
    min_veg_pct           drop tiles below this vegetation fraction, EXCEPT a
                          random `keep_empty_frac` of them, kept as hard
                          negatives - bare soil and stubble is what a seedling
                          detector gets wrong, so it has to see some.
    max_tiles_per_plot    0 = all. Otherwise a random sample per plot.
    per_plot_dirs         write each plot's tiles into its own subfolder
                          out_dir/<plot stem>/ instead of one flat folder.

    Returns (n_tiles, n_plots, index_path).
    """
    import rasterio
    from PIL import Image

    ext = {"png": ".png", "jpeg": ".jpg", "jpg": ".jpg", "tif": ".tif"}.get(fmt.lower(), ".png")
    os.makedirs(out_dir, exist_ok=True)
    files = sorted(f for f in os.listdir(src_dir)
                   if f.lower().endswith((".tif", ".tiff")))
    if not files:
        log(f"[!] no GeoTIFF crops in {src_dir}")
        return 0, 0, ""
    stride = max(1, tile - max(0, overlap))
    rnd = random.Random(seed)
    rows_out = []
    log(f"Tiling {len(files)} plot(s): {tile}px tiles, {overlap}px overlap "
        f"(stride {stride}), {'row band only' if rows_only else 'whole crop'}, "
        f"{ext[1:].upper()} output.")

    for fi, fn in enumerate(files):
        if cancel is not None and cancel.is_set():
            log("Tiling cancelled.")
            break
        stem = os.path.splitext(fn)[0]
        path = os.path.join(src_dir, fn)
        try:
            with rasterio.open(path) as s:
                W, H = s.width, s.height
                # rows run along the long axis; 'across' is the short one
                across_axis = 0 if H < W else 1
                # Ground-referenced tiling. Ground sampling varies plot to plot
                # (0.274-0.470 mm/px on OZ Barley), so a fixed pixel tile covers a
                # different amount of ground on every plot - 2.2 rows on one, 3.8 on
                # another. Reading a window sized in GROUND units and resampling it
                # to `tile` px makes every tile the same real size, which is what a
                # "cover ~3 rows" requirement actually asks for.
                src_mm = _plot_mm_px(s) if target_mm_px > 0 else None
                win_px = tile
                if src_mm and src_mm > 0:
                    win_px = max(8, int(round(tile * target_mm_px / src_mm)))
                if row_pitch_mm > 0 and plot_width_mm > 0:
                    across_px = H if across_axis == 0 else W
                    pitch_px = row_pitch_mm / (plot_width_mm / across_px)
                else:
                    pitch_px = (H if across_axis == 0 else W) / 5.5
                lo, hi = 0, (H if across_axis == 0 else W)
                if rows_only:
                    dec = max(1, int(round(min(W, H) / 900)))
                    small = s.read(out_shape=(3, max(1, H // dec), max(1, W // dec)))
                    l2, h2 = row_band(small, across_axis, pitch_px / dec)
                    lo, hi = l2 * dec, min(hi, h2 * dec)
                    if hi - lo < win_px:    # band thinner than a tile: centre it
                        c = (lo + hi) // 2
                        lo = max(0, c - win_px // 2)
                        hi = min((H if across_axis == 0 else W), lo + win_px)

                y0b, y1b = (lo, hi) if across_axis == 0 else (0, H)
                x0b, x1b = (0, W) if across_axis == 0 else (lo, hi)
                # grid is laid out in SOURCE pixels; win_px == tile unless the
                # tiles are ground-referenced, in which case it is the source
                # window that resamples down to `tile`.
                src_stride = max(1, int(round(stride * win_px / tile)))
                cand = []
                for r_i, y in enumerate(_starts(y0b, y1b, win_px, src_stride)):
                    for c_i, x in enumerate(_starts(x0b, x1b, win_px, src_stride)):
                        cand.append((c_i, r_i, x, y))
                if not cand:
                    cand = [(0, 0, x0b, y0b)]
                if max_tiles_per_plot and len(cand) > max_tiles_per_plot:
                    cand = rnd.sample(cand, max_tiles_per_plot)

                kept = 0
                for (c_i, r_i, x, y) in cand:
                    w = min(win_px, W - x); h = min(win_px, H - y)
                    if w < win_px // 2 or h < win_px // 2:
                        continue
                    a = s.read(window=rasterio.windows.Window(x, y, w, h))
                    if a.shape[0] < 3:
                        continue
                    rgb = a[:3]
                    valid = rgb.max(axis=0) > 4
                    if valid.mean() < 0.5:            # mostly nodata
                        continue
                    R, G, B = (rgb[i].astype(np.float32) for i in range(3))
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
                    im = Image.fromarray(np.transpose(rgb, (1, 2, 0)).astype(np.uint8))
                    if win_px != tile:
                        # proportional, so a clipped window keeps the same scale
                        im = im.resize((max(1, int(round(w * tile / win_px))),
                                        max(1, int(round(h * tile / win_px)))),
                                       Image.LANCZOS)
                    if ext == ".jpg":
                        im.save(dest, quality=95, subsampling=0)
                    else:
                        im.save(dest)
                    rows_out.append(dict(tile=name, rel_path=rel, plot_stem=stem, col=c_i, row=r_i,
                                         x0=x, y0=y, w=w, h=h,
                                         out_px=tile,
                                         src_mm_px=round(src_mm, 4) if src_mm else "",
                                         tile_mm_px=round(target_mm_px, 4) if (src_mm and target_mm_px) else "",
                                         ground_mm=round(w * (src_mm or 0), 1) if src_mm else "",
                                         band_lo=lo, band_hi=hi,
                                         across_axis=across_axis,
                                         src_w=W, src_h=H, veg_pct=round(veg, 2),
                                         src=fn))
                    kept += 1
        except Exception as e:
            log(f"[!] {fn}: {e}")
        if progress:
            progress(fi + 1, len(files))
        if (fi + 1) % 25 == 0:
            log(f"   {fi+1}/{len(files)} plots, {len(rows_out)} tiles")

    idx = os.path.join(out_dir, "tile_index.csv")
    keys = ["tile", "rel_path", "plot_stem", "col", "row", "x0", "y0", "w", "h", "out_px",
            "src_mm_px", "tile_mm_px", "ground_mm",
            "band_lo", "band_hi", "across_axis", "src_w", "src_h", "veg_pct", "src"]
    with open(idx, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        for rr in rows_out:
            wr.writerow(rr)
    plots = len({r["plot_stem"] for r in rows_out})
    log(f"Wrote {len(rows_out)} tiles from {plots} plot(s) -> {out_dir}")
    log(f"Index: {idx}  (needed to map a tile back to its plot)")
    if rows_out:
        v = np.array([r["veg_pct"] for r in rows_out])
        log(f"Vegetation per tile: median {np.median(v):.2f}%  "
            f"empty (<0.5%) {100*(v<0.5).mean():.1f}%")
    return len(rows_out), plots, idx


# ---------------------------------------------------------------------------
# stitching back
# ---------------------------------------------------------------------------
def build_tiles_tab(parent, defaults=None):
    """The 'Training tiles' tab: cut tiles, and put labelled tiles back together.

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
        "tile": tk.StringVar(value=str(d.get("tile", 1024))),
        "overlap": tk.StringVar(value=str(d.get("overlap", 256))),
        "rows_only": tk.BooleanVar(value=bool(d.get("rows_only", True))),
        "pitch": tk.StringVar(value=str(d.get("pitch", ""))),
        "plotw": tk.StringVar(value=str(d.get("plotw", ""))),
        "minveg": tk.StringVar(value=str(d.get("minveg", 0.5))),
        "mmpx": tk.StringVar(value=str(d.get("mmpx", ""))),
        "subdirs": tk.BooleanVar(value=bool(d.get("subdirs", False))),
        "keepempty": tk.StringVar(value=str(d.get("keepempty", 15))),
        "maxtiles": tk.StringVar(value=str(d.get("maxtiles", 0))),
        "fmt": tk.StringVar(value=d.get("fmt", "png")),
        "seed": tk.StringVar(value=str(d.get("seed", 7))),
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

    ttk.Label(frm, text="Cut plot crops into tiles",
              font=("", 10, "bold")).grid(row=r, column=0, columnspan=3, sticky="w")
    r += 1
    path_row("Plot crops folder *", V["src"], "dir",
             "the RawCrops folder of a rectified extraction")
    path_row("Tiles output folder *", V["out"], "dir",
             "tiles are written flat here, plus tile_index.csv")

    o1 = ttk.Frame(frm)
    o1.grid(row=r, column=0, columnspan=3, sticky="w", pady=4)
    r += 1
    for text, var, w in (("Tile size (px):", V["tile"], 6),
                         ("Overlap (px):", V["overlap"], 6),
                         ("Max tiles/plot (0=all):", V["maxtiles"], 5),
                         ("Seed:", V["seed"], 4)):
        ttk.Label(o1, text=text).pack(side="left", padx=(0, 4))
        ttk.Entry(o1, textvariable=var, width=w).pack(side="left", padx=(0, 12))
    ttk.Label(o1, text="Format:").pack(side="left", padx=(0, 4))
    ttk.Combobox(o1, textvariable=V["fmt"], width=6, state="readonly",
                 values=["png", "tif", "jpeg"]).pack(side="left")

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

    o2b = ttk.Frame(frm)
    o2b.grid(row=r, column=0, columnspan=3, sticky="w", pady=4)
    r += 1
    ttk.Label(o2b, text="Fix tile ground scale (mm/px):").pack(side="left", padx=(0, 4))
    ttk.Entry(o2b, textvariable=V["mmpx"], width=7).pack(side="left")
    ttk.Label(frm, text="blank = tile is N source pixels, so it covers different "
                        "ground on each plot. Set it (e.g. 0.40) and every tile "
                        "covers the same real area: tile px x mm/px.",
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
                    variable=V["subdirs"]).pack(side="left")

    btn1 = ttk.Frame(frm)
    btn1.grid(row=r, column=0, columnspan=3, sticky="w", pady=6)
    r += 1
    make_btn = ttk.Button(btn1, text="Make tiles")
    make_btn.pack(side="left")
    cancel_btn = ttk.Button(btn1, text="Cancel", state="disabled")
    cancel_btn.pack(side="left", padx=8)
    prog = ttk.Progressbar(btn1, mode="determinate", length=260)
    prog.pack(side="left", padx=8)

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
        t = int(_f(V["tile"], 1024)); ov = int(_f(V["overlap"], 0))
        if ov >= t:
            messagebox.showwarning("Make tiles",
                                   "Overlap must be smaller than the tile size.")
            return
        start(lambda: make_tiles(
            src, out, tile=t, overlap=ov,
            rows_only=bool(V["rows_only"].get()),
            row_pitch_mm=_f(V["pitch"], 0.0), plot_width_mm=_f(V["plotw"], 0.0),
            fmt=V["fmt"].get(), min_veg_pct=_f(V["minveg"], 0.0),
            keep_empty_frac=max(0.0, min(1.0, _f(V["keepempty"], 15.0) / 100.0)),
            max_tiles_per_plot=int(_f(V["maxtiles"], 0)),
            seed=int(_f(V["seed"], 7)), target_mm_px=_f(V["mmpx"], 0.0),
            per_plot_dirs=bool(V["subdirs"].get()), log=put,
            progress=lambda a, b: q.put(("prog", (a, b))), cancel=cancel_ev))

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
                with Image.open(src) as im:
                    tim = im.convert("RGB")
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
