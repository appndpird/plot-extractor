"""row_review.py - tab 1c "Review crop rows": make, review and edit the crop-row lines
inside every plot of the grid (after the PhenoFuse 4.0.6/4.0.7 crop-row editor).

  Make rows      Fixed spacing: a known number of rows per plot, evenly spaced (plot width /
                 rows, or a given spacing), the same lines tab 1a writes at generation.
                 Detect: the rows found in the RGB orthomosaic, plot by plot, with the
                 PhenoFuse row method (row_detect.py) - the actual sowing, tilt included.
  Review         plots list with the row count of every plot (plots whose count differs from
                 the usual count are marked), Prev / Next plot, zoom to a plot. Zoomed in, the
                 view reads the orthomosaic at screen resolution, so single rows are sharp.
  Edit           click a row line to select it (Shift adds), drag it across its plot; arrow
                 keys nudge 1 cm (Shift 5 cm); Add row (A, then click in a plot); Delete (Del);
                 Even spacing (E) re-spaces the selected rows' plot between its outer rows;
                 Undo (Ctrl+Z).
  Save           <grid>_rows.shp (+ .csv: rows per plot). The first save over an existing
                 file keeps it as <rows>_before_edit.shp.
Colours: cyan = fixed spacing, yellow = detected, orange = moved / added, magenta = selected.
"""
from __future__ import annotations

import math
import os
import shutil
import threading

import numpy as np

import crop_rows as cr

COL = {"count": "#33ddff", "detect": "#ffd23f", "edit": "#ff9933", "sel": "#ff44ff"}
DETAIL_MAX_M = 40.0          # visible width (m) below which the orthomosaic is read at screen res


def _read_ortho(path, grid_crs, bounds, out_w, out_h, geographic):
    """RGB uint8 (out_h, out_w, 3) of the orthomosaic inside bounds (grid CRS)."""
    import rasterio
    from rasterio.crs import CRS
    from rasterio.enums import Resampling
    from rasterio.vrt import WarpedVRT
    from rasterio.windows import from_bounds
    want = None
    if grid_crs:
        for fn in (CRS.from_wkt, CRS.from_user_input):
            try:
                want = fn(grid_crs)
                break
            except Exception:
                pass
    if want is None and geographic:
        want = CRS.from_epsg(4326)
    with rasterio.open(path) as src:
        same = (want is None or src.crs is None or src.crs == want
                or (want.to_epsg() is not None and want.to_epsg() == src.crs.to_epsg()))
        vrt = None if same else WarpedVRT(src, crs=want, resampling=Resampling.bilinear)
        ds = vrt or src
        try:
            w = from_bounds(*bounds, ds.transform)
            bands = [1, 2, 3] if ds.count >= 3 else [1, 1, 1]
            a = ds.read(bands, window=w, out_shape=(3, out_h, out_w), boundless=vrt is None,
                        fill_value=0, resampling=Resampling.bilinear)
        finally:
            if vrt is not None:
                vrt.close()
    a = a.astype(np.float32)
    if a.max() > 255:
        hi = float(np.percentile(a[a > 0], 99.5)) if (a > 0).any() else 1.0
        a = a / max(hi, 1.0) * 255.0
    return np.ascontiguousarray(np.clip(a, 0, 255).astype(np.uint8).transpose(1, 2, 0))


def build_row_review(parent, get_grid=lambda: "", get_ortho=lambda: "", log=print):
    """The tab. get_grid / get_ortho give the current grid and the orthomosaic of tabs
    1a / 1b (used when the fields are empty or on 'Use current grid')."""
    import tkinter as tk
    from tkinter import ttk, messagebox
    from PIL import Image, ImageTk
    import plot_extractor as pe

    frm = ttk.Frame(parent)
    top = ttk.Frame(frm, padding=(8, 6, 8, 2))
    top.pack(fill="x")
    top.columnconfigure(1, weight=1)
    V = dict(grid=tk.StringVar(value=get_grid() or ""), ortho=tk.StringVar(value=get_ortho() or ""),
             rows=tk.StringVar(value=""), method=tk.StringVar(value="count"),
             n=tk.StringVar(value="5"), pitch=tk.StringVar(value=""),
             along=tk.StringVar(value="long side"))

    def browse(var, title, types, save=False):
        from tkinter import filedialog
        f = (filedialog.asksaveasfilename if save else filedialog.askopenfilename)(
            title=title, filetypes=types)
        if f:
            var.set(os.path.normpath(f))

    def path_row(r, label, var, title, types, extra=None):
        ttk.Label(top, text=label).grid(row=r, column=0, sticky="w", pady=1)
        ttk.Entry(top, textvariable=var).grid(row=r, column=1, sticky="ew", padx=6)
        ttk.Button(top, text="Browse...", command=lambda: browse(var, title, types)).grid(
            row=r, column=2)
        if extra:
            ttk.Button(top, text=extra[0], command=extra[1]).grid(row=r, column=3, padx=(6, 0))

    path_row(0, "Plot grid *", V["grid"], "Plot grid", [("Plots", "*.shp *.geojson *.json")],
             ("Use current grid", lambda: (V["grid"].set(get_grid() or V["grid"].get()),
                                           load_all())))
    path_row(1, "RGB orthomosaic", V["ortho"], "RGB orthomosaic", [("GeoTIFF", "*.tif *.tiff")])
    path_row(2, "Crop rows file", V["rows"], "Crop rows", [("Rows", "*_rows.shp *.shp")],
             ("Load", lambda: load_all(rows_only=True)))
    mk = ttk.Frame(top)
    mk.grid(row=3, column=0, columnspan=3, sticky="w", pady=(3, 1))
    ttk.Label(mk, text="Make rows:").pack(side="left")
    ttk.Radiobutton(mk, text="Fixed spacing (row count, as at generation)", value="count",
                    variable=V["method"]).pack(side="left", padx=4)
    ttk.Radiobutton(mk, text="Detect in the orthomosaic (PhenoFuse row method)",
                    value="detect", variable=V["method"]).pack(side="left", padx=4)
    ttk.Label(mk, text="  rows per plot").pack(side="left")
    ttk.Spinbox(mk, from_=0, to=60, width=4, textvariable=V["n"]).pack(side="left")
    ttk.Label(mk, text="(detect: 0 = as found)", foreground="#666").pack(side="left")
    ttk.Label(mk, text="  spacing cm").pack(side="left")
    ttk.Entry(mk, textvariable=V["pitch"], width=6).pack(side="left")
    ttk.Label(mk, text="(blank = width / rows)", foreground="#666").pack(side="left")
    ttk.Label(mk, text="  rows run along the plot's").pack(side="left")
    ttk.Combobox(mk, textvariable=V["along"], values=["long side", "short side"], width=10,
                 state="readonly").pack(side="left")
    mk_btn = ttk.Button(top, text="Make rows")
    mk_btn.grid(row=3, column=3, sticky="e", padx=(6, 0))

    tb = ttk.Frame(frm, padding=(8, 0))
    tb.pack(fill="x")
    body = ttk.Panedwindow(frm, orient="horizontal")
    body.pack(fill="both", expand=True, padx=8, pady=4)
    left = ttk.Frame(body)
    lst = tk.Listbox(left, width=26, exportselection=False, font=("Consolas", 9))
    sb = ttk.Scrollbar(left, command=lst.yview)
    lst.configure(yscrollcommand=sb.set)
    lst.pack(side="left", fill="both", expand=True)
    sb.pack(side="left", fill="y")
    body.add(left, weight=0)
    cv = tk.Canvas(body, bg="#181818", highlightthickness=0, takefocus=1)
    body.add(cv, weight=1)
    status = ttk.Label(frm, anchor="w", padding=(8, 2))
    status.pack(fill="x")
    pbar = ttk.Progressbar(frm, mode="determinate")
    pbar.pack(fill="x", padx=8, pady=(0, 4))

    S = dict(rings=[], ids=[], crs=None, M=None, frames=[], rows=[], report={}, sel=set(),
             undo=[], add=False, scale=1.0, ox=0.0, oy=0.0, base=None, detail=None,
             photo=None, dphoto=None, dkey=None, drag=None, dirty=False, cur=None, busy=False)

    # ---------------- coordinates ----------------
    def m2s(q):
        return q[0] * S["scale"] + S["ox"], S["oy"] - q[1] * S["scale"]

    def s2m(sx, sy):
        return np.array([(sx - S["ox"]) / S["scale"], (S["oy"] - sy) / S["scale"]])

    def w2s(p):
        return m2s(S["M"].to_m(p))

    def view_bounds_m():
        w, h = max(cv.winfo_width(), 50), max(cv.winfo_height(), 50)
        a, b = s2m(0, h), s2m(w, 0)
        return a[0], a[1], b[0], b[1]

    def fit(box=None):
        if not S["rings"]:
            return
        if box is None:
            pts = np.array([S["M"].to_m(p) for r in S["rings"] for p in r])
            x0, y0 = pts.min(0)
            x1, y1 = pts.max(0)
        else:
            x0, y0, x1, y1 = box
        w, h = max(cv.winfo_width(), 400), max(cv.winfo_height(), 300)
        S["scale"] = min((w - 40) / max(x1 - x0, 1e-6), (h - 40) / max(y1 - y0, 1e-6))
        S["ox"] = w / 2 - (x0 + x1) / 2 * S["scale"]
        S["oy"] = h / 2 + (y0 + y1) / 2 * S["scale"]
        redraw()

    # ---------------- loading ----------------
    def along():
        return "short" if V["along"].get().startswith("short") else "long"

    def load_grid():
        g = V["grid"].get().strip()
        if not (g and os.path.exists(g)):
            raise RuntimeError("Choose the plot grid.")
        plots, crs, fields = pe.read_plot_file(g)
        plots = [p for p in plots if len(p["ring"]) >= 3]
        idf = pe.pick_id_field(fields)
        S["rings"] = [p["ring"] for p in plots]
        S["ids"] = [pe.safe_name(p["props"].get(idf, i + 1)) if idf else str(i + 1)
                    for i, p in enumerate(plots)]
        S["crs"] = crs
        S["M"], S["frames"] = cr.grid_frames(S["rings"], along())
        S["fidx"] = {pid: i for i, pid in enumerate(S["ids"])}

    def attach(rows):
        """v_off_m / angle_deg of loaded rows from their lines, in this grid's frames."""
        keep, lost = [], 0
        for r in rows:
            i = S["fidx"].get(str(r["Plot_ID"]))
            if i is None:
                lost += 1
                continue
            c, a_u, a_v, hl, hw = S["frames"][i]
            p0, p1 = (S["M"].to_m(p) for p in r["line"])
            mid, d = (p0 + p1) / 2 - c, p1 - p0
            if float(d @ a_u) < 0:
                d = -d
            r["v_off_m"] = round(float(mid @ a_v), 4)
            r["angle_deg"] = round(math.degrees(math.atan2(float(d @ a_v), float(d @ a_u))), 3)
            keep.append(r)
        if lost:
            log(f"[!] {lost} row line(s) belong to plots not in this grid - left out")
        cr.renumber(keep)
        return keep

    def load_all(rows_only=False):
        try:
            if not rows_only or not S["rings"]:
                load_grid()
            rp = V["rows"].get().strip()
            if not rp or (not rows_only and not os.path.exists(rp)):
                rp = cr.rows_path_for(V["grid"].get().strip())
                V["rows"].set(rp)
            S["rows"], S["report"], S["sel"], S["undo"] = [], {}, set(), []
            if os.path.exists(rp):
                S["rows"] = attach(cr.read_rows(rp))
                say(f"Loaded {len(S['rows'])} row lines from {os.path.basename(rp)}")
            else:
                say("No rows file yet: choose a method and press 'Make rows'.")
            S["dirty"] = False
            load_base()
            fill_list()
            fit()
        except Exception as e:
            messagebox.showerror("Review crop rows", str(e))

    def load_base():
        S["base"] = S["detail"] = None
        o = V["ortho"].get().strip()
        if not (o and os.path.exists(o) and S["rings"]):
            return
        pts = np.array([p for r in S["rings"] for p in r], float)
        x0, y0 = pts.min(0)
        x1, y1 = pts.max(0)
        mx, my = 0.08 * (x1 - x0), 0.08 * (y1 - y0)
        b = (x0 - mx, y0 - my, x1 + mx, y1 + my)
        M = S["M"]
        wm, hm = (b[2] - b[0]) * M.mx, (b[3] - b[1]) * M.my
        f = 2400.0 / max(wm, hm)
        ow, oh = max(1, int(wm * f)), max(1, int(hm * f))

        def job():
            try:
                img = _read_ortho(o, S["crs"], b, ow, oh, M.geographic)
                S["base"] = (Image.fromarray(img), b)
                frm.after(0, redraw)
            except Exception as e:
                frm.after(0, lambda m=str(e): say(f"[!] orthomosaic: {m}"))
        threading.Thread(target=job, daemon=True).start()

    # ---------------- detail layer ----------------
    def schedule_detail():
        if S.get("dtimer"):
            frm.after_cancel(S["dtimer"])
        S["dtimer"] = frm.after(220, update_detail)

    def update_detail():
        S["dtimer"] = None
        o = V["ortho"].get().strip()
        if not (o and os.path.exists(o) and S["M"]):
            return
        x0, y0, x1, y1 = view_bounds_m()
        if x1 - x0 > DETAIL_MAX_M:
            S["detail"] = None
            return
        M = S["M"]
        b = (x0 / M.mx, y0 / M.my, x1 / M.mx, y1 / M.my)
        key = tuple(round(v, 9) for v in b)
        if S["dkey"] == key:
            return
        S["dkey"] = key
        w, h = max(cv.winfo_width(), 50), max(cv.winfo_height(), 50)

        def job():
            try:
                img = _read_ortho(o, S["crs"], b, w, h, M.geographic)
                if S["dkey"] == key:
                    S["detail"] = (Image.fromarray(img), b)
                    frm.after(0, redraw)
            except Exception:
                pass
        threading.Thread(target=job, daemon=True).start()

    # ---------------- drawing ----------------
    def redraw():
        cv.delete("all")
        if not S["rings"]:
            cv.create_text(20, 20, anchor="nw", fill="#bbb", text=(
                "Choose the plot grid (and the RGB orthomosaic) and press 'Use current grid', "
                "or load a rows file."))
            return
        cw, ch = max(cv.winfo_width(), 50), max(cv.winfo_height(), 50)
        layer = S["detail"] if S["detail"] is not None else S["base"]
        if layer is not None:
            img, b = layer
            M = S["M"]
            sx0, sy0 = m2s((b[0] * M.mx, b[3] * M.my))
            sx1, sy1 = m2s((b[2] * M.mx, b[1] * M.my))
            fw, fh = sx1 - sx0, sy1 - sy0
            iw, ih = img.size
            vx0, vy0, vx1, vy1 = max(sx0, 0), max(sy0, 0), min(sx1, cw), min(sy1, ch)
            if fw > 1 and fh > 1 and vx1 > vx0 and vy1 > vy0:
                box = (int((vx0 - sx0) / fw * iw), int((vy0 - sy0) / fh * ih),
                       max(int(math.ceil((vx1 - sx0) / fw * iw)), 1),
                       max(int(math.ceil((vy1 - sy0) / fh * ih)), 1))
                if box[2] > box[0] and box[3] > box[1]:
                    dx0, dy0 = sx0 + box[0] / iw * fw, sy0 + box[1] / ih * fh
                    dx1, dy1 = sx0 + box[2] / iw * fw, sy0 + box[3] / ih * fh
                    dw, dh = max(1, int(round(dx1 - dx0))), max(1, int(round(dy1 - dy0)))
                    if dw < 6 * cw and dh < 6 * ch:
                        S["photo"] = ImageTk.PhotoImage(img.crop(box).resize((dw, dh),
                                                                            Image.BILINEAR))
                        cv.create_image(int(round(dx0)), int(round(dy0)), anchor="nw",
                                        image=S["photo"])
        for i, r in enumerate(S["rings"]):
            pts = [c for p in r for c in w2s(p)]
            cur = S["cur"] == i
            cv.create_polygon(pts, outline="#ffffff" if cur else "#5fae3f", fill="",
                              width=2 if cur else 1)
            if S["scale"] > 8:
                c = S["frames"][i][0]
                sx, sy = m2s(c)
                cv.create_text(sx, sy, text=S["ids"][i], fill="#dddddd", font=("Segoe UI", 8))
        wid = 3 if S["scale"] > 150 else 2 if S["scale"] > 40 else 1
        for k, r in enumerate(S["rows"]):
            (a, b) = (w2s(p) for p in r["line"])
            col = COL["sel"] if k in S["sel"] else COL.get(r["source"], COL["edit"])
            cv.create_line(a[0], a[1], b[0], b[1], fill=col,
                           width=wid + (1 if k in S["sel"] else 0))
        if S["add"]:
            cv.create_text(10, 10, anchor="nw", fill="#80ff80",
                           text="ADD ROW: click inside a plot (A or Esc to stop)")
        status.config(text=status_text())

    def status_text():
        n = len(S["rows"])
        return (f"{len(S['rings'])} plots, {n} row lines, {len(S['sel'])} selected"
                + ("  |  unsaved changes" if S["dirty"] else "")
                + "  |  wheel = zoom, right/middle drag = pan, click = select row, drag = move, "
                  "A = add, Del = delete, E = even spacing, Ctrl+Z = undo, PgUp/PgDn = plot")

    def say(m):
        status.config(text=m)
        log(m)

    # ---------------- plots list ----------------
    def counts():
        c = {pid: 0 for pid in S["ids"]}
        for r in S["rows"]:
            c[r["Plot_ID"]] = c.get(r["Plot_ID"], 0) + 1
        return c

    def fill_list():
        lst.delete(0, "end")
        c = counts()
        vals = [v for v in c.values() if v]
        usual = max(set(vals), key=vals.count) if vals else 0
        for pid in S["ids"]:
            n = c.get(pid, 0)
            note = S["report"].get(pid, "")
            flag = "  <-" if (n != usual or note) else ""
            lst.insert("end", f"{pid:>8}  {n:2d} rows{flag}")
            if flag:
                lst.itemconfig("end", foreground="#c04000")

    def zoom_plot(i):
        if not (0 <= i < len(S["rings"])):
            return
        S["cur"] = i
        pts = np.array([S["M"].to_m(p) for p in S["rings"][i]])
        (x0, y0), (x1, y1) = pts.min(0), pts.max(0)
        m = 0.25 * max(x1 - x0, y1 - y0)
        fit((x0 - m, y0 - m, x1 + m, y1 + m))
        lst.selection_clear(0, "end")
        lst.selection_set(i)
        lst.see(i)
        schedule_detail()

    lst.bind("<<ListboxSelect>>", lambda e: zoom_plot(lst.curselection()[0])
             if lst.curselection() else None)

    def step_plot(d):
        zoom_plot(((S["cur"] if S["cur"] is not None else -1) + d) % max(len(S["rings"]), 1))

    # ---------------- editing ----------------
    def snapshot():
        S["undo"].append([dict(r, line=list(r["line"])) for r in S["rows"]])
        S["undo"] = S["undo"][-40:]

    def undo():
        if S["undo"]:
            S["rows"] = S["undo"].pop()
            S["sel"].clear()
            S["dirty"] = True
            fill_list()
            redraw()

    def frame_of(r):
        return S["frames"][S["fidx"][r["Plot_ID"]]]

    def relayout(rows):
        for r in rows:
            cr.relayout_plot([r], S["M"], frame_of(r))

    def hit_row(sx, sy, tol=7):
        best, bd = None, tol
        for k, r in enumerate(S["rows"]):
            (ax, ay), (bx, by) = (w2s(p) for p in r["line"])
            dx, dy = bx - ax, by - ay
            L2 = dx * dx + dy * dy or 1e-9
            t = max(0.0, min(1.0, ((sx - ax) * dx + (sy - ay) * dy) / L2))
            d = math.hypot(sx - (ax + t * dx), sy - (ay + t * dy))
            if d < bd:
                best, bd = k, d
        return best

    def plot_at(q):
        for i, (c, a_u, a_v, hl, hw) in enumerate(S["frames"]):
            d = q - c
            if abs(d @ a_u) <= hl and abs(d @ a_v) <= hw:
                return i
        return None

    def move_sel(dv):
        for k in S["sel"]:
            r = S["rows"][k]
            r["v_off_m"] = round(r["v_off_m"] + dv, 4)
            if r["source"] != "edit":
                r["source"] = "edit"
        relayout([S["rows"][k] for k in S["sel"]])
        S["dirty"] = True

    def on_press(e):
        cv.focus_set()
        q = s2m(e.x, e.y)
        if S["add"]:
            i = plot_at(q)
            if i is None:
                return
            snapshot()
            c, a_u, a_v, hl, hw = S["frames"][i]
            pid = S["ids"][i]
            same = [r for r in S["rows"] if r["Plot_ID"] == pid]
            ang = float(np.median([r["angle_deg"] for r in same])) if same else 0.0
            r = dict(Plot_ID=pid, row_idx=0, v_off_m=round(float((q - c) @ a_v), 4),
                     length_m=round(2 * hl, 3), angle_deg=ang, source="edit",
                     contrast=float("nan"), line=[])
            cr.relayout_plot([r], S["M"], S["frames"][i])
            S["rows"].append(r)
            cr.renumber(S["rows"])
            S["dirty"] = True
            fill_list()
            redraw()
            return
        k = hit_row(e.x, e.y)
        if k is None:
            if not (e.state & 0x1):
                S["sel"].clear()
            redraw()
            return
        if e.state & 0x1:
            S["sel"] ^= {k}
        elif k not in S["sel"]:
            S["sel"] = {k}
        snapshot()
        S["drag"] = dict(q=q, moved=False)
        redraw()

    def on_motion(e):
        d = S["drag"]
        if not d or not S["sel"]:
            return
        q = s2m(e.x, e.y)
        r0 = S["rows"][next(iter(S["sel"]))]
        a_v = frame_of(r0)[2]
        dv = float((q - d["q"]) @ a_v)
        d["q"] = q
        d["moved"] = True
        move_sel(dv)
        redraw()

    def on_release(_e):
        d, S["drag"] = S["drag"], None
        if d and not d["moved"] and S["undo"]:
            S["undo"].pop()                 # a click, not a move
        elif d:
            cr.renumber(S["rows"])

    def nudge(dv):
        if not S["sel"]:
            return False
        snapshot()
        move_sel(dv)
        cr.renumber(S["rows"])
        redraw()
        return True

    def delete_sel():
        if not S["sel"]:
            return
        snapshot()
        S["rows"] = [r for k, r in enumerate(S["rows"]) if k not in S["sel"]]
        S["sel"].clear()
        cr.renumber(S["rows"])
        S["dirty"] = True
        fill_list()
        redraw()

    def even():
        pids = {S["rows"][k]["Plot_ID"] for k in S["sel"]}
        if not pids and S["cur"] is not None:
            pids = {S["ids"][S["cur"]]}
        if not pids:
            say("Select a row (or a plot in the list) first.")
            return
        snapshot()
        for pid in pids:
            rr = [r for r in S["rows"] if r["Plot_ID"] == pid]
            cr.even_spacing(rr)
            relayout(rr)
        cr.renumber(S["rows"])
        S["dirty"] = True
        redraw()

    def toggle_add(v=None):
        S["add"] = (not S["add"]) if v is None else v
        add_btn.state(["pressed"] if S["add"] else ["!pressed"])
        redraw()

    # ---------------- make / save ----------------
    def make_rows():
        if S["busy"]:
            return
        try:
            load_grid()
        except Exception as e:
            messagebox.showerror("Make rows", str(e))
            return
        if S["rows"] and S["dirty"] and not messagebox.askyesno(
                "Make rows", "Replace the rows shown (unsaved changes are lost)?"):
            return
        try:
            n = int(V["n"].get() or 0)
            pitch = float(V["pitch"].get()) / 100.0 if V["pitch"].get().strip() else None
        except ValueError:
            messagebox.showerror("Make rows", "Rows per plot and spacing must be numbers.")
            return
        if V["method"].get() == "count":
            if n <= 0:
                messagebox.showerror("Make rows", "Set the number of rows per plot.")
                return
            rows, pitches = cr.rows_from_count(S["rings"], S["ids"], n, pitch, along())
            S["rows"], S["report"], S["sel"], S["undo"] = rows, {}, set(), []
            S["dirty"] = True
            p = np.nanmedian(pitches) * 100
            say(f"{len(rows)} rows: {n} per plot, spacing {p:.2f} cm "
                f"({'given' if pitch else 'plot width / rows'}). Review, then Save rows.")
            if not V["rows"].get().strip():
                V["rows"].set(cr.rows_path_for(V["grid"].get().strip()))
            fill_list()
            redraw()
            return
        o = V["ortho"].get().strip()
        if not (o and os.path.exists(o)):
            messagebox.showerror("Make rows", "Detecting rows needs the RGB orthomosaic.")
            return
        S["busy"] = True
        mk_btn.state(["disabled"])
        pbar.configure(value=0, maximum=len(S["rings"]))
        rings, ids, crs, al = S["rings"], S["ids"], S["crs"], along()

        def job():
            try:
                rows, rep = cr.detect_rows_in_grid(
                    rings, ids, o, crs, rows_expected=n, pitch_m=pitch, along=al,
                    log=lambda m: frm.after(0, lambda: say(m)),
                    progress=lambda a, b: frm.after(0, lambda: pbar.configure(value=a)))
                frm.after(0, lambda: done(rows, rep))
            except Exception as e:
                frm.after(0, lambda m=str(e): (messagebox.showerror("Detect rows", m), done(None, None)))
        threading.Thread(target=job, daemon=True).start()

    def done(rows, rep):
        S["busy"] = False
        mk_btn.state(["!disabled"])
        if rows is None:
            return
        S["rows"], S["sel"], S["undo"] = rows, set(), []
        S["report"] = {d["Plot_ID"]: d["reason"] for d in rep if d["reason"]}
        S["last_report"] = rep
        S["dirty"] = True
        ok = [d for d in rep if d["ok"]]
        pit = [d["pitch_m"] for d in ok if d["pitch_m"] != ""]
        say(f"Rows detected in {len(ok)}/{len(rep)} plots, {len(rows)} lines, median spacing "
            f"{np.median(pit) * 100:.1f} cm" if pit else "No rows detected.")
        if not V["rows"].get().strip():
            V["rows"].set(cr.rows_path_for(V["grid"].get().strip()))
        fill_list()
        redraw()

    def save():
        p = V["rows"].get().strip() or cr.rows_path_for(V["grid"].get().strip())
        if not S["rows"]:
            messagebox.showinfo("Save rows", "No rows to save.")
            return
        try:
            if os.path.exists(p) and not S.get("backed_up"):
                b = os.path.splitext(p)[0] + "_before_edit"
                if not os.path.exists(b + ".shp"):
                    for ext in (".shp", ".shx", ".dbf", ".prj", ".csv"):
                        if os.path.exists(os.path.splitext(p)[0] + ext):
                            shutil.copy2(os.path.splitext(p)[0] + ext, b + ext)
                S["backed_up"] = True
            cr.renumber(S["rows"])
            cr.write_rows(p, S["rows"], prj_from=V["grid"].get().strip(),
                          report=S.get("last_report"))
            V["rows"].set(p)
            S["dirty"] = False
            say(f"Saved {len(S['rows'])} row lines: {p} (+ {os.path.basename(os.path.splitext(p)[0])}.csv)")
        except Exception as e:
            messagebox.showerror("Save rows", str(e))
        redraw()

    mk_btn.configure(command=make_rows)

    def tool(text, cmd, w=None):
        b = ttk.Button(tb, text=text, command=cmd, width=w)
        b.pack(side="left", padx=2, pady=2)
        return b
    tool("Fit", lambda: fit())
    tool("< Prev plot", lambda: step_plot(-1))
    tool("Next plot >", lambda: step_plot(1))
    add_btn = tool("Add row (A)", lambda: toggle_add())
    tool("Delete (Del)", delete_sel)
    tool("Even spacing (E)", even)
    tool("Undo (Ctrl+Z)", undo)
    tool("Save rows", save)

    # ---------------- bindings ----------------
    pan = {}

    def pan_start(e):
        pan["xy"] = (e.x, e.y)

    def pan_move(e):
        if "xy" in pan:
            S["ox"] += e.x - pan["xy"][0]
            S["oy"] += e.y - pan["xy"][1]
            pan["xy"] = (e.x, e.y)
            redraw()

    def pan_end(_e):
        pan.pop("xy", None)
        schedule_detail()

    def wheel(e):
        f = 1.25 if e.delta > 0 else 1 / 1.25
        q = s2m(e.x, e.y)
        S["scale"] *= f
        S["ox"] = e.x - q[0] * S["scale"]
        S["oy"] = e.y + q[1] * S["scale"]
        redraw()
        schedule_detail()

    cv.bind("<ButtonPress-1>", on_press)
    cv.bind("<B1-Motion>", on_motion)
    cv.bind("<ButtonRelease-1>", on_release)
    for b in ("2", "3"):
        cv.bind(f"<ButtonPress-{b}>", pan_start)
        cv.bind(f"<B{b}-Motion>", pan_move)
        cv.bind(f"<ButtonRelease-{b}>", pan_end)
    cv.bind("<MouseWheel>", wheel)
    cv.bind("<Enter>", lambda e: cv.focus_set())
    cv.bind("<Configure>", lambda e: redraw())
    cv.bind("<Delete>", lambda e: delete_sel())
    cv.bind("<BackSpace>", lambda e: delete_sel())
    for k in ("a", "A"):
        cv.bind(f"<KeyPress-{k}>", lambda e: toggle_add())
    for k in ("e", "E"):
        cv.bind(f"<KeyPress-{k}>", lambda e: even())
    cv.bind("<Escape>", lambda e: (toggle_add(False), S["sel"].clear(), redraw()))
    cv.bind("<Control-z>", lambda e: undo())
    cv.bind("<Prior>", lambda e: step_plot(-1))
    cv.bind("<Next>", lambda e: step_plot(1))

    def arrow(e, sgn, horiz):
        # nudge the selected rows across their plot (1 cm, Shift 5 cm); else pan
        d = (0.05 if e.state & 0x1 else 0.01) * sgn
        if not nudge(d):
            if horiz:
                S["ox"] -= 60 * sgn
            else:
                S["oy"] += 60 * sgn
            redraw()
            schedule_detail()
    cv.bind("<Left>", lambda e: arrow(e, -1, True))
    cv.bind("<Right>", lambda e: arrow(e, 1, True))
    cv.bind("<Up>", lambda e: arrow(e, 1, False))
    cv.bind("<Down>", lambda e: arrow(e, -1, False))

    frm.api = dict(V=V, S=S, load=load_all, make=make_rows, save=save, even=even,
                   delete=delete_sel, undo=undo, zoom_plot=zoom_plot, nudge=nudge)
    frm.after(400, lambda: load_all() if V["grid"].get().strip() and
              os.path.exists(V["grid"].get().strip()) else redraw())
    return frm
