"""
Grid Designer — click the 4 corners of a field trial on an orthomosaic and
write a plot-wise grid shapefile.

Pure tkinter + rasterio + pyshp + numpy + PIL: no geopandas, no shapely, no
PyQt, so it drops straight into the PlotExtractor app without new deps.

The UI works in the RASTER's own CRS as world coordinates, so nothing is ever
reprojected — a corner clicked on screen is converted with the raster transform
and written out with the raster's WKT as the .prj sidecar.

Public API
----------
build_grid_designer(parent, on_saved=None, log=print, ortho_path='',
                    crs_hint='')   -> ttk.Frame (embeddable, e.g. Notebook tab)
open_grid_designer(parent, on_saved=None, log=print, ortho_path='',
                   crs_hint='')    -> Toplevel wrapping the same frame
"""

import json
import math
import os

# Label sequence for each supported click order. The dict keys are exactly the
# strings shown in the toolbar dropdown.
CLICK_ORDERS = {
    'BL -> BR -> TR -> TL': ('BL', 'BR', 'TR', 'TL'),
    'TL -> TR -> BR -> BL': ('TL', 'TR', 'BR', 'BL'),
}

# Anything below this many seconds for the display read is considered normal;
# slower usually means a striped GeoTIFF with no overviews.
SLOW_READ_S = 6.0


# --------------------------------------------------------------------------
# Grid maths (ported from phenoapp/core/grid_gen.py, shapely-free)
# --------------------------------------------------------------------------
def _slot_edges(n, size, gap):
    """Fractional [start, end] positions of each of `n` plots along [0, 1],
    given a repeating (plot `size`, `gap`) pattern. The alleys between plots
    keep their proportion of the total span, which is what lets fit-mode
    stretch to the clicked quad while still honouring the plot:gap ratio."""
    n = max(int(n), 1)
    size = max(float(size), 1e-9)
    gap = max(float(gap), 0.0)
    total = n * size + (n - 1) * gap
    if total <= 0:
        total = 1.0
    edges, pos = [], 0.0
    for _ in range(n):
        edges.append((pos / total, (pos + size) / total))
        pos += size + gap
    return edges


def _bilinear(bl, br, tr, tl, u, v):
    """Point at parametric (u, v) inside the quad; u runs BL->BR, v BL->TL."""
    return ((1 - u) * (1 - v) * bl[0] + u * (1 - v) * br[0]
            + u * v * tr[0] + (1 - u) * v * tl[0],
            (1 - u) * (1 - v) * bl[1] + u * (1 - v) * br[1]
            + u * v * tr[1] + (1 - u) * v * tl[1])


def _attrs(b, r):
    """Attribute dict for the plot at range index b, row index r (0-based)."""
    return {
        'Plot_ID': (b + 1) * 1000 + (r + 1),
        'B/R': "Rg%dRw%d" % (b + 1, r + 1),
        'Range': b + 1,
        'Row': r + 1,
        'Bank': b + 1,          # kept equal to Range for back-compatibility
    }


def generate_grid_from_corners(corners, n_ranges, n_rows, plot_w, plot_l,
                               gap_row, gap_range, fit_to_corners=False,
                               mpx=1.0, mpy=1.0):
    """Build the plot grid for a trial whose 4 corners are `corners`, given in
    BL, BR, TR, TL order and in world (raster CRS) coordinates.

    Returns a list of dicts: {'ring': [(x, y) x4], 'props': {...}} in
    range-major order, with props exactly Plot_ID / 'B/R' / Range / Row / Bank.

    mpx / mpy are metres per world unit on the x / y axis. They are 1.0 for a
    projected CRS; for a geographic one, pass the local metres-per-degree so
    the metre dimensions below stay true on the ground. The whole construction
    happens in that metric space and is mapped back at the end, which keeps the
    rotation a true rotation instead of a shear.

    fit_to_corners=True stretches the whole grid to fill the clicked quad
    (plot_w / plot_l / gaps then act only as ratios); False builds exact metre
    dimensions rotated onto the BL->BR direction and anchored at BL.
    """
    if corners is None or len(corners) != 4:
        raise ValueError("Need exactly 4 corners (BL, BR, TR, TL)")
    n_ranges = max(int(n_ranges), 1)
    n_rows = max(int(n_rows), 1)
    mpx = float(mpx) or 1.0
    mpy = float(mpy) or 1.0

    # Everything is computed in metres, then divided back into world units.
    bl, br, tr, tl = [(float(c[0]) * mpx, float(c[1]) * mpy) for c in corners]

    if fit_to_corners:
        u_edges = _slot_edges(n_rows, plot_w, gap_row)      # along BL->BR
        v_edges = _slot_edges(n_ranges, plot_l, gap_range)  # along BL->TL
        out = []
        for b, (v0, v1) in enumerate(v_edges):
            for r, (u0, u1) in enumerate(u_edges):
                ring = [_bilinear(bl, br, tr, tl, u0, v0),
                        _bilinear(bl, br, tr, tl, u1, v0),
                        _bilinear(bl, br, tr, tl, u1, v1),
                        _bilinear(bl, br, tr, tl, u0, v1)]
                out.append({'ring': [(p[0] / mpx, p[1] / mpy) for p in ring],
                            'props': _attrs(b, r)})
        return out

    # Exact-size mode: orientation comes from the bottom edge BL -> BR.
    rvx, rvy = br[0] - bl[0], br[1] - bl[1]
    row_len = math.hypot(rvx, rvy)
    if row_len <= 0:
        raise ValueError("BL and BR are the same point")
    bvx, bvy = tl[0] - bl[0], tl[1] - bl[1]
    if math.hypot(bvx, bvy) <= 0:
        raise ValueError("BL and TL are the same point")
    ca, sa = rvx / row_len, rvy / row_len       # cos/sin of the row bearing
    plot_w = float(plot_w)
    plot_l = float(plot_l)
    gap_row = float(gap_row)
    gap_range = float(gap_range)

    out = []
    for b in range(n_ranges):
        for r in range(n_rows):
            x0 = r * (plot_w + gap_row)
            y0 = b * (plot_l + gap_range)
            local = [(x0, y0), (x0 + plot_w, y0),
                     (x0 + plot_w, y0 + plot_l), (x0, y0 + plot_l)]
            ring = []
            for lx, ly in local:
                mx = bl[0] + lx * ca - ly * sa
                my = bl[1] + lx * sa + ly * ca
                ring.append((mx / mpx, my / mpy))
            out.append({'ring': ring, 'props': _attrs(b, r)})
    return out


def measure_plot_size(ring, mpx=1.0, mpy=1.0):
    """(width, length) of one plot ring in METRES: width across rows is the
    p00->p10 edge, length along the range is the p00->p01 edge. Averaged with
    the opposite edge so a slightly sheared fit-mode plot still reads sanely."""
    p = [(x * mpx, y * mpy) for x, y in ring[:4]]
    if len(p) < 4:
        return 0.0, 0.0

    def d(a, b):
        return math.hypot(p[b][0] - p[a][0], p[b][1] - p[a][1])
    return (d(0, 1) + d(3, 2)) / 2.0, (d(0, 3) + d(1, 2)) / 2.0


def metres_per_unit(crs, centre_lat):
    """Metres per CRS unit on (x, y). 1.0 / 1.0 for a projected CRS; for a
    geographic one, the local degree lengths at `centre_lat`."""
    try:
        geographic = bool(crs is not None and crs.is_geographic)
    except Exception:
        geographic = False
    if not geographic:
        return 1.0, 1.0
    lat = max(-89.5, min(89.5, float(centre_lat)))
    return 111320.0 * math.cos(math.radians(lat)), 110540.0


# --------------------------------------------------------------------------
# Writers
# --------------------------------------------------------------------------
def write_grid_shapefile(path, plots, crs=None):
    """Write `plots` (as returned by generate_grid_from_corners) to a polygon
    shapefile plus a .prj holding the CRS WKT. Rings are closed on write."""
    import shapefile
    w = shapefile.Writer(path, shapeType=shapefile.POLYGON)
    try:
        w.field('Plot_ID', 'N', 10, 0)
        w.field('B/R', 'C', 24)
        w.field('Range', 'N', 6, 0)
        w.field('Row', 'N', 6, 0)
        w.field('Bank', 'N', 6, 0)
        for pl in plots:
            ring = [list(pt) for pt in pl['ring']]
            # The shapefile spec wants outer rings CLOCKWISE; a counter-clockwise
            # ring is legally a hole, and some readers do treat it that way.
            s = 0.0
            for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1]):
                s += x0 * y1 - x1 * y0
            if s > 0:                       # positive = counter-clockwise
                ring.reverse()
            if ring[0] != ring[-1]:
                ring.append(list(ring[0]))
            w.poly([ring])
            pr = pl['props']
            w.record(pr['Plot_ID'], pr['B/R'], pr['Range'],
                     pr['Row'], pr['Bank'])
    finally:
        w.close()
    if crs is not None:
        try:
            # Prefer ESRI-dialect WKT1: ArcGIS chokes on WKT2, and it matches the
            # .prj files this trial's other tools already write. But the ESRI
            # dialect renames datums ("D_WGS_1984"), and for a GEOGRAPHIC CRS
            # that no longer reads back as the same CRS - so only use it when it
            # actually round-trips, else keep the CRS's own WKT.
            wkt = crs.to_wkt()
            try:
                esri = crs.to_wkt(version="WKT1_ESRI")
                if esri and type(crs).from_wkt(esri) == crs:
                    wkt = esri
            except Exception:
                pass
            with open(os.path.splitext(path)[0] + '.prj', 'w') as f:
                f.write(wkt)
        except Exception:
            pass          # the .shp is still usable without a .prj


def write_grid_geojson(path, plots, crs=None):
    """Write the same plots as plain GeoJSON (no fiona/geopandas needed)."""
    feats = []
    for pl in plots:
        ring = [[float(x), float(y)] for x, y in pl['ring']]
        if ring[0] != ring[-1]:
            ring.append(list(ring[0]))
        feats.append({'type': 'Feature', 'properties': dict(pl['props']),
                      'geometry': {'type': 'Polygon', 'coordinates': [ring]}})
    doc = {'type': 'FeatureCollection', 'features': feats}
    if crs is not None:
        try:
            epsg = crs.to_epsg()
            name = 'urn:ogc:def:crs:EPSG::%d' % epsg if epsg else crs.to_wkt()
        except Exception:
            name = None
        if name:
            doc['crs'] = {'type': 'name', 'properties': {'name': name}}
    with open(path, 'w') as f:
        json.dump(doc, f)


def _to_uint8(arr):
    """Stretch any dtype to display bytes. 16-bit / float orthos would come out
    black under a plain astype, so scale off the 2-98 percentile."""
    import numpy as np
    a = np.asarray(arr)
    if a.dtype == np.uint8:
        return a
    a = a.astype('float32')
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return np.zeros(a.shape, dtype='uint8')
    lo, hi = np.percentile(finite, 2), np.percentile(finite, 98)
    if not (hi > lo):
        lo, hi = float(finite.min()), float(finite.max())
    if not (hi > lo):
        hi = lo + 1.0
    a = np.clip((a - lo) / (hi - lo), 0.0, 1.0) * 255.0
    return np.nan_to_num(a).astype('uint8')


class _Designer(object):
    """Handle on a built designer UI: state containers plus the callables the
    toolbar uses, so a host app (or a test) can drive it without a mouse.
    Supports attribute and dict-style access."""

    def __init__(self, **kw):
        self.__dict__.update(kw)

    def __getitem__(self, key):
        return self.__dict__[key]

    def __contains__(self, key):
        return key in self.__dict__


# --------------------------------------------------------------------------
# GUI
# --------------------------------------------------------------------------
def build_grid_designer(parent, on_saved=None, log=print, ortho_path='',
                        crs_hint=''):
    """Build the grid designer as a ttk.Frame child of `parent` and return it.

    The frame is NOT packed or gridded — the caller places it, so the same UI
    works as a Notebook tab or inside a Toplevel. Every key binding is on the
    frame's own canvas (which grabs focus on mouse-enter) rather than on the
    root window, so the shortcuts cannot collide with the host app's bindings.

    parent      : any tk container widget
    on_saved    : optional callback, called with the saved path by "Save & use"
    log         : message sink, defaults to print
    ortho_path  : optional GeoTIFF to load as soon as the frame is realised
    crs_hint    : shown in the status bar when the raster carries no CRS

    Returns the frame with the controller attached as `frame.designer`.
    """
    import tkinter as tk
    from tkinter import ttk, messagebox, filedialog

    frame = ttk.Frame(parent)
    bar = ttk.Frame(frame, padding=5)
    bar.pack(fill='x')
    params = ttk.Frame(frame, padding=(6, 0, 6, 4))
    params.pack(fill='x')
    # takefocus so the canvas can own the keyboard while the pointer is over it
    cv = tk.Canvas(frame, bg='#181818', highlightthickness=0, takefocus=1)
    cv.pack(fill='both', expand=True)
    status = ttk.Label(frame, text="", anchor='w', padding=4)
    status.pack(fill='x')

    T = {'scale': 1.0, 'ox': 0.0, 'oy': 0.0}
    # k squashes the y axis on screen so a geographic ortho is not stretched.
    k = {'v': 1.0}
    bm = {'img': None, 'bounds': None, 'photo': None, 'key': None}
    src = {'path': '', 'crs': None, 'mpx': 1.0, 'mpy': 1.0,
           'geographic': False}
    corners = []                     # world coords, in click order
    grid = []                        # last previewed plots
    drag = {'kind': None, 'last': None}
    space_held = {'v': False}
    # close_cb lets "Save & use" dismiss a hosting Toplevel; embedded in a
    # Notebook there is nothing to close, so it stays None.
    st = {'msg': '', 'suppress_dialogs': False, 'close_cb': None}

    order_var = tk.StringVar(value='BL -> BR -> TR -> TL')
    pan_mode = tk.BooleanVar(value=False)
    fit_var = tk.BooleanVar(value=True)
    v_ranges = tk.StringVar(value='10')
    v_rows = tk.StringVar(value='20')
    v_pw = tk.StringVar(value='1.00')
    v_pl = tk.StringVar(value='4.00')
    v_grow = tk.StringVar(value='0.20')
    v_grange = tk.StringVar(value='0.50')

    def _err(title, msg):
        """Report a failure through both sinks; dialogs are skippable so the
        UI can be driven headlessly."""
        log("%s: %s" % (title, msg))
        if not st['suppress_dialogs']:
            try:
                messagebox.showerror(title, msg, parent=frame)
            except Exception:
                pass

    def _info(title, msg):
        log("%s: %s" % (title, msg))
        if not st['suppress_dialogs']:
            try:
                messagebox.showinfo(title, msg, parent=frame)
            except Exception:
                pass

    def _num(var, default, cast=float):
        try:
            return cast(str(var.get()).strip())
        except (ValueError, TypeError):
            return default

    def labels():
        return CLICK_ORDERS.get(order_var.get(), ('BL', 'BR', 'TR', 'TL'))

    def ordered_corners():
        """The 4 clicked points remapped to BL, BR, TR, TL, or None."""
        if len(corners) != 4:
            return None
        by = dict(zip(labels(), corners))
        try:
            return [by['BL'], by['BR'], by['TR'], by['TL']]
        except KeyError:
            return None

    # ---- world <-> screen -------------------------------------------------
    def w2s(x, y):
        return x * T['scale'] + T['ox'], T['oy'] - y * T['scale'] * k['v']

    def s2w(sx, sy):
        return ((sx - T['ox']) / T['scale'],
                (T['oy'] - sy) / (T['scale'] * k['v']))

    def bounds():
        if bm['bounds'] is not None:
            return bm['bounds']
        pts = list(corners) + [p for pl in grid for p in pl['ring']]
        if not pts:
            return 0.0, 0.0, 1.0, 1.0
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        return (min(xs), min(ys), max(xs) + 1e-9, max(ys) + 1e-9)

    def fit():
        x0, y0, x1, y1 = bounds()
        w = cv.winfo_width() or 1100
        h = cv.winfo_height() or 660
        pad = 40
        T['scale'] = min((w - 2 * pad) / max(x1 - x0, 1e-12),
                         (h - 2 * pad) / max((y1 - y0) * k['v'], 1e-12))
        T['ox'] = pad - x0 * T['scale']
        T['oy'] = pad + y1 * T['scale'] * k['v']

    def zoom_at(sx, sy, f):
        wx, wy = s2w(sx, sy)
        T['scale'] *= f
        T['ox'] = sx - wx * T['scale']
        T['oy'] = sy + wy * T['scale'] * k['v']
        redraw()

    def pan_by(dx, dy):
        T['ox'] += dx
        T['oy'] += dy
        redraw()

    # ---- drawing ---------------------------------------------------------
    def _draw_basemap():
        """Render ONLY the canvas-visible window of the raster, at canvas
        resolution. Scaling the whole ortho to the zoom level explodes (10x on
        a 5000 px image = 50000 px wide -> MemoryError)."""
        if bm['img'] is None or bm['bounds'] is None:
            return
        from PIL import Image as _I, ImageTk
        L, B, R, Tp = bm['bounds']
        sxL, syT = w2s(L, Tp)
        sxR, syB = w2s(R, B)
        fw, fh = sxR - sxL, syB - syT          # whole basemap at this zoom
        iw, ih = bm['img'].size
        cw = cv.winfo_width() or 1100
        chv = cv.winfo_height() or 660
        vx0, vy0 = max(sxL, 0.0), max(syT, 0.0)
        vx1, vy1 = min(sxR, float(cw)), min(syB, float(chv))
        if not (fw >= 1 and fh >= 1 and vx1 - vx0 >= 1 and vy1 - vy0 >= 1):
            return
        box = (max(0, int(math.floor((vx0 - sxL) / fw * iw))),
               max(0, int(math.floor((vy0 - syT) / fh * ih))),
               min(iw, int(math.ceil((vx1 - sxL) / fw * iw))),
               min(ih, int(math.ceil((vy1 - syT) / fh * ih))))
        if not (box[2] > box[0] and box[3] > box[1]):
            return
        # Destination rect comes from the SNAPPED source box, so the image
        # stays exactly registered to the grid outlines.
        dx0, dy0 = sxL + box[0] / iw * fw, syT + box[1] / ih * fh
        dx1, dy1 = sxL + box[2] / iw * fw, syT + box[3] / ih * fh
        dw = max(1, min(int(round(dx1 - dx0)), 4 * cw))
        dh = max(1, min(int(round(dy1 - dy0)), 4 * chv))
        if bm['key'] != (box, dw, dh):
            try:
                bm['photo'] = ImageTk.PhotoImage(
                    bm['img'].crop(box).resize((dw, dh), _I.BILINEAR))
                bm['key'] = (box, dw, dh)
            except (MemoryError, OSError, ValueError):
                bm['photo'] = None
                bm['key'] = None
        if bm['photo'] is not None:
            cv.create_image(int(round(dx0)), int(round(dy0)),
                            anchor='nw', image=bm['photo'])

    def redraw():
        try:
            cv.delete('all')
            _draw_basemap()
            for pl in grid:
                pts = []
                for x, y in pl['ring']:
                    sx, sy = w2s(x, y)
                    pts += [sx, sy]
                cv.create_polygon(pts, outline='#33ddee', fill='', width=1)
            lab = labels()
            if len(corners) >= 2:
                pts = []
                for x, y in corners:
                    sx, sy = w2s(x, y)
                    pts += [sx, sy]
                if len(corners) >= 3:
                    cv.create_polygon(pts, outline='#ffcc33', fill='',
                                      dash=(5, 3))
                else:
                    cv.create_line(pts, fill='#ffcc33', dash=(5, 3))
            for i, (x, y) in enumerate(corners):
                sx, sy = w2s(x, y)
                # Constant on-screen marker size: drawn in screen units, so it
                # never balloons when zoomed in.
                cv.create_oval(sx - 7, sy - 7, sx + 7, sy + 7,
                               outline='#ffcc33', fill='#7a4a00', width=2)
                cv.create_text(sx, sy, text=str(i + 1), fill='#ffffff',
                               font=('Segoe UI', 8, 'bold'))
                cv.create_text(sx + 14, sy - 12,
                               text=lab[i] if i < len(lab) else '',
                               fill='#ffcc33', anchor='w',
                               font=('Segoe UI', 8))
            status.config(text=status_text())
        except Exception as e:                 # never let a paint kill the UI
            log("Grid designer redraw failed: %s" % e)

    def status_text():
        who = os.path.basename(src['path']) if src['path'] else 'no ortho'
        if src['crs'] is not None:
            crs_s = str(src['crs'])
        else:
            crs_s = crs_hint or 'CRS unknown'
        units = ('geographic CRS (degrees) — metre sizes converted locally'
                 if src['geographic'] else 'projected CRS (metres)')
        base = ("%s   |   %s   |   %s   |   corners %d/4 (%s)   |   "
                "plots %d   |   zoom %.3g"
                % (who, crs_s, units, len(corners),
                   order_var.get(), len(grid), T['scale']))
        if st['msg']:
            base = st['msg'] + "   |   " + base
        return base + ("   |   wheel=zoom · Space/MMB/arrows=pan · "
                       "left-click=corner · right-click=undo corner")

    def say(msg):
        st['msg'] = msg
        if msg:
            log(msg)
        status.config(text=status_text())

    # ---- ortho loading ---------------------------------------------------
    def load_ortho(path=None):
        """Load a GeoTIFF as the basemap and adopt its CRS as world space."""
        if not path:
            path = filedialog.askopenfilename(
                parent=frame, title="Orthomosaic (georeferenced GeoTIFF)",
                filetypes=[("GeoTIFF", "*.tif *.tiff"), ("All files", "*.*")])
        if not path:
            return False
        try:
            import time
            import numpy as np
            import rasterio
            from rasterio.enums import Resampling
            from PIL import Image as _I
            t0 = time.time()
            with rasterio.open(path) as ds:
                dec = max(1, int(max(ds.width, ds.height) / 2500))
                ow = max(1, ds.width // dec)
                oh = max(1, ds.height // dec)
                if ds.count >= 3:
                    idx = [1, 2, 3]           # first 3 bands = RGB
                elif ds.count == 2:
                    idx = [1, 1, 1]           # band 2 is an alpha, ignore it
                else:
                    idx = [1, 1, 1]           # single band -> greyscale
                arr = ds.read(idx, out_shape=(3, oh, ow),
                              resampling=Resampling.average)
                b = ds.bounds
                crs = ds.crs
            dt = time.time() - t0
            rgb = _to_uint8(arr)
            bm['img'] = _I.fromarray(np.transpose(rgb, (1, 2, 0)), 'RGB')
            bm['bounds'] = (b.left, b.bottom, b.right, b.top)
            bm['key'] = None
            bm['photo'] = None
            src['path'] = path
            src['crs'] = crs
            try:
                src['geographic'] = bool(crs is not None and crs.is_geographic)
            except Exception:
                src['geographic'] = False
            clat = (b.bottom + b.top) / 2.0
            src['mpx'], src['mpy'] = metres_per_unit(crs, clat)
            k['v'] = src['mpy'] / src['mpx'] if src['mpx'] else 1.0
            del corners[:]
            del grid[:]
            if dt > SLOW_READ_S:
                log("Hint: that ortho took %.1fs to read for display — it is "
                    "probably striped with no overviews. Building overviews "
                    "(gdaladdo) or a tiled copy will make this snappier." % dt)
            log("Ortho loaded: %s  (%s)" % (os.path.basename(path),
                                            crs if crs else 'no CRS'))
            if crs is None:
                log("Warning: the raster has no CRS; the saved grid will have "
                    "no .prj sidecar.")
            fit()
            say("Ortho loaded — click the 4 trial corners in %s order"
                % order_var.get())
            redraw()
            return True
        except Exception as e:
            _err("Load orthomosaic", "Could not read the orthomosaic:\n%s" % e)
            return False

    # ---- corner picking --------------------------------------------------
    def add_corner(wx, wy):
        if len(corners) >= 4:
            say("Already have 4 corners — clear or undo one first")
            return
        corners.append((wx, wy))
        lab = labels()
        say("Corner %d (%s) placed" % (len(corners), lab[len(corners) - 1]))
        redraw()

    def undo_corner():
        if corners:
            corners.pop()
            del grid[:]
            say("Removed the last corner")
        redraw()

    def clear_corners():
        del corners[:]
        del grid[:]
        say("Corners cleared")
        redraw()

    # ---- grid build ------------------------------------------------------
    def build_grid():
        """Generate plots from the current corners + parameters, or None."""
        quad = ordered_corners()
        if quad is None:
            _err("Preview grid", "Click all 4 trial corners first.")
            return None
        try:
            plots = generate_grid_from_corners(
                quad,
                _num(v_ranges, 10, int), _num(v_rows, 20, int),
                _num(v_pw, 1.0), _num(v_pl, 4.0),
                _num(v_grow, 0.2), _num(v_grange, 0.5),
                fit_to_corners=bool(fit_var.get()),
                mpx=src['mpx'], mpy=src['mpy'])
        except Exception as e:
            _err("Preview grid", "Could not build the grid:\n%s" % e)
            return None
        return plots

    def preview():
        plots = build_grid()
        if not plots:
            return False
        del grid[:]
        grid.extend(plots)
        ws, ls = [], []
        for pl in plots:
            w, l = measure_plot_size(pl['ring'], src['mpx'], src['mpy'])
            ws.append(w)
            ls.append(l)
        say("%d plots — mean plot %.3f m wide x %.3f m long (%s)"
            % (len(plots), sum(ws) / len(ws), sum(ls) / len(ls),
               'stretched to corners' if fit_var.get() else 'exact metres'))
        redraw()
        return True

    # ---- saving ----------------------------------------------------------
    def default_out():
        base = os.path.dirname(src['path']) if src['path'] else os.getcwd()
        return os.path.join(base, 'plot_grid.shp')

    def save(and_use=False, path=None):
        """Write the grid; with and_use, also call on_saved(path) and close any
        hosting window (a Notebook tab has nothing to close)."""
        plots = grid if grid else build_grid()
        if not plots:
            return None
        if not path:
            d = default_out()
            path = filedialog.asksaveasfilename(
                parent=frame, title="Save plot grid",
                initialdir=os.path.dirname(d),
                initialfile=os.path.basename(d), defaultextension='.shp',
                filetypes=[("Shapefile", "*.shp"), ("GeoJSON", "*.geojson")])
        if not path:
            return None
        try:
            if os.path.splitext(path)[1].lower() in ('.geojson', '.json'):
                write_grid_geojson(path, plots, src['crs'])
            else:
                if not path.lower().endswith('.shp'):
                    path += '.shp'
                write_grid_shapefile(path, plots, src['crs'])
        except Exception as e:
            _err("Save grid", "Could not save the grid:\n%s" % e)
            return None
        del grid[:]
        grid.extend(plots)
        log("Grid saved: %s  (%d plots)" % (path, len(plots)))
        if and_use:
            if on_saved is not None:
                try:
                    on_saved(path)
                except Exception as e:
                    _err("Save & use", "The grid was saved but the caller "
                                       "could not use it:\n%s" % e)
            if st['close_cb'] is not None:
                try:
                    st['close_cb']()
                except Exception:
                    pass
            else:
                say("Saved %d plots to %s"
                    % (len(plots), os.path.basename(path)))
            return path
        _info("Save grid", "Saved %d plots:\n%s" % (len(plots), path))
        say("Saved %d plots to %s" % (len(plots), os.path.basename(path)))
        return path

    # ---- mouse / keys ----------------------------------------------------
    def on_press(e):
        cv.focus_set()
        if pan_mode.get() or space_held['v']:
            drag['kind'] = 'pan'
            drag['last'] = (e.x, e.y)
            return
        wx, wy = s2w(e.x, e.y)
        add_corner(wx, wy)

    def on_motion(e):
        if drag['kind'] == 'pan' and drag['last']:
            lx, ly = drag['last']
            T['ox'] += e.x - lx
            T['oy'] += e.y - ly
            drag['last'] = (e.x, e.y)
            redraw()

    def on_release(_e):
        drag['kind'] = None

    def on_mmb(e):
        drag['kind'] = 'pan'
        drag['last'] = (e.x, e.y)

    def on_wheel(e):
        zoom_at(e.x, e.y, 1.2 if e.delta > 0 else 1 / 1.2)

    # ---- toolbar (each tool shows its keyboard shortcut underneath) ------
    btns = {}

    def tool(label, short, cmd, w=13):
        f = ttk.Frame(bar)
        f.pack(side='left', padx=3)
        b = ttk.Button(f, text=label, command=cmd, width=w)
        b.pack()
        ttk.Label(f, text=short, foreground='#888',
                  font=('Segoe UI', 7)).pack()
        btns[label] = b
        return b

    tool("Load orthomosaic...", "Ctrl+O", lambda: load_ortho(), w=18)
    of = ttk.Frame(bar)
    of.pack(side='left', padx=(8, 3))
    ttk.Combobox(of, textvariable=order_var, width=20, state='readonly',
                 values=list(CLICK_ORDERS.keys())).pack()
    ttk.Label(of, text="click order", foreground='#888',
              font=('Segoe UI', 7)).pack()
    pf = ttk.Frame(bar)
    pf.pack(side='left', padx=3)
    ttk.Checkbutton(pf, text="Pan mode (else place corners)",
                    variable=pan_mode).pack()
    ttk.Label(pf, text="P · Space / MMB", foreground='#888',
              font=('Segoe UI', 7)).pack()
    tool("Undo corner", "Ctrl+Z", undo_corner, w=11)
    tool("Clear corners", "Esc", clear_corners, w=12)
    tool("Fit view", "F", lambda: (fit(), redraw()), w=8)
    tool("Preview grid", "G", preview, w=11)
    tool("Save shapefile", "Ctrl+S", lambda: save(False), w=13)
    tool("Save & use", "S", lambda: save(True), w=10)

    # ---- parameters panel ------------------------------------------------
    def field(label, var, width=6):
        f = ttk.Frame(params)
        f.pack(side='left', padx=(0, 10))
        ttk.Label(f, text=label, font=('Segoe UI', 8)).pack(anchor='w')
        ttk.Spinbox(f, textvariable=var, from_=0, to=100000, increment=0.1,
                    width=width).pack(anchor='w')

    field("Ranges (count)", v_ranges)
    field("Rows per range (count)", v_rows)
    field("Plot width across rows (m)", v_pw)
    field("Plot length along range (m)", v_pl)
    field("Gap between rows (m)", v_grow)
    field("Gap between ranges (m)", v_grange)
    ttk.Checkbutton(params, text="Stretch grid to fit the 4 clicked corners",
                    variable=fit_var).pack(side='left', padx=(4, 0))

    # ---- bindings (canvas-local, so an embedded tab cannot steal the host
    # window's own shortcuts; the canvas takes focus on mouse-enter) -------
    cv.bind('<Enter>', lambda e: cv.focus_set())
    cv.bind('<ButtonPress-1>', on_press)
    cv.bind('<B1-Motion>', on_motion)
    cv.bind('<ButtonRelease-1>', on_release)
    cv.bind('<ButtonPress-3>', lambda e: undo_corner())
    cv.bind('<ButtonPress-2>', on_mmb)
    cv.bind('<B2-Motion>', on_motion)
    cv.bind('<ButtonRelease-2>', on_release)
    cv.bind('<MouseWheel>', on_wheel)
    cv.bind('<Configure>', lambda e: redraw())
    cv.bind('<Control-o>', lambda e: load_ortho())
    cv.bind('<Control-O>', lambda e: load_ortho())
    cv.bind('<Control-s>', lambda e: save(False))
    cv.bind('<Control-S>', lambda e: save(False))
    cv.bind('<Control-z>', lambda e: undo_corner())
    cv.bind('<Control-Z>', lambda e: undo_corner())
    cv.bind('<Escape>', lambda e: clear_corners())
    cv.bind('<KeyPress-space>', lambda e: space_held.__setitem__('v', True))
    cv.bind('<KeyRelease-space>', lambda e: space_held.__setitem__('v', False))
    for kk in ('<KeyPress-f>', '<KeyPress-F>'):
        cv.bind(kk, lambda e: (fit(), redraw()))
    for kk in ('<KeyPress-g>', '<KeyPress-G>'):
        cv.bind(kk, lambda e: preview())
    for kk in ('<KeyPress-s>', '<KeyPress-S>'):
        cv.bind(kk, lambda e: save(True))
    for kk in ('<KeyPress-p>', '<KeyPress-P>'):
        cv.bind(kk, lambda e: pan_mode.set(not pan_mode.get()))
    cv.bind('<Left>', lambda e: pan_by(60, 0))
    cv.bind('<Right>', lambda e: pan_by(-60, 0))
    cv.bind('<Up>', lambda e: pan_by(0, 60))
    cv.bind('<Down>', lambda e: pan_by(0, -60))

    designer = _Designer(
        frame=frame, canvas=cv, status=status, state=st, transform=T,
        corners=corners, grid=grid, src=src, basemap=bm, buttons=btns,
        vars={'order': order_var, 'pan_mode': pan_mode, 'fit': fit_var,
              'ranges': v_ranges, 'rows': v_rows, 'plot_w': v_pw,
              'plot_l': v_pl, 'gap_row': v_grow, 'gap_range': v_grange},
        load_ortho=load_ortho, preview=preview, save=save,
        build_grid=build_grid, add_corner=add_corner,
        clear_corners=clear_corners, undo_corner=undo_corner,
        fit=fit, redraw=redraw, w2s=w2s, s2w=s2w, zoom_at=zoom_at,
        ordered_corners=ordered_corners, log=log)
    frame.designer = designer

    say("Load an orthomosaic, then click the 4 trial corners")
    if ortho_path:
        frame.after(60, lambda: load_ortho(ortho_path))
    else:
        frame.after(60, lambda: (fit(), redraw()))
    return frame


def open_grid_designer(parent, on_saved=None, log=print, ortho_path='',
                       crs_hint=''):
    """Open the grid designer in its own window: click the 4 corners of a trial
    on an orthomosaic, set the plot layout, and save a plot-wise grid shapefile
    in the raster's own CRS.

    Thin wrapper over build_grid_designer. Returns the Toplevel;
    `win.designer` / `win.gd_state` are the same controller the frame carries.
    """
    import tkinter as tk

    win = tk.Toplevel(parent)
    win.title("Generate grid — click 4 trial corners on the orthomosaic")
    win.geometry("1340x880")
    frame = build_grid_designer(win, on_saved=on_saved, log=log,
                                ortho_path=ortho_path, crs_hint=crs_hint)
    frame.pack(fill='both', expand=True)
    # "Save & use" should dismiss the window it created, but not a Notebook tab
    frame.designer.state['close_cb'] = win.destroy
    win.designer = frame.designer
    win.gd_state = frame.designer
    win.grid_frame = frame
    return win


if __name__ == '__main__':
    import tkinter as tk

    _root = tk.Tk()
    _root.withdraw()
    _w = open_grid_designer(_root, on_saved=lambda p: print("saved:", p))
    _w.protocol('WM_DELETE_WINDOW', _root.destroy)
    _root.mainloop()
