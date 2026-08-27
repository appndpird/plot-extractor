"""PlotExtractor — simple desktop tool to extract one raw image (and optionally a
deskewed orthomosaic crop) per plot from a Metashape project.

Inputs: Metashape project (.psx), plot boundaries (.shp or .geojson), optional raw
image folder (fallback when the project's photo paths are stale), output folder,
margin, images per plot, limit.

Runs on any Python 3.9+ with: Metashape (standalone wheel), Pillow, pyshp.
Launch with PlotExtractor.bat or:  python plot_extractor.py
"""
import os, sys, re, json, math, csv, threading, queue, traceback

# Metashape standalone module license (floating license server config lives here)
os.environ.setdefault('agisoft_LICENSE', r"C:\Program Files\Agisoft\Licensing\licenses")

# When frozen (PyInstaller .exe) config.json lives NEXT TO the exe (writable);
# from source it lives next to this script.
if getattr(sys, 'frozen', False):
    APP_DIR = os.path.dirname(sys.executable)
else:
    APP_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.join(APP_DIR, "config.json")
RAW_EXTS = ['.jpg', '.jpeg', '.tif', '.tiff', '.png']
AUTO_ID = "(auto: plot id field or 1..N)"


# --------------------------------------------------------------------------
# Plot boundary file reading (no Metashape required)
# --------------------------------------------------------------------------
def _ring_area(ring):
    s = 0.0
    for (x0, y0), (x1, y1) in zip(ring, ring[1:] + ring[:1]):
        s += x0 * y1 - x1 * y0
    return abs(s) / 2


def _ring_area_m2(ring_wgs):
    """Area (m^2) of a WGS84 lon/lat ring via local metres-per-degree shoelace."""
    if len(ring_wgs) < 3:
        return 0.0
    cy = sum(q[1] for q in ring_wgs) / len(ring_wgs)
    mlon = 111320 * math.cos(math.radians(cy))
    mlat = 110540
    pts = [(lon * mlon, lat * mlat) for lon, lat in ring_wgs]
    s = 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:] + pts[:1]):
        s += x0 * y1 - x1 * y0
    return abs(s) / 2


def _point_in_ring(x, y, ring):
    """Ray-cast point-in-polygon for an open ring of (x,y)."""
    inside = False
    n = len(ring)
    for i in range(n):
        x0, y0 = ring[i]
        x1, y1 = ring[(i + 1) % n]
        if (y0 > y) != (y1 > y):
            xint = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
            if x < xint:
                inside = not inside
    return inside


def _clip_poly_rect(poly, w, h):
    """Sutherland-Hodgman clip of a polygon [(x,y), ...] to the rectangle
    [0,w] x [0,h]. Returns the clipped polygon (possibly empty)."""
    def clip(pts, inside, inter):
        out = []
        n = len(pts)
        for i in range(n):
            a = pts[i]
            b = pts[(i + 1) % n]
            ia, ib = inside(a), inside(b)
            if ia:
                out.append(a)
                if not ib:
                    out.append(inter(a, b))
            elif ib:
                out.append(inter(a, b))
        return out

    def ix(a, b, xc):                                   # crossing of a-b with x=xc
        t = (xc - a[0]) / (b[0] - a[0]) if b[0] != a[0] else 0.0
        return (xc, a[1] + t * (b[1] - a[1]))

    def iy(a, b, yc):                                   # crossing of a-b with y=yc
        t = (yc - a[1]) / (b[1] - a[1]) if b[1] != a[1] else 0.0
        return (a[0] + t * (b[0] - a[0]), yc)

    pts = list(poly)
    for inside, inter in (
            (lambda p: p[0] >= 0.0, lambda a, b: ix(a, b, 0.0)),
            (lambda p: p[0] <= w,   lambda a, b: ix(a, b, w)),
            (lambda p: p[1] >= 0.0, lambda a, b: iy(a, b, 0.0)),
            (lambda p: p[1] <= h,   lambda a, b: iy(a, b, h))):
        if not pts:
            return []
        pts = clip(pts, inside, inter)
    return pts


def frame_coverage(quad_px, w, h):
    """Fraction (0..1) of a plot quad, projected into a w x h frame, that is
    actually inside the frame. Exact polygon clip - far more informative than
    counting how many of the 4 corners land in the image (the old measure),
    which is 2/4 for both a frame holding 51% of the plot and one holding 95%."""
    a0 = _ring_area(quad_px)
    if a0 <= 1e-9:
        return 0.0
    cl = _clip_poly_rect(quad_px, float(w), float(h))
    if len(cl) < 3:
        return 0.0
    return max(0.0, min(1.0, _ring_area(cl) / a0))


# --------------------------------------------------------------------------
# Exact rectification (projected mesh) + minimum-error seam
# --------------------------------------------------------------------------
MESH_STEP_PX = 96           # target output-pixel size of one mesh cell
MESH_MAX_CELLS = 512        # per axis


def mesh_dims(wpx, hpx, step=MESH_STEP_PX):
    """Number of mesh cells across/down for an output crop of wpx x hpx.

    The cell size matters because the plot surface is NOT smooth at plot scale:
    the seed ridges are ~1-3 cm high and ~20 cm apart, so a mesh coarser than
    the ridge spacing cannot follow the relief and the two frames of a stitch
    still disagree over a ridge. Measured cross-frame misregistration on OZ
    Barley York: 15-24 px mean with one quad, 4-6 px at 256 px cells,
    ~1-2 px at 64-96 px cells."""
    nx = max(2, min(MESH_MAX_CELLS, int(math.ceil(wpx / float(step)))))
    ny = max(2, min(MESH_MAX_CELLS, int(math.ceil(hpx / float(step)))))
    return nx, ny


def mesh_uv_nodes(ring, nx, ny):
    """Ground positions of the (nx+1) x (ny+1) mesh nodes, bilinearly spread over
    the plot quad `ring` = [A, B, C, D]. u runs A->B (output x), v runs A->D
    (output y) - the same convention the QUAD warp used for its 4 corners."""
    A, B, C, D = ring
    out = []
    for j in range(ny + 1):
        v = j / float(ny)
        for i in range(nx + 1):
            u = i / float(nx)
            out.append([(1 - u) * (1 - v) * A[k] + u * (1 - v) * B[k]
                        + u * v * C[k] + (1 - u) * v * D[k] for k in (0, 1)])
    return out


def build_pil_mesh(src_nodes, nx, ny, wpx, hpx):
    """Turn a grid of source-pixel coordinates into PIL's MESH transform data.

    Each cell becomes (destination box, source quad). Because every frame is
    sampled at the SAME ground nodes, two frames rectified through their own
    meshes land on the ground identically - which is what stops a plant from
    appearing twice across a stitch boundary. Cells with a node that failed to
    project are dropped (they fall outside the frame anyway and would be black).
    """
    mesh = []
    for j in range(ny):
        y0 = int(round(hpx * j / float(ny)))
        y1 = int(round(hpx * (j + 1) / float(ny)))
        if y1 <= y0:
            continue
        for i in range(nx):
            x0 = int(round(wpx * i / float(nx)))
            x1 = int(round(wpx * (i + 1) / float(nx)))
            if x1 <= x0:
                continue
            ul = src_nodes[j * (nx + 1) + i]
            ur = src_nodes[j * (nx + 1) + i + 1]
            ll = src_nodes[(j + 1) * (nx + 1) + i]
            lr = src_nodes[(j + 1) * (nx + 1) + i + 1]
            if ul is None or ur is None or ll is None or lr is None:
                continue
            # PIL quad order: upper-left, lower-left, lower-right, upper-right
            mesh.append(((x0, y0, x1, y1),
                         (ul[0], ul[1], ll[0], ll[1], lr[0], lr[1], ur[0], ur[1])))
    return mesh


def min_error_seam(np, g_old, g_new, filled, vmask, slack=300, bias=0.05,
                   maxstep=4):
    """Choose where two frames should hand over to each other.

    `filled` = pixels the composite already holds, `vmask` = pixels the new frame
    can supply. Letting the new frame fill only the holes puts the boundary
    exactly on the old frame's footprint edge, wherever that happens to land -
    including straight through a seedling, which is then drawn twice because its
    top is displaced differently in the two views. So the boundary is allowed to
    move up to `slack` pixels back into the overlap and is routed along the path
    where the two frames look most alike, i.e. over bare soil. Classic
    minimum-error boundary cut, solved with dynamic programming.

    `slack` is deliberately small (300 px is ~10 cm here, several times a
    seedling) and `bias` pulls the path back towards the footprint edge: the
    point is to step around a plant, NOT to hand half the plot to a second
    frame. Letting the seam settle mid-overlap looks tidy but throws away the
    best frame's pixels over a large area.

    Returns a dict (take-mask, path, orientation) or None when the geometry
    isn't a clean band, in which case the caller falls back to plain hole
    filling.
    """
    ov = vmask & filled
    need = vmask & (~filled)
    if not need.any() or ov.sum() < 5000:
        return None
    # Which way does the boundary run? Count membership flips between vertical
    # vs horizontal neighbours; the seam is a function of the other axis.
    m = need.astype(np.int8) - ov.astype(np.int8)
    flips_v = int(np.count_nonzero(np.diff(m, axis=0)))     # changes down a column
    flips_h = int(np.count_nonzero(np.diff(m, axis=1)))     # changes along a row
    transposed = flips_h > flips_v
    if transposed:
        g_old = g_old.T; g_new = g_new.T
        ov = ov.T; need = need.T; vmask = vmask.T
    H, W = ov.shape
    idx = np.arange(H)[:, None]
    ov_any = ov.any(0)
    need_any = need.any(0)
    active = ov_any & need_any
    if active.sum() < 0.2 * W:
        return None
    # per-column overlap band and which side of it the new pixels sit on
    big = H + 10
    lo = np.where(ov_any, np.where(ov, idx, big).min(0), 0)
    hi = np.where(ov_any, np.where(ov, idx, -1).max(0), -1)
    need_c = np.where(need_any, np.where(need, idx, 0).sum(0) /
                      np.maximum(need.sum(0), 1), 0.0)
    ov_c = np.where(ov_any, np.where(ov, idx, 0).sum(0) /
                    np.maximum(ov.sum(0), 1), 0.0)
    side_col = np.sign(need_c - ov_c)[active]
    if side_col.size == 0:
        return None
    side = 1.0 if side_col.sum() >= 0 else -1.0        # +1: new frame is BELOW
    # The plain hole-fill boundary, per column: the last overlap row on the side
    # the new pixels come from. The seam may retreat up to `slack` from it.
    edge = hi if side > 0 else lo
    slack = max(8, int(slack))
    if side > 0:
        Lo = int(max(0, edge[active].min() - slack))
        Hi = int(min(H - 1, edge[active].max()))
    else:
        Lo = int(max(0, edge[active].min()))
        Hi = int(min(H - 1, edge[active].max() + slack))
    if Hi - Lo < 2:
        return None
    band = slice(Lo, Hi + 1)
    BH = Hi - Lo + 1
    if BH * W > 120_000_000:                            # keep the DP affordable
        return None
    cost = np.abs(g_old[band].astype(np.float32) - g_new[band].astype(np.float32))
    rows = np.arange(Lo, Hi + 1, dtype=np.float32)[:, None]
    # legal: inside the overlap AND within `slack` of this column's own edge
    eg = edge[None, :].astype(np.float32)
    if side > 0:
        legal = ov[band] & (rows <= eg) & (rows >= eg - slack)
    else:
        legal = ov[band] & (rows >= eg) & (rows <= eg + slack)
    PEN = 1e4
    cost = np.where(legal, cost, PEN)
    if bias:                    # pull back towards the footprint edge
        cost = cost + bias * np.abs(rows - eg)
    # DP left->right. The step limit has to be more than 1 row per column: a
    # seedling is only ~30-90 px across, and at +/-1 the path needs ~50 columns
    # to descend past it and 50 more to climb back, so it could not detour around
    # anything narrower than ~100 px - exactly the obstacles that matter here.
    ms = max(1, int(maxstep))
    offs = list(range(-ms, ms + 1))                 # back value = prev row - cur row
    acc = np.empty((BH, W), np.float32)
    back = np.zeros((BH, W), np.int8)
    acc[:, 0] = cost[:, 0]
    stk = np.empty((len(offs), BH), np.float32)
    ar = np.arange(BH)
    offa = np.asarray(offs, np.int8)
    for c in range(1, W):
        prev = acc[:, c - 1]
        for i, d in enumerate(offs):
            col = stk[i]
            if d == 0:
                col[:] = prev
            elif d < 0:                             # came from |d| rows above
                col[:-d] = np.inf
                col[-d:] = prev[:d]
            else:                                   # came from d rows below
                col[-d:] = np.inf
                col[:-d] = prev[d:]
        k = np.argmin(stk, 0)
        acc[:, c] = stk[k, ar] + cost[:, c]
        back[:, c] = offa[k]
    path = np.empty(W, np.int32)
    path[W - 1] = int(np.argmin(acc[:, W - 1]))
    for c in range(W - 1, 0, -1):
        path[c - 1] = path[c] + int(back[path[c], c])
    path = np.clip(path + Lo, 0, H - 1)
    # columns with no overlap keep the plain hole-fill boundary
    take = need.copy()
    if side > 0:
        seam_take = idx > path[None, :]
    else:
        seam_take = idx < path[None, :]
    take |= (ov & seam_take & active[None, :])
    take &= vmask
    return dict(take=(take.T if transposed else take), path=path,
                transposed=transposed, side=side, active=active,
                cut=float(np.mean(np.abs(np.diff(path)))))


def feather_seam(np, out, ra, seam, width=600, halfband=24, block=256):
    """Erase the residual brightness/colour LINE along a stitch boundary WITHOUT
    mixing pixels from two frames.

    The step across the seam is measured column-by-column in a thin band on
    either side, smoothed hard along the seam, and then added to the incoming
    frame as an offset that fades to zero `width` pixels away. Only a smooth
    low-frequency offset is applied - every output pixel keeps 100% of one
    frame's texture, so nothing is averaged and nothing is duplicated.

    `ra` (float32) is corrected IN PLACE and clipped to 0..255. Everything is
    done in row blocks: a 13000x3400 crop would otherwise need ~1 GB of
    full-size float temporaries per worker, times however many cores are busy.
    """
    path = seam['path']
    side = seam['side']
    active = seam['active']
    # transpose() returns a view, so in-place writes reach the caller's array
    o = out.transpose(1, 0, 2) if seam['transposed'] else out
    r = ra.transpose(1, 0, 2) if seam['transposed'] else ra
    H, W = o.shape[:2]
    pf = path[None, :].astype(np.float32)
    sgn = np.float32(1.0 if side >= 0 else -1.0)

    def rel_of(y0, y1):
        rows = np.arange(y0, y1, dtype=np.float32)[:, None]
        return (rows - pf) * sgn                    # <0 = old side, >0 = new side

    # --- measure the step across the seam, per column, in thin bands ---
    no = np.zeros(W, np.float32); nn = np.zeros(W, np.float32)
    so = np.zeros((3, W), np.float32); sn = np.zeros((3, W), np.float32)
    for y0 in range(0, H, block):
        y1 = min(H, y0 + block)
        rel = rel_of(y0, y1)
        ob = (rel < 0) & (rel >= -halfband)
        nb = (rel >= 0) & (rel <= halfband)
        if not (ob.any() or nb.any()):
            continue
        no += ob.sum(0); nn += nb.sum(0)
        for c in range(3):
            so[c] += (o[y0:y1, :, c] * ob).sum(0)
            sn[c] += (r[y0:y1, :, c] * nb).sum(0)
    good = (no > 4) & (nn > 4) & active
    step = np.zeros((3, W), np.float32)
    if good.any():
        k = max(3, int(W / 12) | 1)
        ker = np.ones(k, np.float32)
        den = np.convolve(good.astype(np.float32), ker, 'same')
        for c in range(3):
            d = np.zeros(W, np.float32)
            d[good] = so[c][good] / no[good] - sn[c][good] / nn[good]
            # smooth hard along the seam, and carry the value into dead columns
            num = np.convolve(np.where(good, d, 0.0), ker, 'same')
            step[c] = np.clip(np.where(den > 0, num / np.maximum(den, 1e-6), 0.0),
                              -60.0, 60.0)

    # --- apply it, fading out away from the seam ---
    for y0 in range(0, H, block):
        y1 = min(H, y0 + block)
        ramp = np.clip(1.0 - np.maximum(rel_of(y0, y1), 0.0) / float(width),
                       0.0, 1.0)
        for c in range(3):
            r[y0:y1, :, c] += step[c][None, :] * ramp
    np.clip(ra, 0, 255, out=ra)
    return ra


def _tilt_from_vertical(pu):
    """Acute angle (deg, 0-90) between a plot's long edge and the image vertical,
    from its 4 corners [A,B,C,D] in pixel coords. 0 = already upright."""
    if len(pu) < 4:
        return 0.0
    ax, ay = pu[0]
    dx, dy = pu[3]
    a = abs(math.degrees(math.atan2(dx - ax, dy - ay))) % 180.0
    return min(a, 180.0 - a)


def _geojson_crs(data):
    """Extract a Metashape-parseable CRS spec from a GeoJSON 'crs' member, if any.
    Handles 'EPSG:xxxx', 'EPSG::xxxx' and 'urn:ogc:def:crs:EPSG::xxxx'."""
    name = (((data.get('crs') or {}).get('properties')) or {}).get('name', '')
    m = re.search(r'EPSG:*(\d+)', str(name), re.I)
    if m:
        return f"EPSG::{m.group(1)}"
    return None


def read_plot_file(path):
    """Return (plots, crs_def, field_names).
    plots = [{ring: [[x,y],...] (open, exterior only), props: {...}}, ...]
    crs_def = WKT string or 'EPSG::xxxx'; None means 'assume WGS84 lon/lat'."""
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.geojson', '.json'):
        with open(path, encoding='utf-8') as f:
            data = json.load(f)
        plots, fields = [], []
        for ft in data.get('features', []):
            geom = ft.get('geometry') or {}
            gtype, coords = geom.get('type'), geom.get('coordinates')
            if gtype == 'Polygon':
                ring = coords[0]
            elif gtype == 'MultiPolygon':
                ring = max(coords, key=lambda poly: _ring_area(poly[0]))[0]
            else:
                continue
            ring = [list(p[:2]) for p in ring]
            if len(ring) > 1 and ring[0] == ring[-1]:
                ring = ring[:-1]
            props = ft.get('properties') or {}
            for k in props:
                if k not in fields:
                    fields.append(k)
            plots.append(dict(ring=ring, props=props))
        return plots, _geojson_crs(data), fields

    if ext == '.shp':
        import shapefile
        rdr = shapefile.Reader(path)
        fields = [f[0] for f in rdr.fields[1:]]
        plots = []
        for sr in rdr.iterShapeRecords():
            shp = sr.shape
            if 'POLYGON' not in shp.shapeTypeName.upper():
                continue
            parts = list(shp.parts) + [len(shp.points)]
            rings = [shp.points[parts[i]:parts[i + 1]] for i in range(len(parts) - 1)]
            ring = max(rings, key=_ring_area)          # exterior ring
            ring = [list(p[:2]) for p in ring]
            if len(ring) > 1 and ring[0] == ring[-1]:
                ring = ring[:-1]
            plots.append(dict(ring=ring, props=dict(zip(fields, sr.record))))
        wkt = None
        prj = os.path.splitext(path)[0] + '.prj'
        if os.path.exists(prj):
            with open(prj) as f:
                wkt = f.read().strip()
        return plots, wkt, fields

    raise ValueError("Unsupported plot boundary file (use .shp or .geojson): " + path)


def pick_id_field(fields):
    """Best guess for the plot-id attribute."""
    low = {f.lower(): f for f in fields}
    for cand in ('plot_id', 'plotid', 'plot', 'id'):
        if cand in low:
            return low[cand]
    return None


def safe_name(val):
    s = str(val).strip()
    return "".join(c if (c.isalnum() or c in '-_.') else '_' for c in s) or "unnamed"


# --------------------------------------------------------------------------
# Extraction engine (worker thread)
# --------------------------------------------------------------------------
def _periodicity(prof):
    """Strength of the dominant ridge/row oscillation in a 1-D profile."""
    import numpy as np
    n = len(prof)
    if n < 32:
        return 0.0
    x = np.arange(n, dtype='float64')
    p = prof - np.polyval(np.polyfit(x, prof, 3), x)      # drop illumination trend
    f = np.abs(np.fft.rfft((p - p.mean()) * np.hanning(n)))
    lo, hi = max(2, n // 120), max(3, n // 5)
    if hi <= lo:
        return 0.0
    band = f[lo:hi]
    return float(band.max() / (band.mean() + 1e-9))


def _across_row_axis(im, hint=None):
    """Which image axis runs ACROSS the planted rows: 'x' (rows vertical) or
    'y' (rows horizontal). Rows run along the plot's LONG axis, so the
    across-row axis is the SHORT image axis - which side that is depends on the
    shapefile's vertex order, so it must be measured, not assumed. Near-square
    crops fall back to whichever axis carries the stronger row periodicity."""
    if hint in ('x', 'y'):
        return hint
    W, H = im.size
    ar = W / float(H or 1)
    if ar <= 0.8:
        return 'x'
    if ar >= 1.25:
        return 'y'
    import numpy as np
    a = np.asarray(im.resize((min(W, 700), min(H, 700))).convert('L')).astype('float32')
    return 'x' if _periodicity(a.mean(0)) >= _periodicity(a.mean(1)) else 'y'


def _refit_gcps_box(gcps, x0, y0, x1, y1):
    """Rebuild 4 corner GCPs for a sub-window (x0,y0,x1,y1) of the current frame.

    Same bilinear refit used when auto-fit trims a crop: fit
    v = a + b*row + c*col + d*row*col through the existing corner GCPs, then
    evaluate at the new corners. Filtering GCPs by whether they fall inside the
    window is NOT equivalent - the corners always sit outside it, so filtering
    leaves none at all.
    """
    if len(gcps) < 4:
        return [(r - y0, c - x0, lon, lat, z) for (r, c, lon, lat, z) in gcps]
    A = [[1.0, g[0], g[1], g[0] * g[1]] for g in gcps]
    coefs = []
    for k in (2, 3, 4):
        coefs.append(_lstsq4(A, [g[k] for g in gcps]))
    out = []
    for r_src, r_dst in ((y0, 0.0), (y1, float(y1 - y0))):
        for c_src, c_dst in ((x0, 0.0), (x1, float(x1 - x0))):
            vals = [cf[0] + cf[1] * r_src + cf[2] * c_src + cf[3] * r_src * c_src
                    for cf in coefs]
            out.append((r_dst, c_dst, vals[0], vals[1], vals[2]))
    return out


def _lstsq4(A, b):
    """Least-squares solve of a 4-column design matrix without numpy.

    Used to refit GCPs when auto-fit trims the crop (see _apply_auto_fit). Solves
    the 4x4 normal equations by Gaussian elimination with partial pivoting; that
    is plenty for 4 corner points and keeps this module free of a numpy import at
    module scope.
    """
    n = 4
    M = [[sum(A[k][i] * A[k][j] for k in range(len(A))) for j in range(n)]
         + [sum(A[k][i] * b[k] for k in range(len(A)))] for i in range(n)]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(M[r][i]))
        if abs(M[p][i]) < 1e-12:
            raise ValueError("singular GCP design matrix")
        M[i], M[p] = M[p], M[i]
        inv = 1.0 / M[i][i]
        for j in range(i, n + 1):
            M[i][j] *= inv
        for r in range(n):
            if r == i:
                continue
            f = M[r][i]
            if f:
                for j in range(i, n + 1):
                    M[r][j] -= f * M[i][j]
    return [M[i][n] for i in range(n)]


def detect_row_fit(im, rows_target=0, mode='width', buffer_frac=0.5, across=None):
    """Find the planted crop rows in a rectified plot image and return a pixel
    crop box that tightly encloses them. Uses the ROW/RIDGE periodicity
    (detrended-brightness oscillation) reinforced by greenness, so it works even
    before emergence and doesn't chase stray edge grass.

    Because the caller renders WIDER than the grid boundary, this can INCLUDE rows
    the grid clipped, then crop to exactly `rows_target` rows (or all detected).

    mode: 'width' (band width, keep original centre) | 'width_shift' (band as
    detected, recentre on rows) | 'width_length' (also trim along-row ends).
    `across` is the image axis running across the rows ('x' = rows vertical,
    'y' = rows horizontal); None decides from the crop shape.
    Returns dict(box=(x0,x1,y0,y1), rows, spacing_px, peaks, across) or None.
    """
    import numpy as np
    from PIL import Image
    if _across_row_axis(im, across) == 'y':
        # rows are horizontal bands here: detect on the transposed image (which
        # puts them upright) and map the box back, so the trimming happens
        # ACROSS the rows instead of along the plot's length
        r = detect_row_fit(im.transpose(Image.TRANSPOSE), rows_target=rows_target,
                           mode=mode, buffer_frac=buffer_frac, across='x')
        if r:
            bx0, bx1, by0, by1 = r['box']
            r['box'] = (by0, by1, bx0, bx1)
            r['across'] = 'y'
        return r
    W0, H0 = im.size
    ds = max(1, int(max(W0, H0) / 1800))
    a = np.asarray(im.resize((max(1, W0 // ds), max(1, H0 // ds))).convert('RGB')
                   ).astype('float32')
    H, W, _ = a.shape

    def sm(v, k):
        return np.convolve(v, np.ones(2 * k + 1) / (2 * k + 1), 'same') if k > 0 else v
    lum = a.mean(2)
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    exg = np.clip(2 * g - r - b, 0, None)
    # ignore uncovered pixels (BLACK or WHITE nodata) so column means reflect
    # only real imagery
    valid = ((lum > 12) & (a.min(2) < 245)).astype('float32')
    vcount = valid.sum(0)
    denom = np.maximum(vcount, 1.0)
    lum_col = (lum * valid).sum(0) / denom
    exg_col = (exg * valid).sum(0) / denom
    lum_col[vcount < 0.2 * H] = np.median(lum_col)   # mostly-black columns -> neutral
    trend = sm(lum_col, max(3, W // 15))
    osc = lum_col - trend
    # row spacing via FFT of the cross-row oscillation
    f = np.abs(np.fft.rfft((osc - osc.mean()) * np.hanning(W)))
    f[:max(2, W // 200)] = 0
    lo, hi = max(2, W // 120), max(3, W // 5)
    k = lo + int(np.argmax(f[lo:hi])) if hi > lo else 0
    sp = (W / k) if k > 0 else (W / 8.0)
    # row centres from greenness (fallback: ridge crests)
    sig = sm(exg_col, max(1, int(sp * 0.25)))
    if sig.max() < 1e-3:
        sig = sm(np.clip(osc, 0, None), max(1, int(sp * 0.25)))
    sign = sig / (sig.max() + 1e-6)
    sep = max(2, int(sp * 0.55))
    thr = max(0.12, float(np.percentile(sign, 55)))
    peaks = []
    i = 0
    while i < W:
        if sign[i] >= thr:
            j = min(W, i + sep)
            c = i + int(np.argmax(sign[i:j]))
            if not peaks or c - peaks[-1] >= sep:
                peaks.append(c)
            i = c + sep
        else:
            i += 1
    if len(peaks) < 2:
        return None
    if rows_target and rows_target < len(peaks):
        best = None
        for s in range(0, len(peaks) - rows_target + 1):
            grp = peaks[s:s + rows_target]
            score = abs((grp[0] + grp[-1]) / 2.0 - W / 2.0)
            if best is None or score < best[0]:
                best = (score, grp)
        grp = best[1]
    else:
        grp = peaks
    half = sp * buffer_frac
    bx0 = max(0.0, grp[0] - half)
    bx1 = min(float(W), grp[-1] + half)
    if mode == 'width':
        bw = bx1 - bx0
        cc = W / 2.0
        x0 = max(0.0, cc - bw / 2.0)
        x1 = min(float(W), cc + bw / 2.0)
    else:
        x0, x1 = bx0, bx1
    y0, y1 = 0, H
    if mode == 'width_length':
        ry = sm(exg.mean(1), max(2, H // 60))
        ryn = ry / (np.percentile(ry, 95) + 1e-6)
        yon = np.where(ryn > 0.20)[0]
        if len(yon) > H * 0.3:
            y0, y1 = int(yon[0]), int(yon[-1])
    return dict(box=(int(x0 * ds), int(x1 * ds), int(y0 * ds), int(y1 * ds)),
                rows=len(grp), spacing_px=sp * ds, across='x',
                peaks=[int(p * ds) for p in grp])


def read_psx_chunks(psx_path):
    """List chunk labels in a .psx by reading its project zips directly - no
    Metashape/licence needed. Returns (labels, active_label); ([],None) on error."""
    import zipfile
    try:
        base = psx_path[:-4] if psx_path.lower().endswith('.psx') else psx_path
        files_dir = base + '.files'
        with zipfile.ZipFile(os.path.join(files_dir, 'project.zip')) as z:
            docxml = z.read('doc.xml').decode('utf-8', 'replace')
        m = re.search(r'active_id="(\d+)"', docxml)
        active_id = m.group(1) if m else None
        entries = re.findall(r'<chunk id="(\d+)" path="([^"]+)"', docxml)
        labels, active_label = [], None
        for cid, cpath in entries:
            lbl = f"chunk {cid}"
            try:
                with zipfile.ZipFile(os.path.join(files_dir, cpath)) as cz:
                    cxml = cz.read('doc.xml').decode('utf-8', 'replace')
                mm = re.search(r'label="([^"]*)"', cxml)
                if mm and mm.group(1):
                    lbl = mm.group(1)
            except Exception:
                pass
            labels.append(lbl)
            if cid == active_id:
                active_label = lbl
        return labels, active_label
    except Exception:
        return [], None


def run_extraction(p, log, set_progress, cancel):
    """p: dict of parameters. log(str), set_progress(done,total), cancel: threading.Event"""
    import Metashape as ms
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None

    lic = ms.License()
    if not lic.valid:
        raise RuntimeError("No valid Metashape license found (license server unreachable?).")

    # --- use the GPU for Metashape compute (ortho raster export etc.) if present ---
    try:
        gpus = list(ms.app.enumGPUDevices())
        if gpus:
            def _gname(g):
                try:
                    return g['name'] if isinstance(g, dict) else getattr(g, 'name', str(g))
                except Exception:
                    return 'GPU'
            ms.app.gpu_mask = (1 << len(gpus)) - 1      # enable every device
            ms.app.cpu_enable = True                    # keep CPU helping alongside GPU
            log(f"GPU acceleration ON: {len(gpus)} device(s) - "
                + ", ".join(_gname(g) for g in gpus))
        else:
            log("GPU acceleration: no CUDA/OpenCL GPU detected - running on CPU")
    except Exception as e:
        log(f"GPU setup skipped ({e}) - running on CPU")

    plots, wkt, fields = read_plot_file(p['plots'])
    if not plots:
        raise RuntimeError("No polygons found in " + p['plots'])
    log(f"Loaded {len(plots)} plot polygons from {os.path.basename(p['plots'])}")

    wgs84 = ms.CoordinateSystem("EPSG::4326")
    override = (p.get('crs_override') or '').strip()
    if override and override.lower() != 'auto':
        spec = override if not override.isdigit() else f"EPSG::{override}"
        spec = re.sub(r'^EPSG:+', 'EPSG::', spec, flags=re.I)
        try:
            src_crs = ms.CoordinateSystem(spec)
        except Exception:
            raise RuntimeError(f"Could not parse CRS override '{override}'. "
                               "Use an EPSG code like 28350, or leave it on auto.")
        log(f"Plot CRS (user override): {src_crs.name}")
    elif wkt:
        try:
            src_crs = ms.CoordinateSystem(wkt)
        except Exception:
            raise RuntimeError("Could not parse the coordinate system stored in the plot "
                               "file. Set the 'Plot CRS' EPSG code manually.")
        log(f"Plot CRS (from file): {src_crs.name}")
    else:
        src_crs = wgs84
        xs = [abs(v) for pl in plots[:5] for v in pl['ring'][0]]
        if xs and max(xs) > 360:
            raise RuntimeError("The plot file carries no CRS information and its coordinates "
                               "do not look like lon/lat degrees. Set the 'Plot CRS' EPSG "
                               "code (e.g. 28350 for GDA2020 / MGA zone 50).")
        log("Plot CRS: assumed WGS84 lon/lat (file carries no CRS)")

    # plot ids
    id_field = p.get('id_field')
    if id_field in (None, '', AUTO_ID):
        id_field = pick_id_field(fields)
    for i, pl in enumerate(plots, 1):
        pl['id'] = safe_name(pl['props'].get(id_field, i)) if id_field else str(i)
    if id_field:
        log(f"Plot IDs from field '{id_field}'")
    else:
        log("No plot-id field found; using sequential IDs 1..N")

    only = p.get('plot_filter')            # optional: process only these plot ids
    if only:
        want = set(str(x) for x in only)
        plots = [pl for pl in plots if str(pl['id']) in want]
        log(f"Plot filter: {len(plots)} of the requested {len(want)} plot id(s)")
    limit = int(p.get('limit') or 0)
    if limit:
        plots = plots[:limit]
        log(f"Limit: first {limit} plots only (test run)")

    log("Opening Metashape project ...")
    doc = ms.Document()
    doc.open(p['psx'], read_only=True)
    chunks = [c for c in doc.chunks if c is not None]
    if not chunks:
        raise RuntimeError("Project has no chunks.")
    log("Project chunks: " + ", ".join(f"'{c.label}'" for c in chunks))
    want = (p.get('chunk') or '').strip()
    ch = None
    if want and want.lower() not in ('active', '(active)', 'auto', ''):
        for c in chunks:                                   # match by label first
            if (c.label or '').strip() == want:
                ch = c
                break
        if ch is None:                                     # then by index
            try:
                i = int(want)
                if 0 <= i < len(chunks):
                    ch = chunks[i]
            except ValueError:
                pass
        if ch is None:
            raise RuntimeError(
                f"Chunk '{want}' not found. Available: "
                + ", ".join(f"'{c.label}'" for c in chunks))
    if ch is None:
        ch = doc.chunk or chunks[0]                        # fall back to active
    if ch is None:
        raise RuntimeError("Project has no chunk.")
    log(f"Using chunk: '{ch.label}'  ({len(ch.cameras)} cameras)")
    if ch.crs is None or ch.transform.matrix is None:
        raise RuntimeError("Project chunk is not georeferenced (no CRS / transform).")
    T = ch.transform.matrix
    Tinv = T.inv()
    dem = ch.elevation
    cams = [c for c in ch.cameras if c.transform]
    log(f"Chunk '{ch.label}': {len(cams)} aligned cameras, CRS: {ch.crs.name}, "
        f"DEM: {'yes' if dem else 'NO'}, orthomosaic: {'yes' if ch.orthomosaic else 'no'}")
    if not cams:
        raise RuntimeError("No aligned cameras in the project.")

    def tr(v, a, b):
        return ms.CoordinateSystem.transform(v, a, b)

    # fallback ground height = chunk region centre height (chunk CRS)
    centre_chunk = ch.crs.project(T.mulp(ch.region.center))
    z_fallback = centre_chunk.z

    # sanity check: does the first plot actually sit near the project?
    p0 = plots[0]['ring'][0]
    p0w = tr(ms.Vector([p0[0], p0[1], 0]), src_crs, wgs84)
    c0w = tr(ms.Vector([centre_chunk.x, centre_chunk.y, centre_chunk.z]), ch.crs, wgs84)
    dist_km = math.hypot((p0w.x - c0w.x) * 111.32 * math.cos(math.radians(c0w.y)),
                         (p0w.y - c0w.y) * 110.54)
    if dist_km > 10:
        log(f"[!] WARNING: plots are ~{dist_km:.0f} km away from the project centre - "
            "wrong CRS or wrong project? Continuing anyway ...")

    def ground_z(pt_chunk_2d):
        if dem is not None:
            z = dem.altitude(pt_chunk_2d)
            if z is not None and not math.isnan(z):
                return z
        return z_fallback

    margin = float(p.get('margin') or 0.0)

    def expand_wgs(ring, m):
        """Grow the plot ring outward by m metres on ALL sides (uniform perpendicular
        buffer), so every edge — including the long sides where the outer plant rows
        sit — gains m metres of extra ground. Done in a local metric frame, then
        converted back to WGS84 degrees.

        (The old method pushed each corner radially from the centroid, which for a
        long thin plot is just a scale: it added ~m at the short ends but almost
        nothing on the long sides, clipping the edge rows. A true buffer fixes that.)"""
        if m <= 0:
            return ring
        cx = sum(q[0] for q in ring) / len(ring)
        cy = sum(q[1] for q in ring) / len(ring)
        mlon = 111320 * math.cos(math.radians(cy))
        mlat = 110540
        # ring -> local metres relative to centroid
        pts_m = [((lon - cx) * mlon, (lat - cy) * mlat) for lon, lat in ring]
        try:
            from shapely.geometry import Polygon
            poly = Polygon(pts_m)
            if not poly.is_valid:
                poly = poly.buffer(0)
            grown = poly.buffer(m, join_style='mitre', mitre_limit=10.0)
            ext = list(grown.exterior.coords)[:-1]          # drop closing duplicate
            if len(ring) == 4 and len(ext) == 4:            # keep A,B,C,D order/count
                # re-order buffered corners to match each original corner (nearest)
                ext = [min(ext, key=lambda e: (e[0]-px)**2 + (e[1]-py)**2)
                       for px, py in pts_m]
            return [[cx + x / mlon, cy + y / mlat] for x, y in ext]
        except Exception:
            # fallback: radial push from centroid (legacy behaviour)
            out = []
            for lon, lat in ring:
                dx = (lon - cx) * mlon
                dy = (lat - cy) * mlat
                d = math.hypot(dx, dy) or 1.0
                out.append([lon + dx / d * m / mlon, lat + dy / d * m / mlat])
            return out

    # raw image fallback index (basename and stem, lowercase)
    raw_index = {}
    if p.get('rawdir') and os.path.isdir(p['rawdir']):
        n = 0
        for root, _dirs, files in os.walk(p['rawdir']):
            for fn in files:
                if os.path.splitext(fn)[1].lower() in RAW_EXTS:
                    full = os.path.join(root, fn)
                    raw_index.setdefault(fn.lower(), full)
                    raw_index.setdefault(os.path.splitext(fn)[0].lower(), full)
                    n += 1
        log(f"Raw image folder indexed: {n} image files "
            "(preferred over the project's own photo paths)")

    def resolve_photo(cam):
        # the user-supplied raw folder wins over the project's own photo path,
        # so extraction always uses the rendition the user considers original
        path = cam.photo.path
        base = os.path.basename(path).lower()
        stem = os.path.splitext(base)[0]
        for key in (base, stem, cam.label.lower()):
            if key in raw_index:
                return raw_index[key]
        for e in RAW_EXTS:
            key = (cam.label + e).lower()
            if key in raw_index:
                return raw_index[key]
        if os.path.exists(path):
            return path
        return None

    do_raw = bool(p.get('do_raw', True))
    do_ortho = bool(p.get('do_ortho', False))
    raw_mode = (p.get('raw_mode') or 'auto').lower()   # auto|bbox|mask|rectify|native|best
    smart_seam = bool(p.get('smart_seam', True))       # min-error seam + feathering
    GAPFILL_N = max(2, int(p.get('gapfill_frames') or 8))   # candidates offered to fill
    mesh_step = max(16, int(p.get('mesh_step') or MESH_STEP_PX))
    seam_map = bool(p.get('seam_map', False))   # QC: dump a per-pixel source map
    seam_slack = max(8, int(p.get('seam_slack_px') or 300))  # how far a seam may
    #        retreat into the overlap to step around a plant (px of the output crop)
    CAND_REPORT_N = 12                                 # rows per plot in candidates.csv
    want_cand_csv = bool(p.get('candidates_csv', True))
    candidates_csv = [] if (do_raw and want_cand_csv) else None
    _cand_lock = threading.Lock()
    to_srgb = bool(p.get('to_srgb', False))            # convert raw crops to sRGB
    auto_fit = bool(p.get('auto_fit', False))          # snap crop to detected rows
    fit_rows = int(p.get('rows_per_plot') or 0)        # target rows/plot (0=all found)
    fit_mode = (p.get('fit_mode') or 'width').lower().replace(' ', '_')
    fit_probe_m = float(p.get('fit_probe_m') or 0.35)  # extra width rendered each side
    if do_raw and auto_fit:
        log(f"Auto-fit to rows: ON  (mode={fit_mode}, "
            f"{'target ' + str(fit_rows) + ' rows' if fit_rows else 'all detected rows'}, "
            f"probe ±{fit_probe_m:.2f} m). Rows clipped by the grid are recovered by "
            f"rendering wider, then the crop is snapped to the rows.")
        if raw_mode not in ('rectify', 'auto', 'best'):
            log("[!] Auto-fit needs an upright crop - forcing raw style to 'Rectified'.")
            raw_mode = 'rectify'
    if do_raw and to_srgb:
        log("Colour: converting raw crops from their embedded ICC profile "
            "(e.g. Adobe RGB) to sRGB")
    if do_raw:
        log({'native': "Raw crop style: native bounding box (no resample, no black; "
                       "full plot in frame, bit-identical pixels)",
             'auto': "Raw crop style: auto-straighten (bit-exact if upright, "
                     "rectify only plots tilted >3 deg)",
             'bbox': "Raw crop style: bounding box (whole tilted plot + surrounds)",
             'mask': "Raw crop style: masked to plot polygon (outside set to black)",
             'rectify': "Raw crop style: rectified to an upright plot rectangle",
             'best': "Raw crop style: BEST single frame per plot - the one image "
                     "covering most of the plot, rectified, NO stitching"}
            .get(raw_mode, f"Raw crop style: {raw_mode}"))
        log("Rectification: exact projected mesh (full camera model incl. lens "
            f"distortion, ~{mesh_step}px cells) - all frames share one ground "
            "grid, so they register to each other")
        if raw_mode in ('rectify', 'auto'):
            log("Stitching: " + ("minimum-error seam through the overlap + "
                                 "feathered radiometric match (no pixel blending)"
                                 if smart_seam else
                                 "plain gap fill at the frame footprint boundary"))
    if do_ortho and ch.orthomosaic is None:
        log("[!] Orthomosaic crops requested but the project has no orthomosaic - skipping those.")
        do_ortho = False
    raw_dir = os.path.join(p['outdir'], "RawCrops")
    ortho_dir = os.path.join(p['outdir'], "OrthoCrops")
    if do_raw:
        os.makedirs(raw_dir, exist_ok=True)
    if do_ortho:
        os.makedirs(ortho_dir, exist_ok=True)
    quality = int(p.get('quality') or 92)
    n_per_plot = max(1, int(p.get('n_per_plot') or 1))
    fmt = (p.get('format') or 'JPEG').upper()
    geotag = fmt.startswith('GEO')
    # CRS for the GeoTIFF ground control points. Coordinates are transformed into
    # it, not relabelled; blank keeps the internal EPSG:4326.
    try:
        _GCP_EPSG[0] = int(str(p.get('gcp_epsg') or '').strip() or 0) or None
    except (TypeError, ValueError):
        _GCP_EPSG[0] = None
    if geotag and _GCP_EPSG[0]:
        log(f"GCP CRS: writing ground control points in EPSG:{_GCP_EPSG[0]}")
    if geotag:
        ext, save_kw = '.tif', dict(compression='tiff_lzw')        # lossless GeoTIFF
        import rasterio                                            # fail early if missing
        log("Output format: GeoTIFF (lossless, GCP-geotagged, source EXIF copied)")
    elif fmt.startswith('PNG'):
        ext, save_kw = '.png', {}                                  # lossless
    elif fmt.startswith('TIF'):
        ext, save_kw = '.tif', dict(compression='tiff_lzw')        # lossless
    else:
        ext, save_kw = '.jpg', dict(quality=quality)
    if not geotag:
        log(f"Output format: {ext[1:].upper()}"
            + (f" quality {quality} (lossy)" if ext == '.jpg' else " (lossless)"))

    # local units-per-metre in chunk CRS (for ortho export), measured at first plot
    def units_per_metre(lon, lat):
        dlon = 1.0 / (111320 * math.cos(math.radians(lat)))
        dlat = 1.0 / 110540
        p0 = tr(ms.Vector([lon, lat, 0]), wgs84, ch.crs)
        pe = tr(ms.Vector([lon + dlon, lat, 0]), wgs84, ch.crs)
        pn = tr(ms.Vector([lon, lat + dlat, 0]), wgs84, ch.crs)
        return math.hypot(pe.x - p0.x, pe.y - p0.y), math.hypot(pn.x - p0.x, pn.y - p0.y)

    # one-image cache: consecutive plots often crop from the same raw photo
    # how many CPU cores to use for the image work (decode/crop/convert/save)
    workers = int(p.get('workers') or 0)
    if workers <= 0:
        # each in-flight frame is a ~370 MB decoded image, so cap the auto
        # default for a safe memory footprint; the user can raise it in the GUI.
        workers = min(8, (os.cpu_count() or 4))
    workers = max(1, workers)

    # thread-safe LRU cache of decoded source frames + their ICC profile.
    # Several plots often crop from the same frame, so each frame decodes once;
    # decoding happens OUTSIDE the lock so different frames decode in parallel.
    from collections import OrderedDict
    _img_lru = OrderedDict()          # path -> (RGB image, icc_profile_bytes)
    _img_lock = threading.Lock()
    # Cap defaults to workers+1, which assumes one frame in flight per worker.
    # That holds for plain crops, but a rectified crop needs every frame the plot
    # touches, and auto-fit widens the render by 2*fit_probe_m so it touches more
    # of them. workers*frames_per_plot then exceeds the cap, frames get evicted
    # and re-read from disk, and throughput collapses (measured: 2.5 MB/s against
    # 75.8 MB/s single-threaded on the same array). Override with PE_FRAME_CACHE
    # when the source frames are large: each cached frame costs w*h*3 bytes
    # (~366 MB for a 122 MP iXM-GS120 frame).
    try:
        _cache_cap = max(2, int(os.environ.get('PE_FRAME_CACHE', workers + 1)))
    except (TypeError, ValueError):
        _cache_cap = max(2, workers + 1)

    def load_image(path):
        """Return (decoded RGB image, icc_profile_bytes). Thread-safe."""
        with _img_lock:
            hit = _img_lru.get(path)
            if hit is not None:
                _img_lru.move_to_end(path)
                return hit
        # Slurp the whole frame in ONE sequential read, then decode from memory.
        # PIL streams a JPEG in small blocks as it decodes; with several worker
        # threads interleaving those blocks across different files, a USB/spinning
        # array degenerates into a seek storm and nothing finishes. Measured on the
        # QNAP TR-004 that holds these frames: 8 threads reading whole files reach
        # 109 MB/s, the same 8 threads letting PIL stream reach 2-5 MB/s. The blob
        # is transient (one per in-flight frame, ~150 MB) and is freed once decoded.
        import io as _io
        with open(path, 'rb', buffering=0) as _fh:
            _blob = _fh.read()
        im = Image.open(_io.BytesIO(_blob))
        icc = im.info.get('icc_profile')          # grab before convert
        if im.mode != 'RGB':
            im = im.convert('RGB')
        im.load()                                 # force full decode -> crop() is read-only
        with _img_lock:
            _img_lru[path] = (im, icc)
            _img_lru.move_to_end(path)
            while len(_img_lru) > _cache_cap:
                _img_lru.popitem(last=False)
        return (im, icc)

    _srgb_dst = [None]                            # sRGB target profile, built once
    _srgb_lock = threading.Lock()

    def to_srgb_crop(im, icc):
        """Convert a crop from its embedded ICC profile (e.g. Adobe RGB) to sRGB.
        No-op if conversion is off or the source has no profile. Thread-safe."""
        if not to_srgb or not icc:
            return im
        try:
            import io
            from PIL import ImageCms
            with _srgb_lock:
                if _srgb_dst[0] is None:
                    _srgb_dst[0] = ImageCms.createProfile('sRGB')
            src = ImageCms.ImageCmsProfile(io.BytesIO(icc))
            return ImageCms.profileToProfile(im, src, _srgb_dst[0], outputMode='RGB')
        except Exception as e:
            log(f"[!] sRGB conversion failed: {e}")
            return im

    exif_cache = {}
    _exif_lock = threading.Lock()

    def exif_tags(path):
        """Source image EXIF (incl. Exif/GPS sub-IFDs) as flat str->str dict.
        Thread-safe; keeps one entry per unique frame."""
        with _exif_lock:
            cached = exif_cache.get(path)
        if cached is not None:
            return cached
        from PIL import ExifTags
        tags = {}
        try:
            ex = Image.open(path).getexif()
            ifds = [(ex, ExifTags.TAGS)]
            for ifd_id, names in ((ExifTags.IFD.Exif, ExifTags.TAGS),
                                  (ExifTags.IFD.GPSInfo, ExifTags.GPSTAGS)):
                try:
                    ifds.append((ex.get_ifd(ifd_id), names))
                except Exception:
                    pass
            for ifd, names in ifds:
                for k, v in ifd.items():
                    if isinstance(v, bytes) or isinstance(v, dict):
                        continue
                    name = names.get(k, f"Tag{k}")
                    tags.setdefault(f"EXIF_{name}", str(v)[:500])
        except Exception:
            pass
        tags['SOURCE_IMAGE'] = path
        with _exif_lock:
            exif_cache[path] = tags
        return tags

    _white = {}
    _white_lock = threading.Lock()

    def _white_full(size):
        """Cached all-white L image the size of a source frame (read-only)."""
        with _white_lock:
            im = _white.get(size)
            if im is None:
                im = Image.new('L', size, 255)
                _white[size] = im
            return im

    RENDER_MAX_ANISO = 1.6            # frames more oblique than this smear when warped

    def _rectify(job, layer, W, H):
        """Warp ONE source frame onto the plot rectangle.

        The mapping comes from the real camera model: a grid of ground points
        spanning the plot was projected with cam.project() (full perspective +
        lens distortion) in phase A, and PIL replays it as a piecewise mesh.
        The old code fitted a single bilinear QUAD through only the plot's 4
        projected corners, which on this data is wrong by 20-95 px in the middle
        of a plot - so each frame bent the rows differently and two frames could
        not agree at a stitch boundary. Sampling every frame at the SAME ground
        nodes removes that error and makes the frames register to each other.
        Returns (rgb uint8 HxWx3, valid mask HxW)."""
        import numpy as np
        im, icc = load_image(layer['path'])
        nodes = layer.get('nodes')
        if nodes:
            nx, ny = layer['nx'], layer['ny']
            data = build_pil_mesh(nodes, nx, ny, W, H)
            if not data:
                return (np.zeros((H, W, 3), np.uint8), np.zeros((H, W), bool))
            rect = im.transform((W, H), Image.MESH, data, Image.BICUBIC)
            vm = _white_full(im.size).transform((W, H), Image.MESH, data,
                                                Image.NEAREST)
        else:                                   # fallback: legacy 4-corner quad
            A, B, C, D = layer['corners']
            quad = [A[0], A[1], D[0], D[1], C[0], C[1], B[0], B[1]]
            rect = im.transform((W, H), Image.QUAD, quad, Image.BICUBIC)
            vm = _white_full(im.size).transform((W, H), Image.QUAD, quad,
                                                Image.NEAREST)
        rect = to_srgb_crop(rect, icc)
        return (np.array(rect, dtype=np.uint8),      # writable copy
                np.asarray(vm) > 127)

    def _render_fill(job):
        """Rectify mode with gap-fill. Every candidate frame is warped onto the
        SAME ground grid (see _rectify), so they overlay pixel-for-pixel; the
        near-nadir frames do the bulk and oblique frames only patch what is left.

        Each output pixel is native raw from EXACTLY ONE frame (one resample) -
        nothing is averaged, so a plant can never be blended into a ghost. Two
        further steps kill the stitch artefacts:
          * the hand-over boundary is routed by a minimum-error cut through the
            overlap, so it runs over bare soil rather than slicing a seedling
            (a plant sitting on the old footprint edge used to appear twice,
            because its top is displaced differently in the two views);
          * the leftover brightness/colour step is removed with a per-column,
            heavily smoothed offset that fades out away from the seam - a
            low-frequency correction only, the texture stays 100% single-frame.
        """
        import numpy as np
        W, H = job['wpx'], job['hpx']
        layers = job['layers']
        good = [i for i, l in enumerate(layers)
                if l.get('aniso', 99.0) <= RENDER_MAX_ANISO]
        rest = [i for i, l in enumerate(layers)
                if l.get('aniso', 99.0) > RENDER_MAX_ANISO]
        order = good + rest              # good (clean) frames first, oblique last
        out = None                       # uint8 HxWx3 accumulator
        filled = None                    # bool HxW: pixel already has data
        srcid = None                     # uint8 HxW: which layer filled each pixel
        used = 0
        seams = []
        for pos, i in enumerate(order):
            if pos >= len(good) and filled is not None and filled.all():
                break                    # holes closed -> stop before oblique frames
            layer = layers[i]
            rgb, vmask = _rectify(job, layer, W, H)
            if out is None:
                out = rgb
                filled = vmask.copy()
                srcid = np.where(vmask, np.uint8(used), np.uint8(255))
                used += 1
                if filled.all():
                    break
                continue
            need = vmask & (~filled)                    # still-empty pixels
            if not need.any():
                continue
            ra = rgb.astype(np.float32)
            # radiometric match: make this frame match the already-placed pixels
            # in the overlap. A per-channel AFFINE match (mean AND contrast/std),
            # so a seam is removed even when the two frames differ in view angle
            # (oblique frames have deeper furrow shadows -> higher contrast). The
            # std ratio is clamped so noise isn't amplified / detail not crushed.
            ov = vmask & filled
            if ov.sum() > 500:
                for c in range(3):
                    oc = out[..., c][ov].astype(np.float32); rc = ra[..., c][ov]
                    m_out, m_rect = float(oc.mean()), float(rc.mean())
                    s_out, s_rect = float(oc.std()), float(rc.std())
                    s_ratio = min(1.6, max(0.6, s_out / s_rect)) if s_rect > 1e-3 else 1.0
                    ra[..., c] = (ra[..., c] - m_rect) * s_ratio + m_out
                np.clip(ra, 0, 255, out=ra)
            take = need
            if job.get('smart_seam', True):
                try:
                    # green channel only: enough to score how alike two frames
                    # look, and it avoids two full-size float temporaries
                    seam = min_error_seam(np, out[..., 1], ra[..., 1],
                                          filled, vmask,
                                          slack=job.get('seam_slack', 300))
                except Exception as e:
                    seam = None
                    log(f"plot {job['pid']}: seam search failed ({e}); "
                        "falling back to plain gap fill")
                if seam is not None:
                    take = seam['take']
                    seams.append(seam)
                    try:
                        ra = feather_seam(np, out, ra, seam)
                    except Exception as e:
                        log(f"plot {job['pid']}: seam feather failed ({e})")
            out[take] = ra[take].astype(np.uint8)
            srcid[take] = np.uint8(used)
            filled |= vmask
            used += 1
            if filled.all():
                break
        if out is None:
            out = np.zeros((H, W, 3), np.uint8)
            filled = np.zeros((H, W), bool)
            srcid = np.full((H, W), 255, np.uint8)
        job['filled_pct'] = 100.0 * float(filled.mean())
        job['layers_used'] = used
        # seam location, recorded so a reviewer knows where to spot-check. These
        # long plots ALWAYS need 2 frames along their length (one PhaseOne frame
        # covers ~4.2 m of a 4.35 m plot), so a seam is normal; what matters is
        # that it no longer cuts through plants.
        job['seam_y_pct'] = ''
        if used > 1:
            chg = np.zeros((H, W), bool)
            chg[1:, :] = (srcid[1:, :] != srcid[:-1, :]) & filled[1:, :] & filled[:-1, :]
            rows = np.where(chg.any(1))[0]
            if len(rows):
                job['seam_y_pct'] = f"{100.0 * float(np.median(rows)) / H:.0f}"
            else:
                chg = np.zeros((H, W), bool)
                chg[:, 1:] = (srcid[:, 1:] != srcid[:, :-1]) & filled[:, 1:] & filled[:, :-1]
                cols = np.where(chg.any(0))[0]
                if len(cols):
                    job['seam_y_pct'] = f"x{100.0 * float(np.median(cols)) / W:.0f}"
        job['seam_smart'] = len(seams)
        if job.get('seam_map') and used > 1:
            # QC aid: which frame each pixel came from (0,1,2..., 255 = no data)
            try:
                sp = os.path.splitext(job['out'])[0] + "_seam.png"
                sd = os.path.dirname(sp)
                os.makedirs(sd, exist_ok=True)
                vis = np.where(srcid == 255, np.uint8(0),
                               (60 + srcid.astype(np.uint16) * 70) % 256).astype(np.uint8)
                Image.fromarray(vis, 'L').save(sp)
            except Exception as e:
                log(f"plot {job['pid']}: seam map not written ({e})")
        return Image.fromarray(out, 'RGB')

    def _unprobe(crop, job, why):
        """Auto-fit could not fit: trim the probe margin so what ships is the PLOT
        rectangle, not the widened render.

        Returning the widened crop - which is what this did originally - ships a crop
        2*fit_probe_m wider than the plot, containing the NEIGHBOURING plots' rows.
        That is worse than never enabling auto-fit. Measured on OZ Barley York: the
        4 plots with too little emergence for rows to be detected came out 1.82 m
        across against a 1.12 m plot, and the extra 0.70 m held a neighbour's rows
        on each side.
        """
        job['fit_note'] = why
        pm = float(job.get('fit_probe_used') or 0.0)
        if pm <= 0:
            return crop
        W, H = crop.size
        # rows run along the plot's long axis, so 'across rows' is the shorter side
        if W <= H:
            a_px, mpp, axis = W, float(job.get('mpp_x') or 0.0), 'x'
        else:
            a_px, mpp, axis = H, float(job.get('mpp_y') or 0.0), 'y'
        across_m = a_px * mpp
        if mpp <= 0 or across_m <= 2.5 * pm:
            return crop                       # cannot tell probe from plot: leave it
        cut = int(round(pm / across_m * a_px))
        if cut <= 0 or 2 * cut >= a_px:
            return crop
        box = (cut, 0, W - cut, H) if axis == 'x' else (0, cut, W, H - cut)
        out = crop.crop(box)
        job['fit_note'] = f"{why}; probe trimmed back to the plot rectangle"
        if job.get('geotag') and job.get('gcps'):
            job['gcps'] = _refit_gcps_box(job['gcps'], box[0], box[1], box[2], box[3])
        return out

    def _apply_auto_fit(crop, job):
        """Detect the planted rows in the (wider-than-grid) rectified crop and
        crop back to exactly the rows. Records fitted width/length/area (m) and
        the row count, and shifts geotag GCPs into the new frame. On any doubt it
        keeps the full crop (never cuts rows), and REFUSES to fit a crop with
        significant black (incomplete coverage) since the area would be wrong."""
        import numpy as np
        # black fraction: incomplete-coverage crops give garbage row detection
        sc = max(crop.width, crop.height) / 1200.0
        small = (crop.resize((max(1, int(crop.width / sc)), max(1, int(crop.height / sc))))
                 if sc > 1 else crop)
        arr = np.asarray(small.convert('RGB'))
        black_pct = 100.0 * float(np.mean(arr.sum(2) < 24))
        if black_pct > 10.0:
            return _unprobe(crop, job,
                            f'low coverage ({black_pct:.0f}% black) - area NOT fitted')
        try:
            fb = detect_row_fit(crop, rows_target=job.get('fit_rows', 0),
                                 mode=job.get('fit_mode', 'width'),
                                 across=job.get('fit_across'))
        except Exception as e:
            log(f"plot {job['pid']}: auto-fit detect failed ({e}); keeping full crop")
            fb = None
        if not fb:
            return _unprobe(crop, job, 'rows not detected')
        if job.get('fit_rows') and fb['rows'] < job['fit_rows']:
            return _unprobe(crop, job,
                            f"only {fb['rows']}/{job['fit_rows']} rows found")
        x0, x1, y0, y1 = fb['box']
        W, H = crop.size
        x0 = max(0, min(x0, W - 2)); x1 = max(x0 + 1, min(x1, W))
        y0 = max(0, min(y0, H - 2)); y1 = max(y0 + 1, min(y1, H))
        crop = crop.crop((x0, y0, x1, y1))
        if job['geotag'] and job.get('gcps'):
            # REFIT the GCPs onto the trimmed frame; do NOT filter to the ones that
            # happen to fall inside it. For a rectified crop the GCPs sit at the
            # corners of the rendered rectangle (rows 0 and H), and auto-fit cuts
            # the top and bottom away - so filtering drops ALL of them and the
            # GeoTIFF ships with no georeferencing at all. That is exactly what
            # happened on the 2026-08-12 autofit run: 691 of 696 crops came out
            # with 0 GCPs and save_geotiff logged "no GCPs landed inside the crop"
            # for every one. Instead fit v = a + b*row + c*col + d*row*col (exactly
            # determined by 4 corner GCPs, which is what these are) for lon, lat
            # and z, then evaluate it at the corners of the new frame.
            g = job['gcps']
            ng = []
            if len(g) >= 4:
                try:
                    A = [[1.0, gg[0], gg[1], gg[0] * gg[1]] for gg in g]
                    coefs = []
                    for k in (2, 3, 4):                      # lon, lat, z
                        b = [gg[k] for gg in g]
                        coefs.append(_lstsq4(A, b))
                    for r_src, r_dst in ((y0, 0.0), (y1, float(y1 - y0))):
                        for c_src, c_dst in ((x0, 0.0), (x1, float(x1 - x0))):
                            vals = [cf[0] + cf[1] * r_src + cf[2] * c_src
                                    + cf[3] * r_src * c_src for cf in coefs]
                            ng.append((r_dst, c_dst, vals[0], vals[1], vals[2]))
                except Exception as e:
                    log(f"plot {job['pid']}: GCP refit failed ({e}); shifting instead")
                    ng = []
            if not ng:                                       # <4 GCPs, or refit failed
                ng = [(r - y0, c - x0, lon, lat, z) for (r, c, lon, lat, z) in g]
            job['gcps'] = ng
        mppx, mppy = job.get('mpp_x', 0.0), job.get('mpp_y', 0.0)
        job['fit_rows_found'] = fb['rows']
        # width = ACROSS the rows, length = ALONG them, whichever image axis
        # each happens to be in this trial's shapefile vertex order
        if fb.get('across', 'x') == 'y':
            job['fit_w_m'] = (y1 - y0) * mppy if mppy else 0.0
            job['fit_l_m'] = (x1 - x0) * mppx if mppx else 0.0
        else:
            job['fit_w_m'] = (x1 - x0) * mppx if mppx else 0.0
            job['fit_l_m'] = (y1 - y0) * mppy if mppy else 0.0
        job['fit_area_m2'] = job['fit_w_m'] * job['fit_l_m']
        job['cw'], job['ch'] = crop.size
        return crop

    def _render(job):
        """Do the pixel work for one raw crop (runs in a worker thread):
        decode/crop/mask/rectify/sRGB/save. Returns (job, error-or-None)."""
        try:
            if job.get('layers'):
                crop = _render_fill(job)
            elif job['eff'] == 'rectify' and job.get('self_layer'):
                # single-frame rectify: same exact projected-mesh warp as the
                # gap-fill path, so the plot's internal geometry is right (rows
                # stay straight) instead of bent by a 4-corner bilinear fit
                import numpy as np
                rgb, vmask = _rectify(job, job['self_layer'], job['wpx'], job['hpx'])
                job['filled_pct'] = 100.0 * float(vmask.mean())
                job['layers_used'] = 1
                job['seam_y_pct'] = ''
                crop = Image.fromarray(rgb, 'RGB')
            else:
                im, icc = load_image(job['path'])
                crop = im.crop((job['x0'], job['y0'], job['x1'], job['y1']))
                eff, plot_uv = job['eff'], job['plot_uv']
                if eff == 'mask' and len(plot_uv) >= 3:
                    from PIL import ImageDraw
                    mk = Image.new('L', crop.size, 0)
                    ImageDraw.Draw(mk).polygon([tuple(q) for q in plot_uv], fill=255)
                    crop = Image.composite(crop, Image.new('RGB', crop.size, (0, 0, 0)), mk)
                elif eff == 'rectify' and len(plot_uv) == 4:
                    A, B, C, D = plot_uv
                    quad = [A[0], A[1], D[0], D[1], C[0], C[1], B[0], B[1]]
                    crop = crop.transform((job['wpx'], job['hpx']),
                                          Image.QUAD, quad, Image.BICUBIC)
                crop = to_srgb_crop(crop, icc)
            if job.get('auto_fit'):
                crop = _apply_auto_fit(crop, job)
            if job['geotag']:
                save_geotiff(job['out'], crop, job['gcps'], exif_tags(job['path']), log)
            else:
                crop.save(job['out'], **save_kw)
            return (job, None)
        except Exception as e:
            return (job, str(e))

    manifest = []
    ok_raw = ok_ortho = full_raw = 0
    total = len(plots)

    # Raw crops render on a live thread pool: each job is submitted the moment
    # its geometry is ready, so the heavy image work runs in parallel WHILE the
    # main thread moves on (incl. the sequential ortho export). Crops therefore
    # appear in RawCrops\ progressively during the run, not all at the end.
    from concurrent.futures import ThreadPoolExecutor
    raw_pool = ThreadPoolExecutor(max_workers=workers) if do_raw else None
    raw_futures = []
    _rawcount = {'ok': 0, 'full': 0, 'done': 0}
    _rawcount_lock = threading.Lock()
    _fallback = []            # plots whose raw was too incomplete -> use ortho

    def _on_raw_done(fut):
        # runs in the worker thread as each crop finishes -> live progress log
        try:
            job, err = fut.result()
        except Exception as e:
            log(f"raw render callback error: {e}")
            return
        fp = job.get('filled_pct')
        if not err and job['is_rank1'] and fp is not None:
            # gap-filled crop: report the real post-fill coverage in the manifest
            try:
                area = float(job['row'].get('plot_area_m2') or 0)
            except ValueError:
                area = 0.0
            job['row']['raw_covered_pct'] = f"{fp:.1f}"
            job['row']['raw_covered_m2'] = f"{area * fp / 100:.3f}"
            job['cov'] = fp
            job['full'] = fp >= 99.5
            job['flag'] = "" if fp >= 99.5 else (
                f"  [!] {fp:.1f}% of the plot in this frame (rest is black)"
                if job.get('layers_used', 1) == 1 else
                f"  [!] {fp:.0f}% covered even after fill")
            nfr = job.get('layers_used', 1)
            job['cam_label'] = (
                f"{nfr} frames (stitched"
                + (", min-error seam" if job.get('seam_smart') else ", gap-fill")
                + ")" if nfr > 1 else f"{job['cam_label']} (single frame)")
        seam_msg = ''
        if not err and job['is_rank1']:
            nfr = job.get('layers_used', 1)
            job['row']['frames_used'] = str(nfr)
            job['row']['seam_y_pct'] = job.get('seam_y_pct', '') if nfr > 1 else ''
            if nfr > 1 and job.get('seam_y_pct'):
                seam_msg = f"  seam@{job['seam_y_pct']}%"
        fit_msg = ''
        if not err and job['is_rank1'] and job.get('fit_rows_found') is not None:
            job['row']['fit_rows'] = str(job['fit_rows_found'])
            job['row']['fit_width_m'] = f"{job.get('fit_w_m', 0):.3f}"
            job['row']['fit_length_m'] = f"{job.get('fit_l_m', 0):.3f}"
            job['row']['fit_area_m2'] = f"{job.get('fit_area_m2', 0):.3f}"
            fit_msg = (f"  fit {job['fit_rows_found']} rows "
                       f"{job.get('fit_w_m', 0):.2f}×{job.get('fit_l_m', 0):.2f}m "
                       f"={job.get('fit_area_m2', 0):.2f}m²")
        elif not err and job['is_rank1'] and job.get('fit_note'):
            job['row']['note'] = ((job['row'].get('note') or '') +
                                  f" auto-fit: {job['fit_note']}").strip()
            fit_msg = f"  [auto-fit: {job['fit_note']}]"
            # too much black to fit -> queue an ortho fallback (if ortho is on)
            if do_ortho and 'low coverage' in job['fit_note']:
                with _rawcount_lock:
                    _fallback.append((job['pid'], job['row'], job['out']))
        with _rawcount_lock:
            _rawcount['done'] += 1
            if not err and job['is_rank1']:
                _rawcount['ok'] += 1
                if job['full']:
                    _rawcount['full'] += 1
        if err:
            if job['is_rank1']:
                job['row']['note'] = (
                    (job['row'].get('note') or '') + f' render error: {err}').strip()
            log(f"plot {job['pid']}: raw render ERROR {err}")
        else:
            log(f"  raw  plot {job['pid']} <- {job['cam_label']}  "
                f"{job['cw']}x{job['ch']}px  cover {job['cov']:.0f}%"
                f"{job['flag']}{fit_msg}{seam_msg}")

    if do_raw:
        log(f"Raw crops render in parallel across {workers} core(s) "
            f"as each plot's geometry is ready.")
    for idx, pl in enumerate(plots):
        if cancel.is_set():
            log("Cancelled by user.")
            break
        pid = pl['id']
        row = dict(plot_id=pid, plot_area_m2='', margin_m=f"{margin:.2f}",
                   raw_image='', raw_file='', crop_w='', crop_h='',
                   corners_in_frame='', raw_covered_pct='', raw_covered_m2='',
                   cand_frames='', best_frame='', best_cover_pct='', best_cover_m2='',
                   frames_used='', seam_y_pct='',
                   fit_rows='', fit_width_m='', fit_length_m='', fit_area_m2='',
                   ortho_file='', ortho_covered_pct='', ortho_covered_m2='',
                   source_path='', note='')
        try:
            # ring -> WGS84; keep both the bare plot ring (coverage scoring) and
            # the margin-expanded ring (crop extent)
            ring0 = [[v.x, v.y] for v in
                     (tr(ms.Vector([x, y, 0]), src_crs, wgs84) for x, y in pl['ring'])]
            plot_area = _ring_area_m2(ring0)                       # true plot area, m^2
            row['plot_area_m2'] = f"{plot_area:.3f}"
            # Auto-fit: render WIDER than the grid so rows the grid clips are still
            # captured, then the worker crops back to the detected rows. Only the
            # width (short) edges are pushed out; length & the true area are kept.
            fit_probe_used = 0.0
            if auto_fit and len(ring0) == 4:
                latc = sum(q[1] for q in ring0) / 4.0
                _mx = 111320.0 * math.cos(math.radians(latc)); _my = 110540.0

                def _uv_wid(a, bpt):          # unit vector a->b in metres, back to deg
                    dx = (bpt[0] - a[0]) * _mx; dy = (bpt[1] - a[1]) * _my
                    d = math.hypot(dx, dy) or 1.0
                    return (dx / d / _mx, dy / d / _my)

                def _len_m(a, bpt):
                    return math.hypot((bpt[0] - a[0]) * _mx, (bpt[1] - a[1]) * _my)
                A, B, C, D = ring0
                pm = fit_probe_m
                # Push out the SHORT (across-row) edges only: rows run along the
                # plot's long axis, so it's the short axis that has to grow to
                # recover rows the grid clipped. Which vertex pair is short
                # depends on the shapefile's vertex order -> measure it.
                if (_len_m(A, B) + _len_m(C, D)) <= (_len_m(B, C) + _len_m(D, A)):
                    u1 = _uv_wid(A, B); u2 = _uv_wid(D, C)    # A-B / D-C are short
                    ring0 = [[A[0] - u1[0] * pm, A[1] - u1[1] * pm],
                             [B[0] + u1[0] * pm, B[1] + u1[1] * pm],
                             [C[0] + u2[0] * pm, C[1] + u2[1] * pm],
                             [D[0] - u2[0] * pm, D[1] - u2[1] * pm]]
                else:
                    u1 = _uv_wid(A, D); u2 = _uv_wid(B, C)    # A-D / B-C are short
                    ring0 = [[A[0] - u1[0] * pm, A[1] - u1[1] * pm],
                             [B[0] - u2[0] * pm, B[1] - u2[1] * pm],
                             [C[0] + u2[0] * pm, C[1] + u2[1] * pm],
                             [D[0] + u1[0] * pm, D[1] + u1[1] * pm]]
                fit_probe_used = pm
            ring_w = expand_wgs(ring0, margin)
            lon0 = sum(q[0] for q in ring_w) / len(ring_w)
            lat0 = sum(q[1] for q in ring_w) / len(ring_w)

            def to3d(lon, lat):
                v = tr(ms.Vector([lon, lat, 0]), wgs84, ch.crs)
                z = ground_z(ms.Vector([v.x, v.y]))     # per-vertex terrain height
                return Tinv.mulp(ch.crs.unproject(ms.Vector([v.x, v.y, z]))), z

            pts3d = [to3d(lon, lat)[0] for lon, lat in ring_w]      # margin ring (crop)
            pts3d_plot = [to3d(lon, lat)[0] for lon, lat in ring0]  # plot ring (coverage)

            if do_raw:
                # Score cameras for the most ACCURATE native crop, in priority order:
                #   1. most PLOT corners in-frame (coverage of the actual plot)
                #   2. crop extent fully inside the frame with an inset  -> no black
                #      edge in the native bounding-box crop, no resampling needed
                #   3. plot most centred in the frame (closest to nadir/principal
                #      point) -> sharpest, least perspective distortion, most uniform
                #      GSD. This mirrors EasyIDP's sort_img_by_distance heuristic.
                INSET = 8  # px a corner must sit inside the image edge to count as "clear"
                # For a rectified crop the frame is WARPED, so the priority is a
                # near-NADIR view (low foreshortening) — an oblique frame stretches
                # the far edge into streaks. Gap-fill then covers what one frame
                # misses. For native/mask (no warp) corner-count still leads.
                prefer_nadir = raw_mode in ('rectify', 'auto', 'best')
                mpx_g = 111320.0 * math.cos(math.radians(lat0))
                mpy_g = 110540.0

                def _gm(a, b):
                    return math.hypot((a[0] - b[0]) * mpx_g, (a[1] - b[1]) * mpy_g)
                if len(ring0) == 4:
                    egA = _gm(ring0[0], ring0[1]); egB = _gm(ring0[1], ring0[2])
                    egC = _gm(ring0[2], ring0[3]); egD = _gm(ring0[3], ring0[0])
                    gW = (egA + egC) / 2.0 or 1e-6       # width  (edges 0-1, 2-3)
                    gL = (egB + egD) / 2.0 or 1e-6       # length (edges 1-2, 3-0)
                else:
                    egA = egB = egC = egD = None; gW = gL = None
                cand = []
                for cam in cams:
                    uvs_p = [cam.project(q) for q in pts3d_plot]
                    if any(uv is None for uv in uvs_p):
                        continue
                    w, h = cam.sensor.width, cam.sensor.height
                    # cheap reject: cam.project happily extrapolates outside the
                    # sensor, so without this every camera in the flight (3000+)
                    # counted as a "candidate" and got 8 projections spent on it
                    if (max(uv.x for uv in uvs_p) < 0 or min(uv.x for uv in uvs_p) > w or
                            max(uv.y for uv in uvs_p) < 0 or min(uv.y for uv in uvs_p) > h):
                        continue
                    uvs = [cam.project(q) for q in pts3d]
                    if any(uv is None for uv in uvs):
                        continue
                    inside = sum(0 <= uv.x < w and 0 <= uv.y < h for uv in uvs_p)
                    over = max(max(-uv.x, uv.x - w, -uv.y, uv.y - h, 0) for uv in uvs_p)
                    crop_clear = all(INSET <= uv.x <= w - INSET and
                                     INSET <= uv.y <= h - INSET for uv in uvs)
                    cx = sum(uv.x for uv in uvs_p) / len(uvs_p)
                    cy = sum(uv.y for uv in uvs_p) / len(uvs_p)
                    centre_dist = math.hypot(cx - w / 2.0, cy - h / 2.0) / \
                        math.hypot(w / 2.0, h / 2.0)
                    # foreshortening (nadir->1) and resolution from the projected
                    # plot quad's 4 edge scales (frame px per ground metre)
                    aniso, gsd = 99.0, 0.0
                    if len(uvs_p) == 4 and gW:
                        sc = []
                        for (a, b), eg in (((0, 1), egA), ((1, 2), egB),
                                           ((2, 3), egC), ((3, 0), egD)):
                            dd = math.hypot(uvs_p[a].x - uvs_p[b].x,
                                            uvs_p[a].y - uvs_p[b].y)
                            if eg and eg > 1e-9:
                                sc.append(dd / eg)
                        if sc:
                            aniso = max(sc) / max(min(sc), 1e-9)
                            gsd = sorted(sc)[len(sc) // 2]
                    # EXACT share of the plot this one frame holds: the projected
                    # plot quad clipped to the frame rectangle. "Corners in frame"
                    # is 2/4 for a frame holding 51% and for one holding 95% of
                    # the plot, so it cannot pick the best single image; this can.
                    covq = frame_coverage([(uv.x, uv.y) for uv in uvs_p], w, h)
                    if covq <= 0.0:
                        continue                    # frame does not see the plot
                    if raw_mode == 'best':
                        # one image per plot, no stitching: take the frame that
                        # holds the LARGEST part of the plot, then the most nadir
                        # / most centred of those.
                        key = (-round(covq, 4), round(aniso, 3),
                               round(centre_dist, 3))
                    elif prefer_nadir:
                        # bucket foreshortening coarsely so that among comparably
                        # near-nadir frames, the one that COVERS more plot wins.
                        # This avoids picking a barely-clipping nadir frame (which
                        # leaves black) over a slightly-more-oblique full-coverage one.
                        key = (round(aniso / 0.35) * 0.35, -round(covq, 3),
                               round(centre_dist, 3), over)
                    else:
                        # native / mask / bbox: no warp, so what matters is how
                        # much of the plot this frame holds (was corner count,
                        # which is 2/4 for an 81% frame and for a 99% one alike),
                        # then whether the crop box is clear of the image edge
                        # (no black border), then how centred the plot is.
                        key = (-round(covq, 4), -int(crop_clear),
                               round(centre_dist, 4), over)
                    cand.append((key, cam, uvs, uvs_p, w, h, inside, aniso, gsd, covq))
                cand.sort(key=lambda c: c[0])
                row['cand_frames'] = str(len(cand))
                if cand:
                    bestc = max(cand, key=lambda c: c[9])
                    row['best_frame'] = bestc[1].label
                    row['best_cover_pct'] = f"{100 * bestc[9]:.1f}"
                    row['best_cover_m2'] = f"{plot_area * bestc[9]:.3f}"
                    if candidates_csv is not None:
                        chosen = cand[0][1].label
                        with _cand_lock:
                            for crank, c in enumerate(cand[:CAND_REPORT_N], 1):
                                candidates_csv.append(dict(
                                    plot_id=pid, rank=crank, camera=c[1].label,
                                    cover_pct=f"{100 * c[9]:.1f}",
                                    cover_m2=f"{plot_area * c[9]:.3f}",
                                    corners_in_frame=f"{c[6]}/4",
                                    aniso=f"{c[7]:.3f}",
                                    mm_per_px=(f"{1000.0 / c[8]:.3f}" if c[8] else ''),
                                    used=('yes' if c[1].label == chosen else '')))
                if not cand:
                    row['note'] = 'no covering raw image'
                    log(f"plot {pid}: NO covering raw image")
                else:
                    # PHASE A (this main thread): resolve all Metashape geometry and
                    # queue a render job. The heavy pixel work runs later in parallel.
                    mesh_cache = {}       # (nx,ny) -> 3D ground nodes for this plot

                    def _mesh_nodes3d(nx, ny):
                        """Ground points of the rectification mesh, shared by every
                        frame of this plot - that shared grid is what makes the
                        frames land on each other."""
                        key = (nx, ny)
                        got = mesh_cache.get(key)
                        if got is None:
                            got = [to3d(lon, lat)[0]
                                   for lon, lat in mesh_uv_nodes(ring0, nx, ny)]
                            mesh_cache[key] = got
                        return got

                    def _project_nodes(cam_, nodes3d):
                        out_ = []
                        for q in nodes3d:
                            uv = cam_.project(q)
                            out_.append(None if uv is None else (uv.x, uv.y))
                        return out_

                    for rank, (_s, cam, uvs, uvs_p, w, h, inside, aniso, gsd, covq) in \
                            enumerate(cand[:n_per_plot], 1):
                        xs = [uv.x for uv in uvs]
                        ys = [uv.y for uv in uvs]
                        if raw_mode == 'native':
                            # Centre the crop on the projected PLOT centroid and pad
                            # symmetrically, so the plot sits in the middle with equal
                            # margin on every side. Oblique perspective makes the plain
                            # bounding box lopsided (equal ground margin -> unequal pixel
                            # margin); centring fixes that. Still native pixels, no resample.
                            cxp = sum(uv.x for uv in uvs_p) / len(uvs_p)
                            cyp = sum(uv.y for uv in uvs_p) / len(uvs_p)
                            half_w = max(abs(uv.x - cxp) for uv in uvs)
                            half_h = max(abs(uv.y - cyp) for uv in uvs)
                            x0 = max(0, int(math.floor(cxp - half_w)))
                            y0 = max(0, int(math.floor(cyp - half_h)))
                            x1 = min(w, int(math.ceil(cxp + half_w)))
                            y1 = min(h, int(math.ceil(cyp + half_h)))
                        else:
                            x0 = max(0, int(min(xs))); y0 = max(0, int(min(ys)))
                            x1 = min(w, int(max(xs))); y1 = min(h, int(max(ys)))
                        if x1 <= x0 or y1 <= y0:
                            row['note'] = 'degenerate crop'
                            continue
                        path = resolve_photo(cam)
                        if path is None:
                            if rank == 1:
                                row['note'] = f'raw file missing ({cam.label})'
                            log(f"plot {pid}: raw file missing for {cam.label} "
                                f"(set the raw image folder?)")
                            continue
                        fn = f"plot_{pid}{ext}" if rank == 1 else f"plot_{pid}_{rank}{ext}"
                        out = os.path.join(raw_dir, fn)
                        # plot polygon corners in crop-local pixel coords
                        plot_uv = [(uv.x - x0, uv.y - y0) for uv in uvs_p]
                        eff = raw_mode
                        if eff == 'auto':                       # straighten tilted plots only
                            eff = 'rectify' if _tilt_from_vertical(plot_uv) > 3.0 else 'mask'
                        elif eff == 'best':                     # one frame, upright, no stitch
                            eff = 'rectify'
                        rectified = (eff == 'rectify' and len(plot_uv) == 4)
                        # output crop size, computed analytically (no pixels needed yet)
                        if rectified:
                            # Size from the plot's TRUE GROUND dimensions at this frame's
                            # representative (median) resolution -> correct aspect and an
                            # upright strip, immune to the oblique blow-up that made some
                            # crops huge & near-square. Falls back to the pixel-edge
                            # average if ground dims are unavailable.
                            if gW and gL and gsd > 1e-9:
                                wpx = max(1, int(round(gW * gsd)))
                                hpx = max(1, int(round(gL * gsd)))
                            else:
                                A, B, C, D = plot_uv

                                def _d(a, b):
                                    return math.hypot(a[0] - b[0], a[1] - b[1])
                                wpx = max(1, int(round((_d(A, B) + _d(D, C)) / 2)))
                                hpx = max(1, int(round((_d(A, D) + _d(B, C)) / 2)))
                            MAXPX = 20000                    # sane pixel-budget cap
                            wpx = min(wpx, MAXPX); hpx = min(hpx, MAXPX)
                            cw2, ch2 = wpx, hpx
                        else:
                            wpx = hpx = None
                            cw2, ch2 = x1 - x0, y1 - y0
                        gcps = []
                        if geotag:
                            if rectified and len(ring0) == 4:
                                # output rectangle corners map to the bare plot corners
                                for (gx, gy), (lon, lat) in zip(
                                        [(0, 0), (0, ch2), (cw2, ch2), (cw2, 0)],
                                        [ring0[0], ring0[3], ring0[2], ring0[1]]):
                                    _p3, zg = to3d(lon, lat)
                                    gcps.append((gy, gx, lon, lat, zg))
                            else:
                                # 3x3 grid over the margin quad, projected into crop px
                                if len(ring_w) == 4:
                                    A, B, C, D = ring_w
                                    grid = [tuple((1-s)*(1-t)*a + s*(1-t)*b + s*t*c + (1-s)*t*d
                                                  for a, b, c, d in zip(A, B, C, D))
                                            for s in (0.0, 0.5, 1.0) for t in (0.0, 0.5, 1.0)]
                                else:
                                    grid = [tuple(q) for q in ring_w] + [(lon0, lat0)]
                                for lon, lat in grid:
                                    p3, zg = to3d(lon, lat)
                                    uv = cam.project(p3)
                                    if uv is None:
                                        continue
                                    px, py = uv.x - x0, uv.y - y0
                                    if -0.5 <= px <= (x1-x0)+0.5 and -0.5 <= py <= (y1-y0)+0.5:
                                        gcps.append((py, px, lon, lat, zg))
                        # Rectification mesh: the plot's ground grid projected into
                        # this frame with the FULL camera model. Every frame of the
                        # plot uses the same ground nodes, so the frames overlay
                        # each other exactly (see _rectify).
                        nx = ny = 0
                        self_layer = None
                        if rectified:
                            nx, ny = mesh_dims(cw2, ch2, mesh_step)
                            n3d = _mesh_nodes3d(nx, ny)
                            self_layer = dict(path=path, aniso=aniso, nx=nx, ny=ny,
                                              nodes=_project_nodes(cam, n3d),
                                              corners=[(uv.x, uv.y) for uv in uvs_p])
                        # Gap-fill: for a rectified rank-1 crop, gather the top
                        # near-nadir candidate frames (cand is already sorted
                        # nadir-first) so a neighbour can fill what the best single
                        # frame misses. Each layer carries its foreshortening so the
                        # renderer can skip oblique frames that would smear.
                        # In "best single frame" mode there is deliberately no
                        # gap-fill: one untouched frame per plot, zero stitching.
                        layers = None
                        if rank == 1 and rectified and raw_mode != 'best':
                            # Choose the FEWEST frames that still cover the whole
                            # plot (greedy set cover on a coarse ground grid).
                            # Every extra frame is another seam, so taking the
                            # highest-coverage frames in rank order - which is
                            # what the old code did - could pull in 4-5 frames
                            # where 2 suffice. Fewer frames = fewer boundaries.
                            c3d = _mesh_nodes3d(12, 4)
                            pool = cand[:max(2 * GAPFILL_N, 12)]
                            masks = []
                            for c in pool:
                                w2, h2 = c[4], c[5]
                                masks.append([
                                    (q is not None and 0 <= q[0] < w2 and 0 <= q[1] < h2)
                                    for q in _project_nodes(c[1], c3d)])
                            covered = list(masks[0])
                            picked = [0]
                            while (not all(covered)) and len(picked) < GAPFILL_N:
                                bi, bgain = None, 0
                                for i2, m in enumerate(masks):
                                    if i2 in picked:
                                        continue
                                    gain = sum(1 for k, cv in enumerate(covered)
                                               if m[k] and not cv)
                                    if gain > bgain:
                                        bi, bgain = i2, gain
                                if bi is None:
                                    break
                                picked.append(bi)
                                covered = [a or b for a, b in zip(covered, masks[bi])]
                            layers = [self_layer]
                            for i2 in picked[1:]:
                                cam2 = pool[i2][1]
                                p2 = resolve_photo(cam2)
                                if p2 is None:
                                    continue
                                layers.append(dict(
                                    path=p2, aniso=pool[i2][7], nx=nx, ny=ny,
                                    nodes=_project_nodes(cam2, n3d),
                                    corners=[(uv.x, uv.y) for uv in pool[i2][3]]))
                            if len(layers) <= 1:
                                layers = None          # one frame covers it all
                        # ground metres-per-pixel of the rectified output (for
                        # auto-fit area + width in real units)
                        mpp_x = (gW / wpx) if (rectified and gW and wpx) else 0.0
                        mpp_y = (gL / hpx) if (rectified and gL and hpx) else 0.0
                        # image axis running ACROSS the rows = the plot's SHORT
                        # ground axis (x is the gW edge, y is the gL edge)
                        fit_across = ('x' if (not gW or not gL or gW <= gL)
                                      else 'y')
                        job = dict(pid=pid, path=path, x0=x0, y0=y0, x1=x1, y1=y1,
                                   eff=eff, plot_uv=plot_uv, wpx=wpx, hpx=hpx,
                                   out=out, geotag=geotag, gcps=gcps, layers=layers,
                                   is_rank1=(rank == 1), row=row, cam_label=cam.label,
                                   cw=cw2, ch=ch2, full=False, cov=0.0, flag='',
                                   auto_fit=(auto_fit and rectified and rank == 1),
                                   fit_rows=fit_rows, fit_mode=fit_mode,
                                   mpp_x=mpp_x, mpp_y=mpp_y,
                                   fit_across=fit_across, self_layer=self_layer,
                                   smart_seam=smart_seam, seam_map=seam_map,
                                   seam_slack=seam_slack)
                        if rank == 1:
                            # coverage of the plot by THIS single frame, from the
                            # exact clipped-quad area (replaces the old 15x15 point
                            # sampling: same answer, ~200x cheaper, no quantisation)
                            cov = 100.0 * covq
                            # Judge completeness by real coverage, not by how many
                            # of the 4 corners are in frame: this frame holds 99%
                            # of the plot with only 1 corner inside, which the old
                            # "1/4 corners" warning made look like a failure.
                            job['full'] = cov >= 99.5
                            job['cov'] = cov
                            job['flag'] = "" if cov >= 99.5 else \
                                f"  [!] {cov:.1f}% of the plot in this frame " \
                                f"({inside}/{len(uvs)} corners)"
                            row.update(raw_image=cam.label, raw_file=fn,
                                       crop_w=cw2, crop_h=ch2,
                                       corners_in_frame=f"{inside}/{len(uvs)}",
                                       raw_covered_pct=f"{cov:.1f}",
                                       raw_covered_m2=f"{plot_area * cov / 100:.3f}",
                                       source_path=path)
                        fut = raw_pool.submit(_render, job)
                        fut.add_done_callback(_on_raw_done)
                        raw_futures.append(fut)

            if do_ortho:
                res_o = export_ortho_crop(ms, Image, ch, tr, wgs84, units_per_metre,
                                          ring_w, margin, pid, ortho_dir, ext, save_kw,
                                          geotag, dict(PLOT_ID=str(pid)) if geotag else None)
                if res_o:
                    fn, ocov = res_o
                    ok_ortho += 1
                    row['ortho_file'] = fn
                    row['ortho_covered_pct'] = f"{ocov:.1f}"
                    row['ortho_covered_m2'] = f"{plot_area * ocov / 100:.3f}"
                    log(f"ortho  plot {pid}  {idx + 1}/{total}  cover {ocov:.0f}%")
        except Exception as e:
            row['note'] = f'error: {e}'
            log(f"plot {pid}: ERROR {e}")
        manifest.append(row)
        set_progress(idx + 1, total)

    # Wait for the raw-crop renders that were submitted during the loop. Each
    # one already logged itself live via _on_raw_done as it finished.
    if raw_pool is not None:
        nj = len(raw_futures)
        if nj:
            log("")
            log(f"Finishing raw crops ({_rawcount['done']}/{nj} already rendered) ...")
        raw_pool.shutdown(wait=True)
        ok_raw = _rawcount['ok']
        full_raw = _rawcount['full']

    # Ortho fallback: for field-edge plots where NO near-nadir raw frame exists
    # (raw came out mostly black), rebuild the raw crop from the plot's ORTHO
    # crop (full coverage, correct location), auto-fit to the rows. Every ortho
    # pixel is native raw stitched by Metashape - softer than a good raw frame,
    # but a complete, correctly-located crop beats a smeared/black one.
    if _fallback and do_ortho:
        import numpy as np
        log("")
        log(f"Ortho fallback for {len(_fallback)} low-coverage plot(s) "
            f"(no adequate raw frame; rebuilding from the orthomosaic crop) ...")
        for pid, row, rawout in _fallback:
            ofn = row.get('ortho_file')
            opath = os.path.join(ortho_dir, ofn) if ofn else None
            if not (opath and os.path.exists(opath)):
                log(f"plot {pid}: ortho fallback skipped (no ortho crop)")
                continue
            try:
                oim = Image.open(opath).convert('RGB')
                oa = np.asarray(oim).copy()
                oa[oa.min(2) > 235] = 0                # white nodata -> black
                oim = Image.fromarray(oa, 'RGB')
                if auto_fit:
                    fb = detect_row_fit(oim, rows_target=fit_rows, mode=fit_mode)
                    if fb:
                        x0, x1, y0, y1 = fb['box']
                        oim = oim.crop((x0, y0, x1, y1))
                        row['fit_rows'] = str(fb['rows'])
                oim.save(rawout, **save_kw)
                blk = 100.0 * float((np.asarray(oim).sum(2) < 24).mean())
                row['crop_w'], row['crop_h'] = str(oim.width), str(oim.height)
                row['raw_covered_pct'] = f"{100 - blk:.1f}"
                row['source_path'] = opath
                row['note'] = ("ortho-fallback: field-edge, no near-nadir raw frame; "
                               "rebuilt from orthomosaic crop, auto-fit rows "
                               "(nodata=black)")
                log(f"  plot {pid} <- ortho crop  {oim.width}x{oim.height}px  "
                    f"nodata {blk:.0f}%")
            except Exception as e:
                log(f"plot {pid}: ortho fallback FAILED ({e})")

    # per-plot candidate report: every raw frame that sees the plot, how much of
    # the plot it holds, and which one was used. This is the "all the images for
    # this plot / which is the best one" record.
    if candidates_csv:
        os.makedirs(p['outdir'], exist_ok=True)
        cpath = os.path.join(p['outdir'], "candidates.csv")
        try:
            with open(cpath, 'w', newline='', encoding='utf-8') as f:
                wc = csv.DictWriter(f, fieldnames=list(candidates_csv[0].keys()))
                wc.writeheader()
                wc.writerows(candidates_csv)
            log(f"Candidate frames per plot: {cpath} ({len(candidates_csv)} rows)")
        except OSError as e:
            log(f"[!] Could not write candidates.csv: {e}")

    # manifest
    os.makedirs(p['outdir'], exist_ok=True)
    mpath = os.path.join(p['outdir'], "manifest.csv")
    with open(mpath, 'w', newline='', encoding='utf-8') as f:
        wcsv = csv.DictWriter(f, fieldnames=list(manifest[0].keys()) if manifest else ['plot_id'])
        wcsv.writeheader()
        wcsv.writerows(manifest)
    log("")
    log(f"DONE: raw crops {ok_raw}/{len(plots)}"
        + (f" ({full_raw} with the whole plot in frame)" if do_raw else "")
        + (f", ortho crops {ok_ortho}/{len(plots)}" if do_ortho else "")
        + f"  ->  {p['outdir']}")
    log(f"Manifest: {mpath}")
    return p['outdir']


_GCP_TR = {}


_GCP_EPSG = [None]          # set per run from the 'gcp_epsg' parameter


def _gcp_target():
    """Target CRS for the GCPs written into GeoTIFF crops.

    Geometry is normalised to EPSG:4326 internally (see the wgs84 CoordinateSystem
    and the plot-ring transform on load), so GCPs are lon/lat by default. Set
    PE_GCP_EPSG to write them in another CRS instead - the coordinates are
    TRANSFORMED, not merely relabelled. Use this to keep a whole project in one
    system, e.g. PE_GCP_EPSG=7850 for GDA2020 / MGA zone 50.

    Note the default 4326 label is inherited from the Metashape chunk CRS. For an
    Australian project those degrees are really GDA2020 values, whose correct code
    is 7844 - numerically identical, but the label is loose. Setting this variable
    removes the ambiguity.
    """
    v = _GCP_EPSG[0] or os.environ.get('PE_GCP_EPSG')
    if not v:
        return 4326, None
    try:
        tgt = int(v)
    except (TypeError, ValueError):
        return 4326, None
    if tgt == 4326:
        return 4326, None
    if tgt not in _GCP_TR:
        from pyproj import Transformer
        _GCP_TR[tgt] = Transformer.from_crs("EPSG:4326", f"EPSG:{tgt}",
                                            always_xy=True)
    return tgt, _GCP_TR[tgt]


def retag_gcps(folder, epsg, log=print, progress=None):
    """Re-express the ground control points of every GeoTIFF in `folder`.

    Coordinates are TRANSFORMED into the target CRS, not relabelled, and the
    pixels are untouched - only the GCP block is rewritten, about 4 KB per file.

    This exists because the output CRS is decided at extraction time: a set
    produced before it was set ships in the internal EPSG:4326, and a delivery
    can end up half in one system and half in another. That happened once here,
    and this is the repair.

    Returns (converted, already_correct, skipped, failed).
    """
    import glob
    import rasterio
    from rasterio.control import GroundControlPoint
    from rasterio.crs import CRS
    from pyproj import Transformer

    target = int(epsg)
    files = sorted(glob.glob(os.path.join(folder, "*.tif"))
                   + glob.glob(os.path.join(folder, "*.tiff")))
    if not files:
        log(f"[!] no GeoTIFFs found in {folder}")
        return 0, 0, 0, 0
    log(f"Re-expressing GCPs of {len(files)} file(s) in EPSG:{target} - "
        f"pixels are not touched.")
    trs, done, same, skip, bad = {}, 0, 0, 0, 0
    for i, f in enumerate(files):
        try:
            with rasterio.open(f) as s:
                g, c = s.gcps
                src_epsg = c.to_epsg() if c else None
            if not g:
                skip += 1
                continue
            if src_epsg == target:
                same += 1
                continue
            if src_epsg is None:
                log(f"[!] {os.path.basename(f)}: GCPs carry no CRS - skipped")
                skip += 1
                continue
            if src_epsg not in trs:
                trs[src_epsg] = Transformer.from_crs(f"EPSG:{src_epsg}",
                                                     f"EPSG:{target}",
                                                     always_xy=True)
            tf = trs[src_epsg]
            ng = []
            for q in g:
                x, y = tf.transform(q.x, q.y)
                ng.append(GroundControlPoint(row=q.row, col=q.col,
                                             x=x, y=y, z=q.z))
            with rasterio.open(f, "r+") as d:
                d.gcps = (ng, CRS.from_epsg(target))
            done += 1
        except Exception as e:
            log(f"[!] {os.path.basename(f)}: {e}")
            bad += 1
        if progress and (i + 1) % 25 == 0:
            progress(i + 1, len(files))
    log(f"GCP CRS: {done} converted, {same} already EPSG:{target}, "
        f"{skip} skipped, {bad} failed.")
    return done, same, skip, bad


def save_geotiff(path, pil_img, gcps, tags, log):
    """Write an RGB crop as a lossless GeoTIFF with GCPs and the source EXIF
    copied into the GDAL metadata (readable via gdalinfo / rasterio.tags()).

    GCPs are written in EPSG:4326 unless PE_GCP_EPSG asks for another CRS."""
    import numpy as np
    import rasterio
    from rasterio.control import GroundControlPoint
    from rasterio.crs import CRS
    arr = np.asarray(pil_img)
    prof = dict(driver='GTiff', height=arr.shape[0], width=arr.shape[1], count=3,
                dtype='uint8', compress='lzw', photometric='RGB',
                tiled=True, blockxsize=512, blockysize=512, bigtiff='IF_SAFER')
    with rasterio.open(path, 'w', **prof) as dst:
        dst.write(arr.transpose(2, 0, 1))
        if gcps:
            epsg, tf = _gcp_target()
            if tf is not None:
                pts = []
                for r, c, x, y, z in gcps:
                    X, Y = tf.transform(x, y)
                    pts.append(GroundControlPoint(row=r, col=c, x=X, y=Y, z=z))
            else:
                pts = [GroundControlPoint(row=r, col=c, x=x, y=y, z=z)
                       for r, c, x, y, z in gcps]
            dst.gcps = (pts, CRS.from_epsg(epsg))
        else:
            log(f"[!] {os.path.basename(path)}: no GCPs landed inside the crop")
        if tags:
            dst.update_tags(**tags)


def export_ortho_crop(ms, Image, ch, tr, wgs84, units_per_metre,
                      ring_w, margin, pid, outdir, ext, save_kw,
                      geotag=False, tags=None):
    """Export a native-resolution ortho crop of one plot (ring_w already includes
    the margin, WGS84). geotag=True keeps Metashape's own georeferencing (axis-
    aligned GeoTIFF, GDAL-ready); otherwise the plot is deskewed to a tight
    rectangle. Returns filename or None."""
    import tempfile
    res = ch.orthomosaic.resolution                     # metres / pixel
    lon0 = sum(q[0] for q in ring_w) / len(ring_w)
    lat0 = sum(q[1] for q in ring_w) / len(ring_w)
    ux, uy = units_per_metre(lon0, lat0)                # chunk units per metre
    rx, ry = res * ux, res * uy                         # chunk units per pixel
    ring_c = [tr(ms.Vector([lon, lat, 0]), wgs84, ch.crs) for lon, lat in ring_w]
    pts = [(v.x, v.y) for v in ring_c]
    pad = max(2.0 * margin, 0.5)                        # metres of slack around bbox
    x0 = min(q[0] for q in pts) - pad * ux
    x1 = max(q[0] for q in pts) + pad * ux
    y0 = min(q[1] for q in pts) - pad * uy
    y1 = max(q[1] for q in pts) + pad * uy
    box = ms.BBox()
    box.min = ms.Vector([x0, y0])
    box.max = ms.Vector([x1, y1])
    if geotag:
        # axis-aligned GeoTIFF with Metashape's exact geotransform + CRS, whole
        # plot in frame (rotated), everything outside the plot polygon set to
        # nodata so the file stays small. Plot pixels are untouched (lossless).
        import numpy as np
        import rasterio
        from rasterio.features import geometry_mask
        fn = f"plot_{pid}.tif"
        outp = os.path.join(outdir, fn)
        tmpg = os.path.join(tempfile.gettempdir(), f"_plotg_{os.getpid()}.tif")
        ch.exportRaster(path=tmpg, source_data=ms.OrthomosaicData, region=box,
                        resolution_x=rx, resolution_y=ry, save_alpha=False,
                        image_format=ms.ImageFormatTIFF)
        try:
            with rasterio.open(tmpg) as s:
                arr = s.read()
                prof = s.profile
            poly = dict(type='Polygon', coordinates=[ring_w + [ring_w[0]]])
            keep = geometry_mask([poly], out_shape=arr.shape[1:],
                                 transform=prof['transform'], invert=True)
            arr *= keep
            valid = arr.max(axis=0) > 0                     # any channel non-nodata
            kc = int(keep.sum())
            covered = 100.0 * float((valid & keep).sum()) / kc if kc else 0.0
            prof.update(compress='lzw', predictor=2, nodata=0, tiled=True,
                        blockxsize=512, blockysize=512, bigtiff='IF_SAFER')
            with rasterio.open(outp, 'w', **prof) as dst:
                dst.write(arr)
                if tags:
                    dst.update_tags(**tags)
        finally:
            try:
                os.remove(tmpg)
            except OSError:
                pass
        return fn, covered
    tmp = os.path.join(tempfile.gettempdir(), f"_plotx_{os.getpid()}.tif")
    try:
        ch.exportRaster(path=tmp, source_data=ms.OrthomosaicData, region=box,
                        resolution_x=rx, resolution_y=ry, save_alpha=False,
                        image_format=ms.ImageFormatTIFF)
    except Exception as e:
        raise RuntimeError(f"ortho export failed: {str(e)[:80]}")
    try:
        with Image.open(tmp) as im0:
            im = im0.convert("RGB")
            im.load()
    finally:
        try:
            os.remove(tmp)
        except OSError:
            pass

    def to_px(q):
        return ((q[0] - x0) / rx, (y1 - q[1]) / ry)     # row 0 = max northing

    fn = f"plot_{pid}{ext}"
    if len(pts) == 4:
        # deskew the (rotated) quad to an axis-aligned rectangle
        A, B, C, D = [to_px(q) for q in pts]
        def dist(a, b):
            return math.hypot(a[0] - b[0], a[1] - b[1])
        wpx = max(1, int(round((dist(A, B) + dist(D, C)) / 2)))
        hpx = max(1, int(round((dist(A, D) + dist(B, C)) / 2)))
        quad = [A[0], A[1], D[0], D[1], C[0], C[1], B[0], B[1]]
        out = im.transform((wpx, hpx), Image.QUAD, quad, Image.BICUBIC)
    else:
        out = im                                        # non-quad ring: keep bbox crop
    import numpy as np
    oa = np.asarray(out.convert("RGB")).copy()
    # the Metashape orthomosaic fills uncovered ground WHITE; recolor that nodata
    # to BLACK so it reads as a consistent "no data" region (not a bright strip
    # that would be miscounted as soil / skew a whole-plot radiometric analysis)
    nod = oa.min(axis=2) > 235
    valid = (oa.max(axis=2) > 0) & (~nod)               # real imagery = not black/white
    oa[nod] = 0
    out = Image.fromarray(oa, "RGB")
    covered = 100.0 * float(valid.sum()) / (oa.shape[0] * oa.shape[1])
    out.save(os.path.join(outdir, fn), **save_kw)
    return fn, covered


# --------------------------------------------------------------------------
# Ortho-only mode: crop plots straight from a GeoTIFF (no Metashape)
# --------------------------------------------------------------------------
def export_ortho_crop_from_raster(ds, np, rasterio, Image, ring_r, pid, outdir,
                                  ext, save_kw, geotag, tags, pad_m, upmx, upmy):
    """Crop one plot out of an open rasterio orthomosaic. ring_r is the (margin-
    expanded) plot ring in the raster CRS. geotag=True writes an axis-aligned
    GeoTIFF with the source geotransform and everything outside the polygon set to
    nodata (lossless); otherwise the plot is deskewed to a tight upright rectangle.
    Returns (filename, covered_pct, width, height) or None if the plot is off-raster."""
    from rasterio.windows import Window
    from rasterio.features import geometry_mask
    xs = [q[0] for q in ring_r]
    ys = [q[1] for q in ring_r]
    # upmx/upmy = metres per CRS unit (metres-per-degree for a geographic CRS, 1 for
    # a projected one), so a pad in metres becomes CRS units by dividing.
    padx, pady = pad_m / upmx, pad_m / upmy
    minx, maxx = min(xs) - padx, max(xs) + padx
    miny, maxy = min(ys) - pady, max(ys) + pady
    # bounding window in pixel space (clamped to the raster) via the affine inverse
    inv = ~ds.transform
    corners_px = [inv * (cx, cy) for cx in (minx, maxx) for cy in (miny, maxy)]
    cols = [c for c, _ in corners_px]
    rowspx = [rr for _, rr in corners_px]
    col_off = max(0, int(math.floor(min(cols))))
    row_off = max(0, int(math.floor(min(rowspx))))
    col_end = min(ds.width, int(math.ceil(max(cols))))
    row_end = min(ds.height, int(math.ceil(max(rowspx))))
    if col_end <= col_off or row_end <= row_off:
        return None
    win = Window(col_off, row_off, col_end - col_off, row_end - row_off)
    arr = ds.read(window=win)                              # (bands, h, w)
    wt = ds.window_transform(win)
    ring_closed = [tuple(q) for q in ring_r] + [tuple(ring_r[0])]
    poly = dict(type='Polygon', coordinates=[ring_closed])
    fn = f"plot_{pid}{ext}"
    outp = os.path.join(outdir, fn)

    if geotag:
        keep = geometry_mask([poly], out_shape=(arr.shape[1], arr.shape[2]),
                             transform=wt, invert=True)
        band_n = min(3, arr.shape[0])
        arr = arr[:band_n] * keep
        valid = arr.max(axis=0) > 0
        kc = int(keep.sum())
        covered = 100.0 * float((valid & keep).sum()) / kc if kc else 0.0
        prof = ds.profile.copy()
        prof.update(driver='GTiff', height=arr.shape[1], width=arr.shape[2],
                    count=band_n, transform=wt, compress='lzw', predictor=2,
                    nodata=0, tiled=True, blockxsize=512, blockysize=512,
                    bigtiff='IF_SAFER')
        prof.pop('photometric', None)
        with rasterio.open(outp, 'w', **prof) as dst:
            dst.write(arr)
            if tags:
                dst.update_tags(**tags)
        return fn, covered, int(arr.shape[2]), int(arr.shape[1])

    # deskew to an upright rectangle for plain PNG/TIFF/JPEG output
    rgb = arr[:3] if arr.shape[0] >= 3 else np.repeat(arr[:1], 3, axis=0)
    if rgb.dtype != np.uint8:                              # 16-bit ortho -> 8-bit view
        mx = float(rgb.max()) or 1.0
        rgb = np.clip(rgb.astype('float32') * (255.0 / mx), 0, 255).astype('uint8')
    im = Image.fromarray(np.transpose(rgb, (1, 2, 0)), 'RGB')
    inv = ~wt

    def to_px(q):
        c, r = inv * (q[0], q[1])
        return (c, r)

    if len(ring_r) == 4:
        A, B, C, D = [to_px(q) for q in ring_r]
        def dist(a, b):
            return math.hypot(a[0] - b[0], a[1] - b[1])
        wpx = max(1, int(round((dist(A, B) + dist(D, C)) / 2)))
        hpx = max(1, int(round((dist(A, D) + dist(B, C)) / 2)))
        quad = [A[0], A[1], D[0], D[1], C[0], C[1], B[0], B[1]]
        out = im.transform((wpx, hpx), Image.QUAD, quad, Image.BICUBIC)
    else:
        out = im
    oa = np.asarray(out.convert('RGB'))
    covered = 100.0 * float((oa.max(axis=2) > 0).sum()) / (oa.shape[0] * oa.shape[1])
    out.save(outp, **save_kw)
    return fn, covered, out.size[0], out.size[1]


def run_extraction_ortho_tif(p, log, set_progress, cancel):
    """Ortho-only extraction: crop one deskewed orthomosaic tile per plot directly
    from a georeferenced GeoTIFF, WITHOUT opening Metashape (rasterio only)."""
    import numpy as np
    import rasterio
    from rasterio.crs import CRS
    from rasterio.warp import transform as warp_transform
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None

    ortho_path = (p.get('ortho_tif') or '').strip()
    if not (ortho_path and os.path.exists(ortho_path)):
        raise RuntimeError("Select a valid orthomosaic GeoTIFF for ortho-only mode.")

    plots, wkt, fields = read_plot_file(p['plots'])
    if not plots:
        raise RuntimeError("No polygons found in " + p['plots'])
    log(f"Loaded {len(plots)} plot polygons from {os.path.basename(p['plots'])}")

    # resolve the plot source CRS as a rasterio CRS
    override = (p.get('crs_override') or '').strip()
    if override and override.lower() != 'auto':
        spec = re.sub(r'^EPSG:+', '', override, flags=re.I)
        try:
            src_crs = CRS.from_epsg(int(spec)) if spec.isdigit() \
                else CRS.from_user_input(override)
        except Exception:
            raise RuntimeError(f"Could not parse CRS override '{override}'. "
                               "Use an EPSG code like 28350, or leave it on auto.")
        log(f"Plot CRS (user override): {src_crs.to_string()}")
    elif wkt:
        try:
            src_crs = CRS.from_wkt(wkt)
        except Exception:
            try:
                src_crs = CRS.from_user_input(wkt)
            except Exception:
                raise RuntimeError("Could not parse the plot file CRS. "
                                   "Set the 'Plot CRS' EPSG code manually.")
        log(f"Plot CRS (from file): {src_crs.to_string()}")
    else:
        src_crs = CRS.from_epsg(4326)
        xs = [abs(v) for pl in plots[:5] for v in pl['ring'][0]]
        if xs and max(xs) > 360:
            raise RuntimeError("The plot file carries no CRS and its coordinates do "
                               "not look like lon/lat. Set the 'Plot CRS' EPSG code.")
        log("Plot CRS: assumed WGS84 lon/lat (file carries no CRS)")

    id_field = p.get('id_field')
    if id_field in (None, '', AUTO_ID):
        id_field = pick_id_field(fields)
    for i, pl in enumerate(plots, 1):
        pl['id'] = safe_name(pl['props'].get(id_field, i)) if id_field else str(i)
    limit = int(p.get('limit') or 0)
    if limit:
        plots = plots[:limit]
        log(f"Limit: first {limit} plots only (test run)")

    margin = float(p.get('margin') or 0.0)
    fmt = (p.get('format') or 'TIFF').upper()
    geotag = fmt.startswith('GEO')
    # CRS for the GeoTIFF ground control points. Coordinates are transformed into
    # it, not relabelled; blank keeps the internal EPSG:4326.
    try:
        _GCP_EPSG[0] = int(str(p.get('gcp_epsg') or '').strip() or 0) or None
    except (TypeError, ValueError):
        _GCP_EPSG[0] = None
    if geotag and _GCP_EPSG[0]:
        log(f"GCP CRS: writing ground control points in EPSG:{_GCP_EPSG[0]}")
    if geotag or fmt.startswith('TIF'):
        ext, save_kw = '.tif', dict(compression='tiff_lzw')
    elif fmt.startswith('PNG'):
        ext, save_kw = '.png', {}
    else:
        ext, save_kw = '.jpg', dict(quality=int(p.get('quality') or 92))
    log("Output format: " + ("GeoTIFF (lossless, masked to plot)" if geotag
                             else f"{ext[1:].upper()} deskewed rectangle"))

    ortho_dir = os.path.join(p['outdir'], "OrthoCrops")
    os.makedirs(ortho_dir, exist_ok=True)

    wgs84 = CRS.from_epsg(4326)
    manifest = []
    ok = 0
    total = len(plots)
    with rasterio.open(ortho_path) as ds:
        if ds.crs is None:
            raise RuntimeError("The orthomosaic GeoTIFF has no CRS / is not georeferenced.")
        log(f"Orthomosaic: {os.path.basename(ortho_path)}  {ds.width}x{ds.height}px, "
            f"{ds.count} bands, CRS {ds.crs.to_string()}, dtype {ds.dtypes[0]}")
        if not ds.is_tiled:
            log("[!] This ortho is STRIPED (not tiled): per-plot reads will be slow. "
                "Use 'Make tiled copy' once for a much faster tiled/overview version.")
        geographic = ds.crs.is_geographic
        for idx, pl in enumerate(plots):
            if cancel.is_set():
                log("Cancelled by user.")
                break
            pid = pl['id']
            row = dict(plot_id=pid, plot_area_m2='', margin_m=f"{margin:.2f}",
                       ortho_file='', crop_w='', crop_h='',
                       ortho_covered_pct='', ortho_covered_m2='',
                       source_path=ortho_path, note='')
            try:
                xs_s = [q[0] for q in pl['ring']]
                ys_s = [q[1] for q in pl['ring']]
                rx, ry = warp_transform(src_crs, ds.crs, xs_s, ys_s)
                ring_r = [[a, b] for a, b in zip(rx, ry)]
                wlon, wlat = warp_transform(src_crs, wgs84, xs_s, ys_s)
                ring_wgs = [[a, b] for a, b in zip(wlon, wlat)]
                area = _ring_area_m2(ring_wgs)
                row['plot_area_m2'] = f"{area:.3f}"
                cyd = sum(q[1] for q in ring_wgs) / len(ring_wgs)
                if geographic:
                    upmx = 111320 * math.cos(math.radians(cyd))
                    upmy = 110540
                else:
                    upmx = upmy = 1.0
                if margin > 0:                              # push vertices out m metres
                    cx = sum(q[0] for q in ring_r) / len(ring_r)
                    cy = sum(q[1] for q in ring_r) / len(ring_r)
                    grown = []
                    for x, y in ring_r:
                        dx, dy = (x - cx) * upmx, (y - cy) * upmy
                        d = math.hypot(dx, dy) or 1.0
                        grown.append([x + dx / d * margin / upmx,
                                      y + dy / d * margin / upmy])
                    ring_r = grown
                res = export_ortho_crop_from_raster(
                    ds, np, rasterio, Image, ring_r, pid, ortho_dir, ext, save_kw,
                    geotag, dict(PLOT_ID=str(pid)) if geotag else None,
                    max(2.0 * margin, 0.5), upmx, upmy)
                if res:
                    fn, cov, cw, ch2 = res
                    ok += 1
                    row.update(ortho_file=fn, crop_w=cw, crop_h=ch2,
                               ortho_covered_pct=f"{cov:.1f}",
                               ortho_covered_m2=f"{area * cov / 100:.3f}")
                    log(f"plot {pid} -> {fn}  {cw}x{ch2}px  cover {cov:.0f}%")
                else:
                    row['note'] = 'plot outside orthomosaic extent'
                    log(f"plot {pid}: outside the orthomosaic extent")
            except Exception as e:
                row['note'] = f'error: {e}'
                log(f"plot {pid}: ERROR {e}")
            manifest.append(row)
            set_progress(idx + 1, total)

    os.makedirs(p['outdir'], exist_ok=True)
    mpath = os.path.join(p['outdir'], "manifest.csv")
    with open(mpath, 'w', newline='', encoding='utf-8') as f:
        wcsv = csv.DictWriter(f, fieldnames=list(manifest[0].keys())
                              if manifest else ['plot_id'])
        wcsv.writeheader()
        wcsv.writerows(manifest)
    log("")
    log(f"DONE: ortho crops {ok}/{total}  ->  {p['outdir']}")
    log(f"Manifest: {mpath}")
    return p['outdir']


def make_tiled_ortho(src_path, dst_path, log, set_progress, cancel):
    """Rewrite a (typically striped) orthomosaic GeoTIFF as an internally TILED +
    overview'd copy so per-plot windowed reads become fast. Returns the output path,
    or the source path itself if it is already tiled with overviews. One-time cost."""
    import rasterio
    from rasterio.windows import Window
    from rasterio.enums import Resampling
    if not (src_path and os.path.exists(src_path)):
        raise RuntimeError("Select a valid GeoTIFF to tile.")
    with rasterio.open(src_path) as src:
        already = src.is_tiled and bool(src.overviews(1))
        if already:
            log(f"'{os.path.basename(src_path)}' is already tiled with overviews - "
                "no copy needed.")
            return src_path
        if dst_path is None:
            base, _ext = os.path.splitext(src_path)
            dst_path = base + "_tiled.tif"
        W, H, cnt = src.width, src.height, src.count
        prof = src.profile.copy()
        prof.update(driver='GTiff', tiled=True, blockxsize=512, blockysize=512,
                    compress='deflate', predictor=2, bigtiff='YES')
        prof.pop('photometric', None)
        band_rows = 2048                                   # copy this many rows at a time
        n = (H + band_rows - 1) // band_rows
        log(f"Tiling {os.path.basename(src_path)} -> {os.path.basename(dst_path)}  "
            f"({W}x{H}px, {cnt} bands): sequential copy then overviews ...")
        with rasterio.open(dst_path, 'w', NUM_THREADS='ALL_CPUS', **prof) as dst:
            for bi in range(n):
                if cancel.is_set():
                    log("Tiling cancelled by user.")
                    return None
                r0 = bi * band_rows
                rows = min(band_rows, H - r0)
                win = Window(0, r0, W, rows)
                dst.write(src.read(window=win), window=win)
                set_progress(bi + 1, n + 1)
                if bi % 5 == 0 or bi == n - 1:
                    log(f"  copied rows {r0 + rows}/{H} ({100*(r0+rows)//H}%)")
            log("  building internal overviews (2,4,8,16,32) ...")
            dst.build_overviews([2, 4, 8, 16, 32], Resampling.average)
            dst.update_tags(ns='rio_overview', resampling='average')
    set_progress(n + 1, n + 1)
    log(f"DONE: tiled copy written -> {dst_path}")
    return dst_path


# --------------------------------------------------------------------------
# Interactive plot editor (resize / reshape before extraction)
# --------------------------------------------------------------------------
def open_plot_editor(parent, plots_path, on_saved, log=lambda *_: None):
    """Interactive editor: wheel-zoom, pan (Space/middle-mouse/buttons), multi-
    select (click, Shift/Ctrl-click, rubber-band, Select all), move plots by
    dragging, reshape via corner handles, uniform grow/shrink, optional ortho
    basemap. Saves <name>_edited.shp and calls on_saved(new_path)."""
    import tkinter as tk
    from tkinter import ttk, messagebox, filedialog
    try:
        plots, _crs, fields = read_plot_file(plots_path)
    except Exception as e:
        messagebox.showerror("Edit plots", f"Could not read plots:\n{e}")
        return
    plots = [pl for pl in plots if len(pl['ring']) >= 3]
    if not plots:
        messagebox.showerror("Edit plots", "No polygons found in the file.")
        return
    orig = [[list(pt) for pt in pl['ring']] for pl in plots]
    rings = [[list(pt) for pt in pl['ring']] for pl in plots]
    id_field = pick_id_field(fields)
    ids = [safe_name(pl['props'].get(id_field, i + 1)) if id_field else str(i + 1)
           for i, pl in enumerate(plots)]
    allx = [p[0] for r in rings for p in r]
    ally = [p[1] for r in rings for p in r]
    geographic = max(abs(min(allx)), abs(max(allx)),
                     abs(min(ally)), abs(max(ally))) <= 360
    if geographic:
        cyc = sum(ally) / len(ally)
        mpx, mpy = 111320 * math.cos(math.radians(cyc)), 110540
    else:
        mpx = mpy = 1.0
    k = mpy / mpx

    win = tk.Toplevel(parent)
    win.title("Edit plots — zoom · pan · select · move · reshape · width/length/shift")
    win.geometry("1320x860")
    bar = ttk.Frame(win, padding=5)
    bar.pack(fill='x')
    cv = tk.Canvas(win, bg='#181818', highlightthickness=0)
    cv.pack(fill='both', expand=True)
    status = ttk.Label(win, text="", anchor='w', padding=4)
    status.pack(fill='x')

    T = {'scale': 1.0, 'ox': 0.0, 'oy': 0.0}
    sel = set()
    drag = {'kind': None, 'last': None, 'band': None, 'corner': None, 'resize': None}
    bm = {'img': None, 'bounds': None, 'photo': None, 'key': None}
    pan_mode = tk.BooleanVar(value=False)
    space_held = {'v': False}

    def _cfg_get(key, default):
        try:
            with open(CONFIG_PATH) as f:
                return json.load(f).get(key, default)
        except Exception:
            return default

    def _cfg_set(key, val):
        try:
            d = {}
            if os.path.exists(CONFIG_PATH):
                with open(CONFIG_PATH) as f:
                    d = json.load(f)
            d[key] = val
            with open(CONFIG_PATH, 'w') as f:
                json.dump(d, f, indent=1)
        except Exception:
            pass

    # Step default is read from config and remembered whenever it changes, so
    # the value you set becomes the default the next time the editor opens.
    step_var = tk.StringVar(value=str(_cfg_get('edit_step', '0.10')))

    def _step():
        try:
            return float(step_var.get())
        except ValueError:
            return 0.1

    def _persist_step(*_):
        try:
            float(step_var.get())          # only save valid numbers
        except ValueError:
            return
        _cfg_set('edit_step', step_var.get())
    step_var.trace_add('write', _persist_step)

    def w2s(x, y):
        return x * T['scale'] + T['ox'], T['oy'] - y * T['scale'] * k

    def s2w(sx, sy):
        return (sx - T['ox']) / T['scale'], (T['oy'] - sy) / (T['scale'] * k)

    def bounds():
        xs = [p[0] for r in rings for p in r]
        ys = [p[1] for r in rings for p in r]
        return min(xs), min(ys), max(xs), max(ys)

    def fit():
        x0, y0, x1, y1 = bounds()
        w = cv.winfo_width() or 1100
        h = cv.winfo_height() or 660
        pad = 50
        T['scale'] = min((w - 2 * pad) / max(x1 - x0, 1e-9),
                         (h - 2 * pad) / max((y1 - y0) * k, 1e-9))
        T['ox'] = pad - x0 * T['scale']
        T['oy'] = pad + y1 * T['scale'] * k

    def zoom_at(sx, sy, f):
        wx, wy = s2w(sx, sy)
        T['scale'] *= f
        T['ox'] = sx - wx * T['scale']
        T['oy'] = sy + wy * T['scale'] * k
        redraw()

    def redraw():
        cv.delete('all')
        if bm['img'] is not None and bm['bounds'] is not None:
            from PIL import Image as _I, ImageTk
            L, B, R, Tp = bm['bounds']
            sxL, syT = w2s(L, Tp)
            sxR, syB = w2s(R, B)
            fw, fh = sxR - sxL, syB - syT          # whole basemap, at this zoom
            iw, ih = bm['img'].size
            cw = cv.winfo_width() or 1100
            chv = cv.winfo_height() or 660
            # Render ONLY the part of the basemap inside the canvas, at canvas
            # resolution. Scaling the whole raster to the zoom level blows up
            # (a 10x zoom on a 5000px basemap = 50000px wide -> MemoryError).
            vx0, vy0 = max(sxL, 0.0), max(syT, 0.0)
            vx1, vy1 = min(sxR, float(cw)), min(syB, float(chv))
            if fw >= 1 and fh >= 1 and vx1 - vx0 >= 1 and vy1 - vy0 >= 1:
                box = (max(0, int(math.floor((vx0 - sxL) / fw * iw))),
                       max(0, int(math.floor((vy0 - syT) / fh * ih))),
                       min(iw, int(math.ceil((vx1 - sxL) / fw * iw))),
                       min(ih, int(math.ceil((vy1 - syT) / fh * ih))))
                if box[2] > box[0] and box[3] > box[1]:
                    # destination rect from the SNAPPED source box, so the
                    # basemap stays exactly registered to the plot outlines
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
        for i, r in enumerate(rings):
            pts = []
            for x, y in r:
                sx, sy = w2s(x, y)
                pts += [sx, sy]
            q = i in sel
            cv.create_polygon(pts, outline='#ffcc33' if q else '#5fae3f',
                              fill='', width=2 if q else 1)
            cxp = sum(pts[0::2]) / (len(pts) // 2)
            cyp = sum(pts[1::2]) / (len(pts) // 2)
            cv.create_text(cxp, cyp, text=ids[i],
                           fill='#e6e6e6' if q else '#9a9a9a', font=('Segoe UI', 8))
        if len(sel) == 1:
            for x, y in rings[next(iter(sel))]:
                sx, sy = w2s(x, y)
                cv.create_rectangle(sx - 4, sy - 4, sx + 4, sy + 4,
                                    outline='#ffcc33', fill='#3a2a00')
        elif len(sel) >= 2:
            hs = group_handles_screen()
            if hs:
                pts2 = []
                for hx, hy in hs:
                    pts2 += [hx, hy]
                cv.create_polygon(pts2, outline='#66ccff', fill='',
                                  dash=(5, 3), width=1)
                for hx, hy in hs:
                    cv.create_rectangle(hx - 5, hy - 5, hx + 5, hy + 5,
                                        outline='#66ccff', fill='#083048')
        if drag['kind'] == 'band' and drag['band']:
            x0, y0, x1, y1 = drag['band']
            cv.create_rectangle(x0, y0, x1, y1, outline='#66ccff', dash=(4, 3))
        status.config(text=status_text())

    def hit_corner(sx, sy):
        if len(sel) == 1:
            i = next(iter(sel))
            for vi, (x, y) in enumerate(rings[i]):
                px, py = w2s(x, y)
                if abs(px - sx) <= 6 and abs(py - sy) <= 6:
                    return (i, vi)
        return None

    def group_handles_screen():
        """Screen (x,y) of the 4 selection-bounding-box corners when 2+ plots are
        selected, else None. Dragging one scales every selected plot together."""
        if len(sel) < 2:
            return None
        xs = [p[0] for i in sel for p in rings[i]]
        ys = [p[1] for i in sel for p in rings[i]]
        x0, y0, x1, y1 = min(xs), min(ys), max(xs), max(ys)
        return [w2s(x0, y0), w2s(x1, y0), w2s(x1, y1), w2s(x0, y1)]

    def hit_group_handle(sx, sy):
        hs = group_handles_screen()
        if not hs:
            return None
        for gi, (hx, hy) in enumerate(hs):
            if abs(hx - sx) <= 7 and abs(hy - sy) <= 7:
                return gi
        return None

    def plot_at(wx, wy):
        for i in range(len(rings) - 1, -1, -1):
            if _point_in_ring(wx, wy, rings[i]):
                return i
        return None

    def on_press(e):
        cv.focus_set()
        if pan_mode.get() or space_held['v']:
            drag['kind'] = 'pan'
            drag['last'] = (e.x, e.y)
            return
        if not (e.state & 0x0005):               # not Shift/Ctrl -> allow group resize
            gh = hit_group_handle(e.x, e.y)
            if gh is not None:
                hs = group_handles_screen()
                gx = sum(h[0] for h in hs) / 4.0
                gy = sum(h[1] for h in hs) / 4.0
                drag['kind'] = 'resize'
                drag['resize'] = {
                    'gx': gx, 'gy': gy,
                    'd0': math.hypot(e.x - gx, e.y - gy) or 1.0,
                    'base': {i: [list(pt) for pt in rings[i]] for i in sel},
                    'cent': {i: (sum(p[0] for p in rings[i]) / len(rings[i]),
                                 sum(p[1] for p in rings[i]) / len(rings[i]))
                             for i in sel}}
                return
        c = hit_corner(e.x, e.y)
        if c is not None:
            drag['kind'] = 'corner'
            drag['corner'] = c
            return
        wx, wy = s2w(e.x, e.y)
        i = plot_at(wx, wy)
        multi = bool(e.state & 0x0005)          # Shift (0x1) or Ctrl (0x4)
        if i is not None:
            if multi:
                sel.symmetric_difference_update({i})
                drag['kind'] = None
                redraw()
                return
            if i not in sel:
                sel.clear()
                sel.add(i)
            drag['kind'] = 'move'
            drag['last'] = (wx, wy)
            redraw()
            return
        if not multi:
            sel.clear()
        drag['kind'] = 'band'
        drag['band'] = (e.x, e.y, e.x, e.y)
        redraw()

    def on_motion(e):
        kd = drag['kind']
        if kd == 'pan':
            lx, ly = drag['last']
            T['ox'] += e.x - lx
            T['oy'] += e.y - ly
            drag['last'] = (e.x, e.y)
            redraw()
        elif kd == 'corner':
            i, vi = drag['corner']
            rings[i][vi] = list(s2w(e.x, e.y))
            redraw()
        elif kd == 'resize':
            rz = drag['resize']
            f = max(0.05, math.hypot(e.x - rz['gx'], e.y - rz['gy']) / rz['d0'])
            for i in sel:
                cx, cy = rz['cent'][i]
                rings[i] = [[cx + (px - cx) * f, cy + (py - cy) * f]
                            for px, py in rz['base'][i]]
            redraw()
            status.config(text=f"Resizing {len(sel)} plots  ×{f:.2f}   "
                               "(drag the blue handle out = larger, in = smaller)")
        elif kd == 'move':
            wx, wy = s2w(e.x, e.y)
            lx, ly = drag['last']
            dx, dy = wx - lx, wy - ly
            for i in sel:
                for p in rings[i]:
                    p[0] += dx
                    p[1] += dy
            drag['last'] = (wx, wy)
            redraw()
        elif kd == 'band':
            x0, y0, _a, _b = drag['band']
            drag['band'] = (x0, y0, e.x, e.y)
            redraw()

    def on_release(_e):
        if drag['kind'] == 'band' and drag['band']:
            x0, y0, x1, y1 = drag['band']
            if abs(x1 - x0) > 3 and abs(y1 - y0) > 3:
                ax, ay = s2w(x0, y0)
                bx, by = s2w(x1, y1)
                xmn, xmx = min(ax, bx), max(ax, bx)
                ymn, ymx = min(ay, by), max(ay, by)
                for i, r in enumerate(rings):
                    cx = sum(p[0] for p in r) / len(r)
                    cy = sum(p[1] for p in r) / len(r)
                    if xmn <= cx <= xmx and ymn <= cy <= ymx:
                        sel.add(i)
        drag['kind'] = None
        drag['band'] = None
        redraw()

    def on_wheel(e):
        zoom_at(e.x, e.y, 1.2 if e.delta > 0 else 1 / 1.2)

    def on_mmb(e):
        drag['kind'] = 'pan'
        drag['last'] = (e.x, e.y)

    def on_mmb_move(e):
        if drag['kind'] == 'pan':
            lx, ly = drag['last']
            T['ox'] += e.x - lx
            T['oy'] += e.y - ly
            drag['last'] = (e.x, e.y)
            redraw()

    def on_mmb_up(_e):
        drag['kind'] = None

    def select_all():
        sel.clear()
        sel.update(range(len(rings)))
        redraw()

    def clear_sel():
        sel.clear()
        redraw()

    def _targets():
        """Plots to act on: the selection, or ALL plots when nothing is selected."""
        return sel if sel else set(range(len(rings)))

    def do_buffer(d):
        idxs = _targets()
        for i in idxs:
            r = rings[i]
            cx = sum(p[0] for p in r) / len(r)
            cy = sum(p[1] for p in r) / len(r)
            for p in r:
                dx, dy = (p[0] - cx) * mpx, (p[1] - cy) * mpy
                dd = math.hypot(dx, dy) or 1.0
                p[0] += dx / dd * d / mpx
                p[1] += dy / dd * d / mpy
        redraw()
        status.config(text=status_text(
            f"Buffer {'+' if d >= 0 else '−'}{abs(d):.2f} m -> "
            f"{'selection' if sel else 'ALL'} ({len(idxs)})"))

    def _axes_metric(r):
        """Principal axes of a plot in metric space (PCA): returns centroid and
        unit vectors (uL = length/long axis, uW = width/short axis)."""
        cx = sum(p[0] for p in r) / len(r)
        cy = sum(p[1] for p in r) / len(r)
        sxx = sxy = syy = 0.0
        for p in r:
            x = (p[0] - cx) * mpx
            y = (p[1] - cy) * mpy
            sxx += x * x; sxy += x * y; syy += y * y
        t = 0.5 * math.atan2(2 * sxy, sxx - syy)     # major-axis (max-variance) angle
        uL = (math.cos(t), math.sin(t))              # length  (long axis)
        uW = (-math.sin(t), math.cos(t))             # width   (short axis, perpendicular)
        return cx, cy, uL, uW

    def _wl(r):
        """(width, length) of a plot in metres: extent along its short/long axes."""
        cx, cy, uL, uW = _axes_metric(r)
        wc = [(p[0] - cx) * mpx * uW[0] + (p[1] - cy) * mpy * uW[1] for p in r]
        lc = [(p[0] - cx) * mpx * uL[0] + (p[1] - cy) * mpy * uL[1] for p in r]
        return max(wc) - min(wc), max(lc) - min(lc)

    def _wl_text():
        idxs = list(sel) if sel else list(range(len(rings)))
        if not idxs:
            return "W×L: n/a"
        ws, ls = [], []
        for i in idxs:
            w, l = _wl(rings[i]); ws.append(w); ls.append(l)
        if len(idxs) == 1:
            return f"W {ws[0]:.3f} × L {ls[0]:.3f} m"

        def rng(v):
            return (f"{v[0]:.2f}" if max(v) - min(v) < 0.005
                    else f"{min(v):.2f}–{max(v):.2f}")
        return (f"W {rng(ws)} × L {rng(ls)} m  "
                f"(mean {sum(ws)/len(ws):.2f}×{sum(ls)/len(ls):.2f})")

    def status_text(action=''):
        base = (f"{len(sel)} of {len(rings)} selected   |   {_wl_text()}   |   "
                f"zoom {T['scale']:.0f}   |   step {step_var.get()} m")
        if action:
            return f"{action}   |   {base}"
        return (base + "   |   wheel=zoom · Space/MMB=pan · click=select · "
                "Shift+click=add · drag plot=move · 1 sel: corner=reshape · "
                "2+ sel: blue handle=resize")

    def do_axis(which, d):
        """Change WIDTH (short axis) or LENGTH (long axis) of the target plots by
        d metres total (each of the two opposite sides moves d/2), keeping each
        plot's centre, orientation, and the other dimension unchanged."""
        idxs = _targets()
        axis_name = 'width' if which == 'width' else 'length'
        for i in idxs:
            r = rings[i]
            cx, cy, uL, uW = _axes_metric(r)
            ax = uW if which == 'width' else uL
            for p in r:
                x = (p[0] - cx) * mpx
                y = (p[1] - cy) * mpy
                comp = x * ax[0] + y * ax[1]                 # signed offset along axis (m)
                s = 1.0 if comp >= 0 else -1.0
                newmag = max(0.02, abs(comp) + d / 2.0)      # guard against collapse/inversion
                move = s * newmag - comp                     # metres to shift along the axis
                x += move * ax[0]; y += move * ax[1]
                p[0] = cx + x / mpx
                p[1] = cy + y / mpy
        redraw()
        status.config(text=status_text(
            f"{axis_name.capitalize()} {'+' if d >= 0 else '−'}{abs(d):.2f} m -> "
            f"{'selection' if sel else 'ALL'} ({len(idxs)})"))

    def do_shift(which, d):
        """Slide target plots by d metres along the WIDTH (cross-row) or LENGTH
        axis, keeping size/shape/orientation. Use it to move the boundary onto
        the planted rows when the plot sits off to one side. Direction is made
        consistent across plots so every plot moves the same way on screen."""
        idxs = _targets()
        for i in idxs:
            r = rings[i]
            _cx, _cy, uL, uW = _axes_metric(r)
            ax = uW if which == 'width' else uL
            if ax[0] < -1e-9 or (abs(ax[0]) <= 1e-9 and ax[1] < 0):   # +x (width) / +y (length)
                ax = (-ax[0], -ax[1])
            dxm, dym = d * ax[0], d * ax[1]
            for p in r:
                p[0] += dxm / mpx
                p[1] += dym / mpy
        redraw()
        arrow = {('width', True): '►', ('width', False): '◄',
                 ('length', True): '▲', ('length', False): '▼'}[(which, d >= 0)]
        status.config(text=status_text(
            f"Shift {arrow} {abs(d):.2f} m -> {'selection' if sel else 'ALL'} ({len(idxs)})"))

    def pan_by(dx, dy):
        T['ox'] += dx
        T['oy'] += dy
        redraw()

    def reset_all():
        for i in range(len(rings)):
            rings[i] = [list(pt) for pt in orig[i]]
        sel.clear()
        redraw()

    def load_basemap():
        pth = filedialog.askopenfilename(
            title="Ortho basemap (georeferenced GeoTIFF; any CRS - reprojected to the plots)",
            filetypes=[("GeoTIFF", "*.tif *.tiff")])
        if not pth:
            return
        try:
            import numpy as np
            import rasterio
            from rasterio.crs import CRS
            from rasterio.enums import Resampling
            from rasterio.vrt import WarpedVRT
            from rasterio.windows import from_bounds
            from PIL import Image as _I

            # --- plot CRS (the editor works in the shapefile's coordinates) ---
            plot_crs = None
            if _crs:
                for fn in (CRS.from_wkt, CRS.from_user_input):
                    try:
                        plot_crs = fn(_crs)
                        break
                    except Exception:
                        pass
            if plot_crs is None and geographic:
                plot_crs = CRS.from_epsg(4326)

            px0, px1 = min(allx), max(allx)
            py0, py1 = min(ally), max(ally)
            mx, my = 0.15 * max(px1 - px0, 1e-9), 0.15 * max(py1 - py0, 1e-9)

            with rasterio.open(pth) as ds:
                src_crs = ds.crs
                log(f"Basemap CRS: {src_crs.to_string() if src_crs else 'NONE'}; "
                    f"plots CRS: {plot_crs.to_string() if plot_crs else 'unknown'}")
                if plot_crs is not None and src_crs is not None and src_crs != plot_crs:
                    vrt = WarpedVRT(ds, crs=plot_crs, resampling=Resampling.bilinear)
                    log("Basemap reprojected to the plot CRS")
                else:
                    vrt = ds
                    if plot_crs is None or src_crs is None:
                        log("WARNING: could not compare CRSs - basemap used as-is")
                try:
                    b = vrt.bounds
                    # intersect raster with the plot extent (+15 % margin)
                    L = max(b.left, px0 - mx);   R = min(b.right, px1 + mx)
                    B = max(b.bottom, py0 - my); Tp = min(b.top, py1 + my)
                    if R <= L or Tp <= B:
                        raise RuntimeError(
                            "The basemap does not overlap the plots.\n"
                            f"Basemap extent: x {b.left:.2f}..{b.right:.2f}, "
                            f"y {b.bottom:.2f}..{b.top:.2f}\n"
                            f"Plot extent:    x {px0:.2f}..{px1:.2f}, "
                            f"y {py0:.2f}..{py1:.2f}\n"
                            "Check that the shapefile .prj and the GeoTIFF CRS are correct.")
                    win_ = from_bounds(L, B, R, Tp, vrt.transform)
                    ww, wh = win_.width, win_.height
                    dec = max(1.0, max(ww, wh) / 3000.0)
                    ow, oh = max(1, int(ww / dec)), max(1, int(wh / dec))
                    nb = vrt.count
                    bands = [1, 2, 3] if nb >= 3 else [1, 1, 1]
                    arr = vrt.read(bands, window=win_, out_shape=(3, oh, ow),
                                   resampling=Resampling.average).astype('float32')
                    alpha = None
                    if nb == 4 or nb == 2:
                        alpha = vrt.read(nb, window=win_, out_shape=(oh, ow),
                                         resampling=Resampling.nearest)
                    nod = vrt.nodata
                finally:
                    if vrt is not ds:
                        vrt.close()

            valid = np.ones(arr.shape[1:], bool)
            if nod is not None and not (isinstance(nod, float) and math.isnan(nod)):
                valid &= ~np.all(arr == nod, axis=0)
            elif nod is not None:
                valid &= ~np.any(np.isnan(arr), axis=0)
            if alpha is not None:
                valid &= alpha > 0
            arr = np.nan_to_num(arr)
            # scale to 8-bit: 8-bit data passes through, anything else gets a
            # 2-98 percentile stretch (16-bit orthos otherwise come out black)
            if arr[:, valid].size and (arr.max() > 255 or arr.max() <= 1.0):
                lo, hi = np.percentile(arr[:, valid], [2, 98])
                arr = (arr - lo) / max(hi - lo, 1e-9) * 255.0
            arr = np.clip(arr, 0, 255)
            arr[:, ~valid] = 24                       # nodata -> canvas grey
            bm['img'] = _I.fromarray(np.transpose(arr, (1, 2, 0)).astype('uint8'), 'RGB')
            bm['bounds'] = (L, B, R, Tp)
            bm['key'] = None
            log(f"Basemap loaded: {os.path.basename(pth)} ({ow}x{oh} px, "
                f"{100.0 * valid.mean():.0f}% valid data in the plot window)")
            redraw()
        except Exception as e:
            log("Basemap error: " + traceback.format_exc(limit=3))
            messagebox.showerror("Basemap", f"Could not load basemap:\n{e}")

    def save():
        import shapefile
        import shutil
        base = os.path.splitext(plots_path)[0]
        outp = base + "_edited.shp"
        try:
            w = shapefile.Writer(outp, shapeType=shapefile.POLYGON)
            recs = None
            try:
                rdr = shapefile.Reader(plots_path)
                for f in rdr.fields[1:]:
                    w.field(*f)
                recs = list(rdr.records())
            except Exception:
                w.field('Plot_ID', 'C', 40)
            for i, r in enumerate(rings):
                w.poly([[list(p) for p in (r + [r[0]])]])
                if recs is not None and i < len(recs):
                    w.record(*recs[i])
                else:
                    w.record(ids[i])
            w.close()
            prj = base + '.prj'
            if os.path.exists(prj):
                shutil.copyfile(prj, outp[:-4] + '.prj')
        except Exception as e:
            messagebox.showerror("Edit plots", f"Could not save:\n{e}")
            return
        log(f"Edited plots saved: {outp}")
        on_saved(outp)
        messagebox.showinfo("Edit plots",
                            f"Saved adjusted plots:\n{outp}\n\n"
                            "The extractor will now use this edited file.")
        win.destroy()

    # ---- toolbar (each tool shows its keyboard shortcut underneath) ----
    def tool(label, short, cmd, w=11):
        f = ttk.Frame(bar)
        f.pack(side='left', padx=3)
        ttk.Button(f, text=label, command=cmd, width=w).pack()
        ttk.Label(f, text=short, foreground='#888', font=('Segoe UI', 7)).pack()

    tool("Select all", "Ctrl+A", select_all)
    tool("Clear", "Esc", clear_sel, w=7)
    pf = ttk.Frame(bar)
    pf.pack(side='left', padx=3)
    ttk.Checkbutton(pf, text="Pan mode", variable=pan_mode).pack()
    ttk.Label(pf, text="Space / MMB", foreground='#888', font=('Segoe UI', 7)).pack()
    tool("Zoom in", "wheel up", lambda: zoom_at(cv.winfo_width() / 2, cv.winfo_height() / 2, 1.2), w=8)
    tool("Zoom out", "wheel down", lambda: zoom_at(cv.winfo_width() / 2, cv.winfo_height() / 2, 1 / 1.2), w=8)
    tool("Fit", "F", lambda: (fit(), redraw()), w=5)
    sf = ttk.Frame(bar)
    sf.pack(side='left', padx=(10, 3))
    rowf = ttk.Frame(sf)
    rowf.pack()
    ttk.Label(rowf, text="Step m:").pack(side='left')
    ttk.Entry(rowf, textvariable=step_var, width=5).pack(side='left')
    ttk.Label(sf, text="selection · else ALL plots", foreground='#888',
              font=('Segoe UI', 7)).pack()
    tool("Grow +", "+ key", lambda: do_buffer(_step()), w=8)
    tool("Shrink −", "− key", lambda: do_buffer(-_step()), w=8)
    tool("Width +", "] key", lambda: do_axis('width', _step()), w=8)
    tool("Width −", "[ key", lambda: do_axis('width', -_step()), w=8)
    tool("Length +", "' key", lambda: do_axis('length', _step()), w=8)
    tool("Length −", "; key", lambda: do_axis('length', -_step()), w=8)
    tool("Shift ◄", "Shift+←", lambda: do_shift('width', -_step()), w=6)
    tool("Shift ►", "Shift+→", lambda: do_shift('width', _step()), w=6)
    tool("Shift ▲", "Shift+↑", lambda: do_shift('length', _step()), w=6)
    tool("Shift ▼", "Shift+↓", lambda: do_shift('length', -_step()), w=6)
    tool("Basemap", "ortho .tif", load_basemap, w=9)
    tool("Reset all", "Ctrl+Z", reset_all, w=9)
    tool("Save & use", "S", save, w=10)

    # ---- bindings ----
    cv.bind('<ButtonPress-1>', on_press)
    cv.bind('<B1-Motion>', on_motion)
    cv.bind('<ButtonRelease-1>', on_release)
    cv.bind('<ButtonPress-2>', on_mmb)
    cv.bind('<B2-Motion>', on_mmb_move)
    cv.bind('<ButtonRelease-2>', on_mmb_up)
    cv.bind('<MouseWheel>', on_wheel)
    win.bind('<Control-a>', lambda e: select_all())
    win.bind('<Control-A>', lambda e: select_all())
    win.bind('<Escape>', lambda e: clear_sel())
    win.bind('<Control-z>', lambda e: reset_all())
    win.bind('<KeyPress-space>', lambda e: space_held.__setitem__('v', True))
    win.bind('<KeyRelease-space>', lambda e: space_held.__setitem__('v', False))
    for kk in ('<KeyPress-f>', '<KeyPress-F>'):
        win.bind(kk, lambda e: (fit(), redraw()))
    for kk in ('<KeyPress-s>', '<KeyPress-S>'):
        win.bind(kk, lambda e: save())
    for kk in ('+', '=', '<KP_Add>'):
        win.bind(kk, lambda e: do_buffer(_step()))
    for kk in ('-', '<KP_Subtract>'):
        win.bind(kk, lambda e: do_buffer(-_step()))
    win.bind(']', lambda e: do_axis('width', _step()))       # width wider
    win.bind('[', lambda e: do_axis('width', -_step()))      # width narrower
    win.bind("'", lambda e: do_axis('length', _step()))      # length longer
    win.bind(';', lambda e: do_axis('length', -_step()))     # length shorter
    win.bind('<Shift-Left>', lambda e: do_shift('width', -_step()))    # slide onto rows
    win.bind('<Shift-Right>', lambda e: do_shift('width', _step()))
    win.bind('<Shift-Up>', lambda e: do_shift('length', _step()))
    win.bind('<Shift-Down>', lambda e: do_shift('length', -_step()))
    win.bind('<Left>', lambda e: pan_by(60, 0))
    win.bind('<Right>', lambda e: pan_by(-60, 0))
    win.bind('<Up>', lambda e: pan_by(0, 60))
    win.bind('<Down>', lambda e: pan_by(0, -60))
    cv.bind('<Configure>', lambda e: redraw())
    win.after(80, lambda: (fit(), redraw()))


# --------------------------------------------------------------------------
# GUI
# --------------------------------------------------------------------------
def main():
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox

    cfg = {}
    if os.path.exists(CONFIG_PATH):
        try:
            with open(CONFIG_PATH) as f:
                cfg = json.load(f)
        except Exception:
            cfg = {}

    root = tk.Tk()
    root.title("Plot Extractor — raw image per plot")
    root.geometry("860x720")
    root.minsize(760, 620)

    # ---- header banner with partner logos (DPIRD left, APPN right) --------
    def _asset(name):
        base = getattr(sys, '_MEIPASS',
                       os.path.dirname(os.path.abspath(__file__)))
        return os.path.join(base, 'assets', name)
    try:
        from PIL import Image as _HImg, ImageTk as _HTk
        header = tk.Frame(root, bg='white')
        header.pack(side='top', fill='x')
        root._logo_refs = []

        def _logo(name, target_h=52):
            im = _HImg.open(_asset(name)).convert('RGBA')
            w, h = im.size
            im = im.resize((max(1, round(w * target_h / h)), target_h),
                           _HImg.LANCZOS)
            ph = _HTk.PhotoImage(im)
            root._logo_refs.append(ph)      # keep ref so tk doesn't GC it
            return ph

        tk.Label(header, image=_logo('logo_dpird.png'), bg='white'
                 ).pack(side='left', padx=14, pady=8)
        tk.Label(header, image=_logo('logo_appn.png'), bg='white'
                 ).pack(side='right', padx=14, pady=8)
        tk.Label(header, text="Plot Extractor", bg='white', fg='#0c1a2b',
                 font=('Segoe UI', 14, 'bold')).pack(side='top', pady=16)
        ttk.Separator(root, orient='horizontal').pack(side='top', fill='x')
    except Exception:
        pass    # logos are cosmetic — never block the app

    # Tabs: the extractor form, and the plot-grid designer that builds the
    # shapefile the extractor consumes (ortho basemap, no Metashape needed).
    nb = ttk.Notebook(root)
    nb.pack(fill='both', expand=True)
    tab_extract = ttk.Frame(nb)
    tab_grid = ttk.Frame(nb)
    tab_tiles = ttk.Frame(nb)
    nb.add(tab_extract, text="  Extract plots  ")
    nb.add(tab_grid, text="  Generate plot grid  ")
    nb.add(tab_tiles, text="  Training tiles  ")

    frm = ttk.Frame(tab_extract, padding=10)
    frm.pack(fill='both', expand=True)
    frm.columnconfigure(1, weight=1)

    vars_ = {
        'source': tk.StringVar(value=cfg.get('source', 'Metashape project')),
        'psx': tk.StringVar(value=cfg.get('psx', '')),
        'chunk': tk.StringVar(value=cfg.get('chunk', '(active)')),
        'ortho_tif': tk.StringVar(value=cfg.get('ortho_tif', '')),
        'plots': tk.StringVar(value=cfg.get('plots', '')),
        'id_field': tk.StringVar(value=cfg.get('id_field', AUTO_ID)),
        'rawdir': tk.StringVar(value=cfg.get('rawdir', '')),
        'outdir': tk.StringVar(value=cfg.get('outdir', '')),
        'margin': tk.StringVar(value=str(cfg.get('margin', 0.2))),
        'n_per_plot': tk.StringVar(value=str(cfg.get('n_per_plot', 1))),
        'limit': tk.StringVar(value=str(cfg.get('limit', 0))),
        'quality': tk.StringVar(value=str(cfg.get('quality', 92))),
        'do_raw': tk.BooleanVar(value=cfg.get('do_raw', True)),
        'raw_mode': tk.StringVar(value=cfg.get('raw_mode', 'Native (no resample, no black)')),
        'engine': tk.StringVar(value=cfg.get('engine', 'Metashape API')),
        'dem_tif': tk.StringVar(value=cfg.get('dem_tif', '')),
        'do_ortho': tk.BooleanVar(value=cfg.get('do_ortho', False)),
        'crs_override': tk.StringVar(value=cfg.get('crs_override', 'auto')),
        'gcp_epsg': tk.StringVar(value=str(cfg.get('gcp_epsg', ''))),
        'format': tk.StringVar(value=cfg.get('format', 'TIFF (lossless)')),
        'to_srgb': tk.BooleanVar(value=cfg.get('to_srgb', False)),
        'workers': tk.StringVar(value=str(cfg.get('workers', 0))),
        'auto_fit': tk.BooleanVar(value=cfg.get('auto_fit', False)),
        'rows_per_plot': tk.StringVar(value=str(cfg.get('rows_per_plot', 0))),
        'fit_mode': tk.StringVar(value=cfg.get('fit_mode', 'width')),
        'smart_seam': tk.BooleanVar(value=cfg.get('smart_seam', True)),
        'seam_map': tk.BooleanVar(value=cfg.get('seam_map', False)),
        'mesh_step': tk.StringVar(value=str(cfg.get('mesh_step', MESH_STEP_PX))),
    }

    r = 0
    def add_path_row(label, var, kind, filetypes=None, hint=None):
        nonlocal r
        ttk.Label(frm, text=label).grid(row=r, column=0, sticky='w', pady=3)
        ent = ttk.Entry(frm, textvariable=var)
        ent.grid(row=r, column=1, sticky='ew', padx=6)
        def browse():
            if kind == 'file':
                v = filedialog.askopenfilename(filetypes=filetypes)
            else:
                v = filedialog.askdirectory()
            if v:
                var.set(v)
        ttk.Button(frm, text="Browse...", command=browse).grid(row=r, column=2)
        r += 1
        if hint:
            ttk.Label(frm, text=hint, foreground='grey').grid(
                row=r, column=1, sticky='w', padx=6)
            r += 1
        return ent

    ttk.Label(frm, text="Source").grid(row=r, column=0, sticky='w', pady=3)
    src_combo = ttk.Combobox(frm, textvariable=vars_['source'], state='readonly',
                             values=["Metashape project",
                                     "Orthomosaic GeoTIFF (no Metashape)"])
    src_combo.grid(row=r, column=1, sticky='w', padx=6)
    r += 1

    psx_ent = add_path_row("Metashape project (.psx) *", vars_['psx'], 'file',
                           [("Metashape project", "*.psx")])

    ttk.Label(frm, text="Chunk / sub-project").grid(row=r, column=0, sticky='w', pady=3)
    chunk_combo = ttk.Combobox(frm, textvariable=vars_['chunk'], state='readonly',
                               values=["(active)"], width=32)
    chunk_combo.grid(row=r, column=1, sticky='w', padx=6)
    r += 1
    ttk.Label(frm, text="Pick the right chunk when a project holds several "
                        "(e.g. RootPheno_Day1 / RootPheno_Day2)",
              foreground='grey').grid(row=r, column=1, sticky='w', padx=6)
    r += 1

    def populate_chunks(*_):
        pth = vars_['psx'].get().strip()
        if not (pth and os.path.exists(pth)):
            chunk_combo['values'] = ["(active)"]
            return
        labels, active = read_psx_chunks(pth)
        vals = ["(active)"] + labels
        chunk_combo['values'] = vals
        if vars_['chunk'].get() not in vals:
            vars_['chunk'].set("(active)")
    vars_['psx'].trace_add('write', lambda *_: root.after(200, populate_chunks))
    populate_chunks()
    ot_ent = add_path_row("Orthomosaic GeoTIFF *", vars_['ortho_tif'], 'file',
                          [("GeoTIFF", "*.tif *.tiff")],
                          hint="Ortho-only mode: crops plots straight from this raster, "
                               "no Metashape / license needed")
    add_path_row("Plot boundaries (.shp / .geojson) *", vars_['plots'], 'file',
                 [("Plot boundaries", "*.shp *.geojson *.json")])

    idrow = ttk.Frame(frm)
    idrow.grid(row=r, column=1, sticky='w', padx=6)
    ttk.Label(frm, text="Plot ID attribute").grid(row=r, column=0, sticky='w', pady=3)
    id_combo = ttk.Combobox(idrow, textvariable=vars_['id_field'],
                            values=[AUTO_ID], state='readonly', width=30)
    id_combo.pack(side='left')
    ttk.Label(idrow, text="   Plot CRS:").pack(side='left')
    ttk.Entry(idrow, textvariable=vars_['crs_override'], width=10).pack(side='left', padx=4)
    ttk.Label(idrow, text="(auto = read from file; else EPSG code)",
              foreground='grey').pack(side='left')
    r += 1

    # CRS the GeoTIFF ground control points are written in. Geometry is
    # normalised to EPSG:4326 internally, so without this every crop ships in
    # degrees regardless of the grid's own CRS - which is how one delivery ended
    # up half in 4326 and half in 7850. Coordinates are transformed, not relabelled.
    crsrow = ttk.Frame(frm)
    crsrow.grid(row=r, column=1, sticky='w', padx=6)
    ttk.Label(frm, text="Output GCP CRS").grid(row=r, column=0, sticky='w', pady=3)
    ttk.Entry(crsrow, textvariable=vars_['gcp_epsg'], width=10).pack(side='left')
    ttk.Label(crsrow, text="EPSG code for the GeoTIFF ground control points "
                           "(blank = 4326). Set it once for the whole project.",
              foreground='grey').pack(side='left', padx=4)
    r += 1

    add_path_row("Raw images folder (optional)", vars_['rawdir'], 'dir',
                 hint="Only needed if the .psx points to photos that were moved")
    dem_ent = add_path_row("DEM / DSM GeoTIFF (EasyIDP engine)", vars_['dem_tif'], 'file',
                           [("GeoTIFF", "*.tif *.tiff")],
                           hint="Used by the EasyIDP engine to lift plots to 3D. "
                                "Leave blank to auto-export the DEM from the .psx")
    add_path_row("Output folder *", vars_['outdir'], 'dir')

    opts = ttk.Frame(frm)
    opts.grid(row=r, column=0, columnspan=3, sticky='w', pady=6)
    r += 1
    for text, var, width in [("Margin around plot (m):", vars_['margin'], 6),
                             ("Images per plot:", vars_['n_per_plot'], 4),
                             ("Limit (0 = all plots):", vars_['limit'], 6)]:
        ttk.Label(opts, text=text).pack(side='left', padx=(0, 4))
        ttk.Entry(opts, textvariable=var, width=width).pack(side='left', padx=(0, 14))
    ttk.Label(opts, text="Format:").pack(side='left', padx=(0, 4))
    fmt_combo = ttk.Combobox(opts, textvariable=vars_['format'], width=15, state='readonly',
                             values=["GeoTIFF (geotagged, lossless)", "PNG (lossless)",
                                     "TIFF (lossless)", "JPEG (small files)"])
    fmt_combo.pack(side='left', padx=(0, 14))
    quality_lbl = ttk.Label(opts, text="JPEG quality:")
    quality_ent = ttk.Entry(opts, textvariable=vars_['quality'], width=4)

    def toggle_quality(*_):
        if vars_['format'].get().startswith('JPEG'):
            quality_lbl.pack(side='left', padx=(0, 4))
            quality_ent.pack(side='left')
        else:
            quality_lbl.pack_forget()
            quality_ent.pack_forget()
    vars_['format'].trace_add('write', toggle_quality)
    toggle_quality()

    checks = ttk.Frame(frm)
    checks.grid(row=r, column=0, columnspan=3, sticky='w')
    r += 1
    raw_chk = ttk.Checkbutton(checks, text="Raw image crop per plot",
                              variable=vars_['do_raw'])
    raw_chk.pack(side='left', padx=(0, 6))
    ttk.Label(checks, text="style:").pack(side='left', padx=(0, 3))
    raw_mode_combo = ttk.Combobox(checks, textvariable=vars_['raw_mode'], width=30,
                                  state='readonly',
                                  values=["Native (no resample, no black)",
                                          "Auto-straighten tilted", "Bounding box",
                                          "Masked to plot", "Rectified rectangle",
                                          "Best single frame (no stitching)"])
    raw_mode_combo.pack(side='left', padx=(0, 12))
    ttk.Label(checks, text="engine:").pack(side='left', padx=(0, 3))
    engine_combo = ttk.Combobox(checks, textvariable=vars_['engine'], width=14,
                                state='readonly',
                                values=["Metashape API", "EasyIDP"])
    engine_combo.pack(side='left', padx=(0, 16))
    ortho_chk = ttk.Checkbutton(checks, text="Orthomosaic crop per plot (deskewed)",
                                variable=vars_['do_ortho'])
    ortho_chk.pack(side='left')

    checks2 = ttk.Frame(frm)
    checks2.grid(row=r, column=0, columnspan=3, sticky='w', pady=(4, 0))
    r += 1
    srgb_chk = ttk.Checkbutton(
        checks2, text="Convert raw crops to sRGB (from Adobe RGB / embedded ICC)",
        variable=vars_['to_srgb'])
    srgb_chk.pack(side='left')
    _ncores = os.cpu_count() or 4
    ttk.Label(checks2, text="CPU cores:").pack(side='left', padx=(18, 3))
    ttk.Spinbox(checks2, from_=0, to=_ncores, width=4,
                textvariable=vars_['workers']).pack(side='left')
    ttk.Label(checks2, text=f"(0 = auto, max {_ncores})",
              foreground='#888').pack(side='left', padx=(4, 0))

    checks3 = ttk.Frame(frm)
    checks3.grid(row=r, column=0, columnspan=3, sticky='w', pady=(4, 0))
    r += 1
    fit_chk = ttk.Checkbutton(
        checks3, text="Auto-fit crop to detected rows",
        variable=vars_['auto_fit'])
    fit_chk.pack(side='left')
    ttk.Label(checks3, text="rows/plot:").pack(side='left', padx=(14, 3))
    ttk.Spinbox(checks3, from_=0, to=60, width=4,
                textvariable=vars_['rows_per_plot']).pack(side='left')
    ttk.Label(checks3, text="(0 = all found)", foreground='#888').pack(
        side='left', padx=(3, 0))
    ttk.Label(checks3, text="mode:").pack(side='left', padx=(14, 3))
    ttk.Combobox(checks3, textvariable=vars_['fit_mode'], width=14,
                 state='readonly',
                 values=["width", "width_shift", "width_length"]).pack(side='left')
    ttk.Label(checks3, text="— recovers rows the grid clipped by rendering wider",
              foreground='#888').pack(side='left', padx=(8, 0))

    checks4 = ttk.Frame(frm)
    checks4.grid(row=r, column=0, columnspan=3, sticky='w', pady=(4, 0))
    r += 1
    seam_chk = ttk.Checkbutton(
        checks4, text="Seamless stitching (min-error seam + feathered colour match)",
        variable=vars_['smart_seam'])
    seam_chk.pack(side='left')
    ttk.Label(checks4, text="rectify grid:").pack(side='left', padx=(18, 3))
    mesh_combo = ttk.Combobox(checks4, textvariable=vars_['mesh_step'], width=5,
                              state='readonly', values=["48", "96", "192", "256"])
    mesh_combo.pack(side='left')
    ttk.Label(checks4, text="px cells (smaller = more exact, slower)",
              foreground='#888').pack(side='left', padx=(4, 0))
    seammap_chk = ttk.Checkbutton(checks4, text="save seam map (QC)",
                                  variable=vars_['seam_map'])
    seammap_chk.pack(side='left', padx=(18, 0))

    def toggle_source(*_):
        ortho_only = vars_['source'].get().startswith('Orthomosaic')
        psx_ent.configure(state='disabled' if ortho_only else 'normal')
        ot_ent.configure(state='normal' if ortho_only else 'disabled')
        raw_chk.configure(state='disabled' if ortho_only else 'normal')
        raw_mode_combo.configure(state='disabled' if ortho_only else 'readonly')
        srgb_chk.configure(state='disabled' if ortho_only else 'normal')
        for wdg in (seam_chk, mesh_combo, seammap_chk):
            wdg.configure(state='disabled' if ortho_only else
                          ('readonly' if wdg is mesh_combo else 'normal'))
        if ortho_only:
            vars_['do_ortho'].set(True)
    vars_['source'].trace_add('write', toggle_source)

    def toggle_seam(*_):
        # stitching controls only mean something when frames are being stitched
        stitching = vars_['raw_mode'].get().startswith(('Rectified', 'Auto-straighten'))
        if not vars_['source'].get().startswith('Orthomosaic'):
            seam_chk.configure(state='normal' if stitching else 'disabled')
            seammap_chk.configure(state='normal' if stitching else 'disabled')
    vars_['raw_mode'].trace_add('write', toggle_seam)
    toggle_source()
    toggle_seam()

    btns = ttk.Frame(frm)
    btns.grid(row=r, column=0, columnspan=3, sticky='ew', pady=8)
    r += 1
    run_btn = ttk.Button(btns, text="Run extraction")
    run_btn.pack(side='left')
    cancel_btn = ttk.Button(btns, text="Cancel", state='disabled')
    cancel_btn.pack(side='left', padx=8)
    open_btn = ttk.Button(btns, text="Open output folder", state='disabled')
    open_btn.pack(side='left', padx=8)
    edit_btn = ttk.Button(btns, text="Edit plots…")
    edit_btn.pack(side='left', padx=8)
    tile_btn = ttk.Button(btns, text="Make tiled copy")
    tile_btn.pack(side='left', padx=8)
    fixgrid_btn = ttk.Button(btns, text="Correct grid to rows…")
    fixgrid_btn.pack(side='left', padx=8)
    setcrs_btn = ttk.Button(btns, text="Set CRS of existing crops…")
    setcrs_btn.pack(side='left', padx=8)
    prog = ttk.Progressbar(btns, mode='determinate')
    prog.pack(side='left', fill='x', expand=True, padx=8)

    logbox = tk.Text(frm, height=14, wrap='none', state='disabled',
                     font=('Consolas', 9))
    logbox.grid(row=r, column=0, columnspan=3, sticky='nsew')
    frm.rowconfigure(r, weight=1)
    sb = ttk.Scrollbar(frm, command=logbox.yview)
    sb.grid(row=r, column=3, sticky='ns')
    logbox.configure(yscrollcommand=sb.set)

    msg_q = queue.Queue()
    cancel_ev = threading.Event()
    state = {'running': False, 'outdir': None}

    def log(msg):
        msg_q.put(('log', str(msg)))

    def set_progress(done, total):
        msg_q.put(('prog', (done, total)))

    def refresh_id_fields(*_):
        path = vars_['plots'].get().strip()
        if not path or not os.path.exists(path):
            return
        try:
            _plots, _wkt, fields = read_plot_file(path)
            vals = [AUTO_ID] + fields
            id_combo['values'] = vals
            best = pick_id_field(fields)
            if vars_['id_field'].get() not in vals:
                vars_['id_field'].set(best or AUTO_ID)
            elif vars_['id_field'].get() == AUTO_ID and best:
                vars_['id_field'].set(best)
            log(f"Boundary file: {len(_plots)} plots, fields: {', '.join(fields) or '-'}")
            if not vars_['outdir'].get().strip():
                vars_['outdir'].set(os.path.join(os.path.dirname(path), "PlotImages"))
        except Exception as e:
            log(f"[!] Could not read boundary file: {e}")

    vars_['plots'].trace_add('write', lambda *_: root.after(300, refresh_id_fields))

    def save_config():
        out = {k: (v.get() if not isinstance(v, tk.BooleanVar) else bool(v.get()))
               for k, v in vars_.items()}
        try:
            with open(CONFIG_PATH, 'w') as f:
                json.dump(out, f, indent=1)
        except OSError:
            pass

    def worker(params):
        try:
            if params.get('source') == 'ortho_tif':
                outdir = run_extraction_ortho_tif(params, log, set_progress, cancel_ev)
            else:
                outdir = run_extraction(params, log, set_progress, cancel_ev)
            msg_q.put(('done', outdir))
        except Exception:
            log("ERROR:\n" + traceback.format_exc(limit=3))
            msg_q.put(('done', None))

    def on_run():
        ortho_only = vars_['source'].get().startswith('Orthomosaic')
        psx = vars_['psx'].get().strip()
        ortho_tif = vars_['ortho_tif'].get().strip()
        plots = vars_['plots'].get().strip()
        outdir = vars_['outdir'].get().strip()
        if ortho_only:
            if not (ortho_tif and os.path.exists(ortho_tif)):
                messagebox.showerror("Plot Extractor",
                                     "Please select a valid orthomosaic GeoTIFF.")
                return
        elif not (psx and os.path.exists(psx)):
            messagebox.showerror("Plot Extractor", "Please select a valid .psx project.")
            return
        if not (plots and os.path.exists(plots)):
            messagebox.showerror("Plot Extractor", "Please select a valid plot boundary file.")
            return
        if not outdir:
            messagebox.showerror("Plot Extractor", "Please choose an output folder.")
            return
        if not ortho_only and not (vars_['do_raw'].get() or vars_['do_ortho'].get()):
            messagebox.showerror("Plot Extractor", "Select at least one output type.")
            return
        try:
            if ortho_only:
                params = dict(source='ortho_tif', ortho_tif=ortho_tif, plots=plots,
                              outdir=outdir, id_field=vars_['id_field'].get(),
                              margin=float(vars_['margin'].get() or 0),
                              limit=int(vars_['limit'].get() or 0),
                              quality=int(vars_['quality'].get() or 92),
                              crs_override=vars_['crs_override'].get(),
                              gcp_epsg=vars_['gcp_epsg'].get().strip(),
                              format=vars_['format'].get())
            else:
                params = dict(source='metashape', psx=psx, plots=plots, outdir=outdir,
                              chunk=vars_['chunk'].get().strip(),
                              rawdir=vars_['rawdir'].get().strip(),
                              id_field=vars_['id_field'].get(),
                              margin=float(vars_['margin'].get() or 0),
                              n_per_plot=int(vars_['n_per_plot'].get() or 1),
                              limit=int(vars_['limit'].get() or 0),
                              quality=int(vars_['quality'].get() or 92),
                              do_raw=vars_['do_raw'].get(),
                              raw_mode={"Native (no resample, no black)": "native",
                                        "Auto-straighten tilted": "auto",
                                        "Bounding box": "bbox", "Masked to plot": "mask",
                                        "Rectified rectangle": "rectify",
                                        "Best single frame (no stitching)": "best"}.get(
                                            vars_['raw_mode'].get(), "native"),
                              smart_seam=vars_['smart_seam'].get(),
                              seam_map=vars_['seam_map'].get(),
                              mesh_step=int(vars_['mesh_step'].get() or 96),
                              engine={"Metashape API": "metashape",
                                      "EasyIDP": "easyidp"}.get(
                                          vars_['engine'].get(), "metashape"),
                              dem_tif=vars_['dem_tif'].get().strip(),
                              do_ortho=vars_['do_ortho'].get(),
                              to_srgb=vars_['to_srgb'].get(),
                              workers=int(vars_['workers'].get() or 0),
                              auto_fit=vars_['auto_fit'].get(),
                              rows_per_plot=int(vars_['rows_per_plot'].get() or 0),
                              fit_mode=vars_['fit_mode'].get(),
                              crs_override=vars_['crs_override'].get(),
                              gcp_epsg=vars_['gcp_epsg'].get().strip(),
                              format=vars_['format'].get())
        except ValueError:
            messagebox.showerror("Plot Extractor",
                                 "Margin / images per plot / limit / quality must be numbers.")
            return
        save_config()
        cancel_ev.clear()
        state['running'] = True
        state['outdir'] = outdir
        run_btn.configure(state='disabled')
        cancel_btn.configure(state='normal')
        open_btn.configure(state='disabled')
        prog['value'] = 0
        log("=" * 60)
        threading.Thread(target=worker, args=(params,), daemon=True).start()

    def on_cancel():
        cancel_ev.set()
        log("Cancelling after current plot ...")

    def on_open():
        if state['outdir'] and os.path.isdir(state['outdir']):
            os.startfile(state['outdir'])

    def on_edit_plots():
        pth = vars_['plots'].get().strip()
        if not (pth and os.path.exists(pth)):
            messagebox.showerror("Plot Extractor",
                                 "Select a plot boundary file first.")
            return
        open_plot_editor(root, pth, lambda new: vars_['plots'].set(new), log)

    def on_make_tiled():
        if state['running']:
            return
        src = vars_['ortho_tif'].get().strip()
        if not (src and os.path.exists(src)):
            src = filedialog.askopenfilename(
                title="Orthomosaic GeoTIFF to convert to a tiled copy",
                filetypes=[("GeoTIFF", "*.tif *.tiff")])
            if not src:
                return
            vars_['ortho_tif'].set(src)
        cancel_ev.clear()
        state['running'] = True
        run_btn.configure(state='disabled')
        tile_btn.configure(state='disabled')
        cancel_btn.configure(state='normal')
        prog['value'] = 0
        log("=" * 60)

        def tworker():
            try:
                dst = make_tiled_ortho(src, None, log, set_progress, cancel_ev)
                msg_q.put(('tiled_done', dst))
            except Exception:
                log("ERROR:\n" + traceback.format_exc(limit=3))
                msg_q.put(('tiled_done', None))
        threading.Thread(target=tworker, args=(), daemon=True).start()

    run_btn.configure(command=on_run)
    cancel_btn.configure(command=on_cancel)
    open_btn.configure(command=on_open)
    def on_fix_grid():
        """Centre the plot grid on the rows that were actually drilled.

        Measures every plot against the orthomosaic and writes <grid>_refit.shp
        beside the original, with an audit CSV. It refuses when the trial is
        drilled continuously (row pitch == plot width / rows), because then a plot
        boundary leaves no signature in the imagery and 'centring on the rows'
        would just be moving polygons on detection noise.
        """
        shp = vars_['plots'].get().strip()
        if not shp:
            messagebox.showwarning("Correct grid",
                                   "Choose the plot boundaries first.")
            return
        ortho = vars_['ortho_tif'].get().strip()
        if not ortho or not os.path.exists(ortho):
            ortho = filedialog.askopenfilename(
                title="Orthomosaic covering the trial",
                filetypes=[("GeoTIFF", "*.tif *.tiff"), ("All files", "*.*")])
            if not ortho:
                return
            vars_['ortho_tif'].set(ortho)
        try:
            rows = int(vars_['rows_per_plot'].get() or 0)
        except (TypeError, ValueError):
            rows = 0
        out = os.path.splitext(shp)[0] + "_refit.shp"
        fixgrid_btn.configure(state='disabled')
        log("")
        log("Correcting the plot grid against the orthomosaic — this reads the "
            "whole trial, allow a few minutes.")

        def work():
            import io as _io, contextlib as _ctx
            code, buf = 1, _io.StringIO()
            try:
                import grid_refit
                argv = ["--shp", shp, "--ortho", ortho, "--out", out]
                if rows:
                    argv += ["--rows", str(rows)]
                with _ctx.redirect_stdout(buf):
                    code = grid_refit.main(argv)
            except Exception as e:
                buf.write("[!] grid correction failed: " + str(e))
            for line in buf.getvalue().splitlines():
                msg_q.put(('log', line))
            msg_q.put(('gridfix', (code, out)))

        threading.Thread(target=work, daemon=True).start()

    def on_set_crs():
        """Re-express the GCPs of an existing folder of crops in the output CRS.

        For sets extracted before 'Output GCP CRS' was set, so a project does not
        end up with some crops in one system and some in another. Pixels are not
        touched - only the GCP block is rewritten.
        """
        epsg = vars_['gcp_epsg'].get().strip()
        if not epsg.isdigit():
            messagebox.showwarning(
                "Set CRS",
                "Enter the EPSG code in 'Output GCP CRS' first, "
                "then choose the folder of crops.")
            return
        folder = filedialog.askdirectory(
            title="Folder of GeoTIFF crops to re-express in EPSG:" + epsg)
        if not folder:
            return
        setcrs_btn.configure(state='disabled')

        def work():
            try:
                retag_gcps(folder, int(epsg),
                           log=lambda m: msg_q.put(('log', str(m))))
            except Exception as e:
                msg_q.put(('log', "[!] could not set the CRS: " + str(e)))
            msg_q.put(('setcrs', None))

        threading.Thread(target=work, daemon=True).start()

    edit_btn.configure(command=on_edit_plots)
    tile_btn.configure(command=on_make_tiled)
    fixgrid_btn.configure(command=on_fix_grid)
    setcrs_btn.configure(command=on_set_crs)

    NL = chr(10)          # newline for message-box text, kept out of literals

    def poll():
        try:
            while True:
                kind, payload = msg_q.get_nowait()
                if kind == 'log':
                    logbox.configure(state='normal')
                    logbox.insert('end', payload + "\n")
                    logbox.see('end')
                    logbox.configure(state='disabled')
                elif kind == 'prog':
                    done, total = payload
                    prog['maximum'] = total
                    prog['value'] = done
                elif kind == 'done':
                    state['running'] = False
                    run_btn.configure(state='normal')
                    tile_btn.configure(state='normal')
                    cancel_btn.configure(state='disabled')
                    if payload:
                        open_btn.configure(state='normal')
                elif kind == 'setcrs':
                    setcrs_btn.configure(state='normal')
                elif kind == 'gridfix':
                    code, path = payload
                    fixgrid_btn.configure(state='normal')
                    if code == 0 and path and os.path.exists(path):
                        if messagebox.askyesno(
                                "Grid corrected",
                                "A corrected grid was written to:" + NL + NL
                                + path + NL + NL
                                + "Use it for the extraction?"):
                            vars_['plots'].set(path)
                            log(f"Plot boundaries set to the corrected grid: {path}")
                            log("With the grid corrected, leave 'Auto-fit to rows' OFF "
                                "— the polygon is now the plot.")
                    elif code == 3:
                        messagebox.showinfo(
                            "Grid not changed",
                            "The rows in this trial are drilled continuously across "
                            "plots, so a plot boundary has no signature in the "
                            "imagery and correcting the grid to the rows is not "
                            "meaningful." + NL + NL
                            + "The existing grid is the only definition of a plot "
                              "edge — use it as-is with auto-fit OFF.")
                    else:
                        messagebox.showwarning(
                            "Grid not corrected",
                            "The grid could not be corrected. See the log for detail.")
                elif kind == 'tiled_done':
                    state['running'] = False
                    run_btn.configure(state='normal')
                    tile_btn.configure(state='normal')
                    cancel_btn.configure(state='disabled')
                    if payload:
                        vars_['ortho_tif'].set(payload)
                        log(f"Ortho source set to the tiled copy: {payload}")
        except queue.Empty:
            pass
        root.after(120, poll)

    def on_close():
        save_config()
        cancel_ev.set()
        root.destroy()

    # ---- "Generate plot grid" tab -----------------------------------------
    # Builds the plot shapefile by clicking the 4 trial corners on an
    # orthomosaic; "Save & use" feeds it straight back into the form above.
    def _use_grid(path):
        vars_['plots'].set(path)
        log(f"Plot boundaries set to the new grid: {path}")
        nb.select(tab_extract)

    try:
        from grid_designer import build_grid_designer
        gd = build_grid_designer(tab_grid, on_saved=_use_grid, log=log,
                                 ortho_path=vars_['ortho_tif'].get().strip(),
                                 crs_hint=vars_['crs_override'].get().strip())
        gd.pack(fill='both', expand=True)
    except Exception as e:
        msg = ("The plot-grid designer could not be loaded:\n\n"
               f"{e}\n\nIt needs rasterio and pyshp, and grid_designer.py next "
               "to the application.")
        ttk.Label(tab_grid, text=msg, justify='left', padding=20,
                  foreground='#a00').pack(anchor='nw')

    # ---- "Training tiles" tab ---------------------------------------------
    # Cuts the per-plot crops into fixed-size tiles for labelling, and puts
    # labelled tiles back together into one image per plot so a whole plot can be
    # reviewed at once. Tile names follow the group's existing convention
    # (<stem>_tile_c<col>_r<row>), so tiles drop straight into the SAM3 notebook.
    try:
        from tiler import build_tiles_tab
        _out = vars_['outdir'].get().strip()
        tt = build_tiles_tab(tab_tiles, defaults=dict(
            src=os.path.join(_out, "RawCrops") if _out else "",
            out=os.path.join(_out, "Tiles") if _out else "",
            tile=int(cfg.get('tile_size', 1024)),
            overlap=int(cfg.get('tile_overlap', 256)),
            rows_only=bool(cfg.get('tile_rows_only', True)),
            pitch=cfg.get('row_pitch_mm', ''),
            plotw=cfg.get('plot_width_mm', ''),
        ))
        tt.pack(fill='both', expand=True)
    except Exception as e:
        ttk.Label(tab_tiles, justify='left', padding=20, foreground='#a00',
                  text=("The tiling tab could not be loaded:" + NL + NL + str(e)
                        + NL + NL + "It needs rasterio and Pillow, and tiler.py "
                        "next to the application.")).pack(anchor='nw')

    root.protocol("WM_DELETE_WINDOW", on_close)
    if vars_['plots'].get():
        root.after(400, refresh_id_fields)
    poll()
    root.mainloop()


def _selftest():
    """Headless import + license check for the frozen build (no GUI)."""
    ok = True
    for m in ("numpy", "PIL", "PIL.ImageCms", "PIL.ImageTk", "rasterio",
              "rasterio.features", "rasterio.warp", "shapefile", "shapely.geometry"):
        try:
            __import__(m); print(f"  [ok] {m}")
        except Exception as e:
            ok = False; print(f"  [FAIL] {m}: {e}")
    try:
        import Metashape
        lic = Metashape.License()
        print(f"  [ok] Metashape {Metashape.version}  licence_valid={lic.valid}")
    except Exception as e:
        print(f"  [warn] Metashape unavailable ({e}) - ortho-only mode still works")
    print("SELFTEST", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        _selftest()
    main()
