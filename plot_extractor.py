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
def detect_row_fit(im, rows_target=0, mode='width', buffer_frac=0.5):
    """Find the planted crop rows (vertical) in a rectified plot image and return
    a pixel crop box that tightly encloses them. Uses the ROW/RIDGE periodicity
    (detrended-brightness oscillation) reinforced by greenness, so it works even
    before emergence and doesn't chase stray edge grass.

    Because the caller renders WIDER than the grid boundary, this can INCLUDE rows
    the grid clipped, then crop to exactly `rows_target` rows (or all detected).

    mode: 'width' (band width, keep original centre) | 'width_shift' (band as
    detected, recentre on rows) | 'width_length' (also trim along-row ends).
    Returns dict(box=(x0,x1,y0,y1), rows, spacing_px, peaks) or None.
    """
    import numpy as np
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
                rows=len(grp), spacing_px=sp * ds,
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
    raw_mode = (p.get('raw_mode') or 'auto').lower()   # auto | bbox | mask | rectify
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
        if raw_mode not in ('rectify', 'auto'):
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
             'rectify': "Raw crop style: rectified to an upright plot rectangle"}
            .get(raw_mode, f"Raw crop style: {raw_mode}"))
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
    _cache_cap = max(2, workers + 1)

    def load_image(path):
        """Return (decoded RGB image, icc_profile_bytes). Thread-safe."""
        with _img_lock:
            hit = _img_lru.get(path)
            if hit is not None:
                _img_lru.move_to_end(path)
                return hit
        im = Image.open(path)
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

    def _render_fill(job):
        """Rectify mode with gap-fill. Rectify candidate frames to the SAME plot
        rectangle (they align pixel-for-pixel) and composite. Near-nadir frames
        go first and do the bulk; oblique frames only patch remaining holes.

        Each output pixel is native raw from EXACTLY ONE frame (one resample) -
        no blending, no duplicated plants. Before a fill frame is dropped in, its
        RGB is radiometrically matched to the already-placed pixels in the overlap
        zone (per-channel gain), so the frame boundary has no brightness/colour
        step (the 'two shades' seam). A per-pixel source-id map is kept so the
        seam lines can be QC-checked for a plant sitting across a boundary."""
        import numpy as np
        from PIL import ImageChops
        W, H = job['wpx'], job['hpx']
        layers = job['layers']
        good = [i for i, l in enumerate(layers)
                if l.get('aniso', 99.0) <= RENDER_MAX_ANISO]
        rest = [i for i, l in enumerate(layers)
                if l.get('aniso', 99.0) > RENDER_MAX_ANISO]
        order = good + rest              # good (clean) frames first, oblique last
        out = None                       # float32 HxWx3 accumulator
        filled = None                    # bool HxW: pixel already has data
        srcid = None                     # uint8 HxW: which layer filled each pixel
        used = 0
        for pos, i in enumerate(order):
            if pos >= len(good) and filled is not None and filled.all():
                break                    # holes closed -> stop before oblique frames
            layer = layers[i]
            im, icc = load_image(layer['path'])
            A, B, C, D = layer['corners']
            quad = [A[0], A[1], D[0], D[1], C[0], C[1], B[0], B[1]]
            rect = to_srgb_crop(im.transform((W, H), Image.QUAD, quad, Image.BICUBIC),
                                icc)
            vmask = np.asarray(_white_full(im.size).transform(
                (W, H), Image.QUAD, quad, Image.NEAREST)) > 127
            ra = np.asarray(rect, dtype=np.float32)
            if out is None:
                out = ra.copy()
                filled = vmask.copy()
                srcid = np.where(vmask, np.uint8(used), np.uint8(255))
                used += 1
                if filled.all():
                    break
                continue
            need = vmask & (~filled)                    # only still-empty pixels
            if not need.any():
                continue
            # radiometric match: make this frame match the already-placed pixels
            # in the overlap. A per-channel AFFINE match (mean AND contrast/std),
            # so a seam is removed even when the two frames differ in view angle
            # (oblique frames have deeper furrow shadows -> higher contrast). The
            # std ratio is clamped so noise isn't amplified / detail not crushed.
            ov = vmask & filled
            if ov.sum() > 500:
                for c in range(3):
                    oc = out[..., c][ov]; rc = ra[..., c][ov]
                    m_out, m_rect = float(oc.mean()), float(rc.mean())
                    s_out, s_rect = float(oc.std()), float(rc.std())
                    s_ratio = min(1.6, max(0.6, s_out / s_rect)) if s_rect > 1e-3 else 1.0
                    ra[..., c] = (ra[..., c] - m_rect) * s_ratio + m_out
                np.clip(ra, 0, 255, out=ra)
            out[need] = ra[need]
            srcid[need] = np.uint8(used)
            filled |= vmask
            used += 1
            if filled.all():
                break
        if out is None:
            out = np.zeros((H, W, 3), np.float32)
            filled = np.zeros((H, W), bool)
            srcid = np.full((H, W), 255, np.uint8)
        job['filled_pct'] = 100.0 * float(filled.mean())
        job['layers_used'] = used
        # seam location: median row where the source frame changes (the boundary
        # a plant could straddle). Recorded so a reviewer knows where to look;
        # not auto-flagged, since these long plots ALWAYS need 2 frames along
        # their length and the seam is a thin line (parallax tiny at emergence).
        job['seam_y_pct'] = ''
        if used > 1:
            vchg = np.zeros((H, W), bool)
            vchg[1:, :] = (srcid[1:, :] != srcid[:-1, :]) & filled[1:, :] & filled[:-1, :]
            rows = np.where(vchg.any(1))[0]
            if len(rows):
                job['seam_y_pct'] = f"{100.0 * float(np.median(rows)) / H:.0f}"
        return Image.fromarray(out.astype(np.uint8), 'RGB')

    def _apply_auto_fit(crop, job):
        """Detect the planted rows in the (wider-than-grid) rectified crop and
        crop back to exactly the rows. Records fitted width/length/area (m) and
        the row count, and shifts geotag GCPs into the new frame. On any doubt it
        keeps the full crop (never cuts rows), and REFUSES to fit a crop with
        significant black (incomplete coverage) since the area would be wrong."""
        import numpy as np
        # black fraction: incomplete-coverage crops give garbage row detection
        small = crop.resize((min(400, crop.width), min(1200, crop.height)))
        arr = np.asarray(small.convert('RGB'))
        black_pct = 100.0 * float(np.mean(arr.sum(2) < 24))
        if black_pct > 10.0:
            job['fit_note'] = f'low coverage ({black_pct:.0f}% black) - area NOT fitted'
            return crop                      # keep full crop, don't emit a bogus area
        try:
            fb = detect_row_fit(crop, rows_target=job.get('fit_rows', 0),
                                 mode=job.get('fit_mode', 'width'))
        except Exception as e:
            log(f"plot {job['pid']}: auto-fit detect failed ({e}); keeping full crop")
            fb = None
        if not fb:
            job['fit_note'] = 'rows not detected'
            return crop
        # sanity: fitted width must be within 60% of the expected 7-row span
        exp_w = job.get('mpp_x', 0) * (fb['box'][1] - fb['box'][0])
        if job.get('fit_rows') and fb['rows'] < job['fit_rows']:
            job['fit_note'] = f"only {fb['rows']}/{job['fit_rows']} rows found - kept full crop"
            return crop
        x0, x1, y0, y1 = fb['box']
        W, H = crop.size
        x0 = max(0, min(x0, W - 2)); x1 = max(x0 + 1, min(x1, W))
        y0 = max(0, min(y0, H - 2)); y1 = max(y0 + 1, min(y1, H))
        crop = crop.crop((x0, y0, x1, y1))
        if job['geotag'] and job.get('gcps'):
            ng = []
            for (r, c, lon, lat, z) in job['gcps']:
                nc, nr = c - x0, r - y0
                if -0.5 <= nc <= (x1 - x0) + 0.5 and -0.5 <= nr <= (y1 - y0) + 0.5:
                    ng.append((nr, nc, lon, lat, z))
            job['gcps'] = ng
        mppx, mppy = job.get('mpp_x', 0.0), job.get('mpp_y', 0.0)
        job['fit_rows_found'] = fb['rows']
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
            job['flag'] = "" if fp >= 99.5 else \
                f"  [!] {fp:.0f}% covered even after fill"
            nfr = job.get('layers_used', 1)
            job['cam_label'] = (f"{nfr} frame{'s' if nfr != 1 else ''}"
                                + (" (gap-filled, matched)" if nfr > 1 else " (single)"))
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
                A, B, C, D = ring0
                uAB = _uv_wid(A, B); uDC = _uv_wid(D, C)      # width-edge directions
                pm = fit_probe_m
                ring0 = [[A[0] - uAB[0] * pm, A[1] - uAB[1] * pm],
                         [B[0] + uAB[0] * pm, B[1] + uAB[1] * pm],
                         [C[0] + uDC[0] * pm, C[1] + uDC[1] * pm],
                         [D[0] - uDC[0] * pm, D[1] - uDC[1] * pm]]
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
                prefer_nadir = raw_mode in ('rectify', 'auto')
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
                    uvs = [cam.project(q) for q in pts3d]
                    uvs_p = [cam.project(q) for q in pts3d_plot]
                    if any(uv is None for uv in uvs) or any(uv is None for uv in uvs_p):
                        continue
                    w, h = cam.sensor.width, cam.sensor.height
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
                    if prefer_nadir:
                        # bucket foreshortening coarsely so that among comparably
                        # near-nadir frames, the one that COVERS more plot wins.
                        # This avoids picking a barely-clipping nadir frame (which
                        # leaves black) over a slightly-more-oblique full-coverage one.
                        key = (round(aniso / 0.35) * 0.35, -inside,
                               round(centre_dist, 3), over)
                    else:
                        key = (-inside, -int(crop_clear), round(centre_dist, 4), over)
                    cand.append((key, cam, uvs, uvs_p, w, h, inside, aniso, gsd))
                cand.sort(key=lambda c: c[0])
                if not cand:
                    row['note'] = 'no covering raw image'
                    log(f"plot {pid}: NO covering raw image")
                else:
                    # PHASE A (this main thread): resolve all Metashape geometry and
                    # queue a render job. The heavy pixel work runs later in parallel.
                    for rank, (_s, cam, uvs, uvs_p, w, h, inside, aniso, gsd) in \
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
                        # Gap-fill: for a rectified rank-1 crop, gather the top
                        # near-nadir candidate frames (cand is already sorted
                        # nadir-first) so a neighbour can fill what the best single
                        # frame misses. Each layer carries its foreshortening so the
                        # renderer can skip oblique frames that would smear.
                        layers = None
                        if rank == 1 and rectified:
                            layers = []
                            for (_s2, cam2, _u2, uvs_p2, _w2, _h2, _in2, an2, _g2) in cand[:8]:
                                p2 = resolve_photo(cam2)
                                if p2 is None:
                                    continue
                                layers.append(dict(
                                    path=p2, aniso=an2,
                                    corners=[(uv.x, uv.y) for uv in uvs_p2]))
                            if len(layers) <= 1:
                                layers = None          # nothing to fill from
                        # ground metres-per-pixel of the rectified output (for
                        # auto-fit area + width in real units)
                        mpp_x = (gW / wpx) if (rectified and gW and wpx) else 0.0
                        mpp_y = (gL / hpx) if (rectified and gL and hpx) else 0.0
                        job = dict(pid=pid, path=path, x0=x0, y0=y0, x1=x1, y1=y1,
                                   eff=eff, plot_uv=plot_uv, wpx=wpx, hpx=hpx,
                                   out=out, geotag=geotag, gcps=gcps, layers=layers,
                                   is_rank1=(rank == 1), row=row, cam_label=cam.label,
                                   cw=cw2, ch=ch2, full=False, cov=0.0, flag='',
                                   auto_fit=(auto_fit and rectified and rank == 1),
                                   fit_rows=fit_rows, fit_mode=fit_mode,
                                   mpp_x=mpp_x, mpp_y=mpp_y)
                        if rank == 1:
                            # coverage: fraction of the plot polygon imaged in this frame
                            xs_r = [q[0] for q in ring0]; ys_r = [q[1] for q in ring0]
                            gx0, gx1 = min(xs_r), max(xs_r)
                            gy0, gy1 = min(ys_r), max(ys_r)
                            NN = 15; tot = infr = 0
                            for gi in range(NN):
                                for gj in range(NN):
                                    lon = gx0 + (gx1 - gx0) * gi / (NN - 1)
                                    lat = gy0 + (gy1 - gy0) * gj / (NN - 1)
                                    if not _point_in_ring(lon, lat, ring0):
                                        continue
                                    tot += 1
                                    p3c, _z = to3d(lon, lat)
                                    uvc = cam.project(p3c)
                                    if uvc is not None and 0 <= uvc.x < w and 0 <= uvc.y < h:
                                        infr += 1
                            cov = (100.0 * infr / tot) if tot else 0.0
                            job['full'] = (inside == len(uvs))
                            job['cov'] = cov
                            job['flag'] = "" if inside == len(uvs) else \
                                f"  [!] only {inside}/{len(uvs)} corners in-frame"
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


def save_geotiff(path, pil_img, gcps, tags, log):
    """Write an RGB crop as a lossless GeoTIFF with WGS84 GCPs and source EXIF
    copied into the GDAL metadata (readable via gdalinfo / rasterio.tags())."""
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
            dst.gcps = ([GroundControlPoint(row=r, col=c, x=x, y=y, z=z)
                         for r, c, x, y, z in gcps], CRS.from_epsg(4326))
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
            bw, bh = max(1, int(round(sxR - sxL))), max(1, int(round(syB - syT)))
            if bm['key'] != (bw, bh):
                bm['photo'] = ImageTk.PhotoImage(bm['img'].resize((bw, bh), _I.BILINEAR))
                bm['key'] = (bw, bh)
            cv.create_image(sxL, syT, anchor='nw', image=bm['photo'])
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
            title="Ortho basemap (georeferenced GeoTIFF, same CRS as plots)",
            filetypes=[("GeoTIFF", "*.tif *.tiff")])
        if not pth:
            return
        try:
            import numpy as np
            import rasterio
            from rasterio.enums import Resampling
            from PIL import Image as _I
            with rasterio.open(pth) as ds:
                dec = max(1, int(max(ds.width, ds.height) / 2500))
                ow, oh = max(1, ds.width // dec), max(1, ds.height // dec)
                bands = [1, 2, 3] if ds.count >= 3 else [1, 1, 1]
                arr = ds.read(bands, out_shape=(3, oh, ow), resampling=Resampling.average)
                b = ds.bounds
            bm['img'] = _I.fromarray(np.transpose(arr, (1, 2, 0)).astype('uint8'), 'RGB')
            bm['bounds'] = (b.left, b.bottom, b.right, b.top)
            bm['key'] = None
            log(f"Basemap loaded: {os.path.basename(pth)}")
            redraw()
        except Exception as e:
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

    frm = ttk.Frame(root, padding=10)
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
        'format': tk.StringVar(value=cfg.get('format', 'TIFF (lossless)')),
        'to_srgb': tk.BooleanVar(value=cfg.get('to_srgb', False)),
        'workers': tk.StringVar(value=str(cfg.get('workers', 0))),
        'auto_fit': tk.BooleanVar(value=cfg.get('auto_fit', False)),
        'rows_per_plot': tk.StringVar(value=str(cfg.get('rows_per_plot', 0))),
        'fit_mode': tk.StringVar(value=cfg.get('fit_mode', 'width')),
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
    raw_mode_combo = ttk.Combobox(checks, textvariable=vars_['raw_mode'], width=26,
                                  state='readonly',
                                  values=["Native (no resample, no black)",
                                          "Auto-straighten tilted", "Bounding box",
                                          "Masked to plot", "Rectified rectangle"])
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

    def toggle_source(*_):
        ortho_only = vars_['source'].get().startswith('Orthomosaic')
        psx_ent.configure(state='disabled' if ortho_only else 'normal')
        ot_ent.configure(state='normal' if ortho_only else 'disabled')
        raw_chk.configure(state='disabled' if ortho_only else 'normal')
        raw_mode_combo.configure(state='disabled' if ortho_only else 'readonly')
        srgb_chk.configure(state='disabled' if ortho_only else 'normal')
        if ortho_only:
            vars_['do_ortho'].set(True)
    vars_['source'].trace_add('write', toggle_source)
    toggle_source()

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
                                        "Rectified rectangle": "rectify"}.get(
                                            vars_['raw_mode'].get(), "native"),
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
    edit_btn.configure(command=on_edit_plots)
    tile_btn.configure(command=on_make_tiled)

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
