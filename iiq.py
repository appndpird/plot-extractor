"""iiq.py - Phase One IIQ support for PlotExtractor.

  * SDK discovery and the p1tool developer (iiq_tools\\p1tool.exe, built from
    iiq_tools\\src\\p1tool.cpp against the Phase One ImageSDK).
  * Developer: develops whole frames or plot boxes straight from IIQ files. Used by
    run_extraction (raw_source='iiq') and by the QA re-development check.
  * convert_folder + build_convert_tab: the "IIQ conversion" tab (IIQ -> TIFF, GeoTIFF
    or JPEG, 8 or 16 bit).
  * build_iiq_extract_tab: the "IIQ plot extraction" tab - the plot-wise extraction of
    the Extract tab, with the pixels developed from the IIQ files instead of read from
    the JPEG/TIFF renditions.

IIQ files are only ever READ. Nothing is written next to them: settings, job lists and
temporary pixels go to the chosen output folder, and an output folder inside the IIQ
folder is refused. Each run records the size and modification time of every IIQ it
used and checks them again at the end.

Acknowledgements
  The approach of developing each plot's box directly from the IIQ (back-projected plot
  corners -> crop -> Phase One SDK convert) and the HR settings preset come from the
  APPN-UQ PhaseOne_Pipeline (Australian Plant Phenomics Network, The University of
  Queensland; MIT licence, Copyright (c) 2025 APPN-UQ; contact appn@uq.edu.au).
  The Phase One ImageSDK is Copyright Phase One A/S and is used under Phase One's SDK
  licence. See THIRD_PARTY_NOTICES.md. A stand-alone IIQ converter with a GUI is
  available from the DPIRD APPN node: github.com/appndpird/P1ImageSDK-IIQ-Converter.
"""
from __future__ import annotations

import csv
import json
import math
import os
import shutil
import subprocess
import sys
import threading
import time
import uuid

CREDIT = ("IIQ processing approach adapted from the APPN-UQ PhaseOne_Pipeline "
          "(Australian Plant Phenomics Network, The University of Queensland; MIT licence). "
          "Pixels developed by the Phase One ImageSDK (c) Phase One A/S, used under "
          "Phase One's SDK licence.")

IIQ_EXTS = (".iiq",)

# label -> preset file in iiq_tools\presets
PRESETS = {
    "APPN HR_clip - noise reduction 0.8 + sharpening 3.1 (recommended)": "appn_hr.json",
    "Neutral - all processing off (least processed, visibly noisy)": "neutral.json",
    "APPN template - exposure +0.8, highlight/shadow recovery": "appn_template.json",
}
PRESET_DEFAULT = "APPN HR_clip - noise reduction 0.8 + sharpening 3.1 (recommended)"
CUSTOM = "Custom settings file..."


def app_dir():
    return getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))


def tools_dir():
    return os.path.join(app_dir(), "iiq_tools")


def p1tool_path():
    p = os.path.join(tools_dir(), "p1tool.exe")
    return p if os.path.exists(p) else None


def _preset_or_default(label):
    return label if label in PRESETS or label == CUSTOM else PRESET_DEFAULT


def preset_path(label, custom=""):
    if label == CUSTOM:
        return custom
    return os.path.join(tools_dir(), "presets", PRESETS.get(label, PRESETS[PRESET_DEFAULT]))


def _dll_version(path):
    """File version of a DLL (e.g. '4.1.84.0'), or '' if unavailable."""
    try:
        import ctypes
        from ctypes import wintypes
        ver = ctypes.WinDLL("version")
        size = ver.GetFileVersionInfoSizeW(path, None)
        if not size:
            return ""
        buf = ctypes.create_string_buffer(size)
        ver.GetFileVersionInfoW(path, 0, size, buf)
        p = ctypes.c_void_p()
        n = wintypes.UINT()
        ver.VerQueryValueW(buf, "\\", ctypes.byref(p), ctypes.byref(n))

        class FFI(ctypes.Structure):
            _fields_ = [(k, wintypes.DWORD) for k in (
                "sig", "struc", "ms", "ls", "pms", "pls", "mask", "flags", "os", "type",
                "sub", "dms", "dls")]
        f = ctypes.cast(p, ctypes.POINTER(FFI)).contents
        return f"{f.ms >> 16}.{f.ms & 0xFFFF}.{f.ls >> 16}.{f.ls & 0xFFFF}"
    except Exception:
        return ""


def find_sdk(hint=""):
    """Locate a Phase One ImageSDK runtime: a folder with ImageSdkCpp.dll (in it or in
    its lib\\) and SensorProfiles\\ (in it or its parent). Tries the user's folder, then
    the copy shipped with the application (iiq_tools\\sdk). Returns a dict or None."""
    cands = [hint] if hint else []
    cands += [os.path.join(tools_dir(), "sdk"), tools_dir()]
    for d in cands:
        if not d or not os.path.isdir(d):
            continue
        for dll_dir in (d, os.path.join(d, "lib")):
            dll = os.path.join(dll_dir, "ImageSdkCpp.dll")
            if not os.path.exists(dll):
                continue
            for prof in (os.path.join(d, "SensorProfiles"),
                         os.path.join(os.path.dirname(d), "SensorProfiles"),
                         os.path.join(dll_dir, "SensorProfiles")):
                if os.path.isdir(prof):
                    return dict(root=d, dll_dir=dll_dir, profiles=prof,
                                version=_dll_version(dll))
    return None


def sdk_status(hint=""):
    exe = p1tool_path()
    sdk = find_sdk(hint)
    if not exe:
        return False, "p1tool.exe not found in iiq_tools\\ (see iiq_tools\\build_p1tool.bat)"
    if not sdk:
        return False, ("Phase One ImageSDK not found - set the SDK folder "
                       "(the folder with ImageSdkCpp.dll and SensorProfiles)")
    return True, f"Phase One ImageSDK {sdk['version'] or '(version unknown)'}  -  {sdk['root']}"


def _norm(p):
    return os.path.normcase(os.path.realpath(os.path.abspath(p)))


def assert_safe_output(src_dirs, out_dir):
    """Refuse an output folder that is (inside) an IIQ source folder."""
    o = _norm(out_dir)
    for s in src_dirs:
        if not s:
            continue
        sn = _norm(s)
        if o == sn or o.startswith(sn.rstrip("\\/") + os.sep):
            raise RuntimeError(
                "The output folder is inside the IIQ folder. The IIQ folder is never "
                "written to - choose an output folder somewhere else.")


def list_iiq(folder, recursive=True):
    out = []
    if recursive:
        for root, _d, files in os.walk(folder):
            out += [os.path.join(root, f) for f in files if f.lower().endswith(IIQ_EXTS)]
    else:
        out = [os.path.join(folder, f) for f in os.listdir(folder)
               if f.lower().endswith(IIQ_EXTS)]
    return sorted(out)


def snapshot(paths):
    snap = {}
    for p in paths:
        try:
            st = os.stat(p)
            snap[p] = (st.st_size, st.st_mtime_ns)
        except OSError:
            snap[p] = None
    return snap


def snapshot_changed(snap):
    """Paths whose size or modification time changed since the snapshot."""
    now = snapshot(list(snap))
    return [p for p in snap if snap[p] != now[p]]


class Developer:
    """Develops IIQ frames or boxes with p1tool. Thread-safe; at most `parallel`
    p1tool processes at once (they share the GPU)."""

    def __init__(self, settings, work_dir, sdk_hint="", log=print, parallel=2):
        self.exe = p1tool_path()
        self.sdk = find_sdk(sdk_hint)
        ok, msg = sdk_status(sdk_hint)
        if not ok:
            raise RuntimeError(msg)
        if not settings or not os.path.exists(settings):
            raise RuntimeError(f"IIQ settings file not found: {settings}")
        self.settings = settings
        self.work = work_dir
        os.makedirs(work_dir, exist_ok=True)
        self.log = log
        self.sem = threading.Semaphore(max(1, int(parallel)))
        self._snap = {}
        self._snap_lock = threading.Lock()
        self._info = {}

    def _env(self):
        env = dict(os.environ)
        env["PATH"] = self.sdk["dll_dir"] + os.pathsep + env.get("PATH", "")
        return env

    def _remember(self, path):
        with self._snap_lock:
            if path not in self._snap:
                self._snap.update(snapshot([path]))

    def used_unchanged(self):
        """(number of IIQ files used, list of those that changed during the run)."""
        with self._snap_lock:
            snap = dict(self._snap)
        return len(snap), snapshot_changed(snap)

    def run(self, jobs, on_line=None, cancel=None):
        """jobs: list of dicts (mode, iiq, out, x, y, w, h, scale, bits, quality, lzw).
        Returns one result dict per job: ok, w, h, bits, error."""
        res = [dict(ok=False, error="not run") for _ in jobs]
        if not jobs:
            return res
        for j in jobs:
            self._remember(j["iiq"])
        tsv = os.path.join(self.work, f"jobs_{uuid.uuid4().hex}.tsv")
        with open(tsv, "w", encoding="utf-8", newline="\n") as f:
            for j in jobs:
                f.write("\t".join(str(v) for v in (
                    j.get("mode", "raw"), j["iiq"], j.get("out", ""), self.settings,
                    j.get("x", 0), j.get("y", 0), j.get("w", 0), j.get("h", 0),
                    j.get("scale", 1.0), j.get("bits", 8), j.get("quality", 95),
                    1 if j.get("lzw", True) else 0)) + "\n")
        flags = getattr(subprocess, "CREATE_NO_WINDOW", 0)
        with self.sem:
            proc = subprocess.Popen([self.exe, self.sdk["profiles"], tsv], cwd=self.work,
                                    env=self._env(), stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, text=True,
                                    encoding="utf-8", errors="replace",
                                    creationflags=flags)
            try:
                for line in proc.stdout:
                    if cancel is not None and cancel.is_set():
                        proc.kill()
                        break
                    if not line.startswith("@@"):
                        continue
                    parts = line[2:].strip().split(" ", 2)
                    tag = parts[0]
                    try:
                        i = int(parts[1])
                    except (IndexError, ValueError):
                        continue
                    if not 0 <= i < len(jobs):
                        continue
                    if tag == "OK":
                        w, h, b = (int(v) for v in parts[2].split())
                        res[i] = dict(ok=True, w=w, h=h, bits=b, error="")
                    elif tag == "INFO":
                        w, h, rw, rh = (int(v) for v in parts[2].split())
                        res[i] = dict(ok=True, w=w, h=h, raw_w=rw, raw_h=rh, error="")
                    elif tag == "ERR":
                        res[i] = dict(ok=False, error=parts[2] if len(parts) > 2 else "error")
                    if on_line:
                        on_line(i, res[i])
            finally:
                proc.wait()
        try:
            os.remove(tsv)
        except OSError:
            pass
        return res

    def info(self, path):
        """(width, height) of the developed full-size frame."""
        if path not in self._info:
            r = self.run([dict(mode="info", iiq=path)])[0]
            if not r["ok"]:
                raise RuntimeError(f"cannot read {os.path.basename(path)}: {r['error']}")
            self._info[path] = (r["w"], r["h"])
        return self._info[path]

    def develop(self, path, box=None, bits=8):
        """Develop a whole frame (box None) or the box (x, y, w, h) in full-size pixels.
        Returns an HxWx3 numpy array, uint8 (bits=8) or uint16 (bits=16)."""
        import numpy as np
        x, y, w, h = box if box else (0, 0, 0, 0)
        bits = 16 if int(bits) == 16 else 8
        dt = np.uint16 if bits == 16 else np.uint8
        tmp = os.path.join(self.work, f"px_{uuid.uuid4().hex}.rgb")
        try:
            r = self.run([dict(mode="raw", iiq=path, out=tmp, x=x, y=y, w=w, h=h,
                               bits=bits)])[0]
            if not r["ok"]:
                raise RuntimeError(f"IIQ develop failed ({os.path.basename(path)}): "
                                   f"{r['error']}")
            a = np.fromfile(tmp, dt)
            if a.size != r["w"] * r["h"] * 3:
                raise RuntimeError(f"IIQ develop returned {a.size} samples for "
                                   f"{r['w']}x{r['h']}")
            return a.reshape(r["h"], r["w"], 3)
        finally:
            try:
                os.remove(tmp)
            except OSError:
                pass

    def cleanup(self):
        try:
            if os.path.isdir(self.work) and not os.listdir(self.work):
                os.rmdir(self.work)
        except OSError:
            pass


class ScaledCam:
    """A Metashape camera whose pixel coordinates are scaled by k - for a project
    aligned on downscaled renditions (e.g. 0.3x TIFFs) while the pixels are developed
    from the full-size IIQ."""

    class _Sensor:
        def __init__(self, s, k):
            self._s, self.width, self.height = s, int(round(s.width * k)), \
                int(round(s.height * k))

        def __getattr__(self, n):
            return getattr(self._s, n)

    def __init__(self, cam, k):
        self._c, self._k = cam, k
        self.sensor = ScaledCam._Sensor(cam.sensor, k)

    def project(self, p):
        uv = self._c.project(p)
        if uv is None:
            return None
        import Metashape
        return Metashape.Vector([uv.x * self._k, uv.y * self._k])

    def __getattr__(self, n):
        return getattr(self._c, n)


# ---------------------------------------------------------------------------
# frame georeferencing (GeoTIFF conversion)
# ---------------------------------------------------------------------------
def open_chunk(psx, chunk=""):
    import Metashape as ms
    doc = ms.Document()
    doc.open(psx, read_only=True)
    chunks = [c for c in doc.chunks if c is not None]
    ch = None
    for c in chunks:
        if chunk and (c.label or "").strip() == chunk:
            ch = c
    ch = ch or doc.chunk or (chunks[0] if chunks else None)
    if ch is None or ch.crs is None:
        raise RuntimeError("The Metashape project has no georeferenced chunk.")
    return doc, ch


def _epsg_of(crs):
    try:
        a = crs.authority or ""
        if a.upper().startswith("EPSG::"):
            return int(a.split("::")[1])
    except Exception:
        pass
    return None


def _is_geographic(crs):
    """True for a lon/lat chunk CRS (WGS 84, GDA2020 geographic, ...)."""
    try:
        w = crs.wkt or ""
        return w.lstrip().upper().startswith(("GEOGCS", "GEOGCRS", "GEODCRS"))
    except Exception:
        return False


def frame_gcps(ch, cam, out_w, out_h, n=15):
    """Ground control points for a whole developed frame: a ground grid around the
    camera (heights from the DEM) projected with the full camera model (lens
    distortion included) and kept where it lands inside the frame. Returns
    [(row, col, X, Y, Z)] in the chunk CRS."""
    import Metashape as ms
    T = ch.transform.matrix
    Tinv = T.inv()
    dem = ch.elevation
    sw, sh = cam.sensor.width, cam.sensor.height
    kx, ky = out_w / float(sw), out_h / float(sh)
    c = ch.crs.project(T.mulp(cam.center))

    def zat(x, y):
        try:
            z = dem.altitude(ms.Vector([x, y])) if dem is not None else None
        except Exception:
            z = None
        return z if z is not None and z == z else None

    zg = zat(c.x, c.y)
    if zg is None:
        zg = ch.crs.project(T.mulp(ch.region.center)).z
    # ground footprint from the focal length; a point far outside the view can be
    # folded back into the frame by the lens-distortion polynomial, so every point
    # must also lie inside the (undistorted) field of view, checked in camera space
    f = float(cam.sensor.calibration.f or 0) or float(sw)
    hx, hy = 0.5 * sw / f, 0.5 * sh / f               # tan of the half view angles
    H = abs(c.z - zg)
    reach = max(2.0, 1.15 * H * max(hx, hy))           # metres
    # the grid is laid out in metres; a geographic chunk CRS (degrees) gets the
    # local metres-per-degree so the ground grid is the same size on the ground
    if _is_geographic(ch.crs):
        mx, my = 111320.0 * math.cos(math.radians(c.y)), 110540.0
    else:
        mx = my = 1.0
    Cinv = cam.transform.inv()
    pts = []
    for i in range(n):
        for j in range(n):
            x = c.x + (-reach + 2 * reach * i / (n - 1)) / mx
            y = c.y + (-reach + 2 * reach * j / (n - 1)) / my
            z = zat(x, y)
            if z is None:
                z = zg
            p3 = Tinv.mulp(ch.crs.unproject(ms.Vector([x, y, z])))
            pc = Cinv.mulp(p3)
            if pc.z <= 0 or abs(pc.x / pc.z) > 1.05 * hx or abs(pc.y / pc.z) > 1.05 * hy:
                continue
            uv = cam.project(p3)
            if uv is None:
                continue
            col, row = uv.x * kx, uv.y * ky
            if 0 <= col <= out_w and 0 <= row <= out_h:
                pts.append((row, col, x, y, z))
    return pts


def projected_epsg_for(lon, lat, geo_epsg=4326):
    """The UTM / MGA zone EPSG for a lon/lat: MGA2020 for GDA2020 (7844), MGA94 for
    GDA94 (4283), WGS 84 UTM otherwise."""
    zone = int((lon + 180.0) // 6) + 1
    if geo_epsg == 7844 and 46 <= zone <= 59:
        return 7800 + zone
    if geo_epsg == 4283 and 46 <= zone <= 59:
        return 28300 + zone
    return (32700 if lat < 0 else 32600) + zone


def gcps_to_epsg(gcps, src_epsg, dst_epsg):
    """[(row, col, X, Y, Z)] reprojected from src_epsg to dst_epsg (Z kept)."""
    from rasterio.crs import CRS
    from rasterio.warp import transform
    xs, ys = transform(CRS.from_epsg(src_epsg), CRS.from_epsg(dst_epsg),
                       [g[2] for g in gcps], [g[3] for g in gcps])
    return [(g[0], g[1], x, y, g[4]) for g, x, y in zip(gcps, xs, ys)]


def georeference_tiff(path, gcps, epsg, log=print):
    """Add georeferencing to a TIFF written by p1tool (pixels and EXIF untouched): an
    affine fit through the GCPs when it holds to 10 cm, else the GCPs themselves; the
    GCPs are always kept in the PE_GCPS tag. Returns the max affine residual (m)."""
    import rasterio
    from rasterio.control import GroundControlPoint
    from rasterio.crs import CRS
    import plot_extractor as pe
    pts = [GroundControlPoint(row=r, col=c, x=x, y=y, z=z) for r, c, x, y, z in gcps]
    if len(pts) < 3:
        log(f"[!] {os.path.basename(path)}: too few ground points in the frame - "
            "left without georeference")
        return None
    crs = CRS.from_epsg(epsg)
    aff, res = pe._affine_from_gcps(pts)
    with rasterio.open(path, "r+") as d:
        d.update_tags(PE_GCPS=json.dumps([[g.row, g.col, g.x, g.y, g.z] for g in pts]),
                      PE_GCPS_EPSG=str(epsg))
        if aff is not None and res.max() <= 0.10:
            d.transform = aff
            d.crs = crs
            d.update_tags(PE_AFFINE_RESID=f"mean={res.mean() * 1000:.1f}mm "
                                          f"max={res.max() * 1000:.1f}mm n={len(pts)}")
        else:
            d.gcps = (pts, crs)
    return float(res.max()) if res is not None else None


# ---------------------------------------------------------------------------
# IIQ conversion (tab 1)
# ---------------------------------------------------------------------------
def convert_folder(p, log, progress, cancel):
    """p: iiq_dir, recursive, out_dir, fmt ('tiff'|'geotiff'|'jpg'), bits, settings,
    quality, lzw, scale, psx, chunk, sdk."""
    src, out = p["iiq_dir"], p["out_dir"]
    if not os.path.isdir(src):
        raise RuntimeError("IIQ folder not found.")
    if not out:
        raise RuntimeError("Choose an output folder.")
    assert_safe_output([src], out)
    files = list_iiq(src, p.get("recursive", True))
    if not files:
        raise RuntimeError("No .IIQ files in the folder.")
    fmt = p.get("fmt", "tiff")
    bits = 8 if fmt == "jpg" else int(p.get("bits", 8))
    scale = float(p.get("scale") or 1.0)
    os.makedirs(out, exist_ok=True)
    log(CREDIT)
    settings = p["settings"]
    shutil.copy2(settings, os.path.join(out, "IIQ_settings_used.json"))
    dev = Developer(settings, os.path.join(out, "_iiq_work"), p.get("sdk", ""), log)
    log(sdk_status(p.get("sdk", ""))[1])
    log(f"{len(files)} IIQ file(s); output {fmt.upper()} {bits}-bit, scale {scale:g}, "
        f"settings: {os.path.basename(settings)}")
    snap = snapshot(files)

    ch = cams = epsg = geo_epsg = None
    if fmt == "geotiff":
        if not p.get("psx"):
            raise RuntimeError("GeoTIFF needs the Metashape project the frames were "
                               "aligned in (each frame is georeferenced from its aligned "
                               "camera and the DEM). Choose TIFF to keep only the "
                               "camera GPS in the EXIF.")
        _doc, ch = open_chunk(p["psx"], p.get("chunk", ""))
        # ground points come in the chunk CRS; the frames are written in the chosen
        # output CRS (default EPSG:7850), else a geographic chunk goes to its local
        # UTM / MGA zone and a projected one stays as it is
        epsg = _epsg_of(ch.crs) or (4326 if _is_geographic(ch.crs) else None)
        want = str(p.get("epsg") or "").strip()
        if epsg is None:
            log("[!] chunk CRS has no EPSG code - frames written with GCPs in the chunk "
                "coordinates, labelled EPSG:4326")
        elif want.isdigit() and int(want) != epsg:
            geo_epsg, epsg = epsg, int(want)
        elif not want and _is_geographic(ch.crs):
            cc = ch.crs.project(ch.transform.matrix.mulp(ch.region.center))
            geo_epsg = epsg
            epsg = projected_epsg_for(cc.x, cc.y, geo_epsg)
        log(f"GeoTIFF CRS: EPSG:{epsg or 4326} (chunk CRS {ch.crs.name})")
        cams = {}
        for c in ch.cameras:
            if c.transform is None:
                continue
            cams[c.label.lower()] = c
            try:
                cams[os.path.splitext(os.path.basename(c.photo.path))[0].lower()] = c
            except Exception:
                pass
        log(f"GeoTIFF: {len(cams)} aligned cameras in '{ch.label}', CRS {ch.crs.name}")

    ext = ".jpg" if fmt == "jpg" else ".tif"
    jobs = []
    for f in files:
        rel = os.path.relpath(f, src)
        o = os.path.join(out, os.path.splitext(rel)[0] + ext)
        os.makedirs(os.path.dirname(o), exist_ok=True)
        jobs.append(dict(mode="jpg" if fmt == "jpg" else "tiff", iiq=f, out=o, scale=scale,
                         bits=bits, quality=int(p.get("quality", 95)),
                         lzw=bool(p.get("lzw", True))))
    rows = []
    done = [0]
    t0 = time.time()
    batch = 10
    for b0 in range(0, len(jobs), batch):
        if cancel.is_set():
            log("Cancelled.")
            break
        part = jobs[b0:b0 + batch]
        res = dev.run(part, cancel=cancel)
        for j, r in zip(part, res):
            rec = dict(iiq=j["iiq"], output=j["out"], ok=r["ok"], width=r.get("w", ""),
                       height=r.get("h", ""), bits=r.get("bits", ""), error=r["error"],
                       georef="")
            if r["ok"] and fmt == "geotiff":
                stem = os.path.splitext(os.path.basename(j["iiq"]))[0].lower()
                cam = cams.get(stem)
                if cam is None:
                    rec["georef"] = "no aligned camera"
                    log(f"[!] {stem}: not an aligned camera in the project - "
                        "written without georeference")
                else:
                    gc = frame_gcps(ch, cam, r["w"], r["h"])
                    if geo_epsg and gc:
                        gc = gcps_to_epsg(gc, geo_epsg, epsg)
                    mx = georeference_tiff(j["out"], gc, epsg or 4326, log)
                    rec["georef"] = ("affine" if mx is not None and mx <= 0.10 else "gcps") \
                        + (f" max {mx * 100:.1f} cm" if mx is not None else "")
            if not r["ok"]:
                log(f"[!] {os.path.basename(j['iiq'])}: {r['error']}")
            rows.append(rec)
            done[0] += 1
        progress(done[0], len(jobs))
        el = time.time() - t0
        log(f"{done[0]}/{len(jobs)} converted  ({el / max(done[0], 1):.1f} s per frame)")
    if rows:
        with open(os.path.join(out, "conversion_log.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    changed = snapshot_changed(snap)
    log(f"Source check: {len(snap) - len(changed)}/{len(snap)} IIQ files unchanged"
        + (f"  [!] CHANGED: {changed[:5]}" if changed else ""))
    dev.cleanup()
    ok = sum(1 for r in rows if r["ok"])
    log(f"DONE: {ok}/{len(files)} converted -> {out}")
    return ok


# ---------------------------------------------------------------------------
# UI helpers
# ---------------------------------------------------------------------------
_CFG = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0] if getattr(
    sys, "frozen", False) else __file__)), "iiq_config.json")


def _cfg_load():
    try:
        with open(_CFG) as f:
            return json.load(f)
    except Exception:
        return {}


def _cfg_save(d):
    try:
        cur = _cfg_load()
        cur.update(d)
        with open(_CFG, "w") as f:
            json.dump(cur, f, indent=1)
    except Exception:
        pass


def _runner(frm, V, key, log_box, pbar, run_btn, cancel_btn, work):
    """Run work(log, progress, cancel) in a thread; stream its log to log_box."""
    import queue
    q = queue.Queue()
    st = {"busy": False}
    cancel_ev = threading.Event()

    def log(m):
        q.put(("log", time.strftime("%H:%M:%S ") + str(m)))

    def progress(a, b):
        q.put(("prog", (a, b)))

    def poll():
        try:
            while True:
                k, v = q.get_nowait()
                if k == "log":
                    log_box.insert("end", v + "\n")
                    log_box.see("end")
                elif k == "prog":
                    pbar["maximum"] = max(1, v[1])
                    pbar["value"] = v[0]
                elif k == "end":
                    st["busy"] = False
                    run_btn.state(["!disabled"])
                    cancel_btn.state(["disabled"])
        except queue.Empty:
            pass
        frm.after(150, poll)

    def start():
        if st["busy"]:
            return
        _cfg_save({f"{key}_{k}": (v.get()) for k, v in V.items()})
        st["busy"] = True
        cancel_ev.clear()
        run_btn.state(["disabled"])
        cancel_btn.state(["!disabled"])
        log_box.delete("1.0", "end")

        def th():
            try:
                work(log, progress, cancel_ev)
            except Exception as e:
                import traceback
                log(f"ERROR: {e}")
                log(traceback.format_exc())
            q.put(("end", None))
        threading.Thread(target=th, daemon=True).start()

    run_btn.configure(command=start)
    cancel_btn.configure(command=cancel_ev.set)
    cancel_btn.state(["disabled"])
    poll()


def _path_row(frm, r, label, var, kind, hint="", filetypes=None, on_set=None):
    from tkinter import ttk, filedialog
    ttk.Label(frm, text=label).grid(row=r, column=0, sticky="w", pady=2)
    ttk.Entry(frm, textvariable=var).grid(row=r, column=1, sticky="ew", padx=6)

    def browse():
        if kind == "dir":
            v = filedialog.askdirectory()
        else:
            v = filedialog.askopenfilename(filetypes=filetypes or [("All", "*.*")])
        if v:
            var.set(os.path.normpath(v))
            if on_set:
                on_set()
    ttk.Button(frm, text="Browse...", command=browse).grid(row=r, column=2)
    if hint:
        ttk.Label(frm, text=hint, foreground="#666").grid(row=r + 1, column=1, sticky="w",
                                                          padx=6)
        return r + 2
    return r + 1


DEFAULT_EPSG = "7850"      # GDA2020 / MGA zone 50


def _epsg_row(frm, r, V, what):
    from tkinter import ttk
    ttk.Label(frm, text="Output CRS (EPSG)").grid(row=r, column=0, sticky="w", pady=2)
    o = ttk.Frame(frm)
    o.grid(row=r, column=1, columnspan=2, sticky="w", padx=6)
    ttk.Entry(o, textvariable=V["epsg"], width=8).pack(side="left")
    ttk.Label(o, text=f"  {what}; default 7850 = GDA2020 / MGA zone 50, any EPSG code",
              foreground="#666").pack(side="left")
    return r + 1


def _sdk_rows(frm, r, V):
    from tkinter import ttk
    status = ttk.Label(frm, text="", foreground="#555")

    def refresh(*_):
        ok, msg = sdk_status(V["sdk"].get().strip())
        status.configure(text=("OK: " if ok else "[!] ") + msg,
                         foreground="#1a6b1a" if ok else "#a00")
    r = _path_row(frm, r, "Phase One SDK folder", V["sdk"], "dir", on_set=refresh)
    status.grid(row=r, column=1, columnspan=2, sticky="w", padx=6)
    V["sdk"].trace_add("write", refresh)
    refresh()
    return r + 1


def _settings_rows(frm, r, V):
    from tkinter import ttk, filedialog
    ttk.Label(frm, text="Development settings").grid(row=r, column=0, sticky="w", pady=2)
    cb = ttk.Combobox(frm, textvariable=V["preset"], state="readonly",
                      values=list(PRESETS) + [CUSTOM])
    cb.grid(row=r, column=1, sticky="ew", padx=6)

    def browse():
        v = filedialog.askopenfilename(filetypes=[("Phase One settings", "*.json *.txt"),
                                                  ("All", "*.*")])
        if v:
            V["custom"].set(os.path.normpath(v))
            V["preset"].set(CUSTOM)
    ttk.Button(frm, text="File...", command=browse).grid(row=r, column=2)
    lab = ttk.Label(frm, text="", foreground="#666")
    lab.grid(row=r + 1, column=1, columnspan=2, sticky="w", padx=6)

    def show(*_):
        if V["preset"].get() == CUSTOM:
            lab.configure(text=V["custom"].get() or "choose a settings file (Phase One "
                                                     "ConvertConfig JSON)")
        elif V["preset"].get().startswith("Neutral"):
            lab.configure(text="no noise reduction or sharpening: closest to the sensor "
                               "data, with its colour noise; bit-repeatable")
        elif "HR_clip" in V["preset"].get():
            lab.configure(text="APPN's HR_clip settings: noise like Phase One's own JPEGs, "
                               "lossless output (repeatable to +-2/255)")
        else:
            lab.configure(text="APPN's template: brighter, highlight/shadow recovery, "
                               "colour noise reduction only")
    V["preset"].trace_add("write", show)
    show()
    return r + 2


def _log_area(frm, r):
    import tkinter as tk
    from tkinter import ttk
    bar = ttk.Frame(frm)
    bar.grid(row=r, column=0, columnspan=3, sticky="ew", pady=(8, 4))
    run_btn = ttk.Button(bar, text="Run")
    run_btn.pack(side="left")
    cancel_btn = ttk.Button(bar, text="Cancel")
    cancel_btn.pack(side="left", padx=6)
    pbar = ttk.Progressbar(bar, mode="determinate")
    pbar.pack(side="left", fill="x", expand=True, padx=6)
    box = tk.Text(frm, height=10, wrap="word", font=("Consolas", 9))
    box.grid(row=r + 1, column=0, columnspan=3, sticky="nsew")
    frm.rowconfigure(r + 1, weight=1)
    ttk.Label(frm, text=CREDIT, foreground="#666", wraplength=780, justify="left",
              font=("", 8)).grid(row=r + 2, column=0, columnspan=3, sticky="w",
                                 pady=(6, 0))
    return run_btn, cancel_btn, pbar, box


def _chunk_combo(frm, r, V):
    from tkinter import ttk
    ttk.Label(frm, text="Chunk").grid(row=r, column=0, sticky="w", pady=2)
    cb = ttk.Combobox(frm, textvariable=V["chunk"], width=30)
    cb.grid(row=r, column=1, sticky="w", padx=6)

    def fill(*_):
        try:
            import plot_extractor as pe
            labels = pe.read_psx_chunks(V["psx"].get().strip())[0]
            cb["values"] = [c if isinstance(c, str) else str(c) for c in labels]
        except Exception:
            cb["values"] = []
    cb.bind("<Button-1>", fill)
    return r + 1


# ---------------------------------------------------------------------------
# Tab: IIQ conversion
# ---------------------------------------------------------------------------
def build_convert_tab(parent):
    import tkinter as tk
    from tkinter import ttk
    c = _cfg_load()
    g = lambda k, d="": c.get(f"conv_{k}", d)
    frm = ttk.Frame(parent, padding=10)
    frm.columnconfigure(1, weight=1)
    V = {"iiq": tk.StringVar(value=g("iiq")), "rec": tk.BooleanVar(value=g("rec", True)),
         "out": tk.StringVar(value=g("out")),
         "fmt": tk.StringVar(value=g("fmt", "TIFF (lossless)")),
         "bits": tk.StringVar(value=g("bits", "8")),
         "quality": tk.StringVar(value=g("quality", "95")),
         "lzw": tk.BooleanVar(value=g("lzw", True)),
         "scale": tk.StringVar(value=g("scale", "1.0")),
         "preset": tk.StringVar(value=_preset_or_default(g("preset", PRESET_DEFAULT))),
         "custom": tk.StringVar(value=g("custom")),
         "psx": tk.StringVar(value=g("psx")), "chunk": tk.StringVar(value=g("chunk")),
         "epsg": tk.StringVar(value=g("epsg") or DEFAULT_EPSG),
         "sdk": tk.StringVar(value=g("sdk"))}
    r = 0
    ttk.Label(frm, text="Convert Phase One IIQ raw images to TIFF, GeoTIFF or JPEG",
              font=("", 10, "bold")).grid(row=r, column=0, columnspan=3, sticky="w")
    r += 1
    r = _path_row(frm, r, "IIQ folder *", V["iiq"], "dir",
                  "only read - nothing is written into this folder")
    ttk.Checkbutton(frm, text="include subfolders (mirrored in the output)",
                    variable=V["rec"]).grid(row=r, column=1, sticky="w", padx=6)
    r += 1
    r = _path_row(frm, r, "Output folder *", V["out"], "dir")
    o = ttk.Frame(frm)
    o.grid(row=r, column=1, columnspan=2, sticky="w", padx=6, pady=2)
    ttk.Label(frm, text="Output").grid(row=r, column=0, sticky="w")
    ttk.Combobox(o, textvariable=V["fmt"], state="readonly", width=34,
                 values=["TIFF (lossless)", "GeoTIFF (lossless, georeferenced)",
                         "JPEG"]).pack(side="left")
    ttk.Label(o, text="  bit depth").pack(side="left")
    ttk.Combobox(o, textvariable=V["bits"], state="readonly", width=4,
                 values=["8", "16"]).pack(side="left")
    ttk.Label(o, text="  JPEG quality").pack(side="left")
    ttk.Spinbox(o, from_=50, to=100, textvariable=V["quality"], width=5).pack(side="left")
    ttk.Checkbutton(o, text="LZW", variable=V["lzw"]).pack(side="left", padx=6)
    ttk.Label(o, text="scale").pack(side="left")
    ttk.Entry(o, textvariable=V["scale"], width=5).pack(side="left")
    r += 1
    ttk.Label(frm, text="16-bit: TIFF/GeoTIFF only.  Scale 1.0 = full size "
                        "(e.g. 0.3 for alignment-only copies).",
              foreground="#666").grid(row=r, column=1, columnspan=2, sticky="w", padx=6)
    r += 1
    r = _settings_rows(frm, r, V)
    r = _path_row(frm, r, "Metashape project", V["psx"], "file",
                  "GeoTIFF only: the project the frames were aligned in (read only)",
                  [("Metashape project", "*.psx")])
    r = _chunk_combo(frm, r, V)
    r = _epsg_row(frm, r, V, "GeoTIFF only: CRS of the georeferenced frames")
    r = _sdk_rows(frm, r, V)
    run_btn, cancel_btn, pbar, box = _log_area(frm, r)
    run_btn.configure(text="Convert")

    def work(log, progress, cancel):
        fmt = V["fmt"].get()
        convert_folder(dict(
            iiq_dir=V["iiq"].get().strip(), recursive=V["rec"].get(),
            out_dir=V["out"].get().strip(),
            fmt="geotiff" if fmt.startswith("Geo") else "jpg" if fmt.startswith("JPEG")
            else "tiff",
            bits=int(V["bits"].get() or 8), quality=int(V["quality"].get() or 95),
            lzw=V["lzw"].get(), scale=float(V["scale"].get() or 1.0),
            settings=preset_path(V["preset"].get(), V["custom"].get().strip()),
            psx=V["psx"].get().strip(), chunk=V["chunk"].get().strip(),
            epsg=V["epsg"].get().strip(), sdk=V["sdk"].get().strip()), log, progress, cancel)
    _runner(frm, V, "conv", box, pbar, run_btn, cancel_btn, work)
    frm.V, frm.run_btn, frm.log_box = V, run_btn, box
    return frm


# ---------------------------------------------------------------------------
# Tab: IIQ plot extraction
# ---------------------------------------------------------------------------
COMPRESS = ["LZW (lossless, smaller)", "None (uncompressed)"]


def compress_key(label):
    """'lzw' or 'none' from a COMPRESS label (or the key itself)."""
    return "none" if str(label).lower().startswith(("none", "unc")) else "lzw"


STYLES = {"Masked to the plot (developed pixels, outside black)": "mask",
          "Bounding box (developed pixels, no mask)": "bbox",
          "Rectified whole plot (stitched from frames, resampled)": "rectify"}


def run_iiq_extraction(p, log, progress, cancel):
    """The Extract-tab extraction with the pixels developed from the IIQ files."""
    import plot_extractor as pe
    for k, what in (("psx", "Metashape project"), ("plots", "plots file"),
                    ("rawdir", "IIQ folder"), ("outdir", "output folder")):
        if not p.get(k):
            raise RuntimeError(f"Choose the {what}.")
    assert_safe_output([p["rawdir"]], p["outdir"])
    os.makedirs(p["outdir"], exist_ok=True)
    return pe.run_extraction(p, log, progress, cancel)


def build_iiq_extract_tab(parent):
    import tkinter as tk
    from tkinter import ttk
    c = _cfg_load()
    g = lambda k, d="": c.get(f"iqx_{k}", d)
    frm = ttk.Frame(parent, padding=10)
    frm.columnconfigure(1, weight=1)
    V = {"psx": tk.StringVar(value=g("psx")), "chunk": tk.StringVar(value=g("chunk")),
         "plots": tk.StringVar(value=g("plots")), "idf": tk.StringVar(value=g("idf")),
         "iiq": tk.StringVar(value=g("iiq")), "out": tk.StringVar(value=g("out")),
         "preset": tk.StringVar(value=_preset_or_default(g("preset", PRESET_DEFAULT))),
         "custom": tk.StringVar(value=g("custom")),
         "style": tk.StringVar(value=g("style", list(STYLES)[0])),
         "margin": tk.StringVar(value=g("margin", "0.10")),
         "fmt": tk.StringVar(value=g("fmt", "GeoTIFF (geotagged, lossless)")),
         "bits": tk.StringVar(value=g("bits", "8")),
         "compress": tk.StringVar(value=g("compress", COMPRESS[0])),
         "workers": tk.StringVar(value=g("workers", "4")),
         "limit": tk.StringVar(value=g("limit", "0")),
         "epsg": tk.StringVar(value=g("epsg") or DEFAULT_EPSG),
         "sdk": tk.StringVar(value=g("sdk"))}
    r = 0
    ttk.Label(frm, text="Extract plot images directly from the IIQ raw files",
              font=("", 10, "bold")).grid(row=r, column=0, columnspan=3, sticky="w")
    r += 1
    r = _path_row(frm, r, "Metashape project *", V["psx"], "file",
                  "aligned project (opened read only); it may be aligned on JPEGs or on "
                  "downscaled TIFFs - coordinates are scaled to the IIQ size",
                  [("Metashape project", "*.psx")])
    r = _chunk_combo(frm, r, V)
    idcb = ttk.Combobox(frm, textvariable=V["idf"], width=30)

    def fill_ids(*_):
        try:
            import plot_extractor as pe
            _pl, _wkt, fields = pe.read_plot_file(V["plots"].get().strip())
            idcb["values"] = list(fields)
            if not V["idf"].get() and fields:
                V["idf"].set(pe.pick_id_field(fields) or fields[0])
        except Exception:
            pass
    r = _path_row(frm, r, "Plot grid *", V["plots"], "file",
                  "default: the grid made in tab 1 or edited in tab 2; Browse to use another",
                  [("Plots", "*.shp *.geojson *.json *.gpkg"), ("All", "*.*")], fill_ids)
    ttk.Label(frm, text="Plot ID field").grid(row=r, column=0, sticky="w", pady=2)
    idcb.grid(row=r, column=1, sticky="w", padx=6)
    idcb.bind("<Button-1>", fill_ids)
    r += 1
    r = _path_row(frm, r, "IIQ folder *", V["iiq"], "dir",
                  "the raw .IIQ files of the flight (only read); matched to the cameras "
                  "by file name")
    r = _path_row(frm, r, "Output folder *", V["out"], "dir")
    r = _settings_rows(frm, r, V)
    ttk.Label(frm, text="Crop").grid(row=r, column=0, sticky="w", pady=2)
    o = ttk.Frame(frm)
    o.grid(row=r, column=1, columnspan=2, sticky="w", padx=6)
    ttk.Combobox(o, textvariable=V["style"], state="readonly", width=48,
                 values=list(STYLES)).pack(side="left")
    ttk.Label(o, text="  margin (m)").pack(side="left")
    ttk.Entry(o, textvariable=V["margin"], width=6).pack(side="left")
    r += 1
    o2 = ttk.Frame(frm)
    o2.grid(row=r, column=1, columnspan=2, sticky="w", padx=6, pady=2)
    ttk.Label(frm, text="Output").grid(row=r, column=0, sticky="w")
    ttk.Combobox(o2, textvariable=V["fmt"], state="readonly", width=30,
                 values=["GeoTIFF (geotagged, lossless)", "TIFF (lossless)",
                         "PNG (lossless)", "JPEG"]).pack(side="left")
    ttk.Label(o2, text="  bits").pack(side="left")
    ttk.Combobox(o2, textvariable=V["bits"], state="readonly", width=3,
                 values=["8", "16"]).pack(side="left")
    ttk.Label(o2, text="  TIFF compression").pack(side="left")
    ttk.Combobox(o2, textvariable=V["compress"], state="readonly", width=22,
                 values=COMPRESS).pack(side="left")
    r += 1
    o2 = ttk.Frame(frm)
    o2.grid(row=r, column=1, columnspan=2, sticky="w", padx=6, pady=2)
    ttk.Label(o2, text="parallel plots").pack(side="left")
    ttk.Spinbox(o2, from_=1, to=16, textvariable=V["workers"], width=4).pack(side="left")
    ttk.Label(o2, text="  test: first").pack(side="left")
    ttk.Spinbox(o2, from_=0, to=10000, textvariable=V["limit"], width=6).pack(side="left")
    ttk.Label(o2, text="plots (0 = all)").pack(side="left")
    r += 1
    ttk.Label(frm, foreground="#555", justify="left", wraplength=640, text=(
        "Frame per plot: a photo that holds the WHOLE plot (the most nadir/central of "
        "those); if no photo holds the whole plot, the photo holding the largest part of "
        "it. The chosen box is then developed from that photo's IIQ. "
        "Bits: 8 (0-255) or 16 (0-65535, from the same development; Masked or Bounding "
        "box, TIFF/GeoTIFF/PNG). Compression: LZW and uncompressed hold identical pixels; "
        "LZW files are smaller.")).grid(
        row=r, column=1, columnspan=2, sticky="w", padx=6, pady=(2, 4))
    r += 1
    r = _epsg_row(frm, r, V, "CRS of the GeoTIFF georeference")
    r = _sdk_rows(frm, r, V)
    run_btn, cancel_btn, pbar, box = _log_area(frm, r)
    run_btn.configure(text="Extract")

    def work(log, progress, cancel):
        p = dict(psx=V["psx"].get().strip(), chunk=V["chunk"].get().strip(),
                 plots=V["plots"].get().strip(), id_field=V["idf"].get().strip(),
                 rawdir=V["iiq"].get().strip(), outdir=V["out"].get().strip(),
                 margin=V["margin"].get().strip() or "0", n_per_plot="1",
                 limit=V["limit"].get().strip() or "0", do_raw=True, do_ortho=False,
                 crs_override="auto", raw_mode=STYLES.get(V["style"].get(), "mask"),
                 format=V["fmt"].get(), to_srgb=False,
                 bits=V["bits"].get(), compress=compress_key(V["compress"].get()),
                 workers=V["workers"].get().strip() or "4", auto_fit=False,
                 frame_choice="whole_first", raw_source="iiq",
                 gcp_epsg=V["epsg"].get().strip(),
                 iiq_settings=preset_path(V["preset"].get(), V["custom"].get().strip()),
                 iiq_sdk=V["sdk"].get().strip())
        run_iiq_extraction(p, log, progress, cancel)
    _runner(frm, V, "iqx", box, pbar, run_btn, cancel_btn, work)
    frm.V, frm.run_btn, frm.log_box = V, run_btn, box   # main window sets V["plots"]
    return frm


# ---------------------------------------------------------------------------
# Tab: Acknowledgements
# ---------------------------------------------------------------------------
ACK_SECTIONS = [
    ("IIQ plot extraction: approach by the APPN node at The University of Queensland", [
        "The IIQ plot extraction (tab 5) follows the PhaseOne_Pipeline of the Australian Plant "
        "Phenomics Network (APPN), The University of Queensland node (APPN-UQ; "
        "appn@uq.edu.au; MIT licence, Copyright (c) 2025 APPN-UQ). Their method: convert the "
        "IIQ files to 0.3-scale TIFFs, align them in Metashape, back-project the four corners "
        "of every plot into the aligned photos, keep the photos that hold all four corners, "
        "pick the one whose corners are closest to the image centre, and develop that "
        "bounding box at full size from the IIQ with HR_clip (Phase One ImageSDK, their "
        "config_hr settings).",
        "Code taken from APPN-UQ: iiq_tools/src/p1tool.cpp is adapted from their IIQ_Proc "
        "tools (IIQ2TIFF.cpp, IIQHRclip.cpp); the presets 'APPN HR_clip' and 'APPN template' "
        "are their config_hr settings. The full licence is in THIRD_PARTY_NOTICES.md."]),
    ("What we changed in their method", [
        "Partial plots: their pipeline drops a plot when no single photo holds all four "
        "corners. We take a photo holding the whole plot first and, when none does, the photo "
        "holding the largest part of it (the coverage is recorded in the manifest), or "
        "stitch the whole plot from several photos (rectified style).",
        "Plot shape: the whole polygon (all vertices, plus a margin) is projected with the "
        "full camera model, not only four corners; any polygon and any Plot ID field work.",
        "Exact scale: the box is scaled by the true frame ratio (e.g. 12768/3830 = 3.3337, "
        "not 1/0.3), and projects aligned on full-size JPEGs work as well as on 0.3-scale "
        "TIFFs, so no extra conversion is needed for alignment.",
        "Output: masked to the plot, bounding box or rectified; TIFF, GeoTIFF (georeferenced), "
        "PNG or JPEG; files named by Plot ID (no overwriting when plots share a photo); "
        "camera EXIF plus the source photo and development settings in the tags; "
        "manifest.csv and candidates.csv.",
        "p1tool rewritten for batches: one process develops many boxes, Unicode paths, "
        "8- or 16-bit, TIFF / JPEG / raw RGB output, frame info; a Neutral preset (no noise "
        "reduction or sharpening) next to their two presets.",
        "Safety and checks: output inside the IIQ folder is refused, every IIQ used is "
        "re-checked unchanged after the run, and the QA tab re-develops each box from the "
        "IIQ and compares it with the saved crop.",
        "IIQ conversion (tab 3): whole IIQ frames to TIFF / GeoTIFF / JPEG, GeoTIFFs "
        "georeferenced from the aligned cameras."]),
    ("Phase One", [
        "Pixels in tabs 3 and 5 are developed by the Phase One ImageSDK (Copyright Phase One "
        "A/S), used under Phase One's SDK licence. The SDK runtime is not part of the public "
        "source code. The Phase One SDK use at DPIRD follows the approval obtained for the "
        "DPIRD IIQ converter (github.com/appndpird/P1ImageSDK-IIQ-Converter, Bipul Neupane)."]),
    ("Everything else", [
        "All other parts of PlotExtractor - plot grid generation and editing (tabs 1 and 2), "
        "plot extraction from the Phase One JPEGs or the orthomosaic (tab 4), training tiles "
        "and row detection (tab 6, after the PhenoFuse row method) and QA (tab 7) - were "
        "developed by DPIRD for the APPN node in Western Australia."]),
]


def build_ack_tab(parent):
    import tkinter as tk
    from tkinter import ttk
    frm = ttk.Frame(parent, padding=10)
    txt = tk.Text(frm, wrap="word", relief="flat", padx=14, pady=10,
                  font=("Segoe UI", 10), background="white")
    sb = ttk.Scrollbar(frm, command=txt.yview)
    txt.configure(yscrollcommand=sb.set)
    sb.pack(side="right", fill="y")
    txt.pack(side="left", fill="both", expand=True)
    txt.tag_configure("h", font=("Segoe UI", 12, "bold"), foreground="#0c1a2b",
                      spacing1=10, spacing3=4)
    txt.tag_configure("p", spacing3=6, lmargin1=10, lmargin2=10)
    txt.tag_configure("b", spacing3=4, lmargin1=10, lmargin2=26)
    nl = chr(10)
    txt.insert("end", "Acknowledgements" + nl, ("h",))
    for head, paras in ACK_SECTIONS:
        txt.insert("end", head + nl, ("h",))
        bullets = head.startswith("What we changed")
        for para in paras:
            txt.insert("end", ("•  " if bullets else "") + para + nl,
                       ("b" if bullets else "p",))
    txt.configure(state="disabled")
    frm.text = txt
    return frm
