"""Convert GCP-only RawCrops GeoTIFFs into QGIS-ready GeoTIFFs.

For each TIF: least-squares affine (col,row)->(E,N) through the embedded
EPSG:7850 GCPs, then in-place header edit: clear GCPs, set the affine as the
geotransform + CRS EPSG:7850. Pixels are never rewritten. The original GCPs
and the fit residual are preserved as metadata tags (PE_GCPS, PE_AFFINE_RESID).

Usage:  python fix_georef.py <file-or-folder> [more ...]
"""
import glob
import json
import os
import sys

import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()


def fit_affine(gcps):
    A = np.array([[g.GCPPixel, g.GCPLine, 1.0] for g in gcps])
    E = np.array([g.GCPX for g in gcps])
    N = np.array([g.GCPY for g in gcps])
    ce, *_ = np.linalg.lstsq(A, E, rcond=None)
    cn, *_ = np.linalg.lstsq(A, N, rcond=None)
    res = np.hypot(A @ ce - E, A @ cn - N)
    gt = (float(ce[2]), float(ce[0]), float(ce[1]),
          float(cn[2]), float(cn[0]), float(cn[1]))
    return gt, res


def fix(path):
    ds = gdal.Open(path, gdal.GA_Update)
    try:
        n = ds.GetGCPCount()
        gt0 = ds.GetGeoTransform(can_return_null=True)
        if n == 0:
            if gt0:
                return "already-georef"
            return "no-gcps-no-transform"
        gcps = ds.GetGCPs()
        gcp_srs = ds.GetGCPSpatialRef()
        epsg = None
        if gcp_srs is not None:
            code = gcp_srs.GetAuthorityCode(None)
            epsg = int(code) if code else None
        if epsg != 7850:
            return f"gcp-crs-not-7850 ({epsg})"
        if n < 3:
            return f"too-few-gcps ({n})"
        gt, res = fit_affine(gcps)
        # preserve the original GCPs + fit quality as metadata
        ds.SetMetadataItem("PE_GCPS", json.dumps(
            [[g.GCPLine, g.GCPPixel, g.GCPX, g.GCPY, g.GCPZ] for g in gcps]))
        ds.SetMetadataItem("PE_GCPS_EPSG", "7850")
        ds.SetMetadataItem("PE_AFFINE_RESID",
                           f"mean={res.mean()*1000:.2f}mm max={res.max()*1000:.2f}mm n={n}")
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(7850)
        ds.SetGCPs([], None)            # clear GCPs so the affine governs
        ds.SetSpatialRef(srs)
        ds.SetGeoTransform(gt)
        return f"ok resid_mean={res.mean()*1000:.1f}mm max={res.max()*1000:.1f}mm"
    finally:
        ds = None


def main(argv):
    files = []
    for a in argv:
        if os.path.isdir(a):
            files += sorted(glob.glob(os.path.join(a, "*.tif")))
        else:
            files.append(a)
    counts = {}
    worst = 0.0
    for i, f in enumerate(files, 1):
        try:
            r = fix(f)
        except Exception as e:
            r = f"ERROR {e}"
        key = r.split()[0]
        counts[key] = counts.get(key, 0) + 1
        if key != "ok" or len(files) == 1:
            print(f"{os.path.basename(f)}: {r}")
        elif "max=" in r:
            worst = max(worst, float(r.split("max=")[1].split("mm")[0]))
        if i % 50 == 0:
            print(f"  ... {i}/{len(files)}")
    print("\nSummary:", counts, f"worst max-residual {worst:.1f} mm")


if __name__ == "__main__":
    main(sys.argv[1:])
