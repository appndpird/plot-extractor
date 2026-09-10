import sys
import rasterio

path = r"F:\Ibrahim's Workspace2\Plant Emergency Dataset UWA\2025_NUE+RootPhenotyping_I_DPIRD-20263108\RawCrops\plot_1001.tif"
if len(sys.argv) > 1:
    path = sys.argv[1]

with rasterio.open(path) as ds:
    print("driver:", ds.driver)
    print("crs:", ds.crs)
    print("transform:", ds.transform)
    print("bounds:", ds.bounds)
    print("size:", ds.width, "x", ds.height, "bands:", ds.count, ds.dtypes)
    print("profile:", {k: v for k, v in ds.profile.items() if k not in ("crs", "transform")})
    print("dataset tags:", ds.tags())
    print("ns=IMAGE_STRUCTURE tags:", ds.tags(ns="IMAGE_STRUCTURE"))
    print("gcps:", len(ds.gcps[0]) if ds.gcps and ds.gcps[0] else 0, "gcp crs:", ds.gcps[1])
    print("rpcs:", ds.rpcs is not None)
