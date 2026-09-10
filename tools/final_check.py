from osgeo import gdal

gdal.UseExceptions()
p = (r"F:\Ibrahim's Workspace2\Plant Emergency Dataset UWA"
     r"\2025_NUE+RootPhenotyping_I_DPIRD-20263108\RawCrops\plot_7043.tif")
ds = gdal.Open(p)
print("EPSG:", ds.GetSpatialRef().GetAuthorityCode(None))
print("GT:", tuple(round(v, 6) for v in ds.GetGeoTransform()))
print("GCPs:", ds.GetGCPCount())
print("resid:", ds.GetMetadataItem("PE_AFFINE_RESID"))
