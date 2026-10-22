# formats: Readers for Data Formats

*Status: draft.* Part of [`data/`](../), an edge of the architecture.

## 1. Purpose

Datasets arrive in many formats. This module reads them into plain arrays and dictionaries, so that product builders need not know how a file is laid out. Readers use the Python standard library. They may use numpy where it is available, but must not require it, except to deliver pixels: `tiff.py` parses headers and decodes with the standard library, and needs numpy only for the arrays of a window. Writers, which store data fetched in other forms (such as OPeNDAP responses) in a standard format, may require numpy.

## 2. Modules

| Module | Format | Status |
|---|---|---|
| [`netcdf3.py`](netcdf3.py) | netCDF classic (CDF-1) and 64-bit offset (CDF-2): dimensions, attributes, variables. `values()` unpacks `scale_factor` and `add_offset` and turns fill values into NaN; `array()` does the same with numpy. Record variables, and files written with a streaming record count, are supported. `write()` stores fixed-size variables in the 64-bit offset format, with typed attributes (needs numpy). | draft |
| [`tiff.py`](tiff.py) | GeoTIFF and Cloud-Optimised GeoTIFF: windows of any image or overview, read locally or by HTTP range requests; strips or tiles, LZW and Deflate, the horizontal and floating-point predictors, integers of 8–64 bits and floats; GeoKeys (EPSG code, pixel as area or point), tie point and pixel scale. `encode_grid()` writes a float32 grid in EPSG:4326 (GeoTIFF 1.1, Deflate, NaN as GDAL_NODATA) that GDAL and QGIS read. The UTM zones on WGS 84 by the transverse Mercator projection to 5 nm (Karney 2011). Needs numpy for pixels. | draft |
| [`export.py`](export.py) | Products in standard formats for other software: a field as GeoTIFF, netCDF (CF-1.11, packed as stored, with its grid mapping), CSV or JSON, with a box and a stride; observations, events and tracks as GeoJSON (RFC 7946). The API's export endpoint calls it. | draft |
| netCDF-4 / HDF5 | Climate-model output (CMIP), reanalyses | planned |
| FITS | Astronomical images and tables | planned |
| Zarr | Cloud-native arrays | planned |

## 3. References

- Unidata, *NetCDF Classic and 64-bit Offset File Formats*, https://docs.unidata.ucar.edu/netcdf-c/current/file_format_specifications.html
- *NetCDF Climate and Forecast (CF) Metadata Conventions*, version 1.11: standard names, units, `cell_methods`, grid mappings.
- Adobe, *TIFF Revision 6.0* (1992); Adobe, *Photoshop TIFF Technical Note 3* (2005), the floating-point predictor.
- OGC 19-008r4, *GeoTIFF Standard* 1.1 (2019); OGC 21-026, *Cloud Optimized GeoTIFF* (2023).
- Karney, C. F. F. (2011), Transverse Mercator with an accuracy of a few nanometers, *J. Geodesy* 85, 475–485.
- RFC 7946, *The GeoJSON Format* (2016).
