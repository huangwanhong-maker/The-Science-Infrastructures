# formats: Readers for Data Formats

*Status: draft.* Part of [`data/`](../), an edge of the architecture.

## 1. Purpose

Datasets arrive in many formats. This module reads them into plain arrays and dictionaries, so that product builders need not know how a file is laid out. Readers use the Python standard library. They may use numpy where it is available, but must not require it. Writers, which store data fetched in other forms (such as OPeNDAP responses) in a standard format, may require numpy.

## 2. Modules

| Module | Format | Status |
|---|---|---|
| [`netcdf3.py`](netcdf3.py) | netCDF classic (CDF-1) and 64-bit offset (CDF-2): dimensions, attributes, variables. `values()` unpacks `scale_factor` and `add_offset` and turns fill values into NaN; `array()` does the same with numpy. Record variables, and files written with a streaming record count, are supported. `write()` stores fixed-size variables in the 64-bit offset format, with typed attributes (needs numpy). | draft |
| netCDF-4 / HDF5 | Climate-model output (CMIP), reanalyses | planned |
| FITS | Astronomical images and tables | planned |
| Zarr | Cloud-native arrays | planned |

## 3. References

- Unidata, *NetCDF Classic and 64-bit Offset File Formats*, https://docs.unidata.ucar.edu/netcdf-c/current/file_format_specifications.html
- *NetCDF Climate and Forecast (CF) Metadata Conventions*, version 1.11: standard names, units, `cell_methods`.
