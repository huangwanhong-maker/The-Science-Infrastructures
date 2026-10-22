# registry: Dataset Manifests

*Status: draft.* Each external dataset the project uses has one manifest here. A manifest says:
- what the data are;
- where they come from;
- what licence they carry;
- how to acquire them.

Manifests are read by [`../acquisition/fetch.py`](../acquisition/fetch.py). It downloads or queries the data into `data/cache/`, and writes a provenance record (`retrieval.json`) beside each result.

## 1. Kinds of Dataset

| Kind | Meaning | Example |
|---|---|---|
| `files` | Fixed files, downloaded as they are, or taken from a zip archive on the web. Each file's SHA-256 is recorded in the manifest and verified on every retrieval. | `jpl-de440s`, `natural-earth-50m`, `solar-system-maps`, `koppen-geiger`, `pangu-weather`, `pb2002`, `crust1` |
| `tap` | A query to an astronomical archive through the IVOA Table Access Protocol | `gaia-dr3`, `ungc`, `tully-2mass-groups` |
| `griddap` | A subset of a gridded dataset served by ERDDAP: which variables, which times, how fine a grid | `gfs`, `oisst`, `ww3`, `etopo-2022` |
| `http` | One or more web requests built from URL templates, and optionally those the results link to | `worldbank-indicators`, `usgs-earthquakes`, `gibs-daily`, `gibs-static`, `gibs-relief`, `gfs-analysis-0p25`, `gdacs`, `ngl-midas`, `ngl-station`, `gcmt` |
| `ncss` | Groups of variables of a forecast model's collection on a THREDDS server, by its NetCDF Subset Service: levels, times and grid spacing | `gfs-levels`, `aigfs`, `aigefs` |
| `opendap` | One variable from each of several datasets on an OPeNDAP (DAP2) server, every n-th point of the grid, or the column beneath a place, saved as netCDF classic | `terraclimate-normals`, `hycom-column` |

All kinds except `files` declare parameters and profiles, so that a user chooses which part of a dataset to take and how much of it ([ADR 0012](../../development/decisions/0012-selectable-parameterised-datasets.md)).

## 2. Common Fields

| Field | Meaning |
|---|---|
| `id` | Identifier; also the name of the cache folder |
| `kind` | `files` (the default), `tap`, `griddap`, `http`, `ncss` or `opendap` |
| `domain` | The field the data belong to: `astronomy`, `geography`, `oceanography`, `meteorology`, `geophysics`, `climate`. Applications group data by it. |
| `title`, `provider` | What the dataset is and who produced it |
| `credit` | A short attribution, shown on maps |
| `model` | For forecast models: the model's short name (GFS, AIGFS), by which applications offer a choice of model |
| `climatology` | For normals: `{period: "YYYY-MM-DD/YYYY-MM-DD"}`. The product's times become the months of a typical year (§9). |
| `format`, `coverage` | Units, frames, epochs, grids; what the data cover |
| `licence`, `citation` | Terms of use, and how to cite. The project's MIT licence does not cover data. |
| `used_by` | The applications or modules that depend on the dataset |

## 3. Fields of `files` Datasets

`files` is a list of `{url, file, size, sha256}`. A `null` checksum marks a file never verified; run `fetch.py --record` once to record it. Prefer URLs pinned to a release (for example a git tag), so that checksums stay valid. A file may carry further fields that its product builder reads; for example, `solar-system-maps` names the `body` each map shows, and `koppen-geiger` the `period` and `scenario` of each map. A map may also name a tile service for detail at close range (`live`: a WMTS template and its levels), which its product passes on.

A file with a `member` is that member of the zip archive at `url`. The fetcher reads the archive's directory and the member by HTTP range requests, without downloading the archive; the checksum is the member's.

## 4. Parameters and Profiles

These fields are shared by the parameterised kinds (`tap`, `griddap`, `http`, `ncss`, `opendap`).

| Field | Meaning |
|---|---|
| `parameters` | For each parameter: `type`, `default`, optional `min` and `max`, and a `description`. The type is `int`, `float`, `str` or `time`. A `str` parameter may have a `pattern` (a regular expression its whole value must match), as one that names a file on a server must. |
| `profiles` | Named sets of parameter values, each with a `description` |
| `default_profile` | The profile used when none is chosen |

A retrieval chooses a profile with `--profile` and may override any parameter with `--set name=value`.

**Parameters of type `time`** accept two forms:
- an ISO 8601 instant, read as UTC when it has no offset;
- an expression `now`, `today` or `last`, optionally followed by `+` or `−` and a number of hours (`h`) or days (`d`), for example `now-1d` or `last-29d`.

The anchors mean:
- `now` is the present, rounded down to the hour;
- `today` is 00:00 UTC;
- `last` is the dataset's last time (`griddap` only).

A `time` parameter may also carry `align_hours`: the resolved instant is then moved down to a multiple of that many hours from 00:00 UTC, as a forecast's cycle is (00, 06, 12, 18 UTC for 6).

The resolved instants are recorded in the provenance.

**Cache folders.** Results go to `data/cache/<id>/<profile>/`. When parameters are overridden, they go to `data/cache/<id>/custom-<hash>/`, named by a hash of the values.

**Refreshing.** A retrieval whose request is identical to the cached one is not repeated. Relative times make "the latest forecast" a moving target, so fetching again replaces it once a new request results.

## 5. Fields of `tap` Datasets

| Field | Meaning |
|---|---|
| `service` | Base URL of the TAP service |
| `query` | ADQL as a list of lines, with `{parameter}` placeholders |
| `file` | Name of the result file (CSV) |
| `maxrec_parameter` | The parameter passed to the service as `MAXREC` (default `limit`). Without it, many services truncate results silently. |

## 6. Fields of `griddap` Datasets

| Field | Meaning |
|---|---|
| `service`, `dataset` | The ERDDAP server (ending in `/erddap`) and the dataset identifier |
| `variables` | Comma-separated variable names, or a `{parameter}` |
| `axes` | The axes, in the dataset's order. Each is `{name, from, to, stride}`, with values in the axis's units. The time axis instead takes `{from, to, step_hours}`, with `time` parameters. A fixed axis (such as a depth of 0) has `from` equal to `to`. |
| `attributes` | Attributes to add to variables whose source lacks them, such as `units` or `standard_name` |
| `file` | Name of the result (netCDF classic, default `data.nc`) |

**Choosing times.** For the time axis, `fetch.py` reads the dataset's time values within the window. It takes those on the step, aligned to 00:00 UTC, or to the first time when none falls there, as for daily means stamped at 12:00. It then requests them with a single stride. The window is clamped to the dataset's coverage.

## 7. Fields of `http` Datasets

`requests` is a list of `{url, file}`. Each URL is a template with `{parameter}` placeholders. A `time` parameter may carry a format, as in `{start:%Y-%m-%dT%H:%M:%S}`.

`follow` is an optional list of rules for requests that the results name, fetched after `requests`:

| Field | Meaning |
|---|---|
| `from` | A JSON file already fetched |
| `items` | The key of the list in it, such as `features` |
| `where` | Fields an item must have, by dotted paths, such as `{"properties.eventtype": "TC"}` |
| `url`, `file` | Templates over the item, indexed as Python's `str.format` does, such as `{properties[url][geometry]}` |

The names come from fetched content. A file name that would leave the dataset's folder, and a URL that is not HTTPS, are refused. A linked request that fails is skipped and the retrieval continues.

## 8. Fields of `ncss` Datasets

| Field | Meaning |
|---|---|
| `service`, `dataset` | The NCSS endpoint (ending in `/ncss/grid`) and the collection's path, such as `grib/NCEP/GFS/Global_onedeg/Best` |
| `time_step_hours` | The collection's own time step; `step_hours` must be a multiple of it |
| `groups` | Requests, each a list of `variables` sharing a coordinate system, with an optional vertical level (`vertCoord`, in the collection's units, such as Pa) and its name (`level`, such as `500 hPa`); one file per group |
| `attributes` | Standard names and units for the variables, which GRIB-derived collections lack |

Times run from `start` to `end` every `step_hours`, aligned to 00:00 UTC. A product built from the files names each variable with its level, as in `Temperature_isobaric@850hPa`.

## 9. Fields of `opendap` Datasets

| Field | Meaning |
|---|---|
| `service` | The OPeNDAP base URL, such as a THREDDS server's `.../thredds/dodsC/` |
| `datasets` | One `{path, variable, file}` for each variable: the dataset's path under `service`, the variable to take, and the netCDF file to write |
| `attributes` | Standard names, long names, units and cell methods to set on the variables |
| `climatology` | For normals, their period (§2) |

The parameter `stride` takes every n-th point along the horizontal dimensions (`lat`, `lon` and their synonyms), and every point along the others. The subset is a point sample, not an area mean. The fetcher reads the dataset's structure and attributes (DDS and DAS), requests the hyperslab in binary (DAP2), and writes it with the variable's attributes and coordinates.

**Columns.** With a `column` block (`{lat, lon, depth, time, run}`: the names of the coordinates, `run` optional), the fetcher takes instead the column at the grid point nearest the parameters `lat` and `lon` (longitudes compared on the circle). It takes the depths down to `max_depth`, and the times from `start` to `end` on multiples of `step_hours` since 00:00 UTC; CF time units are read from the DAS. Times are requested in runs of at most four, `workers` requests at a time (default 4), and each variable is written as netCDF with dimensions (time, depth), times in seconds since 1970, and the grid point's coordinates as global attributes.

## 10. Requests on Demand

A manifest may offer itself for requests made while an application waits:

| Field | Meaning |
|---|---|
| `on_demand.parameters` | The parameters a request may set; the others keep their defaults |
| `on_demand.builder` | The module of `data/products` whose `build(folder, manifest)` makes the product, returning its path |
| `on_demand.seconds` | How long a request usually takes, to tell the user |
| `on_demand.description` | What a request gives |

See `data/acquisition/on_demand.py` and [ADR 0023](../../development/decisions/0023-ocean-column-and-data-on-demand.md).

## 11. Usage

```
python data/acquisition/fetch.py --list                               # datasets, parameters, profiles
python data/acquisition/fetch.py gaia-dr3 --profile nearby
python data/acquisition/fetch.py gaia-dr3 --set limit=50000 --set g_max=9.5
python data/acquisition/fetch.py gfs --profile glance
python data/acquisition/fetch.py oisst --set start=2026-01-01 --set end=2026-01-31
python data/acquisition/fetch.py worldbank-indicators --set indicator=SP.POP.65UP.TO.ZS   # share aged 65 and over
```
