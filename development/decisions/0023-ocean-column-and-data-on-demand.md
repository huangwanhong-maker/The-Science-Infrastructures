# 0023. The Ocean Beneath a Place, and Data Requested on Demand

- **Status:** Accepted, provisionally: built at the user's request; open to revision on review.
- **Date:** 2026-10-22
- **Decided by:** the user asked whether the Earth Atlas could show, for a clicked place, a record like a ship's current profiler's (eastward and northward current, and echo strength, against depth and time). Claude found the sources and chose the design below.

## 1. Context

A current profiler's record exists only where an instrument is: along research ships' tracks (archived by JASADCP, or streamed by a few ships) and at a few moorings. For any place, the only source is an ocean model. The global models that are openly served and can be subset are these:
- **HYCOM / ESPC-D-V02** (US Navy): 1/12°, 40 levels, every 3 hours, analyses about ten days back and a forecast a week ahead. Public release, through HYCOM.org's THREDDS server.
- **NOAA's RTOFS:** its OPeNDAP service on NOMADS was retired (Service Change Notice 25-81). Only whole files remain.
- **Copernicus Marine:** requires an account.

HYCOM's server opens one file per time step, so a request costs about 10 s per time step and variable, whatever its spatial extent. A strided global subset of all levels timed out. Fetching columns ahead of time is therefore impossible, and a column must be fetched when it is wanted.

The data library fetched only ahead of time, from the command line. Applications read finished products.

## 2. Decision

1. **Requests on demand, through the data library.** A manifest may offer an `on_demand` block, which gives:
   - the parameters a request may set;
   - a product builder;
   - a typical duration.

   `data/acquisition/on_demand.py` runs such requests as jobs:
   - parameters are checked against the manifest's types and bounds;
   - identical requests share a job;
   - two jobs run at once;
   - each job's log is kept;
   - retrievals and products are cached as usual.

   The Earth Atlas's local server exposes this at `/request?dataset=…` and `/request/<job>`. It is bound to 127.0.0.1, as before.
2. **Columns from OPeNDAP.** An `opendap` manifest with a `column` block takes, for one variable of each dataset, the column at the grid point nearest a place:
   - The fetcher reads the coordinates (latitude, longitude on the circle, depth and CF times).
   - It selects times on multiples of a step between two times, and depths down to a limit.
   - It requests the times in short runs, four requests at a time. This is about three times faster on HYCOM's server, which opens one file per step.
   - It writes one netCDF file per variable, with dimensions time × depth.
3. **The `hycom-column` dataset and the `columns` product kind.**
   - The dataset takes eastward and northward current and temperature, from two days back to four ahead, every 6 hours, at all 40 levels. A first request takes about four minutes.
   - `data/products/ocean_column.py` writes a `columns/1` product. It holds the times, depths and values, the grid point, the deepest level with water, and the model's latest run. A run time later than the retrieval is a mislabelled file of the server's collection, and is ignored.
   - Columns are kept in `products/columns/` and are not listed in the index, since there may be any number of them.
4. **The record in the atlas.** For a place at sea, the place card offers the column. While it is being fetched, the card counts the time. When it is ready, the card gives:
   - the current at the surface and at 500 m;
   - the mixed-layer depth;
   - the water temperature at the surface and at 500 m.

   The column itself is drawn over the map as a profiler's record (`src/ui/column.js`):
   - **Panels:** eastward current, northward current and temperature, with isotherms and the mixed layer.
   - **Marks:** a strip of the Sun's height at the place; the start of the forecast; the time shown; the sea floor.
   - **Interaction:** a reading under the pointer; depths 0–1000 m, or down to the sea floor.

   The computations are in `src/core/ocean.js`: the inertial period, the direction of flow, the mixed-layer depth (0.2 °C from 10 m; de Boyer Montégut et al. 2004), and interpolation within the column.

## 3. Alternatives Considered

- **Fetching a coarse three-dimensional product ahead of time.** One global time step at all levels and 2° timed out after five minutes. A week of them for three variables would take hours and download gigabytes, almost all of it never looked at.
- **Fetching from the page directly.** HYCOM's server sends no CORS header. A proxy in the page's server would bypass the data library's provenance and caching. The data library does the fetching, as everywhere else.
- **Echo strength, as in a profiler's record.** No model has it. The caption says what a profiler sees that a model cannot: the plankton that migrate at dawn and dusk.
- **Copernicus Marine's time-series stores,** which answer a point in seconds. They need an account. That remains open, for the user to decide.

## 4. Consequences

- **Any place at sea** can be seen in depth, from two days back to four ahead, without anything fetched in advance.
- **A pattern for other data that are too large to fetch ahead,** such as a station's time series, which the solid-Earth layers will use.
- **The server becomes an edge with state:** jobs, threads and a cache. It remains local.
- **Tests:** the Earth Atlas has 124 checks, adding:
  - the inertial period (one sidereal day at 30°);
  - directions of flow;
  - the mixed layer of a constructed profile;
  - exact interpolation of bilinear fields;
  - nothing beneath the sea floor.

  The fetcher's selection of the grid point (across 0°/360°), depths and times was checked on constructed axes.
- *Caveat.*
  - HYCOM's analyses lag about two days, so the present is already a forecast.
  - The server's record of runs is inconsistent.
  - A column's first request takes minutes.
  - A page served without the laboratory's server cannot request columns.
