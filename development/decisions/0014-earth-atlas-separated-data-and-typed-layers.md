# 0014. The Earth Atlas: Separated Data and Typed Layers

- **Status:** Accepted, provisionally: built while the user was away; open to revision on review. Amended by [0015](0015-earth-atlas-models-imagery-detail.md).
- **Date:** 2026-10-22
- **Decided by:** the user asked for an application for the Earth that integrates geography, oceanography, meteorology and climate, taking Windy as an inspiration but not as a design; asked that geographic data be kept separate, because many kinds of data will be shown on the Earth (climate simulations, climate data, our own results, economic and population data); and asked for high quality. Claude chose the design below.

## 1. Context

Windy shows one weather variable at a time over a map, with moving particles for the wind and a timeline. It is a fine viewer of forecasts, but its design is bound to its subject:

- its layers are named after weather variables;
- its colours are baked into the data it serves;
- its map is for navigation, not for comparing areas;
- it has no place for data about people, or for our own results.

The laboratory's purpose differs. The atlas must show any quantity that has a place on the Earth, from any field, at a stated time and with stated units and sources, and must let new data in without changes to its code. The base map itself (coasts, countries, names) is one dataset among others, and must not carry the values of any other.

## 2. Decision

1. **Name and place.** The application is the *Earth Atlas*, at `applications/geography/earth-atlas/`, a tool.
   - It lives under geography because its own work is cartographic: projecting, composing layers, joining data to places.
   - The data it shows come from every field through the data library.
   - The name pairs it with the Celestial Atlas.
2. **Data kept apart, joined only in the application.** Each kind of data is its own product (`data/products/`), with its own provenance:
   - **Geography** (`geography.py`, from Natural Earth). Geometry, identifiers, names in 26 languages, and cartographic ranks only. Natural Earth's population, GDP and income columns are dropped on purpose.
   - **Fields** (`field_grid.py`). Variables on a regular longitude–latitude grid, over time, in the units of their source, with CF standard names, on a normalised grid layout. Any netCDF file can become a field product, including our own model output (`--input`).
   - **Indicators** (`country_indicators.py`, from the World Bank). One value per economy and year, keyed by codes, with no geometry.
   - **Events** (`event_catalog.py`, from the USGS). Things at a place and an instant.
3. **Explicit joins.** Indicators are joined to countries by key in `src/core/join.js`.
   - Each country carries several systems of identifiers. They are tried in a stated order: World Bank, ISO alpha-3 (Natural Earth's corrected column, then the plain one), Natural Earth's own, then ISO alpha-2.
   - Unmatched records are reported, not hidden. The test suite checks that only Gibraltar, absent from the 1:50 million map, is unmatched.
4. **Styles belong to the application.** Colour maps, ranges and scales are chosen in `src/layers/catalogue.js`:
   - by CF standard name for fields, so that the same quantity from any source looks the same;
   - from the data's own range for a variable the catalogue does not know, so that new data are shown at once.
5. **Layers typed by the shape of their data**, not by subject:
   - a *surface* (one quantity over the surface, a field or a value per country);
   - a *flow* (a vector field, shown by particles);
   - *lines* (isolines with highs and lows);
   - *events* (points in time);
   - a *base* (relief, coasts, boundaries, rivers, graticule, names, day and night).

   One of each can be shown together, for example sea-surface temperature, the wind, isobars and earthquakes.
6. **One instant, stated for every layer.** The atlas shows one instant in UTC.
   - Each layer shows its data nearest to that instant. Fields are interpolated linearly in time, and means are taken at the nearest step.
   - Each layer states which time it shows, and how far that time is from the instant shown.
   - Indicators, which are yearly, have their own year.
7. **The data library grows**, as [ADR 0012](0012-selectable-parameterised-datasets.md) intended:
   - new kinds `griddap` (subsets of ERDDAP grids) and `http` (URL templates);
   - a parameter type `time`, with expressions relative to the present or to the end of a dataset;
   - manifest fields `domain`, `credit` and `attributes`;
   - a reader for netCDF classic files in `data/formats/`, using the standard library only.
8. **One renderer for globe and maps.**
   - A single fragment shader inverts the projection at each pixel (Equal Earth, equirectangular, or the sphere for the globe), then samples the base, the data and the Sun there.
   - It interpolates exactly as `src/core/grid.js` does, so the colour under the cursor and the value in the readout agree.
   - Equal Earth is the default flat map because it is equal-area, as maps of quantities per area require.
9. **Physics stated, not imitated.**
   - **Particles** move on the sphere by the midpoint rule with parallel transport, and their speed-up is stated in the legend.
   - **Wave particles** move at the deep-water group speed g T / 4π.
   - **Day and night** are drawn as the bands of civil, nautical and astronomical twilight, from a low-precision Sun tested against JPL Horizons.
10. **Provisional computation in JavaScript**, as in [ADR 0011](0011-first-application-celestial-atlas.md). It is a stand-in for `core/` and is tested in `tests/run-tests.mjs`.
11. **Shared look.** The atlas shares the Celestial Atlas's colours (plum, with Dusty Orchid #CA8DBE) and serif type, so that the laboratory's tools read as one family. Its main colour can be changed on request.

## 3. Alternatives Considered

- **Windy's model:** one weather layer at a time, a navigation projection, colours rendered on the server. Rejected for the reasons in §1. Its strengths are kept: particles, a timeline, and a readout at a point.
- **Joining indicators to geometry when building products.** Rejected. It would bind every indicator to one map, and hide the choices made in the join.
- **Map tiles (XYZ pyramids)** with a library such as MapLibre GL. Deferred. At the present resolutions (0.25° to 1°) whole grids fit easily in memory, and tiles would add a server and a build step. They become necessary for kilometre-scale data.
- **Mercator** for the flat map. Rejected, because it inflates high latitudes several-fold, which misleads every map of a quantity per area.
- **A new bridge for the Earth system** in `core/domains/` and `applications/`. Not taken without the user's agreement (CLAUDE.md §6.4); the climate bridge was the other candidate.
- **Importing the Celestial Atlas's modules.** Not done, because applications must not depend on each other (CLAUDE.md §3.3). A few small parts are therefore duplicated:
  - the local server;
  - the time-scale arithmetic;
  - the palette.

  Moving shared web code into `visualization/` is proposed for the user's agreement.

## 4. Consequences

- **New data appear without changes to code.**
  - A dataset that becomes a product of a known kind appears in the atlas.
  - A field with a known CF standard name takes that quantity's style. An unknown one gets a style fitted to its range.
- **The data library now reaches beyond astronomy:**
  - ERDDAP: forecasts, analyses and relief;
  - the World Bank;
  - the USGS;
  - Natural Earth.
- **The tests check the core against references:**
  - subsolar points against JPL Horizons, 1955 to 2045;
  - Equal Earth's equal area and axis ratio against its paper;
  - GMST against Meeus;
  - interpolation, contours, rasters, particles, joins and plausibility checks on the real products.
- **The core's building blocks have a later home.** Projections, rasters, joins, isolines and particle advection belong in `core/domains/geography` and physis once the C++ core exists.
- **Limits:**
  - fields are drawn at their own resolution (1° for the default forecast);
  - Natural Earth's boundaries are its de facto depiction;
  - climate-model and reanalysis archives are not registered yet. ERA5 needs a Copernicus account; CMIP data are served by ESGF.
