# data: Acquisition, Registry and Provenance

*Status: draft.* An edge of the architecture. Written in Python first, using the standard library (with numpy in product builders), with C++ readers later where performance matters. Running an AI model needs more (onnxruntime and ecCodes), kept in a separate environment ([ADR 0022](../development/decisions/0022-ai-forecast-cyclones-and-disasters.md)):

```
python -m venv .venv
.venv/Scripts/python -m pip install -r data/products/requirements-models.txt
```

## 1. Purpose

This directory finds, downloads or queries, verifies and records the external data the project uses. It then hands the data to the core and the applications as local files. The core never fetches data itself.

A user chooses **which** data to take and **how much** of it ([ADR 0012](../development/decisions/0012-selectable-parameterised-datasets.md)):

- a dataset is chosen by its manifest;
- the part and the amount are chosen by the manifest's parameters and named profiles.

Data of different kinds are kept apart ([ADR 0014](../development/decisions/0014-earth-atlas-separated-data-and-typed-layers.md)):

- a base map carries geometry, identifiers and names, and no values about the places it draws;
- values come as their own products, and applications join them to places by stated keys.

## 2. Modules

| Module | Scope | Status |
|---|---|---|
| [`registry/`](registry/) | One manifest per dataset: source, licence, checksum or query, parameters, profiles. The schema is in its README. | draft |
| [`acquisition/`](acquisition/) | `fetch.py`: fixed files, catalogue queries, gridded subsets and web requests into `data/cache/`, with `retrieval.json`. `erddap.py`: a client for ERDDAP's griddap. `ncss.py`: a client for THREDDS's NetCDF Subset Service. `opendap.py`: a client for OPeNDAP (DAP2) servers, for strided grids or the column beneath a place. `remotezip.py`: members of zip archives on the web, by range requests. `on_demand.py`: requests made while an application waits (§3.1). | draft |
| [`catalogs/`](catalogs/) | `tap.py`: a minimal client for the IVOA Table Access Protocol, used for ESA Gaia, CDS VizieR and SIMBAD | draft |
| [`formats/`](formats/) | Readers and writers: netCDF classic (`netcdf3.py`). FITS, netCDF-4/HDF5 and Zarr are planned. | draft |
| [`products/`](products/) | Derived artefacts, built reproducibly from registered datasets, and their index `index.json` (§4). `ai_forecast.py` runs an AI weather model to make one. | draft |
| `provenance` | Records of which data and which code produced each result. For now this is `retrieval.json` for each retrieval, and a `source` block in each product. | seed |

## 3. Usage

From the repository root:

```
python data/acquisition/fetch.py --list                          # what can be fetched, with parameters
python data/acquisition/fetch.py                                 # every fixed-file dataset
python data/acquisition/fetch.py gaia-dr3 --profile small        # a catalogue query, by profile
python data/acquisition/fetch.py gaia-dr3 --set g_max=9.5        # ... or with parameters overridden
python data/acquisition/fetch.py gfs                             # the latest weather forecast (ERDDAP)
python data/acquisition/fetch.py worldbank-indicators --profile life-expectancy
python data/products/star_catalog.py                             # build products (§4)
python data/products/field_grid.py --source gfs
python data/acquisition/fetch.py gfs-levels aigfs                # forecast models aloft (THREDDS)
python data/products/imagery.py --source gibs-daily             # yesterday's Earth from space
python data/products/body_maps.py                                # maps of the planets and the sky
python data/acquisition/fetch.py koppen-geiger terraclimate-normals   # climate classes and normals
python data/products/class_raster.py --source koppen-geiger
python data/products/field_grid.py --source terraclimate-normals
python data/acquisition/fetch.py gdacs                           # disasters of the last month, and cyclone tracks
python data/products/disaster_catalog.py
python data/acquisition/fetch.py pangu-weather gfs-analysis-0p25 # an AI model's weights, and an analysis to start it
.venv/Scripts/python data/products/ai_forecast.py --source gfs-analysis-0p25 --hours 168
python data/acquisition/fetch.py pb2002 ngl-midas gcmt crust1          # the solid Earth
python data/products/plate_model.py
python data/products/station_velocities.py
python data/products/mechanism_catalog.py
python data/products/crust_model.py
```

The cache, `data/cache/`, is not under version control.

### 3.1 Requests on Demand

Some data cannot be fetched ahead of time, such as the ocean beneath a place a user has just chosen: one column out of millions. A manifest offers itself for such requests with an `on_demand` block (the parameters a request may set, a product builder, a typical duration). `on_demand.py` runs requests as jobs, checks their parameters against the manifest, shares identical requests, runs two at a time, and caches retrievals and products as any others ([ADR 0023](../development/decisions/0023-ocean-column-and-data-on-demand.md)). A local server exposes them to its pages; the Earth Atlas's does at `/request`:

```
http://127.0.0.1:8790/request?dataset=hycom-column&lat=32&lon=135     starts the request, or finds it
http://127.0.0.1:8790/request/<job id>                                its state; when done, the product's path
```

The same column from the command line:

```
python data/acquisition/fetch.py hycom-column --set lat=32 --set lon=135
python data/products/ocean_column.py --profile custom-<hash>
```

## 4. Products

Products are written to `data/cache/products/` and listed by kind in `index.json`, which applications read.

| Kind | Builder | Content | Format |
|---|---|---|---|
| `stars` | `star_catalog.py` | Stars from HYG or Gaia: position, magnitude, colour, distance | `.bin` (float32 records) and `.json` |
| `galaxies` | `galaxy_catalog.py` | Galaxies of the Local Volume (UNGC) and 2MRS groups, in supergalactic coordinates | `.json`, or `.bin` and `.json` |
| `geography` | `geography.py` | The base map: countries (with several systems of codes and names in 26 languages), land, lakes, glaciers, ice shelves, sea names, coastlines, boundaries, rivers, places. No thematic values. | `.bin` (float32 longitude–latitude pairs) and `.json` (features referring to runs of pairs) |
| `fields` | `field_grid.py`, `ai_forecast.py` | Variables on a regular longitude–latitude grid, over time: forecasts, analyses, relief, climate normals, the laboratory's own AI forecasts, or our own netCDF files. A forecast that has been verified records its errors (`verification`). Rows run north to south; global grids start at 180° W; values in their source units. A climatology's times are the months of its period's middle year, and it records the period. | `.bin` (int16, a linear scale per variable) and `.json` (grid, times, CF metadata, statistics) |
| `classes` | `class_raster.py` | A category per cell of a regular grid, such as climate classes, by period and scenario; the legend (code, name, the authors' colour) and the share of land in each class | `.bin` (uint8, 0 for no value) and `.json` |
| `indicators` | `country_indicators.py` | A value per economy and year, keyed by World Bank codes with ISO alpha-2; aggregates kept apart | `.json` |
| `events` | `event_catalog.py`, `disaster_catalog.py`, `mechanism_catalog.py` | Events at a place and an instant. Earthquakes (`earthquake`): time, longitude, latitude, depth, magnitude. Disasters (`disaster`): start and last update, a point, the kind (tropical cyclone, earthquake, flood, volcano, drought, wildfire), GDACS's alert level, severity, countries. Mechanisms (`mechanism`): the centroid, Mw, the six moment-tensor elements over M0, and the catalogue's nodal planes. | `.json` |
| `plates` | `plate_model.py` | Tectonic plates: their outlines (rings of points) and the steps of their boundaries, each with its class, its left and right plate, the polarity of subduction and the source's own relative motion. Geometry only: the rotations are a model, kept with the code. | `.bin` (float32 longitude–latitude pairs) and `.json` |
| `stations` | `station_velocities.py` | GPS stations: place, height, east, north and up velocity with uncertainties (mm/yr), first and last epoch, in a stated frame | `.json` |
| `series` | `station_series.py` | One station's daily positions, east, north and up (mm from the first day) with their sigmas. Made on demand, in `products/series/`, not listed in the index. | `.json` |
| `crust` | `crust_model.py` | The layers of the crust on a grid: each layer's top, P and S speeds and density, as scaled 16-bit integers; with a field product of the Moho's depth and the crust's and sediments' thicknesses | `.bin` (int16) and `.json` |
| `columns` | `ocean_column.py` | The ocean beneath a place: currents and temperature by depth and time, at a model's grid point, with its latest run and deepest level. Made on demand, kept in `products/columns/`, not listed in the index. | `.json` |
| `tracks` | `disaster_catalog.py` | Moving things, as positions over time. Tropical cyclones (`cyclone`): positions every 6 or 12 hours, each marked observed or forecast by the responsible warning centre. | `.json` |
| `imagery` | `imagery.py`, `body_maps.py` | Equirectangular images with a role. Shaded relief of the land (`relief-shading`), with its grey levels of water and flat land, for detail at close range. The Earth (`day`, `night-lights`): true colour of a day (gaps filled from the day before), the Blue Marble, the Black Marble, with the provider's service for detail. Other bodies (`body`), Saturn's rings (`rings`, a radial strip) and the whole sky in galactic coordinates (`sky`, its orientation checked against bright stars). | `.jpg` or `.png`, and `.json` |

## 5. Datasets

### 5.1 Registered

| Dataset | Kind | Provider | Use |
|---|---|---|---|
| `jpl-de440s` | files | NASA JPL / NAIF | Positions of the Sun, Moon and planets, 1849–2150 |
| `hyg-v41` | files | Hipparcos, Yale, Gliese (compiled by D. Nash) | Stars to about magnitude 7.5 |
| `d3-celestial-constellations` | files | O. Frohn | Constellation figures and names |
| `gaia-dr3` | tap | ESA Gaia (DPAC); distances from Bailer-Jones et al. 2021 | Stars with parallaxes, brightest first |
| `ungc` | tap | Karachentsev et al. 2013, via CDS VizieR | Galaxies of the Local Volume, with tidal indices |
| `tully-2mass-groups` | tap | Tully 2015 and Huchra et al. 2012, via CDS VizieR | Galaxy groups of the 2MASS Redshift Survey |
| `natural-earth-50m` | files | Natural Earth, v5.1.2 | The base map of the Earth |
| `etopo-2022` | griddap | NOAA NCEI, via CoastWatch ERDDAP | Relief of land and sea floor, 15″ source grid |
| `gfs` | griddap | NOAA NCEP, via PacIOOS ERDDAP | Weather: analyses and forecasts of wind, temperature, pressure, precipitation, humidity, radiation |
| `oisst` | griddap | NOAA NCEI (OISST v2.1), via CoastWatch ERDDAP | Daily sea-surface temperature, anomaly and sea ice |
| `ww3` | griddap | NOAA NCEP (WaveWatch III), via PacIOOS ERDDAP | Ocean waves: height, period, direction |
| `worldbank-indicators` | http | The World Bank | Indicators by country and year: population, economy, health, energy, emissions |
| `usgs-earthquakes` | http | USGS ComCat | Earthquakes by time, place, depth and magnitude |
| `gfs-levels` | ncss | NOAA NCEP, via Unidata THREDDS | GFS aloft (850, 500, 250 hPa) and over the column (cloud, water vapour, CAPE) |
| `aigfs` | ncss | NOAA NCEP and EPIC, via Unidata THREDDS | NOAA's AI global forecast model, near the surface and aloft |
| `aigefs` | ncss | NOAA NCEP and EPIC, via Unidata THREDDS | NOAA's AI ensemble: mean and spread |
| `gibs-daily` | http | NASA GIBS | Yesterday's Earth from VIIRS, true colour |
| `gibs-static` | http | NASA GIBS | The Blue Marble by day, the Black Marble by night |
| `gibs-relief` | http | NASA GIBS, from ASTER GDEM (METI and NASA) | Shaded relief of the land from 30 m elevations: a global preview, and detail at close range |
| `solar-system-maps` | files | Solar System Scope (INOVE), from NASA, JPL and ESA data, via Wikimedia Commons | Maps of the Sun, the Moon, the planets, Saturn's rings and the Milky Way |
| `koppen-geiger` | files (members of a remote archive) | Beck et al. 2023, via figshare | Köppen–Geiger climate classes at 0.1°, 1901–2099, past and under four scenarios |
| `terraclimate-normals` | opendap | Abatzoglou et al. 2018 (TerraClimate), via the Northwest Knowledge Network | Monthly normals of 1991–2020: daily maximum and minimum temperature, precipitation |
| `gfs-analysis-0p25` | http | NOAA NCEP, via the NOMADS GRIB filter | The GFS analysis at 0.25° and 13 levels, to start an AI model |
| `pangu-weather` | files | Bi et al. 2023 (Huawei Cloud), via ECMWF | The AI model Pangu-Weather, 24 h and 6 h, as ONNX (CC BY-NC-SA 4.0: never committed, never used commercially) |
| `gdacs` | http (following links) | GDACS (UN and EC JRC) | Disasters of the last month with alert levels; tracks and official forecast tracks of tropical cyclones |
| `hycom-column` | opendap (a column), on demand | US Navy, ESPC-D-V02 (HYCOM), via HYCOM.org | The ocean beneath a place: currents and temperature at 40 levels to 5000 m, two days back and four ahead |
| `pb2002` | files | Bird 2003, via fraxen/tectonicplates (ODC-By) | The outlines of 52 tectonic plates and the 5,819 steps of their boundaries |
| `ngl-midas` | http | Nevada Geodetic Laboratory | GPS velocities of about 22,000 stations (MIDAS), in IGS20 or a plate's frame |
| `ngl-station` | http, on demand | Nevada Geodetic Laboratory | One GPS station's daily positions |
| `gcmt` | http | Global CMT Project (LDEO) | The moment tensors of about 70,000 earthquakes since 1976 |
| `crust1` | files | Laske et al. 2013 (IGPP, UC San Diego) | CRUST1.0: the crust's layers, speeds and densities on 1° cells |

### 5.2 Planned

| Dataset | Provider | Use |
|---|---|---|
| Earth-orientation parameters, leap seconds | IERS | Time scales and reference frames |
| Tide-gauge records and predictions | NOAA CO-OPS | Tidal analysis and validation |
| Cosmicflows-4 distances and peculiar velocities | Tully et al. 2023, via VizieR | Flows of galaxies |
| ERA5 reanalysis | ECMWF / Copernicus (needs an account) | The atmosphere's past, hourly since 1940 |
| CMIP6 climate projections | ESGF | Simulated climates, by scenario |
| Ocean currents (OSCAR or CMEMS) | NASA / Copernicus | Surface currents, as a flow |
| Landsat, Sentinel-2 | USGS / NASA, ESA | Remote sensing of the land surface |
| OpenStreetMap, GHSL, WorldPop | OSM contributors, JRC, WorldPop | Settlement, land use, population on a grid |

## 6. Licences

Each dataset keeps its own licence, which its manifest records. The project's MIT licence does not cover datasets.
