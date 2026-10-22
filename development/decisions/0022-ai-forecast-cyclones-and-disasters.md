# 0022. An AI Weather Model Run in the Laboratory, Tropical Cyclones Tracked, and Disaster Alerts

- **Status:** Accepted, provisionally: built while the user was away; open to revision on review.
- **Date:** 2026-10-22
- **Decided by:** the user suggested weather and disaster prediction, by any algorithm, state-of-the-art AI included if open source. They also asked whether the atlas's predictions were only displays of others' simulations. Claude chose the design below.

## 1. Context

The Earth Atlas showed forecasts from four sources (ADR 0015, ADR 0019):
- NOAA's GFS, AIGFS and AIGEFS, which it displays but does not run;
- the laboratory's barotropic and shallow-water models, which it runs in the page.

The laboratory's own models are real predictions but simple: a single layer, with no moisture, no surface and no baroclinic dynamics. The user's question was fair: the laboratory itself predicted only with its simplest models. Nor did the atlas say anything about hazards. It showed no tropical cyclones as such, no warnings, and no disasters other than earthquakes.

## 2. Decision

1. **A state-of-the-art AI model, run by the laboratory.** Pangu-Weather (Bi et al. 2023) runs on the laboratory's own computer, from an analysis the data library fetches. Its forecast becomes a field product like any other.
   - **Initial state.** The GFS analysis at 0.25°, from NOAA's NOMADS GRIB filter: five variables at Pangu-Weather's 13 levels, and four at the surface. It is the dataset `gfs-analysis-0p25`, whose cycle is a time parameter moved down to 00, 06, 12 or 18 UTC.
   - **Model.** The 24 h and 6 h models, as ONNX files (`pangu-weather`, a `files` dataset with their checksums). They are composed by hierarchical temporal aggregation, which uses the fewest steps ([theory note](../../theory/meteorology/learned-forecast-models.md), Proposition 1). They run on the CPU through onnxruntime at about 85 s a step, so a week takes 28 steps and about 40 minutes.
   - **Output.** `data/products/ai_forecast.py` writes the forecast every 6 hours at 1°: sea-level pressure, 2 m temperature, 10 m wind, 500 hPa height, 850 hPa temperature, and wind at 850, 500 and 250 hPa. With `--verify`, the product also records the root-mean-square error of 500 hPa height over 20–80° N against the GFS's later analyses, with persistence, as for the laboratory's other models.
   - **Environment.** The model needs numpy, onnxruntime and ecCodes, which the rest of the data library does not. They live in a separate environment, `.venv`, pinned in `data/products/requirements-models.txt`. The fetcher and the other builders keep to the standard library and numpy.
   - **Licence.** The weights are CC BY-NC-SA 4.0. They are fetched into the cache and never committed. The forecasts are for learning and research, and their products say so.
   - *Placement, provisional.* The forecast step is pure computation and belongs in `core/inference` (learning) once the C++ core exists. Until then it lives with the product builders, as the laboratory's own models live in the Earth Atlas's `src/core/`.
2. **Tropical cyclones found in every forecast.** `src/core/cyclones.js` detects and links cyclones in sea-level pressure after TempestExtremes (Ullrich and Zarzycki 2017), with its thresholds for grids of 1° ([theory note](../../theory/meteorology/tropical-cyclone-tracking.md)).
   - The warm-core test is replaced by requirements that the track deepen to 1000 hPa and begin over the sea within 25° of the equator. The sea is judged by the relief product.
   - Any field product with sea-level pressure at more than four times is tracked, in the page, in about 0.2 s. A new model's cyclones therefore appear without new code.
   - Each model's tracks take its colour, and take the name of the official cyclone they follow.
3. **Observed and official tracks, and disaster alerts, from GDACS.** The dataset `gdacs` fetches the last month's events, of six kinds:
   - tropical cyclones, earthquakes, floods, volcanoes, droughts and wildfires;
   - each with its alert level, which assesses the potential need for international assistance.

   It also fetches the geometry of each tropical cyclone. `data/products/disaster_catalog.py` builds two products:
   - an `events` product of kind `disaster`, drawn by alert level (colour) and kind (size), from the start of each event until a week after its last update;
   - a product of a new kind, `tracks`: positions every 6 or 12 hours, observed up to the last observation and forecast beyond it by the responsible warning centre, as GDACS relays them. The observed track is drawn solid and the official forecast dashed.
4. **Requests that follow links.** An `http` dataset may declare `follow` rules. A rule names a JSON file already fetched, the items in it, a filter on their fields, and templates for a URL and a file name. The fetcher downloads what the items point to, so the GDACS cyclone geometries are named by the event list.
   - Names come from fetched content, so the fetcher refuses a file name that would leave the dataset's folder, and any URL that is not HTTPS.
   - A failed linked request is skipped and does not fail the retrieval.
5. **Hazards in the Earth Atlas.**
   - A *Cyclones* section offers the official tracks and the tracks found in each model. With several chosen, the forecasts of one cyclone can be compared on the map.
   - At the time shown, each cyclone is labelled once, with its name and each model's central pressure in that model's colour.
   - *Disaster alerts* is a choice beside the earthquakes.
   - The panel's verification table adds the AI forecast's errors beside the laboratory's models for the same start.

## 3. Alternatives Considered

- **Other AI models.** GraphCast (Google DeepMind), FourCastNet (NVIDIA) and Aurora (Microsoft) were considered. Pangu-Weather was chosen first for two reasons:
  - it is distributed as ONNX files that onnxruntime runs on a CPU, with no deep-learning framework;
  - it takes exactly the 13 pressure levels the GFS analysis provides.

  The others can be added the same way, as further models in `ai_forecast.py`, each verified as this one is.
- **Running on the GPU.** The laptop's GPU has 4 GB of memory, and was not tried. The CPU makes a week's forecast in about 40 minutes, which suffices for forecasts made a few times a day.
- **Only displaying NOAA's AI forecasts.** These are already shown. The user asked for prediction, and running a model makes its initial state, its version and its verification the laboratory's own.
- **Tracking cyclones in the product builders.** The tracker runs in the page instead. It then serves every field product, including ones added later, and its thresholds can be changed and the result seen at once. It is fast enough at 1°.
- **IBTrACS, or each warning centre directly.** IBTrACS is the archive of best tracks. Its records of current storms are provisional and updated only every few days, and it carries no forecasts. The warning centres (NHC, JTWC, JMA and others) each publish in their own form for their own ocean. GDACS gathers current cyclones, their official forecasts and other hazards in one interface.

## 4. Consequences

- **Forecasts.** The laboratory now makes a forecast that is close to the state of the art, on its own computer. It is verified against analyses that arrive afterwards, beside the laboratory's simple models and persistence.
  - From the analysis of 2026-10-19 00:00 UTC, the 500 hPa height error over 20–80° N was 6.6 m at +24 h, 12.7 m at +48 h and 18.3 m at +72 h. The shallow-water model's was 43, 74 and 100 m, and persistence's 81, 108 and 126 m.
  - Typhoon Champi's track in that forecast stayed within 49 km of the observed centre at +24 h, 58 km at +48 h and 79 km at +60 h. Its central pressure was about 45 hPa higher than the GFS's on the same grid, as its authors' report of weak cyclones leads one to expect.
- **Cyclones.** The tracks of one cyclone can be compared across the models and the official forecast. This shows the disagreement between models, which is the uncertainty of the forecast.
- **Tests.** The Earth Atlas has 114 checks, adding the tracker's detection, rejection, linking and filters, and interpolation across the date line.
- **Data.**
  - The model's weights are two files of 1.18 GB each, in the cache.
  - An analysis is 57 MB.
  - The GDACS event list and the cyclone geometries are fetched into the cache and not redistributed, since GDACS states no licence for reuse. They are credited on the map.
- *Caveat.*
  - GDACS's cyclone and earthquake alerts are automatic, and do not replace official warnings. The atlas says so in its legend.
  - Intensities at 1° are well short of a cyclone's true intensity, and Pangu-Weather underestimates them further (Bi et al. 2023).
  - The forecasts are for learning and research. They are not for decisions about safety.
