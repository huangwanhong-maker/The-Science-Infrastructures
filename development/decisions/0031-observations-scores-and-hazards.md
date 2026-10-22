# 0031. Observations, Scores and Hazards: the Earth Atlas as a Forecasting System

- **Status:** Accepted, provisionally: built at the user's request; open to revision on review. Extends [ADR 0022](0022-ai-forecast-cyclones-and-disasters.md), [ADR 0027](0027-earth-atlas-plugins.md).
- **Date:** 2026-10-23
- **Decided by:** the user asked for what weather prediction, disaster prediction and geographical analysis still lacked, as Claude had listed it, and for the user's own sensors to serve rendering and prediction. Claude chose the design below. The research arc (an inference foundation, the C++ toolchain, a primitive-equation model with a real assimilation cycle) was deferred by the user's decision.

## 1. Context

The Earth Atlas showed and compared other centres' forecasts well, and ran small models of its own. It held no observations, kept no past forecasts and scored them only against their own analyses: it could not tell how good a forecast had been. Its hazards were alerts issued elsewhere; its geography stopped at countries. Pangu-Weather, its own AI model, has weights that forbid commercial use.

## 2. Decision

1. **Observations at places** are a kind of product of their own (`observations/1`): one row per station and time, SI units under CF standard names, null where missing, values beyond the extremes ever observed nulled and counted, kept as a rolling window merged at each rebuild. METAR (about 5,000 aerodromes) and NDBC buoys are the first sources; IBTrACS gives best tracks of tropical cyclones.
2. **The user's own sensors** have a home in `observation/sensors`: a store that checks readings and converts them to SI on arrival, fed through the API, from files (a logger's, followed as it grows), or by a simulated sensor that is the template for a device. Each source becomes an observations product, drawn and read as any other, and posted readings reach listeners at once.
3. **Forecasts are kept and scored.** Each run of each forecast is archived, named by its start; once its weather has happened it is scored: 500 hPa height against analyses (RMSE and centred anomaly correlation, beside persistence), 2 m temperature, pressure and waves against station reports, cyclones against best tracks. Each score names its truth, its definition and its count of cases.
4. **A schedule** fetches observations hourly and forecasts six-hourly, then archives and scores, one job at a time under an operating-system lock that the page's own runs also take.
5. **ECMWF's open forecasts** join NOAA's: the IFS (the strongest global model), AIFS (its AI model, CC BY 4.0, open for commercial use) and the ensemble's probabilities, beside GFS under every quantity.
6. **Hazards are models in the domains, run as plug-ins**: fire weather (the Canadian FWI System), heat (wet-bulb temperature, heat index, excess over the normal), heavy rain (against the month's normal), tsunami first arrivals over the bathymetry, aftershock probabilities, and the settling of volcanic ash. Each plug-in says it is not a warning service.
7. **Geography gains people and regions**: GHS-POP on a 0.25° grid, countries and provinces as regions on 0.1°, zonal statistics as exact overlap integrals on the sphere, and drainage from elevation. Exposure (people under heat or gales, by country and province) joins hazard to people.
8. **Assimilation of real observations**: the optimal-interpolation plug-in analyses at the time shown, the background interpolated in time, toward station reports or the user's sensors, with gross-error rejection and the diagnostics of an analysis.
9. **The page** draws observations as dots coloured by a variable, the scorecard as charts against lead, and plug-in results by region with outlines. On a public address, its data recipe runs only for the holder of a token.

## 3. Alternatives Considered

- **Observations as events.** Events are things that happen once; a station reports again and again, and its variables are many. A kind of its own keeps both simple.
- **Scores against each model's own analyses only.** Cheap, but it flatters each model; GFS analyses serve all, and the scores say so.
- **ERA5 for climatology and truth.** Better, but behind an account; NCEP/NCAR Reanalysis 1 means serve the anomaly correlation openly.
- **Storm surge, regional models, radar nowcasting, river networks.** Named as next steps: surge needs the tides domain; regional models and radar need sources with clear licences.

## 4. Consequences

- **The archive grows** by about 50 MB a cycle (only the scored variables are kept), 6 GB over thirty days.
- **Scores rest on few cases** until the archive has grown; each says its count, and the page shows sparse points faintly.
- **A failed refresh is silent**: when a source is down, the fetcher keeps the cached retrieval and the job reports success. Flagging stale retrievals is a next step.
- **Height corrections** of station temperatures await a model orography.
