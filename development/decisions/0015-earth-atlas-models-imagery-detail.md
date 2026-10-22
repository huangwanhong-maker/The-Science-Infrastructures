# 0015. Forecast Models, a Laboratory Model, Imagery and Detail in the Earth Atlas

- **Status:** Accepted, provisionally: built while the user was away; open to revision on review. Amends [ADR 0014](0014-earth-atlas-separated-data-and-typed-layers.md). Amended by [ADR 0018](0018-fine-relief-at-close-range.md) (relief detail).
- **Date:** 2026-10-22
- **Decided by:** the user asked for weather prediction models, finer textures on the 3D globe, richer geography, and creativity. Claude chose the design below.

## 1. Context

The first Earth Atlas showed one model (GFS at the surface), relief at 10′, and no imagery. A forecast is better understood as one of several. Comparing models, including one we can run and inspect ourselves, is both a learning activity and a research activity. The globe at close range needs imagery and detail finer than one global texture can hold.

## 2. Decision

1. **A new acquisition kind, `ncss`.** It reads THREDDS NetCDF Subset Service collections (Unidata's forecast models), by groups of variables at stated levels.
   - Registered so far:
     - `gfs-levels`: GFS at 850, 500 and 250 hPa, and cloud, water vapour and CAPE;
     - `aigfs`: NOAA's AI model;
     - `aigefs`: NOAA's AI ensemble, mean and spread.
   - A product built from several files names each variable with its level (`Temperature_isobaric@850hPa`).
2. **Models as a dimension of the catalogue.**
   - Layers of the same quantity at the same level share a key. The panel offers their models as a choice.
   - A layer can be shown minus another model's, as a difference on a diverging scale fitted to the whole forecast.
   - The place card gives every model's value and draws their courses together.
   - Ensemble spread is a quantity of its own.
3. **The laboratory's own model.**
   - **The model.** The barotropic vorticity equation (Charney, Fjørtoft and von Neumann 1950), solved by the spectral transform method at T42. It runs in a Web Worker from the 500 hPa wind of the GFS at the time shown.
   - **Its output.** A computed product, the same interface as a library product. It is drawn, contoured, compared and verified like any other.
   - **Verification.** The panel reports the root-mean-square error of height over 20–80° N, against the GFS's analyses and against persistence.
   - **Tests.** The model is tested against the Rossby–Haurwitz wave, an exact solution.
4. **A new product kind, `imagery`.** Images of the whole Earth, from NASA GIBS (`gibs-daily`, `gibs-static`):
   - yesterday's VIIRS true colour, its gaps filled from the day before;
   - the Blue Marble;
   - the Black Marble at night.

   They are drawn as base styles, and as lights on the night side.
5. **Detail where the user looks.**
   - An imagery product records how to ask its provider for a part of it: its `live` service, a WMS request with the same layer and date.
   - When every pixel of the view falls on the Earth, the atlas requests the region in view at the screen's resolution and draws it over the global image.
   - This amends ADR 0014's rule that applications read only products. An application may also use a *live service declared by a product*, for the same data at another resolution. The registry thus stays the single record of sources.
6. **Forecast clouds.** The GFS total cloud cover can be drawn as a layer of cloud over the base, at the time shown.

## 3. Alternatives Considered

- **ECMWF open data and DWD ICON.** Deferred. They are distributed as GRIB2 with packings that need libraries we do not have (CCSDS). THREDDS serves NOAA's models as netCDF, which our reader handles.
- **A shallow-water model** instead of the barotropic one. Deferred. The barotropic model is the historical first step, has an exact solution to test against, and fits in the page. The shallow-water equations are the natural next step (Williamson et al. 1992 give its test suite).
- **Map tiles for detail.** A single WMS request for the region in view is simpler, and serves the same purpose at the zoom levels a globe needs.
- **Fetching detail through the local server.** Unnecessary for GIBS, which allows browsers to fetch it directly (CORS). It would be needed for relief detail from ERDDAP, which does not; that is left for later.

## 4. Consequences

- **Three model sources and one model of our own.** They can be compared at any place and time, and our model is verified as its builders verified theirs.
- **The tests grow to 78 checks.** They cover the spectral transforms (exact to 10⁻¹⁴), the Rossby–Haurwitz wave (7 × 10⁻¹⁰), conservation of energy and enstrophy, and linear balance.
- **The globe shows imagery at any zoom.** Detail at 250 m per pixel comes from VIIRS when close. Black Marble lights show at night.
- **Providers' terms.**
  - THREDDS can take minutes to cut a subset. The fetcher waits for it and retries.
  - GIBS imagery carries NASA's acknowledgement, given in the atlas's credits.
