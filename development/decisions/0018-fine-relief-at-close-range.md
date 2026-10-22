# 0018. Fine Relief at Close Range, and a Finer Base Map

- **Status:** Accepted, provisionally: built while the user was away; open to revision on review. Amends [ADR 0015](0015-earth-atlas-models-imagery-detail.md), which left relief detail for later.
- **Date:** 2026-10-22
- **Decided by:** the user asked for finer textures on the 3D globe. Claude chose the design below.

## 1. Context

The Earth Atlas painted its relief map at 4096 × 2048 pixels from ETOPO 2022 at 5′ (about 9 km). Closer than a country, mountains became soft blurs. On the steepest slopes the shading, tuned for a 10′ grid, turned black.

ADR 0015 brought detail for imagery from NASA GIBS, but not for relief. ETOPO's ERDDAP server does not allow browsers to fetch from it, and fetching through the local server was left for later.

## 2. Decision

1. **Relief detail from GIBS.**
   - GIBS serves grey shaded relief computed from the ASTER global digital elevation model (1″, about 30 m), with CORS. A new dataset, `gibs-relief`, and the imagery builder make it a product with the role `relief-shading`: a global preview and its `live` WMS service.
   - The builder measures two grey levels: water (the uniform grey of lakes, sea and gaps; 88) and flat land (144).
2. **The shader replaces coarse shading with fine.**
   - On the relief map at close range, the atlas fetches the region in view, at the screen's resolution, as PNG. A pixel of the water's grey is taken as water only if its whole neighbourhood has that grey, so that a shadow of the same grey is not (`waterMask` in `src/core/raster.js`).
   - On land, the shader divides the relief map's own shading out of its colour and applies the fine shading in the same form, 1 + (s − 1) k. Data drawn over the relief are shaded by the fine relief too.
   - The sea keeps the bathymetric shading of ETOPO, which ASTER does not have.
3. **A finer base map, painted in a worker.**
   - The relief map is painted at 8192 × 4096 (about 5 km to a pixel), so coasts, lakes and glaciers are twice as sharp at the scale of a region.
   - The painting, about two seconds, runs in a Web Worker (`src/workers/relief.worker.js`).
4. **Shading tied to its grid.**
   - The hillshade's exaggeration is stated for a 10′ grid and scaled with the grid's spacing (7× on the 5′ grid).
   - The darkening has a floor, so the steepest slopes stay coloured.
5. **A closer view.** The globe can be approached to 1.04 Earth radii from its centre (about 250 km above the surface), where the fine relief shows its valleys.

## 3. Alternatives Considered

- **ETOPO at 15″ through the local server.** Elevations rather than a picture of them would allow our own shading and tints. But ETOPO's server cuts subsets slowly, and the local server would become a proxy. The GIBS shading is ready, fast and in the browser's reach. Our own shading of fine elevations stays possible later.
- **A global texture of the fine shading.** At 8192 pixels it would show little more than the 5′ grid. Larger, it would take hundreds of megabytes of graphics memory. Fine shading matters only at close range, where the detail mechanism already fetches only what is in view.
- **Shading imagery styles too.** Satellite images carry their own illumination. Shading them again would double it.

## 4. Consequences

- **Detail at close range.** At close range the relief map shows 30 m relief. The first load paints a base map of 128 MB of texture.
- **Tests:** the Earth Atlas has 95 checks, including the water mask.
- **Credit.** The credits name ASTER GDEM (METI and NASA) whenever the fine relief is shown.
