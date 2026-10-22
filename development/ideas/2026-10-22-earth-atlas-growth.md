---
title: Directions for the Earth Atlas
status: seed
date: 2026-10-22
---

# Directions for the Earth Atlas

Written while building the first version ([ADR 0014](../decisions/0014-earth-atlas-separated-data-and-typed-layers.md)). None of this is decided.

1. **The climate itself.** Register CMIP6 projections (ESGF) and ERA5 (Copernicus, with an account) as field products. Two things are then needed:
   - a timeline in years, beside the timeline in hours;
   - the comparison of two fields as a difference: a scenario against the historical period, or a forecast against an analysis.
2. **True trajectories.** Integrate air or water parcels through the forecast's own time, rather than along the instantaneous field. Each trajectory starts at a chosen place; this is the Lagrangian view of a flow.
3. **Ocean currents** (OSCAR, or CMEMS with an account) as a second vector field. They test that a new flow appears without changes to code.
4. **Questions across fields**, the reason for one atlas:
   - Does the SST anomaly explain a precipitation anomaly?
   - Where do people live relative to forecast heat?

   This needs gridded population (GHSL or WorldPop) beside country totals. It also needs zonal statistics, a field summarised over each country, which is again a join, this time from a grid to polygons.
5. **Our own models.** The output of the laboratory's own simulations (calx, physis) can be written as netCDF and shown with `field_grid.py --input`. A shallow-water model on the sphere would be a first test. (A shallow-water model now runs in the atlas itself: [ADR 0019](../decisions/0019-shallow-water-model.md). The next step is [a primitive-equation model](2026-10-22-primitive-equation-model.md).)
6. **Shared web modules.** The Celestial Atlas and the Earth Atlas share a server, label placement, palette and time arithmetic. These could move into `visualization/` (renderers and styles) once the user agrees.
7. **Higher resolution.** Tiles or chunked grids (Zarr) for kilometre-scale data such as Sentinel imagery and regional models.
