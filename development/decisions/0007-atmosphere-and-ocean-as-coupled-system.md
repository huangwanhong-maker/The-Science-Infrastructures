# 0007. Atmosphere and Ocean as One Coupled System

- **Status:** Accepted; amended by [0008](0008-independent-libraries.md): the shared physics now lives in physis; amended by [0010](0010-domain-organised-fields-and-applications.md): the field is named meteorology
- **Date:** 2026-10-21
- **Decided by:** the user ("also add atmospheric and meteorology … ocean cannot leave them"). Claude designed the placement.

## 1. Context

The initial scope named astronomy, geometry and oceanography. Climate appeared only as a possible adjacent field. Yet the ocean cannot be understood without the atmosphere:

- winds drive the upper-ocean circulation and the waves;
- heat, freshwater and gases pass through the sea surface;
- phenomena such as ENSO and the monsoons exist only in the coupled system.

Both fluids obey the same equations of rotating, stratified flow. Both are driven by solar radiation, and astronomy sets how that radiation is distributed over the Earth.

## 2. Decision

1. **Scope.** Atmospheric science, including meteorology, is in scope on the same footing as oceanography.
2. **Shared physics.** The parts of physics that ocean and atmosphere share are written once, in `core/physics/`:
   - geophysical fluid dynamics lives in `continuum`;
   - a new `thermodynamics` module provides the common base for seawater and moist air;
   - `radiation` changes from "later" to planned.
3. **New field and bridges.** `core/domains/` gains one field and two bridges:
   - `atmosphere` (field): atmospheric dynamics, moist thermodynamics, clouds and precipitation, the radiation budget, the boundary layer, atmospheric chemistry, and weather analysis and prediction.
   - `air_sea` (bridge: atmosphere, oceanography): fluxes of momentum, heat, freshwater and gases across the sea surface, and wind waves.
   - `climate` (bridge: atmosphere, oceanography, astronomy): the coupled system, energy balance, variability such as ENSO and the monsoons, and orbital forcing (Milankovitch cycles).
4. **Tides.** The `tides` bridge extends to atmospheric tides.

## 3. Alternatives Considered

- **A single merged `geofluids` domain for ocean and atmosphere.** Rejected. The two differ in thermodynamics (seawater versus moist air), in their observations and in their literature. What they share belongs in `core/physics/` and in the bridges, not in one module.
- **Atmosphere as a later, adjacent field.** Rejected by the user: the ocean cannot be studied without the atmosphere.

## 4. Consequences

- The scope statement, the README and the layout are updated.
- New datasets and validation problems are planned. Examples are the ERA5 reanalysis, and the Held–Suarez and Jablonowski–Williamson dynamical-core tests.
- Coupled modelling enters the roadmap. The coupling interface between the ocean and atmosphere components is where the two fields meet, so it must be designed carefully.
