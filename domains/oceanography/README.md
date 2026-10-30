# oceanography

*Status: draft. Field, height 2.* The ocean: its circulation, waves and water, its chemistry and life, and its floor.

## 1. Scope

| Module | Scope | Status |
|---|---|---|
| `physical` | Circulation and waves: the ocean column beneath a place (inertial period, mixed layer, sampling in time and depth); the propagation of wave energy; tsunami travel times, by shortest paths over the bathymetry | draft (JavaScript: `ocean.js`, `waves.js`, `tsunami.js`) |
| `water` | The upper ocean's water as a body in it meets it: its density presumed from the range of surface seawater (1020–1030 kg/m³) and its near-surface current assumed, each realised for a world's seed, as fields (`presumed.hpp`, ADR 0059) | draft (C++) |
| `seawater` | Seawater thermodynamics (TEOS-10), from which `water`'s density will follow where the water's temperature and salinity are known | planned |
| `biogeochemistry` | Chemical and biological oceanography | planned |
| `seafloor` | Bathymetry and the sea floor | planned |

## 2. Dependencies

- **Domains:** geodesy, for the geoid and other reference surfaces.
- **Foundations:** physis (`continuum`, `thermodynamics`, `chemistry`; `analysis` and `fields` for `water`); calx (numerical methods and device kernels).

## 3. Code

`js/index.js` exports `ocean.js`, `waves.js` and `tsunami.js`. The checks (the inertial period, the mixed layer by de Boyer Montégut et al.'s criterion, sampling in time and depth, deep-water group speed after Holthuijsen 2007, tsunami travel times over constructed oceans: great-circle times, walls, closed seas) are run by the command below, which also builds the C++ (`water`, its presumption's interval, its realisation's statistics, its current horizontal, its fields) with CMake and runs it on calx's verification:

```
node domains/oceanography/tests/run-tests.mjs
```

## 4. Notes

- [The ocean column](theory/the-ocean-column.md)
- [Wave energy propagation](theory/wave-energy-propagation.md)
- [Tsunami travel times](theory/tsunami-travel-times.md)

## 5. Used By

air-sea, tides and climate (planned); the Earth Atlas.
