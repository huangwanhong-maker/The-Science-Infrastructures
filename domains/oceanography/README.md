# oceanography

*Status: draft. Field, height 2.* The ocean: its circulation, waves and water, its chemistry and life, and its floor.

## 1. Scope

| Module | Scope | Status |
|---|---|---|
| `physical` | Circulation and waves: the ocean column beneath a place (inertial period, mixed layer, sampling in time and depth); the propagation of wave energy | draft (JavaScript: `ocean.js`, `waves.js`) |
| `seawater` | Seawater thermodynamics (TEOS-10) | planned |
| `biogeochemistry` | Chemical and biological oceanography | planned |
| `seafloor` | Bathymetry and the sea floor | planned |

## 2. Dependencies

- **Domains:** geodesy, for the geoid and other reference surfaces.
- **Foundations:** physis (`continuum`, `thermodynamics`, `chemistry`); calx (numerical methods and device kernels).

## 3. Code

`js/index.js` exports `ocean.js` and `waves.js`. The checks (the inertial period, the mixed layer by de Boyer Montégut et al.'s criterion, sampling in time and depth, deep-water group speed after Holthuijsen 2007) are run by:

```
node domains/oceanography/tests/run-tests.mjs
```

## 4. Notes

- [The ocean column](theory/the-ocean-column.md)
- [Wave energy propagation](theory/wave-energy-propagation.md)

## 5. Used By

air-sea, tides and climate (planned); the Earth Atlas.
