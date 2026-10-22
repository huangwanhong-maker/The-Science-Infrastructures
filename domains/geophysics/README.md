# geophysics

*Status: draft. Field, height 2.* The solid Earth: how it moves, how it breaks, and what lies beneath a place. Proposed as a field in [ADR 0024](../../development/decisions/0024-the-solid-earth.md), and made a domain of its own by [ADR 0026](../../development/decisions/0026-domains-as-libraries.md).

Geophysics is the physics of the Earth itself. It draws on continuum mechanics, elasticity and potential theory (physis), and places its objects in geodesy's frames.

## 1. Scope

| Module | Scope | Status |
|---|---|---|
| `tectonics` | The motion of rigid plates on a sphere: angular velocities (NNR-MORVEL56; ITRF2020's plate motion model), relative motion across boundaries, after Bird's (2003) conventions | draft (JavaScript: `plates.js`) |
| `seismology` | Earthquake mechanisms from moment tensors: principal axes, nodal planes, Frohlich's classes, the radiation of P waves; aftershock rates and probabilities (Reasenberg and Jones 1989) and rupture lengths (Wells and Coppersmith 1994); later waves and the Earth's interior | draft (JavaScript: `moment.js`, `aftershocks.js`) |
| `structure` | The crust beneath a place (CRUST1.0); sections through the Earth along great circles | draft (JavaScript: `crust.js`, `section.js`) |
| `gravity` | Gravity anomalies and isostasy | planned |
| `magnetism` | The geomagnetic field and its variation | planned |

## 2. Dependencies

- **Domains:** geodesy, for frames (ITRF) and the figure of the Earth.
- **Foundations:** physis (`continuum`, `gravitation`, `geometry`); morphe; calx.

## 3. Code

`js/index.js` exports `plates.js`, `moment.js`, `crust.js`, `section.js` and `aftershocks.js`. The checks are run by:

```
node domains/geophysics/tests/run-tests.mjs
```

They compare plate velocities with the EarthScope plate calculator, GPS stations in stable plate interiors with their plates, faults with their double couples, and the nodal planes of all Global CMT events with the catalogue's; they check the crust of constructed columns and of Tibet, the Pacific and the Canadian Shield, and sections against exact spherical distances; and the aftershock counts in closed form against numerical integration. The checks of the plates' boundaries, which need geography's rasters as well, are the Earth Atlas's.

## 4. Notes

- [The motion of rigid plates on a sphere](theory/plate-kinematics.md)
- [Earthquake mechanisms: moment tensors, beach balls and slabs](theory/earthquake-mechanisms.md)
- [Aftershock forecasts](theory/aftershock-forecasts.md)

## 5. Used By

The Earth Atlas.
