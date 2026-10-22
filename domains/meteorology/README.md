# meteorology

*Status: draft. Field, height 2.* Meteorology in the broad sense of atmospheric science. Ocean and atmosphere are studied as one coupled system ([ADR 0007](../../development/decisions/0007-atmosphere-and-ocean-as-coupled-system.md)).

## 1. Scope

| Module | Scope | Status |
|---|---|---|
| `dynamics` | Atmospheric dynamics. Two models of our own, by the spectral transform method: the barotropic vorticity equation of the first numerical forecast (1950), and the shallow-water equations, started by a digital filter | draft (JavaScript: `barotropic.js`, `shallowwater.js`) |
| `prediction` | Weather analysis and prediction: tropical cyclones found and tracked in sea-level pressure (after TempestExtremes); the wind | draft (JavaScript: `cyclones.js`, `wind.js`) |
| `moist` | Moist thermodynamics, clouds and precipitation | planned |
| `radiation_budget` | The radiation budget of the atmosphere | planned |
| `boundary_layer` | The atmospheric boundary layer | planned |
| `chemistry` | Atmospheric chemistry | planned |

## 2. Dependencies

- **Domains:** geodesy, for reference surfaces.
- **Foundations:** physis (`continuum`, `thermodynamics`, `radiation`, `chemistry`; today the spectral transform of its JavaScript stand-in); calx.

The Earth's radius and rotation rate, as Williamson et al. (1992) give them for the standard tests, are this domain's (`EARTH_RADIUS`, `OMEGA` in `barotropic.js`): physis's spectral transform takes the radius as an argument.

## 3. Code

`js/index.js` exports `wind.js`, `barotropic.js`, `shallowwater.js` and `cyclones.js`. The checks (the Rossby–Haurwitz wave; Williamson's test case 2, also across the poles; gravity waves after Lamb 1932; the digital filter; the invariants; cyclones in constructed fields) are run by:

```
node domains/meteorology/tests/run-tests.mjs
```

## 4. Notes

- [The barotropic vorticity model](theory/barotropic-vorticity-model.md)
- [The shallow-water model](theory/shallow-water-model.md)
- [Learned forecast models](theory/learned-forecast-models.md)
- [Tropical cyclone tracking](theory/tropical-cyclone-tracking.md)

## 5. Used By

air-sea, tides and climate (planned); the Earth Atlas.
