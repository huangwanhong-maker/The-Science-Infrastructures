# meteorology

*Status: draft. Field, height 2.* Meteorology in the broad sense of atmospheric science. Ocean and atmosphere are studied as one coupled system ([ADR 0007](../../development/decisions/0007-atmosphere-and-ocean-as-coupled-system.md)).

## 1. Scope

| Module | Scope | Status |
|---|---|---|
| `dynamics` | Atmospheric dynamics. Two models of our own, by the spectral transform method: the barotropic vorticity equation of the first numerical forecast (1950), and the shallow-water equations, started by a digital filter | draft (JavaScript: `barotropic.js`, `shallowwater.js`) |
| `prediction` | Weather analysis and prediction: tropical cyclones found and tracked in sea-level pressure (after TempestExtremes); the wind | draft (JavaScript: `cyclones.js`, `wind.js`) |
| `hazards` | Hazards read from forecasts: fire weather by the Canadian Fire Weather Index System; heat stress (Stull's wet-bulb temperature, the NWS heat index); daily weather from a forecast's samples (local noon, local days, amounts from rates) | draft (JavaScript: `fireweather.js`, `heat.js`) |
| `transport` | Particles in the air: the ICAO standard atmosphere, settling by Stokes' law, the paths of falling parcels in the wind (volcanic ash) | draft (JavaScript: `settling.js`) |
| `moist` | Moist thermodynamics, clouds and precipitation | planned |
| `radiation_budget` | The radiation budget of the atmosphere | planned |
| `boundary_layer` | The atmospheric boundary layer | planned |
| `chemistry` | Atmospheric chemistry | planned |

## 2. Dependencies

- **Domains:** geodesy, for reference surfaces.
- **Foundations:** physis (`continuum`, `thermodynamics`, `radiation`, `chemistry`; today the spectral transform of its JavaScript stand-in); calx.

The Earth's radius and rotation rate, as Williamson et al. (1992) give them for the standard tests, are this domain's (`EARTH_RADIUS`, `OMEGA` in `barotropic.js`): physis's spectral transform takes the radius as an argument.

## 3. Code

`js/index.js` exports `wind.js`, `barotropic.js`, `shallowwater.js`, `cyclones.js`, `fireweather.js`, `heat.js` and `settling.js`. The checks (the Rossby–Haurwitz wave; Williamson's test case 2, also across the poles; gravity waves after Lamb 1932; the digital filter; the invariants; cyclones in constructed fields; the Fire Weather Index System against Van Wagner and Pickett's standard test data; Stull's wet-bulb temperature and the NWS heat index chart; the standard atmosphere and Stokes settling) are run by:

```
node domains/meteorology/tests/run-tests.mjs
```

## 4. Notes

- [The barotropic vorticity model](theory/barotropic-vorticity-model.md)
- [The shallow-water model](theory/shallow-water-model.md)
- [Learned forecast models](theory/learned-forecast-models.md)
- [Tropical cyclone tracking](theory/tropical-cyclone-tracking.md)
- [The Canadian Forest Fire Weather Index System](theory/fire-weather-index.md)
- [Daily weather from a forecast](theory/daily-weather-from-a-forecast.md)
- [Heat stress: the wet-bulb temperature and the heat index](theory/heat-stress.md)
- [Volcanic ash carried by the wind and settling through the air](theory/settling-of-volcanic-ash.md)

## 5. Used By

air-sea, tides and climate (planned); the Earth Atlas.
