# astronomy

*Status: draft. Field, height 2.* The positions and motions of the bodies of the sky, and later their physics.

## 1. Scope

| Module | Scope | Status |
|---|---|---|
| `ephemerides` | An evaluator for the JPL DE format (SPK); the bodies of the Solar System; osculating orbits; small bodies | draft (JavaScript: `spk.js`, `bodies.js`, `kepler.js`) |
| `positional` | Apparent places, the celestial sphere, rising, transit and setting; the Sun seen from the Earth | draft (JavaScript: `sky.js`, `sun.js`) |
| `rotation` | The IAU rotation models of the Sun, the Moon and the planets (Archinal et al. 2018) | draft (JavaScript: `rotation.js`) |
| `galaxies` | The supergalactic frame, tidal indices and groups of galaxies | draft (JavaScript: `galaxies.js`) |
| `astrometry` | Astrometric reduction of images | planned |
| later | Astrophysics and cosmology | planned |

## 2. Dependencies

- **Domains:** metrology (time scales and frames); geodesy (the observer's place).
- **Foundations:** physis (mechanics, gravitation and geometry); calx (numerical methods).

## 3. Code

`js/index.js` exports `bodies.js`, `spk.js`, `kepler.js`, `rotation.js`, `sky.js`, `sun.js` and `galaxies.js`. Reading the ephemeris kernel is the domain's part; fetching it is the application's.

The checks, against JPL Horizons, Meeus, Tully (2015) and Karachentsev et al. (2013), are run by:

```
node domains/astronomy/tests/run-tests.mjs
```

The ephemeris checks need DE440s (`python data/acquisition/fetch.py jpl-de440s`), and those of the galaxies their catalogues (`tully-2mass-groups`, `ungc`). The reference values are in `tests/fixtures/`, fetched again by `fetch_horizons.py`, `fetch_orientation.py` and `fetch_subsolar.py`.

## 4. Notes

- [Relations between galaxies: supergalactic frame, tidal index, Hubble flow](theory/galaxy-relations.md)
- [The subsolar point, the terminator and twilight](theory/subsolar-point-and-twilight.md)

## 5. Open Questions

- The Sun is computed twice: from DE440 in `sky.js`, and by the Astronomical Almanac's low-precision theory in `sun.js`, which needs no ephemeris file. `sun.js` also keeps its own sidereal time, where `sky.js` uses metrology's. One position of the Sun, at a declared accuracy, should serve both.

## 6. Used By

climate (orbital forcing, planned), tides and navigation (planned); the Celestial Atlas and the Earth Atlas.
