# geodesy

*Status: draft. Field, height 1.* The figure of the Earth, its gravity field and rotation, and its reference systems. The other Earth sciences rest on it: geography, oceanography, meteorology and geophysics all place their objects in its frames.

## 1. Scope

| Module | Scope | Status |
|---|---|---|
| `figure` | The figure of the Earth: ellipsoids and the geoid. Geodetic and geocentric latitude on WGS 84; the Earth as a sphere of mean radius for drawing, sampling and motion. | draft (JavaScript: `geodesy.js`, `sphere.js`; C++: WGS 84, GRS 80, the IERS 2010 ellipsoid and spheres; geodetic and Earth-centred coordinates both ways, the inverse by Vermeille 2002 with a Newton step, held to ERFA; the local ENU and NED frames of a place, look angles) |
| `geodesic` | Geodesics on the ellipsoid: the direct and inverse problems, lines sampled for drawing | draft (C++: Karney 2013, ported from GeographicLib 2.1 with its MIT notice; held to GeographicLib to 4 nm over 1040 cases and its own 20 test cases) |
| `projection` | Map projections and grids: the transverse Mercator, UTM and UPS, MGRS | draft (C++: Krüger's series to n⁶ after Karney 2011; UTM's zones with the Norway and Svalbard exceptions, UPS, MGRS with its polar lettering; held to PROJ to 5 nm and to GEOTRANS's references exactly) |
| `geoid` | The geoid's undulation and heights above mean sea level | draft (C++: EGM96's 15′ grid, bilinear or cubic; bilinear held to PROJ's grid to its 3 mm quantum) |
| `directions` | Azimuths in the local horizon, the points of the compass, and their means | draft (JavaScript: `compass.js`) |
| `gravity` | The normal gravity field: the level ellipsoid of WGS 84 and GRS 80 (or any four defining constants, with f or J₂) and its derived constants (b, E, m, q₀, U₀, γₑ, γₚ, k, f*, J₂ to J₈, mean gravity); Somigliana's formula; the exact closed field at any height in ellipsoidal-harmonic coordinates, with its vector in Earth-fixed axes and the local frame; the second-order series in height | draft (C++: `gravity/normal.hpp`, target `geodesy::gravity`; every derived constant NGA and Moritz print, to their last digit; boule's field to 100 km) |
| `rotation` | The rotation of the Earth | planned |
| `reference_systems` | Coordinate reference systems, ITRF and its realisations | planned |

## 2. Dependencies

- **Domains:** metrology, for time scales, frames and Earth orientation.
- **Foundations:** physis (`geometry`: rotations and motions; `gravitation` later); calx (`quadrature`, for the mean of normal gravity; the checks' verification).

## 3. Code

`js/index.js` exports `geodesy.js`, `sphere.js` and `compass.js`. The checks (round trips on the sphere and the ellipsoid, the greatest difference of geodetic and geocentric latitude after Torge and Müller 2012, distances, directions) are run by:

```
node domains/geodesy/tests/run-tests.mjs
```

The C++ (header-only, targets `geodesy::figure`, `geodesy::geodesic`, `geodesy::projection`, `geodesy::geoid` and `geodesy::gravity`; [ADR 0049](../../development/decisions/0049-the-spacetime-of-the-earth-and-the-moon.md)) runs in the same suite. Its references are in `tests/fixtures/`, each made by a script beside it: ERFA's gd2gc and gc2gd (`make_fixtures.py`), GeographicLib's geodesics (`make_geodesic_fixtures.py`), PROJ's UTM and UPS and GEOTRANS's MGRS (`make_projection_fixtures.py`), PROJ's EGM96 (`make_geoid_fixtures.py`), and the derived constants of WGS 84 and GRS 80 as NGA.STND.0036 and Moritz print them, with boule's level ellipsoid and normal gravity (`make_gravity_fixtures.py`). The geoid's checks read the data library's grid (`python data/products/geoid_grid.py`) and are skipped without it.

## 4. Notes

- [The normal gravity of a level ellipsoid](theory/normal-gravity.md)

## 5. Open Questions

- The Earth is a sphere of mean radius (6 371 008.8 m) in `sphere.js` and the WGS 84 ellipsoid in `geodesy.js`. When should each be used, and should one module serve both, with the figure as an argument?

## 6. Used By

astronomy, geography, geophysics, oceanography, meteorology and the bridges above them, the navigation bridge's C++ among them; the shared visualization layer; both atlases; Cosmopoiesis.
