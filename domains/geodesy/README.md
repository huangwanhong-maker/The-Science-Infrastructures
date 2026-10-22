# geodesy

*Status: draft. Field, height 1.* The figure of the Earth, its gravity field and rotation, and its reference systems. The other Earth sciences rest on it: geography, oceanography, meteorology and geophysics all place their objects in its frames.

## 1. Scope

| Module | Scope | Status |
|---|---|---|
| `figure` | The figure of the Earth: ellipsoids and the geoid. Geodetic and geocentric latitude on WGS 84; the Earth as a sphere of mean radius for drawing, sampling and motion. | draft (JavaScript: `geodesy.js`, `sphere.js`) |
| `directions` | Azimuths in the local horizon, the points of the compass, and their means | draft (JavaScript: `compass.js`) |
| `gravity` | The gravity field | planned |
| `rotation` | The rotation of the Earth | planned |
| `reference_systems` | Coordinate reference systems, ITRF and its realisations | planned |

## 2. Dependencies

- **Domains:** metrology, for time scales, frames and Earth orientation.
- **Foundations:** physis (`gravitation`, `geometry`); calx.

## 3. Code

`js/index.js` exports `geodesy.js`, `sphere.js` and `compass.js`. The checks (round trips on the sphere and the ellipsoid, the greatest difference of geodetic and geocentric latitude after Torge and Müller 2012, distances, directions) are run by:

```
node domains/geodesy/tests/run-tests.mjs
```

## 4. Open Questions

- The Earth is a sphere of mean radius (6 371 008.8 m) in `sphere.js` and the WGS 84 ellipsoid in `geodesy.js`. When should each be used, and should one module serve both, with the figure as an argument?

## 5. Used By

astronomy, geography, geophysics, oceanography, meteorology and the bridges above them; the shared visualization layer; both atlases.
