# metrology

*Status: draft. Field, height 0.* Time scales and reference systems: every time value carries a scale, and every coordinate a frame. This domain makes both explicit and provides tested conversions between them.

Units and physical constants are general physics, so they live in physis (`units`). This domain holds only what is particular: the time scales of astronomy and of the Earth, its frames, and its orientation.

## 1. Scope

| Module | Scope | Status |
|---|---|---|
| `time` | Time scales: UTC, TAI, TT, TDB, UT1, GPS. Representations: JD, MJD, two-part JD, ISO 8601. Leap seconds. The Earth rotation angle and sidereal time. | draft (JavaScript) |
| `frames` | ICRS, GCRS, CIRS, TIRS, ITRS, ecliptic, galactic, supergalactic, topocentric and horizon frames. Precession–nutation (IAU 2006/2000A). Explicit transformations between frames. | draft (JavaScript: precession, the ecliptic, galactic and supergalactic frames, the horizon) |
| `astronomical_constants` | IAU and IERS constants, and the defining parameters of Earth ellipsoids such as WGS 84 and GRS 80, each versioned by the standard that defines it | planned |
| `earth_orientation` | Earth-orientation parameters and their interpolation; the data come through `data/` | planned |

## 2. Dependencies

- **Domains:** none. Metrology is the lowest domain.
- **Foundations:** physis (`units`, and rotations in `geometry`); calx (numerical methods).

## 3. Code

`js/index.js` exports `time.js` and `frames.js`. The checks, against the IERS and Meeus (*Astronomical Algorithms*, examples 12.a, 13.b and 21.b), are run by:

```
node domains/metrology/tests/run-tests.mjs
```

The reference for validating the C++ will be IAU SOFA, with its BSD-licensed derivative ERFA.

## 4. Used By

astronomy and geodesy; the Celestial Atlas.
