# metrology: Time Scales and Reference Systems

*Status: planned.* This is layer 1 of `core/`, and the first layer to be built in Phase 1.

## 1. Purpose

Every time value in the laboratory carries a time scale, and every coordinate carries a reference frame. This layer makes both explicit and provides tested conversions between them.

Units and physical constants are general physics, so they live in physis (module `units`). This layer holds only what is specific to astronomy and the Earth.

## 2. Planned Modules

| Module | Scope |
|---|---|
| `time` | Time scales: UTC, TAI, TT, TDB, UT1, GPS. Representations: JD, MJD, two-part JD, ISO 8601. Leap seconds. |
| `frames` | Frames: ICRS, GCRS, CIRS, TIRS, ITRS, ecliptic, topocentric and horizon frames. Precession–nutation (IAU 2006/2000A). Explicit transformations between frames. |
| `astronomical_constants` | IAU and IERS constants, and the defining parameters of Earth ellipsoids such as WGS 84 and GRS 80. Each is versioned by the standard that defines it. |
| `earth_orientation` | Earth-orientation parameters and their interpolation. The data are supplied by `data/`. |

## 3. Dependencies

- physis: units and rotations.
- calx: numerical methods.

## 4. Reference

The reference for validation is IAU SOFA, together with its BSD-licensed derivative ERFA.
