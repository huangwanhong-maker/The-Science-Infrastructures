# metrology

*Status: draft. Field, height 0.* Time scales and reference systems: every time value carries a scale, and every coordinate a frame. This domain makes both explicit and provides tested conversions between them.

Units and physical constants are general physics, so they live in physis (`units`). This domain holds only what is particular: the time scales of astronomy and of the Earth, its frames, and its orientation.

## 1. Scope

| Module | Scope | Status |
|---|---|---|
| `time` | Time scales: UTC, TAI, TT, TDB, TCG, TCB, UT1, GPS. Representations: JD, MJD, two-part JD, ISO 8601 with its scale. Leap seconds. The Earth rotation angle and sidereal time. | draft (JavaScript; C++: every scale and the conversions between them, after SOFA; leap seconds from the IANA's list or a built-in table, UTC's scaled days before 1972; TDB − TT by Fairhead and Bretagnon's 787 terms; UT1 from the IERS's series, or approximated outside it, and saying which) |
| `frames` | ICRS, GCRS, CIRS, TIRS, ITRS, ecliptic, galactic, supergalactic, topocentric and horizon frames. Precession–nutation (IAU 2006/2000A). Explicit transformations between frames. | draft (JavaScript: precession, the ecliptic, galactic and supergalactic frames, the horizon; C++: IAU 2006/2000A in full, its 1365-term nutation transcribed from ERFA, the CIO and the equinox routes, the Earth rotation angle, GMST and GAST, polar motion with s′; the Earth's frames EME2000, MOD, TOD, CIRS, TIRS, ITRS and PEF as edges of a physis frame graph, the equinox route a loop whose holonomy is the celestial pole offsets; GMST of IAU 1982, TEME's) |
| `astronomical_constants` | IAU and IERS constants, and the defining parameters of Earth ellipsoids such as WGS 84 and GRS 80, each versioned by the standard that defines it | planned |
| `earth_orientation` | Earth-orientation parameters and their interpolation; the data come through `data/` | draft (C++: the IERS's finals2000A read, its days measured or predicted, interpolated with the leap-second jump of UT1 − UTC removed) |

## 2. Dependencies

- **Domains:** none. Metrology is the lowest domain.
- **Foundations:** physis (`units`, and rotations in `geometry`); calx (numerical methods).

## 3. Code

`js/index.js` exports `time.js` and `frames.js`. The checks, against the IERS and Meeus (*Astronomical Algorithms*, examples 12.a, 13.b and 21.b), are run by:

```
node domains/metrology/tests/run-tests.mjs
```

The C++ (header-only, targets `metrology::time`, `metrology::earth_orientation` and `metrology::frames`; [ADR 0049](../../development/decisions/0049-the-spacetime-of-the-earth-and-the-moon.md)) is checked against IAU SOFA's own test values and ERFA through pyerfa (`tests/fixtures/erfa-reference.txt`, made by `tests/make_fixtures.py`): the scales and calendars, leap seconds and UTC's scaled days, the nutation, the CIO locator, the sidereal times of 2006 and 1982 and the GCRS-to-ITRS matrix, each to rounding; the Earth's orientation against rows of the IERS's own file. The same suite runs it. Its nutation and CIO series are transcribed from ERFA by `tools/transcribe_erfa.py`, which keeps ERFA's BSD notice.

## 4. Used By

astronomy, geodesy and the navigation bridge; the Celestial Atlas; Cosmopoiesis.
