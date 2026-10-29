# geophysics

*Status: draft. Field, height 2.* The solid Earth: how it moves, how it breaks, what lies beneath a place, and the magnetic field its core sustains. Proposed as a field in [ADR 0024](../../development/decisions/0024-the-solid-earth.md), and made a domain of its own by [ADR 0026](../../development/decisions/0026-domains-as-libraries.md).

Geophysics is the physics of the Earth itself. It draws on continuum mechanics, elasticity and potential theory (physis), and places its objects in geodesy's frames.

## 1. Scope

| Module | Scope | Status |
|---|---|---|
| `tectonics` | The motion of rigid plates on a sphere: angular velocities (NNR-MORVEL56; ITRF2020's plate motion model), relative motion across boundaries, after Bird's (2003) conventions | draft (JavaScript: `plates.js`) |
| `seismology` | Earthquake mechanisms from moment tensors: principal axes, nodal planes, Frohlich's classes, the radiation of P waves; aftershock rates and probabilities (Reasenberg and Jones 1989) and rupture lengths (Wells and Coppersmith 1994); later waves and the Earth's interior | draft (JavaScript: `moment.js`, `aftershocks.js`) |
| `structure` | The crust beneath a place (CRUST1.0); sections through the Earth along great circles | draft (JavaScript: `crust.js`, `section.js`) |
| `gravity` | Gravity anomalies and isostasy | planned |
| `geomagnetism` | The Earth's main magnetic field and its secular variation by the World Magnetic Model 2025 (NOAA NCEI and BGS): at a geodetic place on WGS 84 and an instant, the elements X, Y, Z, H, F, I, D and the polar grid variation with their rates; the field on the Earth-fixed axes; the centred dipole, its geomagnetic poles and McIlwain's L of the dipole; the model's span in time and height, its blackout and caution zones, and its error model; the model as one of physis's fields (`WmmField`, `geomagnetism/field.hpp`: on its base, or continued outward as the core's field beyond 850 km, noted; [ADR 0053](../../development/decisions/0053-fields-and-their-visual-expression.md)) | draft (C++: [`geomagnetism/`](geomagnetism/)) |

## 2. Dependencies

- **Domains:** geodesy, for frames (ITRF) and the figure of the Earth (WGS 84, geodetic coordinates and local frames, in the C++ too); metrology, for the calendar (the decimal years of the geomagnetic model).
- **Foundations:** physis (`continuum`, `gravitation`, `geometry`); morphe; calx (`quadrature`, `random` and `verification` in the C++).
- **Data:** WMM2025's coefficients, through the data library: `python data/acquisition/fetch.py noaa-wmm` (into `data/cache/noaa-wmm/`; public domain).

## 3. Code

`js/index.js` exports `plates.js`, `moment.js`, `crust.js`, `section.js` and `aftershocks.js`. The checks are run by:

```
node domains/geophysics/tests/run-tests.mjs
```

They compare plate velocities with the EarthScope plate calculator, GPS stations in stable plate interiors with their plates, faults with their double couples, and the nodal planes of all Global CMT events with the catalogue's; they check the crust of constructed columns and of Tibet, the Pacific and the Canadian Shield, and sections against exact spherical distances; and the aftershock counts in closed form against numerical integration. The checks of the plates' boundaries, which need geography's rasters as well, are the Earth Atlas's.

The C++ is header-only, built by CMake (`CMakeLists.txt`, the target `geophysics::geomagnetism`), and its checks run in the same suite through calx's verification. A magnetometer at a place:

```cpp
#include <geophysics/geomagnetism/wmm.hpp>
namespace gm = geophysics::geomagnetism;
constexpr double deg = std::numbers::pi / 180;

const gm::Coefficients wmm = gm::read_cof("data/cache/noaa-wmm/WMM.COF");
const double t = gm::decimal_year(2026, 10, 28, 0.5);                   // a UTC date and fraction of its day
const gm::Geodetic place{48.85 * deg, 2.35 * deg, 35.0};                // WGS 84: rad, rad, m
const gm::Field f = gm::field(wmm, place, t);                           // f.elements: T and rad; f.itrs_T: T
const gm::Validity v = gm::validity(wmm, place, t, f);                  // v.valid(), v.zone, v.reasons
```

The checks (`geomagnetism/tests/check_wmm.cpp`) hold it to the WMM2025 report's high-precision example to 10⁻¹⁰ nT, to NOAA's 100 test values and the report's Table 6 to their printed digits (but for X, where NOAA's file departs from the report's equations, and from NOAA's own C library, by up to 7.2 × 10⁻⁴ nT), to the field above the poles and the geomagnetic poles the report gives, and to the mathematics: the field as the gradient of its potential, the rates as derivatives in time, the degree-1 field as the centred dipole's. As a field (`WmmField`) it gives NOAA's points on the Earth-fixed axes, refuses outside its base in height and time with the reason, and has no divergence or curl by central differences (some $2\times10^{-7}$ of $\lvert B\rvert/r$, as the step's truncation predicts), on its base and continued outward. The Schmidt functions are calx's (`special_functions/legendre.hpp`), written here first and moved there, with their checks, since they name no object in nature. Their references are in `tests/fixtures/`, written from the data library's copy of NOAA's files by `python -I domains/geophysics/tests/make_wmm_fixtures.py` (it needs pypdf); without the coefficients in the cache the checks that need them are skipped, and say so.

## 4. Notes

- [The motion of rigid plates on a sphere](theory/plate-kinematics.md)
- [Earthquake mechanisms: moment tensors, beach balls and slabs](theory/earthquake-mechanisms.md)
- [Aftershock forecasts](theory/aftershock-forecasts.md)
- [The geomagnetic field of the World Magnetic Model](theory/geomagnetism.md)

## 5. Used By

The Earth Atlas. The geomagnetism module was written for Cosmopoiesis ([applications/cosmopoiesis-lab](../../applications/cosmopoiesis-lab/)), whose virtual magnetometers it serves.
