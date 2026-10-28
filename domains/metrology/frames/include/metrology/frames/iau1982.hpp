// The Greenwich mean sidereal time of IAU 1982, the sidereal time of the frames older than the
// IAU 2006/2000A system: of SGP4's TEME above all (the navigation bridge, orbits/teme.hpp), whose x
// axis is the mean equinox and which is turned to the Earth by this angle alone.
//
//   θ_GMST82 = 67310.54841 s + (876600 h + 8640184.812866 s) T + 0.093104 s T² − 6.2 × 10⁻⁶ s T³,
//
// T in Julian centuries of UT1 from J2000.0 (Aoki et al. 1982, Astron. Astrophys. 105, 359–361;
// SOFA's iauGmst82; Vallado 2013, Fundamentals of Astrodynamics and Applications, eq. 3-47). The
// checks hold it to ERFA's gmst82. It departs from the IAU 2006 GMST (iau2006.hpp, `gmst`) by
// −15 mas at J2000.0 and 60 mas in 2026, growing by about 2.8 mas a year (ERFA's gmst82 against
// gmst06): the change of the precession model and the definition of the sidereal time by the Earth
// rotation angle.
//
// Two forms: of a two-part Julian date, keeping the parts' precision (the whole turns of the term
// 876600 h T, one turn a day, taken from the date's fractions of a day; as SOFA and skyfield); and
// of a date in one double, exactly as Vallado's gstime computes it, which SGP4 takes for the sidereal
// time of its epoch so that its answers are Vallado's to the last digit.

#pragma once

#include <cmath>
#include <numbers>

#include <metrology/time/julian.hpp>

namespace metrology::frames {

// θ_GMST82 (rad, in [0, 2π)) of a Julian date of UT1 in one double (Vallado's gstime).
inline double gmst1982(double jd_ut1) {
  constexpr double twopi = 2.0 * std::numbers::pi, deg2rad = std::numbers::pi / 180.0;
  const double tut1 = (jd_ut1 - 2451545.0) / 36525.0;
  double temp = -6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 + (876600.0 * 3600 + 8640184.812866) * tut1 + 67310.54841;  // s
  temp = std::fmod(temp * deg2rad / 240.0, twopi);  // 360°/86400 s = 1°/240 s, then radians
  if (temp < 0.0) temp += twopi;
  return temp;
}

// θ_GMST82 (rad, in [0, 2π)) of a two-part Julian date of UT1, to the parts' precision.
inline double gmst1982(const metrology::time::JulianDate& ut1) {
  constexpr double twopi = 2.0 * std::numbers::pi;
  const double t = ut1.since_j2000() / 36525.0;
  const double g = 67310.54841 + (8640184.812866 + (0.093104 + (-6.2e-6) * t) * t) * t;  // s
  double turns = std::fmod(ut1.d1, 1.0) + std::fmod(ut1.d2, 1.0) + std::fmod(g / 86400.0, 1.0);
  turns -= std::floor(turns);
  return turns * twopi;
}

// dθ_GMST82/dt, rad per second of UT1.
inline double gmst1982_rate(const metrology::time::JulianDate& ut1) {
  const double t = ut1.since_j2000() / 36525.0;
  const double dg = 8640184.812866 + (0.093104 * 2.0 + (-6.2e-6 * 3.0) * t) * t;  // s per century
  return 2.0 * std::numbers::pi / 86400.0 * (1.0 + dg / (86400.0 * 36525.0));
}

}  // namespace metrology::frames
