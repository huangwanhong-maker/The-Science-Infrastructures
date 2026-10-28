// TEME, the frame of SGP4, and its rotation to the Earth.
//
// SGP4's positions and velocities are in the "true equator, mean equinox" frame of the instant:
// the z axis the true celestial pole (the true equator of date), the x axis the mean equinox of
// date measured along it — the frame in which the Space Force's theory was written, and which no
// standard defines otherwise (Vallado, Crawford, Hujsak and Kelso 2006, "Revisiting Spacetrack
// Report #3", AIAA 2006-6753, rev. 2, §6 and Appendix C). Because its x axis is the mean equinox,
// TEME is turned to the pseudo Earth-fixed frame (PEF, the Earth's frame before polar motion) by
// the Greenwich mean sidereal time of IAU 1982 alone, not by the apparent one:
//   r_PEF = R3(θ_GMST82(UT1)) r_TEME,
//   v_PEF = R3(θ) v_TEME − ω⊕ × r_PEF,
// and the PEF to the ITRS by polar motion, r_ITRS = R1(−y_p) R2(−x_p) r_PEF (metrology's
// polar_motion, IERS Conventions 2010 eq. 5.3, with s′ neglected as Vallado neglects it; it is
// below 0.1 mas within this century). θ_GMST82 is Aoki et al.'s (1982, Astron. Astrophys. 105,
// 359) as Vallado writes it (2013, Fundamentals of Astrodynamics and Applications, eq. 3-47;
// gstime in sgp4unit.cpp):
//   θ = 67310.54841 s + (876600 h + 8640184.812866 s) T + 0.093104 s T² − 6.2 × 10⁻⁶ s T³,
// T in Julian centuries of UT1 from J2000.0.
//
// The rate of the rotation, for velocities: Vallado takes the Earth's angular velocity, ω⊕ =
// 7.29211514670698 × 10⁻⁵ rad/s × (1 − LOD/86400 s), the IERS's nominal rate corrected by the
// length of day (his teme2ecef; the default here, `TemeRate::earth`, which reproduces his published
// numbers); the derivative of θ_GMST82 itself is larger by the general precession in right
// ascension, 7.1 × 10⁻¹² rad/s (`TemeRate::gmst1982`, as skyfield takes it): 0.06 mm/s at the
// distance of a low orbit, beneath SGP4's own error by many orders.
//
// GMST 1982 is a convention of the Earth's rotation, and so metrology's (frames/iau1982.hpp).

#pragma once

#include <cmath>
#include <numbers>

#include <metrology/frames/iau1982.hpp>
#include <metrology/frames/iau2006.hpp>
#include <metrology/time/julian.hpp>
#include <physis/geometry/kinematics.hpp>
#include <physis/geometry/rotation.hpp>

namespace navigation::orbits {

using physis::geometry::Axis;
using physis::geometry::Mat3;
using physis::geometry::StateTransform;

// ω⊕ of Vallado's teme2ecef: the IERS's nominal angular velocity of the Earth, rad/s.
inline constexpr double EARTH_ANGULAR_VELOCITY = 7.29211514670698e-05;

// θ_GMST82 and its rate are metrology's (frames/iau1982.hpp).
using metrology::frames::gmst1982;
using metrology::frames::gmst1982_rate;

enum class TemeRate { earth, gmst1982 };

// TEME → PEF at an instant of UT1: x_PEF = R3(θ) x_TEME, with the rate of the frames' turning (and
// so the velocity law v_PEF = R3(θ) v_TEME − ω × r_PEF) per second of the scale the velocities are
// counted in, the UT1 rate scaled by 1 − LOD/86400 s. lod in seconds.
inline StateTransform teme_to_pef(const metrology::time::JulianDate& ut1, double lod = 0, TemeRate rate = TemeRate::earth) {
  const double theta = gmst1982(ut1);
  const double omega = (rate == TemeRate::earth ? EARTH_ANGULAR_VELOCITY : gmst1982_rate(ut1)) * (1.0 - lod / 86400.0);
  // TEME turns at −ω about z relative to the PEF (the angular velocity of the source, in target axes).
  return StateTransform::rotating(physis::geometry::frame_rotation(Axis::z, theta), {0, 0, -omega});
}

// PEF → ITRS: polar motion (x_p, y_p in radians), its own slow motion neglected.
inline StateTransform pef_to_itrs(double xp, double yp) {
  return StateTransform::fixed(metrology::frames::polar_motion(xp, yp, 0.0));
}

// TEME → ITRS at an instant of UT1, with the pole's coordinates (rad) and the length of day (s).
inline StateTransform teme_to_itrs(const metrology::time::JulianDate& ut1, double xp = 0, double yp = 0, double lod = 0,
                                   TemeRate rate = TemeRate::earth) {
  return pef_to_itrs(xp, yp) * teme_to_pef(ut1, lod, rate);
}

}  // namespace navigation::orbits
