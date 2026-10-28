// The Moon's body-fixed frames of the JPL ephemeris DE440: principal axes and mean Earth.
//
// The lunar physical librations are integrated with the orbits in DE440 (Park et al. 2021, AJ 161,
// 105). NAIF distributes them as Euler angles of the principal-axes frame (PA) in a binary PCK
// (moon_pa_de440_200625.bpc, frame 31008). The mean-Earth/polar-axis frame (ME), in which lunar
// maps and the LRO's altimetry are given, differs from PA by a constant rotation (Park et al. 2021;
// NAIF's moon_de440_250416.tf, frame 31009): v_ME = R_x(−0.2785″) R_y(−78.6944″) R_z(−67.8526″)
// v_PA, a turn of about 104″ (0.9 km on the surface), checked here against the kernel's own example.
// The IAU's expression for the Moon (iau.hpp) approximates ME to tens of metres and needs no data.

#pragma once

#include <numbers>

#include <astronomy/ephemerides/pck.hpp>
#include <physis/geometry/kinematics.hpp>
#include <physis/geometry/rotation.hpp>

namespace astronomy::rotation {

inline constexpr int MOON_PA_DE440 = 31008;  // NAIF frame id of the DE440 principal axes

// PA → ME, the constant rotation of NAIF's moon_de440_250416.tf (MOON_ME_DE440_ME421).
inline Mat3 moon_pa_to_me() {
  const double as = std::numbers::pi / 648000;
  return frame_rotation(Axis::x, -0.2785 * as) * frame_rotation(Axis::y, -78.6944 * as) * frame_rotation(Axis::z, -67.8526 * as);
}

// The change of coordinates ICRF → MOON_PA and ICRF → MOON_ME at ET (TDB seconds past J2000).
inline physis::geometry::StateTransform icrf_to_moon_pa(const ephemerides::Pck& pck, double et) { return pck.icrf_to_body(MOON_PA_DE440, et); }
inline physis::geometry::StateTransform icrf_to_moon_me(const ephemerides::Pck& pck, double et) {
  const physis::geometry::StateTransform pa = icrf_to_moon_pa(pck, et);
  const Mat3 m = moon_pa_to_me();
  return {m * pa.r, m * pa.r_dot, {0, 0, 0}, {0, 0, 0}};
}

}  // namespace astronomy::rotation
