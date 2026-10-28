// Local level frames at a place: east-north-up, north-east-down, and azimuth, elevation and range.
//
// At a point of geodetic latitude φ and longitude λ, the local level frame has its z axis along
// the ellipsoid's normal and its x, y axes in the tangent plane. ENU (east, north, up) is the
// convention of geodesy and of surveying; NED (north, east, down) that of aircraft and ships,
// whose attitude angles are measured from it (physis rotation.hpp, Euler sequence Z-Y-X). The
// rotation from Earth-centred (ITRS) coordinates to ENU is
//   R = [[−sin λ, cos λ, 0], [−sin φ cos λ, −sin φ sin λ, cos φ], [cos φ cos λ, cos φ sin λ, sin φ]]
// (Hofmann-Wellenhof, Lichtenegger and Wasle 2008, GNSS, §10.2.2.3), and NED's rows are ENU's
// north, east and −up. Azimuth is measured from north through east, elevation above the tangent
// plane (CLAUDE.md §4.1: the origin of azimuth is north, clockwise seen from above).
//
// The frame of a place is fixed to the Earth: its edge in a frame graph (physis frame_graph.hpp)
// is a constant rigid motion from the place's frame to the Earth-centred one.

#pragma once

#include <cmath>
#include <numbers>

#include <geodesy/figure/ellipsoid.hpp>
#include <physis/geometry/kinematics.hpp>

namespace geodesy::figure {

using physis::geometry::Mat3;
using physis::geometry::StateTransform;
using physis::geometry::transpose;
using namespace physis::geometry::operators;

// Earth-centred coordinates to east-north-up at a geodetic latitude and longitude.
inline Mat3 ecef_to_enu(double phi, double lambda) {
  const double sp = std::sin(phi), cp = std::cos(phi), sl = std::sin(lambda), cl = std::cos(lambda);
  return {{{-sl, cl, 0}, {-sp * cl, -sp * sl, cp}, {cp * cl, cp * sl, sp}}};
}
inline Mat3 ecef_to_ned(double phi, double lambda) {
  const Mat3 e = ecef_to_enu(phi, lambda);
  return {{e[1], e[0], {-e[2][0], -e[2][1], -e[2][2]}}};
}
inline Mat3 enu_to_ned() { return {{{0, 1, 0}, {1, 0, 0}, {0, 0, -1}}}; }

// Azimuth (from north through east, [0, 2π)), elevation and range of a vector given in ENU.
struct LookAngles {
  double azimuth = 0, elevation = 0, range = 0;
};
inline LookAngles look_angles(const Vec3& enu) {
  LookAngles l;
  l.range = physis::geometry::norm(enu);
  double az = std::atan2(enu[0], enu[1]);
  if (az < 0) az += 2 * std::numbers::pi;
  l.azimuth = az;
  l.elevation = std::atan2(enu[2], std::hypot(enu[0], enu[1]));
  return l;
}
inline Vec3 enu_from_look(const LookAngles& l) {
  const double ce = std::cos(l.elevation);
  return {l.range * ce * std::sin(l.azimuth), l.range * ce * std::cos(l.azimuth), l.range * std::sin(l.elevation)};
}

// A place on an ellipsoid, with its local level frame.
struct Place {
  Geodetic where;
  Vec3 origin;  // Earth-centred
  Mat3 to_enu;  // Earth-centred directions → ENU

  static Place at(const Ellipsoid& e, const Geodetic& g) { return {g, e.cartesian(g), ecef_to_enu(g.latitude, g.longitude)}; }
  Vec3 enu_of(const Vec3& ecef_point) const { return to_enu * (ecef_point - origin); }
  Vec3 ecef_of(const Vec3& enu) const { return transpose(to_enu) * enu + origin; }
  LookAngles look_at(const Vec3& ecef_point) const { return look_angles(enu_of(ecef_point)); }
  // The constant motion from this place's ENU frame to the Earth-centred frame, as an edge of a
  // frame graph: x_ecef = Rᵀ x_enu + origin.
  StateTransform enu_to_earth() const { return StateTransform::fixed(transpose(to_enu), origin); }
  StateTransform ned_to_earth() const { return StateTransform::fixed(transpose(enu_to_ned() * to_enu), origin); }
};

}  // namespace geodesy::figure
