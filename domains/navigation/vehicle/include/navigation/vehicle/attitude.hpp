// The frames of a vehicle and what navigation reads from its state: attitude by the aerospace Euler
// angles, the body's rates, ground speed, track and flight path, and the air-relative angles.
//
// Conventions (ISO 1151-1:1988; ANSI/AIAA R-004-1992, "Atmospheric and Space Flight Vehicle
// Coordinate Systems"; Stevens, Lewis and Johnson 2016, Aircraft Control and Simulation, 3rd ed.,
// §1.3–1.4):
//   - the body frame has x forward along the vehicle's reference axis, y to the right (starboard
//     wing), z down, completing a right-handed triad;
//   - the local frame is north-east-down (NED) at the vehicle's geodetic position (geodesy);
//   - the attitude is the yaw ψ (heading, from north toward east), pitch θ (nose up positive) and
//     roll φ (right wing down positive) of the 3-2-1 sequence: the body's axes are the NED axes
//     turned by ψ about z, then θ about the new y, then φ about the new x, so that coordinates
//     transform as x_body = R₁(φ) R₂(θ) R₃(ψ) x_NED, with Rᵢ the passive (frame) rotations;
//   - the body rates p, q, r about x, y, z relate to the Euler rates by
//       φ̇ = p + (q sin φ + r cos φ) tan θ,  θ̇ = q cos φ − r sin φ,  ψ̇ = (q sin φ + r cos φ)/cos θ
//     (Stevens et al., eq. 1.4-4), singular at θ = ±90° (gimbal lock), where a quaternion serves;
//   - from a velocity in NED: the ground speed √(v_N² + v_E²), the track (course over ground)
//     atan2(v_E, v_N), the flight path angle γ = atan2(−v_D, ground speed) and the vertical speed −v_D;
//   - from the air-relative velocity in body axes (u, v, w): the airspeed V, the angle of attack
//     α = atan2(w, u) and the sideslip β = asin(v/V) (Stevens et al., eq. 2.3-6).

#pragma once

#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>

#include <geodesy/figure/ellipsoid.hpp>
#include <geodesy/figure/local.hpp>
#include <physis/geometry/kinematics.hpp>
#include <physis/geometry/rotation.hpp>

namespace navigation::vehicle {

using physis::geometry::Mat3;
using physis::geometry::StateTransform;
using physis::geometry::Vec3;
using namespace physis::geometry::operators;

struct Euler {
  double yaw = 0, pitch = 0, roll = 0;  // rad: ψ, θ, φ of the 3-2-1 sequence
};

// The matrix taking NED coordinates to body coordinates: R₁(φ) R₂(θ) R₃(ψ).
inline Mat3 ned_to_body(const Euler& a) {
  using physis::geometry::Axis;
  using physis::geometry::frame_rotation;
  return frame_rotation(Axis::x, a.roll) * frame_rotation(Axis::y, a.pitch) * frame_rotation(Axis::z, a.yaw);
}

// The Euler angles of a matrix taking NED to body coordinates (ψ in (−π, π], θ in [−π/2, π/2]). At
// gimbal lock (|θ| = π/2) the roll is taken as zero and the yaw carries the whole turn about z.
inline Euler euler_of(const Mat3& c) {
  Euler a;
  const double s = -c[0][2];
  if (std::fabs(s) >= 1 - 1e-15) {
    a.pitch = std::copysign(std::numbers::pi / 2, s);
    a.roll = 0;
    a.yaw = std::atan2(-c[1][0], c[1][1]);
    return a;
  }
  a.pitch = std::asin(s);
  a.roll = std::atan2(c[1][2], c[2][2]);
  a.yaw = std::atan2(c[0][1], c[0][0]);
  return a;
}

// The Euler rates of body rates (p, q, r) at an attitude.
inline Euler euler_rates(const Euler& a, const Vec3& pqr) {
  const double sp = std::sin(a.roll), cp = std::cos(a.roll), ct = std::cos(a.pitch), tt = std::tan(a.pitch);
  if (std::fabs(ct) < 1e-12) throw std::domain_error("navigation: Euler rates are undefined at gimbal lock (pitch ±90°)");
  Euler d;
  d.roll = pqr[0] + (pqr[1] * sp + pqr[2] * cp) * tt;
  d.pitch = pqr[1] * cp - pqr[2] * sp;
  d.yaw = (pqr[1] * sp + pqr[2] * cp) / ct;
  return d;
}

// The body rates of Euler rates at an attitude (the inverse relation).
inline Vec3 body_rates(const Euler& a, const Euler& d) {
  const double sp = std::sin(a.roll), cp = std::cos(a.roll), st = std::sin(a.pitch), ct = std::cos(a.pitch);
  return {d.roll - d.yaw * st, d.pitch * cp + d.yaw * ct * sp, -d.pitch * sp + d.yaw * ct * cp};
}

// What navigation reads from a velocity in NED.
struct Motion {
  double ground_speed = 0;   // m/s
  double track = 0;          // rad, course over ground from north toward east, in [0, 2π)
  double flight_path = 0;    // rad, the climb angle γ
  double vertical_speed = 0; // m/s, up positive
};
inline Motion motion_of(const Vec3& v_ned) {
  Motion m;
  m.ground_speed = std::hypot(v_ned[0], v_ned[1]);
  m.track = std::atan2(v_ned[1], v_ned[0]);
  if (m.track < 0) m.track += 2 * std::numbers::pi;
  m.flight_path = std::atan2(-v_ned[2], m.ground_speed);
  m.vertical_speed = -v_ned[2];
  return m;
}

// The air-relative angles of the air-relative velocity in body axes.
struct AirData {
  double airspeed = 0;  // m/s
  double alpha = 0;     // rad, the angle of attack
  double beta = 0;      // rad, the sideslip
};
inline AirData air_data(const Vec3& v_body) {
  AirData d;
  d.airspeed = std::sqrt(v_body[0] * v_body[0] + v_body[1] * v_body[1] + v_body[2] * v_body[2]);
  if (!(d.airspeed > 0)) return d;
  d.alpha = std::atan2(v_body[2], v_body[0]);
  d.beta = std::asin(std::max(-1.0, std::min(1.0, v_body[1] / d.airspeed)));
  return d;
}

// A vehicle at a geodetic place with an attitude: the motion from its body frame to the Earth's frame
// (ITRS for WGS 84), x_earth = R x_body + b, with the body's rotation rate (p, q, r) and velocity in NED.
inline StateTransform body_to_earth(const geodesy::figure::Ellipsoid& e, const geodesy::figure::Geodetic& where, const Euler& attitude,
                                    const Vec3& pqr = {0, 0, 0}, const Vec3& v_ned = {0, 0, 0}) {
  const auto place = geodesy::figure::Place::at(e, where);
  const StateTransform ned = place.ned_to_earth();  // NED coordinates to Earth coordinates
  const Mat3 cbn = ned_to_body(attitude);           // NED to body
  StateTransform t;
  t.r = ned.r * physis::geometry::transpose(cbn);   // body to Earth
  // The body turns at ω_b = (p, q, r) in body axes: Ṙ = R [ω_b]× (the place's frame held still).
  t.r_dot = t.r * physis::geometry::skew(pqr);
  t.b = ned.b;
  t.b_dot = ned.r * v_ned;
  return t;
}

}  // namespace navigation::vehicle
