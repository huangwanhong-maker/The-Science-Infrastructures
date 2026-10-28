// Checks of a vehicle's frames (navigation/vehicle): the aerospace Euler angles, body and Euler rates,
// ground speed, track and flight path, air data, and the body frame on the Earth.
//
// References: closed forms of the 3-2-1 sequence (Stevens, Lewis and Johnson 2016, §1.3–1.4); physis's
// own Euler sequences (an intrinsic Z-Y-X rotation is the body's attitude); the time derivative of the
// attitude matrix, Ċ = −[ω]× C, against central differences; geodesy's local frames.

#include <cmath>
#include <numbers>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <navigation/vehicle/attitude.hpp>

using namespace navigation::vehicle;
using calx::verification::Suite;
namespace pg = physis::geometry;

namespace {
constexpr double pi = std::numbers::pi, deg = pi / 180;
double max_abs(const Vec3& v) { return std::max({std::fabs(v[0]), std::fabs(v[1]), std::fabs(v[2])}); }
double gap(double a, double b) { return std::fabs(std::remainder(a - b, 2 * pi)); }
}  // namespace

int main() {
  Suite s;
  calx::random::Stream rng(20261028, 2, 0, 0);

  // ---- the Euler angles --------------------------------------------------------------------------------
  {
    double round = 0, against_physis = 0;
    for (int k = 0; k < 500; k++) {
      const Euler a{pi * (2 * rng.uniform() - 1), 0.49 * pi * (2 * rng.uniform() - 1), pi * (2 * rng.uniform() - 1)};
      const Euler b = euler_of(ned_to_body(a));
      round = std::max({round, gap(a.yaw, b.yaw), std::fabs(a.pitch - b.pitch), gap(a.roll, b.roll)});
      // The body's attitude as an active rotation of the NED frame: intrinsic z (yaw), y (pitch), x (roll).
      const pg::Rotation r = pg::Rotation::from_euler(pg::EulerSequence::parse("ZYX", pg::EulerConvention::intrinsic), a.yaw, a.pitch, a.roll);
      against_physis = std::max(against_physis, pg::max_difference(r.matrix(), pg::transpose(ned_to_body(a))));
    }
    s.check("Euler angles → attitude matrix → Euler angles, 500 attitudes", round, 1e-12, "rad");
    s.check("the attitude is physis's intrinsic Z-Y-X rotation of the NED axes (to rounding)", against_physis, 4e-15);
    // Heading east, level: the nose points east, the right wing south.
    const Mat3 c = ned_to_body({90 * deg, 0, 0});
    s.check("heading 90°: the nose points east", max_abs(pg::transpose(c) * Vec3{1, 0, 0} - Vec3{0, 1, 0}), 1e-15);
    s.check("heading 90°: the right wing points south", max_abs(pg::transpose(c) * Vec3{0, 1, 0} - Vec3{-1, 0, 0}), 1e-15);
    // Nose up 10°: the nose rises above the horizon (NED's z is down).
    s.check("pitch 10°: the nose rises", (pg::transpose(ned_to_body({0, 10 * deg, 0})) * Vec3{1, 0, 0})[2] + std::sin(10 * deg), 1e-15);
    // Gimbal lock: the yaw carries the turn about z, the roll is zero, and the matrix is kept.
    const Euler locked = euler_of(ned_to_body({30 * deg, 90 * deg, 20 * deg}));
    s.check("at gimbal lock the attitude is kept", pg::max_difference(ned_to_body(locked), ned_to_body({30 * deg, 90 * deg, 20 * deg})), 1e-12);
  }

  // ---- rates ------------------------------------------------------------------------------------------------
  {
    double round = 0, derivative = 0;
    for (int k = 0; k < 200; k++) {
      const Euler a{pi * (2 * rng.uniform() - 1), 0.45 * pi * (2 * rng.uniform() - 1), pi * (2 * rng.uniform() - 1)};
      const Vec3 pqr{0.3 * (2 * rng.uniform() - 1), 0.3 * (2 * rng.uniform() - 1), 0.3 * (2 * rng.uniform() - 1)};
      const Euler d = euler_rates(a, pqr);
      round = std::max(round, max_abs(body_rates(a, d) - pqr));
      // Ċ = −[ω]× C, with C the attitude matrix and ω the body rates in body axes.
      const double h = 1e-5;
      const Euler ap{a.yaw + h * d.yaw, a.pitch + h * d.pitch, a.roll + h * d.roll}, am{a.yaw - h * d.yaw, a.pitch - h * d.pitch, a.roll - h * d.roll};
      const Mat3 cdot = (1 / (2 * h)) * (ned_to_body(ap) - ned_to_body(am));
      derivative = std::max(derivative, pg::max_difference(cdot, -1.0 * (pg::skew(pqr) * ned_to_body(a))));
    }
    s.check("body rates → Euler rates → body rates", round, 1e-14, "rad/s");
    s.check("the attitude matrix turns at the body rates: Ċ = −[ω]× C", derivative, 1e-9, "1/s");
    bool refused = false;
    try { euler_rates({0, 90 * deg, 0}, {0.1, 0, 0}); } catch (const std::domain_error&) { refused = true; }
    s.require("Euler rates are refused at gimbal lock", refused);
  }

  // ---- motion and air data ------------------------------------------------------------------------------------
  {
    // Climbing at 3° along a track of 45° at 100 m/s.
    const double v = 100, g = 3 * deg, t = 45 * deg;
    const Motion m = motion_of({v * std::cos(g) * std::cos(t), v * std::cos(g) * std::sin(t), -v * std::sin(g)});
    s.check("a 3° climb on track 045 at 100 m/s: ground speed", m.ground_speed - v * std::cos(g), 1e-12, "m/s");
    s.check("a 3° climb on track 045: the track", m.track - t, 1e-15, "rad");
    s.check("a 3° climb: the flight path angle", m.flight_path - g, 1e-15, "rad");
    s.check("a 3° climb: the vertical speed", m.vertical_speed - v * std::sin(g), 1e-12, "m/s");
    s.check("a westward track is 270°", motion_of({0, -5, 0}).track - 1.5 * pi, 1e-15, "rad");
    const AirData d = air_data({60 * std::cos(5 * deg) * std::cos(2 * deg), 60 * std::sin(2 * deg), 60 * std::sin(5 * deg) * std::cos(2 * deg)});
    s.check("air data: the airspeed", d.airspeed - 60, 1e-12, "m/s");
    s.check("air data: the angle of attack", d.alpha - 5 * deg, 1e-15, "rad");
    s.check("air data: the sideslip", d.beta - 2 * deg, 1e-15, "rad");
  }

  // ---- the body on the Earth ------------------------------------------------------------------------------------
  {
    const auto e = geodesy::figure::Ellipsoid::wgs84();
    const geodesy::figure::Geodetic here{48.85 * deg, 2.35 * deg, 120};
    const StateTransform t = body_to_earth(e, here, {90 * deg, 0, 0}, {0, 0, 0.01}, {0, 50, 0});
    const auto place = geodesy::figure::Place::at(e, here);
    const Vec3 east = place.enu_to_earth().direction({1, 0, 0});
    s.check("a vehicle heading east has its nose along the local east", max_abs(t.direction({1, 0, 0}) - east), 1e-15);
    s.check("its body origin is its place", max_abs(t.b - e.cartesian(here)), 1e-9, "m");
    s.check("its velocity in the Earth's frame is its NED velocity turned", max_abs(t.b_dot - 50.0 * east), 1e-12, "m/s");
    s.check("it turns about its own z (down) at its yaw rate", max_abs(t.angular_velocity() - 0.01 * t.direction({0, 0, 1})), 1e-17, "rad/s");
  }
  return s.done();
}
