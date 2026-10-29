// Checks of the waypoint follower (control/waypoints.hpp; theory/vehicles.md §4), on the Moon's sphere.
//
// - The heading error's sign and size: a place to the left is positive, one straight ahead zero, one
//   behind ±π, to the precision of a bearing between unit vectors d apart, R ε/d (7.6 × 10⁻¹² rad at 50 m).
// - Both vehicles, given the same follower, reach a waypoint 100 m away that starts behind them to the
//   right, within the arrival radius, turning in place first (the speed is zero while the waypoint is
//   behind).
// - Both drive a lawnmower of three lines and reach every waypoint in order.

#include <robotics/control/waypoints.hpp>
#include <robotics/planning/survey.hpp>
#include <robotics/vehicles/dynamic.hpp>
#include <robotics/vehicles/kinematic.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

using namespace robotics;
using namespace physis::geometry::operators;
using vehicles::Vec3;

int main() {
  calx::verification::Suite s;
  const double pi = std::numbers::pi;
  const vehicles::Sphere moon{1737.4e3, 4.9028e12};
  const auto start = vehicles::SurfacePose::at(0.13, -1.03, 0.0);  // facing north
  const auto ahead = [&](double forward, double left) { return vehicles::exponential(start.place(moon), forward * start.forward() + left * start.left(), moon); };

  s.check("heading: a place straight ahead has error 0", std::fabs(control::WaypointFollower::heading_error(start, ahead(50, 0))), 1e-10, "rad");
  s.check("heading: one to the left, at 45°, has error +π/4", std::fabs(control::WaypointFollower::heading_error(start, ahead(50, 50)) - pi / 4), 1e-9, "rad");
  s.check("heading: one behind has error ±π", std::fabs(std::fabs(control::WaypointFollower::heading_error(start, ahead(-50, 1e-9))) - pi), 1e-9, "rad");

  const vehicles::Limits limits{0.3, 0.5, 0.3};
  const vehicles::KinematicRover kin(moon, limits);
  const vehicles::DynamicRover dyn(moon, {140, 400, 200, 0.5, limits}, {0.5, 0.5, 0.1});
  const control::FollowerGains gains{0.2, 0.8, 0.2, 0.5};
  const double h = 0.05;

  // One waypoint behind, to the right.
  {
    const Vec3 target = ahead(-80, -60);
    control::WaypointFollower fk({target}, moon, gains), fd({target}, moon, gains);
    auto a = kin.start(start, 0);
    auto b = dyn.start(start, 0);
    double moved_while_behind = 0;
    for (int i = 0; i < 40000 && !(fk.done() && fd.done()); i++) {
      const double behind = std::fabs(control::WaypointFollower::heading_error(a.pose, target)) > pi / 2;
      const auto ca = fk.command(a.pose), cb = fd.command(b.pose);
      kin.step(a, ca, h), dyn.step(b, cb, h);
      if (behind) moved_while_behind = std::max(moved_while_behind, std::fabs(ca.speed));
    }
    s.require("follower: the kinematic vehicle reaches a waypoint 100 m away, behind it", fk.done());
    s.require("follower: and so does the dynamic one", fd.done());
    s.check("follower: within the arrival radius (kinematic, m)", vehicles::distance(a.place(), target, moon), gains.arrival + 0.01, "m");
    s.check("follower: within the arrival radius (dynamic, m)", vehicles::distance(b.place(), target, moon), gains.arrival + 0.01, "m");
    s.check("follower: while the waypoint is behind it turns in place (its commanded speed, m/s)", moved_while_behind, 0.0, "m/s");
  }

  // A lawnmower of three lines, 60 m long and 10 m apart, north of the start.
  {
    const auto lines = planning::lawnmower(ahead(40, 0), 0.0, 60, 10, 3, moon);
    control::WaypointFollower fk(lines, moon, gains), fd(lines, moon, gains);
    auto a = kin.start(start, 0);
    auto b = dyn.start(start, 0);
    std::size_t ka = 0, kb = 0;
    bool in_order = true;
    for (int i = 0; i < 200000 && !(fk.done() && fd.done()); i++) {
      kin.step(a, fk.command(a.pose), h), dyn.step(b, fd.command(b.pose), h);
      if (fk.next() < ka || fd.next() < kb) in_order = false;
      ka = fk.next(), kb = fd.next();
    }
    s.require("survey: both vehicles drive the lawnmower to its end, every waypoint in order", fk.done() && fd.done() && in_order,
              std::to_string(fk.next()) + " and " + std::to_string(fd.next()) + " of " + std::to_string(lines.size()));
  }

  return s.done();
}
