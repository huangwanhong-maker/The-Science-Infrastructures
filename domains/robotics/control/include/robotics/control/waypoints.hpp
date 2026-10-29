// Following waypoints on a sphere: the controller that turns a vehicle toward its next place and drives it
// there (theory/vehicles.md §4).
//
// A controller sees what its vehicle believes, not the world: it is given a pose (from the vehicle's own
// estimate) and returns a Command, the same for the kinematic and the dynamic vehicle. Toward a waypoint w
// from the pose (p, e):
//   - the bearing b is the direction at p of the great circle to w (surface.hpp), and the heading error the
//     angle from e to b about the up axis, ψ = atan2((e × b) · û, e · b);
//   - the yaw rate is k_ψ ψ, within the vehicle's limit: a first-order approach of the heading to the
//     bearing, with time constant 1/k_ψ;
//   - the speed is the cruising speed times max(0, cos ψ), so that a vehicle facing away turns in place
//     before it goes, and at most k_d times the distance left, so that it slows as it arrives;
//   - within the arrival radius the waypoint is reached and the next one taken.
// It is a proportional heading controller on the sphere's geometry, a relative of pure pursuit (which
// commands the curvature 2 sin ψ / L to a point at the lookahead L; Coulter 1992, CMU-RI-TR-92-01, §2). The
// heading error decays at the rate k_ψ while the bearing turns more slowly than that, and the distance left
// at the rate k_d near the end, so a vehicle that turns as commanded comes within the arrival radius in a
// finite time; the checks drive both vehicles there.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <robotics/vehicles/kinematic.hpp>
#include <robotics/vehicles/surface.hpp>

namespace robotics::control {

using robotics::vehicles::Command;
using robotics::vehicles::Sphere;
using robotics::vehicles::SurfacePose;
using robotics::vehicles::Vec3;

struct FollowerGains {
  double cruise = 0.1;     // m/s
  double heading = 0.5;    // k_ψ, 1/s
  double approach = 0.2;   // k_d, 1/s
  double arrival = 0.5;    // m
};

class WaypointFollower {
 public:
  WaypointFollower(std::vector<Vec3> waypoints, Sphere sphere, FollowerGains gains) : w_(std::move(waypoints)), s_(sphere), k_(gains) {
    if (!(gains.cruise >= 0 && gains.heading > 0 && gains.approach > 0 && gains.arrival > 0))
      throw std::invalid_argument("robotics::control::WaypointFollower: positive gains and arrival radius");
  }

  const std::vector<Vec3>& waypoints() const { return w_; }
  std::size_t next() const { return next_; }
  bool done() const { return next_ >= w_.size(); }

  // The heading error toward a place, radians in (−π, π], positive when the place is to the left.
  static double heading_error(const SurfacePose& pose, const Vec3& target) {
    const Vec3 up = pose.up();
    const Vec3 b = robotics::vehicles::bearing(up, target);
    const Vec3 e = pose.forward();
    return std::atan2(physis::geometry::dot(physis::geometry::cross(e, b), up), physis::geometry::dot(e, b));
  }

  // The command at a pose; waypoints within the arrival radius are passed first.
  Command command(const SurfacePose& pose) {
    const Vec3 here = pose.place(s_);
    while (!done() && robotics::vehicles::distance(here, w_[next_], s_) <= k_.arrival) next_++;
    if (done()) return {};
    const double left = robotics::vehicles::distance(here, w_[next_], s_);
    const double psi = heading_error(pose, w_[next_]);
    Command c;
    c.yaw_rate = k_.heading * psi;
    c.speed = std::min(k_.cruise * std::max(0.0, std::cos(psi)), k_.approach * left);
    return c;
  }

 private:
  std::vector<Vec3> w_;
  Sphere s_;
  FollowerGains k_;
  std::size_t next_ = 0;
};

}  // namespace robotics::control
