// The kinematic vehicle: a unicycle on the sphere, commanded by its speed and its yaw rate
// (theory/vehicles.md §2).
//
// The fast model: what a vehicle that rolls without slipping, and whose motors keep the speed they are
// given, does. Its state is its pose (surface.hpp) and its forward speed; a command is a wished speed and a
// yaw rate. Each step the speed moves toward the command by at most the acceleration limit, the yaw rate is
// clamped to its limit, and the pose turns by S = exp(h (0, v̄/R, ω)^), v̄ the step's mean speed: exact when
// the speed is constant over the step, of order 2 in h while it slews (the curvature ω/v then varies within
// the step).
//
// Precision. Composed step after step, a full rotation gathers about ε of angle a step, R ε ≈ 4 × 10⁻¹⁰ m of
// place on the Moon. The place is therefore kept as the dynamic vehicle keeps it (dynamic.hpp): a local origin
// o on the sphere and the place x from it, which gathers only the steps' displacements, R A (S e_z − e_z),
// each computed without cancellation as R A (2w (s × e_z) + 2 s × (s × e_z)) from S's quaternion (w, s); the
// up axis is the place's direction and the forward axis is carried by A S and made tangent again; the origin
// moves to the vehicle when it has gone a kilometre. The two vehicles share this state's form, and the
// translation between them is exact.
//
// The dynamic model (dynamic.hpp) is this one's refinement: inertia, the ground's friction and a wheel motor's
// speed loop; this one is its limit as the loop's time constant goes to zero on firm ground.

#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <robotics/vehicles/surface.hpp>

namespace robotics::vehicles {

// What the vehicle is asked to do: a forward speed (m/s, negative backward) and a yaw rate (rad/s, positive
// to the left, about the up axis).
struct Command {
  double speed = 0;
  double yaw_rate = 0;
};

// What its motors allow.
struct Limits {
  double max_speed = 0;         // m/s, forward or backward
  double max_acceleration = 0;  // m/s²
  double max_yaw_rate = 0;      // rad/s
};

struct KinematicState {
  double t = 0;
  Vec3 origin{};        // a place on the sphere near the vehicle: its local origin, m
  Vec3 x{};             // the vehicle's place from the origin, m
  SurfacePose pose;     // its axes: forward, left, up
  double speed = 0;     // m/s, forward
  double odometer = 0;  // m driven, whichever way

  Vec3 place() const {
    using namespace physis::geometry::operators;
    return origin + x;
  }
};

class KinematicRover {
 public:
  KinematicRover(Sphere sphere, Limits limits, double recentre = 1000) : s_(sphere), l_(limits), recentre_(recentre) {
    if (!(sphere.radius > 0)) throw std::invalid_argument("robotics::vehicles::KinematicRover: a sphere of positive radius");
    if (!(limits.max_speed >= 0 && limits.max_acceleration > 0 && limits.max_yaw_rate >= 0))
      throw std::invalid_argument("robotics::vehicles::KinematicRover: limits of zero or more, and an acceleration above zero");
  }

  const Sphere& sphere() const { return s_; }
  const Limits& limits() const { return l_; }

  KinematicState start(const SurfacePose& pose, double speed, double t = 0) const {
    KinematicState x;
    x.t = t;
    x.origin = pose.place(s_);
    x.pose = pose;
    x.speed = speed;
    return x;
  }

  void step(KinematicState& x, const Command& c, double h) const {
    using namespace physis::geometry::operators;
    if (!(h > 0)) throw std::invalid_argument("robotics::vehicles::KinematicRover::step: a positive step");
    const double R = s_.radius;
    const double wished = std::clamp(c.speed, -l_.max_speed, l_.max_speed);
    const double dv = std::clamp(wished - x.speed, -l_.max_acceleration * h, l_.max_acceleration * h);
    const double mean = x.speed + 0.5 * dv;
    const double w = std::clamp(c.yaw_rate, -l_.max_yaw_rate, l_.max_yaw_rate);
    const Rotation A = x.pose.attitude();
    const Rotation S = Rotation::exp(Vec3{0, mean * h / R, w * h});
    const Vec3 s{S.x(), S.y(), S.z()}, ez{0, 0, 1};
    const Vec3 t = 2.0 * physis::geometry::cross(s, ez);
    const Vec3 rise = S.w() * t + physis::geometry::cross(s, t);  // S e_z − e_z
    x.x = x.x + R * A.apply(rise);
    x.pose = SurfacePose::at(x.place(), (A * S).apply(Vec3{1, 0, 0}));
    const double v0 = x.speed, v1 = x.speed + dv;
    // ∫|v| dt for v linear from v0 to v1: |v̄| h, or (v0² + v1²) h / (2 |v1 − v0|) when it changes sign.
    x.odometer += (v0 * v1 >= 0) ? std::fabs(mean) * h : (v0 * v0 + v1 * v1) * h / (2 * std::fabs(dv));
    x.speed = v1;
    x.t += h;
    if (physis::geometry::norm(x.x) > recentre_) {
      const Vec3 o2 = R * x.pose.up();
      x.x = x.x + (x.origin - o2);
      x.origin = o2;
    }
  }

 private:
  Sphere s_;
  Limits l_;
  double recentre_;
};

}  // namespace robotics::vehicles
