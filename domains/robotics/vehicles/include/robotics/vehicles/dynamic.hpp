// The dynamic vehicle: a mass held to the sphere by its weight and the ground, driven by its wheels against
// the ground's friction (theory/vehicles.md §3).
//
// The accurate model. The vehicle is a point mass m at its place q on the sphere, kept there by the
// holonomic constraint g(q) = 0, whose multiplier is the ground's normal force N; gravity, −GM m q/|q|³, acts
// on it, the wheels' drive force F along the heading, and the ground's Coulomb friction against its rolling
// and its sideslip. The motion is calx's RATTLE (symplectic on the constraint, of order 2) for gravity and
// the constraint, with the wheels and the ground split off on either side of it (Strang) as two half kicks,
// each solved exactly for the force held over the step:
//   - forward, v′ = (F − R sign v)/m with R = C_rr N the rolling resistance: the speed moves at that rate,
//     and if it reaches zero within the kick it stops there, staying stopped while |F| ≤ R (stiction) and
//     moving off at (F ∓ R)/m when the drive overcomes it;
//   - sideways, v′ = −μ_lat N sign v / m: the speed falls at that rate until it is zero and stays there.
// The drift between the kicks moves the vehicle at the speed after the first, so a speed that changes at a
// constant rate is integrated by the trapezoidal rule.
//
// The wheels. A motor's speed loop gives the drive force F = m (v* − v)/τ + C_rr N sign(v*), the wished
// speed's error over the loop's time constant τ and the rolling resistance it must overcome, limited by the
// traction the ground allows, μ N, the motors' force and their power over the speed. The heading is turned
// at the commanded yaw rate (a skid-steered vehicle's yaw is much faster than its translation), carried with
// the place by the rotation taking one up to the next (the discrete parallel transport along the step's
// great circle) and made tangent again.
//
// Precision. Far from the sphere's centre (1.7 × 10⁶ m on the Moon) |q|² − R² cannot be met to better than
// ε R² ≈ 10⁻³ m², and a momentum carries the rounding of the place over h (physis mechanics,
// variational.hpp). The vehicle is integrated about a local origin o on the sphere, q = o + x, where the
// constraint is
//   g(x) = (|x|² + 2 o · x)/(2R)   (= |o + x| − R to first order: a height, in metres),   ∇g = (o + x)/R,
// with no cancellation between large numbers; the origin moves to the vehicle when it has gone a
// kilometre.
//
// The kinematic model (kinematic.hpp) is this one's limit as τ → 0 on firm ground (no slip); the
// translation between them (project, lift) carries the pose and the speed, and the normal force is the
// weight less the centripetal term, m (g − v²/R).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <calx/integrators/symplectic.hpp>
#include <robotics/vehicles/kinematic.hpp>
#include <robotics/vehicles/surface.hpp>

namespace robotics::vehicles {

// The ground under the wheels.
struct Ground {
  double traction = 0.5;  // μ: the largest drive or braking force over N
  double lateral = 0.5;   // μ_lat: the largest sideways force over N
  double rolling = 0.1;   // C_rr: the rolling resistance over N
};

// The vehicle's mass and its motors.
struct Drive {
  double mass = 0;           // kg
  double max_force = 0;      // N, the wheels' force together
  double max_power = 0;      // W
  double time_constant = 0;  // s, of the speed loop
  Limits limits;             // speed and yaw rate (its acceleration limit is the kinematic model's)
};

struct DynamicState {
  double t = 0;
  Vec3 origin{};      // a place on the sphere near the vehicle: its local origin, m
  Vec3 x{};           // the vehicle's place from the origin, m
  Vec3 p{};           // its momentum, kg m/s (the sphere frame's axes)
  SurfacePose pose;   // its axes: forward, left, up
  double normal = 0;  // the ground's normal force over the last step, N
  double odometer = 0;

  Vec3 place() const {
    using namespace physis::geometry::operators;
    return origin + x;
  }
};

class DynamicRover {
 public:
  DynamicRover(Sphere sphere, Drive drive, Ground ground, double recentre = 1000) : s_(sphere), d_(drive), g_(ground), recentre_(recentre) {
    if (!(sphere.radius > 0 && sphere.gm > 0)) throw std::invalid_argument("robotics::vehicles::DynamicRover: a sphere of positive radius and GM");
    if (!(drive.mass > 0 && drive.max_force > 0 && drive.max_power > 0 && drive.time_constant > 0))
      throw std::invalid_argument("robotics::vehicles::DynamicRover: a positive mass, force, power and time constant");
  }

  const Sphere& sphere() const { return s_; }
  const Drive& drive() const { return d_; }
  const Ground& ground() const { return g_; }

  // The state at a pose, moving forward at a speed: the lift of a kinematic state.
  DynamicState start(const SurfacePose& pose, double speed, double t = 0) const {
    using namespace physis::geometry::operators;
    DynamicState x;
    x.t = t;
    x.origin = pose.place(s_);
    x.p = (d_.mass * speed) * pose.forward();
    x.pose = pose;
    x.normal = d_.mass * (s_.gravity() - speed * speed / s_.radius);
    return x;
  }
  DynamicState lift(const KinematicState& k) const {
    using namespace physis::geometry::operators;
    DynamicState x;
    x.t = k.t;
    x.origin = k.origin;
    x.x = k.x;
    x.p = (d_.mass * k.speed) * k.pose.forward();
    x.pose = k.pose;
    x.normal = d_.mass * (s_.gravity() - k.speed * k.speed / s_.radius);
    x.odometer = k.odometer;
    return x;
  }
  // Its kinematic state: the pose and the forward speed.
  KinematicState project(const DynamicState& x) const {
    KinematicState k;
    k.t = x.t;
    k.origin = x.origin;
    k.x = x.x;
    k.pose = x.pose;
    k.speed = forward_speed(x);
    k.odometer = x.odometer;
    return k;
  }
  double forward_speed(const DynamicState& x) const { return physis::geometry::dot(x.p, x.pose.forward()) / d_.mass; }
  double lateral_speed(const DynamicState& x) const { return physis::geometry::dot(x.p, x.pose.left()) / d_.mass; }

  // The drive force the speed loop asks for at a forward speed, within traction, force and power.
  double drive_force(double speed, double wished, double normal) const {
    const double m = d_.mass;
    const double feedforward = wished > 0 ? g_.rolling * normal : wished < 0 ? -g_.rolling * normal : 0.0;
    const double F = m * (wished - speed) / d_.time_constant + feedforward;
    const double cap = std::min({g_.traction * std::max(normal, 0.0), d_.max_force, d_.max_power / std::max(std::fabs(speed), 1e-9)});
    return std::clamp(F, -cap, cap);
  }

  void step(DynamicState& x, const Command& c, double h) const {
    using namespace physis::geometry::operators;
    if (!(h > 0)) throw std::invalid_argument("robotics::vehicles::DynamicRover::step: a positive step");
    const double m = d_.mass, R = s_.radius;
    const double wished = std::clamp(c.speed, -d_.limits.max_speed, d_.limits.max_speed);
    const double yaw = std::clamp(c.yaw_rate, -d_.limits.max_yaw_rate, d_.limits.max_yaw_rate);
    const double N = x.normal > 0 ? x.normal : 0.0;
    const double F = drive_force(forward_speed(x), wished, N);

    kick(x, F, N, 0.5 * h);
    const Vec3 o = x.origin;
    const double GMm = s_.gm * m;
    const auto gradient = [&](const calx::integrators::Vec<3>& q, calx::integrators::Vec<3>& g) {
      const Vec3 r = o + Vec3{q[0], q[1], q[2]};
      const double r3 = std::pow(physis::geometry::dot(r, r), 1.5);
      for (int i = 0; i < 3; i++) g[i] = GMm * r[i] / r3;
    };
    const auto constraint = [&](const calx::integrators::Vec<3>& q, calx::integrators::Vec<1>& v, std::array<calx::integrators::Vec<3>, 1>& G) {
      v[0] = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + 2 * (o[0] * q[0] + o[1] * q[1] + o[2] * q[2])) / (2 * R);
      for (int i = 0; i < 3; i++) G[0][i] = (o[i] + q[i]) / R;
    };
    const Vec3 before = x.place(), local = x.x;
    calx::integrators::Vec<3> q{x.x[0], x.x[1], x.x[2]}, p{x.p[0], x.p[1], x.p[2]};
    const double inv = 1 / m;
    const auto r = calx::integrators::rattle<3, 1>({inv, inv, inv}, gradient, constraint, q, p, h, 1e-9);
    x.x = {q[0], q[1], q[2]};
    x.p = {p[0], p[1], p[2]};
    // The constraint's force is −∇g λ over the first half and −∇g μ over the second; ∇g is the outward
    // normal, so the ground pushes outward with N = −(λ + μ)/2.
    x.normal = -0.5 * (r.lambda[0] + r.mu[0]);
    // The heading carried to the new place, turned by the yaw, made tangent again.
    const Vec3 after = x.place();
    const Vec3 u0 = physis::geometry::unit(before), u1 = physis::geometry::unit(after);
    const Vec3 axis = physis::geometry::cross(u0, u1);
    const double angle = std::atan2(physis::geometry::norm(axis), physis::geometry::dot(u0, u1));
    Vec3 forward = x.pose.forward();
    if (angle > 0) forward = Rotation::about(axis, angle).apply(forward);
    forward = Rotation::about(u1, yaw * h).apply(forward);
    x.pose = SurfacePose::at(u1, forward);
    // The arc of the step from its chord in the local coordinates, which carry it without the rounding of R.
    const double chord = physis::geometry::norm(x.x - local);
    x.odometer += 2 * R * std::asin(std::min(1.0, chord / (2 * R)));
    kick(x, F, x.normal > 0 ? x.normal : 0.0, 0.5 * h);
    x.t += h;
    // A new origin when the vehicle has gone far from it.
    if (physis::geometry::norm(x.x) > recentre_) {
      const Vec3 o2 = R * u1;
      x.x = x.x + (x.origin - o2);
      x.origin = o2;
    }
  }

 private:
  // v after a time τ of v′ = (F − R sign v)/m, R ≥ 0 a Coulomb friction, F held: exact, with stiction at zero.
  static double coulomb(double v, double F, double R, double m, double tau) {
    if (v == 0) {
      if (std::fabs(F) <= R) return 0;
      return (F - std::copysign(R, F)) / m * tau;
    }
    const double a = (F - std::copysign(R, v)) / m;
    const double end = v + a * tau;
    if ((end > 0) == (v > 0) && end != 0) return end;  // no stop within τ
    // It stops at t₀ = −v/a, then stays, or moves off the other way if the drive overcomes the friction.
    const double rest = tau + v / a;
    return coulomb(0, F, R, m, rest);
  }
  // The wheels and the ground over a time τ: forward, the drive against rolling resistance; sideways, the
  // lateral friction alone.
  void kick(DynamicState& x, double F, double N, double tau) const {
    using namespace physis::geometry::operators;
    const double m = d_.mass;
    const Vec3 e = x.pose.forward(), n = x.pose.left(), u = x.pose.up();
    const double vf = physis::geometry::dot(x.p, e) / m, vl = physis::geometry::dot(x.p, n) / m, vu = physis::geometry::dot(x.p, u) / m;
    const double f = coulomb(vf, F, g_.rolling * N, m, tau), l = coulomb(vl, 0, g_.lateral * N, m, tau);
    x.p = m * (f * e + l * n + vu * u);
  }

  Sphere s_;
  Drive d_;
  Ground g_;
  double recentre_;
};

}  // namespace robotics::vehicles
