// Checks of the vehicles on a sphere (vehicles/; theory/vehicles.md §1–3), on the Moon's sphere.
//
// The kinematic vehicle:
// - driving straight it follows a great circle: the distance from the start is the speed times the time,
//   and it never leaves the plane of its first heading;
// - turning at constant speed and yaw rate it keeps a constant angle atan(v/(Rω)) to the fixed axis of its
//   rotation in space, and after 2π/|w|, |w| = √(ω² + (v/R)²) (not 2π/ω: the rotation also carries it round
//   the sphere), it is back where it began, facing the same way;
// - its speed slews within the acceleration limit, and the odometer counts the distance driven, both ways.
// The dynamic vehicle:
// - it stays on the sphere and its momentum tangent to it, through 3 km and two changes of origin, on the
//   great circle of its heading;
// - the ground carries its weight, N = m (g − v²/R);
// - coasting, rolling resistance stops it in v₀²/(2 C_rr g), and sideslip dies in v₀²/(2 μ_lat g);
// - its speed loop answers a step exactly as the discrete loop does, v_k = v*(1 − (1 − h/τ)ᵏ), until the
//   ground's traction limits it to the acceleration (μ − C_rr) g.
// The translation between them: lift then project returns the kinematic state, and driving straight the
// dynamic vehicle lags the kinematic one by v*(τ − h): its loop's lag v*τ, less a step (both schemes are
// trapezoidal in the speed).

#include <robotics/vehicles/dynamic.hpp>
#include <robotics/vehicles/kinematic.hpp>
#include <robotics/vehicles/surface.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>

using namespace robotics::vehicles;
using namespace physis::geometry::operators;
namespace pg = physis::geometry;

int main() {
  calx::verification::Suite s;
  const double pi = std::numbers::pi;
  const Sphere moon{1737.4e3, 4.9028e12};  // the mean radius and about DE440's GM of the Moon
  const SurfacePose start = SurfacePose::at(0.13, -1.03, 0.4);  // near Reiner Gamma, heading north-east

  // ---- the kinematic vehicle ----------------------------------------------------------------------------------
  {
    const KinematicRover k(moon, {0.2, 1e3, 0.1});
    KinematicState x = k.start(start, 0.2);
    for (int i = 0; i < 10000; i++) k.step(x, {0.2, 0}, 1.0);
    s.check("kinematic: driving straight 2 km, its distance along the sphere is v t (m)", std::fabs(distance(start.place(moon), x.place(), moon) - 2000), 1e-8, "m");
    s.check("kinematic: on the great circle of its first heading (its height above that plane, m)", std::fabs(pg::dot(x.place(), start.left())), 1e-8, "m");
    s.check("kinematic: its odometer is the distance driven (m)", std::fabs(x.odometer - 2000), 1e-9, "m");

    // A constant turn: the axis in space A w/|w| is fixed, the place keeps the angle atan(v/(Rω)) to it.
    const double v = 0.2, w = 0.01, T = 2 * pi / std::hypot(w, v / moon.radius);
    const Vec3 body{0, v / moon.radius, w};
    const Vec3 axis = pg::unit(start.attitude().apply(body));
    KinematicState c = k.start(start, v);
    double angle = 0;
    const int n = 1000;
    for (int i = 0; i < n; i++) {
      k.step(c, {v, w}, T / n);
      const Vec3 u = c.pose.up();  // the angle by atan2(|u × a|, u · a), which keeps its digits when it is small
      angle = std::max(angle, std::fabs(std::atan2(pg::norm(pg::cross(u, axis)), pg::dot(u, axis)) - std::atan(v / (moon.radius * w))));
    }
    s.check("kinematic: turning, its place keeps the angle atan(v/(Rω)) to the fixed axis of its rotation (rad)", angle, 1e-12, "rad");
    s.check("kinematic: after 2π/|w| it is back at its start (m)", distance(c.place(), start.place(moon), moon), 1e-9, "m");
    s.check("kinematic: facing the way it began", pg::norm(c.pose.forward() - start.forward()), 1e-12);

    // Speed slewing: 0 → 0.2 m/s at 0.01 m/s² covers ½ a t² = 2 m; 0.1 → −0.1 m/s goes 0.5 m out and back.
    const KinematicRover slow(moon, {0.2, 0.01, 0.1});
    KinematicState a = slow.start(start, 0);
    for (int i = 0; i < 20; i++) slow.step(a, {0.2, 0}, 1.0);
    s.check("kinematic: from rest to 0.2 m/s at 0.01 m/s², the speed reached in 20 s", std::fabs(a.speed - 0.2), 1e-15);
    s.check("kinematic: and ½ a t² = 2 m driven (m)", std::fabs(a.odometer - 2) + std::fabs(distance(a.place(), start.place(moon), moon) - 2), 1e-9, "m");
    KinematicState b = slow.start(start, 0.1);
    for (int i = 0; i < 7; i++) slow.step(b, {-0.1, 0}, 20.0 / 7);  // the sign changes within a step
    s.check("kinematic: reversing from 0.1 to −0.1 m/s, the odometer counts 0.5 m out and 0.5 m back (m)", std::fabs(b.odometer - 1.0), 1e-12, "m");
    s.check("kinematic: and it is back where it began (m)", distance(b.place(), start.place(moon), moon), 1e-9, "m");
  }

  // ---- the dynamic vehicle ------------------------------------------------------------------------------------------
  const Drive drive{140, 400, 200, 0.5, {1.0, 1e3, 0.2}};  // a 140 kg rover
  const Ground ground{0.5, 0.5, 0.1};
  {
    const DynamicRover d(moon, drive, ground);
    DynamicState x = d.start(start, 1.0);
    double off = 0, radial = 0, normal = 0;
    for (int i = 0; i < 60000; i++) {
      d.step(x, {1.0, 0}, 0.05);
      off = std::max(off, std::fabs(pg::norm(x.place()) - moon.radius));
      radial = std::max(radial, std::fabs(pg::dot(x.p, x.pose.up())) / (drive.mass * 1.0));
      if (i > 100) normal = std::max(normal, std::fabs(x.normal / (drive.mass * (moon.gravity() - 1.0 / moon.radius)) - 1));
    }
    s.check("dynamic: on the sphere through 3 km and two new origins (|q| − R, m)", off, 1e-8, "m");
    s.check("dynamic: its momentum tangent to it (radial speed over speed)", radial, 1e-9);
    s.check("dynamic: the ground carries it with N = m (g − v²/R) (relative)", normal, 1e-9);
    s.check("dynamic: on the great circle of its heading (its height above that plane, m)", std::fabs(pg::dot(x.place(), start.left())), 1e-6, "m");
    s.check("dynamic: its odometer is the distance along the sphere (relative)", std::fabs(x.odometer / distance(start.place(moon), x.place(), moon) - 1), 1e-9);
  }
  {
    // Coasting: the speed loop slackened (τ = 10⁹ s), the ground's friction alone.
    Drive slack = drive;
    slack.time_constant = 1e9;
    const DynamicRover d(moon, slack, ground);
    const double v0 = 0.2, g = moon.gravity(), h = 1e-3;
    DynamicState x = d.start(start, v0);
    while (d.forward_speed(x) > 1e-12 && x.t < 60) d.step(x, {0, 0}, h);  // the stop is exact: the speed then rounding's
    s.check("dynamic: coasting, rolling resistance stops it in v₀²/(2 C_rr g) (relative)", std::fabs(x.odometer / (v0 * v0 / (2 * ground.rolling * g)) - 1), 1e-3);
    // Sideslip: started moving to its left at 0.1 m/s.
    DynamicState y = d.start(start, 0);
    y.p = (slack.mass * 0.1) * start.left();
    while (std::fabs(d.lateral_speed(y)) > 1e-12 && y.t < 60) d.step(y, {0, 0}, h);
    s.check("dynamic: sideslip dies in v₀²/(2 μ_lat g) (relative)", std::fabs(y.odometer / (0.01 / (2 * ground.lateral * g)) - 1), 1e-3);
  }
  {
    // The speed loop's step response, and the traction limit.
    const DynamicRover d(moon, drive, ground);
    const double h = 0.01, target = 0.2;
    DynamicState x = d.start(start, 0);
    double worst = 0;
    for (int k = 1; k <= 400; k++) {
      d.step(x, {target, 0}, h);
      worst = std::max(worst, std::fabs(d.forward_speed(x) - target * (1 - std::pow(1 - h / drive.time_constant, k))));
    }
    s.check("dynamic: its speed loop answers a step as the discrete loop, v_k = v*(1 − (1 − h/τ)ᵏ) (m/s)", worst, 1e-9, "m/s");
    Drive eager = drive;
    eager.time_constant = 0.01;
    const DynamicRover e(moon, eager, ground);
    DynamicState y = e.start(start, 0);
    for (int k = 0; k < 10; k++) e.step(y, {target, 0}, h);
    s.check("dynamic: asked too much, it accelerates at the traction's (μ − C_rr) g (m/s after 0.1 s)",
            std::fabs(e.forward_speed(y) - (ground.traction - ground.rolling) * moon.gravity() * 0.1), 1e-9, "m/s");
  }

  // ---- between the fidelities -----------------------------------------------------------------------------------------
  {
    const DynamicRover d(moon, drive, ground);
    const KinematicRover k(moon, {1.0, 1e3, 0.2});
    KinematicState k0 = k.start(start, 0.15);
    k0.odometer = 3.0;
    const KinematicState back = d.project(d.lift(k0));
    s.check("translation: lift, then project, returns the kinematic state", std::max({pg::norm(back.pose.up() - k0.pose.up()), pg::norm(back.pose.forward() - k0.pose.forward()),
                                                                                        std::fabs(back.speed - k0.speed), std::fabs(back.odometer - k0.odometer)}), 1e-15);
    const double h = 0.01, target = 0.2;
    KinematicState a = k.start(start, 0);
    DynamicState b = d.start(start, 0);
    for (int i = 0; i < 3000; i++) k.step(a, {target, 0}, h), d.step(b, {target, 0}, h);
    const double lag = distance(a.place(), b.place(), moon);
    s.check("translation: driving straight, the dynamic vehicle lags the kinematic one by v*(τ − h) (m)", std::fabs(lag - target * (drive.time_constant - h)), 1e-6, "m");
  }

  return s.done();
}
