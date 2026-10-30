// Checks of the dog (quadruped/dog.hpp, ADR 0059, theory/quadruped.md).
//
// - The body: its mass; each leg's inverse kinematics against its forward kinematics.
// - Standing on rock: still, level, its weight carried by its four feet, its energy accounted for.
// - Trotting: in place it stays near where it was, level; forward it reaches the speed asked, level and at its
//   height; asked to stop, it stands.
// - The three methods are one motion: the distance trotted in 4 s agrees between them.
// - On ice (μ 0.15), told so by its presumption, it asks no more friction than it believes is there: it gathers
//   speed more gently than on rock, and trots on it level.

#include <robotics/quadruped/dog.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace rq = robotics::quadruped;

namespace {

struct Trot {
  double distance = 0, mean_speed = 0, max_pitch = 0, min_height = 1e9, max_height = 0, residual = 0, cost = 0;
  std::string behaviour;
};

Trot trot(rq::Method method, double seconds, double speed, double friction = 0.8, double presumed = 0.8295) {
  rq::SimulationConfig c;
  c.method = method;
  c.ground_at = [friction](double) { auto r = rq::rock(); r.friction = friction; return r; };
  c.presumption.friction = {presumed, 0.185};
  rq::Simulation sim(rq::dog(), c);
  sim.advance(0.5);
  sim.command({rq::Behaviour::trot, speed});
  Trot r;
  double x0 = 0, t0 = 0;
  bool late = false;
  while (sim.time() < 0.5 + seconds - 1e-9) {
    sim.advance(0.05);
    if (sim.time() > 1.5) r.max_pitch = std::max(r.max_pitch, std::fabs(sim.q()[rq::PITCH])), r.min_height = std::min(r.min_height, sim.q()[rq::Z]),
                          r.max_height = std::max(r.max_height, sim.q()[rq::Z]);
    if (!late && sim.time() >= 0.5 + seconds - 4.0) late = true, x0 = sim.q()[rq::X], t0 = sim.time();
  }
  r.distance = sim.q()[rq::X];
  r.mean_speed = (sim.q()[rq::X] - x0) / (sim.time() - t0);
  r.residual = sim.energy_residual();
  r.cost = sim.cost();
  r.behaviour = rq::behaviour_name(sim.controller().current().behaviour);
  return r;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const auto dog = rq::dog();
  s.check("body: the dog's mass, 10.2 kg (kg)", dog.mass() - 10.2, 1e-12);
  {
    double worst = 0;
    for (const auto& T : std::vector<std::array<double, 2>>{{0.0, -0.30}, {0.08, -0.28}, {-0.10, -0.25}, {0.05, -0.18}}) {
      const auto ik = rq::leg_ik(0.20, 0.20, T);
      if (!ik) { worst = 1; continue; }
      worst = std::max({worst, std::fabs(0.20 * std::cos((*ik)[0]) + 0.20 * std::cos((*ik)[1]) - T[0]), std::fabs(0.20 * std::sin((*ik)[0]) + 0.20 * std::sin((*ik)[1]) - T[1])});
    }
    s.check("legs: each leg's inverse kinematics reaches its target, the knee behind (m)", worst, 1e-15);
  }
  {
    rq::SimulationConfig c;
    rq::Simulation sim(dog, c);
    sim.advance(0.5);
    const double z1 = sim.q()[rq::Z];
    sim.advance(1.0);
    double load = 0;
    for (const auto& f : sim.contacts()) load += f.force[1];
    s.check("standing: on rock it is still (its torso's height over a second, m)", sim.q()[rq::Z] - z1, 1e-4);
    s.check("standing: level (rad)", sim.q()[rq::PITCH], 0.01);
    s.check("standing: its four feet carry its weight (relative)", load / (dog.mass() * 9.81) - 1, 0.01);
    s.check("standing: its energy accounted for (J)", sim.energy_residual(), 0.1);
  }
  {
    const auto r = trot(rq::Method::semi_implicit, 6.0, 0.0);
    s.check("trotting in place: it stays near where it stood (m)", r.distance, 0.15, "m");
    s.check("trotting in place: level within 5° (rad)", r.max_pitch, 5 * rq::DEG);
  }
  {
    const auto r = trot(rq::Method::semi_implicit, 8.0, 0.3);
    s.check("trotting forward: the speed asked, 0.3 m/s, reached (its mean over the last 4 s, m/s)", r.mean_speed - 0.3, 0.05, "m/s");
    s.check("trotting forward: level within 5° (rad)", r.max_pitch, 5 * rq::DEG);
    s.check("trotting forward: its torso held at its height, within 3 cm of 0.30 m (m)", std::max(std::fabs(r.min_height - 0.30), std::fabs(r.max_height - 0.30)), 0.03);
  }
  {
    rq::SimulationConfig c;
    rq::Simulation sim(dog, c);
    sim.advance(0.5);
    sim.command({rq::Behaviour::trot, 0.3});
    sim.advance(5.0);
    sim.command({rq::Behaviour::stop, 0.0});
    sim.advance(3.0);
    s.require("stopping: from a trot at 0.3 m/s it stops and stands", std::string(rq::behaviour_name(sim.controller().current().behaviour)) == "stand",
              rq::behaviour_name(sim.controller().current().behaviour));
    s.check("stopping: still and level when it stands (the larger of its speed and its pitch)", std::max(std::fabs(sim.v()[rq::X]), std::fabs(sim.q()[rq::PITCH])), 0.03);
  }
  {
    const auto e = trot(rq::Method::semi_implicit, 4.0, 0.3), k = trot(rq::Method::runge_kutta, 4.0, 0.3), v = trot(rq::Method::variational, 4.0, 0.3);
    s.check("methods: the distance trotted in 4 s, semi-implicit Euler against Runge–Kutta (relative)", e.distance / k.distance - 1, 0.05);
    s.check("methods: and the variational integrator against Runge–Kutta (relative)", v.distance / k.distance - 1, 0.05);
    s.require("methods: each states its energy's account and its cost",
              std::isfinite(e.residual) && std::isfinite(k.residual) && std::isfinite(v.residual) && e.cost > 0 && k.cost > 0 && v.cost > 0,
              std::to_string(e.residual) + ", " + std::to_string(k.residual) + ", " + std::to_string(v.residual) + " J; " + std::to_string(e.cost) + ", " +
                  std::to_string(k.cost) + ", " + std::to_string(v.cost) + " s/s");
  }
  {
    // A steady trot asks little of the ground; gathering speed asks more, and on ice it asks no more than its
    // belief allows: it gathers speed more gently than on rock.
    const auto ice = trot(rq::Method::semi_implicit, 1.5, 0.3, 0.15, 0.15), rock = trot(rq::Method::semi_implicit, 1.5, 0.3);
    s.require("ice: told it is on ice (μ 0.15), it gathers speed more gently: less distance in its first 1.5 s than on rock", ice.distance < rock.distance,
              std::to_string(ice.distance) + " m against " + std::to_string(rock.distance) + " m");
    const auto r = trot(rq::Method::semi_implicit, 8.0, 0.3, 0.15, 0.15);
    s.check("ice: and trots on it, level within 5° (rad)", r.max_pitch, 5 * rq::DEG);
  }
  return s.done();
}
