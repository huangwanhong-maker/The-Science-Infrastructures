// Checks of the fish (swimming/fish.hpp, ADR 0059, theory/swimming.md).
//
// - The body: its mass and length; each joint's reduced inertia with the water below its diagonal one;
//   holding straight in still water it stays where it is (it is neutrally buoyant).
// - Swimming at 0.15 m/s asked: its speed reached (the mean over the last 8 s of 20), its heading and course
//   held; the wake it leaves pushes the water backward (the mean force the fin and the peduncle give it).
// - Turning left; stopping to a glide.
// - The three methods are one motion: the distance swum in 6 s agrees between them; each states its energy's
//   account and its cost.
// - Holding still in a current across it, it is carried with the water.

#include <robotics/swimming/fish.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <string>

namespace rs = robotics::swimming;

namespace {

struct Swim {
  double distance = 0, mean_speed = 0, max_heading = 0, lateral = 0, residual = 0, cost = 0, wake_back = 0;
};

Swim swim(rs::Method method, double seconds, double speed) {
  rs::SimulationConfig c;
  c.method = method;
  rs::Simulation sim(rs::fish(), c);
  sim.command({rs::Behaviour::swim, speed, 0.0});
  Swim r;
  double x0 = 0, t0 = 0;
  bool late = false;
  while (sim.time() < seconds - 1e-9) {
    sim.advance(0.1);
    if (sim.time() > 4.0) r.max_heading = std::max(r.max_heading, std::fabs(sim.q()[rs::HEADING]));
    if (!late && sim.time() >= seconds - 8.0) late = true, x0 = sim.q()[rs::X], t0 = sim.time();
  }
  r.distance = sim.q()[rs::X];
  r.mean_speed = (sim.q()[rs::X] - x0) / (sim.time() - t0);
  r.lateral = sim.q()[rs::Y];
  r.residual = sim.energy_residual();
  r.cost = sim.cost();
  // The wake: the mean along x of the force the peduncle and the fin give the water, over the frames kept.
  double fx = 0;
  int n = 0;
  for (const auto& f : sim.frames())
    if (f.t > 4.0) {
      for (const auto& w : f.wake) fx += w.fx;
      n++;
    }
  r.wake_back = n ? fx / n : 0;
  return r;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const auto fish = rs::fish();
  s.check("body: the fish's mass, 1.85 kg (kg)", fish.mass() - 1.85, 1e-12);
  s.check("body: and its length, 0.68 m (m)", fish.length() - 0.68, 1e-12);
  {
    rs::SimulationConfig c;
    rs::Simulation sim(fish, c);
    const auto I = sim.joint_inertia();
    bool ok = true;
    for (double x : I) ok = ok && x > 0;
    s.require("body: each joint's reduced inertia with the water's added mass is positive", ok);
    sim.advance(2.0);
    s.check("still: holding straight in still water it stays where it is (neutrally buoyant; m)", std::hypot(sim.q()[rs::X], sim.q()[rs::Y]), 1e-9);
  }
  {
    const auto r = swim(rs::Method::semi_implicit, 20.0, 0.15);
    s.check("swimming: asked 0.15 m/s, it reaches it (the mean over the last 8 s of 20, m/s)", r.mean_speed - 0.15, 0.04, "m/s");
    s.check("swimming: its heading held within 15° once under way (rad)", r.max_heading, 15 * std::numbers::pi / 180);
    s.check("swimming: and its course: it ends within 0.3 m of its line (m)", r.lateral, 0.3, "m");
    s.require("swimming: its wake pushes the water backward (the mean force its tail gives the water, along −x)", r.wake_back < 0, std::to_string(r.wake_back) + " N");
  }
  {
    rs::SimulationConfig c;
    rs::Simulation sim(fish, c);
    sim.command({rs::Behaviour::swim, 0.15, 0.0});
    sim.advance(6.0);
    const double h0 = sim.q()[rs::HEADING];
    sim.command({rs::Behaviour::swim, 0.15, 1.0});
    sim.advance(8.0);
    s.require("turning: asked to turn left, its heading turns more than 45° in 8 s", sim.q()[rs::HEADING] - h0 > std::numbers::pi / 4,
              std::to_string((sim.q()[rs::HEADING] - h0) * 180 / std::numbers::pi) + "°");
    sim.command({rs::Behaviour::stop, 0.0, 0.0});
    sim.advance(10.0);
    s.check("stopping: it straightens and glides to rest (its speed after 10 s, m/s)", std::hypot(sim.v()[rs::X], sim.v()[rs::Y]), 0.06, "m/s");
  }
  {
    const auto e = swim(rs::Method::semi_implicit, 6.0, 0.15), k = swim(rs::Method::runge_kutta, 6.0, 0.15), v = swim(rs::Method::variational, 6.0, 0.15);
    s.check("methods: the distance swum in 6 s, semi-implicit Euler against Runge–Kutta (relative)", e.distance / k.distance - 1, 0.05);
    s.check("methods: and the variational integrator against Runge–Kutta (relative)", v.distance / k.distance - 1, 0.05);
    s.check("methods: Runge–Kutta's energy account, the body's and the water's (J)", k.residual, 0.05, "J");
    s.require("methods: each states its cost, the machine's seconds per second of the fish's", e.cost > 0 && k.cost > 0 && v.cost > 0,
              std::to_string(e.cost) + ", " + std::to_string(k.cost) + ", " + std::to_string(v.cost) + " s/s");
  }
  {
    rs::SimulationConfig c;
    c.current = {0.0, 0.2};
    rs::Simulation sim(fish, c);
    sim.advance(30.0);
    s.check("current: holding still in a current of 0.2 m/s across it, it is carried with the water (its velocity across, m/s)", sim.v()[rs::Y] - 0.2, 0.03, "m/s");
  }
  return s.done();
}
