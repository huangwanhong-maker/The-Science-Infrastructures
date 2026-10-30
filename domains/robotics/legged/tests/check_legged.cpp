// Checks of the rabbit (legged/rabbit.hpp, controller.hpp, simulation.hpp; ADR 0058, theory/legged.md).
//
// - the body: its mass, the legs' inverse kinematics against their forward kinematics, the standing
//   posture's feet on the ground, a motor's torque–speed line, the potential's force against differences;
// - standing on the Moon: still, level, its feet carrying its weight, its energy accounted for;
// - hopping in place: sustained, the pitch held, the apex brought to the height asked, the speed near zero;
// - hopping forward at the speed asked (an integral of the speed's error), and stopping to sit;
// - the three methods are one motion: the first hop's apex agrees between them; the variational integrator
//   keeps the energy's account best and semi-implicit Euler worst (at contact's stiff and nonsmooth instants
//   an explicit method's order is not reached); Runge–Kutta keeps the angular momentum about the centre in
//   flight;
// - on the Earth the hop is refused, with its reason: this build's ankle cannot launch it there.

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include <calx/verification/check.hpp>
#include <physis/analysis/stationary_field.hpp>
#include <robotics/legged/controller.hpp>
#include <robotics/legged/rabbit.hpp>
#include <robotics/legged/simulation.hpp>

using calx::verification::Suite;
namespace rl = robotics::legged;

namespace {

struct Run {
  double max_pitch = 0, min_height = 1e9, apex = 0, speed = 0, residual = 0, flight_momentum = 0, cost = 0;
  int hops = 0;
  std::string phase;
};

Run hop(rl::Method method, double gravity, double seconds_stand, double seconds_hop, double speed, double height) {
  rl::SimulationConfig c;
  c.method = method, c.gravity = gravity;
  rl::Simulation sim(rl::rabbit(), c);
  sim.advance(seconds_stand);
  sim.command({rl::Behaviour::hop, speed, height});
  Run r;
  const double t0 = sim.time();
  double x_start = 0, t_start = 0;
  bool started = false;
  while (sim.time() < t0 + seconds_hop) {
    sim.advance(0.05);
    r.max_pitch = std::max(r.max_pitch, std::fabs(sim.q()[rl::PITCH]));
    if (!started && sim.time() > t0 + seconds_hop / 2) started = true, x_start = sim.q()[rl::X], t_start = sim.time();
  }
  r.speed = (sim.q()[rl::X] - x_start) / (sim.time() - t_start);
  r.apex = sim.controller().last_apex();
  r.hops = sim.controller().hops();
  r.residual = std::fabs(sim.energy_residual());
  r.flight_momentum = sim.flight_momentum_change();
  r.cost = sim.cost();
  r.phase = rl::phase_name(sim.controller().phase());
  return r;
}

}  // namespace

int main() {
  Suite s;
  const auto rabbit = rl::rabbit();
  constexpr double MOON = 1.62;

  // ---- the body -------------------------------------------------------------------------------------------
  s.check("body: the rabbit's mass, 1.69 kg (kg)", rabbit.mass() - 1.69, 1e-12);
  {
    double worst = 0;
    for (double fwd : {+1.0, -1.0})
      for (const auto& T : std::vector<std::array<double, 2>>{{0.05, -0.15}, {-0.08, -0.12}, {0.12, -0.05}, {0.0, -0.19}}) {
        const auto ik = rl::two_link(0.10, 0.12, T, fwd);
        if (!ik) { worst = 1; continue; }
        const double x = 0.10 * std::cos((*ik)[0]) + 0.12 * std::cos((*ik)[1]), z = 0.10 * std::sin((*ik)[0]) + 0.12 * std::sin((*ik)[1]);
        worst = std::max({worst, std::fabs(x - T[0]), std::fabs(z - T[1])});
      }
    s.check("legs: the two-link inverse kinematics reaches its target, either knee (m)", worst, 1e-15);
  }
  {
    const physis::mechanics::PlanarTree<rl::DOF> tree(physis::mechanics::Base::floating, rabbit.links, MOON);
    rl::Controller ctl(rabbit, MOON);
    ctl.attach(tree);
    const auto q = ctl.standing_state(0.0, 0.0);
    const auto st = tree.kinematics(q);
    const auto toe = tree.point(st, rl::FOOT, rabbit.feet[0].point), heel = tree.point(st, rl::FOOT, rabbit.feet[1].point);
    const auto paw = tree.point(st, rl::FOREARM, rabbit.feet[2].point);
    s.check("posture: standing, the toe, the heel and the paw are on the ground, pressed in by their share of the weight (m)",
            std::max({std::fabs(toe[1] + 0.0009), std::fabs(heel[1] + 0.0009), std::fabs(paw[1] + 0.0009)}), 1e-12);
    s.check("posture: the hind toe 0.01 m ahead of the torso's centre (m)", toe[0] - 0.01, 1e-12);
  }
  {
    const auto& m = rabbit.actuators[0];
    const bool line = std::fabs(m.available(100, 0) - m.peak) < 1e-12 && std::fabs(m.available(-100, 0) + m.peak) < 1e-12 &&
                      std::fabs(m.available(100, m.no_load)) < 1e-12 && std::fabs(m.available(100, 0.8 * m.no_load) - 0.2 * m.stall()) < 1e-12 &&
                      std::fabs(m.available(-100, 0.8 * m.no_load) + m.peak) < 1e-12;
    s.require("motor: the torque within the driver's peak and the torque–speed line, none driving at the no-load speed, full braking", line);
    s.check("motor: its electrical power, τω plus the copper loss (τ/k_τ)² R (W)", m.electrical_power(2.0, 10.0) - (20.0 + (2.0 / 0.45) * (2.0 / 0.45) * 0.35), 1e-12);
  }
  {
    rl::Coordinates q = {0.0, 0.15, 0.1, 0.4, -0.1, 3.2, -2.0, 1.0};  // the hip and the ankle beyond their limits
    const auto f = rl::potential_force(rabbit, q);
    double worst = 0;
    for (std::size_t k = 0; k < rl::DOF; k++) {
      auto a = q, b = q;
      a[k] += 1e-6, b[k] -= 1e-6;
      const double d = -(rl::RabbitPotential{&rabbit}(a) - rl::RabbitPotential{&rabbit}(b)) / 2e-6;
      worst = std::max(worst, std::fabs(d - f[k]));
    }
    s.check("potential: the tendon's and the limits' force −∂V/∂q against differences (N m)", worst, 1e-6);
  }

  // ---- standing ----------------------------------------------------------------------------------------------
  {
    rl::SimulationConfig c;
    c.gravity = MOON;
    rl::Simulation sim(rabbit, c);
    sim.advance(1.0);
    const double z1 = sim.tree().centre_of_mass(sim.q())[1];
    sim.advance(1.0);
    const double z2 = sim.tree().centre_of_mass(sim.q())[1];
    double load = 0;
    for (const auto& f : sim.contacts()) load += f.force[1];
    s.check("standing: on the Moon the rabbit is still (the centre's height over a second, m)", z2 - z1, 1e-4);
    s.check("standing: level (rad)", sim.q()[rl::PITCH], 0.01);
    s.check("standing: its feet carry its weight (relative)", load / (rabbit.mass() * MOON) - 1, 0.01);
    s.check("standing: its energy accounted for, the work of its motors, friction and the ground (J)", sim.energy_residual(), 5e-3);
  }

  // ---- hopping ----------------------------------------------------------------------------------------------
  {
    const auto r = hop(rl::Method::semi_implicit, MOON, 1.0, 8.0, 0.0, 0.15);
    s.require("hopping in place: sustained (at least 6 hops in 8 s)", r.hops >= 6, std::to_string(r.hops) + " hops");
    s.check("hopping in place: the pitch held within 15° (rad)", r.max_pitch, 15 * rl::DEG);
    s.check("hopping in place: the apex brought to the 0.15 m asked (m)", r.apex - 0.15, 0.04, "m");
    s.check("hopping in place: the speed near zero (m/s)", r.speed, 0.1, "m/s");
  }
  {
    const auto r = hop(rl::Method::semi_implicit, MOON, 1.0, 10.0, 0.3, 0.15);
    s.check("hopping forward: the speed asked, 0.3 m/s, reached (its mean over the last 5 s, m/s)", r.speed - 0.3, 0.06, "m/s");
    s.check("hopping forward: the pitch held within 15° (rad)", r.max_pitch, 15 * rl::DEG);
  }
  {
    rl::SimulationConfig c;
    c.gravity = MOON;
    rl::Simulation sim(rabbit, c);
    sim.advance(1.0);
    sim.command({rl::Behaviour::hop, 0.2, 0.12});
    sim.advance(6.0);
    sim.command({rl::Behaviour::stop, 0.0, 0.0});
    sim.advance(4.0);
    const bool sits = std::string(rl::phase_name(sim.controller().phase())) == "standing";
    s.require("stopping: from hopping at 0.2 m/s the rabbit lands and sits", sits, rl::phase_name(sim.controller().phase()));
    s.check("stopping: still and level when it sits (rad)", std::max(std::fabs(sim.q()[rl::PITCH]), std::fabs(sim.v()[rl::X])), 0.03);
  }

  // ---- the three methods: one motion ---------------------------------------------------------------------------
  {
    const auto e = hop(rl::Method::semi_implicit, MOON, 1.0, 2.6, 0.0, 0.15);
    const auto k = hop(rl::Method::runge_kutta, MOON, 1.0, 2.6, 0.0, 0.15);
    const auto v = hop(rl::Method::variational, MOON, 1.0, 2.6, 0.0, 0.15);
    s.check("methods: the first hop's apex, Runge–Kutta against the variational integrator (m)", k.apex - v.apex, 0.01, "m");
    s.check("methods: and against semi-implicit Euler (m)", k.apex - e.apex, 0.02, "m");
    s.require("methods: the variational integrator keeps the energy's account best, semi-implicit Euler worst (the residuals)",
              v.residual < k.residual && k.residual < e.residual,
              std::to_string(v.residual) + " < " + std::to_string(k.residual) + " < " + std::to_string(e.residual) + " J");
    s.check("methods: Runge–Kutta keeps the angular momentum about the centre in flight (relative)", k.flight_momentum, 1e-8);
    s.require("methods: each states its cost, the machine's seconds per second of the robot's",
              e.cost > 0 && k.cost > e.cost && v.cost > 0, std::to_string(e.cost) + ", " + std::to_string(k.cost) + ", " + std::to_string(v.cost) + " s/s");
  }

  // ---- the ground it learns ---------------------------------------------------------------------------------------
  {
    // A slippery patch from 0.8 m to 1.6 m (μ 0.3), the presumed regolith elsewhere (μ 0.83); the rabbit is told
    // the presumption (median 0.83, spread 0.185) and hops forward across it.
    rl::SimulationConfig c;
    c.gravity = MOON;
    c.ground_at = [base = rabbit.ground_law](double x) {
      auto law = base;
      law.friction = x > 0.8 && x < 1.6 ? 0.3 : 0.8295;
      return law;
    };
    rl::Simulation sim(rabbit, c);
    sim.advance(1.0);
    sim.command({rl::Behaviour::hop, 0.25, 0.12});
    double max_pitch = 0;
    while (sim.time() < 13.0) {
      sim.advance(0.05);
      max_pitch = std::max(max_pitch, std::fabs(sim.q()[rl::PITCH]));
    }
    s.require("ground: across a slippery patch (μ 0.3) the rabbit keeps hopping (at least 8 hops) and crosses it",
              sim.controller().hops() >= 8 && sim.q()[rl::X] > 1.8, std::to_string(sim.controller().hops()) + " hops, at " + std::to_string(sim.q()[rl::X]) + " m");
    s.check("ground: its pitch held within 15° there (rad)", max_pitch, 15 * rl::DEG);
    std::size_t on_patch = 0, off_patch = 0;
    for (const auto& r : sim.readings())
      if (r.parameter == 0) (r.x > 0.85 && r.x < 1.55 ? on_patch : off_patch)++;
    s.require("ground: its feet read the friction where they slid, on the patch and off it", on_patch > 0 && off_patch > 0,
              std::to_string(on_patch) + " on, " + std::to_string(off_patch) + " off");
    double worst = 0;
    for (const auto& r : sim.readings())
      if (r.parameter == 0 && r.x > 0.85 && r.x < 1.55) worst = std::max(worst, std::fabs(r.value / 0.3 - 1));
    s.check("ground: a sliding foot on the patch reads its friction, 0.3, within its sensor's noise (2 %; the largest relative error, five deviations)", worst, 0.10);
    s.require("ground: and the belief knows its presumption was wrong there: its whitened innovations' mean square far above 1",
              sim.belief().friction.whitened_mean_square() > 10, std::to_string(sim.belief().friction.whitened_mean_square()));
    const auto ahead = sim.belief().friction.at(sim.q()[rl::X] + 5.0);
    s.check("ground: five metres ahead, where no foot has been, it still believes the presumption (relative)", ahead.median / 0.8295 - 1, 1e-3);
  }
  {
    // A ground drawn from the presumption itself (μ and the load at 1 mm log-normal about their medians, varying
    // on a metre's scale): there the belief is calibrated.
    const physis::analysis::StationaryField<1> zm(1.0, 1.0, 5, 1), zl(1.0, 1.0, 5, 2);
    rl::SimulationConfig c;
    c.gravity = MOON;
    c.ground_at = [base = rabbit.ground_law, zm, zl](double x) {
      auto law = base;
      law.friction = 0.8295 * std::exp(0.185 * zm({x}));
      law.stiffness = 9.5 * std::exp(0.5 * zl({x})) / std::pow(1e-3, 1.5);
      return law;
    };
    rl::Simulation sim(rabbit, c);
    sim.advance(1.0);
    sim.command({rl::Behaviour::hop, 0.25, 0.12});
    sim.advance(12.0);
    const auto& b = sim.belief();
    s.require("calibration: on a ground drawn from its presumption, its feet read both the friction and the load",
              b.friction.readings() >= 5 && b.load.readings() >= 5, std::to_string(b.friction.readings()) + " and " + std::to_string(b.load.readings()));
    // The whitened mean square of n readings has mean 1 and deviation √(2/n) when the presumption and the
    // sensors' noise are the truth; the readings are few and correlated by the ground, so a wide band.
    s.check("calibration: its friction readings' whitened mean square is near 1 (|log|)", std::fabs(std::log(b.friction.whitened_mean_square())), std::log(3.0));
    s.check("calibration: and its load readings' (|log|)", std::fabs(std::log(b.load.whitened_mean_square())), std::log(3.0));
    double e = 0;
    int n = 0;
    for (const auto& r : sim.readings())
      if (r.parameter == 0) {
        const auto at = b.friction.at(r.x);
        e += std::pow(std::log(at.median / (0.8295 * std::exp(0.185 * zm({r.x})))), 2), n++;
      }
    s.check("calibration: where its feet slid, it knows the friction within 5 % (the rms of log(belief/truth))", n ? std::sqrt(e / n) : 1.0, 0.05);
  }

  // ---- the Earth ------------------------------------------------------------------------------------------------
  {
    rl::SimulationConfig c;
    c.gravity = 9.81;
    rl::Simulation sim(rabbit, c);
    sim.advance(0.5);
    sim.command({rl::Behaviour::hop, 0.0, 0.05});
    s.require("earth: the hop is refused, the ankle unable to launch it, and says why", !sim.controller().refusal().empty(), sim.controller().refusal());
    sim.advance(1.0);
    double load = 0;
    for (const auto& f : sim.contacts()) load += f.force[1];
    s.check("earth: it stands, its feet carrying its weight (relative)", load / (rabbit.mass() * 9.81) - 1, 0.01);
  }

  return s.done();
}
