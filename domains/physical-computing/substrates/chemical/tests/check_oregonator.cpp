// Checks of the Oregonator and the BZ medium (oregonator.hpp).
//
// - The steady state solves the kinetics (to rounding), with light and without.
// - Its Hopf points in f, found where the trace vanishes, agree with physis's continuation in f (to 1e−6), which
//   finds them as a complex pair crossing.
// - Between them the stirred medium oscillates; beyond the upper one it rests.
// - Beyond the upper one it is excitable: on a chain whose lattice resolves the front, a strong kick runs a pulse 25
//   units without decay and a weak one dies; the pulse's speed converges as the lattice refines (within 5%).
// - Light quenches the oscillation: continuation in φ at an oscillating f finds a Hopf point at some φ > 0.

#include <physical_computing/chemical/oregonator.hpp>

#include <physis/dynamics/continuation.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <tuple>

using namespace physical_computing;

int main() {
  calx::verification::Suite s;
  chemical::Oregonator o;
  for (double light : {0.0, 0.01}) {
    o.light = light;
    const double u = o.steady_u();
    double c[2] = {u, u}, dc[2];
    o.kinetics(c, dc);
    s.check(light == 0 ? "steady state: solves the kinetics, no light" : "steady state: solves the kinetics, with light", std::fabs(dc[0]) * o.epsilon + std::fabs(dc[1]), 1e-13);
  }
  o.light = 0;
  const auto hopf = chemical::hopf_points_in_f(o);
  s.require("Hopf: the trace vanishes at two values of f", hopf.size() == 2);
  if (hopf.size() == 2) {
    auto f = [&](const physis::dynamics::State& x, double ff, physis::dynamics::State& out) {
      chemical::Oregonator p = o;
      p.f = ff;
      double dc[2];
      p.kinetics(x.data(), dc);
      out[0] = dc[0], out[1] = dc[1];
    };
    chemical::Oregonator start = o;
    start.f = 0.5 * hopf[0];
    physis::dynamics::ContinuationOptions co;
    co.step = 0.005, co.max_step = 0.02, co.mu_max = hopf[1] + 0.3, co.max_points = 4000;
    const auto B = physis::dynamics::continue_equilibria(f, {start.steady_u(), start.steady_u()}, start.f, co);
    std::vector<double> found;
    for (const auto& b : B.bifurcations)
      if (b.kind == physis::dynamics::Bifurcation::Kind::hopf) found.push_back(b.mu);
    s.require("Hopf: continuation finds the same two", found.size() == 2);
    if (found.size() == 2) {
      s.check("Hopf: the lower point, by the trace and by continuation", found[0] - hopf[0], 1e-6);
      s.check("Hopf: the upper point, by the trace and by continuation", found[1] - hopf[1], 1e-6);
    }
    auto amplitude = [&](double ff) {
      chemical::Oregonator p = o;
      p.f = ff;
      physis::dynamics::State x{p.steady_u() * 1.2, p.steady_u()};
      physis::dynamics::integrate(p, 0.0, 200.0, x, 1e-3);
      double lo = 1e300, hi = -1e300;
      physis::dynamics::integrate(p, 0.0, 50.0, x, 1e-3, [&](double, const physis::dynamics::State& y) { lo = std::min(lo, y[1]), hi = std::max(hi, y[1]); });
      return hi - lo;
    };
    s.require("between the Hopf points the stirred medium oscillates", amplitude(0.5 * (hopf[0] + hopf[1])) > 0.01);
    s.check("beyond the upper one it rests (the catalyst's range)", amplitude(hopf[1] + 0.4), 1e-6);

    // Excitability on a chain beyond the upper Hopf point, its fronts resolved (spacing 0.1 against a front about
    // √(εD) = 0.2 wide; on a lattice coarser than the front, propagation fails: the lattice pins it).
    chemical::Oregonator ex = o;
    ex.f = hopf[1] + 0.3;
    auto pulse = [&](double kick, double spacing, double& arrival) {
      const std::size_t n = static_cast<std::size_t>(30 / spacing);
      auto m = chemical::BZMedium::chain(ex, n, spacing);
      for (std::size_t i = 0; i < static_cast<std::size_t>(1.0 / spacing); i++) m.excite(i, 0, kick);
      const std::size_t p10 = static_cast<std::size_t>(10 / spacing), p25 = static_cast<std::size_t>(25 / spacing);
      double far = 0, t10 = -1;
      arrival = -1;
      const double h = m.stable_step();
      for (int k = 0; k < 120; k++) {
        m.run(0.05, h);
        const auto u = m.activator();
        far = std::max(far, u[p25]);
        if (t10 < 0 && u[p10] > 0.2) t10 = m.time();
        if (arrival < 0 && u[p25] > 0.2) arrival = m.time();
      }
      if (t10 > 0 && arrival > 0) arrival = 15 / (arrival - t10);  // the speed between 10 and 25
      return far;
    };
    double speed_coarse = 0, speed_fine = 0, dummy = 0;
    s.require("excitable: a strong kick runs a pulse 25 units without decay", pulse(0.8, 0.1, speed_coarse) > 0.5);
    s.require("excitable: a weak kick dies", pulse(ex.steady_u() * 1.05, 0.1, dummy) < 0.01);
    pulse(0.8, 0.05, speed_fine);
    s.check("excitable: the pulse's speed converges as the lattice refines (spacing 0.1 against 0.05, relative)", speed_coarse / speed_fine - 1, 0.05);

    // Light quenches the oscillation.
    chemical::Oregonator osc = o;
    osc.f = 0.5 * (hopf[0] + hopf[1]);
    auto g = [&](const physis::dynamics::State& x, double phi, physis::dynamics::State& out) {
      chemical::Oregonator p = osc;
      p.light = phi;
      double dc[2];
      p.kinetics(x.data(), dc);
      out[0] = dc[0], out[1] = dc[1];
    };
    physis::dynamics::ContinuationOptions lo;
    lo.step = 1e-4, lo.max_step = 2e-3, lo.mu_max = 0.2, lo.max_points = 4000;
    const auto L = physis::dynamics::continue_equilibria(g, {osc.steady_u(), osc.steady_u()}, 0.0, lo);
    double phi_c = -1;
    for (const auto& b : L.bifurcations)
      if (b.kind == physis::dynamics::Bifurcation::Kind::hopf && phi_c < 0) phi_c = b.mu;
    s.require("light: a Hopf point at some φ > 0 quenches the oscillation", phi_c > 0, "φ_c = " + std::to_string(phi_c));
  }
  return s.done();
}
