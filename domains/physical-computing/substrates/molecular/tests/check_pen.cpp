// Checks of the PEN predator–prey ecosystem (pen.hpp).
//
// - The coexistence equilibrium is a zero of the law; its Hopf point in the prey's template K is K_H = 2x* + 1/(a h),
//   and physis's continuation in K finds the same (to 1e−6).
// - Below K_H the ecosystem settles; above it, it oscillates.
// - The prey alone invades as a pulled front at 2√(D r) less Bramson's correction (within 2%).

#include <physical_computing/molecular/pen.hpp>

#include <physis/dynamics/continuation.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>

using namespace physical_computing;

int main() {
  calx::verification::Suite s;
  molecular::PenPredatorPrey pen;
  {
    auto& L = pen.law;
    physis::dynamics::State x{L.prey_at_coexistence(), L.predator_at_coexistence()}, dx(2);
    L(0.0, x, dx);
    s.check("equilibrium: the coexistence state is a zero of the law", std::hypot(dx[0], dx[1]), 1e-14);
    auto f = [&](const physis::dynamics::State& y, double K, physis::dynamics::State& out) {
      auto law = L;
      law.K = K;
      law(0.0, y, out);
    };
    auto start = L;
    start.K = 1.6;
    physis::dynamics::ContinuationOptions o;
    o.step = 0.01, o.max_step = 0.05, o.mu_max = 4.0, o.max_points = 2000;
    const auto B = physis::dynamics::continue_equilibria(f, {start.prey_at_coexistence(), start.predator_at_coexistence()}, 1.6, o);
    double found = -1;
    for (const auto& b : B.bifurcations)
      if (b.kind == physis::dynamics::Bifurcation::Kind::hopf && found < 0) found = b.mu;
    s.check("Hopf: continuation finds K_H = 2x* + 1/(a h)", found - pen.hopf_template(), 1e-6);
  }
  auto amplitude = [&](double K) {
    auto L = pen.law;
    L.K = K;
    physis::dynamics::State x{L.prey_at_coexistence() * 1.1, L.predator_at_coexistence()};
    physis::dynamics::integrate(L, 0.0, 800.0, x, 0.01);
    double lo = 1e300, hi = -1e300;
    physis::dynamics::integrate(L, 0.0, 200.0, x, 0.01, [&](double, const physis::dynamics::State& y) { lo = std::min(lo, y[0]), hi = std::max(hi, y[0]); });
    return hi - lo;
  };
  s.check("below K_H the ecosystem settles (the prey's range)", amplitude(0.9 * pen.hopf_template()), 1e-6);
  s.require("above K_H it oscillates", amplitude(1.3 * pen.hopf_template()) > 0.1);
  {
    const double D = 1.0, spacing = 0.1;
    auto F = pen.prey_front(2000, spacing, D);
    for (std::size_t i = 0; i < 20; i++) F.state()[i] = pen.law.K;
    auto front = [&]() {
      for (std::size_t i = 0; i < F.vertices(); i++)
        if (F.at(i, 0) < 0.5 * pen.law.K) return i * spacing;
      return F.vertices() * spacing;
    };
    const double dt = 0.2 * F.diffusive_step();
    F.run(30.0, dt);
    const double x1 = front(), t1 = F.time();
    F.run(30.0, dt);
    const double t2 = F.time(), speed = (front() - x1) / (t2 - t1);
    const double r = pen.law.r, bramson = 2 * std::sqrt(D * r) - 1.5 * std::sqrt(D / r) * std::log(t2 / t1) / (t2 - t1);
    s.check("waves: the prey's front at 2√(D r) less Bramson's correction (relative)", speed / bramson - 1, 0.02);
  }
  return s.done();
}
