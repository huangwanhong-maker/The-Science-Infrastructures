// The Solow–Swan model of growth (Solow 1956; Swan 1956; Barro and Sala-i-Martin 2004, ch. 1).
//
// Capital per effective worker k = K / (A L) obeys
//
//   dk/dt = s f(k) − (n + g + δ) k,                                           (B&S 1.13)
//
// s the saving rate, n the growth of labour, g of technology, δ depreciation. The technology is
// Cobb–Douglas, f(k) = k^α, or of constant elasticity of substitution σ (CES),
// f(k) = [α k^ρ + (1 − α)]^(1/ρ), ρ = (σ − 1)/σ (Arrow, Chenery, Minhas and Solow 1961; σ = 1 is the
// Cobb–Douglas limit).
//
//   the steady state k*:  s f(k*) = (n + g + δ) k*; for Cobb–Douglas k* = (s / (n + g + δ))^(1/(1−α));
//   the golden rule:      f'(k_g) = n + g + δ, the k that makes consumption c = f(k) − (n+g+δ)k largest;
//                         for Cobb–Douglas s_g = α (Phelps 1961);
//   the speed of convergence near k*: λ = (n + g + δ)(1 − ε(k*)), ε = k f'(k) / f(k) the capital share
//                         (B&S 1.45); for Cobb–Douglas λ = (1 − α)(n + g + δ), half-life ln 2 / λ.
//
// The path from k(0): for Cobb–Douglas in closed form, since k^(1−α) is linear in its own equation,
//
//   k(t)^(1−α) = k*^(1−α) + (k(0)^(1−α) − k*^(1−α)) e^(−(1−α)(n+g+δ) t);       (B&S 1.A)
//
// for CES by calx's fourth-order Runge–Kutta, steps of a hundredth of a year by default.
// Output per worker grows at g in the steady state; the path reports output, consumption and
// investment per effective worker, and the growth rate of output per effective worker.

#pragma once

#include <calx/integrators/runge_kutta.hpp>
#include <calx/solvers/roots.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace economics::growth {

struct Solow {
  double saving = 0.2;           // s
  double labour_growth = 0.01;   // n, a year
  double technology_growth = 0.02;  // g
  double depreciation = 0.05;    // δ
  double alpha = 1.0 / 3.0;      // the capital share parameter
  double sigma = 1.0;            // the elasticity of substitution; 1 is Cobb–Douglas

  double effective_depreciation() const { return labour_growth + technology_growth + depreciation; }
  bool cobb_douglas() const { return std::fabs(sigma - 1) < 1e-12; }
  void check() const {
    if (!(saving > 0 && saving < 1)) throw std::invalid_argument("economics::growth::Solow: the saving rate lies in (0, 1)");
    if (!(alpha > 0 && alpha < 1)) throw std::invalid_argument("economics::growth::Solow: α lies in (0, 1)");
    if (!(effective_depreciation() > 0)) throw std::invalid_argument("economics::growth::Solow: n + g + δ is positive");
    if (!(sigma > 0)) throw std::invalid_argument("economics::growth::Solow: σ is positive");
  }
  // Output per effective worker and its derivative.
  double f(double k) const {
    if (cobb_douglas()) return std::pow(k, alpha);
    const double rho = (sigma - 1) / sigma;
    return std::pow(alpha * std::pow(k, rho) + (1 - alpha), 1 / rho);
  }
  double f_prime(double k) const {
    if (cobb_douglas()) return alpha * std::pow(k, alpha - 1);
    const double rho = (sigma - 1) / sigma;
    return alpha * std::pow(k, rho - 1) * std::pow(alpha * std::pow(k, rho) + (1 - alpha), 1 / rho - 1);
  }
  double motion(double k) const { return saving * f(k) - effective_depreciation() * k; }
};

struct SteadyState {
  double k = 0, y = 0, c = 0, investment = 0;
  double capital_share = 0;      // ε(k*) = k f'(k) / f(k)
  double convergence = 0;        // λ, a year
  double half_life = 0;          // ln 2 / λ, years
  double golden_k = 0, golden_saving = 0, golden_c = 0;
};

namespace detail {
// A bracket [a, b] with a sign change of g, widened geometrically from a guess.
template <class G>
inline double root(G&& g, double guess) {
  double a = guess, b = guess;
  for (int i = 0; i < 400; i++) {
    a *= 0.5;
    b *= 2;
    if ((g(a) > 0) != (g(b) > 0)) return calx::solvers::brent(g, a, b, 1e-14).x;
  }
  throw std::domain_error("economics::growth: no steady state (the saving cannot sustain the effective depreciation, or always outruns it)");
}
}  // namespace detail

inline SteadyState steady_state(const Solow& m) {
  m.check();
  const double d = m.effective_depreciation();
  SteadyState S;
  if (m.cobb_douglas()) {
    S.k = std::pow(m.saving / d, 1 / (1 - m.alpha));
    S.golden_k = std::pow(m.alpha / d, 1 / (1 - m.alpha));
  } else {
    S.k = detail::root([&](double k) { return m.saving * m.f(k) / k - d; }, 1.0);
    S.golden_k = detail::root([&](double k) { return m.f_prime(k) - d; }, 1.0);
  }
  S.y = m.f(S.k);
  S.investment = m.saving * S.y;
  S.c = S.y - S.investment;
  S.capital_share = S.k * m.f_prime(S.k) / S.y;
  S.convergence = d * (1 - S.capital_share);
  S.half_life = std::log(2.0) / S.convergence;
  S.golden_c = m.f(S.golden_k) - d * S.golden_k;
  S.golden_saving = d * S.golden_k / m.f(S.golden_k);
  return S;
}

struct SolowPath {
  std::vector<double> t, k, y, c, investment, growth;  // growth: d ln y / dt, a year
};

// The path from k0 over the years [0, T], sampled `samples` + 1 times.
inline SolowPath path(const Solow& m, double k0, double T, std::size_t samples = 200, double step = 0.01) {
  m.check();
  if (!(k0 > 0 && T > 0 && samples > 0)) throw std::invalid_argument("economics::growth::path: k0 > 0, T > 0 and at least one sample");
  const double d = m.effective_depreciation();
  SolowPath P;
  auto record = [&](double t, double k) {
    const double y = m.f(k);
    P.t.push_back(t);
    P.k.push_back(k);
    P.y.push_back(y);
    P.investment.push_back(m.saving * y);
    P.c.push_back((1 - m.saving) * y);
    P.growth.push_back(m.f_prime(k) * m.motion(k) / y);  // d ln y/dt = f'(k) k̇ / f(k)
  };
  if (m.cobb_douglas()) {
    const double ks = std::pow(m.saving / d, 1 / (1 - m.alpha)), a = 1 - m.alpha;
    for (std::size_t i = 0; i <= samples; i++) {
      const double t = T * static_cast<double>(i) / static_cast<double>(samples);
      record(t, std::pow(std::pow(ks, a) + (std::pow(k0, a) - std::pow(ks, a)) * std::exp(-a * d * t), 1 / a));
    }
    return P;
  }
  calx::integrators::Vec<1> x{k0};
  const auto rhs = [&](const calx::integrators::Vec<1>& s, calx::integrators::Vec<1>& dk) { dk[0] = m.motion(s[0]); };
  double t = 0;
  record(0, k0);
  for (std::size_t i = 1; i <= samples; i++) {
    const double target = T * static_cast<double>(i) / static_cast<double>(samples);
    while (t < target - 1e-12) {
      const double h = std::min(step, target - t);
      calx::integrators::rk4<1>(rhs, x, h);
      t += h;
    }
    record(target, x[0]);
  }
  return P;
}

}  // namespace economics::growth
