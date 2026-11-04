// Checks of growth (solow.hpp, accounting.hpp) and macro (islm.hpp).
//
// Solow: the Cobb–Douglas steady state k* = (s/(n+g+δ))^(1/(1−α)) holds its equation of motion, the
// golden rule's saving is α and gives more consumption than any other saving rate, the closed-form
// path agrees with calx's Runge–Kutta, the gap ln k − ln k* halves in about ln 2 / λ near k*, and a
// CES technology with σ → 1 tends to the Cobb–Douglas answers while σ ≠ 1 holds its own steady state.
// Growth accounting: on data made from a known A, K and L the residual is A's growth, the
// contributions add up to output's growth exactly, and the perpetual inventory reproduces a stock
// built by its own recursion. IS–LM: the equilibrium lies on both curves, the multipliers equal
// finite differences, and the Keynesian cross and the classical case are the limits h → ∞ and h → 0.

#include <economics/growth/accounting.hpp>
#include <economics/growth/solow.hpp>
#include <economics/macro/islm.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace economics;

int main() {
  calx::verification::Suite s;

  // ---- Solow ---------------------------------------------------------------------------------------
  growth::Solow m;  // s 0.2, n 0.01, g 0.02, δ 0.05, α 1/3
  const auto S = growth::steady_state(m);
  const double d = m.effective_depreciation();
  s.check("Solow: k* = (s/(n+g+δ))^(1/(1−α)) holds s f(k*) = (n+g+δ) k*", (m.saving * m.f(S.k) - d * S.k) / S.k, 1e-14);
  s.check("... k* = 2.5^1.5 for s = 0.2, n+g+δ = 0.08, α = 1/3", S.k - std::pow(2.5, 1.5), 1e-12);
  s.check("... its speed of convergence (1 − α)(n+g+δ)", S.convergence - (1 - m.alpha) * d, 1e-12);
  s.check("Solow: the golden rule's saving is α", S.golden_saving - m.alpha, 1e-12);
  double best = 0;
  for (double sv = 0.01; sv < 0.99; sv += 0.01) {
    growth::Solow q = m;
    q.saving = sv;
    best = std::max(best, growth::steady_state(q).c);
  }
  s.require("... and no saving rate on a grid of 0.01 gives more steady consumption", best <= S.golden_c + 1e-12);
  const auto closed = growth::path(m, 0.5 * S.k, 60, 60);
  double worst = 0;
  {
    calx::integrators::Vec<1> x{0.5 * S.k};
    const auto rhs = [&](const calx::integrators::Vec<1>& q, calx::integrators::Vec<1>& dk) { dk[0] = m.motion(q[0]); };
    for (std::size_t i = 1; i <= 60; i++) {
      for (int j = 0; j < 100; j++) calx::integrators::rk4<1>(rhs, x, 0.01);
      worst = std::max(worst, std::fabs(x[0] - closed.k[i]) / closed.k[i]);
    }
  }
  s.check("Solow: the closed-form path agrees with calx's Runge–Kutta (relative)", worst, 1e-10);
  // Near k*, the log gap halves in about the half-life.
  const auto near = growth::path(m, 0.98 * S.k, 4 * S.half_life, 4000);
  std::size_t half = 0;
  const double gap0 = std::log(S.k / near.k[0]);
  while (half < near.k.size() && std::log(S.k / near.k[half]) > 0.5 * gap0) half++;
  s.check("Solow: near k*, the log gap halves in ln 2 / λ (relative, the linearisation's error)", (near.t[half] - S.half_life) / S.half_life, 0.02);
  growth::Solow almost = m;
  almost.sigma = 1.0001;
  const auto SA = growth::steady_state(almost);
  s.check("Solow: CES with σ = 1.0001 is near Cobb–Douglas (k*, relative)", (SA.k - S.k) / S.k, 1e-3);
  growth::Solow sub = m;
  sub.sigma = 0.5;
  const auto SS = growth::steady_state(sub);
  s.check("Solow: CES with σ = 0.5 holds its steady state", (sub.saving * sub.f(SS.k) - d * SS.k) / SS.k, 1e-12);
  const auto P = growth::path(sub, 0.3 * SS.k, 200, 200);
  s.check("... and its Runge–Kutta path reaches it (relative)", (P.k.back() - SS.k) / SS.k, 1e-3);

  // ---- growth accounting ------------------------------------------------------------------------------
  std::vector<double> Y, K, L, A;
  for (int t = 0; t < 12; t++) {
    A.push_back(std::exp(0.013 * t + 0.002 * std::sin(t)));
    K.push_back(100 * std::exp(0.04 * t));
    L.push_back(50 * std::exp(0.008 * t + 0.003 * std::cos(t)));
    Y.push_back(A.back() * std::pow(K.back(), 0.35) * std::pow(L.back(), 0.65));
  }
  const auto G = growth::growth_accounts(Y, K, L, {0.35});
  double residual = 0, adds = 0;
  for (std::size_t i = 0; i < G.tfp.size(); i++) {
    residual = std::max(residual, std::fabs(G.tfp[i] - std::log(A[i + 1] / A[i])));
    adds = std::max(adds, std::fabs(G.capital[i] + G.labour[i] + G.tfp[i] - G.output[i]));
  }
  s.check("growth accounting: the residual is A's growth when Y = A K^0.35 L^0.65", residual, 1e-12);
  s.check("... and the contributions add up to output's growth", adds, 1e-15);
  std::vector<double> I{10, 11, 12.5, 13, 14, 15.5, 16};
  const auto Kp = growth::perpetual_inventory(I, 0.06, 0.05);
  double rec = std::fabs(Kp[0] - 10 / 0.11);
  for (std::size_t t = 1; t < I.size(); t++) rec = std::max(rec, std::fabs(Kp[t] - (0.94 * Kp[t - 1] + I[t - 1])));
  s.check("perpetual inventory: K₀ = I₀/(g+δ) and K_t = (1−δ)K_{t−1} + I_{t−1}", rec, 1e-12);

  // ---- IS–LM -------------------------------------------------------------------------------------------
  macro::ISLM is;
  const auto E = macro::equilibrium(is);
  s.check("IS–LM: the equilibrium lies on the IS curve", E.rate - is.is_rate(E.output), 1e-12);
  s.check("... and on the LM curve", E.rate - is.lm_rate(E.output), 1e-12);
  s.check("... goods market: Y = C + I + G", E.output - (E.consumption + E.investment + is.government), 1e-9);
  auto output_with = [&](auto change) { macro::ISLM q = is; change(q); return macro::equilibrium(q).output; };
  const double h = 1e-3;
  s.check("IS–LM: the fiscal multiplier is dY/dG (finite difference)", E.fiscal_multiplier - (output_with([&](macro::ISLM& q) { q.government += h; }) - E.output) / h, 1e-6);
  s.check("... the tax multiplier dY/dT", E.tax_multiplier - (output_with([&](macro::ISLM& q) { q.taxes += h; }) - E.output) / h, 1e-6);
  s.check("... the money multiplier dY/d(M/P)", E.money_multiplier - (output_with([&](macro::ISLM& q) { q.money += h * q.price; }) - E.output) / h, 1e-6);
  macro::ISLM trap = is;
  trap.h = 1e12;
  s.check("IS–LM: with money demand infinitely elastic, the fiscal multiplier is the cross's 1/(1−c₁)", macro::equilibrium(trap).fiscal_multiplier - 4, 1e-9);
  macro::ISLM classical = is;
  classical.h = 1e-12;
  s.check("... and with it inelastic, fiscal policy moves nothing", macro::equilibrium(classical).fiscal_multiplier, 1e-9);
  const auto C = macro::curves(is, 500, 2000, 1, 4, 31);
  s.require("IS–LM: the IS curve falls and the LM curve rises; aggregate demand falls in P",
            C.is.front() > C.is.back() && C.lm.front() < C.lm.back() && C.demand.front() > C.demand.back());
  return s.done();
}
