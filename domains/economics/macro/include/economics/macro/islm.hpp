// The IS–LM model of a closed economy in the short run (Hicks 1937), in its linear textbook form
// (Mankiw, Macroeconomics, ch. 11–12; Blanchard, Macroeconomics, ch. 5):
//
//   goods (IS):   Y = C + I + G,  C = c₀ + c₁ (Y − T),  I = i₀ − i₁ r
//                 ⇒  r = (A − (1 − c₁) Y) / i₁,   A = c₀ − c₁ T + i₀ + G   (autonomous spending)
//   money (LM):   M / P = k Y − h r
//                 ⇒  r = (k Y − M/P) / h
//
// with 0 < c₁ < 1, i₁, k, h > 0. Their crossing is the equilibrium,
//
//   Y* = (h A + i₁ M/P) / D,   r* = (k A − (1 − c₁) M/P) / D,   D = h (1 − c₁) + i₁ k,
//
// and its comparative statics the multipliers dY*/dG = h / D, dY*/dT = −c₁ h / D,
// dY*/d(M/P) = i₁ / D. The Keynesian cross is the limit of a fixed rate (h → ∞): dY/dG = 1/(1 − c₁);
// the classical one, inelastic money demand (h → 0): dY/dG = 0. Letting the price level vary traces
// the aggregate demand curve, Y*(P) = (h A + i₁ M / P) / D.

#pragma once

#include <stdexcept>
#include <vector>

namespace economics::macro {

struct ISLM {
  double c0 = 200, c1 = 0.75;   // autonomous consumption, the marginal propensity to consume
  double taxes = 100;           // T
  double i0 = 200, i1 = 25;     // autonomous investment; its fall per point of the interest rate
  double government = 100;      // G
  double money = 1000, price = 2;  // M and P
  double k = 0.5, h = 50;       // money demand: k Y − h r

  void check() const {
    if (!(c1 > 0 && c1 < 1)) throw std::invalid_argument("economics::macro::ISLM: the marginal propensity to consume lies in (0, 1)");
    if (!(i1 > 0 && k > 0 && h > 0 && price > 0)) throw std::invalid_argument("economics::macro::ISLM: i₁, k, h and P are positive");
  }
  double autonomous() const { return c0 - c1 * taxes + i0 + government; }
  double real_money() const { return money / price; }
  double is_rate(double Y) const { return (autonomous() - (1 - c1) * Y) / i1; }   // r on the IS curve
  double lm_rate(double Y) const { return (k * Y - real_money()) / h; }           // r on the LM curve
};

struct ISLMEquilibrium {
  double output = 0, rate = 0, consumption = 0, investment = 0;
  double fiscal_multiplier = 0, tax_multiplier = 0, money_multiplier = 0, cross_multiplier = 0;
};

inline ISLMEquilibrium equilibrium(const ISLM& m) {
  m.check();
  const double D = m.h * (1 - m.c1) + m.i1 * m.k, A = m.autonomous(), Mp = m.real_money();
  ISLMEquilibrium E;
  E.output = (m.h * A + m.i1 * Mp) / D;
  E.rate = (m.k * A - (1 - m.c1) * Mp) / D;
  E.consumption = m.c0 + m.c1 * (E.output - m.taxes);
  E.investment = m.i0 - m.i1 * E.rate;
  E.fiscal_multiplier = m.h / D;
  E.tax_multiplier = -m.c1 * m.h / D;
  E.money_multiplier = m.i1 / D;
  E.cross_multiplier = 1 / (1 - m.c1);
  return E;
}

struct ISLMCurves {
  std::vector<double> output, is, lm;     // r on each curve at each Y
  std::vector<double> price, demand;      // the aggregate demand curve Y*(P)
};

// The two curves over [Y_lo, Y_hi], and aggregate demand over prices [P_lo, P_hi].
inline ISLMCurves curves(const ISLM& m, double y_lo, double y_hi, double p_lo, double p_hi, std::size_t points = 60) {
  m.check();
  if (!(y_hi > y_lo && p_hi > p_lo && p_lo > 0 && points > 1)) throw std::invalid_argument("economics::macro::curves: increasing ranges, positive prices");
  ISLMCurves C;
  const double D = m.h * (1 - m.c1) + m.i1 * m.k, A = m.autonomous();
  for (std::size_t i = 0; i < points; i++) {
    const double u = static_cast<double>(i) / static_cast<double>(points - 1);
    const double Y = y_lo + u * (y_hi - y_lo), P = p_lo + u * (p_hi - p_lo);
    C.output.push_back(Y);
    C.is.push_back(m.is_rate(Y));
    C.lm.push_back(m.lm_rate(Y));
    C.price.push_back(P);
    C.demand.push_back((m.h * A + m.i1 * m.money / P) / D);
  }
  return C;
}

}  // namespace economics::macro
