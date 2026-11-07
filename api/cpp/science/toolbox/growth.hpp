// The toolbox's growth and macroeconomics: the economics domain's Solow model, IS-LM and growth
// accounting.

#pragma once

#include <economics/growth/accounting.hpp>
#include <economics/growth/solow.hpp>
#include <economics/macro/islm.hpp>

#include "answer.hpp"

namespace science::toolbox {

// The Solow–Swan model: its steady state, golden rule and path from a share of the steady capital.
inline Json solow(const Request& r) {
  economics::growth::Solow m;
  m.saving = r.number("saving", 0.2);
  m.labour_growth = r.number("labour_growth", 0.01);
  m.technology_growth = r.number("technology_growth", 0.02);
  m.depreciation = r.number("depreciation", 0.05);
  m.alpha = r.number("alpha", 1.0 / 3);
  m.sigma = r.number("sigma", 1.0);
  const auto S = economics::growth::steady_state(m);
  const double years = r.number("years", 80), start = r.number("start", 0.5);
  const auto P = economics::growth::path(m, start * S.k, years, 160);
  // The same economy saving at the golden rule's rate, for comparison.
  economics::growth::Solow g = m;
  g.saving = std::min(0.99, std::max(0.01, S.golden_saving));
  const auto PG = economics::growth::path(g, start * S.k, years, 160);
  Arr blocks;
  blocks.push_back(figures("The steady state, per effective worker", {
      {"capital k*", S.k}, {"output y*", S.y}, {"consumption c*", S.c}, {"capital share at k*", S.capital_share},
      {"speed of convergence λ, a year", S.convergence}, {"half-life, years", S.half_life},
      {"golden-rule saving", S.golden_saving}, {"golden-rule consumption", S.golden_c}}));
  blocks.push_back(lines("Capital per effective worker", P.t, "years", "k", {{"saving at s", P.k, {}, {}, false}, {"at the golden rule's s", PG.k, {}, {}, true}}));
  blocks.push_back(lines("Output and consumption per effective worker", P.t, "years", "", {{"output y", P.y, {}, {}, false}, {"consumption c", P.c, {}, {}, false}}));
  std::vector<double> growth_pct;
  for (double x : P.growth) growth_pct.push_back(100 * x);
  blocks.push_back(lines("Growth of output per effective worker, % a year (per worker: add g)", P.t, "years", "%", {{"d ln y / dt", growth_pct, {}, {}, false}}));
  blocks.push_back(note(m.cobb_douglas() ? "Cobb–Douglas: the path in closed form (Barro and Sala-i-Martin 2004, ch. 1)."
                                         : "CES: the path by calx's fourth-order Runge–Kutta, the steady state by Brent's method."));
  return answer("The Solow–Swan model", std::move(blocks));
}

// IS–LM: the equilibrium, the multipliers, and the two curves with a policy's shift.
inline Json islm(const Request& r) {
  economics::macro::ISLM m;
  m.c0 = r.number("c0", m.c0);
  m.c1 = r.number("c1", m.c1);
  m.taxes = r.number("taxes", m.taxes);
  m.i0 = r.number("i0", m.i0);
  m.i1 = r.number("i1", m.i1);
  m.government = r.number("government", m.government);
  m.money = r.number("money", m.money);
  m.price = r.number("price", m.price);
  m.k = r.number("k", m.k);
  m.h = r.number("h", m.h);
  const auto E = economics::macro::equilibrium(m);
  economics::macro::ISLM shifted = m;
  shifted.government += r.number("delta_government", 0);
  shifted.money += r.number("delta_money", 0);
  shifted.taxes += r.number("delta_taxes", 0);
  const auto E2 = economics::macro::equilibrium(shifted);
  const double lo = 0.5 * std::min(E.output, E2.output), hi = 1.5 * std::max(E.output, E2.output);
  const auto C = economics::macro::curves(m, lo, hi, 0.5 * m.price, 2 * m.price, 60);
  const auto C2 = economics::macro::curves(shifted, lo, hi, 0.5 * m.price, 2 * m.price, 60);
  Arr blocks;
  blocks.push_back(figures("The equilibrium", {{"output Y*", E.output}, {"interest rate r*", E.rate}, {"consumption", E.consumption}, {"investment", E.investment},
                                               {"fiscal multiplier dY/dG", E.fiscal_multiplier}, {"tax multiplier dY/dT", E.tax_multiplier},
                                               {"money multiplier dY/d(M/P)", E.money_multiplier}, {"Keynesian cross 1/(1−c₁)", E.cross_multiplier}}));
  if (E2.output != E.output || E2.rate != E.rate)
    blocks.push_back(figures("After the policy", {{"output Y*", E2.output}, {"interest rate r*", E2.rate}, {"change in output", E2.output - E.output},
                                                  {"change in the rate", E2.rate - E.rate}}));
  // A curve after the policy only where the policy moves it: fiscal policy the IS, money the LM.
  const bool fiscal = shifted.government != m.government || shifted.taxes != m.taxes, monetary = shifted.money != m.money;
  std::vector<Series> isl = {{"IS", C.is, {}, {}, false}, {"LM", C.lm, {}, {}, false}};
  if (fiscal) isl.push_back({"IS after", C2.is, {}, {}, true});
  if (monetary) isl.push_back({"LM after", C2.lm, {}, {}, true});
  blocks.push_back(lines("The IS and LM curves", C.output, "output Y", "interest rate r", isl));
  std::vector<Series> ad = {{"AD", C.demand, {}, {}, false}};
  if (fiscal || monetary) ad.push_back({"AD after", C2.demand, {}, {}, true});
  blocks.push_back(lines("Aggregate demand: output at each price level", C.price, "price level P", "output Y", ad));
  return answer("IS–LM", std::move(blocks));
}

// Growth accounting from a table of output, capital (or investment) and labour over time.
inline Json growth_accounting(const Request& r) {
  const auto time = r.numbers(r.one("time"));
  std::vector<double> Y = r.numbers(r.one("output")), L = r.numbers(r.one("labour"));
  const bool from_investment = !r.role("investment").empty() && r.role("capital").empty();
  std::vector<double> K = from_investment ? r.numbers(r.one("investment")) : r.numbers(r.one("capital"));
  auto C = complete({time, Y, K, L});
  if (C.rows.size() < 3) throw std::invalid_argument("growth accounting needs at least three complete periods");
  // In the order of time.
  std::vector<std::size_t> order(C.rows.size());
  for (std::size_t i = 0; i < order.size(); i++) order[i] = i;
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return C.columns[0][a] < C.columns[0][b]; });
  std::vector<double> t, y, k, l;
  for (std::size_t i : order) {
    t.push_back(C.columns[0][i]);
    y.push_back(C.columns[1][i]);
    k.push_back(C.columns[2][i]);
    l.push_back(C.columns[3][i]);
  }
  if (from_investment) k = economics::growth::perpetual_inventory(k, r.number("depreciation", 0.05));
  const auto G = economics::growth::growth_accounts(y, k, l, {r.number("capital_share", 1.0 / 3)});
  std::vector<std::string> labels;
  std::vector<double> cap, lab, tfp;
  Arr rows;
  for (std::size_t i = 0; i < G.output.size(); i++) {
    labels.push_back(Json(t[i + 1]).dump());
    cap.push_back(100 * G.capital[i]);
    lab.push_back(100 * G.labour[i]);
    tfp.push_back(100 * G.tfp[i]);
    rows.push_back(Json(Arr{Json(t[i + 1]), num(100 * G.output[i]), num(100 * G.capital[i]), num(100 * G.labour[i]), num(100 * G.tfp[i])}));
  }
  Arr blocks;
  blocks.push_back(figures("A period, on average (log points × 100)", {{"output", 100 * G.mean_output}, {"capital", 100 * G.mean_capital},
                                                                       {"labour", 100 * G.mean_labour}, {"total factor productivity", 100 * G.mean_tfp},
                                                                       {"TFP's share of output growth", G.tfp_share}}));
  blocks.push_back(bars("The growth of output, period by period, by its sources (log points × 100)", labels,
                        {{"capital", cap}, {"labour", lab}, {"TFP (the residual)", tfp}}, true));
  blocks.push_back(table("Growth accounts", {"period", "output", "capital", "labour", "TFP"}, std::move(rows),
                         from_investment ? "Capital by the perpetual inventory of the investment column." : ""));
  return answer("Growth accounting", std::move(blocks));
}

}  // namespace science::toolbox
