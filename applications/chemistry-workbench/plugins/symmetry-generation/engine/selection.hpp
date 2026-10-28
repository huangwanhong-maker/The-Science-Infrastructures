// Symmetry selection in the manuscript's §7 and App. I: Frank's autocatalytic rule family as mass
// action (with the enantiomeric excess m and its equation) and as a Markov chain of molecule
// numbers (its exactly symmetric law of outcomes, its generator commuting with the mirror, sample
// paths), the cubic normal form's pitchfork, and the Landau potential along the field's modes.

#pragma once

#include <cmath>
#include <string>
#include <vector>

#include <calx/integrators/runge_kutta.hpp>
#include <calx/random/stream.hpp>
#include <physis/chemistry/kinetics.hpp>

#include "molecules.hpp"

namespace chem {

namespace pc = physis::chemistry;

struct FrankParameters {
  double k = 1, ki = 1, k0 = 0;     // autocatalysis, mutual inhibition, uncatalysed production
  double a0 = 1, l0 = 0.011, r0 = 0.01, t_end = 30;
  long long precursor = 30, seed_l = 1, seed_r = 1;  // the Markov chain's start
};

// A + L → 2L, A + R → 2R (k), L + R → D (k_i), and optionally A → L, A → R (k₀).
inline pc::Network frank_network(const Rational& k, const Rational& ki, const Rational& k0) {
  pc::Network net;
  net.species = {"A", "L", "R", "D"};
  net.add("A + L → 2L", {{"A", 1}, {"L", 1}}, {{"L", 2}}, k);
  net.add("A + R → 2R", {{"A", 1}, {"R", 1}}, {{"R", 2}}, k);
  net.add("L + R → D", {{"L", 1}, {"R", 1}}, {{"D", 1}}, ki);
  if (k0 != Rational(0)) {
    net.add("A → L", {{"A", 1}}, {{"L", 1}}, k0);
    net.add("A → R", {{"A", 1}}, {{"R", 1}}, k0);
  }
  return net;
}

inline Rational rational_of(double x) {
  // Rates are given to at most 4 decimals.
  return Rational(Integer(static_cast<long long>(std::llround(x * 10000))), Integer(10000));
}

inline Json frank(const FrankParameters& P) {
  const auto net = frank_network(rational_of(P.k), rational_of(P.ki), rational_of(P.k0));
  Json out = Json::object({{"reactions", Json::array()}});
  Json rx = Json::array();
  for (const auto& r : net.reactions) rx.push(r.name);
  out.set("reactions", rx);
  // Conservation: a + ℓ + r + 2d.
  Json laws = Json::array();
  for (const auto& y : pc::conservation_laws(net)) {
    Json l = Json::array();
    for (const auto& c : y) l.push(static_cast<double>(c.to_int64()));
    laws.push(l);
  }
  out.set("conservation", laws);
  out.set("mirrorSymmetric", pc::symmetry_of(net, {0, 2, 1, 3}));

  // ---- mass action ----
  using V4 = calx::integrators::Vec<4>;
  V4 c{P.a0, P.l0, P.r0, 0};
  const auto f = [&](const V4& x, V4& dx) { pc::mass_action(net, x, dx); };
  Json t = Json::array(), a = Json::array(), l = Json::array(), r = Json::array(), d = Json::array(), m = Json::array();
  double worst_identity = 0, worst_conservation = 0, sub = 0;
  const double total0 = c[0] + c[1] + c[2] + 2 * c[3];
  const int steps = 600;
  const double h = P.t_end / steps;
  for (int i = 0; i <= steps; i++) {
    const double s = c[1] + c[2];
    t.push(i * h), a.push(c[0]), l.push(c[1]), r.push(c[2]), d.push(c[3]), m.push(s > 0 ? (c[1] - c[2]) / s : 0.0);
    if (s > 1e-12) {
      // ṁ from the vector field, against the manuscript's (k_i s / 2) m (1 − m²), with the term
      // −(2 k₀ a / s) m that uncatalysed production adds.
      V4 dx;
      f(c, dx);
      const double mm = (c[1] - c[2]) / s;
      const double mdot = ((dx[1] - dx[2]) * s - (c[1] - c[2]) * (dx[1] + dx[2])) / (s * s);
      const double formula = P.ki * s / 2 * mm * (1 - mm * mm) - 2 * P.k0 * c[0] / s * mm;
      worst_identity = std::max(worst_identity, std::fabs(mdot - formula) / std::max(1.0, std::fabs(formula)));
    }
    worst_conservation = std::max(worst_conservation, std::fabs(c[0] + c[1] + c[2] + 2 * c[3] - total0));
    if (i < steps) calx::integrators::dopri5<4>(f, c, h, {1, 1, 1, 1}, {1e-12, 1e-10}, sub);
  }
  out.set("ode", Json::object({{"t", t}, {"a", a}, {"l", l}, {"r", r}, {"d", d}, {"m", m}, {"identityError", worst_identity},
                               {"conservationError", worst_conservation}, {"mFinal", (c[1] - c[2]) / (c[1] + c[2])}}));

  // ---- the Markov chain from a symmetric start ----
  const pc::State start{P.precursor, P.seed_l, P.seed_r, 0};
  const auto law = pc::absorption<double>(net, start);
  std::map<long long, double> excess;  // n_L − n_R at absorption
  double extinct = 0, total = 0;
  for (const auto& [state, p] : law) {
    total += p;
    if (state[1] == 0 && state[2] == 0) extinct += p;
    else excess[state[1] - state[2]] += p;
  }
  double asym = 0;
  for (const auto& [k, p] : excess) asym = std::max(asym, std::fabs(p - (excess.count(-k) ? excess.at(-k) : 0.0)));
  Json hist = Json::array();
  for (const auto& [k, p] : excess) hist.push(Json::array({static_cast<double>(k), p}));
  // The generator commutes with the mirror L ↔ R on the reachable states, exactly.
  const auto S = pc::reachable(net, start);
  const auto Q = pc::generator(net, S);
  bool commute = true;
  const pc::SpeciesPermutation mirror{0, 2, 1, 3};
  for (std::size_t x = 0; x < S.states.size() && commute; x++)
    for (const auto& [y, rate] : Q[x]) {
      const std::size_t sx = S.find(pc::permute(mirror, S.states[x])), sy = S.find(pc::permute(mirror, S.states[y]));
      if (sx == S.states.size()) { commute = false; break; }
      const auto it = Q[sx].find(sy);
      commute = commute && it != Q[sx].end() && it->second == rate;
    }
  // Sample paths of m by Gillespie's method.
  Json paths = Json::array();
  calx::random::Stream rng(20261028, 1, 0, 0);
  for (int k = 0; k < 16; k++) {
    std::vector<double> times;
    std::vector<pc::State> states;
    pc::gillespie(net, start, 1e9, rng, times, states);
    Json pt = Json::array(), pm = Json::array();
    for (std::size_t e = 0; e < states.size(); e++) {
      const double s = static_cast<double>(states[e][1] + states[e][2]);
      pt.push(times[e]);
      pm.push(s > 0 ? (states[e][1] - states[e][2]) / s : 0.0);
    }
    paths.push(Json::object({{"t", pt}, {"m", pm}}));
  }
  out.set("chain", Json::object({{"start", Json(std::vector<double>{static_cast<double>(start[0]), static_cast<double>(start[1]), static_cast<double>(start[2]), 0})},
                                 {"states", static_cast<double>(S.states.size())}, {"excess", hist}, {"extinct", extinct}, {"total", total},
                                 {"asymmetry", asym}, {"commutes", commute}, {"paths", paths}}));
  return out;
}

// The exact symmetric law, in rationals, for a small chain (one of the checks).
inline bool frank_law_exactly_symmetric(long long precursor) {
  const auto net = frank_network(1, 1, 0);
  const auto law = pc::absorption<Rational>(net, {precursor, 1, 1, 0});
  for (const auto& [state, p] : law) {
    const pc::State mirror{state[0], state[2], state[1], state[3]};
    const auto it = law.find(mirror);
    if (it == law.end() || it->second != p) return false;
  }
  return true;
}

// ṁ = a m − b m³: equilibria and their stability for a range of a.
inline Json pitchfork(double b) {
  Json branches = Json::array();
  for (int i = 0; i <= 80; i++) {
    const double a = -1 + 2.0 * i / 80;
    Json eq = Json::array();
    eq.push(Json::object({{"m", 0.0}, {"stable", a < 0}}));
    if (a > 0)
      for (double s : {1.0, -1.0}) eq.push(Json::object({{"m", s * std::sqrt(a / b)}, {"stable", a - 3 * b * (a / b) < 0}}));
    branches.push(Json::object({{"a", a}, {"equilibria", eq}}));
  }
  return Json::object({{"b", b}, {"branches", branches}});
}

// The action restricted to t·χ for a character χ of ℤ₂³: (μ/2)·8 t² + λ V[χ] t⁴; V[χ] is 0 along the
// manuscript's (−1)^a and 64 along (−1)^b.
inline Json landau(double mu, double lambda) {
  Json curves = Json::array();
  for (const auto& [name, V] : std::vector<std::pair<std::string, double>>{{"(−1)^a, the manuscript's η", 0}, {"(−1)^b", 64}, {"(−1)^{a+b+c}", 64}}) {
    Json t = Json::array(), S = Json::array();
    for (int i = 0; i <= 120; i++) {
      const double x = -1.5 + 3.0 * i / 120;
      t.push(x), S.push(mu / 2 * 8 * x * x + lambda * V * x * x * x * x);
    }
    double minimum = 0;
    if (V > 0 && mu < 0) minimum = std::sqrt(-4 * mu / (2 * 64 * lambda));  // α = 4μ, β = 64λ: η* = √(−α/2β)
    curves.push(Json::object({{"mode", name}, {"V", V}, {"t", t}, {"S", S}, {"minimum", minimum}}));
  }
  return Json::object({{"mu", mu}, {"lambda", lambda}, {"curves", curves}});
}

}  // namespace chem
