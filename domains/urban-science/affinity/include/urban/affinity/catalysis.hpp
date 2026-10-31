// Catalysis in space and time: the bridal store and the ring store (theory/street-oceanography.md,
// Representation 5.2; Huang, Street Oceanography, §11.3 third mechanism and Claim 11.1). A representation in
// the sense of §24: it shows what such a coupling would imply, and is offered for exclusion.
//
// A visitor walks a street of frontages 0…n−1 (spacing Δs). Its state is (frontage, need) with need 0 (a
// wedding dress, the errand it came with) or 1 (rings, a need the bridal store raises), and two absorbing
// states, bought (rings bought) and left (the street left). A Markov jump process (physis
// `statistical_mechanics`):
//   - walking: to each neighbouring frontage at D/Δs² (D the wandering diffusivity, m²/s); with need 1 and a
//     search speed v toward the ring store, the hop toward it gains v/(2Δs) and the hop away loses it (a
//     drift v: the store seen or known — the symbolic terrain of §9); a frontage holding a store is entered
//     instead, its rates out 1/(2τ_store) each way (retention, §10.1: a store is a place where the walk stops);
//   - the ends of the street lead off it (to left), and fatigue leaves from anywhere at the rate ℓ;
//   - catalysis: in the bridal store B, while it is open, need 0 → 1 at the rate κ; B is not consumed;
//   - purchase: in the ring store J, while it is open, need 1 → bought at the rate μ.
// The yield Y is the probability of absorption in bought from (s_B, 0): physis's absorption probabilities when
// both are open, and when hours differ, the forward equation stepped by the propagators e^{QΔt} of the four
// configurations (calx's expm), arrivals spread evenly over B's open hours.
//
// Derived, exactly on the representation:
//   - the coupling length ξ_c: the e-folding distance of Y(d) beyond the stores' neighbourhood (by least
//     squares on ln Y over the farther half of the separations, J away from the ends). For a pure walk it is the decay of a killed walk's
//     hitting probability, ξ_c = −Δs / ln λ with λ + 1/λ = 2 + ℓΔs²/D (λ < 1), which tends to √(D/ℓ) as
//     Δs → 0; a search speed lengthens it toward v/ℓ. The half-length (Y at half its value with J adjacent)
//     is reported too: it is set mostly by the step from B straight into a neighbouring J;
//   - the coupling window: Y against an offset of J's hours from B's — J closing before B cuts the pathway at
//     any separation;
//   - creation without capture (Claim 11.1): B's contribution to J's takings is Y v_J per visitor of B, which
//     an allocation of frontage by B's own takings p_B v_B does not register;
//   - the dwell: the visitor's time on the street, phase-type (physis first passage).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/expm.hpp>
#include <physis/statistical_mechanics/first_passage.hpp>
#include <physis/statistical_mechanics/master_equation.hpp>

namespace urban::affinity {

namespace sm = physis::statistical_mechanics;

struct Pathway {
  double length = 600;            // m
  double spacing = 10;            // Δs, m
  double diffusivity = 1.0;       // D, m²/s
  double search_speed = 0;        // v toward J with need 1, m/s
  double leave_rate = 1.0 / 1800; // ℓ, 1/s
  double store_residence = 900;   // τ_store, s
  double catalysis = 1.0 / 1200;  // κ, 1/s
  double purchase = 1.0 / 1800;   // μ, 1/s
  double bridal = 300;            // s_B, m
  double ring = 320;              // s_J, m
  bool ends_exit = true;
  std::size_t frontages() const { return static_cast<std::size_t>(std::lround(length / spacing)) + 1; }
};

// The rates of the process: transient states (frontage i, need m) at 2i + m, then bought, then left.
// b_open, j_open: whether the catalysis and the purchase act.
inline calx::linalg::Matrix pathway_rates(const Pathway& p, bool b_open = true, bool j_open = true) {
  const std::size_t n = p.frontages();
  if (n < 3 || !(p.spacing > 0 && p.diffusivity > 0 && p.store_residence > 0)) throw std::invalid_argument("urban::affinity::pathway_rates: a street of three frontages or more, positive rates");
  const std::size_t ib = static_cast<std::size_t>(std::lround(std::clamp(p.bridal, 0.0, p.length) / p.spacing));
  const std::size_t ij = static_cast<std::size_t>(std::lround(std::clamp(p.ring, 0.0, p.length) / p.spacing));
  const std::size_t bought = 2 * n, left = 2 * n + 1;
  calx::linalg::Matrix k(2 * n + 2);
  const double hop = p.diffusivity / (p.spacing * p.spacing), drift = p.search_speed / (2 * p.spacing);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t m = 0; m < 2; m++) {
      const std::size_t s = 2 * i + m;
      const bool store = i == ib || i == ij;
      for (int dir : {-1, 1}) {
        double rate = store ? 1.0 / (2 * p.store_residence) : hop;
        if (m == 1 && i != ij && !store) {
          const bool toward = (ij > i) == (dir > 0);
          rate = std::max(0.0, rate + (toward ? drift : -drift));
        }
        if (dir < 0 && i == 0) { if (p.ends_exit) k(s, left) += rate; continue; }
        if (dir > 0 && i + 1 == n) { if (p.ends_exit) k(s, left) += rate; continue; }
        const std::size_t t = 2 * (dir < 0 ? i - 1 : i + 1) + m;
        k(s, t) += rate;
      }
      k(s, left) += p.leave_rate;
      if (i == ib && m == 0 && b_open) k(s, 2 * i + 1) += p.catalysis;
      if (i == ij && m == 1 && j_open) k(s, bought) += p.purchase;
    }
  return k;
}

inline std::vector<std::size_t> pathway_transient(const Pathway& p) {
  std::vector<std::size_t> t(2 * p.frontages());
  for (std::size_t i = 0; i < t.size(); i++) t[i] = i;
  return t;
}

inline std::size_t start_state(const Pathway& p) { return 2 * static_cast<std::size_t>(std::lround(std::clamp(p.bridal, 0.0, p.length) / p.spacing)); }

// Y with both stores open: the probability of reaching bought from (s_B, need 0).
inline double pathway_yield(const Pathway& p) {
  const auto k = pathway_rates(p);
  const std::size_t n = p.frontages();
  const auto transient = pathway_transient(p);
  const std::vector<std::size_t> absorbing{2 * n, 2 * n + 1};
  const auto B = sm::absorption_probabilities(k, transient, absorbing);
  return B(start_state(p), 0);
}

// The e-folding length of a killed walk's hitting probability on frontages Δs apart.
inline double walk_coupling_length(double diffusivity, double leave_rate, double spacing) {
  const double c = 2 + leave_rate * spacing * spacing / diffusivity;
  const double lambda = (c - std::sqrt(c * c - 4)) / 2;
  return -spacing / std::log(lambda);
}

struct Coupling {
  std::vector<double> separation, yield;   // Y(d), J placed at s_B + d (or s_B − d if beyond the street)
  double adjacent = 0;                     // Y with J one frontage from B
  double coupling_length = 0;              // ξ_c, m: the e-folding distance of Y beyond three frontages (NaN if too few)
  double half_length = 0;                  // m: where Y falls to half its value with J adjacent (NaN if it never does)
  double dwell_mean = 0, dwell_median = 0, dwell_p90 = 0;   // s, the visitor's time on the street (J adjacent)
};

inline Coupling coupling(Pathway p, const std::vector<double>& separations) {
  Coupling out;
  for (double d : separations) {
    Pathway q = p;
    q.ring = p.bridal + d <= p.length ? p.bridal + d : p.bridal - d;
    if (std::lround(q.ring / q.spacing) == std::lround(q.bridal / q.spacing)) q.ring = q.bridal + q.spacing;
    out.separation.push_back(d);
    out.yield.push_back(pathway_yield(q));
  }
  Pathway a = p;
  a.ring = p.bridal + p.spacing;
  out.adjacent = pathway_yield(a);
  out.half_length = std::nan("");
  for (std::size_t i = 1; i < out.yield.size(); i++)
    if (out.yield[i] <= 0.5 * out.adjacent && out.yield[i - 1] > 0.5 * out.adjacent) {
      const double f = (out.yield[i - 1] - 0.5 * out.adjacent) / (out.yield[i - 1] - out.yield[i]);
      out.half_length = out.separation[i - 1] + f * (out.separation[i] - out.separation[i - 1]);
      break;
    }
  // ξ_c: the slope of ln Y against d over the farther half of the separations of three frontages or more with
  // J five frontages or more from the street's ends (the near field adds terms of order λ^{2d/Δs}; an end near
  // J changes the chance of buying once there).
  std::vector<std::size_t> usable;
  for (std::size_t i = 0; i < out.yield.size(); i++) {
    const double ring = p.bridal + out.separation[i] <= p.length ? p.bridal + out.separation[i] : p.bridal - out.separation[i];
    if (out.separation[i] >= 3 * p.spacing - 1e-9 && out.yield[i] > 0 && ring >= 5 * p.spacing && ring <= p.length - 5 * p.spacing) usable.push_back(i);
  }
  double sx = 0, sy = 0, sxx = 0, sxy = 0, m = 0;
  for (std::size_t q = usable.size() / 2; q < usable.size(); q++) {
    const std::size_t i = usable[q];
    const double x = out.separation[i], y = std::log(out.yield[i]);
    sx += x, sy += y, sxx += x * x, sxy += x * y, m += 1;
  }
  const double slope = m >= 2 ? (m * sxy - sx * sy) / (m * sxx - sx * sx) : 0;
  out.coupling_length = slope < 0 ? -1 / slope : std::nan("");
  // The dwell: phase-type from (s_B, 0), with J adjacent.
  const auto k = pathway_rates(a);
  const auto transient = pathway_transient(a);
  std::vector<double> alpha(transient.size(), 0);
  alpha[start_state(a)] = 1;
  const auto ph = sm::phase_type(k, transient, alpha);
  out.dwell_mean = ph.moment(1);
  out.dwell_median = ph.quantile(0.5);
  out.dwell_p90 = ph.quantile(0.9);
  return out;
}

struct Window {
  std::vector<double> offset, yield;   // hours by which J's day is shifted from B's; Y averaged over arrivals at B
  double b_open = 0, b_close = 0, j_open = 0, j_close = 0;   // hours of the day
};

// Y over a day when B is open [b_open, b_close) and J [j_open + offset, j_close + offset); arrivals at B every
// `arrival_step` hours from b_open to b_close; the forward equation stepped every `step` seconds to the day's
// end, mass still on the street then counted as left.
inline Window coupling_window(const Pathway& p, double b_open, double b_close, double j_open, double j_close, const std::vector<double>& offsets,
                              double step = 300, double arrival_step = 0.5) {
  const std::size_t S = 2 * p.frontages() + 2, bought = S - 2;
  std::array<calx::linalg::Matrix, 4> E;
  for (int c = 0; c < 4; c++) {
    auto Q = sm::generator(pathway_rates(p, c & 1, c & 2));
    for (auto& v : Q.a) v *= step;
    E[static_cast<std::size_t>(c)] = calx::linalg::expm(Q);
  }
  Window out;
  out.b_open = b_open, out.b_close = b_close, out.j_open = j_open, out.j_close = j_close;
  for (double off : offsets) {
    double total = 0;
    std::size_t arrivals = 0;
    for (double t0 = b_open; t0 < b_close - 1e-9; t0 += arrival_step) {
      std::vector<double> prob(S, 0), next(S);
      prob[start_state(p)] = 1;
      for (double t = t0 * 3600; t < 24 * 3600; t += step) {
        const double h = t / 3600;
        const bool b = h >= b_open && h < b_close, j = h >= j_open + off && h < j_close + off;
        const auto& M = E[static_cast<std::size_t>((b ? 1 : 0) | (j ? 2 : 0))];
        std::fill(next.begin(), next.end(), 0);
        for (std::size_t i = 0; i < S; i++) {
          if (prob[i] == 0) continue;
          for (std::size_t q = 0; q < S; q++) next[q] += prob[i] * M(i, q);
        }
        prob.swap(next);
      }
      total += prob[bought];
      arrivals++;
    }
    out.offset.push_back(off);
    out.yield.push_back(arrivals ? total / static_cast<double>(arrivals) : 0);
  }
  return out;
}

// The plume: the law of a visitor's state at times after it leaves B's door with need 0, both stores open. It
// gives P(at frontage i, need m), P(bought) and P(left). The forward equation is stepped by the propagator
// e^{Q h}, h = horizon / steps (calx's expm). A store's frontage holds whoever is inside it. P(bought) at the
// horizon tends to the yield Y as the horizon grows. The masses sum to one at every time, to rounding.
struct Plume {
  std::vector<double> times;                      // s since leaving B
  std::vector<double> position;                   // m, the frontages along the street
  std::vector<std::vector<double>> errand, need;  // per time, per frontage: P(there with need 0), P(there with need 1)
  std::vector<double> bought, left;               // per time
};

inline Plume pathway_plume(const Pathway& p, double horizon = 3600, std::size_t steps = 60) {
  if (!(horizon > 0) || steps == 0) throw std::invalid_argument("urban::affinity::pathway_plume: a positive horizon and at least one step");
  const std::size_t n = p.frontages(), S = 2 * n + 2;
  auto Q = sm::generator(pathway_rates(p));
  const double h = horizon / static_cast<double>(steps);
  for (auto& v : Q.a) v *= h;
  const auto E = calx::linalg::expm(Q);
  Plume out;
  for (std::size_t i = 0; i < n; i++) out.position.push_back(static_cast<double>(i) * p.spacing);
  std::vector<double> prob(S, 0), next(S);
  prob[start_state(p)] = 1;
  for (std::size_t k = 0; k <= steps; k++) {
    out.times.push_back(static_cast<double>(k) * h);
    std::vector<double> e(n), m(n);
    for (std::size_t i = 0; i < n; i++) e[i] = prob[2 * i], m[i] = prob[2 * i + 1];
    out.errand.push_back(std::move(e)), out.need.push_back(std::move(m));
    out.bought.push_back(prob[S - 2]), out.left.push_back(prob[S - 1]);
    if (k == steps) break;
    std::fill(next.begin(), next.end(), 0);
    for (std::size_t i = 0; i < S; i++) {
      if (prob[i] == 0) continue;
      for (std::size_t q = 0; q < S; q++) next[q] += prob[i] * E(i, q);
    }
    prob.swap(next);
  }
  return out;
}

struct Capture {
  double own = 0;          // B's takings per visitor: p_B v_B
  double created = 0;      // B's contribution to J's takings per visitor of B: Y v_J
  double ratio = 0;        // created / own
};

inline Capture capture(double yield, double p_bridal, double v_bridal, double v_ring) {
  Capture c;
  c.own = p_bridal * v_bridal;
  c.created = yield * v_ring;
  c.ratio = c.own > 0 ? c.created / c.own : 0;
  return c;
}

}  // namespace urban::affinity
