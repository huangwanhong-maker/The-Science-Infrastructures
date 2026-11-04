// The structured search of study (2): the value read hour by hour, and a frontage's exact best response (README §3).
//
// THE STRUCTURE. The shared problem's value is a sum over the week's hours of the margin captured in that hour, which
// depends only on the frontages open in it, their kinds and the passers then, less costs that are each a frontage's
// own (staffing by its hours, outgoings and refit by its kind) or a link's:
//
//   V(plan) = Σ_h M_h(open set at h, kinds, P(h)) − Σ_f cost_f(k_f, hours_f) − Σ_links cost.
//
// The flow model makes the passers affine in the openings, P(h) = B(h) + c(h) Σ_g w_{k_g} [g open at h] R_{·g}, so
// a frontage's opening at h changes M_h and nothing else. Holding every other frontage fixed, the value of frontage f
// holding kind k with weekday hours [s, e) and weekend hours [s', e') is therefore
//
//   Σ_{t ∈ [s, e)} Σ_{weekdays d} gain_k(24 d + t) + Σ_{t ∈ [s', e')} Σ_{weekend d} gain_k(24 d + t) − cost_f(k, hours),
//
// gain_k(h) = M_h(f open as k) − M_h(f closed): exact, not a linearisation. Its maximum over the intervals is found by
// enumerating the 300 intervals of a day, once a day-type, so that f's best response over its 6 kinds and every pair of
// hours costs 6 × 168 margins of an hour, the work of six evaluations of the value. Sweeping the frontages so, each in
// turn taking its best response, is block-coordinate ascent: the value never falls, and it stops where no frontage alone
// can do better.
//
// THE BUDGET. Work is counted in margins of an hour: one evaluation of a plan's value is 168 of them, so a best response
// costs 6 evaluations and a sweep of the 61 frontages 366.

#pragma once

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include <calx/random/stream.hpp>

#include "economy.hpp"
#include "flow.hpp"

namespace structured {

using economy::Plan;
inline constexpr std::size_t H = economy::week;

class Engine {
 public:
  const economy::Problem& P;
  const flow::Model& M;
  mutable std::atomic<long long> margins{0};  // margins of an hour computed

  Engine(const economy::Problem& p, const flow::Model& m) : P(p), M(m) {}

  double evaluations() const { return static_cast<double>(margins.load()) / H; }

  // The margin captured in hour h (A$), as economy::Problem::evaluate computes it, given the frontages open then, the
  // kinds and the passers then. `p` is scratch of one entry a frontage.
  double hour(const Plan& plan, const economy::Distances& D, std::size_t h, const std::vector<char>& open, const double* passers, std::vector<double>& p) const {
    margins++;
    const std::size_t n = P.size();
    const auto& Q = P.P;
    for (std::size_t f = 0; f < n; f++) {
      p[f] = 0;
      if (!open[f] || !economy::trades(plan.kind[f])) continue;
      const auto k = static_cast<std::size_t>(plan.kind[f]);
      double kappa = 0;
      for (const auto& [g, w] : D.near[f])
        if (open[g]) kappa += w * Q.affinity[static_cast<std::size_t>(plan.kind[g])][k];
      p[f] = Q.propensity[k] * Q.period[k][economy::band(h)] * std::clamp(1 + Q.catalysis * kappa, Q.catalysis_floor, Q.catalysis_ceiling);
    }
    double margin = 0;
    for (std::size_t f = 0; f < n; f++) {
      if (p[f] <= 0) continue;
      double none = 1 - p[f], sum = p[f];
      for (auto g : P.on_segment[f])
        if (p[g] > 0 && plan.kind[g] == plan.kind[f]) none *= 1 - p[g], sum += p[g];
      const auto k = static_cast<std::size_t>(plan.kind[f]);
      margin += passers[f] * p[f] * ((1 - none) / sum) * Q.spend[k] * Q.margin[k];
    }
    return margin;
  }

  // Frontage f's own costs (staffing, outgoings, refit) with kind k and the plan's hours for it.
  double frontage_cost(const Plan& plan, std::size_t f) const {
    const auto& Q = P.P;
    const auto k = static_cast<std::size_t>(plan.kind[f]);
    if (plan.kind[f] == economy::vacant) return 0;
    double c = Q.staff[k] * Q.wage * (plan.weekly_hours(f) + plan.days_open(f)) + Q.outgoings[k] * P.frontages[f].trading_area;
    if (plan.kind[f] != P.frontages[f].kind0) c += Q.refit_weekly(P.frontages[f].trading_area);
    return c;
  }
  double link_cost(const Plan& plan) const {
    double c = 0;
    for (std::size_t k = 0; k < P.links.size(); k++)
      if (plan.links[k]) c += P.P.link_weekly(P.links[k].length);
    return c;
  }

  // The hour's passers of every frontage, for the plan's openings.
  void hour_passers(const Plan& plan, const flow::Response& R, std::size_t h, const std::vector<char>& open, std::vector<double>& out) const {
    const std::size_t n = P.size();
    const double c = M.x[flow::visit][h] / M.W_ref;
    for (std::size_t f = 0; f < n; f++) out[f] = R.B[f * H + h];
    if (c == 0) return;
    for (std::size_t g = 0; g < n; g++) {
      if (!open[g]) continue;
      const double w = M.o.attraction[static_cast<std::size_t>(plan.kind[g])];
      if (w == 0) continue;
      for (std::size_t f = 0; f < n; f++) out[f] += c * w * R.R[f * n + g];
    }
  }

  // The value of a plan, hour by hour (equal to economy::Problem::evaluate's value; the checks hold them equal).
  double value(const Plan& plan, const flow::Response& R) const {
    const std::size_t n = P.size();
    std::vector<char> open(n);
    std::vector<double> pass(n), p(n);
    double v = 0;
    for (std::size_t h = 0; h < H; h++) {
      for (std::size_t f = 0; f < n; f++) open[f] = plan.open(f, h);
      hour_passers(plan, R, h, open, pass);
      v += hour(plan, *R.distances, h, open, pass.data(), p);
    }
    for (std::size_t f = 0; f < n; f++) v -= frontage_cost(plan, f);
    return v - link_cost(plan);
  }

  struct Best {
    double value = -std::numeric_limits<double>::infinity();  // relative to f closed at every hour
    std::int64_t kind = economy::vacant;
    std::array<std::int64_t, 4> hours{0, 0, 0, 0};
  };

  // The best interval [s, e) of a day for gains g(t) less `per_hour` an hour open and `per_day` if open at all.
  static std::pair<double, std::array<std::int64_t, 2>> best_interval(const std::array<double, 24>& g, double per_hour, double per_day) {
    double best = 0;
    std::array<std::int64_t, 2> se{0, 0};
    for (std::int64_t s = 0; s < 24; s++) {
      double run = -per_day;
      for (std::int64_t e = s + 1; e <= 24; e++) {
        run += g[static_cast<std::size_t>(e - 1)] - per_hour;
        if (run > best + 1e-9) best = run, se = {s, e - s};
      }
    }
    return {best, se};
  }

  struct Gains {
    std::array<std::array<double, H>, economy::kinds> gain{};  // M_h(f open as k) − M_h(f closed)
  };

  // gain_k(h) for frontage f, every kind, the rest of the plan fixed.
  Gains gains(const Plan& plan, const flow::Response& R, std::size_t f) const {
    const std::size_t n = P.size();
    Gains G;
    Plan q = plan;
    std::vector<char> open(n);
    std::vector<double> base(n), pass(n), p(n);
    for (std::size_t h = 0; h < H; h++) {
      for (std::size_t g = 0; g < n; g++) open[g] = plan.open(g, h);
      open[f] = 0;
      hour_passers(plan, R, h, open, base);
      const double closed = hour(plan, *R.distances, h, open, base.data(), p);
      const double c = M.x[flow::visit][h] / M.W_ref;
      open[f] = 1;
      for (std::size_t k = 1; k < economy::kinds; k++) {
        q.kind[f] = static_cast<std::int64_t>(k);
        const double w = M.o.attraction[k];
        for (std::size_t g = 0; g < n; g++) pass[g] = base[g] + c * w * R.R[g * n + f];
        G.gain[k][h] = hour(q, *R.distances, h, open, pass.data(), p) - closed;
      }
    }
    return G;
  }

  // The best (kind, hours) of frontage f alone, from its gains, and the value of each kind with its best hours.
  std::pair<Best, std::array<Best, economy::kinds>> choose(const Plan& plan, std::size_t f, const Gains& G) const {
    const auto& Q = P.P;
    std::array<Best, economy::kinds> per;
    per[economy::vacant] = {0.0, economy::vacant, {plan.hours[4 * f], 0, plan.hours[4 * f + 2], 0}};
    Best best = per[economy::vacant];
    for (std::size_t k = 1; k < economy::kinds; k++) {
      const double hourly = Q.staff[k] * Q.wage;
      std::array<double, 24> wd{}, we{};
      for (std::size_t t = 0; t < 24; t++) {
        for (std::size_t d = 0; d < 5; d++) wd[t] += G.gain[k][24 * d + t];
        for (std::size_t d = 5; d < 7; d++) we[t] += G.gain[k][24 * d + t];
      }
      const auto [vw, sw] = best_interval(wd, 5 * hourly, 5 * hourly);
      const auto [ve, se] = best_interval(we, 2 * hourly, 2 * hourly);
      double v = vw + ve - Q.outgoings[k] * P.frontages[f].trading_area;
      if (static_cast<std::int64_t>(k) != P.frontages[f].kind0) v -= Q.refit_weekly(P.frontages[f].trading_area);
      per[k] = {v, static_cast<std::int64_t>(k), {sw[0], sw[1], se[0], se[1]}};
      if (v > best.value + 1e-9) best = per[k];
    }
    return {best, per};
  }
  std::pair<Best, std::array<Best, economy::kinds>> best_response(const Plan& plan, const flow::Response& R, std::size_t f) const { return choose(plan, f, gains(plan, R, f)); }

  // The present (kind, hours) of f valued on the same gains, relative to f closed.
  double present(const Plan& plan, std::size_t f, const Gains& G) const {
    if (plan.kind[f] == economy::vacant) return 0;
    double v = 0;
    for (std::size_t h = 0; h < H; h++)
      if (plan.open(f, h)) v += G.gain[static_cast<std::size_t>(plan.kind[f])][h];
    return v - frontage_cost(plan, f);
  }

  // One sweep of best responses in the given order; returns the number of frontages changed.
  std::size_t sweep(Plan& plan, const flow::Response& R, const std::vector<std::size_t>& order) const {
    std::size_t changed = 0;
    for (auto f : order) {
      const auto G = gains(plan, R, f);
      const auto [best, per] = choose(plan, f, G);
      (void)per;
      if (best.value > present(plan, f, G) + 1e-6) {
        plan.kind[f] = best.kind;
        for (std::size_t j = 0; j < 4; j++) plan.hours[4 * f + j] = best.hours[j];
        changed++;
      }
    }
    return changed;
  }
};

}  // namespace structured
