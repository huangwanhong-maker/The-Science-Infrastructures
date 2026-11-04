// The people without a model of their flow: the study's own provisional code (README §2).
//
// The passers at a frontage are read from the data alone, with no conservation law, no network of flow and no
// circulation:
//
// 1. Interpolation (kriging). At each hour h of the mean week, y_i(h) = ln(1 + c_i(h)) at the City's counters i (their
//    mean weeks, results/sensors.json of the shared problem) is a Gaussian random field in the plane: a mean μ(h), the
//    hour's mean over the counters, and a deviation of covariance s² [(1 − η) exp(−d/ℓ) + η δ] (an exponential,
//    Matérn ν = ½, covariance with a nugget η for the difference between the two footpaths of one street and the
//    counters' own scatter). The deviation at a frontage's anchor is its posterior mean given the counters, by tyche's
//    Gaussian conditioning (tyche::spatial::GaussianConditioning: simple kriging, Cressie 1993 §3.4.1), the same
//    weights at every hour (one factorisation; each hour's innovations whitened by it). The passers are the posterior
//    median of the count, exp(μ̂) − 1. ℓ and η are chosen by the marginal likelihood, the 168 hours taken as
//    independent replicates of the field and s² profiled out (an approximation: the hours are not independent), on a
//    grid; leave-one-out prediction of every counter's week measures the interpolation's error.
//
// 2. How a plan changes them (an empirical rule, stated as such). Across the counters, ln(1 + c_i(h)) is regressed on
//    ln(1 + A_i(h)), A_i(h) the number of CLUE's 2024 establishments of the trading kinds (retail, food and drink,
//    services, culture) within r of the counter that are open at h by their kind's typical hours, with the hours'
//    168 means absorbed (Frisch–Waugh–Lovell) and errors clustered by counter (tyche's least squares). The slope ε is
//    a cross-sectional elasticity: where more trading establishments are open, more people are counted. A plan then
//    scales a frontage's passers by ((1 + A_f(h; plan)) / (1 + A_f(h; as built)))^ε, clamped to [½, 2], where
//    A_f(h; plan) is CLUE's count about the frontage with the area's ground floors replaced by the plan's: the as built
//    trading ground floors within r of f open at h taken out and the plan's put in.
//
// What the rule cannot see (README §6.6). It conserves nobody: people drawn to an open frontage come from nowhere. It
// knows no network: a link made changes nobody's way, so links enter the value only through their cost and the
// shared problem's access. And it reads a correlation as a response: places with more open shops have more people,
// whichever causes which.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/linalg/rectangular.hpp>
#include <tyche/regression/linear.hpp>
#include <tyche/spatial/conditioning.hpp>

#include "economy.hpp"

namespace without {

using economy::Json;
using economy::Point;
inline constexpr std::size_t H = economy::week;

struct Counter {
  long long id = 0;
  std::string name, way;
  Point at{};
  std::array<double, H> week{};  // people an hour, the mean week
  std::array<double, H> y{};     // ln(1 + week)
};

inline std::vector<Counter> counters_of(const Json& sensors) {
  std::vector<Counter> out;
  for (const auto& s : sensors["sensors"].as_array()) {
    Counter c;
    c.id = static_cast<long long>(s["id"].as_number()), c.name = s["name"].text(), c.way = s["way"]["kind"].text();
    c.at = {s["x"].as_number(), s["y"].as_number()};
    const auto& w = s["mean_week"].as_array();
    for (std::size_t h = 0; h < H; h++) c.week[h] = std::max(0.0, w[h].as_number()), c.y[h] = std::log1p(c.week[h]);
    out.push_back(std::move(c));
  }
  return out;
}

// ---- 1. kriging ---------------------------------------------------------------------------------------------------

class Kriging {
 public:
  double length = 100, nugget = 0.1, scale = 1;  // ℓ (m), η, s²
  std::array<double, H> mean{};                  // μ(h)
  double log_likelihood = -std::numeric_limits<double>::infinity();

  Kriging() = default;
  Kriging(const std::vector<Counter>& c, std::vector<std::size_t> use, double length_m, double nugget_share) : length(length_m), nugget(nugget_share), use_(std::move(use)) {
    for (std::size_t h = 0; h < H; h++) {
      double m = 0;
      for (auto i : use_) m += c[i].y[h];
      mean[h] = m / static_cast<double>(use_.size());
    }
    for (std::size_t a = 0; a < use_.size(); a++) {
      std::vector<double> k;
      for (std::size_t b = 0; b < a; b++) k.push_back(signal(dist(c[use_[a]].at, c[use_[b]].at)));
      gc_.observe(k, 1 - nugget, nugget, 0.0);
      at_.push_back(c[use_[a]].at);
    }
    const double log_det_half = -gc_.log_likelihood() - 0.5 * static_cast<double>(use_.size()) * std::log(2 * std::numbers::pi);
    double q = 0;
    z_.assign(H, {});
    for (std::size_t h = 0; h < H; h++) {
      std::vector<double> d;
      for (auto i : use_) d.push_back(c[i].y[h] - mean[h]);
      z_[h] = gc_.whiten(d);
      for (double v : z_[h]) q += v * v;
    }
    const double n = static_cast<double>(use_.size()) * H;
    scale = q / n;
    log_likelihood = -0.5 * n * (std::log(2 * std::numbers::pi * scale) + 1) - static_cast<double>(H) * log_det_half;
  }

  // The posterior mean of ln(1 + count) at a point, hour by hour, and the posterior variance (the smooth field's).
  std::array<double, H> predict(const Point& q, double* variance = nullptr) const {
    std::vector<double> k;
    for (const auto& p : at_) k.push_back(signal(dist(q, p)));
    const auto w = gc_.whiten(k);
    std::array<double, H> out{};
    for (std::size_t h = 0; h < H; h++) {
      double m = mean[h];
      for (std::size_t i = 0; i < w.size(); i++) m += w[i] * z_[h][i];
      out[h] = m;
    }
    if (variance) {
      double v = 1 - nugget;
      for (double x : w) v -= x * x;
      *variance = scale * std::max(0.0, v);
    }
    return out;
  }

  static double dist(const Point& a, const Point& b) { return std::hypot(a[0] - b[0], a[1] - b[1]); }

 private:
  double signal(double d) const { return (1 - nugget) * std::exp(-d / length); }
  std::vector<std::size_t> use_;
  std::vector<Point> at_;
  tyche::spatial::GaussianConditioning gc_;
  std::vector<std::vector<double>> z_;
};

struct KrigingFit {
  Kriging best;
  std::vector<std::array<double, 3>> grid;  // ℓ, η, log-likelihood
  // Leave-one-out: each counter's week predicted from the others.
  std::vector<std::array<double, H>> loo;    // ln(1 + count)
  double loo_rmse = 0, loo_r2 = 0, loo_weekly_median_ratio = 0;  // in ln(1 + count); R² over all hours and counters; median of predicted/observed weekly totals
};

inline KrigingFit fit_kriging(const std::vector<Counter>& c) {
  KrigingFit F;
  std::vector<std::size_t> all(c.size());
  std::iota(all.begin(), all.end(), std::size_t{0});
  for (double l : {20.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0, 400.0, 600.0, 900.0})
    for (double eta : {0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6}) {
      Kriging k;
      try {
        k = Kriging(c, all, l, eta);
      } catch (const std::exception&) {
        continue;  // a correlation matrix that is not positive definite (η = 0 with coincident counters)
      }
      F.grid.push_back({l, eta, k.log_likelihood});
      if (k.log_likelihood > F.best.log_likelihood) F.best = k;
    }
  // Leave one out at the chosen ℓ and η.
  double se = 0, ss = 0, n = 0, mean = 0;
  for (const auto& x : c)
    for (double v : x.y) mean += v, n += 1;
  mean /= n;
  std::vector<double> ratios;
  for (std::size_t i = 0; i < c.size(); i++) {
    std::vector<std::size_t> others;
    for (std::size_t j = 0; j < c.size(); j++)
      if (j != i) others.push_back(j);
    const Kriging k(c, others, F.best.length, F.best.nugget);
    const auto p = k.predict(c[i].at);
    F.loo.push_back(p);
    double po = 0, ob = 0;
    for (std::size_t h = 0; h < H; h++) {
      se += (p[h] - c[i].y[h]) * (p[h] - c[i].y[h]), ss += (c[i].y[h] - mean) * (c[i].y[h] - mean);
      po += std::expm1(p[h]), ob += c[i].week[h];
    }
    ratios.push_back(po / ob);
  }
  F.loo_rmse = std::sqrt(se / n), F.loo_r2 = 1 - se / ss;
  std::sort(ratios.begin(), ratios.end());
  F.loo_weekly_median_ratio = ratios[ratios.size() / 2];
  return F;
}

// ---- 2. the empirical rule -------------------------------------------------------------------------------------------

// The trading establishments open at hour h within r of a point, by their kind's typical hours.
struct Establishments {
  std::vector<Point> at;
  std::vector<std::int64_t> kind;
};

inline Establishments establishments_of(const Json& inputs) {
  Establishments E;
  for (const auto& e : inputs["establishments"]["rows"].as_array()) {
    const auto k = static_cast<std::int64_t>(e[2].as_number());
    if (!economy::trades(k)) continue;
    E.at.push_back({e[0].as_number(), e[1].as_number()}), E.kind.push_back(k);
  }
  return E;
}

// The typical hours (start, length on a weekday and at the weekend) by kind: the shared problem's hours as built.
inline std::array<std::array<std::int64_t, 4>, economy::kinds> typical_hours(const economy::Problem& P) {
  std::array<std::array<std::int64_t, 4>, economy::kinds> t{};
  std::array<bool, economy::kinds> seen{};
  for (const auto& f : P.frontages)
    if (!seen[static_cast<std::size_t>(f.kind0)]) t[static_cast<std::size_t>(f.kind0)] = f.hours0, seen[static_cast<std::size_t>(f.kind0)] = true;
  return t;
}

inline bool open_by(const std::array<std::int64_t, 4>& hrs, std::size_t h) {
  const std::size_t o = economy::is_weekend(h) ? 2 : 0;
  return economy::open_in(hrs[o], hrs[o + 1], economy::hour_of_day(h));
}

inline std::array<double, H> active_near(const Establishments& E, const std::array<std::array<std::int64_t, 4>, economy::kinds>& hours, const Point& q, double r) {
  std::array<double, H> a{};
  for (std::size_t e = 0; e < E.at.size(); e++) {
    if (Kriging::dist(E.at[e], q) > r) continue;
    for (std::size_t h = 0; h < H; h++) a[h] += open_by(hours[static_cast<std::size_t>(E.kind[e])], h);
  }
  return a;
}

struct Elasticity {
  double radius = 50, slope = 0, se = 0, t = 0, r2_within = 0;
  std::size_t rows = 0, clusters = 0;
};

inline Elasticity fit_elasticity(const std::vector<Counter>& c, const Establishments& E, const std::array<std::array<std::int64_t, 4>, economy::kinds>& hours, double r) {
  std::vector<double> y, x;
  std::vector<std::int64_t> cl;
  std::vector<std::size_t> hour;
  for (std::size_t i = 0; i < c.size(); i++) {
    const auto a = active_near(E, hours, c[i].at, r);
    for (std::size_t h = 0; h < H; h++) y.push_back(c[i].y[h]), x.push_back(std::log1p(a[h])), cl.push_back(static_cast<std::int64_t>(i)), hour.push_back(h);
  }
  // The hours' means absorbed.
  const auto demean = [&](std::vector<double>& v) {
    std::array<double, H> s{}, n{};
    for (std::size_t j = 0; j < v.size(); j++) s[hour[j]] += v[j], n[hour[j]] += 1;
    for (std::size_t j = 0; j < v.size(); j++) v[j] -= s[hour[j]] / n[hour[j]];
  };
  demean(y), demean(x);
  tyche::regression::OlsOptions o;
  o.covariance = tyche::regression::Covariance::cluster;
  o.clusters = cl;
  const auto fit = tyche::regression::ols(calx::linalg::Rectangular::from_columns({x}), y, o);
  Elasticity e;
  e.radius = r, e.slope = fit.beta[0], e.se = fit.se[0], e.t = fit.t[0], e.r2_within = fit.r_squared, e.rows = y.size(), e.clusters = fit.clusters;
  return e;
}

// ---- the flows of a plan ---------------------------------------------------------------------------------------------

class Flows {
 public:
  std::vector<double> base;                 // passers as built, f * 168 + h
  std::vector<std::array<double, H>> clue;  // CLUE's trading establishments open within r of each frontage, as built
  std::vector<std::vector<std::size_t>> near;  // frontages within r of each frontage's anchor (itself included)
  Elasticity rule;
  double floor = 0.5, ceiling = 2.0;

  Flows(const economy::Problem& P, const Kriging& K, const Establishments& E, const Elasticity& e) : rule(e), P_(&P), built_(P.as_built()) {
    const std::size_t n = P.size();
    const auto hours = typical_hours(P);
    base.assign(n * H, 0.0);
    for (std::size_t f = 0; f < n; f++) {
      const auto m = K.predict(P.frontages[f].anchor);
      for (std::size_t h = 0; h < H; h++) base[f * H + h] = std::expm1(m[h]);
      clue.push_back(active_near(E, hours, P.frontages[f].anchor, e.radius));
      std::vector<std::size_t> nb;
      for (std::size_t g = 0; g < n; g++)
        if (Kriging::dist(P.frontages[f].anchor, P.frontages[g].anchor) <= e.radius) nb.push_back(g);
      near.push_back(std::move(nb));
    }
  }

  // The passers of a plan: the rule's multiplier on the passers as built.
  std::vector<double> operator()(const economy::Plan& plan) const {
    const std::size_t n = P_->size();
    std::vector<double> out(base);
    for (std::size_t f = 0; f < n; f++)
      for (std::size_t h = 0; h < H; h++) {
        double delta = 0;
        for (auto g : near[f]) delta += active(plan, g, h) - active(built_, g, h);
        if (delta == 0) continue;
        const double a0 = clue[f][h], a1 = std::max(0.0, a0 + delta);
        out[f * H + h] *= std::clamp(std::pow((1 + a1) / (1 + a0), rule.slope), floor, ceiling);
      }
    return out;
  }

 private:
  static double active(const economy::Plan& p, std::size_t g, std::size_t h) { return economy::trades(p.kind[g]) && p.open(g, h) ? 1.0 : 0.0; }
  const economy::Problem* P_;
  economy::Plan built_;
};

}  // namespace without
