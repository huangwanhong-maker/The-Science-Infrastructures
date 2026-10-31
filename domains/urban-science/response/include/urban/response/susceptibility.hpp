// Forcing and response (theory/street-oceanography.md §4; Huang, Street Oceanography, §13): how much of an
// irregular or dated forcing a street converts, and the comparison across a repeating alteration.
//
// The response function (Definition 4.1). On the record's hours with a value and every forcing's lags,
//   ln(1 + y_t) = b(hour of the civil week) + a₁ cos 2πθ_Y + b₁ sin 2πθ_Y + a₂ cos 4πθ_Y + b₂ sin 4πθ_Y + c t
//                 + Σ_f Σ_ℓ β_{fℓ} F_f(t − ℓ) + ε_t.
// The 168 levels b are absorbed by demeaning within the hours of the week (Frisch–Waugh–Lovell: the slopes
// and their residuals are those of the full regression; Davidson and MacKinnon 2004, §2.4), and the rest is
// fitted by tyche's least squares with standard errors clustered by civil day (the residuals of one day are
// correlated; Cameron and Miller 2015, J. Hum. Resour. 50, 317). A forcing's susceptibility is Σ_ℓ β_{fℓ},
// the proportional response to a unit of it sustained over its lags, with variance 1ᵀ C 1 over its block.
// Caveat: the degrees of freedom absorbed by the 168 levels are not subtracted from the clustered
// inference, which counts days (G − 1); with years of hours the difference is negligible.
//
// The comparison across a repeating alteration (Definition 4.3). Each day of configuration A is paired with
// the nearest day of configuration B of the same weekday, within a fortnight; the effect is the mean of the
// paired differences of ln(1 + daily total), its significance by flipping signs within pairs (randomisation
// inference; Fisher 1935, §21), Monte Carlo with calx's seeded stream.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/linalg/rectangular.hpp>
#include <calx/random/stream.hpp>
#include <tyche/regression/linear.hpp>
#include <urban/rhythm/calendar.hpp>
#include <urban/rhythm/record.hpp>

namespace urban::response {

using urban::rhythm::Record;

struct Forcing {
  std::string name, unit;
  std::vector<double> value;     // one per hour of the record (NaN where unknown)
  std::vector<int> lags;         // hours, each ≥ 0
};

struct Susceptibility {
  std::string name, unit;
  std::vector<int> lags;
  std::vector<double> beta, se;  // per lag
  double total = 0, total_se = 0, z = 0;
  double percent = 0;            // 100 (e^{total} − 1): the proportional change in counts per unit sustained
};

struct ResponseFit {
  std::vector<Susceptibility> forcings;
  std::size_t observations = 0, days = 0;
  double r2_within = 0;          // R² of the demeaned regression: the share of the within-hour variance explained
  double sigma = 0;
  std::vector<double> annual;    // a₁, b₁, a₂, b₂
  double trend = 0;              // per year
};

inline ResponseFit respond(const Record& r, const std::vector<Forcing>& forcings, bool annual = true, bool trend = true) {
  for (const auto& f : forcings)
    if (f.value.size() != r.size() || f.lags.empty()) throw std::invalid_argument("urban::response::respond: each forcing has a value for every hour and a lag");
  std::size_t max_lag = 0;
  for (const auto& f : forcings)
    for (int l : f.lags) {
      if (l < 0) throw std::invalid_argument("urban::response::respond: lags are not negative");
      max_lag = std::max(max_lag, static_cast<std::size_t>(l));
    }
  // The rows: hours with a value and every lagged forcing known.
  std::vector<std::size_t> rows;
  for (std::size_t i = max_lag; i < r.size(); i++) {
    if (!std::isfinite(r.y[i])) continue;
    bool ok = true;
    for (const auto& f : forcings)
      for (int l : f.lags) ok = ok && std::isfinite(f.value[i - static_cast<std::size_t>(l)]);
    if (ok) rows.push_back(i);
  }
  if (rows.size() < 500) throw std::invalid_argument("urban::response::respond: fewer than 500 hours with a value and every forcing");
  // Columns: the forcings' lags, then the annual harmonics and the trend.
  std::vector<std::vector<double>> cols;
  for (const auto& f : forcings)
    for (int l : f.lags) {
      std::vector<double> c;
      for (auto i : rows) c.push_back(f.value[i - static_cast<std::size_t>(l)]);
      cols.push_back(std::move(c));
    }
  const std::size_t n_forcing = cols.size();
  if (annual)
    for (int h = 1; h <= 2; h++)
      for (int s = 0; s < 2; s++) {
        std::vector<double> c;
        for (auto i : rows) {
          const double a = 2 * std::numbers::pi * h * urban::rhythm::year_angle(r.time(i));
          c.push_back(s == 0 ? std::cos(a) : std::sin(a));
        }
        cols.push_back(std::move(c));
      }
  if (trend) {
    std::vector<double> c;
    for (auto i : rows) c.push_back((r.time(i) - r.t0) / (365.25 * 86400.0));
    cols.push_back(std::move(c));
  }
  std::vector<double> y;
  std::vector<int> hour;
  std::vector<std::int64_t> day;
  for (auto i : rows) {
    y.push_back(std::log1p(std::max(0.0, r.y[i])));
    const double t = r.time(i);
    hour.push_back(static_cast<int>(std::floor(urban::rhythm::week_angle(r.civil(t)) * 168.0 + 1e-9)) % 168);
    day.push_back(r.civil_day(t));
  }
  // Demean within the hours of the week.
  auto demean = [&](std::vector<double>& v) {
    std::vector<double> sum(168, 0), count(168, 0);
    for (std::size_t j = 0; j < v.size(); j++) sum[hour[j]] += v[j], count[hour[j]] += 1;
    for (std::size_t j = 0; j < v.size(); j++) v[j] -= sum[hour[j]] / count[hour[j]];
  };
  demean(y);
  for (auto& c : cols) demean(c);
  const auto X = calx::linalg::Rectangular::from_columns(cols);
  tyche::regression::OlsOptions o;
  o.covariance = tyche::regression::Covariance::cluster;
  o.clusters = day;
  o.drop_collinear = true;
  const auto fit = tyche::regression::ols(X, y, o);
  ResponseFit out;
  out.observations = rows.size();
  out.days = fit.clusters;
  out.r2_within = fit.tss > 0 ? 1 - fit.rss / fit.tss : 0;
  out.sigma = fit.sigma;
  std::size_t c = 0;
  for (const auto& f : forcings) {
    Susceptibility s;
    s.name = f.name, s.unit = f.unit, s.lags = f.lags;
    const std::size_t first = c;
    for (std::size_t l = 0; l < f.lags.size(); l++, c++) {
      s.beta.push_back(fit.beta[c]);
      s.se.push_back(fit.se[c]);
      if (std::isfinite(fit.beta[c])) s.total += fit.beta[c];
    }
    double v = 0;
    for (std::size_t a = first; a < c; a++)
      for (std::size_t b = first; b < c; b++)
        if (std::isfinite(fit.covariance(a, b))) v += fit.covariance(a, b);
    s.total_se = std::sqrt(std::max(0.0, v));
    s.z = s.total_se > 0 ? s.total / s.total_se : 0;
    s.percent = 100 * std::expm1(s.total);
    out.forcings.push_back(std::move(s));
  }
  if (annual)
    for (int k = 0; k < 4; k++) out.annual.push_back(fit.beta[n_forcing + static_cast<std::size_t>(k)]);
  if (trend) out.trend = fit.beta[cols.size() - 1];
  return out;
}

// ---- the comparison across a repeating alteration ----------------------------------------------------------

struct PairedEffect {
  std::size_t pairs = 0;
  double effect = 0;             // mean of ln(1 + y_A) − ln(1 + y_B), daily totals
  double se = 0;
  double percent = 0;            // 100 (e^{effect} − 1)
  double p = 1;                  // two-sided, sign flips
  std::size_t permutations = 0;
  std::vector<double> differences;
  std::vector<long long> a_days, b_days;
};

// Daily totals by civil day (days with every hour known only).
inline std::map<long long, double> daily_totals(const Record& r) {
  std::map<long long, double> sum;
  std::map<long long, int> hours;
  for (std::size_t i = 0; i < r.size(); i++) {
    if (!std::isfinite(r.y[i])) continue;
    const long long d = r.civil_day(r.time(i));
    sum[d] += r.y[i];
    hours[d] += 1;
  }
  std::map<long long, double> out;
  for (const auto& [d, s] : sum)
    if (hours[d] >= 23) out[d] = s;   // 23: the day the clocks go forward
  return out;
}

// A-days (civil days since 1970) paired with the nearest B-day of the same weekday within 14 days: B-days are
// the days not in A (and not in `exclude`).
inline PairedEffect paired_effect(const Record& r, const std::vector<long long>& a_days, const std::vector<long long>& exclude = {},
                                  std::size_t permutations = 9999, std::uint64_t seed = 1) {
  const auto totals = daily_totals(r);
  auto banned = [&](long long d) {
    return std::find(a_days.begin(), a_days.end(), d) != a_days.end() || std::find(exclude.begin(), exclude.end(), d) != exclude.end();
  };
  PairedEffect out;
  for (long long a : a_days) {
    const auto ia = totals.find(a);
    if (ia == totals.end()) continue;
    for (int w : {7, -7, 14, -14}) {
      const long long b = a + w;
      const auto ib = totals.find(b);
      if (ib == totals.end() || banned(b)) continue;
      out.differences.push_back(std::log1p(ia->second) - std::log1p(ib->second));
      out.a_days.push_back(a), out.b_days.push_back(b);
      break;
    }
  }
  out.pairs = out.differences.size();
  if (out.pairs < 2) return out;
  double mean = 0;
  for (double d : out.differences) mean += d;
  mean /= static_cast<double>(out.pairs);
  double ss = 0;
  for (double d : out.differences) ss += (d - mean) * (d - mean);
  out.effect = mean;
  out.se = std::sqrt(ss / static_cast<double>(out.pairs - 1) / static_cast<double>(out.pairs));
  out.percent = 100 * std::expm1(mean);
  calx::random::Stream rng(seed, 0xfa, 0x17, 0);
  std::size_t extreme = 0;
  for (std::size_t k = 0; k < permutations; k++) {
    double s = 0;
    for (double d : out.differences) s += rng.bernoulli(0.5) ? d : -d;
    if (std::fabs(s / static_cast<double>(out.pairs)) >= std::fabs(mean) - 1e-15) extreme++;
  }
  out.permutations = permutations;
  out.p = (1.0 + static_cast<double>(extreme)) / (1.0 + static_cast<double>(permutations));
  return out;
}

}  // namespace urban::response
