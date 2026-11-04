// The regimes of economies, perceived from their growth: a hidden Markov model of the business cycle fitted to many
// economies at once, and each economy's belief about its regime, year by year.
//
// An economy's annual growth of real GDP is read as the emission of a hidden regime, a Markov chain on K states
// (recession, normal growth, boom, for K = 3), each regime's growth normal (Hamilton 1989, "A new approach to the
// economic analysis of nonstationary time series and the business cycle", Econometrica 57, 357–384, for the
// quarterly growth of the United States; here annual, and pooled). The model is fitted jointly to the economies' series
// by tyche's EM (one transition matrix and one law of growth a regime, shared), the regimes numbered by their mean
// growth, so that state 0 is the worst. Each economy's perception of its regime is the filtered law (what the years so
// far tell), the smoothed law (what the whole record tells) and the most probable sequence of regimes.
//
// Standardised. Pooled on raw growth, the regimes separate kinds of economy (steady, volatile, fast-growing), not the
// phases of any one economy's cycle: the United States then sits in one regime through 1982, 2009 and 2020. Each
// series is therefore read against its own normal, z = (g − median) / (1.4826 MAD), held to ±6, before the pooled fit;
// the regimes are then "well below this economy's normal", "normal" and "well above", in its own deviations.
// With a deviation shared by the regimes (Hamilton's model switches only the mean), the regimes differ by their mean
// growth alone: without it, the worst regime takes the turbulent years, deep falls and their rebounds alike.
//
// Caveats. Sharing one model across economies treats a contraction in a small, volatile economy and in a large, steady
// one alike, in units of each one's spread. The World Bank's growth series begin at different years and have gaps;
// missing years only predict.

#pragma once

#include <tyche/processes/hidden_markov.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace economics {

struct GrowthSeries {
  std::string economy;      // a code (ISO 3166-1 alpha-3)
  int first_year = 0;
  std::vector<double> growth;  // per cent a year, NaN where missing
};

struct RegimeReading {
  std::string economy;
  int first_year = 0;
  std::vector<double> growth;          // per cent
  double centre = 0, scale = 1;        // the economy's median growth and robust deviation (standardised fits)
  std::vector<double> z;               // what the model read: (growth − centre) / scale, or growth
  std::vector<std::vector<double>> filtered, smoothed;
  std::vector<std::size_t> path;
  double log_likelihood = 0;
};

struct Regimes {
  tyche::processes::GaussianHmm model;
  double log_likelihood = 0;
  int iterations = 0;
  bool converged = false;
  std::vector<RegimeReading> readings;
  std::size_t observations = 0;
  bool standardised = false;           // the regimes' units: each economy's deviations, or per cent

  // Mean years in regime k: 1 / (1 − A_kk).
  double duration(std::size_t k) const { return 1 / (1 - model.a(k, k)); }
  // The stationary law of the regimes (the long-run share of years in each).
  std::vector<double> stationary() const {
    const std::size_t K = model.states;
    std::vector<double> p(K, 1.0 / static_cast<double>(K)), q(K);
    for (int it = 0; it < 10000; it++) {
      double change = 0;
      for (std::size_t j = 0; j < K; j++) {
        q[j] = 0;
        for (std::size_t i = 0; i < K; i++) q[j] += p[i] * model.a(i, j);
      }
      for (std::size_t j = 0; j < K; j++) change = std::max(change, std::fabs(q[j] - p[j]));
      p.swap(q);
      if (change < 1e-15) break;
    }
    return p;
  }
};

namespace detail {
inline double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const std::size_t n = v.size();
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
}  // namespace detail

// The pooled fit and every economy's reading. Series with fewer than `least` observed years are left out.
inline Regimes perceive_regimes(const std::vector<GrowthSeries>& series, std::size_t K = 3, std::size_t least = 10, bool standardise = false,
                                 bool shared_sd = false) {
  std::vector<std::vector<double>> seqs;
  std::vector<RegimeReading> readings;
  Regimes R;
  R.standardised = standardise;
  for (const auto& s : series) {
    std::vector<double> seen;
    for (const double g : s.growth)
      if (!std::isnan(g)) seen.push_back(g);
    if (seen.size() < least) continue;
    RegimeReading r;
    r.economy = s.economy, r.first_year = s.first_year, r.growth = s.growth;
    if (standardise) {
      r.centre = detail::median(seen);
      std::vector<double> dev;
      for (const double g : seen) dev.push_back(std::fabs(g - r.centre));
      r.scale = std::max(1e-6, 1.4826 * detail::median(dev));
    }
    for (const double g : s.growth) r.z.push_back(std::isnan(g) || !standardise ? g : std::clamp((g - r.centre) / r.scale, -6.0, 6.0));
    seqs.push_back(r.z);
    readings.push_back(std::move(r));
    R.observations += seen.size();
  }
  const auto F = tyche::processes::fit_hmm(seqs, K, 5, 500, 1e-9, 20261104, shared_sd);
  R.model = F.model, R.log_likelihood = F.log_likelihood, R.iterations = F.iterations, R.converged = F.converged;
  for (auto& r : readings) {
    const auto D = tyche::processes::decode(R.model, r.z);
    r.filtered = D.filtered, r.smoothed = D.smoothed, r.path = D.path, r.log_likelihood = D.log_likelihood;
  }
  R.readings = std::move(readings);
  return R;
}

}  // namespace economics
