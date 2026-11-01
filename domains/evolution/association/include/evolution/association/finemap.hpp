// Fine-mapping: which variants of an associated region carry its effects, as posterior inclusion
// probabilities and credible sets of variants.
//
// - One causal variant (Wakefield 2009, Genetic Epidemiology 33, 79–86; Maller et al. 2012, Nature
//   Genetics 44, 1294–1301). From an estimate β̂ with variance V and a prior β ~ N(0, W), the approximate
//   Bayes factor of association against none is ABF = √(V/(V + W)) exp(z² W/(2(V + W))), z = β̂/√V
//   (the reciprocal of Wakefield's eq. 2, which is H₀ against H₁). If exactly one variant of the region
//   is causal, with a uniform prior, its posterior is ABF_j/Σ_k ABF_k, and the 95% credible set is the
//   smallest set of variants whose posteriors sum to 0.95. W = (0.15 sd y)² for a quantitative trait and
//   0.2² on the log-odds scale for a case–control one are coloc's defaults (Giambartolomei et al. 2014),
//   from Wakefield's choice of a prior that puts 95% of effects below a relative risk of about 1.5.
// - Several causal variants: SuSiE (tyche regression's susie, Wang et al. 2020) on the region's
//   individual dosages with the trait, both residualised on the covariates first (the practice of
//   susieR's documentation for covariates); or SuSiE-RSS (tyche's susie_rss, Zou et al. 2022) on the
//   region's z-scores and the LD correlation of a panel. The sets of columns become credible sets of
//   variants (site indices), each holding one effect with 95% probability.
// - Against the truth (a simulation's causal variants): a set covers when it holds a causal variant;
//   the coverage is the fraction of sets that cover, the power the fraction of causal variants in some
//   set, as Wang et al. (2020, §5) measure them.
//
// COMPLEXITY. ABF O(p); SuSiE O(n p²) to form XᵀX, then O(L p²) an iteration (tyche).

#pragma once

#include <evolution/association/ld.hpp>
#include <evolution/association/scan.hpp>
#include <evolution/population/genotypes.hpp>

#include <tyche/regression/susie.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

namespace evolution::association {

// log ABF of association against none (natural logarithm).
inline double log_abf(double beta, double variance, double prior_variance) {
  if (!(variance > 0) || !(prior_variance > 0)) return 0;
  const double z2 = beta * beta / variance, r = prior_variance / (variance + prior_variance);
  return 0.5 * std::log1p(-r) + 0.5 * z2 * r;
}

struct SingleCausal {
  std::vector<double> log_abf;      // per variant
  std::vector<double> posterior;    // per variant, summing to 1
  std::vector<std::size_t> credible;  // the smallest set reaching `coverage`, by decreasing posterior
  double coverage = 0;              // its posterior mass
};

inline SingleCausal finemap_abf(std::span<const double> beta, std::span<const double> se, double prior_sd, double coverage = 0.95) {
  const std::size_t p = beta.size();
  if (se.size() != p || p == 0) throw std::invalid_argument("evolution::association::finemap_abf: one standard error per estimate");
  SingleCausal f;
  f.log_abf.resize(p);
  double top = -INFINITY;
  for (std::size_t j = 0; j < p; j++) top = std::max(top, f.log_abf[j] = log_abf(beta[j], se[j] * se[j], prior_sd * prior_sd));
  double z = 0;
  for (const double l : f.log_abf) z += std::exp(l - top);
  f.posterior.resize(p);
  for (std::size_t j = 0; j < p; j++) f.posterior[j] = std::exp(f.log_abf[j] - top) / z;
  std::vector<std::size_t> o(p);
  std::iota(o.begin(), o.end(), 0);
  std::stable_sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return f.posterior[a] > f.posterior[b]; });
  for (const std::size_t j : o) {
    f.credible.push_back(j);
    f.coverage += f.posterior[j];
    if (f.coverage >= coverage) break;
  }
  return f;
}

// SuSiE on a region's individual data: the trait (one value per listed sample) and the sites listed,
// both residualised on C. The fit's columns are the sites in the order given.
inline tyche::regression::SusieFit finemap_susie(const population::Genotypes& g, std::span<const std::size_t> samples, std::span<const std::size_t> sites,
                                                 std::span<const double> y, const Rectangular& C, const tyche::regression::SusieOptions& o = {},
                                                 calx::scheduler::ThreadPool* pool = nullptr) {
  const std::size_t n = samples.size(), p = sites.size();
  if (y.size() != n || C.rows != n || p == 0) throw std::invalid_argument("evolution::association::finemap_susie: one y and one row of C per sample, some sites");
  const Residualizer M(C);
  Rectangular X(n, p);
  for (std::size_t c = 0; c < p; c++) {
    auto d = sample_dosages(g, sites[c], samples);
    M.apply(d);
    for (std::size_t i = 0; i < n; i++) X(i, c) = d[i];
  }
  const auto yt = M.residuals(y);
  return tyche::regression::susie(X, yt, o, pool);
}

// SuSiE-RSS on a region's z-scores, the LD correlation R of a panel and the sample size n (σ² fixed, as
// susieR's susie_rss has it).
inline tyche::regression::SusieFit finemap_rss(std::span<const double> z, const calx::linalg::Matrix& R, std::size_t n, tyche::regression::SusieOptions o = {}) {
  o.estimate_residual_variance = false;
  return tyche::regression::susie_rss(z, R, n, o);
}

struct FinemapTruth {
  std::size_t sets = 0, covering = 0, causal = 0, found = 0;
  double coverage() const { return sets ? static_cast<double>(covering) / static_cast<double>(sets) : 0; }
  double power() const { return causal ? static_cast<double>(found) / static_cast<double>(causal) : 0; }
};

// Credible sets (positions into `sites`) against the causal sites among them.
inline FinemapTruth against_truth(std::span<const std::vector<std::size_t>> sets, std::span<const std::size_t> sites, std::span<const std::size_t> causal) {
  FinemapTruth t;
  const auto is_causal = [&](std::size_t site) { return std::find(causal.begin(), causal.end(), site) != causal.end(); };
  for (const auto& s : sets) {
    t.sets++;
    bool hit = false;
    for (const std::size_t k : s) hit = hit || is_causal(sites[k]);
    t.covering += hit;
  }
  for (const std::size_t c : causal) {
    if (std::find(sites.begin(), sites.end(), c) == sites.end()) continue;
    t.causal++;
    bool in = false;
    for (const auto& s : sets)
      for (const std::size_t k : s) in = in || sites[k] == c;
    t.found += in;
  }
  return t;
}

}  // namespace evolution::association
