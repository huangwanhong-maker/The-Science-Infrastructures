// Colocalisation: whether two traits measured at one region (a disease and a gene's expression, say)
// share a causal variant.
//
// THE MODEL (Giambartolomei et al. 2014, PLoS Genetics 10, e1004383, the method of R's coloc.abf). Each
// trait has at most one causal variant in the region; a configuration is a pair (j, k) of the two traits'
// causal variants (or none). Five hypotheses partition the configurations: H₀ neither trait has one; H₁
// only the first; H₂ only the second; H₃ both, at different variants; H₄ both, at the same one. With
// per-variant priors p₁, p₂ (one trait) and p₁₂ (both at that variant), and Wakefield's approximate Bayes
// factors ABF¹_j and ABF²_k of each trait (finemap.hpp), the posterior of each hypothesis is proportional to
//   H₀: 1,  H₁: p₁ Σ_j ABF¹_j,  H₂: p₂ Σ_k ABF²_k,
//   H₃: p₁p₂ Σ_{j≠k} ABF¹_j ABF²_k = p₁p₂ (Σ_j ABF¹_j Σ_k ABF²_k − Σ_j ABF¹_j ABF²_j),
//   H₄: p₁₂ Σ_j ABF¹_j ABF²_j
// (their eqs. 1–4 and Table 1, in coloc's combine.abf form, computed on the log scale with log-sum-exp
// and log-difference), and the posterior that variant j is the shared one, given H₄, is
// ABF¹_j ABF²_j/Σ_k ABF¹_k ABF²_k. The defaults p₁ = p₂ = 10⁻⁴, p₁₂ = 10⁻⁵ are coloc's.
//
// SEVERAL CAUSAL VARIANTS (Wallace 2021, PLoS Genetics 17, e1009440, coloc.susie): SuSiE decomposes each
// trait's signal into single effects; each pair of effects (one of each trait) is colocalised as above
// from their per-variant log Bayes factors (the SER's lbf_variable), relaxing the one-variant assumption.
//
// COMPLEXITY. O(p) for p variants; coloc.susie O(L₁ L₂ p).

#pragma once

#include <evolution/association/finemap.hpp>

#include <tyche/regression/susie.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace evolution::association {

struct ColocPriors {
  double p1 = 1e-4, p2 = 1e-4, p12 = 1e-5;
};

struct Coloc {
  std::array<double, 5> posterior{};      // PP.H0 … PP.H4
  std::array<double, 5> log_bayes{};      // the unnormalised log posteriors (lH0 = 0)
  std::vector<double> shared;             // per variant: P(it is the shared one | H₄)
  std::size_t variants = 0;
};

namespace detail {

inline double log_sum(std::span<const double> x) {
  double top = -INFINITY;
  for (const double v : x) top = std::max(top, v);
  if (!std::isfinite(top)) return top;
  double s = 0;
  for (const double v : x) s += std::exp(v - top);
  return top + std::log(s);
}

// log(eᵃ − eᵇ), a > b.
inline double log_diff(double a, double b) { return a + std::log1p(-std::exp(b - a)); }

}  // namespace detail

// Colocalisation from the two traits' per-variant log Bayes factors.
inline Coloc coloc_from_lbf(std::span<const double> l1, std::span<const double> l2, const ColocPriors& pr = {}) {
  const std::size_t p = l1.size();
  if (l2.size() != p || p == 0) throw std::invalid_argument("evolution::association::coloc: both traits at the same variants");
  std::vector<double> sum(p);
  for (std::size_t j = 0; j < p; j++) sum[j] = l1[j] + l2[j];
  const double s1 = detail::log_sum(l1), s2 = detail::log_sum(l2), s12 = detail::log_sum(sum);
  Coloc c;
  c.variants = p;
  c.log_bayes = {0.0, std::log(pr.p1) + s1, std::log(pr.p2) + s2,
                 std::log(pr.p1) + std::log(pr.p2) + (s1 + s2 > s12 ? detail::log_diff(s1 + s2, s12) : -INFINITY), std::log(pr.p12) + s12};
  const double denom = detail::log_sum(c.log_bayes);
  for (std::size_t h = 0; h < 5; h++) c.posterior[h] = std::exp(c.log_bayes[h] - denom);
  c.shared.resize(p);
  for (std::size_t j = 0; j < p; j++) c.shared[j] = std::exp(sum[j] - s12);
  return c;
}

// Colocalisation from the two traits' estimates and standard errors at the same variants, with the
// priors' standard deviations of an effect (0.15 sd y for a quantitative trait, 0.2 for a case–control one).
inline Coloc coloc_abf(std::span<const double> beta1, std::span<const double> se1, double prior_sd1, std::span<const double> beta2,
                       std::span<const double> se2, double prior_sd2, const ColocPriors& pr = {}) {
  const std::size_t p = beta1.size();
  if (se1.size() != p || beta2.size() != p || se2.size() != p) throw std::invalid_argument("evolution::association::coloc_abf: four vectors of one length");
  std::vector<double> l1(p), l2(p);
  for (std::size_t j = 0; j < p; j++) {
    l1[j] = log_abf(beta1[j], se1[j] * se1[j], prior_sd1 * prior_sd1);
    l2[j] = log_abf(beta2[j], se2[j] * se2[j], prior_sd2 * prior_sd2);
  }
  return coloc_from_lbf(l1, l2, pr);
}

struct ColocPair {
  std::size_t effect1 = 0, effect2 = 0;   // the SuSiE effects
  Coloc coloc;
};

// coloc.susie: every pair of the two fits' kept effects (σ₀² above the fits' tolerance) with a credible set.
inline std::vector<ColocPair> coloc_susie(const tyche::regression::SusieFit& a, const tyche::regression::SusieFit& b, const ColocPriors& pr = {}) {
  if (a.p != b.p) throw std::invalid_argument("evolution::association::coloc_susie: two fits on the same variants");
  std::vector<ColocPair> out;
  for (const auto& sa : a.sets)
    for (const auto& sb : b.sets) {
      ColocPair c;
      c.effect1 = sa.effect;
      c.effect2 = sb.effect;
      c.coloc = coloc_from_lbf(std::span<const double>(&a.lbf_variable[sa.effect * a.p], a.p), std::span<const double>(&b.lbf_variable[sb.effect * b.p], b.p), pr);
      out.push_back(std::move(c));
    }
  return out;
}

}  // namespace evolution::association
