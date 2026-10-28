// Genotype imputation: the alleles of a target at the sites it was not typed at, from a
// reference panel typed at all of them, as the posterior predictive of the copying model.
//
// The untyped sites of the target are missing data (−1) of the copying chain (copying.hpp),
// emitted with probability 1, so the chain runs over every site of the panel and its posterior
// at an untyped site j carries the evidence of the typed sites on both sides. Given the copied
// haplotype the allele is independent of everything else, so
//
//   P(x_j = 1 | the typed alleles) = Σ_k γ_j(k) P(1 | h_k(j)),   P(1 | h) = 1 − ε if h = 1, else ε,
//
// the dosage of a haplotype; and for a diploid's unphased genotypes, from the pair model's
// posterior (phasing.hpp),
//
//   P(g_j = m | the typed genotypes) = Σ_{a,b} γ_j(a, b) P(m | a, b),   m = 0, 1, 2,
//
// with the dosage E[g_j] = Σ_m m P(g_j = m). Marchini et al. (2007, Nature Genetics 39,
// 906–913) imputed genotypes in this way, from a model over pairs of panel haplotypes. A typed
// site keeps its observed value.
//
// The quality of an imputed site across targets is reported, as is usual, by the squared
// correlation r² between the imputed dosages and the true alleles or genotypes (dosage_r2).
//
// Complexity: O(N K) for a haplotype (one forward–backward pass of the copying chain),
// O(N K²) time and O(N K²) memory for a diploid. The posterior predictive is exact under the
// model; its accuracy against the truth depends on how closely the panel is related to the
// target, which the checks measure on the coalescent's samples.

#pragma once

#include <evolution/haplotypes/copying.hpp>
#include <evolution/haplotypes/phasing.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace evolution::haplotypes {

struct Imputation {
  std::vector<double> dosage;    // E[x_j] (a haplotype) or E[g_j] (a diploid) at each site
  std::vector<double> genotype;  // a diploid's P(g_j = 0, 1, 2), sites × 3; empty for a haplotype
  double log_likelihood = NEG_INF;  // of the typed alleles
  double at(std::size_t j, int m) const { return genotype[3 * j + static_cast<std::size_t>(m)]; }
};

// The dosages of a haplotype typed at the sites where target_j ≥ 0.
inline Imputation impute_haplotype(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> target) {
  const Posterior p = posterior(panel, c, target);
  Imputation out;
  out.log_likelihood = p.log_likelihood;
  out.dosage.assign(c.sites, std::numeric_limits<double>::quiet_NaN());
  if (p.log_likelihood == NEG_INF) return out;
  const std::size_t K = c.states;
  for (std::size_t j = 0; j < c.sites; j++) {
    if (target[j] >= 0) {
      out.dosage[j] = target[j];
      continue;
    }
    const std::int8_t* h = panel.alleles.data() + j * K;
    double d = 0;
    for (std::size_t k = 0; k < K; k++) d += p.gamma[j * K + k] * (h[k] ? 1 - c.epsilon : c.epsilon);
    out.dosage[j] = d;
  }
  return out;
}

// The genotype probabilities and dosages of a diploid typed (unphased) where g_j ≥ 0.
inline Imputation impute_genotypes(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> g) {
  const std::size_t K = c.states, N = c.sites;
  Imputation out;
  out.dosage.assign(N, std::numeric_limits<double>::quiet_NaN());
  out.genotype.assign(3 * N, std::numeric_limits<double>::quiet_NaN());
  std::vector<double> p1(K);
  out.log_likelihood = for_each_pair_posterior(panel, c, g, [&](std::size_t j, std::span<const double> gamma) {
    double* q = &out.genotype[3 * j];
    if (g[j] >= 0) {
      q[0] = q[1] = q[2] = 0;
      q[g[j]] = 1;
      out.dosage[j] = g[j];
      return;
    }
    detail::allele_one(panel, j, c.epsilon, p1.data());
    // Σ γ(a, b) P(m | a, b) with P(1 | a) = p_a: the marginals of the two copies suffice for E[g],
    // the pairs for the genotype law.
    double g0 = 0, g1 = 0, g2 = 0;
    for (std::size_t a = 0; a < K; a++)
      for (std::size_t b = 0; b < K; b++) {
        const double w = gamma[a * K + b], p = p1[a], r = p1[b];
        g0 += w * (1 - p) * (1 - r);
        g1 += w * (p * (1 - r) + (1 - p) * r);
        g2 += w * p * r;
      }
    q[0] = g0;
    q[1] = g1;
    q[2] = g2;
    out.dosage[j] = g1 + 2 * g2;
  });
  return out;
}

// The squared Pearson correlation of imputed dosages with the true values across targets at one
// site; NaN when either is constant.
inline double dosage_r2(std::span<const double> dosage, std::span<const double> truth) {
  const std::size_t n = std::min(dosage.size(), truth.size());
  double sx = 0, sy = 0;
  for (std::size_t i = 0; i < n; i++) {
    sx += dosage[i];
    sy += truth[i];
  }
  const double mx = sx / static_cast<double>(n), my = sy / static_cast<double>(n);
  double sxx = 0, syy = 0, sxy = 0;
  for (std::size_t i = 0; i < n; i++) {
    const double dx = dosage[i] - mx, dy = truth[i] - my;
    sxx += dx * dx;
    syy += dy * dy;
    sxy += dx * dy;
  }
  if (!(sxx > 0) || !(syy > 0)) return std::numeric_limits<double>::quiet_NaN();
  return sxy * sxy / (sxx * syy);
}

}  // namespace evolution::haplotypes
