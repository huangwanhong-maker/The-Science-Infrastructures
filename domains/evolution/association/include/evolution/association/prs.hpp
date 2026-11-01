// Polygenic scores: a trait predicted from many variants' estimated effects, and their accuracy in
// samples the effects were not estimated from.
//
// - Clumping and thresholding (C+T; Purcell et al. 2009, Nature 460, 748–752; Choi, Mak and O'Reilly
//   2020, Nature Protocols 15, 2759–2772): the scan's sites are clumped (ld.hpp, r² > 0.1 within 250 kb,
//   p1 = p2 = 1), and for each threshold p_T the score is Σ_j β̂_j g_ij over the index sites with p < p_T.
//   The threshold is chosen on a validation sample and the accuracy reported on a test sample apart.
// - LDpred-inf (Vilhjálmsson et al. 2015, AJHG 97, 576–592, "LDpred-inf"): under the infinitesimal model
//   (every standardised effect N(0, h²/M)), the posterior mean of the joint standardised effects given
//   the marginal ones β̃ = Xᵀy/(N − 1) (standardised X and y) is E(β | β̃, D) ≈ (D + (M/(N h²)) I)⁻¹ β̃ with
//   D the LD correlation, computed within blocks of consecutive sites (their LD windows; D is
//   block-diagonal beyond them). With the training sample's own D this is ridge regression, the BLUP
//   (checks). Effects are returned per allele (divided by the sites' standard deviations).
// - Accuracy: the squared correlation of the score with the trait adjusted for covariates in the test
//   sample. Daetwyler, Villanueva and Woolliams (2008, PLoS ONE 3, e3395) expect about
//   h² · N h²/(N h² + M_e) for M_e effectively independent segments, a guide, not a check.
//
// COMPLEXITY. C+T as clumping; LDpred-inf O(N m b) for blocks of b sites plus O(m b²).

#pragma once

#include <evolution/association/ld.hpp>
#include <evolution/association/scan.hpp>
#include <evolution/population/genotypes.hpp>

#include <calx/linalg/cholesky.hpp>
#include <calx/linalg/dense.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace evolution::association {

// Σ_j w_j g_ij for the samples (per-allele weights at the sites; missing calls at the site's mean).
inline std::vector<double> polygenic_score(const population::Genotypes& g, std::span<const std::size_t> samples, std::span<const std::size_t> sites,
                                           std::span<const double> weights) {
  if (weights.size() != sites.size()) throw std::invalid_argument("evolution::association::polygenic_score: one weight per site");
  std::vector<double> s(samples.size(), 0.0);
  for (std::size_t c = 0; c < sites.size(); c++) {
    if (weights[c] == 0) continue;
    const auto d = sample_dosages(g, sites[c], samples);
    for (std::size_t i = 0; i < d.size(); i++) s[i] += weights[c] * d[i];
  }
  return s;
}

// The squared correlation of a score with y, both adjusted for C.
inline double score_r2(std::span<const double> score, std::span<const double> y, const Rectangular& C) {
  const Residualizer M(C);
  const auto a = M.residuals(score), b = M.residuals(y);
  double ab = 0, aa = 0, bb = 0;
  for (std::size_t i = 0; i < a.size(); i++) {
    ab += a[i] * b[i];
    aa += a[i] * a[i];
    bb += b[i] * b[i];
  }
  return aa > 0 && bb > 0 ? ab * ab / (aa * bb) : 0;
}

struct ThresholdScore {
  double threshold = 1;
  std::vector<std::size_t> sites;   // the index sites below it
  std::vector<double> weights;      // their β̂
};

// C+T: the clumped index sites of a scan (LD from `ld_samples`) at each threshold.
inline std::vector<ThresholdScore> clump_and_threshold(const population::Genotypes& g, std::span<const std::size_t> ld_samples, const Scan& scan,
                                                       std::span<const double> thresholds, ClumpOptions o = {1.0, 1.0, 0.1, 250'000}) {
  std::vector<double> p = scan.p;
  for (double& v : p)
    if (!(v >= 0 && v <= 1)) v = 1;
  o.p1 = std::nextafter(o.p1, 2.0);
  o.p2 = std::nextafter(o.p2, 2.0);
  const auto clumps = clump(g, ld_samples, p, o);
  std::vector<ThresholdScore> out;
  for (const double t : thresholds) {
    ThresholdScore s;
    s.threshold = t;
    for (const auto& c : clumps)
      if (p[c.index] < t && std::isfinite(scan.beta[c.index])) {
        s.sites.push_back(c.index);
        s.weights.push_back(scan.beta[c.index]);
      }
    std::vector<std::size_t> order(s.sites.size());
    for (std::size_t k = 0; k < order.size(); k++) order[k] = k;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return s.sites[a] < s.sites[b]; });
    ThresholdScore sorted;
    sorted.threshold = t;
    for (const std::size_t k : order) {
      sorted.sites.push_back(s.sites[k]);
      sorted.weights.push_back(s.weights[k]);
    }
    out.push_back(std::move(sorted));
  }
  return out;
}

struct LdpredInf {
  std::vector<std::size_t> sites;
  std::vector<double> standardised;  // the posterior mean effects on standardised genotypes
  std::vector<double> weights;       // per allele
  double shrinkage = 0;              // M/(N h²)
};

// LDpred-inf from the training samples' trait (adjusted for C) and dosages, in blocks of `block` sites
// of one chromosome each.
inline LdpredInf ldpred_inf(const population::Genotypes& g, std::span<const std::size_t> samples, std::span<const double> y, const Rectangular& C,
                            double h2, std::size_t block = 200, calx::scheduler::ThreadPool* pool = nullptr) {
  const std::size_t n = samples.size(), m = g.sites;
  if (y.size() != n || C.rows != n) throw std::invalid_argument("evolution::association::ldpred_inf: one y and one row of C per sample");
  if (!(h2 > 0 && h2 <= 1) || block == 0) throw std::invalid_argument("evolution::association::ldpred_inf: 0 < h² ≤ 1 and blocks of sites");
  const Residualizer M(C);
  auto yt = M.residuals(y);
  {
    double ss = 0;
    for (const double v : yt) ss += v * v;
    const double sd = std::sqrt(ss / static_cast<double>(n - 1));
    for (double& v : yt) v /= sd;
  }
  LdpredInf out;
  out.shrinkage = static_cast<double>(m) / (static_cast<double>(n) * h2);
  out.sites.resize(m);
  for (std::size_t j = 0; j < m; j++) out.sites[j] = j;
  out.standardised.assign(m, 0.0);
  out.weights.assign(m, 0.0);
  // Blocks of consecutive sites within a chromosome.
  std::vector<std::pair<std::size_t, std::size_t>> blocks;
  for (std::size_t c = 0; c < g.chromosomes.size(); c++) {
    const auto [lo, hi] = g.chromosome_sites(static_cast<std::uint32_t>(c));
    for (std::size_t b = lo; b < hi; b += block) blocks.emplace_back(b, std::min(hi, b + block));
  }
  std::vector<double> sd(m, 0.0);
  const auto body = [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t q = begin; q < end; q++) {
      const auto [lo, hi] = blocks[q];
      const std::size_t b = hi - lo;
      std::vector<std::size_t> idx(b);
      for (std::size_t k = 0; k < b; k++) idx[k] = lo + k;
      Rectangular Z = standardised_dosages(g, samples, idx);
      for (std::size_t k = 0; k < b; k++) {
        double mean = 0;
        const auto d = sample_dosages(g, idx[k], samples, &mean);
        double ss = 0;
        for (const double v : d) ss += (v - mean) * (v - mean);
        sd[idx[k]] = std::sqrt(ss / static_cast<double>(n - 1));
      }
      calx::linalg::Matrix D(b);
      std::vector<double> marginal(b, 0.0);
      for (std::size_t a = 0; a < b; a++) {
        for (std::size_t i = 0; i < n; i++) marginal[a] += Z(i, a) * yt[i];
        marginal[a] /= static_cast<double>(n - 1);
        for (std::size_t c = 0; c <= a; c++) {
          double s = 0;
          for (std::size_t i = 0; i < n; i++) s += Z(i, a) * Z(i, c);
          D(a, c) = D(c, a) = s / static_cast<double>(n - 1);
        }
        D(a, a) += out.shrinkage;
      }
      const auto ch = calx::linalg::cholesky(D);
      if (!ch.positive_definite) continue;
      const auto beta = calx::linalg::cholesky_solve(ch, marginal);
      for (std::size_t k = 0; k < b; k++) {
        out.standardised[idx[k]] = beta[k];
        out.weights[idx[k]] = sd[idx[k]] > 0 ? beta[k] / sd[idx[k]] : 0;
      }
    }
  };
  if (pool) calx::scheduler::parallel_for(*pool, blocks.size(), 1, body);
  else body(0, blocks.size(), 0);
  return out;
}

}  // namespace evolution::association
