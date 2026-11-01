// Expression predicted from genotypes, and its association with a trait: PrediXcan, S-PrediXcan and
// the transcriptome-wide association study (TWAS).
//
// - PrediXcan (Gamazon et al. 2015, Nature Genetics 47, 1091–1098). For each gene, its expression in a
//   reference panel (adjusted for covariates) is regressed on the dosages of the variants within a
//   window of its TSS by the elastic net with α = 0.5, λ chosen by 10-fold cross-validation (tyche
//   regression's cross_validate, glmnet's algorithm; Friedman, Hastie and Tibshirani 2010; a path of 20
//   penalties down to 0.01 λ_max by default, glmnet's ratio for n < p, against PrediXcan's 100). The model's
//   accuracy is estimated by nested cross-validation, as PredictDB does (Barbeira et al. 2018, Nature
//   Communications 9, 1825, Methods): in each of 5 outer folds the whole inner procedure is run on the
//   other four and the held-out fold predicted; ρ̄ is the mean of the folds' correlations of predicted
//   and observed expression, R² the held-out coefficient of determination, and the folds' correlations
//   are combined as Stouffer's z of their Fisher transforms, z_k = atanh(ρ_k)√(n_k − 3); a gene's model
//   is kept when ρ̄ > 0.1 and that z's two-sided p < 0.05 (PredictDB's criteria). The final weights are
//   those of the inner procedure on the whole panel, per allele.
// - TWAS from individual data: the predicted expression Ê_g = Σ_l w_l g_l of each sample of the
//   association study is regressed on the trait with the covariates (scan.hpp's Frisch–Waugh–Lovell
//   regression), its t and p the gene's association.
// - S-PrediXcan (Barbeira et al. 2018, eq. 2): from the GWAS's summary statistics alone,
//   z_g ≈ Σ_l w_l (σ̂_l/σ̂_g)(β̂_l/se_l), σ̂_l the sd of variant l's dosage and σ̂_g² = wᵀΓw the variance of
//   the predicted expression, Γ the covariance of the dosages in a reference panel. FUSION's statistic
//   (Gusev et al. 2016, Nature Genetics 48, 245–252, its eq.: z = w̃ᵀz/√(w̃ᵀRw̃) with w̃ = σ̂ ∘ w and R the
//   LD correlation) is the same number, since Γ = diag σ̂ · R · diag σ̂. With in-sample LD it equals the
//   individual-level z up to the approximations of β̂/se ≈ the marginal correlation's z (checks).
//
// COMPLEXITY. Per gene, O(F · I · L · n p) coordinate descent passes for F outer and I inner folds and a
// path of L penalties (warm-started), genes in parallel; prediction O(n p); S-PrediXcan O(p²).

#pragma once

#include <evolution/association/ld.hpp>
#include <evolution/association/scan.hpp>
#include <evolution/population/genotypes.hpp>

#include <tyche/probability/distributions.hpp>
#include <tyche/regression/penalized.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace evolution::association {

struct PredictOptions {
  double alpha = 0.5;               // the elastic net's mixing (PrediXcan's)
  std::size_t folds = 10, outer_folds = 5;
  std::uint64_t window = 1'000'000; // ± of the TSS
  std::size_t lambda_count = 20;
  double lambda_min_ratio = 0.01;   // λ_min/λ_max: glmnet's ratio when n < p, as in a cis window
  double tolerance = 1e-7;          // glmnet's default threshold
  double min_rho = 0.1, max_p = 0.05;  // PredictDB's criteria for keeping a model
  std::uint64_t seed = 1;
};

struct GeneModel {
  std::size_t gene = 0;
  std::vector<std::size_t> sites;   // the window's sites
  std::vector<double> weights;      // per allele, one per site (mostly zero)
  double intercept = 0;
  double lambda = 0;
  std::size_t nonzero = 0;
  double rho_mean = 0, r2_nested = 0, z = 0, p = 1;  // nested cross-validation
  bool kept = false;
};

namespace detail {

inline double correlation(std::span<const double> a, std::span<const double> b) {
  const std::size_t n = a.size();
  double ma = 0, mb = 0;
  for (std::size_t i = 0; i < n; i++) {
    ma += a[i];
    mb += b[i];
  }
  ma /= static_cast<double>(n);
  mb /= static_cast<double>(n);
  double ab = 0, aa = 0, bb = 0;
  for (std::size_t i = 0; i < n; i++) {
    ab += (a[i] - ma) * (b[i] - mb);
    aa += (a[i] - ma) * (a[i] - ma);
    bb += (b[i] - mb) * (b[i] - mb);
  }
  return aa > 0 && bb > 0 ? ab / std::sqrt(aa * bb) : 0;
}

inline Rectangular rows_of(const Rectangular& X, std::span<const std::size_t> rows) {
  Rectangular out(rows.size(), X.cols);
  for (std::size_t k = 0; k < rows.size(); k++)
    for (std::size_t j = 0; j < X.cols; j++) out(k, j) = X(rows[k], j);
  return out;
}

// The inner procedure: cv.glmnet's λ_min fit.
inline tyche::regression::PenalizedFit fit_inner(const Rectangular& X, std::span<const double> y, const PredictOptions& o, std::uint64_t seed) {
  tyche::regression::PenalizedOptions po;
  po.alpha = o.alpha;
  po.lambda_count = o.lambda_count;
  po.lambda_min_ratio = o.lambda_min_ratio;
  po.tolerance = o.tolerance;
  tyche::regression::CvOptions cv;
  cv.folds = std::min(o.folds, y.size());
  cv.seed = seed;
  const auto c = tyche::regression::cross_validate(X, y, po, cv);
  return c.path.fits[c.index_min];
}

}  // namespace detail

// The PrediXcan model of one gene: the window's dosages X (samples × sites, per allele) and its
// expression y (adjusted for covariates).
inline GeneModel train_expression_model(const Rectangular& X, std::span<const double> y, const PredictOptions& o = {}) {
  const std::size_t n = X.rows, p = X.cols;
  GeneModel m;
  m.weights.assign(p, 0.0);
  if (p == 0 || n < 2 * o.outer_folds) return m;
  // Nested cross-validation.
  const auto outer = tyche::regression::make_folds(n, o.outer_folds, o.seed);
  std::vector<double> predicted(n, 0.0);
  double z_sum = 0, rho_sum = 0;
  for (std::size_t f = 0; f < o.outer_folds; f++) {
    std::vector<std::size_t> train, test;
    for (std::size_t i = 0; i < n; i++) (outer[i] == f ? test : train).push_back(i);
    std::vector<double> yt(train.size()), ytest(test.size());
    for (std::size_t k = 0; k < train.size(); k++) yt[k] = y[train[k]];
    for (std::size_t k = 0; k < test.size(); k++) ytest[k] = y[test[k]];
    const auto fit = detail::fit_inner(detail::rows_of(X, train), yt, o, o.seed + 1 + f);
    const auto pr = tyche::regression::predict(fit, detail::rows_of(X, test));
    for (std::size_t k = 0; k < test.size(); k++) predicted[test[k]] = pr[k];
    const double rho = detail::correlation(pr, ytest);
    rho_sum += rho;
    z_sum += std::atanh(std::clamp(rho, -0.999999, 0.999999)) * std::sqrt(std::max(1.0, static_cast<double>(test.size()) - 3));
  }
  m.rho_mean = rho_sum / static_cast<double>(o.outer_folds);
  m.z = z_sum / std::sqrt(static_cast<double>(o.outer_folds));
  m.p = std::min(1.0, 2 * tyche::probability::Normal{}.survival(std::fabs(m.z)));
  double my = 0;
  for (const double v : y) my += v;
  my /= static_cast<double>(n);
  double rss = 0, tss = 0;
  for (std::size_t i = 0; i < n; i++) {
    rss += (y[i] - predicted[i]) * (y[i] - predicted[i]);
    tss += (y[i] - my) * (y[i] - my);
  }
  m.r2_nested = tss > 0 ? 1 - rss / tss : 0;
  // The final model on the whole panel.
  const auto fit = detail::fit_inner(X, y, o, o.seed);
  m.weights = fit.beta;
  m.intercept = fit.intercept;
  m.lambda = fit.lambda;
  for (const double w : m.weights) m.nonzero += w != 0;
  m.kept = m.nonzero > 0 && m.rho_mean > o.min_rho && m.p < o.max_p;
  return m;
}

// PrediXcan models of every gene from the reference panel's genotypes and expression (each gene's
// values one per listed sample), adjusted for C; in parallel over genes.
template <class Locus>
inline std::vector<GeneModel> train_expression_models(const population::Genotypes& g, std::span<const std::size_t> samples,
                                                      const std::vector<std::vector<double>>& expression, std::span<const Locus> loci, const Rectangular& C,
                                                      const PredictOptions& o = {}, calx::scheduler::ThreadPool* pool = nullptr) {
  const std::size_t G = expression.size(), n = samples.size();
  if (loci.size() != G || C.rows != n) throw std::invalid_argument("evolution::association::train_expression_models: one locus per gene, one row of C per sample");
  const Residualizer M(C);
  std::vector<GeneModel> out(G);
  const auto body = [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t gi = begin; gi < end; gi++) {
      const std::uint64_t tss = loci[gi].tss, from = tss > o.window ? tss - o.window : 0;
      const auto sites = sites_between(g, loci[gi].chromosome, from, tss + o.window + 1);
      std::vector<std::size_t> kept;
      for (const std::size_t j : sites) {
        const auto [alt, total] = population::count_site(g.site(j), samples);
        if (alt > 0 && alt < total) kept.push_back(j);
      }
      Rectangular X(n, kept.size());
      for (std::size_t c = 0; c < kept.size(); c++) {
        const auto d = sample_dosages(g, kept[c], samples);
        for (std::size_t i = 0; i < n; i++) X(i, c) = d[i];
      }
      const auto y = M.residuals(expression[gi]);
      PredictOptions og = o;
      og.seed = o.seed + 7919 * gi;
      GeneModel m = train_expression_model(X, y, og);
      m.gene = gi;
      m.sites = std::move(kept);
      out[gi] = std::move(m);
    }
  };
  if (pool) calx::scheduler::parallel_for(*pool, G, 1, body);
  else body(0, G, 0);
  return out;
}

// The predicted expression Ê = intercept + Σ w g of the samples listed.
inline std::vector<double> predict_expression(const GeneModel& m, const population::Genotypes& g, std::span<const std::size_t> samples) {
  std::vector<double> e(samples.size(), m.intercept);
  for (std::size_t c = 0; c < m.sites.size(); c++) {
    if (m.weights[c] == 0) continue;
    const auto d = sample_dosages(g, m.sites[c], samples);
    for (std::size_t i = 0; i < d.size(); i++) e[i] += m.weights[c] * d[i];
  }
  return e;
}

struct Association {
  double beta = 0, se = 0, z = 0, p = 1;
};

// The trait (one per sample) regressed on a predictor with the covariates C (the TWAS of individual data).
inline Association associate(std::span<const double> x, std::span<const double> y, const Rectangular& C) {
  const std::size_t n = x.size();
  if (y.size() != n || C.rows != n) throw std::invalid_argument("evolution::association::associate: one x, y and row of C per sample");
  const Residualizer M(C);
  const auto xt = M.residuals(x), yt = M.residuals(y);
  double xx = 0, xy = 0, yy = 0;
  for (std::size_t i = 0; i < n; i++) {
    xx += xt[i] * xt[i];
    xy += xt[i] * yt[i];
    yy += yt[i] * yt[i];
  }
  Association a;
  if (!(xx > 0)) return a;
  const double df = static_cast<double>(n - M.rank() - 1);
  a.beta = xy / xx;
  a.se = std::sqrt(std::max(0.0, yy - a.beta * xy) / df / xx);
  a.z = a.beta / a.se;
  a.p = std::min(1.0, 2 * tyche::probability::StudentT{df}.survival(std::fabs(a.z)));
  return a;
}

// S-PrediXcan's z from a GWAS's per-site z-scores (β̂/se, one per site of g) and the LD of a panel.
inline Association s_predixcan(const GeneModel& m, std::span<const double> z, const population::Genotypes& g, std::span<const std::size_t> panel) {
  std::vector<std::size_t> used;
  std::vector<double> w;
  for (std::size_t c = 0; c < m.sites.size(); c++)
    if (m.weights[c] != 0 && std::isfinite(z[m.sites[c]])) {
      used.push_back(m.sites[c]);
      w.push_back(m.weights[c]);
    }
  Association a;
  if (used.empty()) return a;
  const auto R = ld_correlation(g, panel, used);
  std::vector<double> sd(used.size());
  for (std::size_t k = 0; k < used.size(); k++) {
    double mean = 0;
    const auto d = sample_dosages(g, used[k], panel, &mean);
    double ss = 0;
    for (const double v : d) ss += (v - mean) * (v - mean);
    sd[k] = std::sqrt(ss / static_cast<double>(panel.size() - 1));
  }
  // σ̂_g² = wᵀΓw = (w∘σ)ᵀ R (w∘σ); z_g = Σ w_l σ_l z_l / σ̂_g.
  double num = 0, var = 0;
  for (std::size_t a1 = 0; a1 < used.size(); a1++) {
    num += w[a1] * sd[a1] * z[used[a1]];
    for (std::size_t b = 0; b < used.size(); b++) var += w[a1] * sd[a1] * R(a1, b) * w[b] * sd[b];
  }
  if (!(var > 0)) return a;
  a.z = num / std::sqrt(var);
  a.p = std::min(1.0, 2 * tyche::probability::Normal{}.survival(std::fabs(a.z)));
  return a;
}

}  // namespace evolution::association
