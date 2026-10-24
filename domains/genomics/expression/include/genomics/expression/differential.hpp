// Differential expression from counts: size factors by the median of ratios, a negative
// binomial generalised linear model for each gene fitted by IRLS, dispersions by Cox–Reid
// adjusted maximum likelihood shrunk toward a trend in the mean, Wald tests, and the false
// discovery rate of Benjamini and Hochberg. A documented simplification of DESeq2 (Love, Huber
// and Anders 2014, Genome Biology 15, 550, Methods).
//
// - Size factors (Anders and Huber 2010, Genome Biology 11, R106, eq. 5):
//     ŝ_j = median over genes i of k_ij / (Π_v k_iv)^{1/m},
//   over the genes with no zero count; a sample's depth relative to a pseudo-reference sample,
//   the geometric mean of all, robust to a minority of genes that change.
// - The model (Love et al. 2014, Methods): K_ij ~ NB(μ_ij, α_i), mean μ_ij = s_j q_ij and
//   variance μ + α μ², log q_ij = Σ_r x_jr β_ir. The negative binomial is the gamma–Poisson
//   mixture: λ ~ Γ(shape 1/α, scale α μ), K | λ ~ Poisson(λ) (theory note, Proposition 6.1).
// - IRLS (McCullagh and Nelder 1989, Generalized Linear Models, 2nd ed., §2.5): with the log
//   link and offset log s_j, the working weights w_j = μ_j / (1 + α μ_j) and working response
//   z_j = log(μ_j / s_j) + (y_j − μ_j)/μ_j, β ← (Xᵀ W X + λ I)⁻¹ Xᵀ W z until the deviance changes
//   by less than 10⁻⁸ of itself; λ = 10⁻⁶ and μ held at 0.5 or above in the weights, as DESeq2
//   does. Standard errors from (Xᵀ W X + λ I)⁻¹ at the fit. Linear solves by calx::linalg.
// - Dispersion (Cox and Reid 1987, J. R. Statist. Soc. B 49, 1–39; McCarthy, Chen and Smyth
//   2012, Nucleic Acids Research 40, 4288–4297): at the fitted means, the adjusted profile
//   log-likelihood ℓ(α) − ½ log det(Xᵀ W X) is maximised over log α in [log 10⁻⁸, log max(10, m)]:
//   a grid of 25 points, then Brent's method (calx::solvers) on its analytic derivative
//   (digamma, calx::special_functions) in the bracket of the grid's best point. The means are
//   refitted at the estimate and the dispersion estimated again (two passes; DESeq2 one).
// - The trend (Love et al. 2014, Methods): α_tr(μ̄) = a₀ + a₁/μ̄ in the mean of normalised counts,
//   fitted by a gamma-family GLM with identity link to the gene-wise estimates above 100 × 10⁻⁸,
//   iteratively dropping genes whose ratio to the trend leaves (10⁻⁴, 15), from the start
//   (0.1, 1) and for at most eleven fits, as DESeq2's parametric fit does. If a coefficient is
//   not positive, the trend is the mean estimate (DESeq2 falls back to a local fit).
// - Shrinkage (Love et al. 2014, Methods): a normal prior on log α centred on log α_tr, of
//   variance σ_d² = max(s²_lr − ψ′((m − p)/2), 0.25), s_lr the median absolute deviation (×
//   1.4826) of the log residuals of the genes in the fit, ψ′ the trigamma function; the
//   maximum a posteriori of the adjusted likelihood and the prior is the final dispersion,
//   except for genes whose own estimate lies more than 2 s_lr above the trend, which keep it.
//   *Simplification.* DESeq2's outlier replacement (Cook's distances), independent filtering
//   and the shrinkage of fold changes are not done; the trend has no local fallback.
// - The Wald test (Love et al. 2014, Methods): z = β̂_r / SE(β̂_r) against N(0, 1), two-sided,
//   P = Q(½, z²/2) (the χ²₁ tail, calx::special_functions). Fold changes are reported as log₂.
// - Benjamini and Hochberg (1995, J. R. Statist. Soc. B 57, 289–300, §3): with the m p-values
//   sorted, p_(1) ≤ … ≤ p_(m), the adjusted value of p_(i) is min_{k ≥ i} min(1, m p_(k) / k);
//   rejecting where it is at most q controls the false discovery rate at q for independent
//   tests. NaN p-values (genes with no counts) are left out of m and stay NaN.
//
// Cost: O(m p² + p³) per IRLS iteration and per evaluation of the adjusted likelihood; genes in
// parallel, each writing its own slot, so the results are the same on any number of threads.
//
// *Remark.* The negative binomial GLM, the dispersion's shrinkage, the tests and the FDR are
// general statistics; they move to the inference foundation when it exists (CLAUDE.md §8, open
// decision 3).

#pragma once

#include <genomics/expression/transcriptome.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/random/stream.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/solvers/roots.hpp>
#include <calx/special_functions/gamma.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace genomics::expression {

struct CountMatrix {
  std::size_t genes = 0, samples = 0;
  std::vector<double> counts;  // genes × samples, row-major

  CountMatrix() = default;
  CountMatrix(std::size_t g, std::size_t m) : genes(g), samples(m), counts(g * m, 0.0) {}
  double operator()(std::size_t g, std::size_t j) const { return counts[g * samples + j]; }
  double& at(std::size_t g, std::size_t j) { return counts[g * samples + j]; }
  const double* row(std::size_t g) const { return counts.data() + g * samples; }
};

// The median of a list (the mean of the two middle values for an even count); NaN when empty.
inline double median(std::vector<double> v) {
  if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
  const std::size_t h = v.size() / 2;
  std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(h), v.end());
  const double upper = v[h];
  if (v.size() % 2) return upper;
  return 0.5 * (upper + *std::max_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(h)));
}

// Anders and Huber's size factors (eq. 5).
inline std::vector<double> size_factors(const CountMatrix& k) {
  std::vector<double> logmean(k.genes, 0.0);
  std::vector<bool> keep(k.genes, true);
  for (std::size_t g = 0; g < k.genes; g++) {
    double s = 0;
    for (std::size_t j = 0; j < k.samples; j++) {
      const double x = k(g, j);
      if (!(x > 0)) { keep[g] = false; break; }
      s += std::log(x);
    }
    logmean[g] = s / static_cast<double>(k.samples);
  }
  std::vector<double> out(k.samples, 1.0), ratios;
  for (std::size_t j = 0; j < k.samples; j++) {
    ratios.clear();
    for (std::size_t g = 0; g < k.genes; g++)
      if (keep[g]) ratios.push_back(std::log(k(g, j)) - logmean[g]);
    if (ratios.empty()) throw std::invalid_argument("genomics::expression::size_factors: no gene is counted in every sample");
    out[j] = std::exp(median(ratios));
  }
  return out;
}

// ---- the design ------------------------------------------------------------------------------------
struct Design {
  std::size_t samples = 0, p = 0;
  std::vector<double> x;  // samples × p, row-major

  double operator()(std::size_t j, std::size_t r) const { return x[j * p + r]; }
  // An intercept and the indicator of group 1: β₁ is the log fold change of group 1 over group 0.
  static Design two_groups(const std::vector<int>& group) {
    Design d;
    d.samples = group.size();
    d.p = 2;
    for (int g : group) {
      if (g != 0 && g != 1) throw std::invalid_argument("genomics::expression::Design::two_groups: groups are 0 and 1");
      d.x.push_back(1.0);
      d.x.push_back(g == 1 ? 1.0 : 0.0);
    }
    return d;
  }
};

// ---- the negative binomial ------------------------------------------------------------------------------------
namespace detail {

// log Γ(y + r) − log Γ(r): a sum for small integer y, log Γ otherwise.
inline double log_gamma_ratio(double y, double r) {
  if (y == std::floor(y) && y < 40) {
    double s = 0;
    for (int i = 0; i < static_cast<int>(y); i++) s += std::log(r + i);
    return s;
  }
  return calx::special_functions::log_gamma(y + r) - calx::special_functions::log_gamma(r);
}

inline constexpr double MIN_MU = 0.5, RIDGE = 1e-6, MIN_DISPERSION = 1e-8;

// A = Xᵀ diag(w) X + λ I.
inline calx::linalg::Matrix gram(const Design& X, const std::vector<double>& w, double ridge) {
  calx::linalg::Matrix A(X.p);
  for (std::size_t j = 0; j < X.samples; j++)
    for (std::size_t a = 0; a < X.p; a++)
      for (std::size_t b = 0; b < X.p; b++) A(a, b) += w[j] * X(j, a) * X(j, b);
  for (std::size_t a = 0; a < X.p; a++) A(a, a) += ridge;
  return A;
}

// log det of a symmetric positive definite matrix by elimination (no pivoting is needed).
inline double log_det(calx::linalg::Matrix A) {
  const std::size_t n = A.n;
  double s = 0;
  for (std::size_t c = 0; c < n; c++) {
    if (!(A(c, c) > 0)) return -INFINITY;
    s += std::log(A(c, c));
    for (std::size_t r = c + 1; r < n; r++) {
      const double f = A(r, c) / A(c, c);
      for (std::size_t k = c; k < n; k++) A(r, k) -= f * A(c, k);
    }
  }
  return s;
}

}  // namespace detail

// Σ_j log P(y_j) under NB(μ_j, α).
inline double nb_log_likelihood(const double* y, const double* mu, std::size_t m, double alpha) {
  const double r = 1.0 / alpha;
  double ll = 0;
  for (std::size_t j = 0; j < m; j++) {
    const double l1p = std::log1p(alpha * mu[j]);  // log((r + μ)/r)
    ll += detail::log_gamma_ratio(y[j], r) - calx::special_functions::log_gamma(y[j] + 1) - r * l1p + y[j] * (std::log(mu[j]) - std::log(r + mu[j]));
  }
  return ll;
}

// A draw of NB(μ, α) as a gamma–Poisson mixture.
inline std::int64_t nb_draw(calx::random::Stream& r, double mu, double alpha) {
  if (!(mu > 0)) return 0;
  if (!(alpha > 0)) return r.poisson(mu);
  return r.poisson(r.gamma(1.0 / alpha, alpha * mu));
}

struct NbFit {
  std::vector<double> beta, se;  // natural log
  std::vector<double> mu;        // the fitted means
  double log_likelihood = 0;
  int iterations = 0;
  bool converged = false;
};

// The NB GLM of one gene's counts y (one per sample) at dispersion α, by IRLS.
inline NbFit fit_glm(const double* y, const std::vector<double>& s, const Design& X, double alpha, int max_iterations = 100) {
  const std::size_t m = X.samples, p = X.p;
  NbFit fit;
  fit.beta.assign(p, 0.0);
  fit.se.assign(p, 0.0);
  fit.mu.assign(m, 0.0);
  // Start: least squares of log(y/s + 0.1) on X.
  {
    std::vector<double> ones(m, 1.0);
    calx::linalg::Matrix A = detail::gram(X, ones, detail::RIDGE), B(p);
    for (std::size_t j = 0; j < m; j++)
      for (std::size_t a = 0; a < p; a++) B(a, 0) += X(j, a) * std::log(y[j] / s[j] + 0.1);
    const auto b = calx::linalg::solve(A, B);
    for (std::size_t a = 0; a < p; a++) fit.beta[a] = b(a, 0);
  }
  std::vector<double> w(m), z(m);
  double dev_old = INFINITY;
  calx::linalg::Matrix A(p);
  for (int it = 1; it <= max_iterations; it++) {
    for (std::size_t j = 0; j < m; j++) {
      double eta = 0;
      for (std::size_t a = 0; a < p; a++) eta += X(j, a) * fit.beta[a];
      const double mu = std::max(detail::MIN_MU, s[j] * std::exp(std::clamp(eta, -30.0, 30.0)));
      fit.mu[j] = mu;
      w[j] = mu / (1 + alpha * mu);
      z[j] = std::log(mu / s[j]) + (y[j] - mu) / mu;
    }
    A = detail::gram(X, w, detail::RIDGE);
    calx::linalg::Matrix B(p);
    for (std::size_t j = 0; j < m; j++)
      for (std::size_t a = 0; a < p; a++) B(a, 0) += X(j, a) * w[j] * z[j];
    const auto b = calx::linalg::solve(A, B);
    for (std::size_t a = 0; a < p; a++) fit.beta[a] = std::clamp(b(a, 0), -30.0, 30.0);
    for (std::size_t j = 0; j < m; j++) {
      double eta = 0;
      for (std::size_t a = 0; a < p; a++) eta += X(j, a) * fit.beta[a];
      fit.mu[j] = s[j] * std::exp(std::clamp(eta, -30.0, 30.0));
    }
    const double dev = -2 * nb_log_likelihood(y, fit.mu.data(), m, alpha);
    fit.iterations = it;
    if (std::fabs(dev - dev_old) / (std::fabs(dev) + 0.1) < 1e-8) {
      fit.converged = true;
      break;
    }
    dev_old = dev;
  }
  for (std::size_t j = 0; j < m; j++) w[j] = std::max(detail::MIN_MU, fit.mu[j]) / (1 + alpha * std::max(detail::MIN_MU, fit.mu[j]));
  A = detail::gram(X, w, detail::RIDGE);
  const auto inv = calx::linalg::solve(A, calx::linalg::Matrix::identity(p));
  for (std::size_t a = 0; a < p; a++) fit.se[a] = std::sqrt(std::max(0.0, inv(a, a)));
  fit.log_likelihood = nb_log_likelihood(y, fit.mu.data(), m, alpha);
  return fit;
}

// ---- the dispersion ------------------------------------------------------------------------------------
struct DispersionPrior {
  double mean_log = 0;   // log α_tr
  double variance = 0;   // σ_d²; 0 for none
};

// The objective in θ = log α at fixed means: ℓ(α) [− ½ log det(Xᵀ W X)] [+ the log prior].
inline double dispersion_objective(const double* y, const std::vector<double>& mu, const Design& X, double theta, bool cox_reid, const DispersionPrior& prior = {}) {
  const double alpha = std::exp(theta);
  double v = nb_log_likelihood(y, mu.data(), X.samples, alpha);
  if (cox_reid) {
    std::vector<double> w(X.samples);
    for (std::size_t j = 0; j < X.samples; j++) w[j] = mu[j] / (1 + alpha * mu[j]);
    v -= 0.5 * detail::log_det(detail::gram(X, w, 0.0));
  }
  if (prior.variance > 0) v -= (theta - prior.mean_log) * (theta - prior.mean_log) / (2 * prior.variance);
  return v;
}

// Its derivative in θ.
inline double dispersion_score(const double* y, const std::vector<double>& mu, const Design& X, double theta, bool cox_reid, const DispersionPrior& prior = {}) {
  const double alpha = std::exp(theta), r = 1.0 / alpha;
  double dr = 0;  // dℓ/dr
  for (std::size_t j = 0; j < X.samples; j++)
    dr += calx::special_functions::digamma(y[j] + r) - calx::special_functions::digamma(r) - std::log1p(mu[j] / r) + (mu[j] - y[j]) / (r + mu[j]);
  double d = -r * dr;
  if (cox_reid) {
    const std::size_t p = X.p;
    std::vector<double> w(X.samples), v(X.samples);
    for (std::size_t j = 0; j < X.samples; j++) {
      w[j] = mu[j] / (1 + alpha * mu[j]);
      v[j] = w[j] * w[j];  // μ²/(1 + αμ)²
    }
    const auto inv = calx::linalg::solve(detail::gram(X, w, 0.0), calx::linalg::Matrix::identity(p));
    const auto Bm = detail::gram(X, v, 0.0);
    double tr = 0;
    for (std::size_t a = 0; a < p; a++)
      for (std::size_t b = 0; b < p; b++) tr += inv(a, b) * Bm(b, a);
    d += 0.5 * alpha * tr;
  }
  if (prior.variance > 0) d -= (theta - prior.mean_log) / prior.variance;
  return d;
}

// The α maximising the objective over [10⁻⁸, max(10, m)]: a grid of 25, then Brent's method on
// the derivative inside the best point's bracket.
inline double estimate_dispersion(const double* y, const std::vector<double>& mu, const Design& X, bool cox_reid, const DispersionPrior& prior = {}) {
  const double lo = std::log(detail::MIN_DISPERSION), hi = std::log(std::max(10.0, static_cast<double>(X.samples)));
  constexpr int G = 25;
  double best = -INFINITY;
  int at = 0;
  for (int k = 0; k < G; k++) {
    const double th = lo + (hi - lo) * k / (G - 1);
    const double v = dispersion_objective(y, mu, X, th, cox_reid, prior);
    if (v > best) { best = v; at = k; }
  }
  const double step = (hi - lo) / (G - 1);
  const double a = at > 0 ? lo + step * (at - 1) : lo, b = at < G - 1 ? lo + step * (at + 1) : hi;
  const auto f = [&](double th) { return dispersion_score(y, mu, X, th, cox_reid, prior); };
  const double fa = f(a), fb = f(b);
  double theta = lo + step * at;
  if (fa > 0 && fb < 0) theta = calx::solvers::brent(f, a, b, 1e-10).x;
  else if (at == 0 && fa <= 0) theta = lo;
  else if (at == G - 1 && fb >= 0) theta = hi;
  return std::exp(theta);
}

// The method-of-moments start: (var − mean · mean(1/s)) / mean² of the normalised counts.
inline double moments_dispersion(const double* y, const std::vector<double>& s) {
  const std::size_t m = s.size();
  double mean = 0, inv = 0;
  for (std::size_t j = 0; j < m; j++) {
    mean += y[j] / s[j];
    inv += 1 / s[j];
  }
  mean /= static_cast<double>(m);
  inv /= static_cast<double>(m);
  double var = 0;
  for (std::size_t j = 0; j < m; j++) var += (y[j] / s[j] - mean) * (y[j] / s[j] - mean);
  var /= static_cast<double>(m > 1 ? m - 1 : 1);
  return mean > 0 ? std::clamp((var - mean * inv) / (mean * mean), detail::MIN_DISPERSION, 10.0) : detail::MIN_DISPERSION;
}

// ---- Benjamini–Hochberg ------------------------------------------------------------------------------------
inline std::vector<double> benjamini_hochberg(const std::vector<double>& p) {
  std::vector<std::size_t> order;
  for (std::size_t i = 0; i < p.size(); i++)
    if (!std::isnan(p[i])) order.push_back(i);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return p[a] < p[b]; });
  std::vector<double> q(p.size(), std::numeric_limits<double>::quiet_NaN());
  const double m = static_cast<double>(order.size());
  double running = 1.0;
  for (std::size_t k = order.size(); k-- > 0;) {
    running = std::min(running, std::min(1.0, m * p[order[k]] / static_cast<double>(k + 1)));
    q[order[k]] = running;
  }
  return q;
}

// ---- the whole analysis ------------------------------------------------------------------------------------
struct DifferentialOptions {
  bool cox_reid = true;
  bool shrink = true;            // the MAP toward the trend; false keeps the gene-wise estimates
  std::size_t coefficient = 1;   // the coefficient tested
  double outlier_sd = 2;         // genes this many σ_d above the trend keep their own estimate
  std::vector<double> size_factors;  // empty: estimated
};

struct GeneResult {
  double base_mean = 0;                // the mean of the normalised counts
  double log2_fold_change = 0, lfc_se = 0, statistic = 0;
  double p_value = std::numeric_limits<double>::quiet_NaN(), padj = std::numeric_limits<double>::quiet_NaN();
  double dispersion_gene = 0, dispersion_trend = 0, dispersion = 0;
  bool dispersion_outlier = false, converged = false;
};

struct Differential {
  std::vector<double> size_factors;
  std::vector<GeneResult> genes;
  double trend_a0 = 0, trend_a1 = 0;   // α_tr(μ̄) = a₀ + a₁ / μ̄
  bool trend_parametric = false;       // false: the trend is the mean estimate
  double log_dispersion_variance = 0;  // s²_lr, the squared MAD of the log residuals
  double prior_variance = 0;           // σ_d² = max(s²_lr − ψ′((m − p)/2), 0.25)
};

namespace detail {

// The gamma GLM with identity link of α on (1, 1/μ̄), by IRLS (weights 1/fitted²).
inline bool gamma_identity_fit(const std::vector<double>& alpha, const std::vector<double>& mean, const std::vector<std::size_t>& use, double& a0, double& a1) {
  for (int it = 0; it < 50; it++) {
    double s00 = 0, s01 = 0, s11 = 0, t0 = 0, t1 = 0;
    for (std::size_t g : use) {
      const double x = 1.0 / mean[g], fitted = std::max(a0 + a1 * x, 1e-12), w = 1.0 / (fitted * fitted);
      s00 += w;
      s01 += w * x;
      s11 += w * x * x;
      t0 += w * alpha[g];
      t1 += w * x * alpha[g];
    }
    const double det = s00 * s11 - s01 * s01;
    if (!(std::fabs(det) > 0)) return false;
    const double n0 = (s11 * t0 - s01 * t1) / det, n1 = (s00 * t1 - s01 * t0) / det;
    const bool done = std::fabs(n0 - a0) <= 1e-10 * std::fabs(n0) && std::fabs(n1 - a1) <= 1e-10 * std::fabs(n1);
    a0 = n0;
    a1 = n1;
    if (done) break;
  }
  return std::isfinite(a0) && std::isfinite(a1);
}

}  // namespace detail

inline Differential differential(const CountMatrix& k, const Design& X, const DifferentialOptions& o, calx::scheduler::ThreadPool& pool) {
  if (X.samples != k.samples || X.p == 0 || X.p >= k.samples) throw std::invalid_argument("genomics::expression::differential: a design of p < m columns for the samples");
  if (o.coefficient >= X.p) throw std::invalid_argument("genomics::expression::differential: no such coefficient");
  Differential out;
  out.size_factors = o.size_factors.empty() ? size_factors(k) : o.size_factors;
  if (out.size_factors.size() != k.samples) throw std::invalid_argument("genomics::expression::differential: one size factor per sample");
  const auto& s = out.size_factors;
  const std::size_t G = k.genes, m = k.samples;
  out.genes.assign(G, GeneResult{});
  std::vector<std::vector<double>> mu(G);
  std::vector<bool> counted(G, false);
  // Gene-wise dispersions: two passes of fit and estimate.
  calx::scheduler::parallel_for(pool, G, 64, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t g = b; g < e; g++) {
      const double* y = k.row(g);
      GeneResult& r = out.genes[g];
      double total = 0;
      for (std::size_t j = 0; j < m; j++) {
        total += y[j];
        r.base_mean += y[j] / s[j] / static_cast<double>(m);
      }
      if (!(total > 0)) continue;
      counted[g] = true;
      double alpha = moments_dispersion(y, s);
      for (int pass = 0; pass < 2; pass++) {
        const NbFit f = fit_glm(y, s, X, alpha);
        mu[g] = f.mu;
        for (double& x : mu[g]) x = std::max(x, detail::MIN_MU);
        alpha = estimate_dispersion(y, mu[g], X, o.cox_reid);
      }
      r.dispersion_gene = alpha;
    }
  });
  // The trend.
  std::vector<double> gene_alpha(G), mean(G);
  std::vector<std::size_t> use;
  for (std::size_t g = 0; g < G; g++) {
    gene_alpha[g] = out.genes[g].dispersion_gene;
    mean[g] = out.genes[g].base_mean;
    if (counted[g] && gene_alpha[g] >= 100 * detail::MIN_DISPERSION && mean[g] > 0) use.push_back(g);
  }
  double a0 = 0.1, a1 = 1.0;  // DESeq2's start
  bool parametric = !use.empty();
  for (int iter = 0; iter <= 10 && parametric; iter++) {
    std::vector<std::size_t> kept;
    for (std::size_t g : use) {
      const double ratio = gene_alpha[g] / (a0 + a1 / mean[g]);
      if (ratio > 1e-4 && ratio < 15) kept.push_back(g);
    }
    double n0 = a0, n1 = a1;
    if (kept.size() < 3 || !detail::gamma_identity_fit(gene_alpha, mean, kept, n0, n1) || !(n0 > 0) || !(n1 > 0)) { parametric = false; break; }
    const double change = std::pow(std::log(n0 / a0), 2) + std::pow(std::log(n1 / a1), 2);
    a0 = n0;
    a1 = n1;
    if (change < 1e-6) break;
  }
  if (!parametric) {
    double sum = 0;
    for (std::size_t g : use) sum += gene_alpha[g];
    a0 = use.empty() ? 0.1 : sum / static_cast<double>(use.size());
    a1 = 0;
  }
  out.trend_a0 = a0;
  out.trend_a1 = a1;
  out.trend_parametric = parametric;
  // The prior's variance: the MAD of the log residuals, less the expected sampling variance.
  std::vector<double> residual;
  for (std::size_t g : use) residual.push_back(std::log(gene_alpha[g]) - std::log(a0 + a1 / mean[g]));
  if (!residual.empty()) {
    const double med = median(residual);
    std::vector<double> dev;
    for (double x : residual) dev.push_back(std::fabs(x - med));
    const double mad = 1.4826 * median(dev);
    out.log_dispersion_variance = mad * mad;
    out.prior_variance = std::max(mad * mad - calx::special_functions::trigamma(static_cast<double>(m - X.p) / 2), 0.25);
  } else {
    out.prior_variance = 0.25;
  }
  // MAP dispersions, the final fits and the Wald tests.
  calx::scheduler::parallel_for(pool, G, 64, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t g = b; g < e; g++) {
      if (!counted[g]) continue;
      GeneResult& r = out.genes[g];
      const double* y = k.row(g);
      r.dispersion_trend = a0 + a1 / std::max(r.base_mean, 1e-12);
      const double log_trend = std::log(r.dispersion_trend);
      double alpha = r.dispersion_gene;
      if (o.shrink) {
        r.dispersion_outlier = std::log(r.dispersion_gene) > log_trend + o.outlier_sd * std::sqrt(out.log_dispersion_variance);
        if (!r.dispersion_outlier) alpha = estimate_dispersion(y, mu[g], X, o.cox_reid, {log_trend, out.prior_variance});
      }
      r.dispersion = alpha;
      const NbFit f = fit_glm(y, s, X, alpha);
      r.converged = f.converged;
      const double beta = f.beta[o.coefficient], se = f.se[o.coefficient];
      r.log2_fold_change = beta / std::log(2.0);
      r.lfc_se = se / std::log(2.0);
      r.statistic = se > 0 ? beta / se : 0;
      r.p_value = calx::special_functions::gamma_q(0.5, r.statistic * r.statistic / 2);
    }
  });
  std::vector<double> p(G);
  for (std::size_t g = 0; g < G; g++) p[g] = out.genes[g].p_value;
  const auto q = benjamini_hochberg(p);
  for (std::size_t g = 0; g < G; g++) out.genes[g].padj = q[g];
  return out;
}

// Counts of genes at base means q_g (per unit size factor), dispersions α_g and log fold
// changes β_g (natural log) of group 1 over group 0, in samples of the given groups and size
// factors; gene g from its own stream.
inline CountMatrix simulate_counts(const std::vector<double>& base, const std::vector<double>& dispersion, const std::vector<double>& log_fold,
                                   const std::vector<int>& group, const std::vector<double>& size, std::uint64_t seed) {
  if (dispersion.size() != base.size() || log_fold.size() != base.size() || size.size() != group.size())
    throw std::invalid_argument("genomics::expression::simulate_counts: inconsistent sizes");
  CountMatrix k(base.size(), group.size());
  for (std::size_t g = 0; g < base.size(); g++) {
    calx::random::Stream r(seed, COUNTS, static_cast<std::uint32_t>(g), 0);
    for (std::size_t j = 0; j < group.size(); j++) {
      const double mu = size[j] * base[g] * (group[j] == 1 ? std::exp(log_fold[g]) : 1.0);
      k.at(g, j) = static_cast<double>(nb_draw(r, mu, dispersion[g]));
    }
  }
  return k;
}

}  // namespace genomics::expression
