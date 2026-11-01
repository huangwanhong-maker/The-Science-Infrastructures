// Association scans: every variant's effect on a trait, adjusted for covariates, as genome-wide
// association studies estimate them; genomic control; the quantile–quantile plot's points and band.
//
// THE LINEAR SCAN. For each site j, the regression y = β g_j + Cγ + e, with C the covariates (an
// intercept, sex, age, principal components; Price et al. 2006, Nature Genetics 38, 904–909, for the
// components), fitted by least squares. By the Frisch–Waugh–Lovell theorem (Frisch and Waugh 1933;
// Lovell 1963, JASA 58, 993–1010) β̂ and its standard error are those of the simple regression of
// ỹ = M y on g̃ = M g_j, M = I − Q Qᵀ the projection off C's columns (Q orthonormal, from calx's
// pivoted QR of C with unit-norm columns): β̂ = g̃ᵀỹ/g̃ᵀg̃, σ̂² = (ỹᵀỹ − β̂² g̃ᵀg̃)/(n − r − 1),
// se = σ̂/√(g̃ᵀg̃), with r = rank C; t = β̂/se on n − r − 1 degrees of freedom (exact under Gaussian
// errors). So ỹ is formed once and each site costs O(n r): the vectorised scan of PLINK 2's --glm and of
// REGENIE's second step, equal to tyche's ols on [C, g_j] (checks).
//
// THE LOGISTIC SCAN. A case–control trait is fitted by the logistic model logit P(y = 1) = β g + Cγ.
// The null model (β = 0) is fitted once by tyche's glm; each site then gets the score test
// (Rao 1948): U = g̃ᵀ(y − μ̂₀), its efficient information V = gᵀWg − gᵀWC(CᵀWC)⁻¹CᵀWg with W = diag
// μ̂₀(1 − μ̂₀) (formed as ‖√W g‖² − ‖Q_Wᵀ√W g‖², Q_W orthonormal for √W C), χ² = U²/V on one degree of
// freedom; the one-step estimate β̂ = U/V with standard error 1/√V. This is the score test SAIGE and
// REGENIE use for scans (Zhou et al. 2018, Nature Genetics 50, 1335–1341). The sites whose score p
// falls below wald_below are refitted by tyche's glm and their Wald β̂, standard error and p kept apart.
//
// GENOMIC CONTROL (Devlin and Roeder 1999, Biometrics 55, 997–1004): λ = median χ²/F⁻¹_{χ²₁}(½), the
// χ² read from each p-value (F⁻¹_{χ²₁}(½) = 0.4549…), near 1 under the null. Under polygenic
// inheritance λ exceeds 1 without any confounding (Yang et al. 2011, EJHG 19, 807–812), which LD score
// regression separates (heritability.hpp).
//
// THE Q–Q PLOT. The i-th smallest of m independent uniform p-values is Beta(i, m − i + 1): expected
// −log₁₀ at i/(m + 1), and a pointwise 95% band from that beta's quantiles. Points are thinned for
// drawing: every one with p below thin_below, an even sample by rank of the rest.
//
// MISSING CALLS are replaced by the site's mean dosage over the samples (the scan's own imputation),
// which keeps g̃ orthogonal to C as required. A site monomorphic in the samples gets NaN.
//
// COMPLEXITY. O(n m r) for m sites, in parallel over sites (each writes its own slot: the same on any
// number of threads); the logistic Wald refits O(n r²) each for the few that need them.

#pragma once

#include <evolution/population/genotypes.hpp>

#include <tyche/probability/distributions.hpp>
#include <tyche/regression/glm.hpp>
#include <tyche/regression/linear.hpp>

#include <calx/linalg/qr.hpp>
#include <calx/linalg/rectangular.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace evolution::association {

using calx::linalg::Rectangular;

// The projection off a covariate matrix's columns: M v = v − Q Qᵀ v.
class Residualizer {
 public:
  Residualizer() = default;
  explicit Residualizer(const Rectangular& C) : n_(C.rows) {
    if (C.cols == 0) return;
    Rectangular A = C;
    for (std::size_t j = 0; j < A.cols; j++) {
      double ss = 0;
      for (std::size_t i = 0; i < n_; i++) ss += A(i, j) * A(i, j);
      const double norm = std::sqrt(ss);
      if (norm > 0)
        for (std::size_t i = 0; i < n_; i++) A(i, j) /= norm;
    }
    const auto f = calx::linalg::qr(std::move(A), true);
    q_ = calx::linalg::thin_q(f, f.rank);
    rank_ = f.rank;
  }
  std::size_t rank() const { return rank_; }
  std::size_t rows() const { return n_; }
  // v ← M v.
  void apply(std::span<double> v) const {
    for (std::size_t a = 0; a < rank_; a++) {
      double s = 0;
      for (std::size_t i = 0; i < n_; i++) s += q_(i, a) * v[i];
      for (std::size_t i = 0; i < n_; i++) v[i] -= s * q_(i, a);
    }
  }
  std::vector<double> residuals(std::span<const double> v) const {
    std::vector<double> r(v.begin(), v.end());
    apply(r);
    return r;
  }

 private:
  std::size_t n_ = 0, rank_ = 0;
  Rectangular q_;
};

// The dosages of one site over the samples listed, missing calls set to the samples' mean.
inline std::vector<double> sample_dosages(const population::Genotypes& g, std::size_t site, std::span<const std::size_t> samples, double* mean = nullptr) {
  const auto row = g.site(site);
  std::vector<double> d(samples.size());
  double s = 0, typed = 0;
  for (std::size_t k = 0; k < samples.size(); k++) {
    const std::int8_t x = row[samples[k]];
    d[k] = x;
    if (x >= 0) {
      s += x;
      typed++;
    }
  }
  const double m = typed > 0 ? s / typed : 0;
  for (std::size_t k = 0; k < samples.size(); k++)
    if (row[samples[k]] < 0) d[k] = m;
  if (mean) *mean = m;
  return d;
}

struct Scan {
  std::size_t samples = 0, sites = 0, covariates = 0;  // covariates: the rank of C
  double df = 0;                                       // n − r − 1 (linear)
  bool logistic = false;
  std::vector<double> beta, se, statistic, p;          // per site: β̂, its se, t (linear) or the signed √χ² (score), p
  std::vector<double> frequency;                       // the alternate allele's frequency in the samples
  std::vector<double> wald_beta, wald_se, wald_p;      // logistic: the glm's Wald refits (NaN where not refitted)
};

// The linear scan of y (one value per listed sample) on every site, adjusted for C (samples × k).
inline Scan linear_scan(const population::Genotypes& g, std::span<const std::size_t> samples, std::span<const double> y, const Rectangular& C,
                        calx::scheduler::ThreadPool* pool = nullptr) {
  const std::size_t n = samples.size(), m = g.sites;
  if (y.size() != n || C.rows != n) throw std::invalid_argument("evolution::association::linear_scan: one y and one row of C per sample");
  const Residualizer M(C);
  const std::size_t r = M.rank();
  if (n < r + 3) throw std::invalid_argument("evolution::association::linear_scan: more samples than covariates");
  const auto yt = M.residuals(y);
  double yy = 0;
  for (const double v : yt) yy += v * v;
  Scan s;
  s.samples = n;
  s.sites = m;
  s.covariates = r;
  s.df = static_cast<double>(n - r - 1);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  s.beta.assign(m, nan);
  s.se = s.statistic = s.p = s.beta;
  s.frequency.assign(m, nan);
  const tyche::probability::StudentT t{s.df};
  const auto body = [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t j = begin; j < end; j++) {
      double mean = 0;
      auto d = sample_dosages(g, j, samples, &mean);
      s.frequency[j] = mean / 2;
      M.apply(d);
      double gg = 0, gy = 0;
      for (std::size_t i = 0; i < n; i++) {
        gg += d[i] * d[i];
        gy += d[i] * yt[i];
      }
      if (!(gg > 1e-12 * static_cast<double>(n))) continue;
      const double b = gy / gg, rss = std::max(0.0, yy - b * gy), se = std::sqrt(rss / s.df / gg);
      s.beta[j] = b;
      s.se[j] = se;
      s.statistic[j] = b / se;
      s.p[j] = std::min(1.0, 2 * t.survival(std::fabs(b / se)));
    }
  };
  if (pool) calx::scheduler::parallel_for(*pool, m, 64, body);
  else body(0, m, 0);
  return s;
}

// The logistic scan of a 0/1 trait: the score test at every site, the Wald test where its p < wald_below.
inline Scan logistic_scan(const population::Genotypes& g, std::span<const std::size_t> samples, std::span<const double> cases, const Rectangular& C,
                          calx::scheduler::ThreadPool* pool = nullptr, double wald_below = 1e-4) {
  const std::size_t n = samples.size(), m = g.sites, k = C.cols;
  if (cases.size() != n || C.rows != n) throw std::invalid_argument("evolution::association::logistic_scan: one case label and one row of C per sample");
  tyche::regression::GlmOptions o;
  o.family = tyche::regression::Family::binomial;
  o.drop_collinear = true;
  const auto null = tyche::regression::glm(C, cases, o);
  std::vector<double> resid(n), sw(n);
  for (std::size_t i = 0; i < n; i++) {
    resid[i] = cases[i] - null.fitted[i];
    sw[i] = std::sqrt(null.fitted[i] * (1 - null.fitted[i]));
  }
  Rectangular WC(n, k);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t a = 0; a < k; a++) WC(i, a) = sw[i] * C(i, a);
  const Residualizer MW(WC);
  Scan s;
  s.samples = n;
  s.sites = m;
  s.covariates = MW.rank();
  s.df = static_cast<double>(n - MW.rank() - 1);
  s.logistic = true;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  s.beta.assign(m, nan);
  s.se = s.statistic = s.p = s.wald_beta = s.wald_se = s.wald_p = s.beta;
  s.frequency.assign(m, nan);
  const tyche::probability::ChiSquared chi{1};
  const auto body = [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t j = begin; j < end; j++) {
      double mean = 0;
      auto d = sample_dosages(g, j, samples, &mean);
      s.frequency[j] = mean / 2;
      double U = 0;
      std::vector<double> wd(n);
      for (std::size_t i = 0; i < n; i++) {
        U += d[i] * resid[i];
        wd[i] = sw[i] * d[i];
      }
      MW.apply(wd);
      double V = 0;
      for (const double v : wd) V += v * v;
      if (!(V > 1e-12 * static_cast<double>(n))) continue;
      const double x2 = U * U / V;
      s.beta[j] = U / V;
      s.se[j] = 1 / std::sqrt(V);
      s.statistic[j] = (U >= 0 ? 1 : -1) * std::sqrt(x2);
      s.p[j] = chi.survival(x2);
    }
  };
  if (pool) calx::scheduler::parallel_for(*pool, m, 64, body);
  else body(0, m, 0);
  // The Wald refits of the strongest sites.
  std::vector<std::size_t> refit;
  for (std::size_t j = 0; j < m; j++)
    if (s.p[j] < wald_below) refit.push_back(j);
  const auto wald = [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t q = begin; q < end; q++) {
      const std::size_t j = refit[q];
      const auto d = sample_dosages(g, j, samples);
      Rectangular X(n, k + 1);
      for (std::size_t i = 0; i < n; i++) {
        for (std::size_t a = 0; a < k; a++) X(i, a) = C(i, a);
        X(i, k) = d[i];
      }
      try {
        const auto fit = tyche::regression::glm(X, cases, o);
        s.wald_beta[j] = fit.beta[k];
        s.wald_se[j] = fit.se[k];
        s.wald_p[j] = fit.p[k];
      } catch (const std::exception&) {
      }
    }
  };
  if (pool) calx::scheduler::parallel_for(*pool, refit.size(), 4, wald);
  else wald(0, refit.size(), 0);
  return s;
}

// λ = median χ²₁ / F⁻¹(½), the χ² of each finite p-value.
inline double genomic_control(std::span<const double> p) {
  std::vector<double> x;
  const tyche::probability::ChiSquared chi{1};
  for (const double v : p)
    if (v >= 0 && v <= 1) x.push_back(chi.upper_quantile(v));
  if (x.empty()) return std::numeric_limits<double>::quiet_NaN();
  const std::size_t h = x.size() / 2;
  std::nth_element(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(h), x.end());
  double med = x[h];
  if (x.size() % 2 == 0) {
    const double lower = *std::max_element(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(h));
    med = (med + lower) / 2;
  }
  return med / chi.quantile(0.5);
}

struct QqPoint {
  double expected = 0, observed = 0;  // −log₁₀ p
  double lower = 0, upper = 0;        // the pointwise 95% band, −log₁₀
};

// The Q–Q plot's points (thinned: all with p < thin_below, about `keep` of the others evenly by rank).
inline std::vector<QqPoint> qq_points(std::span<const double> p, double thin_below = 1e-2, std::size_t keep = 400) {
  std::vector<double> x;
  for (const double v : p)
    if (v >= 0 && v <= 1) x.push_back(std::max(v, 1e-300));
  std::sort(x.begin(), x.end());
  const std::size_t m = x.size();
  std::size_t rest = 0;
  for (const double v : x) rest += v >= thin_below;
  const std::size_t step = std::max<std::size_t>(1, rest / std::max<std::size_t>(1, keep));
  std::vector<QqPoint> out;
  std::size_t seen = 0;
  for (std::size_t i = 0; i < m; i++) {
    if (x[i] >= thin_below && (seen++ % step) != 0 && i + 1 != m) continue;
    const double a = static_cast<double>(i + 1), b = static_cast<double>(m - i);
    const tyche::probability::Beta law{a, b};
    QqPoint q;
    q.expected = -std::log10(a / static_cast<double>(m + 1));
    q.observed = -std::log10(x[i]);
    q.lower = -std::log10(law.upper_quantile(0.025));
    q.upper = -std::log10(law.quantile(0.025));
    out.push_back(q);
  }
  return out;
}

}  // namespace evolution::association
