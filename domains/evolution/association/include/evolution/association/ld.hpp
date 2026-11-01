// Linkage disequilibrium as association methods read it: the correlation matrix of a region's
// dosages, LD scores, clumping of association signals, and principal components of the genotypes.
//
// - The correlation r_jk of two sites' dosages over a sample (Pearson's, missing calls at the site's
//   mean) is what summary-statistic methods take as LD: SuSiE-RSS, S-PrediXcan, LD score regression,
//   LDpred. For phased haplotypes it is close to the haplotypic r of population/linkage.hpp (equal in
//   expectation under Hardy–Weinberg equilibrium); signed, so that z-scores of correlated sites align.
// - The LD score of site j (Bulik-Sullivan et al. 2015, Nature Genetics 47, 291–295, eq. 1 and its
//   Methods) is ℓ_j = Σ_k r²_jk over the sites within `window` base pairs (itself included), each r²
//   corrected for its sampling bias, r̃² = r² − (1 − r²)/(n − 2) (their Methods, after Hill and
//   Robertson's E r̂² ≈ r² + (1 − r²)/n).
// - Clumping (PLINK's --clump; Purcell et al. 2007, AJHG 81, 559–575): sites with p below p1 in
//   increasing order of p become index sites; each takes into its clump the sites within `window`
//   base pairs with r² above r2 and p below p2 not yet clumped. The defaults are PLINK's (p1 = 10⁻⁴,
//   p2 = 10⁻², r² = 0.5, 250 kb); polygenic scores clump with p1 = p2 = 1 and r² = 0.1.
// - Principal components of the genotypes for adjusting scans (Price et al. 2006): the leading
//   eigenvectors of the samples' kinship M Mᵀ/m (population/structure.hpp's pca, Patterson et al. 2006)
//   are obtained here through the m′ × m′ covariance MᵀM/n of a subset of m′ sites (every site, or an
//   evenly spaced subset, LD pruning by thinning), by the duality of the singular value decomposition:
//   if MᵀM v = λ v then M v/√λ is a unit eigenvector of M Mᵀ with the same eigenvalue. For n samples and
//   m′ sites this costs O(n m′² + m′³) instead of O(n² m′ + n³), the difference between a cohort of
//   thousands and one of hundreds; the axes are population/structure.hpp's (checks), on Patterson's
//   scale, with its sign convention (the largest entry positive).
//
// COMPLEXITY. A region's correlation matrix O(n p²); LD scores O(n m w) for w sites a window, in
// parallel over sites; clumping O(I w n) for I index sites.

#pragma once

#include <evolution/association/scan.hpp>
#include <evolution/population/genotypes.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigen.hpp>
#include <calx/linalg/rectangular.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

namespace evolution::association {

// The sites' dosages over the samples, centred and scaled to unit variance (divisor n − 1), as an
// n × p matrix (a constant site is all zero).
inline Rectangular standardised_dosages(const population::Genotypes& g, std::span<const std::size_t> samples, std::span<const std::size_t> sites) {
  const std::size_t n = samples.size(), p = sites.size();
  Rectangular Z(n, p);
  for (std::size_t c = 0; c < p; c++) {
    double mean = 0;
    const auto d = sample_dosages(g, sites[c], samples, &mean);
    double ss = 0;
    for (const double v : d) ss += (v - mean) * (v - mean);
    const double sd = std::sqrt(ss / static_cast<double>(n - 1));
    for (std::size_t i = 0; i < n; i++) Z(i, c) = sd > 0 ? (d[i] - mean) / sd : 0;
  }
  return Z;
}

// The correlation matrix of the sites' dosages over the samples (p × p).
inline calx::linalg::Matrix ld_correlation(const population::Genotypes& g, std::span<const std::size_t> samples, std::span<const std::size_t> sites,
                                           calx::scheduler::ThreadPool* pool = nullptr) {
  const Rectangular Z = standardised_dosages(g, samples, sites);
  const std::size_t n = Z.rows, p = Z.cols;
  // Column-major copy for the inner products.
  std::vector<double> col(n * p);
  for (std::size_t c = 0; c < p; c++)
    for (std::size_t i = 0; i < n; i++) col[c * n + i] = Z(i, c);
  calx::linalg::Matrix R(p);
  const auto rows = [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t a = begin; a < end; a++)
      for (std::size_t b = 0; b <= a; b++) {
        double s = 0;
        for (std::size_t i = 0; i < n; i++) s += col[a * n + i] * col[b * n + i];
        R(a, b) = s / static_cast<double>(n - 1);
      }
  };
  if (pool) calx::scheduler::parallel_for(*pool, p, 8, rows);
  else rows(0, p, 0);
  for (std::size_t a = 0; a < p; a++) {
    if (R(a, a) == 0) R(a, a) = 1;  // a constant site: uncorrelated with the rest, unit on the diagonal
    for (std::size_t b = 0; b < a; b++) R(b, a) = R(a, b);
  }
  return R;
}

// The sites of a chromosome within [from, to) base pairs, in order.
inline std::vector<std::size_t> sites_between(const population::Genotypes& g, std::uint32_t chromosome, std::uint64_t from, std::uint64_t to) {
  const auto [lo, hi] = g.chromosome_sites(chromosome);
  std::vector<std::size_t> out;
  for (std::size_t j = lo; j < hi; j++)
    if (g.positions[j] >= from && g.positions[j] < to) out.push_back(j);
  return out;
}

// LD scores ℓ_j = Σ_k r̃²_jk over the sites within `window` base pairs (bias-corrected r²).
inline std::vector<double> ld_scores(const population::Genotypes& g, std::span<const std::size_t> samples, std::uint64_t window,
                                     calx::scheduler::ThreadPool* pool = nullptr, bool corrected = true) {
  const std::size_t n = samples.size(), m = g.sites;
  // Standardised dosages of every site, site-major.
  std::vector<double> z(m * n);
  for (std::size_t j = 0; j < m; j++) {
    double mean = 0;
    const auto d = sample_dosages(g, j, samples, &mean);
    double ss = 0;
    for (const double v : d) ss += (v - mean) * (v - mean);
    const double sd = std::sqrt(ss / static_cast<double>(n - 1));
    for (std::size_t i = 0; i < n; i++) z[j * n + i] = sd > 0 ? (d[i] - mean) / sd : 0;
  }
  std::vector<double> l(m, 0.0);
  const double dn = static_cast<double>(n);
  const auto body = [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t j = begin; j < end; j++) {
      const auto [lo, hi] = g.chromosome_sites(g.chromosome[j]);
      double s = 0;
      for (std::size_t k = lo; k < hi; k++) {
        const auto d = g.positions[k] > g.positions[j] ? g.positions[k] - g.positions[j] : g.positions[j] - g.positions[k];
        if (d > window) continue;
        double c = 0;
        for (std::size_t i = 0; i < n; i++) c += z[j * n + i] * z[k * n + i];
        const double r = c / (dn - 1), r2 = r * r;
        s += corrected && k != j ? r2 - (1 - r2) / (dn - 2) : r2;
      }
      l[j] = s;
    }
  };
  if (pool) calx::scheduler::parallel_for(*pool, m, 16, body);
  else body(0, m, 0);
  return l;
}

struct ClumpOptions {
  double p1 = 1e-4, p2 = 1e-2, r2 = 0.5;
  std::uint64_t window = 250'000;
};

struct Clump {
  std::size_t index = 0;               // the index site
  std::vector<std::size_t> members;    // the sites it took in (not itself), increasing
};

// Clumps of the scan's p-values (one per site of g), r² from the samples' dosages.
inline std::vector<Clump> clump(const population::Genotypes& g, std::span<const std::size_t> samples, std::span<const double> p, const ClumpOptions& o = {}) {
  const std::size_t m = g.sites, n = samples.size();
  if (p.size() != m) throw std::invalid_argument("evolution::association::clump: one p-value per site");
  std::vector<std::size_t> order;
  for (std::size_t j = 0; j < m; j++)
    if (p[j] < o.p1) order.push_back(j);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return p[a] < p[b]; });
  std::vector<char> taken(m, 0);
  std::vector<Clump> out;
  const auto standardised = [&](std::size_t j) {
    double mean = 0;
    auto d = sample_dosages(g, j, samples, &mean);
    double ss = 0;
    for (const double v : d) ss += (v - mean) * (v - mean);
    const double sd = std::sqrt(ss / static_cast<double>(n - 1));
    for (double& v : d) v = sd > 0 ? (v - mean) / sd : 0;
    return d;
  };
  for (const std::size_t j : order) {
    if (taken[j]) continue;
    taken[j] = 1;
    Clump c;
    c.index = j;
    const auto zj = standardised(j);
    const auto [lo, hi] = g.chromosome_sites(g.chromosome[j]);
    for (std::size_t k = lo; k < hi; k++) {
      if (taken[k] || !(p[k] < o.p2)) continue;
      const auto d = g.positions[k] > g.positions[j] ? g.positions[k] - g.positions[j] : g.positions[j] - g.positions[k];
      if (d > o.window) continue;
      const auto zk = standardised(k);
      double s = 0;
      for (std::size_t i = 0; i < n; i++) s += zj[i] * zk[i];
      const double r = s / static_cast<double>(n - 1);
      if (r * r > o.r2) {
        taken[k] = 1;
        c.members.push_back(k);
      }
    }
    out.push_back(std::move(c));
  }
  return out;
}

struct Components {
  std::size_t samples = 0, components = 0, sites_used = 0;
  std::vector<double> values;   // the leading eigenvalues on Patterson's scale (all n of them sum to n − 1)
  std::vector<double> vectors;  // samples × components, unit columns
};

// The leading principal components of the samples' genotypes at the sites listed (empty: every
// site), standardised as population::kinship standardises them, through the m′ × m′ covariance.
inline Components principal_components(const population::Genotypes& g, std::span<const std::size_t> samples, std::span<const std::size_t> sites,
                                       std::size_t components) {
  std::vector<std::size_t> use(sites.begin(), sites.end());
  if (use.empty()) {
    use.resize(g.sites);
    std::iota(use.begin(), use.end(), 0);
  }
  const std::size_t n = samples.size();
  // Patterson's standardisation: centred by the mean dosage, scaled by √(p(1 − p)), missing → 0.
  std::vector<std::size_t> kept;
  std::vector<double> M;  // site-major, kept × n
  for (const std::size_t j : use) {
    const auto row = g.site(j);
    double s = 0, typed = 0;
    for (const std::size_t i : samples)
      if (row[i] >= 0) {
        s += row[i];
        typed++;
      }
    if (typed == 0) continue;
    const double mu = s / typed, p = mu / 2;
    if (!(p > 0 && p < 1)) continue;
    const double scale = 1 / std::sqrt(p * (1 - p));
    kept.push_back(j);
    for (const std::size_t i : samples) M.push_back(row[i] >= 0 ? (row[i] - mu) * scale : 0.0);
  }
  const std::size_t m = kept.size();
  if (m == 0 || n < 2) throw std::invalid_argument("evolution::association::principal_components: polymorphic sites and two samples");
  calx::linalg::Matrix S(m);  // MᵀM (m × m)
  for (std::size_t a = 0; a < m; a++)
    for (std::size_t b = 0; b <= a; b++) {
      double s = 0;
      for (std::size_t i = 0; i < n; i++) s += M[a * n + i] * M[b * n + i];
      S(a, b) = S(b, a) = s;
    }
  const auto e = calx::linalg::symmetric_eigen(S);
  Components out;
  out.samples = n;
  out.sites_used = m;
  components = std::min({components, n, m});
  out.components = components;
  // The kinship X = M Mᵀ/m has the same non-zero eigenvalues as MᵀM/m; Patterson's scale makes the n
  // eigenvalues of X sum to n − 1, so it multiplies by (n − 1)/tr X, tr X = tr(MᵀM)/m.
  double trace = 0;
  for (std::size_t a = 0; a < m; a++) trace += S(a, a);
  const double scale = trace > 0 ? static_cast<double>(n - 1) / trace : 0;
  out.vectors.assign(n * components, 0.0);
  for (std::size_t c = 0; c < components; c++) {
    const std::size_t col = m - 1 - c;
    const double lambda = std::max(0.0, e.values[col]);
    out.values.push_back(lambda * scale);
    std::vector<double> u(n, 0.0);
    for (std::size_t a = 0; a < m; a++) {
      const double v = e.vectors(a, col);
      for (std::size_t i = 0; i < n; i++) u[i] += M[a * n + i] * v;
    }
    double norm = 0;
    for (const double x : u) norm += x * x;
    norm = std::sqrt(norm);
    double big = 0;
    for (const double x : u)
      if (std::fabs(x) > std::fabs(big)) big = x;
    const double sign = big < 0 ? -1 : 1;
    for (std::size_t i = 0; i < n; i++) out.vectors[i * components + c] = norm > 0 ? sign * u[i] / norm : 0;
  }
  return out;
}

}  // namespace evolution::association
