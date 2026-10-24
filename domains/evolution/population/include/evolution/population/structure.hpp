// Population structure: principal components of a genotype matrix, ancestry proportions by
// maximum likelihood, and a simulator of structured and admixed samples to check them on.
//
// - Principal components (Patterson, Price and Reich 2006, PLoS Genetics 2, e190; their m
//   individuals and n markers are our n samples and m sites). Each site j is centred by its
//   mean dosage μ_j and scaled by √(p_j(1 − p_j)), p_j = μ_j/2 (their eq. 3); a missing call is
//   set to 0 after centring (their §"Missing data"); monomorphic sites are dropped. The samples'
//   covariance X = M Mᵀ / m, an n × n kinship matrix, is decomposed by calx's symmetric
//   eigensolver; its leading eigenvectors are the axes of variation, and the samples'
//   coordinates are the eigenvectors' entries. The eigenvalues are reported on Patterson's
//   scale, summing to n − 1 (one eigenvalue is a structural zero, since every centred column
//   sums to zero). The Tracy–Widom statistic of each component (their eqs. 5–10) is
//   x = (L − μ)/σ with L = M λ_k / Σ_{i≥k} λ_i over the M = n − k eigenvalues from λ_k on,
//   μ = (√(n′ − 1) + √M)²/n′, σ = ((√(n′ − 1) + √M)/n′)(1/√(n′ − 1) + 1/√M)^{1/3}, and the
//   effective number of markers n′ = (M + 2)(Σλ)² / (MΣλ² − (Σλ)²) (their eq. 10, whose m
//   individuals give the M = m − 1 eigenvalues of the centred matrix); without structure x is
//   approximately Tracy–Widom (TW₁) distributed, and x above about 3 marks a significant axis
//   (the p-value needs the TW₁ distribution function, which calx does not have yet).
// - Ancestry proportions. Each individual i draws each of its 2 alleles at site j from one of K
//   ancestral populations, population k with probability q_ik, and the allele is alternate with
//   probability f_kj; so g_ij ~ Binomial(2, π_ij), π_ij = Σ_k q_ik f_kj. The likelihood is
//   maximised by the EM algorithm of Tang, Peng, Wang and Risch (2005, Genetic Epidemiology 28,
//   289–301, FRAPPE): the latent population of each allele is the missing datum; given the
//   current (q, f), an alternate allele came from k with probability q_ik f_kj / π_ij and a
//   reference allele with q_ik (1 − f_kj)/(1 − π_ij); q is re-estimated as the fraction of i's
//   alleles from each k, f as the fraction of alternate alleles among those from k. Each step
//   raises the likelihood. The EM is accelerated by SQUAREM (Varadhan and Roland 2008, Scand.
//   J. Statist. 35, 335–353, scheme S3), the role ADMIXTURE's quasi-Newton steps play
//   (Alexander, Novembre and Lange 2009, Genome Research 19, 1655–1664): from θ₀, two EM steps
//   θ₁, θ₂, r = θ₁ − θ₀, v = θ₂ − 2θ₁ + θ₀, α = −‖r‖/‖v‖ (at most −1), θ′ = θ₀ − 2αr + α²v,
//   projected onto the parameter space and followed by one EM step; when that loses likelihood,
//   θ₂ is taken instead, so the likelihood never falls. The starting point is random (seeded);
//   the likelihood has local maxima and label switching, so `restarts` runs from several seeds
//   keep the best.
// - A simulator of structure: ancestral frequencies p_j uniform on [low, high]; K populations
//   diverged from it by the Balding–Nichols model (Balding and Nichols 1995, Genetica 96, 3–12),
//   f_kj ~ Beta(p_j(1 − F_k)/F_k, (1 − p_j)(1 − F_k)/F_k), so that E[f_kj] = p_j and
//   Var[f_kj] = F_k p_j(1 − p_j); unadmixed individuals of each population, and admixed ones with
//   ancestry q_i ~ Dirichlet(α, …, α); genotypes g_ij ~ Binomial(2, Σ_k q_ik f_kj). The draws are
//   gamma, beta and binomial variates from calx streams keyed by site and by individual.
//
// Complexity: the kinship matrix O(n² m), in tiles of sites on vector lanes and in parallel
// over rows (each entry summed by one thread in a fixed order, so the same on any number of
// threads); the eigenproblem O(n³). One EM step O(n m K), in parallel over fixed shards of
// sites whose partial sums are merged in shard order; SQUAREM spends three steps per cycle.

#pragma once

#include <evolution/population/genotypes.hpp>

#include <calx/cpu/lanes.hpp>
#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigen.hpp>
#include <calx/random/stream.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace evolution::population {

inline constexpr std::uint32_t ADMIXTURE = 0x41444d58, STRUCTURE = 0x42414c44;

// ---- the kinship matrix and its principal components -------------------------------------------------------

struct Kinship {
  calx::linalg::Matrix x;    // M Mᵀ / m, n × n
  std::size_t sites_used = 0;  // polymorphic sites
};

namespace detail {

// Σ_t a[t] b[t] over T values (T a multiple of the lane width), in a fixed order.
template <class V>
inline void dot4(const double* a, const double* b0, const double* b1, const double* b2, const double* b3, std::size_t T, double out[4]) {
  V s0 = V::zero(), s1 = V::zero(), s2 = V::zero(), s3 = V::zero();
  for (std::size_t t = 0; t < T; t += V::width) {
    const V x = V::load(a + t);
    s0 = s0 + x * V::load(b0 + t);
    s1 = s1 + x * V::load(b1 + t);
    s2 = s2 + x * V::load(b2 + t);
    s3 = s3 + x * V::load(b3 + t);
  }
  out[0] = s0.hsum();
  out[1] = s1.hsum();
  out[2] = s2.hsum();
  out[3] = s3.hsum();
}

}  // namespace detail

inline Kinship kinship(const Genotypes& g, calx::scheduler::ThreadPool* pool = nullptr) {
  using V = calx::cpu::Lanes<double>;
  const std::size_t n = g.samples, T = 256;  // a tile of sites: n × T doubles
  static_assert(256 % V::width == 0);
  Kinship k;
  k.x = calx::linalg::Matrix(n);
  std::vector<double> tile(n * T, 0.0);
  std::vector<double> X(n * n, 0.0);  // lower triangle accumulated, row-major
  std::size_t filled = 0;
  const auto flush = [&] {
    if (filled == 0) return;
    // Columns beyond `filled` are zero, so they add nothing.
    const auto rows = [&](std::size_t begin, std::size_t end, std::size_t) {
      double d[4];
      for (std::size_t i = begin; i < end; i++) {
        const double* a = &tile[i * T];
        std::size_t c = 0;
        for (; c + 4 <= i + 1; c += 4) {
          detail::dot4<V>(a, &tile[c * T], &tile[(c + 1) * T], &tile[(c + 2) * T], &tile[(c + 3) * T], T, d);
          for (std::size_t u = 0; u < 4; u++) X[i * n + c + u] += d[u];
        }
        for (; c <= i; c++) {
          detail::dot4<V>(a, &tile[c * T], &tile[c * T], &tile[c * T], &tile[c * T], T, d);
          X[i * n + c] += d[0];
        }
      }
    };
    if (pool) calx::scheduler::parallel_for(*pool, n, 8, rows);
    else rows(0, n, 0);
    std::fill(tile.begin(), tile.end(), 0.0);
    filled = 0;
  };
  for (std::size_t j = 0; j < g.sites; j++) {
    const auto row = g.site(j);
    double sum = 0, typed = 0;
    for (const std::int8_t x : row)
      if (x >= 0) {
        sum += x;
        typed++;
      }
    if (typed == 0) continue;
    const double mu = sum / typed, p = mu / 2;
    if (!(p > 0 && p < 1)) continue;
    const double scale = 1 / std::sqrt(p * (1 - p));
    for (std::size_t i = 0; i < n; i++) tile[i * T + filled] = row[i] >= 0 ? (row[i] - mu) * scale : 0.0;
    k.sites_used++;
    if (++filled == T) flush();
  }
  flush();
  const double m = static_cast<double>(std::max<std::size_t>(1, k.sites_used));
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t c = 0; c <= i; c++) k.x(i, c) = k.x(c, i) = X[i * n + c] / m;
  return k;
}

struct Pca {
  std::size_t samples = 0, components = 0, sites_used = 0;
  std::vector<double> values;        // the leading eigenvalues, descending, on Patterson's scale (all n sum to n − 1)
  std::vector<double> vectors;       // samples × components: each sample's coordinate on each axis (unit eigenvectors)
  std::vector<double> all_values;    // every eigenvalue, descending, on the same scale
  std::vector<double> tracy_widom;   // x of each leading component
  std::vector<double> effective_markers;  // n′ of each
  double trace = 0;                  // of X = M Mᵀ / m
};

// Patterson's statistics for the eigenvalues λ_k ≥ … of `values` (descending), from k on.
inline void tracy_widom_statistic(const std::vector<double>& values, std::size_t k, double& x, double& effective) {
  // The structural zero (the last eigenvalue) is left out: M = n − 1 − k eigenvalues from λ_k.
  const std::size_t count = values.size() - 1;
  if (k + 2 > count) {
    x = effective = std::numeric_limits<double>::quiet_NaN();
    return;
  }
  double s = 0, s2 = 0;
  for (std::size_t i = k; i < count; i++) {
    s += values[i];
    s2 += values[i] * values[i];
  }
  // Eq. 10 with m the number of individuals the M eigenvalues come from, M + 1 (one dimension
  // is lost to centring): n′ = (m + 1)(Σλ)² / ((m − 1)Σλ² − (Σλ)²), the moment estimate that
  // makes Σλ²/(Σλ)² that of a Wishart matrix with n′ degrees of freedom.
  const double M = static_cast<double>(count - k);
  effective = (M + 2) * s * s / (M * s2 - s * s);
  const double L = M * values[k] / s, a = std::sqrt(effective - 1), b = std::sqrt(M);
  const double mu = (a + b) * (a + b) / effective, sigma = (a + b) / effective * std::cbrt(1 / a + 1 / b);
  x = (L - mu) / sigma;
}

inline Pca pca(const Genotypes& g, std::size_t components, calx::scheduler::ThreadPool* pool = nullptr) {
  const std::size_t n = g.samples;
  if (n < 2) throw std::invalid_argument("evolution::population::pca: at least two samples");
  components = std::min(components, n);
  const Kinship k = kinship(g, pool);
  const calx::linalg::SymmetricEigen e = calx::linalg::symmetric_eigen(k.x);
  Pca out;
  out.samples = n;
  out.components = components;
  out.sites_used = k.sites_used;
  for (std::size_t i = 0; i < n; i++) out.trace += k.x(i, i);
  const double scale = out.trace > 0 ? static_cast<double>(n - 1) / out.trace : 0;
  out.all_values.resize(n);
  for (std::size_t i = 0; i < n; i++) out.all_values[i] = std::max(0.0, e.values[n - 1 - i]) * scale;
  out.values.assign(out.all_values.begin(), out.all_values.begin() + static_cast<std::ptrdiff_t>(components));
  out.vectors.resize(n * components);
  for (std::size_t c = 0; c < components; c++) {
    // A sign fixed by convention: the largest entry positive.
    const std::size_t col = n - 1 - c;
    double big = 0;
    for (std::size_t i = 0; i < n; i++)
      if (std::fabs(e.vectors(i, col)) > std::fabs(big)) big = e.vectors(i, col);
    const double sign = big < 0 ? -1 : 1;
    for (std::size_t i = 0; i < n; i++) out.vectors[i * components + c] = sign * e.vectors(i, col);
  }
  out.tracy_widom.resize(components);
  out.effective_markers.resize(components);
  for (std::size_t c = 0; c < components; c++) tracy_widom_statistic(out.all_values, c, out.tracy_widom[c], out.effective_markers[c]);
  return out;
}

// ---- ancestry by maximum likelihood --------------------------------------------------------------------------

struct AdmixtureOptions {
  std::size_t populations = 2;       // K
  std::size_t max_iterations = 2000; // EM steps
  double tolerance = 1e-4;           // stop when a cycle raises the log-likelihood by less
  bool accelerate = true;            // SQUAREM
  std::size_t restarts = 1;          // random starts (seeds seed, seed + 1, …); the best is kept
  std::uint64_t seed = 1;
};

struct Admixture {
  std::size_t samples = 0, populations = 0, sites = 0;
  std::vector<double> q;          // samples × K: ancestry proportions, each row summing to 1
  std::vector<double> f;          // sites × K: the alternate allele's frequency in each population
  double log_likelihood = 0;
  std::size_t iterations = 0;     // EM steps taken
  bool converged = false;
  std::vector<double> history;    // the log-likelihood after each cycle
};

namespace detail {

inline constexpr double FREQUENCY_FLOOR = 1e-6, ANCESTRY_FLOOR = 1e-9;

// One EM step from (q, f) to (q_next, f_next); returns the log-likelihood of (q, f). Shards of
// 256 sites sum their own contributions to q, merged in shard order.
inline double em_step(const Genotypes& g, std::size_t K, const std::vector<double>& q, const std::vector<double>& f, std::vector<double>& q_next,
                      std::vector<double>& f_next, const std::vector<double>& typed, std::vector<double>& partial, std::vector<double>& ll_partial,
                      calx::scheduler::ThreadPool* pool) {
  const std::size_t n = g.samples, m = g.sites, shard = 256, shards = (m + shard - 1) / shard;
  partial.assign(shards * n * K, 0.0);
  ll_partial.assign(shards, 0.0);
  const auto work = [&](std::size_t begin, std::size_t end, std::size_t) {
    std::vector<double> a(K), b(K), A(K), B(K);
    for (std::size_t s = begin; s < end; s++) {
      double* Q = &partial[s * n * K];
      double ll = 0;
      for (std::size_t j = s * shard; j < std::min(m, (s + 1) * shard); j++) {
        const auto row = g.site(j);
        const double* fj = &f[j * K];
        std::fill(A.begin(), A.end(), 0.0);
        std::fill(B.begin(), B.end(), 0.0);
        for (std::size_t i = 0; i < n; i++) {
          const std::int8_t x = row[i];
          if (x < 0) continue;
          const double* qi = &q[i * K];
          double pi = 0;
          for (std::size_t k = 0; k < K; k++) pi += qi[k] * fj[k];
          pi = std::clamp(pi, FREQUENCY_FLOOR, 1 - FREQUENCY_FLOOR);
          // One logarithm and one division per genotype.
          const double alt = x, ref = 2 - x, inv = 1 / (pi * (1 - pi));
          ll += x == 0 ? 2 * std::log1p(-pi) : x == 2 ? 2 * std::log(pi) : std::log(pi * (1 - pi));
          const double ua = alt * (1 - pi) * inv, ub = ref * pi * inv;
          double* Qi = &Q[i * K];
          for (std::size_t k = 0; k < K; k++) {
            const double ak = ua * qi[k] * fj[k], bk = ub * qi[k] * (1 - fj[k]);
            A[k] += ak;
            B[k] += bk;
            Qi[k] += ak + bk;
          }
        }
        for (std::size_t k = 0; k < K; k++)
          f_next[j * K + k] = A[k] + B[k] > 0 ? std::clamp(A[k] / (A[k] + B[k]), FREQUENCY_FLOOR, 1 - FREQUENCY_FLOOR) : fj[k];
      }
      ll_partial[s] = ll;
    }
  };
  if (pool) calx::scheduler::parallel_for(*pool, shards, 1, work);
  else work(0, shards, 0);
  double ll = 0;
  std::fill(q_next.begin(), q_next.end(), 0.0);
  for (std::size_t s = 0; s < shards; s++) {
    ll += ll_partial[s];
    const double* Q = &partial[s * n * K];
    for (std::size_t x = 0; x < n * K; x++) q_next[x] += Q[x];
  }
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t k = 0; k < K; k++) q_next[i * K + k] = typed[i] > 0 ? q_next[i * K + k] / (2 * typed[i]) : q[i * K + k];
  return ll;
}

// Onto the parameter space: frequencies in [floor, 1 − floor], ancestries on the simplex.
inline void project(std::size_t K, std::vector<double>& q, std::vector<double>& f) {
  for (double& x : f) x = std::clamp(x, FREQUENCY_FLOOR, 1 - FREQUENCY_FLOOR);
  for (std::size_t i = 0; i < q.size() / K; i++) {
    double s = 0;
    for (std::size_t k = 0; k < K; k++) s += q[i * K + k] = std::max(ANCESTRY_FLOOR, q[i * K + k]);
    for (std::size_t k = 0; k < K; k++) q[i * K + k] /= s;
  }
}

inline Admixture admixture_from(const Genotypes& g, const AdmixtureOptions& o, std::uint64_t seed, calx::scheduler::ThreadPool* pool) {
  const std::size_t n = g.samples, m = g.sites, K = o.populations;
  Admixture out;
  out.samples = n;
  out.populations = K;
  out.sites = m;
  std::vector<double> q(n * K), f(m * K);
  for (std::size_t i = 0; i < n; i++) {
    calx::random::Stream r(seed, ADMIXTURE, 0, static_cast<std::uint32_t>(i));
    double s = 0;
    for (std::size_t k = 0; k < K; k++) s += q[i * K + k] = r.gamma(1);
    for (std::size_t k = 0; k < K; k++) q[i * K + k] /= s;
  }
  for (std::size_t j = 0; j < m; j++) {
    calx::random::Stream r(seed, ADMIXTURE, 1, static_cast<std::uint32_t>(j));
    for (std::size_t k = 0; k < K; k++) f[j * K + k] = r.uniform(0.05, 0.95);
  }
  std::vector<double> typed(n, 0.0);
  for (std::size_t j = 0; j < m; j++) {
    const auto row = g.site(j);
    for (std::size_t i = 0; i < n; i++) typed[i] += row[i] >= 0;
  }
  std::vector<double> q1(n * K), f1(m * K), q2(n * K), f2(m * K), q3(n * K), f3(m * K), partial, ll_partial;
  const auto step = [&](const std::vector<double>& qa, const std::vector<double>& fa, std::vector<double>& qb, std::vector<double>& fb) {
    out.iterations++;
    return em_step(g, K, qa, fa, qb, fb, typed, partial, ll_partial, pool);
  };
  double previous = -std::numeric_limits<double>::infinity();
  while (out.iterations < o.max_iterations) {
    double ll0;
    if (!o.accelerate) {
      ll0 = step(q, f, q1, f1);
      q.swap(q1);
      f.swap(f1);
    } else {
      ll0 = step(q, f, q1, f1);
      step(q1, f1, q2, f2);
      // r = θ₁ − θ₀, v = θ₂ − 2θ₁ + θ₀ over both blocks of parameters.
      double rr = 0, vv = 0;
      for (std::size_t x = 0; x < q.size(); x++) {
        const double r = q1[x] - q[x], v = q2[x] - 2 * q1[x] + q[x];
        rr += r * r;
        vv += v * v;
      }
      for (std::size_t x = 0; x < f.size(); x++) {
        const double r = f1[x] - f[x], v = f2[x] - 2 * f1[x] + f[x];
        rr += r * r;
        vv += v * v;
      }
      double alpha = vv > 0 ? -std::sqrt(rr / vv) : -1;
      alpha = std::min(alpha, -1.0);
      for (std::size_t x = 0; x < q.size(); x++) q3[x] = q[x] - 2 * alpha * (q1[x] - q[x]) + alpha * alpha * (q2[x] - 2 * q1[x] + q[x]);
      for (std::size_t x = 0; x < f.size(); x++) f3[x] = f[x] - 2 * alpha * (f1[x] - f[x]) + alpha * alpha * (f2[x] - 2 * f1[x] + f[x]);
      project(K, q3, f3);
      // The stabilising step: θ′ → F(θ′), whose likelihood is that of θ′.
      const double ll3 = step(q3, f3, q1, f1);
      if (std::isfinite(ll3) && ll3 >= ll0) {
        q.swap(q1);
        f.swap(f1);
      } else {
        q.swap(q2);
        f.swap(f2);
      }
    }
    out.history.push_back(ll0);
    if (ll0 - previous < o.tolerance) {
      out.converged = true;
      break;
    }
    previous = ll0;
  }
  out.log_likelihood = em_step(g, K, q, f, q1, f1, typed, partial, ll_partial, pool);
  out.q = std::move(q);
  out.f = std::move(f);
  return out;
}

}  // namespace detail

// The log-likelihood Σ_ij [g log π + (2 − g) log(1 − π)] of ancestries q (samples × K) and
// frequencies f (sites × K).
inline double admixture_log_likelihood(const Genotypes& g, std::size_t K, const std::vector<double>& q, const std::vector<double>& f) {
  double ll = 0;
  for (std::size_t j = 0; j < g.sites; j++) {
    const auto row = g.site(j);
    for (std::size_t i = 0; i < g.samples; i++) {
      if (row[i] < 0) continue;
      double pi = 0;
      for (std::size_t k = 0; k < K; k++) pi += q[i * K + k] * f[j * K + k];
      pi = std::clamp(pi, detail::FREQUENCY_FLOOR, 1 - detail::FREQUENCY_FLOOR);
      ll += row[i] * std::log(pi) + (2 - row[i]) * std::log1p(-pi);
    }
  }
  return ll;
}

inline Admixture admixture(const Genotypes& g, const AdmixtureOptions& o, calx::scheduler::ThreadPool* pool = nullptr) {
  if (o.populations < 1 || g.samples == 0 || g.sites == 0) throw std::invalid_argument("evolution::population::admixture: K >= 1 and a non-empty matrix");
  Admixture best;
  for (std::size_t r = 0; r < std::max<std::size_t>(1, o.restarts); r++) {
    Admixture a = detail::admixture_from(g, o, o.seed + r, pool);
    if (r == 0 || a.log_likelihood > best.log_likelihood) best = std::move(a);
  }
  return best;
}

// ---- structured samples -------------------------------------------------------------------------------

struct StructureSpec {
  std::size_t populations = 3;          // K
  std::size_t sites = 5000;
  std::vector<double> fst = {0.1};      // Balding–Nichols F of each population (one value serves all)
  std::vector<std::size_t> unadmixed = {20, 20, 20};  // individuals wholly of each population
  std::size_t admixed = 0;              // individuals with Dirichlet ancestry
  double alpha = 1;                     // its concentration in each component
  double low = 0.05, high = 0.95;       // the ancestral frequency, uniform between
  std::uint64_t spacing = 1000;         // bases between consecutive sites, on one chromosome
};

struct StructuredSample {
  Genotypes genotypes;
  std::vector<double> ancestry;            // samples × K
  std::vector<double> frequencies;         // sites × K
  std::vector<double> ancestral;           // of each site
  std::vector<std::uint32_t> population;   // of each unadmixed individual; of an admixed one, its largest component
  std::vector<std::uint8_t> admixed;       // whether each individual is admixed
};

inline StructuredSample simulate_structure(const StructureSpec& spec, std::uint64_t seed) {
  const std::size_t K = spec.populations, m = spec.sites;
  if (K == 0 || spec.fst.empty() || (spec.fst.size() != 1 && spec.fst.size() != K) || spec.unadmixed.size() > K || !(spec.low >= 0) ||
      !(spec.high <= 1) || !(spec.low <= spec.high) || !(spec.alpha > 0))
    throw std::invalid_argument("evolution::population::simulate_structure: K populations, an F for each (or one), frequencies in [0, 1]");
  StructuredSample out;
  std::size_t n = spec.admixed;
  for (const std::size_t u : spec.unadmixed) n += u;
  out.ancestry.assign(n * K, 0.0);
  out.population.resize(n);
  out.admixed.assign(n, 0);
  std::size_t i = 0;
  for (std::size_t k = 0; k < spec.unadmixed.size(); k++)
    for (std::size_t c = 0; c < spec.unadmixed[k]; c++, i++) {
      out.ancestry[i * K + k] = 1;
      out.population[i] = static_cast<std::uint32_t>(k);
    }
  for (; i < n; i++) {
    calx::random::Stream r(seed, STRUCTURE, 1, static_cast<std::uint32_t>(i));
    double s = 0;
    for (std::size_t k = 0; k < K; k++) s += out.ancestry[i * K + k] = r.gamma(spec.alpha);
    std::size_t top = 0;
    for (std::size_t k = 0; k < K; k++) {
      out.ancestry[i * K + k] /= s;
      if (out.ancestry[i * K + k] > out.ancestry[i * K + top]) top = k;
    }
    out.population[i] = static_cast<std::uint32_t>(top);
    out.admixed[i] = 1;
  }
  out.ancestral.resize(m);
  out.frequencies.resize(m * K);
  std::vector<std::int8_t> dosage(m * n);
  std::vector<std::uint64_t> positions(m);
  for (std::size_t j = 0; j < m; j++) {
    calx::random::Stream r(seed, STRUCTURE, 0, static_cast<std::uint32_t>(j));
    const double p = r.uniform(spec.low, spec.high);
    out.ancestral[j] = p;
    for (std::size_t k = 0; k < K; k++) {
      const double F = spec.fst.size() == 1 ? spec.fst[0] : spec.fst[k];
      out.frequencies[j * K + k] = F > 0 ? r.beta(p * (1 - F) / F, (1 - p) * (1 - F) / F) : p;
    }
    calx::random::Stream d(seed, STRUCTURE, 2, static_cast<std::uint32_t>(j));
    for (std::size_t x = 0; x < n; x++) {
      double pi = 0;
      for (std::size_t k = 0; k < K; k++) pi += out.ancestry[x * K + k] * out.frequencies[j * K + k];
      dosage[j * n + x] = static_cast<std::int8_t>(d.binomial(2, pi));
    }
    positions[j] = j * spec.spacing;
  }
  out.genotypes = genotypes_from_dosages(n, std::move(dosage), std::move(positions));
  return out;
}

}  // namespace evolution::population
