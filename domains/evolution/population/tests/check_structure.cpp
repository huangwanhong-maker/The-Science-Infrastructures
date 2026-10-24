// Checks of population structure: the simulator's Balding–Nichols frequencies have mean p and
// variance F p(1 − p); the kinship matrix equals M Mᵀ/m computed naively (Patterson et al. 2006,
// eq. 3, missing calls set to 0), and is the same on 1 and 8 threads; its eigenpairs satisfy
// X v = λ v, sum to n − 1 on Patterson's scale, and for K populations K − 1 of them stand out,
// with Tracy–Widom statistics above 3 while the next is not, and none stands out without
// structure; the leading components separate the populations (every sample nearest its own
// population's centroid); the EM of ancestry (Tang et al. 2005, accelerated by SQUAREM)
// recovers simulated proportions (the mean absolute error is stated), reaches a likelihood at
// least that of the truth, takes fewer steps accelerated than plain, and is the same on 1 and
// 8 threads.

#include <evolution/population/structure.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

using namespace evolution::population;

namespace {

std::string fmt(double x) { return calx::verification::exponential(x); }

// The mean absolute error of estimated ancestries against the truth, under the best matching of
// the estimated populations to the true ones.
double ancestry_error(const std::vector<double>& truth, const std::vector<double>& estimate, std::size_t n, std::size_t K) {
  std::vector<std::size_t> perm(K);
  std::iota(perm.begin(), perm.end(), 0);
  double best = INFINITY;
  do {
    double e = 0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t k = 0; k < K; k++) e += std::fabs(truth[i * K + k] - estimate[i * K + perm[k]]);
    best = std::min(best, e / static_cast<double>(n * K));
  } while (std::next_permutation(perm.begin(), perm.end()));
  return best;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::scheduler::ThreadPool one(1), eight(8);

  // ---- the simulator --------------------------------------------------------------------------------------
  {
    StructureSpec spec;
    spec.populations = 2;
    spec.sites = 20000;
    spec.fst = {0.1, 0.3};
    spec.unadmixed = {1, 1};
    const StructuredSample x = simulate_structure(spec, 1);
    for (std::size_t k = 0; k < 2; k++) {
      double mean = 0, ratio = 0, ratio2 = 0;
      for (std::size_t j = 0; j < spec.sites; j++) {
        const double p = x.ancestral[j], f = x.frequencies[j * 2 + k], r = (f - p) * (f - p) / (p * (1 - p));
        mean += f - p;
        ratio += r;
        ratio2 += r * r;
      }
      const double m = static_cast<double>(spec.sites), F = spec.fst[k];
      const double se = std::sqrt((ratio2 / m - (ratio / m) * (ratio / m)) / m);
      s.check("Balding–Nichols frequencies vary about p with variance F p(1 − p), F = " + fmt(F) + " (z-score over 20000 sites)", (ratio / m - F) / se, 4,
              "SE");
      s.check("and have mean p (mean of f − p)", mean / m, 0.005);
    }
    s.require("a structured sample is a function of its seed", simulate_structure(spec, 1).genotypes.dosage == x.genotypes.dosage &&
                                                                   simulate_structure(spec, 2).genotypes.dosage != x.genotypes.dosage);
  }

  // ---- the kinship matrix ----------------------------------------------------------------------------------
  {
    StructureSpec spec;
    spec.sites = 700;
    spec.unadmixed = {10, 10, 11};
    StructuredSample x = simulate_structure(spec, 3);
    Genotypes& g = x.genotypes;
    for (std::size_t j = 0; j < g.sites; j += 7) g.dosage[j * g.samples + (j % g.samples)] = MISSING;
    for (std::size_t i = 0; i < g.samples; i++) g.dosage[5 * g.samples + i] = 1;  // everyone heterozygous: p = ½, kept
    for (std::size_t i = 0; i < g.samples; i++) g.dosage[9 * g.samples + i] = 0;  // a monomorphic site, dropped
    // Naively: M (samples × sites), X = M Mᵀ / m.
    const std::size_t n = g.samples;
    std::vector<std::vector<double>> M;
    for (std::size_t j = 0; j < g.sites; j++) {
      double sum = 0, typed = 0;
      for (std::size_t i = 0; i < n; i++)
        if (g(i, j) >= 0) {
          sum += g(i, j);
          typed++;
        }
      const double mu = sum / typed, p = mu / 2;
      if (p <= 0 || p >= 1) continue;
      std::vector<double> col(n);
      for (std::size_t i = 0; i < n; i++) col[i] = g(i, j) >= 0 ? (g(i, j) - mu) / std::sqrt(p * (1 - p)) : 0;
      M.push_back(col);
    }
    const Kinship k1 = kinship(g, &one), k8 = kinship(g, &eight);
    double worst = 0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t c = 0; c < n; c++) {
        double v = 0;
        for (const auto& col : M) v += col[i] * col[c];
        worst = std::max(worst, std::fabs(v / static_cast<double>(M.size()) - k1.x(i, c)));
      }
    s.require("the kinship matrix uses the sites with 0 < p < 1", k1.sites_used == M.size(), std::to_string(k1.sites_used) + " of 700");
    s.check("it equals M Mᵀ/m computed naively, with missing calls as 0 (largest difference)", worst, 1e-12);
    s.require("and is the same on 1 and 8 threads", k1.x.a == k8.x.a);
  }

  // ---- principal components ----------------------------------------------------------------------------
  {
    StructureSpec spec;
    spec.sites = 5000;
    spec.fst = {0.05};
    spec.unadmixed = {40, 40, 40};
    const StructuredSample x = simulate_structure(spec, 4);
    const auto t0 = std::chrono::steady_clock::now();
    const Pca p = pca(x.genotypes, 4, &eight);
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    const std::size_t n = p.samples;
    double total = 0;
    for (const double v : p.all_values) total += v;
    s.check("the eigenvalues sum to n − 1 on Patterson's scale", total - static_cast<double>(n - 1), 1e-9);
    // X v = λ v, on X's own scale.
    const Kinship k = kinship(x.genotypes, &eight);
    double residual = 0;
    for (std::size_t c = 0; c < p.components; c++) {
      const double lambda = p.values[c] * p.trace / static_cast<double>(n - 1);
      for (std::size_t i = 0; i < n; i++) {
        double xv = 0;
        for (std::size_t u = 0; u < n; u++) xv += k.x(i, u) * p.vectors[u * p.components + c];
        residual = std::max(residual, std::fabs(xv - lambda * p.vectors[i * p.components + c]));
      }
    }
    s.check("the components are eigenvectors of the kinship matrix (largest |Xv − λv|)", residual, 1e-10);
    s.require("three populations give two eigenvalues that stand out (λ₂ > 2λ₃)", p.values[1] > 2 * p.values[2],
              fmt(p.values[0]) + ", " + fmt(p.values[1]) + ", " + fmt(p.values[2]) + ", " + fmt(p.values[3]));
    s.require("their Tracy–Widom statistics are above 3, and the third's is not (Patterson et al. 2006)",
              p.tracy_widom[0] > 3 && p.tracy_widom[1] > 3 && p.tracy_widom[2] < 3,
              fmt(p.tracy_widom[0]) + ", " + fmt(p.tracy_widom[1]) + ", " + fmt(p.tracy_widom[2]));
    // Nearest centroid on the first two components.
    double cx[3] = {0, 0, 0}, cy[3] = {0, 0, 0}, count[3] = {0, 0, 0};
    for (std::size_t i = 0; i < n; i++) {
      cx[x.population[i]] += p.vectors[i * 4];
      cy[x.population[i]] += p.vectors[i * 4 + 1];
      count[x.population[i]]++;
    }
    std::size_t right = 0;
    for (std::size_t i = 0; i < n; i++) {
      std::size_t best = 0;
      double bd = INFINITY;
      for (std::size_t c = 0; c < 3; c++) {
        const double dx = p.vectors[i * 4] - cx[c] / count[c], dy = p.vectors[i * 4 + 1] - cy[c] / count[c];
        if (dx * dx + dy * dy < bd) {
          bd = dx * dx + dy * dy;
          best = c;
        }
      }
      right += best == x.population[i];
    }
    s.require("the first two components separate the populations: every sample nearest its population's centroid", right == n,
              std::to_string(right) + " of " + std::to_string(n) + "; 120 samples × 5000 sites in " + fmt(ms) + " ms");

    StructureSpec none;
    none.populations = 1;
    none.sites = 5000;
    none.unadmixed = {120};
    const Pca q = pca(simulate_structure(none, 5).genotypes, 2, &eight);
    s.require("without structure no component stands out (Tracy–Widom statistic below 3)", q.tracy_widom[0] < 3, fmt(q.tracy_widom[0]));
    // Under no structure x is about TW₁: mean −1.2065, standard deviation 1.2680, 5% above
    // 0.9793 (Tracy and Widom 1996; the tables of Patterson et al. 2006, Table 1). Over 200
    // unstructured samples of 60 individuals and 2000 sites:
    double sum = 0;
    std::size_t above = 0;
    for (std::uint64_t r = 0; r < 200; r++) {
      StructureSpec small;
      small.populations = 1;
      small.sites = 2000;
      small.unadmixed = {60};
      const double x0 = pca(simulate_structure(small, 100 + r).genotypes, 1).tracy_widom[0];
      sum += x0;
      above += x0 > 0.9793;
    }
    s.check("the statistic's mean without structure is TW₁'s, −1.2065 (absolute error)", sum / 200 + 1.2065, 0.4);
    s.require("and it exceeds TW₁'s 95% point in about 5% of samples (at most 10%)", above <= 20, std::to_string(above) + " of 200");
  }

  // ---- ancestry by EM -------------------------------------------------------------------------------------
  {
    StructureSpec spec;
    spec.sites = 4000;
    spec.fst = {0.15};
    spec.unadmixed = {15, 15, 15};
    spec.admixed = 90;
    spec.alpha = 0.5;
    const StructuredSample x = simulate_structure(spec, 6);
    const std::size_t n = x.genotypes.samples, K = 3;
    AdmixtureOptions o;
    o.populations = K;
    o.seed = 7;
    const auto t0 = std::chrono::steady_clock::now();
    const Admixture a = admixture(x.genotypes, o, &eight);
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    const double error = ancestry_error(x.ancestry, a.q, n, K);
    const double truth = admixture_log_likelihood(x.genotypes, K, x.ancestry, x.frequencies);
    s.require("the accelerated EM converges", a.converged, std::to_string(a.iterations) + " EM steps in " + fmt(ms) + " ms");
    s.check("it recovers the simulated ancestry proportions (mean absolute error, 135 samples, 4000 sites, K = 3)", error, 0.03);
    s.require("its likelihood is at least the truth's", a.log_likelihood >= truth, fmt(a.log_likelihood) + " ≥ " + fmt(truth));
    bool simplex = true;
    for (std::size_t i = 0; i < n; i++) {
      double sum = 0;
      for (std::size_t k = 0; k < K; k++) {
        sum += a.q[i * K + k];
        simplex = simplex && a.q[i * K + k] >= 0;
      }
      simplex = simplex && std::fabs(sum - 1) < 1e-9;
    }
    bool rising = true;
    for (std::size_t c = 1; c < a.history.size(); c++) rising = rising && a.history[c] >= a.history[c - 1] - 1e-9;
    s.require("ancestries stay on the simplex, and the log-likelihood never falls", simplex && rising);
    AdmixtureOptions plain = o;
    plain.accelerate = false;
    const Admixture b = admixture(x.genotypes, plain, &eight);
    s.require("SQUAREM takes fewer EM steps than the plain EM to the same tolerance", a.iterations < b.iterations,
              std::to_string(a.iterations) + " < " + std::to_string(b.iterations) + "; plain error " + fmt(ancestry_error(x.ancestry, b.q, n, K)));
    // Thirty steps on one thread and on eight (the whole run on one thread would take seconds).
    AdmixtureOptions short_run = o;
    short_run.max_iterations = 30;
    const Admixture c = admixture(x.genotypes, short_run, &one), d = admixture(x.genotypes, short_run, &eight);
    s.require("the estimate is the same on 1 and 8 threads", c.q == d.q && c.f == d.f && c.log_likelihood == d.log_likelihood && c.iterations == d.iterations);
  }
  return s.done();
}
