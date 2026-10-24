// Checks of differential expression: Benjamini and Hochberg's adjustment on a hand-worked list
// and on the fifteen p-values of their own example (Benjamini and Hochberg 1995, §4, from
// Neuhaus et al. 1992: four rejections at q = 0.05); the median of ratios recovering exact and
// noisy known scalings; the negative binomial's probabilities summing to 1 and equal to the
// closed form C(y + r − 1, y) p^r (1 − p)^y at integer r; its gamma–Poisson draws with mean μ and
// variance μ + αμ²; IRLS against the closed form of a two-group design at equal size factors,
// and its score equations at zero with unequal ones; Cox–Reid dispersion estimates at a known
// dispersion (the median of 400 genes, reported) and a stationary point of their objective;
// the trend's coefficients; Wald p-values uniform under the null (Kolmogorov–Smirnov distance
// reported, limit stated) and simulated two-fold changes detected (power and the observed false
// discovery rate, reported); the same results on 1 and 8 threads.

#include <genomics/expression/differential.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

using namespace genomics::expression;

namespace {

// sup |F_n(x) − x| of a sample against the uniform law on (0, 1).
double ks_uniform(std::vector<double> p) {
  std::sort(p.begin(), p.end());
  const double n = static_cast<double>(p.size());
  double d = 0;
  for (std::size_t i = 0; i < p.size(); i++) d = std::max({d, static_cast<double>(i + 1) / n - p[i], p[i] - static_cast<double>(i) / n});
  return d;
}

struct Simulation {
  CountMatrix counts;
  std::vector<double> base, dispersion, lfc, size;
  std::vector<int> group;
};

// G genes, m/2 samples in each group; a fraction `de` of the genes changed by ±log 2.
Simulation simulate(std::size_t G, std::size_t m, double de, std::uint64_t seed) {
  Simulation sim;
  calx::random::Stream r(seed, 0x44455331, 0, 0);
  for (std::size_t j = 0; j < m; j++) {
    sim.group.push_back(j < m / 2 ? 0 : 1);
    sim.size.push_back(0.6 + r.uniform());
  }
  for (std::size_t g = 0; g < G; g++) {
    const double b = std::max(5.0, r.lognormal(std::log(150.0), 1.2));
    sim.base.push_back(b);
    sim.dispersion.push_back((0.03 + 3.0 / b) * r.lognormal(0, 0.3));
    sim.lfc.push_back(r.uniform() < de ? (r.bernoulli(0.5) ? 1 : -1) * std::log(2.0) : 0.0);
  }
  sim.counts = simulate_counts(sim.base, sim.dispersion, sim.lfc, sim.group, sim.size, seed);
  return sim;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x44494646, 0, 0);
  calx::scheduler::ThreadPool one(1), eight(8);

  // ---- Benjamini–Hochberg ------------------------------------------------------------------------------------
  {
    // Sorted 0.005, 0.01, 0.03, 0.04, 0.2: m p/i = 0.025, 0.025, 0.05, 0.05, 0.2, already monotone.
    const auto q = benjamini_hochberg({0.01, 0.04, 0.03, 0.005, 0.20});
    const std::vector<double> want = {0.025, 0.05, 0.05, 0.025, 0.2};
    double err = 0;
    for (std::size_t i = 0; i < q.size(); i++) err = std::max(err, std::fabs(q[i] - want[i]));
    s.check("BH on 0.01, 0.04, 0.03, 0.005, 0.2 gives 0.025, 0.05, 0.05, 0.025, 0.2", err, 1e-15);
    const auto qn = benjamini_hochberg({0.01, std::numeric_limits<double>::quiet_NaN(), 0.04});
    s.require("a NaN p-value is left out of m and stays NaN: 0.01, NaN, 0.04 gives 0.02, NaN, 0.04",
              std::fabs(qn[0] - 0.02) < 1e-15 && std::isnan(qn[1]) && std::fabs(qn[2] - 0.04) < 1e-15);
    const std::vector<double> bh = {0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298, 0.0344, 0.0459, 0.3240, 0.4262, 0.5719, 0.6528, 0.7590, 1.000};
    const auto qb = benjamini_hochberg(bh);
    std::size_t rejected = 0;
    for (double x : qb) rejected += x <= 0.05;
    s.require("on Benjamini and Hochberg's fifteen p-values the procedure at q = 0.05 rejects four hypotheses, as their §4 reports", rejected == 4);
    s.check("the adjusted value of the fourth is 15 · 0.0095 / 4 = 0.035625", qb[3] - 0.035625, 1e-15);
    s.check("of the fifth, min over k ≥ 5 of 15 p_(k)/k = 15 · 0.0201 / 5 = 0.0603", qb[4] - 0.0603, 1e-15);
  }

  // ---- size factors ------------------------------------------------------------------------------------
  {
    const std::vector<double> scale = {1, 2, 0.5, 4, 1.25};
    CountMatrix k(300, scale.size());
    for (std::size_t g = 0; g < 300; g++) {
      const double base = 4.0 * static_cast<double>(1 + r.below(500));
      for (std::size_t j = 0; j < scale.size(); j++) k.at(g, j) = base * scale[j];
    }
    const auto sf = size_factors(k);
    double geo = 0;
    for (double x : scale) geo += std::log(x) / static_cast<double>(scale.size());
    double err = 0;
    for (std::size_t j = 0; j < scale.size(); j++) err = std::max(err, std::fabs(sf[j] / (scale[j] / std::exp(geo)) - 1));
    s.check("counts scaled exactly by s give the size factors s over their geometric mean", err, 1e-12);
    const Simulation sim = simulate(2000, 8, 0.1, 31);
    const auto est = size_factors(sim.counts);
    double ge = 0, gs = 0;
    for (std::size_t j = 0; j < 8; j++) { ge += std::log(est[j]) / 8; gs += std::log(sim.size[j]) / 8; }
    double worst = 0;
    for (std::size_t j = 0; j < 8; j++) worst = std::max(worst, std::fabs(est[j] / std::exp(ge) / (sim.size[j] / std::exp(gs)) - 1));
    s.check("from negative binomial counts of 2000 genes (10% changed), the size factors relative to their geometric mean", worst, 0.03);
  }

  // ---- the negative binomial ------------------------------------------------------------------------------------
  {
    const double mu = 7.3, alpha = 0.4;
    double total = 0;
    for (int y = 0; y < 3000; y++) {
      const double yy = y;
      total += std::exp(nb_log_likelihood(&yy, &mu, 1, alpha));
    }
    s.check("the probabilities of NB(7.3, 0.4) sum to 1", total - 1, 1e-10);
    // r = 3: P(y) = C(y + 2, 2) p³ (1 − p)^y, p = r / (r + μ).
    const double m2 = 5.0, r3 = 3.0, pp = r3 / (r3 + m2);
    double worst = 0;
    for (int y = 0; y <= 12; y++) {
      const double yy = y, closed = (y + 1.0) * (y + 2.0) / 2 * std::pow(pp, 3) * std::pow(1 - pp, y);
      worst = std::max(worst, std::fabs(std::exp(nb_log_likelihood(&yy, &m2, 1, 1 / r3)) / closed - 1));
    }
    s.check("at α = 1/3 the probabilities are C(y + 2, 2) p³ (1 − p)^y (y = 0 to 12)", worst, 1e-12);
    calx::random::Stream d(5, 0x4e424452, 0, 0);
    double sum = 0, sq = 0;
    const int n = 200000;
    for (int i = 0; i < n; i++) {
      const double x = static_cast<double>(nb_draw(d, 20, 0.25));
      sum += x;
      sq += x * x;
    }
    const double mean = sum / n, var = sq / n - mean * mean;
    s.check("gamma–Poisson draws of NB(20, 0.25) have mean 20", mean / 20 - 1, 0.005);
    s.check("and variance μ + αμ² = 120", var / 120 - 1, 0.03);
  }

  // ---- IRLS ------------------------------------------------------------------------------------
  {
    const Design X = Design::two_groups({0, 0, 0, 0, 1, 1, 1, 1});
    const std::vector<double> y = {12, 30, 18, 25, 60, 41, 77, 52}, ones(8, 1.0);
    const NbFit f = fit_glm(y.data(), ones, X, 0.1);
    const double a = (12 + 30 + 18 + 25) / 4.0, b = (60 + 41 + 77 + 52) / 4.0;
    s.check("at equal size factors IRLS gives β₀ = log of group 0's mean", f.beta[0] - std::log(a), 1e-6);
    s.check("and β₁ = the log of the ratio of the groups' means", f.beta[1] - std::log(b / a), 1e-6);
    const std::vector<double> sz = {0.5, 1.2, 0.8, 1.7, 0.9, 1.3, 0.6, 2.0};
    const NbFit g = fit_glm(y.data(), sz, X, 0.1);
    double worst = 0;
    for (std::size_t c = 0; c < 2; c++) {
      double score = 0;
      for (std::size_t j = 0; j < 8; j++) score += X(j, c) * (y[j] - g.mu[j]) / (1 + 0.1 * g.mu[j]);
      worst = std::max(worst, std::fabs(score));
    }
    s.check("with unequal size factors the score equations Σ x (y − μ)/(1 + αμ) vanish at the fit", worst, 1e-5);
    s.require("and it converged", g.converged && f.converged, std::to_string(g.iterations) + " iterations");
  }

  // ---- dispersions ------------------------------------------------------------------------------------
  {
    const std::size_t G = 400, m = 10;
    const double alpha = 0.1;
    std::vector<int> group;
    for (std::size_t j = 0; j < m; j++) group.push_back(j < m / 2 ? 0 : 1);
    const CountMatrix k = simulate_counts(std::vector<double>(G, 100.0), std::vector<double>(G, alpha), std::vector<double>(G, 0.0), group, std::vector<double>(m, 1.0), 41);
    const Design X = Design::two_groups(group);
    const std::vector<double> ones(m, 1.0);
    std::vector<double> est, est_ml;
    bool stationary = true;
    for (std::size_t g = 0; g < G; g++) {
      const double* y = k.row(g);
      NbFit f = fit_glm(y, ones, X, alpha);
      const double a = estimate_dispersion(y, f.mu, X, true);
      est.push_back(a);
      est_ml.push_back(estimate_dispersion(y, f.mu, X, false));
      const double th = std::log(a);
      if (a > 2e-8 && a < 9) {
        const double v = dispersion_objective(y, f.mu, X, th, true);
        stationary = stationary && std::fabs(dispersion_score(y, f.mu, X, th, true)) < 1e-6 && v >= dispersion_objective(y, f.mu, X, th + 0.01, true) &&
                     v >= dispersion_objective(y, f.mu, X, th - 0.01, true);
      }
    }
    const double med = median(est), med_ml = median(est_ml);
    std::printf("dispersion 0.1, 5 + 5 samples, mean 100: median Cox–Reid estimate %.4f, median maximum likelihood %.4f (400 genes)\n", med, med_ml);
    s.check("the median Cox–Reid estimate of 400 genes recovers the dispersion 0.1: relative error", med / alpha - 1, 0.15);
    s.require("the Cox–Reid adjustment raises the estimates, which maximum likelihood biases down", med > med_ml);
    s.require("each estimate is a zero of the derivative and a local maximum of the adjusted likelihood", stationary);
  }

  // ---- the whole analysis: the null, power, threads ------------------------------------------------------------------------------------
  {
    const Simulation null = simulate(4000, 10, 0.0, 51);
    const Design X = Design::two_groups(null.group);
    DifferentialOptions o;
    const Differential d = differential(null.counts, X, o, eight);
    std::vector<double> p;
    std::size_t small = 0;
    for (const auto& g : d.genes)
      if (!std::isnan(g.p_value)) {
        p.push_back(g.p_value);
        small += g.p_value < 0.05;
      }
    const double ks = ks_uniform(p);
    std::printf("null of 4000 genes, 5 + 5 samples: trend %.4f + %.3f / mean (simulated 0.03 + 3 / mean, scattered), sigma_d^2 %.3f, KS distance %.4f, P(p < 0.05) = %.4f\n",
                d.trend_a0, d.trend_a1, d.prior_variance, ks, static_cast<double>(small) / static_cast<double>(p.size()));
    s.check("under the null the Wald p-values are uniform: the Kolmogorov–Smirnov distance (1% critical value 0.026 for 4000)", ks, 0.04);
    // A Wald test with estimated dispersions at 5 + 5 samples is a little liberal (here 0.064 at
    // 0.05); the limit allows for it.
    s.check("the fraction of p-values below 0.05 is 0.05, within 0.02", static_cast<double>(small) / static_cast<double>(p.size()) - 0.05, 0.02);
    s.check("the trend's asymptotic dispersion a₀ near 0.03 (scattered by a log-normal of sd 0.3, mean factor 1.046)", d.trend_a0 / 0.03 - 1, 0.35);
    s.check("its coefficient a₁ near 3", d.trend_a1 / 3 - 1, 0.35);

    const Simulation alt = simulate(4000, 10, 0.1, 61);
    const Differential e = differential(alt.counts, X, o, eight);
    std::size_t changed = 0, found = 0, discoveries = 0, false_discoveries = 0;
    for (std::size_t g = 0; g < e.genes.size(); g++) {
      const bool de = alt.lfc[g] != 0;
      const bool call = !std::isnan(e.genes[g].padj) && e.genes[g].padj < 0.1;
      changed += de;
      found += de && call;
      discoveries += call;
      false_discoveries += call && !de;
    }
    const double power = static_cast<double>(found) / static_cast<double>(changed), fdr = static_cast<double>(false_discoveries) / static_cast<double>(std::max<std::size_t>(1, discoveries));
    std::printf("400 of 4000 genes changed two-fold, 5 + 5 samples, padj < 0.1: %zu discoveries, power %.3f, observed FDR %.3f\n", discoveries, power, fdr);
    s.check("the observed false discovery rate at padj < 0.1 stays near its bound: FDR − 0.1", std::max(0.0, fdr - 0.1), 0.05);
    s.require("two-fold changes are detected with power above one half", power > 0.5, std::to_string(power));
    const Differential e1 = differential(alt.counts, X, o, one);
    bool same = true;
    for (std::size_t g = 0; g < e.genes.size(); g++)
      same = same && (e.genes[g].p_value == e1.genes[g].p_value || (std::isnan(e.genes[g].p_value) && std::isnan(e1.genes[g].p_value))) &&
             e.genes[g].dispersion == e1.genes[g].dispersion;
    s.require("the analysis is the same on 1 and 8 threads, gene for gene", same);
  }
  return s.done();
}
