// Checks of the coalescent without recombination: a genealogy is a function of its seed and
// replicate; under constant size E[T_MRCA] = 2(1 − 1/n), Var[T_MRCA] = Σ_{k=2}^{n} 4/(k(k − 1))²
// and E[L] = 2a_n (Kingman 1982; Hudson 1990, in units of 2N generations); with mutation
// E[S] = θa_n and Var[S] = θa_n + θ²b_n (Watterson 1975), E[π] = θ and Var[π] = b₁θ + b₂θ²
// (Tajima 1983), Var[π − S/a₁] = c₁θ + c₂θ² (Tajima 1989: the variance his D divides by),
// E[ξ_i] = θ/i (Fu 1995), E[θ_H] = θ (Fay and Wu 2000), and Tajima's D near 0; one tree, whose
// sites never show four gametes; under exponential growth and a bottleneck P(T₂ > t) =
// e^{−Λ(t)} (time rescaling), and after an expansion Tajima's D is negative; replicates in
// parallel are the same on 1 and 8 threads. Each expectation is checked by simulation within
// four standard errors (the z-score is reported).

#include <evolution/population/coalescent.hpp>
#include <evolution/population/diversity.hpp>

#include <calx/verification/check.hpp>

#include <chrono>
#include <cmath>
#include <string>
#include <vector>

using namespace evolution::population;

namespace {

struct Moments {
  double n = 0, sum = 0, sum2 = 0, sum3 = 0, sum4 = 0;
  void add(double x) {
    n++;
    sum += x;
    sum2 += x * x;
    sum3 += x * x * x;
    sum4 += x * x * x * x;
  }
  double mean() const { return sum / n; }
  double variance() const { return (sum2 - sum * sum / n) / (n - 1); }
  double se() const { return std::sqrt(variance() / n); }
  // The standard error of the sample variance: √((μ₄ − σ⁴)/n), with the central moments estimated.
  double se_variance() const {
    const double m = mean(), m2 = sum2 / n - m * m;
    const double m4 = sum4 / n - 4 * m * sum3 / n + 6 * m * m * sum2 / n - 3 * m * m * m * m;
    return std::sqrt((m4 - m2 * m2) / n);
  }
  double z(double expected) const { return (mean() - expected) / se(); }
  double z_variance(double expected) const { return (variance() - expected) / se_variance(); }
};

std::string fmt(double x) { return calx::verification::exponential(x); }

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- the genealogy is a function of the seed --------------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 12;
    spec.theta = 4;
    const Genealogy a = coalescent(spec, 7, 3), b = coalescent(spec, 7, 3), c = coalescent(spec, 7, 4);
    const bool same = a.time == b.time && a.edges.size() == b.edges.size() && a.mutations.size() == b.mutations.size();
    s.require("a genealogy is a function of its seed and replicate", same && a.time != c.time);
    s.require("n samples coalesce in n − 1 events into 2n − 1 nodes", a.nodes() == 23 && a.common_ancestors == 11 && a.edges.size() == 22);
    spec.theta = 0;
    s.require("the genealogy does not depend on θ (mutations have their own stream)", coalescent(spec, 7, 3).time == a.time);
  }

  // ---- constant size, no recombination ------------------------------------------------------------------
  const std::size_t n = 10;
  const double a_n = watterson_a(n), b_n = watterson_b(n);
  {
    CoalescentSpec spec;
    spec.samples = n;
    Moments tmrca, length;
    bool one_tree = true;
    for (std::uint32_t rep = 0; rep < 20000; rep++) {
      const Genealogy g = coalescent(spec, 1, rep);
      std::size_t trees = 0;
      for_each_tree(g, [&](const TreeView& v) {
        trees++;
        tmrca.add(v.tmrca());
        length.add(v.length);
      });
      one_tree = one_tree && trees == 1;
    }
    double var_tmrca = 0;
    for (std::size_t k = 2; k <= n; k++) var_tmrca += 4 / std::pow(static_cast<double>(k * (k - 1)), 2);
    s.require("without recombination the genealogy is one tree", one_tree);
    s.check("E[T_MRCA] = 2(1 − 1/n) = 1.8 for n = 10 (z-score, 20000 genealogies)", tmrca.z(expected_tmrca(n)), 4, "SE");
    s.check("Var[T_MRCA] = Σ 4/(k(k − 1))² (z-score)", tmrca.z_variance(var_tmrca), 4, "SE");
    s.check("E[L] = 2a_n, the total branch length (z-score)", length.z(expected_tree_length(n)), 4, "SE");
  }

  // ---- mutation: the estimators of θ and their variances ---------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = n;
    spec.theta = 5;
    const double theta = spec.theta;
    const TajimaConstants tc = tajima_constants(n);
    Moments S, pi, d, th, D;
    std::vector<Moments> xi(n);
    bool segregating = true, four_gametes = false;
    for (std::uint32_t rep = 0; rep < 20000; rep++) {
      const Genealogy g = coalescent(spec, 2, rep);
      const Haplotypes h = haplotypes(g);
      const AlleleCounts c = allele_counts(h);
      const double sites = static_cast<double>(segregating_sites(c));
      segregating = segregating && sites == static_cast<double>(h.sites);
      const double p = nucleotide_diversity(c);
      S.add(sites);
      pi.add(p);
      d.add(p - sites / tc.a1);
      th.add(theta_h(c));
      if (sites > 0) D.add(tajima_d(n, sites, p));
      std::vector<double> spectrum = site_frequency_spectrum(c, n);
      for (std::size_t i = 1; i < n; i++) xi[i].add(spectrum[i]);
      if (rep < 200)
        for (std::size_t x = 0; x < h.sites && !four_gametes; x++)
          for (std::size_t y = x + 1; y < h.sites && !four_gametes; y++) {
            int seen = 0;
            for (std::size_t k = 0; k < n; k++) seen |= 1 << (2 * h(k, x) + h(k, y));
            four_gametes = seen == 15;
          }
    }
    s.require("every infinite-sites mutation makes a segregating site", segregating);
    s.require("on one tree no two sites show all four gametes (the four-gamete test)", !four_gametes);
    s.check("E[S] = θa_n (Watterson 1975; z-score)", S.z(theta * a_n), 4, "SE");
    s.check("Var[S] = θa_n + θ²b_n (Watterson 1975; z-score)", S.z_variance(theta * a_n + theta * theta * b_n), 4, "SE");
    s.check("E[π] = θ (Tajima 1983; z-score)", pi.z(theta), 4, "SE");
    s.check("Var[π] = b₁θ + b₂θ² (Tajima 1983; z-score)", pi.z_variance(tc.b1 * theta + tc.b2 * theta * theta), 4, "SE");
    s.check("Var[π − S/a₁] = c₁θ + c₂θ², the variance of Tajima's numerator (z-score)", d.z_variance(tc.c1 * theta + tc.c2 * theta * theta), 4, "SE");
    s.check("E[θ_H] = θ, so E[H] = 0 (Fay and Wu 2000; z-score)", th.z(theta), 4, "SE");
    double worst = 0;
    std::size_t at = 0;
    for (std::size_t i = 1; i < n; i++) {
      const double z = xi[i].z(theta / static_cast<double>(i));
      if (std::fabs(z) > std::fabs(worst)) {
        worst = z;
        at = i;
      }
    }
    s.check("E[ξ_i] = θ/i for i = 1 … 9 (Fu 1995; the largest z-score, at i = " + std::to_string(at) + ")", worst, 4, "SE");
    s.check("Tajima's D is near 0 under the neutral coalescent (its mean)", D.mean(), 0.1);
    s.require("and its variance is below 1 (D's denominator estimates the variance; D is not standardised exactly)",
              D.variance() < 1 && D.variance() > 0.5, fmt(D.variance()));
  }

  // ---- changing size: time rescaling ---------------------------------------------------------------------
  {
    const auto survival = [&](const Demography& demography, std::uint64_t seed, const std::string& what) {
      CoalescentSpec spec;
      spec.samples = 2;
      spec.demography = demography;
      const std::size_t R = 40000;
      std::vector<double> times(R);
      for (std::uint32_t rep = 0; rep < R; rep++) times[rep] = coalescent(spec, seed, rep).time.back();
      double worst = 0;
      for (const double t : {0.05, 0.15, 0.3, 0.6}) {
        const double p = std::exp(-demography.intensity(t));
        double k = 0;
        for (const double x : times) k += x > t;
        const double z = (k / R - p) / std::sqrt(p * (1 - p) / R);
        if (std::fabs(z) > std::fabs(worst)) worst = z;
      }
      s.check("P(T₂ > t) = exp(−Λ(t)) at t = 0.05, 0.15, 0.3, 0.6 under " + what + " (largest z-score)", worst, 4, "SE");
    };
    survival(Demography::exponential(5), 3, "exponential growth, α = 5");
    survival(Demography::bottleneck(0.1, 0.2, 0.1), 4, "a bottleneck to 0.1 between 0.1 and 0.3");
    const Demography d({{0, 1, 0}, {0.2, 0.5, 2}, {0.5, 0.3, 0}});
    s.check("Λ⁻¹ inverts Λ across epochs of constant and exponential size", std::fabs(d.intensity(d.advance(0.1, 0.9)) - d.intensity(0.1) - 0.9), 1e-12);
    survival(d, 5, "three epochs, one growing");

    CoalescentSpec spec;
    spec.samples = 20;
    spec.theta = 20;
    spec.demography = Demography::exponential(20);
    Moments D;
    for (std::uint32_t rep = 0; rep < 2000; rep++) {
      const AlleleCounts c = allele_counts(haplotypes(coalescent(spec, 6, rep)));
      if (segregating_sites(c) > 0) D.add(tajima_d(c));
    }
    s.require("after an expansion (α = 20) Tajima's D is negative: an excess of rare alleles", D.mean() < -0.8, fmt(D.mean()));
  }

  // ---- replicates in parallel ---------------------------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 30;
    spec.theta = 10;
    spec.rho = 5;
    const auto summary = [](const Genealogy& g) { return std::pair<double, std::size_t>(g.time.back(), g.mutations.size()); };
    calx::scheduler::ThreadPool one(1), eight(8);
    const auto a = coalescent_replicates(spec, 9, 500, one, summary);
    const auto b = coalescent_replicates(spec, 9, 500, eight, summary);
    s.require("replicates in parallel are the same on 1 and 8 threads", a == b);
  }

  // ---- speed --------------------------------------------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 1000;
    spec.theta = 100;
    const auto t0 = std::chrono::steady_clock::now();
    const Genealogy g = coalescent(spec, 10);
    const Haplotypes h = haplotypes(g);
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    s.require("a genealogy of 1000 samples with θ = 100 and its haplotypes in under 200 ms", ms < 200,
              fmt(ms) + " ms, " + std::to_string(h.sites) + " sites");
  }
  return s.done();
}
