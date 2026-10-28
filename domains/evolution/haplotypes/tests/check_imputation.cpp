// Checks of imputation by the copying model: on the coalescent's samples, the alleles of held-out
// haplotypes and the genotypes of held-out diploids at the sites they were not typed at,
// imputed from a panel and compared with the truth and with the panel's allele frequency
// (concordance of the best guess, the squared error of the dosage, r² at common sites); and the
// posterior predictive's exactness under the model, by the identity E[x_j] = E[d_j] at hidden
// sites of targets drawn from the model (d_j the imputed dosage, the expectation of x_j given
// the typed alleles).

#include <evolution/haplotypes/imputation.hpp>
#include <evolution/population/coalescent.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace evolution::haplotypes;
using evolution::population::CoalescentSpec;
using evolution::population::coalescent;

namespace {

std::string fmt(double x) { return calx::verification::exponential(x); }

struct Split {
  Haplotypes panel;
  std::vector<std::vector<std::int8_t>> targets;  // held-out haplotypes
  std::vector<double> frequency;                  // of allele 1 in the panel
};
Split split(const Haplotypes& all, std::size_t K) {
  Split s;
  s.panel.sequences = K;
  s.panel.sites = all.sites;
  s.panel.positions = all.positions;
  s.targets.assign(all.sequences - K, std::vector<std::int8_t>(all.sites));
  s.frequency.resize(all.sites);
  for (std::size_t j = 0; j < all.sites; j++) {
    const auto row = all.site(j);
    s.panel.alleles.insert(s.panel.alleles.end(), row.begin(), row.begin() + static_cast<std::ptrdiff_t>(K));
    s.frequency[j] = static_cast<double>(std::count(row.begin(), row.begin() + static_cast<std::ptrdiff_t>(K), 1)) / static_cast<double>(K);
    for (std::size_t i = K; i < all.sequences; i++) s.targets[i - K][j] = row[i];
  }
  return s;
}

struct Tally {
  double calls = 0, agree = 0, baseline = 0, sq = 0, sq_baseline = 0;
  double common_calls = 0, common_agree = 0, common_baseline = 0;
  double r2 = 0, r2_sites = 0;
  double concordance() const { return agree / calls; }
};

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- haplotypes ------------------------------------------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 220;
    spec.theta = 100;
    spec.rho = 100;
    Tally t;
    bool typed_kept = true;
    for (std::uint32_t rep = 0; rep < 5; rep++) {
      const Split sp = split(evolution::population::haplotypes(coalescent(spec, 31, rep)), 200);
      const Chain c = copying_chain(sp.panel, {100, -1, -1});
      const std::size_t N = sp.panel.sites, T = sp.targets.size();
      // A quarter of the sites typed, the same for every target, as on an array.
      calx::random::Stream r(31, 1, rep, 0);
      std::vector<std::uint8_t> typed(N);
      for (auto& x : typed) x = r.bernoulli(0.25);
      std::vector<std::vector<double>> dosage(T);
      for (std::size_t i = 0; i < T; i++) {
        std::vector<std::int8_t> x = sp.targets[i];
        for (std::size_t j = 0; j < N; j++)
          if (!typed[j]) x[j] = -1;
        dosage[i] = impute_haplotype(sp.panel, c, x).dosage;
        for (std::size_t j = 0; j < N; j++) typed_kept = typed_kept && (!typed[j] || dosage[i][j] == sp.targets[i][j]);
      }
      for (std::size_t j = 0; j < N; j++) {
        if (typed[j]) continue;
        const double p = sp.frequency[j];
        const int major = p > 0.5;
        const bool common = std::min(p, 1 - p) >= 0.05;
        std::vector<double> d(T), truth(T);
        for (std::size_t i = 0; i < T; i++) {
          const int x = sp.targets[i][j];
          const int guess = dosage[i][j] >= 0.5;
          t.calls++;
          t.agree += guess == x;
          t.baseline += major == x;
          t.sq += (dosage[i][j] - x) * (dosage[i][j] - x);
          t.sq_baseline += (p - x) * (p - x);
          if (common) {
            t.common_calls++;
            t.common_agree += guess == x;
            t.common_baseline += major == x;
          }
          d[i] = dosage[i][j];
          truth[i] = x;
        }
        const double r2 = dosage_r2(d, truth);
        if (common && std::isfinite(r2)) {
          t.r2 += r2;
          t.r2_sites++;
        }
      }
    }
    s.require("a typed site keeps its allele as its dosage", typed_kept);
    s.require("haplotypes: untyped sites imputed with a concordance above the panel's major allele's (5 × 20 haplotypes, panels of 200, a quarter of the sites typed)",
              t.concordance() > t.baseline / t.calls, fmt(t.concordance()) + " against " + fmt(t.baseline / t.calls) + " over " + fmt(t.calls) + " calls");
    s.require("at common sites (panel MAF ≥ 0.05) the concordance is well above the major allele's", t.common_agree / t.common_calls > t.common_baseline / t.common_calls + 0.05,
              fmt(t.common_agree / t.common_calls) + " against " + fmt(t.common_baseline / t.common_calls));
    s.require("the imputed dosage's squared error is below the panel frequency's (as a dosage)", t.sq < 0.5 * t.sq_baseline,
              fmt(t.sq / t.calls) + " against " + fmt(t.sq_baseline / t.calls) + " per call");
    s.require("the mean r² of dosage and truth at common untyped sites is above 0.7", t.r2 / t.r2_sites > 0.7, fmt(t.r2 / t.r2_sites) + " over " + fmt(t.r2_sites) + " sites");
  }

  // ---- diploids ------------------------------------------------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 120;
    spec.theta = 50;
    spec.rho = 50;
    double calls = 0, agree = 0, baseline = 0, sq = 0, sq_baseline = 0;
    bool laws = true;
    for (std::uint32_t rep = 0; rep < 4; rep++) {
      const Split sp = split(evolution::population::haplotypes(coalescent(spec, 32, rep)), 100);
      const Chain c = copying_chain(sp.panel, {50, -1, -1});
      const std::size_t N = sp.panel.sites;
      calx::random::Stream r(32, 1, rep, 0);
      std::vector<std::uint8_t> typed(N);
      for (auto& x : typed) x = r.bernoulli(0.25);
      for (std::size_t i = 0; i + 1 < sp.targets.size(); i += 2) {
        std::vector<std::int8_t> g(N);
        for (std::size_t j = 0; j < N; j++) g[j] = typed[j] ? static_cast<std::int8_t>(sp.targets[i][j] + sp.targets[i + 1][j]) : -1;
        const Imputation im = impute_genotypes(sp.panel, c, g);
        for (std::size_t j = 0; j < N; j++) {
          laws = laws && std::fabs(im.at(j, 0) + im.at(j, 1) + im.at(j, 2) - 1) < 1e-12 && std::fabs(im.dosage[j] - im.at(j, 1) - 2 * im.at(j, 2)) < 1e-12;
          if (typed[j]) continue;
          const int truth = sp.targets[i][j] + sp.targets[i + 1][j];
          const int guess = static_cast<int>(std::max_element(&im.genotype[3 * j], &im.genotype[3 * j] + 3) - &im.genotype[3 * j]);
          // The baseline: the most probable genotype under Hardy–Weinberg at the panel's frequency.
          const double p = sp.frequency[j];
          const double hw[3] = {(1 - p) * (1 - p), 2 * p * (1 - p), p * p};
          const int mode = static_cast<int>(std::max_element(hw, hw + 3) - hw);
          calls++;
          agree += guess == truth;
          baseline += mode == truth;
          sq += (im.dosage[j] - truth) * (im.dosage[j] - truth);
          sq_baseline += (2 * p - truth) * (2 * p - truth);
        }
      }
    }
    s.require("diploids: each site's genotype law sums to one, and the dosage is its mean", laws);
    s.require("diploids: untyped genotypes imputed from unphased ones with a concordance above Hardy–Weinberg's mode (4 × 10 diploids, panels of 100)", agree > baseline,
              fmt(agree / calls) + " against " + fmt(baseline / calls) + " over " + fmt(calls) + " calls");
    s.require("and the dosage's squared error below that of twice the panel frequency", sq < 0.5 * sq_baseline, fmt(sq / calls) + " against " + fmt(sq_baseline / calls) + " per call");
  }

  // ---- exactness under the model -----------------------------------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 40;
    spec.theta = 30;
    spec.rho = 30;
    const Haplotypes panel = evolution::population::haplotypes(coalescent(spec, 33));
    CopyingModel cm;
    cm.rho = 30;
    cm.epsilon = 0.02;
    const Chain c = copying_chain(panel, cm);
    double n = 0, sum = 0, sum2 = 0, gn = 0, gsum = 0, gsum2 = 0;
    for (std::uint32_t rep = 0; rep < 300; rep++) {
      calx::random::Stream r(33, MOSAIC, rep, 0);
      const Mosaic a = sample_mosaic(panel, c, r), b = sample_mosaic(panel, c, r);
      std::vector<std::int8_t> x = a.alleles, g(a.alleles.size());
      std::vector<std::uint8_t> hidden(x.size());
      for (std::size_t j = 0; j < x.size(); j++) {
        hidden[j] = r.bernoulli(0.5);
        if (hidden[j]) x[j] = -1;
        g[j] = hidden[j] ? -1 : static_cast<std::int8_t>(a.alleles[j] + b.alleles[j]);
      }
      const Imputation hx = impute_haplotype(panel, c, x);
      double d = 0, k = 0;
      for (std::size_t j = 0; j < x.size(); j++)
        if (hidden[j]) {
          d += a.alleles[j] - hx.dosage[j];
          k++;
        }
      d /= k;
      n++;
      sum += d;
      sum2 += d * d;
      if (rep < 60) {
        const Imputation hg = impute_genotypes(panel, c, g);
        double e = 0;
        for (std::size_t j = 0; j < g.size(); j++)
          if (hidden[j]) e += a.alleles[j] + b.alleles[j] - hg.dosage[j];
        e /= k;
        gn++;
        gsum += e;
        gsum2 += e * e;
      }
    }
    const double mean = sum / n, se = std::sqrt((sum2 - sum * sum / n) / (n - 1) / n);
    const double gmean = gsum / gn, gse = std::sqrt((gsum2 - gsum * gsum / gn) / (gn - 1) / gn);
    s.check("E[x_j − d_j] = 0 at hidden sites of 300 haplotypes drawn from the model (z-score)", mean / se, 4, "SE");
    s.check("E[g_j − d_j] = 0 at hidden sites of 60 diploids drawn from the model (z-score)", gmean / gse, 4, "SE");
  }
  return s.done();
}
