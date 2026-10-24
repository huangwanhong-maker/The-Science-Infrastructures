// Checks of the Hardy–Weinberg tests: the exact law of the heterozygote count against
// Wigginton, Cutler and Abecasis (2005, AJHG 76, 887–893), Table 1 (N = 100 individuals and 21
// copies of the minor allele: each outcome's probability and the P values P_HWE, P_high and
// P_low, printed to six decimals), against brute-force enumeration of every arrangement of the
// alleles for small samples, and against the closed form (their eq. 1) for a million
// individuals; the paper's statements about Table 1 (the exact test rejects at α = 0.05 when 13
// or fewer heterozygotes are seen, an event of probability 0.010; χ² when 15 or fewer, of
// probability 0.070); the exact test never exceeds its nominal size; χ² is N F̂²; and on
// simulated samples the exact test keeps its size in one population and detects the Wahlund
// effect in two pooled.

#include <evolution/population/hardy_weinberg.hpp>
#include <evolution/population/structure.hpp>

#include <calx/special_functions/gamma.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace evolution::population;

namespace {

std::string fmt(double x) { return calx::verification::exponential(x); }

// N individuals, n_minor copies of the minor allele, n_het heterozygotes; the minor allele is
// the alternate here.
GenotypeCounts counts(std::uint64_t N, std::uint64_t n_minor, std::uint64_t het) {
  GenotypeCounts c;
  c.het = het;
  c.hom_alt = (n_minor - het) / 2;
  c.hom_ref = N - het - c.hom_alt;
  return c;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- Table 1 of Wigginton et al. 2005 -----------------------------------------------------------------
  {
    // n_AB, P(n_AB | N = 100, n_A = 21), P_HWE, P_high, P_low; −1 stands for "<.000001".
    // Source: Wigginton, Cutler and Abecasis (2005), AJHG 76, 887–893, Table 1 (PMC1199378,
    // retrieved 2026-10-24).
    const double table[9][5] = {
        {5, -1, -1, 1.000000, -1},
        {7, .000001, .000001, 1.000000, .000001},
        {9, .000047, .000048, .999999, .000048},
        {11, .000870, .000919, .999952, .000919},
        {13, .009375, .010293, .999081, .010293},
        {15, .059283, .069576, .989707, .069576},
        {17, .214465, .284042, .930424, .284042},
        {19, .406355, 1.000000, .715958, .690396},
        {21, .309604, .593645, .309604, 1.000000},
    };
    double worst = 0;
    bool small = true;
    for (const auto& row : table) {
      const ExactTest t = hwe_exact(counts(100, 21, static_cast<std::uint64_t>(row[0])));
      const double got[4] = {t.probability, t.p_hwe, t.p_high, t.p_low};
      for (int k = 0; k < 4; k++) {
        if (row[k + 1] < 0) small = small && got[k] < 0.0000005;
        else worst = std::max(worst, std::fabs(got[k] - row[k + 1]));
      }
    }
    s.check("Table 1 of Wigginton et al. (2005): every probability and P value to six decimals (largest difference)", worst, 5e-7);
    s.require("and those printed as <.000001 round to 0 at six decimals", small);
    const HeterozygoteLaw law = heterozygote_law(100, 21);
    s.require("the heterozygote counts possible with 21 minor copies are 1, 3, …, 21", law.first == 1 && law.last() == 21 && law.probability.size() == 11);
    // The paper: at α = 0.05 the exact test rejects for ≤ 13 heterozygotes, an event of
    // probability 0.010; the χ² test for ≤ 15 (P ≤ .045), of probability 0.070.
    std::uint64_t exact_cut = 0, chi_cut = 0;
    for (std::uint64_t h = 1; h <= 21; h += 2) {
      if (hwe_exact(counts(100, 21, h)).p_hwe <= 0.05) exact_cut = h;
      if (hwe_chi_square(counts(100, 21, h)).p <= 0.05) chi_cut = h;
    }
    double exact_size = 0, chi_size = 0;
    for (std::uint64_t h = 1; h <= 21; h += 2) {
      if (h <= exact_cut) exact_size += law.at(h);
      if (h <= chi_cut) chi_size += law.at(h);
    }
    s.require("at α = 0.05 the exact test rejects for at most 13 heterozygotes, the χ² test for at most 15 (as the paper says)", exact_cut == 13 && chi_cut == 15);
    s.check("the exact test's true size there is 0.010 (the paper's figure; absolute error)", exact_size - 0.010, 0.0005);
    s.check("the χ² test's is 0.070, above its nominal 0.05 (absolute error)", chi_size - 0.070, 0.0005);
  }

  // ---- brute force: every arrangement of the alleles ----------------------------------------------------------
  {
    double worst = 0;
    for (const auto& [N, nA] : std::vector<std::pair<int, int>>{{6, 5}, {8, 5}, {9, 8}, {10, 10}}) {
      // Under HWE given the allele counts, the 2N alleles are a uniformly random arrangement:
      // choose the positions of the n_A minor alleles among 2N, each choice equally likely;
      // individual i holds positions 2i and 2i + 1.
      std::vector<double> count(static_cast<std::size_t>(nA) + 1, 0);
      double total = 0;
      for (std::uint32_t mask = 0; mask < (1u << (2 * N)); mask++) {
        if (__builtin_popcount(mask) != nA) continue;
        int het = 0;
        for (int i = 0; i < N; i++) het += ((mask >> (2 * i)) & 1) != ((mask >> (2 * i + 1)) & 1);
        count[static_cast<std::size_t>(het)]++;
        total++;
      }
      const HeterozygoteLaw law = heterozygote_law(static_cast<std::uint64_t>(N), static_cast<std::uint64_t>(nA));
      for (std::size_t h = 0; h < count.size(); h++) worst = std::max(worst, std::fabs(count[h] / total - law.at(h)));
    }
    s.check("the law of heterozygotes equals the enumeration of every arrangement of the alleles (N ≤ 10; largest difference)", worst, 1e-14);
  }

  // ---- a million individuals: the closed form ------------------------------------------------------------
  {
    const std::uint64_t N = 1000000, nA = 300000, nB = 2 * N - nA;
    const HeterozygoteLaw law = heterozygote_law(N, nA);
    double sum = 0, best = 0;
    std::uint64_t mode = 0;
    for (std::size_t k = 0; k < law.probability.size(); k++) {
      sum += law.probability[k];
      if (law.probability[k] > best) {
        best = law.probability[k];
        mode = law.first + 2 * k;
      }
    }
    // Eq. 1: log P = n_AB log 2 + log N! − log n_AA! − log n_AB! − log n_BB! + log n_A! + log n_B! − log (2N)!
    using calx::special_functions::log_gamma;
    const auto lf = [](double x) { return log_gamma(x + 1); };
    const double het = static_cast<double>(mode), aa = (nA - het) / 2, bb = (nB - het) / 2;
    const double closed = std::exp(het * std::log(2.0) + lf(N) - lf(aa) - lf(het) - lf(bb) + lf(nA) + lf(nB) - lf(2.0 * N));
    s.check("for a million individuals the probabilities sum to 1", sum - 1, 1e-12);
    s.check("and the mode's probability is eq. 1's (relative error)", (best - closed) / closed, 1e-6);
    s.require("the mode is near n_A n_B / (2N − 1), the expected count", std::fabs(het - 300000.0 * 1700000.0 / 1999999.0) < 2, std::to_string(mode));
  }

  // ---- the exact test never exceeds its size ---------------------------------------------------------------
  {
    double largest = 0;
    for (std::uint64_t nA = 1; nA <= 100; nA++) {
      const HeterozygoteLaw law = heterozygote_law(100, nA);
      double size = 0;
      for (std::size_t k = 0; k < law.probability.size(); k++) {
        const std::uint64_t h = law.first + 2 * k;
        GenotypeCounts c;
        c.het = h;
        c.hom_alt = (nA - h) / 2;
        c.hom_ref = 100 - h - c.hom_alt;
        if (hwe_exact(c).p_hwe <= 0.05) size += law.probability[k];
      }
      largest = std::max(largest, size);
    }
    s.require("the exact test's size never exceeds α = 0.05, for N = 100 and every minor count", largest <= 0.05, fmt(largest));
  }

  // ---- χ² -----------------------------------------------------------------------------------------------
  {
    const GenotypeCounts c{298, 489, 213};
    const ChiSquare x = hwe_chi_square(c);
    s.check("χ² equals N F̂², F̂ = 1 − H_obs/H_exp (relative error)", (x.statistic - 1000 * x.f * x.f) / x.statistic, 1e-12);
    s.check("and its P value is the χ²₁ tail, 2(1 − Φ(√χ²)) (absolute error)", x.p - std::erfc(std::sqrt(x.statistic / 2)), 1e-12);
  }

  // ---- simulated samples -------------------------------------------------------------------------------
  {
    StructureSpec spec;
    spec.populations = 2;
    spec.sites = 4000;
    spec.fst = {0.1};
    spec.unadmixed = {200, 200};
    const StructuredSample x = simulate_structure(spec, 21);
    std::vector<std::size_t> first(200);
    for (std::size_t i = 0; i < 200; i++) first[i] = i;
    double one = 0, pooled = 0, mean_f = 0;
    for (std::size_t j = 0; j < spec.sites; j++) {
      one += hwe_exact(genotype_counts(x.genotypes.site(j), first)).p_hwe <= 0.05;
      const GenotypeCounts c = genotype_counts(x.genotypes.site(j));
      pooled += hwe_exact(c).p_low <= 0.05;
      mean_f += hwe_chi_square(c).f;
    }
    s.require("in one population at equilibrium the exact test rejects at most 5% of sites (α = 0.05)", one / 4000 <= 0.05, fmt(one / 4000));
    // Two populations pooled: a deficit of heterozygotes, F̂ ≈ F_ST between halves ≈ F/2 (Wahlund 1928).
    s.require("two populations pooled show the Wahlund effect: F̂ > 0 on average and heterozygote deficits often significant",
              mean_f / 4000 > 0.03 && pooled / 4000 > 0.1, "mean F̂ " + fmt(mean_f / 4000) + ", " + fmt(pooled / 4000) + " of sites with P_low ≤ 0.05");
  }
  return s.done();
}
