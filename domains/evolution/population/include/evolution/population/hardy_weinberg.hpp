// Tests of Hardy–Weinberg equilibrium at a biallelic site: Pearson's χ² test and the exact
// test of Wigginton, Cutler and Abecasis.
//
// Under random union of gametes a sample of N diploids with allele frequency p has genotype
// frequencies p², 2pq, q² (Hardy 1908; Weinberg 1908).
//
// - χ²: the goodness of fit of the counts (n_AA, n_AB, n_BB) to N p̂², 2N p̂q̂, N q̂² with p̂ the
//   sample frequency, one degree of freedom. It equals N F̂² with F̂ = 1 − H_obs/H_exp the
//   inbreeding coefficient (the 2 × 2 table of the two alleles of an individual), and its tail
//   is calx's χ² tail. It is asymptotic, and anticonservative when the minor allele is rare.
//
// - The exact test conditions on the allele counts: given N and the n_A copies of the rarer
//   allele (n_B = 2N − n_A), the number of heterozygotes has the law (Levene 1949; Haldane 1954;
//   Wigginton et al. 2005, AJHG 76, 887–893, eq. 1)
//
//     P(n_AB | N, n_A) = 2^{n_AB} N! / (n_AA! n_AB! n_BB!) · n_A! n_B! / (2N)!,
//
//   with n_AB ≡ n_A (mod 2), n_AA = (n_A − n_AB)/2, n_BB = (n_B − n_AB)/2. The probabilities are
//   computed by the recurrences of their eq. 2 (Guo and Thompson 1992),
//
//     P(n_AB − 2) = P(n_AB) · n_AB (n_AB − 1) / (4 (n_AA + 1)(n_BB + 1)),
//     P(n_AB + 2) = P(n_AB) · 4 n_AA n_BB / ((n_AB + 2)(n_AB + 1)),
//
//   started at 1 near the mode, n_AB ≈ n_A n_B / (2N), and normalised; so nothing overflows or
//   underflows near the mass, and terms far in the tails underflow harmlessly to 0. The p-values
//   are their P_HWE = Σ P(k) over the k with P(k) ≤ P(n_AB) (the probability of a sample
//   no more likely than the one observed, as Fisher's exact test), and the one-sided
//   P_low = P(N_AB ≤ n_AB), P_high = P(N_AB ≥ n_AB). Probabilities that differ by a relative
//   1e-9 count as equal, so that ties of the exact law are not split by rounding.
//
// Complexity: O(n_A) for the exact test, O(1) for χ². Accuracy: the recurrence loses about
// one ulp per step, so the probabilities near the mode carry a relative error below 1e-12 for
// samples of a million; the p-values are sums of them, never 1 − a sum.

#pragma once

#include <evolution/population/genotypes.hpp>

#include <calx/special_functions/gamma.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace evolution::population {

// The genotypes of a site: homozygous for the reference allele, heterozygous, homozygous for
// the alternate.
struct GenotypeCounts {
  std::uint64_t hom_ref = 0, het = 0, hom_alt = 0;
  std::uint64_t individuals() const { return hom_ref + het + hom_alt; }
  std::uint64_t alt_copies() const { return 2 * hom_alt + het; }
};

inline GenotypeCounts genotype_counts(std::span<const std::int8_t> dosages, std::span<const std::size_t> which = {}) {
  GenotypeCounts c;
  const auto add = [&](std::int8_t x) {
    if (x == 0) c.hom_ref++;
    else if (x == 1) c.het++;
    else if (x == 2) c.hom_alt++;
  };
  if (which.empty())
    for (const std::int8_t x : dosages) add(x);
  else
    for (const std::size_t i : which) add(dosages[i]);
  return c;
}

struct ChiSquare {
  double statistic = 0, p = 1;
  double f = 0;  // 1 − H_obs / H_exp
};

inline ChiSquare hwe_chi_square(const GenotypeCounts& c) {
  ChiSquare r;
  const double n = static_cast<double>(c.individuals());
  if (n == 0) return r;
  const double p = static_cast<double>(c.alt_copies()) / (2 * n), q = 1 - p;
  if (p == 0 || q == 0) return r;
  const double expected[3] = {n * q * q, 2 * n * p * q, n * p * p};
  const double observed[3] = {static_cast<double>(c.hom_ref), static_cast<double>(c.het), static_cast<double>(c.hom_alt)};
  for (int k = 0; k < 3; k++) r.statistic += (observed[k] - expected[k]) * (observed[k] - expected[k]) / expected[k];
  r.f = 1 - observed[1] / expected[1];
  r.p = calx::special_functions::chi_square_q(r.statistic, 1);
  return r;
}

// The law of the number of heterozygotes among N diploids carrying n_A copies of one allele:
// probability[k] = P(n_AB = first + 2k), summing to 1.
struct HeterozygoteLaw {
  std::uint64_t first = 0;
  std::vector<double> probability;
  std::uint64_t last() const { return first + 2 * (probability.size() - 1); }
  double at(std::uint64_t het) const {
    if (het < first || het > last() || (het - first) % 2) return 0;
    return probability[(het - first) / 2];
  }
};

inline HeterozygoteLaw heterozygote_law(std::uint64_t individuals, std::uint64_t copies) {
  const std::uint64_t two_n = 2 * individuals;
  if (copies > two_n) throw std::invalid_argument("evolution::population::heterozygote_law: at most 2N copies");
  const std::uint64_t n_a = std::min(copies, two_n - copies), n_b = two_n - n_a;
  HeterozygoteLaw law;
  law.first = n_a % 2;
  const std::size_t states = static_cast<std::size_t>((n_a - law.first) / 2 + 1);
  law.probability.assign(states, 0.0);
  if (individuals == 0) {
    law.probability[0] = 1;
    return law;
  }
  // Start near the mode, with the parity of n_A.
  std::uint64_t mid = static_cast<std::uint64_t>(static_cast<double>(n_a) * static_cast<double>(n_b) / static_cast<double>(two_n));
  if ((mid & 1) != (n_a & 1)) mid++;
  mid = std::min(mid, n_a);
  const std::size_t m = static_cast<std::size_t>((mid - law.first) / 2);
  std::vector<double>& P = law.probability;
  P[m] = 1;
  double sum = 1;
  // Downwards: from het h to h − 2, with hom_a = (n_a − h)/2, hom_b = (n_b − h)/2 at h.
  for (std::size_t k = m; k > 0; k--) {
    const double h = static_cast<double>(law.first + 2 * k);
    const double hom_a = (static_cast<double>(n_a) - h) / 2, hom_b = (static_cast<double>(n_b) - h) / 2;
    P[k - 1] = P[k] * h * (h - 1) / (4 * (hom_a + 1) * (hom_b + 1));
    sum += P[k - 1];
  }
  // Upwards: from h to h + 2.
  for (std::size_t k = m; k + 1 < states; k++) {
    const double h = static_cast<double>(law.first + 2 * k);
    const double hom_a = (static_cast<double>(n_a) - h) / 2, hom_b = (static_cast<double>(n_b) - h) / 2;
    P[k + 1] = P[k] * 4 * hom_a * hom_b / ((h + 2) * (h + 1));
    sum += P[k + 1];
  }
  for (double& x : P) x /= sum;
  return law;
}

struct ExactTest {
  double p_hwe = 1;        // Σ P(k) over P(k) ≤ P(observed)
  double p_low = 1;        // P(N_AB ≤ observed): a deficit of heterozygotes
  double p_high = 1;       // P(N_AB ≥ observed): an excess
  double probability = 1;  // P(observed)
};

// The exact test of one site's genotype counts.
inline ExactTest hwe_exact(const GenotypeCounts& c) {
  ExactTest r;
  const std::uint64_t n = c.individuals();
  if (n == 0) return r;
  const HeterozygoteLaw law = heterozygote_law(n, c.alt_copies());
  const double observed = law.at(c.het);
  r.probability = observed;
  const double tie = observed * (1 + 1e-9);
  double p_hwe = 0, low = 0, high = 0;
  for (std::size_t k = 0; k < law.probability.size(); k++) {
    const std::uint64_t h = law.first + 2 * k;
    const double x = law.probability[k];
    if (x <= tie) p_hwe += x;
    if (h <= c.het) low += x;
    if (h >= c.het) high += x;
  }
  r.p_hwe = std::min(1.0, p_hwe);
  r.p_low = std::min(1.0, low);
  r.p_high = std::min(1.0, high);
  return r;
}

// The exact test at every site of a genotype matrix (over the samples listed, or all).
inline std::vector<ExactTest> hwe_exact(const Genotypes& g, std::span<const std::size_t> which = {}) {
  std::vector<ExactTest> r(g.sites);
  for (std::size_t j = 0; j < g.sites; j++) r[j] = hwe_exact(genotype_counts(g.site(j), which));
  return r;
}

}  // namespace evolution::population
