// Genotype likelihoods of a diploid from the reads at a site, the posterior under a prior of
// heterozygosity, genotype and site qualities, and Fisher's exact test of strand bias.
//
// - Bases (Li 2011, "A statistical framework for SNP calling, mutation discovery, association
//   mapping and population genetical parameter estimation from sequencing data",
//   Bioinformatics 27, 2987–2993, eq. 2 and §2.2). At a biallelic site with reference allele
//   R and alternate A, a diploid genotype with r copies of R (and 2 − r of A) and reads whose
//   bases are R or A, the base j read with error probability ε_j = 10^(−Q_j/10) (its Phred
//   quality) has
//
//     P(b_j = R | r) = [r (1 − ε_j) + (2 − r) ε_j] / 2,   P(b_j = A | r) = [(2 − r)(1 − ε_j) + r ε_j] / 2,
//
//   and the reads are independent given the genotype: L(r) = ∏_j P(b_j | r). The genotype is
//   named here, as VCF names it, by its count g = 2 − r of alternate alleles (0/0, 0/1, 1/1).
//   *Remark.* The model has ε, not ε/3, for reading the other allele: it is Li's biallelic
//   simplification, every error landing on the one other allele considered. With n alleles
//   the same rule, P(b | {a₁, a₂}) = [p(b | a₁) + p(b | a₂)] / 2 with p(b | a) = 1 − ε when b = a
//   and ε otherwise, gives the likelihoods of the n(n + 1)/2 diploid genotypes, which reduce to
//   eq. 2 for n = 2; bases that are none of the alleles are left out.
// - Indels are counted, not scored by base: of the d reads spanning a site, k carry the indel,
//   and each read is taken to show it (or not) wrongly with a fixed probability ε (0.01 by
//   default, Phred 20: an indel error in reads with an indel rate near 10⁻⁴ a base is far rarer,
//   and the rate stands for mapping and alignment artefacts), in the same formula with k
//   alternate and d − k reference reads.
// - Prior (the infinite-sites prior of a sample of two chromosomes, P(k alternate copies) = θ/k;
//   DePristo et al. 2011, Nature Genetics 43, 491–498, the heterozygosity of GATK's caller):
//   π(0/0) = 1 − 3θ/2, π(0/a) = θ, π(a/a) = θ/2, and for two different alternate alleles
//   π(a/b) = θ² (two mutations at one site), normalised; θ = 10⁻³ for SNVs and 1.25 · 10⁻⁴ for
//   indels (GATK's defaults). The posterior is P(g | D) ∝ L(g) π(g).
// - Qualities: GQ = −10 log₁₀ P(g ≠ ĝ | D) for the genotype ĝ of largest posterior, held to 99;
//   QUAL = −10 log₁₀ P(0/0 | D), the Phred-scaled probability that the site is not variant; PL
//   the likelihoods as −10 log₁₀ L(g), less their minimum, rounded. All in log space, so that
//   nothing underflows at any depth.
// - Strand bias: Fisher's exact test of the 2 × 2 table (reference and alternate reads by
//   strand), two-sided, the sum of the hypergeometric probabilities of every table with the
//   margins of the observed one and no greater probability than it (Fisher 1935, The Design of
//   Experiments, ch. 2; Agresti 1992, Statistical Science 7, 131–153), each from log Γ (calx's
//   log_choose); FS = −10 log₁₀ p, GATK's FisherStrand. O(min(margins)).
//   *Remark.* A hypergeometric distribution (pmf, cdf) in calx::special_functions would serve
//   here and elsewhere (enrichment tests); this file builds it from log_choose.
//
// Complexity: O(n_g d) for d bases and n_g genotypes. Probabilities in double precision; the
// likelihoods of thousands of reads are sums of logarithms, exact to rounding.

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <calx/special_functions/gamma.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace genomics::variant {

// A base of the pileup reduced to what genotyping needs.
struct BaseCall {
  sequence::Code base = 0;
  std::uint8_t quality = 0;  // Phred
};

// The index of genotype j/k (j ≤ k) in VCF's order: k(k + 1)/2 + j.
constexpr std::size_t genotype_index(std::size_t j, std::size_t k) { return j <= k ? k * (k + 1) / 2 + j : j * (j + 1) / 2 + k; }
constexpr std::size_t genotype_count(std::size_t alleles) { return alleles * (alleles + 1) / 2; }

// log₁₀ L of each diploid genotype over the alleles (codes), in VCF's order; bases that are
// none of the alleles are left out.
inline std::vector<double> log10_likelihoods(std::span<const BaseCall> bases, std::span<const sequence::Code> alleles) {
  const std::size_t n = alleles.size();
  if (n < 1) throw std::invalid_argument("genomics::variant::log10_likelihoods: at least one allele");
  std::vector<double> l(genotype_count(n), 0.0);
  for (const auto& b : bases) {
    std::size_t which = n;
    for (std::size_t a = 0; a < n; a++)
      if (alleles[a] == b.base) { which = a; break; }
    if (which == n) continue;
    const double e = std::pow(10.0, -static_cast<double>(b.quality) / 10.0);
    for (std::size_t k = 0; k < n; k++)
      for (std::size_t j = 0; j <= k; j++) {
        const double pj = j == which ? 1 - e : e, pk = k == which ? 1 - e : e;
        l[genotype_index(j, k)] += std::log10((pj + pk) / 2);
      }
  }
  return l;
}

// The biallelic likelihoods of eq. 2 from counts: k alternate and d − k reference observations,
// each wrong with probability ε.
inline std::array<double, 3> log10_likelihoods(std::size_t reference_count, std::size_t alternate_count, double error) {
  if (!(error > 0 && error < 1)) throw std::invalid_argument("genomics::variant::log10_likelihoods: 0 < error < 1");
  std::array<double, 3> l{};
  for (int g = 0; g <= 2; g++) {
    const double r = 2 - g;
    const double pr = (r * (1 - error) + (2 - r) * error) / 2, pa = ((2 - r) * (1 - error) + r * error) / 2;
    l[static_cast<std::size_t>(g)] = static_cast<double>(reference_count) * std::log10(pr) + static_cast<double>(alternate_count) * std::log10(pa);
  }
  return l;
}

// The prior of each genotype over n alleles (allele 0 the reference), in VCF's order, normalised.
inline std::vector<double> genotype_prior(std::size_t alleles, double theta) {
  if (!(theta > 0 && theta < 0.5)) throw std::invalid_argument("genomics::variant::genotype_prior: 0 < θ < 1/2");
  std::vector<double> p(genotype_count(alleles));
  for (std::size_t k = 0; k < alleles; k++)
    for (std::size_t j = 0; j <= k; j++)
      p[genotype_index(j, k)] = (j == 0 && k == 0) ? 1 - 1.5 * theta : j == 0 ? theta : j == k ? theta / 2 : theta * theta;
  double s = 0;
  for (double x : p) s += x;
  for (double& x : p) x /= s;
  return p;
}

struct GenotypeCall {
  std::size_t best = 0;              // the index of the genotype of largest posterior
  std::size_t allele_a = 0, allele_b = 0;  // its alleles, a ≤ b
  std::vector<double> log10_posterior;
  int gq = 0;                        // genotype quality, at most 99
  double qual = 0;                   // −10 log₁₀ P(0/0 | D)
  std::vector<int> pl;               // Phred-scaled likelihoods, the least 0
};

// The posterior and qualities from log₁₀ likelihoods (VCF order) and the prior.
inline GenotypeCall posterior(std::span<const double> log10_l, std::span<const double> prior) {
  if (log10_l.size() != prior.size() || log10_l.empty()) throw std::invalid_argument("genomics::variant::posterior: one prior for each genotype");
  GenotypeCall c;
  const std::size_t n = log10_l.size();
  std::vector<double> lp(n);
  double top = -std::numeric_limits<double>::infinity();
  for (std::size_t g = 0; g < n; g++) {
    lp[g] = log10_l[g] + std::log10(prior[g]);
    top = std::max(top, lp[g]);
  }
  double sum = 0;
  for (std::size_t g = 0; g < n; g++) sum += std::pow(10.0, lp[g] - top);
  const double log_total = top + std::log10(sum);
  c.log10_posterior.resize(n);
  for (std::size_t g = 0; g < n; g++) c.log10_posterior[g] = lp[g] - log_total;
  for (std::size_t g = 1; g < n; g++)
    if (c.log10_posterior[g] > c.log10_posterior[c.best]) c.best = g;
  // P(g ≠ best | D), summed from the others so that it does not cancel.
  double other = 0, other_top = -std::numeric_limits<double>::infinity();
  for (std::size_t g = 0; g < n; g++)
    if (g != c.best) other_top = std::max(other_top, c.log10_posterior[g]);
  if (std::isfinite(other_top)) {
    for (std::size_t g = 0; g < n; g++)
      if (g != c.best) other += std::pow(10.0, c.log10_posterior[g] - other_top);
    const double log_other = other_top + std::log10(other);
    c.gq = static_cast<int>(std::clamp(std::floor(-10 * log_other), 0.0, 99.0));
  } else {
    c.gq = 99;
  }
  c.qual = -10 * c.log10_posterior[0];
  if (c.qual < 0) c.qual = 0;  // −0 and rounding
  double lmax = -std::numeric_limits<double>::infinity();
  for (double x : log10_l) lmax = std::max(lmax, x);
  c.pl.resize(n);
  for (std::size_t g = 0; g < n; g++) c.pl[g] = static_cast<int>(std::lround(-10 * (log10_l[g] - lmax)));
  // The alleles of genotype `best`: the k with k(k + 1)/2 ≤ best, then j.
  std::size_t k = 0;
  while ((k + 1) * (k + 2) / 2 <= c.best) k++;
  c.allele_b = k;
  c.allele_a = c.best - k * (k + 1) / 2;
  return c;
}

// ---- Fisher's exact test -------------------------------------------------------------------------------
// log P of the table [[a, b], [c, d]] under the hypergeometric law of its margins.
inline double log_hypergeometric(std::uint64_t a, std::uint64_t b, std::uint64_t c, std::uint64_t d) {
  using calx::special_functions::log_choose;
  const double r1 = static_cast<double>(a + b), r2 = static_cast<double>(c + d), c1 = static_cast<double>(a + c);
  return log_choose(r1, static_cast<double>(a)) + log_choose(r2, static_cast<double>(c)) - log_choose(r1 + r2, c1);
}

// The two-sided p-value of Fisher's exact test of [[a, b], [c, d]].
inline double fisher_exact(std::uint64_t a, std::uint64_t b, std::uint64_t c, std::uint64_t d) {
  const std::uint64_t r1 = a + b, r2 = c + d, c1 = a + c;
  const std::uint64_t lo = c1 > r2 ? c1 - r2 : 0, hi = std::min(r1, c1);
  const double observed = log_hypergeometric(a, b, c, d);
  // Tables of probability no greater than the observed one, with a relative tolerance as R's fisher.test.
  const double tolerance = 1e-7;
  double p = 0;
  for (std::uint64_t x = lo; x <= hi; x++) {
    const double lx = log_hypergeometric(x, r1 - x, c1 - x, r2 - (c1 - x));
    if (lx <= observed + std::log1p(tolerance)) p += std::exp(lx);
  }
  return std::min(1.0, p);
}

// GATK's FisherStrand: −10 log₁₀ of Fisher's p for reference and alternate reads by strand.
inline double strand_bias(std::uint64_t ref_forward, std::uint64_t ref_reverse, std::uint64_t alt_forward, std::uint64_t alt_reverse) {
  const double p = fisher_exact(ref_forward, ref_reverse, alt_forward, alt_reverse);
  return p > 0 ? std::max(0.0, -10 * std::log10(p)) : 3000.0;
}

}  // namespace genomics::variant
