// Checks of genotyping: Li's biallelic likelihoods (eq. 2) against their closed forms on small
// pileups; the counting form equal to the per-base form at equal qualities; the n-allele form
// reducing to eq. 2; VCF's order of genotypes; the posterior, GQ, QUAL and PL against direct
// computation; the prior; and Fisher's exact test against Fisher's tea-tasting experiment and
// against exact enumeration in integers.

#include <genomics/sequence/alphabet.hpp>
#include <genomics/variant/genotype.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace genomics;
using variant::BaseCall;

namespace {

constexpr std::uint32_t CHECK = 0x47454e4f;

// C(n, k) exactly, for small n.
std::uint64_t choose(std::uint64_t n, std::uint64_t k) {
  if (k > n) return 0;
  std::uint64_t c = 1;
  for (std::uint64_t i = 1; i <= k; i++) c = c * (n - k + i) / i;
  return c;
}

// Fisher's two-sided p by enumeration in integers: Σ over tables with P ≤ P(observed).
double fisher_exact_integers(std::uint64_t a, std::uint64_t b, std::uint64_t c, std::uint64_t d) {
  const std::uint64_t r1 = a + b, r2 = c + d, c1 = a + c, n = r1 + r2;
  const auto weight = [&](std::uint64_t x) { return choose(r1, x) * choose(r2, c1 - x); };
  const std::uint64_t observed = weight(a), total = choose(n, c1);
  std::uint64_t sum = 0;
  for (std::uint64_t x = (c1 > r2 ? c1 - r2 : 0); x <= std::min(r1, c1); x++)
    if (weight(x) <= observed) sum += weight(x);
  return static_cast<double>(sum) / static_cast<double>(total);
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, CHECK, 0, 0);
  using sequence::A;
  using sequence::C;
  using sequence::G;

  // ---- Li's eq. 2 ---------------------------------------------------------------------------------------
  {
    // Three reference bases at Q20 (ε = 0.01) and one alternate at Q30 (ε = 0.001):
    //   L(0/0) = 0.99³ · 0.001,  L(0/1) = (1/2)⁴,  L(1/1) = 0.01³ · 0.999.
    const std::vector<BaseCall> bases{{A, 20}, {A, 20}, {A, 20}, {G, 30}};
    const std::vector<sequence::Code> alleles{A, G};
    const auto l = variant::log10_likelihoods(bases, alleles);
    const double e0 = std::log10(std::pow(0.99, 3) * 0.001), e1 = std::log10(std::pow(0.5, 4)), e2 = std::log10(std::pow(0.01, 3) * 0.999);
    const double err = std::max({std::fabs(l[0] - e0), std::fabs(l[1] - e1), std::fabs(l[2] - e2)});
    s.check("eq. 2 on three reference bases at Q20 and one alternate at Q30, against its closed form", err, 1e-12, "log10");
    // Mixed qualities: one reference base at Q10, two alternates at Q40 and Q15, and a base of a third allele, left out.
    const std::vector<BaseCall> mixed{{A, 10}, {G, 40}, {G, 15}, {C, 30}};
    const auto m = variant::log10_likelihoods(mixed, alleles);
    const double a = 0.1, b = 1e-4, c = std::pow(10.0, -1.5);
    const double f0 = std::log10((1 - a) * b * c), f1 = std::log10(0.125), f2 = std::log10(a * (1 - b) * (1 - c));
    s.check("eq. 2 on mixed qualities, a base of a third allele left out", std::max({std::fabs(m[0] - f0), std::fabs(m[1] - f1), std::fabs(m[2] - f2)}), 1e-12, "log10");
    // The counting form at equal qualities.
    std::vector<BaseCall> equal;
    for (int i = 0; i < 7; i++) equal.push_back({A, 20});
    for (int i = 0; i < 5; i++) equal.push_back({G, 20});
    const auto lb = variant::log10_likelihoods(equal, alleles);
    const auto lc = variant::log10_likelihoods(7, 5, 0.01);
    s.check("the counting form equals eq. 2 at equal qualities (7 reference, 5 alternate, Q20)", std::max({std::fabs(lb[0] - lc[0]), std::fabs(lb[1] - lc[1]), std::fabs(lb[2] - lc[2])}),
            1e-12, "log10");
    // Three alleles: the genotypes over {A, G} are eq. 2's when no base is the third allele.
    const std::vector<sequence::Code> three{A, G, C};
    const auto l3 = variant::log10_likelihoods(bases, three);
    s.check("over three alleles, the genotypes 0/0, 0/1, 1/1 are eq. 2's", std::max({std::fabs(l3[0] - l[0]), std::fabs(l3[1] - l[1]), std::fabs(l3[2] - l[2])}), 1e-12, "log10");
    s.require("VCF's order of genotypes: F(j/k) = k(k + 1)/2 + j (0/0, 0/1, 1/1, 0/2, 1/2, 2/2)",
              variant::genotype_index(0, 0) == 0 && variant::genotype_index(0, 1) == 1 && variant::genotype_index(1, 1) == 2 && variant::genotype_index(0, 2) == 3 &&
                  variant::genotype_index(1, 2) == 4 && variant::genotype_index(2, 2) == 5 && variant::genotype_index(2, 1) == 4 && l3.size() == 6);
    // L(1/2) for one G base at Q30 between alleles G and C: (0.999 + 0.001)/2 = 1/2; and the four bases: 0.01³ · 0.5.
    s.check("L(1/2) = ∏ [p(b|G) + p(b|C)]/2", std::fabs(l3[4] - std::log10(std::pow(0.01, 3) * 0.5)), 1e-12, "log10");
  }

  // ---- prior, posterior and qualities ------------------------------------------------------------------
  {
    const double theta = 1e-3;
    const auto prior = variant::genotype_prior(2, theta);
    s.check("the prior: π(0/0) = 1 − 3θ/2, π(0/1) = θ, π(1/1) = θ/2", std::max({std::fabs(prior[0] - (1 - 1.5 * theta)), std::fabs(prior[1] - theta), std::fabs(prior[2] - theta / 2)}),
            1e-15);
    // Ten alternate reads at Q30, against direct computation in double precision.
    const auto l = variant::log10_likelihoods(0, 10, 1e-3);
    const auto g = variant::posterior(l, prior);
    double z = 0;
    double post[3];
    for (int k = 0; k < 3; k++) z += std::pow(10.0, l[k]) * prior[k];
    for (int k = 0; k < 3; k++) post[k] = std::pow(10.0, l[k]) * prior[k] / z;
    const double qual = -10 * std::log10(post[0]), gq = std::floor(-10 * std::log10(post[0] + post[1]));
    s.require("ten alternate reads call 1/1", g.best == 2 && g.allele_a == 1 && g.allele_b == 1);
    s.check("QUAL = −10 log₁₀ P(0/0 | D) against direct computation", std::fabs(g.qual - qual) / qual, 1e-12);
    s.require("GQ = −10 log₁₀ P(g ≠ 1/1 | D), rounded down, against direct computation", g.gq == static_cast<int>(std::min(99.0, gq)), std::to_string(g.gq) + " vs " + std::to_string(gq));
    double perr = 0;
    for (int k = 0; k < 3; k++) perr = std::max(perr, std::fabs(std::pow(10.0, g.log10_posterior[k]) - post[k]));
    s.check("the posterior against direct computation", perr, 1e-15);
    s.require("PL is −10 log₁₀ L less its minimum, rounded", g.pl[2] == 0 && g.pl[1] == static_cast<int>(std::lround(-10 * (l[1] - l[2]))) && g.pl[0] == static_cast<int>(std::lround(-10 * (l[0] - l[2]))));
    // 1000 reads, half and half: the likelihoods underflow in linear space, not here.
    const auto deep = variant::posterior(variant::log10_likelihoods(500, 500, 1e-3), prior);
    s.require("1000 reads, half alternate, call 0/1 with GQ 99 and a finite QUAL", deep.best == 1 && deep.gq == 99 && std::isfinite(deep.qual) && deep.qual > 1000, std::to_string(deep.qual));
    const auto ref = variant::posterior(variant::log10_likelihoods(30, 0, 1e-3), prior);
    s.require("thirty reference reads call 0/0 with QUAL near 0", ref.best == 0 && ref.qual < 1e-3 && ref.gq >= 60, std::to_string(ref.gq));
    const auto prior3 = variant::genotype_prior(3, theta);
    double sum3 = 0;
    for (double x : prior3) sum3 += x;
    s.check("the prior over three alleles sums to 1, π(1/2) = θ² of the unnormalised", std::fabs(sum3 - 1) + std::fabs(prior3[4] / prior3[1] - theta), 1e-15);
  }

  // ---- Fisher's exact test --------------------------------------------------------------------------------------
  {
    // The lady tasting tea (Fisher 1935, The Design of Experiments, ch. 2): eight cups, four of
    // each kind, three of each named rightly: P = 17/70 one-sided, 34/70 two-sided.
    s.check("Fisher's tea-tasting table [[3, 1], [1, 3]]: p = 34/70", std::fabs(variant::fisher_exact(3, 1, 1, 3) - 34.0 / 70.0), 1e-12);
    s.check("the table of a perfect taster [[4, 0], [0, 4]]: p = 2/70", std::fabs(variant::fisher_exact(4, 0, 0, 4) - 2.0 / 70.0), 1e-12);
    double worst = 0;
    for (int trial = 0; trial < 400; trial++) {
      const auto a = r.below(15), b = r.below(15), c = r.below(15), d = r.below(15);
      worst = std::max(worst, std::fabs(variant::fisher_exact(a, b, c, d) - fisher_exact_integers(a, b, c, d)));
    }
    s.check("Fisher's p against exact enumeration in integers (400 tables with entries below 15)", worst, 1e-10);
    const double fs = variant::strand_bias(15, 15, 12, 0), fs_exact = -10 * std::log10(fisher_exact_integers(15, 15, 12, 0));
    s.require("FS is 0 for a balanced table, −10 log₁₀ p of the exact test for one-sided support, finite when p underflows",
              variant::strand_bias(10, 10, 10, 10) < 1e-9 && std::fabs(fs - fs_exact) < 1e-9 && fs > 20 && std::isfinite(variant::strand_bias(1000, 1000, 1000, 0)),
              "FS(15, 15, 12, 0) = " + std::to_string(fs));
  }

  return s.done();
}
