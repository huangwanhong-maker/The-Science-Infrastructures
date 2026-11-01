// Checks of the genotypes of structural variants: SVTyper's bayes_gt reproduced (its own code run
// on 26 pairs of supports, tests/fixtures/structural-svtyper.tsv, made by make_structural_fixtures.py):
// the three log₁₀ likelihoods to 10⁻⁹, the genotype, GQ and SQ; the expected fractions of a tandem
// duplication of k extra copies, k/(k + 2) heterozygous and k/(k + 1) homozygous, against the
// copies of the junction in a constructed pileup; real-valued supports give the same as integers
// where they are integers.

#include <genomics/structural/genotype.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace genomics::structural;

int main() {
  calx::verification::Suite s;
  {
    std::ifstream f(std::string(GENOMICS_FIXTURES) + "/structural-svtyper.tsv");
    std::string line;
    std::size_t rows = 0, likelihoods = 0, genotypes = 0, gqs = 0, sqs = 0;
    double worst = 0;
    while (std::getline(f, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream in(line);
      double ref, alt, l0, l1, l2, sq;
      int dup, gt, gq;
      in >> ref >> alt >> dup >> l0 >> l1 >> l2 >> gt >> gq >> sq;
      // SVTyper's fractions: (10⁻³, 1/2, 0.9), and for a duplication (10⁻², 0.2, 1/3).
      const SvGenotype g = dup ? genotype_of(ref, alt, 0.2, 1 / 3.0, 1e-2) : genotype_of(ref, alt, 0.5, 0.9, 1e-3);
      rows++;
      const double e = std::max({std::fabs(g.log10_likelihood[0] - l0), std::fabs(g.log10_likelihood[1] - l1), std::fabs(g.log10_likelihood[2] - l2)});
      worst = std::max(worst, e);
      likelihoods += e < 1e-9;
      genotypes += g.genotype == gt;
      gqs += g.gq == gq;
      sqs += std::fabs(g.qual - sq) < 1e-6 * std::max(1.0, sq);
    }
    s.require("SVTyper's fixture is read (26 rows)", rows == 26, std::to_string(rows));
    s.check("the log₁₀ likelihoods are bayes_gt's", worst, 1e-9);
    s.require("and so are the genotypes, GQs and SQs of bayesian_genotype", genotypes == rows && gqs == rows && sqs == rows,
              std::to_string(genotypes) + " " + std::to_string(gqs) + " " + std::to_string(sqs));
  }
  {
    // A tandem duplication of [a, b) with k extra copies on one haplotype: reads across the junction
    // b−1|a come only from the k copies' junctions; reads across the reference junction a−1|a come
    // from both haplotypes. With equal coverage per copy, alt/(alt + ref) = k/(k + 2).
    bool fractions = true;
    for (int k = 1; k <= 4; k++) {
      const double alt = 10.0 * k, ref = 20.0;  // 10 reads per copy of each junction
      fractions = fractions && std::fabs(alt / (alt + ref) - k / (k + 2.0)) < 1e-12;
      const SvGenotype het = genotype_of(ref, alt, k / (k + 2.0), k / (k + 1.0));
      // Homozygous: 2k junction copies against the 2 reference junctions.
      const SvGenotype hom = genotype_of(20.0, 20.0 * k, k / (k + 2.0), k / (k + 1.0));
      fractions = fractions && het.genotype == 1 && hom.genotype == 2;
    }
    s.require("a duplication of k extra copies: a heterozygous pileup at k/(k + 2) is called 0/1, a homozygous one at k/(k + 1) is called 1/1 (k = 1 … 4)", fractions);
    const SvGenotype a = genotype_of(17, 9, 0.5, 0.999), b = genotype_of(17.0, 9.0, 0.5, 0.999);
    s.require("real-valued supports are the integers' where they are integers", a.log10_likelihood == b.log10_likelihood);
    const SvGenotype frac = genotype_of(12.4, 11.7, 0.5, 0.999);
    s.require("weighted supports (Σ of 1 − 10^(−Q/10)) are genotyped: 12.4 against 11.7 is 0/1", frac.genotype == 1 && frac.gq > 20, std::to_string(frac.gq));
  }
  return s.done();
}
