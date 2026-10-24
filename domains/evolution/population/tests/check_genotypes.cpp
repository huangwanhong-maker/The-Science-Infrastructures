// Checks of genotype matrices: the three constructors (site by site, the samples × sites matrix,
// phased haplotypes) agree; malformed input is refused; allele counts, frequencies, minor-allele
// frequencies and missingness on a matrix worked by hand; the filters by frequency and
// missingness; selecting sites and samples keeps the haplotypes with the dosages; and diploids
// paired from the coalescent's haplotypes.

#include <evolution/population/genotypes.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

using namespace evolution::population;

int main() {
  calx::verification::Suite s;

  // A matrix by hand: 3 samples (rows) × 4 sites (columns).
  //   sample 0:  0  1  2  −1
  //   sample 1:  1  1  0   0
  //   sample 2:  2 −1  0   0
  const std::vector<std::int8_t> matrix = {0, 1, 2, -1, 1, 1, 0, 0, 2, -1, 0, 0};
  const Genotypes g = genotypes_from_matrix(3, 4, matrix, {10, 20, 30, 40});
  {
    const Genotypes h = genotypes_from_dosages(3, {0, 1, 2, 1, 1, -1, 2, 0, 0, -1, 0, 0}, {10, 20, 30, 40});
    s.require("the matrix and its sites agree: G(i, j) is the dosage of sample i at site j", g.dosage == h.dosage && g(2, 0) == 2 && g(0, 3) == MISSING);
    // The same as phased haplotypes (sample i's copies 2i, 2i + 1); a missing copy makes a missing dosage.
    const Genotypes p = genotypes_from_haplotypes(3, {0, 0, 1, 0, 1, 1, /**/ 0, 1, 0, 1, -1, 1, /**/ 1, 1, 0, 0, 0, 0, /**/ -1, 0, 0, 0, 0, 0},
                                                  {10, 20, 30, 40});
    s.require("phased haplotypes give the same dosages, and are kept", p.dosage == g.dosage && p.phased() && p.haplotypes_at(1)[3] == 1);
    s.require("one chromosome named \"1\" by default", g.chromosomes.size() == 1 && g.chromosomes[0] == "1" && g.chromosome_sites(0).second == 4);
  }
  {
    int refused = 0;
    const auto refuses = [&](auto&& make) {
      try {
        make();
      } catch (const std::invalid_argument&) {
        refused++;
      }
    };
    refuses([] { genotypes_from_dosages(2, {0, 3}, {1}); });                      // a dosage of 3
    refuses([] { genotypes_from_dosages(2, {0, 1, 0, 1}, {5, 4}); });             // positions out of order
    refuses([] { genotypes_from_dosages(2, {0, 1, 0}, {5}); });                   // not samples × sites
    refuses([] { genotypes_from_dosages(1, {0, 1}, {1, 2}, {1, 0}); });           // chromosomes out of order
    refuses([] { genotypes_from_haplotypes(1, {0, 2}, {1}); });                   // an allele of 2
    s.require("malformed matrices are refused with invalid_argument", refused == 5);
  }
  {
    const AlleleCounts c = allele_counts(g);
    s.require("allele counts at each site: alternate copies and chromosomes typed", c.alt == std::vector<std::uint32_t>{3, 2, 2, 0} &&
                                                                                      c.total == std::vector<std::uint32_t>{6, 4, 6, 4});
    const std::vector<std::size_t> two = {1, 2};
    const AlleleCounts d = allele_counts(g, two);
    s.require("and over a subset of the samples", d.alt == std::vector<std::uint32_t>{3, 1, 0, 0} && d.total == std::vector<std::uint32_t>{4, 2, 4, 4});
    const auto f = frequencies(c), maf = minor_allele_frequencies(c);
    s.require("frequencies and minor-allele frequencies", f[0] == 0.5 && f[2] == 1.0 / 3 && maf[2] == 1.0 / 3 && maf[3] == 0);
    const AlleleCounts none = allele_counts(select_samples(g, std::vector<std::size_t>{0}));
    s.require("a site with nothing typed has frequency NaN", std::isnan(frequencies(none)[3]));
    const auto ms = missing_by_site(g), mi = missing_by_sample(g);
    s.require("missingness by site and by sample", ms == std::vector<double>{0, 1.0 / 3, 0, 1.0 / 3} && mi == std::vector<double>{0.25, 0, 0.25});
    s.require("filters: minor-allele frequency at least 0.4", passing_sites(g, 0.4) == std::vector<std::size_t>{0, 1});
    s.require("and at most 20% missing", passing_sites(g, 0.4, 0.2) == std::vector<std::size_t>{0});
    s.require("monomorphic sites fail any positive frequency filter", passing_sites(g, 1e-9) == std::vector<std::size_t>{0, 1, 2});
  }
  {
    const Genotypes p = genotypes_from_haplotypes(3, {0, 0, 1, 0, 1, 1, /**/ 0, 1, 0, 1, -1, 1, /**/ 1, 1, 0, 0, 0, 0, /**/ -1, 0, 0, 0, 0, 0},
                                                  {10, 20, 30, 40});
    const Genotypes sites = select_sites(p, std::vector<std::size_t>{1, 3});
    const Genotypes samples = select_samples(p, std::vector<std::size_t>{2, 0});
    s.require("selecting sites keeps their dosages, haplotypes and positions",
              sites.sites == 2 && sites(0, 0) == 1 && sites.haplotypes_at(1)[0] == -1 && sites.positions == std::vector<std::uint64_t>{20, 40});
    s.require("selecting samples reorders them, with their haplotypes", samples.samples == 2 && samples(0, 0) == 2 && samples(1, 0) == 0 &&
                                                                         samples.haplotypes_at(1)[0] == -1 && samples.haplotypes_at(1)[1] == 1);
  }
  {
    Haplotypes h;
    h.sequences = 4;
    h.sites = 2;
    h.alleles = {1, 0, 0, 1, /**/ 1, 1, 0, 0};
    h.positions = {0.25, 0.75};
    const Genotypes d = diploids(h, 1000);
    s.require("diploids pair sequences 2i and 2i + 1, positions ⌊x · L⌋", d.samples == 2 && d(0, 0) == 1 && d(1, 0) == 1 && d(0, 1) == 2 && d(1, 1) == 0 &&
                                                                              d.positions == std::vector<std::uint64_t>{250, 750} && d.phased());
    const Haplotypes back = haplotypes_of(d);
    s.require("and their haplotypes are the sequences again", back.alleles == h.alleles && back.sequences == 4);
  }
  return s.done();
}
