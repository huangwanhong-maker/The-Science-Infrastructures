// Genotype matrices from VCF: the records of the genomics domain's VCF reader
// (genomics/variant/vcf.hpp) turned into the Genotypes every statistic of this module reads.
//
// - Sites. A record is a site when it is biallelic: one ALT allele (not '.', not several).
//   Records with no ALT or with several are skipped and counted; splitting multi-allelic records
//   into biallelic ones (bcftools norm -m −) is left to the caller of this function, because the
//   statistics here read a site as one alternate allele against the rest.
// - Dosages. Each sample's GT (VCF 4.2, §1.6.2) gives its copies of the ALT allele: 0/0 → 0,
//   0/1 and 1/0 → 1, 1/1 → 2, phased or not. A call with a missing allele (./., 0/., .|1), with
//   another ploidy than two (a haploid 1, a triploid 0/0/1), naming an allele the record does
//   not have (2 at a biallelic site), or absent (no GT, or a missing sample column) is missing
//   (−1), and the last three are counted.
// - Positions are POS − 1, 0-based, as the reader keeps them; chromosomes are indexed in the order
//   they first appear, with their names, and the sites sorted by (chromosome, position) when the
//   records are not.
// - Phase. A record is phased when every call in it is a phased diploid call (0|1, 1|1, .|1);
//   when every site kept is phased, the haplotypes are kept too (copies 2i and 2i + 1 of sample
//   i, −1 for a missing allele), since a Genotypes holds haplotypes for all its sites or none.
//
// Complexity: O(n m) for n samples and m records.

#pragma once

#include <evolution/population/genotypes.hpp>

#include <genomics/variant/vcf.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

namespace evolution::population {

// What the conversion read, kept and could not use.
struct VcfConversion {
  std::size_t records = 0, sites = 0;
  std::size_t not_biallelic = 0;   // records skipped: no ALT, or several
  std::size_t other_ploidy = 0;    // calls of another ploidy than two, made missing
  std::size_t unknown_allele = 0;  // calls naming an allele the record does not have, made missing
  std::size_t absent = 0;          // calls without a GT, made missing
  std::size_t phased_records = 0;  // sites every call of which is phased
  bool phased = false;             // the haplotypes were kept
};

// The genotypes of `samples` (sample columns in their order) from VCF records.
inline Genotypes genotypes_from_vcf(const std::vector<genomics::variant::VcfRecord>& records, const std::vector<std::string>& samples, VcfConversion* report = nullptr) {
  VcfConversion r;
  r.records = records.size();
  const std::size_t n = samples.size();
  std::vector<std::string> names;          // chromosomes, by first appearance
  std::vector<std::uint32_t> chromosome;   // of each site kept
  std::vector<std::uint64_t> position;
  std::vector<std::int8_t> dosage, haplotype;
  std::size_t phased_sites = 0;
  for (const auto& rec : records) {
    if (rec.alts.size() != 1) { r.not_biallelic++; continue; }
    auto it = std::find(names.begin(), names.end(), rec.chrom);
    if (it == names.end()) it = names.insert(names.end(), rec.chrom);
    chromosome.push_back(static_cast<std::uint32_t>(it - names.begin()));
    position.push_back(rec.position);
    bool phased = true;
    for (std::size_t i = 0; i < n; i++) {
      const std::string_view gt = i < rec.samples.size() ? rec.sample_value(i, "GT") : std::string_view{};
      std::int8_t a = MISSING, b = MISSING;
      if (gt.empty()) {
        r.absent++;
        phased = false;
      } else {
        const genomics::variant::VcfGenotype g = genomics::variant::parse_genotype(gt);
        if (g.alleles.size() != 2) {
          r.other_ploidy++;
          phased = false;
        } else if (g.alleles[0] > 1 || g.alleles[1] > 1) {
          r.unknown_allele++;
          phased = phased && g.phased;
        } else {
          a = static_cast<std::int8_t>(g.alleles[0]);
          b = static_cast<std::int8_t>(g.alleles[1]);
          phased = phased && g.phased;
        }
      }
      dosage.push_back(a < 0 || b < 0 ? MISSING : static_cast<std::int8_t>(a + b));
      haplotype.push_back(a);
      haplotype.push_back(b);
    }
    phased_sites += phased;
  }
  const std::size_t m = position.size();
  r.sites = m;
  r.phased_records = phased_sites;
  r.phased = m > 0 && n > 0 && phased_sites == m;
  // Sorted by (chromosome, position), stably, when the records were not.
  std::vector<std::size_t> order(m);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::stable_sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) {
    return chromosome[x] != chromosome[y] ? chromosome[x] < chromosome[y] : position[x] < position[y];
  });
  std::vector<std::int8_t> d(m * n), h(r.phased ? m * 2 * n : 0);
  std::vector<std::uint64_t> pos(m);
  std::vector<std::uint32_t> chr(m);
  for (std::size_t j = 0; j < m; j++) {
    const std::size_t s = order[j];
    std::copy_n(dosage.begin() + static_cast<std::ptrdiff_t>(s * n), n, d.begin() + static_cast<std::ptrdiff_t>(j * n));
    if (r.phased) std::copy_n(haplotype.begin() + static_cast<std::ptrdiff_t>(s * 2 * n), 2 * n, h.begin() + static_cast<std::ptrdiff_t>(j * 2 * n));
    pos[j] = position[s];
    chr[j] = chromosome[s];
  }
  Genotypes g = genotypes_from_dosages(n, std::move(d), std::move(pos), std::move(chr), std::move(names));
  g.haplotypes = std::move(h);
  g.names = samples;
  if (report) *report = r;
  return g;
}

// The genotypes of every sample of a VCF, named by its header.
inline Genotypes genotypes_from_vcf(const genomics::variant::Vcf& vcf, VcfConversion* report = nullptr) {
  return genotypes_from_vcf(vcf.records, vcf.header.samples, report);
}

}  // namespace evolution::population
