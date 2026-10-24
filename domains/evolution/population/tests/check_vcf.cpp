// Checks of genotypes from VCF: the VCF 4.2 specification's example (§1.1) gives its two
// biallelic sites with their dosages, unphased, the multi-allelic and monomorphic records
// skipped; a VCF written here from known genotypes on two chromosomes, phased, with missing
// calls, a haploid call, an allele the record lacks and a multi-allelic record, reads back to
// the same dosages, haplotypes, positions, chromosomes and names; unsorted records are sorted;
// and the allele counts of the result are those of the genotypes written.

#include <evolution/population/genotypes.hpp>
#include <evolution/population/vcf.hpp>

#include <genomics/variant/vcf.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <string>
#include <vector>

using namespace evolution::population;
namespace variant = genomics::variant;

int main() {
  calx::verification::Suite s;

  // ---- the specification's example ---------------------------------------------------------------------
  {
    const std::string text =
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003\n"
        "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.\n"
        "20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3\n"
        "20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4\n"
        "20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2\n"
        "20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3\n";
    // The records of the VCF 4.2 specification (samtools/hts-specs, §1.1), with one of its meta lines.
    const auto vcf = variant::read_vcf(text);
    VcfConversion report;
    const Genotypes g = genotypes_from_vcf(vcf, &report);
    s.require("the specification's example gives its two biallelic sites, the other three skipped", g.sites == 2 && report.records == 5 && report.not_biallelic == 3 && g.samples == 3);
    s.require("dosages from GT: 0|0 1|0 1/1 → 0 1 2; 0|0 0|1 0/0 → 0 1 0", g.site(0)[0] == 0 && g.site(0)[1] == 1 && g.site(0)[2] == 2 && g.site(1)[1] == 1 && g.site(1)[2] == 0);
    s.require("positions POS − 1, chromosome 20, the samples named", g.positions == std::vector<std::uint64_t>{14369, 17329} && g.chromosomes == std::vector<std::string>{"20"} &&
                                                                       g.names == std::vector<std::string>{"NA00001", "NA00002", "NA00003"});
    s.require("partly phased records keep no haplotypes", !g.phased() && report.phased_records == 0);
  }

  // ---- genotypes written and read back -------------------------------------------------------------------
  {
    calx::random::Stream r(20261024, 0x56434647, 0, 0);
    const std::size_t n = 7, m = 40;
    // Phased haplotypes on two chromosomes, some alleles missing.
    std::vector<std::int8_t> hap(m * 2 * n);
    for (auto& x : hap) x = r.uniform() < 0.03 ? MISSING : static_cast<std::int8_t>(r.bernoulli(0.3));
    std::vector<std::uint64_t> pos(m);
    std::vector<std::uint32_t> chr(m);
    for (std::size_t j = 0; j < m; j++) {
      chr[j] = j < 25 ? 0 : 1;
      pos[j] = (j < 25 ? j : j - 25) * 1000 + 17;
    }
    Genotypes truth = genotypes_from_haplotypes(n, hap, pos, chr, {"chrA", "chrB"});
    std::vector<std::string> names;
    for (std::size_t i = 0; i < n; i++) names.push_back("s" + std::to_string(i));
    truth.names = names;
    variant::VcfHeader h;
    h.samples = names;
    std::vector<variant::VcfRecord> records;
    for (std::size_t j = 0; j < m; j++) {
      variant::VcfRecord rec;
      rec.chrom = truth.chromosomes[truth.chromosome[j]];
      rec.position = truth.positions[j];
      rec.ref = "A";
      rec.alts = {"G"};
      rec.format = {"GT"};
      for (std::size_t i = 0; i < n; i++) {
        const auto a = truth.haplotypes_at(j)[2 * i], b = truth.haplotypes_at(j)[2 * i + 1];
        rec.samples.push_back({variant::genotype_string({{a, b}, true})});
      }
      records.push_back(std::move(rec));
    }
    // A multi-allelic record between them, skipped.
    variant::VcfRecord multi = records[3];
    multi.position += 1;
    multi.alts = {"G", "T"};
    records.insert(records.begin() + 4, multi);
    const std::string text = variant::write_vcf(h, records);
    std::vector<genomics::sequence::Problem> problems;
    const auto vcf = variant::read_vcf(text, &problems);
    VcfConversion report;
    const Genotypes g = genotypes_from_vcf(vcf, &report);
    s.require("a phased VCF written from known genotypes reads back to them: dosages, haplotypes, positions, chromosomes, names",
              problems.empty() && g.sites == m && g.dosage == truth.dosage && g.phased() && g.haplotypes == truth.haplotypes && g.positions == truth.positions &&
                  g.chromosome == truth.chromosome && g.chromosomes == truth.chromosomes && g.names == truth.names,
              std::to_string(g.sites) + " sites");
    s.require("the multi-allelic record is skipped and counted", report.not_biallelic == 1 && report.phased_records == m);
    const auto c1 = allele_counts(g), c2 = allele_counts(truth);
    s.require("the allele counts are those of the genotypes written", c1.alt == c2.alt && c1.total == c2.total);

    // Unphased calls, a haploid call, an allele the record lacks, a missing column, records out of order.
    auto mixed = records;
    mixed.erase(mixed.begin() + 4);
    mixed[0].samples[0] = {"0/1"};
    mixed[1].samples[2] = {"1"};
    mixed[2].samples[3] = {"0|2"};
    mixed[5].samples.pop_back();
    std::swap(mixed[6], mixed[7]);
    std::swap(mixed[10], mixed[30]);  // across the chromosomes
    VcfConversion rep2;
    const Genotypes u = genotypes_from_vcf(mixed, names, &rep2);
    bool sorted = true;
    for (std::size_t j = 1; j < u.sites; j++)
      sorted = sorted && (u.chromosome[j] > u.chromosome[j - 1] || (u.chromosome[j] == u.chromosome[j - 1] && u.positions[j] > u.positions[j - 1]));
    s.require("records out of order are sorted by chromosome and position", sorted && u.positions == truth.positions);
    s.require("a haploid call, an allele the record lacks and a missing sample column are missing and counted",
              rep2.other_ploidy == 1 && rep2.unknown_allele == 1 && rep2.absent == 1 && u(2, 1) == MISSING && u(3, 2) == MISSING && u(n - 1, 5) == MISSING);
    s.require("an unphased call makes its site unphased, and the matrix keeps no haplotypes", !u.phased() && rep2.phased_records == m - 3 && u(0, 0) == 1,
              std::to_string(rep2.phased_records) + " phased sites");
  }

  return s.done();
}
