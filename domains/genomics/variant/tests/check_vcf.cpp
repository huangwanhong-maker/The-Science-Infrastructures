// Checks of VCF: the example of the VCF 4.2 specification (§1.1) read field by field and
// written back to the same text; genotypes parsed with their phasing, missing alleles and
// ploidy; calls written as records and read back; and lines that are not records reported as
// problems with their line numbers while the others are read.

#include <genomics/variant/caller.hpp>
#include <genomics/variant/vcf.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace genomics;
using variant::VcfGenotype;

int main() {
  calx::verification::Suite s;

  // ---- the specification's example ------------------------------------------------------------------------
  {
    // The Variant Call Format Specification, VCFv4.2 (samtools/hts-specs), §1.1, "An example",
    // columns separated by tabs as the format requires.
    const std::string text =
        "##fileformat=VCFv4.2\n"
        "##fileDate=20090805\n"
        "##source=myImputationProgramV3.1\n"
        "##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta\n"
        "##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>\n"
        "##phasing=partial\n"
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
        "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
        "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n"
        "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">\n"
        "##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">\n"
        "##FILTER=<ID=q10,Description=\"Quality below 10\">\n"
        "##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n"
        "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003\n"
        "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.\n"
        "20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3\n"
        "20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4\n"
        "20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2\n"
        "20\t1234567\tmicrosat1\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3\n";
    std::vector<sequence::Problem> problems;
    const auto vcf = variant::read_vcf(text, &problems);
    s.require("the specification's example reads without problems: 5 records, 3 samples, 17 meta lines, contig 20 of 62 435 964 bases",
              problems.empty() && vcf.records.size() == 5 && vcf.header.samples.size() == 3 && vcf.header.meta.size() == 17 && vcf.header.contigs.size() == 1 &&
                  vcf.header.contigs[0].first == "20" && vcf.header.contigs[0].second == 62435964,
              std::to_string(problems.size()) + " problems");
    if (vcf.records.size() == 5) {
      const auto& r = vcf.records;
      s.require("POS is read 1-based into 0-based positions; IDs, REF, ALT ('.' none, 'G,T' two)",
                r[0].position == 14369 && r[0].id == "rs6054257" && r[2].alts == std::vector<std::string>{"G", "T"} && r[3].alts.empty() && r[4].ref == "GTC");
      s.require("QUAL, FILTER and INFO, flags included (DB, H2)", r[1].qual == 3 && r[1].filters == std::vector<std::string>{"q10"} && r[0].passed() && r[0].info_value("AF") == "0.5" &&
                                                                       r[0].info.size() == 5 && r[0].info[3].first == "DB" && r[0].info[3].second.empty());
      s.require("genotypes with their phasing: 0|0 phased, 1/1 not, 1|2 between two alternates",
                r[0].genotype(0) == VcfGenotype{{0, 0}, true} && r[0].genotype(2) == VcfGenotype{{1, 1}, false} && r[2].genotype(0) == VcfGenotype{{1, 2}, true} &&
                    r[4].genotype(1) == VcfGenotype{{0, 2}, false});
      s.require("trailing FORMAT fields may be dropped: NA00003's HQ at 17330 is absent", r[1].sample_value(2, "HQ").empty() && r[1].sample_value(2, "DP") == "3");
    }
    const std::string back = variant::write_vcf(vcf.header, vcf.records);
    s.require("the example written back is the same text", back == text);
  }

  // ---- genotypes --------------------------------------------------------------------------------------
  {
    bool ok = true, all = true;
    all = all && variant::parse_genotype("./.", &ok).missing() && ok;
    all = all && variant::parse_genotype(".", &ok).alleles == std::vector<int>{-1} && ok;
    all = all && variant::parse_genotype("1", &ok) == VcfGenotype{{1}, false} && ok;
    all = all && variant::parse_genotype("0|.", &ok) == VcfGenotype{{0, -1}, true} && ok;
    all = all && variant::parse_genotype("0/1|2", &ok) == VcfGenotype{{0, 1, 2}, false} && ok;
    all = all && variant::parse_genotype("10|3", &ok) == VcfGenotype{{10, 3}, true} && ok;
    s.require("GT values: missing, haploid, partly missing, triploid, alleles above 9", all);
    variant::parse_genotype("0/", &ok);
    const bool trailing = !ok;
    variant::parse_genotype("A/T", &ok);
    s.require("a GT that is not one is refused (0/, A/T)", trailing && !ok);
    s.require("genotypes are written as read", variant::genotype_string({{0, 1}, true}) == "0|1" && variant::genotype_string({{-1, -1}, false}) == "./." &&
                                                    variant::genotype_string({}) == ".");
  }

  // ---- calls as records ----------------------------------------------------------------------------------
  {
    variant::Call snv;
    snv.contig = 1;
    snv.position = 999;
    snv.ref = "A";
    snv.alts = {"G"};
    snv.qual = 123.456;
    snv.depth = 31;
    snv.allele_depths = {15, 16};
    snv.genotype = {0, 1};
    snv.gq = 99;
    snv.strand_bias = 1.25;
    variant::Call del = snv;
    del.position = 1999;
    del.ref = "ACA";
    del.alts = {"A"};
    del.genotype = {1, 1};
    del.indel = true;
    del.filters = {"LowDepth", "StrandBias"};
    del.allele_depths = {0, 4};
    del.depth = 4;
    const std::vector<std::string> names{"chr1", "chr2"};
    const auto header = variant::caller_header(names, {5000, 8000}, "sample1");
    const std::string text = variant::write_vcf(header, {variant::vcf_record(snv, names), variant::vcf_record(del, names)});
    std::vector<sequence::Problem> problems;
    const auto vcf = variant::read_vcf(text, &problems);
    bool same = problems.empty() && vcf.records.size() == 2;
    if (same) {
      const auto& a = vcf.records[0];
      const auto& b = vcf.records[1];
      same = a.chrom == "chr2" && a.position == 999 && a.ref == "A" && a.alts[0] == "G" && std::fabs(a.qual - 123.46) < 1e-9 && a.passed() && a.info_value("DP") == "31" &&
             a.info_value("FS") == "1.25" && a.format == std::vector<std::string>{"GT", "GQ", "DP", "AD"} && a.sample_value(0, "AD") == "15,16" &&
             a.genotype(0) == VcfGenotype{{0, 1}, false} && a.sample_value(0, "GQ") == "99" && b.position == 1999 && b.ref == "ACA" &&
             b.filters == std::vector<std::string>{"LowDepth", "StrandBias"} && b.genotype(0) == VcfGenotype{{1, 1}, false};
    }
    s.require("calls written as VCF read back: CHROM, POS, REF, ALT, QUAL, FILTER, INFO DP and FS, GT:GQ:DP:AD", same, std::to_string(problems.size()) + " problems");
    s.require("the caller's header names its contigs, INFO, FILTER and FORMAT fields and its sample",
              vcf.header.contigs.size() == 2 && vcf.header.contigs[1].second == 8000 && vcf.header.samples == std::vector<std::string>{"sample1"} &&
                  text.find("##FORMAT=<ID=AD,Number=R") != std::string::npos);
  }

  // ---- problems ----------------------------------------------------------------------------------------
  {
    const std::string text =
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
        "1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
        "1\tzero\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
        "1\t0\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
        "1\t200\t.\tA\tG\n"
        "1\t300\t.\tA\tG\thigh\tPASS\t.\tGT\t0/1\n"
        "1\t400\t.\tA\tG\t.\t.\t.\tGT\t0/x\n"
        "1\t500\t.\tA\tG\t.\t.\t.\tGT\t0/1\t1/1\n"
        "1\t600\t.\t\tG\t.\t.\t.\tGT\t0/1\n";
    std::vector<sequence::Problem> problems;
    const auto vcf = variant::read_vcf(text, &problems);
    std::vector<std::size_t> lines;
    for (const auto& p : problems) lines.push_back(p.line);
    s.require("lines that are not records are problems with their line numbers (4, 5 POS; 6 columns; 7 QUAL; 8 GT; 9 samples; 10 REF)",
              lines == std::vector<std::size_t>{4, 5, 6, 7, 8, 9, 10}, std::to_string(problems.size()) + " problems");
    s.require("the other records are read, never throwing; a missing QUAL is NaN", vcf.records.size() == 3 && !vcf.records[1].has_qual() && vcf.records[2].samples.size() == 2);
  }

  return s.done();
}
