// Checks of the caller on a simulated diploid: a 500 kb reference with interspersed repeats and
// tandem arrays, mutated with SNVs and short indels (heterozygous and homozygous), sequenced as
// 30× of 150-base reads with 1% substitutions, mapped by the mapper, piled up and called. SNVs
// and indels against the truth after normalisation: precision, recall, F₁ and genotype
// concordance; the pileup's depth and its left-aligned indels; the calls the same on one and
// eight threads; the calls written as VCF and read back.

#include <genomics/mapping/mapper.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/mutate.hpp>
#include <genomics/simulate/reads.hpp>
#include <genomics/variant/caller.hpp>
#include <genomics/variant/evaluate.hpp>
#include <genomics/variant/pileup.hpp>
#include <genomics/variant/vcf.hpp>

#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <chrono>
#include <cstdio>
#include <string>
#include <vector>

using namespace genomics;
using sequence::Codes;

namespace {

double seconds_since(std::chrono::steady_clock::time_point t0) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }

std::string summary(const variant::Counts& c) {
  char buf[200];
  std::snprintf(buf, sizeof buf, "TP %zu FP %zu FN %zu: precision %.4f recall %.4f F1 %.4f, genotypes %.4f", c.tp, c.fp, c.fn, c.precision(), c.recall(), c.f1(),
                c.genotype_concordance());
  return buf;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::scheduler::ThreadPool one(1), eight(8);

  simulate::GenomeSpec spec;
  spec.length = 500000;
  spec.gc = 0.42;
  spec.repeats = {6, 4, 1000, 0.03};
  spec.tandem = {6, 40, 10, 0.02};
  const Codes reference = simulate::make_genome(spec, 21).sequence;
  simulate::MutationRates rates;
  rates.snv = 1e-3;
  rates.indel = 4e-4;
  rates.indel_mean = 2.5;
  rates.indel_max = 12;
  const simulate::Diploid diploid = simulate::mutate_diploid(reference, rates, 22);
  const std::vector<Codes> contigs{reference};
  const auto source = simulate::Source::diploid({diploid});
  simulate::ShortReadModel model;
  model.error_first = 0.005;
  model.error_last = 0.015;
  const std::size_t count = 30 * spec.length / model.length;
  const auto reads = simulate::short_reads(source, count, model, 23);

  auto t0 = std::chrono::steady_clock::now();
  const mapping::Mapper mapper(eight, contigs, mapping::MapperOptions::short_reads());
  const auto mapped = mapper.map_all(eight, reads);
  const double map_time = seconds_since(t0);
  const auto views = variant::read_views(reads);

  t0 = std::chrono::steady_clock::now();
  const auto calls = variant::call_variants(eight, contigs, mapped, views);
  const double call_time = seconds_since(t0);
  const auto calls1 = variant::call_variants(one, contigs, mapped, views);

  const auto truth = variant::truth_set(diploid.truth);
  const auto cmp = variant::compare(calls, truth, contigs);
  std::size_t truth_snv = 0, truth_indel = 0;
  for (const auto& v : diploid.truth) (v.snv() ? truth_snv : truth_indel)++;

  char buf[256];
  std::snprintf(buf, sizeof buf, "%zu SNVs and %zu indels in the truth; %zu reads mapped in %.2f s, called in %.2f s on 8 threads", truth_snv, truth_indel, reads.size(), map_time, call_time);
  s.require("a diploid of 500 kb at 30×", truth_snv > 300 && truth_indel > 120, buf);
  s.require("SNVs: precision above 0.99 — " + summary(cmp.snv), cmp.snv.precision() > 0.99);
  s.require("SNVs: recall above 0.98", cmp.snv.recall() > 0.98);
  s.require("SNVs: genotypes concordant above 0.98", cmp.snv.genotype_concordance() > 0.98);
  s.require("insertions: " + summary(cmp.insertion) + "; precision and recall above 0.9", cmp.insertion.precision() > 0.9 && cmp.insertion.recall() > 0.9);
  s.require("deletions: " + summary(cmp.deletion) + "; precision and recall above 0.9", cmp.deletion.precision() > 0.9 && cmp.deletion.recall() > 0.9);
  s.require("indels: genotypes concordant above 0.9", cmp.insertion.genotype_concordance() > 0.9 && cmp.deletion.genotype_concordance() > 0.9);
  std::string missed;
  for (std::size_t i = 0; i < cmp.false_negatives.size() && i < 6; i++)
    missed += " " + std::to_string(cmp.false_negatives[i].position + 1) + ":" + cmp.false_negatives[i].ref + ">" + cmp.false_negatives[i].alt;
  std::string extra;
  for (std::size_t i = 0; i < cmp.false_positives.size() && i < 6; i++)
    extra += " " + std::to_string(cmp.false_positives[i].position + 1) + ":" + cmp.false_positives[i].ref + ">" + cmp.false_positives[i].alt;
  s.require("all variants: " + summary(cmp.all), cmp.all.f1() > 0.98, "missed" + (missed.empty() ? std::string(" none") : missed) + "; extra" + (extra.empty() ? std::string(" none") : extra));

  bool same = calls.size() == calls1.size();
  for (std::size_t i = 0; same && i < calls.size(); i++)
    same = calls[i].position == calls1[i].position && calls[i].ref == calls1[i].ref && calls[i].alts == calls1[i].alts && calls[i].genotype == calls1[i].genotype &&
           calls[i].qual == calls1[i].qual && calls[i].gq == calls1[i].gq;
  s.require("the calls are the same on one and eight threads", same, std::to_string(calls.size()) + " calls");

  // The pileup: depth near 30 away from variants, and the reads of a truth indel reporting it at its normal place.
  {
    std::vector<mapping::Mapping> flat;
    for (const auto& ms : mapped)
      for (const auto& m : ms) flat.push_back(m);
    const auto pu = variant::pileup(reference, 0, 200000, 210000, flat, views);
    double depth = 0;
    for (std::size_t p = 200000; p < 210000; p++) depth += static_cast<double>(pu.depth(p));
    depth /= 10000;
    std::snprintf(buf, sizeof buf, "mean depth %.1f", depth);
    s.require("the pileup's depth is near the 30× sampled (bases of quality ≥ 13, MAPQ ≥ 20)", depth > 24 && depth < 34, buf);
    std::size_t checked = 0, at_place = 0;
    for (const auto& v : diploid.truth) {
      if (v.snv() || v.position < 1000 || v.position + 1000 > reference.size()) continue;
      const auto n = variant::normalize(v.position, v.ref, v.alt, reference);
      const auto local = variant::pileup(reference, 0, v.position > 200 ? v.position - 200 : 0, v.position + 200, flat, views);
      std::size_t here = 0, elsewhere = 0;
      for (const auto& e : local.indels) {
        const int length = static_cast<int>(n.alts[0].size()) - static_cast<int>(n.ref.size());
        if (e.length != length) continue;
        (e.position == n.position ? here : elsewhere)++;
      }
      checked++;
      at_place += here > 0 && here >= 4 * elsewhere;
    }
    std::snprintf(buf, sizeof buf, "%zu of %zu indels", at_place, checked);
    s.require("the reads of a truth indel report it at its normalised place (left-aligned in the pileup)", checked > 30 && at_place >= checked * 9 / 10, buf);
  }

  // The calls as VCF and back.
  {
    const auto header = variant::caller_header({"sim"}, {reference.size()}, "diploid");
    std::vector<variant::VcfRecord> records;
    for (const auto& c : calls) records.push_back(variant::vcf_record(c, {"sim"}));
    std::vector<sequence::Problem> problems;
    const auto vcf = variant::read_vcf(variant::write_vcf(header, records), &problems);
    bool back = problems.empty() && vcf.records.size() == calls.size();
    for (std::size_t i = 0; back && i < calls.size(); i++) {
      const auto& r = vcf.records[i];
      const auto g = r.genotype(0);
      back = r.position == calls[i].position && r.ref == calls[i].ref && r.alts == calls[i].alts && g.alleles.size() == 2 && g.alleles[0] == calls[i].genotype[0] &&
             g.alleles[1] == calls[i].genotype[1] && r.passed() == calls[i].passed() && r.sample_value(0, "DP") == std::to_string(calls[i].depth);
    }
    s.require("the calls written as VCF read back record for record", back, std::to_string(vcf.records.size()) + " records");
  }

  return s.done();
}
