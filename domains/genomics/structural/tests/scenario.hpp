// A simulated sample run end to end for the checks: the individual made, its pairs drawn, mapped and
// detected, and the detection compared with the truth, with a report printed as comment lines.

#pragma once

#include <genomics/mapping/mapper.hpp>
#include <genomics/sequence/fasta.hpp>
#include <genomics/sequence/genbank.hpp>
#include <genomics/simulate/pairs.hpp>
#include <genomics/simulate/structural.hpp>
#include <genomics/structural/detect.hpp>
#include <genomics/structural/evaluate.hpp>

#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace scenario {

using namespace genomics;
using sequence::Codes;

inline double seconds_since(std::chrono::steady_clock::time_point t0) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }

inline std::string tally(const structural::Tally& t) {
  char buf[200];
  std::snprintf(buf, sizeof buf, "TP %zu FP %zu FN %zu: precision %.3f recall %.3f, genotypes %.3f", t.tp_call, t.fp, t.fn, t.precision(), t.recall(), t.genotype_concordance());
  return buf;
}

// The data library's cache (SCIENCE_CACHE when set).
inline std::string cache() {
  const char* e = std::getenv("SCIENCE_CACHE");
  return e && *e ? e : GENOMICS_CACHE;
}

inline std::string read_file(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  if (!f) return {};
  std::ostringstream s;
  s << f.rdbuf();
  return s.str();
}

struct Run {
  simulate::StructuralSample sample;
  std::vector<simulate::ReadPair> pairs;
  structural::MateMappings mates;
  structural::Detection detection;
  std::vector<structural::TruthEvent> truth;
  structural::SvComparison comparison;
  structural::CopyConcordance copy;
  double seconds_simulate = 0, seconds_index = 0, seconds_map = 0, seconds_detect = 0;
};

inline Run run(calx::scheduler::ThreadPool& pool, const std::vector<Codes>& contigs, const std::vector<simulate::StructuralVariant>& svs, double coverage,
               const simulate::PairModel& pm, std::uint64_t seed, const structural::DetectOptions& options = {}) {
  Run r;
  auto t0 = std::chrono::steady_clock::now();
  r.sample = simulate::make_structural_sample(contigs, svs, simulate::MutationRates{}, seed);
  r.pairs = simulate::read_pairs(r.sample.chromosomes, simulate::pairs_for_coverage(r.sample.chromosomes, coverage, pm.read.length), pm, seed + 1);
  r.seconds_simulate = seconds_since(t0);
  t0 = std::chrono::steady_clock::now();
  const mapping::Mapper mapper(pool, contigs, mapping::MapperOptions::short_reads());
  r.seconds_index = seconds_since(t0);
  t0 = std::chrono::steady_clock::now();
  r.mates = structural::map_pairs(pool, mapper, r.pairs);
  r.seconds_map = seconds_since(t0);
  t0 = std::chrono::steady_clock::now();
  r.detection = structural::detect(pool, contigs, r.mates, r.pairs, options);
  r.seconds_detect = seconds_since(t0);
  r.truth = structural::truth_events(r.sample.truth);
  r.comparison = structural::compare_sv(contigs, r.detection.calls, r.truth);
  r.copy = structural::copy_concordance(r.detection.depth, r.detection.segments, [&](std::uint32_t c, std::size_t p) { const auto k = r.sample.copies_at(c, p); return k[0] + k[1]; });
  return r;
}

// Σ (n − e)²/e per bin over the mappable bins whose true copy number is 2 throughout.
inline double dispersion_at_two(const Run& r) {
  double sum = 0;
  std::size_t n = 0;
  for (const auto& b : r.detection.depth.bins) {
    if (!b.mappable || b.exposure <= 0) continue;
    bool two = true;
    for (std::size_t p = b.begin; p < b.end && two; p += 100) { const auto k = r.sample.copies_at(b.contig, p); two = k[0] + k[1] == 2; }
    if (!two) continue;
    sum += (b.count - b.exposure) * (b.count - b.exposure) / b.exposure;
    n++;
  }
  return n ? sum / static_cast<double>(n) : 0;
}

inline void report(const Run& r, const char* title) {
  const auto& d = r.detection;
  const auto& cmp = r.comparison;
  std::printf("# %s: %zu pairs simulated in %.2f s, index %.2f s, mapped in %.2f s, detected in %.2f s (library %.2f, evidence %.2f, cluster %.2f, call %.2f, depth %.2f, segment %.2f, snv %.2f, allelic %.2f)\n",
              title, r.pairs.size(), r.seconds_simulate, r.seconds_index, r.seconds_map, r.seconds_detect, d.seconds.library, d.seconds.evidence, d.seconds.cluster, d.seconds.call,
              d.seconds.depth, d.seconds.segment, d.seconds.snv, d.seconds.allelic);
  std::printf("# library: median %.1f MAD %.1f, normal [%lld, %lld], cut-offs [%lld, %lld], read length %zu\n", d.library.median, d.library.mad, static_cast<long long>(d.library.min_normal),
              static_cast<long long>(d.library.max_normal), static_cast<long long>(d.library.min_cutoff), static_cast<long long>(d.library.max_cutoff), d.library.read_length);
  std::printf("# pairs:");
  for (int k = 0; k < 9; k++) std::printf(" %s %zu,", structural::pair_class_name(static_cast<structural::PairClass>(k)), d.counts.by_class[k]);
  std::printf("\n# evidence %zu, clips %zu, breakpoints %zu, calls %zu\n", d.evidence.size(), d.clips.size(), d.breakpoints.size(), d.calls.size());
  std::printf("# overall: %s\n", tally(cmp.overall).c_str());
  for (const auto& [k, t] : cmp.by_kind) std::printf("#   %s: %s\n", k.c_str(), tally(t).c_str());
  for (std::size_t k = 0; k < 5; k++) std::printf("#   size %s: %s\n", structural::size_bin_name(k), tally(cmp.by_size[k]).c_str());
  if (std::getenv("STRUCTURAL_VERBOSE")) {
    for (std::size_t j = 0; j < r.truth.size(); j++) {
      const auto& tr = r.truth[j];
      const structural::SvMatch* m = nullptr;
      for (const auto& x : cmp.matches) if (x.truth == j) m = &x;
      std::printf("#   truth %s %u:[%zu, %zu) length %zu gt %d partner %u:%zu inv %d", structural::sv_kind_name(tr.kind), tr.contig, tr.begin, tr.end, tr.length, tr.genotype, tr.partner, tr.partner_position, tr.inverted);
      if (m) {
        const auto& c = d.calls[m->call];
        std::printf(" <- %s [%zu, %zu) gt %d (ref %.1f alt %.1f) copies %.2f precise %d error %zu partner %u:%zu", structural::sv_kind_name(c.kind), c.begin, c.end, c.genotype.genotype, c.genotype.reference,
                    c.genotype.alternate, c.copies, m->precise, m->breakpoint_error, c.partner, c.partner_position);
      }
      std::printf("\n");
    }
  }
  for (const auto j : cmp.false_negatives) {
    const auto& tr = r.truth[j];
    std::printf("#   missed %s %u:[%zu, %zu) length %zu gt %d partner %u:%zu\n", structural::sv_kind_name(tr.kind), tr.contig, tr.begin, tr.end, tr.length, tr.genotype, tr.partner, tr.partner_position);
  }
  for (const auto i : cmp.false_positives) {
    const auto& c = d.calls[i];
    std::printf("#   false %s %u:[%zu, %zu) length %zu qual %.1f pairs %zu splits %zu clips %zu partner %u:%zu copies %.2f\n", structural::sv_kind_name(c.kind), c.contig, c.begin, c.end, c.length, c.qual,
                c.pairs, c.splits, c.clips, c.partner, c.partner_position, c.copies);
  }
  std::size_t exact = 0, worst = 0;
  for (const auto e : cmp.breakpoint_errors) { exact += e == 0; worst = std::max(worst, e); }
  std::printf("# precise matches %zu: exact %zu, worst %zu bases\n", cmp.breakpoint_errors.size(), exact, worst);
  std::printf("# segments:");
  for (const auto& g : d.segments)
    if (g.copies != 2) std::printf(" [%u:%zu–%zu cn %d (%.2f)]", g.contig, g.begin, g.end, g.copies, g.estimate());
  std::printf("\n# aneuploidy:");
  for (const auto& a : d.aneuploidy) std::printf(" [%u: %d copies over %.3f]", a.contig, a.copies, a.fraction);
  std::printf("\n# allelic:");
  for (const auto& a : d.allelic) std::printf(" [%u:%zu–%zu %s cn %d minor %d sites %zu]", a.contig, a.begin, a.end, structural::allelic_state_name(a.state), a.copies, a.minor, a.sites);
  std::printf("\n# truth beyond junctions:");
  for (const auto& v : r.sample.truth)
    if (v.type == simulate::SvType::loss_of_heterozygosity || v.type == simulate::SvType::aneuploidy || v.length() >= 50000)
      std::printf(" [%s %u:%zu–%zu on %d%d]", simulate::sv_type_name(v.type), v.contig, v.begin, v.end, v.on[0], v.on[1]);
  std::printf("\n# copy number concordance %.4f over %zu bins; dispersion %.3f per bin at two copies; small variants %zu, sites %zu, runs %zu\n", r.copy.fraction(), r.copy.bins,
              dispersion_at_two(r), d.small.size(), d.sites.size(), d.runs.size());
}

// The bases of the truth's copy-neutral LOH and of the called one, as (contig, begin, end).
inline void loh_intervals(const Run& r, std::vector<std::array<std::size_t, 3>>& called, std::vector<std::array<std::size_t, 3>>& truth) {
  for (const auto& v : r.sample.truth)
    if (v.type == simulate::SvType::loss_of_heterozygosity) truth.push_back({v.contig, v.begin, v.end});
  for (const auto& a : r.detection.allelic)
    if (a.state == structural::AllelicState::copy_neutral_loh) called.push_back({a.contig, a.begin, a.end});
}

}  // namespace scenario
