// Checks of the mapper on a simulated genome of two contigs, with diverged interspersed repeats,
// tandem arrays and three exact copies of a 2.5 kb segment (one reverse-complemented): short
// reads (150 bases, 1% substitutions), seeded by minimizers and by SMEMs, and long reads (8%
// errors, indel-dominated) map to their
// true places; reads inside the exact copies get mapping quality near 0; every CIGAR is
// consistent with its read and its interval and rescores to its AS, its NM the path's edit
// distance; the mappings are the same on one and eight threads; the accuracy curve by MAPQ.

#include <genomics/align/types.hpp>
#include <genomics/mapping/evaluate.hpp>
#include <genomics/mapping/mapper.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

using namespace genomics;
using mapping::Mapper;
using mapping::MapperOptions;
using mapping::Mapping;
using sequence::Codes;

namespace {

double seconds_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

bool same_mapping(const Mapping& a, const Mapping& b) {
  return a.read == b.read && a.contig == b.contig && a.target_begin == b.target_begin && a.target_end == b.target_end && a.reverse == b.reverse && a.mapq == b.mapq &&
         a.primary == b.primary && a.supplementary == b.supplementary && a.alignment.score == b.alignment.score && a.alignment.cigar == b.alignment.cigar &&
         a.alignment.query_begin == b.alignment.query_begin && a.nm == b.nm && a.chain_score == b.chain_score;
}

bool same_all(const std::vector<std::vector<Mapping>>& a, const std::vector<std::vector<Mapping>>& b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); i++) {
    if (a[i].size() != b[i].size()) return false;
    for (std::size_t j = 0; j < a[i].size(); j++)
      if (!same_mapping(a[i][j], b[i][j])) return false;
  }
  return true;
}

// Every mapping's path agrees with its read and its interval, and rescores to its claims.
std::string consistency(const Mapper& mapper, const std::vector<simulate::Read>& reads, const std::vector<std::vector<Mapping>>& all) {
  std::size_t checked = 0;
  for (std::size_t i = 0; i < reads.size(); i++) {
    const Codes rc = sequence::reverse_complement(reads[i].bases);
    for (const auto& m : all[i]) {
      const Codes& q = m.reverse ? rc : reads[i].bases;
      const Codes& t = mapper.contigs()[m.contig];
      const auto& a = m.alignment;
      if (align::query_length(m.cigar()) != reads[i].bases.size() || m.query_length != reads[i].bases.size()) return reads[i].name + ": the CIGAR's query length";
      if (align::query_length(a.cigar) != a.query_end - a.query_begin) return reads[i].name + ": the path's query length";
      if (align::target_length(a.cigar) != m.target_end - m.target_begin || m.target_end > t.size()) return reads[i].name + ": the path's target length";
      if (align::rescore(a, q, t, mapper.options().scoring) != a.score) return reads[i].name + ": AS is not the path's score";
      if (align::edit_distance(a, q, t) != m.nm) return reads[i].name + ": NM is not the path's edit distance";
      for (const auto& o : a.cigar)
        if (o.op == 'S' || o.length == 0) return reads[i].name + ": a clip or an empty operation inside the path";
      checked++;
    }
  }
  return "every one of " + std::to_string(checked) + " mappings";
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::scheduler::ThreadPool one(1), eight(8);

  // ---- the genome --------------------------------------------------------------------------------------
  simulate::GenomeSpec spec;
  spec.length = 600000;
  spec.gc = 0.45;
  spec.repeats = {8, 6, 1200, 0.02};
  spec.tandem = {4, 60, 20, 0.01};
  simulate::Genome g = simulate::make_genome(spec, 11);
  // Three exact copies of a 2.5 kb segment, the third reverse-complemented, where no repeat lies.
  const std::size_t copy_length = 2500;
  std::vector<std::pair<std::size_t, std::size_t>> exact;
  const auto free_at = [&](std::size_t b) {
    for (const auto& c : g.repeats)
      if (b < c.end + 200 && c.begin < b + copy_length + 200) return false;
    for (const auto& [x, y] : exact)
      if (b < y + 200 && x < b + copy_length + 200) return false;
    return true;
  };
  for (std::size_t b = 60000; exact.size() < 3 && b + copy_length < spec.length; b += 37003)
    if (free_at(b)) exact.push_back({b, b + copy_length});
  const Codes segment(g.sequence.begin() + static_cast<std::ptrdiff_t>(exact[0].first), g.sequence.begin() + static_cast<std::ptrdiff_t>(exact[0].second));
  std::copy(segment.begin(), segment.end(), g.sequence.begin() + static_cast<std::ptrdiff_t>(exact[1].first));
  const Codes rseg = sequence::reverse_complement(segment);
  std::copy(rseg.begin(), rseg.end(), g.sequence.begin() + static_cast<std::ptrdiff_t>(exact[2].first));
  const std::vector<Codes> contigs{g.sequence, simulate::random_sequence(150000, 0.5, 12, 1)};
  const auto source = simulate::Source::haploid(contigs);

  const auto in_repeat = [&](const simulate::Read& r) {
    if (r.contig != 0) return false;
    for (const auto& c : g.repeats)
      if (r.begin < c.end && c.begin < r.end) return true;
    for (const auto& [x, y] : exact)
      if (r.begin < y && x < r.end) return true;
    return false;
  };
  const auto inside_exact = [&](const simulate::Read& r) {
    if (r.contig != 0) return false;
    for (const auto& [x, y] : exact)
      if (x <= r.begin && r.end <= y) return true;
    return false;
  };

  // ---- short reads, seeded by minimizers and by SMEMs ------------------------------------------------------
  for (const bool smem : {false, true}) {
    const std::string seeded = smem ? " (seeded by SMEMs)" : " (seeded by minimizers)";
    const auto t0 = std::chrono::steady_clock::now();
    const Mapper mapper(eight, contigs, smem ? MapperOptions::short_reads_smem() : MapperOptions::short_reads());
    const double build = seconds_since(t0);
    simulate::ShortReadModel model;
    model.length = 150;
    model.error_first = 0.005;
    model.error_last = 0.015;
    const auto reads = simulate::short_reads(source, 20000, model, 101);
    std::size_t bases = 0;
    for (const auto& r : reads) bases += r.bases.size();

    auto t1 = std::chrono::steady_clock::now();
    const auto on8 = mapper.map_all(eight, reads);
    const double s8 = seconds_since(t1);
    t1 = std::chrono::steady_clock::now();
    const auto on1 = mapper.map_all(one, reads);
    const double s1 = seconds_since(t1);

    std::size_t outside = 0, outside_ok = 0, exact_reads = 0, exact_low = 0, exact_q = 0;
    for (std::size_t i = 0; i < reads.size(); i++) {
      const Mapping* m = mapping::primary_of(on8[i]);
      if (!in_repeat(reads[i])) {
        outside++;
        outside_ok += m && mapping::correct(*m, mapping::truth_of(reads[i]));
      }
      if (inside_exact(reads[i])) {
        exact_reads++;
        exact_low += m && m->mapq <= 3;
        exact_q += m ? static_cast<std::size_t>(m->mapq) : 0;
      }
    }
    const auto acc = mapping::accuracy_of(on8, reads);
    const double fo = static_cast<double>(outside_ok) / static_cast<double>(outside);
    char buf[256];
    std::snprintf(buf, sizeof buf, "%.4f of %zu", fo, outside);
    s.require("short reads (150 b, 1% substitutions) outside repeats map to their true place: above 0.995" + seeded, fo > 0.995, buf);
    std::snprintf(buf, sizeof buf, "%.4f mapped, %.4f correct of %zu", acc.fraction_mapped(), acc.fraction_correct(), acc.reads);
    s.require("short reads across the whole genome, repeats included, map correctly: above 0.98" + seeded, acc.fraction_correct() > 0.98, buf);
    std::snprintf(buf, sizeof buf, "%zu of %zu reads, mean MAPQ %.2f", exact_low, exact_reads, exact_reads ? static_cast<double>(exact_q) / static_cast<double>(exact_reads) : 0.0);
    s.require("short reads inside the exact copies have MAPQ at most 3" + seeded, exact_reads > 100 && exact_low == exact_reads, buf);
    const auto& q20 = acc.curve[4];
    std::snprintf(buf, sizeof buf, "MAPQ >= %d: %.4f of reads, error rate %.2e", q20.mapq, q20.fraction_mapped, q20.error_rate);
    s.require("short reads of MAPQ ≥ 20 are wrong less than once in a thousand" + seeded, q20.error_rate < 1e-3 && q20.fraction_mapped > 0.97, buf);
    std::string curve;
    for (const auto& p : acc.curve) {
      std::snprintf(buf, sizeof buf, " q%d:%.4f/%.1e", p.mapq, p.fraction_mapped, p.error_rate);
      curve += buf;
    }
    s.require("the accuracy curve of short reads by MAPQ (fraction mapped / error rate)" + seeded, true, curve);
    s.require("short-read mappings are consistent with their reads and rescore to AS and NM" + seeded, consistency(mapper, reads, on8).rfind("every", 0) == 0, consistency(mapper, reads, on8));
    s.require("short-read mappings are the same on one and eight threads" + seeded, same_all(on1, on8));
    const mapping::Throughput t8{reads.size(), bases, s8}, t1x{reads.size(), bases, s1};
    std::snprintf(buf, sizeof buf, "index %.2f s; %.0f reads/s on 1 thread, %.0f on 8 (%.1f Mb/s)", build, t1x.reads_per_second(), t8.reads_per_second(), t8.bases_per_second() / 1e6);
    s.require("short-read throughput" + seeded, true, buf);
  }

  // ---- long reads --------------------------------------------------------------------------------------
  {
    const Mapper mapper(eight, contigs, MapperOptions::long_reads());
    simulate::LongReadModel model;
    model.mean_length = 6000;
    model.error = 0.08;
    const auto reads = simulate::long_reads(source, 300, model, 202);
    std::size_t bases = 0;
    for (const auto& r : reads) bases += r.bases.size();
    auto t1 = std::chrono::steady_clock::now();
    const auto on8 = mapper.map_all(eight, reads, 1);
    const double s8 = seconds_since(t1);
    t1 = std::chrono::steady_clock::now();
    const auto on1 = mapper.map_all(one, reads, 1);
    const double s1 = seconds_since(t1);
    const auto acc = mapping::accuracy_of(on8, reads);
    char buf[256];
    std::snprintf(buf, sizeof buf, "%.4f mapped, %.4f correct of %zu", acc.fraction_mapped(), acc.fraction_correct(), acc.reads);
    s.require("long reads (8% errors) map to their true place: above 0.97", acc.fraction_correct() > 0.97, buf);
    double clipped = 0;
    std::size_t primaries = 0;
    for (std::size_t i = 0; i < reads.size(); i++)
      if (const Mapping* m = mapping::primary_of(on8[i])) {
        clipped += static_cast<double>(m->clip_left() + m->clip_right()) / static_cast<double>(m->query_length);
        primaries++;
      }
    std::snprintf(buf, sizeof buf, "%.4f of a read clipped on average", primaries ? clipped / static_cast<double>(primaries) : 0.0);
    s.require("long reads are aligned nearly end to end: under 2% of their bases clipped on average", primaries && clipped / static_cast<double>(primaries) < 0.02, buf);
    s.require("long-read mappings are consistent with their reads and rescore to AS and NM", consistency(mapper, reads, on8).rfind("every", 0) == 0, consistency(mapper, reads, on8));
    s.require("long-read mappings are the same on one and eight threads", same_all(on1, on8));
    std::snprintf(buf, sizeof buf, "%.1f reads/s on 1 thread, %.1f on 8 (%.2f Mb/s)", static_cast<double>(reads.size()) / s1, static_cast<double>(reads.size()) / s8,
                  static_cast<double>(bases) / s8 / 1e6);
    s.require("long-read throughput", true, buf);
  }

  return s.done();
}
