// The operations of the Reads view: reads mapped to their genome, and the variants of a diploid
// called from them.
//
// - `reads.map` is a job (cluster.py): `.plan` cuts n reads into shards; `.shard` simulates its
//   reads (read i from its own stream, so a shard is the same on any worker), maps them with the
//   worker's index of the genome (built once and kept) and counts, by mapping quality, the reads
//   placed and those placed wrongly against the simulator's truth; `.merge` adds the counts and
//   draws Li's curve (2018, Fig. 1): the fraction of reads mapped at a MAPQ threshold against the
//   error rate among them.
// - `reads.call` makes a diploid of the genome (SNVs and indels with their truth), reads it at a
//   coverage, maps the reads, calls the variants and compares the calls with the truth; the run
//   is kept by the worker, and `reads.window` serves the reads and calls of a stretch of it for
//   the browser.
//
// Every computation is the domain's (genomics::simulate, genomics::mapping, genomics::variant):
// the engine only gathers, times and reports.

#pragma once

#include "genomes.hpp"
#include "shared.hpp"
#include "shared_mapping.hpp"

#include <genomics/mapping/evaluate.hpp>
#include <genomics/mapping/mapper.hpp>
#include <genomics/mapping/sam.hpp>
#include <genomics/variant/caller.hpp>
#include <genomics/variant/evaluate.hpp>
#include <genomics/variant/vcf.hpp>

#include <algorithm>
#include <chrono>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace workbench {

namespace gvar = genomics::variant;

namespace detail {

// The counts of a shard: reads, bases, unmapped, and the primaries and wrong ones at each MAPQ.
struct MapSums {
  std::size_t reads = 0, bases = 0, unmapped = 0;
  std::vector<double> mapped = std::vector<double>(61, 0), wrong = std::vector<double>(61, 0);
  double map_seconds = 0;
  void add(const MapSums& o) {
    reads += o.reads; bases += o.bases; unmapped += o.unmapped; map_seconds += o.map_seconds;
    for (std::size_t q = 0; q < 61; q++) { mapped[q] += o.mapped[q]; wrong[q] += o.wrong[q]; }
  }
};

}  // namespace detail

BENCH_OP("reads.map.plan") {
  const auto g = worker.genome(request);
  const auto n = static_cast<std::size_t>(bounded(request, "reads", 200000, 1, 20'000'000));
  const auto shards = static_cast<std::size_t>(bounded(request, "shards", 16, 1, 4096));
  Json list = Json::array();
  for (std::size_t k = 0; k < shards; k++) {
    const auto r = calx::scheduler::part(n, shards, k);
    if (r.size()) list.push(Json::object({{"first", r.begin}, {"count", r.size()}, {"size", r.size()}}));
  }
  return Json::object({{"shards", list}, {"reads", n}, {"genome_length", g->length()}});
}

BENCH_OP("reads.map.shard") {
  const auto g = worker.genome(request);
  const Json& s = request["shard"];
  const auto first = static_cast<std::size_t>(integer(s, "first", 0)), count = static_cast<std::size_t>(integer(s, "count", 0));
  const bool long_reads = detail::long_reads_of(request);
  double index_seconds = 0;
  const auto mapper = detail::mapper_of(worker, *g, long_reads, &index_seconds);
  const auto source = detail::source_of(worker, *g, request);
  const sim::Source& src = *source;
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
  const auto reads = long_reads ? sim::long_reads(src, count, detail::long_model(request), seed, first) : sim::short_reads(src, count, detail::short_model(request), seed, first);
  const auto t0 = std::chrono::steady_clock::now();
  const auto mappings = mapper->map_all(worker.pool, reads, long_reads ? 1 : 64);
  detail::MapSums sums;
  sums.map_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  for (std::size_t i = 0; i < reads.size(); i++) {
    sums.reads++;
    sums.bases += reads[i].bases.size();
    const gmap::Mapping* m = gmap::primary_of(mappings[i]);
    if (!m) { sums.unmapped++; continue; }
    const auto q = static_cast<std::size_t>(std::clamp(m->mapq, 0, 60));
    sums.mapped[q]++;
    if (!gmap::correct(*m, gmap::truth_of(reads[i]))) sums.wrong[q]++;
  }
  return Json::object({{"reads", sums.reads}, {"bases", sums.bases}, {"unmapped", sums.unmapped}, {"mapped", Json(sums.mapped)}, {"wrong", Json(sums.wrong)},
                       {"map_seconds", sums.map_seconds}, {"index_seconds", index_seconds}, {"threads", worker.pool.size()}});
}

BENCH_OP("reads.map.merge") {
  detail::MapSums total;
  double index_seconds = 0;
  for (const auto& p : request["parts"].as_array()) {
    detail::MapSums s;
    s.reads = static_cast<std::size_t>(integer(p, "reads", 0));
    s.bases = static_cast<std::size_t>(integer(p, "bases", 0));
    s.unmapped = static_cast<std::size_t>(integer(p, "unmapped", 0));
    s.map_seconds = num(p, "map_seconds", 0);
    const auto& m = p["mapped"].as_array();
    const auto& w = p["wrong"].as_array();
    for (std::size_t q = 0; q < 61 && q < m.size(); q++) { s.mapped[q] = m[q].as_number(); s.wrong[q] = w[q].as_number(); }
    total.add(s);
    index_seconds = std::max(index_seconds, num(p, "index_seconds", 0));
  }
  // Li's curve: at each threshold, the primaries at or above it and the wrong among them.
  Json curve = Json::array();
  double mapped = 0, wrong = 0;
  for (int q = 60; q >= 0; q--) {
    mapped += total.mapped[static_cast<std::size_t>(q)];
    wrong += total.wrong[static_cast<std::size_t>(q)];
    curve.push(Json::array({q, mapped, wrong, total.reads ? mapped / static_cast<double>(total.reads) : 0, mapped ? wrong / mapped : 0}));
  }
  return Json::object({{"reads", total.reads}, {"bases", total.bases}, {"unmapped", total.unmapped}, {"mapped", mapped}, {"wrong", wrong},
                       {"by_mapq", Json(total.mapped)}, {"wrong_by_mapq", Json(total.wrong)}, {"curve", curve},
                       {"map_seconds", total.map_seconds}, {"index_seconds", index_seconds}});
}

namespace detail {

// A diploid read, mapped and called: kept by the worker for the browser's windows.
struct CallRun {
  std::vector<seq::Codes> reference;
  std::vector<std::string> names;
  std::vector<sim::Diploid> diploids;
  std::vector<sim::Read> reads;
  std::vector<std::vector<gmap::Mapping>> mappings;
  std::vector<std::uint32_t> by_position;  // primaries, by (contig, target_begin): indices of reads
  std::vector<gvar::Call> calls;
  std::vector<gvar::TruthVariant> truth;
  gvar::Comparison comparison;
  gmap::Accuracy accuracy;
  std::set<std::tuple<std::uint32_t, std::size_t, std::string, std::string>> false_positives;
  double seconds_simulate = 0, seconds_index = 0, seconds_map = 0, seconds_call = 0, seconds_compare = 0;
  std::size_t mean_length = 0;
};

inline std::string call_key(const Json& r) {
  std::ostringstream k;
  k << "call:" << str(r, "kind", "short") << ":" << num(r, "coverage", 30) << ":" << num(r, "snv", 1e-3) << ":" << num(r, "indel", 1e-4) << ":" << integer(r, "seed", 1)
    << ":" << num(r, "error_first", 0.001) << ":" << num(r, "error_last", 0.01) << ":" << num(r, "long_error", 0.08);
  return k.str();
}

inline std::shared_ptr<CallRun> call_run(Worker& worker, const Genome& g, const Json& request, double* seconds = nullptr) {
  return worker.derived<CallRun>(g, call_key(request), [&] {
    using clock = std::chrono::steady_clock;
    const auto since = [](clock::time_point t) { return std::chrono::duration<double>(clock::now() - t).count(); };
    CallRun run;
    run.reference = contigs_of(g);
    for (const auto& c : g.contigs) run.names.push_back(c.name);
    std::size_t L = 0;
    for (const auto& c : run.reference) L += c.size();
    if (L > 12'000'000) throw Fault("the Reads view calls genomes of up to 12 Mb");
    const bool long_reads = long_reads_of(request);
    const double coverage = num(request, "coverage", 30);
    if (!(coverage >= 2 && coverage <= 200)) throw Fault("a coverage is between 2 and 200");
    sim::MutationRates rates;
    rates.snv = num(request, "snv", 1e-3);
    rates.indel = num(request, "indel", 1e-4);
    if (!(rates.snv >= 0 && rates.snv <= 0.02 && rates.indel >= 0 && rates.indel <= 0.005)) throw Fault("rates of SNVs up to 0.02 and of indels up to 0.005 a base");
    const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
    auto t = clock::now();
    for (std::size_t k = 0; k < run.reference.size(); k++) run.diploids.push_back(sim::mutate_diploid(run.reference[k], rates, seed, static_cast<std::uint32_t>(k)));
    const sim::Source src = sim::Source::diploid(run.diploids);
    if (long_reads) {
      const auto model = long_model(request);
      const auto n = static_cast<std::size_t>(std::ceil(coverage * static_cast<double>(L) / model.mean_length));
      run.reads = sim::long_reads(src, n, model, seed + 7);
    } else {
      const auto model = short_model(request);
      const auto n = static_cast<std::size_t>(std::ceil(coverage * static_cast<double>(L) / static_cast<double>(model.length)));
      run.reads = short_reads(worker, src, n, model, seed + 7);
    }
    run.seconds_simulate = since(t);
    std::size_t bases = 0;
    for (const auto& r : run.reads) bases += r.bases.size();
    run.mean_length = run.reads.empty() ? 0 : bases / run.reads.size();
    const auto mapper = mapper_of(worker, g, long_reads, &run.seconds_index);
    t = clock::now();
    run.mappings = mapper->map_all(worker.pool, run.reads, long_reads ? 1 : 64);
    run.seconds_map = since(t);
    std::vector<gmap::Truth> truth;
    for (const auto& r : run.reads) truth.push_back(gmap::truth_of(r));
    run.accuracy = gmap::accuracy(run.mappings, truth);
    for (std::uint32_t i = 0; i < run.mappings.size(); i++)
      if (gmap::primary_of(run.mappings[i])) run.by_position.push_back(i);
    std::sort(run.by_position.begin(), run.by_position.end(), [&](std::uint32_t a, std::uint32_t b) {
      const auto& x = run.mappings[a].front();
      const auto& y = run.mappings[b].front();
      return std::tie(x.contig, x.target_begin, a) < std::tie(y.contig, y.target_begin, b);
    });
    t = clock::now();
    const auto views = gvar::read_views(run.reads);
    run.calls = gvar::call_variants(worker.pool, run.reference, run.mappings, views);
    run.seconds_call = since(t);
    t = clock::now();
    for (std::uint32_t k = 0; k < run.diploids.size(); k++)
      for (const auto& v : gvar::truth_set(run.diploids[k].truth, k)) run.truth.push_back(v);
    run.comparison = gvar::compare(run.calls, run.truth, run.reference);
    for (const auto& v : run.comparison.false_positives) run.false_positives.insert({v.contig, v.position, v.ref, v.alt});
    run.seconds_compare = since(t);
    return run;
  }, seconds);
}

inline Json counts_json(const gvar::Counts& c) {
  return Json::object({{"tp", c.tp}, {"fp", c.fp}, {"fn", c.fn}, {"precision", std::round(c.precision() * 1e5) / 1e5}, {"recall", std::round(c.recall() * 1e5) / 1e5},
                       {"f1", std::round(c.f1() * 1e5) / 1e5}, {"genotype", std::round(c.genotype_concordance() * 1e5) / 1e5}});
}

inline Json call_json(const CallRun& run, const gvar::Call& c) {
  // A call is false when the comparison found its (normalised) allele nowhere in the truth.
  bool wrong = false;
  for (const auto& a : c.alts) wrong = wrong || run.false_positives.count({c.contig, c.position, c.ref, a});
  return Json::object({{"contig", c.contig}, {"position", c.position}, {"ref", c.ref}, {"alts", c.alts}, {"qual", std::round(c.qual * 10) / 10},
                       {"genotype", std::to_string(c.genotype[0]) + "/" + std::to_string(c.genotype[1])}, {"gq", c.gq}, {"depth", c.depth},
                       {"ad", Json(std::vector<double>(c.allele_depths.begin(), c.allele_depths.end()))}, {"passed", c.passed()},
                       {"filters", c.filters}, {"indel", c.indel}, {"true", c.passed() ? Json(!wrong) : Json()}});
}

inline Json truth_json(const gvar::TruthVariant& v) {
  return Json::object({{"contig", v.contig}, {"position", v.position}, {"ref", v.ref}, {"alt", v.alt}, {"dosage", v.dosage}});
}

}  // namespace detail

// A diploid of the genome, read, mapped and called; the calls against the truth.
BENCH_OP("reads.call") {
  const auto g = worker.genome(request);
  double seconds = 0;
  const auto run = detail::call_run(worker, *g, request, &seconds);
  const auto& cmp = run->comparison;
  Json types = Json::object({{"snv", detail::counts_json(cmp.snv)}, {"insertion", detail::counts_json(cmp.insertion)}, {"deletion", detail::counts_json(cmp.deletion)},
                             {"all", detail::counts_json(cmp.all)}});
  Json calls = Json::array();
  for (const auto& c : run->calls) {
    if (calls.size() >= 2000) break;
    calls.push(detail::call_json(*run, c));
  }
  Json missed = Json::array();
  for (const auto& v : cmp.false_negatives) {
    if (missed.size() >= 500) break;
    missed.push(detail::truth_json(v));
  }
  // The first records of the SAM and of the VCF, as they would be written.
  std::vector<std::size_t> lengths;
  for (const auto& c : run->reference) lengths.push_back(c.size());
  std::ostringstream sam;
  sam << gmap::sam_header(run->names, lengths, "genome-workbench", "reads.call");
  for (std::size_t i = 0, n = 0; i < run->reads.size() && n < 12; i++)
    if (const auto* m = gmap::primary_of(run->mappings[i])) {
      sam << gmap::sam_line(gmap::sam_record(*m, run->reads[i].name, run->reads[i].bases, run->reads[i].quality, run->names));  // the line ends itself
      n++;
    }
  std::ostringstream vcf;
  {
    std::vector<gvar::VcfRecord> recs;
    for (std::size_t i = 0; i < run->calls.size() && i < 12; i++) recs.push_back(gvar::vcf_record(run->calls[i], run->names));
    vcf << gvar::write_vcf(gvar::caller_header(run->names, lengths, "diploid"), recs);
  }
  std::size_t covered = 0, depth = 0;
  for (const auto i : run->by_position) {
    const auto& m = run->mappings[i].front();
    depth += m.target_end - m.target_begin;
  }
  for (const auto& c : run->reference) covered += c.size();
  const auto& acc = run->accuracy;
  return Json::object({{"reads", run->reads.size()}, {"mean_length", run->mean_length}, {"length", covered},
                       {"mean_depth", covered ? static_cast<double>(depth) / static_cast<double>(covered) : 0},
                       {"mapped", acc.mapped}, {"correct", acc.correct}, {"unmapped", acc.unmapped},
                       {"truth", run->truth.size()}, {"called", run->calls.size()}, {"types", types}, {"calls", calls}, {"missed", missed},
                       {"sam", sam.str()}, {"vcf", vcf.str()}, {"contigs", run->names},
                       {"seconds", Json::object({{"simulate", run->seconds_simulate}, {"index", run->seconds_index}, {"map", run->seconds_map},
                                                 {"call", run->seconds_call}, {"compare", run->seconds_compare}, {"built", seconds}})},
                       {"threads", worker.pool.size()}});
}

// The reads of a stretch of the called diploid, with what differs from the reference, and the
// calls and the truth there.
BENCH_OP("reads.window") {
  const auto g = worker.genome(request);
  const auto run = detail::call_run(worker, *g, request);
  const auto contig = static_cast<std::uint32_t>(integer(request, "contig", 0));
  if (contig >= run->reference.size()) throw Fault("no such contig");
  const auto& ref = run->reference[contig];
  const auto begin = static_cast<std::size_t>(std::min<long long>(integer(request, "begin", 0), static_cast<long long>(ref.size())));
  const auto end = static_cast<std::size_t>(std::min<long long>(integer(request, "end", static_cast<long long>(begin) + 300), static_cast<long long>(ref.size())));
  if (end <= begin || end - begin > 5000) throw Fault("a window of 1 to 5000 bases");
  // The primaries that overlap the window: those that begin before its end and end after its start.
  const std::size_t reach = run->mean_length * 4 + 2000;
  const auto lo = std::lower_bound(run->by_position.begin(), run->by_position.end(), std::make_pair(contig, begin > reach ? begin - reach : 0), [&](std::uint32_t i, const auto& key) {
    const auto& m = run->mappings[i].front();
    return std::make_pair(m.contig, m.target_begin) < key;
  });
  Json reads = Json::array();
  const std::size_t cap = static_cast<std::size_t>(bounded(request, "max_reads", 1500, 1, 20000));
  std::size_t total = 0;
  for (auto it = lo; it != run->by_position.end(); ++it) {
    const auto& m = run->mappings[*it].front();
    if (m.contig != contig || m.target_begin >= end) break;
    if (m.target_end <= begin) continue;
    total++;
    if (reads.size() >= cap) continue;
    const auto& read = run->reads[*it];
    const seq::Codes bases = m.reverse ? seq::reverse_complement(read.bases) : read.bases;
    Json mismatches = Json::array(), insertions = Json::array(), deletions = Json::array();
    std::size_t q = m.alignment.query_begin, t = m.alignment.target_begin;
    for (const auto& o : m.alignment.cigar) {
      if (o.op == 'M' || o.op == '=' || o.op == 'X') {
        for (std::uint32_t k = 0; k < o.length; k++, q++, t++)
          if (t >= begin && t < end && bases[q] != ref[t]) mismatches.push(Json::array({t, std::string(1, seq::dna_letter(bases[q]))}));
      } else if (o.op == 'I') {
        if (t >= begin && t <= end) insertions.push(Json::array({t, o.length}));
        q += o.length;
      } else if (o.op == 'D') {
        if (t + o.length > begin && t < end) deletions.push(Json::array({t, o.length}));
        t += o.length;
      }
    }
    reads.push(Json::array({m.target_begin, m.target_end, m.reverse, m.mapq, mismatches, insertions, deletions, Json::array({m.clip_left(), m.clip_right()}), read.name}));
  }
  Json calls = Json::array(), truth = Json::array();
  for (const auto& c : run->calls)
    if (c.contig == contig && c.position + c.ref.size() > begin && c.position < end) calls.push(detail::call_json(*run, c));
  for (const auto& v : run->truth)
    if (v.contig == contig && v.position + v.ref.size() > begin && v.position < end) truth.push(detail::truth_json(v));
  return Json::object({{"contig", contig}, {"begin", begin}, {"end", end}, {"reference", seq::decode_dna(seq::Codes(ref.begin() + static_cast<std::ptrdiff_t>(begin), ref.begin() + static_cast<std::ptrdiff_t>(end)))},
                       {"reads", reads}, {"total", total}, {"calls", calls}, {"truth", truth}});
}

}  // namespace workbench
