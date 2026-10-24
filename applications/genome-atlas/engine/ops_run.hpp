// The operations of the Run view: a job cut into shards (simulated reads, with the profile of
// their errors and qualities), and the scaling of a computation with the number of threads,
// traced task by task on calx's task graph.
//
// A job's three parts (cluster.py): `reads.sample.plan` cuts n reads into shards of reads
// [first, first + count); `reads.sample.shard` simulates one shard (read i from its own stream,
// so a shard is the same on any worker) and sums what is measured; `reads.sample.merge` adds the
// shards' sums in order. The shard's own work is a task graph of chunks on the worker's threads,
// whose trace travels back for the timeline.

#pragma once

#include "genomes.hpp"
#include "shared.hpp"

#include <genomics/sequence/kmer.hpp>
#include <genomics/simulate/mutate.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/scheduler/task_graph.hpp>

#include <algorithm>
#include <chrono>
#include <string>
#include <vector>

namespace atlas {

namespace detail {

// What a shard of reads measures, summed.
struct ReadSums {
  std::size_t reads = 0, bases = 0, errors = 0, minimizers = 0;
  std::vector<double> quality, wrong;  // by position in the read: Phred sums, and errors (substitutions vs the source)
  std::array<std::size_t, 5> composition{};
  void add(const ReadSums& o) {
    reads += o.reads; bases += o.bases; errors += o.errors; minimizers += o.minimizers;
    if (quality.size() < o.quality.size()) { quality.resize(o.quality.size()); wrong.resize(o.wrong.size()); }
    for (std::size_t i = 0; i < o.quality.size(); i++) { quality[i] += o.quality[i]; wrong[i] += o.wrong[i]; }
    for (int k = 0; k < 5; k++) composition[static_cast<std::size_t>(k)] += o.composition[static_cast<std::size_t>(k)];
  }
};

inline Json sums_json(const ReadSums& s) {
  return Json::object({{"reads", s.reads}, {"bases", s.bases}, {"errors", s.errors}, {"minimizers", s.minimizers}, {"quality", numbers(s.quality, 1)},
                       {"wrong", numbers(s.wrong, 0)}, {"composition", Json(std::vector<double>(s.composition.begin(), s.composition.end()))}});
}
inline ReadSums sums_of(const Json& j) {
  ReadSums s;
  s.reads = static_cast<std::size_t>(integer(j, "reads", 0));
  s.bases = static_cast<std::size_t>(integer(j, "bases", 0));
  s.errors = static_cast<std::size_t>(integer(j, "errors", 0));
  s.minimizers = static_cast<std::size_t>(integer(j, "minimizers", 0));
  for (const auto& v : j["quality"].as_array()) s.quality.push_back(v.is_number() ? v.as_number() : 0);
  for (const auto& v : j["wrong"].as_array()) s.wrong.push_back(v.is_number() ? v.as_number() : 0);
  const auto& c = j["composition"].as_array();
  for (std::size_t k = 0; k < 5 && k < c.size(); k++) s.composition[k] = static_cast<std::size_t>(c[k].as_number());
  return s;
}

// Reads [first, first + count) simulated and measured, as a task graph of chunks on the pool.
inline std::pair<ReadSums, calx::scheduler::Trace> sample_reads(const sim::Source& src, const Json& request, std::size_t first, std::size_t count,
                                                                calx::scheduler::ThreadPool& pool) {
  const auto model = short_model(request);
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
  // Four chunks a thread, so that the threads stay busy to the end; at least 200 reads a chunk.
  const std::size_t chunk = std::max<std::size_t>(200, (count + 4 * pool.size() - 1) / (4 * pool.size()));
  const std::size_t chunks = (count + chunk - 1) / chunk;
  std::vector<ReadSums> part(chunks);
  calx::scheduler::TaskGraph graph;
  for (std::size_t k = 0; k < chunks; k++) {
    graph.add("reads", [&, k](std::size_t) {
      const std::size_t a = first + k * chunk, n = std::min(chunk, first + count - a);
      ReadSums s;
      s.quality.assign(model.length, 0);
      s.wrong.assign(model.length, 0);
      for (const auto& r : sim::short_reads(src, n, model, seed, a)) {
        s.reads++;
        s.bases += r.bases.size();
        s.errors += r.errors;
        for (std::size_t i = 0; i < r.bases.size() && i < model.length; i++) {
          s.quality[i] += seq::phred_of(r.quality[i]);
          s.composition[r.bases[i] < 4 ? r.bases[i] : 4]++;
        }
        // Errors by position, against the template when the read has no indel.
        const auto& hap = src.copies[static_cast<std::size_t>(r.haplotype)][r.contig];
        if (r.end - r.begin == r.bases.size() && !flag(request, "diploid", false)) {
          seq::Codes t(hap.begin() + static_cast<std::ptrdiff_t>(r.begin), hap.begin() + static_cast<std::ptrdiff_t>(r.end));
          if (r.reverse) t = seq::reverse_complement(t);
          for (std::size_t i = 0; i < t.size() && i < model.length; i++) s.wrong[i] += t[i] != r.bases[i];
        }
        s.minimizers += seq::minimizers(r.bases, 15, 10).size();
      }
      part[k] = std::move(s);
    }, {}, static_cast<double>(std::min(chunk, count - k * chunk)));
  }
  const auto trace = graph.run(pool);
  ReadSums total;
  total.quality.assign(model.length, 0);
  total.wrong.assign(model.length, 0);
  for (const auto& p : part) total.add(p);
  return {total, trace};
}

}  // namespace detail

ATLAS_OP("reads.sample.plan") {
  const auto g = worker.genome(request);
  const auto n = static_cast<std::size_t>(bounded(request, "reads", 100000, 1, 50'000'000));
  const auto shards = static_cast<std::size_t>(bounded(request, "shards", 16, 1, 4096));
  Json list = Json::array();
  for (std::size_t k = 0; k < shards; k++) {
    const auto r = calx::scheduler::part(n, shards, k);
    if (r.size()) list.push(Json::object({{"first", r.begin}, {"count", r.size()}, {"size", r.size()}}));
  }
  return Json::object({{"shards", list}, {"reads", n}, {"genome_length", g->length()}});
}

ATLAS_OP("reads.sample.shard") {
  const auto g = worker.genome(request);
  const Json& s = request["shard"];
  const auto [sums, trace] = detail::sample_reads(*detail::source_of(worker, *g, request), request, static_cast<std::size_t>(integer(s, "first", 0)), static_cast<std::size_t>(integer(s, "count", 0)), worker.pool);
  Json out = detail::sums_json(sums);
  out.set("_trace", detail::trace_json(trace));
  return out;
}

ATLAS_OP("reads.sample.merge") {
  detail::ReadSums total;
  for (const auto& p : request["parts"].as_array()) total.add(detail::sums_of(p));
  const double reads = std::max<double>(1, static_cast<double>(total.reads));
  std::vector<double> quality(total.quality.size()), wrong(total.wrong.size());
  for (std::size_t i = 0; i < quality.size(); i++) quality[i] = total.quality[i] / reads;
  // Per-position error rates among the reads compared (haploid sources).
  for (std::size_t i = 0; i < wrong.size(); i++) wrong[i] = total.wrong[i] / reads;
  Json out = detail::sums_json(total);
  out.set("mean_quality", numbers(quality, 2));
  out.set("error_rate", numbers(wrong, 5));
  out.set("error_per_base", std::round(1e6 * static_cast<double>(total.errors) / std::max<double>(1, static_cast<double>(total.bases))) / 1e6);
  return out;
}

// The same computation on 1, 2, 4, ... threads: seconds and throughput, and the trace of each.
ATLAS_OP("run.scaling") {
  const auto g = worker.genome(request);
  const auto n = static_cast<std::size_t>(bounded(request, "reads", 40000, 100, 5'000'000));
  std::vector<std::size_t> threads;
  for (std::size_t t = 1; t <= std::max<std::size_t>(1, worker.pool.size()); t *= 2) threads.push_back(t);
  if (threads.back() != worker.pool.size()) threads.push_back(worker.pool.size());
  Json runs = Json::array();
  double one = 0;
  const auto source = detail::source_of(worker, *g, request);
  for (std::size_t t : threads) {
    calx::scheduler::ThreadPool pool(t);
    const auto t0 = std::chrono::steady_clock::now();
    const auto [sums, trace] = detail::sample_reads(*source, request, 0, n, pool);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (t == 1) one = seconds;
    runs.push(Json::object({{"threads", t}, {"seconds", std::round(seconds * 1e5) / 1e5}, {"reads_per_second", std::round(static_cast<double>(n) / seconds)},
                            {"speedup", one > 0 ? std::round(100 * one / seconds) / 100 : 1}, {"trace", detail::trace_json(trace)}, {"reads", sums.reads}}));
  }
  return Json::object({{"reads", n}, {"runs", runs}, {"work", "simulate short reads and sketch each by its (10, 15)-minimizers"}});
}

}  // namespace atlas
