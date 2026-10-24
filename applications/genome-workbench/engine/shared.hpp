// What several of the engine's operations share: a gene found by its name, the source of
// simulated reads (a genome, or a diploid made from it, kept by the worker), the model of short
// reads a request describes, and the trace of a task graph as the page draws it. The operations
// themselves are in ops_*.hpp, each compiled in its own unit (ops_*.cpp), so that they include
// this and not each other.

#pragma once

#include "genomes.hpp"

#include <genomics/simulate/mutate.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/scheduler/task_graph.hpp>

#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace workbench {

namespace detail {

// A gene by its name (its gene or locus_tag qualifier): the contig and the feature, a CDS before
// any other feature of the name.
struct GeneAt {
  std::size_t contig = 0, feature = 0;
};
inline std::optional<GeneAt> find_gene(const Genome& g, const std::string& name) {
  for (int pass = 0; pass < 2; pass++)
    for (std::size_t ci = 0; ci < g.contigs.size(); ci++) {
      const auto& c = g.contigs[ci];
      for (std::size_t i = 0; i < c.features.size(); i++) {
        const auto& f = c.features[i];
        if (f.key != "CDS" && pass == 0) continue;
        if (f.qualifier("gene") == name || f.qualifier("locus_tag") == name) return GeneAt{ci, i};
      }
    }
  return std::nullopt;
}

// The source of simulated reads: the genome's contigs, or a diploid individual mutated from
// them by a seed.
inline sim::Source read_source(const Genome& g, const Json& request) {
  std::vector<seq::Codes> contigs;
  for (const auto& c : g.contigs) contigs.push_back(c.codes);
  if (!flag(request, "diploid", false)) return sim::Source::haploid(contigs);
  sim::MutationRates rates;
  rates.snv = num(request, "snv", 1e-3);
  rates.indel = num(request, "indel", 1e-4);
  std::vector<sim::Diploid> ds;
  for (std::size_t k = 0; k < contigs.size(); k++)
    ds.push_back(sim::mutate_diploid(contigs[k], rates, static_cast<std::uint64_t>(integer(request, "seed", 1)), static_cast<std::uint32_t>(k)));
  return sim::Source::diploid(ds);
}

// The same, made once per worker and kept with the genome: every shard of a job on this worker
// draws from it, and a diploid of a bacterial genome takes a second to make.
inline std::shared_ptr<sim::Source> source_of(Worker& worker, const Genome& g, const Json& request) {
  std::string key = "source:haploid";
  if (flag(request, "diploid", false))
    key = "source:diploid:" + std::to_string(num(request, "snv", 1e-3)) + ":" + std::to_string(num(request, "indel", 1e-4)) + ":" + std::to_string(integer(request, "seed", 1));
  return worker.derived<sim::Source>(g, key, [&] { return read_source(g, request); });
}

inline sim::ShortReadModel short_model(const Json& r) {
  sim::ShortReadModel m;
  m.length = static_cast<std::size_t>(bounded(r, "length", 150, 30, 1000));
  m.error_first = num(r, "error_first", 0.001);
  m.error_last = num(r, "error_last", 0.01);
  m.indel = num(r, "read_indel", 0.0002);
  return m;
}

// The trace of a task graph as the page draws it: [task, thread, start, end] in seconds.
inline Json trace_json(const calx::scheduler::Trace& t) {
  Json rows = Json::array();
  for (const auto& r : t.records) rows.push(Json::array({r.task, r.thread, std::round(r.start * 1e6) / 1e6, std::round(r.end * 1e6) / 1e6}));
  return Json::object({{"threads", t.threads}, {"seconds", std::round(t.seconds * 1e6) / 1e6}, {"utilisation", std::round(t.utilisation() * 1e4) / 1e4}, {"tasks", rows}});
}

}  // namespace detail

}  // namespace workbench
