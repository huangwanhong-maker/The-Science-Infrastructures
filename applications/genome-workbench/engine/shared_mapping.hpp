// What the Reads view (reads.map, reads.call) and the toolbox (tool.reads.map, tool.reads.call)
// share of mapping: the kind of reads a request names, the model of simulated long reads, a
// genome's contigs as codes, and the worker's mapper of a genome for a kind of read, built once
// (its minimizer index) and kept while the genome is.

#pragma once

#include "genomes.hpp"

#include <genomics/mapping/mapper.hpp>
#include <genomics/simulate/reads.hpp>

#include <memory>
#include <string>
#include <vector>

namespace workbench {

namespace gmap = genomics::mapping;

namespace detail {

inline bool long_reads_of(const Json& request) {
  const std::string kind = str(request, "kind", "short");
  if (kind != "short" && kind != "long") throw Fault("reads are short or long");
  return kind == "long";
}

inline sim::LongReadModel long_model(const Json& r) {
  sim::LongReadModel m;
  m.mean_length = num(r, "mean_length", 8000);
  m.min_length = static_cast<std::size_t>(bounded(r, "min_length", 1000, 100, 100000));
  m.error = num(r, "long_error", 0.08);
  if (!(m.error >= 0 && m.error <= 0.25)) throw Fault("a long-read error rate is between 0 and 0.25");
  return m;
}

// The reference's contigs as codes.
inline std::vector<seq::Codes> contigs_of(const Genome& g) {
  std::vector<seq::Codes> out;
  for (const auto& c : g.contigs) out.push_back(c.codes);
  return out;
}

inline std::shared_ptr<gmap::Mapper> mapper_of(Worker& worker, const Genome& g, bool long_reads, double* seconds = nullptr) {
  return worker.derived<gmap::Mapper>(g, long_reads ? "mapper:long" : "mapper:short", [&] {
    return gmap::Mapper(worker.pool, contigs_of(g), long_reads ? gmap::MapperOptions::long_reads() : gmap::MapperOptions::short_reads());
  }, seconds);
}

}  // namespace detail

}  // namespace workbench
