// What the Assembly view (assembly.run, on simulated reads) and the toolbox (tool.assembly.reads, on
// the user's reads) report alike of an assembly on the de Bruijn graph: the spectrum of k-mer counts
// as far as its fit reaches, the error threshold and the peak, the k-mer coverage and the size of
// the genome the spectrum implies, the steps of cleaning, the contigs' statistics, the time of each
// stage and the memory.

#pragma once

#include "genomes.hpp"

#include <genomics/assembly/assemble.hpp>

#include <algorithm>
#include <cmath>
#include <string>

namespace workbench {

namespace asmb = genomics::assembly;

namespace detail {

inline void assembly_report(const asmb::Assembly& a, Json& out) {
  // The spectrum: distinct k-mers by their count, as far as the fit's interest reaches.
  Json spectrum = Json::array();
  const std::size_t last = std::min<std::size_t>(a.spectrum.histogram.size(), std::max<std::size_t>(24, 4 * static_cast<std::size_t>(a.fit.peak) + 4));
  for (std::size_t c = 1; c < last; c++)
    if (a.spectrum.at(c)) spectrum.push(Json::array({c, a.spectrum.at(c)}));
  out.set("spectrum", spectrum);
  out.set("threshold", a.threshold);
  out.set("peak", a.fit.peak);
  out.set("kmer_coverage", std::round(a.fit.coverage * 100) / 100);
  out.set("size_estimate", std::round(a.fit.genome_size));
  out.set("error_kmers", a.fit.error_kmers);
  out.set("solid_kmers", a.solid_kmers);
  Json steps = Json::array();
  for (const auto& s : a.cleaning.steps)
    steps.push(Json::object({{"name", std::string(asmb::step_name(s.kind)) + ", round " + std::to_string(s.round + 1)}, {"kind", asmb::step_name(s.kind)},
                             {"removed", s.removed.size()}, {"before", s.unitigs_before}, {"after", s.unitigs_after}}));
  out.set("steps", steps);
  out.set("raw_unitigs", a.raw.unitigs.size());
  out.set("stats", Json::object({{"contigs", a.stats.contigs}, {"total", a.stats.total}, {"largest", a.stats.largest}, {"n50", a.stats.n50}, {"l50", a.stats.l50},
                                 {"ng50", a.stats.ng50}, {"lg50", a.stats.lg50}}));
  out.set("memory", Json::object({{"count", a.count_bytes}, {"graph", a.graph_bytes}}));
}

inline Json assembly_seconds(const asmb::Assembly& a) {
  return Json::object({{"count", a.seconds.count}, {"spectrum", a.seconds.spectrum}, {"graph", a.seconds.graph}, {"clean", a.seconds.clean}});
}

}  // namespace detail

}  // namespace workbench
