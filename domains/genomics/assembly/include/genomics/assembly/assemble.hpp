// De novo assembly of short reads on the de Bruijn graph, from reads (or counts) to contigs,
// and the statistics by which assemblies are compared.
//
// The pipeline: count the canonical k-mers of the reads (kmers/counter.hpp); read the error
// threshold from the spectrum's first minimum (kmers/spectrum.hpp) unless one is given; keep
// the solid k-mers, those counted at least that often; build and compact their de Bruijn graph
// (debruijn.hpp); clean it of tips, bubbles and low-coverage unitigs (clean.hpp); and report
// the cleaned graph's unitigs as the contigs, longest first. This is the classical short-read
// pipeline (Pevzner, Tang and Waterman 2001, PNAS 98, 9748–9753; Zerbino and Birney 2008;
// Compeau, Pevzner and Tesler 2011) without the use of read pairs or of read threading, so
// contigs end wherever the graph branches: at every repeat longer than k − 1 bases, and at
// the sites the cleaning could not resolve.
//
// - Statistics (Lander and Waterman 1988 for the model; Miller, Koren and Sutton 2010,
//   Genomics 95, 315–327, and QUAST, Gurevich et al. 2013, Bioinformatics 29, 1072–1075, for
//   the definitions): total length, the largest contig, N50 the length ℓ such that contigs of
//   length at least ℓ hold at least half the total, L50 the fewest contigs that do; NG50 and
//   LG50 the same against half a given genome size (0 when the contigs do not reach it).
// - Lander and Waterman (1988, "Genomic mapping by fingerprinting random clones: a
//   mathematical analysis", Genomics 2, 231–239): N reads of length L placed uniformly on G
//   bases, at coverage c = NL/G, two consecutive reads joining when they overlap by at least
//   T = θL bases, make N e^{−cσ} islands (contigs) in expectation, σ = 1 − θ, of mean length
//   L((e^{cσ} − 1)/c + 1 − σ). On the node-centric de Bruijn graph two reads join when they
//   share k − 1 bases, so T = k − 1, for error-free reads of a genome without repeats of k − 1
//   bases.
//
// Cost: the counting's O(bases); the graph's O(solid k-mers); the cleaning's O(rounds ·
// unitigs + total length). Times are measured and reported; memory is the counter's peak and
// the graph build's.

#pragma once

#include <genomics/assembly/clean.hpp>
#include <genomics/assembly/debruijn.hpp>
#include <genomics/kmers/counter.hpp>
#include <genomics/kmers/spectrum.hpp>

#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace genomics::assembly {

struct AssemblyStats {
  std::size_t contigs = 0, total = 0, largest = 0, n50 = 0, l50 = 0, ng50 = 0, lg50 = 0;
};

// The statistics of contigs of the given lengths; NG50 and LG50 against genome_size when it
// is given.
inline AssemblyStats statistics(std::vector<std::size_t> lengths, std::size_t genome_size = 0) {
  AssemblyStats s;
  std::sort(lengths.begin(), lengths.end(), std::greater<>());
  s.contigs = lengths.size();
  s.total = std::accumulate(lengths.begin(), lengths.end(), std::size_t{0});
  s.largest = lengths.empty() ? 0 : lengths.front();
  const auto at = [&](std::size_t half, std::size_t& n, std::size_t& l) {
    std::size_t sum = 0;
    for (std::size_t i = 0; i < lengths.size(); i++) {
      sum += lengths[i];
      if (2 * sum >= half) { n = lengths[i]; l = i + 1; return; }
    }
  };
  if (s.total) at(s.total, s.n50, s.l50);
  if (genome_size) at(genome_size, s.ng50, s.lg50);
  return s;
}
inline AssemblyStats statistics(const std::vector<Codes>& contigs, std::size_t genome_size = 0) {
  std::vector<std::size_t> lengths;
  lengths.reserve(contigs.size());
  for (const auto& c : contigs) lengths.push_back(c.size());
  return statistics(std::move(lengths), genome_size);
}

// Lander and Waterman's expected number of islands, and their mean length in bases.
inline double lander_waterman_contigs(double reads, double length, double genome, double overlap) {
  const double c = reads * length / genome, sigma = 1 - overlap / length;
  return reads * std::exp(-c * sigma);
}
inline double lander_waterman_contig_length(double reads, double length, double genome, double overlap) {
  const double c = reads * length / genome, sigma = 1 - overlap / length;
  return length * ((std::exp(c * sigma) - 1) / c + 1 - sigma);
}

struct AssemblyOptions {
  int k = 31;
  std::uint32_t min_count = 0;    // solid k-mers are counted at least this often; 0: the spectrum's first minimum
  kmers::CounterOptions counter;  // its k is set to k
  CleaningOptions cleaning;
  bool clean = true;
  std::size_t min_contig = 0;     // shorter unitigs are not reported as contigs
};

struct Assembly {
  kmers::Spectrum spectrum;
  kmers::SpectrumFit fit;
  std::uint32_t threshold = 0;
  std::size_t solid_kmers = 0;
  Graph raw;           // the compacted graph of the solid k-mers
  Cleaning cleaning;   // its steps, and the cleaned graph
  std::vector<Codes> contigs;  // longest first
  AssemblyStats stats;
  struct Seconds {
    double count = 0, spectrum = 0, graph = 0, clean = 0;
    double total() const { return count + spectrum + graph + clean; }
  } seconds;
  std::size_t count_bytes = 0, graph_bytes = 0;  // the counter's peak, and the graph build's
};

namespace detail {
inline double since(std::chrono::steady_clock::time_point t) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count(); }
}  // namespace detail

// The assembly of the k-mers counted (the counter's k must be the options' k).
inline Assembly assemble(const kmers::KmerCounter& counts, calx::scheduler::ThreadPool& pool, const AssemblyOptions& o = {}, std::size_t genome_size = 0) {
  check_odd_k(o.k, "assemble");
  if (counts.k() != o.k) throw std::invalid_argument("genomics::assembly::assemble: the counter's k is not the options' k");
  Assembly a;
  a.count_bytes = counts.peak_bytes();
  auto t = std::chrono::steady_clock::now();
  a.spectrum = kmers::spectrum(counts, pool);
  kmers::FitOptions fo;
  fo.threshold = o.min_count;
  a.fit = kmers::fit_spectrum(a.spectrum, fo);
  a.threshold = a.fit.threshold;
  a.seconds.spectrum = detail::since(t);
  t = std::chrono::steady_clock::now();
  const std::vector<KmerCount> solid = counts.sorted(a.threshold, pool);
  a.solid_kmers = solid.size();
  BuildReport report;
  a.raw = build_graph(solid, o.k, pool, &report);
  a.graph_bytes = report.peak_bytes + solid.size() * sizeof(KmerCount);
  a.seconds.graph = detail::since(t);
  t = std::chrono::steady_clock::now();
  if (o.clean) a.cleaning = clean(a.raw, o.cleaning);
  else a.cleaning.graph = a.raw;
  a.seconds.clean = detail::since(t);
  const Graph& g = a.cleaning.graph;
  std::vector<std::uint32_t> order;
  for (std::uint32_t u = 0; u < g.size(); u++)
    if (g.length(u) >= o.min_contig) order.push_back(u);
  std::stable_sort(order.begin(), order.end(), [&](std::uint32_t x, std::uint32_t y) { return g.length(x) > g.length(y); });
  a.contigs.reserve(order.size());
  for (auto u : order) a.contigs.push_back(g.unitigs[u].sequence);
  a.stats = statistics(a.contigs, genome_size);
  return a;
}

// The assembly of reads 0, …, n − 1, get(i) giving the i-th as a const Codes& (for
// simulate::Read, its bases).
template <class Get>
Assembly assemble(std::size_t n, Get&& get, calx::scheduler::ThreadPool& pool, const AssemblyOptions& o = {}, std::size_t genome_size = 0) {
  kmers::CounterOptions co = o.counter;
  co.k = o.k;
  const auto t = std::chrono::steady_clock::now();
  kmers::KmerCounter counts(co);
  counts.add(n, get, pool);
  counts.release_buffers();
  const double counting = detail::since(t);
  Assembly a = assemble(counts, pool, o, genome_size);
  a.seconds.count = counting;
  return a;
}
inline Assembly assemble(const std::vector<Codes>& reads, calx::scheduler::ThreadPool& pool, const AssemblyOptions& o = {}, std::size_t genome_size = 0) {
  return assemble(reads.size(), [&](std::size_t i) -> const Codes& { return reads[i]; }, pool, o, genome_size);
}

}  // namespace genomics::assembly
