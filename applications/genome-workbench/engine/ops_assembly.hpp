// The operations of the Assembly view: reads simulated from a genome, counted, and assembled on
// the de Bruijn graph, with what the page draws of it (the spectrum of counts and its fit, the
// graph before and after every step of its cleaning, the contigs placed on the genome beside its
// repeats); and, for learning, the directed graph of a short text's k-mers, an Eulerian trail
// through it, and how many texts share its composition (the BEST theorem, exactly).
//
// Every computation is the domain's (genomics::assembly, genomics::kmers): the engine only
// simulates the input, places the contigs for the drawing and reports.

#pragma once

#include "genomes.hpp"
#include "shared_assembly.hpp"

#include <genomics/annotation/repeats.hpp>
#include <genomics/assembly/assemble.hpp>
#include <genomics/assembly/eulerian.hpp>
#include <genomics/sequence/kmer.hpp>
#include <genomics/simulate/reads.hpp>

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace workbench {

namespace detail {

// A graph as the page draws it: each unitig's length and coverage, the links between unitigs
// (unoriented, once each), and which unitigs the next step of cleaning removes.
inline Json assembly_frame(const asmb::Graph& g, const std::vector<asmb::Removal>* next) {
  Json unitigs = Json::array();
  for (std::uint32_t u = 0; u < g.unitigs.size(); u++) unitigs.push(Json::array({g.length(u), std::round(g.coverage(u) * 100) / 100}));
  std::set<std::pair<std::uint32_t, std::uint32_t>> seen;
  Json links = Json::array();
  for (const auto& l : g.links()) {
    const auto key = std::minmax(l.from, l.to);
    if (seen.insert({key.first, key.second}).second) links.push(Json::array({key.first, key.second}));
  }
  Json removed = Json::array();
  if (next)
    for (const auto& r : *next) removed.push(r.unitig);
  return Json::object({{"unitigs", unitigs}, {"links", links}, {"removed", removed}});
}

// Where each contig lies on the genome, by the first of its k-mers that occurs once in the genome
// (on either strand): its interval and strand. Contigs with no such k-mer are not placed.
inline Json place_contigs(const seq::Codes& genome, const std::vector<seq::Codes>& contigs, int k, double& covered) {
  std::unordered_map<seq::Kmer, std::int64_t> at;  // −1 where a k-mer occurs more than once
  seq::for_each_kmer(genome, k, [&](std::size_t i, seq::Kmer f, seq::Kmer) {
    auto [it, fresh] = at.emplace(f, static_cast<std::int64_t>(i));
    if (!fresh) it->second = -1;
  });
  const std::int64_t L = static_cast<std::int64_t>(genome.size());
  std::vector<std::pair<std::int64_t, std::int64_t>> spans;
  Json out = Json::array();
  for (std::size_t c = 0; c < contigs.size(); c++) {
    const auto& x = contigs[c];
    std::optional<Json> placed;
    seq::for_each_kmer(x, k, [&](std::size_t i, seq::Kmer f, seq::Kmer r) {
      if (placed) return;
      const auto n = static_cast<std::int64_t>(x.size()), off = static_cast<std::int64_t>(i);
      if (auto it = at.find(f); it != at.end() && it->second >= 0) {
        const std::int64_t b = it->second - off;
        placed = Json::array({std::max<std::int64_t>(0, b), std::min(L, b + n), false});
      } else if (auto jt = at.find(r); jt != at.end() && jt->second >= 0) {
        const std::int64_t e = jt->second + k + off;
        placed = Json::array({std::max<std::int64_t>(0, e - n), std::min(L, e), true});
      }
    });
    if (!placed) continue;
    const auto b = static_cast<std::int64_t>((*placed)[0].as_number()), e = static_cast<std::int64_t>((*placed)[1].as_number());
    spans.push_back({b, e});
    out.push(Json::object({{"contig", c}, {"length", x.size()}, {"begin", b}, {"end", e}, {"reverse", (*placed)[2].as_bool()}}));
  }
  // The share of the genome under some contig: the union of the intervals.
  std::sort(spans.begin(), spans.end());
  std::int64_t total = 0, reach = 0;
  for (const auto& [b, e] : spans) {
    if (e <= reach) continue;
    total += e - std::max(b, reach);
    reach = e;
  }
  covered = L ? static_cast<double>(total) / static_cast<double>(L) : 0;
  return out;
}

}  // namespace detail

// Reads of a genome at a coverage, with errors, counted and assembled; the steps of cleaning.
BENCH_OP("assembly.run") {
  const auto g = worker.genome(request);
  const auto& contig = g->contig(integer(request, "contig", 0));
  const seq::Codes& genome = contig.codes;
  if (genome.size() > 6'000'000) throw Fault("the Assembly view assembles genomes of up to 6 Mb");
  if (seq::count_ambiguous(genome)) throw Fault("this genome has bases other than A, C, G and T; the simulated reads would carry them");
  const int k = static_cast<int>(bounded(request, "k", 31, 3, 31));
  if (k % 2 == 0) throw Fault("k is odd, so that no k-mer is its own reverse complement");
  const double coverage = num(request, "coverage", 30), error = num(request, "error", 0.005);
  if (!(coverage >= 1 && coverage <= 200)) throw Fault("a coverage is between 1 and 200");
  if (!(error >= 0 && error <= 0.05)) throw Fault("an error rate is between 0 and 0.05");
  sim::ShortReadModel model;
  model.length = static_cast<std::size_t>(bounded(request, "length", 150, 50, 300));
  model.error_first = model.error_last = error;  // substitutions, the same along the read
  model.indel = error / 25;
  const auto count = static_cast<std::size_t>(std::ceil(coverage * static_cast<double>(genome.size()) / static_cast<double>(model.length)));
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
  const auto t0 = std::chrono::steady_clock::now();
  const auto reads = short_reads(worker, sim::Source::haploid({genome}), count, model, seed);
  const double simulating = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  asmb::AssemblyOptions o;
  o.k = k;
  const bool frames = genome.size() <= 400'000;
  o.cleaning.keep_graphs = frames;
  const auto a = asmb::assemble(reads.size(), [&](std::size_t i) -> const seq::Codes& { return reads[i].bases; }, worker.pool, o, genome.size());

  Json out = Json::object();
  out.set("length", genome.size());
  out.set("reads", reads.size());
  out.set("k", k);
  detail::assembly_report(a, out);
  // The graph before cleaning and after each step, each with what the next step removes.
  Json graphs = Json::array();
  const auto& S = a.cleaning.steps;
  if (frames && a.raw.unitigs.size() <= 4000) {
    graphs.push(detail::assembly_frame(a.raw, S.empty() ? nullptr : &S[0].removed));
    for (std::size_t i = 0; i < S.size(); i++) graphs.push(detail::assembly_frame(S[i].after, i + 1 < S.size() ? &S[i + 1].removed : nullptr));
  }
  out.set("graphs", graphs);
  double covered = 0;
  out.set("placements", detail::place_contigs(genome, a.contigs, k, covered));
  out.set("covered", covered);
  // The genome's repeats of k bases or more: a simulated genome's own record of them, or the
  // maximal repeats on both strands found by the suffix array.
  Json repeats = Json::array();
  if (!g->repeats.empty()) {
    for (const auto& r : g->repeats) repeats.push(Json::array({r.begin, r.end, r.tandem}));
  } else {
    genomics::annotation::RepeatOptions ro;
    ro.min_length = static_cast<std::size_t>(k);
    ro.both_strands = true;
    for (const auto& r : genomics::annotation::maximal_repeats(genome, ro)) {
      for (const auto& oc : r.occurrences) repeats.push(Json::array({oc.position, oc.position + r.length, false}));
      if (repeats.size() > 4000) break;
    }
  }
  out.set("repeats", repeats);
  Json seconds = detail::assembly_seconds(a);
  seconds.set("simulate", simulating);
  out.set("seconds", seconds);
  out.set("threads", worker.pool.size());
  return out;
}

// The directed de Bruijn graph of a text's k-mers: its (k − 1)-mers, its k-mers as edges, an
// Eulerian trail and the text it spells, and the number of trails and of texts (BEST).
BENCH_OP("assembly.eulerian") {
  const std::string text = str(request, "text");
  const int k = static_cast<int>(bounded(request, "k", 3, 2, 16));
  if (text.size() < static_cast<std::size_t>(k) + 1 || text.size() > 400) throw Fault("a text of k + 1 to 400 letters");
  const seq::Codes s = seq::encode_dna(text);
  if (seq::count_ambiguous(s)) throw Fault("a text is written in A, C, G and T");
  const auto kg = asmb::kmer_graph(asmb::composition(s, k));
  Json nodes = Json::array(), edges = Json::array();
  for (const auto x : kg.nodes) nodes.push(seq::kmer_string(x, k - 1));
  for (std::size_t e = 0; e < kg.graph.edges.size(); e++)
    edges.push(Json::array({kg.graph.edges[e][0], kg.graph.edges[e][1], seq::kmer_string(kg.edge_kmers[e], k)}));
  Json out = Json::object({{"text", text}, {"k", k}, {"nodes", nodes}, {"edges", edges}});
  if (const auto trail = morphe::combinatorics::eulerian_trail(kg.graph)) {
    out.set("trail", Json(std::vector<double>(trail->begin(), trail->end())));
    if (const auto spelled = asmb::spell(kg, *trail)) out.set("spelled", seq::decode_dna(*spelled));
  }
  out.set("trails", morphe::combinatorics::count_eulerian_trails(kg.graph).str());
  out.set("texts", asmb::count_strings(kg).str());
  return out;
}

}  // namespace workbench
