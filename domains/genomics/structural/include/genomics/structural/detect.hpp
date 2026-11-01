// The detection of an abnormal genome from its paired reads, every step of the module in order:
// the library, the classes of the pairs, the evidence of junctions, breakpoints, events and their
// genotypes; the depth, its segments and copy numbers, aneuploidy; the small variants of the
// domain's caller, their B-allele fractions, the runs of homozygosity, and the alleles read with
// the copy number (loss of heterozygosity told from deletion). Each step is its header's; this one
// only runs them in order and keeps what each made, and the time it took.
//
// The mates' mappings (structural/pairs.hpp, map_pairs) must outlive the detection: its pair
// records point into them. Deterministic: every step is serial in a fixed order or parallel into
// its own slots (the depth's chunks, the caller's windows).

#pragma once

#include <genomics/structural/allelic.hpp>
#include <genomics/structural/call.hpp>
#include <genomics/structural/cluster.hpp>
#include <genomics/structural/copy_number.hpp>
#include <genomics/structural/depth.hpp>
#include <genomics/structural/evidence.hpp>
#include <genomics/structural/genotype.hpp>
#include <genomics/structural/pairs.hpp>
#include <genomics/variant/caller.hpp>
#include <genomics/variant/pileup.hpp>

#include <calx/scheduler/thread_pool.hpp>

#include <chrono>
#include <span>
#include <string_view>
#include <vector>

namespace genomics::structural {

struct DetectOptions {
  LibraryOptions library;
  int min_mapq = 20;
  EvidenceOptions evidence;
  ClusterOptions cluster;
  CallOptions call;
  DepthOptions depth;
  CopyNumberOptions copy_number;
  bool small_variants = true;  // call SNVs for B-allele fractions and runs of homozygosity
  variant::CallerOptions caller;
  RohOptions roh;
  std::size_t min_copy_call = 10000;  // copy-number segments called on their own from this length
};

struct Timings {
  double library = 0, evidence = 0, cluster = 0, call = 0, depth = 0, segment = 0, snv = 0, allelic = 0;
};

struct Detection {
  Library library;
  PairCounts counts;
  std::vector<PairRecord> pairs;
  std::vector<Evidence> evidence;
  std::vector<Clip> clips;
  std::vector<Breakpoint> breakpoints;
  std::vector<SvCall> calls;          // junction events, then copy-number calls, by position
  Depth depth;
  std::vector<CopySegmentCall> segments;
  std::vector<Aneuploidy> aneuploidy;
  std::vector<variant::Call> small;   // the caller's SNVs and indels
  std::vector<AlleleSite> sites;
  std::vector<Run> runs;
  std::vector<AllelicSegment> allelic;
  Timings seconds;
};

// Everything from the mates' mappings (2N, mate m of pair i at 2i + m) of `pairs` (anything with
// mates[m].bases and mates[m].quality: simulate::ReadPair).
template <class Pair>
Detection detect(calx::scheduler::ThreadPool& pool, const std::vector<Codes>& contigs, const MateMappings& mates, const std::vector<Pair>& pairs, const DetectOptions& o = {}) {
  using clock = std::chrono::steady_clock;
  const auto since = [](clock::time_point t) { return std::chrono::duration<double>(clock::now() - t).count(); };
  Detection d;
  auto t = clock::now();
  d.library = estimate_library(mates, o.library);
  d.pairs = classify_pairs(mates, d.library, o.min_mapq, &d.counts);
  d.seconds.library = since(t);

  t = clock::now();
  const auto bases = [&](std::size_t r) -> const Codes& { return pairs[r / 2].mates[r % 2].bases; };
  gather_evidence(contigs, mates, d.pairs, d.library, bases, d.evidence, d.clips, o.evidence);
  d.seconds.evidence = since(t);

  t = clock::now();
  d.breakpoints = cluster(contigs, d.evidence, o.cluster);
  d.seconds.cluster = since(t);

  t = clock::now();
  d.depth = read_depth(pool, contigs, mates, d.pairs, d.library, o.depth);
  d.seconds.depth = since(t);

  t = clock::now();
  d.segments = segment_copy_number(d.depth, o.copy_number);
  d.aneuploidy = find_aneuploidy(d.depth, d.segments, o.copy_number);
  d.seconds.segment = since(t);

  t = clock::now();
  const ReadIndex index(mates, d.pairs, d.library, contigs.size(), o.min_mapq);
  d.calls = interpret(contigs, d.evidence, d.breakpoints, d.clips, d.pairs, d.library, index, d.depth, o.call);
  auto cnv = copy_number_calls(d.segments, d.aneuploidy, d.calls, o.min_copy_call);
  d.calls.insert(d.calls.end(), cnv.begin(), cnv.end());
  std::stable_sort(d.calls.begin(), d.calls.end(), [](const SvCall& x, const SvCall& y) { return x.contig != y.contig ? x.contig < y.contig : x.begin < y.begin; });
  d.seconds.call = since(t);

  if (o.small_variants) {
    t = clock::now();
    std::vector<variant::ReadView> views;
    views.reserve(2 * pairs.size());
    for (const auto& p : pairs)
      for (int m = 0; m < 2; m++) views.push_back({std::span<const sequence::Code>(p.mates[m].bases), std::string_view(p.mates[m].quality)});
    d.small = variant::call_variants(pool, contigs, mates, views, o.caller);
    d.seconds.snv = since(t);
    t = clock::now();
    d.sites = allele_sites(d.small);
    d.runs = runs_of_homozygosity(d.sites, o.roh);
    d.allelic = allelic_segments(d.sites, d.runs, d.segments);
    d.seconds.allelic = since(t);
  }
  return d;
}

}  // namespace genomics::structural
