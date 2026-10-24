// Checks of assembly/debruijn.hpp: a genome without repeats of k − 1 bases is one unitig; a
// genome a R b R c with a repeat R longer than k makes the four unitigs and four links the
// theory predicts, R's at coverage 2; on error-free reads of a 200 kb genome with repeats
// longer and shorter than k, every unitig is a substring of the genome or of its reverse
// complement, the unitigs hold every solid k-mer exactly once, every unitig is maximal, every
// link is an overlap of k − 1 bases, coverage is the mean count; the graph is the same on 1
// and 8 threads and in canonical form; compaction after removing unitigs gives the graph built
// from the remaining k-mers; and GFA has a segment per unitig and a line per link.

#include "truth.hpp"

#include <genomics/assembly/debruijn.hpp>
#include <genomics/kmers/counter.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

using namespace genomics;
using assembly::Graph;
using sequence::Codes;
using sequence::Kmer;

namespace {

bool same_graph(const Graph& a, const Graph& b) {
  if (a.k != b.k || a.size() != b.size()) return false;
  for (std::size_t u = 0; u < a.size(); u++)
    if (a.unitigs[u].sequence != b.unitigs[u].sequence || a.unitigs[u].count_sum != b.unitigs[u].count_sum) return false;
  for (std::size_t o = 0; o < a.next.size(); o++) {
    auto x = a.next[o].to, y = b.next[o].to;
    std::sort(x.begin(), x.begin() + a.next[o].size);
    std::sort(y.begin(), y.begin() + b.next[o].size);
    if (a.next[o].size != b.next[o].size || x != y) return false;
  }
  return true;
}

// The canonical k-mers of the unitigs, each with how many times the unitigs hold it.
std::unordered_map<Kmer, int> unitig_kmers(const Graph& g) {
  std::unordered_map<Kmer, int> m;
  for (const auto& u : g.unitigs) sequence::for_each_kmer(u.sequence, g.k, [&](std::size_t, Kmer f, Kmer r) { m[std::min(f, r)]++; });
  return m;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::scheduler::ThreadPool one(1), eight(8);
  const int k = 31;

  // A genome with no repeated (k − 1)-mer: one unitig, the genome on one of its strands.
  {
    const Codes genome = simulate::random_sequence(5000, 0.5, 1);
    const auto counts = kmers::count_kmers({genome}, eight);
    const Graph g = assembly::build_graph(counts.sorted(1, eight), k, eight);
    s.require("de Bruijn: a random 5 kb genome is one unitig, spelled as the genome or its reverse complement, with no links",
              g.size() == 1 && (g.unitigs[0].sequence == genome || g.unitigs[0].sequence == sequence::reverse_complement(genome)) && g.links().empty());
  }

  // a R b R c: the repeat R (100 bases) is a unitig of coverage 2 between three others.
  {
    const Codes a = simulate::random_sequence(300, 0.5, 2), R = simulate::random_sequence(100, 0.5, 3), b = simulate::random_sequence(300, 0.5, 4),
                c = simulate::random_sequence(300, 0.5, 5);
    Codes x = a;
    for (const Codes* part : {&R, &b, &R, &c}) x.insert(x.end(), part->begin(), part->end());
    const auto counts = kmers::count_kmers({x}, eight);
    const Graph g = assembly::build_graph(counts.sorted(1, eight), k, eight);
    std::size_t repeat = assembly::NONE;
    for (std::uint32_t u = 0; u < g.size(); u++)
      if (g.unitigs[u].sequence == R || g.unitigs[u].sequence == sequence::reverse_complement(R)) repeat = u;
    s.require("de Bruijn: a R b R c makes four unitigs and four links, R one of them", g.size() == 4 && g.links().size() == 4 && repeat != assembly::NONE,
              std::to_string(g.size()) + " unitigs, " + std::to_string(g.links().size()) + " links");
    s.require("de Bruijn: the repeat's unitig has coverage 2, two arcs in and two out", repeat != assembly::NONE && g.coverage(static_cast<std::uint32_t>(repeat)) == 2 &&
                                                                                       g.in_degree(assembly::oriented(static_cast<std::uint32_t>(repeat), false)) == 2 &&
                                                                                       g.out_degree(assembly::oriented(static_cast<std::uint32_t>(repeat), false)) == 2);
  }

  // Error-free reads of a 200 kb genome with repeats longer than k (400 to 1000 bases) and
  // shorter (20 bases).
  simulate::GenomeSpec spec;
  spec.length = 200'000;
  spec.repeats = {6, 3, 700, 0.0};
  auto genome = simulate::make_genome(spec, 10);
  {
    calx::random::Stream r(11, 0x53484f52, 0, 0);
    for (int e = 0; e < 5; e++) {
      const Codes unit = simulate::random_sequence(20, 0.5, 12, static_cast<std::uint32_t>(e));
      for (int c = 0; c < 8; c++) {
        const std::size_t at = r.below(spec.length - 20);
        std::copy(unit.begin(), unit.end(), genome.sequence.begin() + static_cast<std::ptrdiff_t>(at));
      }
    }
  }
  const auto source = simulate::Source::haploid({genome.sequence});
  const auto reads = truth::reads(source, 30 * spec.length / 150, {150, 0, 0, 0}, 13, eight);
  const auto counts = kmers::count_kmers(reads, eight);
  const auto solid = counts.sorted(1, eight);
  const Graph g = assembly::build_graph(solid, k, eight);
  const truth::GenomeIndex index(genome.sequence, k);
  {
    std::size_t outside = 0;
    for (const auto& u : g.unitigs) outside += !index.contains(u.sequence);
    s.require("de Bruijn: every unitig of error-free reads is a substring of the genome or of its reverse complement", outside == 0,
              std::to_string(g.size()) + " unitigs, " + std::to_string(outside) + " outside");
    const auto held = unitig_kmers(g);
    bool once = held.size() == solid.size();
    for (const auto& e : solid) {
      const auto it = held.find(e.kmer);
      once = once && it != held.end() && it->second == 1;
    }
    s.require("de Bruijn: the unitigs hold every solid k-mer exactly once", once, std::to_string(solid.size()) + " k-mers");
    bool maximal = true;
    for (std::uint32_t o = 0; o < 2 * g.size(); o++) {
      if (g.out_degree(o) != 1) continue;
      const std::uint32_t p = g.successors(o)[0];
      // A unitig may continue into one that does not branch only when that would close a cycle on itself.
      maximal = maximal && (g.in_degree(p) != 1 || assembly::element(p) == assembly::element(o));
    }
    s.require("de Bruijn: every unitig is maximal: each end branches, or its neighbour has another arc in", maximal);
    bool overlaps = true;
    for (std::uint32_t o = 0; o < 2 * g.size(); o++)
      for (auto p : g.successors(o)) {
        const Codes x = g.sequence(o), y = g.sequence(p);
        overlaps = overlaps && std::equal(x.end() - (k - 1), x.end(), y.begin());
      }
    s.require("de Bruijn: every link joins the last k − 1 bases of one oriented unitig to the first k − 1 of the next", overlaps, std::to_string(g.links().size()) + " links");
    bool coverage = true, canonical = true;
    Kmer previous = 0;
    for (std::uint32_t u = 0; u < g.size(); u++) {
      std::uint64_t sum = 0;
      Kmer least = ~Kmer{0};
      sequence::for_each_kmer(g.unitigs[u].sequence, k, [&](std::size_t, Kmer f, Kmer r) {
        sum += counts.count(f);
        least = std::min(least, std::min(f, r));
      });
      coverage = coverage && sum == g.unitigs[u].count_sum;
      const Codes& q = g.unitigs[u].sequence;
      canonical = canonical && sequence::encode_kmer(q.data(), k) < sequence::reverse_complement(sequence::encode_kmer(q.data() + q.size() - k, k), k);
      if (u > 0) canonical = canonical && previous < least;
      previous = least;
    }
    s.require("de Bruijn: a unitig's coverage is the mean count of its k-mers", coverage);
    s.require("de Bruijn: canonical form: each unitig spelled from its smaller end, numbered by its least k-mer", canonical);
    s.require("de Bruijn: the same graph on 1 and 8 threads", same_graph(g, assembly::build_graph(solid, k, one)));
    std::vector<std::vector<std::uint32_t>> merged;
    s.require("de Bruijn: compaction with nothing removed is the identity", same_graph(g, assembly::compact(g, std::vector<char>(g.size(), 0), &merged)));
    // Remove every seventh unitig: the compaction equals the graph of the k-mers that remain.
    std::vector<char> removed(g.size(), 0);
    for (std::size_t u = 0; u < g.size(); u += 7) removed[u] = 1;
    const Graph h = assembly::compact(g, removed, &merged);
    std::vector<kmers::KmerCount> rest;
    for (std::uint32_t u = 0; u < g.size(); u++)
      if (!removed[u]) sequence::for_each_kmer(g.unitigs[u].sequence, k, [&](std::size_t, Kmer f, Kmer r) { rest.push_back({std::min(f, r), counts.count(f)}); });
    std::sort(rest.begin(), rest.end(), [](const auto& x, const auto& y) { return x.kmer < y.kmer; });
    s.require("de Bruijn: compaction after removing every seventh unitig is the graph built from the k-mers left", same_graph(h, assembly::build_graph(rest, k, eight)),
              std::to_string(g.size()) + " unitigs, then " + std::to_string(h.size()));
    std::size_t parts = 0;
    for (const auto& m : merged) parts += m.size();
    s.require("de Bruijn: the merged unitigs are made of the unitigs kept, each once", parts == g.size() - static_cast<std::size_t>(std::count(removed.begin(), removed.end(), 1)));
    const std::string gfa = g.gfa();
    s.require("de Bruijn: GFA has a segment per unitig and a link line per bidirected edge",
              static_cast<std::size_t>(std::count(gfa.begin(), gfa.end(), '\n')) == 1 + g.size() + g.links().size() && gfa.rfind("H\tVN:Z:1.0", 0) == 0);
  }
  return s.done();
}
