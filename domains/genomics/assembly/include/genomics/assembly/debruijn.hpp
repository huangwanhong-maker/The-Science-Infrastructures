// The de Bruijn graph of a set of canonical k-mers as a bidirected graph, and its compaction
// into unitigs, the maximal non-branching paths: the compacted graph, an object an
// application can draw (unitig sequences, lengths, coverage, links with their orientations).
//
// - The graph of a double strand (Medvedev, Georgiou, Myers and Brudno 2007, "Computability of
//   models for sequence assembly", WABI, LNBI 4645, 289–301, §2; Compeau, Pevzner and Tesler
//   2011, "How to apply de Bruijn graphs to genome assembly", Nature Biotechnology 29,
//   987–991): a node is a canonical k-mer x, standing for x and its reverse complement x̄. It
//   is read in two orientations, (x, +) spelling x and (x, −) spelling x̄; an oriented node o
//   has its reverse o⁻ ((x, ±) ↦ (x, ∓)). There is an arc o → o′ when the last k − 1 bases o
//   spells are the first k − 1 bases o′ spells, and with it the arc o′⁻ → o⁻ that reads the
//   same overlap on the other strand: the arcs come in pairs, and a pair is one bidirected
//   edge, joining a side of each node. k is odd, so that no k-mer is its own reverse
//   complement and the two orientations of a node always differ. The graph is node-centric:
//   every (k − 1)-overlap of two k-mers of the set is an edge, whether or not a read spans it.
//   The structure is a skew-symmetric graph (Goldberg and Karzanov 1996, "Path problems in
//   skew-symmetric graphs", Combinatorica 16, 353–382): a digraph on the oriented nodes with
//   the involution o ↦ o⁻ reversing its arcs.
// - Unitigs (Chikhi, Limasset and Medvedev 2016, "Compacting de Bruijn graphs from sequencing
//   data quickly and in low memory", Bioinformatics 32, i201–i208, BCALM 2, §2, for the
//   definition): a unitig is a maximal path o₁ → … → o_m of distinct nodes in which every arc
//   oᵢ → oᵢ₊₁ is the only arc out of oᵢ and the only arc into oᵢ₊₁. The unitigs partition the
//   nodes, their spelling has m + k − 1 bases, and the arcs between them join only their ends,
//   so the compacted graph is a bidirected graph on the unitigs with the same paths spelled.
// - The algorithm is the direct one: each node's two sides find their neighbours by four hash
//   lookups each (in parallel over the nodes; deterministic, each node writes its own slots);
//   then the non-branching paths, walked in parallel from their ends (detail::chains); then the
//   arcs between unitig ends. O(n) expected time and O(n) memory for n nodes: about 70 bytes a
//   node during the build (the index, 24 to 48 bytes a node at a load between 1/4 and 1/2; the
//   degrees and successors of the sides, 10; the chains, about 12), 357 MB for the 4.9 million
//   nodes of a 5 Mb genome, besides the 16 bytes a node of the sorted k-mers given.
// - Canonical form: each unitig is spelled in the orientation whose first k-mer is the smaller
//   integer (as sequence/kmer.hpp orders k-mers), and unitigs are numbered in increasing order
//   of their least canonical k-mer, so that the compacted graph is a function of the k-mer set
//   alone, however it was reached (compact() after cleaning gives the graph that building from
//   the remaining k-mers gives). Coverage is the mean count of a unitig's k-mers, BCALM's km.
// - The links are written as GFA 1 writes them (the GFA Format Specification Working Group
//   2016): L from ± to ± with overlap (k − 1)M; a bidirected edge is written once.

#pragma once

#include <genomics/kmers/counter.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/kmer.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace genomics::assembly {

using kmers::KmerCount;
using kmers::mix64;
using sequence::Code;
using sequence::Codes;
using sequence::Kmer;

inline constexpr std::uint32_t NONE = UINT32_MAX;

// An oriented element (a node or a unitig): 2i + 1 for its reverse complement, 2i for itself.
constexpr std::uint32_t oriented(std::uint32_t i, bool reverse) { return 2 * i + (reverse ? 1u : 0u); }
constexpr std::uint32_t element(std::uint32_t o) { return o >> 1; }
constexpr bool is_reverse(std::uint32_t o) { return o & 1; }
constexpr std::uint32_t flip(std::uint32_t o) { return o ^ 1u; }

inline void check_odd_k(int k, const char* who) {
  if (k < 3 || k > 31 || k % 2 == 0) throw std::invalid_argument(std::string("genomics::assembly::") + who + ": k odd, 3 <= k <= 31");
}

// The successors of one side: at most four, one for each next base.
struct Adjacency {
  std::array<std::uint32_t, 4> to{NONE, NONE, NONE, NONE};
  std::uint8_t size = 0;
  std::span<const std::uint32_t> view() const { return {to.data(), size}; }
  void push(std::uint32_t o) {
    if (size < 4) to[size++] = o;
  }
};

struct Unitig {
  Codes sequence;              // in its canonical orientation
  std::uint64_t count_sum = 0; // the counts of its k-mers, summed
  std::size_t length() const { return sequence.size(); }
  std::size_t kmers(int k) const { return sequence.size() - static_cast<std::size_t>(k) + 1; }
  double coverage(int k) const { return static_cast<double>(count_sum) / static_cast<double>(kmers(k)); }
};

// A bidirected edge as GFA writes it: the end of `from` (read as its reverse complement when
// from_reverse) overlaps the start of `to` (likewise) by k − 1 bases.
struct Link {
  std::uint32_t from = 0, to = 0;
  bool from_reverse = false, to_reverse = false;
  friend bool operator==(const Link&, const Link&) = default;
};

// The compacted de Bruijn graph: unitigs, and for each oriented unitig its successors.
class Graph {
 public:
  int k = 0;
  std::vector<Unitig> unitigs;
  std::vector<Adjacency> next;  // per oriented unitig

  std::size_t size() const { return unitigs.size(); }
  std::span<const std::uint32_t> successors(std::uint32_t o) const { return next[o].view(); }
  // The predecessors of o: the reverses of the successors of o⁻.
  Adjacency predecessors(std::uint32_t o) const {
    Adjacency a;
    for (auto p : next[flip(o)].view()) a.push(flip(p));
    return a;
  }
  std::size_t out_degree(std::uint32_t o) const { return next[o].size; }
  std::size_t in_degree(std::uint32_t o) const { return next[flip(o)].size; }
  std::size_t length(std::uint32_t u) const { return unitigs[u].length(); }
  double coverage(std::uint32_t u) const { return unitigs[u].coverage(k); }
  std::size_t total_length() const {
    std::size_t t = 0;
    for (const auto& u : unitigs) t += u.length();
    return t;
  }
  // The sequence an oriented unitig spells.
  Codes sequence(std::uint32_t o) const {
    const Codes& s = unitigs[element(o)].sequence;
    return is_reverse(o) ? sequence::reverse_complement(s) : s;
  }

  // Every bidirected edge once: o → o′ is written when (o, o′) ≤ (o′⁻, o⁻).
  std::vector<Link> links() const {
    std::vector<Link> out;
    for (std::uint32_t o = 0; o < next.size(); o++)
      for (auto p : next[o].view())
        if (std::pair(o, p) <= std::pair(flip(p), flip(o))) out.push_back({element(o), element(p), is_reverse(o), is_reverse(p)});
    return out;
  }

  // GFA 1: a segment per unitig with its length, k-mer count and coverage, a link per edge.
  std::string gfa() const {
    std::string s = "H\tVN:Z:1.0\n";
    for (std::size_t u = 0; u < unitigs.size(); u++) {
      s += "S\t" + std::to_string(u) + "\t" + sequence::decode_dna(unitigs[u].sequence) + "\tLN:i:" + std::to_string(unitigs[u].length()) +
           "\tKC:i:" + std::to_string(unitigs[u].count_sum) + "\tkm:f:" + std::to_string(unitigs[u].coverage(k)) + "\n";
    }
    for (const auto& l : links())
      s += "L\t" + std::to_string(l.from) + "\t" + (l.from_reverse ? "-" : "+") + "\t" + std::to_string(l.to) + "\t" + (l.to_reverse ? "-" : "+") + "\t" +
           std::to_string(k - 1) + "M\n";
    return s;
  }
};

namespace detail {

// The first k-mer of a sequence, as an integer.
inline Kmer first_kmer(const Codes& s, int k) { return sequence::encode_kmer(s.data(), k); }
inline Kmer last_kmer(const Codes& s, int k) { return sequence::encode_kmer(s.data() + s.size() - static_cast<std::size_t>(k), k); }

// Whether a spelling should be reversed into its canonical orientation: its first k-mer is
// larger than the first k-mer of its reverse complement.
inline bool spelled_backward(const Codes& s, int k) { return first_kmer(s, k) > sequence::reverse_complement(last_kmer(s, k), k); }

// The least canonical k-mer of a sequence (all bases below 4).
inline Kmer least_kmer(const Codes& s, int k) {
  Kmer m = ~Kmer{0};
  sequence::for_each_kmer(s, k, [&](std::size_t, Kmer f, Kmer r) { m = std::min(m, std::min(f, r)); });
  return m;
}

// The maximal non-branching chains of n elements, each a sequence of oriented elements.
// degree(o) is the number of successors of oriented o, unique(o) the successor when there is
// one; alive(i) excludes elements. The chains come in order of their least element, and
// chain_of receives each element's chain.
//
// An arc o → p extends a chain when it is the only arc out of o and the only one into p, and
// joins two different elements. These arcs form paths and cycles on the oriented elements,
// each path P paired with its reverse P⁻ (the arcs come in pairs, and no path is its own
// reverse, since its middle arc would join o to o⁻). A path begins at a head, an oriented
// element with no extending arc into it, and every path has two heads, its first element and
// the reverse of its last. The paths are walked from all heads in parallel, each kept from the
// head that is the smaller of its two; the elements left over lie on cycles, which are walked
// one at a time from their least element in the + orientation. The result is the same on any
// number of threads. O(n) work, each path walked twice.
template <class Alive, class Degree, class Unique>
std::vector<std::vector<std::uint32_t>> chains(std::size_t n, Alive&& alive, Degree&& degree, Unique&& unique, std::vector<std::uint32_t>& chain_of,
                                               calx::scheduler::ThreadPool& pool) {
  std::vector<std::uint8_t> extends(2 * n, 0);
  calx::scheduler::parallel_for(pool, n, 4096, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t i = b; i < e; i++) {
      if (!alive(static_cast<std::uint32_t>(i))) continue;
      for (std::uint32_t r = 0; r < 2; r++) {
        const std::uint32_t o = oriented(static_cast<std::uint32_t>(i), r);
        if (degree(o) != 1) continue;
        const std::uint32_t p = unique(o);
        extends[o] = degree(flip(p)) == 1 && element(p) != i;
      }
    }
  });
  struct Found {
    std::uint32_t least = 0;
    std::vector<std::uint32_t> chain;
  };
  constexpr std::size_t GRAIN = 4096;
  const std::size_t parts = (n + GRAIN - 1) / GRAIN;
  std::vector<std::vector<Found>> found(parts);
  calx::scheduler::parallel_for(pool, parts, 1, [&](std::size_t b, std::size_t e, std::size_t) {
    std::vector<std::uint32_t> walk;
    for (std::size_t q = b; q < e; q++)
      for (std::size_t i = q * GRAIN; i < std::min(n, (q + 1) * GRAIN); i++) {
        if (!alive(static_cast<std::uint32_t>(i))) continue;
        for (std::uint32_t r = 0; r < 2; r++) {
          const std::uint32_t h = oriented(static_cast<std::uint32_t>(i), r);
          if (extends[flip(h)]) continue;  // an arc extends into h: not a head
          walk.assign(1, h);
          std::uint32_t least = static_cast<std::uint32_t>(i);
          for (std::uint32_t o = h; extends[o];) {
            o = unique(o);
            walk.push_back(o);
            least = std::min(least, element(o));
          }
          if (h < flip(walk.back())) found[q].push_back({least, walk});
        }
      }
  });
  std::vector<Found> all;
  for (auto& f : found)
    for (auto& x : f) all.push_back(std::move(x));
  chain_of.assign(n, NONE);
  calx::scheduler::parallel_for(pool, all.size(), 256, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t c = b; c < e; c++)
      for (auto o : all[c].chain) chain_of[element(o)] = static_cast<std::uint32_t>(c);
  });
  // The cycles, each from its least element.
  for (std::uint32_t i = 0; i < n; i++) {
    if (!alive(i) || chain_of[i] != NONE) continue;
    Found f{i, {oriented(i, false)}};
    chain_of[i] = static_cast<std::uint32_t>(all.size());
    for (std::uint32_t o = unique(oriented(i, false)); o != oriented(i, false); o = unique(o)) {
      f.chain.push_back(o);
      chain_of[element(o)] = static_cast<std::uint32_t>(all.size());
    }
    all.push_back(std::move(f));
  }
  std::sort(all.begin(), all.end(), [](const Found& x, const Found& y) { return x.least < y.least; });
  std::vector<std::vector<std::uint32_t>> out(all.size());
  calx::scheduler::parallel_for(pool, all.size(), 256, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t c = b; c < e; c++) {
      for (auto o : all[c].chain) chain_of[element(o)] = static_cast<std::uint32_t>(c);
      out[c] = std::move(all[c].chain);
    }
  });
  return out;
}

// A hash index from canonical k-mers to node numbers, partitioned by the high bits of mix64
// and built one partition per task.
class NodeIndex {
 public:
  static constexpr unsigned BITS = 10;
  NodeIndex(const std::vector<KmerCount>& nodes, calx::scheduler::ThreadPool& pool) : tables_(std::size_t{1} << BITS) {
    const std::size_t P = tables_.size();
    std::vector<std::vector<std::uint32_t>> members(P);
    for (std::uint32_t i = 0; i < nodes.size(); i++) members[mix64(nodes[i].kmer) >> (64 - BITS)].push_back(i);
    calx::scheduler::parallel_for(pool, P, 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t p = b; p < e; p++) {
        tables_[p] = kmers::KmerMap(2 * members[p].size());
        for (auto i : members[p]) tables_[p].at(nodes[i].kmer, mix64(nodes[i].kmer)) = i;
      }
    });
  }
  std::uint32_t find(Kmer x) const {
    const std::uint64_t h = mix64(x);
    const std::uint32_t* v = tables_[h >> (64 - BITS)].find(x, h);
    return v ? *v : NONE;
  }
  void prefetch(Kmer x) const {
    const std::uint64_t h = mix64(x);
    tables_[h >> (64 - BITS)].prefetch(h);
  }
  std::size_t memory_bytes() const {
    std::size_t b = 0;
    for (const auto& t : tables_) b += t.memory_bytes();
    return b;
  }

 private:
  std::vector<kmers::KmerMap> tables_;
};

}  // namespace detail

// The successors of the oriented k-mer o = (node, reverse) among the nodes: up to four
// oriented nodes, in the order of the next base A, C, G, T.
inline Adjacency kmer_successors(const std::vector<KmerCount>& nodes, const detail::NodeIndex& index, int k, std::uint32_t o) {
  const Kmer x = nodes[element(o)].kmer;
  const Kmer v = is_reverse(o) ? sequence::reverse_complement(x, k) : x, rv = is_reverse(o) ? x : sequence::reverse_complement(x, k);
  const Kmer mask = sequence::kmer_mask(k);
  const int top = 2 * (k - 1);
  std::array<Kmer, 4> w, c;
  for (unsigned b = 0; b < 4; b++) {
    w[b] = ((v << 2) | b) & mask;
    const Kmer rw = (rv >> 2) | (static_cast<Kmer>(3 - b) << top);
    c[b] = std::min(w[b], rw);
    index.prefetch(c[b]);
  }
  Adjacency a;
  for (unsigned b = 0; b < 4; b++) {
    const std::uint32_t j = index.find(c[b]);
    if (j != NONE) a.push(oriented(j, w[b] != c[b]));
  }
  return a;
}

struct BuildReport {
  std::size_t nodes = 0, index_bytes = 0, peak_bytes = 0;
};

// The compacted de Bruijn graph of the k-mers `nodes` (canonical, distinct, sorted by k-mer,
// as KmerCounter::sorted gives them; their counts give the coverage).
inline Graph build_graph(const std::vector<KmerCount>& nodes, int k, calx::scheduler::ThreadPool& pool, BuildReport* report = nullptr) {
  check_odd_k(k, "build_graph");
  for (std::size_t i = 1; i < nodes.size(); i++)
    if (!(nodes[i - 1].kmer < nodes[i].kmer)) throw std::invalid_argument("genomics::assembly::build_graph: nodes distinct and sorted by k-mer");
  if (nodes.size() >= NONE / 2) throw std::invalid_argument("genomics::assembly::build_graph: fewer than 2^31 nodes");
  const std::size_t n = nodes.size();
  const detail::NodeIndex index(nodes, pool);
  // Each side's degree and, when it is one, its successor.
  std::vector<std::uint8_t> degree(2 * n);
  std::vector<std::uint32_t> unique(2 * n, NONE);
  calx::scheduler::parallel_for(pool, n, 4096, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t i = b; i < e; i++)
      for (std::uint32_t r = 0; r < 2; r++) {
        const std::uint32_t o = oriented(static_cast<std::uint32_t>(i), r);
        const Adjacency a = kmer_successors(nodes, index, k, o);
        degree[o] = a.size;
        if (a.size == 1) unique[o] = a.to[0];
      }
  });
  std::vector<std::uint32_t> unitig_of;
  auto chains = detail::chains(
      n, [](std::uint32_t) { return true; }, [&](std::uint32_t o) { return degree[o]; }, [&](std::uint32_t o) { return unique[o]; }, unitig_of, pool);
  // Spell each chain in its canonical orientation; chains come in order of their least node,
  // which is their least k-mer, so this numbering is already the canonical one.
  Graph g;
  g.k = k;
  g.unitigs.resize(chains.size());
  std::vector<std::uint32_t> first(chains.size()), last(chains.size());  // oriented nodes, in the unitig's orientation
  calx::scheduler::parallel_for(pool, chains.size(), 256, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t u = b; u < e; u++) {
      auto& c = chains[u];
      const auto spell = [&](std::uint32_t o) {
        const Kmer x = nodes[element(o)].kmer;
        return is_reverse(o) ? sequence::reverse_complement(x, k) : x;
      };
      if (spell(c.front()) > sequence::reverse_complement(spell(c.back()), k)) {
        std::reverse(c.begin(), c.end());
        for (auto& o : c) o = flip(o);
      }
      Unitig& t = g.unitigs[u];
      t.sequence.resize(c.size() + static_cast<std::size_t>(k) - 1);
      const Kmer x0 = spell(c.front());
      for (int i = 0; i < k; i++) t.sequence[static_cast<std::size_t>(i)] = static_cast<Code>((x0 >> (2 * (k - 1 - i))) & 3);
      for (std::size_t j = 1; j < c.size(); j++) t.sequence[j + static_cast<std::size_t>(k) - 1] = static_cast<Code>(spell(c[j]) & 3);
      for (auto o : c) t.count_sum += nodes[element(o)].count;
      first[u] = c.front();
      last[u] = c.back();
    }
  });
  // The arcs between unitig ends.
  g.next.resize(2 * g.unitigs.size());
  const auto resolve = [&](std::uint32_t q) {
    const std::uint32_t w = unitig_of[element(q)];
    if (first[w] == q) return oriented(w, false);
    if (last[w] == flip(q)) return oriented(w, true);
    throw std::logic_error("genomics::assembly::build_graph: an arc into the middle of a unitig");
  };
  calx::scheduler::parallel_for(pool, g.unitigs.size(), 256, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t u = b; u < e; u++)
      for (std::uint32_t r = 0; r < 2; r++) {
        // Out of (u, +) is out of its last node; out of (u, −) is out of its first node reversed.
        const std::uint32_t end = r ? flip(first[u]) : last[u];
        Adjacency a;
        if (degree[end] == 1) a.push(unique[end]);
        else if (degree[end] > 1) a = kmer_successors(nodes, index, k, end);
        Adjacency& out = g.next[oriented(static_cast<std::uint32_t>(u), r)];
        for (auto q : a.view()) out.push(resolve(q));
      }
  });
  if (report) {
    report->nodes = n;
    report->index_bytes = index.memory_bytes();
    // The index, the degrees and successors of the sides, the chain of each node and the chains.
    report->peak_bytes = index.memory_bytes() + n * (2 * sizeof(std::uint8_t) + 2 * sizeof(std::uint32_t) + 2 * sizeof(std::uint32_t)) + g.total_length();
  }
  return g;
}

// The graph with the unitigs `removed` taken out and the paths that no longer branch merged:
// the compacted graph of the remaining k-mers, in canonical form. merged_from, when given,
// receives for each new unitig the oriented old unitigs it is made of, in its orientation.
inline Graph compact(const Graph& g, const std::vector<char>& removed, std::vector<std::vector<std::uint32_t>>* merged_from = nullptr) {
  const std::size_t n = g.size();
  if (removed.size() != n) throw std::invalid_argument("genomics::assembly::compact: one flag per unitig");
  // The successors among the unitigs kept.
  std::vector<Adjacency> next(2 * n);
  for (std::uint32_t o = 0; o < 2 * n; o++) {
    if (removed[element(o)]) continue;
    for (auto p : g.next[o].view())
      if (!removed[element(p)]) next[o].push(p);
  }
  // A graph of unitigs is small: it is compacted on the calling thread.
  calx::scheduler::ThreadPool serial(1);
  std::vector<std::uint32_t> chain_of;
  auto chains = detail::chains(
      n, [&](std::uint32_t i) { return !removed[i]; }, [&](std::uint32_t o) { return next[o].size; }, [&](std::uint32_t o) { return next[o].to[0]; }, chain_of, serial);
  const int k = g.k;
  struct Built {
    Unitig unitig;
    std::vector<std::uint32_t> parts;
    Kmer least = 0;
  };
  std::vector<Built> built(chains.size());
  for (std::size_t c = 0; c < chains.size(); c++) {
    auto& parts = chains[c];
    Codes s = g.sequence(parts.front());
    for (std::size_t j = 1; j < parts.size(); j++) {
      const Codes t = g.sequence(parts[j]);
      s.insert(s.end(), t.begin() + (k - 1), t.end());
    }
    if (detail::spelled_backward(s, k)) {
      s = sequence::reverse_complement(s);
      std::reverse(parts.begin(), parts.end());
      for (auto& o : parts) o = flip(o);
    }
    Built& b = built[c];
    for (auto o : parts) b.unitig.count_sum += g.unitigs[element(o)].count_sum;
    b.least = detail::least_kmer(s, k);
    b.unitig.sequence = std::move(s);
    b.parts = std::move(parts);
  }
  // Number the new unitigs by their least k-mer.
  std::vector<std::uint32_t> order(built.size());
  for (std::uint32_t i = 0; i < order.size(); i++) order[i] = i;
  std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) { return built[a].least < built[b].least; });
  std::vector<std::uint32_t> rank(built.size());
  for (std::uint32_t r = 0; r < order.size(); r++) rank[order[r]] = r;
  // The new oriented unitig that each old oriented unitig begins, where it begins one.
  std::vector<std::uint32_t> begins(2 * n, NONE);
  for (std::uint32_t c = 0; c < built.size(); c++) {
    begins[built[c].parts.front()] = oriented(rank[c], false);
    begins[flip(built[c].parts.back())] = oriented(rank[c], true);
  }
  Graph h;
  h.k = k;
  h.unitigs.resize(built.size());
  h.next.resize(2 * built.size());
  if (merged_from) merged_from->assign(built.size(), {});
  for (std::uint32_t c = 0; c < built.size(); c++) {
    const std::uint32_t u = rank[c];
    for (std::uint32_t r = 0; r < 2; r++) {
      const std::uint32_t old_end = r ? flip(built[c].parts.front()) : built[c].parts.back();
      for (auto p : next[old_end].view()) {
        if (begins[p] == NONE) throw std::logic_error("genomics::assembly::compact: an arc into the middle of a unitig");
        h.next[oriented(u, r)].push(begins[p]);
      }
    }
    if (merged_from) (*merged_from)[u] = built[c].parts;
    h.unitigs[u] = std::move(built[c].unitig);
  }
  return h;
}

}  // namespace genomics::assembly
