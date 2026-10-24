// Assembly as an Eulerian path (learning mode): the k-mer composition of one strand, its
// directed de Bruijn graph, the strings it spells, how many there are, and the de Bruijn
// sequences.
//
// - The k-mer composition of a string s is the multiset of its |s| − k + 1 words of length k.
//   Its de Bruijn graph has a node for each distinct (k − 1)-mer of the composition and an
//   edge, from prefix to suffix, for each k-mer, with its multiplicity (Pevzner 1989,
//   "l-tuple DNA sequencing: computer analysis", J. Biomol. Struct. Dyn. 7, 63–73; Idury and
//   Waterman 1995, J. Comput. Biol. 2, 291–306; Compeau, Pevzner and Tesler 2011). It is the
//   line-graph construction turned around: the de Bruijn graph B(σ, n) of all words is the line
//   digraph of B(σ, n − 1), and the graph of a composition is a subgraph of B(4, k − 1).
// - Strings with the composition are the Eulerian trails, read by their edge labels: a trail
//   spells its first k-mer and then the last base of each next one. Parallel edges carry the
//   same k-mer (a prefix and a suffix of k − 1 bases determine it), so permuting the copies of
//   a repeated k-mer permutes trails without changing the string, and the map from trails to
//   strings is ∏_x m(x)!-to-one, m(x) the multiplicity of the k-mer x. Hence
//
//     #{t : composition(t) = composition(s)} = trails(G) / ∏_x m(x)!,
//
//   with trails(G) counted as edge sequences by the BEST theorem (morphe/combinatorics/
//   eulerian.hpp): ec(G + (t → s)) when the graph has a trail from s ≠ t, and m · ec(G) when it
//   is balanced (s begins and ends with the same (k − 1)-mer, and every rotation of a circuit
//   is a string). Exact, in morphe's integers.
// - Reconstruction by Hierholzer's algorithm (1873), the unused edge of least index first:
//   it returns the string when there is exactly one, and one of them when there are several.
//   The ambiguity is the subject: every repeat of k − 1 bases or more is a vertex the path
//   passes twice, where it may leave by either exit.
// - De Bruijn sequences of order n over σ letters, the cyclic words in which every word of
//   length n occurs once (de Bruijn 1946, "A combinatorial problem", Proc. Ned. Akad. Wet. 49,
//   758–764; Flye Sainte-Marie 1894 for σ = 2), are the Eulerian circuits of B(σ, n − 1); the
//   BEST theorem counts them as (σ!)^{σ^{n−1}} / σ^n (van Aardenne-Ehrenfest and de Bruijn
//   1951), 2^{2^{n−1} − n} for σ = 2.
//
// Cost: the graph in O(m log m) for m k-mers; a trail in O(m); the count is a determinant of
// order n − 1 for n nodes, O(n³) operations on exact integers.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/kmer.hpp>

#include <morphe/arithmetic/integer.hpp>
#include <morphe/combinatorics/eulerian.hpp>
#include <morphe/combinatorics/graph.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

namespace genomics::assembly {

using morphe::arithmetic::Integer;
using morphe::combinatorics::Multigraph;
using sequence::Code;
using sequence::Codes;
using sequence::Kmer;

// The k-mers of one strand of s, with multiplicity, sorted; s free of N.
struct Composition {
  int k = 0;
  std::vector<Kmer> kmers;
  friend bool operator==(const Composition&, const Composition&) = default;
};

inline Composition composition(const Codes& s, int k) {
  sequence::check_k(k, "assembly::composition");
  for (Code c : s)
    if (c > 3) throw std::invalid_argument("genomics::assembly::composition: a sequence of A, C, G and T");
  Composition c;
  c.k = k;
  sequence::for_each_kmer(s, k, [&](std::size_t, Kmer f, Kmer) { c.kmers.push_back(f); });
  std::sort(c.kmers.begin(), c.kmers.end());
  return c;
}

// The directed de Bruijn graph of a composition: node i the (k − 1)-mer nodes[i] (sorted),
// edge e from the prefix of edge_kmers[e] to its suffix, the edges in the order of the
// composition.
struct KmerGraph {
  int k = 0;
  std::vector<Kmer> nodes, edge_kmers;
  Multigraph graph;
  std::size_t node(Kmer x) const {
    const auto it = std::lower_bound(nodes.begin(), nodes.end(), x);
    return it != nodes.end() && *it == x ? static_cast<std::size_t>(it - nodes.begin()) : SIZE_MAX;
  }
};

inline KmerGraph kmer_graph(const Composition& c) {
  if (c.k < 2) throw std::invalid_argument("genomics::assembly::kmer_graph: k >= 2");
  KmerGraph g;
  g.k = c.k;
  const Kmer low = sequence::kmer_mask(c.k - 1);
  for (Kmer x : c.kmers) {
    g.nodes.push_back(x >> 2);
    g.nodes.push_back(x & low);
  }
  std::sort(g.nodes.begin(), g.nodes.end());
  g.nodes.erase(std::unique(g.nodes.begin(), g.nodes.end()), g.nodes.end());
  g.graph.vertices = g.nodes.size();
  g.edge_kmers = c.kmers;
  for (Kmer x : c.kmers) g.graph.edges.push_back({g.node(x >> 2), g.node(x & low)});
  return g;
}

// The string a walk spells: its first k-mer, then the last base of each next one; none if the
// edges do not form a walk.
inline std::optional<Codes> spell(const KmerGraph& g, const std::vector<std::size_t>& walk) {
  Codes s;
  if (walk.empty()) return s;
  for (std::size_t i = 0; i + 1 < walk.size(); i++)
    if (g.graph.edges[walk[i]][1] != g.graph.edges[walk[i + 1]][0]) return std::nullopt;
  const Kmer x = g.edge_kmers[walk[0]];
  for (int i = g.k - 1; i >= 0; i--) s.push_back(static_cast<Code>((x >> (2 * i)) & 3));
  for (std::size_t i = 1; i < walk.size(); i++) s.push_back(static_cast<Code>(g.edge_kmers[walk[i]] & 3));
  return s;
}

// A string with the composition, by Hierholzer's algorithm; none if there is none.
inline std::optional<Codes> reconstruct(const KmerGraph& g) {
  const auto trail = morphe::combinatorics::eulerian_trail(g.graph);
  if (!trail) return std::nullopt;
  return spell(g, *trail);
}

// The number of distinct strings whose k-mer composition is the graph's.
inline Integer count_strings(const KmerGraph& g) {
  Integer trails = morphe::combinatorics::count_eulerian_trails(g.graph);
  if (trails.is_zero()) return trails;
  // ∏ m(x)! over the distinct k-mers; the edge k-mers are sorted.
  Integer copies = 1;
  for (std::size_t i = 0; i < g.edge_kmers.size();) {
    std::size_t j = i;
    while (j < g.edge_kmers.size() && g.edge_kmers[j] == g.edge_kmers[i]) j++;
    if (j - i > 1) copies *= morphe::arithmetic::factorial(j - i);
    i = j;
  }
  return morphe::arithmetic::exact_div(trails, copies);
}
inline Integer count_strings(const Codes& s, int k) { return count_strings(kmer_graph(composition(s, k))); }

// The de Bruijn graph B(σ, n): the words of length n over σ letters, word w = Σ wᵢ σ^{n−1−i};
// an edge x → y for each word of length n + 1 with prefix x and suffix y, edge e = x σ + a
// for the word x·a.
inline Multigraph complete_de_bruijn(unsigned sigma, unsigned n) {
  if (sigma < 1) throw std::invalid_argument("genomics::assembly::complete_de_bruijn: at least one letter");
  std::size_t V = 1;
  for (unsigned i = 0; i < n; i++) V *= sigma;
  Multigraph G;
  G.vertices = V;
  G.edges.reserve(V * sigma);
  for (std::size_t x = 0; x < V; x++)
    for (unsigned a = 0; a < sigma; a++) G.edges.push_back({x, (x * sigma + a) % V});
  return G;
}

// The number of de Bruijn sequences of order n ≥ 1 over σ letters, by the BEST theorem on
// B(σ, n − 1); and the closed form (σ!)^{σ^{n−1}} / σ^n.
inline Integer de_bruijn_sequences(unsigned sigma, unsigned n) {
  if (n < 1) throw std::invalid_argument("genomics::assembly::de_bruijn_sequences: n >= 1");
  return morphe::combinatorics::count_eulerian_circuits(complete_de_bruijn(sigma, n - 1));
}
inline Integer de_bruijn_sequences_formula(unsigned sigma, unsigned n) {
  if (n < 1) throw std::invalid_argument("genomics::assembly::de_bruijn_sequences_formula: n >= 1");
  unsigned long long power = 1;
  for (unsigned i = 1; i < n; i++) power *= sigma;
  return morphe::arithmetic::exact_div(morphe::arithmetic::pow(morphe::arithmetic::factorial(sigma), power),
                                       morphe::arithmetic::pow(Integer(static_cast<long long>(sigma)), n));
}

// One de Bruijn sequence of order n over letters 0, …, σ − 1, as the cyclic word of length σ^n
// read along a Hierholzer circuit of B(σ, n − 1) from edge 0 (the word 0^n).
inline std::vector<unsigned> de_bruijn_sequence(unsigned sigma, unsigned n) {
  if (n < 1 || sigma < 1) throw std::invalid_argument("genomics::assembly::de_bruijn_sequence: n >= 1, sigma >= 1");
  const auto c = morphe::combinatorics::eulerian_circuit(complete_de_bruijn(sigma, n - 1), 0);
  std::vector<unsigned> w;
  for (auto e : *c) w.push_back(static_cast<unsigned>(e % sigma));
  return w;
}

}  // namespace genomics::assembly
