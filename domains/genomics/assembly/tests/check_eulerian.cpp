// Checks of assembly/eulerian.hpp: from the perfect k-spectrum of a random string, Hierholzer
// reconstructs it when it is the only string with that composition; the number of strings with
// a composition, by the BEST theorem and the multiplicities of repeated k-mers, against every
// string of length 9 enumerated, and against the distinct spellings of every Eulerian trail
// of the textbook example TAATGCCATGGGATGTT; a reconstruction always has the composition; the
// de Bruijn sequences of order n, 2^(2^(n−1) − n) binary ones (de Bruijn 1946) and
// (σ!)^(σ^(n−1)) / σ^n over σ letters, exactly; and one of them, spelled by Hierholzer, holds
// every word of length n once.

#include <genomics/assembly/eulerian.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <morphe/combinatorics/eulerian.hpp>

#include <algorithm>
#include <array>
#include <set>
#include <string>
#include <vector>

using namespace genomics;
using assembly::Integer;
using sequence::Codes;

namespace {

Codes from_index(std::size_t x, std::size_t length) {
  Codes s(length);
  for (std::size_t i = length; i-- > 0; x >>= 2) s[i] = static_cast<sequence::Code>(x & 3);
  return s;
}

// The strings spelled by every Eulerian trail: the circuits of the graph with an edge added from
// the sink to the source, cut there (or, for a balanced graph, every rotation of every circuit).
std::set<Codes> spellings(const assembly::KmerGraph& g) {
  std::set<Codes> out;
  const auto st = morphe::combinatorics::eulerian_status(g.graph);
  if (st.kind == morphe::combinatorics::EulerianKind::trail) {
    morphe::combinatorics::Multigraph h = g.graph;
    h.edges.insert(h.edges.begin(), std::array<std::size_t, 2>{st.end, st.start});
    morphe::combinatorics::for_each_eulerian_circuit(h, [&](const std::vector<std::size_t>& c) {
      std::vector<std::size_t> walk;
      for (std::size_t i = 1; i < c.size(); i++) walk.push_back(c[i] - 1);  // edge 0 is the added one
      out.insert(*assembly::spell(g, walk));
      return true;
    });
  } else if (st.kind == morphe::combinatorics::EulerianKind::circuit) {
    morphe::combinatorics::for_each_eulerian_circuit(g.graph, [&](const std::vector<std::size_t>& c) {
      for (std::size_t r = 0; r < c.size(); r++) {
        std::vector<std::size_t> walk(c.begin() + static_cast<std::ptrdiff_t>(r), c.end());
        walk.insert(walk.end(), c.begin(), c.begin() + static_cast<std::ptrdiff_t>(r));
        out.insert(*assembly::spell(g, walk));
      }
      return true;
    });
  }
  return out;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // Random strings of 2000 bases, k = 21: no repeated 20-mer, one string, reconstructed.
  {
    int bad = 0;
    for (std::uint32_t t = 0; t < 20; t++) {
      const Codes x = simulate::random_sequence(2000, 0.5, 70, t);
      const auto g = assembly::kmer_graph(assembly::composition(x, 21));
      const auto y = assembly::reconstruct(g);
      if (assembly::count_strings(g) != Integer(1) || !y || *y != x) bad++;
    }
    s.require("Eulerian: Hierholzer reconstructs 20 random strings of 2000 bases from their 21-mers, each the only string with its composition", bad == 0, std::to_string(bad) + " failures");
  }

  // Short strings with repeats: the count against all 4^9 strings of length 9.
  {
    const std::size_t L = 9, strings = std::size_t{1} << (2 * L);
    calx::random::Stream r(71, 0x45554c45, 0, 0);
    int bad = 0, reconstructed = 0, unique = 0, several = 0;
    std::size_t largest = 0;
    for (int k = 2; k <= 5; k++) {
      // The composition of every string of length 9, once for this k.
      std::vector<assembly::Composition> all(strings);
      for (std::size_t x = 0; x < strings; x++) all[x] = assembly::composition(from_index(x, L), k);
      for (int t = 0; t < 10; t++) {
        // Strings over two letters repeat their k-mers; over four they rarely do.
        Codes target(L);
        for (auto& c : target) c = static_cast<sequence::Code>(t % 2 ? r.below(2) : r.below(4));
        const auto comp = assembly::composition(target, k);
        const std::size_t brute = static_cast<std::size_t>(std::count(all.begin(), all.end(), comp));
        const auto g = assembly::kmer_graph(comp);
        if (assembly::count_strings(g) != Integer(static_cast<unsigned long long>(brute))) bad++;
        const auto y = assembly::reconstruct(g);
        if (y && assembly::composition(*y, k) == comp) reconstructed++;
        if (brute == 1) {
          unique++;
          if (!y || *y != target) bad++;
        } else {
          several++;
        }
        largest = std::max(largest, brute);
      }
    }
    s.require("Eulerian: the strings with a composition, trails / ∏ m(x)!, against all 4^9 strings of length 9, k = 2 to 5, 40 targets; the unique ones reconstructed", bad == 0, std::to_string(bad) + " failures");
    s.require("Eulerian: every reconstruction has the composition it was made from", reconstructed == 40,
              std::to_string(unique) + " targets with one string, " + std::to_string(several) + " with several, at most " + std::to_string(largest));
  }

  // The textbook example (Compeau and Pevzner 2015, Bioinformatics Algorithms, ch. 3): the
  // 3-mers of TAATGCCATGGGATGTT, among them ATG three times.
  {
    const Codes x = sequence::encode_dna("TAATGCCATGGGATGTT");
    const auto g = assembly::kmer_graph(assembly::composition(x, 3));
    const auto spelled = spellings(g);
    const Integer n = assembly::count_strings(g);
    s.require("Eulerian: TAATGCCATGGGATGTT has as many strings with its 3-mers as distinct spellings of its Eulerian trails", n == Integer(static_cast<long long>(spelled.size())) && spelled.count(x),
              n.str() + " strings, among them TAATGGGATGCCATGTT: " + (spelled.count(sequence::encode_dna("TAATGGGATGCCATGTT")) ? "yes" : "no"));
    s.require("Eulerian: its graph has 11 nodes and 15 edges, a trail from TA to TT", g.nodes.size() == 11 && g.graph.edges.size() == 15 &&
              morphe::combinatorics::eulerian_status(g.graph).kind == morphe::combinatorics::EulerianKind::trail);
  }

  // A balanced composition: a string that begins and ends with the same (k − 1)-mer.
  {
    const Codes x = sequence::encode_dna("ACGTTACGAC");  // AC … AC, k = 3
    const auto g = assembly::kmer_graph(assembly::composition(x, 3));
    std::size_t brute = 0;
    for (std::size_t i = 0; i < (std::size_t{1} << 20); i++) brute += assembly::composition(from_index(i, 10), 3) == assembly::composition(x, 3);
    s.require("Eulerian: a balanced composition counts every rotation of every circuit, against all 4^10 strings", assembly::count_strings(g) == Integer(static_cast<long long>(brute)) &&
              spellings(g).size() == brute, std::to_string(brute) + " strings");
  }

  // De Bruijn sequences.
  {
    bool ok = true;
    for (unsigned n = 1; n <= 6; n++)
      ok = ok && assembly::de_bruijn_sequences(2, n) == morphe::arithmetic::pow(Integer(2), (1ull << (n - 1)) - n);
    s.require("de Bruijn: the BEST theorem on B(2, n − 1) gives 2^(2^(n−1) − n) binary sequences of order n = 1 to 6 (de Bruijn 1946)", ok,
              "order 6: " + assembly::de_bruijn_sequences(2, 6).str());
    ok = true;
    for (unsigned sigma = 2; sigma <= 4; sigma++)
      for (unsigned n = 1; n <= 4; n++) ok = ok && assembly::de_bruijn_sequences(sigma, n) == assembly::de_bruijn_sequences_formula(sigma, n);
    s.require("de Bruijn: the BEST theorem on B(σ, n − 1) gives (σ!)^(σ^(n−1)) / σ^n for σ = 2, 3, 4 and n = 1 to 4", ok,
              "σ = 4, n = 3: " + assembly::de_bruijn_sequences(4, 3).str());
    // 24^16 / 64 = 12116574790945106558976 / 64, by hand.
    s.require("de Bruijn: the DNA de Bruijn sequences of order 3, (4!)^16 / 4^3 = 189321481108517289984", assembly::de_bruijn_sequences(4, 3).str() == "189321481108517289984", assembly::de_bruijn_sequences(4, 3).str());
    ok = true;
    for (unsigned sigma : {2u, 4u})
      for (unsigned n = 1; n <= (sigma == 2 ? 12u : 6u); n++) {
        const auto w = assembly::de_bruijn_sequence(sigma, n);
        std::size_t words = 1;
        for (unsigned i = 0; i < n; i++) words *= sigma;
        std::vector<int> seen(words, 0);
        for (std::size_t i = 0; i < w.size(); i++) {
          std::size_t x = 0;
          for (unsigned j = 0; j < n; j++) x = x * sigma + w[(i + j) % w.size()];
          seen[x]++;
        }
        ok = ok && w.size() == words && std::all_of(seen.begin(), seen.end(), [](int c) { return c == 1; });
      }
    s.require("de Bruijn: the sequence Hierholzer spells holds every word of length n once, cyclically (σ = 2, n ≤ 12; σ = 4, n ≤ 6)", ok);
  }
  return s.done();
}
