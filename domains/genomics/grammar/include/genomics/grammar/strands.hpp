// The symmetries of DNA words: the Klein group of the bases (nucleotides.hpp) and reversal, which
// together make a group of order 8 acting on words, and its strand subgroup {1, κ, ρ, κρ}:
// complement, reversal, and reverse complement, the map from a word to the other strand read in
// its own direction. Each element acts on words by morphe's GroupAction, reversing words when its
// sign is −1 (theory/group-grammars.md §4).
//
// Canonical k-mers. The words of length k up to a group are counted exactly by Burnside's lemma
// with the fixed-word counts of morphe/algebra/action.hpp. Under the reverse complement alone
// this gives the familiar (4^k + 4^{k/2})/2 for even k (the reverse-complement palindromes
// 4^{k/2} are fixed) and 4^k/2 for odd k; under the whole strand group it gives
// (4^k + 4^⌈k/2⌉ + [k even] 4^{k/2})/4, since no word is its own complement and reversal fixes
// 4^⌈k/2⌉.
//
// Which symmetries a genome keeps. A group element g acts on the frequencies f of a genome's
// k-mers (one strand, as the record gives it) by (g·f)(w) = f(g⁻¹·w); the total variation
// distance TV(f, g·f) = ½ Σ_w |f(w) − f(g⁻¹·w)| is 0 exactly when the statistics are invariant
// under g. The answer for real genomes is an empirical fact, not a structural one: Chargaff's
// second parity rule, that a single strand's k-mer counts nearly equal its reverse complements'
// for long genomes (Rudner, Karkas and Chargaff 1968, PNAS 60, 921–922, for the bases), says
// that of the eight elements only the identity and the reverse complement nearly keep them. This
// is the draft's question Q2 (group discovery), asked of data: the symmetry group a genome's
// statistics have is measured here, element by element.

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <morphe/algebra/action.hpp>
#include <morphe/algebra/finite_group.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "nucleotides.hpp"

namespace genomics::grammar {

struct WordGroup {
  ma::FiniteGroup group;
  ma::GroupAction words;              // on {A, C, G, T}, signs reversing words
  std::vector<std::string> names;     // of the elements
  std::uint32_t complement = 0, reverse = 0, reverse_complement = 0;
};

// The group of order 8, V₄ × ⟨ρ⟩: element 2g + r is the base symmetry g of klein(), followed by
// reversal when r = 1.
inline const WordGroup& base_and_strand_group() {
  static const WordGroup w = [] {
    WordGroup out;
    const auto& k = klein();
    out.group = ma::direct_product(k.group, ma::cyclic(2));
    std::vector<std::vector<std::uint32_t>> images(8, std::vector<std::uint32_t>(4));
    std::vector<std::int8_t> sign(8);
    const std::array<const char*, 4> base_names{"1", "τ (transition)", "κ (complement)", "σ (A↔C, G↔T)"};
    for (std::uint32_t e = 0; e < 8; e++) {
      const std::uint32_t g = e / 2, r = e % 2;
      for (std::uint32_t x = 0; x < 4; x++) images[e][x] = k.bases.act(g, x);
      sign[e] = r ? -1 : 1;
      out.names.push_back(r ? (g == 0 ? std::string("ρ (reversal)") : std::string(base_names[g]) + " ∘ ρ") : std::string(base_names[g]));
    }
    out.names[2 * k.complement + 1] = "κ ∘ ρ (reverse complement)";
    out.words = ma::GroupAction::from_images(out.group, 4, images, sign);
    out.complement = 2 * k.complement;
    out.reverse = 1;
    out.reverse_complement = 2 * k.complement + 1;
    return out;
  }();
  return w;
}

// The strand group {1, κ, ρ, κρ} ≅ C₂ × C₂: element 2a + b is the complement when a = 1 followed
// by reversal when b = 1.
inline const WordGroup& strand_group() {
  static const WordGroup w = [] {
    WordGroup out;
    const auto& k = klein();
    out.group = ma::direct_product(ma::cyclic(2), ma::cyclic(2));
    std::vector<std::vector<std::uint32_t>> images(4, std::vector<std::uint32_t>(4));
    std::vector<std::int8_t> sign(4);
    for (std::uint32_t e = 0; e < 4; e++) {
      for (std::uint32_t x = 0; x < 4; x++) images[e][x] = (e / 2) ? k.bases.act(k.complement, x) : x;
      sign[e] = (e % 2) ? -1 : 1;
    }
    out.words = ma::GroupAction::from_images(out.group, 4, images, sign);
    out.names = {"1", "ρ (reversal)", "κ (complement)", "κ ∘ ρ (reverse complement)"};
    out.complement = 2;
    out.reverse = 1;
    out.reverse_complement = 3;
    return out;
  }();
  return w;
}

// The reverse-complement group {1, κρ} alone.
inline const WordGroup& reverse_complement_group() {
  static const WordGroup w = [] {
    WordGroup out;
    const auto& k = klein();
    out.group = ma::cyclic(2);
    std::vector<std::vector<std::uint32_t>> images(2, std::vector<std::uint32_t>(4));
    for (std::uint32_t x = 0; x < 4; x++) { images[0][x] = x; images[1][x] = k.bases.act(k.complement, x); }
    out.words = ma::GroupAction::from_images(out.group, 4, images, {1, -1});
    out.names = {"1", "κ ∘ ρ (reverse complement)"};
    out.reverse_complement = 1;
    return out;
  }();
  return w;
}

// ---- k-mer statistics and their symmetry ----------------------------------------------------------

// The counts of the k-mers of a sequence (one strand), indexed by the base-4 number of the word,
// first base most significant; windows with an N are skipped.
inline std::vector<double> kmer_frequencies(const sequence::Codes& x, std::size_t k) {
  if (k == 0 || k > 12) throw std::invalid_argument("genomics::grammar::kmer_frequencies: k must lie in 1..12");
  std::vector<double> f(std::size_t{1} << (2 * k), 0.0);
  std::uint64_t code = 0;
  std::size_t valid = 0;
  const std::uint64_t mask = (std::uint64_t{1} << (2 * k)) - 1;
  double total = 0;
  for (Code b : x) {
    if (b > 3) { valid = 0; code = 0; continue; }
    code = ((code << 2) | b) & mask;
    if (++valid >= k) { f[code] += 1; total += 1; }
  }
  if (total > 0)
    for (auto& v : f) v /= total;
  return f;
}

// The index of g·w for the k-mer of index w.
inline std::size_t act_on_kmer(const ma::GroupAction& a, std::uint32_t g, std::size_t w, std::size_t k) {
  std::vector<std::uint8_t> word(k);
  for (std::size_t i = 0; i < k; i++) word[k - 1 - i] = static_cast<std::uint8_t>((w >> (2 * i)) & 3);
  const auto image = a.act_on_word(g, word);
  std::size_t out = 0;
  for (auto b : image) out = (out << 2) | b;
  return out;
}

struct Invariance {
  std::uint32_t element = 0;
  std::string name;
  double total_variation = 0;  // ½ Σ_w |f(w) − f(g⁻¹·w)|
};

inline std::vector<Invariance> kmer_invariance(const sequence::Codes& x, std::size_t k, const WordGroup& w) {
  const auto f = kmer_frequencies(x, k);
  std::vector<Invariance> out;
  for (std::uint32_t g = 0; g < w.group.order(); g++) {
    const std::uint32_t gi = w.group.inv(g);
    double tv = 0;
    for (std::size_t u = 0; u < f.size(); u++) tv += std::fabs(f[u] - f[act_on_kmer(w.words, gi, u, k)]);
    out.push_back({g, w.names[g], tv / 2});
  }
  return out;
}

}  // namespace genomics::grammar
