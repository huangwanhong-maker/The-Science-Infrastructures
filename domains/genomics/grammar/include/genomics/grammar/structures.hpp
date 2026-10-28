// Grammars of nucleic-acid structure: the hairpin (a stem of paired bases closing a loop), as a
// context-free grammar over {A, C, G, T}, its symmetry, and a grammar that finds one anywhere in
// a sequence (theory/group-grammars.md §5).
//
// The hairpin language H(s, l₁, l₂) is the set of words u v u′ with |u| = |u′| ≥ s, the pairs
// (u_i, u′_{|u|+1−i}) all allowed, and l₁ ≤ |v| ≤ l₂. With Watson–Crick pairs only, u′ is the
// reverse complement of u; with the wobble pairs G·T and T·G (G·U in RNA, where T reads as U)
// it is not. The grammar emits the stem from the outside in,
//
//   Stem_i → a Stem_{i+1} b            (i < s),          one rule per allowed pair (a, b),
//   Stem_s → a Stem_s b | a Loop_m b   (l₁ ≤ m ≤ l₂),
//   Loop_m → x Loop_{m−2} y  (m ≥ 3),  Loop_2 → x y,  Loop_1 → x,  for all letters x, y,
//
// so that every rule is symmetric under reversal (a pair read backwards is the reversed pair,
// and the loop is built outside in too). Its symmetry group is then exact and readable from the
// rules (formal/equivariant.hpp's invariant): with Watson–Crick pairs the grammar is invariant
// under the whole strand group {1, κ, ρ, κρ}, since (κa, κb) is a Watson–Crick pair when (a, b)
// is; with wobble pairs it keeps reversal but loses the complement, since κ(G·T) = C·A pairs in
// neither. The chemistry of the G·U pair breaks the symmetry, and the grammar records which part.
//
// A hairpin anywhere: Region → Hairpin | X Hairpin | Hairpin X | X Hairpin X, X → x X | x,
// ε-free so that both Earley's parser and the inside algorithm take it. The cost on a window of
// length n is O(n³) in general; windows of a few tens of bases are the use here, not whole
// genomes. Its images under the breaking maps of nucleotides.hpp (to R/Y, W/S or M/K) are the
// coarse grammars of graduated breaking: parsed first, they mask the fine parse exactly
// (formal/breaking.hpp's coarse_mask).

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <morphe/formal/grammar.hpp>

#include "nucleotides.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace genomics::grammar {

struct HairpinSpec {
  std::size_t min_stem = 4;
  std::size_t min_loop = 3, max_loop = 8;
  bool wobble = false;  // allow G·T and T·G pairs (G·U in RNA)
};

inline std::vector<std::pair<Code, Code>> allowed_pairs(bool wobble) {
  std::vector<std::pair<Code, Code>> p{{sequence::A, sequence::T}, {sequence::T, sequence::A}, {sequence::C, sequence::G}, {sequence::G, sequence::C}};
  if (wobble) { p.push_back({sequence::G, sequence::T}); p.push_back({sequence::T, sequence::G}); }
  return p;
}

// The hairpin grammar, its start symbol the outermost stem nonterminal.
inline mf::Grammar hairpin_grammar(const HairpinSpec& s) {
  if (s.min_stem < 1 || s.min_loop < 1 || s.max_loop < s.min_loop) throw std::invalid_argument("genomics::grammar::hairpin_grammar: a stem and a loop of at least one base");
  mf::Grammar g;
  g.terminals = 4;
  std::vector<mf::Symbol> stem;
  for (std::size_t i = 1; i <= s.min_stem; i++) stem.push_back(g.add_nonterminal("Stem" + std::to_string(i)));
  std::vector<mf::Symbol> loop(s.max_loop + 1, 0);
  for (std::size_t m = 1; m <= s.max_loop; m++) loop[m] = g.add_nonterminal("Loop" + std::to_string(m));
  const auto pairs = allowed_pairs(s.wobble);
  for (std::size_t i = 0; i < s.min_stem; i++)
    for (const auto& [a, b] : pairs) {
      if (i + 1 < s.min_stem) g.add(stem[i], mf::Form{a, stem[i + 1], b});
      else {
        g.add(stem[i], mf::Form{a, stem[i], b});
        for (std::size_t m = s.min_loop; m <= s.max_loop; m++) g.add(stem[i], mf::Form{a, loop[m], b});
      }
    }
  for (std::size_t m = 1; m <= s.max_loop; m++)
    for (Code x = 0; x < 4; x++) {
      if (m == 1) { g.add(loop[m], mf::Form{x}); continue; }
      for (Code y = 0; y < 4; y++) {
        if (m == 2) g.add(loop[m], mf::Form{x, y});
        else g.add(loop[m], mf::Form{x, loop[m - 2], y});
      }
    }
  return g;
}

// The same with flanks of any letters on both sides: the words containing a hairpin.
inline mf::Grammar contains_hairpin_grammar(const HairpinSpec& s) {
  mf::Grammar h = hairpin_grammar(s);
  mf::Grammar g;
  g.terminals = 4;
  const mf::Symbol region = g.add_nonterminal("Region");
  const mf::Symbol flank = g.add_nonterminal("Flank");
  // Copy the hairpin's nonterminals after these two.
  const std::size_t offset = g.nonterminals.size();
  for (const auto& name : h.nonterminals) g.add_nonterminal(name);
  const auto shift = [&](const mf::Form& f) {
    mf::Form out;
    for (mf::Symbol x : f) out.push_back(h.terminal(x) ? x : g.nonterminal(offset + h.index(x)));
    return out;
  };
  for (const auto& r : h.rules) g.add(shift(r.lhs), shift(r.rhs));
  // ε-free, so that the inside algorithm takes it as well as Earley's parser: a flank is one
  // letter or more, and the region has a flank on neither, either or both sides.
  const mf::Symbol hp = g.nonterminal(offset);
  g.add(region, mf::Form{hp});
  g.add(region, mf::Form{flank, hp});
  g.add(region, mf::Form{hp, flank});
  g.add(region, mf::Form{flank, hp, flank});
  for (Code x = 0; x < 4; x++) {
    g.add(flank, mf::Form{x, flank});
    g.add(flank, mf::Form{x});
  }
  return g;
}

// The definition, by brute force: whether w is in H(s, l₁, l₂).
inline bool is_hairpin(const std::vector<mf::Symbol>& w, const HairpinSpec& s) {
  const auto pairs = allowed_pairs(s.wobble);
  const auto pairs_with = [&](mf::Symbol a, mf::Symbol b) {
    for (const auto& [x, y] : pairs) if (x == a && y == b) return true;
    return false;
  };
  const std::size_t n = w.size();
  for (std::size_t stem = s.min_stem; 2 * stem + s.min_loop <= n; stem++) {
    const std::size_t loop = n - 2 * stem;
    if (loop < s.min_loop || loop > s.max_loop) continue;
    bool ok = true;
    for (std::size_t i = 0; i < stem && ok; i++) ok = pairs_with(w[i], w[n - 1 - i]);
    if (ok) return true;
  }
  return false;
}

inline bool contains_hairpin(const std::vector<mf::Symbol>& w, const HairpinSpec& s) {
  for (std::size_t i = 0; i < w.size(); i++)
    for (std::size_t j = i + 1; j <= w.size(); j++)
      if (is_hairpin(std::vector<mf::Symbol>(w.begin() + static_cast<std::ptrdiff_t>(i), w.begin() + static_cast<std::ptrdiff_t>(j)), s)) return true;
  return false;
}

}  // namespace genomics::grammar
