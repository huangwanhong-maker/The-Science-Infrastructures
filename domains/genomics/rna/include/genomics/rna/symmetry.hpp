// The symmetries a structure grammar can be tied by: subgroups of the group of order 8 generated
// by the Klein group of the bases and reversal (genomics/grammar/strands.hpp), acting on the
// grammar's terminals and fixing its nonterminals (theory/rna-structure.md §3.4 and §5).
//
// On a pair rule X → a F b an element acting by a letter permutation g gives X → (g a) F (g b);
// a reversing one gives X → (g b) F (g a). So
//
//   ρ    (reversal)            ties (a, b) with (b, a):    A·U ~ U·A, G·C ~ C·G, G·U ~ U·G;
//   κ    (complement)          ties (a, b) with (κa, κb):  A·U ~ U·A, G·C ~ C·G, but G·U ~ C·A;
//   κρ   (reverse complement)  ties (a, b) with (κb, κa):  fixes A·U, U·A, G·C, C·G; G·U ~ A·C;
//
// and the strand group {1, κ, ρ, κρ} all three at once. The complement also ties the unpaired
// bases A ~ U and C ~ G, and reversal ties nothing but pairs. Knudsen and Hein's grammar admits
// only the letter permutations (κ); the reversible one admits every element.

#pragma once

#include <genomics/grammar/strands.hpp>

#include <morphe/algebra/action.hpp>
#include <morphe/algebra/finite_group.hpp>
#include <morphe/formal/equivariant.hpp>

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "scfg.hpp"
#include "training.hpp"

namespace genomics::rna {

namespace ma = morphe::algebra;

struct Symmetry {
  std::string name;
  const ma::GroupAction* letters = nullptr;  // on {A, C, G, U}, reversing words where its sign is −1
};

namespace detail {
inline const ma::FiniteGroup& order_two() {
  static const ma::FiniteGroup g = ma::cyclic(2);
  return g;
}
}  // namespace detail

// {1, ρ}: letters fixed, words reversed.
inline const ma::GroupAction& reversal_action() {
  static const ma::GroupAction a = ma::GroupAction::from_images(detail::order_two(), 4, {{0, 1, 2, 3}, {0, 1, 2, 3}}, {1, -1});
  return a;
}
// {1, κ}: A ↔ U, C ↔ G on letters, words kept in order.
inline const ma::GroupAction& complement_action() {
  static const ma::GroupAction a = ma::GroupAction::from_images(detail::order_two(), 4, {{0, 1, 2, 3}, {3, 2, 1, 0}}, {1, 1});
  return a;
}

// The symmetries tested here, smallest first.
inline std::vector<Symmetry> pair_symmetries() {
  return {{"ρ (reversal)", &reversal_action()},
          {"κ (complement)", &complement_action()},
          {"κρ (reverse complement)", &grammar::reverse_complement_group().words},
          {"{1, κ, ρ, κρ} (strand group)", &grammar::strand_group().words}};
}

inline mf::GrammarAction grammar_action(const Symmetry& s) { return {s.letters, {}}; }

// Whether the grammar's rule set is invariant under the symmetry.
inline bool admits(const StructureGrammar& sg, const Symmetry& s) { return mf::invariant(sg.grammar, grammar_action(s)); }

// The image of a pair under an element: (g·a, g·b), or (g·b, g·a) when g reverses.
inline std::pair<Code, Code> image(const Symmetry& s, std::uint32_t g, Code a, Code b) {
  const auto x = static_cast<Code>(s.letters->act(g, a)), y = static_cast<Code>(s.letters->act(g, b));
  return s.letters->reverses(g) ? std::pair<Code, Code>{y, x} : std::pair<Code, Code>{x, y};
}

// The least set of pairs containing p and closed under the symmetry: the pair set of the
// smallest invariant grammar containing a grammar of pair set p (morphe's closure, on pairs).
inline PairSet closure(const Symmetry& s, PairSet p) {
  PairSet out = 0;
  for (Code a = 0; a < 4; a++)
    for (Code b = 0; b < 4; b++)
      if (contains(p, a, b))
        for (std::uint32_t g = 0; g < s.letters->order(); g++) {
          const auto [x, y] = image(s, g, a, b);
          out = static_cast<PairSet>(out | (1u << (4 * x + y)));
        }
  return out;
}

// The tie along the symmetry's rule orbits.
inline Tie orbit_tie(const StructureGrammar& sg, const Symmetry& s) { return orbit_tie(sg, grammar_action(s), s.name); }

// g·x for an element of the symmetry's group.
inline Codes act(const Symmetry& s, std::uint32_t g, const Codes& x) { return s.letters->act_on_word(g, x); }

}  // namespace genomics::rna
