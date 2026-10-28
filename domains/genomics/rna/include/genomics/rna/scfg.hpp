// Stochastic context-free grammars of RNA secondary structure, built as morphe's grammars and
// computed with morphe's inside–outside: the probability of a sequence, the Viterbi structure,
// base-pair probabilities, the derivation of a given structure, and sampling
// (theory/rna-structure.md §3–§4).
//
// The grammars. Knudsen and Hein's grammar (Knudsen and Hein 1999, RNA secondary structure
// prediction using stochastic context-free grammars and evolutionary history, Bioinformatics 15,
// 446–454; the grammar G6 of Dowell and Eddy 2004, Evaluation of several lightweight stochastic
// context-free grammars for RNA secondary structure prediction, BMC Bioinformatics 5, 71) is
//
//   S → L S | L,     L → s | d F d,     F → d F d | L S,
//
// with s an unpaired base and d F d a pair; here each emission is its own rule, L → a for the
// four bases and L → a F b, F → a F b for every allowed pair (a, b), so that a pair rule's
// probability is the joint probability of pairing and of its two bases. S is a list of one or
// more components L, each an unpaired base or a pair closing F; F is the inside of a pair: a
// further pair stacked on it, or a list of two or more components. Three choices, each kept
// explicit:
//
//   the shape     Knudsen and Hein's lists are built from the left. Their rules S → L S and
//                 F → L S are not reversal-closed (reversed they read S → S L, F → S L), so
//                 reversal is no symmetry of the grammar in morphe's sense
//                 (formal/equivariant.hpp). The reversible shape builds lists from both ends,
//                   S → L | L L | L S L,   F → d F d | L L | L S L,
//                 derives the same structures and contains Knudsen and Hein's distributions as
//                 a sub-model (theory §3.3, as_reversible below); its rules are invariant under
//                 reversal as well as under the letters' permutations.
//   lonely pairs  allowed (Knudsen and Hein), or forbidden: L → s | d H d, H → d F d, so that
//                 every helix has at least two stacked pairs (theory §3.2).
//   the pairs     all sixteen, or a subset: the six canonical ones (structure.hpp), or their
//                 closure under a symmetry (symmetry.hpp). A pair outside the set has no rule.
//
// Every one of these grammars derives each of its structures by exactly one derivation (theory
// §3.2): the secondary structures with every pair allowed, every pair enclosing at least two
// bases (a hairpin loop of two or more) and, when lonely pairs are forbidden, every pair stacked
// on a neighbour. Its symmetries are the elements of the group of order 8 (the Klein group of the
// bases and reversal, genomics::grammar's base_and_strand_group) under which its rule set is
// invariant: the letter permutations that keep the pair set, and with the reversible shape the
// reversing elements that keep it.
//
// Probabilities are morphe's Arithmetic<double>: exact up to rounding while they stay in the
// normal range of doubles. A base costs about 1.5–2 nats, so sequences of up to about 350 bases
// are safe and longer ones underflow (refused, std::range_error); morphe has no scaled or
// logarithmic inside algorithm yet. Time O(n³ |R_binary|) for every computation on a sequence of
// n bases (|R_binary| = 66 for Knudsen and Hein's grammar with sixteen pairs, 42 for the
// reversible one with the six canonical pairs and no lonely ones); space O(n² |V_binary|).

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <morphe/formal/breaking.hpp>
#include <morphe/formal/grammar.hpp>
#include <morphe/formal/inside_outside.hpp>
#include <morphe/formal/semiring.hpp>

#include <calx/random/stream.hpp>

#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "pairing.hpp"
#include "structure.hpp"

namespace genomics::rna {

namespace mf = morphe::formal;

enum class Shape { knudsen_hein, reversible };

// A set of pairs (a, b), bit 4a + b.
using PairSet = std::uint16_t;
inline constexpr PairSet ALL_PAIRS = 0xFFFF;
inline constexpr bool contains(PairSet s, Code a, Code b) { return (s >> (4 * a + b)) & 1u; }
inline PairSet canonical_pairs(bool with_wobble = true) {
  PairSet s = 0;
  for (Code a = 0; a < 4; a++)
    for (Code b = 0; b < 4; b++)
      if (canonical(a, b, with_wobble)) s = static_cast<PairSet>(s | (1u << (4 * a + b)));
  return s;
}
inline std::string pair_names(PairSet s) {
  std::string out;
  for (Code a = 0; a < 4; a++)
    for (Code b = 0; b < 4; b++)
      if (contains(s, a, b)) {
        if (!out.empty()) out += ' ';
        out += "ACGU"[a];
        out += "ACGU"[b];
      }
  return out;
}

// What a rule of a structure grammar does: list structure, an unpaired base a, or a pair (a, b)
// around the nonterminal `inner`.
struct RuleInfo {
  enum Kind : std::uint8_t { structural, unpaired, pair } kind = structural;
  std::size_t lhs = 0;    // the left side's nonterminal index (S 0, L 1, F 2, H 3)
  Code a = 0, b = 0;
  std::size_t inner = 0;  // a pair rule's middle nonterminal
};

struct StructureGrammar {
  Shape shape = Shape::knudsen_hein;
  bool lonely_pairs = true;  // false: every helix has two pairs or more (the nonterminal H)
  PairSet pairs = ALL_PAIRS;
  std::string name;
  mf::Grammar grammar;
  mf::BinaryGrammar binary;
  std::vector<RuleInfo> info;  // per rule of grammar
  // The rules, by name (indices into grammar.rules; npos when absent).
  static constexpr std::size_t npos = SIZE_MAX;
  std::size_t S_LS = npos, S_L = npos, S_LL = npos, S_LSL = npos;
  std::size_t F_LS = npos, F_LL = npos, F_LSL = npos;
  std::size_t L_unpaired[4]{}, L_pair[4][4]{}, H_pair[4][4]{}, F_pair[4][4]{};

  std::size_t rules() const { return grammar.rules.size(); }
  std::size_t nonterminals() const { return grammar.nonterminals.size(); }
};

namespace detail {

inline StructureGrammar build(Shape shape, bool lonely_pairs, PairSet pairs) {
  if (!pairs) throw std::invalid_argument("genomics::rna::structure_grammar: a grammar needs at least one pair");
  StructureGrammar sg;
  sg.shape = shape;
  sg.lonely_pairs = lonely_pairs;
  sg.pairs = pairs;
  sg.name = std::string(shape == Shape::knudsen_hein ? "Knudsen–Hein" : "Knudsen–Hein, reversal-closed") + (lonely_pairs ? "" : ", no lonely pairs") +
            (pairs == ALL_PAIRS ? "" : ", pairs " + pair_names(pairs));
  auto& g = sg.grammar;
  g.terminals = 4;
  const mf::Symbol S = g.add_nonterminal("S"), L = g.add_nonterminal("L"), F = g.add_nonterminal("F");
  const mf::Symbol H = lonely_pairs ? F : g.add_nonterminal("H");
  const auto add = [&](mf::Symbol lhs, mf::Form rhs, RuleInfo info) {
    info.lhs = g.index(lhs);
    if (info.kind == RuleInfo::pair) info.inner = g.index(rhs[1]);
    g.add(lhs, std::move(rhs));
    sg.info.push_back(info);
    return g.rules.size() - 1;
  };
  for (Code a = 0; a < 4; a++)
    for (Code b = 0; b < 4; b++) sg.L_pair[a][b] = sg.H_pair[a][b] = sg.F_pair[a][b] = StructureGrammar::npos;
  if (shape == Shape::knudsen_hein) {
    sg.S_LS = add(S, {L, S}, {});
    sg.S_L = add(S, {L}, {});
  } else {
    sg.S_L = add(S, {L}, {});
    sg.S_LL = add(S, {L, L}, {});
    sg.S_LSL = add(S, {L, S, L}, {});
  }
  for (Code a = 0; a < 4; a++) sg.L_unpaired[a] = add(L, {a}, {RuleInfo::unpaired, 0, a, 0, 0});
  const auto pair_rules = [&](mf::Symbol lhs, mf::Symbol inner, std::size_t (&index)[4][4]) {
    for (Code a = 0; a < 4; a++)
      for (Code b = 0; b < 4; b++)
        if (contains(pairs, a, b)) index[a][b] = add(lhs, {a, inner, b}, {RuleInfo::pair, 0, a, b, 0});
  };
  pair_rules(L, H, sg.L_pair);
  if (!lonely_pairs) pair_rules(H, F, sg.H_pair);
  pair_rules(F, F, sg.F_pair);
  if (shape == Shape::knudsen_hein) sg.F_LS = add(F, {L, S}, {});
  else {
    sg.F_LL = add(F, {L, L}, {});
    sg.F_LSL = add(F, {L, S, L}, {});
  }
  g.validate();
  sg.binary = mf::BinaryGrammar::from(g);
  return sg;
}

inline std::vector<mf::Symbol> word(const Codes& x) {
  std::vector<mf::Symbol> w(x.size());
  for (std::size_t i = 0; i < x.size(); i++) {
    if (x[i] > 3) throw std::invalid_argument("genomics::rna: a base other than A, C, G, U");
    w[i] = x[i];
  }
  return w;
}

}  // namespace detail

// The grammar of a shape, a choice on lonely pairs and a pair set; built once and kept.
inline const StructureGrammar& structure_grammar(Shape shape, bool lonely_pairs = true, PairSet pairs = ALL_PAIRS) {
  static std::mutex mutex;
  static std::map<std::tuple<int, bool, PairSet>, std::unique_ptr<StructureGrammar>> built;
  std::lock_guard<std::mutex> lock(mutex);
  auto& slot = built[{static_cast<int>(shape), lonely_pairs, pairs}];
  if (!slot) slot = std::make_unique<StructureGrammar>(detail::build(shape, lonely_pairs, pairs));
  return *slot;
}
// Knudsen and Hein's grammar as published: lonely pairs allowed, every pair a rule.
inline const StructureGrammar& knudsen_hein() { return structure_grammar(Shape::knudsen_hein); }
inline const StructureGrammar& knudsen_hein_reversible() { return structure_grammar(Shape::reversible); }

// ---- parameters ---------------------------------------------------------------------------------------
// The rule probabilities: θ_r ≥ 0, summing to 1 over each left side.
using Parameters = std::vector<double>;

// Throws std::invalid_argument unless θ is a probability vector on each left side.
inline void validate(const StructureGrammar& sg, const Parameters& theta, double tolerance = 1e-9) {
  if (theta.size() != sg.rules()) throw std::invalid_argument("genomics::rna::validate: one probability per rule is needed");
  std::vector<double> sum(sg.nonterminals(), 0);
  for (std::size_t r = 0; r < theta.size(); r++) {
    if (!(theta[r] >= 0) || !std::isfinite(theta[r])) throw std::invalid_argument("genomics::rna::validate: a probability that is negative or not finite");
    sum[sg.info[r].lhs] += theta[r];
  }
  for (double s : sum)
    if (std::fabs(s - 1) > tolerance) throw std::invalid_argument("genomics::rna::validate: the rules of a left side do not sum to 1");
}

// A starting point for EM: every base equally likely unpaired; the pairs of each left side in
// proportion to 1, or to canonical_weight for the canonical ones (1: a start that knows nothing
// of pairing); P(L → a pair) = 1 − unpaired, P(F → a stacked pair) = stacked, H all pairs; the
// list rules of each left side sharing the rest equally. The canonical pairs are a
// reversal-closed set, so the start is invariant under every symmetry the grammar admits whose
// orbits keep canonical and other pairs apart, and EM is started from its projection onto a tie
// otherwise (training.hpp).
struct Start {
  double unpaired = 0.6, stacked = 0.5;
  double canonical_weight = 1;
};

inline Parameters starting_parameters(const StructureGrammar& sg, const Start& s = {}) {
  if (!(s.unpaired > 0 && s.unpaired < 1 && s.stacked > 0 && s.stacked < 1 && s.canonical_weight > 0))
    throw std::invalid_argument("genomics::rna::starting_parameters: probabilities in (0, 1) and a positive weight");
  Parameters t(sg.rules(), 0);
  std::vector<std::size_t> lists(sg.nonterminals(), 0);
  std::vector<double> pair_weight(sg.nonterminals(), 0);
  const auto weight = [&](const RuleInfo& i) { return canonical(i.a, i.b) ? s.canonical_weight : 1.0; };
  for (const auto& i : sg.info) {
    if (i.kind == RuleInfo::structural) lists[i.lhs]++;
    if (i.kind == RuleInfo::pair) pair_weight[i.lhs] += weight(i);
  }
  for (std::size_t r = 0; r < sg.rules(); r++) {
    const auto& i = sg.info[r];
    if (i.kind == RuleInfo::structural) t[r] = i.lhs == 0 ? 1.0 / static_cast<double>(lists[0]) : (1 - s.stacked) / static_cast<double>(lists[i.lhs]);
    else if (i.kind == RuleInfo::unpaired) t[r] = s.unpaired / 4;
    else t[r] = (i.lhs == 1 ? 1 - s.unpaired : i.lhs == 2 ? s.stacked : 1.0) * weight(i) / pair_weight[i.lhs];
  }
  validate(sg, t);
  return t;
}

// Knudsen and Hein's parameters carried to the reversible grammar of the same lonely-pair choice
// and pair set, giving every sequence and structure the same probability (theory §3.3): with
// p = θ(S → L S) and q = θ(F → L S),
//   S → L : 1 − p,   S → L L : p (1 − p),   S → L S L : p²,
//   F → L L : q (1 − p),   F → L S L : q p,
// and the emissions unchanged.
inline Parameters as_reversible(const StructureGrammar& from, const Parameters& kh) {
  if (from.shape != Shape::knudsen_hein) throw std::invalid_argument("genomics::rna::as_reversible: the parameters of a Knudsen–Hein grammar are needed");
  validate(from, kh);
  const auto& b = structure_grammar(Shape::reversible, from.lonely_pairs, from.pairs);
  Parameters t(b.rules(), 0);
  const double p = kh[from.S_LS], q = kh[from.F_LS];
  t[b.S_L] = 1 - p;
  t[b.S_LL] = p * (1 - p);
  t[b.S_LSL] = p * p;
  t[b.F_LL] = q * (1 - p);
  t[b.F_LSL] = q * p;
  for (Code x = 0; x < 4; x++) {
    t[b.L_unpaired[x]] = kh[from.L_unpaired[x]];
    for (Code y = 0; y < 4; y++) {
      if (!contains(from.pairs, x, y)) continue;
      t[b.L_pair[x][y]] = kh[from.L_pair[x][y]];
      t[b.F_pair[x][y]] = kh[from.F_pair[x][y]];
      if (!from.lonely_pairs) t[b.H_pair[x][y]] = kh[from.H_pair[x][y]];
    }
  }
  return t;
}

inline Parameters logarithms(const Parameters& theta) {
  Parameters l(theta.size());
  for (std::size_t r = 0; r < theta.size(); r++) l[r] = theta[r] > 0 ? std::log(theta[r]) : -std::numeric_limits<double>::infinity();
  return l;
}

// ---- the derivation of a structure ------------------------------------------------------------------------
// Whether a structure is one the grammar derives: every pair in its set and enclosing at least
// two bases, and (without lonely pairs) every pair stacked on a neighbour.
inline bool derivable(const StructureGrammar& sg, const Codes& x, const Structure& y) {
  if (x.size() != y.size() || x.empty() || !violation(y, 2).empty()) return false;
  for (std::size_t i = 0; i < y.size(); i++) {
    if (x[i] > 3) return false;
    if (y.partner[i] <= static_cast<std::int32_t>(i)) continue;
    const auto j = static_cast<std::size_t>(y.partner[i]);
    if (!contains(sg.pairs, x[i], x[j])) return false;
    if (!sg.lonely_pairs && !y.has(i + 1, j - 1) && !(i > 0 && y.has(i - 1, j + 1))) return false;
  }
  return true;
}

// The rules of the one derivation of (x, y), in leftmost (preorder) order; throws
// std::invalid_argument when the grammar does not derive the structure (derivable).
inline std::vector<std::size_t> derivation(const StructureGrammar& sg, const Codes& x, const Structure& y) {
  if (!derivable(sg, x, y)) throw std::invalid_argument("genomics::rna::derivation: the grammar does not derive this structure of this sequence");
  const std::size_t n = x.size();
  std::vector<std::size_t> out;
  // The components of [i, j): their first positions, each an unpaired base or a pair.
  const auto components = [&](std::size_t i, std::size_t j) {
    std::vector<std::size_t> c;
    for (std::size_t k = i; k < j; k = y.paired(k) ? static_cast<std::size_t>(y.partner[k]) + 1 : k + 1) c.push_back(k);
    return c;
  };
  std::function<void(std::size_t)> component;
  std::function<void(const std::vector<std::size_t>&, std::size_t, std::size_t)> list;
  std::function<void(std::size_t, std::size_t)> inside;  // F over the interior [i, j) of a pair
  component = [&](std::size_t k) {
    if (!y.paired(k)) { out.push_back(sg.L_unpaired[x[k]]); return; }
    const auto e = static_cast<std::size_t>(y.partner[k]);
    out.push_back(sg.L_pair[x[k]][x[e]]);
    if (sg.lonely_pairs) { inside(k + 1, e); return; }
    out.push_back(sg.H_pair[x[k + 1]][x[e - 1]]);
    inside(k + 2, e - 1);
  };
  // S over the components c[from, to).
  list = [&](const std::vector<std::size_t>& c, std::size_t from, std::size_t to) {
    const std::size_t m = to - from;
    if (sg.shape == Shape::knudsen_hein) {
      for (std::size_t t = from; t < to; t++) {
        out.push_back(t + 1 < to ? sg.S_LS : sg.S_L);
        component(c[t]);
      }
      return;
    }
    if (m == 1) { out.push_back(sg.S_L); component(c[from]); }
    else if (m == 2) { out.push_back(sg.S_LL); component(c[from]); component(c[from + 1]); }
    else { out.push_back(sg.S_LSL); component(c[from]); list(c, from + 1, to - 1); component(c[to - 1]); }
  };
  inside = [&](std::size_t i, std::size_t j) {
    if (j > i && y.has(i, j - 1)) {
      out.push_back(sg.F_pair[x[i]][x[j - 1]]);
      inside(i + 1, j - 1);
      return;
    }
    const auto c = components(i, j);
    if (sg.shape == Shape::knudsen_hein) { out.push_back(sg.F_LS); component(c[0]); list(c, 1, c.size()); }
    else if (c.size() == 2) { out.push_back(sg.F_LL); component(c[0]); component(c[1]); }
    else { out.push_back(sg.F_LSL); component(c[0]); list(c, 1, c.size() - 1); component(c.back()); }
  };
  const auto top = components(0, n);
  list(top, 0, top.size());
  return out;
}

// log P(x, y) under θ (−∞ when a rule of its derivation has probability 0).
inline double log_joint(const StructureGrammar& sg, const Parameters& theta, const Codes& x, const Structure& y) {
  double s = 0;
  for (std::size_t r : derivation(sg, x, y)) s += theta[r] > 0 ? std::log(theta[r]) : -std::numeric_limits<double>::infinity();
  return s;
}

// The count of each rule in a derivation: the complete-data statistics of a known structure.
inline std::vector<double> rule_counts(const StructureGrammar& sg, const Codes& x, const Structure& y) {
  std::vector<double> c(sg.rules(), 0);
  for (std::size_t r : derivation(sg, x, y)) c[r] += 1;
  return c;
}

// ---- the probability of a sequence -------------------------------------------------------------------------
namespace detail {
inline double checked_log(double z) {
  if (!(z >= DBL_MIN)) throw std::range_error("genomics::rna: the inside probability is zero or below the normal range of doubles");
  return std::log(z);
}
}  // namespace detail

// log P(x) = log Σ_y P(x, y), by the inside algorithm.
inline double log_likelihood(const StructureGrammar& sg, const Parameters& theta, const Codes& x) {
  const auto w = detail::word(x);
  if (w.empty()) throw std::invalid_argument("genomics::rna::log_likelihood: an empty sequence");
  const auto in = mf::inside<mf::Arithmetic<double>>(sg.binary, theta, w);
  return detail::checked_log(in.at(0, 0, w.size()));
}

// The expected number of uses of every rule in the derivations of x (the E-step), and log P(x).
inline std::vector<double> expected_rule_counts(const StructureGrammar& sg, const Parameters& theta, const Codes& x, double* log_p = nullptr) {
  const auto w = detail::word(x);
  if (w.empty()) throw std::invalid_argument("genomics::rna::expected_rule_counts: an empty sequence");
  double z = 0;
  auto c = mf::expected_counts<mf::Arithmetic<double>>(sg.grammar, sg.binary, theta, w, &z);
  const double l = detail::checked_log(z);
  if (log_p) *log_p = l;
  return c;
}

// ---- folding -------------------------------------------------------------------------------------------------
struct Folding {
  Structure structure;
  double log_probability = 0;  // log P(x, y) of the structure found
};

// The Viterbi structure: argmax_y P(x, y), by morphe's best parse over (max, +) on log θ.
inline Folding viterbi(const StructureGrammar& sg, const Parameters& theta, const Codes& x) {
  const auto w = detail::word(x);
  if (w.empty()) throw std::invalid_argument("genomics::rna::viterbi: an empty sequence");
  double best = 0;
  const auto tree = mf::best_parse<mf::MaxPlus<double>>(sg.grammar, sg.binary, logarithms(theta), w, &best);
  Folding f{Structure(w.size()), best};
  std::vector<const mf::ParseTree*> todo{&tree};
  while (!todo.empty()) {
    const mf::ParseTree* t = todo.back();
    todo.pop_back();
    if (t->rule < 0) continue;
    if (sg.info[static_cast<std::size_t>(t->rule)].kind == RuleInfo::pair) f.structure.add(t->begin, t->end - 1);
    for (const auto& c : t->children) todo.push_back(&c);
  }
  return f;
}

// The probability of every pair, from inside and outside: for i < j,
//   P(i·j | x) = Σ_{pair rules r = X → x_i Y x_j} O(X, i, j + 1) θ_r I(Y, i + 1, j) / P(x),
// the probability that a derivation uses a pair rule over the span [i, j] (theory §4.1), the
// grammar's analogue of McCaskill's partition-function pair probabilities (McCaskill 1990,
// The equilibrium partition function and base pair binding probabilities for RNA secondary
// structure, Biopolymers 29, 1105–1119). The original nonterminals keep their indices in
// morphe's binary form, so their inside and outside values are read from its charts directly.
inline PairProbabilities pair_probabilities(const StructureGrammar& sg, const Parameters& theta, const Codes& x) {
  const auto w = detail::word(x);
  const std::size_t n = w.size();
  if (n == 0) throw std::invalid_argument("genomics::rna::pair_probabilities: an empty sequence");
  using R = mf::Arithmetic<double>;
  const auto in = mf::inside<R>(sg.binary, theta, w);
  const auto out = mf::outside<R>(sg.binary, theta, w, in);
  PairProbabilities bp;
  bp.n = n;
  bp.p.assign(n * n, 0);
  const double z = in.at(0, 0, n);
  bp.log_likelihood = detail::checked_log(z);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 3; j < n; j++) {
      double v = 0;
      for (const std::size_t r : {sg.L_pair[w[i]][w[j]], sg.H_pair[w[i]][w[j]], sg.F_pair[w[i]][w[j]]}) {
        if (r == StructureGrammar::npos || theta[r] == 0) continue;
        const auto& info = sg.info[r];
        v += out.at(info.lhs, i, j + 1) * theta[r] * in.at(info.inner, i + 1, j);
      }
      bp.p[i * n + j] = bp.p[j * n + i] = v / z;
    }
  return bp;
}

// The MEA structure of x under the grammar (pairing.hpp), with its pair probabilities.
inline BestStructure<double> mea_fold(const StructureGrammar& sg, const Parameters& theta, const Codes& x, double gamma = 1.0,
                                      PairProbabilities* probabilities = nullptr) {
  auto bp = pair_probabilities(sg, theta, x);
  auto best = maximum_expected_accuracy(bp, gamma, 2);
  if (probabilities) *probabilities = std::move(bp);
  return best;
}

// ---- sampling --------------------------------------------------------------------------------------------------
struct Sample {
  Codes sequence;
  Structure structure;
};

// A sequence and its structure drawn from the grammar: morphe's leftmost derivation with each
// rule chosen by its probability, then replayed to place the pairs. Draws until the length lies
// in [min_length, max_length] (at most `attempts` times; then std::runtime_error).
inline Sample sample(const StructureGrammar& sg, const Parameters& theta, calx::random::Stream& rng, std::size_t min_length = 1,
                     std::size_t max_length = 1000, std::size_t attempts = 10000) {
  validate(sg, theta);
  const auto& g = sg.grammar;
  struct TooLong {};
  for (std::size_t attempt = 0; attempt < attempts; attempt++) {
    std::vector<std::uint32_t> chosen;
    std::vector<mf::Symbol> w;
    try {
      w = mf::derive(g, [&](mf::Symbol, const std::vector<std::uint32_t>& options) {
            if (chosen.size() > 4 * max_length) throw TooLong{};
            const double u = rng.uniform();
            double acc = 0;
            std::size_t pick = options.size();
            for (std::size_t k = 0; k < options.size() && pick == options.size(); k++) {
              acc += theta[options[k]];
              if (u < acc) pick = k;
            }
            // Rounding left u above the sum: the last rule of positive probability.
            if (pick == options.size())
              for (std::size_t k = options.size(); k-- > 0 && pick == options.size();)
                if (theta[options[k]] > 0) pick = k;
            chosen.push_back(options[pick]);
            return pick;
          }, 8 * max_length + 64);
    } catch (const TooLong&) {
      continue;
    } catch (const std::runtime_error&) {
      continue;
    }
    if (w.size() < min_length || w.size() > max_length) continue;
    Sample s;
    s.sequence.assign(w.begin(), w.end());
    s.structure = Structure(w.size());
    std::size_t next = 0, pos = 0;
    std::function<void(mf::Symbol)> replay = [&](mf::Symbol x) {
      if (g.terminal(x)) { pos++; return; }
      const std::size_t r = chosen.at(next++);
      const std::size_t start = pos;
      for (mf::Symbol c : g.rules[r].rhs) replay(c);
      if (sg.info[r].kind == RuleInfo::pair) s.structure.add(start, pos - 1);
    };
    replay(g.start());
    return s;
  }
  throw std::runtime_error("genomics::rna::sample: no derivation of a length in the range was drawn");
}

}  // namespace genomics::rna
