// The workbench's operations on group grammars (src/views/grammar.js and src/tools/grammar.js): the
// user's draft "Towards a group context-sensitive grammar for genome analysis" (GE-CSG), made exact
// in domains/genomics/theory/group-grammars.md, whose sections the comments below cite.
//
//   grammar.nucleotides   the Klein group of the bases, its subgroups (IUPAC's codes) and characters
//                         (the classifications), and the chains of breakings (§3)
//   grammar.characters    a finite group's exact character table, its Frobenius–Schur indicators and
//                         the Wedderburn blocks of ℂ[G] and ℚ[G] (§7)
//   grammar.parse         a grammar and a group acting on its letters: invariance, rule orbits by
//                         Burnside's lemma, Earley's parse, the equivariant parser beside the draft's
//                         orbit-only scan, spans shared across orbits, and the canonical search of a
//                         monotone grammar (§1, §2)
//   grammar.hairpins      every hairpin of a short sequence by the inside algorithm over all spans,
//                         the grammar's symmetry, and graduated breaking over R/Y, W/S and M/K (§2.5, §5)
//   grammar.symmetry      the eight symmetries of a genome's k-mer statistics and its eight Wedderburn
//                         channels (§4, §7)
//   grammar.code          a genetic code as a function on (C₂)⁶: its invariance profile, stabiliser and
//                         power spectrum (§6)
//
// Every computation is morphe's (algebra/action, characters, group_algebra, symmetric_representations;
// formal/grammar, earley, inside_outside, equivariant, breaking) or the genomics domain's (grammar/):
// the engine reads the request, holds it within bounds, and reports. The bounds keep a worker
// responsive: grammars of at most 300 rules (3,000 in a closure) under groups of at most 128 elements,
// words of at most 64 letters (14 for the exponential search of a monotone grammar, which also stops
// after two million forms), hairpins of at most 200 bases in a stretch of at most 400 (some seconds at
// worst), groups of at most 400 elements for their characters (their Wedderburn blocks measured by
// exact elimination up to 64).

#pragma once

#include "genomes.hpp"
#include "inputs.hpp"
#include "protocol.hpp"

#include <genomics/grammar/channels.hpp>
#include <genomics/grammar/code_symmetry.hpp>
#include <genomics/grammar/nucleotides.hpp>
#include <genomics/grammar/strands.hpp>
#include <genomics/grammar/structures.hpp>
#include <genomics/sequence/genetic_code.hpp>

#include <morphe/algebra/action.hpp>
#include <morphe/algebra/characters.hpp>
#include <morphe/algebra/finite_group.hpp>
#include <morphe/algebra/group_algebra.hpp>
#include <morphe/algebra/symmetric_representations.hpp>
#include <morphe/formal/breaking.hpp>
#include <morphe/formal/earley.hpp>
#include <morphe/formal/equivariant.hpp>
#include <morphe/formal/grammar.hpp>
#include <morphe/formal/inside_outside.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace workbench {

namespace grammar_ops {

namespace ma = morphe::algebra;
namespace mf = morphe::formal;
namespace gg = genomics::grammar;

using RuleSet = std::set<mf::Rule, mf::detail::RuleLess>;

inline double since(std::chrono::steady_clock::time_point t0) {
  return std::round(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() * 1e6) / 1e6;
}

// ---- groups acting on letters ---------------------------------------------------------------------
// A group acting on t letters with signs (morphe's GroupAction), held where its pointers stay valid.
// It is the closure of generators given as permutations of the letters and of two points more, which
// a reversing element exchanges: the sign is then part of the permutation and closes with it, so
// that the closure of (letters, ε) pairs is computed by morphe's closure of permutations.
struct Acting {
  ma::FiniteGroup group;
  ma::GroupAction action;
  mf::GrammarAction on_grammar;
  std::vector<ma::Permutation> elements;  // of the t + 2 points
  std::size_t letters = 0;
};

inline ma::Permutation encode(const std::vector<std::uint32_t>& letters, bool reverse) {
  ma::Permutation p(letters.begin(), letters.end());
  const auto t = static_cast<std::uint32_t>(letters.size());
  p.push_back(reverse ? t + 1 : t);
  p.push_back(reverse ? t : t + 1);
  return p;
}

inline std::unique_ptr<Acting> acting(const std::vector<ma::Permutation>& generators, std::size_t t, std::size_t limit = 128) {
  auto a = std::make_unique<Acting>();
  a->letters = t;
  ma::PermutationGroup P;
  try {
    P = ma::permutation_group(generators, t + 2, "G", limit);
  } catch (const std::length_error&) {
    throw Fault("the generators make a group of more than " + std::to_string(limit) + " elements, more than this takes");
  }
  a->group = std::move(P.group);
  a->elements = std::move(P.elements);
  std::vector<std::vector<std::uint32_t>> images;
  std::vector<std::int8_t> sign;
  for (const auto& p : a->elements) {
    images.emplace_back(p.begin(), p.begin() + static_cast<std::ptrdiff_t>(t));
    sign.push_back(p[t] == t ? 1 : -1);
  }
  a->action = ma::GroupAction::from_images(a->group, t, images, sign);
  a->on_grammar = mf::GrammarAction{&a->action, {}};
  return a;
}

// The subgroup of some of a group's elements (closed under the product), as a group of its own.
inline std::unique_ptr<Acting> subgroup(const Acting& a, const std::vector<std::uint32_t>& members) {
  std::vector<ma::Permutation> generators;
  for (auto g : members) generators.push_back(a.elements[g]);
  return acting(generators, a.letters, a.elements.size());
}

// An element's name: 1; ρ for reversal; κ, τ, σ for the involutions of the bases; else its cycles
// (a 2-cycle written a↔b), followed by ρ when it reverses.
inline std::string element_name(const Acting& a, std::uint32_t g, const std::string& alphabet) {
  const auto& p = a.elements[g];
  const std::size_t t = a.letters;
  const bool reverse = p[t] != t;
  bool identity = true;
  for (std::size_t i = 0; i < t; i++) identity = identity && p[i] == i;
  std::string letters;
  if (!identity && alphabet == "ACGT") {
    static const std::array<std::pair<const char*, std::array<std::uint32_t, 4>>, 3> known{{{"κ", {3, 2, 1, 0}}, {"τ", {2, 3, 0, 1}}, {"σ", {1, 0, 3, 2}}}};
    for (const auto& [name, q] : known)
      if (std::equal(q.begin(), q.end(), p.begin())) letters = name;
  }
  if (!identity && letters.empty()) {
    std::vector<char> seen(t, 0);
    for (std::size_t i = 0; i < t; i++) {
      if (seen[i] || p[i] == i) continue;
      std::vector<std::size_t> cycle;
      for (std::size_t j = i; !seen[j]; j = p[j]) { seen[j] = 1; cycle.push_back(j); }
      if (!letters.empty()) letters += " ";
      if (cycle.size() == 2) letters += std::string(1, alphabet[cycle[0]]) + "↔" + alphabet[cycle[1]];
      else {
        letters += "(";
        for (std::size_t k = 0; k < cycle.size(); k++) letters += (k ? " " : "") + std::string(1, alphabet[cycle[k]]);
        letters += ")";
      }
    }
  }
  if (identity) return reverse ? "ρ" : "1";
  if (!reverse) return letters;
  return letters == "κ" || letters == "τ" || letters == "σ" ? letters + "ρ" : letters + " ∘ ρ";
}

inline Json group_json(const Acting& a, const std::string& alphabet) {
  Json elements = Json::array();
  for (std::uint32_t g = 0; g < a.group.order(); g++) {
    std::string image;
    for (std::size_t i = 0; i < a.letters; i++) image += alphabet[a.action.act(g, static_cast<std::uint32_t>(i))];
    elements.push(Json::object({{"name", element_name(a, g, alphabet)}, {"letters", image}, {"reverses", a.action.reverses(g)},
                                {"order", a.group.element_order(g)}}));
  }
  return Json::object({{"order", a.group.order()}, {"elements", elements}, {"abelian", a.group.is_abelian()}});
}

// ---- grammars written as text ----------------------------------------------------------------------
struct Written {
  mf::Grammar grammar;
  std::string alphabet;
};

inline std::string symbol_name(const mf::Grammar& g, const std::string& alphabet, mf::Symbol s) {
  return g.terminal(s) ? std::string(1, alphabet[s]) : g.nonterminals[g.index(s)];
}

inline std::string form_text(const mf::Grammar& g, const std::string& alphabet, const mf::Form& f) {
  if (f.empty()) return "ε";
  std::string out;
  for (std::size_t i = 0; i < f.size(); i++) out += (i ? " " : "") + symbol_name(g, alphabet, f[i]);
  return out;
}

inline std::string rule_text(const mf::Grammar& g, const std::string& alphabet, const mf::Rule& r) {
  return form_text(g, alphabet, r.lhs) + " → " + form_text(g, alphabet, r.rhs);
}

// The tokens of one side of a rule: whitespace separates symbols; ε (or eps) is the empty word; a
// one-character token is a terminal when it is in the alphabet and a nonterminal otherwise; letters
// run together are split when they hold a terminal and are otherwise terminals or capitals (aSb,
// cB), and are a nonterminal's name otherwise (Stem, Loop2).
inline std::vector<std::string> tokens(const Json& side, const std::string& alphabet) {
  std::vector<std::string> raw;
  if (side.is_array()) {
    for (const auto& x : side.as_array()) raw.push_back(x.text());
  } else if (side.is_string()) {
    std::string cur;
    for (char c : side.as_string()) {
      if (std::isspace(static_cast<unsigned char>(c))) {
        if (!cur.empty()) raw.push_back(cur);
        cur.clear();
      } else cur += c;
    }
    if (!cur.empty()) raw.push_back(cur);
  } else if (!side.is_null()) {
    throw Fault("a side of a rule is a string of symbols");
  }
  std::vector<std::string> out;
  for (const auto& tok : raw) {
    if (tok.empty() || tok == "ε" || tok == "eps" || tok == "epsilon") continue;
    if (tok.size() == 1) { out.push_back(tok); continue; }
    bool has_terminal = false, splittable = true;
    for (unsigned char c : tok) {
      if (c < 128 && alphabet.find(static_cast<char>(c)) != std::string::npos) has_terminal = true;
      else if (!(c >= 'A' && c <= 'Z')) splittable = false;
    }
    if (has_terminal && splittable)
      for (char c : tok) out.push_back(std::string(1, c));
    else out.push_back(tok);
  }
  return out;
}

inline std::string alphabet_of(const Json& request) {
  const std::string a = str(request, "alphabet");
  if (a.empty() || a.size() > 16) throw Fault("an alphabet of 1 to 16 letters is needed");
  for (std::size_t i = 0; i < a.size(); i++) {
    const unsigned char c = static_cast<unsigned char>(a[i]);
    if (c < 33 || c > 126 || c == '|') throw Fault("the alphabet's letters are printable ASCII characters other than |");
    if (a.find(a[i]) != i) throw Fault(std::string("the letter ") + a[i] + " is in the alphabet twice");
  }
  return a;
}

inline Written written_grammar(const Json& request) {
  Written w;
  w.alphabet = alphabet_of(request);
  const std::string& alphabet = w.alphabet;
  const Json& rules = request["rules"];
  if (!rules.is_array() || rules.size() == 0) throw Fault("give the rules as [left, right] pairs");
  if (rules.size() > 300) throw Fault("a grammar here has at most 300 rules");
  const auto terminal = [&](const std::string& tok) { return tok.size() == 1 && alphabet.find(tok[0]) != std::string::npos; };
  std::vector<std::pair<std::vector<std::string>, std::vector<std::string>>> sides;
  std::vector<std::string> names;
  const auto meet = [&](const std::vector<std::string>& ts) {
    for (const auto& tok : ts)
      if (!terminal(tok) && std::find(names.begin(), names.end(), tok) == names.end()) names.push_back(tok);
  };
  for (const auto& r : rules.as_array()) {
    if (!r.is_array() || r.size() != 2) throw Fault("each rule is a pair [left, right]");
    sides.emplace_back(tokens(r[0], alphabet), tokens(r[1], alphabet));
    if (sides.back().first.empty()) throw Fault("a rule has an empty left side");
    meet(sides.back().first);
    meet(sides.back().second);
  }
  std::string start = str(request, "start");
  if (start.empty()) {
    for (const auto& tok : sides[0].first)
      if (!terminal(tok)) { start = tok; break; }
  }
  if (start.empty() || std::find(names.begin(), names.end(), start) == names.end()) throw Fault("the start symbol " + (start.empty() ? std::string("is missing") : start + " occurs in no rule"));
  w.grammar.terminals = alphabet.size();
  w.grammar.add_nonterminal(start);
  for (const auto& n : names)
    if (n != start) w.grammar.add_nonterminal(n);
  const auto symbol = [&](const std::string& tok) -> mf::Symbol {
    if (terminal(tok)) return static_cast<mf::Symbol>(alphabet.find(tok[0]));
    const auto at = std::find(w.grammar.nonterminals.begin(), w.grammar.nonterminals.end(), tok) - w.grammar.nonterminals.begin();
    return w.grammar.nonterminal(static_cast<std::size_t>(at));
  };
  for (const auto& [l, r] : sides) {
    mf::Form lhs, rhs;
    for (const auto& tok : l) lhs.push_back(symbol(tok));
    for (const auto& tok : r) rhs.push_back(symbol(tok));
    w.grammar.add(std::move(lhs), std::move(rhs));
  }
  try {
    w.grammar.validate();
  } catch (const std::invalid_argument& e) {
    throw Fault(detail::message_of(e));
  }
  return w;
}

inline std::vector<mf::Symbol> word_of(const std::string& text, const std::string& alphabet) {
  std::vector<mf::Symbol> w;
  for (char c : text) {
    if (std::isspace(static_cast<unsigned char>(c))) continue;
    const auto at = alphabet.find(c);
    if (at == std::string::npos) throw Fault(std::string("the word has a letter, ") + c + ", outside the alphabet");
    w.push_back(static_cast<mf::Symbol>(at));
  }
  return w;
}

inline std::string word_text(const std::vector<mf::Symbol>& w, const std::string& alphabet) {
  std::string out;
  for (auto s : w) out += alphabet[s];
  return out;
}

// The generators of a request's group: [letters, reverses] pairs (or {letters, reverses}), each
// letters the images of the alphabet's letters in order.
inline std::vector<ma::Permutation> generators_of(const Json& group, const std::string& alphabet) {
  std::vector<ma::Permutation> out;
  if (group.is_null()) return out;
  if (!group.is_array()) throw Fault("a group is a list of generators");
  if (group.size() > 8) throw Fault("a group here has at most 8 generators");
  for (const auto& g : group.as_array()) {
    std::string letters;
    bool reverse = false;
    if (g.is_array()) {
      letters = g[0].text();
      reverse = g[1].is_bool() && g[1].as_bool();
    } else if (g.is_object()) {
      letters = g["letters"].text();
      reverse = flag(g, "reverses", false);
    } else if (g.is_string()) {
      letters = g.as_string();
    } else {
      throw Fault("a generator is [letters, reverses]");
    }
    if (letters.empty()) letters = alphabet;
    if (letters.size() != alphabet.size()) throw Fault("a generator, " + letters + ", must give one image for each letter of " + alphabet);
    std::vector<std::uint32_t> p;
    std::set<char> seen;
    for (char c : letters) {
      const auto at = alphabet.find(c);
      if (at == std::string::npos) throw Fault(std::string("a generator has a letter, ") + c + ", outside the alphabet");
      if (!seen.insert(c).second) throw Fault("a generator, " + letters + ", is not a permutation of " + alphabet);
      p.push_back(static_cast<std::uint32_t>(at));
    }
    out.push_back(encode(p, reverse));
  }
  return out;
}

// Whether an element carries the rules into themselves.
inline bool keeps(const mf::Grammar& g, const mf::GrammarAction& a, std::uint32_t element, const RuleSet& R) {
  for (const auto& r : g.rules)
    if (!R.count(a.act(g, element, r))) return false;
  return true;
}

// A parse tree as nested objects: a node is a rule over a span, a leaf a letter at a position.
inline Json tree_json(const mf::Grammar& g, const std::string& alphabet, const mf::ParseTree& t, const std::vector<mf::Symbol>& w, std::size_t offset = 0) {
  if (t.rule < 0) return Json::object({{"leaf", true}, {"at", t.begin + offset}, {"letter", std::string(1, alphabet[w[t.begin]])}});
  const auto& r = g.rules[static_cast<std::size_t>(t.rule)];
  Json children = Json::array();
  for (const auto& c : t.children) children.push(tree_json(g, alphabet, c, w, offset));
  return Json::object({{"rule", t.rule}, {"lhs", form_text(g, alphabet, r.lhs)}, {"text", rule_text(g, alphabet, r)}, {"begin", t.begin + offset},
                       {"end", t.end + offset}, {"children", children}});
}

struct Budget {};

// ---- groups named by a phrase ----------------------------------------------------------------------
inline std::string subscript(std::size_t n) {
  static const char* digits[] = {"₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"};
  const std::string s = std::to_string(n);
  std::string out;
  for (char c : s) out += digits[c - '0'];
  return out;
}

struct Named {
  ma::FiniteGroup group;
  std::string name;
  long symmetric = -1;  // n, when the group is S_n alone
};

// One factor: cyclic n, klein, dihedral n (order 2n), quaternion, symmetric n, alternating n; the
// short forms C4, V4, D4, Q8, S4, A5 too.
inline Named factor_group(std::string s) {
  std::string t;
  for (char c : s)
    if (!std::isspace(static_cast<unsigned char>(c))) t += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  std::size_t cut = t.size();
  while (cut > 0 && std::isdigit(static_cast<unsigned char>(t[cut - 1]))) cut--;
  const std::string name = t.substr(0, cut);
  const long n = cut < t.size() ? std::stol(t.substr(cut)) : -1;
  const auto need = [&](long lo, long hi, const char* what) {
    if (n < lo || n > hi) throw Fault(std::string(what) + " takes n from " + std::to_string(lo) + " to " + std::to_string(hi));
  };
  Named out;
  if (name == "klein" || (name == "v" && n == 4)) {
    out.group = ma::direct_product(ma::cyclic(2), ma::cyclic(2));
    out.name = "V₄";
  } else if (name == "cyclic" || name == "c" || name == "z") {
    need(1, 400, "a cyclic group");
    out.group = ma::cyclic(static_cast<std::size_t>(n));
    out.name = "C" + subscript(static_cast<std::size_t>(n));
  } else if (name == "dihedral" || name == "d") {
    need(3, 200, "a dihedral group (of order 2n)");
    out.group = ma::dihedral(static_cast<std::size_t>(n));
    out.name = "D" + subscript(static_cast<std::size_t>(n));
  } else if (name == "quaternion" || name == "q" || name == "quaternions") {
    if (n != -1 && n != 8) throw Fault("the quaternion group here is Q8");
    out.group = ma::quaternion();
    out.name = "Q₈";
  } else if (name == "symmetric" || name == "s") {
    need(1, 5, "a symmetric group");
    out.group = ma::symmetric(static_cast<std::size_t>(n));
    out.name = "S" + subscript(static_cast<std::size_t>(n));
    out.symmetric = n;
  } else if (name == "alternating" || name == "a") {
    need(3, 5, "an alternating group");
    out.group = ma::alternating(static_cast<std::size_t>(n));
    out.name = "A" + subscript(static_cast<std::size_t>(n));
  } else {
    throw Fault("no group \"" + s + "\": write cyclic n, klein, dihedral n, quaternion, symmetric n or alternating n, or products of them with ×");
  }
  return out;
}

inline Named named_group(const std::string& spec) {
  std::vector<std::string> parts;
  std::string cur;
  for (std::size_t i = 0; i < spec.size(); i++) {
    if (spec.compare(i, 2, "×") == 0) { parts.push_back(cur); cur.clear(); i += 1; continue; }
    if (spec[i] == '*' || ((spec[i] == 'x' || spec[i] == 'X') && i > 0 && i + 1 < spec.size() && std::isspace(static_cast<unsigned char>(spec[i - 1])) && std::isspace(static_cast<unsigned char>(spec[i + 1])))) {
      parts.push_back(cur);
      cur.clear();
      continue;
    }
    cur += spec[i];
  }
  parts.push_back(cur);
  if (parts.size() > 4) throw Fault("a product of at most four groups");
  Named out = factor_group(parts[0]);
  for (std::size_t k = 1; k < parts.size(); k++) {
    Named f = factor_group(parts[k]);
    if (out.group.order() * f.group.order() > 400) throw Fault("a group here has at most 400 elements");
    out.group = ma::direct_product(out.group, f.group);
    out.name += " × " + f.name;
    out.symmetric = -1;
  }
  if (out.group.order() > 400) throw Fault("a group here has at most 400 elements");
  return out;
}

// ---- the bases ---------------------------------------------------------------------------------------
inline const std::array<const char*, 4>& klein_names() {
  static const std::array<const char*, 4> n{"1", "τ", "κ", "σ"};  // klein()'s numbering
  return n;
}

inline std::string bases_of_mask(std::uint8_t m) {
  std::string s;
  for (seq::Code x = 0; x < 4; x++)
    if ((m >> x) & 1) s += seq::dna_letter(x);
  return s;
}

// The name of an element of V₄³ (the codon group) and of a character of it: per position, the
// involution (τ, κ, σ) or the classification whose character it is (RY, MK, WS), with the position.
inline std::string codon_element_name(std::uint32_t m) {
  static const std::array<const char*, 4> by_bits{"", "σ", "τ", "κ"};  // the 2-bit masks of nucleotides.hpp
  std::string s;
  for (int p = 0; p < 3; p++) {
    const auto b = (m >> (2 * (2 - p))) & 3u;
    if (b) s += std::string(by_bits[b]) + subscript(static_cast<std::size_t>(p + 1));
  }
  return s.empty() ? "1" : s;
}
inline std::string codon_character_name(std::uint32_t u) {
  static const std::array<const char*, 4> by_bits{"", "RY", "MK", "WS"};
  std::string s;
  for (int p = 0; p < 3; p++) {
    const auto b = (u >> (2 * (2 - p))) & 3u;
    if (b) s += (s.empty() ? "" : " ") + std::string(by_bits[b]) + subscript(static_cast<std::size_t>(p + 1));
  }
  return s.empty() ? "trivial" : s;
}

inline std::string strand_short(std::uint32_t e) {
  const std::uint32_t g = e / 2, r = e % 2;
  if (!r) return klein_names()[g];
  return g == 0 ? std::string("ρ") : std::string(klein_names()[g]) + "ρ";
}

}  // namespace grammar_ops

// ---- the operations ------------------------------------------------------------------------------------

// The Klein group of the bases (§3): its elements as permutations of A, C, G, T; its subgroups,
// whose cosets are IUPAC's codes; its characters, the three classifications, on the bases at a base
// point; and the chains of breakings from N to the bases.
BENCH_OP("grammar.nucleotides") {
  using namespace grammar_ops;
  const std::string bp_text = str(request, "base_point", "A");
  if (bp_text.size() != 1 || std::string("ACGT").find(bp_text[0]) == std::string::npos) throw Fault("the base point is A, C, G or T");
  const seq::Code bp = seq::dna_code(bp_text[0]);
  const auto& k = gg::klein();
  static const std::array<const char*, 4> about{"the identity", "transition: A↔G, C↔T, within the purines and within the pyrimidines",
                                                "Watson–Crick complementarity: A↔T, C↔G", "the third involution, κτ: A↔C, G↔T"};
  Json elements = Json::array(), product = Json::array(), torsor = Json::array();
  for (std::uint32_t g = 0; g < 4; g++) {
    std::string image;
    for (seq::Code x = 0; x < 4; x++) image += seq::dna_letter(static_cast<seq::Code>(k.bases.act(g, x)));
    elements.push(Json::object({{"name", klein_names()[g]}, {"about", about[g]}, {"letters", image}, {"xor", k.bases.act(g, 0)}}));
    Json row = Json::array();
    for (std::uint32_t h = 0; h < 4; h++) row.push(klein_names()[k.group.mul(g, h)]);
    product.push(row);
  }
  // The trivialisation at the base point: the base g·σ₀ for each element g.
  for (std::uint32_t g = 0; g < 4; g++) torsor.push(Json::object({{"element", klein_names()[g]}, {"base", std::string(1, seq::dna_letter(static_cast<seq::Code>(k.bases.act(g, bp))))}}));
  Json classes = Json::array();
  std::array<std::array<int, 4>, 3> values{};
  const auto cl = gg::classifications(bp);
  for (std::size_t c = 0; c < cl.size(); c++) {
    Json v = Json::array();
    for (seq::Code x = 0; x < 4; x++) {
      values[c][x] = gg::character_value(cl[c].kernel, bp, x);
      v.push(values[c][x]);
    }
    Json on_group = Json::array();
    for (std::uint32_t g = 0; g < 4; g++) on_group.push(g == 0 || g == cl[c].kernel ? 1 : -1);
    classes.push(Json::object({{"name", cl[c].name}, {"kernel", klein_names()[cl[c].kernel]},
                               {"letters", Json::array({std::string(1, cl[c].letters[0]), std::string(1, cl[c].letters[1])})},
                               {"cosets", Json::array({bases_of_mask(cl[c].masks[0]), bases_of_mask(cl[c].masks[1])})}, {"on_bases", v}, {"on_group", on_group}}));
  }
  bool product_rule = true;  // χ_RY χ_WS = χ_MK on the bases (Proposition 3.3)
  for (seq::Code x = 0; x < 4; x++) product_rule = product_rule && values[0][x] * values[1][x] == values[2][x];
  // morphe's character table of V₄, the classifications found in it by their kernels.
  const auto T = ma::character_table(k.group);
  Json table = Json::array();
  for (std::size_t chi = 0; chi < T.size(); chi++) {
    Json row = Json::array();
    // The kernel: the elements where χ is 1, V₄ itself for the trivial character.
    std::vector<std::uint32_t> ones;
    for (std::uint32_t g = 0; g < 4; g++) {
      row.push(T.value(chi, g).str());
      if (T.value(chi, g) == morphe::arithmetic::Cyclotomic(1)) ones.push_back(g);
    }
    const std::string kernel = ones.size() == 4 ? "V₄" : ones.size() == 2 ? klein_names()[ones[1]] : "?";
    std::string named = kernel == "V₄" ? "trivial" : "";
    for (const auto& c : cl)
      if (klein_names()[c.kernel] == kernel) named = c.name;
    table.push(Json::object({{"values", row}, {"kernel", kernel}, {"classification", named}}));
  }
  // The subgroups and their cosets, as IUPAC letters.
  const std::vector<std::pair<std::string, std::vector<std::uint32_t>>> subgroups{
      {"1", {0}}, {"⟨τ⟩", {0, k.transition}}, {"⟨κ⟩", {0, k.complement}}, {"⟨σ⟩", {0, k.transversion}}, {"V₄", {0, 1, 2, 3}}};
  Json subs = Json::array();
  std::set<int> all_cosets;
  for (const auto& [name, h] : subgroups) {
    Json cosets = Json::array();
    std::set<int> seen;
    for (seq::Code x = 0; x < 4; x++) {
      int m = 0;
      for (auto g : h) m |= 1 << k.bases.act(g, x);
      all_cosets.insert(m);
      if (!seen.insert(m).second) continue;
      cosets.push(Json::object({{"letter", std::string(1, seq::iupac_letter(static_cast<std::uint8_t>(m)))}, {"bases", bases_of_mask(static_cast<std::uint8_t>(m))}}));
    }
    Json names = Json::array();
    for (auto g : h) names.push(klein_names()[g]);
    subs.push(Json::object({{"name", name}, {"elements", names}, {"cosets", cosets}}));
  }
  bool three_letter_none = true;  // B, D, H, V are cosets of no subgroup (Proposition 3.2)
  for (char c : std::string("BDHV")) three_letter_none = three_letter_none && !all_cosets.count(seq::iupac_mask(c));
  // The chains of breakings N → H's cosets → the bases (grammar::quotient at each step).
  Json chains = Json::array();
  for (const auto& c : cl) {
    const auto mid = gg::quotient({0, c.kernel});
    const auto top = gg::quotient({0, 1, 2, 3});
    Json parent = Json::array(), mid_letters = Json::array();
    for (std::uint32_t q = 0; q < mid.to; q++) {
      for (seq::Code x = 0; x < 4; x++)
        if (mid.image[x] == q) { mid_letters.push(std::string(1, seq::iupac_letter(gg::coset_mask(c.kernel, x)))); break; }
    }
    for (seq::Code x = 0; x < 4; x++) parent.push(mid.image[x]);
    chains.push(Json::object({{"name", c.name}, {"kernel", klein_names()[c.kernel]}, {"top", top.to}, {"middle", mid_letters}, {"parent", parent},
                              {"residual", "V₄/⟨" + std::string(klein_names()[c.kernel]) + "⟩"}}));
  }
  return Json::object({{"base_point", bp_text}, {"elements", elements}, {"product", product}, {"torsor", torsor}, {"classifications", classes},
                       {"product_rule", product_rule}, {"table", table}, {"subgroups", subs}, {"three_letter_codes_are_cosets", !three_letter_none},
                       {"coset_count", all_cosets.size()}, {"chains", chains}});
}

// A finite group's exact character table (Dixon–Schneider, in cyclotomic numbers), its
// Frobenius–Schur indicators, and the Wedderburn blocks of ℂ[G] and ℚ[G], their dimensions measured
// by exact elimination (§7). For S_n alone, the degrees are held to the hook length formula.
BENCH_OP("grammar.characters") {
  using namespace grammar_ops;
  const auto t0 = std::chrono::steady_clock::now();
  const Named N = named_group(str(request, "group", "dihedral 4"));
  const auto& G = N.group;
  const auto T = ma::character_table(G);
  const double table_seconds = since(t0);
  // Class labels in the style of the ATLAS: the order of the elements, then a letter.
  Json classes = Json::array();
  std::map<std::size_t, int> seen_order;
  for (std::size_t j = 0; j < T.classes.size(); j++) {
    const std::size_t o = T.representative_orders[j];
    const int letter = seen_order[o]++;
    classes.push(Json::object({{"size", T.class_sizes[j]}, {"order", o}, {"label", std::to_string(o) + std::string(1, static_cast<char>('a' + letter % 26))}}));
  }
  Json characters = Json::array();
  std::vector<std::size_t> degrees;
  std::size_t squares = 0;
  for (std::size_t chi = 0; chi < T.size(); chi++) {
    Json values = Json::array(), complex = Json::array();
    bool real = true;
    for (std::size_t j = 0; j < T.classes.size(); j++) {
      values.push(T(chi, j).str());
      const auto z = T(chi, j).to_complex();
      complex.push(Json::array({std::round(z.real() * 1e9) / 1e9, std::round(z.imag() * 1e9) / 1e9}));
      real = real && T(chi, j) == T(chi, j).conj();
    }
    const std::size_t d = T.degree(chi);
    degrees.push_back(d);
    squares += d * d;
    const auto fs = ma::frobenius_schur(G, T, chi);
    characters.push(Json::object({{"degree", d}, {"values", values}, {"complex", complex}, {"real", real}, {"indicator", fs.to_double()}}));
  }
  // Row orthogonality, ⟨χ, ψ⟩ = δ_χψ, exactly.
  bool orthonormal = true;
  for (std::size_t a = 0; a < T.size(); a++)
    for (std::size_t b = 0; b < T.size(); b++)
      orthonormal = orthonormal && ma::inner_product(T, T.values[a], T.values[b]) == morphe::arithmetic::Cyclotomic(a == b ? 1 : 0);
  Json out = Json::object({{"name", N.name}, {"order", G.order()}, {"abelian", G.is_abelian()}, {"exponent", T.exponent}, {"prime", T.prime},
                           {"classes", classes}, {"characters", characters}, {"degrees", degrees}, {"sum_of_squares", squares}, {"orthonormal", orthonormal},
                           {"table_seconds", table_seconds}});
  // The blocks of ℂ[G] ≅ ⊕ M_χ(1)(ℂ), each the rank of x ↦ x e_χ, and of ℚ[G] (Galois orbits).
  const std::size_t cap = static_cast<std::size_t>(bounded(request, "measure_up_to", 64, 1, 128));
  Json blocks = Json::array();
  if (G.order() <= cap) {
    const auto t1 = std::chrono::steady_clock::now();
    const ma::GroupAlgebra<morphe::arithmetic::Cyclotomic> A(G);
    std::size_t total = 0;
    for (const auto& B : ma::wedderburn(A, T)) {
      blocks.push(Json::object({{"character", B.character}, {"degree", B.degree}, {"dimension", B.dimension}}));
      total += B.dimension;
    }
    const ma::GroupAlgebra<ma::Rational> Q(G);
    Json rational = Json::array();
    for (const auto& B : ma::rational_wedderburn(Q, T)) rational.push(Json::object({{"characters", B.characters}, {"dimension", B.dimension}}));
    out.set("wedderburn", Json::object({{"measured", true}, {"blocks", blocks}, {"total", total}, {"rational", rational}, {"seconds", since(t1)}}));
  } else {
    for (std::size_t chi = 0; chi < T.size(); chi++) blocks.push(Json::object({{"character", chi}, {"degree", degrees[chi]}, {"dimension", degrees[chi] * degrees[chi]}}));
    out.set("wedderburn", Json::object({{"measured", false}, {"blocks", blocks}, {"total", squares}}));
  }
  if (N.symmetric > 0) {
    Json parts = Json::array();
    std::vector<std::size_t> hooks;
    for (const auto& lambda : ma::partitions(static_cast<int>(N.symmetric))) {
      const auto d = static_cast<std::size_t>(ma::hook_length_dimension(lambda).to_int64());
      hooks.push_back(d);
      parts.push(Json::object({{"partition", std::vector<int>(lambda.begin(), lambda.end())}, {"dimension", d}}));
    }
    std::vector<std::size_t> a = degrees, b = hooks;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    out.set("symmetric", Json::object({{"n", N.symmetric}, {"partitions", parts}, {"degrees_agree", a == b}}));
  }
  return out;
}

// A grammar and a group acting on its letters (§1, §2): whether the grammar is invariant, the
// stabiliser of its rules and their closure; the rule orbits of the invariant grammar used, by
// Burnside's lemma; Earley's parse and its tree; the equivariant parser and the draft's orbit-only
// scan (on the subgroup that does not reverse words), with a survey of short words that shows the
// scan's language; spans shared by orbit and by content; and, for a monotone grammar that is not
// context-free, the search of its forms with and without the canonicaliser.
BENCH_OP("grammar.parse") {
  using namespace grammar_ops;
  const auto t0 = std::chrono::steady_clock::now();
  const Written wr = written_grammar(request);
  const mf::Grammar& given = wr.grammar;
  const std::string& alphabet = wr.alphabet;
  const std::size_t t = given.terminals;
  const auto G = acting(generators_of(request["group"], alphabet), t);
  const RuleSet R(given.rules.begin(), given.rules.end());
  std::vector<std::uint32_t> stab;
  Json keeps_json = Json::array();
  for (std::uint32_t g = 0; g < G->group.order(); g++) {
    const bool k = keeps(given, G->on_grammar, g, R);
    if (k) stab.push_back(g);
    keeps_json.push(k);
  }
  const bool invariant = stab.size() == G->group.order();
  const mf::Grammar closed = mf::closure(given, G->on_grammar);
  const std::string use = str(request, "use", "stabiliser");
  const bool on_closure = !invariant && use == "closure";
  if (on_closure && closed.rules.size() > 3000) throw Fault("the closure has " + std::to_string(closed.rules.size()) + " rules, more than this takes");
  std::unique_ptr<Acting> S = invariant || on_closure ? nullptr : subgroup(*G, stab);
  const mf::Grammar& used = on_closure ? closed : given;
  const Acting& A = S ? *S : *G;

  Json rules = Json::array();
  for (const auto& r : given.rules) rules.push(rule_text(given, alphabet, r));
  Json added = Json::array();
  for (const auto& r : closed.rules)
    if (!R.count(r) && added.size() < 200) added.push(rule_text(closed, alphabet, r));
  Json stab_names = Json::array();
  for (auto g : stab) stab_names.push(element_name(*G, g, alphabet));
  Json out = Json::object({{"alphabet", alphabet}, {"nonterminals", given.nonterminals}, {"rules", rules}, {"start", given.nonterminals[0]},
                           {"context_free", given.context_free()}, {"monotone", given.monotone()}, {"group", group_json(*G, alphabet)},
                           {"invariant", invariant}, {"keeps", keeps_json}, {"stabiliser", stab_names},
                           {"closure", Json::object({{"rules", closed.rules.size()}, {"added", added}, {"added_count", closed.rules.size() - given.rules.size()}})},
                           {"used", on_closure ? "closure" : "given"}, {"acting", group_json(A, alphabet)}});

  // The rule orbits of the grammar used, under the group that keeps it (Proposition 1.4).
  const auto orbits = mf::rule_orbits(used, A.on_grammar);
  Json orbit_list = Json::array();
  for (std::uint32_t o = 0; o < orbits.representatives.size() && o < 200; o++) {
    Json members = Json::array();
    for (std::uint32_t r = 0; r < used.rules.size(); r++)
      if (orbits.orbit_of[r] == o && members.size() < 16) members.push(rule_text(used, alphabet, used.rules[r]));
    orbit_list.push(Json::object({{"representative", rule_text(used, alphabet, used.rules[orbits.representatives[o]])}, {"size", orbits.size[o]}, {"members", members}}));
  }
  out.set("orbits", Json::object({{"count", orbits.representatives.size()}, {"burnside", orbits.burnside}, {"rules", used.rules.size()},
                                  {"group_order", A.group.order()}, {"ratio", static_cast<double>(used.rules.size()) / static_cast<double>(A.group.order())},
                                  {"list", orbit_list}}));

  const std::vector<mf::Symbol> w = word_of(str(request, "word"), alphabet);
  const bool cf = used.context_free();
  const bool monotone_only = !cf && used.monotone();
  if (cf && w.size() > 64) throw Fault("a word here has at most 64 letters");
  if (monotone_only && w.size() > 14) throw Fault("the search of a monotone grammar's forms is exponential: a word of at most 14 letters");
  out.set("word", word_text(w, alphabet));
  {
    std::set<std::vector<mf::Symbol>> images;
    for (std::uint32_t g = 0; g < A.group.order(); g++) images.insert(A.action.act_on_word(g, w));
    out.set("word_orbit", Json::object({{"canonical", word_text(A.action.canonical(w), alphabet)}, {"size", images.size()}}));
  }

  if (cf) {
    const auto t1 = std::chrono::steady_clock::now();
    const mf::Earley earley(used);
    const auto parsed = earley.parse(w, true);
    Json p = Json::object({{"accepts", parsed.accepts}, {"items", parsed.items}, {"seconds", since(t1)}});
    if (parsed.accepts) {
      p.set("tree", tree_json(used, alphabet, parsed.tree, w));
      p.set("yield_checked", mf::yield(used, parsed.tree, w) == w);
    }
    if (on_closure) p.set("given_accepts", mf::Earley(given).accepts(w));
    out.set("earley", p);

    // The equivariant parser binds each rule's element; the orbit-only scan does not (§2.1–2.2).
    // Both need an action by letters alone: the subgroup of the elements that do not reverse.
    std::vector<std::uint32_t> letters_only;
    for (std::uint32_t g = 0; g < A.group.order(); g++)
      if (!A.action.reverses(g)) letters_only.push_back(g);
    std::unique_ptr<Acting> Hown = letters_only.size() == A.group.order() ? nullptr : subgroup(A, letters_only);
    const Acting& H = Hown ? *Hown : A;
    mf::EquivariantEarley equivariant(used, H.on_grammar, true), scan(used, H.on_grammar, false);
    const auto t2 = std::chrono::steady_clock::now();
    const auto e = equivariant.parse(w);
    const auto s = scan.parse(w);
    Json hn = Json::array();
    for (std::uint32_t g = 0; g < H.group.order(); g++) hn.push(element_name(H, g, alphabet));
    Json eq = Json::object({{"group", hn}, {"order", H.group.order()}, {"orbits", equivariant.orbits().representatives.size()},
                            {"consistent", Json::object({{"accepts", e.accepts}, {"items", e.items}, {"plain_equivalent", e.plain_equivalent}})},
                            {"scan", Json::object({{"accepts", s.accepts}, {"items", s.items}})},
                            {"chart_agrees", e.plain_equivalent == parsed.items}, {"seconds", since(t2)}});
    // Short words, all of them: where the scan accepts outside the language (Proposition 2.1), and
    // whether the equivariant parser agrees with Earley's on every one (Proposition 2.3).
    if (H.group.order() > 1) {
      std::size_t length = 0, total = 0, power = 1;
      while (length < 10) {
        power *= t;
        if (total + power > 4000) break;
        total += power;
        length++;
      }
      std::size_t in_language = 0, by_scan = 0, outside = 0;
      bool agree = true;
      Json examples = Json::array();
      std::vector<mf::Symbol> x;
      for (std::size_t n = 1; n <= length; n++) {
        x.assign(n, 0);
        for (;;) {
          const bool a = earley.accepts(x), b = scan.parse(x).accepts;
          agree = agree && equivariant.parse(x).accepts == a;
          in_language += a;
          by_scan += b;
          if (b && !a) {
            outside++;
            if (examples.size() < 8) examples.push(word_text(x, alphabet));
          }
          std::size_t i = n;
          while (i > 0 && ++x[i - 1] == t) x[--i] = 0;
          if (i == 0) break;
        }
      }
      eq.set("survey", Json::object({{"up_to", length}, {"words", total}, {"in_language", in_language}, {"by_scan", by_scan}, {"outside", outside},
                                     {"examples", examples}, {"equivariant_agrees", agree}}));
    }
    out.set("equivariant", eq);

    // Spans decided once per class (§2.3): under the group, and under content alone.
    try {
      const auto trivial = acting({}, t);
      const auto t3 = std::chrono::steady_clock::now();
      const auto by_group = mf::span_recognition(used, A.on_grammar, w);
      const auto by_content = mf::span_recognition(used, trivial->on_grammar, w);
      const auto pack = [](const mf::SpanRecognition& r) {
        return Json::object({{"accepts", r.accepts}, {"spans", r.spans}, {"reused", r.reused}, {"classes", r.classes}});
      };
      out.set("spans", Json::object({{"by_group", pack(by_group)}, {"by_content", pack(by_content)}, {"seconds", since(t3)}}));
    } catch (const std::invalid_argument& e) {
      out.set("spans", Json::object({{"skipped", detail::message_of(e)}}));
    }
  } else if (monotone_only) {
    // Membership by searching the forms no longer than the word (§2.4), with and without the
    // canonicaliser that replaces each form by the least of its images.
    const std::size_t budget = 2'000'000;
    const auto canon = mf::form_canonicaliser(used, A.on_grammar);
    const auto search = [&](bool classes) {
      std::size_t calls = 0;
      const std::function<mf::Form(const mf::Form&)> counted = [&](const mf::Form& f) {
        if (++calls > budget) throw Budget{};
        return classes ? canon(f) : f;
      };
      const auto t4 = std::chrono::steady_clock::now();
      try {
        const auto r = mf::derives_monotone(used, w, counted);
        return Json::object({{"derives", r.derives}, {"visited", r.visited}, {"seconds", since(t4)}});
      } catch (const Budget&) {
        return Json::object({{"stopped", true}, {"budget", budget}, {"seconds", since(t4)}});
      }
    };
    out.set("monotone_search", Json::object({{"classes", search(true)}, {"forms", search(false)}}));
  } else {
    out.set("note", "the grammar is neither context-free nor monotone: membership is undecidable in general, and nothing is parsed");
  }
  out.set("seconds", since(t0));
  return out;
}

// Every hairpin of a short sequence (§5): the spans the hairpin grammar's start derives, from one
// inside pass over every span no longer than max_span (Boolean semiring), held to the definition
// by brute force; the maximal ones with their stems and loops, and the Viterbi parse (the most
// pairs) of the longest; the grammar's symmetry under the eight symmetries of words; and graduated
// breaking (§2.5): the grammar of a hairpin anywhere parsed first over R/Y, W/S or M/K, the spans
// each keeps open, and that the fine pass over them finds every hairpin.
BENCH_OP("grammar.hairpins") {
  using namespace grammar_ops;
  const auto t0 = std::chrono::steady_clock::now();
  detail::Bounds b;
  b.bytes = std::size_t{64} << 20;
  detail::Reading reading;
  const auto first = first_record(request, "fasta", b, reading);
  require_dna(first.record, "a hairpin search");
  const seq::Codes all = seq::encode_dna(first.record.sequence);
  const auto begin = static_cast<std::size_t>(bounded(request, "begin", 0, 0, static_cast<long long>(all.size())));
  const auto end = static_cast<std::size_t>(bounded(request, "end", static_cast<long long>(all.size()), static_cast<long long>(begin), static_cast<long long>(all.size())));
  if (end - begin > 400) throw Fault(first.record.id + " has " + std::to_string(end - begin) + " bases; hairpins are searched here in up to 400 (give begin and end, or a stretch such as a tRNA gene)");
  if (end == begin) throw Fault("no bases to search");
  gg::HairpinSpec spec;
  spec.min_stem = static_cast<std::size_t>(bounded(request, "min_stem", 4, 2, 12));
  spec.min_loop = static_cast<std::size_t>(bounded(request, "min_loop", 3, 1, 12));
  spec.max_loop = static_cast<std::size_t>(bounded(request, "max_loop", 8, static_cast<long long>(spec.min_loop), 20));
  spec.wobble = flag(request, "wobble", false);
  const auto max_span = static_cast<std::size_t>(bounded(request, "max_span", 80, 10, 200));
  const std::vector<mf::Symbol> w(all.begin() + static_cast<std::ptrdiff_t>(begin), all.begin() + static_cast<std::ptrdiff_t>(end));
  const std::size_t n = w.size();
  const std::string letters = seq::decode_dna(seq::Codes(w.begin(), w.end()));

  const mf::Grammar g = gg::hairpin_grammar(spec);
  const auto bin = mf::BinaryGrammar::from(g);
  const std::vector<bool> ones(g.rules.size(), true);
  std::vector<char> mask((n + 1) * (n + 1), 0);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 1; j <= n && j - i <= max_span; j++) mask[i * (n + 1) + j] = 1;
  const auto t1 = std::chrono::steady_clock::now();
  const auto in = mf::inside<mf::Boolean>(bin, ones, w, &mask);
  const double inside_seconds = since(t1);
  std::vector<std::pair<std::size_t, std::size_t>> spans;
  for (std::size_t len = 1; len <= std::min(n, max_span); len++)
    for (std::size_t i = 0; i + len <= n; i++)
      if (in.at(0, i, i + len)) spans.emplace_back(i, i + len);
  // The definition, by brute force over the same spans.
  std::size_t by_definition = 0;
  bool agrees = true;
  for (std::size_t len = 1; len <= std::min(n, max_span); len++)
    for (std::size_t i = 0; i + len <= n; i++) {
      const bool h = gg::is_hairpin(std::vector<mf::Symbol>(w.begin() + static_cast<std::ptrdiff_t>(i), w.begin() + static_cast<std::ptrdiff_t>(i + len)), spec);
      by_definition += h;
      agrees = agrees && h == static_cast<bool>(in.at(0, i, i + len));
    }
  // The maximal spans, longest first: those inside no other.
  std::vector<std::pair<std::size_t, std::size_t>> maximal;
  for (auto it = spans.rbegin(); it != spans.rend(); ++it) {
    bool inside_other = false;
    for (const auto& m : maximal) inside_other = inside_other || (m.first <= it->first && it->second <= m.second);
    if (!inside_other) maximal.push_back(*it);
  }
  std::sort(maximal.begin(), maximal.end());
  const auto pairs_with = [&](mf::Symbol a, mf::Symbol c) {
    for (const auto& [x, y] : gg::allowed_pairs(spec.wobble))
      if (x == a && y == c) return true;
    return false;
  };
  Json found = Json::array();
  std::size_t longest = 0, longest_at = 0;
  for (std::size_t m = 0; m < maximal.size() && m < 300; m++) {
    const auto [i, j] = maximal[m];
    const std::size_t len = j - i;
    // The longest stem with a loop in range: the parse with the most pairs.
    std::size_t stem = 0;
    for (std::size_t s = spec.min_stem; 2 * s + spec.min_loop <= len; s++) {
      bool ok = true;
      for (std::size_t q = 0; q < s && ok; q++) ok = pairs_with(w[i + q], w[j - 1 - q]);
      if (!ok) break;
      const std::size_t loop = len - 2 * s;
      if (loop >= spec.min_loop && loop <= spec.max_loop) stem = s;
    }
    std::string dots(len, '.');
    Json pairs = Json::array();
    for (std::size_t q = 0; q < stem; q++) {
      dots[q] = '(';
      dots[len - 1 - q] = ')';
      pairs.push(Json::array({i + q + begin, j - 1 - q + begin}));
    }
    found.push(Json::object({{"begin", i + begin}, {"end", j + begin}, {"length", len}, {"stem", stem}, {"loop", len - 2 * stem},
                             {"structure", dots}, {"letters", letters.substr(i, len)}, {"pairs", pairs}}));
    if (len > longest) { longest = len; longest_at = m; }
  }
  Json out = Json::object({{"name", first.record.id}, {"length", all.size()}, {"begin", begin}, {"end", end}, {"letters", letters},
                           {"spec", Json::object({{"min_stem", spec.min_stem}, {"min_loop", spec.min_loop}, {"max_loop", spec.max_loop}, {"wobble", spec.wobble}, {"max_span", max_span}})},
                           {"spans", spans.size()}, {"by_definition", by_definition}, {"agrees_with_definition", agrees}, {"maximal", found},
                           {"maximal_count", maximal.size()}, {"inside_seconds", inside_seconds}, {"reading", reading.json()}});
  // The Viterbi parse of the longest: stem rules weigh 1, so the best parse has the most pairs.
  if (!maximal.empty()) {
    const auto [i, j] = maximal[longest_at];
    const std::vector<mf::Symbol> sub(w.begin() + static_cast<std::ptrdiff_t>(i), w.begin() + static_cast<std::ptrdiff_t>(j));
    std::vector<long long> weights(g.rules.size(), 0);
    for (std::size_t r = 0; r < g.rules.size(); r++)
      if (g.index(g.rules[r].lhs[0]) < spec.min_stem) weights[r] = 1;
    long long best = 0;
    const auto tree = mf::best_parse<mf::MaxPlus<long long>>(g, bin, weights, sub, &best);
    out.set("longest", longest_at);
    out.set("tree", tree_json(g, "ACGT", tree, sub, i + begin));
    out.set("tree_pairs", best);
  }

  // The grammar's symmetry under the eight symmetries of words: which keep its rules, the subgroup
  // they form, and the rule orbits under it (§5, Proposition 5.1).
  {
    const auto& w8 = gg::base_and_strand_group();
    const mf::GrammarAction on8{&w8.words, {}};
    const RuleSet R(g.rules.begin(), g.rules.end());
    Json by_element = Json::array();
    std::vector<ma::Permutation> kept;
    std::vector<std::string> kept_names;
    for (std::uint32_t e = 0; e < 8; e++) {
      const bool k = keeps(g, on8, e, R);
      by_element.push(Json::object({{"name", strand_short(e)}, {"keeps", k}}));
      if (k) {
        std::vector<std::uint32_t> p;
        for (std::uint32_t x = 0; x < 4; x++) p.push_back(w8.words.act(e, x));
        kept.push_back(encode(p, w8.words.reverses(e)));
        kept_names.push_back(strand_short(e));
      }
    }
    const auto S = acting(kept, 4, 8);
    const auto orbits = mf::rule_orbits(g, S->on_grammar);
    bool strand = true;
    for (std::size_t e : {std::size_t{0}, std::size_t{1}, std::size_t{4}, std::size_t{5}}) strand = strand && by_element[e]["keeps"].as_bool();
    Json rules = Json::array();
    for (std::size_t r = 0; r < g.rules.size() && r < 400; r++) rules.push(rule_text(g, "ACGT", g.rules[r]));
    out.set("grammar", Json::object({{"rules", g.rules.size()}, {"nonterminals", g.nonterminals.size()}, {"rule_texts", rules}, {"by_element", by_element},
                                     {"stabiliser", kept_names}, {"strand_invariant", strand},
                                     {"orbits", Json::object({{"count", orbits.representatives.size()}, {"burnside", orbits.burnside}, {"group_order", S->group.order()}})}}));
  }

  // Graduated breaking on a window of at most 120 bases (the whole when it is short, else centred
  // on the longest hairpin), with the grammar of a hairpin anywhere.
  {
    const std::size_t W = std::min<std::size_t>(n, 120);
    std::size_t at = 0;
    if (n > W && !maximal.empty()) {
      const auto [i, j] = maximal[longest_at];
      const std::size_t mid = (i + j) / 2;
      at = mid > W / 2 ? std::min(mid - W / 2, n - W) : 0;
    }
    const std::vector<mf::Symbol> ww(w.begin() + static_cast<std::ptrdiff_t>(at), w.begin() + static_cast<std::ptrdiff_t>(at + W));
    if (std::any_of(ww.begin(), ww.end(), [](mf::Symbol x) { return x > 3; })) {
      out.set("coarse", Json::object({{"skipped", "the window holds a base that is not A, C, G or T"}}));
    } else {
      const auto t2 = std::chrono::steady_clock::now();
      const mf::Grammar fine = gg::contains_hairpin_grammar(spec);
      const auto bfine = mf::BinaryGrammar::from(fine);
      const std::vector<bool> fine_ones(fine.rules.size(), true);
      const auto full = mf::inside<mf::Boolean>(bfine, fine_ones, ww);
      const std::size_t hp = 2;  // the hairpin's start, after Region and Flank
      const auto hairpin_spans = [&](const mf::Chart<mf::Boolean>& c) {
        std::vector<std::pair<std::size_t, std::size_t>> s;
        for (std::size_t i = 0; i < W; i++)
          for (std::size_t j = i + 1; j <= W; j++)
            if (c.at(hp, i, j)) s.emplace_back(i, j);
        return s;
      };
      const auto reference = hairpin_spans(full);
      const auto& k = gg::klein();
      Json by = Json::array();
      bool exact = true;
      std::size_t total = 0;
      const std::array<std::pair<const char*, std::uint32_t>, 3> levels{{{"R/Y", k.transition}, {"W/S", k.complement}, {"M/K", k.transversion}}};
      for (const auto& [name, h] : levels) {
        mf::CoarseToFine st;
        const auto cm = mf::coarse_mask(fine, gg::quotient({0, h}), ww, &st);
        const auto pruned = mf::inside<mf::Boolean>(bfine, fine_ones, ww, &cm);
        const bool same = static_cast<bool>(pruned.at(0, 0, W)) == static_cast<bool>(full.at(0, 0, W)) && hairpin_spans(pruned) == reference;
        exact = exact && same;
        total = st.spans;
        by.push(Json::object({{"name", name}, {"kept", st.kept}, {"spans", st.spans}, {"fraction", st.spans ? static_cast<double>(st.kept) / static_cast<double>(st.spans) : 0.0},
                              {"exact", same}}));
      }
      out.set("coarse", Json::object({{"begin", at + begin}, {"end", at + begin + W}, {"spans", total}, {"hairpins", reference.size()}, {"by", by},
                                      {"exact", exact}, {"seconds", since(t2)}}));
    }
  }
  out.set("seconds", since(t0));
  return out;
}

// The symmetries a genome's k-mer statistics keep (§4): the total variation between the frequencies
// and their image under each of the eight symmetries of words, and the eight isotypic channels of the
// frequency vector under that group (§7), with the canonical k-mers counted by Burnside's lemma.
BENCH_OP("grammar.symmetry") {
  using namespace grammar_ops;
  const auto t0 = std::chrono::steady_clock::now();
  const auto k = static_cast<std::size_t>(bounded(request, "k", 4, 1, 6));
  seq::Codes x;
  std::string name;
  Json reading_json;
  std::size_t records = 0;
  if (request.has("genome")) {
    const auto g = worker.genome(request);
    for (const auto& c : g->contigs) {
      if (!x.empty()) x.push_back(seq::N);  // no k-mer spans two contigs
      x.insert(x.end(), c.codes.begin(), c.codes.end());
      records++;
    }
    name = g->title;
  } else {
    detail::Bounds b;
    b.bytes = std::size_t{256} << 20;
    b.whole = false;
    const auto reading = each_record(request, "fasta", b, [&](seq::Record&& r) {
      require_dna(r, "the symmetry of k-mers");
      if (!x.empty()) x.push_back(seq::N);
      const auto codes = seq::encode_dna(r.sequence);
      x.insert(x.end(), codes.begin(), codes.end());
      if (name.empty()) name = r.id;
      records++;
    });
    reading_json = reading.json();
    if (records > 1) name += " and " + std::to_string(records - 1) + " more";
  }
  std::size_t acgt = 0;
  for (auto c : x) acgt += c < 4;
  if (acgt < k) throw Fault("too few bases for k-mers of " + std::to_string(k));
  const auto& w8 = gg::base_and_strand_group();
  const auto inv = gg::kmer_invariance(x, k, w8);
  Json elements = Json::array();
  for (const auto& v : inv)
    elements.push(Json::object({{"element", v.element}, {"name", v.name}, {"short", strand_short(v.element)}, {"reverses", w8.words.reverses(v.element)},
                                {"total_variation", v.total_variation}}));
  double norm = 0, sum = 0, odd_rc = 0, odd_complement = 0, odd_reverse = 0;
  const auto channels = gg::kmer_channels(x, k, w8, &norm);
  Json ch = Json::array();
  for (const auto& c : channels) {
    sum += c.energy;
    Json o = Json::object({{"character", c.character}, {"signs", c.signs}, {"energy", c.energy}, {"fraction", c.fraction}});
    if (!c.signs.empty()) {
      o.set("kappa", c.signs[w8.complement]);
      o.set("tau", c.signs[2 * gg::klein().transition]);
      o.set("rho", c.signs[w8.reverse]);
      o.set("rc", c.signs[w8.reverse_complement]);
      if (c.signs[w8.reverse_complement] < 0) odd_rc += c.fraction;
      if (c.signs[w8.complement] < 0) odd_complement += c.fraction;
      if (c.signs[w8.reverse] < 0) odd_reverse += c.fraction;
    }
    ch.push(o);
  }
  Json names = Json::array();
  for (std::uint32_t e = 0; e < 8; e++) names.push(strand_short(e));
  const auto orbits = [&](const gg::WordGroup& g) { return g.words.word_orbits(k).str(); };
  Json out = Json::object({{"name", name}, {"records", records}, {"bases", acgt}, {"k", k}, {"kmers", std::size_t{1} << (2 * k)}, {"element_names", names},
                           {"elements", elements}, {"channels", ch}, {"norm", norm}, {"channel_sum", sum}, {"odd_rc", odd_rc}, {"odd_complement", odd_complement},
                           {"odd_reverse", odd_reverse},
                           {"canonical", Json::object({{"reverse_complement", orbits(gg::reverse_complement_group())}, {"strand", orbits(gg::strand_group())},
                                                       {"all_eight", orbits(w8)}})},
                           {"seconds", since(t0)}});
  if (!reading_json.is_null()) out.set("reading", reading_json);
  return out;
}

// A genetic code of NCBI's as a function on (C₂)⁶ (§6): the 64 codons as six-bit numbers, the
// invariance profile I(m), the stabiliser, the codons where a chosen element fails, the power
// spectrum P(u) with the identity P = Ŵ I checked exactly, and the share of the spectrum in the
// characters that ignore the third base, and the second and third.
BENCH_OP("grammar.code") {
  using namespace grammar_ops;
  const auto& code = genetic_code_of(request, 1);
  const auto aa = gg::code_on_group(code);
  const auto profile = gg::invariance_profile(aa);
  const auto stab = gg::stabiliser(aa);
  const auto power = gg::code_power_spectrum(aa);
  const std::uint32_t tau3 = gg::at_position(gg::klein().transition, 2);
  const auto m = static_cast<std::uint32_t>(bounded(request, "element", tau3, 0, 63));
  Json codons = Json::array();
  for (std::uint32_t c = 0; c < 64; c++) {
    std::string letters;
    for (auto b : gg::codon_bases(c)) letters += seq::dna_letter(b);
    codons.push(Json::object({{"bits", c}, {"codon", letters}, {"aa", std::string(1, aa[c])}}));
  }
  Json elements = Json::array(), characters = Json::array();
  for (std::uint32_t u = 0; u < 64; u++) {
    elements.push(Json::object({{"bits", u}, {"name", codon_element_name(u)}, {"keeps", profile[u]}}));
    characters.push(Json::object({{"bits", u}, {"name", codon_character_name(u)}, {"power", power[u]}}));
  }
  Json stab_names = Json::array();
  for (auto s : stab) stab_names.push(Json::object({{"bits", s}, {"name", codon_element_name(s)}}));
  // P = Ŵ I exactly (the Wiener–Khinchin identity on the group), and Proposition 6.2: m keeps the
  // code exactly when P vanishes on every character u with u·m odd.
  const auto wh = ma::walsh_hadamard(std::vector<long long>(profile.begin(), profile.end()));
  bool wiener_khinchin = true;
  for (std::size_t u = 0; u < 64; u++) wiener_khinchin = wiener_khinchin && wh[u] == power[u];
  bool duality = true;
  for (std::uint32_t e = 0; e < 64; e++) {
    bool vanishes = true;
    for (std::uint32_t u = 0; u < 64; u++)
      if (__builtin_popcount(u & e) % 2) vanishes = vanishes && power[u] == 0;
    duality = duality && (profile[e] == 64) == vanishes;
  }
  long long all = 0, first_two = 0, first = 0, third = 0;
  for (std::uint32_t u = 1; u < 64; u++) {
    all += power[u];
    if (gg::ignores_positions(u, {2})) first_two += power[u];
    if (gg::ignores_positions(u, {1, 2})) first += power[u];
    if (gg::ignores_positions(u, {0, 1})) third += power[u];
  }
  const auto share = [&](long long x) { return all ? static_cast<double>(x) / static_cast<double>(all) : 0.0; };
  return Json::object({{"table", code.id}, {"name", std::string(code.name)}, {"codons", codons}, {"elements", elements}, {"characters", characters},
                       {"stabiliser", stab_names}, {"element", Json::object({{"bits", m}, {"name", codon_element_name(m)}, {"keeps", profile[m]}, {"failures", gg::failures(aa, m)}})},
                       {"tau3", tau3}, {"wiener_khinchin", wiener_khinchin}, {"duality", duality},
                       {"energy", Json::object({{"total", all}, {"first_two", share(first_two)}, {"first", share(first)}, {"third", share(third)}})}});
}

}  // namespace workbench
