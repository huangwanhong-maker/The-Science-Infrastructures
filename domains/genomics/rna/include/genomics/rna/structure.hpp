// Secondary structures of RNA: the pairs of a sequence as a pair table, their dot-bracket
// notation, the rules a structure obeys (pairs not crossing, a minimum hairpin loop, the pairs
// the chemistry allows), the number of structures of a length, and their enumeration for checks
// (theory/rna-structure.md §1).
//
// A secondary structure on the positions 0, …, n − 1 is a set P of pairs (i, j), i < j, each
// position in at most one pair, no two pairs crossing: for (i, j), (k, l) ∈ P with i < k, either
// l < j (nested) or j < k (side by side) (Nussinov, Pieczenik, Griggs and Kleitman 1978;
// Durbin, Eddy, Krogh and Mitchison 1998, ch. 10). Its pair table gives each position its
// partner or −1. A hairpin loop is the interior of a pair that encloses no pair; a structure has
// minimum loop m when j − i − 1 ≥ m for every pair. Dot-bracket writes a pair as ( … ) and an
// unpaired base as '.'; a string of balanced brackets is exactly a structure without crossings.
//
// The canonical pairs are the Watson–Crick pairs A·U and G·C, both ways round, and the wobble
// pair G·U, both ways round (the domain's codes: A = 0, C = 1, G = 2, U = T = 3). Whether a pair
// is canonical is chemistry, stated here once; the grammars of scfg.hpp allow every pair and
// let their probabilities say which occur.
//
// The number of structures. With every pair allowed and minimum loop m, the structures on n
// positions number S_m(0) = 1 and S_m(n) = S_m(n − 1) + Σ_{k = 0}^{n − m − 2} S_m(k) S_m(n − k − 2)
// (the last position unpaired, or paired with position k, its interior of length n − k − 2 ≥ m
// and the prefix of length k independent of it); exact, in morphe's integers. Enumeration visits
// every structure, at least S_m(n) of them, so it serves only short sequences.

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <morphe/arithmetic/integer.hpp>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace genomics::rna {

using sequence::Code;
using sequence::Codes;

// A secondary structure as its pair table: partner[i] = j when (i, j) is a pair, −1 when i is
// unpaired.
struct Structure {
  std::vector<std::int32_t> partner;

  Structure() = default;
  explicit Structure(std::size_t n) : partner(n, -1) {}
  std::size_t size() const { return partner.size(); }
  bool paired(std::size_t i) const { return partner[i] >= 0; }
  void add(std::size_t i, std::size_t j) {
    if (i >= j || j >= size()) throw std::invalid_argument("genomics::rna::Structure::add: a pair (i, j) needs i < j < n");
    if (partner[i] >= 0 || partner[j] >= 0) throw std::invalid_argument("genomics::rna::Structure::add: a position already paired");
    partner[i] = static_cast<std::int32_t>(j);
    partner[j] = static_cast<std::int32_t>(i);
  }
  // The pairs (i, j), i < j, in increasing order of i.
  std::vector<std::pair<std::size_t, std::size_t>> pairs() const {
    std::vector<std::pair<std::size_t, std::size_t>> out;
    for (std::size_t i = 0; i < size(); i++)
      if (partner[i] > static_cast<std::int32_t>(i)) out.emplace_back(i, static_cast<std::size_t>(partner[i]));
    return out;
  }
  std::size_t pair_count() const {
    std::size_t c = 0;
    for (std::size_t i = 0; i < size(); i++) c += partner[i] > static_cast<std::int32_t>(i);
    return c;
  }
  bool has(std::size_t i, std::size_t j) const { return i < size() && partner[i] == static_cast<std::int32_t>(j); }
  friend bool operator==(const Structure&, const Structure&) = default;
};

// ---- canonical pairs ----------------------------------------------------------------------------
inline bool watson_crick(Code a, Code b) { return a < 4 && b < 4 && a + b == 3; }  // A·U, U·A, C·G, G·C
inline bool wobble(Code a, Code b) { return (a == sequence::G && b == sequence::T) || (a == sequence::T && b == sequence::G); }
inline bool canonical(Code a, Code b, bool with_wobble = true) { return watson_crick(a, b) || (with_wobble && wobble(a, b)); }

// ---- validity -----------------------------------------------------------------------------------
// The first rule the structure breaks, or "" when it is a secondary structure of the sequence
// with minimum loop min_loop (and, when allowed is given, only allowed pairs).
inline std::string violation(const Structure& s, std::size_t min_loop = 0, const Codes* x = nullptr, bool with_wobble = true) {
  const std::size_t n = s.size();
  if (x && x->size() != n) return "the structure and the sequence differ in length";
  for (std::size_t i = 0; i < n; i++) {
    const std::int32_t j = s.partner[i];
    if (j < 0) continue;
    if (static_cast<std::size_t>(j) >= n || j == static_cast<std::int32_t>(i)) return "a partner outside the sequence";
    if (s.partner[static_cast<std::size_t>(j)] != static_cast<std::int32_t>(i)) return "a pair not recorded at both ends";
  }
  // Crossings: a stack of open pairs, scanning left to right.
  std::vector<std::size_t> open;
  for (std::size_t i = 0; i < n; i++) {
    const std::int32_t j = s.partner[i];
    if (j < 0) continue;
    if (static_cast<std::size_t>(j) > i) {
      if (static_cast<std::size_t>(j) - i - 1 < min_loop) return "a loop shorter than the minimum";
      if (x && !canonical((*x)[i], (*x)[static_cast<std::size_t>(j)], with_wobble)) return "a pair that is not canonical";
      open.push_back(i);
    } else {
      if (open.empty() || open.back() != static_cast<std::size_t>(j)) return "two pairs cross";
      open.pop_back();
    }
  }
  return "";
}

// ---- dot-bracket ----------------------------------------------------------------------------------
inline std::string dot_bracket(const Structure& s) {
  std::string out(s.size(), '.');
  for (std::size_t i = 0; i < s.size(); i++)
    if (s.partner[i] >= 0) out[i] = s.partner[i] > static_cast<std::int32_t>(i) ? '(' : ')';
  return out;
}

// Reads '(' ')' and '.' (also '-', ',', ':', '_' and '~' as unpaired, as several tools write
// them); throws std::invalid_argument on an unbalanced bracket or another character.
inline Structure parse_dot_bracket(std::string_view text) {
  Structure s(text.size());
  std::vector<std::size_t> open;
  for (std::size_t i = 0; i < text.size(); i++) {
    const char c = text[i];
    if (c == '(') open.push_back(i);
    else if (c == ')') {
      if (open.empty()) throw std::invalid_argument("genomics::rna::parse_dot_bracket: a ')' without its '(' at " + std::to_string(i));
      s.add(open.back(), i);
      open.pop_back();
    } else if (c != '.' && c != '-' && c != ',' && c != ':' && c != '_' && c != '~') {
      throw std::invalid_argument(std::string("genomics::rna::parse_dot_bracket: the character '") + c + "' at " + std::to_string(i));
    }
  }
  if (!open.empty()) throw std::invalid_argument("genomics::rna::parse_dot_bracket: a '(' without its ')' at " + std::to_string(open.back()));
  return s;
}

// The letters of an RNA: A, C, G, U (codes above 3 as N).
inline std::string rna_letters(const Codes& x) {
  std::string s(x.size(), 'N');
  for (std::size_t i = 0; i < x.size(); i++) s[i] = x[i] < 4 ? "ACGU"[x[i]] : 'N';
  return s;
}

// The base-pair distance: the number of pairs in one structure and not the other.
inline std::size_t base_pair_distance(const Structure& a, const Structure& b) {
  if (a.size() != b.size()) throw std::invalid_argument("genomics::rna::base_pair_distance: structures of different lengths");
  std::size_t d = 0;
  for (std::size_t i = 0; i < a.size(); i++)
    if (a.partner[i] != b.partner[i]) {
      d += a.partner[i] > static_cast<std::int32_t>(i);
      d += b.partner[i] > static_cast<std::int32_t>(i);
    }
  return d;
}

// ---- counting and enumeration ------------------------------------------------------------------------
// S_m(n), exactly: the structures on n positions with every pair allowed and minimum loop m.
inline morphe::arithmetic::Integer structure_count(std::size_t n, std::size_t min_loop) {
  using morphe::arithmetic::Integer;
  std::vector<Integer> c(n + 1, Integer(0));
  c[0] = Integer(1);
  for (std::size_t len = 1; len <= n; len++) {
    Integer t = c[len - 1];
    for (std::size_t k = 0; k + min_loop + 2 <= len; k++) t = t + c[k] * c[len - k - 2];
    c[len] = t;
  }
  return c[n];
}

// Calls visit on every structure of n positions whose pairs satisfy allowed(i, j) and
// j − i − 1 ≥ min_loop. Exponential in n: for checks and examples of a few tens of bases.
inline void enumerate_structures(std::size_t n, std::size_t min_loop, const std::function<bool(std::size_t, std::size_t)>& allowed,
                                 const std::function<void(const Structure&)>& visit) {
  Structure s(n);
  // Decides position i, given that every pair opened before it closes after it or is closed.
  std::function<void(std::size_t, std::vector<std::size_t>&)> go = [&](std::size_t i, std::vector<std::size_t>& open) {
    if (i == n) {
      if (open.empty()) visit(s);
      return;
    }
    if (open.size() > n - i) return;  // too many pairs left open
    // i unpaired
    go(i + 1, open);
    // i closes the innermost open pair
    if (!open.empty()) {
      const std::size_t k = open.back();
      if (i - k - 1 >= min_loop && allowed(k, i)) {
        open.pop_back();
        s.partner[k] = static_cast<std::int32_t>(i);
        s.partner[i] = static_cast<std::int32_t>(k);
        go(i + 1, open);
        s.partner[k] = s.partner[i] = -1;
        open.push_back(k);
      }
    }
    // i opens a pair
    open.push_back(i);
    go(i + 1, open);
    open.pop_back();
  };
  std::vector<std::size_t> open;
  go(0, open);
}

}  // namespace genomics::rna
