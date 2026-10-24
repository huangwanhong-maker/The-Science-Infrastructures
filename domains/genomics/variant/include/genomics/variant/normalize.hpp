// The normalisation of variants: Tan, Abecasis and Kang's left-alignment and trimming of the
// alleles of a variant (2015, "Unified representation of genetic variants", Bioinformatics 31,
// 2202–2204, Algorithm 1), and the left-alignment of the gaps of an alignment, its counterpart
// on reads.
//
// - A variant is a position p (0-based, of REF's first base) and alleles: REF, the reference's
//   bases there, and the ALTs. Its effect is the haplotype it makes of the reference: REF's
//   bases replaced by an ALT's. Two variants are equivalent when they make the same haplotype;
//   an indel in a repeat has many equivalent representations (a CA deleted from CACACA can be
//   written at four positions, with longer or shorter context). The normalised representation is
//   the one Algorithm 1 reaches:
//
//     repeat
//       if every allele ends with the same base, drop that base from each;
//       if some allele is empty, prepend the reference base before p to each, and p ← p − 1;
//     until nothing changes;
//     while every allele starts with the same base and has at least two, drop that base.
//
//   It is *parsimonious* (no shorter alleles represent the variant) and *left-aligned* (no
//   equivalent representation has a smaller position), and it is unique: the theory note
//   (mapping-and-variants.md, §7) proves that it is the canonical form of the class of
//   equivalent representations. At p = 0 nothing precedes the first base; an allele left empty
//   there takes the base after it, as VCF 4.2 pads an event at position 1 (§1.6.1, REF). Letters
//   are read in upper case. O(ℓ + s) for alleles of total length ℓ shifted by s bases.
// - Gaps of an alignment are left-aligned without changing its score: a deletion of target
//   bases t[j, j + L) after a pair (q[i − 1], t[j − 1]) moves one column left when
//   t[j − 1] = t[j + L − 1], the pair then joining q[i − 1] to t[j + L − 1], which scores as
//   t[j − 1] did; an insertion of q[i, i + L) moves left when q[i − 1] = q[i + L − 1]. Gaps
//   stop at the alignment's start and at another gap (and merge with one of the same kind).
//   The path's score never falls, and is unchanged when no gaps merge. This is what makes reads
//   that carry one indel in a repeat report it at one place, whatever the aligner chose (GATK's
//   LeftAlignIndels does the same; DePristo et al. 2011, Nature Genetics 43, 491–498). O(n s)
//   for a path of n columns and shifts of s.

#pragma once

#include <genomics/align/types.hpp>
#include <genomics/sequence/alphabet.hpp>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::variant {

using sequence::Code;
using sequence::Codes;

struct Normalized {
  std::size_t position = 0;
  std::string ref;
  std::vector<std::string> alts;
  bool changed = false;      // the representation was not already normal
  bool ref_matches = true;   // REF agreed with the reference where it was given
};

namespace detail {

inline std::string upper(std::string s) {
  for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

}  // namespace detail

// Whether REF's letters are the reference's bases at p.
inline bool ref_agrees(std::size_t p, std::string_view ref, const Codes& reference) {
  if (p + ref.size() > reference.size()) return false;
  for (std::size_t i = 0; i < ref.size(); i++)
    if (sequence::dna_code(ref[i]) != reference[p + i]) return false;
  return true;
}

// Tan, Abecasis and Kang's Algorithm 1 on a variant of any number of alleles.
inline Normalized normalize(std::size_t position, std::string ref, std::vector<std::string> alts, const Codes& reference) {
  if (alts.empty()) throw std::invalid_argument("genomics::variant::normalize: a variant has at least one ALT");
  Normalized out;
  out.ref_matches = ref_agrees(position, detail::upper(ref), reference);
  std::vector<std::string> a;
  a.reserve(alts.size() + 1);
  a.push_back(detail::upper(std::move(ref)));
  for (auto& x : alts) a.push_back(detail::upper(std::move(x)));
  std::size_t p = position;
  if (p > reference.size()) throw std::invalid_argument("genomics::variant::normalize: the position lies beyond the reference");
  const auto all_end_same = [&] {
    for (const auto& x : a)
      if (x.empty() || x.back() != a[0].back()) return false;
    return true;
  };
  const auto some_empty = [&] {
    for (const auto& x : a)
      if (x.empty()) return true;
    return false;
  };
  bool padded_right = false;
  for (bool change = true; change;) {
    change = false;
    if (all_end_same()) {
      for (auto& x : a) x.pop_back();
      change = true;
    }
    if (some_empty()) {
      if (p == 0) {
        // Nothing before the first base: pad with the base after REF, once, and stop.
        const std::size_t after = a[0].size();
        if (after < reference.size()) {
          const char b = sequence::dna_letter(reference[after]);
          for (auto& x : a) x.push_back(b);
          padded_right = true;
        }
        break;
      }
      const char b = sequence::dna_letter(reference[p - 1]);
      for (auto& x : a) x.insert(x.begin(), b);
      p--;
      change = true;
    }
  }
  if (!padded_right) {
    const auto all_start_same = [&] {
      for (const auto& x : a)
        if (x.size() < 2 || x.front() != a[0].front()) return false;
      return true;
    };
    while (all_start_same()) {
      for (auto& x : a) x.erase(x.begin());
      p++;
    }
  }
  out.position = p;
  out.ref = std::move(a[0]);
  out.alts.assign(std::make_move_iterator(a.begin() + 1), std::make_move_iterator(a.end()));
  return out;
}

inline Normalized normalize(std::size_t position, std::string ref, std::string alt, const Codes& reference) {
  std::vector<std::string> alts{std::move(alt)};
  Normalized n = normalize(position, std::move(ref), std::move(alts), reference);
  return n;
}

// The haplotype a variant makes of the reference: REF's bases at p replaced by ALT's.
inline Codes apply(const Codes& reference, std::size_t position, std::string_view ref, std::string_view alt) {
  if (position + ref.size() > reference.size()) throw std::invalid_argument("genomics::variant::apply: REF runs past the reference");
  Codes h(reference.begin(), reference.begin() + static_cast<std::ptrdiff_t>(position));
  for (char c : alt) h.push_back(sequence::dna_code(c));
  h.insert(h.end(), reference.begin() + static_cast<std::ptrdiff_t>(position + ref.size()), reference.end());
  return h;
}

// ---- the gaps of an alignment ---------------------------------------------------------------------
// The path of an alignment (without clipping) with every gap moved as far left as it goes
// without lowering the score.
inline align::Cigar left_align_gaps(const align::Alignment& a, const Codes& query, const Codes& target) {
  // One operation a column: M (a pair), I or D.
  std::string col;
  for (const auto& o : a.cigar) {
    const char op = (o.op == '=' || o.op == 'X') ? 'M' : o.op;
    if (op != 'M' && op != 'I' && op != 'D') throw std::invalid_argument("genomics::variant::left_align_gaps: a path of M, =, X, I and D");
    col.append(o.length, op);
  }
  // q and t positions at the start of each column.
  std::vector<std::size_t> qi(col.size() + 1), tj(col.size() + 1);
  const auto positions = [&] {
    std::size_t i = a.query_begin, j = a.target_begin;
    for (std::size_t c = 0; c < col.size(); c++) {
      qi[c] = i;
      tj[c] = j;
      if (col[c] != 'D') i++;
      if (col[c] != 'I') j++;
    }
    qi[col.size()] = i;
    tj[col.size()] = j;
  };
  positions();
  std::size_t c = 0;
  while (c < col.size()) {
    if (col[c] == 'M') { c++; continue; }
    const char g = col[c];
    std::size_t e = c;
    while (e < col.size() && col[e] == g) e++;
    const std::size_t L = e - c;
    // Shift [c, e) left while the column before is a pair and the condition holds.
    std::size_t start = c;
    while (start > 0 && col[start - 1] == 'M') {
      bool can;
      if (g == 'D') can = target[tj[start] - 1] == target[tj[start] + L - 1];
      else can = query[qi[start] - 1] == query[qi[start] + L - 1];
      if (!can) break;
      // The pair moves after the gap.
      col[start - 1] = g;
      col[start - 1 + L] = 'M';
      // Positions of the moved columns: the gap now starts one column earlier.
      start--;
      // Recompute the two positions that changed (column start + L holds the pair now).
      std::size_t i = qi[start], j = tj[start];
      for (std::size_t k = start; k <= start + L; k++) {
        qi[k] = i;
        tj[k] = j;
        if (col[k] != 'D') i++;
        if (col[k] != 'I') j++;
      }
    }
    c = std::max(start + L, c + 1);
    // Skip past the gap's (possibly merged) run.
    while (c < col.size() && col[c] == g) c++;
  }
  align::Cigar out;
  for (const char op : col) align::push(out, op);
  return out;
}

}  // namespace genomics::variant
