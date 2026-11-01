// Breakends and adjacencies: the language in which every structural variant is written
// (theory/structural-variants.md, §1).
//
// - A breakend b = (c, p, o) is a side of a base of the reference: o = + is the junction after
//   base p of contig c, the reference to its left kept; o = − the junction before base p, the
//   reference to its right kept. A read mapped forward points to a + breakend downstream of it,
//   one mapped reverse to a − breakend upstream of it. These are VCF's breakend orientations
//   (The Variant Call Format Specification v4.2, §5.4: "t[p[" and "t]p]" join after t, "]p]t"
//   and "[p[t" before it).
// - An adjacency {b₁, b₂} is a junction of a rearranged genome: reading the derivative chromosome
//   across it one leaves the reference at b₁ and enters it at b₂, with possibly some novel bases
//   between. Read on the other strand the same junction leaves at b₂ and enters at b₁: an
//   adjacency is an unordered pair, written here with a ≤ b in the order (contig, position,
//   orientation). Its signature names the rearrangement that makes it alone: on one contig,
//   (+, −) a deletion (a < b), (−, +) a tandem duplication, (+, +) or (−, −) one of the two
//   junctions of an inversion; on two contigs a translocation.
// - Reading across a junction. With the derivative's bases indexed by k, k < 0 before the junction
//   and k ≥ 0 after it, the side left at a breakend (p, +) is ref[p + 1 + k], at (p, −) the
//   complement of ref[p − 1 − k]; the side entered at (q, −) is ref[q + k], at (q, +) the
//   complement of ref[q − k]. Both formulas extend to every k, which is what a shift of the
//   junction does: moved m bases later in the derivative, (p, +) becomes (p + m, +), (p, −)
//   becomes (p − m, −), (q, −) becomes (q + m, −), (q, +) becomes (q − m, +); the derivative is
//   unchanged while the leaving side's extension equals the entering side's bases. The bases over
//   which it can slide are the junction's microhomology, and an adjacency without novel bases is
//   written canonically at the placement that puts its first breakend leftmost on the reference,
//   as VCF left-aligns an indel (Tan et al. 2015; the same lattice argument, theory note §1.3).
//
// Complexity: homology is O(h) for h bases of homology, bounded by `limit`.

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

namespace genomics::structural {

using sequence::Code;
using sequence::Codes;

struct Breakend {
  std::uint32_t contig = 0;
  std::size_t position = 0;
  bool plus = true;  // + : the reference to the left of the position kept (junction after it)
  friend bool operator==(const Breakend&, const Breakend&) = default;
  friend bool operator<(const Breakend& x, const Breakend& y) { return std::tie(x.contig, x.position, x.plus) < std::tie(y.contig, y.position, y.plus); }
};

// A junction: left at `a`, entered at `b` (as read in one direction), with novel bases between.
struct Adjacency {
  Breakend a, b;
  std::string inserted;  // letters, read in the direction that leaves at a
};

enum class Signature : std::uint8_t { deletion, duplication, inversion_plus, inversion_minus, translocation };

inline const char* signature_name(Signature s) {
  switch (s) {
    case Signature::deletion: return "deletion";
    case Signature::duplication: return "duplication";
    case Signature::inversion_plus: return "inversion(+,+)";
    case Signature::inversion_minus: return "inversion(−,−)";
    case Signature::translocation: return "translocation";
  }
  return "?";
}

// The signature of an adjacency written with a ≤ b.
inline Signature signature(const Adjacency& j) {
  if (j.a.contig != j.b.contig) return Signature::translocation;
  if (j.a.plus && !j.b.plus) return Signature::deletion;
  if (!j.a.plus && j.b.plus) return Signature::duplication;
  return j.a.plus ? Signature::inversion_plus : Signature::inversion_minus;
}

// The complement of a letter of the alphabet ACGTN.
inline char complement_letter(char c) { return sequence::dna_letter(sequence::complement(sequence::dna_code(c))); }

// The same junction read in the other direction: it leaves at b and enters at a, the novel bases
// reverse-complemented.
inline Adjacency reversed(const Adjacency& j) {
  Adjacency r{j.b, j.a, {}};
  for (auto it = j.inserted.rbegin(); it != j.inserted.rend(); ++it) r.inserted += complement_letter(*it);
  return r;
}

// Written with a ≤ b.
inline Adjacency ordered(const Adjacency& j) { return j.b < j.a ? reversed(j) : j; }

namespace detail {

// The base at derivative offset k of the side left at e (k < 0 its own bases, k ≥ 0 its
// extension), or N outside the contig.
inline Code leaving(const std::vector<Codes>& contigs, const Breakend& e, std::int64_t k) {
  const Codes& x = contigs[e.contig];
  const std::int64_t i = e.plus ? static_cast<std::int64_t>(e.position) + 1 + k : static_cast<std::int64_t>(e.position) - 1 - k;
  if (i < 0 || i >= static_cast<std::int64_t>(x.size())) return sequence::N;
  const Code c = x[static_cast<std::size_t>(i)];
  return e.plus ? c : sequence::complement(c);
}
// The base at derivative offset k of the side entered at e (k ≥ 0 its own bases, k < 0 its extension).
inline Code entering(const std::vector<Codes>& contigs, const Breakend& e, std::int64_t k) {
  const Codes& x = contigs[e.contig];
  const std::int64_t i = e.plus ? static_cast<std::int64_t>(e.position) - k : static_cast<std::int64_t>(e.position) + k;
  if (i < 0 || i >= static_cast<std::int64_t>(x.size())) return sequence::N;
  const Code c = x[static_cast<std::size_t>(i)];
  return e.plus ? sequence::complement(c) : c;
}

inline Breakend moved(Breakend e, std::int64_t m, bool leaving_side) {
  // m > 0: the junction m bases later in the derivative.
  const bool up = leaving_side ? e.plus : !e.plus;
  const std::int64_t p = static_cast<std::int64_t>(e.position) + (up ? m : -m);
  e.position = static_cast<std::size_t>(std::max<std::int64_t>(0, p));
  return e;
}

}  // namespace detail

// The junction moved m bases later (m < 0: earlier) in the derivative read leaving at a.
inline Adjacency shifted(const Adjacency& j, std::int64_t m) {
  Adjacency r = j;
  r.a = detail::moved(j.a, m, true);
  r.b = detail::moved(j.b, m, false);
  return r;
}

// How far the junction can slide without changing the derivative: `left` bases earlier and
// `right` later (0 when it carries novel bases), each at most `limit`.
struct Homology {
  std::size_t left = 0, right = 0;
  std::size_t length() const { return left + right; }
};

inline Homology homology(const std::vector<Codes>& contigs, const Adjacency& j, std::size_t limit = 1000) {
  Homology h;
  if (!j.inserted.empty()) return h;
  for (std::int64_t k = 0; h.right < limit; k++) {
    const Code x = detail::leaving(contigs, j.a, k), y = detail::entering(contigs, j.b, k);
    if (x > 3 || x != y) break;
    h.right++;
  }
  for (std::int64_t k = -1; h.left < limit; k--) {
    const Code x = detail::leaving(contigs, j.a, k), y = detail::entering(contigs, j.b, k);
    if (x > 3 || x != y) break;
    h.left++;
  }
  return h;
}

// The canonical form: written with a ≤ b, then slid as far as the derivative allows in the direction
// that lowers a's position (earlier in the derivative when a is +, later when it is −), so that a
// deletion and a tandem duplication are left-aligned on the reference as VCF left-aligns them; an
// adjacency with novel bases is only ordered.
inline Adjacency canonical(const std::vector<Codes>& contigs, const Adjacency& j, std::size_t limit = 1000) {
  Adjacency x = ordered(j);
  if (!x.inserted.empty()) return x;
  const Homology h = homology(contigs, x, limit);
  const std::int64_t m = x.a.plus ? -static_cast<std::int64_t>(std::min(h.left, x.a.position)) : static_cast<std::int64_t>(std::min(h.right, x.a.position));
  return m ? ordered(shifted(x, m)) : x;
}

// The derivative's bases around the junction, `flank` each side, novel bases between (letters).
inline std::string junction_sequence(const std::vector<Codes>& contigs, const Adjacency& j, std::size_t flank) {
  std::string s;
  const auto f = static_cast<std::int64_t>(flank);
  for (std::int64_t k = -f; k < 0; k++) s += sequence::dna_letter(detail::leaving(contigs, j.a, k));
  s += j.inserted;
  for (std::int64_t k = 0; k < f; k++) s += sequence::dna_letter(detail::entering(contigs, j.b, k));
  return s;
}

}  // namespace genomics::structural
