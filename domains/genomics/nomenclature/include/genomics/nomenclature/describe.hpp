// A change to a reference described as the HGVS Nomenclature describes it, with the class of its
// equivalent representations and the canonical choices three standards make from it.
//
// - A change replaces the bases [p, p + |D|) of a reference r by I; it makes the haplotype
//   h = r[0, p) · I · r[p + |D|, n). Two changes are equivalent when they make the same haplotype.
//   After trimming the bases D and I share at their ends, either both are non-empty and the change
//   is the only minimal one of its class, or one is empty (a pure deletion or insertion of L bases)
//   and its minimal representations are the placements p of an interval [p_min, p_max]
//   (theory/nomenclature.md §3, Theorem 1, with its closed form by the longest common prefix and
//   suffix of r and h). The interval is found here by sliding: a deletion of r[p, p + L) moves right
//   while r[p] = r[p + L] and left while r[p − 1] = r[p + L − 1]; an insertion of u at p moves right
//   while r[p] = u₀ (u rotating left) and left while r[p − 1] = u_{L−1}. O(|D| + |I| + p_max − p_min).
// - Three canonical forms are three choices from the class (Theorem 2): VCF's left-aligned and
//   parsimonious form, Tan, Abecasis and Kang's Algorithm 1 (variant/normalize.hpp, reused here),
//   is p_min with the base before it; the HGVS Nomenclature's 3′ rule ("the most 3′ position
//   possible of the reference sequence is arbitrarily assigned to have been changed", general.md) is
//   p_max; NCBI's VOCA (Holmes et al. 2020, Bioinformatics 36, 1902–1907, "Correcting overprecision")
//   expands the change to the whole interval, deleting r[p_min, p_max + |D|) — the hull of the class.
// - The HGVS description by the standard's order of preference ("(1) substitution, (2) deletion,
//   (3) inversion, (4) duplication, (5) insertion", general.md), delins for the rest:
//   one base by one is a substitution (DNA/substitution.md); a pure deletion is del at p_max; an
//   equal-length replacement by the reverse complement, of more than one base, is inv
//   (DNA/inversion.md); a pure insertion of u at p_max is a dup of r[p_max − L, p_max) when u is
//   that copy, which happens exactly when p_max − p_min ≥ L ("a copy … inserted directly 3′ of the
//   original copy", DNA/duplication.md); when u is k ≥ 2 copies of its primitive root v and v ends
//   the reference before p_max, it is a repeated sequence over the reference's run of v
//   ("when more than one additional copies are inserted directly 3′ … the format for Repeated
//   sequences", DNA/duplication.md; Mutalyzer 3 writes 4_12CAG[5] and 13_15T[6] so); otherwise ins,
//   between p_max and p_max + 1 (1-based), or, when asked, as a copy of a range of the reference,
//   inverted or not (DNA/insertion.md: g.17_18ins5_16, g.234_235ins123_234inv), the nearest
//   occurrence taken. Everything else is delins of the trimmed window.
//
// Coordinates here are 0-based and half-open; the descriptions written are 1-based, from `first`
// (the number of the reference's first base, 1 for a whole sequence; a region of a chromosome passes
// its own). A change at the very ends that HGVS has no flank for (an insertion before the first or
// after the last base) is written as a delins of the end base. Letters are read in capitals.

#pragma once

#include <genomics/nomenclature/hgvs.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/variant/normalize.hpp>

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::nomenclature {

using sequence::Code;
using sequence::Codes;

// The bases [position, position + deleted.size()) of the reference replaced by `inserted`.
struct Change {
  std::size_t position = 0;
  std::string deleted, inserted;
  friend bool operator==(const Change&, const Change&) = default;
};

// The class of a change's equivalent representations.
struct Placements {
  Change minimal;               // trimmed of the bases its two sides share at their ends
  bool identity = false;        // the change makes the reference
  bool indel = false;           // a pure deletion or insertion: its placements form an interval
  bool deletion = false;        // of the indels, a deletion
  std::size_t length = 0;       // L, the bases deleted or inserted
  std::size_t first = 0, last = 0;            // p_min and p_max (the minimal change's position)
  std::string inserted_first, inserted_last;  // an insertion's sequence at p_min and at p_max (its rotations)
  std::size_t hull_begin = 0, hull_end = 0;   // [p_min, p_max + |D|): the bases any equivalent minimal change touches
  Change at_first() const { return deletion || !indel ? Change{first, minimal.deleted, minimal.inserted} : Change{first, "", inserted_first}; }
  Change at_last() const {
    if (!indel) return minimal;
    return deletion ? Change{last, minimal.deleted, ""} : Change{last, "", inserted_last};
  }
};

namespace detail {

inline std::string upper(std::string_view s) {
  std::string out(s);
  for (char& c : out) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return out;
}
inline char letter(const Codes& r, std::size_t i) { return sequence::dna_letter(r[i]); }

inline std::string slice(const Codes& r, std::size_t b, std::size_t e) {
  std::string s;
  s.reserve(e - b);
  for (std::size_t i = b; i < e; i++) s += letter(r, i);
  return s;
}

// The shortest v with u = v^k.
inline std::string primitive_root(std::string_view u) {
  const std::size_t n = u.size();
  for (std::size_t d = 1; d <= n; d++) {
    if (n % d) continue;
    bool ok = true;
    for (std::size_t i = d; i < n && ok; i++) ok = u[i] == u[i - d];
    if (ok) return std::string(u.substr(0, d));
  }
  return std::string(u);
}

}  // namespace detail

// Whether the change's deleted bases are the reference's.
inline bool deleted_agrees(const Codes& reference, const Change& c) {
  if (c.position + c.deleted.size() > reference.size()) return false;
  for (std::size_t i = 0; i < c.deleted.size(); i++)
    if (sequence::dna_code(c.deleted[i]) != reference[c.position + i]) return false;
  return true;
}

// The haplotype a change makes.
inline Codes apply_change(const Codes& reference, const Change& c) {
  if (!deleted_agrees(reference, c)) throw std::invalid_argument("genomics::nomenclature::apply_change: the deleted bases are not the reference's");
  Codes h(reference.begin(), reference.begin() + static_cast<std::ptrdiff_t>(c.position));
  for (char x : c.inserted) h.push_back(sequence::dna_code(x));
  h.insert(h.end(), reference.begin() + static_cast<std::ptrdiff_t>(c.position + c.deleted.size()), reference.end());
  return h;
}

inline Placements placements(const Codes& reference, Change change) {
  change.deleted = detail::upper(change.deleted);
  change.inserted = detail::upper(change.inserted);
  if (!deleted_agrees(reference, change)) throw std::invalid_argument("genomics::nomenclature::placements: the deleted bases are not the reference's at the position given");
  // Trim the common suffix, then the common prefix (VOCA's shrinking; Tan et al.'s order).
  std::string& d = change.deleted;
  std::string& u = change.inserted;
  while (!d.empty() && !u.empty() && d.back() == u.back()) { d.pop_back(); u.pop_back(); }
  std::size_t k = 0;
  while (k < d.size() && k < u.size() && d[k] == u[k]) k++;
  d.erase(0, k);
  u.erase(0, k);
  change.position += k;
  Placements out;
  out.minimal = change;
  out.first = out.last = change.position;
  const std::size_t n = reference.size();
  if (d.empty() && u.empty()) {
    out.identity = true;
    out.hull_begin = out.hull_end = change.position;
    return out;
  }
  if (!d.empty() && !u.empty()) {
    out.hull_begin = change.position;
    out.hull_end = change.position + d.size();
    return out;
  }
  out.indel = true;
  if (!d.empty()) {
    out.deletion = true;
    const std::size_t L = d.size();
    out.length = L;
    std::size_t lo = change.position, hi = change.position;
    while (hi + L < n && reference[hi] == reference[hi + L]) hi++;
    while (lo > 0 && reference[lo - 1] == reference[lo + L - 1]) lo--;
    out.first = lo;
    out.last = hi;
    out.minimal = Change{lo, detail::slice(reference, lo, lo + L), ""};
    out.hull_begin = lo;
    out.hull_end = hi + L;
    return out;
  }
  const std::size_t L = u.size();
  out.length = L;
  std::string right = u, left = u;
  std::size_t hi = change.position, lo = change.position;
  while (hi < n && sequence::dna_code(right.front()) == reference[hi]) {
    std::rotate(right.begin(), right.begin() + 1, right.end());
    hi++;
  }
  while (lo > 0 && sequence::dna_code(left.back()) == reference[lo - 1]) {
    std::rotate(left.rbegin(), left.rbegin() + 1, left.rend());
    lo--;
  }
  out.first = lo;
  out.last = hi;
  out.inserted_first = left;
  out.inserted_last = right;
  out.minimal = Change{lo, "", left};
  out.hull_begin = lo;
  out.hull_end = hi;
  return out;
}

// ---- the three canonical forms ------------------------------------------------------------------------
// VCF's: Tan et al.'s normal form, POS 0-based of REF's first base, with the padding base.
inline variant::Normalized vcf_form(const Codes& reference, const Change& change) {
  const Placements pl = placements(reference, change);
  Change c = pl.at_first();
  std::string ref = c.deleted, alt = c.inserted;
  std::size_t pos = c.position;
  if (ref.empty() || alt.empty()) {
    if (pos > 0) {
      const char b = detail::letter(reference, pos - 1);
      ref.insert(ref.begin(), b);
      alt.insert(alt.begin(), b);
      pos--;
    } else if (ref.size() < reference.size()) {
      const char b = detail::letter(reference, ref.size());
      ref += b;
      alt += b;
    }
  }
  if (ref.empty()) ref = alt;  // an empty reference: nothing to normalise
  return variant::normalize(pos, ref, alt, reference);
}

// SPDI's contextual allele (VOCA): the class's hull deleted and the haplotype over it inserted.
struct Spdi {
  std::string sequence;     // the reference's accession and version
  std::size_t position = 0; // 0-based interbase: the bases to advance before deleting
  std::string deleted, inserted;
  std::string text() const { return sequence + ":" + std::to_string(position) + ":" + deleted + ":" + inserted; }
  friend bool operator==(const Spdi&, const Spdi&) = default;
};

inline Spdi spdi_form(const Codes& reference, const Change& change, std::string sequence_id = {}) {
  const Placements pl = placements(reference, change);
  Spdi s;
  s.sequence = std::move(sequence_id);
  if (pl.identity) {
    s.position = change.position;
    s.deleted = s.inserted = detail::upper(change.deleted);
    return s;
  }
  s.position = pl.hull_begin;
  s.deleted = detail::slice(reference, pl.hull_begin, pl.hull_end);
  if (!pl.indel) {
    s.inserted = pl.minimal.inserted;
  } else if (pl.deletion) {
    s.inserted = detail::slice(reference, pl.first + pl.length, pl.hull_end);
  } else {
    s.inserted = pl.inserted_first + s.deleted;
  }
  return s;
}

// ---- the HGVS description ---------------------------------------------------------------------------------
struct DescribeOptions {
  std::int64_t first = 1;           // the number of the reference's first base
  bool repeats = true;              // write insertions of two or more tandem copies as repeated sequences
  bool unit_multiple_of_three = false;  // only units whose length is a multiple of three (a coding sequence, DNA/repeated.md)
  std::size_t copy_length = 0;      // write an inserted sequence of at least this length as a copy of a reference range, when it is one (0: never)
  std::size_t copy_window = 100000; // how far from the insertion to look for the copy
  bool open_start = false, open_end = false;  // the reference is a window: bases lie beyond its start or its end
};

struct Described {
  Variant variant;        // in g. numbering from options.first
  Placements placements;
  Change change;          // the placement described (p_max for an indel)
  std::string text() const { return format_variant(variant); }
};

namespace detail {

inline Position number_position(std::int64_t n) {
  Position p;
  p.number = n;
  return p;
}
inline Location<Position> bases(std::int64_t a, std::int64_t b) {  // 1-based inclusive
  Location<Position> l;
  l.start.first = l.start.last = number_position(a);
  l.end.first = l.end.last = number_position(b);
  l.range = a != b;
  if (!l.range) l.end = l.start;
  return l;
}

// The nearest occurrence of u in r within the window, forward or reverse-complemented: [begin, end) and whether inverted.
inline bool find_copy(const Codes& r, const std::string& u, std::size_t at, std::size_t window, std::size_t& begin, bool& inverted) {
  const std::size_t L = u.size();
  if (L == 0 || L > r.size()) return false;
  const Codes fw = sequence::encode_dna(u), rc = sequence::reverse_complement(fw);
  const std::size_t lo = at > window ? at - window : 0, hi = std::min(r.size() - L, at + window);
  std::size_t best = SIZE_MAX, best_distance = SIZE_MAX;
  bool best_inv = false;
  for (std::size_t i = lo; i <= hi; i++) {
    for (int k = 0; k < 2; k++) {
      const Codes& w = k ? rc : fw;
      if (!std::equal(w.begin(), w.end(), r.begin() + static_cast<std::ptrdiff_t>(i))) continue;
      const std::size_t e = i + L;
      const std::size_t dist = e <= at ? at - e : i >= at ? i - at : 0;
      if (dist < best_distance || (dist == best_distance && k == 0 && best_inv)) { best = i; best_distance = dist; best_inv = k; }
    }
  }
  if (best == SIZE_MAX) return false;
  begin = best;
  inverted = best_inv;
  return true;
}

}  // namespace detail

inline Described describe(const Codes& reference, const Change& change, const DescribeOptions& o = {}) {
  Described out;
  out.placements = placements(reference, change);
  const Placements& pl = out.placements;
  out.change = pl.at_last();
  const Change& c = out.change;
  Variant& v = out.variant;
  const std::int64_t base = o.first;  // number of reference[0]
  const auto num = [&](std::size_t i) { return base + static_cast<std::int64_t>(i); };  // 0-based → numbered
  const std::size_t n = reference.size();
  if (pl.identity) {
    v.kind = Kind::identity;
    const std::size_t L = std::max<std::size_t>(1, change.deleted.size());
    v.location = detail::bases(num(change.position), num(change.position + L - 1));
    return out;
  }
  if (!pl.indel) {
    const std::size_t D = c.deleted.size();
    v.location = detail::bases(num(c.position), num(c.position + D - 1));
    if (D == 1 && c.inserted.size() == 1) {
      v.kind = Kind::substitution;
      v.was = c.deleted;
      v.now = c.inserted;
    } else if (D == c.inserted.size() && D >= 2 && sequence::reverse_complement(c.deleted) == c.inserted) {
      v.kind = Kind::inversion;
    } else {
      v.kind = Kind::delins;
      v.inserted.push_back(Piece{Piece::Kind::sequence, c.inserted, {}, {}, false, {}, Coordinate::g});
    }
    return out;
  }
  const std::size_t L = pl.length;
  if (pl.deletion) {
    v.kind = Kind::deletion;
    v.location = detail::bases(num(c.position), num(c.position + L - 1));
    return out;
  }
  const std::string& u = c.inserted;
  const std::size_t p = c.position;  // interbase: between p − 1 and p
  // A repeated sequence: k ≥ 2 copies of the primitive root v, v ending the reference before p.
  if (o.repeats) {
    const std::string root = detail::primitive_root(u);
    const std::size_t m = root.size(), k = L / m;
    const bool unit_ok = !o.unit_multiple_of_three || m % 3 == 0;
    if (k >= 2 && unit_ok && p >= m && detail::slice(reference, p - m, p) == root) {
      std::size_t s = p, j = 0;
      while (s >= m && detail::slice(reference, s - m, s) == root) { s -= m; j++; }
      v.kind = Kind::repeat;
      v.location = detail::bases(num(s), num(p - 1));
      v.units.push_back(RepeatUnit{root, Copies{Copies::Kind::exact, static_cast<std::int64_t>(j + k), static_cast<std::int64_t>(j + k), false}});
      return out;
    }
  }
  if (pl.last - pl.first >= L) {
    v.kind = Kind::duplication;
    v.location = detail::bases(num(p - L), num(p - 1));
    return out;
  }
  if ((p == 0 && !o.open_start) || (p == n && !o.open_end)) {
    // No flank on one side: a delins of the end base.
    v.kind = Kind::delins;
    const std::size_t at = p == 0 ? 0 : n - 1;
    v.location = detail::bases(num(at), num(at));
    const std::string end(1, detail::letter(reference, at));
    v.inserted.push_back(Piece{Piece::Kind::sequence, p == 0 ? u + end : end + u, {}, {}, false, {}, Coordinate::g});
    out.change = Change{at, end, p == 0 ? u + end : end + u};
    return out;
  }
  v.kind = Kind::insertion;
  v.location = detail::bases(num(p - 1), num(p));
  v.location.range = true;
  std::size_t from = 0;
  bool inverted = false;
  if (o.copy_length && L >= o.copy_length && detail::find_copy(reference, u, p, o.copy_window, from, inverted)) {
    Piece piece;
    piece.kind = Piece::Kind::location;
    piece.location = detail::bases(num(from), num(from + L - 1));
    piece.inverted = inverted;
    v.inserted.push_back(piece);
  } else {
    v.inserted.push_back(Piece{Piece::Kind::sequence, u, {}, {}, false, {}, Coordinate::g});
  }
  return out;
}

// ---- the meaning of a description ----------------------------------------------------------------------------------
// The haplotype an exact g.-numbered variant makes of the reference (theory/nomenclature.md §3: a
// description denotes a haplotype): its location replaced by what the edit says. Nothing when the
// variant is not exact (uncertain or unknown positions, counts of unknown bases, other references)
// or does not fit the reference (a substitution whose `was` is not the reference's base).
inline std::optional<Change> change_of(const Codes& reference, const Variant& v, std::int64_t first = 1) {
  if (v.protein) return std::nullopt;
  const auto& l = v.location;
  if (l.start.uncertain || l.end.uncertain || l.start.first.kind != Position::Kind::number || l.end.first.kind != Position::Kind::number ||
      l.start.first.offset_sign || l.end.first.offset_sign || l.start.first.utr || l.end.first.utr)
    return std::nullopt;
  const std::int64_t a = l.start.first.number - first, b = l.end.first.number - first;  // 0-based, inclusive
  if (a < 0 || b < a || b >= static_cast<std::int64_t>(reference.size())) return std::nullopt;
  const auto ua = static_cast<std::size_t>(a), ub = static_cast<std::size_t>(b);
  const std::string named = detail::slice(reference, ua, ub + 1);
  const auto pieces = [&]() -> std::optional<std::string> {
    std::string s;
    for (const auto& p : v.inserted) {
      if (p.kind == Piece::Kind::remote) return std::nullopt;
      if (p.kind == Piece::Kind::sequence) {
        if (p.count) return std::nullopt;
        s += detail::upper(p.sequence);
        continue;
      }
      const auto c = change_of(reference, Variant{Kind::identity, false, p.location}, first);
      if (!c) return std::nullopt;
      s += p.inverted ? sequence::reverse_complement(c->deleted) : c->deleted;
    }
    return s;
  };
  switch (v.kind) {
    case Kind::substitution:
      if (l.range || detail::upper(v.was) != named) return std::nullopt;
      return Change{ua, named, detail::upper(v.now)};
    case Kind::identity: return Change{ua, named, named};
    case Kind::deletion: return Change{ua, named, ""};
    case Kind::duplication: return Change{ub + 1, "", named};
    case Kind::inversion: return Change{ua, named, sequence::reverse_complement(named)};
    case Kind::insertion: {
      if (!l.range || ub != ua + 1) return std::nullopt;
      const auto s = pieces();
      if (!s) return std::nullopt;
      return Change{ub, "", *s};
    }
    case Kind::delins: {
      const auto s = pieces();
      if (!s) return std::nullopt;
      return Change{ua, named, *s};
    }
    case Kind::repeat: {
      if (v.units.empty()) return std::nullopt;
      std::string s;
      for (const auto& u : v.units) {
        if (u.copies.kind != Copies::Kind::exact) return std::nullopt;
        for (std::int64_t k = 0; k < u.copies.low; k++) s += detail::upper(u.sequence);
      }
      // The located range is the repeat's reference extent (from its first position alone, the run of the unit there).
      if (l.range) return Change{ua, named, s};
      const std::string unit = detail::upper(v.units.front().sequence);
      std::size_t e = ua;
      while (e + unit.size() <= reference.size() && detail::slice(reference, e, e + unit.size()) == unit) e += unit.size();
      return Change{ua, detail::slice(reference, ua, e), s};
    }
    default: return std::nullopt;
  }
}

// The description with its reference: NC_000023.11:g.33344591del.
inline std::string describe_text(const Codes& reference, const Change& change, const std::string& accession, Coordinate type = Coordinate::g,
                                 const DescribeOptions& o = {}) {
  const Described d = describe(reference, change, o);
  return (accession.empty() ? std::string() : accession + ":") + static_cast<char>(type) + "." + format_variant(d.variant);
}

}  // namespace genomics::nomenclature
