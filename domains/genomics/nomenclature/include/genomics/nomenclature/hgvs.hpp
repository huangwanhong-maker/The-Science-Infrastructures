// HGVS descriptions as typed values: read from text by the grammar (grammar.hpp), written back in
// the standard's canonical form, checked for the conditions of meaning that a context-free grammar
// does not express, and read out in plain words.
//
// - The value (Hgvs) follows the grammar's structure (theory/nomenclature.md §2): alternatives
//   `^` of descriptions; a description is a reference (accession, version, selector), a coordinate
//   type and a list of variants; the list is alternatives `^` of groups joined by `(;)` (unknown
//   phase); a group is alleles in trans `[…];[…]`, or one chain; an allele is chains separated by
//   `;` (cis), `(;)` or `,` (transcripts of one allele), possibly predicted as a whole; a chain is
//   items joined by `::`; an item is a variant, or a predicted choice of variants `(a^b)`.
// - A variant (Variant) is a location and an edit. Nucleotide positions keep everything written
//   (c.-14, c.*32, c.88+1, c.4072-?, pter); a location's ends may be uncertain, (a_b) or (a); residue
//   positions are a residue and a number (Lys23, Ter110). Edits are the standard's kinds; inserted
//   material is a list of pieces (a sequence with a count, a range of the same reference, inverted or
//   not, a range of another reference).
// - The canonical form (format) is the one the standard recommends where it gives a choice
//   (general.md, DNA/deletion.md, protein/substitution.md): three-letter amino acids with Ter, the
//   deleted or duplicated sequence not repeated after del or dup (it is kept when it carries
//   information: a count, N[15], or an uncertain location), no spaces. Style asks for one-letter
//   codes or '*' instead.
// - Conditions of meaning (check): a substitution names one position and one base; a range names
//   two different positions, 5′ to 3′; an insertion is between two adjacent positions; an inversion
//   covers more than one; c., n. and r. positions name their offsets on an exon's boundary only in
//   the sense numbering.md gives; three-letter codes are preferred to one-letter ones, Xaa to X
//   (background/standards.md). Each finding is a Diagnostic with its severity and where it is.
//
// Readers never throw on content: parse_hgvs returns whether the text is in the language, the
// named tree, the furthest viable prefix and what may follow it (morphe's Earley diagnostics), the
// value and the findings. Time O(n²) in the length of the text on this unambiguous grammar.

#pragma once

#include <genomics/nomenclature/grammar.hpp>

#include <morphe/formal/bnf.hpp>

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace genomics::nomenclature {

using morphe::formal::ParseNode;

enum class Coordinate : char { g = 'g', m = 'm', o = 'o', c = 'c', n = 'n', r = 'r', p = 'p' };

inline bool nucleotide(Coordinate c) { return c != Coordinate::p; }
inline bool genomic(Coordinate c) { return c == Coordinate::g || c == Coordinate::m || c == Coordinate::o; }

struct Reference {
  std::string accession, version;                    // NM_004006, 2
  bool has_selector = false;
  std::string selector_accession, selector_version;  // the transcript or protein selected: NG_012232.1(NM_004006.2)
  std::string id() const { return version.empty() ? accession : accession + "." + version; }
  std::string selector() const { return selector_version.empty() ? selector_accession : selector_accession + "." + selector_version; }
  std::string text() const { return id() + (has_selector ? "(" + selector() + ")" : ""); }
  friend bool operator==(const Reference&, const Reference&) = default;
};

// A nucleotide position as written (numbering.md): g.1234, c.-14 (5′ UTR), c.*32 (3′ UTR), c.88+1 and
// c.89-2 (intronic offsets), c.4072-? (an unknown offset), ?, pter, qter, cen.
struct Position {
  enum class Kind : std::uint8_t { number, unknown, pter, qter, cen };
  Kind kind = Kind::number;
  char utr = 0;               // '-', '*' or 0
  std::int64_t number = 0;    // as written, ≥ 0
  char offset_sign = 0;       // '+', '-' or 0 (no offset)
  std::int64_t offset = 0;    // its magnitude
  bool offset_unknown = false;
  std::int64_t signed_offset() const { return offset_sign == '-' ? -offset : offset; }
  friend bool operator==(const Position&, const Position&) = default;
};

// A residue position: a residue (one-letter code, '*' for the stop) and its number, Lys23, Ter110.
struct Residue {
  char residue = 'X';
  std::int64_t number = 0;
  friend bool operator==(const Residue&, const Residue&) = default;
};

// A point of a location: exact (a), somewhere in (a_b), or a parenthesised (a).
template <class P>
struct Point {
  P first, last;
  bool uncertain = false;  // written in parentheses
  bool span = false;       // (first_last); otherwise one position
  friend bool operator==(const Point&, const Point&) = default;
};

template <class P>
struct Location {
  Point<P> start, end;
  bool range = false;  // start_end; otherwise the start alone
  friend bool operator==(const Location&, const Location&) = default;
};

// A number of copies: [23], [?], [(50_60)], [(60_?)].
struct Copies {
  enum class Kind : std::uint8_t { exact, unknown, range };
  Kind kind = Kind::exact;
  std::int64_t low = 0, high = 0;
  bool high_unknown = false;
  friend bool operator==(const Copies&, const Copies&) = default;
};

// A piece of inserted material.
struct Piece {
  enum class Kind : std::uint8_t { sequence, location, remote };
  Kind kind = Kind::sequence;
  std::string sequence;                // in the coordinate type's letters
  std::optional<Copies> count;         // A[26], N[15]
  Location<Position> location;         // a range of this or of the remote reference
  bool inverted = false;
  Reference reference;                 // remote
  Coordinate coordinate = Coordinate::g;
  friend bool operator==(const Piece&, const Piece&) = default;
};

struct RepeatUnit {
  std::string sequence;
  Copies copies;
  friend bool operator==(const RepeatUnit&, const RepeatUnit&) = default;
};

enum class Kind : std::uint8_t {
  substitution, deletion, duplication, insertion, delins, inversion, identity, repeat, supernumerary, methylation,
  frameshift, extension, n_extension, unknown
};

inline const char* kind_name(Kind k) {
  static constexpr const char* names[] = {"substitution", "deletion", "duplication", "insertion", "delins", "inversion", "identity",
                                          "repeat", "supernumerary", "methylation", "frameshift", "extension", "n_extension", "unknown"};
  return names[static_cast<std::size_t>(k)];
}

// One variant: a location and an edit.
struct Variant {
  Kind kind = Kind::substitution;
  bool protein = false;
  Location<Position> location;  // nucleotides
  Location<Residue> residues;   // proteins
  std::uint8_t mixture = 0;     // 1: mosaic "=/", 2: chimeric "=//" (some cells unchanged at the location)
  // substitution: nucleotide `was` > `now`; protein `now` holds the alternatives, one-letter (Gly56Ala^Ser^Cys → "ASC")
  std::string was, now;
  // deletion and duplication: the sequence written after them, if any (N[15] for an unsequenced one)
  std::string stated;
  std::optional<Copies> stated_count;
  // insertion and delins (nucleotides)
  std::vector<Piece> inserted;
  bool bracketed = false;
  // insertion and delins (proteins): a sequence (one-letter) with a count, a stop at a length, or a range
  std::string residues_inserted;
  std::optional<Copies> inserted_count;
  std::int64_t stop_length = 0;  // insTer12
  bool inserted_range = false;
  Location<Residue> inserted_from;
  // repeat: units with their copies, or the located segment's copies; the other allele's counts
  std::vector<RepeatUnit> units;
  std::optional<Copies> copies;
  std::vector<Copies> others;
  // methylation: gom, lom, met=
  std::string methylation;
  // frameshift (fs_new: the new residue, or its alternatives) and extensions: the new stop
  std::string fs_new;
  bool has_stop = false, stop_unknown = false;
  std::int64_t stop_at = 0;
  std::int64_t upstream = 0;  // Met1ext-5: 5
  friend bool operator==(const Variant&, const Variant&) = default;
};

struct Item {
  std::vector<Variant> choices;  // one, or the alternatives of (a^b)
  bool predicted = false;
  std::uint8_t mixed = 0;        // "=/v": the whole sequence unchanged in some cells, v in others
  friend bool operator==(const Item&, const Item&) = default;
};

struct Chain {
  std::vector<Item> items;  // joined by "::"
  friend bool operator==(const Chain&, const Chain&) = default;
};

struct Allele {
  std::vector<Chain> chains;
  std::string separators;  // one per gap: ';' in cis, '?' unknown phase "(;)", ',' transcripts of one allele
  bool predicted = false;  // [(a;b)]
  std::string whole;       // [=], [?], [0]
  bool supernumerary = false;
  friend bool operator==(const Allele&, const Allele&) = default;
};

struct Group {
  bool bracketed = false;
  std::vector<Allele> alleles;  // in trans
  Chain chain;                  // when not bracketed
  friend bool operator==(const Group&, const Group&) = default;
};

struct Alternative {
  std::vector<Group> groups;  // joined by "(;)"
  friend bool operator==(const Alternative&, const Alternative&) = default;
};

struct VariantList {
  Coordinate coordinate = Coordinate::g;
  std::vector<Alternative> alternatives;  // joined by "^"
  std::string whole;                      // c.=, p.?, p.0, r.spl, ...
  friend bool operator==(const VariantList&, const VariantList&) = default;
};

struct AdjoinedPart {
  Reference reference;
  Point<Position> first, last;
  friend bool operator==(const AdjoinedPart&, const AdjoinedPart&) = default;
};

struct Description {
  bool has_reference = false;
  Reference reference;
  VariantList variants;
  bool adjoined = false;
  AdjoinedPart five, three;
  std::string linker;
  friend bool operator==(const Description&, const Description&) = default;
};

struct Hgvs {
  std::vector<Description> alternatives;            // joined by "^"
  bool genome_alleles = false;
  std::vector<std::vector<Description>> alleles;    // [d;d];[d]
  friend bool operator==(const Hgvs&, const Hgvs&) = default;
};

// ---- amino acids ------------------------------------------------------------------------------------
// IUPAC-IUB's three-letter codes and HGVS's Ter (background/standards.md), with Pyl, which that
// table leaves out.
inline const char* three_letter(char one) {
  switch (one) {
    case 'A': return "Ala"; case 'R': return "Arg"; case 'N': return "Asn"; case 'D': return "Asp"; case 'C': return "Cys";
    case 'Q': return "Gln"; case 'E': return "Glu"; case 'G': return "Gly"; case 'H': return "His"; case 'I': return "Ile";
    case 'L': return "Leu"; case 'K': return "Lys"; case 'M': return "Met"; case 'F': return "Phe"; case 'P': return "Pro";
    case 'S': return "Ser"; case 'T': return "Thr"; case 'W': return "Trp"; case 'Y': return "Tyr"; case 'V': return "Val";
    case 'U': return "Sec"; case 'O': return "Pyl"; case 'B': return "Asx"; case 'Z': return "Glx"; case 'J': return "Xle";
    case 'X': return "Xaa"; case '*': return "Ter";
    default: return "Xaa";
  }
}
inline char one_letter(std::string_view three) {
  static constexpr std::string_view letters = "ARNDCQEGHILKMFPSTWYVUOBZJX*";
  for (char c : letters)
    if (three == three_letter(c)) return c;
  return 0;
}

struct Style {
  bool one_letter = false;   // p.W24C rather than p.Trp24Cys
  bool star = false;         // '*' rather than Ter in three-letter descriptions
  bool keep_stated = false;  // keep a deleted or duplicated sequence written after del or dup
};

// ---- reading the named tree --------------------------------------------------------------------------
namespace detail {

// A node's name without its family prefix: "c_deletion" → "deletion".
inline std::string_view bare(const ParseNode& n) {
  const std::string_view s = n.name;
  return s.size() > 2 && s[1] == '_' && (s[0] == 'g' || s[0] == 'c' || s[0] == 'n' || s[0] == 'r' || s[0] == 'p') ? s.substr(2) : s;
}
inline const ParseNode* child(const ParseNode& n, std::string_view name) {
  for (const auto& c : n.children)
    if (bare(c) == name) return &c;
  return nullptr;
}
inline std::vector<const ParseNode*> children(const ParseNode& n, std::string_view name) {
  std::vector<const ParseNode*> out;
  for (const auto& c : n.children)
    if (bare(c) == name) out.push_back(&c);
  return out;
}
inline std::int64_t integer(std::string_view s) {
  std::int64_t v = 0;
  for (char c : s) v = v * 10 + (c - '0');
  return v;
}

inline Reference reference_of(const ParseNode& n) {
  Reference r;
  const auto accessions = children(n, "accession");
  const auto read = [](const ParseNode& a, std::string& acc, std::string& ver) {
    acc = child(a, "identifier") ? child(a, "identifier")->text : a.text;
    ver = child(a, "version") ? child(a, "version")->text : "";
  };
  if (!accessions.empty()) read(*accessions[0], r.accession, r.version);
  if (const ParseNode* s = child(n, "selector")) {
    r.has_selector = true;
    if (const ParseNode* a = child(*s, "accession")) read(*a, r.selector_accession, r.selector_version);
  }
  return r;
}

inline Position position_of(const ParseNode& n) {
  Position p;
  const std::string& t = n.text;
  if (t == "?") { p.kind = Position::Kind::unknown; return p; }
  if (t == "pter") { p.kind = Position::Kind::pter; return p; }
  if (t == "qter") { p.kind = Position::Kind::qter; return p; }
  if (t == "cen") { p.kind = Position::Kind::cen; return p; }
  if (const ParseNode* u = child(n, "utr")) p.utr = u->text[0];
  if (const ParseNode* k = child(n, "number")) p.number = integer(k->text);
  if (const ParseNode* o = child(n, "offset")) {
    p.offset_sign = o->text[0];
    const ParseNode* amount = child(*o, "offset_amount");
    if (amount && amount->text == "?") p.offset_unknown = true;
    else if (amount) p.offset = integer(amount->text);
  }
  return p;
}

inline Residue residue_of(const ParseNode& n) {  // a p_pos
  Residue r;
  if (const ParseNode* x = child(n, "residue")) r.residue = x->text.size() == 1 ? x->text[0] : one_letter(x->text);
  if (const ParseNode* k = child(n, "number")) r.number = integer(k->text);
  return r;
}

template <class P, class F>
Point<P> point_of(const ParseNode& upos, F&& pos) {  // an *_upos node
  Point<P> pt;
  if (const ParseNode* u = child(upos, "uncertain")) {
    pt.uncertain = true;
    const auto ps = children(*u, "pos");
    if (!ps.empty()) pt.first = pt.last = pos(*ps[0]);
    if (ps.size() > 1) { pt.last = pos(*ps[1]); pt.span = true; }
  } else if (const ParseNode* p = child(upos, "pos")) {
    pt.first = pt.last = pos(*p);
  }
  return pt;
}

template <class P, class F>
Location<P> location_of(const ParseNode& loc, F&& pos) {
  Location<P> l;
  const auto ups = children(loc, "upos");
  if (!ups.empty()) l.start = l.end = point_of<P>(*ups[0], pos);
  if (ups.size() > 1) { l.end = point_of<P>(*ups[1], pos); l.range = true; }
  return l;
}

inline Copies copies_of(const ParseNode& count) {  // a count node: "[" copies "]"
  Copies c;
  const ParseNode* k = child(count, "copies");
  if (!k) return c;
  if (k->text == "?") { c.kind = Copies::Kind::unknown; return c; }
  const auto nums = children(*k, "number");
  if (k->text.front() == '(') {
    c.kind = Copies::Kind::range;
    if (!nums.empty()) c.low = integer(nums[0]->text);
    const ParseNode* up = child(*k, "copies_upper");
    if (up && up->text == "?") c.high_unknown = true;
    else if (up) c.high = integer(up->text);
    return c;
  }
  c.low = c.high = integer(k->text);
  return c;
}

inline std::string protein_sequence(const ParseNode& n) {  // a p_sequence: one-letter codes
  std::string out;
  const std::string& t = n.text;
  if (!t.empty() && t[0] >= 'A' && t[0] <= 'Z' && (t.size() < 2 || !(t[1] >= 'a' && t[1] <= 'z'))) {
    return t;  // already one-letter
  }
  for (std::size_t i = 0; i < t.size();) {
    if (t[i] == '*') { out += '*'; i++; continue; }
    out += one_letter(std::string_view(t).substr(i, 3));
    i += 3;
  }
  return out;
}

inline char residue_letter(const ParseNode& residue) { return residue.text.size() == 1 ? residue.text[0] : one_letter(residue.text); }

inline Coordinate coordinate_of_remote(const ParseNode& remote) {
  for (const auto& c : remote.children) {
    if (c.name == "g_type") return static_cast<Coordinate>(c.text[0]);
    if (c.name == "c_loc") return Coordinate::c;
    if (c.name == "n_loc") return Coordinate::n;
    if (c.name == "r_loc") return Coordinate::r;
  }
  return Coordinate::g;
}

inline Variant nucleotide_variant(const ParseNode& edit) {
  Variant v;
  const ParseNode* loc = child(edit, "loc");
  if (loc) v.location = location_of<Position>(*loc, position_of);
  if (const ParseNode* m = child(edit, "mixture")) v.mixture = m->text == "/" ? 1 : 2;
  const ParseNode* change = child(edit, "change");
  if (!change || change->children.empty()) return v;
  const ParseNode& c = change->children.front();
  const std::string_view k = bare(c);
  if (k == "substitution") {
    v.kind = Kind::substitution;
    if (const ParseNode* w = child(c, "was")) v.was = w->text;
    if (const ParseNode* w = child(c, "now")) v.now = w->text;
  } else if (k == "deletion" || k == "duplication") {
    v.kind = k == "deletion" ? Kind::deletion : Kind::duplication;
    if (const ParseNode* s = child(c, "stated")) {
      if (const ParseNode* q = child(*s, "sequence")) v.stated = q->text;
      if (const ParseNode* n = child(*s, "count")) v.stated_count = copies_of(*n);
    }
  } else if (k == "insertion" || k == "delins") {
    v.kind = k == "insertion" ? Kind::insertion : Kind::delins;
    if (const ParseNode* ins = child(c, "inserted")) {
      v.bracketed = !ins->text.empty() && ins->text.front() == '[';
      std::vector<const ParseNode*> pieces;
      for (const auto& x : ins->children) {
        if (bare(x) == "piece") pieces.push_back(&x);
        else if (bare(x) == "bracketed" && !x.children.empty()) pieces.push_back(&x);
      }
      for (const ParseNode* q : pieces) {
        Piece piece;
        const ParseNode* p = bare(*q) == "bracketed" && child(*q, "piece") ? child(*q, "piece") : q;
        if (const ParseNode* r = child(*q, "remote")) {
          piece.kind = Piece::Kind::remote;
          if (const ParseNode* ref = child(*r, "reference")) piece.reference = reference_of(*ref);
          piece.coordinate = coordinate_of_remote(*r);
          for (const auto& rc : r->children)
            if (rc.name == "g_loc" || rc.name == "c_loc" || rc.name == "n_loc" || rc.name == "r_loc") piece.location = location_of<Position>(rc, position_of);
          piece.inverted = child(*r, "inverted") != nullptr;
        } else if (const ParseNode* l = child(*p, "loc")) {
          piece.kind = Piece::Kind::location;
          piece.location = location_of<Position>(*l, position_of);
          piece.inverted = child(*p, "inverted") != nullptr;
        } else {
          piece.kind = Piece::Kind::sequence;
          if (const ParseNode* q = child(*p, "sequence")) piece.sequence = q->text;
          if (const ParseNode* n = child(*p, "count")) piece.count = copies_of(*n);
        }
        v.inserted.push_back(std::move(piece));
      }
    }
  } else if (k == "inversion") {
    v.kind = Kind::inversion;
  } else if (k == "identity") {
    v.kind = Kind::identity;
    if (const ParseNode* w = child(c, "was")) v.was = w->text;
  } else if (k == "repeat") {
    v.kind = Kind::repeat;
    for (const ParseNode* u : children(c, "unit")) {
      RepeatUnit unit;
      if (const ParseNode* q = child(*u, "sequence")) unit.sequence = q->text;
      if (const ParseNode* n = child(*u, "count")) unit.copies = copies_of(*n);
      v.units.push_back(std::move(unit));
    }
    if (const ParseNode* n = child(c, "count")) v.copies = copies_of(*n);
    for (const ParseNode* o : children(c, "other_allele"))
      if (const ParseNode* n = child(*o, "count")) v.others.push_back(copies_of(*n));
  } else if (k == "supernumerary") {
    v.kind = Kind::supernumerary;
  } else if (k == "methylation") {
    v.kind = Kind::methylation;
    if (const ParseNode* s = child(c, "methylation_state")) v.methylation = s->text;
  }
  return v;
}

inline Variant protein_variant(const ParseNode& edit) {
  Variant v;
  v.protein = true;
  const auto rpos = [](const ParseNode& n) { return residue_of(n); };
  if (const ParseNode* loc = child(edit, "loc")) v.residues = location_of<Residue>(*loc, rpos);
  if (const ParseNode* m = child(edit, "mixture")) v.mixture = m->text == "/" ? 1 : 2;
  const ParseNode* change = child(edit, "change");
  if (!change || change->children.empty()) return v;
  const ParseNode& c = change->children.front();
  const std::string_view k = bare(c);
  const auto stop_of = [&](const ParseNode& node) {
    if (const ParseNode* s = child(node, "new_stop")) {
      v.has_stop = true;
      const ParseNode* at = child(*s, "stop_at");
      if (at && at->text == "?") v.stop_unknown = true;
      else if (at) v.stop_at = integer(at->text);
    }
  };
  if (k == "substitution") {
    v.kind = Kind::substitution;
    for (const ParseNode* r : children(c, "residue")) v.now += residue_letter(*r);
  } else if (k == "deletion") {
    v.kind = Kind::deletion;
  } else if (k == "duplication") {
    v.kind = Kind::duplication;
  } else if (k == "insertion" || k == "delins") {
    v.kind = k == "insertion" ? Kind::insertion : Kind::delins;
    if (const ParseNode* ins = child(c, "inserted")) {
      if (const ParseNode* s = child(*ins, "sequence")) v.residues_inserted = protein_sequence(*s);
      if (const ParseNode* n = child(*ins, "count")) v.inserted_count = copies_of(*n);
      if (const ParseNode* sl = child(*ins, "stop_length"))
        if (const ParseNode* n = child(*sl, "number")) v.stop_length = integer(n->text);
      if (const ParseNode* r = child(*ins, "inserted_range")) {
        v.inserted_range = true;
        const auto ps = children(*r, "pos");
        if (ps.size() == 2) {
          v.inserted_from.range = true;
          v.inserted_from.start.first = v.inserted_from.start.last = residue_of(*ps[0]);
          v.inserted_from.end.first = v.inserted_from.end.last = residue_of(*ps[1]);
        }
      }
    }
  } else if (k == "frameshift") {
    v.kind = Kind::frameshift;
    if (const ParseNode* nw = child(c, "fs_new"))
      for (const ParseNode* r : children(*nw, "residue")) v.fs_new += residue_letter(*r);
    stop_of(c);
  } else if (k == "extension") {
    v.kind = Kind::extension;
    if (const ParseNode* r = child(c, "residue")) v.now = std::string(1, residue_letter(*r));
    stop_of(c);
  } else if (k == "n_extension") {
    v.kind = Kind::n_extension;
    if (const ParseNode* n = child(c, "number")) v.upstream = integer(n->text);
  } else if (k == "identity") {
    v.kind = Kind::identity;
  } else if (k == "unknown") {
    v.kind = Kind::unknown;
  } else if (k == "repeat") {
    v.kind = Kind::repeat;
    if (const ParseNode* n = child(c, "count")) v.copies = copies_of(*n);
    for (const ParseNode* o : children(c, "other_allele"))
      if (const ParseNode* n = child(*o, "count")) v.others.push_back(copies_of(*n));
  }
  return v;
}

inline Variant variant_of(const ParseNode& edit) { return edit.name[0] == 'p' ? protein_variant(edit) : nucleotide_variant(edit); }

inline Item item_of(const ParseNode& item) {
  Item it;
  if (const ParseNode* e = child(item, "edit")) {
    it.choices.push_back(variant_of(*e));
  } else if (const ParseNode* p = child(item, "predicted")) {
    it.predicted = true;
    for (const ParseNode* e : children(*p, "edit")) it.choices.push_back(variant_of(*e));
    if (const ParseNode* m = child(*p, "mixed")) {
      if (const ParseNode* x = child(*m, "mixture")) it.mixed = x->text == "/" ? 1 : 2;
      if (const ParseNode* e = child(*m, "edit")) it.choices.push_back(variant_of(*e));
    }
  } else if (const ParseNode* m = child(item, "mixed")) {
    if (const ParseNode* x = child(*m, "mixture")) it.mixed = x->text == "/" ? 1 : 2;
    if (const ParseNode* e = child(*m, "edit")) it.choices.push_back(variant_of(*e));
  }
  return it;
}

inline Chain chain_of(const ParseNode& chain) {
  Chain c;
  for (const ParseNode* i : children(chain, "item")) c.items.push_back(item_of(*i));
  return c;
}

inline char separator_char(const ParseNode& s) { return s.text == ";" ? ';' : s.text == "," ? ',' : '?'; }

inline Allele allele_of(const ParseNode& allele) {
  Allele a;
  a.supernumerary = child(allele, "supernumerary") != nullptr;
  const ParseNode* body = child(allele, "body");
  if (!body) return a;
  const ParseNode* list = child(*body, "predicted_list");
  const ParseNode& holder = list ? *list : *body;
  a.predicted = list != nullptr;
  for (const auto& c : holder.children) {
    if (bare(c) == "chain") a.chains.push_back(chain_of(c));
    else if (c.name == "separator") a.separators += separator_char(c);
    else if (c.name == "whole" || c.name == "r_whole") a.whole = c.text;
  }
  return a;
}

inline VariantList list_of(const ParseNode& variants) {  // an *_variants node
  VariantList l;
  if (variants.name == "g_variants") l.coordinate = static_cast<Coordinate>(variants.text[0]);
  else l.coordinate = static_cast<Coordinate>(variants.name[0]);
  const ParseNode* list = child(variants, "list");
  if (!list) return l;
  for (const auto& c : list->children) {
    if (c.name == "whole" || c.name == "r_whole") { l.whole = c.text; continue; }
    if (bare(c) != "alternative") continue;
    Alternative alt;
    for (const ParseNode* g : children(c, "group")) {
      Group grp;
      const auto alleles = children(*g, "allele");
      if (!alleles.empty()) {
        grp.bracketed = true;
        for (const ParseNode* a : alleles) grp.alleles.push_back(allele_of(*a));
      } else if (const ParseNode* ch = child(*g, "chain")) {
        grp.chain = chain_of(*ch);
      }
      alt.groups.push_back(std::move(grp));
    }
    l.alternatives.push_back(std::move(alt));
  }
  return l;
}

inline AdjoinedPart adjoined_part_of(const ParseNode& part) {
  AdjoinedPart a;
  if (const ParseNode* r = child(part, "reference")) a.reference = reference_of(*r);
  const auto ups = children(part, "upos");
  if (ups.size() == 2) {
    a.first = point_of<Position>(*ups[0], position_of);
    a.last = point_of<Position>(*ups[1], position_of);
  }
  return a;
}

inline Description description_of(const ParseNode& d) {
  Description out;
  if (const ParseNode* adj = child(d, "adjoined")) {
    out.adjoined = true;
    out.variants.coordinate = Coordinate::r;
    const auto parts = children(*adj, "adjoined_part");
    if (parts.size() == 2) {
      out.five = adjoined_part_of(*parts[0]);
      out.three = adjoined_part_of(*parts[1]);
    }
    if (const ParseNode* l = child(*adj, "linker")) out.linker = l->text;
    return out;
  }
  if (const ParseNode* r = child(d, "reference")) {
    out.has_reference = true;
    out.reference = reference_of(*r);
  }
  for (const auto& c : d.children)
    if (bare(c) == "variants") out.variants = list_of(c);
  return out;
}

}  // namespace detail

// The value of a named tree from the start rule.
inline Hgvs hgvs_of(const ParseNode& tree) {
  Hgvs h;
  if (const ParseNode* ga = detail::child(tree, "genome_alleles")) {
    h.genome_alleles = true;
    for (const ParseNode* a : detail::children(*ga, "genome_allele")) {
      std::vector<Description> ds;
      for (const ParseNode* d : detail::children(*a, "description")) ds.push_back(detail::description_of(*d));
      h.alleles.push_back(std::move(ds));
    }
    return h;
  }
  for (const ParseNode* d : detail::children(tree, "description")) h.alternatives.push_back(detail::description_of(*d));
  return h;
}

// ---- writing ---------------------------------------------------------------------------------------------
inline std::string format_position(const Position& p) {
  switch (p.kind) {
    case Position::Kind::unknown: return "?";
    case Position::Kind::pter: return "pter";
    case Position::Kind::qter: return "qter";
    case Position::Kind::cen: return "cen";
    case Position::Kind::number: break;
  }
  std::string s;
  if (p.utr) s += p.utr;
  s += std::to_string(p.number);
  if (p.offset_sign) {
    s += p.offset_sign;
    s += p.offset_unknown ? "?" : std::to_string(p.offset);
  }
  return s;
}

inline std::string format_residue(char r, const Style& st) {
  if (st.one_letter) return std::string(1, r);
  if (r == '*') return st.star ? "*" : "Ter";
  return three_letter(r);
}
inline std::string format_residue_position(const Residue& r, const Style& st) { return format_residue(r.residue, st) + std::to_string(r.number); }

template <class P, class F>
std::string format_point(const Point<P>& p, F&& f) {
  if (!p.uncertain) return f(p.first);
  return "(" + f(p.first) + (p.span ? "_" + f(p.last) : std::string()) + ")";
}
template <class P, class F>
std::string format_location(const Location<P>& l, F&& f) {
  return format_point(l.start, f) + (l.range ? "_" + format_point(l.end, f) : std::string());
}

inline std::string format_copies(const Copies& c) {
  switch (c.kind) {
    case Copies::Kind::unknown: return "[?]";
    case Copies::Kind::range: return "[(" + std::to_string(c.low) + "_" + (c.high_unknown ? std::string("?") : std::to_string(c.high)) + ")]";
    case Copies::Kind::exact: break;
  }
  return "[" + std::to_string(c.low) + "]";
}

inline std::string format_protein_sequence(std::string_view s, const Style& st) {
  std::string out;
  for (char c : s) out += format_residue(c, st);
  return out;
}

inline std::string format_variant(const Variant& v, const Style& st = {}) {
  std::string s;
  if (v.protein) {
    s = format_location(v.residues, [&](const Residue& r) { return format_residue_position(r, st); });
  } else {
    s = format_location(v.location, format_position);
  }
  if (v.mixture) s += v.mixture == 1 ? "=/" : "=//";
  const auto stop = [&] {
    if (!v.has_stop) return std::string();
    return format_residue('*', st) + (v.stop_unknown ? std::string("?") : std::to_string(v.stop_at));
  };
  switch (v.kind) {
    case Kind::substitution:
      if (v.protein) {
        for (std::size_t i = 0; i < v.now.size(); i++) s += (i ? "^" : "") + format_residue(v.now[i], st);
      } else {
        s += v.was + ">" + v.now;
      }
      break;
    case Kind::deletion:
    case Kind::duplication: {
      s += v.kind == Kind::deletion ? "del" : "dup";
      const bool informative = v.stated_count.has_value() || v.location.start.uncertain || v.location.end.uncertain;
      if (!v.protein && (st.keep_stated || informative)) {
        s += v.stated;
        if (v.stated_count) s += format_copies(*v.stated_count);
      }
      break;
    }
    case Kind::insertion:
    case Kind::delins: {
      s += v.kind == Kind::insertion ? "ins" : "delins";
      if (v.protein) {
        if (v.stop_length) s += format_residue('*', st) + std::to_string(v.stop_length);
        else if (v.inserted_range)
          s += format_location(v.inserted_from, [&](const Residue& r) { return format_residue_position(r, st); });
        else {
          s += format_protein_sequence(v.residues_inserted, st);
          if (v.inserted_count) s += format_copies(*v.inserted_count);
        }
        break;
      }
      const bool brackets = v.bracketed || v.inserted.size() > 1 || (v.inserted.size() == 1 && v.inserted[0].kind == Piece::Kind::remote);
      if (brackets) s += "[";
      for (std::size_t i = 0; i < v.inserted.size(); i++) {
        if (i) s += ";";
        const Piece& p = v.inserted[i];
        if (p.kind == Piece::Kind::sequence) {
          s += p.sequence;
          if (p.count) s += format_copies(*p.count);
        } else {
          if (p.kind == Piece::Kind::remote) s += p.reference.text() + ":" + static_cast<char>(p.coordinate) + ".";
          s += format_location(p.location, format_position);
          if (p.inverted) s += "inv";
        }
      }
      if (brackets) s += "]";
      break;
    }
    case Kind::inversion: s += "inv"; break;
    case Kind::identity: s += (v.protein ? std::string() : v.was) + "="; break;
    case Kind::repeat:
      for (const auto& u : v.units) s += u.sequence + format_copies(u.copies);
      if (v.copies) s += format_copies(*v.copies);
      for (const auto& o : v.others) s += ";" + format_copies(o);
      break;
    case Kind::supernumerary: s += "sup"; break;
    case Kind::methylation: s += "|" + v.methylation; break;
    case Kind::frameshift:
      if (v.fs_new.size() == 1) s += format_residue(v.fs_new[0], st);
      else if (v.fs_new.size() > 1) {
        s += "(";
        for (std::size_t i = 0; i < v.fs_new.size(); i++) s += (i ? "^" : "") + format_residue(v.fs_new[i], st);
        s += ")";
      }
      s += "fs" + stop();
      break;
    case Kind::extension: s += format_residue(v.now.empty() ? 'X' : v.now[0], st) + "ext" + stop(); break;
    case Kind::n_extension: s += "ext-" + std::to_string(v.upstream); break;
    case Kind::unknown: s += "?"; break;
  }
  return s;
}

inline std::string format_item(const Item& it, const Style& st = {}) {
  std::string s;
  if (it.mixed) s += it.mixed == 1 ? "=/" : "=//";
  if (it.mixed && it.predicted) return "(" + s + format_variant(it.choices.front(), st) + ")";
  for (std::size_t i = 0; i < it.choices.size(); i++) s += (i ? "^" : "") + format_variant(it.choices[i], st);
  return it.predicted ? "(" + s + ")" : s;
}
inline std::string format_chain(const Chain& c, const Style& st = {}) {
  std::string s;
  for (std::size_t i = 0; i < c.items.size(); i++) s += (i ? "::" : "") + format_item(c.items[i], st);
  return s;
}
inline std::string format_allele(const Allele& a, const Style& st = {}) {
  std::string s;
  if (!a.whole.empty()) s = a.whole;
  for (std::size_t i = 0; i < a.chains.size(); i++) {
    if (i) {
      const char sep = i - 1 < a.separators.size() ? a.separators[i - 1] : ';';
      s += sep == ';' ? ";" : sep == ',' ? "," : "(;)";
    }
    s += format_chain(a.chains[i], st);
  }
  if (a.predicted) s = "(" + s + ")";
  return "[" + s + "]" + (a.supernumerary ? "sup" : "");
}
inline std::string format_list(const VariantList& l, const Style& st = {}) {
  if (!l.whole.empty()) return l.whole;
  std::string s;
  for (std::size_t i = 0; i < l.alternatives.size(); i++) {
    if (i) s += "^";
    const auto& alt = l.alternatives[i];
    for (std::size_t j = 0; j < alt.groups.size(); j++) {
      if (j) s += "(;)";
      const Group& g = alt.groups[j];
      if (g.bracketed) {
        for (std::size_t k = 0; k < g.alleles.size(); k++) s += (k ? ";" : "") + format_allele(g.alleles[k], st);
      } else {
        s += format_chain(g.chain, st);
      }
    }
  }
  return s;
}
inline std::string format_description(const Description& d, const Style& st = {}) {
  if (d.adjoined) {
    const auto part = [](const AdjoinedPart& p) {
      return p.reference.text() + ":r." + format_point(p.first, format_position) + "_" + format_point(p.last, format_position);
    };
    return part(d.five) + "::" + (d.linker.empty() ? std::string() : d.linker + "::") + part(d.three);
  }
  return (d.has_reference ? d.reference.text() + ":" : std::string()) + static_cast<char>(d.variants.coordinate) + "." + format_list(d.variants, st);
}
inline std::string format(const Hgvs& h, const Style& st = {}) {
  std::string s;
  if (h.genome_alleles) {
    for (std::size_t i = 0; i < h.alleles.size(); i++) {
      if (i) s += ";";
      s += "[";
      for (std::size_t j = 0; j < h.alleles[i].size(); j++) s += (j ? ";" : "") + format_description(h.alleles[i][j], st);
      s += "]";
    }
    return s;
  }
  for (std::size_t i = 0; i < h.alternatives.size(); i++) s += (i ? "^" : "") + format_description(h.alternatives[i], st);
  return s;
}

// ---- conditions of meaning --------------------------------------------------------------------------------
struct Diagnostic {
  enum class Severity : std::uint8_t { error, warning, note };
  Severity severity = Severity::error;
  std::string message;
  std::string where;  // the variant, as written in canonical form
};

namespace detail {

// The order of two positions of one coordinate type, when it is decided by what is written:
// −1, 0, +1, or 2 when it is not (unknown positions, offsets on different anchors are compared by
// their anchors first).
inline int compare(const Position& a, const Position& b) {
  if (a.kind != Position::Kind::number || b.kind != Position::Kind::number) {
    if (a.kind == Position::Kind::pter && b.kind != Position::Kind::pter) return -1;
    if (b.kind == Position::Kind::qter && a.kind != Position::Kind::qter) return -1;
    if (a.kind == b.kind && a.kind != Position::Kind::unknown) return 0;
    return 2;
  }
  const auto anchor = [](const Position& p) -> std::pair<int, std::int64_t> {
    return {p.utr == '-' ? 0 : p.utr == '*' ? 2 : 1, p.utr == '-' ? -p.number : p.number};
  };
  const auto x = anchor(a), y = anchor(b);
  if (x != y) return x < y ? -1 : 1;
  if (a.offset_unknown || b.offset_unknown) return a.signed_offset() == b.signed_offset() && a.offset_unknown == b.offset_unknown ? 0 : 2;
  const std::int64_t u = a.signed_offset(), w = b.signed_offset();
  return u == w ? 0 : u < w ? -1 : 1;
}

// Whether b is the position right after a, as far as the text decides it; nothing when it takes the
// transcript to decide (across an intron's middle, c.88+5_89-5; the CDS's last base and c.*1).
inline std::optional<bool> next_position(const Position& a, const Position& b) {
  if (a.kind != Position::Kind::number || b.kind != Position::Kind::number || a.offset_unknown || b.offset_unknown) return std::nullopt;
  if (a.utr == b.utr && a.number == b.number) return b.signed_offset() == a.signed_offset() + 1;
  if (a.offset_sign == '+' && b.offset_sign == '-') return std::nullopt;
  if (a.offset_sign || b.offset_sign) return false;
  if (!a.utr && b.utr == '*' && b.number == 1) return std::nullopt;
  const auto linear = [](const Position& p) { return p.utr == '-' ? -p.number : p.utr == '*' ? (std::int64_t{1} << 40) + p.number : p.number; };
  const std::int64_t x = linear(a), y = linear(b);
  if (a.utr == '-' && !b.utr) return a.number == 1 && b.number == 1;  // c.-1 and c.1: there is no c.0
  return y == x + 1;
}

inline void check_variant(const Variant& v, Coordinate type, std::vector<Diagnostic>& out) {
  const std::string where = format_variant(v);
  const auto add = [&](Diagnostic::Severity s, std::string m) { out.push_back({s, std::move(m), where}); };
  using S = Diagnostic::Severity;
  if (v.protein) {
    const auto& l = v.residues;
    if (l.range && !l.start.uncertain && !l.end.uncertain) {
      if (l.start.first.number == l.end.first.number) add(S::error, "a range names two different residues (Cys76_Glu79, not Cys76_Cys76)");
      else if (l.start.first.number > l.end.first.number) add(S::error, "a range is written from the N- to the C-terminus (Cys76_Glu79, not Glu79_Cys76)");
    }
    if (v.kind == Kind::insertion && l.range && !l.start.uncertain && !l.end.uncertain && l.end.first.number != l.start.first.number + 1)
      add(S::error, "an insertion lies between two adjacent residues (Lys23_Leu24)");
    if (v.kind == Kind::insertion && !l.range && !l.start.uncertain)
      add(S::error, "an insertion names its two flanking residues (p.Lys23_Leu24insAsp, not p.Lys23insAsp)");
    if (v.kind == Kind::substitution && l.range && !l.start.uncertain)
      add(S::error, "a substitution changes one residue; two or more are a delins (p.Trp24_Val25delinsCysArg)");
    if (v.kind == Kind::frameshift && v.fs_new.size() == 1 && v.fs_new[0] == '*')
      add(S::error, "a stop at the first changed residue is a nonsense substitution (p.Tyr4Ter, not p.Tyr4TerfsTer1)");
    if (v.kind == Kind::frameshift && v.has_stop && !v.stop_unknown && v.stop_at < 2)
      add(S::error, "the shortest frameshift is fsTer2: a stop at the first changed residue is a substitution");
    for (char c : v.now + v.fs_new + v.residues_inserted)
      if (c == 'X') add(S::warning, "an unknown residue is written Xaa (X once meant a stop; background/standards.md)");
    if (v.kind == Kind::substitution && !l.range && !l.start.uncertain && l.start.first.residue == 'M' && l.start.first.number == 1 && v.now != "M")
      add(S::error, "a change of the initiating Met1 is not a substitution: write p.0, p.0? or p.(Met1?) (protein/substitution.md)");
    if (v.kind == Kind::frameshift && v.fs_new.size() == 1 && !l.start.uncertain && v.fs_new[0] == l.start.first.residue)
      add(S::error, "a frameshift is written from the first residue it changes (p.Gln151Thrfs*9, not p.His150Hisfs*10)");
    if (v.kind == Kind::delins && !v.residues_inserted.empty() && !l.start.uncertain && !l.end.uncertain &&
        (v.residues_inserted.front() == l.start.first.residue || v.residues_inserted.back() == l.end.first.residue))
      add(S::error, "the delins is not minimal: a residue at an end of the range is replaced by itself (p.(Val559_Glu561del), not p.(Val559_Glu562delinsGlu))");
    return;
  }
  const auto& l = v.location;
  const bool exact = !l.start.uncertain && !l.end.uncertain;
  if (genomic(type)) {
    for (const Position* p : {&l.start.first, &l.start.last, &l.end.first, &l.end.last})
      if (p->utr || p->offset_sign) add(S::error, std::string("a ") + static_cast<char>(type) + ". position takes no '-', '*' or offset (numbering.md)");
  }
  if (l.range && exact) {
    const int c = compare(l.start.first, l.end.first);
    if (c == 0) add(S::error, "a range names two different positions (123_126, not 123_123)");
    else if (c == 1 && type != Coordinate::o && type != Coordinate::m)
      add(S::error, "a range is written 5′ to 3′ (123_126, not 126_123)");
  }
  switch (v.kind) {
    case Kind::substitution:
      if (l.range && exact) add(S::error, "a substitution replaces one nucleotide; two or more consecutive are a delins (c.79_80delinsTT)");
      if (!v.mixture && v.was == v.now) add(S::error, "a substitution by the same nucleotide is written '=' (c.123=, not c.123C>C)");
      break;
    case Kind::insertion:
      if (!l.range && !l.start.uncertain) add(S::error, "an insertion names its two flanking positions (g.123_124insG, not g.123insG)");
      if (l.range && exact) {
        const auto adjacent = next_position(l.start.first, l.end.first);
        if (adjacent && !*adjacent) add(S::error, "an insertion lies between two adjacent positions (123_124, not 123_125)");
        if (v.inserted.size() == 1 && v.inserted[0].kind == Piece::Kind::location && !v.inserted[0].inverted && v.inserted[0].location.range &&
            !v.inserted[0].location.end.uncertain && compare(v.inserted[0].location.end.first, l.start.first) == 0)
          add(S::error, "a copy inserted directly 3′ of its original is a duplication (g.123_456dup, not g.456_457ins123_456)");
      }
      break;
    case Kind::inversion:
      if (!l.range) add(S::error, "an inversion covers more than one nucleotide; one is a substitution (g.234inv is not allowed)");
      break;
    case Kind::deletion:
    case Kind::duplication:
      if (!v.stated.empty() && !v.stated_count && exact)
        add(S::note, std::string("the ") + (v.kind == Kind::deletion ? "deleted" : "duplicated") + " sequence need not be written (g.123_125del, not g.123_125delAGC)");
      break;
    case Kind::repeat:
      if (type == Coordinate::r && l.range && !v.units.empty())
        add(S::warning, "RNA/repeated.md calls a range with a unit redundant (r.-125_-123cug[4]), while the same page's proposal allows only that form (r.123_191cag[23])");
      if ((type == Coordinate::c || type == Coordinate::r) && !v.units.empty()) {
        // A coding repeat's unit is a multiple of three (DNA/repeated.md); outside the CDS this is not known from the text alone.
        bool in_cds = true;
        for (const Position* p : {&l.start.first, &l.end.first}) in_cds = in_cds && !p->utr && !p->offset_sign;
        for (const auto& u : v.units)
          if (in_cds && u.sequence.size() % 3 != 0)
            add(S::error, "in a coding sequence a repeat's unit is a multiple of three nucleotides; write dup or ins (c.2692_2693dup, not c.2686A[10])");
      }
      break;
    default: break;
  }
}

// The bases (or residues) a variant names, as a closed interval of a linear numbering, when the text
// decides it: positions without offsets or uncertainty. `between`: an insertion's two flanks.
struct Span {
  std::int64_t low = 0, high = 0;
  bool between = false;
};
inline std::optional<Span> span_of(const Variant& v) {
  const auto linear = [](const Position& p) -> std::optional<std::int64_t> {
    if (p.kind != Position::Kind::number || p.offset_sign) return std::nullopt;
    return p.utr == '-' ? -p.number : p.utr == '*' ? (std::int64_t{1} << 40) + p.number : p.number;
  };
  Span s;
  if (v.protein) {
    if (v.residues.start.uncertain || v.residues.end.uncertain) return std::nullopt;
    s.low = v.residues.start.first.number;
    s.high = v.residues.end.first.number;
  } else {
    if (v.location.start.uncertain || v.location.end.uncertain) return std::nullopt;
    const auto a = linear(v.location.start.first), b = linear(v.location.end.first);
    if (!a || !b) return std::nullopt;
    s.low = *a;
    s.high = *b;
  }
  s.between = v.kind == Kind::insertion;
  return s;
}
inline bool touch(const Span& x, const Span& y) {
  if (x.between && y.between) return x.low == y.low && x.high == y.high;
  if (x.between || y.between) return y.low <= x.high && y.high >= x.low;
  return y.low <= x.high + 1 && y.high + 1 >= x.low;
}

// Conditions on an allele's variants together (DNA/alleles.md, DNA/substitution.md).
inline void check_allele(const Allele& a, Coordinate type, std::vector<Diagnostic>& out) {
  std::vector<const Variant*> cis;
  for (std::size_t i = 0; i < a.chains.size(); i++) {
    if (i > 0 && i - 1 < a.separators.size() && a.separators[i - 1] != ';') { cis.clear(); }
    for (const auto& it : a.chains[i].items)
      if (it.choices.size() == 1 && !it.mixed) {
        // Small changes only: the standard keeps an inversion or a repeat apart from a change at its
        // breakpoint (DNA/complex.md: g.[776788_93191545inv;93191546T>C]).
        const Kind k = it.choices[0].kind;
        if (k == Kind::substitution || k == Kind::deletion || k == Kind::insertion || k == Kind::delins || k == Kind::duplication) cis.push_back(&it.choices[0]);
      }
  }
  for (std::size_t i = 0; i < cis.size(); i++)
    for (std::size_t j = i + 1; j < cis.size(); j++) {
      const auto x = span_of(*cis[i]), y = span_of(*cis[j]);
      if (!x || !y) continue;
      const std::string where = format_variant(*cis[i]) + ";" + format_variant(*cis[j]);
      if (touch(*x, *y)) {
        out.push_back({Diagnostic::Severity::error, "variants in cis that touch are one delins (c.79_80delinsTT, not c.[79G>T;80C>T])", where});
      } else if (type == Coordinate::c && cis[i]->kind == Kind::substitution && cis[j]->kind == Kind::substitution && !x->between && x->low == x->high &&
                 y->low == y->high && x->low > 0 && y->low > 0 && x->low < (std::int64_t{1} << 40) && y->low < (std::int64_t{1} << 40) &&
                 std::abs(x->low - y->low) == 2 && (x->low - 1) / 3 == (y->low - 1) / 3) {
        out.push_back({Diagnostic::Severity::error, "two changes in one codon, one nucleotide apart, are one delins (c.145_147delinsTGG, not c.[145C>T;147C>G])", where});
      }
    }
}

}  // namespace detail

inline std::vector<Diagnostic> check(const Hgvs& h) {
  std::vector<Diagnostic> out;
  const auto list = [&](const Description& d) {
    if (d.adjoined) return;
    const Coordinate t = d.variants.coordinate;
    if ((t == Coordinate::c || t == Coordinate::r) && d.has_reference && !d.reference.has_selector) {
      // Intronic positions need a genomic reference (numbering.md).
      bool intronic = false;
      const auto scan = [&](const Chain& ch) {
        for (const auto& it : ch.items)
          for (const auto& v : it.choices)
            for (const Position* p : {&v.location.start.first, &v.location.end.first, &v.location.start.last, &v.location.end.last}) intronic = intronic || p->offset_sign;
      };
      for (const auto& a : d.variants.alternatives)
        for (const auto& g : a.groups) {
          scan(g.chain);
          for (const auto& al : g.alleles)
            for (const auto& ch : al.chains) scan(ch);
        }
      const std::string acc = d.reference.accession;
      if (intronic && (acc.rfind("NM_", 0) == 0 || acc.rfind("NR_", 0) == 0 || acc.rfind("XM_", 0) == 0 || acc.rfind("ENST", 0) == 0))
        out.push_back({Diagnostic::Severity::warning, "a transcript holds no introns: intronic positions need a genomic reference, as NC_000023.10(NM_004006.2):c.357+1G>A (numbering.md; the standard's own repeated.md writes NM_000492.3:c.1210-12_1210-6T[7])", d.reference.text()});
    }
    if (t == Coordinate::p) {
      // One-letter codes are allowed; three-letter ones are preferred (general.md).
    }
    const auto chain = [&](const Chain& ch) {
      for (const auto& it : ch.items)
        for (const auto& v : it.choices) detail::check_variant(v, t, out);
    };
    for (const auto& a : d.variants.alternatives) {
      for (const auto& g : a.groups) {
        chain(g.chain);
        for (const auto& al : g.alleles) {
          for (const auto& ch : al.chains) chain(ch);
          detail::check_allele(al, t, out);
        }
        // "=" is written on an allele only when it holds no variant (DNA/alleles.md).
        if (g.alleles.size() >= 2) {
          bool every = true, some_identity = false;
          for (const auto& al : g.alleles) {
            bool variant = false;
            for (const auto& ch : al.chains)
              for (const auto& it : ch.items)
                for (const auto& v : it.choices) {
                  if (v.kind == Kind::identity) some_identity = true;
                  else variant = true;
                }
            every = every && variant;
          }
          if (every && some_identity)
            out.push_back({Diagnostic::Severity::error, "no change is written on an allele only when it holds no variant (c.[2376G>C];[3103del], not c.[2376G>C;3103=];[2376=;3103del])", format_list(d.variants)});
        }
      }
      if (a.groups.size() == 1 && a.groups[0].bracketed && a.groups[0].alleles.size() == 1 && a.groups[0].alleles[0].whole.empty() &&
          a.groups[0].alleles[0].chains.size() == 1 && d.variants.alternatives.size() == 1)
        out.push_back({Diagnostic::Severity::warning, "one allele is named; write the other too (c.[76A>C];[76=]), or leave out the brackets", format_list(d.variants)});
    }
  };
  for (const auto& d : h.alternatives) list(d);
  for (const auto& a : h.alleles)
    for (const auto& d : a) list(d);
  return out;
}

// ---- parsing ----------------------------------------------------------------------------------------------------
struct HgvsParse {
  bool accepted = false;
  Hgvs value;
  ParseNode tree;
  std::size_t furthest = 0;
  std::vector<std::string> expected;
  std::string error;
  std::vector<Diagnostic> diagnostics;  // conditions of meaning, when accepted
  std::string canonical;                // format(value), when accepted
  bool ok() const {
    if (!accepted) return false;
    for (const auto& d : diagnostics)
      if (d.severity == Diagnostic::Severity::error) return false;
    return true;
  }
};

inline HgvsParse parse_hgvs(std::string_view text) {
  HgvsParse out;
  const auto p = hgvs_grammar().parse(text);
  out.accepted = p.accepted;
  out.furthest = p.furthest;
  out.expected = p.expected;
  out.error = p.error;
  if (!p.accepted) return out;
  out.tree = p.tree;
  out.value = hgvs_of(p.tree);
  out.diagnostics = check(out.value);
  out.canonical = format(out.value);
  return out;
}

// ---- in plain words ------------------------------------------------------------------------------------------------
namespace detail {

inline std::string plural(std::int64_t n, const char* one, const char* many) { return std::to_string(n) + " " + (n == 1 ? one : many); }

inline std::string place(const Variant& v, Coordinate t) {
  const std::string type(1, static_cast<char>(t));
  if (v.protein) {
    const auto& l = v.residues;
    const auto at = [](const Point<Residue>& p) {
      if (!p.uncertain) return std::string(three_letter(p.first.residue)) + std::to_string(p.first.number);
      return std::string("somewhere in ") + three_letter(p.first.residue) + std::to_string(p.first.number) + (p.span ? std::string("–") + three_letter(p.last.residue) + std::to_string(p.last.number) : "");
    };
    return l.range ? at(l.start) + " to " + at(l.end) : at(l.start);
  }
  const auto& l = v.location;
  const auto at = [&](const Point<Position>& p) {
    if (!p.uncertain) return type + "." + format_position(p.first);
    if (!p.span) return "about " + type + "." + format_position(p.first);
    const bool u1 = p.first.kind == Position::Kind::unknown, u2 = p.last.kind == Position::Kind::unknown;
    if (u1 && u2) return std::string("an unknown position");
    if (u1) return "somewhere before " + type + "." + format_position(p.last);
    if (u2) return "somewhere after " + type + "." + format_position(p.first);
    return "somewhere between " + type + "." + format_position(p.first) + " and " + type + "." + format_position(p.last);
  };
  return l.range ? at(l.start) + " to " + at(l.end) : at(l.start);
}

inline std::int64_t span_length(const Variant& v) {
  const auto& l = v.location;
  if (!l.range) return 1;
  const Position &a = l.start.first, &b = l.end.first;
  if (l.start.uncertain || l.end.uncertain || a.kind != Position::Kind::number || b.kind != Position::Kind::number) return -1;
  if (a.offset_sign || b.offset_sign) {
    if (a.number == b.number && a.utr == b.utr) return b.signed_offset() - a.signed_offset() + 1;
    return -1;
  }
  if (a.utr != b.utr) return -1;
  return (a.utr == '-' ? a.number - b.number : b.number - a.number) + 1;
}

inline std::string reading_of(const Variant& v, Coordinate t) {
  const std::string where = place(v, t);
  const std::int64_t n = v.protein ? (v.residues.range ? v.residues.end.first.number - v.residues.start.first.number + 1 : 1) : span_length(v);
  const char* unit = v.protein ? "residue" : "base";
  const char* units = v.protein ? "residues" : "bases";
  const auto count = [&] { return n > 0 ? plural(n, unit, units) + ", " : std::string(); };
  std::string s;
  switch (v.kind) {
    case Kind::substitution:
      if (v.protein) {
        std::string alts;
        for (std::size_t i = 0; i < v.now.size(); i++) alts += (i ? " or " : "") + std::string(v.now[i] == '*' ? "a stop" : three_letter(v.now[i]));
        s = (v.now == "*" ? "a stop in place of " : "a substitution of ") + where + (v.now == "*" ? "" : " by " + alts);
      } else {
        s = "a substitution of " + v.was + " by " + v.now + " at " + where;
      }
      break;
    case Kind::deletion: s = "a deletion of " + count() + where; break;
    case Kind::duplication: s = "a duplication of " + count() + where; break;
    case Kind::insertion:
    case Kind::delins: {
      std::string what;
      if (v.protein) {
        what = v.stop_length ? plural(v.stop_length - 1, "residue", "residues") + " and a stop"
               : v.inserted_range ? "a copy of " + std::string(three_letter(v.inserted_from.start.first.residue)) + std::to_string(v.inserted_from.start.first.number) + " to " + three_letter(v.inserted_from.end.first.residue) + std::to_string(v.inserted_from.end.first.number)
               : v.inserted_count ? "unknown residues" : plural(static_cast<std::int64_t>(v.residues_inserted.size()), "residue", "residues");
      } else {
        for (std::size_t i = 0; i < v.inserted.size(); i++) {
          const Piece& p = v.inserted[i];
          if (i) what += ", then ";
          if (p.kind == Piece::Kind::sequence) {
            if (p.count && p.count->kind == Copies::Kind::exact) what += plural(p.count->low * static_cast<std::int64_t>(p.sequence.size()), p.sequence == "N" ? "unknown base" : "base", p.sequence == "N" ? "unknown bases" : "bases");
            else if (p.count) what += "an uncertain number of " + std::string(p.sequence == "N" ? "unknown " : "") + "bases";
            else what += p.sequence.size() <= 12 ? p.sequence : plural(static_cast<std::int64_t>(p.sequence.size()), "base", "bases");
          } else {
            what += (p.inverted ? "an inverted copy of " : "a copy of ") + std::string(p.kind == Piece::Kind::remote ? p.reference.text() + " " : "") +
                    format_location(p.location, format_position);
          }
        }
      }
      s = v.kind == Kind::insertion ? "an insertion of " + what + (v.protein || v.location.range ? " between " : " at ") + where
                                    : "a deletion of " + count() + where + ", replaced by " + what;
      break;
    }
    case Kind::inversion: s = "an inversion of " + count() + where; break;
    case Kind::identity: s = "no change at " + where; break;
    case Kind::repeat: {
      std::string units_text;
      for (const auto& u : v.units) units_text += (units_text.empty() ? "" : ", then ") + u.sequence + " " + format_copies(u.copies);
      s = v.units.empty() ? "the segment " + where + " in " + format_copies(*v.copies) + " copies" : "a repeated sequence at " + where + ": " + units_text;
      for (const auto& o : v.others) s += "; on the other allele " + format_copies(o);
      break;
    }
    case Kind::supernumerary: s = "a supernumerary copy of " + where; break;
    case Kind::methylation: s = (v.methylation == "gom" ? "a gain of methylation over " : v.methylation == "lom" ? "a loss of methylation over " : "normal methylation over ") + where; break;
    case Kind::frameshift:
      s = "a frameshift from " + where + (v.fs_new.empty() ? std::string() : std::string(", which becomes ") + three_letter(v.fs_new[0])) +
          (v.has_stop ? (v.stop_unknown ? ", with no new stop" : ", with a new stop at codon " + std::to_string(v.stop_at) + " of the new frame") : "");
      break;
    case Kind::extension:
      s = "the stop at " + where + " lost, read as " + three_letter(v.now.empty() ? 'X' : v.now[0]) +
          (v.stop_unknown ? ", the protein extended with no new stop" : ", the protein extended to a new stop " + std::to_string(v.stop_at) + " codons on");
      break;
    case Kind::n_extension: s = "a new start " + std::to_string(v.upstream) + " codons upstream of " + where + ", extending the protein's N-terminus"; break;
    case Kind::unknown: s = "a change at " + where + " whose consequence cannot be predicted"; break;
  }
  if (v.mixture) s += v.mixture == 1 ? ", in some cells of the individual (mosaic)" : ", in cells of another zygote (chimeric)";
  return s;
}

}  // namespace detail

// The description read out in plain words, one clause per variant: "a deletion of 3 bases, c.1521 to
// c.1523 of NM_000492.4".
inline std::string reading(const Hgvs& h) {
  std::string out;
  const auto desc = [&](const Description& d) {
    std::string s;
    if (d.adjoined) return "an adjoined transcript: " + d.five.reference.text() + " joined to " + d.three.reference.text();
    const Coordinate t = d.variants.coordinate;
    const std::string of = d.has_reference ? " of " + d.reference.text() : std::string();
    if (!d.variants.whole.empty()) {
      const std::string& w = d.variants.whole;
      s = w == "=" ? "no change" : w == "?" ? "an unknown effect" : w == "0" ? "no product" : w == "0?" ? "probably no product" : w == "(=)" ? "no change expected" : w == "spl" ? "splicing affected" : w == "spl?" ? "splicing possibly affected" : "no change expected beyond the DNA's";
      return s + of;
    }
    for (const auto& alt : d.variants.alternatives) {
      if (!s.empty()) s += "; or ";
      for (std::size_t gi = 0; gi < alt.groups.size(); gi++) {
        const Group& g = alt.groups[gi];
        if (gi) s += "; phase unknown with ";
        const auto chain = [&](const Chain& c, bool predicted) {
          std::string x;
          for (const auto& it : c.items) {
            if (!x.empty()) x += " joined to ";
            for (std::size_t k = 0; k < it.choices.size(); k++) x += (k ? " or " : "") + detail::reading_of(it.choices[k], t);
            if (it.predicted || predicted) x += " (predicted)";
          }
          return x;
        };
        if (!g.bracketed) {
          s += chain(g.chain, false);
          continue;
        }
        for (std::size_t ai = 0; ai < g.alleles.size(); ai++) {
          const Allele& a = g.alleles[ai];
          s += g.alleles.size() > 1 ? (ai ? "; and on another allele: " : "on one allele: ") : "on one allele: ";
          if (!a.whole.empty()) { s += a.whole == "=" ? "no change" : a.whole == "?" ? "a variant not identified" : a.whole == "0" ? "no allele" : a.whole; continue; }
          for (std::size_t k = 0; k < a.chains.size(); k++) {
            if (k) s += a.separators[k - 1] == ',' ? ", or in another transcript " : a.separators[k - 1] == '?' ? ", phase unknown, " : ", and ";
            s += chain(a.chains[k], a.predicted);
          }
        }
      }
    }
    return s + of;
  };
  for (const auto& d : h.alternatives) out += (out.empty() ? "" : "; or ") + desc(d);
  for (std::size_t i = 0; i < h.alleles.size(); i++)
    for (const auto& d : h.alleles[i]) out += (out.empty() ? "" : "; ") + desc(d);
  return out;
}

}  // namespace genomics::nomenclature
