// Results described in ISCN: copy-number segments in coordinates as a microarray description
// (arr[GRCh38] 1q21.1q21.2(147,105,904_147,917,509)x1), and the junctions of a genome found by
// sequencing as a sequencing description in the short and the detailed systems (seq[GRCh38]
// t(1;8)(q23.3;q12.2), and 1pter→1q23.3(164,789,896)::8q12.2(61,760,650)→8qter;…) with each
// junction as Ordulu et al. write a rearrangement sequence (1q23.3(+)(164,789,896)::8q12.2(+)
// (61,760,650)).
//
// - Microarrays (ISCN 2024 ch. 8; ISCN 2020 ch. 14, as summarised by Miller and Madan 2021, and
//   its examples in Conlin et al. 2025): each departure from the normal copy number is written as
//   its chromosome, the bands at its two ends (one when they are the same), its first and last
//   bases (1-based, with commas grouping thousands, separated by an underscore) and its copies
//   (x1, x3, x2~3 for a mosaic range, [0.3] for a proportion of cells); hmz marks a run of
//   homozygosity; a whole chromosome is (13)x3, chromosomes with the same copies are grouped,
//   (1-20,22)x3; the sex chromosomes come first; a normal result is arr(X,1-22)x2 for a female and
//   arr(X,Y)x1,(1-22)x2 for a male. The bands are those holding the first and the last base.
// - Sequencing (Ordulu et al. 2014, adopted by ISCN 2016–2024): an adjacency joins the piece of one
//   chromosome on one side of a break to the piece of another (or the same) on one side of
//   another. The pieces the breaks cut the chromosomes into, joined by the adjacencies and, where a
//   break is joined on neither side, by the reference itself, form a graph of degree at most two;
//   its walks from the telomeres (and its cycles, which are rings) are the derivative chromosomes.
//   A derivative is named by the centromere it holds and read from pter to qter of its centric
//   segment. A balanced reciprocal translocation, an interstitial deletion, a tandem duplication,
//   an inversion and an unbalanced derivative are recognised and named; anything else is a complex
//   rearrangement, (2,3,8)cx, written chromosome by chromosome.
//
// Positions inside are 0-based boundaries between bases: a break at p separates base p − 1 from
// base p (0-based), so the left piece ends at base p and the right piece begins at base p + 1 in
// ISCN's 1-based numbering.
//
// Complexity: O(k log k) in the segments or the adjacencies.

#pragma once

#include <genomics/karyotype/bands.hpp>
#include <genomics/karyotype/format.hpp>
#include <genomics/karyotype/interpret.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace genomics::karyotype {

// ---- microarrays ---------------------------------------------------------------------------------------

struct CopySegment {
  int chromosome = -1;
  std::uint64_t start = 0, end = 0;  // 0-based, half-open
  int copies = 2, copies_hi = -1;    // a range of copies (x2~3) when copies_hi > copies
  double proportion = -1;            // the proportion of cells, when known
  bool homozygous = false;           // a copy-neutral run of homozygosity (hmz)
  std::string inheritance;           // dn, mat, pat, dmat, dpat, …
  CopySegment() = default;
  CopySegment(int c, std::uint64_t s, std::uint64_t e, int k, int k_hi = -1, double p = -1)
      : chromosome(c), start(s), end(e), copies(k), copies_hi(k_hi), proportion(p) {}
};

struct ArrOptions {
  std::string build = "GRCh38";
  bool commas = true;
  bool unicode = false;
  std::optional<bool> male;          // the sex of the sample, for the normal copies of X and Y
};

inline int normal_copies(int c, std::optional<bool> male) {
  if (c == chr_x) return male ? (*male ? 1 : 2) : -1;
  if (c == chr_y) return male ? (*male ? 1 : 0) : -1;
  return 2;
}

// Writes copy-number segments as a microarray description. Segments at their normal copy number
// (and not homozygous) are left out; a segment covering a whole chromosome is written as one.
inline std::string describe_arr(std::vector<CopySegment> segs, const Ideogram& g, const ArrOptions& o = {}) {
  std::vector<ArrEntry> entries;
  std::map<std::tuple<int, int, std::string, std::string>, std::vector<int>> wholes;  // (copies, hi, flags, inheritance) → chromosomes
  for (const auto& s : segs) {
    if (!g.has(s.chromosome) || s.end <= s.start) continue;
    const int normal = normal_copies(s.chromosome, o.male);
    if (s.copies == normal && s.copies_hi <= s.copies && !s.homozygous) continue;
    const bool whole = s.start == 0 && s.end >= g.length[s.chromosome];
    if (whole && !s.homozygous && s.proportion < 0) {
      wholes[{s.copies, s.copies_hi, "", s.inheritance}].push_back(s.chromosome);
      continue;
    }
    ArrEntry e;
    e.kind = ArrEntry::Kind::region;
    e.chromosome = s.chromosome;
    const BandRange r = g.bands_over(s.chromosome, s.start, std::min(s.end, g.length[s.chromosome]));
    e.arm1 = g.bands[r.first].arm;
    e.band1 = g.bands[r.first].digits;
    if (r.last != r.first) {
      e.arm2 = g.bands[r.last].arm;
      e.band2 = g.bands[r.last].digits;
    }
    e.start = s.start + 1;
    e.end = s.end;
    e.copies = s.copies;
    e.copies_hi = s.copies_hi;
    e.proportion = s.proportion;
    if (s.homozygous) e.flags.push_back("hmz");
    if (!s.inheritance.empty()) e.flags.push_back(s.inheritance);
    entries.push_back(std::move(e));
  }
  for (auto& [key, cs] : wholes) {
    std::sort(cs.begin(), cs.end(), [](int a, int b) { return iscn_rank(a) < iscn_rank(b); });
    ArrEntry e;
    e.kind = ArrEntry::Kind::whole;
    for (std::size_t i = 0; i < cs.size();) {
      std::size_t j = i;
      while (j + 1 < cs.size() && !is_sex(cs[j]) && cs[j + 1] == cs[j] + 1) j++;
      e.wholes.push_back({cs[i], cs[j]});
      i = j + 1;
    }
    e.copies = std::get<0>(key);
    e.copies_hi = std::get<1>(key);
    if (!std::get<3>(key).empty()) e.flags.push_back(std::get<3>(key));
    entries.push_back(std::move(e));
  }
  Technique t;
  t.name = "arr";
  t.build = o.build;
  if (entries.empty()) {
    // A normal result.
    ArrEntry e;
    e.kind = ArrEntry::Kind::whole;
    e.copies = 2;
    if (o.male && *o.male) {
      ArrEntry xy;
      xy.kind = ArrEntry::Kind::whole;
      xy.wholes = {{chr_x, chr_x}, {chr_y, chr_y}};
      xy.copies = 1;
      e.wholes = {{0, 21}};
      t.arr = {xy, e};
    } else if (o.male) {
      e.wholes = {{chr_x, chr_x}, {0, 21}};
      t.arr = {e};
    } else {
      e.wholes = {{0, 21}};
      t.arr = {e};
    }
    t.build.clear();  // ISCN writes a normal result without its build
  } else {
    std::stable_sort(entries.begin(), entries.end(), arr_before);
    t.arr = std::move(entries);
  }
  FormatOptions f;
  f.commas = o.commas;
  f.unicode = o.unicode;
  f.order = false;
  return format(t, f);
}

// ---- sequencing ----------------------------------------------------------------------------------------

struct Adjacency {
  // The piece of chromosome a ending at break_a ([0, break_a), a_left) or starting there
  // ([break_a, L)), joined to the piece of chromosome b starting at break_b ([break_b, L),
  // b_right) or ending there ([0, break_b)); nontemplated bases between them.
  int chromosome_a = -1;
  std::uint64_t break_a = 0;
  bool a_left = true;
  int chromosome_b = -1;
  std::uint64_t break_b = 0;
  bool b_right = true;
  std::string bases;
  Adjacency() = default;
  Adjacency(int a, std::uint64_t at_a, bool left, int b, std::uint64_t at_b, bool right, std::string nt = {})
      : chromosome_a(a), break_a(at_a), a_left(left), chromosome_b(b), break_b(at_b), b_right(right), bases(std::move(nt)) {}
};

struct SeqOptions {
  std::string build = "GRCh38";
  bool unicode = true;  // → and −, as ISCN prints them
};

struct SeqDescription {
  std::string kind;        // t, der, del, dup, inv, cx, or empty when nothing departs from the reference
  std::string short_form;  // seq[GRCh38] t(1;8)(q23.3;q12.2)
  std::string detailed;    // seq[GRCh38] t(1;8)(1pter→1q23.3(164,789,896)::…)
  std::vector<std::string> junctions;  // 1q23.3(+)(164,789,896)::8q12.2(+)(61,760,650)
  std::vector<Chromosome> derivatives; // the walks with a centromere and a junction, read from pter
  std::vector<Chromosome> acentric;    // pieces left without a centromere
  std::vector<Diagnostic> diagnostics;
};

namespace detail {

inline std::string band_of_base(const Ideogram& g, int c, std::uint64_t base0) {
  return g.bands[g.band_at(c, std::min(base0, g.length[c] ? g.length[c] - 1 : 0))].name();
}

}  // namespace detail

// One rearrangement: junctions that belong together (describe_seq groups them).
inline SeqDescription describe_seq_event(const std::vector<Adjacency>& adjs, const Ideogram& g, const SeqOptions& o = {}) {
  SeqDescription out;
  std::vector<Diagnostic> sink;
  detail::Interpreter in(g, sink);
  const std::string minus = o.unicode ? "\xE2\x88\x92" : "-";

  const std::string head = "seq" + (o.build.empty() ? std::string() : "[" + o.build + "]") + " ";
  LongFormOptions lo1;
  lo1.unicode = o.unicode;
  lo1.chromosome_names = false;
  const auto point_on = [&](int c, std::uint64_t pos, bool lower) {
    if (pos == 0) return in.pter(c);
    if (pos == g.length[c]) return in.qter(c);
    Point p;
    p.kind = Point::Kind::band;
    p.chromosome = c;
    p.at = p.lo = p.hi = pos;
    p.written_position = true;
    p.label = detail::band_of_base(g, c, lower ? pos : pos - 1);
    return p;
  };
  // Rearrangements within one chromosome whose pieces a graph of single copies cannot hold (a tandem
  // duplication uses its piece twice) or that read most plainly as their own words.
  if (!adjs.empty() && g.has(adjs[0].chromosome_a)) {
    const int c = adjs[0].chromosome_a;
    const bool one = std::all_of(adjs.begin(), adjs.end(), [&](const Adjacency& a) { return a.chromosome_a == c && a.chromosome_b == c; });
    std::string kind;
    detail::Word w;
    if (one && adjs.size() == 1 && adjs[0].a_left && adjs[0].b_right && adjs[0].break_a != adjs[0].break_b) {
      const std::uint64_t x = adjs[0].break_a, y = adjs[0].break_b;
      kind = y > x ? "del" : "dup";
      w = {in.ref(c, point_on(c, 0, true), point_on(c, x, false)), in.ref(c, point_on(c, y, true), point_on(c, g.length[c], false))};
      if (kind == "del") out.short_form = head + "del(" + chromosome_name(c) + ")(" + w[0].to.label + w[1].from.label + ")";
      else out.short_form = head + "dup(" + chromosome_name(c) + ")(" + w[1].from.label + w[0].to.label + ")";
    } else if (one && adjs.size() == 2) {
      // inv(c)(x y): x left joined to y left, and x right to y right.
      const Adjacency* l = nullptr;
      const Adjacency* r = nullptr;
      for (const auto& a : adjs) {
        if (a.a_left && !a.b_right) l = &a;
        if (!a.a_left && a.b_right) r = &a;
      }
      if (l && r && l->break_a == r->break_a && l->break_b == r->break_b) {
        const std::uint64_t x = std::min(l->break_a, l->break_b), y = std::max(l->break_a, l->break_b);
        kind = "inv";
        w = {in.ref(c, point_on(c, 0, true), point_on(c, x, false)), in.ref(c, point_on(c, y, false), point_on(c, x, true)),
             in.ref(c, point_on(c, y, true), point_on(c, g.length[c], false))};
        out.short_form = head + "inv(" + chromosome_name(c) + ")(" + w[0].to.label + w[2].from.label + ")";
      }
    }
    if (!kind.empty()) {
      Chromosome ch;
      ch.segments = w;
      ch.centromere_of = c;
      ch.abnormal = true;
      ch.label = "der(" + chromosome_name(c) + ")";
      out.kind = kind;
      out.detailed = head + kind + "(" + chromosome_name(c) + ")(" + long_form(ch, lo1) + ")";
      out.derivatives.push_back(ch);
      for (const auto& a : adjs) {
        const auto end = [&](int cc, std::uint64_t brk, bool plus, bool is_a) {
          const std::uint64_t base1 = (is_a == plus) ? brk : brk + 1;
          return chromosome_name(cc) + detail::band_of_base(g, cc, base1 - 1) + "(" + (plus ? "+" : minus) + ")(" + with_commas(base1) + ")";
        };
        out.junctions.push_back(end(a.chromosome_a, a.break_a, a.a_left, true) + "::" + end(a.chromosome_b, a.break_b, a.b_right, false));
      }
      return out;
    }
  }

  // The breaks of each chromosome, and its pieces.
  std::map<int, std::set<std::uint64_t>> breaks;
  for (const auto& a : adjs) {
    if (!g.has(a.chromosome_a) || !g.has(a.chromosome_b)) {
      out.diagnostics.push_back({Diagnostic::Level::error, 0, 0, "an adjacency names a chromosome the ideogram does not hold"});
      return out;
    }
    if (a.break_a == 0 || a.break_a >= g.length[a.chromosome_a] || a.break_b == 0 || a.break_b >= g.length[a.chromosome_b]) {
      out.diagnostics.push_back({Diagnostic::Level::error, 0, 0, "an adjacency's break lies at or beyond a telomere"});
      return out;
    }
    breaks[a.chromosome_a].insert(a.break_a);
    breaks[a.chromosome_b].insert(a.break_b);
  }
  struct Piece { int c; std::uint64_t start, end; };
  std::vector<Piece> pieces;
  std::map<std::pair<int, std::uint64_t>, std::size_t> starting, ending;  // (c, position) → piece
  for (auto& [c, bs] : breaks) {
    std::uint64_t at = 0;
    for (std::uint64_t b : bs) {
      ending[{c, b}] = pieces.size();
      starting[{c, at}] = pieces.size();
      pieces.push_back({c, at, b});
      at = b;
    }
    ending[{c, g.length[c]}] = pieces.size();
    starting[{c, at}] = pieces.size();
    pieces.push_back({c, at, g.length[c]});
  }
  // Ends: 2i is the left end of piece i, 2i + 1 its right end.
  std::vector<long> join(pieces.size() * 2, -1);
  std::vector<const Adjacency*> join_adj(pieces.size() * 2, nullptr);
  for (const auto& a : adjs) {
    const long ea = a.a_left ? static_cast<long>(2 * ending[{a.chromosome_a, a.break_a}] + 1) : static_cast<long>(2 * starting[{a.chromosome_a, a.break_a}]);
    const long eb = a.b_right ? static_cast<long>(2 * starting[{a.chromosome_b, a.break_b}]) : static_cast<long>(2 * ending[{a.chromosome_b, a.break_b}] + 1);
    if (join[static_cast<std::size_t>(ea)] >= 0 || join[static_cast<std::size_t>(eb)] >= 0)
      out.diagnostics.push_back({Diagnostic::Level::warning, 0, 0, "an end is joined twice; the later junction is kept"});
    join[static_cast<std::size_t>(ea)] = eb;
    join[static_cast<std::size_t>(eb)] = ea;
    join_adj[static_cast<std::size_t>(ea)] = join_adj[static_cast<std::size_t>(eb)] = &a;
  }
  // The reference joins a break where neither of its sides is joined.
  std::vector<char> free_end(pieces.size() * 2, 0);
  for (auto& [c, bs] : breaks)
    for (std::uint64_t b : bs) {
      const std::size_t l = 2 * ending[{c, b}] + 1, r = 2 * starting[{c, b}];
      if (join[l] < 0 && join[r] < 0) { join[l] = static_cast<long>(r); join[r] = static_cast<long>(l); }
      else { if (join[l] < 0) free_end[l] = 1; if (join[r] < 0) free_end[r] = 1; }
    }

  // Walks.
  const auto point_at = [&](int c, std::uint64_t pos, bool lower) {
    Point p;
    if (pos == 0) return in.pter(c);
    if (pos == g.length[c]) return in.qter(c);
    p.kind = Point::Kind::band;
    p.chromosome = c;
    p.at = p.lo = p.hi = pos;
    p.written_position = true;
    p.label = detail::band_of_base(g, c, lower ? pos : pos - 1);
    return p;
  };
  std::vector<char> used(pieces.size(), 0);
  std::vector<Chromosome> walks;
  std::vector<int> junction_count;
  const auto walk_from = [&](std::size_t end, bool cycle) {
    detail::Word w;
    int junctions_here = 0;
    std::size_t e = end;
    Chromosome ch;
    ch.ring = cycle;
    if (free_end[e]) ch.open = ":";
    for (;;) {
      const std::size_t piece = e / 2;
      if (used[piece]) break;
      used[piece] = 1;
      const Piece& p = pieces[piece];
      const bool forward = e % 2 == 0;
      const Point a = point_at(p.c, forward ? p.start : p.end, forward);
      const Point b = point_at(p.c, forward ? p.end : p.start, !forward);
      w.push_back(in.ref(p.c, a, b));
      const std::size_t other = forward ? e + 1 : e - 1;
      const long next = join[other];
      if (next < 0) {
        if (free_end[other]) ch.close = ":";
        break;
      }
      if (join_adj[other]) {
        junctions_here++;
        if (!join_adj[other]->bases.empty()) {
          Segment s;
          s.kind = Segment::Kind::bases;
          s.bases = join_adj[other]->bases;
          w.push_back(s);
        }
      }
      e = static_cast<std::size_t>(next);
    }
    ch.segments = std::move(w);
    walks.push_back(std::move(ch));
    junction_count.push_back(junctions_here);
  };
  for (std::size_t i = 0; i < pieces.size(); i++)
    if (!used[i] && pieces[i].start == 0) walk_from(2 * i, false);
  for (std::size_t i = 0; i < pieces.size(); i++)
    if (!used[i] && pieces[i].end == g.length[pieces[i].c]) walk_from(2 * i + 1, false);
  for (std::size_t i = 0; i < pieces.size(); i++)
    if (!used[i] && (free_end[2 * i] || free_end[2 * i + 1])) walk_from(free_end[2 * i] ? 2 * i : 2 * i + 1, false);
  for (std::size_t i = 0; i < pieces.size(); i++)
    if (!used[i]) walk_from(2 * i, true);

  // Derivatives: walks holding a centromere and a junction, oriented and named.
  struct Derivative { Chromosome ch; int named; int junctions; };
  std::vector<Derivative> ders;
  for (std::size_t k = 0; k < walks.size(); k++) {
    Chromosome& ch = walks[k];
    int named = -1;
    for (const auto& s : ch.segments)
      if (s.kind == Segment::Kind::reference && s.start() < g.centromere[s.chromosome] && g.centromere[s.chromosome] < s.end()) { named = s.chromosome; break; }
    // A walk without a junction is the reference, or a remnant of a chromosome whose other
    // junction was not seen; only an acentric remnant is reported.
    if (junction_count[k] == 0 && !ch.ring) {
      bool centric = false;
      for (const auto& s : ch.segments)
        if (s.kind == Segment::Kind::reference && s.start() < g.centromere[s.chromosome] && g.centromere[s.chromosome] < s.end()) centric = true;
      if (!centric && (!ch.open.empty() || !ch.close.empty())) out.acentric.push_back(ch);
      continue;
    }
    detail::Word oriented = in.orient(ch.segments, named);
    const auto same = [](const Segment& x, const Segment& y) {
      return x.kind == y.kind && x.chromosome == y.chromosome && x.from.at == y.from.at && x.to.at == y.to.at;
    };
    const bool flip = !oriented.empty() && !same(oriented.front(), ch.segments.front());
    if (flip) {
      ch.segments = std::move(oriented);
      std::swap(ch.open, ch.close);
    }
    ch.centromere_of = named;
    ch.abnormal = true;
    if (named < 0) { out.acentric.push_back(ch); continue; }
    ders.push_back({ch, named, junction_count[k]});
  }
  std::stable_sort(ders.begin(), ders.end(), [](const Derivative& a, const Derivative& b) { return iscn_rank(a.named) < iscn_rank(b.named); });

  // The junctions as rearrangement sequences.
  for (const auto& a : adjs) {
    const auto end = [&](int c, std::uint64_t brk, bool plus, bool is_a) {
      // The a end reads toward its break, the b end away from it.
      const std::uint64_t base1 = (is_a == plus) ? brk : brk + 1;  // 1-based
      return chromosome_name(c) + detail::band_of_base(g, c, base1 - 1) + "(" + (plus ? "+" : minus) + ")(" + with_commas(base1) + ")";
    };
    out.junctions.push_back(end(a.chromosome_a, a.break_a, a.a_left, true) + "::" + (a.bases.empty() ? "" : a.bases + "::") +
                            end(a.chromosome_b, a.break_b, a.b_right, false));
  }

  LongFormOptions lo;
  lo.unicode = o.unicode;
  const auto seg_of = [](const Chromosome& ch, std::size_t i) -> const Segment& { return ch.segments[i]; };
  const auto band_at_break = [&](const Segment& s, bool from) {
    const Point& p = from ? s.from : s.to;
    return p.label;
  };
  (void)seg_of;
  if (ders.empty() && out.acentric.empty()) return out;

  // Recognised patterns.
  const auto refs = [](const Chromosome& ch) {
    std::vector<const Segment*> r;
    for (const auto& s : ch.segments)
      if (s.kind == Segment::Kind::reference) r.push_back(&s);
    return r;
  };
  if (ders.size() == 2 && ders[0].named != ders[1].named && ders[0].junctions == 1 && ders[1].junctions == 1) {
    const auto r0 = refs(ders[0].ch), r1 = refs(ders[1].ch);
    if (r0.size() == 2 && r1.size() == 2) {
      // der(a) = a's centric piece and b's; der(b) the reverse.
      const int a = ders[0].named, b = ders[1].named;
      const auto own = [&](const std::vector<const Segment*>& r, int c) { return r[0]->chromosome == c ? r[0] : r[1]; };
      const auto other = [&](const std::vector<const Segment*>& r, int c) { return r[0]->chromosome == c ? r[1] : r[0]; };
      if (other(r0, a)->chromosome == b && other(r1, b)->chromosome == a) {
        const Segment* sa = own(r0, a);
        const Segment* sb = own(r1, b);
        // The break on each chromosome: the end of its centric piece away from its telomere.
        const auto brk = [&](const Segment* s) { return s->from.kind == Point::Kind::pter || s->from.kind == Point::Kind::qter ? band_at_break(*s, false) : band_at_break(*s, true); };
        out.kind = "t";
        const std::string bands = "(" + brk(sa) + ";" + brk(sb) + ")";
        out.short_form = head + "t(" + chromosome_name(a) + ";" + chromosome_name(b) + ")" + bands;
        out.detailed = head + "t(" + chromosome_name(a) + ";" + chromosome_name(b) + ")(" + long_form(ders[0].ch, lo) + ";" + long_form(ders[1].ch, lo) + ")";
      }
    }
  }
  if (out.kind.empty() && ders.size() == 1 && out.acentric.size() <= 1) {
    const Derivative& d = ders[0];
    const auto r = refs(d.ch);
    const int c = d.named;
    const bool one = std::all_of(r.begin(), r.end(), [&](const Segment* s) { return s->chromosome == c; });
    lo.chromosome_names = !one;
    if (one && d.junctions == 1 && r.size() == 2 && !r[0]->reversed() && !r[1]->reversed()) {
      if (r[0]->end() < r[1]->start()) {
        out.kind = "del";
        out.short_form = head + "del(" + chromosome_name(c) + ")(" + r[0]->to.label + r[1]->from.label + ")";
      } else if (r[1]->start() < r[0]->end()) {
        out.kind = "dup";
        out.short_form = head + "dup(" + chromosome_name(c) + ")(" + r[1]->from.label + r[0]->to.label + ")";
      }
    } else if (one && d.junctions == 2 && r.size() == 3 && r[1]->reversed() && !r[0]->reversed() && !r[2]->reversed()) {
      out.kind = "inv";
      out.short_form = head + "inv(" + chromosome_name(c) + ")(" + r[0]->to.label + r[2]->from.label + ")";
    } else if (!one && d.junctions == 1 && r.size() == 2) {
      const Segment* mine = r[0]->chromosome == c ? r[0] : r[1];
      const Segment* theirs = r[0]->chromosome == c ? r[1] : r[0];
      const auto brk_of = [&](const Segment* s) { return (s->from.kind == Point::Kind::pter || s->from.kind == Point::Kind::qter) ? s->to.label : s->from.label; };
      out.kind = "der";
      FormatOptions f;
      Abnormality t;
      t.type = Type::t;
      t.chromosomes = {ChromosomeRef{c}, ChromosomeRef{theirs->chromosome}};
      Breakpoint x, y;
      x.kind = y.kind = Breakpoint::Kind::band;
      x.arm = brk_of(mine)[0];
      x.digits = *band_digits(brk_of(mine).substr(1));
      y.arm = brk_of(theirs)[0];
      y.digits = *band_digits(brk_of(theirs).substr(1));
      t.has_groups = true;
      t.groups = {{x}, {y}};
      out.short_form = head + "der(" + chromosome_name(c) + ")" + format(t, f);
    }
    if (!out.kind.empty()) out.detailed = head + (out.kind == "der" ? "der" : out.kind) + "(" + chromosome_name(c) + ")(" + long_form(d.ch, lo) + ")";
  }
  if (out.kind.empty()) {
    std::set<int> cs;
    for (const auto& d : ders)
      for (const auto& s : d.ch.segments)
        if (s.kind == Segment::Kind::reference) cs.insert(s.chromosome);
    std::vector<int> order(cs.begin(), cs.end());
    std::sort(order.begin(), order.end(), [](int a, int b) { return iscn_rank(a) < iscn_rank(b); });
    std::string list;
    for (std::size_t i = 0; i < order.size(); i++) list += (i ? "," : "") + chromosome_name(order[i]);
    out.kind = "cx";
    out.short_form = head + "(" + list + ")cx";
    out.detailed = head + "(" + list + ")cx";
    for (const auto& d : ders) out.detailed += ",der(" + chromosome_name(d.named) + ")(" + long_form(d.ch, lo) + ")";
  }
  for (auto& d : ders) {
    d.ch.label = "der(" + chromosome_name(d.named) + ")";
    out.derivatives.push_back(d.ch);
  }
  return out;
}

// A junction read from its other end is the same junction: the form kept joins a left piece or
// ends at a right one (a_left or b_right), and within one chromosome a junction of two left or two
// right pieces names its lower break first.
inline Adjacency normalised(Adjacency a) {
  if (!a.a_left && !a.b_right) {
    std::string rc;
    for (auto it = a.bases.rbegin(); it != a.bases.rend(); ++it)
      rc += *it == 'A' ? 'T' : *it == 'C' ? 'G' : *it == 'G' ? 'C' : *it == 'T' ? 'A' : 'N';  // read from the other strand
    a = Adjacency(a.chromosome_b, a.break_b, true, a.chromosome_a, a.break_a, true, rc);
  }
  if (a.chromosome_a == a.chromosome_b && a.a_left != a.b_right && a.break_a > a.break_b) {
    std::swap(a.break_a, a.break_b);
  }
  return a;
}

// Junctions of a sample described: grouped into rearrangements (junctions whose breaks lie within
// `tolerance` bases of one another on a chromosome belong to one: the two junctions of a reciprocal
// translocation or an inversion, with the small losses at their breaks), each described, and the
// descriptions joined in ISCN's order (seq[GRCh38] t(21;22)(…),del(21)(…)).
inline SeqDescription describe_seq(const std::vector<Adjacency>& raw, const Ideogram& g, const SeqOptions& o = {}, std::uint64_t tolerance = 1000) {
  std::vector<Adjacency> adjs;
  for (const auto& a : raw) adjs.push_back(normalised(a));
  const std::size_t n = adjs.size();
  std::vector<std::size_t> parent(n);
  for (std::size_t i = 0; i < n; i++) parent[i] = i;
  const auto root = [&](std::size_t i) {
    while (parent[i] != i) i = parent[i] = parent[parent[i]];
    return i;
  };
  const auto near = [&](int c, std::uint64_t x, int d, std::uint64_t y) { return c == d && (x > y ? x - y : y - x) <= tolerance; };
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 1; j < n; j++) {
      const Adjacency &a = adjs[i], &b = adjs[j];
      if (near(a.chromosome_a, a.break_a, b.chromosome_a, b.break_a) || near(a.chromosome_a, a.break_a, b.chromosome_b, b.break_b) ||
          near(a.chromosome_b, a.break_b, b.chromosome_a, b.break_a) || near(a.chromosome_b, a.break_b, b.chromosome_b, b.break_b))
        parent[root(i)] = root(j);
    }
  std::map<std::size_t, std::vector<Adjacency>> groups;
  for (std::size_t i = 0; i < n; i++) groups[root(i)].push_back(adjs[i]);
  if (groups.size() <= 1) return describe_seq_event(adjs, g, o);
  std::vector<SeqDescription> parts;
  for (auto& [k, v] : groups) parts.push_back(describe_seq_event(v, g, o));
  const auto first = [](const SeqDescription& d) {
    int best = 1000;
    for (const auto& ch : d.derivatives)
      for (const auto& s : ch.segments)
        if (s.kind == Segment::Kind::reference) best = std::min(best, iscn_rank(s.chromosome));
    return best;
  };
  std::stable_sort(parts.begin(), parts.end(), [&](const SeqDescription& a, const SeqDescription& b) { return first(a) < first(b); });
  SeqDescription out;
  const std::string head = "seq" + (o.build.empty() ? std::string() : "[" + o.build + "]") + " ";
  for (const auto& p : parts) {
    if (p.kind.empty()) continue;
    out.kind += (out.kind.empty() ? "" : ",") + p.kind;
    out.short_form += (out.short_form.empty() ? head : ",") + p.short_form.substr(std::min(p.short_form.size(), head.size()));
    out.detailed += (out.detailed.empty() ? head : ",") + p.detailed.substr(std::min(p.detailed.size(), head.size()));
    out.junctions.insert(out.junctions.end(), p.junctions.begin(), p.junctions.end());
    out.derivatives.insert(out.derivatives.end(), p.derivatives.begin(), p.derivatives.end());
    out.acentric.insert(out.acentric.end(), p.acentric.begin(), p.acentric.end());
    out.diagnostics.insert(out.diagnostics.end(), p.diagnostics.begin(), p.diagnostics.end());
  }
  return out;
}

}  // namespace genomics::karyotype
