// Karyotypes interpreted: every chromosome of a cell line as a word over oriented segments of the
// normal chromosomes, its copy number at every band against the normal complement, its losses,
// gains and fusions, and the detailed system written from the short (theory/karyotypes.md §3–5).
//
// - Points. A breakpoint written as a band (q21) lies somewhere in that band: it is a point with an
//   interval of uncertainty [lo, hi) (the band's extent on GRCh38, at the finest level of the
//   ideogram that holds it) and a representative position `at` within it (the middle; breakpoints
//   of one abnormality that share a band are spread through it in the order written). pter, qter
//   and the centromere (cen, p10, q10) are exact; a nucleotide position of a sequencing result is
//   exact. Breakpoints written once are one point wherever they recur: the two derivatives of
//   t(9;22) share the breaks at 9q34 and 22q11.2, which is what makes a reciprocal translocation
//   balanced exactly.
// - Words. A chromosome of the cell is a word of segments a→b of reference chromosomes (reversed
//   when b < a), of material of unknown origin (?), of homogeneously staining regions and of
//   nucleotides that align nowhere; a ring is a circular word. Each structural abnormality is an
//   operation on words made of the pieces a breakpoint cuts a chromosome into: the centric piece
//   prox(c, x) ending at x and the acentric piece dist(c, x) starting from it (at the centromere,
//   p10 and q10 name the arm that is the piece). A translocation of n chromosomes joins
//   der(c_k) = prox(c_k, x_k) · dist(c_{k−1}, x_{k−1}); a dicentric prox(a, x) · prox(b, y)⁻¹; an
//   inversion reverses a factor; a deletion removes one; a duplication repeats one, directly or
//   inverted; an insertion moves one; a derivative applies its abnormalities in order to the word
//   of its chromosome (ISCN 2024 §5; the detailed system's examples of ISCN 2009 quoted by
//   Warrender et al. 2019; Ordulu et al. 2014 for the orientation).
// - Orientation. A derivative is read so that the segment holding its centromere runs from pter
//   toward qter (Ordulu et al. 2014, "Strand direction"); when no segment holds a centromere
//   strictly (whole-arm rearrangements), from a pter if it has one, else from the chromosome named
//   first (ISCN 2009's t(1;3)(p10;p10), CytoGPS's der(2;13)(q10;q10)).
// - The cell. The normal complement is the ploidy (written as <3n>, or the level whose range
//   23k ± 11 holds the count, ISCN 1995 as quoted by Bakhshi et al. 2004) times the autosomes, and
//   the sex chromosomes written. A structural abnormality replaces the normal chromosomes it is
//   made from (one 9 and one 22 for t(9;22), one 22 for der(22)t(9;22), one 13 and one 21 for
//   der(13;21)); with + it is additional, with − it is missing; numerical changes add or remove
//   normal chromosomes. The sex chromosomes written are the normal ones present: an abnormal X is
//   written as such (46,X,i(X)(q10)) and replaces none of them, and −X, −Y of an acquired loss are
//   already reflected (45,X,−Y), unless the sex chromosomes are marked constitutional (46,XXYc,−X)
//   (ISCN 2009 pp. 56–57). The reference is the normal diploid of the constitutional sex (male
//   when a Y is written anywhere), times the ploidy.
// - Copy number. The coverage of each reference chromosome by the segments of all the cell's
//   chromosomes, a step function on the representative positions; its runs that depart from the
//   reference are the losses and gains, each with the bands at its ends as written, and three
//   extents: the representative one, the inner one (certainly within it: from the far side of its
//   first breakpoint's band to the near side of its last's) and the outer one (possibly within
//   it). A band's copy number is the set of values of the runs whose outer extents meet it: one
//   value when no breakpoint can lie in it.
// - Loss, gain and fusion (CytoGPS's LGF model, Abrams et al. 2015, 2019): for each abnormality
//   the bands where its chromosomes hold fewer copies than the chromosomes it replaced (loss),
//   more (gain), and the bands of its junctions (fusion), over the outer extents; the cell line's
//   model is their union. It is a model of events, not of net copy number: in
//   46,XX,+2,der(2;13)(q10;q10) the band 2p25.3 is both lost (by the derivative) and gained (by +2).
//
// Complexity: linear in the size of the description and the bands of the chromosomes involved.

#pragma once

#include <genomics/karyotype/bands.hpp>
#include <genomics/karyotype/format.hpp>
#include <genomics/karyotype/iscn.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace genomics::karyotype {

struct Point {
  enum class Kind : std::uint8_t { pter, qter, cen, band, unknown };
  Kind kind = Kind::band;
  int chromosome = -1;
  std::uint64_t at = 0, lo = 0, hi = 0;  // representative position; the interval it may lie in
  std::string label;                     // as written: q21, p36.33, q2?, pter, p10
  bool written_position = false;         // a nucleotide position of a sequencing result
  std::uint32_t id = 0;
  bool exact() const { return lo == hi; }
};

struct Segment {
  enum class Kind : std::uint8_t { reference, unknown, hsr, bases };
  Kind kind = Kind::reference;
  int chromosome = -1;
  Point from, to;
  std::string bases;
  bool reversed() const { return kind == Kind::reference && to.at < from.at; }
  std::uint64_t start() const { return std::min(from.at, to.at); }
  std::uint64_t end() const { return std::max(from.at, to.at); }
};

struct Chromosome {
  std::string key;     // unique in its cell line
  std::string label;   // "9", "der(22)t(9;22)(q34;q11.2)", "mar1"
  int centromere_of = -1;
  bool abnormal = false, ring = false, marker = false, counted = true;
  int source = -1;     // its abnormality in LineInterpretation::abnormalities; −1 for the normal complement
  std::vector<Segment> segments;
  std::string open, close;  // ":" where it begins or ends at a break that joined nothing
};

// A junction (::) of a chromosome of the cell: the end of one segment joined to the start of the
// next, or the closing of a ring.
struct Fusion {
  std::size_t chromosome = 0;
  Point a, b;
  bool ring = false;
};

struct BandCopies {
  std::size_t band = 0;
  int min = 0, max = 0, normal = 0;
};

// A run of the copy number that departs from the reference.
struct Change {
  int chromosome = -1;
  int copies = 0, normal = 0;
  std::uint64_t start = 0, end = 0;
  std::uint64_t inner_start = 0, inner_end = 0;  // empty when inner_start ≥ inner_end
  std::uint64_t outer_start = 0, outer_end = 0;
  std::string from_band, to_band;               // the bands at its ends, as written (pter and qter as the terminal bands)
  BandRange bands;                              // the ideogram's bands over the outer extent
  bool loss() const { return copies < normal; }
};

struct Lgf {
  std::vector<std::size_t> loss, gain, fusion;  // bands of the ideogram
  std::vector<std::string> fusion_labels;       // the bands of the junctions as written (2q10, 9q34)
};

struct EventResult {
  std::size_t abnormality = 0;
  std::vector<std::size_t> products;  // its chromosomes in the cell
  Lgf lgf;
};

struct LineInterpretation {
  std::size_t line = 0;
  int ploidy = 2;
  bool male = false, sex_known = true;
  std::array<int, chromosome_count> normal{};       // the reference's copies
  std::vector<Abnormality> abnormalities;      // with idem, sl and sdl expanded and omitted breakpoints filled in
  std::vector<Chromosome> chromosomes;
  std::vector<Fusion> fusions;
  std::vector<BandCopies> bands;               // every band of the ideogram, in its order
  std::vector<Change> changes;
  std::vector<EventResult> events;
  Lgf lgf;
  int counted = 0;                             // chromosomes counted (double minutes are not)
  int unknown = 0;                             // chromosomes or pieces of unknown origin
  bool balanced = true;                        // no change of copy number, unknown material aside
  std::vector<Diagnostic> diagnostics;
};

// A microarray's or a sequencing result's reading: the changes an arr (or seq region) states, and
// for seq the chromosomes its abnormalities make (on autosomes and the sex chromosomes as written,
// the sex being unknown).
struct TechniqueInterpretation {
  std::size_t technique = 0;
  std::string assembly;               // GRCh38 or GRCh37, from the build written
  std::vector<Change> changes;        // of the arr entries, at their coordinates
  bool has_line = false;
  LineInterpretation line;            // of the seq result's abnormalities
  std::vector<Diagnostic> diagnostics;
};

// The assembly of a build as ISCN writes it: GRCh38, hg38 → GRCh38; GRCh37, hg19, GRCh37/hg19 → GRCh37.
inline std::string assembly_of(std::string_view build) {
  if (build.find("38") != std::string_view::npos) return "GRCh38";
  if (build.find("37") != std::string_view::npos || build.find("19") != std::string_view::npos) return "GRCh37";
  if (build.find("36") != std::string_view::npos || build.find("18") != std::string_view::npos) return "NCBI36";
  return std::string(build);
}

struct Interpretation {
  std::vector<LineInterpretation> lines;
  std::vector<TechniqueInterpretation> techniques;
  std::vector<Diagnostic> diagnostics;
  bool ok() const {
    for (const auto& d : diagnostics)
      if (d.level == Diagnostic::Level::error) return false;
    for (const auto& l : lines)
      for (const auto& d : l.diagnostics)
        if (d.level == Diagnostic::Level::error) return false;
    return true;
  }
};

namespace detail {

using Word = std::vector<Segment>;

inline Word reversed(const Word& w) {
  Word out(w.rbegin(), w.rend());
  for (auto& s : out) std::swap(s.from, s.to);
  return out;
}
inline Word operator+(Word a, const Word& b) {
  a.insert(a.end(), b.begin(), b.end());
  return a;
}

class Interpreter {
 public:
  Interpreter(const Ideogram& g, std::vector<Diagnostic>& diag) : g_(g), diag_(&diag) {}
  void set_diagnostics(std::vector<Diagnostic>& d) { diag_ = &d; }

  void error(const Span& s, std::string m) { diag_->push_back({Diagnostic::Level::error, s.begin, s.end, std::move(m)}); }
  void warning(const Span& s, std::string m) { diag_->push_back({Diagnostic::Level::warning, s.begin, s.end, std::move(m)}); }
  void note(const Span& s, std::string m) { diag_->push_back({Diagnostic::Level::note, s.begin, s.end, std::move(m)}); }

  // ---- points --------------------------------------------------------------------------------------
  Point exact(int c, Point::Kind k, std::uint64_t at, std::string label) {
    Point p;
    p.kind = k;
    p.chromosome = c;
    p.at = p.lo = p.hi = at;
    p.label = std::move(label);
    p.id = ++ids_;
    return p;
  }
  Point pter(int c) { return exact(c, Point::Kind::pter, 0, "pter"); }
  Point qter(int c) { return exact(c, Point::Kind::qter, g_.length[c], "qter"); }
  Point cen(int c, char arm) { return exact(c, Point::Kind::cen, g_.centromere[c], arm == 'p' ? "p10" : "q10"); }
  std::uint64_t centromere(int c) const { return g_.centromere[c]; }

  Point resolve(const Breakpoint& b, int c, const Span& where) {
    if (b.chromosome >= 0) c = b.chromosome;
    Point p;
    p.chromosome = c;
    p.label = b.label();
    p.id = ++ids_;
    if (c < 0 || !g_.has(c) || b.kind == Breakpoint::Kind::unknown) {
      p.kind = Point::Kind::unknown;
      return p;
    }
    switch (b.kind) {
      case Breakpoint::Kind::pter: return pter(c);
      case Breakpoint::Kind::qter: return qter(c);
      case Breakpoint::Kind::cen: return cen(c, 'q');
      default: break;
    }
    if (b.centromeric()) return cen(c, b.arm);
    p.kind = Point::Kind::band;
    std::optional<BandRange> r = g_.find(c, b.arm, b.digits);
    if (!r) {
      error(where, "there is no band " + chromosome_name(c) + band_label(b.arm, b.digits) + " in the ideogram (" + g_.assembly + ", " +
                       std::to_string(g_.level) + " bands)");
      r = g_.arm(c, b.arm);
      if (!r) { p.kind = Point::Kind::unknown; return p; }
    }
    std::uint64_t lo = r->start, hi = r->end;
    if (b.interval()) {
      const auto r2 = g_.find(c, b.to_arm ? b.to_arm : b.arm, b.to_digits);
      if (!r2) error(where, "there is no band " + chromosome_name(c) + band_label(b.to_arm ? b.to_arm : b.arm, b.to_digits) + " in the ideogram");
      else { lo = std::min(lo, r2->start); hi = std::max(hi, r2->end); }
    }
    p.lo = lo;
    p.hi = hi;
    p.at = lo + (hi - lo) / 2;
    return p;
  }

  // A breakpoint given by a nucleotide position (1-based) at the lower or the upper end of the
  // segment it bounds.
  void position(Point& p, const Breakpoint& b, bool lower_end) {
    if (!b.has_position) return;
    const auto boundary = [&](std::uint64_t v) { return lower_end ? (v > 0 ? v - 1 : 0) : v; };
    std::uint64_t x = boundary(b.position), y = b.position_to ? boundary(b.position_to) : x;
    if (y < x) std::swap(x, y);
    p.lo = x;
    p.hi = y;
    p.at = x + (y - x) / 2;
    p.written_position = true;
    if (p.kind == Point::Kind::band || p.kind == Point::Kind::unknown) p.kind = Point::Kind::band;
  }

  // ---- pieces --------------------------------------------------------------------------------------
  Segment ref(int c, const Point& a, const Point& b) {
    Segment s;
    s.kind = Segment::Kind::reference;
    s.chromosome = c;
    s.from = a;
    s.to = b;
    return s;
  }
  Segment unknown_material() {
    Segment s;
    s.kind = Segment::Kind::unknown;
    return s;
  }
  Segment hsr_material() {
    Segment s;
    s.kind = Segment::Kind::hsr;
    return s;
  }
  Word normal(int c) { return {ref(c, pter(c), qter(c))}; }
  static bool centromeric(const Point& x) { return x.kind == Point::Kind::cen; }
  char arm_of(const Point& x) const {
    if (centromeric(x)) return x.label == "p10" ? 'p' : 'q';
    return x.at < centromere(x.chromosome) ? 'p' : 'q';
  }
  // The centric piece at x, ending at x; at the centromere, the arm named.
  Word prox(int c, const Point& x) {
    if (c < 0 || x.kind == Point::Kind::unknown) return {unknown_material()};
    if (centromeric(x)) return arm_of(x) == 'p' ? Word{ref(c, pter(c), cen(c, 'p'))} : Word{ref(c, qter(c), cen(c, 'q'))};
    return arm_of(x) == 'q' ? Word{ref(c, pter(c), x)} : Word{ref(c, qter(c), x)};
  }
  // The acentric piece at x, starting from x; at the centromere, the arm named.
  Word dist(int c, const Point& x) {
    if (c < 0 || x.kind == Point::Kind::unknown) return {unknown_material()};
    if (centromeric(x)) return arm_of(x) == 'p' ? Word{ref(c, cen(c, 'p'), pter(c))} : Word{ref(c, cen(c, 'q'), qter(c))};
    return arm_of(x) == 'q' ? Word{ref(c, x, qter(c))} : Word{ref(c, x, pter(c))};
  }
  // At the centromere, the arm not named, starting from the centromere.
  Word other_arm(int c, const Point& x) {
    return arm_of(x) == 'p' ? Word{ref(c, cen(c, 'q'), qter(c))} : Word{ref(c, cen(c, 'p'), pter(c))};
  }

  // The word read with its centric segment from pter toward qter; else from a pter; else as made.
  Word orient(const Word& w, int named) const {
    // A word that begins or ends at the pter of its own chromosome is read from it (ISCN 2009's
    // inv(3)(pter→p13::q21→p13::q21→qter), its centric segment inverted).
    if (named >= 0 && !w.empty()) {
      const auto own_pter = [&](const Point& p, const Segment& s) { return s.kind == Segment::Kind::reference && s.chromosome == named && p.kind == Point::Kind::pter; };
      const bool starts = own_pter(w.front().from, w.front()), ends = own_pter(w.back().to, w.back());
      if (starts) return w;
      if (ends) return reversed(w);
    }
    if (named >= 0) {
      const std::uint64_t c = centromere(named);
      for (const auto& s : w)
        if (s.kind == Segment::Kind::reference && s.chromosome == named && s.start() < c && c < s.end())
          return s.reversed() ? reversed(w) : w;
    }
    for (const auto& s : w) {
      if (s.kind != Segment::Kind::reference) continue;
      if (s.chromosome != named && named >= 0) continue;
      if (s.start() < centromere(s.chromosome) && centromere(s.chromosome) < s.end()) return s.reversed() ? reversed(w) : w;
    }
    if (!w.empty() && w.back().kind == Segment::Kind::reference && w.back().to.kind == Point::Kind::pter &&
        !(w.front().kind == Segment::Kind::reference && w.front().from.kind == Point::Kind::pter))
      return reversed(w);
    return w;
  }

  // ---- operations on words (components of derivatives) -------------------------------------------
  // The place of x in w: the index of the reference segment holding it, or −1.
  long locate(const Word& w, const Point& x) const {
    for (std::size_t i = 0; i < w.size(); i++) {
      const Segment& s = w[i];
      if (s.kind == Segment::Kind::reference && s.chromosome == x.chromosome && s.start() <= x.at && x.at <= s.end()) return static_cast<long>(i);
    }
    return -1;
  }
  // Cuts w at x; returns the index k of the first segment after the cut.
  std::size_t cut(Word& w, std::size_t i, const Point& x) {
    Segment& s = w[i];
    if (x.at == s.from.at) return i;
    if (x.at == s.to.at) return i + 1;
    Segment right = s;
    s.to = x;
    right.from = x;
    w.insert(w.begin() + static_cast<std::ptrdiff_t>(i) + 1, right);
    return i + 1;
  }
  // Whether the centromere of `named` lies in w[0, k) (true) or in w[k, …) (false); `x` decides
  // a cut at the centromere itself (the arm it names is kept).
  bool centromere_left(const Word& w, std::size_t k, int named, const Point& x) const {
    if (named < 0) return true;
    const std::uint64_t c = centromere(named);
    for (std::size_t i = 0; i < w.size(); i++) {
      const Segment& s = w[i];
      if (s.kind != Segment::Kind::reference || s.chromosome != named) continue;
      if (s.start() < c && c < s.end()) return i < k;
    }
    if (centromeric(x) && k > 0 && k < w.size()) {
      const Segment& left = w[k - 1];
      const char arm = arm_of(x);
      const bool left_on_arm = left.kind == Segment::Kind::reference && (arm == 'p' ? left.start() < c : left.end() > c);
      return left_on_arm;
    }
    return true;
  }
  Word middle(Word& w, const Point& x, const Point& y, std::size_t& a, std::size_t& b, const Span& where, bool& ok) {
    ok = false;
    const long ix = locate(w, x);
    if (ix < 0) { warning(where, "the breakpoint " + chromosome_name(x.chromosome) + x.label + " is not on this derivative"); return {}; }
    a = cut(w, static_cast<std::size_t>(ix), x);
    const long iy = locate(w, y);
    if (iy < 0) { warning(where, "the breakpoint " + chromosome_name(y.chromosome) + y.label + " is not on this derivative"); return {}; }
    const std::size_t before = w.size();
    b = cut(w, static_cast<std::size_t>(iy), y);
    // A split of a segment before the first cut moves that cut one place on.
    if (w.size() > before && static_cast<std::size_t>(iy) + 1 <= a) a++;
    if (b < a) std::swap(a, b);
    ok = true;
    return Word(w.begin() + static_cast<std::ptrdiff_t>(a), w.begin() + static_cast<std::ptrdiff_t>(b));
  }

  // ---- the point table of one abnormality -----------------------------------------------------------
  struct Table {
    std::map<const Breakpoint*, Point> points;
  };

  // Resolves every breakpoint of an abnormality (head, components, detailed form): a label that
  // recurs in another group of the same abnormality is the same point; breakpoints of one band are
  // spread through it in the order written.
  void resolve_all(const Abnormality& a, Table& t) {
    std::vector<std::pair<const Breakpoint*, std::pair<int, int>>> order;  // breakpoint, (chromosome, group)
    std::map<const Breakpoint*, bool> homologue;
    int group = 0;
    const auto add_groups = [&](const Abnormality& x) {
      for (std::size_t g = 0; g < x.groups.size(); g++) {
        const int c = g < x.chromosomes.size() ? (x.chromosomes[g].unknown ? -1 : x.chromosomes[g].index)
                                               : (x.chromosomes.empty() ? -1 : x.chromosomes[0].index);
        const bool h = g < x.chromosomes.size() && x.chromosomes[g].homologue;
        for (const auto& b : x.groups[g]) { order.push_back({&b, {c, group}}); homologue[&b] = h; }
        group++;
      }
    };
    add_groups(a);
    for (const auto& c : a.components) add_groups(c);
    // The first breakpoint of each (chromosome, label) and its group; a later one in another group
    // repeats it.
    std::map<std::pair<int, std::string>, std::pair<const Breakpoint*, int>> first;
    std::vector<std::pair<const Breakpoint*, const Breakpoint*>> repeats;
    std::vector<const Breakpoint*> fresh;
    for (const auto& [b, cg] : order) {
      const int c = b->chromosome >= 0 ? b->chromosome : cg.first;
      // 7 and 7* are two homologues: the same band on each is two breaks.
      const auto key = std::make_pair(c * 2 + (homologue[b] ? 1 : 0), b->label());
      const auto it = first.find(key);
      if (it != first.end() && it->second.second != cg.second && b->kind == Breakpoint::Kind::band) {
        repeats.push_back({b, it->second.first});
        continue;
      }
      t.points[b] = resolve(*b, c, b->span);
      fresh.push_back(b);
      if (it == first.end()) first[key] = {b, cg.second};
    }
    // Spread the band points that share an interval, in the order written.
    std::map<std::tuple<int, std::uint64_t, std::uint64_t>, std::vector<const Breakpoint*>> share;
    for (const Breakpoint* b : fresh) {
      const Point& p = t.points[b];
      if (p.kind == Point::Kind::band && !p.exact()) share[{p.chromosome, p.lo, p.hi}].push_back(b);
    }
    for (auto& [k, list] : share) {
      if (list.size() < 2) continue;
      const std::uint64_t lo = std::get<1>(k), hi = std::get<2>(k);
      for (std::size_t i = 0; i < list.size(); i++) t.points[list[i]].at = lo + (hi - lo) * (i + 1) / (list.size() + 1);
    }
    for (const auto& [b, src] : repeats) t.points[b] = t.points[src];
  }

  // Joins consecutive segments that continue one another through the same point in the same
  // direction (pter→x::x→y is pter→y).
  static Word merged(const Word& w) {
    Word out;
    for (const auto& s : w) {
      if (!out.empty()) {
        Segment& p = out.back();
        if (p.kind == Segment::Kind::reference && s.kind == Segment::Kind::reference && p.chromosome == s.chromosome && p.to.id == s.from.id &&
            p.to.at == s.from.at && (p.from.at <= p.to.at) == (s.from.at <= s.to.at) && p.from.at != p.to.at && s.from.at != s.to.at) {
          p.to = s.to;
          continue;
        }
      }
      out.push_back(s);
    }
    return out;
  }

  // ---- the products of one abnormality ----------------------------------------------------------------
  struct Products {
    std::vector<Chromosome> chromosomes;
    std::vector<int> replaced;  // normal chromosomes it replaces, with multiplicity
  };

  Chromosome make(Word w, int named, std::string label, bool ring = false) {
    Chromosome c;
    c.segments = orient(w, named);
    c.centromere_of = named;
    c.label = std::move(label);
    c.abnormal = true;
    c.ring = ring;
    return c;
  }

  static std::string product_label(const Abnormality& a, std::size_t k) {
    Abnormality bare = a;
    bare.sign = 0;
    bare.multiplier = 1;
    bare.constitutional = false;
    bare.inheritance.clear();
    bare.alternatives.clear();
    FormatOptions o;
    o.order = false;
    const std::string text = format(bare, o);
    if ((a.type == Type::t || a.type == Type::ins) && a.chromosomes.size() > 1 && k < a.chromosomes.size())
      return "der(" + a.chromosomes[k].name() + ")" + text;
    return text;
  }

  Word words_of_long(const LongChromosome& lc, int default_chromosome, Table& t, std::map<std::pair<int, std::string>, Point>& seen,
                     const Span& where, bool& ring, std::string& open, std::string& close) {
    Word w;
    ring = lc.open == "::" && lc.close == "::";
    open = lc.open == ":" ? ":" : "";
    close = lc.close == ":" ? ":" : "";
    (void)t;
    const auto point = [&](const Breakpoint& b, bool lower) {
      const int c = b.chromosome >= 0 ? b.chromosome : default_chromosome;
      const auto key = std::make_pair(c, b.label());
      Point p;
      const auto it = seen.find(key);
      if (it != seen.end() && b.kind == Breakpoint::Kind::band && !b.has_position) p = it->second;
      else {
        p = resolve(b, c, where);
        if (b.kind == Breakpoint::Kind::band && !b.has_position) seen[key] = p;
      }
      position(p, b, lower);
      return p;
    };
    for (const auto& part : lc.parts) {
      switch (part.kind) {
        case LongPart::Kind::segment: {
          const int c = part.from.chromosome >= 0 ? part.from.chromosome : default_chromosome;
          // Which end is lower decides how a written position becomes a boundary.
          Point a = point(part.from, true), b = point(part.to, false);
          if (a.kind != Point::Kind::unknown && b.kind != Point::Kind::unknown && b.at < a.at) {
            a = point(part.from, false);
            b = point(part.to, true);
          }
          if (c < 0 || a.kind == Point::Kind::unknown || b.kind == Point::Kind::unknown) { w.push_back(unknown_material()); break; }
          if (part.to.chromosome >= 0 && part.to.chromosome != c) error(where, "a segment runs from one chromosome to another");
          w.push_back(ref(c, a, b));
          break;
        }
        case LongPart::Kind::unknown: w.push_back(unknown_material()); break;
        case LongPart::Kind::bases: {
          Segment s;
          s.kind = Segment::Kind::bases;
          s.bases = part.bases;
          w.push_back(s);
          break;
        }
        case LongPart::Kind::region: {
          const int c = part.from.chromosome >= 0 ? part.from.chromosome : default_chromosome;
          Point a = resolve(part.from, c, where), b = resolve(part.to, c, where);
          position(a, part.from, true);
          position(b, part.to, false);
          for (int k = 0; k < std::max(1, part.multiplier); k++) w.push_back(ref(c, a, b));
          break;
        }
      }
    }
    return w;
  }

  // The word of a derivative after one of its abnormalities.
  void apply(Word& w, int named, const Abnormality& comp, Table& t, bool& ring, const std::vector<int>& partners = {}) {
    const auto P = [&](std::size_t g, std::size_t i) -> const Point& { return t.points.at(&comp.groups[g][i]); };
    const auto has = [&](std::size_t g, std::size_t n) { return comp.has_groups && g < comp.groups.size() && comp.groups[g].size() >= n; };
    switch (comp.type) {
      case Type::t: {
        if (!comp.has_groups) { error(comp.span, "the breakpoints of " + comp.text + " are given nowhere"); return; }
        const std::size_t n = comp.chromosomes.size();
        for (std::size_t k = 0; k < n; k++) {
          if (comp.chromosomes[k].unknown || !has(k, 1)) continue;
          const Point& x = P(k, 0);
          const long i = locate(w, x);
          if (i < 0) continue;
          const std::size_t prev = (k + n - 1) % n;
          const int pc = comp.chromosomes[prev].unknown ? -1 : comp.chromosomes[prev].index;
          // A derivative of two chromosomes, der(5;7), holds both centromeres: the second's piece
          // is its centric one.
          bool centric = false;
          if (pc >= 0 && std::find(partners.begin(), partners.end(), pc) != partners.end()) {
            centric = true;
            for (const auto& s : w)
              if (s.kind == Segment::Kind::reference && s.chromosome == pc) centric = false;
          }
          const Word piece = !has(prev, 1) ? Word{unknown_material()} : centric ? reversed(prox(pc, P(prev, 0))) : dist(pc, P(prev, 0));
          const std::size_t at = cut(w, static_cast<std::size_t>(i), x);
          if (centromere_left(w, at, named, x)) {
            w.erase(w.begin() + static_cast<std::ptrdiff_t>(at), w.end());
            w = w + piece;
          } else {
            w.erase(w.begin(), w.begin() + static_cast<std::ptrdiff_t>(at));
            w = reversed(piece) + w;
          }
          return;
        }
        // The derivative of an unknown chromosome receives the piece of the other.
        for (std::size_t k = 0; k < n; k++) {
          if (!comp.chromosomes[k].unknown) continue;
          const std::size_t prev = (k + n - 1) % n;
          if (comp.chromosomes[prev].unknown || !has(prev, 1)) continue;
          w = w + dist(comp.chromosomes[prev].index, P(prev, 0));
          return;
        }
        warning(comp.span, "none of the breakpoints of " + comp.text + " is on this derivative");
        return;
      }
      case Type::del: {
        if (has(0, 2)) {
          std::size_t a = 0, b = 0;
          bool ok = false;
          middle(w, P(0, 0), P(0, 1), a, b, comp.span, ok);
          if (ok) w.erase(w.begin() + static_cast<std::ptrdiff_t>(a), w.begin() + static_cast<std::ptrdiff_t>(b));
        } else if (has(0, 1)) {
          const Point& x = P(0, 0);
          const long i = locate(w, x);
          if (i < 0) { warning(comp.span, "the breakpoint " + x.label + " is not on this derivative"); return; }
          const std::size_t at = cut(w, static_cast<std::size_t>(i), x);
          if (centromere_left(w, at, named, x)) w.erase(w.begin() + static_cast<std::ptrdiff_t>(at), w.end());
          else w.erase(w.begin(), w.begin() + static_cast<std::ptrdiff_t>(at));
        }
        return;
      }
      case Type::dup: case Type::trp: case Type::qdp: case Type::inv: {
        if (!has(0, 2)) return;
        std::size_t a = 0, b = 0;
        bool ok = false;
        const Word m = middle(w, P(0, 0), P(0, 1), a, b, comp.span, ok);
        if (!ok) return;
        const bool inverted = P(0, 0).at > P(0, 1).at || comp.has_prefix("inv");
        if (comp.type == Type::inv) {
          const Word r = reversed(m);
          std::copy(r.begin(), r.end(), w.begin() + static_cast<std::ptrdiff_t>(a));
          return;
        }
        const int extra = comp.type == Type::dup ? 1 : comp.type == Type::trp ? 2 : 3;
        Word add;
        for (int k = 0; k < extra; k++) add = add + ((inverted && k % 2 == 0) ? reversed(m) : m);
        w.insert(w.begin() + static_cast<std::ptrdiff_t>(b), add.begin(), add.end());
        return;
      }
      case Type::ins: {
        const bool two = comp.chromosomes.size() == 2;
        if (two && has(0, 1)) {
          const Point& x = P(0, 0);
          const long i = locate(w, x);
          if (i >= 0) {
            Word seg;
            if (comp.chromosomes[1].unknown || !has(1, 2)) seg = {unknown_material()};
            else seg = {ref(comp.chromosomes[1].index, P(1, 0), P(1, 1))};
            if (w[static_cast<std::size_t>(i)].reversed()) seg = reversed(seg);
            const std::size_t at = cut(w, static_cast<std::size_t>(i), x);
            w.insert(w.begin() + static_cast<std::ptrdiff_t>(at), seg.begin(), seg.end());
            return;
          }
          if (has(1, 2)) {
            std::size_t a = 0, b = 0;
            bool ok = false;
            middle(w, P(1, 0), P(1, 1), a, b, comp.span, ok);
            if (ok) { w.erase(w.begin() + static_cast<std::ptrdiff_t>(a), w.begin() + static_cast<std::ptrdiff_t>(b)); return; }
          }
        }
        if (two && comp.chromosomes[0].unknown) {
          // ins(?;7)(?;q11.2q32) on a derivative of unknown origin: into its unknown material.
          for (std::size_t i = 0; i < w.size(); i++)
            if (w[i].kind == Segment::Kind::unknown) {
              const Word seg = (comp.chromosomes[1].unknown || !has(1, 2)) ? Word{unknown_material()} : Word{ref(comp.chromosomes[1].index, P(1, 0), P(1, 1))};
              w.insert(w.begin() + static_cast<std::ptrdiff_t>(i) + 1, seg.begin(), seg.end());
              return;
            }
        }
        warning(comp.span, comp.text + " is not interpreted on this derivative");
        return;
      }
      case Type::add: case Type::hsr: {
        if (comp.type == Type::hsr && (!has(0, 1) || comp.chromosomes.empty() || comp.chromosomes[0].unknown)) {
          w.push_back(hsr_material());  // hsr(?): amplified material of unknown origin
          return;
        }
        if (!has(0, 1)) return;
        const Point& x = P(0, 0);
        const long i = locate(w, x);
        if (i < 0) { warning(comp.span, "the breakpoint " + x.label + " is not on this derivative"); return; }
        const std::size_t at = cut(w, static_cast<std::size_t>(i), x);
        if (comp.type == Type::hsr) {
          w.insert(w.begin() + static_cast<std::ptrdiff_t>(at), hsr_material());
        } else if (centromere_left(w, at, named, x)) {
          w.erase(w.begin() + static_cast<std::ptrdiff_t>(at), w.end());
          w.push_back(unknown_material());
        } else {
          w.erase(w.begin(), w.begin() + static_cast<std::ptrdiff_t>(at));
          w.insert(w.begin(), unknown_material());
        }
        return;
      }
      case Type::r: {
        Word rw;
        for (std::size_t g = 0; g < comp.groups.size(); g++) {
          const int c = g < comp.chromosomes.size() && !comp.chromosomes[g].unknown ? comp.chromosomes[g].index : -1;
          if (c < 0 || !has(g, 2)) rw.push_back(unknown_material());
          else rw.push_back(ref(c, P(g, 0), P(g, 1)));
        }
        if (!rw.empty()) { w = rw; ring = true; }
        return;
      }
      default:
        warning(comp.span, comp.text + " is not interpreted inside a derivative chromosome");
        return;
    }
  }

  Products realise(const Abnormality& a) {
    Products out;
    Table t;
    resolve_all(a, t);
    const auto P = [&](std::size_t g, std::size_t i) -> const Point& { return t.points.at(&a.groups[g][i]); };
    const auto has = [&](std::size_t g, std::size_t n) { return a.has_groups && g < a.groups.size() && a.groups[g].size() >= n; };
    const auto chrom = [&](std::size_t k) { return k < a.chromosomes.size() && !a.chromosomes[k].unknown ? a.chromosomes[k].index : -1; };
    const auto replace_all = [&]() {
      for (const auto& c : a.chromosomes)
        if (!c.unknown && c.index >= 0) out.replaced.push_back(c.index);
    };
    const int c0 = chrom(0);
    const auto single = [&](Word w, bool ring = false) {
      out.chromosomes.push_back(make(std::move(w), c0, product_label(a, 0), ring));
      if (c0 >= 0) out.replaced.push_back(c0);
    };

    // The detailed system: one word per chromosome written.
    if (!a.long_form.empty()) {
      std::map<std::pair<int, std::string>, Point> seen;
      for (std::size_t k = 0; k < a.long_form.size(); k++) {
        bool ring = false;
        std::string open, close;
        const int named = a.long_form.size() > 1 ? chrom(k) : c0;
        Word w = words_of_long(a.long_form[k], a.chromosomes.size() == 1 ? c0 : -1, t, seen, a.span, ring, open, close);
        Chromosome ch = make(std::move(w), named, product_label(a, k), ring);
        ch.open = open;
        ch.close = close;
        out.chromosomes.push_back(std::move(ch));
      }
      if (a.type == Type::der && c0 >= 0) out.replaced.push_back(c0);
      else replace_all();
      return out;
    }

    switch (a.type) {
      case Type::t: {
        const std::size_t n = a.chromosomes.size();
        if (!a.has_groups || a.groups.size() != n) { error(a.span, "the breakpoints of " + a.text + " are given nowhere"); return out; }
        const bool pairs = std::all_of(a.groups.begin(), a.groups.end(), [](const auto& g) { return g.size() == 2; });
        if (pairs) {
          // Segments exchanged: der(c_k) has its segment replaced by that of c_{k−1}, in the
          // orientation written.
          for (std::size_t k = 0; k < n; k++) {
            const std::size_t prev = (k + n - 1) % n;
            const int c = chrom(k), pc = chrom(prev);
            if (c < 0) { out.chromosomes.push_back(make({unknown_material()}, -1, product_label(a, k))); continue; }
            const Point& x1 = P(k, 0);
            const Point& x2 = P(k, 1);
            const Point& lo = x1.at <= x2.at ? x1 : x2;
            const Point& hi = x1.at <= x2.at ? x2 : x1;
            const Word seg = pc < 0 ? Word{unknown_material()} : Word{ref(pc, P(prev, 0), P(prev, 1))};
            out.chromosomes.push_back(make(Word{ref(c, pter(c), lo)} + seg + Word{ref(c, hi, qter(c))}, c, product_label(a, k)));
          }
        } else {
          const bool whole_arm = n == 2 && has(0, 1) && has(1, 1) && centromeric(P(0, 0)) && centromeric(P(1, 0));
          for (std::size_t k = 0; k < n; k++) {
            const std::size_t prev = (k + n - 1) % n;
            const int c = chrom(k), pc = chrom(prev);
            Word w;
            if (whole_arm && k == 1) {
              w = reversed(other_arm(chrom(0), P(0, 0))) + other_arm(chrom(1), P(1, 0));
            } else if (whole_arm && k == 0) {
              w = prox(c, P(0, 0)) + dist(chrom(1), P(1, 0));
            } else {
              w = prox(c, P(k, 0)) + dist(pc, P(prev, 0));
            }
            out.chromosomes.push_back(make(w, c, product_label(a, k)));
          }
        }
        replace_all();
        return out;
      }
      case Type::del: {
        if (c0 < 0) return out;
        if (has(0, 2)) {
          const Point& x = P(0, 0).at <= P(0, 1).at ? P(0, 0) : P(0, 1);
          const Point& y = P(0, 0).at <= P(0, 1).at ? P(0, 1) : P(0, 0);
          single(Word{ref(c0, pter(c0), x), ref(c0, y, qter(c0))});
        } else if (has(0, 1)) {
          const Point& x = P(0, 0);
          if (arm_of(x) == 'q') { single(Word{ref(c0, pter(c0), x)}); out.chromosomes.back().close = ":"; }
          else { single(Word{ref(c0, x, qter(c0))}); out.chromosomes.back().open = ":"; }
        }
        return out;
      }
      case Type::dup: case Type::trp: case Type::qdp: {
        if (c0 < 0 || !has(0, 2)) return out;
        const Point& w1 = P(0, 0);
        const Point& w2 = P(0, 1);
        const bool inverted = w2.at < w1.at;
        const Point& x = inverted ? w2 : w1;
        const Point& y = inverted ? w1 : w2;
        // pter→x, then the segment n + 1 times (alternately reversed when inverted), then y→qter;
        // merged, this is ISCN 2009's pter→q25::q22→qter for dup(1)(q22q25), pter→q32::q21→q32::
        // q21→qter for trp(1)(q21q32), and pter→q25::q25→q22::q25→qter for dup(1)(q25q22) (p. 69,
        // the first of its two descriptions).
        const int copies = a.type == Type::dup ? 2 : a.type == Type::trp ? 3 : 4;
        Word w{ref(c0, pter(c0), x)};
        for (int k = 0; k < copies; k++) w.push_back(inverted && k % 2 == 1 ? ref(c0, y, x) : ref(c0, x, y));
        w.push_back(ref(c0, y, qter(c0)));
        single(merged(w));
        return out;
      }
      case Type::inv: {
        if (c0 < 0 || !has(0, 2)) return out;
        const Point& x = P(0, 0).at <= P(0, 1).at ? P(0, 0) : P(0, 1);
        const Point& y = P(0, 0).at <= P(0, 1).at ? P(0, 1) : P(0, 0);
        single(Word{ref(c0, pter(c0), x), ref(c0, y, x), ref(c0, y, qter(c0))});
        return out;
      }
      case Type::ins: {
        if (a.chromosomes.size() == 2 && has(0, 1)) {
          const int b = chrom(1);
          const Point& x = P(0, 0);
          const Word seg = (b < 0 || !has(1, 2)) ? Word{unknown_material()} : Word{ref(b, P(1, 0), P(1, 1))};
          if (c0 >= 0) {
            out.chromosomes.push_back(make(Word{ref(c0, pter(c0), x)} + seg + Word{ref(c0, x, qter(c0))}, c0, product_label(a, 0)));
            out.replaced.push_back(c0);
          }
          if (b >= 0 && has(1, 2)) {
            const Point& y1 = P(1, 0).at <= P(1, 1).at ? P(1, 0) : P(1, 1);
            const Point& y2 = P(1, 0).at <= P(1, 1).at ? P(1, 1) : P(1, 0);
            out.chromosomes.push_back(make(Word{ref(b, pter(b), y1), ref(b, y2, qter(b))}, b, product_label(a, 1)));
            out.replaced.push_back(b);
          }
        } else if (a.chromosomes.size() == 1 && has(0, 3) && c0 >= 0) {
          const Point& x = P(0, 0);
          const Point& y1 = P(0, 1);
          const Point& y2 = P(0, 2);
          const Point& m = y1.at <= y2.at ? y1 : y2;
          const Point& M = y1.at <= y2.at ? y2 : y1;
          const Segment seg = ref(c0, y1, y2);
          if (x.at < m.at) single(Word{ref(c0, pter(c0), x), seg, ref(c0, x, m), ref(c0, M, qter(c0))});
          else if (x.at > M.at) single(Word{ref(c0, pter(c0), m), ref(c0, M, x), seg, ref(c0, x, qter(c0))});
          else error(a.span, "the point of insertion lies within the inserted segment");
        }
        return out;
      }
      case Type::add: {
        if (c0 < 0 || !has(0, 1)) return out;
        const Point& x = P(0, 0);
        if (arm_of(x) == 'q') single(Word{ref(c0, pter(c0), x), unknown_material()});
        else single(Word{unknown_material(), ref(c0, x, qter(c0))});
        return out;
      }
      case Type::hsr: {
        if (c0 < 0 || !has(0, 1)) return out;
        const Point& x = P(0, 0);
        single(Word{ref(c0, pter(c0), x), hsr_material(), ref(c0, x, qter(c0))});
        return out;
      }
      case Type::fra: case Type::neo: {
        if (c0 < 0) return out;
        single(normal(c0));
        return out;
      }
      case Type::i: case Type::idic: {
        if (c0 < 0 || !has(0, 1)) return out;
        const Point& x = P(0, 0);
        if (a.type == Type::i && centromeric(x)) {
          const Word armw = arm_of(x) == 'p' ? Word{ref(c0, pter(c0), cen(c0, 'p'))} : Word{ref(c0, qter(c0), cen(c0, 'q'))};
          single(armw + reversed(armw));
        } else {
          const Word p = prox(c0, x);
          single(p + reversed(p));
        }
        return out;
      }
      case Type::fis: {
        if (c0 < 0 || !has(0, 1)) return out;
        // Each product of a centric fission is read from its telomere to its broken centromere
        // (ISCN 2009 p. 69: fis(10)(pter→p10:), fis(10)(qter→q10:)).
        const Point& x = P(0, 0);
        single(arm_of(x) == 'p' ? Word{ref(c0, pter(c0), cen(c0, 'p'))} : Word{ref(c0, qter(c0), cen(c0, 'q'))});
        out.chromosomes.back().segments = arm_of(x) == 'p' ? Word{ref(c0, pter(c0), cen(c0, 'p'))} : Word{ref(c0, qter(c0), cen(c0, 'q'))};
        out.chromosomes.back().close = ":";
        return out;
      }
      case Type::r: {
        Word w;
        if (!a.has_groups) {
          // A ring without breakpoints loses at least its terminal bands (as CytoGPS reads it).
          if (c0 < 0) return out;
          const Band& tp = g_.bands[g_.pter_band(c0)];
          const Band& tq = g_.bands[g_.qter_band(c0)];
          Point x = exact(c0, Point::Kind::band, tp.start + tp.length() / 2, tp.name());
          x.lo = tp.start;
          x.hi = tp.end;
          Point y = exact(c0, Point::Kind::band, tq.start + tq.length() / 2, tq.name());
          y.lo = tq.start;
          y.hi = tq.end;
          w.push_back(ref(c0, x, y));
          note(a.span, a.text + " gives no breakpoints; they are taken in its terminal bands, " + chromosome_name(c0) + tp.name() + " and " +
                           chromosome_name(c0) + tq.name());
        } else {
          for (std::size_t g = 0; g < a.groups.size(); g++) {
            const int c = chrom(g < a.chromosomes.size() ? g : 0);
            if (c < 0 || !has(g, 2)) w.push_back(unknown_material());
            else w.push_back(ref(c, P(g, 0), P(g, 1)));
          }
        }
        out.chromosomes.push_back(make(w, c0, product_label(a, 0), true));
        replace_all();
        return out;
      }
      case Type::dic: case Type::rob: {
        if (a.chromosomes.size() < 2 || !has(0, 1) || !has(1, 1)) return out;
        out.chromosomes.push_back(make(prox(c0, P(0, 0)) + reversed(prox(chrom(1), P(1, 0))), c0, product_label(a, 0)));
        replace_all();
        return out;
      }
      case Type::trc: {
        if (a.chromosomes.size() != 3 || !has(0, 1) || !has(1, 2) || !has(2, 1)) return out;
        const int b = chrom(1);
        const Word mid = b < 0 ? Word{unknown_material()} : Word{ref(b, P(1, 0), P(1, 1))};
        out.chromosomes.push_back(make(prox(c0, P(0, 0)) + mid + reversed(prox(chrom(2), P(2, 0))), c0, product_label(a, 0)));
        replace_all();
        return out;
      }
      case Type::tas: {
        // A telomeric association joins chromosomes end to end without a proven fusion (ISCN 2020
        // §9.2.16): each stays a chromosome, counted, its content whole.
        // The sex chromosomes stay in the sex chromosomes written.
        for (std::size_t k = 0; k < a.chromosomes.size(); k++) {
          const int c = chrom(k);
          if (c < 0 || is_sex(c)) continue;
          out.chromosomes.push_back(make(normal(c), c, product_label(a, 0)));
          out.replaced.push_back(c);
        }
        return out;
      }
      case Type::der: case Type::ider: case Type::rec: {
        Word w;
        int named = c0;
        bool ring = false;
        if (a.type == Type::rec) {
          if (c0 < 0) return out;
          // The recombinants of a pericentric inversion (theory §4.6): rec(c)dup(cq)inv(c)(x y) is
          // qter→y::x→qter; dup(cp), pter→y::x→pter.
          const Abnormality* inv = nullptr;
          for (const auto& comp : a.components)
            if (comp.type == Type::inv && comp.groups.size() == 1 && comp.groups[0].size() == 2) inv = &comp;
          if (!inv) { warning(a.span, "a recombinant chromosome is interpreted from the inversion it comes from"); return out; }
          const Point& p1 = t.points.at(&inv->groups[0][0]);
          const Point& p2 = t.points.at(&inv->groups[0][1]);
          const Point& x = p1.at <= p2.at ? p1 : p2;
          const Point& y = p1.at <= p2.at ? p2 : p1;
          if (!(arm_of(x) == 'p' && arm_of(y) == 'q')) {
            warning(a.span, "the recombinants of a paracentric inversion are dicentric or acentric and are not interpreted");
            return out;
          }
          const bool dup_q = (a.rec_dup && a.rec_arm == 'q') || (!a.rec_dup && a.rec_arm == 'p');
          w = dup_q ? Word{ref(c0, qter(c0), y), ref(c0, x, qter(c0))} : Word{ref(c0, pter(c0), y), ref(c0, x, pter(c0))};
          single(w);
          return out;
        }
        if (a.has_prefix("der") && a.type != Type::der) {
          // der dic(1;20)(…)ins(…): the head is made first.
          Abnormality head = a;
          head.components.clear();
          head.prefixes.clear();
          Products hp = realise(head);
          if (hp.chromosomes.empty()) return out;
          w = hp.chromosomes[0].segments;
          out.replaced = hp.replaced;
        } else if (a.chromosomes.size() >= 2 && a.has_groups && a.groups.size() >= 2 && has(0, 1) && has(1, 1)) {
          w = prox(c0, P(0, 0)) + reversed(prox(chrom(1), P(1, 0)));  // der(1;3)(p10;q10), as a dicentric
          replace_all();
        } else if (a.chromosomes.size() >= 2 && c0 >= 0) {
          // der(5;7)t(5;7)(q22;p13)…: made by its abnormalities from the first, the others' centric
          // pieces brought in by them.
          w = normal(c0);
          replace_all();
        } else if (c0 >= 0) {
          w = normal(c0);
          out.replaced.push_back(c0);
        } else {
          w = {unknown_material()};
          named = -1;
        }
        std::vector<int> partners;
        for (std::size_t k = 1; k < a.chromosomes.size(); k++)
          if (!a.chromosomes[k].unknown) partners.push_back(a.chromosomes[k].index);
        for (const auto& comp : a.components) apply(w, named, comp, t, ring, partners);
        if (a.type == Type::ider) {
          if (!has(0, 1) || c0 < 0) return out;
          w = orient(w, c0);
          const char armx = arm_of(P(0, 0));
          const std::uint64_t c = centromere(c0);
          long i = -1;
          for (std::size_t k = 0; k < w.size(); k++)
            if (w[k].kind == Segment::Kind::reference && w[k].chromosome == c0 && w[k].start() <= c && c <= w[k].end()) { i = static_cast<long>(k); break; }
          if (i < 0) { warning(a.span, "the derivative has no centromere of " + chromosome_name(c0) + " to double an arm at"); return out; }
          const std::size_t at = cut(w, static_cast<std::size_t>(i), cen(c0, armx));
          Word half = armx == 'q' ? Word(w.begin() + static_cast<std::ptrdiff_t>(at), w.end()) : Word(w.begin(), w.begin() + static_cast<std::ptrdiff_t>(at));
          // Label the point of fusion by the arm doubled.
          w = armx == 'q' ? reversed(half) + half : half + reversed(half);
        }
        out.chromosomes.push_back(make(std::move(w), named, product_label(a, 0), ring));
        return out;
      }
      default:
        return out;
    }
  }

  const Ideogram& g_;
  std::vector<Diagnostic>* diag_;
  std::uint32_t ids_ = 0;
};

// The coverage of chromosome c by the segments of a set of chromosomes, as runs over the
// representative positions; each boundary keeps the widest uncertainty of the points there.
struct Run {
  std::uint64_t start = 0, end = 0;
  int copies = 0;
  Point first, last;  // the points at its ends
};

inline std::vector<Run> coverage(const std::vector<const Chromosome*>& cs, int c, const Ideogram& g, int base = 0) {
  struct Event { std::uint64_t at; int delta; Point p; };
  std::vector<Event> ev;
  for (const Chromosome* ch : cs)
    for (const Segment& s : ch->segments) {
      if (s.kind != Segment::Kind::reference || s.chromosome != c || s.start() == s.end()) continue;
      const Point& a = s.from.at <= s.to.at ? s.from : s.to;
      const Point& b = s.from.at <= s.to.at ? s.to : s.from;
      ev.push_back({a.at, +1, a});
      ev.push_back({b.at, -1, b});
    }
  Point pter_p;
  pter_p.kind = Point::Kind::pter;
  pter_p.chromosome = c;
  pter_p.label = "pter";
  Point qter_p = pter_p;
  qter_p.kind = Point::Kind::qter;
  qter_p.label = "qter";
  qter_p.at = qter_p.lo = qter_p.hi = g.length[c];
  ev.push_back({0, 0, pter_p});
  ev.push_back({g.length[c], 0, qter_p});
  std::sort(ev.begin(), ev.end(), [](const Event& x, const Event& y) { return x.at < y.at; });
  std::vector<Run> runs;
  int level = base;
  std::size_t i = 0;
  std::optional<Point> at_point;
  std::uint64_t pos = 0;
  Point last_point = pter_p;
  while (i < ev.size()) {
    const std::uint64_t x = ev[i].at;
    Point here = ev[i].p;
    int delta = 0;
    for (; i < ev.size() && ev[i].at == x; i++) {
      delta += ev[i].delta;
      const Point& q = ev[i].p;
      // The widest uncertainty, and a written label over pter or qter.
      if (q.hi - q.lo > here.hi - here.lo || (here.kind != Point::Kind::band && q.kind == Point::Kind::band && q.hi - q.lo >= here.hi - here.lo)) here = q;
    }
    if (x > pos) {
      if (!runs.empty() && runs.back().copies == level && runs.back().end == pos) {
        runs.back().end = x;
        runs.back().last = here;
      } else {
        runs.push_back({pos, x, level, last_point, here});
      }
    }
    level += delta;
    pos = x;
    last_point = here;
  }
  return runs;
}

inline std::string end_band(const Point& p, int c, const Ideogram& g, bool start) {
  if (p.kind == Point::Kind::pter || (start && p.at == 0)) return g.bands[g.pter_band(c)].name();
  if (p.kind == Point::Kind::qter || (!start && p.at == g.length[c])) return g.bands[g.qter_band(c)].name();
  if (p.kind == Point::Kind::cen) return g.bands[g.band_at(c, start ? p.at : (p.at > 0 ? p.at - 1 : 0))].name();
  if (p.written_position) return g.bands[g.band_at(c, p.at > 0 && !start ? p.at - 1 : p.at)].name();
  return p.label;
}

inline void mark(std::vector<char>& v, const Ideogram& g, int c, std::uint64_t lo, std::uint64_t hi) {
  if (!g.has(c)) return;
  if (hi <= lo) {
    v[g.band_at(c, std::min(lo, g.length[c] > 0 ? g.length[c] - 1 : 0))] = 1;
    return;
  }
  const BandRange r = g.bands_over(c, lo, hi);
  for (std::size_t b = r.first; b <= r.last; b++) v[b] = 1;
}

inline std::vector<std::size_t> marked(const std::vector<char>& v) {
  std::vector<std::size_t> out;
  for (std::size_t i = 0; i < v.size(); i++)
    if (v[i]) out.push_back(i);
  return out;
}

// The fusions of one chromosome: consecutive segments, and the closing of a ring.
inline void junctions(const Chromosome& ch, std::size_t index, std::vector<Fusion>& out) {
  const auto& s = ch.segments;
  for (std::size_t i = 0; i + 1 < s.size(); i++)
    if (s[i].kind == Segment::Kind::reference && s[i + 1].kind == Segment::Kind::reference) out.push_back({index, s[i].to, s[i + 1].from, false});
    else if (s[i].kind == Segment::Kind::reference && s[i + 1].kind != Segment::Kind::reference) {
      Point none;
      none.kind = Point::Kind::unknown;
      out.push_back({index, s[i].to, none, false});
    } else if (s[i].kind != Segment::Kind::reference && s[i + 1].kind == Segment::Kind::reference) {
      Point none;
      none.kind = Point::Kind::unknown;
      out.push_back({index, none, s[i + 1].from, false});
    }
  if (ch.ring && !s.empty() && s.back().kind == Segment::Kind::reference && s.front().kind == Segment::Kind::reference)
    out.push_back({index, s.back().to, s.front().from, true});
}

inline std::string fusion_label(const Point& p) { return chromosome_name(p.chromosome) + p.label; }

// Marks the bands of a junction's point: its interval of uncertainty, or the band holding it; a
// centromeric point is named (2q10) but holds no band.
inline void mark_point(std::vector<char>& v, std::vector<std::string>& labels, const Point& p, const Ideogram& g) {
  if (p.kind == Point::Kind::unknown || p.chromosome < 0) return;
  if (p.kind == Point::Kind::cen) {
    labels.push_back(fusion_label(p));
    return;
  }
  labels.push_back(chromosome_name(p.chromosome) + (p.kind == Point::Kind::band ? p.label : end_band(p, p.chromosome, g, p.kind == Point::Kind::pter)));
  if (p.exact()) {
    const std::uint64_t x = p.at;
    // An exact point between two bases belongs to the band holding the base after it, or the
    // last band at the end.
    mark(v, g, p.chromosome, x, x);
  } else {
    mark(v, g, p.chromosome, p.lo, p.hi);
  }
}

}  // namespace detail

// Interprets every cell line of a description on an ideogram.
inline Interpretation interpret(const Description& d, const Ideogram& g) {
  Interpretation out;
  detail::Interpreter in(g, out.diagnostics);

  // The abnormalities of each line, with idem, sl and sdl expanded.
  std::vector<std::vector<Abnormality>> expanded(d.lines.size());
  std::vector<const CellLine*> sex_of(d.lines.size(), nullptr);
  for (std::size_t i = 0; i < d.lines.size(); i++) {
    const CellLine& l = d.lines[i];
    std::vector<Abnormality> list;
    sex_of[i] = &l;
    if (l.reference != CellLine::Reference::none && i > 0) {
      std::size_t from = 0;
      if (l.reference == CellLine::Reference::sdl) from = l.sdl_number > 0 ? std::min<std::size_t>(static_cast<std::size_t>(l.sdl_number), i - 1) : i - 1;
      list = expanded[from];
      sex_of[i] = sex_of[from];
    }
    for (const auto& a : l.abnormalities) {
      if (a.sign < 0 && !a.numerical() && a.type != Type::marker) {
        // −dic(1;20): the line lacks an abnormal chromosome of the line it refers to, named by
        // its type and chromosomes.
        // The chromosomes it replaced stay missing; only the abnormal chromosome itself is gone.
        bool removed = false;
        for (auto it = list.begin(); it != list.end(); ++it) {
          if (it->type == a.type && it->sign >= 0 && !it->lost && it->chromosomes.size() == a.chromosomes.size() &&
              std::equal(it->chromosomes.begin(), it->chromosomes.end(), a.chromosomes.begin(),
                         [](const ChromosomeRef& x, const ChromosomeRef& y) { return x.index == y.index && x.unknown == y.unknown; })) {
            it->lost = true;
            removed = true;
            break;
          }
        }
        if (removed) continue;
      }
      list.push_back(a);
    }
    expanded[i] = std::move(list);
  }

  // Breakpoints omitted where the same rearrangement was written with them: +der(22)t(9;22).
  const auto fill = [&](Abnormality& a, std::size_t line) {
    const auto find = [&](const Abnormality& x) -> const Abnormality* {
      const auto same = [&](const Abnormality& y) {
        return y.type == x.type && y.has_groups && y.chromosomes.size() == x.chromosomes.size() &&
               std::equal(y.chromosomes.begin(), y.chromosomes.end(), x.chromosomes.begin(),
                          [](const ChromosomeRef& p, const ChromosomeRef& q) { return p.index == q.index && p.unknown == q.unknown; });
      };
      for (std::size_t k = 0; k < d.lines.size(); k++) {
        const std::size_t li = (line + k) % d.lines.size();
        for (const auto& y : expanded[li]) {
          if (same(y)) return &y;
          for (const auto& c : y.components)
            if (same(c)) return &c;
        }
      }
      return nullptr;
    };
    const auto fill_one = [&](Abnormality& x) {
      if (x.has_groups || !x.long_form.empty() || x.numerical() || x.type == Type::marker || x.type == Type::ring_marker || x.type == Type::dmin) return;
      if (x.type != Type::t && x.type != Type::ins && x.type != Type::dic && x.type != Type::inv) return;
      if (const Abnormality* src = find(x)) {
        x.groups = src->groups;
        x.has_groups = true;
      }
    };
    fill_one(a);
    for (auto& c : a.components) fill_one(c);
  };

  for (std::size_t li = 0; li < d.lines.size(); li++) {
    LineInterpretation L;
    L.line = li;
    in.set_diagnostics(L.diagnostics);
    const CellLine& line = d.lines[li];
    const CellLine& sexline = *sex_of[li];
    for (auto& a : expanded[li]) fill(a, li);
    L.abnormalities = expanded[li];

    // The ploidy, and the constitutional sex.
    int count = line.count_unknown ? 46 : line.count;
    if (line.count < 0 && sexline.count > 0) count = sexline.count;
    L.ploidy = line.ploidy > 0 ? line.ploidy : std::max(1, (count + 11) / 23);
    if (line.ploidy == 0 && line.reference != CellLine::Reference::none && sexline.ploidy > 0) L.ploidy = sexline.ploidy;
    const int k = L.ploidy;
    L.sex_known = !sexline.undisclosed;
    bool y_seen = sexline.sex.find('Y') != std::string::npos;
    for (const auto& a : L.abnormalities)
      for (const auto& c : a.chromosomes)
        if (c.index == chr_y) y_seen = true;
    L.male = y_seen;
    for (int c = 0; c < 22; c++) L.normal[c] = k;
    if (L.male) { L.normal[chr_x] = k - k / 2; L.normal[chr_y] = k / 2; }
    else { L.normal[chr_x] = k; L.normal[chr_y] = 0; }
    if (sexline.undisclosed) { L.sex_known = false; L.normal[chr_x] = L.normal[chr_y] = -1; }

    // The normal complement present before the abnormalities.
    std::vector<Chromosome>& cs = L.chromosomes;
    std::array<int, chromosome_count> serial{};
    const auto add_normal = [&](int c) {
      Chromosome ch;
      ch.segments = in.normal(c);
      ch.centromere_of = c;
      ch.label = chromosome_name(c);
      ch.key = chromosome_name(c) + "#" + std::to_string(++serial[c]);
      cs.push_back(std::move(ch));
    };
    for (int c = 0; c < 22; c++)
      for (int j = 0; j < k; j++) add_normal(c);
    std::string sex = sexline.sex;
    for (char s : sex) add_normal(s == 'X' ? chr_x : chr_y);
    const auto remove_normal = [&](int c) -> bool {
      for (auto it = cs.rbegin(); it != cs.rend(); ++it)
        if (!it->abnormal && it->centromere_of == c) {
          cs.erase(std::next(it).base());
          return true;
        }
      return false;
    };

    int abnormal_serial = 0;
    for (std::size_t ai = 0; ai < L.abnormalities.size(); ai++) {
      const Abnormality& a = L.abnormalities[ai];
      EventResult ev;
      ev.abnormality = ai;
      std::vector<Chromosome> products;
      std::vector<int> replaced;
      const int m = std::max(1, a.multiplier);
      if (a.numerical()) {
        const int c = a.chromosome();
        if (c < 0) continue;
        if (a.sign > 0) {
          for (int j = 0; j < m; j++) {
            add_normal(c);
            ev.products.push_back(cs.size() - 1);
          }
        } else {
          // −X and −Y of an acquired loss are already reflected in the sex chromosomes written,
          // unless they are marked constitutional (46,XXYc,−X) or inherited from another clone.
          const bool implicit = is_sex(c) && line.has_sex && !line.sex_constitutional;
          if (!implicit && !remove_normal(c))
            in.warning(a.span, "there is no normal chromosome " + chromosome_name(c) + " left to lose");
        }
        // The event's loss or gain: the whole chromosome.
        std::vector<char> v(g.bands.size(), 0);
        for (std::size_t b = g.begin_of(c); b < g.end_of(c); b++) v[b] = 1;
        (a.sign > 0 ? ev.lgf.gain : ev.lgf.loss) = detail::marked(v);
        L.events.push_back(std::move(ev));
        continue;
      }
      if (a.type == Type::marker || a.type == Type::ring_marker || a.type == Type::dmin) {
        const int n = (a.type == Type::dmin ? a.count : a.count) * m;
        for (int j = 0; j < n; j++) {
          Chromosome ch;
          ch.segments = {in.unknown_material()};
          ch.abnormal = true;
          ch.marker = true;
          ch.ring = a.type == Type::ring_marker;
          ch.counted = a.type != Type::dmin;
          Abnormality bare = a;
          bare.sign = 0;
          bare.count = bare.count_hi = 1;
          bare.multiplier = 1;
          ch.label = format(bare);
          ch.source = static_cast<int>(ai);
          ch.key = ch.label + "#" + std::to_string(++abnormal_serial);
          cs.push_back(std::move(ch));
          ev.products.push_back(cs.size() - 1);
        }
        L.unknown += n;
        L.events.push_back(std::move(ev));
        continue;
      }
      if (a.sign < 0) {
        // A missing abnormal chromosome not found in the line it refers to.
        in.warning(a.span, a.text + " is not among the chromosomes of the clone it refers to");
        L.events.push_back(std::move(ev));
        continue;
      }
      detail::Interpreter::Products p = in.realise(a);
      std::vector<Chromosome> made;
      for (int j = 0; j < m; j++)
        for (auto ch : p.chromosomes) {
          ch.source = static_cast<int>(ai);
          ch.key = ch.label + "#" + std::to_string(++abnormal_serial);
          made.push_back(std::move(ch));
        }
      std::vector<int> repl;
      if (a.sign <= 0)  // with + it is additional and replaces nothing
        for (int j = 0; j < m; j++) repl.insert(repl.end(), p.replaced.begin(), p.replaced.end());
      if (a.lost) made.clear();
      // The event's own loss, gain and fusions, against the normal chromosomes it replaced.
      {
        std::vector<char> lv(g.bands.size(), 0), gv(g.bands.size(), 0), fv(g.bands.size(), 0);
        std::vector<Chromosome> norms;
        for (int c : repl) {
          Chromosome n;
          n.segments = in.normal(c);
          norms.push_back(std::move(n));
        }
        std::set<int> touched(repl.begin(), repl.end());
        for (const auto& ch : made)
          for (const auto& s : ch.segments)
            if (s.kind == Segment::Kind::reference) touched.insert(s.chromosome);
        for (int c : touched) {
          std::vector<const Chromosome*> mine, theirs;
          for (const auto& ch : made) mine.push_back(&ch);
          for (const auto& ch : norms) theirs.push_back(&ch);
          const auto r1 = detail::coverage(mine, c, g), r0 = detail::coverage(theirs, c, g);
          // Compare piecewise on the union of the boundaries.
          std::vector<std::pair<std::uint64_t, std::pair<std::uint64_t, std::uint64_t>>> cuts;  // position, its uncertainty
          for (const auto& r : r1) { cuts.push_back({r.start, {r.first.lo, r.first.hi}}); cuts.push_back({r.end, {r.last.lo, r.last.hi}}); }
          for (const auto& r : r0) { cuts.push_back({r.start, {r.first.lo, r.first.hi}}); cuts.push_back({r.end, {r.last.lo, r.last.hi}}); }
          std::sort(cuts.begin(), cuts.end());
          const auto value = [](const std::vector<detail::Run>& rs, std::uint64_t x) {
            for (const auto& r : rs)
              if (r.start <= x && x < r.end) return r.copies;
            return 0;
          };
          for (std::size_t q = 0; q + 1 < cuts.size(); q++) {
            const std::uint64_t s0 = cuts[q].first, s1 = cuts[q + 1].first;
            if (s1 <= s0) continue;
            const int diff = value(r1, s0) - value(r0, s0);
            if (diff == 0) continue;
            const std::uint64_t lo = std::min(s0, cuts[q].second.first), hi = std::max(s1, cuts[q + 1].second.second);
            detail::mark(diff < 0 ? lv : gv, g, c, lo, hi);
          }
        }
        std::vector<Fusion> fs;
        for (std::size_t q = 0; q < made.size(); q++) detail::junctions(made[q], q, fs);
        for (const auto& f : fs) {
          detail::mark_point(fv, ev.lgf.fusion_labels, f.a, g);
          detail::mark_point(fv, ev.lgf.fusion_labels, f.b, g);
        }
        ev.lgf.loss = detail::marked(lv);
        ev.lgf.gain = detail::marked(gv);
        ev.lgf.fusion = detail::marked(fv);
      }
      for (int c : repl) {
        if (is_sex(c)) continue;  // the sex chromosomes written are the normal ones present
        if (!remove_normal(c))
          in.warning(a.span, "there is no normal chromosome " + chromosome_name(c) + " left for " + a.text + " to replace");
      }
      if (a.lost) made.clear();
      for (auto& ch : made) {
        cs.push_back(std::move(ch));
        ev.products.push_back(cs.size() - 1);
        for (const auto& s : cs.back().segments)
          if (s.kind != Segment::Kind::reference) L.unknown++;
      }
      L.events.push_back(std::move(ev));
    }

    // The count.
    L.counted = 0;
    for (const auto& ch : cs) L.counted += ch.counted ? 1 : 0;
    if (!line.count_unknown && line.count > 0) {
      const int lo = line.count, hi = std::max(line.count, line.count_hi);
      if (L.counted < lo || L.counted > hi) {
        const std::string written = line.count_hi > line.count ? std::to_string(lo) + "~" + std::to_string(hi) : std::to_string(lo);
        const std::string msg = "the abnormalities give " + std::to_string(L.counted) + " chromosomes, and the count is " + written;
        if (line.composite || line.incomplete) in.note(line.span, msg + " (a composite or incomplete karyotype need not agree)");
        else in.warning(line.span, msg);
      }
    }

    // Fusions, copy number, changes, bands.
    for (std::size_t q = 0; q < cs.size(); q++)
      if (cs[q].abnormal) detail::junctions(cs[q], q, L.fusions);
    L.bands.resize(g.bands.size());
    for (int c = 0; c < chromosome_count; c++) {
      if (!g.has(c)) continue;
      std::vector<const Chromosome*> all;
      for (const auto& ch : cs) all.push_back(&ch);
      const auto runs = detail::coverage(all, c, g);
      const int normal = L.normal[c];
      for (std::size_t b = g.begin_of(c); b < g.end_of(c); b++) {
        L.bands[b].band = b;
        L.bands[b].min = 1 << 30;
        L.bands[b].max = -1;
        L.bands[b].normal = normal;
      }
      for (const auto& r : runs) {
        const std::uint64_t olo = std::min(r.start, r.first.lo), ohi = std::max(r.end, r.last.hi);
        const BandRange br = g.bands_over(c, olo, ohi > olo ? ohi : olo + 1);
        for (std::size_t b = br.first; b <= br.last; b++) {
          // A run only meets a band that it overlaps in its outer extent.
          if (g.bands[b].end <= olo || g.bands[b].start >= ohi) continue;
          L.bands[b].min = std::min(L.bands[b].min, r.copies);
          L.bands[b].max = std::max(L.bands[b].max, r.copies);
        }
        if (normal < 0 || r.copies == normal) continue;
        Change ch;
        ch.chromosome = c;
        ch.copies = r.copies;
        ch.normal = normal;
        ch.start = r.start;
        ch.end = r.end;
        ch.outer_start = olo;
        ch.outer_end = ohi;
        ch.inner_start = std::max(r.start, r.first.hi);
        ch.inner_end = std::min(r.end, r.last.lo);
        if (ch.inner_end < ch.inner_start) ch.inner_end = ch.inner_start;
        ch.from_band = detail::end_band(r.first, c, g, true);
        ch.to_band = detail::end_band(r.last, c, g, false);
        ch.bands = g.bands_over(c, olo, ohi);
        L.changes.push_back(std::move(ch));
        L.balanced = false;
      }
      for (std::size_t b = g.begin_of(c); b < g.end_of(c); b++)
        if (L.bands[b].max < 0) L.bands[b].min = L.bands[b].max = 0;
    }

    // The line's LGF: the union of its events'.
    {
      std::vector<char> lv(g.bands.size(), 0), gv(g.bands.size(), 0), fv(g.bands.size(), 0);
      std::set<std::string> labels;
      for (const auto& e : L.events) {
        for (auto b : e.lgf.loss) lv[b] = 1;
        for (auto b : e.lgf.gain) gv[b] = 1;
        for (auto b : e.lgf.fusion) fv[b] = 1;
        labels.insert(e.lgf.fusion_labels.begin(), e.lgf.fusion_labels.end());
      }
      L.lgf.loss = detail::marked(lv);
      L.lgf.gain = detail::marked(gv);
      L.lgf.fusion = detail::marked(fv);
      L.lgf.fusion_labels.assign(labels.begin(), labels.end());
    }
    out.lines.push_back(std::move(L));
  }

  // The other techniques.
  for (std::size_t ti = 0; ti < d.techniques.size(); ti++) {
    const Technique& tq = d.techniques[ti];
    if (tq.name != "arr" && tq.name != "sseq" && tq.name != "seq") continue;
    TechniqueInterpretation T;
    T.technique = ti;
    T.assembly = tq.build.empty() ? g.assembly : assembly_of(tq.build);
    if (T.assembly != g.assembly)
      T.diagnostics.push_back({Diagnostic::Level::warning, tq.span.begin, tq.span.end,
                               "the result is on " + T.assembly + " and the ideogram on " + g.assembly + ": its bands are not checked"});
    for (const auto& e : tq.arr) {
      if (e.kind == ArrEntry::Kind::complex || e.copies < 0) continue;
      const auto add = [&](int c, std::uint64_t start, std::uint64_t end, const std::string& from, const std::string& to) {
        Change ch;
        ch.chromosome = c;
        ch.copies = e.copies;
        ch.normal = is_sex(c) ? -1 : 2;
        ch.start = ch.inner_start = ch.outer_start = start;
        ch.end = ch.inner_end = ch.outer_end = end;
        ch.from_band = from;
        ch.to_band = to;
        if (g.has(c) && end <= g.length[c] && end > start) ch.bands = g.bands_over(c, start, end);
        T.changes.push_back(ch);
      };
      if (e.kind == ArrEntry::Kind::whole) {
        for (const auto& [x, y] : e.wholes)
          for (int c = x; c <= y && c >= 0; c++)
            if (g.has(c)) add(c, 0, g.length[c], g.bands[g.pter_band(c)].name(), g.bands[g.qter_band(c)].name());
        continue;
      }
      if (!g.has(e.chromosome) || e.start == 0) continue;
      const std::uint64_t start = e.start - 1, end = e.end;
      add(e.chromosome, start, end, band_label(e.arm1, e.band1), e.arm2 ? band_label(e.arm2, e.band2) : band_label(e.arm1, e.band1));
      // The bands written hold the coordinates (ISCN 2024 ch. 8), on the ideogram of the build.
      if (T.assembly == g.assembly && end <= g.length[e.chromosome]) {
        const auto r1 = g.find(e.chromosome, e.arm1, e.band1);
        const auto r2 = e.arm2 ? g.find(e.chromosome, e.arm2, e.band2) : r1;
        const bool holds = r1 && r2 && r1->start <= start && start < r1->end && r2->start < end && end <= r2->end;
        if (!holds)
          T.diagnostics.push_back({Diagnostic::Level::warning, e.span.begin, e.span.end,
                                   "the bands written, " + chromosome_name(e.chromosome) + e.band_name() + ", do not hold " + with_commas(e.start) + "-" +
                                       with_commas(e.end) + ", which lie in " + chromosome_name(e.chromosome) + g.range_name(e.chromosome, start, end)});
      }
    }
    if (!tq.structural.empty()) {
      Description sd;
      CellLine line;
      line.count_unknown = true;
      line.undisclosed = true;
      line.has_sex = false;
      line.abnormalities = tq.structural;
      line.span = tq.span;
      sd.lines.push_back(line);
      Interpretation sub = interpret(sd, g);
      if (!sub.lines.empty()) {
        T.has_line = true;
        T.line = std::move(sub.lines[0]);
      }
    }
    out.techniques.push_back(std::move(T));
  }
  return out;
}

// ---- the detailed system written from the interpretation ----------------------------------------------

struct LongFormOptions {
  bool unicode = true;        // → (ISCN) or ->
  bool chromosome_names = true;
};

// One chromosome of the cell from end to end: 22pter→22q11.2::9q34→9qter.
inline std::string long_form(const Chromosome& ch, const LongFormOptions& o = {}) {
  const std::string arrow = o.unicode ? "\xE2\x86\x92" : "->";
  const auto label = [&](const Segment& s, const Point& p, bool lower) {
    std::string out = o.chromosome_names ? chromosome_name(s.chromosome) : "";
    if (p.kind == Point::Kind::cen) {
      // The centromere is written by the arm the segment lies on.
      const bool on_p = std::min(s.from.at, s.to.at) < p.at;
      return out + (on_p ? "p10" : "q10");
    }
    out += p.label;
    if (p.written_position) out += "(" + with_commas(lower ? p.at + 1 : p.at) + ")";
    return out;
  };
  std::string s = ch.ring ? "::" : ch.open;
  for (std::size_t i = 0; i < ch.segments.size(); i++) {
    if (i) s += "::";
    const Segment& x = ch.segments[i];
    switch (x.kind) {
      case Segment::Kind::reference: {
        const bool fwd = x.from.at <= x.to.at;
        s += label(x, x.from, fwd) + arrow + label(x, x.to, !fwd);
        break;
      }
      case Segment::Kind::unknown: s += "?"; break;
      case Segment::Kind::hsr: s += "hsr"; break;
      case Segment::Kind::bases: s += x.bases; break;
    }
  }
  s += ch.ring ? "::" : ch.close;
  return s;
}

// The cell line in the detailed system: each structural abnormality's breakpoints replaced by its
// chromosomes written from end to end (ISCN 2024 §4.3; its chromosome numbers left out when it
// involves one chromosome only).
inline std::string detailed(const LineInterpretation& L, const CellLine& line, const LongFormOptions& o = {}) {
  FormatOptions f;
  f.unicode = o.unicode;
  std::string s = line.count_unknown ? "?" : std::to_string(line.count);
  if (line.count_hi > line.count) s += "~" + std::to_string(line.count_hi);
  if (line.ploidy > 0) s += "<" + std::to_string(line.ploidy) + "n>";
  if (line.has_sex) s += "," + (line.undisclosed ? std::string("U") : line.sex + (line.sex_uncertain ? "?" : "") + (line.sex_constitutional ? "c" : ""));
  for (std::size_t ai = 0; ai < L.abnormalities.size(); ai++) {
    const Abnormality& a = L.abnormalities[ai];
    s += ",";
    if (a.numerical() || a.type == Type::marker || a.type == Type::ring_marker || a.type == Type::dmin || L.events[ai].products.empty()) {
      s += format(a, f);
      continue;
    }
    if (a.sign > 0) s += "+";
    for (const auto& p : a.prefixes) s += p + " ";
    s += std::string(type_name(a.type)) + "(";
    for (std::size_t i = 0; i < a.chromosomes.size(); i++) s += (i ? ";" : "") + a.chromosomes[i].name();
    s += ")(";
    std::set<int> refs;
    const auto& prods = L.events[ai].products;
    const std::size_t per = prods.size() / static_cast<std::size_t>(std::max(1, a.multiplier));
    for (std::size_t q = 0; q < per; q++)
      for (const auto& seg : L.chromosomes[prods[q]].segments)
        if (seg.kind == Segment::Kind::reference) refs.insert(seg.chromosome);
    LongFormOptions lo = o;
    lo.chromosome_names = refs.size() > 1;
    for (std::size_t q = 0; q < per; q++) s += (q ? ";" : "") + long_form(L.chromosomes[prods[q]], lo);
    s += ")";
    if (a.multiplier > 1) s += (o.unicode ? "\xC3\x97" : "x") + std::to_string(a.multiplier);
    if (a.constitutional) s += "c";
    s += a.inheritance;
  }
  if (line.incomplete) s += ",inc";
  if (line.cells >= 0) s += "[" + std::string(line.composite ? "cp" : "") + std::to_string(line.cells) + "]";
  return s;
}

}  // namespace genomics::karyotype
