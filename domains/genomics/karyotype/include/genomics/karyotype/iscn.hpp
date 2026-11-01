// ISCN descriptions as typed values: the parse tree of grammar.hpp read into cell lines,
// abnormalities, breakpoints, detailed forms, microarray and sequencing results, with every
// fault of the text reported where it is.
//
// - Syntax. morphe's BNF reader decides whether the text is in the language and, when it is not,
//   says where the longest viable prefix ends and what could have followed it.
// - Sense. A text in the language may still not describe anything: t(9;22) with one group of
//   breakpoints, del(5)(q13q21q31) with three, a chromosome 23, components after a del. These are
//   errors with their spans; the description is still returned, for what it says.
// - Editions. ISCN 2020 lists the breakpoints of a rearrangement within one chromosome from pter to
//   qter, and the inserted segment of an insertion with the band nearer the recipient's pter first;
//   ISCN 2016 and earlier listed them from the centromere outward (Miller and Madan 2021, on
//   ISCN 2020 §4.3.1.1, §4.3.1.2 and §9.2.5, §9.2.9). The two editions agree on the q arm and differ
//   on the p arm, which matters for duplications (direct or inverted) and insertions (direct or
//   inverted); a text read as ISCN 2016 is turned into ISCN 2020's order here, so that every value
//   below means what ISCN 2024 means.
//
// Complexity: linear in the size of the parse tree.

#pragma once

#include <genomics/karyotype/bands.hpp>
#include <genomics/karyotype/grammar.hpp>

#include <morphe/formal/bnf.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace genomics::karyotype {

struct Span {
  std::size_t begin = 0, end = 0;
};

struct Diagnostic {
  enum class Level : std::uint8_t { error, warning, note };
  Level level = Level::error;
  std::size_t begin = 0, end = 0;
  std::string message;
};
inline std::string_view level_name(Diagnostic::Level l) {
  return l == Diagnostic::Level::error ? "error" : l == Diagnostic::Level::warning ? "warning" : "note";
}

enum class Edition : std::uint8_t { iscn2020, iscn2016 };

// A breakpoint as written: a band at any resolution (or an interval of bands), pter, qter, cen, or
// unknown (?); with the chromosome when the detailed system names it, and with nucleotide
// positions in sequencing results.
struct Breakpoint {
  enum class Kind : std::uint8_t { band, pter, qter, cen, unknown };
  Kind kind = Kind::band;
  int chromosome = -1;    // −1: the abnormality's chromosome
  char arm = 0;           // 'p' or 'q'
  std::string digits;     // "3633"; empty with `uncertain`: somewhere on the arm
  bool uncertain = false; // q2? (a band of region 2), q? (somewhere on q)
  char to_arm = 0;        // an interval of bands, q21~q24 or q11.2-13.1: its far end
  std::string to_digits;
  char to_sep = '~';      // as written: '~' or '-' (ISCN 2020's hyphen of a band at low resolution)
  bool to_arm_written = false;
  bool has_position = false, position_open = false;  // seq: (164,789,896), (128,534,000∼128,546,000), (158,573,186∼)
  std::uint64_t position = 0, position_to = 0;       // 1-based, as written
  Span span;

  bool centromeric() const { return kind == Kind::band && digits == "10" && !uncertain && to_digits.empty(); }
  bool interval() const { return !to_digits.empty(); }
  std::string label() const {
    switch (kind) {
      case Kind::pter: return "pter";
      case Kind::qter: return "qter";
      case Kind::cen: return "cen";
      case Kind::unknown: return "?";
      case Kind::band: break;
    }
    std::string s = band_label(arm, digits);
    if (uncertain) s += "?";
    if (interval()) s += std::string(1, to_sep) + (to_arm_written ? band_label(to_arm, to_digits) : band_label(arm, to_digits).substr(1));
    return s;
  }
};

enum class Type : std::uint8_t {
  gain, loss,  // numerical: +21, −7
  add, del, der, dic, dup, fis, fra, hsr, i, idic, ider, ins, inv, neo, qdp, r, rec, rob, t, tas, trc, trp,
  marker, ring_marker, dmin
};
inline std::string_view type_name(Type t) {
  static constexpr std::string_view names[] = {"+", "-", "add", "del", "der", "dic", "dup", "fis", "fra", "hsr", "i", "idic", "ider", "ins", "inv",
                                               "neo", "qdp", "r", "rec", "rob", "t", "tas", "trc", "trp", "mar", "r", "dmin"};
  return names[static_cast<int>(t)];
}
inline std::optional<Type> type_from(std::string_view a) {
  static constexpr std::pair<std::string_view, Type> table[] = {
      {"add", Type::add}, {"del", Type::del}, {"der", Type::der}, {"dic", Type::dic}, {"dup", Type::dup}, {"fis", Type::fis},
      {"fra", Type::fra}, {"hsr", Type::hsr}, {"i", Type::i}, {"idic", Type::idic}, {"ider", Type::ider}, {"ins", Type::ins},
      {"inv", Type::inv}, {"neo", Type::neo}, {"qdp", Type::qdp}, {"r", Type::r}, {"rec", Type::rec}, {"rob", Type::rob},
      {"t", Type::t}, {"tas", Type::tas}, {"trc", Type::trc}, {"trp", Type::trp}};
  for (const auto& [n, t] : table)
    if (n == a) return t;
  return std::nullopt;
}

struct ChromosomeRef {
  int index = -1;
  bool unknown = false;    // ?
  bool homologue = false;  // 7*: the homologue of a chromosome named before (ISCN's underlining)
  std::string name() const { return unknown ? "?" : chromosome_name(index) + (homologue ? "*" : ""); }
};

// A part of a chromosome in the detailed system: a segment a→b, material of unknown origin (?),
// nucleotides that align nowhere (ACG…), or a region with its copies (q11.21q24.23(…)×2).
struct LongPart {
  enum class Kind : std::uint8_t { segment, unknown, bases, region };
  Kind kind = Kind::segment;
  Breakpoint from, to;
  std::string bases;
  int multiplier = 1;
};
// A chromosome in the detailed system: its parts joined by ::, opened and closed by "", ":" (a
// break) or "::" (a ring).
struct LongChromosome {
  std::string open, close;
  std::vector<LongPart> parts;
};

struct Abnormality {
  Type type = Type::gain;
  int sign = 0;                        // +1: an additional chromosome; −1: a missing one
  bool uncertain = false;              // ?del(…)
  std::vector<std::string> prefixes;   // psu, inv (inv dup), dir, dic (dic r), trc, der (der dic)
  std::vector<ChromosomeRef> chromosomes;
  bool has_groups = false;
  std::vector<std::vector<Breakpoint>> groups;  // one group per chromosome, in its order
  std::vector<LongChromosome> long_form;        // when written in the detailed system
  std::vector<Abnormality> components;          // of der, ider, rec
  bool rec_dup = true;                          // rec(2)dup(2q): the duplicated arm
  int rec_chromosome = -1;
  char rec_arm = 0;
  int multiplier = 1;                           // x2
  int count = 1, count_hi = 1;                  // +2mar; 3~9dmin
  int marker_number = 0;                        // mar1
  bool constitutional = false;                  // c
  std::string inheritance;                      // dn, mat, pat, inh, dmat, dpat, dinh
  std::vector<Abnormality> alternatives;        // … or …
  bool lost = false;                            // in a subclone: this abnormal chromosome of the clone it refers to is missing (−dic(1;20))
  Span span;
  std::string text;

  bool numerical() const { return type == Type::gain || type == Type::loss; }
  // The chromosome by which the abnormality is ordered and named: its first.
  int chromosome() const {
    for (const auto& c : chromosomes)
      if (!c.unknown) return c.index;
    return -1;
  }
  bool has_prefix(std::string_view p) const { return std::find(prefixes.begin(), prefixes.end(), p) != prefixes.end(); }
};

struct ArrEntry {
  enum class Kind : std::uint8_t { region, whole, complex };
  Kind kind = Kind::region;
  int chromosome = -1;
  char arm1 = 0, arm2 = 0;
  std::string band1, band2;                // words of digits
  std::uint64_t start = 0, end = 0;        // 1-based, closed, as written
  bool extended = false;                   // with the normal flanks: (a×2,b_c×1,d×2)
  std::uint64_t flank_start = 0, flank_end = 0;
  int flank_start_copies = -1, flank_end_copies = -1;
  int copies = -1, copies_hi = -1;         // x1, x2~3
  double proportion = -1;                  // [0.5]
  std::vector<std::pair<int, int>> wholes; // (X,1-22): ranges of chromosomes
  std::vector<int> complex_chromosomes;    // (2,3,8)cx
  std::vector<std::string> flags;          // hmz, htz, cth, cha, dn, mat, pat, …
  Span span;
  std::string band_name() const { return (arm1 ? band_label(arm1, band1) : "") + (arm2 ? band_label(arm2, band2) : ""); }
};

// A junction of a sequencing result: 1q23.3(+)(164,789,896)::8q12.2(+)(61,760,650).
struct Junction {
  struct End {
    int chromosome = -1;
    char arm = 0;
    std::string digits;
    int strand = 0;  // +1, −1, 0 when not written
    std::uint64_t position = 0;  // 1-based
  };
  End a, b;
  std::string bases;
  Span span;
};

struct Technique {
  std::string name;   // arr, sseq, seq, ish, nuc ish, rsa
  std::string build;  // GRCh38, GRCh37/hg19, …
  std::vector<ArrEntry> arr;
  std::vector<Abnormality> structural;
  std::vector<Junction> junctions;
  std::string text;
  Span span;
};

struct CellLine {
  bool count_unknown = false;
  int count = -1, count_hi = -1;
  char count_range = 0;   // '~' or '-'
  int ploidy = 0;         // <3n>
  bool has_sex = false;
  std::string sex;        // the sex chromosomes written: XX, XY, X, XXY, …
  bool undisclosed = false, sex_constitutional = false, sex_uncertain = false;
  enum class Reference : std::uint8_t { none, idem, sl, sdl };
  Reference reference = Reference::none;
  int sdl_number = 0;
  std::vector<Abnormality> abnormalities;
  bool incomplete = false;
  int cells = -1;
  bool composite = false;
  std::string separator;  // before this line: "", "/" or "//"
  Span span;
};

struct Description {
  std::string mosaic;  // mos, chi or empty
  std::vector<CellLine> lines;
  std::vector<Technique> techniques;
};

struct Parsed {
  morphe::formal::BnfParse syntax;
  Description description;
  std::vector<Diagnostic> diagnostics;
  bool accepted() const { return syntax.accepted; }
  bool ok() const {
    if (!syntax.accepted) return false;
    for (const auto& d : diagnostics)
      if (d.level == Diagnostic::Level::error) return false;
    return true;
  }
  std::string errors() const {
    std::string s;
    for (const auto& d : diagnostics)
      if (d.level == Diagnostic::Level::error) s += (s.empty() ? "" : "; ") + d.message;
    return s;
  }
};

namespace detail {

using Node = morphe::formal::ParseNode;

inline int chromosome_of(const Node& n) { return chromosome_index(n.text); }

inline std::uint64_t coord_value(std::string_view s) {
  std::uint64_t v = 0;
  for (char c : s)
    if (c >= '0' && c <= '9') v = v * 10 + static_cast<std::uint64_t>(c - '0');
  return v;
}
inline int small_number(std::string_view s) {
  long v = 0;
  for (char c : s) {
    if (c < '0' || c > '9') break;
    v = v * 10 + (c - '0');
    if (v > 1000000) return 1000000;
  }
  return static_cast<int>(v);
}

class Builder {
 public:
  Builder(std::string_view text, Edition edition, std::vector<Diagnostic>& out) : text_(text), edition_(edition), out_(out) {}

  void error(std::size_t b, std::size_t e, std::string m) { out_.push_back({Diagnostic::Level::error, b, e, std::move(m)}); }
  void warning(std::size_t b, std::size_t e, std::string m) { out_.push_back({Diagnostic::Level::warning, b, e, std::move(m)}); }
  void note(std::size_t b, std::size_t e, std::string m) { out_.push_back({Diagnostic::Level::note, b, e, std::move(m)}); }

  int chromosome(const Node& n) {
    const int c = chromosome_of(n);
    if (c < 0) error(n.begin, n.end, "there is no chromosome " + n.text + " in the human karyotype (1–22, X, Y)");
    return c;
  }

  Breakpoint breakpoint(const Node& n) {
    Breakpoint b;
    b.span = {n.begin, n.end};
    if (const Node* t = n.child("terminal")) {
      b.kind = t->text == "pter" ? Breakpoint::Kind::pter : t->text == "qter" ? Breakpoint::Kind::qter : Breakpoint::Kind::cen;
    } else {
      b.kind = Breakpoint::Kind::band;
      if (const Node* a = n.child("arm")) b.arm = a->text[0];
      if (const Node* d = n.child("band_digits")) b.digits = digits(*d);
      b.uncertain = n.child("uncertain") != nullptr;
      if (const Node* e = n.child("band_end")) {
        const Node* ea = e->child("arm");
        b.to_arm = ea ? ea->text[0] : b.arm;
        b.to_arm_written = ea != nullptr;
        if (const Node* rs = n.child("range_sep")) b.to_sep = rs->text == "-" ? '-' : '~';
        if (const Node* ed = e->child("band_digits")) b.to_digits = digits(*ed);
        if (b.to_arm != b.arm) error(n.begin, n.end, "the interval " + n.text + " runs from one arm to the other");
      }
    }
    if (const Node* c = n.child("chromosome")) b.chromosome = chromosome(*c);
    if (const Node* lc = n.child("lf_coord")) coordinates(*lc, b);
    return b;
  }

  std::string digits(const Node& d) {
    const auto w = band_digits(d.text);
    if (!w) {
      error(d.begin, d.end, "not a band: " + d.text + " (a band is written with its region and band, and sub-bands after a point: 36.33)");
      return {};
    }
    return *w;
  }

  void coordinates(const Node& lc, Breakpoint& b) {
    const auto cs = lc.all("coord");
    if (cs.empty()) return;
    b.has_position = true;
    b.position = coord_value(cs[0]->text);
    if (cs.size() > 1) {
      b.position_to = coord_value(cs[1]->text);
    } else if (lc.text.size() > cs[0]->text.size() + 2) {
      b.position_open = true;  // (158,573,186∼)
    }
    for (const Node* c : cs) check_grouping(*c);
  }

  void check_grouping(const Node& c) {
    // A number written with commas groups its digits by three (ISCN 2024 ch. 8).
    if (c.text.find(',') != std::string::npos && c.text.find(',') > 3)
      error(c.begin, c.end, "the commas of " + c.text + " do not group its digits by three");
  }

  LongChromosome long_chromosome(const Node& n) {
    LongChromosome lc;
    const auto ends = n.all("lf_end");
    for (const Node* e : ends) {
      if (e->begin == n.begin) lc.open = e->text;
      else lc.close = e->text;
    }
    for (const auto& c : n.children) {
      if (c.name == "lf_segment") {
        LongPart p;
        p.kind = LongPart::Kind::segment;
        const auto pts = c.all("lf_point");
        if (pts.size() == 2) {
          p.from = breakpoint(*pts[0]);
          p.to = breakpoint(*pts[1]);
        }
        lc.parts.push_back(std::move(p));
      } else if (c.name == "lf_unknown") {
        LongPart p;
        p.kind = LongPart::Kind::unknown;
        lc.parts.push_back(std::move(p));
      } else if (c.name == "lf_bases") {
        LongPart p;
        p.kind = LongPart::Kind::bases;
        p.bases = c.text;
        lc.parts.push_back(std::move(p));
      } else if (c.name == "lf_region") {
        LongPart p;
        p.kind = LongPart::Kind::region;
        const auto arms = c.all("arm");
        const auto ds = c.all("band_digits");
        int chrom = -1;
        if (const Node* ch = c.child("chromosome")) chrom = chromosome(*ch);
        if (arms.size() == 2 && ds.size() == 2) {
          p.from.kind = p.to.kind = Breakpoint::Kind::band;
          p.from.arm = arms[0]->text[0];
          p.to.arm = arms[1]->text[0];
          p.from.digits = digits(*ds[0]);
          p.to.digits = digits(*ds[1]);
          p.from.chromosome = p.to.chromosome = chrom;
        }
        if (const Node* lco = c.child("lf_coord")) {
          const auto cs = lco->all("coord");
          if (cs.size() >= 2) {
            p.from.has_position = p.to.has_position = true;
            p.from.position = coord_value(cs[0]->text);
            p.to.position = coord_value(cs[1]->text);
          }
        }
        if (const Node* m = c.child("multiplier")) p.multiplier = small_number(m->text.substr(m->text.find_first_of("0123456789")));
        lc.parts.push_back(std::move(p));
      }
    }
    return lc;
  }

  // The parts of a head or component: abbreviation, chromosomes, breakpoints or detailed form.
  void part(const Node& n, Abnormality& a) {
    if (const Node* ab = n.child("abbr")) {
      const auto t = type_from(ab->text);
      if (t) a.type = *t;
    }
    if (const Node* ch = n.child("chroms")) {
      for (const Node* r : ch->all("chrom_ref")) {
        ChromosomeRef ref;
        if (r->child("unknown")) ref.unknown = true;
        else if (const Node* c = r->child("chromosome")) ref.index = chromosome(*c);
        ref.homologue = r->child("homologue") != nullptr;
        a.chromosomes.push_back(ref);
      }
    }
    for (const auto& c : n.children) {
      if (c.name == "group") {
        a.has_groups = true;
        std::vector<Breakpoint> g;
        if (c.child("unknown")) {
          Breakpoint b;
          b.kind = Breakpoint::Kind::unknown;
          b.span = {c.begin, c.end};
          g.push_back(b);
        }
        for (const Node* bp : c.all("breakpoint")) g.push_back(breakpoint(*bp));
        a.groups.push_back(std::move(g));
      } else if (c.name == "long_body") {
        for (const Node* lc : c.all("long_chromosome")) a.long_form.push_back(long_chromosome(*lc));
      }
    }
  }

  Abnormality abnormality(const Node& n) {
    Abnormality a;
    a.span = {n.begin, n.end};
    a.text = n.text;
    if (n.name == "numeric") {
      const Node* s = n.child("sign");
      a.sign = s && s->text == "+" ? 1 : -1;
      a.type = a.sign > 0 ? Type::gain : Type::loss;
      ChromosomeRef ref;
      if (const Node* c = n.child("chromosome")) ref.index = chromosome(*c);
      a.chromosomes.push_back(ref);
    } else if (n.name == "structural") {
      if (const Node* s = n.child("sign")) a.sign = s->text == "+" ? 1 : -1;
      a.uncertain = n.child("uncertain") != nullptr;
      for (const Node* p : n.all("prefix")) {
        std::string w = p->text;
        while (!w.empty() && w.back() == ' ') w.pop_back();
        a.prefixes.push_back(w);
      }
      if (const Node* h = n.child("head")) part(*h, a);
      for (const Node* c : n.all("component")) {
        if (const Node* rp = c->child("rec_part")) {
          a.rec_dup = rp->text.substr(0, 3) == "dup";
          if (const Node* ch = rp->child("chromosome")) a.rec_chromosome = chromosome(*ch);
          if (const Node* ar = rp->child("arm")) a.rec_arm = ar->text[0];
          continue;
        }
        Abnormality comp;
        comp.span = {c->begin, c->end};
        comp.text = c->text;
        part(*c, comp);
        a.components.push_back(std::move(comp));
      }
      if (const Node* m = n.child("multiplier")) a.multiplier = std::max(1, small_number(m->text.substr(m->text.find_first_of("0123456789"))));
    } else if (n.name == "marker" || n.name == "ring_marker") {
      a.type = n.name == "marker" ? Type::marker : Type::ring_marker;
      if (const Node* s = n.child("sign")) a.sign = s->text == "+" ? 1 : -1;
      if (const Node* c = n.child("marker_count")) a.count = a.count_hi = std::max(1, small_number(c->text));
      if (const Node* c = n.child("marker_number")) a.marker_number = small_number(c->text);
      if (const Node* m = n.child("multiplier")) a.multiplier = std::max(1, small_number(m->text.substr(m->text.find_first_of("0123456789"))));
    } else if (n.name == "dmin") {
      a.type = Type::dmin;
      const auto ns = n.all("number");
      a.count = a.count_hi = ns.empty() ? 1 : small_number(ns[0]->text);
      if (ns.size() > 1) a.count_hi = small_number(ns[1]->text);
    }
    for (const Node* f : n.all("flag")) {
      std::string w = f->text;
      w.erase(0, w.find_first_not_of(' '));
      if (w == "c") a.constitutional = true;
      else if (!a.inheritance.empty()) error(f->begin, f->end, "two inheritances, " + a.inheritance + " and " + w);
      else a.inheritance = w;
    }
    return a;
  }

  // The rules of sense of one structural abnormality.
  void validate(Abnormality& a, bool component = false) {
    if (a.numerical() || a.type == Type::marker || a.type == Type::ring_marker || a.type == Type::dmin) return;
    const std::size_t nc = a.chromosomes.size();
    const auto fail = [&](const std::string& m) { error(a.span.begin, a.span.end, std::string(type_name(a.type)) + ": " + m); };
    const bool der_like = a.type == Type::der || a.type == Type::ider || a.type == Type::rec || a.has_prefix("der");
    if (!a.components.empty() && !der_like) fail("only a derivative chromosome (der, ider, rec) is followed by the abnormalities that made it");
    if (!a.long_form.empty()) {
      if (a.has_groups) fail("both breakpoints and a detailed description");
      return;
    }
    if (a.has_groups && a.groups.size() != nc && !(nc == 1 && a.type != Type::t))
      fail(std::to_string(nc) + " chromosome" + (nc == 1 ? "" : "s") + " but " + std::to_string(a.groups.size()) + " group" +
           (a.groups.size() == 1 ? "" : "s") + " of breakpoints; each chromosome has its group, separated by semicolons");
    const auto count_in = [&](std::size_t g) { return g < a.groups.size() ? a.groups[g].size() : 0; };
    std::size_t total = 0;
    for (const auto& g : a.groups) total += g.size();
    switch (a.type) {
      case Type::t:
        if (nc < 2) fail("a translocation involves two or more chromosomes");
        if (a.has_groups) {
          for (std::size_t g = 0; g < a.groups.size(); g++)
            if (count_in(g) != 1 && count_in(g) != 2) fail("each chromosome of a translocation has one breakpoint, or two for an exchanged segment");
          for (std::size_t g = 1; g < a.groups.size(); g++)
            if ((count_in(g) == 2) != (count_in(0) == 2)) fail("a translocation of segments has two breakpoints on every chromosome");
        } else if (!component && a.sign == 0) {
          // t(9;22) without breakpoints is allowed where the same translocation was written before.
          note(a.span.begin, a.span.end, "the breakpoints of " + a.text + " are taken from where the same rearrangement is written with them");
        }
        break;
      case Type::del:
        if (nc != 1) fail("a deletion is of one chromosome");
        if (a.has_groups && total != 1 && total != 2) fail("a deletion has one breakpoint (terminal) or two (interstitial)");
        if (!a.has_groups) fail("a deletion needs its breakpoints");
        break;
      case Type::dup: case Type::trp: case Type::qdp: case Type::inv:
        if (nc != 1) fail("it is of one chromosome");
        if (!a.has_groups || total != 2) fail("it needs two breakpoints, the ends of its segment");
        break;
      case Type::ins:
        if (nc == 1) {
          if (!a.has_groups || total != 3) fail("an insertion within a chromosome has three breakpoints: the point of insertion and the ends of the inserted segment");
        } else if (nc == 2) {
          if (!a.has_groups || count_in(0) != 1 || count_in(1) != 2)
            if (!(a.has_groups && count_in(1) == 1 && a.groups[1][0].kind == Breakpoint::Kind::unknown))
              fail("an insertion between two chromosomes has the point of insertion on the first and the ends of the segment on the second");
        } else {
          fail("an insertion involves one chromosome or two");
        }
        break;
      case Type::i: case Type::idic: case Type::add: case Type::hsr: case Type::fra: case Type::fis: case Type::neo:
        if (a.type == Type::hsr && nc == 1 && a.chromosomes[0].unknown && !a.has_groups) break;  // hsr(?)
        if (nc != 1 && !(a.type == Type::hsr && nc == 2)) fail("it is of one chromosome");
        if (!a.has_groups || (total != 1 && !(a.type == Type::hsr && total == 2))) fail("it needs one breakpoint");
        if (a.type == Type::i && a.has_groups && total == 1 && !a.groups[0][0].centromeric() && a.groups[0][0].kind == Breakpoint::Kind::band)
          warning(a.span.begin, a.span.end, "an isochromosome breaks in its centromere, p10 or q10; " + a.groups[0][0].label() + " is read as its point of fusion");
        if (a.type == Type::fis && a.has_groups && total == 1 && !a.groups[0][0].centromeric()) fail("a centric fission breaks at p10 or q10");
        break;
      case Type::r:
        if (a.has_groups && (a.groups.size() != nc || std::any_of(a.groups.begin(), a.groups.end(), [](const auto& g) { return g.size() != 2 && !(g.size() == 1 && g[0].kind == Breakpoint::Kind::unknown); }))) fail("each chromosome of a ring has two breakpoints, the ends of its segment");
        break;
      case Type::dic: case Type::rob: case Type::tas:
        if (nc < 2) fail("it joins two chromosomes or more");
        if (!a.has_groups) fail("it needs its breakpoints");
        break;
      case Type::trc:
        if (nc != 3) fail("a tricentric chromosome joins three chromosomes");
        if (!a.has_groups) fail("it needs its breakpoints");
        break;
      case Type::der: case Type::ider:
        if (nc == 0) fail("a derivative names its chromosome");
        if (a.components.empty() && !a.has_groups) fail("a derivative chromosome is described by its breakpoints or by the abnormalities that made it");
        if (a.type == Type::ider && !a.has_groups) fail("an isoderivative names the arm doubled, p10 or q10");
        for (auto& c : a.components) validate(c, true);
        break;
      case Type::rec:
        if (a.rec_chromosome < 0) fail("a recombinant chromosome names its duplicated or deleted arm, as dup(2q)");
        for (auto& c : a.components) validate(c, true);
        break;
      default:
        break;
    }
    if (a.has_prefix("der"))
      for (auto& c : a.components) validate(c, true);
  }

  // ISCN 2016's order (from the centromere) turned into ISCN 2020's (from pter): duplications and
  // the inserted segments of insertions on the short arm are reversed.
  void to_edition_2020(Abnormality& a) {
    if (edition_ != Edition::iscn2016) return;
    if ((a.type == Type::dup || a.type == Type::trp || a.type == Type::qdp) && a.has_groups && a.groups.size() == 1 && a.groups[0].size() == 2) {
      auto& g = a.groups[0];
      if (g[0].arm == 'p' && g[1].arm == 'p') std::swap(g[0], g[1]);
    }
    if (a.type == Type::ins && a.has_groups) {
      // ins(5;2)(p14;q22q32) or ins(2)(p13q21q31): the recipient's breakpoint decides.
      std::vector<Breakpoint>* seg = nullptr;
      const Breakpoint* at = nullptr;
      if (a.groups.size() == 2 && a.groups[0].size() == 1 && a.groups[1].size() == 2) { at = &a.groups[0][0]; seg = &a.groups[1]; }
      if (a.groups.size() == 1 && a.groups[0].size() == 3) { at = &a.groups[0][0]; seg = &a.groups[0]; }
      if (at && seg && at->arm == 'p') {
        auto& s = *seg;
        std::swap(s[s.size() - 2], s[s.size() - 1]);
      }
    }
    for (auto& c : a.components) to_edition_2020(c);
  }

  CellLine cell_line(const Node& n) {
    CellLine line;
    line.span = {n.begin, n.end};
    if (const Node* c = n.child("count")) {
      if (c->text == "?") {
        line.count_unknown = true;
      } else {
        const auto ns = c->all("number");
        if (!ns.empty()) line.count = line.count_hi = small_number(ns[0]->text);
        if (ns.size() > 1) {
          line.count_hi = small_number(ns[1]->text);
          line.count_range = c->text.find('-') != std::string::npos ? '-' : '~';
          if (line.count_hi < line.count) error(c->begin, c->end, "the range of the chromosome count " + c->text + " runs downward");
        }
        if (const Node* p = c->child("ploidy")) {
          if (const Node* pn = p->child("number")) line.ploidy = small_number(pn->text);
        }
      }
    }
    if (const Node* s = n.child("sex")) {
      line.has_sex = true;
      if (s->text == "U") {
        line.undisclosed = true;
      } else {
        if (const Node* sc = s->child("sex_chromosomes")) line.sex = sc->text;
        if (const Node* f = s->child("sex_flag")) {
          line.sex_constitutional = f->text.find('c') != std::string::npos;
          line.sex_uncertain = f->text.find('?') != std::string::npos;
        }
      }
    }
    if (const Node* st = n.child("stemline")) {
      if (st->text == "idem") line.reference = CellLine::Reference::idem;
      else if (st->text == "sl") line.reference = CellLine::Reference::sl;
      else {
        line.reference = CellLine::Reference::sdl;
        if (const Node* sn = st->child("number")) line.sdl_number = small_number(sn->text);
      }
    }
    for (const Node* it : n.all("item")) {
      const Node* first = nullptr;
      std::vector<const Node*> alts;
      for (const auto& c : it->children) {
        if (c.name == "alternative") {
          for (const auto& d : c.children)
            if (d.name == "numeric" || d.name == "structural" || d.name == "marker" || d.name == "ring_marker" || d.name == "dmin" || d.name == "inc") alts.push_back(&d);
        } else if (!first) {
          first = &c;
        }
      }
      if (!first) continue;
      if (first->name == "inc") { line.incomplete = true; continue; }
      Abnormality a = abnormality(*first);
      for (const Node* alt : alts) {
        if (alt->name == "inc") continue;
        Abnormality b = abnormality(*alt);
        validate(b);
        to_edition_2020(b);
        a.alternatives.push_back(std::move(b));
      }
      validate(a);
      to_edition_2020(a);
      line.abnormalities.push_back(std::move(a));
    }
    if (const Node* cl = n.child("cells")) {
      line.composite = cl->child("composite") != nullptr;
      if (const Node* cn = cl->child("number")) line.cells = small_number(cn->text);
    }
    return line;
  }

  ArrEntry arr_entry(const Node& n) {
    ArrEntry e;
    e.span = {n.begin, n.end};
    for (const auto& c : n.children) {
      if (c.name == "arr_region") {
        e.kind = ArrEntry::Kind::region;
        if (const Node* ch = c.child("chromosome")) e.chromosome = chromosome(*ch);
        const auto arms = c.all("arm");
        const auto ds = c.all("band_digits");
        if (!arms.empty() && !ds.empty()) { e.arm1 = arms[0]->text[0]; e.band1 = digits(*ds[0]); }
        if (arms.size() > 1 && ds.size() > 1) { e.arm2 = arms[1]->text[0]; e.band2 = digits(*ds[1]); }
        if (const Node* sp = c.child("arr_span")) {
          const auto cs = sp->all("coord");
          if (cs.size() == 2) { e.start = coord_value(cs[0]->text); e.end = coord_value(cs[1]->text); check_grouping(*cs[0]); check_grouping(*cs[1]); }
        } else if (const Node* ex = c.child("arr_extended")) {
          e.extended = true;
          const auto cs = ex->all("coord");
          const auto cps = ex->all("copies");
          if (cs.size() == 4 && cps.size() == 3) {
            e.flank_start = coord_value(cs[0]->text);
            e.start = coord_value(cs[1]->text);
            e.end = coord_value(cs[2]->text);
            e.flank_end = coord_value(cs[3]->text);
            int lo, hi;
            copies(*cps[0], e.flank_start_copies, hi, nullptr);
            copies(*cps[1], e.copies, e.copies_hi, &e.proportion);
            copies(*cps[2], e.flank_end_copies, lo, nullptr);
          }
        }
        if (const Node* cp = c.child("copies")) copies(*cp, e.copies, e.copies_hi, &e.proportion);
        if (e.start > e.end) error(c.begin, c.end, "the region " + c.text + " ends before it begins");
      } else if (c.name == "arr_whole") {
        e.kind = ArrEntry::Kind::whole;
        for (const Node* w : c.all("whole_ref")) {
          if (w->text == "XY") { e.wholes.push_back({chr_x, chr_x}); e.wholes.push_back({chr_y, chr_y}); continue; }
          const auto chs = w->all("chromosome");
          if (chs.size() == 1) { const int x = chromosome(*chs[0]); e.wholes.push_back({x, x}); }
          else if (chs.size() == 2) {
            const int x = chromosome(*chs[0]), y = chromosome(*chs[1]);
            if (y < x) error(w->begin, w->end, "the range of chromosomes " + w->text + " runs downward");
            e.wholes.push_back({x, y});
          }
        }
        if (const Node* cp = c.child("copies")) copies(*cp, e.copies, e.copies_hi, &e.proportion);
      } else if (c.name == "arr_complex") {
        e.kind = ArrEntry::Kind::complex;
        for (const Node* ch : c.all("chromosome")) e.complex_chromosomes.push_back(chromosome(*ch));
      } else if (c.name == "arr_flag") {
        std::string w = c.text;
        w.erase(0, w.find_first_not_of(' '));
        e.flags.push_back(w);
      }
    }
    return e;
  }

  void copies(const Node& c, int& lo, int& hi, double* proportion) {
    // " × 2", "x2~3", "x3[0.5]": the numbers after the sign.
    std::string_view s = c.text;
    std::size_t i = 0;
    while (i < s.size() && !(s[i] >= '0' && s[i] <= '9')) i++;
    std::size_t j = i;
    while (j < s.size() && s[j] >= '0' && s[j] <= '9') j++;
    lo = hi = small_number(s.substr(i, j - i));
    const std::size_t br = s.find('[');
    std::size_t k = j;
    while (k < s.size() && k < br && !(s[k] >= '0' && s[k] <= '9')) k++;
    if (k < s.size() && k < br) hi = small_number(s.substr(k));
    if (proportion) {
      if (const Node* p = c.child("proportion")) {
        double v = 0, scale = 1;
        bool point = false;
        for (char ch : p->text) {
          if (ch == '.') point = true;
          else if (ch >= '0' && ch <= '9') { if (point) { scale /= 10; v += (ch - '0') * scale; } else v = v * 10 + (ch - '0'); }
        }
        *proportion = v;
        if (v > 1) error(p->begin, p->end, "a proportion of cells lies between 0 and 1");
      }
    }
  }

  Junction junction(const Node& n) {
    Junction j;
    j.span = {n.begin, n.end};
    const auto ends = n.all("junction_end");
    for (std::size_t k = 0; k < ends.size() && k < 2; k++) {
      Junction::End& e = k == 0 ? j.a : j.b;
      const Node& en = *ends[k];
      if (const Node* ch = en.child("chromosome")) e.chromosome = chromosome(*ch);
      if (const Node* a = en.child("arm")) e.arm = a->text[0];
      if (const Node* d = en.child("band_digits")) e.digits = digits(*d);
      if (const Node* s = en.child("strand")) e.strand = s->text == "(+)" ? 1 : -1;
      if (const Node* lc = en.child("lf_coord")) {
        const auto cs = lc->all("coord");
        if (!cs.empty()) e.position = coord_value(cs[0]->text);
      }
    }
    if (const Node* b = n.child("lf_bases")) j.bases = b->text;
    return j;
  }

  Technique technique(const Node& n) {
    Technique t;
    t.span = {n.begin, n.end};
    t.text = n.text;
    if (n.name == "arr") t.name = n.text.substr(0, 4) == "sseq" ? "sseq" : "arr";
    else if (n.name == "seq") t.name = "seq";
    else if (n.name == "rsa") t.name = "rsa";
    else t.name = n.text.substr(0, 3) == "nuc" ? "nuc ish" : "ish";
    if (const Node* b = n.child("build")) t.build = b->text.substr(1, b->text.size() - 2);
    for (const auto& c : n.children) {
      if (c.name == "arr_entry") t.arr.push_back(arr_entry(c));
      else if (c.name == "structural") {
        Abnormality a = abnormality(c);
        validate(a);
        to_edition_2020(a);
        t.structural.push_back(std::move(a));
      } else if (c.name == "junction") {
        t.junctions.push_back(junction(c));
      }
    }
    return t;
  }

  Description description(const Node& root) {
    Description d;
    for (const auto& c : root.children) {
      if (c.name == "karyotype") {
        if (const Node* m = c.child("mosaic")) d.mosaic = m->text;
        std::string sep;
        for (const auto& k : c.children) {
          if (k.name == "line_sep") {
            sep = k.text;
            sep.erase(std::remove(sep.begin(), sep.end(), ' '), sep.end());
          } else if (k.name == "cell_line") {
            CellLine line = cell_line(k);
            line.separator = sep;
            d.lines.push_back(std::move(line));
          }
        }
      } else if (c.name == "arr" || c.name == "seq" || c.name == "ish" || c.name == "rsa") {
        d.techniques.push_back(technique(c));
      }
    }
    // A reference to an earlier clone needs one.
    for (std::size_t i = 0; i < d.lines.size(); i++) {
      const CellLine& l = d.lines[i];
      if (l.reference != CellLine::Reference::none && i == 0)
        error(l.span.begin, l.span.end, "idem, sl and sdl refer to an earlier clone, and this is the first");
      if (l.reference == CellLine::Reference::sdl && l.sdl_number > 0 && static_cast<std::size_t>(l.sdl_number) >= i)
        error(l.span.begin, l.span.end, "sdl" + std::to_string(l.sdl_number) + " names a sideline that is not written before");
      if (i > 0 && l.separator == "//" && d.mosaic == "mos")
        warning(l.span.begin, l.span.end, "// separates the cell lines of a chimera, and the description says mos");
    }
    return d;
  }

 private:
  std::string_view text_;
  Edition edition_;
  std::vector<Diagnostic>& out_;
};

}  // namespace detail

// Parses an ISCN description. The syntax decides acceptance; the diagnostics say what an accepted
// text fails to mean, and the description is filled in either case for what could be read.
inline Parsed parse(std::string_view text, Edition edition = Edition::iscn2020) {
  Parsed out;
  out.syntax = iscn_bnf().parse(text);
  if (!out.syntax.accepted) {
    out.diagnostics.push_back({Diagnostic::Level::error, out.syntax.furthest, std::min(out.syntax.furthest + 1, text.size()), out.syntax.error});
    return out;
  }
  detail::Builder b(text, edition, out.diagnostics);
  out.description = b.description(out.syntax.tree);
  return out;
}

}  // namespace genomics::karyotype
