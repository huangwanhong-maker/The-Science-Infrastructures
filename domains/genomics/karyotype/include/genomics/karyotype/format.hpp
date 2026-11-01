// ISCN descriptions written back in a canonical form: the order of the abnormalities, of the
// chromosomes within a translocation and of the breakpoints within one chromosome as ISCN 2024
// sets them, its spacing, and either its signs (×, →, −) or their ASCII stand-ins (x, ->, -).
//
// The order (ISCN 2024 §4 and ch. 5–6, as summarised by Hastings et al. 2025; ISCN 2020 ch. 6; the
// rules below are those its examples follow, quoted from ISCN 2009 and 2020 by Warrender et al.
// 2019 and Miller and Madan 2021, and from the GenQA assessment of 2021, Hastings, Tabiner and
// Deans):
// 1. sex chromosomes first (X before Y), then the autosomes in ascending order, each abnormality
//    filed under its chromosome (a translocation under its first, a derivative under its own);
//    abnormalities of unknown chromosomes last: der(?), then rings of unknown origin, markers and
//    double minutes, then inc;
// 2. for one chromosome, numerical changes before structural ones, a constitutional change (c)
//    before an acquired one, a gain before a loss;
// 3. several structural abnormalities of one chromosome in the alphabetical order of their
//    abbreviations, one present in the stemline before an additional (+) one;
// 4. within t(a;b) the sex chromosome or the lower-numbered chromosome first; within a
//    translocation of three chromosomes or more the same chromosome first, the order of the others
//    kept (it is the order in which the segments move);
// 5. the breakpoints of an inversion or a deletion from pter to qter (ISCN 2020 §4.3.1.1).
// The microarray's entries: sex chromosomes first, then autosomes, each from pter to qter
// (ISCN 2020 §14.2.2, reversing ISCN 2016).
//
// Spacing (ISCN 2024): a space after the technique and its build (arr[GRCh38] 1p36.33…), around
// mos, chi and or, between two abbreviations not separated otherwise (x1 dn, mar c); none between a
// parenthesis and the abbreviation that follows (…(q34;q11.2)dn).
//
// Complexity: O(k log k) for k abnormalities.

#pragma once

#include <genomics/karyotype/iscn.hpp>

#include <algorithm>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

namespace genomics::karyotype {

struct FormatOptions {
  bool unicode = false;   // ×, →, − instead of x, ->, -
  bool order = true;      // put the abnormalities and the chromosomes of translocations in ISCN's order
  bool commas = true;     // 21,340,001 rather than 21340001 in coordinates
};

// The rank of a chromosome in ISCN's order: X, Y, 1, …, 22; unknown last.
inline int iscn_rank(int c) {
  if (c == chr_x) return 0;
  if (c == chr_y) return 1;
  if (c < 0) return 100;
  return c + 2;
}

inline std::string with_commas(std::uint64_t v) {
  std::string s = std::to_string(v);
  for (int i = static_cast<int>(s.size()) - 3; i > 0; i -= 3) s.insert(static_cast<std::size_t>(i), ",");
  return s;
}

// The order of two breakpoints of one chromosome from pter to qter: −1, +1, or 0 when it cannot be
// told from the names alone (nested bands, unknowns).
inline int compare_by_name(const Breakpoint& a, const Breakpoint& b) {
  const auto key = [](const Breakpoint& x) -> int {
    switch (x.kind) {
      case Breakpoint::Kind::pter: return 0;
      case Breakpoint::Kind::cen: return 2;
      case Breakpoint::Kind::qter: return 4;
      case Breakpoint::Kind::unknown: return -1;
      case Breakpoint::Kind::band: return x.arm == 'p' ? 1 : 3;
    }
    return -1;
  };
  const int ka = key(a), kb = key(b);
  if (ka < 0 || kb < 0) return 0;
  if (ka != kb) return ka < kb ? -1 : 1;
  if (ka != 1 && ka != 3) return 0;
  if (a.centromeric() || b.centromeric()) {
    if (a.centromeric() && b.centromeric()) return 0;
    // p10 is the arm's centromeric end: after every band of p, before every band of q.
    if (ka == 1) return a.centromeric() ? 1 : -1;
    return a.centromeric() ? -1 : 1;
  }
  const std::string& x = a.digits;
  const std::string& y = b.digits;
  const std::size_t n = std::min(x.size(), y.size());
  const int c = x.compare(0, n, y, 0, n);
  if (c == 0) return 0;  // one band holds the other
  // Bands are numbered outward from the centromere: on p the higher comes first.
  return (ka == 1) == (c > 0) ? -1 : 1;
}

namespace detail {

inline std::string times(const FormatOptions& o) { return o.unicode ? "\xC3\x97" : "x"; }
inline std::string minus(const FormatOptions& o) { return o.unicode ? "\xE2\x88\x92" : "-"; }
inline std::string arrow(const FormatOptions& o) { return o.unicode ? "\xE2\x86\x92" : "->"; }

inline std::string format_breakpoint(const Breakpoint& b, const FormatOptions& o) {
  std::string s;
  s = b.label();
  if (b.has_position) {
    s += "(" + (o.commas ? with_commas(b.position) : std::to_string(b.position));
    if (b.position_to) s += "~" + (o.commas ? with_commas(b.position_to) : std::to_string(b.position_to));
    else if (b.position_open) s += "~";
    s += ")";
  }
  return s;
}

// Puts a two-way translocation's chromosomes in order, and rotates one of more chromosomes so
// that its first is the lowest; the groups follow their chromosomes.
inline void order_translocation(Abnormality& a) {
  if (a.type != Type::t || a.chromosomes.size() < 2) return;
  for (const auto& c : a.chromosomes)
    if (c.unknown) return;  // t(?;3): the order of what is not known is kept
  const std::size_t n = a.chromosomes.size();
  std::size_t best = 0;
  for (std::size_t i = 1; i < n; i++)
    if (iscn_rank(a.chromosomes[i].index) < iscn_rank(a.chromosomes[best].index) && !a.chromosomes[i].unknown) best = i;
  if (a.chromosomes[best].unknown) return;
  if (n == 2) {
    if (best == 1 && iscn_rank(a.chromosomes[1].index) != iscn_rank(a.chromosomes[0].index)) {
      std::swap(a.chromosomes[0], a.chromosomes[1]);
      if (a.groups.size() == 2) std::swap(a.groups[0], a.groups[1]);
      if (a.long_form.size() == 2) std::swap(a.long_form[0], a.long_form[1]);
    }
    return;
  }
  if (best == 0) return;
  std::rotate(a.chromosomes.begin(), a.chromosomes.begin() + static_cast<std::ptrdiff_t>(best), a.chromosomes.end());
  if (a.groups.size() == n) std::rotate(a.groups.begin(), a.groups.begin() + static_cast<std::ptrdiff_t>(best), a.groups.end());
  if (a.long_form.size() == n) std::rotate(a.long_form.begin(), a.long_form.begin() + static_cast<std::ptrdiff_t>(best), a.long_form.end());
}

inline void order_breakpoints(Abnormality& a) {
  if ((a.type == Type::inv || a.type == Type::del) && a.groups.size() == 1 && a.groups[0].size() == 2)
    if (compare_by_name(a.groups[0][0], a.groups[0][1]) > 0) std::swap(a.groups[0][0], a.groups[0][1]);
}

inline void canonical(Abnormality& a, const FormatOptions& o) {
  if (!o.order) return;
  order_translocation(a);
  order_breakpoints(a);
  for (auto& c : a.components) canonical(c, o);
  for (auto& c : a.alternatives) canonical(c, o);
}

inline std::string format_long_point(const Breakpoint& b, bool with_chromosome, const FormatOptions& o) {
  std::string s;
  if (with_chromosome && b.chromosome >= 0) s += chromosome_name(b.chromosome);
  return s + format_breakpoint(b, o);
}

inline std::string format_long(const std::vector<LongChromosome>& lf, const FormatOptions& o) {
  std::string out;
  for (std::size_t k = 0; k < lf.size(); k++) {
    if (k) out += ";";
    const LongChromosome& c = lf[k];
    out += c.open;
    for (std::size_t i = 0; i < c.parts.size(); i++) {
      if (i) out += "::";
      const LongPart& p = c.parts[i];
      switch (p.kind) {
        case LongPart::Kind::segment:
          out += format_long_point(p.from, true, o) + arrow(o) + format_long_point(p.to, true, o);
          break;
        case LongPart::Kind::unknown: out += "?"; break;
        case LongPart::Kind::bases: out += p.bases; break;
        case LongPart::Kind::region:
          out += (p.from.chromosome >= 0 ? chromosome_name(p.from.chromosome) : "") + band_label(p.from.arm, p.from.digits) +
                 band_label(p.to.arm, p.to.digits) + "(" + (o.commas ? with_commas(p.from.position) : std::to_string(p.from.position)) + "-" +
                 (o.commas ? with_commas(p.to.position) : std::to_string(p.to.position)) + ")" + times(o) + std::to_string(p.multiplier);
          break;
      }
    }
    out += c.close;
  }
  return out;
}

inline std::string format_head(const Abnormality& a, const FormatOptions& o) {
  std::string s(type_name(a.type));
  s += "(";
  for (std::size_t i = 0; i < a.chromosomes.size(); i++) s += (i ? ";" : "") + a.chromosomes[i].name();
  s += ")";
  if (!a.long_form.empty()) {
    s += "(" + format_long(a.long_form, o) + ")";
  } else if (a.has_groups) {
    s += "(";
    for (std::size_t g = 0; g < a.groups.size(); g++) {
      if (g) s += ";";
      for (const auto& b : a.groups[g]) s += format_breakpoint(b, o);
    }
    s += ")";
  }
  return s;
}

inline std::string format_flags(const Abnormality& a, bool after_letter) {
  std::string s;
  if (a.constitutional) s += (after_letter ? " c" : "c");
  if (!a.inheritance.empty()) s += ((after_letter || a.constitutional) ? " " : "") + a.inheritance;
  return s;
}

}  // namespace detail

// One abnormality in the short system (or, when it was written so, the detailed one).
inline std::string format(const Abnormality& in, const FormatOptions& o = {}) {
  Abnormality a = in;
  detail::canonical(a, o);
  std::string s;
  if (a.sign > 0) s += "+";
  if (a.sign < 0) s += detail::minus(o);
  switch (a.type) {
    case Type::gain:
    case Type::loss:
      s += a.chromosomes.empty() ? "?" : a.chromosomes[0].name();
      s += detail::format_flags(a, false);
      break;
    case Type::marker:
      if (a.count > 1) s += std::to_string(a.count);
      s += "mar";
      if (a.marker_number > 0) s += std::to_string(a.marker_number);
      if (a.multiplier > 1) s += detail::times(o) + std::to_string(a.multiplier);
      s += detail::format_flags(a, true);
      break;
    case Type::ring_marker:
      if (a.count > 1) s += std::to_string(a.count);
      s += "r";
      if (a.multiplier > 1) s += detail::times(o) + std::to_string(a.multiplier);
      s += detail::format_flags(a, true);
      break;
    case Type::dmin:
      if (a.count_hi > a.count) s += std::to_string(a.count) + "~" + std::to_string(a.count_hi);
      else if (a.count > 1) s += std::to_string(a.count);
      s += "dmin";
      break;
    default: {
      if (a.uncertain) s += "?";
      for (const auto& p : a.prefixes) s += p + " ";
      s += detail::format_head(a, o);
      if (a.type == Type::rec && a.rec_chromosome >= 0)
        s += std::string(a.rec_dup ? "dup(" : "del(") + chromosome_name(a.rec_chromosome) + std::string(1, a.rec_arm) + ")";
      for (const auto& c : a.components) s += detail::format_head(c, o);
      if (a.multiplier > 1) s += detail::times(o) + std::to_string(a.multiplier);
      s += detail::format_flags(a, false);
      break;
    }
  }
  for (const auto& alt : a.alternatives) s += " or " + format(alt, o);
  return s;
}

namespace detail {

// The place of an abnormality in ISCN's order (rule 1–3 above).
inline std::tuple<int, int, int, std::string, int> order_key(const Abnormality& a) {
  int group = 0, kind = 0, sub = 0;
  switch (a.type) {
    case Type::ring_marker: group = 201; break;
    case Type::marker: group = 202; break;
    case Type::dmin: group = 203; break;
    default: group = iscn_rank(a.chromosome()); break;
  }
  if (a.numerical()) {
    kind = 0;
    sub = (a.constitutional ? 0 : 2) + (a.sign > 0 ? 0 : 1);
  } else {
    kind = 1;
    sub = a.sign > 0 ? 1 : 0;
  }
  std::string abbr = a.numerical() ? std::string() : std::string(type_name(a.type));
  return {group, kind, 0, abbr, sub};
}

}  // namespace detail

// A cell line: count, sex chromosomes, abnormalities in ISCN's order, cells.
inline std::string format(const CellLine& line, const FormatOptions& o = {}) {
  std::string s;
  if (line.count_unknown) s = "?";
  else {
    s = std::to_string(line.count);
    if (line.count_hi > line.count) s += "~" + std::to_string(line.count_hi);
  }
  if (line.ploidy > 0) s += "<" + std::to_string(line.ploidy) + "n>";
  if (line.reference == CellLine::Reference::idem) s += ",idem";
  else if (line.reference == CellLine::Reference::sl) s += ",sl";
  else if (line.reference == CellLine::Reference::sdl) s += ",sdl" + (line.sdl_number > 0 ? std::to_string(line.sdl_number) : std::string());
  else if (line.has_sex) {
    if (line.undisclosed) s += ",U";
    else {
      s += "," + line.sex;
      if (line.sex_uncertain) s += "?";
      if (line.sex_constitutional) s += "c";
    }
  }
  std::vector<const Abnormality*> list;
  for (const auto& a : line.abnormalities) list.push_back(&a);
  if (o.order) {
    // Order by the canonical forms, so that the chromosome of a translocation is its first after
    // the chromosomes are put in order.
    std::vector<std::pair<std::tuple<int, int, int, std::string, int>, const Abnormality*>> keyed;
    for (const Abnormality* a : list) {
      Abnormality c = *a;
      detail::canonical(c, o);
      keyed.push_back({detail::order_key(c), a});
    }
    std::stable_sort(keyed.begin(), keyed.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
    list.clear();
    for (const auto& k : keyed) list.push_back(k.second);
  }
  for (const Abnormality* a : list) s += "," + format(*a, o);
  if (line.incomplete) s += ",inc";
  if (line.cells >= 0) s += "[" + std::string(line.composite ? "cp" : "") + std::to_string(line.cells) + "]";
  return s;
}

inline std::string format(const ArrEntry& e, const FormatOptions& o = {}) {
  std::string s;
  const auto num = [&](std::uint64_t v) { return o.commas ? with_commas(v) : std::to_string(v); };
  const auto cps = [&](int lo, int hi) {
    std::string t = detail::times(o) + std::to_string(lo);
    if (hi > lo) t += "~" + std::to_string(hi);
    return t;
  };
  switch (e.kind) {
    case ArrEntry::Kind::region:
      s = chromosome_name(e.chromosome) + e.band_name() + "(";
      if (e.extended)
        s += std::to_string(e.flank_start) + detail::times(o) + std::to_string(e.flank_start_copies) + "," + std::to_string(e.start) + "_" +
             std::to_string(e.end) + detail::times(o) + std::to_string(e.copies) + "," + std::to_string(e.flank_end) + detail::times(o) +
             std::to_string(e.flank_end_copies) + ")";
      else {
        s += num(e.start) + "_" + num(e.end) + ")";
        if (e.copies >= 0) s += cps(e.copies, e.copies_hi);
      }
      break;
    case ArrEntry::Kind::whole: {
      s = "(";
      for (std::size_t i = 0; i < e.wholes.size(); i++) {
        if (i) s += ",";
        s += chromosome_name(e.wholes[i].first);
        if (e.wholes[i].second != e.wholes[i].first) s += "-" + chromosome_name(e.wholes[i].second);
      }
      s += ")" + cps(e.copies, e.copies_hi);
      break;
    }
    case ArrEntry::Kind::complex:
      s = "(";
      for (std::size_t i = 0; i < e.complex_chromosomes.size(); i++) s += (i ? "," : "") + chromosome_name(e.complex_chromosomes[i]);
      s += ")cx";
      break;
  }
  if (e.proportion >= 0) {
    std::string p = std::to_string(e.proportion);
    while (p.size() > 1 && p.back() == '0') p.pop_back();
    if (!p.empty() && p.back() == '.') p.pop_back();
    s += "[" + p + "]";
  }
  for (const auto& f : e.flags) s += (s.back() == ')' ? "" : " ") + f;
  return s;
}

// The order of microarray entries: sex chromosomes first, then autosomes, each from pter to qter;
// whole chromosomes before regions of the same chromosome.
inline bool arr_before(const ArrEntry& a, const ArrEntry& b) {
  const auto first_chrom = [](const ArrEntry& e) {
    if (e.kind == ArrEntry::Kind::whole && !e.wholes.empty()) return e.wholes[0].first;
    if (e.kind == ArrEntry::Kind::complex && !e.complex_chromosomes.empty()) return e.complex_chromosomes[0];
    return e.chromosome;
  };
  const int ra = iscn_rank(first_chrom(a)), rb = iscn_rank(first_chrom(b));
  if (ra != rb) return ra < rb;
  if ((a.kind == ArrEntry::Kind::whole) != (b.kind == ArrEntry::Kind::whole)) return a.kind == ArrEntry::Kind::whole;
  return a.start < b.start;
}

inline std::string format(const Technique& t, const FormatOptions& o = {}) {
  if (t.name == "ish" || t.name == "nuc ish" || t.name == "rsa") return t.text;
  std::string s = t.name;
  if (!t.build.empty()) s += "[" + t.build + "]";
  s += " ";
  std::vector<std::string> parts;
  if (t.name == "arr" || t.name == "sseq") {
    std::vector<ArrEntry> entries = t.arr;
    if (o.order) std::stable_sort(entries.begin(), entries.end(), arr_before);
    for (const auto& e : entries) parts.push_back(format(e, o));
  } else {
    for (const auto& a : t.structural) parts.push_back(format(a, o));
    for (const auto& e : t.arr) parts.push_back(format(e, o));
    for (const auto& j : t.junctions) {
      const auto end = [&](const Junction::End& e) {
        std::string x = chromosome_name(e.chromosome) + band_label(e.arm, e.digits);
        if (e.strand) x += e.strand > 0 ? "(+)" : "(" + detail::minus(o) + ")";
        return x + "(" + (o.commas ? with_commas(e.position) : std::to_string(e.position)) + ")";
      };
      parts.push_back(end(j.a) + "::" + (j.bases.empty() ? "" : j.bases + "::") + end(j.b));
    }
  }
  for (std::size_t i = 0; i < parts.size(); i++) s += (i ? "," : "") + parts[i];
  return s;
}

inline std::string format(const Description& d, const FormatOptions& o = {}) {
  std::string s;
  if (!d.mosaic.empty()) s += d.mosaic + " ";
  for (std::size_t i = 0; i < d.lines.size(); i++) {
    if (i) s += d.lines[i].separator.empty() ? "/" : d.lines[i].separator;
    s += format(d.lines[i], o);
  }
  for (std::size_t i = 0; i < d.techniques.size(); i++) {
    if (!s.empty()) s += ".";
    s += format(d.techniques[i], o);
  }
  return s;
}

}  // namespace genomics::karyotype
