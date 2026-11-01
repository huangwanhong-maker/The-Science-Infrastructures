// Structural variants as breakends, adjacencies and copy-number segments: read from and written to
// VCF's symbolic and breakend records (VCF 4.3, §1.4.5 symbolic alleles, §3 the INFO keys for
// structural variants, §5.3 "Encoding Structural Variants", §5.4 "Specifying complex
// rearrangements with breakends"; the same sections of VCF 4.2), and described in HGVS
// (DNA/deletion.md, duplication.md, inversion.md, insertion.md, uncertain.md, complex.md).
//
// - A breakend is a base of a contig and the side of it a junction is on: strand '+' when the
//   sequence kept runs up to the base from the left and the junction follows it, '-' when the
//   sequence kept runs on to the right of the base and the junction precedes it. An adjacency joins
//   two breakends, with any bases inserted between them, read from the first breakend's kept
//   sequence across the junction into the second's; the same adjacency read the other way has its
//   breakends swapped and the inserted bases reverse-complemented. Its canonical form puts the
//   smaller breakend first.
// - VCF's four bracket forms (§5.4) are these strands: for a record at POS with REF s,
//     t[p[  POS '+', mate p '-'        t]p]  POS '+', mate p '+'
//     ]p]t  POS '-', mate p '+'        [p[t  POS '-', mate p '-'
//   with t = s·x for the first two (x the inserted bases, as read from POS to the mate) and
//   t = x′·s for the last two (x′ the reverse complement of the bases read from POS to the mate).
//   Single breakends are "s." and ".s" (with inserted bases "TCC.", ".TGCA"); telomeric breakends
//   stand at POS 0 and N + 1 (§5.4.4), which variant/vcf.hpp reads only from POS 1 on, so they are
//   read here.
// - Symbolic records (§5.3): <DEL> deletes POS+1..END (POS is the base before the event, required
//   with symbolic alleles, §1.6.1), so it is the adjacency (POS,'+')–(END+1,'-'); <DUP> and
//   <DUP:TANDEM> repeat POS+1..END, the adjacency (END,'+')–(POS+1,'-') when the copy is in tandem;
//   <INV> inverts POS+1..END, two adjacencies (POS,'+')–(END,'+') and (POS+1,'-')–(END+1,'-'), the
//   spec's own INV0 written as four breakends; <INS> inserts SVLEN bases after POS; <CNV> is a
//   segment of its copy number (INFO or FORMAT CN). CIPOS and CIEND bound POS and END; IMPRECISE,
//   MATEID, EVENT, SVLEN are kept. SVLEN is signed in VCF 4.2 and 4.3 (negative for deletions) and
//   unsigned from 4.4; its absolute value is read.
// - HGVS descriptions: a deletion g.x_ydel, imprecise g.(a_b)_(c_d)del with B_C the minimal and A_D
//   the maximal extent (uncertain.md); a tandem duplication dup; a duplication with no evidence of
//   tandem position, as HGVS asks, an insertion of a copy at an unknown place, g.?_?ins[ref:g.x_y],
//   or, with its copy number, g.x_y[n] (duplication.md, "three copies … were detected"); an
//   inversion inv; an insertion insN[n]; a copy-number loss del; a translocation as the derivative
//   chromosome's delins, the part replaced and the part inserted, named by the chromosome whose
//   centromere the derivative keeps (complex.md: NC_000002.12:g.pter_8247756delins[NC_000011.10:g.pter_15825266]),
//   an inverted piece written with its ends 5′ to 3′ and inv (the standard writes both orders in one
//   example; theory/nomenclature.md §6).
//
// Coordinates are 0-based inside (a breakend's base), 1-based in VCF and HGVS.

#pragma once

#include <genomics/nomenclature/hgvs.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/variant/vcf.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace genomics::nomenclature {

struct Breakend {
  std::string contig;
  std::int64_t position = 0;  // 0-based; −1 and N for VCF's telomeric breakends at POS 0 and N + 1
  char strand = '+';          // '+': the junction follows the base; '-': it precedes it
  friend bool operator==(const Breakend&, const Breakend&) = default;
  friend auto operator<=>(const Breakend& a, const Breakend& b) { return std::tie(a.contig, a.position, a.strand) <=> std::tie(b.contig, b.position, b.strand); }
};

struct Adjacency {
  Breakend a, b;
  std::string inserted;       // read from a across the junction into b
  bool single = false;        // a single breakend: b unknown
  friend bool operator==(const Adjacency&, const Adjacency&) = default;
};

// The same adjacency read from b to a.
inline Adjacency reversed(const Adjacency& x) {
  if (x.single) return x;
  return Adjacency{x.b, x.a, sequence::reverse_complement(x.inserted), false};
}
inline Adjacency canonical(const Adjacency& x) { return !x.single && x.b < x.a ? reversed(x) : x; }

struct CopyNumberSegment {
  std::string contig;
  std::int64_t begin = 0, end = 0;  // 0-based, [begin, end)
  int copy_number = 2;
  int ploidy = 2;
};

enum class SvKind : std::uint8_t { deletion, duplication, tandem_duplication, inversion, insertion, copy_number, breakend, single_breakend, other };

inline const char* sv_kind_name(SvKind k) {
  static constexpr const char* names[] = {"deletion", "duplication", "tandem duplication", "inversion", "insertion", "copy-number variant", "breakend", "single breakend", "other"};
  return names[static_cast<std::size_t>(k)];
}

// A structural-variant record of a VCF.
struct StructuralVariant {
  SvKind kind = SvKind::other;
  std::string id, chrom, alt, subtype;  // subtype: ME:ALU of <DEL:ME:ALU>
  std::int64_t pos = 0, end = 0;         // VCF's POS and END, 1-based
  std::int64_t length = 0;               // |SVLEN|
  std::array<std::int64_t, 2> cipos{0, 0}, ciend{0, 0};
  bool imprecise = false;
  std::optional<int> copy_number;        // INFO CN, else the first sample's FORMAT CN
  std::vector<std::string> mates;        // MATEID
  std::string event;                     // EVENT
  Adjacency adjacency;                   // breakends: as the record reads it, from POS
  std::string ref;
  std::size_t line = 0;
};

// ---- the bracket notation ---------------------------------------------------------------------------------
struct BreakendAlt {
  bool ok = false;
  Adjacency adjacency;  // from the record's own breakend
  std::string error;
};

// Reads a breakend replacement (§5.4) of a record at `chrom`:`pos` (1-based) with REF `ref`.
inline BreakendAlt parse_breakend(std::string_view alt, std::string_view chrom, std::int64_t pos, std::string_view ref) {
  BreakendAlt out;
  Adjacency& adj = out.adjacency;
  adj.a.contig = std::string(chrom);
  adj.a.position = pos - 1;
  const std::string up = [&] { std::string s(alt); for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c))); return s; }();
  const std::string s(ref.empty() ? "N" : std::string(1, static_cast<char>(std::toupper(static_cast<unsigned char>(ref[0])))));
  // Single breakends: "s." or ".s", with inserted bases.
  if (!up.empty() && up.find('[') == std::string::npos && up.find(']') == std::string::npos) {
    if (up.back() == '.' && up.size() >= 2) {
      adj.single = true;
      adj.a.strand = '+';
      adj.inserted = up.substr(1, up.size() - 2);
      out.ok = up.substr(0, 1) == s || s == "N";
      if (!out.ok) out.error = "a single breakend's first base is not REF's";
      return out;
    }
    if (up.front() == '.' && up.size() >= 2) {
      adj.single = true;
      adj.a.strand = '-';
      adj.inserted = sequence::reverse_complement(up.substr(1, up.size() - 2));
      out.ok = up.substr(up.size() - 1) == s || s == "N";
      if (!out.ok) out.error = "a single breakend's last base is not REF's";
      return out;
    }
    out.error = "not a breakend replacement";
    return out;
  }
  const std::size_t open = up.find_first_of("[]");
  const std::size_t close = up.find_first_of("[]", open + 1);
  if (open == std::string::npos || close == std::string::npos || up[open] != up[close]) { out.error = "a breakend's mate is written between two brackets of one kind"; return out; }
  const char bracket = up[open];
  const std::string mate(alt.substr(open + 1, close - open - 1));
  const std::size_t colon = mate.rfind(':');
  if (colon == std::string::npos) { out.error = "the mate is written chrom:pos"; return out; }
  std::int64_t mpos = 0;
  const std::string_view ps = std::string_view(mate).substr(colon + 1);
  if (std::from_chars(ps.data(), ps.data() + ps.size(), mpos).ptr != ps.data() + ps.size()) { out.error = "the mate's position is not a number"; return out; }
  adj.b.contig = mate.substr(0, colon);
  if (!adj.b.contig.empty() && adj.b.contig.front() == '<' && adj.b.contig.back() == '>') adj.b.contig = adj.b.contig.substr(1, adj.b.contig.size() - 2);
  adj.b.position = mpos - 1;
  const std::string before = up.substr(0, open), after = up.substr(close + 1);
  if (!before.empty() && after.empty()) {
    // t[p[ or t]p]: the junction after POS.
    adj.a.strand = '+';
    adj.b.strand = bracket == '[' ? '-' : '+';
    if (before.substr(0, 1) != s && s != "N" && before != ".") { out.error = "t does not begin with REF's base"; return out; }
    adj.inserted = before == "." ? std::string() : before.substr(1);
  } else if (before.empty() && !after.empty()) {
    // ]p]t or [p[t: the junction before POS.
    adj.a.strand = '-';
    adj.b.strand = bracket == ']' ? '+' : '-';
    if (after.substr(after.size() - 1) != s && s != "N" && after != ".") { out.error = "t does not end with REF's base"; return out; }
    adj.inserted = after == "." ? std::string() : sequence::reverse_complement(after.substr(0, after.size() - 1));
  } else {
    out.error = "a breakend replacement has t on one side of its mate";
    return out;
  }
  out.ok = true;
  return out;
}

// Writes the breakend replacement of an adjacency at its first breakend, REF's base given.
inline std::string breakend_alt(const Adjacency& x, char ref_base) {
  const std::string s = x.a.position < 0 ? std::string(".") : std::string(1, ref_base);  // a telomere at POS 0 (§5.4.4)
  if (x.single) return x.a.strand == '+' ? s + x.inserted + "." : "." + sequence::reverse_complement(x.inserted) + s;
  const std::string p = x.b.contig + ":" + std::to_string(x.b.position + 1);
  const char br = x.a.strand == '+' ? (x.b.strand == '-' ? '[' : ']') : (x.b.strand == '+' ? ']' : '[');
  const std::string mate = std::string(1, br) + p + std::string(1, br);
  return x.a.strand == '+' ? s + x.inserted + mate : mate + sequence::reverse_complement(x.inserted) + s;
}

// ---- reading -------------------------------------------------------------------------------------------------------
namespace detail {

inline std::int64_t to_int(std::string_view s, bool& ok) {
  std::int64_t v = 0;
  ok = !s.empty() && std::from_chars(s.data(), s.data() + s.size(), v).ptr == s.data() + s.size();
  return v;
}
inline std::array<std::int64_t, 2> interval(std::string_view s, bool& ok) {
  ok = false;
  const std::size_t c = s.find(',');
  if (c == std::string_view::npos) return {0, 0};
  bool a = false, b = false;
  const std::int64_t x = to_int(s.substr(0, c), a), y = to_int(s.substr(c + 1), b);
  ok = a && b;
  return {x, y};
}

}  // namespace detail

// The structural-variant records of a VCF text: symbolic alleles and breakends (records of small
// variants are left out). Telomeric breakends at POS 0 are read too. Never throws on content.
inline std::vector<StructuralVariant> read_structural(std::string_view text, std::vector<sequence::Problem>* problems = nullptr) {
  const auto problem = [&](std::size_t line, std::string m) { if (problems) problems->push_back({line, std::move(m)}); };
  // POS 0 (a telomeric breakend) read as POS 1 and noted.
  std::string patched;
  patched.reserve(text.size());
  std::vector<std::size_t> telomeric;
  {
    std::size_t number = 0, i = 0;
    while (i < text.size()) {
      std::size_t j = text.find('\n', i);
      if (j == std::string_view::npos) j = text.size();
      std::string line(text.substr(i, j - i));
      number++;
      if (!line.empty() && line[0] != '#') {
        const std::size_t t1 = line.find('\t'), t2 = t1 == std::string::npos ? std::string::npos : line.find('\t', t1 + 1);
        if (t2 != std::string::npos && line.substr(t1 + 1, t2 - t1 - 1) == "0") {
          line.replace(t1 + 1, 1, "1");
          telomeric.push_back(number);
        }
      }
      patched += line;
      patched += '\n';
      i = j + 1;
    }
  }
  std::vector<sequence::Problem> own;
  const variant::Vcf vcf = variant::read_vcf(patched, &own);
  for (auto& p : own) problem(p.line, p.message);
  std::vector<StructuralVariant> out;
  for (const auto& r : vcf.records) {
    const bool telo = std::find(telomeric.begin(), telomeric.end(), r.line) != telomeric.end();
    for (const auto& alt : r.alts) {
      StructuralVariant sv;
      sv.id = r.id;
      sv.chrom = r.chrom;
      sv.pos = telo ? 0 : static_cast<std::int64_t>(r.position) + 1;
      sv.ref = r.ref;
      sv.alt = alt;
      sv.line = r.line;
      sv.end = sv.pos + static_cast<std::int64_t>(r.ref.size()) - 1;
      bool ok = false;
      const auto info = [&](std::string_view key) { return r.info_value(key); };
      for (const auto& [k, v] : r.info) if (k == "IMPRECISE") sv.imprecise = true;
      if (const auto e = info("END"); !e.empty()) { const auto x = detail::to_int(e, ok); if (ok) sv.end = x; else problem(r.line, "END is not an integer"); }
      if (const auto e = info("SVLEN"); !e.empty()) {
        const auto first = e.substr(0, e.find(','));
        const auto x = detail::to_int(first, ok);
        if (ok) sv.length = std::llabs(x); else problem(r.line, "SVLEN is not an integer");
      }
      if (const auto e = info("CIPOS"); !e.empty()) { sv.cipos = detail::interval(e, ok); if (!ok) problem(r.line, "CIPOS is not two integers"); }
      if (const auto e = info("CIEND"); !e.empty()) { sv.ciend = detail::interval(e, ok); if (!ok) problem(r.line, "CIEND is not two integers"); }
      if (const auto e = info("CN"); !e.empty()) { const auto x = detail::to_int(e, ok); if (ok) sv.copy_number = static_cast<int>(x); }
      if (!sv.copy_number) {
        const auto cn = r.sample_value(0, "CN");
        if (!cn.empty()) { const auto x = detail::to_int(cn, ok); if (ok) sv.copy_number = static_cast<int>(x); }
      }
      if (const auto e = info("MATEID"); !e.empty())
        for (const auto m : variant::detail::fields(e, ',')) sv.mates.emplace_back(m);
      sv.event = std::string(info("EVENT"));
      if (alt.size() > 2 && alt.front() == '<' && alt.back() == '>') {
        const std::string id = alt.substr(1, alt.size() - 2);
        const std::string type = id.substr(0, id.find(':'));
        sv.subtype = id.find(':') == std::string::npos ? "" : id.substr(id.find(':') + 1);
        if (type == "DEL") sv.kind = SvKind::deletion;
        else if (type == "DUP") sv.kind = sv.subtype == "TANDEM" ? SvKind::tandem_duplication : SvKind::duplication;
        else if (type == "INV") sv.kind = SvKind::inversion;
        else if (type == "INS") sv.kind = SvKind::insertion;
        else if (type == "CNV") sv.kind = SvKind::copy_number;
        else sv.kind = SvKind::other;
        if (sv.kind == SvKind::insertion && info("END").empty()) sv.end = sv.pos;
        if (sv.kind != SvKind::insertion && sv.kind != SvKind::other && info("END").empty() && sv.length > 0) sv.end = sv.pos + sv.length;
        if (sv.kind != SvKind::insertion && sv.kind != SvKind::other && sv.length == 0) sv.length = sv.end - sv.pos;
        out.push_back(std::move(sv));
        continue;
      }
      if (alt.find_first_of("[]") != std::string::npos || (alt.size() >= 2 && (alt.front() == '.' || alt.back() == '.'))) {
        const BreakendAlt b = parse_breakend(alt, r.chrom, sv.pos, r.ref);
        if (!b.ok) { problem(r.line, "a breakend that cannot be read (" + alt + "): " + b.error); continue; }
        sv.kind = b.adjacency.single ? SvKind::single_breakend : SvKind::breakend;
        sv.adjacency = b.adjacency;
        out.push_back(std::move(sv));
      }
    }
  }
  return out;
}

// The adjacencies the records assert, each once, in canonical form; a breakend and its mate assert one.
inline std::vector<Adjacency> adjacencies(const std::vector<StructuralVariant>& svs) {
  std::vector<Adjacency> out;
  const auto add = [&](const Adjacency& x) {
    const Adjacency c = canonical(x);
    if (std::find(out.begin(), out.end(), c) == out.end()) out.push_back(c);
  };
  for (const auto& sv : svs) {
    const std::string& c = sv.chrom;
    switch (sv.kind) {
      case SvKind::deletion: add({{c, sv.pos - 1, '+'}, {c, sv.end, '-'}, "", false}); break;
      case SvKind::tandem_duplication: add({{c, sv.end - 1, '+'}, {c, sv.pos, '-'}, "", false}); break;
      case SvKind::inversion:
        add({{c, sv.pos - 1, '+'}, {c, sv.end - 1, '+'}, "", false});
        add({{c, sv.pos, '-'}, {c, sv.end, '-'}, "", false});
        break;
      case SvKind::insertion: add({{c, sv.pos - 1, '+'}, {c, sv.pos, '-'}, std::string(static_cast<std::size_t>(std::max<std::int64_t>(sv.length, 0)), 'N'), false}); break;
      case SvKind::breakend:
      case SvKind::single_breakend: add(sv.adjacency); break;
      default: break;
    }
  }
  return out;
}

// The copy-number segments the records assert (with their copy number when known).
inline std::vector<CopyNumberSegment> copy_number_segments(const std::vector<StructuralVariant>& svs, int ploidy = 2) {
  std::vector<CopyNumberSegment> out;
  for (const auto& sv : svs) {
    if (sv.kind != SvKind::copy_number && sv.kind != SvKind::deletion && sv.kind != SvKind::duplication && sv.kind != SvKind::tandem_duplication) continue;
    if (!sv.copy_number && sv.kind == SvKind::copy_number) continue;
    CopyNumberSegment s;
    s.contig = sv.chrom;
    s.begin = sv.pos;  // POS+1 (1-based) as 0-based
    s.end = sv.end;
    s.ploidy = ploidy;
    s.copy_number = sv.copy_number ? *sv.copy_number : sv.kind == SvKind::deletion ? ploidy - 1 : ploidy + 1;
    out.push_back(s);
  }
  return out;
}

// ---- writing ---------------------------------------------------------------------------------------------------------
// REF's base at a breakend; 'N' when the reference is not given.
using BaseAt = std::function<char(const std::string& contig, std::int64_t position)>;

// The two breakend records of an adjacency (one for a single breakend), with MATEID; ids `id`_1 and `id`_2.
inline std::vector<variant::VcfRecord> breakend_records(const Adjacency& x, const std::string& id, const BaseAt& base = {}) {
  const auto ref_at = [&](const Breakend& b) { return base ? base(b.contig, b.position) : 'N'; };
  std::vector<variant::VcfRecord> out;
  const auto record = [&](const Adjacency& y, const std::string& own, const std::string& mate) {
    variant::VcfRecord r;
    r.chrom = y.a.contig;
    r.position = static_cast<std::size_t>(std::max<std::int64_t>(y.a.position, 0));
    r.id = own;
    const char s = y.a.position < 0 ? 'N' : ref_at(y.a);
    r.ref = std::string(1, s);
    r.alts = {breakend_alt(y, s)};
    r.filters = {"PASS"};
    r.info.push_back({"SVTYPE", "BND"});
    if (!mate.empty()) r.info.push_back({"MATEID", mate});
    out.push_back(std::move(r));
  };
  if (x.single) {
    record(x, id, "");
    return out;
  }
  record(x, id + "_1", id + "_2");
  record(reversed(x), id + "_2", id + "_1");
  return out;
}

// A symbolic record (<DEL>, <DUP>, <DUP:TANDEM>, <INV>, <INS>, <CNV>) of a 1-based interval POS+1..END.
inline variant::VcfRecord symbolic_record(SvKind kind, const std::string& contig, std::int64_t pos, std::int64_t end, char ref_base = 'N',
                                          std::optional<int> copy_number = std::nullopt, std::int64_t inserted_length = 0) {
  variant::VcfRecord r;
  r.chrom = contig;
  r.position = static_cast<std::size_t>(pos - 1);
  r.ref = std::string(1, ref_base);
  const char* alt = kind == SvKind::deletion ? "<DEL>" : kind == SvKind::duplication ? "<DUP>" : kind == SvKind::tandem_duplication ? "<DUP:TANDEM>"
                    : kind == SvKind::inversion ? "<INV>" : kind == SvKind::insertion ? "<INS>" : "<CNV>";
  r.alts = {alt};
  r.filters = {"PASS"};
  const char* type = kind == SvKind::deletion ? "DEL" : kind == SvKind::insertion ? "INS" : kind == SvKind::inversion ? "INV" : kind == SvKind::copy_number ? "CNV" : "DUP";
  r.info.push_back({"SVTYPE", type});
  r.info.push_back({"END", std::to_string(kind == SvKind::insertion ? pos : end)});
  const std::int64_t len = kind == SvKind::insertion ? inserted_length : end - pos;
  r.info.push_back({"SVLEN", std::to_string(kind == SvKind::deletion ? -len : len)});
  if (copy_number) r.info.push_back({"CN", std::to_string(*copy_number)});
  return r;
}

// ---- HGVS -------------------------------------------------------------------------------------------------------------
struct SvHgvsOptions {
  std::function<std::string(const std::string& contig)> accession;          // chr2 → NC_000002.12; the contig as written by default
  std::function<std::optional<std::int64_t>(const std::string& contig)> centromere;  // a 0-based base within the centromere, when known
  std::function<std::int64_t(const std::string& contig)> length;             // the contig's length, for pter and qter (unknown: none)
  bool assume_tandem = false;  // describe <DUP> as a tandem duplication
};

namespace detail {

inline std::string acc(const SvHgvsOptions& o, const std::string& contig) { return o.accession ? o.accession(contig) : contig; }

// An end of a range, exact or (a_b) when its confidence interval is not a point.
inline std::string uncertain_end(std::int64_t x, const std::array<std::int64_t, 2>& ci) {
  if (ci[0] == 0 && ci[1] == 0) return std::to_string(x);
  return "(" + std::to_string(x + ci[0]) + "_" + std::to_string(x + ci[1]) + ")";
}

}  // namespace detail

// The HGVS description of a structural-variant record (g. on the contig's accession).
inline std::string hgvs_of(const StructuralVariant& sv, const SvHgvsOptions& o = {}) {
  const std::string ref = detail::acc(o, sv.chrom) + ":g.";
  const std::string a = detail::uncertain_end(sv.pos + 1, sv.cipos), b = detail::uncertain_end(sv.end, sv.ciend);
  const std::string range = sv.pos + 1 == sv.end && sv.cipos == sv.ciend ? a : a + "_" + b;
  switch (sv.kind) {
    case SvKind::deletion: return ref + range + "del";
    case SvKind::tandem_duplication: return ref + range + "dup";
    case SvKind::duplication:
      if (o.assume_tandem) return ref + range + "dup";
      if (sv.copy_number) return ref + range + "[" + std::to_string(*sv.copy_number) + "]";
      return ref + "?_?ins[" + detail::acc(o, sv.chrom) + ":g." + range + "]";
    case SvKind::inversion: return ref + range + "inv";
    case SvKind::insertion: {
      const std::string len = sv.length > 0 ? "N[" + std::to_string(sv.length) + "]" : std::string("N[?]");
      if (sv.cipos[0] == 0 && sv.cipos[1] == 0) return ref + std::to_string(sv.pos) + "_" + std::to_string(sv.pos + 1) + "ins" + len;
      return ref + "(" + std::to_string(sv.pos + sv.cipos[0]) + "_" + std::to_string(sv.pos + 1 + sv.cipos[1]) + ")ins" + len;
    }
    case SvKind::copy_number:
      if (!sv.copy_number) return ref + range + "[?]";
      if (*sv.copy_number < 2) return ref + range + "del";
      return ref + range + "[" + std::to_string(*sv.copy_number) + "]";
    default: return {};
  }
}

inline std::string hgvs_of(const CopyNumberSegment& s, const SvHgvsOptions& o = {}) {
  const std::string range = detail::acc(o, s.contig) + ":g." + std::to_string(s.begin + 1) + "_" + std::to_string(s.end);
  if (s.copy_number == 0) {
    std::string out = detail::acc(o, s.contig) + ":g.";
    for (int i = 0; i < s.ploidy; i++) out += (i ? ";[" : "[") + std::to_string(s.begin + 1) + "_" + std::to_string(s.end) + "del]";
    return out;
  }
  if (s.copy_number < s.ploidy) return range + "del";
  if (s.copy_number == s.ploidy) return range + "=";
  return range + "[" + std::to_string(s.copy_number) + "]";
}

// The HGVS description of an adjacency: a deletion, a tandem duplication or a delins on one contig
// when its breakends face each other so; otherwise the derivative chromosome's delins, named by the
// contig whose centromere the derivative keeps (the first breakend's when it is not known).
inline std::string hgvs_of(const Adjacency& x0, const SvHgvsOptions& o = {}) {
  if (x0.single) return {};
  Adjacency x = x0;
  // Which breakend's contig names the derivative.
  const auto keeps_centromere = [&](const Breakend& b) -> std::optional<bool> {
    if (!o.centromere) return std::nullopt;
    const auto c = o.centromere(b.contig);
    if (!c) return std::nullopt;
    return b.strand == '+' ? *c <= b.position : *c >= b.position;
  };
  if (x.a.contig == x.b.contig && x.a.strand != x.b.strand) {
    // Facing breakends on one contig.
    const Breakend& left = x.a.strand == '+' ? x.a : x.b;   // the junction follows it
    const Breakend& right = x.a.strand == '+' ? x.b : x.a;  // the junction precedes it
    const std::string ins = x.a.strand == '+' ? x.inserted : sequence::reverse_complement(x.inserted);
    const std::string ref = detail::acc(o, x.a.contig) + ":g.";
    if (right.position > left.position + 1) {
      const std::string range = right.position - left.position - 1 == 1 ? std::to_string(left.position + 2) : std::to_string(left.position + 2) + "_" + std::to_string(right.position);
      return ref + range + (ins.empty() ? "del" : "delins" + ins);
    }
    if (right.position <= left.position && ins.empty()) {
      const std::string range = right.position == left.position ? std::to_string(left.position + 1) : std::to_string(right.position + 1) + "_" + std::to_string(left.position + 1);
      return ref + range + "dup";
    }
    if (right.position == left.position + 1) return ref + std::to_string(left.position + 1) + "_" + std::to_string(right.position + 1) + "ins" + (ins.empty() ? std::string("N[?]") : ins);
  }
  const auto ka = keeps_centromere(x.a), kb = keeps_centromere(x.b);
  if ((ka && !*ka && kb && *kb) || (!ka && kb && *kb)) x = reversed(x);
  // The derivative read from a's kept part: [a's part] inserted [b's part].
  const std::string A = detail::acc(o, x.a.contig), B = detail::acc(o, x.b.contig);
  const auto piece_b = [&](bool forward_of_a) {
    // b's part as read from the junction (forward_of_a: the derivative read in a's forward direction).
    const std::string p1 = std::to_string(x.b.position + 1);
    const bool b_runs_right = x.b.strand == '-';
    const bool inverted = forward_of_a ? !b_runs_right : b_runs_right;
    const std::string range = b_runs_right ? p1 + "_qter" : "pter_" + p1;
    return B + ":g." + range + (inverted ? "inv" : "");
  };
  if (x.a.strand == '+') {
    // A[pter..a] kept; A[a+1..qter] replaced by [inserted; b's part].
    std::string inner = x.inserted.empty() ? piece_b(true) : x.inserted + ";" + piece_b(true);
    return A + ":g." + std::to_string(x.a.position + 2) + "_qterdelins[" + inner + "]";
  }
  // A[a..qter] kept; A[pter..a−1] replaced by [b's part, read toward the junction; inserted].
  const std::string ins = sequence::reverse_complement(x.inserted);
  std::string inner = ins.empty() ? piece_b(false) : piece_b(false) + ";" + ins;
  return A + ":g.pter_" + std::to_string(x.a.position) + "delins[" + inner + "]";
}

// Descriptions of a set of adjacencies: the two junctions of an inversion, (x,'+')–(y,'+') and
// (x+1,'-')–(y+1,'-'), described together as g.(x+2)_(y+1)inv; every other adjacency by itself.
inline std::vector<std::string> hgvs_of(const std::vector<Adjacency>& xs, const SvHgvsOptions& o = {}) {
  std::vector<std::string> out;
  std::vector<char> used(xs.size(), 0);
  for (std::size_t i = 0; i < xs.size(); i++) {
    if (used[i]) continue;
    const Adjacency x = canonical(xs[i]);
    if (!x.single && x.a.contig == x.b.contig && x.a.strand == '+' && x.b.strand == '+' && x.inserted.empty()) {
      for (std::size_t j = 0; j < xs.size(); j++) {
        if (j == i || used[j]) continue;
        const Adjacency y = canonical(xs[j]);
        if (!y.single && y.a.contig == x.a.contig && y.a.strand == '-' && y.b.strand == '-' && y.inserted.empty() && y.a.position == x.a.position + 1 &&
            y.b.position == x.b.position + 1) {
          used[i] = used[j] = 1;
          out.push_back(detail::acc(o, x.a.contig) + ":g." + std::to_string(x.a.position + 2) + "_" + std::to_string(x.b.position + 1) + "inv");
          break;
        }
      }
      if (used[i]) continue;
    }
    used[i] = 1;
    out.push_back(hgvs_of(xs[i], o));
  }
  return out;
}

}  // namespace genomics::nomenclature
