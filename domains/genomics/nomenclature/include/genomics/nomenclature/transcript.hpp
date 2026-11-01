// Transcripts and the numbering of their positions (HGVS Nomenclature 21.1, background/numbering.md),
// and the projection of a genomic change onto a transcript as its c. or n. description.
//
// - A transcript is its exons in the order it reads them, each a segment of a genomic sequence, on
//   one strand; a coding transcript has its CDS, from the A of the start codon to the last base of
//   the stop codon, in transcript coordinates. Models come from GenBank (an mRNA's join and its CDS
//   on a genomic record, both strands; a CDS alone where no mRNA is annotated, as in mitochondrial
//   genomes; or an mRNA record, which is its own transcript) and from GFF3 (an mRNA or transcript
//   with its exon and CDS children).
// - Numbering (numbering.md): c.1 is the A of ATG; c.-1 the base before it, c.*1 the base after the
//   stop; there is no c.0. An intronic base is numbered from the nearer exon, c.87+1 after exon
//   end c.87, c.88-1 before exon start c.88, and in an intron of odd length the central base from
//   the upstream exon with '+' ("the central nucleotide is numbered relative to the upstream exon
//   followed by a +"). n. numbers a non-coding transcript from n.1 with the same offsets. Bases
//   beyond the transcript's ends are numbered on by extrapolation (c.-128354, c.*21568), which the
//   standard allows only with a genomic reference (NC_000023.10(NM_004006.1):c.-128354C>T); they
//   are marked `outside`.
// - The 3′ rule on a transcript is the 3′ rule in the transcript's direction, so a change is placed
//   on the strand the transcript reads (general.md: "the 3'rule applies to ALL descriptions …
//   given variant"; its example, g.32343183del on the plus strand is c.5697del on a minus-strand
//   transcript). The change is described on the reverse complement of the genome for a
//   minus-strand transcript, and the positions mapped back. Placing it on the genome makes the
//   standard's exception at exon–exon junctions (numbering.md: c.3921del, not c.3922del) hold by
//   itself, since the intron stands between the exons; on a transcript-only reference (an NM_
//   record) the exception is applied by not letting the change move past the end of its exon.
// - In a coding sequence a repeated sequence is written only for units of a length divisible by
//   three (DNA/repeated.md, "exception"); elsewhere (UTRs, introns) as on the genome.
//
// Coordinates are 0-based and half-open on the genomic sequence (the record's), 0-based along the
// transcript; HGVS positions are written 1-based. `first` gives the genomic number of the record's
// first base (a region of a chromosome).

#pragma once

#include <genomics/annotation/features.hpp>
#include <genomics/nomenclature/describe.hpp>
#include <genomics/nomenclature/hgvs.hpp>
#include <genomics/sequence/genbank.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::nomenclature {

struct Exon {
  std::size_t begin = 0, end = 0;  // genomic, [begin, end)
  friend bool operator==(const Exon&, const Exon&) = default;
};

struct Transcript {
  std::string id, gene, protein_id, product;
  bool reverse = false;
  std::vector<Exon> exons;   // in the transcript's order, 5′ to 3′
  bool coding = false;
  std::size_t cds_begin = 0, cds_end = 0;  // along the transcript, [cds_begin, cds_end), the stop codon included
  int genetic_code = 1;
  std::string translation;   // the record's /translation, when given
  bool genomic = true;       // false: the record is the transcript itself (an mRNA record)

  std::size_t length() const {
    std::size_t n = 0;
    for (const auto& e : exons) n += e.end - e.begin;
    return n;
  }
  std::size_t lowest() const {
    std::size_t m = SIZE_MAX;
    for (const auto& e : exons) m = std::min(m, e.begin);
    return exons.empty() ? 0 : m;
  }
  std::size_t highest() const {  // one past the highest base
    std::size_t m = 0;
    for (const auto& e : exons) m = std::max(m, e.end);
    return m;
  }
  Coordinate coordinate() const { return coding ? Coordinate::c : Coordinate::n; }
};

// ---- models from records ----------------------------------------------------------------------------------
namespace detail {

inline std::vector<Exon> exons_of(const sequence::Location& l) {
  std::vector<Exon> out;
  for (const auto& p : l.parts) out.push_back({p.begin, p.end});
  return out;
}

// The transcript coordinate of a genomic base lying in an exon, or SIZE_MAX.
inline std::size_t along(const std::vector<Exon>& exons, bool reverse, std::size_t g) {
  std::size_t t = 0;
  for (const auto& e : exons) {
    if (g >= e.begin && g < e.end) return t + (reverse ? e.end - 1 - g : g - e.begin);
    t += e.end - e.begin;
  }
  return SIZE_MAX;
}

// Whether a CDS's segments are those of the transcript between its ends: every junction inside the
// CDS one of the transcript's.
inline bool compatible(const Transcript& t, const sequence::Location& cds) {
  if (cds.parts.empty() || cds.reverse() != t.reverse) return false;
  const std::vector<Exon> c = exons_of(cds);
  const auto first_base = [&](const Exon& e) { return t.reverse ? e.end - 1 : e.begin; };
  const auto last_base = [&](const Exon& e) { return t.reverse ? e.begin : e.end - 1; };
  const std::size_t a = along(t.exons, t.reverse, first_base(c.front())), b = along(t.exons, t.reverse, last_base(c.back()));
  if (a == SIZE_MAX || b == SIZE_MAX || b < a || b - a + 1 != cds.length()) return false;
  for (const auto& e : c)
    if (along(t.exons, t.reverse, first_base(e)) == SIZE_MAX || along(t.exons, t.reverse, last_base(e)) == SIZE_MAX) return false;
  return true;
}

inline std::string variant_number(std::string_view text) {  // "transcript variant 12" → "12"
  const auto at = text.find("transcript variant ");
  if (at == std::string_view::npos) return {};
  std::size_t i = at + 19, j = i;
  while (j < text.size() && text[j] != ';' && text[j] != ',' && text[j] != ' ' && text[j] != '.') j++;
  return std::string(text.substr(i, j - i));
}

inline void set_cds(Transcript& t, const sequence::Feature& cds) {
  const auto& l = cds.location;
  const std::size_t g0 = t.reverse ? l.parts.front().end - 1 : l.parts.front().begin;
  t.coding = true;
  t.cds_begin = along(t.exons, t.reverse, g0);
  t.cds_end = t.cds_begin + l.length();
  t.protein_id = cds.qualifier("protein_id");
  t.translation = cds.qualifier("translation");
  const std::string table = cds.qualifier("transl_table");
  if (!table.empty()) t.genetic_code = std::stoi(table);
  if (t.gene.empty()) t.gene = cds.qualifier("gene");
  if (t.product.empty()) t.product = cds.qualifier("product");
}

}  // namespace detail

// The transcripts of a GenBank record. An mRNA record (molecule mRNA) is one transcript of its own
// sequence; on a genomic record every mRNA or non-coding RNA with a /transcript_id is one, paired
// with its CDS (a CDS fits an mRNA when its segments are the mRNA's between its ends; RefSeq's
// "transcript variant N" in the mRNA's /product and the CDS's /note decides between isoforms), and
// every CDS that fits no mRNA is a transcript of its own (named by /gene, else /locus_tag).
inline std::vector<Transcript> transcripts(const sequence::GenBank& record, std::vector<sequence::Problem>* problems = nullptr) {
  std::vector<Transcript> out;
  const auto problem = [&](std::size_t line, std::string m) { if (problems) problems->push_back({line, std::move(m)}); };
  const auto prefixed = [&](const char* x) { return record.accession.rfind(x, 0) == 0; };
  const bool mrna_record = record.molecule.find("mRNA") != std::string::npos || prefixed("NM_") || prefixed("NR_") || prefixed("XM_") || prefixed("XR_");
  if (mrna_record) {
    Transcript t;
    t.id = record.version.empty() ? record.accession : record.version;
    t.genomic = false;
    for (const auto& f : record.features)
      if (f.key == "exon" && !f.location.parts.empty()) t.exons.push_back({f.location.parts.front().begin, f.location.parts.back().end});
    std::sort(t.exons.begin(), t.exons.end(), [](const Exon& a, const Exon& b) { return a.begin < b.begin; });
    if (t.exons.empty()) t.exons.push_back({0, record.sequence.size()});
    for (const auto& f : record.features) {
      if (f.key == "gene" && t.gene.empty()) t.gene = f.qualifier("gene");
      if (f.key == "CDS" && !t.coding && !f.location.parts.empty() && !f.location.reverse()) detail::set_cds(t, f);
    }
    out.push_back(std::move(t));
    return out;
  }
  static const std::vector<std::string> rna_keys = {"mRNA", "ncRNA", "misc_RNA", "precursor_RNA", "tRNA", "rRNA", "transcript"};
  std::vector<std::size_t> rna_feature;
  for (std::size_t i = 0; i < record.features.size(); i++) {
    const auto& f = record.features[i];
    if (std::find(rna_keys.begin(), rna_keys.end(), f.key) == rna_keys.end() || f.location.parts.empty()) continue;
    // An RNA with no /transcript_id (a mitochondrial tRNA or rRNA) is named by its gene.
    const std::string id = !f.qualifier("transcript_id").empty() ? f.qualifier("transcript_id") : f.key != "mRNA" ? f.qualifier("gene") : std::string();
    if (id.empty()) continue;
    if (f.location.remote) { problem(f.line, f.key + " " + id + ": a part lies in another record; left out"); continue; }
    Transcript t;
    t.id = id;
    t.gene = f.qualifier("gene");
    t.product = f.qualifier("product");
    t.reverse = f.location.reverse();
    t.exons = detail::exons_of(f.location);
    out.push_back(std::move(t));
    rna_feature.push_back(i);
  }
  std::vector<char> used(out.size(), 0);
  for (const auto& f : record.features) {
    if (f.key != "CDS" || f.location.parts.empty() || f.location.remote) continue;
    const std::string variant = detail::variant_number(f.qualifier("note"));
    long best = -1;
    for (std::size_t k = 0; k < out.size(); k++) {
      if (used[k] || record.features[rna_feature[k]].key != "mRNA" || out[k].gene != f.qualifier("gene") || !detail::compatible(out[k], f.location)) continue;
      const bool named = !variant.empty() && detail::variant_number(out[k].product) == variant;
      if (named) { best = static_cast<long>(k); break; }
      if (best < 0) best = static_cast<long>(k);
    }
    if (best >= 0) {
      detail::set_cds(out[static_cast<std::size_t>(best)], f);
      used[static_cast<std::size_t>(best)] = 1;
      continue;
    }
    // A CDS of its own (mitochondrial genes; prokaryotes).
    Transcript t;
    t.gene = f.qualifier("gene");
    t.id = !t.gene.empty() ? t.gene : f.qualifier("locus_tag");
    if (t.id.empty()) t.id = f.qualifier("protein_id");
    t.reverse = f.location.reverse();
    t.exons = detail::exons_of(f.location);
    detail::set_cds(t, f);
    out.push_back(std::move(t));
    used.push_back(1);
  }
  return out;
}

// The transcripts of a GFF3 landmark: entities of type mRNA, transcript or a non-coding RNA, their
// exon children as exons (their CDS children's extent as the CDS).
inline std::vector<Transcript> transcripts(const annotation::Gff& gff, std::string_view seqid) {
  std::vector<Transcript> out;
  static const std::vector<std::string> types = {"mRNA", "transcript", "ncRNA", "lnc_RNA", "lncRNA", "primary_transcript", "snRNA", "snoRNA", "miRNA", "rRNA", "tRNA"};
  for (const auto& e : gff.entities) {
    if (e.seqid != seqid || std::find(types.begin(), types.end(), e.type) == types.end()) continue;
    const auto& head = gff.features[e.lines.front()];
    Transcript t;
    t.id = head.id;
    if (const auto* v = head.attribute("transcript_id")) t.id = v->front();
    else if (!head.name.empty() && head.id.empty()) t.id = head.name;
    if (const auto* g = head.attribute("gene")) t.gene = g->front();
    t.reverse = head.strand == '-';
    std::size_t cds_lo = SIZE_MAX, cds_hi = 0;
    std::string protein;
    for (std::size_t c : e.children) {
      const auto& child = gff.entities[c];
      for (std::size_t line : child.lines) {
        const auto& f = gff.features[line];
        if (child.type == "exon") t.exons.push_back({f.begin, f.end});
        if (child.type == "CDS") {
          cds_lo = std::min(cds_lo, f.begin);
          cds_hi = std::max(cds_hi, f.end);
          if (const auto* p = f.attribute("protein_id")) protein = p->front();
        }
      }
    }
    if (t.exons.empty()) t.exons.push_back({e.begin, e.end});
    std::sort(t.exons.begin(), t.exons.end(), [](const Exon& a, const Exon& b) { return a.begin < b.begin; });
    if (t.reverse) std::reverse(t.exons.begin(), t.exons.end());
    if (cds_lo != SIZE_MAX) {
      const std::size_t a = detail::along(t.exons, t.reverse, t.reverse ? cds_hi - 1 : cds_lo);
      const std::size_t b = detail::along(t.exons, t.reverse, t.reverse ? cds_lo : cds_hi - 1);
      if (a != SIZE_MAX && b != SIZE_MAX && b >= a) {
        t.coding = true;
        t.cds_begin = a;
        t.cds_end = b + 1;
        t.protein_id = protein;
      }
    }
    out.push_back(std::move(t));
  }
  return out;
}

// ---- numbering ------------------------------------------------------------------------------------------------
struct Mapped {
  enum class Region : std::uint8_t { exon, intron, upstream, downstream };
  Position position;          // as c. or n. writes it
  Region region = Region::exon;
  std::size_t exon = 0;       // the exon it lies in, or the intron after exon `exon`
  bool outside() const { return region == Region::upstream || region == Region::downstream; }
};

namespace detail {

// A transcript coordinate (possibly before 0 or past the end) as c. or n. numbering, with an offset.
inline Position number_of(const Transcript& t, std::int64_t along, char sign, std::int64_t offset) {
  Position p;
  if (t.coding) {
    const auto cb = static_cast<std::int64_t>(t.cds_begin), ce = static_cast<std::int64_t>(t.cds_end);
    if (along < cb) { p.utr = '-'; p.number = cb - along; }
    else if (along >= ce) { p.utr = '*'; p.number = along - ce + 1; }
    else p.number = along - cb + 1;
  } else {
    p.number = along + 1;  // n. upstream of the transcript has no number; the caller marks it outside
    if (along < 0) { p.utr = '-'; p.number = -along; }
  }
  if (sign && offset) {
    p.offset_sign = sign;
    p.offset = offset;
  }
  return p;
}

}  // namespace detail

// The c. or n. numbering of a genomic base (0-based).
inline Mapped map_genomic(const Transcript& t, std::size_t g) {
  Mapped m;
  const auto gi = static_cast<std::int64_t>(g);
  std::int64_t start = 0;
  for (std::size_t k = 0; k < t.exons.size(); k++) {
    const Exon& e = t.exons[k];
    const auto len = static_cast<std::int64_t>(e.end - e.begin);
    if (g >= e.begin && g < e.end) {
      m.region = Mapped::Region::exon;
      m.exon = k;
      m.position = detail::number_of(t, start + (t.reverse ? static_cast<std::int64_t>(e.end) - 1 - gi : gi - static_cast<std::int64_t>(e.begin)), 0, 0);
      return m;
    }
    if (k + 1 < t.exons.size()) {
      const Exon& f = t.exons[k + 1];
      // The intron between e and f, in genomic terms.
      const std::int64_t last = t.reverse ? static_cast<std::int64_t>(e.begin) : static_cast<std::int64_t>(e.end) - 1;     // e's last base
      const std::int64_t next = t.reverse ? static_cast<std::int64_t>(f.end) - 1 : static_cast<std::int64_t>(f.begin);     // f's first base
      const bool inside = t.reverse ? (gi < last && gi > next) : (gi > last && gi < next);
      if (inside) {
        const std::int64_t d5 = t.reverse ? last - gi : gi - last, d3 = t.reverse ? gi - next : next - gi;
        m.region = Mapped::Region::intron;
        m.exon = k;
        if (d5 <= d3) m.position = detail::number_of(t, start + len - 1, '+', d5);
        else m.position = detail::number_of(t, start + len, '-', d3);
        return m;
      }
    }
    start += len;
  }
  // Beyond the ends.
  const Exon& a = t.exons.front();
  const Exon& b = t.exons.back();
  const std::int64_t L = start;
  const bool before = t.reverse ? gi >= static_cast<std::int64_t>(a.end) : gi < static_cast<std::int64_t>(a.begin);
  if (before) {
    m.region = Mapped::Region::upstream;
    m.position = detail::number_of(t, t.reverse ? static_cast<std::int64_t>(a.end) - 1 - gi : gi - static_cast<std::int64_t>(a.begin), 0, 0);
  } else {
    m.region = Mapped::Region::downstream;
    m.exon = t.exons.size() - 1;
    m.position = detail::number_of(t, L - 1 + (t.reverse ? static_cast<std::int64_t>(b.begin) - gi : gi - (static_cast<std::int64_t>(b.end) - 1)), 0, 0);
  }
  return m;
}

// The genomic base (0-based) a c. or n. position names, or nothing with the reason: an offset whose
// anchor is not at an exon's edge facing the intron, an offset beyond its intron, an unknown position.
inline std::optional<std::size_t> map_transcript(const Transcript& t, const Position& p, std::string* why = nullptr) {
  const auto fail = [&](std::string m) -> std::optional<std::size_t> { if (why) *why = std::move(m); return std::nullopt; };
  if (p.kind != Position::Kind::number) return fail("only numbered positions map to the genome");
  if (p.offset_unknown) return fail("an unknown offset maps to no base");
  std::int64_t along;
  if (t.coding) {
    if (p.utr == '-') along = static_cast<std::int64_t>(t.cds_begin) - p.number;
    else if (p.utr == '*') along = static_cast<std::int64_t>(t.cds_end) + p.number - 1;
    else along = static_cast<std::int64_t>(t.cds_begin) + p.number - 1;
  } else {
    if (p.utr) return fail("n. numbering has no '-' or '*'");
    along = p.number - 1;
  }
  const auto L = static_cast<std::int64_t>(t.length());
  const auto genomic_of = [&](std::int64_t a) -> std::int64_t {  // along → genomic, extrapolated beyond the ends
    if (a < 0) {
      const Exon& e = t.exons.front();
      return t.reverse ? static_cast<std::int64_t>(e.end) - 1 - a : static_cast<std::int64_t>(e.begin) + a;
    }
    if (a >= L) {
      const Exon& e = t.exons.back();
      const std::int64_t beyond = a - (L - 1);
      return t.reverse ? static_cast<std::int64_t>(e.begin) - beyond : static_cast<std::int64_t>(e.end) - 1 + beyond;
    }
    std::int64_t start = 0;
    for (const auto& e : t.exons) {
      const auto len = static_cast<std::int64_t>(e.end - e.begin);
      if (a < start + len) return t.reverse ? static_cast<std::int64_t>(e.end) - 1 - (a - start) : static_cast<std::int64_t>(e.begin) + (a - start);
      start += len;
    }
    return -1;
  };
  std::int64_t g = genomic_of(along);
  if (p.offset_sign) {
    // The anchor must be an exon's last base (for '+') or first base (for '-'), facing an intron.
    std::int64_t start = 0;
    bool ok = false;
    std::int64_t intron = 0;
    for (std::size_t k = 0; k < t.exons.size(); k++) {
      const auto len = static_cast<std::int64_t>(t.exons[k].end - t.exons[k].begin);
      if (p.offset_sign == '+' && along == start + len - 1) {
        ok = true;
        if (k + 1 < t.exons.size()) {
          const Exon &e = t.exons[k], &f = t.exons[k + 1];
          intron = t.reverse ? static_cast<std::int64_t>(e.begin) - static_cast<std::int64_t>(f.end) : static_cast<std::int64_t>(f.begin) - static_cast<std::int64_t>(e.end);
        } else {
          intron = INT64_MAX;  // past the last exon: downstream of the transcript
        }
      }
      if (p.offset_sign == '-' && along == start) {
        ok = true;
        if (k > 0) {
          const Exon &e = t.exons[k - 1], &f = t.exons[k];
          intron = t.reverse ? static_cast<std::int64_t>(e.begin) - static_cast<std::int64_t>(f.end) : static_cast<std::int64_t>(f.begin) - static_cast<std::int64_t>(e.end);
        } else {
          intron = INT64_MAX;
        }
      }
      start += len;
    }
    if (!ok) return fail("an offset is counted from an exon's last base (+) or first base (-); " + format_position(p) + " is not at one");
    if (t.genomic && p.offset > intron) return fail("the offset " + format_position(p) + " runs past its intron of " + std::to_string(intron) + " bases");
    if (!t.genomic) return fail("a transcript record holds no introns: " + format_position(p) + " needs a genomic reference");
    const std::int64_t step = (p.offset_sign == '+') != t.reverse ? 1 : -1;
    g += step * p.offset;
  }
  if (g < 0) return fail(format_position(p) + " lies before the record's first base");
  return static_cast<std::size_t>(g);
}

// ---- projection --------------------------------------------------------------------------------------------------
struct ProjectOptions {
  bool repeats = true;
};

struct Projection {
  std::string transcript;
  Coordinate coordinate = Coordinate::c;
  Variant variant;
  bool outside = false;         // some position lies beyond the transcript's ends
  bool intronic = false;        // some position is intronic
  bool exonic = false;          // some base the change names is in an exon
  Change change;                // the genomic change at the placement the transcript's description names (0-based, forward strand)
  std::string text() const { return transcript + ":" + static_cast<char>(coordinate) + "." + format_variant(variant); }
};

namespace detail {

inline Position map_position(const Transcript& t, const Position& window_position, std::size_t w0, std::size_t w1, bool reverse, Projection& pr) {
  if (window_position.kind != Position::Kind::number) return window_position;
  const std::size_t i = static_cast<std::size_t>(window_position.number - 1);
  const std::size_t g = reverse ? w1 - 1 - i : w0 + i;
  const Mapped m = map_genomic(t, g);
  if (m.outside()) pr.outside = true;
  if (m.region == Mapped::Region::intron) pr.intronic = true;
  if (m.region == Mapped::Region::exon) pr.exonic = true;
  return m.position;
}

}  // namespace detail

// The change described on the transcript: placed by the 3′ rule in the transcript's direction, its
// positions in c. or n. numbering.
inline Projection project(const Codes& genome, const Change& change, const Transcript& t, const ProjectOptions& o = {}) {
  Projection pr;
  pr.transcript = t.id;
  pr.coordinate = t.coordinate();
  const Placements pl = placements(genome, change);
  // The window: the class's hull and a margin; on a transcript record, cut at the end of the exon
  // holding the 5′-most placement (numbering.md's exception at exon–exon junctions).
  const std::size_t margin = 64;
  std::size_t w0 = pl.hull_begin > margin ? pl.hull_begin - margin : 0;
  std::size_t w1 = std::min(genome.size(), pl.hull_end + margin);
  if (!t.genomic && t.exons.size() > 1 && !t.reverse) {
    const std::size_t anchor = pl.indel && !pl.deletion ? (pl.first > 0 ? pl.first - 1 : 0) : pl.first;
    for (const auto& e : t.exons)
      if (anchor >= e.begin && anchor < e.end) w1 = std::max(std::min(w1, e.end), std::min(genome.size(), pl.first + pl.minimal.deleted.size()));
  }
  Codes window(genome.begin() + static_cast<std::ptrdiff_t>(w0), genome.begin() + static_cast<std::ptrdiff_t>(w1));
  Change local{pl.minimal.position - w0, pl.minimal.deleted, pl.minimal.inserted};
  if (t.reverse) {
    window = sequence::reverse_complement(window);
    local = Change{window.size() - (local.position + local.deleted.size()), sequence::reverse_complement(local.deleted), sequence::reverse_complement(local.inserted)};
  }
  DescribeOptions d;
  d.repeats = o.repeats;
  d.open_start = t.reverse ? w1 < genome.size() : w0 > 0;
  d.open_end = t.reverse ? w0 > 0 : w1 < genome.size();
  // A repeat inside the CDS is written only for units of a length divisible by three.
  if (t.coding) {
    const Mapped a = map_genomic(t, pl.hull_begin), b = map_genomic(t, pl.hull_end > 0 ? pl.hull_end - 1 : 0);
    const bool coding_a = a.region == Mapped::Region::exon && !a.position.utr, coding_b = b.region == Mapped::Region::exon && !b.position.utr;
    d.unit_multiple_of_three = coding_a || coding_b;
  }
  const Described described = describe(window, local, d);
  pr.variant = described.variant;
  const auto map = [&](Position& p) { p = detail::map_position(t, p, w0, w1, t.reverse, pr); };
  for (Position* p : {&pr.variant.location.start.first, &pr.variant.location.start.last, &pr.variant.location.end.first, &pr.variant.location.end.last}) map(*p);
  // The genomic change at the described placement.
  const Change& c = described.change;
  if (!t.reverse) pr.change = Change{w0 + c.position, c.deleted, c.inserted};
  else pr.change = Change{w1 - (c.position + c.deleted.size()), sequence::reverse_complement(c.deleted), sequence::reverse_complement(c.inserted)};
  return pr;
}

}  // namespace genomics::nomenclature
