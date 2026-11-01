// The calls said in the international notations: every detected abnormality as a VCF record (its
// symbolic allele or its breakends), as an HGVS description of the change on the chromosome's
// accession (nomenclature/describe.hpp: the 3′ rule, del, dup, inv, ins, the copy of a range), and
// in ISCN with the ideogram's bands (karyotype/describe.hpp): the copy numbers and runs of
// homozygosity as a microarray description, the junctions as a sequencing description with its
// derivative chromosomes ("from where to where is lacking, which part swapped to which part").
//
// - Where a contig is. The analysis may run on regions of chromosomes; each contig is placed on
//   its chromosome by its accession (NC_000022.11), its karyotype index (21 for chromosome 22) and
//   the offset of its first base, and every position is written in the chromosome's numbering.
// - HGVS (HGVS Nomenclature, DNA/deletion, duplication, inversion, insertion). A precise deletion,
//   tandem duplication or resolved insertion of up to `exact_limit` bases is described by the
//   nomenclature module from the change itself, so that the 3′ rule and the choice of dup or ins
//   are its; longer ones are written g.a_bdel, g.a_bdup; an inversion g.a_binv; a dispersed copy
//   as an insertion of the copied range (ins[NC_…:g.a_b], inverted when it is); an imprecise
//   event with its confidence intervals as uncertain positions, g.(a₁_a₂)_(b₁_b₂)del; a
//   copy-number call from depth alone as g.(?_a)_(b_?)del or dup; an unresolved insertion as
//   ins(?). Translocations and other breakends are not given an HGVS form here (the HGVS
//   recommendations write them with "::" between references, which the nomenclature module does not
//   yet produce); their VCF breakends and ISCN sequencing description say them.
// - ISCN. The copy-number segments other than two copies, and the copy-neutral runs of
//   homozygosity (hmz), are an arr description; the junctions of the passing calls are a seq
//   description, the karyotype module reading the derivative chromosomes from them (an inversion's
//   or a translocation's two junctions put on one cut per chromosome). Each call is also described
//   alone (t(21;22)(q22.11;q12.3), del(21)(q22.11q22.11)); the whole sample's junctions together
//   are one description, which the karyotype module names complex, (21,22)cx, when several events
//   share chromosomes.
//
// Complexity: linear in the calls and segments, and the nomenclature module's for each change.

#pragma once

#include <genomics/karyotype/bands.hpp>
#include <genomics/karyotype/describe.hpp>
#include <genomics/nomenclature/describe.hpp>
#include <genomics/structural/allelic.hpp>
#include <genomics/structural/call.hpp>
#include <genomics/structural/copy_number.hpp>
#include <genomics/structural/vcf.hpp>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace genomics::structural {

// Where a contig of the analysis lies: its chromosome's accession, its karyotype index (−1 when it
// is not a chromosome of the ideogram) and the chromosome position of its first base (0-based).
struct Placement {
  std::string accession;
  int chromosome = -1;
  std::size_t offset = 0;
};

struct CallDescription {
  std::string hgvs;              // empty when HGVS is not written for its kind
  std::string iscn;              // arr for a copy-number change, seq for a junction event
  std::vector<std::string> vcf;  // its records' lines
};

struct SampleDescription {
  std::vector<CallDescription> calls;
  std::string arr;                 // the copy numbers and runs of homozygosity
  karyotype::SeqDescription seq;   // the junctions and the derivative chromosomes
};

struct DescribeOptions {
  std::size_t exact_limit = 100000;  // precise changes up to this length go to the nomenclature module
  std::string build = "GRCh38";
};

namespace detail {

inline std::string number(std::size_t x) { return std::to_string(x); }

// The karyotype adjacency of a junction, in chromosome coordinates (boundaries between bases).
inline karyotype::Adjacency karyotype_adjacency(const Adjacency& junction, const std::vector<Placement>& at) {
  // Read from a piece's right end when the junction allows it (a tandem duplication as [0, b)::[a, L)).
  const Adjacency j = !junction.a.plus && junction.b.plus ? reversed(junction) : junction;
  karyotype::Adjacency k;
  const Placement &pa = at[j.a.contig], &pb = at[j.b.contig];
  k.chromosome_a = pa.chromosome;
  k.a_left = j.a.plus;
  k.break_a = pa.offset + (j.a.plus ? j.a.position + 1 : j.a.position);
  k.chromosome_b = pb.chromosome;
  k.b_right = !j.b.plus;
  k.break_b = pb.offset + (j.b.plus ? j.b.position + 1 : j.b.position);
  k.bases = j.inserted;
  return k;
}

// The junctions of a call, an inversion's or a translocation's put on one cut per side: their
// breakends may lie a few bases apart (microhomology, small losses at the breaks), and ISCN describes
// one break of each chromosome.
inline std::vector<Adjacency> cut_junctions(const SvCall& c) {
  std::vector<Adjacency> out = c.junctions;
  if (c.kind != SvKind::translocation && c.kind != SvKind::inversion) return out;
  const auto d = [](std::size_t x, std::size_t y) { return x > y ? x - y : y - x; };
  const auto snap = [&](Breakend& e) {
    // The cut a breakend belongs to: its contig's, or for an inversion the nearer of its two.
    const std::size_t cut = c.kind == SvKind::translocation ? (e.contig == c.contig ? c.begin : c.partner_position)
                                                             : (d(e.position, c.begin) <= d(e.position, c.end) ? c.begin : c.end);
    e.position = e.plus ? cut - 1 : cut;
  };
  for (auto& j : out) {
    snap(j.a);
    snap(j.b);
  }
  return out;
}

}  // namespace detail

// The HGVS description of one call (see the header), or "" for kinds it does not write.
inline std::string hgvs_of(const std::vector<Codes>& contigs, const std::vector<Placement>& at, const SvCall& c, const DescribeOptions& o = {}) {
  const Placement& p = at[c.contig];
  const std::string acc = p.accession + ":g.";
  const Codes& ref = contigs[c.contig];
  const std::size_t b1 = p.offset + c.begin + 1, e1 = p.offset + c.end;  // 1-based first and last base
  nomenclature::DescribeOptions no;
  no.first = static_cast<std::int64_t>(p.offset) + 1;
  no.open_start = p.offset > 0;
  no.open_end = true;
  const auto letters = [&](std::size_t a, std::size_t b) {
    std::string s;
    for (std::size_t i = a; i < b; i++) s += sequence::dna_letter(ref[i]);
    return s;
  };
  const auto uncertain = [&](const char* edit) {
    // CIPOS on the first changed base and CIEND on the last, 1-based.
    return acc + "(" + detail::number(p.offset + c.ci_begin[0] + 2) + "_" + detail::number(p.offset + c.ci_begin[1] + 2) + ")_(" + detail::number(p.offset + c.ci_end[0]) + "_" +
           detail::number(p.offset + c.ci_end[1]) + ")" + edit;
  };
  switch (c.kind) {
    case SvKind::deletion:
      if (!c.precise) return uncertain("del");
      if (c.length <= o.exact_limit) return nomenclature::describe_text(ref, {c.begin, letters(c.begin, c.end), c.sequence}, p.accession, nomenclature::Coordinate::g, no);
      return acc + detail::number(b1) + "_" + detail::number(e1) + "del";
    case SvKind::duplication:
      if (!c.precise) return uncertain("dup");
      if (c.length <= o.exact_limit) return nomenclature::describe_text(ref, {c.end, "", letters(c.begin, c.end)}, p.accession, nomenclature::Coordinate::g, no);
      return acc + detail::number(b1) + "_" + detail::number(e1) + "dup";
    case SvKind::inversion:
      if (!c.precise) return uncertain("inv");
      return acc + detail::number(b1) + "_" + detail::number(e1) + "inv";
    case SvKind::insertion:
      if (c.resolved && !c.sequence.empty()) return nomenclature::describe_text(ref, {c.begin, "", c.sequence}, p.accession, nomenclature::Coordinate::g, no);
      return acc + detail::number(p.offset + c.begin) + "_" + detail::number(p.offset + c.begin + 1) + "ins(?)";
    case SvKind::dispersed_duplication: {
      const Placement& q = at[c.partner];
      const std::size_t site = q.offset + c.partner_position;  // 1-based: the bases site and site + 1 flank it
      return q.accession + ":g." + detail::number(site) + "_" + detail::number(site + 1) + "ins[" + p.accession + ":g." + detail::number(b1) + "_" + detail::number(e1) + (c.inverted ? "inv" : "") + "]";
    }
    case SvKind::copy_loss:
    case SvKind::copy_gain:
      return acc + "(?_" + detail::number(b1) + ")_(" + detail::number(e1) + "_?)" + (c.kind == SvKind::copy_loss ? "del" : "dup");
    default: return "";
  }
}

// Every call described, and the sample's ISCN descriptions; `ideogram` may be null (no ISCN then).
inline SampleDescription describe_sample(const std::vector<Codes>& contigs, const std::vector<Placement>& at, const std::vector<SvCall>& calls,
                                         const std::vector<CopySegmentCall>& segments, const std::vector<AllelicSegment>& allelic, const karyotype::Ideogram* ideogram,
                                         const DescribeOptions& o = {}) {
  SampleDescription out;
  // VCF, in the chromosomes' numbering.
  std::vector<std::string> names;
  for (const auto& p : at) names.push_back(p.accession);
  auto records = sv_records(contigs, names, calls);
  for (auto& r : records)
    for (std::size_t k = 0; k < at.size(); k++)
      if (r.chrom == names[k]) r.position += at[k].offset;
  karyotype::SeqOptions so;
  so.build = o.build;
  karyotype::ArrOptions ao;
  ao.build = o.build;
  const bool iscn = ideogram != nullptr;
  for (std::size_t i = 0; i < calls.size(); i++) {
    const SvCall& c = calls[i];
    CallDescription d;
    d.hgvs = hgvs_of(contigs, at, c, o);
    const std::string id = "sv" + std::to_string(i + 1);
    for (const auto& r : records)
      if (r.id == id || r.id.rfind(id + "_", 0) == 0) d.vcf.push_back(variant::vcf_line(r));
    if (iscn && at[c.contig].chromosome >= 0) {
      if (c.kind == SvKind::copy_loss || c.kind == SvKind::copy_gain) {
        karyotype::CopySegment s;
        s.chromosome = at[c.contig].chromosome;
        s.start = at[c.contig].offset + c.begin;
        s.end = at[c.contig].offset + c.end;
        s.copies = static_cast<int>(std::lround(c.copies));
        d.iscn = karyotype::describe_arr({s}, *ideogram, ao);
      } else if (!c.junctions.empty()) {
        std::vector<karyotype::Adjacency> adj;
        bool placed = true;
        for (const auto& j : detail::cut_junctions(c)) {
          placed = placed && at[j.a.contig].chromosome >= 0 && at[j.b.contig].chromosome >= 0;
          adj.push_back(detail::karyotype_adjacency(j, at));
        }
        if (placed) d.iscn = karyotype::describe_seq(adj, *ideogram, so).short_form;
      }
    }
    out.calls.push_back(std::move(d));
  }
  if (iscn) {
    std::vector<karyotype::CopySegment> segs;
    for (const auto& s : segments) {
      if (s.copies == 2 || at[s.contig].chromosome < 0) continue;
      karyotype::CopySegment k;
      k.chromosome = at[s.contig].chromosome;
      k.start = at[s.contig].offset + s.begin;
      k.end = at[s.contig].offset + s.end;
      k.copies = s.copies;
      segs.push_back(k);
    }
    for (const auto& a : allelic) {
      if (a.state != AllelicState::copy_neutral_loh || at[a.contig].chromosome < 0) continue;
      karyotype::CopySegment k;
      k.chromosome = at[a.contig].chromosome;
      k.start = at[a.contig].offset + a.begin;
      k.end = at[a.contig].offset + a.end;
      k.copies = 2;
      k.homozygous = true;
      segs.push_back(k);
    }
    out.arr = karyotype::describe_arr(segs, *ideogram, ao);
    std::vector<karyotype::Adjacency> adj;
    for (const auto& c : calls) {
      if (!c.passed()) continue;
      for (const auto& j : detail::cut_junctions(c))
        if (at[j.a.contig].chromosome >= 0 && at[j.b.contig].chromosome >= 0) adj.push_back(detail::karyotype_adjacency(j, at));
    }
    out.seq = karyotype::describe_seq(adj, *ideogram, so);
  }
  return out;
}

}  // namespace genomics::structural
