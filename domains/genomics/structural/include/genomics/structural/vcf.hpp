// Structural variant calls as VCF 4.2 records (The Variant Call Format Specification, VCFv4.2,
// samtools/hts-specs, §1.4.5 "Alternative allele field format", §3 "INFO keys used for structural
// variants", §5 "Representing variation in VCF records", §5.3 symbolic alleles and §5.4 breakends),
// on variant/vcf.hpp's records and writer.
//
// - Symbolic alleles. A deletion, tandem duplication or inversion is written at the base before
//   it (POS, 1-based, with that base as REF) with ALT <DEL>, <DUP:TANDEM> or <INV>, INFO
//   SVTYPE, END (the last base changed, 1-based), SVLEN (negative for a deletion), and, when
//   imprecise, IMPRECISE with CIPOS and CIEND (the confidence intervals relative to POS and END);
//   a precise one carries HOMLEN, the junction's microhomology. A copy-number call is <CNV> with
//   its copy number in FORMAT CN. A resolved insertion is written with its bases (REF the base
//   before, ALT that base and the inserted ones); an unresolved one as <INS> with SVLEN a lower
//   bound and IMPRECISE.
// - Breakends (§5.4). A translocation, a dispersed duplication and an unpaired junction between
//   contigs are written as breakend records, two for each junction (one at each of its breakends,
//   mates of each other by MATEID, the records of one event sharing EVENT). For a breakend
//   (c, p, +) with base t = ref[p], joined to (c′, q, −) the ALT is t[c′:q[ ("the piece extending
//   to the right of q is joined after t"), to (c′, q, +) it is t]c′:q] (the reverse complement of
//   the piece extending to the left of q); for (c, p, −), ]c′:q]t and [c′:q[t mirror them. Novel
//   bases at the junction stand between t and the bracket, read in the record's own direction.
//
// Coordinates are 0-based inside and 1-based in the records, as everywhere in the domain.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/structural/breakend.hpp>
#include <genomics/structural/call.hpp>
#include <genomics/variant/vcf.hpp>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

namespace genomics::structural {

// The header of the structural calls' VCF: contigs, ALT, INFO and FORMAT definitions, one sample.
inline variant::VcfHeader sv_header(const std::vector<std::string>& names, const std::vector<std::size_t>& lengths, std::string sample, std::string source = "genomics::structural") {
  variant::VcfHeader h;
  h.meta.push_back("source=" + source);
  for (std::size_t i = 0; i < names.size(); i++) {
    h.meta.push_back("contig=<ID=" + names[i] + ",length=" + std::to_string(lengths[i]) + ">");
    h.contigs.push_back({names[i], lengths[i]});
  }
  for (const char* line : {
           "ALT=<ID=DEL,Description=\"Deletion relative to the reference\">",
           "ALT=<ID=DUP:TANDEM,Description=\"Tandem duplication\">",
           "ALT=<ID=INV,Description=\"Inversion of reference sequence\">",
           "ALT=<ID=INS,Description=\"Insertion of novel sequence\">",
           "ALT=<ID=CNV,Description=\"Copy number variable region\">",
           "INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
           "INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
           "INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
           "INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">",
           "INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">",
           "INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">",
           "INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">",
           "INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">",
           "INFO=<ID=EVENT,Number=1,Type=String,Description=\"ID of event associated to breakend\">",
           "INFO=<ID=PE,Number=1,Type=Integer,Description=\"Discordant pairs supporting the variant\">",
           "INFO=<ID=SR,Number=1,Type=Integer,Description=\"Split reads and gapped alignments supporting the variant\">",
           "INFO=<ID=SC,Number=1,Type=Integer,Description=\"Soft-clipped reads supporting the variant\">",
           "FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
           "FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">",
           "FORMAT=<ID=CN,Number=1,Type=Float,Description=\"Copy number measured by read depth\">",
           "FORMAT=<ID=RO,Number=1,Type=Float,Description=\"Reference support, weighted by placement\">",
           "FORMAT=<ID=AO,Number=1,Type=Float,Description=\"Alternate support, weighted by placement\">",
       })
    h.meta.push_back(line);
  h.samples.push_back(std::move(sample));
  return h;
}

namespace detail {

inline std::string fixed(double x, int digits) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.*f", digits, x);
  return buf;
}

inline std::string genotype_text(int g) { return g == 2 ? "1/1" : g == 1 ? "0/1" : "0/0"; }

inline void sample_columns(variant::VcfRecord& r, const SvCall& c) {
  r.format = {"GT", "GQ", "CN", "RO", "AO"};
  r.samples = {{genotype_text(c.genotype.genotype), std::to_string(c.genotype.gq), c.copies >= 0 ? fixed(c.copies, 2) : ".", fixed(c.genotype.reference, 1), fixed(c.genotype.alternate, 1)}};
}

// The ALT of the breakend record at `a` of the junction leaving at a and entering at b.
inline std::string breakend_alt(const std::vector<Codes>& contigs, const std::vector<std::string>& names, const Adjacency& j) {
  const char t = sequence::dna_letter(contigs[j.a.contig][j.a.position]);
  const std::string mate = names[j.b.contig] + ":" + std::to_string(j.b.position + 1);
  if (j.a.plus) {
    const std::string tail = j.b.plus ? "]" + mate + "]" : "[" + mate + "[";
    return std::string(1, t) + j.inserted + tail;
  }
  const std::string head = j.b.plus ? "]" + mate + "]" : "[" + mate + "[";
  return head + sequence::reverse_complement(j.inserted) + std::string(1, t);
}

}  // namespace detail

// The records of the calls, in order; `names` names the contigs.
inline std::vector<variant::VcfRecord> sv_records(const std::vector<Codes>& contigs, const std::vector<std::string>& names, const std::vector<SvCall>& calls) {
  std::vector<variant::VcfRecord> out;
  for (std::size_t k = 0; k < calls.size(); k++) {
    const SvCall& c = calls[k];
    const std::string id = "sv" + std::to_string(k + 1);
    const auto info_common = [&](variant::VcfRecord& r) {
      if (c.pairs) r.info.push_back({"PE", std::to_string(c.pairs)});
      if (c.splits) r.info.push_back({"SR", std::to_string(c.splits)});
      if (c.clips) r.info.push_back({"SC", std::to_string(c.clips)});
    };
    const bool breakends = c.kind == SvKind::translocation || c.kind == SvKind::dispersed_duplication || c.kind == SvKind::breakend;
    if (breakends) {
      for (std::size_t j = 0; j < c.junctions.size(); j++) {
        const Adjacency& x = c.junctions[j];
        for (int side = 0; side < 2; side++) {
          const Adjacency y = side == 0 ? x : reversed(x);
          variant::VcfRecord r;
          r.chrom = names[y.a.contig];
          r.position = y.a.position;
          r.id = id + "_" + std::to_string(j + 1) + (side ? "b" : "a");
          r.ref = std::string(1, sequence::dna_letter(contigs[y.a.contig][y.a.position]));
          r.alts = {detail::breakend_alt(contigs, names, y)};
          r.qual = std::round(c.qual * 100) / 100;
          r.filters = c.passed() ? std::vector<std::string>{"PASS"} : c.filters;
          r.info = {{"SVTYPE", "BND"}, {"MATEID", id + "_" + std::to_string(j + 1) + (side ? "a" : "b")}, {"EVENT", id}};
          if (!c.precise) r.info.push_back({"IMPRECISE", ""});
          info_common(r);
          detail::sample_columns(r, c);
          out.push_back(std::move(r));
        }
      }
      continue;
    }
    variant::VcfRecord r;
    r.chrom = names[c.contig];
    r.id = id;
    r.qual = std::round(c.qual * 100) / 100;
    r.filters = c.passed() ? std::vector<std::string>{"PASS"} : c.filters;
    const std::size_t before = c.begin > 0 ? c.begin - 1 : 0;
    r.position = before;
    r.ref = std::string(1, sequence::dna_letter(contigs[c.contig][before]));
    const auto symbolic = [&](const char* alt, const char* type, long long svlen) {
      r.alts = {alt};
      r.info = {{"SVTYPE", type}, {"END", std::to_string(c.end)}, {"SVLEN", std::to_string(svlen)}};
      if (!c.precise) {
        r.info.push_back({"IMPRECISE", ""});
        const auto rel = [](std::size_t lo, std::size_t hi, std::size_t at) {
          return std::to_string(static_cast<long long>(lo) - static_cast<long long>(at)) + "," + std::to_string(static_cast<long long>(hi) - static_cast<long long>(at));
        };
        r.info.push_back({"CIPOS", rel(c.ci_begin[0], c.ci_begin[1], c.ci_begin[0] <= before && before <= c.ci_begin[1] ? before : c.ci_begin[0])});
        r.info.push_back({"CIEND", rel(c.ci_end[0], c.ci_end[1], c.ci_end[0] <= c.end && c.end <= c.ci_end[1] ? c.end : c.ci_end[0])});
      } else if (!c.junctions.empty()) {
        r.info.push_back({"HOMLEN", std::to_string(homology(contigs, c.junctions.front()).length())});
      }
    };
    switch (c.kind) {
      case SvKind::deletion: symbolic("<DEL>", "DEL", -static_cast<long long>(c.length)); break;
      case SvKind::duplication: symbolic("<DUP:TANDEM>", "DUP", static_cast<long long>(c.length)); break;
      case SvKind::inversion: symbolic("<INV>", "INV", static_cast<long long>(c.length)); break;
      case SvKind::copy_loss:
      case SvKind::copy_gain: symbolic("<CNV>", "CNV", static_cast<long long>(c.length)); break;
      case SvKind::insertion:
        if (c.resolved) {
          r.alts = {r.ref + c.sequence};
          r.info = {{"SVTYPE", "INS"}, {"END", std::to_string(c.begin)}, {"SVLEN", std::to_string(c.sequence.size())}};
        } else {
          r.alts = {"<INS>"};
          r.info = {{"SVTYPE", "INS"}, {"END", std::to_string(c.begin)}, {"SVLEN", std::to_string(c.length)}, {"IMPRECISE", ""}};
        }
        break;
      default: break;
    }
    info_common(r);
    detail::sample_columns(r, c);
    out.push_back(std::move(r));
  }
  return out;
}

}  // namespace genomics::structural
