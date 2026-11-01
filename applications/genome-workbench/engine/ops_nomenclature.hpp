// The operations of the nomenclature: a description read by the HGVS grammar (its value, its named
// tree for drawing the derivation, its diagnostics); a change on a real reference described in every
// form (HGVS by the 3′ rule, VCF's left form, SPDI's hull, the interval of placements) and projected
// through every transcript to its protein and consequences; VCF's structural-variant and breakend
// records as adjacencies, segments and HGVS descriptions.
//
// Every computation is the domain's (genomics::nomenclature); the engine reads the inputs and writes
// the answers.

#pragma once

#include "genomes.hpp"
#include "inputs.hpp"

#include <genomics/nomenclature/hgvs.hpp>
#include <genomics/nomenclature/protein.hpp>
#include <genomics/nomenclature/spdi.hpp>
#include <genomics/nomenclature/structural.hpp>
#include <genomics/nomenclature/transcript.hpp>

#include <algorithm>
#include <string>
#include <vector>

namespace workbench {

namespace gnm = genomics::nomenclature;

namespace detail {

inline Json tree_json(const morphe::formal::ParseNode& n) {
  Json c = Json::array();
  for (const auto& x : n.children) c.push(tree_json(x));
  return Json::object({{"name", n.name}, {"begin", n.begin}, {"end", n.end}, {"text", n.text}, {"children", c}});
}

inline Json position_json(const gnm::Position& p) {
  static const char* kinds[] = {"number", "unknown", "pter", "qter", "cen"};
  Json j = Json::object({{"text", gnm::format_position(p)}, {"kind", kinds[static_cast<int>(p.kind)]}, {"number", p.number}});
  if (p.utr) j.set("utr", std::string(1, p.utr));
  if (p.offset_sign) j.set("offset", p.offset_unknown ? Json("?") : Json(p.signed_offset()));
  return j;
}

template <class P, class F>
Json point_json(const gnm::Point<P>& p, F&& f) {
  Json j = Json::object({{"first", f(p.first)}, {"uncertain", p.uncertain}});
  if (p.span) j.set("last", f(p.last));
  return j;
}

inline Json variant_json(const gnm::Variant& v) {
  Json j = Json::object({{"kind", gnm::kind_name(v.kind)}, {"protein", v.protein}, {"text", gnm::format_variant(v)}});
  if (v.protein) {
    const auto f = [](const gnm::Residue& r) { return Json::object({{"residue", std::string(gnm::three_letter(r.residue))}, {"number", r.number}}); };
    j.set("start", point_json(v.residues.start, f));
    if (v.residues.range) j.set("end", point_json(v.residues.end, f));
  } else {
    j.set("start", point_json(v.location.start, position_json));
    if (v.location.range) j.set("end", point_json(v.location.end, position_json));
  }
  if (v.mixture) j.set("mixture", v.mixture == 1 ? "mosaic" : "chimeric");
  if (!v.was.empty()) j.set("was", v.was);
  if (!v.now.empty()) j.set("now", v.now);
  if (!v.stated.empty()) j.set("stated", v.stated);
  if (!v.inserted.empty()) {
    Json ins = Json::array();
    for (const auto& p : v.inserted) {
      static const char* kinds[] = {"sequence", "location", "remote"};
      Json x = Json::object({{"kind", kinds[static_cast<int>(p.kind)]}, {"inverted", p.inverted}});
      if (p.kind == gnm::Piece::Kind::sequence) x.set("sequence", p.sequence);
      else x.set("location", gnm::format_location(p.location, gnm::format_position));
      if (p.kind == gnm::Piece::Kind::remote) x.set("reference", p.reference.text());
      if (p.count) x.set("count", gnm::format_copies(*p.count));
      ins.push(x);
    }
    j.set("inserted", ins);
  }
  if (!v.residues_inserted.empty()) j.set("residues_inserted", v.residues_inserted);
  if (!v.units.empty()) {
    Json u = Json::array();
    for (const auto& x : v.units) u.push(Json::object({{"unit", x.sequence}, {"copies", gnm::format_copies(x.copies)}}));
    j.set("units", u);
  }
  if (v.copies) j.set("copies", gnm::format_copies(*v.copies));
  if (!v.fs_new.empty()) j.set("new_residue", v.fs_new);
  if (v.has_stop) j.set("new_stop", v.stop_unknown ? Json("?") : Json(v.stop_at));
  return j;
}

inline Json chain_json(const gnm::Chain& c) {
  Json items = Json::array();
  for (const auto& it : c.items) {
    Json ch = Json::array();
    for (const auto& v : it.choices) ch.push(variant_json(v));
    items.push(Json::object({{"predicted", it.predicted}, {"mixed", it.mixed ? Json(it.mixed == 1 ? "mosaic" : "chimeric") : Json()}, {"choices", ch}}));
  }
  return items;
}

inline Json description_json(const gnm::Description& d) {
  Json j = Json::object({{"coordinate", std::string(1, static_cast<char>(d.variants.coordinate))}, {"text", gnm::format_description(d)}});
  if (d.has_reference)
    j.set("reference", Json::object({{"accession", d.reference.accession}, {"version", d.reference.version}, {"selector", d.reference.has_selector ? Json(d.reference.selector()) : Json()}}));
  if (d.adjoined) {
    j.set("adjoined", Json::object({{"five", d.five.reference.text()}, {"three", d.three.reference.text()}, {"linker", d.linker}}));
    return j;
  }
  if (!d.variants.whole.empty()) j.set("whole", d.variants.whole);
  Json alts = Json::array();
  for (const auto& a : d.variants.alternatives) {
    Json groups = Json::array();
    for (const auto& g : a.groups) {
      Json gj = Json::object({{"bracketed", g.bracketed}});
      if (g.bracketed) {
        Json alleles = Json::array();
        for (const auto& al : g.alleles) {
          Json chains = Json::array();
          for (const auto& ch : al.chains) chains.push(chain_json(ch));
          std::string seps;
          for (char c : al.separators) seps += c == '?' ? "(;)" : std::string(1, c);
          alleles.push(Json::object({{"predicted", al.predicted}, {"whole", al.whole}, {"separators", seps}, {"supernumerary", al.supernumerary}, {"chains", chains}}));
        }
        gj.set("alleles", alleles);
      } else {
        gj.set("chain", chain_json(g.chain));
      }
      groups.push(gj);
    }
    alts.push(Json::object({{"groups", groups}}));
  }
  j.set("alternatives", alts);
  return j;
}

inline const seq::Codes& contig_codes(const Genome& g, const Json& request, std::string& accession) {
  const Json& c = request["contig"];
  if (c.is_string()) {
    for (const auto& x : g.contigs)
      if (x.accession == c.text() || x.name == c.text()) { accession = x.accession; return x.codes; }
    throw Fault("no contig " + c.text());
  }
  const Contig& x = g.contig(integer(request, "contig", 0));
  accession = x.accession;
  return x.codes;
}

inline seq::GenBank record_of(const Genome& g, const Json& request) {
  const Json& c = request["contig"];
  const Contig* x = nullptr;
  if (c.is_string()) {
    for (const auto& y : g.contigs)
      if (y.accession == c.text() || y.name == c.text()) x = &y;
  } else {
    x = &g.contig(integer(request, "contig", 0));
  }
  if (!x) throw Fault("no contig");
  seq::GenBank r;
  r.locus = x->name;
  r.accession = r.version = x->accession;
  r.definition = x->definition;
  r.features = x->features;
  return r;
}

}  // namespace detail

// {"text": "NM_004006.2:c.4375C>T"} → {"accepted", "tree" (morphe's named tree), "furthest", "expected",
// "error", "canonical", "value" (descriptions), "reading", "diagnostics": [{severity, message, where}]}.
BENCH_OP("nomenclature.parse") {
  const std::string text = str(request, "text");
  if (text.size() > 4000) throw Fault("a description of at most 4000 characters");
  const gnm::HgvsParse p = gnm::parse_hgvs(text);
  Json expected = Json::array();
  for (const auto& e : p.expected) expected.push(e);
  Json out = Json::object({{"accepted", p.accepted}, {"furthest", p.furthest}, {"expected", expected}, {"error", p.error}, {"length", text.size()}});
  if (!p.accepted) return out;
  out.set("tree", detail::tree_json(p.tree));
  out.set("canonical", p.canonical);
  out.set("reading", gnm::reading(p.value));
  Json ds = Json::array();
  for (const auto& d : p.value.alternatives) ds.push(detail::description_json(d));
  for (const auto& a : p.value.alleles)
    for (const auto& d : a) ds.push(detail::description_json(d));
  out.set("value", Json::object({{"genome_alleles", p.value.genome_alleles}, {"descriptions", ds}}));
  Json diag = Json::array();
  static const char* sev[] = {"error", "warning", "note"};
  for (const auto& d : p.diagnostics) diag.push(Json::object({{"severity", sev[static_cast<int>(d.severity)]}, {"message", d.message}, {"where", d.where}}));
  out.set("diagnostics", diag);
  out.set("ok", p.ok());
  return out;
}

// {"genome": {"file": …GenBank of the data library…}, "contig": 0, "first": number of the record's first base
// (a region's offset + 1), and either "description": "NC_000011.10:g.5227002T>A" or "position" (numbered
// from first), "deleted", "inserted"; "transcripts": [ids] (default: those within 5 kb)} →
// {"accession", "coordinate", "hgvs", "change": {position, deleted, inserted}, "placements": {first, last,
// hull: [begin, end]} (numbered from first; insertions between bases), "vcf": {pos, ref, alt}, "spdi",
// "transcripts": [{id, gene, protein, strand, c, p, r, consequences, impact, note, intronic, outside}]}.
BENCH_OP("nomenclature.describe") {
  const auto genome = worker.genome(request);
  std::string accession;
  const seq::Codes& codes = detail::contig_codes(*genome, request, accession);
  const std::int64_t first = integer(request, "first", 1);
  gnm::Change change;
  if (request.has("description")) {
    const auto p = gnm::parse_hgvs(str(request, "description"));
    if (!p.accepted) throw Fault(p.error);
    const auto& alts = p.value.alternatives;
    if (alts.empty() || alts[0].variants.alternatives.empty() || alts[0].variants.alternatives[0].groups.empty() ||
        alts[0].variants.alternatives[0].groups[0].bracketed || alts[0].variants.alternatives[0].groups[0].chain.items.size() != 1)
      throw Fault("one exact g. or m. variant");
    if (!gnm::genomic(alts[0].variants.coordinate)) throw Fault("a g., m. or o. description (c. needs its transcript's mapping first)");
    const auto c = gnm::change_of(codes, alts[0].variants.alternatives[0].groups[0].chain.items[0].choices[0], first);
    if (!c) throw Fault("the description is not exact, or its reference bases are not this record's");
    change = *c;
  } else {
    const std::int64_t pos = integer(request, "position", 0) - first;
    if (pos < 0 || pos > static_cast<std::int64_t>(codes.size())) throw Fault("the position lies outside the record");
    change = gnm::Change{static_cast<std::size_t>(pos), str(request, "deleted"), str(request, "inserted")};
    if (!gnm::deleted_agrees(codes, change)) throw Fault("the deleted bases are not the reference's at that position");
  }
  // m. for a mitochondrial genome, o. for another circular one, g. otherwise (numbering.md).
  const Contig& head = genome->contigs.front();
  const char type = head.definition.find("mitochondri") != std::string::npos ? 'm' : head.circular ? 'o' : 'g';
  gnm::DescribeOptions o;
  o.first = first;
  const gnm::Described d = gnm::describe(codes, change, o);
  const gnm::Placements& pl = d.placements;
  const auto vcf = gnm::vcf_form(codes, change);
  gnm::Spdi sp = gnm::spdi_form(codes, change, accession);
  sp.position += static_cast<std::size_t>(first - 1);
  const auto num = [&](std::size_t i) { return static_cast<long long>(i) + first; };
  Json out = Json::object({{"accession", accession}, {"coordinate", std::string(1, type)},
                           {"hgvs", accession + ":" + type + "." + d.text()},
                           {"kind", gnm::kind_name(d.variant.kind)},
                           {"change", Json::object({{"position", num(d.change.position)}, {"deleted", d.change.deleted}, {"inserted", d.change.inserted}})},
                           {"placements", Json::object({{"indel", pl.indel}, {"first", num(pl.first)}, {"last", num(pl.last)},
                                                        {"hull", Json::array({num(pl.hull_begin), num(pl.hull_end) - 1})}, {"length", pl.length}})},
                           {"vcf", Json::object({{"pos", num(vcf.position)}, {"ref", vcf.ref}, {"alt", vcf.alts.empty() ? std::string() : vcf.alts[0]}})},
                           {"spdi", sp.text()}});
  // The letters about the class, for drawing its placements: about 20 bases on each side of the hull.
  {
    const std::size_t b = pl.hull_begin > 20 ? pl.hull_begin - 20 : 0, e = std::min(codes.size(), pl.hull_end + 20);
    std::string letters;
    for (std::size_t i = b; i < e; i++) letters += seq::dna_letter(codes[i]);
    out.set("context", Json::object({{"begin", num(b)}, {"letters", letters}}));
  }
  // Transcripts.
  const auto record = detail::record_of(*genome, request);
  const auto models = gnm::transcripts(record);
  std::vector<std::string> wanted;
  if (request["transcripts"].is_array())
    for (const auto& x : request["transcripts"].as_array()) wanted.push_back(x.text());
  Json ts = Json::array();
  for (const auto& t : models) {
    if (!wanted.empty() && std::find(wanted.begin(), wanted.end(), t.id) == wanted.end()) continue;
    if (wanted.empty() && (pl.hull_end + 5000 < t.lowest() || pl.hull_begin > t.highest() + 5000)) continue;
    const gnm::ProteinPrediction pp = gnm::predict(codes, change, t);
    Json terms = Json::array();
    for (const auto& x : pp.consequences) terms.push(x);
    Json exons = Json::array();  // numbered from first, [begin, end] inclusive, in the transcript's order
    for (const auto& e : t.exons) exons.push(Json::array({num(e.begin), num(e.end) - 1}));
    Json cds;
    if (t.coding) {
      // The CDS's genomic ends: its first and last bases along the transcript.
      std::size_t along = 0, lo = SIZE_MAX, hi = 0;
      for (const auto& e : t.exons) {
        for (std::size_t k = 0; k < e.end - e.begin; k++, along++) {
          if (along != t.cds_begin && along + 1 != t.cds_end) continue;
          const std::size_t g = t.reverse ? e.end - 1 - k : e.begin + k;
          lo = std::min(lo, g);
          hi = std::max(hi, g);
        }
      }
      if (lo != SIZE_MAX) cds = Json::array({num(lo), num(hi)});
    }
    ts.push(Json::object({{"id", t.id}, {"gene", t.gene}, {"protein", t.protein_id}, {"strand", t.reverse ? "-" : "+"}, {"coding", t.coding},
                          {"c", pp.projection.text()}, {"p", pp.text()}, {"r", pp.rna}, {"consequences", terms}, {"impact", pp.impact()},
                          {"note", pp.note}, {"intronic", pp.projection.intronic}, {"outside", pp.projection.outside},
                          {"reference_protein_length", pp.reference.size()}, {"observed_protein_length", pp.observed.size()},
                          {"exons", exons}, {"cds", cds}}));
    if (ts.size() >= 60) break;
  }
  out.set("transcripts", ts);
  return out;
}

// {"vcf_file": path, "accessions": {"chr2": "NC_000002.12", …}, "centromeres": {"chr2": 93000000, …}} →
// {"records": [{line, id, chrom, pos, end, kind, alt, length, cipos, ciend, copy_number, hgvs}],
//  "adjacencies": [{a: {contig, position, strand}, b, inserted, alt_a, alt_b, hgvs}], "inversions": [hgvs],
//  "segments": [{contig, begin, end, copy_number, hgvs}], "reading"}.
BENCH_OP("nomenclature.structural") {
  detail::Reading reading;
  const std::string text = vcf_text(request, 64u << 20, reading);
  const auto svs = gnm::read_structural(text, &reading.problems);
  gnm::SvHgvsOptions o;
  const Json& accs = request["accessions"];
  const Json& cens = request["centromeres"];
  if (accs.is_object()) o.accession = [&accs](const std::string& c) { return accs.has(c) ? accs[c].text() : c; };
  if (cens.is_object())
    o.centromere = [&cens](const std::string& c) -> std::optional<std::int64_t> {
      if (!cens.has(c)) return std::nullopt;
      return static_cast<std::int64_t>(cens[c].as_number());
    };
  Json recs = Json::array();
  for (const auto& sv : svs) {
    if (recs.size() >= 2000) break;
    Json r = Json::object({{"line", sv.line}, {"id", sv.id}, {"chrom", sv.chrom}, {"pos", sv.pos}, {"end", sv.end}, {"kind", gnm::sv_kind_name(sv.kind)},
                           {"alt", sv.alt}, {"length", sv.length}, {"imprecise", sv.imprecise},
                           {"cipos", Json::array({sv.cipos[0], sv.cipos[1]})}, {"ciend", Json::array({sv.ciend[0], sv.ciend[1]})}});
    if (sv.copy_number) r.set("copy_number", *sv.copy_number);
    const std::string h = sv.kind == gnm::SvKind::breakend ? gnm::hgvs_of(sv.adjacency, o) : gnm::hgvs_of(sv, o);
    if (!h.empty()) r.set("hgvs", h);
    recs.push(r);
  }
  const auto adj = gnm::adjacencies(svs);
  Json as = Json::array();
  const auto bk = [](const gnm::Breakend& b) { return Json::object({{"contig", b.contig}, {"position", b.position + 1}, {"strand", std::string(1, b.strand)}}); };
  for (const auto& x : adj) {
    if (as.size() >= 2000) break;
    Json j = Json::object({{"a", bk(x.a)}, {"inserted", x.inserted}, {"single", x.single}, {"alt_a", gnm::breakend_alt(x, 'N')}});
    if (!x.single) {
      j.set("b", bk(x.b));
      j.set("alt_b", gnm::breakend_alt(gnm::reversed(x), 'N'));
      j.set("hgvs", gnm::hgvs_of(x, o));
    }
    as.push(j);
  }
  Json inv = Json::array();
  for (const auto& h : gnm::hgvs_of(adj, o)) inv.push(h);
  Json segs = Json::array();
  for (const auto& sgm : gnm::copy_number_segments(svs))
    segs.push(Json::object({{"contig", sgm.contig}, {"begin", sgm.begin + 1}, {"end", sgm.end}, {"copy_number", sgm.copy_number}, {"hgvs", gnm::hgvs_of(sgm, o)}}));
  return Json::object({{"records", recs}, {"adjacencies", as}, {"descriptions", inv}, {"segments", segs}, {"reading", reading.json()}});
}

}  // namespace workbench
