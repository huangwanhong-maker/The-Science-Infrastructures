// The operations of the karyotype view: the human ideogram for drawing, ISCN descriptions parsed
// (the named tree, the diagnostics, the canonical form), interpreted (every chromosome of every cell
// line as segments of the normal chromosomes, the copy number at every band, the losses and gains
// from where to where, the fusions, the detailed system), and results described (arr from
// copy-number segments, seq from junctions).
//
// Every computation is the domain's (genomics::karyotype); the engine reads the ideogram from the
// data library's cache (data/cache/ncbi-ideogram, the directory the server names) and keeps it.

#pragma once

#include "protocol.hpp"

#include <genomics/karyotype/bands.hpp>
#include <genomics/karyotype/describe.hpp>
#include <genomics/karyotype/format.hpp>
#include <genomics/karyotype/interpret.hpp>
#include <genomics/karyotype/iscn.hpp>

#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

namespace workbench {

namespace kary {

namespace gk = genomics::karyotype;

// The ideogram of an assembly at a level, read once from the cache directory the request names.
inline std::shared_ptr<const gk::Ideogram> ideogram(const Json& request, int level_default = 850) {
  const std::string dir = str(request, "dir");
  if (dir.empty()) throw Fault("the request names no ideogram directory (dir: data/cache/ncbi-ideogram)");
  const std::string assembly = str(request, "assembly", "GRCh38");
  const int level = static_cast<int>(integer(request, "level", level_default));
  if (level != 400 && level != 550 && level != 850) throw Fault("an ideogram's level is 400, 550 or 850 bands");
  std::string file;
  bool ucsc = false;
  if (assembly == "GRCh38") file = "ideogram_9606_GCF_000001305.16_" + std::to_string(level) + "_V1.tsv";
  else if (assembly == "GRCh37") {
    if (level == 850) { file = "ucsc-hg19-cytoBand.txt"; ucsc = true; }  // NCBI's GRCh37 850 level lacks densities
    else file = "ideogram_9606_GCF_000001305.13_" + std::to_string(level) + "_V1.tsv";
  } else {
    throw Fault("the ideograms are of GRCh38 and GRCh37");
  }
  static std::mutex lock;
  static std::map<std::string, std::shared_ptr<const gk::Ideogram>> kept;
  const std::string key = dir + "/" + file;
  std::lock_guard<std::mutex> hold(lock);
  if (auto it = kept.find(key); it != kept.end()) return it->second;
  std::ifstream f(key, std::ios::binary);
  if (!f) throw Fault("cannot read " + key + " (python data/acquisition/fetch.py ncbi-ideogram)");
  std::stringstream s;
  s << f.rdbuf();
  auto g = std::make_shared<gk::Ideogram>(ucsc ? gk::read_cytoband(s.str(), level, nullptr, assembly) : gk::read_ideogram(s.str(), level, nullptr, assembly));
  kept[key] = g;
  return g;
}

inline Json tree(const morphe::formal::ParseNode& n) {
  Json kids = Json::array();
  for (const auto& c : n.children) kids.push(tree(c));
  return Json::object({{"name", n.name}, {"begin", n.begin}, {"end", n.end}, {"children", std::move(kids)}});
}

inline Json diagnostics(const std::vector<gk::Diagnostic>& ds) {
  Json out = Json::array();
  for (const auto& d : ds)
    out.push(Json::object({{"level", std::string(gk::level_name(d.level))}, {"begin", d.begin}, {"end", d.end}, {"message", d.message}}));
  return out;
}

inline gk::Edition edition(const Json& request) { return str(request, "edition", "2020") == "2016" ? gk::Edition::iscn2016 : gk::Edition::iscn2020; }

inline std::string point_band(const gk::Point& p, const gk::Ideogram& g, bool lower) {
  if (p.chromosome < 0 || !g.has(p.chromosome)) return "?";
  if (p.kind == gk::Point::Kind::band && !p.written_position) return gk::chromosome_name(p.chromosome) + p.label;
  const std::uint64_t x = lower ? p.at : (p.at > 0 ? p.at - 1 : 0);
  return g.bands[g.band_at(p.chromosome, std::min(x, g.length[p.chromosome] - 1))].full_name();
}

inline Json chromosome(const gk::Chromosome& ch, const gk::Ideogram& g) {
  Json segs = Json::array();
  for (const auto& s : ch.segments) {
    if (s.kind != gk::Segment::Kind::reference) {
      segs.push(Json::object({{"kind", s.kind == gk::Segment::Kind::hsr ? "hsr" : s.kind == gk::Segment::Kind::bases ? "bases" : "unknown"},
                              {"bases", s.bases}}));
      continue;
    }
    const bool rev = s.reversed();
    const gk::Point& lo = rev ? s.to : s.from;
    const gk::Point& hi = rev ? s.from : s.to;
    segs.push(Json::object({{"kind", "reference"}, {"chrom", gk::chromosome_name(s.chromosome)}, {"start", s.start()}, {"end", s.end()},
                            {"reverse", rev}, {"from_band", point_band(s.from, g, !rev)}, {"to_band", point_band(s.to, g, rev)},
                            {"start_uncertain", Json::array({lo.lo, lo.hi})}, {"end_uncertain", Json::array({hi.lo, hi.hi})}}));
  }
  gk::LongFormOptions o;
  return Json::object({{"key", ch.key}, {"label", ch.label}, {"centromere_of", ch.centromere_of >= 0 ? Json(gk::chromosome_name(ch.centromere_of)) : Json()},
                       {"abnormal", ch.abnormal}, {"ring", ch.ring}, {"marker", ch.marker}, {"counted", ch.counted},
                       {"long_form", ch.abnormal && !ch.marker ? Json(gk::long_form(ch, o)) : Json()}, {"segments", std::move(segs)}});
}

inline Json change(const gk::Change& c) {
  const std::string n = gk::chromosome_name(c.chromosome);
  return Json::object({{"chrom", n}, {"start", c.start}, {"end", c.end}, {"copies", c.copies}, {"normal", c.normal},
                       {"inner", Json::array({c.inner_start, c.inner_end})}, {"outer", Json::array({c.outer_start, c.outer_end})},
                       {"from_band", n + c.from_band}, {"to_band", n + c.to_band}});
}

inline Json band_names(const std::vector<std::size_t>& bs, const gk::Ideogram& g) {
  Json out = Json::array();
  for (auto b : bs) out.push(g.bands[b].full_name());
  return out;
}

inline Json line(const gk::LineInterpretation& L, const gk::CellLine* cl, const gk::Ideogram& g) {
  Json chroms = Json::array(), changes = Json::array(), fusions = Json::array(), copies = Json::array();
  for (const auto& ch : L.chromosomes) chroms.push(chromosome(ch, g));
  for (const auto& c : L.changes) changes.push(change(c));
  for (const auto& f : L.fusions) {
    const auto end = [&](const gk::Point& p) {
      return p.kind == gk::Point::Kind::unknown ? Json() : Json::object({{"chrom", gk::chromosome_name(p.chromosome)}, {"at", p.at}, {"label", p.label}});
    };
    fusions.push(Json::object({{"chromosome", f.chromosome}, {"a", end(f.a)}, {"b", end(f.b)}, {"ring", f.ring}}));
  }
  // Every band of the ideogram: [min, max, normal], in the ideogram's order.
  for (const auto& b : L.bands) copies.push(Json::array({b.min, b.max, b.normal}));
  Json events = Json::array();
  for (const auto& e : L.events) {
    Json prods = Json::array();
    for (auto p : e.products) prods.push(p);
    events.push(Json::object({{"abnormality", e.abnormality}, {"text", L.abnormalities[e.abnormality].text}, {"products", std::move(prods)},
                              {"loss", band_names(e.lgf.loss, g)}, {"gain", band_names(e.lgf.gain, g)}, {"fusion", band_names(e.lgf.fusion, g)}}));
  }
  Json labels = Json::array();
  for (const auto& s : L.lgf.fusion_labels) labels.push(s);
  return Json::object({{"line", L.line}, {"cells", cl ? Json(cl->cells) : Json()}, {"composite", cl ? cl->composite : false}, {"counted", L.counted},
                       {"ploidy", L.ploidy}, {"male", L.male}, {"sex_known", L.sex_known}, {"balanced", L.balanced}, {"unknown", L.unknown},
                       {"chromosomes", std::move(chroms)}, {"changes", std::move(changes)}, {"fusions", std::move(fusions)}, {"band_copies", std::move(copies)},
                       {"events", std::move(events)},
                       {"lgf", Json::object({{"loss", band_names(L.lgf.loss, g)}, {"gain", band_names(L.lgf.gain, g)}, {"fusion", band_names(L.lgf.fusion, g)},
                                             {"fusion_labels", std::move(labels)}})},
                       {"detailed", cl ? Json(gk::detailed(L, *cl)) : Json()}, {"diagnostics", diagnostics(L.diagnostics)}});
}

inline Json parsed(const gk::Parsed& p) {
  Json expected = Json::array();
  for (const auto& e : p.syntax.expected) expected.push(e);
  return Json::object({{"accepted", p.syntax.accepted}, {"ok", p.ok()}, {"tree", p.syntax.accepted ? tree(p.syntax.tree) : Json()},
                       {"furthest", p.syntax.furthest}, {"expected", std::move(expected)}, {"error", p.errors()},
                       {"canonical", p.syntax.accepted ? Json(gk::format(p.description)) : Json()}, {"diagnostics", diagnostics(p.diagnostics)}});
}

inline std::string text_of(const Json& request) {
  const std::string text = str(request, "text");
  if (text.empty()) throw Fault("the request gives no ISCN text");
  if (text.size() > 8000) throw Fault("an ISCN text of at most 8000 characters");
  return text;
}

}  // namespace kary

// The ideogram for drawing: every chromosome's bands with their coordinates, stains and ISCN lengths.
// {"dir", "assembly": "GRCh38", "level": 850}
BENCH_OP("karyotype.ideogram") {
  const auto g = kary::ideogram(request);
  Json chroms = Json::array();
  for (int c = 0; c < kary::gk::chromosome_count; c++) {
    if (!g->has(c)) continue;
    Json bands = Json::array();
    for (std::size_t b = g->begin_of(c); b < g->end_of(c); b++) {
      const auto& x = g->bands[b];
      bands.push(Json::array({x.name(), x.start, x.end, std::string(kary::gk::stain_name(x.stain)), x.iscn_start, x.iscn_end}));
    }
    chroms.push(Json::object({{"name", kary::gk::chromosome_name(c)}, {"length", g->length[c]}, {"centromere", g->centromere[c]},
                              {"iscn_length", g->iscn_length[c]}, {"bands", std::move(bands)}}));
  }
  return Json::object({{"assembly", g->assembly}, {"level", g->level}, {"columns", Json::array({"name", "start", "end", "stain", "iscn_start", "iscn_end"})},
                       {"chromosomes", std::move(chroms)}});
}

// An ISCN text parsed: acceptance, the named tree, the furthest viable prefix and what could follow,
// the diagnostics of sense, the canonical form. {"text", "edition": "2020" | "2016"}
BENCH_OP("karyotype.parse") {
  const std::string text = kary::text_of(request);
  return kary::parsed(kary::gk::parse(text, kary::edition(request)));
}

// An ISCN text interpreted on the ideogram: the parse, and for every cell line and technique its
// chromosomes, changes, band copies, fusions, LGF and detailed system.
// {"text", "edition", "dir", "assembly", "level"}
BENCH_OP("karyotype.interpret") {
  const std::string text = kary::text_of(request);
  const auto p = kary::gk::parse(text, kary::edition(request));
  Json out = kary::parsed(p);
  const auto g = kary::ideogram(request);
  Json lines = Json::array(), techniques = Json::array();
  if (p.syntax.accepted) {
    const auto in = kary::gk::interpret(p.description, *g);
    for (std::size_t i = 0; i < in.lines.size(); i++) lines.push(kary::line(in.lines[i], &p.description.lines[i], *g));
    for (const auto& t : in.techniques) {
      Json changes = Json::array();
      for (const auto& c : t.changes) changes.push(kary::change(c));
      techniques.push(Json::object({{"technique", p.description.techniques[t.technique].name}, {"build", p.description.techniques[t.technique].build},
                                    {"assembly", t.assembly}, {"changes", std::move(changes)},
                                    {"line", t.has_line ? kary::line(t.line, nullptr, *g) : Json()}, {"diagnostics", kary::diagnostics(t.diagnostics)}}));
    }
  }
  out.set("assembly", g->assembly);
  out.set("level", g->level);
  out.set("lines", std::move(lines));
  out.set("techniques", std::move(techniques));
  return out;
}

// Results described in ISCN: copy-number segments as arr, junctions as seq.
// {"segments": [{"chrom", "start", "end", "copies", "copies_hi", "proportion", "hmz", "inheritance"}],
//  "adjacencies": [{"a", "break_a", "a_left", "b", "break_b", "b_right", "bases"}], "male", "build", "dir", "assembly"}
BENCH_OP("karyotype.describe") {
  const auto g = kary::ideogram(request);
  std::vector<kary::gk::CopySegment> segs;
  if (request["segments"].is_array())
    for (const auto& s : request["segments"].as_array()) {
      kary::gk::CopySegment c;
      c.chromosome = kary::gk::chromosome_index(str(s, "chrom"));
      if (c.chromosome < 0) throw Fault("a segment's chromosome is 1-22, X or Y");
      c.start = static_cast<std::uint64_t>(integer(s, "start", 0));
      c.end = static_cast<std::uint64_t>(integer(s, "end", 0));
      c.copies = static_cast<int>(integer(s, "copies", 2));
      c.copies_hi = static_cast<int>(integer(s, "copies_hi", -1));
      c.proportion = num(s, "proportion", -1);
      c.homozygous = flag(s, "hmz", false);
      c.inheritance = str(s, "inheritance");
      segs.push_back(c);
    }
  std::vector<kary::gk::Adjacency> adjs;
  if (request["adjacencies"].is_array())
    for (const auto& a : request["adjacencies"].as_array()) {
      kary::gk::Adjacency x;
      x.chromosome_a = kary::gk::chromosome_index(str(a, "a"));
      x.chromosome_b = kary::gk::chromosome_index(str(a, "b"));
      if (x.chromosome_a < 0 || x.chromosome_b < 0) throw Fault("an adjacency's chromosomes are 1-22, X or Y");
      x.break_a = static_cast<std::uint64_t>(integer(a, "break_a", 0));
      x.break_b = static_cast<std::uint64_t>(integer(a, "break_b", 0));
      x.a_left = flag(a, "a_left", true);
      x.b_right = flag(a, "b_right", true);
      x.bases = str(a, "bases");
      adjs.push_back(x);
    }
  const std::string build = str(request, "build", g->assembly);
  kary::gk::ArrOptions ao;
  ao.build = build;
  if (request["male"].is_bool()) ao.male = request["male"].as_bool();
  kary::gk::SeqOptions so;
  so.build = build;
  const auto seq = kary::gk::describe_seq(adjs, *g, so);
  Json ders = Json::array(), junctions = Json::array();
  for (const auto& d : seq.derivatives) ders.push(kary::chromosome(d, *g));
  for (const auto& j : seq.junctions) junctions.push(j);
  return Json::object({{"arr", segs.empty() && !request["male"].is_bool() ? Json() : Json(kary::gk::describe_arr(segs, *g, ao))},
                       {"seq", Json::object({{"kind", seq.kind}, {"short", seq.short_form}, {"detailed", seq.detailed}, {"junctions", std::move(junctions)},
                                             {"derivatives", std::move(ders)}, {"diagnostics", kary::diagnostics(seq.diagnostics)}})}});
}

// Examples the page may offer, each with its source.
BENCH_OP("karyotype.examples") {
  struct E { const char* text; const char* kind; const char* source; const char* edition; };
  static const E examples[] = {
      {"47,XX,+21", "constitutional, numerical", "trisomy 21; ISCN 2009 p. 57", "2020"},
      {"45,X", "constitutional, sex chromosomes", "Turner syndrome; ISCN 2009 p. 56", "2020"},
      {"46,XY,t(9;22)(q34;q11.2)", "acquired, balanced translocation", "the Philadelphia translocation; ISCN 2009 pp. 64-66", "2020"},
      {"46,XY,der(6)t(3;6)(p13;q21)", "unbalanced derivative", "Conlin et al. 2025, Clin. Chem. 71, Table 1", "2020"},
      {"46,XX,del(5)(q13q33)", "interstitial deletion", "ISCN 2009 p. 61", "2020"},
      {"46,XX,ins(5;2)(p14;q32q22)", "insertion between chromosomes", "ISCN 2020 §9.2.9 (Miller and Madan 2021)", "2020"},
      {"45,XX,der(13;21)(q10;q10)", "Robertsonian translocation", "ISCN 2009 p. 81", "2020"},
      {"46,X,i(X)(q10)", "isochromosome", "ISCN 2020 §9.2.11 (Miller and Madan 2021)", "2020"},
      {"46,XX,r(7)(p22q36)", "ring chromosome", "ISCN 2009 p. 75", "2020"},
      {"46,XX,der(1)t(1;3)(p32;q21)t(1;11)(q25;q13)", "derivative of two translocations", "ISCN 2009 p. 63", "2020"},
      {"46,XX,t(2;7;5)(p21;q22;q23)", "three-way translocation", "ISCN 2009 p. 78", "2020"},
      {"46,XY,rec(2)dup(2q)inv(2)(p23q32)dpat", "recombinant of a parental inversion", "Conlin et al. 2025, Table 1", "2020"},
      {"47,XY,+12[10]/49,idem,+18,+19[5]/46,XY[20]", "clonal evolution in cancer", "after CytoGPS's example (Abrams et al. 2019)", "2020"},
      {"mos 47,XXY[10]/46,XY[20]", "mosaic", "ISCN 2009 p. 56", "2020"},
      {"49~55,XY,-3,+6,-7,+10,+17,+18,+21c,+21,+r,+mar1,+mar2x2,3~9dmin[cp10]/47,XY,+21c[10]", "composite hyperdiploid karyotype with constitutional trisomy 21",
       "GenQA's ISCN assessment (Hastings, Tabiner and Deans, ESHG 2022)", "2020"},
      {"arr[GRCh38] 7q11.23(73,265,388_74,726,130)x1", "microarray, Williams syndrome region", "Conlin et al. 2025", "2020"},
      {"arr[GRCh38] 1q21.1q21.2(147,105,904_147,917,509)x1", "microarray deletion", "Conlin et al. 2025, Table 1", "2020"},
      {"seq[GRCh38] t(9;10)(q34.3;q24.2)", "sequencing, balanced translocation", "Conlin et al. 2025", "2020"},
      {"seq[GRCh37/hg19] t(1;8)(1pter→1q23.3(164,789,896)::8q12.2(61,760,650)→8qter;8pter→8q12.2(61,760,649)::1q23.3(164,789,897)→1qter)dn",
       "sequencing, detailed system (GRCh37: use the GRCh37 ideogram)", "Ordulu et al. 2014, Am. J. Hum. Genet. 94, 695-709", "2020"},
  };
  Json out = Json::array();
  for (const auto& e : examples)
    out.push(Json::object({{"text", e.text}, {"kind", e.kind}, {"source", e.source}, {"edition", e.edition}}));
  return Json::object({{"examples", std::move(out)}});
}

}  // namespace workbench
