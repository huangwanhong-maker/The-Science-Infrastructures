// The operations of structural variants and copy number: a simulated individual made on a real
// reference, read, mapped, and everything the structural module detects, with its truth and the
// evaluation, in the shapes a page draws.
//
// - `structural.sample` makes the individual (genomics::simulate, structural.hpp: the variants
//   given, or drawn at random, on top of a diploid's small variants), draws its paired reads
//   (simulate/pairs.hpp), maps both mates with the worker's short-read mapper of the reference,
//   runs genomics::structural::detect (library, pairs, evidence, breakpoints, events and their
//   genotypes, depth, segments, aneuploidy, small variants, B-allele fractions, runs of
//   homozygosity, allelic states) and compares the calls with the truth (Truvari's matching, copy
//   number bin by bin, LOH base by base). The reference is a genome of the data library ("genome",
//   as every view names one), or several ("genomes", one contig each joined), or regions of them
//   ("regions": parts of their contigs taken as the sample's contigs, their offsets kept so that
//   positions are written in the chromosome's own numbering). The run is kept by the worker, and
//   `structural.window` serves the mates of a stretch of it (or around a call) for a browser.
//   Every call is described (genomics::structural::describe): its VCF records, its HGVS g. form on
//   the chromosome's accession, and, with an ideogram ("ideogram": {"file": the data library's
//   ncbi-ideogram TSV}), its ISCN form and the sample's ISCN arr and seq descriptions; a contig is
//   placed on a human chromosome by its RefSeq accession (NC_000001 … NC_000024).
//
// The sample is simulated and labelled so in every answer ("simulated": true); only the reference
// is real. Every computation is the domains'; the engine gathers, times and reports.

#pragma once

#include "genomes.hpp"
#include "shared.hpp"
#include "shared_mapping.hpp"

#include <genomics/simulate/pairs.hpp>
#include <genomics/simulate/structural.hpp>
#include <genomics/karyotype/bands.hpp>
#include <genomics/structural/describe.hpp>
#include <genomics/structural/detect.hpp>
#include <genomics/structural/evaluate.hpp>
#include <genomics/structural/vcf.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace workbench {

namespace gsv = genomics::structural;

namespace detail {

// The sample's contigs: each a region of a contig of a reference genome.
struct SvContig {
  std::string name, accession;
  std::size_t offset = 0;  // its first base's position on the chromosome (0-based)
  seq::Codes codes;
};

struct StructuralRun {
  std::vector<SvContig> contigs;
  std::vector<seq::Codes> reference;
  std::vector<std::string> names;
  sim::StructuralSample sample;
  std::vector<sim::ReadPair> pairs;
  gsv::MateMappings mates;
  gsv::Detection detection;
  std::vector<gsv::TruthEvent> truth;
  gsv::SvComparison comparison;
  gsv::CopyConcordance copy;
  gsv::OverlapTally loh;
  double seconds_simulate = 0, seconds_index = 0, seconds_map = 0, seconds_detect = 0, seconds_compare = 0;
};

inline sim::SvType sv_type_of(const std::string& t) {
  if (t == "DEL") return sim::SvType::deletion;
  if (t == "DUP" || t == "DUP:TANDEM") return sim::SvType::tandem_duplication;
  if (t == "DUP:DISPERSED") return sim::SvType::dispersed_duplication;
  if (t == "INV") return sim::SvType::inversion;
  if (t == "INS") return sim::SvType::insertion;
  if (t == "TRA") return sim::SvType::translocation;
  if (t == "ANEUPLOIDY") return sim::SvType::aneuploidy;
  if (t == "LOH") return sim::SvType::loss_of_heterozygosity;
  throw Fault("a structural variant is DEL, DUP, DUP:DISPERSED, INV, INS, TRA, ANEUPLOIDY or LOH, not " + t);
}

// The variants a request gives: {"type", "contig", "begin", "end", "zygosity": "het"|"hom",
// "haplotype": 0|1, "copies", "length", "sequence", "partner", "partner_position", "inverted"},
// positions on the chromosome's numbering (0-based) when the contigs are regions.
inline std::vector<sim::StructuralVariant> variants_of(const Json& list, const std::vector<SvContig>& contigs) {
  std::vector<sim::StructuralVariant> out;
  for (const auto& v : list.as_array()) {
    sim::StructuralVariant x;
    x.type = sv_type_of(str(v, "type", "DEL"));
    x.contig = static_cast<std::uint32_t>(bounded(v, "contig", 0, 0, static_cast<long long>(contigs.size()) - 1));
    const auto local = [&](std::uint32_t c, long long p) -> std::size_t {
      const long long q = p - static_cast<long long>(contigs[c].offset);
      if (q < 0 || q > static_cast<long long>(contigs[c].codes.size())) throw Fault("a position outside its contig's region");
      return static_cast<std::size_t>(q);
    };
    x.begin = x.type == sim::SvType::aneuploidy ? 0 : local(x.contig, integer(v, "begin", 0));
    x.end = x.type == sim::SvType::aneuploidy ? contigs[x.contig].codes.size() : local(x.contig, integer(v, "end", static_cast<long long>(contigs[x.contig].offset + x.begin)));
    if (x.type == sim::SvType::insertion || x.type == sim::SvType::translocation) x.end = x.begin;
    x.copies = static_cast<int>(integer(v, "copies", 1));
    x.insert_length = static_cast<std::size_t>(bounded(v, "length", 100, 1, 100000));
    x.sequence = str(v, "sequence");
    x.inverted = flag(v, "inverted", false);
    if (x.type == sim::SvType::dispersed_duplication || x.type == sim::SvType::translocation) {
      x.partner = static_cast<std::uint32_t>(bounded(v, "partner", 0, 0, static_cast<long long>(contigs.size()) - 1));
      x.partner_position = local(x.partner, integer(v, "partner_position", 0));
    }
    const int h = static_cast<int>(bounded(v, "haplotype", 0, 0, 1));
    x.on = str(v, "zygosity", "het") == "hom" && x.type != sim::SvType::loss_of_heterozygosity && x.type != sim::SvType::aneuploidy ? std::array<bool, 2>{true, true}
                                                                                                                                : std::array<bool, 2>{h == 0, h == 1};
    out.push_back(x);
  }
  return out;
}

inline sim::StructuralSpec random_spec(const Json& r) {
  sim::StructuralSpec s;
  const auto n = [&](const char* k, std::size_t d) { return static_cast<std::size_t>(bounded(r, k, static_cast<long long>(d), 0, 200)); };
  s.deletions = n("deletions", 8);
  s.tandem_duplications = n("tandem_duplications", 3);
  s.dispersed_duplications = n("dispersed_duplications", 1);
  s.inversions = n("inversions", 2);
  s.insertions = n("insertions", 3);
  s.translocations = n("translocations", 0);
  s.gains = n("gains", 1);
  s.losses = n("losses", 1);
  s.losses_of_heterozygosity = n("losses_of_heterozygosity", 1);
  s.homozygous = num(r, "homozygous", 0.3);
  return s;
}

inline std::string run_key(const Json& r) {
  Json k = Json::object({{"genome", r["genome"]}, {"genomes", r["genomes"]}, {"regions", r["regions"]}, {"variants", r["variants"]}, {"random", r["random"]},
                         {"coverage", r["coverage"]}, {"seed", r["seed"]}, {"fragment", r["fragment"]}, {"read_length", r["read_length"]}, {"bin", r["bin"]}, {"snv", r["snv"]}});
  return "structural:" + k.dump();
}

// The genomes of the request (one, or several joined) and the regions of their contigs taken.
inline std::vector<SvContig> sample_contigs(Worker& worker, const Json& request) {
  std::vector<std::shared_ptr<Genome>> genomes;
  if (request["genomes"].is_array())
    for (const auto& spec : request["genomes"].as_array()) genomes.push_back(worker.genome(Json::object({{"genome", spec}})));
  else genomes.push_back(worker.genome(request));
  std::vector<const Contig*> all;
  for (const auto& g : genomes)
    for (const auto& c : g->contigs) all.push_back(&c);
  std::vector<SvContig> out;
  const auto take = [&](const Contig& c, std::size_t b, std::size_t e) {
    SvContig x;
    x.name = c.name;
    x.accession = c.accession.empty() ? c.name : c.accession;
    x.offset = b;
    x.codes.assign(c.codes.begin() + static_cast<std::ptrdiff_t>(b), c.codes.begin() + static_cast<std::ptrdiff_t>(e));
    out.push_back(std::move(x));
  };
  if (request["regions"].is_array()) {
    for (const auto& r : request["regions"].as_array()) {
      const auto k = static_cast<std::size_t>(bounded(r, "contig", 0, 0, static_cast<long long>(all.size()) - 1));
      const auto b = static_cast<std::size_t>(bounded(r, "begin", 0, 0, static_cast<long long>(all[k]->codes.size())));
      const auto e = static_cast<std::size_t>(bounded(r, "end", static_cast<long long>(all[k]->codes.size()), static_cast<long long>(b) + 1, static_cast<long long>(all[k]->codes.size())));
      take(*all[k], b, e);
    }
  } else {
    for (const Contig* c : all) take(*c, 0, c->codes.size());
  }
  std::size_t total = 0;
  for (const auto& c : out) total += c.codes.size();
  if (total > 12'000'000) throw Fault("a structural sample is made on up to 12 Mb of reference: take regions of a chromosome");
  return out;
}

inline std::shared_ptr<StructuralRun> structural_run(Worker& worker, const Json& request) {
  // Kept by the worker under the first genome's key and the request's.
  const auto g = request["genomes"].is_array() ? worker.genome(Json::object({{"genome", request["genomes"][0]}})) : worker.genome(request);
  return worker.derived<StructuralRun>(*g, run_key(request), [&] {
    using clock = std::chrono::steady_clock;
    const auto since = [](clock::time_point t) { return std::chrono::duration<double>(clock::now() - t).count(); };
    StructuralRun run;
    run.contigs = sample_contigs(worker, request);
    for (const auto& c : run.contigs) {
      run.reference.push_back(c.codes);
      run.names.push_back(c.accession);
    }
    const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
    const double coverage = num(request, "coverage", 30);
    if (!(coverage >= 2 && coverage <= 100)) throw Fault("a coverage is between 2 and 100");
    auto t = clock::now();
    std::vector<sim::StructuralVariant> svs;
    if (request["variants"].is_array()) svs = variants_of(request["variants"], run.contigs);
    else svs = sim::random_structural(run.reference, random_spec(request["random"]), seed);
    sim::MutationRates rates;
    rates.snv = num(request, "snv", 1e-3);
    rates.indel = num(request, "indel", 1e-4);
    try {
      run.sample = sim::make_structural_sample(run.reference, svs, rates, seed, run.names);
    } catch (const std::invalid_argument& e) {
      throw Fault(e.what());
    }
    sim::PairModel pm;
    pm.read.length = static_cast<std::size_t>(bounded(request, "read_length", 150, 50, 300));
    const Json& f = request["fragment"];
    pm.fragment.mean = num(f, "mean", 450);
    pm.fragment.sd = num(f, "sd", 45);
    pm.fragment.gc_bias = num(f, "gc_bias", 0.5);
    pm.fragment.gc_optimum = num(f, "gc_optimum", 0.45);
    pm.fragment.gc_width = num(f, "gc_width", 0.08);
    if (!(pm.fragment.mean >= static_cast<double>(pm.read.length) + 8 && pm.fragment.sd > 0 && pm.fragment.mean <= 2000)) throw Fault("fragments longer than the reads, of at most 2000 bases");
    const std::size_t n = sim::pairs_for_coverage(run.sample.chromosomes, coverage, pm.read.length);
    // Pairs in chunks on the worker's threads, each pair from its own stream.
    constexpr std::size_t CHUNK = 8192;
    const std::size_t chunks = (n + CHUNK - 1) / CHUNK;
    std::vector<std::vector<sim::ReadPair>> parts(chunks);
    calx::scheduler::parallel_for(worker.pool, chunks, 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t c = b; c < e; c++) parts[c] = sim::read_pairs(run.sample.chromosomes, std::min(CHUNK, n - c * CHUNK), pm, seed + 1, c * CHUNK);
    });
    for (auto& p : parts) std::move(p.begin(), p.end(), std::back_inserter(run.pairs));
    run.seconds_simulate = since(t);
    t = clock::now();
    const gmap::Mapper mapper(worker.pool, run.reference, gmap::MapperOptions::short_reads());
    run.seconds_index = since(t);
    t = clock::now();
    run.mates = gsv::map_pairs(worker.pool, mapper, run.pairs);
    run.seconds_map = since(t);
    t = clock::now();
    gsv::DetectOptions o;
    o.depth.bin = static_cast<std::size_t>(bounded(request, "bin", 1000, 100, 100000));
    run.detection = gsv::detect(worker.pool, run.reference, run.mates, run.pairs, o);
    run.seconds_detect = since(t);
    t = clock::now();
    run.truth = gsv::truth_events(run.sample.truth);
    run.comparison = gsv::compare_sv(run.reference, run.detection.calls, run.truth);
    run.copy = gsv::copy_concordance(run.detection.depth, run.detection.segments, [&](std::uint32_t c, std::size_t p) { const auto k = run.sample.copies_at(c, p); return k[0] + k[1]; });
    std::vector<std::array<std::size_t, 3>> called, truth;
    for (const auto& v : run.sample.truth)
      if (v.type == sim::SvType::loss_of_heterozygosity) truth.push_back({v.contig, v.begin, v.end});
    for (const auto& a : run.detection.allelic)
      if (a.state == gsv::AllelicState::copy_neutral_loh) called.push_back({a.contig, a.begin, a.end});
    run.loh = gsv::base_overlap(called, truth);
    run.seconds_compare = since(t);
    return run;
  });
}

inline double r4(double x) { return std::isfinite(x) ? std::round(x * 1e4) / 1e4 : 0; }

// A contig placed on its chromosome: the karyotype index of a human RefSeq chromosome
// (NC_000001 … NC_000022, NC_000023 X, NC_000024 Y), or −1.
inline std::vector<gsv::Placement> placements_of(const StructuralRun& run) {
  std::vector<gsv::Placement> out;
  for (const auto& c : run.contigs) {
    int chromosome = -1;
    if (c.accession.rfind("NC_0000", 0) == 0 && c.accession.size() >= 9) {
      const int n = std::atoi(c.accession.substr(7, 2).c_str());
      if (n >= 1 && n <= 24) chromosome = n - 1;
    }
    out.push_back({c.accession, chromosome, c.offset});
  }
  return out;
}

inline Json tally_json(const gsv::Tally& t) {
  return Json::object({{"tp", t.tp_call}, {"fp", t.fp}, {"fn", t.fn}, {"precision", r4(t.precision())}, {"recall", r4(t.recall())}, {"f1", r4(t.f1())}, {"genotype", r4(t.genotype_concordance())}});
}

inline Json breakend_json(const StructuralRun& run, const gsv::Breakend& e) {
  return Json::object({{"contig", e.contig}, {"position", run.contigs[e.contig].offset + e.position}, {"plus", e.plus}});
}

inline Json call_json(const StructuralRun& run, const gsv::SvCall& c, std::size_t index, const std::vector<genomics::variant::VcfRecord>& records) {
  const auto& off = run.contigs;
  Json junctions = Json::array();
  for (const auto& j : c.junctions)
    junctions.push(Json::object({{"a", breakend_json(run, j.a)}, {"b", breakend_json(run, j.b)}, {"inserted", j.inserted.size() > 200 ? j.inserted.substr(0, 200) + "…" : j.inserted},
                                 {"homology", gsv::homology(run.reference, j).length()}}));
  // The evidence of its breakpoints, for drawing: pairs as their two mates, splits as their parts.
  Json pairs = Json::array(), splits = Json::array();
  for (const auto b : c.breakpoints)
    for (const auto e : run.detection.breakpoints[b].evidence) {
      const auto& ev = run.detection.evidence[e];
      if (ev.source == gsv::Source::pair && pairs.size() < 80) {
        const auto& rec = run.detection.pairs[ev.read];
        Json mates = Json::array();
        for (int m = 0; m < 2; m++) {
          const auto* x = rec.mate[m];
          mates.push(Json::array({x->contig, off[x->contig].offset + x->target_begin, off[x->contig].offset + x->target_end, x->reverse, x->mapq}));
        }
        pairs.push(mates);
      } else if (ev.source != gsv::Source::pair && splits.size() < 80) {
        Json parts = Json::array();
        for (const auto& m : run.mates[ev.read])
          if (m.primary) parts.push(Json::array({m.contig, off[m.contig].offset + m.target_begin, off[m.contig].offset + m.target_end, m.reverse, m.mapq, m.clip_left(), m.clip_right()}));
        splits.push(Json::object({{"source", gsv::source_name(ev.source)}, {"parts", parts}}));
      }
    }
  int truth = -1;
  std::size_t error = 0;
  for (const auto& m : run.comparison.matches)
    if (m.call == index) { truth = static_cast<int>(m.truth); error = m.breakpoint_error; }
  Json vcf = Json::array();
  for (const auto& r : records)
    if (r.id == "sv" + std::to_string(index + 1) || r.id.rfind("sv" + std::to_string(index + 1) + "_", 0) == 0) vcf.push(genomics::variant::vcf_line(r));
  return Json::object({{"id", "sv" + std::to_string(index + 1)}, {"kind", gsv::sv_kind_name(c.kind)}, {"contig", c.contig}, {"begin", off[c.contig].offset + c.begin},
                       {"end", off[c.contig].offset + c.end}, {"partner", c.partner}, {"partner_position", c.kind == gsv::SvKind::translocation || c.kind == gsv::SvKind::dispersed_duplication || c.kind == gsv::SvKind::breakend ? Json(off[c.partner].offset + c.partner_position) : Json()},
                       {"inverted", c.inverted}, {"length", c.length}, {"precise", c.precise},
                       {"ci_begin", Json::array({off[c.contig].offset + c.ci_begin[0], off[c.contig].offset + c.ci_begin[1]})},
                       {"ci_end", Json::array({off[c.contig].offset + c.ci_end[0], off[c.contig].offset + c.ci_end[1]})},
                       {"support", Json::object({{"pairs", c.pairs}, {"splits", c.splits}, {"clips", c.clips}, {"anchored", c.anchored}})},
                       {"copies", c.copies >= 0 ? Json(r4(c.copies)) : Json()},
                       {"genotype", c.genotype.genotype == 2 ? "1/1" : c.genotype.genotype == 1 ? "0/1" : "0/0"}, {"gq", c.genotype.gq}, {"qual", std::round(c.qual * 10) / 10},
                       {"reference_support", r4(c.genotype.reference)}, {"alternate_support", r4(c.genotype.alternate)}, {"passed", c.passed()}, {"filters", c.filters},
                       {"sequence", c.sequence.size() > 300 ? c.sequence.substr(0, 300) + "…" : c.sequence}, {"resolved", c.resolved}, {"junctions", junctions},
                       {"evidence", Json::object({{"pairs", pairs}, {"splits", splits}})}, {"truth", truth >= 0 ? Json(truth) : Json()},
                       {"breakpoint_error", truth >= 0 && c.precise ? Json(error) : Json()}, {"vcf", vcf}});
}

}  // namespace detail

// A simulated individual on a real reference, read, mapped and detected; everything against its truth.
BENCH_OP("structural.sample") {
  const auto run = detail::structural_run(worker, request);
  const auto& d = run->detection;
  const auto& off = run->contigs;
  Json contigs = Json::array();
  for (const auto& c : off) contigs.push(Json::object({{"name", c.name}, {"accession", c.accession}, {"offset", c.offset}, {"length", c.codes.size()}}));
  // The library: its histogram of templates by 5 bases, from the survival function.
  Json histogram = Json::array();
  for (std::size_t x = 0; x + 5 < d.library.survival.size(); x += 5) {
    const double p = d.library.survival[x] - d.library.survival[x + 5];
    if (p > 0) histogram.push(Json::array({x, detail::r4(p)}));
  }
  Json pairs = Json::object();
  for (int k = 0; k < 9; k++) pairs.set(gsv::pair_class_name(static_cast<gsv::PairClass>(k)), d.counts.by_class[k]);
  // The depth in columns.
  std::vector<double> bc, bb, bn, be, br, bm, bg;
  for (const auto& b : d.depth.bins) {
    bc.push_back(b.contig);
    bb.push_back(static_cast<double>(off[b.contig].offset + b.begin));
    bn.push_back(b.count);
    be.push_back(detail::r4(b.exposure));
    br.push_back(detail::r4(2 * b.ratio()));
    bm.push_back(b.mappable ? 1 : 0);
    bg.push_back(detail::r4(b.gc));
  }
  Json segments = Json::array();
  for (const auto& g : d.segments)
    segments.push(Json::object({{"contig", g.contig}, {"begin", off[g.contig].offset + g.begin}, {"end", off[g.contig].offset + g.end}, {"copies", g.copies}, {"estimate", detail::r4(g.estimate())},
                                {"quality", std::round(g.quality * 10) / 10}, {"bins", g.bins}}));
  Json aneuploidy = Json::array();
  for (const auto& a : d.aneuploidy) aneuploidy.push(Json::object({{"contig", a.contig}, {"copies", a.copies}, {"fraction", detail::r4(a.fraction)}}));
  std::vector<double> sc, sp, sb, sg, sd;
  for (const auto& x : d.sites) {
    sc.push_back(x.contig);
    sp.push_back(static_cast<double>(off[x.contig].offset + x.position));
    sb.push_back(detail::r4(x.baf()));
    sg.push_back(x.genotype);
    sd.push_back(static_cast<double>(x.reference + x.alternate));
  }
  Json runs = Json::array();
  for (const auto& r : d.runs) runs.push(Json::object({{"contig", r.contig}, {"begin", off[r.contig].offset + r.begin}, {"end", off[r.contig].offset + r.end}, {"sites", r.sites}, {"quality", std::round(r.quality * 10) / 10}}));
  Json allelic = Json::array();
  for (const auto& a : d.allelic)
    allelic.push(Json::object({{"contig", a.contig}, {"begin", off[a.contig].offset + a.begin}, {"end", off[a.contig].offset + a.end}, {"state", gsv::allelic_state_name(a.state)}, {"copies", a.copies},
                               {"minor", a.minor}, {"sites", a.sites}, {"heterozygous", a.heterozygous}}));
  // VCF of the calls, in the chromosomes' numbering when the contigs are regions: records written on
  // the regions, positions shifted.
  auto records = gsv::sv_records(run->reference, run->names, d.calls);
  for (auto& r : records)
    for (std::size_t k = 0; k < off.size(); k++)
      if (r.chrom == run->names[k]) r.position += off[k].offset;
  std::vector<std::size_t> lengths;
  for (const auto& c : off) lengths.push_back(c.offset + c.codes.size());
  const std::string vcf = genomics::variant::write_vcf(gsv::sv_header(run->names, lengths, "SIMULATED"), records);
  // The descriptions: HGVS always, ISCN with an ideogram.
  std::unique_ptr<genomics::karyotype::Ideogram> ideogram;
  if (request["ideogram"].is_object() && request["ideogram"].has("file"))
    ideogram = std::make_unique<genomics::karyotype::Ideogram>(genomics::karyotype::read_ideogram(detail::read_text(request["ideogram"]["file"].text()), 850));
  const auto described = gsv::describe_sample(run->reference, detail::placements_of(*run), d.calls, d.segments, d.allelic, ideogram.get());
  Json calls = Json::array();
  for (std::size_t i = 0; i < d.calls.size() && i < 2000; i++) {
    Json c = detail::call_json(*run, d.calls[i], i, records);
    c.set("hgvs", described.calls[i].hgvs);
    c.set("iscn", described.calls[i].iscn);
    calls.push(std::move(c));
  }
  Json descriptions = Json::object({{"arr", described.arr}, {"seq", described.seq.short_form}, {"seq_detailed", described.seq.detailed},
                                    {"seq_junctions", Json(described.seq.junctions)}, {"ideogram", ideogram != nullptr}});
  // The truth.
  Json truth = Json::array();
  for (const auto& v : run->sample.truth)
    truth.push(Json::object({{"type", sim::sv_type_name(v.type)}, {"contig", v.contig}, {"begin", off[v.contig].offset + v.begin}, {"end", off[v.contig].offset + v.end},
                             {"partner", v.partner}, {"partner_position", v.type == sim::SvType::translocation || v.type == sim::SvType::dispersed_duplication ? Json(off[v.partner].offset + v.partner_position) : Json()},
                             {"inverted", v.inverted}, {"copies", v.copies}, {"length", v.length()}, {"zygosity", v.heterozygous() ? "het" : "hom"}, {"haplotypes", Json::array({v.on[0], v.on[1]})},
                             {"sequence", v.sequence.size() > 300 ? v.sequence.substr(0, 300) + "…" : v.sequence}}));
  Json copy_truth = Json::array();
  for (const auto& s : run->sample.copy_number)
    copy_truth.push(Json::array({s.contig, off[s.contig].offset + s.begin, off[s.contig].offset + s.end, s.copies[0], s.copies[1]}));
  Json chromosomes = Json::array();
  for (const auto& ch : run->sample.chromosomes) {
    Json pieces = Json::array();
    for (const auto& p : ch.pieces)
      pieces.push(p.novel ? Json::array({"novel", p.end - p.begin}) : Json::array({p.contig, off[p.contig].offset + p.begin, off[p.contig].offset + p.end, p.reverse, p.haplotype}));
    chromosomes.push(Json::object({{"name", ch.name}, {"length", ch.sequence.size()}, {"pieces", pieces}}));
  }
  const auto& cmp = run->comparison;
  Json by_kind = Json::object();
  for (const auto& [k, t] : cmp.by_kind) by_kind.set(k, detail::tally_json(t));
  Json by_size = Json::array();
  for (std::size_t k = 0; k < 5; k++) by_size.push(Json::object({{"bin", gsv::size_bin_name(k)}, {"tally", detail::tally_json(cmp.by_size[k])}}));
  Json confusion = Json::array();
  for (const auto& [tc, n] : run->copy.confusion) confusion.push(Json::array({tc.first, tc.second, n}));
  std::vector<double> errors(cmp.breakpoint_errors.begin(), cmp.breakpoint_errors.end());
  Json evaluation = Json::object({{"overall", detail::tally_json(cmp.overall)}, {"by_kind", by_kind}, {"by_size", by_size}, {"breakpoint_errors", Json(errors)},
                                  {"false_positives", Json(std::vector<double>(cmp.false_positives.begin(), cmp.false_positives.end()))},
                                  {"false_negatives", Json(std::vector<double>(cmp.false_negatives.begin(), cmp.false_negatives.end()))},
                                  {"copy_number", Json::object({{"bins", run->copy.bins}, {"agree", run->copy.agree}, {"fraction", detail::r4(run->copy.fraction())}, {"confusion", confusion}})},
                                  {"loh", Json::object({{"truth", run->loh.truth}, {"called", run->loh.called}, {"recall", detail::r4(run->loh.recall())}, {"precision", detail::r4(run->loh.precision())}})}});
  const auto& L = d.library;
  return Json::object({
      {"simulated", true},
      {"contigs", contigs},
      {"pairs_simulated", run->pairs.size()},
      {"library", Json::object({{"median", L.median}, {"mad", L.mad}, {"sd", detail::r4(L.sd())}, {"normal", Json::array({L.min_normal, L.max_normal})},
                                {"cutoffs", Json::array({L.min_cutoff, L.max_cutoff})}, {"read_length", L.read_length}, {"pairs", L.pairs}, {"histogram", histogram}})},
      {"pairs", pairs},
      {"depth", Json::object({{"bin", d.depth.bins.empty() ? 0 : d.depth.bins.front().end - d.depth.bins.front().begin}, {"fragment", d.depth.fragment}, {"contig", Json(bc)}, {"begin", Json(bb)},
                              {"count", Json(bn)}, {"exposure", Json(be)}, {"copies", Json(br)}, {"mappable", Json(bm)}, {"gc", Json(bg)}, {"gc_rate", numbers(d.depth.rate, 6)}})},
      {"segments", segments},
      {"aneuploidy", aneuploidy},
      {"baf", Json::object({{"contig", Json(sc)}, {"position", Json(sp)}, {"baf", Json(sb)}, {"genotype", Json(sg)}, {"depth", Json(sd)}})},
      {"runs", runs},
      {"allelic", allelic},
      {"calls", calls},
      {"descriptions", descriptions},
      {"vcf", vcf},
      {"truth", Json::object({{"variants", truth}, {"copy_number", copy_truth}, {"chromosomes", chromosomes}})},
      {"evaluation", evaluation},
      {"seconds", Json::object({{"simulate", detail::r4(run->seconds_simulate)}, {"index", detail::r4(run->seconds_index)}, {"map", detail::r4(run->seconds_map)},
                                {"detect", detail::r4(run->seconds_detect)}, {"compare", detail::r4(run->seconds_compare)}, {"library", detail::r4(d.seconds.library)},
                                {"evidence", detail::r4(d.seconds.evidence)}, {"cluster", detail::r4(d.seconds.cluster)}, {"call", detail::r4(d.seconds.call)},
                                {"depth", detail::r4(d.seconds.depth)}, {"segment", detail::r4(d.seconds.segment)}, {"snv", detail::r4(d.seconds.snv)}, {"allelic", detail::r4(d.seconds.allelic)}})},
      {"threads", worker.pool.size()},
  });
}

// The mates of a stretch of a structural run (the same request as structural.sample, with
// "window": {"contig", "begin", "end"} in the chromosome's numbering, at most 50 kb): each mate's
// interval, strand, mapping quality, clips and pair class, and its mate's place.
BENCH_OP("structural.window") {
  const auto run = detail::structural_run(worker, request);
  const Json& w = request["window"];
  const auto c = static_cast<std::uint32_t>(bounded(w, "contig", 0, 0, static_cast<long long>(run->contigs.size()) - 1));
  const auto off = run->contigs[c].offset;
  const long long b = integer(w, "begin", static_cast<long long>(off)) - static_cast<long long>(off), e = integer(w, "end", b + 5000 + static_cast<long long>(off)) - static_cast<long long>(off);
  if (b < 0 || e <= b || e - b > 50000) throw Fault("a window of at most 50 kb inside the contig's region");
  Json mates = Json::array();
  const auto& d = run->detection;
  for (std::size_t i = 0; i < d.pairs.size() && mates.size() < 4000; i++)
    for (int m = 0; m < 2; m++) {
      const auto* x = d.pairs[i].mate[m];
      if (!x || x->contig != c || static_cast<long long>(x->target_end) <= b || static_cast<long long>(x->target_begin) >= e) continue;
      const auto* y = d.pairs[i].mate[1 - m];
      mates.push(Json::array({i, m, off + x->target_begin, off + x->target_end, x->reverse, x->mapq, x->clip_left(), x->clip_right(), gsv::pair_class_name(d.pairs[i].type),
                              y ? Json(y->contig) : Json(), y ? Json(run->contigs[y->contig].offset + y->target_begin) : Json(), run->mates[2 * i + m].size() > 1}));
    }
  return Json::object({{"simulated", true}, {"contig", c}, {"begin", off + static_cast<std::size_t>(b)}, {"end", off + static_cast<std::size_t>(e)},
                       {"columns", Json::array({"pair", "mate", "begin", "end", "reverse", "mapq", "clip_left", "clip_right", "class", "mate_contig", "mate_begin", "split"})}, {"mates", mates}});
}

}  // namespace workbench
