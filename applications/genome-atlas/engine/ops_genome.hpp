// The operations of the Genome view: the worker's target, a genome's summary, its tracks over
// a window (composition in bins, the cumulative skew, features, ORFs, the letters themselves
// when the window is short), one feature in full, and the search of a pattern on both strands.

#pragma once

#include "genomes.hpp"
#include "shared.hpp"

#include <genomics/annotation/genes.hpp>
#include <genomics/annotation/repeats.hpp>
#include <genomics/hmm/cpg.hpp>
#include <genomics/sequence/composition.hpp>
#include <genomics/sequence/genetic_code.hpp>
#include <genomics/sequence/kmer.hpp>

#include <calx/cpu/target.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <optional>
#include <string>
#include <vector>

namespace atlas {

namespace detail {

inline std::string feature_name(const seq::Feature& f) {
  for (const char* q : {"gene", "locus_tag", "product", "note", "label", "rpt_family", "mobile_element_type"}) {
    const auto v = f.qualifier(q);
    if (!v.empty()) return v.size() > 60 ? v.substr(0, 57) + "..." : v;
  }
  return f.key;
}

// The features to draw: every feature but a gene that a CDS or RNA of the same locus repeats.
inline std::vector<std::size_t> drawable(const Contig& c) {
  std::set<std::string> covered;
  for (const auto& f : c.features)
    if (f.key == "CDS" || f.key.find("RNA") != std::string::npos) {
      const auto t = f.qualifier("locus_tag");
      covered.insert(t.empty() ? f.qualifier("gene") : t);
    }
  std::vector<std::size_t> out;
  for (std::size_t i = 0; i < c.features.size(); i++) {
    const auto& f = c.features[i];
    if (f.key == "gene") {
      const auto t = f.qualifier("locus_tag");
      if (covered.count(t.empty() ? f.qualifier("gene") : t)) continue;
    }
    out.push_back(i);
  }
  return out;
}

inline Json feature_json(const seq::Feature& f, std::size_t index) {
  Json parts = Json::array();
  for (const auto& p : f.location.parts) parts.push(Json::array({p.begin, p.end}));
  Json j = Json::object({{"i", index}, {"key", f.key}, {"name", feature_name(f)}, {"strand", f.location.reverse() ? "-" : "+"}, {"parts", parts}});
  for (const char* q : {"product", "locus_tag", "gene"}) {
    const auto v = f.qualifier(q);
    if (!v.empty()) j.set(q, v);
  }
  if (f.has("pseudo") || f.has("pseudogene")) j.set("pseudo", true);
  return j;
}

}  // namespace detail

ATLAS_OP("target") {
  const auto& t = calx::cpu::host();
  Json ops = Json::array();
  for (const auto& [name, op] : operations()) ops.push(name);
  return Json::object({{"architecture", t.architecture}, {"brand", t.brand}, {"threads", worker.pool.size()}, {"hardware_threads", t.threads},
                       {"built_for", calx::cpu::ISA}, {"lane_bytes", calx::cpu::LANE_BYTES}, {"best", t.best()}, {"widest", t.widest()},
                       {"avx2", t.avx2}, {"avx512", t.avx512bw}, {"describe", calx::cpu::describe(t)}, {"operations", ops}});
}

ATLAS_OP("genome.summary") {
  const auto g = worker.genome(request);
  Json contigs = Json::array();
  for (const auto& c : g->contigs) {
    std::map<std::string, int> kinds;
    for (const auto& f : c.features) kinds[f.key]++;
    Json k = Json::object();
    for (const auto& [key, n] : kinds) k.set(key, n);
    const auto counts = seq::base_counts(c.codes);
    contigs.push(Json::object({{"name", c.name}, {"accession", c.accession}, {"definition", c.definition}, {"length", c.codes.size()},
                               {"circular", c.circular}, {"gc", std::round(seq::gc_fraction(c.codes) * 1e5) / 1e5}, {"ambiguous", counts[4]},
                               {"features", k}}));
  }
  Json problems = Json::array();
  for (const auto& p : g->problems) problems.push(Json::object({{"line", p.line}, {"message", p.message}}));
  Json repeats = Json::array();
  for (const auto& r : g->repeats) repeats.push(Json::object({{"element", r.element}, {"begin", r.begin}, {"end", r.end}, {"reverse", r.reverse}, {"tandem", r.tandem}}));
  return Json::object({{"title", g->title}, {"organism", g->organism}, {"lineage", Json(g->lineage)}, {"source", g->source}, {"length", g->length()},
                       {"table", g->table}, {"table_name", std::string(seq::genetic_code(g->table).name)}, {"contigs", contigs},
                       {"repeats", repeats}, {"problems", problems}});
}

// Tracks over [begin, end) of one contig, in `bins` bins.
ATLAS_OP("genome.tracks") {
  const auto g = worker.genome(request);
  const Contig& c = g->contig(integer(request, "contig", 0));
  const auto n = static_cast<long long>(c.codes.size());
  const long long begin = std::clamp<long long>(integer(request, "begin", 0), 0, n);
  const long long end = std::clamp<long long>(integer(request, "end", n), begin, n);
  const long long bins = bounded(request, "bins", 1000, 1, 20000);
  const auto span = static_cast<double>(end - begin);
  const long long nb = std::max<long long>(1, std::min<long long>(bins, end - begin));
  const auto edge = [&](long long b) { return begin + static_cast<long long>(std::llround(span * static_cast<double>(b) / static_cast<double>(nb))); };

  std::vector<double> gc(static_cast<std::size_t>(nb)), skew(static_cast<std::size_t>(nb)), cpg(static_cast<std::size_t>(nb)), entropy(static_cast<std::size_t>(nb)),
      ambiguous(static_cast<std::size_t>(nb));
  for (long long b = 0; b < nb; b++) {
    const auto lo = static_cast<std::size_t>(edge(b)), hi = static_cast<std::size_t>(edge(b + 1));
    const auto k = seq::base_counts(c.codes, lo, hi);
    const double acgt = static_cast<double>(k[0] + k[1] + k[2] + k[3]);
    const auto i = static_cast<std::size_t>(b);
    gc[i] = acgt > 0 ? static_cast<double>(k[1] + k[2]) / acgt : NAN;
    skew[i] = k[1] + k[2] ? (static_cast<double>(k[2]) - static_cast<double>(k[1])) / static_cast<double>(k[1] + k[2]) : 0;
    cpg[i] = seq::cpg_observed_expected(c.codes, lo, hi);
    entropy[i] = hi - lo >= 32 ? seq::kmer_entropy(c.codes, 3, lo, hi) : NAN;
    ambiguous[i] = hi > lo ? static_cast<double>(k[4]) / static_cast<double>(hi - lo) : 0;
  }
  // The cumulative skew of the whole contig, read at the bins' edges.
  std::vector<double> cumulative(static_cast<std::size_t>(nb) + 1);
  {
    double s = 0;
    long long b = 0;
    for (long long i = 0; i < n && b <= nb; i++) {
      while (b <= nb && edge(b) == i) cumulative[static_cast<std::size_t>(b++)] = s;
      s += (c.codes[static_cast<std::size_t>(i)] == seq::G) - (c.codes[static_cast<std::size_t>(i)] == seq::C);
    }
    while (b <= nb) cumulative[static_cast<std::size_t>(b++)] = s;
  }
  const auto whole = seq::cumulative_gc_skew(c.codes, std::max<std::size_t>(1, c.codes.size() / 2000));

  Json out = Json::object({{"contig", integer(request, "contig", 0)}, {"begin", begin}, {"end", end}, {"bins", nb}, {"length", n}});
  out.set("gc", numbers(gc));
  out.set("skew", numbers(skew));
  out.set("cpg", numbers(cpg));
  out.set("entropy", numbers(entropy, 3));
  out.set("ambiguous", numbers(ambiguous));
  out.set("cumulative", numbers(cumulative, 0));
  out.set("origin", whole.minimum);
  out.set("terminus", whole.maximum);
  out.set("mean_gc", std::round(seq::gc_fraction(c.codes) * 1e5) / 1e5);

  // Features overlapping the window.
  Json features = Json::array();
  std::size_t shown = 0, hidden = 0;
  const std::size_t cap = static_cast<std::size_t>(bounded(request, "max_features", 12000, 0, 100000));
  for (std::size_t i : detail::drawable(c)) {
    const auto& f = c.features[i];
    bool overlaps = false;
    for (const auto& p : f.location.parts) overlaps = overlaps || (static_cast<long long>(p.end) > begin && static_cast<long long>(p.begin) < end);
    if (!overlaps) continue;
    if (shown >= cap) { hidden++; continue; }
    features.push(detail::feature_json(f, i));
    shown++;
  }
  out.set("features", features);
  out.set("features_hidden", hidden);

  // ORFs of at least `orf_min` codons overlapping the window, by the genome's code.
  if (flag(request, "orfs", true)) {
    const auto min_codons = static_cast<std::size_t>(bounded(request, "orf_min", 100, 10, 100000));
    Json orfs = Json::array();
    std::size_t count = 0;
    for (const auto& o : seq::open_reading_frames(c.codes, seq::genetic_code(g->table), min_codons)) {
      if (static_cast<long long>(o.end) <= begin || static_cast<long long>(o.begin) >= end) continue;
      if (++count > 8000) break;
      // The frame on the forward strand's coordinates, for drawing the six frames: 0-2 forward, 3-5 reverse.
      const long long frame = o.reverse ? 3 + static_cast<long long>((n - static_cast<long long>(o.end)) % 3) : static_cast<long long>(o.begin % 3);
      orfs.push(Json::array({o.begin, o.end, frame}));
    }
    out.set("orfs", orfs);
  }
  if (end - begin <= bounded(request, "letters_max", 4000, 0, 50000)) out.set("letters", seq::decode_dna(seq::Codes(c.codes.begin() + begin, c.codes.begin() + end)));
  return out;
}

// One feature in full: its qualifiers, its bases' composition, its protein and its codons.
ATLAS_OP("genome.feature") {
  const auto g = worker.genome(request);
  const Contig& c = g->contig(integer(request, "contig", 0));
  const long long i = integer(request, "index", -1);
  if (i < 0 || static_cast<std::size_t>(i) >= c.features.size()) throw Fault("no feature " + std::to_string(i));
  const auto& f = c.features[static_cast<std::size_t>(i)];
  Json q = Json::array();
  for (const auto& [k, v] : f.qualifiers) q.push(Json::array({k, v}));
  const seq::Codes bases = seq::extract(c.codes, f.location);
  Json out = detail::feature_json(f, static_cast<std::size_t>(i));
  out.set("location", f.location.text);
  out.set("qualifiers", q);
  out.set("length", bases.size());
  out.set("gc", std::round(seq::gc_fraction(bases) * 1e4) / 1e4);
  if (f.key == "CDS") {
    const std::string mine = seq::translate_cds(c.codes, f, g->table), theirs = f.qualifier("translation");
    out.set("translation", mine);
    out.set("agrees", theirs.empty() ? Json() : Json(mine == theirs));
    const auto counts = seq::codon_counts({bases});
    const auto& code = seq::genetic_code(f.qualifier("transl_table").empty() ? g->table : std::stoi(f.qualifier("transl_table")));
    Json codons = Json::array();
    for (int k = 0; k < 64; k++) {
      const auto cod = seq::codon_of(k);
      codons.push(Json::array({seq::decode_dna({cod[0], cod[1], cod[2]}), std::string(1, code.amino_acids[static_cast<std::size_t>(k)]), counts[static_cast<std::size_t>(k)]}));
    }
    out.set("codons", codons);
    out.set("nc", std::round(seq::effective_number_of_codons(counts, code) * 100) / 100);
    // GC at the third positions of codons.
    std::size_t gc3 = 0, n3 = 0;
    for (std::size_t k = 2; k < bases.size(); k += 3)
      if (bases[k] < 4) { n3++; gc3 += bases[k] == seq::C || bases[k] == seq::G; }
    out.set("gc3", n3 ? std::round(1e4 * static_cast<double>(gc3) / static_cast<double>(n3)) / 1e4 : 0);
  }
  if (bases.size() <= 20000) out.set("sequence", seq::decode_dna(bases));
  return out;
}

// A gene by its name (/gene, then /locus_tag): its bases and, for a CDS, its protein.
ATLAS_OP("genome.gene") {
  const auto g = worker.genome(request);
  const std::string name = str(request, "gene");
  const auto at = detail::find_gene(*g, name);
  if (!at) throw Fault("no gene " + name + " in " + g->title);
  const auto& c = g->contigs[at->contig];
  const auto& f = c.features[at->feature];
  Json out = detail::feature_json(f, at->feature);
  out.set("contig", at->contig);
  out.set("sequence", seq::decode_dna(seq::extract(c.codes, f.location)));
  if (f.key == "CDS") out.set("translation", seq::translate_cds(c.codes, f, g->table));
  out.set("organism", g->organism);
  return out;
}

// Annotation computed from the sequence alone, for one contig: genes predicted by the self-trained
// gene finder (compared with the record's CDSs: a CDS is found when a prediction ends at its stop
// on its strand), CpG islands by the HMM of Durbin et al., or maximal repeats on both strands by
// the suffix array. Each is computed once per worker and kept with the genome.
ATLAS_OP("genome.annotate") {
  const auto g = worker.genome(request);
  const long long ci = integer(request, "contig", 0);
  const Contig& c = g->contig(ci);
  const std::string what = str(request, "what", "genes");
  double seconds = 0;
  Json out = Json::object({{"what", what}, {"contig", ci}});
  if (what == "genes") {
    struct Found { std::vector<genomics::annotation::Gene> genes; genomics::annotation::GeneReport report; };
    const auto f = worker.derived<Found>(*g, "genes:" + std::to_string(ci), [&] {
      Found x;
      genomics::annotation::GeneOptions o;
      o.table = g->table;
      x.genes = genomics::annotation::find_genes(c.codes, o, &x.report);
      return x;
    }, &seconds);
    Json genes = Json::array();
    std::set<std::pair<std::size_t, bool>> stops;  // (stop end, strand) of the predictions
    for (const auto& x : f->genes) {
      genes.push(Json::array({x.begin, x.end, x.reverse, std::round(x.score * 10) / 10, x.start}));
      stops.insert({x.reverse ? x.begin : x.end, x.reverse});
    }
    std::size_t cds = 0, found = 0;
    std::set<std::pair<std::size_t, bool>> annotated;
    for (const auto& feat : c.features) {
      if (feat.key != "CDS" || feat.location.parts.size() != 1 || feat.has("pseudo")) continue;
      const auto& part = feat.location.parts.front();
      const std::pair<std::size_t, bool> stop{part.reverse ? part.begin : part.end, part.reverse};
      cds++;
      annotated.insert(stop);
      found += stops.count(stop);
    }
    std::size_t right = 0;
    for (const auto& st : stops) right += annotated.count(st);
    out.set("genes", genes);
    out.set("rounds", f->report.rounds);
    out.set("annotated", cds);
    out.set("sensitivity", cds ? Json(static_cast<double>(found) / static_cast<double>(cds)) : Json());
    out.set("specificity", stops.empty() || !cds ? Json() : Json(static_cast<double>(right) / static_cast<double>(stops.size())));
  } else if (what == "islands") {
    const auto f = worker.derived<std::vector<genomics::hmm::Island>>(*g, "islands:" + std::to_string(ci), [&] {
      genomics::hmm::IslandOptions o;
      o.min_length = 200;
      return genomics::hmm::cpg_islands(c.codes, o);
    }, &seconds);
    Json list = Json::array();
    for (const auto& x : *f) list.push(Json::array({x.begin, x.end, std::round(x.log_odds * 10) / 10, std::round(x.gc * 1000) / 1000, std::round(x.cpg_observed_expected * 1000) / 1000}));
    out.set("islands", list);
  } else if (what == "repeats") {
    const auto min_length = static_cast<std::size_t>(bounded(request, "min_length", 100, 20, 100000));
    const auto f = worker.derived<std::vector<genomics::annotation::Repeat>>(*g, "repeats:" + std::to_string(ci) + ":" + std::to_string(min_length), [&] {
      genomics::annotation::RepeatOptions o;
      o.min_length = min_length;
      o.both_strands = true;
      return genomics::annotation::maximal_repeats(c.codes, o);
    }, &seconds);
    Json list = Json::array();
    for (const auto& r : *f) {
      if (list.size() >= 3000) break;
      Json occ = Json::array();
      for (const auto& o : r.occurrences) occ.push(Json::array({o.position, o.reverse}));
      list.push(Json::object({{"length", r.length}, {"supermaximal", r.supermaximal}, {"at", occ}}));
    }
    out.set("repeats", list);
    out.set("count", f->size());
    out.set("min_length", min_length);
  } else {
    throw Fault("annotate genes, islands or repeats");
  }
  out.set("seconds", seconds);
  return out;
}

// The exact occurrences of a pattern on both strands of every contig.
ATLAS_OP("genome.search") {
  const auto g = worker.genome(request);
  const std::string text = str(request, "pattern");
  if (text.size() < 4 || text.size() > 1000) throw Fault("a pattern has 4 to 1000 bases");
  const seq::Codes p = seq::encode_dna(text);
  if (seq::count_ambiguous(p)) throw Fault("a pattern is written in A, C, G and T");
  const seq::Codes rc = seq::reverse_complement(p);
  const std::size_t cap = static_cast<std::size_t>(bounded(request, "max_hits", 2000, 1, 100000));
  Json hits = Json::array();
  std::size_t total = 0;
  for (std::size_t ci = 0; ci < g->contigs.size(); ci++) {
    const auto& codes = g->contigs[ci].codes;
    for (const auto* q : {&p, &rc}) {
      if (q == &rc && rc == p) break;  // a palindrome is found once
      auto it = codes.begin();
      while ((it = std::search(it, codes.end(), q->begin(), q->end())) != codes.end()) {
        total++;
        if (hits.size() < cap) hits.push(Json::array({ci, static_cast<std::size_t>(it - codes.begin()), q == &rc}));
        ++it;
      }
    }
  }
  return Json::object({{"pattern", text}, {"total", total}, {"hits", hits}, {"method", "scan"}});
}

}  // namespace atlas
