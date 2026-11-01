// The operations of expression predicted from sequence (ADR 0064): codon usage indices of a real
// genome's genes with their places on it, the agreement of the indices with measured protein
// abundance; σ⁷⁰ promoters learnt from E. coli's start sites and found along a region; and a model
// learnt from measured promoter activity, with its held-out predictions, learning curve and filters.
//
// Every computation is the domain's (genomics::codon, genomics::regulation) and the foundations'
// (tyche::learning); the inputs are real data from the data library, whose files the request names:
//   codon.indices      {"genome": {"file": ".../NC_000913.3.gb"}, "reference": "ribosomal" | "all", "host": {"file": ...}?,
//                       "abundance": ".../511145-WHOLE_ORGANISM-integrated.txt"?, "s_permutations": 199}
//   regulation.sigma70 {"tss": ".../ettwiller-2016-tss.gtf", "genome_fasta": ".../U00096.2.fasta", "sites": 2000}
//   regulation.promoters {"tss": …, "genome_fasta": …, "genome": {"file": …}, "begin": 0, "end": 20000, "threshold": 8}
//   regulation.learn   {"data": ".../filtered_test_data_with_MAUDE_expression.txt", "motifs": ".../JASPAR2024_….txt"?,
//                       "train": 8500, "epochs": 8, "filters": 32, "width": 15, "seed": 7}
//   regulation.predict {"model": <the key regulation.learn returned>, "genome": {"file": …}, "begin": …, "end": …}
// The σ⁷⁰ model and the learnt network are kept by the worker (keyed by their requests), so that
// scanning and predicting do not learn again.

#pragma once

#include "genomes.hpp"

#include <genomics/codon/abundance.hpp>
#include <genomics/codon/adaptation.hpp>
#include <genomics/codon/trna.hpp>
#include <genomics/regulation/learn.hpp>
#include <genomics/regulation/sigma70.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace workbench {

namespace reg {

namespace gc = genomics::codon;
namespace gr = genomics::regulation;

inline double reg_round(double x, double scale = 1e4) { return std::isfinite(x) ? std::round(x * scale) / scale : 0; }
inline double since(std::chrono::steady_clock::time_point t) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count(); }

// The first GenBank record of a request's genome, read again for its features (the worker's Genome
// keeps contigs and features, not the record's header).
inline genomics::sequence::GenBank record_of(const Json& spec) {
  if (!spec.is_object() || !spec.has("file")) throw Fault("codon usage needs a genome given as a GenBank file");
  auto recs = genomics::sequence::read_genbank(detail::read_text(spec["file"].text()));
  if (recs.empty()) throw Fault("no GenBank record in " + spec["file"].text());
  return std::move(recs.front());
}

// Models learnt once per worker, by the request that made them.
struct Store {
  std::mutex lock;
  std::map<std::string, std::shared_ptr<gr::Sigma70>> sigma;
  std::map<std::string, std::shared_ptr<gr::SequenceCnn>> cnn;
};
inline Store& store() {
  static Store s;
  return s;
}

inline std::shared_ptr<gr::Sigma70> sigma70(const Json& request, Worker& worker, Json* about = nullptr) {
  const std::string tss = str(request, "tss"), fasta = str(request, "genome_fasta");
  if (tss.empty() || fasta.empty()) throw Fault("σ70 needs the start sites (tss) and their genome (genome_fasta) of data/cache/ecoli-tss-cappable");
  const auto n = static_cast<std::size_t>(bounded(request, "sites", 2000, 100, 8000));
  const std::string key = tss + "|" + fasta + "|" + std::to_string(n);
  {
    std::lock_guard<std::mutex> g(store().lock);
    if (auto it = store().sigma.find(key); it != store().sigma.end()) return it->second;
  }
  auto sites = gr::read_start_sites(detail::read_text(tss));
  const auto records = genomics::sequence::read_fasta(detail::read_text(fasta));
  if (records.empty()) throw Fault("no sequence in " + fasta);
  const auto genome = genomics::sequence::encode_dna(records[0].sequence);
  std::sort(sites.begin(), sites.end(), [](const gr::StartSite& a, const gr::StartSite& b) { return a.strength > b.strength; });
  std::vector<genomics::sequence::Codes> windows;
  for (std::size_t i = 0; i < sites.size() && windows.size() < n; i++) {
    auto w = gr::site_window(genome, sites[i]);
    if (!w.empty()) windows.push_back(std::move(w));
  }
  const auto t = std::chrono::steady_clock::now();
  auto m = std::make_shared<gr::Sigma70>(gr::learn_sigma70(windows, {}, &worker.pool));
  if (about) *about = Json::object({{"windows", windows.size()}, {"seconds", since(t)}});
  std::lock_guard<std::mutex> g(store().lock);
  store().sigma[key] = m;
  return m;
}

inline Json logo_of(const genomics::annotation::Pwm& p) {
  Json cols = Json::array();
  for (const auto& c : genomics::annotation::logo(p)) cols.push(Json::object({{"bits", c.information}, {"heights", numbers(c.height, 3)}}));
  return cols;
}

inline Json pwm_json(const genomics::annotation::Pwm& p) {
  return Json::object({{"consensus", gr::consensus(p)}, {"probability", numbers(p.probability, 4)}, {"logo", logo_of(p)}});
}

}  // namespace reg

// Codon usage of a genome's genes: per gene its place, CAI, tAI, Fop, N_c (Wright, Sun), GC3s, MILC,
// Karlin and Mrázek's E and PHX; the reference set and its weights; the S value; and, given PaxDb's
// abundances, the rank correlation of each index with them.
BENCH_OP("codon.indices") {
  using namespace reg;
  const auto t0 = std::chrono::steady_clock::now();
  const auto record = record_of(request["genome"]);
  const auto set = gc::coding_genes(record, {static_cast<int>(integer(request, "table", 11)), 0});
  if (set.genes.empty()) throw Fault("the genome has no coding gene codon usage can count");
  const int table = set.genes.front().table;
  const auto& code = genomics::sequence::genetic_code(table);
  std::vector<gc::CodonCounts> counts;
  for (const auto& g : set.genes) counts.push_back(gc::codon_counts(g));
  // The reference set: the genome's ribosomal proteins, all its genes, or a host's ribosomal proteins.
  std::string reference = str(request, "reference", "ribosomal");
  gc::CodonCounts ref{};
  std::size_t ref_genes = 0;
  std::string ref_source = record.version;
  if (request["host"].is_object()) {
    const auto host = gc::coding_genes(record_of(request["host"]), {11, 0});
    const auto rp = gc::genes_by_product(host.genes, gc::RIBOSOMAL_PROTEINS);
    ref = gc::pooled_counts(host.genes, rp);
    ref_genes = rp.size();
    ref_source = host.accession;
    reference = "host ribosomal proteins";
  } else if (reference == "ribosomal") {
    const auto rp = gc::genes_by_product(set.genes, gc::RIBOSOMAL_PROTEINS);
    if (rp.empty()) { reference = "all"; }
    else { ref = gc::pooled_counts(set.genes, rp); ref_genes = rp.size(); }
  }
  if (reference == "all") { ref = gc::pooled_counts(set.genes); ref_genes = set.genes.size(); }
  const auto w = gc::relative_adaptiveness(ref, code);
  const auto optimal = gc::optimal_by_adaptiveness(w);
  const auto all = gc::pooled_counts(set.genes);
  // tAI: only where the wobble rules hold (not in mitochondria) and the genome has tRNA genes.
  const auto trnas = gc::trna_genes(record);
  const auto pool = gc::trna_pool(trnas);
  const bool mitochondrial = table == 2 || table == 3 || table == 4 || table == 5 || table == 9 || table == 13 || table == 14 || table == 21 || table == 24 || table == 33;
  std::unique_ptr<gc::TrnaAdaptiveness> ta;
  std::string tai_note;
  if (mitochondrial) tai_note = "tAI is not computed: mitochondria decode with superwobbling U34, which dos Reis et al.'s rules do not describe";
  else if (pool.empty()) tai_note = "tAI is not computed: the genome has no tRNA genes";
  else ta = std::make_unique<gc::TrnaAdaptiveness>(gc::trna_adaptiveness(pool, code));
  // Karlin and Mrázek's standards, from E. coli's gene names (where present).
  gc::ExpressionStandards st;
  st.all = all;
  st.ribosomal = ref;
  st.chaperones = gc::pooled_counts(set.genes, gc::genes_by_name(set.genes, gc::ecoli_chaperones()));
  st.factors = gc::pooled_counts(set.genes, gc::genes_by_name(set.genes, gc::ecoli_translation_factors()));
  const bool karlin = !gc::genes_by_name(set.genes, gc::ecoli_chaperones()).empty();
  Json genes = Json::array();
  std::vector<double> v_cai, v_tai, v_fop, v_nc, v_milc, v_e;
  for (std::size_t i = 0; i < set.genes.size(); i++) {
    const auto& g = set.genes[i];
    const auto& n = counts[i];
    const double c = gc::cai(n, w), t = ta ? gc::tai(n, *ta) : 0, f = gc::fop(n, optimal, code), nc = gc::enc_wright(n, code);
    const double milc = gc::milc(n, ref, code);
    const auto p = karlin ? gc::predicted_expression(n, st, code) : gc::PredictedExpression{};
    v_cai.push_back(c); v_tai.push_back(t); v_fop.push_back(f); v_nc.push_back(nc); v_milc.push_back(-milc); v_e.push_back(p.e);
    genes.push(Json::array({g.locus_tag, g.gene, g.product, g.start, g.stop, g.reverse, g.codons(), reg_round(c), reg_round(t), reg_round(f), reg_round(nc),
                            reg_round(gc::enc_sun(n, code)), reg_round(gc::third_positions(n, code).gc3s), reg_round(milc), reg_round(p.e), karlin && p.phx}));
  }
  Json weights = Json::array();
  for (int k = 0; k < 64; k++)
    if (code.amino_acids[static_cast<std::size_t>(k)] != '*') {
      const auto b = genomics::sequence::codon_of(k);
      const std::string codon = {genomics::sequence::dna_letter(b[0]), genomics::sequence::dna_letter(b[1]), genomics::sequence::dna_letter(b[2])};
      weights.push(Json::array({codon, std::string(1, code.amino_acids[static_cast<std::size_t>(k)]), reg_round(w.w[static_cast<std::size_t>(k)]),
                                ta ? Json(reg_round(ta->w[static_cast<std::size_t>(k)])) : Json(), all[static_cast<std::size_t>(k)], ref[static_cast<std::size_t>(k)]}));
    }
  Json result = Json::object({{"accession", set.accession}, {"length", set.length}, {"circular", set.circular}, {"table", table},
                              {"columns", Json::array({"locus_tag", "gene", "product", "start", "stop", "reverse", "codons", "cai", "tai", "fop", "nc_wright", "nc_sun", "gc3s", "milc", "e", "phx"})},
                              {"genes", genes}, {"skipped", set.skipped.total()},
                              {"reference", Json::object({{"kind", reference}, {"genes", ref_genes}, {"source", ref_source}})},
                              {"weights", weights}, {"trna", Json::object({{"genes", trnas.genes.size()}, {"isoacceptors", pool.size()}, {"note", tai_note}})}});
  if (ta) {
    const auto sv = gc::s_value(counts, *ta, code, {}, static_cast<std::size_t>(bounded(request, "s_permutations", 199, 0, 9999)), 1, &worker.pool);
    result.set("s_value", Json::object({{"s", reg_round(sv.s)}, {"p", sv.p_value}, {"permutations", sv.permutations},
                                        {"psi", numbers(sv.psi, 3)}, {"tai", numbers(sv.tai, 4)}, {"gc3s", numbers(sv.gc3s, 4)}, {"nc", numbers(sv.nc, 3)}}));
  }
  const std::string ab = str(request, "abundance");
  if (!ab.empty()) {
    const auto table_ab = gc::read_paxdb(detail::read_text(ab));
    Json agree = Json::object();
    for (const auto& [name, v] : std::vector<std::pair<std::string, const std::vector<double>*>>{{"cai", &v_cai}, {"tai", &v_tai}, {"fop", &v_fop}, {"nc_wright", &v_nc}, {"milc_negated", &v_milc}, {"e", &v_e}}) {
      const auto a = gc::agreement(set.genes, *v, table_ab);
      agree.set(name, Json::object({{"spearman", reg_round(a.spearman)}, {"pearson_log", reg_round(a.pearson_log)}, {"matched", a.matched}}));
    }
    const auto a = gc::agreement(set.genes, v_cai, table_ab);
    Json pairs = Json::array();
    for (std::size_t k = 0; k < a.gene.size(); k++) pairs.push(Json::array({a.gene[k], reg_round(v_cai[a.gene[k]]), reg_round(v_tai[a.gene[k]]), a.abundance[k]}));
    result.set("abundance", Json::object({{"source", table_ab.name}, {"agreement", agree}, {"points", pairs}, {"columns", Json::array({"gene", "cai", "tai", "ppm"})}}));
  }
  result.set("seconds", since(t0));
  return result;
}

// The σ⁷⁰ model learnt from E. coli's strongest start sites: its boxes as logos, its spacer and
// offset laws, the EM's trace.
BENCH_OP("regulation.sigma70") {
  using namespace reg;
  Json about = Json::object();
  const auto m = sigma70(request, worker, &about);
  return Json::object({{"minus10", pwm_json(m->minus10)}, {"minus35", pwm_json(m->minus35)}, {"spacers", Json::object({{"from", m->s_min}, {"p", numbers(m->pi_s, 4)}})},
                       {"offsets", Json::object({{"from", m->r_min}, {"p", numbers(m->pi_r, 4)}})}, {"gamma", reg_round(m->gamma)}, {"trace", numbers(m->trace, 3)},
                       {"iterations", m->iterations}, {"learnt", about}});
}

// σ⁷⁰ promoters along a region of a genome, both strands, above a threshold in bits.
BENCH_OP("regulation.promoters") {
  using namespace reg;
  const auto m = sigma70(request, worker);
  const auto g = worker.genome(request);
  const auto& contig = g->contig(integer(request, "contig", 0));
  const auto begin = static_cast<std::size_t>(std::max<long long>(0, integer(request, "begin", 0)));
  const auto end = std::min(contig.codes.size(), static_cast<std::size_t>(integer(request, "end", static_cast<long long>(begin + 20000))));
  if (end <= begin || end - begin > 2'000'000) throw Fault("a region of at most 2,000,000 bases");
  const double threshold = num(request, "threshold", 8.0);
  const genomics::sequence::Codes x(contig.codes.begin() + static_cast<std::ptrdiff_t>(begin), contig.codes.begin() + static_cast<std::ptrdiff_t>(end));
  const auto hits = gr::scan_promoters(*m, x, threshold);
  Json out = Json::array();
  for (const auto& h : hits)
    out.push(Json::array({begin + h.minus35, begin + h.minus10, h.spacer, h.reverse, reg_round(h.score), reg_round(h.score10), reg_round(h.score35), begin + h.start}));
  return Json::object({{"contig", contig.name}, {"begin", begin}, {"end", end}, {"threshold", threshold},
                       {"columns", Json::array({"minus35", "minus10", "spacer", "reverse", "bits", "bits10", "bits35", "start"})}, {"hits", out}});
}

// A sequence-to-expression model learnt from measured promoter activity (a stated subset, families
// held out): the k-mer baseline and the CNN, held-out predictions, learning curve, filters as logos
// matched to known motifs.
BENCH_OP("regulation.learn") {
  using namespace reg;
  const std::string data = str(request, "data");
  if (data.empty()) throw Fault("regulation.learn needs measured activities (data: data/cache/dream-promoters/…)");
  const auto t0 = std::chrono::steady_clock::now();
  const auto all = gr::read_measured(detail::read_text(data));
  if (all.size() < 1000) throw Fault("fewer than 1000 measured sequences in " + data);
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 7));
  const auto parts = gr::grouped_split(all, {{17, 20}, {77, 20}}, 0.02, 0.05, seed);
  const auto n_train = static_cast<std::size_t>(bounded(request, "train", 8500, 500, 60000));
  std::vector<std::size_t> rows(parts.train.begin(), parts.train.begin() + static_cast<std::ptrdiff_t>(std::min(parts.train.size(), n_train)));
  const auto tr = gr::subset(all, rows), va = gr::subset(all, parts.validation), te = gr::subset(all, parts.test);
  gr::KmerBaselineOptions bo;
  bo.from = 17;
  bo.length = 80;
  bo.gapped_k = 0;
  bo.folds = 3;
  bo.lambda_count = 12;
  const auto t1 = std::chrono::steady_clock::now();
  const auto base = gr::fit_kmer_baseline(tr, bo);
  const double base_seconds = since(t1);
  const auto base_pred = gr::predict(base, te.sequences);
  gr::CnnOptions co;
  co.epochs = static_cast<std::size_t>(bounded(request, "epochs", 8, 1, 60));
  co.filters = static_cast<std::size_t>(bounded(request, "filters", 32, 4, 256));
  co.width = static_cast<std::size_t>(bounded(request, "width", 15, 4, 31));
  co.seed = seed;
  const auto t2 = std::chrono::steady_clock::now();
  auto cnn = std::make_shared<gr::SequenceCnn>(gr::train_cnn(tr, va, co, &worker.pool));
  const double cnn_seconds = since(t2);
  const auto cnn_pred = gr::predict(*cnn, te.sequences, &worker.pool);
  const auto ab = gr::agreement(base_pred, te.values), ac = gr::agreement(cnn_pred, te.values);
  Json curve = Json::array();
  for (const auto& e : cnn->training.curve) curve.push(Json::array({e.epoch, reg_round(e.train_loss), reg_round(e.validation_loss), reg_round(e.seconds)}));
  // Filters as motifs, matched to known ones when a motif file is given.
  std::vector<gr::KnownMotif> known;
  if (!str(request, "motifs").empty()) known = gr::read_meme_motifs(detail::read_text(str(request, "motifs")));
  Json filters = Json::array();
  for (const auto& f : gr::filter_motifs(*cnn, te.sequences)) {
    Json match;
    const auto [a, b] = gr::informative_core(f.pwm);
    if (!known.empty() && f.sites >= 20 && b - a >= 5) {
      std::vector<double> core(f.pwm.probability.begin() + static_cast<std::ptrdiff_t>(4 * a), f.pwm.probability.begin() + static_cast<std::ptrdiff_t>(4 * b));
      const auto mm = gr::match_motif(core, b - a, known, 5, 100, f.filter + 1);
      match = Json::object({{"id", known[mm.motif].id}, {"name", known[mm.motif].name}, {"similarity", reg_round(mm.similarity)}, {"p", mm.p_value}, {"reverse", mm.reverse}, {"offset", mm.offset}});
    }
    filters.push(Json::object({{"filter", f.filter}, {"sites", f.sites}, {"bits", reg_round(f.information)}, {"consensus", gr::consensus(f.pwm)}, {"logo", logo_of(f.pwm)}, {"core", Json::array({a, b})}, {"match", match}}));
  }
  // The baseline's strongest words.
  std::vector<std::size_t> order(base.fit.beta.size());
  for (std::size_t k = 0; k < order.size(); k++) order[k] = k;
  std::sort(order.begin(), order.end(), [&](std::size_t p, std::size_t q) { return std::fabs(base.fit.beta[p]) > std::fabs(base.fit.beta[q]); });
  Json words = Json::array();
  for (std::size_t k = 0; k < std::min<std::size_t>(30, order.size()) && base.fit.beta[order[k]] != 0; k++) words.push(Json::array({base.names[order[k]], reg_round(base.fit.beta[order[k]])}));
  Json held = Json::array();
  for (std::size_t i = 0; i < te.size(); i++) held.push(Json::array({reg_round(te.values[i]), reg_round(base_pred[i]), reg_round(cnn_pred[i])}));
  const std::string key = data + "|" + std::to_string(n_train) + "|" + std::to_string(co.epochs) + "|" + std::to_string(co.filters) + "|" + std::to_string(co.width) + "|" + std::to_string(seed);
  {
    std::lock_guard<std::mutex> g(store().lock);
    store().cnn[key] = cnn;
  }
  return Json::object({{"model", key}, {"sizes", Json::object({{"train", tr.size()}, {"validation", va.size()}, {"test", te.size()}, {"all", all.size()}})},
                       {"baseline", Json::object({{"pearson", reg_round(ab.pearson)}, {"spearman", reg_round(ab.spearman)}, {"r2", reg_round(ab.r2)}, {"nonzero", base.nonzero},
                                                  {"columns", base.names.size()}, {"lambda", base.lambda}, {"seconds", base_seconds}, {"words", words}})},
                       {"cnn", Json::object({{"pearson", reg_round(ac.pearson)}, {"spearman", reg_round(ac.spearman)}, {"r2", reg_round(ac.r2)}, {"architecture", cnn->net.describe()},
                                             {"parameters", cnn->net.size()}, {"best_epoch", cnn->training.best_epoch}, {"seconds", cnn_seconds}, {"curve", curve}})},
                       {"held_out", held}, {"held_out_columns", Json::array({"measured", "baseline", "cnn"})}, {"filters", filters}, {"seconds", since(t0)}});
}

// A learnt network's prediction along a region of a real genome: each window of the model's length,
// every `step` bases, on the forward strand. (The model is a yeast promoter model: on another genome
// the numbers are a model's reading, not a measurement.)
BENCH_OP("regulation.predict") {
  using namespace reg;
  std::shared_ptr<gr::SequenceCnn> m;
  {
    std::lock_guard<std::mutex> g(store().lock);
    if (auto it = store().cnn.find(str(request, "model")); it != store().cnn.end()) m = it->second;
  }
  if (!m) throw Fault("no learnt model by that key on this worker: run regulation.learn first");
  const auto g = worker.genome(request);
  const auto& contig = g->contig(integer(request, "contig", 0));
  const auto begin = static_cast<std::size_t>(std::max<long long>(0, integer(request, "begin", 0)));
  const auto end = std::min(contig.codes.size(), static_cast<std::size_t>(integer(request, "end", static_cast<long long>(begin + 5000))));
  const auto step = static_cast<std::size_t>(bounded(request, "step", 10, 1, 1000));
  if (end < begin + m->length || (end - begin) / step > 200000) throw Fault("a region longer than the model's window, at most 200,000 windows");
  std::vector<genomics::sequence::Codes> windows;
  std::vector<std::size_t> at;
  for (std::size_t p = begin; p + m->length <= end; p += step) {
    windows.emplace_back(contig.codes.begin() + static_cast<std::ptrdiff_t>(p), contig.codes.begin() + static_cast<std::ptrdiff_t>(p + m->length));
    at.push_back(p);
  }
  const auto y = gr::predict(*m, windows, &worker.pool);
  return Json::object({{"contig", contig.name}, {"window", m->length}, {"step", step}, {"positions", Json(at)}, {"values", numbers(y, 4)}});
}

}  // namespace workbench
