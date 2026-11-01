// The operations of the somatic genome: a tumour and its normal simulated on a stretch of a real
// reference, sequenced, mapped and analysed as a cancer genome is (genomics::somatic), with the
// simulator's truth beside every answer; and a cohort's catalogues for de novo signatures.
//
// - `somatic.tumour`: the reference stretch [begin, begin + length) of one contig of the genome
//   (a GenBank file of the data library, e.g. NC_000913.3, or a simulated genome); a normal diploid
//   mutated from it; a tumour of a trunk and nested subclones, each with its own synthetic
//   signature, a random allele-specific copy-number profile (with or without whole-genome
//   doubling), microsatellites that slip; reads of both; the mapper, the normal's caller, MuTect's
//   tests, ASPCF and ASCAT, cancer cell fractions and their clones, the trees the sum rule allows,
//   the catalogue refitted and bootstrapped, MSIsensor's test, the burden; and the evaluation of
//   each against the truth. Everything simulated is labelled so ("simulated": true); the rates are
//   raised far above a real tumour's so that a stretch of a bacterial genome shows clones and
//   signatures.
// - `somatic.signatures`: a cohort of catalogues drawn from the synthetic signatures with random
//   exposures (simulate_catalogues), the NMF stability of ranks 1 to max_rank, the signatures
//   extracted at the chosen rank and matched to the truth.
//
// Every computation is the domain's (genomics::simulate, genomics::mapping, genomics::variant,
// genomics::somatic) and the foundations' (tyche, calx); the engine gathers, times and reports.

#pragma once

#include "genomes.hpp"

#include <genomics/mapping/mapper.hpp>
#include <genomics/sequence/contexts.hpp>
#include <genomics/simulate/mutate.hpp>
#include <genomics/simulate/reads.hpp>
#include <genomics/simulate/tumour.hpp>
#include <genomics/somatic/catalogue.hpp>
#include <genomics/somatic/clonality.hpp>
#include <genomics/somatic/copy_number.hpp>
#include <genomics/somatic/msi.hpp>
#include <genomics/somatic/mutect.hpp>
#include <genomics/somatic/signatures.hpp>
#include <genomics/variant/caller.hpp>
#include <genomics/variant/evaluate.hpp>

#include <calx/random/stream.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace workbench {

namespace gsom = genomics::somatic;

namespace detail {

inline double since_s(std::chrono::steady_clock::time_point t) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count(); }
inline double r4(double x) { return std::isfinite(x) ? std::round(x * 1e4) / 1e4 : 0; }

// The synthetic signatures as JSON.
inline Json signatures_json(const std::vector<sim::NamedSignature>& s) {
  Json a = Json::array();
  for (const auto& x : s) a.push(Json::object({{"name", x.name}, {"profile", numbers(x.profile, 6)}}));
  return a;
}

inline Json class_names() {
  Json sbs = Json::array(), dbs = Json::array();
  for (std::size_t k = 0; k < seq::SBS96; k++) sbs.push(seq::sbs96_name(k));
  for (std::size_t k = 0; k < seq::DBS78; k++) dbs.push(seq::dbs78_name(k));
  return Json::object({{"sbs96", sbs}, {"dbs78", dbs}});
}

inline Json stability_json(const tyche::multivariate::RankStability& r) {
  return Json::object({{"rank", r.rank}, {"stability", r4(r.average_stability)}, {"minimum", r4(r.minimum_stability)}, {"relative_error", r4(r.relative_error)},
                       {"stable", r.stable}, {"cluster_stability", numbers(r.stability)}});
}

}  // namespace detail

BENCH_OP("somatic.tumour") {
  using clock = std::chrono::steady_clock;
  Json timing = Json::object();
  auto t = clock::now();
  const auto g = worker.genome(request);
  const auto& contig = g->contig(integer(request, "contig", 0));
  const auto begin = static_cast<std::size_t>(bounded(request, "begin", 0, 0, 1'000'000'000));
  const auto length = static_cast<std::size_t>(bounded(request, "length", 600000, 50000, 3'000'000));
  if (begin + length > contig.codes.size()) throw Fault("the stretch [begin, begin + length) lies beyond the contig (" + std::to_string(contig.codes.size()) + " bases)");
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
  seq::Codes ref(contig.codes.begin() + static_cast<std::ptrdiff_t>(begin), contig.codes.begin() + static_cast<std::ptrdiff_t>(begin + length));
  for (auto& c : ref)
    if (c > 3) c = seq::N;
  const std::vector<seq::Codes> contigs{ref};

  // The microsatellites of the stretch: homopolymers of 7 or more by default, below MSIsensor's 10,
  // because a bacterium has almost none of 10 (E. coli K-12: one in 4.6 Mb, 1239 of 7 or more).
  gsom::MicrosatelliteOptions mo;
  mo.min_homopolymer = static_cast<std::size_t>(bounded(request, "min_homopolymer", 7, 5, 50));
  mo.min_repeats = static_cast<std::size_t>(bounded(request, "min_repeats", 5, 3, 40));
  const auto microsatellites = gsom::find_microsatellites(contigs, mo);

  // The normal and the tumour.
  sim::MutationRates rates;
  rates.snv = num(request, "germline_snv", 1.5e-3);
  rates.indel = num(request, "germline_indel", 1e-4);
  if (!(rates.snv > 0 && rates.snv < 0.01 && rates.indel >= 0 && rates.indel < 0.002)) throw Fault("germline rates: SNVs in (0, 0.01), indels in [0, 0.002)");
  const sim::Diploid normal = sim::mutate_diploid(ref, rates, seed);
  const auto synthetic = sim::synthetic_signatures();
  sim::TumourSpec spec;
  spec.purity = num(request, "purity", 0.7);
  if (!(spec.purity >= 0.2 && spec.purity <= 1)) throw Fault("a purity between 0.2 and 1");
  const double per_mb = num(request, "mutations_per_mb", 2000);
  if (!(per_mb >= 10 && per_mb <= 20000)) throw Fault("mutations per megabase between 10 and 20000");
  const double mb = static_cast<double>(length) / 1e6;
  const auto n_of = [&](double share) { return static_cast<std::size_t>(std::llround(share * per_mb * mb)); };
  const double ccf1 = num(request, "subclone_ccf", 0.55), ccf2 = num(request, "subsubclone_ccf", 0.3);
  if (!(ccf1 > 0 && ccf1 < 1 && ccf2 > 0 && ccf2 < ccf1)) throw Fault("CCFs 0 < subsubclone < subclone < 1");
  spec.signatures = {synthetic[0].profile, synthetic[4].profile, synthetic[1].profile, synthetic[3].profile};
  spec.clones = {{-1, 1.0, n_of(0.5), n_of(0.025), {0.7, 0.3, 0, 0}}, {0, ccf1, n_of(0.3), n_of(0.015), {0.2, 0, 0.8, 0}}, {1, ccf2, n_of(0.2), n_of(0.01), {0, 0, 0, 1}}};
  sim::CopyNumberSpec cn;
  cn.segments = static_cast<std::size_t>(bounded(request, "segments", 7, 1, 30));
  cn.min_length = std::max<std::size_t>(20000, length / (4 * cn.segments));
  cn.whole_genome_doubling = flag(request, "whole_genome_doubling", false);
  spec.copy_number = sim::random_copy_number({length}, cn, seed + 1);
  spec.whole_genome_doubling = cn.whole_genome_doubling;
  for (const auto& m : microsatellites) spec.msi.sites.push_back({m.contig, m.begin, m.unit, m.copies});
  spec.msi.rate = num(request, "msi_rate", 0.3);
  if (!(spec.msi.rate >= 0 && spec.msi.rate <= 1)) throw Fault("an MSI rate between 0 and 1");
  const sim::Tumour tumour = sim::make_tumour(contigs, {normal}, spec, seed + 2);
  timing.set("simulate", detail::since_s(t));

  // Reads, mapped.
  t = clock::now();
  sim::ShortReadModel model;
  const double tdepth = num(request, "tumour_depth", 60), ndepth = num(request, "normal_depth", 35);
  if (!(tdepth >= 10 && tdepth <= 200 && ndepth >= 10 && ndepth <= 200)) throw Fault("depths between 10 and 200");
  const auto treads = sim::tumour_reads(tumour, sim::reads_for_depth(tumour, tdepth, model.length), model, seed + 3);
  const auto nreads = sim::short_reads(sim::normal_source(tumour), sim::reads_for_depth(tumour, ndepth, model.length), model, seed + 4);
  timing.set("reads", detail::since_s(t));
  t = clock::now();
  const genomics::mapping::Mapper mapper(worker.pool, contigs, genomics::mapping::MapperOptions::short_reads());
  const auto tmap = mapper.map_all(worker.pool, treads), nmap = mapper.map_all(worker.pool, nreads);
  timing.set("map", detail::since_s(t));
  const auto tviews = genomics::variant::read_views(treads), nviews = genomics::variant::read_views(nreads);
  const gsom::MappedReads T(tmap, tviews), N(nmap, nviews);

  // Somatic calls.
  t = clock::now();
  const auto calls = gsom::call_somatic(worker.pool, contigs, T, N);
  timing.set("call", detail::since_s(t));

  // Copy number.
  t = clock::now();
  const auto ncalls = genomics::variant::call_variants(worker.pool, contigs, nmap, nviews);
  const auto profile = gsom::copy_number(worker.pool, contigs, T, N, ncalls);
  timing.set("copy_number", detail::since_s(t));

  // Clonality.
  t = clock::now();
  auto ccf = profile.fit.found ? gsom::cancer_cell_fractions(calls.calls, profile.segments, profile.fit.purity) : std::vector<gsom::MutationCcf>{};
  const auto clusters = gsom::cluster_mutations(ccf, 6);
  std::vector<double> clone_ccf;
  for (const auto& c : clusters.clones) clone_ccf.push_back(c.ccf);
  const auto trees = gsom::clone_trees(clone_ccf);
  timing.set("clonality", detail::since_s(t));

  // Catalogue and signatures.
  t = clock::now();
  const auto cat = gsom::catalogue(calls.calls, contigs);
  const auto fit = gsom::refit(cat.sbs96, spec.signatures);
  const auto boot = gsom::bootstrap_exposures(cat.sbs96, spec.signatures, static_cast<std::size_t>(bounded(request, "bootstrap", 100, 0, 1000)), seed + 5);
  timing.set("signatures", detail::since_s(t));

  // MSI.
  t = clock::now();
  const auto msi = gsom::msi(worker.pool, microsatellites, T, N);
  timing.set("msi", detail::since_s(t));
  const auto burden = gsom::mutational_burden(calls);

  // ---- the truth and the evaluation --------------------------------------------------------------
  std::map<std::tuple<std::size_t, std::string, std::string>, const sim::SomaticVariant*> truth_at;
  for (const auto& v : tumour.truth) truth_at[{v.position, v.ref, v.alt}] = &v;
  std::vector<genomics::variant::Call> as;
  for (const auto& c : calls.calls) as.push_back(gsom::as_call(c));
  std::vector<genomics::variant::TruthVariant> tv;
  for (const auto& v : tumour.truth) tv.push_back({v.contig, v.position, v.ref, v.alt, 1});
  const auto cmp = genomics::variant::compare(as, tv, contigs);
  // Normalised keys of the truth for the per-call labels.
  std::map<std::tuple<std::size_t, std::string, std::string>, const sim::SomaticVariant*> truth_norm;
  for (const auto& v : tumour.truth) {
    const auto n = genomics::variant::normalize(v.position, v.ref, v.alt, ref);
    truth_norm[{n.position, n.ref, n.alts[0]}] = &v;
  }
  const auto counts = [](const genomics::variant::Counts& c) {
    return Json::object({{"tp", c.tp}, {"fp", c.fp}, {"fn", c.fn}, {"precision", detail::r4(c.precision())}, {"recall", detail::r4(c.recall())}, {"f1", detail::r4(c.f1())}});
  };
  genomics::variant::Counts indel;
  indel.tp = cmp.insertion.tp + cmp.deletion.tp + cmp.other.tp;
  indel.fp = cmp.insertion.fp + cmp.deletion.fp + cmp.other.fp;
  indel.fn = cmp.insertion.fn + cmp.deletion.fn + cmp.other.fn;
  // Recall of substitutions by expected VAF.
  std::set<std::tuple<std::size_t, std::string, std::string>> passed;
  for (const auto& c : calls.calls)
    if (c.passed()) passed.insert({c.position, c.ref, c.alt});
  const std::vector<double> edges{0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1.01};
  Json by_vaf = Json::array();
  for (std::size_t b = 0; b + 1 < edges.size(); b++) {
    std::size_t n = 0, f = 0;
    for (const auto& v : tumour.truth)
      if (v.snv() && v.vaf >= edges[b] && v.vaf < edges[b + 1]) { n++; f += passed.count({v.position, v.ref, v.alt}) > 0; }
    by_vaf.push(Json::object({{"low", edges[b]}, {"high", std::min(1.0, edges[b + 1])}, {"truth", n}, {"found", f}}));
  }
  // Copy number by base.
  std::size_t right = 0, covered = 0;
  for (const auto& s : profile.segments)
    for (std::size_t p = s.begin; p < s.end; p += 100) {
      const auto& x = tumour.segment_at(0, p);
      covered += 100;
      right += 100 * (s.major() == std::max(x.copies[0], x.copies[1]) && s.minor() == std::min(x.copies[0], x.copies[1]));
    }

  // ---- the answer --------------------------------------------------------------------------------
  Json truth_clones = Json::array();
  for (std::size_t j = 0; j < tumour.clones.size(); j++) {
    std::size_t snv = 0, ind = 0;
    for (const auto& v : tumour.truth)
      if (v.clone == static_cast<int>(j) && !v.microsatellite) (v.snv() ? snv : ind)++;
    truth_clones.push(Json::object({{"parent", tumour.clones[j].parent}, {"ccf", tumour.clones[j].ccf}, {"fraction", detail::r4(tumour.fraction[j])}, {"snvs", snv}, {"indels", ind},
                                    {"exposures", numbers(tumour.clones[j].exposures)}}));
  }
  Json truth_segments = Json::array();
  for (const auto& s : tumour.copy_number)
    truth_segments.push(Json::object({{"begin", s.begin}, {"end", s.end}, {"a", s.copies[0]}, {"b", s.copies[1]}, {"major", std::max(s.copies[0], s.copies[1])},
                                      {"minor", std::min(s.copies[0], s.copies[1])}}));
  Json used_signatures = Json::array();
  for (const int k : {0, 4, 1, 3}) used_signatures.push(Json::object({{"name", synthetic[static_cast<std::size_t>(k)].name}, {"profile", numbers(synthetic[static_cast<std::size_t>(k)].profile, 6)}}));

  // Bins, at most 2000 points for the page.
  Json bins_begin = Json::array(), bins_logr = Json::array(), bins_seg = Json::array(), bins_tumour = Json::array(), bins_normal = Json::array();
  const std::size_t stride = std::max<std::size_t>(1, profile.bins.size() / 2000);
  for (std::size_t i = 0; i < profile.bins.size(); i += stride) {
    const auto& b = profile.bins[i];
    bins_begin.push(b.begin);
    bins_logr.push(std::isfinite(b.logr) ? Json(detail::r4(b.logr)) : Json());
    bins_seg.push(b.segment < profile.segments.size() ? Json(detail::r4(profile.segments[b.segment].logr)) : Json());
    bins_tumour.push(detail::r4(b.tumour));
    bins_normal.push(detail::r4(b.normal));
  }
  Json sites_pos = Json::array(), sites_baf = Json::array(), sites_logr = Json::array(), sites_seg = Json::array();
  for (const auto& h : profile.sites) {
    sites_pos.push(h.position);
    sites_baf.push(detail::r4(h.baf));
    sites_logr.push(detail::r4(h.logr));
    sites_seg.push(h.segment);
  }
  Json segments = Json::array();
  for (const auto& s : profile.segments)
    segments.push(Json::object({{"begin", s.begin}, {"end", s.end}, {"sites", s.sites()}, {"logr", detail::r4(s.logr)}, {"baf", detail::r4(s.baf)}, {"n_a", s.n_a}, {"n_b", s.n_b},
                                {"major", s.major()}, {"minor", s.minor()}, {"n_a_raw", detail::r4(s.n_a_raw)}, {"n_b_raw", detail::r4(s.n_b_raw)}}));
  Json optima = Json::array();
  for (const auto& o : profile.fit.optima) optima.push(Json::object({{"psi", o.psi}, {"rho", o.rho}, {"ploidy", detail::r4(o.ploidy)}, {"goodness", detail::r4(o.goodness)}}));
  Json sunrise = Json::object({{"psi", numbers(profile.fit.psi_grid, 2)}, {"rho", numbers(profile.fit.rho_grid, 2)}, {"distance", numbers(profile.fit.distances, 3)}, {"optima", optima}});

  // Calls with their CCFs and their truth.
  std::map<std::size_t, const gsom::MutationCcf*> ccf_of;
  for (const auto& m : ccf) ccf_of[m.call] = &m;
  Json call_list = Json::array();
  for (std::size_t i = 0; i < calls.calls.size(); i++) {
    const auto& c = calls.calls[i];
    Json f = Json::array();
    for (const auto& x : c.filters) f.push(x);
    Json o = Json::object({{"position", c.position}, {"ref", c.ref}, {"alt", c.alt}, {"indel", c.indel}, {"vaf", detail::r4(c.vaf())}, {"tumour_depth", c.tumour_depth},
                           {"tumour_alt", c.tumour_alt}, {"normal_depth", c.normal_depth}, {"normal_alt", c.normal_alt}, {"lod_tumour", detail::r4(c.lod_tumour)},
                           {"lod_normal", detail::r4(c.lod_normal)}, {"filters", f}});
    if (auto it = ccf_of.find(i); it != ccf_of.end()) {
      const auto& m = *it->second;
      o.set("ccf", detail::r4(m.ccf));
      o.set("ccf_low", m.ccf_low);
      o.set("ccf_high", m.ccf_high);
      o.set("clonal", m.clonal);
      o.set("multiplicity", m.multiplicity);
      o.set("cluster", m.cluster);
    }
    const auto n = genomics::variant::normalize(c.position, c.ref, c.alt, ref);
    if (auto it = truth_norm.find({n.position, n.ref, n.alts[0]}); it != truth_norm.end()) {
      const auto& v = *it->second;
      o.set("truth", Json::object({{"clone", v.clone}, {"ccf", v.ccf}, {"multiplicity", v.multiplicity}, {"vaf", detail::r4(v.vaf)}, {"signature", v.signature},
                                   {"microsatellite", v.microsatellite}, {"early", v.early}}));
    } else {
      o.set("truth", Json());
    }
    call_list.push(std::move(o));
  }
  Json clones = Json::array();
  for (const auto& c : clusters.clones) clones.push(Json::object({{"ccf", detail::r4(c.ccf)}, {"weight", detail::r4(c.weight)}, {"mutations", c.mutations}}));
  Json tree_list = Json::array();
  for (std::size_t k = 0; k < trees.trees.size() && k < 50; k++) tree_list.push(Json::object({{"parent", Json(trees.trees[k].parent)}, {"slack", detail::r4(trees.trees[k].slack)}}));
  Json forced = Json::array();
  for (const auto& [a, b] : trees.forced) forced.push(Json::array({a, b}));
  Json boot_list = Json::array();
  for (const auto& e : boot.exposures) boot_list.push(Json::object({{"estimate", detail::r4(e.estimate)}, {"mean", detail::r4(e.mean)}, {"sd", detail::r4(e.sd)}, {"low", detail::r4(e.low)}, {"high", detail::r4(e.high)}}));
  // The truth of the exposures: the passing true substitutions by the signature that made them.
  std::vector<double> want(spec.signatures.size(), 0);
  for (const auto& c : calls.calls)
    if (c.passed() && !c.indel)
      if (auto it = truth_at.find({c.position, c.ref, c.alt}); it != truth_at.end() && it->second->signature >= 0) want[static_cast<std::size_t>(it->second->signature)]++;
  // MSI sites.
  std::set<std::size_t> slipped(tumour.slipped.begin(), tumour.slipped.end());
  Json msi_sites = Json::array();
  std::size_t mtp = 0, mfp = 0, mfn = 0;
  for (const auto& x : msi.sites) {
    const auto& m = microsatellites[x.site];
    const bool slip = slipped.count(x.site) > 0;
    if (x.tested) {
      if (x.somatic && slip) mtp++;
      else if (x.somatic) mfp++;
      else if (slip) mfn++;
    }
    Json nd = Json::object(), td = Json::object();
    for (const auto& [k, v] : x.normal) nd.set(std::to_string(k), v);
    for (const auto& [k, v] : x.tumour) td.set(std::to_string(k), v);
    msi_sites.push(Json::object({{"position", m.begin}, {"unit", m.unit}, {"copies", m.copies}, {"motif", seq::decode_dna(m.motif)}, {"tested", x.tested}, {"p", x.test.p},
                                 {"fdr", detail::r4(x.fdr)}, {"somatic", x.somatic}, {"slipped", slip}, {"normal", nd}, {"tumour", td}}));
  }

  return Json::object({
      {"simulated", true},
      {"reference", Json::object({{"title", g->title}, {"contig", contig.name}, {"accession", contig.accession}, {"begin", begin}, {"length", length}, {"source", g->source}})},
      {"truth", Json::object({{"purity", tumour.purity}, {"ploidy", detail::r4(tumour.ploidy())}, {"whole_genome_doubling", tumour.whole_genome_doubling}, {"clones", truth_clones},
                              {"segments", truth_segments}, {"variants", tumour.truth.size()}, {"msi_sites", tumour.msi_sites}, {"msi_slipped", tumour.slipped.size()},
                              {"signatures", used_signatures}})},
      {"copy_number", Json::object({{"bin", 1000}, {"bins", Json::object({{"begin", bins_begin}, {"logr", bins_logr}, {"segment_logr", bins_seg}, {"tumour_depth", bins_tumour}, {"normal_depth", bins_normal}})},
                                    {"sites", Json::object({{"position", sites_pos}, {"baf", sites_baf}, {"logr", sites_logr}, {"segment", sites_seg}})},
                                    {"segments", segments},
                                    {"fit", Json::object({{"found", profile.fit.found}, {"purity", profile.fit.purity}, {"psi", profile.fit.psi}, {"ploidy", detail::r4(profile.ploidy)},
                                                          {"goodness", detail::r4(profile.fit.goodness)}, {"pass", profile.fit.pass}, {"nonaberrant", profile.fit.nonaberrant}})},
                                    {"sunrise", sunrise}})},
      {"calls", call_list},
      {"clonality", Json::object({{"clones", clones}, {"bic", numbers(clusters.bic, 2)}, {"trees", tree_list}, {"trees_consistent", trees.trees.size()}, {"trees_examined", trees.examined},
                                  {"best_tree", trees.best}, {"tree_nodes_ccf", numbers(trees.ccf)}, {"tree_nodes_clone", Json(trees.clone)}, {"forced", forced}})},
      {"catalogue", Json::object({{"sbs96", numbers(cat.sbs96, 0)}, {"dbs78", numbers(cat.dbs78, 0)}, {"substitutions", cat.substitutions}, {"doublets", cat.doublets},
                                  {"classes", detail::class_names()}})},
      {"signatures", Json::object({{"refit", Json::object({{"exposures", numbers(fit.exposures, 2)}, {"cosine", detail::r4(fit.cosine)}, {"reconstruction", numbers(fit.reconstruction, 3)}})},
                                   {"bootstrap", boot_list}, {"bootstrap_replicates", boot.replicates}, {"truth_exposures", numbers(want, 0)}})},
      {"msi", Json::object({{"sites", msi_sites}, {"tested", msi.tested}, {"somatic", msi.somatic}, {"score", detail::r4(msi.score)}, {"high", msi.high},
                            {"min_homopolymer", mo.min_homopolymer}, {"min_repeats", mo.min_repeats}})},
      {"burden", Json::object({{"per_megabase", detail::r4(burden.per_megabase)}, {"substitutions", burden.substitutions}, {"indels", burden.indels}, {"callable", burden.callable}})},
      {"evaluation", Json::object({{"snv", counts(cmp.snv)}, {"indel", counts(indel)}, {"recall_by_vaf", by_vaf},
                                   {"copy_number_bases_right", covered ? detail::r4(static_cast<double>(right) / static_cast<double>(covered)) : 0.0},
                                   {"purity_error", detail::r4(profile.fit.purity - tumour.purity)}, {"ploidy_error", detail::r4(profile.ploidy - tumour.ploidy())},
                                   {"msi", Json::object({{"tp", mtp}, {"fp", mfp}, {"fn", mfn}})}})},
      {"reads", Json::object({{"tumour", treads.size()}, {"normal", nreads.size()}, {"tumour_depth", tdepth}, {"normal_depth", ndepth}, {"length", model.length}})},
      {"timing", timing}});
}

BENCH_OP("somatic.signatures") {
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
  const auto n = static_cast<std::size_t>(bounded(request, "tumours", 40, 5, 300));
  const auto max_rank = static_cast<std::size_t>(bounded(request, "max_rank", 6, 1, 10));
  const double mean = num(request, "mean_mutations", 400);
  if (!(mean >= 20 && mean <= 100000)) throw Fault("mean mutations per tumour and signature between 20 and 100000");
  const auto synthetic = sim::synthetic_signatures();
  std::vector<std::size_t> use;
  if (request["signatures"].is_array())
    for (const auto& x : request["signatures"].as_array()) use.push_back(static_cast<std::size_t>(x.as_number()));
  else
    use = {0, 1, 2, 3};
  std::vector<std::vector<double>> truth;
  for (auto k : use) {
    if (k >= synthetic.size()) throw Fault("synthetic signatures are numbered 0 to " + std::to_string(synthetic.size() - 1));
    truth.push_back(synthetic[k].profile);
  }
  calx::random::Stream r(seed, 0x434f484fu /* "COHO" */, 0, 0);
  std::vector<std::vector<double>> exposures;
  for (std::size_t i = 0; i < n; i++) {
    std::vector<double> e(truth.size());
    for (auto& x : e) x = r.bernoulli(0.7) ? r.gamma(2, mean / 2) : 0;
    exposures.push_back(e);
  }
  const auto cohort = sim::simulate_catalogues(truth, exposures, seed + 1);
  tyche::multivariate::StabilityOptions so;
  so.replicates = static_cast<std::size_t>(bounded(request, "replicates", 20, 2, 200));
  so.nmf.max_iterations = 4000;
  so.nmf.tolerance = 1e-9;
  so.seed = seed + 2;
  const auto t = std::chrono::steady_clock::now();
  const auto ex = gsom::extract_signatures(cohort, 1, max_rank, so);
  const double seconds = detail::since_s(t);
  const auto extracted = gsom::columns_of(ex.solution().signatures);
  Json ranks = Json::array();
  for (const auto& k : ex.ranks) ranks.push(detail::stability_json(k));
  Json sigs = Json::array();
  for (const auto& s : extracted) sigs.push(numbers(s, 6));
  Json match = Json::array();
  if (extracted.size() <= truth.size()) {
    const auto m = gsom::match_signatures(extracted, truth);
    for (std::size_t i = 0; i < extracted.size(); i++) match.push(Json::object({{"truth", use[m.reference[i]]}, {"name", synthetic[use[m.reference[i]]].name}, {"cosine", detail::r4(m.cosine[i])}}));
  }
  Json cats = Json::array(), exps = Json::array();
  for (std::size_t i = 0; i < n; i++) {
    cats.push(numbers(cohort[i], 0));
    exps.push(numbers(exposures[i], 1));
  }
  std::vector<sim::NamedSignature> used;
  for (auto k : use) used.push_back(synthetic[k]);
  return Json::object({{"simulated", true}, {"truth", detail::signatures_json(used)}, {"exposures", exps}, {"catalogues", cats}, {"ranks", ranks},
                       {"chosen_rank", ex.solution().rank}, {"stable", ex.stable}, {"extracted", sigs}, {"match", match},
                       {"exposures_extracted", numbers(ex.solution().exposures.a, 1)}, {"classes", detail::class_names()}, {"seconds", seconds}});
}

}  // namespace workbench
