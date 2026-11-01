// The operations of association: a simulated cohort (clearly simulated) with genotypes in linkage, a
// trait of known heritability and genes under cis-eQTLs, analysed as a study would be — the scan with
// principal components, genomic control, a locus fine-mapped, eQTLs across genes, colocalisation,
// expression predicted from genotypes and its association (PrediXcan, S-PrediXcan), heritability and
// polygenic scores — every result beside the truth that made it.
//
// Every computation is the evolution domain's (evolution::association) and tyche's; the engine sizes
// the study, times the steps, and thins what a page draws.

#pragma once

#include "protocol.hpp"
#include "genomes.hpp"

#include <evolution/association/coloc.hpp>
#include <evolution/association/cohort.hpp>
#include <evolution/association/eqtl.hpp>
#include <evolution/association/finemap.hpp>
#include <evolution/association/heritability.hpp>
#include <evolution/association/ld.hpp>
#include <evolution/association/predict.hpp>
#include <evolution/association/prs.hpp>
#include <evolution/association/scan.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

namespace workbench {

namespace assoc = evolution::association;

namespace detail_association {

inline double round_to(double x, double scale = 1e4) { return std::isfinite(x) ? std::round(x * scale) / scale : 0; }
inline double mlog10(double p) { return std::isfinite(p) ? -std::log10(std::max(p, 1e-300)) : 0; }

inline calx::linalg::Rectangular covariates(const assoc::Study& st, std::span<const std::size_t> who, const assoc::Components* pcs, std::size_t k,
                                            const std::vector<std::vector<double>>& extra = {}) {
  calx::linalg::Rectangular C(who.size(), 3 + (pcs ? k : 0) + extra.size());
  for (std::size_t i = 0; i < who.size(); i++) {
    C(i, 0) = 1;
    C(i, 1) = st.sex[who[i]];
    C(i, 2) = st.age[who[i]];
    for (std::size_t c = 0; pcs && c < k; c++) C(i, 3 + c) = pcs->vectors[i * pcs->components + c];
    for (std::size_t e = 0; e < extra.size(); e++) C(i, 3 + (pcs ? k : 0) + e) = extra[e][i];
  }
  return C;
}

inline std::vector<double> gather(std::span<const double> x, std::span<const std::size_t> who) {
  std::vector<double> out(who.size());
  for (std::size_t i = 0; i < who.size(); i++) out[i] = x[who[i]];
  return out;
}

inline bool is_causal(const assoc::Study& st, std::size_t site) { return std::binary_search(st.trait_causal.begin(), st.trait_causal.end(), site); }

}  // namespace detail_association

// A simulated study, analysed end to end.
BENCH_OP("association.study") {
  namespace da = detail_association;
  using clock = std::chrono::steady_clock;
  const auto since = [](clock::time_point t) { return std::chrono::duration<double>(clock::now() - t).count(); };
  assoc::StudySpec spec;
  spec.cohort.individuals = static_cast<std::size_t>(bounded(request, "individuals", 3000, 600, 8000));
  spec.cohort.regions = static_cast<std::size_t>(bounded(request, "regions", 4, 1, 8));
  spec.cohort.populations = static_cast<std::size_t>(bounded(request, "populations", 2, 1, 4));
  spec.heritability = num(request, "heritability", 0.4);
  if (!(spec.heritability > 0.05 && spec.heritability < 0.9)) throw Fault("a heritability between 0.05 and 0.9");
  spec.causal = static_cast<std::size_t>(bounded(request, "causal", 20, 1, 500));
  spec.mediated = std::min(spec.heritability / 2, num(request, "mediated", 0.05));
  spec.genes = static_cast<std::size_t>(bounded(request, "genes", 30, 4, 120));
  spec.reference = static_cast<std::size_t>(bounded(request, "reference", 500, 200, 2000));
  spec.population_effect = num(request, "confounding", 0.25);
  spec.prevalence = num(request, "prevalence", 0.2);
  if (!(spec.prevalence > 0.01 && spec.prevalence < 0.6)) throw Fault("a prevalence between 0.01 and 0.6");
  if (spec.reference + 300 > spec.cohort.individuals) throw Fault("the cohort must hold the expression panel and an association sample of at least 300");
  const auto permutations = static_cast<std::size_t>(bounded(request, "permutations", 300, 50, 5000));
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
  Json seconds = Json::object();

  auto t = clock::now();
  const auto st = assoc::simulate_study(spec, seed, &worker.pool);
  seconds.set("simulate", since(t));
  const auto& g = st.genotypes;
  Json out = Json::object({{"simulated", true}, {"seed", static_cast<double>(seed)}, {"individuals", st.individuals}, {"sample", st.sample.size()},
                           {"reference", st.reference.size()}, {"sites", g.sites}, {"populations", st.populations},
                           {"heritability_true", da::round_to(st.heritability)}, {"causal", st.causal.size()},
                           {"trait_causal", st.trait_causal.size()}});
  {
    Json regions = Json::array();
    for (std::size_t c = 0; c < g.chromosomes.size(); c++) {
      const auto [lo, hi] = g.chromosome_sites(static_cast<std::uint32_t>(c));
      regions.push(Json::object({{"name", g.chromosomes[c]}, {"length", static_cast<double>(spec.cohort.region_length)}, {"sites", hi - lo},
                                 {"founder_sites", st.founder_sites[c]}}));
    }
    out.set("regions", regions);
  }

  // ---- structure and the scan ----
  t = clock::now();
  std::vector<std::size_t> thin;
  for (std::size_t j = 0; j < g.sites; j += std::max<std::size_t>(1, g.sites / 600)) thin.push_back(j);
  const auto pcs = assoc::principal_components(g, st.sample, thin, 4);
  const auto C = da::covariates(st, st.sample, &pcs, 4);
  const auto C0 = da::covariates(st, st.sample, nullptr, 0);
  const auto y = da::gather(st.trait, st.sample);
  const auto scan = assoc::linear_scan(g, st.sample, y, C, &worker.pool);
  const auto scan0 = assoc::linear_scan(g, st.sample, y, C0, &worker.pool);
  const auto cases = da::gather(st.cases, st.sample);
  const auto lscan = assoc::logistic_scan(g, st.sample, cases, C, &worker.pool, 1e-5);
  seconds.set("scan", since(t));
  {
    Json points = Json::array();
    const std::size_t step = std::max<std::size_t>(1, st.sample.size() / 1500);
    for (std::size_t i = 0; i < st.sample.size(); i += step)
      points.push(Json::array({da::round_to(pcs.vectors[i * 4]), da::round_to(pcs.vectors[i * 4 + 1]), st.population[st.sample[i]]}));
    out.set("pcs", Json::object({{"values", numbers(pcs.values, 3)}, {"points", points}}));
    Json manhattan = Json::array();
    for (std::size_t j = 0; j < g.sites; j++)
      manhattan.push(Json::array({g.chromosome[j], static_cast<double>(g.positions[j]), da::round_to(da::mlog10(scan.p[j]), 1e3), da::is_causal(st, j) ? 1 : 0,
                                  da::round_to(da::mlog10(lscan.p[j]), 1e3)}));
    Json qq = Json::array();
    for (const auto& q : assoc::qq_points(scan.p)) qq.push(Json::array({da::round_to(q.expected, 1e3), da::round_to(q.observed, 1e3), da::round_to(q.lower, 1e3), da::round_to(q.upper, 1e3)}));
    Json qq0 = Json::array();
    for (const auto& q : assoc::qq_points(scan0.p, 1e-2, 200)) qq0.push(Json::array({da::round_to(q.expected, 1e3), da::round_to(q.observed, 1e3)}));
    std::size_t found = 0;
    for (const std::size_t c : st.causal) found += scan.p[c] < 5e-8;
    out.set("gwas", Json::object({{"lambda", da::round_to(assoc::genomic_control(scan.p))}, {"lambda_no_pcs", da::round_to(assoc::genomic_control(scan0.p))},
                                  {"lambda_logistic", da::round_to(assoc::genomic_control(lscan.p))}, {"manhattan", manhattan}, {"qq", qq}, {"qq_no_pcs", qq0},
                                  {"causal_found", found}, {"df", scan.df}}));
  }

  // ---- the strongest locus, fine-mapped ----
  t = clock::now();
  assoc::ClumpOptions co;
  co.p1 = 5e-8;
  co.r2 = 0.1;
  const auto clumps = assoc::clump(g, st.sample, scan.p, co);
  {
    Json hits = Json::array();
    for (const auto& c : clumps)
      hits.push(Json::object({{"site", c.index}, {"chromosome", g.chromosome[c.index]}, {"position", static_cast<double>(g.positions[c.index])},
                              {"mlog10p", da::round_to(da::mlog10(scan.p[c.index]), 1e3)}, {"beta", da::round_to(scan.beta[c.index])}, {"members", c.members.size()},
                              {"causal", da::is_causal(st, c.index)}}));
    out.set("hits", hits);
  }
  if (!clumps.empty()) {
    const std::size_t lead = clumps.front().index;
    const std::uint64_t pos = g.positions[lead], from = pos > 150'000 ? pos - 150'000 : 0, to = pos + 150'000;
    const auto sites = assoc::sites_between(g, g.chromosome[lead], from, to);
    const auto fit = assoc::finemap_susie(g, st.sample, sites, y, C, {}, &worker.pool);
    std::vector<double> b, se, z;
    for (const std::size_t j : sites) {
      b.push_back(scan.beta[j]);
      se.push_back(scan.se[j]);
      z.push_back(scan.statistic[j]);
    }
    double sdy = 0;
    {
      const double m = std::accumulate(y.begin(), y.end(), 0.0) / static_cast<double>(y.size());
      for (const double v : y) sdy += (v - m) * (v - m);
      sdy = std::sqrt(sdy / static_cast<double>(y.size() - 1));
    }
    const auto abf = assoc::finemap_abf(b, se, 0.15 * sdy);
    const auto R = assoc::ld_correlation(g, st.sample, sites, &worker.pool);
    const auto rss = assoc::finemap_rss(z, R, st.sample.size());
    const std::size_t lead_k = static_cast<std::size_t>(std::find(sites.begin(), sites.end(), lead) - sites.begin());
    Json rows = Json::array();
    for (std::size_t k = 0; k < sites.size(); k++)
      rows.push(Json::array({static_cast<double>(g.positions[sites[k]]), da::round_to(da::mlog10(scan.p[sites[k]]), 1e3), da::round_to(R(k, lead_k) * R(k, lead_k), 1e3),
                             da::round_to(fit.pip[k], 1e3), da::round_to(abf.posterior[k], 1e3), da::round_to(rss.pip[k], 1e3), da::is_causal(st, sites[k]) ? 1 : 0}));
    const auto sets_json = [&](const tyche::regression::SusieFit& f) {
      Json sets = Json::array();
      std::vector<std::vector<std::size_t>> cs;
      for (const auto& set : f.sets) {
        Json members = Json::array();
        bool covers = false;
        for (const std::size_t k : set.variables) {
          members.push(static_cast<double>(g.positions[sites[k]]));
          covers = covers || da::is_causal(st, sites[k]);
        }
        cs.push_back(set.variables);
        sets.push(Json::object({{"effect", set.effect}, {"positions", members}, {"purity", da::round_to(set.min_abs_corr)}, {"coverage", da::round_to(set.coverage)},
                                {"log_bf", da::round_to(set.log_bayes_factor)}, {"covers", covers}}));
      }
      const auto truth = assoc::against_truth(cs, sites, st.trait_causal);
      return Json::object({{"sets", sets}, {"covering", truth.covering}, {"causal_in_region", truth.causal}, {"causal_found", truth.found}});
    };
    Json abf_set = Json::array();
    for (const std::size_t k : abf.credible) abf_set.push(static_cast<double>(g.positions[sites[k]]));
    out.set("locus", Json::object({{"chromosome", g.chromosome[lead]}, {"from", static_cast<double>(from)}, {"to", static_cast<double>(to)},
                                   {"lead", static_cast<double>(pos)}, {"columns", Json::array({"position", "mlog10p", "r2_lead", "pip_susie", "pip_abf", "pip_rss", "causal"})},
                                   {"sites", rows}, {"susie", sets_json(fit)}, {"rss", sets_json(rss)}, {"abf_set", abf_set},
                                   {"iterations", fit.iterations}, {"residual_variance", da::round_to(fit.residual_variance)}}));
  }
  seconds.set("finemap", since(t));

  // ---- cis-eQTLs, colocalisation ----
  t = clock::now();
  std::vector<std::vector<double>> expr_ref(st.genes.size());
  for (std::size_t q = 0; q < st.genes.size(); q++) expr_ref[q] = da::gather(st.gene_expression(q), st.reference);
  const auto factors = assoc::expression_components(expr_ref, 3);
  const auto pcs_ref = assoc::principal_components(g, st.reference, thin, 3);
  const auto Cref = da::covariates(st, st.reference, &pcs_ref, 3, factors);
  std::vector<assoc::GeneLocus> loci;
  for (const auto& gene : st.genes) loci.push_back({gene.chromosome, gene.tss});
  assoc::CisOptions eo;
  eo.window = 200'000;
  eo.permutations = permutations;
  eo.seed = seed;
  const auto map = assoc::cis_eqtl(g, st.reference, expr_ref, loci, Cref, eo, &worker.pool);
  seconds.set("eqtl", since(t));
  {
    Json genes = Json::array(), colocs = Json::array();
    for (std::size_t q = 0; q < st.genes.size(); q++) {
      const auto& gene = st.genes[q];
      const auto& cg = map.genes[q];
      const bool has = !cg.sites.empty();
      bool lead_causal = false;
      if (has) lead_causal = std::find(gene.causal.begin(), gene.causal.end(), cg.sites[cg.lead]) != gene.causal.end();
      genes.push(Json::object({{"name", gene.name}, {"chromosome", gene.chromosome}, {"tss", static_cast<double>(gene.tss)}, {"causal", gene.causal.size()},
                               {"cis_h2", da::round_to(gene.cis_h2)}, {"mediating", gene.mediating}, {"variants", cg.sites.size()},
                               {"p_nominal", has ? Json(cg.p_nominal) : Json()}, {"p_beta", has ? Json(cg.p_beta) : Json()}, {"p_permutation", has ? Json(cg.p_permutation) : Json()},
                               {"q", has ? Json(da::round_to(cg.q, 1e6)) : Json()}, {"egene", cg.egene}, {"beta_shape", Json::array({da::round_to(cg.beta_a), da::round_to(cg.beta_b, 100)})},
                               {"df_effective", da::round_to(cg.df_effective, 100)},
                               {"lead", has ? Json(static_cast<double>(g.positions[cg.sites[cg.lead]])) : Json()}, {"lead_causal", lead_causal}}));
      if (!has) continue;
      // Colocalisation of the eQTL with the trait over the gene's window.
      std::vector<double> se_e;
      for (const double v : cg.se) se_e.push_back(std::fabs(v));
      std::vector<double> bt, st_;
      for (const std::size_t j : cg.sites) {
        bt.push_back(scan.beta[j]);
        st_.push_back(scan.se[j]);
      }
      const auto c = assoc::coloc_abf(bt, st_, 0.15, cg.beta, se_e, 0.15);
      colocs.push(Json::object({{"gene", gene.name}, {"pp", numbers(c.posterior, 4)}, {"mediating", gene.mediating}, {"egene", cg.egene}}));
    }
    out.set("eqtl", Json::object({{"genes", genes}, {"pi0", da::round_to(map.pi0)}, {"egenes", map.egenes}, {"p_threshold", map.p_threshold},
                                  {"permutations", permutations}, {"window", static_cast<double>(eo.window)}}));
    out.set("coloc", colocs);
  }

  // ---- PrediXcan and the TWAS ----
  t = clock::now();
  assoc::PredictOptions po;
  // The nested cross-validation is the cost: (outer · (inner + 1) + inner + 1) paths a gene. PredictDB's
  // 5 × 10 folds and 100 penalties are cut to 3 × 4 folds, 12 penalties and a looser threshold (about a
  // sixth of this unit's earlier cost); request "predixcan_outer", "predixcan_folds", "predixcan_lambdas" to raise them.
  po.window = 100'000;
  po.outer_folds = static_cast<std::size_t>(bounded(request, "predixcan_outer", 3, 2, 10));
  po.folds = static_cast<std::size_t>(bounded(request, "predixcan_folds", 4, 3, 10));
  po.lambda_count = static_cast<std::size_t>(bounded(request, "predixcan_lambdas", 12, 5, 100));
  po.tolerance = 1e-5;
  po.seed = seed;
  const auto models = assoc::train_expression_models(g, st.reference, expr_ref, std::span<const assoc::GeneLocus>(loci), Cref, po, &worker.pool);
  seconds.set("predixcan", since(t));
  {
    Json genes = Json::array();
    Json scatter;
    double best_r2 = -1;
    for (std::size_t q = 0; q < st.genes.size(); q++) {
      const auto& m = models[q];
      const auto& gene = st.genes[q];
      Json row = Json::object({{"name", gene.name}, {"cis_h2", da::round_to(gene.cis_h2)}, {"mediating", gene.mediating}, {"kept", m.kept}, {"nonzero", m.nonzero},
                               {"r2_nested", da::round_to(m.r2_nested)}, {"rho", da::round_to(m.rho_mean)}, {"p_model", m.p}});
      if (m.kept) {
        const auto pred = assoc::predict_expression(m, g, st.sample);
        const auto truth = da::gather(std::span<const double>(&st.expression_genetic[q * st.individuals], st.individuals), st.sample);
        const auto observed = da::gather(st.gene_expression(q), st.sample);
        const double r = assoc::detail::correlation(pred, truth), ro = assoc::detail::correlation(pred, observed);
        const auto twas = assoc::associate(pred, y, C);
        const auto sp = assoc::s_predixcan(m, scan.statistic, g, st.sample);
        row.set("r2_genetic", da::round_to(r * r));
        row.set("r2_observed", da::round_to(ro * ro));
        row.set("twas_z", da::round_to(twas.z, 1e3));
        row.set("twas_p", twas.p);
        row.set("spredixcan_z", da::round_to(sp.z, 1e3));
        if (ro * ro > best_r2) {
          best_r2 = ro * ro;
          Json pts = Json::array();
          const std::size_t step = std::max<std::size_t>(1, pred.size() / 800);
          for (std::size_t i = 0; i < pred.size(); i += step) pts.push(Json::array({da::round_to(pred[i], 1e3), da::round_to(observed[i], 1e3), da::round_to(truth[i], 1e3)}));
          scatter = Json::object({{"gene", gene.name}, {"r2_observed", da::round_to(ro * ro)}, {"r2_genetic", da::round_to(r * r)}, {"cis_h2", da::round_to(gene.cis_h2)},
                                  {"has_eqtl", !gene.causal.empty()}, {"columns", Json::array({"predicted", "observed", "true_genetic"})}, {"points", pts}});
        }
      }
      genes.push(row);
    }
    out.set("predixcan", Json::object({{"genes", genes}, {"scatter", scatter}, {"alpha", po.alpha}, {"window", static_cast<double>(po.window)},
                                       {"outer_folds", po.outer_folds}, {"folds", po.folds}, {"lambdas", po.lambda_count}}));
  }

  // ---- heritability ----
  t = clock::now();
  {
    const std::size_t nh = std::min<std::size_t>(st.sample.size(), static_cast<std::size_t>(bounded(request, "heritability_samples", 1000, 200, 2500)));
    std::vector<std::size_t> who(st.sample.begin(), st.sample.begin() + static_cast<std::ptrdiff_t>(nh));
    const auto pw = assoc::principal_components(g, who, thin, 4);
    const auto Cw = da::covariates(st, who, &pw, 4);
    const auto h = assoc::snp_heritability(g, who, da::gather(st.trait, who), Cw, &worker.pool);
    const auto cc = da::gather(st.cases, who);
    const auto hc = assoc::snp_heritability(g, who, cc, Cw, &worker.pool);
    const double frac = std::accumulate(cc.begin(), cc.end(), 0.0) / static_cast<double>(cc.size());
    const double liab = assoc::liability_scale(hc.greml.h2, spec.prevalence, frac);
    Json profile = Json::array();
    for (std::size_t k = 0; k < h.greml.profile_h2.size(); k += 2)
      profile.push(Json::array({da::round_to(h.greml.profile_h2[k]), std::isfinite(h.greml.profile_log_likelihood[k]) ? Json(da::round_to(h.greml.profile_log_likelihood[k], 1e3)) : Json()}));
    out.set("heritability", Json::object({{"true", da::round_to(st.heritability)}, {"samples", nh},
                                          {"greml", da::round_to(h.greml.h2)}, {"greml_se", da::round_to(h.greml.h2_se)},
                                          {"he", da::round_to(h.he.h2)}, {"he_se", da::round_to(h.he.se)}, {"profile", profile},
                                          {"case_control", Json::object({{"observed", da::round_to(hc.greml.h2)}, {"liability", da::round_to(liab)},
                                                                         {"liability_se", da::round_to(hc.greml.h2_se * liab / std::max(1e-12, hc.greml.h2))},
                                                                         {"prevalence", spec.prevalence}})}}));
  }
  seconds.set("heritability", since(t));

  // ---- polygenic scores ----
  t = clock::now();
  {
    const std::size_t half = st.sample.size() * 2 / 3;
    std::vector<std::size_t> train(st.sample.begin(), st.sample.begin() + static_cast<std::ptrdiff_t>(half)), test(st.sample.begin() + static_cast<std::ptrdiff_t>(half), st.sample.end());
    const auto pt = assoc::principal_components(g, train, thin, 4), pv = assoc::principal_components(g, test, thin, 4);
    const auto Ct = da::covariates(st, train, &pt, 4), Cv = da::covariates(st, test, &pv, 4);
    const auto yt = da::gather(st.trait, train), yv = da::gather(st.trait, test);
    const auto tscan = assoc::linear_scan(g, train, yt, Ct, &worker.pool);
    Json ct = Json::array();
    const std::vector<double> thresholds{5e-8, 1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.5, 1};
    for (const auto& s : assoc::clump_and_threshold(g, train, tscan, thresholds))
      ct.push(Json::object({{"threshold", s.threshold}, {"sites", s.sites.size()}, {"r2", da::round_to(assoc::score_r2(assoc::polygenic_score(g, test, s.sites, s.weights), yv, Cv))}}));
    const auto ld = assoc::ldpred_inf(g, train, yt, Ct, st.heritability, 300, &worker.pool);
    const double r2_ld = assoc::score_r2(assoc::polygenic_score(g, test, ld.sites, ld.weights), yv, Cv);
    // Daetwyler et al. (2008)'s guide, with M_e the sites.
    const double Nh = static_cast<double>(train.size()) * st.heritability;
    out.set("prs", Json::object({{"train", train.size()}, {"test", test.size()}, {"ct", ct}, {"ldpred_inf", da::round_to(r2_ld)},
                                 {"guide", da::round_to(st.heritability * Nh / (Nh + static_cast<double>(g.sites)))}}));
  }
  seconds.set("prs", since(t));
  out.set("seconds", seconds);
  out.set("threads", worker.pool.size());
  return out;
}

}  // namespace workbench
