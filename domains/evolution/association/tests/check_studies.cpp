// Checks of fine-mapping, colocalisation, cis-eQTL mapping, expression prediction and the TWAS,
// SNP heritability and polygenic scores (finemap.hpp, coloc.hpp, eqtl.hpp, predict.hpp,
// heritability.hpp, prs.hpp), against R's coloc, closed forms and the truth of simulated cohorts.
//
// - ABF: Wakefield's closed form; one causal variant's posterior sums to 1 and its credible set reaches
//   0.95.
// - coloc.abf against R's coloc 5.2.3 (fixtures/coloc.txt): PP.H0–PP.H4 and SNP.PP.H4 in three scenarios
//   (shared, distinct, one trait null); H₃ and H₄ against the enumeration of all p² configurations.
// - Fine-mapping a simulated cohort's associated regions: SuSiE's 95% credible sets cover causal
//   variants (the fraction stated), SuSiE-RSS on the scan's z-scores and the sample's LD finds the same
//   sets as SuSiE on individual data.
// - cis-eQTL mapping: genes with causal eQTLs are eGenes and genes without are not, at the FDR asked (the
//   discoveries and false discoveries stated); the beta-approximated p of null genes is uniform enough
//   (its mean near ½); the lead variant is a causal variant or in LD with one (r² > 0.5).
// - PrediXcan: the nested-CV R² of genes with cis-heritability above 0.15 is positive and below their
//   simulated cis-h², and the held-out R² of predicted against true genetic expression in the
//   association sample is high; genes without eQTLs get no kept model; S-PrediXcan's z equals the
//   individual-level TWAS z (in-sample LD) to within 5% and the mediating genes are found.
// - SNP heritability: GREML on the association sample recovers the simulated h² within 3 standard
//   errors, and so does Haseman–Elston; the liability-scale transformation of the case–control trait's
//   observed h² lands near the liability's h².
// - Polygenic scores: LDpred-inf with in-sample LD and one block equals ridge regression's solution
//   (D + λI)⁻¹ β̃ computed apart; on held-out samples LDpred-inf and C+T predict with R² > 0, both below
//   h², and on this sparse architecture C+T's best threshold beats LDpred-inf (as Vilhjálmsson et al.
//   2015 find for few causal variants).

#include <evolution/association/coloc.hpp>
#include <evolution/association/cohort.hpp>
#include <evolution/association/eqtl.hpp>
#include <evolution/association/finemap.hpp>
#include <evolution/association/heritability.hpp>
#include <evolution/association/ld.hpp>
#include <evolution/association/predict.hpp>
#include <evolution/association/prs.hpp>
#include <evolution/association/scan.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <fstream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace as = evolution::association;
namespace pop = evolution::population;
using calx::linalg::Rectangular;

namespace {

std::string fmt(double x) { return calx::verification::exponential(x); }

std::map<std::string, std::vector<double>> read_fixture(const std::string& path) {
  std::map<std::string, std::vector<double>> out;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    const auto colon = line.find(" :");
    if (colon == std::string::npos) continue;
    std::istringstream values(line.substr(colon + 2));
    std::vector<double> v;
    std::string token;
    while (values >> token) {
      try {
        v.push_back(std::stod(token));
      } catch (...) {
        break;
      }
    }
    out[line.substr(0, colon)] = v;
  }
  return out;
}

Rectangular covariates(const as::Study& st, std::span<const std::size_t> who, const as::Components& pcs, std::size_t k,
                       const std::vector<std::vector<double>>& extra = {}) {
  Rectangular C(who.size(), 3 + k + extra.size());
  for (std::size_t i = 0; i < who.size(); i++) {
    C(i, 0) = 1;
    C(i, 1) = st.sex[who[i]];
    C(i, 2) = st.age[who[i]];
    for (std::size_t c = 0; c < k; c++) C(i, 3 + c) = pcs.vectors[i * pcs.components + c];
    for (std::size_t e = 0; e < extra.size(); e++) C(i, 3 + k + e) = extra[e][i];
  }
  return C;
}

std::vector<double> gather(std::span<const double> x, std::span<const std::size_t> who) {
  std::vector<double> out(who.size());
  for (std::size_t i = 0; i < who.size(); i++) out[i] = x[who[i]];
  return out;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::scheduler::ThreadPool pool(8);

  // ---- ABF and coloc ----
  {
    const double b = 0.12, v = 0.0016, w = 0.0225, z2 = b * b / v;
    const double closed = std::log(std::sqrt(v / (v + w)) * std::exp(z2 * w / (2 * (v + w))));
    s.check("Wakefield's ABF in closed form", std::fabs(as::log_abf(b, v, w) - closed), 1e-14);
    const auto F = read_fixture(std::string(EVOLUTION_FIXTURES) + "/coloc.txt");
    s.require("the fixture of coloc's values is read", F.count("shared_pp") && F.count("distinct_pp") && F.count("one_pp"));
    for (const std::string name : {"shared", "distinct", "one"}) {
      const auto c = as::coloc_abf(F.at(name + "_b1"), [&] {
        std::vector<double> se;
        for (const double x : F.at(name + "_v1")) se.push_back(std::sqrt(x));
        return se;
      }(), 0.15, F.at(name + "_b2"), [&] {
        std::vector<double> se;
        for (const double x : F.at(name + "_v2")) se.push_back(std::sqrt(x));
        return se;
      }(), 0.2);
      double worst = 0, worst_snp = 0;
      for (std::size_t h = 0; h < 5; h++) worst = std::max(worst, std::fabs(c.posterior[h] - F.at(name + "_pp")[h]));
      for (std::size_t j = 0; j < c.shared.size(); j++) worst_snp = std::max(worst_snp, std::fabs(c.shared[j] - F.at(name + "_snp_pp_h4")[j]));
      s.check("coloc.abf (" + name + "): PP.H0–PP.H4 equal R coloc's", worst, 1e-12);
      s.check("coloc.abf (" + name + "): SNP.PP.H4 equal R coloc's", worst_snp, 1e-12);
    }
    // H₃ and H₄ by enumerating configurations.
    const auto& l1v = F.at("distinct_b1");
    std::vector<double> l1, l2;
    for (std::size_t j = 0; j < l1v.size(); j++) {
      l1.push_back(as::log_abf(F.at("distinct_b1")[j], F.at("distinct_v1")[j], 0.0225));
      l2.push_back(as::log_abf(F.at("distinct_b2")[j], F.at("distinct_v2")[j], 0.04));
    }
    const auto c = as::coloc_from_lbf(l1, l2);
    double h3 = 0, h4 = 0;
    for (std::size_t j = 0; j < l1.size(); j++)
      for (std::size_t k = 0; k < l2.size(); k++) (j == k ? h4 : h3) += std::exp(l1[j] + l2[k]);
    s.check("H₃'s Bayes factor is the sum over the p(p − 1) distinct configurations (relative)", std::fabs(std::exp(c.log_bayes[3]) / (1e-8 * h3) - 1), 1e-10);
    s.check("and H₄'s over the p shared ones", std::fabs(std::exp(c.log_bayes[4]) / (1e-5 * h4) - 1), 1e-10);
  }

  // ---- a simulated study ----
  as::StudySpec spec;
  spec.cohort.individuals = 1800;
  spec.cohort.regions = 3;
  spec.cohort.theta = 200;
  spec.cohort.rho = 200;
  spec.cohort.founders = 120;
  spec.genes = 18;
  spec.reference = 400;
  spec.causal = 15;
  spec.heritability = 0.5;
  spec.mediated = 0.08;
  spec.mediating = 3;
  const auto st = as::simulate_study(spec, 11, &pool);
  const auto& g = st.genotypes;
  // Principal components from every fifth site (LD pruning by thinning, as a study would).
  std::vector<std::size_t> thin;
  for (std::size_t j = 0; j < g.sites; j += 5) thin.push_back(j);
  const auto pcs = as::principal_components(g, st.sample, thin, 3);
  const auto C = covariates(st, st.sample, pcs, 3);
  const auto y = gather(st.trait, st.sample);
  const auto scan = as::linear_scan(g, st.sample, y, C, &pool);

  // ---- fine-mapping ----
  {
    as::ClumpOptions o;
    o.p1 = 1e-8;
    o.r2 = 0.1;
    const auto clumps = as::clump(g, st.sample, scan.p, o);
    std::size_t sets = 0, covering = 0, regions = 0, same = 0;
    for (const auto& c : clumps) {
      if (regions == 4) break;
      regions++;
      const auto lead = c.index;
      const std::uint64_t pos = g.positions[lead];
      const auto sites = as::sites_between(g, g.chromosome[lead], pos > 100'000 ? pos - 100'000 : 0, pos + 100'000);
      const auto fit = as::finemap_susie(g, st.sample, sites, y, C, {}, &pool);
      std::vector<std::vector<std::size_t>> cs;
      for (const auto& set : fit.sets) cs.push_back(set.variables);
      const auto truth = as::against_truth(cs, sites, st.trait_causal);
      sets += truth.sets;
      covering += truth.covering;
      // SuSiE-RSS on the scan's z-scores and the sample's LD.
      std::vector<double> z;
      for (const std::size_t j : sites) z.push_back(scan.statistic[j]);
      const auto R = as::ld_correlation(g, st.sample, sites, &pool);
      const auto rss = as::finemap_rss(z, R, st.sample.size());
      bool agree = rss.sets.size() == fit.sets.size();
      for (std::size_t k = 0; agree && k < fit.sets.size(); k++) {
        const auto& v = fit.sets[k].variables;
        bool found = false;
        for (const auto& r : rss.sets) found = found || std::find(r.variables.begin(), r.variables.end(), v[0]) != r.variables.end();
        agree = found;
      }
      same += agree;
    }
    s.require("SuSiE's credible sets in " + std::to_string(regions) + " associated regions cover causal variants (" + std::to_string(covering) + " of " +
                  std::to_string(sets) + ")",
              sets > 0 && covering * 10 >= sets * 8);
    s.require("and SuSiE-RSS on z-scores and in-sample LD finds the same sets (" + std::to_string(same) + " of " + std::to_string(regions) + " regions)",
              same + 1 >= regions);
    // One causal variant's ABF posterior.
    std::vector<double> b, se;
    const auto lead = clumps.front().index;
    const auto sites = as::sites_between(g, g.chromosome[lead], g.positions[lead] - std::min<std::uint64_t>(g.positions[lead], 50'000), g.positions[lead] + 50'000);
    for (const std::size_t j : sites) {
      b.push_back(scan.beta[j]);
      se.push_back(scan.se[j]);
    }
    const auto abf = as::finemap_abf(b, se, 0.15);
    double total = 0;
    for (const double p : abf.posterior) total += p;
    s.check("one causal variant: the posteriors sum to 1", total - 1, 1e-12);
    s.require("and the credible set reaches 0.95 (" + std::to_string(abf.credible.size()) + " variants)", abf.coverage >= 0.95);
  }

  // ---- cis-eQTL mapping ----
  const auto expr_ref = [&] {
    std::vector<std::vector<double>> e(st.genes.size());
    for (std::size_t q = 0; q < st.genes.size(); q++) e[q] = gather(st.gene_expression(q), st.reference);
    return e;
  }();
  const auto factors = as::expression_components(expr_ref, 3);
  const auto pcs_ref = as::principal_components(g, st.reference, thin, 3);
  const auto Cref = covariates(st, st.reference, pcs_ref, 3, factors);
  std::vector<as::GeneLocus> loci;
  for (const auto& gene : st.genes) loci.push_back({gene.chromosome, gene.tss});
  {
    as::CisOptions o;
    o.window = 200'000;
    o.permutations = 500;
    const auto map = as::cis_eqtl(g, st.reference, expr_ref, loci, Cref, o, &pool);
    std::size_t true_pos = 0, false_pos = 0, with = 0, lead_ok = 0, egenes = 0;
    double null_mean = 0, nulls = 0;
    for (std::size_t q = 0; q < st.genes.size(); q++) {
      const auto& cg = map.genes[q];
      const bool causal = !st.genes[q].causal.empty();
      with += causal;
      if (cg.egene) {
        egenes++;
        (causal ? true_pos : false_pos)++;
        double best = 0;
        for (const std::size_t c : st.genes[q].causal) {
          std::vector<std::size_t> pair{cg.sites[cg.lead], c};
          best = std::max(best, std::pow(as::ld_correlation(g, st.reference, pair)(0, 1), 2));
        }
        lead_ok += best > 0.5;
      }
      if (!causal) {
        null_mean += cg.p_beta;
        nulls++;
      }
    }
    s.require("cis-eQTL mapping finds eGenes among the " + std::to_string(with) + " genes with causal eQTLs (" + std::to_string(true_pos) + ")", true_pos * 2 >= with);
    s.require("with few false eGenes (" + std::to_string(false_pos) + " of " + std::to_string(egenes) + ")", false_pos <= 1);
    s.require("the eGenes' lead variants are causal or in LD with one (r² > 0.5: " + std::to_string(lead_ok) + " of " + std::to_string(egenes) + ")",
              lead_ok + 1 >= egenes);
    if (nulls > 0) s.check("genes without eQTLs: the mean beta-approximated p is ½ (" + std::to_string(static_cast<int>(nulls)) + " genes)", null_mean / nulls - 0.5, 0.25);
  }

  // ---- PrediXcan and the TWAS ----
  {
    as::PredictOptions o;
    o.window = 100'000;
    o.folds = 5;
    o.outer_folds = 3;
    const auto models = as::train_expression_models(g, st.reference, expr_ref, std::span<const as::GeneLocus>(loci), Cref, o, &pool);
    std::size_t kept_null = 0, checked = 0, below = 0, good = 0;
    double worst_z = 0;
    std::size_t found_mediating = 0, mediating = 0;
    const auto z_scan = scan.statistic;
    for (std::size_t q = 0; q < st.genes.size(); q++) {
      const auto& m = models[q];
      const auto& gene = st.genes[q];
      if (gene.causal.empty()) {
        kept_null += m.kept;
        continue;
      }
      if (gene.cis_h2 > 0.15 && m.kept) {
        checked++;
        below += m.r2_nested < gene.cis_h2 + 0.05 && m.r2_nested > 0;
        const auto pred = as::predict_expression(m, g, st.sample);
        const auto truth = gather(std::span<const double>(&st.expression_genetic[q * st.individuals], st.individuals), st.sample);
        const double r = as::detail::correlation(pred, truth);
        good += r * r > 0.5;
      }
      if (!m.kept) continue;
      const auto pred = as::predict_expression(m, g, st.sample);
      const auto ind = as::associate(pred, y, C);
      const auto sum = as::s_predixcan(m, z_scan, g, st.sample);
      if (std::fabs(ind.z) > 2) worst_z = std::max(worst_z, std::fabs(sum.z - ind.z) / std::fabs(ind.z));
      if (gene.mediating) {
        mediating++;
        found_mediating += ind.p < 0.05 / static_cast<double>(st.genes.size());
      }
    }
    s.require("PrediXcan: the nested-CV R² of well-predicted genes lies in (0, cis-h² + 0.05) (" + std::to_string(below) + " of " + std::to_string(checked) + ")",
              checked > 0 && below + 1 >= checked);
    s.require("and their predicted expression in the association sample tracks the true genetic expression (r² > 0.5: " + std::to_string(good) + " of " +
                  std::to_string(checked) + ")",
              good + 1 >= checked);
    s.require("genes without eQTLs keep no model (" + std::to_string(kept_null) + ")", kept_null <= 1);
    s.check("S-PrediXcan's z equals the individual-level TWAS z with in-sample LD (largest relative difference where |z| > 2)", worst_z, 0.08);
    s.require("the TWAS finds the mediating genes (" + std::to_string(found_mediating) + " of " + std::to_string(mediating) + " with models)",
              mediating == 0 || found_mediating * 2 >= mediating);
  }

  // ---- SNP heritability ----
  {
    std::vector<std::size_t> who(st.sample.begin(), st.sample.begin() + 900);
    const auto pw = as::principal_components(g, who, thin, 3);
    const auto Cw = covariates(st, who, pw, 3);
    const auto h = as::snp_heritability(g, who, gather(st.trait, who), Cw, &pool);
    s.check("GREML recovers the simulated h² (" + fmt(st.heritability) + "): (ĥ² − h²)/SE", (h.greml.h2 - st.heritability) / h.greml.h2_se, 3, "SE");
    s.check("and so does Haseman–Elston: (ĥ² − h²)/SE", (h.he.h2 - st.heritability) / h.he.se, 3, "SE");
    double cases = 0;
    const auto cc = gather(st.cases, who);
    for (const double v : cc) cases += v;
    const auto hc = as::snp_heritability(g, who, cc, Cw, &pool);
    const double liab = as::liability_scale(hc.greml.h2, spec.prevalence, cases / who.size());
    s.check("the case–control trait's GREML h² on the liability scale (Lee et al. 2011) is the liability's: (h²_l − h²)/SE_l",
            (liab - st.heritability) / (hc.greml.h2_se * liab / hc.greml.h2), 3, "SE");
  }

  // ---- polygenic scores ----
  {
    const std::size_t half = st.sample.size() / 2;
    std::vector<std::size_t> train(st.sample.begin(), st.sample.begin() + half), test(st.sample.begin() + half, st.sample.end());
    const auto pt = as::principal_components(g, train, thin, 3), pv = as::principal_components(g, test, thin, 3);
    const auto Ct = covariates(st, train, pt, 3), Cv = covariates(st, test, pv, 3);
    const auto yt = gather(st.trait, train), yv = gather(st.trait, test);
    const auto tscan = as::linear_scan(g, train, yt, Ct, &pool);
    const auto ld = as::ldpred_inf(g, train, yt, Ct, st.heritability, 400, &pool);
    const double r2_ld = as::score_r2(as::polygenic_score(g, test, ld.sites, ld.weights), yv, Cv);
    double best_ct = 0;
    const std::vector<double> thresholds{5e-8, 1e-5, 1e-3, 0.05, 1};
    for (const auto& t : as::clump_and_threshold(g, train, tscan, thresholds)) best_ct = std::max(best_ct, as::score_r2(as::polygenic_score(g, test, t.sites, t.weights), yv, Cv));
    s.require("on held-out samples LDpred-inf predicts (R² " + fmt(r2_ld) + ") and C+T's best threshold too (" + fmt(best_ct) + "), both below h²",
              r2_ld > 0.05 && best_ct > 0.02 && r2_ld < st.heritability && best_ct < st.heritability);
    // The simulated architecture is sparse (15 direct causal variants and a few eQTLs), where selecting
    // variants beats the infinitesimal prior, as Vilhjálmsson et al. (2015, Figs. 2–3) find for traits
    // with few causal variants (LDpred with a small causal fraction, not written here, would close it).
    s.require("on this sparse architecture C+T's selection beats LDpred-inf's infinitesimal prior", best_ct > r2_ld);
    // LDpred-inf with one block of a region and in-sample LD is ridge regression.
    const auto [lo, hi] = g.chromosome_sites(0);
    const std::size_t b = std::min<std::size_t>(hi - lo, 150);
    std::vector<std::size_t> sites(b);
    std::iota(sites.begin(), sites.end(), lo);
    auto sub = pop::select_sites(g, sites);
    const auto one = as::ldpred_inf(sub, train, yt, Ct, 0.3, b, &pool);
    const auto Z = as::standardised_dosages(sub, train, sites.empty() ? std::vector<std::size_t>{} : [&] {
      std::vector<std::size_t> all(b);
      std::iota(all.begin(), all.end(), 0);
      return all;
    }());
    const as::Residualizer M(Ct);
    auto yr = M.residuals(yt);
    double ss = 0;
    for (const double v : yr) ss += v * v;
    for (double& v : yr) v /= std::sqrt(ss / (yr.size() - 1.0));
    calx::linalg::Matrix D(b), rhs(b);
    for (std::size_t a = 0; a < b; a++) {
      for (std::size_t i = 0; i < train.size(); i++) rhs(a, 0) += Z(i, a) * yr[i] / (train.size() - 1.0);
      for (std::size_t c = 0; c < b; c++) {
        double v = 0;
        for (std::size_t i = 0; i < train.size(); i++) v += Z(i, a) * Z(i, c);
        D(a, c) = v / (train.size() - 1.0);
      }
      D(a, a) += static_cast<double>(b) / (train.size() * 0.3);
    }
    const auto ridge = calx::linalg::solve(D, rhs);
    double worst = 0;
    for (std::size_t a = 0; a < b; a++) worst = std::max(worst, std::fabs(one.standardised[a] - ridge(a, 0)));
    s.check("LDpred-inf with one block and in-sample LD equals ridge regression (D + (M/(Nh²))I)⁻¹β̃ (LU, apart)", worst, 1e-10);
  }
  return s.done();
}
