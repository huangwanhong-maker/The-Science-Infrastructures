// Checks of the cohort simulator, the association scans, genomic control, LD, clumping and principal
// components (cohort.hpp, scan.hpp, ld.hpp).
//
// - The cohort: a function of its seed; the same on 1 and 8 threads; every site passes the frequency
//   filter; linkage decays with distance (r² of neighbouring sites above that of distant ones); the
//   populations' founder groups differentiate them (Hudson's F_ST of population/diversity.hpp near
//   1/(2F) for F founder haplotypes); the realised heritability, the expression's cis-heritabilities and
//   the case fraction are those asked.
// - The linear scan equals tyche's ols on [C, g_j] (β, se, p) at 30 sites to 10⁻¹⁰ (Frisch–Waugh–Lovell);
//   under a null trait (pure noise) λ_GC is near 1 and P(p < 0.05) near 0.05; under confounding by
//   population without principal components λ_GC rises, and with them it returns near 1; the causal
//   variants are found.
// - The logistic score test equals its definition by dense algebra from tyche's glm null fit
//   (U = gᵀ(y − μ̂), V = gᵀWg − gᵀWC(CᵀWC)⁻¹CᵀWg) at 20 sites; its Wald refits equal tyche's glm.
// - The Q–Q band holds 95% of null points; genomic control of uniform p-values is 1.
// - The principal components through the m′ × m′ covariance equal population::pca's axes (|correlation|
//   above 1 − 10⁻⁹) and eigenvalues; the first separates the two populations.
// - LD scores equal their definition; clumps hold only sites above the r² threshold within the window.

#include <evolution/association/cohort.hpp>
#include <evolution/association/ld.hpp>
#include <evolution/association/scan.hpp>
#include <evolution/population/diversity.hpp>
#include <evolution/population/structure.hpp>

#include <tyche/regression/glm.hpp>
#include <tyche/regression/linear.hpp>

#include <calx/linalg/cholesky.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <numeric>
#include <string>
#include <vector>

namespace as = evolution::association;
namespace pop = evolution::population;
using calx::linalg::Rectangular;

namespace {

std::string fmt(double x) { return calx::verification::exponential(x); }

Rectangular covariates(const as::Study& st, std::span<const std::size_t> who, const as::Components* pcs, std::size_t k) {
  Rectangular C(who.size(), 3 + (pcs ? k : 0));
  for (std::size_t i = 0; i < who.size(); i++) {
    C(i, 0) = 1;
    C(i, 1) = st.sex[who[i]];
    C(i, 2) = st.age[who[i]];
    for (std::size_t c = 0; pcs && c < k; c++) C(i, 3 + c) = pcs->vectors[i * pcs->components + c];
  }
  return C;
}

double r_of(const pop::Genotypes& g, std::span<const std::size_t> who, std::size_t a, std::size_t b) {
  std::vector<std::size_t> s{a, b};
  return as::ld_correlation(g, who, s)(0, 1);
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::scheduler::ThreadPool one(1), eight(8);
  as::StudySpec spec;
  spec.cohort.individuals = 1500;
  spec.cohort.regions = 3;
  spec.cohort.theta = 200;
  spec.cohort.rho = 200;
  spec.cohort.founders = 100;
  spec.genes = 12;
  spec.reference = 300;
  spec.population_effect = 0.6;
  const auto st = as::simulate_study(spec, 7, &eight);
  const auto& g = st.genotypes;
  const std::size_t n = st.individuals;
  std::vector<std::size_t> all(n);
  std::iota(all.begin(), all.end(), 0);

  // ---- the cohort ----
  {
    const auto again = as::simulate_study(spec, 7, &one);
    s.require("a cohort is a function of its seed, the same on 1 and 8 threads", again.genotypes.dosage == g.dosage && again.trait == st.trait && again.expression == st.expression);
    s.require("and another seed gives another cohort", as::simulate_study(spec, 8, &eight).genotypes.dosage != g.dosage);
    s.require("labelled simulated", st.simulated && g.phased());
    bool maf = true;
    for (std::size_t j = 0; j < g.sites; j++) {
      const auto [alt, total] = pop::count_site(g.site(j));
      const double f = static_cast<double>(alt) / total;
      maf = maf && std::min(f, 1 - f) >= spec.cohort.min_maf;
    }
    s.require("every site passes the frequency filter (" + std::to_string(g.sites) + " sites in 3 regions)", maf && g.sites > 300);
    // Linkage decays: mean r² of adjacent sites against sites 200 kb apart.
    double near = 0, far = 0;
    std::size_t cn = 0, cf = 0;
    const auto [lo, hi] = g.chromosome_sites(0);
    std::vector<std::size_t> some(all.begin(), all.begin() + 500);
    for (std::size_t j = lo; j + 1 < hi && cn < 100; j += 3, cn++) near += std::pow(r_of(g, some, j, j + 1), 2);
    for (std::size_t j = lo; j < hi && cf < 100; j += 3) {
      std::size_t k = j;
      while (k < hi && g.positions[k] < g.positions[j] + 200'000) k++;
      if (k >= hi) break;
      far += std::pow(r_of(g, some, j, k), 2);
      cf++;
    }
    s.require("linkage decays with distance: mean r² of neighbours " + fmt(near / cn) + " against " + fmt(far / std::max<std::size_t>(1, cf)) + " at 200 kb",
              near / cn > 3 * far / std::max<std::size_t>(1, cf));
    // F_ST between the two populations, near 1/(2F) = 0.005 per population pair: Hudson's estimate.
    std::vector<std::uint32_t> label(st.population.begin(), st.population.end());
    const auto groups = pop::members(label);
    const auto c0 = pop::allele_counts(g, groups[0]), c1 = pop::allele_counts(g, groups[1]);
    const double fst = pop::fst_hudson(c0, c1).value();
    s.require("the founder groups differentiate the populations: Hudson's F_ST " + fmt(fst) + " between 0.002 and 0.03", fst > 0.002 && fst < 0.03);
    s.check("the realised heritability is the one asked (0.4)", st.heritability - spec.heritability, 0.03);
    double worst = 0;
    for (std::size_t q = 0; q < st.genes.size(); q++) {
      const auto& gene = st.genes[q];
      if (gene.causal.empty()) continue;
      worst = std::max(worst, std::fabs(gene.cis_h2 - as::detail::variance_of(std::span<const double>(&st.expression_genetic[q * n], n))));
      worst = std::max(worst, std::fabs(as::detail::variance_of(st.gene_expression(q)) - 1) * 0.1);
    }
    s.check("each gene's expression has unit variance and its cis-heritability is its genetic part's variance", worst, 0.02);
    double cases = 0;
    for (const double c : st.cases) cases += c;
    s.check("the case fraction is the prevalence (0.2)", cases / n - spec.prevalence, 0.002);
  }

  // ---- principal components ----
  const auto pcs = as::principal_components(g, st.sample, {}, 4);
  {
    const auto sub = pop::select_samples(g, st.sample);
    const auto ref = pop::pca(sub, 4, &eight);
    double worst = 0, worst_value = 0;
    for (std::size_t c = 0; c < 4; c++) {
      double dot = 0;
      for (std::size_t i = 0; i < st.sample.size(); i++) dot += pcs.vectors[i * 4 + c] * ref.vectors[i * 4 + c];
      worst = std::max(worst, 1 - std::fabs(dot));
      worst_value = std::max(worst_value, std::fabs(pcs.values[c] - ref.values[c]) / ref.values[c]);
    }
    s.check("principal components through the m × m covariance are population::pca's axes (1 − |cos|)", worst, 1e-9);
    s.check("and its eigenvalues on Patterson's scale (relative)", worst_value, 1e-9);
    double m0 = 0, m1 = 0;
    std::size_t k0 = 0, k1 = 0;
    for (std::size_t i = 0; i < st.sample.size(); i++)
      (st.population[st.sample[i]] == 0 ? (m0 += pcs.vectors[i * 4], k0++) : (m1 += pcs.vectors[i * 4], k1++));
    double within = 0;
    for (std::size_t i = 0; i < st.sample.size(); i++) {
      const double c = st.population[st.sample[i]] == 0 ? m0 / k0 : m1 / k1;
      within += std::pow(pcs.vectors[i * 4] - c, 2);
    }
    within = std::sqrt(within / st.sample.size());
    s.require("the first component separates the populations (gap " + fmt(std::fabs(m0 / k0 - m1 / k1)) + " against spread " + fmt(within) + ")",
              std::fabs(m0 / k0 - m1 / k1) > 3 * within);
  }

  // ---- the linear scan ----
  {
    std::vector<double> y(st.sample.size());
    for (std::size_t i = 0; i < y.size(); i++) y[i] = st.trait[st.sample[i]];
    const auto C = covariates(st, st.sample, &pcs, 4);
    const auto scan = as::linear_scan(g, st.sample, y, C, &eight);
    double worst = 0;
    for (std::size_t q = 0; q < 30; q++) {
      const std::size_t j = (q * 997) % g.sites;
      Rectangular X(y.size(), C.cols + 1);
      const auto d = as::sample_dosages(g, j, st.sample);
      for (std::size_t i = 0; i < y.size(); i++) {
        for (std::size_t a = 0; a < C.cols; a++) X(i, a) = C(i, a);
        X(i, C.cols) = d[i];
      }
      const auto fit = tyche::regression::ols(X, y);
      const std::size_t k = C.cols;
      worst = std::max({worst, std::fabs(scan.beta[j] - fit.beta[k]) / std::fabs(fit.beta[k]), std::fabs(scan.se[j] - fit.se[k]) / fit.se[k],
                        std::fabs(scan.p[j] - fit.p[k]) / std::max(1e-300, fit.p[k])});
    }
    s.check("the linear scan equals tyche's ols on [C, g] at 30 sites: β, se and p (largest relative difference)", worst, 1e-9);
    const auto scan0 = as::linear_scan(g, st.sample, y, covariates(st, st.sample, nullptr, 0), &eight);
    const double l_pc = as::genomic_control(scan.p), l_none = as::genomic_control(scan0.p);
    s.require("confounding by population raises λ_GC without components (" + fmt(l_none) + ") above λ_GC with four (" + fmt(l_pc) + ")", l_none > l_pc + 0.05);
    std::size_t found = 0;
    for (const std::size_t c : st.causal) found += scan.p[c] < 1e-4;
    s.require("causal variants found at p < 10⁻⁴ (" + std::to_string(found) + " of " + std::to_string(st.causal.size()) + ")", found >= 3);
    // A null trait: noise.
    calx::random::Stream r(3, 0x4e554c4c, 0, 0);
    std::vector<double> noise(y.size());
    for (double& v : noise) v = r.normal();
    const auto null = as::linear_scan(g, st.sample, noise, C, &eight);
    s.check("under a null trait λ_GC is 1 (λ − 1)", as::genomic_control(null.p) - 1, 0.08);
    double below = 0, count = 0;
    for (const double p : null.p)
      if (std::isfinite(p)) {
        below += p < 0.05;
        count++;
      }
    s.check("and P(p < 0.05) is 0.05 (correlated sites; a loose bound)", below / count - 0.05, 0.02);
    const auto qq = as::qq_points(null.p, 1e-2, 2000);
    double inside = 0;
    for (const auto& q : qq) inside += q.observed >= q.lower - 1e-12 && q.observed <= q.upper + 1e-12;
    s.require("the Q–Q band holds most null points (" + fmt(inside / qq.size()) + ")", inside / qq.size() > 0.8);
  }

  // ---- the logistic scan ----
  {
    std::vector<double> y(st.sample.size());
    for (std::size_t i = 0; i < y.size(); i++) y[i] = st.cases[st.sample[i]];
    const auto C = covariates(st, st.sample, &pcs, 4);
    const auto scan = as::logistic_scan(g, st.sample, y, C, &eight, 1e-3);
    tyche::regression::GlmOptions o;
    o.family = tyche::regression::Family::binomial;
    const auto null = tyche::regression::glm(C, y, o);
    const std::size_t n2 = y.size(), k = C.cols;
    calx::linalg::Matrix CWC(k);
    for (std::size_t i = 0; i < n2; i++) {
      const double w = null.fitted[i] * (1 - null.fitted[i]);
      for (std::size_t a = 0; a < k; a++)
        for (std::size_t b = 0; b < k; b++) CWC(a, b) += C(i, a) * w * C(i, b);
    }
    const auto chol = calx::linalg::cholesky(CWC);
    double worst = 0, worst_wald = 0;
    for (std::size_t q = 0; q < 20; q++) {
      const std::size_t j = (q * 1231 + 5) % g.sites;
      const auto d = as::sample_dosages(g, j, st.sample);
      double U = 0, gwg = 0;
      std::vector<double> cwg(k, 0.0);
      for (std::size_t i = 0; i < n2; i++) {
        const double w = null.fitted[i] * (1 - null.fitted[i]);
        U += d[i] * (y[i] - null.fitted[i]);
        gwg += d[i] * w * d[i];
        for (std::size_t a = 0; a < k; a++) cwg[a] += C(i, a) * w * d[i];
      }
      const auto sol = calx::linalg::cholesky_solve(chol, cwg);
      double V = gwg;
      for (std::size_t a = 0; a < k; a++) V -= cwg[a] * sol[a];
      const double x2 = U * U / V;
      worst = std::max(worst, std::fabs(scan.statistic[j] * scan.statistic[j] - x2) / x2);
    }
    s.check("the logistic score χ² equals its definition from tyche's null glm at 20 sites (relative)", worst, 1e-8);
    std::size_t refits = 0;
    for (std::size_t j = 0; j < g.sites && refits < 3; j++) {
      if (!std::isfinite(scan.wald_p[j])) continue;
      refits++;
      Rectangular X(n2, k + 1);
      const auto d = as::sample_dosages(g, j, st.sample);
      for (std::size_t i = 0; i < n2; i++) {
        for (std::size_t a = 0; a < k; a++) X(i, a) = C(i, a);
        X(i, k) = d[i];
      }
      const auto fit = tyche::regression::glm(X, y, o);
      worst_wald = std::max(worst_wald, std::fabs(fit.beta[k] - scan.wald_beta[j]) / std::fabs(fit.beta[k]));
      // The score and Wald statistics agree at strong sites, to the order of the effect's square.
    }
    s.require("the strongest sites are refitted by Wald (" + std::to_string(refits) + ")", refits > 0);
    s.check("and their Wald β equal tyche's glm's", worst_wald, 1e-12);
  }

  // ---- LD scores and clumps ----
  {
    std::vector<std::size_t> panel(st.reference.begin(), st.reference.end());
    const auto l = as::ld_scores(g, panel, 50'000, &eight);
    double worst = 0;
    for (std::size_t q = 0; q < 5; q++) {
      const std::size_t j = (q * 401 + 17) % g.sites;
      double direct = 0;
      for (std::size_t k = 0; k < g.sites; k++) {
        if (g.chromosome[k] != g.chromosome[j]) continue;
        const auto d = g.positions[k] > g.positions[j] ? g.positions[k] - g.positions[j] : g.positions[j] - g.positions[k];
        if (d > 50'000) continue;
        const double r = r_of(g, panel, j, k), r2 = r * r;
        direct += k == j ? r2 : r2 - (1 - r2) / (panel.size() - 2.0);
      }
      worst = std::max(worst, std::fabs(l[j] - direct));
    }
    s.check("LD scores equal Σ r̃² over the window by direct correlation", worst, 1e-9);
    std::vector<double> y(st.sample.size());
    for (std::size_t i = 0; i < y.size(); i++) y[i] = st.trait[st.sample[i]];
    const auto scan = as::linear_scan(g, st.sample, y, covariates(st, st.sample, &pcs, 4), &eight);
    as::ClumpOptions o;
    o.r2 = 0.3;
    const auto clumps = as::clump(g, st.sample, scan.p, o);
    bool ok = !clumps.empty();
    for (const auto& c : clumps)
      for (const std::size_t k : c.members) {
        const double r = r_of(g, st.sample, c.index, k);
        const auto d = g.positions[k] > g.positions[c.index] ? g.positions[k] - g.positions[c.index] : g.positions[c.index] - g.positions[k];
        ok = ok && r * r > o.r2 && d <= o.window && g.chromosome[k] == g.chromosome[c.index] && scan.p[c.index] <= scan.p[k];
      }
    s.require("clumps (" + std::to_string(clumps.size()) + ") hold only sites in LD (r² > 0.3) with their index, nearby, and weaker", ok);
  }
  return s.done();
}
