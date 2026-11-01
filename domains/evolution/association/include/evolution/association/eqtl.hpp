// Expression quantitative trait loci in cis: for each gene, the association of its expression with
// every variant within a window of its transcription start site, the gene-level significance of the
// strongest by permutations and their beta approximation, and the false discovery rate over genes.
//
// THE NOMINAL PASS (Matrix eQTL; Shabalin 2012, Bioinformatics 28, 1353–1358). Expression and
// dosages are residualised on the covariates C (genotype principal components, sex, factors of the
// expression itself) and each pair's partial correlation r is formed; t = r √(ν/(1 − r²)) on
// ν = n − rank C − 1 degrees of freedom is exactly the t of the regression of expression on the dosage
// and C (Frisch–Waugh–Lovell, as in scan.hpp), β̂ = r ‖ẽ‖/‖g̃‖ the effect per allele, and p is
// Student's two-sided tail.
//
// THE PERMUTATION PASS (FastQTL; Ongen, Buil, Brown, Dermitzakis and Delaneau 2016, Bioinformatics 32,
// 1479–1485). A gene is summarised by its strongest variant, and the null law of that minimum depends
// on the window's linkage, so the residualised expression is permuted B times (a calx stream per gene,
// Fisher–Yates) and the largest r² over the window recorded each time. The empirical p is
// (b + 1)/(B + 1) (tyche testing's permutation_p_value); the permutation minima, as p-values on the
// effective degrees of freedom (learnt as FastQTL and tensorQTL learn them, tyche's effective_df), are
// fitted by a beta law (tyche's fit_beta), whose distribution function at the observed minimum is the
// adjusted p_beta, resolved far below 1/B.
//
// THE FALSE DISCOVERY RATE over genes is Storey's q-values of p_beta (tyche's qvalues; GTEx
// Consortium 2017, Nature 550, 204–213, Methods), eGenes those with q ≤ fdr. Each eGene's nominal
// threshold is the beta law's quantile at the p_beta that lies at the FDR boundary (the mean of the
// largest p_beta with q ≤ fdr and the least above, as GTEx and tensorQTL take it): its variants with
// nominal p below it are its significant pairs.
//
// COMPLEXITY. For a gene with p variants in its window and n samples, O(n p) nominal and O(B n p) for
// the permutations; genes in parallel, each writing its own slot, so the same on any number of threads.

#pragma once

#include <evolution/association/scan.hpp>
#include <evolution/population/genotypes.hpp>

#include <tyche/probability/distributions.hpp>
#include <tyche/testing/multiple.hpp>
#include <tyche/testing/permutation.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/eigen.hpp>
#include <calx/random/stream.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace evolution::association {

inline constexpr std::uint32_t PERMUTATION = 0x5045524d;

struct GeneLocus {
  std::uint32_t chromosome = 0;  // an index into the genotypes' chromosomes
  std::uint64_t tss = 0;         // 0-based
};

struct CisOptions {
  std::uint64_t window = 1'000'000;  // ± of the TSS (FastQTL's and GTEx's 1 Mb)
  std::size_t permutations = 1000;
  bool learn_df = true;              // FastQTL's effective degrees of freedom
  std::uint64_t seed = 1;
  double fdr = 0.05;
};

struct CisGene {
  std::size_t gene = 0;
  std::vector<std::size_t> sites;           // in the window, increasing
  std::vector<double> r, beta, se, t, p;    // nominal, per site of the window
  std::size_t lead = 0;                     // the position in `sites` of the least p
  double p_nominal = std::numeric_limits<double>::quiet_NaN();  // the lead's
  double df = 0, df_effective = 0;
  double beta_a = 1, beta_b = 1;            // the beta law of the permutation minima
  double p_beta = 1, p_permutation = 1;
  double q = 1;
  double threshold = 0;                     // the nominal threshold of an eGene (0 otherwise)
  bool egene = false;
};

struct CisMapping {
  std::vector<CisGene> genes;
  double pi0 = 1;                // Storey's π̂₀ of the genes' p_beta
  double p_threshold = 0;        // the p_beta at the FDR boundary
  std::size_t egenes = 0;
  std::size_t covariates = 0;    // rank of C
};

// cis-eQTL mapping of the genes' expression (each one value per listed sample) at the genotypes' sites.
inline CisMapping cis_eqtl(const population::Genotypes& g, std::span<const std::size_t> samples, const std::vector<std::vector<double>>& expression,
                           std::span<const GeneLocus> loci, const Rectangular& C, const CisOptions& o = {}, calx::scheduler::ThreadPool* pool = nullptr) {
  const std::size_t n = samples.size(), G = expression.size();
  if (loci.size() != G || C.rows != n) throw std::invalid_argument("evolution::association::cis_eqtl: one locus per gene and one row of C per sample");
  const Residualizer M(C);
  CisMapping out;
  out.covariates = M.rank();
  if (n < M.rank() + 4) throw std::invalid_argument("evolution::association::cis_eqtl: more samples than covariates");
  const double df = static_cast<double>(n - M.rank() - 1);
  out.genes.resize(G);
  const tyche::probability::StudentT tdist{df};
  const auto body = [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t gi = begin; gi < end; gi++) {
      CisGene& cg = out.genes[gi];
      cg.gene = gi;
      cg.df = cg.df_effective = df;
      if (expression[gi].size() != n) throw std::invalid_argument("evolution::association::cis_eqtl: one expression value per sample");
      auto e = M.residuals(expression[gi]);
      double ee = 0;
      for (const double v : e) ee += v * v;
      const double enorm = std::sqrt(ee);
      if (!(enorm > 0)) continue;
      for (double& v : e) v /= enorm;
      const std::uint64_t from = loci[gi].tss > o.window ? loci[gi].tss - o.window : 0;
      const auto [lo, hi] = g.chromosome_sites(loci[gi].chromosome);
      std::vector<double> Z;  // the window's residualised, unit-norm dosages, site-major
      std::vector<double> gnorm;
      for (std::size_t j = lo; j < hi; j++) {
        if (g.positions[j] < from || g.positions[j] > loci[gi].tss + o.window) continue;
        auto d = sample_dosages(g, j, samples);
        M.apply(d);
        double gg = 0;
        for (const double v : d) gg += v * v;
        if (!(gg > 1e-12 * static_cast<double>(n))) continue;
        const double norm = std::sqrt(gg);
        for (double& v : d) v /= norm;
        cg.sites.push_back(j);
        gnorm.push_back(norm);
        Z.insert(Z.end(), d.begin(), d.end());
      }
      const std::size_t p = cg.sites.size();
      if (p == 0) continue;
      double best = -1;
      for (std::size_t k = 0; k < p; k++) {
        double r = 0;
        for (std::size_t i = 0; i < n; i++) r += Z[k * n + i] * e[i];
        const double t = r * std::sqrt(df / std::max(1e-300, 1 - r * r));
        cg.r.push_back(r);
        cg.t.push_back(t);
        cg.beta.push_back(r * enorm / gnorm[k]);
        cg.se.push_back(r * enorm / gnorm[k] / t);
        cg.p.push_back(std::min(1.0, 2 * tdist.survival(std::fabs(t))));
        if (r * r > best) {
          best = r * r;
          cg.lead = k;
        }
      }
      cg.p_nominal = cg.p[cg.lead];
      if (o.permutations == 0) continue;
      // The permutations: the largest r² over the window, each time.
      calx::random::Stream rng(o.seed, PERMUTATION, static_cast<std::uint32_t>(gi), 0);
      std::vector<double> perm = e, maxima(o.permutations);
      std::size_t extreme = 0;
      for (std::size_t b = 0; b < o.permutations; b++) {
        for (std::size_t i = n; i > 1; i--) std::swap(perm[i - 1], perm[rng.below(i)]);
        double top = 0;
        for (std::size_t k = 0; k < p; k++) {
          const double* z = &Z[k * n];
          double r = 0;
          for (std::size_t i = 0; i < n; i++) r += z[i] * perm[i];
          top = std::max(top, r * r);
        }
        maxima[b] = std::min(top, 1 - 1e-15);
        extreme += top >= best;
      }
      cg.p_permutation = tyche::testing::permutation_p_value(extreme, o.permutations);
      cg.df_effective = o.learn_df ? tyche::testing::effective_df(maxima, df) : df;
      std::vector<double> pv(o.permutations);
      for (std::size_t b = 0; b < o.permutations; b++) pv[b] = std::clamp(tyche::testing::correlation_p_value(maxima[b], cg.df_effective), 1e-300, 1 - 1e-16);
      try {
        const auto fit = tyche::testing::fit_beta(pv);
        cg.beta_a = fit.a;
        cg.beta_b = fit.b;
        cg.p_beta = fit.cdf(tyche::testing::correlation_p_value(std::min(best, 1 - 1e-15), cg.df_effective));
      } catch (const std::exception&) {
        cg.p_beta = cg.p_permutation;
      }
    }
  };
  if (pool) calx::scheduler::parallel_for(*pool, G, 1, body);
  else body(0, G, 0);
  // q-values over the genes with variants in their windows.
  std::vector<double> pb;
  std::vector<std::size_t> which;
  for (std::size_t gi = 0; gi < G; gi++)
    if (!out.genes[gi].sites.empty() && o.permutations > 0) {
      pb.push_back(std::clamp(out.genes[gi].p_beta, 0.0, 1.0));
      which.push_back(gi);
    }
  if (!pb.empty()) {
    tyche::testing::Pi0Options po;
    if (pb.size() < 50) po.method = tyche::testing::Pi0Method::fixed;  // the smoother needs many tests; λ = ½ for few
    tyche::testing::QValues q;
    try {
      q = tyche::testing::qvalues(pb, po);
    } catch (const std::exception&) {
      q = tyche::testing::qvalues(pb, po, 1.0);
    }
    out.pi0 = q.pi0.pi0;
    double lower = -1, upper = 2;
    for (std::size_t k = 0; k < which.size(); k++) {
      CisGene& cg = out.genes[which[k]];
      cg.q = q.q[k];
      cg.egene = cg.q <= o.fdr;
      out.egenes += cg.egene;
      if (cg.egene) lower = std::max(lower, cg.p_beta);
      else upper = std::min(upper, cg.p_beta);
    }
    out.p_threshold = lower < 0 ? 0 : upper > 1 ? lower : (lower + upper) / 2;
    for (const std::size_t gi : which) {
      CisGene& cg = out.genes[gi];
      if (cg.egene) cg.threshold = tyche::probability::Beta{cg.beta_a, cg.beta_b}.quantile(std::min(out.p_threshold, 1.0));
    }
  }
  return out;
}

// The leading principal components of expression (genes × samples, each gene standardised), as
// covariates for the hidden factors that move many genes at once (in place of PEER factors).
inline std::vector<std::vector<double>> expression_components(const std::vector<std::vector<double>>& expression, std::size_t k) {
  const std::size_t G = expression.size();
  if (G == 0) return {};
  const std::size_t n = expression[0].size();
  std::vector<std::vector<double>> z(G, std::vector<double>(n));
  for (std::size_t g = 0; g < G; g++) {
    double m = 0, ss = 0;
    for (const double v : expression[g]) m += v;
    m /= static_cast<double>(n);
    for (const double v : expression[g]) ss += (v - m) * (v - m);
    const double sd = std::sqrt(ss / static_cast<double>(n - 1));
    for (std::size_t i = 0; i < n; i++) z[g][i] = sd > 0 ? (expression[g][i] - m) / sd : 0;
  }
  calx::linalg::Matrix S(G);
  for (std::size_t a = 0; a < G; a++)
    for (std::size_t b = 0; b <= a; b++) {
      double s = 0;
      for (std::size_t i = 0; i < n; i++) s += z[a][i] * z[b][i];
      S(a, b) = S(b, a) = s;
    }
  const auto e = calx::linalg::symmetric_eigen(S);
  k = std::min(k, G);
  std::vector<std::vector<double>> out(k, std::vector<double>(n, 0.0));
  for (std::size_t c = 0; c < k; c++) {
    const std::size_t col = G - 1 - c;
    for (std::size_t g = 0; g < G; g++)
      for (std::size_t i = 0; i < n; i++) out[c][i] += e.vectors(g, col) * z[g][i];
  }
  return out;
}

}  // namespace evolution::association
