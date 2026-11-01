// Heritability from genotypes: the share of a trait's variance that the variants of a sample explain
// together (SNP heritability), by restricted maximum likelihood on the genetic relationship matrix and by
// Haseman–Elston regression.
//
// - The genetic relationship matrix (GRM; Yang, Lee, Goddard and Visscher 2011, AJHG 88, 76–82, eq. 3):
//   A_ik = (1/m) Σ_j (g_ij − 2p_j)(g_kj − 2p_j)/(2p_j(1 − p_j)). It is half population/structure.hpp's
//   kinship M Mᵀ/m (Patterson's standardisation divides by p(1 − p), not 2p(1 − p)), which is reused.
// - GREML: y = Cγ + u + e, u ~ N(0, σ_g² A), e ~ N(0, σ_e² I), h² = σ_g²/(σ_g² + σ_e²), by tyche's
//   reml_spectral (exact, one decomposition) or reml (AI-REML, as GCTA). Under the model of
//   cohort.hpp, with every causal variant among the sites, h² is the simulated heritability; the
//   estimator is unbiased when the effects are standardised (the α = −1 model the simulation uses).
// - Haseman–Elston: tyche's haseman_elston on the GRM, a method of moments.
// - A case–control trait's heritability on the observed 0/1 scale transforms to the liability scale by
//   Lee, Wray, Goddard and Visscher (2011, AJHG 88, 294–305, eq. 23): h²_l = h²_o K²(1 − K)²/(z²P(1 − P)),
//   K the prevalence, P the sample's fraction of cases and z the standard normal density at the
//   threshold Φ⁻¹(1 − K) (Robertson's 1950 transformation when P = K).
//
// COMPLEXITY. The GRM O(n² m) (population's kinship, on vector lanes and in parallel); REML O(n³).

#pragma once

#include <evolution/population/genotypes.hpp>
#include <evolution/population/structure.hpp>

#include <tyche/probability/distributions.hpp>
#include <tyche/regression/variance_components.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/rectangular.hpp>

#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

namespace evolution::association {

// The GRM of the samples at the sites listed (empty: every site), n × n.
inline calx::linalg::Matrix genetic_relationship(const population::Genotypes& g, std::span<const std::size_t> samples, std::span<const std::size_t> sites = {},
                                                 calx::scheduler::ThreadPool* pool = nullptr) {
  population::Genotypes sub = population::select_samples(g, samples);
  if (!sites.empty()) sub = population::select_sites(sub, sites);
  auto k = population::kinship(sub, pool).x;
  for (double& v : k.a) v /= 2;
  return k;
}

// h² on the liability scale from h² on the observed 0/1 scale (Lee et al. 2011, eq. 23).
inline double liability_scale(double h2_observed, double prevalence, double case_fraction) {
  if (!(prevalence > 0 && prevalence < 1) || !(case_fraction > 0 && case_fraction < 1))
    throw std::invalid_argument("evolution::association::liability_scale: a prevalence and a case fraction in (0, 1)");
  const double t = tyche::probability::Normal{}.upper_quantile(prevalence);
  const double z = std::exp(-0.5 * t * t) / std::sqrt(2 * std::numbers::pi);
  const double K = prevalence, P = case_fraction;
  return h2_observed * K * K * (1 - K) * (1 - K) / (z * z * P * (1 - P));
}

struct Heritability {
  tyche::regression::VarianceComponents greml;
  tyche::regression::HasemanElston he;
  std::size_t samples = 0, sites = 0;
};

// GREML (spectral) and Haseman–Elston on the GRM of the samples, the trait adjusted for C.
inline Heritability snp_heritability(const population::Genotypes& g, std::span<const std::size_t> samples, std::span<const double> y,
                                     const calx::linalg::Rectangular& C, calx::scheduler::ThreadPool* pool = nullptr) {
  if (y.size() != samples.size() || C.rows != samples.size()) throw std::invalid_argument("evolution::association::snp_heritability: one y and one row of C per sample");
  const auto A = genetic_relationship(g, samples, {}, pool);
  Heritability h;
  h.samples = samples.size();
  h.sites = g.sites;
  h.greml = tyche::regression::reml_spectral(C, y, A);
  h.he = tyche::regression::haseman_elston(A, y, C);
  return h;
}

}  // namespace evolution::association
