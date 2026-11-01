// Thermodynamic models of transcription: the probability that RNA polymerase occupies a promoter,
// in equilibrium with activators and repressors bound at sites near it, and the fold change in
// expression that regulators make.
//
// The model (Shea and Ackers 1985, Journal of Molecular Biology 181, 211–230, for λ's O_R; Bintu,
// Buchler, Garcia, Gerland, Hwa, Kondev and Phillips 2005, Current Opinion in Genetics & Development
// 15, 116–124, "Transcriptional regulation by the numbers: models", and 125–135, "… applications").
// A promoter region has sites; each site is bound by one species (RNA polymerase, an activator, a
// repressor) or empty. Each species X has n_X copies held on N_NS non-specific sites of the genome,
// and a specific site binds it with energy Δε (specific minus non-specific, in k_BT); a bound site
// then has the statistical weight x = (n_X/N_NS) e^{−Δε} (Bintu et al., eqs. 1–4 in the limit
// N_NS ≫ n_X). Sites that overlap cannot be bound together; two bound sites may interact with energy
// ε_ij (negative: they attract), a factor ω_ij = e^{−ε_ij}. The weight of a state (a set of bound
// sites, no two overlapping) is the product of its sites' weights and its pairs' factors; the
// partition function Z sums the states; the probability of a site's occupancy is the sum over the
// states that bind it, over Z. Here the states are enumerated exactly (2^k for k sites, k ≤ 20).
//
// Results the enumeration reproduces (the checks): with RNA polymerase alone, p_bound = 1/(1 +
// (N_NS/P) e^{Δε_pd}) (eq. 4); with regulators, p_bound = 1/(1 + (N_NS/(P F_reg)) e^{Δε_pd}) (eq. 5),
// the regulation factor F_reg = (1 + a f)/(1 + a) for an activator recruiting polymerase with f =
// e^{−ε_ap} (eq. 8), 1/(1 + r) for a repressor excluding it (eq. 9), a = (A/N_NS) e^{−Δε_ad} and r
// likewise; a helper raising an activator's effective number by (1 + h ω)/(1 + h) (eq. 10); and for
// two repressors that both exclude polymerase, 1/(1 + r₁ + r₂ + ω r₁ r₂), which is 1/((1 + r₁)(1 +
// r₂)) when they do not interact (Bintu et al. 2005, Table 1). In the weak-promoter limit (p_bound ≪
// 1) the fold change p_bound(regulated)/p_bound(unregulated) is F_reg (the companion paper's
// §"Quantitative Characteristics").
//
// Energies from sequence. The log-odds of a site under a position weight matrix approximates
// −Δε/k_BT up to a scale and an offset (Berg and von Hippel 1987, Journal of Molecular Biology 193,
// 723–750: the frequency of a base at a position in natural sites is about proportional to e^{−ε};
// Stormo 2000, Bioinformatics 16, 16–23): site_energy() takes −Δε = λ · score(bits) · ln 2 + μ with a
// caller's λ (k_BT per nat of log-odds; 1 in the ideal case) and μ.
//
// Limits: equilibrium is assumed (the open complex forms slowly beside binding and unbinding);
// expression is taken proportional to p_bound; DNA looping enters only as an interaction factor
// between the two sites a loop joins. Complexity: O(2^k k²) for k sites.

#pragma once

#include <genomics/annotation/motifs.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace genomics::regulation {

struct Species {
  std::string name;
  double copies = 0;            // n_X, molecules available for binding
  bool polymerase = false;      // RNA polymerase, whose occupancy makes expression
};

struct BindingSite {
  std::size_t species = 0;
  double energy = 0;            // Δε, specific minus non-specific, k_BT
  long begin = 0, end = 0;      // its footprint on the DNA, [begin, end): overlapping sites exclude each other
  std::string name;
};

struct SiteInteraction {
  std::size_t a = 0, b = 0;     // two sites
  double energy = 0;            // ε_ab, k_BT; negative attracts (f or ω = e^{−ε} > 1)
};

struct Architecture {
  std::vector<Species> species;
  std::vector<BindingSite> sites;
  std::vector<SiteInteraction> interactions;
  double nonspecific = 5e6;     // N_NS, about the genome's length in base pairs
};

struct Occupancy {
  double z = 0;                       // the partition function, the empty state's weight 1
  double p_bound = 0;                 // P(some polymerase site bound)
  std::vector<double> site;           // P(site bound)
  std::vector<std::uint32_t> states;  // the allowed states (bit masks of bound sites)
  std::vector<double> weights;        // their weights
};

inline Occupancy occupancy(const Architecture& a) {
  const std::size_t k = a.sites.size();
  if (k > 20) throw std::invalid_argument("genomics::regulation::occupancy: at most 20 sites");
  std::vector<double> w(k);
  for (std::size_t i = 0; i < k; i++) {
    const auto& s = a.sites[i];
    if (s.species >= a.species.size()) throw std::invalid_argument("genomics::regulation::occupancy: a site of an unknown species");
    w[i] = a.species[s.species].copies / a.nonspecific * std::exp(-s.energy);
  }
  std::vector<std::uint32_t> clash(k, 0);
  for (std::size_t i = 0; i < k; i++)
    for (std::size_t j = 0; j < k; j++)
      if (i != j && a.sites[i].begin < a.sites[j].end && a.sites[j].begin < a.sites[i].end) clash[i] |= 1u << j;
  std::vector<std::vector<double>> omega(k, std::vector<double>(k, 1.0));
  for (const auto& it : a.interactions) {
    if (it.a >= k || it.b >= k) throw std::invalid_argument("genomics::regulation::occupancy: an interaction of unknown sites");
    omega[it.a][it.b] = omega[it.b][it.a] = omega[it.a][it.b] * std::exp(-it.energy);
  }
  Occupancy o;
  o.site.assign(k, 0);
  double bound = 0;
  for (std::uint32_t m = 0; m < (1u << k); m++) {
    bool ok = true;
    double v = 1;
    for (std::size_t i = 0; i < k && ok; i++) {
      if (!(m >> i & 1u)) continue;
      if (clash[i] & m) { ok = false; break; }
      v *= w[i];
      for (std::size_t j = i + 1; j < k; j++)
        if (m >> j & 1u) v *= omega[i][j];
    }
    if (!ok) continue;
    o.states.push_back(m);
    o.weights.push_back(v);
    o.z += v;
    bool poly = false;
    for (std::size_t i = 0; i < k; i++)
      if (m >> i & 1u) {
        o.site[i] += v;
        poly = poly || a.species[a.sites[i].species].polymerase;
      }
    if (poly) bound += v;
  }
  for (auto& p : o.site) p /= o.z;
  o.p_bound = bound / o.z;
  return o;
}

// p_bound with the regulators over p_bound without them (their copies set to zero).
inline double fold_change(const Architecture& a) {
  Architecture bare = a;
  for (auto& s : bare.species)
    if (!s.polymerase) s.copies = 0;
  return occupancy(a).p_bound / occupancy(bare).p_bound;
}

// Bintu et al.'s closed forms, for the checks and for callers who want them.
inline double p_bound_basal(double P, double nonspecific, double energy) { return 1 / (1 + nonspecific / P * std::exp(energy)); }                // eq. 4
inline double p_bound_regulated(double P, double nonspecific, double energy, double F) { return 1 / (1 + nonspecific / (P * F) * std::exp(energy)); }  // eq. 5
inline double regulation_activation(double a, double f) { return (1 + a * f) / (1 + a); }                                                         // eq. 8
inline double regulation_repression(double r) { return 1 / (1 + r); }                                                                             // eq. 9

// −Δε/k_BT of a site from its PWM score in bits: λ ln 2 · score + μ.
inline double site_energy(double score_bits, double lambda = 1, double mu = 0) { return -(lambda * std::numbers::ln2 * score_bits + mu); }

}  // namespace genomics::regulation
