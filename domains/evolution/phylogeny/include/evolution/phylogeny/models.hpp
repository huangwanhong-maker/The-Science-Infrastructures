// Nucleotide substitution models: reversible continuous-time Markov chains on {A, C, G, T}, their
// transition probabilities P(t) = exp(Qt) with the first two derivatives in t, and the
// distribution of rates across sites (discrete Γ, invariant sites).
//
// A reversible model is a stationary distribution π and symmetric exchangeabilities s_ij, with
// rates Q_ij = s_ij π_j (i ≠ j) and Q_ii = −Σ_{j≠i} Q_ij (Tavaré 1986, "Some probabilistic and
// statistical problems in the analysis of DNA sequences", Lectures on Mathematics in the Life
// Sciences 17, 57–86; Yang 2014, Molecular Evolution, ch. 1). The named models are its special
// cases, s in the order AC, AG, AT, CG, CT, GT:
//
//   JC69   π uniform, s = 1                     Jukes and Cantor 1969
//   K80    π uniform, s_AG = s_CT = κ           Kimura 1980, J. Mol. Evol. 16, 111–120
//   F81    π free, s = 1                         Felsenstein 1981, J. Mol. Evol. 17, 368–376
//   HKY85  π free, s_AG = s_CT = κ              Hasegawa, Kishino and Yano 1985, J. Mol. Evol.
//                                                22, 160–174
//   TN93   π free, s_AG = κ₁, s_CT = κ₂         Tamura and Nei 1993, Mol. Biol. Evol. 10, 512–526
//   GTR    π free, s free                        Tavaré 1986
//
// Q is normalised so that the expected rate at stationarity, −Σ_i π_i Q_ii, is 1: branch
// lengths are expected substitutions per site.
//
// Reversibility, π_i Q_ij = π_j Q_ji, makes S = Π^{1/2} Q Π^{−1/2} symmetric, S_ij = s_ij
// √(π_i π_j); with S = U Λ Uᵀ (calx::linalg::symmetric_eigen),
//
//   P(t) = Π^{−1/2} U e^{Λt} Uᵀ Π^{1/2},  P_ij(t) = Σ_k A_ik e^{λ_k t} B_kj,
//
// A = Π^{−1/2} U, B = Uᵀ Π^{1/2}, computed once; then P(t), dP/dt = A Λ e^{Λt} B and
// d²P/dt² = A Λ² e^{Λt} B cost four exponentials and 48 products each, for every t (Yang 2014,
// ch. 1). The eigenvalue 0 and its vector √π are set exactly, so that P(t) → Π (every row π)
// as t → ∞ to rounding. Entries of P below zero by rounding (for t near 0) are set to zero.
//
// Rates across sites (Yang 1994, "Maximum likelihood phylogenetic estimation from DNA
// sequences with variable rates over sites: approximate methods", J. Mol. Evol. 39, 306–314):
// a site's rate is drawn from Γ(α, α), mean 1, approximated by K equiprobable categories whose
// rates are the means of the distribution in each (the cut points b_k the quantiles k/K, and
// r_k = K [P(α + 1, α b_k) − P(α + 1, α b_{k−1})], P the regularised incomplete gamma
// function; Yang 1994), or their medians rescaled to mean 1. A proportion p of
// invariant sites (+I) is one more category, of rate 0 and weight p, the others' rates divided
// by 1 − p so that the mean stays 1 (Gu, Fu and Li 1995, Mol. Biol. Evol. 12, 546–557). The
// quantiles are calx::special_functions::gamma_p_inverse's, to a relative 1e-13.

#pragma once

#include <calx/linalg/eigen.hpp>
#include <calx/special_functions/gamma.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace evolution::phylogeny {

enum class ModelKind { jc69, k80, f81, hky85, tn93, gtr };

// Exchangeabilities in the order AC, AG, AT, CG, CT, GT.
inline constexpr int AC = 0, AG = 1, AT = 2, CG = 3, CT = 4, GT = 5;

struct Model {
  ModelKind kind = ModelKind::jc69;
  std::array<double, 4> pi{0.25, 0.25, 0.25, 0.25};
  std::array<double, 6> exchange{1, 1, 1, 1, 1, 1};  // as given; Q is their normalisation
  // Derived by prepare():
  std::array<double, 16> Q{};
  std::array<double, 4> lambda{};  // eigenvalues of Q, ascending, the last 0
  std::array<double, 16> left{};   // A_ik
  std::array<double, 16> right{};  // B_kj
  double scale = 1;                // the factor that normalised the rates

  // The exchangeability of the pair (i, j), i ≠ j.
  static int pair(int i, int j) {
    if (i > j) std::swap(i, j);
    return i == 0 ? j - 1 : i == 1 ? j + 1 : 5;
  }

  void prepare() {
    double total = 0;
    for (double p : pi) {
      if (!(p > 0)) throw std::invalid_argument("evolution::phylogeny::Model: every base has a frequency above zero");
      total += p;
    }
    for (double& p : pi) p /= total;
    for (double s : exchange)
      if (!(s > 0) || !std::isfinite(s)) throw std::invalid_argument("evolution::phylogeny::Model: exchangeabilities are positive and finite");
    Q.fill(0);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        if (i != j) {
          Q[i * 4 + j] = exchange[static_cast<std::size_t>(pair(i, j))] * pi[static_cast<std::size_t>(j)];
          Q[i * 4 + i] -= Q[i * 4 + j];
        }
    double rate = 0;
    for (int i = 0; i < 4; i++) rate -= pi[static_cast<std::size_t>(i)] * Q[i * 4 + i];
    scale = 1 / rate;
    for (double& q : Q) q *= scale;
    calx::linalg::Matrix S(4);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) S(i, j) = i == j ? Q[i * 4 + i] : Q[i * 4 + j] * std::sqrt(pi[static_cast<std::size_t>(i)] / pi[static_cast<std::size_t>(j)]);
    auto e = calx::linalg::symmetric_eigen(S);
    // The largest eigenvalue is 0, with the vector √π: set them exactly.
    e.values[3] = 0;
    for (int i = 0; i < 4; i++) e.vectors(i, 3) = std::sqrt(pi[static_cast<std::size_t>(i)]);
    for (int k = 0; k < 4; k++) lambda[static_cast<std::size_t>(k)] = e.values[static_cast<std::size_t>(k)];
    for (int i = 0; i < 4; i++)
      for (int k = 0; k < 4; k++) {
        left[static_cast<std::size_t>(i * 4 + k)] = e.vectors(i, k) / std::sqrt(pi[static_cast<std::size_t>(i)]);
        right[static_cast<std::size_t>(k * 4 + i)] = e.vectors(i, k) * std::sqrt(pi[static_cast<std::size_t>(i)]);
      }
  }

  // P(t), row-major: P[4i + j] = P(j at time t | i at time 0).
  void transition(double t, double* P) const {
    double e[4];
    for (int k = 0; k < 4; k++) e[k] = std::exp(lambda[static_cast<std::size_t>(k)] * t);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        double s = 0;
        for (int k = 0; k < 4; k++) s += left[static_cast<std::size_t>(i * 4 + k)] * e[k] * right[static_cast<std::size_t>(k * 4 + j)];
        P[i * 4 + j] = s > 0 ? s : 0;
      }
  }
  // P(t) and its first two derivatives in t.
  void transition(double t, double* P, double* dP, double* d2P) const {
    double e[4], l[4];
    for (int k = 0; k < 4; k++) {
      l[k] = lambda[static_cast<std::size_t>(k)];
      e[k] = std::exp(l[k] * t);
    }
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        double s = 0, s1 = 0, s2 = 0;
        for (int k = 0; k < 4; k++) {
          const double c = left[static_cast<std::size_t>(i * 4 + k)] * e[k] * right[static_cast<std::size_t>(k * 4 + j)];
          s += c;
          s1 += c * l[k];
          s2 += c * l[k] * l[k];
        }
        P[i * 4 + j] = s > 0 ? s : 0;
        dP[i * 4 + j] = s1;
        d2P[i * 4 + j] = s2;
      }
  }
};

inline Model make_model(ModelKind kind, const std::array<double, 4>& pi, const std::array<double, 6>& exchange) {
  Model m;
  m.kind = kind;
  m.pi = pi;
  m.exchange = exchange;
  m.prepare();
  return m;
}
inline Model jc69() { return make_model(ModelKind::jc69, {0.25, 0.25, 0.25, 0.25}, {1, 1, 1, 1, 1, 1}); }
inline Model k80(double kappa) { return make_model(ModelKind::k80, {0.25, 0.25, 0.25, 0.25}, {1, kappa, 1, 1, kappa, 1}); }
inline Model f81(const std::array<double, 4>& pi) { return make_model(ModelKind::f81, pi, {1, 1, 1, 1, 1, 1}); }
inline Model hky85(double kappa, const std::array<double, 4>& pi) { return make_model(ModelKind::hky85, pi, {1, kappa, 1, 1, kappa, 1}); }
inline Model tn93(double kappa_purines, double kappa_pyrimidines, const std::array<double, 4>& pi) {
  return make_model(ModelKind::tn93, pi, {1, kappa_purines, 1, 1, kappa_pyrimidines, 1});
}
inline Model gtr(const std::array<double, 6>& exchange, const std::array<double, 4>& pi) { return make_model(ModelKind::gtr, pi, exchange); }

// Rates across sites: categories with their rates and probabilities (weights sum to 1).
struct Rates {
  std::vector<double> rate{1}, weight{1};
  std::size_t size() const { return rate.size(); }
};

// The K rates of Yang's discrete Γ with shape α: the category means, or the medians rescaled
// to mean 1.
inline std::vector<double> discrete_gamma(double alpha, int K, bool median = false) {
  if (!(alpha > 0) || K < 1) throw std::invalid_argument("evolution::phylogeny::discrete_gamma: alpha > 0 and K >= 1");
  std::vector<double> r(static_cast<std::size_t>(K));
  if (K == 1) {
    r[0] = 1;
    return r;
  }
  namespace sf = calx::special_functions;
  if (median) {
    double total = 0;
    for (int k = 0; k < K; k++) {
      r[static_cast<std::size_t>(k)] = sf::gamma_p_inverse(alpha, (2.0 * k + 1) / (2.0 * K)) / alpha;
      total += r[static_cast<std::size_t>(k)];
    }
    for (double& x : r) x *= K / total;
    return r;
  }
  double previous = 0;  // P(α + 1, α b_{k−1})
  for (int k = 0; k < K; k++) {
    const double next = k + 1 == K ? 1.0 : sf::gamma_p(alpha + 1, sf::gamma_p_inverse(alpha, static_cast<double>(k + 1) / K));
    r[static_cast<std::size_t>(k)] = K * (next - previous);
    previous = next;
  }
  return r;
}

// Discrete Γ with K categories (K = 1: equal rates) and a proportion of invariant sites.
inline Rates gamma_rates(double alpha, int K, double invariant = 0, bool median = false) {
  if (invariant < 0 || invariant >= 1) throw std::invalid_argument("evolution::phylogeny::gamma_rates: 0 <= invariant < 1");
  Rates r;
  r.rate = discrete_gamma(alpha, K, median);
  r.weight.assign(r.rate.size(), (1 - invariant) / static_cast<double>(K));
  if (invariant > 0) {
    for (double& x : r.rate) x /= 1 - invariant;
    r.rate.push_back(0);
    r.weight.push_back(invariant);
  }
  return r;
}

// How rates across sites are modelled, as parameters to estimate.
struct RateSpec {
  int categories = 1;      // 1: equal rates; K > 1: discrete Γ with K categories
  double alpha = 1;        // the Γ shape
  double invariant = 0;    // the proportion of invariant sites
  bool estimate_invariant = false;
  bool median = false;
  Rates rates() const { return gamma_rates(alpha, categories, invariant, median); }
};

}  // namespace evolution::phylogeny
