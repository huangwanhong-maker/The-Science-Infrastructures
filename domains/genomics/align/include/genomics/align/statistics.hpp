// Karlin–Altschul statistics of local alignment: λ, K and H of a scoring under background
// frequencies, bit scores and E-values, and NCBI BLAST's table of the gapped parameters.
//
// Ungapped theory (Karlin and Altschul 1990, PNAS 87, 2264–2268; Karlin and Dembo 1992, Adv. Appl.
// Prob. 24, 113–140). Let the score of one aligned pair take the value s with probability
// P(s) = Σ_{a,b} p_a q_b [s_ab = s], with negative mean and some positive score. Then:
//
//   λ   the unique positive root of Σ_s P(s) e^{λs} = 1 (the function is convex, 0 at λ = 0 with
//       slope E[s] < 0, and grows without bound), found by Brent's method (calx::solvers);
//   H   the relative entropy of the target frequencies p_a q_b e^{λ s_ab} to the background,
//       H = λ Σ_s s P(s) e^{λs}, in nats per aligned pair;
//   K   K = λ δ e^{−2σ} / (H (1 − e^{−λδ})), δ the greatest common divisor of the scores of positive
//       probability and σ = Σ_{k≥1} (1/k) [E(e^{λ S_k}; S_k < 0) + P(S_k ≥ 0)] for S_k the sum of k
//       independent scores (Karlin and Altschul 1990, the Appendix; the series as NCBI BLAST sums it
//       in BlastKarlinLHtoK of blast_stat.c, the NCBI C++ Toolkit, retrieved 2026-10-24, each P_k by
//       convolution, here to a term below 10⁻¹² rather than BLAST's 10⁻⁴). When the reduced scores
//       reach only −1 or only +1 the series has closed forms, which BLAST uses and so does this
//       header: both −1 and +1, K = (P(−1) − P(1))² / P(−1); low −1, K = (μ² λ / H)(1 − e^{−λ});
//       high +1, K = (H / λ)(1 − e^{−λ}), in reduced units (μ the mean score).
//
// The number of distinct local alignments of two random sequences of lengths m and n with score at
// least S is then asymptotically Poisson with mean E = K m n e^{−λS}; the bit score is
// S′ = (λS − ln K) / ln 2, so that E = m n 2^{−S′}. The check holds BLOSUM62 under Robinson and
// Robinson's frequencies (scoring.hpp) to λ = 0.3176, K = 0.134, H = 0.4012, NCBI BLAST's ungapped
// values, and the seven other matrices (BLOSUM80 in the half-bit form BLAST uses) to theirs.
//
// Gapped alignment has no closed theory; its λ and K are estimated by simulation. The table below
// is NCBI BLAST's (blast_stat.c, the arrays <matrix>_values; Altschul and Gish 1996, Methods
// Enzymol. 266, 460–480, and Altschul, Bundschuh, Olsen and Hwa 2001, Nucleic Acids Res. 29,
// 351–361, for α and β, the slope and intercept of the mean length of an alignment against its
// score; C, α_v and σ for Spouge's finite-size correction), transcribed by script, every row: gap
// open and extend (32767 for no gaps, the ungapped row), λ, K, H, α, β, C, α_v, σ. The rows are
// for BLAST's built-in matrices, which equal NCBI's FTP files (scoring.hpp) on the twenty
// residues except BLOSUM80: BLAST's is at half-bit scale, the file at third-bit, so BLAST's
// BLOSUM80 rows are filed under "BLOSUM80-BLAST", the matrix they belong to, and gapped_table
// ("BLOSUM80") is empty rather than wrong.
//
// K is computed to a term of 10⁻¹² by default; karlin_K_blast truncates as BLAST's run-time code
// does (a term of at most 10⁻⁴, at most 100 terms), which omits positive terms of σ and so leaves K
// larger by up to a few parts in 10³. BLAST's tabulated ungapped K lie between the two and agree
// with the full series to within one unit of their last digit (the check).
//
// Accuracy: λ to 10⁻¹² relative; K to the truncation of its series (relative 10⁻⁹ or better); the
// theory itself is asymptotic in m and n, without BLAST's edge corrections (effective lengths),
// which matter when ln(K m n) / H is not small against m and n.

#pragma once

#include <genomics/align/types.hpp>

#include <calx/solvers/roots.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <numeric>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::align {

// The distribution of the score of one aligned pair: p[s − low] = P(s) for s in [low, high].
struct ScoreDistribution {
  int low = 0, high = 0;
  std::vector<double> p;
  double probability(int s) const { return s < low || s > high ? 0.0 : p[static_cast<std::size_t>(s - low)]; }
  double mean() const {
    double m = 0;
    for (int s = low; s <= high; s++) m += s * probability(s);
    return m;
  }
};

// P(s) = Σ p_a q_b [s_ab = s] over the codes of the matrix with positive frequency in both; the
// frequencies are normalised, so counts or percentages will do.
inline ScoreDistribution score_distribution(const Matrix& matrix, std::span<const double> p_query, std::span<const double> p_target) {
  if (p_query.size() < static_cast<std::size_t>(matrix.size) || p_target.size() < static_cast<std::size_t>(matrix.size))
    throw std::invalid_argument("genomics::align::score_distribution: a frequency for every code of the matrix");
  ScoreDistribution d;
  bool any = false;
  for (int a = 0; a < matrix.size; a++)
    for (int b = 0; b < matrix.size; b++)
      if (p_query[static_cast<std::size_t>(a)] > 0 && p_target[static_cast<std::size_t>(b)] > 0) {
        const int s = matrix(static_cast<Code>(a), static_cast<Code>(b));
        if (!any) { d.low = d.high = s; any = true; }
        d.low = std::min(d.low, s);
        d.high = std::max(d.high, s);
      }
  if (!any) throw std::invalid_argument("genomics::align::score_distribution: no code has positive frequency");
  d.p.assign(static_cast<std::size_t>(d.high - d.low + 1), 0.0);
  double total = 0;
  for (int a = 0; a < matrix.size; a++)
    for (int b = 0; b < matrix.size; b++) {
      const double w = p_query[static_cast<std::size_t>(a)] * p_target[static_cast<std::size_t>(b)];
      if (w > 0) {
        d.p[static_cast<std::size_t>(matrix(static_cast<Code>(a), static_cast<Code>(b)) - d.low)] += w;
        total += w;
      }
    }
  for (auto& x : d.p) x /= total;
  return d;
}

// A distribution given directly: scores low .. low + p.size() − 1 with weights p (normalised).
inline ScoreDistribution score_distribution(int low, std::vector<double> p) {
  ScoreDistribution d;
  d.low = low;
  d.high = low + static_cast<int>(p.size()) - 1;
  double total = 0;
  for (double x : p) total += x;
  for (auto& x : p) x /= total;
  d.p = std::move(p);
  return d;
}

struct KarlinAltschul {
  double lambda = 0, K = 0, H = 0;  // λ per score unit, K, H in nats per aligned pair
};

namespace detail {

inline void require_local_regime(const ScoreDistribution& d, const char* who) {
  if (!(d.mean() < 0)) throw std::invalid_argument(std::string(who) + ": the expected score must be negative");
  if (d.high <= 0 || d.probability(d.high) <= 0) throw std::invalid_argument(std::string(who) + ": some positive score must have positive probability");
}

inline int score_divisor(const ScoreDistribution& d) {
  int delta = 0;
  for (int s = d.low; s <= d.high; s++)
    if (d.probability(s) > 0) delta = std::gcd(delta, s < 0 ? -s : s);
  return delta;
}

}  // namespace detail

// λ: the positive root of Σ P(s) e^{λs} = 1.
inline double karlin_lambda(const ScoreDistribution& d) {
  detail::require_local_regime(d, "genomics::align::karlin_lambda");
  const auto f = [&](double lambda) {
    double v = 0;
    for (int s = d.low; s <= d.high; s++) v += d.probability(s) * std::exp(lambda * s);
    return v - 1;
  };
  double hi = 1;
  while (f(hi) <= 0) hi *= 2;
  double lo = hi;
  while (f(lo) >= 0) lo /= 2;
  return calx::solvers::brent(f, lo, hi, 1e-15).x;
}

// H = λ Σ s P(s) e^{λs}, in nats.
inline double karlin_H(const ScoreDistribution& d, double lambda) {
  double h = 0;
  for (int s = d.low; s <= d.high; s++) h += s * d.probability(s) * std::exp(lambda * s);
  return lambda * h;
}

// K by the series of Karlin and Altschul, always (without the closed forms): for checking them.
inline double karlin_K_series(const ScoreDistribution& d, double lambda, double H, double tolerance = 1e-12, int max_terms = 4000) {
  detail::require_local_regime(d, "genomics::align::karlin_K");
  const int delta = detail::score_divisor(d);
  const int low = d.low / delta, high = d.high / delta, range = high - low;
  const double l = lambda * delta;
  std::vector<double> step(static_cast<std::size_t>(range + 1), 0.0);
  for (int s = low; s <= high; s++) step[static_cast<std::size_t>(s - low)] = d.probability(s * delta);
  // P_k on [k low, k high], by convolution with the one-step distribution.
  std::vector<double> pk{1.0}, next;
  double sigma = 0;
  for (int k = 1; k <= max_terms; k++) {
    next.assign(pk.size() + static_cast<std::size_t>(range), 0.0);
    for (std::size_t a = 0; a < pk.size(); a++)
      if (pk[a] != 0)
        for (std::size_t b = 0; b < step.size(); b++) next[a + b] += pk[a] * step[b];
    pk.swap(next);
    const int base = k * low;  // the score of pk[0]
    double term = 0;
    for (std::size_t a = 0; a < pk.size(); a++) {
      const int s = base + static_cast<int>(a);
      term += s < 0 ? pk[a] * std::exp(l * s) : pk[a];
    }
    sigma += term / k;
    if (term / k <= tolerance) break;  // BLAST stops once a term is at most its sumlimit
  }
  return l * std::exp(-2 * sigma) / (H * -std::expm1(-l));
}

// K, by the closed forms where they apply and the series otherwise (to `tolerance`).
inline double karlin_K(const ScoreDistribution& d, double lambda, double H, double tolerance = 1e-12, int max_terms = 4000) {
  detail::require_local_regime(d, "genomics::align::karlin_K");
  const int delta = detail::score_divisor(d);
  const int low = d.low / delta, high = d.high / delta;
  const double l = lambda * delta;
  if (low == -1 && high == 1) {
    const double pm = d.probability(-delta), pp = d.probability(delta);
    return (pm - pp) * (pm - pp) / pm;
  }
  if (low == -1 || high == 1) {
    double first = H / l;
    if (high != 1) {
      const double mu = d.mean() / delta;
      first = mu * mu / first;
    }
    return first * -std::expm1(-l);
  }
  return karlin_K_series(d, lambda, H, tolerance, max_terms);
}

// K as NCBI BLAST computes it: BlastKarlinLHtoK's sumlimit 10⁻⁴ and at most 100 terms.
inline double karlin_K_blast(const ScoreDistribution& d, double lambda, double H) { return karlin_K(d, lambda, H, 1e-4, 100); }

inline KarlinAltschul karlin_altschul(const ScoreDistribution& d) {
  KarlinAltschul k;
  k.lambda = karlin_lambda(d);
  k.H = karlin_H(d, k.lambda);
  k.K = karlin_K(d, k.lambda, k.H);
  return k;
}
inline KarlinAltschul karlin_altschul(const Matrix& matrix, std::span<const double> p_query, std::span<const double> p_target) {
  return karlin_altschul(score_distribution(matrix, p_query, p_target));
}

// The bit score (λS − ln K) / ln 2, the E-value K m n e^{−λS}, and the probability 1 − e^{−E} of
// at least one alignment scoring S or more.
inline double bit_score(double raw, const KarlinAltschul& k) { return (k.lambda * raw - std::log(k.K)) / std::numbers::ln2; }
inline double evalue(double raw, double m, double n, const KarlinAltschul& k) { return k.K * m * n * std::exp(-k.lambda * raw); }
inline double pvalue(double e) { return -std::expm1(-e); }

// ---- NCBI BLAST's gapped parameters ------------------------------------------------------------------
struct GappedParameters {
  std::string_view matrix;
  int gap_open, gap_extend;  // 32767 for both: the ungapped row
  double lambda, K, H, alpha, beta, C, alpha_v, sigma;
  bool ungapped() const { return gap_open == 32767; }
  KarlinAltschul statistics() const { return {lambda, K, H}; }
};

namespace detail {

// blast_stat.c (NCBI C++ Toolkit, src/algo/blast/core, retrieved 2026-10-24), the arrays
// <matrix>_values, in their order; the third column (decline to align, ignored by BLAST) dropped.
inline constexpr GappedParameters GAPPED[] = {
    // BLOSUM45 (blosum45_values)
    {"BLOSUM45", 32767, 32767, 0.2291, 0.0924, 0.2514, 0.9113, -5.7, 0.641318, 9.611060, 9.611060},
    {"BLOSUM45", 13, 3, 0.207, 0.049, 0.14, 1.5, -22, 0.671128, 35.855900, 35.963900},
    {"BLOSUM45", 12, 3, 0.199, 0.039, 0.11, 1.8, -34, 0.691530, 45.693600, 45.851700},
    {"BLOSUM45", 11, 3, 0.190, 0.031, 0.095, 2.0, -38, 0.691181, 62.874100, 63.103700},
    {"BLOSUM45", 10, 3, 0.179, 0.023, 0.075, 2.4, -51, 0.710529, 88.286800, 88.639100},
    {"BLOSUM45", 16, 2, 0.210, 0.051, 0.14, 1.5, -24, 0.666680, 36.279800, 36.452400},
    {"BLOSUM45", 15, 2, 0.203, 0.041, 0.12, 1.7, -31, 0.673871, 44.825700, 45.060400},
    {"BLOSUM45", 14, 2, 0.195, 0.032, 0.10, 1.9, -36, 0.685753, 60.736200, 61.102300},
    {"BLOSUM45", 13, 2, 0.185, 0.024, 0.084, 2.2, -45, 0.698480, 85.148100, 85.689400},
    {"BLOSUM45", 12, 2, 0.171, 0.016, 0.061, 2.8, -65, 0.713429, 127.758000, 128.582000},
    {"BLOSUM45", 19, 1, 0.205, 0.040, 0.11, 1.9, -43, 0.672302, 53.071400, 53.828200},
    {"BLOSUM45", 18, 1, 0.198, 0.032, 0.10, 2.0, -43, 0.682580, 72.342400, 73.403900},
    {"BLOSUM45", 17, 1, 0.189, 0.024, 0.079, 2.4, -57, 0.695035, 103.055000, 104.721000},
    {"BLOSUM45", 16, 1, 0.176, 0.016, 0.063, 2.8, -67, 0.712966, 170.100000, 173.003000},
    // BLOSUM50 (blosum50_values)
    {"BLOSUM50", 32767, 32767, 0.2318, 0.112, 0.3362, 0.6895, -4.0, 0.609639, 5.388310, 5.388310},
    {"BLOSUM50", 13, 3, 0.212, 0.063, 0.19, 1.1, -16, 0.639287, 18.113800, 18.202800},
    {"BLOSUM50", 12, 3, 0.206, 0.055, 0.17, 1.2, -18, 0.644715, 22.654600, 22.777700},
    {"BLOSUM50", 11, 3, 0.197, 0.042, 0.14, 1.4, -25, 0.656327, 29.861100, 30.045700},
    {"BLOSUM50", 10, 3, 0.186, 0.031, 0.11, 1.7, -34, 0.671150, 42.393800, 42.674000},
    {"BLOSUM50", 9, 3, 0.172, 0.022, 0.082, 2.1, -48, 0.694326, 66.069600, 66.516400},
    {"BLOSUM50", 16, 2, 0.215, 0.066, 0.20, 1.05, -15, 0.633899, 17.951800, 18.092100},
    {"BLOSUM50", 15, 2, 0.210, 0.058, 0.17, 1.2, -20, 0.641985, 21.940100, 22.141800},
    {"BLOSUM50", 14, 2, 0.202, 0.045, 0.14, 1.4, -27, 0.650682, 28.681200, 28.961900},
    {"BLOSUM50", 13, 2, 0.193, 0.035, 0.12, 1.6, -32, 0.660984, 42.059500, 42.471600},
    {"BLOSUM50", 12, 2, 0.181, 0.025, 0.095, 1.9, -41, 0.678090, 63.747600, 64.397300},
    {"BLOSUM50", 19, 1, 0.212, 0.057, 0.18, 1.2, -21, 0.635714, 26.311200, 26.923300},
    {"BLOSUM50", 18, 1, 0.207, 0.050, 0.15, 1.4, -28, 0.643523, 34.903700, 35.734800},
    {"BLOSUM50", 17, 1, 0.198, 0.037, 0.12, 1.6, -33, 0.654504, 48.895800, 50.148600},
    {"BLOSUM50", 16, 1, 0.186, 0.025, 0.10, 1.9, -42, 0.667750, 76.469100, 78.443000},
    {"BLOSUM50", 15, 1, 0.171, 0.015, 0.063, 2.7, -76, 0.694575, 140.053000, 144.160000},
    // BLOSUM62 (blosum62_values)
    {"BLOSUM62", 32767, 32767, 0.3176, 0.134, 0.4012, 0.7916, -3.2, 0.623757, 4.964660, 4.964660},
    {"BLOSUM62", 11, 2, 0.297, 0.082, 0.27, 1.1, -10, 0.641766, 12.673800, 12.757600},
    {"BLOSUM62", 10, 2, 0.291, 0.075, 0.23, 1.3, -15, 0.649362, 16.474000, 16.602600},
    {"BLOSUM62", 9, 2, 0.279, 0.058, 0.19, 1.5, -19, 0.659245, 22.751900, 22.950000},
    {"BLOSUM62", 8, 2, 0.264, 0.045, 0.15, 1.8, -26, 0.672692, 35.483800, 35.821300},
    {"BLOSUM62", 7, 2, 0.239, 0.027, 0.10, 2.5, -46, 0.702056, 61.238300, 61.886000},
    {"BLOSUM62", 6, 2, 0.201, 0.012, 0.061, 3.3, -58, 0.740802, 140.417000, 141.882000},
    {"BLOSUM62", 13, 1, 0.292, 0.071, 0.23, 1.2, -11, 0.647715, 19.506300, 19.893100},
    {"BLOSUM62", 12, 1, 0.283, 0.059, 0.19, 1.5, -19, 0.656391, 27.856200, 28.469900},
    {"BLOSUM62", 11, 1, 0.267, 0.041, 0.14, 1.9, -30, 0.669720, 42.602800, 43.636200},
    {"BLOSUM62", 10, 1, 0.243, 0.024, 0.10, 2.5, -44, 0.693267, 83.178700, 85.065600},
    {"BLOSUM62", 9, 1, 0.206, 0.010, 0.052, 4.0, -87, 0.731887, 210.333000, 214.842000},
    // BLOSUM80 (blosum80_values), for BLAST's half-bit BLOSUM80 (scoring.hpp: "BLOSUM80-BLAST")
    {"BLOSUM80-BLAST", 32767, 32767, 0.3430, 0.177, 0.6568, 0.5222, -1.6, 0.564057, 1.918130, 1.918130},
    {"BLOSUM80-BLAST", 25, 2, 0.342, 0.17, 0.66, 0.52, -1.6, 0.563956, 1.731000, 1.731300},
    {"BLOSUM80-BLAST", 13, 2, 0.336, 0.15, 0.57, 0.59, -3, 0.570979, 2.673470, 2.692300},
    {"BLOSUM80-BLAST", 9, 2, 0.319, 0.11, 0.42, 0.76, -6, 0.587837, 5.576090, 5.667860},
    {"BLOSUM80-BLAST", 8, 2, 0.308, 0.090, 0.35, 0.89, -9, 0.597556, 7.536950, 7.686230},
    {"BLOSUM80-BLAST", 7, 2, 0.293, 0.070, 0.27, 1.1, -14, 0.615254, 11.586600, 11.840400},
    {"BLOSUM80-BLAST", 6, 2, 0.268, 0.045, 0.19, 1.4, -19, 0.644054, 19.958100, 20.441200},
    {"BLOSUM80-BLAST", 11, 1, 0.314, 0.095, 0.35, 0.90, -9, 0.590702, 8.808610, 9.223320},
    {"BLOSUM80-BLAST", 10, 1, 0.299, 0.071, 0.27, 1.1, -14, 0.609620, 13.833800, 14.533400},
    {"BLOSUM80-BLAST", 9, 1, 0.279, 0.048, 0.20, 1.4, -19, 0.623800, 24.252000, 25.490400},
    // BLOSUM90 (blosum90_values)
    {"BLOSUM90", 32767, 32767, 0.3346, 0.190, 0.7547, 0.4434, -1.4, 0.544178, 1.377760, 1.377760},
    {"BLOSUM90", 9, 2, 0.310, 0.12, 0.46, 0.67, -6, 0.570267, 4.232290, 4.334170},
    {"BLOSUM90", 8, 2, 0.300, 0.099, 0.39, 0.76, -7, 0.581580, 5.797020, 5.961420},
    {"BLOSUM90", 7, 2, 0.283, 0.072, 0.30, 0.93, -11, 0.600024, 9.040880, 9.321600},
    {"BLOSUM90", 6, 2, 0.259, 0.048, 0.22, 1.2, -16, 0.629344, 16.024400, 16.531600},
    {"BLOSUM90", 11, 1, 0.302, 0.093, 0.39, 0.78, -8, 0.576919, 7.143250, 7.619190},
    {"BLOSUM90", 10, 1, 0.290, 0.075, 0.28, 1.04, -15, 0.591366, 11.483900, 12.269800},
    {"BLOSUM90", 9, 1, 0.265, 0.044, 0.20, 1.3, -19, 0.613013, 21.408300, 22.840900},
    // PAM30 (pam30_values)
    {"PAM30", 32767, 32767, 0.3400, 0.283, 1.754, 0.1938, -0.3, 0.436164, 0.161818, 0.161818},
    {"PAM30", 7, 2, 0.305, 0.15, 0.87, 0.35, -3, 0.479087, 1.014010, 1.162730},
    {"PAM30", 6, 2, 0.287, 0.11, 0.68, 0.42, -4, 0.499980, 1.688060, 1.951430},
    {"PAM30", 5, 2, 0.264, 0.079, 0.45, 0.59, -7, 0.533009, 3.377010, 3.871950},
    {"PAM30", 10, 1, 0.309, 0.15, 0.88, 0.35, -3, 0.474741, 1.372050, 1.788770},
    {"PAM30", 9, 1, 0.294, 0.11, 0.61, 0.48, -6, 0.492716, 2.463920, 3.186150},
    {"PAM30", 8, 1, 0.270, 0.072, 0.40, 0.68, -10, 0.521286, 5.368130, 6.763480},
    {"PAM30", 15, 3, 0.339, 0.28, 1.70, 0.20, -0.5, 0.437688, 0.157089, 0.155299},
    {"PAM30", 14, 2, 0.337, 0.27, 1.62, 0.21, -0.8, 0.440010, 0.206970, 0.198524},
    {"PAM30", 14, 1, 0.333, 0.27, 1.43, 0.23, -1.4, 0.444817, 0.436301, 0.361947},
    {"PAM30", 13, 3, 0.338, 0.27, 1.69, 0.20, -0.5, 0.439086, 0.178973, 0.175436},
    // PAM70 (pam70_values)
    {"PAM70", 32767, 32767, 0.3345, 0.229, 1.029, 0.3250, -0.7, 0.511296, 0.633439, 0.633439},
    {"PAM70", 8, 2, 0.301, 0.12, 0.54, 0.56, -5, 0.549019, 2.881650, 3.025710},
    {"PAM70", 7, 2, 0.286, 0.093, 0.43, 0.67, -7, 0.565659, 4.534540, 4.785780},
    {"PAM70", 6, 2, 0.264, 0.064, 0.29, 0.90, -12, 0.596330, 7.942630, 8.402720},
    {"PAM70", 11, 1, 0.305, 0.12, 0.52, 0.59, -6, 0.543514, 3.681400, 4.108020},
    {"PAM70", 10, 1, 0.291, 0.091, 0.41, 0.71, -9, 0.560723, 6.002970, 6.716570},
    {"PAM70", 9, 1, 0.270, 0.060, 0.28, 0.97, -14, 0.585186, 11.360800, 12.636700},
    {"PAM70", 11, 2, 0.323, 0.186, 0.80, 1.32, -27, 0.524062, 1.321301, 1.281671},
    {"PAM70", 12, 3, 0.330, 0.219, 0.93, 0.82, -16, 0.516845, 0.818768, 0.811240},
    // PAM250 (pam250_values)
    {"PAM250", 32767, 32767, 0.2252, 0.0868, 0.2223, 0.98, -5.0, 0.660059, 11.754300, 11.754300},
    {"PAM250", 15, 3, 0.205, 0.049, 0.13, 1.6, -23, 0.687656, 34.578400, 34.928000},
    {"PAM250", 14, 3, 0.200, 0.043, 0.12, 1.7, -26, 0.689768, 43.353000, 43.443800},
    {"PAM250", 13, 3, 0.194, 0.036, 0.10, 1.9, -31, 0.697431, 50.948500, 51.081700},
    {"PAM250", 12, 3, 0.186, 0.029, 0.085, 2.2, -41, 0.704565, 69.606500, 69.793600},
    {"PAM250", 11, 3, 0.174, 0.020, 0.070, 2.5, -48, 0.722438, 98.653500, 98.927100},
    {"PAM250", 17, 2, 0.204, 0.047, 0.12, 1.7, -28, 0.684799, 41.583800, 41.735800},
    {"PAM250", 16, 2, 0.198, 0.038, 0.11, 1.8, -29, 0.691098, 51.635200, 51.843900},
    {"PAM250", 15, 2, 0.191, 0.031, 0.087, 2.2, -44, 0.699051, 67.256700, 67.558500},
    {"PAM250", 14, 2, 0.182, 0.024, 0.073, 2.5, -53, 0.714103, 96.315100, 96.756800},
    {"PAM250", 13, 2, 0.171, 0.017, 0.059, 2.9, -64, 0.728738, 135.653000, 136.339000},
    {"PAM250", 21, 1, 0.205, 0.045, 0.11, 1.8, -34, 0.683265, 48.728200, 49.218800},
    {"PAM250", 20, 1, 0.199, 0.037, 0.10, 1.9, -35, 0.689380, 60.832000, 61.514100},
    {"PAM250", 19, 1, 0.192, 0.029, 0.083, 2.3, -52, 0.696344, 84.019700, 84.985600},
    {"PAM250", 18, 1, 0.183, 0.021, 0.070, 2.6, -60, 0.710525, 113.829000, 115.184000},
    {"PAM250", 17, 1, 0.171, 0.014, 0.052, 3.3, -86, 0.727000, 175.071000, 177.196000},
};

}  // namespace detail

// The rows of NCBI's table for a matrix (the ungapped row first); empty for a matrix it lacks.
inline std::vector<GappedParameters> gapped_table(std::string_view matrix) {
  std::vector<GappedParameters> rows;
  for (const auto& g : detail::GAPPED)
    if (g.matrix == matrix) rows.push_back(g);
  return rows;
}

// The parameters of a matrix with affine gaps (open, extend), as BLAST tabulates them.
inline std::optional<GappedParameters> gapped_parameters(std::string_view matrix, int gap_open, int gap_extend) {
  for (const auto& g : detail::GAPPED)
    if (g.matrix == matrix && g.gap_open == gap_open && g.gap_extend == gap_extend) return g;
  return std::nullopt;
}

}  // namespace genomics::align
