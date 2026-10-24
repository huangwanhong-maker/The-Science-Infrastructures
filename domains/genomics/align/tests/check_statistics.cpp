// Checks of statistics.hpp: λ, K and H of the eight NCBI matrices under Robinson and Robinson's
// frequencies against NCBI BLAST's ungapped values (BLOSUM62: λ = 0.3176, K = 0.134, H = 0.4012);
// the series for K against its closed forms on walks that reach only −1 or only +1, and λ of the
// ±1 walk against ln(q/p); the target frequencies e^{λs} p q summing to 1; E = m n 2^{−S′}; the
// gapped table of BLAST (BLASTP's default BLOSUM62 11/1: λ = 0.267, K = 0.041, H = 0.14); and the
// theory against simulation: the maximal ungapped scores of random protein pairs follow the
// extreme-value law with these λ and K.

#include <genomics/align/scoring.hpp>
#include <genomics/align/statistics.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::align;

int main() {
  calx::verification::Suite s;
  const auto p = robinson_frequencies();

  // ---- the ungapped parameters of the NCBI matrices --------------------------------------------------
  {
    // NCBI BLAST's ungapped rows (blast_stat.c, the first row of each <matrix>_values, retrieved
    // 2026-10-24): λ, K, H, which BLAST computes from these matrices and frequencies.
    // λ and H with half a unit in the last digit BLAST gives as the limit; K, whose series BLAST's
    // table summed differently (its values lie between this full series and the 10⁻⁴ truncation
    // of today's BlastKarlinLHtoK), within one unit. BLAST's BLOSUM80 is its built-in half-bit
    // matrix, BLOSUM80-BLAST here (scoring.hpp).
    struct Row {
      const char* name;
      double lambda, K, H, dK, dH;
    };
    const Row rows[] = {{"BLOSUM45", 0.2291, 0.0924, 0.2514, 5e-5, 5e-5},     {"BLOSUM50", 0.2318, 0.112, 0.3362, 5e-4, 5e-5}, {"BLOSUM62", 0.3176, 0.134, 0.4012, 5e-4, 5e-5},
                        {"BLOSUM80-BLAST", 0.3430, 0.177, 0.6568, 5e-4, 5e-5}, {"BLOSUM90", 0.3346, 0.190, 0.7547, 5e-4, 5e-5}, {"PAM30", 0.3400, 0.283, 1.754, 5e-4, 5e-4},
                        {"PAM70", 0.3345, 0.229, 1.029, 5e-4, 5e-4},           {"PAM250", 0.2252, 0.0868, 0.2223, 5e-5, 5e-5}};
    for (const auto& row : rows) {
      const auto d = score_distribution(matrix(row.name), p, p);
      const KarlinAltschul k = karlin_altschul(d);
      s.check(std::string(row.name) + ": λ against NCBI BLAST's", k.lambda - row.lambda, 5e-5);
      s.check(std::string(row.name) + ": K within one unit of the last digit of NCBI BLAST's", k.K - row.K, 2 * row.dK);
      // BLAST's run-time truncation can only add to σ's omitted tail, so it never lowers K.
      s.require(std::string(row.name) + ": K truncated as BLAST truncates it is at least the full series, by under 0.5%",
                karlin_K_blast(d, k.lambda, k.H) >= k.K && karlin_K_blast(d, k.lambda, k.H) < 1.005 * k.K);
      s.check(std::string(row.name) + ": H (nats) against NCBI BLAST's", k.H - row.H, row.dH);
    }
    // The FTP file's BLOSUM80 is at a third of a bit: its λ is near ln 2 / 3, not BLAST's 0.343.
    const KarlinAltschul third = karlin_altschul(matrix("BLOSUM80"), p, p);
    s.check("BLOSUM80 (NCBI's file, 1/3 bit): λ near ln 2 / 3 = 0.2310, relative", third.lambda / (std::numbers::ln2 / 3) - 1, 0.05);
    s.require("BLOSUM80: BLAST's table is filed under the matrix it belongs to, BLOSUM80-BLAST", gapped_table("BLOSUM80").empty() && gapped_table("BLOSUM80-BLAST").size() == 10);
    const KarlinAltschul k = karlin_altschul(matrix("BLOSUM62"), p, p);
    const auto table = gapped_table("BLOSUM62");
    s.require("BLOSUM62: the table's ungapped row is λ 0.3176, K 0.134, H 0.4012", !table.empty() && table[0].ungapped() && table[0].lambda == 0.3176 && table[0].K == 0.134 && table[0].H == 0.4012);
    // The target frequencies implied by λ (Altschul 1991) are a distribution.
    const Matrix B = matrix("BLOSUM62");
    double total = 0;
    for (int a = 0; a < 24; a++)
      for (int b = 0; b < 24; b++) total += p[static_cast<std::size_t>(a)] * p[static_cast<std::size_t>(b)] * std::exp(k.lambda * B(static_cast<Code>(a), static_cast<Code>(b)));
    s.check("BLOSUM62: the target frequencies p_a p_b e^{λ s_ab} sum to 1", total - 1, 1e-12);
    const double bits = bit_score(50, k);
    s.check("E-values: E = K m n e^{−λS} = m n 2^{−S′} for S = 50, m = 300, n = 10⁶", evalue(50, 300, 1e6, k) / (300 * 1e6 * std::pow(2.0, -bits)) - 1, 1e-12);
  }

  // ---- closed forms ---------------------------------------------------------------------------------------
  {
    // The ±1 walk with P(+1) = 0.3, P(−1) = 0.7: λ = ln(0.7/0.3), and K = (q − p)²/q.
    const auto walk = score_distribution(-1, {0.7, 0.0, 0.3});
    const double lambda = karlin_lambda(walk), H = karlin_H(walk, lambda);
    s.check("±1 walk: λ = ln(q/p)", lambda - std::log(0.7 / 0.3), 1e-13);
    s.check("±1 walk: the series for K equals the closed form (q − p)²/q", karlin_K_series(walk, lambda, H) / ((0.7 - 0.3) * (0.7 - 0.3) / 0.7) - 1, 1e-8);
    // Highest score +1: {−2: 0.5, 0: 0.25, +1: 0.25}.
    const auto up = score_distribution(-2, {0.5, 0.0, 0.25, 0.25});
    const double lu = karlin_lambda(up), hu = karlin_H(up, lu);
    s.check("highest score +1: the series equals the closed form (H/λ)(1 − e^{−λ})", karlin_K_series(up, lu, hu) / karlin_K(up, lu, hu) - 1, 1e-8);
    // Lowest score −1: {−1: 0.7, 0: 0.05, +2: 0.25}.
    const auto down = score_distribution(-1, {0.7, 0.05, 0.0, 0.25});
    const double ld = karlin_lambda(down), hd = karlin_H(down, ld);
    s.check("lowest score −1: the series equals the closed form (μ²λ/H)(1 − e^{−λ})", karlin_K_series(down, ld, hd) / karlin_K(down, ld, hd) - 1, 1e-8);
    // Scores with a common divisor: multiplying every score by 4 divides λ by 4 and leaves K and H.
    const auto doubled = score_distribution(-4, {0.7, 0, 0, 0, 0.05, 0, 0, 0, 0, 0, 0, 0, 0.25});  // the scores of `down` times 4
    const auto kd = karlin_altschul(doubled), k1 = karlin_altschul(down);
    s.check("scaling: multiplying every score by 4 divides λ by 4 and keeps K and H", std::max({std::fabs(kd.lambda * 4 - k1.lambda), std::fabs(kd.K - k1.K), std::fabs(kd.H - k1.H)}), 1e-10);
    bool refused = false;
    try {
      karlin_lambda(score_distribution(-1, {0.4, 0.0, 0.6}));
    } catch (const std::invalid_argument&) { refused = true; }
    s.require("a scoring of positive mean is refused", refused);
  }

  // ---- BLAST's gapped table ---------------------------------------------------------------------------------
  {
    const auto g = gapped_parameters("BLOSUM62", 11, 1);
    s.require("gapped: BLOSUM62 with open 11, extend 1 (BLASTP's default) has λ 0.267, K 0.041, H 0.14, α 1.9, β −30",
              g && g->lambda == 0.267 && g->K == 0.041 && g->H == 0.14 && g->alpha == 1.9 && g->beta == -30);
    bool sizes = gapped_table("BLOSUM45").size() == 14 && gapped_table("BLOSUM50").size() == 16 && gapped_table("BLOSUM62").size() == 12 && gapped_table("BLOSUM80-BLAST").size() == 10 &&
                 gapped_table("BLOSUM90").size() == 8 && gapped_table("PAM30").size() == 11 && gapped_table("PAM70").size() == 9 && gapped_table("PAM250").size() == 16;
    s.require("gapped: the table has BLAST's 14, 16, 12, 10, 8, 11, 9 and 16 rows for its eight matrices", sizes);
    s.require("gapped: a pair of costs BLAST does not tabulate is absent", !gapped_parameters("BLOSUM62", 5, 5).has_value());
  }

  // ---- against simulation -------------------------------------------------------------------------------------
  {
    // The best ungapped local score S of two random proteins of length 400 (Robinson and Robinson's
    // frequencies, BLOSUM62). Karlin and Altschul: P(S ≥ x) ≈ 1 − exp(−K m′ n′ e^{−λx}), with the
    // edge-corrected lengths m′ = m − ℓ, ℓ = ln(K m n)/H (Altschul and Gish 1996). Compared: the
    // fraction of 3000 pairs with S ≥ x at the x where the prediction is 0.1, and the mean of S
    // against the Gumbel mean (ln(K m′ n′) + γ)/λ, less 1/2 for integer scores (the sum
    // Σ_{x≥1} P(S ≥ x) is the integral of the continuous tail from 1/2, by the midpoint rule).
    const Matrix B = matrix("BLOSUM62");
    const KarlinAltschul k = karlin_altschul(B, p, p);
    calx::random::Stream r(20261024, 0x6b61726c, 0, 0);
    std::vector<double> cumulative(20);
    double acc = 0;
    for (int a = 0; a < 20; a++) cumulative[static_cast<std::size_t>(a)] = acc += p[static_cast<std::size_t>(a)];
    const auto draw = [&](std::size_t n) {
      Codes x(n);
      for (auto& c : x) {
        const double u = r.uniform() * acc;
        c = static_cast<Code>(std::lower_bound(cumulative.begin(), cumulative.end(), u) - cumulative.begin());
      }
      return x;
    };
    const std::size_t m = 400, n = 400, pairs = 3000;
    const double ell = std::log(k.K * m * n) / k.H, me = m - ell, ne = n - ell;
    const double x01 = std::log(k.K * me * ne / -std::log(1 - 0.1)) / k.lambda;  // P(S ≥ x) = 0.1
    const double threshold = std::ceil(x01);
    const double predicted = pvalue(evalue(threshold, me, ne, k));
    std::size_t above = 0;
    double mean = 0;
    for (std::size_t t = 0; t < pairs; t++) {
      const Codes x = draw(m), y = draw(n);
      std::int64_t best = 0;
      for (std::int64_t d = -static_cast<std::int64_t>(m) + 1; d < static_cast<std::int64_t>(n); d++) {
        std::int64_t run = 0;
        for (std::size_t i = d < 0 ? static_cast<std::size_t>(-d) : 0, j = i + static_cast<std::size_t>(d); i < m && j < n; i++, j++) {
          run = std::max<std::int64_t>(0, run + B(x[i], y[j]));
          best = std::max(best, run);
        }
      }
      above += best >= threshold;
      mean += static_cast<double>(best) / pairs;
    }
    const double observed = static_cast<double>(above) / pairs;
    const double sd = std::sqrt(predicted * (1 - predicted) / pairs);
    s.check("simulation: P(S ≥ " + std::to_string(static_cast<int>(threshold)) + ") of 3000 random pairs of 400 residues against Karlin–Altschul, in standard errors",
            (observed - predicted) / sd, 3.0, "σ");
    const double gumbel = (std::log(k.K * me * ne) + 0.5772156649015329) / k.lambda - 0.5;
    s.check("simulation: the mean best ungapped score against the Gumbel mean (ln(K m′ n′) + γ)/λ − 1/2, relative", (mean - gumbel) / gumbel, 0.02);
  }
  return s.done();
}
