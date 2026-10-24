// CpG islands: the two Markov chains of Durbin, Eddy, Krogh and Mitchison (1998, Biological
// Sequence Analysis, §3.1), their log-odds, and the eight-state hidden Markov model that finds
// islands along a sequence (§3.2).
//
// In vertebrate genomes the dinucleotide CG is rare, because a methylated C in CG mutates to T
// at a high rate; around promoters methylation is suppressed and CG keeps its expected
// frequency, in stretches of a few hundred to a few thousand bases (Bird 1987, Trends in
// Genetics 3, 342–347; Gardiner-Garden and Frommer 1987).
//
// - The chains. a⁺_st and a⁻_st, the probabilities that base t follows base s inside islands
//   (+) and outside them (−), estimated by Durbin et al. from putative islands in human
//   sequence and the regions between them (the tables on p. 50, transcribed below; their rows
//   of three decimals sum to 1 ± 0.001 and are normalised here). Transcription checked 2026-10-24 against Clemens Gröpl's
//   lecture notes (FU Berlin, Hidden Markov Models, 2012, §9.4, "Source: DEMK98, p 50") and
//   the Stanford Stats 366 notes (web.stanford.edu/class/stats366/exs/HMM1.html).
// - The log-odds of a window, S(x) = Σ_i log₂ (a⁺_{x_{i−1}x_i} / a⁻_{x_{i−1}x_i}) = Σ_i β_{x_{i−1}x_i}
//   in bits (§3.1; the table of β on p. 51), and S(x)/L, the length-normalised score of the
//   book's Figure 3.2. Pairs with an N are skipped.
// - The HMM. States A⁺ C⁺ G⁺ T⁺ (0–3) and A⁻ C⁻ G⁻ T⁻ (4–7), each emitting its own base. Within
//   a group, transitions are those of its chain scaled by 1 − p (inside) or 1 − q (outside);
//   leaving, the label changes first and the next base is drawn from the new chain:
//   a(s⁺ → t⁻) = p a⁻_st and a(s⁻ → t⁺) = q a⁺_st. Then the label is a Markov chain with mean
//   sojourns 1/p in islands and 1/q between them, and, given the labels, the bases are the
//   chain of each label. The book leaves p and q to the modeller (§3.2); the defaults are an
//   island of 1 kb on average and 100 kb between islands, the order of magnitude of the human
//   genome, whose 3 Gb hold some tens of thousands of islands (Lander et al. 2001, Nature 409,
//   860–921). The
//   initial distribution is the label's stationary distribution times each chain's.
// - Decoding. Islands are the maximal runs of + states of the Viterbi path, or of positions
//   whose posterior probability of + exceeds a threshold (posterior decoding, §3.2); each is
//   reported with its log-odds, GC content and CpG observed/expected.
//
// Complexity: O(L) for the log-odds; O(64 L) for decoding (eight states), with O(8 L) bytes for
// Viterbi's back-pointers or O(64 L) bytes for the posterior.
//
// Limits: the chains were estimated from a small set of human sequences; they describe human
// islands, and are a teaching model rather than a genome annotator (which would add a minimum
// length and the GC and o/e criteria of Takai and Jones 2002, PNAS 99, 3740–3745).

#pragma once

#include <genomics/hmm/model.hpp>
#include <genomics/sequence/composition.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace genomics::hmm {

// A first-order Markov chain on A, C, G, T: a[4 s + t] = P(t | s).
struct MarkovChain {
  std::array<double, 16> a{};
  double operator()(Code s, Code t) const { return a[4 * s + t]; }
  // Rows scaled to sum to one.
  MarkovChain normalised() const {
    MarkovChain c = *this;
    for (int s = 0; s < 4; s++) {
      double r = 0;
      for (int t = 0; t < 4; t++) r += c.a[4 * s + t];
      for (int t = 0; t < 4; t++) c.a[4 * s + t] /= r;
    }
    return c;
  }
  // The stationary distribution, by power iteration (the chains here are primitive).
  std::array<double, 4> stationary() const {
    const MarkovChain c = normalised();
    std::array<double, 4> p{0.25, 0.25, 0.25, 0.25};
    for (int it = 0; it < 1000; it++) {
      std::array<double, 4> q{};
      for (int s = 0; s < 4; s++)
        for (int t = 0; t < 4; t++) q[t] += p[s] * c.a[4 * s + t];
      double d = 0;
      for (int t = 0; t < 4; t++) d += std::fabs(q[t] - p[t]);
      p = q;
      if (d < 1e-15) break;
    }
    return p;
  }
};

// Durbin et al. 1998, p. 50: the + model (inside islands) and the − model, rows A, C, G, T.
inline constexpr MarkovChain CPG_PLUS{{0.180, 0.274, 0.426, 0.120,  //
                                       0.171, 0.368, 0.274, 0.188,  //
                                       0.161, 0.339, 0.375, 0.125,  //
                                       0.079, 0.355, 0.384, 0.182}};
inline constexpr MarkovChain CPG_MINUS{{0.300, 0.205, 0.285, 0.210,  //
                                        0.322, 0.298, 0.078, 0.302,  //
                                        0.248, 0.246, 0.298, 0.208,  //
                                        0.177, 0.239, 0.292, 0.292}};

// β_st = log₂ (a⁺_st / a⁻_st), in bits, of the normalised chains.
inline std::array<double, 16> log_odds_table(const MarkovChain& plus = CPG_PLUS, const MarkovChain& minus = CPG_MINUS) {
  const MarkovChain p = plus.normalised(), m = minus.normalised();
  std::array<double, 16> b{};
  for (int i = 0; i < 16; i++) b[i] = std::log2(p.a[i] / m.a[i]);
  return b;
}

// S(x) over the pairs (x_{i−1}, x_i) with begin < i < end, in bits.
inline double log_odds(const Codes& x, std::size_t begin = 0, std::size_t end = SIZE_MAX, const std::array<double, 16>& beta = log_odds_table()) {
  end = std::min(end, x.size());
  double s = 0;
  for (std::size_t i = begin + 1; i < end; i++)
    if (x[i - 1] < 4 && x[i] < 4) s += beta[4 * x[i - 1] + x[i]];
  return s;
}

// S(x)/L in windows, bits per base.
inline std::vector<sequence::Window> log_odds_windows(const Codes& x, std::size_t size, std::size_t step) {
  const auto beta = log_odds_table();
  std::vector<sequence::Window> out;
  if (size == 0 || step == 0) return out;
  for (std::size_t b = 0; b < x.size(); b += step) {
    const std::size_t e = std::min(x.size(), b + size);
    out.push_back({b, e, e > b ? log_odds(x, b, e, beta) / static_cast<double>(e - b) : 0});
    if (e == x.size()) break;
  }
  return out;
}

// The maximum-likelihood chain of a set of sequences, a_st = (c_st + r) / Σ_t' (c_st' + r)
// (Durbin et al. §3.1), pairs with an N skipped.
inline MarkovChain estimate_chain(const std::vector<Codes>& xs, double pseudocount = 0) {
  std::array<double, 16> c{};
  for (const auto& x : xs)
    for (std::size_t i = 1; i < x.size(); i++)
      if (x[i - 1] < 4 && x[i] < 4) c[4 * x[i - 1] + x[i]] += 1;
  MarkovChain m;
  for (int s = 0; s < 4; s++) {
    double r = 0;
    for (int t = 0; t < 4; t++) r += c[4 * s + t] + pseudocount;
    for (int t = 0; t < 4; t++) m.a[4 * s + t] = r > 0 ? (c[4 * s + t] + pseudocount) / r : 0.25;
  }
  return m;
}

// The eight-state model: states 0–3 are A⁺ … T⁺, 4–7 A⁻ … T⁻.
inline Model cpg_model(double mean_island = 1000, double mean_background = 100000, const MarkovChain& plus = CPG_PLUS, const MarkovChain& minus = CPG_MINUS) {
  if (!(mean_island >= 1) || !(mean_background >= 1)) throw std::invalid_argument("genomics::hmm::cpg_model: mean lengths of at least one base");
  const double p = 1 / mean_island, q = 1 / mean_background;
  const MarkovChain P = plus.normalised(), M = minus.normalised();
  Model m(8, 4);
  const auto sp = P.stationary(), sm = M.stationary();
  const double in = q / (p + q);  // the stationary probability of the + label
  for (int s = 0; s < 4; s++) {
    m.initial[s] = in * sp[s];
    m.initial[4 + s] = (1 - in) * sm[s];
    m.e(s, s) = 1;
    m.e(4 + s, s) = 1;
    for (int t = 0; t < 4; t++) {
      m.a(s, t) = (1 - p) * P(s, t);
      m.a(s, 4 + t) = p * M(s, t);
      m.a(4 + s, 4 + t) = (1 - q) * M(s, t);
      m.a(4 + s, t) = q * P(s, t);
    }
  }
  return m;
}

struct Island {
  std::size_t begin = 0, end = 0;
  double log_odds = 0;   // S(x) of the island, bits
  double gc = 0, cpg_observed_expected = 0;
};

struct IslandOptions {
  double mean_island = 1000, mean_background = 100000;
  bool posterior = false;    // posterior decoding instead of Viterbi
  double threshold = 0.5;    // the posterior probability of + that calls a base, when posterior
  std::size_t min_length = 1;
};

// The maximal runs of positions where `inside` holds, as islands with their statistics.
template <class Inside>
std::vector<Island> runs(const Codes& x, std::size_t n, Inside&& inside, std::size_t min_length) {
  std::vector<Island> out;
  const auto beta = log_odds_table();
  std::size_t i = 0;
  while (i < n) {
    if (!inside(i)) { i++; continue; }
    std::size_t j = i;
    while (j < n && inside(j)) j++;
    if (j - i >= min_length) {
      Island is;
      is.begin = i;
      is.end = j;
      is.log_odds = log_odds(x, i, j, beta);
      is.gc = sequence::gc_fraction(x, i, j);
      is.cpg_observed_expected = sequence::cpg_observed_expected(x, i, j);
      out.push_back(is);
    }
    i = j;
  }
  return out;
}

// The islands of x decoded by the eight-state model.
inline std::vector<Island> cpg_islands(const Codes& x, const IslandOptions& o = {}) {
  const Model m = cpg_model(o.mean_island, o.mean_background);
  if (o.posterior) {
    const Posterior p = posterior(m, x);
    if (p.log_likelihood == NEG_INF) return {};
    std::vector<double> plus(x.size());
    for (std::size_t t = 0; t < x.size(); t++) plus[t] = p.at(t, 0) + p.at(t, 1) + p.at(t, 2) + p.at(t, 3);
    return runs(x, x.size(), [&](std::size_t t) { return plus[t] > o.threshold; }, o.min_length);
  }
  const Path v = viterbi(m, x);
  if (v.states.empty()) return {};
  return runs(x, x.size(), [&](std::size_t t) { return v.states[t] < 4; }, o.min_length);
}

}  // namespace genomics::hmm
