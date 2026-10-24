// Motifs of DNA: position weight matrices, their information, the exact distribution of
// their scores, scanning both strands, and the discovery of a motif shared by a set of
// sequences, by Gibbs sampling and by expectation maximisation.
//
// - A PWM of width W from aligned sites: counts n_{i,b} at each position i, probabilities
//   p_{i,b} = (n_{i,b} + α q_b)/(N + α) with a total pseudocount α spread by the background q,
//   and weights w_{i,b} = log₂(p_{i,b}/q_b) in bits, so that a word's score Σ_i w_{i,x_i} is
//   the log-odds of the motif against the background (Stormo 2000, Bioinformatics 16, 16–23).
//   A window with an N has no score and is skipped; on the reverse strand the score is that
//   of the window's reverse complement.
// - Information. The relative entropy of each column, I_i = Σ_b p_{i,b} log₂(p_{i,b}/q_b), the
//   mean score of a site (Stormo 2000); and Schneider's R_sequence = H_g − (H_i + e(n)) from the
//   counts, with H_g the background's entropy, H_i the column's, and e(n) the small-sample
//   correction, the bias E[H_g − H_n] of the entropy of n bases drawn from the background
//   (Schneider, Stormo, Gold and Ehrenfeucht 1986, J. Mol. Biol. 188, 415–431): exact by summing
//   over the compositions of n for n ≤ 150, and (s − 1)/(2 n ln 2) above, Miller's (1955)
//   first-order bias. The heights of a sequence logo are f_{i,b} · R_sequence(i) (Schneider and
//   Stephens 1990, Nucleic Acids Res. 18, 6097–6100). A column of one base has 2 bits.
// - The exact distribution of scores of a random word drawn from the background, by dynamic
//   programming over a lattice: each weight rounded to an integer multiple of a granularity ε,
//   the distribution of the rounded score built column by column, P_i(s) = Σ_b q_b P_{i−1}(s −
//   W_{i,b}) (Staden 1989, Comput. Appl. Biosci. 5, 89–96; Touzet and Varré 2007, Algorithms
//   Mol. Biol. 2, 15, whose analysis bounds the rounding). The p-value of
//   a word's rounded score is exact for the rounded matrix; a real score t is bracketed by
//   P(S' ≥ t + Wε/2) ≤ P(S ≥ t) ≤ P(S' ≥ t − Wε/2), since every word's rounding error is at most
//   Wε/2.
// - The Gibbs motif sampler (Lawrence, Altschul, Boguski, Liu, Neuwald and Wootton 1993, Science
//   262, 208–214): one site in each sequence (OOPS); in turn each sequence's site is removed, the
//   PWM q_{i,b} = (c_{i,b} + B p_b)/(N − 1 + B) of the others is built (pseudocounts B = √(N − 1)
//   as the paper's √N), and a new site is drawn with probability proportional to the odds
//   Π_i q_{i,x_i}/p_{x_i} of each window; every few sweeps a phase shift moves all sites by δ,
//   drawn in proportion to the likelihood ratio of the shifted alignment, so that the sampler
//   does not stay on a motif seen a few bases off. The alignment of largest
//   F = Σ_i Σ_b c_{i,b} log₂(q_{i,b}/p_b) over all sweeps and restarts is returned.
// - MEME's expectation maximisation (Bailey and Elkan 1994, Proc. ISMB 2, 28–36; 1995, Machine
//   Learning 21, 51–80) for the one-occurrence-per-sequence (OOPS) and zero-or-one (ZOOPS) models:
//   the E-step gives each window the posterior Z_ij of being the site, Z_ij ∝ R_ij = Π_k
//   p_{k,x}/p_0(x) for OOPS, and Z_ij = λ R_ij / ((1 − γ) + λ Σ_j' R_ij') with λ = γ/m_i for ZOOPS;
//   the M-step re-estimates p from the expected counts (with pseudocounts β q_b) and γ as the mean
//   of the Z; with the pseudocounts, a Dirichlet prior, EM climbs the log posterior, which the
//   trace records. Starting points are MEME's: every distinct word of width W in the data, as a matrix
//   giving its letter probability X = 0.5 and the others (1 − X)/3; one iteration from each, and
//   EM to convergence from the start of largest likelihood.
//
// Both discovery methods are deterministic: the sampler's draws come from calx::random::Stream
// keyed by the seed and the restart, restarts run on a ThreadPool as shards and are compared in
// order, and EM has no randomness.
//
// Complexity: scanning O(L W) per strand; the score distribution O(W R) for a lattice of R
// points (R ≈ W · (max − min weight)/ε); a Gibbs sweep O(Σ L_i W); an EM iteration O(Σ L_i W),
// MEME's starting points O(U Σ L_i W) for U distinct words.
//
// Limits: one motif at a time (MEME's erasing of found sites for the next motif and the TCM
// model are not here); the background is of order zero; e(n) assumes the counts are of
// independent sites. *Caveat.* The maximum-likelihood ZOOPS fit calls weak sites in sequences
// that have none (the checks show its likelihood above that of the planted truth); MEME itself
// chooses the number of sites by the E-value of the motif, which is not here.

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/special_functions/gamma.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

namespace genomics::annotation {

using sequence::Code;
using sequence::Codes;

inline constexpr std::uint32_t GIBBS = 0x47494253;  // "GIBS", the sampler's streams
inline constexpr double MOTIF_NEG_INF = -std::numeric_limits<double>::infinity();

using Background = std::array<double, 4>;
inline constexpr Background UNIFORM{0.25, 0.25, 0.25, 0.25};

// The base frequencies of sequences (N skipped), with one pseudocount each.
inline Background base_frequencies(const std::vector<Codes>& xs) {
  Background q{1, 1, 1, 1};
  for (const auto& x : xs)
    for (Code c : x)
      if (c < 4) q[c] += 1;
  const double s = q[0] + q[1] + q[2] + q[3];
  for (auto& v : q) v /= s;
  return q;
}

// ---- position weight matrices ------------------------------------------------------------------------
struct Pwm {
  std::size_t width = 0;
  std::vector<double> counts;       // n_{i,b}, W × 4
  std::vector<double> probability;  // p_{i,b}
  std::vector<double> weight;       // w_{i,b} = log₂(p_{i,b}/q_b), bits
  Background background = UNIFORM;
  double sites = 0;                 // N

  double w(std::size_t i, Code b) const { return weight[4 * i + b]; }
  double p(std::size_t i, Code b) const { return probability[4 * i + b]; }
  double max_score() const {
    double s = 0;
    for (std::size_t i = 0; i < width; i++) s += *std::max_element(&weight[4 * i], &weight[4 * i] + 4);
    return s;
  }
  double min_score() const {
    double s = 0;
    for (std::size_t i = 0; i < width; i++) s += *std::min_element(&weight[4 * i], &weight[4 * i] + 4);
    return s;
  }
  // The score of the window at `at` (forward-strand coordinates), or of its reverse complement;
  // −∞ when the window has an N or leaves the sequence.
  double score(const Codes& x, std::size_t at, bool reverse = false) const {
    if (at + width > x.size()) return MOTIF_NEG_INF;
    double s = 0;
    for (std::size_t i = 0; i < width; i++) {
      const Code c = reverse ? sequence::complement(x[at + width - 1 - i]) : x[at + i];
      if (c > 3) return MOTIF_NEG_INF;
      s += weight[4 * i + c];
    }
    return s;
  }
};

// A PWM from counts (W × 4), with a total pseudocount α spread by the background.
inline Pwm pwm_from_counts(const std::vector<double>& counts, double pseudocount = 1, const Background& q = UNIFORM) {
  if (counts.empty() || counts.size() % 4) throw std::invalid_argument("genomics::annotation::pwm_from_counts: counts of four bases at each position");
  Pwm m;
  m.width = counts.size() / 4;
  m.counts = counts;
  m.background = q;
  m.probability.assign(counts.size(), 0);
  m.weight.assign(counts.size(), 0);
  for (std::size_t i = 0; i < m.width; i++) {
    double n = 0;
    for (int b = 0; b < 4; b++) n += counts[4 * i + b];
    m.sites = std::max(m.sites, n);
    for (int b = 0; b < 4; b++) {
      const double p = n + pseudocount > 0 ? (counts[4 * i + b] + pseudocount * q[b]) / (n + pseudocount) : q[b];
      m.probability[4 * i + b] = p;
      m.weight[4 * i + b] = p > 0 ? std::log2(p / q[b]) : MOTIF_NEG_INF;
    }
  }
  return m;
}

// A PWM from aligned sites of one width (an N adds nothing to its column).
inline Pwm pwm_from_sites(const std::vector<Codes>& sites, double pseudocount = 1, const Background& q = UNIFORM) {
  if (sites.empty()) throw std::invalid_argument("genomics::annotation::pwm_from_sites: at least one site");
  const std::size_t W = sites[0].size();
  std::vector<double> c(4 * W, 0);
  for (const auto& s : sites) {
    if (s.size() != W) throw std::invalid_argument("genomics::annotation::pwm_from_sites: sites of one width");
    for (std::size_t i = 0; i < W; i++)
      if (s[i] < 4) c[4 * i + s[i]] += 1;
  }
  return pwm_from_counts(c, pseudocount, q);
}

// The relative entropy of each column, bits.
inline std::vector<double> information_content(const Pwm& m) {
  std::vector<double> out(m.width, 0);
  for (std::size_t i = 0; i < m.width; i++)
    for (int b = 0; b < 4; b++) {
      const double p = m.probability[4 * i + b];
      if (p > 0) out[i] += p * std::log2(p / m.background[b]);
    }
  return out;
}

inline double entropy(const Background& q) {
  double h = 0;
  for (double p : q)
    if (p > 0) h -= p * std::log2(p);
  return h;
}

// e(n) = H_g − E[H_n], the expected shortfall of the entropy of n bases drawn from q.
inline double small_sample_correction(std::size_t n, const Background& q = UNIFORM) {
  if (n == 0) return entropy(q);
  if (n > 150) return 3.0 / (2.0 * static_cast<double>(n) * std::numbers::ln2);
  const double N = static_cast<double>(n), lgn = calx::special_functions::log_gamma(N + 1);
  std::array<double, 4> lq{};
  for (int b = 0; b < 4; b++) lq[b] = std::log(q[b]);
  double expected = 0;
  for (std::size_t a = 0; a <= n; a++)
    for (std::size_t c = 0; a + c <= n; c++)
      for (std::size_t g = 0; a + c + g <= n; g++) {
        const std::size_t t = n - a - c - g;
        const std::size_t k[4] = {a, c, g, t};
        double lp = lgn, h = 0;
        for (int b = 0; b < 4; b++) {
          const double kb = static_cast<double>(k[b]);
          lp += kb * lq[b] - calx::special_functions::log_gamma(kb + 1);
          if (k[b]) h -= kb / N * std::log2(kb / N);
        }
        expected += std::exp(lp) * h;
      }
  return entropy(q) - expected;
}

struct LogoColumn {
  std::array<double, 4> height{};  // bits, by base
  double information = 0;          // R_sequence, bits
};

// The heights of a sequence logo from the matrix's counts.
inline std::vector<LogoColumn> logo(const Pwm& m, bool correct = true) {
  std::vector<LogoColumn> out(m.width);
  const double hg = entropy(m.background);
  for (std::size_t i = 0; i < m.width; i++) {
    double n = 0;
    for (int b = 0; b < 4; b++) n += m.counts[4 * i + b];
    if (!(n > 0)) continue;
    double h = 0;
    for (int b = 0; b < 4; b++) {
      const double f = m.counts[4 * i + b] / n;
      if (f > 0) h -= f * std::log2(f);
    }
    const double e = correct ? small_sample_correction(static_cast<std::size_t>(std::lround(n)), m.background) : 0;
    out[i].information = std::max(0.0, hg - (h + e));
    for (int b = 0; b < 4; b++) out[i].height[b] = m.counts[4 * i + b] / n * out[i].information;
  }
  return out;
}

// ---- the distribution of scores ------------------------------------------------------------------------
struct ScoreDistribution {
  double granularity = 0;            // ε, bits
  std::size_t width = 0;
  std::int64_t lowest = 0;           // the lattice index of mass[0]
  std::vector<double> mass;          // P(S' = (lowest + j) ε)
  std::vector<double> at_least;      // P(S' ≥ (lowest + j) ε)
  std::vector<std::int64_t> lattice; // the rounded weights W_{i,b}, W × 4; INT64_MIN for −∞
  double impossible = 0;             // P(S = −∞)

  static constexpr std::int64_t NONE = std::numeric_limits<std::int64_t>::min();

  // P(S' ≥ k ε).
  double tail_index(std::int64_t k) const {
    if (k <= lowest) return at_least.empty() ? 0 : at_least.front();
    const std::int64_t j = k - lowest;
    return j >= static_cast<std::int64_t>(at_least.size()) ? 0 : at_least[static_cast<std::size_t>(j)];
  }
  // P(S' ≥ t) for a real t: the lattice points at or above it.
  double tail(double t) const { return tail_index(static_cast<std::int64_t>(std::ceil(t / granularity - 1e-9))); }
  // Bounds on P(S ≥ t) for the unrounded score.
  std::pair<double, double> bounds(double t) const {
    const double slack = 0.5 * granularity * static_cast<double>(width);
    return {tail(t + slack), tail(t - slack)};
  }
  // The rounded score of a word (forward), and its exact p-value under the rounded matrix.
  std::int64_t rounded(const Codes& x, std::size_t at, bool reverse = false) const {
    std::int64_t s = 0;
    for (std::size_t i = 0; i < width; i++) {
      const Code c = reverse ? sequence::complement(x[at + width - 1 - i]) : x[at + i];
      if (c > 3) return NONE;
      const std::int64_t w = lattice[4 * i + c];
      if (w == NONE) return NONE;
      s += w;
    }
    return s;
  }
  double p_value(const Codes& x, std::size_t at, bool reverse = false) const {
    const std::int64_t s = rounded(x, at, reverse);
    return s == NONE ? 1.0 : tail_index(s);
  }
  // The smallest score t on the lattice with P(S' ≥ t) ≤ p.
  double threshold(double p) const {
    for (std::size_t j = 0; j < at_least.size(); j++)
      if (at_least[j] <= p) return static_cast<double>(lowest + static_cast<std::int64_t>(j)) * granularity;
    return static_cast<double>(lowest + static_cast<std::int64_t>(at_least.size())) * granularity;
  }
};

inline ScoreDistribution score_distribution(const Pwm& m, double granularity = 1e-3) {
  if (!(granularity > 0)) throw std::invalid_argument("genomics::annotation::score_distribution: a positive granularity");
  ScoreDistribution d;
  d.granularity = granularity;
  d.width = m.width;
  d.lattice.assign(4 * m.width, ScoreDistribution::NONE);
  for (std::size_t i = 0; i < 4 * m.width; i++)
    if (std::isfinite(m.weight[i])) d.lattice[i] = static_cast<std::int64_t>(std::llround(m.weight[i] / granularity));
  std::vector<double> cur{1.0}, next;
  std::int64_t lo = 0;
  for (std::size_t i = 0; i < m.width; i++) {
    std::int64_t wmin = std::numeric_limits<std::int64_t>::max(), wmax = std::numeric_limits<std::int64_t>::min();
    double lost = 0;
    for (int b = 0; b < 4; b++) {
      const std::int64_t w = d.lattice[4 * i + b];
      if (w == ScoreDistribution::NONE) { lost += m.background[b]; continue; }
      wmin = std::min(wmin, w);
      wmax = std::max(wmax, w);
    }
    if (wmin > wmax) { cur.assign(1, 0); d.impossible = 1; break; }
    next.assign(cur.size() + static_cast<std::size_t>(wmax - wmin), 0);
    for (int b = 0; b < 4; b++) {
      const std::int64_t w = d.lattice[4 * i + b];
      if (w == ScoreDistribution::NONE) continue;
      const double q = m.background[b];
      double* out = next.data() + (w - wmin);
      for (std::size_t s = 0; s < cur.size(); s++) out[s] += q * cur[s];
    }
    d.impossible += (1 - d.impossible) * lost;
    lo += wmin;
    std::swap(cur, next);
  }
  d.lowest = lo;
  d.mass = cur;
  d.at_least.assign(cur.size(), 0);
  double acc = 0;
  for (std::size_t j = cur.size(); j-- > 0;) d.at_least[j] = acc += cur[j];
  return d;
}

// ---- scanning ------------------------------------------------------------------------------------------
struct Hit {
  std::size_t position = 0;  // the window's first base, forward strand
  bool reverse = false;
  double score = 0;          // bits
  double p_value = 1;        // when a distribution is given
};

// The windows scoring at least `threshold` bits, by position then strand.
inline std::vector<Hit> scan(const Pwm& m, const Codes& x, double threshold, bool both_strands = true, const ScoreDistribution* d = nullptr) {
  std::vector<Hit> out;
  if (x.size() < m.width) return out;
  for (std::size_t j = 0; j + m.width <= x.size(); j++)
    for (int strand = 0; strand < (both_strands ? 2 : 1); strand++) {
      const double s = m.score(x, j, strand == 1);
      if (s >= threshold) out.push_back({j, strand == 1, s, d ? d->p_value(x, j, strand == 1) : 1.0});
    }
  return out;
}

// ---- discovery -----------------------------------------------------------------------------------------
struct Site {
  std::size_t sequence = 0, position = 0;  // position: the first base, forward strand
  bool reverse = false;
  double probability = 1;                  // MEME: the posterior Z of the site
};

struct Motif {
  Pwm pwm;
  std::vector<Site> sites;  // OOPS: one per sequence; ZOOPS: those with Z > 0.5
  double score = 0;         // Gibbs: F, bits; MEME: the log-likelihood ratio against the background, bits
  double gamma = 1;         // ZOOPS: the fraction of sequences with a site
  std::size_t iterations = 0;
  std::vector<double> trace;  // MEME: the log posterior at each iteration (the likelihood ratio plus Σ β q_b log p), bits: EM never lowers it
};

namespace detail {

// The windows of the sequences: for each, its valid starting positions on each strand.
struct Windows {
  std::size_t width = 0;
  bool both = false;
  std::vector<std::vector<std::uint32_t>> valid;  // positions without N, per sequence
};

inline Windows windows(const std::vector<Codes>& xs, std::size_t W, bool both) {
  Windows w;
  w.width = W;
  w.both = both;
  w.valid.resize(xs.size());
  for (std::size_t s = 0; s < xs.size(); s++) {
    const Codes& x = xs[s];
    std::size_t run = 0;
    for (std::size_t i = 0; i < x.size(); i++) {
      run = x[i] < 4 ? run + 1 : 0;
      if (run >= W) w.valid[s].push_back(static_cast<std::uint32_t>(i + 1 - W));
    }
  }
  return w;
}

inline Code letter(const Codes& x, std::size_t at, std::size_t i, std::size_t W, bool reverse) {
  return reverse ? sequence::complement(x[at + W - 1 - i]) : x[at + i];
}

inline void add_site(std::vector<double>& c, const Codes& x, std::size_t at, bool reverse, std::size_t W, double weight) {
  for (std::size_t i = 0; i < W; i++) c[4 * i + letter(x, at, i, W, reverse)] += weight;
}

// log(q_{i,b}/p_b) of counts c with pseudocounts B p_b over `n` sites.
inline void log_ratios(const std::vector<double>& c, double n, double B, const Background& p, std::vector<double>& lr) {
  const std::size_t W = c.size() / 4;
  lr.resize(c.size());
  for (std::size_t i = 0; i < W; i++)
    for (int b = 0; b < 4; b++) lr[4 * i + b] = std::log((c[4 * i + b] + B * p[b]) / (n + B) / p[b]);
}

inline double window_score(const std::vector<double>& lr, const Codes& x, std::size_t at, bool reverse, std::size_t W) {
  double s = 0;
  for (std::size_t i = 0; i < W; i++) s += lr[4 * i + letter(x, at, i, W, reverse)];
  return s;
}

}  // namespace detail

struct GibbsOptions {
  std::size_t width = 10;
  std::size_t sweeps = 200;        // each sweep resamples every sequence's site once
  std::size_t restarts = 8;
  std::size_t shift_every = 5;     // sweeps between phase shifts (0: none)
  std::size_t max_shift = 0;       // 0: W/2
  double pseudocount = -1;         // B; negative: √(N − 1)
  bool both_strands = false;
  std::uint64_t seed = 1;
};

inline Motif gibbs_motif(const std::vector<Codes>& xs, const GibbsOptions& o = {}, calx::scheduler::ThreadPool* pool = nullptr) {
  const std::size_t n = xs.size(), W = o.width;
  if (n < 2 || W == 0) throw std::invalid_argument("genomics::annotation::gibbs_motif: two sequences and a positive width");
  const detail::Windows win = detail::windows(xs, W, o.both_strands);
  for (const auto& v : win.valid)
    if (v.empty()) throw std::invalid_argument("genomics::annotation::gibbs_motif: every sequence needs a window of the width without N");
  const Background p = base_frequencies(xs);
  const double B = o.pseudocount >= 0 ? o.pseudocount : std::sqrt(static_cast<double>(n - 1));
  const int strands = o.both_strands ? 2 : 1;
  const std::size_t max_shift = o.max_shift ? o.max_shift : W / 2;

  struct Result {
    std::vector<Site> sites;
    double F = MOTIF_NEG_INF;
    std::size_t sweeps = 0;
  };
  std::vector<Result> results(std::max<std::size_t>(1, o.restarts));
  // F of an alignment, bits (all N sites, pseudocounts B).
  const auto F_of = [&](const std::vector<Site>& sites) {
    std::vector<double> c(4 * W, 0);
    for (const auto& s : sites) detail::add_site(c, xs[s.sequence], s.position, s.reverse, W, 1);
    double f = 0;
    for (std::size_t i = 0; i < W; i++)
      for (int b = 0; b < 4; b++) {
        const double q = (c[4 * i + b] + B * p[b]) / (static_cast<double>(n) + B);
        f += c[4 * i + b] * std::log2(q / p[b]);
      }
    return f;
  };
  const auto run = [&](std::size_t restart) {
    calx::random::Stream r(o.seed, GIBBS, static_cast<std::uint32_t>(restart), 0);
    std::vector<Site> sites(n);
    for (std::size_t s = 0; s < n; s++) {
      const auto& v = win.valid[s];
      sites[s] = {s, v[r.below(v.size())], strands == 2 && r.bernoulli(0.5), 1};
    }
    std::vector<double> c(4 * W, 0), lr, logA;
    for (const auto& st : sites) detail::add_site(c, xs[st.sequence], st.position, st.reverse, W, 1);
    Result best;
    best.sites = sites;
    best.F = F_of(sites);
    for (std::size_t sweep = 1; sweep <= o.sweeps; sweep++) {
      for (std::size_t z = 0; z < n; z++) {
        detail::add_site(c, xs[z], sites[z].position, sites[z].reverse, W, -1);
        detail::log_ratios(c, static_cast<double>(n - 1), B, p, lr);
        const auto& v = win.valid[z];
        logA.resize(v.size() * static_cast<std::size_t>(strands));
        double top = MOTIF_NEG_INF;
        for (std::size_t j = 0; j < v.size(); j++)
          for (int st = 0; st < strands; st++) top = std::max(top, logA[j * strands + st] = detail::window_score(lr, xs[z], v[j], st == 1, W));
        double total = 0;
        for (auto& a : logA) total += a = std::exp(a - top);
        double u = r.uniform() * total;
        std::size_t pick = logA.size() - 1;
        for (std::size_t k = 0; k < logA.size(); k++) {
          if (u < logA[k]) { pick = k; break; }
          u -= logA[k];
        }
        sites[z].position = v[pick / static_cast<std::size_t>(strands)];
        sites[z].reverse = pick % static_cast<std::size_t>(strands) == 1;
        detail::add_site(c, xs[z], sites[z].position, sites[z].reverse, W, 1);
      }
      if (o.shift_every && sweep % o.shift_every == 0) {
        // Shift every site by δ along its own strand; δ drawn in proportion to 2^F.
        std::vector<double> F(2 * max_shift + 1, MOTIF_NEG_INF);
        std::vector<std::vector<Site>> moved(F.size());
        for (std::size_t k = 0; k < F.size(); k++) {
          const long delta = static_cast<long>(k) - static_cast<long>(max_shift);
          std::vector<Site> m = sites;
          bool ok = true;
          for (auto& st : m) {
            const long at = static_cast<long>(st.position) + (st.reverse ? -delta : delta);
            const auto& v = win.valid[st.sequence];
            if (at < 0 || !std::binary_search(v.begin(), v.end(), static_cast<std::uint32_t>(at))) { ok = false; break; }
            st.position = static_cast<std::size_t>(at);
          }
          if (!ok) continue;
          F[k] = F_of(m);
          moved[k] = std::move(m);
        }
        const double top = *std::max_element(F.begin(), F.end());
        double total = 0;
        for (auto& f : F) total += f = f == MOTIF_NEG_INF ? 0 : std::exp2(f - top);
        double u = r.uniform() * total;
        for (std::size_t k = 0; k < F.size(); k++) {
          if (F[k] > 0 && u < F[k]) {
            sites = moved[k];
            break;
          }
          u -= F[k];
        }
        std::fill(c.begin(), c.end(), 0);
        for (const auto& st : sites) detail::add_site(c, xs[st.sequence], st.position, st.reverse, W, 1);
      }
      const double f = F_of(sites);
      if (f > best.F) { best.F = f; best.sites = sites; best.sweeps = sweep; }
    }
    results[restart] = std::move(best);
  };
  if (pool) calx::scheduler::parallel_for(*pool, results.size(), 1, [&](std::size_t b, std::size_t e, std::size_t) { for (std::size_t k = b; k < e; k++) run(k); });
  else
    for (std::size_t k = 0; k < results.size(); k++) run(k);
  std::size_t best = 0;
  for (std::size_t k = 1; k < results.size(); k++)
    if (results[k].F > results[best].F) best = k;
  Motif out;
  out.sites = results[best].sites;
  out.score = results[best].F;
  out.iterations = results[best].sweeps;
  std::vector<double> c(4 * W, 0);
  for (const auto& s : out.sites) detail::add_site(c, xs[s.sequence], s.position, s.reverse, W, 1);
  out.pwm = pwm_from_counts(c, B, p);
  return out;
}

enum class MotifModel { oops, zoops };

struct MemeOptions {
  std::size_t width = 10;
  MotifModel model = MotifModel::oops;
  std::size_t iterations = 500;
  double tolerance = 1e-6;   // stop when the log posterior rises by less (nats)
  double pseudocount = 1;    // β, spread by the background
  bool both_strands = false;
  double start_probability = 0.5;  // X, the probability of a start word's letter
  std::size_t max_starts = 0;      // 0: every distinct word; otherwise every k-th of them in order
  double gamma = 0.5;              // ZOOPS: the starting fraction of sequences with a site
};

namespace detail {

struct MemeState {
  std::vector<double> lr;  // log(p_{k,b}/p_0(b)), W × 4
  double gamma = 1;
};

// One E-step: the log-likelihood ratio (nats) and, if `counts`, the expected counts and Σ Z.
inline double meme_e_step(const std::vector<Codes>& xs, const Windows& win, const MemeState& st, MotifModel model, std::vector<double>* counts, double* z_sum,
                          std::vector<std::vector<double>>* z_out = nullptr) {
  const std::size_t W = win.width;
  const int strands = win.both ? 2 : 1;
  double ll = 0;
  std::vector<double> lz;
  for (std::size_t s = 0; s < xs.size(); s++) {
    const auto& v = win.valid[s];
    lz.resize(v.size() * static_cast<std::size_t>(strands));
    double top = MOTIF_NEG_INF;
    for (std::size_t j = 0; j < v.size(); j++)
      for (int t = 0; t < strands; t++) top = std::max(top, lz[j * strands + t] = window_score(st.lr, xs[s], v[j], t == 1, W));
    double sum = 0;
    for (double a : lz) sum += std::exp(a - top);
    const double m = static_cast<double>(lz.size());
    // log Σ_j R_ij = top + log sum.
    const double lsum = top + std::log(sum);
    double scale;  // Z_ij = exp(lz − top) · scale
    if (model == MotifModel::oops) {
      ll += lsum - std::log(m);
      scale = 1 / sum;
    } else {
      const double lam = st.gamma / m;
      // log((1 − γ) + λ Σ R), stably.
      const double a = std::log1p(-st.gamma), b = std::log(lam) + lsum;
      const double hi = std::max(a, b), l = hi + std::log(std::exp(a - hi) + std::exp(b - hi));
      ll += l;
      scale = std::exp(std::log(lam) + top - l);
    }
    if (z_out) (*z_out)[s].assign(lz.size(), 0);
    if (counts || z_out)
      for (std::size_t j = 0; j < v.size(); j++)
        for (int t = 0; t < strands; t++) {
          const double z = std::exp(lz[j * strands + t] - top) * scale;
          if (z_out) (*z_out)[s][j * strands + t] = z;
          if (counts && z > 1e-12) add_site(*counts, xs[s], v[j], t == 1, W, z);
          if (z_sum) *z_sum += z;
        }
  }
  return ll;
}

inline void meme_m_step(const std::vector<double>& counts, double z_sum, std::size_t sequences, const Background& p, double beta, MotifModel model, MemeState& st) {
  const std::size_t W = counts.size() / 4;
  st.lr.resize(counts.size());
  for (std::size_t i = 0; i < W; i++) {
    double n = beta;
    for (int b = 0; b < 4; b++) n += counts[4 * i + b];
    for (int b = 0; b < 4; b++) st.lr[4 * i + b] = std::log((counts[4 * i + b] + beta * p[b]) / n / p[b]);
  }
  if (model == MotifModel::zoops) st.gamma = std::clamp(z_sum / static_cast<double>(sequences), 1e-6, 1 - 1e-6);
}

}  // namespace detail

// The log-likelihood ratio of sequences under a motif model against the background (the PWM's
// own), bits: Σ_i log((1/m_i) Σ_j R_ij) for OOPS, Σ_i log((1 − γ) + (γ/m_i) Σ_j R_ij) for ZOOPS,
// R_ij the odds of window j; the quantity MEME's EM climbs.
inline double motif_log_likelihood(const std::vector<Codes>& xs, const Pwm& m, MotifModel model, double gamma = 1, bool both_strands = false) {
  const detail::Windows win = detail::windows(xs, m.width, both_strands);
  detail::MemeState st;
  st.gamma = gamma;
  st.lr.resize(4 * m.width);
  for (std::size_t i = 0; i < 4 * m.width; i++) st.lr[i] = std::log(m.probability[i] / m.background[i % 4]);
  return detail::meme_e_step(xs, win, st, model, nullptr, nullptr) / std::numbers::ln2;
}

inline Motif meme_motif(const std::vector<Codes>& xs, const MemeOptions& o = {}, calx::scheduler::ThreadPool* pool = nullptr) {
  const std::size_t n = xs.size(), W = o.width;
  if (n < 1 || W == 0 || W > 32) throw std::invalid_argument("genomics::annotation::meme_motif: sequences and a width from 1 to 32");
  const detail::Windows win = detail::windows(xs, W, o.both_strands);
  for (const auto& v : win.valid)
    if (v.empty()) throw std::invalid_argument("genomics::annotation::meme_motif: every sequence needs a window of the width without N");
  const Background p = base_frequencies(xs);
  // The distinct words of width W, as integers (two bits a base), in increasing order.
  std::vector<std::uint64_t> words;
  for (std::size_t s = 0; s < n; s++)
    for (std::uint32_t j : win.valid[s]) {
      std::uint64_t w = 0;
      for (std::size_t i = 0; i < W; i++) w = (w << 2) | xs[s][j + i];
      words.push_back(w);
    }
  std::sort(words.begin(), words.end());
  words.erase(std::unique(words.begin(), words.end()), words.end());
  if (o.max_starts && words.size() > o.max_starts) {
    const std::size_t step = (words.size() + o.max_starts - 1) / o.max_starts;
    std::vector<std::uint64_t> kept;
    for (std::size_t k = 0; k < words.size(); k += step) kept.push_back(words[k]);
    words.swap(kept);
  }
  const auto start_state = [&](std::uint64_t w) {
    detail::MemeState st;
    st.gamma = o.model == MotifModel::zoops ? o.gamma : 1;
    st.lr.resize(4 * W);
    for (std::size_t i = 0; i < W; i++) {
      const Code letter = static_cast<Code>((w >> (2 * (W - 1 - i))) & 3);
      for (int b = 0; b < 4; b++) st.lr[4 * i + b] = std::log((b == letter ? o.start_probability : (1 - o.start_probability) / 3) / p[b]);
    }
    return st;
  };
  // One iteration from each start: the likelihood after it.
  std::vector<double> after(words.size(), MOTIF_NEG_INF);
  const auto try_start = [&](std::size_t k) {
    detail::MemeState st = start_state(words[k]);
    std::vector<double> c(4 * W, 0);
    double zs = 0;
    detail::meme_e_step(xs, win, st, o.model, &c, &zs);
    detail::meme_m_step(c, zs, n, p, o.pseudocount, o.model, st);
    after[k] = detail::meme_e_step(xs, win, st, o.model, nullptr, nullptr);
  };
  if (pool) calx::scheduler::parallel_for(*pool, words.size(), 16, [&](std::size_t b, std::size_t e, std::size_t) { for (std::size_t k = b; k < e; k++) try_start(k); });
  else
    for (std::size_t k = 0; k < words.size(); k++) try_start(k);
  std::size_t best = 0;
  for (std::size_t k = 1; k < words.size(); k++)
    if (after[k] > after[best]) best = k;
  // EM to convergence from the best start.
  detail::MemeState st = start_state(words[best]);
  Motif out;
  double last = MOTIF_NEG_INF;
  std::size_t it = 0;
  std::vector<double> c(4 * W);
  for (; it < o.iterations; it++) {
    std::fill(c.begin(), c.end(), 0);
    double zs = 0;
    const double ll = detail::meme_e_step(xs, win, st, o.model, &c, &zs);
    // The pseudocounts are a Dirichlet prior (β q_b + 1): EM climbs the log posterior.
    double prior = 0;
    for (std::size_t i = 0; i < 4 * W; i++) prior += o.pseudocount * p[i % 4] * (st.lr[i] + std::log(p[i % 4]));
    out.trace.push_back((ll + prior) / std::numbers::ln2);
    if (it > 0 && ll + prior - last < o.tolerance) break;
    last = ll + prior;
    detail::meme_m_step(c, zs, n, p, o.pseudocount, o.model, st);
  }
  out.iterations = it;
  // The sites by the posteriors of the converged parameters; the matrix and γ of one more M-step,
  // which at convergence differ from them by less than the tolerance.
  std::vector<std::vector<double>> Z(n);
  std::fill(c.begin(), c.end(), 0);
  double zs = 0;
  detail::meme_e_step(xs, win, st, o.model, &c, &zs, &Z);
  detail::meme_m_step(c, zs, n, p, o.pseudocount, o.model, st);
  out.gamma = st.gamma;
  const int strands = o.both_strands ? 2 : 1;
  for (std::size_t s = 0; s < n; s++) {
    const auto& z = Z[s];
    const std::size_t k = static_cast<std::size_t>(std::max_element(z.begin(), z.end()) - z.begin());
    if (o.model == MotifModel::zoops && z[k] <= 0.5) continue;
    out.sites.push_back({s, win.valid[s][k / static_cast<std::size_t>(strands)], k % static_cast<std::size_t>(strands) == 1, z[k]});
  }
  out.pwm = pwm_from_counts(c, o.pseudocount, p);
  out.score = motif_log_likelihood(xs, out.pwm, o.model, out.gamma, o.both_strands);
  return out;
}

}  // namespace genomics::annotation
