// Promoters of the bacterial σ⁷⁰ holoenzyme: a model of two boxes and the spacer between them,
// learnt from transcription start sites by expectation maximisation, and scanned along a genome on
// both strands.
//
// The biology. RNA polymerase with σ⁷⁰ (σᴰ, RpoD) recognises two hexamers upstream of the start of
// transcription (+1): the −35 box, consensus TTGACA, and the −10 box, consensus TATAAT, separated by
// a spacer of 15 to 19 bases, 17 the commonest and most active (Hawley and McClure 1983, Nucleic
// Acids Research 11, 2237–2255; Harley and Reynolds 1987, Nucleic Acids Research 15, 2343–2361, the
// consensus from 263 promoters); the −10 box ends about 5 to 9 bases before +1. An extended −10
// (TGn just upstream of the −10 box) and an AT-rich UP element upstream of the −35 box (Ross et al.
// 1993, Science 262, 1407–1413) add to some promoters' strength; the model below has the two boxes
// and the spacer only.
//
// The model. A window x of the region −U … +D about a start site (U = 60 upstream, D = 10 from +1)
// either holds a promoter, with probability γ, or is background, independent bases of composition q.
// A promoter places its −10 box (width w₁₀ = 6) at offset r from +1 (r ∈ [r_min, r_max]) and its −35
// box (width w₃₅ = 6) a spacer s (s ∈ [15, 19]) before it, with probabilities π_r and π_s; the boxes'
// bases are drawn from position probability matrices P₁₀, P₃₅, the rest from q. The likelihood
// ratio of configuration z = (r, s) against the background is
//   R(x, z) = Π_k P₁₀[k, x_{a+k}]/q(x_{a+k}) · Π_k P₃₅[k, x_{b+k}]/q(x_{b+k}),  a = U + r, b = a − s − w₃₅,
// and that of the window Λ(x) = (1 − γ) + γ Σ_z π_r π_s R(x, z) (a zero-or-one-occurrence model with
// two boxes and a variable gap, as BioProspector's two-block motifs; Liu, Brutlag and Liu 2001, Pac.
// Symp. Biocomput. 6, 127–138).
//
// Learning, by expectation maximisation (Dempster, Laird and Rubin 1977, JRSS B 39, 1–38): the E-step
// gives each window the posterior of each configuration, γ π_r π_s R(x, z)/Λ(x), and of none; the
// M-step re-estimates P₁₀ and P₃₅ from the posterior-weighted bases (with pseudocounts β q_b), π_r,
// π_s (with a pseudocount of one each) and γ. With the pseudocounts a Dirichlet prior, each iteration
// raises the log posterior, which the trace records (Dempster et al., Theorem 1). Starting point: the
// −10 box found by MEME's EM (annotation::meme_motif, zero-or-one occurrence) in the windows'
// −10 region, then the −35 box found the same way in the region the spacers allow upstream of each
// window's best −10 box, each moved by a few bases either way (MEME's sites may lie a few bases
// off their box), EM run from every pair and the most probable model kept: the model learns its
// boxes from the data, the consensus is not given.
//
// Scanning. Along a sequence, on both strands, each −10 position is scored by the log-odds, in
// bits, of the best −35 box and spacer with it: S = log₂ R₁₀ + max_s (log₂ R₃₅ + log₂ π_s/max π);
// hits above a threshold are kept. A threshold for a rate of false hits is found by scanning random
// sequence of the background's composition (Monte Carlo on a calx stream). The promoter score of a
// window about a known start site is log₂ Σ_z π_r π_s R(x, z), its likelihood ratio of holding a
// promoter.
//
// Limits: the boxes have fixed widths; the extended −10, the UP element and the discriminator are
// not modelled, nor other σ factors (σ³², σ⁵⁴, σ³⁸ promoters are background here); the background
// is of order zero. Complexity: an EM iteration O(N · |r| · |s| · (w₁₀ + w₃₅)); a scan O(L · |s| · w)
// per strand.

#pragma once

#include <genomics/annotation/motifs.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/fasta.hpp>

#include <calx/random/stream.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::regulation {

using sequence::Code;
using sequence::Codes;

inline constexpr std::uint32_t PROMOTER_NULL = 0x504e554c;  // "PNUL": random sequence for thresholds

// ---- transcription start sites -----------------------------------------------------------------------
struct StartSite {
  std::size_t position = 0;  // 0-based, the base of +1
  bool reverse = false;
  double strength = 0;       // the relative read score (RRS) of Cappable-seq
  std::size_t line = 0;
};

// Ettwiller et al. (2016)'s sites: a GTF line per site, start = end = +1 (1-based), the RRS in the
// score column, the strand in the seventh.
inline std::vector<StartSite> read_start_sites(std::string_view text, std::vector<sequence::Problem>* problems = nullptr) {
  std::vector<StartSite> out;
  sequence::detail::for_each_line(text, [&](std::string_view line, std::size_t number) {
    if (line.empty() || line[0] == '#') return;
    std::vector<std::string_view> col;
    std::size_t from = 0;
    for (;;) {
      const auto tab = line.find('\t', from);
      col.push_back(line.substr(from, tab == std::string_view::npos ? std::string_view::npos : tab - from));
      if (tab == std::string_view::npos) break;
      from = tab + 1;
    }
    if (col.size() < 7) { if (problems) problems->push_back({number, "a site line has fewer than seven columns"}); return; }
    StartSite s;
    const std::string pos(col[3]), score(col[5]);
    char* e1 = nullptr;
    char* e2 = nullptr;
    const unsigned long long p = std::strtoull(pos.c_str(), &e1, 10);
    s.strength = std::strtod(score.c_str(), &e2);
    if (pos.empty() || *e1 || p < 1 || score.empty() || *e2 || (col[6] != "+" && col[6] != "-")) {
      if (problems) problems->push_back({number, "a site's position, score or strand cannot be read"});
      return;
    }
    s.position = static_cast<std::size_t>(p - 1);
    s.reverse = col[6] == "-";
    s.line = number;
    out.push_back(s);
  });
  return out;
}

// The window −upstream … +downstream about a site, read on the site's strand (+1 at index
// `upstream`); empty when it leaves the genome (a circular genome is not wrapped here).
inline Codes site_window(const Codes& genome, const StartSite& s, std::size_t upstream = 60, std::size_t downstream = 10) {
  if (!s.reverse) {
    if (s.position < upstream || s.position + downstream > genome.size()) return {};
    return Codes(genome.begin() + static_cast<std::ptrdiff_t>(s.position - upstream), genome.begin() + static_cast<std::ptrdiff_t>(s.position + downstream));
  }
  if (s.position + upstream >= genome.size() || s.position + 1 < downstream) return {};
  const Codes forward(genome.begin() + static_cast<std::ptrdiff_t>(s.position + 1 - downstream), genome.begin() + static_cast<std::ptrdiff_t>(s.position + upstream + 1));
  return sequence::reverse_complement(forward);
}

// ---- the model ---------------------------------------------------------------------------------------------
struct Sigma70 {
  std::size_t upstream = 60;                // +1 at this index of a window
  int r_min = -17, r_max = -9;              // the −10 box's first base, relative to +1
  std::size_t s_min = 15, s_max = 19;       // spacers
  annotation::Pwm minus10, minus35;         // probabilities and log₂-odds against the background
  std::vector<double> pi_r, pi_s;           // over r_min … r_max and s_min … s_max
  double gamma = 0.5;                       // the fraction of windows with a promoter
  annotation::Background background = annotation::UNIFORM;
  std::vector<double> trace;                // the log posterior (bits) at each EM iteration
  std::size_t iterations = 0;

  std::size_t w10() const { return minus10.width; }
  std::size_t w35() const { return minus35.width; }
  std::size_t spacers() const { return s_max - s_min + 1; }
  std::size_t offsets() const { return static_cast<std::size_t>(r_max - r_min + 1); }
  std::size_t spacer_mode() const { return s_min + static_cast<std::size_t>(std::max_element(pi_s.begin(), pi_s.end()) - pi_s.begin()); }
  int offset_mode() const { return r_min + static_cast<int>(std::max_element(pi_r.begin(), pi_r.end()) - pi_r.begin()); }
};

namespace detail {

// log of P/q for each column and base, natural.
inline std::vector<double> log_ratios(const annotation::Pwm& m, const annotation::Background& q) {
  std::vector<double> lr(4 * m.width);
  for (std::size_t i = 0; i < m.width; i++)
    for (int b = 0; b < 4; b++) lr[4 * i + b] = std::log(m.probability[4 * i + b] / q[b]);
  return lr;
}

inline double box(const std::vector<double>& lr, std::size_t width, const Codes& x, std::size_t at) {
  double s = 0;
  for (std::size_t k = 0; k < width; k++) {
    const Code c = x[at + k];
    if (c > 3) return -std::numeric_limits<double>::infinity();
    s += lr[4 * k + c];
  }
  return s;
}

inline std::string consensus(const annotation::Pwm& m) {
  std::string s;
  for (std::size_t i = 0; i < m.width; i++) {
    int best = 0;
    for (int b = 1; b < 4; b++)
      if (m.probability[4 * i + b] > m.probability[4 * i + best]) best = b;
    s += "ACGT"[best];
  }
  return s;
}

inline annotation::Pwm pwm_from(const std::vector<double>& counts, double beta, const annotation::Background& q) { return annotation::pwm_from_counts(counts, beta, q); }

}  // namespace detail

inline std::string consensus(const annotation::Pwm& m) { return detail::consensus(m); }

struct WindowPosterior {
  double log_ratio = 0;            // log Λ(x)/(1 − γ)… see promoter_log_odds; here log Σ_z π π R (nats)
  double present = 0;              // P(promoter | x)
  int best_r = 0;
  std::size_t best_s = 0;
  double best = 0;                 // the posterior of the best configuration
};

// The posterior of one window under the model (natural logarithms inside).
inline WindowPosterior window_posterior(const Sigma70& m, const std::vector<double>& lr10, const std::vector<double>& lr35, const Codes& x, std::vector<double>* z = nullptr) {
  const std::size_t R = m.offsets(), S = m.spacers();
  std::vector<double> lz(R * S, -std::numeric_limits<double>::infinity());
  double top = -std::numeric_limits<double>::infinity();
  for (std::size_t ri = 0; ri < R; ri++) {
    const long a = static_cast<long>(m.upstream) + m.r_min + static_cast<long>(ri);
    if (a < 0 || static_cast<std::size_t>(a) + m.w10() > x.size()) continue;
    const double s10 = detail::box(lr10, m.w10(), x, static_cast<std::size_t>(a));
    for (std::size_t si = 0; si < S; si++) {
      const long b = a - static_cast<long>(m.s_min + si) - static_cast<long>(m.w35());
      if (b < 0) continue;
      const double v = std::log(m.pi_r[ri]) + std::log(m.pi_s[si]) + s10 + detail::box(lr35, m.w35(), x, static_cast<std::size_t>(b));
      lz[ri * S + si] = v;
      top = std::max(top, v);
    }
  }
  WindowPosterior w;
  if (!std::isfinite(top)) { w.log_ratio = -std::numeric_limits<double>::infinity(); if (z) z->assign(R * S, 0); return w; }
  double sum = 0;
  for (double v : lz) sum += std::exp(v - top);
  w.log_ratio = top + std::log(sum);  // log Σ π π R
  // P(present) = γ e^{ℓ} / ((1 − γ) + γ e^{ℓ}).
  const double a = std::log(m.gamma) + w.log_ratio, b = std::log1p(-m.gamma), hi = std::max(a, b);
  const double lnorm = hi + std::log(std::exp(a - hi) + std::exp(b - hi));
  w.present = std::exp(a - lnorm);
  std::size_t arg = 0;
  for (std::size_t k = 0; k < lz.size(); k++)
    if (lz[k] > lz[arg]) arg = k;
  w.best_r = m.r_min + static_cast<int>(arg / S);
  w.best_s = m.s_min + arg % S;
  w.best = std::exp(lz[arg] - top) / sum * w.present;
  if (z) {
    z->resize(R * S);
    for (std::size_t k = 0; k < lz.size(); k++) (*z)[k] = std::exp(lz[k] - top) / sum * w.present;
  }
  return w;
}

// The log₂ likelihood ratio of a window holding a promoter, log₂ Σ_z π_r π_s R(x, z).
inline double promoter_score(const Sigma70& m, const Codes& window) {
  const auto lr10 = detail::log_ratios(m.minus10, m.background), lr35 = detail::log_ratios(m.minus35, m.background);
  return window_posterior(m, lr10, lr35, window).log_ratio / std::numbers::ln2;
}

struct Sigma70Options {
  std::size_t width10 = 6, width35 = 6;
  int r_min = -17, r_max = -9;
  std::size_t s_min = 15, s_max = 19;
  std::size_t upstream = 60;
  double pseudocount = 1;          // β, spread by the background
  std::size_t iterations = 200;
  double tolerance = 1e-6;         // nats of log posterior per window, to stop
  std::size_t meme_starts = 400;   // MEME's starting words, at most
  std::size_t shift = 2;           // the phase shifts of MEME's sites tried, ±
};

namespace detail {

// Counts of a box's bases at MEME's sites moved by δ within their regions (sites that leave a
// region are dropped).
inline std::vector<double> shifted_counts(const std::vector<Codes>& regions, const annotation::Motif& motif, long delta, std::size_t width) {
  std::vector<double> c(4 * width, 0);
  for (const auto& s : motif.sites) {
    const long at = static_cast<long>(s.position) + delta;
    const Codes& x = regions[s.sequence];
    if (at < 0 || static_cast<std::size_t>(at) + width > x.size()) continue;
    for (std::size_t k = 0; k < width; k++)
      if (x[static_cast<std::size_t>(at) + k] < 4) c[4 * k + x[static_cast<std::size_t>(at) + k]] += s.probability;
  }
  return c;
}

// EM from a model's starting boxes (uniform π, γ = ½), to convergence.
inline void em_sigma70(const std::vector<Codes>& windows, Sigma70& m, const Sigma70Options& o) {
  const auto& q = m.background;
  m.pi_r.assign(m.offsets(), 1.0 / static_cast<double>(m.offsets()));
  m.pi_s.assign(m.spacers(), 1.0 / static_cast<double>(m.spacers()));
  m.gamma = 0.5;
  m.trace.clear();
  const std::size_t R = m.offsets(), S = m.spacers(), W10 = m.w10(), W35 = m.w35();
  double last = -std::numeric_limits<double>::infinity();
  std::vector<double> z;
  for (std::size_t it = 0; it < o.iterations; it++) {
    const auto l10 = log_ratios(m.minus10, q), l35 = log_ratios(m.minus35, q);
    std::vector<double> c10(4 * W10, 0), c35(4 * W35, 0), nr(R, 0), ns(S, 0);
    double present = 0, loglik = 0;
    for (const auto& x : windows) {
      const auto w = window_posterior(m, l10, l35, x, &z);
      // log Λ = log((1 − γ) + γ e^{ℓ})
      const double a = std::log(m.gamma) + w.log_ratio, b = std::log1p(-m.gamma), hi = std::max(a, b);
      loglik += std::isfinite(hi) ? hi + std::log(std::exp(a - hi) + std::exp(b - hi)) : b;
      present += w.present;
      for (std::size_t ri = 0; ri < R; ri++)
        for (std::size_t si = 0; si < S; si++) {
          const double p = z[ri * S + si];
          if (p < 1e-12) continue;
          nr[ri] += p;
          ns[si] += p;
          const std::size_t at = m.upstream + static_cast<std::size_t>(m.r_min + static_cast<int>(ri));
          const std::size_t bt = at - (m.s_min + si) - W35;
          for (std::size_t k = 0; k < W10; k++)
            if (x[at + k] < 4) c10[4 * k + x[at + k]] += p;
          for (std::size_t k = 0; k < W35; k++)
            if (x[bt + k] < 4) c35[4 * k + x[bt + k]] += p;
        }
    }
    // The Dirichlet prior's term, so that the trace is the log posterior EM climbs.
    double prior = 0;
    for (std::size_t i = 0; i < 4 * W10; i++) prior += o.pseudocount * q[i % 4] * std::log(m.minus10.probability[i]);
    for (std::size_t i = 0; i < 4 * W35; i++) prior += o.pseudocount * q[i % 4] * std::log(m.minus35.probability[i]);
    for (double p : m.pi_r) prior += std::log(p);
    for (double p : m.pi_s) prior += std::log(p);
    const double objective = loglik + prior;
    m.trace.push_back(objective / std::numbers::ln2);
    m.iterations = it + 1;
    if (it > 0 && objective - last < o.tolerance * static_cast<double>(windows.size())) break;
    last = objective;
    m.minus10 = annotation::pwm_from_counts(c10, o.pseudocount, q);
    m.minus35 = annotation::pwm_from_counts(c35, o.pseudocount, q);
    double tr = 0, ts = 0;
    for (std::size_t ri = 0; ri < R; ri++) tr += nr[ri] + 1;
    for (std::size_t si = 0; si < S; si++) ts += ns[si] + 1;
    for (std::size_t ri = 0; ri < R; ri++) m.pi_r[ri] = (nr[ri] + 1) / tr;
    for (std::size_t si = 0; si < S; si++) m.pi_s[si] = (ns[si] + 1) / ts;
    m.gamma = std::clamp(present / static_cast<double>(windows.size()), 1e-6, 1 - 1e-6);
  }
}

}  // namespace detail

// Learns the model from windows about start sites (each `upstream` bases before +1 and some after).
// MEME's sites may be found a few bases off the box they belong to (a phase shift, the trap Lawrence
// et al.'s sampler shifts out of); so the starting boxes are MEME's sites moved by δ₁₀ and δ₃₅ in
// [−shift, shift], EM runs from each pair, and the model of largest log posterior is kept.
inline Sigma70 learn_sigma70(const std::vector<Codes>& windows, const Sigma70Options& o = {}, calx::scheduler::ThreadPool* pool = nullptr) {
  if (windows.size() < 10) throw std::invalid_argument("genomics::regulation::learn_sigma70: at least ten windows");
  Sigma70 base;
  base.upstream = o.upstream;
  base.r_min = o.r_min;
  base.r_max = o.r_max;
  base.s_min = o.s_min;
  base.s_max = o.s_max;
  base.background = annotation::base_frequencies(windows);
  const auto& q = base.background;
  // The −10 region of each window; MEME's −10 box in it.
  std::vector<Codes> region10;
  const long a10 = static_cast<long>(o.upstream) + o.r_min, b10 = static_cast<long>(o.upstream) + o.r_max + static_cast<long>(o.width10);
  for (const auto& x : windows) {
    if (a10 < static_cast<long>(o.s_max + o.width35) || static_cast<std::size_t>(b10) > x.size()) throw std::invalid_argument("genomics::regulation::learn_sigma70: windows shorter than the model");
    region10.emplace_back(x.begin() + a10, x.begin() + b10);
  }
  annotation::MemeOptions mo;
  mo.width = o.width10;
  mo.model = annotation::MotifModel::zoops;
  mo.max_starts = o.meme_starts;
  const auto m10 = annotation::meme_motif(region10, mo, pool);
  Sigma70 best;
  double best_value = -std::numeric_limits<double>::infinity();
  const long shift = static_cast<long>(o.shift);
  for (long d10 = -shift; d10 <= shift; d10++) {
    Sigma70 start = base;
    start.minus10 = annotation::pwm_from_counts(detail::shifted_counts(region10, m10, d10, o.width10), o.pseudocount, q);
    // The −35 region each window's best −10 box allows; MEME's −35 box in it.
    const auto lr10 = detail::log_ratios(start.minus10, q);
    std::vector<Codes> region35;
    for (const auto& x : windows) {
      long best_a = a10;
      double top = -std::numeric_limits<double>::infinity();
      for (int r = o.r_min; r <= o.r_max; r++) {
        const long a = static_cast<long>(o.upstream) + r;
        const double v = detail::box(lr10, o.width10, x, static_cast<std::size_t>(a));
        if (v > top) { top = v; best_a = a; }
      }
      const long lo = best_a - static_cast<long>(o.s_max + o.width35), hi = best_a - static_cast<long>(o.s_min);
      region35.emplace_back(x.begin() + lo, x.begin() + hi);
    }
    annotation::MemeOptions m35o = mo;
    m35o.width = o.width35;
    const auto m35 = annotation::meme_motif(region35, m35o, pool);
    for (long d35 = -shift; d35 <= shift; d35++) {
      Sigma70 m = start;
      m.minus35 = annotation::pwm_from_counts(detail::shifted_counts(region35, m35, d35, o.width35), o.pseudocount, q);
      detail::em_sigma70(windows, m, o);
      if (m.trace.back() > best_value) {
        best_value = m.trace.back();
        best = std::move(m);
      }
    }
  }
  return best;
}

// ---- scanning ------------------------------------------------------------------------------------------------
struct PromoterHit {
  std::size_t minus10 = 0, minus35 = 0;  // the boxes' first bases, forward-strand coordinates
  std::size_t spacer = 0;
  bool reverse = false;
  double score = 0;                       // bits
  double score10 = 0, score35 = 0;        // bits, the boxes alone
  std::size_t start = 0;                  // the start site the most likely offset predicts (+1), forward strand
};

// Every −10 position on both strands scoring at least `threshold` bits, with its best −35 box and
// spacer, by position then strand.
inline std::vector<PromoterHit> scan_promoters(const Sigma70& m, const Codes& x, double threshold) {
  std::vector<PromoterHit> out;
  const auto l10 = detail::log_ratios(m.minus10, m.background), l35 = detail::log_ratios(m.minus35, m.background);
  const double top_s = *std::max_element(m.pi_s.begin(), m.pi_s.end());
  const std::size_t W10 = m.w10(), W35 = m.w35(), L = x.size();
  const int mode = m.offset_mode();
  for (int strand = 0; strand < 2; strand++) {
    const Codes y = strand ? sequence::reverse_complement(x) : Codes{};
    const Codes& s = strand ? y : x;
    for (std::size_t a = m.s_min + W35; a + W10 <= L; a++) {
      const double s10 = detail::box(l10, W10, s, a);
      if (!std::isfinite(s10)) continue;
      double best = -std::numeric_limits<double>::infinity(), best35 = 0;
      std::size_t best_sp = m.s_min;
      for (std::size_t sp = m.s_min; sp <= m.s_max; sp++) {
        if (a < sp + W35) continue;
        const double v35 = detail::box(l35, W35, s, a - sp - W35);
        const double v = v35 + std::log(m.pi_s[sp - m.s_min] / top_s);
        if (v > best) { best = v; best35 = v35; best_sp = sp; }
      }
      const double score = (s10 + best) / std::numbers::ln2;
      if (!(score >= threshold)) continue;
      PromoterHit h;
      h.reverse = strand == 1;
      h.spacer = best_sp;
      h.score = score;
      h.score10 = s10 / std::numbers::ln2;
      h.score35 = best35 / std::numbers::ln2;
      const std::size_t b = a - best_sp - W35;
      const long tss = static_cast<long>(a) - mode;  // +1 lies −r bases after the −10 box's first base
      if (!h.reverse) {
        h.minus10 = a;
        h.minus35 = b;
        h.start = static_cast<std::size_t>(std::clamp<long>(tss, 0, static_cast<long>(L) - 1));
      } else {
        // Position p on the reverse strand is L − 1 − p forward; a box [p, p + w) is [L − p − w, L − p).
        h.minus10 = L - a - W10;
        h.minus35 = L - b - W35;
        h.start = L - 1 - static_cast<std::size_t>(std::clamp<long>(tss, 0, static_cast<long>(L) - 1));
      }
      out.push_back(h);
    }
  }
  std::sort(out.begin(), out.end(), [](const PromoterHit& p, const PromoterHit& q) { return p.minus10 != q.minus10 ? p.minus10 < q.minus10 : p.reverse < q.reverse; });
  return out;
}

// The score above which random sequence of the background's composition gives `rate` hits per
// position and strand, from `length` simulated bases.
inline double promoter_threshold(const Sigma70& m, double rate, std::size_t length = 2'000'000, std::uint64_t seed = 1) {
  calx::random::Stream r(seed, PROMOTER_NULL, 0, 0);
  Codes x(length);
  const auto& q = m.background;
  for (auto& c : x) {
    const double u = r.uniform();
    c = u < q[0] ? 0 : u < q[0] + q[1] ? 1 : u < q[0] + q[1] + q[2] ? 2 : 3;
  }
  auto hits = scan_promoters(m, x, -1e9);
  std::vector<double> scores;
  scores.reserve(hits.size());
  for (const auto& h : hits) scores.push_back(h.score);
  const std::size_t k = std::min(scores.size() - 1, static_cast<std::size_t>(rate * 2.0 * static_cast<double>(length)));
  std::nth_element(scores.begin(), scores.begin() + static_cast<std::ptrdiff_t>(k), scores.end(), std::greater<double>());
  return scores[k];
}

}  // namespace genomics::regulation
