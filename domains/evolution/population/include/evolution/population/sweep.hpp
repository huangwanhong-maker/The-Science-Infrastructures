// A selective sweep forward in time: a beneficial allele arising on one chromosome of a
// Wright–Fisher population whose standing variation is a coalescent sample, carried up in
// frequency by selection together with the chromosome around it, recombination breaking that
// chromosome up; and a sample of the population at the end, as haplotypes.
//
// The model, stated so that its every step can be checked:
//
// - The population is 2N chromosomes (N diploid individuals), each the interval [0, 1) of the
//   coalescent. Its standing variation at the start is a coalescent sample of all 2N
//   chromosomes, with θ = 4Nμ and ρ = 4Nr for the whole chromosome (coalescent.hpp): the
//   neutral population's state in the coalescent's limit, an approximation to Wright–Fisher's
//   stationary state that is good when N is large.
// - At generation 0 the beneficial allele arises at `position` on one chromosome drawn
//   uniformly: a hard sweep from a single new copy.
// - Each chromosome of the next generation is formed from two chromosomes of the present one,
//   each drawn independently with probability proportional to its fitness, 1 + s with the
//   beneficial allele and 1 without (genic selection with the chromosome as the unit). With
//   k ~ Poisson(R) crossovers at uniform positions, R = ρ/(4N) the length of the chromosome in
//   Morgans, the child copies the first chromosome up to the first crossover, the second up to
//   the next, and so on. Each child carries the beneficial allele with probability
//   ψ(p) = p(1 + s)/(1 + sp), independently, so the frequency is Binomial(2N, ψ(p))/(2N): the
//   Wright–Fisher chain of wright_fisher.hpp with the haploid ψ. Its diffusion has drift
//   s x(1 − x) and variance x(1 − x)/(2N), so the probability of fixing from one copy is
//   Kimura's u(1/(2N)) = (1 − e^{−2s})/(1 − e^{−4Ns}) (kimura_fixation(N, s, 1/(2N))).
// - No new mutation arises during the sweep. The deterministic recursion p′ = ψ(p) multiplies
//   the odds p/(1 − p) by 1 + s each generation, so it goes from 1/(2N) to 1 − 1/(2N) in
//   2 ln(2N − 1)/ln(1 + s) generations, about (2/s) ln(2N): short against the 2N generations of a
//   unit of coalescent time when 2Ns is large.
// - An attempt whose allele is lost starts again from the same standing variation with new
//   draws, until one reaches the stop frequency (1 by default: fixation; less for a partial
//   sweep). The attempts are counted, so 1/attempts estimates the probability of reaching it.
// - At the end n chromosomes are sampled without replacement; the sites polymorphic in the
//   sample are kept, the selected site among them while it segregates, 0 the ancestral allele
//   and 1 the derived (the beneficial one at the selected site).
//
// Draws: the standing variation is the coalescent's replicate (seed, replicate); attempt a
// uses the stream (seed, SWEEP, replicate, a), the sample (seed, SWEEP_SAMPLE, replicate), so
// a replicate is the same on any thread.
//
// Complexity: O(2N W) per generation, W = ⌈(S + 1)/64⌉ words per chromosome for S standing
// sites (each child is a few masked copies of words), O(G 2N W) for G generations; plus the
// coalescent of 2N chromosomes.
//
// Limits: one population of constant size; no dominance (the chromosome is the unit of
// selection) and no new mutation; standing variation from the coalescent's limit, not from a
// forward burn-in. Soft sweeps from standing variation and the structured coalescent of a
// sweep (Kaplan, Hudson and Langley 1989) are not written.

#pragma once

#include <evolution/population/coalescent.hpp>
#include <evolution/population/genotypes.hpp>

#include <calx/random/stream.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace evolution::population {

inline constexpr std::uint32_t SWEEP = 0x53574550, SWEEP_SAMPLE = 0x53574153;

struct SweepSpec {
  std::size_t individuals = 500;  // N: 2N chromosomes
  std::size_t samples = 100;      // n chromosomes sampled at the end
  double theta = 50, rho = 100;   // 4Nμ and 4Nr for the whole chromosome [0, 1)
  double s = 0.05;                // a chromosome with the beneficial allele is a parent with weight 1 + s
  double position = 0.5;          // of the selected site
  double stop = 1;                // the population frequency at which the sweep is sampled (1: fixation)
  std::size_t max_attempts = 1000000;
  std::size_t max_generations = 1000000;  // of one attempt
};

struct Sweep {
  static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();
  Haplotypes haplotypes;            // the sample at its polymorphic sites
  std::size_t selected = npos;      // the selected site among them, npos when it does not segregate in the sample
  double frequency = 0;             // of the beneficial allele in the population at the end
  std::vector<double> trajectory;   // its frequency at generations 0, 1, … of the attempt that succeeded
  std::size_t generations = 0, attempts = 0, standing_sites = 0;
};

namespace detail {

// Copies bits [lo, hi) of src into dst.
inline void copy_bits(std::uint64_t* dst, const std::uint64_t* src, std::size_t lo, std::size_t hi) {
  if (lo >= hi) return;
  const std::size_t wl = lo / 64, wh = (hi - 1) / 64;
  const std::uint64_t ml = ~0ull << (lo % 64), mh = hi % 64 ? (1ull << (hi % 64)) - 1 : ~0ull;
  if (wl == wh) {
    const std::uint64_t m = ml & mh;
    dst[wl] = (dst[wl] & ~m) | (src[wl] & m);
    return;
  }
  dst[wl] = (dst[wl] & ~ml) | (src[wl] & ml);
  for (std::size_t w = wl + 1; w < wh; w++) dst[w] = src[w];
  dst[wh] = (dst[wh] & ~mh) | (src[wh] & mh);
}

inline bool bit(const std::uint64_t* x, std::size_t b) { return (x[b / 64] >> (b % 64)) & 1u; }

}  // namespace detail

inline Sweep sweep(const SweepSpec& spec, std::uint64_t seed, std::uint32_t replicate = 0) {
  if (spec.individuals == 0 || spec.samples == 0 || spec.samples > 2 * spec.individuals)
    throw std::invalid_argument("evolution::population::sweep: a sample of 1 … 2N chromosomes");
  if (!(spec.s > -1) || !(spec.stop > 0 && spec.stop <= 1) || !(spec.position >= 0 && spec.position < 1) || !(spec.rho >= 0) || !(spec.theta >= 0))
    throw std::invalid_argument("evolution::population::sweep: s > −1, a stop frequency in (0, 1], a position on [0, 1), θ and ρ >= 0");
  const std::size_t M = 2 * spec.individuals;

  // The standing variation, and the selected site inserted among its sites by position.
  CoalescentSpec cs;
  cs.samples = M;
  cs.theta = spec.theta;
  cs.rho = spec.rho;
  const Haplotypes standing = haplotypes(coalescent(cs, seed, replicate));
  const std::size_t S = standing.sites + 1;
  const std::size_t at = static_cast<std::size_t>(std::lower_bound(standing.positions.begin(), standing.positions.end(), spec.position) - standing.positions.begin());
  std::vector<double> positions(S);
  for (std::size_t b = 0; b < S; b++) positions[b] = b < at ? standing.positions[b] : b == at ? spec.position : standing.positions[b - 1];
  const std::size_t W = (S + 63) / 64;
  std::vector<std::uint64_t> initial(M * W, 0);
  for (std::size_t b = 0; b < S; b++) {
    if (b == at) continue;
    const auto row = standing.site(b < at ? b : b - 1);
    for (std::size_t i = 0; i < M; i++)
      if (row[i] == 1) initial[i * W + b / 64] |= 1ull << (b % 64);
  }

  Sweep out;
  out.standing_sites = standing.sites;
  const double R = spec.rho / (4.0 * static_cast<double>(spec.individuals));
  const std::size_t stop_count = static_cast<std::size_t>(std::ceil(spec.stop * static_cast<double>(M) - 1e-9));
  std::vector<std::uint64_t> now(M * W), next(M * W);
  std::vector<std::uint32_t> carriers, others;
  std::vector<double> breaks;
  bool done = false;
  for (std::size_t attempt = 0; attempt < spec.max_attempts && !done; attempt++) {
    calx::random::Stream r(seed, SWEEP, replicate, static_cast<std::uint32_t>(attempt));
    now = initial;
    const std::size_t founder = static_cast<std::size_t>(r.below(M));
    now[founder * W + at / 64] |= 1ull << (at % 64);
    std::size_t count = 1;
    out.trajectory.assign(1, 1.0 / static_cast<double>(M));
    out.attempts = attempt + 1;
    for (std::size_t gen = 0; gen < spec.max_generations; gen++) {
      if (count >= stop_count) {
        done = true;
        out.generations = gen;
        break;
      }
      if (count == 0) break;
      carriers.clear();
      others.clear();
      for (std::size_t i = 0; i < M; i++) (detail::bit(&now[i * W], at) ? carriers : others).push_back(static_cast<std::uint32_t>(i));
      const double wb = static_cast<double>(carriers.size()) * (1 + spec.s), ww = static_cast<double>(others.size());
      const double pb = wb / (wb + ww);
      const auto parent = [&]() -> std::size_t {
        if (others.empty() || (!carriers.empty() && r.uniform() < pb)) return carriers[r.below(carriers.size())];
        return others[r.below(others.size())];
      };
      count = 0;
      for (std::size_t c = 0; c < M; c++) {
        const std::size_t p1 = parent(), p2 = parent();
        std::uint64_t* child = &next[c * W];
        const std::int64_t k = R > 0 ? r.poisson(R) : 0;
        if (k == 0) {
          std::copy(&now[p1 * W], &now[p1 * W] + W, child);
        } else {
          breaks.resize(static_cast<std::size_t>(k));
          for (double& y : breaks) y = r.uniform();
          std::sort(breaks.begin(), breaks.end());
          const std::uint64_t* src[2] = {&now[p1 * W], &now[p2 * W]};
          std::size_t start = 0, which = 0;
          for (const double y : breaks) {
            const std::size_t b = static_cast<std::size_t>(std::lower_bound(positions.begin(), positions.end(), y) - positions.begin());
            detail::copy_bits(child, src[which], start, b);
            start = std::max(start, b);
            which ^= 1;
          }
          detail::copy_bits(child, src[which], start, S);
        }
        count += detail::bit(child, at);
      }
      now.swap(next);
      out.trajectory.push_back(static_cast<double>(count) / static_cast<double>(M));
    }
  }
  if (!done) throw std::runtime_error("evolution::population::sweep: no attempt reached the stop frequency within the limits");
  out.frequency = out.trajectory.back();

  // The sample, and its polymorphic sites.
  calx::random::Stream r(seed, SWEEP_SAMPLE, replicate, 0);
  std::vector<std::uint32_t> order(M);
  for (std::size_t i = 0; i < M; i++) order[i] = static_cast<std::uint32_t>(i);
  for (std::size_t i = 0; i < spec.samples; i++) std::swap(order[i], order[i + r.below(M - i)]);
  Haplotypes& h = out.haplotypes;
  h.sequences = spec.samples;
  std::vector<std::int8_t> row(spec.samples);
  for (std::size_t b = 0; b < S; b++) {
    std::size_t ones = 0;
    for (std::size_t i = 0; i < spec.samples; i++) ones += row[i] = static_cast<std::int8_t>(detail::bit(&now[order[i] * W], b));
    if (ones == 0 || ones == spec.samples) continue;
    if (b == at) out.selected = h.sites;
    h.alleles.insert(h.alleles.end(), row.begin(), row.end());
    h.positions.push_back(positions[b]);
    h.sites++;
  }
  return out;
}

}  // namespace evolution::population
