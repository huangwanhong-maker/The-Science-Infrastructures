// Linkage disequilibrium between two biallelic sites: D, D′ and r², from phased haplotypes or
// from unphased genotypes, and its decay with distance.
//
// With A and B the alternate alleles of sites 1 and 2, p_A and p_B their frequencies and p_AB
// that of the haplotype carrying both:
//
//   D = p_AB − p_A p_B                                (Lewontin and Kojima 1960)
//   D′ = D / D_max, D_max = min(p_A q_B, q_A p_B) for D > 0, min(p_A p_B, q_A q_B) for D < 0
//                                                     (Lewontin 1964, Genetics 49, 49–67)
//   r² = D² / (p_A q_A p_B q_B)                       (Hill and Robertson 1968, TAG 38, 226–231)
//
// D′ keeps the sign of D, so it lies in [−1, 1]; r² is the squared correlation of the two
// alleles on a haplotype. Both are undefined (NaN) when either site is monomorphic among the
// haplotypes typed at both.
//
// From phased haplotypes the four counts are read directly; on bit-packed haplotypes
// (PackedHaplotypes) a pair of sites costs four popcounts per 64 haplotypes, which is what the
// decay of r² over all pairs within a distance needs.
//
// From unphased genotypes only the double heterozygote is ambiguous (AB/ab or Ab/aB). The
// haplotype frequencies are estimated by maximum likelihood with the EM algorithm (Excoffier
// and Slatkin 1995, Mol. Biol. Evol. 12, 921–927; for two sites the iteration of Hill 1974,
// Heredity 33, 229–239): the double heterozygotes are split in the proportion
// φ = x_AB x_ab / (x_AB x_ab + x_Ab x_aB) and the frequencies re-counted. The allele frequencies
// are the sample's at every step, so the likelihood is a function of x_AB alone on
// [max(0, p_A + p_B − 1), min(p_A, p_B)]; its stationary points are the roots of Hill's cubic,
// at most two of them maxima, so the EM is started at linkage equilibrium and near both ends
// of the interval, and the start of highest likelihood is kept.
//
// Complexity: O(n) for one pair (O(n/64) packed); the EM O(1) per iteration after the 3 × 3
// table is counted. The decay over m sites with w pairs each within the distance is O(m w n/64)
// phased, O(m w n) unphased, run in parallel in fixed shards, so its sums are the same on any
// number of threads.

#pragma once

#include <evolution/population/genotypes.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace evolution::population {

struct Linkage {
  double p_a = 0, p_b = 0;  // the alternate alleles' frequencies at the two sites
  double p_ab = 0;          // the frequency of the haplotype carrying both
  double d = 0, d_prime = 0, r2 = 0;
  std::size_t n = 0;        // haplotypes (phased) or individuals (unphased) typed at both
};

inline Linkage linkage_from_frequencies(double p_a, double p_b, double p_ab, std::size_t n = 0) {
  Linkage l;
  l.p_a = p_a;
  l.p_b = p_b;
  l.p_ab = p_ab;
  l.n = n;
  l.d = p_ab - p_a * p_b;
  const double var = p_a * (1 - p_a) * p_b * (1 - p_b);
  if (!(var > 0)) {
    l.d_prime = l.r2 = std::numeric_limits<double>::quiet_NaN();
    return l;
  }
  const double dmax = l.d > 0 ? std::min(p_a * (1 - p_b), (1 - p_a) * p_b) : std::min(p_a * p_b, (1 - p_a) * (1 - p_b));
  l.d_prime = dmax > 0 ? l.d / dmax : 0;
  l.r2 = l.d * l.d / var;
  return l;
}

// From the phased alleles of two sites (0, 1, −1 missing), over the haplotypes typed at both.
inline Linkage linkage_phased(std::span<const std::int8_t> a, std::span<const std::int8_t> b) {
  if (a.size() != b.size()) throw std::invalid_argument("evolution::population::linkage_phased: the same haplotypes at both sites");
  std::size_t n = 0, na = 0, nb = 0, nab = 0;
  for (std::size_t k = 0; k < a.size(); k++) {
    if (a[k] < 0 || b[k] < 0) continue;
    n++;
    na += a[k];
    nb += b[k];
    nab += a[k] & b[k];
  }
  if (n == 0) return linkage_from_frequencies(0, 0, 0, 0);
  const double x = static_cast<double>(n);
  return linkage_from_frequencies(na / x, nb / x, nab / x, n);
}

// Haplotypes packed 64 to a word, site by site: bit k of alt is set where haplotype k carries
// the alternate allele, bit k of typed where it is not missing.
struct PackedHaplotypes {
  std::size_t sequences = 0, sites = 0, words = 0;
  std::vector<std::uint64_t> alt, typed;  // sites × words

  PackedHaplotypes() = default;
  explicit PackedHaplotypes(const Haplotypes& h) : PackedHaplotypes(h.sequences, h.sites, h.alleles) {}
  PackedHaplotypes(std::size_t n, std::size_t m, std::span<const std::int8_t> by_site) : sequences(n), sites(m), words((n + 63) / 64) {
    if (by_site.size() != n * m) throw std::invalid_argument("evolution::population::PackedHaplotypes: sites × sequences alleles");
    alt.assign(m * words, 0);
    typed.assign(m * words, 0);
    for (std::size_t j = 0; j < m; j++)
      for (std::size_t k = 0; k < n; k++) {
        const std::int8_t x = by_site[j * n + k];
        if (x < 0) continue;
        typed[j * words + k / 64] |= std::uint64_t{1} << (k % 64);
        if (x > 0) alt[j * words + k / 64] |= std::uint64_t{1} << (k % 64);
      }
  }

  Linkage linkage(std::size_t i, std::size_t j) const {
    const std::uint64_t *ai = &alt[i * words], *aj = &alt[j * words], *ti = &typed[i * words], *tj = &typed[j * words];
    std::uint64_t n = 0, na = 0, nb = 0, nab = 0;
    for (std::size_t w = 0; w < words; w++) {
      n += static_cast<std::uint64_t>(std::popcount(ti[w] & tj[w]));
      na += static_cast<std::uint64_t>(std::popcount(ai[w] & tj[w]));
      nb += static_cast<std::uint64_t>(std::popcount(aj[w] & ti[w]));
      nab += static_cast<std::uint64_t>(std::popcount(ai[w] & aj[w]));
    }
    if (n == 0) return linkage_from_frequencies(0, 0, 0, 0);
    const double x = static_cast<double>(n);
    return linkage_from_frequencies(static_cast<double>(na) / x, static_cast<double>(nb) / x, static_cast<double>(nab) / x, n);
  }
};

// ---- unphased genotypes ---------------------------------------------------------------------------------

// The 3 × 3 table of two sites' dosages over the individuals typed at both: n[3 g₁ + g₂].
struct GenotypeTable {
  std::array<std::uint64_t, 9> n{};
  std::uint64_t operator()(int g1, int g2) const { return n[static_cast<std::size_t>(3 * g1 + g2)]; }
  std::uint64_t individuals() const {
    std::uint64_t s = 0;
    for (const auto x : n) s += x;
    return s;
  }
};

inline GenotypeTable genotype_table(std::span<const std::int8_t> g1, std::span<const std::int8_t> g2) {
  if (g1.size() != g2.size()) throw std::invalid_argument("evolution::population::genotype_table: the same individuals at both sites");
  GenotypeTable t;
  for (std::size_t i = 0; i < g1.size(); i++)
    if (g1[i] >= 0 && g2[i] >= 0) t.n[static_cast<std::size_t>(3 * g1[i] + g2[i])]++;
  return t;
}

// The log-likelihood of haplotype frequencies (x_AB, x_Ab, x_aB, x_ab) for a table, the
// genotypes being unions of two haplotypes drawn independently.
inline double haplotype_log_likelihood(const GenotypeTable& t, const std::array<double, 4>& x) {
  const double AB = x[0], Ab = x[1], aB = x[2], ab = x[3];
  const double P[9] = {ab * ab,     2 * aB * ab,                aB * aB,
                       2 * Ab * ab, 2 * (AB * ab + Ab * aB),    2 * AB * aB,
                       Ab * Ab,     2 * AB * Ab,                AB * AB};
  // P is indexed by 3 g₁ + g₂; g = 0 is homozygous for the reference allele (a or b).
  double l = 0;
  for (std::size_t k = 0; k < 9; k++)
    if (t.n[k]) l += static_cast<double>(t.n[k]) * (P[k] > 0 ? std::log(P[k]) : -1e300);
  return l;
}

struct HaplotypeEstimate {
  Linkage ld;
  std::array<double, 4> frequencies{};  // x_AB, x_Ab, x_aB, x_ab
  double log_likelihood = 0;
  std::size_t iterations = 0;
  bool converged = false;
};

namespace detail {

inline HaplotypeEstimate em_from(const GenotypeTable& t, double x_ab_start, double tolerance, std::size_t max_iterations) {
  const double N2 = 2.0 * static_cast<double>(t.individuals());
  // Haplotypes known from unambiguous genotypes (those of ab follow from the others).
  const double cAB = 2.0 * t(2, 2) + t(2, 1) + t(1, 2), cAb = 2.0 * t(2, 0) + t(2, 1) + t(1, 0);
  const double caB = 2.0 * t(0, 2) + t(1, 2) + t(0, 1);
  const double d = t(1, 1);
  const double pA = (cAB + cAb + d) / N2, pB = (cAB + caB + d) / N2;
  HaplotypeEstimate e;
  double AB = x_ab_start;
  std::size_t it = 0;
  bool converged = false;
  for (; it < max_iterations; it++) {
    const double Ab = pA - AB, aB = pB - AB, ab = 1 - pA - pB + AB;
    const double cis = AB * ab, trans = Ab * aB;
    const double phi = cis + trans > 0 ? cis / (cis + trans) : 0.5;
    const double next = (cAB + d * phi) / N2;
    const double change = std::fabs(next - AB);
    AB = next;
    if (change < tolerance) {
      converged = true;
      it++;
      break;
    }
  }
  e.frequencies = {AB, pA - AB, pB - AB, 1 - pA - pB + AB};
  for (double& x : e.frequencies) x = std::max(0.0, x);
  e.ld = linkage_from_frequencies(pA, pB, AB, static_cast<std::size_t>(t.individuals()));
  e.log_likelihood = haplotype_log_likelihood(t, e.frequencies);
  e.iterations = it;
  e.converged = converged;
  return e;
}

}  // namespace detail

// The maximum-likelihood haplotype frequencies of a table, and the linkage they imply.
inline HaplotypeEstimate linkage_unphased(const GenotypeTable& t, double tolerance = 1e-12, std::size_t max_iterations = 10000) {
  const std::uint64_t n = t.individuals();
  if (n == 0) return {};
  const double N2 = 2.0 * static_cast<double>(n);
  const double pA = (2.0 * (t(2, 0) + t(2, 1) + t(2, 2)) + t(1, 0) + t(1, 1) + t(1, 2)) / N2;
  const double pB = (2.0 * (t(0, 2) + t(1, 2) + t(2, 2)) + t(0, 1) + t(1, 1) + t(2, 1)) / N2;
  const double lo = std::max(0.0, pA + pB - 1), hi = std::min(pA, pB);
  HaplotypeEstimate best = detail::em_from(t, pA * pB, tolerance, max_iterations);
  if (t(1, 1) > 0 && hi > lo) {
    for (const double start : {lo + 1e-3 * (hi - lo), hi - 1e-3 * (hi - lo)}) {
      HaplotypeEstimate e = detail::em_from(t, start, tolerance, max_iterations);
      if (e.log_likelihood > best.log_likelihood + 1e-12) best = e;
    }
  }
  return best;
}

inline HaplotypeEstimate linkage_unphased(std::span<const std::int8_t> g1, std::span<const std::int8_t> g2) {
  return linkage_unphased(genotype_table(g1, g2));
}

// ---- decay with distance --------------------------------------------------------------------------------

// Mean r² of the pairs whose distance lies in [begin, end).
struct DecayBin {
  double begin = 0, end = 0;
  double mean_r2 = 0;
  std::size_t pairs = 0;
};

namespace detail {

// r²(i, j) for sites i < j whose positions differ by less than max_distance and whose groups
// (chromosomes) agree, binned by distance. Shards of 64 first sites each sum their own bins,
// merged in shard order.
template <class R2>
std::vector<DecayBin> decay(std::span<const double> positions, std::span<const std::uint32_t> group, double max_distance, std::size_t bins, R2&& r2,
                            calx::scheduler::ThreadPool* pool) {
  if (bins == 0 || !(max_distance > 0)) throw std::invalid_argument("evolution::population::ld_decay: bins and a positive distance");
  const std::size_t m = positions.size(), shard = 64, shards = (m + shard - 1) / shard;
  const double width = max_distance / static_cast<double>(bins);
  std::vector<double> sum(shards * bins, 0.0);
  std::vector<std::size_t> count(shards * bins, 0);
  const auto work = [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t s = begin; s < end; s++)
      for (std::size_t i = s * shard; i < std::min(m, (s + 1) * shard); i++)
        for (std::size_t j = i + 1; j < m; j++) {
          if (!group.empty() && group[j] != group[i]) break;
          const double dist = positions[j] - positions[i];
          if (dist >= max_distance) break;
          const double x = r2(i, j);
          if (!std::isfinite(x)) continue;
          const std::size_t b = std::min(bins - 1, static_cast<std::size_t>(dist / width));
          sum[s * bins + b] += x;
          count[s * bins + b]++;
        }
  };
  if (pool) calx::scheduler::parallel_for(*pool, shards, 1, work);
  else work(0, shards, 0);
  std::vector<DecayBin> out(bins);
  for (std::size_t b = 0; b < bins; b++) {
    out[b].begin = static_cast<double>(b) * width;
    out[b].end = static_cast<double>(b + 1) * width;
    double total = 0;
    std::size_t pairs = 0;
    for (std::size_t s = 0; s < shards; s++) {
      total += sum[s * bins + b];
      pairs += count[s * bins + b];
    }
    out[b].pairs = pairs;
    out[b].mean_r2 = pairs ? total / static_cast<double>(pairs) : std::numeric_limits<double>::quiet_NaN();
  }
  return out;
}

}  // namespace detail

// The decay of r² with distance among phased haplotypes on one chromosome.
inline std::vector<DecayBin> ld_decay(const Haplotypes& h, double max_distance, std::size_t bins, calx::scheduler::ThreadPool* pool = nullptr) {
  const PackedHaplotypes packed(h);
  return detail::decay(h.positions, {}, max_distance, bins, [&](std::size_t i, std::size_t j) { return packed.linkage(i, j).r2; }, pool);
}

// The decay of r² with distance (in base pairs) within each chromosome of a genotype matrix:
// from the haplotypes when phased, by the EM of haplotype frequencies otherwise.
inline std::vector<DecayBin> ld_decay(const Genotypes& g, double max_distance, std::size_t bins, calx::scheduler::ThreadPool* pool = nullptr) {
  const std::vector<double> pos(g.positions.begin(), g.positions.end());
  if (g.phased()) {
    const PackedHaplotypes packed(2 * g.samples, g.sites, g.haplotypes);
    return detail::decay(pos, g.chromosome, max_distance, bins, [&](std::size_t i, std::size_t j) { return packed.linkage(i, j).r2; }, pool);
  }
  return detail::decay(pos, g.chromosome, max_distance, bins, [&](std::size_t i, std::size_t j) { return linkage_unphased(g.site(i), g.site(j)).ld.r2; }, pool);
}

// The r² of every pair of the sites [first, last) of phased haplotypes, row-major, for drawing
// (NaN on the diagonal of a monomorphic site).
inline std::vector<double> ld_matrix(const PackedHaplotypes& packed, std::size_t first, std::size_t last) {
  if (first > last || last > packed.sites) throw std::invalid_argument("evolution::population::ld_matrix: a range of sites");
  const std::size_t k = last - first;
  std::vector<double> r(k * k);
  for (std::size_t a = 0; a < k; a++)
    for (std::size_t b = a; b < k; b++) r[a * k + b] = r[b * k + a] = packed.linkage(first + a, first + b).r2;
  return r;
}

}  // namespace evolution::population
