// The k-mer spectrum of a sample of reads (the histogram of k-mer counts) and what it says
// about the genome sequenced: the threshold below which k-mers are errors, the coverage of a
// k-mer, the size of the genome and its heterozygosity.
//
// - The spectrum h(c) is the number of distinct canonical k-mers counted c times. Reads drawn
//   uniformly from a genome cover a k-mer of it by a number of reads that is Poisson
//   distributed (Lander and Waterman 1988, Genomics 2, 231–239): a single-copy k-mer of a
//   haploid genome by Poisson(λ), with the k-mer coverage λ = c (L − k + 1)/L of reads of
//   length L at base coverage c, thinned by the chance that the read has no error inside the
//   k-mer. Sequencing errors make k-mers that occur once or a few times. The spectrum is
//   therefore a steep descent of errors from c = 1, a first minimum, and a peak at λ, with
//   smaller peaks at 2λ, 3λ, … from repeats and, in a diploid genome, a peak at λ/2 from the
//   k-mers that cover a heterozygous site, each seen on one haplotype only.
// - The error threshold is the first minimum: the least c ≥ 1 with h(c + 1) ≥ h(c).
// - The peak is the mode of h at or above the threshold.
// - The model is a simplification of GenomeScope's (Vurture et al. 2017, "GenomeScope:
//   fast reference-free genome profiling from short reads", Bioinformatics 33, 2202–2204),
//   stated as such: above the threshold, h is fitted as a mixture of Poisson laws with means
//   λ/2, λ, 2λ, 3λ and 4λ (GenomeScope fits negative binomials whose weights are functions of
//   the heterozygosity and the duplication levels; here the five weights are free, the means
//   share one λ, and the law is Poisson), by EM (Dempster, Laird and Rubin 1977, J. R. Stat.
//   Soc. B 39, 1–38): the E-step gives each count its responsibilities r_j(c), the M-step
//   w_j = Σ h r_j / Σ h and λ = Σ c h / Σ_c h Σ_j m_j r_j, which is the maximiser of the
//   expected log-likelihood in λ for means m_j λ. The fit is over [threshold, 5 · peak]; the
//   truncation of the laws at the threshold is ignored, which is exact enough when the
//   threshold lies well below λ/2, as it does at the coverages assembly needs.
//   A spectrum with peaks at x and 2x is a diploid genome with heterozygous k-mers at x or a
//   haploid one with duplications at 2x: the data cannot tell them apart. Two fits are made
//   over the same range, from λ₀ = the peak and from λ₀ = twice the peak, and the more likely
//   is kept, the first on a tie; a caller who knows the coverage gives it
//   (FitOptions::coverage), and the range is then [threshold, 5λ₀].
// - Genome size: the k-mers at and above the threshold, counted with multiplicity, divided
//   by λ (Li and Waterman 2003, "Estimating the repeat structure and length of DNA sequences
//   using ℓ-tuples", Genome Research 13, 1916–1922; Liu et al. 2013, GCE, arXiv:1308.2012):
//   G ≈ Σ_{c ≥ t} c h(c) / λ, the haploid length; a heterozygous window contributes two
//   k-mers at λ/2, which is λ in all.
// - Heterozygosity: a window of k bases is heterozygous when one of its bases is, which for
//   independent sites of rate r has probability 1 − (1 − r)^k; it gives two distinct k-mers
//   at λ/2, a homozygous window one at λ (or one at jλ shared by j copies). With W_j the
//   distinct k-mers of the component at m_j λ, the fraction of heterozygous windows is
//   f = (W_{1/2}/2) / (W_{1/2}/2 + Σ_{j ≥ 1} j W_j), and r = 1 − (1 − f)^{1/k}. Indels count
//   as heterozygous sites, as they do in GenomeScope.
//
// Limits: Poisson, not overdispersed, coverage, so on real reads λ and the weights are those
// of the best Poisson fit, not of GenomeScope's; no ploidy above two; repeats of five copies
// or more are counted in the genome size but not fitted. The EM costs O(I · range) for I
// iterations (at most 2000; convergence in λ to 10⁻¹⁰).

#pragma once

#include <genomics/kmers/counter.hpp>

#include <calx/scheduler/parallel_for.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace genomics::kmers {

struct Spectrum {
  int k = 0;
  // histogram[c]: the distinct k-mers counted c times, for 1 ≤ c < size − 1; the last bin
  // holds every count at or above its index; histogram[0] = 0.
  std::vector<std::uint64_t> histogram;
  std::uint64_t distinct = 0, total = 0;  // distinct k-mers, and their occurrences

  std::uint64_t at(std::size_t c) const { return c < histogram.size() ? histogram[c] : 0; }
  // Σ_{c < t} c h(c): occurrences of the k-mers counted fewer than t times (t ≤ the last bin).
  std::uint64_t occurrences_below(std::size_t t) const {
    std::uint64_t s = 0;
    for (std::size_t c = 1; c < t && c + 1 < histogram.size(); c++) s += c * histogram[c];
    return s;
  }
  std::uint64_t distinct_below(std::size_t t) const {
    std::uint64_t s = 0;
    for (std::size_t c = 1; c < t && c < histogram.size(); c++) s += histogram[c];
    return s;
  }
};

// The spectrum of the counts, with exact bins up to max_count − 1.
inline Spectrum spectrum(const KmerCounter& counts, calx::scheduler::ThreadPool& pool, std::uint32_t max_count = 10000) {
  if (max_count < 2) throw std::invalid_argument("genomics::kmers::spectrum: max_count >= 2");
  const std::size_t P = counts.partitions(), H = max_count + 1;
  std::vector<std::vector<std::uint64_t>> part(P);
  calx::scheduler::parallel_for(pool, P, 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t p = b; p < e; p++) {
      auto& h = part[p];
      h.assign(H, 0);
      counts.table(p).for_each([&](Kmer, std::uint32_t c) { h[std::min<std::size_t>(c, max_count)]++; });
    }
  });
  Spectrum s;
  s.k = counts.k();
  s.histogram.assign(H, 0);
  for (const auto& h : part)
    for (std::size_t c = 0; c < H; c++) s.histogram[c] += h[c];
  s.distinct = counts.distinct();
  s.total = counts.total();
  return s;
}

// The first minimum of the spectrum: the least c ≥ 1 with h(c + 1) ≥ h(c); 1 if there is none.
inline std::uint32_t error_threshold(const Spectrum& s) {
  for (std::size_t c = 1; c + 2 < s.histogram.size(); c++)
    if (s.histogram[c + 1] >= s.histogram[c]) return static_cast<std::uint32_t>(c);
  return 1;
}

// The mode of the spectrum at or above `from` (below the last bin).
inline std::uint32_t coverage_peak(const Spectrum& s, std::uint32_t from) {
  std::size_t best = from;
  for (std::size_t c = from; c + 1 < s.histogram.size(); c++)
    if (s.histogram[c] > s.histogram[best]) best = c;
  return static_cast<std::uint32_t>(best);
}

// The k-mer coverage of reads of length L at base coverage c, and back.
constexpr double kmer_coverage(double base_coverage, double read_length, int k) { return base_coverage * (read_length - k + 1) / read_length; }
constexpr double base_coverage(double kmer_cov, double read_length, int k) { return kmer_cov * read_length / (read_length - k + 1); }

struct FitOptions {
  std::uint32_t threshold = 0;   // 0: the first minimum
  double coverage = 0;           // the k-mer coverage if known (λ₀); 0: from the peak
  int max_iterations = 2000;
};

struct SpectrumFit {
  static constexpr std::array<double, 5> MULTIPLES = {0.5, 1, 2, 3, 4};

  std::uint32_t threshold = 0;          // counts below are errors
  std::uint32_t peak = 0;               // the mode at or above the threshold
  double coverage = 0;                  // λ, the k-mer coverage of a single-copy homozygous k-mer
  std::array<double, 5> kmers{};        // distinct k-mers of the components at λ/2, λ, 2λ, 3λ, 4λ
  double log_likelihood = 0;            // of the fitted counts
  double genome_size = 0;               // Σ_{c ≥ threshold} c h(c) / λ
  double genome_size_at_peak = 0;       // the same over the peak
  double heterozygous_windows = 0;      // f
  double heterozygosity = 0;            // r = 1 − (1 − f)^{1/k}, per base
  std::uint64_t error_kmers = 0;        // distinct k-mers below the threshold
  std::uint64_t solid_kmers = 0;        // distinct k-mers at or above it
  int iterations = 0;
  bool converged = false;
};

namespace detail {

struct Em {
  double lambda = 0, log_likelihood = -INFINITY;
  std::array<double, 5> w{};
  int iterations = 0;
  bool converged = false;
};

inline Em spectrum_em(const Spectrum& s, std::size_t lo, std::size_t hi, double lambda0, int max_iterations) {
  constexpr auto& m = SpectrumFit::MULTIPLES;
  Em em;
  em.lambda = lambda0;
  em.w = {0.1, 0.6, 0.1, 0.1, 0.1};
  // log c! for the range.
  std::vector<double> log_factorial(hi + 1, 0);
  for (std::size_t c = 2; c <= hi; c++) log_factorial[c] = log_factorial[c - 1] + std::log(static_cast<double>(c));
  double n = 0;
  for (std::size_t c = lo; c <= hi; c++) n += static_cast<double>(s.at(c));
  if (n == 0 || !(lambda0 > 0)) return em;
  std::array<double, 5> r{};
  for (int it = 0; it < max_iterations; it++) {
    std::array<double, 5> sw{};
    double sc = 0, sm = 0, ll = 0;
    for (std::size_t c = lo; c <= hi; c++) {
      const double h = static_cast<double>(s.at(c));
      if (h == 0) continue;
      const double x = static_cast<double>(c);
      double top = -INFINITY;
      for (std::size_t j = 0; j < 5; j++) {
        const double mu = m[j] * em.lambda;
        r[j] = em.w[j] > 0 ? std::log(em.w[j]) + x * std::log(mu) - mu - log_factorial[c] : -INFINITY;
        top = std::max(top, r[j]);
      }
      double z = 0;
      for (auto& v : r) z += (v = std::exp(v - top));
      ll += h * (top + std::log(z));
      for (std::size_t j = 0; j < 5; j++) {
        const double q = h * r[j] / z;
        sw[j] += q;
        sm += q * m[j];
      }
      sc += h * x;
    }
    const double lambda = sc / sm;
    double change = std::fabs(lambda - em.lambda) / em.lambda;
    for (std::size_t j = 0; j < 5; j++) {
      const double w = sw[j] / n;
      change = std::max(change, std::fabs(w - em.w[j]));
      em.w[j] = w;
    }
    em.lambda = lambda;
    em.log_likelihood = ll;
    em.iterations = it + 1;
    if (change < 1e-10) { em.converged = true; break; }
  }
  return em;
}

}  // namespace detail

inline SpectrumFit fit_spectrum(const Spectrum& s, const FitOptions& options = {}) {
  SpectrumFit f;
  f.threshold = options.threshold ? options.threshold : error_threshold(s);
  f.peak = coverage_peak(s, f.threshold);
  f.error_kmers = s.distinct_below(f.threshold);
  f.solid_kmers = s.distinct - f.error_kmers;
  const std::size_t last = s.histogram.size() - 2;  // the last exact bin
  // One range for every fit, so that their likelihoods compare.
  const double base = options.coverage > 0 ? options.coverage : static_cast<double>(f.peak);
  const std::size_t hi = std::max<std::size_t>(f.threshold, std::min<std::size_t>(last, static_cast<std::size_t>(std::ceil(5 * base))));
  const auto fit = [&](double lambda0) { return detail::spectrum_em(s, f.threshold, hi, lambda0, options.max_iterations); };
  detail::Em em;
  if (options.coverage > 0) {
    em = fit(options.coverage);
  } else {
    em = fit(f.peak);
    const detail::Em twice = fit(2.0 * f.peak);
    if (twice.log_likelihood > em.log_likelihood + 1e-9 * std::fabs(em.log_likelihood)) em = twice;
  }
  f.coverage = em.lambda;
  f.log_likelihood = em.log_likelihood;
  f.iterations = em.iterations;
  f.converged = em.converged;
  double fitted = 0;
  for (std::size_t c = f.threshold; c <= hi; c++) fitted += static_cast<double>(s.at(c));
  for (std::size_t j = 0; j < 5; j++) f.kmers[j] = em.w[j] * fitted;
  const double occurrences = static_cast<double>(s.total - s.occurrences_below(f.threshold));
  if (f.coverage > 0) f.genome_size = occurrences / f.coverage;
  if (f.peak > 0) f.genome_size_at_peak = occurrences / f.peak;
  const double het = f.kmers[0] / 2, hom = f.kmers[1] + 2 * f.kmers[2] + 3 * f.kmers[3] + 4 * f.kmers[4];
  if (het + hom > 0) {
    f.heterozygous_windows = het / (het + hom);
    f.heterozygosity = 1 - std::pow(1 - f.heterozygous_windows, 1.0 / s.k);
  }
  return f;
}

}  // namespace genomics::kmers
