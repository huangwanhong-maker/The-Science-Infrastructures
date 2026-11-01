// Read depth in bins, with mappability and the GC content of fragments: the observations of copy
// number (CNVnator; Abyzov et al. 2011, Genome Research 21, 974–984: depth in bins, corrected for
// GC, segmented).
//
// - Counts. A fragment is counted once, in the bin of its start (the unclipped start of the
//   forward mate), when its pair is proper and both mates have mapping quality ≥ `min_mapq`: the
//   observation of Benjamini and Speed's model below. Mates of lower mapping quality are counted
//   apart: the fraction of a bin's mates that are ambiguous is its mappability, read from the data
//   as CNVnator reads q0 (the fraction of reads of mapping quality 0); a bin with more than
//   `max_ambiguous` (0.5) of them, or with more than `max_n` (0.1) of its bases N, is not
//   mappable and is left out of segmentation; a mappable bin's exposure is scaled by its share of
//   confidently placed mates (1 − the ambiguous fraction), since only those are counted.
// - GC (Benjamini and Speed 2012, "Summarizing and correcting the GC content bias in
//   high-throughput sequencing", Nucleic Acids Research 40, e72, §2: the rate of fragments at a
//   position is a function of the GC content of the fragment that would start there). With F the
//   library's median template and g(p) the GC stratum (a hundredth) of reference bases
//   [p, p + F), the rate of stratum g is λ_g = N_g / P_g: the fragments counted whose start has
//   stratum g over the mappable positions of stratum g (strata with fewer than `min_positions`
//   positions pooled with their neighbours until they have them). The expected count of a bin is
//   its exposure e = Σ_{p ∈ bin} λ_{g(p)} over its positions whose fragment has no N; without the
//   GC model, λ̄ = ΣN/ΣP for every position (the correction of Yoon et al. 2009, Genome Research 19,
//   1586–1592, by the median of bins of one GC, is the same idea at the scale of bins).
// - Scale. The exposures are scaled so that the median ratio count/exposure over mappable bins is
//   1: the typical copy number of the genome is taken as its baseline (2 in a diploid), which
//   assumes that most of the genome is at it (*Caveat*: a genome mostly gained or lost, a triploid
//   one, needs a ploidy of its own; the somatic module estimates it).
// - Under the model, a bin's count is Poisson with mean (c/2) e for copy number c, which
//   copy_number.hpp segments and the checks hold to: the dispersion Σ (n − e)²/e over mappable
//   bins at copy number 2 is near its Poisson value, the number of bins, once GC is corrected.
//
// Complexity: O(n) for a reference of n bases (one sliding window of GC), in parallel over chunks
// of bins whose histograms are added in order, so the result does not depend on the threads.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/structural/pairs.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace genomics::structural {

struct DepthOptions {
  std::size_t bin = 1000;
  int min_mapq = 20;
  double max_ambiguous = 0.5;
  double max_n = 0.1;
  std::size_t fragment = 0;         // the fragment of the GC model; 0: the library's median
  bool gc = true;
  double min_positions = 2000;      // a stratum's positions before it stands alone
};

struct DepthBin {
  std::uint32_t contig = 0;
  std::size_t begin = 0, end = 0;
  double count = 0;      // fragments starting in it
  double exposure = 0;   // expected fragments at the baseline copy number
  double gc = 0;         // of its bases
  double ambiguous = 0;  // the fraction of its mates of low mapping quality
  double n = 0;          // the fraction of its bases that are N
  bool mappable = true;
  double ratio() const { return exposure > 0 ? count / exposure : 0; }
};

struct Depth {
  std::vector<DepthBin> bins;          // by contig and position
  std::vector<std::size_t> first;      // the first bin of each contig, and the total at the end
  std::vector<double> rate;            // λ_g, g = 0 .. 100 (after pooling), scaled
  std::vector<double> positions, fragments;  // P_g and N_g
  std::size_t fragment = 0;
  double scale = 1;                    // the factor that set the median ratio to 1
  double dispersion(std::size_t from = 0, std::size_t to = SIZE_MAX) const {
    double d = 0;
    for (std::size_t i = from; i < std::min(to, bins.size()); i++)
      if (bins[i].mappable && bins[i].exposure > 0) d += (bins[i].count - bins[i].exposure) * (bins[i].count - bins[i].exposure) / bins[i].exposure;
    return d;
  }
};

namespace detail {

inline std::size_t stratum(std::size_t gc, std::size_t F) { return (gc * 100 + F / 2) / F; }

}  // namespace detail

inline Depth read_depth(calx::scheduler::ThreadPool& pool, const std::vector<Codes>& contigs, const MateMappings& mates, const std::vector<PairRecord>& pairs, const Library& lib,
                        const DepthOptions& o = {}) {
  Depth d;
  const std::size_t W = std::max<std::size_t>(1, o.bin);
  d.fragment = o.fragment ? o.fragment : static_cast<std::size_t>(std::max(1.0, lib.median));
  const std::size_t F = d.fragment;
  for (std::uint32_t c = 0; c < contigs.size(); c++) {
    d.first.push_back(d.bins.size());
    for (std::size_t b = 0; b < contigs[c].size(); b += W) d.bins.push_back({c, b, std::min(contigs[c].size(), b + W)});
  }
  d.first.push_back(d.bins.size());
  const auto bin_of = [&](std::uint32_t c, std::size_t p) { return d.first[c] + std::min(p / W, d.first[c + 1] - d.first[c] - 1); };
  // Fragments and mates per bin, and the strata of the fragments' starts.
  std::vector<double> total(d.bins.size(), 0), low(d.bins.size(), 0);
  for (const auto& ms : mates) {
    const Mapping* m = mapping::primary_of(ms);
    if (!m) continue;
    const std::size_t b = bin_of(m->contig, m->target_begin);
    total[b]++;
    low[b] += m->mapq < o.min_mapq;
  }
  std::vector<std::pair<std::uint32_t, std::size_t>> starts;
  for (const auto& p : pairs) {
    if (p.type != PairClass::proper || p.mate[0]->mapq < o.min_mapq || p.mate[1]->mapq < o.min_mapq) continue;
    const Mapping& f = p.mate[0]->reverse ? *p.mate[1] : *p.mate[0];
    const std::int64_t s = std::max<std::int64_t>(0, unclipped_begin(f));
    starts.push_back({f.contig, static_cast<std::size_t>(s)});
    d.bins[bin_of(f.contig, static_cast<std::size_t>(s))].count++;
  }
  // Per bin: GC, N, and the histogram of the strata of its clean positions (computed in chunks).
  const std::size_t S = 101;
  const std::size_t B = d.bins.size(), chunk = 64, chunks = (B + chunk - 1) / chunk;
  std::vector<double> hist(chunks * S, 0);
  calx::scheduler::parallel_for(pool, chunks, 1, [&](std::size_t lo, std::size_t hi, std::size_t) {
    for (std::size_t k = lo; k < hi; k++)
      for (std::size_t i = k * chunk; i < std::min(B, (k + 1) * chunk); i++) {
        DepthBin& bin = d.bins[i];
        const Codes& x = contigs[bin.contig];
        std::size_t gc = 0, n = 0;
        for (std::size_t p = bin.begin; p < bin.end; p++) {
          gc += x[p] == sequence::C || x[p] == sequence::G;
          n += x[p] > 3;
        }
        const double len = static_cast<double>(bin.end - bin.begin);
        bin.gc = n < bin.end - bin.begin ? static_cast<double>(gc) / (len - static_cast<double>(n)) : 0;
        bin.n = static_cast<double>(n) / len;
        bin.ambiguous = total[i] > 0 ? low[i] / total[i] : 0;
        bin.mappable = bin.n <= o.max_n && bin.ambiguous <= o.max_ambiguous;
        if (!bin.mappable || bin.begin + F > x.size()) continue;
        // The sliding window [p, p + F).
        std::size_t wgc = 0, wn = 0;
        for (std::size_t p = bin.begin; p < bin.begin + F; p++) { wgc += x[p] == sequence::C || x[p] == sequence::G; wn += x[p] > 3; }
        for (std::size_t p = bin.begin; p < bin.end && p + F <= x.size(); p++) {
          if (p > bin.begin) {
            const Code out = x[p - 1], in = x[p + F - 1];
            wgc += (in == sequence::C || in == sequence::G);
            wgc -= (out == sequence::C || out == sequence::G);
            wn += in > 3;
            wn -= out > 3;
          }
          if (wn == 0) hist[k * S + detail::stratum(wgc, F)]++;
        }
      }
  });
  d.positions.assign(S, 0);
  for (std::size_t k = 0; k < chunks; k++)
    for (std::size_t g = 0; g < S; g++) d.positions[g] += hist[k * S + g];
  d.fragments.assign(S, 0);
  for (const auto& [c, s] : starts) {
    const Codes& x = contigs[c];
    if (s + F > x.size() || !d.bins[bin_of(c, s)].mappable) continue;
    std::size_t gc = 0, n = 0;
    for (std::size_t p = s; p < s + F; p++) { gc += x[p] == sequence::C || x[p] == sequence::G; n += x[p] > 3; }
    if (n == 0) d.fragments[detail::stratum(gc, F)]++;
  }
  // λ_g, pooling strata with their neighbours until they hold `min_positions` positions.
  double all_n = 0, all_p = 0;
  for (std::size_t g = 0; g < S; g++) { all_n += d.fragments[g]; all_p += d.positions[g]; }
  const double mean_rate = all_p > 0 ? all_n / all_p : 0;
  d.rate.assign(S, mean_rate);
  if (o.gc)
    for (std::size_t g = 0; g < S; g++) {
      double n = d.fragments[g], p = d.positions[g];
      std::size_t r = 0;
      while (p < o.min_positions && r < S) {
        r++;
        if (g >= r) { n += d.fragments[g - r]; p += d.positions[g - r]; }
        if (g + r < S) { n += d.fragments[g + r]; p += d.positions[g + r]; }
      }
      d.rate[g] = p > 0 ? n / p : mean_rate;
    }
  // Exposures.
  calx::scheduler::parallel_for(pool, chunks, 1, [&](std::size_t lo, std::size_t hi, std::size_t) {
    for (std::size_t k = lo; k < hi; k++)
      for (std::size_t i = k * chunk; i < std::min(B, (k + 1) * chunk); i++) {
        DepthBin& bin = d.bins[i];
        const Codes& x = contigs[bin.contig];
        bin.exposure = 0;
        if (bin.begin + F > x.size()) continue;
        std::size_t wgc = 0, wn = 0;
        for (std::size_t p = bin.begin; p < bin.begin + F; p++) { wgc += x[p] == sequence::C || x[p] == sequence::G; wn += x[p] > 3; }
        for (std::size_t p = bin.begin; p < bin.end && p + F <= x.size(); p++) {
          if (p > bin.begin) {
            const Code out = x[p - 1], in = x[p + F - 1];
            wgc += (in == sequence::C || in == sequence::G);
            wgc -= (out == sequence::C || out == sequence::G);
            wn += in > 3;
            wn -= out > 3;
          }
          if (wn == 0) bin.exposure += d.rate[detail::stratum(wgc, F)];
        }
        // Mappability from the data: only the confidently placed share of a bin's fragments is counted.
        bin.exposure *= 1 - bin.ambiguous;
        if (bin.exposure <= 0) bin.mappable = false;
      }
  });
  // The scale: the median ratio of mappable bins set to 1.
  std::vector<double> ratios;
  for (const auto& b : d.bins)
    if (b.mappable && b.exposure > 0) ratios.push_back(b.count / b.exposure);
  if (!ratios.empty()) {
    std::nth_element(ratios.begin(), ratios.begin() + static_cast<std::ptrdiff_t>(ratios.size() / 2), ratios.end());
    const double m = ratios[ratios.size() / 2];
    if (m > 0) {
      d.scale = m;
      for (auto& b : d.bins) b.exposure *= m;
      for (auto& r : d.rate) r *= m;
    }
  }
  return d;
}

}  // namespace genomics::structural
