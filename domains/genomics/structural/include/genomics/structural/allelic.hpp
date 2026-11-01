// The alleles beside the depth: B-allele fractions at the sample's variant sites, runs of
// homozygosity by bcftools roh's hidden Markov model, and the reading of both with the copy number,
// which tells a copy-neutral loss of heterozygosity from a deletion (theory/structural-variants.md, §5).
//
// - Sites. The passing biallelic SNVs of the domain's caller (variant/caller.hpp, Li 2011's
//   genotype likelihoods): at each, the reads of the reference and of the alternate allele, the
//   B-allele fraction BAF = alt/(ref + alt), the genotype, and the likelihoods L(RR), L(RA), L(AA)
//   from the PL field (L = 10^(−PL/10), relative).
// - Runs of homozygosity (Narasimhan et al. 2016, "BCFtools/RoH: a hidden Markov model approach for
//   detecting autozygosity from next-generation sequencing data", Bioinformatics 32, 1749–1751;
//   the HMM of bcftools' vcfroh.c at bc63cf2 and its manual). Two states over the sites in order,
//   AZ (autozygous: one haplotype's alleles twice) and HW (Hardy–Weinberg), with f the frequency of
//   the alternate allele:
//     P(D | AZ) = (1 − f) L(RR) + f L(AA),
//     P(D | HW) = (1 − f)² L(RR) + 2 f (1 − f) L(RA) + f² L(AA),
//   and transitions P(AZ | HW) = 6.7 · 10⁻⁸ and P(HW | AZ) = 5 · 10⁻⁹ per site (bcftools' -a and -H
//   defaults), multiplied by the distance to the next site times a recombination rate when one is
//   given (bcftools' -M; capped at 1). Without a population, f is one value for every site:
//   bcftools' --AF-dflt, 0.4 in its manual's example. The states are decoded by Viterbi
//   (hmm/chain.hpp), and each run's quality is the mean posterior of its state, Phred-scaled.
// - Reading with depth. A run of homozygosity where the copy number is 2 is a copy-neutral loss of
//   heterozygosity (or autozygosity), where it is 1 a hemizygous deletion, where it is ≥ 3 a gain of
//   one haplotype; outside runs, at copy number c ≥ 2, the minor allele's copies m are estimated from
//   the heterozygous sites by maximum likelihood, each site's alternate reads a ~ ½ Bin(n, m/c) +
//   ½ Bin(n, 1 − m/c) (which haplotype carries the alternate is unknown), m = 1 … ⌊c/2⌋: an allelic
//   imbalance (2 + 1 in a trisomy) is seen as BAF at 1/3 and 2/3.
//
// Complexity O(s) for s sites. Limits: the caller emits variant sites only, so a run is read from
// the homozygous sites it holds and its edges are known to the nearest site; the density of sites
// is not used (a deletion also thins them); with a population's allele frequencies (a VCF of
// them) the model is bcftools' own.

#pragma once

#include <genomics/hmm/chain.hpp>
#include <genomics/structural/copy_number.hpp>
#include <genomics/variant/caller.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace genomics::structural {

struct AlleleSite {
  std::uint32_t contig = 0;
  std::size_t position = 0;
  std::size_t reference = 0, alternate = 0;  // reads
  int genotype = 0;                           // alternate copies called
  std::array<double, 3> log10_likelihood{};   // RR, RA, AA (from PL)
  double baf() const { return reference + alternate ? static_cast<double>(alternate) / static_cast<double>(reference + alternate) : 0; }
};

// The sites of the caller's passing biallelic SNVs, in order.
inline std::vector<AlleleSite> allele_sites(const std::vector<variant::Call>& calls) {
  std::vector<AlleleSite> out;
  for (const auto& c : calls) {
    if (!c.passed() || c.indel || c.alts.size() != 1 || c.pl.size() != 3 || c.allele_depths.size() != 2) continue;
    AlleleSite s;
    s.contig = c.contig;
    s.position = c.position;
    s.reference = c.allele_depths[0];
    s.alternate = c.allele_depths[1];
    s.genotype = c.genotype[0] + c.genotype[1];
    for (int g = 0; g < 3; g++) s.log10_likelihood[static_cast<std::size_t>(g)] = -c.pl[static_cast<std::size_t>(g)] / 10.0;
    out.push_back(s);
  }
  return out;
}

struct RohOptions {
  double allele_frequency = 0.4;  // f, bcftools --AF-dflt as in its manual's example
  double hw_to_az = 6.7e-8;       // bcftools -a
  double az_to_hw = 5e-9;         // bcftools -H
  double recombination = 0;       // per base (bcftools -M); 0: the transitions per site
  std::size_t min_sites = 10;     // a run reported must hold this many sites
};

struct Run {
  std::uint32_t contig = 0;
  std::size_t begin = 0, end = 0;     // from its first site to one past its last
  std::size_t first = 0, last = 0;    // its sites [first, last] in the sites' order
  std::size_t sites = 0, heterozygous = 0;
  double quality = 0;                  // the mean posterior of AZ, Phred-scaled
};

// bcftools roh over the sites (all contigs, each its own chain): the runs in the AZ state.
inline std::vector<Run> runs_of_homozygosity(const std::vector<AlleleSite>& sites, const RohOptions& o = {}) {
  std::vector<Run> out;
  const double f = o.allele_frequency, ln10 = std::log(10.0);
  for (std::size_t lo = 0; lo < sites.size();) {
    std::size_t hi = lo;
    while (hi < sites.size() && sites[hi].contig == sites[lo].contig) hi++;
    const std::size_t T = hi - lo;
    hmm::Chain c(2, T, o.recombination <= 0);  // state 0 HW, 1 AZ
    c.log_initial = {std::log(0.5), std::log(0.5)};
    for (std::size_t t = 0; t < T; t++) {
      const auto& L = sites[lo + t].log10_likelihood;
      const double top = std::max({L[0], L[1], L[2]});
      const double rr = std::pow(10.0, L[0] - top), ra = std::pow(10.0, L[1] - top), aa = std::pow(10.0, L[2] - top);
      c.e(t, 0) = std::log((1 - f) * (1 - f) * rr + 2 * f * (1 - f) * ra + f * f * aa) + top * ln10;
      c.e(t, 1) = std::log((1 - f) * rr + f * aa) + top * ln10;
    }
    const auto set = [&](std::size_t t, double scale) {
      const double to_az = std::min(1.0, o.hw_to_az * scale), to_hw = std::min(1.0, o.az_to_hw * scale);
      c.a(t, 0, 0) = std::log1p(-to_az);
      c.a(t, 0, 1) = std::log(to_az);
      c.a(t, 1, 0) = std::log(to_hw);
      c.a(t, 1, 1) = std::log1p(-to_hw);
    };
    if (o.recombination <= 0) set(0, 1);
    else
      for (std::size_t t = 0; t + 1 < T; t++) set(t, static_cast<double>(sites[lo + t + 1].position - sites[lo + t].position) * o.recombination);
    const hmm::Path path = hmm::viterbi(c);
    const hmm::Posterior post = hmm::posterior(c);
    for (std::size_t t = 0; t < T && !path.states.empty();) {
      if (path.states[t] != 1) { t++; continue; }
      std::size_t u = t;
      Run r;
      r.contig = sites[lo].contig;
      double q = 0;
      while (u < T && path.states[u] == 1) {
        r.heterozygous += sites[lo + u].genotype == 1;
        q += post.at(u, 1);
        u++;
      }
      r.first = lo + t;
      r.last = lo + u - 1;
      r.sites = u - t;
      r.begin = sites[r.first].position;
      r.end = sites[r.last].position + 1;
      const double mean = q / static_cast<double>(r.sites);
      r.quality = mean >= 1 ? 99 : std::min(99.0, -10 * std::log10(1 - mean));
      if (r.sites >= o.min_sites) out.push_back(r);
      t = u;
    }
    lo = hi;
  }
  return out;
}

enum class AllelicState : std::uint8_t { heterozygous, copy_neutral_loh, hemizygous_deletion, homozygous_deletion, loh_gain, imbalance };

inline const char* allelic_state_name(AllelicState s) {
  switch (s) {
    case AllelicState::heterozygous: return "heterozygous";
    case AllelicState::copy_neutral_loh: return "copy-neutral LOH";
    case AllelicState::hemizygous_deletion: return "hemizygous deletion";
    case AllelicState::homozygous_deletion: return "homozygous deletion";
    case AllelicState::loh_gain: return "LOH with gain";
    case AllelicState::imbalance: return "allelic imbalance";
  }
  return "?";
}

// A stretch of a contig with its copy number and its alleles: total copies c, minor copies m.
struct AllelicSegment {
  std::uint32_t contig = 0;
  std::size_t begin = 0, end = 0;
  int copies = 2, minor = 1;
  std::size_t sites = 0, heterozygous = 0;
  double quality = 0;
  AllelicState state = AllelicState::heterozygous;
};

// The minor copies of a stretch at total copies c from its heterozygous sites, by maximum likelihood.
inline int minor_copies(const std::vector<AlleleSite>& sites, std::size_t first, std::size_t last, int c, double error = 0.01) {
  if (c < 2) return 0;
  int best = 1;
  double best_l = -INFINITY;
  for (int m = 1; m <= c / 2; m++) {
    const double p = std::clamp(static_cast<double>(m) / c, error, 1 - error);
    double l = 0;
    for (std::size_t i = first; i <= last && i < sites.size(); i++) {
      const auto& s = sites[i];
      if (s.genotype != 1) continue;
      const double a = static_cast<double>(s.alternate), r = static_cast<double>(s.reference);
      const double x = a * std::log(p) + r * std::log(1 - p), y = a * std::log(1 - p) + r * std::log(p);
      l += std::max(x, y) + std::log1p(std::exp(-std::fabs(x - y))) - std::log(2.0);
    }
    if (l > best_l) { best_l = l; best = m; }
  }
  return best;
}

// The runs of homozygosity and the copy-number segments read together.
inline std::vector<AllelicSegment> allelic_segments(const std::vector<AlleleSite>& sites, const std::vector<Run>& runs, const std::vector<CopySegmentCall>& segments) {
  std::vector<AllelicSegment> out;
  // Runs, classified by the copy number over them (the most common over their sites).
  for (const auto& r : runs) {
    AllelicSegment a;
    a.contig = r.contig;
    a.begin = r.begin;
    a.end = r.end;
    a.sites = r.sites;
    a.heterozygous = r.heterozygous;
    a.quality = r.quality;
    a.minor = 0;
    std::vector<std::size_t> by(16, 0);
    for (std::size_t i = r.first; i <= r.last; i++) by[static_cast<std::size_t>(std::min(15, copies_at(segments, sites[i].contig, sites[i].position)))]++;
    a.copies = static_cast<int>(std::max_element(by.begin(), by.end()) - by.begin());
    a.state = a.copies == 2 ? AllelicState::copy_neutral_loh : a.copies == 1 ? AllelicState::hemizygous_deletion : a.copies == 0 ? AllelicState::homozygous_deletion : AllelicState::loh_gain;
    out.push_back(a);
  }
  // Copy-number segments other than the baseline outside runs: their minor copies.
  for (const auto& s : segments) {
    if (s.copies == 2) continue;
    std::size_t first = SIZE_MAX, last = 0, n = 0, het = 0;
    for (std::size_t i = 0; i < sites.size(); i++) {
      if (sites[i].contig != s.contig || sites[i].position < s.begin || sites[i].position >= s.end) continue;
      bool in_run = false;
      for (const auto& r : runs) in_run = in_run || (r.contig == s.contig && sites[i].position >= r.begin && sites[i].position < r.end);
      if (in_run) continue;
      first = std::min(first, i);
      last = i;
      n++;
      het += sites[i].genotype == 1;
    }
    if (n == 0 || s.copies < 2) continue;
    AllelicSegment a;
    a.contig = s.contig;
    a.begin = s.begin;
    a.end = s.end;
    a.copies = s.copies;
    a.sites = n;
    a.heterozygous = het;
    a.minor = minor_copies(sites, first, last, s.copies);
    a.quality = s.quality;
    a.state = 2 * a.minor == a.copies ? AllelicState::heterozygous : AllelicState::imbalance;
    out.push_back(a);
  }
  std::sort(out.begin(), out.end(), [](const AllelicSegment& x, const AllelicSegment& y) { return x.contig != y.contig ? x.contig < y.contig : x.begin < y.begin; });
  return out;
}

}  // namespace genomics::structural
