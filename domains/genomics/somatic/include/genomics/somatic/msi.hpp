// Microsatellite instability from a tumour and its normal, as MSIsensor measures it (Niu et al.
// 2014, "MSIsensor: microsatellite instability detection using paired tumor-normal sequence
// data", Bioinformatics 30, 1015–1016; the program ding-lab/msisensor, whose scan, homo.cpp,
// chi.cpp and sample.cpp give the rules followed here).
//
// - Sites (find_microsatellites): runs of a primitive unit of u = 1 to 5 bases repeated, maximal
//   on both sides; kept for the analysis when a homopolymer has 10 to 50 copies (MSIsensor's
//   -p and -m) and a longer unit 5 to 40 (-s and -w), with 5 bases of flank on each side (its
//   ContextLength) free of N and of other sites. A unit is primitive when it is no power of a
//   shorter word (ATAT is (AT)²); a run of a non-primitive unit is the run of its root.
// - Lengths read (read_lengths): a read spanning a site, as mapped (its reverse complement when
//   mapped to the reverse strand), contributes the number of whole units between an occurrence of
//   the left flank and the right flank that follows them at once (homo.cpp's scan of the front
//   and end k-mers); reads in which no such reading exists contribute nothing. Mappings of
//   quality ≥ 20, primaries and supplementaries.
// - The test (chi.cpp, X2BetweenTwo): Pearson's χ² of the 2 × K table of the normal's and the
//   tumour's counts over the K lengths read in either, with K − 1 degrees of freedom, p = Q((K −
//   1)/2, χ²/2) (calx's regularised upper incomplete gamma; MSIsensor sums the series of the lower
//   one to 200 terms, the same quantity to its rounding); p = 1 when K = 1 or χ² = 0. A site is
//   tested when both samples have at least `min_coverage` such reads (20 for exomes, 15 for
//   genomes: its -c).
// - Somatic sites (sample.cpp, calculateFDR): the N tested sites by increasing p, the site of rank
//   k is somatic when p_k N / k ≤ 0.05. *Remark.* This is not Benjamini and Hochberg's step-up
//   procedure, which also rejects every site ranked before the last one passing: MSIsensor
//   rejects a subset of BH's rejections. The rule is MSIsensor's, so that its score is reproduced.
// - The score: the percentage of tested sites that are somatic; MSI-high at ≥ 3.5% (the cut-off
//   of the paper's endometrial validation: 70 of 71 MSI tumours above it, 165 of 168 stable ones
//   below).
//
// Complexity: O(n · 5) for the scan of n bases; O(d (ℓ + s)) for a site of span s read by d reads
// of length ℓ. Limits: MSIsensor's coverage normalisation (-z) is off, as its default; no
// tumour-only mode (MSIsensor's entropy, MSIsensor-pro's polymerase slippage model); a read is
// read as sequenced, without its qualities.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/somatic/coverage.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/special_functions/gamma.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <span>
#include <string>
#include <vector>

namespace genomics::somatic {

struct MicrosatelliteOptions {
  std::size_t min_homopolymer = 10, max_homopolymer = 50;  // copies of a one-base unit
  std::size_t min_repeats = 5, max_repeats = 40;           // copies of a unit of 2–5 bases
  std::size_t max_unit = 5;
  std::size_t flank = 5;
};

struct Microsatellite {
  std::uint32_t contig = 0;
  std::size_t begin = 0;     // the first base of the repeat
  std::size_t unit = 1, copies = 0;
  Codes motif;               // the unit
  Codes left, right;         // the flanks
  std::size_t end() const { return begin + unit * copies; }
};

namespace detail {

// Whether w[a, a + u) is a power of a shorter word.
inline bool primitive(const Codes& x, std::size_t a, std::size_t u) {
  for (std::size_t d = 1; d < u; d++) {
    if (u % d) continue;
    bool power = true;
    for (std::size_t i = d; i < u && power; i++) power = x[a + i] == x[a + i % d];
    if (power) return false;
  }
  return true;
}

}  // namespace detail

// The microsatellites of the contigs, by contig and position.
inline std::vector<Microsatellite> find_microsatellites(const std::vector<Codes>& contigs, const MicrosatelliteOptions& o = {}) {
  std::vector<Microsatellite> out;
  for (std::uint32_t c = 0; c < contigs.size(); c++) {
    const Codes& x = contigs[c];
    const std::size_t n = x.size();
    std::vector<Microsatellite> found;
    for (std::size_t u = 1; u <= o.max_unit; u++) {
      // Maximal runs with x[i] = x[i + u]: a run [a, b) of agreement makes the repeat [a, b + u).
      std::size_t i = 0;
      while (i + u < n) {
        if (x[i] > 3 || x[i] != x[i + u]) { i++; continue; }
        std::size_t j = i;
        while (j + u < n && x[j] <= 3 && x[j] == x[j + u]) j++;
        const std::size_t begin = i, span = j - i + u, copies = span / u;
        const std::size_t lo = u == 1 ? o.min_homopolymer : o.min_repeats, hi = u == 1 ? o.max_homopolymer : o.max_repeats;
        if (copies >= lo && copies <= hi && detail::primitive(x, begin, u)) {
          Microsatellite m;
          m.contig = c;
          m.begin = begin;
          m.unit = u;
          m.copies = copies;
          m.motif.assign(x.begin() + static_cast<std::ptrdiff_t>(begin), x.begin() + static_cast<std::ptrdiff_t>(begin + u));
          found.push_back(std::move(m));
        }
        i = j + 1;
      }
    }
    std::sort(found.begin(), found.end(), [](const Microsatellite& a, const Microsatellite& b) { return a.begin != b.begin ? a.begin < b.begin : a.unit < b.unit; });
    // Flanks free of N and of other sites (a site's whole span with its flanks); a site spans at
    // most `reach` bases, so only the sites starting within it of each other can meet.
    const std::size_t reach = std::max(o.max_homopolymer, o.max_repeats * o.max_unit) + o.flank;
    for (std::size_t k = 0; k < found.size(); k++) {
      Microsatellite& m = found[k];
      if (m.begin < o.flank || m.end() + o.flank > n) continue;
      bool clear = true;
      for (std::size_t q = m.begin - o.flank; q < m.end() + o.flank && clear; q++) clear = x[q] <= 3;
      const auto meets = [&](std::size_t q) { return found[q].begin < m.end() + o.flank && m.begin < found[q].end() + o.flank; };
      for (std::size_t q = k + 1; q < found.size() && found[q].begin < m.end() + o.flank && clear; q++) clear = !meets(q);
      for (std::size_t q = k; q-- > 0 && found[q].begin + reach > m.begin && clear;) clear = !meets(q);
      if (!clear) continue;
      m.left.assign(x.begin() + static_cast<std::ptrdiff_t>(m.begin - o.flank), x.begin() + static_cast<std::ptrdiff_t>(m.begin));
      m.right.assign(x.begin() + static_cast<std::ptrdiff_t>(m.end()), x.begin() + static_cast<std::ptrdiff_t>(m.end() + o.flank));
      out.push_back(m);
    }
  }
  return out;
}

// The copies of the site's unit a read shows between the flanks, or −1.
inline long long copies_in_read(std::span<const Code> read, const Microsatellite& m) {
  const std::size_t f = m.left.size(), u = m.unit;
  if (read.size() < 2 * f + u) return -1;
  for (std::size_t s = 0; s + f <= read.size(); s++) {
    if (!std::equal(m.left.begin(), m.left.end(), read.begin() + static_cast<std::ptrdiff_t>(s))) continue;
    std::size_t at = s + f, k = 0;
    while (at + u <= read.size() && std::equal(m.motif.begin(), m.motif.end(), read.begin() + static_cast<std::ptrdiff_t>(at))) { at += u; k++; }
    if (at + m.right.size() <= read.size() && std::equal(m.right.begin(), m.right.end(), read.begin() + static_cast<std::ptrdiff_t>(at))) return static_cast<long long>(k);
  }
  return -1;
}

// The distribution of lengths (copies → reads) read at each site by a sample.
inline std::vector<std::map<std::size_t, std::size_t>> read_lengths(calx::scheduler::ThreadPool& pool, const std::vector<Microsatellite>& sites, const MappedReads& sample,
                                                                    int min_mapq = 20) {
  std::vector<std::map<std::size_t, std::size_t>> out(sites.size());
  calx::scheduler::parallel_for(pool, sites.size(), 8, [&](std::size_t lo, std::size_t hi, std::size_t) {
    Codes rc;
    for (std::size_t i = lo; i < hi; i++) {
      const Microsatellite& m = sites[i];
      const std::size_t a = m.begin - m.left.size(), b = m.end() + m.right.size();
      for (const auto& x : sample.overlapping(m.contig, a, b)) {
        if (x.mapq < min_mapq || x.target_begin > a || x.target_end < b || x.read >= sample.reads().size()) continue;
        const auto bases = sample.reads()[x.read].bases;
        long long k;
        if (x.reverse) {
          rc.assign(bases.begin(), bases.end());
          rc = sequence::reverse_complement(rc);
          k = copies_in_read(rc, m);
        } else {
          k = copies_in_read(bases, m);
        }
        if (k >= 0) out[i][static_cast<std::size_t>(k)]++;
      }
    }
  });
  return out;
}

struct ChiSquare {
  double statistic = 0, p = 1;
  std::size_t dof = 0;
};

// Pearson's χ² of two distributions over the same lengths (MSIsensor's X2BetweenTwo).
inline ChiSquare chi_square_two(const std::map<std::size_t, std::size_t>& first, const std::map<std::size_t, std::size_t>& second) {
  std::map<std::size_t, std::pair<double, double>> table;
  double s1 = 0, s2 = 0;
  for (const auto& [k, v] : first) { table[k].first += static_cast<double>(v); s1 += static_cast<double>(v); }
  for (const auto& [k, v] : second) { table[k].second += static_cast<double>(v); s2 += static_cast<double>(v); }
  ChiSquare r;
  const double s = s1 + s2;
  std::size_t K = 0;
  for (const auto& [k, v] : table) {
    if (v.first + v.second <= 0) continue;
    K++;
    const double e1 = (v.first + v.second) * s1 / s, e2 = (v.first + v.second) * s2 / s;
    if (e1 > 0) r.statistic += (v.first - e1) * (v.first - e1) / e1;
    if (e2 > 0) r.statistic += (v.second - e2) * (v.second - e2) / e2;
  }
  if (K <= 1 || r.statistic == 0) return {r.statistic, 1, K ? K - 1 : 0};
  r.dof = K - 1;
  r.p = calx::special_functions::chi_square_q(r.statistic, static_cast<double>(r.dof));
  return r;
}

struct MsiSite {
  std::size_t site = 0;                          // in the list of microsatellites
  std::map<std::size_t, std::size_t> normal, tumour;
  std::size_t normal_reads = 0, tumour_reads = 0;
  bool tested = false;
  ChiSquare test;
  double fdr = 1;
  bool somatic = false;
};

struct MsiOptions {
  std::size_t min_coverage = 15;  // MSIsensor's -c for genomes (20 for exomes)
  double fdr = 0.05;
  double high = 3.5;              // MSI-high at a score of at least this (%)
  int min_mapq = 20;
};

struct MsiResult {
  std::vector<MsiSite> sites;
  std::size_t tested = 0, somatic = 0;
  double score = 0;   // %
  bool high = false;
};

// MSIsensor's rule: the sites by increasing p; rank k somatic when p_k · N / k ≤ α.
inline void msisensor_fdr(std::vector<MsiSite>& sites, double alpha) {
  std::vector<std::size_t> order;
  for (std::size_t i = 0; i < sites.size(); i++)
    if (sites[i].tested) order.push_back(i);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return sites[a].test.p < sites[b].test.p; });
  const double N = static_cast<double>(order.size());
  for (std::size_t r = 0; r < order.size(); r++) {
    MsiSite& s = sites[order[r]];
    s.fdr = s.test.p * N / static_cast<double>(r + 1);
    s.somatic = s.fdr <= alpha;
  }
}

inline MsiResult msi(calx::scheduler::ThreadPool& pool, const std::vector<Microsatellite>& sites, const MappedReads& tumour, const MappedReads& normal, const MsiOptions& o = {}) {
  MsiResult r;
  const auto t = read_lengths(pool, sites, tumour, o.min_mapq), n = read_lengths(pool, sites, normal, o.min_mapq);
  r.sites.resize(sites.size());
  for (std::size_t i = 0; i < sites.size(); i++) {
    MsiSite& s = r.sites[i];
    s.site = i;
    s.normal = n[i];
    s.tumour = t[i];
    for (const auto& [k, v] : s.normal) s.normal_reads += v;
    for (const auto& [k, v] : s.tumour) s.tumour_reads += v;
    s.tested = s.normal_reads >= o.min_coverage && s.tumour_reads >= o.min_coverage;
    if (s.tested) {
      s.test = chi_square_two(s.normal, s.tumour);
      r.tested++;
    }
  }
  msisensor_fdr(r.sites, o.fdr);
  for (const auto& s : r.sites) r.somatic += s.somatic;
  r.score = r.tested ? 100.0 * static_cast<double>(r.somatic) / static_cast<double>(r.tested) : 0;
  r.high = r.score >= o.high;
  return r;
}

}  // namespace genomics::somatic
