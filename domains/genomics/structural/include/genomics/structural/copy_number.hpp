// Copy number as a step function: the depth of each contig (depth.hpp) segmented by tyche's change
// points, an integer copy number for each segment by maximum likelihood, and the aneuploidy of
// whole contigs (theory/structural-variants.md, §4).
//
// - The model. The fragments n_b counted in mappable bin b of a contig are independent Poisson
//   variables of mean (c_b/2) e_b, e_b the bin's exposure at the baseline (depth.hpp), c_b the copy
//   number, a step function of the bins. Its change points are found by PELT (Killick, Fearnhead and
//   Eckley 2012, JASA 107, 1590–1598) on tyche's Poisson cost with exposures, C = 2(Y − Y log(Y/E))
//   for a segment's sums Y = Σ n and E = Σ e (tyche/time_series/changepoint.hpp), with the BIC
//   penalty β = 2 log n per change point by default and segments of at least `min_bins` bins.
// - Integer copy number. On a segment, ℓ(c) = Y log μ_c − μ_c E (up to terms free of c) with
//   μ_c = c/2, and μ_0 = `floor` (the rate of mismapped fragments where nothing is, 0.02 of the
//   baseline); ℓ is concave in μ, so its maximum over the integers 0 … `max_copies` is at ⌊ĉ⌋ or ⌈ĉ⌉
//   of the continuous ĉ = 2Y/E, or at 0, and all are tried. Neighbouring segments of one integer are
//   merged (the segmentation may cut a level the integers cannot tell apart: overdispersion, a
//   bin's error of the GC model). The quality of a segment's copy number is the Phred-scaled
//   likelihood ratio of its best integer against its second, Q = 10 log₁₀ e · (ℓ(ĉ) − ℓ(c₂)).
// - Aneuploidy. A contig whose mappable bins are, for at least `aneuploidy_fraction` (0.9) of
//   them, in segments of one copy number c ≠ the baseline, has c copies: a monosomy (1), trisomy
//   (3), tetrasomy (4).
//
// Complexity: PELT's, expected O(n) for n bins; then O(n). Limits: a segment shorter than
// `min_bins` bins is not seen (a deletion of a few kb is the pairs' and splits' to find).

#pragma once

#include <genomics/structural/depth.hpp>

#include <tyche/time_series/changepoint.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace genomics::structural {

struct CopyNumberOptions {
  tyche::time_series::Penalty penalty = tyche::time_series::Penalty::bic;
  double beta = 0;               // for Penalty::manual
  std::size_t min_bins = 2;
  int max_copies = 8;
  double floor = 0.02;
  int baseline = 2;
  double aneuploidy_fraction = 0.9;
};

struct CopySegmentCall {
  std::uint32_t contig = 0;
  std::size_t begin = 0, end = 0;          // reference bases
  std::size_t first_bin = 0, last_bin = 0; // bins [first, last] of Depth::bins (mappable or not)
  std::size_t bins = 0;                    // mappable bins in it
  double count = 0, exposure = 0;
  int copies = 2;
  double quality = 0;
  double ratio() const { return exposure > 0 ? count / exposure : 0; }
  double estimate() const { return 2 * ratio(); }  // the continuous copy number 2Y/E
};

struct Aneuploidy {
  std::uint32_t contig = 0;
  int copies = 2;
  double fraction = 0;  // of the contig's mappable bins at that copy number
};

namespace detail {

inline double copy_log_likelihood(double Y, double E, int c, double floor) {
  const double mu = c == 0 ? floor : c / 2.0;
  return (Y > 0 ? Y * std::log(mu) : 0) - mu * E;
}

inline void integer_copies(CopySegmentCall& s, const CopyNumberOptions& o) {
  if (s.exposure <= 0) { s.copies = o.baseline; s.quality = 0; return; }
  const double ch = 2 * s.count / s.exposure;
  std::vector<int> candidates{0, static_cast<int>(std::floor(ch)), static_cast<int>(std::ceil(ch)), static_cast<int>(std::floor(ch)) - 1, static_cast<int>(std::ceil(ch)) + 1};
  double best = -INFINITY, second = -INFINITY;
  int arg = o.baseline;
  for (int c : candidates) {
    if (c < 0 || c > o.max_copies) continue;
    const double l = copy_log_likelihood(s.count, s.exposure, c, o.floor);
    if (l > best) { second = best; best = l; arg = c; }
    else if (l > second && c != arg) second = l;
  }
  s.copies = arg;
  s.quality = std::isfinite(second) ? 10 * std::log10(std::exp(1.0)) * (best - second) : 99;
}

}  // namespace detail

// The copy-number segments of every contig, in order, covering each contig's mappable bins.
inline std::vector<CopySegmentCall> segment_copy_number(const Depth& d, const CopyNumberOptions& o = {}) {
  std::vector<CopySegmentCall> out;
  for (std::uint32_t c = 0; c + 1 < d.first.size(); c++) {
    std::vector<std::size_t> idx;
    std::vector<double> y, e;
    for (std::size_t b = d.first[c]; b < d.first[c + 1]; b++)
      if (d.bins[b].mappable && d.bins[b].exposure > 0) {
        idx.push_back(b);
        y.push_back(d.bins[b].count);
        e.push_back(d.bins[b].exposure);
      }
    if (idx.empty()) continue;
    tyche::time_series::ChangepointOptions co;
    co.cost.model = tyche::time_series::SegmentModel::poisson;
    co.cost.exposure = e;
    co.penalty = o.penalty;
    co.beta = o.beta;
    co.min_size = std::min(o.min_bins, idx.size());
    const auto seg = tyche::time_series::pelt(y, co);
    std::vector<CopySegmentCall> mine;
    for (const auto& g : seg.segments) {
      CopySegmentCall s;
      s.contig = c;
      s.first_bin = idx[g.begin];
      s.last_bin = idx[g.end - 1];
      s.bins = g.end - g.begin;
      for (std::size_t k = g.begin; k < g.end; k++) { s.count += y[k]; s.exposure += e[k]; }
      detail::integer_copies(s, o);
      if (!mine.empty() && mine.back().copies == s.copies) {
        CopySegmentCall& m = mine.back();
        m.last_bin = s.last_bin;
        m.bins += s.bins;
        m.count += s.count;
        m.exposure += s.exposure;
        detail::integer_copies(m, o);
        if (m.copies != s.copies) {  // merging moved the estimate across an integer: keep them apart
          m.last_bin = idx[g.begin - 1];
          m.bins -= s.bins;
          m.count -= s.count;
          m.exposure -= s.exposure;
          detail::integer_copies(m, o);
          mine.push_back(s);
        }
      } else {
        mine.push_back(s);
      }
    }
    // Reference coordinates: a segment runs from its first bin to the bin before the next segment.
    for (std::size_t k = 0; k < mine.size(); k++) {
      mine[k].begin = k == 0 ? d.bins[d.first[c]].begin : d.bins[mine[k].first_bin].begin;
      mine[k].end = k + 1 < mine.size() ? d.bins[mine[k + 1].first_bin].begin : d.bins[d.first[c + 1] - 1].end;
    }
    out.insert(out.end(), mine.begin(), mine.end());
  }
  return out;
}

// The contigs at one copy number other than the baseline over (nearly) all their length.
inline std::vector<Aneuploidy> find_aneuploidy(const Depth& d, const std::vector<CopySegmentCall>& segments, const CopyNumberOptions& o = {}) {
  std::vector<Aneuploidy> out;
  for (std::uint32_t c = 0; c + 1 < d.first.size(); c++) {
    std::vector<double> by(static_cast<std::size_t>(o.max_copies) + 1, 0);
    double total = 0;
    for (const auto& s : segments)
      if (s.contig == c) { by[static_cast<std::size_t>(s.copies)] += static_cast<double>(s.bins); total += static_cast<double>(s.bins); }
    if (total == 0) continue;
    const auto best = static_cast<int>(std::max_element(by.begin(), by.end()) - by.begin());
    const double f = by[static_cast<std::size_t>(best)] / total;
    if (best != o.baseline && f >= o.aneuploidy_fraction) out.push_back({c, best, f});
  }
  return out;
}

// The copy number at a reference position (the baseline outside every segment).
inline int copies_at(const std::vector<CopySegmentCall>& segments, std::uint32_t contig, std::size_t p, int baseline = 2) {
  for (const auto& s : segments)
    if (s.contig == contig && p >= s.begin && p < s.end) return s.copies;
  return baseline;
}

}  // namespace genomics::structural
