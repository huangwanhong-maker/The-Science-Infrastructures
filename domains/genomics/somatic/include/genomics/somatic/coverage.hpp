// What a sample's mapped reads say at a place: the mappings sorted by position, the depth of bins,
// and the reads of each base at chosen sites (what ASCAT's alleleCounter reads; Van Loo et al.
// 2010; ascat.prepareHTS).
//
// - MappedReads keeps a sample's primary and supplementary mappings (any mapping quality) sorted by
//   contig and start, with the longest span, so that the mappings overlapping an interval are a
//   contiguous range found by two binary searches.
// - Depth in bins: the mean depth of [b, b + w) is the sum over mappings of mapping quality ≥ q of
//   their overlap with the bin, over w. The tumour's depth over the normal's, each relative to its
//   genome-wide mean, is the depth ratio whose log₂ is ASCAT's logR (copy_number.hpp).
// - Allele counts at sites: the bases (quality ≥ 13 and mapping quality ≥ 20, the pileup's
//   defaults) read at each site, by the variant module's pileup of [p, p + 1), in parallel over the
//   sites, each site's counts in its own slot.
//
// Complexity: O(m log m) to sort m mappings; O(m + bins) for the bins; O(d ℓ) per site for its d
// reads of length ℓ.

#pragma once

#include <genomics/mapping/mapper.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/variant/pileup.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace genomics::somatic {

using sequence::Code;
using sequence::Codes;

class MappedReads {
 public:
  MappedReads() = default;
  MappedReads(const std::vector<std::vector<mapping::Mapping>>& mappings, std::span<const variant::ReadView> reads) : reads_(reads) {
    for (const auto& ms : mappings)
      for (const auto& m : ms)
        if (m.primary) flat_.push_back(m);
    std::sort(flat_.begin(), flat_.end(), [](const mapping::Mapping& a, const mapping::Mapping& b) {
      return std::tie(a.contig, a.target_begin, a.read) < std::tie(b.contig, b.target_begin, b.read);
    });
    for (const auto& m : flat_) longest_ = std::max(longest_, m.target_end - m.target_begin);
  }
  // The mappings that may overlap [b, e) of a contig (every one that does, and a few that end before b).
  std::span<const mapping::Mapping> overlapping(std::uint32_t contig, std::size_t b, std::size_t e) const {
    const auto key = [&](std::size_t pos) {
      return std::lower_bound(flat_.begin(), flat_.end(), std::pair{contig, pos}, [](const mapping::Mapping& m, const std::pair<std::uint32_t, std::size_t>& k) {
        return m.contig != k.first ? m.contig < k.first : m.target_begin < k.second;
      });
    };
    const auto lo = key(b > longest_ ? b - longest_ : 0), hi = key(e);
    return {flat_.data() + (lo - flat_.begin()), static_cast<std::size_t>(hi - lo)};
  }
  std::span<const variant::ReadView> reads() const { return reads_; }
  const std::vector<mapping::Mapping>& mappings() const { return flat_; }
  std::size_t longest() const { return longest_; }

 private:
  std::vector<mapping::Mapping> flat_;
  std::span<const variant::ReadView> reads_;
  std::size_t longest_ = 0;
};

struct DepthBin {
  std::uint32_t contig = 0;
  std::size_t begin = 0, end = 0;
  double depth = 0;  // mean depth
};

// The mean depth of every bin of `width` bases of every contig, from mappings of quality ≥ min_mapq.
inline std::vector<DepthBin> depth_bins(const std::vector<Codes>& contigs, const MappedReads& sample, std::size_t width, int min_mapq = 20) {
  if (width == 0) throw std::invalid_argument("genomics::somatic::depth_bins: bins of at least one base");
  std::vector<std::size_t> first(contigs.size() + 1, 0);
  for (std::size_t c = 0; c < contigs.size(); c++) first[c + 1] = first[c] + (contigs[c].size() + width - 1) / width;
  std::vector<double> bases(first.back(), 0);
  for (const auto& m : sample.mappings()) {
    if (m.contig >= contigs.size() || m.mapq < min_mapq) continue;
    for (std::size_t b = m.target_begin / width; b * width < m.target_end; b++) {
      const std::size_t lo = std::max(m.target_begin, b * width), hi = std::min(m.target_end, (b + 1) * width);
      if (hi > lo) bases[first[m.contig] + b] += static_cast<double>(hi - lo);
    }
  }
  std::vector<DepthBin> out;
  out.reserve(bases.size());
  for (std::uint32_t c = 0; c < contigs.size(); c++)
    for (std::size_t b = 0; first[c] + b < first[c + 1]; b++) {
      const std::size_t lo = b * width, hi = std::min(contigs[c].size(), lo + width);
      out.push_back({c, lo, hi, bases[first[c] + b] / static_cast<double>(hi - lo)});
    }
  return out;
}

struct Site {
  std::uint32_t contig = 0;
  std::size_t position = 0;
  bool operator<(const Site& o) const { return contig != o.contig ? contig < o.contig : position < o.position; }
};

// The reads of A, C, G, T at each site (any order; the counts in the sites' order).
inline std::vector<std::array<std::uint32_t, 4>> allele_counts(calx::scheduler::ThreadPool& pool, const std::vector<Codes>& contigs, std::span<const Site> sites,
                                                               const MappedReads& sample, const variant::PileupOptions& o = {}) {
  std::vector<std::array<std::uint32_t, 4>> out(sites.size(), {0, 0, 0, 0});
  calx::scheduler::parallel_for(pool, sites.size(), 16, [&](std::size_t lo, std::size_t hi, std::size_t) {
    for (std::size_t i = lo; i < hi; i++) {
      const Site& s = sites[i];
      if (s.contig >= contigs.size() || s.position >= contigs[s.contig].size()) continue;
      const auto ms = sample.overlapping(s.contig, s.position, s.position + 1);
      const variant::Pileup pu = variant::pileup(contigs[s.contig], s.contig, s.position, s.position + 1, ms, sample.reads(), o);
      for (const auto& b : pu.at(s.position))
        if (b.base < 4) out[i][b.base]++;
    }
  });
  return out;
}

}  // namespace genomics::somatic
