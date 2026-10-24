// The accuracy of a mapper against the truth of simulated reads, by mapping quality, and its
// throughput.
//
// - A read is mapped correctly when its primary mapping lies on the contig and strand it came
//   from and overlaps the reference interval it truly covers (simulate/reads.hpp's begin and
//   end). minimap2's own evaluation (paftools.js mapeval) asks for an overlap of a fraction of
//   the intervals, 0.1 by default; any overlap differs from it only for mappings that touch the
//   end of the true interval.
// - The curve of Li 2018, Fig. 1: for each threshold q of mapping quality, the fraction of all
//   reads whose primary mapping has MAPQ ≥ q, and the fraction of those that are wrong. A mapper
//   whose qualities mean what they say has an error rate that falls as q rises, near 10^(−q/10).
// - Throughput: reads and bases a second, from the time the caller measured.
//
// Complexity: O(R) for R reads. The truth is any type with contig, begin, end and reverse
// members (simulate::Read), so that the mapper does not depend on the simulator.

#pragma once

#include <genomics/mapping/mapper.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace genomics::mapping {

struct Truth {
  std::uint32_t contig = 0;
  std::size_t begin = 0, end = 0;
  bool reverse = false;
};

template <class R>
Truth truth_of(const R& r) {
  return {static_cast<std::uint32_t>(r.contig), static_cast<std::size_t>(r.begin), static_cast<std::size_t>(r.end), static_cast<bool>(r.reverse)};
}

inline bool correct(const Mapping& m, const Truth& t) {
  return m.contig == t.contig && m.reverse == t.reverse && m.target_begin < t.end && t.begin < m.target_end;
}

// The primary mapping of a read (the first), or nullptr when it is unmapped.
inline const Mapping* primary_of(const std::vector<Mapping>& mappings) {
  return mappings.empty() || !mappings.front().primary ? nullptr : &mappings.front();
}

struct AccuracyPoint {
  int mapq = 0;
  std::size_t mapped = 0, wrong = 0;     // primaries with MAPQ ≥ mapq, and the wrong among them
  double fraction_mapped = 0, error_rate = 0;
};

struct Accuracy {
  std::size_t reads = 0, mapped = 0, correct = 0, unmapped = 0;
  std::vector<AccuracyPoint> curve;  // thresholds in decreasing order
  double fraction_mapped() const { return reads ? static_cast<double>(mapped) / static_cast<double>(reads) : 0; }
  double fraction_correct() const { return reads ? static_cast<double>(correct) / static_cast<double>(reads) : 0; }
  double error_rate() const { return mapped ? static_cast<double>(mapped - correct) / static_cast<double>(mapped) : 0; }
};

inline Accuracy accuracy(const std::vector<std::vector<Mapping>>& mappings, const std::vector<Truth>& truth,
                         const std::vector<int>& thresholds = {60, 50, 40, 30, 20, 10, 5, 3, 1, 0}) {
  if (mappings.size() != truth.size()) throw std::invalid_argument("genomics::mapping::accuracy: one truth for each read");
  Accuracy a;
  a.reads = mappings.size();
  std::vector<std::size_t> mapped_at(62, 0), wrong_at(62, 0);  // by MAPQ, 0 .. 60
  for (std::size_t i = 0; i < mappings.size(); i++) {
    const Mapping* m = primary_of(mappings[i]);
    if (!m) { a.unmapped++; continue; }
    a.mapped++;
    const bool ok = correct(*m, truth[i]);
    a.correct += ok;
    const auto mq = static_cast<std::size_t>(std::clamp(m->mapq, 0, 60));
    mapped_at[mq]++;
    wrong_at[mq] += !ok;
  }
  for (const int t : thresholds) {
    AccuracyPoint p;
    p.mapq = t;
    for (int v = std::clamp(t, 0, 60); v <= 60; v++) {
      p.mapped += mapped_at[static_cast<std::size_t>(v)];
      p.wrong += wrong_at[static_cast<std::size_t>(v)];
    }
    p.fraction_mapped = a.reads ? static_cast<double>(p.mapped) / static_cast<double>(a.reads) : 0;
    p.error_rate = p.mapped ? static_cast<double>(p.wrong) / static_cast<double>(p.mapped) : 0;
    a.curve.push_back(p);
  }
  return a;
}

// The accuracy against simulated reads (anything truth_of reads).
template <class R>
Accuracy accuracy_of(const std::vector<std::vector<Mapping>>& mappings, const std::vector<R>& reads, const std::vector<int>& thresholds = {60, 50, 40, 30, 20, 10, 5, 3, 1, 0}) {
  std::vector<Truth> t;
  t.reserve(reads.size());
  for (const auto& r : reads) t.push_back(truth_of(r));
  return accuracy(mappings, t, thresholds);
}

struct Throughput {
  std::size_t reads = 0, bases = 0;
  double seconds = 0;
  double reads_per_second() const { return seconds > 0 ? static_cast<double>(reads) / seconds : 0; }
  double bases_per_second() const { return seconds > 0 ? static_cast<double>(bases) / seconds : 0; }
};

}  // namespace genomics::mapping
