// Paired reads on a reference: both mates mapped by the domain's mapper (mapping/mapper.hpp), the
// library's distribution of fragment lengths estimated robustly from the pairs that look
// unrearranged, and every pair classified by orientation and distance.
//
// - Mapping. The mates are mapped independently, as single reads, in parallel (each into its own
//   slot): mate m of pair i is read 2i + m. A chimeric mate, one that crosses a junction, comes
//   back with a primary and a supplementary part (the mapper's second primary chain), which is
//   the split read of evidence.hpp.
// - The template of an FR pair on one contig: from the unclipped start of the forward mate to the
//   unclipped end of the reverse mate (soft clips counted, so that a mate clipped at a junction
//   keeps the fragment's length).
// - The library (DELLY; Rausch et al. 2012, Bioinformatics 28, i333–i339, and its implementation,
//   delly src/util.h at 88ea887, getLibraryParams): over pairs with both mates of mapping quality
//   ≥ 20 on one contig in FR orientation, the median m of the template lengths and their median
//   absolute deviation MAD (unscaled); normal pairs lie within m ± 5 MAD (madNormalCutoff), and a
//   pair is discordant by distance beyond m + 9 MAD (DELLY's -s 9), that cut-off at least 500 and
//   twice the read length (DELLY's sanity floors); σ̂ = 1.4826 MAD estimates the standard
//   deviation of a normal law. Pairs beyond 100 kb are left out of the estimate. The survival
//   function S(x) = P(F ≥ x) of the template length F is the empirical one of the normal pairs,
//   which LUMPY's probabilities of breakpoints need (Layer et al. 2014, Genome Biology 15, R84:
//   "D(j), the proportion of fragments longer than j").
// - Classes. With both primaries of mapping quality ≥ `min_mapq`:
//     proper        one contig, FR (→ ←), template within the normal range
//     long          FR, template beyond the cut-off: the signature of a deletion (+, −)
//     short         FR, template below m − 9 MAD (and ≥ 0)
//     everted       one contig, RF (← →): a tandem duplication (−, +)
//     same strand   → → or ← ←: one of an inversion's junctions, (+, +) or (−, −)
//     other contig  a translocation or an insertion of a copy from another contig
//   and with one mate placed and the other unmapped, one-end anchored (the signature of an
//   insertion longer than the reads, or of a sequence absent from the reference). A pair with a
//   mate of low mapping quality is `ambiguous` and is not evidence.
//
// Complexity: O(N) to classify N pairs, O(N log N) to estimate the library.

#pragma once

#include <genomics/mapping/evaluate.hpp>
#include <genomics/mapping/mapper.hpp>
#include <genomics/sequence/alphabet.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace genomics::structural {

using mapping::Mapping;
using sequence::Code;
using sequence::Codes;

// The mappings of 2N mates (mate m of pair i at 2i + m), each read's primary first.
using MateMappings = std::vector<std::vector<Mapping>>;

// Both mates of every pair mapped, in pair order, the same on any number of threads. `Pair` is
// anything with `mates[0..1].bases` (simulate::ReadPair).
template <class Pair>
MateMappings map_pairs(calx::scheduler::ThreadPool& pool, const mapping::Mapper& mapper, const std::vector<Pair>& pairs, std::size_t grain = 32) {
  MateMappings out(2 * pairs.size());
  calx::scheduler::parallel_for(pool, out.size(), grain, [&](std::size_t lo, std::size_t hi, std::size_t) {
    for (std::size_t r = lo; r < hi; r++) out[r] = mapper.map(pairs[r / 2].mates[r % 2].bases, static_cast<std::uint32_t>(r));
  });
  return out;
}

// The unclipped interval of a mapping on the reference: its aligned bases with the soft clips
// carried on.
inline std::int64_t unclipped_begin(const Mapping& m) { return static_cast<std::int64_t>(m.target_begin) - static_cast<std::int64_t>(m.clip_left()); }
inline std::int64_t unclipped_end(const Mapping& m) { return static_cast<std::int64_t>(m.target_end + m.clip_right()); }

// The probability that a placement of mapping quality Q is right, 1 − 10^(−Q/10) (SVTyper's
// prob_mapq, svtyper/utils.py).
inline double placement_probability(int mapq) { return 1 - std::pow(10.0, -static_cast<double>(mapq) / 10.0); }

struct LibraryOptions {
  int min_mapq = 20;
  std::int64_t max_template = 100000;
  double normal_mads = 5, discordant_mads = 9;
  std::int64_t min_cutoff = 500;
  double tail = 1e-3;  // the survival function is kept while above this
};

struct Library {
  double median = 0, mad = 0;
  std::int64_t min_normal = 0, max_normal = 0, min_cutoff = 0, max_cutoff = 0;
  std::size_t read_length = 0;      // the median mate length
  std::size_t pairs = 0;            // pairs the estimate used
  std::vector<double> survival;     // S(x) = P(F ≥ x) for x = 0 .. survival.size() − 1, then 0
  double sd() const { return 1.4826 * mad; }
  double survival_at(std::int64_t x) const {
    if (x <= 0) return 1;
    return static_cast<std::size_t>(x) < survival.size() ? survival[static_cast<std::size_t>(x)] : 0;
  }
  // The longest template the survival function admits.
  std::int64_t longest() const { return static_cast<std::int64_t>(survival.size()) - 1; }
};

namespace detail {

inline const Mapping* placed(const std::vector<Mapping>& ms) { return mapping::primary_of(ms); }

// The template length of an FR pair (forward mate f, reverse mate r) on one contig, or −1.
inline std::int64_t fr_template(const Mapping& f, const Mapping& r) {
  if (f.contig != r.contig || f.reverse || !r.reverse) return -1;
  const std::int64_t b = unclipped_begin(f), e = unclipped_end(r);
  return e > b && b <= unclipped_begin(r) ? e - b : -1;
}

inline double median_of(std::vector<double> v) {
  if (v.empty()) return 0;
  const std::size_t k = v.size() / 2;
  std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(k), v.end());
  return v[k];
}

}  // namespace detail

// The library of a sample of 2N mates' mappings.
inline Library estimate_library(const MateMappings& mates, const LibraryOptions& o = {}) {
  if (mates.size() % 2) throw std::invalid_argument("genomics::structural::estimate_library: two mates for each pair");
  std::vector<double> t;
  std::vector<double> lengths;
  t.reserve(mates.size() / 2);
  for (std::size_t i = 0; i + 1 < mates.size(); i += 2) {
    const Mapping *x = detail::placed(mates[i]), *y = detail::placed(mates[i + 1]);
    if (x) lengths.push_back(static_cast<double>(x->query_length));
    if (!x || !y || x->mapq < o.min_mapq || y->mapq < o.min_mapq) continue;
    const Mapping& f = x->reverse ? *y : *x;
    const Mapping& r = x->reverse ? *x : *y;
    const std::int64_t len = detail::fr_template(f, r);
    if (len > 0 && len <= o.max_template) t.push_back(static_cast<double>(len));
  }
  Library lib;
  lib.read_length = static_cast<std::size_t>(detail::median_of(lengths));
  if (t.empty()) return lib;
  // DELLY's median (the upper middle element) and MAD.
  std::sort(t.begin(), t.end());
  lib.median = t[t.size() / 2];
  std::vector<double> dev(t.size());
  for (std::size_t i = 0; i < t.size(); i++) dev[i] = std::fabs(t[i] - lib.median);
  std::sort(dev.begin(), dev.end());
  lib.mad = std::max(1.0, dev[dev.size() / 2]);
  lib.min_normal = std::max<std::int64_t>(0, static_cast<std::int64_t>(lib.median - o.normal_mads * lib.mad));
  lib.max_normal = static_cast<std::int64_t>(lib.median + o.normal_mads * lib.mad);
  lib.min_cutoff = std::max<std::int64_t>(0, static_cast<std::int64_t>(lib.median - o.discordant_mads * lib.mad));
  lib.max_cutoff = std::max({static_cast<std::int64_t>(lib.median + o.discordant_mads * lib.mad), o.min_cutoff, static_cast<std::int64_t>(2 * lib.read_length)});
  // The empirical survival of the normal templates.
  std::vector<double> normal;
  for (const double x : t)
    if (x >= static_cast<double>(lib.min_normal) && x <= static_cast<double>(lib.max_normal)) normal.push_back(x);
  lib.pairs = normal.size();
  const auto top = static_cast<std::size_t>(lib.max_normal) + 1;
  std::vector<double> count(top + 1, 0);
  for (const double x : normal) count[static_cast<std::size_t>(x)]++;
  lib.survival.assign(top + 1, 0);
  double above = 0;
  for (std::size_t x = top + 1; x-- > 0;) {
    above += count[x];
    lib.survival[x] = above / static_cast<double>(normal.size());
  }
  // Kept while above the tail.
  std::size_t keep = lib.survival.size();
  while (keep > 1 && lib.survival[keep - 1] < o.tail) keep--;
  lib.survival.resize(keep);
  return lib;
}

enum class PairClass : std::uint8_t { unmapped, one_end, ambiguous, proper, long_template, short_template, everted, same_strand, other_contig };

inline const char* pair_class_name(PairClass c) {
  switch (c) {
    case PairClass::unmapped: return "unmapped";
    case PairClass::one_end: return "one-end anchored";
    case PairClass::ambiguous: return "ambiguous";
    case PairClass::proper: return "proper";
    case PairClass::long_template: return "long";
    case PairClass::short_template: return "short";
    case PairClass::everted: return "everted";
    case PairClass::same_strand: return "same strand";
    case PairClass::other_contig: return "other contig";
  }
  return "?";
}

struct PairRecord {
  std::uint32_t pair = 0;
  PairClass type = PairClass::unmapped;
  std::int64_t length = -1;  // the template length of an FR pair on one contig
  bool normal = false;       // a proper pair within the normal range (m ± 5 MAD)
  const Mapping* mate[2]{nullptr, nullptr};  // the primaries (null when unmapped)
  bool discordant() const {
    return type == PairClass::long_template || type == PairClass::everted || type == PairClass::same_strand || type == PairClass::other_contig;
  }
};

struct PairCounts {
  std::size_t by_class[9]{};
  std::size_t of(PairClass c) const { return by_class[static_cast<std::size_t>(c)]; }
};

// The class of every pair.
inline std::vector<PairRecord> classify_pairs(const MateMappings& mates, const Library& lib, int min_mapq = 20, PairCounts* counts = nullptr) {
  if (mates.size() % 2) throw std::invalid_argument("genomics::structural::classify_pairs: two mates for each pair");
  std::vector<PairRecord> out(mates.size() / 2);
  for (std::size_t i = 0; i < out.size(); i++) {
    PairRecord& p = out[i];
    p.pair = static_cast<std::uint32_t>(i);
    const Mapping *x = detail::placed(mates[2 * i]), *y = detail::placed(mates[2 * i + 1]);
    p.mate[0] = x;
    p.mate[1] = y;
    if (!x && !y) p.type = PairClass::unmapped;
    else if (!x || !y) p.type = (x ? x : y)->mapq >= min_mapq ? PairClass::one_end : PairClass::ambiguous;
    else if (x->mapq < min_mapq || y->mapq < min_mapq) p.type = PairClass::ambiguous;
    else if (x->contig != y->contig) p.type = PairClass::other_contig;
    else if (x->reverse == y->reverse) p.type = PairClass::same_strand;
    else {
      const Mapping& f = x->reverse ? *y : *x;
      const Mapping& r = x->reverse ? *x : *y;
      p.length = detail::fr_template(f, r);
      if (p.length < 0) p.type = PairClass::everted;
      else if (p.length > lib.max_cutoff) p.type = PairClass::long_template;
      else if (p.length < lib.min_cutoff) p.type = PairClass::short_template;
      else {
        p.type = PairClass::proper;
        p.normal = p.length >= lib.min_normal && p.length <= lib.max_normal;
      }
    }
    if (counts) counts->by_class[static_cast<std::size_t>(p.type)]++;
  }
  return out;
}

}  // namespace genomics::structural
