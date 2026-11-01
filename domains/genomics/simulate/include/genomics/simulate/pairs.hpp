// Paired-end reads simulated from the chromosomes of an individual (structural.hpp), with the
// truth of both mates: a fragment of a length drawn from the library's distribution is cut from
// a chromosome chosen in proportion to the fragments it holds, and read from both ends inwards,
// mate 1 on the fragment's strand and mate 2 on the other (Illumina's forward–reverse, FR,
// orientation: on the reference an unrearranged pair maps as →…←).
//
// - Fragment lengths: normal of the given mean and standard deviation, held to [length + 8,
//   mean + 6 sd] (a fragment shorter than a read, read through into its adapter, is not drawn).
// - GC bias of the library (Benjamini and Speed 2012, Nucleic Acids Research 40, e72: the rate
//   of fragments depends on the GC content of the whole fragment, unimodally): a fragment of GC
//   fraction g is kept with probability (1 − β) + β exp(−½ ((g − g₀)/w)²), else a new one is
//   drawn, up to `max_attempts` times; β = 0 is an unbiased library. A fragment touching an N is
//   never kept (no read comes from an unknown base).
// - Each mate is a short read of the model of reads.hpp (detail::sequence_short: errors rising
//   along the read, a few indels, qualities that track the errors), from the first length + 8
//   bases of its strand of the fragment.
// - Truth: the chromosome, the fragment's interval on it and its strand; for each mate its
//   interval on the chromosome, its strand there, the loci of its two ends on the reference (the
//   chromosome's liftover, which knows inversions, copies, other contigs and novel bases), and
//   whether it crosses a junction (a split read).
//
// Pair i comes from a stream keyed by the seed, PAIRS and i, so that pairs a to b are the same
// made alone or among all of them. Complexity O(count · (length + mean)).

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/simulate/reads.hpp>
#include <genomics/simulate/structural.hpp>

#include <calx/random/stream.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace genomics::simulate {

inline constexpr std::uint32_t PAIRS = 0x50414952;  // "PAIR"

struct FragmentModel {
  double mean = 450, sd = 45;
  double gc_bias = 0;          // β: 0 none, 1 the full curve
  double gc_optimum = 0.45;    // g₀
  double gc_width = 0.12;      // w
  int max_attempts = 64;
};

struct PairModel {
  ShortReadModel read;
  FragmentModel fragment;
};

struct Mate {
  Codes bases;
  std::string quality;
  std::size_t begin = 0, end = 0;  // the chromosome bases it covers
  bool reverse = false;            // read from the reverse strand of the chromosome
  std::size_t errors = 0;
  Locus first, last;               // the loci of chromosome bases begin and end − 1
  bool junction = false;           // it crosses a junction of its chromosome
};

struct ReadPair {
  std::string name;
  std::uint32_t chromosome = 0;
  std::size_t fragment_begin = 0, fragment_length = 0;
  bool reverse = false;            // the fragment read from the chromosome's reverse strand (mate 1 reverse)
  std::array<Mate, 2> mates;
};

// The number of pairs giving `coverage` reads over the sample's chromosomes.
inline std::size_t pairs_for_coverage(const std::vector<Chromosome>& chromosomes, double coverage, std::size_t read_length) {
  std::size_t n = 0;
  for (const auto& c : chromosomes) n += c.sequence.size();
  return static_cast<std::size_t>(std::llround(coverage * static_cast<double>(n) / (2.0 * static_cast<double>(read_length))));
}

// The fraction of G and C among the bases of x[a, b), or −1 when one of them is N.
inline double gc_or_n(const Codes& x, std::size_t a, std::size_t b) {
  std::size_t gc = 0;
  for (std::size_t i = a; i < b; i++) {
    if (x[i] > 3) return -1;
    gc += x[i] == sequence::C || x[i] == sequence::G;
  }
  return b > a ? static_cast<double>(gc) / static_cast<double>(b - a) : 0;
}

// Pairs first, ..., first + count − 1 of a sample.
inline std::vector<ReadPair> read_pairs(const std::vector<Chromosome>& chromosomes, std::size_t count, const PairModel& model, std::uint64_t seed, std::size_t first = 0) {
  std::vector<ReadPair> out;
  out.reserve(count);
  const std::size_t L = model.read.length;
  const auto& fm = model.fragment;
  const double fmin = static_cast<double>(L + 8), fmax = fm.mean + 6 * fm.sd;
  for (std::size_t i = first; i < first + count; i++) {
    calx::random::Stream r(seed, PAIRS, static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(i >> 32));
    ReadPair pair;
    pair.name = "pair" + std::to_string(i);
    bool placed = false;
    std::size_t F = 0, start = 0, which = 0;
    for (int attempt = 0; attempt < std::max(1, fm.max_attempts) && !placed; attempt++) {
      F = static_cast<std::size_t>(std::llround(std::clamp(r.normal(fm.mean, fm.sd), fmin, std::max(fmin, fmax))));
      std::size_t total = 0;
      for (const auto& c : chromosomes) total += c.sequence.size() >= F ? c.sequence.size() - F + 1 : 0;
      if (total == 0) break;
      std::size_t u = r.below(total);
      for (which = 0; which < chromosomes.size(); which++) {
        const std::size_t room = chromosomes[which].sequence.size() >= F ? chromosomes[which].sequence.size() - F + 1 : 0;
        if (u < room) { start = u; break; }
        u -= room;
      }
      const double g = gc_or_n(chromosomes[which].sequence, start, start + F);
      if (g < 0) continue;
      const double keep = (1 - fm.gc_bias) + fm.gc_bias * std::exp(-0.5 * std::pow((g - fm.gc_optimum) / fm.gc_width, 2));
      placed = r.uniform() < keep;
    }
    if (!placed) continue;
    const Chromosome& ch = chromosomes[which];
    pair.chromosome = static_cast<std::uint32_t>(which);
    pair.fragment_begin = start;
    pair.fragment_length = F;
    pair.reverse = r.bernoulli(0.5);
    const Codes forward(ch.sequence.begin() + static_cast<std::ptrdiff_t>(start), ch.sequence.begin() + static_cast<std::ptrdiff_t>(start + F));
    const Codes backward = sequence::reverse_complement(forward);
    const std::size_t span = std::min(F, L + 8);
    for (int m = 0; m < 2; m++) {
      Mate& mate = pair.mates[m];
      // Mate 1 reads the fragment's strand, mate 2 the other.
      mate.reverse = (m == 0) == pair.reverse;
      const Codes& strand = mate.reverse ? backward : forward;
      const Codes t(strand.begin(), strand.begin() + static_cast<std::ptrdiff_t>(span));
      const std::size_t used = detail::sequence_short(r, t, model.read, mate.bases, mate.quality, mate.errors);
      mate.begin = mate.reverse ? start + F - used : start;
      mate.end = mate.reverse ? start + F : start + used;
      mate.first = ch.map.locate(mate.begin);
      mate.last = ch.map.locate(mate.end - 1);
      const auto j = std::upper_bound(ch.junctions.begin(), ch.junctions.end(), mate.begin);
      mate.junction = j != ch.junctions.end() && *j < mate.end;
    }
    out.push_back(std::move(pair));
  }
  return out;
}

}  // namespace genomics::simulate
