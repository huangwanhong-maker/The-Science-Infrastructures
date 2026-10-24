// Genomes simulated from a seed, with what is known about them: random sequence of a given GC
// content, interspersed repeats (copies of a few elements, each copy diverged by substitutions
// at a given rate) and tandem repeats (a unit repeated, with point changes). They are what the
// assembler, the mapper and the indexes are checked on, because their truth is known exactly:
// a repeat longer than a read, or than k, is where reads become ambiguous and de Bruijn graphs
// tangle.
//
// Every draw comes from a calx::random::Stream keyed by the seed and by what it is for
// (the constants below), so that a genome is the same on any machine, thread or shard.

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <calx/random/stream.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace genomics::simulate {

using sequence::Code;
using sequence::Codes;

// The purposes of the streams of this module (the first identity word).
inline constexpr std::uint32_t GENOME = 0x47454e4f, REPEATS = 0x52455053, MUTATIONS = 0x4d555441, READS = 0x52454144, LONG_READS = 0x4c4f4e47,
                               FAMILY = 0x46414d49;

// A base drawn with P(G) = P(C) = gc/2, P(A) = P(T) = (1 − gc)/2.
inline Code base(calx::random::Stream& r, double gc) {
  const double u = r.uniform();
  return u < gc / 2 ? sequence::G : u < gc ? sequence::C : u < (1 + gc) / 2 ? sequence::A : sequence::T;
}

// A base other than b, uniformly.
inline Code substitute(calx::random::Stream& r, Code b) {
  const Code x = static_cast<Code>(r.below(3));
  return b > 3 ? x : static_cast<Code>(x >= b ? x + 1 : x);
}

inline Codes random_sequence(std::size_t length, double gc, std::uint64_t seed, std::uint32_t object = 0) {
  if (gc < 0 || gc > 1) throw std::invalid_argument("genomics::simulate::random_sequence: 0 <= gc <= 1");
  calx::random::Stream r(seed, GENOME, object, 0);
  Codes x(length);
  for (auto& c : x) c = base(r, gc);
  return x;
}

struct RepeatSpec {
  std::size_t elements = 0;    // distinct repeat elements
  std::size_t copies = 0;      // copies of each
  std::size_t length = 0;      // bases of an element
  double divergence = 0;       // the fraction of a copy's bases substituted
};
struct TandemSpec {
  std::size_t arrays = 0;      // tandem arrays
  std::size_t unit = 0;        // bases of the unit
  std::size_t copies = 0;      // units in an array
  double divergence = 0;
};
struct GenomeSpec {
  std::size_t length = 100'000;
  double gc = 0.5;
  RepeatSpec repeats;
  TandemSpec tandem;
};

// Where a copy of a repeat lies.
struct RepeatCopy {
  std::size_t element = 0, begin = 0, end = 0;
  bool reverse = false;   // inserted as the reverse complement
  bool tandem = false;
};

struct Genome {
  Codes sequence;
  std::vector<RepeatCopy> repeats;  // in order of position, never overlapping
};

// A genome of the spec's length: random background, then the copies written over it at
// random positions that do not overlap, each copy on a random strand.
inline Genome make_genome(const GenomeSpec& spec, std::uint64_t seed) {
  Genome g;
  g.sequence = random_sequence(spec.length, spec.gc, seed);
  calx::random::Stream r(seed, REPEATS, 0, 0);
  // The intervals taken, to keep copies apart.
  std::vector<std::pair<std::size_t, std::size_t>> taken;
  const auto place = [&](std::size_t len) -> std::size_t {
    if (len + 2 > spec.length) return SIZE_MAX;
    for (int attempt = 0; attempt < 1000; attempt++) {
      const std::size_t b = r.below(spec.length - len);
      bool free = true;
      for (const auto& [x, y] : taken)
        if (b < y + 1 && x < b + len + 1) { free = false; break; }
      if (free) { taken.push_back({b, b + len}); return b; }
    }
    return SIZE_MAX;
  };
  for (std::size_t e = 0; e < spec.repeats.elements; e++) {
    Codes element(spec.repeats.length);
    for (auto& c : element) c = base(r, spec.gc);
    for (std::size_t c = 0; c < spec.repeats.copies; c++) {
      const std::size_t b = place(element.size());
      if (b == SIZE_MAX) continue;
      const bool rev = r.bernoulli(0.5);
      Codes copy = rev ? sequence::reverse_complement(element) : element;
      for (auto& x : copy)
        if (r.uniform() < spec.repeats.divergence) x = substitute(r, x);
      std::copy(copy.begin(), copy.end(), g.sequence.begin() + static_cast<std::ptrdiff_t>(b));
      g.repeats.push_back({e, b, b + copy.size(), rev, false});
    }
  }
  for (std::size_t a = 0; a < spec.tandem.arrays; a++) {
    Codes unit(spec.tandem.unit);
    for (auto& c : unit) c = base(r, spec.gc);
    const std::size_t len = spec.tandem.unit * spec.tandem.copies;
    const std::size_t b = place(len);
    if (b == SIZE_MAX) continue;
    for (std::size_t i = 0; i < len; i++) {
      Code x = unit[i % unit.size()];
      if (r.uniform() < spec.tandem.divergence) x = substitute(r, x);
      g.sequence[b + i] = x;
    }
    g.repeats.push_back({spec.repeats.elements + a, b, b + len, false, true});
  }
  std::sort(g.repeats.begin(), g.repeats.end(), [](const RepeatCopy& x, const RepeatCopy& y) { return x.begin < y.begin; });
  return g;
}

}  // namespace genomics::simulate
