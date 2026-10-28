// Strobemers: seeds of n short k-mers ("strobes"), the first at the seed's position and each other
// chosen from a window downstream, so that a seed spans more sequence than its n ℓ bases and
// survives an insertion or deletion between its strobes, where a k-mer of the same span does not
// (Sahlin 2021, "Effective sequence similarity detection with strobemers", Genome Research 31,
// 2080–2094: minstrobes and randstrobes).
//
// - Windows. A strobemer of order n with strobes of ℓ bases and window bounds w_min ≤ w_max starts
//   at position i: strobe 1 is the ℓ-mer at i, and strobe j ≥ 2 is chosen among the ℓ-mers that
//   start in W_j(i) = [i + w_min + (j − 2) w_max, i + (j − 1) w_max], windows of
//   w = w_max − w_min + 1 positions one after another at fixed offsets from i. This is the layout
//   followed here, as Sahlin's paper places the window of order 2. A strobemer exists here only
//   where all its windows fit in the sequence; none is made with shortened windows at its end.
// - Minstrobes: strobe j is the ℓ-mer of least hash in W_j(i), the leftmost on ties, chosen
//   independently of the other strobes. The windows of successive i slide, so one sliding minimum
//   over the ℓ-mers' hashes serves every j: O(L) for a sequence of L bases.
// - Randstrobes: strobe j depends on the strobes before it, through their combined hash H_{j−1}
//   (H_1 = h(strobe 1)): it is the ℓ-mer p of W_j(i) minimising g(H_{j−1}, h(p)), and
//   H_j = g(H_{j−1}, h(strobe j)). Here g(H, h) = mix64(H ⊕ h), SplitMix64's finaliser (a
//   bijection of 64-bit words), so that a randstrobe minimises its own next combined hash; the
//   choice of g is ours, as Sahlin's is a different function with the same purpose: to make the
//   selection of later strobes pseudo-random given the earlier ones. O(L n w).
// - The hash of an ℓ-mer is minimap2's invertible hash64 of the forward ℓ-mer (sequence/kmer.hpp),
//   so distinct ℓ-mers never tie; a seed's hash is H_n. Strobemers are of the forward strand: the
//   reverse strand's are those of the reverse complement. ℓ-mers containing N are never chosen; a
//   strobemer exists at i when the ℓ-mer at i is free of N and every window has an ℓ-mer free of N.
//
// A strobemer of order 2 survives a one-base indel between its first strobe and its window
// (theory, Proposition 3.4): the window then holds the same ℓ-mers but one, which gives way to its
// neighbour outside the window, and the second strobe is the same ℓ-mer unless the one lost was
// the window's choice or the one gained beats it. A k-mer spanning the indel never survives.

#pragma once

#include <genomics/kmers/counter.hpp>
#include <genomics/sequence/kmer.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace genomics::kmers {

using sequence::Codes;

struct StrobemerShape {
  int order = 2;      // n, 2 to 4
  int length = 20;    // ℓ, the bases of each strobe, 1 to 31 (below 32, so no hash is UINT64_MAX)
  int w_min = 25;     // the first window's offset
  int w_max = 50;     // the first window's end; later windows follow every w_max
};

struct Strobemer {
  std::uint64_t hash = 0;              // H_n
  std::uint32_t position = 0;          // i, the first strobe's start
  std::array<std::uint32_t, 4> strobes{};  // the strobes' starts (the first n)
  int order = 0;
  std::uint32_t end(int length) const { return strobes[static_cast<std::size_t>(order - 1)] + static_cast<std::uint32_t>(length); }
};

namespace detail {

inline void check_shape(const StrobemerShape& s) {
  if (s.length < 1 || s.length > 31) throw std::invalid_argument("genomics::kmers::strobemers: strobes of 1 to 31 bases");
  if (s.order < 2 || s.order > 4) throw std::invalid_argument("genomics::kmers::strobemers: an order from 2 to 4");
  if (s.w_min < 1 || s.w_max < s.w_min) throw std::invalid_argument("genomics::kmers::strobemers: 1 <= w_min <= w_max");
}

// h[p] = hash64 of the forward ℓ-mer at p, UINT64_MAX where it holds N or does not fit.
inline std::vector<std::uint64_t> lmer_hashes(const Codes& x, int length) {
  std::vector<std::uint64_t> h(x.size(), UINT64_MAX);
  const sequence::Kmer mask = sequence::kmer_mask(length);
  sequence::for_each_kmer(x, length, [&](std::size_t p, sequence::Kmer f, sequence::Kmer) { h[p] = sequence::hash64(f, mask); });
  return h;
}

inline std::uint64_t combine(std::uint64_t H, std::uint64_t h) { return mix64(H ^ h); }

}  // namespace detail

// The positions i at which a strobemer exists: its last window fits, i + (n − 1) w_max + ℓ ≤ L.
inline std::size_t strobemer_positions(std::size_t L, const StrobemerShape& s) {
  const std::size_t span = static_cast<std::size_t>((s.order - 1) * s.w_max + s.length);
  return L >= span ? L - span + 1 : 0;
}

// The minstrobes of x, in order of position.
inline std::vector<Strobemer> minstrobes(const Codes& x, const StrobemerShape& s) {
  detail::check_shape(s);
  std::vector<Strobemer> out;
  const std::size_t positions = strobemer_positions(x.size(), s);
  if (positions == 0) return out;
  const std::vector<std::uint64_t> h = detail::lmer_hashes(x, s.length);
  // M[p]: the leftmost position of least hash in [p, p + w), by a monotone queue.
  const std::size_t w = static_cast<std::size_t>(s.w_max - s.w_min + 1);
  const std::size_t last = x.size() - static_cast<std::size_t>(s.length);  // the last ℓ-mer that fits
  std::vector<std::uint32_t> M(last + 1 >= w ? last + 2 - w : 0);
  {
    std::vector<std::uint32_t> q(last + 1);
    std::size_t head = 0, tail = 0;
    for (std::size_t p = 0; p <= last; p++) {
      while (tail > head && h[q[tail - 1]] > h[p]) tail--;
      q[tail++] = static_cast<std::uint32_t>(p);
      if (p + 1 >= w) {
        const std::size_t first = p + 1 - w;
        while (q[head] < first) head++;
        M[first] = q[head];
      }
    }
  }
  out.reserve(positions);
  for (std::size_t i = 0; i < positions; i++) {
    if (h[i] == UINT64_MAX) continue;
    Strobemer m;
    m.position = static_cast<std::uint32_t>(i);
    m.order = s.order;
    m.strobes[0] = static_cast<std::uint32_t>(i);
    std::uint64_t H = h[i];
    bool ok = true;
    for (int j = 2; j <= s.order && ok; j++) {
      const std::size_t start = i + static_cast<std::size_t>(s.w_min + (j - 2) * s.w_max);
      const std::uint32_t p = M[start];
      ok = h[p] != UINT64_MAX;
      m.strobes[static_cast<std::size_t>(j - 1)] = p;
      H = detail::combine(H, h[p]);
    }
    if (!ok) continue;
    m.hash = H;
    out.push_back(m);
  }
  return out;
}

// The randstrobes of x, in order of position.
inline std::vector<Strobemer> randstrobes(const Codes& x, const StrobemerShape& s) {
  detail::check_shape(s);
  std::vector<Strobemer> out;
  const std::size_t positions = strobemer_positions(x.size(), s);
  if (positions == 0) return out;
  const std::vector<std::uint64_t> h = detail::lmer_hashes(x, s.length);
  const std::size_t w = static_cast<std::size_t>(s.w_max - s.w_min + 1);
  out.reserve(positions);
  for (std::size_t i = 0; i < positions; i++) {
    if (h[i] == UINT64_MAX) continue;
    Strobemer m;
    m.position = static_cast<std::uint32_t>(i);
    m.order = s.order;
    m.strobes[0] = static_cast<std::uint32_t>(i);
    std::uint64_t H = h[i];
    bool ok = true;
    for (int j = 2; j <= s.order && ok; j++) {
      const std::size_t start = i + static_cast<std::size_t>(s.w_min + (j - 2) * s.w_max);
      std::uint64_t best = UINT64_MAX;
      std::size_t at = start;
      bool any = false;
      for (std::size_t p = start; p < start + w; p++) {
        if (h[p] == UINT64_MAX) continue;
        const std::uint64_t v = detail::combine(H, h[p]);
        if (!any || v < best) { best = v; at = p; any = true; }
      }
      ok = any;
      m.strobes[static_cast<std::size_t>(j - 1)] = static_cast<std::uint32_t>(at);
      H = best;
    }
    if (!ok) continue;
    m.hash = H;
    out.push_back(m);
  }
  return out;
}

}  // namespace genomics::kmers
