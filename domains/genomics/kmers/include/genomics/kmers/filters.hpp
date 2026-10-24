// Approximate membership and approximate counts in fixed memory: the Bloom filter and the
// count-min sketch, keyed by 64-bit words (canonical k-mers, or any key).
//
// - Bloom filter (Bloom 1970, "Space/time trade-offs in hash coding with allowable errors",
//   Comm. ACM 13, 422–426): m bits and h hash functions; a key is inserted by setting its h
//   bits and is reported present when all are set, so there are no false negatives. The h
//   positions come from two hashes by double hashing, gᵢ(x) = a(x) + i b(x) mod m, i < h,
//   which has the same asymptotic false-positive rate as h independent hashes (Kirsch and
//   Mitzenmacher 2006, "Less hashing, same performance: building a better Bloom filter",
//   ESA, LNCS 4168, 456–467). After n insertions a bit is still 0 with
//   probability (1 − 1/m)^{hn} ≈ e^{−hn/m}, so the false-positive rate is
//
//     p ≈ (1 − e^{−hn/m})^h,
//
//   least at h = (m/n) ln 2, where p = 2^{−h} and m = −n ln p / (ln 2)² bits (Broder and
//   Mitzenmacher 2004, Internet Mathematics 1, 485–509, §2.1). The formula is an
//   approximation that slightly underestimates the exact rate at small m (Bose et al. 2008,
//   Inform. Process. Lett. 108, 210–213); the rate of one filter, given the fraction ρ of its
//   bits set, is ρ^h.
// - Count-min sketch (Cormode and Muthukrishnan 2005, "An improved data stream summary: the
//   count-min sketch and its applications", J. Algorithms 55, 58–75): d rows of w counters,
//   row j hashed by a pairwise independent hⱼ; a key's count is added to one counter of each
//   row, and its estimate is the least of them. The estimate never falls below the true count
//   a(x), and with w = ⌈e/ε⌉ and d = ⌈ln(1/δ)⌉
//
//     P[â(x) > a(x) + ε ‖a‖₁] ≤ δ
//
//   (their Theorem 1: each row's excess has mean at most ‖a‖₁/w, Markov's inequality bounds
//   one row by 1/e, and the rows are independent). The hashes are the Carter–Wegman family
//   hⱼ(x) = ((αⱼ x_lo + βⱼ x_hi + γⱼ) mod p) mod w over p = 2^61 − 1 with the key split into
//   32-bit halves (Carter and Wegman 1979, J. Comput. Syst. Sci. 18, 143–154), pairwise
//   independent up to the final reduction mod w, whose bias is below w/p. Counters are 32
//   bits and saturate.
//
// Both are deterministic functions of their seed (one calx::random::Stream each), so that
// shards built with the same seed merge: the union of Bloom filters by bitwise or, the sum of
// sketches counter by counter. Insertion and query cost O(h) and O(d).

#pragma once

#include <genomics/kmers/counter.hpp>

#include <calx/random/stream.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace genomics::kmers {

// The purposes of this header's streams.
inline constexpr std::uint32_t BLOOM = 0x424c4f4d, COUNT_MIN = 0x434d494e;

// (1 − e^{−hn/m})^h, the false-positive rate after n insertions into m bits with h hashes.
inline double bloom_false_positive_rate(double n, double m, double h) { return std::pow(1 - std::exp(-h * n / m), h); }
// The number of hashes that minimises it, (m/n) ln 2, rounded, at least 1.
inline unsigned bloom_optimal_hashes(double n, double m) { return static_cast<unsigned>(std::max(1.0, std::round(m / n * std::numbers::ln2))); }
// The bits that give rate p to n keys at the optimal number of hashes.
inline std::size_t bloom_bits(double n, double p) { return static_cast<std::size_t>(std::ceil(-n * std::log(p) / (std::numbers::ln2 * std::numbers::ln2))); }

class BloomFilter {
 public:
  BloomFilter(std::size_t bits, unsigned hashes, std::uint64_t seed = 0) : m_(bits), h_(hashes), words_((bits + 63) / 64, 0) {
    if (bits == 0 || hashes == 0) throw std::invalid_argument("genomics::kmers::BloomFilter: at least one bit and one hash");
    calx::random::Stream r(seed, BLOOM, 0, 0);
    s1_ = r.next64();
    s2_ = r.next64();
  }
  // A filter for n keys at false-positive rate p.
  static BloomFilter for_capacity(std::size_t n, double p, std::uint64_t seed = 0) {
    if (n == 0 || !(p > 0 && p < 1)) throw std::invalid_argument("genomics::kmers::BloomFilter::for_capacity: n >= 1, 0 < p < 1");
    const std::size_t m = bloom_bits(static_cast<double>(n), p);
    return BloomFilter(m, bloom_optimal_hashes(static_cast<double>(n), static_cast<double>(m)), seed);
  }

  std::size_t bits() const { return m_; }
  unsigned hashes() const { return h_; }
  std::size_t memory_bytes() const { return words_.size() * 8; }
  std::size_t set_bits() const {
    std::size_t c = 0;
    for (auto w : words_) c += static_cast<std::size_t>(std::popcount(w));
    return c;
  }
  double fill() const { return static_cast<double>(set_bits()) / static_cast<double>(m_); }

  void insert(std::uint64_t key) {
    std::uint64_t g, b;
    start(key, g, b);
    for (unsigned i = 0; i < h_; i++) {
      words_[g >> 6] |= 1ULL << (g & 63);
      step(g, b);
    }
  }
  bool contains(std::uint64_t key) const {
    std::uint64_t g, b;
    start(key, g, b);
    for (unsigned i = 0; i < h_; i++) {
      if (!((words_[g >> 6] >> (g & 63)) & 1)) return false;
      step(g, b);
    }
    return true;
  }
  // The union with a filter of the same bits, hashes and seed.
  void merge(const BloomFilter& o) {
    if (o.m_ != m_ || o.h_ != h_ || o.s1_ != s1_ || o.s2_ != s2_) throw std::invalid_argument("genomics::kmers::BloomFilter::merge: filters of one shape and seed");
    for (std::size_t i = 0; i < words_.size(); i++) words_[i] |= o.words_[i];
  }

 private:
  // a(x) and b(x) in [0, m), b ≠ 0, by Lemire's multiply-shift reduction of two mixes.
  void start(std::uint64_t key, std::uint64_t& g, std::uint64_t& b) const {
    const std::uint64_t x = mix64(key ^ s1_), y = mix64(key ^ s2_);
    g = static_cast<std::uint64_t>((static_cast<__uint128_t>(x) * m_) >> 64);
    b = static_cast<std::uint64_t>((static_cast<__uint128_t>(y) * m_) >> 64);
    if (b == 0) b = 1;
  }
  void step(std::uint64_t& g, std::uint64_t b) const {
    g += b;
    if (g >= m_) g -= m_;
  }

  std::size_t m_;
  unsigned h_;
  std::vector<std::uint64_t> words_;
  std::uint64_t s1_ = 0, s2_ = 0;
};

class CountMinSketch {
 public:
  static constexpr std::uint64_t P = (1ULL << 61) - 1;

  CountMinSketch(std::size_t width, std::size_t depth, std::uint64_t seed = 0) : w_(width), d_(depth), cells_(width * depth, 0), rows_(depth) {
    if (width == 0 || depth == 0) throw std::invalid_argument("genomics::kmers::CountMinSketch: width and depth at least 1");
    for (std::size_t j = 0; j < depth; j++) {
      calx::random::Stream r(seed, COUNT_MIN, static_cast<std::uint32_t>(j), 0);
      rows_[j] = {1 + r.below(P - 1), 1 + r.below(P - 1), r.below(P)};
    }
  }
  // w = ⌈e/ε⌉, d = ⌈ln(1/δ)⌉.
  static CountMinSketch for_error(double epsilon, double delta, std::uint64_t seed = 0) {
    if (!(epsilon > 0 && delta > 0 && delta < 1)) throw std::invalid_argument("genomics::kmers::CountMinSketch::for_error: ε > 0, 0 < δ < 1");
    return CountMinSketch(static_cast<std::size_t>(std::ceil(std::numbers::e / epsilon)), static_cast<std::size_t>(std::ceil(std::log(1 / delta))), seed);
  }

  std::size_t width() const { return w_; }
  std::size_t depth() const { return d_; }
  double epsilon() const { return std::numbers::e / static_cast<double>(w_); }
  double delta() const { return std::exp(-static_cast<double>(d_)); }
  std::uint64_t total() const { return total_; }
  std::size_t memory_bytes() const { return cells_.size() * sizeof(std::uint32_t); }

  void add(std::uint64_t key, std::uint32_t count = 1) {
    total_ += count;
    for (std::size_t j = 0; j < d_; j++) {
      std::uint32_t& c = cells_[j * w_ + hash(j, key)];
      c = c > UINT32_MAX - count ? UINT32_MAX : c + count;
    }
  }
  std::uint32_t estimate(std::uint64_t key) const {
    std::uint32_t m = UINT32_MAX;
    for (std::size_t j = 0; j < d_; j++) m = std::min(m, cells_[j * w_ + hash(j, key)]);
    return m;
  }
  // The sum with a sketch of the same shape and seed.
  void merge(const CountMinSketch& o) {
    if (o.w_ != w_ || o.d_ != d_ || o.rows_.front().a != rows_.front().a) throw std::invalid_argument("genomics::kmers::CountMinSketch::merge: sketches of one shape and seed");
    for (std::size_t i = 0; i < cells_.size(); i++) cells_[i] = cells_[i] > UINT32_MAX - o.cells_[i] ? UINT32_MAX : cells_[i] + o.cells_[i];
    total_ += o.total_;
  }

 private:
  struct Row {
    std::uint64_t a, b, c;
  };
  // x mod p for x < 2^122.
  static std::uint64_t reduce(__uint128_t x) {
    std::uint64_t r = static_cast<std::uint64_t>(x & P) + static_cast<std::uint64_t>(x >> 61);
    r = (r & P) + (r >> 61);
    return r >= P ? r - P : r;
  }
  std::size_t hash(std::size_t j, std::uint64_t key) const {
    const Row& r = rows_[j];
    const __uint128_t s = static_cast<__uint128_t>(r.a) * (key & 0xffffffffULL) + static_cast<__uint128_t>(r.b) * (key >> 32) + r.c;
    return static_cast<std::size_t>(reduce(s) % w_);
  }

  std::size_t w_, d_;
  std::vector<std::uint32_t> cells_;
  std::vector<Row> rows_;
  std::uint64_t total_ = 0;
};

}  // namespace genomics::kmers
