// k-mers and their samples: the words of length k of a sequence, as integers; the hashes that
// order them; and the subsets that stand for them (minimizers and syncmers), which every
// index, sketch, mapper and assembler of the domain is built on.
//
// - A k-mer, k ≤ 32, is a 64-bit integer with its first base most significant (two bits a
//   base), so that the order of integers is the lexicographic order of the words. Its reverse
//   complement is computed in O(1) word operations; the canonical k-mer is the smaller of the
//   two, so that a k-mer and its reverse complement, which are the same piece of a double
//   strand, have one name. k-mers containing N are skipped.
// - hash64 is the invertible integer hash minimap2 uses (Thomas Wang's 64-bit mix, restricted
//   to 2k bits; Li 2018, Bioinformatics 34, 3094–3100, §2.1.1): a bijection of [0, 4^k), so
//   that distinct k-mers never collide and a random order of them costs nothing.
//   *Remark.* For k ≤ 15 (2k ≤ 31 bits) its last step, key + (key << 31), is the identity
//   modulo 2^{2k}, and the mix is weaker: the density of (w, 15)-minimizers measured on random
//   sequence exceeds 2/(w + 1) by 2 to 3 per cent (w = 5 to 19), where a strong mix (SplitMix64's
//   finaliser) and hash64 at k ≥ 16 both give it within 0.2 per cent (the checks, and a
//   measurement on 2 Mb, 2026-10-24). It is kept for k-mers' names being invertible.
// - ntHash (Mohamadi, Chu, Vandervalk and Birol 2016, Bioinformatics 32, 3492–3494) is a
//   rolling hash of k-mers of any length: h(s_1..s_k) = ⊕_i rol^{k−i} f(s_i), updated in O(1)
//   per base by h' = rol(h) ⊕ rol^k f(s_out) ⊕ f(s_in), with its reverse-complement twin; the
//   canonical value is the smaller of the two (the paper's NTC64).
// - (w, k)-minimizers (Roberts, Hayes, Hunt, Mount and Yorke 2004, Bioinformatics 20,
//   3363–3369; Schleimer, Wilkerson and Aiken 2003, winnowing): of every w consecutive
//   k-mers, the one of smallest hash; ties keep every minimal position, as minimap2 does, and
//   k-mers equal to their reverse complement are skipped, their strand being undefined; a run
//   between Ns shorter than one window gives the minimum of what it has. A
//   random order samples a fraction 2/(w + 1) of the positions (Schleimer et al. 2003, Lemma
//   1), the density the checks measure.
// - Syncmers (Edgar 2021, PeerJ 9, e10805): a k-mer is a closed syncmer when the smallest of
//   its s-mers (by hash) is its first or its last, and an open syncmer when it is at offset t.
//   Unlike minimizers they are chosen by the k-mer alone, not by its neighbours, so a
//   mutation outside a k-mer never changes whether it is chosen. Density 2/(k − s + 1) for
//   closed syncmers, 1/(k − s + 1) for open ones.

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace genomics::sequence {

using Kmer = std::uint64_t;

constexpr Kmer kmer_mask(int k) { return k >= 32 ? ~0ULL : (1ULL << (2 * k)) - 1; }

inline void check_k(int k, const char* who) {
  if (k < 1 || k > 32) throw std::invalid_argument(std::string("genomics::sequence::") + who + ": 1 <= k <= 32");
}

// The k-mer of k codes from p (all below 4).
inline Kmer encode_kmer(const Code* p, int k) {
  Kmer x = 0;
  for (int i = 0; i < k; i++) x = (x << 2) | p[i];
  return x;
}
inline Kmer encode_kmer(std::string_view s) {
  check_k(static_cast<int>(s.size()), "encode_kmer");
  Kmer x = 0;
  for (char c : s) x = (x << 2) | (dna_code(c) & 3);
  return x;
}
inline std::string kmer_string(Kmer x, int k) {
  std::string s(static_cast<std::size_t>(k), 'A');
  for (int i = k - 1; i >= 0; i--, x >>= 2) s[static_cast<std::size_t>(i)] = DNA_LETTERS[x & 3];
  return s;
}

// The reverse complement of a k-mer: reverse the 2-bit groups, complement them, align.
inline Kmer reverse_complement(Kmer x, int k) {
  x = ~x;
  x = ((x >> 2) & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) << 2);
  x = ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((x & 0x0F0F0F0F0F0F0F0FULL) << 4);
  x = ((x >> 8) & 0x00FF00FF00FF00FFULL) | ((x & 0x00FF00FF00FF00FFULL) << 8);
  x = ((x >> 16) & 0x0000FFFF0000FFFFULL) | ((x & 0x0000FFFF0000FFFFULL) << 16);
  x = (x >> 32) | (x << 32);
  return x >> (64 - 2 * k);
}
inline Kmer canonical(Kmer x, int k) {
  const Kmer r = reverse_complement(x, k);
  return x < r ? x : r;
}

// Every k-mer free of N, in order: f(position, forward, reverse complement); positions count from
// the start of s, which may be any stretch of a longer sequence.
template <class F>
void for_each_kmer(std::span<const Code> s, int k, F&& f) {
  check_k(k, "for_each_kmer");
  const Kmer mask = kmer_mask(k);
  const int shift = 2 * (k - 1);
  Kmer fwd = 0, rev = 0;
  int valid = 0;
  for (std::size_t i = 0; i < s.size(); i++) {
    const Code c = s[i];
    if (c > 3) { valid = 0; fwd = rev = 0; continue; }
    fwd = ((fwd << 2) | c) & mask;
    rev = (rev >> 2) | (static_cast<Kmer>(3 - c) << shift);
    if (++valid >= k) f(i + 1 - static_cast<std::size_t>(k), fwd, rev);
  }
}
template <class F>
void for_each_kmer(const Codes& s, int k, F&& f) {
  for_each_kmer(std::span<const Code>(s.data(), s.size()), k, std::forward<F>(f));
}

// minimap2's invertible hash of a key below mask + 1 (mask = 4^k − 1).
constexpr std::uint64_t hash64(std::uint64_t key, std::uint64_t mask) {
  key = (~key + (key << 21)) & mask;
  key = key ^ key >> 24;
  key = ((key + (key << 3)) + (key << 8)) & mask;
  key = key ^ key >> 14;
  key = ((key + (key << 2)) + (key << 4)) & mask;
  key = key ^ key >> 28;
  key = (key + (key << 31)) & mask;
  return key;
}

// ---- ntHash ------------------------------------------------------------------------------------
namespace detail {
inline constexpr std::uint64_t NT_SEED[4] = {0x3c8bfbb395c60474ULL, 0x3193c18562a02b4cULL, 0x20323ed082572324ULL, 0x295549f54be24456ULL};
constexpr std::uint64_t rol(std::uint64_t x, unsigned r) { r &= 63; return r ? (x << r) | (x >> (64 - r)) : x; }
constexpr std::uint64_t ror(std::uint64_t x, unsigned r) { r &= 63; return r ? (x >> r) | (x << (64 - r)) : x; }
}  // namespace detail

// The forward and reverse-complement ntHash of the k bases from p, computed directly.
inline std::pair<std::uint64_t, std::uint64_t> nthash(const Code* p, unsigned k) {
  std::uint64_t f = 0, r = 0;
  for (unsigned i = 0; i < k; i++) {
    f ^= detail::rol(detail::NT_SEED[p[i] & 3], k - 1 - i);
    r ^= detail::rol(detail::NT_SEED[3 - (p[i] & 3)], i);
  }
  return {f, r};
}

// Every k-mer free of N with its ntHash values, rolling: f(position, forward, reverse, canonical).
template <class F>
void for_each_nthash(const Codes& s, unsigned k, F&& f) {
  if (k == 0) throw std::invalid_argument("genomics::sequence::for_each_nthash: k >= 1");
  std::size_t start = 0;  // the first base of the current run free of N
  std::uint64_t fh = 0, rh = 0;
  for (std::size_t i = 0; i < s.size(); i++) {
    if (s[i] > 3) { start = i + 1; continue; }
    if (i + 1 - start < k) continue;
    const std::size_t pos = i + 1 - k;
    if (pos == start) {
      std::tie(fh, rh) = nthash(&s[pos], k);
    } else {
      const Code out = s[pos - 1], in = s[i];
      fh = detail::rol(fh, 1) ^ detail::rol(detail::NT_SEED[out], k) ^ detail::NT_SEED[in];
      rh = detail::ror(rh, 1) ^ detail::ror(detail::NT_SEED[3 - out], 1) ^ detail::rol(detail::NT_SEED[3 - in], k - 1);
    }
    f(pos, fh, rh, fh < rh ? fh : rh);
  }
}

// ---- minimizers ------------------------------------------------------------------------------------
namespace detail {

// The monotone queue of a sliding window, in a ring of fixed size: at most `capacity` items, with
// none of std::deque's allocations (a window of w k-mers holds at most w + 1 candidates at once).
template <class T>
class Ring {
 public:
  explicit Ring(std::size_t capacity) {
    std::size_t n = 1;
    while (n < capacity) n <<= 1;
    buf_.resize(n);
    mask_ = n - 1;
  }
  bool empty() const { return head_ == tail_; }
  std::size_t size() const { return tail_ - head_; }
  T& front() { return buf_[head_ & mask_]; }
  T& back() { return buf_[(tail_ - 1) & mask_]; }
  const T& operator[](std::size_t i) const { return buf_[(head_ + i) & mask_]; }
  void push_back(const T& x) { buf_[tail_++ & mask_] = x; }
  void pop_back() { --tail_; }
  void pop_front() { ++head_; }
  void clear() { head_ = tail_ = 0; }

 private:
  std::vector<T> buf_;
  std::size_t mask_ = 0, head_ = 0, tail_ = 0;
};

}  // namespace detail

struct Minimizer {
  std::uint64_t hash = 0;     // hash64 of the canonical k-mer
  std::uint32_t position = 0; // the k-mer's first base
  bool reverse = false;       // the canonical k-mer is the reverse complement of the forward one
  friend bool operator==(const Minimizer&, const Minimizer&) = default;
};

// The (w, k)-minimizers of a sequence, in order of position, each once.
inline std::vector<Minimizer> minimizers(const Codes& s, int k, int w) {
  check_k(k, "minimizers");
  if (w < 1) throw std::invalid_argument("genomics::sequence::minimizers: w >= 1");
  const Kmer mask = kmer_mask(k);
  std::vector<Minimizer> out;
  out.reserve(2 * s.size() / static_cast<std::size_t>(w + 1) + 16);
  // The candidates of the current window, by increasing position, with non-decreasing hashes.
  detail::Ring<Minimizer> q(static_cast<std::size_t>(w) + 1);
  std::size_t seen = 0;        // k-mers in the current run free of N
  std::int64_t last = -1;      // the last position written
  const auto emit_window = [&] {
    const std::uint64_t m = q.front().hash;
    for (std::size_t j = 0; j < q.size(); j++) {
      const Minimizer& c = q[j];
      if (c.hash != m) break;
      if (static_cast<std::int64_t>(c.position) > last) { out.push_back(c); last = c.position; }
    }
  };
  Kmer fwd = 0, rev = 0;
  int valid = 0;
  for (std::size_t i = 0; i < s.size(); i++) {
    const Code c = s[i];
    if (c > 3) {
      // A window may not span an N: the run ends, with its last partial window if it had no full one.
      if (seen > 0 && seen < static_cast<std::size_t>(w) && !q.empty()) emit_window();
      q.clear();
      seen = 0;
      valid = 0;
      fwd = rev = 0;
      continue;
    }
    fwd = ((fwd << 2) | c) & mask;
    rev = (rev >> 2) | (static_cast<Kmer>(3 - c) << (2 * (k - 1)));
    if (++valid < k) continue;
    const auto pos = static_cast<std::uint32_t>(i + 1 - static_cast<std::size_t>(k));
    seen++;
    if (fwd != rev) {
      const Minimizer m{hash64(fwd < rev ? fwd : rev, mask), pos, rev < fwd};
      while (!q.empty() && q.back().hash > m.hash) q.pop_back();
      q.push_back(m);
    }
    // Drop candidates that have left the window [pos − w + 1, pos].
    while (!q.empty() && q.front().position + static_cast<std::uint32_t>(w) <= pos) q.pop_front();
    if (seen >= static_cast<std::size_t>(w) && !q.empty()) emit_window();
  }
  if (seen > 0 && seen < static_cast<std::size_t>(w) && !q.empty()) emit_window();
  return out;
}

// The positions of the closed (t < 0) or open (offset t) syncmers among the k-mers, by the
// hash64 order of canonical s-mers.
inline std::vector<std::uint32_t> syncmers(const Codes& seq, int k, int s, int t = -1) {
  check_k(k, "syncmers");
  check_k(s, "syncmers");
  if (s > k || t >= k - s + 1) throw std::invalid_argument("genomics::sequence::syncmers: s <= k and t <= k - s");
  // The hash of every s-mer free of N (UINT64_MAX where there is none).
  std::vector<std::uint64_t> h(seq.size(), UINT64_MAX);
  const Kmer mask = kmer_mask(s);
  for_each_kmer(seq, s, [&](std::size_t p, Kmer f, Kmer r) { h[p] = hash64(f < r ? f : r, mask); });
  std::vector<std::uint32_t> out;
  const int m = k - s + 1;  // s-mers in a k-mer
  detail::Ring<std::size_t> q(static_cast<std::size_t>(m) + 1);  // positions of s-mers, increasing hash, in the window of m
  std::size_t run = 0;
  for (std::size_t j = 0; j + static_cast<std::size_t>(s) <= seq.size(); j++) {
    if (h[j] == UINT64_MAX) { q.clear(); run = 0; continue; }
    while (!q.empty() && h[q.back()] > h[j]) q.pop_back();
    q.push_back(j);
    if (++run < static_cast<std::size_t>(m)) continue;
    const std::size_t first = j + 1 - static_cast<std::size_t>(m);  // the k-mer's start
    while (q.front() < first) q.pop_front();
    const std::size_t at = q.front() - first;  // the first minimal s-mer's offset
    const bool chosen = t < 0 ? (at == 0 || at == static_cast<std::size_t>(m - 1)) : at == static_cast<std::size_t>(t);
    if (chosen) out.push_back(static_cast<std::uint32_t>(first));
  }
  return out;
}

}  // namespace genomics::sequence
