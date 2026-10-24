// DNA packed two bits to a base: thirty-two bases in a 64-bit word, base i in bits
// 2(i mod 32) and 2(i mod 32) + 1 of word i / 32, with A = 0, C = 1, G = 2, T = 3
// (alphabet.hpp). The positions that are not A, C, G or T are kept beside the words as
// sorted runs, and read back as N; inside the words they hold A.
//
// A genome of n bases takes n / 4 bytes (Escherichia coli's 4.6 Mb in 1.2 MB), and the
// operations a genome needs run a word at a time:
//
// - the reverse complement: the complement of a word is its bitwise negation (3 − x is x
//   xor 3), and the reversal of the 2-bit groups of a word is five swaps of halves (Warren,
//   Hacker's Delight, 2nd ed. 2013, §7-1), followed by one shift of the whole array by the
//   unused bases of the last word;
// - a k-mer (k ≤ 32) at any position, from at most two words;
// - counts of each base, by popcount on masks of the two bit planes.

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace genomics::sequence {

class PackedDna {
 public:
  PackedDna() = default;
  explicit PackedDna(const Codes& x) : size_(x.size()), words_((x.size() + 31) / 32, 0) {
    for (std::size_t i = 0; i < x.size(); i++) {
      const Code c = x[i];
      if (c < 4) words_[i >> 5] |= static_cast<std::uint64_t>(c) << (2 * (i & 31));
      else if (!runs_.empty() && runs_.back().second == i) runs_.back().second = i + 1;
      else runs_.push_back({i, i + 1});
    }
  }

  std::size_t size() const { return size_; }
  bool empty() const { return size_ == 0; }
  const std::vector<std::uint64_t>& words() const { return words_; }
  // The runs [begin, end) of positions that read as N, in order.
  const std::vector<std::pair<std::size_t, std::size_t>>& ambiguous() const { return runs_; }
  std::size_t bytes() const { return words_.size() * 8 + runs_.size() * sizeof(runs_[0]); }

  bool is_ambiguous(std::size_t i) const {
    auto it = std::upper_bound(runs_.begin(), runs_.end(), std::pair<std::size_t, std::size_t>{i, SIZE_MAX});
    return it != runs_.begin() && std::prev(it)->second > i;
  }
  // The two bits at i, without the N runs.
  Code bits(std::size_t i) const { return static_cast<Code>((words_[i >> 5] >> (2 * (i & 31))) & 3); }
  Code operator[](std::size_t i) const { return is_ambiguous(i) ? N : bits(i); }

  Codes codes(std::size_t begin = 0, std::size_t end = SIZE_MAX) const {
    end = std::min(end, size_);
    Codes out;
    if (begin >= end) return out;
    out.resize(end - begin);
    for (std::size_t i = begin; i < end; i++) out[i - begin] = bits(i);
    auto it = std::upper_bound(runs_.begin(), runs_.end(), std::pair<std::size_t, std::size_t>{begin, SIZE_MAX});
    if (it != runs_.begin()) --it;
    for (; it != runs_.end() && it->first < end; ++it)
      for (std::size_t i = std::max(begin, it->first); i < std::min(end, it->second); i++) out[i - begin] = N;
    return out;
  }

  // The k bases from i (k ≤ 32) as a k-mer, the first base most significant (kmer.hpp's
  // convention); the N runs are not consulted.
  std::uint64_t kmer(std::size_t i, int k) const {
    if (k < 1 || k > 32 || i + static_cast<std::size_t>(k) > size_) throw std::invalid_argument("genomics::sequence::PackedDna::kmer: 1 <= k <= 32, inside the sequence");
    // The 2k bits from bit 2i of the stream, little-endian: base i lowest.
    const std::size_t bit = 2 * i, w = bit >> 6, off = bit & 63;
    std::uint64_t lo = words_[w] >> off;
    if (off && w + 1 < words_.size()) lo |= words_[w + 1] << (64 - off);
    const std::uint64_t mask = k == 32 ? ~0ULL : (1ULL << (2 * k)) - 1;
    lo &= mask;
    // Reverse the order of the k groups, so that base i is the most significant.
    return reverse_groups(lo) >> (64 - 2 * k);
  }

  // The number of A, C, G, T and N, by popcount over the words.
  std::array<std::size_t, 5> counts() const {
    std::array<std::size_t, 5> c{};
    for (std::size_t w = 0; w < words_.size(); w++) {
      const std::size_t bases = std::min<std::size_t>(32, size_ - 32 * w);
      const std::uint64_t valid = bases == 32 ? ~0ULL : (1ULL << (2 * bases)) - 1;
      const std::uint64_t x = words_[w];
      const std::uint64_t lo = x & 0x5555555555555555ULL & valid, hi = (x >> 1) & 0x5555555555555555ULL & valid;
      const std::uint64_t ones = valid & 0x5555555555555555ULL;
      c[T] += std::popcount(lo & hi);
      c[G] += std::popcount(~lo & hi & ones);
      c[C] += std::popcount(lo & ~hi & ones);
      c[A] += std::popcount(~lo & ~hi & ones);
    }
    // The N runs hold A inside the words.
    for (const auto& [b, e] : runs_) {
      c[A] -= e - b;
      c[N] += e - b;
    }
    return c;
  }

  // The reverse complement, a word at a time.
  PackedDna reverse_complement() const {
    PackedDna r;
    r.size_ = size_;
    const std::size_t W = words_.size();
    r.words_.assign(W, 0);
    if (W == 0) return r;
    for (std::size_t w = 0; w < W; w++) r.words_[W - 1 - w] = ~reverse_groups(words_[w]);
    // The array now holds 32W bases, the first `pad` of them from the unused end of the last word.
    const std::size_t pad = 32 * W - size_;
    if (pad) {
      const unsigned s = static_cast<unsigned>(2 * pad);
      for (std::size_t w = 0; w + 1 < W; w++) r.words_[w] = (r.words_[w] >> s) | (r.words_[w + 1] << (64 - s));
      r.words_[W - 1] >>= s;
    }
    // Clear the bits beyond the end, as the constructor leaves them.
    const std::size_t tail = size_ & 31;
    if (tail) r.words_[W - 1] &= (1ULL << (2 * tail)) - 1;
    // The N runs, mirrored; their bits are complemented A, that is T, and are set back to A.
    for (auto it = runs_.rbegin(); it != runs_.rend(); ++it) {
      const std::size_t b = size_ - it->second, e = size_ - it->first;
      r.runs_.push_back({b, e});
      for (std::size_t i = b; i < e; i++) r.words_[i >> 5] &= ~(3ULL << (2 * (i & 31)));
    }
    return r;
  }

  friend bool operator==(const PackedDna& a, const PackedDna& b) { return a.size_ == b.size_ && a.words_ == b.words_ && a.runs_ == b.runs_; }

  // The 2-bit groups of a word in reverse order.
  static std::uint64_t reverse_groups(std::uint64_t x) {
    x = ((x >> 2) & 0x3333333333333333ULL) | ((x & 0x3333333333333333ULL) << 2);
    x = ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((x & 0x0F0F0F0F0F0F0F0FULL) << 4);
    x = ((x >> 8) & 0x00FF00FF00FF00FFULL) | ((x & 0x00FF00FF00FF00FFULL) << 8);
    x = ((x >> 16) & 0x0000FFFF0000FFFFULL) | ((x & 0x0000FFFF0000FFFFULL) << 16);
    return (x >> 32) | (x << 32);
  }

 private:
  std::size_t size_ = 0;
  std::vector<std::uint64_t> words_;
  std::vector<std::pair<std::size_t, std::size_t>> runs_;
};

}  // namespace genomics::sequence
