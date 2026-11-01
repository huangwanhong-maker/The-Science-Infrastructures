// Sequences as the inputs of models of regulation: one-hot matrices, and the counts of k-mers and
// of gapped k-mers.
//
// - One-hot: a sequence of length L as a 4 × L matrix stored base by base (the channels A, C, G, T,
//   each a row of L), 1 where the base is, an N as ¼ in every row (as DeepBind reads it; Alipanahi
//   et al. 2015, Supplementary Notes §1.1); the shape tyche::learning's convolutions take.
// - k-mer counts: every word of length k in [k_min, k_max] counted along the sequence, words with an
//   N skipped; with `canonical`, a word and its reverse complement are one feature (the strand of a
//   regulatory element is unknown).
// - Gapped k-mers (Ghandi, Lee, Mohammad-Noori and Beer 2014, PLoS Computational Biology 10,
//   e1003711, §"Gapped k-mer features"): in every window of length l, the bases at k chosen positions
//   (the informative ones; l − k gaps), so that a feature is a mask of k of the l positions with a
//   word of k bases: C(l, k) 4^k features, which tolerate a mismatch at a gap, as binding sites do.
//   The masks here have their first and last positions informative (a gap at an end is a shorter
//   l-mer's feature), giving C(l − 2, k − 2) masks; with `canonical`, (mask m, word w) and
//   (reverse m, reverse complement w) are one feature. gkm-SVM puts a kernel on these counts and an
//   SVM on top; here they are the columns of a penalised linear model (tyche's elastic net), which
//   is that model's linear analogue, not gkm-SVM itself.
//
// Complexity: one-hot O(L); k-mer counts O(L (k_max − k_min + 1)); gapped k-mers O(L C(l − 2, k − 2)).

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace genomics::regulation {

using sequence::Code;
using sequence::Codes;

// ---- one-hot --------------------------------------------------------------------------------------------
template <class T>
void one_hot_into(const Codes& x, T* out) {
  const std::size_t L = x.size();
  std::fill(out, out + 4 * L, T(0));
  for (std::size_t t = 0; t < L; t++) {
    if (x[t] < 4) out[x[t] * L + t] = T(1);
    else
      for (int b = 0; b < 4; b++) out[b * L + t] = T(0.25);
  }
}
template <class T = double>
std::vector<T> one_hot(const Codes& x) {
  std::vector<T> out(4 * x.size());
  one_hot_into(x, out.data());
  return out;
}

// ---- k-mers -----------------------------------------------------------------------------------------------
namespace detail {

inline std::uint32_t reverse_complement_word(std::uint32_t w, std::size_t k) {
  std::uint32_t r = 0;
  for (std::size_t i = 0; i < k; i++) {
    r = (r << 2) | (3u - (w & 3u));
    w >>= 2;
  }
  return r;
}

inline std::string word_letters(std::uint32_t w, std::size_t k) {
  std::string s(k, 'A');
  for (std::size_t i = k; i-- > 0;) {
    s[i] = "ACGT"[w & 3u];
    w >>= 2;
  }
  return s;
}

}  // namespace detail

class KmerFeatures {
 public:
  KmerFeatures(std::size_t k_min, std::size_t k_max, bool canonical = true) : k_min_(k_min), k_max_(k_max), canonical_(canonical) {
    if (k_min < 1 || k_max < k_min || k_max > 10) throw std::invalid_argument("genomics::regulation::KmerFeatures: 1 ≤ k_min ≤ k_max ≤ 10");
    for (std::size_t k = k_min; k <= k_max; k++) {
      std::vector<std::int32_t> col(std::size_t{1} << (2 * k), -1);
      offsets_.push_back(names_.size());
      for (std::uint32_t w = 0; w < col.size(); w++) {
        const std::uint32_t c = canonical ? std::min(w, detail::reverse_complement_word(w, k)) : w;
        if (col[c] < 0) {
          col[c] = static_cast<std::int32_t>(names_.size() - offsets_.back());
          names_.push_back(detail::word_letters(c, k));
        }
        col[w] = col[c];
      }
      column_.push_back(std::move(col));
    }
  }
  std::size_t columns() const { return names_.size(); }
  const std::vector<std::string>& names() const { return names_; }
  // Adds the counts of x's k-mers into row (of `columns()` values).
  void count(const Codes& x, double* row) const {
    for (std::size_t k = k_min_; k <= k_max_; k++) {
      const auto& col = column_[k - k_min_];
      const std::size_t base = offsets_[k - k_min_];
      const std::uint32_t mask = static_cast<std::uint32_t>((std::size_t{1} << (2 * k)) - 1);
      std::uint32_t w = 0;
      std::size_t valid = 0;
      for (Code c : x) {
        if (c > 3) { valid = 0; continue; }
        w = ((w << 2) | c) & mask;
        if (++valid >= k) row[base + static_cast<std::size_t>(col[w])] += 1;
      }
    }
  }

 private:
  std::size_t k_min_, k_max_;
  bool canonical_;
  std::vector<std::vector<std::int32_t>> column_;
  std::vector<std::size_t> offsets_;
  std::vector<std::string> names_;
};

class GappedKmerFeatures {
 public:
  GappedKmerFeatures(std::size_t l, std::size_t k, bool canonical = true) : l_(l), k_(k) {
    if (k < 2 || l < k || l > 16 || k > 8) throw std::invalid_argument("genomics::regulation::GappedKmerFeatures: 2 ≤ k ≤ l ≤ 16, k ≤ 8");
    // Masks over l positions with the first and last informative.
    for (std::uint32_t m = 0; m < (1u << l); m++)
      if ((m & 1u) && (m >> (l - 1) & 1u) && static_cast<std::size_t>(__builtin_popcount(m)) == k) masks_.push_back(m);
    const std::size_t words = std::size_t{1} << (2 * k);
    column_.assign(masks_.size() * words, -1);
    for (std::size_t mi = 0; mi < masks_.size(); mi++) {
      std::uint32_t rm = 0;  // the mask read backwards
      for (std::size_t i = 0; i < l; i++)
        if (masks_[mi] >> i & 1u) rm |= 1u << (l - 1 - i);
      const std::size_t rmi = static_cast<std::size_t>(std::find(masks_.begin(), masks_.end(), rm) - masks_.begin());
      for (std::uint32_t w = 0; w < words; w++) {
        const std::size_t self = mi * words + w, twin = rmi * words + detail::reverse_complement_word(w, k);
        const std::size_t c = canonical ? std::min(self, twin) : self;
        if (column_[c] < 0) {
          column_[c] = static_cast<std::int32_t>(names_.size());
          names_.push_back(name(masks_[c / words], static_cast<std::uint32_t>(c % words)));
        }
        column_[self] = column_[c];
      }
    }
  }
  std::size_t columns() const { return names_.size(); }
  std::size_t masks() const { return masks_.size(); }
  const std::vector<std::string>& names() const { return names_; }
  void count(const Codes& x, double* row) const {
    if (x.size() < l_) return;
    const std::size_t words = std::size_t{1} << (2 * k_);
    for (std::size_t s = 0; s + l_ <= x.size(); s++)
      for (std::size_t mi = 0; mi < masks_.size(); mi++) {
        std::uint32_t w = 0;
        bool ok = true;
        for (std::size_t i = 0; i < l_ && ok; i++)
          if (masks_[mi] >> i & 1u) {
            const Code c = x[s + i];
            ok = c < 4;
            w = (w << 2) | (c & 3u);
          }
        if (ok) row[static_cast<std::size_t>(column_[mi * words + w])] += 1;
      }
  }

 private:
  std::string name(std::uint32_t mask, std::uint32_t w) const {
    std::string s(l_, 'N');
    std::size_t j = k_;
    for (std::size_t i = l_; i-- > 0;)
      if (mask >> i & 1u) {
        j--;
        s[i] = "ACGT"[(w >> (2 * (k_ - 1 - j))) & 3u];
      }
    return s;
  }
  std::size_t l_, k_;
  std::vector<std::uint32_t> masks_;
  std::vector<std::int32_t> column_;
  std::vector<std::string> names_;
};

}  // namespace genomics::regulation
