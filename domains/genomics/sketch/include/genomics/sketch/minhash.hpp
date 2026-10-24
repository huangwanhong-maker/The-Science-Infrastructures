// MinHash sketches of the sets of canonical k-mers of genomes: the Jaccard index, the Mash distance
// and its p-value from bottom-s sketches, and containment and ANI from FracMinHash sketches.
//
// - Bottom-s MinHash (Broder 1997, "On the resemblance and containment of documents",
//   Compression and Complexity of Sequences, 21–29): the s smallest hashes of the set A of
//   canonical k-mers, S(A). For two sets, the s smallest hashes of S(A) ∪ S(B) are the s smallest
//   of A ∪ B, a uniform sample of it, and the fraction of them in both sketches estimates
//   J(A, B) = |A ∩ B| / |A ∪ B| (Ondov et al. 2016, "Mash: fast genome and metagenome distance
//   estimation using MinHash", Genome Biology 17, 132, eq. 3), with the binomial's standard error
//   √(J(1 − J)/s), independent of the genomes' sizes. The sketch is kept as s values bounded by a
//   threshold, refreshed when a buffer of candidates fills: O(n + s log s log(n/s)) expected for
//   n k-mers (Ondov et al., Methods). Sketches of different sizes compare by the smaller.
// - The Mash distance (eq. 4): D = −(1/k) ln(2j/(1 + j)), the mutation rate under a Poisson model
//   in which each k-mer survives with probability e^{−kD} (2j/(1 + j) is the fraction of shared
//   k-mers for sets of equal size). D = 1 when j = 0. Under the simple model of independent
//   substitutions at rate p a k-mer survives with probability (1 − p)^k, so D estimates
//   −ln(1 − p) > p: Mash overestimates the substitution rate by about p²/2 (Hera, Pierce-Ward
//   and Koslicki 2023, Fig. 3A). Ondov et al.'s Table 1 bounds the sketch's error at 0.99
//   probability: for k = 21 and s = 1000, 0.0068, 0.0158, 0.0323, 0.0630 at D = 0.05, 0.10,
//   0.15, 0.20.
// - The p-value (eqs. 1, 5 and 8): a random k-mer lies in a random genome of n distinct k-mers with
//   probability P = 1 − (1 − 1/K)^n, K the number of possible k-mers; two random genomes have the
//   expected Jaccard index r = P_X P_Y / (P_X + P_Y − P_X P_Y); and the chance of x or more shared
//   hashes among s is the binomial tail P(X ≥ x) = I_r(x, s − x + 1) (the regularised incomplete
//   beta function, calx). n is estimated from the sketch as s · 2^64 / v, v its largest hash
//   (Mash, Methods), and is exact when the sketch is not full.
//   *Caveat.* Mash takes K = |Σ|^k = 4^k (eq. 1) though it sketches canonical k-mers, of which
//   there are only K = (4^k + 4^{k/2})/2 (k even) or 4^k/2 (k odd). With 4^k the chance matches
//   expected of unrelated genomes are about half those observed (28 against 63 at k = 12 for two
//   random megabases, in the checks), and the p-value is too small. compare() therefore uses the
//   canonical K, and mash_p_value takes K as an argument; mash_p_value(…, std::pow(4.0, k))
//   reproduces Mash.
// - FracMinHash (Irber et al. 2022, "Lightweight compositional analysis of metagenomes with
//   FracMinHash and minimum metagenome covers", bioRxiv 2022.01.11.475838; Hera et al. 2023,
//   "Deriving confidence intervals for mutation rates across a wide range of evolutionary
//   distances using FracMinHash", Genome Research 33, 1061–1068): every hash at most H·s, with
//   s = 1/scaled, so that the sketch's size grows with the set's, |FRAC_s(A)| ~ Binom(|A|, s), and
//   |A| ≈ |FRAC_s(A)|/s. Containment C(A, B) = |A ∩ B|/|A| is estimated by
//   Ĉ = |FRAC(A) ∩ FRAC(B)| / |FRAC(A)| (Hera et al., eq. 2), whose mean is
//   C(A, B)(1 − (1 − s)^{|A|}) (Theorem 1), so that dividing by that factor debiases it (eq. 3),
//   with variance ≈ m n (1 − s)/(s (m + n)³) for n = |A ∩ B|, m = |A ∖ B| (Theorem 3). Under the
//   simple mutation model E[C] = (1 − p)^k (Theorem 5), and ANI = 1 − p is estimated by C^{1/k}.
//
// Exactness: the sketches are exact functions of the k-mer sets; only the estimates they give
// are random, with the errors stated. Hashes are those of hash.hpp, a bijection of 64-bit words.

#pragma once

#include <genomics/sketch/hash.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/special_functions/gamma.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace genomics::sketch {

// ---- bottom-s MinHash ------------------------------------------------------------------------------------
struct MinHash {
  int k = 21;
  std::size_t size = 1000;       // s, the sketch's capacity
  std::uint64_t seed = 42;
  std::vector<std::uint64_t> hashes;  // the smallest distinct hashes, increasing, at most s
  std::uint64_t kmers = 0;            // k-mers seen, with repeats

  bool full() const { return hashes.size() >= size; }
  // The number of distinct k-mers: exact when the sketch is not full, else s · 2^64 / v.
  double cardinality() const {
    if (!full()) return static_cast<double>(hashes.size());
    return static_cast<double>(size) * 0x1.0p64 / (static_cast<double>(hashes.back()) + 1);
  }
};

namespace detail {

inline void check_sketch(int k, std::size_t size, const char* who) {
  sequence::check_k(k, who);
  if (size == 0) throw std::invalid_argument(std::string("genomics::sketch::") + who + ": a sketch of at least one hash");
}

// Keeps the s smallest distinct values of a stream: candidates below the current threshold are
// buffered, and the buffer is cut back to s whenever it reaches 2s.
class BottomS {
 public:
  explicit BottomS(std::size_t s) : s_(s) { buffer_.reserve(2 * s + 1); }
  void add(std::uint64_t h) {
    if (h >= threshold_) return;
    buffer_.push_back(h);
    if (buffer_.size() >= 2 * s_) trim();
  }
  std::vector<std::uint64_t> take() {
    trim();
    return std::move(buffer_);
  }

 private:
  void trim() {
    std::sort(buffer_.begin(), buffer_.end());
    buffer_.erase(std::unique(buffer_.begin(), buffer_.end()), buffer_.end());
    if (buffer_.size() >= s_) {
      buffer_.resize(s_);
      threshold_ = buffer_.back();  // a value equal to the largest kept is a repeat
    }
  }
  std::size_t s_;
  std::uint64_t threshold_ = std::numeric_limits<std::uint64_t>::max();
  std::vector<std::uint64_t> buffer_;
};

}  // namespace detail

inline MinHash minhash(const std::vector<Codes>& contigs, int k, std::size_t size, std::uint64_t seed = 42) {
  detail::check_sketch(k, size, "minhash");
  MinHash m{k, size, seed, {}, 0};
  detail::BottomS bottom(size);
  for_each_kmer_hash(contigs, k, seed, [&](std::uint64_t h) {
    m.kmers++;
    bottom.add(h);
  });
  m.hashes = bottom.take();
  return m;
}
inline MinHash minhash(const Codes& s, int k, std::size_t size, std::uint64_t seed = 42) {
  return minhash(std::vector<Codes>{s}, k, size, seed);
}

// The sketches of many genomes on the pool, each in its own slot.
inline std::vector<MinHash> minhash_all(calx::scheduler::ThreadPool& pool, const std::vector<std::vector<Codes>>& genomes, int k, std::size_t size,
                                        std::uint64_t seed = 42) {
  std::vector<MinHash> out(genomes.size());
  calx::scheduler::parallel_for(pool, genomes.size(), 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t i = b; i < e; i++) out[i] = minhash(genomes[i], k, size, seed);
  });
  return out;
}

// D = −(1/k) ln(2j/(1 + j)), and 1 when j = 0 (Ondov et al. 2016, eq. 4).
inline double mash_distance(double j, int k) {
  if (!(j > 0)) return 1;
  if (j >= 1) return 0;
  return -std::log(2 * j / (1 + j)) / k;
}

// P(X ≥ x) for X binomial with s trials of probability r: I_r(x, s − x + 1); 1 when x = 0.
inline double binomial_upper_tail(std::size_t x, std::size_t s, double r) {
  if (x == 0) return 1;
  if (x > s) return 0;
  if (!(r > 0)) return 0;
  if (r >= 1) return 1;
  return calx::special_functions::beta_i(static_cast<double>(x), static_cast<double>(s - x + 1), r);
}

// The number of canonical k-mers: (4^k + 4^{k/2})/2 for even k (the palindromes are their own
// reverse complements), 4^k/2 for odd k.
inline double canonical_kmer_space(int k) {
  const double all = std::pow(4.0, k);
  return k % 2 == 0 ? (all + std::pow(4.0, k / 2)) / 2 : all / 2;
}

// The probability that a random k-mer, of K possible, lies in a random genome of n distinct
// k-mers (eq. 1): 1 − (1 − 1/K)^n.
inline double random_kmer_probability(double n, double space) {
  return -std::expm1(n * std::log1p(-1 / space));
}

// The p-value of x shared hashes among s for genomes of na and nb distinct k-mers drawn from a
// space of K (eqs. 1, 5, 8).
inline double mash_p_value(std::size_t x, std::size_t s, double na, double nb, double space) {
  const double px = random_kmer_probability(na, space), py = random_kmer_probability(nb, space);
  const double r = px * py / (px + py - px * py);
  return binomial_upper_tail(x, s, r);
}

struct MashComparison {
  std::size_t shared = 0, considered = 0;  // x and s′ of the merge
  double jaccard = 0, distance = 1, p_value = 1;
};

// The merge of two sketches (eq. 3): the smallest s′ = min(s_a, s_b) distinct hashes of their
// union, and how many of them lie in both. O(s).
inline MashComparison compare(const MinHash& a, const MinHash& b) {
  if (a.k != b.k || a.seed != b.seed) throw std::invalid_argument("genomics::sketch::compare: sketches of the same k and seed");
  const std::size_t s = std::min(a.size, b.size);
  MashComparison c;
  std::size_t i = 0, j = 0;
  while (c.considered < s && (i < a.hashes.size() || j < b.hashes.size())) {
    if (j == b.hashes.size() || (i < a.hashes.size() && a.hashes[i] < b.hashes[j])) i++;
    else if (i == a.hashes.size() || b.hashes[j] < a.hashes[i]) j++;
    else { c.shared++; i++; j++; }
    c.considered++;
  }
  c.jaccard = c.considered ? static_cast<double>(c.shared) / static_cast<double>(c.considered) : 0;
  c.distance = mash_distance(c.jaccard, a.k);
  c.p_value = mash_p_value(c.shared, c.considered, a.cardinality(), b.cardinality(), canonical_kmer_space(a.k));
  return c;
}

// ---- FracMinHash ------------------------------------------------------------------------------------
struct FracMinHash {
  int k = 21;
  std::uint64_t scaled = 1000;   // 1/s
  std::uint64_t seed = 42;
  std::vector<std::uint64_t> hashes;  // every distinct hash at most max_hash(), increasing

  std::uint64_t max_hash() const { return std::numeric_limits<std::uint64_t>::max() / scaled; }
  double fraction() const { return 1.0 / static_cast<double>(scaled); }
  // |A| estimated as |FRAC_s(A)| / s (Hera et al. 2023, Supplement A.6).
  double cardinality() const { return static_cast<double>(hashes.size()) * static_cast<double>(scaled); }
};

inline FracMinHash frac_minhash(const std::vector<Codes>& contigs, int k, std::uint64_t scaled, std::uint64_t seed = 42) {
  sequence::check_k(k, "frac_minhash");
  if (scaled == 0) throw std::invalid_argument("genomics::sketch::frac_minhash: scaled >= 1");
  FracMinHash f{k, scaled, seed, {}};
  const std::uint64_t max = f.max_hash();
  for_each_kmer_hash(contigs, k, seed, [&](std::uint64_t h) {
    if (h <= max) f.hashes.push_back(h);
  });
  std::sort(f.hashes.begin(), f.hashes.end());
  f.hashes.erase(std::unique(f.hashes.begin(), f.hashes.end()), f.hashes.end());
  return f;
}
inline FracMinHash frac_minhash(const Codes& s, int k, std::uint64_t scaled, std::uint64_t seed = 42) {
  return frac_minhash(std::vector<Codes>{s}, k, scaled, seed);
}

// |FRAC(A) ∩ FRAC(B)| by a merge, O(|FRAC(A)| + |FRAC(B)|).
inline std::size_t intersection(const FracMinHash& a, const FracMinHash& b) {
  if (a.k != b.k || a.seed != b.seed || a.scaled != b.scaled) throw std::invalid_argument("genomics::sketch::intersection: sketches of the same k, seed and scale");
  std::size_t i = 0, j = 0, n = 0;
  while (i < a.hashes.size() && j < b.hashes.size()) {
    if (a.hashes[i] < b.hashes[j]) i++;
    else if (b.hashes[j] < a.hashes[i]) j++;
    else { n++; i++; j++; }
  }
  return n;
}

// Ĉ(A, B) = |FRAC(A) ∩ FRAC(B)| / |FRAC(A)| (Hera et al. 2023, eq. 2); 0 for an empty sketch.
inline double containment(const FracMinHash& a, const FracMinHash& b) {
  return a.hashes.empty() ? 0 : static_cast<double>(intersection(a, b)) / static_cast<double>(a.hashes.size());
}
// The debiased containment for |A| = size_a (eq. 3): Ĉ / (1 − (1 − s)^{|A|}).
inline double containment_debiased(const FracMinHash& a, const FracMinHash& b, double size_a) {
  const double factor = -std::expm1(size_a * std::log1p(-a.fraction()));
  return factor > 0 ? containment(a, b) / factor : 0;
}
// |FRAC(A) ∩ FRAC(B)| / |FRAC(A) ∪ FRAC(B)|.
inline double frac_jaccard(const FracMinHash& a, const FracMinHash& b) {
  const std::size_t n = intersection(a, b), u = a.hashes.size() + b.hashes.size() - n;
  return u ? static_cast<double>(n) / static_cast<double>(u) : 0;
}
// ANI = C^{1/k} (Hera et al. 2023, after Theorem 6).
inline double ani(double containment, int k) { return containment > 0 ? std::pow(containment, 1.0 / k) : 0; }
// The variance of Ĉ for n = |A ∩ B| and m = |A ∖ B| at fraction s (Theorem 3, first order).
inline double containment_variance(double n, double m, double s) {
  const double t = m + n;
  return t > 0 ? m * n * (1 - s) / (s * t * t * t) : 0;
}

}  // namespace genomics::sketch
