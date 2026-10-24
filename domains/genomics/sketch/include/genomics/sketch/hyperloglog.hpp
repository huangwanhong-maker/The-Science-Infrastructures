// HyperLogLog: the number of distinct k-mers of a genome, and of a union of genomes, in 2^p bytes,
// with a relative standard error of about 1.04/√m for m = 2^p registers (Flajolet, Fusy,
// Gandouet and Meunier 2007, "HyperLogLog: the analysis of a near-optimal cardinality estimation
// algorithm", AofA, DMTCS Proceedings AH, 127–146).
//
// - The sketch. A 64-bit hash h is split into its top p bits, the register j, and the remaining
//   q = 64 − p bits w; the register keeps the largest ρ(w) seen, the position of the first one bit
//   of w (q + 1 when w = 0). Inserting is O(1); inserting the same element twice changes nothing;
//   the sketch of a union is the register-wise maximum of the sketches, exactly (a join
//   semilattice homomorphism from sets under ∪ to registers under max).
// - The classical estimate (Flajolet et al., Fig. 3): E = α_m m² / Σ_j 2^{−M_j}, with
//   α_m = 0.7213/(1 + 1.079/m) (and the constants for m = 16, 32, 64), replaced by linear
//   counting m ln(m/V) when E ≤ 5m/2 and V registers are zero. With a 64-bit hash no large-range
//   correction is needed (Heule, Nunkesser and Hall 2013, "HyperLogLog in practice", EDBT,
//   683–692, §5.2, the improvement adopted here); their empirical bias tables are not, as the
//   next estimator makes them unnecessary.
// - The improved estimate (Ertl 2017, "New cardinality estimation algorithms for HyperLogLog
//   sketches", arXiv:1702.01284, eq. 32 and Algorithm 6), the default:
//     n̂ = α_∞ m² / (m σ(C₀/m) + Σ_{k=1}^{q} C_k 2^{−k} + m τ(1 − C_{q+1}/m) 2^{−q}),
//   with C_k the number of registers equal to k, α_∞ = 1/(2 ln 2),
//   σ(x) = x + Σ_{k≥1} x^{2^k} 2^{k−1} and τ(x) = (1/3)(1 − x − Σ_{k≥1} (1 − x^{2^{−k}})² 2^{−k}).
//   It is nearly unbiased over the whole range, from empty to 2^64, with no switch between
//   estimators and no tables; σ and τ converge to double precision in at most a few dozen steps.
// - Jaccard by inclusion and exclusion (Baker and Langmead 2019, "Dashing: fast and accurate
//   genomic distances with HyperLogLog", Genome Biology 20, 265, Methods): |A ∩ B| = |A| + |B| −
//   |A ∪ B|, each term estimated, so J = (n̂_A + n̂_B − n̂_{A∪B}) / n̂_{A∪B}, clamped to [0, 1]. Its
//   error is that of three estimates, larger than MinHash's for small J.
//
// The estimates assume the hashes of distinct elements are independent and uniform (hash.hpp).
// p from 4 to 18.

#pragma once

#include <genomics/sketch/hash.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace genomics::sketch {

class HyperLogLog {
 public:
  explicit HyperLogLog(int p = 12, std::uint64_t seed = 42) : p_(p), seed_(seed) {
    if (p < 4 || p > 18) throw std::invalid_argument("genomics::sketch::HyperLogLog: 4 <= p <= 18");
    registers_.assign(std::size_t{1} << p, 0);
  }

  int precision() const { return p_; }
  std::uint64_t seed() const { return seed_; }
  std::size_t size() const { return registers_.size(); }
  const std::vector<std::uint8_t>& registers() const { return registers_; }
  // The relative standard error of the estimate, 1.04/√m.
  double standard_error() const { return 1.04 / std::sqrt(static_cast<double>(registers_.size())); }

  void insert(std::uint64_t hash) {
    const std::size_t j = static_cast<std::size_t>(hash >> (64 - p_));
    const std::uint64_t w = hash << p_;
    const int q = 64 - p_;
    const auto rho = static_cast<std::uint8_t>(w == 0 ? q + 1 : std::min(std::countl_zero(w), q) + 1);
    if (rho > registers_[j]) registers_[j] = rho;
  }
  // The canonical k-mers of a sequence or of a genome's contigs.
  void add(const Codes& s, int k) {
    for_each_kmer_hash(s, k, seed_, [&](std::uint64_t h) { insert(h); });
  }
  void add(const std::vector<Codes>& contigs, int k) {
    for (const auto& c : contigs) add(c, k);
  }

  // The union: register-wise maximum.
  void merge(const HyperLogLog& o) {
    if (o.p_ != p_ || o.seed_ != seed_) throw std::invalid_argument("genomics::sketch::HyperLogLog::merge: sketches of the same precision and seed");
    for (std::size_t j = 0; j < registers_.size(); j++) registers_[j] = std::max(registers_[j], o.registers_[j]);
  }

  // Ertl's improved estimate.
  double estimate() const {
    const int q = 64 - p_;
    std::vector<double> C(static_cast<std::size_t>(q) + 2, 0);
    for (auto r : registers_) C[r]++;
    const double m = static_cast<double>(registers_.size());
    double z = m * tau(1 - C[static_cast<std::size_t>(q) + 1] / m);
    for (int k = q; k >= 1; k--) z = 0.5 * (z + C[static_cast<std::size_t>(k)]);
    z += m * sigma(C[0] / m);
    if (!std::isfinite(z)) return 0;  // every register zero
    return m * m / (2 * std::log(2.0)) / z;
  }

  // Flajolet et al.'s estimate, with linear counting below 5m/2.
  double estimate_classic() const {
    const double m = static_cast<double>(registers_.size());
    double sum = 0;
    std::size_t zeros = 0;
    for (auto r : registers_) {
      sum += std::ldexp(1.0, -static_cast<int>(r));
      zeros += r == 0;
    }
    const double alpha = registers_.size() == 16 ? 0.673 : registers_.size() == 32 ? 0.697 : registers_.size() == 64 ? 0.709 : 0.7213 / (1 + 1.079 / m);
    const double e = alpha * m * m / sum;
    if (e <= 2.5 * m && zeros > 0) return m * std::log(m / static_cast<double>(zeros));
    return e;
  }

 private:
  // σ(x) = x + Σ_{k≥1} x^{2^k} 2^{k−1} (Ertl 2017, eq. 33); ∞ at x = 1.
  static double sigma(double x) {
    if (x == 1) return INFINITY;
    double y = 1, z = x, previous;
    do {
      x *= x;
      previous = z;
      z += x * y;
      y += y;
    } while (z != previous);
    return z;
  }
  // τ(x) = (1/3)(1 − x − Σ_{k≥1} (1 − x^{2^{−k}})² 2^{−k}) (eq. 37); 0 at x = 0 and x = 1.
  static double tau(double x) {
    if (x == 0 || x == 1) return 0;
    double y = 1, z = 1 - x, previous;
    do {
      x = std::sqrt(x);
      previous = z;
      y *= 0.5;
      z -= (1 - x) * (1 - x) * y;
    } while (z != previous);
    return z / 3;
  }

  int p_;
  std::uint64_t seed_;
  std::vector<std::uint8_t> registers_;
};

// The sketch of A ∪ B.
inline HyperLogLog union_of(const HyperLogLog& a, const HyperLogLog& b) {
  HyperLogLog u = a;
  u.merge(b);
  return u;
}

// J(A, B) by inclusion and exclusion of three estimates, clamped to [0, 1].
inline double jaccard(const HyperLogLog& a, const HyperLogLog& b) {
  const double na = a.estimate(), nb = b.estimate(), nu = union_of(a, b).estimate();
  if (!(nu > 0)) return 0;
  return std::clamp((na + nb - nu) / nu, 0.0, 1.0);
}

}  // namespace genomics::sketch
