// Distances between aligned nucleotide sequences: the expected number of substitutions per
// site since their common ancestor, estimated from the differences between them under a model
// of substitution, with Γ-distributed rates across sites where a closed form exists.
//
// Gaps and N are removed pair by pair (pairwise deletion): two sequences are compared at the
// sites where both have a base. Let a pair have L such sites, P1 the fraction of A↔G
// differences, P2 of C↔T, P = P1 + P2 the transitions, Q the transversions, p = P + Q, and π
// the base frequencies (of the pair, as Tamura and Nei 1993 take them, or of the alignment).
// With α the shape of the Γ distribution of rates, each −ln x below becomes α(x^{−1/α} − 1):
//
//   p        the proportion of sites that differ;
//   JC69     d = −(3/4) ln(1 − 4p/3)                            (Jukes and Cantor 1969, in
//            Mammalian Protein Metabolism III, 21–132); Γ: d = (3/4) α [(1 − 4p/3)^{−1/α} − 1]
//            (Nei and Gojobori 1986, Mol. Biol. Evol. 3, 418–426, eq. 4 inverted);
//   K80      d = −(1/2) ln(1 − 2P − Q) − (1/4) ln(1 − 2Q)       (Kimura 1980, J. Mol. Evol. 16,
//            111–120, eq. 10); Γ: d = (α/2)[(1 − 2P − Q)^{−1/α} + (1/2)(1 − 2Q)^{−1/α} − 3/2]
//            (Jin and Nei 1990, Mol. Biol. Evol. 7, 82–102, eq. A4);
//   F81      d = −B ln(1 − p/B), B = 1 − Σ π_i²                (Felsenstein 1981, J. Mol. Evol.
//            17, 368–376; Tajima and Nei 1984, Mol. Biol. Evol. 1, 269–285);
//            Γ: d = B α [(1 − p/B)^{−1/α} − 1];
//   TN93     d = −(2π_Aπ_G/π_R) ln(1 − π_R P1/(2π_Aπ_G) − Q/(2π_R))
//              − (2π_Cπ_T/π_Y) ln(1 − π_Y P2/(2π_Cπ_T) − Q/(2π_Y))
//              − 2(π_Rπ_Y − π_Aπ_Gπ_Y/π_R − π_Cπ_Tπ_R/π_Y) ln(1 − Q/(2π_Rπ_Y))
//            (Tamura and Nei 1993, Mol. Biol. Evol. 10, 512–526, eq. 7; Γ: eq. 15);
//   LogDet   d = −(1/4)[ln det F − (1/2) ln(det Π_x det Π_y)], F the 4 × 4 matrix of joint
//            frequencies and Π the diagonal matrices of each sequence's base frequencies:
//            the paralinear distance (Lake 1994, PNAS 91, 1455–1459), the LogDet of Lockhart,
//            Steel, Hendy and Penny (1994, Mol. Biol. Evol. 11, 605–612) normalised so that it
//            equals the JC69 distance under JC69. It is additive on the tree under any (also
//            non-stationary) Markov model with rates equal across sites, so that tree building
//            on it is consistent, but under a stationary model it is −tr(Q) t / 4 (det P(t) =
//            e^{tr(Q) t}), not t, unless the frequencies are equal; it has no Γ form.
//
// Saturation: when an argument of a logarithm is not positive (the differences more than the
// model allows), or a pair shares no site, the distance is not defined; the pair is reported
// saturated and given the cap of the options (10 substitutions per site by default, far beyond
// any estimate that means something), never a NaN or an infinity, so that tree building
// proceeds and the caller sees which pairs to distrust. A TN93 term whose two bases are absent
// is zero, its limit.
//
// Complexity: O(n² · patterns) for n sequences, in parallel over rows on a pool when given.
// The estimates are the moment estimators of their papers; their variances are not computed.

#pragma once

#include <evolution/phylogeny/alignment.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/scheduler/parallel_for.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace evolution::phylogeny {

enum class DistanceModel { p, jc69, k80, f81, tn93, logdet };

struct DistanceOptions {
  DistanceModel model = DistanceModel::jc69;
  double alpha = INFINITY;          // the Γ shape; infinite for rates equal across sites
  bool pair_frequencies = true;     // π of the pair (Tamura and Nei 1993), else of the alignment
  double cap = 10;                  // the distance given a saturated pair
};

struct Distance {
  double d = 0;
  bool saturated = false;
};

// A symmetric matrix of distances with the names of its rows.
struct DistanceMatrix {
  std::vector<std::string> names;
  std::size_t n = 0;
  std::vector<double> d;                 // n × n
  std::vector<std::uint8_t> saturated;   // n × n
  std::vector<double> sites;             // n × n, the sites compared

  DistanceMatrix() = default;
  DistanceMatrix(std::vector<std::string> labels, std::vector<double> values) : names(std::move(labels)), n(names.size()), d(std::move(values)) {
    if (d.size() != n * n) throw std::invalid_argument("evolution::phylogeny::DistanceMatrix: n × n values");
    saturated.assign(n * n, 0);
    sites.assign(n * n, 0);
  }
  double operator()(std::size_t i, std::size_t j) const { return d[i * n + j]; }
  std::size_t saturated_pairs() const {
    std::size_t k = 0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = i + 1; j < n; j++) k += saturated[i * n + j];
    return k;
  }
};

namespace detail {

// −ln x, or α(x^{−1/α} − 1) with Γ rates; false when x ≤ 0.
inline bool minus_log(double x, double alpha, double& out) {
  if (!(x > 0)) return false;
  out = std::isinf(alpha) ? -std::log(x) : alpha * (std::pow(x, -1 / alpha) - 1);
  return true;
}


}  // namespace detail

inline Distance p_distance(double p) { return {p, false}; }

inline Distance jc69_distance(double p, double alpha = INFINITY) {
  Distance r;
  r.saturated = !detail::minus_log(1 - 4 * p / 3, alpha, r.d);
  r.d *= 0.75;
  return r;
}

// P: transitions per site; Q: transversions per site.
inline Distance k80_distance(double P, double Q, double alpha = INFINITY) {
  double a = 0, b = 0;
  Distance r;
  if (!detail::minus_log(1 - 2 * P - Q, alpha, a) || !detail::minus_log(1 - 2 * Q, alpha, b)) {
    r.saturated = true;
    return r;
  }
  r.d = a / 2 + b / 4;
  return r;
}

inline Distance f81_distance(double p, const std::array<double, 4>& pi, double alpha = INFINITY) {
  const double B = 1 - (pi[0] * pi[0] + pi[1] * pi[1] + pi[2] * pi[2] + pi[3] * pi[3]);
  Distance r;
  if (!(B > 0)) {
    r.saturated = p > 0;
    return r;
  }
  r.saturated = !detail::minus_log(1 - p / B, alpha, r.d);
  r.d *= B;
  return r;
}

// P1: A↔G differences per site; P2: C↔T; Q: transversions; π in the order A, C, G, T.
inline Distance tn93_distance(double P1, double P2, double Q, const std::array<double, 4>& pi, double alpha = INFINITY) {
  const double pA = pi[0], pC = pi[1], pG = pi[2], pT = pi[3];
  const double pR = pA + pG, pY = pC + pT;
  Distance r;
  double x = 0, d = 0;
  if (pA * pG > 0) {
    if (!detail::minus_log(1 - pR * P1 / (2 * pA * pG) - Q / (2 * pR), alpha, x)) r.saturated = true;
    d += 2 * pA * pG / pR * x;
  } else if (P1 > 0) {
    r.saturated = true;
  }
  if (pC * pT > 0) {
    if (!detail::minus_log(1 - pY * P2 / (2 * pC * pT) - Q / (2 * pY), alpha, x)) r.saturated = true;
    d += 2 * pC * pT / pY * x;
  } else if (P2 > 0) {
    r.saturated = true;
  }
  if (pR > 0 && pY > 0) {
    if (!detail::minus_log(1 - Q / (2 * pR * pY), alpha, x)) r.saturated = true;
    d += 2 * (pR * pY - pA * pG * pY / pR - pC * pT * pR / pY) * x;
  } else if (Q > 0) {
    r.saturated = true;
  }
  r.d = r.saturated ? 0 : d;
  return r;
}

// F: joint counts or frequencies of the pair's bases, F[4x + y] for x in the first, y in the second.
inline Distance logdet_distance(const std::array<double, 16>& F) {
  double total = 0;
  for (double v : F) total += v;
  Distance r;
  if (!(total > 0)) {
    r.saturated = true;
    return r;
  }
  std::array<double, 16> J{};
  std::array<double, 4> fx{}, fy{};
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 4; y++) {
      J[x * 4 + y] = F[x * 4 + y] / total;
      fx[x] += J[x * 4 + y];
      fy[y] += J[x * 4 + y];
    }
  calx::linalg::Matrix Jm(4);
  std::copy(J.begin(), J.end(), Jm.a.begin());
  const double det = calx::linalg::determinant(Jm);
  const double dx = fx[0] * fx[1] * fx[2] * fx[3], dy = fy[0] * fy[1] * fy[2] * fy[3];
  if (!(det > 0) || !(dx > 0) || !(dy > 0)) {
    r.saturated = true;
    return r;
  }
  r.d = -0.25 * (std::log(det) - 0.5 * (std::log(dx) + std::log(dy)));
  return r;
}

// The joint counts of the bases of rows i and j at the sites where both are known.
inline std::array<double, 16> pair_counts(const Patterns& p, std::size_t i, std::size_t j) {
  std::array<double, 16> F{};
  for (std::size_t q = 0; q < p.count; q++) {
    const Code x = p.states[q * p.taxa + i], y = p.states[q * p.taxa + j];
    if (x < 4 && y < 4) F[x * 4 + y] += p.weight[q];
  }
  return F;
}

// The distance of a pair from its joint counts.
inline Distance distance_from_counts(const std::array<double, 16>& F, const DistanceOptions& o, const std::array<double, 4>& alignment_pi) {
  double L = 0;
  for (double v : F) L += v;
  if (!(L > 0)) return {o.cap, true};
  using genomics::sequence::A;
  using genomics::sequence::C;
  using genomics::sequence::G;
  using genomics::sequence::T;
  const double P1 = (F[A * 4 + G] + F[G * 4 + A]) / L, P2 = (F[C * 4 + T] + F[T * 4 + C]) / L;
  double same = 0;
  for (int x = 0; x < 4; x++) same += F[x * 4 + x];
  const double p = 1 - same / L, Q = std::max(0.0, p - P1 - P2);
  std::array<double, 4> pi = alignment_pi;
  if (o.pair_frequencies) {
    for (int x = 0; x < 4; x++) {
      double s = 0;
      for (int y = 0; y < 4; y++) s += F[x * 4 + y] + F[y * 4 + x];
      pi[x] = s / (2 * L);
    }
  }
  Distance r;
  switch (o.model) {
    case DistanceModel::p: r = p_distance(p); break;
    case DistanceModel::jc69: r = jc69_distance(p, o.alpha); break;
    case DistanceModel::k80: r = k80_distance(P1 + P2, Q, o.alpha); break;
    case DistanceModel::f81: r = f81_distance(p, pi, o.alpha); break;
    case DistanceModel::tn93: r = tn93_distance(P1, P2, Q, pi, o.alpha); break;
    case DistanceModel::logdet: r = logdet_distance(F); break;
  }
  if (r.saturated || !std::isfinite(r.d)) return {o.cap, true};
  return r;
}

// The matrix of distances between every two rows of the patterns.
inline DistanceMatrix distances(const Patterns& p, const DistanceOptions& o = {}, calx::scheduler::ThreadPool* pool = nullptr) {
  DistanceMatrix m;
  m.names = p.names;
  m.n = p.taxa;
  const std::size_t n = m.n;
  m.d.assign(n * n, 0);
  m.saturated.assign(n * n, 0);
  m.sites.assign(n * n, 0);
  const auto pi = base_frequencies(p);
  const auto row = [&](std::size_t i) {
    for (std::size_t j = i + 1; j < n; j++) {
      const auto F = pair_counts(p, i, j);
      double L = 0;
      for (double v : F) L += v;
      const Distance r = distance_from_counts(F, o, pi);
      m.d[i * n + j] = m.d[j * n + i] = r.d;
      m.saturated[i * n + j] = m.saturated[j * n + i] = r.saturated;
      m.sites[i * n + j] = m.sites[j * n + i] = L;
    }
  };
  if (pool && pool->size() > 1) {
    calx::scheduler::parallel_for(*pool, n, 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t i = b; i < e; i++) row(i);
    });
  } else {
    for (std::size_t i = 0; i < n; i++) row(i);
  }
  return m;
}

inline DistanceMatrix distances(const Alignment& a, const DistanceOptions& o = {}, calx::scheduler::ThreadPool* pool = nullptr) {
  return distances(compress(a), o, pool);
}

}  // namespace evolution::phylogeny
