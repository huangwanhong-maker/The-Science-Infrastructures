// The best secondary structure for additive weights, and its two uses here: Nussinov's maximal
// base pairing and the maximum expected accuracy structure (theory/rna-structure.md §2 and §4).
//
// One recursion. Given a weight u_i for leaving position i unpaired and a weight w(i, j) for the
// pair (i, j) (or none, when the pair is not allowed), the structure maximising
// Σ_{i unpaired} u_i + Σ_{(i, j) paired} w(i, j) over the secondary structures with minimum loop m
// is found by
//
//   M(i, i) = 0,
//   M(i, j) = max( u_i + M(i + 1, j),
//                  max_{k : i + m + 1 ≤ k < j, w(i, k) allowed} w(i, k) + M(i + 1, k) + M(k + 1, j) )
//
// over the half-open intervals [i, j): position i is unpaired, or paired with some k, and then
// the interior (i + 1, k) and the rest (k + 1, j) are independent because pairs do not cross.
// Each structure has exactly one such decomposition, so the maximum is exact; time O(n³) and
// space O(n²). The traceback prefers the unpaired case, then the least k, on ties.
//
// Nussinov's maximal base pairing (Nussinov, Pieczenik, Griggs and Kleitman 1978, Algorithms
// for loop matchings, SIAM J. Appl. Math. 35, 68–82; Nussinov and Jacobson 1980, Fast algorithm
// for predicting the secondary structure of single-stranded RNA, PNAS 77, 6309–6313): u = 0,
// w = 1 on the canonical pairs (Watson–Crick, and G·U when asked), over exact integers. It
// counts pairs and nothing else: no stacking, no loop costs, so its structures are rarely the
// real ones, and it is here as the simplest exact case of the recursion.
//
// The maximum expected accuracy (MEA) structure from base-pair probabilities p_ij (a grammar's,
// scfg.hpp, or any model's): with q_i = 1 − Σ_j p_ij the probability that i is unpaired, the
// expected number of positions whose state (unpaired, or paired with the right partner) a
// structure predicts correctly is Σ_{i unpaired} q_i + Σ_{(i, j)} 2 p_ij; weighting the pairs by
// γ > 0 gives the objective Σ_{i unpaired} q_i + 2γ Σ_{(i, j)} p_ij of Do, Woods and
// Batzoglou (2006, CONTRAfold: RNA secondary structure prediction without physics-based
// models, Bioinformatics 22, e90–e98), larger γ trading specificity for sensitivity. It is the
// recursion with u_i = q_i and w(i, j) = 2γ p_ij on the pairs with p_ij > 0. As γ grows, the
// expected number of correct paired positions of the optimum, Σ 2 p_ij, never decreases
// (theory, Proposition 4.4).

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

#include "structure.hpp"

namespace genomics::rna {

template <class T>
struct BestStructure {
  Structure structure;
  T score{};
};

// The maximising structure of n positions. unpaired(i) is u_i; pair(i, j) is w(i, j), or nothing
// when the pair is not allowed.
template <class T>
BestStructure<T> best_structure(std::size_t n, const std::function<T(std::size_t)>& unpaired,
                                const std::function<std::optional<T>(std::size_t, std::size_t)>& pair, std::size_t min_loop = 0) {
  // M[i][j] for 0 ≤ i ≤ j ≤ n, row-major in (n + 1)²; pair weights cached the same way.
  const std::size_t m = n + 1;
  std::vector<T> M(m * m, T{});
  std::vector<std::optional<T>> w(n * n);
  std::vector<T> u(n);
  for (std::size_t i = 0; i < n; i++) {
    u[i] = unpaired(i);
    for (std::size_t k = i + min_loop + 1; k < n; k++) w[i * n + k] = pair(i, k);
  }
  for (std::size_t len = 1; len <= n; len++)
    for (std::size_t i = 0; i + len <= n; i++) {
      const std::size_t j = i + len;
      T best = u[i] + M[(i + 1) * m + j];
      for (std::size_t k = i + min_loop + 1; k < j; k++) {
        const auto& wk = w[i * n + k];
        if (!wk) continue;
        const T v = *wk + M[(i + 1) * m + k] + M[(k + 1) * m + j];
        if (best < v) best = v;
      }
      M[i * m + j] = best;
    }
  BestStructure<T> out;
  out.structure = Structure(n);
  out.score = n ? M[0 * m + n] : T{};
  // Traceback over a stack of intervals.
  std::vector<std::pair<std::size_t, std::size_t>> todo{{0, n}};
  while (!todo.empty()) {
    const auto [i, j] = todo.back();
    todo.pop_back();
    if (i >= j) continue;
    const T target = M[i * m + j];
    if (u[i] + M[(i + 1) * m + j] == target) { todo.push_back({i + 1, j}); continue; }
    bool found = false;
    for (std::size_t k = i + min_loop + 1; k < j && !found; k++) {
      const auto& wk = w[i * n + k];
      if (!wk) continue;
      if (*wk + M[(i + 1) * m + k] + M[(k + 1) * m + j] == target) {
        out.structure.add(i, k);
        todo.push_back({i + 1, k});
        todo.push_back({k + 1, j});
        found = true;
      }
    }
    if (!found) throw std::logic_error("genomics::rna::best_structure: the traceback lost its way");
  }
  return out;
}

// ---- Nussinov's maximal base pairing ------------------------------------------------------------------
struct NussinovOptions {
  std::size_t min_loop = 3;  // the fewest unpaired bases a hairpin encloses
  bool wobble = true;        // allow G·U pairs beside the Watson–Crick ones
};

// A structure with the most canonical pairs; score is their number.
inline BestStructure<std::int64_t> nussinov(const Codes& x, const NussinovOptions& o = {}) {
  return best_structure<std::int64_t>(
      x.size(), [](std::size_t) { return std::int64_t{0}; },
      [&](std::size_t i, std::size_t j) -> std::optional<std::int64_t> {
        if (canonical(x[i], x[j], o.wobble)) return std::int64_t{1};
        return std::nullopt;
      },
      o.min_loop);
}

// ---- maximum expected accuracy ----------------------------------------------------------------------------
// Base-pair probabilities of a sequence of n positions: p[i * n + j] = p[j * n + i] = P(i·j).
struct PairProbabilities {
  std::size_t n = 0;
  std::vector<double> p;
  double log_likelihood = 0;  // of the sequence, under the model that gave them (natural log)

  double at(std::size_t i, std::size_t j) const { return p[i * n + j]; }
  // q_i = 1 − Σ_j p_ij, clamped at 0 against rounding.
  double unpaired(std::size_t i) const {
    double s = 0;
    for (std::size_t j = 0; j < n; j++) s += p[i * n + j];
    return s < 1 ? 1 - s : 0.0;
  }
  // The expected number of pairs, Σ_{i<j} p_ij.
  double expected_pairs() const {
    double s = 0;
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = i + 1; j < n; j++) s += p[i * n + j];
    return s;
  }
};

// The MEA structure for weight γ on the pairs; its score is the expected accuracy
// Σ_{i unpaired} q_i + 2γ Σ_{(i, j)} p_ij. Pairs with p_ij = 0 or a loop shorter than min_loop are
// not considered.
inline BestStructure<double> maximum_expected_accuracy(const PairProbabilities& bp, double gamma = 1.0, std::size_t min_loop = 0) {
  if (!(gamma > 0)) throw std::invalid_argument("genomics::rna::maximum_expected_accuracy: γ must be positive");
  std::vector<double> q(bp.n);
  for (std::size_t i = 0; i < bp.n; i++) q[i] = bp.unpaired(i);
  return best_structure<double>(
      bp.n, [&](std::size_t i) { return q[i]; },
      [&](std::size_t i, std::size_t j) -> std::optional<double> {
        const double p = bp.at(i, j);
        if (!(p > 0)) return std::nullopt;
        return 2 * gamma * p;
      },
      min_loop);
}

// The expected accuracy of any structure, by the same objective.
inline double expected_accuracy(const PairProbabilities& bp, const Structure& s, double gamma = 1.0) {
  if (s.size() != bp.n) throw std::invalid_argument("genomics::rna::expected_accuracy: a structure of another length");
  double a = 0;
  for (std::size_t i = 0; i < bp.n; i++) {
    if (!s.paired(i)) a += bp.unpaired(i);
    else if (s.partner[i] > static_cast<std::int32_t>(i)) a += 2 * gamma * bp.at(i, static_cast<std::size_t>(s.partner[i]));
  }
  return a;
}

}  // namespace genomics::rna
