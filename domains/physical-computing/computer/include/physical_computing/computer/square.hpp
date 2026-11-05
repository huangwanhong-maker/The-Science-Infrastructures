// The commuting square: when a substrate computes (theory note §1.1; ADR 0076 §2.2).
//
// A physical computer (S, Φ_u, M, C, E, R, ε) computes C on a set K when
//
//   d(R(Φ_u(E(m))), C(m)) ≤ ε   for every m ∈ K,
//
// and it is non-trivial when E and R are of a class fixed in advance and cheaper than C (Horsman, Stepney,
// Wagner and Kendon 2014). commutes() evaluates the square on K and returns the largest error, where it was
// reached, and whether the tolerance holds; the class of E and R is declared beside it, not inferred, since
// no program can tell whether a readout secretly does the computation. For chaotic flows the cases of K are
// statistics of runs (invariants, theory note §1.4), not trajectories: the square is the same, with Φ a run
// and R an estimator.
#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>

namespace physical_computing::computer {

// The class of an encoding or a readout, declared.
enum class Map {
  identity,        // the state is the result
  rescaling,       // a fixed change of units or of variables (analogue computing by similarity)
  linear,          // a fixed linear map
  trained_linear,  // a linear map fitted to examples (reservoirs)
  threshold,       // a comparison with a level (bits, the sign of an amplitude, a phase binarised)
  estimator,       // a statistic of a run (invariants of chaotic flows)
};

struct Verdict {
  double error = 0;            // the largest d over K
  double tolerance = 0;        // ε
  bool holds = false;          // error ≤ ε (and finite)
  std::size_t cases = 0;       // |K|
  std::size_t worst = 0;       // where the largest error was reached
  Map encoding = Map::identity, readout = Map::identity;
};

// K: a range of abstract states; E(m) → substrate state; Φ(s) → evolved state; R(s) → abstract state;
// C(m) → the model's answer; d(a, b) → error.
template <class Cases, class Encode, class Evolve, class Read, class Model, class Metric>
Verdict commutes(const Cases& K, Encode&& E, Evolve&& Phi, Read&& R, Model&& C, Metric&& d, double tolerance,
                 Map encoding = Map::identity, Map readout = Map::identity) {
  Verdict v;
  v.tolerance = tolerance, v.encoding = encoding, v.readout = readout;
  std::size_t i = 0;
  for (const auto& m : K) {
    const double e = d(R(Phi(E(m))), C(m));
    if (!std::isfinite(e)) {
      v.error = std::numeric_limits<double>::infinity(), v.worst = i;
    } else if (e > v.error) {
      v.error = e, v.worst = i;
    }
    i++;
  }
  v.cases = i;
  v.holds = std::isfinite(v.error) && v.error <= tolerance && i > 0;
  return v;
}

}  // namespace physical_computing::computer
