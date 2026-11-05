// A task: what is to be computed, on no particular substrate (ADR 0076 §2.3; tasks/README.md).
//
// A task names its kind, the metric by which its result is judged, and the tolerance it must be met to.
// Its set K, its model C and its check are the caller's, passed to computer::commutes (square.hpp); until
// the intermediate representation exists (ir/, optional and future), a task is this value and those
// callables. New kinds are added at the end of the enumeration, so that a kind's number never changes.
#pragma once

#include <array>
#include <string>
#include <string_view>

namespace physical_computing::tasks {

enum class Kind {
  filter,      // a causal map from input series to output series, with fading memory
  flow,        // a dynamical system to realise, judged by its invariants or over a horizon
  energy,      // an Ising, QUBO or Potts problem, to minimise or to sample
  function,    // a function on counts or bits, or a linear map
  phenomenon,  // a behaviour to produce, judged against its normal form
  learning,    // a filter the substrate learns by its own rule
};

inline constexpr std::array<Kind, 6> all_kinds{Kind::filter, Kind::flow, Kind::energy, Kind::function, Kind::phenomenon, Kind::learning};

inline std::string_view name(Kind k) {
  switch (k) {
    case Kind::filter: return "filter";
    case Kind::flow: return "flow";
    case Kind::energy: return "energy";
    case Kind::function: return "function";
    case Kind::phenomenon: return "phenomenon";
    case Kind::learning: return "learning";
  }
  return "unknown";
}

// The tolerance of a task: the largest error d(R(Φ(E(m))), C(m)) allowed, in the metric's unit.
struct Tolerance {
  double value = 0;
  std::string metric;  // e.g. "absolute error", "relative error of the Lyapunov spectrum", "NMSE"
  std::string unit;
};

struct Task {
  std::string name;
  Kind kind = Kind::function;
  Tolerance tolerance;
  std::string check;  // how the result is verified, in words; the check itself is a callable
};

}  // namespace physical_computing::tasks
