// A schedule: when inputs are applied and states read, kept apart from the program (ADR 0076 §2.3;
// schedule/README.md), as Halide keeps a schedule apart from its algorithm.
//
// The first schedule is a sample-and-hold: each input vector is held for `hold` seconds of the substrate's
// own time and the substrate is read at the end of each hold; the first `washout` holds are run but not
// counted, so that the initial state is forgotten (the echo state property). Others (annealing, virtual
// nodes, repetitions) are added as fields with defaults, so that a schedule written before still means the
// same.
#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace physical_computing::schedule {

struct Schedule {
  double hold = 1;            // s of the substrate's time per input
  std::size_t washout = 0;    // holds run before readings count
  std::size_t repetitions = 1;  // runs, each from its own seed
  std::uint64_t seed = 1;     // the first run's seed; run r uses seed + r

  void validate() const {
    if (!(hold > 0) || repetitions == 0) throw std::invalid_argument("physical_computing::schedule::Schedule: a positive hold and at least one repetition");
  }
  // The substrate's time a schedule of n inputs takes, every repetition included.
  double duration(std::size_t inputs) const { return hold * static_cast<double>(inputs) * static_cast<double>(repetitions); }
};

}  // namespace physical_computing::schedule
