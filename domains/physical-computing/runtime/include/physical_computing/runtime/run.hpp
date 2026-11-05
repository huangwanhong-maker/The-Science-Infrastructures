// Programs driven on backends under a schedule (ADR 0076 §2.3; runtime/README.md).
//
// drive() applies an input series to a backend under a sample-and-hold schedule and returns the readings after
// the washout, run by run, with each run's provenance, the substrate time taken and, where the backend declares
// its power, the energy of operation. It knows nothing of any substrate: it is the same for a tank of water, a
// chemical medium or a device behind observation/.
#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

#include <physical_computing/runtime/backend.hpp>
#include <physical_computing/schedule/schedule.hpp>

namespace physical_computing::runtime {

struct Run {
  std::vector<std::vector<double>> readings;  // per counted step, one value per output channel
  Provenance provenance;
  double duration = 0;  // s of substrate time
};

struct Drive {
  std::vector<Run> runs;
  double energy = 0;          // J of operation, power × time; 0 when the power is not declared
  bool energy_known = false;
};

inline Drive drive(Backend& backend, const std::vector<std::vector<double>>& inputs, const schedule::Schedule& s) {
  s.validate();
  const auto target = backend.describe();
  for (const auto& u : inputs)
    if (u.size() != target.inputs.size()) throw std::invalid_argument("physical_computing::runtime::drive: one input per input channel");
  Drive out;
  for (std::size_t r = 0; r < s.repetitions; r++) {
    backend.reset(s.seed + r);
    const double t0 = backend.time();
    Run run;
    for (std::size_t k = 0; k < inputs.size(); k++) {
      backend.apply(inputs[k], s.hold);
      if (k >= s.washout) run.readings.push_back(backend.read());
    }
    run.duration = backend.time() - t0;
    run.provenance = backend.provenance();
    out.runs.push_back(std::move(run));
  }
  if (target.power > 0) {
    out.energy_known = true;
    for (const auto& r : out.runs) out.energy += target.power * r.duration;
  }
  return out;
}

}  // namespace physical_computing::runtime
