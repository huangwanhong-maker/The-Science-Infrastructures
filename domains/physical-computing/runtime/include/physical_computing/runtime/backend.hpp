// The backend interface: what every simulator, device and surrogate answers (ADR 0076 §2.3; runtime/README.md).
//
// A backend is a substrate the runtime can drive. It describes itself (the target description the compiler
// reads), is reset under a seed, has inputs applied for a duration of its own time, and is read. Every reading
// carries its provenance, so that a reading from a device is never taken for one from a simulator. A family's
// simulator is the reference its devices are held to; a device is an instrument behind observation/, driven by
// a work, never by the domain, whose core stays pure.
//
// Extension. A new substrate implements Backend (or adapts its own interface to it); the runtime, the schedules
// and the checks of the square then apply unchanged. Fields are added to the descriptions with defaults, so that
// a backend written before still compiles and means the same.
#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <physical_computing/computer/materials.hpp>
#include <physical_computing/tasks/task.hpp>

namespace physical_computing::runtime {

enum class BackendKind { simulator, device, surrogate };

struct Channel {
  std::string name, unit;
};

struct Provenance {
  BackendKind kind = BackendKind::simulator;
  std::string family, substrate;  // "fluids", "convection cell"
  std::string model, version;     // the model and its source, the version of the code
  std::uint64_t seed = 0;
  std::vector<std::pair<std::string, double>> conditions;  // what the run was under: drive, temperature
  bool observed() const { return kind == BackendKind::device; }
};

struct TargetDescription {
  std::string family, substrate;
  std::vector<tasks::Kind> accepts;      // the kinds of task the compiler may lower onto it
  std::vector<Channel> inputs, outputs;
  double time_scale = 1;                 // s: the substrate's characteristic time
  double power = 0;                      // W while running, declared (0: not known)
  computer::Marks materials;
  std::string model;                     // the model and its sources
  std::vector<std::string> limits;       // what the model leaves out, in words
  bool accepts_kind(tasks::Kind k) const {
    for (auto a : accepts)
      if (a == k) return true;
    return false;
  }
};

class Backend {
 public:
  virtual ~Backend() = default;
  virtual TargetDescription describe() const = 0;
  virtual void reset(std::uint64_t seed) = 0;
  // Holds the inputs (one per input channel) for `duration` seconds of the substrate's time.
  virtual void apply(std::span<const double> inputs, double duration) = 0;
  // One value per output channel, now.
  virtual std::vector<double> read() const = 0;
  virtual double time() const = 0;
  virtual Provenance provenance() const = 0;
};

}  // namespace physical_computing::runtime
