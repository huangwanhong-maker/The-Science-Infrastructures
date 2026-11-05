// Checks of the runtime (backend.hpp, run.hpp) on a substrate whose answer is known in closed form: a leaky
// integrator, dx/dt = (u − x)/τ, exactly x(t + h) = u + (x(t) − u) e^{−h/τ}.
//
// - drive() under a sample-and-hold schedule returns the closed form at every counted step (1e−15).
// - The washout drops the first readings; repetitions are run from their own seeds, and the provenance says
//   which; the substrate time and the energy (declared power × time) are reported.
// - A schedule with a non-positive hold, or inputs of the wrong width, are refused.

#include <physical_computing/runtime/run.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

using namespace physical_computing;

class LeakyIntegrator final : public runtime::Backend {
 public:
  explicit LeakyIntegrator(double tau) : tau_(tau) {}
  runtime::TargetDescription describe() const override {
    runtime::TargetDescription d;
    d.family = "test", d.substrate = "leaky integrator";
    d.accepts = {tasks::Kind::filter};
    d.inputs = {{"u", "1"}}, d.outputs = {{"x", "1"}};
    d.time_scale = tau_, d.power = 2.0;
    d.model = "dx/dt = (u − x)/τ, exactly";
    return d;
  }
  void reset(std::uint64_t seed) override { x_ = 0, t_ = 0, seed_ = seed; }
  void apply(std::span<const double> u, double h) override {
    x_ = u[0] + (x_ - u[0]) * std::exp(-h / tau_);
    t_ += h;
  }
  std::vector<double> read() const override { return {x_}; }
  double time() const override { return t_; }
  runtime::Provenance provenance() const override {
    runtime::Provenance p;
    p.family = "test", p.substrate = "leaky integrator", p.model = "closed form", p.seed = seed_;
    return p;
  }

 private:
  double tau_, x_ = 0, t_ = 0;
  std::uint64_t seed_ = 0;
};

int main() {
  calx::verification::Suite s;
  LeakyIntegrator b(0.5);
  const std::vector<std::vector<double>> u{{1}, {1}, {0}, {2}, {-1}};
  schedule::Schedule sch;
  sch.hold = 0.2, sch.washout = 2, sch.repetitions = 3, sch.seed = 40;
  const auto d = runtime::drive(b, u, sch);
  double x = 0, worst = 0;
  std::vector<double> expected;
  for (const auto& v : u) x = v[0] + (x - v[0]) * std::exp(-0.2 / 0.5), expected.push_back(x);
  for (const auto& run : d.runs)
    for (std::size_t k = 0; k < run.readings.size(); k++) worst = std::max(worst, std::fabs(run.readings[k][0] - expected[k + 2]));
  s.check("drive: the readings are the closed form at every counted step", worst, 1e-15);
  s.require("drive: the washout drops the first two readings", d.runs.size() == 3 && d.runs[0].readings.size() == 3);
  s.require("drive: each repetition from its own seed", d.runs[0].provenance.seed == 40 && d.runs[2].provenance.seed == 42);
  s.require("drive: a simulator's readings are not observed", !d.runs[0].provenance.observed());
  s.check("drive: substrate time 5 × 0.2 s a run (s)", d.runs[1].duration - 1.0, 1e-15, "s");
  s.check("drive: energy 2 W × 3 s (J)", d.energy - 6.0, 1e-12, "J");
  bool refused = false;
  try {
    schedule::Schedule bad;
    bad.hold = 0;
    runtime::drive(b, u, bad);
  } catch (const std::invalid_argument&) {
    refused = true;
  }
  s.require("drive: a non-positive hold is refused", refused);
  refused = false;
  try {
    runtime::drive(b, {{1, 2}}, sch);
  } catch (const std::invalid_argument&) {
    refused = true;
  }
  s.require("drive: inputs of the wrong width are refused", refused);
  return s.done();
}
