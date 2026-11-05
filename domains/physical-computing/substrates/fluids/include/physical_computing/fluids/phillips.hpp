// Phillips's hydraulic machine of the economy: water as a computer of the circular flow of income
// (Phillips 1950) (substrates/fluids/README.md, simulator 9).
//
// The analogue, kept to its income–expenditure core. A tank of area A holds water at level h; it drains
// through a linear valve of resistance R, so its outflow Q = h/R stands for income Y. A share c of that flow
// (the propensity to consume) is pumped back in, and injections I (investment and government spending, an
// input) are added:
//
//   A dh/dt = c h/R + I − h/R.
//
// The map Y = h/R carries it exactly onto the income model dY/dt = (c Y + I − Y)/τ, τ = A R, a conjugacy of
// the hydraulic flow and the economic one (theory note §1.2, the analogue mode). Its steady income is the
// multiplier's, Y* = I/(1 − c), reached with the time constant A R/(1 − c). Each hold is solved exactly, since
// the system is linear.
//
// Limits (declared): the machine's monetary sector (its balances and interest rate), its lags and its
// cams are left out; a linear valve. Checks: tests/check_phillips.cpp.
#pragma once

#include <cmath>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

#include <physical_computing/fluids/fluid.hpp>
#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::fluids {

struct PhillipsOptions {
  double area = 0.05;          // A, m²
  double resistance = 2.0e3;   // R, s/m² (Q = h/R, m³/s per m of level)
  double propensity = 0.75;    // c, the share returned as consumption
};

class PhillipsMachine final : public runtime::Backend {
 public:
  explicit PhillipsMachine(const PhillipsOptions& o) : o_(o) {
    if (!(o.area > 0) || !(o.resistance > 0) || !(o.propensity >= 0 && o.propensity < 1)) throw std::invalid_argument("physical_computing::fluids::PhillipsMachine: positive area and resistance, 0 ≤ c < 1");
  }
  double income() const { return h_ / o_.resistance; }
  double level() const { return h_; }
  double steady_income(double injections) const { return injections / (1 - o_.propensity); }
  double time_constant() const { return o_.area * o_.resistance / (1 - o_.propensity); }

  runtime::TargetDescription describe() const override {
    runtime::TargetDescription d;
    d.family = "fluids", d.substrate = "Phillips's machine";
    d.accepts = {tasks::Kind::flow, tasks::Kind::filter};
    d.inputs = {{"injections I", "m³/s"}};
    d.outputs = {{"income Y = h/R", "m³/s"}, {"consumption c Y", "m³/s"}, {"level h", "m"}};
    d.time_scale = time_constant();
    d.materials.materials = {{"water", "money"}, {"acrylic (PMMA)", "the tanks"}, {"copper", "the pumps' windings"}};
    d.materials.end_of_life = computer::EndOfLife::recyclable;
    d.model = "a tank, a linear valve and a return of the share c, conjugate to the income–expenditure model (Phillips 1950)";
    d.limits = {"the monetary sector, lags and cams left out", "a linear valve"};
    return d;
  }
  void reset(std::uint64_t seed) override { seed_ = seed, h_ = 0, t_ = 0; }
  // Holds the injections I for `duration` s, exactly.
  void apply(std::span<const double> inputs, double duration) override {
    if (inputs.size() != 1) throw std::invalid_argument("physical_computing::fluids::PhillipsMachine::apply: one input, the injections");
    const double hs = o_.resistance * steady_income(inputs[0]);
    h_ = hs + (h_ - hs) * std::exp(-duration / time_constant());
    t_ += duration;
  }
  std::vector<double> read() const override { return {income(), o_.propensity * income(), h_}; }
  double time() const override { return t_; }
  runtime::Provenance provenance() const override {
    auto p = simulator_provenance("Phillips's machine", "a linear tank conjugate to the income model, solved exactly", seed_);
    p.conditions = {{"c", o_.propensity}};
    return p;
  }

 private:
  PhillipsOptions o_;
  double h_ = 0, t_ = 0;
  std::uint64_t seed_ = 0;
};

}  // namespace physical_computing::fluids
