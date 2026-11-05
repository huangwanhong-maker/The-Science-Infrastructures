// Electrochemical oscillators: an array of nickel electrodes dissolving in acid, each an oscillator, coupled
// through a common resistance (Kiss, Zhai and Hudson 2002) (substrates/chemical/README.md, simulator 4).
//
// The model is the array's phase reduction: each electrode's oscillation is a phase with its own natural
// frequency, and the common resistance couples them all to all, so that the array is Kuramoto's network
// (physis::dynamics::Kuramoto). Kiss, Zhai and Hudson saw its transition to partial synchrony on 64 electrodes as
// the coupling grew. The electrodes' frequencies are drawn from a declared Lorentzian of half-width γ about ω₀,
// for which the transition is at K_c = 2γ and the order parameter r = √(1 − 2γ/K) above it.
//
// Limits (declared, provisional): phases only, not the electrodes' current and potential (a two-variable model
// of nickel's electrodissolution is to be chosen from Hudson's group's papers); a sinusoidal coupling; the map from
// the resistance to K is not modelled, K is the input. Checks: tests/check_electrochemical.cpp.
#pragma once

#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <physis/dynamics/flows.hpp>
#include <physis/dynamics/kuramoto.hpp>

#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::chemical {

class ElectrodeArray {
 public:
  ElectrodeArray(std::size_t electrodes, double centre_frequency, double half_width)
      : net_{physis::dynamics::lorentzian_frequencies(electrodes, centre_frequency, half_width), 0.0, {}}, gamma_(half_width), phase_(electrodes) {
    if (electrodes < 2 || !(half_width > 0)) throw std::invalid_argument("physical_computing::chemical::ElectrodeArray: two electrodes or more, a positive spread");
    for (std::size_t i = 0; i < electrodes; i++) phase_[i] = 2 * 3.141592653589793 * std::fmod(0.6180339887498949 * static_cast<double>(i), 1.0);
  }
  double critical_coupling() const { return 2 * gamma_; }
  double predicted_order(double K) const { return physis::dynamics::kuramoto_lorentzian_order(K, gamma_); }
  // Runs at coupling K for `duration` and returns the mean order parameter over the last `average` of it.
  double run(double K, double duration, double average, double step = 0.01) {
    net_.coupling = K;
    physis::dynamics::integrate(net_, 0.0, duration - average, phase_, step);
    double r = 0;
    const int m = 100;
    for (int k = 0; k < m; k++) {
      physis::dynamics::integrate(net_, 0.0, average / m, phase_, step);
      r += std::abs(physis::dynamics::order_parameter(phase_));
    }
    return r / m;
  }
  const physis::dynamics::State& phases() const { return phase_; }

  runtime::TargetDescription describe() const {
    runtime::TargetDescription d;
    d.family = "chemical", d.substrate = "electrochemical oscillators";
    d.accepts = {tasks::Kind::phenomenon, tasks::Kind::energy};
    d.inputs = {{"coupling K (through the common resistance)", "1/s"}};
    d.outputs = {{"electrodes' phases", "rad"}};
    d.materials.materials = {{"nickel", "the electrodes"}, {"sulfuric acid", "the electrolyte"}, {"platinum group metals", "the counter electrode"}};
    d.materials.hazards = {"sulfuric acid: corrosive", "nickel salts: sensitising, toxic"};
    d.materials.reuse = computer::Reuse::degrades_with_use;
    d.model = "Kuramoto's network as the electrodes' phase reduction (physis::dynamics; Kiss, Zhai and Hudson 2002)";
    d.limits = {"phases only", "sinusoidal coupling", "K given, not derived from the resistance"};
    return d;
  }

 private:
  physis::dynamics::Kuramoto net_;
  double gamma_;
  physis::dynamics::State phase_;
};

}  // namespace physical_computing::chemical
