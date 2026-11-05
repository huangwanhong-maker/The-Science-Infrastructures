// Checks of seesaw gates (seesaw.hpp; Qian and Winfree 2011).
//
// - Threshold 0.6 is OR: its output is above 0.6 when either input is on and below 0.05 when both are off.
// - Threshold 1.2 is AND: above 0.6 when both are on, below 0.05 otherwise.

#include <physical_computing/molecular/seesaw.hpp>

#include <calx/verification/check.hpp>

using namespace physical_computing;

int main() {
  calx::verification::Suite s;
  const auto OR = molecular::two_input_gate(0.6), AND = molecular::two_input_gate(1.2);
  auto y = [](const molecular::SeesawGate& g, bool a, bool b) { return molecular::run_gate(g, a, b); };
  s.require("OR: off when both are off", y(OR, false, false) < 0.05);
  s.require("OR: on with x₁", y(OR, true, false) > 0.6);
  s.require("OR: on with x₂", y(OR, false, true) > 0.6);
  s.require("OR: on with both", y(OR, true, true) > 0.6);
  s.require("AND: off when both are off", y(AND, false, false) < 0.05);
  s.require("AND: off with x₁ alone", y(AND, true, false) < 0.05);
  s.require("AND: off with x₂ alone", y(AND, false, true) < 0.05);
  s.require("AND: on with both", y(AND, true, true) > 0.6);
  return s.done();
}
