// Checks of bubble logic (bubble_logic.hpp).
//
// - The AND–OR gate's truth table: no bubbles, none out; one bubble, at the OR outlet; two together, one at each.
// - Bubbles are conserved: injected = exited + in flight.

#include <physical_computing/fluids/bubble_logic.hpp>

#include <calx/verification/check.hpp>

using namespace physical_computing;

int main() {
  calx::verification::Suite s;
  using fluids::BubbleCircuit;
  const auto g00 = BubbleCircuit::evaluate(false, false), g10 = BubbleCircuit::evaluate(true, false), g01 = BubbleCircuit::evaluate(false, true),
             g11 = BubbleCircuit::evaluate(true, true);
  s.require("gate: no bubbles in, none out", g00.and_out == 0 && g00.or_out == 0);
  s.require("gate: A alone leaves at OR", g10.or_out == 1 && g10.and_out == 0);
  s.require("gate: B alone leaves at OR", g01.or_out == 1 && g01.and_out == 0);
  s.require("gate: A and B together, one at OR and one at AND", g11.or_out == 1 && g11.and_out == 1);
  auto g = BubbleCircuit::and_or_gate();
  g.inject(0, 0), g.inject(1, 1), g.inject(0, 0);
  g.run(0.01, 1e-3);
  s.require("conservation mid-run: injected = exited + in flight", g.injected() == g.exited() + g.in_flight());
  g.run(10, 1e-3);
  s.require("conservation: all have left", g.exited() == 3 && g.in_flight() == 0);
  return s.done();
}
