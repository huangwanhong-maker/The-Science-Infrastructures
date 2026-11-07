// Checks of exposure to a hazard on a grid (exposure/exposure.hpp).
//
// - A grid of known depths, its rows from the south and its columns from the west: a place reads its cell, and outside it reads nothing.
// - Assets on a wet half and a dry half: the count and value reached, those beyond the grid apart, each band and zone adding up to the whole.
// - A threshold: at it an asset is reached, below it not.
// - A declared damage curve: linear between points, constant beyond them, the damage the asset's value times its ratio.
// - The curve read between its points, and refusals of a grid of the wrong size.

#include <urban/exposure/exposure.hpp>

#include <calx/verification/check.hpp>

#include <stdexcept>
#include <vector>

namespace ue = urban::exposure;

int main() {
  calx::verification::Suite s;
  ue::Grid g;
  g.south = 0, g.north = 2, g.west = 10, g.east = 14, g.width = 4, g.height = 2, g.unit = 0.5;
  // Rows from the south: the first row 0, 2, 4, 6 (0, 1, 2, 3 m), the second row all dry.
  g.values = {0, 2, 4, 6, 0, 0, 0, 0};
  s.check("a place reads its cell, rows from the south", g.at(0.5, 12.5) - 2, 0);
  s.check("the next row is another", g.at(1.5, 12.5), 0);
  s.require("outside the grid it reads nothing", g.at(3, 12) < 0 && g.at(1, 9) < 0);
  const std::vector<ue::Asset> A = {{0.5, 10.5, 100, 0}, {0.5, 11.5, 200, 0}, {0.5, 12.5, 400, 1}, {0.5, 13.5, 800, 1}, {1.5, 12.0, 1600, 1}, {5.0, 12.0, 3200, 0}};
  const auto E = ue::exposed(g, A, 0.5, {1.5, 2.5}, {0, 2}, {0, 1}, 2);
  s.require("those reached are counted", E.reached.count == 3, std::to_string(E.reached.count));
  s.check("and their value", E.reached.value - 1400, 0);
  s.require("those with a reading, and those beyond the grid apart", E.all.count == 5 && E.outside == 1);
  s.check("all with a reading", E.all.value - 3100, 0);
  s.require("the bands add up to those reached", E.by_band[0].count + E.by_band[1].count + E.by_band[2].count == 3 && E.by_band[0].count == 1 && E.by_band[1].count == 1 && E.by_band[2].count == 1);
  s.require("and the zones", E.zones[0].count == 1 && E.zones[1].count == 2);
  s.check("the damage is the value times the curve: 1 m at half, 2 m and 3 m whole", E.reached.damage - (200 * 0.5 + 400 * 1 + 800 * 1), 1e-12);
  const auto T = ue::exposed(g, A, 1.0);
  s.require("at the threshold an asset is reached", T.reached.count == 3);
  s.require("above it, it is not", ue::exposed(g, A, 1.5).reached.count == 2);
  s.check("the curve between its points", ue::interpolate({0, 1, 3}, {0, 0.5, 1}, 2) - 0.75, 1e-12);
  s.check("and constant beyond them", ue::interpolate({0, 1, 3}, {0.1, 0.5, 1}, -1) - 0.1, 0);
  bool refused = false;
  g.values.pop_back();
  try { ue::exposed(g, A, 0); } catch (const std::invalid_argument&) { refused = true; }
  s.require("a grid of the wrong size is refused", refused);
  return s.done();
}
