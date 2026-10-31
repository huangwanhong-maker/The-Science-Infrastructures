// Checks of the massing (massing.hpp): the buildings of a district as one regular 3-complex.
// - Two boxes side by side, 10 m and 20 m high: two 3-cells, regular; their shared wall one 2-cell with two cofaces,
//   the taller's wall above it its own; ∂∂ = 0 in every degree; the union a ball (Betti 1 0 0 0).
// - A third box along both their backs, meeting the corner between them in the middle of its edge (a T-junction):
//   its edge cut there, so that it shares a wall with each; still regular, still a ball.
// - A bow tie and a clockwise footprint skipped, with their reasons; a footprint with a repeated corner too.

#include <cmath>
#include <string>
#include <vector>

#include <architecture/space/massing.hpp>
#include <calx/verification/check.hpp>

using calx::verification::Suite;
namespace as = architecture::space;
namespace mt = morphe::topology;

namespace {

std::string text(const std::vector<std::size_t>& v) {
  std::string s;
  for (std::size_t x : v) s += (s.empty() ? "" : " ") + std::to_string(x);
  return s;
}

// ∂_k ∂_{k+1} = 0 for every k, over the integers (morphe's validation of a chain complex).
bool boundary_squares_to_zero(const mt::RegularComplex& X) {
  try { mt::chain_complex(X).validate(); return true; } catch (const std::invalid_argument&) { return false; }
}

std::size_t shared_faces(const as::Massing& M) {
  const auto up = M.complex.cofaces();
  std::size_t n = 0;
  for (std::size_t f = 0; f < M.complex.size(); f++)
    if (M.complex.cells[f].dim == 2 && up[f].size() == 2) n++;
  return n;
}

}  // namespace

int main() {
  Suite s;
  {
    const auto M = as::massing({{"a", "A", {{0, 0}, {10, 0}, {10, 10}, {0, 10}}, 0, 10}, {"b", "B", {{10, 0}, {20, 0}, {20, 10}, {10, 10}}, 0, 20}});
    s.require("two boxes side by side: two 3-cells, regular", M.cell.size() == 2 && mt::admissible(M.complex).regular, mt::admissible(M.complex).reason);
    s.require("one wall between them, shared; the taller's wall above it its own", shared_faces(M) == 1);
    s.require("∂∂ = 0 in every degree", boundary_squares_to_zero(M.complex));
    s.require("their union a ball: Betti 1 0 0 0", text(morphe::algebra::betti(mt::chain_complex(M.complex))) == "1 0 0 0");
    s.check("volumes", M.volume[0] + M.volume[1] - (1000 + 2000), 1e-9, "m³");
  }
  {
    const auto M = as::massing({{"a", "A", {{0, 0}, {10, 0}, {10, 10}, {0, 10}}, 0, 10},
                                {"b", "B", {{10, 0}, {20, 0}, {20, 10}, {10, 10}}, 0, 20},
                                {"c", "C", {{0, 10}, {20, 10}, {20, 20}, {0, 20}}, 0, 15}});
    s.require("a third along both their backs: its edge cut at the corner between them", M.ring[2].size() == 5, std::to_string(M.ring[2].size()));
    s.require("still regular, three 3-cells", M.cell.size() == 3 && mt::admissible(M.complex).regular, mt::admissible(M.complex).reason);
    s.require("three shared walls: A with B, C with A, C with B", shared_faces(M) == 3, std::to_string(shared_faces(M)));
    s.require("∂∂ = 0, and the union still a ball", boundary_squares_to_zero(M.complex) && text(morphe::algebra::betti(mt::chain_complex(M.complex))) == "1 0 0 0");
  }
  {
    const auto M = as::massing({{"bow", "", {{0, 0}, {10, 10}, {10, 0}, {0, 10}}, 0, 10},
                                {"cw", "", {{0, 0}, {0, 10}, {10, 10}, {10, 0}}, 0, 10},
                                {"twice", "", {{0, 0}, {10, 0}, {5, 5}, {10, 10}, {0, 10}, {5, 5}}, 0, 10},
                                {"ok", "", {{30, 0}, {40, 0}, {40, 10}, {30, 10}}, 0, 10}});
    s.require("a bow tie, a clockwise footprint and one with a repeated corner skipped; the rest kept", M.cell.size() == 1 && M.skipped.size() == 3, std::to_string(M.skipped.size()));
  }
  return s.done();
}
