// Checks of the space module (rooms.hpp): a building of 4 × 2 rooms on 3 floors.
// - a ball (Betti 1 0 0 0), 24 rooms, regular;
// - a stair opened: the slab relabelled, nothing else made or removed;
// - two rooms merged across a wall: 23 rooms, regular, the building's heat unchanged; parted back: 24 rooms,
//   isomorphic to the building before, the heat unchanged again, each part at the merged temperature;
// - from the entrance: one room, then the whole ground floor once its doors are opened, then every floor once a
//   stair joins each two.

#include <cmath>
#include <string>
#include <vector>

#include <architecture/space/rooms.hpp>
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

}  // namespace

int main() {
  Suite s;
  auto Bd = as::building({{0, 0}, 0, 20, 12, 3.5, 4, 2, 3});
  s.require("4 × 2 rooms on 3 floors: 24 rooms, a ball, regular",
            Bd.capacity.size() == 24 && text(morphe::algebra::betti(mt::chain_complex(Bd.complex))) == "1 0 0 0" && mt::admissible(Bd.complex).regular);
  // Rooms at different temperatures, so that merging has something to average.
  {
    double t = 12;
    for (auto& [room, q] : Bd.heat) q = Bd.capacity.at(room) * (t += 0.5);
  }
  const std::uint64_t slab = [&] { for (const auto& c : Bd.complex.cells) if (c.type == "slab") return c.id; return std::uint64_t{0}; }();
  const double before = Bd.total_heat();
  const auto st = as::relabel_face(Bd, slab, "stair");
  s.require("a stair opened: the slab relabelled, nothing else", Bd.complex.cells[Bd.complex.find(slab)].type == "stair" && st.created.empty() && st.deleted.empty());
  const std::uint64_t wall = [&] { for (const auto& c : Bd.complex.cells) if (c.type == "wall" && Bd.rooms_of(c.id).size() == 2) return c.id; return std::uint64_t{0}; }();
  const auto pair = Bd.rooms_of(wall);
  const auto first = Bd.units.at(pair[0]);
  const double Tm = (Bd.heat.at(pair[0]) + Bd.heat.at(pair[1])) / (Bd.capacity.at(pair[0]) + Bd.capacity.at(pair[1]));
  const auto original = Bd.complex;
  const auto m = as::merge_rooms(Bd, wall);
  const std::uint64_t merged = m.created[0];
  s.require("two rooms merged: 23 rooms, regular", Bd.capacity.size() == 23 && mt::admissible(Bd.complex).regular);
  s.check("the building's heat unchanged by the merge (relative)", (Bd.total_heat() - before) / before, 1e-14);
  s.check("the merged room at the temperature the capacities weigh", Bd.temperature(merged) - Tm, 1e-12, "K");
  const auto back = as::partition_room(Bd, merged, first);
  s.require("parted back: 24 rooms, isomorphic to the building before", Bd.capacity.size() == 24 && mt::isomorphic(Bd.complex, original));
  s.check("the heat unchanged again (relative)", (Bd.total_heat() - before) / before, 1e-14);
  s.check("each part at the merged temperature", std::fabs(Bd.temperature(back.created[1]) - Tm) + std::fabs(Bd.temperature(back.created[2]) - Tm), 1e-12, "K");

  auto E = as::building({{0, 0}, 0, 20, 12, 3.5, 4, 2, 2});
  const std::uint64_t facade = [&] { for (const auto& c : E.complex.cells) if (c.type == "facade" && E.corner.at(c.id)[2] == 0 && E.corner.at(c.id)[1] == 0 && E.corner.at(c.id)[0] == 0 && E.span.at(c.id)[0]) return c.id; return std::uint64_t{0}; }();
  as::relabel_face(E, facade, "entrance");
  const std::size_t alone = as::reachable_rooms(E).size();
  for (const auto& c : std::vector<mt::RegularCell>(E.complex.cells))
    if (c.type == "wall" && E.corner.at(c.id)[2] == 0) as::relabel_face(E, c.id, "door");
  const std::size_t ground = as::reachable_rooms(E).size();
  for (const auto& c : std::vector<mt::RegularCell>(E.complex.cells))
    if (c.type == "wall") as::relabel_face(E, c.id, "door");
  for (const auto& c : std::vector<mt::RegularCell>(E.complex.cells))
    if (c.type == "slab" && E.corner.at(c.id)[0] == 0 && E.corner.at(c.id)[1] == 0) as::relabel_face(E, c.id, "stair");
  const std::size_t all = as::reachable_rooms(E).size();
  s.require("from the entrance: one room, then the ground floor once its doors are opened, then both floors by a stair",
            alone == 1 && ground == 8 && all == 16, std::to_string(alone) + " " + std::to_string(ground) + " " + std::to_string(all));
  return s.done();
}
