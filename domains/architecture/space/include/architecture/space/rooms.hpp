// A building as a complex of rooms (the infrastructure's ADR 0062; the urban-science domain's
// theory/landscape-dynamics.md §7): a box of rooms fitted to a building, its rooms the 3-cells, the slabs, walls and
// facades between them the 2-cells, typed; what each room carries; and the events that change it — a stair or a door
// opened, two rooms merged and parted again — with the rooms one can reach from its entrance.
//
// The complex. A box of nx × ny rooms on each of `floors` floors, the cubical grid of morphe regular_complex.hpp,
// its 2-cells typed by where they lie: horizontal ones "ground" at the bottom, "roof" at the top, "slab" between
// floors; vertical ones "facade" on the outside, "wall" between rooms. It is a ball: b₀ = 1, the rest 0. Each face's
// fabric (its transmittance U by its type) is part of the building, read by the environment module.
//
// What a room carries (the schema): its heat Q (J) and capacity C (J/K), both extensive, its temperature Q/C
// intensive. Events carry them as the schema says: summed when rooms merge, shared in proportion when a room is
// parted, so that the building's heat is unchanged by every event.
//
// Events. A stair opened is a slab relabelled "stair" (an opening conducts as air moves: U = 25 W/m²K), a door a
// wall relabelled "door": the complex keeps its cells, only what they carry and mean changes. Two rooms merged
// across a wall is merge_cells: one 3-cell with the facets of both but the wall, its heat the sum of theirs and its
// capacity the sum; parted back along the first room's facets it has them again, each with its share of capacity
// at the merged temperature. Every derived complex is regular. The rooms reached from the entrance are found by a
// search through doors and stairs. Checks: tests/check_space.cpp.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <morphe/formal/traces.hpp>
#include <morphe/topology/complex_rewriting.hpp>
#include <morphe/topology/regular_complex.hpp>

namespace architecture::space {

namespace mt = morphe::topology;
using Point = std::array<double, 2>;

// One rewrite of a building: its kind, what it matched, the cells it made and removed, what it read and wrote, and
// whether it describes the same building more finely (a refinement) or changes it (an intervention).
struct Step {
  std::string kind, match;
  std::vector<std::uint64_t> created, deleted;
  morphe::formal::traces::Footprint footprint;
  bool refinement = false;
};

struct BuildingSpec {
  Point origin{};      // local east, north of the corner (0, 0) of the box (m)
  double heading = 0;  // the box's first axis, radians counter-clockwise from east
  double width = 20, depth = 12, storey = 3.5;  // m
  int nx = 4, ny = 2, floors = 3;
};

struct Fabric {
  std::map<std::string, double> u{{"facade", 1.2}, {"roof", 0.6}, {"ground", 0.4}, {"wall", 2.5}, {"slab", 1.8}, {"stair", 25.0}, {"door", 8.0}, {"entrance", 8.0}};
  double capacity_per_m3 = 6.0e4;  // J/(m³ K)
  double absorptance = 0.5;         // of the opaque fabric's outer surface
  double outer_coefficient = 25;    // h_o, W/(m² K): the outer surface to the air
  double glazing = 0.3;             // the glazed share of a facade
  double gain_coefficient = 0.6;    // the glass's solar heat gain coefficient
};

struct Building {
  BuildingSpec spec;
  Fabric fabric;
  mt::RegularComplex complex;
  std::map<std::uint64_t, std::array<int, 3>> corner, span;     // each cell of the grid as built
  std::map<std::uint64_t, std::vector<std::array<int, 3>>> units;  // each room's unit cubes (more than one once merged)
  std::map<std::uint64_t, double> heat, capacity;                // each room's Q (J) and C (J/K)
  std::map<std::uint64_t, std::vector<std::uint64_t>> merged_from; // a merged room ↦ the facets of its first room

  double dx() const { return spec.width / spec.nx; }
  double dy() const { return spec.depth / spec.ny; }
  // A face's area (m²) from its span.
  double area(std::uint64_t face) const {
    const auto& s = span.at(face);
    return (s[0] ? dx() : 1) * (s[1] ? dy() : 1) * (s[2] ? spec.storey : 1);
  }
  // The rooms on each side of a face (one for an outer face).
  std::vector<std::uint64_t> rooms_of(std::uint64_t face) const {
    std::vector<std::uint64_t> out;
    const auto up = complex.cofaces()[complex.find(face)];
    for (std::size_t c : up) out.push_back(complex.cells[c].id);
    return out;
  }
  double temperature(std::uint64_t room) const { return heat.at(room) / capacity.at(room); }
  double total_heat() const {
    double s = 0;
    for (const auto& [r, q] : heat) s += q;
    return s;
  }
  // A face's outward normal in the local frame (east, north, up), for an outer face; its world centre.
  std::array<double, 3> normal(std::uint64_t face) const {
    const auto& c = corner.at(face);
    const auto& s = span.at(face);
    std::array<double, 3> n{0, 0, 0};  // in the box's axes
    if (!s[0]) n[0] = c[0] == 0 ? -1 : 1;
    else if (!s[1]) n[1] = c[1] == 0 ? -1 : 1;
    else n[2] = c[2] == 0 ? -1 : 1;
    const double ch = std::cos(spec.heading), sh = std::sin(spec.heading);
    return {n[0] * ch - n[1] * sh, n[0] * sh + n[1] * ch, n[2]};
  }
  std::array<double, 3> world(double x, double y, double z) const {
    const double ch = std::cos(spec.heading), sh = std::sin(spec.heading);
    return {spec.origin[0] + x * ch - y * sh, spec.origin[1] + x * sh + y * ch, z};
  }
};

inline Building building(const BuildingSpec& spec, const Fabric& fabric = {}) {
  if (spec.nx < 1 || spec.ny < 1 || spec.floors < 1) throw std::invalid_argument("urban::landscape::building: at least one room");
  if (spec.nx * spec.ny * spec.floors > 160) throw std::invalid_argument("urban::landscape::building: at most 160 rooms");
  Building B;
  B.spec = spec, B.fabric = fabric;
  const auto G = mt::cubical_grid(spec.nx, spec.ny, spec.floors, "building");
  B.complex = G.complex;
  for (std::size_t i = 0; i < B.complex.size(); i++) {
    auto& c = B.complex.cells[i];
    const auto co = G.corner[i], sp = G.span[i];
    B.corner[c.id] = co, B.span[c.id] = sp;
    if (c.dim == 3) {
      c.type = "room";
      B.units[c.id] = {co};
      const double C = fabric.capacity_per_m3 * B.dx() * B.dy() * spec.storey;
      B.capacity[c.id] = C;
      B.heat[c.id] = C * 20.0;
    } else if (c.dim == 2) {
      if (!sp[2]) c.type = co[2] == 0 ? "ground" : co[2] == spec.floors ? "roof" : "slab";
      else if (!sp[0]) c.type = (co[0] == 0 || co[0] == spec.nx) ? "facade" : "wall";
      else c.type = (co[1] == 0 || co[1] == spec.ny) ? "facade" : "wall";
    } else if (c.dim == 1) c.type = "edge";
    else c.type = "corner";
  }
  return B;
}

// ---- events ----------------------------------------------------------------------------------------------------

// A face relabelled: a slab as a stair, a wall as a door, a facade as the entrance.
inline Step relabel_face(Building& B, std::uint64_t face, const std::string& type) {
  const auto& c = B.complex.cells.at(B.complex.find(face));
  const std::map<std::string, std::string> allowed{{"stair", "slab"}, {"door", "wall"}, {"entrance", "facade"}};
  const auto it = allowed.find(type);
  if (it == allowed.end() || c.type != it->second) throw std::invalid_argument("urban::landscape::relabel_face: a " + type + " is opened in a " + (it == allowed.end() ? std::string("?") : it->second) + ", not a " + c.type);
  const auto r = mt::rewrite(B.complex, mt::relabel(B.complex, B.complex.find(face), type));
  if (!r.applied) throw std::invalid_argument("urban::landscape::relabel_face: " + r.reason);
  B.complex = r.result;
  return {"relabel", "a " + it->second + " opened as a " + type, r.created, r.deleted, r.footprint};
}

// Two rooms merged across an interior face (a wall or a door): heat and capacity summed.
inline Step merge_rooms(Building& B, std::uint64_t face) {
  const auto rooms = B.rooms_of(face);
  if (rooms.size() != 2) throw std::invalid_argument("urban::landscape::merge_rooms: the face does not lie between two rooms");
  const auto type = B.complex.cells.at(B.complex.find(face)).type;
  if (type != "wall" && type != "door") throw std::invalid_argument("urban::landscape::merge_rooms: rooms are merged across a wall or a door, not a " + type);
  std::vector<std::uint64_t> first;
  for (std::size_t f : B.complex.cells[B.complex.find(rooms[0])].facets)
    if (B.complex.cells[f].id != face) first.push_back(B.complex.cells[f].id);
  const auto r = mt::rewrite(B.complex, mt::merge_cells(B.complex, B.complex.find(face)));
  if (!r.applied) throw std::invalid_argument("urban::landscape::merge_rooms: " + r.reason);
  B.complex = r.result;
  const std::uint64_t m = r.created[0];
  B.heat[m] = B.heat.at(rooms[0]) + B.heat.at(rooms[1]);
  B.capacity[m] = B.capacity.at(rooms[0]) + B.capacity.at(rooms[1]);
  auto u = B.units.at(rooms[0]);
  const auto& u2 = B.units.at(rooms[1]);
  u.insert(u.end(), u2.begin(), u2.end());
  B.units[m] = u;
  B.merged_from[m] = first;
  for (auto id : rooms) B.heat.erase(id), B.capacity.erase(id), B.units.erase(id);
  return {"merge", "two rooms, across a " + type, r.created, r.deleted, r.footprint};
}

// A merged room partitioned back along its first room's facets: each part keeps the merged temperature, with its
// share of capacity by volume (so the heat is unchanged).
inline Step partition_room(Building& B, std::uint64_t room, const std::vector<std::array<int, 3>>& first_units) {
  const auto from = B.merged_from.find(room);
  if (from == B.merged_from.end()) throw std::invalid_argument("urban::landscape::partition_room: the room was not merged");
  std::vector<std::size_t> part;
  for (auto id : from->second) part.push_back(B.complex.find(id));
  const double T = B.temperature(room), C = B.capacity.at(room);
  const auto all = B.units.at(room);
  const auto r = mt::rewrite(B.complex, mt::partition_cell(B.complex, B.complex.find(room), part, "wall"));
  if (!r.applied) throw std::invalid_argument("urban::landscape::partition_room: " + r.reason);
  B.complex = r.result;
  const std::uint64_t wall = r.created[0], r1 = r.created[1], r2 = r.created[2];
  std::vector<std::array<int, 3>> second;
  for (const auto& x : all)
    if (std::find(first_units.begin(), first_units.end(), x) == first_units.end()) second.push_back(x);
  const double share = static_cast<double>(first_units.size()) / all.size();
  B.capacity[r1] = C * share, B.capacity[r2] = C * (1 - share);
  B.heat[r1] = B.capacity[r1] * T, B.heat[r2] = B.capacity[r2] * T;
  B.units[r1] = first_units, B.units[r2] = second;
  // The new wall lies where the old one did: the face between the two unit sets, recovered from the first room's
  // unit cube next to the second's.
  for (const auto& a : first_units)
    for (const auto& b : second) {
      int diff = 0, axis = -1;
      for (int k = 0; k < 3; k++)
        if (a[k] != b[k]) diff += std::abs(a[k] - b[k]), axis = k;
      if (diff == 1) {
        std::array<int, 3> c = a, s{1, 1, 1};
        c[axis] = std::max(a[axis], b[axis]);
        s[axis] = 0;
        B.corner[wall] = c, B.span[wall] = s;
      }
    }
  for (auto id : {r1, r2}) B.merged_from.erase(id);
  B.merged_from.erase(room);
  B.heat.erase(room), B.capacity.erase(room), B.units.erase(room);
  return {"partition", "a room, back into two", r.created, r.deleted, r.footprint};
}

// The rooms reachable from the entrance through doors, stairs and merged rooms.
inline std::set<std::uint64_t> reachable_rooms(const Building& B) {
  std::set<std::uint64_t> seen;
  std::queue<std::uint64_t> q;
  const auto up = B.complex.cofaces();
  for (std::size_t i = 0; i < B.complex.size(); i++)
    if (B.complex.cells[i].type == "entrance" && up[i].size() == 1) {
      const auto r = B.complex.cells[up[i][0]].id;
      if (seen.insert(r).second) q.push(r);
    }
  while (!q.empty()) {
    const auto r = q.front();
    q.pop();
    for (std::size_t f : B.complex.cells[B.complex.find(r)].facets) {
      const auto& fc = B.complex.cells[f];
      if (fc.type != "door" && fc.type != "stair") continue;
      for (std::size_t c : up[f]) {
        const auto other = B.complex.cells[c].id;
        if (seen.insert(other).second) q.push(other);
      }
    }
  }
  return seen;
}

}  // namespace architecture::space
