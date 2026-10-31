// The buildings of a district as one regular 3-complex (the infrastructure's ADR 0063): each building a 3-cell
// bounded by its ground face, its roof and its walls; two buildings standing against each other share the wall
// between them, one 2-cell with two cofaces, so that the buildings of a street are joined into terraces as the
// streets join blocks. Buildings are space, as the streets are.
//
// The construction. Corners are pooled (footprints given in a local frame, a corner within `snap` of another the
// same corner), and a corner lying on another footprint's edge is inserted into it (a T-junction), so that two
// footprints sharing a stretch of edge share its pieces exactly. Each building i stands from its base b_i to its top
// t_i. On a pooled edge e the heights at which something begins or ends are L_e = ∪ {b_i, t_i : i uses e}; at a
// corner v, L_v = ∪_{e ∋ v} L_e. The cells:
//   0-cells   (v, z), z ∈ L_v, made when a face needs them;
//   1-cells   horizontal (e, z), z ∈ L_e, oriented from the edge's first corner to its second; vertical (v, z, z′) for
//             consecutive z < z′ in L_v, oriented upward;
//   2-cells   for each building its ground face at b_i and its roof at t_i, walked counter-clockwise from above (their
//             normals up); for each pooled edge e = (a → b) and consecutive z < z′ in L_e its wall piece, walked
//             (a,z) → (b,z) → up b → (b,z′) → (a,z′) → down a, so that its normal points to the right of a → b;
//   3-cells   each building, its faces with the signs that orient them outward: the roof +1, the ground −1, a wall
//             piece +1 where the building's counter-clockwise ring runs a → b (its outside to the right), −1 where it
//             runs b → a. A wall piece two buildings share is +1 for one and −1 for the other.
// A footprint that is not a simple polygon of three or more corners after pooling is skipped with the reason.
// Courtyards (holes) are not carved: the ground face would be an annulus, not a cell.
//
// Checks: tests/check_massing.cpp.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <morphe/topology/regular_complex.hpp>

namespace architecture::space {

namespace mt = morphe::topology;

struct Footprint {
  std::string key, name;
  std::vector<std::array<double, 2>> ring;  // counter-clockwise, local east and north (m), not closed
  double base = 0, top = 10;                // m
};

struct Massing {
  mt::RegularComplex complex;
  std::map<std::uint64_t, std::array<double, 3>> at;   // each 0-cell's position
  std::vector<std::size_t> building;                   // input index of each building kept, in 3-cell order
  std::map<std::uint64_t, std::size_t> of_cell;        // a 3-cell id ↦ its index in `building`
  std::vector<std::uint64_t> cell;                     // the 3-cell id of each building kept
  std::vector<double> volume, footprint_area;
  std::vector<std::vector<std::array<double, 2>>> ring;  // each building's ring after pooling and T-junctions
  // Each 2-cell's kind ("ground", "roof", "wall"), area (m²), unit normal (as walked) and, for walls, its edge's
  // two corners and heights.
  struct Face {
    std::string kind;
    double area = 0;
    std::array<double, 3> normal{};
    std::array<double, 2> a{}, b{};
    double z0 = 0, z1 = 0;
  };
  std::map<std::uint64_t, Face> face;
  std::vector<std::pair<std::size_t, std::string>> skipped;  // input index, reason
};

namespace detail {

inline double ring_area(const std::vector<std::array<double, 2>>& r) {
  double s = 0;
  for (std::size_t i = 0; i < r.size(); i++) {
    const auto& p = r[i];
    const auto& q = r[(i + 1) % r.size()];
    s += p[0] * q[1] - q[0] * p[1];
  }
  return s / 2;
}

}  // namespace detail

inline Massing massing(const std::vector<Footprint>& F, double snap = 0.05) {
  Massing M;
  // Pooled corners.
  std::vector<std::array<double, 2>> corner;
  std::map<std::pair<long long, long long>, std::size_t> pool;
  const auto pooled = [&](std::array<double, 2> p) {
    const long long kx = std::llround(p[0] / snap), ky = std::llround(p[1] / snap);
    for (long long dx = -1; dx <= 1; dx++)
      for (long long dy = -1; dy <= 1; dy++) {
        const auto it = pool.find({kx + dx, ky + dy});
        if (it != pool.end() && std::hypot(corner[it->second][0] - p[0], corner[it->second][1] - p[1]) <= snap) return it->second;
      }
    pool[{kx, ky}] = corner.size();
    corner.push_back(p);
    return corner.size() - 1;
  };
  std::vector<std::vector<std::size_t>> rings(F.size());
  std::vector<char> keep(F.size(), 1);
  for (std::size_t i = 0; i < F.size(); i++) {
    if (!(F[i].top > F[i].base)) { M.skipped.push_back({i, "its top is not above its base"}); keep[i] = 0; continue; }
    std::vector<std::size_t> r;
    for (const auto& p : F[i].ring) {
      const auto k = pooled(p);
      if (r.empty() || r.back() != k) r.push_back(k);
    }
    while (r.size() > 1 && r.front() == r.back()) r.pop_back();
    rings[i] = r;
  }
  // T-junctions: corners lying on another footprint's edge are inserted into it.
  const double tol = std::max(snap, 0.02);
  for (std::size_t i = 0; i < F.size(); i++) {
    if (!keep[i]) continue;
    std::vector<std::size_t> out;
    const auto& r = rings[i];
    for (std::size_t k = 0; k < r.size(); k++) {
      const std::size_t a = r[k], b = r[(k + 1) % r.size()];
      out.push_back(a);
      const auto& A = corner[a];
      const auto& B = corner[b];
      const double dx = B[0] - A[0], dy = B[1] - A[1], L2 = dx * dx + dy * dy;
      if (L2 <= 0) continue;
      std::vector<std::pair<double, std::size_t>> on;
      for (std::size_t c = 0; c < corner.size(); c++) {
        if (c == a || c == b) continue;
        const auto& C = corner[c];
        if (C[0] < std::min(A[0], B[0]) - tol || C[0] > std::max(A[0], B[0]) + tol || C[1] < std::min(A[1], B[1]) - tol || C[1] > std::max(A[1], B[1]) + tol) continue;
        const double t = ((C[0] - A[0]) * dx + (C[1] - A[1]) * dy) / L2;
        if (t <= 1e-9 || t >= 1 - 1e-9) continue;
        const double px = A[0] + t * dx - C[0], py = A[1] + t * dy - C[1];
        if (px * px + py * py <= tol * tol) on.push_back({t, c});
      }
      std::sort(on.begin(), on.end());
      for (const auto& [t, c] : on) out.push_back(c);
    }
    rings[i] = out;
    // Simple: no corner twice, at least three, a positive area.
    std::set<std::size_t> seen(out.begin(), out.end());
    std::vector<std::array<double, 2>> pts;
    for (auto c : out) pts.push_back(corner[c]);
    if (out.size() < 3 || seen.size() != out.size()) { M.skipped.push_back({i, "its footprint is not a simple polygon"}); keep[i] = 0; continue; }
    if (!(detail::ring_area(pts) > 0.5)) { M.skipped.push_back({i, "its footprint is not counter-clockwise, or has no area"}); keep[i] = 0; continue; }
  }
  // Heights on each pooled edge and at each corner.
  using Edge = std::pair<std::size_t, std::size_t>;  // (smaller, larger) corner index
  const auto edge_of = [](std::size_t a, std::size_t b) { return a < b ? Edge{a, b} : Edge{b, a}; };
  std::map<Edge, std::set<double>> L_e;
  std::map<std::size_t, std::set<double>> L_v;
  for (std::size_t i = 0; i < F.size(); i++) {
    if (!keep[i]) continue;
    const auto& r = rings[i];
    for (std::size_t k = 0; k < r.size(); k++) {
      const auto e = edge_of(r[k], r[(k + 1) % r.size()]);
      L_e[e].insert(F[i].base), L_e[e].insert(F[i].top);
    }
  }
  for (const auto& [e, zs] : L_e)
    for (double z : zs) L_v[e.first].insert(z), L_v[e.second].insert(z);
  // Cells made on demand.
  auto& X = M.complex;
  X.name = "massing";
  std::map<std::pair<std::size_t, double>, std::size_t> vcell;
  const auto vertex = [&](std::size_t v, double z) {
    const auto key = std::pair{v, z};
    const auto it = vcell.find(key);
    if (it != vcell.end()) return it->second;
    const std::size_t c = X.add_vertex("corner");
    M.at[X.cells[c].id] = {corner[v][0], corner[v][1], z};
    return vcell[key] = c;
  };
  std::map<std::pair<Edge, double>, std::size_t> hcell;
  const auto horizontal = [&](Edge e, double z) {
    const auto key = std::pair{e, z};
    const auto it = hcell.find(key);
    if (it != hcell.end()) return it->second;
    return hcell[key] = X.add_edge(vertex(e.first, z), vertex(e.second, z), "horizontal");
  };
  std::map<std::tuple<std::size_t, double, double>, std::size_t> vertical_cell;
  // The vertical pieces at v from z to z′ (z < z′), bottom to top.
  const auto vertical = [&](std::size_t v, double z, double z2) {
    std::vector<std::size_t> out;
    const auto& zs = L_v.at(v);
    auto it = zs.find(z);
    if (it == zs.end()) throw std::logic_error("architecture::space::massing: a height missing at a corner");
    while (*it < z2) {
      const double lo = *it, hi = *std::next(it);
      const auto key = std::tuple{v, lo, hi};
      auto f = vertical_cell.find(key);
      if (f == vertical_cell.end()) f = vertical_cell.emplace(key, X.add_edge(vertex(v, lo), vertex(v, hi), "vertical")).first;
      out.push_back(f->second);
      ++it;
    }
    return out;
  };
  std::map<std::pair<Edge, double>, std::size_t> wall_cell;  // (edge, lower height) ↦ its wall piece
  const auto wall = [&](Edge e, double z, double z2) {
    const auto key = std::pair{e, z};
    const auto it = wall_cell.find(key);
    if (it != wall_cell.end()) return it->second;
    const auto [a, b] = e;
    std::vector<std::size_t> facets{horizontal(e, z)};
    std::vector<int> signs{1};
    for (auto c : vertical(b, z, z2)) facets.push_back(c), signs.push_back(1);
    facets.push_back(horizontal(e, z2)), signs.push_back(-1);
    for (auto c : vertical(a, z, z2)) facets.push_back(c), signs.push_back(-1);
    const std::size_t w = X.add(2, "wall", facets, signs);
    Massing::Face f;
    f.kind = "wall", f.a = corner[a], f.b = corner[b], f.z0 = z, f.z1 = z2;
    const double dx = corner[b][0] - corner[a][0], dy = corner[b][1] - corner[a][1], L = std::hypot(dx, dy);
    f.area = L * (z2 - z), f.normal = {dy / L, -dx / L, 0};
    M.face[X.cells[w].id] = f;
    return wall_cell[key] = w;
  };
  const auto flat = [&](const std::vector<std::size_t>& r, double z, const std::string& kind) {
    std::vector<std::size_t> facets;
    std::vector<int> signs;
    std::vector<std::array<double, 2>> pts;
    for (std::size_t k = 0; k < r.size(); k++) {
      const std::size_t a = r[k], b = r[(k + 1) % r.size()];
      facets.push_back(horizontal(edge_of(a, b), z)), signs.push_back(a < b ? 1 : -1);
      pts.push_back(corner[a]);
    }
    const std::size_t c = X.add(2, kind, facets, signs);
    Massing::Face f;
    f.kind = kind, f.area = detail::ring_area(pts), f.normal = {0, 0, 1}, f.z0 = f.z1 = z;
    M.face[X.cells[c].id] = f;
    return c;
  };
  for (std::size_t i = 0; i < F.size(); i++) {
    if (!keep[i]) continue;
    const auto& r = rings[i];
    std::vector<std::size_t> facets{flat(r, F[i].base, "ground"), flat(r, F[i].top, "roof")};
    std::vector<int> signs{-1, 1};
    for (std::size_t k = 0; k < r.size(); k++) {
      const std::size_t a = r[k], b = r[(k + 1) % r.size()];
      const auto e = edge_of(a, b);
      const auto& zs = L_e.at(e);
      for (auto it = zs.find(F[i].base); *it < F[i].top; ++it) facets.push_back(wall(e, *it, *std::next(it))), signs.push_back(a < b ? 1 : -1);
    }
    const std::size_t c = X.add(3, "building", facets, signs);
    std::vector<std::array<double, 2>> pts;
    for (auto v : r) pts.push_back(corner[v]);
    const double area = detail::ring_area(pts);
    M.of_cell[X.cells[c].id] = M.building.size();
    M.building.push_back(i), M.cell.push_back(X.cells[c].id);
    M.footprint_area.push_back(area), M.volume.push_back(area * (F[i].top - F[i].base)), M.ring.push_back(pts);
  }
  return M;
}

}  // namespace architecture::space
