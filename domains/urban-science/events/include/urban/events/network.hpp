// A district's volumes and interfaces (the design, §1.3; ADR 0067 §2.3): the spin network of the moment of a city,
// the dual of its decomposition of space. Couplings are volumes, primitives the interfaces between two of them.
//
//   couplings   each building of the massing (architecture's `space/massing.hpp`, a 3-cell); each street space, one
//               for each segment of the district (urban-science's `landscape/district.hpp`, a 1-cell, its kind kept:
//               street, laneway, mall, lane); each junction space (a 0-cell).
//   primitives  party walls: the wall pieces two buildings share (2-cells with two cofaces), one primitive for each
//               pair of buildings, its area their sum and its normal their area-weighted outward normal from the first
//               building toward the second;
//               frontages: the outer wall pieces of a building (2-cells with one coface) that face a street space:
//               the nearest segment to the piece's middle within `frontage_within` metres, toward which the piece's
//               outward normal points (the cosine of the angle at least `facing`); one primitive for each building
//               and street, its area the pieces' sum;
//               cross-sections: a street space and each of its two junctions.
// The roofs, the ground and the facades that face no street are interfaces with the outside, which is not yet a
// coupling: no event of the record reaches it, and it would join every building to every other.
//
// Each primitive keeps its shape for drawing: the base lines of its wall pieces, or the junction's point.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <architecture/space/massing.hpp>
#include <physis/computational_geometry/planar.hpp>
#include <urban/events/record.hpp>
#include <urban/landscape/district.hpp>

namespace urban::events {

namespace as = architecture::space;
namespace cg = physis::computational_geometry;

struct NetworkOptions {
  double frontage_within = 25;  // m
  double facing = 0.25;         // the least cosine between a facade's outward normal and the way to its street
};

struct CityNetwork {
  std::vector<Coupling> couplings;
  std::vector<Primitive> primitives;
  std::vector<std::vector<std::array<double, 3>>> shape;  // each primitive's: wall pieces' base lines in pairs, or a point
  std::map<std::string, std::size_t> building;            // a footprint's key ↦ its coupling
  std::map<std::uint64_t, std::size_t> street, junction;  // a district cell's id ↦ its coupling
  std::size_t buildings = 0, streets = 0, junctions = 0;
  std::size_t party_walls = 0, frontages = 0, cross_sections = 0;
  std::size_t facades = 0, facades_fronting = 0;           // outer wall pieces, and those assigned to a street
};

namespace detail {

inline std::array<double, 2> centroid(const std::vector<std::array<double, 2>>& r) {
  double a = 0, x = 0, y = 0;
  for (std::size_t i = 0; i < r.size(); i++) {
    const auto& p = r[i];
    const auto& q = r[(i + 1) % r.size()];
    const double w = p[0] * q[1] - q[0] * p[1];
    a += w, x += (p[0] + q[0]) * w, y += (p[1] + q[1]) * w;
  }
  if (std::fabs(a) < 1e-12) {
    std::array<double, 2> m{};
    for (const auto& p : r) m[0] += p[0] / static_cast<double>(r.size()), m[1] += p[1] / static_cast<double>(r.size());
    return m;
  }
  return {x / (3 * a), y / (3 * a)};
}

inline std::vector<cg::Point2> line_of(const landscape::District& D, std::uint64_t e) {
  const auto it = D.line.find(e);
  if (it != D.line.end() && it->second.size() >= 2) return it->second;
  const auto& c = D.cell(e);
  const auto ends = D.complex.ends(D.at(e));
  (void)c;
  return {D.point.at(D.complex.cells[ends[0]].id), D.point.at(D.complex.cells[ends[1]].id)};
}

}  // namespace detail

inline CityNetwork city_network(const as::Massing& M, const std::vector<as::Footprint>& F, const landscape::District& D, NetworkOptions o = {}) {
  CityNetwork N;
  // The buildings.
  std::map<std::uint64_t, std::size_t> of_volume;  // a massing 3-cell id ↦ its coupling
  for (std::size_t k = 0; k < M.cell.size(); k++) {
    const auto& f = F.at(M.building[k]);
    Coupling c;
    c.id = "building:" + f.key, c.type = "building", c.name = f.name;
    c.at = detail::centroid(M.ring[k].empty() ? f.ring : M.ring[k]);
    c.base = f.base, c.top = f.top;
    of_volume[M.cell[k]] = N.couplings.size();
    N.building[f.key] = N.couplings.size();
    N.couplings.push_back(c);
  }
  N.buildings = N.couplings.size();
  // The junction spaces and the street spaces.
  for (const auto& cell : D.complex.cells) {
    if (cell.dim != 0) continue;
    Coupling c;
    c.id = "junction:" + std::to_string(cell.id), c.type = "junction";
    c.at = D.point.at(cell.id);
    N.junction[cell.id] = N.couplings.size();
    N.couplings.push_back(c);
  }
  N.junctions = N.junction.size();
  std::vector<std::pair<std::uint64_t, std::vector<cg::Point2>>> lines;
  for (const auto& cell : D.complex.cells) {
    if (cell.dim != 1) continue;
    const auto line = detail::line_of(D, cell.id);
    Coupling c;
    c.id = "street:" + std::to_string(cell.id), c.type = cell.type;
    const auto nm = D.name.find(cell.id);
    c.name = nm == D.name.end() ? std::string() : nm->second;
    c.at = cg::point_at(line, cg::polyline_length(line) / 2);
    N.street[cell.id] = N.couplings.size();
    N.couplings.push_back(c);
    lines.push_back({cell.id, line});
  }
  N.streets = N.street.size();
  // The cross-sections.
  for (const auto& [e, line] : lines) {
    const auto ends = D.complex.ends(D.at(e));
    for (std::size_t k = 0; k < 2; k++) {
      const auto j = D.complex.cells[ends[k]].id;
      Primitive p;
      p.id = "cross:" + std::to_string(e) + ":" + std::to_string(j), p.type = "cross-section";
      p.a = N.street.at(e), p.b = N.junction.at(j);
      const auto& at = D.point.at(j);
      N.primitives.push_back(p);
      N.shape.push_back({{at[0], at[1], 0}});
    }
  }
  N.cross_sections = N.primitives.size();
  // The massing's walls: shared ones are party walls, outer ones may front a street.
  const auto cof = M.complex.cofaces();
  struct Pool {
    double area = 0;
    std::array<double, 3> normal{};
    std::vector<std::array<double, 3>> shape;
  };
  std::map<std::pair<std::size_t, std::size_t>, Pool> party, front;
  const auto sign_in = [&](std::size_t volume, std::size_t face) {
    const auto& c = M.complex.cells[volume];
    for (std::size_t q = 0; q < c.facets.size(); q++)
      if (c.facets[q] == face) return c.signs[q] >= 0 ? 1.0 : -1.0;
    return 1.0;
  };
  for (std::size_t i = 0; i < M.complex.cells.size(); i++) {
    const auto& cell = M.complex.cells[i];
    if (cell.dim != 2) continue;
    const auto fit = M.face.find(cell.id);
    if (fit == M.face.end() || fit->second.kind != "wall") continue;
    const auto& w = fit->second;
    std::vector<std::size_t> vols;
    for (auto c : cof[i])
      if (M.complex.cells[c].dim == 3) vols.push_back(c);
    const double z = (w.z0 + w.z1) / 2;
    if (vols.size() == 2) {
      std::size_t x = of_volume.at(M.complex.cells[vols[0]].id), y = of_volume.at(M.complex.cells[vols[1]].id);
      double s = sign_in(vols[0], i);
      if (x > y) std::swap(x, y), s = sign_in(vols[1], i);
      auto& P = party[{x, y}];
      P.area += w.area;
      for (int d = 0; d < 3; d++) P.normal[d] += s * w.normal[d] * w.area;
      P.shape.push_back({w.a[0], w.a[1], z}), P.shape.push_back({w.b[0], w.b[1], z});
    } else if (vols.size() == 1) {
      N.facades++;
      const double s = sign_in(vols[0], i);
      const double nx = s * w.normal[0], ny = s * w.normal[1], nh = std::hypot(nx, ny);
      if (nh < 1e-9) continue;
      const cg::Point2 m{(w.a[0] + w.b[0]) / 2, (w.a[1] + w.b[1]) / 2};
      double best = o.frontage_within;
      long which = -1;
      for (std::size_t q = 0; q < lines.size(); q++) {
        const auto [at, dist] = cg::nearest_on_polyline(lines[q].second, m);
        if (dist > best || dist < 1e-9) continue;
        const auto p = cg::point_at(lines[q].second, at);
        const double dx = p[0] - m[0], dy = p[1] - m[1];
        if ((nx * dx + ny * dy) / (nh * dist) < o.facing) continue;
        best = dist, which = static_cast<long>(q);
      }
      if (which < 0) continue;
      N.facades_fronting++;
      const std::size_t x = of_volume.at(M.complex.cells[vols[0]].id), y = N.street.at(lines[static_cast<std::size_t>(which)].first);
      auto& P = front[{x, y}];
      P.area += w.area;
      for (int d = 0; d < 3; d++) P.normal[d] += s * w.normal[d] * w.area;
      P.shape.push_back({w.a[0], w.a[1], z}), P.shape.push_back({w.b[0], w.b[1], z});
    }
  }
  const auto add = [&](const std::map<std::pair<std::size_t, std::size_t>, Pool>& pools, const std::string& type) {
    for (const auto& [ab, P] : pools) {
      Primitive p;
      p.type = type, p.a = ab.first, p.b = ab.second, p.area = P.area;
      p.id = (type == "party-wall" ? "party:" : "front:") + N.couplings[ab.first].id + "|" + N.couplings[ab.second].id;
      const double n = std::sqrt(P.normal[0] * P.normal[0] + P.normal[1] * P.normal[1] + P.normal[2] * P.normal[2]);
      if (n > 0) p.normal = {P.normal[0] / n, P.normal[1] / n, P.normal[2] / n};
      N.primitives.push_back(p);
      N.shape.push_back(P.shape);
    }
  };
  add(party, "party-wall");
  N.party_walls = party.size();
  add(front, "frontage");
  N.frontages = front.size();
  return N;
}

}  // namespace urban::events
