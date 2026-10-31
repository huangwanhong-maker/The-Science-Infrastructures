// A district as a typed cellular configuration (theory/landscape-dynamics.md §2): its walking network as a
// regular cell complex — junctions the 0-cells, the segments of streets and lanes the 1-cells, the blocks they
// enclose the 2-cells — with the geometry of each cell, the stores attached to the segments they face, and the
// state each kind of cell carries by its schema; and a lane through a block, built and closed as rewrites of the
// complex with the transfer of state they need.
//
// The complex. The segments are a plane graph (OpenStreetMap's ways split at their shared nodes, in the local
// frame of the city); its faces (physis computational_geometry, planar.hpp) are traced, and every bounded face
// whose boundary is a circuit becomes a 2-cell, a block, oriented counter-clockwise. A face with a dead end or a
// bridge inside it is not a disc bounded by a circuit and stays a hole, counted, so that the complex is regular
// (morphe regular_complex.hpp) and its homology says how many holes there are: b₁ of the district is the number
// of bounded faces left open. Euler's relation V − E + F = 2c over the traced faces is checked, so a network that
// is not plane (ways crossing without a node) is refused rather than given a false complex.
//
// The schema (§2.3). A cell's type says what it carries: a junction its walkers (persons, extensive), a segment
// its flow (persons an hour, a 1-cochain, intensive along the segment), a block its temperature (°C, intensive)
// and heat capacity (J/K, extensive). A rewrite carries the state with it by the kind of each fibre: an
// extensive quantity is split in proportion (area, length) or summed, an intensive one copied, or averaged with
// the weights that keep its extensive partner (heat = capacity × temperature) conserved; a cell that is new and
// carries what nothing gives it starts at zero (a new junction has no walkers yet).
//
// The lane (§3.4). A passage through a block between two of its streets, at the points facing two stores, is
// three rewrites: each street subdivided at its point (subdivide_edge, the cofaces keeping e ↦ e₁ + e₂), and
// the block split between the two new junctions (split_face). The two subdivisions are independent and the
// split depends on both: the record's causal order (morphe traces.hpp) has two linear orders. Closing the lane
// is the inverse: the two blocks merged across it (merge_cells), then the streets merged at the two junctions
// (merge_at_vertex). Every derived complex is tested for regularity, and a lane that would leave its block, or
// two stores that do not face one block, are refused with the reason. Checks: tests/check_landscape.cpp.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <morphe/formal/traces.hpp>
#include <morphe/topology/complex_rewriting.hpp>
#include <morphe/topology/regular_complex.hpp>
#include <physis/computational_geometry/planar.hpp>

namespace urban::landscape {

namespace mt = morphe::topology;
namespace cg = physis::computational_geometry;
using Point = cg::Point2;

struct Segment {
  std::size_t u = 0, v = 0;
  std::vector<Point> line;  // from junction u to junction v; empty for straight
  std::string kind = "street";
  std::string name;
};

struct Network {
  std::vector<Point> junctions;  // local east, north (m)
  std::vector<Segment> segments;
};

struct Fibre {
  std::string name, unit;
  bool extensive = false;
};

struct Store {
  std::string name, use;
  Point at{};
  std::uint64_t segment = 0;  // the 1-cell it faces
  double s = 0;               // its arc length along that segment, from the tail
  double offset = 0;          // its distance from the segment (m)
};

// A structural event of the record: what was rewritten, where, the rewrite's footprint, and whether it only
// describes the same place more finely (a refinement: a street cut at a new junction is the same street) or changes
// what is there (an intervention: a passage that was not there, two blocks where there was one).
struct Step {
  std::string kind, match;
  std::vector<std::uint64_t> created, deleted;
  morphe::formal::traces::Footprint footprint;
  bool refinement = false;
};

inline std::map<std::string, std::vector<Fibre>> default_schema() {
  return {
      {"junction", {{"walkers", "persons", true}}},
      {"portal", {{"walkers", "persons", true}}},
      {"street", {{"flow", "persons/h", false}}},
      {"footway", {{"flow", "persons/h", false}}},
      {"lane", {{"flow", "persons/h", false}}},
      {"block", {{"temperature", "°C", false}, {"capacity", "J/K", true}}},
  };
}

struct District {
  mt::RegularComplex complex;
  std::map<std::uint64_t, Point> point;              // each 0-cell's position
  std::map<std::uint64_t, std::vector<Point>> line;  // each 1-cell's polyline, tail to head
  std::map<std::uint64_t, std::string> name;         // each 1-cell's street name
  std::vector<Store> stores;
  std::map<std::string, std::vector<Fibre>> schema = default_schema();
  std::map<std::uint64_t, std::vector<double>> state;  // cell id ↦ its values, in its type's fibre order
  std::size_t holes = 0;                               // bounded faces left open
  double hole_area = 0;

  std::size_t at(std::uint64_t id) const {
    const std::size_t i = complex.find(id);
    if (i == mt::no_cell) throw std::invalid_argument("urban::landscape::District: no cell #" + std::to_string(id));
    return i;
  }
  const mt::RegularCell& cell(std::uint64_t id) const { return complex.cells[at(id)]; }
  double length(std::uint64_t edge) const { return cg::polyline_length(line.at(edge)); }
  // The boundary of a 2-cell walked in the direction of its orientation: (edge id, +1 when walked tail → head).
  std::vector<std::pair<std::uint64_t, int>> walk(std::uint64_t face) const {
    const auto& c = cell(face);
    const std::size_t m = c.facets.size();
    const auto start_of = [&](std::size_t j) { const auto e = complex.ends(c.facets[j]); return c.signs[j] > 0 ? e[0] : e[1]; };
    const auto end_of = [&](std::size_t j) { const auto e = complex.ends(c.facets[j]); return c.signs[j] > 0 ? e[1] : e[0]; };
    std::vector<char> used(m, 0);
    std::vector<std::pair<std::uint64_t, int>> out{{complex.cells[c.facets[0]].id, c.signs[0]}};
    used[0] = 1;
    std::size_t here = end_of(0);
    for (std::size_t step = 1; step < m; step++) {
      std::size_t j = m;
      for (std::size_t q = 0; q < m; q++)
        if (!used[q] && start_of(q) == here) { j = q; break; }
      if (j == m) throw std::logic_error("urban::landscape::District::walk: a block's boundary is not a circuit");
      used[j] = 1;
      out.push_back({complex.cells[c.facets[j]].id, c.signs[j]});
      here = end_of(j);
    }
    return out;
  }
  // The ring of a 2-cell, counter-clockwise, from its edges' polylines.
  std::vector<Point> ring(std::uint64_t face) const {
    std::vector<Point> r;
    for (const auto& [e, sign] : walk(face)) {
      auto l = line.at(e);
      if (sign < 0) std::reverse(l.begin(), l.end());
      r.insert(r.end(), l.begin(), l.end() - 1);
    }
    return r;
  }
  double area(std::uint64_t face) const { return cg::ring_area(ring(face)); }
  std::size_t fibre(const std::string& type, const std::string& fibre_name) const {
    const auto it = schema.find(type);
    if (it != schema.end())
      for (std::size_t k = 0; k < it->second.size(); k++)
        if (it->second[k].name == fibre_name) return k;
    throw std::invalid_argument("urban::landscape::District: a " + type + " carries no " + fibre_name);
  }
  double value(std::uint64_t id, const std::string& f) const { return state.at(id).at(fibre(cell(id).type, f)); }
  void set(std::uint64_t id, const std::string& f, double v) {
    auto& s = state[id];
    const std::size_t k = fibre(cell(id).type, f);
    if (s.size() <= k) s.resize(schema.at(cell(id).type).size(), 0.0);
    s[k] = v;
  }
  // Every cell of a type, by id, in the complex's order.
  std::vector<std::uint64_t> of(std::size_t dim) const {
    std::vector<std::uint64_t> out;
    for (const auto& c : complex.cells)
      if (c.dim == dim) out.push_back(c.id);
    return out;
  }
};

// A plane network from one that may not be: of every two segments that cross without a junction (a ramp over a
// lane, ways at two levels), the shorter is dropped, until none cross; the names of those dropped are reported.
// The district then says which streets it left out, rather than inventing a junction where there is none.
inline Network planarised(Network N, std::vector<std::string>* dropped = nullptr) {
  const auto pieces = [&](const Segment& s) { return s.line.empty() ? std::vector<Point>{N.junctions[s.u], N.junctions[s.v]} : s.line; };
  while (true) {
    std::vector<std::vector<Point>> lines;
    std::vector<std::array<double, 4>> box;
    for (const auto& sg : N.segments) {
      lines.push_back(pieces(sg));
      std::array<double, 4> b{1e300, 1e300, -1e300, -1e300};
      for (const auto& q : lines.back()) b = {std::min(b[0], q[0]), std::min(b[1], q[1]), std::max(b[2], q[0]), std::max(b[3], q[1])};
      box.push_back(b);
    }
    std::vector<std::size_t> crossings(N.segments.size(), 0);
    std::size_t worst = N.segments.size();
    for (std::size_t i = 0; i < lines.size(); i++)
      for (std::size_t j = i + 1; j < lines.size(); j++) {
        if (box[i][2] < box[j][0] || box[j][2] < box[i][0] || box[i][3] < box[j][1] || box[j][3] < box[i][1]) continue;
        bool cross = false;
        for (std::size_t a = 0; a + 1 < lines[i].size() && !cross; a++)
          for (std::size_t b = 0; b + 1 < lines[j].size() && !cross; b++) cross = cg::segments_cross(lines[i][a], lines[i][a + 1], lines[j][b], lines[j][b + 1]);
        if (!cross) continue;
        const std::size_t shorter = cg::polyline_length(lines[i]) <= cg::polyline_length(lines[j]) ? i : j;
        crossings[shorter]++;
        if (worst == N.segments.size() || crossings[shorter] > crossings[worst]) worst = shorter;
      }
    if (worst == N.segments.size()) return N;
    if (dropped) dropped->push_back(N.segments[worst].name.empty() ? N.segments[worst].kind : N.segments[worst].name);
    N.segments.erase(N.segments.begin() + static_cast<std::ptrdiff_t>(worst));
  }
}

// The district of a walking network; blocks start at temperature t0 with heat capacity per square metre cap.
inline District district(const Network& N, double t0 = 20, double cap_per_m2 = 2.0e5) {
  std::vector<cg::PlaneEdge> E;
  for (const auto& s : N.segments) E.push_back({s.u, s.v, s.line});
  const auto faces = cg::plane_faces(N.junctions, E);
  // Components over the vertices that some segment touches.
  std::vector<std::size_t> parent(N.junctions.size());
  std::iota(parent.begin(), parent.end(), 0);
  const auto root = [&](std::size_t x) { while (parent[x] != x) x = parent[x] = parent[parent[x]]; return x; };
  std::vector<char> touched(N.junctions.size(), 0);
  for (const auto& s : N.segments) parent[root(s.u)] = root(s.v), touched[s.u] = touched[s.v] = 1;
  std::set<std::size_t> comps;
  std::size_t V = 0;
  for (std::size_t v = 0; v < N.junctions.size(); v++)
    if (touched[v]) comps.insert(root(v)), V++;
  if (static_cast<long long>(V) - static_cast<long long>(N.segments.size()) + static_cast<long long>(faces.size()) != 2 * static_cast<long long>(comps.size()))
    throw std::invalid_argument("urban::landscape::district: the network is not plane (segments cross without a junction)");
  // Euler's relation cannot see two trees crossing; every pair of pieces is tested, boxes first.
  {
    std::vector<std::vector<Point>> lines;
    std::vector<std::array<double, 4>> box;
    for (const auto& sg : N.segments) {
      lines.push_back(sg.line.empty() ? std::vector<Point>{N.junctions[sg.u], N.junctions[sg.v]} : sg.line);
      std::array<double, 4> b{1e300, 1e300, -1e300, -1e300};
      for (const auto& q : lines.back()) b = {std::min(b[0], q[0]), std::min(b[1], q[1]), std::max(b[2], q[0]), std::max(b[3], q[1])};
      box.push_back(b);
    }
    for (std::size_t i = 0; i < lines.size(); i++)
      for (std::size_t j = i + 1; j < lines.size(); j++) {
        if (box[i][2] < box[j][0] || box[j][2] < box[i][0] || box[i][3] < box[j][1] || box[j][3] < box[i][1]) continue;
        for (std::size_t a = 0; a + 1 < lines[i].size(); a++)
          for (std::size_t b = 0; b + 1 < lines[j].size(); b++)
            if (cg::segments_cross(lines[i][a], lines[i][a + 1], lines[j][b], lines[j][b + 1]))
              throw std::invalid_argument("urban::landscape::district: the network is not plane: " + (N.segments[i].name.empty() ? std::string("a segment") : N.segments[i].name) + " crosses " + (N.segments[j].name.empty() ? std::string("another") : N.segments[j].name) + " without a junction");
      }
  }
  District D;
  D.complex.name = "district";
  std::vector<std::size_t> vcell(N.junctions.size(), mt::no_cell);
  for (std::size_t v = 0; v < N.junctions.size(); v++) {
    if (!touched[v]) continue;
    vcell[v] = D.complex.add_vertex("junction");
    const auto id = D.complex.cells[vcell[v]].id;
    D.point[id] = N.junctions[v];
    D.state[id] = {0.0};
  }
  std::vector<std::size_t> ecell(N.segments.size());
  for (std::size_t e = 0; e < N.segments.size(); e++) {
    const auto& s = N.segments[e];
    if (!D.schema.count(s.kind)) D.schema[s.kind] = {{"flow", "persons/h", false}};
    ecell[e] = D.complex.add_edge(vcell[s.u], vcell[s.v], s.kind);
    const auto id = D.complex.cells[ecell[e]].id;
    D.line[id] = s.line.empty() ? std::vector<Point>{N.junctions[s.u], N.junctions[s.v]} : s.line;
    D.name[id] = s.name;
    D.state[id] = {0.0};
  }
  for (const auto& f : faces) {
    if (!(f.area > 0)) continue;
    if (!f.simple) { D.holes++, D.hole_area += f.area; continue; }
    std::vector<std::size_t> facets;
    std::vector<int> signs;
    for (const auto& [e, sgn] : f.edges) facets.push_back(ecell[e]), signs.push_back(sgn);
    const std::size_t c = D.complex.add(2, "block", facets, signs);
    const auto id = D.complex.cells[c].id;
    D.state[id] = {t0, cap_per_m2 * f.area};
  }
  const auto a = mt::admissible(D.complex);
  if (!a.regular) throw std::logic_error("urban::landscape::district: the traced complex is not regular: " + a.reason);
  return D;
}

// Attaches a store to the segment nearest to it.
inline std::size_t attach_store(District& D, Store s) {
  double best = std::numeric_limits<double>::infinity();
  for (const auto& [id, l] : D.line) {
    const auto [at, dist] = cg::nearest_on_polyline(l, s.at);
    if (dist < best) best = dist, s.segment = id, s.s = at, s.offset = dist;
  }
  if (!std::isfinite(best)) throw std::invalid_argument("urban::landscape::attach_store: the district has no segments");
  D.stores.push_back(std::move(s));
  return D.stores.size() - 1;
}

// ---- the lane --------------------------------------------------------------------------------------------------

struct LanePlan {
  std::uint64_t block = 0, edge_a = 0, edge_b = 0;
  double s_a = 0, s_b = 0;
  Point a{}, b{};
  std::string refused;
};

// Where a lane between two stores would run: from the point of each store's segment that faces it, through the
// block both segments bound.
inline LanePlan plan_lane(const District& D, std::size_t store_a, std::size_t store_b) {
  LanePlan P;
  const auto& A = D.stores.at(store_a);
  const auto& B = D.stores.at(store_b);
  if (A.segment == B.segment) { P.refused = "the two stores face one segment: they need no lane"; return P; }
  const auto up = D.complex.cofaces();
  std::set<std::uint64_t> fa, fb;
  for (std::size_t c : up[D.at(A.segment)]) fa.insert(D.complex.cells[c].id);
  for (std::size_t c : up[D.at(B.segment)]) fb.insert(D.complex.cells[c].id);
  std::vector<std::uint64_t> common;
  for (auto f : fa)
    if (fb.count(f)) common.push_back(f);
  if (common.empty()) { P.refused = "the two stores do not face one block"; return P; }
  P.edge_a = A.segment, P.edge_b = B.segment;
  const double la = D.length(A.segment), lb = D.length(B.segment);
  P.s_a = std::clamp(A.s, 0.08 * la, 0.92 * la), P.s_b = std::clamp(B.s, 0.08 * lb, 0.92 * lb);
  P.a = cg::point_at(D.line.at(A.segment), P.s_a), P.b = cg::point_at(D.line.at(B.segment), P.s_b);
  // The block whose interior holds the chord.
  for (auto f : common) {
    const auto r = D.ring(f);
    const Point mid{(P.a[0] + P.b[0]) / 2, (P.a[1] + P.b[1]) / 2};
    // The chord starts and ends on the ring: the ring's edges through its ends are not crossings, whatever the
    // rounding of where on them it starts.
    const auto through = [](Point p, Point c, Point d) { return cg::nearest_on_polyline(std::vector<Point>{c, d}, p).second < 1e-6; };
    bool crosses = false;
    for (std::size_t k = 0; k < r.size() && !crosses; k++) {
      const Point c = r[k], d = r[(k + 1) % r.size()];
      if (through(P.a, c, d) || through(P.b, c, d)) continue;
      crosses = cg::segments_cross(P.a, P.b, c, d);
    }
    if (!crosses && cg::inside_ring(r, mid)) { P.block = f; return P; }
  }
  P.refused = "a straight lane between the two stores would leave their block";
  return P;
}

namespace detail {

// Subdivides a segment at arc length s: geometry, stores and state carried; returns the step.
inline Step subdivide(District& D, std::uint64_t edge, double s) {
  const auto q = mt::subdivide_edge(D.complex, D.at(edge), "portal");
  const auto r = mt::rewrite(D.complex, q);
  if (!r.applied) throw std::invalid_argument("urban::landscape: the street cannot be subdivided: " + r.reason);
  const std::uint64_t m = r.created[0], e1 = r.created[1], e2 = r.created[2];
  const auto [first, second] = cg::split_polyline(D.line.at(edge), s);
  D.complex = r.result;
  D.point[m] = first.back();
  D.line[e1] = first, D.line[e2] = second;
  D.name[e1] = D.name[e2] = D.name[edge];
  D.state[m] = {0.0};                                       // a new junction: no walkers yet
  D.state[e1] = D.state[e2] = D.state.at(edge);             // the flow along the segment, on both pieces
  for (auto& st : D.stores)
    if (st.segment == edge) {
      if (st.s <= s) st.segment = e1;
      else st.segment = e2, st.s -= s;
    }
  D.line.erase(edge), D.name.erase(edge), D.state.erase(edge);
  return {"subdivide", D.name[e1].empty() ? "a street" : D.name[e1], r.created, r.deleted, r.footprint, true};
}

}  // namespace detail

// Builds the lane: three rewrites, applied to D in order, or none (D unchanged) when one is refused.
inline std::vector<Step> build_lane(District& D, const LanePlan& P) {
  if (!P.refused.empty()) throw std::invalid_argument("urban::landscape::build_lane: " + P.refused);
  for (const auto id : {P.edge_a, P.edge_b, P.block})
    if (D.complex.find(id) == mt::no_cell) throw std::invalid_argument("urban::landscape::build_lane: the plan was made for another configuration of the cells (no cell #" + std::to_string(id) + "): plan the lane again");
  District before = D;
  try {
    std::vector<Step> steps;
    steps.push_back(detail::subdivide(D, P.edge_a, P.s_a));
    const std::uint64_t ma = steps.back().created[0];
    steps.push_back(detail::subdivide(D, P.edge_b, P.s_b));
    const std::uint64_t mb = steps.back().created[0];
    const double t = D.value(P.block, "temperature"), C = D.value(P.block, "capacity"), A = D.area(P.block);
    const auto q = mt::split_face(D.complex, D.at(P.block), D.at(ma), D.at(mb), "lane");
    const auto r = mt::rewrite(D.complex, q);
    if (!r.applied) throw std::invalid_argument("the block cannot be split: " + r.reason);
    D.complex = r.result;
    const std::uint64_t lane = r.created[0], f1 = r.created[1], f2 = r.created[2];
    const auto ends = D.complex.ends(D.at(lane));
    D.line[lane] = {D.point.at(D.complex.cells[ends[0]].id), D.point.at(D.complex.cells[ends[1]].id)};
    D.name[lane] = "a lane";
    D.state[lane] = {0.0};
    // Temperature copied; heat capacity split by area, so that heat = C T is kept.
    const double a1 = D.area(f1), a2 = D.area(f2);
    D.state[f1] = {t, C * a1 / A}, D.state[f2] = {t, C * a2 / A};
    D.state.erase(P.block);
    steps.push_back({"split", "a block, by a lane", r.created, r.deleted, r.footprint});
    return steps;
  } catch (...) {
    D = std::move(before);
    throw;
  }
}

// Closes a lane: the two blocks merged across it, then each street merged at the lane's junction if that
// junction is left on exactly two streets.
inline std::vector<Step> close_lane(District& D, std::uint64_t lane) {
  District before = D;
  try {
    std::vector<Step> steps;
    const auto ends = D.complex.ends(D.at(lane));
    const std::uint64_t pa = D.complex.cells[ends[0]].id, pb = D.complex.cells[ends[1]].id;
    const auto up = D.complex.cofaces()[D.at(lane)];
    if (up.size() != 2) throw std::invalid_argument("the lane does not run between two blocks");
    const std::uint64_t f1 = D.complex.cells[up[0]].id, f2 = D.complex.cells[up[1]].id;
    const double C1 = D.value(f1, "capacity"), C2 = D.value(f2, "capacity"), T1 = D.value(f1, "temperature"), T2 = D.value(f2, "temperature");
    const auto r = mt::rewrite(D.complex, mt::merge_cells(D.complex, D.at(lane)));
    if (!r.applied) throw std::invalid_argument("the blocks cannot be merged: " + r.reason);
    D.complex = r.result;
    const std::uint64_t merged = r.created[0];
    D.state[merged] = {(C1 * T1 + C2 * T2) / (C1 + C2), C1 + C2};  // heat and capacity summed
    for (auto id : {f1, f2, lane}) D.state.erase(id);
    D.line.erase(lane), D.name.erase(lane);
    steps.push_back({"merge", "two blocks, across a lane", r.created, r.deleted, r.footprint});
    for (const std::uint64_t m : {pa, pb}) {
      const std::size_t mi = D.at(m);
      const auto upm = D.complex.cofaces()[mi];
      if (upm.size() != 2 || D.complex.cells[mi].type != "portal") continue;
      const std::uint64_t e1 = D.complex.cells[upm[0]].id, e2 = D.complex.cells[upm[1]].id;
      const auto other = [&](std::uint64_t e) { const auto x = D.complex.ends(D.at(e)); return D.complex.cells[x[0]].id == m ? D.complex.cells[x[1]].id : D.complex.cells[x[0]].id; };
      const std::uint64_t a = other(e1), b = other(e2);
      auto l1 = D.line.at(e1), l2 = D.line.at(e2);
      if (D.complex.cells[D.complex.ends(D.at(e1))[0]].id != a) std::reverse(l1.begin(), l1.end());  // a → m
      if (D.complex.cells[D.complex.ends(D.at(e2))[0]].id != m) std::reverse(l2.begin(), l2.end());  // m → b
      const double walkers = D.value(m, "walkers"), flow = 0.5 * (D.value(e1, "flow") + D.value(e2, "flow"));
      const auto q = mt::merge_at_vertex(D.complex, mi);
      const auto rr = mt::rewrite(D.complex, q);
      if (!rr.applied) throw std::invalid_argument("the street cannot be merged at the lane's junction: " + rr.reason);
      D.complex = rr.result;
      const std::uint64_t e = rr.created[0];
      std::vector<Point> joined = l1;
      joined.insert(joined.end(), l2.begin() + 1, l2.end());
      {
        // The junction's point, inserted on a straight piece when the street was cut, is no vertex of the street:
        // dropped, the street has its polyline again.
        const std::size_t k = l1.size() - 1;
        const Point &p = joined[k - 1], &q = joined[k], &w = joined[k + 1];
        const double cross = (q[0] - p[0]) * (w[1] - p[1]) - (q[1] - p[1]) * (w[0] - p[0]);
        const double scale = std::hypot(q[0] - p[0], q[1] - p[1]) * std::hypot(w[0] - q[0], w[1] - q[1]);
        if (std::fabs(cross) <= 1e-9 * scale || scale == 0) joined.erase(joined.begin() + static_cast<std::ptrdiff_t>(k));
      }
      D.line[e] = joined;
      D.name[e] = D.name[e1];
      D.state[e] = {flow};
      // The junction's walkers go to the two ends, half each.
      D.set(a, "walkers", D.value(a, "walkers") + walkers / 2);
      D.set(b, "walkers", D.value(b, "walkers") + walkers / 2);
      for (auto& st : D.stores)
        if (st.segment == e1 || st.segment == e2) {
          const Point p = st.at;
          const auto [s, dist] = cg::nearest_on_polyline(joined, p);
          st.segment = e, st.s = s, st.offset = dist;
        }
      for (auto id : {m, e1, e2}) D.state.erase(id), D.line.erase(id), D.name.erase(id), D.point.erase(id);
      steps.push_back({"merge", D.name[e].empty() ? "a street" : D.name[e], rr.created, rr.deleted, rr.footprint, true});
    }
    return steps;
  } catch (...) {
    D = std::move(before);
    throw;
  }
}

}  // namespace urban::landscape
