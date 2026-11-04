// The measures of a part: its surface's topology, its area, the solid it bounds (volume, centre, inertia,
// principal axes and the extents along them), a compound's measures composed from its elements', and the
// parts that touch (ADR 0066; theory/anatomy.md §4).
//
// Nothing here is new mathematics; it is the foundations' applied to a part of the body:
// - the topology of a mesh is morphe's (topology/triangulated.hpp): its components, each classified as a
//   surface or not, and whether the triangles' own orientation is consistent. The mesh is welded first
//   (physis's welded, at 10⁻⁴ mm), since a seam whose points are repeated is not a surface;
// - its area and sections are physis's (computational_geometry/surfaces.hpp);
// - the solid it bounds is physis's solid_polyhedron (mechanics/rigid_body.hpp), by the divergence theorem,
//   which needs every component closed and consistently oriented. A consistent mesh oriented inward gives a
//   negative volume, and is read reversed (`inward`); a mesh that is not closed has no volume, and its centre
//   is the centre of its bounds;
// - the principal axes are physis's principal_axes, the eigenvectors of the inertia of the solid of unit
//   density, the first of least moment: the axis along which the part is longest. The extents are the
//   lengths of the part's projection on each axis, its length, breadth and thickness;
// - a compound's solid is the sum of its elements' (physis's combine, the parallel axis theorem), its copies
//   left out. Overlapping elements are counted twice: BodyParts3D's elements partition a compound only where
//   its tree does;
// - parts touch when some of their points come within a distance (physis's close_pairs): a measure of
//   contiguity, not of shared boundary.
// Volumes are in mm³ and areas in mm², as the meshes' millimetres give them; the database's published
// volumes are in cm³. Checks: tests/check_measures.cpp.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <morphe/topology/triangulated.hpp>
#include <physis/computational_geometry/surfaces.hpp>
#include <physis/mechanics/rigid_body.hpp>

namespace anatomy::measures {

namespace cg = physis::computational_geometry;
namespace pm = physis::mechanics;
using Vec3 = std::array<double, 3>;
using Mat3 = std::array<std::array<double, 3>, 3>;

constexpr double weld_tolerance = 1e-4;  // mm

struct Topology {
  std::size_t components = 0;  // connected parts of the mesh
  std::size_t surfaces = 0;    // of which surfaces (manifolds, possibly with boundary)
  std::size_t closed = 0;      // of which closed
  bool consistent = false;     // every component a closed surface with its triangles consistently oriented
  long long euler = 0;         // V − E + F of the whole mesh
  std::size_t genus = 0;       // summed over the closed orientable components
  std::size_t boundary_circles = 0;
  std::size_t welded = 0;      // points made one by welding
  std::string reason;          // why the first component that is not a surface is not
};

struct Measures {
  Topology topology;
  double area = 0;                                           // mm²
  double volume = std::numeric_limits<double>::quiet_NaN();  // mm³, when the mesh bounds a solid
  bool inward = false;                                       // read reversed
  Vec3 centre{};                                             // of the solid, or of the bounds
  Mat3 inertia{};                                            // of the solid at unit density, about its centre (mm⁵)
  std::array<Vec3, 3> axes{};                                // principal axes, the first of least moment
  Vec3 extents{};                                            // along the axes, mm
  cg::Box3 bounds;
  bool solid() const { return std::isfinite(volume); }
};

inline Topology topology_of(const cg::TriangleSurface& W, std::size_t welded) {
  Topology T;
  T.welded = welded;
  const auto parts = morphe::topology::classify_triangles(W.points.size(), {W.triangles.begin(), W.triangles.end()});
  T.components = parts.size();
  bool all = !parts.empty();
  for (const auto& p : parts) {
    T.euler += static_cast<long long>(p.vertices) - static_cast<long long>(p.edges) + static_cast<long long>(p.triangles);
    if (!p.type.surface) {
      if (T.reason.empty()) T.reason = p.type.reason;
      all = false;
      continue;
    }
    T.surfaces++;
    T.boundary_circles += p.type.boundary_components;
    if (p.type.closed) {
      T.closed++;
      if (p.type.orientable) T.genus += p.type.genus;
    }
    all = all && p.type.closed && p.consistent;
  }
  T.consistent = all;
  return T;
}

// Extents of the points along three unit axes.
inline Vec3 extents_along(const std::vector<cg::Point3>& P, const std::array<Vec3, 3>& axes) {
  Vec3 out{};
  for (int a = 0; a < 3; a++) {
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (const auto& p : P) {
      const double s = p[0] * axes[a][0] + p[1] * axes[a][1] + p[2] * axes[a][2];
      lo = std::min(lo, s), hi = std::max(hi, s);
    }
    out[a] = P.empty() ? 0 : hi - lo;
  }
  return out;
}

inline void set_axes(Measures& M) {
  const auto pa = pm::principal_axes(M.inertia);
  const auto R = pa.axes.matrix();
  for (int a = 0; a < 3; a++) M.axes[a] = {R[0][a], R[1][a], R[2][a]};
}

// The measures of one mesh.
inline Measures measure(const cg::TriangleSurface& S) {
  Measures M;
  std::vector<std::size_t> map;
  cg::TriangleSurface W = cg::welded(S, weld_tolerance, &map);
  M.topology = topology_of(W, S.points.size() - W.points.size());
  M.area = cg::surface_area(W);
  M.bounds = cg::bounds(W.points);
  M.centre = {(M.bounds.lower[0] + M.bounds.upper[0]) / 2, (M.bounds.lower[1] + M.bounds.upper[1]) / 2, (M.bounds.lower[2] + M.bounds.upper[2]) / 2};
  M.axes = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
  if (M.topology.consistent) {
    std::vector<pm::Vec3> P(W.points.begin(), W.points.end());
    std::vector<std::array<std::size_t, 3>> T(W.triangles.begin(), W.triangles.end());
    // physis refuses a surface enclosing no positive volume: a consistent one is then oriented inward.
    auto solid = [&]() -> bool {
      try {
        const auto mp = pm::solid_polyhedron(P, T, 1);
        M.volume = mp.mass;
        M.centre = mp.centre;
        M.inertia = mp.inertia;
        set_axes(M);
        return true;
      } catch (const std::invalid_argument&) {
        return false;
      }
    };
    if (!solid()) {
      for (auto& t : T) std::swap(t[1], t[2]);
      M.inward = solid();
    }
  }
  M.extents = extents_along(W.points, M.axes);
  return M;
}

// The measures of a compound from its elements' (copies left out by the caller): volumes add, the solid is
// combined, the bounds joined, the extents taken over all the elements' points.
inline Measures combine(const std::vector<const Measures*>& parts, const std::vector<const cg::TriangleSurface*>& surfaces) {
  Measures M;
  std::vector<pm::MassProperties> solids;
  bool all_solid = !parts.empty();
  for (const auto* p : parts) {
    M.area += p->area;
    for (int k = 0; k < 3; k++) M.bounds.lower[k] = std::min(M.bounds.lower[k], p->bounds.lower[k]), M.bounds.upper[k] = std::max(M.bounds.upper[k], p->bounds.upper[k]);
    M.topology.components += p->topology.components;
    M.topology.surfaces += p->topology.surfaces;
    M.topology.closed += p->topology.closed;
    M.topology.genus += p->topology.genus;
    M.topology.boundary_circles += p->topology.boundary_circles;
    M.topology.euler += p->topology.euler;
    if (p->solid()) solids.push_back({p->volume, p->centre, p->inertia});
    else all_solid = false;
  }
  M.topology.consistent = all_solid;
  M.centre = {(M.bounds.lower[0] + M.bounds.upper[0]) / 2, (M.bounds.lower[1] + M.bounds.upper[1]) / 2, (M.bounds.lower[2] + M.bounds.upper[2]) / 2};
  M.axes = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
  if (!solids.empty()) {
    const auto B = pm::combine(solids);
    if (all_solid) M.volume = B.mass;
    M.centre = B.centre;
    M.inertia = B.inertia;
    set_axes(M);
  }
  std::vector<cg::Point3> P;
  for (const auto* s : surfaces) P.insert(P.end(), s->points.begin(), s->points.end());
  M.extents = extents_along(P, M.axes);
  return M;
}

// Pairs of parts some of whose points lie within `distance` mm, with their least distance.
inline std::vector<cg::ClosePair> touching(const std::vector<const cg::TriangleSurface*>& surfaces, double distance) {
  std::vector<cg::Point3> P;
  std::vector<std::size_t> label;
  for (std::size_t k = 0; k < surfaces.size(); k++) {
    if (!surfaces[k]) continue;
    for (const auto& p : surfaces[k]->points) P.push_back(p), label.push_back(k);
  }
  return cg::close_pairs(P, label, distance);
}

}  // namespace anatomy::measures
