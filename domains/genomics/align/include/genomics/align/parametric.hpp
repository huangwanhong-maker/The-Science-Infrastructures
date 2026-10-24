// Parametric alignment of two DNA or protein sequences: every optimal global alignment under the
// scores match 1, mismatch −α, gap −β (per column, linear gaps), for all (α, β) at once, by morphe's
// polytope propagation (morphe/formal/polytope.hpp).
//
// An alignment's score is matches − α mismatches − β gaps. In global alignment
// 2 matches + 2 mismatches + gaps = m + n, so with x the mismatches and g the gap columns
//
//   score = (m + n)/2 − (1 + α) x − (1/2 + β) g,
//
// a function of the point (x, g) alone, and the optimal alignments at (α, β) are those whose point
// minimises ⟨c, ·⟩ for c = (1 + α, 1/2 + β) over the alignment polygon, the convex hull of all
// the points (Pachter and Sturmfels 2004, PNAS 101, 16138–16143; Waterman, Eggert and Lander 1992,
// PNAS 89, 6090–6093; Gusfield, Balasubramanian and Naor 1994, Algorithmica 12, 312–326). A vertex
// v is optimal exactly when −c lies in its outer normal cone N(v), so its region in the (α, β)
// plane is the cone with apex (−1, −1/2) spanned by the negated generators of N(v): the plane is
// cut into cones from one point, one for each vertex, and the parameters on a ray between two
// cones are those at which both vertices (and every alignment on the edge between them) are
// optimal. The quadrant α, β ≥ 0 that biology uses is the part of this fan beyond (0, 0).
//
// Each vertex carries a witness, a global alignment whose point it is (the tropical optimum at an
// interior direction of its cone, which is unique there), as a CIGAR. Letters are equal when
// same(a, b, limit) of types.hpp holds: A, C, G, T for DNA (limit 4, N matches nothing); the
// twenty amino acids for protein (limit 20, B, Z, X and * match nothing).
//
// Complexity: O(mn V log V) for polygons of at most V = O((m + n)^{2/3}) vertices, and O(V mn)
// for the witnesses; memory O(n V) for the propagation and O(mn) for one traceback. Sequences of a
// few hundred letters take a fraction of a second.

#pragma once

#include <genomics/align/types.hpp>

#include <morphe/formal/alignment.hpp>
#include <morphe/formal/polytope.hpp>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace genomics::align {

struct ParametricVertex {
  std::int64_t mismatches = 0, gaps = 0, matches = 0;
  Alignment alignment;  // a witness: a global alignment with these counts
  // The region in the (α, β) plane where the vertex is optimal: the cone with apex (−1, −1/2)
  // from `first` counterclockwise to `second` (directions in the (α, β) plane), or all of the
  // plane (a polygon of one point), or the half-plane {(α, β) : ⟨(1 + α, 1/2 + β), first⟩ ≥ 0}
  // (a segment).
  morphe::formal::NormalCone::Kind kind = morphe::formal::NormalCone::Kind::plane;
  morphe::formal::LatticePoint first, second;
  // Whether the vertex is optimal at (α, β).
  bool optimal_at(double alpha, double beta) const {
    const double a = 1 + alpha, b = 0.5 + beta;
    const auto cross = [](double ux, double uy, double vx, double vy) { return ux * vy - uy * vx; };
    switch (kind) {
      case morphe::formal::NormalCone::Kind::plane: return true;
      case morphe::formal::NormalCone::Kind::half_plane: return a * static_cast<double>(first.x) + b * static_cast<double>(first.y) >= -1e-12;
      case morphe::formal::NormalCone::Kind::sector:
        return cross(static_cast<double>(first.x), static_cast<double>(first.y), a, b) >= -1e-12 && cross(a, b, static_cast<double>(second.x), static_cast<double>(second.y)) >= -1e-12;
    }
    return false;
  }
};

struct Parametric {
  std::size_t query_length = 0, target_length = 0;
  std::vector<ParametricVertex> vertices;  // counterclockwise around the polygon of (mismatches, gaps)
  // The score of a vertex's alignments at (α, β).
  double score(std::size_t vertex, double alpha, double beta) const {
    const auto& v = vertices[vertex];
    return static_cast<double>(v.matches) - alpha * static_cast<double>(v.mismatches) - beta * static_cast<double>(v.gaps);
  }
  // A vertex optimal at (α, β): the highest score, the first in the vertex order on ties.
  std::size_t optimal(double alpha, double beta) const {
    std::size_t best = 0;
    for (std::size_t k = 1; k < vertices.size(); k++)
      if (score(k, alpha, beta) > score(best, alpha, beta)) best = k;
    return best;
  }
  const Alignment& optimal_alignment(double alpha, double beta) const { return vertices[optimal(alpha, beta)].alignment; }
};

// The parametric global alignment of the query and the target (DNA codes when `dna`, protein codes
// otherwise).
inline Parametric parametric(const Codes& query, const Codes& target, bool dna = true) {
  const Code limit = dna ? 4 : 20;
  const auto equal = [limit](Code a, Code b) { return same(a, b, limit); };  // types.hpp
  const auto poly = morphe::formal::alignment_polygon(query, target, equal, morphe::formal::Mode::global);
  Parametric out;
  out.query_length = query.size();
  out.target_length = target.size();
  const std::int64_t total = static_cast<std::int64_t>(query.size() + target.size());
  for (std::size_t k = 0; k < poly.polygon.size(); k++) {
    const auto p = poly.polygon.vertices()[k];
    const auto& cone = poly.fan[k];
    ParametricVertex v;
    v.mismatches = p.x;
    v.gaps = p.y;
    v.matches = (total - p.y) / 2 - p.x;
    v.kind = cone.kind;
    // The minimising direction c ranges over −N(v): negate the generators (a rotation by π keeps
    // the counterclockwise order of a sector).
    v.first = {-cone.first.x, -cone.first.y};
    v.second = {-cone.second.x, -cone.second.y};
    const auto& w = poly.witness[k];
    v.alignment.query_begin = v.alignment.target_begin = 0;
    v.alignment.query_end = query.size();
    v.alignment.target_end = target.size();
    for (const auto& r : w.runs()) push(v.alignment.cigar, r.op, static_cast<std::uint32_t>(r.length));
    v.alignment.score = v.matches;  // the score at (α, β) = (0, 0)
    out.vertices.push_back(std::move(v));
  }
  return out;
}

}  // namespace genomics::align
