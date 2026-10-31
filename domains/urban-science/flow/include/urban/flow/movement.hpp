// Movement among the places of a district as a Markov jump process (theory/street-oceanography.md,
// Definition 6.2 and Theorem 6.3): its transit times, and its circulation, which only the non-exact part of
// its rates can drive.
//
// The district is a grid of streets bounding blocks: a 2-complex (morphe's CW complex) whose vertices are
// intersections, whose edges are street segments and whose faces are blocks; a block left out is a hole (a
// park, a station, a river). Movement along an edge e = (u → v) has rates k_uv = c e^{φ_e/2}, k_vu = c e^{−φ_e/2},
// so that φ = ln(k_uv / k_vu) is a 1-cochain, built here as
//   φ = dV + Σ_blocks Γ_b ∂-orientation/4 + Σ_holes Γ_h ∂-orientation/|∂h|,
// an attraction V on the intersections (movement toward higher V favoured), a circulation Γ_b round each
// block, and one round each hole.
//
// Theorem 6.3: the process satisfies detailed balance exactly when φ is exact (Kolmogorov 1936); its cycle
// affinities are the periods of φ's non-exact part, whatever V; its entropy production σ = Σ_C J(C) A(C) ≥ 0
// (Schnakenberg 1976). physis computes σ and the affinities (statistical_mechanics) and the Hodge
// decomposition φ = dα + δγ + h (discretisation): α recovers V (up to a constant), δγ the blocks' circulation,
// h the holes'.
//
// Transit (Definition 6.2): exits at chosen intersections lead off the district at a rate; the time from an
// entry to leaving is phase-type (physis first_passage).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <morphe/topology/cw_complex.hpp>
#include <physis/discretisation/hodge.hpp>
#include <physis/statistical_mechanics/first_passage.hpp>
#include <physis/statistical_mechanics/master_equation.hpp>

namespace urban::flow {

namespace sm = physis::statistical_mechanics;
namespace dz = physis::discretisation;

struct District {
  std::size_t nx = 4, ny = 3;                     // blocks across and up
  std::vector<std::size_t> holes;                 // blocks (j * nx + i) left out of the complex
  std::vector<double> attraction;                 // V per intersection ((nx+1)(ny+1)); empty: 0
  std::vector<double> block_circulation;          // Γ_b per block (nx·ny; holes ignored); empty: 0
  std::vector<double> hole_circulation;           // Γ_h per hole, in the order of `holes`
  double rate = 1.0 / 60;                         // c, 1/s (a minute to walk a segment)
  std::vector<std::size_t> exits;                 // intersections that lead off the district
  double exit_rate = 1.0 / 60;                    // 1/s
  std::size_t entry = 0;                          // where transit starts
};

inline std::size_t vertex(const District& d, std::size_t i, std::size_t j) { return j * (d.nx + 1) + i; }
inline std::size_t h_edge(const District& d, std::size_t i, std::size_t j) { return j * d.nx + i; }                     // (i,j)→(i+1,j)
inline std::size_t v_edge(const District& d, std::size_t i, std::size_t j) { return (d.ny + 1) * d.nx + j * (d.nx + 1) + i; }  // (i,j)→(i,j+1)

// The block's boundary, counterclockwise from its lower-left corner.
inline morphe::topology::Path block_boundary(const District& d, std::size_t i, std::size_t j) {
  return {{h_edge(d, i, j), +1}, {v_edge(d, i + 1, j), +1}, {h_edge(d, i, j + 1), -1}, {v_edge(d, i, j), -1}};
}

inline morphe::topology::CWComplex district_complex(const District& d) {
  morphe::topology::CWComplex X;
  X.name = "district";
  X.vertices = (d.nx + 1) * (d.ny + 1);
  for (std::size_t j = 0; j <= d.ny; j++)
    for (std::size_t i = 0; i < d.nx; i++) X.edges.push_back({vertex(d, i, j), vertex(d, i + 1, j)});
  for (std::size_t j = 0; j < d.ny; j++)
    for (std::size_t i = 0; i <= d.nx; i++) X.edges.push_back({vertex(d, i, j), vertex(d, i, j + 1)});
  for (std::size_t j = 0; j < d.ny; j++)
    for (std::size_t i = 0; i < d.nx; i++) {
      const std::size_t b = j * d.nx + i;
      if (std::find(d.holes.begin(), d.holes.end(), b) != d.holes.end()) continue;
      X.faces.push_back({vertex(d, i, j), block_boundary(d, i, j)});
    }
  return X;
}

// φ on the edges.
inline dz::Cochain log_ratios(const District& d) {
  const auto X = district_complex(d);
  dz::Cochain phi(X.edges.size(), 0);
  if (!d.attraction.empty()) {
    if (d.attraction.size() != X.vertices) throw std::invalid_argument("urban::flow::log_ratios: an attraction for each intersection");
    for (std::size_t e = 0; e < X.edges.size(); e++) phi[e] += d.attraction[X.edges[e][1]] - d.attraction[X.edges[e][0]];
  }
  for (std::size_t j = 0; j < d.ny; j++)
    for (std::size_t i = 0; i < d.nx; i++) {
      const std::size_t b = j * d.nx + i;
      const auto hole = std::find(d.holes.begin(), d.holes.end(), b);
      double gamma = 0, share = 4;
      if (hole == d.holes.end()) gamma = b < d.block_circulation.size() ? d.block_circulation[b] : 0;
      else {
        const std::size_t h = static_cast<std::size_t>(hole - d.holes.begin());
        gamma = h < d.hole_circulation.size() ? d.hole_circulation[h] : 0;
      }
      for (const auto& s : block_boundary(d, i, j)) phi[s.edge] += s.sign * gamma / share;
    }
  return phi;
}

// The rates among the intersections (closed district), and with an absorbing "away" state for the exits.
inline calx::linalg::Matrix movement_rates(const District& d, bool with_exits) {
  const auto X = district_complex(d);
  const auto phi = log_ratios(d);
  const std::size_t n = X.vertices + (with_exits ? 1 : 0);
  calx::linalg::Matrix k(n);
  for (std::size_t e = 0; e < X.edges.size(); e++) {
    const auto [u, v] = X.edges[e];
    k(u, v) += d.rate * std::exp(phi[e] / 2);
    k(v, u) += d.rate * std::exp(-phi[e] / 2);
  }
  if (with_exits)
    for (auto x : d.exits) k(x, X.vertices) += d.exit_rate;
  return k;
}

struct Circulation {
  morphe::topology::CWComplex complex;
  dz::Cochain phi;
  dz::HodgeDecomposition hodge;
  sm::Irreversibility irreversibility;
  std::vector<double> current;           // J on each edge (tail → head), the stationary current
  double attraction_free_affinity_gap = 0;   // max |A(C) − A₀(C)| between the district and itself without V
};

inline Circulation circulation(const District& d) {
  Circulation c;
  c.complex = district_complex(d);
  c.phi = log_ratios(d);
  c.hodge = dz::hodge_decomposition(c.complex, c.phi);
  const auto k = movement_rates(d, false);
  c.irreversibility = sm::irreversibility(k);
  for (const auto& e : c.complex.edges) c.current.push_back(c.irreversibility.current(e[0], e[1]));
  District flat = d;
  flat.attraction.clear();
  const auto i0 = sm::irreversibility(movement_rates(flat, false));
  for (std::size_t q = 0; q < std::min(i0.cycles.size(), c.irreversibility.cycles.size()); q++)
    c.attraction_free_affinity_gap = std::max(c.attraction_free_affinity_gap, std::fabs(i0.cycles[q].affinity - c.irreversibility.cycles[q].affinity));
  return c;
}

struct Transit {
  std::vector<double> time, density, survival;   // s
  double mean = 0, median = 0, p90 = 0;           // s
};

inline Transit transit(const District& d, double t_max, std::size_t points = 120) {
  if (d.exits.empty()) throw std::invalid_argument("urban::flow::transit: a district with exits");
  const auto k = movement_rates(d, true);
  const std::size_t n = k.n - 1;
  std::vector<std::size_t> transient(n);
  for (std::size_t i = 0; i < n; i++) transient[i] = i;
  std::vector<double> alpha(n, 0);
  alpha.at(d.entry) = 1;
  const auto ph = sm::phase_type(k, transient, alpha);
  Transit t;
  for (std::size_t q = 0; q <= points; q++) {
    const double tau = t_max * static_cast<double>(q) / static_cast<double>(points);
    t.time.push_back(tau);
    t.density.push_back(ph.density(tau));
    t.survival.push_back(ph.survival(tau));
  }
  t.mean = ph.moment(1);
  t.median = ph.quantile(0.5);
  t.p90 = ph.quantile(0.9);
  return t;
}

}  // namespace urban::flow
