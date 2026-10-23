// The Complex Atlas's engine: fields and state sums on complexes, by physis's gauge theory
// and quantum gravity (README §4, the protocol).
//
// - gauge: a lattice gauge field of ℤ_N, U(1), SU(2) or SU(N) on the complex's edges, its
//   2-cells the plaquettes, run by heat bath (with overrelaxation) or Metropolis from a seed;
//   the answer has the final configuration (class angles of links and holonomies), the
//   history of the mean plaquette with its error, Wilson loops along given closed paths, the
//   topological charge of U(1) on a closed orientable surface, and the exact mean plaquette
//   where it is known (Migdal–Witten on closed orientable surfaces, the area law on discs).
// - turaev_viro: the Turaev–Viro invariant at level r of a closed triangulated 3-manifold,
//   given by its simplices or by a 4-coloured graph.

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <physis/gauge/run.hpp>
#include <physis/quantum_gravity/turaev_viro.hpp>
#include <physis/symmetry/quantum.hpp>

#include "complexes.hpp"

namespace atlas {

inline Json estimate_json(const physis::field_theory::Estimate& e) {
  return Json::object({{"mean", e.mean}, {"error", e.error}, {"tau", e.tau}, {"n", static_cast<double>(e.n)}});
}

inline Json gauge(const Model& m, const Json& request) {
  const auto& X = m.X;
  if (X.edges.empty() || X.faces.empty()) throw std::invalid_argument("a gauge field needs edges (links) and 2-cells (plaquettes)");
  const Json& g = request["group"];
  physis::gauge::GroupSpec spec;
  spec.kind = g["kind"].text("U1");
  spec.n = static_cast<std::uint32_t>(integer(g, "n", spec.kind == "SU" ? 3 : spec.kind == "Z" ? 2 : 0, 0, 64));
  physis::gauge::RunOptions o;
  o.beta = number(request, "beta", 2.0);
  if (!(o.beta >= 0 && o.beta <= 1000)) throw std::invalid_argument("β must lie in 0..1000");
  o.sweeps = static_cast<std::size_t>(integer(request, "sweeps", 200, 1, 100000));
  o.thermalise = static_cast<std::size_t>(integer(request, "thermalise", 100, 0, 100000));
  o.seed = static_cast<std::uint64_t>(integer(request, "seed", 1, 0, 1ll << 52));
  o.update = request["update"].text("heatbath");
  if (o.update != "heatbath" && o.update != "metropolis") throw std::invalid_argument("the update is \"heatbath\" or \"metropolis\"");
  o.overrelax = static_cast<std::size_t>(integer(request, "overrelax", 1, 0, 10));
  o.hot = request["hot"].is_bool() && request["hot"].as_bool();
  // The work allowed in one request: link updates, weighted by the cost of the group.
  const double weight = spec.kind == "SU" ? static_cast<double>(spec.n * spec.n * spec.n) : spec.kind == "SU2" ? 4.0 : 1.0;
  if (static_cast<double>(X.edges.size()) * static_cast<double>(o.sweeps + o.thermalise) * weight > 2e8)
    throw std::invalid_argument("too much work for one run: fewer sweeps, a smaller complex or a smaller group");
  if (const Json* loops = request.find("loops"); loops && loops->is_array())
    for (const auto& l : loops->as_array()) {
      mt::Path p;
      for (const auto& x : l.as_array()) {
        const double v = x.as_number();
        const std::size_t e = static_cast<std::size_t>(std::fabs(v)) - 1;
        if (v == 0 || e >= X.edges.size()) throw std::invalid_argument("a loop names an edge that is not there");
        p.push_back({e, v > 0 ? 1 : -1});
      }
      if (!p.empty()) {
        const std::size_t start = p[0].sign > 0 ? X.edges[p[0].edge][0] : X.edges[p[0].edge][1];
        if (X.walk(start, p) != start) throw std::invalid_argument("a Wilson loop must be a closed path");
      }
      o.loops.push_back(p);
    }
  const auto R = physis::gauge::run(X, spec, o);
  Json links = Json::array();
  for (const auto& l : R.links) links.push(Json(l));
  Json loops = Json::array();
  for (std::size_t i = 0; i < R.loops.size(); i++) {
    Json e = estimate_json(R.loops[i]);
    e.set("loop", word_json(o.loops[i]));
    loops.push(e);
  }
  Json out = Json::object({{"group", R.group}, {"beta", o.beta}, {"link_angle", Json(R.link_angle)}, {"links", links},
                           {"plaquettes", Json(R.plaquettes)}, {"plaquette_angle", Json(R.plaquette_angle)},
                           {"history", Json::object({{"plaquette", Json(R.history)}})}});
  Json mean = estimate_json(R.plaquette);
  mean.set("plaquette", R.plaquette.mean);
  out.set("mean", mean);
  out.set("exact", R.exact_plaquette ? Json::object({{"plaquette", *R.exact_plaquette}}) : Json());
  out.set("loops", loops);
  out.set("charge", R.charge ? Json(*R.charge) : Json());
  out.set("acceptance", R.acceptance);
  return out;
}

inline Json turaev_viro(const Model& m, const Json& request) {
  const int r = static_cast<int>(integer(request, "r", 4, 3, 9));
  physis::quantum_gravity::Triangulation3 T;
  Json manifold;
  if (m.simplicial) {
    const auto M = mt::manifold_type(*m.simplicial);
    if (M.dimension != 3 || !M.manifold || !M.closed)
      throw std::invalid_argument("the Turaev–Viro invariant needs a closed triangulated 3-manifold" + (M.reason.empty() ? std::string() : ": " + M.reason));
    T = physis::quantum_gravity::triangulation(*m.simplicial);
    manifold = Json::object({{"vertices", static_cast<double>(m.simplicial->count(0))}, {"tetrahedra", static_cast<double>(m.simplicial->count(3))}});
  } else if (m.graph && m.graph->colours == 4) {
    T = physis::quantum_gravity::triangulation(*m.graph);
    manifold = Json::object({{"coloured_graph", true}});
  } else throw std::invalid_argument("the Turaev–Viro invariant needs a triangulated 3-manifold (its simplices) or a 4-coloured graph");
  if (T.edges > 48) throw std::invalid_argument("more than 48 edges: the state sum would take too long");
  const auto S = physis::quantum_gravity::turaev_viro(T, r);
  const double w2 = physis::symmetry::total_quantum_dimension_squared(r);
  return Json::object({{"value", S.value}, {"colourings", static_cast<double>(S.colourings)}, {"r", r}, {"edges", static_cast<double>(T.edges)},
                       {"w2", w2}, {"sphere", 1 / w2}, {"triangulation", manifold}});
}

}  // namespace atlas
