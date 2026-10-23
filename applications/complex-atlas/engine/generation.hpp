// The Complex Atlas's engine: the generation of spacetime by field theories (physis decision
// 0006; the infrastructure's ADR 0040). Theories and their Feynman expansions, the complexes
// of their graphs, their large-N solutions and condensates, random geometries of the
// continuum and their dimensions, and spin foams.

#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#include <optional>

#include <physis/quantum_gravity/condensate.hpp>
#include <physis/quantum_gravity/random_geometry.hpp>
#include <physis/quantum_gravity/spin_foam.hpp>
#include <physis/spacetime/dimensions.hpp>

#include "complexes.hpp"

namespace atlas {

// ---- coloured graphs in JSON: {"colours", "p", "sigma"} --------------------------------------------------------------

inline Json coloured_graph_json(const mt::ColouredGraph& G) {
  Json sig = Json::array();
  for (const auto& s : G.sigma) {
    Json a = Json::array();
    for (auto w : s) a.push(static_cast<double>(w));
    sig.push(a);
  }
  return Json::object({{"colours", static_cast<double>(G.colours)}, {"p", static_cast<double>(G.p)}, {"sigma", sig}});
}

inline mt::ColouredGraph coloured_graph_from_json(const Json& g) {
  if (!g.is_object()) throw std::invalid_argument("a coloured graph is {\"colours\", \"p\", \"sigma\"}");
  mt::ColouredGraph G;
  G.colours = static_cast<std::size_t>(number(g, "colours", 0));
  G.p = static_cast<std::size_t>(number(g, "p", 0));
  if (G.colours > 8 || G.p > 4096) throw std::invalid_argument("the coloured graph is too large");
  for (const auto& s : g["sigma"].as_array()) {
    std::vector<std::size_t> perm;
    for (const auto& w : s.as_array()) perm.push_back(index(w, G.p, "a permutation's value"));
    G.sigma.push_back(perm);
  }
  G.validate();
  return G;
}

// The 2-complex (vertices, edges and bicoloured cycles) or the pseudomanifold of a coloured graph.
inline Json graph_complex(const Json& request) {
  const auto G = coloured_graph_from_json(request["coloured_graph"]);
  const std::string view = request["view"].text("2-complex");
  Model m;
  if (view == "pseudomanifold") {
    if (G.colours < 3 || G.colours > 5) throw std::invalid_argument("the pseudomanifold of a graph with 3 to 5 colours");
    m.X = G.dual_complex();
  } else if (view == "2-complex") m.X = G.two_complex();
  else throw std::invalid_argument("the view is \"2-complex\" or \"pseudomanifold\"");
  m.graph = G;
  m.title = m.X.name;
  return Json::object({{"complex", to_json(m)}});
}

// ---- the condensate: cosmology of a group field theory's condensed phase ----------------------------------------------

inline Json condensate(const Json& request) {
  namespace qg = physis::quantum_gravity;
  qg::Universe U;
  U.G = number(request, "G", 1.0);
  if (!(U.G > 0 && U.G < 1e3)) throw std::invalid_argument("G must be positive");
  U.m2 = number(request, "m2", -1.0);
  U.v0 = number(request, "v0", 1.0);
  for (const auto& s : request["spins"].is_array() ? request["spins"].as_array() : Json::Array{}) {
    qg::CondensateMode md;
    md.j = static_cast<double>(integer(s, "two_j", 1, 1, 200)) / 2;
    md.Q = number(s, "Q", 1.0);
    md.E = number(s, "E", 0.0);
    md.phi_bounce = number(s, "phi_bounce", 0.0);
    U.modes.push_back(md);
  }
  if (U.modes.empty() || U.modes.size() > 16) throw std::invalid_argument("one to sixteen condensate modes (spins)");
  double lo = -1.5, hi = 1.5;
  if (request["phi"].is_array() && request["phi"].size() == 2) lo = request["phi"][0].as_number(), hi = request["phi"][1].as_number();
  const auto points = static_cast<std::size_t>(integer(request, "points", 401, 2, 20001));
  const auto C = qg::evolve(U, lo, hi, points, request["integrate"].is_bool() && request["integrate"].as_bool());
  Json bounce = C.bounces ? Json::object({{"phi", C.bounce_phi}, {"V", C.bounce_volume}, {"rho", C.bounce_density}, {"rho_c", C.critical_density}}) : Json();
  const std::vector<double> classical_h2(C.phi.size(), C.classical_rate / 9);
  return Json::object({{"phi", Json(C.phi)}, {"V", Json(C.volume)}, {"H", Json(C.hubble)}, {"rho", Json(C.density)}, {"acceleration", Json(C.acceleration)},
                       {"classical", Json::object({{"V", Json(C.classical_expanding)}, {"V_contracting", Json(C.classical_contracting)}})},
                       {"bounce", bounce}, {"pi_phi", C.pi_phi}, {"rho_c", C.critical_density}, {"m2", U.mass2()},
                       {"friedmann", Json::object({{"rho", Json(C.density)}, {"H2", Json(C.hubble_squared)}, {"effective", Json(C.friedmann)},
                                                   {"classical_H2", Json(classical_h2)}})}});
}

// ---- emergent dimensions of a complex -----------------------------------------------------------------------------------

inline Json dimensions(const Model& m, const Json& request) {
  namespace ps = physis::spacetime;
  const std::string which = request["graph"].text("skeleton");
  ps::Graph g;
  if (which == "skeleton") g = ps::skeleton_graph(m.X);
  else if (which == "dual") g = ps::dual_graph(m.X);
  else if (which == "coloured" && m.graph) g = ps::coloured_graph(*m.graph);
  else throw std::invalid_argument("the graph is \"skeleton\", \"dual\" or, for a coloured graph, \"coloured\"");
  if (g.vertices() < 2) throw std::invalid_argument("the graph has fewer than two vertices");
  ps::SpectralOptions so;
  so.steps = static_cast<std::size_t>(integer(request, "steps", 200, 2, 20000));
  so.starts = static_cast<std::size_t>(integer(request, "starts", 16, 0, 4096));
  so.seed = static_cast<std::uint64_t>(integer(request, "seed", 1, 0, 1ll << 52));
  if (static_cast<double>(so.steps) * static_cast<double>(so.starts ? so.starts : g.vertices()) * static_cast<double>(g.vertices()) > 4e9)
    throw std::invalid_argument("too much diffusion for one request: fewer steps or starts");
  ps::HausdorffOptions ho;
  ho.centres = so.starts;
  ho.seed = so.seed;
  const auto S = ps::spectral_dimension(g, so);
  const auto H = ps::hausdorff_dimension(g, ho);
  // The curves where the estimators are defined (their centred differences need room at the ends).
  std::vector<double> t, P, ds, r, ball, dh;
  for (std::size_t i = 0; i < S.t.size(); i++)
    if (std::isfinite(S.ds[i])) t.push_back(S.t[i]), P.push_back(S.P[i]), ds.push_back(S.ds[i]);
  for (std::size_t i = 0; i < H.r.size(); i++)
    if (std::isfinite(H.dh[i])) r.push_back(H.r[i]), ball.push_back(H.volume[i]), dh.push_back(H.dh[i]);
  return Json::object({{"graph", which}, {"vertices", static_cast<double>(g.vertices())},
                       {"spectral", Json::object({{"t", Json(t)}, {"P", Json(P)}, {"ds", Json(ds)}})},
                       {"hausdorff", Json::object({{"r", Json(r)}, {"ball", Json(ball)}, {"dh", Json(dh)}})}});
}

// ---- random geometries of the continuum ---------------------------------------------------------------------------------

inline Json random_geometry(const Json& request) {
  namespace qg = physis::quantum_gravity;
  const std::string kind = request["kind"].text("quadrangulation");
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1, 0, 1ll << 52));
  Model m;
  Json extra = Json::object();
  if (kind == "quadrangulation") {
    const auto n = static_cast<std::size_t>(integer(request, "size", 500, 1, 50000));
    const auto Q = qg::random_quadrangulation(n, seed);
    m.X = Q.complex;
    m.X.name = "random planar quadrangulation of " + std::to_string(n) + " faces";
    Json d = Json::array();
    for (int l : Q.label) d.push(l);
    extra.set("distance", d);
    extra.set("root", static_cast<double>(Q.extra));
  } else if (kind == "melonic") {
    const auto D = static_cast<std::size_t>(integer(request, "D", 3, 2, 5));
    const auto p = static_cast<std::size_t>(integer(request, "size", 20, 1, 4000));
    m.graph = qg::random_melonic_graph(D, p, seed);
    m.X = request["view"].text("2-complex") == "pseudomanifold" ? m.graph->dual_complex() : m.graph->two_complex();
    m.X.name = "random melonic " + std::to_string(D + 1) + "-coloured graph of " + std::to_string(2 * p) + " vertices";
  } else throw std::invalid_argument("the kind is \"quadrangulation\" or \"melonic\"");
  m.title = m.X.name;
  Json out = Json::object({{"complex", to_json(m)}});
  for (const auto& [k, v] : extra.as_object()) out.set(k, v);
  return out;
}

// ---- spin foams --------------------------------------------------------------------------------------------------------------

inline Json spin_foam(const Model& m, const Json& request) {
  namespace qg = physis::quantum_gravity;
  std::optional<qg::SpinFoam> built;
  if (m.simplicial) {
    const auto M = mt::manifold_type(*m.simplicial);
    if (M.dimension != 3 || !M.manifold || !M.closed)
      throw std::invalid_argument("a spin foam needs a closed triangulated 3-manifold" + (M.reason.empty() ? std::string() : ": " + M.reason));
    built = qg::spin_foam(*m.simplicial);
  } else if (m.graph && m.graph->colours == 4) built = qg::spin_foam(*m.graph);
  else throw std::invalid_argument("a spin foam needs a closed triangulated 3-manifold (its simplices) or a 4-coloured graph");
  const qg::SpinFoam& F = *built;
  const int r = static_cast<int>(integer(request, "r", 5, 0, 12));
  const std::string mode = request["mode"].text("sample");
  std::vector<int> spins(F.faces(), 0);
  Json statistics;
  if (mode == "sample") {
    if (r < 3) throw std::invalid_argument("sampling needs a level r >= 3, where the sum is finite");
    qg::SampleOptions o;
    o.r = r;
    o.sweeps = static_cast<std::size_t>(integer(request, "sweeps", 2000, 1, 200000));
    o.thermalise = static_cast<std::size_t>(integer(request, "thermalise", 200, 0, 200000));
    o.seed = static_cast<std::uint64_t>(integer(request, "seed", 1, 0, 1ll << 52));
    if (static_cast<double>(F.faces()) * static_cast<double>(o.sweeps + o.thermalise) > 5e7) throw std::invalid_argument("too many sweeps for this foam");
    const auto S = qg::sample(F, o);
    spins = S.last;
    Json mean = Json::array();
    for (double x : S.mean_spin) mean.push(2 * x);
    Json hist = Json::array();
    for (const auto& h : S.histogram) hist.push(Json(h));
    statistics = Json::object({{"mean_two_j", mean}, {"histogram", hist}, {"sign_mean", S.mean_sign}, {"positive", S.positive},
                               {"acceptance", S.acceptance}, {"samples", static_cast<double>(o.sweeps)}, {"log_amplitude", Json(S.history)}});
  } else if (mode == "uniform") {
    std::fill(spins.begin(), spins.end(), static_cast<int>(integer(request, "two_j", 0, 0, 200)));
  } else if (mode == "given") {
    const Json& g = request["spins"];
    if (!g.is_array() || g.size() != F.faces()) throw std::invalid_argument("give one two_j per face of the foam");
    for (std::size_t f = 0; f < F.faces(); f++) spins[f] = static_cast<int>(g[f].as_number());
  } else throw std::invalid_argument("the mode is \"sample\", \"uniform\" or \"given\"");
  const auto A = qg::amplitude(F, spins, r);
  if (!A.admissible) throw std::invalid_argument("these spins are not admissible at this level (every triangle needs a triad with j1 + j2 + j3 an integer, within r − 2)");
  Json vertices = Json::array();
  for (std::size_t v = 0; v < F.vertices(); v++) {
    const auto Q = qg::vertex_geometry(F, spins, v, r);
    Json faces = Json::array(), two_j = Json::array(), links = Json::array();
    for (int k = 0; k < 6; k++) {
      faces.push(static_cast<double>(F.vertex_faces[v][static_cast<std::size_t>(k)]));
      two_j.push(Q.two_j[static_cast<std::size_t>(k)]);
      links.push(Json::array({Q.links[static_cast<std::size_t>(k)][0], Q.links[static_cast<std::size_t>(k)][1]}));
    }
    Json lengths = Json::array(), dihedral = Json::array();
    for (int k = 0; k < 6; k++) lengths.push(Q.lengths[static_cast<std::size_t>(k)]), dihedral.push(Q.dihedral[static_cast<std::size_t>(k)]);
    vertices.push(Json::object({{"faces", faces}, {"spins", two_j}, {"links", links}, {"amplitude", Q.six_j}, {"asymptotic", Q.ponzano_regge}, {"envelope", Q.envelope},
                                {"tetrahedron", Json::object({{"lengths", lengths}, {"volume", Q.volume}, {"euclidean", Q.euclidean}, {"dihedral", dihedral}})}}));
  }
  Json face_edge = Json::array();
  for (std::size_t f = 0; f < F.faces(); f++) {
    const std::size_t e = F.face_edge[f];
    if (m.simplicial && e < m.simplicial->count(1)) {
      const auto& s = m.simplicial->simplices[1][e];
      face_edge.push(Json::array({static_cast<double>(s[0]), static_cast<double>(s[1])}));
    } else face_edge.push(static_cast<double>(e));
  }
  const auto R = qg::regge_geometry(F, spins);
  Json regge;
  if (R.complete) {
    Json deficits = Json::array();
    for (std::size_t f = 0; f < F.faces(); f++) deficits.push(R.deficit[F.face_edge[f]]);
    regge = Json::object({{"deficits", deficits}, {"action", R.action}});
  }
  Json tv;
  if (r >= 3 && F.triangulation.edges <= 30) tv = qg::turaev_viro(F.triangulation, r).value;
  Model foam;
  foam.X = F.complex;
  foam.X.name = "spin foam of " + (m.X.name.empty() ? std::string("the complex") : m.X.name);
  foam.title = foam.X.name;
  foam.graph = m.graph;
  Json sp = Json::array();
  for (int x : spins) sp.push(x);
  return Json::object({{"foam", to_json(foam)}, {"r", r}, {"face_spins", sp}, {"face_edge", face_edge}, {"amplitude", A.value},
                       {"sign", A.value > 0 ? 1 : A.value < 0 ? -1 : 0}, {"vertices", vertices}, {"statistics", statistics}, {"tv", tv}, {"regge", regge},
                       {"euclidean_tetrahedra", static_cast<double>(R.euclidean)}});
}

}  // namespace atlas
