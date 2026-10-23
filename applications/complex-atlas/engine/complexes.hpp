// The Complex Atlas's engine: complexes to and from JSON, the catalogue of families, and the
// analysis of a complex, all by morphe's topology (README §4, the protocol).
//
// A complex in JSON: {"name", "title", "vertices", "edges": [[tail, head]], "faces": [{"base",
// "boundary": [±(edge + 1)]}], "cells": [[{"boundary": [[index, coefficient]]}] for k = 3, 4,
// ...], and optionally "positions", "labels", "simplices" (the facets, when the complex is
// simplicial) and "coloured_graph" ({"colours", "p", "sigma"}, when it comes from one).

#pragma once

#include <cmath>
#include <cstddef>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <morphe/topology/collapse.hpp>
#include <morphe/topology/coloured_graph.hpp>
#include <morphe/topology/constructions.hpp>
#include <morphe/topology/dual.hpp>
#include <morphe/topology/invariants.hpp>
#include <morphe/topology/manifolds.hpp>
#include <morphe/topology/surfaces.hpp>
#include <science/json.hpp>

namespace atlas {

using science::Json;
namespace mt = morphe::topology;
namespace ma = morphe::algebra;

// A complex with what it came from.
struct Model {
  mt::CWComplex X;
  std::string title;
  std::vector<std::array<double, 3>> positions;     // empty, or one per vertex
  std::optional<mt::SimplicialComplex> simplicial;  // when the complex is a triangulation
  std::optional<mt::ColouredGraph> graph;           // when it comes from a coloured graph
};

// ---- JSON -------------------------------------------------------------------------------------------------------------

inline double number(const Json& j, const char* key, double fallback) {
  const Json* v = j.find(key);
  return v && v->is_number() ? v->as_number() : fallback;
}
inline long long integer(const Json& j, const char* key, long long fallback, long long lo, long long hi) {
  const double x = number(j, key, static_cast<double>(fallback));
  if (!std::isfinite(x) || x != std::floor(x)) throw std::invalid_argument(std::string("the parameter ") + key + " must be an integer");
  if (x < static_cast<double>(lo) || x > static_cast<double>(hi))
    throw std::invalid_argument(std::string("the parameter ") + key + " must lie in " + std::to_string(lo) + ".." + std::to_string(hi));
  return static_cast<long long>(x);
}
inline std::size_t index(const Json& j, std::size_t bound, const char* what) {
  if (!j.is_number()) throw std::invalid_argument(std::string(what) + " must be a number");
  const double x = j.as_number();
  if (x < 0 || x != std::floor(x) || x >= static_cast<double>(bound)) throw std::invalid_argument(std::string(what) + " out of range");
  return static_cast<std::size_t>(x);
}

inline Json word_json(const mt::Path& p) {
  Json w = Json::array();
  for (const auto& s : p) w.push(s.sign * static_cast<long long>(s.edge + 1));
  return w;
}
inline Json letters_json(const ma::Word& w) {
  Json a = Json::array();
  for (int x : w) a.push(x);
  return a;
}
inline Json chain_json(const ma::Chain& c) {
  Json a = Json::array();
  for (const auto& e : c) a.push(Json::array({static_cast<double>(e.index), e.value.to_double()}));
  return a;
}

inline Json to_json(const Model& m) {
  const auto& X = m.X;
  Json j = Json::object({{"name", X.name}, {"title", m.title.empty() ? X.name : m.title}, {"vertices", static_cast<double>(X.vertices)}});
  Json edges = Json::array();
  for (const auto& e : X.edges) edges.push(Json::array({static_cast<double>(e[0]), static_cast<double>(e[1])}));
  j.set("edges", edges);
  Json faces = Json::array();
  for (const auto& f : X.faces) faces.push(Json::object({{"base", static_cast<double>(f.base)}, {"boundary", word_json(f.boundary)}}));
  j.set("faces", faces);
  Json cells = Json::array();
  for (const auto& level : X.higher) {
    Json l = Json::array();
    for (const auto& c : level) {
      Json b = Json::array();
      for (const auto& e : c) b.push(Json::array({static_cast<double>(e.index), e.value.to_double()}));
      l.push(Json::object({{"boundary", b}}));
    }
    cells.push(l);
  }
  j.set("cells", cells);
  if (!m.positions.empty()) {
    Json p = Json::array();
    for (const auto& x : m.positions) p.push(Json::array({x[0], x[1], x[2]}));
    j.set("positions", p);
  }
  if (!X.edge_names.empty()) j.set("labels", Json::object({{"edges", Json(X.edge_names)}}));
  if (m.simplicial) {
    Json s = Json::array();
    for (const auto& f : m.simplicial->facets()) {
      Json t = Json::array();
      for (auto v : f) t.push(static_cast<double>(v));
      s.push(t);
    }
    j.set("simplices", s);
  }
  if (m.graph) {
    Json sig = Json::array();
    for (const auto& s : m.graph->sigma) {
      Json a = Json::array();
      for (auto w : s) a.push(static_cast<double>(w));
      sig.push(a);
    }
    j.set("coloured_graph", Json::object({{"colours", static_cast<double>(m.graph->colours)}, {"p", static_cast<double>(m.graph->p)}, {"sigma", sig}}));
  }
  return j;
}

inline Model from_json(const Json& j) {
  if (!j.is_object()) throw std::invalid_argument("a complex must be a JSON object");
  Model m;
  auto& X = m.X;
  X.name = j["name"].text("your complex");
  m.title = j["title"].text(X.name);
  const double nv = number(j, "vertices", -1);
  if (!(nv >= 0) || nv != std::floor(nv) || nv > 5e6) throw std::invalid_argument("\"vertices\" must be a count of vertices");
  X.vertices = static_cast<std::size_t>(nv);
  for (const auto& e : j["edges"].is_array() ? j["edges"].as_array() : Json::Array{}) {
    if (e.size() != 2) throw std::invalid_argument("an edge is [tail, head]");
    X.edges.push_back({index(e[0], X.vertices, "an edge's tail"), index(e[1], X.vertices, "an edge's head")});
  }
  for (const auto& f : j["faces"].is_array() ? j["faces"].as_array() : Json::Array{}) {
    mt::Face face;
    face.base = f.find("base") ? index(f["base"], X.vertices, "a face's base") : 0;
    for (const auto& x : f["boundary"].is_array() ? f["boundary"].as_array() : Json::Array{}) {
      if (!x.is_number() || x.as_number() == 0 || x.as_number() != std::floor(x.as_number())) throw std::invalid_argument("a word's letters are ±(edge + 1)");
      const double v = x.as_number();
      const std::size_t e = static_cast<std::size_t>(std::fabs(v)) - 1;
      if (e >= X.edges.size()) throw std::invalid_argument("a word names an edge that is not there");
      face.boundary.push_back({e, v > 0 ? 1 : -1});
    }
    // A face without a base starts where its first step starts.
    if (!f.find("base") && !face.boundary.empty()) face.base = face.boundary[0].sign > 0 ? X.edges[face.boundary[0].edge][0] : X.edges[face.boundary[0].edge][1];
    X.faces.push_back(face);
  }
  for (const auto& level : j["cells"].is_array() ? j["cells"].as_array() : Json::Array{}) {
    std::vector<ma::Vector> cells;
    const std::size_t k = X.higher.size() + 3;
    const std::size_t below = X.cells(k - 1);
    for (const auto& c : level.is_array() ? level.as_array() : Json::Array{}) {
      ma::Vector v;
      for (const auto& e : c["boundary"].is_array() ? c["boundary"].as_array() : Json::Array{}) {
        if (e.size() != 2 || !e[1].is_number()) throw std::invalid_argument("a cell's boundary is [[index, coefficient]]");
        v = ma::axpy(v, ma::Integer(static_cast<long long>(e[1].as_number())), ma::Vector{{index(e[0], below, "a boundary cell"), 1}});
      }
      cells.push_back(v);
    }
    X.higher.push_back(cells);
  }
  if (const Json* labels = j.find("labels"))
    if ((*labels)["edges"].is_array())
      for (const auto& n : (*labels)["edges"].as_array()) X.edge_names.push_back(n.text());
  if (const Json* p = j.find("positions"); p && p->is_array() && p->size() == X.vertices)
    for (const auto& x : p->as_array()) m.positions.push_back({x[0].as_number(), x[1].as_number(), x[2].as_number()});
  if (const Json* s = j.find("simplices"); s && s->is_array() && s->size()) {
    std::vector<mt::Simplex> facets;
    for (const auto& f : s->as_array()) {
      mt::Simplex t;
      for (const auto& v : f.as_array()) t.push_back(index(v, X.vertices, "a simplex's vertex"));
      facets.push_back(t);
    }
    m.simplicial = mt::SimplicialComplex::from_facets(facets, X.vertices, X.name);
  }
  if (const Json* g = j.find("coloured_graph"); g && g->is_object()) {
    mt::ColouredGraph G;
    G.colours = static_cast<std::size_t>(number(*g, "colours", 0));
    G.p = static_cast<std::size_t>(number(*g, "p", 0));
    for (const auto& s : (*g)["sigma"].as_array()) {
      std::vector<std::size_t> perm;
      for (const auto& w : s.as_array()) perm.push_back(index(w, G.p, "a permutation's value"));
      G.sigma.push_back(perm);
    }
    G.validate();
    m.graph = G;
  }
  X.validate();
  return m;
}

// ---- the catalogue -------------------------------------------------------------------------------------------------------

inline Json param(const char* name, const char* type, double lo, double hi, Json fallback) {
  Json p = Json::object({{"name", name}, {"type", type}, {"default", fallback}});
  if (std::string(type) != "text") p.set("min", lo), p.set("max", hi);
  return p;
}

inline Json catalogue() {
  struct Family { const char* id; const char* title; const char* group; const char* description; Json params; };
  const Family families[] = {
      {"surface", "Orientable surface", "Surfaces", "The closed orientable surface of genus g: one 4g-gon glued by a₁ b₁ a₁⁻¹ b₁⁻¹ ⋯ (genus 0: the sphere as a disc glued by a a⁻¹).",
       Json::array({param("genus", "int", 0, 8, 1)})},
      {"nonorientable", "Non-orientable surface", "Surfaces", "The closed surface with k crosscaps: one 2k-gon glued by c₁ c₁ ⋯ c_k c_k (k = 1, the projective plane; k = 2, the Klein bottle).",
       Json::array({param("crosscaps", "int", 1, 8, 2)})},
      {"klein-bottle", "Klein bottle", "Surfaces", "A square glued by a b a⁻¹ b: one vertex, two edges, one 2-cell.", Json::array()},
      {"words", "Polygons glued by words", "Surfaces", "Any polygons, each a word: a lowercase letter is a side along the polygon, its capital a side against it; equal letters are one edge.",
       Json::array({param("words", "text", 0, 0, "abAB")})},
      {"dunce-hat", "Dunce hat", "Homotopy", "A triangle glued by a a a⁻¹: contractible, yet without a free edge, so it does not collapse (Zeeman 1964).", Json::array()},
      {"sphere", "Sphere", "Cell structures", "The n-sphere as a point with one n-cell.", Json::array({param("n", "int", 1, 8, 2)})},
      {"real-projective-space", "Real projective space", "Cell structures", "ℝPⁿ with one cell in each dimension, ∂e_k = (1 + (−1)^k) e_{k−1}.", Json::array({param("n", "int", 1, 8, 3)})},
      {"complex-projective-space", "Complex projective space", "Cell structures", "ℂPⁿ with one cell in each even dimension.", Json::array({param("n", "int", 1, 4, 2)})},
      {"lens-space", "Lens space", "Cell structures", "L(p, q) with one cell in each dimension, the 2-cell glued by aᵖ: π₁ = ℤ/p.", Json::array({param("p", "int", 1, 12, 5)})},
      {"presentation", "Presentation complex", "Groups", "One vertex, a loop for each letter, a 2-cell for each relator: its π₁ is the presented group. Relators are words, separated by spaces.",
       Json::array({param("relators", "text", 0, 0, "aa bbb ababababab")})},
      {"lattice-torus", "Lattice torus", "Lattices", "The periodic n × m square lattice: the torus as a product of circles, with its links and plaquettes, the arena of lattice gauge theory.",
       Json::array({param("n", "int", 1, 64, 8), param("m", "int", 1, 64, 8)})},
      {"lattice-disc", "Lattice disc", "Lattices", "The n × m square grid: a disc, on which the plaquettes of a gauge field are independent and Wilson loops obey an exact area law.",
       Json::array({param("n", "int", 1, 64, 8), param("m", "int", 1, 64, 8)})},
      {"lattice-cylinder", "Lattice cylinder", "Lattices", "A circle of n links times an interval of m: an annulus.", Json::array({param("n", "int", 1, 64, 12), param("m", "int", 1, 64, 4)})},
      {"lattice-3-torus", "Lattice 3-torus", "Lattices", "The periodic n × n × n cubic lattice, with its cubes.", Json::array({param("n", "int", 1, 10, 3)})},
      {"boundary-simplex", "Boundary of a simplex", "Triangulations", "∂Δⁿ, the (n−1)-sphere on n + 1 vertices.", Json::array({param("n", "int", 2, 7, 3)})},
      {"octahedron", "Octahedron", "Triangulations", "The boundary of the cross-polytope: a 2-sphere on 6 vertices.", Json::array()},
      {"torus-7", "Seven-vertex torus", "Triangulations", "The minimal triangulation of the torus (Möbius 1886; Császár 1949).", Json::array()},
      {"projective-plane-6", "Six-vertex projective plane", "Triangulations", "The minimal triangulation of ℝP², the icosahedron with antipodes identified.", Json::array()},
      {"three-sphere", "Suspended 3-sphere", "Triangulations", "The suspension of ∂Δ³: a 3-sphere on 6 vertices.", Json::array()},
      {"facets", "Simplicial complex from facets", "Triangulations", "Any simplicial complex, its facets as lists of vertices separated by semicolons.",
       Json::array({param("facets", "text", 0, 0, "0 1 2; 0 2 3; 0 3 1; 1 2 3")})},
      {"dual-simplex-boundary", "Dual 2-complex of ∂Δⁿ", "Spin foams", "The 2-complex dual to the boundary of the n-simplex: a vertex per facet, an edge per ridge, a 2-cell around each codimension-2 face. For n = 4, the spin-foam 2-complex of the 3-sphere.",
       Json::array({param("n", "int", 3, 6, 4)})},
      {"dual-torus-7", "Dual of the seven-vertex torus", "Spin foams", "The hexagonal torus dual to the minimal triangulation.", Json::array()},
      {"melon", "Melonic graph", "Coloured graphs", "A melonic (D+1)-coloured graph: the elementary melon with k melons inserted at random. Its 2-complex is the dual 2-skeleton of a D-sphere; the leading order of a coloured tensor model.",
       Json::array({param("D", "int", 2, 5, 3), param("insertions", "int", 0, 40, 3), param("seed", "int", 1, 999999, 1)})},
      {"coloured-graph", "Random coloured graph", "Coloured graphs", "A random (D+1)-coloured bipartite graph with 2p vertices: a Feynman graph of a group field theory, shown as its 2-complex of bicoloured cycles.",
       Json::array({param("D", "int", 2, 5, 3), param("p", "int", 1, 12, 3), param("seed", "int", 1, 999999, 1)})},
      {"coloured-pseudomanifold", "Pseudomanifold of a coloured graph", "Coloured graphs", "The D-dimensional pseudomanifold a random coloured graph encodes, as a Δ-complex: k-cells are the components of the subgraphs on D − k colours.",
       Json::array({param("D", "int", 2, 4, 3), param("p", "int", 1, 10, 3), param("seed", "int", 1, 999999, 1)})},
  };
  Json out = Json::array();
  for (const auto& f : families) out.push(Json::object({{"id", f.id}, {"title", f.title}, {"group", f.group}, {"description", f.description}, {"params", f.params}}));
  return Json::object({{"families", out}});
}

// ---- building --------------------------------------------------------------------------------------------------------------

inline Model simplicial_model(const mt::SimplicialComplex& K, std::string title) {
  Model m;
  m.X = K.cw();
  m.X.name = title;
  m.title = std::move(title);
  m.simplicial = K;
  return m;
}

inline mt::SimplicialComplex parse_facets(const std::string& text) {
  std::vector<mt::Simplex> facets;
  mt::Simplex cur;
  std::string num;
  auto flush_num = [&] {
    if (!num.empty()) cur.push_back(std::stoul(num)), num.clear();
  };
  for (char c : text + ";") {
    if (c >= '0' && c <= '9') num += c;
    else {
      flush_num();
      if (c == ';' || c == '\n') {
        if (!cur.empty()) facets.push_back(cur), cur.clear();
      } else if (c != ' ' && c != ',' && c != '\t' && c != '\r')
        throw std::invalid_argument("facets are lists of vertex numbers separated by semicolons");
    }
  }
  if (facets.empty()) throw std::invalid_argument("no facets");
  return mt::SimplicialComplex::from_facets(facets, 0, "your simplicial complex");
}

inline Model build(const std::string& family, const Json& p) {
  Model m;
  auto& X = m.X;
  constexpr double tau = 2 * std::numbers::pi;
  if (family == "surface") X = mt::orientable_surface(static_cast<std::size_t>(integer(p, "genus", 1, 0, 8)));
  else if (family == "nonorientable") X = mt::nonorientable_surface(static_cast<std::size_t>(integer(p, "crosscaps", 2, 1, 8)));
  else if (family == "klein-bottle") X = mt::klein_bottle();
  else if (family == "words") {
    const std::string w = p["words"].text("abAB");
    if (w.size() > 2000) throw std::invalid_argument("the words are too long");
    X = mt::polygon_complex(w);
  } else if (family == "dunce-hat") X = mt::dunce_hat();
  else if (family == "sphere") X = mt::sphere(static_cast<std::size_t>(integer(p, "n", 2, 1, 8)));
  else if (family == "real-projective-space") X = mt::real_projective_space(static_cast<std::size_t>(integer(p, "n", 3, 1, 8)));
  else if (family == "complex-projective-space") X = mt::complex_projective_space(static_cast<std::size_t>(integer(p, "n", 2, 1, 4)));
  else if (family == "lens-space") X = mt::lens_space(static_cast<std::size_t>(integer(p, "p", 5, 1, 12)));
  else if (family == "presentation") {
    const auto polys = mt::parse_polygons(p["relators"].text("aa bbb ababababab"));
    ma::Presentation P;
    std::map<std::string, int> gen;
    for (const auto& poly : polys) {
      ma::Word r;
      for (const auto& s : poly) {
        auto [it, fresh] = gen.emplace(s.label, static_cast<int>(gen.size()) + 1);
        if (fresh) P.names.push_back(s.label);
        r.push_back(s.sign * it->second);
      }
      P.relators.push_back(r);
    }
    P.generators = gen.size();
    X = mt::presentation_complex(P);
  } else if (family == "lattice-torus" || family == "lattice-disc" || family == "lattice-cylinder") {
    const auto n = static_cast<std::size_t>(integer(p, "n", 8, 1, 64)), k = static_cast<std::size_t>(integer(p, "m", 8, 1, 64));
    const bool torus = family == "lattice-torus", disc = family == "lattice-disc";
    X = torus ? mt::torus_lattice({n, k}) : disc ? mt::product(mt::interval(n), mt::interval(k)) : mt::product(mt::circle(n), mt::interval(k));
    const std::size_t ny = torus ? k : k + 1;
    for (std::size_t v = 0; v < X.vertices; v++) {
      const double a = static_cast<double>(v / ny), b = static_cast<double>(v % ny);
      if (torus) {
        const double th = tau * a / static_cast<double>(n), ph = tau * b / static_cast<double>(k);
        m.positions.push_back({(2 + 0.8 * std::cos(ph)) * std::cos(th), (2 + 0.8 * std::cos(ph)) * std::sin(th), 0.8 * std::sin(ph)});
      } else if (disc) {
        const double s = 4.0 / static_cast<double>(std::max(n, k));
        m.positions.push_back({(a - static_cast<double>(n) / 2) * s, (b - static_cast<double>(k) / 2) * s, 0});
      } else {
        const double th = tau * a / static_cast<double>(n);
        m.positions.push_back({1.5 * std::cos(th), 1.5 * std::sin(th), (b - static_cast<double>(k) / 2) * 3.0 / static_cast<double>(std::max<std::size_t>(k, 1))});
      }
    }
    X.name = torus ? "lattice torus " + std::to_string(n) + " × " + std::to_string(k) : disc ? "lattice disc " + std::to_string(n) + " × " + std::to_string(k)
                                                                                      : "lattice cylinder " + std::to_string(n) + " × " + std::to_string(k);
  } else if (family == "lattice-3-torus") {
    const auto n = static_cast<std::size_t>(integer(p, "n", 3, 1, 10));
    X = mt::torus_lattice({n, n, n});
  } else if (family == "boundary-simplex") {
    const auto n = static_cast<std::size_t>(integer(p, "n", 3, 2, 7));
    m = simplicial_model(mt::simplex_boundary(n), "boundary of the " + std::to_string(n) + "-simplex");
    if (n == 3) m.positions = {{{1, 1, 1}}, {{1, -1, -1}}, {{-1, 1, -1}}, {{-1, -1, 1}}};
    if (n == 2) m.positions = {{{1, 0, 0}}, {{-0.5, 0.866, 0}}, {{-0.5, -0.866, 0}}};
  } else if (family == "octahedron") {
    m = simplicial_model(mt::octahedron(), "octahedron");
    m.positions = {{{1.5, 0, 0}}, {{-1.5, 0, 0}}, {{0, 1.5, 0}}, {{0, -1.5, 0}}, {{0, 0, 1.5}}, {{0, 0, -1.5}}};
  } else if (family == "torus-7") m = simplicial_model(mt::torus_7(), "seven-vertex torus");
  else if (family == "projective-plane-6") m = simplicial_model(mt::projective_plane_6(), "six-vertex projective plane");
  else if (family == "three-sphere") m = simplicial_model(mt::suspension(mt::simplex_boundary(3)), "suspended 3-sphere");
  else if (family == "facets") {
    const std::string text = p["facets"].text("0 1 2; 0 2 3; 0 3 1; 1 2 3");
    if (text.size() > 20000) throw std::invalid_argument("too many facets");
    m = simplicial_model(parse_facets(text), "your simplicial complex");
  } else if (family == "dual-simplex-boundary") {
    const auto n = static_cast<std::size_t>(integer(p, "n", 4, 3, 6));
    X = mt::dual_two_complex(mt::simplex_boundary(n)).complex;
    X.name = "dual 2-complex of ∂Δ" + std::to_string(n);
  } else if (family == "dual-torus-7") X = mt::dual_two_complex(mt::torus_7()).complex;
  else if (family == "melon" || family == "coloured-graph" || family == "coloured-pseudomanifold") {
    const auto D = static_cast<std::size_t>(integer(p, "D", 3, 2, 5));
    mt::ColouredGraph G;
    if (family == "melon") {
      const auto k = static_cast<std::size_t>(integer(p, "insertions", 3, 0, 40));
      std::mt19937_64 g(static_cast<std::uint64_t>(integer(p, "seed", 1, 1, 999999)));
      G = mt::elementary_melon(D);
      for (std::size_t i = 0; i < k; i++) G = mt::insert_melon(G, g() % (D + 1), g() % G.p);
    } else {
      const auto pp = static_cast<std::size_t>(integer(p, "p", 3, 1, family == "coloured-graph" ? 12 : 10));
      const auto seed = static_cast<std::uint64_t>(integer(p, "seed", 1, 1, 999999));
      if (family == "coloured-pseudomanifold" && D > 4) throw std::invalid_argument("D at most 4 for the pseudomanifold");
      G = mt::random_coloured_graph(D, pp, seed);
    }
    X = family == "coloured-pseudomanifold" ? G.dual_complex() : G.two_complex();
    m.graph = G;
  } else throw std::invalid_argument("no family \"" + family + "\"");
  if (m.title.empty()) m.title = X.name;
  X.validate();
  return m;
}

// ---- analysis ------------------------------------------------------------------------------------------------------------

inline Json analyse(const Model& m, bool generators = true) {
  const auto& X = m.X;
  const std::size_t cells = X.vertices + X.edges.size() + X.faces.size();
  if (cells > 400000) throw std::invalid_argument("the complex has more than 400000 cells of dimension at most 2");
  Json a = Json::object({{"dimension", static_cast<double>(X.dimension())}, {"euler", static_cast<double>(X.euler_characteristic())}});
  a.set("counts", Json(X.counts()));
  const auto F = morphe::combinatorics::spanning_forest(X.one_skeleton());
  a.set("components", static_cast<double>(F.components()));
  const auto H = mt::homology(X, generators && cells <= 60000);
  Json hj = Json::array();
  for (const auto& g : H) {
    Json gens = Json::array(), tors = Json::array();
    for (const auto& t : g.torsion) tors.push(t.to_double());
    for (const auto& z : g.free_generators) gens.push(Json::object({{"order", 0}, {"chain", chain_json(z)}}));
    for (std::size_t i = 0; i < g.torsion_generators.size(); i++)
      gens.push(Json::object({{"order", g.torsion[i].to_double()}, {"chain", chain_json(g.torsion_generators[i])}}));
    hj.push(Json::object({{"dimension", static_cast<double>(g.dimension)}, {"rank", static_cast<double>(g.rank)}, {"torsion", tors}, {"generators", gens}, {"text", g.str()}}));
  }
  a.set("homology", hj);
  Json co = Json::array();
  for (const auto& c : ma::cohomology(H)) {
    Json tors = Json::array();
    for (const auto& t : c.torsion) tors.push(t.to_double());
    co.push(Json::object({{"dimension", static_cast<double>(c.dimension)}, {"rank", static_cast<double>(c.rank)}, {"torsion", tors}}));
  }
  a.set("cohomology", co);
  a.set("betti", Json::object({{"2", Json(mt::betti(X, 2))}, {"3", Json(mt::betti(X, 3))}}));
  // The fundamental group of the component of vertex 0.
  if (X.vertices > 0 && X.edges.size() <= 20000) {
    const auto G = mt::fundamental_group(X, 0);
    auto present = [](const ma::Presentation& P, const std::vector<mt::Path>& loops) {
      Json gens = Json::array(), rels = Json::array(), ls = Json::array();
      for (std::size_t i = 0; i < P.generators; i++) gens.push(P.name(i));
      for (const auto& r : P.relators) rels.push(letters_json(r));
      for (const auto& l : loops) ls.push(word_json(l));
      return Json::object({{"generators", gens}, {"relators", rels}, {"loops", ls}, {"text", P.str()}});
    };
    Json fg = present(G.presentation, G.loops);
    fg.set("simplified", present(G.simplified.presentation, G.simplified_loops));
    Json tors = Json::array();
    for (const auto& t : G.abelianisation.torsion) tors.push(t.to_double());
    fg.set("abelianisation", Json::object({{"rank", static_cast<double>(G.abelianisation.rank)}, {"torsion", tors}, {"text", G.abelianisation.str()}}));
    a.set("fundamental_group", fg);
  } else a.set("fundamental_group", nullptr);
  if (X.dimension() == 2) {
    const auto S = mt::classify_surface(X);
    if (S.surface) {
      Json o = Json::array();
      for (int x : S.orientation) o.push(x);
      a.set("surface", Json::object({{"closed", S.closed}, {"orientable", S.orientable}, {"genus", static_cast<double>(S.genus)}, {"crosscaps", static_cast<double>(S.crosscaps)},
                                     {"boundary_components", static_cast<double>(S.boundary_components)}, {"euler", static_cast<double>(S.euler)}, {"name", S.name},
                                     {"normal_form", letters_json(S.normal_form)}, {"orientation", o}}));
    } else {
      a.set("surface", nullptr);
      a.set("surface_reason", S.reason);
    }
  } else a.set("surface", nullptr);
  if (m.simplicial) {
    const auto M = mt::manifold_type(*m.simplicial);
    a.set("manifold", Json::object({{"pseudomanifold", M.pseudomanifold}, {"manifold", M.manifold}, {"exact", M.exact}, {"closed", M.closed},
                                    {"dimension", static_cast<double>(M.dimension)}, {"reason", M.reason}}));
  } else a.set("manifold", nullptr);
  const auto C = mt::collapse(X);
  if (C.attempted) a.set("collapse", Json::object({{"collapsible", C.collapsible}, {"remaining", Json(C.remaining)}, {"steps", static_cast<double>(C.steps)}}));
  else a.set("collapse", nullptr);
  if (m.graph) {
    const auto& G = *m.graph;
    const bool connected = G.connected();
    a.set("coloured_graph", Json::object({{"colours", static_cast<double>(G.colours)}, {"vertices", static_cast<double>(G.vertices())},
                                          {"faces", static_cast<double>(G.face_count())}, {"bubbles", static_cast<double>(G.bubbles())},
                                          {"connected", connected}, {"jackets", Json(G.jacket_genera())},
                                          {"degree", connected ? Json(static_cast<double>(G.degree())) : Json()}, {"melonic", connected && G.colours >= 4 && G.degree() == 0}}));
  } else a.set("coloured_graph", nullptr);
  return a;
}

}  // namespace atlas
