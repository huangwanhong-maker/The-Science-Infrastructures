// Formalism II of the manuscript on physis's stranded field theories and morphe's complexes: the
// benzene cycle and its filled 2-complex, the two-vertex complex of Case Study II, the tetrahedral
// interaction on fields over finite groups with its symmetries and order parameters, and the Wick
// expansion whose Feynman graphs generate the 2-complexes.

#pragma once

#include <map>
#include <string>
#include <vector>

#include <physis/quantum_gravity/stranded.hpp>

#include "molecules.hpp"

namespace chem {

namespace pq = physis::quantum_gravity;
namespace ma = morphe::algebra;
namespace mt = morphe::topology;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

inline const ma::FiniteGroup& group_named(const std::string& name) {
  static const ma::FiniteGroup one = ma::cyclic(1), Z2 = ma::cyclic(2), Z3 = ma::cyclic(3), S3 = ma::symmetric(3);
  if (name == "1") return one;
  if (name == "Z2") return Z2;
  if (name == "Z3") return Z3;
  if (name == "S3") return S3;
  throw std::invalid_argument("no group " + name + " (1, Z2, Z3, S3)");
}

// A complex as JSON: cells, attaching words, integer boundary matrices, homology.
inline Json complex_json(const mt::CWComplex& X) {
  Json edges = Json::array(), faces = Json::array();
  for (const auto& e : X.edges) edges.push(Json::array({static_cast<double>(e[0]), static_cast<double>(e[1])}));
  for (const auto& f : X.faces) {
    Json word = Json::array();
    for (const auto& s : f.boundary) word.push(Json::array({static_cast<double>(s.edge), s.sign}));
    faces.push(Json::object({{"base", static_cast<double>(f.base)}, {"word", word}}));
  }
  const auto C = X.chain_complex();
  Json boundaries = Json::array();
  bool zero = true;
  for (std::size_t k = 1; k < C.boundary.size(); k++) {
    const auto& B = C.boundary[k];
    Json rows = Json::array();
    std::vector<std::vector<long long>> dense(B.rows, std::vector<long long>(B.columns.size(), 0));
    for (std::size_t j = 0; j < B.columns.size(); j++)
      for (const auto& e : B.columns[j]) dense[e.index][j] = e.value.to_int64();
    for (const auto& r : dense) {
      Json row = Json::array();
      for (auto v : r) row.push(static_cast<double>(v));
      rows.push(row);
    }
    boundaries.push(rows);
  }
  // ∂₁∂₂ = 0, entry by entry.
  if (C.boundary.size() > 2) {
    const auto P = C.boundary[1] * C.boundary[2];
    for (const auto& col : P.columns)
      for (const auto& e : col) zero = zero && e.value.is_zero();
  }
  Json homology = Json::array();
  for (const auto& h : mt::homology(X, false)) homology.push(Json::object({{"rank", static_cast<double>(h.rank)}, {"group", h.str()}}));
  Json counts = Json::array();
  for (auto c : X.counts()) counts.push(static_cast<double>(c));
  return Json::object({{"name", X.name}, {"vertices", static_cast<double>(X.vertices)}, {"edges", edges}, {"faces", faces}, {"counts", counts},
                       {"euler", static_cast<double>(X.euler_characteristic())}, {"boundaries", boundaries}, {"boundaryOfBoundary", zero},
                       {"homology", homology}});
}

inline mt::CWComplex hexagon(bool filled) {
  mt::CWComplex X;
  X.name = filled ? "the six-cycle with one 2-cell" : "the six-cycle";
  X.vertices = 6;
  for (std::size_t i = 0; i < 6; i++) X.edges.push_back({i, (i + 1) % 6});
  if (filled) {
    mt::Face f;
    f.base = 0;
    for (std::size_t i = 0; i < 6; i++) f.boundary.push_back({i, 1});
    X.faces.push_back(f);
  }
  X.validate();
  return X;
}

// The identity gluing of two tetrahedral vertices: factor i of copy 0 to factor i of copy 1.
inline pq::StrandedGraph dipole() {
  pq::Pairing p(8);
  for (std::size_t i = 0; i < 4; i++) p[i] = 4 + i, p[4 + i] = i;
  auto S = pq::stranded_graph(pq::tetrahedral_vertex(), 2, p);
  S.complex.name = "the two-vertex complex of Case Study II";
  return S;
}

inline Json complexes() {
  const auto D = dipole();
  return Json::object({{"ring", complex_json(hexagon(false))}, {"filled", complex_json(hexagon(true))}, {"dipole", complex_json(D.complex)}});
}

// ---- the field -------------------------------------------------------------------------------------------------------------------

inline std::string point_name(std::size_t p, std::size_t q) {
  const auto x = pq::coordinates(p, q, 3);
  return "(" + std::to_string(x[0]) + "," + std::to_string(x[1]) + "," + std::to_string(x[2]) + ")";
}

// The tetrahedral interaction over G: its value on a given field (or the paper's two test fields),
// its slot classes and translation symmetries, and, over ℤ₂, its values along the eight
// characters, the order parameters η_s and the effect of the gauge projection.
inline Json field(const std::string& group, const Json& values) {
  const auto& G = group_named(group);
  const auto v = pq::tetrahedral_vertex();
  const std::size_t q = G.order(), size = pq::field_size(q, 3);
  std::vector<Rational> psi(size, Rational(1));
  if (values.is_array() && values.size() == size)
    for (std::size_t p = 0; p < size; p++) psi[p] = Rational(static_cast<long long>(std::llround(values[p].as_number())));
  Json out = Json::object({{"group", G.name}, {"order", static_cast<double>(q)}, {"vertex", Json(v.words)}});
  Json pts = Json::array(), vals = Json::array();
  for (std::size_t p = 0; p < size; p++) pts.push(point_name(p, q)), vals.push(psi[p].to_double());
  out.set("points", pts);
  out.set("psi", vals);
  out.set("value", pq::interaction(v, G, psi).to_double());
  out.set("indices", Json(std::vector<std::string>{"a", "b", "c", "d", "e", "f"}));
  Json classes = Json::array();
  for (const auto& c : pq::slot_classes(v)) {
    Json cj = Json::array();
    for (auto s : c) cj.push(static_cast<double>(s + 1));
    classes.push(cj);
  }
  out.set("slotClasses", classes);
  // Invariance under every slot translation (g₁, g₂, g₃), decided on a generic field (fixed values
  // with no symmetry of their own), so that the table shows the symmetries of V, not of the field.
  Json table = Json::array();
  std::vector<Rational> generic(size);
  for (std::size_t p = 0; p < size; p++) generic[p] = Rational(static_cast<long long>((p * 7 + 3) % 11) - 5);
  const Rational v0 = pq::interaction(v, G, generic);
  for (std::uint32_t a = 0; a < q; a++)
    for (std::uint32_t b = 0; b < q; b++)
      for (std::uint32_t c = 0; c < q; c++)
        table.push(Json::object({{"g", Json::array({static_cast<double>(a), static_cast<double>(b), static_cast<double>(c)})},
                                 {"invariant", pq::interaction(v, G, pq::translate(generic, G, {a, b, c})) == v0}}));
  out.set("translations", table);
  if (q == 2) {
    // The characters (−1)^{s·x}: V along each, the order parameter η_s = 2⁻³ Σ (−1)^{s·x} ψ(x), its
    // parity under the simultaneous translation, and η_s after the gauge projection.
    const auto Ppsi = pq::gauge_project(psi, G, 3);
    Json modes = Json::array();
    for (std::size_t s = 0; s < 8; s++) {
      const auto sv = pq::coordinates(s, 2, 3);
      std::vector<Rational> chi(8);
      Rational eta = 0, eta_p = 0;
      for (std::size_t p = 0; p < 8; p++) {
        const auto x = pq::coordinates(p, 2, 3);
        chi[p] = (sv[0] * x[0] + sv[1] * x[1] + sv[2] * x[2]) % 2 ? Rational(-1) : Rational(1);
        eta = eta + chi[p] * psi[p];
        eta_p = eta_p + chi[p] * Ppsi[p];
      }
      modes.push(Json::object({{"s", Json::array({static_cast<double>(sv[0]), static_cast<double>(sv[1]), static_cast<double>(sv[2])})},
                               {"odd", (sv[0] + sv[1] + sv[2]) % 2 == 1}, {"V", pq::interaction(v, G, chi).to_double()},
                               {"eta", (eta / Rational(8)).to_double()}, {"etaProjected", (eta_p / Rational(8)).to_double()}}));
    }
    out.set("modes", modes);
    Json projected = Json::array();
    for (const auto& x : Ppsi) projected.push(x.to_double());
    out.set("projected", projected);
    // A field on which V is negative: V is even and quartic, so unbounded below along it.
    const std::vector<long long> neg{-2, -2, -2, -2, 2, -2, 1, -2};
    std::vector<Rational> nf(8);
    for (int p = 0; p < 8; p++) nf[p] = Rational(neg[p]);
    out.set("negative", Json::object({{"psi", Json(std::vector<double>(neg.begin(), neg.end()))}, {"V", pq::interaction(v, G, nf).to_double()}}));
  }
  return out;
}

// ---- the Wick expansion ------------------------------------------------------------------------------------------------------------

inline Json wick(std::size_t n, const std::string& group) {
  if (n < 1 || n > 3) throw std::invalid_argument("wick: n is 1, 2 or 3");
  const auto& G = group_named(group);
  const auto v = pq::tetrahedral_vertex();
  const auto classes = pq::feynman_classes(v, n);
  Json cl = Json::array();
  Integer sum = 0;
  const Integer q(static_cast<long long>(G.order()));
  for (const auto& c : classes) {
    const auto& X = c.graph.complex;
    const Integer amp = morphe::arithmetic::pow(q, X.faces.size());
    sum += amp * Integer(static_cast<long long>(c.pairings));
    Json betti = Json::array();
    for (auto b : mt::betti(X)) betti.push(static_cast<double>(b));
    Json rep = Json::array();
    for (auto p : c.representative) rep.push(static_cast<double>(p));
    cl.push(Json::object({{"pairings", static_cast<double>(c.pairings)}, {"automorphisms", static_cast<double>(c.automorphisms)},
                          {"faces", static_cast<double>(X.faces.size())}, {"connected", c.graph.connected}, {"betti", betti},
                          {"amplitude", amp.str()}, {"bf", pq::gft_amplitude(X, G, true).str()}, {"pairing", rep}, {"complex", complex_json(X)}}));
  }
  Json out = Json::object({{"n", static_cast<double>(n)}, {"group", G.name}, {"classes", cl}, {"sum", sum.str()}});
  // The brute-force moment where it is affordable: |G|^{6n} index assignments.
  double work = 1;
  for (std::size_t k = 0; k < 6 * n; k++) work *= static_cast<double>(G.order());
  if (work <= 3e6) out.set("moment", pq::gaussian_moment(v, G, n).str());
  if (n == 2) {
    // The 24 gluings factor to factor.
    std::map<std::size_t, std::size_t> by;
    std::vector<std::size_t> sigma{0, 1, 2, 3};
    Json glue = Json::array();
    do {
      pq::Pairing p(8);
      for (std::size_t i = 0; i < 4; i++) p[i] = 4 + sigma[i], p[4 + sigma[i]] = i;
      const auto S = pq::stranded_graph(v, 2, p);
      by[S.faces()]++;
      glue.push(Json::object({{"sigma", Json::array({static_cast<double>(sigma[0]), static_cast<double>(sigma[1]), static_cast<double>(sigma[2]), static_cast<double>(sigma[3])})},
                              {"faces", static_cast<double>(S.faces())}}));
    } while (std::next_permutation(sigma.begin(), sigma.end()));
    out.set("gluings", glue);
    Json summary = Json::object();
    for (const auto& [f, k] : by) summary.set(std::to_string(f), static_cast<double>(k));
    out.set("gluingFaces", summary);
  }
  out.set("vertexAutomorphisms", static_cast<double>(pq::vertex_automorphisms(v).size()));
  return out;
}

}  // namespace chem
