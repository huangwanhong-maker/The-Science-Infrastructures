// Formalism I of the manuscript as executable grammars on morphe's graph rewriting: the ethanol
// derivation of its Fig. grammar-syntax, the labelled glycine backbones of Figs. grammar-growth and
// grammar-derivation-tree, the tetrahedral stereocentre of Case Study I, and the acid–base edit of
// Fig. proton-graph-edit; with their derivation spaces, derivation counts, orbits, equivariance,
// and material accounting.

#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>

#include <morphe/formal/graph_rewriting.hpp>
#include <physis/chemistry/kinetics.hpp>

#include "molecules.hpp"

namespace chem {

namespace rw = morphe::formal::rewriting;

inline Json graph_json(const rw::Graph& G) {
  Json v = Json::array(), e = Json::array();
  for (const auto& l : G.vertex) v.push(l);
  for (const auto& ed : G.edges) e.push(Json::array({static_cast<double>(ed.u), static_cast<double>(ed.v), ed.label}));
  return Json::object({{"vertices", v}, {"edges", e}, {"directed", G.directed}});
}

inline Json rule_json(const rw::Rule& q) {
  Json kv = Json::array(), ke = Json::array();
  for (const auto& [l, r] : q.kept_vertices) kv.push(Json::array({static_cast<double>(l), static_cast<double>(r)}));
  for (const auto& [l, r] : q.kept_edges) ke.push(Json::array({static_cast<double>(l), static_cast<double>(r)}));
  return Json::object({{"name", q.name}, {"left", graph_json(q.left)}, {"right", graph_json(q.right)}, {"keptVertices", kv},
                       {"keptEdges", ke}, {"automorphismsOfL", static_cast<double>(rw::automorphisms(q.left))}, {"conditional", static_cast<bool>(q.condition)}});
}

// A derivation space as JSON: its classes, arcs and the derivations reaching each class.
inline Json space_json(const rw::DerivationSpace& D, const std::function<bool(const rw::Graph&)>& terminal) {
  Json classes = Json::array(), arcs = Json::array();
  for (std::size_t c = 0; c < D.graphs.size(); c++) {
    Json ways = Json::array();
    for (auto w : D.ways[c]) ways.push(static_cast<double>(w));
    classes.push(Json::object({{"graph", graph_json(D.graphs[c])}, {"depth", static_cast<double>(D.depth[c])}, {"terminal", terminal(D.graphs[c])},
                               {"derivations", static_cast<double>(D.derivations_to(c))}, {"ways", ways}}));
  }
  for (const auto& a : D.arcs)
    arcs.push(Json::object({{"from", static_cast<double>(a.from)}, {"to", static_cast<double>(a.to)}, {"rule", static_cast<double>(a.rule)},
                            {"multiplicity", static_cast<double>(a.multiplicity)}}));
  return Json::object({{"classes", classes}, {"arcs", arcs}});
}

inline bool has_label_prefix(const rw::Graph& G, const std::string& prefix) {
  for (const auto& l : G.vertex)
    if (l.rfind(prefix, 0) == 0) return true;
  return false;
}

// ---- ethanol (Fig. grammar-syntax) -------------------------------------------------------------------------------------------

// Nonterminals are vertices labelled A0 and Tail1, a port an edge labelled "port". Terminal
// fragments abbreviate their hydrogens (CH3, CH2, OH), as the figure does.
inline std::vector<rw::Rule> ethanol_rules() {
  rw::Rule p1;
  p1.name = "p1: A0() → CH3–Tail1";
  p1.left.add_vertex("A0");
  const std::size_t m = p1.right.add_vertex("CH3"), t = p1.right.add_vertex("Tail1");
  p1.right.add_edge(m, t, "port");
  p1.validate();
  rw::Rule p2;
  p2.name = "p2: Tail1(u) → u–CH2–OH";
  const std::size_t lu = p2.left.add_vertex(rw::any_label), lt = p2.left.add_vertex("Tail1");
  p2.left.add_edge(lu, lt, "port");
  const std::size_t ru = p2.right.add_vertex(rw::any_label), c = p2.right.add_vertex("CH2"), o = p2.right.add_vertex("OH");
  p2.right.add_edge(ru, c, "single");
  p2.right.add_edge(c, o, "single");
  p2.kept_vertices = {{lu, ru}};
  p2.validate();
  return {p1, p2};
}

inline Json ethanol() {
  const auto rules = ethanol_rules();
  rw::Graph S;
  S.add_vertex("A0");
  const auto D = rw::explore(rules, S, 4);
  const auto terminal = [](const rw::Graph& G) { return !has_label_prefix(G, "A0") && !has_label_prefix(G, "Tail"); };
  Json rj = Json::array();
  for (const auto& q : rules) rj.push(rule_json(q));
  return Json::object({{"title", "Ethanol from a start symbol (Fig. grammar-syntax)"}, {"rules", rj}, {"space", space_json(D, terminal)}});
}

// ---- glycine backbones (Figs. grammar-growth and grammar-derivation-tree) -------------------------------------------------------

// Each unit is a vertex Gly_i with an amine terminus N• and a carboxyl terminus C• (• free). The
// coupling rule joins a free C• to a free N• of another chain into an amide C–N, both termini now
// bound, and creates one water: the net bookkeeping of the figure.
inline rw::Rule coupling_rule() {
  rw::Rule q;
  q.name = "couple: C• + N• → C–N + H2O (different chains)";
  q.left.directed = q.right.directed = true;
  const std::size_t lc = q.left.add_vertex("C•"), ln = q.left.add_vertex("N•");
  const std::size_t rc = q.right.add_vertex("C"), rn = q.right.add_vertex("N");
  q.right.add_edge(rc, rn, "amide");
  q.right.add_vertex("H2O");
  q.kept_vertices = {{lc, rc}, {ln, rn}};
  q.condition = [](const rw::Graph& H, const rw::Match& m) {
    const auto c = H.components();
    return c[m.vertex[0]] != c[m.vertex[1]];
  };
  q.validate();
  return q;
}

inline rw::Graph glycines(std::size_t n, bool labelled) {
  rw::Graph G;
  G.directed = true;
  for (std::size_t i = 0; i < n; i++) {
    const std::size_t u = G.add_vertex(labelled ? "Gly" + std::to_string(i + 1) : "Gly");
    G.add_edge(G.add_vertex("N•"), u, "amine");
    G.add_edge(u, G.add_vertex("C•"), "carboxyl");
  }
  return G;
}

inline Json glycine(std::size_t n = 3) {
  const std::vector<rw::Rule> rules{coupling_rule()};
  const auto terminal = [](const rw::Graph& G) {
    const auto c = G.components();
    std::set<std::size_t> chains;
    for (std::size_t v = 0; v < G.size(); v++)
      if (G.vertex[v] != "H2O") chains.insert(c[v]);
    return chains.size() == 1;
  };
  Json out = Json::object({{"title", "Labelled glycine backbones (Figs. grammar-growth, grammar-derivation-tree)"}, {"units", static_cast<double>(n)}});
  Json rj = Json::array();
  rj.push(rule_json(rules[0]));
  out.set("rules", rj);
  for (bool labelled : {true, false}) {
    const auto D = rw::explore(rules, glycines(n, labelled), n - 1);
    Json part = Json::object({{"space", space_json(D, terminal)}});
    // Terminal classes, their orbits under S_n (labelled), equivariance at every class.
    std::vector<rw::Graph> ends;
    for (std::size_t c = 0; c < D.graphs.size(); c++)
      if (terminal(D.graphs[c])) ends.push_back(D.graphs[c]);
    unsigned long long derivations = 0;
    for (std::size_t c = 0; c < D.graphs.size(); c++)
      if (terminal(D.graphs[c])) derivations += D.derivations_to(c);
    part.set("terminalClasses", static_cast<double>(ends.size()));
    part.set("terminalDerivations", static_cast<double>(derivations));
    if (labelled) {
      std::vector<rw::Relabelling> gens{{{"Gly1", "Gly2"}, {"Gly2", "Gly1"}}, {}};
      for (std::size_t i = 0; i < n; i++) gens[1]["Gly" + std::to_string(i + 1)] = "Gly" + std::to_string((i + 1) % n + 1);
      const auto orb = rw::orbits(ends, gens);
      part.set("orbits", static_cast<double>(*std::max_element(orb.begin(), orb.end()) + 1));
      bool eq = true;
      for (std::size_t c = 0; c < D.graphs.size(); c++)
        for (const auto& g : gens) eq = eq && rw::equivariant(rules, D.graphs[c], g);
      part.set("equivariant", eq);
    }
    // Water created along every derivation of a full chain: n − 1.
    std::size_t water = 0;
    for (const auto& G : ends) water = std::max<std::size_t>(water, static_cast<std::size_t>(std::count(G.vertex.begin(), G.vertex.end(), "H2O")));
    part.set("water", static_cast<double>(water));
    out.set(labelled ? "labelled" : "identical", part);
  }
  // The net reaction n Gly → Gly_n + (n − 1) H2O balanced by element, with the species' formulas
  // (glycine C2H5NO2; the chain C(2n) H(3n+2) N(n) O(n+1)).
  physis::chemistry::Network net;
  net.species = {"Gly", "chain", "H2O"};
  net.add("condense", {{"Gly", static_cast<unsigned>(n)}}, {{"chain", 1}, {"H2O", static_cast<unsigned>(n - 1)}}, 1);
  morphe::algebra::Matrix A(4, 3);  // C, H, N, O
  const long long N = static_cast<long long>(n);
  const long long gly[4] = {2, 5, 1, 2}, chain[4] = {2 * N, 3 * N + 2, N, N + 1}, water[4] = {0, 2, 0, 1};
  for (int e = 0; e < 4; e++) A(e, 0) = gly[e], A(e, 1) = chain[e], A(e, 2) = water[e];
  out.set("balanced", physis::chemistry::balance(net, A, {0, 0, 0}).balanced());
  return out;
}

// ---- the tetrahedral stereocentre (Case Study I) --------------------------------------------------------------------------------

// The start graph: the nonterminal A_tet joined by ports 1..4 to the ligand atoms H, F, Cl, Br.
// p± relabel A_tet to the carbon with its orientation label, and each port to a single bond.
inline rw::Rule tet_rule(int sign) {
  rw::Rule q;
  q.name = sign > 0 ? "p+: A_tet(u1..u4) → C(τ = +1)" : "p−: A_tet(u1..u4) → C(τ = −1)";
  const std::size_t la = q.left.add_vertex("A_tet"), ra = q.right.add_vertex(sign > 0 ? "C@+" : "C@-");
  q.kept_vertices.push_back({la, ra});
  for (int i = 1; i <= 4; i++) {
    const std::size_t lu = q.left.add_vertex(rw::any_label), ru = q.right.add_vertex(rw::any_label);
    q.left.add_edge(la, lu, "port" + std::to_string(i));
    q.right.add_edge(ra, ru, "single");
    q.kept_vertices.push_back({lu, ru});
    q.kept_edges.push_back({static_cast<std::size_t>(i - 1), static_cast<std::size_t>(i - 1)});
  }
  q.validate();
  return q;
}

inline Json stereocentre() {
  const std::vector<rw::Rule> rules{tet_rule(+1), tet_rule(-1)};
  rw::Graph S;
  const std::size_t a = S.add_vertex("A_tet");
  const std::vector<std::string> ligands{"H", "F", "Cl", "Br"};
  for (int i = 0; i < 4; i++) S.add_edge(a, S.add_vertex(ligands[i]), "port" + std::to_string(i + 1));
  const auto D = rw::explore(rules, S, 2);
  const auto terminal = [](const rw::Graph& G) { return !has_label_prefix(G, "A_tet"); };
  const rw::Relabelling p{{"C@+", "C@-"}, {"C@-", "C@+"}};
  std::vector<rw::Graph> ends;
  for (std::size_t c = 0; c < D.graphs.size(); c++)
    if (terminal(D.graphs[c])) ends.push_back(D.graphs[c]);
  const auto orb = rw::orbits(ends, {p});
  // The constitutional projection T_conn forgets the orientation label.
  const rw::Relabelling forget{{"C@+", "C"}, {"C@-", "C"}};
  std::set<std::string> constitutional;
  for (const auto& G : ends) {
    rw::Graph H = G;
    for (auto& l : H.vertex) l = (l == "C@+" || l == "C@-") ? "C" : l;
    constitutional.insert(rw::code(H));
  }
  // Geometry: the two configurations, their signs, descriptors and congruences.
  const Molecule plus = chbrclf(+1), minus = chbrclf(-1);
  const std::array<std::size_t, 4> lig{1, 2, 3, 4};
  const auto between = physis::symmetry::congruences(plus.points(), minus.points());
  const auto self = physis::symmetry::point_group(plus.points());
  Json geometry = Json::object({{"plus", molecule_json(plus)}, {"minus", molecule_json(minus)}, {"tauPlus", tau(plus, lig)}, {"tauMinus", tau(minus, lig)},
                                {"cipPlus", cip(plus, lig)}, {"cipMinus", cip(minus, lig)},
                                {"congruences", static_cast<double>(between.size())},
                                {"congruenceDet", between.empty() ? 0 : between[0].det()}, {"pointGroup", self.schoenflies}});
  // An odd relabelling of the ligands flips τ as a reflection does.
  Molecule swapped = plus;
  std::swap(swapped.atoms[2].element, swapped.atoms[3].element);  // F ↔ Cl at their places
  geometry.set("tauAfterSwap", tau(swapped, {1, 3, 2, 4}));
  Json rj = Json::array();
  for (const auto& q : rules) rj.push(rule_json(q));
  bool closed = rw::closed_under(rules, p), eq = true;
  for (std::size_t c = 0; c < D.graphs.size(); c++) eq = eq && rw::equivariant(rules, D.graphs[c], p);
  const std::vector<rw::Rule> half{tet_rule(+1)};
  return Json::object({{"title", "Bromochlorofluoromethane: a group-equivariant grammar (Case Study I)"}, {"rules", rj},
                       {"space", space_json(D, terminal)}, {"terminalClasses", static_cast<double>(ends.size())},
                       {"orbits", static_cast<double>(orb.empty() ? 0 : *std::max_element(orb.begin(), orb.end()) + 1)},
                       {"constitutionalClasses", static_cast<double>(constitutional.size())}, {"closed", closed}, {"equivariant", eq},
                       {"halfClosed", rw::closed_under(half, p)}, {"halfEquivariant", rw::equivariant(half, S, p)}, {"geometry", geometry}});
}

// ---- acid–base edit (Fig. proton-graph-edit) ---------------------------------------------------------------------------------

// Explicit atoms; charges in the labels (O-, O+). The rule matches the acid's C(=O)–O–H* and a
// water O with its two H, deletes O–H*, creates O(water)–H*, and relabels the two oxygens: every
// atom kept, its identity traced.
inline rw::Rule proton_rule() {
  rw::Rule q;
  q.name = "q: CH3C(=O)O–H* + H2O → CH3C(=O)O⁻ + H2O–H*⁺";
  auto& L = q.left;
  auto& R = q.right;
  const std::size_t c = L.add_vertex("C"), oc = L.add_vertex("O"), oa = L.add_vertex("O"), h = L.add_vertex("H"), ow = L.add_vertex("O"),
                    w1 = L.add_vertex("H"), w2 = L.add_vertex("H");
  L.add_edge(c, oc, "double"), L.add_edge(c, oa, "single"), L.add_edge(oa, h, "single"), L.add_edge(ow, w1, "single"), L.add_edge(ow, w2, "single");
  const std::size_t rc = R.add_vertex("C"), roc = R.add_vertex("O"), roa = R.add_vertex("O-"), rh = R.add_vertex("H"), row = R.add_vertex("O+"),
                    rw1 = R.add_vertex("H"), rw2 = R.add_vertex("H");
  R.add_edge(rc, roc, "double"), R.add_edge(rc, roa, "single"), R.add_edge(row, rw1, "single"), R.add_edge(row, rw2, "single"), R.add_edge(row, rh, "single");
  q.kept_vertices = {{c, rc}, {oc, roc}, {oa, roa}, {h, rh}, {ow, row}, {w1, rw1}, {w2, rw2}};
  q.kept_edges = {{0, 0}, {1, 1}, {3, 2}, {4, 3}};
  q.validate();
  return q;
}

inline rw::Graph acid_and_water() {
  rw::Graph G;
  const std::size_t cm = G.add_vertex("C"), c = G.add_vertex("C"), oc = G.add_vertex("O"), oa = G.add_vertex("O"), h = G.add_vertex("H");
  for (int i = 0; i < 3; i++) G.add_edge(cm, G.add_vertex("H"), "single");
  G.add_edge(cm, c, "single"), G.add_edge(c, oc, "double"), G.add_edge(c, oa, "single"), G.add_edge(oa, h, "single");
  const std::size_t ow = G.add_vertex("O");
  G.add_edge(ow, G.add_vertex("H"), "single"), G.add_edge(ow, G.add_vertex("H"), "single");
  return G;
}

inline std::map<std::string, int> ledger(const rw::Graph& G) {
  std::map<std::string, int> count;
  int charge = 0;
  for (const auto& l : G.vertex) {
    std::string e = l;
    if (!e.empty() && (e.back() == '+' || e.back() == '-')) charge += e.back() == '+' ? 1 : -1, e.pop_back();
    count[e]++;
  }
  count["charge"] = charge;
  return count;
}

inline Json proton_transfer() {
  const auto q = proton_rule();
  const auto H = acid_and_water();
  const auto ds = rw::derivations({q}, H);
  std::set<std::string> classes;
  for (const auto& d : ds) classes.insert(rw::code(d.result));
  Json before = Json::object(), after = Json::object();
  for (const auto& [k, v] : ledger(H)) before.set(k, v);
  if (!ds.empty())
    for (const auto& [k, v] : ledger(ds[0].result)) after.set(k, v);
  Json trace = Json::array();
  if (!ds.empty())
    for (auto t : ds[0].trace) trace.push(static_cast<double>(t));
  // The species-level balance: CH3COOH + H2O → CH3COO⁻ + H3O⁺, AN = 0 and zᵀN = 0.
  physis::chemistry::Network net;
  net.species = {"CH3COOH", "H2O", "CH3COO-", "H3O+"};
  net.add("transfer", {{"CH3COOH", 1}, {"H2O", 1}}, {{"CH3COO-", 1}, {"H3O+", 1}}, 1);
  morphe::algebra::Matrix A(3, 4);  // C, H, O
  const long long comp[3][4] = {{2, 0, 2, 0}, {4, 2, 3, 3}, {2, 1, 2, 1}};
  for (int e = 0; e < 3; e++)
    for (int s = 0; s < 4; s++) A(e, s) = comp[e][s];
  return Json::object({{"title", "Acid–base edit with atom and charge accounting (Fig. proton-graph-edit)"}, {"rules", Json::array({rule_json(q)})},
                       {"host", graph_json(H)}, {"result", ds.empty() ? Json() : graph_json(ds[0].result)}, {"trace", trace},
                       {"matches", static_cast<double>(ds.size())}, {"automorphismsOfL", static_cast<double>(rw::automorphisms(q.left))},
                       {"productClasses", static_cast<double>(classes.size())}, {"before", before}, {"after", after},
                       {"balanced", physis::chemistry::balance(net, A, {0, 0, -1, 1}).balanced()},
                       {"labelPreservingSpanPossible", false}});
}

}  // namespace chem
