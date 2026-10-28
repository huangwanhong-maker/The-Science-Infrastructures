// The chemistry the plug-in needs, provisionally here until a chemistry domain holds it (the user's
// decision of 2026-10-28): a few elements, idealized geometries of the example molecules, the
// orientation sign of a tetrahedral centre and its CIP descriptor, and their point groups through
// physis's symmetry.
//
// Geometries are idealized: VSEPR angles (tetrahedral 109.47°, trigonal 120°, octahedral 90°) and
// bond lengths taken as sums of covalent radii (Cordero et al. 2008, Dalton Trans. 2832, Table 2),
// not measured structures. They fix the symmetry; their lengths only draw the molecules.

#pragma once

#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <physis/symmetry/point_group.hpp>
// The exact character tables come from morphe's characters (its decision 0006) where they are present.
#if __has_include(<morphe/algebra/characters.hpp>)
#include <morphe/algebra/characters.hpp>
#define CHEM_CHARACTER_TABLES 1
#endif
#include <science/json.hpp>

namespace chem {

using physis::symmetry::Vec3;
using science::Json;

struct Element {
  std::string symbol;
  int number;     // Z, which also orders CIP priorities among single atoms
  double radius;  // covalent radius, Å (Cordero et al. 2008)
};

inline const std::vector<Element>& elements() {
  static const std::vector<Element> table{
      {"H", 1, 0.31}, {"B", 5, 0.84}, {"C", 6, 0.76}, {"N", 7, 0.71}, {"O", 8, 0.66}, {"F", 9, 0.57},
      {"S", 16, 1.05}, {"Cl", 17, 1.02}, {"Br", 35, 1.20},
  };
  return table;
}
inline const Element& element(const std::string& s) {
  for (const auto& e : elements())
    if (e.symbol == s) return e;
  throw std::invalid_argument("chemistry: no element " + s);
}
inline double bond_length(const std::string& a, const std::string& b) { return element(a).radius + element(b).radius; }

struct Atom {
  std::string element;
  Vec3 x{};
  int charge = 0;
};
struct Bond {
  std::size_t a = 0, b = 0;
  int order = 1;
};
struct Molecule {
  std::string id, name, formula, note;
  std::vector<Atom> atoms;
  std::vector<Bond> bonds;

  std::size_t add(const std::string& e, Vec3 x) {
    atoms.push_back({e, x, 0});
    return atoms.size() - 1;
  }
  void bond(std::size_t a, std::size_t b, int order = 1) { bonds.push_back({a, b, order}); }
  physis::symmetry::LabelledPoints points() const {
    physis::symmetry::LabelledPoints P;
    for (const auto& a : atoms) P.add(a.element, a.x);
    return P;
  }
};

inline Vec3 add(Vec3 a, Vec3 b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
inline Vec3 scaled(double s, Vec3 a) { return {s * a[0], s * a[1], s * a[2]}; }
inline Vec3 unit(Vec3 a) {
  const double n = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  return {a[0] / n, a[1] / n, a[2] / n};
}
// The atom of element e bonded to atom `from` in direction d, at the bond length.
inline std::size_t substituent(Molecule& m, std::size_t from, const std::string& e, Vec3 d, int order = 1) {
  const std::size_t i = m.add(e, add(m.atoms[from].x, scaled(bond_length(m.atoms[from].element, e), unit(d))));
  m.bond(from, i, order);
  return i;
}

// The four tetrahedral directions.
inline std::array<Vec3, 4> tetrahedron() { return {{{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}}}; }

inline double det3(Vec3 a, Vec3 b, Vec3 c) {
  return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
}
inline Vec3 minus(Vec3 a, Vec3 b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }

// The manuscript's orientation sign, Eq. (case-chiral-sign): τ = sgn det[r_F − r_H, r_Cl − r_H, r_Br − r_H]
// for a centre bearing H, F, Cl and Br (the ligands in this order, by index into the molecule).
inline int tau(const Molecule& m, const std::array<std::size_t, 4>& ligands) {
  const Vec3 h = m.atoms[ligands[0]].x;
  const double d = det3(minus(m.atoms[ligands[1]].x, h), minus(m.atoms[ligands[2]].x, h), minus(m.atoms[ligands[3]].x, h));
  return d > 0 ? 1 : d < 0 ? -1 : 0;
}

// The CIP descriptor of a centre whose four ligands are single atoms of distinct atomic numbers
// (the first sphere decides; IUPAC 2013 Recommendations P-92.1.4, Rule 1a): R when, seen with the
// lowest-priority ligand d pointing away, the others a > b > c run clockwise, which is
// det[a − d, b − d, c − d] < 0.
inline std::string cip(const Molecule& m, std::array<std::size_t, 4> ligands) {
  std::sort(ligands.begin(), ligands.end(), [&](std::size_t i, std::size_t j) { return element(m.atoms[i].element).number > element(m.atoms[j].element).number; });
  for (int k = 0; k + 1 < 4; k++)
    if (m.atoms[ligands[k]].element == m.atoms[ligands[k + 1]].element) return "";
  const Vec3 d = m.atoms[ligands[3]].x;
  const double s = det3(minus(m.atoms[ligands[0]].x, d), minus(m.atoms[ligands[1]].x, d), minus(m.atoms[ligands[2]].x, d));
  return s < 0 ? "R" : "S";
}

// CHBrClF with the ligands H, F, Cl, Br placed so that τ = sign.
inline Molecule chbrclf(int sign) {
  Molecule m;
  m.id = sign > 0 ? "chbrclf-plus" : "chbrclf-minus";
  m.formula = "CHBrClF";
  const std::size_t c = m.add("C", {0, 0, 0});
  const auto t = tetrahedron();
  std::array<std::size_t, 4> lig{};
  const std::array<std::string, 4> names{"H", "F", "Cl", "Br"};
  for (int i = 0; i < 4; i++) lig[i] = substituent(m, c, names[i], t[i]);
  if (tau(m, lig) != sign) {
    // Reflect in the plane x = 0: the improper map that carries one configuration to the other.
    for (auto& a : m.atoms) a.x[0] = -a.x[0];
  }
  m.name = std::string("bromochlorofluoromethane, X") + (sign > 0 ? "+" : "−") + " (" + cip(m, lig) + ")";
  return m;
}

// The idealized molecules of the symmetry view.
inline std::vector<Molecule> catalogue() {
  std::vector<Molecule> out;
  const double pi = std::acos(-1.0);
  const auto tet = tetrahedron();
  {
    out.push_back(chbrclf(+1));
    out.push_back(chbrclf(-1));
  }
  for (bool cis : {true, false}) {
    Molecule m;
    m.id = cis ? "dichloroethene-z" : "dichloroethene-e";
    m.name = cis ? "(Z)-1,2-dichloroethene" : "(E)-1,2-dichloroethene";
    m.formula = "C2H2Cl2";
    const double half = bond_length("C", "C") / 2 * 0.87;  // a double bond, about 13% shorter than a single
    const std::size_t c1 = m.add("C", {-half, 0, 0}), c2 = m.add("C", {half, 0, 0});
    m.bond(c1, c2, 2);
    const double a = 2 * pi / 3;
    substituent(m, c1, "Cl", {std::cos(a), std::sin(a), 0});
    substituent(m, c1, "H", {std::cos(a), -std::sin(a), 0});
    substituent(m, c2, cis ? "Cl" : "H", {-std::cos(a), std::sin(a), 0});
    substituent(m, c2, cis ? "H" : "Cl", {-std::cos(a), -std::sin(a), 0});
    m.note = cis ? "Both chlorines on one side of the C=C bond." : "The chlorines on opposite sides: centrosymmetric.";
    out.push_back(m);
  }
  {
    Molecule m;
    m.id = "benzene";
    m.name = "benzene";
    m.formula = "C6H6";
    const double r = 1.39;
    std::vector<std::size_t> c;
    for (int k = 0; k < 6; k++) c.push_back(m.add("C", {r * std::cos(k * pi / 3), r * std::sin(k * pi / 3), 0}));
    for (int k = 0; k < 6; k++) m.bond(c[k], c[(k + 1) % 6], k % 2 ? 1 : 2);
    for (int k = 0; k < 6; k++) substituent(m, c[k], "H", {std::cos(k * pi / 3), std::sin(k * pi / 3), 0});
    m.note = "Drawn with alternating bonds (one Kekulé structure); the symmetry is that of the nuclei.";
    out.push_back(m);
  }
  {
    Molecule m;
    m.id = "water", m.name = "water", m.formula = "H2O";
    const std::size_t o = m.add("O", {0, 0, 0});
    const double a = 104.5 / 2 * pi / 180;
    substituent(m, o, "H", {std::sin(a), 0, std::cos(a)});
    substituent(m, o, "H", {-std::sin(a), 0, std::cos(a)});
    out.push_back(m);
  }
  {
    Molecule m;
    m.id = "ammonia", m.name = "ammonia", m.formula = "NH3";
    const std::size_t n = m.add("N", {0, 0, 0});
    for (int k = 1; k < 4; k++) substituent(m, n, "H", tet[k]);  // the lone pair along tet[0]
    out.push_back(m);
  }
  {
    Molecule m;
    m.id = "methane", m.name = "methane", m.formula = "CH4";
    const std::size_t c = m.add("C", {0, 0, 0});
    for (const auto& d : tet) substituent(m, c, "H", d);
    out.push_back(m);
  }
  {
    Molecule m;
    m.id = "ethylene", m.name = "ethylene", m.formula = "C2H4";
    const double half = bond_length("C", "C") / 2 * 0.87;
    const std::size_t c1 = m.add("C", {-half, 0, 0}), c2 = m.add("C", {half, 0, 0});
    m.bond(c1, c2, 2);
    const double a = 2 * pi / 3;
    for (double s : {1.0, -1.0}) {
      substituent(m, c1, "H", {std::cos(a), s * std::sin(a), 0});
      substituent(m, c2, "H", {-std::cos(a), s * std::sin(a), 0});
    }
    out.push_back(m);
  }
  {
    Molecule m;
    m.id = "boron-trifluoride", m.name = "boron trifluoride", m.formula = "BF3";
    const std::size_t b = m.add("B", {0, 0, 0});
    for (int k = 0; k < 3; k++) substituent(m, b, "F", {std::cos(2 * pi * k / 3), std::sin(2 * pi * k / 3), 0});
    out.push_back(m);
  }
  {
    Molecule m;
    m.id = "sulfur-hexafluoride", m.name = "sulfur hexafluoride", m.formula = "SF6";
    const std::size_t s = m.add("S", {0, 0, 0});
    for (Vec3 d : std::vector<Vec3>{{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}}) substituent(m, s, "F", d);
    out.push_back(m);
  }
  {
    // Staggered ethane: the methyls about the C–C axis (z), turned 60° from each other.
    Molecule m;
    m.id = "ethane", m.name = "ethane (staggered)", m.formula = "C2H6";
    const double half = bond_length("C", "C") / 2;
    const std::size_t c1 = m.add("C", {0, 0, -half}), c2 = m.add("C", {0, 0, half});
    m.bond(c1, c2);
    const double s = std::sin(std::acos(-1.0 / 3)), z = 1.0 / 3;
    for (int k = 0; k < 3; k++) {
      const double a = 2 * pi * k / 3;
      substituent(m, c1, "H", {s * std::cos(a), s * std::sin(a), -z});
      substituent(m, c2, "H", {s * std::cos(a + pi / 3), s * std::sin(a + pi / 3), z});
    }
    out.push_back(m);
  }
  {
    // Hydrogen peroxide, skew: dihedral angle 111.5° (H–O–O–H), O–O–H 95°.
    Molecule m;
    m.id = "hydrogen-peroxide", m.name = "hydrogen peroxide (skew)", m.formula = "H2O2";
    const double half = bond_length("O", "O") / 2;
    const std::size_t o1 = m.add("O", {-half, 0, 0}), o2 = m.add("O", {half, 0, 0});
    m.bond(o1, o2);
    const double tilt = (180 - 95) * pi / 180, phi = 111.5 / 2 * pi / 180;
    substituent(m, o1, "H", {-std::cos(tilt), std::sin(tilt) * std::cos(phi), std::sin(tilt) * std::sin(phi)});
    substituent(m, o2, "H", {std::cos(tilt), std::sin(tilt) * std::cos(phi), -std::sin(tilt) * std::sin(phi)});
    m.note = "Chiral in any fixed skew conformation; the two mirror forms interconvert by rotation about O–O.";
    out.push_back(m);
  }
  {
    // Allene: the two CH2 planes perpendicular.
    Molecule m;
    m.id = "allene", m.name = "allene", m.formula = "C3H4";
    const double l = bond_length("C", "C") * 0.87;
    const std::size_t c0 = m.add("C", {0, 0, 0}), c1 = m.add("C", {0, 0, -l}), c2 = m.add("C", {0, 0, l});
    m.bond(c0, c1, 2), m.bond(c0, c2, 2);
    const double a = pi / 3;
    substituent(m, c1, "H", {std::sin(a), 0, -std::cos(a)});
    substituent(m, c1, "H", {-std::sin(a), 0, -std::cos(a)});
    substituent(m, c2, "H", {0, std::sin(a), std::cos(a)});
    substituent(m, c2, "H", {0, -std::sin(a), std::cos(a)});
    out.push_back(m);
  }
  {
    Molecule m;
    m.id = "carbon-dioxide", m.name = "carbon dioxide", m.formula = "CO2";
    const std::size_t c = m.add("C", {0, 0, 0});
    substituent(m, c, "O", {0, 0, 1}, 2), substituent(m, c, "O", {0, 0, -1}, 2);
    out.push_back(m);
  }
  for (auto& m : out)
    if (m.name.empty()) m.name = m.id;
  return out;
}

// ---- JSON ----------------------------------------------------------------------------------------------------------------------

inline Json vec(const Vec3& v) { return Json::array({v[0], v[1], v[2]}); }

inline Json molecule_json(const Molecule& m) {
  Json atoms = Json::array(), bonds = Json::array();
  for (const auto& a : m.atoms) atoms.push(Json::object({{"element", a.element}, {"x", vec(a.x)}, {"charge", a.charge}, {"radius", element(a.element).radius}}));
  for (const auto& b : m.bonds) bonds.push(Json::array({static_cast<double>(b.a), static_cast<double>(b.b), b.order}));
  return Json::object({{"id", m.id}, {"name", m.name}, {"formula", m.formula}, {"note", m.note}, {"atoms", atoms}, {"bonds", bonds}});
}

// The geometric kind of an isometry: E, Cn^k, σ, i, Sn^k.
inline std::string operation_name(const physis::symmetry::Isometry& g) {
  using namespace physis::symmetry::pg_detail;
  const double pi = std::acos(-1.0);
  if (g.proper()) {
    const double theta = axis_angle(g.Q).second;
    if (theta < 1e-6) return "E";
    // θ = 2πk/n with k/n in lowest terms.
    for (int n = 2; n <= 12; n++)
      for (int k = 1; k < n; k++)
        if (std::fabs(theta - 2 * pi * k / n) < 1e-6 && std::gcd(k, n) == 1) return "C" + std::to_string(n) + (k > 1 ? "^" + std::to_string(k) : "");
    return "C?";
  }
  if (std::fabs(g.trace() + 3) < 1e-6) return "i";
  if (std::fabs(g.trace() - 1) < 1e-6) return "σ";
  const int k = element_order(g.Q, 1e-6);
  return "S" + std::to_string(k);
}

// The point group of a molecule: its symbol, order, elements for drawing, the responses it allows,
// and its character table (morphe, exactly).
inline Json point_group_json(const Molecule& m) {
  const auto G = physis::symmetry::point_group(m.points(), 1e-6);
  Json out = Json::object({{"schoenflies", G.schoenflies}, {"linear", G.linear}, {"order", static_cast<double>(G.order())}});
  out.set("centroid", vec(G.centroid));
  if (G.linear) {
    out.set("chiral", false);
    out.set("dipole", G.schoenflies == "C∞v" ? 1 : 0);
    out.set("axial", 1);
    return out;
  }
  out.set("chiral", G.chiral());
  out.set("dipole", static_cast<double>(G.fixed_vectors(false)));
  out.set("axial", static_cast<double>(G.fixed_vectors(true)));
  Json els = Json::array();
  for (const auto& e : G.symmetry_elements) {
    const char* kind = e.kind == physis::symmetry::SymmetryElement::Kind::rotation ? "rotation"
                       : e.kind == physis::symmetry::SymmetryElement::Kind::rotoreflection ? "rotoreflection"
                       : e.kind == physis::symmetry::SymmetryElement::Kind::mirror ? "mirror" : "inversion";
    els.push(Json::object({{"kind", kind}, {"order", e.order}, {"axis", vec(e.axis)}}));
  }
  out.set("elements", els);
#ifdef CHEM_CHARACTER_TABLES
  // The character table, its classes named by their geometric kind.
  const auto T = morphe::algebra::character_table(G.group);
  Json classes = Json::array(), rows = Json::array();
  for (std::size_t j = 0; j < T.classes.size(); j++)
    classes.push(Json::object({{"size", static_cast<double>(T.class_sizes[j])}, {"operation", operation_name(G.elements[T.representatives[j]])}}));
  for (std::size_t chi = 0; chi < T.size(); chi++) {
    Json row = Json::array();
    for (std::size_t j = 0; j < T.classes.size(); j++) {
      const auto z = T(chi, j).to_complex();
      row.push(Json::object({{"exact", T(chi, j).str()}, {"re", z.real()}, {"im", z.imag()}}));
    }
    rows.push(row);
  }
  // The vector and axial representations' characters on the classes.
  Json vecchi = Json::array(), axchi = Json::array();
  for (std::size_t j = 0; j < T.classes.size(); j++) {
    const auto& g = G.elements[T.representatives[j]];
    vecchi.push(g.trace());
    axchi.push(g.det() * g.trace());
  }
  out.set("characters", Json::object({{"classes", classes}, {"rows", rows}, {"vector", vecchi}, {"axial", axchi}}));
#endif
  return out;
}

}  // namespace chem
