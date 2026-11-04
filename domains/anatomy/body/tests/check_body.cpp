// Checks of the body read from the data library's product (body.hpp), skipped when it has not been built.
//
// - The product holds BodyParts3D 4.0 whole: 2,234 elements, 4,209,773 points and 6,681,030 triangles; the
//   is-a relation one tree under "anatomical entity", the part-of relation one tree under "human body".
// - Every rule of the systems' table names a concept of the release, and every element falls under a rule:
//   no part of the body is left without a system. Named parts fall where anatomy puts them.
// - The woman of the Human Reference Atlas, when her product is built: her counts, one part-of tree, a system for
//   every element by the Atlas's layers and names, her height, and her organs on their sides in LPS.
// - BodyParts3D's trees are spanning trees of FMA's relations: for every compound, the elements whose
//   concepts lie below it in the tree are among those its table lists. The lists go beyond the tree where
//   FMA's part-of graph gives a concept other parents (the body wall's muscles, also parts of their regions),
//   and in is-a where an element without a concept of its own is listed under two incomparable compounds.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <set>
#include <tuple>
#include <string>
#include <vector>

#include <anatomy/body/body.hpp>
#include <calx/verification/check.hpp>

using calx::verification::Suite;
namespace ab = anatomy::body;
namespace ap = anatomy::parts;

// The woman of the Human Reference Atlas (body-hra-f; data/products/body_female.py), when built.
void female(Suite& s, const std::filesystem::path& dir) {
  if (!std::filesystem::exists(dir / "body-hra-f-elements.tsv")) {
    s.skip("the female body: the product is not built (python data/acquisition/fetch.py hra-female; python data/products/body_female.py)");
    return;
  }
  const auto B = ab::read_body(dir, "body-hra-f");
  const auto& A = B.atlas;
  std::size_t triangles = 0;
  for (const auto& e : B.elements) triangles += e.triangles;
  s.require("the woman of the Human Reference Atlas, v1.8: 911 elements, 7,060,019 triangles", B.elements.size() == 911 && triangles == 7060019,
            std::to_string(B.elements.size()) + " elements, " + std::to_string(triangles) + " triangles");
  const auto rp = A.partof.roots();
  s.require("her part-of tree is one tree under the united model, read by the Atlas's rules", rp.size() == 1 && A.concepts[rp[0]].id == "VH_F" && ap::is_reference_atlas(A));
  std::map<std::string, std::size_t> count;
  std::size_t without = 0;
  for (std::size_t e = 0; e < B.elements.size(); e++) {
    const auto k = ap::system_of_element(A, e);
    if (k == ap::none) without++;
    else count[std::string(ap::SYSTEMS[k].key)]++;
  }
  std::string tally;
  for (const auto& sys : ap::SYSTEMS) tally += std::string(sys.key) + " " + std::to_string(count[std::string(sys.key)]) + ", ";
  std::printf("      the woman: %s\n", tally.c_str());
  s.require("every element of hers has a system", without == 0, std::to_string(without) + " without");
  s.require("her vessels told apart by name: 55 arteries and 55 veins, the heart's chambers and valves the heart",
            count["arteries"] == 55 && count["veins"] == 55 && count["heart"] == 14);
  // an organ's centre: the middle of the box of the elements below its node
  const auto centre = [&](const std::string& id) {
    const std::size_t c = A.find(id);
    std::array<double, 3> lo{INFINITY, INFINITY, INFINITY}, hi{-INFINITY, -INFINITY, -INFINITY};
    for (const auto& e : B.elements)
      if (c != ap::none && A.forest(ap::Relation::partof).below(e.concept_, c))
        for (int a = 0; a < 3; a++) lo[a] = std::min(lo[a], e.lower[a]), hi[a] = std::max(hi[a], e.upper[a]);
    return std::array<double, 3>{(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2};
  };
  std::printf("      spleen x %.0f mm, liver x %.0f mm, kidneys y %.0f mm against the liver's %.0f mm\n", centre("VH_F_spleen")[0], centre("VH_F_liver")[0],
              centre("VH_F_kidney")[1], centre("VH_F_liver")[1]);
  double height = 0;
  for (const auto& e : B.elements)
    if (e.id == "VH_F_skin") height = e.upper[2] - e.lower[2];
  s.check("her height, soles to crown, the Atlas's 1,657.8 mm", height - 1657.8, 0.5, "mm");
  s.require("in LPS: the spleen on her left (+x), the liver on her right, the kidneys behind the liver (+y)",
            centre("VH_F_spleen")[0] > 0 && centre("VH_F_liver")[0] < 0 && centre("VH_F_kidney")[1] > centre("VH_F_liver")[1]);
}

int main() {
  Suite s;
  const std::filesystem::path dir = ANATOMY_PRODUCTS;
  if (!std::filesystem::exists(dir / "body-bp3d-elements.tsv")) {
    s.skip("the body: the product is not built (python data/acquisition/fetch.py bodyparts3d; python data/products/body_anatomy.py)");
    return s.done();
  }
  const auto B = ab::read_body(dir);
  const auto& A = B.atlas;
  {
    std::size_t points = 0, triangles = 0;
    for (const auto& e : B.elements) points += e.vertices, triangles += e.triangles;
    s.require("BodyParts3D 4.0 whole: 2,234 elements, 4,209,773 points, 6,681,030 triangles",
              B.elements.size() == 2234 && points == 4209773 && triangles == 6681030, std::to_string(B.elements.size()) + " elements");
    const auto ri = A.isa.roots(), rp = A.partof.roots();
    s.require("is-a is one tree under the anatomical entity", ri.size() == 1 && A.concepts[ri[0]].name == "anatomical entity");
    s.require("part-of is one tree under the human body", rp.size() == 1 && A.concepts[rp[0]].name == "human body");
    const auto S = B.surface(0);
    s.require("an element's surface: its points and triangles as counted", S.points.size() == B.elements[0].vertices && S.triangles.size() == B.elements[0].triangles);
  }
  {
    s.require("every rule of the systems' table names a concept of this release", ap::missing_rules(A).empty());
    std::map<std::string, std::size_t> count;
    std::size_t without = 0;
    for (std::size_t e = 0; e < B.elements.size(); e++) {
      const auto k = ap::system_of_element(A, e);
      if (k == ap::none) without++;
      else count[std::string(ap::SYSTEMS[k].key)]++;
    }
    s.require("every element falls under a rule: no part without a system", without == 0, std::to_string(without) + " without");
    std::string tally;
    for (const auto& sys : ap::SYSTEMS) tally += std::string(sys.key) + " " + std::to_string(count[std::string(sys.key)]) + ", ";
    std::printf("      %s\n", tally.c_str());
    const std::map<std::string, std::string> expected{
        {"Left kidney", "urinary"}, {"Right great saphenous vein", "veins"}, {"Left superior pulmonary vein", "veins"},
        {"Manubrium", "skeleton"}, {"Cavity of left ventricle", "heart"}, {"Anterior leaflet of mitral valve", "heart"},
        {"Trachea", "respiratory"}, {"Thyroid cartilage", "respiratory"}, {"Left lacrimal gland", "eye"}, {"Spleen", "lymphoid"},
        {"Left testis", "genital"}, {"Pituitary gland", "endocrine"}, {"Skin", "skin"}, {"Corpus callosum", "brain"},
        {"Diaphragm", "muscles"}, {"Left scapula", "skeleton"}, {"Inferior vena cava", "veins"}, {"Basilar artery", "arteries"},
        {"Left optic nerve", "nerves"}, {"Esophagus", "digestive"}, {"Left iliotibial tract", "fascia"}, {"External ear", "ear"}};
    std::string wrong;
    std::size_t found = 0;
    for (std::size_t e = 0; e < B.elements.size(); e++) {
      const auto it = expected.find(B.elements[e].name);
      if (it == expected.end()) continue;
      found++;
      const auto k = ap::system_of_element(A, e);
      if (k == ap::none || ap::SYSTEMS[k].key != it->second) wrong += B.elements[e].name + "; ";
    }
    s.require("named parts fall where anatomy puts them (a kidney urinary, a pulmonary vein a vein, the mitral valve the heart, ...)",
              found >= expected.size() && wrong.empty(), wrong);
  }
  {
    for (auto [relation, table, label] : {std::tuple{ap::Relation::isa, &B.compounds_isa, "is-a"}, std::tuple{ap::Relation::partof, &B.compounds_partof, "part-of"}}) {
      std::map<std::size_t, std::set<std::size_t>> given;
      for (const auto& [c, e] : *table) given[c].insert(e);
      std::size_t within = 0, equal = 0;
      std::string first;
      for (const auto& [c, elements] : given) {
        std::set<std::size_t> derived;
        for (std::size_t e = 0; e < B.elements.size(); e++)
          if (A.forest(relation).below(B.elements[e].concept_, c)) derived.insert(e);
        if (std::includes(elements.begin(), elements.end(), derived.begin(), derived.end())) within++;
        else if (first.empty()) first = A.concepts[c].name;
        equal += derived == elements;
      }
      s.require(std::string("the ") + label + " tree is a spanning tree of its relation: every compound's closure in the tree lies within its list",
                within == given.size(), first);
      std::printf("      %s: %zu of %zu compounds' lists are the tree's closure exactly\n", label, equal, given.size());
    }
  }
  female(s, dir);
  return s.done();
}
