// Checks of the parts of the body (parts.hpp), on a small atlas built here.
//
// - A forest: ancestors as a chain, the order (below) reflexive and transitive, the meet of two concepts,
//   the down-set; a second parent and a cycle refused.
// - Elements of a compound are the closure: the elements whose concepts lie below it.
// - Systems by the ordered table: a femur is skeleton by is-a; a kidney urinary by part-of; a vein of the
//   lung is a vein, not respiratory, because tissue kinds come before organ systems; a cartilage of the
//   larynx is respiratory, because that regional rule comes before cartilage; a concept below no rule has no
//   system.
// - Sides read from names as words: "left kidney", "flexor retinaculum of right wrist"; neither in
//   "leftover" nor where both are named.

#include <string>
#include <vector>

#include <anatomy/parts/parts.hpp>
#include <calx/verification/check.hpp>

using calx::verification::Suite;
namespace ap = anatomy::parts;

int main() {
  Suite s;
  ap::Atlas A;
  auto add = [&](const char* id, const char* name) { return A.add_concept({id, name, true, true}); };
  // is-a
  const auto entity = add("FMA62955", "anatomical entity");
  const auto bone = add("FMA5018", "bone organ");
  const auto femur = add("FMA9611", "femur");
  const auto vein_kind = add("FMA86188", "segment of venous tree organ");
  const auto pulmonary_vein = add("FMA49913", "left inferior pulmonary vein");
  const auto cartilage = add("FMA55107", "cartilage organ");
  const auto laryngeal = add("FMA55108", "laryngeal cartilage");
  const auto thyroid_cartilage = add("FMA9615", "thyroid cartilage");
  // part-of
  const auto body = add("FMA20394", "human body");
  const auto urinary = add("FMA7159", "urinary system");
  const auto kidney = add("FMA7205", "left kidney");
  const auto respiratory = add("FMA7158", "respiratory system");
  const auto lung = add("FMA7310", "left lung");
  const auto orphan = add("FMA99999", "a concept below no rule");
  A.isa.relate(bone, entity);
  A.isa.relate(femur, bone);
  A.isa.relate(vein_kind, entity);
  A.isa.relate(pulmonary_vein, vein_kind);
  A.isa.relate(cartilage, entity);
  A.isa.relate(laryngeal, cartilage);
  A.isa.relate(thyroid_cartilage, laryngeal);
  A.partof.relate(urinary, body);
  A.partof.relate(kidney, urinary);
  A.partof.relate(respiratory, body);
  A.partof.relate(lung, respiratory);
  A.partof.relate(pulmonary_vein, lung);
  A.partof.relate(thyroid_cartilage, respiratory);
  {
    s.require("ancestors as a chain, the concept first", A.isa.ancestors(femur) == std::vector<std::size_t>{femur, bone, entity});
    s.require("the order is reflexive and transitive", A.isa.below(femur, femur) && A.isa.below(femur, entity) && !A.isa.below(entity, femur));
    s.require("the meet of a femur and a pulmonary vein is the anatomical entity", A.isa.meet(femur, pulmonary_vein) == entity);
    s.require("the meet of a kidney and a lung is the human body; of concepts in different trees, none",
              A.partof.meet(kidney, lung) == body && A.partof.meet(kidney, femur) == ap::none);
    const auto d = A.partof.down_set(respiratory);
    s.require("the down-set of the respiratory system: itself, the lung, its vein, the thyroid cartilage", d.size() == 4 && d[0] == respiratory);
    bool refused_parent = false, refused_cycle = false;
    try { A.isa.relate(femur, cartilage); } catch (const std::invalid_argument&) { refused_parent = true; }
    try { A.isa.relate(entity, femur); } catch (const std::invalid_argument&) { refused_cycle = true; }
    s.require("a second parent is refused, and so is a cycle", refused_parent && refused_cycle);
  }
  {
    A.element_concept = {femur, kidney, pulmonary_vein, thyroid_cartilage, orphan};
    s.require("the elements of the respiratory system by part-of: the vein and the cartilage",
              A.elements_of(respiratory, ap::Relation::partof) == std::vector<std::size_t>{2, 3});
    s.require("the elements of the anatomical entity by is-a: all four below it", A.elements_of(entity, ap::Relation::isa) == std::vector<std::size_t>{0, 2, 3});
    // FMA's part-of is a graph: the source lists the femur's element under the urinary system too (a second
    // parent the spanning tree leaves out); listed membership counts beside the tree.
    A.listed_partof.assign(A.element_concept.size(), {});
    A.listed_partof[0].push_back(urinary);
    s.require("an element listed under a compound lies below it and its ancestors, beside its place in the tree",
              A.element_below(0, urinary, ap::Relation::partof) && A.element_below(0, body, ap::Relation::partof) &&
                  !A.element_below(0, respiratory, ap::Relation::partof) && A.elements_of(urinary, ap::Relation::partof) == std::vector<std::size_t>{0, 1});
  }
  {
    auto sys = [&](std::size_t c) { const auto k = ap::system_of_concept(A, c); return k == ap::none ? std::string("none") : std::string(ap::SYSTEMS[k].key); };
    s.require("a femur is skeleton, by is-a", sys(femur) == "skeleton");
    s.require("a kidney is urinary, by part-of", sys(kidney) == "urinary");
    s.require("a vein of the lung is a vein: tissue kinds come before organ systems", sys(pulmonary_vein) == "veins");
    s.require("a cartilage of the larynx is respiratory: that regional rule comes before cartilage", sys(thyroid_cartilage) == "respiratory");
    s.require("a concept below no rule has no system", sys(orphan) == "none");
    bool known = true;
    for (const auto& r : ap::SYSTEM_RULES) known = known && ap::system_index(r.system) != ap::none;
    s.require("every rule names one of the seventeen systems", known && ap::SYSTEMS.size() == 17);
  }
  {
    using ap::Side;
    s.require("sides read from names as words",
              ap::side_by_name("Left kidney") == Side::left && ap::side_by_name("Flexor retinaculum of right wrist") == Side::right &&
                  ap::side_by_name("leftover tissue") == Side::none && ap::side_by_name("Left and right x") == Side::none && ap::side_by_name("Aorta") == Side::none);
  }
  return s.done();
}
