// The parts of the body: concepts of the Foundational Model of Anatomy, the is-a and part-of relations
// between them, the elements that represent them, and the systems they belong to (ADR 0066;
// theory/anatomy.md §3).
//
// Concepts and relations. A concept is an FMA class (FMA7205, "left kidney"). The FMA relates classes by
// is-a (a left kidney is a kidney, is a parenchymatous organ, ...) and by part-of (the left kidney is part of
// the left upper urinary tract, of the urinary system, of the human body) (Rosse and Mejino 2003, J. Biomed.
// Inform. 36, 478, §3). Each relation, as BodyParts3D gives it, is a forest: every concept has at most one
// parent. A Forest holds it and answers the questions of a partial order: the ancestors of a concept (its
// principal up-set read as a chain), whether one concept lies below another (the reflexive transitive
// closure), their meet (the least common ancestor), and the down-set of a concept.
//
// Elements. An element is one mesh, the representation of one concept; a compound concept is represented by
// the elements below it. FMA's part-of is not a forest but a directed acyclic graph (a muscle of the body
// wall is part of the body wall and of its region), and BodyParts3D publishes a spanning tree of it, one
// parent each, beside lists of the elements composing each compound, which follow the whole graph. So an
// element lies below a concept when its own concept does in the tree, or when a compound it is listed under
// does (element_below); elements_of is that closure. The checks hold the tree to the lists: for every
// compound the tree's closure lies within its list (tests/check_body.cpp).
//
// Systems. The systems of the body cut across both relations: a vein is a vein by is-a wherever it lies,
// and a kidney is urinary by part-of. A system is assigned by an ordered table of rules, each "below the
// concept C in relation R"; the first rule that holds names the element's system. Tissue kinds come before
// organ systems, so that the vessels of the lung are vessels and not respiratory, and a few regional rules
// (the larynx's cartilages, the teeth) come before the kinds they would otherwise fall under. The table is
// the domain's; theory/anatomy.md §3.3 gives the reasons for its order. Checks: tests/check_parts.cpp.

#pragma once

#include <algorithm>
#include <array>
#include <cctype>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace anatomy::parts {

constexpr std::size_t none = static_cast<std::size_t>(-1);

enum class Relation { isa, partof };

// A relation as a forest on concepts 0..n−1.
struct Forest {
  std::vector<std::size_t> parent;                 // none at a root
  std::vector<std::vector<std::size_t>> children;  // in the order given

  explicit Forest(std::size_t n = 0) : parent(n, none), children(n) {}

  // Adds child → parent; a second parent, or an edge that closes a cycle, is refused.
  void relate(std::size_t child, std::size_t parent_) {
    if (child >= parent.size() || parent_ >= parent.size()) throw std::out_of_range("anatomy::parts::Forest::relate: a concept beyond the forest");
    if (parent[child] != none) {
      if (parent[child] == parent_) return;
      throw std::invalid_argument("anatomy::parts::Forest::relate: a concept with two parents");
    }
    for (std::size_t a = parent_; a != none; a = parent[a])
      if (a == child) throw std::invalid_argument("anatomy::parts::Forest::relate: a cycle");
    parent[child] = parent_;
    children[parent_].push_back(child);
  }
  // The chain from c up to its root, c first.
  std::vector<std::size_t> ancestors(std::size_t c) const {
    std::vector<std::size_t> out;
    for (std::size_t a = c; a != none; a = parent[a]) out.push_back(a);
    return out;
  }
  // c ≤ a: c lies in the down-set of a (reflexively).
  bool below(std::size_t c, std::size_t a) const {
    for (std::size_t x = c; x != none; x = parent[x])
      if (x == a) return true;
    return false;
  }
  std::size_t depth(std::size_t c) const { return ancestors(c).size() - 1; }
  // The meet of a and b: their least common ancestor, or none when they lie in different trees.
  std::size_t meet(std::size_t a, std::size_t b) const {
    const auto A = ancestors(a);
    for (std::size_t x = b; x != none; x = parent[x])
      if (std::find(A.begin(), A.end(), x) != A.end()) return x;
    return none;
  }
  // The down-set of c, c first, depth-first.
  std::vector<std::size_t> down_set(std::size_t c) const {
    std::vector<std::size_t> out, stack{c};
    while (!stack.empty()) {
      const std::size_t x = stack.back();
      stack.pop_back();
      out.push_back(x);
      for (auto it = children[x].rbegin(); it != children[x].rend(); ++it) stack.push_back(*it);
    }
    return out;
  }
  std::vector<std::size_t> roots() const {
    std::vector<std::size_t> out;
    for (std::size_t c = 0; c < parent.size(); c++)
      if (parent[c] == none && !children[c].empty()) out.push_back(c);
    return out;
  }
};

struct Concept {
  std::string id;    // FMA identifier, "FMA7205"
  std::string name;  // FMA's preferred name
  bool in_isa = false, in_partof = false;
};

// The concepts, the two forests, and the concept each element represents.
struct Atlas {
  std::vector<Concept> concepts;
  std::unordered_map<std::string, std::size_t> index;  // FMA id → concept
  Forest isa, partof;
  std::vector<std::size_t> element_concept;  // per element
  // The compounds each element is listed under in the source's tables, per relation (empty when none are given).
  std::vector<std::vector<std::size_t>> listed_isa, listed_partof;

  std::size_t add_concept(Concept c) {
    const auto it = index.find(c.id);
    if (it != index.end()) return it->second;
    const std::size_t k = concepts.size();
    index.emplace(c.id, k);
    concepts.push_back(std::move(c));
    isa.parent.push_back(none), isa.children.emplace_back();
    partof.parent.push_back(none), partof.children.emplace_back();
    return k;
  }
  std::size_t find(std::string_view id) const {
    const auto it = index.find(std::string(id));
    return it == index.end() ? none : it->second;
  }
  const Forest& forest(Relation r) const { return r == Relation::isa ? isa : partof; }
  Forest& forest(Relation r) { return r == Relation::isa ? isa : partof; }

  const std::vector<std::vector<std::size_t>>& listed(Relation r) const { return r == Relation::isa ? listed_isa : listed_partof; }

  // Element e lies below concept c in relation r: its own concept, or a compound it is listed under, lies below
  // c in the forest.
  bool element_below(std::size_t e, std::size_t c, Relation r) const {
    const auto& F = forest(r);
    if (F.below(element_concept.at(e), c)) return true;
    const auto& L = listed(r);
    if (e < L.size())
      for (std::size_t k : L[e])
        if (F.below(k, c)) return true;
    return false;
  }
  // The elements below concept c in relation r.
  std::vector<std::size_t> elements_of(std::size_t c, Relation r) const {
    std::vector<std::size_t> out;
    for (std::size_t e = 0; e < element_concept.size(); e++)
      if (element_below(e, c, r)) out.push_back(e);
    return out;
  }
};

// ---- systems ---------------------------------------------------------------------------------------------------

struct System {
  std::string_view key, name;
};

// The systems, in the order a reader meets them, from the surface inward.
inline constexpr std::array<System, 17> SYSTEMS{{
    {"skin", "Skin and hair"},
    {"skeleton", "Skeleton and joints"},
    {"muscles", "Muscles and tendons"},
    {"fascia", "Fascia"},
    {"heart", "Heart"},
    {"arteries", "Arteries"},
    {"veins", "Veins"},
    {"brain", "Brain and spinal cord"},
    {"nerves", "Nerves"},
    {"eye", "Eyes"},
    {"ear", "Ears"},
    {"respiratory", "Respiratory system"},
    {"digestive", "Digestive system"},
    {"urinary", "Urinary system"},
    {"genital", "Genital system"},
    {"endocrine", "Endocrine glands"},
    {"lymphoid", "Lymphoid organs"},
}};

struct SystemRule {
  std::string_view system;
  Relation relation;
  std::string_view concept_id;  // FMA id
  std::string_view concept_name;  // FMA's name, for the reader and the checks
};

// The ordered table (theory/anatomy.md §3.3): the first rule whose concept lies above an element's concept
// in its relation names the element's system.
inline constexpr SystemRule SYSTEM_RULES[] = {
    // Regional rules that come before the kinds they would fall under.
    {"digestive", Relation::isa, "FMA12516", "tooth"},
    {"respiratory", Relation::isa, "FMA55108", "laryngeal cartilage"},
    {"respiratory", Relation::isa, "FMA55250", "conus elasticus"},
    {"respiratory", Relation::isa, "FMA55244", "vocal ligament"},
    {"respiratory", Relation::isa, "FMA55233", "cricothyroid ligament"},
    {"respiratory", Relation::isa, "FMA55132", "thyrohyoid membrane"},
    {"respiratory", Relation::isa, "FMA55130", "epiglottis"},
    // Vessels, wherever they lie.
    {"arteries", Relation::isa, "FMA86187", "segment of arterial tree organ"},
    {"arteries", Relation::isa, "FMA63812", "set of arteries"},
    {"arteries", Relation::isa, "FMA3728", "arterial anastomosis"},
    {"arteries", Relation::isa, "FMA66332", "zone of artery"},
    {"veins", Relation::isa, "FMA86188", "segment of venous tree organ"},
    {"veins", Relation::isa, "FMA14284", "venous tree organ"},
    {"veins", Relation::isa, "FMA63814", "set of veins"},
    {"veins", Relation::isa, "FMA3729", "venous anastomosis"},
    // The nervous system.
    {"nerves", Relation::isa, "FMA11195", "segment of neural tree organ"},
    {"nerves", Relation::isa, "FMA61284", "neural tree organ"},
    {"nerves", Relation::isa, "FMA45638", "subdivision of nervous system"},
    {"brain", Relation::isa, "FMA256237", "segment of neuraxis"},
    {"brain", Relation::isa, "FMA83143", "cell part cluster of neuraxis"},
    {"brain", Relation::isa, "FMA83153", "organ component of neuraxis"},
    {"brain", Relation::isa, "FMA242770", "region of wall of ventricular system of neuraxis"},
    {"brain", Relation::isa, "FMA71235", "region of dura mater"},
    {"brain", Relation::isa, "FMA84881", "subdivision of subarachnoid space"},
    {"brain", Relation::isa, "FMA9338", "anatomical conduit space"},
    // Bone, cartilage and the ligaments of joints.
    {"skeleton", Relation::isa, "FMA5018", "bone organ"},
    {"skeleton", Relation::partof, "FMA7485", "sternum"},
    {"skeleton", Relation::isa, "FMA55107", "cartilage organ"},
    {"skeleton", Relation::isa, "FMA7591", "costal cartilage"},
    {"skeleton", Relation::isa, "FMA7538", "cartilage organ component"},
    {"skeleton", Relation::isa, "FMA21496", "ligament organ"},
    {"skeleton", Relation::isa, "FMA54839", "interosseous membrane"},
    // Muscle and tendon.
    {"muscles", Relation::isa, "FMA5022", "muscle organ"},
    {"muscles", Relation::isa, "FMA85453", "head of muscle organ"},
    {"muscles", Relation::isa, "FMA10474", "zone of muscle organ"},
    {"muscles", Relation::isa, "FMA32558", "musculature"},
    {"muscles", Relation::isa, "FMA9721", "tendon"},
    {"muscles", Relation::isa, "FMA32555", "ascending part of trapezius"},
    {"muscles", Relation::isa, "FMA32557", "descending part of trapezius"},
    {"muscles", Relation::isa, "FMA11336", "linea alba"},
    {"fascia", Relation::partof, "FMA79063", "deep fascial system"},
    {"fascia", Relation::isa, "FMA7646", "retinaculum"},
    // Organs, by their part-of systems and a few kinds the trees leave outside them.
    {"heart", Relation::isa, "FMA86212", "region of wall of heart"},
    {"heart", Relation::partof, "FMA7088", "heart"},
    {"heart", Relation::partof, "FMA7165", "right side of heart"},
    {"heart", Relation::partof, "FMA7166", "left side of heart"},
    {"eye", Relation::partof, "FMA54450", "left eye"},
    {"eye", Relation::partof, "FMA54449", "right eye"},
    {"eye", Relation::isa, "FMA225818", "region of layer of wall of eyeball"},
    {"eye", Relation::isa, "FMA58071", "chamber of eyeball"},
    {"eye", Relation::isa, "FMA49143", "check ligament of lateral rectus"},
    {"eye", Relation::isa, "FMA49146", "check ligament of medial rectus"},
    {"eye", Relation::isa, "FMA49071", "common tendinous ring"},
    {"ear", Relation::isa, "FMA52781", "external ear"},
    {"respiratory", Relation::isa, "FMA12224", "segment of tracheobronchial tree"},
    {"respiratory", Relation::partof, "FMA7158", "respiratory system"},
    {"digestive", Relation::isa, "FMA9597", "salivary gland"},
    {"digestive", Relation::isa, "FMA59816", "lip"},
    {"digestive", Relation::isa, "FMA55077", "pharyngeal raphe"},
    {"digestive", Relation::isa, "FMA55618", "pterygomandibular raphe"},
    {"digestive", Relation::isa, "FMA15737", "hepatovenous subsector"},
    {"digestive", Relation::isa, "FMA71856", "segment of biliary tree"},
    {"digestive", Relation::partof, "FMA7152", "alimentary system"},
    {"urinary", Relation::partof, "FMA7159", "urinary system"},
    {"genital", Relation::partof, "FMA7160", "genital system"},
    {"genital", Relation::isa, "FMA7210", "testis"},
    {"genital", Relation::isa, "FMA18255", "epididymis"},
    {"genital", Relation::isa, "FMA19234", "deferent duct"},
    {"genital", Relation::isa, "FMA19386", "seminal vesicle"},
    {"genital", Relation::isa, "FMA20194", "cavernous organ"},
    {"genital", Relation::isa, "FMA20221", "zone of cavernous organ"},
    {"lymphoid", Relation::isa, "FMA7196", "spleen"},
    {"lymphoid", Relation::isa, "FMA71193", "lobe of thymus"},
    {"lymphoid", Relation::partof, "FMA9607", "thymus"},
    {"endocrine", Relation::partof, "FMA9668", "endocrine system"},
    {"skin", Relation::partof, "FMA74657", "integument"},
    {"skin", Relation::isa, "FMA53667", "hair"},
    {"skin", Relation::isa, "FMA70752", "set of hairs"},
    // The rest of the nervous system by part-of, last.
    {"brain", Relation::partof, "FMA7157", "nervous system"},
};

inline std::size_t system_index(std::string_view key) {
  for (std::size_t k = 0; k < SYSTEMS.size(); k++)
    if (SYSTEMS[k].key == key) return k;
  return none;
}

// The system of a concept: the first rule that holds, or none.
inline std::size_t system_of_concept(const Atlas& A, std::size_t c) {
  for (const auto& r : SYSTEM_RULES) {
    const std::size_t k = A.find(r.concept_id);
    if (k != none && A.forest(r.relation).below(c, k)) return system_index(r.system);
  }
  return none;
}

// The Human Reference Atlas's bodies (the woman of its United, Female set; data/products/body_female.py) name their
// structures by UBERON and give no is-a tree: their concepts are the model's nodes, grouped by system layers
// (VH_F_circulatory_system, ...). Their systems are read from the layer an element lies in by part-of, and, where
// a layer holds more than one of the seventeen systems, from the element's name: in the circulatory layer veins
// and arteries by their names, the rest the heart; in the integumentary layer the breasts' parts (mammary gland,
// nipple, areola, lactiferous ducts and sinuses, the breast's suspensory ligaments), accessory organs of
// reproduction, with the genital system, the rest the skin. The ordered table: the first rule that holds.
struct LayerRule {
  std::string_view system;
  std::string_view layer;          // a concept (a node of the model) above the element by part-of
  std::string_view name_contains;  // empty: any name
};

inline constexpr LayerRule HRA_SYSTEM_RULES[] = {
    {"eye", "VH_F_eyes", ""},
    {"brain", "VH_F_nervous_system", ""},
    {"genital", "VH_F_integumentary_system", "mammary"},
    {"genital", "VH_F_integumentary_system", "nipple"},
    {"genital", "VH_F_integumentary_system", "areola"},
    {"genital", "VH_F_integumentary_system", "lactiferous"},
    {"genital", "VH_F_integumentary_system", "breast"},
    {"skin", "VH_F_integumentary_system", ""},
    {"muscles", "VH_F_muscular_system", ""},
    {"genital", "VH_F_reproductive_system", ""},
    {"digestive", "VH_F_digestive_system", ""},
    {"urinary", "VH_F_urinary_system", ""},
    {"veins", "VH_F_circulatory_system", "vein"},
    {"veins", "VH_F_circulatory_system", "vena"},
    {"veins", "VH_F_circulatory_system", "venous"},
    {"veins", "VH_F_circulatory_system", "sinus"},
    {"arteries", "VH_F_circulatory_system", "arter"},
    {"arteries", "VH_F_circulatory_system", "aorta"},
    {"arteries", "VH_F_circulatory_system", "trunk"},
    {"heart", "VH_F_circulatory_system", ""},
    {"respiratory", "VH_F_respiratory_system", ""},
    {"lymphoid", "VH_F_lymphatic_system", ""},
    {"skeleton", "VH_F_skeletal_system", ""},
};

// Whether an atlas is one of the Human Reference Atlas's (its root node VH_F, or VH_M for a man).
inline bool is_reference_atlas(const Atlas& A) { return A.find("VH_F") != none || A.find("VH_M") != none; }

inline std::size_t system_of_reference_element(const Atlas& A, std::size_t e) {
  const std::size_t c = A.element_concept[e];
  const std::string& name = A.concepts[c].name;
  for (const auto& r : HRA_SYSTEM_RULES) {
    const std::size_t k = A.find(r.layer);
    if (k == none || !A.forest(Relation::partof).below(c, k)) continue;
    if (!r.name_contains.empty() && name.find(r.name_contains) == std::string::npos) continue;
    return system_index(r.system);
  }
  return none;
}

// The system of an element: the first rule under which it lies, by its concept or its listed compounds.
inline std::size_t system_of_element(const Atlas& A, std::size_t e) {
  if (is_reference_atlas(A)) return system_of_reference_element(A, e);
  for (const auto& r : SYSTEM_RULES) {
    const std::size_t k = A.find(r.concept_id);
    if (k != none && A.element_below(e, k, r.relation)) return system_index(r.system);
  }
  return none;
}

// The rules of the table whose concepts the atlas lacks (a table written against another release).
inline std::vector<std::string_view> missing_rules(const Atlas& A) {
  std::vector<std::string_view> out;
  for (const auto& r : SYSTEM_RULES)
    if (A.find(r.concept_id) == none) out.push_back(r.concept_id);
  return out;
}

// ---- sides -----------------------------------------------------------------------------------------------------

enum class Side { none, left, right };

// The side a name gives: FMA writes it as a word, "left kidney", "flexor retinaculum of right wrist".
inline Side side_by_name(std::string_view name) {
  std::string s(name);
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
  auto word = [&](std::string_view w) {
    for (std::size_t at = s.find(w); at != std::string::npos; at = s.find(w, at + 1)) {
      const bool start = at == 0 || !std::isalpha(static_cast<unsigned char>(s[at - 1]));
      const bool end = at + w.size() == s.size() || !std::isalpha(static_cast<unsigned char>(s[at + w.size()]));
      if (start && end) return true;
    }
    return false;
  };
  const bool l = word("left"), r = word("right");
  return l && !r ? Side::left : r && !l ? Side::right : Side::none;
}

}  // namespace anatomy::parts
