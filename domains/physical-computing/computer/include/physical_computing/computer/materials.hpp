// The marks of a substrate's materials (theory note §1.5; ADR 0076 §2.3).
//
// Five marks: supply (the critical raw materials it needs), toxicity and persistence, reuse, the conditions of
// operation, and the end of life. Materials are named as the raw materials the supply mark is judged on, so
// that "silicon metal" (critical) is not confused with silicone (a polymer that needs no silicon metal in the
// substrate itself).
//
// The list of critical raw materials is that of Regulation (EU) 2024/1252 (the Critical Raw Materials Act),
// Annex II, Section 1 (Section 2 is the method of its calculation): 34 entries, which include every strategic raw
// material of Annex I, copper and battery-grade nickel among them. Its count and first entries agree with EIT
// RawMaterials' copy of the Annex ("bauxite/alumina/aluminium" is entered here as "aluminium"); the whole list is
// to be checked against the Official Journal before it is relied on. *Caveat.* The list measures economic
// importance and the risk to supply, not scarcity or harm; it is a mark among five, not a verdict.
#pragma once

#include <algorithm>
#include <array>
#include <string>
#include <string_view>
#include <vector>

namespace physical_computing::computer {

inline constexpr std::array<std::string_view, 34> critical_raw_materials{
    "aluminium", "antimony", "arsenic", "baryte", "beryllium", "bismuth", "boron", "cobalt", "coking coal", "copper",
    "feldspar", "fluorspar", "gallium", "germanium", "hafnium", "helium", "heavy rare earth elements",
    "light rare earth elements", "lithium", "magnesium", "manganese", "natural graphite", "nickel (battery grade)",
    "niobium", "phosphate rock", "phosphorus", "platinum group metals", "scandium", "silicon metal", "strontium",
    "tantalum", "titanium metal", "tungsten", "vanadium"};

inline bool is_critical(std::string_view raw_material) {
  return std::find(critical_raw_materials.begin(), critical_raw_materials.end(), raw_material) != critical_raw_materials.end();
}

enum class Reuse { reusable, consumed_per_run, degrades_with_use };
enum class EndOfLife { biodegradable, recyclable, neither, unknown };

struct Material {
  std::string raw_material;  // as the supply mark names it: "copper", "water", "silicone (PDMS)"
  std::string role;          // what it does in the substrate
};

struct Marks {
  std::vector<Material> materials;
  std::vector<std::string> hazards;  // toxicity and persistence, each with its source
  Reuse reuse = Reuse::reusable;
  double operating_temperature = 293.15;  // K
  bool vacuum = false;
  EndOfLife end_of_life = EndOfLife::unknown;

  // The critical raw materials among the materials.
  std::vector<std::string> critical() const {
    std::vector<std::string> out;
    for (const auto& m : materials)
      if (is_critical(m.raw_material)) out.push_back(m.raw_material);
    return out;
  }
  // Room temperature and pressure: no cryogenics, no furnace, no vacuum.
  bool ambient() const { return !vacuum && operating_temperature > 250 && operating_temperature < 350; }
};

}  // namespace physical_computing::computer
