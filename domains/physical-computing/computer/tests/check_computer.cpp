// Checks of the formalism's contract and accounts (square.hpp, cost.hpp, materials.hpp; ADR 0076).
//
// - The square: a substrate that doubles by a physical rescaling (a lever) computes x ↦ 2x on K within its
//   tolerance, and the verdict names the worst case; a readout off by a bias breaks it; a non-finite reading is
//   never within tolerance; an empty K proves nothing.
// - Cost: fixed costs spread over a lifetime of results, the per-result terms added; a missing stage is
//   reported, not taken as zero; the declared share; the Landauer floor of a byte at 300 K, 8 k_B T ln 2;
//   an ADC of 8 effective bits at 5 fJ per step costs 1.28 pJ a sample.
// - Materials: silicon metal, copper and the platinum group metals are critical, water and silicone are not;
//   the list has 34 entries; ambient conditions.

#include <physical_computing/computer/cost.hpp>
#include <physical_computing/computer/materials.hpp>
#include <physical_computing/computer/square.hpp>
#include <physical_computing/tasks/task.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <limits>
#include <vector>

using namespace physical_computing;

int main() {
  calx::verification::Suite s;
  {
    const std::vector<double> K{0.1, 0.5, 1.0, 2.5};
    auto lever = [](double x) { return 2 * x * (1 + 1e-6 * x); };  // a slightly nonlinear lever
    const auto v = computer::commutes(K, [](double m) { return m; }, lever, [](double y) { return y; }, [](double m) { return 2 * m; },
                                      [](double a, double b) { return std::fabs(a - b); }, 1e-4, computer::Map::rescaling, computer::Map::identity);
    s.require("square: a lever computes x ↦ 2x within 1e−4 on K", v.holds && v.cases == 4);
    s.require("square: the worst case is the largest x", v.worst == 3);
    s.check("square: its error 2 × 1e−6 × 2.5²", v.error - 2e-6 * 6.25, 1e-15);
    const auto biased = computer::commutes(K, [](double m) { return m; }, lever, [](double y) { return y + 1e-3; }, [](double m) { return 2 * m; },
                                           [](double a, double b) { return std::fabs(a - b); }, 1e-4);
    s.require("square: a biased readout breaks it", !biased.holds);
    const auto nan = computer::commutes(K, [](double m) { return m; }, [](double) { return std::numeric_limits<double>::quiet_NaN(); },
                                        [](double y) { return y; }, [](double m) { return m; }, [](double a, double b) { return std::fabs(a - b); }, 1.0);
    s.require("square: a non-finite reading is never within tolerance", !nan.holds);
    const std::vector<double> none;
    const auto empty = computer::commutes(none, [](double m) { return m; }, [](double m) { return m; }, [](double m) { return m; },
                                          [](double m) { return m; }, [](double a, double b) { return std::fabs(a - b); }, 1.0);
    s.require("square: an empty set proves nothing", !empty.holds);
  }
  {
    computer::CostAccount a;
    a.add({"fabrication of the tank", computer::Stage::fabrication, 5.0e6, computer::Basis::declared, "an estimate"});
    a.add({"pumps and heaters", computer::Stage::operation, 120.0, computer::Basis::derived, "power × time"});
    a.add({"camera's conversions", computer::Stage::conversion, 3.0, computer::Basis::declared, "sensor's datasheet"});
    s.check("cost: fixed costs spread over 10⁶ results plus the per-result terms (J)", a.per_verified_result(1e6) - (5.0 + 123.0), 1e-12, "J");
    const auto missing = a.missing_stages();
    s.require("cost: the stages without a term are reported (materials, verification)", missing.size() == 2 && missing[0] == computer::Stage::materials && missing[1] == computer::Stage::verification);
    s.check("cost: the declared share (8/128)", a.declared_share(1e6) - 8.0 / 128, 1e-15);
    s.check("cost: the Landauer floor of a byte at 300 K, 8 k_B T ln 2 (J)", computer::landauer_floor(300, 8) - 8 * 1.380649e-23 * 300 * std::log(2.0), 1e-36, "J");
    s.check("cost: an ADC of 8 bits at 5 fJ per step, 1.28 pJ a sample (J)", computer::conversion_energy(8, 5e-15) - 1.28e-12, 1e-27, "J");
  }
  {
    s.require("materials: silicon metal, copper and the platinum group metals are critical",
              computer::is_critical("silicon metal") && computer::is_critical("copper") && computer::is_critical("platinum group metals"));
    s.require("materials: water and silicone (PDMS) are not", !computer::is_critical("water") && !computer::is_critical("silicone (PDMS)"));
    s.require("materials: the list has 34 entries", computer::critical_raw_materials.size() == 34);
    computer::Marks m;
    m.materials = {{"water", "the medium"}, {"copper", "the motors' windings"}, {"acrylic (PMMA)", "the tank"}};
    s.require("materials: a tank's critical materials are its copper alone", m.critical().size() == 1 && m.critical()[0] == "copper");
    s.require("materials: at room temperature without vacuum the conditions are ambient", m.ambient());
    m.operating_temperature = 4.2;
    s.require("materials: at 4.2 K they are not", !m.ambient());
  }
  s.require("tasks: six kinds, each named", tasks::all_kinds.size() == 6 && tasks::name(tasks::Kind::phenomenon) == "phenomenon");
  return s.done();
}
