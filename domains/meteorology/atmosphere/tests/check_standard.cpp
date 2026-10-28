// Checks of the U.S. Standard Atmosphere 1976 (meteorology/atmosphere/standard.hpp).
//
// References (tests/fixtures/atmosphere-reference.txt, made by tests/make_atmosphere_fixtures.py):
// the document's printed Tables I-III at thirty heights from -1 km to 86 km, each held to half a unit
// in its last printed digit (at 85 km corrected by Table 8 as the document's Sec. 1.2.4 directs);
// fluids' ATMOSPHERE_1976, which shares the document's constants; ambiance's ICAO 1993 atmosphere to
// 80 km, which shares its temperature profile but not all its constants. And by construction: the
// continuity of temperature and pressure at the bases of the layers, the hydrostatic equation
// dP/dZ = -rho g by finite differences, the inverse of the pressure, Table 8 and T7 (eq. 25), the
// bounds of the model, and the columns of levels.

#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/verification/check.hpp>
#include <meteorology/atmosphere/standard.hpp>

using namespace meteorology::atmosphere;
using calx::verification::relative;
using calx::verification::Suite;

namespace {

struct TableEntry {
  std::string quantity;
  double z, value, half_unit;
};

// The model's value of a quantity named as the fixture names it.
double quantity(const State& s, const std::string& q) {
  if (q == "T") return s.temperature_K;
  if (q == "TM") return s.molecular_temperature_K;
  if (q == "P") return s.pressure_Pa;
  if (q == "rho") return s.density_kg_m3;
  if (q == "g") return s.gravity_m_s2;
  if (q == "Hp") return s.pressure_scale_height_m;
  if (q == "n") return s.number_density_per_m3;
  if (q == "V") return s.mean_particle_speed_m_s;
  if (q == "nu") return s.collision_frequency_per_s;
  if (q == "L") return s.mean_free_path_m;
  if (q == "Cs") return s.speed_of_sound_m_s;
  if (q == "mu") return s.dynamic_viscosity_Pa_s;
  if (q == "eta") return s.kinematic_viscosity_m2_s;
  throw std::invalid_argument("unknown quantity " + q);
}

}  // namespace

int main() {
  Suite s;
  using namespace us1976;

  // ---- read the references ------------------------------------------------------------------------------
  std::vector<TableEntry> table, table_geopotential;
  std::vector<std::vector<double>> fluids, ambiance;
  {
    std::ifstream in(METEOROLOGY_FIXTURES "/atmosphere-reference.txt");
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream f(line);
      std::string kind;
      f >> kind;
      if (kind == "table" || kind == "tableH") {
        TableEntry e;
        int page;
        f >> e.quantity >> e.z >> e.value >> e.half_unit >> page;
        (kind == "table" ? table : table_geopotential).push_back(e);
      } else {
        std::vector<double> row;
        double v;
        while (f >> v) row.push_back(v);
        (kind == "fluids" ? fluids : ambiance).push_back(row);
      }
    }
  }
  s.require("the references were read (272 + 24 table entries, fluids and ambiance rows)",
            table.size() == 272 && table_geopotential.size() == 24 && fluids.size() > 300 && ambiance.size() > 250);

  // ---- heights and gravity (eqs. 17-19) ----------------------------------------------------------------------
  {
    double round_trip = 0;
    for (double z = -5000; z <= 86000; z += 37.5) round_trip = std::max(round_trip, std::fabs(geometric_height(geopotential_height(z)) - z));
    s.check("geometric to geopotential height and back (eqs. 18, 19)", round_trip, 1e-8, "m");
    // Table 8 prints its pairs to 0.1 m but does not round them consistently (Z = 81 km: 79 980.8 m' for
    // 79 980.86; H = 79 500 m': 80 506.9 m for 80 506.85), so they are held to one unit of 0.1 m.
    double pairs = 0;
    for (const auto& [z, h] : std::vector<std::pair<double, double>>{{80500, 79493.3}, {81000, 79980.8}, {82000, 80955.7}, {84000, 82904.4}, {86000, 84852.0}})
      pairs = std::max(pairs, std::fabs(geopotential_height(z) - h));
    for (const auto& [h, z] : std::vector<std::pair<double, double>>{{79500, 80506.9}, {80000, 81019.6}, {83000, 84098.0}, {84000, 85124.8}})
      pairs = std::max(pairs, std::fabs(geometric_height(h) - z));
    s.check("Table 8's nine pairs of geometric and geopotential heights (eqs. 18, 19)", pairs, 0.1, "m");
    s.check("gravity at sea level is g0 (eq. 17)", gravity(0) - g0_m_s2, 1e-15, "m/s2");
  }

  // ---- the layers (Table 4): continuity at their bases --------------------------------------------------------
  {
    const auto& t = layers();
    double dT = 0, dP = 0, jumpT = 0, jumpP = 0;
    for (int b = 1; b < 7; b++) {
      const double H = t[b].base_geopotential_m;
      dT = std::max(dT, std::fabs(molecular_temperature_in_layer(b - 1, H) - t[b].base_temperature_K));
      dP = std::max(dP, relative(pressure_in_layer(b - 1, H), t[b].base_pressure_Pa));
      const auto below = at_geopotential_height(H - 1e-6), above = at_geopotential_height(H + 1e-6);
      jumpT = std::max(jumpT, std::fabs(above->molecular_temperature_K - below->molecular_temperature_K));
      jumpP = std::max(jumpP, relative(above->pressure_Pa, below->pressure_Pa));
    }
    s.check("T_M of each layer meets the next layer's base (eq. 23, six bases)", dT, 1e-11, "K");
    s.check("P of each layer meets the next layer's base (eqs. 33a, 33b)", dP, 1e-14);
    s.check("no jump of T_M across a base (at H_b -/+ 1 um)", jumpT, 1e-8, "K");
    s.check("no jump of P across a base (relative)", jumpP, 1e-9);
    s.check("T_M7 = 186.946 K at H7 = 84.852 km' (p. 10: T7 + 0.0787 K)", molecular_temperature_in_layer(6, H7_geopotential_m) - 186.946, 1e-9, "K");
  }

  // ---- sea level ----------------------------------------------------------------------------------------------
  {
    const State z = *at_geometric_height(0);
    s.check("sea level: T0 = 288.15 K and P0 = 101 325 Pa", std::fabs(z.temperature_K - T0_K) + std::fabs(z.pressure_Pa - P0_Pa), 1e-12);
    s.require("sea level: layer 0, the troposphere", z.layer == 0 && z.layer_name == "troposphere");
    // The 1976's R* is "not exactly consistent with the cited values of k and N_A" (p. 3): P = N k T holds to 2.3e-6.
    s.check("the ideal-gas law P = N k T with Table 2's k, to the inconsistency p. 3 admits", z.number_density_per_m3 * k_J_K * z.temperature_K / z.pressure_Pa - 1, 3e-6);
    s.check("the mole volume times the density is the molecular weight (eqs. 42, 43)", relative(z.mole_volume_m3_mol * z.density_kg_m3, z.molecular_weight_kg_mol), 1e-14);
  }

  // ---- the document's Tables I-III ------------------------------------------------------------------------------
  {
    // Table I's pressure column (mb) is truncated, not rounded: the equations' value less the printed one
    // lies between 0 and one unit of the last digit at all 36 heights (226.9996 mb printed 226.99 at
    // Z = 11 km; 845.560 printed 845.55 at 1500 m', beside a P/P0 of 0.834 50 that is 845.560 mb).
    double low = 0, high = 0;
    int pressures = 0;
    for (const auto* entries : {&table, &table_geopotential})
      for (const auto& e : *entries) {
        if (e.quantity != "P") continue;
        const State st = entries == &table ? *at_geometric_height(e.z) : *at_geopotential_height(e.z);
        const double units = (st.pressure_Pa - e.value) / (2 * e.half_unit);
        low = std::min(low, units);
        high = std::max(high, units);
        pressures++;
      }
    s.check("Table I: the printed pressure is the equations' truncated, at " + std::to_string(pressures) + " heights (excess beyond [0, 1] unit of the last digit)",
            std::max({0.0, -low, high - 1.0}), 0.05, "unit");
    // The other columns are rounded. On the geopotential pages they agree with the equations to half a unit
    // of the last printed digit; ties (9.806 65 printed 9.8066) reach it exactly.
    std::map<std::string, double> worstH;
    for (const auto& e : table_geopotential) {
      const State st = *at_geopotential_height(e.z);
      worstH[e.quantity] = std::max(worstH[e.quantity], std::fabs(quantity(st, e.quantity) - e.value) / e.half_unit);
    }
    s.check("Table I, geopotential pages: temperature at 8 heights, in half-units of the last printed digit", worstH["T"], 1.0 + 1e-6);
    s.check("Table I, geopotential pages: density at 8 heights, in half-units", worstH["rho"], 1.0 + 1e-6);
    // On the geometric pages some derived columns depart by up to one unit (two half-units): the density at
    // 2 and 85 km, the number density and collision frequency at 50-80 km. The document does not say why.
    std::map<std::string, double> worst;
    std::map<std::string, int> count;
    for (const auto& e : table) {
      if (e.quantity == "P") continue;
      const State st = *at_geometric_height(e.z);
      worst[e.quantity] = std::max(worst[e.quantity], std::fabs(quantity(st, e.quantity) - e.value) / e.half_unit);
      count[e.quantity]++;
    }
    const std::map<std::string, std::string> names{
        {"T", "Table I: kinetic temperature T"}, {"TM", "Table I: molecular-scale temperature T_M at 86 km"}, {"rho", "Table I: density"}, {"g", "Table II: gravity"}, {"Hp", "Table II: pressure scale height"}, {"n", "Table II: number density"},
        {"V", "Table II: mean particle speed"}, {"nu", "Table II: collision frequency"}, {"L", "Table II: mean free path"},
        {"Cs", "Table III: speed of sound"}, {"mu", "Table III: dynamic viscosity"}, {"eta", "Table III: kinematic viscosity"}};
    for (const auto& [q, name] : names)
      s.check(name + " at " + std::to_string(count[q]) + " geometric heights, in half-units of the last printed digit", worst[q], 2.0 + 1e-6);
  }

  // ---- fluids' ATMOSPHERE_1976 ------------------------------------------------------------------------------------
  {
    // fluids: Z T_M P rho C_s mu k g. Its T is T_M (it has no Table 8); its viscosity and conductivity are
    // compared below 80 km, where T = T_M. At Z = 86 km exactly (5 cm above H7) it switches to its
    // isothermal layer 7, which is not this model; that one row is left out.
    double TM = 0, P = 0, rho = 0, cs = 0, mu = 0, k = 0, g = 0;
    int n = 0;
    for (const auto& r : fluids) {
      if (geopotential_height(r[0]) > H7_geopotential_m) continue;
      const State st = *at_geometric_height(r[0]);
      TM = std::max(TM, relative(st.molecular_temperature_K, r[1]));
      P = std::max(P, relative(st.pressure_Pa, r[2]));
      rho = std::max(rho, relative(st.density_kg_m3, r[3]));
      cs = std::max(cs, relative(st.speed_of_sound_m_s, r[4]));
      if (r[0] <= 80000) {
        mu = std::max(mu, relative(st.dynamic_viscosity_Pa_s, r[5]));
        k = std::max(k, relative(st.thermal_conductivity_W_m_K, r[6]));
      }
      g = std::max(g, relative(st.gravity_m_s2, r[7]));
      n++;
    }
    const std::string at = " (fluids, " + std::to_string(n) + " heights to 86 km, relative)";
    s.check("molecular-scale temperature T_M" + at, TM, 1e-13);
    s.check("pressure" + at, P, 1e-12);
    s.check("density" + at, rho, 1e-12);
    s.check("speed of sound" + at, cs, 1e-13);
    s.check("gravity" + at, g, 1e-13);
    s.check("dynamic viscosity (fluids, to 80 km, relative)", mu, 1e-13);
    s.check("thermal conductivity, eq. 53 (fluids, to 80 km, relative)", k, 1e-13);
  }

  // ---- ambiance's ICAO 1993 atmosphere, to 80 km ---------------------------------------------------------------
  {
    // ambiance: Z H T P rho C_s mu eta n L V nu g H_P. ICAO 1993 rounds its base pressures to six figures
    // (Table D: 2.26320e4 Pa at 11 km', against the 1976's 22 632.06), so pressure, density and the
    // kinematic viscosity differ by up to 1e-5; its R = 287.052 87 J kg-1 K-1 against R*/M0 = 287.053 07
    // moves the speeds and the scale height by 4e-7; its N_A = 6.022 57e23 mol-1 against the 1976's
    // 6.022 169e23 moves the number density, the mean free path and the collision frequency by 6.7e-5.
    double H = 0, T = 0, P = 0, rho = 0, cs = 0, mu = 0, eta = 0, num = 0, L = 0, V = 0, nu = 0, g = 0, Hp = 0;
    for (const auto& r : ambiance) {
      const State st = *at_geometric_height(r[0]);
      H = std::max(H, std::fabs(st.geopotential_height_m - r[1]));
      T = std::max(T, std::fabs(st.temperature_K - r[2]));
      P = std::max(P, relative(st.pressure_Pa, r[3]));
      rho = std::max(rho, relative(st.density_kg_m3, r[4]));
      cs = std::max(cs, relative(st.speed_of_sound_m_s, r[5]));
      mu = std::max(mu, relative(st.dynamic_viscosity_Pa_s, r[6]));
      eta = std::max(eta, relative(st.kinematic_viscosity_m2_s, r[7]));
      num = std::max(num, relative(st.number_density_per_m3, r[8]));
      L = std::max(L, relative(st.mean_free_path_m, r[9]));
      V = std::max(V, relative(st.mean_particle_speed_m_s, r[10]));
      nu = std::max(nu, relative(st.collision_frequency_per_s, r[11]));
      g = std::max(g, relative(st.gravity_m_s2, r[12]));
      Hp = std::max(Hp, relative(st.pressure_scale_height_m, r[13]));
    }
    const std::string at = " (ambiance, ICAO 1993, " + std::to_string(ambiance.size()) + " heights to 80 km)";
    s.check("geopotential height" + at, H, 1e-8, "m");
    s.check("temperature" + at, T, 1e-10, "K");
    s.check("gravity, relative" + at, g, 1e-13);
    s.check("dynamic viscosity, relative" + at, mu, 1e-13);
    s.check("pressure, relative: ICAO's six-figure base pressures" + at, P, 1.5e-5);
    s.check("density, relative: the same" + at, rho, 1.5e-5);
    s.check("kinematic viscosity, relative: the same" + at, eta, 1.5e-5);
    s.check("speed of sound, relative: ICAO's R" + at, cs, 1e-6);
    s.check("mean particle speed, relative: ICAO's R" + at, V, 1e-6);
    s.check("pressure scale height, relative: ICAO's R" + at, Hp, 1e-6);
    s.check("number density, relative: ICAO's N_A" + at, num, 8e-5);
    s.check("mean free path, relative: ICAO's N_A" + at, L, 8e-5);
    s.check("collision frequency, relative: ICAO's N_A" + at, nu, 8e-5);
  }

  // ---- the hydrostatic equation, dP/dZ = -rho g ---------------------------------------------------------------
  {
    std::mt19937_64 rng(20261028);
    std::uniform_real_distribution<double> height(-4990.0, 85990.0);
    double worst = 0;
    const double d = 0.25;
    for (int i = 0; i < 2000; i++) {
      const double z = height(rng);
      const double dPdZ = (at_geometric_height(z + d)->pressure_Pa - at_geometric_height(z - d)->pressure_Pa) / (2 * d);
      const State st = *at_geometric_height(z);
      worst = std::max(worst, std::fabs(dPdZ + st.density_kg_m3 * st.gravity_m_s2) / (st.density_kg_m3 * st.gravity_m_s2));
    }
    s.check("hydrostatic equilibrium dP/dZ = -rho g, by central differences at 2000 heights (relative)", worst, 1e-7);
  }

  // ---- between 80 and 86 km: Table 8 and T7 ------------------------------------------------------------------
  {
    double nodes = 0;
    bool falling = true;
    for (int i = 0; i <= 12; i++) {
      nodes = std::max(nodes, std::fabs(molecular_weight_ratio(80000 + 500.0 * i) - molecular_weight_ratios[i]));
      if (i > 0 && !(molecular_weight_ratios[i] < molecular_weight_ratios[i - 1])) falling = false;
    }
    s.check("M/M0 at Table 8's nodes, 80.0 to 86.0 km", nodes, 1e-15);
    s.require("M/M0 falls from 80 to 86 km, and is 1 below 80 km", falling && molecular_weight_ratio(79999.0) == 1.0 && molecular_weight_ratio(-3000) == 1.0);
    const State top = *at_geometric_height(86000);
    s.check("T7 = 186.8673 K at 86 km (eq. 25)", top.temperature_K - 186.8673, 5e-4, "K");
    s.check("and 0.0787 K below T_M there (p. 10)", top.molecular_temperature_K - top.temperature_K - 0.0787, 5e-4, "K");
    s.check("M = 28.9522 kg/kmol at 86 km (p. 9)", top.molecular_weight_kg_mol * 1e3 - 28.9522, 5e-5, "kg/kmol");
    s.check("T is continuous at 80 km, where the correction begins", at_geometric_height(80000.001)->temperature_K - at_geometric_height(79999.999)->temperature_K, 1e-5, "K");
  }

  // ---- the bounds of the model and the names of its layers ------------------------------------------------------
  {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    s.require("no value above 86 km (86 000.001 m, 90 km, 500 km): the upper model is not implemented",
              !at_geometric_height(86000.001) && !at_geometric_height(90000) && !at_geometric_height(5e5) && !at_geopotential_height(84853));
    s.require("no value below -5 km' or for NaN", !at_geopotential_height(-5000.001) && !at_geometric_height(-5005) && !at_geometric_height(nan) && !layer_index(nan));
    s.require("a value at 86 km and at -5 km'", at_geometric_height(86000).has_value() && at_geopotential_height(-5000).has_value() && at_geopotential_height(top_geopotential_m()).has_value());
    bool names = true;
    for (const auto& [z, name] : std::vector<std::pair<double, std::string>>{{-1000, "troposphere"}, {5000, "troposphere"}, {15000, "tropopause"}, {30000, "stratosphere"},
                                                                            {40000, "stratosphere"}, {49000, "stratopause"}, {60000, "mesosphere"}, {80000, "mesosphere"}})
      names = names && at_geometric_height(z)->layer_name == name;
    s.require("the layers' names: troposphere, tropopause, stratosphere, stratopause, mesosphere", names);
    s.require("the base of a layer belongs to it (11 km' is layer 1, 47 km' layer 4)", *layer_index(11000) == 1 && *layer_index(10999.999) == 0 && *layer_index(47000) == 4);
  }

  // ---- the pressure altitude ---------------------------------------------------------------------------------------
  {
    double worst = 0;
    for (double H = -5000; H <= top_geopotential_m(); H += 113.7) worst = std::max(worst, std::fabs(*geopotential_height_at_pressure(at_geopotential_height(H)->pressure_Pa) - H));
    for (const auto& l : layers()) worst = std::max(worst, std::fabs(*geopotential_height_at_pressure(l.base_pressure_Pa) - l.base_geopotential_m));
    s.check("the pressure altitude inverts eqs. 33a and 33b (-5 to 84.852 km', and at every base)", worst, 1e-6, "m'");
    s.require("no pressure altitude beyond the model (0.3 Pa, 2 bar)", !geopotential_height_at_pressure(0.3) && !geopotential_height_at_pressure(2e5));
  }

  // ---- columns of levels -------------------------------------------------------------------------------------------
  {
    const auto up = profile(0, 30000, 100);
    bool falling = true;
    for (std::size_t i = 1; i < up.size(); i++) falling = falling && up[i].pressure_Pa < up[i - 1].pressure_Pa;
    s.require("a rising column: 0 to 30 km every 100 m is 301 levels, pressure falling", up.size() == 301 && falling && up.back().geometric_height_m == 30000);
    const auto down = profile(86000, 80000, 700);
    s.require("a falling column: 86 to 80 km every 700 m ends at 80 km (10 levels)", down.size() == 10 && down.back().geometric_height_m == 80000 && down[1].geometric_height_m == 85300);
    bool threw_domain = false, threw_step = false;
    try { profile(0, 90000, 1000); } catch (const std::domain_error&) { threw_domain = true; }
    try { profile(0, 1000, 0); } catch (const std::invalid_argument&) { threw_step = true; }
    s.require("a column beyond 86 km, or with no step, is refused", threw_domain && threw_step);
  }
  return s.done();
}
