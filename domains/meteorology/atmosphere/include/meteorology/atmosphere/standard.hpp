// The U.S. Standard Atmosphere 1976 from −5 km′ to 86 km: temperature, pressure, density and the
// properties of the air derived from them, at a geometric or a geopotential height, and columns of
// levels for an instrument that rises or falls through it.
//
// Source: NOAA, NASA and USAF (1976), U.S. Standard Atmosphere, 1976, NOAA-S/T 76-1562,
// Washington, D.C. (NASA-TM-X-74335, NTRS 19770009539), Part 1, "Defining constants and
// equations". Equation, table and page numbers below are the document's own.
//
// The model. Below 86 km the air is a mixed ideal gas of fixed composition in hydrostatic
// equilibrium. Two changes of variable make the hydrostatic equation integrable in closed form
// (§1.2.3–1.2.5): the geopotential height H, in standard geopotential metres m′, absorbs the fall of
// gravity with height (eqs. 15–19, with the 1976's own Earth radius r₀ = 6 356 766 m and
// g₀ = 9.806 65 m s⁻², eq. 17: g = g₀ (r₀/(r₀ + Z))²); and the molecular-scale temperature
// T_M = T M₀/M absorbs the change of the mean molecular weight M (eq. 22). T_M is linear in H in
// seven layers b = 0 … 6 with the bases H_b and gradients L_M,b of Table 4 (eq. 23), from
// T_M,0 = T₀ = 288.15 K; the pressure follows by eq. 33a where L_M,b ≠ 0 and eq. 33b where the
// layer is isothermal, each layer from the pressure at its base, P₀ = 101 325 Pa and P_b the value
// of the layer below at H_b (p. 12). Density, number density, the speed of sound and the rest
// follow from T_M, T and P by eqs. 41–53.
//
// Between 80 and 86 km. The 1976 lets M fall from M₀ = 28.9644 kg kmol⁻¹ at 80 km to
// 28.9522 kg kmol⁻¹ at 86 km, by the ratios M/M₀ of Table 8, tabulated at every 0.5 geometric
// kilometre (p. 9); the kinetic temperature is then T = T_M (M/M₀) (eq. 22), and at 86 km
// T₇ = 186.8673 K (eq. 25). The document's own tables were computed without these ratios below
// 86 km (§1.2.4, p. 9: "these values of M/M₀ were not included in the program used for computing
// the tables"), and it says how to correct them: T, M and the mean free path are multiplied by
// M/M₀, the number density and the collision frequency divided by it, and the viscosities and
// conductivity recomputed from the corrected T. This module does what the document directs: it
// interpolates Table 8 linearly in Z between its nodes (the document names no interpolation) and
// computes every property from the corrected T. Pressure, density, T_M, the speed of sound, the
// mean particle speed and the pressure scale height depend on T_M and M₀ alone (eqs. 33, 42, 44,
// 46, 50) and need no correction.
//
// Above 86 km: no value. There the 1976 is another model: the kinetic temperature is given in
// four segments of geometric height (eqs. 24–31), each species N₂, O, O₂, Ar, He and H diffuses
// under its own scale height (eqs. 6–8, 35–40), and the pressure is the sum of the partial
// pressures (eq. 33c). That model is not implemented; the functions here return std::nullopt.
// Below −5 km′ there is no value either: the document extends layer 0 that far (p. 12) and no
// further.
//
// Constants (Table 2 and §1.2.4), in SI with the mole, not the 1976's kilomole:
//   g₀ = 9.806 65 m s⁻²; r₀ = 6 356 766 m; R* = 8.314 32 J mol⁻¹ K⁻¹ (the 1962 value, which the
//   1976 keeps, p. 3); M₀ = 28.9644 × 10⁻³ kg mol⁻¹ (eq. 21, from Table 3);
//   N_A = 6.022 169 × 10²³ mol⁻¹; k = 1.380 622 × 10⁻²³ J K⁻¹; P₀ = 101 325 Pa; T₀ = 288.15 K;
//   β = 1.458 × 10⁻⁶ kg s⁻¹ m⁻¹ K⁻¹ᐟ²; S = 110.4 K; γ = 1.40; σ = 3.65 × 10⁻¹⁰ m.
// Remark. Table 2 and the discussion of the constants on p. 4 print the Sutherland constant as
// "110 K"; the text of eq. 51 (p. 19) gives 110.4 K, and the tables are computed with 110.4 K
// (Table III's 1.7894 × 10⁻⁵ N s m⁻² at sea level). This module uses 110.4 K.
//
// The names of the layers. The document numbers its layers b = 0 … 6 and does not name them. They
// are named here by the thermal structure each idealises: the troposphere where T_M falls
// (b = 0, from −5 to 11 km′), the tropopause where it is constant (b = 1, 11–20 km′), the
// stratosphere where it rises (b = 2, 3, 20–47 km′), the stratopause where it is constant again
// (b = 4, 47–51 km′) and the mesosphere where it falls (b = 5, 6, 51–84.852 km′).
//
// Accuracy. The functions reproduce the document's equations to rounding; the checks hold them to
// two independent implementations (to 1e-12 where the constants are shared), to the hydrostatic
// equation, and to the printed Tables I–III at thirty-eight heights: their rounded columns to half a
// unit of the last printed digit, or one unit for a few derived columns that depart from the
// document's own equations by that much; and Table I's pressure in millibars, which is truncated, not
// rounded (226.9996 mb printed 226.99 at Z = 11 km), to between nought and one unit below the
// equations' value. The model is a standard, not a
// forecast: the real atmosphere departs from it by tens of kelvin and by tens of per cent in
// density, by latitude, season and weather (the document's Part 2 gives the variability).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace meteorology::atmosphere {

// ---- the constants (Table 2; M₀ from §1.2.4) ------------------------------------------------------------
namespace us1976 {
inline constexpr double g0_m_s2 = 9.80665;               // sea-level gravity; also g₀′ in m² s⁻² m′⁻¹
inline constexpr double r0_m = 6356766.0;                // effective Earth radius (p. 4; eq. 17)
inline constexpr double R_star_J_mol_K = 8.31432;        // gas constant (Table 2: 8.314 32 × 10³ J kmol⁻¹ K⁻¹)
inline constexpr double M0_kg_mol = 28.9644e-3;          // sea-level mean molecular weight (eq. 21)
inline constexpr double N_A_per_mol = 6.022169e23;       // Avogadro (Table 2: 6.022 169 × 10²⁶ kmol⁻¹)
inline constexpr double k_J_K = 1.380622e-23;            // Boltzmann
inline constexpr double P0_Pa = 101325.0;                // sea-level pressure
inline constexpr double T0_K = 288.15;                   // sea-level temperature
inline constexpr double beta_kg_s_m_K12 = 1.458e-6;      // Sutherland's β (eq. 51)
inline constexpr double S_K = 110.4;                     // Sutherland's constant (eq. 51; see the remark above)
inline constexpr double gamma = 1.40;                    // ratio of specific heats
inline constexpr double sigma_m = 3.65e-10;              // mean effective collision diameter
inline constexpr double bottom_geopotential_m = -5000.0; // the lowest height of the model (p. 12), m′
inline constexpr double top_geometric_m = 86000.0;       // the top of the lower model (Z₇)
inline constexpr double H7_geopotential_m = 84852.0;     // Table 4: H₇ = 84.8520 km′
}  // namespace us1976

// ---- geopotential and geometric height, gravity (eqs. 17–19) ---------------------------------------------
// The geopotential height H (m′) of a geometric height Z (m), eq. 18 with Γ = g₀/g₀′ = 1 m′ m⁻¹.
inline double geopotential_height(double geometric_height_m) {
  return us1976::r0_m * geometric_height_m / (us1976::r0_m + geometric_height_m);
}
// The geometric height Z (m) of a geopotential height H (m′), eq. 19.
inline double geometric_height(double geopotential_height_m) {
  return us1976::r0_m * geopotential_height_m / (us1976::r0_m - geopotential_height_m);
}
// The 1976's gravity g (m s⁻²) at a geometric height, eq. 17.
inline double gravity(double geometric_height_m) {
  const double s = us1976::r0_m / (us1976::r0_m + geometric_height_m);
  return us1976::g0_m_s2 * s * s;
}

// ---- the seven layers (Table 4; base pressures by eqs. 33a/33b, p. 12) -----------------------------------
struct Layer {
  int b = 0;                         // the subscript of Table 4
  std::string_view name;             // its name, by thermal structure (see above)
  double base_geopotential_m = 0;    // H_b (m′)
  double top_geopotential_m = 0;     // H_{b+1} (m′)
  double gradient_K_per_m = 0;       // L_M,b (K m′⁻¹)
  double base_temperature_K = 0;     // T_M,b
  double base_pressure_Pa = 0;       // P_b
};

namespace detail {
inline double exponent(double gradient_K_per_m) {
  return us1976::g0_m_s2 * us1976::M0_kg_mol / (us1976::R_star_J_mol_K * gradient_K_per_m);
}
// Eqs. 33a and 33b in a layer of base H_b, T_M,b, P_b and gradient L_M,b, at H.
inline double pressure(double base_H, double base_T, double base_P, double L, double H) {
  if (L == 0) return base_P * std::exp(-us1976::g0_m_s2 * us1976::M0_kg_mol * (H - base_H) / (us1976::R_star_J_mol_K * base_T));
  return base_P * std::pow(base_T / (base_T + L * (H - base_H)), exponent(L));
}
}  // namespace detail

// The layers b = 0 … 6, with their base temperatures and pressures carried up from T₀ and P₀.
inline const std::array<Layer, 7>& layers() {
  static const std::array<Layer, 7> table = [] {
    constexpr std::array<double, 8> H{0, 11000, 20000, 32000, 47000, 51000, 71000, us1976::H7_geopotential_m};
    constexpr std::array<double, 7> L{-6.5e-3, 0.0, 1.0e-3, 2.8e-3, 0.0, -2.8e-3, -2.0e-3};
    constexpr std::array<std::string_view, 7> names{"troposphere", "tropopause", "stratosphere", "stratosphere",
                                                    "stratopause", "mesosphere", "mesosphere"};
    std::array<Layer, 7> t{};
    double T = us1976::T0_K, P = us1976::P0_Pa;
    for (int b = 0; b < 7; b++) {
      t[b] = {b, names[b], H[b], H[b + 1], L[b], T, P};
      P = detail::pressure(H[b], T, P, L[b], H[b + 1]);
      T = T + L[b] * (H[b + 1] - H[b]);
    }
    return t;
  }();
  return table;
}

// The highest geopotential height of the model, H(86 km) = 84 852.05 m′ (5 cm above H₇; layer 6
// is carried across that step).
inline double top_geopotential_m() { return geopotential_height(us1976::top_geometric_m); }

// The layer b in which a geopotential height lies (H_b ≤ H < H_{b+1}, layer 6 to the top), or
// none outside [−5 km′, H(86 km)].
inline std::optional<int> layer_index(double geopotential_height_m) {
  const double H = geopotential_height_m;
  if (!(H >= us1976::bottom_geopotential_m && H <= top_geopotential_m())) return std::nullopt;
  const auto& t = layers();
  int b = 0;
  while (b < 6 && H >= t[b + 1].base_geopotential_m) b++;
  return b;
}

// A layer's own formulas, carried beyond the layer if asked (for the checks of continuity):
// T_M by eq. 23 and P by eq. 33a or 33b.
inline double molecular_temperature_in_layer(int b, double geopotential_height_m) {
  const Layer& l = layers().at(static_cast<std::size_t>(b));
  return l.base_temperature_K + l.gradient_K_per_m * (geopotential_height_m - l.base_geopotential_m);
}
inline double pressure_in_layer(int b, double geopotential_height_m) {
  const Layer& l = layers().at(static_cast<std::size_t>(b));
  return detail::pressure(l.base_geopotential_m, l.base_temperature_K, l.base_pressure_Pa, l.gradient_K_per_m, geopotential_height_m);
}

// ---- the molecular weight between 80 and 86 km (Table 8) --------------------------------------------------
// M/M₀ at Z = 80.0, 80.5, …, 86.0 km, from the right-hand half of Table 8 (p. 9).
inline constexpr std::array<double, 13> molecular_weight_ratios{1.000000, 0.999996, 0.999989, 0.999971, 0.999941, 0.999909, 0.999870,
                                                                0.999829, 0.999786, 0.999741, 0.999694, 0.999641, 0.999579};

// M/M₀ at a geometric height: 1 below 80 km, Table 8 interpolated linearly in Z from 80 to 86 km.
inline double molecular_weight_ratio(double geometric_height_m) {
  if (geometric_height_m <= 80000.0) return 1.0;
  const double x = std::min((geometric_height_m - 80000.0) / 500.0, 12.0);
  const int i = std::min(11, static_cast<int>(x));
  return molecular_weight_ratios[i] + (molecular_weight_ratios[i + 1] - molecular_weight_ratios[i]) * (x - i);
}

// ---- the state of the air at a height ----------------------------------------------------------------------
struct State {
  double geometric_height_m = 0;           // Z
  double geopotential_height_m = 0;        // H (m′)
  int layer = 0;                           // b of Table 4
  std::string_view layer_name;
  double temperature_K = 0;                // kinetic temperature T = T_M M/M₀ (eq. 22)
  double molecular_temperature_K = 0;      // T_M (eq. 23)
  double molecular_weight_ratio = 1;       // M/M₀ (Table 8)
  double molecular_weight_kg_mol = 0;      // M
  double pressure_Pa = 0;                  // P (eqs. 33a, 33b)
  double density_kg_m3 = 0;                // ρ = P M₀/(R* T_M) (eq. 42)
  double gravity_m_s2 = 0;                 // g (eq. 17)
  double pressure_scale_height_m = 0;      // H_P = R* T_M/(g M₀) (eq. 44)
  double number_density_per_m3 = 0;        // N = N_A P/(R* T) (eq. 41)
  double mole_volume_m3_mol = 0;           // v_m = R* T/P (eq. 43)
  double mean_particle_speed_m_s = 0;      // V = (8 R* T_M/(π M₀))^½ (eq. 46)
  double mean_free_path_m = 0;             // L = √2 R* T/(2π N_A σ² P) (eq. 47)
  double collision_frequency_per_s = 0;    // ν = V/L (eqs. 48, 49)
  double speed_of_sound_m_s = 0;           // C_s = (γ R* T_M/M₀)^½ (eq. 50)
  double dynamic_viscosity_Pa_s = 0;       // μ = β T^{3/2}/(T + S) (eq. 51)
  double kinematic_viscosity_m2_s = 0;     // η = μ/ρ (eq. 52)
  double thermal_conductivity_W_m_K = 0;   // k_t = 2.646 38 × 10⁻³ T^{3/2}/(T + 245.4 × 10^{−12/T}) (eq. 53)
};

namespace detail {
inline State state(double Z, double H, int b) {
  using namespace us1976;
  constexpr double pi = 3.14159265358979323846;
  State s;
  s.geometric_height_m = Z;
  s.geopotential_height_m = H;
  s.layer = b;
  s.layer_name = layers()[b].name;
  s.molecular_temperature_K = molecular_temperature_in_layer(b, H);
  s.molecular_weight_ratio = molecular_weight_ratio(Z);
  s.molecular_weight_kg_mol = M0_kg_mol * s.molecular_weight_ratio;
  s.temperature_K = s.molecular_temperature_K * s.molecular_weight_ratio;
  s.pressure_Pa = pressure_in_layer(b, H);
  const double TM = s.molecular_temperature_K, T = s.temperature_K, P = s.pressure_Pa;
  s.density_kg_m3 = P * M0_kg_mol / (R_star_J_mol_K * TM);
  s.gravity_m_s2 = gravity(Z);
  s.pressure_scale_height_m = R_star_J_mol_K * TM / (s.gravity_m_s2 * M0_kg_mol);
  s.number_density_per_m3 = N_A_per_mol * P / (R_star_J_mol_K * T);
  s.mole_volume_m3_mol = R_star_J_mol_K * T / P;
  s.mean_particle_speed_m_s = std::sqrt(8 * R_star_J_mol_K * TM / (pi * M0_kg_mol));
  s.mean_free_path_m = std::sqrt(2.0) * R_star_J_mol_K * T / (2 * pi * N_A_per_mol * sigma_m * sigma_m * P);
  s.collision_frequency_per_s = s.mean_particle_speed_m_s / s.mean_free_path_m;
  s.speed_of_sound_m_s = std::sqrt(gamma * R_star_J_mol_K * TM / M0_kg_mol);
  s.dynamic_viscosity_Pa_s = beta_kg_s_m_K12 * T * std::sqrt(T) / (T + S_K);
  s.kinematic_viscosity_m2_s = s.dynamic_viscosity_Pa_s / s.density_kg_m3;
  s.thermal_conductivity_W_m_K = 2.64638e-3 * T * std::sqrt(T) / (T + 245.4 * std::pow(10.0, -12.0 / T));
  return s;
}
}  // namespace detail

// The air at a geometric height Z (m), from Z(−5 km′) = r₀H/(r₀ − H) = −4 996.1 m to 86 km; none outside.
inline std::optional<State> at_geometric_height(double geometric_height_m) {
  const double Z = geometric_height_m;
  if (!(Z <= us1976::top_geometric_m && Z > -us1976::r0_m)) return std::nullopt;
  const double H = geopotential_height(Z);
  const auto b = layer_index(H);
  if (!b) return std::nullopt;
  return detail::state(Z, H, *b);
}

// The air at a geopotential height H (m′), from −5 km′ to H(86 km); none outside.
inline std::optional<State> at_geopotential_height(double geopotential_height_m) {
  const auto b = layer_index(geopotential_height_m);
  if (!b) return std::nullopt;
  return detail::state(std::min(geometric_height(geopotential_height_m), us1976::top_geometric_m), geopotential_height_m, *b);
}

// The pressure altitude: the geopotential height (m′) at which the standard pressure is p (Pa),
// by eqs. 33a and 33b solved for H in the layer whose pressures bracket p; none for pressures
// outside [P(86 km), P(−5 km′)].
inline std::optional<double> geopotential_height_at_pressure(double pressure_Pa) {
  const double p = pressure_Pa;
  const double p_top = pressure_in_layer(6, top_geopotential_m()), p_bottom = pressure_in_layer(0, us1976::bottom_geopotential_m);
  if (!(p >= p_top && p <= p_bottom)) return std::nullopt;
  const auto& t = layers();
  int b = 0;
  while (b < 6 && p <= t[b + 1].base_pressure_Pa) b++;
  const Layer& l = t[b];
  if (l.gradient_K_per_m == 0)
    return l.base_geopotential_m - us1976::R_star_J_mol_K * l.base_temperature_K / (us1976::g0_m_s2 * us1976::M0_kg_mol) * std::log(p / l.base_pressure_Pa);
  return l.base_geopotential_m + l.base_temperature_K / l.gradient_K_per_m * (std::pow(p / l.base_pressure_Pa, -1 / detail::exponent(l.gradient_K_per_m)) - 1);
}

// A column of levels, for an instrument that rises or falls through the standard atmosphere: the
// air at geometric heights from `from` to `to` (m), every `step` metres (step > 0, in the
// direction from `from` to `to`), with `to` itself as the last level. Throws std::invalid_argument
// for a step that is not positive and finite, and std::domain_error when either end lies outside
// the model (above 86 km or below −5 km′).
inline std::vector<State> profile(double from_geometric_m, double to_geometric_m, double step_m) {
  if (!(step_m > 0) || !std::isfinite(step_m)) throw std::invalid_argument("meteorology: a profile needs a positive, finite step");
  const auto first = at_geometric_height(from_geometric_m), last = at_geometric_height(to_geometric_m);
  if (!first || !last) throw std::domain_error("meteorology: the U.S. Standard Atmosphere 1976 is implemented from -5 km' to 86 km only");
  const double span = to_geometric_m - from_geometric_m, direction = span >= 0 ? 1.0 : -1.0;
  const auto steps = static_cast<long long>(std::floor(std::fabs(span) / step_m * (1 + 1e-12)));
  std::vector<State> column;
  column.reserve(static_cast<std::size_t>(steps) + 2);
  for (long long i = 0; i <= steps; i++) {
    double z = from_geometric_m + direction * step_m * static_cast<double>(i);
    if (direction * (z - to_geometric_m) > 0) z = to_geometric_m;  // rounding never steps past the end
    column.push_back(*at_geometric_height(z));
  }
  if (std::fabs(column.back().geometric_height_m - to_geometric_m) > 1e-9 * std::max(1.0, std::fabs(to_geometric_m))) column.push_back(*last);
  return column;
}

}  // namespace meteorology::atmosphere
