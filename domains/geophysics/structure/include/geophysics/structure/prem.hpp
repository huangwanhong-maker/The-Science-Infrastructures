// The Preliminary Reference Earth Model (PREM): the Earth as a sphere of shells, its density and the
// speeds of P and S waves as polynomials of the radius in each (Dziewonski and Anderson 1981, Table I).
//
// PREM was fitted to the Earth's free oscillations, to the travel times of body waves, and to its mass
// and moment of inertia. This is its isotropic form, at the reference period of 1 s, with its ocean of
// 3 km: values in shell k are ρ, vp, vs = Σ c_j x^j with x = r / a, a = 6371 km, in g/cm³ and km/s
// (here converted to SI). At a discontinuity a query takes the shell below unless asked for the one
// above. The 80 km of the anisotropic model's lid and low-velocity zone has no jump in the isotropic
// one, whose two shells share their polynomials.
//
// The mass within a radius is the polynomial density integrated exactly, ∫ 4π ρ(r) r² dr; gravity is
// G M(r) / r², which is exact for a spherically symmetric Earth. Limits: PREM is one radial profile for
// the whole Earth; it has no lateral variation, its crust is a global average of two layers, and its
// values are reference values at 1 s (the anelastic dispersion of other periods is not applied).

#pragma once

#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace geophysics::structure {

// PREM's radius (m) and the gravitational constant (CODATA 2018, m³ kg⁻¹ s⁻²).
inline constexpr double PREM_RADIUS_M = 6371.0e3;
inline constexpr double GRAVITATIONAL_CONSTANT = 6.67430e-11;

// One shell: its radii (km), name, region, and the polynomials in x = r / a of density (g/cm³) and of
// the speeds of P and S waves (km/s).
struct PremShell {
  double inner_km, outer_km;
  std::string_view name, region;
  std::array<double, 4> rho, vp, vs;
};

// Dziewonski and Anderson (1981), Table I, isotropic, from the centre outward.
inline constexpr std::array<PremShell, 13> PREM_SHELLS{{
    {0.0, 1221.5, "inner core", "inner core", {13.0885, 0, -8.8381, 0}, {11.2622, 0, -6.3640, 0}, {3.6678, 0, -4.4475, 0}},
    {1221.5, 3480.0, "outer core", "outer core", {12.5815, -1.2638, -3.6426, -5.5281}, {11.0487, -4.0362, 4.8023, -13.5732}, {0, 0, 0, 0}},
    {3480.0, 3630.0, "D″", "lower mantle", {7.9565, -6.4761, 5.5283, -3.0807}, {15.3891, -5.3181, 5.5242, -2.5514}, {6.9254, 1.4672, -2.0834, 0.9783}},
    {3630.0, 5600.0, "lower mantle", "lower mantle", {7.9565, -6.4761, 5.5283, -3.0807}, {24.9520, -40.4673, 51.4832, -26.6419}, {11.1671, -13.7818, 17.4575, -9.2777}},
    {5600.0, 5701.0, "top of the lower mantle", "lower mantle", {7.9565, -6.4761, 5.5283, -3.0807}, {29.2766, -23.6027, 5.5242, -2.5514}, {22.3459, -17.2473, -2.0834, 0.9783}},
    {5701.0, 5771.0, "transition zone (600–670 km)", "transition zone", {5.3197, -1.4836, 0, 0}, {19.0957, -9.8672, 0, 0}, {9.9839, -4.9324, 0, 0}},
    {5771.0, 5971.0, "transition zone (400–600 km)", "transition zone", {11.2494, -8.0298, 0, 0}, {39.7027, -32.6166, 0, 0}, {22.3512, -18.5856, 0, 0}},
    {5971.0, 6151.0, "upper mantle (220–400 km)", "upper mantle", {7.1089, -3.8045, 0, 0}, {20.3926, -12.2569, 0, 0}, {8.9496, -4.4597, 0, 0}},
    {6151.0, 6291.0, "low-velocity zone", "upper mantle", {2.6910, 0.6924, 0, 0}, {4.1875, 3.9382, 0, 0}, {2.1519, 2.3481, 0, 0}},
    {6291.0, 6346.6, "lid", "upper mantle", {2.6910, 0.6924, 0, 0}, {4.1875, 3.9382, 0, 0}, {2.1519, 2.3481, 0, 0}},
    {6346.6, 6356.0, "lower crust", "crust", {2.900, 0, 0, 0}, {6.800, 0, 0, 0}, {3.900, 0, 0, 0}},
    {6356.0, 6368.0, "upper crust", "crust", {2.600, 0, 0, 0}, {5.800, 0, 0, 0}, {3.200, 0, 0, 0}},
    {6368.0, 6371.0, "ocean", "ocean", {1.020, 0, 0, 0}, {1.450, 0, 0, 0}, {0, 0, 0, 0}},
}};

// The discontinuities of PREM, by depth below its surface (m), with their usual names.
struct PremDiscontinuity {
  double depth_m;
  std::string_view name;
};
inline constexpr std::array<PremDiscontinuity, 8> PREM_DISCONTINUITIES{{
    {3.0e3, "the sea floor"},
    {15.0e3, "the middle crust"},
    {24.4e3, "the Moho"},
    {220.0e3, "the 220 km discontinuity"},
    {400.0e3, "the 400 km discontinuity"},
    {670.0e3, "the 670 km discontinuity"},
    {2891.0e3, "the core–mantle boundary"},
    {5149.5e3, "the inner core boundary"},
}};

// What PREM gives at a radius: SI values, and the shell's name and region.
struct PremValues {
  double density_kg_m3 = 0, vp_m_s = 0, vs_m_s = 0;
  std::string_view shell, region;
  int index = -1;
};

namespace prem_detail {
inline double poly(const std::array<double, 4>& c, double x) { return c[0] + x * (c[1] + x * (c[2] + x * c[3])); }
}  // namespace prem_detail

// The index of the shell holding radius r (m), 0 ≤ r ≤ a; at a boundary the shell below, or with
// `above` the one above.
inline int prem_shell_index(double r_m, bool above = false) {
  if (!(r_m >= 0 && r_m <= PREM_RADIUS_M)) throw std::invalid_argument("geophysics::structure: PREM holds radii from 0 to 6371 km");
  const double r = r_m / 1e3;
  for (int k = 0; k < int(PREM_SHELLS.size()); k++) {
    const auto& s = PREM_SHELLS[k];
    if (above ? (r >= s.inner_km && r < s.outer_km) : (r > s.inner_km && r <= s.outer_km)) return k;
  }
  return above ? int(PREM_SHELLS.size()) - 1 : 0;
}

inline PremValues prem_at_radius(double r_m, bool above = false) {
  const int k = prem_shell_index(r_m, above);
  const auto& s = PREM_SHELLS[k];
  const double x = r_m / PREM_RADIUS_M;
  return {prem_detail::poly(s.rho, x) * 1e3, prem_detail::poly(s.vp, x) * 1e3, prem_detail::poly(s.vs, x) * 1e3, s.name, s.region, k};
}

// At a depth below PREM's surface (m), 0 to 6371 km; at a discontinuity the side below, or with `above`
// the side above (toward the surface).
inline PremValues prem_at_depth(double depth_m, bool above = false) { return prem_at_radius(PREM_RADIUS_M - depth_m, above); }

// The mass within radius r (kg): Σ over shells of 4π a³ Σ_j c_j (x₁^{j+3} − x₀^{j+3}) / (j + 3), exact.
inline double prem_mass_within(double r_m) {
  if (!(r_m >= 0)) throw std::invalid_argument("geophysics::structure: a radius of at least 0");
  const double a = PREM_RADIUS_M, xr = std::min(r_m, a) / a;
  double m = 0;
  for (const auto& s : PREM_SHELLS) {
    const double x0 = s.inner_km * 1e3 / a, x1 = std::min(s.outer_km * 1e3 / a, xr);
    if (x1 <= x0) break;
    for (int j = 0; j < 4; j++) m += s.rho[j] * 1e3 * (std::pow(x1, j + 3) - std::pow(x0, j + 3)) / (j + 3);
  }
  return 4 * std::numbers::pi * a * a * a * m;
}

// Gravity at radius r (m/s²), G M(r) / r², for the spherically symmetric Earth; 0 at the centre.
inline double prem_gravity(double r_m) { return r_m <= 0 ? 0.0 : GRAVITATIONAL_CONSTANT * prem_mass_within(r_m) / (r_m * r_m); }

// The moment of inertia about an axis, (8π/3) ∫ ρ r⁴ dr, exact (kg m²).
inline double prem_moment_of_inertia() {
  const double a = PREM_RADIUS_M;
  double m = 0;
  for (const auto& s : PREM_SHELLS) {
    const double x0 = s.inner_km * 1e3 / a, x1 = s.outer_km * 1e3 / a;
    for (int j = 0; j < 4; j++) m += s.rho[j] * 1e3 * (std::pow(x1, j + 5) - std::pow(x0, j + 5)) / (j + 5);
  }
  return 8 * std::numbers::pi / 3 * std::pow(a, 5) * m;
}

}  // namespace geophysics::structure
