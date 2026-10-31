// A building in its environment (the infrastructure's ADR 0062; the urban-science domain's
// theory/landscape-dynamics.md §7): the heat of its rooms conducted through their faces, which is the Hodge
// Laplacian Δ₃ of the building's complex with the fabric in its stars; the Sun on its facades and roof; the
// building's thermal modes; and each hour propagated exactly.
//
// Heat (physis hodge_laplacian.hpp). With Q_c = C_c T_c the heat of room c (C_c its capacity, the volume times an
// effective volumetric capacity of air and fabric) and G_f = U_f A_f the conductance of face f (U by its type, in
// W/m²K), conduction is dQ/dt = −d G dᵀ T with d the coboundary from faces to rooms, and that is exactly
// −Δ₃ Q for the stars w_f = 1/G_f and w_c = 1/C_c: the constitutive law is the Hodge star. An outer face bounds one
// room and conducts to the outside: dQ/dt = −Δ₃ Q + g T_out + S, g_c = Σ_{outer f ∋ c} G_f, S the gains. Δ₃ is
// positive definite (no harmonic 3-cochain on a ball), so each mode θ_k decays at its own rate λ_k: the
// building's time constants 1/λ_k. Each hour is propagated exactly, the forcing held over it (the exponential of
// the augmented matrix, Van Loan 1978).
//
// The Sun (a clear sky; the caller gives its altitude and azimuth, from north through east). The direct normal
// irradiance is Meinel's I = 1353 · 0.7^{m^0.678} W/m², m = 1/sin(altitude) the air mass (Meinel and Meinel 1976,
// Applied Solar Energy, §2), on a face of outward normal n the flux I max(0, s·n), s the direction of the Sun. What
// reaches the room behind it is the sol-air gain (ASHRAE Handbook — Fundamentals 2017, ch. 18): sunlight absorbed
// on opaque fabric warms its outer surface, which loses most of it to the air outside, so that through a face of
// transmittance U only U α I cos/h_o enters, h_o ≈ 25 W/m²K the outer surface's coefficient; through the glazed
// share g of a facade the glass's solar heat gain coefficient times the light. The ground faces none. Diffuse light
// and shading by neighbours are left out — a stated limit. Checks: tests/check_environment.cpp.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numbers>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/expm.hpp>
#include <physis/discretisation/hodge_laplacian.hpp>
#include <architecture/space/rooms.hpp>

namespace architecture::environment {

namespace mt = morphe::topology;
using space::Building;

// The building's conduction: Δ₃ with stars 1/G on the faces and 1/C on the rooms, and the rooms' indices.
struct Conduction {
  physis::discretisation::HodgeLaplacian H;
  std::vector<std::uint64_t> rooms;   // room of each row
  std::vector<double> outer;          // g_c = Σ_{outer f ∋ c} G_f (W/K)
};
inline Conduction conduction(const Building& B) {
  const auto C = mt::chain_complex(B.complex);
  std::vector<std::vector<double>> w(4);
  Conduction K;
  const auto up = B.complex.cofaces();
  std::map<std::uint64_t, std::size_t> row;
  for (const auto& c : B.complex.cells)
    if (c.dim == 3) row[c.id] = K.rooms.size(), K.rooms.push_back(c.id), w[3].push_back(1 / B.capacity.at(c.id));
  K.outer.assign(K.rooms.size(), 0.0);
  for (std::size_t i = 0; i < B.complex.size(); i++) {
    const auto& c = B.complex.cells[i];
    if (c.dim != 2) continue;
    const auto it = B.fabric.u.find(c.type);
    if (it == B.fabric.u.end()) throw std::invalid_argument("urban::landscape::conduction: no U value for a face of type " + c.type);
    const double G = it->second * B.area(c.id);
    w[2].push_back(1 / G);
    if (up[i].size() == 1) K.outer[row.at(B.complex.cells[up[i][0]].id)] += G;
  }
  K.H = physis::discretisation::hodge_laplacian(C, 3, w);
  return K;
}

// The heat each room gains from the Sun (W): altitude and azimuth in radians (azimuth from north through east).
inline std::map<std::uint64_t, double> solar_gains(const Building& B, double altitude, double azimuth) {
  std::map<std::uint64_t, double> S;
  for (const auto& [r, q] : B.heat) S[r] = 0;
  if (!(altitude > 0)) return S;
  const double m = 1 / std::sin(altitude), I = 1353 * std::pow(0.7, std::pow(m, 0.678));
  const std::array<double, 3> sun{std::cos(altitude) * std::sin(azimuth), std::cos(altitude) * std::cos(azimuth), std::sin(altitude)};
  const auto up = B.complex.cofaces();
  for (std::size_t i = 0; i < B.complex.size(); i++) {
    const auto& c = B.complex.cells[i];
    if (c.dim != 2 || up[i].size() != 1 || c.type == "ground") continue;
    const auto n = B.normal(c.id);
    const double cosine = sun[0] * n[0] + sun[1] * n[1] + sun[2] * n[2];
    if (!(cosine > 0)) continue;
    const auto& F = B.fabric;
    const double glass = c.type == "facade" ? F.glazing : 0, opaque = F.u.at(c.type) * F.absorptance / F.outer_coefficient;
    S[B.complex.cells[up[i][0]].id] += B.area(c.id) * I * cosine * (glass * F.gain_coefficient + (1 - glass) * opaque);
  }
  return S;
}

// Advances the rooms' heat by h seconds, the outside at t_out and the gains held: Q ← e^{−Δh} Q + Φ (g T_out + S),
// Φ = ∫₀ʰ e^{−Δs} ds, both from one exponential of the augmented matrix [[−Δ, I], [0, 0]] h.
inline void advance_heat(Building& B, const Conduction& K, double h, double t_out, const std::map<std::uint64_t, double>& gains) {
  const std::size_t n = K.rooms.size();
  calx::linalg::Matrix M(2 * n);
  for (std::size_t i = 0; i < n; i++) {
    for (std::size_t j = 0; j < n; j++) M(i, j) = -K.H.laplacian(i, j) * h;
    M(i, n + i) = h;
  }
  const auto E = calx::linalg::expm(M);
  std::vector<double> Q(n), b(n), next(n, 0.0);
  for (std::size_t i = 0; i < n; i++) {
    Q[i] = B.heat.at(K.rooms[i]);
    const auto g = gains.find(K.rooms[i]);
    b[i] = K.outer[i] * t_out + (g == gains.end() ? 0.0 : g->second);
  }
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) next[i] += E(i, j) * Q[j] + E(i, n + j) * b[j];
  for (std::size_t i = 0; i < n; i++) B.heat[K.rooms[i]] = next[i];
}

}  // namespace architecture::environment
