// The heat of a district's buildings (the infrastructure's ADR 0063): each building a 3-cell of the massing
// (space/massing.hpp) holding heat Q = C T, conducted through the walls it shares with its neighbours and lost
// through its facades, roof and underside to the air and through its ground face to the ground; warmed by the Sun
// on its facades and roof. On the massing, oriented outward, this is the Hodge Laplacian Δ₃ with stars 1/G on the
// faces and 1/C on the buildings, as for a building's rooms (heat.hpp): C Ṫ = −Δ T + g_air T_air + g_ground T_ground
// + S, with Δ the weighted Laplacian of the buildings through their shared faces plus g_air + g_ground on the
// diagonal. Each step is Crank–Nicolson (calx's θ-method), conjugate gradients on the sparse rows.
//
// Conductances U A by the face: a wall shared by two buildings "party", a wall with one "facade", a roof "roof", a
// ground face at the ground "ground" (to the ground's temperature), a ground face raised above it (an underside)
// "roof" to the air; capacities the volume times an effective volumetric capacity. The Sun's gains are sol-air, as
// in heat.hpp: through glass its gain coefficient times the light, through opaque fabric U α I cos θ / h_o.
// Checks: tests/check_environment.cpp.

#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <span>
#include <stdexcept>
#include <vector>

#include <architecture/space/massing.hpp>
#include <architecture/space/rooms.hpp>
#include <calx/integrators/theta.hpp>
#include <calx/linalg/sparse.hpp>

namespace architecture::environment {

struct MassingHeat {
  std::size_t n = 0;
  std::vector<double> capacity;                                   // J/K, by building
  std::vector<double> to_air, to_ground;                          // W/K
  struct Party { std::size_t a = 0, b = 0; double G = 0; };       // W/K through the walls two buildings share
  std::vector<Party> party;
  calx::linalg::Csr<std::uint32_t> A;
  std::vector<double> a;                                          // A's values: the Laplacian plus the losses
};

inline MassingHeat massing_heat(const space::Massing& M, const space::Fabric& fabric = {}, double capacity_per_m3 = 6.0e4) {
  MassingHeat H;
  const auto& X = M.complex;
  H.n = M.cell.size();
  H.to_air.assign(H.n, 0.0), H.to_ground.assign(H.n, 0.0);
  for (double v : M.volume) H.capacity.push_back(capacity_per_m3 * v);
  const auto U = [&](const char* k) { const auto it = fabric.u.find(k); return it == fabric.u.end() ? 1.0 : it->second; };
  const auto up = X.cofaces();
  std::map<std::pair<std::size_t, std::size_t>, double> shared;
  for (std::size_t f = 0; f < X.size(); f++) {
    const auto& c = X.cells[f];
    if (c.dim != 2) continue;
    const auto& face = M.face.at(c.id);
    std::vector<std::size_t> owners;
    for (auto k : up[f]) owners.push_back(M.of_cell.at(X.cells[k].id));
    if (owners.size() == 2) {
      const double G = (face.kind == "wall" ? U("party") : U("slab")) * face.area;
      shared[{std::min(owners[0], owners[1]), std::max(owners[0], owners[1])}] += G;
    } else if (owners.size() == 1) {
      const std::size_t i = owners[0];
      if (face.kind == "wall") H.to_air[i] += U("facade") * face.area;
      else if (face.kind == "roof") H.to_air[i] += U("roof") * face.area;
      else if (face.z0 <= 1e-9) H.to_ground[i] += U("ground") * face.area;
      else H.to_air[i] += U("roof") * face.area;  // an underside, raised above the ground
    }
  }
  std::vector<std::uint32_t> ri, ci;
  std::vector<double> av;
  std::vector<double> diagonal(H.n, 0.0);
  for (const auto& [ij, G] : shared) {
    H.party.push_back({ij.first, ij.second, G});
    for (const auto& [p, q] : {std::pair{ij.first, ij.second}, std::pair{ij.second, ij.first}})
      ri.push_back(static_cast<std::uint32_t>(p)), ci.push_back(static_cast<std::uint32_t>(q)), av.push_back(-G);
    diagonal[ij.first] += G, diagonal[ij.second] += G;
  }
  for (std::size_t i = 0; i < H.n; i++)
    ri.push_back(static_cast<std::uint32_t>(i)), ci.push_back(static_cast<std::uint32_t>(i)), av.push_back(diagonal[i] + H.to_air[i] + H.to_ground[i]);
  std::vector<std::uint64_t> order;
  H.A = calx::linalg::compress<std::uint32_t>(H.n, H.n, ri, ci, &order);
  H.a.resize(order.size());
  for (std::size_t k = 0; k < order.size(); k++) H.a[k] = av[order[k]];
  return H;
}

// The Sun's gains (W) on each building at an altitude and azimuth (radians, azimuth from north through east).
inline std::vector<double> massing_gains(const space::Massing& M, double altitude, double azimuth, const space::Fabric& fabric = {}) {
  std::vector<double> S(M.cell.size(), 0.0);
  if (!(altitude > 0)) return S;
  const double I = 1353 * std::pow(0.7, std::pow(1 / std::sin(altitude), 0.678));
  const std::array<double, 3> sun{std::cos(altitude) * std::sin(azimuth), std::cos(altitude) * std::cos(azimuth), std::sin(altitude)};
  const auto& X = M.complex;
  const auto up = X.cofaces();
  const auto U = [&](const std::string& k) { const auto it = fabric.u.find(k); return it == fabric.u.end() ? 1.0 : it->second; };
  for (std::size_t b = 0; b < M.cell.size(); b++) {
    const auto& c = X.cells[X.find(M.cell[b])];
    for (std::size_t k = 0; k < c.facets.size(); k++) {
      const std::size_t f = c.facets[k];
      if (up[f].size() != 1) continue;  // shared with a neighbour: no sky
      const auto& face = M.face.at(X.cells[f].id);
      if (face.kind == "ground") continue;
      const auto& n = face.normal;
      const double s = c.signs[k];
      const double cosine = s * (sun[0] * n[0] + sun[1] * n[1] + sun[2] * n[2]);
      if (!(cosine > 0)) continue;
      const bool facade = face.kind == "wall";
      const double glass = facade ? fabric.glazing : 0, opaque = U(facade ? "facade" : "roof") * fabric.absorptance / fabric.outer_coefficient;
      S[b] += face.area * I * cosine * (glass * fabric.gain_coefficient + (1 - glass) * opaque);
    }
  }
  return S;
}

// Temperatures advanced by h seconds, the air, the ground and the gains held over the step.
inline std::vector<double> massing_step(const MassingHeat& H, const std::vector<double>& T, double air, double ground, const std::vector<double>& gains, double h) {
  if (T.size() != H.n || gains.size() != H.n) throw std::invalid_argument("architecture::environment::massing_step: a temperature and a gain for each building");
  std::vector<double> c(H.n);
  for (std::size_t i = 0; i < H.n; i++) c[i] = H.to_air[i] * air + H.to_ground[i] * ground + gains[i];
  return calx::integrators::theta_step(H.A, std::span<const double>(H.a), std::span<const double>(H.capacity), std::span<const double>(T), std::span<const double>(c), h, 0.5).y;
}

}  // namespace architecture::environment
