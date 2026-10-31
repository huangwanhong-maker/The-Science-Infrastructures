// Checks of the environment module (heat.hpp): a building of 4 × 2 rooms on 3 floors.
// - Δ₃ positive definite: no harmonic 3-cochain, every mode decays; with no Sun every room relaxes to the outside;
// - the Sun at 50° gives gains, a fraction of the light on the outer faces (the sol-air gain); an hour's change of
//   heat equals the flux through the outer faces and the gains integrated over it, propagated in one step or in
//   3600 alike;
// - the noon Sun, due north at Melbourne's latitude, warms the north facade and not the south;
// - a stair opened raises the conduction: the slowest mode is faster.

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

#include <architecture/environment/heat.hpp>
#include <architecture/space/rooms.hpp>
#include <calx/verification/check.hpp>

using calx::verification::Suite;
namespace as = architecture::space;
namespace ae = architecture::environment;

int main() {
  Suite s;
  auto Bd = as::building({{0, 0}, 0, 20, 12, 3.5, 4, 2, 3});
  auto K = ae::conduction(Bd);
  double lo = 1e300;
  for (double v : K.H.values) lo = std::min(lo, v);
  s.require("Δ₃ positive definite: no harmonic 3-cochain, every mode decays", lo > 0 && K.H.harmonic() == 0);
  ae::advance_heat(Bd, K, 3600.0 * 2000, 15, {});  // fifty times the slowest time constant, in one exact step
  double worst = 0;
  for (const auto& [room, q] : Bd.heat) worst = std::max(worst, std::fabs(Bd.temperature(room) - 15));
  s.check("with no Sun every room relaxes to the outside", worst, 1e-6, "K");
  // An hour with the Sun: the change of heat against the integrated flux, by fine sub-steps.
  const double alt = 50 * std::numbers::pi / 180, az = 20 * std::numbers::pi / 180;
  const auto gains = ae::solar_gains(Bd, alt, az);
  double sun = 0;
  for (const auto& [r, g] : gains) sun += g;
  s.require("the Sun at 50° gives gains", sun > 1000);
  {
    // Never more than the light falling on the outer faces, and only a few per cent of it through opaque fabric.
    double light = 0;
    const auto up = Bd.complex.cofaces();
    const double I = 1353 * std::pow(0.7, std::pow(1 / std::sin(alt), 0.678));
    for (std::size_t i = 0; i < Bd.complex.size(); i++) {
      const auto& c = Bd.complex.cells[i];
      if (c.dim != 2 || up[i].size() != 1 || c.type == "ground") continue;
      const auto n = Bd.normal(c.id);
      light += Bd.area(c.id) * I * std::max(0.0, std::cos(alt) * std::sin(az) * n[0] + std::cos(alt) * std::cos(az) * n[1] + std::sin(alt) * n[2]);
    }
    s.require("the sol-air gain: a fraction of the light on the outer faces (between 2% and 30%)", sun > 0.02 * light && sun < 0.3 * light, std::to_string(sun / light));
  }
  const double Q0 = Bd.total_heat();
  double flux = 0;
  auto Bs = Bd;
  const int sub = 3600;
  const auto outer_flux = [&](const as::Building& X) {
    double f = 0;
    for (std::size_t i = 0; i < K.rooms.size(); i++) f += K.outer[i] * (15 - X.temperature(K.rooms[i]));
    return f;
  };
  for (int q = 0; q < sub; q++) {
    const double f0 = outer_flux(Bs);
    ae::advance_heat(Bs, K, 3600.0 / sub, 15, gains);
    flux += 0.5 * (f0 + outer_flux(Bs)) * 3600.0 / sub;
  }
  ae::advance_heat(Bd, K, 3600, 15, gains);
  s.check("an hour's change of heat is the flux through the outer faces and the gains (relative)", (Bd.total_heat() - Q0 - flux - sun * 3600) / (sun * 3600), 1e-6);
  s.check("propagated in one step or in 3600 (relative)", (Bd.total_heat() - Bs.total_heat()) / (sun * 3600), 1e-9);
  // The noon Sun, due north in Melbourne: on the north facade, not on the south.
  const auto noon = ae::solar_gains(Bd, 52 * std::numbers::pi / 180, 0);
  double north = 0, south = 0;
  for (const auto& [room, g] : noon) {
    const auto u = Bd.units.at(room)[0];
    if (u[2] != 1) continue;  // the middle floor: no roof
    if (u[1] == 1) north += g;
    if (u[1] == 0) south += g;
  }
  s.require("the noon Sun in the southern hemisphere warms the north facade, not the south", north > 0 && south == 0);
  // A stair opened raises the conduction.
  const std::uint64_t slab = [&] { for (const auto& c : Bd.complex.cells) if (c.type == "slab") return c.id; return std::uint64_t{0}; }();
  as::relabel_face(Bd, slab, "stair");
  const auto K2 = ae::conduction(Bd);
  s.require("a stair opened raises the conduction: the slowest mode faster", K2.H.values.front() >= K.H.values.front() - 1e-12);
  return s.done();
}
