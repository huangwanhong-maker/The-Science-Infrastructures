// What the fluids family shares: the properties of its liquids, and the provenance of its simulators
// (substrates/fluids/README.md; ADR 0076).
//
// Properties are declared at 20 °C and 1 atm from standard tables (water: the CRC Handbook of Chemistry and
// Physics, 97th ed., §6; a silicone oil of 20 cSt: its makers' sheets, the walkers' liquid of Couder and Fort
// 2006), with their sources, and may be replaced by a caller's measured values.
#pragma once

#include <cstdint>
#include <string>

#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::fluids {

inline constexpr double standard_gravity = 9.80665;  // m/s², by definition (CGPM 1901)

struct Fluid {
  std::string name;
  double density = 0;              // ρ, kg/m³
  double kinematic_viscosity = 0;  // ν, m²/s
  double thermal_diffusivity = 0;  // κ, m²/s
  double thermal_expansion = 0;    // α, 1/K
  double surface_tension = 0;      // σ, N/m
  std::string source;

  static Fluid water() {
    return {"water at 20 °C", 998.2, 1.004e-6, 1.43e-7, 2.07e-4, 0.0728, "CRC Handbook of Chemistry and Physics, 97th ed., §6 (declared)"};
  }
  static Fluid silicone_oil_20cSt() {
    return {"silicone oil, 20 cSt, at 20 °C", 950.0, 20e-6, 1.0e-7, 1.07e-3, 0.0206, "makers' sheets for polydimethylsiloxane of 20 cSt (declared)"};
  }
};

inline runtime::Provenance simulator_provenance(std::string substrate, std::string model, std::uint64_t seed = 0) {
  runtime::Provenance p;
  p.kind = runtime::BackendKind::simulator;
  p.family = "fluids", p.substrate = std::move(substrate), p.model = std::move(model), p.version = "0.1 (ADR 0076)", p.seed = seed;
  return p;
}

}  // namespace physical_computing::fluids
