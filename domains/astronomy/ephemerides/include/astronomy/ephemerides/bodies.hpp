// The bodies of the Solar System as the ephemerides name them: NAIF identifiers, gravitational
// parameters and mean radii.
//
// Identifiers: NAIF, "NAIF Integer ID codes" (0 the solar-system barycentre, 1–9 the barycentres
// of the planetary systems, 10 the Sun, 199 Mercury, 299 Venus, 399 the Earth, 301 the Moon,
// 499 Mars, ...). DE440 holds the planets beyond Mars as the barycentres of their systems.
// GM: Park et al. 2021 (DE440), in km³ s⁻², here m³ s⁻². Mean radii: Archinal et al. 2018 (IAU
// WGCCRE 2015), here in metres. The astronomy domain's JavaScript (bodies.js) holds the same.

#pragma once

#include <array>
#include <optional>
#include <string_view>

namespace astronomy::ephemerides {

struct Body {
  std::string_view key, name;
  int naif;           // the body's own identifier
  int ephemeris;      // the identifier DE440 gives it (a barycentre for the outer planets)
  double gm;          // m³ s⁻²
  double radius;      // mean radius, m
  int parent;         // NAIF id of the body it orbits (0: the barycentre)
};

inline constexpr std::array<Body, 11> BODIES{{
    {"sun", "Sun", 10, 10, 132712440041.279419e9, 695700e3, 0},
    {"mercury", "Mercury", 199, 1, 22031.868551e9, 2439.4e3, 10},
    {"venus", "Venus", 299, 2, 324858.592e9, 6051.8e3, 10},
    {"earth", "Earth", 399, 399, 398600.435507e9, 6371.0084e3, 10},
    {"moon", "Moon", 301, 301, 4902.800118e9, 1737.4e3, 399},
    {"mars", "Mars", 499, 4, 42828.375816e9, 3389.5e3, 10},
    {"jupiter", "Jupiter", 599, 5, 126712764.1e9, 69911e3, 10},
    {"saturn", "Saturn", 699, 6, 37940584.8418e9, 58232e3, 10},
    {"uranus", "Uranus", 799, 7, 5794556.4e9, 25362e3, 10},
    {"neptune", "Neptune", 899, 8, 6836527.10058e9, 24622e3, 10},
    {"pluto", "Pluto", 999, 9, 975.5e9, 1188.3e3, 10},
}};

inline std::optional<Body> body(std::string_view key) {
  for (const Body& b : BODIES)
    if (b.key == key) return b;
  return std::nullopt;
}
inline std::optional<Body> body_by_naif(int id) {
  for (const Body& b : BODIES)
    if (b.naif == id || b.ephemeris == id) return b;
  return std::nullopt;
}

}  // namespace astronomy::ephemerides
