// The galactocentric frame: the centre of the Milky Way as origin, x from the Sun's projection on
// the galactic midplane toward the centre, z toward the north galactic pole, the Sun at a height
// z☉ above the plane: the frame of the Galaxy's dynamics, and the root of a world larger than the
// Solar System.
//
// The definition is astropy's Galactocentric (astropy.coordinates, its "v4.0" parameters, the
// default since astropy 4.0), so that positions and velocities agree with the community's tool:
//
//   x_gc = A x_icrs + o,   v_gc = A v_icrs + v☉,
//   A = R_y(−asin(z☉/d)) R_x(η − roll) R_y(−δ_GC) R_z(α_GC),   o = −R_y(−asin(z☉/d)) (d, 0, 0)ᵀ,
//
// with R_x, R_y, R_z the frame (passive) rotations of SOFA's iauRx, iauRy, iauRz, x_icrs barycentric
// (BCRS) and:
//   (α_GC, δ_GC) = (266.4051°, −28.936175°)  Sgr A* in the ICRS (Reid and Brunthaler 2004, ApJ 616, 872);
//   d = 8.122 kpc                            the Sun's distance from it (GRAVITY Collaboration 2018, A&A 615, L15);
//   z☉ = 20.8 pc                             the Sun's height above the plane (Bennett and Bovy 2019, MNRAS 482, 1417);
//   v☉ = (12.9, 245.6, 7.78) km/s            the Sun's velocity about the centre (Drimmel and Poggio 2018,
//                                            RNAAS 2, 210; from Reid and Brunthaler 2004 and GRAVITY 2018);
//   η = 58.5986320306°                       astropy's roll, measured so that the plane of galactic longitude
//                                            0 and the frame's x–z plane agree; roll = 0.
//
// Remark. Aligning the north galactic pole of SOFA's galactic frame (galactic axes of the
// Hipparcos realisation) with the x–z plane instead of using η turns the frame by 1.3 × 10⁻⁶ rad
// (0.27″): two conventions, one tenth of a light-year apart at the Galaxy's centre. These are the
// parameters of a scenario, not constants of nature (Cosmopoiesis formalism §9): `Galactocentric`
// holds them, and a world may choose others.
//
// Caveat. The translation is that of the epoch at which a position is given, as astropy's is: the
// Sun's motion v☉ enters the velocities only, so that the origin does not drift with the
// ephemeris's time (15 au over its three centuries, 3 × 10⁻⁹ of the distance to the centre).

#pragma once

#include <cmath>
#include <numbers>

#include <physis/geometry/frame_graph.hpp>

namespace navigation::frames {

inline constexpr double AU_M = 149597870700.0;                               // IAU 2012 Resolution B2, exact
inline constexpr double PARSEC_M = AU_M * 648000.0 / std::numbers::pi;       // IAU 2015 Resolution B2
inline constexpr double LIGHT_YEAR_M = 299792458.0 * 365.25 * 86400.0;      // c times the Julian year (IAU)

struct Galactocentric {
  double ra = 266.4051 * std::numbers::pi / 180;
  double dec = -28.936175 * std::numbers::pi / 180;
  double distance = 8.122e3 * PARSEC_M;
  double z_sun = 20.8 * PARSEC_M;
  physis::geometry::Vec3 v_sun{12.9e3, 245.6e3, 7.78e3};  // m/s, in the galactocentric axes
  double roll = 0;                                        // rad
  static constexpr double roll0 = 58.5986320306 * std::numbers::pi / 180;

  // The frame (passive) rotations of SOFA (iauRx, iauRy, iauRz).
  static physis::geometry::Mat3 rx(double a) { const double c = std::cos(a), s = std::sin(a); return {{{1, 0, 0}, {0, c, s}, {0, -s, c}}}; }
  static physis::geometry::Mat3 ry(double a) { const double c = std::cos(a), s = std::sin(a); return {{{c, 0, -s}, {0, 1, 0}, {s, 0, c}}}; }
  static physis::geometry::Mat3 rz(double a) { const double c = std::cos(a), s = std::sin(a); return {{{c, s, 0}, {-s, c, 0}, {0, 0, 1}}}; }

  physis::geometry::Mat3 matrix() const {
    using namespace physis::geometry::operators;
    return ry(-std::asin(z_sun / distance)) * rx(roll0 - roll) * ry(-dec) * rz(ra);
  }
  physis::geometry::Vec3 offset() const {
    using namespace physis::geometry::operators;
    return -1.0 * (ry(-std::asin(z_sun / distance)) * physis::geometry::Vec3{distance, 0, 0});
  }
  // The edge of the frame graph, galactocentric → ICRS: x_icrs = Aᵀ (x_gc − o), v_icrs = Aᵀ (v_gc − v☉).
  physis::geometry::StateTransform to_icrs() const {
    using namespace physis::geometry::operators;
    const auto at = physis::geometry::transpose(matrix());
    return {at, physis::geometry::zero3(), -1.0 * (at * offset()), -1.0 * (at * v_sun)};
  }
};

}  // namespace navigation::frames
