// The Universal Transverse Mercator (UTM) and Universal Polar Stereographic (UPS) grids on WGS 84
// (NGA.SIG.0012_2.0.0_UTMUPS, 2014: "The Universal Grids and the Transverse Mercator and Polar
// Stereographic Map Projections").
//
// UTM divides the Earth from 80° S to 84° N into sixty zones of 6° of longitude, zone z centred on
// λ₀ = 6z − 183°, each a transverse Mercator with k₀ = 0.9996, a false easting of 500 000 m and a
// false northing of 0 in the north and 10 000 000 m in the south. Two regions depart from the rule
// (NGA.SIG.0012, §2.3): zone 32 is widened to 3°–12° E between 56° and 64° N (south-west Norway),
// and between 72° and 84° N only the odd zones 31, 33, 35 and 37 are used, over 0°–9°, 9°–21°,
// 21°–33° and 33°–42° E (Svalbard). The latitude bands C to X (I and O omitted) are 8° tall from
// 80° S, X 12° tall. UPS covers the polar caps beyond (north of 84° N, south of 80° S) by a polar
// stereographic projection with k₀ = 0.994 at the pole and false easting and northing 2 000 000 m.
//
// The polar stereographic projection of the ellipsoid (Snyder 1987, "Map Projections — A Working
// Manual", USGS PP 1395, eqs. 15-9, 21-33, 21-34; Karney's form with the conformal latitude):
//   ρ = 2 k₀ a t / √((1 + e)^(1+e) (1 − e)^(1−e)),   t = tan(π/4 − φ/2) ((1 + e sin φ)/(1 − e sin φ))^(e/2),
// in the north x = ρ sin λ, y = −ρ cos λ; in the south x = ρ sin λ, y = ρ cos λ (with t of −φ). The
// point scale is k = ρ / (a m), m = cos φ / √(1 − e² sin² φ); the convergence is λ in the north and
// −λ in the south.

#pragma once

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>

#include <geodesy/figure/ellipsoid.hpp>
#include <geodesy/projection/transverse_mercator.hpp>

namespace geodesy::projection {

inline constexpr double UTM_K0 = 0.9996;
inline constexpr double UTM_FALSE_EASTING = 500000.0;
inline constexpr double UTM_FALSE_NORTHING_SOUTH = 10000000.0;
inline constexpr double UPS_K0 = 0.994;
inline constexpr double UPS_FALSE = 2000000.0;
inline constexpr double DEG = std::numbers::pi / 180;

// A position on the universal grids: a UTM zone (1–60) or UPS (zone 0), its hemisphere, easting and
// northing (with their false origins), and the convergence and scale there.
struct UtmUps {
  int zone = 0;            // 1..60 for UTM, 0 for UPS
  bool north = true;
  double easting = 0, northing = 0;  // m
  double convergence = 0;  // rad: the bearing of grid north, clockwise from true north
  double scale = 1;
  char band = 0;           // the latitude band (UTM: C..X; UPS: A, B in the south, Y, Z in the north)
};

// The latitude band of a latitude (rad) in UTM's range: C (80°–72° S) to X (72°–84° N).
inline char utm_band(double lat) {
  static const char* bands = "CDEFGHJKLMNPQRSTUVWX";
  const double d = lat / DEG;
  if (d < -80 || d > 84) return 0;
  int i = int(std::floor((d + 80) / 8));
  if (i > 19) i = 19;  // X runs to 84°
  return bands[i];
}

// The standard zone of a point (rad): UTM 1–60 with the Norway and Svalbard exceptions, or 0 for UPS.
inline int standard_zone(double lat, double lon) {
  const double d = lat / DEG;
  double l = std::remainder(lon / DEG, 360.0);
  if (l >= 180) l -= 360;
  if (d < -80 || d >= 84) return 0;
  int zone = int(std::floor((l + 180) / 6)) + 1;
  if (zone > 60) zone = 60;
  if (d >= 56 && d < 64 && l >= 3 && l < 12) zone = 32;  // Norway
  if (d >= 72 && d < 84 && l >= 0 && l < 42) zone = l < 9 ? 31 : l < 21 ? 33 : l < 33 ? 35 : 37;  // Svalbard
  return zone;
}

inline double central_meridian(int zone) { return (6.0 * zone - 183.0) * DEG; }

class UniversalGrids {
 public:
  explicit UniversalGrids(const Ellipsoid& e = Ellipsoid::wgs84()) : e_(e), tm_(e, UTM_K0) {
    const double ecc = std::sqrt(e.e2());
    c_ = std::sqrt(std::pow(1 + ecc, 1 + ecc) * std::pow(1 - ecc, 1 - ecc));
    ecc_ = ecc;
  }

  // A point (rad) on the grids: in the zone given (1..60, 0 for UPS) or, by default (−1), its standard zone.
  UtmUps forward(double lat, double lon, int zone = -1) const {
    if (zone < 0) zone = standard_zone(lat, lon);
    UtmUps r;
    r.zone = zone;
    r.north = lat >= 0;
    if (zone == 0) {
      const PolarPoint p = polar_forward(lat, lon, r.north);
      r.easting = UPS_FALSE + p.x;
      r.northing = UPS_FALSE + p.y;
      r.convergence = p.convergence;
      r.scale = p.scale;
      double l = std::remainder(lon / DEG, 360.0);
      r.band = r.north ? (l < 0 ? 'Y' : 'Z') : (l < 0 ? 'A' : 'B');
      return r;
    }
    if (zone < 1 || zone > 60) throw std::invalid_argument("geodesy: a UTM zone is 1 to 60");
    const double dl = std::remainder(lon - central_meridian(zone), 2 * std::numbers::pi);
    const GridPoint g = tm_.forward(lat, dl);
    r.easting = UTM_FALSE_EASTING + g.easting;
    r.northing = (r.north ? 0 : UTM_FALSE_NORTHING_SOUTH) + g.northing;
    r.convergence = g.convergence;
    r.scale = g.scale;
    r.band = utm_band(lat);
    return r;
  }

  // The latitude and longitude (rad) of a point of the grids.
  geodesy::figure::Geodetic reverse(int zone, bool north, double easting, double northing, double* convergence = nullptr, double* scale = nullptr) const {
    geodesy::figure::Geodetic g;
    if (zone == 0) {
      const GeoPoint p = polar_reverse(easting - UPS_FALSE, northing - UPS_FALSE, north);
      g.latitude = p.latitude;
      g.longitude = p.longitude;
      if (convergence) *convergence = p.convergence;
      if (scale) *scale = p.scale;
      return g;
    }
    if (zone < 1 || zone > 60) throw std::invalid_argument("geodesy: a UTM zone is 1 to 60");
    const GeoPoint p = tm_.reverse(easting - UTM_FALSE_EASTING, northing - (north ? 0 : UTM_FALSE_NORTHING_SOUTH));
    g.latitude = p.latitude;
    g.longitude = std::remainder(p.longitude + central_meridian(zone), 2 * std::numbers::pi);
    if (convergence) *convergence = p.convergence;
    if (scale) *scale = p.scale;
    return g;
  }

  const TransverseMercator& transverse_mercator() const { return tm_; }
  const Ellipsoid& ellipsoid() const { return e_; }

 private:
  struct PolarPoint { double x, y, convergence, scale; };

  // t of Snyder eq. 15-9 for the latitude measured toward the projection's pole.
  double t_of(double phi) const {
    const double s = std::sin(phi);
    return std::tan(std::numbers::pi / 4 - phi / 2) * std::pow((1 + ecc_ * s) / (1 - ecc_ * s), ecc_ / 2);
  }

  PolarPoint polar_forward(double lat, double lon, bool north) const {
    const double phi = north ? lat : -lat;
    const double t = t_of(phi);
    const double rho = 2 * UPS_K0 * e_.a() * t / c_;
    PolarPoint p;
    p.x = rho * std::sin(lon);
    p.y = north ? -rho * std::cos(lon) : rho * std::cos(lon);
    const double s = std::sin(phi), m = std::cos(phi) / std::sqrt(1 - e_.e2() * s * s);
    p.scale = m > 0 ? rho / (e_.a() * m) : UPS_K0;
    p.convergence = north ? lon : -lon;
    return p;
  }

  GeoPoint polar_reverse(double x, double y, bool north) const {
    const double rho = std::hypot(x, y);
    const double t = rho * c_ / (2 * UPS_K0 * e_.a());
    // φ from t by iteration (Snyder eq. 7-9): χ = π/2 − 2 atan t, then φ = π/2 − 2 atan(t ((1 − e sin φ)/(1 + e sin φ))^(e/2)).
    double phi = std::numbers::pi / 2 - 2 * std::atan(t);
    for (int i = 0; i < 20; i++) {
      const double s = std::sin(phi);
      const double next = std::numbers::pi / 2 - 2 * std::atan(t * std::pow((1 - ecc_ * s) / (1 + ecc_ * s), ecc_ / 2));
      if (std::fabs(next - phi) < 1e-15) { phi = next; break; }
      phi = next;
    }
    GeoPoint g;
    g.latitude = north ? phi : -phi;
    g.longitude = rho == 0 ? 0.0 : (north ? std::atan2(x, -y) : std::atan2(x, y));
    const double s = std::sin(phi), m = std::cos(phi) / std::sqrt(1 - e_.e2() * s * s);
    g.scale = m > 0 ? rho / (e_.a() * m) : UPS_K0;
    g.convergence = north ? g.longitude : -g.longitude;
    return g;
  }

  Ellipsoid e_;
  TransverseMercator tm_;
  double c_ = 1, ecc_ = 0;
};

}  // namespace geodesy::projection
