// The World Magnetic Model: the Earth's main magnetic field and its secular variation at a place
// and an instant, as NOAA NCEI and the British Geological Survey define it for 2025 to 2030.
//
// The main field of internal origin is −∇V, V a sum of spherical harmonics to degree N = 12 on
// the reference sphere a = 6 371 200 m, in geocentric spherical coordinates (λ, φ′, r), whose
// Gauss coefficients gₙᵐ, hₙᵐ (nT) vary linearly in time from the epoch t₀ by ġₙᵐ, ḣₙᵐ (nT per
// year). This header follows the WMM2025 technical report step by step (Chulliat et al. 2025,
// The US/UK World Magnetic Model for 2025–2030: Technical Report, NCEI, NOAA,
// doi:10.25923/prbc-s316, §1.2): a geodetic place on WGS 84 to (λ, φ′, r) (eqs. 7, 8); the
// coefficients at the instant (eq. 9); the geocentric components X′, Y′, Z′ and their rates
// (eqs. 10–15), with the Schmidt functions of schmidt.hpp (eq. 5); the rotation into the
// geodetic north, east and down (eqs. 17, 18); the elements H, F, I, D (eq. 19) and their
// secular variation (eq. 20); the grid variation in the polar regions (eq. 1). It adds the field
// on the Earth-fixed axes, the centred dipole of the first three coefficients (§1.6) with the
// McIlwain L of that dipole, the model's validity (§1.1, §1.8, §1.9) and its error model
// (§3.4). The theory, with the derivations, is domains/geophysics/theory/geomagnetism.md.
//
// Units. Within the infrastructure SI: tesla, radians, metres. The coefficient file and the
// report are in nanotesla, and the coefficients are kept so (their names say nT); every field
// returned is in tesla (×10⁻⁹ from nT). Rates are per year of the decimal-year scale (a calendar
// year, 365 or 366 days), as the model's coefficients are; the names say per_year.
//
// Frames. A place is geodetic on WGS 84 (latitude φ and longitude λ east of Greenwich in
// radians, height h along the normal in metres; geodesy::figure). X, Y, Z are along the
// geodetic north, east and down of the place: the ellipsoid's normal, not the radius, is the
// vertical (report §1.1.1). D is measured from true north through east, I from the horizontal,
// positive down. The Earth-fixed vector is on the axes of WGS 84, which at this model's
// accuracy are the ITRS's.
//
// Accuracy. The computation meets the report's high-precision example (Table 3b) to about
// 10⁻¹⁰ nT, and NOAA's 100 test values to their printed digits but for X and the elements made
// from it, where the file departs by up to 7.2 × 10⁻⁴ nT from the report's equations and, alike,
// from NOAA's own C library (the checks, and the theory note's §9). Such differences are far
// below the 0.1 nT the report allows an implementation (§1.10). The model itself is far
// coarser: its uncertainties are about 140 nT in X and Z and 0.3° in D away from the poles
// (§3.4); it leaves out the crustal field at
// wavelengths below about 3 200 km, which can reach thousands of nT locally, and all fields of
// external origin (ionosphere, magnetosphere) with the currents they induce, which magnetic
// storms raise to hundreds of nT; its secular variation is a prediction, linear in time.

#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <geodesy/figure/ellipsoid.hpp>
#include <geodesy/figure/local.hpp>
#include <geophysics/geomagnetism/schmidt.hpp>
#include <metrology/time/julian.hpp>
#include <physis/geometry/rotation.hpp>

namespace geophysics::geomagnetism {

using geodesy::figure::Geodetic;
using physis::geometry::Vec3;

// ---- units, constants and the model's validity ------------------------------------------------------

inline constexpr double TESLA_PER_NANOTESLA = 1e-9;
// The magnetic constant μ₀ (N A⁻²), the CODATA 2022 recommended value (as scipy.constants 1.18
// gives it); since the SI of 2019 it is measured, and differs from 4π × 10⁻⁷ by 5 × 10⁻¹⁰.
inline constexpr double MAGNETIC_CONSTANT = 1.25663706127e-6;

// The report's statements of validity. A model is valid from its epoch t₀ to t₀ + 5 years,
// WMM2025 "from 2025.0 to 2030.0" (§1.1), "to be used from January 1, 2025 to December 31, 2029"
// (§1.9), and NOAA's software accepts dates up to t₀ + 5 inclusive. It gives the field "from
// 1 km below the WGS 84 ellipsoid surface to at least 850 km above it" (§1.1, after the military
// specification MIL-PRF-89500B); §3.5 finds it within that specification higher in quiet times
// (Table 17), which is not assumed here.
inline constexpr double LIFETIME_YEARS = 5;
inline constexpr double LOWEST_HEIGHT_M = -1000;
inline constexpr double HIGHEST_HEIGHT_M = 850e3;
// Where the horizontal intensity H is weak the compass, and the model's D, fail: the blackout
// zone where H < 2000 nT, and the caution zone around it where 2000 nT ≤ H < 6000 nT (§1.8).
inline constexpr double BLACKOUT_HORIZONTAL_T = 2000e-9;
inline constexpr double CAUTION_HORIZONTAL_T = 6000e-9;
// The grid variation is defined poleward of 55° (§1.1.2, eq. 1; NOAA's software takes |φ| ≥ 55°).
inline constexpr double GRID_VARIATION_LATITUDE_RAD = 55 * std::numbers::pi / 180;

// ---- the coefficients ------------------------------------------------------------------------------

struct Coefficients {
  std::string model;         // the file's name of the model, "WMM-2025"
  std::string release_date;  // as the file writes it, month/day/year: "11/13/2024"
  double epoch_year = 0;     // t₀, a decimal year
  int degree = 0;            // N
  double reference_radius_m = 6371200.0;  // a, the geomagnetic reference radius (report §1.2)
  // gₙᵐ and hₙᵐ at t₀ (nT) and their secular variation (nT per year), at SchmidtFunctions::index(n, m);
  // the n = 0 entry is 0, as is h for m = 0.
  std::vector<double> g_nT, h_nT, g_nT_per_year, h_nT_per_year;

  static constexpr std::size_t index(int n, int m) { return SchmidtFunctions::index(n, m); }
  // The coefficients at a decimal year t: gₙᵐ(t) = gₙᵐ(t₀) + (t − t₀) ġₙᵐ(t₀) (report eq. 9), nT.
  double g(int n, int m, double year) const { return g_nT[index(n, m)] + (year - epoch_year) * g_nT_per_year[index(n, m)]; }
  double h(int n, int m, double year) const { return h_nT[index(n, m)] + (year - epoch_year) * h_nT_per_year[index(n, m)]; }
  double valid_from_year() const { return epoch_year; }
  double valid_until_year() const { return epoch_year + LIFETIME_YEARS; }

  // The model cut to a lower degree (degree 1 is its centred dipole).
  Coefficients truncated(int to_degree) const {
    if (to_degree < 1 || to_degree > degree) throw std::invalid_argument("geophysics::geomagnetism: no such degree to truncate to");
    Coefficients c = *this;
    c.degree = to_degree;
    const std::size_t size = index(to_degree, to_degree) + 1;
    for (auto* v : {&c.g_nT, &c.h_nT, &c.g_nT_per_year, &c.h_nT_per_year}) v->resize(size);
    return c;
  }
};

// The coefficients from the text of a coefficient file in NOAA's format (WMM.COF): the header
// line "epoch model release-date", rows "n m g h ġ ḣ", and a line of 9s that ends them.
inline Coefficients parse_cof(const std::string& text) {
  std::istringstream in(text);
  std::string line;
  Coefficients c;
  bool header = false, ended = false;
  struct Row { int n, m; double g, h, gd, hd; };
  std::vector<Row> rows;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.find_first_not_of(" \t") == std::string::npos) continue;
    std::istringstream f(line);
    if (!header) {
      if (!(f >> c.epoch_year >> c.model >> c.release_date)) throw std::invalid_argument("geophysics::geomagnetism: the coefficient file's header is not 'epoch model date'");
      header = true;
      continue;
    }
    if (line.find_first_not_of(" \t") == line.find("9999")) { ended = true; break; }
    Row r{};
    if (!(f >> r.n >> r.m >> r.g >> r.h >> r.gd >> r.hd)) throw std::invalid_argument("geophysics::geomagnetism: a coefficient row is not 'n m g h g-dot h-dot': " + line);
    if (r.n < 1 || r.m < 0 || r.m > r.n) throw std::invalid_argument("geophysics::geomagnetism: a coefficient row has no valid degree and order: " + line);
    rows.push_back(r);
  }
  if (!header) throw std::invalid_argument("geophysics::geomagnetism: the coefficient file is empty");
  if (!ended) throw std::invalid_argument("geophysics::geomagnetism: the coefficient file has no line of 9s at its end");
  for (const Row& r : rows) c.degree = std::max(c.degree, r.n);
  if (c.degree < 1) throw std::invalid_argument("geophysics::geomagnetism: the coefficient file holds no coefficients");
  const std::size_t size = Coefficients::index(c.degree, c.degree) + 1;
  c.g_nT.assign(size, 0.0);
  c.h_nT.assign(size, 0.0);
  c.g_nT_per_year.assign(size, 0.0);
  c.h_nT_per_year.assign(size, 0.0);
  std::vector<int> seen(size, 0);
  for (const Row& r : rows) {
    const std::size_t k = Coefficients::index(r.n, r.m);
    if (seen[k]++) throw std::invalid_argument("geophysics::geomagnetism: degree " + std::to_string(r.n) + " order " + std::to_string(r.m) + " appears twice");
    c.g_nT[k] = r.g;
    c.h_nT[k] = r.h;
    c.g_nT_per_year[k] = r.gd;
    c.h_nT_per_year[k] = r.hd;
  }
  for (int n = 1; n <= c.degree; n++)
    for (int m = 0; m <= n; m++)
      if (!seen[Coefficients::index(n, m)]) throw std::invalid_argument("geophysics::geomagnetism: degree " + std::to_string(n) + " order " + std::to_string(m) + " is missing");
  return c;
}

// The coefficients of a coefficient file (the data library's noaa-wmm/WMM.COF), whose path is passed in.
inline Coefficients read_cof(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::invalid_argument("geophysics::geomagnetism: cannot read the coefficient file " + path);
  std::stringstream b;
  b << in.rdbuf();
  return parse_cof(b.str());
}

// ---- time ---------------------------------------------------------------------------------------------

// The decimal year of a Gregorian date and a fraction of its day from 0h: the year plus the days
// elapsed since its 1 January, 0h, over the days of that year (365 or 366). This is the convention
// of NOAA's software (GeomagnetismLibrary.c, MAG_DateToYear, which takes 0h of the day), extended
// to a time of day. The date is UTC, the civil date; the scale does not matter at this model's
// accuracy: TT − UTC = 69.184 s is 2.2 × 10⁻⁶ year, over which the field changes by less than
// 5 × 10⁻⁴ nT (the secular variation stays below 230 nT per year at the surface; the checks).
// The convention does matter at the level of 0.1 nT: decimal years of Julian years (365.25 days
// from J2000.0) stray from these by up to half a day in 2025–2030, and the field by up to 0.3 nT.
inline double decimal_year(int year, int month, int day, double fraction_of_day = 0) {
  if (day < 1 || day > metrology::time::days_in_month(year, month)) throw std::invalid_argument("geophysics::geomagnetism: no such day in that month");
  int before = 0;
  for (int m = 1; m < month; m++) before += metrology::time::days_in_month(year, m);
  const int days = metrology::time::is_leap_year(year) ? 366 : 365;
  return year + (before + (day - 1) + fraction_of_day) / days;
}
// The decimal year of a Julian date in UTC (metrology::time's two-part Julian date).
inline double decimal_year(const metrology::time::JulianDate& utc) {
  const metrology::time::CalendarDate c = metrology::time::jd_to_calendar(utc);
  return decimal_year(c.year, c.month, c.day, c.fraction);
}

// ---- positions --------------------------------------------------------------------------------------------

struct GeocentricPosition {
  double longitude_rad = 0;  // λ, east of Greenwich: the geodetic longitude
  double latitude_rad = 0;   // φ′, geocentric
  double radius_m = 0;       // r, from the Earth's centre
};

inline const geodesy::figure::Ellipsoid& wgs84() {
  static const geodesy::figure::Ellipsoid e = geodesy::figure::Ellipsoid::wgs84();
  return e;
}

// A geodetic place on WGS 84 in geocentric spherical coordinates (report eqs. 7, 8):
// p = (R_c + h) cos φ, z = (R_c (1 − e²) + h) sin φ, r = √(p² + z²), φ′ = arcsin(z/r), with
// R_c = A/√(1 − e² sin² φ) and e² = f(2 − f); φ′ is taken as atan2(z, p), the same angle without
// the arcsine's loss of precision near the poles.
inline GeocentricPosition geocentric_position(const Geodetic& place) {
  const geodesy::figure::Ellipsoid& e = wgs84();
  const double rc = e.prime_vertical_radius(place.latitude);
  const double p = (rc + place.height) * std::cos(place.latitude), z = (rc * (1 - e.e2()) + place.height) * std::sin(place.latitude);
  return {place.longitude, std::atan2(z, p), std::hypot(p, z)};
}
// An Earth-centred position (m, ITRS axes) in geocentric spherical coordinates.
inline GeocentricPosition geocentric_position(const Vec3& itrs_m) {
  const double p = std::hypot(itrs_m[0], itrs_m[1]);
  return {p > 0 ? std::atan2(itrs_m[1], itrs_m[0]) : 0.0, std::atan2(itrs_m[2], p), std::hypot(p, itrs_m[2])};
}

// ---- the field in geocentric spherical components -----------------------------------------------------

struct GeocentricField {
  double north_T = 0, east_T = 0, down_T = 0;                                // X′, Y′, Z′
  double north_T_per_year = 0, east_T_per_year = 0, down_T_per_year = 0;     // Ẋ′, Ẏ′, Ż′
};

// X′, Y′, Z′ and their rates at a geocentric position and a decimal year (report eqs. 10–15):
//   X′ = −Σₙ (a/r)ⁿ⁺² Σₘ (gₙᵐ cos mλ + hₙᵐ sin mλ) dP̆ₙᵐ/dφ′,
//   Y′ = (1/cos φ′) Σₙ (a/r)ⁿ⁺² Σₘ m (gₙᵐ sin mλ − hₙᵐ cos mλ) P̆ₙᵐ,
//   Z′ = −Σₙ (n + 1)(a/r)ⁿ⁺² Σₘ (gₙᵐ cos mλ + hₙᵐ sin mλ) P̆ₙᵐ,
// and the same with ġ, ḣ for the rates. Y′ takes P̆ₙᵐ/cos φ′ directly, so that it holds at the
// poles, where the north and east of a point are those of the meridian λ given (report §1.4).
inline GeocentricField geocentric_field(const Coefficients& c, const GeocentricPosition& at, double year) {
  if (!(at.radius_m > 0)) throw std::invalid_argument("geophysics::geomagnetism: the field is not defined at the Earth's centre");
  const SchmidtFunctions s = schmidt_functions(c.degree, at.latitude_rad);
  const double dt = year - c.epoch_year, ratio = c.reference_radius_m / at.radius_m;
  double x = 0, y = 0, z = 0, xd = 0, yd = 0, zd = 0;
  double power = ratio * ratio;  // (a/r)ⁿ⁺², from n = 1
  for (int n = 1; n <= c.degree; n++) {
    power *= ratio;
    for (int m = 0; m <= n; m++) {
      const std::size_t k = Coefficients::index(n, m);
      const double cm = std::cos(m * at.longitude_rad), sm = std::sin(m * at.longitude_rad);
      const double g = c.g_nT[k] + dt * c.g_nT_per_year[k], h = c.h_nT[k] + dt * c.h_nT_per_year[k];
      const double gd = c.g_nT_per_year[k], hd = c.h_nT_per_year[k];
      const double cos_part = g * cm + h * sm, cos_rate = gd * cm + hd * sm;
      x -= power * cos_part * s.dp[k];
      xd -= power * cos_rate * s.dp[k];
      z -= power * (n + 1) * cos_part * s.p[k];
      zd -= power * (n + 1) * cos_rate * s.p[k];
      if (m > 0) {
        y += power * m * (g * sm - h * cm) * s.p_over_cos[k];
        yd += power * m * (gd * sm - hd * cm) * s.p_over_cos[k];
      }
    }
  }
  const double t = TESLA_PER_NANOTESLA;
  return {t * x, t * y, t * z, t * xd, t * yd, t * zd};
}

// The potential V (T m) at a geocentric position and a decimal year, whose negative gradient is
// the field (report eqs. 3, 4): V = a Σₙ (a/r)ⁿ⁺¹ Σₘ (gₙᵐ cos mλ + hₙᵐ sin mλ) P̆ₙᵐ(sin φ′).
inline double potential_T_m(const Coefficients& c, const GeocentricPosition& at, double year) {
  const SchmidtFunctions s = schmidt_functions(c.degree, at.latitude_rad);
  const double ratio = c.reference_radius_m / at.radius_m;
  double v = 0, power = ratio;
  for (int n = 1; n <= c.degree; n++) {
    power *= ratio;
    for (int m = 0; m <= n; m++)
      v += power * (c.g(n, m, year) * std::cos(m * at.longitude_rad) + c.h(n, m, year) * std::sin(m * at.longitude_rad)) * s.p[Coefficients::index(n, m)];
  }
  return c.reference_radius_m * v * TESLA_PER_NANOTESLA;
}

// ---- the elements at a place ----------------------------------------------------------------------------------

struct MagneticElements {
  double north_T = 0, east_T = 0, down_T = 0;  // X, Y, Z: along the geodetic north, east and down
  double horizontal_T = 0, total_T = 0;        // H, F
  double inclination_rad = 0;                  // I, from the horizontal, positive down; in [−π/2, π/2]
  double declination_rad = 0;                  // D, from true north through east; in (−π, π]
  // GV = D − λ poleward of 55° N, D + λ poleward of 55° S (report eq. 1), in (−π, π]; NaN elsewhere.
  double grid_variation_rad = std::numeric_limits<double>::quiet_NaN();
};

struct SecularVariation {
  double north_T_per_year = 0, east_T_per_year = 0, down_T_per_year = 0;  // Ẋ, Ẏ, Ż
  double horizontal_T_per_year = 0, total_T_per_year = 0;                 // Ḣ, Ḟ
  double inclination_rad_per_year = 0, declination_rad_per_year = 0;      // İ, Ḋ
  double grid_variation_rad_per_year = std::numeric_limits<double>::quiet_NaN();  // equal to Ḋ where GV is defined
};

struct Field {
  MagneticElements elements;
  SecularVariation secular_variation;
  GeocentricPosition position;  // (λ, φ′, r) of the place
  GeocentricField geocentric;   // X′, Y′, Z′ and their rates
  Vec3 itrs_T{};                // B on the Earth-fixed axes (ITRS, as WGS 84 realises them)
  Vec3 itrs_T_per_year{};       // dB/dt on the same axes
};

namespace detail {
inline double wrap_pi(double a) {
  double w = std::remainder(a, 2 * std::numbers::pi);
  if (w <= -std::numbers::pi) w += 2 * std::numbers::pi;
  return w;
}
}  // namespace detail

// The field at a geodetic place on WGS 84 and a decimal year.
inline Field field(const Coefficients& c, const Geodetic& place, double year) {
  Field f;
  f.position = geocentric_position(place);
  f.geocentric = geocentric_field(c, f.position, year);
  const GeocentricField& g = f.geocentric;
  // From the sphere's north and down to the ellipsoid's (report eqs. 17, 18), ψ = φ′ − φ.
  const double psi = f.position.latitude_rad - place.latitude, cp = std::cos(psi), sp = std::sin(psi);
  MagneticElements& e = f.elements;
  SecularVariation& v = f.secular_variation;
  e.north_T = g.north_T * cp - g.down_T * sp;
  e.east_T = g.east_T;
  e.down_T = g.north_T * sp + g.down_T * cp;
  v.north_T_per_year = g.north_T_per_year * cp - g.down_T_per_year * sp;
  v.east_T_per_year = g.east_T_per_year;
  v.down_T_per_year = g.north_T_per_year * sp + g.down_T_per_year * cp;
  // The elements (report eq. 19) and their rates (eq. 20); at H = 0, a dip pole, D and the
  // rates that divide by H are undefined.
  const double X = e.north_T, Y = e.east_T, Z = e.down_T;
  const double Xd = v.north_T_per_year, Yd = v.east_T_per_year, Zd = v.down_T_per_year;
  e.horizontal_T = std::hypot(X, Y);
  e.total_T = std::hypot(e.horizontal_T, Z);
  e.inclination_rad = std::atan2(Z, e.horizontal_T);
  e.declination_rad = std::atan2(Y, X);
  const double H = e.horizontal_T, F = e.total_T;
  v.horizontal_T_per_year = (X * Xd + Y * Yd) / H;
  v.total_T_per_year = (X * Xd + Y * Yd + Z * Zd) / F;
  v.inclination_rad_per_year = (H * Zd - Z * v.horizontal_T_per_year) / (F * F);
  v.declination_rad_per_year = (X * Yd - Y * Xd) / (H * H);
  if (place.latitude >= GRID_VARIATION_LATITUDE_RAD || place.latitude <= -GRID_VARIATION_LATITUDE_RAD) {
    e.grid_variation_rad = detail::wrap_pi(place.latitude > 0 ? e.declination_rad - place.longitude : e.declination_rad + place.longitude);
    v.grid_variation_rad_per_year = v.declination_rad_per_year;
  }
  // North, east and down of the place are the rows of geodesy's Earth-fixed → NED rotation.
  using namespace physis::geometry::operators;
  const physis::geometry::Mat3 to_ned = geodesy::figure::ecef_to_ned(place.latitude, place.longitude);
  const physis::geometry::Mat3 from_ned = physis::geometry::transpose(to_ned);
  f.itrs_T = from_ned * Vec3{X, Y, Z};
  f.itrs_T_per_year = from_ned * Vec3{Xd, Yd, Zd};
  return f;
}

// The field at an Earth-centred position (m, ITRS axes), its geodetic place found on WGS 84
// (geodesy::figure, valid outside about 45 km of the centre).
inline Field field(const Coefficients& c, const Vec3& itrs_position_m, double year) {
  return field(c, wgs84().geodetic(itrs_position_m), year);
}

// ---- validity ----------------------------------------------------------------------------------------------

enum class CompassZone { reliable, caution, blackout };

struct Validity {
  bool in_time = true;    // t₀ ≤ t ≤ t₀ + 5
  bool in_height = true;  // −1 km ≤ h ≤ 850 km
  CompassZone zone = CompassZone::reliable;
  std::string reasons;    // why not, in words; empty when valid and the compass reliable
  bool valid() const { return in_time && in_height; }
  bool declination_reliable() const { return valid() && zone == CompassZone::reliable; }
};

// Whether a place and an instant are within the model's validity, and whether its declination
// can be trusted there. The zones are judged by H at the place, as NOAA's calculators warn; the
// report draws them on the ellipsoid (§1.8), which is the same at h = 0.
inline Validity validity(const Coefficients& c, const Geodetic& place, double year, const Field& at_place) {
  Validity v;
  std::ostringstream why;
  auto fixed = [](double x, int digits) { std::ostringstream s; s.setf(std::ios::fixed); s.precision(digits); s << x; return s.str(); };
  auto add = [&](const std::string& s) { why << (why.tellp() > 0 ? "; " : "") << s; };
  if (!(year >= c.valid_from_year() && year <= c.valid_until_year())) {
    v.in_time = false;
    add("the time " + fixed(year, 4) + " is outside " + c.model + "'s span, " + fixed(c.valid_from_year(), 1) + " to " + fixed(c.valid_until_year(), 1));
  }
  if (!(place.height >= LOWEST_HEIGHT_M && place.height <= HIGHEST_HEIGHT_M)) {
    v.in_height = false;
    add("the height " + fixed(place.height / 1000, 3) + " km is outside the range the model is specified for, 1 km below to 850 km above the WGS 84 ellipsoid");
  }
  const double H = at_place.elements.horizontal_T;
  if (H < BLACKOUT_HORIZONTAL_T) {
    v.zone = CompassZone::blackout;
    add("H = " + fixed(H / TESLA_PER_NANOTESLA, 1) + " nT: the blackout zone about a magnetic pole (H < 2000 nT), where the declination is unreliable and a compass cannot be trusted");
  } else if (H < CAUTION_HORIZONTAL_T) {
    v.zone = CompassZone::caution;
    add("H = " + fixed(H / TESLA_PER_NANOTESLA, 1) + " nT: the caution zone about a magnetic pole (2000 nT <= H < 6000 nT), where the declination may be degraded");
  }
  v.reasons = why.str();
  return v;
}
inline Validity validity(const Coefficients& c, const Geodetic& place, double year) { return validity(c, place, year, field(c, place, year)); }

// ---- the error model ---------------------------------------------------------------------------------------

struct Uncertainty {
  double north_T = 0, east_T = 0, down_T = 0, horizontal_T = 0, total_T = 0;
  double inclination_rad = 0, declination_rad = 0;
};

// WMM2025's error model (report §3.4): constant uncertainties of X, Y, Z (eqs. 27–29: 137, 89 and
// 141 nT), of H, F and I (eqs. 37–39: 133 nT, 138 nT, 0.20°), and of D growing as 1/H about the
// magnetic poles, δD = √(0.26² + (5417/H)²) degrees with H in nT (eq. 43). They estimate the
// model's departure from the true field, crustal and disturbance fields included; they are
// WMM2025's, and do not describe another model's coefficients.
inline Uncertainty wmm2025_uncertainty(double horizontal_T) {
  const double nT = TESLA_PER_NANOTESLA, deg = std::numbers::pi / 180;
  const double h_nT = horizontal_T / nT;
  return {137 * nT, 89 * nT, 141 * nT, 133 * nT, 138 * nT, 0.20 * deg, std::hypot(0.26, 5417 / h_nT) * deg};
}

// ---- the centred dipole --------------------------------------------------------------------------------------

// The degree-1 part of the potential, V₁ = a³ (g · r̂)/r² with g = (g₁¹, h₁¹, g₁⁰) on the
// Earth-fixed axes, is the potential μ₀ (m · r̂)/(4π r²) of a dipole at the centre with moment
// m = 4π a³ g/μ₀. Its axis meets the sphere at the geomagnetic poles, the northern one in the
// direction −g (report §1.6). The pole's geodetic latitude is that of the point of the WGS 84
// ellipsoid on the axis, tan φ = tan φ′/(1 − e²). The south geomagnetic pole is antipodal.
struct Dipole {
  double year = 0;
  double g10_T = 0, g11_T = 0, h11_T = 0;      // g₁⁰, g₁¹, h₁¹ at the instant
  double equatorial_field_T = 0;               // B₀ = √(g₁⁰² + g₁¹² + h₁¹²): |B| on the dipole's equator at r = a
  double moment_A_m2 = 0;                      // |m| = 4π a³ B₀/μ₀
  Vec3 north_axis{};                           // unit vector from the centre to the geomagnetic north pole; m points along −north_axis
  double north_pole_latitude_rad = 0;          // geocentric
  double north_pole_geodetic_latitude_rad = 0;
  double north_pole_longitude_rad = 0;
  double tilt_rad = 0;                         // the angle between the dipole's axis and the rotation axis
  double reference_radius_m = 0;               // a

  Vec3 moment_vector_A_m2() const { return {-moment_A_m2 * north_axis[0], -moment_A_m2 * north_axis[1], -moment_A_m2 * north_axis[2]}; }
  // The dipole's field at an Earth-centred position (T, ITRS axes):
  // B = (a/r)³ (3 (g · r̂) r̂ − g) = μ₀/(4π r³) (3 (m · r̂) r̂ − m).
  Vec3 field_T(const Vec3& itrs_position_m) const {
    const double r = physis::geometry::norm(itrs_position_m);
    const Vec3 u{itrs_position_m[0] / r, itrs_position_m[1] / r, itrs_position_m[2] / r};
    const Vec3 g{-equatorial_field_T * north_axis[0], -equatorial_field_T * north_axis[1], -equatorial_field_T * north_axis[2]};
    const double k = std::pow(reference_radius_m / r, 3), gu = physis::geometry::dot(g, u);
    return {k * (3 * gu * u[0] - g[0]), k * (3 * gu * u[1] - g[1]), k * (3 * gu * u[2] - g[2])};
  }
};

inline Dipole dipole(const Coefficients& c, double year) {
  Dipole d;
  d.year = year;
  d.reference_radius_m = c.reference_radius_m;
  d.g10_T = c.g(1, 0, year) * TESLA_PER_NANOTESLA;
  d.g11_T = c.g(1, 1, year) * TESLA_PER_NANOTESLA;
  d.h11_T = c.h(1, 1, year) * TESLA_PER_NANOTESLA;
  d.equatorial_field_T = std::sqrt(d.g10_T * d.g10_T + d.g11_T * d.g11_T + d.h11_T * d.h11_T);
  if (!(d.equatorial_field_T > 0)) throw std::domain_error("geophysics::geomagnetism: the model has no dipole");
  const double a = c.reference_radius_m;
  d.moment_A_m2 = 4 * std::numbers::pi * a * a * a * d.equatorial_field_T / MAGNETIC_CONSTANT;
  d.north_axis = {-d.g11_T / d.equatorial_field_T, -d.h11_T / d.equatorial_field_T, -d.g10_T / d.equatorial_field_T};
  d.north_pole_latitude_rad = std::atan2(d.north_axis[2], std::hypot(d.north_axis[0], d.north_axis[1]));
  d.north_pole_longitude_rad = std::atan2(d.north_axis[1], d.north_axis[0]);
  d.tilt_rad = std::acos(std::min(1.0, std::max(-1.0, d.north_axis[2])));
  d.north_pole_geodetic_latitude_rad = std::atan2(std::sin(d.north_pole_latitude_rad), (1 - wgs84().e2()) * std::cos(d.north_pole_latitude_rad));
  return d;
}

// The geomagnetic latitude of an Earth-centred position: its latitude about the centred dipole's axis.
inline double geomagnetic_latitude_rad(const Dipole& d, const Vec3& itrs_position_m) {
  const double r = physis::geometry::norm(itrs_position_m);
  return std::asin(std::min(1.0, std::max(-1.0, physis::geometry::dot(itrs_position_m, d.north_axis) / r)));
}

// McIlwain's L of an Earth-centred position in the dipole approximation: the field line of the
// centred dipole through the point is r = L a cos² λₘ (λₘ the geomagnetic latitude), so
// L = r/(a cos² λₘ), in units of the reference radius a. The real field's L (McIlwain 1961), by
// the adiabatic invariant along its drift shell, departs from this by several per cent over the
// South Atlantic and at high L.
inline double mcilwain_l(const Dipole& d, const Vec3& itrs_position_m) {
  const double r = physis::geometry::norm(itrs_position_m), lm = geomagnetic_latitude_rad(d, itrs_position_m);
  return r / (d.reference_radius_m * std::cos(lm) * std::cos(lm));
}

}  // namespace geophysics::geomagnetism
