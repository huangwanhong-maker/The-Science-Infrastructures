// Checks of the geomagnetism module (geophysics/geomagnetism): the World Magnetic Model 2025.
//
// References, all NOAA's (tests/fixtures/, written by tests/make_wmm_fixtures.py from the files
// the data library fetches, python data/acquisition/fetch.py noaa-wmm):
// - the technical report's high-precision example (Chulliat et al. 2025, Tables 3a, 3b, 3d), its
//   field above the geographic poles (§1.4), its geomagnetic poles (§1.6) and the construction of
//   its error model (§3.4): wmm2025-report.txt;
// - the 100 test points distributed with the coefficients: wmm2025-test-values.txt;
// - the 12 test points of the report's Table 6, with the grid variation: wmm2025-table6.txt.
// And the mathematics: the field as the gradient of the potential; the rates as derivatives in time; the Earth-fixed vector against the geocentric basis; the degree-1 model
// as the field of the centred dipole; McIlwain's L constant along the dipole's field lines.
// The Schmidt functions are calx's and are checked there (special_functions/tests/check_legendre.cpp).
// The coefficients come from the data library's cache (noaa-wmm/WMM.COF); without them only the
// checks of the mathematics run.

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <geophysics/geomagnetism/wmm.hpp>

using namespace geophysics::geomagnetism;
using calx::verification::Suite;
using physis::geometry::Mat3;

namespace {

constexpr double pi = std::numbers::pi, deg = pi / 180, nT = 1e-9;

std::string cache_path(const std::string& relative) {
  const char* env = std::getenv("SCIENCE_CACHE");
  return (env ? std::string(env) : std::string(GEOPHYSICS_ROOT "/../../data/cache")) + "/" + relative;
}

// Rows of numbers of a fixture, comment lines left out; "NaN" read as NaN.
std::vector<std::vector<double>> numeric_rows(const std::string& path) {
  std::vector<std::vector<double>> out;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream f(line);
    std::vector<double> row;
    std::string x;
    while (f >> x) row.push_back(x == "NaN" ? NAN : std::strtod(x.c_str(), nullptr));
    if (!row.empty()) out.push_back(row);
  }
  return out;
}

// The report's values, by key.
std::map<std::string, double> report_values(const std::string& path) {
  std::map<std::string, double> out;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream f(line);
    std::string key, value;
    if (f >> key >> value) out[key] = std::strtod(value.c_str(), nullptr);
  }
  return out;
}

double angle_difference(double a, double b) { return std::remainder(a - b, 2 * pi); }

}  // namespace

int main() {
  Suite s;

  // ---- decimal years -------------------------------------------------------------------------------------------------
  {
    s.require("decimal year: 1 January 0h is the year itself", decimal_year(2025, 1, 1) == 2025.0 && decimal_year(2028, 1, 1) == 2028.0);
    s.require("decimal year: NOAA's convention, days elapsed over the days of the year",
              decimal_year(2025, 7, 2) == 2025 + 182.0 / 365 && decimal_year(2028, 7, 2) == 2028 + 183.0 / 366 && decimal_year(2029, 12, 31, 1.0) == 2030.0);
    double worst = 0;
    calx::random::Stream r(20261028, 1, 0, 0);
    for (int i = 0; i < 200; i++) {
      const int y = 2024 + int(r.below(8)), m = 1 + int(r.below(12)), d = 1 + int(r.below(std::uint64_t(metrology::time::days_in_month(y, m))));
      const double f = r.uniform();
      const metrology::time::JulianDate jd = metrology::time::calendar_to_jd(y, m, d).plus_days(f);
      worst = std::max(worst, std::fabs(decimal_year(jd) - decimal_year(y, m, d, f)));
    }
    s.check("decimal year: from a Julian date as from its calendar date (200 instants)", worst, 1e-12, "yr");
    // Decimal years of Julian years, 2000 + (JD − 2451545)/365.25, stray from NOAA's by up to half a
    // day over the model's span (the theory note's §5.3), every 6 hours from 2025 to 2030.
    double stray = 0;
    const metrology::time::JulianDate start = metrology::time::calendar_to_jd(2025, 1, 1);
    for (int q = 0; q < 4 * 1826; q++) {
      const metrology::time::JulianDate t = start.plus_days(q * 0.25);
      const double julian = 2000 + (t.value() - metrology::time::J2000) / 365.25;
      stray = std::max(stray, std::fabs(julian - decimal_year(t)) * 365.25);
    }
    s.check("decimal year: Julian years from J2000.0 stray from NOAA's convention by up to half a day, 2025-2030", stray, 0.5000001, "day");
  }

  // ---- the coefficients ----------------------------------------------------------------------------------------------
  const std::string cof = cache_path("noaa-wmm/WMM.COF");
  if (!std::filesystem::exists(cof)) {
    s.skip("WMM2025's coefficients are not in the data library's cache (python data/acquisition/fetch.py noaa-wmm)");
    return s.done();
  }
  const Coefficients c = read_cof(cof);
  s.require("the coefficient file: WMM-2025, epoch 2025.0, released 11/13/2024, to degree 12",
            c.model == "WMM-2025" && c.epoch_year == 2025.0 && c.release_date == "11/13/2024" && c.degree == 12 && c.g_nT.size() == 91);
  {
    bool threw_missing = false, threw_end = false;
    std::ifstream in(cof);
    std::stringstream b;
    b << in.rdbuf();
    const std::string text = b.str();
    const std::size_t row = text.find("\n  3  2 ");
    try { parse_cof(text.substr(0, row) + text.substr(text.find('\n', row + 1))); } catch (const std::invalid_argument&) { threw_missing = true; }
    try { parse_cof(text.substr(0, text.find("9999"))); } catch (const std::invalid_argument&) { threw_end = true; }
    s.require("the coefficient file: a missing row, or a file without its line of 9s, is refused", threw_missing && threw_end);
  }

  const auto R = report_values(GEOPHYSICS_FIXTURES "/wmm2025-report.txt");
  s.require("the report's values were read", R.size() >= 70);
  auto rv = [&](const std::string& key) {
    const auto it = R.find(key);
    return it == R.end() ? NAN : it->second;
  };

  // ---- the report's high-precision example (Tables 3a, 3b) ----------------------------------------------------------
  {
    const double year = rv("table3a.time");
    const Geodetic place{rv("table3a.latitude") * deg, rv("table3a.longitude") * deg, rv("table3a.height") * 1000};
    const Field f = field(c, place, year);
    double coeffs = 0;
    for (int n = 1; n <= 2; n++)
      for (int m = 0; m <= n; m++) {
        const std::string nm = std::to_string(n) + "," + std::to_string(m);
        coeffs = std::max({coeffs, std::fabs(c.g(n, m, year) - rv("table3b.g(" + nm + ",t)")), std::fabs(c.h(n, m, year) - rv("table3b.h(" + nm + ",t)"))});
      }
    s.check("Table 3b: the coefficients at 2027.5 (eq. 9)", coeffs, 1e-9, "nT");
    s.check("Table 3b: the geocentric latitude φ′ (eqs. 7, 8)", f.position.latitude_rad - rv("table3b.phi-prime"), 1e-10, "rad");
    s.check("Table 3b: the radius r", f.position.radius_m - rv("table3b.r"), 1e-6, "m");
    const GeocentricField& g = f.geocentric;
    s.check("Table 3b: X′, Y′, Z′ (eqs. 10-12)",
            std::max({std::fabs(g.north_T / nT - rv("table3b.Xprime")), std::fabs(g.east_T / nT - rv("table3b.Yprime")), std::fabs(g.down_T / nT - rv("table3b.Zprime"))}), 1e-8, "nT");
    s.check("Table 3b: their rates (eqs. 13-15)",
            std::max({std::fabs(g.north_T_per_year / nT - rv("table3b.Xprime-dot")), std::fabs(g.east_T_per_year / nT - rv("table3b.Yprime-dot")),
                      std::fabs(g.down_T_per_year / nT - rv("table3b.Zprime-dot"))}), 1e-8, "nT/yr");
    const MagneticElements& e = f.elements;
    const SecularVariation& v = f.secular_variation;
    s.check("Table 3b: X, Y, Z, H, F (eqs. 17, 19)",
            std::max({std::fabs(e.north_T / nT - rv("table3b.X")), std::fabs(e.east_T / nT - rv("table3b.Y")), std::fabs(e.down_T / nT - rv("table3b.Z")),
                      std::fabs(e.horizontal_T / nT - rv("table3b.H")), std::fabs(e.total_T / nT - rv("table3b.F"))}), 1e-8, "nT");
    s.check("Table 3b: D and I (eq. 19)", std::max(std::fabs(e.declination_rad - rv("table3b.D")), std::fabs(e.inclination_rad - rv("table3b.I"))), 1e-10, "rad");
    s.check("Table 3b: the rates of X, Y, Z, H, F (eqs. 18, 20)",
            std::max({std::fabs(v.north_T_per_year / nT - rv("table3b.Xdot")), std::fabs(v.east_T_per_year / nT - rv("table3b.Ydot")),
                      std::fabs(v.down_T_per_year / nT - rv("table3b.Zdot")), std::fabs(v.horizontal_T_per_year / nT - rv("table3b.Hdot")),
                      std::fabs(v.total_T_per_year / nT - rv("table3b.Fdot"))}), 1e-8, "nT/yr");
    s.check("Table 3b: the rates of D and I (eq. 20)",
            std::max(std::fabs(v.declination_rad_per_year - rv("table3b.Ddot")), std::fabs(v.inclination_rad_per_year - rv("table3b.Idot"))), 1e-10, "rad/yr");
    // Table 3d: D (TrueN-to-MagN) and the grid variation of the south polar grid, D + λ (eq. 1).
    s.check("Table 3d: D in degrees", e.declination_rad / deg - rv("table3d.TrueN-to-MagN.UPS-south"), 1e-8, "°");
    s.check("Table 3d: the grid variation south of 55° S, D + λ, to a whole turn", angle_difference(e.grid_variation_rad, rv("table3d.GridN-to-MagN.UPS-south") * deg) / deg, 1e-8, "°");
  }

  // ---- NOAA's 100 test points ------------------------------------------------------------------------------------------
  const auto points = numeric_rows(GEOPHYSICS_FIXTURES "/wmm2025-test-values.txt");
  s.require("the 100 test points were read", points.size() == 100 && points[0].size() == 18);
  {
    // Fields: year, height (km), latitude, longitude (°), D, I (°), H, X, Y, Z, F (nT),
    // dD/dt, dI/dt (°/yr), dH/dt, dX/dt, dY/dt, dZ/dt, dF/dt (nT/yr).
    double angles = 0, y = 0, rates_yz = 0, rates_di = 0, xzhf = 0, rates_xhf = 0, round_trip = 0, basis = 0;
    int zones_agree = 0, blackout = 0, caution = 0;
    for (const auto& p : points) {
      const Geodetic place{p[2] * deg, p[3] * deg, p[1] * 1000};
      const Field f = field(c, place, p[0]);
      const MagneticElements& e = f.elements;
      const SecularVariation& v = f.secular_variation;
      angles = std::max({angles, std::fabs(e.declination_rad / deg - p[4]), std::fabs(e.inclination_rad / deg - p[5])});
      y = std::max(y, std::fabs(e.east_T / nT - p[8]));
      rates_yz = std::max({rates_yz, std::fabs(v.east_T_per_year / nT - p[15]), std::fabs(v.down_T_per_year / nT - p[16])});
      rates_di = std::max({rates_di, std::fabs(v.declination_rad_per_year / deg - p[11]), std::fabs(v.inclination_rad_per_year / deg - p[12])});
      xzhf = std::max({xzhf, std::fabs(e.north_T / nT - p[7]), std::fabs(e.down_T / nT - p[9]), std::fabs(e.horizontal_T / nT - p[6]), std::fabs(e.total_T / nT - p[10])});
      rates_xhf = std::max({rates_xhf, std::fabs(v.north_T_per_year / nT - p[14]), std::fabs(v.horizontal_T_per_year / nT - p[13]), std::fabs(v.total_T_per_year / nT - p[17])});
      // The Earth-fixed vector: back to north, east and down; and built from the sphere's basis.
      const physis::geometry::Vec3 ned = [&] {
        using namespace physis::geometry::operators;
        return geodesy::figure::ecef_to_ned(place.latitude, place.longitude) * f.itrs_T;
      }();
      round_trip = std::max({round_trip, std::fabs(ned[0] - e.north_T) / nT, std::fabs(ned[1] - e.east_T) / nT, std::fabs(ned[2] - e.down_T) / nT,
                             std::fabs(physis::geometry::norm(f.itrs_T) - e.total_T) / nT});
      const double lam = f.position.longitude_rad, phc = f.position.latitude_rad;
      const physis::geometry::Vec3 north{-std::sin(phc) * std::cos(lam), -std::sin(phc) * std::sin(lam), std::cos(phc)}, east{-std::sin(lam), std::cos(lam), 0},
          down{-std::cos(phc) * std::cos(lam), -std::cos(phc) * std::sin(lam), -std::sin(phc)};
      for (int i = 0; i < 3; i++)
        basis = std::max(basis, std::fabs(f.geocentric.north_T * north[i] + f.geocentric.east_T * east[i] + f.geocentric.down_T * down[i] - f.itrs_T[i]) / nT);
      // The compass zones, by H.
      const Validity val = validity(c, place, p[0], f);
      const CompassZone expected = p[6] < 2000 ? CompassZone::blackout : p[6] < 6000 ? CompassZone::caution : CompassZone::reliable;
      zones_agree += val.zone == expected && val.valid();
      blackout += expected == CompassZone::blackout;
      caution += expected == CompassZone::caution;
    }
    s.check("100 test points: D and I, printed to 0.01°", angles, 0.0050001, "°");
    s.check("100 test points: Y, to a unit of the file's last place", y, 1e-6, "nT");
    s.check("100 test points: the rates of Y and Z, to a unit of the file's last place", rates_yz, 1e-6, "nT/yr");
    s.check("100 test points: the rates of D and I, to a unit of the file's last place", rates_di, 1e-6, "°/yr");
    // NOAA's file departs from the report's equations in X by up to 7.2e-4 nT (Table 3b is met to
    // 1e-10 nT above), and so in H, F and, through the rotation by φ′ − φ, in Z. NOAA's own C
    // library (GeomagnetismLibrary.c of wmm2025_Linux) departs from the file alike, and agrees
    // with this module to the file's last place (theory note §9); the report allows
    // implementations 0.1 nT (§1.10).
    s.check("100 test points: X, Z, H, F (the file departs from the equations by up to 7.2e-4 nT in X)", xzhf, 1e-3, "nT");
    s.check("100 test points: the rates of X, H, F", rates_xhf, 2e-6, "nT/yr");
    s.check("100 test points: the Earth-fixed vector back to north, east and down, and its length F", round_trip, 1e-9, "nT");
    s.check("100 test points: the Earth-fixed vector is X′, Y′, Z′ on the sphere's basis (eq. 17 inverted)", basis, 1e-8, "nT");
    s.require("100 test points: the compass zones by H, " + std::to_string(blackout) + " in the blackout zone and " + std::to_string(caution) + " in the caution zone",
              zones_agree == 100 && blackout > 0 && caution > 0, std::to_string(zones_agree) + " of 100 agree");
  }

  // ---- the report's Table 6, with the grid variation ------------------------------------------------------------------
  {
    const auto rows = numeric_rows(GEOPHYSICS_FIXTURES "/wmm2025-table6.txt");
    s.require("Table 6: 12 points read", rows.size() == 12 && rows[0].size() == 19);
    // Fields: year, height (km), latitude, longitude, X, Y, Z, H, F (nT), I, D, GV (°),
    // Xdot, Ydot, Zdot, Hdot, Fdot (nT/yr), Idot, Ddot (°/yr).
    double field_nT = 0, angles = 0, rates = 0, angle_rates = 0;
    bool gv_defined_alike = true;
    for (const auto& p : rows) {
      const Field f = field(c, Geodetic{p[2] * deg, p[3] * deg, p[1] * 1000}, p[0]);
      const MagneticElements& e = f.elements;
      const SecularVariation& v = f.secular_variation;
      field_nT = std::max({field_nT, std::fabs(e.north_T / nT - p[4]), std::fabs(e.east_T / nT - p[5]), std::fabs(e.down_T / nT - p[6]),
                           std::fabs(e.horizontal_T / nT - p[7]), std::fabs(e.total_T / nT - p[8])});
      angles = std::max({angles, std::fabs(e.inclination_rad / deg - p[9]), std::fabs(e.declination_rad / deg - p[10])});
      if (std::isnan(p[11]) != std::isnan(e.grid_variation_rad)) gv_defined_alike = false;
      if (!std::isnan(p[11])) angles = std::max(angles, std::fabs(angle_difference(e.grid_variation_rad, p[11] * deg) / deg));
      rates = std::max({rates, std::fabs(v.north_T_per_year / nT - p[12]), std::fabs(v.east_T_per_year / nT - p[13]), std::fabs(v.down_T_per_year / nT - p[14]),
                        std::fabs(v.horizontal_T_per_year / nT - p[15]), std::fabs(v.total_T_per_year / nT - p[16])});
      angle_rates = std::max({angle_rates, std::fabs(v.inclination_rad_per_year / deg - p[17]), std::fabs(v.declination_rad_per_year / deg - p[18])});
    }
    s.check("Table 6: X, Y, Z, H, F, printed to 0.1 nT", field_nT, 0.0500001, "nT");
    s.check("Table 6: I, D and the grid variation, printed to 0.01°", angles, 0.0050001, "°");
    s.require("Table 6: the grid variation is defined poleward of 55° and nowhere else", gv_defined_alike);
    s.check("Table 6: the rates of X, Y, Z, H, F, printed to 0.1 nT/yr", rates, 0.0500001, "nT/yr");
    s.check("Table 6: the rates of I and D, printed to 0.01°/yr", angle_rates, 0.0050001, "°/yr");
  }

  // ---- above the geographic poles (report §1.4) ---------------------------------------------------------------------------
  {
    const double year = rv("poles.year"), r = rv("poles.radius");
    double comps = 0;
    std::array<double, 2> meridian{};
    std::array<double, 2> implied{};
    for (int i = 0; i < 2; i++) {
      const std::string pole = i == 0 ? "north" : "south";
      const double sign = i == 0 ? 1 : -1;
      const GeocentricField g = geocentric_field(c, {0.0, sign * pi / 2, r}, year);
      comps = std::max({comps, std::fabs(g.north_T / nT - rv("poles." + pole + ".Xprime")), std::fabs(g.east_T / nT - rv("poles." + pole + ".Yprime")),
                        std::fabs(g.down_T / nT - rv("poles." + pole + ".Zprime"))});
      // With the pole given λ = 0, its north is −x̂ at the north pole and +x̂ at the south, its east +ŷ:
      // the horizontal field points to the meridian atan2(Y′, ∓X′).
      meridian[i] = std::atan2(g.east_T, -sign * g.north_T) / deg;
      implied[i] = std::atan2(rv("poles." + pole + ".Yprime"), -sign * rv("poles." + pole + ".Xprime")) / deg;
    }
    s.check("above the poles at r = a, 2025.0: X′, Y′, Z′ as the report gives them (to its 0.1 nT)", comps, 0.0500001, "nT");
    s.check("above the poles: the horizontal field's meridians agree with those of the report's components", std::max(std::fabs(meridian[0] - implied[0]), std::fabs(meridian[1] - implied[1])), 0.01, "°");
    s.check("above the south pole: the field in the half-plane of the 31.62° W meridian", meridian[1] - rv("poles.south.meridian"), 0.0050001, "°");
    // The report prints 166.27° W; its own X′ = 1709.5 nT and Y′ = 417.8 nT, in the frame it defines, put the field east of Greenwich.
    s.check("above the north pole: the 166.27° meridian, on the east as the report's own components put it (it prints W)", std::fabs(meridian[0]) - std::fabs(rv("poles.north.meridian")), 0.0050001, "°");
    // Continuity: the east component just off the pole tends to its value at the pole.
    const GeocentricField at = geocentric_field(c, {0.3, pi / 2, r}, year), near = geocentric_field(c, {0.3, pi / 2 - 1e-8, r}, year);
    s.check("at the pole the east component is the limit of its values beside it", (near.east_T - at.east_T) / nT, 1e-3, "nT");
    const Field top = field(c, Geodetic{pi / 2, 1.0, 0.0}, 2026.0);
    s.require("a geodetic place at the pole has a finite field", std::isfinite(top.elements.declination_rad) && std::isfinite(top.elements.north_T) && std::isfinite(top.itrs_T[0]));
  }

  // ---- the field as the gradient of its potential, and the rates as derivatives in time -----------------------------------
  {
    calx::random::Stream rnd(20261028, 2, 0, 0);
    double grad = 0, time = 0, elements = 0, angle_rates = 0;
    for (int i = 0; i < 20; i++) {
      const Geodetic place{std::asin(rnd.uniform(-0.98, 0.98)), rnd.uniform(-pi, pi), rnd.uniform(-1000, 850e3)};
      const double year = rnd.uniform(2025, 2030);
      const Field f = field(c, place, year);
      const physis::geometry::Vec3 x = wgs84().cartesian(place);
      const double h = 10;  // m
      for (int k = 0; k < 3; k++) {
        physis::geometry::Vec3 a = x, b = x;
        a[k] += h;
        b[k] -= h;
        const double dv = (potential_T_m(c, geocentric_position(a), year) - potential_T_m(c, geocentric_position(b), year)) / (2 * h);
        grad = std::max(grad, std::fabs(-dv - f.itrs_T[k]) / nT);
      }
      const double dt = 1e-3;
      const Field later = field(c, place, year + dt), earlier = field(c, place, year - dt);
      for (int k = 0; k < 3; k++) time = std::max(time, std::fabs((later.itrs_T[k] - earlier.itrs_T[k]) / (2 * dt) - f.itrs_T_per_year[k]) / nT);
      const auto rate = [&](auto get) { return (get(later.elements) - get(earlier.elements)) / (2 * dt); };
      elements = std::max({elements, std::fabs(rate([](const MagneticElements& e) { return e.horizontal_T; }) - f.secular_variation.horizontal_T_per_year) / nT,
                           std::fabs(rate([](const MagneticElements& e) { return e.total_T; }) - f.secular_variation.total_T_per_year) / nT});
      angle_rates = std::max({angle_rates, std::fabs(rate([](const MagneticElements& e) { return e.declination_rad; }) - f.secular_variation.declination_rad_per_year),
                              std::fabs(rate([](const MagneticElements& e) { return e.inclination_rad; }) - f.secular_variation.inclination_rad_per_year)});
    }
    s.check("the field is −∇V, by central differences of the potential (20 places, 1 km below to 850 km above)", grad, 1e-3, "nT");
    s.check("the rate of the Earth-fixed vector is its derivative in time", time, 1e-6, "nT/yr");
    s.check("the rates of H and F (eq. 20) are their derivatives in time", elements, 1e-5, "nT/yr");
    s.check("the rates of D and I (eq. 20) are their derivatives in time", angle_rates, 1e-9, "rad/yr");
    // The greatest secular variation at the surface, on which the note's bound on the time scale's
    // effect rests (§5.3); the rates of the components do not depend on the time.
    double fastest = 0;
    for (int la = -89; la <= 89; la++)
      for (int lo = -180; lo < 180; lo++) {
        const physis::geometry::Vec3 v = field(c, Geodetic{la * deg, lo * deg, 0.0}, 2025.0).itrs_T_per_year;
        fastest = std::max(fastest, physis::geometry::norm(v) / nT);
      }
    s.check("the secular variation on the ellipsoid stays below 230 nT/yr (its greatest on a 1° grid)", fastest, 230, "nT/yr");
  }

  // ---- the centred dipole (report §1.6), and L ------------------------------------------------------------------------------
  {
    const Dipole d = dipole(c, rv("dipole.year"));
    const double north_lon = d.north_pole_longitude_rad / deg, north_lat = d.north_pole_latitude_rad / deg, north_geodetic = d.north_pole_geodetic_latitude_rad / deg;
    s.check("the geomagnetic north pole at 2025.0: 72.76° W, 80.79° N geocentric, 80.85° N geodetic",
            std::max({std::fabs(north_lon - rv("dipole.north.longitude")), std::fabs(north_lat - rv("dipole.north.geocentric_latitude")),
                      std::fabs(north_geodetic - rv("dipole.north.geodetic_latitude"))}), 0.0050001, "°");
    // The south pole is the antipode.
    s.check("the geomagnetic south pole, the antipode: 107.24° E, 80.79° S geocentric, 80.85° S geodetic",
            std::max({std::fabs(angle_difference((north_lon + 180) * deg, rv("dipole.south.longitude") * deg) / deg), std::fabs(-north_lat - rv("dipole.south.geocentric_latitude")),
                      std::fabs(-north_geodetic - rv("dipole.south.geodetic_latitude"))}), 0.0050001, "°");
    s.check("the dipole's axis is inclined at 9.21° to the rotation axis", d.tilt_rad / deg - rv("dipole.tilt"), 0.0050001, "°");
    // The degree-1 model is the field of the moment m at the centre: μ₀/(4π r³) (3 (m · r̂) r̂ − m).
    const Coefficients one = c.truncated(1);
    calx::random::Stream rnd(20261028, 3, 0, 0);
    double worst = 0, constancy = 0;
    const physis::geometry::Vec3 m = d.moment_vector_A_m2();
    for (int i = 0; i < 20; i++) {
      const Geodetic place{std::asin(rnd.uniform(-0.98, 0.98)), rnd.uniform(-pi, pi), rnd.uniform(0, 2e7)};
      const physis::geometry::Vec3 x = wgs84().cartesian(place);
      const Field f = field(one, place, 2025.0);
      const double r = physis::geometry::norm(x), k = MAGNETIC_CONSTANT / (4 * pi * r * r * r);
      const double mu = physis::geometry::dot(m, x) / r;
      for (int j = 0; j < 3; j++) worst = std::max(worst, std::fabs(k * (3 * mu * x[j] / r - m[j]) - f.itrs_T[j]) / nT);
      // L is constant along the dipole's field lines: B · ∇L = 0.
      const physis::geometry::Vec3 b = d.field_T(x);
      physis::geometry::Vec3 gradient{};
      for (int j = 0; j < 3; j++) {
        physis::geometry::Vec3 up = x, down = x;
        up[j] += 10;
        down[j] -= 10;
        gradient[j] = (mcilwain_l(d, up) - mcilwain_l(d, down)) / 20;
      }
      constancy = std::max(constancy, std::fabs(physis::geometry::dot(b, gradient)) / (physis::geometry::norm(b) * physis::geometry::norm(gradient)));
    }
    s.check("the degree-1 model is the field of the centred dipole's moment 4π a³ g/μ₀ (20 places to 20 000 km)", worst, 1e-8, "nT");
    s.check("McIlwain's L of the centred dipole is constant along its field lines (B · ∇L = 0)", constancy, 1e-6);
    const physis::geometry::Vec3 east_of_axis = physis::geometry::unit(physis::geometry::cross(d.north_axis, physis::geometry::Vec3{0, 0, 1}));
    const physis::geometry::Vec3 on_equator{2 * c.reference_radius_m * east_of_axis[0], 2 * c.reference_radius_m * east_of_axis[1], 2 * c.reference_radius_m * east_of_axis[2]};
    s.check("L is 2 on the geomagnetic equator at r = 2a", mcilwain_l(d, on_equator) - 2, 1e-12);
    s.check("the dipole field on its equator at r = a is B₀, the moment's μ₀ M/(4π a³)",
            std::max(std::fabs(physis::geometry::norm(d.field_T(on_equator)) * 8 - d.equatorial_field_T),
                     std::fabs(MAGNETIC_CONSTANT * d.moment_A_m2 / (4 * pi * std::pow(c.reference_radius_m, 3)) - d.equatorial_field_T)) / nT, 1e-9, "nT");
  }

  // ---- validity, the compass zones and the error model ---------------------------------------------------------------------
  {
    const Geodetic paris{48.85 * deg, 2.35 * deg, 35.0};
    const Validity ok = validity(c, paris, 2026.5);
    s.require("valid in Paris in 2026.5, the compass reliable, nothing to say", ok.valid() && ok.declination_reliable() && ok.reasons.empty(), ok.reasons);
    s.require("valid from 2025.0 to 2030.0 inclusive, and not beyond",
              validity(c, paris, 2025.0).in_time && validity(c, paris, 2030.0).in_time && !validity(c, paris, 2024.99).in_time && !validity(c, paris, 2030.01).in_time);
    s.require("specified from 1 km below to 850 km above the ellipsoid",
              validity(c, Geodetic{paris.latitude, paris.longitude, -1000}, 2026).in_height && validity(c, Geodetic{paris.latitude, paris.longitude, 850e3}, 2026).in_height &&
                  !validity(c, Geodetic{paris.latitude, paris.longitude, -1001}, 2026).in_height && !validity(c, Geodetic{paris.latitude, paris.longitude, 850001}, 2026).in_height);
    const Validity late_high = validity(c, Geodetic{paris.latitude, paris.longitude, 9e5}, 2031);
    s.require("and says why not", !late_high.valid() && late_high.reasons.find("2031.0000") != std::string::npos && late_high.reasons.find("900.000 km") != std::string::npos, late_high.reasons);

    // The error model's δD where H is greatest on the surface, which the report states; and that
    // greatest H itself, found over the ellipsoid at 2025.0 (a 1° grid, then a local climb).
    double best = 0, blat = 0, blon = 0;
    for (int la = -60; la <= 60; la++)
      for (int lo = -180; lo < 180; lo++) {
        const double hh = field(c, Geodetic{la * deg, lo * deg, 0.0}, 2025.0).elements.horizontal_T;
        if (hh > best) { best = hh; blat = la * deg; blon = lo * deg; }
      }
    for (double step = 0.5 * deg; step > 1e-7; step /= 2)
      for (bool moved = true; moved;) {
        moved = false;
        for (const auto& [dla, dlo] : {std::pair{step, 0.0}, {-step, 0.0}, {0.0, step}, {0.0, -step}}) {
          const double hh = field(c, Geodetic{blat + dla, blon + dlo, 0.0}, 2025.0).elements.horizontal_T;
          if (hh > best) { best = hh; blat += dla; blon += dlo; moved = true; }
        }
      }
    s.check("the greatest H on the ellipsoid at 2025.0 is the report's 41875 nT (§3.4)", best / nT - rv("uncertainty.horizontal"), 0.5000001, "nT");
    s.check("WMM2025's error model gives δD = 0.29° there (eq. 43)", wmm2025_uncertainty(rv("uncertainty.horizontal") * nT).declination_rad / deg - rv("uncertainty.declination"), 0.0050001, "°");
    s.require("δD grows without bound towards a dip pole", wmm2025_uncertainty(100 * nT).declination_rad > 50 * deg);
  }
  return s.done();
}
