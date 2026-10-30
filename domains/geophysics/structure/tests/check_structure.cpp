// Checks of the structure module (geophysics/structure): the crust (CRUST1.0), PREM, the interior
// they compose, and sections along geodesics.
//
// References:
// - PREM as ObsPy's TauP tabulates it (tests/fixtures/prem-obspy.txt, written by
//   tests/make_prem_fixtures.py from obspy/taup/data/prem.nd at a fixed commit);
// - the Earth's mass, GM of WGS 84 (NGA.STND.0036, Table 3.1) over G (CODATA 2018), which PREM was
//   fitted to;
// - the same checks of the crust as the domain's JavaScript (tests/run-tests.mjs): a constructed
//   oceanic column, the cells' corners, Tibet, the central Pacific and the Canadian Shield from the
//   data library's product crust-crust1 (skipped, and said, when it is not in the cache);
// - sections against geodesy's direct problem: a place built at a known distance along and across.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

#include <calx/verification/check.hpp>
#include <geodesy/geodesic/geodesic.hpp>
#include <geodesy/gravity/normal.hpp>
#include <geophysics/structure/crust.hpp>
#include <geophysics/structure/interior.hpp>
#include <geophysics/structure/prem.hpp>
#include <geophysics/structure/section.hpp>

using namespace geophysics::structure;
using calx::verification::Suite;

namespace {

constexpr double pi = std::numbers::pi, deg = pi / 180;

std::string cache_path(const std::string& relative) {
  const char* env = std::getenv("SCIENCE_CACHE");
  return (env ? std::string(env) : std::string(GEOPHYSICS_ROOT "/../../data/cache")) + "/" + relative;
}

struct PremRow {
  double depth_km, vp, vs, rho;
};

std::vector<PremRow> prem_rows(const std::string& path) {
  std::vector<PremRow> out;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream f(line);
    PremRow r{};
    if (f >> r.depth_km >> r.vp >> r.vs >> r.rho) out.push_back(r);
  }
  return out;
}

}  // namespace

int main() {
  Suite s;

  // ---- PREM ---------------------------------------------------------------------------------------------------------
  {
    const auto rows = prem_rows(GEOPHYSICS_FIXTURES "/prem-obspy.txt");
    s.require("PREM: ObsPy's table read (88 rows)", rows.size() == 88, std::to_string(rows.size()));
    // ObsPy's table is not linear in the radius between 24.4 and 220 km, where PREM's isotropic lid and
    // low-velocity zone are: there the two depart by up to 1.2·10⁻³ km/s (they agree at both ends), and
    // are compared apart.
    double dv = 0, dv_upper = 0, dr = 0;
    int compared = 0, sides = 0, matched = 0;
    for (std::size_t i = 0; i < rows.size(); i++) {
      const PremRow& r = rows[i];
      if (r.depth_km < 3) continue;  // the table's crust reaches the surface, where PREM has its ocean
      // Two rows of one depth: the upper side first.
      const bool upper = i + 1 < rows.size() && rows[i + 1].depth_km == r.depth_km;
      if (upper) {
        sides++;
        for (const auto& d : PREM_DISCONTINUITIES)
          if (std::fabs(d.depth_m - r.depth_km * 1e3) < 1) matched++;
      }
      const PremValues p = prem_at_depth(r.depth_km * 1e3, upper);
      const double e = std::max(std::fabs(p.vp_m_s / 1e3 - r.vp), std::fabs(p.vs_m_s / 1e3 - r.vs));
      if (r.depth_km > 24.4 && r.depth_km < 220) dv_upper = std::max(dv_upper, e);
      else dv = std::max(dv, e);
      dr = std::max(dr, std::fabs(p.density_kg_m3 / 1e3 - r.rho));
      compared++;
    }
    s.require("PREM: rows compared from 3 km to the centre", compared == 87, std::to_string(compared));
    s.check("PREM: vp and vs against ObsPy's table, both sides of each discontinuity, outside 24.4–220 km (PREM's coefficients have four decimals)", dv, 2e-4, "km/s");
    s.check("PREM: vp and vs in the lid and the low-velocity zone, 24.4–220 km, where the table is not linear in the radius", dv_upper, 1.5e-3, "km/s");
    s.check("PREM: density against ObsPy's table", dr, 2e-4, "g/cm³");
    s.require("PREM: every discontinuity of the table is one of PREM's (15, 24.4, 220, 400, 670, 2891 and 5149.5 km)", sides == 7 && matched == 7,
              std::to_string(sides) + " sides, " + std::to_string(matched) + " matched");
    // The boundaries at 771 and 2741 km are not discontinuities: the polynomials meet there.
    double jump = 0;
    for (const double r_km : {5600.0, 3630.0}) {
      const PremValues a = prem_at_radius(r_km * 1e3, true), b = prem_at_radius(r_km * 1e3, false);
      jump = std::max({jump, std::fabs(a.vp_m_s - b.vp_m_s), std::fabs(a.vs_m_s - b.vs_m_s), std::fabs(a.density_kg_m3 - b.density_kg_m3)});
    }
    s.check("PREM: continuous at 771 and 2741 km, where its shells meet without a jump", jump, 1.0, "m/s or kg/m³");
    s.require("PREM: the outer core is fluid (no S waves), the inner core solid",
              prem_at_depth(2891e3).vs_m_s == 0 && prem_at_depth(3500e3).vs_m_s == 0 && prem_at_depth(5149.5e3, true).vs_m_s == 0 &&
                  prem_at_depth(5149.5e3).vs_m_s > 3000 && prem_at_depth(2891e3, true).vs_m_s > 7000);
    // The mass PREM was fitted to: GM of WGS 84 over G of CODATA 2018 (G has changed by 3·10⁻⁴ since 1981).
    const double gm = geodesy::gravity::LevelEllipsoid::wgs84().GM(), mass = prem_mass_within(PREM_RADIUS_M);
    s.check("PREM: its mass, the densities integrated exactly, against GM / G (relative)", std::fabs(mass - gm / GRAVITATIONAL_CONSTANT) / mass, 1e-3);
    s.require("PREM: the mass within a radius grows to the whole and vanishes at the centre",
              prem_mass_within(0) == 0 && prem_mass_within(3480e3) < prem_mass_within(6000e3) && prem_mass_within(7000e3) == mass);
    // Gravity peaks at the core–mantle boundary: the dense core is a third of the mass in a sixth of the volume.
    double best_r = 0, best_g = 0;
    for (double r = 1e3; r <= PREM_RADIUS_M; r += 1e3)
      if (prem_gravity(r) > best_g) best_g = prem_gravity(r), best_r = r;
    s.require("PREM: gravity is greatest at the core–mantle boundary (3480 km from the centre)", std::fabs(best_r - 3480e3) <= 1e3,
              std::to_string(best_r / 1e3) + " km, " + std::to_string(best_g) + " m/s²");
    s.check("PREM: gravity at its surface, G M / a², against GM / a² (relative)", std::fabs(prem_gravity(PREM_RADIUS_M) - gm / (PREM_RADIUS_M * PREM_RADIUS_M)) / prem_gravity(PREM_RADIUS_M), 1e-3);
    const double factor = prem_moment_of_inertia() / (mass * PREM_RADIUS_M * PREM_RADIUS_M);
    s.require("PREM: its moment of inertia below a uniform sphere's 0.4 M a² (the mass is concentrated toward the centre)", factor > 0.3 && factor < 0.4,
              std::to_string(factor));
  }

  // ---- the crust ----------------------------------------------------------------------------------------------------
  {
    // A constructed oceanic column: 4 km of water, 0.5 km of sediment, 7 km of crust (no ice, no middle or
    // lower sediments, no middle crust), as the JavaScript's check.
    const double tops[9] = {0, -4, -4, -4.5, -4.5, -4.5, -6, -6, -11.5}, vp[9] = {1.5, 0, 2, 0, 0, 5, 0, 7, 8.1};
    const CrustColumn c = crust_column([&](CrustQuantity q, int k) { return q == CrustQuantity::top ? tops[k] : q == CrustQuantity::vp ? vp[k] : 0.0; });
    std::string names;
    for (const auto& l : c.layers) names += std::string(l.name) + ",";
    s.require("crust: a constructed oceanic column, layer by layer",
              names == "water,upper sediments,upper crystalline crust,lower crystalline crust," && std::fabs(c.moho_depth_m - 11.5e3) < 1e-9 &&
                  std::fabs(c.crust_m - 7.5e3) < 1e-9 && std::fabs(c.sediments_m - 0.5e3) < 1e-9 && c.oceanic && c.mantle.vp_m_s == 8100,
              names);
    s.require("crust: the layer that holds a height: above the sea, in the water, at a boundary the layer above, below the Moho",
              c.layer_at(10) == -1 && c.layer_at(-1000) == 0 && c.layer_at(-4000) == 0 && c.layer_at(-4200) == 1 && c.layer_at(-11500) == 3 &&
                  c.layer_at(-12000) == 4);
    const CrustGrid g;
    s.require("crust: the 1° cell containing a place (its corners and the antimeridian)",
              g.cell(90, -180) == 0 && g.cell(-89.99, 179.99) == 360 * 180 - 1 && g.cell(0.5, 0.5) == 89 * 360 + 180 && g.cell(0, 180) == g.cell(0, -180));
    s.require("crust: a cell's centre", g.centre(0)[0] == 89.5 && g.centre(0)[1] == -179.5 && g.centre(89 * 360 + 180)[0] == 0.5 && g.centre(89 * 360 + 180)[1] == 0.5);

    const std::string product = cache_path("products/crust-crust1.bin");
    if (std::filesystem::exists(product)) {
      const auto m = std::make_shared<const Crust1>(Crust1::read_product(product));
      const CrustColumn tibet = m->column(29.65, 91.1), pacific = m->column(10, -150), shield = m->column(55, -100);
      s.require("crust: Tibet over 60 km thick; the central Pacific oceanic, its Moho 10 to 14 km below sea level; the Canadian Shield 35 to 45 km",
                tibet.crust_m > 60e3 && pacific.oceanic && pacific.moho_depth_m > 10e3 && pacific.moho_depth_m < 14e3 && shield.crust_m > 35e3 &&
                    shield.crust_m < 45e3 && !shield.oceanic,
                std::to_string(tibet.crust_m) + " " + std::to_string(pacific.moho_depth_m) + " " + std::to_string(shield.crust_m));
      s.require("crust: P waves quicken to 7.6-8.4 km/s beneath the Moho everywhere checked (Pn)",
                tibet.mantle.vp_m_s > 7600 && tibet.mantle.vp_m_s < 8400 && pacific.mantle.vp_m_s > 7600 && pacific.mantle.vp_m_s < 8400 &&
                    shield.mantle.vp_m_s > 7600 && shield.mantle.vp_m_s < 8400);
      // Every cell: its layers contiguous from its top to its Moho, each thicker than nothing, the crust
      // between 4 and 80 km, and P waves faster beneath the Moho than in any layer above it.
      int bad = 0;
      double thinnest = 1e9, thickest = 0;
      for (int k = 0; k < Crust1::cells; k++) {
        const CrustColumn x = m->column_of_cell(k);
        bool ok = !x.layers.empty() && x.layers.front().top_m == x.surface_m && x.layers.back().bottom_m == -x.moho_depth_m;
        for (std::size_t i = 0; i < x.layers.size(); i++) {
          ok = ok && x.layers[i].thickness_m() > 0 && (i == 0 || x.layers[i].top_m == x.layers[i - 1].bottom_m) && x.layers[i].vp_m_s < x.mantle.vp_m_s;
        }
        thinnest = std::min(thinnest, x.crust_m), thickest = std::max(thickest, x.crust_m);
        if (!ok) bad++;
      }
      s.require("crust: all 64,800 cells contiguous from the top to the Moho, slower than the mantle beneath", bad == 0, std::to_string(bad) + " cells not");
      s.require("crust: its thickness from 4 to 80 km in every cell", thinnest > 4e3 && thickest <= 80e3,
                std::to_string(thinnest / 1e3) + " to " + std::to_string(thickest / 1e3) + " km");

      // The interior they compose, beneath the central Pacific.
      const InteriorValues air = interior_at(pacific, 100), water = interior_at(pacific, -1000), lid = interior_at(pacific, -pacific.moho_depth_m - 100);
      const double x = (PREM_RADIUS_M - pacific.moho_depth_m - 100) / PREM_RADIUS_M;
      s.require("interior: nothing in the air; CRUST1.0's water in the sea; PREM's lid beneath the shallow Moho, its polynomial continued upward",
                air.layer.empty() && water.crust1 && water.layer == "water" && !lid.crust1 && lid.layer == "lid (PREM, continued upward)" &&
                    std::fabs(lid.vp_m_s - (4.1875 + 3.9382 * x) * 1e3) < 1e-9 && std::fabs(lid.vs_m_s - (2.1519 + 2.3481 * x) * 1e3) < 1e-9);
      const auto profile = interior_profile(tibet);
      bool sorted = true;
      for (std::size_t i = 1; i < profile.size(); i++) sorted = sorted && profile[i].depth_m >= profile[i - 1].depth_m;
      bool cmb = false;
      for (std::size_t i = 1; i < profile.size(); i++)
        if (std::fabs(profile[i].depth_m - 2891e3) < 1 && std::fabs(profile[i - 1].depth_m - 2891e3) < 1) cmb = profile[i - 1].vs_m_s > 7000 && profile[i].vs_m_s == 0;
      s.require("interior: the profile beneath Tibet from its top to the centre, in order, the core–mantle boundary's two sides in it",
                sorted && profile.front().depth_m == -tibet.surface_m && profile.back().depth_m == PREM_RADIUS_M && cmb &&
                    std::fabs(profile.back().density_kg_m3 - 13088.5) < 1e-6,
                std::to_string(profile.size()) + " samples");
      const InteriorField rho(m, InteriorQuantity::density);
      const auto up = rho.at({0, {10 * deg, -150 * deg, 100}}), down = rho.at({0, {10 * deg, -150 * deg, -3000e3}}), beyond = rho.at({0, {10 * deg, -150 * deg, -6400e3}});
      s.require("interior as a field: refused in the air and beyond the centre, with the reasons; the outer core's density at 3000 km",
                !up.ok() && up.reason == "above the ground or the sea: in the air" && !beyond.ok() && beyond.reason == "beyond the centre of the Earth" && down.ok() &&
                    std::fabs(down.value - prem_at_depth(3000e3).density_kg_m3) < 1e-9 && down.note == "PREM, below CRUST1.0's Moho");
    } else {
      s.skip("crust: crust-crust1 is not in the cache (python data/acquisition/fetch.py crust1; python data/products/crust_model.py)");
    }
  }

  // ---- sections ----------------------------------------------------------------------------------------------------
  {
    const auto g = geodesy::geodesic::Geodesic::wgs84();
    const double lat0 = 38 * deg, lon0 = 140 * deg, az = 100 * deg;
    const SectionLine line(g, lat0, lon0, az);
    // A place 80 km to the right of the line 400 km along it, by the direct problem at a right angle.
    const SectionPoint p1 = line.at(400e3);
    const auto right = g.direct(p1.latitude, p1.longitude, p1.azimuth + pi / 2, 80e3);
    const SectionOffset r1 = line.offset(right.lat2, right.lon2);
    // One 250 km behind the start, on the line; one 120 km to the left 1200 km along.
    const SectionPoint behind = line.at(-250e3);
    const SectionOffset r2 = line.offset(behind.latitude, behind.longitude);
    const SectionPoint p3 = line.at(1200e3);
    const auto left = g.direct(p3.latitude, p3.longitude, p3.azimuth - pi / 2, 120e3);
    const SectionOffset r3 = line.offset(left.lat2, left.lon2);
    s.check("sections: along and across a place 80 km to the right 400 km along the geodesic (WGS 84)", std::fabs(r1.along_m - 400e3) + std::fabs(r1.across_m - 80e3), 1e-6, "m");
    s.check("sections: a place on the line 250 km behind its centre", std::fabs(r2.along_m + 250e3) + std::fabs(r2.across_m), 1e-6, "m");
    s.check("sections: a place 120 km to the left 1200 km along", std::fabs(r3.along_m - 1200e3) + std::fabs(r3.across_m + 120e3), 1e-6, "m");
    const Section sec = line.sample(500e3, 5e3);
    double step_error = 0;
    for (std::size_t i = 1; i < sec.points.size(); i++) {
      const auto& a = sec.points[i - 1];
      const auto& b = sec.points[i];
      step_error = std::max(step_error, std::fabs(g.inverse(a.latitude, a.longitude, b.latitude, b.longitude).s12 - 5e3));
    }
    s.require("sections: 201 samples from −500 to +500 km, the centre the place itself",
              sec.points.size() == 201 && sec.points[100].s_m == 0 && sec.points[100].latitude == lat0 && std::fabs(sec.points[100].longitude - lon0) < 1e-15 &&
                  sec.points.front().s_m == -500e3 && sec.points.back().s_m == 500e3);
    s.check("sections: successive samples a step apart by the inverse problem", step_error, 1e-6, "m");
    // Half the equator, a geodesic of the ellipsoid: from 90° W to 90° E.
    const SectionLine equator(g, 0, 0, pi / 2);
    const Section half = equator.sample(pi * g.a() / 2, pi * g.a() / 2);
    s.require("sections: half the equator from 90° W to 90° E", half.points.size() == 3 && std::fabs(half.points[0].longitude + pi / 2) < 1e-12 &&
                                                                     std::fabs(half.points[2].longitude - pi / 2) < 1e-12 && std::fabs(half.points[2].latitude) < 1e-12);
  }
  return s.done();
}
