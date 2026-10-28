// Checks of the geoid (geodesy/geoid): EGM96's grid read and interpolated.
//
// References: PROJ 9.8's bilinear interpolation of NGA's EGM96 15′ grid (tests/fixtures/
// geoid-reference.txt, made by tests/make_geoid_fixtures.py), which the bilinear interpolation of
// GeographicLib's copy of the same grid must meet to the copy's 3 mm quantum; the grid's own nodes;
// the continuity of the cubic interpolation across the prime meridian and at the poles. The grid is
// the data library's product geoid/egm96-15.pgm (data/products/geoid_grid.py), read from the cache;
// without it the checks of the grid are skipped.

#include <cmath>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>

#include <calx/verification/check.hpp>
#include <geodesy/geoid/geoid.hpp>

using namespace geodesy::geoid;
using calx::verification::Suite;

int main() {
  Suite s;
  const double deg = std::numbers::pi / 180;
  std::ifstream probe(GEODESY_GEOID_GRID);
  if (!probe) {
    s.skip("the EGM96 grid is not in the cache (python data/acquisition/fetch.py egm96-geoid; python data/products/geoid_grid.py)");
    return s.done();
  }
  const GeoidGrid g = GeoidGrid::read_pgm(GEODESY_GEOID_GRID);
  s.require("the grid is 1440 × 721 nodes at 15′", g.width() == 1440 && g.height() == 721 && std::fabs(g.step_degrees() - 0.25) < 1e-12);

  // ---- against PROJ -------------------------------------------------------------------------------------
  {
    std::ifstream in(GEODESY_FIXTURES "/geoid-reference.txt");
    std::string line;
    double worst = 0, cubic_gap = 0;
    int n = 0;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream f(line);
      double lon, lat, N;
      f >> lon >> lat >> N;
      worst = std::max(worst, std::fabs(g.undulation(lat * deg, lon * deg, Interpolation::bilinear) - N));
      cubic_gap = std::max(cubic_gap, std::fabs(g.undulation(lat * deg, lon * deg, Interpolation::cubic) - N));
      n++;
    }
    s.require("the reference file was read", n > 400);
    s.check("EGM96 bilinear against PROJ's (408 points; the copy's quantum is 3 mm)", worst, 2e-3, "m");
    s.check("EGM96 cubic stays near the bilinear (the two differ between nodes only)", cubic_gap, 1.0, "m");
    s.check("EGM96 at 0° N 0° E: 17.16 m", g.undulation(0, 0) - 17.16, 0.01, "m");
  }

  // ---- the nodes, the seams and the poles ------------------------------------------------------------------
  {
    double node = 0;
    for (const auto& [r, c] : {std::pair{100, 200}, {360, 0}, {500, 1439}, {1, 720}}) {
      const double lat = (90 - 0.25 * r) * deg, lon = 0.25 * c * deg;
      node = std::max({node, std::fabs(g.undulation(lat, lon, Interpolation::bilinear) - g.node(r, c)),
                       std::fabs(g.undulation(lat, lon, Interpolation::cubic) - g.node(r, c))});
    }
    s.check("both interpolations pass through the nodes", node, 1e-9, "m");
    const double east = g.undulation(12.3 * deg, 1e-9 * deg), west = g.undulation(12.3 * deg, -1e-9 * deg);
    s.check("the cubic is continuous across the prime meridian", east - west, 1e-6, "m");
    double spread = 0;
    const double at0 = g.undulation(89.9999 * deg, 0);
    for (double lon : {45.0, 133.0, -77.0, 179.0}) spread = std::max(spread, std::fabs(g.undulation(89.9999 * deg, lon * deg) - at0));
    s.check("near the north pole the undulation hardly depends on the longitude", spread, 0.01, "m");
    s.check("the orthometric height is h − N", g.orthometric_height(0.5, 1.0, 100.0) - (100.0 - g.undulation(0.5, 1.0)), 1e-12, "m");
  }
  return s.done();
}
