// A grid read as float32 from C++, and its mean: GFS's temperature at 2 m over a
// box at an instant, its mean over the points, and its mean over the area.
//
// The two differ because the points of a regular longitude-latitude grid are
// not equally spaced on the Earth: a cell of Δλ by Δφ at latitude φ covers
// R² Δλ (sin φ₂ − sin φ₁) ≈ R² Δλ Δφ cos φ of a sphere of radius R, so its share
// of the area goes as cos φ. The mean over the area weights each point by cos φ;
// over a box far from the equator, or a whole globe, the mean over the points
// gives the cold high latitudes too much weight.
//
//   api/build/cpp-grid_mean [variable] [west south east north] [time] [url]
//   (by default tmp2m over Japan, 122 24 150 46, now)

#include "science.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv) {
  const std::string var = argc > 1 ? argv[1] : "tmp2m";
  science::Box box{122, 24, 150, 46};
  if (argc > 5) box = {std::atof(argv[2]), std::atof(argv[3]), std::atof(argv[4]), std::atof(argv[5])};
  const std::string time = argc > 6 ? argv[6] : "now";
  try {
    science::Client api(argc > 7 ? argv[7] : "");
    // Rows north to south, row-major, NaN where there is no value; between two
    // stored times, the linear blend of the two, as the Earth Atlas draws it.
    science::Grid g = api.field("gfs-forecast").grid(var, time, box);
    const double rad = std::acos(-1.0) / 180.0;
    double sum = 0, weighted = 0, weights = 0, lowest = INFINITY, highest = -INFINITY;
    int points = 0;
    for (int row = 0; row < g.nlat; row++) {
      double w = std::cos(g.lat(row) * rad);
      for (int col = 0; col < g.nlon; col++) {
        double x = g.at(row, col);
        if (std::isnan(x)) continue;
        sum += x;
        points++;
        weighted += w * x;
        weights += w;
        lowest = std::fmin(lowest, x);
        highest = std::fmax(highest, x);
      }
    }
    std::printf("%s over %g..%g E, %g..%g N at %s: %d x %d points, %d with a value\n", var.c_str(), box.west, box.east, box.south,
                box.north, g.time.c_str(), g.nlon, g.nlat, points);
    if (points == 0) return 0;
    std::printf("  mean over the points  %.3f %s\n", sum / points, g.units.c_str());
    std::printf("  mean over the area    %.3f %s   (weights cos φ)\n", weighted / weights, g.units.c_str());
    std::printf("  from %.3f to %.3f %s\n", lowest, highest, g.units.c_str());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
  return 0;
}
