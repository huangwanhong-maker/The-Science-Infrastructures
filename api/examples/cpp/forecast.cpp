// The forecast at a place over time, from C++: GFS's temperature at 2 m and its
// wind at 10 m at every time of the product, with the direction the wind comes
// from given by the meteorology domain's own function, called through the
// service (one tested implementation for every language).
//
//   make -C api && api/build/cpp-forecast [lon lat] [url]     (the service: $SCIENCE_API_URL or http://127.0.0.1:8800)

#include "science.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv) {
  double lon = argc > 2 ? std::atof(argv[1]) : 139.69;  // degrees east, WGS 84 (Tokyo by default)
  double lat = argc > 2 ? std::atof(argv[2]) : 35.69;
  try {
    science::Client api(argc > 3 ? argv[3] : "");
    science::Field gfs = api.field("gfs-forecast");
    // A series is the value at the place at every stored time: bilinear between
    // the four grid points around it, as the Earth Atlas samples it.
    science::Series t = gfs.series("tmp2m", lon, lat);    // K
    science::Series u = gfs.series("ugrd10m", lon, lat);  // m s-1, toward the east
    science::Series v = gfs.series("vgrd10m", lon, lat);  // m s-1, toward the north
    std::printf("GFS at %.2f, %.2f (times UTC)\n", lon, lat);
    for (std::size_t k = 0; k < t.times.size(); k++) {
      double speed = std::hypot(u.values[k], v.values[k]);
      // Meteorology's convention: the direction the wind blows from, clockwise from north.
      double from = api.call("meteorology", "windFromDirection", u.values[k], v.values[k]).as_number();
      std::printf("  %.10s %.5s  %6.1f °C   wind %4.1f m/s from %3.0f°\n", t.times[k].c_str(), t.times[k].c_str() + 11,
                  t.values[k] - 273.15, speed, from);
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
  return 0;
}
