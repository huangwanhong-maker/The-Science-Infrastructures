// A domain's functions from C++: where the Sun stands, now, computed by the
// astronomy domain (domains/astronomy/js/sun.js) through the service's compute
// bridge, so that C++ runs the algorithm the pages run.
//
//   sunPosition(jd)            the Sun's apparent right ascension and declination
//                              (degrees, equinox of date), ecliptic longitude and
//                              distance (au), from the Astronomical Almanac's
//                              low-precision formulae: 0.01° from 1950 to 2050;
//   subsolarPoint(ms)          where the Sun is at the zenith: longitude RA - GMST,
//                              latitude the declination (UT1 taken as UTC);
//   sunAltAz(lon, lat, sub)    its altitude and azimuth (from north through east)
//                              at a place, parallax neglected;
//   daylight(altitude)         day, civil, nautical or astronomical twilight, or night.
//
// The Julian date is the domain's julianDate(ms), of UTC: JD = ms / 86 400 000 + 2 440 587.5.
//
//   api/build/cpp-sun [lon lat] [url]

#include "science.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv) {
  double lon = argc > 2 ? std::atof(argv[1]) : 139.69, lat = argc > 2 ? std::atof(argv[2]) : 35.69;
  try {
    science::Client api(argc > 3 ? argv[3] : "");
    auto now = std::chrono::system_clock::now();
    double ms = std::chrono::duration<double, std::milli>(now.time_since_epoch()).count();  // UTC, since 1970

    double jd = api.call("astronomy", "julianDate", ms).as_number();
    science::Json sun = api.call("astronomy", "sunPosition", jd);
    science::Json sub = api.call("astronomy", "subsolarPoint", ms);           // {lon, lat}
    science::Json sky = api.call("astronomy", "sunAltAz", lon, lat, sub);      // {altitude, azimuth}
    std::string light = api.call("astronomy", "daylight", sky["altitude"]).text();

    std::printf("%s (JD %.5f UTC)\n", science::utc_text(now).c_str(), jd);
    std::printf("  the Sun: right ascension %.3f°, declination %+.3f° (equinox of date), %.5f au away\n",
                sun["ra"].as_number(), sun["dec"].as_number(), sun["distanceAu"].as_number());
    std::printf("  overhead at %.2f°, %+.2f°\n", sub["lon"].as_number(), sub["lat"].as_number());
    std::printf("  at %.2f°, %+.2f°: altitude %+.1f°, azimuth %.1f° (%s)\n", lon, lat, sky["altitude"].as_number(),
                sky["azimuth"].as_number(), light.c_str());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
  return 0;
}
