// Air parcels in C++: a Trajectories plug-in (api/cpp/plugin.hpp) that the server
// compiles and runs, the counterpart of trajectories.js. The parcels start on a
// ring around the chosen place (or, with no place chosen or when asked, on a
// lattice over the globe) and are carried through the wind chosen, at its own
// times, by the header's fourth-order Runge–Kutta step on the sphere,
//
//   dλ/dt = u / (a cos φ),   dφ/dt = v / a.
//
// The result is a sketch: each path, and a point on it that moves with the time
// shown; with a place, a marker and a label there; and a note of how far the
// parcels went. A point of the ring at distance s and azimuth θ (from north,
// clockwise) from the place (λ1, φ1), with δ = s / a, is
//
//   φ = asin(sin φ1 cos δ + cos φ1 sin δ cos θ),
//   λ = λ1 + atan2(sin θ sin δ cos φ1, cos δ − sin φ1 sin φ),
//
// (Snyder 1987, eqs. 5-5 and 5-6, the direct problem on the sphere).
//
// Check (applications/earth-atlas/py/check.py): a parcel in a solid-body rotation
// u = U cos φ keeps its latitude and turns by U t / a.

#include "plugin.hpp"

using plugin::DEG;
using plugin::EARTH_RADIUS_M;

// The point at distance s (m) and azimuth θ (radians from north) from (lon, lat), in degrees.
static plugin::LonLat destination(double lon, double lat, double s, double azimuth) {
  const double d = s / EARTH_RADIUS_M, p1 = lat * DEG;
  const double p = std::asin(std::sin(p1) * std::cos(d) + std::cos(p1) * std::sin(d) * std::cos(azimuth));
  const double l = lon * DEG + std::atan2(std::sin(azimuth) * std::sin(d) * std::cos(p1), std::cos(d) - std::sin(p1) * std::sin(p));
  return {plugin::wrap_longitude(l / DEG), p / DEG};
}

// The great-circle distance between two places, in metres (the haversine form, as domains/geodesy).
static double distance(double lon1, double lat1, double lon2, double lat2) {
  const double dp = (lat2 - lat1) * DEG, dl = (lon2 - lon1) * DEG;
  const double h = std::sin(dp / 2) * std::sin(dp / 2) + std::cos(lat1 * DEG) * std::cos(lat2 * DEG) * std::sin(dl / 2) * std::sin(dl / 2);
  return 2 * EARTH_RADIUS_M * std::asin(std::min(1.0, std::sqrt(h)));
}

struct Parcels : plugin::Trajectories {
  Parcels() {
    id = "parcels-cpp";
    title = "Air parcels (C++)";
    description = "Parcels carried through the wind chosen from a ring around the chosen place, or from a lattice over the globe, "
                  "by a Runge-Kutta step on the sphere in a C++ program in the server.";
    needs = {{"wind", "Wind", {"eastward_wind", "northward_wind"}, "", 2}};
    parameters = {{"hours", "Hours ahead", 48, 6, 240, 6, "h"},
                  {"step", "Time step", 1, 0.25, 6, 0.25, "h"},
                  {"count", "Parcels", 24, 1, 200, 1, ""},
                  {"radius", "Ring radius", 300, 0, 3000, 10, "km"},
                  plugin::choice("start", "Start from", "ring", {{"ring", "a ring around the chosen place"}, {"lattice", "a lattice over the globe"}})};
  }

  std::vector<plugin::Seed> seeds(const plugin::Inputs&, const plugin::Params& params, const plugin::Context& context) override {
    const int n = std::max(1, int(std::lround(params.number("count", 24))));
    std::vector<plugin::Seed> out;
    if (!context.pinned || params.choice("start", "ring") == "lattice") {
      for (const auto& p : plugin::lattice(n)) out.push_back({p.lon, p.lat, ""});
      return out;
    }
    const double radius = params.number("radius", 300) * 1000;
    for (int k = 0; k < n; k++) {
      const auto p = destination(context.pinned->lon, context.pinned->lat, radius, 2 * plugin::PI * k / n);
      out.push_back({p.lon, p.lat, ""});
    }
    return out;
  }

  plugin::Results run(const plugin::Inputs& inputs, const plugin::Params& params, const plugin::Context& context) override {
    plugin::Results r = Trajectories::run(inputs, params, context);
    plugin::Sketch& sketch = r.sketches.front();
    if (context.pinned) {
      sketch.markers.push_back({context.pinned->lon, context.pinned->lat, "#f6c47e", 4});
      sketch.labels.push_back({context.pinned->lon, context.pinned->lat, "start"});
    }
    double total = 0, most = 0;
    for (const auto& path : sketch.polylines) {
      const auto& c = path.coords;
      const double d = distance(c[0], c[1], c[c.size() - 2], c[c.size() - 1]);
      total += d;
      most = std::max(most, d);
    }
    if (!sketch.polylines.empty()) {
      char text[200];
      std::snprintf(text, sizeof text, "From start to end the parcels moved %.0f km on average, and at most %.0f km (great-circle distances).",
                    total / sketch.polylines.size() / 1000, most / 1000);
      r.notes.push_back(text);
    }
    return r;
  }
};

PLUGIN_MAIN(Parcels)
