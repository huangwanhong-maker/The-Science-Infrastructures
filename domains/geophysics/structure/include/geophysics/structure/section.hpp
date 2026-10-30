// A vertical section through the Earth along a geodesic of the ellipsoid: the places along it, and
// where other places (earthquakes) lie relative to it. The C++ of the domain's `structure` beneath
// js/section.js, on the ellipsoid where the JavaScript is on a sphere.
//
// The section is the geodesic through a centre at an azimuth, sampled every `step` from −half to
// +half (geodesy's geodesics, Karney 2013: the direct problem along one line). A place lies at
// distance `along` on the line from the centre and `across` from it, positive to the right of the
// direction of travel: the foot of the shortest geodesic from the place to the line, found by Brent's
// minimisation of the distance along the line (calx), refined by the right triangle at the foot,
// the sign from the azimuths there. Near
// the line (within some thousand kilometres) the foot is unique and the geodesic from it meets the
// line at a right angle; farther from the line, beyond a quarter of a meridian, it need not be, and
// a place is not asked for there.

#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <calx/optimisation/minimise.hpp>
#include <geodesy/geodesic/geodesic.hpp>

namespace geophysics::structure {

struct SectionPoint {
  double s_m = 0;                       // distance along the section from its centre (m), negative behind it
  double latitude = 0, longitude = 0;  // geodetic, WGS 84 (rad)
  double azimuth = 0;                  // of the section there (rad, clockwise from north)
};

struct Section {
  double latitude = 0, longitude = 0, azimuth = 0;  // the centre and the azimuth there (rad)
  double half_m = 0, step_m = 0;
  std::vector<SectionPoint> points;  // from −half to +half
};

// Where a place lies relative to a section: along it (m, from the centre) and across it (m, to the
// right positive), and the azimuth of the line at the foot (rad).
struct SectionOffset {
  double along_m = 0, across_m = 0, azimuth_at_foot = 0;
};

class SectionLine {
 public:
  // The section through (latitude, longitude) at `azimuth` (rad) on the ellipsoid of `g`.
  SectionLine(const geodesy::geodesic::Geodesic& g, double latitude, double longitude, double azimuth)
      : g_(g), line_(g.line(latitude, longitude, azimuth)), latitude_(latitude), longitude_(longitude), azimuth_(azimuth) {}

  // The point at distance s (m) along the line from the centre, negative behind it.
  SectionPoint at(double s_m) const {
    const auto d = line_.position(s_m);
    return {s_m, d.lat2, d.lon2, d.azi2};
  }

  // Samples every `step_m` from −half_m to +half_m; the number of steps each side is half / step
  // rounded, so the ends are ±half when half is a multiple of the step.
  Section sample(double half_m, double step_m) const {
    if (!(half_m > 0 && step_m > 0)) throw std::invalid_argument("geophysics::structure: a section needs a positive half-length and step");
    const int n = int(std::lround(half_m / step_m));
    if (n > 100000) throw std::invalid_argument("geophysics::structure: more than 200,000 samples along a section");
    Section s{latitude_, longitude_, azimuth_, half_m, step_m, {}};
    s.points.reserve(std::size_t(2 * n + 1));
    for (int k = -n; k <= n; k++) s.points.push_back(at(k == -n ? -half_m : k == n ? half_m : k * step_m));
    return s;
  }

  // Where a place (rad) lies relative to the line, searching for the foot within `search_m` of
  // `guess_m` along it.
  SectionOffset offset(double latitude, double longitude, double guess_m, double search_m) const {
    const auto distance = [&](double s) {
      const auto p = line_.position(s);
      return g_.inverse(p.lat2, p.lon2, latitude, longitude).s12;
    };
    // Brent finds the foot to about √ε |s| (some millimetres at 1000 km); two steps of the sphere's
    // right triangle, tan δs = tan(d / R) cos Δα, take it to rounding, since each leaves an error of
    // the ellipsoid's flattening times the step before.
    double s = calx::optimisation::brent_minimise(distance, guess_m - search_m, guess_m + search_m, 1e-9).x;
    const double R = g_.a() * (1 - g_.flattening() / 3);
    auto foot = line_.position(s);
    auto to = g_.inverse(foot.lat2, foot.lon2, latitude, longitude);
    for (int k = 0; k < 2 && to.s12 > 0; k++) {
      s += R * std::atan(std::tan(to.s12 / R) * std::cos(to.azi1 - foot.azi2));
      foot = line_.position(s);
      to = g_.inverse(foot.lat2, foot.lon2, latitude, longitude);
    }
    const double turn = std::sin(to.azi1 - foot.azi2);
    return {s, (turn >= 0 ? 1.0 : -1.0) * to.s12, foot.azi2};
  }

  // As offset(), with the guess from the sphere of the ellipsoid's mean radius (js/section.js's
  // relativeTo) and a search wide enough for the ellipsoid's departure from it.
  SectionOffset offset(double latitude, double longitude) const {
    const double R = g_.a() * (1 - g_.flattening() / 3);
    const double p1 = latitude_, p2 = latitude, dl = longitude - longitude_;
    const double h = std::pow(std::sin((p2 - p1) / 2), 2) + std::cos(p1) * std::cos(p2) * std::pow(std::sin(dl / 2), 2);
    const double delta = 2 * std::asin(std::min(1.0, std::sqrt(h)));
    const double theta = std::atan2(std::sin(dl) * std::cos(p2), std::cos(p1) * std::sin(p2) - std::sin(p1) * std::cos(p2) * std::cos(dl));
    const double along = std::atan2(std::sin(delta) * std::cos(theta - azimuth_), std::cos(delta)) * R;
    return offset(latitude, longitude, along, 0.01 * std::fabs(along) + 0.01 * delta * R + 2000.0);
  }

  double latitude() const { return latitude_; }
  double longitude() const { return longitude_; }
  double azimuth() const { return azimuth_; }

 private:
  const geodesy::geodesic::Geodesic& g_;
  geodesy::geodesic::GeodesicLine line_;
  double latitude_, longitude_, azimuth_;
};

}  // namespace geophysics::structure
