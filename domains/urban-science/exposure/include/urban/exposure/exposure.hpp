// What a hazard on the ground reaches: buildings at their places with their values, against a map of a hazard's intensity (a flood's
// depth, say) on a regular grid. The idea of an exposure analysis such as RAPIDA's or a catastrophe model's first stage; written
// here from that statement. It names no vulnerability: how much a depth damages is a further, declared curve (`damage`).
//
// The grid is on longitude and latitude, degrees, rows from the south, one byte a cell with a factor to the hazard's units (a depth code,
// as the data library's rasters hold it). A place outside the grid has no reading, and is counted apart rather than as dry.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace urban::exposure {

struct Grid {
  double south = 0, west = 0, north = 0, east = 0;
  std::size_t width = 0, height = 0;
  double unit = 1;                     // the hazard's units per byte
  std::vector<std::uint8_t> values;    // width * height, rows from the south

  void check() const {
    if (width == 0 || height == 0 || values.size() != width * height) throw std::invalid_argument("exposure: a grid has width x height values");
    if (!(north > south) || !(east > west)) throw std::invalid_argument("exposure: a grid's box is not empty");
  }
  bool covers(double lat, double lon) const { return lat >= south && lat < north && lon >= west && lon < east; }
  // The intensity of the cell holding a place, or a negative number outside the grid.
  double at(double lat, double lon) const {
    if (!covers(lat, lon)) return -1;
    const auto i = std::min(width - 1, static_cast<std::size_t>((lon - west) / (east - west) * static_cast<double>(width)));
    const auto j = std::min(height - 1, static_cast<std::size_t>((lat - south) / (north - south) * static_cast<double>(height)));
    return values[j * width + i] * unit;
  }
};

struct Asset {
  double lat = 0, lon = 0, value = 0;
  int zone = -1;                       // the unit it falls in, or -1
};

struct Tally {
  std::size_t count = 0;
  double value = 0, damage = 0;
};

struct Exposed {
  Tally reached;                       // those at or above the threshold
  Tally all;                           // those with a reading
  std::size_t outside = 0;             // those beyond the grid
  std::vector<Tally> zones;            // reached, by zone
  std::vector<double> bands;           // the upper edges of the bands of intensity
  std::vector<Tally> by_band;          // reached, by band: band k holds intensity in [bands[k-1], bands[k]), the last open above
};

// The damage ratio at an intensity from a curve of points, linearly between them, and constant beyond its ends.
inline double interpolate(const std::vector<double>& x, const std::vector<double>& y, double v) {
  if (x.size() != y.size() || x.empty()) throw std::invalid_argument("exposure: a curve has a value at each point");
  if (v <= x.front()) return y.front();
  if (v >= x.back()) return y.back();
  const auto hi = static_cast<std::size_t>(std::upper_bound(x.begin(), x.end(), v) - x.begin());
  const double t = (v - x[hi - 1]) / (x[hi] - x[hi - 1]);
  return y[hi - 1] + t * (y[hi] - y[hi - 1]);
}

inline Exposed exposed(const Grid& grid, const std::vector<Asset>& assets, double threshold, const std::vector<double>& bands = {}, const std::vector<double>& curve_x = {}, const std::vector<double>& curve_y = {}, std::size_t zones = 0) {
  grid.check();
  Exposed out;
  out.bands = bands;
  out.by_band.assign(bands.size() + 1, Tally{});
  out.zones.assign(zones, Tally{});
  const bool damaged = !curve_x.empty();
  for (const Asset& a : assets) {
    const double v = grid.at(a.lat, a.lon);
    if (v < 0) { out.outside++; continue; }
    out.all.count++;
    out.all.value += a.value;
    if (v < threshold || v <= 0) continue;
    const double d = damaged ? a.value * interpolate(curve_x, curve_y, v) : 0;
    const auto add = [&](Tally& t) { t.count++; t.value += a.value; t.damage += d; };
    add(out.reached);
    add(out.by_band[static_cast<std::size_t>(std::upper_bound(bands.begin(), bands.end(), v) - bands.begin())]);
    if (a.zone >= 0 && static_cast<std::size_t>(a.zone) < zones) add(out.zones[static_cast<std::size_t>(a.zone)]);
  }
  return out;
}

}  // namespace urban::exposure
