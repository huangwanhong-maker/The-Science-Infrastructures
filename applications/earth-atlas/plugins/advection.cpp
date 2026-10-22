// A temperature carried by the wind, in C++: a ForecastModel (api/cpp/plugin.hpp)
// that the server compiles and runs. The field chosen (the 2 m temperature by
// default, or the temperature at any level the library holds) is advected by the
// wind chosen (the 10 m wind by default), the wind taken from the source's own
// times as they pass, by a semi-Lagrangian scheme on the sphere:
//
//   φ(x_a, t + Δt) = φ(x_d, t),
//
// x_a a grid point, x_d the point the parcel arriving there left a step before
// (Staniforth and Côté 1991; Robert 1981). The departure point is found in
// three-dimensional Cartesian coordinates, so that nothing is singular at the
// poles, by the midpoint rule of Ritchie (1987): with r_a the unit vector of x_a,
// V the wind as a Cartesian vector tangent to the sphere and a the Earth's radius,
//
//   r_m ← (r_a − (Δt / 2a) V(r_m, t + Δt/2)) / |·|,   iterated three times from r_m = r_a;
//   r_d = 2 (r_a · r_m) r_m − r_a,                     the point beyond r_m on the great circle from r_a.
//
// φ(x_d) is interpolated by cubic Lagrange polynomials in longitude and then in
// latitude (16 points), the rows beyond a pole taken from across it, 180° away;
// "clipped" bounds the result by the four nearest values, which makes the scheme
// quasi-monotone (Bermejo and Staniforth 1992); "linear" is bilinear. Nothing is
// added or lost but by interpolation: no sources, no diffusion, no diurnal cycle.
// The base class (in the host) verifies the forecast against the source's own
// later times beside persistence.
//
// Check (applications/earth-atlas/py/check.py): the cosine bell of Williamson et
// al. (1992), test case 1, carried by a solid-body rotation over the poles
// (α = π/2 − 0.05) for one revolution of 12 days, returns to itself; and a uniform
// field stays uniform.
//
// Caveat. A plain lat-lon grid is finer in longitude near the poles than at the
// equator; the scheme is stable at any step (it is semi-Lagrangian), but the
// interpolation there is no more accurate than the rows allow.

#include "plugin.hpp"

using plugin::Grid;

struct Tracer {
  Grid grid;
  std::vector<double> phi, next;
  std::vector<std::array<double, 3>> arrival;     // the unit vectors of the grid points
  const plugin::Input* wind = nullptr;
  double time = 0;                                // ms
  std::string interpolation;
};

// The cubic Lagrange weights at nodes −1, 0, 1, 2 for a point at f in [0, 1).
static void cubic_weights(double f, double w[4]) {
  w[0] = -f * (f - 1) * (f - 2) / 6;
  w[1] = (f + 1) * (f - 1) * (f - 2) / 2;
  w[2] = -(f + 1) * f * (f - 2) / 2;
  w[3] = (f + 1) * f * (f - 1) / 6;
}

// A row beyond a pole, reflected across it: the row and whether to look 180° away;
// a row that has no reflection on the grid is clamped.
static int resolve_row(const Grid& g, int row, bool& across) {
  across = false;
  if (row >= 0 && row < g.nlat) return row;
  if (g.periodic) {
    const double m = (row < 0 ? 2 * g.lat0 - 180 : 2 * g.lat0 + 180) / g.dlat;      // row' = m − row
    const long k = std::lround(m);
    if (std::abs(m - k) < 1e-6 && k - row >= 0 && k - row < g.nlat) { across = true; return int(k - row); }
  }
  return std::clamp(row, 0, g.nlat - 1);
}

// The field at (lon, lat) by cubic interpolation; clipped to the four nearest values when asked.
static double cubic(const std::vector<double>& phi, const Grid& g, double lon, double lat, bool clip, const std::vector<float>& fallback) {
  double x = (lon - g.lon0) / g.dlon;
  const double y = (g.lat0 - lat) / g.dlat;
  if (!g.periodic) return plugin::bilinear(fallback, g, lon, lat);
  x = std::fmod(std::fmod(x, g.nlon) + g.nlon, g.nlon);
  const int j = int(std::floor(y));
  double wy[4];
  cubic_weights(y - j, wy);
  double value = 0, lo = INFINITY, hi = -INFINITY;
  for (int q = 0; q < 4; q++) {
    bool across;
    const int row = resolve_row(g, j - 1 + q, across);
    double xr = across ? x + 180 / g.dlon : x;
    xr = std::fmod(xr, g.nlon);
    const int i = int(std::floor(xr));
    double wx[4];
    cubic_weights(xr - i, wx);
    double s = 0;
    for (int p = 0; p < 4; p++) {
      const double v = phi[size_t(row) * g.nlon + ((i - 1 + p) % g.nlon + g.nlon) % g.nlon];
      if (std::isnan(v)) return plugin::bilinear(fallback, g, lon, lat);
      s += wx[p] * v;
      if ((q == 1 || q == 2) && (p == 1 || p == 2)) { lo = std::min(lo, v); hi = std::max(hi, v); }
    }
    value += wy[q] * s;
  }
  return clip ? std::clamp(value, lo, hi) : value;
}

static std::array<double, 3> unit(double lon, double lat) {
  const double l = lon * plugin::DEG, p = lat * plugin::DEG;
  return {std::cos(p) * std::cos(l), std::cos(p) * std::sin(l), std::sin(p)};
}

struct Advection : plugin::ForecastModel<Tracer> {
  Grid grid_;

  Advection() {
    id = "advection-cpp";
    title = "Advection (C++)";
    description = "The temperature chosen, carried by the wind chosen through the source's own times, by a semi-Lagrangian scheme on the "
                  "sphere: no sources, no mixing, only the motion. Verified against the source's later times beside persistence.";
    dt = 3600;
    every = 6;
    needs = {{"scalar", "Carried", {"air_temperature"}, "", 2}, {"wind", "By the wind", {"eastward_wind", "northward_wind"}, "", 2}};
    parameters = {{"hours", "Forecast length", 48, 6, 360, 6, "h"},
                  plugin::choice("interpolation", "Interpolation", "clipped",
                                 {{"clipped", "cubic, clipped (quasi-monotone)"}, {"cubic", "cubic"}, {"linear", "bilinear"}})};
    outputs = {{"t", "Temperature carried by the wind", "air_temperature", "K", "", false}};
    show = {{"surface", "t"}};
  }

  Grid grid() const override { return grid_; }
  plugin::Verify verify() const override { return {"t", "scalar", "temperature", "K"}; }

  Tracer initial(const plugin::Inputs& inputs, int index, const plugin::Params& params, const plugin::Context&) override {
    const plugin::Input& scalar = inputs.at("scalar");
    Tracer s;
    s.grid = grid_ = scalar.grid;
    s.wind = &inputs.at("wind");
    s.time = scalar.times_ms.empty() ? 0 : scalar.times_ms[index];
    s.interpolation = params.choice("interpolation", "clipped");
    const float* start = scalar.slice(0, index);
    s.phi.assign(start, start + s.grid.size());
    s.next.resize(s.phi.size());
    for (int r = 0; r < s.grid.nlat; r++)
      for (int c = 0; c < s.grid.nlon; c++) s.arrival.push_back(unit(s.grid.lon(c), s.grid.lat(r)));
    // The output is the quantity carried, at its level and in its units.
    outputs[0].level = scalar.level;
    if (!scalar.variables.empty() && !scalar.variables[0].units.empty()) outputs[0].units = scalar.variables[0].units;
    outputs[0].long_name = (scalar.level.empty() ? std::string("Temperature") : "Temperature at " + scalar.level) + ", carried by the wind";
    return s;
  }

  void advance(Tracer& s, double seconds, const plugin::Params&) override {
    const Grid& g = s.grid;
    const double tm = s.time + 500 * seconds, half = seconds / 2 / plugin::EARTH_RADIUS_M;
    const bool clip = s.interpolation == "clipped", linear = s.interpolation == "linear";
    std::vector<float> fallback(s.phi.begin(), s.phi.end());
    for (size_t i = 0; i < s.arrival.size(); i++) {
      const auto& a = s.arrival[i];
      std::array<double, 3> m = a;
      for (int iteration = 0; iteration < 3; iteration++) {
        const double lon = std::atan2(m[1], m[0]) / plugin::DEG, lat = std::asin(std::clamp(m[2], -1.0, 1.0)) / plugin::DEG;
        double u = s.wind->sample(0, tm, lon, lat), v = s.wind->sample(1, tm, lon, lat);
        if (std::isnan(u) || std::isnan(v)) u = v = 0;      // no wind there: the parcel stays
        // The wind as a Cartesian vector: u along the eastward unit vector, v along the northward one.
        const double l = lon * plugin::DEG, p = lat * plugin::DEG, sl = std::sin(l), cl = std::cos(l), sp = std::sin(p), cp = std::cos(p);
        const double V[3] = {-u * sl - v * sp * cl, u * cl - v * sp * sl, v * cp};
        double n = 0;
        for (int k = 0; k < 3; k++) { m[k] = a[k] - half * V[k]; n += m[k] * m[k]; }
        n = std::sqrt(n);
        for (int k = 0; k < 3; k++) m[k] /= n;
      }
      const double am = a[0] * m[0] + a[1] * m[1] + a[2] * m[2];
      const double d[3] = {2 * am * m[0] - a[0], 2 * am * m[1] - a[1], 2 * am * m[2] - a[2]};
      const double lon = std::atan2(d[1], d[0]) / plugin::DEG, lat = std::asin(std::clamp(d[2], -1.0, 1.0)) / plugin::DEG;
      s.next[i] = linear ? plugin::bilinear(fallback, g, lon, lat) : cubic(s.phi, g, lon, lat, clip, fallback);
    }
    std::swap(s.phi, s.next);
    s.time += seconds * 1000;
  }

  std::map<std::string, std::vector<float>> record(const Tracer& s, const plugin::Params&) override {
    return {{"t", std::vector<float>(s.phi.begin(), s.phi.end())}};
  }
};

PLUGIN_MAIN(Advection)
