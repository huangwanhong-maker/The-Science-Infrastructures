// The 2 m temperature analysed toward station reports by successive corrections
// (Cressman 1959), in C++: a DataAssimilation plug-in (api/cpp/plugin.hpp) that
// the server compiles and runs. The background is the model chosen, at the
// instant shown (the base class interpolates it in time); the observations are
// each station's report nearest that instant within a window, from the
// observation set chosen: the METAR reports of aerodromes, the buoys, or our own
// sensors. With innovations d_k = y_k − H x at the observations (H bilinear),
//
//   x(i) ← x(i) + Σ_k w_ik d_k / (Σ_k w_ik + ε²),
//   w_ik = (R² − r_ik²) / (R² + r_ik²) for r_ik < R, and 0 beyond,
//
// r_ik the great-circle distance from grid point i to observation k. Each pass
// takes the analysis of the one before as its background, its radius half of the
// last (Cressman 1959; Kalnay 2003, §5.3). With ε² = 0, Cressman's own form, a
// lone observation is drawn exactly: the analysis takes its value wherever it is
// the only one within R. Reports whose first innovation exceeds a limit are
// rejected as gross errors.
//
// Check (applications/earth-atlas/py/check.py): with one observation on a grid
// point, the analysis there equals the observation, and the background is kept
// beyond the radius.
//
// Caveat. No correction is made for the height of a station above the model's
// ground, and the scheme has no notion of the errors of either; optimal
// interpolation (optimal-interpolation.js) weighs them.

#include "plugin.hpp"

using plugin::DEG;
using plugin::EARTH_RADIUS_M;

struct Cressman : plugin::DataAssimilation {
  // What the last analysis found, for its table.
  int used = 0, rejected = 0;
  double omb = 0, oma = 0;
  std::string source;

  Cressman() {
    id = "cressman-cpp";
    title = "Successive correction (C++)";
    description = "The 2 m temperature analysed toward weather stations, buoys or your own sensors by Cressman's successive corrections, "
                  "at the instant shown, in a C++ program in the server. The increment is shown, and the analysis joins the models of the quantity.";
    needs = {{"background", "Background", {"air_temperature"}, "2 m", 0}, plugin::observed("stations", "Observed by", "air_temperature")};
    parameters = {{"radius", "First radius R", 800, 50, 3000, 50, "km"},
                  {"passes", "Passes (R halved each)", 3, 1, 6, 1, ""},
                  {"epsilon", "ε² (0: Cressman's)", 0, 0, 10, 0.05, ""},
                  {"window", "Reports within", 1.5, 0.25, 12, 0.25, "h"},
                  {"reject", "Reject innovations beyond", 10, 0.5, 50, 0.5, "K"}};
  }

  // Each station's report nearest t within the window, with a value.
  std::vector<plugin::Observation> observations(const plugin::Inputs& inputs, double t, const plugin::Params& params, const plugin::Context&) override {
    const plugin::Input& stations = inputs.at("stations");
    const plugin::Observations& set = stations.observations;
    const int id_ = set.column("id"), lon = set.column("lon"), lat = set.column("lat"), time = set.column("time"), value = set.column(stations.variable);
    if (lon < 0 || lat < 0 || value < 0) throw std::runtime_error(stations.product + " has no lon, lat or " + stations.variable);
    const double window = params.number("window", 1.5) * 3600000;
    std::map<std::string, std::pair<double, plugin::Observation>> best;
    for (size_t r = 0; r < set.rows.size(); r++) {
      const double y = set.number(r, value), when = time < 0 ? t : set.time(r, time);
      if (!std::isfinite(y) || !std::isfinite(when) || std::abs(when - t) > window) continue;
      const std::string key = id_ < 0 ? std::to_string(r) : set.text(r, id_);
      auto it = best.find(key);
      if (it == best.end() || std::abs(when - t) < it->second.first)
        best[key] = {std::abs(when - t), {set.number(r, lon), set.number(r, lat), y, plugin::NaN, key}};
    }
    // Places on the sphere only: some reports carry a placeholder, such as (−99.99, −99.99).
    std::vector<plugin::Observation> out;
    for (auto& [key, entry] : best) {
      const plugin::Observation& o = entry.second;
      if (std::abs(o.lat) <= 90 && o.lon >= -180 && o.lon <= 360) out.push_back(o);
    }
    char text[200];
    std::snprintf(text, sizeof text, "%zu reports of %s within %g h", out.size(), stations.model.c_str(), params.number("window", 1.5));
    source = text;
    return out;
  }

  std::vector<float> analyse(const std::vector<float>& background, const plugin::Grid& g, const std::vector<plugin::Observation>& all,
                             const plugin::Params& params, const plugin::Context& context) override {
    const double R0 = params.number("radius", 800) * 1000, eps = params.number("epsilon", 0), limit = params.number("reject", 10);
    const int passes = std::max(1, int(std::lround(params.number("passes", 3))));
    // The observations kept, by latitude, with their unit vectors.
    struct Obs { double lat, value, x, y, z; plugin::LonLat at; };
    std::vector<Obs> obs;
    rejected = 0;
    double sum_b = 0;
    for (const auto& o : all) {
      const double d = o.value - plugin::bilinear(background, g, o.lon, o.lat);
      if (std::isnan(d)) continue;
      if (std::abs(d) > limit) { rejected++; continue; }
      sum_b += d * d;
      const double l = o.lon * DEG, p = o.lat * DEG;
      obs.push_back({o.lat, o.value, std::cos(p) * std::cos(l), std::cos(p) * std::sin(l), std::sin(p), {o.lon, o.lat}});
    }
    std::sort(obs.begin(), obs.end(), [](const Obs& a, const Obs& b) { return a.lat < b.lat; });
    used = int(obs.size());
    omb = used ? std::sqrt(sum_b / used) : 0;
    std::vector<float> x = background;
    std::vector<double> d(obs.size());
    for (int pass = 0; pass < passes; pass++) {
      const double R = R0 / std::pow(2.0, pass), reach = std::cos(std::min(plugin::PI, R / EARTH_RADIUS_M)), band = R / EARTH_RADIUS_M / DEG;
      for (size_t k = 0; k < obs.size(); k++) d[k] = obs[k].value - plugin::bilinear(x, g, obs[k].at.lon, obs[k].at.lat);
      std::vector<float> next = x;
      for (int r = 0; r < g.nlat; r++) {
        const double lat = g.lat(r), p = lat * DEG;
        // The observations within the band of latitudes R can reach.
        const auto lo = std::lower_bound(obs.begin(), obs.end(), lat - band, [](const Obs& o, double v) { return o.lat < v; });
        const auto hi = std::upper_bound(obs.begin(), obs.end(), lat + band, [](double v, const Obs& o) { return v < o.lat; });
        if (lo == hi) continue;
        for (int c = 0; c < g.nlon; c++) {
          const size_t i = size_t(r) * g.nlon + c;
          if (std::isnan(x[i])) continue;
          const double l = g.lon(c) * DEG, px = std::cos(p) * std::cos(l), py = std::cos(p) * std::sin(l), pz = std::sin(p);
          double sw = 0, swd = 0;
          for (auto o = lo; o != hi; ++o) {
            const double cs = px * o->x + py * o->y + pz * o->z;
            if (cs <= reach) continue;
            const double rr = EARTH_RADIUS_M * std::acos(std::min(1.0, cs)), w = (R * R - rr * rr) / (R * R + rr * rr);
            sw += w;
            swd += w * d[size_t(o - obs.begin())];
          }
          if (sw > 0) next[i] = float(x[i] + swd / (sw + eps));
        }
      }
      x = std::move(next);
      context.progress(double(pass + 1) / passes, "pass " + std::to_string(pass + 1) + " of " + std::to_string(passes) + ", R = " +
                       std::to_string(long(std::lround(R / 1000))) + " km");
    }
    double sum_a = 0;
    int n = 0;
    for (const auto& o : obs) {
      const double e = o.value - plugin::bilinear(x, g, o.at.lon, o.at.lat);
      if (!std::isnan(e)) { sum_a += e * e; n++; }
    }
    oma = n ? std::sqrt(sum_a / n) : 0;
    return x;
  }

  plugin::Results run(const plugin::Inputs& inputs, const plugin::Params& params, const plugin::Context& context) override {
    plugin::Results r = DataAssimilation::run(inputs, params, context);
    const plugin::Input& background = inputs.at("background");
    const std::string unit = background.variables.empty() ? "" : background.variables[0].units;
    const auto value = [&](double v) { char b[40]; std::snprintf(b, sizeof b, "%.2f %s", v, unit.c_str()); return std::string(b); };
    r.tables.push_back({"Successive correction of " + background.model + ": " + source, {"", "Value"},
                        {{"Observations used", std::to_string(used)},
                         {"Rejected (gross error)", std::to_string(rejected)},
                         {"RMS observation − background", value(omb)},
                         {"RMS observation − analysis", value(oma)}},
                        "The analysis draws to the reports within each radius, and keeps the background beyond them."});
    return r;
  }
};

PLUGIN_MAIN(Cressman)
