// The toolbox's space: tyche's spatial autocorrelation over places read as longitudes and latitudes. A tool's rows are
// places when the page gives their coordinates; neighbours are found on the sphere, each place's k nearest. Residuals are
// tested in the same way, since standard errors assume what neighbouring places commonly break.

#pragma once

#include <tyche/spatial/autocorrelation.hpp>

#include <numbers>
#include <set>
#include <utility>

#include "answer.hpp"

namespace science::toolbox {

namespace tspatial = tyche::spatial;

// Places as unit vectors: chords are monotone in the angle, so that neighbours are the sphere's, across the antimeridian too.
struct Sphere {
  std::vector<double> x, y, z;
};
inline Sphere sphere_of(const std::vector<double>& lon, const std::vector<double>& lat) {
  constexpr double deg = std::numbers::pi / 180.0;
  Sphere S;
  for (std::size_t i = 0; i < lon.size(); i++) {
    S.x.push_back(std::cos(lat[i] * deg) * std::cos(lon[i] * deg));
    S.y.push_back(std::cos(lat[i] * deg) * std::sin(lon[i] * deg));
    S.z.push_back(std::sin(lat[i] * deg));
  }
  return S;
}

inline const char* quadrant_name(int q) { return q == 1 ? "high among high" : q == 2 ? "low among high" : q == 3 ? "low among low" : q == 4 ? "high among low" : "no neighbours"; }

inline Json moran(const Request& r) {
  const auto values = r.numbers(r.one("value")), lon = r.numbers(r.one("lon")), lat = r.numbers(r.one("lat"));
  const auto C = complete({values, lon, lat});
  const std::size_t n = C.rows.size();
  if (n < 8) throw std::invalid_argument("too few places with a value and a position: at least eight");
  std::set<std::pair<double, double>> seen;
  for (std::size_t i = 0; i < n; i++) seen.insert({C.columns[1][i], C.columns[2][i]});
  if (seen.size() != n) throw std::invalid_argument("a place appears more than once: keep one row for each place (one year, say)");
  const auto k = static_cast<std::size_t>(std::clamp(r.number("k", 6), 1.0, static_cast<double>(n - 1)));
  const auto permutations = static_cast<std::size_t>(std::clamp(r.number("permutations", 499), 0.0, 9999.0));
  const auto S = sphere_of(C.columns[1], C.columns[2]);
  const auto W = tspatial::k_nearest_weights(S.x, S.y, S.z, k);
  const auto M = tspatial::morans_i(C.columns[0], W, permutations, 1);
  const auto L = tspatial::local_morans_i(C.columns[0], W, permutations, 2);

  Arr blocks;
  blocks.push_back(figures("Moran's I of " + r.one("value"), {{"I", M.I}, {"its expectation with no association", M.expected}, {"z, under random assignment", M.z_randomisation},
                                                             {"p", M.p_randomisation}, {"p, by permutation", M.p_permutation}, {"places", static_cast<double>(n)}, {"neighbours each", static_cast<double>(k)}}));
  // Moran's scatter: a place against the mean of its neighbours, whose slope on row-standardised weights is I.
  double mean = 0, var = 0;
  for (double v : C.columns[0]) mean += v / static_cast<double>(n);
  for (double v : C.columns[0]) var += (v - mean) * (v - mean) / static_cast<double>(n);
  const double sd = std::sqrt(var);
  std::vector<double> px(n), py(n);
  for (std::size_t i = 0; i < n; i++) {
    double lag = 0;
    for (std::size_t e = W.row_start[i]; e < W.row_start[i + 1]; e++) lag += W.weight[e] * (C.columns[0][W.column[e]] - mean);
    px[i] = (C.columns[0][i] - mean) / sd;
    py[i] = lag / sd;
  }
  const auto [lo, hi] = std::minmax_element(px.begin(), px.end());
  blocks.push_back(scatter("Moran's scatter: each place against its neighbours, in standard deviations; the line's slope is I", px, py, "the place", "its neighbours", {*lo, *hi}, {M.I * *lo, M.I * *hi}));

  // The clusters: +1 a hot spot, −1 a cold spot, ±½ an outlier, 0 what is not told apart from chance at 5%, uncorrected.
  std::vector<double> local(n), cluster(n);
  std::map<std::string, double> counts;
  for (std::size_t i = 0; i < n; i++) {
    local[i] = L.I[i];
    const bool told = permutations == 0 || L.p[i] <= 0.05;
    const int q = told ? L.quadrant[i] : 0;
    cluster[i] = q == 1 ? 1.0 : q == 3 ? -1.0 : q == 4 ? 0.5 : q == 2 ? -0.5 : 0.0;
    counts[told ? quadrant_name(L.quadrant[i]) : "not told apart from chance"] += 1;
  }
  std::vector<std::pair<std::string, double>> items(counts.begin(), counts.end());
  blocks.push_back(figures("The places, at 5% by conditional permutation, uncorrected for testing each", items));
  blocks.push_back(rows_of("Clusters: +1 high among high, −1 low among low, ±½ outliers, 0 not told apart from chance", C.rows, cluster, "diverging"));
  blocks.push_back(rows_of("Local Moran's I", C.rows, local, "diverging"));
  if (!r.role("id").empty()) {
    const auto ids = r.texts(r.one("id"));
    std::map<std::string, double> by;
    for (std::size_t i = 0; i < n; i++) by[ids[C.rows[i]]] = cluster[i];
    blocks.push_back(map_values("Clusters of " + r.one("value") + ": +1 high among high, −1 low among low, ±½ outliers, 0 not told apart from chance", by, "diverging"));
  }
  blocks.push_back(note("Moran (1950); Cliff and Ord (1981); local after Anselin (1995). Each place's " + std::to_string(k) +
                        " nearest neighbours on the sphere. Where a common trend runs across the map, I measures the trend: a regression's residuals are the thing to test. tyche::spatial, held to PySAL's esda."));
  return answer("Spatial autocorrelation", std::move(blocks));
}

// A regression's residuals, tested where the rows are distinct places: blocks to add to its answer, and its map's.
// The page names the rows' places in "places": {"lon": [...], "lat": [...]}, along the whole table.
inline void residuals_tested(const Request& r, const std::vector<std::size_t>& rows, const std::vector<double>& residuals, Arr& blocks) {
  const Json* places = r.json.find("places");
  if (!places || !places->is_object()) return;
  const Json* jlon = places->find("lon");
  const Json* jlat = places->find("lat");
  if (!jlon || !jlat || !jlon->is_array() || !jlat->is_array()) return;
  std::vector<double> lon, lat, e;
  const auto& alon = jlon->as_array();
  const auto& alat = jlat->as_array();
  std::set<std::pair<double, double>> seen;
  for (std::size_t m = 0; m < rows.size(); m++) {
    const std::size_t i = rows[m];
    if (i >= alon.size() || i >= alat.size() || !alon[i].is_number() || !alat[i].is_number() || !std::isfinite(residuals[m])) continue;
    if (!seen.insert({alon[i].as_number(), alat[i].as_number()}).second) {
      blocks.push_back(note("A place appears more than once here, so that its residuals are not tested across space: keep one row for each place to test them."));
      return;
    }
    lon.push_back(alon[i].as_number()), lat.push_back(alat[i].as_number()), e.push_back(residuals[m]);
  }
  if (e.size() < 8) return;
  const std::size_t k = std::min<std::size_t>(6, e.size() - 1);
  const auto S = sphere_of(lon, lat);
  double spread = 0;
  for (double v : e) spread += v * v;
  if (!(spread > 0)) return;
  const auto M = tspatial::morans_i(e, tspatial::k_nearest_weights(S.x, S.y, S.z, k), 499, 3);
  blocks.push_back(figures("Do neighbours' residuals resemble each other?", {{"Moran's I of the residuals", M.I}, {"p, by permutation", M.p_permutation}, {"neighbours each", static_cast<double>(k)}}));
  blocks.push_back(note(M.p_permutation <= 0.05 && M.I > 0
                            ? "They do. The standard errors above treat places as independent, so that they are too small here, and something shared by neighbours is missing from the model."
                            : "They are not told apart from chance here (p above 5% or I below zero), so that this much of the standard errors' assumption stands."));
}

// What the page asks for by row: a regression's residuals and fitted values, in the rows it kept.
inline bool wants_rows(const Request& r) {
  const Json* b = r.json.find("by_row");
  return b && b->is_bool() && b->as_bool();
}
inline void residual_blocks(const Request& r, const std::vector<std::size_t>& rows, const std::vector<double>& residuals, const std::vector<double>& fitted, const std::string& outcome, Arr& blocks) {
  if (!wants_rows(r)) return;
  residuals_tested(r, rows, residuals, blocks);
  blocks.push_back(rows_of("Residuals: what " + outcome + " holds beyond the fit", rows, residuals, "diverging"));
  blocks.push_back(rows_of("Fitted values of " + outcome, rows, fitted, "sequential"));
}

}  // namespace science::toolbox
