// The gravity of trade: flows between economies against their sizes and the distance between them
// (Tinbergen 1962), fitted on every ordered pair, zeros included.
//
// THE MODEL. With x_ij the exports of i to j and d_ij the great-circle distance between them,
//   E(x_ij) = exp(β ln d_ij + α_i + γ_j)                       (exporter and importer effects)
//   E(x_ij) = exp(β₀ + β₁ ln Y_i + β₂ ln Y_j + β ln d_ij)      (sizes)
// The effects absorb every property of the exporter and of the importer, among them their sizes
// and their multilateral resistance, the average remoteness that a bilateral distance is measured
// against (Anderson and van Wincoop 2003; Feenstra 2016, ch. 5); the model with sizes omits that
// resistance and its distance elasticity is biased by it. Both are fitted by Poisson
// pseudo-maximum likelihood, which is consistent for the conditional mean whatever the
// distribution of x, keeps the pairs that do not trade, and is not biased by the
// heteroskedasticity that biases least squares on logarithms (Santos Silva and Tenreyro 2006).
// The effects model uses nomos's Poisson with two absorbed effects, the sizes model tyche's GLM;
// least squares on the logarithms of the positive flows is given beside them for comparison.
//
// Distances are great-circle angles (physis) times the Earth's mean radius R₁ = 6371.0088 km
// (IUGG; Moritz 2000), between places given in radians; the elasticity β does not depend on the
// unit. A pair of places at one point has no logarithm of distance and is left out.

#pragma once

#include <economics/trade/partners.hpp>
#include <nomos/econometrics/poisson_fixed_effects.hpp>
#include <physis/computational_geometry/sphere.hpp>
#include <tyche/regression/glm.hpp>
#include <tyche/regression/linear.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace economics::trade {

inline constexpr double earth_mean_radius_km = 6371.0088;

// The great-circle distance in km between places at (latitude, longitude) in radians.
inline double great_circle_km(double lat1, double lon1, double lat2, double lon2) {
  namespace cg = physis::computational_geometry;
  return earth_mean_radius_km * cg::angle(cg::unit_vector(lat1, lon1), cg::unit_vector(lat2, lon2));
}

enum class GravityMethod { ppml_effects, ppml_sizes, ols_sizes };

struct Gravity {
  GravityMethod method = GravityMethod::ppml_effects;
  std::vector<std::string> names;                       // the coefficients' names
  std::vector<double> beta, se, p, ci_lower, ci_upper;  // in that order
  double distance_elasticity = std::numeric_limits<double>::quiet_NaN();
  std::size_t pairs = 0, zeros = 0;                     // the pairs fitted, and those without trade
  std::vector<std::size_t> exporter, importer;          // each pair fitted
  std::vector<double> flow, fitted, distance_km;        // its flow, the model's expectation, its distance
  std::vector<double> exporter_effect, importer_effect; // α̂_i, γ̂_j (effects model; NaN otherwise)
  double pseudo_r_squared = std::numeric_limits<double>::quiet_NaN();  // the squared correlation of flows and fitted values
  int iterations = 0;
  bool converged = true;
};

namespace detail {
struct Pairs {
  std::vector<std::size_t> i, j;
  std::vector<double> x, log_d, d;
};
inline Pairs pairs(const Rectangular& x, const std::vector<double>& latitude, const std::vector<double>& longitude) {
  require_square(x);
  if (latitude.size() != x.rows || longitude.size() != x.rows) throw std::invalid_argument("economics::trade::gravity: a position for every economy");
  Pairs P;
  for (std::size_t i = 0; i < x.rows; i++)
    for (std::size_t j = 0; j < x.cols; j++) {
      if (i == j) continue;
      const double d = great_circle_km(latitude[i], longitude[i], latitude[j], longitude[j]);
      if (!(d > 0)) continue;
      P.i.push_back(i);
      P.j.push_back(j);
      P.x.push_back(x(i, j));
      P.d.push_back(d);
      P.log_d.push_back(std::log(d));
    }
  return P;
}
inline double squared_correlation(const std::vector<double>& a, const std::vector<double>& b) {
  const double n = static_cast<double>(a.size());
  double ma = 0, mb = 0;
  for (std::size_t k = 0; k < a.size(); k++) ma += a[k] / n, mb += b[k] / n;
  double sab = 0, saa = 0, sbb = 0;
  for (std::size_t k = 0; k < a.size(); k++) {
    sab += (a[k] - ma) * (b[k] - mb);
    saa += (a[k] - ma) * (a[k] - ma);
    sbb += (b[k] - mb) * (b[k] - mb);
  }
  return saa > 0 && sbb > 0 ? sab * sab / (saa * sbb) : std::numeric_limits<double>::quiet_NaN();
}
}  // namespace detail

// PPML with exporter and importer effects; standard errors clustered by pair (each unordered pair
// one cluster, as its two directions share what is not modelled) unless `by_cell`.
inline Gravity gravity_effects(const Rectangular& x, const std::vector<double>& latitude, const std::vector<double>& longitude, bool by_cell = false) {
  const auto P = detail::pairs(x, latitude, longitude);
  Rectangular X(P.x.size(), 1);
  for (std::size_t k = 0; k < P.x.size(); k++) X(k, 0) = P.log_d[k];
  nomos::econometrics::PoissonTwoWayOptions o;
  if (!by_cell) {
    o.covariance = nomos::econometrics::PoissonCovariance::cluster;
    for (std::size_t k = 0; k < P.x.size(); k++) {
      const auto a = std::min(P.i[k], P.j[k]), b = std::max(P.i[k], P.j[k]);
      o.clusters.push_back(static_cast<std::int64_t>(a * x.rows + b));
    }
  }
  const auto f = nomos::econometrics::poisson_two_way(X, P.x, P.i, P.j, o);
  Gravity g;
  g.method = GravityMethod::ppml_effects;
  g.names = {"ln distance"};
  g.beta = f.beta;
  g.se = f.se;
  g.p = f.p;
  g.ci_lower = f.ci_lower;
  g.ci_upper = f.ci_upper;
  g.distance_elasticity = f.beta[0];
  g.exporter = P.i;
  g.importer = P.j;
  g.flow = P.x;
  g.fitted = f.fitted;
  g.distance_km = P.d;
  g.exporter_effect = f.row_effects;
  g.importer_effect = f.column_effects;
  g.pairs = P.x.size();
  for (double v : P.x) g.zeros += v == 0;
  g.iterations = f.iterations;
  g.converged = f.converged;
  g.pseudo_r_squared = detail::squared_correlation(P.x, f.fitted);
  return g;
}

// With sizes: by PPML (tyche's Poisson GLM, robust errors) or by least squares on the logarithms
// of the positive flows (robust HC1). Pairs where a size is not positive are left out.
inline Gravity gravity_sizes(const Rectangular& x, const std::vector<double>& latitude, const std::vector<double>& longitude, const std::vector<double>& size,
                             GravityMethod method = GravityMethod::ppml_sizes) {
  if (size.size() != x.rows) throw std::invalid_argument("economics::trade::gravity: a size for every economy");
  const auto P = detail::pairs(x, latitude, longitude);
  std::vector<double> lyi, lyj, ld, y;
  std::vector<std::size_t> ei, ej, kept;
  for (std::size_t k = 0; k < P.x.size(); k++) {
    if (!(size[P.i[k]] > 0) || !(size[P.j[k]] > 0) || !std::isfinite(size[P.i[k]]) || !std::isfinite(size[P.j[k]])) continue;
    if (method == GravityMethod::ols_sizes && !(P.x[k] > 0)) continue;
    kept.push_back(k);
    lyi.push_back(std::log(size[P.i[k]]));
    lyj.push_back(std::log(size[P.j[k]]));
    ld.push_back(P.log_d[k]);
    y.push_back(method == GravityMethod::ols_sizes ? std::log(P.x[k]) : P.x[k]);
  }
  if (kept.size() < 5) throw std::invalid_argument("economics::trade::gravity: too few pairs with sizes");
  const auto X = tyche::regression::design_matrix({lyi, lyj, ld}, true);
  Gravity g;
  g.method = method;
  g.names = {"constant", "ln size of the exporter", "ln size of the importer", "ln distance"};
  std::vector<double> fitted;
  if (method == GravityMethod::ols_sizes) {
    tyche::regression::OlsOptions o;
    o.covariance = tyche::regression::Covariance::hc1;
    const auto f = tyche::regression::ols(X, y, o);
    g.beta = f.beta, g.se = f.se, g.p = f.p, g.ci_lower = f.ci_lower, g.ci_upper = f.ci_upper;
    for (double v : f.fitted) fitted.push_back(std::exp(v));
  } else {
    tyche::regression::GlmOptions o;
    o.family = tyche::regression::Family::poisson;
    o.covariance = tyche::regression::GlmCovariance::hc0;
    const auto f = tyche::regression::glm(X, y, o);
    g.beta = f.beta, g.se = f.se, g.p = f.p, g.ci_lower = f.ci_lower, g.ci_upper = f.ci_upper;
    fitted = f.fitted;
    g.iterations = f.iterations;
    g.converged = f.converged;
  }
  g.distance_elasticity = g.beta[3];
  for (std::size_t m = 0; m < kept.size(); m++) {
    const auto k = kept[m];
    g.exporter.push_back(P.i[k]);
    g.importer.push_back(P.j[k]);
    g.flow.push_back(P.x[k]);
    g.distance_km.push_back(P.d[k]);
    g.zeros += P.x[k] == 0;
  }
  g.fitted = fitted;
  g.pairs = kept.size();
  g.pseudo_r_squared = detail::squared_correlation(g.flow, g.fitted);
  return g;
}

}  // namespace economics::trade
