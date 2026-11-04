// The toolbox's time series: tyche's autocorrelation, unit roots and cointegration, ARIMA, the
// Hodrick–Prescott filter and vector autoregressions, on columns of a table read in the order of
// their time column.

#pragma once

#include <tyche/time_series/arima.hpp>
#include <tyche/time_series/autocorrelation.hpp>
#include <tyche/time_series/filters.hpp>
#include <tyche/time_series/unit_root.hpp>
#include <tyche/time_series/var.hpp>

#include <numeric>

#include "econometrics.hpp"

namespace tools {

namespace ts = tyche::time_series;

// Series in the order of time, in the rows where the time and every series are numbers; logged
// (where every value is positive) and differenced once if asked.
struct TimeSeries {
  std::vector<double> t;
  std::vector<std::vector<double>> v;
  std::vector<std::string> names;
};
inline TimeSeries series_of(const Request& r, const std::vector<std::string>& columns, bool allow_difference) {
  std::vector<std::vector<double>> cols{r.numbers(r.one("time"))};
  for (const auto& c : columns) cols.push_back(r.numbers(c));
  auto C = complete(cols);
  std::vector<std::size_t> order(C.rows.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return C.columns[0][a] < C.columns[0][b]; });
  TimeSeries S;
  for (std::size_t i : order) S.t.push_back(C.columns[0][i]);
  const bool log = r.flag("log", true), diff = allow_difference && r.flag("difference", false);
  for (std::size_t j = 0; j < columns.size(); j++) {
    std::vector<double> v;
    for (std::size_t i : order) v.push_back(C.columns[j + 1][i]);
    std::string name = columns[j];
    if (log) {
      for (double x : v)
        if (!(x > 0)) throw std::invalid_argument(columns[j] + " has values that are not positive: leave its logarithm out");
      for (double& x : v) x = std::log(x);
      name = "ln " + name;
    }
    S.v.push_back(v);
    S.names.push_back(name);
  }
  if (diff) {
    for (auto& v : S.v) {
      std::vector<double> d;
      for (std::size_t i = 1; i < v.size(); i++) d.push_back(v[i] - v[i - 1]);
      v = d;
    }
    S.t.erase(S.t.begin());
    for (auto& n : S.names) n = "Δ " + n;
  }
  if (S.t.size() < 10) throw std::invalid_argument("a series needs at least ten complete points");
  return S;
}

inline std::vector<double> lag_axis(std::size_t lags) {
  std::vector<double> k;
  for (std::size_t i = 1; i <= lags; i++) k.push_back(static_cast<double>(i));
  return k;
}
inline std::vector<std::string> lag_labels(std::size_t lags) {
  std::vector<std::string> k;
  for (std::size_t i = 1; i <= lags; i++) k.push_back(std::to_string(i));
  return k;
}

inline Json autocorrelation(const Request& r) {
  const auto S = series_of(r, {r.one("value")}, true);
  const auto& x = S.v[0];
  const auto lags = static_cast<std::size_t>(std::clamp(r.number("lags", 24), 1.0, static_cast<double>(x.size() / 2)));
  const auto a = ts::acf(x, lags), p = ts::pacf(x, lags);
  const auto lb = ts::ljung_box(x, lags);
  std::vector<double> av(a.begin() + 1, a.end()), pv(p.begin(), p.end());
  if (pv.size() > lags) pv.erase(pv.begin());
  const double band = 1.96 / std::sqrt(static_cast<double>(x.size()));
  Arr blocks;
  blocks.push_back(figures(S.names[0], {{"points", static_cast<double>(x.size())}, {"mean", tyche::descriptive::mean(x)}, {"std. dev.", tyche::descriptive::standard_deviation(x)},
                                        {"Ljung–Box Q(" + std::to_string(lags) + ")", lb.statistic}, {"its p", lb.p_value}}));
  blocks.push_back(lines("The series", S.t, "time", S.names[0], {{S.names[0], x, {}, {}, false}}));
  blocks.push_back(bars("Autocorrelation; the band ±1.96/√n", lag_labels(lags), {{"ACF", av}}, false, band));
  blocks.push_back(bars("Partial autocorrelation, by Durbin–Levinson", lag_labels(lags), {{"PACF", pv}}, false, band));
  blocks.push_back(note(lb.p_value < 0.05 ? "Ljung–Box rejects white noise (p < 0.05): the series remembers its past." : "Ljung–Box does not reject white noise."));
  return answer("Autocorrelation", std::move(blocks));
}

inline Json unit_root(const Request& r) {
  const auto S = series_of(r, {r.one("value")}, false);
  const auto& y = S.v[0];
  const auto trend = r.text("trend", "constant");
  ts::AdfOptions o;
  o.deterministic = trend == "none" ? ts::Deterministic::none : trend == "trend" ? ts::Deterministic::trend : ts::Deterministic::constant;
  const Json* lags = r.params().find("lags");
  if (lags && lags->is_number()) {
    o.selection = ts::LagSelection::fixed;
    o.lags = static_cast<int>(lags->as_number());
  }
  const auto a = ts::adf(y, o);
  const auto k = ts::kpss(y, trend == "trend");
  Arr blocks;
  blocks.push_back(figures("Augmented Dickey–Fuller: H₀ a unit root", {{"τ", a.statistic}, {"p (MacKinnon 1994)", a.p_value}, {"lags", static_cast<double>(a.lags)},
                                                                       {"critical 1%", a.critical[0]}, {"5%", a.critical[1]}, {"10%", a.critical[2]}}));
  blocks.push_back(figures("KPSS: H₀ stationary", {{"η", k.statistic}, {"p (within 0.01–0.10)", k.p_value}, {"lags", static_cast<double>(k.lags)},
                                                   {"critical 10%", k.critical[0]}, {"5%", k.critical[1]}, {"1%", k.critical[3]}}));
  const bool adf_rejects = a.p_value < 0.05, kpss_rejects = k.p_value < 0.05;
  blocks.push_back(words("Together", {{"reading", adf_rejects && !kpss_rejects ? "stationary" : !adf_rejects && kpss_rejects ? "a unit root" : adf_rejects ? "both reject: perhaps a break, or long memory" : "neither rejects: too short to tell"}}));
  blocks.push_back(lines("The series", S.t, "time", S.names[0], {{S.names[0], y, {}, {}, false}}));
  blocks.push_back(note("Dickey and Fuller (1979), lags by AIC after Schwert's maximum unless given; Kwiatkowski, Phillips, Schmidt and Shin (1992). tyche::time_series."));
  return answer("Unit roots", std::move(blocks));
}

inline Json arima_tool(const Request& r) {
  const auto S = series_of(r, {r.one("value")}, false);
  const auto& x = S.v[0];
  ts::ArimaOrder order;
  order.p = static_cast<std::size_t>(std::max(0.0, r.number("p", 1)));
  order.d = static_cast<std::size_t>(std::max(0.0, r.number("d", 1)));
  order.q = static_cast<std::size_t>(std::max(0.0, r.number("q", 1)));
  const auto f = ts::arima(x, order);
  const auto h = static_cast<std::size_t>(std::clamp(r.number("horizon", 12), 1.0, 240.0));
  const auto fc = ts::forecast(f, h);
  std::vector<std::string> names;
  if (f.has_constant) names.push_back("μ, the mean of the differenced series");
  for (std::size_t i = 1; i <= order.p; i++) names.push_back("φ" + std::to_string(i));
  for (std::size_t i = 1; i <= order.q; i++) names.push_back("θ" + std::to_string(i));
  std::vector<double> beta;
  if (f.has_constant) beta.push_back(f.mu);
  beta.insert(beta.end(), f.phi.begin(), f.phi.end());
  beta.insert(beta.end(), f.theta.begin(), f.theta.end());
  std::vector<double> z, p, lo, hi;
  const tyche::probability::Normal N;
  for (std::size_t j = 0; j < beta.size(); j++) {
    const double se = j < f.se.size() ? f.se[j] : NaN;
    z.push_back(beta[j] / se);
    p.push_back(2 * N.survival(std::fabs(beta[j] / se)));
    lo.push_back(beta[j] - 1.96 * se);
    hi.push_back(beta[j] + 1.96 * se);
  }
  const auto lb = ts::ljung_box(f.residuals, std::min<std::size_t>(24, f.residuals.size() / 4), order.p + order.q);
  Arr blocks;
  const std::string label = "ARIMA(" + std::to_string(order.p) + ", " + std::to_string(order.d) + ", " + std::to_string(order.q) + ")";
  blocks.push_back(figures(label + " of " + S.names[0], {{"innovation variance σ²", f.sigma2}, {"log-likelihood", f.log_likelihood}, {"AIC", f.aic}, {"BIC", f.bic},
                                                         {"Ljung–Box p of the residuals", lb.p_value}}));
  blocks.push_back(coefficients("Coefficients, by conditional least squares", names, beta, std::vector<double>(f.se.begin(), f.se.end()), z, p, lo, hi, "z"));
  if (!f.is_stationary || !f.is_invertible) blocks.push_back(note(std::string("The fit is ") + (!f.is_stationary ? "not stationary" : "") + (!f.is_stationary && !f.is_invertible ? " and " : "") + (!f.is_invertible ? "not invertible" : "") + ": try another order."));
  // The history's last stretch and the forecast, on one axis of time.
  const std::size_t keep = std::min(x.size(), std::max<std::size_t>(3 * h, 36));
  const double step = (S.t.back() - S.t.front()) / static_cast<double>(S.t.size() - 1);
  std::vector<double> t, hist, mean, lower, upper;
  for (std::size_t i = x.size() - keep; i < x.size(); i++) {
    t.push_back(S.t[i]);
    hist.push_back(x[i]);
    mean.push_back(i + 1 == x.size() ? x[i] : NaN);
    lower.push_back(NaN);
    upper.push_back(NaN);
  }
  for (std::size_t k = 0; k < h; k++) {
    t.push_back(S.t.back() + step * static_cast<double>(k + 1));
    hist.push_back(NaN);
    mean.push_back(fc.mean[k]);
    lower.push_back(fc.lower[k]);
    upper.push_back(fc.upper[k]);
  }
  blocks.push_back(lines("Forecasts, with 95% intervals", t, "time", S.names[0], {{"observed", hist, {}, {}, false}, {"forecast", mean, lower, upper, true}}));
  blocks.push_back(note("Box and Jenkins; conditional least squares (Box, Jenkins and Reinsel 2008, ch. 7); intervals from the ψ weights. tyche::time_series."));
  return answer(label, std::move(blocks));
}

inline Json hp_filter(const Request& r) {
  const auto S = series_of(r, {r.one("value")}, false);
  const double lambda = r.number("lambda", 129600);
  const auto f = ts::hodrick_prescott(S.v[0], lambda);
  Arr blocks;
  blocks.push_back(lines("The series and its trend (λ = " + fixed(lambda, 0) + ")", S.t, "time", S.names[0], {{"series", S.v[0], {}, {}, false}, {"trend", f.trend, {}, {}, false}}));
  blocks.push_back(lines("The cycle, the series less its trend", S.t, "time", "cycle", {{"cycle", f.cycle, {}, {}, false}}));
  blocks.push_back(figures("The cycle", {{"std. dev.", tyche::descriptive::standard_deviation(f.cycle)}, {"largest", *std::max_element(f.cycle.begin(), f.cycle.end())},
                                         {"smallest", *std::min_element(f.cycle.begin(), f.cycle.end())}}));
  blocks.push_back(note("Hodrick and Prescott (1997); λ 1600 for quarterly data, 6.25 yearly and 129 600 monthly (Ravn and Uhlig 2002). The filter is two-sided: its last values move as data arrive. tyche::time_series."));
  return answer("Trend and cycle", std::move(blocks));
}

inline Json var_tool(const Request& r) {
  const auto names = r.role("values");
  if (names.size() < 2) throw std::invalid_argument("a VAR needs two series or more");
  const auto S = series_of(r, names, true);
  const auto Y = calx::linalg::Rectangular::from_columns(S.v);
  const Json* pj = r.params().find("p");
  const std::size_t pmax = std::min<std::size_t>(12, S.t.size() / (3 * names.size()) + 1);
  std::size_t p = pj && pj->is_number() ? static_cast<std::size_t>(std::max(1.0, pj->as_number())) : std::max<std::size_t>(1, ts::select_order(Y, pmax).aic);
  const auto v = ts::var(Y, p);
  const auto h = static_cast<std::size_t>(std::clamp(r.number("horizon", 24), 1.0, 120.0));
  const auto irf = ts::impulse_responses(v, h);
  Arr blocks;
  blocks.push_back(figures("VAR(" + std::to_string(p) + ") of " + std::to_string(names.size()) + " series", {
      {"lag order" + std::string(pj && pj->is_number() ? "" : ", by AIC"), static_cast<double>(p)}, {"observations used", static_cast<double>(v.T_effective)},
      {"AIC", v.aic}, {"BIC", v.bic}, {"spectral radius of the companion", v.spectral_radius}}));
  if (!v.is_stable) blocks.push_back(note("The VAR is not stable (a root on or outside the unit circle): difference the series, or take fewer lags."));
  std::vector<double> steps;
  for (std::size_t k = 0; k <= h; k++) steps.push_back(static_cast<double>(k));
  for (std::size_t j = 0; j < S.names.size(); j++) {
    std::vector<Series> resp;
    for (std::size_t i = 0; i < S.names.size(); i++) {
      std::vector<double> y;
      for (std::size_t k = 0; k <= h && k < irf.size(); k++) y.push_back(irf[k](i, j));
      resp.push_back({S.names[i], y, {}, {}, false});
    }
    blocks.push_back(lines("The responses to a shock of one standard deviation in " + S.names[j], steps, "steps", "response", resp));
  }
  blocks.push_back(note("Lütkepohl (2005); orthogonalised by the Cholesky factor of Σ̂_u in the order of the columns chosen, which is an assumption about which moves first. tyche::time_series."));
  return answer("Vector autoregression", std::move(blocks));
}

inline Json cointegration(const Request& r) {
  std::vector<std::string> cols{r.one("y")};
  for (const auto& x : r.role("x")) cols.push_back(x);
  if (cols.size() < 2) throw std::invalid_argument("choose y and at least one x");
  const auto S = series_of(r, cols, false);
  std::vector<std::vector<double>> xs(S.v.begin() + 1, S.v.end());
  const auto e = ts::engle_granger(S.v[0], calx::linalg::Rectangular::from_columns(xs));
  std::vector<std::string> names{"constant"};
  names.insert(names.end(), S.names.begin() + 1, S.names.end());
  Arr blocks;
  blocks.push_back(figures("Engle and Granger: H₀ no cointegration", {{"τ of the residuals", e.statistic}, {"p (MacKinnon)", e.p_value}, {"lags", static_cast<double>(e.lags)},
                                                                    {"critical 1%", e.critical[0]}, {"5%", e.critical[1]}, {"10%", e.critical[2]}}));
  blocks.push_back(coefficients("The long-run regression of " + S.names[0] + " (its standard errors are not valid for inference)", names, e.cointegrating_regression));
  blocks.push_back(lines("The residuals: the departure from the long run", S.t, "time", "residual", {{"residual", e.cointegrating_regression.residuals, {}, {}, false}}));
  blocks.push_back(note(std::string(e.p_value < 0.05 ? "Cointegrated (p < 0.05): the series share a stochastic trend and return to their long-run relation." : "No cointegration found at 5%.") +
                        " Engle and Granger (1987); MacKinnon (2010)'s critical values. tyche::time_series."));
  return answer("Cointegration", std::move(blocks));
}

}  // namespace tools
