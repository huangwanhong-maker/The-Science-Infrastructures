// The toolbox's statistics and econometrics: tyche's description of samples and its regressions
// (least squares with robust covariances, generalised linear models), and nomos's econometrics
// (marginal effects, panels, instruments, differences in differences).

#pragma once

#include <nomos/econometrics/difference_in_differences.hpp>
#include <nomos/econometrics/discrete_choice.hpp>
#include <nomos/econometrics/instrumental.hpp>
#include <nomos/econometrics/panel.hpp>
#include <tyche/descriptive/rank.hpp>
#include <tyche/descriptive/summary.hpp>
#include <tyche/regression/glm.hpp>
#include <tyche/regression/linear.hpp>

#include <cstdint>
#include <numeric>
#include <set>

#include "answer.hpp"
#include "spatial.hpp"

namespace science::toolbox {

namespace tr = tyche::regression;
namespace ne = nomos::econometrics;

// ---- shared ----------------------------------------------------------------------------------------------

// The logarithms the form asked for: of the outcome, of the regressors, or both.
inline bool log_outcome(const Request& r) { const auto l = r.text("logs", "none"); return l == "y" || l == "both"; }
inline bool log_regressors(const Request& r) { const auto l = r.text("logs", "none"); return l == "x" || l == "both"; }
inline std::vector<double> logged(std::vector<double> v) {
  for (double& x : v) x = x > 0 ? std::log(x) : NaN;
  return v;
}
inline std::string named(const std::string& name, bool log) { return log ? "ln " + name : name; }

// Labels as integers for the estimators: each distinct text its own.
inline std::vector<std::int64_t> labels_of(const std::vector<std::string>& texts) {
  const auto G = groups(texts);
  return {G.of.begin(), G.of.end()};
}

// A table of coefficients.
inline Json coefficients(const std::string& title, const std::vector<std::string>& names, const std::vector<double>& beta, const std::vector<double>& se,
                         const std::vector<double>& stat, const std::vector<double>& p, const std::vector<double>& lo, const std::vector<double>& hi,
                         const std::string& statistic = "t", const std::string& note = "") {
  Arr rows;
  for (std::size_t j = 0; j < beta.size(); j++)
    rows.push_back(Json(Arr{Json(j < names.size() ? names[j] : "x" + std::to_string(j)), num(beta[j]), num(se[j]), num(stat[j]), num(p[j]), num(lo[j]), num(hi[j])}));
  return table(title, {"", "estimate", "std. error", statistic, "p", "95% from", "to"}, std::move(rows), note);
}
inline Json coefficients(const std::string& title, const std::vector<std::string>& names, const tr::LinearFit& f, const std::string& note = "") {
  return coefficients(title, names, f.beta, f.se, f.t, f.p, f.ci_lower, f.ci_upper, "t", note);
}

// At most `cap` points of a scatter, evenly spread over the rows, so that the answer stays small.
inline std::vector<std::size_t> thinned(std::size_t n, std::size_t cap = 2500) {
  std::vector<std::size_t> out;
  const double step = n > cap ? static_cast<double>(n) / static_cast<double>(cap) : 1.0;
  for (double k = 0; k < static_cast<double>(n); k += step) out.push_back(static_cast<std::size_t>(k));
  return out;
}

// The outcome and the regressors, logged as asked, in the rows where all are numbers.
struct Design {
  std::vector<double> y;
  std::vector<std::vector<double>> x;
  std::vector<std::string> y_name, x_names;
  std::vector<std::size_t> rows;
};
inline Design design(const Request& r, const char* y_role, const char* x_role, std::vector<std::vector<double>> extra = {}) {
  Design D;
  const auto yname = r.one(y_role);
  auto y = r.numbers(yname);
  if (log_outcome(r)) y = logged(y);
  std::vector<std::vector<double>> cols{y};
  for (const auto& name : r.role(x_role)) {
    cols.push_back(log_regressors(r) ? logged(r.numbers(name)) : r.numbers(name));
    D.x_names.push_back(named(name, log_regressors(r)));
  }
  if (D.x_names.empty()) throw std::invalid_argument("choose at least one regressor");
  const std::size_t base = cols.size();
  for (auto& e : extra) cols.push_back(std::move(e));
  auto C = complete(cols);
  if (C.rows.size() <= D.x_names.size() + 1) throw std::invalid_argument("too few complete rows for the regressors");
  D.y = C.columns[0];
  for (std::size_t j = 1; j < base; j++) D.x.push_back(C.columns[j]);
  for (std::size_t j = base; j < C.columns.size(); j++) D.x.push_back(C.columns[j]);  // the extra columns follow
  D.y_name = {named(yname, log_outcome(r))};
  D.rows = C.rows;
  return D;
}

// ---- statistics ----------------------------------------------------------------------------------------------

inline Json summary(const Request& r) {
  const auto names = r.role("columns");
  if (names.empty()) throw std::invalid_argument("choose the columns");
  Arr rows;
  std::vector<std::vector<double>> all;
  for (const auto& name : names) {
    const auto v = r.numbers(name);
    all.push_back(v);
    std::vector<double> x;
    for (double d : v)
      if (std::isfinite(d)) x.push_back(d);
    if (x.size() < 2) {
      rows.push_back(Json(Arr{Json(name), num(static_cast<double>(x.size()))}));
      continue;
    }
    const auto s = tyche::descriptive::summarise(x);
    rows.push_back(Json(Arr{Json(name), num(static_cast<double>(s.n)), num(s.mean), num(s.standard_deviation), num(s.min), num(s.q1), num(s.median), num(s.q3), num(s.max),
                            num(s.skewness), num(s.kurtosis)}));
  }
  Arr blocks;
  blocks.push_back(table("The columns", {"", "n", "mean", "std. dev.", "min", "Q1", "median", "Q3", "max", "skewness", "excess kurtosis"}, std::move(rows),
                         "Quartiles of type 7 (Hyndman and Fan 1996); skewness g₁ and excess kurtosis g₂. tyche::descriptive."));
  if (names.size() >= 2) {
    const auto C = complete(all);
    if (C.rows.size() >= 3) {
      const std::size_t k = names.size();
      std::vector<std::vector<double>> pearson(k, std::vector<double>(k)), spearman(k, std::vector<double>(k));
      for (std::size_t a = 0; a < k; a++)
        for (std::size_t b = 0; b < k; b++) {
          pearson[a][b] = tyche::descriptive::correlation(C.columns[a], C.columns[b]);
          spearman[a][b] = tyche::descriptive::spearman(C.columns[a], C.columns[b]);
        }
      blocks.push_back(matrix("Pearson's correlations, on the " + std::to_string(C.rows.size()) + " rows where every column is a number", names, names, pearson));
      blocks.push_back(matrix("Spearman's correlations of ranks", names, names, spearman));
      if (k == 2) blocks.push_back(scatter(names[1] + " against " + names[0], C.columns[0], C.columns[1], names[0], names[1]));
    }
  }
  return answer("Summary and correlation", std::move(blocks));
}

// ---- regression ----------------------------------------------------------------------------------------------

inline tr::Covariance covariance_of(const std::string& name) {
  if (name == "classical") return tr::Covariance::classical;
  if (name == "HC0") return tr::Covariance::hc0;
  if (name == "HC2") return tr::Covariance::hc2;
  if (name == "HC3") return tr::Covariance::hc3;
  if (name == "cluster") return tr::Covariance::cluster;
  return tr::Covariance::hc1;
}

inline Json ols(const Request& r) {
  std::vector<std::vector<double>> extra;
  const bool weighted = !r.role("weights").empty();
  if (weighted) extra.push_back(r.numbers(r.one("weights")));
  auto D = design(r, "y", "x", extra);
  tr::OlsOptions o;
  o.covariance = covariance_of(r.text("covariance", "HC1"));
  o.drop_collinear = true;
  std::vector<std::vector<double>> xs(D.x.begin(), D.x.begin() + static_cast<std::ptrdiff_t>(D.x_names.size()));
  if (weighted) o.weights = D.x.back();
  if (o.covariance == tr::Covariance::cluster) {
    if (r.role("cluster").empty()) throw std::invalid_argument("choose the column of clusters, or another kind of standard errors");
    const auto all = r.texts(r.one("cluster"));
    std::vector<std::string> mine;
    for (std::size_t i : D.rows) mine.push_back(all[i]);
    o.clusters = labels_of(mine);
  }
  const bool constant = r.flag("constant", true);
  const auto X = tr::design_matrix(xs, constant);
  const auto f = tr::ols(X, D.y, o);
  std::vector<std::string> names;
  if (constant) names.push_back("constant");
  for (const auto& n : D.x_names) names.push_back(n);
  const std::string cov = r.text("covariance", "HC1");
  Arr blocks;
  blocks.push_back(figures("The fit of " + D.y_name[0], {{"observations", static_cast<double>(f.n)}, {"R²", f.r_squared}, {"adjusted R²", f.adj_r_squared},
                                                         {"residual std. error σ̂", f.sigma}, {"F, all slopes zero", f.f_statistic}, {"its p", f.f_p_value},
                                                         {"AIC", f.aic}, {"BIC", f.bic}}));
  blocks.push_back(coefficients("Coefficients, " + (cov == "classical" ? std::string("classical") : cov == "cluster" ? "clustered (" + std::to_string(f.clusters) + " clusters)" : "robust " + cov) + " standard errors",
                                names, f, weighted ? "Weighted least squares." : ""));
  // The fit seen: against its one regressor, or the outcome against its fitted values.
  const auto pick = thinned(D.y.size());
  std::vector<double> px, py;
  if (xs.size() == 1) {
    for (std::size_t i : pick) px.push_back(xs[0][i]), py.push_back(D.y[i]);
    const auto [lo, hi] = std::minmax_element(xs[0].begin(), xs[0].end());
    const double b0 = constant ? f.beta[0] : 0, b1 = f.beta[constant ? 1 : 0];
    blocks.push_back(scatter(D.y_name[0] + " against " + D.x_names[0], px, py, D.x_names[0], D.y_name[0], {*lo, *hi}, {b0 + b1 * *lo, b0 + b1 * *hi}));
  } else {
    for (std::size_t i : pick) px.push_back(f.fitted[i]), py.push_back(D.y[i]);
    const auto [lo, hi] = std::minmax_element(f.fitted.begin(), f.fitted.end());
    blocks.push_back(scatter("The outcome against the fit (the diagonal: a perfect fit)", px, py, "fitted", D.y_name[0], {*lo, *hi}, {*lo, *hi}));
  }
  blocks.push_back(note("β̂ = (XᵀX)⁻¹Xᵀy by Householder QR; robust covariances after White (1980) and MacKinnon and White (1985), clustered after Cameron and Miller (2015). tyche::regression."));
  residual_blocks(r, D.rows, f.residuals, f.fitted, D.y_name[0], blocks);
  return answer("Linear regression", std::move(blocks));
}

inline Json glm(const Request& r) {
  const std::string family = r.text("family", "logit");
  auto D = design(r, "y", "x");
  if (family != "poisson") {
    const Json* above = r.params().find("above");
    for (double& v : D.y) {
      if (above && above->is_number()) v = v > above->as_number() ? 1.0 : 0.0;
      else if (v != 0 && v != 1) throw std::invalid_argument("a binary model needs y of 0 and 1, or a threshold above which y is 1");
    }
  } else {
    for (double v : D.y)
      if (v < 0) throw std::invalid_argument("a Poisson model needs y ≥ 0");
  }
  const auto X = tr::design_matrix(D.x, true);
  tr::GlmOptions o;
  o.family = family == "poisson" ? tr::Family::poisson : tr::Family::binomial;
  o.link = family == "probit" ? tr::Link::probit : family == "poisson" ? tr::Link::log : tr::Link::logit;
  o.covariance = r.flag("robust", true) ? tr::GlmCovariance::hc0 : tr::GlmCovariance::model;
  o.drop_collinear = true;
  const auto f = tr::glm(X, D.y, o);
  std::vector<std::string> names{"constant"};
  for (const auto& n : D.x_names) names.push_back(n);
  Arr blocks;
  double ones = 0;
  for (double v : D.y) ones += v;
  std::vector<std::pair<std::string, double>> fig{{"observations", static_cast<double>(f.n)}, {"deviance", f.deviance}, {"null deviance", f.null_deviance},
                                                   {"McFadden's pseudo R²", f.pseudo_r_squared}, {"log-likelihood", f.log_likelihood}, {"AIC", f.aic}, {"iterations", static_cast<double>(f.iterations)}};
  if (family != "poisson") fig.insert(fig.begin() + 1, {"share with y = 1", ones / static_cast<double>(f.n)});
  blocks.push_back(figures(family == "poisson" ? "Poisson (PPML)" : family == "probit" ? "Probit" : "Logit", fig));
  if (!f.converged) blocks.push_back(note("The iterations did not converge: the estimates are not to be trusted (perfect separation, perhaps)."));
  blocks.push_back(coefficients(std::string("Coefficients, ") + (r.flag("robust", true) ? "robust (HC0)" : "model-based") + " standard errors", names, f.beta, f.se, f.z, f.p,
                                f.ci_lower, f.ci_upper, f.uses_t ? "t" : "z",
                                family == "poisson" ? "e^β is the ratio of expected outcomes for a unit more of the regressor; with logs, β an elasticity." : ""));
  if (family != "poisson") {
    const auto m = ne::average_marginal_effects(f, X);
    std::vector<std::string> mn;
    std::vector<double> e, se, z, p, lo, hi;
    for (std::size_t j = 1; j < names.size(); j++) {
      mn.push_back(names[j]);
      e.push_back(m.effect[j]), se.push_back(m.se[j]), z.push_back(m.z[j]), p.push_back(m.p[j]), lo.push_back(m.ci_lower[j]), hi.push_back(m.ci_upper[j]);
    }
    blocks.push_back(coefficients("Average marginal effects on the probability that y = 1", mn, e, se, z, p, lo, hi, "z", "The delta method on the coefficients' covariance. nomos::econometrics."));
  }
  blocks.push_back(note("IRLS (McCullagh and Nelder 1989); PPML after Santos Silva and Tenreyro (2006). tyche::regression."));
  {
    std::vector<double> e(D.y.size());
    for (std::size_t i = 0; i < e.size(); i++) e[i] = D.y[i] - f.fitted[i];
    residual_blocks(r, D.rows, e, f.fitted, D.y_name[0], blocks);
  }
  return answer(family == "poisson" ? "Poisson regression" : "Binary choice", std::move(blocks));
}

// ---- panels ------------------------------------------------------------------------------------------------------

inline Json panel(const Request& r) {
  const auto unit_text = r.texts(r.one("unit"));
  const auto time = r.numbers(r.one("time"));
  auto D = design(r, "y", "x", {time});
  std::vector<std::string> units;
  std::vector<std::int64_t> period;
  std::vector<double> t = D.x.back();
  D.x.pop_back();
  for (std::size_t m = 0; m < D.rows.size(); m++) {
    units.push_back(unit_text[D.rows[m]]);
    period.push_back(static_cast<std::int64_t>(std::llround(t[m])));
  }
  const auto unit = labels_of(units);
  const auto X = calx::linalg::Rectangular::from_columns(D.x);
  ne::PanelOptions o;
  const bool two = r.text("effects", "two-way") == "two-way";
  o.effects = two ? ne::PanelEffects::two_way : ne::PanelEffects::unit;
  const auto fe = ne::within(X, D.y, unit, period, o);
  const auto G = groups(units);
  std::set<std::int64_t> periods(period.begin(), period.end());
  Arr blocks;
  blocks.push_back(figures(std::string("Fixed effects: ") + (two ? "units and periods" : "units"), {
      {"observations", static_cast<double>(D.y.size())}, {"units", static_cast<double>(G.names.size())}, {"periods", static_cast<double>(periods.size())},
      {"within R²", fe.fit.r_squared}, {"effects absorbed", static_cast<double>(fe.absorbed)}}));
  blocks.push_back(coefficients("The within estimator, standard errors clustered by unit (" + std::to_string(fe.fit.clusters) + ")", D.x_names, fe.fit,
                                "A regressor constant within units (or periods) is absorbed by the effects and has no estimate."));
  if (r.flag("random", true)) {
    const auto re = ne::random_effects(X, D.y, unit);
    ne::PanelOptions classical;
    classical.covariance = ne::PanelCovariance::classical;
    const auto h = ne::hausman(ne::within(X, D.y, unit, {}, classical), re);
    std::vector<std::string> names{"constant"};
    for (const auto& n : D.x_names) names.push_back(n);
    blocks.push_back(coefficients("Random effects (Swamy and Arora), classical standard errors", names, re.fit));
    blocks.push_back(figures("The variance of the effects, and Hausman's test", {{"σ̂ of the unit effects", re.sigma_u}, {"σ̂ of the error", re.sigma_e},
                                                                                 {"ρ, the effects' share of the variance", re.rho}, {"Hausman's H", h.statistic},
                                                                                 {"degrees of freedom", h.df}, {"p", h.p_value}}));
    blocks.push_back(note(h.p_value < 0.05 ? "Hausman's test rejects random effects (p < 0.05): the unit effects are correlated with the regressors, and only the within estimator is consistent."
                                           : "Hausman's test does not reject random effects: both estimators are consistent if the effects are independent of the regressors, and random effects is the more efficient."));
  }
  blocks.push_back(note("Wooldridge (2010), ch. 10; Baltagi (2013), chs 2–3; the two-way effects absorbed by alternating projections. nomos::econometrics."));
  if (wants_rows(r)) {
    // Each unit's effect, on each of its rows: what a place holds beyond its regressors, through the panel.
    std::map<std::int64_t, double> effect;
    for (std::size_t k = 0; k < fe.unit_labels.size() && k < fe.unit_effects.size(); k++) effect[fe.unit_labels[k]] = fe.unit_effects[k];
    double mean = 0;
    for (const auto& [label, a] : effect) mean += a / static_cast<double>(effect.size());
    std::vector<double> by_row(D.rows.size());
    for (std::size_t m = 0; m < by_row.size(); m++) by_row[m] = effect[unit[m]] - mean;
    blocks.push_back(rows_of("Each unit's effect, about their mean: what it holds beyond its regressors", D.rows, by_row, "diverging"));
    blocks.push_back(rows_of("Residuals: what " + D.y_name[0] + " holds beyond the fit and the effects", D.rows, fe.fit.residuals, "diverging"));
  }
  return answer("Panel regression", std::move(blocks));
}

inline Json iv(const Request& r) {
  const auto yname = r.one("y");
  auto y = r.numbers(yname);
  if (log_outcome(r)) y = logged(y);
  const bool lx = log_regressors(r);
  auto take = [&](const char* role) {
    std::vector<std::vector<double>> out;
    for (const auto& n : r.role(role)) out.push_back(lx ? logged(r.numbers(n)) : r.numbers(n));
    return out;
  };
  const auto W = take("exogenous"), E = take("endogenous"), Z = take("instruments");
  if (E.empty() || Z.empty()) throw std::invalid_argument("choose the endogenous regressors and their instruments");
  if (Z.size() < E.size()) throw std::invalid_argument("at least as many instruments as endogenous regressors");
  std::vector<std::vector<double>> cols{y};
  for (const auto* g : {&W, &E, &Z}) cols.insert(cols.end(), g->begin(), g->end());
  const auto C = complete(cols);
  if (C.rows.size() < W.size() + Z.size() + 3) throw std::invalid_argument("too few complete rows");
  std::size_t at = 1;
  auto slice = [&](std::size_t k) {
    std::vector<std::vector<double>> v(C.columns.begin() + static_cast<std::ptrdiff_t>(at), C.columns.begin() + static_cast<std::ptrdiff_t>(at + k));
    at += k;
    return v;
  };
  const auto w = slice(W.size()), e = slice(E.size()), z = slice(Z.size());
  ne::IvOptions o;
  o.covariance = r.flag("robust", true) ? ne::IvCovariance::hc1 : ne::IvCovariance::classical;
  // The exogenous regressors with a constant (the constant alone when there are none).
  calx::linalg::Rectangular Wx(C.rows.size(), w.size() + 1);
  for (std::size_t i = 0; i < C.rows.size(); i++) {
    Wx(i, 0) = 1;
    for (std::size_t j = 0; j < w.size(); j++) Wx(i, j + 1) = w[j][i];
  }
  const auto f = ne::two_stage_least_squares(C.columns[0], Wx, calx::linalg::Rectangular::from_columns(e), calx::linalg::Rectangular::from_columns(z), o);
  std::vector<std::string> names{"constant"};
  for (const auto& n : r.role("exogenous")) names.push_back(named(n, lx));
  std::vector<std::string> endo;
  for (const auto& n : r.role("endogenous")) endo.push_back(named(n, lx)), names.push_back(named(n, lx) + " (endogenous)");
  Arr blocks;
  blocks.push_back(coefficients(std::string("Two-stage least squares, ") + (r.flag("robust", true) ? "robust (HC1)" : "classical") + " standard errors", names, f.beta, f.se, f.t, f.p,
                                f.ci_lower, f.ci_upper));
  // Least squares on the same regressors, whose bias the instruments are there to remove.
  {
    calx::linalg::Rectangular Xo(C.rows.size(), Wx.cols + e.size());
    for (std::size_t i = 0; i < C.rows.size(); i++) {
      for (std::size_t j = 0; j < Wx.cols; j++) Xo(i, j) = Wx(i, j);
      for (std::size_t j = 0; j < e.size(); j++) Xo(i, Wx.cols + j) = e[j][i];
    }
    tr::OlsOptions oo;
    oo.covariance = tr::Covariance::hc1;
    blocks.push_back(coefficients("Least squares, for comparison: biased when a regressor moves with the error", names, tr::ols(Xo, C.columns[0], oo)));
  }
  Arr fs;
  for (std::size_t k = 0; k < f.first_stage.size(); k++)
    fs.push_back(Json(Arr{Json(endo[k]), num(f.first_stage[k].f_statistic), num(f.first_stage[k].df1), num(f.first_stage[k].df2), num(f.first_stage[k].p_value),
                          num(f.first_stage[k].partial_r_squared)}));
  blocks.push_back(table("The first stage: the excluded instruments jointly zero", {"endogenous", "F", "df₁", "df₂", "p", "partial R²"}, std::move(fs),
                         "An F below about 10 marks weak instruments (Staiger and Stock 1997; Stock and Yogo 2005), and 2SLS is then biased towards OLS."));
  if (f.overidentification_df > 0)
    blocks.push_back(figures("Overidentifying restrictions (" + std::to_string(static_cast<int>(f.overidentification_df)) + " more instruments than needed)",
                             {{"Sargan's statistic", f.sargan}, {"its p", f.sargan_p_value}, {"Hansen's J", f.hansen_j}, {"its p", f.hansen_p_value}}));
  blocks.push_back(note("Wooldridge (2010), ch. 5. A test can reject instruments; none can prove them valid, which is an argument about the world. nomos::econometrics."));
  return answer("Instrumental variables", std::move(blocks));
}

inline Json did(const Request& r) {
  const auto y = r.numbers(r.one("y")), tr_ = r.numbers(r.one("treated")), post = r.numbers(r.one("post"));
  const bool fe = !r.role("unit").empty() && !r.role("time").empty();
  const auto unit_text = r.role("unit").empty() ? std::vector<std::string>(y.size(), "") : r.texts(r.one("unit"));
  const auto time_text = r.role("time").empty() ? std::vector<std::string>(y.size(), "") : r.texts(r.one("time"));
  std::vector<double> yy, treatment;
  std::vector<int> g, t;
  std::vector<std::string> u, p;
  for (std::size_t i = 0; i < y.size(); i++) {
    if (!std::isfinite(y[i]) || (tr_[i] != 0 && tr_[i] != 1) || (post[i] != 0 && post[i] != 1)) continue;
    if (!r.role("unit").empty() && unit_text[i].empty()) continue;
    yy.push_back(y[i]);
    g.push_back(static_cast<int>(tr_[i]));
    t.push_back(static_cast<int>(post[i]));
    treatment.push_back(tr_[i] * post[i]);
    u.push_back(unit_text[i]);
    p.push_back(time_text[i]);
  }
  if (yy.size() < 8) throw std::invalid_argument("too few rows with an outcome and 0/1 for treated and after");
  const auto units = r.role("unit").empty() ? std::vector<std::int64_t>{} : labels_of(u);
  const auto d = ne::difference_in_differences(yy, g, t, units);
  Arr blocks;
  blocks.push_back(figures("The effect of the treatment", {{"δ̂, difference in differences", d.estimate}, {"std. error" + std::string(units.empty() ? " (robust)" : " (clustered by unit)"), d.se},
                                                          {"t", d.t}, {"p", d.p}, {"95% from", d.ci_lower}, {"to", d.ci_upper}}));
  Arr rows;
  for (int a = 1; a >= 0; a--)
    rows.push_back(Json(Arr{Json(a ? "treated" : "control"), num(d.mean[a][0]), num(d.mean[a][1]), num(d.mean[a][1] - d.mean[a][0]), num(d.count[a][0]), num(d.count[a][1])}));
  blocks.push_back(table("The four means", {"group", "before", "after", "change", "rows before", "rows after"}, std::move(rows),
                         "δ̂ is the treated group's change less the control group's."));
  blocks.push_back(bars("The means before and after", {"before", "after"}, {{"treated", {d.mean[1][0], d.mean[1][1]}}, {"control", {d.mean[0][0], d.mean[0][1]}}}));
  if (fe) {
    const auto w = ne::two_way_fixed_effects(yy, treatment, labels_of(u), labels_of(p));
    blocks.push_back(coefficients("Two-way fixed effects, clustered by unit", {"treated × after"}, w.fit));
  }
  blocks.push_back(note("Angrist and Pischke (2009), ch. 5. The estimate is causal only if, without the treatment, both groups would have changed alike (parallel trends): an assumption to argue, not a result. nomos::econometrics."));
  return answer("Difference in differences", std::move(blocks));
}

// β-convergence over a span, and σ-convergence year by year (Barro and Sala-i-Martin 1992).
inline Json convergence(const Request& r) {
  const auto unit = r.texts(r.one("unit"));
  const auto time = r.numbers(r.one("time")), value = r.numbers(r.one("value"));
  const int from = static_cast<int>(r.number("from", 2000)), to = static_cast<int>(r.number("to", 2019));
  if (to <= from) throw std::invalid_argument("the span must run forward");
  std::map<std::string, std::map<int, double>> byunit;
  for (std::size_t i = 0; i < unit.size(); i++)
    if (!unit[i].empty() && std::isfinite(time[i]) && std::isfinite(value[i]) && value[i] > 0) byunit[unit[i]][static_cast<int>(std::llround(time[i]))] = value[i];
  std::vector<std::string> who;
  std::vector<double> l0, g;
  for (const auto& [u, s] : byunit) {
    const auto a = s.find(from), b = s.find(to);
    if (a == s.end() || b == s.end()) continue;
    who.push_back(u);
    l0.push_back(std::log(a->second));
    g.push_back((std::log(b->second) - std::log(a->second)) / (to - from));
  }
  if (who.size() < 5) throw std::invalid_argument("too few units with a value in both " + std::to_string(from) + " and " + std::to_string(to));
  tr::OlsOptions o;
  o.covariance = tr::Covariance::hc1;
  const auto f = tr::ols(tr::design_matrix({l0}, true), g, o);
  const double T = to - from, beta = f.beta[1];
  const double lambda = 1 + beta * T > 0 ? -std::log(1 + beta * T) / T : NaN;
  Arr blocks;
  blocks.push_back(figures("β-convergence, " + std::to_string(from) + "–" + std::to_string(to), {{"β, growth a year on ln initial income", beta}, {"its std. error (HC1)", f.se[1]},
                                                                                                {"p", f.p[1]}, {"implied speed λ, a year", lambda}, {"half-life, years", std::log(2) / lambda},
                                                                                                {"units", static_cast<double>(who.size())}}));
  std::vector<double> gp;
  for (double v : g) gp.push_back(100 * v);
  const auto [lo, hi] = std::minmax_element(l0.begin(), l0.end());
  blocks.push_back(scatter("Growth against the initial level", l0, gp, "ln income in " + std::to_string(from), "growth, % a year", {*lo, *hi},
                           {100 * (f.beta[0] + beta * *lo), 100 * (f.beta[0] + beta * *hi)}));
  // σ-convergence: the spread of ln income across the same units, year by year.
  std::vector<double> years, sd;
  for (int y = from; y <= to; y++) {
    std::vector<double> v;
    for (const auto& u : who)
      if (const auto it = byunit[u].find(y); it != byunit[u].end()) v.push_back(std::log(it->second));
    if (v.size() * 10 < who.size() * 9) continue;  // a year missing more than a tenth of the units
    years.push_back(y);
    sd.push_back(tyche::descriptive::standard_deviation(v));
  }
  if (years.size() >= 2) blocks.push_back(lines("σ-convergence: the spread of ln income across the units", years, "year", "std. dev. of ln income", {{"σ", sd, {}, {}, false}}));
  std::map<std::string, double> values;
  for (std::size_t k = 0; k < who.size(); k++) values[who[k]] = 100 * g[k];
  blocks.push_back(map_values("Growth a year, " + std::to_string(from) + "–" + std::to_string(to) + ", %", values, "diverging", "%"));
  blocks.push_back(note("β < 0: poorer units grew faster. λ from β = −(1 − e^(−λT))/T (Barro and Sala-i-Martin 1992). Unconditional convergence: other differences between the units are not held fixed. tyche::regression."));
  return answer("Convergence of incomes", std::move(blocks));
}

}  // namespace science::toolbox
