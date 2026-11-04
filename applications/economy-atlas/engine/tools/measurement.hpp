// The toolbox's inequality and index numbers: nomos's measurement (inequality, poverty,
// concentration, decompositions) and indices (bilateral formulas, direct and chained series).

#pragma once

#include <nomos/indices/series.hpp>
#include <nomos/measurement/concentration.hpp>
#include <nomos/measurement/decomposition.hpp>
#include <nomos/measurement/inequality.hpp>
#include <nomos/measurement/poverty.hpp>

#include <numeric>

#include "answer.hpp"

namespace tools {

namespace nm = nomos::measurement;

// An amount and its optional weight and groups, in the rows where the amount and the weight are numbers.
struct Amounts {
  std::vector<double> x, w;
  std::vector<std::string> group;
  std::vector<std::size_t> rows;
};
inline Amounts amounts(const Request& r) {
  const auto value = r.numbers(r.one("value"));
  const bool weighted = !r.role("weight").empty(), grouped = !r.role("group").empty();
  const auto weight = weighted ? r.numbers(r.one("weight")) : std::vector<double>(value.size(), 1.0);
  const auto group = grouped ? r.texts(r.one("group")) : std::vector<std::string>(value.size());
  Amounts A;
  for (std::size_t i = 0; i < value.size(); i++) {
    if (!std::isfinite(value[i]) || !std::isfinite(weight[i]) || weight[i] <= 0 || value[i] < 0) continue;
    if (grouped && group[i].empty()) continue;
    A.x.push_back(value[i]);
    A.w.push_back(weight[i]);
    A.group.push_back(group[i]);
    A.rows.push_back(i);
  }
  if (A.x.size() < 2) throw std::invalid_argument("at least two rows with a non-negative amount and a positive weight");
  return A;
}

// The weighted median, the least x at which the cumulative weight reaches half.
inline double weighted_median(const std::vector<double>& x, const std::vector<double>& w) {
  std::vector<std::size_t> o(x.size());
  std::iota(o.begin(), o.end(), 0);
  std::sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return x[a] < x[b]; });
  const double half = 0.5 * std::accumulate(w.begin(), w.end(), 0.0);
  double c = 0;
  for (std::size_t i : o)
    if ((c += w[i]) >= half) return x[i];
  return x[o.back()];
}

// Inequality: the indices, the Lorenz curve, and Theil's split between and within groups.
inline Json inequality(const Request& r) {
  const auto A = amounts(r);
  const double eps = r.number("epsilon", 1);
  const auto L = nm::lorenz(A.x, A.w);
  Arr blocks;
  blocks.push_back(figures("The indices", {
      {"Gini", nm::gini(L)}, {"Theil T, GE(1)", nm::theil_t(A.x, A.w)}, {"Theil L, mean log deviation", nm::theil_l(A.x, A.w)},
      {"Atkinson A(" + fixed(eps, 2) + ")", nm::atkinson(A.x, eps, A.w)}, {"Palma ratio, top 10% / bottom 40%", nm::palma_ratio(L)},
      {"share of the top 10%", nm::top_share(L, 0.1)}, {"share of the bottom 40%", L.at(0.4)},
      {"coefficient of variation", nm::coefficient_of_variation(A.x, A.w)}, {"mean", L.mean}}));
  std::vector<double> p, l, equal;
  for (std::size_t k = 0; k < L.p.size(); k++) {
    p.push_back(100 * L.p[k]);
    l.push_back(100 * L.L[k]);
    equal.push_back(100 * L.p[k]);
  }
  blocks.push_back(lines("The Lorenz curve: the share of the amount held by the poorest share of the weight", p, "poorest share, %", "share of the amount, %",
                         {{"Lorenz curve", l, {}, {}, false}, {"equality", equal, {}, {}, true}}));
  if (!r.role("group").empty()) {
    const auto G = groups(A.group);
    const auto D = nm::decompose_generalised_entropy(A.x, G.of, 1.0, A.w);
    blocks.push_back(figures("Theil T split by group", {{"total", D.total}, {"within groups", D.within}, {"between groups", D.between},
                                                         {"share between", D.total > 0 ? D.between / D.total : NaN}}));
    Arr rows;
    std::vector<std::size_t> order(G.names.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return D.mean[a] > D.mean[b]; });
    for (std::size_t g : order)
      rows.push_back(Json(Arr{Json(G.names[g]), num(100 * D.population_share[g]), num(100 * D.amount_share[g]), num(D.mean[g]), num(D.index[g])}));
    blocks.push_back(table("The groups", {"group", "weight, %", "amount, %", "mean", "Theil T within"}, std::move(rows),
                           "Theil T = Σ_g s_g T_g + Σ_g s_g ln(μ_g/μ): within plus between (Shorrocks 1980)."));
  }
  blocks.push_back(note(std::to_string(A.x.size()) + " rows" + (r.role("weight").empty() ? "" : ", weighted") + ". nomos::measurement; Cowell (2011)."));
  return answer("Inequality", std::move(blocks));
}

// Poverty: the FGT family below a line, and the headcount as the line moves.
inline Json poverty(const Request& r) {
  const auto A = amounts(r);
  const double median = weighted_median(A.x, A.w);
  const Json* lp = r.params().find("line");
  const double z = lp && lp->is_number() ? lp->as_number() : 0.5 * median;
  if (!(z > 0)) throw std::invalid_argument("a positive poverty line");
  Arr blocks;
  blocks.push_back(figures("Below the line", {{"the line z", z}, {"headcount ratio, FGT(0)", nm::headcount_ratio(A.x, z, A.w)},
                                              {"poverty gap, FGT(1)", nm::poverty_gap_index(A.x, z, A.w)}, {"severity, FGT(2)", nm::poverty_severity(A.x, z, A.w)},
                                              {"Watts index", nm::watts_index(A.x, z, A.w)}, {"the median", median}}));
  std::vector<double> lines_z, h, gap;
  for (int k = 1; k <= 60; k++) {
    const double zz = 2 * z * k / 60;
    lines_z.push_back(zz);
    h.push_back(100 * nm::headcount_ratio(A.x, zz, A.w));
    gap.push_back(100 * nm::poverty_gap_index(A.x, zz, A.w));
  }
  blocks.push_back(lines("As the line moves: the poverty incidence and gap curves", lines_z, "the line", "%",
                         {{"headcount", h, {}, {}, false}, {"gap", gap, {}, {}, false}}));
  blocks.push_back(note("Foster, Greer and Thorbecke (1984); Watts (1968). With no line given, half the weighted median."));
  return answer("Poverty", std::move(blocks));
}

// Concentration of sizes among holders.
inline Json concentration(const Request& r) {
  const auto value = r.numbers(r.one("value"));
  const auto names = r.role("label").empty() ? std::vector<std::string>{} : r.texts(r.one("label"));
  std::vector<double> x;
  std::vector<std::string> who;
  for (std::size_t i = 0; i < value.size(); i++)
    if (std::isfinite(value[i]) && value[i] >= 0) {
      x.push_back(value[i]);
      who.push_back(names.empty() ? "#" + std::to_string(i + 1) : names[i]);
    }
  if (x.size() < 2) throw std::invalid_argument("at least two sizes");
  const auto k = static_cast<std::size_t>(std::max(1.0, r.number("k", 4)));
  const double hhi = nm::herfindahl_hirschman(x);
  Arr blocks;
  blocks.push_back(figures("The indices", {{"Herfindahl–Hirschman, Σ s²", hhi}, {"the same × 10 000", 1e4 * hhi}, {"normalised", nm::herfindahl_hirschman_normalised(x)},
                                           {"numbers equivalent, 1/H", 1 / hhi}, {"CR" + std::to_string(k) + ", the largest " + std::to_string(k) + "'s share", nm::concentration_ratio(x, std::min(k, x.size()))},
                                           {"entropy of the shares", nm::entropy_index(x)}, {"holders", static_cast<double>(x.size())}}));
  const auto s = nm::shares(x);
  std::vector<std::size_t> o(x.size());
  std::iota(o.begin(), o.end(), 0);
  std::sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return s[a] > s[b]; });
  std::vector<std::string> labels;
  std::vector<double> top;
  for (std::size_t n = 0; n < std::min<std::size_t>(15, o.size()); n++) {
    labels.push_back(who[o[n]]);
    top.push_back(100 * s[o[n]]);
  }
  blocks.push_back(bars("The largest shares, %", labels, {{"share", top}}));
  blocks.push_back(note("The US merger guidelines (2023) call a market above 1 800 on the 10 000 scale highly concentrated. nomos::measurement."));
  return answer("Concentration", std::move(blocks));
}

// Index numbers from a long table of periods, items, prices and quantities.
inline Json price_index(const Request& r) {
  const auto period = r.texts(r.one("period")), item = r.texts(r.one("item"));
  const auto price = r.numbers(r.one("price")), quantity = r.numbers(r.one("quantity"));
  // Periods in order: by number when every label is one, else as they come.
  auto P = groups(period);
  bool numeric = true;
  for (const auto& n : P.names) numeric = numeric && !numbers_in(n).empty() && numbers_in(n).size() == 1;
  std::vector<std::size_t> order(P.names.size());
  std::iota(order.begin(), order.end(), 0);
  if (numeric) std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return numbers_in(P.names[a])[0] < numbers_in(P.names[b])[0]; });
  std::vector<std::size_t> rank(order.size());
  for (std::size_t k = 0; k < order.size(); k++) rank[order[k]] = k;
  const auto I = groups(item);
  const std::size_t T = P.names.size(), N = I.names.size();
  if (T < 2 || N < 1) throw std::invalid_argument("at least two periods and one item");
  nomos::indices::Periods p(T, std::vector<double>(N, NaN)), q = p;
  for (std::size_t i = 0; i < period.size(); i++) {
    p[rank[P.of[i]]][I.of[i]] = price[i];
    q[rank[P.of[i]]][I.of[i]] = quantity[i];
  }
  for (std::size_t t = 0; t < T; t++)
    for (std::size_t n = 0; n < N; n++)
      if (!std::isfinite(p[t][n]) || !std::isfinite(q[t][n]))
        throw std::invalid_argument("every item needs a price and a quantity in every period: " + I.names[n] + " in " + P.names[order[t]]);
  const bool chained = r.flag("chained", false);
  using F = nomos::indices::Formula;
  const std::vector<std::pair<std::string, F>> formulas = {{"Laspeyres", F::laspeyres}, {"Paasche", F::paasche}, {"Fisher", F::fisher}, {"Törnqvist", F::tornqvist}};
  auto series = [&](F f, const nomos::indices::Periods& a, const nomos::indices::Periods& b) {
    return chained ? nomos::indices::chained_series(f, a, b) : nomos::indices::direct_series(f, a, b);
  };
  std::vector<double> x;
  std::vector<std::string> names;
  for (std::size_t t = 0; t < T; t++) {
    names.push_back(P.names[order[t]]);
    x.push_back(numeric ? numbers_in(P.names[order[t]])[0] : static_cast<double>(t));
  }
  std::vector<Series> ps, qs;
  Arr rows;
  std::vector<std::vector<double>> pv, qv;
  for (const auto& [name, f] : formulas) {
    auto s = series(f, p, q);
    auto v = series(f, q, p);  // a quantity index is the price formula with prices and quantities exchanged
    for (double& y : s) y *= 100;
    for (double& y : v) y *= 100;
    ps.push_back({name, s, {}, {}, false});
    qs.push_back({name, v, {}, {}, false});
    pv.push_back(s);
    qv.push_back(v);
  }
  for (std::size_t t = 0; t < T; t++) {
    Arr row{Json(names[t])};
    for (const auto& s : pv) row.push_back(num(s[t]));
    for (const auto& s : qv) row.push_back(num(s[t]));
    rows.push_back(Json(std::move(row)));
  }
  Arr blocks;
  blocks.push_back(lines(std::string("Price indices, the first period = 100") + (chained ? ", chained" : ", direct"), x, numeric ? "period" : "period (in order)", "", ps));
  blocks.push_back(lines(std::string("Quantity indices, the first period = 100") + (chained ? ", chained" : ", direct"), x, numeric ? "period" : "period (in order)", "", qs));
  blocks.push_back(table("The indices", {"period", "P Laspeyres", "P Paasche", "P Fisher", "P Törnqvist", "Q Laspeyres", "Q Paasche", "Q Fisher", "Q Törnqvist"}, std::move(rows),
                         "Fisher's price index times Fisher's quantity index is the ratio of values (the factor reversal test). nomos::indices."));
  return answer("Index numbers", std::move(blocks));
}

}  // namespace tools
