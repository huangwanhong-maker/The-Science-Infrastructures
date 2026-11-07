// The toolbox's trade: the economics domain's partners and gravity of trade (economics::trade) on a year of a
// product of flows between places, the library's world trade or one's own; its answers go to the
// map as layers.

#pragma once

#include <economics/trade/gravity.hpp>
#include <economics/trade/partners.hpp>

#include <fstream>
#include <numeric>
#include <sstream>

#include "econometrics.hpp"

namespace science::toolbox {

// A year of flows as a matrix of exports, with its places.
struct Bilateral {
  std::vector<std::string> ids, names;
  std::vector<double> latitude, longitude;
  economics::trade::Rectangular x;
  int year = 0;
  std::string unit;
};

inline Json read_json(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error(path + ": cannot be read");
  std::ostringstream s;
  s << in.rdbuf();
  return Json::parse(s.str());
}

inline Bilateral bilateral(const Request& r) {
  const Json* f = r.json.find("flows");
  if (!f || !f->is_string()) throw std::invalid_argument("this tool reads flows between places: the server names their product");
  const Json product = read_json(f->as_string());
  const auto& years = product["years"].as_array();
  if (years.empty()) throw std::invalid_argument("the flows have no years");
  Bilateral B;
  B.year = static_cast<int>(r.number("year", years.back().as_number()));
  const Json& flows = product["flows"][std::to_string(B.year)];
  if (!flows.is_array()) throw std::invalid_argument("no flows for " + std::to_string(B.year));
  for (const auto& e : product["economies"].as_array()) {
    B.ids.push_back(e["iso3"].as_string());
    B.names.push_back(e.has("name") ? e["name"].as_string() : e["iso3"].as_string());
    B.latitude.push_back(e["latitude"].as_number());
    B.longitude.push_back(e["longitude"].as_number());
  }
  if (const Json* u = product.find("unit"); u && u->is_string()) B.unit = u->as_string();
  const std::size_t n = B.ids.size();
  B.x = economics::trade::Rectangular(n, n);
  for (const auto& l : flows.as_array()) {
    const auto i = static_cast<std::size_t>(l[0].as_number()), j = static_cast<std::size_t>(l[1].as_number());
    if (i < n && j < n && i != j) B.x(i, j) += l[2].as_number();
  }
  return B;
}

inline std::size_t place(const Bilateral& B, const std::string& id) {
  for (std::size_t i = 0; i < B.ids.size(); i++)
    if (B.ids[i] == id || B.names[i] == id) return i;
  std::string upper = id;
  for (char& c : upper) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  for (std::size_t i = 0; i < B.ids.size(); i++)
    if (B.ids[i] == upper) return i;
  throw std::invalid_argument("no place " + id + " in these flows");
}

// For each economy, how concentrated its exports are on a few partners.
inline Json trade_concentration(const Request& r) {
  const auto B = bilateral(r);
  const auto P = economics::trade::partner_concentration(B.x);
  const auto X = economics::trade::exports(B.x);
  std::map<std::string, double> values;
  std::vector<std::size_t> order;
  for (std::size_t i = 0; i < B.ids.size(); i++)
    if (std::isfinite(P.hhi[i])) {
      values[B.ids[i]] = P.hhi[i];
      order.push_back(i);
    }
  if (order.empty()) throw std::invalid_argument("no economy exports in " + std::to_string(B.year));
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return P.hhi[a] > P.hhi[b]; });
  std::vector<double> h;
  for (std::size_t i : order) h.push_back(P.hhi[i]);
  std::vector<double> sorted = h;
  std::sort(sorted.begin(), sorted.end());
  // The world's exports, weighted: Σ_i (X_i / W) H_i.
  double W = 0, weighted = 0;
  for (std::size_t i : order) W += X[i];
  for (std::size_t i : order) weighted += X[i] / W * P.hhi[i];
  Arr blocks;
  blocks.push_back(map_values("Partner concentration of exports, " + std::to_string(B.year) + " (Herfindahl–Hirschman)", values, "sequential"));
  blocks.push_back(figures("The economies, " + std::to_string(B.year), {{"median H", sorted[sorted.size() / 2]}, {"mean H, weighted by exports", weighted},
                                                                       {"median numbers equivalent", 1 / sorted[sorted.size() / 2]}, {"economies", static_cast<double>(order.size())}}));
  Arr rows;
  auto row = [&](std::size_t i) {
    return Json(Arr{Json(B.names[i]), num(P.hhi[i]), num(P.equivalent[i]), Json(P.first[i] < B.names.size() ? B.names[P.first[i]] : "—"), num(100 * P.first_share[i])});
  };
  for (std::size_t k = 0; k < std::min<std::size_t>(15, order.size()); k++) rows.push_back(row(order[k]));
  blocks.push_back(table("Most concentrated", {"economy", "H", "equal partners", "largest partner", "its share, %"}, std::move(rows)));
  Arr least;
  for (std::size_t k = 0; k < std::min<std::size_t>(15, order.size()); k++) least.push_back(row(order[order.size() - 1 - k]));
  blocks.push_back(table("Most spread", {"economy", "H", "equal partners", "largest partner", "its share, %"}, std::move(least)));
  blocks.push_back(note("H_i = Σ_j s_ij², s_ij j's share of i's exports (Hirschman 1945); 1/H the number of equal partners giving the same H. economics::trade on nomos::measurement."));
  return answer("Concentration of partners", std::move(blocks));
}

// How much more an economy trades with each partner than the partner's size predicts.
inline Json trade_intensity(const Request& r) {
  const auto B = bilateral(r);
  const std::size_t e = place(B, r.text("economy", "DEU"));
  const auto I = economics::trade::intensity(B.x);
  const auto G = economics::trade::two_way(B.x);
  const auto X = economics::trade::exports(B.x), M = economics::trade::imports(B.x);
  double W = 0;
  for (double v : X) W += v;
  if (!(X[e] > 0)) throw std::invalid_argument(B.names[e] + " exports nothing in " + std::to_string(B.year));
  std::map<std::string, double> out, in;
  std::vector<std::size_t> partners;
  for (std::size_t j = 0; j < B.ids.size(); j++) {
    if (j == e) continue;
    if (std::isfinite(I(e, j)) && I(e, j) > 0) out[B.ids[j]] = std::log2(I(e, j));
    if (std::isfinite(I(j, e)) && I(j, e) > 0) in[B.ids[j]] = std::log2(I(j, e));
    if (B.x(e, j) > 0) partners.push_back(j);
  }
  std::sort(partners.begin(), partners.end(), [&](std::size_t a, std::size_t b) { return B.x(e, a) > B.x(e, b); });
  Arr rows;
  for (std::size_t k = 0; k < std::min<std::size_t>(20, partners.size()); k++) {
    const std::size_t j = partners[k];
    rows.push_back(Json(Arr{Json(B.names[j]), num(100 * B.x(e, j) / X[e]), num(100 * M[j] / W), num(I(e, j)), num(I(j, e)), num(G(e, j))}));
  }
  double both = 0, gross = 0;
  for (std::size_t j = 0; j < B.ids.size(); j++)
    if (j != e) {
      gross += B.x(e, j) + B.x(j, e);
      both += 2 * std::min(B.x(e, j), B.x(j, e));
    }
  const auto P = economics::trade::partner_concentration(B.x);
  const std::string who = B.names[e], year = std::to_string(B.year);
  Arr blocks;
  blocks.push_back(map_values("Intensity of " + who + "'s exports to each partner, " + year + " (log₂: +1 twice what the partner's size predicts)", out, "diverging", "log₂"));
  blocks.push_back(map_values("Intensity of each partner's exports to " + who + ", " + year + " (log₂)", in, "diverging", "log₂"));
  blocks.push_back(figures(who + ", " + year, {{"partners it exports to", static_cast<double>(partners.size())}, {"partner concentration H", P.hhi[e]},
                                               {"its trade that goes both ways, Grubel–Lloyd", gross > 0 ? both / gross : NaN},
                                               {"the world's, for comparison", economics::trade::reciprocal_share(B.x)}}));
  blocks.push_back(table("Its largest partners", {"partner", "share of its exports, %", "partner's share of world imports, %", "intensity out", "intensity in", "two-way GL"},
                         std::move(rows), "Intensity I_ij = (x_ij/X_i)/(M_j/W) (Brown 1949; Kojima 1964); above 1, more than the partner's size predicts."));
  return answer("Trade intensity", std::move(blocks));
}

// Each economy's GDP in a year, from the toolbox's table of economies (GDP per head times
// population, the World Bank's), where the server named it; empty where it did not.
inline std::vector<double> gdp_of(const Request& r, const Bilateral& B) {
  const Json* path = r.json.find("sizes");
  if (!path || !path->is_string()) return {};
  const Json t = read_json(path->as_string());
  const Json& d = t["data"];
  if (!d.has("gdp_per_capita") || !d.has("population")) return {};
  std::map<std::string, double> gdp;
  const auto& iso = d["iso3"].as_array();
  const auto& year = d["year"].as_array();
  const auto& pc = d["gdp_per_capita"].as_array();
  const auto& pop = d["population"].as_array();
  for (std::size_t i = 0; i < iso.size(); i++)
    if (year[i].is_number() && static_cast<int>(year[i].as_number()) == B.year && pc[i].is_number() && pop[i].is_number() && iso[i].is_string())
      gdp[iso[i].as_string()] = pc[i].as_number() * pop[i].as_number();
  std::vector<double> out;
  for (const auto& id : B.ids) out.push_back(gdp.count(id) ? gdp[id] : NaN);
  return out;
}

// Flows against sizes and distance: PPML with exporter and importer effects, or with sizes.
inline Json trade_gravity(const Request& r) {
  const auto B = bilateral(r);
  constexpr double deg = 3.14159265358979323846 / 180;
  std::vector<double> lat, lon;
  for (std::size_t i = 0; i < B.ids.size(); i++) lat.push_back(B.latitude[i] * deg), lon.push_back(B.longitude[i] * deg);
  const std::string method = r.text("method", "ppml");
  const std::string year = std::to_string(B.year);
  namespace et = economics::trade;
  et::Gravity g;
  std::string sizes_note;
  if (method == "ppml") {
    g = et::gravity_effects(B.x, lat, lon);
  } else {
    auto size = gdp_of(r, B);
    std::size_t known = 0;
    for (double v : size) known += v > 0;
    if (known < B.ids.size() / 2) {
      const auto X = et::exports(B.x), M = et::imports(B.x);
      size.clear();
      for (std::size_t i = 0; i < X.size(); i++) size.push_back(X[i] + M[i]);
      sizes_note = "Sizes are each economy's total trade: the World Bank's GDP is not in the cache for " + year + ".";
    } else {
      sizes_note = "Sizes are GDP, the World Bank's GDP per head times population in " + year + "; economies without it are left out.";
    }
    g = et::gravity_sizes(B.x, lat, lon, size, method == "ols" ? et::GravityMethod::ols_sizes : et::GravityMethod::ppml_sizes);
  }
  std::vector<double> z;
  for (std::size_t j = 0; j < g.beta.size(); j++) z.push_back(g.beta[j] / g.se[j]);
  Arr blocks;
  blocks.push_back(figures("The gravity of trade in goods, " + year, {{"distance elasticity", g.distance_elasticity}, {"its std. error", g.se.back()},
                                                                     {"pairs", static_cast<double>(g.pairs)}, {"pairs without trade", static_cast<double>(g.zeros)},
                                                                     {"squared correlation of flows and fit", g.pseudo_r_squared}}));
  blocks.push_back(coefficients(method == "ppml" ? "PPML with exporter and importer effects, errors clustered by pair" : method == "ols" ? "Least squares on logarithms, robust (HC1)" : "PPML with sizes, robust (HC0)",
                                g.names, g.beta, g.se, z, g.p, g.ci_lower, g.ci_upper, method == "ols" ? "t" : "z"));
  // The flows against the model, on logarithms, for the pairs that trade.
  std::vector<double> lf, lm;
  for (std::size_t k = 0; k < g.pairs; k++)
    if (g.flow[k] > 0 && g.fitted[k] > 0) lf.push_back(std::log(g.flow[k])), lm.push_back(std::log(g.fitted[k]));
  std::vector<double> px, py;
  for (std::size_t k : thinned(lf.size())) px.push_back(lm[k]), py.push_back(lf[k]);
  if (!lm.empty()) {
    const auto [lo, hi] = std::minmax_element(lm.begin(), lm.end());
    blocks.push_back(scatter("Flows against the model (logarithms; the diagonal: the model's)", px, py, "ln expected flow", "ln flow", {*lo, *hi}, {*lo, *hi}));
  }
  if (method == "ppml") {
    std::map<std::string, double> out, in;
    double mo = 0, mi = 0, no = 0, ni = 0;
    for (std::size_t i = 0; i < B.ids.size(); i++) {
      if (std::isfinite(g.exporter_effect[i])) mo += g.exporter_effect[i], no++;
      if (std::isfinite(g.importer_effect[i])) mi += g.importer_effect[i], ni++;
    }
    for (std::size_t i = 0; i < B.ids.size(); i++) {
      if (std::isfinite(g.exporter_effect[i])) out[B.ids[i]] = g.exporter_effect[i] - mo / no;
      if (std::isfinite(g.importer_effect[i])) in[B.ids[i]] = g.importer_effect[i] - mi / ni;
    }
    blocks.push_back(map_values("Exporter effects, " + year + ": what each economy sends beyond its distances (logarithm, about the mean)", out, "diverging"));
    blocks.push_back(map_values("Importer effects, " + year + ": what each economy draws beyond its distances", in, "diverging"));
  }
  blocks.push_back(note((sizes_note.empty() ? std::string() : sizes_note + " ") +
                        "Tinbergen (1962); PPML after Santos Silva and Tenreyro (2006); with the effects, each economy's size and multilateral resistance absorbed (Anderson and van Wincoop 2003). economics::trade on nomos and tyche."));
  return answer("The gravity of trade", std::move(blocks));
}

}  // namespace science::toolbox
