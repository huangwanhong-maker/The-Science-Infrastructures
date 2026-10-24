// The engine of the Economy Atlas: one request in JSON on standard input, one answer in JSON on
// standard output (the application's README, §3). It composes the economics domain (transfer
// networks, circulation, production, the spacetime complex, places) with physis's discrete
// exterior calculus and morphe's complexes, and reads only the products whose paths the request
// names. Errors are answered as {"error": "..."} with exit status 1.
//
//   {"op": "map", "flows": path, "year": y}           the world's trade in a year: reciprocal exchange,
//                                                      and the net flows' driven, local and global parts
//   {"op": "spacetime", "record": path, "first": "YYYY-MM", "last": "YYYY-MM"}
//                                                      currencies held and exchanged: curvature, prices
//   {"op": "production", "table": {...}}               a table's Leontief model and its circular flow
//   {"op": "circuits", "network": {...}}               a transfer network's deficits, prices and circulation

#include <economics/circulation/flows.hpp>
#include <economics/places/map.hpp>
#include <economics/production/input_output.hpp>
#include <economics/spacetime/complex.hpp>
#include <economics/transfer/network.hpp>
#include <physis/discretisation/hodge.hpp>
#include <science/json.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

namespace {

using science::Json;
using Arr = Json::Array;
using Obj = Json::Object;
namespace pd = physis::discretisation;

constexpr double degree = std::numbers::pi / 180;

Json read_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error(path + ": cannot be read");
  std::ostringstream s;
  s << in.rdbuf();
  return Json::parse(s.str());
}

std::string text(const Json& request, const char* key) {
  const Json* v = request.find(key);
  if (!v || !v->is_string() || v->as_string().empty()) throw std::invalid_argument(std::string("the request has no \"") + key + "\"");
  return v->as_string();
}

double round_to(double x, double unit) { return std::round(x / unit) * unit; }

// ---- the map ---------------------------------------------------------------------------------------

Json map(const Json& request) {
  const Json product = read_file(text(request, "flows"));
  const auto& years = product["years"].as_array();
  const int year = request.has("year") ? static_cast<int>(request["year"].as_number()) : static_cast<int>(years.back().as_number());
  const Json& flows = product["flows"][std::to_string(year)];
  if (!flows.is_array()) throw std::invalid_argument("no flows for " + std::to_string(year));
  const auto& economies = product["economies"].as_array();
  const std::size_t n = economies.size();
  economics::TransferNetwork trade;
  std::vector<economics::Place> places;
  for (const auto& e : economies) {
    trade.positions.push_back(e["iso3"].as_string());
    places.push_back({e["iso3"].as_string(), e["latitude"].as_number() * degree, e["longitude"].as_number() * degree});
  }
  std::vector<double> exports(n, 0.0), imports(n, 0.0);
  for (const auto& f : flows.as_array()) {
    const auto i = static_cast<std::size_t>(f[0].as_number()), j = static_cast<std::size_t>(f[1].as_number());
    const double v = f[2].as_number();
    trade.legs.push_back({i, j, 1.0, v, "goods"});
    exports[i] += v;
    imports[j] += v;
  }
  // Reciprocal exchange, and the net flows on the map's faces. Each link's conductance is its gross
  // trade, so the driven net flow between two economies is their gross trade times their
  // difference in attraction (α is then a share of gross trade, without unit); neighbours with no
  // recorded trade are joined at the resolution of the record, the least gross trade recorded.
  const auto net = economics::net_transfer(trade);
  economics::MapFaces faces;
  const auto filled = economics::with_map_faces(net.network, places, &faces);
  double least = 0;
  for (double g : net.gross)
    if (g > 0 && (least == 0 || g < least)) least = g;
  economics::Cochain conductance = net.gross;
  conductance.resize(filled.legs.size(), least);
  const auto C = economics::circulation(filled, conductance);
  // Each economy's circulating share: of the driven and circulating flows on its links, taken in
  // absolute value, the part that circulates (bounded by 1, unlike their ratio to the net flow,
  // which the two parts may each exceed when they run against each other).
  std::vector<double> circulating(n, 0.0), through(n, 0.0);
  Arr legs;
  for (std::size_t e = 0; e < net.network.legs.size(); e++) {
    const auto& l = filled.legs[e];
    const double driven = C.parts.driven[e], local = C.parts.local[e], global = C.parts.global[e];
    for (std::size_t x : {l.from, l.to}) {
      through[x] += std::fabs(driven) + std::fabs(local + global);
      circulating[x] += std::fabs(local + global);
    }
    if (l.flow <= 0) continue;  // a link whose two directions balance has no net flow to draw
    legs.push_back(Json(Arr{Json(l.from), Json(l.to), Json(round_to(l.flow, 0.001)), Json(round_to(driven, 0.001)), Json(round_to(local, 0.001)), Json(round_to(global, 0.001))}));
  }
  Arr out;
  for (std::size_t i = 0; i < n; i++) {
    const auto& e = economies[i];
    out.push_back(Json(Obj{{"iso3", e["iso3"]}, {"name", e["name"]}, {"continent", e["continent"]}, {"latitude", e["latitude"]}, {"longitude", e["longitude"]},
                           {"exports", Json(round_to(exports[i], 0.001))}, {"imports", Json(round_to(imports[i], 0.001))},
                           {"balance", Json(round_to(C.net_inflow[i], 0.001))}, {"attraction", Json(round_to(C.attraction[i], 1e-6))},
                           {"circulating", Json(through[i] > 0 ? circulating[i] / through[i] : 0.0)}}));
  }
  Arr triangles;
  for (std::size_t f = 0; f < faces.triangles.size(); f++) {
    const auto& t = faces.triangles[f];
    triangles.push_back(Json(Obj{{"places", Json(Arr{Json(t[0]), Json(t[1]), Json(t[2])})}, {"curl", Json(round_to(C.parts.circulation[f], 0.001))}}));
  }
  return Json(Obj{{"year", Json(year)},
                  {"years", product["years"]},
                  {"unit", product["unit"]},
                  {"economies", Json(std::move(out))},
                  {"legs", Json(std::move(legs))},
                  {"triangles", Json(std::move(triangles))},
                  {"totals", Json(Obj{{"total", Json(net.total)}, {"reciprocal", Json(net.reciprocal)}, {"net", Json(net.net)},
                                      {"driven", Json(C.driven)}, {"local", Json(C.local)}, {"global", Json(C.global)}})},
                  {"faces", Json(Obj{{"triangles", Json(faces.triangles.size())}, {"neutral_legs", Json(faces.neutral_legs)}, {"covers_sphere", Json(faces.covers_sphere)},
                                     {"neutral_conductance", Json(least)}})},
                  {"method", Json(C.parts.drop.sparse ? "sparse" : "dense")},
                  {"source", product["source"]}});
}

// ---- spacetime -------------------------------------------------------------------------------------

Json spacetime(const Json& request) {
  const Json product = read_file(text(request, "record"));
  const auto& months = product["months"].as_array();
  const std::string first = request.has("first") ? request["first"].as_string() : months.front().as_string();
  const std::string last = request.has("last") ? request["last"].as_string() : months.back().as_string();
  std::vector<std::size_t> span;
  for (std::size_t m = 0; m < months.size(); m++)
    if (months[m].as_string() >= first && months[m].as_string() <= last) span.push_back(m);
  if (span.size() < 2) throw std::invalid_argument("a span of at least two months");
  // The currencies with a complete record over the span (the euro first): a rate each month and a
  // policy rate each month.
  std::vector<std::string> codes, names;
  std::vector<std::string> wanted;
  if (const Json* w = request.find("currencies"); w && w->is_array())
    for (const auto& c : w->as_array()) wanted.push_back(c.as_string());
  for (const auto& c : product["currencies"].as_array()) {
    const std::string code = c["code"].as_string();
    if (!wanted.empty() && code != "EUR" && std::find(wanted.begin(), wanted.end(), code) == wanted.end()) continue;
    bool complete = true;
    for (std::size_t m : span) {
      complete = complete && !product["policy"][code][m].is_null();
      if (code != "EUR") complete = complete && !product["rates"][code][m].is_null();
    }
    if (complete || code == "EUR") {
      if (!complete) throw std::invalid_argument("the euro's policy rate is missing in the span");
      codes.push_back(code);
      names.push_back(c["name"].as_string());
    }
  }
  economics::Transactions T;
  T.agents = codes;
  std::vector<double> times;
  for (std::size_t m : span) {
    const std::string ym = months[m].as_string();
    const double t = std::stod(ym.substr(0, 4)) + (std::stod(ym.substr(5, 2)) - 0.5) / 12;
    times.push_back(t);
    economics::Transactions::Event ev;
    ev.time = t;
    ev.label = ym;
    for (std::size_t a = 0; a < codes.size(); a++) ev.present.push_back(a);
    T.events.push_back(ev);
  }
  T.holding_rates.assign(codes.size(), std::vector<double>(span.size(), 0.0));
  for (std::size_t a = 0; a < codes.size(); a++)
    for (std::size_t k = 0; k < span.size(); k++) T.holding_rates[a][k] = std::log1p(product["policy"][codes[a]][span[k]].as_number() / 100);
  for (std::size_t a = 1; a < codes.size(); a++)
    for (std::size_t k = 0; k < span.size(); k++) T.exchanges.push_back({k, 0, a, product["rates"][codes[a]][span[k]].as_number(), 0.0});
  const auto S = economics::spacetime(T);
  // The same record with nothing earned by holding: its curvature is the exchange rates' alone,
  // κ₀ = ln S_q − ln S_p, against which the rates of holding are measured.
  economics::Transactions T0 = T;
  for (auto& r : T0.holding_rates) std::fill(r.begin(), r.end(), 0.0);
  const auto S0 = economics::spacetime(T0);
  // Lens by lens: the pair's other currency, the months, the curvature (EUR held and exchanged
  // later against exchanged at once and held: the euro's excess return over the currency).
  Arr lenses;
  std::vector<std::vector<double>> cumulative(codes.size()), ks(codes.size()), k0s(codes.size());
  std::vector<double> years(codes.size(), 0.0);
  for (std::size_t f = 0; f < S.lenses; f++) {
    const auto& L = S.worldlines.lenses[f];
    const std::size_t x = S.encounters.pairs[L.pair].b;  // the euro is agent 0, every pair's a
    const double k = S.curvature[f];
    lenses.push_back(Json(Arr{Json(x), Json(L.from), Json(L.to), Json(k)}));
    ks[x].push_back(k);
    k0s[x].push_back(S0.curvature[f]);
    years[x] += times[L.to] - times[L.from];
    cumulative[x].push_back((cumulative[x].empty() ? 0.0 : cumulative[x].back()) + k);
  }
  // Annualised: the mean excess return a year, Σ κ / Σ Δt, and the volatility, the sample standard
  // deviation of the curvatures over the square root of their mean interval. The curvature is
  // invariant under every re-description of value (a change of a currency's unit adds a gradient to
  // the log transports), so these, and the share of the exchange rates' curvature that the rates of
  // holding remove, 1 − Σ κ² / Σ κ₀², are properties of the record.
  Arr stats;
  double sum2 = 0, sum02 = 0;
  std::size_t count = 0;
  for (std::size_t x = 1; x < codes.size(); x++) {
    const std::size_t m = ks[x].size();
    double total = 0, k2 = 0, k02 = 0;
    for (std::size_t i = 0; i < m; i++) {
      total += ks[x][i];
      k2 += ks[x][i] * ks[x][i];
      k02 += k0s[x][i] * k0s[x][i];
    }
    sum2 += k2, sum02 += k02, count += m;
    const double mean_k = m ? total / static_cast<double>(m) : 0;
    double ss = 0;
    for (double k : ks[x]) ss += (k - mean_k) * (k - mean_k);
    const double sd = m > 1 ? std::sqrt(ss / static_cast<double>(m - 1)) : 0;
    const double dt = m ? years[x] / static_cast<double>(m) : 1;
    stats.push_back(Json(Obj{{"currency", Json(x)}, {"mean", Json(years[x] > 0 ? total / years[x] : 0.0)}, {"volatility", Json(sd / std::sqrt(dt))},
                             {"removed", Json(k02 > 0 ? 1 - k2 / k02 : 0.0)}, {"cumulative", Json(cumulative[x])}}));
  }
  // The nearest price system: the exact part of the transports, its potential the log price of each
  // currency in each month (up to one constant), its distance from the record the norm of the rest,
  // which re-description does not change (the exact part's own norm does, and is no measure).
  const auto H = pd::hodge_decomposition(S.network.complex(), S.network.transport());
  Arr presences;
  for (std::size_t p = 0; p < S.refined.presences.size(); p++)
    presences.push_back(Json(Arr{Json(S.refined.presences[p][0]), Json(S.refined.presences[p][1]), Json(H.potential[p])}));
  const double n = static_cast<double>(std::max<std::size_t>(count, 1));
  Arr month_labels, policy;
  for (std::size_t m : span) month_labels.push_back(months[m]);
  for (std::size_t a = 0; a < codes.size(); a++) {
    Arr r;
    for (std::size_t m : span) r.push_back(product["policy"][codes[a]][m]);
    policy.push_back(Json(std::move(r)));
  }
  Arr rates;
  for (std::size_t a = 0; a < codes.size(); a++) {
    Arr r;
    for (std::size_t m : span) r.push_back(a == 0 ? Json(1.0) : product["rates"][codes[a]][m]);
    rates.push_back(Json(std::move(r)));
  }
  return Json(Obj{{"currencies", Json(codes)}, {"names", Json(names)}, {"months", Json(std::move(month_labels))}, {"times", Json(times)},
                  {"lenses", Json(std::move(lenses))}, {"stats", Json(std::move(stats))}, {"presences", Json(std::move(presences))},
                  {"policy", Json(std::move(policy))}, {"rates", Json(std::move(rates))},
                  {"curvature", Json(Obj{{"rms", Json(std::sqrt(sum2 / n))}, {"rms_exchange", Json(std::sqrt(sum02 / n))},
                                         {"removed", Json(sum02 > 0 ? 1 - sum2 / sum02 : 0.0)},
                                         {"distance", Json(std::sqrt(H.coexact_norm2 + H.harmonic_norm2))}})},
                  {"cells", Json(Obj{{"presences", Json(S.refined.presences.size())}, {"edges", Json(S.refined.cw.edges.size())}, {"lenses", Json(S.lenses)},
                                     {"triangles", Json(S.curvature.size() - S.lenses)}})},
                  {"method", Json(H.sparse ? "sparse" : "dense")}, {"iterations", Json(H.iterations)}, {"source", product["source"]}});
}

// ---- production ------------------------------------------------------------------------------------

Json production(const Json& request) {
  const Json& t = request["table"];
  economics::InputOutputTable T;
  for (const auto& s : t["sectors"].as_array()) T.sectors.push_back(s.as_string());
  const std::size_t n = T.sectors.size();
  T.z = calx::linalg::Matrix(n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) T.z(i, j) = t["z"][i][j].as_number();
  for (std::size_t i = 0; i < n; i++) {
    T.final_demand.push_back(t["final_demand"][i].as_number());
    T.value_added.push_back(t["value_added"][i].as_number());
  }
  const auto M = economics::leontief(T);
  auto matrix = [&](const calx::linalg::Matrix& A) {
    Arr rows;
    for (std::size_t i = 0; i < n; i++) {
      Arr row;
      for (std::size_t j = 0; j < n; j++) row.push_back(Json(A(i, j)));
      rows.push_back(Json(std::move(row)));
    }
    return Json(std::move(rows));
  };
  // The net deliveries, each pair of sectors conducting in proportion to its gross deliveries (as
  // the map's links do).
  const auto net = economics::net_transfer(economics::delivery_network(T));
  const auto C = economics::circulation(net.network, net.gross);
  Arr legs;
  for (std::size_t e = 0; e < net.network.legs.size(); e++) {
    const auto& l = net.network.legs[e];
    if (l.flow <= 0) continue;
    legs.push_back(Json(Arr{Json(l.from), Json(l.to), Json(l.flow), Json(C.parts.driven[e]), Json(C.parts.local[e]), Json(C.parts.global[e])}));
  }
  Json answer(Obj{{"sectors", Json(T.sectors)}, {"A", matrix(M.A)}, {"L", matrix(M.L)}, {"output", Json(M.output)}, {"multipliers", Json(M.multipliers)},
                  {"value_added_coefficients", Json(M.value_added_coefficients)}, {"prices", Json(M.prices)}, {"backward", Json(M.backward)},
                  {"forward", Json(M.forward)}, {"productive", Json(M.productive)}, {"imbalance", Json(T.imbalance())},
                  {"deliveries", Json(Obj{{"legs", Json(std::move(legs))}, {"attraction", Json(C.attraction)}, {"reciprocal", Json(net.reciprocal)},
                                          {"net", Json(net.net)}, {"driven", Json(C.driven)}, {"local", Json(C.local)}, {"global", Json(C.global)}})}});
  if (const Json* shock = request.find("value_added_coefficients"); shock && shock->is_array()) {
    std::vector<double> c;
    for (const auto& x : shock->as_array()) c.push_back(x.as_number());
    answer.set("shocked_prices", Json(economics::prices_for(M, c)));
  }
  if (const Json* y = request.find("final_demand"); y && y->is_array()) {
    std::vector<double> d;
    for (const auto& x : y->as_array()) d.push_back(x.as_number());
    answer.set("shocked_output", Json(economics::output_for(M, d)));
  }
  return answer;
}

// ---- circuits --------------------------------------------------------------------------------------

economics::TransferNetwork network_of(const Json& j) {
  economics::TransferNetwork N;
  for (const auto& p : j["positions"].as_array()) N.positions.push_back(p.as_string());
  for (const auto& l : j["legs"].as_array())
    N.legs.push_back({static_cast<std::size_t>(l["from"].as_number()), static_cast<std::size_t>(l["to"].as_number()), l.has("ratio") ? l["ratio"].as_number() : 1.0,
                      l.has("flow") ? l["flow"].as_number() : 0.0, l["kind"].text()});
  if (const Json* faces = j.find("faces"); faces && faces->is_array())
    for (const auto& f : faces->as_array()) {
      economics::Path p;
      for (const auto& s : f.as_array()) p.push_back({static_cast<std::size_t>(s[0].as_number()), static_cast<int>(s[1].as_number())});
      N.faces.push_back(std::move(p));
    }
  return N;
}

Json path_json(const economics::Path& p) {
  Arr a;
  for (const auto& s : p) a.push_back(Json(Arr{Json(s.edge), Json(s.sign)}));
  return Json(std::move(a));
}

Json circuits(const Json& request) {
  const auto N = network_of(request["network"]);
  const auto C = economics::circuits(N);
  Arr cycles;
  for (std::size_t k = 0; k < C.chords.size(); k++) {
    Json c(Obj{{"chord", Json(C.chords[k])}, {"path", path_json(C.cycles[k])}, {"holonomy", Json(C.holonomy[k])}, {"deficit", Json(C.deficit[k])}});
    // Retention, where the cycle runs every leg forwards with positive flows.
    try {
      c.set("retention", Json(economics::retention(N, C.cycles[k])));
    } catch (const std::invalid_argument&) {
    }
    cycles.push_back(std::move(c));
  }
  const auto F = economics::circulation(N);
  const auto P = economics::pairing(N);
  Arr curvature;
  if (!N.faces.empty())
    for (double k : pd::d1(N.complex(), N.transport())) curvature.push_back(Json(k));
  return Json(Obj{{"positions", Json(N.positions)}, {"cycles", Json(std::move(cycles))}, {"log_price", Json(C.log_price)}, {"defect", Json(C.defect)},
                  {"curvature", Json(std::move(curvature))},
                  {"flows", Json(Obj{{"attraction", Json(F.attraction)}, {"net_inflow", Json(F.net_inflow)}, {"driven", Json(F.driven)}, {"local", Json(F.local)},
                                     {"global", Json(F.global)}, {"parts", Json(Obj{{"driven", Json(F.parts.driven)}, {"local", Json(F.parts.local)}, {"global", Json(F.parts.global)}})}})},
                  {"pairing", Json(Obj{{"total", Json(P.total)}, {"prices", Json(P.prices)}, {"local", Json(P.local)}, {"global", Json(P.global)}})}});
}

}  // namespace

int main() {
  std::ios::sync_with_stdio(false);
  const std::string input((std::istreambuf_iterator<char>(std::cin)), std::istreambuf_iterator<char>());
  try {
    const Json request = Json::parse(input);
    const std::string op = text(request, "op");
    Json answer;
    if (op == "map") answer = map(request);
    else if (op == "spacetime") answer = spacetime(request);
    else if (op == "production") answer = production(request);
    else if (op == "circuits") answer = circuits(request);
    else throw std::invalid_argument("no operation \"" + op + "\"");
    std::cout << answer.dump() << "\n";
    return 0;
  } catch (const std::exception& e) {
    std::cout << Json(Obj{{"error", Json(e.what())}}).dump() << "\n";
    return 1;
  }
}
