// The operations of the plug-in "Landscape Dynamics": the urban-science domain's landscape module on a district
// the page sends (its walking network from the city model, two stores), on a building it sends (its box), and on
// declared parameters. Every number is computed here; the page draws.

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <tuple>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/integrators/theta.hpp>
#include <calx/linalg/sparse.hpp>
#include <morphe/formal/traces.hpp>
#include <morphe/topology/complex_rewriting.hpp>
#include <physis/chemistry/doi_peliti.hpp>
#include <physis/discretisation/hodge_laplacian.hpp>
#include <physis/dynamics/elimination.hpp>
#include <physis/statistical_mechanics/jump_process.hpp>
#include <science/json.hpp>
#include <architecture/environment/heat.hpp>
#include <architecture/space/rooms.hpp>
#include <urban/landscape/district.hpp>
#include <urban/landscape/process.hpp>
#include <urban/landscape/skeleton.hpp>

namespace ld {

using science::Json;
namespace ul = urban::landscape;
namespace as = architecture::space;
namespace ae = architecture::environment;
namespace mt = morphe::topology;
namespace tr = morphe::formal::traces;
namespace sm = physis::statistical_mechanics;

inline double num(const Json& j, const char* key, double fallback) { return j[key].is_number() ? j[key].as_number() : fallback; }
inline Json round(double x, int digits = 4) {
  if (!std::isfinite(x)) return Json();
  const double p = std::pow(10.0, digits);
  return Json(std::round(x * p) / p);
}
inline Json point(const ul::Point& p) { return Json::array({round(p[0], 2), round(p[1], 2)}); }
inline Json line(const std::vector<ul::Point>& l) {
  Json a = Json::array();
  for (const auto& p : l) a.push(point(p));
  return a;
}
inline Json numbers(const std::vector<double>& v, int digits = 5) {
  Json a = Json::array();
  for (double x : v) a.push(round(x, digits));
  return a;
}
inline Json counts(const std::vector<std::size_t>& v) {
  Json a = Json::array();
  for (auto x : v) a.push(Json(static_cast<double>(x)));
  return a;
}
inline std::string id_text(std::uint64_t id) { return std::to_string(id); }

// ---- the district ----------------------------------------------------------------------------------------------

// The district of a request; with "planarise", of every two segments that cross without a junction the shorter is
// dropped (and named in `dropped`).
inline ul::District district_of(const Json& r, std::vector<std::string>* dropped = nullptr) {
  ul::Network N;
  for (const auto& p : r["junctions"].as_array()) N.junctions.push_back({p[0].as_number(), p[1].as_number()});
  for (const auto& s : r["segments"].as_array()) {
    ul::Segment g;
    g.u = static_cast<std::size_t>(s["u"].as_number()), g.v = static_cast<std::size_t>(s["v"].as_number());
    if (s["line"].is_array())
      for (const auto& p : s["line"].as_array()) g.line.push_back({p[0].as_number(), p[1].as_number()});
    g.kind = s["kind"].is_string() ? s["kind"].text() : "street";
    if (g.kind.find('|') != std::string::npos || g.kind.empty()) g.kind = "street";
    g.name = s["name"].is_string() ? s["name"].text() : "";
    N.segments.push_back(std::move(g));
  }
  if (N.segments.size() > 2500) throw std::invalid_argument("at most 2500 segments in a district");
  if (r["planarise"].is_bool() && r["planarise"].as_bool()) N = ul::planarised(std::move(N), dropped);
  auto D = ul::district(N, num(r, "temperature", 18));
  if (r["stores"].is_array())
    for (const auto& s : r["stores"].as_array())
      ul::attach_store(D, {s["name"].is_string() ? s["name"].text() : "a store", s["use"].is_string() ? s["use"].text() : "", {s["at"][0].as_number(), s["at"][1].as_number()}});
  return D;
}

// The cells of a district with their geometry.
inline Json cells(const ul::District& D) {
  Json v = Json::array(), e = Json::array(), f = Json::array();
  const auto up = D.complex.cofaces();
  for (std::size_t i = 0; i < D.complex.size(); i++) {
    const auto& c = D.complex.cells[i];
    if (c.dim == 0) v.push(Json::object({{"id", id_text(c.id)}, {"type", c.type}, {"at", point(D.point.at(c.id))}, {"degree", static_cast<double>(up[i].size())}}));
    if (c.dim == 1) {
      const auto ends = D.complex.ends(i);
      e.push(Json::object({{"id", id_text(c.id)}, {"type", c.type}, {"name", D.name.count(c.id) ? D.name.at(c.id) : ""}, {"line", line(D.line.at(c.id))},
                           {"tail", id_text(D.complex.cells[ends[0]].id)}, {"head", id_text(D.complex.cells[ends[1]].id)}, {"blocks", static_cast<double>(up[i].size())},
                           {"length", round(D.length(c.id), 1)}}));
    }
    if (c.dim == 2)
      f.push(Json::object({{"id", id_text(c.id)}, {"type", c.type}, {"ring", line(D.ring(c.id))}, {"area", round(D.area(c.id), 1)},
                           {"edges", static_cast<double>(c.facets.size())}, {"temperature", round(D.value(c.id, "temperature"), 3)}}));
  }
  return Json::object({{"vertices", v}, {"edges", e}, {"blocks", f}});
}

inline Json topology(const ul::District& D) {
  const auto a = mt::admissible(D.complex);
  return Json::object({{"counts", counts(D.complex.counts())}, {"betti", counts(morphe::algebra::betti(mt::chain_complex(D.complex)))},
                       {"euler", static_cast<double>(D.complex.euler_characteristic())}, {"regular", a.regular}, {"reason", a.reason},
                       {"holes", static_cast<double>(D.holes)}, {"holeArea", round(D.hole_area, 0)}});
}

inline Json stores_json(const ul::District& D) {
  Json s = Json::array();
  for (const auto& st : D.stores)
    s.push(Json::object({{"name", st.name}, {"use", st.use}, {"at", point(st.at)}, {"segment", id_text(st.segment)}, {"street", D.name.count(st.segment) ? D.name.at(st.segment) : ""},
                         {"s", round(st.s, 1)}, {"offset", round(st.offset, 1)}}));
  return s;
}

inline Json plan_json(const ul::LanePlan& P) {
  if (!P.refused.empty()) return Json::object({{"refused", P.refused}});
  return Json::object({{"block", id_text(P.block)}, {"from", point(P.a)}, {"to", point(P.b)}, {"length", round(std::hypot(P.b[0] - P.a[0], P.b[1] - P.a[1]), 1)}});
}

inline Json op_district(const Json& r) {
  std::vector<std::string> dropped;
  const auto D = district_of(r, &dropped);
  Json gone = Json::array();
  for (const auto& d : dropped) gone.push(d);
  Json out = Json::object({{"cells", cells(D)}, {"topology", topology(D)}, {"stores", stores_json(D)}, {"dropped", gone}});
  if (D.stores.size() >= 2) out.set("lane", plan_json(ul::plan_lane(D, 0, 1)));
  return out;
}

// Of the stores sent (candidates), the pairs a lane could join: on different segments of one block, the chord inside
// it, nearest first; at most `limit`.
inline Json op_pairs(const Json& r) {
  auto D = district_of(r);
  const std::size_t n = D.stores.size(), limit = static_cast<std::size_t>(num(r, "limit", 12));
  if (n > 600) throw std::invalid_argument("at most 600 candidate stores");
  std::vector<std::tuple<double, std::size_t, std::size_t, ul::LanePlan>> found;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 1; j < n; j++) {
      const double d = std::hypot(D.stores[i].at[0] - D.stores[j].at[0], D.stores[i].at[1] - D.stores[j].at[1]);
      if (d > num(r, "within", 160) || D.stores[i].offset > 25 || D.stores[j].offset > 25) continue;
      auto plan = ul::plan_lane(D, i, j);
      // A lane at least `least` metres long: shorter ones only cut a block's corner.
      if (plan.refused.empty() && std::hypot(plan.a[0] - plan.b[0], plan.a[1] - plan.b[1]) >= num(r, "least", 0)) found.push_back({d, i, j, plan});
    }
  std::sort(found.begin(), found.end(), [](const auto& a, const auto& b) { return std::get<0>(a) < std::get<0>(b); });
  Json out = Json::array();
  for (std::size_t k = 0; k < found.size() && k < limit; k++) {
    const auto& [d, i, j, plan] = found[k];
    out.push(Json::object({{"a", static_cast<double>(i)}, {"b", static_cast<double>(j)}, {"distance", round(d, 1)}, {"lane", plan_json(plan)}}));
  }
  return Json::object({{"pairs", out}, {"candidates", static_cast<double>(n)}});
}

inline Json footprint_json(const tr::Footprint& f) {
  Json r = Json::array(), w = Json::array();
  for (const auto& x : f.reads) r.push(x);
  for (const auto& x : f.writes) w.push(x);
  return Json::object({{"reads", r}, {"writes", w}});
}
template <class Step>
inline Json step_json(const Step& s) {
  Json c = Json::array(), d = Json::array();
  for (auto id : s.created) c.push(id_text(id));
  for (auto id : s.deleted) d.push(id_text(id));
  return Json::object({{"kind", s.kind}, {"match", s.match}, {"created", c}, {"deleted", d}, {"refinement", s.refinement},
                       {"reads", static_cast<double>(s.footprint.reads.size())}, {"writes", static_cast<double>(s.footprint.writes.size())}});
}
inline Json order_json(const tr::CausalOrder& P, bool count = true) {
  Json covers = Json::array(), steps = Json::array();
  for (const auto& [a, b] : P.covers) covers.push(Json::array({static_cast<double>(a), static_cast<double>(b)}));
  for (const auto& s : P.steps) {
    Json l = Json::array();
    for (auto e : s) l.push(Json(static_cast<double>(e)));
    steps.push(l);
  }
  Json out = Json::object({{"covers", covers}, {"steps", steps}});
  if (count && P.n <= 40) out.set("records", tr::linearisations(P).str());
  return out;
}

// Builds the lane between the first two stores, and closes it again.
inline Json op_lane(const Json& r) {
  auto D = district_of(r);
  if (D.stores.size() < 2) throw std::invalid_argument("two stores");
  const auto plan = ul::plan_lane(D, 0, 1);
  if (!plan.refused.empty()) return Json::object({{"lane", plan_json(plan)}});
  const auto before = topology(D);
  const ul::District original = D;
  const double C0 = [&] { double s = 0; for (auto id : D.of(2)) s += D.value(id, "capacity"); return s; }();
  const auto steps = ul::build_lane(D, plan);
  const auto after = topology(D);
  const double C1 = [&] { double s = 0; for (auto id : D.of(2)) s += D.value(id, "capacity"); return s; }();
  Json st = Json::array();
  std::vector<tr::Footprint> fps;
  for (const auto& s : steps) st.push(step_json(s)), fps.push_back(s.footprint);
  const std::uint64_t lane = steps.back().created[0];
  auto closed = D;
  const auto undo = ul::close_lane(closed, lane);
  Json us = Json::array();
  std::vector<tr::Footprint> ufps;
  for (const auto& s : undo) us.push(step_json(s)), ufps.push_back(s.footprint);
  return Json::object({{"lane", plan_json(plan)}, {"laneId", id_text(lane)}, {"steps", st}, {"order", order_json(tr::causal_order(fps))},
                       {"before", before}, {"after", after}, {"cells", cells(D)}, {"stores", stores_json(D)},
                       {"capacityKept", round((C1 - C0) / C0, 15)},
                       {"undo", Json::object({{"steps", us}, {"order", order_json(tr::causal_order(ufps))}, {"isomorphic", mt::isomorphic(closed.complex, original.complex)},
                                                {"again", plan_json(ul::plan_lane(closed, 0, 1))}, {"stores", stores_json(closed)}})}});
}

// ---- the landscape running -------------------------------------------------------------------------------------

inline Json op_run(const Json& r) {
  auto D = district_of(r);
  if (D.stores.size() < 2) throw std::invalid_argument("two stores");
  ul::LandscapeParams P;
  const Json& q = r["params"];
  P.alpha = num(q, "alpha", P.alpha), P.beta = num(q, "beta", P.beta), P.eta = num(q, "eta", P.eta), P.xi = num(q, "xi", P.xi);
  P.open_from = num(q, "openFrom", P.open_from), P.open_to = num(q, "openTo", P.open_to);
  P.permit_rate = num(q, "permit", P.permit_rate), P.build_rate = num(q, "build", P.build_rate), P.remove_rate = num(q, "remove", P.remove_rate);
  P.recognise_rate = num(q, "recognise", P.recognise_rate), P.recognise_reference = num(q, "recogniseReference", P.recognise_reference);
  P.attraction = num(q, "attraction", P.attraction), P.stay = num(q, "stay", P.stay), P.air = num(q, "air", P.air), P.sun_gain = num(q, "sunGain", P.sun_gain);
  P.seed = static_cast<std::uint64_t>(num(q, "seed", 1));
  if (r["entries"].is_array())
    for (const auto& x : r["entries"].as_array()) P.entries.push_back(x.as_number());
  if (r["altitude"].is_array())
    for (const auto& x : r["altitude"].as_array()) P.altitude.push_back(x.as_number());
  if (!P.entries.empty() && P.entries.size() != 24) throw std::invalid_argument("entries: 24 values, one each hour");
  if (!P.altitude.empty() && P.altitude.size() != 24) throw std::invalid_argument("altitude: 24 values, one each hour");
  for (const auto id : D.of(0)) D.set(id, "walkers", 0.0);
  const auto run = ul::run_landscape(D, P);
  Json events = Json::array();
  std::vector<tr::Footprint> fps;
  for (const auto& e : run.events) {
    Json rw = Json::array();
    for (const auto& s : e.rewrites) rw.push(step_json(s));
    events.push(Json::object({{"time", round(e.time, 4)}, {"mechanism", e.mechanism}, {"match", e.match}, {"forced", e.forced}, {"rewrites", rw}, {"footprint", footprint_json(e.footprint)}}));
    fps.push_back(e.footprint);
  }
  // The configurations: the strata (one for each interval between events) with the same cells are one
  // configuration, sent once with the intervals it held; a frame names its configuration.
  std::map<std::string, std::size_t> config;
  std::vector<std::size_t> config_of, first;
  std::vector<std::vector<std::pair<double, double>>> held;
  for (std::size_t k = 0; k < run.strata.size(); k++) {
    const auto& S = run.strata[k];
    std::string key;
    for (const auto& c : S.district.complex.cells) key += id_text(c.id) + ",";
    const auto [it, fresh] = config.emplace(key, first.size());
    if (fresh) first.push_back(k), held.emplace_back();
    config_of.push_back(it->second);
    auto& h = held[it->second];
    if (!h.empty() && std::fabs(h.back().second - S.from) < 1e-9) h.back().second = S.to;  // contiguous: one interval
    else h.push_back({S.from, S.to});
  }
  Json frames = Json::array();
  for (const auto& F : run.frames) {
    Json w = Json::object(), f = Json::object(), t = Json::object();
    for (const auto& [id, n] : F.walkers) w.set(id_text(id), round(n, 2));
    for (const auto& [id, x] : F.flow) f.set(id_text(id), round(x, 3));
    Json fw = Json::object(), bw = Json::object();
    for (const auto& [id, x] : F.forth) fw.set(id_text(id), round(x, 2));
    for (const auto& [id, x] : F.back) bw.set(id_text(id), round(x, 2));
    for (const auto& [id, x] : F.temperature) t.set(id_text(id), round(x, 3));
    Json en = Json::array();
    for (const auto& k : F.enabled) en.push(k);
    frames.push(Json::object({{"time", round(F.time, 4)}, {"walkers", w}, {"flow", f}, {"temperature", t}, {"a", F.a}, {"b", F.b}, {"permit", F.permit},
                              {"recognised", F.recognised}, {"lane", F.lane ? id_text(F.lane) : ""}, {"stratum", static_cast<double>(config_of[F.stratum])}, {"enabled", en},
                              {"forth", fw}, {"back", bw}}));
  }
  Json strata = Json::array();
  for (std::size_t c = 0; c < first.size(); c++) {
    const auto& S = run.strata[first[c]];
    Json iv = Json::array();
    for (const auto& [a, b] : held[c]) iv.push(Json::array({round(a, 4), round(b, 4)}));
    strata.push(Json::object({{"from", round(S.from, 4)}, {"to", round(S.to, 4)}, {"intervals", iv}, {"lane", S.lane ? id_text(S.lane) : ""}, {"betti", counts(S.betti)},
                              {"counts", counts(S.district.complex.counts())}, {"regular", mt::admissible(S.district.complex).regular}, {"cells", cells(S.district)},
                              {"stores", stores_json(S.district)}}));
  }
  const auto order = tr::causal_order(fps);
  return Json::object({{"events", events}, {"frames", frames}, {"strata", strata}, {"order", order_json(order, run.events.size() <= 40)},
                       {"balance", Json::object({{"start", round(run.walkers0, 3)}, {"end", round(run.walkers1, 3)}, {"entered", round(run.entered, 3)}, {"left", round(run.left, 3)},
                                                 {"residual", round(run.walkers1 - run.walkers0 - (run.entered - run.left), 9)}})},
                       {"laneRefused", run.lane_refused}});
}

// ---- the landscape's event skeleton -------------------------------------------------------------------------

inline ul::SkeletonParams skeleton_of(const Json& r) {
  ul::SkeletonParams q;
  q.alpha_a = num(r, "alphaA", q.alpha_a), q.alpha_b = num(r, "alphaB", q.alpha_b), q.beta_a = num(r, "betaA", q.beta_a), q.beta_b = num(r, "betaB", q.beta_b);
  q.eta_a = num(r, "etaA", q.eta_a), q.eta_b = num(r, "etaB", q.eta_b), q.xi_a = num(r, "xiA", q.xi_a), q.xi_b = num(r, "xiB", q.xi_b);
  q.grant = num(r, "grant", q.grant), q.revoke = num(r, "revoke", q.revoke), q.build = num(r, "build", q.build), q.remove = num(r, "remove", q.remove);
  q.recognise = num(r, "recognise", q.recognise), q.lapse = num(r, "lapse", q.lapse), q.use = num(r, "use", q.use), q.draw = num(r, "draw", q.draw);
  for (double x : {q.alpha_a, q.alpha_b, q.beta_a, q.beta_b, q.grant, q.revoke, q.build, q.remove, q.recognise, q.lapse, q.use})
    if (!(x > 0)) throw std::invalid_argument("every rate of the skeleton, and the lane's use, must be positive: each change needs its way back");
  for (double x : {q.eta_a, q.eta_b, q.xi_a, q.xi_b, q.draw})
    if (!(x >= 0)) throw std::invalid_argument("how much a store answers the other or the lane, and how much it draws to the lane, is ≥ 0");
  return q;
}

inline Json op_skeleton(const Json& r) {
  const auto q = skeleton_of(r);
  const auto a = ul::analyse_skeleton(q);
  Json states = Json::array();
  for (std::size_t i = 0; i < a.system.states.size(); i++) {
    const auto& z = a.system.states[i];
    Json st = Json::object();
    for (const auto& part : ul::skeleton_parts()) st.set(part, static_cast<double>(z.at(part)));
    st.set("p", round(a.pi[i], 7));
    states.push(st);
  }
  Json moves = Json::array();
  for (const auto& m : a.system.moves) moves.push(Json::object({{"mechanism", m.mechanism}, {"from", static_cast<double>(m.from)}, {"to", static_cast<double>(m.to)}, {"rate", round(m.rate, 6)}}));
  Json cycles = Json::array();
  for (const auto& c : a.cycles) {
    Json st = Json::array();
    for (auto x : c.states) st.push(Json(static_cast<double>(x)));
    cycles.push(Json::object({{"states", st}, {"affinity", round(c.affinity, 6)}, {"current", round(c.current, 8)}}));
  }
  Json views = Json::array();
  for (const auto& v : a.views)
    views.push(Json::object({{"name", v.name}, {"blocks", static_cast<double>(v.blocks)}, {"eventwise", v.eventwise}, {"whole", v.whole}, {"discrepancy", round(v.discrepancy, 8)},
                             {"horizon", std::isfinite(v.horizon) ? round(v.horizon, 4) : Json()}, {"reason", v.reason}, {"mechanism", v.mechanism},
                             {"state", static_cast<double>(v.state)}, {"other", static_cast<double>(v.other)}, {"by", round(v.by, 6)}}));
  // The time both stores stay shut, under the description that keeps the stores alone.
  const auto alone = ul::skeleton_views()[0];
  const auto blk = ul::view_blocks(a.system, alone);
  std::size_t shut = 0;
  for (std::size_t i = 0; i < blk.size(); i++)
    if (a.system.states[i].at("A") == 0 && a.system.states[i].at("B") == 0) shut = blk[i];
  const auto dwell = ul::skeleton_dwell(a, alone, shut);
  const double mean = dwell.moment(1);
  std::vector<double> tau, density, exponential;
  for (double t = 0; t <= 4 * mean + 1e-12; t += mean / 25) tau.push_back(t), density.push_back(dwell.density(t)), exponential.push_back(std::exp(-t / mean) / mean);
  // Interventions: promote A, promote B, grant permits faster; on three outcomes.
  const double gain = num(r, "gain", 0.5);
  Json interventions = Json::object();
  for (const std::string outcome : {"both", "open", "recognised"}) interventions.set(outcome, numbers(ul::intervention_interactions(q, gain, outcome), 9));
  // Evidence: a record of the skeleton, the posterior over a common coupling η.
  const double hours = num(r, "hours", 500);
  const auto rec = ul::simulate_skeleton(q, hours, static_cast<std::uint64_t>(num(r, "seed", 5)));
  const std::vector<double> etas{0, 0.25, 0.5, 1, 1.5, 2, 3};
  const auto post = ul::coupling_posterior(q, rec, etas, std::vector<double>(etas.size(), 1.0 / etas.size()));
  return Json::object({
      {"states", states}, {"moves", moves}, {"coherent", a.coherent}, {"reversibleRates", a.reversible_rates},
      {"entropyProduction", std::isfinite(a.entropy_production) ? round(a.entropy_production, 8) : Json()}, {"detailedBalance", a.detailed_balance}, {"cycles", cycles},
      {"views", views},
      {"dwell", Json::object({{"tau", numbers(tau, 4)}, {"density", numbers(density, 6)}, {"exponential", numbers(exponential, 6)}, {"mean", round(mean, 5)},
                              {"cv2", round(dwell.moment(2) / (mean * mean) - 1, 6)}})},
      {"interventions", interventions}, {"gain", gain},
      {"evidence", Json::object({{"hours", hours}, {"events", static_cast<double>(rec.changes.size())}, {"worlds", numbers(etas, 2)}, {"posterior", numbers(post.probability, 6)},
                                 {"coupled", round(post.of([](double e) { return e > 0; }), 6)}, {"strong", round(post.of([](double e) { return e >= 1; }), 6)}})}});
}

// ---- a building ------------------------------------------------------------------------------------------------

inline Json op_building(const Json& r) {
  as::BuildingSpec spec;
  spec.origin = {num(r, "x", 0), num(r, "y", 0)};
  spec.heading = num(r, "heading", 0);
  spec.width = num(r, "width", 20), spec.depth = num(r, "depth", 12), spec.storey = num(r, "storey", 3.5);
  spec.nx = static_cast<int>(num(r, "nx", 4)), spec.ny = static_cast<int>(num(r, "ny", 2)), spec.floors = static_cast<int>(num(r, "floors", 3));
  auto B = as::building(spec);
  const double air = num(r, "air", 18);
  for (auto& [room, q] : B.heat) q = B.capacity.at(room) * air;
  // The entrance: the first ground-floor facade on the south side (y = 0).
  for (const auto& c : std::vector<mt::RegularCell>(B.complex.cells))
    if (c.type == "facade" && B.corner.at(c.id)[2] == 0 && B.corner.at(c.id)[1] == 0 && B.span.at(c.id)[0]) { as::relabel_face(B, c.id, "entrance"); break; }
  Json events = Json::array();
  // Events at hours: {hour, kind: "stairs" | "doors" | "merge"}.
  std::vector<std::pair<double, std::string>> plan;
  if (r["events"].is_array())
    for (const auto& e : r["events"].as_array()) plan.push_back({e["hour"].as_number(), e["kind"].text()});
  std::sort(plan.begin(), plan.end());
  std::vector<double> sunAlt, sunAz, airT;
  for (const auto& x : r["altitude"].as_array()) sunAlt.push_back(x.as_number());
  for (const auto& x : r["azimuth"].as_array()) sunAz.push_back(x.as_number());
  if (r["airByHour"].is_array())
    for (const auto& x : r["airByHour"].as_array()) airT.push_back(x.as_number());
  if (sunAlt.size() != 24 || sunAz.size() != 24) throw std::invalid_argument("the Sun's altitude and azimuth in each hour (24 values each)");
  auto K = ae::conduction(B);
  const auto spectrum = K.H.values;
  Json hours = Json::array();
  std::size_t next = 0;
  std::uint64_t merged = 0;
  std::vector<std::array<int, 3>> merged_first;
  for (int h = 0; h < 24; h++) {
    while (next < plan.size() && plan[next].first <= h) {
      const auto kind = plan[next].second;
      const double before = B.total_heat();
      Json done = Json::array();
      if (kind == "stairs") {
        // A stair in each room column at x = 0, y = 0, between every two floors.
        for (const auto& c : std::vector<mt::RegularCell>(B.complex.cells))
          if (c.type == "slab" && B.corner.at(c.id)[0] == 0 && B.corner.at(c.id)[1] == 0) done.push(step_json(as::relabel_face(B, c.id, "stair")));
      } else if (kind == "doors") {
        for (const auto& c : std::vector<mt::RegularCell>(B.complex.cells))
          if (c.type == "wall") done.push(step_json(as::relabel_face(B, c.id, "door")));
      } else if (kind == "merge" && !merged) {
        for (const auto& c : std::vector<mt::RegularCell>(B.complex.cells))
          if ((c.type == "wall" || c.type == "door") && B.rooms_of(c.id).size() == 2 && B.corner.at(c.id)[2] == 0) {
            const auto rooms = B.rooms_of(c.id);
            merged_first = B.units.at(rooms[0]);
            const auto st = as::merge_rooms(B, c.id);
            merged = st.created[0];
            done.push(step_json(st));
            break;
          }
      } else if (kind == "partition" && merged) {
        done.push(step_json(as::partition_room(B, merged, merged_first)));
        merged = 0;
      } else throw std::invalid_argument("a building's event is stairs, doors, merge or partition");
      K = ae::conduction(B);
      events.push(Json::object({{"hour", plan[next].first}, {"kind", kind}, {"steps", done}, {"heatChange", round(B.total_heat() - before, 6)},
                                {"rooms", static_cast<double>(B.capacity.size())}, {"regular", mt::admissible(B.complex).regular},
                                {"reachable", static_cast<double>(as::reachable_rooms(B).size())}}));
      next++;
    }
    // The hour's room temperatures, then the hour advanced.
    Json rooms = Json::array();
    for (const auto& [room, q] : B.heat) {
      const auto& u = B.units.at(room);
      Json cubes = Json::array();
      for (const auto& c : u) cubes.push(Json::array({static_cast<double>(c[0]), static_cast<double>(c[1]), static_cast<double>(c[2])}));
      rooms.push(Json::object({{"id", id_text(room)}, {"units", cubes}, {"temperature", round(B.temperature(room), 3)}}));
    }
    const auto gains = ae::solar_gains(B, sunAlt[h], sunAz[h]);
    double sun = 0;
    for (const auto& [room, g] : gains) sun += g;
    hours.push(Json::object({{"hour", static_cast<double>(h)}, {"rooms", rooms}, {"sun", round(sun, 1)}, {"heat", round(B.total_heat(), 1)}}));
    ae::advance_heat(B, K, 3600, airT.size() == 24 ? airT[h] : air, gains);
  }
  Json faces = Json::array();
  for (const auto& c : B.complex.cells)
    if (c.dim == 2 && c.type != "wall" && c.type != "slab" && c.type != "facade" && c.type != "roof" && c.type != "ground")
      faces.push(Json::object({{"id", id_text(c.id)}, {"type", c.type}, {"corner", Json::array({static_cast<double>(B.corner.at(c.id)[0]), static_cast<double>(B.corner.at(c.id)[1]), static_cast<double>(B.corner.at(c.id)[2])})},
                               {"span", Json::array({static_cast<double>(B.span.at(c.id)[0]), static_cast<double>(B.span.at(c.id)[1]), static_cast<double>(B.span.at(c.id)[2])})}}));
  std::vector<double> constants;
  for (std::size_t i = 0; i < std::min<std::size_t>(spectrum.size(), 12); i++) constants.push_back(1 / spectrum[i] / 3600);
  return Json::object({{"spec", Json::object({{"x", spec.origin[0]}, {"y", spec.origin[1]}, {"heading", spec.heading}, {"width", spec.width}, {"depth", spec.depth},
                                              {"storey", spec.storey}, {"nx", spec.nx}, {"ny", spec.ny}, {"floors", spec.floors}})},
                       {"counts", counts(B.complex.counts())}, {"betti", counts(morphe::algebra::betti(mt::chain_complex(B.complex)))},
                       {"timeConstants", numbers(constants, 3)}, {"harmonic", static_cast<double>(K.H.harmonic())}, {"hours", hours}, {"events", events}, {"openings", faces}});
}

// ---- scales and memory ---------------------------------------------------------------------------------------

// The district's scales. Its Hodge Laplacians with the walking stars — a junction weight 1, a segment its
// walking rate (speed over length, per hour), a block 1 — so that −Δ₀ is the generator of one walker
// hopping along the streets: its spectrum the rates at which the district's differences of crowding relax,
// its slowest mode the district's natural division; the harmonic 1-cochains, the flows round its holes.
// Memory: the walkers at the junctions kept and those inside the shops eliminated exactly, each shop opening
// on its nearest junction; a passer-by enters a shop at `visit` per hour, stays its mean `stay` hours, and
// comes out where it went in. Shops on one junction staying alike are one compartment, exactly (their sum
// obeys the same linear law). The particles and the flow: walkers hopping among three junctions as a
// reaction network, its Doi–Peliti operator against the master equation, and the first moment against the
// rate equation (exact for hops, which are of the first order).
inline Json op_scales(const Json& r) {
  const auto D = district_of(r);
  const auto C = mt::chain_complex(D.complex);
  const double speed = num(r, "speed", 1.3) * 3600;  // metres an hour
  const auto J = D.of(0), E = D.of(1);
  std::map<std::uint64_t, std::size_t> jx;
  for (std::size_t i = 0; i < J.size(); i++) jx[J[i]] = i;
  std::vector<double> conductance;
  for (const auto e : E) conductance.push_back(speed / std::max(1.0, D.length(e)));
  const std::vector<std::vector<double>> stars{std::vector<double>(J.size(), 1.0), conductance, {}};
  Json out = Json::object(), spectra = Json::array();
  for (std::size_t p = 0; p < C.n.size(); p++) {
    const auto H = physis::discretisation::hodge_laplacian(C, p, stars);
    std::vector<double> head(H.values.begin(), H.values.begin() + std::min<std::size_t>(H.values.size(), 24));
    spectra.push(Json::object({{"degree", static_cast<double>(p)}, {"values", numbers(head, 6)}, {"harmonic", static_cast<double>(H.harmonic())}, {"cells", static_cast<double>(C.n[p])}}));
    if (p == 0 && H.values.size() > H.harmonic()) {
      // The slowest mode of the district proper: the first non-harmonic eigenvector carried mostly by the largest
      // piece of the network (the box may cut off small pieces, each with modes of its own).
      std::vector<std::size_t> root(J.size());
      for (std::size_t i = 0; i < J.size(); i++) root[i] = i;
      const auto find = [&](std::size_t i) { while (root[i] != i) i = root[i] = root[root[i]]; return i; };
      for (const auto e : E) {
        const auto ends = D.complex.ends(D.at(e));
        root[find(jx.at(D.complex.cells[ends[0]].id))] = find(jx.at(D.complex.cells[ends[1]].id));
      }
      std::map<std::size_t, std::size_t> size;
      for (std::size_t i = 0; i < J.size(); i++) size[find(i)]++;
      const auto largest = std::max_element(size.begin(), size.end(), [](const auto& a, const auto& b) { return a.second < b.second; })->first;
      std::size_t pick = H.harmonic();
      for (std::size_t q = H.harmonic(); q < H.values.size(); q++) {
        double mass = 0;
        for (std::size_t i = 0; i < J.size(); i++)
          if (find(i) == largest) mass += H.vectors(i, q) * H.vectors(i, q);
        if (mass > 0.5) { pick = q; break; }
      }
      Json m = Json::object();
      for (std::size_t i = 0; i < J.size(); i++) m.set(id_text(J[i]), round(find(i) == largest ? H.vectors(i, pick) : 0.0, 6));
      out.set("slowest", Json::object({{"rate", round(H.values[pick], 6)}, {"mode", m}, {"pieces", static_cast<double>(size.size())},
                                       {"largest", static_cast<double>(size.at(largest))}, {"junctions", static_cast<double>(J.size())}}));
    }
    if (p == 1) {
      // The harmonic flows, as fluxes G φ along the segments (divergence-free at every junction, circulation-free
      // round every block): one for each hole.
      Json modes = Json::array();
      for (std::size_t k = 0; k < H.harmonic(); k++) {
        Json m = Json::object();
        double top = 0;
        for (std::size_t i = 0; i < E.size(); i++) top = std::max(top, std::fabs(conductance[i] * H.vectors(i, k)));
        for (std::size_t i = 0; i < E.size(); i++) m.set(id_text(E[i]), round(top > 0 ? conductance[i] * H.vectors(i, k) / top : 0, 5));
        modes.push(Json::object({{"flow", m}}));
      }
      out.set("harmonic", modes);
    }
  }
  out.set("spectra", spectra);
  out.set("betti", counts(morphe::algebra::betti(C)));

  // ---- memory ----
  const double visit = num(r, "visit", 0.6);
  struct Group { std::size_t junction = 0; double stay = 0, rate = 0; std::size_t shops = 0; };
  std::map<std::pair<std::size_t, long long>, Group> groups;
  std::vector<std::size_t> shops_at(J.size(), 0);
  if (r["shops"].is_array())
    for (const auto& x : r["shops"].as_array()) {
      const ul::Point at{x["at"][0].as_number(), x["at"][1].as_number()};
      const double stay = num(x, "stay", 0.5);
      if (!(stay > 0)) throw std::invalid_argument("a shop's stay is > 0 hours");
      std::size_t best = 0;
      double d = std::numeric_limits<double>::infinity();
      for (std::size_t i = 0; i < J.size(); i++) {
        const auto& q = D.point.at(J[i]);
        const double e = std::hypot(q[0] - at[0], q[1] - at[1]);
        if (e < d) d = e, best = i;
      }
      auto& g = groups[{best, std::llround(stay * 1000)}];
      g.junction = best, g.stay = stay, g.rate += visit, g.shops++;
      shops_at[best]++;
    }
  const std::size_t m = J.size(), k = groups.size(), n = m + k;
  // The hops between junctions, at the walking rates.
  std::map<std::pair<std::size_t, std::size_t>, double> hop;
  for (std::size_t q = 0; q < E.size(); q++) {
    const auto ends = D.complex.ends(D.at(E[q]));
    const std::size_t a = jx.at(D.complex.cells[ends[0]].id), b = jx.at(D.complex.cells[ends[1]].id);
    if (a != b) hop[{a, b}] += conductance[q], hop[{b, a}] += conductance[q];
  }
  // Walkers at the junctions and in the shops' compartments obey a master equation in detailed balance (a visitor
  // enters at `rate` and leaves at 1/stay; π = 1 at a junction, rate · stay in a compartment). In the potential
  // u = p/π it is Π u̇ = −A u with A = −LΠ symmetric positive semidefinite, which calx's θ-method propagates
  // (Crank–Nicolson); the streets alone, the shops returning their visitors at once, are the walking alone.
  std::vector<double> pull(m, 0.0), delay(m, 0.0), weight(n, 1.0);
  const auto laplacian = [&](bool with_shops) {
    const std::size_t size = with_shops ? n : m;
    std::vector<std::uint32_t> ri, ci;
    std::vector<double> av, diagonal(size, 0.0);
    for (const auto& [ab, g] : hop) ri.push_back(static_cast<std::uint32_t>(ab.first)), ci.push_back(static_cast<std::uint32_t>(ab.second)), av.push_back(-g), diagonal[ab.first] += g;
    if (with_shops) {
      std::size_t c = m;
      for (const auto& [key, G] : groups) {
        for (const auto [i, j] : {std::pair{G.junction, c}, std::pair{c, G.junction}}) ri.push_back(static_cast<std::uint32_t>(i)), ci.push_back(static_cast<std::uint32_t>(j)), av.push_back(-G.rate);
        diagonal[G.junction] += G.rate, diagonal[c] += G.rate;
        c++;
      }
    }
    for (std::size_t i = 0; i < size; i++) ri.push_back(static_cast<std::uint32_t>(i)), ci.push_back(static_cast<std::uint32_t>(i)), av.push_back(diagonal[i]);
    std::vector<std::uint64_t> order;
    auto A = calx::linalg::compress<std::uint32_t>(size, size, ri, ci, &order);
    std::vector<double> values(order.size());
    for (std::size_t q = 0; q < order.size(); q++) values[q] = av[order[q]];
    return std::pair{std::move(A), std::move(values)};
  };
  {
    std::size_t c = m;
    for (const auto& [key, G] : groups) weight[c++] = G.rate * G.stay, pull[G.junction] += G.rate, delay[G.junction] += G.rate * G.stay;
  }
  for (std::size_t i = 0; i < m; i++) delay[i] = pull[i] > 0 ? delay[i] / pull[i] : 0;
  // The release: walkers set down at one junction (the one asked for, else the one with most shops).
  std::size_t from = 0;
  for (std::size_t i = 0; i < m; i++)
    if (shops_at[i] > shops_at[from]) from = i;
  if (r["from"].is_string() || r["from"].is_number()) {
    const std::uint64_t want = r["from"].is_number() ? static_cast<std::uint64_t>(r["from"].as_number()) : std::stoull(r["from"].text());
    if (jx.count(want)) from = jx.at(want);
  }
  const double released = 100, h = 0.0025;
  const auto [Af, af] = laplacian(true);
  const auto [Aw, aw] = laplacian(false);
  std::vector<double> u(n, 0.0), w(m, 0.0), zero_n(n, 0.0), zero_m(m, 0.0), ones_m(m, 1.0);
  u[from] = released, w[from] = released;
  std::vector<double> tau, onStreet, onStreetMarkov, here, hereMarkov;
  for (int step = 0; step <= 1200; step++) {
    if (step % 20 == 0) {
      double s1 = 0, s2 = 0;
      for (std::size_t i = 0; i < m; i++) s1 += u[i], s2 += w[i];
      tau.push_back(step * h), onStreet.push_back(s1), onStreetMarkov.push_back(s2), here.push_back(u[from]), hereMarkov.push_back(w[from]);
    }
    if (step == 1200) break;
    u = calx::integrators::theta_step(Af, std::span<const double>(af), std::span<const double>(weight), std::span<const double>(u), std::span<const double>(zero_n), h, 0.5).y;
    w = calx::integrators::theta_step(Aw, std::span<const double>(aw), std::span<const double>(ones_m), std::span<const double>(w), std::span<const double>(zero_m), h, 0.5).y;
  }
  // The memory at the release junction: its shops eliminated exactly (physis elimination), K(τ) = B e^{Dτ} C.
  std::vector<double> kernel;
  {
    std::vector<const Group*> mine;
    for (const auto& [key, G] : groups)
      if (G.junction == from) mine.push_back(&G);
    physis::dynamics::LinearElimination L;
    const std::size_t kf = std::max<std::size_t>(mine.size(), 1);
    L.A = calx::linalg::Matrix(1), L.D = calx::linalg::Matrix(kf), L.B = calx::linalg::Rectangular(1, kf), L.C = calx::linalg::Rectangular(kf, 1);
    if (mine.empty()) L.D(0, 0) = -1;
    for (std::size_t g = 0; g < mine.size(); g++) L.B(0, g) = 1 / mine[g]->stay, L.C(g, 0) = mine[g]->rate, L.D(g, g) = -1 / mine[g]->stay;
    for (double t : tau) kernel.push_back(L.kernel(t)(0, 0));
  }
  Json junctions = Json::array();
  for (std::size_t i = 0; i < m; i++)
    junctions.push(Json::object({{"id", id_text(J[i])}, {"shops", static_cast<double>(shops_at[i])}, {"pull", round(pull[i], 4)}, {"delay", round(delay[i], 4)}}));
  double dip = 0, when = 0;
  for (std::size_t q = 0; q < tau.size(); q++)
    if (released - onStreet[q] > dip) dip = released - onStreet[q], when = tau[q];
  double kept = 0;
  for (std::size_t i = 0; i < n; i++) kept += weight[i] * u[i];
  out.set("memory", Json::object({{"tau", numbers(tau, 4)}, {"junctions", junctions}, {"from", id_text(J[from])}, {"released", released},
                                  {"onStreet", numbers(onStreet, 5)}, {"onStreetMarkov", numbers(onStreetMarkov, 5)}, {"here", numbers(here, 5)},
                                  {"hereMarkov", numbers(hereMarkov, 5)}, {"kernel", numbers(kernel, 5)}, {"inside", round(dip, 4)}, {"insideAt", round(when, 3)},
                                  {"kept", round(kept - released, 9)}, {"compartments", static_cast<double>(k)}, {"visit", visit}}));

  // ---- the particles and the flow: three junctions joined in the district, hops at their walking rates ----
  {
    std::vector<std::size_t> three{from};
    for (std::size_t round_ = 0; round_ < 2 && three.size() < 3; round_++)
      for (std::size_t q = 0; q < E.size() && three.size() < 3; q++) {
        const auto ends = D.complex.ends(D.at(E[q]));
        const std::size_t a = jx.at(D.complex.cells[ends[0]].id), b = jx.at(D.complex.cells[ends[1]].id);
        const bool ha = std::find(three.begin(), three.end(), a) != three.end(), hb = std::find(three.begin(), three.end(), b) != three.end();
        if (ha != hb) three.push_back(ha ? b : a);
      }
    namespace pc = physis::chemistry;
    pc::Network net;
    for (std::size_t i = 0; i < three.size(); i++) net.species.push_back("X" + std::to_string(i));
    for (std::size_t i = 0; i < three.size(); i++)
      for (std::size_t j = 0; j < three.size(); j++) {
        const auto it = hop.find({three[i], three[j]});
        if (i == j || it == hop.end()) continue;
        const long long hundredths = std::max<long long>(1, std::llround(100 * it->second));
        net.add("X" + std::to_string(i) + "→X" + std::to_string(j), {{"X" + std::to_string(i), 1}}, {{"X" + std::to_string(j), 1}}, pc::Rational(hundredths, 100));
      }
    const std::size_t cap = 4;
    const pc::Fock Fk{three.size(), cap};
    const auto Lop = pc::network_operator(Fk, net);
    pc::State n0(three.size(), 0);
    n0[0] = static_cast<long long>(cap);
    const auto S = pc::reachable(net, n0);
    const auto Q = pc::generator<double>(net, S);
    const auto idx = [&](const pc::State& x) { return Fk.index(x); };
    double worst = 0;
    for (std::size_t x = 0; x < S.states.size(); x++)
      for (std::size_t y = 0; y < S.states.size(); y++) {
        if (y == x) continue;
        const auto it = Q[x].find(y);
        worst = std::max(worst, std::fabs(Lop(idx(S.states[y]), idx(S.states[x])) - (it == Q[x].end() ? 0.0 : it->second)));
      }
    Json names = Json::array();
    for (auto i : three) names.push(id_text(J[i]));
    out.set("operators", Json::object({{"junctions", names}, {"walkers", static_cast<double>(cap)}, {"states", static_cast<double>(S.states.size())},
                                       {"fock", static_cast<double>(Fk.size())}, {"deviation", round(worst, 15)}, {"hops", static_cast<double>(net.reactions.size())}}));
  }
  return out;
}

}  // namespace ld
