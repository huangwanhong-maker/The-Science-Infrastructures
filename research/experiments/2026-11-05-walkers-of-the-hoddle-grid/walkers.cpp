// The walkers of the Hoddle Grid (README): urban science's walkers as agents (domains/urban-science/movement) on the
// district of the shared problem "A district's economy", whose trips the flow study calibrated to the counters each way.
//
//   week      the walkers of the mean week: their mean against the flow study's logit flow (exact in expectation, by the
//             hour of departure), their counts against the counters, their time and their frequencies;
//   closure   Little Bourke Street between Swanston and Russell Streets closed to the walkers passing through, as the
//             City's permit for the Lunar New Year of 2017 closed it (27–29 January): how they redistribute;
//   festival  the festival's crowd on Sunday 29 January 2017, fitted to the counters' excess over the Sundays about it,
//             arriving, dwelling and leaving: densities and queues;
//   regime    the street's regime decided online (POMCP on the walkers) against the regimes held all day.
//
// Everything is written to the results directory; the geometry (an extract of OpenStreetMap, ODbL) to maps.json, which
// is not committed.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <tyche/decision/online.hpp>
#include <urban/events/record.hpp>
#include <urban/movement/movement.hpp>
#include <urban/rhythm/record.hpp>
#include <urban/rhythm/tide.hpp>

#include "economy.hpp"
#include "flow.hpp"

namespace {

using economy::Json;
namespace um = urban::movement;
namespace hv = heuresis::variation;
namespace ur = urban::rhythm;
using Clock = std::chrono::steady_clock;
constexpr std::size_t H = economy::week;
constexpr std::size_t none = um::none;

std::string read_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read " + path);
  std::stringstream b;
  b << in.rdbuf();
  return b.str();
}
void write_file(const std::string& path, const std::string& text) {
  std::ofstream f(path, std::ios::binary);
  f << text;
  if (!f) throw std::runtime_error("cannot write " + path);
}
double rnd(double x, double d) { return std::isfinite(x) ? std::round(x * d) / d : std::nan(""); }
Json jr(double x, double d) { return std::isfinite(x) ? Json(rnd(x, d)) : Json(); }
template <class V>
Json arr(const V& v, double d) {
  Json a = Json::array();
  for (const auto x : v) a.push(jr(static_cast<double>(x), d));
  return a;
}
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

// ---- the flow study's network as the walkers' -------------------------------------------------------------------------

struct Bridge {
  std::shared_ptr<um::Walkway> W;
  std::vector<std::size_t> mine, theirs;   // their edge → mine (none for a candidate link's), mine → theirs
  std::size_t state(std::size_t their_state) const { return 2 * mine.at(their_state / 2) + their_state % 2; }
};

// The network as built (no candidate link made): the flow study's segments, their widths and parent lengths kept, so
// that the walkers' law is its law exactly. Signals: at every junction of three carriageways or more (an assumption:
// the data library holds no record of the signals; the Hoddle Grid's road intersections are signalled).
Bridge bridge(const flow::Model& M) {
  Bridge B;
  B.W = std::make_shared<um::Walkway>();
  auto& W = *B.W;
  for (const auto& p : M.net.at) W.at.push_back({p[0], p[1]});
  B.mine.assign(M.net.edges.size(), none);
  for (std::size_t e = 0; e < M.net.edges.size(); e++) {
    const auto& x = M.net.edges[e];
    if (x.link >= 0) continue;
    um::Segment s;
    s.u = x.u, s.v = x.v, s.length = x.length, s.width = x.width, s.parent = x.parent, s.kind = x.kind, s.name = x.name;
    for (const auto& p : x.line) s.line.push_back({p[0], p[1]});
    B.mine[e] = W.segments.size();
    B.theirs.push_back(e);
    W.segments.push_back(std::move(s));
  }
  W.index();
  for (std::size_t j = 0; j < W.junctions(); j++) W.signalled[j] = W.carriageways[j] >= 3;
  return B;
}

double along_of(const um::Walkway& W, std::size_t e, const um::Point& p) {
  const auto& l = W.segments[e].line;
  double best = 1e300, at = 0, run = 0;
  for (std::size_t k = 0; k + 1 < l.size(); k++) {
    const double ex = l[k + 1][0] - l[k][0], ey = l[k + 1][1] - l[k][1], L2 = ex * ex + ey * ey, L = std::sqrt(L2);
    const double t = L2 > 0 ? std::clamp(((p[0] - l[k][0]) * ex + (p[1] - l[k][1]) * ey) / L2, 0.0, 1.0) : 0.0;
    const double d = std::hypot(p[0] - l[k][0] - t * ex, p[1] - l[k][1] - t * ey);
    if (d < best) best = d, at = run + t * L;
    run += L;
  }
  return at;
}

// The counters of the flow study's calibration as the walkers' counting lines, with their fitted handedness.
std::vector<um::Counter> counters_of(const flow::Model& M, const Bridge& B) {
  std::vector<um::Counter> out;
  for (const auto& c : M.counters) {
    um::Counter k;
    k.name = c.name;
    k.edge = B.mine.at(c.edge);
    k.along = along_of(*B.W, k.edge, {c.at[0], c.at[1]});
    k.share = c.share;
    if (c.directed) k.hand = {c.hand, -c.hand}, k.way = {B.state(c.state[0]), B.state(c.state[1])};
    else k.hand = {0, 0}, k.way = {2 * k.edge, 2 * k.edge + 1};
    out.push_back(k);
  }
  return out;
}

// ---- the flow study's trips as streams --------------------------------------------------------------------------------

// The patterns of the flow study (its README §2.3) over the week's hours [h0, h0 + hours): each a fixed law of origins
// and targets times an hourly rate, so that the walkers' mean is the flow's mean_flows at every hour.
std::vector<um::Stream> streams_of(const flow::Model& M, const economy::Problem& P, um::Targets& T, const um::Walkway& W, std::size_t h0, std::size_t hours) {
  std::vector<um::Stream> out;
  const auto week = [&](std::size_t k) { return (h0 + k) % H; };
  const auto draw = M.frontage_draw(P.as_built());
  const std::size_t E = M.entrances.size();
  const auto junction = [&](std::size_t j) { return T.junction(W, j); };
  for (const std::size_t m : {flow::arrive_work, flow::arrive_home, flow::leave_work, flow::leave_home}) {
    double any = 0;
    for (double v : M.x[m]) any += v;
    if (!(any > 0)) continue;
    const auto& mass = (m == flow::arrive_work || m == flow::leave_work) ? M.work : M.home;
    double total = 0;
    for (double v : mass) total += v;
    const bool arrive = m == flow::arrive_work || m == flow::arrive_home;
    for (std::size_t i = 0; i < E; i++) {
      um::Stream s;
      s.name = std::string(flow::mode_name(m)) + ", " + M.entrances[i].street;
      s.mode = static_cast<int>(m);
      for (const std::size_t b : M.building_nodes) {
        if (!(mass[b] > 0) || b == M.entrances[i].node) continue;
        s.origin.push_back(arrive ? M.entrances[i].node : b);
        s.target.push_back(junction(arrive ? b : M.entrances[i].node));
        s.weight.push_back(mass[b] / total);
      }
      // The pair whose origin is its target makes no trip (the flow study's trip(o, o) is zero): the law keeps the
      // others' shares, so the rate is the amplitude times the gross times their sum.
      double kept = 0;
      for (double w : s.weight) kept += w;
      s.rate.resize(hours);
      for (std::size_t k = 0; k < hours; k++) s.rate[k] = M.x[m][week(k)] * M.entrances[i].gross[week(k)] * kept;
      out.push_back(std::move(s));
    }
  }
  for (std::size_t i = 0; i < E; i++)
    for (std::size_t j = 0; j < E; j++) {
      if (i == j) continue;
      um::Stream s;
      s.name = "through, " + M.entrances[i].street + " to " + M.entrances[j].street;
      s.mode = flow::through;
      s.origin = {M.entrances[i].node}, s.target = {junction(M.entrances[j].node)}, s.weight = {1};
      s.rate.resize(hours);
      for (std::size_t k = 0; k < hours; k++) {
        const std::size_t h = week(k);
        s.rate[k] = M.gsum[h] > 0 ? M.x[flow::through][h] * M.entrances[i].gross[h] * M.entrances[j].gross[h] / M.gsum[h] : 0.0;
      }
      out.push_back(std::move(s));
    }
  {
    um::Stream s;
    s.name = "local", s.mode = flow::local;
    double norm = 0;
    for (const std::size_t a : M.building_nodes)
      for (const std::size_t b : M.building_nodes) {
        if (a == b) continue;
        const double t = (M.work[a] + M.home[a]) * (M.work[b] + M.home[b]) * std::exp(-std::hypot(M.net.at[a][0] - M.net.at[b][0], M.net.at[a][1] - M.net.at[b][1]) / M.o.local_length);
        s.origin.push_back(a), s.target.push_back(junction(b)), s.weight.push_back(t), norm += t;
      }
    for (double& w : s.weight) w /= norm;
    s.rate.resize(hours);
    for (std::size_t k = 0; k < hours; k++) s.rate[k] = M.x[flow::local][week(k)];
    out.push_back(std::move(s));
  }
  // Visits: to each open attractor from the origins π and back, per unit of weight.
  const auto visits = [&](const std::vector<std::pair<std::size_t, double>>& nodes, um::Stream& s) {
    double total = 0;
    for (const auto& [j, a] : nodes) {
      if (!(a > 0)) continue;
      for (const std::size_t o : M.origin_nodes) {
        if (o == j) continue;
        const double w = M.origin[o] * a;
        s.origin.push_back(o), s.target.push_back(junction(j)), s.weight.push_back(w);
        s.origin.push_back(j), s.target.push_back(junction(o)), s.weight.push_back(w);
        total += 2 * w;
      }
    }
    return total;
  };
  for (std::size_t q = 0; q < economy::kinds; q++) {
    if (!(M.o.attraction[q] > 0)) continue;
    um::Stream s;
    s.name = std::string("visit, ") + economy::kind_name(static_cast<std::int64_t>(q));
    s.mode = flow::visit;
    std::vector<std::pair<std::size_t, double>> nodes;
    for (const std::size_t j : M.attractor_nodes) nodes.push_back({j, M.attract[q][j] * M.o.attraction[q]});
    const double total = visits(nodes, s);
    if (!(total > 0)) continue;
    s.rate.resize(hours);
    for (std::size_t k = 0; k < hours; k++) s.rate[k] = M.open_kind[q][week(k)] ? M.x[flow::visit][week(k)] / M.W_ref * total : 0.0;
    out.push_back(std::move(s));
  }
  for (std::size_t g = 0; g < P.size(); g++) {
    um::Stream s;
    s.name = "visit, frontage " + std::to_string(g);
    s.mode = flow::visit;
    std::vector<std::pair<std::size_t, double>> nodes;
    for (const auto& [j, wg] : M.frontage_nodes[g]) nodes.push_back({j, wg});
    const double total = visits(nodes, s);
    if (!(total > 0)) continue;
    s.rate.resize(hours);
    for (std::size_t k = 0; k < hours; k++) s.rate[k] = M.x[flow::visit][week(k)] / M.W_ref * draw[g * H + week(k)] * total;
    out.push_back(std::move(s));
  }
  return out;
}

// The trips of an hour of the run: every pair at its share of its stream's rate (of one route class, or all: −1).
std::vector<um::Trip> trips_at(const std::vector<um::Stream>& streams, std::size_t k, int cls = -1) {
  std::vector<um::Trip> t;
  for (const auto& s : streams) {
    if (cls >= 0 && s.route_class != cls) continue;
    if (k >= s.rate.size() || !(s.rate[k] > 0)) continue;
    double total = 0;
    for (double w : s.weight) total += w;
    for (std::size_t p = 0; p < s.origin.size(); p++)
      if (s.weight[p] > 0) t.push_back({s.origin[p], s.target[p], s.rate[k] * s.weight[p] / total});
  }
  return t;
}

// ---- the week's lines ------------------------------------------------------------------------------------------------

Json lines_of(const std::array<double, H>& x) {
  Json o = Json::object();
  for (const auto& [name, m] : std::vector<std::pair<const char*, std::size_t>>{{"W1", 1}, {"D1-W1", 6}, {"D1", 7}, {"D1+W1", 8}, {"D2", 14}, {"D3", 21}}) {
    const auto c = flow::line(std::span<const double>(x.data(), H), m);
    o.set(name, Json::object({{"amplitude", jr(flow::amplitude(c), 10)}, {"crest_h", jr(flow::crest_hour(c, m), 100)}}));
  }
  double mean = 0;
  for (double v : x) mean += v / H;
  o.set("mean", jr(mean, 10));
  return o;
}

double r2_of(const std::vector<double>& y, const std::vector<double>& p) {
  double m = 0;
  for (double v : y) m += v / static_cast<double>(y.size());
  double ss = 0, st = 0;
  for (std::size_t i = 0; i < y.size(); i++) ss += (y[i] - p[i]) * (y[i] - p[i]), st += (y[i] - m) * (y[i] - m);
  return 1 - ss / st;
}

// ---- the festival block ------------------------------------------------------------------------------------------------

// Little Bourke Street between Swanston Street and Russell Street: the path of its segments between the junctions it
// shares with the two.
um::Street festival_block(const um::Walkway& W) {
  const auto meets = [&](std::size_t x, const std::string& name) {
    for (const std::size_t s : W.leaving[x])
      if (W.segments[s / 2].name == name) return true;
    return false;
  };
  std::size_t west = none, east = none;
  for (std::size_t x = 0; x < W.junctions(); x++)
    if (meets(x, "Little Bourke Street")) {
      if (meets(x, "Swanston Street")) west = x;
      if (meets(x, "Russell Street")) east = x;
    }
  if (west == none || east == none) throw std::runtime_error("the festival block: Little Bourke Street does not meet Swanston and Russell Streets");
  std::vector<std::size_t> prev(W.junctions(), none), via(W.junctions(), none);
  std::vector<std::size_t> frontier{west};
  prev[west] = west;
  for (std::size_t q = 0; q < frontier.size(); q++) {
    const std::size_t x = frontier[q];
    for (const std::size_t s : W.leaving[x]) {
      if (W.segments[s / 2].name != "Little Bourke Street") continue;
      const std::size_t y = W.head(s);
      if (prev[y] != none) continue;
      prev[y] = x, via[y] = s / 2, frontier.push_back(y);
    }
  }
  if (prev[east] == none) throw std::runtime_error("the festival block: no path along Little Bourke Street");
  um::Street S;
  S.name = "Little Bourke Street, from Swanston Street past Russell Street to the edge of the district";
  S.entry = west;
  std::set<std::size_t> on_path{west};
  for (std::size_t y = east; y != west; y = prev[y]) S.edges.push_back(via[y]), on_path.insert(y);
  std::reverse(S.edges.begin(), S.edges.end());
  // On past Russell Street, eastward, to the box's edge: the furthest junction along Little Bourke Street from Russell
  // Street that is reached without passing back along the block.
  std::vector<std::size_t> back(W.junctions(), none), by(W.junctions(), none);
  std::vector<double> dist(W.junctions(), -1);
  std::vector<std::size_t> front{east};
  dist[east] = 0, back[east] = east;
  for (std::size_t q = 0; q < front.size(); q++) {
    const std::size_t x = front[q];
    for (const std::size_t st : W.leaving[x]) {
      if (W.segments[st / 2].name != "Little Bourke Street") continue;
      const std::size_t y = W.head(st);
      if (on_path.count(y) || dist[y] >= 0) continue;
      dist[y] = dist[x] + W.segments[st / 2].length, back[y] = x, by[y] = st / 2, front.push_back(y);
    }
  }
  std::size_t far = east;
  for (std::size_t x = 0; x < W.junctions(); x++)
    if (dist[x] > dist[far]) far = x;
  std::vector<std::size_t> beyond;
  for (std::size_t y = far; y != east; y = back[y]) beyond.push_back(by[y]);
  std::reverse(beyond.begin(), beyond.end());
  S.edges.insert(S.edges.end(), beyond.begin(), beyond.end());
  return S;
}

// Each stream's pairs that start or end on the street, or cannot avoid it (no way between them with the street closed:
// the lanes behind it), apart, as walkers of the street's own class (1): under a regime closed to the walkers passing
// through, they still reach it. The law of the whole is kept (Poisson splitting): each part's rate is the stream's
// times its pairs' share of the weight.
std::vector<um::Stream> with_access(const std::vector<um::Stream>& streams, const um::Street& S, const um::Walkway& W, const um::Targets& T, const um::Values& closed) {
  std::set<std::size_t> on;
  for (const std::size_t e : S.edges) on.insert(W.segments[e].u), on.insert(W.segments[e].v);
  std::vector<um::Stream> out;
  for (const auto& st : streams) {
    um::Stream keep = st, access = st;
    keep.origin.clear(), keep.target.clear(), keep.weight.clear(), keep.back.clear();
    access.origin.clear(), access.target.clear(), access.weight.clear(), access.back.clear();
    access.route_class = 1, access.name = st.name + " (to or from the street)";
    double wk = 0, wa = 0;
    for (std::size_t p = 0; p < st.origin.size(); p++) {
      const std::size_t d = T.junction_of(st.target[p]);
      const bool touches = on.count(st.origin[p]) || (d != none && on.count(d)) || !(closed.routes().total(st.origin[p], closed[st.target[p]]) > 0);
      um::Stream& to = touches ? access : keep;
      to.origin.push_back(st.origin[p]), to.target.push_back(st.target[p]), to.weight.push_back(st.weight[p]);
      if (!st.back.empty()) to.back.push_back(st.back[p]);
      (touches ? wa : wk) += st.weight[p];
    }
    const double total = wk + wa;
    if (wk > 0) {
      for (double& r : keep.rate) r *= wk / total;
      out.push_back(std::move(keep));
    }
    if (wa > 0) {
      for (double& r : access.rate) r *= wa / total;
      out.push_back(std::move(access));
    }
  }
  return out;
}

// The counts of a sensor (the data library's street-counts product) as a record in UTC.
ur::Record counts_record(const std::string& path) {
  const Json j = Json::parse(read_file(path));
  ur::Record r;
  r.id = path;
  r.t0 = ur::parse_utc(j["t0_utc"].text());
  r.step = j["step_s"].as_number();
  for (const auto& o : j["offsets"].as_array()) r.offsets.push_back({ur::parse_utc(o[0].text()), o[1].as_number()});
  for (const auto& v : j["counts"].as_array()) r.y.push_back(v.is_null() ? std::nan("") : v.as_number());
  return r;
}

// The share of a dwell of mean D hours, begun uniformly within an hour, that ends d hours later (d = 0, 1, …).
double back_kernel(double D, int d) {
  double s = 0;
  const int n = 400;
  for (int i = 0; i < n; i++) {
    const double u = (i + 0.5) / n;
    const double a = std::max(0.0, d - u), b = std::max(0.0, d + 1 - u);
    s += std::exp(-a / D) - std::exp(-b / D);
  }
  return s / n;
}

// A run's summary as JSON.
Json totals_json(const um::Totals& t) {
  return Json::object({{"departed", t.departed}, {"arrived", t.arrived}, {"failed", t.failed}, {"person_hours", jr(t.present / 3600, 10)}, {"waiting_hours", jr(t.waiting / 3600, 10)},
                       {"los_e_hours", jr(t.los_e / 3600, 10)}, {"los_f_hours", jr(t.los_f / 3600, 10)}, {"mean_travel_s", jr(t.arrived > 0 ? t.travel / t.arrived : 0, 10)},
                       {"mean_delay_s", jr(t.arrived > 0 ? t.delay / t.arrived : 0, 100)}});
}

Json positions_json(const um::Simulation& sim, const std::vector<double>* local = nullptr) {
  Json a = Json::array();
  for (const auto& p : sim.positions()) {
    Json q = Json::array({jr(p.at[0], 10), jr(p.at[1], 10), p.mode, static_cast<int>(p.status)});
    (void)local;
    a.push(q);
  }
  return a;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 5) {
    std::cerr << "usage: walkers <problem-inputs.json> <problem.json> <sensors.json> <results directory> [--threads N] [--task week|closure|festival|regime|all] [--days N] [--regime-cases N]\n";
    return 2;
  }
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  std::size_t threads = 4, replicates = 2, regime_cases = static_cast<std::size_t>(-1);
  std::string task = "all";
  for (int i = 5; i < argc; i++) {
    if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) threads = static_cast<std::size_t>(std::stoul(argv[++i]));
    else if (std::strcmp(argv[i], "--task") == 0 && i + 1 < argc) task = argv[++i];
    else if (std::strcmp(argv[i], "--days") == 0 && i + 1 < argc) replicates = static_cast<std::size_t>(std::stoul(argv[++i]));
    else if (std::strcmp(argv[i], "--regime-cases") == 0 && i + 1 < argc) regime_cases = static_cast<std::size_t>(std::stoul(argv[++i]));
  }
  // "sunday": the closure, the festival and the regime, without the week.
  const auto wants = [&](const char* t) { return task == "all" || task == t || (task == "sunday" && std::string(t) != "week"); };
  const std::string out = argv[4];
  const auto T0 = Clock::now();
  const Json in = Json::parse(read_file(argv[1]));
  const Json pj = Json::parse(read_file(argv[2]));
  const Json sj = Json::parse(read_file(argv[3]));
  const economy::Problem P(in, &pj);
  flow::Options fo;   // θ = 0.02 per metre, as the flow study chose
  flow::Model M(P, in, sj, fo);
  const hv::Bits none_made(P.links.size(), 0);
  // The flow study's calibration, as it ran it: the homes' modes left out, the visits set at ν = 0.1 of the trips.
  std::vector<std::size_t> cal_edges;
  for (const auto& c : M.counters) cal_edges.push_back(c.edge);
  for (auto e : M.fronted_edges()) cal_edges.push_back(e);
  std::sort(cal_edges.begin(), cal_edges.end());
  cal_edges.erase(std::unique(cal_edges.begin(), cal_edges.end()), cal_edges.end());
  const auto Tr = M.trips(none_made);
  const auto fit = M.calibrate(M.terms(Tr, cal_edges), true, {1, 0, 1, 0, 1, 1, 1});
  std::printf("the flow study's calibration again: gross R2 %.4f (leave-one-counter-out %.4f), net R2 %.4f; %zu counters, %zu entrances  [%.0f s]\n", fit.r2, fit.loo_r2,
              fit.net_r2, M.counters.size(), M.entrances.size(), seconds_since(T0));
  double amplitude_gap = 0;
  {
    const std::string cal = std::string(argv[4]) + "/../../2026-11-05-a-district-economy-as-a-flow/results/calibration.json";
    std::ifstream probe(cal);
    if (probe) {
      const Json C = Json::parse(read_file(cal));
      for (std::size_t m = 0; m < flow::modes; m++) {
        const auto& a = C["amplitudes"][flow::mode_name(m)].as_array();
        for (std::size_t h = 0; h < H; h++) amplitude_gap = std::max(amplitude_gap, std::fabs(a[h].as_number() - M.x[m][h]) / std::max(1e-9, std::fabs(a[h].as_number()) + 1e-6));
      }
      std::printf("  amplitudes against the flow study's results/calibration.json: largest relative difference %.2e (it rounds to 4 significant figures)\n", amplitude_gap);
    }
  }

  const Bridge B = bridge(M);
  const um::Walkway& W = *B.W;
  const std::size_t S = W.states(), E = W.edges();
  const um::RouteLaw law{.theta = fo.theta, .reference_width = fo.widths.reference};
  auto routes = std::make_shared<const um::Routes>(B.W, law);
  calx::scheduler::ThreadPool pool(threads);
  const auto counters = counters_of(M, B);
  const um::Street block = festival_block(W);
  std::printf("the walkers' network: %zu junctions, %zu segments (%zu directed), spectral radius %.4f; %zu signalled junctions; the festival block: %zu segments, %.0f m\n",
              W.junctions(), E, S, routes->spectral_radius(), static_cast<std::size_t>(std::count(W.signalled.begin(), W.signalled.end(), 1)), block.edges.size(),
              [&] { double L = 0; for (auto e : block.edges) L += W.segments[e].length; return L; }());

  // The geometry, for the figures (not committed: ODbL).
  Json maps = Json::object();
  {
    Json segs = Json::array();
    for (std::size_t e = 0; e < E; e++) {
      Json line = Json::array();
      for (const auto& p : W.segments[e].line) line.push(Json::array({jr(p[0], 100), jr(p[1], 100)}));
      segs.push(Json::object({{"edge", e}, {"name", W.segments[e].name}, {"kind", W.segments[e].kind}, {"width", W.segments[e].width}, {"line", line}}));
    }
    maps.set("segments", segs);
    Json cs = Json::array();
    for (std::size_t k = 0; k < counters.size(); k++) {
      const auto p = W.point(2 * counters[k].edge, counters[k].along);
      cs.push(Json::object({{"name", counters[k].name}, {"id", static_cast<double>(M.counters[k].id)}, {"x", jr(p[0], 10)}, {"y", jr(p[1], 10)}}));
    }
    maps.set("counters", cs);
    Json ent = Json::array();
    for (const auto& e : M.entrances) ent.push(Json::array({jr(W.at[e.node][0], 10), jr(W.at[e.node][1], 10)}));
    maps.set("entrances", ent);
    maps.set("block", Json(std::vector<double>(block.edges.begin(), block.edges.end())));
    const auto& box = in["district"]["district_box"];
    maps.set("box", Json::array({box["x0"].as_number(), box["y0"].as_number(), box["x1"].as_number(), box["y1"].as_number()}));
  }

  // ===================================================================================================================
  // 1. The mean week
  // ===================================================================================================================
  um::Targets T;
  auto week_streams = streams_of(M, P, T, W, 0, H);
  // The mean Sunday (its pairs that start or end on the festival's street, or cannot avoid it, apart as the street's own
  // walkers once the regimes' routes are solved, below).
  std::vector<um::Stream> sunday_streams = streams_of(M, P, T, W, 144, 24);
  // The festival's visitors (task 3): from the entrances and the buildings, each to one of the street's segments (by its
  // length), dwelling there, and back.
  std::vector<std::pair<std::size_t, double>> along_street;
  for (const std::size_t e : block.edges) along_street.push_back({T.segments({e}), W.segments[e].length});
  for (const std::size_t o : M.origin_nodes) T.junction(W, o);
  for (const auto& E : M.entrances) T.junction(W, E.node);
  auto Tp = std::make_shared<const um::Targets>(T);
  auto V = std::make_shared<const um::Values>(routes, T, &pool);
  std::printf("%zu streams in the week, %zu targets, their walk sums solved  [%.0f s]\n", week_streams.size(), T.size(), seconds_since(T0));
  // The regimes of the street: [regime][class]; class 1, the street's own walkers.
  auto RR = std::make_shared<const um::RegimeRoutes>(um::regime_routes(B.W, law, block, T, &pool));
  sunday_streams = with_access(sunday_streams, block, W, T, *RR->values[2][0]);
  std::printf("the street's regimes solved; the Sunday's %zu streams, its own walkers apart  [%.0f s]\n", sunday_streams.size(), seconds_since(T0));

  if (task == "diag") {
    // Single trips: the module's visits against the flow study's Trips::visits, and the walk sums solved tighter.
    const std::vector<std::size_t> nodes(M.building_nodes.begin(), M.building_nodes.begin() + std::min<std::size_t>(12, M.building_nodes.size()));
    double worst = 0, worst_tight = 0;
    for (const std::size_t d : nodes) {
      const auto& A = T.states(T.junction(W, d));
      const auto z = routes->value(A), zt = routes->value(A, 1e-16, 2000000);
      for (const std::size_t o : nodes) {
        if (o == d) continue;
        const auto f = routes->visits(o, A, z), ft = routes->visits(o, A, zt);
        double top = 0, g = 0, gt = 0;
        for (std::size_t c = 0; c < S; c++) {
          const std::size_t e = B.theirs[c / 2];
          const double v = Tr.visits(o, d, 2 * e + c % 2);
          top = std::max(top, v), g = std::max(g, std::fabs(f[c] - v)), gt = std::max(gt, std::fabs(ft[c] - v));
        }
        if (g / top > worst) {
          std::size_t wc = 0;
          double gg = 0;
          for (std::size_t c = 0; c < S; c++) {
            const double v = Tr.visits(o, d, 2 * B.theirs[c / 2] + c % 2);
            if (std::fabs(f[c] - v) > gg) gg = std::fabs(f[c] - v), wc = c;
          }
          const std::size_t tc = 2 * B.theirs[wc / 2] + wc % 2;
          std::printf("  pair %zu → %zu: worst state %zu (%s, %zu → %zu), visits %.9e against %.9e; z %.9e σ %.9e; |leaving o| %zu, |A| %zu, o's degree %zu, d's %zu\n", o, d, wc,
                      W.segments[wc / 2].name.c_str(), W.tail(wc), W.head(wc), f[wc], Tr.visits(o, d, tc), z[wc], Tr.sigma(Tr.destination(d), tc), W.leaving[o].size(), A.size(),
                      W.leaving[o].size(), W.leaving[d].size());
        }
        worst = std::max(worst, g / top), worst_tight = std::max(worst_tight, gt / top);
      }
    }
    std::printf("single trips, %zu destinations: the module's visits against the flow study's, largest difference %.3e of the largest visit (%.3e with the walk sums solved to 1e-16)\n",
                nodes.size(), worst, worst_tight);
    // One destination: z against the flow study's σ, and π against its π, state by state.
    {
      const std::size_t d = nodes[0], o = nodes[1];
      const auto& A = T.states(T.junction(W, d));
      const auto z = routes->value(A);
      const auto& D = Tr.destination(d);
      double dz = 0, zt = 0;
      std::size_t wz = 0;
      for (std::size_t c = 0; c < S; c++) {
        const double sg = Tr.sigma(D, 2 * B.theirs[c / 2] + c % 2);
        zt = std::max(zt, sg);
        if (std::fabs(z[c] - sg) > dz) dz = std::fabs(z[c] - sg), wz = c;
      }
      std::printf("z against σ for one destination: largest difference %.3e (of the largest %.3e) at state %zu (%s), z %.6e σ %.6e; |A| %zu against %zu\n", dz, zt, wz,
                  W.segments[wz / 2].name.c_str(), z[wz], Tr.sigma(D, 2 * B.theirs[wz / 2] + wz % 2), A.size(), D.A.size());
      const auto yq = Tr.yq(o, D);
      double Wt = 0;
      for (double v : yq) Wt += v;
      std::printf("  W_od: the module's %.12e, the flow study's %.12e\n", routes->total(o, z), Wt);
      const auto f = routes->visits(o, A, z);
      double g = 0;
      std::size_t wc = 0;
      for (std::size_t c = 0; c < S; c++) {
        const double v = Tr.visits(o, d, 2 * B.theirs[c / 2] + c % 2);
        if (std::fabs(f[c] - v) > g) g = std::fabs(f[c] - v), wc = c;
      }
      const std::size_t tc = 2 * B.theirs[wc / 2] + wc % 2;
      const double yoc = Tr.Y(o, tc);
      std::printf("  worst state %zu (%s, tail %zu head %zu; origin %zu, destination %zu): visits %.9e against %.9e; z %.6e σ %.6e; Y %.6e, their π %.6e; in A %d\n", wc,
                  W.segments[wc / 2].name.c_str(), W.tail(wc), W.head(wc), o, d, f[wc], Tr.visits(o, d, tc), z[wc], Tr.sigma(D, tc), yoc, Tr.pi(D, yq, tc, yoc),
                  static_cast<int>(std::find(A.begin(), A.end(), wc) != A.end()));
      // The same visits by the flow study's formulas on the module's own dense Z.
      calx::linalg::Matrix Am(S);
      for (std::size_t a = 0; a < S; a++) {
        Am(a, a) += 1;
        for (const auto& arc : routes->out(a)) Am(a, arc.to) -= arc.w;
      }
      const auto Zm = calx::linalg::solve(Am, calx::linalg::Matrix::identity(S));
      double Ym = 0;
      for (const std::size_t b : W.leaving[o]) Ym += routes->weight(b) * Zm(b, wc);
      std::printf("  Y from the module's own Z: %.6e (the flow study's Z: %.6e)\n", Ym, yoc);
    }
    // The worst pair seen (15 → 16, a dead end): π by Gauss–Seidel to 1e-16, and by the module's own dense Z.
    {
      const std::size_t o = 15, d = 16;
      const auto& A = T.states(T.junction(W, d));
      const auto z = routes->value(A, 1e-16, 2000000);
      std::vector<char> in_t(S, 0);
      for (auto c : A) in_t[c] = 1;
      std::vector<double> src(S, 0.0);
      const double Wm = routes->total(o, z);
      for (const std::size_t b : W.leaving[o]) src[b] = routes->weight(b) / Wm;
      const auto pi = routes->forward(src, in_t, 1e-16, 5000000);
      double sumA = 0, sumA_default = 0;
      const auto f = routes->visits(o, A, z);
      for (auto c : A) sumA += pi[c] * z[c], sumA_default += f[c];
      std::printf("pair 15 → 16: Σ over the target of the visits %.12f (π to 1e-16), %.12f (default); |A| %zu; leaving 15: %zu; W %.6e\n", sumA, sumA_default, A.size(),
                  W.leaving[o].size(), Wm);
      // Sweeps' history: the change per sweep of the forward iteration.
      std::vector<double> p(src);
      for (int sweep = 0; sweep < 4000; sweep++) {
        double change = 0, top = 0;
        for (std::size_t b = 0; b < S; b++) {
          double v = src[b];
          for (const auto& arc : routes->in(b))
            if (!in_t[arc.to]) v += p[arc.to] * arc.w;
          change = std::max(change, std::fabs(v - p[b])), top = std::max(top, std::fabs(v));
          p[b] = v;
        }
        if (sweep % 400 == 0 || sweep < 5) {
          double s2 = 0;
          for (auto c : A) s2 += p[c];
          std::printf("  sweep %4d: change %.3e, top %.3e, Σ over the target %.12f\n", sweep, change, top, s2);
        }
      }
    }
    // The transitions, entry by entry, and the walk sums Z of the flow study against the module's.
    const auto Mt = flow::transitions(M.net, M.omega, none_made);
    double dm = 0;
    std::size_t where_a = 0, where_b = 0, extra = 0, missing = 0;
    for (std::size_t a = 0; a < S; a++) {
      const std::size_t ta = 2 * B.theirs[a / 2] + a % 2;
      std::vector<double> row(S, 0.0);
      for (const auto& arc : routes->out(a)) row[arc.to] = arc.w;
      for (std::size_t b = 0; b < S; b++) {
        const std::size_t tb = 2 * B.theirs[b / 2] + b % 2;
        const double d = std::fabs(row[b] - Mt(ta, tb));
        if (row[b] > 0 && Mt(ta, tb) == 0) extra++;
        if (row[b] == 0 && Mt(ta, tb) > 0) missing++;
        if (d > dm) dm = d, where_a = a, where_b = b;
      }
    }
    std::printf("transitions: largest |ΔM| %.3e at %zu → %zu (%s → %s); %zu entries only in the module's, %zu only in the flow study's\n", dm, where_a, where_b,
                W.segments[where_a / 2].name.c_str(), W.segments[where_b / 2].name.c_str(), extra, missing);
    // The transitions into link states and out of the network as built (they should carry nothing).
    double into_links = 0;
    for (std::size_t a = 0; a < M.net.states(); a++)
      for (std::size_t b = 0; b < M.net.states(); b++)
        if (M.net.edges[b / 2].link >= 0) into_links += Mt(a, b);
    std::printf("transitions into the candidate links' states: %.3e\n", into_links);
    return 0;
  }
  if (wants("week")) {
    // The macroscopic reference, twice: the flow study's own mode_flows on every segment, and the module's expected flow.
    std::vector<std::size_t> all_edges;
    for (std::size_t e = 0; e < E; e++) all_edges.push_back(B.theirs[e]);
    const auto terms = M.terms(Tr, all_edges);
    const auto draw = M.frontage_draw(P.as_built());
    std::vector<double> theirs(H * S, 0.0), mine(H * S, 0.0);
    std::vector<std::array<double, flow::modes>> by_mode(H);
    for (std::size_t h = 0; h < H; h++) {
      const auto fl = M.mode_flows(terms, draw, h);
      for (std::size_t m = 0; m < flow::modes; m++)
        for (std::size_t k = 0; k < all_edges.size(); k++) {
          theirs[h * S + 2 * k] += fl[m][k][0], theirs[h * S + 2 * k + 1] += fl[m][k][1];
          by_mode[h][m] += fl[m][k][0] + fl[m][k][1];
        }
    }
    calx::scheduler::parallel_for(pool, H, 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t h = b; h < e; h++) {
        const auto f = um::expected_flow(*V, T, trips_at(week_streams, h));
        std::copy(f.begin(), f.end(), mine.begin() + static_cast<std::ptrdiff_t>(h * S));
      }
    });
    double gap = 0, top = 0, tot_theirs = 0;
    for (std::size_t i = 0; i < H * S; i++) gap = std::max(gap, std::fabs(mine[i] - theirs[i])), top = std::max(top, theirs[i]), tot_theirs += theirs[i];
    std::printf("the module's expected flow against the flow study's mode_flows: largest difference %.3e of the largest flow (%.0f walker-segments an hour); %.0f walker-segments a week  [%.0f s]\n",
                gap / top, top, tot_theirs, seconds_since(T0));

    // Mode by mode at one hour (Monday 13:00), where the two differ.
    {
      const std::size_t h = 13;
      const auto fl = M.mode_flows(terms, draw, h);
      for (std::size_t m = 0; m < flow::modes; m++) {
        std::vector<um::Stream> only;
        for (const auto& st : week_streams)
          if (static_cast<std::size_t>(st.mode) == m) only.push_back(st);
        const auto f = um::expected_flow(*V, T, trips_at(only, h));
        double g = 0, tp = 0, sm = 0, st = 0;
        for (std::size_t k = 0; k < all_edges.size(); k++)
          for (int d = 0; d < 2; d++) g = std::max(g, std::fabs(f[2 * k + d] - fl[m][k][d])), tp = std::max(tp, fl[m][k][d]), sm += f[2 * k + d], st += fl[m][k][d];
        std::printf("  Monday 13:00, %-14s largest difference %.3e of its largest flow %.1f; totals %.3f and %.3f\n", flow::mode_name(m), tp > 0 ? g / tp : 0.0, tp, sm, st);
      }
    }
    // The walkers.
    um::Options o;
    o.seed = 20261105, o.hours = H + 1, o.minute_bins = 60;
    um::Simulation sim(Tp, {V}, week_streams, counters, o);
    Json snapshots = Json::object();
    const auto t1 = Clock::now();
    std::size_t peak_active = 0;
    for (std::size_t h = 1; h <= H + 1; h++) {
      if (h == 57) {   // Wednesday 08:30
        sim.run(56.5 * 3600, &pool);
        snapshots.set("Wednesday 08:30", positions_json(sim));
      }
      if (h == 134) {  // Saturday 13:30
        sim.run(133.5 * 3600, &pool);
        snapshots.set("Saturday 13:30", positions_json(sim));
      }
      sim.run(static_cast<double>(h) * 3600, &pool);
      peak_active = std::max(peak_active, sim.active());
      if (h % 6 == 0) std::printf("  hour %3zu: %zu walking, %.0f trips so far  [%.0f s]\n", h, sim.active(), sim.totals().arrived, seconds_since(t1));
    }
    maps.set("snapshots", snapshots);
    const double run_s = seconds_since(t1);
    const auto& t = sim.tallies();
    // Micro against macro: entries by the hour of departure against the flow, (n − μ)/√μ over every state-hour with μ ≥ 5.
    double chi = 0, zsum = 0, zsq = 0, z3 = 0, used = 0, abm_total = 0;
    std::vector<double> zhist(40, 0.0);   // (n − μ)/√μ in bins of 0.25 from −5 to 5
    for (std::size_t h = 0; h < H; h++)
      for (std::size_t s = 0; s < S; s++) {
        const double mu = theirs[h * S + s], n = t.entries_by_departure[h * S + s];
        abm_total += n;
        if (mu < 5) continue;
        const double z = (n - mu) / std::sqrt(mu);
        chi += z * z, zsum += z, zsq += z * z, z3 += std::fabs(z) > 3, used++;
        if (z >= -5 && z < 5) zhist[static_cast<std::size_t>((z + 5) / 0.25)] += 1;
      }
    // The dynamic lag: entries by the hour they happen against the quasi-static flow.
    std::vector<double> lag_flow, lag_abm;
    for (std::size_t h = 0; h < H; h++) {
      double f = 0, a = 0;
      for (std::size_t s = 0; s < S; s++) f += theirs[h * S + s], a += t.entries[h * S + s];
      lag_flow.push_back(f), lag_abm.push_back(a);
    }
    std::vector<double> by_dep(H, 0.0);
    for (std::size_t h = 0; h < H; h++)
      for (std::size_t s = 0; s < S; s++) by_dep[h] += t.entries_by_departure[h * S + s];
    std::printf("the walkers: %.0f trips in %.0f s on %zu threads (at most %zu walking at once); %.0f entries against the flow's %.0f; over %.0f state-hours (μ ≥ 5) the standardised deviations have mean %.3f and variance %.3f, %.2f%% beyond 3; χ²/dof %.3f\n",
                sim.totals().arrived, run_s, threads, peak_active, abm_total, tot_theirs, used, zsum / used, zsq / used - (zsum / used) * (zsum / used), 100 * z3 / used, chi / used);
    // Per segment over the week (both ways).
    std::ostringstream seg;
    seg << "edge,name,kind,width_m,length_m,flow_week,walkers_week_by_departure,walkers_week_by_entry,mean_walkers,los_e_hours\n";
    for (std::size_t e = 0; e < E; e++) {
      double f = 0, a = 0, b = 0, occ = 0, le = 0;
      for (std::size_t h = 0; h < H; h++) {
        f += theirs[h * S + 2 * e] + theirs[h * S + 2 * e + 1];
        a += t.entries_by_departure[h * S + 2 * e] + t.entries_by_departure[h * S + 2 * e + 1];
        b += t.entries[h * S + 2 * e] + t.entries[h * S + 2 * e + 1];
        occ += t.occupancy[h * E + e], le += t.los_e[h * E + e];
      }
      seg << e << ",\"" << W.segments[e].name << "\"," << W.segments[e].kind << "," << W.segments[e].width << "," << rnd(W.segments[e].length, 10) << "," << rnd(f, 10) << "," << a << "," << b
          << "," << rnd(occ / (H * 3600.0), 1000) << "," << rnd(le / 3600, 100) << "\n";
    }
    write_file(out + "/week-segments.csv", seg.str());
    // The counters.
    Json cj = Json::array();
    std::vector<double> y_obs, y_abm, y_flow;
    for (std::size_t k = 0; k < counters.size(); k++) {
      const auto& c = M.counters[k];
      std::array<double, H> obs{}, abm{}, fl{};
      std::array<std::array<double, H>, 2> obs_w{}, abm_w{}, fl_w{};
      for (std::size_t h = 0; h < H; h++) {
        const auto rd = sim.reading(k, h);
        for (int w = 0; w < 2; w++) {
          abm_w[w][h] = rd[w];
          fl_w[w][h] = counters[k].share * (1 + counters[k].hand[w]) * theirs[h * S + counters[k].way[w]];
          obs_w[w][h] = c.directed ? c.way[w][h] : std::nan("");
        }
        abm[h] = rd[0] + rd[1], fl[h] = fl_w[0][h] + fl_w[1][h], obs[h] = c.gross()[h];
        y_obs.push_back(obs[h]), y_abm.push_back(abm[h]), y_flow.push_back(fl[h]);
      }
      std::vector<double> yo(obs.begin(), obs.end()), ya(abm.begin(), abm.end()), yf(fl.begin(), fl.end());
      cj.push(Json::object({{"id", static_cast<double>(c.id)}, {"name", c.name}, {"street", c.street}, {"directed", c.directed}, {"toward", Json::array({c.toward[0], c.toward[1]})},
                            {"observed", arr(obs, 10)}, {"flow", arr(fl, 10)}, {"walkers", arr(abm, 10)},
                            {"observed_way", Json::array({arr(obs_w[0], 10), arr(obs_w[1], 10)})}, {"flow_way", Json::array({arr(fl_w[0], 10), arr(fl_w[1], 10)})},
                            {"walkers_way", Json::array({arr(abm_w[0], 10), arr(abm_w[1], 10)})}, {"r2_walkers_observed", jr(r2_of(yo, ya), 1e4)},
                            {"r2_flow_observed", jr(r2_of(yo, yf), 1e4)}, {"r2_walkers_flow", jr(r2_of(yf, ya), 1e4)},
                            {"lines", Json::object({{"observed", lines_of(obs)}, {"flow", lines_of(fl)}, {"walkers", lines_of(abm)}})}}));
    }
    const double r2_abm = r2_of(y_obs, y_abm), r2_flow = r2_of(y_obs, y_flow), r2_af = r2_of(y_flow, y_abm);
    std::printf("the counters, 16 x 168 hours: R2 of the walkers against the counts %.4f (the flow's %.4f); of the walkers against the flow %.4f\n", r2_abm, r2_flow, r2_af);
    // The district's walkers minute by minute, and the hours' trips.
    std::vector<double> present(H * 60, 0.0);
    for (std::size_t b = 0; b < H * 60; b++)
      for (std::size_t e = 0; e < E; e++) present[b] += t.minutes[b * E + e] / 60.0;
    std::vector<double> trips(H, 0.0), travel(H, 0.0), delay(H, 0.0);
    for (std::size_t h = 0; h < H; h++)
      for (std::size_t s = 0; s < sim.streams().size(); s++) {
        const std::size_t i = h * sim.streams().size() + s;
        trips[h] += t.trips[i], travel[h] += t.travel[i], delay[h] += t.delay[i];
      }
    for (std::size_t h = 0; h < H; h++) travel[h] = trips[h] > 0 ? travel[h] / trips[h] : 0, delay[h] = trips[h] > 0 ? delay[h] / trips[h] : 0;
    // The walkers by mode, with the flow's.
    std::array<double, flow::modes> mode_trips{}, mode_flow{};
    for (std::size_t s = 0; s < sim.streams().size(); s++)
      for (std::size_t h = 0; h < H; h++) mode_trips[static_cast<std::size_t>(sim.streams()[s].mode)] += t.trips[h * sim.streams().size() + s];
    for (std::size_t h = 0; h < H; h++)
      for (std::size_t m = 0; m < flow::modes; m++) mode_flow[m] += by_mode[h][m];
    Json modes = Json::array();
    for (std::size_t m = 0; m < flow::modes; m++)
      modes.push(Json::object({{"mode", flow::mode_name(m)}, {"trips", mode_trips[m]}, {"flow_walker_segments", jr(mode_flow[m], 1)}}));
    write_file(out + "/week.json",
               Json::object({{"seed", 20261105.0}, {"step_s", o.kinematics.step}, {"threads", static_cast<double>(threads)}, {"run_s", jr(run_s, 10)},
                             {"calibration", Json::object({{"r2", jr(fit.r2, 1e4)}, {"loo_r2", jr(fit.loo_r2, 1e4)}, {"net_r2", jr(fit.net_r2, 1e4)}, {"amplitude_gap", jr(amplitude_gap, 1e8)}})},
                             {"macro_check", Json::object({{"largest_difference_of_largest_flow", jr(gap / top, 1e18)}, {"largest_flow", jr(top, 10)}, {"week_walker_segments", jr(tot_theirs, 1)}})},
                             {"totals", totals_json(sim.totals())}, {"peak_walking", static_cast<double>(peak_active)},
                             {"micro_macro", Json::object({{"state_hours", used}, {"mean_z", jr(zsum / used, 1e4)}, {"var_z", jr(zsq / used - (zsum / used) * (zsum / used), 1e4)},
                                                           {"beyond_3", jr(z3 / used, 1e5)}, {"chi2_per_dof", jr(chi / used, 1e4)}, {"walker_entries", abm_total}, {"z_histogram", arr(zhist, 1)},
                                                           {"flow_entries", jr(tot_theirs, 1)}})},
                             {"entries_by_hour", Json::object({{"flow", arr(lag_flow, 10)}, {"walkers_by_entry", arr(lag_abm, 1)}, {"walkers_by_departure", arr(by_dep, 1)}})},
                             {"r2", Json::object({{"walkers_observed", jr(r2_abm, 1e4)}, {"flow_observed", jr(r2_flow, 1e4)}, {"walkers_flow", jr(r2_af, 1e4)}})},
                             {"counters", cj}, {"present_by_minute", arr(present, 10)}, {"trips_by_hour", arr(trips, 1)}, {"travel_s_by_hour", arr(travel, 10)},
                             {"delay_s_by_hour", arr(delay, 100)}, {"modes", modes}})
                   .dump());
    std::printf("week done  [%.0f s]\n", seconds_since(T0));
  }

  // ===================================================================================================================
  // 2. A street closed
  // ===================================================================================================================

  // The ordinary Sunday, open: shared by the closure and the festival as their baseline (common random numbers).
  um::Options day;
  day.seed = 20170129, day.hours = 25, day.minute_bins = 60;
  um::Tallies sunday_open_tallies;
  um::Totals sunday_open_totals;
  std::vector<std::array<std::array<double, 2>, 24>> sunday_open_readings(counters.size());
  if (wants("closure") || wants("festival")) {
    um::Simulation open(Tp, {RR->values[0][0], RR->values[0][1]}, sunday_streams, counters, day);
    open.run(25 * 3600, &pool);
    sunday_open_tallies = open.tallies(), sunday_open_totals = open.totals();
    for (std::size_t k = 0; k < counters.size(); k++)
      for (std::size_t h = 0; h < 24; h++) sunday_open_readings[k][h] = open.reading(k, h);
    std::printf("the ordinary Sunday, open: %.0f trips, mean %.0f s, %.1f person-hours  [%.0f s]\n", open.totals().arrived, open.totals().travel / open.totals().arrived,
                open.totals().present / 3600, seconds_since(T0));
  }
  if (wants("closure")) {
    um::Simulation closed(Tp, {RR->values[2][0], RR->values[2][1]}, sunday_streams, counters, day);
    closed.run(25 * 3600, &pool);
    const auto& to = sunday_open_tallies;
    const auto& tc = closed.tallies();
    // The exact flows under both, hour by hour.
    std::vector<double> f_open(24 * S, 0.0), f_closed(24 * S, 0.0);
    calx::scheduler::parallel_for(pool, 24, 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t h = b; h < e; h++) {
        // Each class by its own routes: closed, the walkers passing through avoid the street; its own walkers keep it.
        const auto t0 = trips_at(sunday_streams, h, 0), t1 = trips_at(sunday_streams, h, 1);
        auto a = um::expected_flow(*RR->values[0][0], T, t0), c = um::expected_flow(*RR->values[2][0], T, t0);
        const auto a1 = um::expected_flow(*RR->values[0][1], T, t1), c1 = um::expected_flow(*RR->values[2][1], T, t1);
        for (std::size_t q = 0; q < a.size(); q++) a[q] += a1[q], c[q] += c1[q];
        std::copy(a.begin(), a.end(), f_open.begin() + static_cast<std::ptrdiff_t>(h * S));
        std::copy(c.begin(), c.end(), f_closed.begin() + static_cast<std::ptrdiff_t>(h * S));
      }
    });
    std::ostringstream seg;
    seg << "edge,name,flow_open,flow_closed,walkers_open,walkers_closed\n";
    double moved = 0, onto_block = 0;
    for (std::size_t e = 0; e < E; e++) {
      double fa = 0, fc = 0, wa = 0, wc = 0;
      for (std::size_t h = 0; h < 24; h++) {
        fa += f_open[h * S + 2 * e] + f_open[h * S + 2 * e + 1], fc += f_closed[h * S + 2 * e] + f_closed[h * S + 2 * e + 1];
        wa += to.entries[h * S + 2 * e] + to.entries[h * S + 2 * e + 1], wc += tc.entries[h * S + 2 * e] + tc.entries[h * S + 2 * e + 1];
      }
      if (std::find(block.edges.begin(), block.edges.end(), e) != block.edges.end()) onto_block += fa;
      else moved += std::max(0.0, fc - fa);
      seg << e << ",\"" << W.segments[e].name << "\"," << rnd(fa, 10) << "," << rnd(fc, 10) << "," << wa << "," << wc << "\n";
    }
    write_file(out + "/closure-segments.csv", seg.str());
    std::vector<double> ph_open(24, 0.0), ph_closed(24, 0.0), dl_open(24, 0.0), dl_closed(24, 0.0), tr_open(24, 0.0), tr_closed(24, 0.0);
    for (std::size_t h = 0; h < 24; h++) {
      for (std::size_t e = 0; e < E; e++) ph_open[h] += to.occupancy[h * E + e] / 3600, ph_closed[h] += tc.occupancy[h * E + e] / 3600;
      for (std::size_t s = 0; s < sunday_streams.size(); s++) {
        const std::size_t i = h * sunday_streams.size() + s;
        tr_open[h] += to.trips[i], tr_closed[h] += tc.trips[i], dl_open[h] += to.travel[i], dl_closed[h] += tc.travel[i];
      }
    }
    Json cj = Json::array();
    for (std::size_t k = 0; k < counters.size(); k++) {
      std::vector<double> a(24), c(24);
      for (std::size_t h = 0; h < 24; h++) a[h] = sunday_open_readings[k][h][0] + sunday_open_readings[k][h][1], c[h] = closed.reading(k, h)[0] + closed.reading(k, h)[1];
      cj.push(Json::object({{"name", counters[k].name}, {"open", arr(a, 10)}, {"closed", arr(c, 10)}}));
    }
    std::printf("closed to walkers passing through: %.0f trips (%.0f failed), mean %.1f s against %.1f s open; %.1f person-hours against %.1f; the block carried %.0f walker-passes, %.0f moved elsewhere  [%.0f s]\n",
                closed.totals().arrived, closed.totals().failed, closed.totals().travel / closed.totals().arrived, sunday_open_totals.travel / sunday_open_totals.arrived,
                closed.totals().present / 3600, sunday_open_totals.present / 3600, onto_block, moved, seconds_since(T0));
    write_file(out + "/closure.json",
               Json::object({{"street", block.name}, {"edges", Json(std::vector<double>(block.edges.begin(), block.edges.end()))}, {"open", totals_json(sunday_open_totals)},
                             {"closed", totals_json(closed.totals())}, {"flow_on_block", jr(onto_block, 10)}, {"flow_moved_elsewhere", jr(moved, 10)},
                             {"person_hours_by_hour", Json::object({{"open", arr(ph_open, 100)}, {"closed", arr(ph_closed, 100)}})},
                             {"trips_by_hour", Json::object({{"open", arr(tr_open, 1)}, {"closed", arr(tr_closed, 1)}})},
                             {"travel_s_by_hour", Json::object({{"open", arr(dl_open, 1)}, {"closed", arr(dl_closed, 1)}})}, {"counters", cj}})
                   .dump());
  }

  // ===================================================================================================================
  // 3. The festival's crowd
  // ===================================================================================================================
  std::vector<std::size_t> visitors;
  std::vector<double> arrivals(24, 0.0);
  double dwell_h = 1;
  std::vector<um::Stream> festival_streams = sunday_streams;
  const double festival_width = 8.0;   // m: Little Bourke Street's 10.1 m reserve given to walkers, less 2 m of stalls (an assumption)
  std::vector<std::size_t> fit_counters;
  if (wants("festival") || wants("regime")) {
    // The counters' excess on Sunday 29 January 2017 over the Sundays about it.
    const long long festival_day = ur::days_from_civil(2017, 1, 29);
    std::vector<long long> ordinary;
    for (const auto& [m, d] : std::vector<std::pair<unsigned, unsigned>>{{1, 8}, {1, 15}, {1, 22}, {2, 5}, {2, 12}, {2, 19}}) ordinary.push_back(ur::days_from_civil(2017, m, d));
    std::map<long long, std::string> series;
    for (const auto& s : in["sensors"].as_array()) series[static_cast<long long>(s["id"].as_number())] = s["series"].text();
    std::vector<std::array<double, 24>> excess(counters.size()), festival_counts(counters.size()), ordinary_mean(counters.size());
    for (std::size_t k = 0; k < counters.size(); k++) {
      excess[k].fill(std::nan("")), festival_counts[k].fill(std::nan("")), ordinary_mean[k].fill(std::nan(""));
      const auto it = series.find(M.counters[k].id);
      if (it == series.end()) continue;
      const auto r = counts_record(it->second);
      const auto f = ur::observed_day(r, festival_day);
      std::array<double, 24> m{};
      std::array<int, 24> n{};
      for (const auto d : ordinary) {
        const auto o = ur::observed_day(r, d);
        for (std::size_t h = 0; h < 24; h++)
          if (std::isfinite(o[h])) m[h] += o[h], n[h]++;
      }
      bool complete = true;
      for (std::size_t h = 0; h < 24; h++) {
        if (!std::isfinite(f[h]) || n[h] < 4) { complete = false; continue; }
        festival_counts[k][h] = f[h], ordinary_mean[k][h] = m[h] / n[h], excess[k][h] = f[h] - m[h] / n[h];
      }
      if (complete) fit_counters.push_back(k);
    }
    std::printf("the festival's counters with 2017's counts whole: %zu of %zu\n", fit_counters.size(), counters.size());
    // Two patterns of visitors, each with its own arrivals by hour: from the district's entrances (by their flow) and
    // from its buildings (by their floor), each bound for a junction along the street by the length about it, dwelling,
    // and walking back to where it came from. Their unit passes at each counter (out, and back) by the exact flow.
    std::array<um::Stream, 2> vs;
    {
      std::vector<std::pair<std::size_t, double>> from_entrances, from_buildings;
      for (const auto& E : M.entrances) {
        double g = 0;
        for (double x : E.gross) g += x;
        from_entrances.push_back({E.node, g});
      }
      for (const std::size_t b : M.building_nodes) from_buildings.push_back({b, M.work[b] + M.home[b]});
      const std::array<const std::vector<std::pair<std::size_t, double>>*, 2> from{&from_entrances, &from_buildings};
      const char* names[] = {"the festival's visitors from the entrances", "the festival's visitors from the buildings"};
      for (int q = 0; q < 2; q++) {
        um::Stream& v = vs[static_cast<std::size_t>(q)];
        v.name = names[q], v.mode = 7 + q, v.route_class = 1, v.rate.assign(24, 0.0), v.dwell = 3600;
        for (const auto& [o, wo] : *from[static_cast<std::size_t>(q)])
          for (const auto& [target, l] : along_street) v.origin.push_back(o), v.target.push_back(target), v.back.push_back(T.junction(W, o)), v.weight.push_back(wo * l);
      }
    }
    const std::array<um::UnitCounts, 2> unit{um::unit_counts(*RR->values[0][1], T, vs[0], counters), um::unit_counts(*RR->values[0][1], T, vs[1], counters)};
    const auto gross = [&](const std::array<double, 2>& u, std::size_t k) { return counters[k].share * ((1 + counters[k].hand[0]) * u[0] + (1 + counters[k].hand[1]) * u[1]); };
    for (const std::size_t k : fit_counters)
      std::printf("  %-34s a visitor from the entrances passes it %.3f times out, %.3f back; from the buildings %.3f, %.3f\n", counters[k].name.c_str(), gross(unit[0].out[k], k),
                  gross(unit[0].back[k], k), gross(unit[1].out[k], k), gross(unit[1].back[k], k));
    // The arrivals a_qk ≥ 0 (hours 8 to 22, each pattern) and the mean dwell D, by NNLS for each D of a grid, on the
    // hours 8 to 23.
    double best = 1e300;
    Json grid = Json::array();
    std::array<std::vector<double>, 2> pattern_arrivals;
    for (const double D : {0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0}) {
      std::vector<std::vector<double>> rows;
      std::vector<double> b;
      for (const std::size_t k : fit_counters)
        for (std::size_t h = 8; h < 24; h++) {
          std::vector<double> row;
          for (int q = 0; q < 2; q++)
            for (std::size_t a = 8; a < 23; a++) {
              double x = a == h ? gross(unit[static_cast<std::size_t>(q)].out[k], k) : 0.0;
              if (h >= a) x += back_kernel(D, static_cast<int>(h - a)) * gross(unit[static_cast<std::size_t>(q)].back[k], k);
              row.push_back(x);
            }
          rows.push_back(row), b.push_back(excess[k][h]);
        }
      const auto sol = tyche::regression::nnls(calx::linalg::Rectangular::from_rows(rows), b);
      double ss = 0, st = 0, mean = 0;
      for (double y : b) mean += y / static_cast<double>(b.size());
      for (std::size_t i = 0; i < rows.size(); i++) {
        double p = 0;
        for (std::size_t a = 0; a < sol.x.size(); a++) p += rows[i][a] * sol.x[a];
        ss += (b[i] - p) * (b[i] - p), st += (b[i] - mean) * (b[i] - mean);
      }
      double total = 0;
      for (double x : sol.x) total += x;
      grid.push(Json::object({{"dwell_h", D}, {"r2", jr(1 - ss / st, 1e4)}, {"visitors", jr(total, 1)}}));
      std::printf("  dwell %.2f h: R2 of the excess %.4f, %.0f visitors\n", D, 1 - ss / st, total);
      if (ss < best) {
        best = ss, dwell_h = D;
        arrivals.assign(24, 0.0);
        for (int q = 0; q < 2; q++) {
          pattern_arrivals[static_cast<std::size_t>(q)].assign(24, 0.0);
          for (std::size_t a = 8; a < 23; a++) {
            pattern_arrivals[static_cast<std::size_t>(q)][a] = sol.x[static_cast<std::size_t>(q) * 15 + a - 8];
            arrivals[a] += sol.x[static_cast<std::size_t>(q) * 15 + a - 8];
          }
        }
      }
    }
    for (int q = 0; q < 2; q++) {
      vs[static_cast<std::size_t>(q)].dwell = dwell_h * 3600, vs[static_cast<std::size_t>(q)].rate = pattern_arrivals[static_cast<std::size_t>(q)];
      festival_streams.push_back(vs[static_cast<std::size_t>(q)]);
      visitors.push_back(festival_streams.size() - 1);
    }
    double total = 0;
    for (double a : arrivals) total += a;
    std::printf("the festival: dwell %.2f h, %.0f visitors, arrivals by hour fitted  [%.0f s]\n", dwell_h, total, seconds_since(T0));

    if (wants("festival")) {
      um::Simulation fest(Tp, {RR->values[0][0], RR->values[0][1]}, festival_streams, counters, day);
      for (const std::size_t e : block.edges) fest.set_width(e, festival_width);
      Json peak_positions;
      double peak = -1;
      std::vector<double> block_walkers(24 * 60, 0.0);
      for (std::size_t minute = 1; minute <= 25 * 60; minute++) {
        fest.run(60.0 * static_cast<double>(minute), &pool);
        double on = 0;
        for (const auto& w : fest.walkers())
          if (w.state != static_cast<std::uint32_t>(none) && std::find(block.edges.begin(), block.edges.end(), w.state / 2) != block.edges.end()) on += 1;
        if (minute <= 24 * 60) block_walkers[minute - 1] = on;
        if (on > peak && minute <= 24 * 60) {
          peak = on;
          Json pos = Json::array();
          for (const auto& p : fest.positions()) pos.push(Json::array({jr(p.at[0], 10), jr(p.at[1], 10), p.mode, static_cast<int>(p.status)}));
          peak_positions = Json::object({{"minute", static_cast<double>(minute)}, {"positions", pos}});
        }
      }
      maps.set("festival_peak", peak_positions);
      const auto& tf = fest.tallies();
      // Space–time along the block, minute by minute: walkers and those at Fruin's level E on each of its segments; the
      // queues at the mouths of every segment touching it.
      Json along = Json::array();
      for (const std::size_t e : block.edges) {
        std::vector<double> m(24 * 60), me(24 * 60);
        for (std::size_t b = 0; b < 24 * 60; b++) m[b] = tf.minutes[b * E + e] / 60.0, me[b] = tf.minutes_e[b * E + e] / 60.0;
        along.push(Json::object({{"edge", static_cast<double>(e)}, {"length", jr(W.segments[e].length, 10)}, {"walkers", arr(m, 10)}, {"at_level_e", arr(me, 10)}}));
      }
      std::set<std::size_t> near;
      for (const std::size_t e : block.edges)
        for (const std::size_t x : {W.segments[e].u, W.segments[e].v})
          for (const std::size_t s : W.entering[x]) near.insert(s);
      Json queues = Json::array();
      for (std::size_t s = 0; s < S; s++) {
        std::vector<double> q(24), qm(24);
        double any = 0;
        for (std::size_t h = 0; h < 24; h++) q[h] = tf.queue[h * S + s] / 60, qm[h] = tf.queue_max[h * S + s], any += qm[h];
        if (any >= 50 || (near.count(s) && any > 0))
          queues.push(Json::object({{"state", static_cast<double>(s)}, {"name", W.segments[s / 2].name}, {"walker_minutes", arr(q, 10)}, {"longest", arr(qm, 1)}}));
      }
      Json cj = Json::array();
      for (std::size_t k = 0; k < counters.size(); k++) {
        std::vector<double> sim_excess(24), obs_excess(24), ord(24), fes(24);
        for (std::size_t h = 0; h < 24; h++) {
          const auto r = fest.reading(k, h);
          sim_excess[h] = r[0] + r[1] - sunday_open_readings[k][h][0] - sunday_open_readings[k][h][1];
          obs_excess[h] = excess[k][h], ord[h] = ordinary_mean[k][h], fes[h] = festival_counts[k][h];
        }
        cj.push(Json::object({{"name", counters[k].name}, {"id", static_cast<double>(M.counters[k].id)}, {"fitted", std::find(fit_counters.begin(), fit_counters.end(), k) != fit_counters.end()},
                              {"observed_2017", arr(fes, 10)}, {"ordinary_sundays_2017", arr(ord, 10)}, {"observed_excess", arr(obs_excess, 10)}, {"walkers_excess", arr(sim_excess, 10)}}));
      }
      // The record of the festival day, as simulated events.
      const auto rec = um::simulated_record(fest, {{block, um::Regime::open, 0, 0}}, {.mjd0 = urban::events::mjd(2017, 1, 28) + 13.0 / 24});   // Sunday 00:00 AEDT
      const auto problems = urban::events::problems(rec);
      Json ev = Json::array();
      std::map<std::string, int> kinds;
      for (const auto& e : rec.events) {
        kinds[e.kind]++;
        if (e.kind == "count") continue;
        Json cs = Json::array();
        for (auto c : e.couplings) cs.push(rec.couplings[c].id + "|" + rec.couplings[c].name);
        ev.push(Json::object({{"id", e.id}, {"type", e.type}, {"kind", e.kind}, {"from_mjd", jr(e.when.from, 1e6)}, {"to_mjd", jr(e.when.to, 1e6)}, {"couplings", cs},
                              {"source", e.source}, {"provenance", e.provenance}}));
      }
      Json kj = Json::object();
      for (const auto& [k, n] : kinds) kj.set(k, n);
      std::printf("the festival day: %.0f trips (%.0f failed), %.0f person-hours (%.1f at level E, %.1f at F), at most %.0f walkers on the block (minute %.0f); the record: %zu events, %s  [%.0f s]\n",
                  fest.totals().arrived, fest.totals().failed, fest.totals().present / 3600, fest.totals().los_e / 3600, fest.totals().los_f / 3600, peak, peak_positions["minute"].as_number(),
                  rec.events.size(), problems.empty() ? "sound" : problems[0].c_str(), seconds_since(T0));
      write_file(out + "/festival.json",
                 Json::object({{"day", "Sunday 29 January 2017"}, {"permit", "City of Melbourne event permits, 2017-01-27 to 2017-01-29: Little Bourke Street, Little Collins Street (Public Event, non-ticketed)"},
                               {"block", block.name}, {"block_width_m", festival_width}, {"dwell_grid", grid}, {"dwell_h", dwell_h}, {"arrivals_by_hour", arr(arrivals, 10)},
                               {"arrivals_from_entrances", arr(festival_streams[visitors[0]].rate, 10)}, {"arrivals_from_buildings", arr(festival_streams[visitors[1]].rate, 10)},
                               {"fit_counters", static_cast<double>(fit_counters.size())}, {"totals", totals_json(fest.totals())}, {"ordinary", totals_json(sunday_open_totals)},
                               {"block_walkers_by_minute", arr(block_walkers, 1)}, {"along", along}, {"queues", queues}, {"counters", cj},
                               {"record", Json::object({{"events", static_cast<double>(rec.events.size())}, {"kinds", kj}, {"sound", problems.empty()}, {"list", ev}})}})
                     .dump());
    }
  }

  // ===================================================================================================================
  // 4. The regime decided online
  // ===================================================================================================================
  if (wants("regime")) {
    // Sensed: the counters on Little Bourke Street at Swanston Street (19) and on Russell Street (21).
    std::vector<std::size_t> sensed;
    for (std::size_t k = 0; k < counters.size(); k++)
      if (M.counters[k].id == 19 || M.counters[k].id == 21) sensed.push_back(k);
    const std::vector<double> scales{0.6, 1.0, 1.6}, prior{0.25, 0.5, 0.25};
    const double epoch = 1800, start_h = 9;
    const std::size_t epochs = 16;   // 09:00 to 17:00
    // The bins of an epoch's reading at a sensed counter: from an open festival day at the expected scale, the readings'
    // terciles over the decision hours.
    std::vector<double> bins;
    {
      um::Options o = day;
      o.light = true;
      um::Simulation probe(Tp, {RR->values[0][0], RR->values[0][1]}, festival_streams, counters, o);
      for (const std::size_t e : block.edges) probe.set_width(e, festival_width);
      probe.run(start_h * 3600, &pool);
      std::vector<double> reads;
      for (std::size_t q = 0; q < epochs; q++) {
        const auto before = probe.passes();
        probe.run(probe.time() + epoch, &pool);
        for (const std::size_t k : sensed) {
          const auto& c = counters[k];
          reads.push_back(c.share * ((1 + c.hand[0]) * (probe.passes()[2 * k] - before[2 * k]) + (1 + c.hand[1]) * (probe.passes()[2 * k + 1] - before[2 * k + 1])));
        }
      }
      std::sort(reads.begin(), reads.end());
      bins = {reads[reads.size() / 3], reads[2 * reads.size() / 3]};
    }
    std::printf("regimes: sensed %zu counters, bins of a half hour's reading %.0f and %.0f\n", sensed.size(), bins[0], bins[1]);
    struct Case {
      std::size_t scale = 0, rep = 0;
      int policy = 0;   // 0 open, 1 one-way, 2 closed, 3 POMCP
      std::vector<int> actions;
      std::vector<double> rewards, person_hours, los_f;
      double value = 0, seconds = 0;
      int restarts = 0;
      um::Totals totals;
    };
    std::vector<Case> cases;
    for (std::size_t sc = 0; sc < scales.size(); sc++)
      for (std::size_t rep = 0; rep < replicates; rep++)
        for (int pol = 0; pol < 4; pol++) cases.push_back({sc, rep, pol, {}, {}, {}, {}, 0, 0, 0, {}});
    um::RegimeReward reward{1, 10};
    const auto world_of = [&](std::size_t sc, std::size_t rep) {
      um::Options o = day;
      o.seed = 20170129 + 1000 * (rep + 1), o.light = true;
      um::Simulation w(Tp, {RR->values[0][0], RR->values[0][1]}, festival_streams, counters, o);
      for (const std::size_t e : block.edges) w.set_width(e, festival_width);
      for (const std::size_t v : visitors) w.set_scale(v, scales[sc]);
      return w;
    };
    const auto policy_name = [](int p) { return std::string(p == 3 ? "POMCP" : um::regime_name(static_cast<um::Regime>(p))); };
    const auto case_json = [&](const Case& c) {
      return Json::object({{"scale", scales[c.scale]}, {"day", static_cast<double>(c.rep)}, {"policy", policy_name(c.policy)},
                           {"actions", Json(std::vector<double>(c.actions.begin(), c.actions.end()))}, {"rewards", arr(c.rewards, 100)}, {"person_hours", arr(c.person_hours, 100)},
                           {"los_f_hours", arr(c.los_f, 1000)}, {"value", jr(c.value, 100)}, {"restarts", c.restarts}, {"seconds", jr(c.seconds, 10)}, {"totals", totals_json(c.totals)}});
    };
    // Resumable: each case is appended to regime-cases.jsonl as it ends, and a case found there is not run again. A case
    // depends only on its scale, its day and its policy (and its index, POMCP's seed), so a resumed run gives the cases
    // the whole run would.
    const std::string journal = out + "/regime-cases.jsonl";
    std::vector<Json> done(cases.size());
    {
      std::ifstream in(journal);
      std::string line;
      while (std::getline(in, line)) {
        if (line.empty()) continue;
        const Json j = Json::parse(line);
        for (std::size_t ci = 0; ci < cases.size(); ci++)
          if (std::fabs(scales[cases[ci].scale] - j["scale"].as_number()) < 1e-9 && static_cast<double>(cases[ci].rep) == j["day"].as_number() &&
              policy_name(cases[ci].policy) == j["policy"].text())
            done[ci] = j;
      }
    }
    // --regime-cases N runs at most N of the cases not yet journalled, in order (scale, then day, then policy).
    std::vector<std::size_t> pending;
    for (std::size_t ci = 0; ci < cases.size(); ci++)
      if (done[ci].is_null() && pending.size() < regime_cases) pending.push_back(ci);
    std::printf("regimes: %zu of %zu cases already in %s\n", cases.size() - pending.size(), cases.size(), "regime-cases.jsonl");
    std::mutex journal_lock;
    calx::scheduler::parallel_for(pool, pending.size(), 1, [&](std::size_t p0, std::size_t p1, std::size_t) {
      for (std::size_t pi = p0; pi < p1; pi++) {
        const std::size_t ci = pending[pi];
        Case& c = cases[ci];
        const auto t0 = Clock::now();
        um::Simulation world = world_of(c.scale, c.rep);
        world.run(start_h * 3600);
        um::RegimeModel::Setup setup;
        setup.routes = RR, setup.epoch = epoch, setup.epochs = epochs, setup.sensed = sensed, setup.bins = bins, setup.scaled = visitors;
        setup.scales = scales, setup.prior = prior, setup.reward = reward, setup.discount = 0.95;
        um::RegimeModel model(world.lightened(4.0), setup);
        std::unique_ptr<tyche::decision::Pomcp<um::RegimeModel>> planner;
        const auto fresh = [&](const um::Simulation& from) {
          model = um::RegimeModel(from.lightened(4.0), setup);
          planner = std::make_unique<tyche::decision::Pomcp<um::RegimeModel>>(
              model, tyche::decision::PomcpOptions{.simulations = 12, .reward_span = 300, .depth = 2, .particles = 6, .refill_attempts = 12, .seed = 77 + ci});
          planner->set_leaf([](const um::RegimeModel::State&) { return 0.0; });
        };
        if (c.policy == 3) fresh(world);
        const um::Totals start = world.totals();
        for (std::size_t q = 0; q < epochs; q++) {
          const std::size_t a = c.policy < 3 ? static_cast<std::size_t>(c.policy) : planner->plan();
          const um::Totals before = world.totals();
          const auto [obs, r] = model.act(world, a);
          c.actions.push_back(static_cast<int>(a)), c.rewards.push_back(r), c.value += r;
          c.person_hours.push_back((world.totals().present - before.present) / 3600), c.los_f.push_back((world.totals().los_f - before.los_f) / 3600);
          if (c.policy == 3 && q + 1 < epochs && !planner->advance(a, obs)) c.restarts++, fresh(world);
        }
        c.totals = world.totals();
        c.totals.present -= start.present, c.totals.los_f -= start.los_f, c.totals.los_e -= start.los_e, c.totals.waiting -= start.waiting;
        c.seconds = seconds_since(t0);
        std::printf("  scale %.1f, day %zu, %-8s value %.1f (%.1f person-hours, %.2f at level F) in %.0f s; regimes", scales[c.scale], c.rep,
                    c.policy == 3 ? "POMCP" : um::regime_name(static_cast<um::Regime>(c.policy)), c.value, c.totals.present / 3600, c.totals.los_f / 3600, c.seconds);
        for (int a : c.actions) std::printf(" %d", a);
        std::printf("%s\n", c.restarts ? (" (" + std::to_string(c.restarts) + " restarts)").c_str() : "");
        std::fflush(stdout);
        const std::lock_guard<std::mutex> hold(journal_lock);
        std::ofstream j(journal, std::ios::app);
        j << case_json(c).dump() << "\n";
      }
    });
    // regime.json holds the cases journalled or run now; a case not yet run is left out, and named.
    Json cj = Json::array(), missing = Json::array();
    for (std::size_t ci = 0; ci < cases.size(); ci++) {
      const bool ran = std::find(pending.begin(), pending.end(), ci) != pending.end();
      if (!done[ci].is_null()) cj.push(done[ci]);
      else if (ran) cj.push(case_json(cases[ci]));
      else missing.push(Json::object({{"scale", scales[cases[ci].scale]}, {"day", static_cast<double>(cases[ci].rep)}, {"policy", policy_name(cases[ci].policy)}}));
    }
    write_file(out + "/regime.json",
               Json::object({{"street", block.name}, {"epoch_s", epoch}, {"epochs", static_cast<double>(epochs)}, {"start_h", start_h}, {"scales", Json(scales)}, {"prior", Json(prior)},
                             {"bins", Json(bins)}, {"reward", Json::object({{"per_person_hour", reward.time}, {"per_person_hour_at_level_f", reward.crush}})},
                             {"pomcp", Json::object({{"simulations", 12}, {"depth", 2}, {"particles", 6}, {"refill_attempts", 12}, {"model_step_s", 4.0}, {"exploration", 300}})},
                             {"cases", cj}, {"not_run", missing}})
                   .dump());
    std::printf("regime done  [%.0f s]\n", seconds_since(T0));
  }
  // Kept with what earlier runs of other tasks drew (the week's snapshots, the festival's peak).
  {
    std::ifstream old(out + "/maps.json");
    if (old) {
      std::stringstream b;
      b << old.rdbuf();
      Json merged = Json::parse(b.str());
      for (const auto& [key, value] : maps.as_object()) merged.set(key, value);
      maps = std::move(merged);
    }
  }
  write_file(out + "/maps.json", maps.dump());
  std::printf("all done  [%.0f s]\n", seconds_since(T0));
  return 0;
}
