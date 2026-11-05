// The district's walkers prepared: the walking network the page sends (Landscape Dynamics' district, cut from the city
// model by src/district.js), the City's counters on it, the entrances where its streets leave its box, the floor of its
// buildings and CLUE's establishments, the streams of trips, and their departures calibrated to the counters each way.
//
// The patterns and their calibration are the district-economy flow study's (research/experiments/2026-11-05-a-
// district-economy-as-a-flow, README §§2–3), on the urban-science domain's movement module, whose route law is that
// study's (θ = 0.02 per metre, widths by Hoddle's classes of street, w_ref = 12 m):
//   - entrances: the network's dead ends within 10 m of its box, each with the gross flow, both ways, of the nearest
//     counter outside the box on its own street within 400 m; else, for a street or a mall, any counter within 120 m,
//     scaled by width; else the median of the street counters outside, scaled by width;
//   - buildings: each building standing in the box carries its floor (footprint × storeys of 3.2 m) on the two ends of
//     its nearest segment, by its place along it;
//   - attractors: CLUE's establishments of the trading kinds (retail, food and drink, services, culture) in the box,
//     open by their kind's typical hours (the problem's: OpenStreetMap's declared hours, median), drawing visits with
//     weights 1, 1, ½, 1;
//   - streams: through trips from each entrance to the others; arrivals from each entrance to the buildings by floor;
//     departures from the buildings to each entrance; local trips between buildings ∝ m_a m_b exp(−d/200 m); visits from
//     the entrances and the buildings (half and half) to the open attractors, and back;
//   - calibration: the amplitudes of through, arrive, leave and local fitted hour by hour, each way, by tyche's NNLS
//     (urban::movement::calibrate); the visits set to ν = 0.1 of the hour's trips (a fixed point, four rounds); each
//     footpath's handedness ρ fitted by least squares over the week, alternating with the amplitudes (three rounds).
// A stream's law over its pairs is fixed (movement's Stream), so where the study's law varies by the hour the week's mean
// stands for it: the through trips from entrance i go to the others in proportion to their mean flows, at i's own rate
// at each hour (the study weights them by the hour's flows). This is the one departure from the study.
//
// Signals: the junctions where two of Hoddle's main streets cross (an assumption: the City's signal timings are not in
// the data), each crossed in two phases of a 120 s cycle (movement's walkers.hpp).

#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <science/json.hpp>
#include <urban/landscape/district.hpp>
#include <urban/movement/movement.hpp>

#include "data.hpp"

namespace hf {

namespace um = urban::movement;
namespace ul = urban::landscape;
inline constexpr std::size_t none = um::none;
inline constexpr std::size_t H = 168;

// The patterns of trips (the flow study's modes, homes' left out as the study left them) and an event's crowd.
enum Mode : int { through = 0, arrive = 1, leave = 2, local = 3, visit = 4, crowd = 5, modes = 6 };
inline const char* mode_name(int m) {
  static const char* n[] = {"passing through", "arriving at work", "leaving work", "between buildings", "visiting shops", "an event's crowd"};
  return n[m];
}

// The trading kinds and their weights (the study's w_k) and typical hours (the problem's typical_hours: weekday and
// weekend opening and closing, civil hours).
inline constexpr int kinds = 4;
inline const char* kind_name(int k) {
  static const char* n[] = {"retail", "food and drink", "services", "culture"};
  return n[k];
}
inline constexpr std::array<double, kinds> attraction{1, 1, 0.5, 1};
inline constexpr std::array<std::array<double, 4>, kinds> typical{{{9, 18, 10, 18}, {11, 22, 11, 22}, {9, 18, 9, 15}, {10, 22, 10, 22}}};
inline bool open_at(int k, std::size_t h) {
  const bool weekend = h / 24 >= 5;
  const double hd = static_cast<double>(h % 24);
  const auto& t = typical[static_cast<std::size_t>(k)];
  return weekend ? hd >= t[2] && hd < t[3] : hd >= t[0] && hd < t[1];
}
// CLUE's uses (the data library's terrain classes) to the trading kinds; −1 for the others.
inline int kind_of_use(const std::string& u) {
  if (u == "clothing-and-fashion" || u == "jewellery" || u == "bridal-and-wedding" || u == "other-retail") return 0;
  if (u == "food-and-drink") return 1;
  if (u == "personal-services" || u == "finance" || u == "health") return 2;
  if (u == "entertainment-and-culture") return 3;
  return -1;
}

inline double side_share(const std::string& kind) { return kind == "laneway" ? 1.0 : 0.5; }

struct Options {
  double entrance_edge = 10, same_street = 400, near_counter = 120, local_length = 200;   // m
  double visit_share = 0.1;
  double storey = 3.2;
  double hand_limit = 0.6;
};

struct CounterInfo {
  long long id = 0;
  std::string name, street;
  um::Point at{};
  std::size_t edge = 0;
  double along = 0, cosine = 0, distance = 0;
  std::array<std::string, 2> toward;
  std::array<std::array<double, H>, 2> way{};   // the mean week each way (the current meter)
  std::array<double, H> gross{};                // the tide's mean week
};

struct Entrance {
  std::size_t node = 0;
  std::string street, kind, source;
  long long counter = -1;
  double distance = 0, width = 0, mean = 0;
  std::array<double, H> gross{};
};

struct NamedStreet {
  std::string name;
  std::vector<std::size_t> edges;
  std::vector<std::size_t> ends;   // junctions of the street with one of its segments: its ends (entries of one way)
  double length = 0;
  std::size_t target = none;
};

struct Calibration {
  double r2_way = 0, r2_gross = 0, visit_share = 0;
  std::size_t rows = 0;
  std::vector<std::vector<double>> amplitude;                 // [group][hour]
  std::vector<std::array<std::array<double, H>, 2>> fitted;   // [counter][way][hour]
  std::vector<double> counter_r2;
  double seconds = 0;
};

class District {
 public:
  std::string key;
  Options o;
  std::array<double, 4> box{};   // x0, y0, x1, y1
  std::shared_ptr<um::Walkway> walkway;
  std::shared_ptr<const um::Walkway> W;
  std::vector<std::string> dropped;
  std::vector<char> main;
  um::RouteLaw law;
  std::shared_ptr<um::Targets> targets;
  std::shared_ptr<const um::Targets> T;
  std::vector<std::size_t> junction_target;
  std::vector<NamedStreet> streets;
  std::shared_ptr<const um::Routes> R;
  std::shared_ptr<const um::Values> V;
  double radius = 0;
  std::vector<CounterInfo> info;
  std::vector<um::Counter> counters;
  std::vector<std::string> left_out;     // counters inside the box not used, with the reason
  std::vector<Entrance> entrances;
  std::vector<double> floor;             // m² per junction
  std::array<std::vector<double>, kinds> attract;   // establishments per junction, by kind
  std::size_t establishments = 0;
  std::vector<double> origin;            // the visits' origins, a law over junctions
  struct Placed { long long id; um::Point at; std::string way; };
  std::vector<Placed> placed;            // the sensors the page placed, with their nearest way
  std::vector<um::Stream> week;          // the streams over the 168 hours of the week, calibrated
  std::vector<std::vector<std::size_t>> groups;
  Calibration fit;
  double seconds = 0;

  calx::scheduler::ThreadPool& pool() const {
    static calx::scheduler::ThreadPool p(std::clamp<std::size_t>(std::thread::hardware_concurrency(), 1, 6));
    return p;
  }

  bool inside(const um::Point& p) const { return p[0] >= box[0] && p[0] <= box[2] && p[1] >= box[1] && p[1] <= box[3]; }

  // The nearest segment on the largest component (of a named street, when given).
  um::Nearest nearest(const um::Point& p, const std::string& street = "") const {
    um::Nearest best;
    for (std::size_t e = 0; e < W->edges(); e++) {
      if (!main[W->segments[e].u]) continue;
      if (!street.empty() && W->segments[e].name != street) continue;
      const auto& l = W->segments[e].line;
      double run = 0;
      for (std::size_t k = 0; k + 1 < l.size(); k++) {
        const double ex = l[k + 1][0] - l[k][0], ey = l[k + 1][1] - l[k][1], L2 = ex * ex + ey * ey, L = std::sqrt(L2);
        const double t = L2 > 0 ? std::clamp(((p[0] - l[k][0]) * ex + (p[1] - l[k][1]) * ey) / L2, 0.0, 1.0) : 0.0;
        const double d = std::hypot(p[0] - l[k][0] - t * ex, p[1] - l[k][1] - t * ey);
        if (d < best.distance) best = {e, run + t * L, d};
        run += L;
      }
    }
    return best;
  }
  // A mass at a point carried onto the two ends of its nearest segment, by its place along it.
  void attach(const um::Point& p, double mass, std::vector<double>& to) const {
    const auto n = nearest(p);
    if (n.edge == none) return;
    const auto& e = W->segments[n.edge];
    const double t = e.length > 0 ? std::clamp(n.along / e.length, 0.0, 1.0) : 0.5;
    to[e.u] += mass * (1 - t), to[e.v] += mass * t;
  }

  explicit District(const Json& r, std::string k) : key(std::move(k)) {
    const auto t0 = std::chrono::steady_clock::now();
    read_network(r);
    read_counters(r);
    read_entrances(r);
    read_masses(r);
    make_streams();
    calibrate();
    seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }

 private:
  void read_network(const Json& r) {
    const Json& b = r["box"];
    box = {b["x0"].as_number(), b["y0"].as_number(), b["x1"].as_number(), b["y1"].as_number()};
    ul::Network N;
    for (const auto& p : r["junctions"].as_array()) N.junctions.push_back({p[0].as_number(), p[1].as_number()});
    for (const auto& s : r["segments"].as_array()) {
      ul::Segment g;
      g.u = static_cast<std::size_t>(s["u"].as_number()), g.v = static_cast<std::size_t>(s["v"].as_number());
      if (s["line"].is_array())
        for (const auto& p : s["line"].as_array()) g.line.push_back({p[0].as_number(), p[1].as_number()});
      g.kind = s["kind"].text("street"), g.name = s["name"].text();
      if (g.kind != "street" && g.kind != "mall" && g.kind != "laneway") g.kind = "street";
      N.segments.push_back(std::move(g));
    }
    if (N.segments.empty() || N.segments.size() > 3000) throw std::invalid_argument("a district of 1 to 3000 segments");
    N = ul::planarised(std::move(N), &dropped);
    walkway = std::make_shared<um::Walkway>(um::walkway(N));
    // Signals where two of Hoddle's main streets cross (an assumption).
    for (std::size_t x = 0; x < walkway->junctions(); x++) {
      std::set<std::string> names;
      for (const std::size_t s : walkway->leaving[x]) {
        const auto& g = walkway->segments[s / 2];
        if (g.kind == "street" && (um::main_street(g.name) || g.name == "Swanston Street")) names.insert(g.name);
      }
      walkway->signalled[x] = walkway->carriageways[x] >= 3 && names.size() >= 2;
    }
    W = walkway;
    main = um::largest_component(*W);
    R = std::make_shared<const um::Routes>(W, law);
    radius = R->spectral_radius(1200);
    if (!(radius < 1)) throw std::runtime_error("the walk sums diverge on this network (spectral radius " + std::to_string(radius) + ")");
    // Targets: every junction of the largest component, and each named street's segments there.
    targets = std::make_shared<um::Targets>();
    junction_target.assign(W->junctions(), none);
    for (std::size_t j = 0; j < W->junctions(); j++)
      if (main[j] && !W->entering[j].empty()) junction_target[j] = targets->junction(*W, j);
    std::map<std::string, NamedStreet> by;
    for (std::size_t e = 0; e < W->edges(); e++) {
      const auto& g = W->segments[e];
      if (g.name.empty() || !main[g.u]) continue;
      auto& s = by[g.name];
      s.name = g.name, s.edges.push_back(e), s.length += g.length;
    }
    for (auto& [name, s] : by) {
      std::map<std::size_t, int> degree;
      for (const std::size_t e : s.edges) degree[W->segments[e].u]++, degree[W->segments[e].v]++;
      for (const auto& [j, d] : degree)
        if (d == 1) s.ends.push_back(j);
      s.target = targets->segments(s.edges);
      streets.push_back(std::move(s));
    }
    std::sort(streets.begin(), streets.end(), [](const NamedStreet& a, const NamedStreet& b) { return a.length > b.length; });
    T = targets;
    V = std::make_shared<const um::Values>(R, *T, &pool());
  }

  struct Raw {
    long long id;
    std::string name, street, kind;
    um::Point at;
    const Weeks* weeks;
  };
  std::vector<Raw> raw_;

  void read_counters(const Json& r) {
    const auto& S = sensors();
    for (const auto& s : r["sensors"].as_array()) {
      const long long id = static_cast<long long>(s["id"].as_number());
      placed.push_back({id, {s["at"][0].as_number(), s["at"][1].as_number()}, s["way"]["name"].text()});
      const auto it = S.find(id);
      if (it == S.end() || it->second.location_type != "Outdoor") continue;
      Raw x{id, it->second.name, s["way"]["name"].text(), s["way"]["kind"].text("street"), {s["at"][0].as_number(), s["at"][1].as_number()}, nullptr};
      const bool crossing = x.name.find("Crossing") != std::string::npos;
      const bool near = x.at[0] > box[0] - 450 && x.at[0] < box[2] + 450 && x.at[1] > box[1] - 450 && x.at[1] < box[3] + 450;
      if (!near) continue;
      if (crossing) {
        if (inside(x.at)) left_out.push_back(x.name + ": it counts people crossing a street");
        continue;
      }
      x.weeks = &weeks(id);
      if (!x.weeks->usable) {
        if (inside(x.at)) left_out.push_back(x.name + ": it covers " + std::to_string(static_cast<int>(std::round(100 * x.weeks->covered))) + "% of the reference year");
        continue;
      }
      raw_.push_back(x);
    }
    const auto compass = [](const std::string& c) -> um::Point {
      if (c == "North") return {0, 1};
      if (c == "South") return {0, -1};
      if (c == "East") return {1, 0};
      return {-1, 0};
    };
    for (const auto& x : raw_) {
      if (!inside(x.at)) continue;
      auto n = nearest(x.at, x.street);
      if (n.edge == none || n.distance > 40) n = nearest(x.at);
      if (n.edge == none) continue;
      if (!x.weeks->directed) {
        left_out.push_back(x.name + ": it does not count each way");
        continue;
      }
      const auto& seg = W->segments[n.edge];
      const um::Point a = W->point(2 * n.edge, std::max(0.0, n.along - 3)), b = W->point(2 * n.edge, std::min(seg.length, n.along + 3));
      const double tx = b[0] - a[0], ty = b[1] - a[1], tn = std::hypot(tx, ty);
      const um::Point cv = compass(x.weeks->toward[0]);
      const double cosine = tn > 0 ? (tx * cv[0] + ty * cv[1]) / tn : 0;
      if (std::fabs(cosine) < 0.5) {
        left_out.push_back(x.name + ": its directions run across its segment");
        continue;
      }
      CounterInfo c;
      c.id = x.id, c.name = x.name, c.street = seg.name, c.at = x.at, c.edge = n.edge, c.along = n.along, c.cosine = cosine, c.distance = n.distance;
      c.toward = x.weeks->toward, c.way = x.weeks->way, c.gross = x.weeks->gross;
      um::Counter k;
      k.name = x.name, k.edge = n.edge, k.along = n.along, k.share = side_share(seg.kind), k.hand = {0, 0};
      k.way = cosine > 0 ? std::array<std::size_t, 2>{2 * n.edge, 2 * n.edge + 1} : std::array<std::size_t, 2>{2 * n.edge + 1, 2 * n.edge};
      info.push_back(std::move(c));
      counters.push_back(k);
    }
  }

  void read_entrances(const Json&) {
    std::array<double, H> median{};
    {
      std::vector<std::array<double, H>> ext;
      for (const auto& x : raw_)
        if (!inside(x.at) && x.kind == "street") {
          std::array<double, H> g{};
          for (std::size_t h = 0; h < H; h++) g[h] = x.weeks->gross[h] / side_share("street");
          ext.push_back(g);
        }
      for (std::size_t h = 0; h < H; h++) {
        std::vector<double> v;
        for (const auto& g : ext) v.push_back(g[h]);
        std::sort(v.begin(), v.end());
        median[h] = v.empty() ? 0 : (v.size() % 2 ? v[v.size() / 2] : 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]));
      }
    }
    const um::Widths widths;
    for (std::size_t j = 0; j < W->junctions(); j++) {
      if (!main[j] || W->leaving[j].size() != 1) continue;
      const auto& p = W->at[j];
      const double dd = std::min({std::fabs(p[0] - box[0]), std::fabs(p[0] - box[2]), std::fabs(p[1] - box[1]), std::fabs(p[1] - box[3])});
      if (dd > o.entrance_edge) continue;
      Entrance E;
      E.node = j;
      const auto& seg = W->segments[W->leaving[j][0] / 2];
      E.street = seg.name, E.kind = seg.kind, E.width = seg.width;
      const Raw* best = nullptr;
      double bd = std::numeric_limits<double>::infinity(), scale = 1;
      for (const auto& x : raw_)
        if (!inside(x.at) && !E.street.empty() && x.street == E.street) {
          const double d = std::hypot(x.at[0] - p[0], x.at[1] - p[1]);
          if (d < bd && d <= o.same_street) bd = d, best = &x;
        }
      if (best) E.source = "same street";
      else if (E.kind != "laneway")
        for (const auto& x : raw_)
          if (!inside(x.at)) {
            const double d = std::hypot(x.at[0] - p[0], x.at[1] - p[1]);
            if (d < bd && d <= o.near_counter) bd = d, best = &x, E.source = "near";
          }
      if (best) {
        const std::string ck = best->kind == "laneway" ? "laneway" : (best->kind == "mall" ? "mall" : "street");
        if (E.source == "near") scale = E.width / um::width_of(widths, ck, best->street);
        for (std::size_t h = 0; h < H; h++) E.gross[h] = scale * best->weeks->gross[h] / side_share(ck);
        E.counter = best->id, E.distance = bd;
      } else {
        E.source = "default";
        for (std::size_t h = 0; h < H; h++) E.gross[h] = median[h] * E.width / widths.main;
      }
      for (double g : E.gross) E.mean += g / H;
      entrances.push_back(std::move(E));
    }
    if (entrances.empty()) throw std::runtime_error("the district has no entrances: no dead end of its network lies on its box");
  }

  void read_masses(const Json& r) {
    const std::size_t n = W->junctions();
    floor.assign(n, 0.0);
    for (const auto& b : r["buildings"].as_array()) {
      const auto& ring = b["ring"].as_array();
      if (ring.size() < 3) continue;
      double area = 0, cx = 0, cy = 0;
      for (std::size_t i = 0; i < ring.size(); i++) {
        const double x0 = ring[i][0].as_number(), y0 = ring[i][1].as_number(), x1 = ring[(i + 1) % ring.size()][0].as_number(), y1 = ring[(i + 1) % ring.size()][1].as_number();
        const double c = x0 * y1 - x1 * y0;
        area += c / 2, cx += (x0 + x1) * c, cy += (y0 + y1) * c;
      }
      if (std::fabs(area) < 1) continue;
      cx /= 6 * area, cy /= 6 * area;
      const double height = num(b, "top", 10) - num(b, "base", 0);
      const double storeys = std::max(1.0, std::round(height / o.storey));
      attach({cx, cy}, std::fabs(area) * storeys, floor);
    }
    for (auto& a : attract) a.assign(n, 0.0);
    std::set<std::tuple<long long, long long, std::string, std::string>> seen;
    const auto& tr = terrain();
    for (const auto& s : r["sensors"].as_array()) {
      const long long id = static_cast<long long>(s["id"].as_number());
      const auto it = tr.by_sensor.find(id);
      if (it == tr.by_sensor.end()) continue;
      const um::Point at{s["at"][0].as_number(), s["at"][1].as_number()};
      for (const auto& e : it->second) {
        const int k = kind_of_use(e.use);
        const um::Point p{at[0] + e.x, at[1] + e.y};
        if (k < 0 || !inside(p)) continue;
        // CLUE places every establishment of a property at its point: one establishment is one (place, name, class).
        if (!seen.insert({std::llround(p[0]), std::llround(p[1]), e.name, e.anzsic}).second) continue;
        attach(p, 1.0, attract[static_cast<std::size_t>(k)]);
        establishments++;
      }
    }
    // The visits' origins: half the entrances by their mean flow, half the buildings by floor.
    origin.assign(n, 0.0);
    double gs = 0, ms = 0;
    for (const auto& E : entrances) gs += E.mean;
    for (double f : floor) ms += f;
    for (const auto& E : entrances) origin[E.node] += 0.5 * E.mean / gs;
    for (std::size_t j = 0; j < n; j++) origin[j] += ms > 0 ? 0.5 * floor[j] / ms : 0;
  }

  void make_streams() {
    const std::size_t n = W->junctions();
    std::array<double, H> gsum{};
    for (const auto& E : entrances)
      for (std::size_t h = 0; h < H; h++) gsum[h] += E.gross[h];
    double gmean = 0;
    for (double g : gsum) gmean += g / H;
    std::vector<std::size_t> buildings;
    for (std::size_t j = 0; j < n; j++)
      if (floor[j] > 0 && junction_target[j] != none) buildings.push_back(j);
    std::vector<std::size_t> g_through, g_arrive, g_leave, g_local;
    for (std::size_t i = 0; i < entrances.size(); i++) {
      const auto& Ei = entrances[i];
      um::Stream th;
      th.name = "through from " + (Ei.street.empty() ? std::string("an entrance") : Ei.street), th.mode = through;
      for (std::size_t j = 0; j < entrances.size(); j++)
        if (j != i && entrances[j].mean > 0) th.origin.push_back(Ei.node), th.target.push_back(junction_target[entrances[j].node]), th.weight.push_back(entrances[j].mean);
      th.rate.resize(H);
      for (std::size_t h = 0; h < H; h++) th.rate[h] = gsum[h] > 0 ? Ei.gross[h] * (gsum[h] - Ei.gross[h]) / gsum[h] : 0;
      if (!th.origin.empty()) g_through.push_back(week.size()), week.push_back(std::move(th));
      um::Stream ar;
      ar.name = "arriving from " + (Ei.street.empty() ? std::string("an entrance") : Ei.street), ar.mode = arrive;
      for (const std::size_t b : buildings) ar.origin.push_back(Ei.node), ar.target.push_back(junction_target[b]), ar.weight.push_back(floor[b]);
      ar.rate.assign(Ei.gross.begin(), Ei.gross.end());
      if (!ar.origin.empty()) g_arrive.push_back(week.size()), week.push_back(std::move(ar));
      um::Stream lv;
      lv.name = "leaving by " + (Ei.street.empty() ? std::string("an entrance") : Ei.street), lv.mode = leave;
      for (const std::size_t b : buildings) lv.origin.push_back(b), lv.target.push_back(junction_target[Ei.node]), lv.weight.push_back(floor[b]);
      lv.rate.assign(Ei.gross.begin(), Ei.gross.end());
      if (!lv.origin.empty()) g_leave.push_back(week.size()), week.push_back(std::move(lv));
    }
    {
      um::Stream lo;
      lo.name = "between buildings", lo.mode = local;
      for (const std::size_t a : buildings)
        for (const std::size_t b : buildings) {
          if (a == b) continue;
          const double w = floor[a] * floor[b] * std::exp(-std::hypot(W->at[a][0] - W->at[b][0], W->at[a][1] - W->at[b][1]) / o.local_length);
          if (w > 0) lo.origin.push_back(a), lo.target.push_back(junction_target[b]), lo.weight.push_back(w);
        }
      lo.rate.assign(H, gmean);
      if (!lo.origin.empty()) g_local.push_back(week.size()), week.push_back(std::move(lo));
    }
    // Visits to each kind's open establishments, and the walks back (their rates set in the calibration).
    for (int k = 0; k < kinds; k++) {
      um::Stream out, back;
      out.name = std::string("visiting ") + kind_name(k), out.mode = visit;
      back.name = std::string("back from ") + kind_name(k), back.mode = visit;
      for (std::size_t a = 0; a < n; a++) {
        const double pull = attract[static_cast<std::size_t>(k)][a];
        if (!(pull > 0) || junction_target[a] == none) continue;
        for (std::size_t q = 0; q < n; q++) {
          if (!(origin[q] > 0) || q == a || junction_target[q] == none) continue;
          out.origin.push_back(q), out.target.push_back(junction_target[a]), out.weight.push_back(origin[q] * pull);
          back.origin.push_back(a), back.target.push_back(junction_target[q]), back.weight.push_back(origin[q] * pull);
        }
      }
      out.rate.assign(H, 0.0), back.rate.assign(H, 0.0);
      if (!out.origin.empty()) visit_streams_[static_cast<std::size_t>(k)] = {week.size(), week.size() + 1}, week.push_back(std::move(out)), week.push_back(std::move(back));
    }
    groups = {g_through, g_arrive, g_leave, g_local};
  }

  std::array<std::array<std::size_t, 2>, kinds> visit_streams_{{{none, none}, {none, none}, {none, none}, {none, none}}};

  void calibrate() {
    const auto t0 = std::chrono::steady_clock::now();
    const std::size_t K = counters.size(), S = week.size();
    if (K == 0) throw std::runtime_error("no counter in the district counts each way over the reference year");
    std::vector<um::UnitCounts> unit(S);
    calx::scheduler::parallel_for(pool(), S, 1, [&](std::size_t begin, std::size_t end, std::size_t) {
      for (std::size_t s = begin; s < end; s++) unit[s] = um::unit_counts(*V, *T, week[s], counters);
    });
    std::vector<std::vector<std::array<double, 2>>> observed(K, std::vector<std::array<double, 2>>(H));
    for (std::size_t k = 0; k < K; k++)
      for (std::size_t h = 0; h < H; h++) observed[k][h] = {info[k].way[0][h], info[k].way[1][h]};
    std::vector<char> grouped(S, 0);
    for (const auto& g : groups)
      for (const std::size_t s : g) grouped[s] = 1;
    // The visits' open weight by kind and hour.
    std::array<std::array<double, H>, kinds> open{};
    for (int k = 0; k < kinds; k++) {
      double n = 0;
      for (double a : attract[static_cast<std::size_t>(k)]) n += a;
      for (std::size_t h = 0; h < H; h++) open[static_cast<std::size_t>(k)][h] = open_at(k, h) ? attraction[static_cast<std::size_t>(k)] * n : 0;
    }
    um::Calibration C;
    for (int round = 0; round < 3; round++) {
      for (int vround = 0; vround < 4; vround++) {
        C = um::calibrate(week, unit, counters, observed, groups);
        // The visits: ν of the hour's trips, out and back, among the kinds open.
        for (std::size_t h = 0; h < H; h++) {
          double others = 0, total_open = 0;
          for (std::size_t g = 0; g < groups.size(); g++)
            for (const std::size_t s : groups[g]) others += C.amplitude[g][h] * week[s].rate[h];
          for (int k = 0; k < kinds; k++) total_open += open[static_cast<std::size_t>(k)][h];
          const double visits = o.visit_share / (1 - o.visit_share) * others;
          for (int k = 0; k < kinds; k++) {
            const auto [a, b] = visit_streams_[static_cast<std::size_t>(k)];
            if (a == none) continue;
            const double r = total_open > 0 ? 0.5 * visits * open[static_cast<std::size_t>(k)][h] / total_open : 0;
            week[a].rate[h] = r, week[b].rate[h] = r;
          }
        }
      }
      C = um::calibrate(week, unit, counters, observed, groups);
      // Each footpath's handedness, by least squares over the week given the amplitudes.
      for (std::size_t k = 0; k < K; k++) {
        double num = 0, den = 0;
        for (std::size_t h = 0; h < H; h++) {
          const double p0 = C.fitted[k][h][0] / (1 + counters[k].hand[0]), p1 = C.fitted[k][h][1] / (1 + counters[k].hand[1]);
          num += (observed[k][h][0] - p0) * p0 - (observed[k][h][1] - p1) * p1, den += p0 * p0 + p1 * p1;
        }
        const double rho = den > 0 ? std::clamp(num / den, -o.hand_limit, o.hand_limit) : 0;
        counters[k].hand = {rho, -rho};
      }
    }
    C = um::calibrate(week, unit, counters, observed, groups);
    fit.amplitude = C.amplitude, fit.r2_way = C.r2, fit.rows = C.rows;
    fit.fitted.assign(K, {});
    double ss = 0, st = 0, mean = 0, cnt = 0;
    for (std::size_t k = 0; k < K; k++)
      for (std::size_t h = 0; h < H; h++) {
        fit.fitted[k][0][h] = C.fitted[k][h][0], fit.fitted[k][1][h] = C.fitted[k][h][1];
        mean += observed[k][h][0] + observed[k][h][1], cnt++;
      }
    mean /= std::max(1.0, cnt);
    fit.counter_r2.assign(K, 0.0);
    for (std::size_t k = 0; k < K; k++) {
      double a = 0, b = 0, m = 0;
      for (std::size_t h = 0; h < H; h++) m += (observed[k][h][0] + observed[k][h][1]) / H;
      for (std::size_t h = 0; h < H; h++) {
        const double y = observed[k][h][0] + observed[k][h][1], p = C.fitted[k][h][0] + C.fitted[k][h][1];
        ss += (y - p) * (y - p), st += (y - mean) * (y - mean), a += (y - p) * (y - p), b += (y - m) * (y - m);
      }
      fit.counter_r2[k] = b > 0 ? 1 - a / b : 0;
    }
    fit.r2_gross = st > 0 ? 1 - ss / st : 0;
    week = um::calibrated(week, C, groups);
    double vt = 0, at = 0;
    for (const auto& s : week)
      for (double r : s.rate) (s.mode == visit ? vt : at) += r;
    fit.visit_share = vt / std::max(1e-9, vt + at);
    fit.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  }
};

// The prepared districts, by the key of what the page sent (at most a few kept).
inline std::shared_ptr<const District> district_of(const Json& r) {
  static std::map<std::string, std::shared_ptr<const District>> cache;
  static std::vector<std::string> order;
  static std::mutex lock;
  // The key: FNV-1a of what defines the district.
  std::uint64_t h = 1469598103934665603ull;
  const auto mix = [&](const std::string& s) { for (const unsigned char c : s) h = (h ^ c) * 1099511628211ull; };
  // A request may name a district already prepared by its key alone.
  if (r["district"].is_string() && !r["junctions"].is_array()) {
    std::lock_guard<std::mutex> g(lock);
    if (auto it = cache.find(r["district"].text()); it != cache.end()) return it->second;
    throw std::invalid_argument("unknown district: send it whole");
  }
  for (const char* k : {"box", "junctions", "segments", "buildings", "sensors"}) mix(r[k].dump());
  const std::string key = std::to_string(h);
  {
    std::lock_guard<std::mutex> g(lock);
    if (auto it = cache.find(key); it != cache.end()) return it->second;
  }
  if (!r["junctions"].is_array() || !r["segments"].is_array() || !r["box"].is_object()) throw std::invalid_argument("a district: its box, junctions and segments, its buildings and sensors");
  auto D = std::make_shared<const District>(r, key);
  std::lock_guard<std::mutex> g(lock);
  cache[key] = D, order.push_back(key);
  if (order.size() > 3) cache.erase(order.front()), order.erase(order.begin());
  return D;
}

}  // namespace hf
