// A district's economy: the problem three case studies share (README.md). They differ only in how much structure
// their model of the city's people has; this header holds everything else, so that the three evaluate one value of
// one decision on one district.
//
// The district. The blocks between Swanston and Russell Streets, Flinders Lane and Bourke Street, Melbourne's Hoddle
// Grid: the case study "Designing a district" (research/experiments/2026-11-05-designing-a-district), its walking
// network, its 61 parcels in the same order (along the area's long axis) and its eleven candidate links (eight blind
// lanes carried through, three passages through blocks). Each parcel has one ground-floor frontage: its facade on the
// street it fronts most, read from the urban-science domain's city network (urban::events::city_network, on
// architecture's massing of the parcels and their 524 neighbours as built), which also gives the walls the buildings
// share (party walls). The frontages' relations are the frontage graph: the party walls, the frontages that face one
// street segment, and the frontages within reach of one another along the network.
//
// The decision (Plan). For each frontage f: its kind k_f ∈ {vacant, retail, food and drink, services, office lobby,
// culture}; its hours on a weekday and at the weekend, each an integer start s ∈ {0, …, 23} and a length
// ℓ ∈ {0, …, 24} (closed when ℓ = 0; open over [s, min(24, s + ℓ)) civil hours, every weekday alike and both weekend
// days alike); and for each candidate link whether it is made. Every value of the encoding is a valid plan.
//
// The value (README §4), A$ per week of the mean week, given the people passing each frontage at each hour of the
// week, P_f(h) (the study's model of flow supplies them):
//
//   V = Σ_f Σ_h [open] P_f(h) p_f(h) c_f(h) s_k m_k − Σ_f (staffing + outgoings + refit) − Σ_links (their cost)
//
//   p_f(h) = π_k φ_k(band(h)) · clamp(1 + γ Σ_{g open, g ≠ f, k_g ≠ vacant} w_fg a(k_g → k_f), ½, 2): the chance that one
//            passer stops, the kind's propensity times its multiplier for the day, the evening or the night, times
//            the catalysis of the open neighbours; a(i → j) = ln CLQ_{i→j} of CLUE's establishments, clamped to
//            [−1, 1]; w_fg = exp(−d_fg/ξ) for walking distances d_fg ≤ 3ξ, ξ the coupling length of the domain's
//            catalysis (urban::affinity::walk_coupling_length);
//   c_f(h) = [1 − Π_{g ∈ G}(1 − p_g(h))] / Σ_{g ∈ G} p_g(h): competition, G the open frontages of f's kind on f's
//            segment, each passer an independent chance at each, the stops shared in proportion to the chances
//            (c = 1 when f is alone);
//   s_k, m_k the spend of a visit and its gross margin; staffing = staff_k × wage × (open hours + one hour for each day
//            opened); outgoings per m² of trading floor; refit, a change of kind (to any but vacant) at so much per m²
//            of trading floor, amortised; a link's cost per metre, amortised.
//
// The other measures: the spending's spread over the week (Parseval on the week's 168 hours: the mean line's share
// of the power, 1/(1 + CV²); the daily lines' share; the week's other lines' share), its evening and weekend shares,
// the mix H(K | block) of the frontages' kinds (vacant left out) weighted by their length
// (physis::statistical_mechanics::mixing; the urban-science domain's Definition 2.6), and access: the walk from where
// the people are to the nearest open frontage of each daily need (retail, food and drink, services), averaged over
// the needs and over the week's hours, people-weighted, capped where none is open.
//
// What the study supplies: P_f(h) for the plan, from its model of flow. Everything here is the same for all three.
//
// Reused from the designing study as it is (its own provisional code, gaps it reported in its README §8.2):
// lanes.hpp (blind lanes carried through) and walking.hpp (Dijkstra's shortest walks).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numbers>
#include <numeric>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <architecture/space/massing.hpp>
#include <calx/linalg/rectangular.hpp>
#include <calx/random/stream.hpp>
#include <heuresis/variation/binary.hpp>
#include <heuresis/variation/integer.hpp>
#include <physis/computational_geometry/planar.hpp>
#include <physis/statistical_mechanics/mixing.hpp>
#include <science/json.hpp>
#include <urban/affinity/catalysis.hpp>
#include <urban/events/network.hpp>
#include <urban/landscape/district.hpp>

#include "../2026-11-05-designing-a-district/lanes.hpp"
#include "../2026-11-05-designing-a-district/walking.hpp"

namespace economy {

namespace ul = urban::landscape;
namespace as = architecture::space;
namespace ue = urban::events;
namespace hv = heuresis::variation;
namespace cg = physis::computational_geometry;
using science::Json;
using Point = study::Point;
using Ring = study::Ring;

// ---- kinds, hours, plans --------------------------------------------------------------------------------------------

inline constexpr std::size_t kinds = 6;
enum Kind : std::int64_t { vacant = 0, retail = 1, food = 2, services = 3, lobby = 4, culture = 5 };
inline const char* kind_name(std::int64_t k) {
  static const char* names[] = {"vacant", "retail", "food and drink", "services", "office lobby", "culture"};
  return names[std::clamp<std::int64_t>(k, 0, 5)];
}
inline bool trades(std::int64_t k) { return k == retail || k == food || k == services || k == culture; }
inline constexpr std::array<std::int64_t, 3> needs{retail, food, services};  // the daily needs of access

inline constexpr std::size_t week = 168;  // hours of the week, Monday 00:00 civil first
inline bool is_weekend(std::size_t h) { return h / 24 >= 5; }
inline std::size_t hour_of_day(std::size_t h) { return h % 24; }
// The bands of the day: 0 day (07–18), 1 evening (18–24), 2 night (00–07).
inline std::size_t band(std::size_t h) {
  const std::size_t t = hour_of_day(h);
  return t >= 18 ? 1 : (t >= 7 ? 0 : 2);
}

inline constexpr std::size_t hour_genes = 4;  // weekday start, weekday length, weekend start, weekend length
inline std::int64_t close_of(std::int64_t start, std::int64_t length) { return std::min<std::int64_t>(24, start + length); }
inline bool open_in(std::int64_t start, std::int64_t length, std::size_t t) {
  return length > 0 && static_cast<std::int64_t>(t) >= start && static_cast<std::int64_t>(t) < close_of(start, length);
}

struct Plan {
  hv::Integers kind;   // per frontage, a Kind
  hv::Integers hours;  // hour_genes per frontage
  hv::Bits links;      // per candidate link
  bool operator==(const Plan&) const = default;
  bool open(std::size_t f, std::size_t h) const {
    if (kind[f] == vacant) return false;
    const std::size_t o = hour_genes * f + (is_weekend(h) ? 2 : 0);
    return open_in(hours[o], hours[o + 1], hour_of_day(h));
  }
  double span(std::size_t f, bool weekend) const {
    const std::size_t o = hour_genes * f + (weekend ? 2 : 0);
    return hours[o + 1] > 0 ? static_cast<double>(close_of(hours[o], hours[o + 1]) - hours[o]) : 0.0;
  }
  double weekly_hours(std::size_t f) const { return kind[f] == vacant ? 0.0 : 5 * span(f, false) + 2 * span(f, true); }
  double days_open(std::size_t f) const { return kind[f] == vacant ? 0.0 : 5.0 * (span(f, false) > 0) + 2.0 * (span(f, true) > 0); }
};

// ---- parameters (README §4.3: each sourced or declared an assumption) --------------------------------------------------

inline double annuity_weekly(double capital, double years, double rate) { return capital * rate / (1 - std::pow(1 + rate, -years)) / (365.25 / 7); }

struct Parameters {
  double wage = 35;                                                         // A$ an hour a staff member, with on-costs (assumption)
  std::array<double, kinds> staff{0, 2, 3, 2, 0, 2};                        // staff while open (assumption)
  std::array<double, kinds> propensity{0, 0.006, 0.010, 0.0025, 0, 0.0015};  // π_k: a passer's chance to stop by day (assumption)
  std::array<std::array<double, 3>, kinds> period{{{0, 0, 0}, {1, 0.7, 0.3}, {1, 1.3, 0.6}, {1, 0.4, 0.1}, {0, 0, 0}, {1, 2.0, 0.5}}};  // φ_k (assumption)
  std::array<double, kinds> spend{0, 45, 18, 55, 0, 30};                    // s_k, A$ a visit (assumption)
  std::array<double, kinds> margin{0, 0.40, 0.65, 0.70, 0, 0.75};           // m_k, the share of spending kept (assumption)
  std::array<double, kinds> outgoings{0, 2, 2, 2, 0, 2};                    // A$ a week per m² of trading floor (assumption)
  double refit = 1500;                                                      // A$ per m² of trading floor for a change of kind (assumption)
  double refit_years = 10, link_years = 30, rate = 0.05;                    // amortisation (assumption)
  double link = 20000;                                                      // A$ a metre of link, an arcade through a building (assumption)
  double depth = 20;                                                        // m: trading floor at most its length × this (assumption)
  double catalysis = 0.5;                                                   // γ (assumption)
  // ξ: the coupling length of the urban-science domain's catalysis at its defaults (D = 1 m²/s, ℓ = 1/1800 s⁻¹,
  // Δs = 10 m; urban/affinity/catalysis.hpp), 42.5 m: an assumption of that representation, taken as it stands.
  double coupling = urban::affinity::walk_coupling_length(1.0, 1.0 / 1800, 10);
  double catalysis_floor = 0.5, catalysis_ceiling = 2.0;                    // the clamp of 1 + κ (assumption)
  double access_cap = 400;                                                  // m: the walk when no frontage of a need is open (assumption)
  std::array<std::array<double, kinds>, kinds> affinity{};                  // a(i → j) = clamp(ln CLQ_{i→j}, −1, 1) (CLUE 2024; calibrate.cpp)

  double refit_weekly(double m2) const { return annuity_weekly(refit * m2, refit_years, rate); }
  double link_weekly(double metres) const { return annuity_weekly(link * metres, link_years, rate); }
};

// ---- the district ---------------------------------------------------------------------------------------------------

struct Frontage {
  std::string parcel, name, clue_property, street, way;  // way: the kind of its segment (street, laneway, mall)
  Ring footprint;
  double footprint_area = 0, height = 0;
  Point centre{}, anchor{};      // the parcel's centroid; the middle of its main facade
  std::uint64_t segment = 0;     // the base district's 1-cell it fronts
  double length = 0, facade = 0; // m: its main facade's length, and all its facades on streets
  double trading_area = 0;       // m²
  bool fronts = false;           // false: no facade faces a street (its anchor then the parcel's centre)
  int block = -1;                // the place it stands in (blocks, then holes)
  double along = 0;              // position along the area's long axis (the encoding's order)
  std::int64_t kind0 = vacant;
  std::array<std::int64_t, 4> hours0{};
  std::string read0;             // how the kind as built was read: joined (CLUE's base property), inside (CLUE's points), none
};

struct Link {
  std::string kind;  // "lane carried through" or "passage"
  Point a{}, b{};
  double length = 0;
  std::string from, to;
  std::size_t through = 0;
};

// Walking distances between the frontages with a set of links made, and each frontage's neighbours for catalysis.
struct Distances {
  std::size_t n = 0;
  std::vector<double> d;                                          // n × n, m
  std::vector<std::vector<std::pair<std::size_t, double>>> near;  // (g, w_fg) for d_fg ≤ 3ξ, g ≠ f
  double at(std::size_t f, std::size_t g) const { return d[f * n + g]; }
};

// What a plan does under given flows.
struct Outcome {
  double value = 0, spending = 0, margin = 0, visits = 0;          // A$/week, A$/week, A$/week, visits/week
  double staffing = 0, outgoings = 0, refits = 0, link_costs = 0;  // A$/week
  double footfall = 0;                                             // people passing frontages, summed over the week
  double spread = 0, effective_hours = 0;                          // P₀/ΣP of the spending's 168 hours; 168 × that
  double daily_share = 0, weekly_share = 0;                        // power of the daily lines (m = 7j) and of the others, of ΣP
  double footfall_spread = 0;
  double evening_share = 0, weekend_share = 0;                     // of the spending: 18–24 h; Saturday and Sunday
  double mix = 0, mix_entropy = 0;                                 // nats: H(K | block), H(K)
  double access = 0;                                               // m
  std::array<double, 3> access_by_need{};
  std::size_t open_frontages = 0, links = 0, changed = 0;
  std::array<double, kinds> spending_by_kind{}, frontages_by_kind{};
  std::array<double, week> spending_by_hour{}, visits_by_hour{}, footfall_by_hour{};
  std::vector<double> spending_by_frontage, value_by_frontage, visits_by_frontage;
  std::vector<double> power;                                       // |Ŝ(m)|², m = 0…84, folded (m and 168 − m)
};

class Problem {
 public:
  Parameters P;
  std::vector<Frontage> frontages;
  std::vector<Link> links;
  std::vector<study::Through> throughs;
  ul::Network plane;
  ul::District base;
  std::vector<std::string> dropped;
  std::vector<std::uint64_t> places;
  std::size_t place_count = 0;
  Ring quad;
  // The relations of the buildings (README §2.3).
  std::vector<std::array<std::size_t, 2>> party;  // frontage pairs whose parcels share a wall
  std::vector<double> party_area;                 // m², each pair's
  std::size_t party_with_context = 0;             // walls a parcel shares with a building of the context
  std::size_t context_buildings = 0, massed = 0;
  std::vector<std::vector<std::size_t>> on_segment;  // per frontage: the other frontages on its segment
  std::vector<std::pair<std::uint64_t, std::vector<std::size_t>>> segments;  // the segments fronted, each with its frontages

  explicit Problem(const Json& inputs, const Json* problem = nullptr);

  std::size_t size() const { return frontages.size(); }
  Plan as_built() const {
    Plan p;
    for (const auto& f : frontages) {
      p.kind.push_back(f.kind0);
      for (auto h : f.hours0) p.hours.push_back(h);
    }
    p.links.assign(links.size(), 0);
    return p;
  }
  hv::IntegerBox kind_box() const { return hv::IntegerBox::cube(size(), 0, kinds - 1); }
  hv::IntegerBox hour_box() const {
    hv::IntegerBox b;
    for (std::size_t f = 0; f < size(); f++)
      for (std::size_t j = 0; j < hour_genes; j++) b.lower.push_back(0), b.upper.push_back(j % 2 == 0 ? 23 : 24);
    return b;
  }
  void validate(const Plan& p) const {
    if (p.kind.size() != size() || p.hours.size() != hour_genes * size() || p.links.size() != links.size()) throw std::invalid_argument("economy::Problem: a plan of another district");
    for (auto k : p.kind)
      if (k < 0 || k >= static_cast<std::int64_t>(kinds)) throw std::invalid_argument("economy::Problem: a kind out of range");
    for (std::size_t i = 0; i < p.hours.size(); i++)
      if (p.hours[i] < 0 || p.hours[i] > (i % 2 == 0 ? 23 : 24)) throw std::invalid_argument("economy::Problem: an hour out of range");
  }

  // The walking distances with the plan's links made (cached by the set of links).
  std::shared_ptr<const Distances> distances(const hv::Bits& made) const;
  ul::District with_links(const hv::Bits& made) const;

  // The value and the measures of a plan given the people passing each frontage at each hour, passers[f * 168 + h].
  Outcome evaluate(const Plan& plan, std::span<const double> passers) const;

  // The canonical relaxation (README §3.4): a point of [0, 1]^dimension() to a plan, and a plan to its cell's centre.
  std::size_t dimension() const { return size() * (kinds + hour_genes) + links.size(); }
  Plan decode(std::span<const double> u) const;
  std::vector<double> encode(const Plan& p) const;

 private:
  mutable std::mutex mutex_;
  mutable std::map<std::uint64_t, std::shared_ptr<const Distances>> cache_;
};

// ---- construction -----------------------------------------------------------------------------------------------------

namespace detail {
inline Ring ring_of(const Json& j) {
  Ring r;
  for (const auto& p : j.as_array()) r.push_back({p[0].as_number(), p[1].as_number()});
  return r;
}
inline Point centroid(const Ring& r) {
  double a = 0, cx = 0, cy = 0;
  for (std::size_t k = 0; k < r.size(); k++) {
    const auto& p = r[k];
    const auto& q = r[(k + 1) % r.size()];
    const double w = p[0] * q[1] - q[0] * p[1];
    a += w, cx += (p[0] + q[0]) * w, cy += (p[1] + q[1]) * w;
  }
  return std::fabs(a) > 1e-9 ? Point{cx / (3 * a), cy / (3 * a)} : r[0];
}
inline Ring counter_clockwise(Ring r) {
  if (cg::ring_area(r) < 0) std::reverse(r.begin(), r.end());
  return r;
}
}  // namespace detail

inline Problem::Problem(const Json& in, const Json* problem) {
  const Json& D = in["district"];
  quad = detail::ring_of(D["redesign_quad"]);
  // The area's long axis, as the designing study takes it: from its southern side's middle to its northern side's.
  const Point s0{(quad[0][0] + quad[1][0]) / 2, (quad[0][1] + quad[1][1]) / 2}, s1{(quad[2][0] + quad[3][0]) / 2, (quad[2][1] + quad[3][1]) / 2};
  const double L = std::hypot(s1[0] - s0[0], s1[1] - s0[1]);
  const Point axis{(s1[0] - s0[0]) / L, (s1[1] - s0[1]) / L};

  // Kinds and typical hours as built.
  std::map<std::string, std::int64_t> kind_of_parcel;
  std::map<std::string, std::string> read_of_parcel;
  for (const auto& p : in["parcels"].as_array()) kind_of_parcel[p["id"].text()] = static_cast<std::int64_t>(p["kind"].as_number()), read_of_parcel[p["id"].text()] = p["read"].text();
  std::array<std::array<std::int64_t, 4>, kinds> typical{};
  for (std::size_t k = 1; k < kinds; k++) {
    const Json& t = in["typical_hours"][kind_name(static_cast<std::int64_t>(k))];
    const auto pair = [&](const Json& w, std::int64_t& s, std::int64_t& l) {
      if (w.is_null()) { s = 0, l = 0; return; }
      s = static_cast<std::int64_t>(w[0].as_number()), l = std::max<std::int64_t>(0, static_cast<std::int64_t>(w[1].as_number()) - s);
    };
    pair(t["weekday"], typical[k][0], typical[k][1]);
    pair(t["weekend"], typical[k][2], typical[k][3]);
  }

  for (const auto& j : D["parcels"].as_array()) {
    Frontage f;
    f.parcel = j["id"].text(), f.name = j["name"].text(), f.clue_property = j["clue_property"].text();
    f.footprint = detail::counter_clockwise(detail::ring_of(j["footprint"]));
    f.footprint_area = cg::ring_area(f.footprint);
    f.height = j["height_m"].as_number();
    f.centre = detail::centroid(f.footprint);
    f.along = (f.centre[0] - s0[0]) * axis[0] + (f.centre[1] - s0[1]) * axis[1];
    f.kind0 = kind_of_parcel.count(f.parcel) ? kind_of_parcel.at(f.parcel) : vacant;
    f.read0 = read_of_parcel.count(f.parcel) ? read_of_parcel.at(f.parcel) : "none";
    if (f.kind0 != vacant) f.hours0 = typical[static_cast<std::size_t>(f.kind0)];
    frontages.push_back(std::move(f));
  }
  std::stable_sort(frontages.begin(), frontages.end(), [](const Frontage& a, const Frontage& b) { return a.along < b.along; });

  // The district as the urban-science domain builds it (planarised, as the workbench's engine does).
  ul::Network N;
  for (const auto& p : D["network"]["junctions"].as_array()) N.junctions.push_back({p[0].as_number(), p[1].as_number()});
  for (const auto& s : D["network"]["segments"].as_array()) {
    ul::Segment g;
    g.u = static_cast<std::size_t>(s["u"].as_number()), g.v = static_cast<std::size_t>(s["v"].as_number());
    for (const auto& p : s["line"].as_array()) g.line.push_back({p[0].as_number(), p[1].as_number()});
    g.kind = s["kind"].text("street"), g.name = s["name"].text();
    N.segments.push_back(std::move(g));
  }
  plane = ul::planarised(std::move(N), &dropped);
  base = ul::district(plane);
  for (const auto id : base.of(2)) places.push_back(id);
  place_count = places.size() + base.hole_rings.size();
  for (auto& f : frontages) {
    for (std::size_t k = 0; k < places.size() && f.block < 0; k++)
      if (cg::inside_ring(base.ring(places[k]), f.centre)) f.block = static_cast<int>(k);
    for (std::size_t k = 0; k < base.hole_rings.size() && f.block < 0; k++)
      if (cg::inside_ring(base.hole_rings[k], f.centre)) f.block = static_cast<int>(places.size() + k);
    if (f.block < 0) f.block = static_cast<int>(place_count++);
  }

  // The candidate links, exactly as the designing study makes them (its district.hpp, Problem's constructor): every
  // blind lane of the area carried through (at least 15 m, at most 120 m), then one passage across the middle of each
  // block (a 2-cell of at least 1,000 m² inside the area) between its long sides, by the domain's plan_passage.
  throughs = study::lanes_through(plane, quad, 15, 120);
  for (std::size_t k = 0; k < throughs.size(); k++) links.push_back({"lane carried through", throughs[k].from, throughs[k].to, throughs[k].length, throughs[k].lane, throughs[k].street, k});
  for (const auto id : base.of(2)) {
    const auto r = base.ring(id);
    const Point c = detail::centroid(r);
    if (!cg::inside_ring(quad, c) || !(base.area(id) >= 1000)) continue;
    double sxx = 0, syy = 0, sxy = 0;
    for (std::size_t k = 0; k < r.size(); k++) {
      const auto& a = r[k];
      const auto& b = r[(k + 1) % r.size()];
      const double l = std::hypot(b[0] - a[0], b[1] - a[1]);
      for (double t = 0; t < l; t += 1) {
        const double x = a[0] + (b[0] - a[0]) * t / l - c[0], y = a[1] + (b[1] - a[1]) * t / l - c[1];
        sxx += x * x, syy += y * y, sxy += x * y;
      }
    }
    const double theta = 0.5 * std::atan2(2 * sxy, sxx - syy);
    const Point across{-std::sin(theta), std::cos(theta)};
    const Point a{c[0] + 500 * across[0], c[1] + 500 * across[1]}, b{c[0] - 500 * across[0], c[1] - 500 * across[1]};
    const auto plan = ul::plan_passage(base, id, a, b);
    if (!plan.refused.empty()) continue;
    const double len = std::hypot(plan.a[0] - plan.b[0], plan.a[1] - plan.b[1]);
    if (len < 15) continue;
    links.push_back({"passage", plan.a, plan.b, len, base.name.at(plan.edge_a), base.name.at(plan.edge_b), 0});
  }

  // The buildings as one complex (architecture's massing): the parcels and their context as built; the domain's city
  // network gives their party walls and the facades that front a street (urban::events::city_network).
  std::vector<as::Footprint> F;
  for (const auto& f : frontages) F.push_back({f.parcel, f.name, f.footprint, 0.0, f.height});
  for (const auto& j : D["context"].as_array()) {
    F.push_back({j["id"].text(), j["name"].text(), detail::counter_clockwise(detail::ring_of(j["footprint"])), 0.0, j["height_m"].as_number()});
    context_buildings++;
  }
  const auto M = as::massing(F);
  massed = M.building.size();
  const auto C = ue::city_network(M, F, base);
  std::map<std::string, std::size_t> index;
  for (std::size_t i = 0; i < frontages.size(); i++) index[frontages[i].parcel] = i;
  std::vector<std::size_t> coupling_frontage(C.couplings.size(), std::numeric_limits<std::size_t>::max());
  for (const auto& [key, c] : C.building)
    if (index.count(key)) coupling_frontage[c] = index.at(key);
  std::map<std::size_t, std::uint64_t> street_cell;  // a street coupling ↦ its 1-cell
  for (const auto& [cell, c] : C.street) street_cell[c] = cell;
  // Each frontage's main facade: the street (or mall) it fronts most, if it fronts one by at least 3 m; else the
  // laneway it fronts most. A shop faces the street; a building's back on a service lane is not its frontage.
  struct Face { double area = 0, length = 0; std::uint64_t segment = 0; Point anchor{}; };
  std::vector<Face> best_street(frontages.size()), best_lane(frontages.size());
  for (std::size_t q = 0; q < C.primitives.size(); q++) {
    const auto& p = C.primitives[q];
    if (p.type == "party-wall") {
      const auto fa = coupling_frontage[p.a], fb = coupling_frontage[p.b];
      const bool a_in = fa != std::numeric_limits<std::size_t>::max(), b_in = fb != std::numeric_limits<std::size_t>::max();
      if (a_in && b_in) party.push_back({std::min(fa, fb), std::max(fa, fb)}), party_area.push_back(p.area);
      else if (a_in || b_in) party_with_context++;
    } else if (p.type == "frontage") {
      const auto f = coupling_frontage[p.a];
      if (f == std::numeric_limits<std::size_t>::max()) continue;
      auto& x = frontages[f];
      const double length = p.area / std::max(1e-9, x.height);
      x.facade += length;
      const auto seg = street_cell.at(p.b);
      auto& slot = base.cell(seg).type == "laneway" ? best_lane[f] : best_street[f];
      if (p.area > slot.area) {
        // The facade's middle: its wall pieces' base lines, each weighted by its length.
        double w = 0, mx = 0, my = 0;
        const auto& sh = C.shape[q];
        for (std::size_t k = 0; k + 1 < sh.size(); k += 2) {
          const double l = std::hypot(sh[k + 1][0] - sh[k][0], sh[k + 1][1] - sh[k][1]);
          w += l, mx += l * (sh[k][0] + sh[k + 1][0]) / 2, my += l * (sh[k][1] + sh[k + 1][1]) / 2;
        }
        slot = {p.area, length, seg, w > 0 ? Point{mx / w, my / w} : x.centre};
      }
    }
  }
  for (std::size_t f = 0; f < frontages.size(); f++) {
    const Face& use = best_street[f].length >= 3 ? best_street[f] : best_lane[f];
    if (!(use.area > 0)) continue;
    auto& x = frontages[f];
    x.fronts = true, x.length = use.length, x.segment = use.segment, x.anchor = use.anchor;
  }
  for (auto& f : frontages) {
    if (!f.fronts) {
      // No facade faces a street: the frontage sits on the segment nearest its centre, with no length.
      ul::District copy = base;
      const auto k = ul::attach_store(copy, {f.parcel, "", f.centre});
      f.segment = copy.stores[k].segment, f.anchor = f.centre;
    }
    f.street = base.name.count(f.segment) ? base.name.at(f.segment) : "";
    f.way = base.cell(f.segment).type;
    f.trading_area = std::min(f.footprint_area, f.length * P.depth);
  }
  std::map<std::uint64_t, std::vector<std::size_t>> by_segment;
  for (std::size_t i = 0; i < frontages.size(); i++) by_segment[frontages[i].segment].push_back(i);
  on_segment.assign(frontages.size(), {});
  for (const auto& [seg, fs] : by_segment) {
    segments.push_back({seg, fs});
    for (auto f : fs)
      for (auto g : fs)
        if (f != g) on_segment[f].push_back(g);
  }

  // The parameters a calibration made: the affinities of the kinds.
  if (problem && problem->has("affinity")) {
    const Json& A = (*problem)["affinity"]["a"];
    for (std::size_t i = 0; i < kinds; i++)
      for (std::size_t j = 0; j < kinds; j++) P.affinity[i][j] = A[i][j].is_null() ? 0.0 : A[i][j].as_number();
  }
}

// ---- the walking network with links --------------------------------------------------------------------------------

inline ul::District Problem::with_links(const hv::Bits& made) const {
  std::vector<study::Through> chosen;
  for (std::size_t k = 0; k < links.size(); k++)
    if (made[k] && links[k].kind == "lane carried through") chosen.push_back(throughs[links[k].through]);
  ul::District D = chosen.empty() ? base : ul::district(study::carry_through(plane, chosen, 120));
  for (std::size_t k = 0; k < links.size(); k++) {
    if (!made[k] || links[k].kind != "passage") continue;
    const Point mid{(links[k].a[0] + links[k].b[0]) / 2, (links[k].a[1] + links[k].b[1]) / 2};
    std::uint64_t block = 0;
    for (const auto id : D.of(2))
      if (cg::inside_ring(D.ring(id), mid)) { block = id; break; }
    if (!block) throw std::invalid_argument("economy::Problem::with_links: no block holds the passage between " + links[k].from + " and " + links[k].to);
    ul::open_passage(D, block, links[k].a, links[k].b);
  }
  return D;
}

inline std::shared_ptr<const Distances> Problem::distances(const hv::Bits& made) const {
  std::uint64_t key = 0;
  for (std::size_t k = 0; k < made.size(); k++)
    if (made[k]) key |= std::uint64_t{1} << k;
  {
    std::lock_guard<std::mutex> lock(mutex_);
    const auto it = cache_.find(key);
    if (it != cache_.end()) return it->second;
  }
  ul::District D = with_links(made);
  study::WalkGraph G;
  std::map<std::size_t, std::size_t> node;
  std::map<std::uint64_t, std::size_t> edge;
  for (std::size_t i = 0; i < D.complex.size(); i++)
    if (D.complex.cells[i].dim == 0) node[i] = G.nodes++;
  for (std::size_t i = 0; i < D.complex.size(); i++) {
    const auto& c = D.complex.cells[i];
    if (c.dim != 1) continue;
    const auto ends = D.complex.ends(i);
    edge[c.id] = G.edges.size();
    G.edges.push_back({node.at(ends[0]), node.at(ends[1]), D.length(c.id)});
  }
  G.finish();
  // Each frontage on the nearest segment of the streets that were there: a link is a new way past the frontages, not
  // a frontage of theirs (the domain's attach_store would take whichever segment is nearest, the new one included).
  const auto made_by_link = [&](std::uint64_t id) {
    const auto it = D.name.find(id);
    if (it == D.name.end()) return false;
    const std::string& nm = it->second;
    const std::string tail = " carried through";
    return nm == "a passage" || (nm.size() >= tail.size() && nm.compare(nm.size() - tail.size(), tail.size(), tail) == 0);
  };
  std::vector<study::OnSegment> at;
  for (const auto& f : frontages) {
    double best = std::numeric_limits<double>::infinity(), s = 0;
    std::uint64_t seg = 0;
    for (const auto& [id, l] : D.line) {
      if (made_by_link(id)) continue;
      const auto [along, dist] = cg::nearest_on_polyline(l, f.anchor);
      if (dist < best) best = dist, seg = id, s = along;
    }
    at.push_back({edge.at(seg), s});
  }
  auto out = std::make_shared<Distances>();
  const std::size_t n = frontages.size();
  out->n = n, out->d.assign(n * n, 0.0), out->near.assign(n, {});
  for (std::size_t f = 0; f < n; f++) {
    const std::vector<study::OnSegment> src{at[f]};
    const auto d = study::from_sources(G, src);
    for (std::size_t g = 0; g < n; g++) out->d[f * n + g] = f == g ? 0.0 : study::to_target(G, d, src, at[g]);
  }
  for (std::size_t f = 0; f < n; f++)
    for (std::size_t g = 0; g < f; g++) {  // symmetric to rounding: the mean of the two directions
      const double m = 0.5 * (out->d[f * n + g] + out->d[g * n + f]);
      out->d[f * n + g] = out->d[g * n + f] = m;
    }
  for (std::size_t f = 0; f < n; f++)
    for (std::size_t g = 0; g < n; g++)
      if (g != f && out->d[f * n + g] <= 3 * P.coupling) out->near[f].push_back({g, std::exp(-out->d[f * n + g] / P.coupling)});
  std::lock_guard<std::mutex> lock(mutex_);
  return cache_.emplace(key, std::move(out)).first->second;
}

// ---- the value ------------------------------------------------------------------------------------------------------

// The week's lines of a series of 168 hours: |X̂(m)|² folded, m = 0…84, X̂(m) = (1/168) Σ_h x_h e^{−2πimh/168}; by
// Parseval Σ_m (folded) = mean of x².
inline std::vector<double> week_power(std::span<const double> x) {
  static const auto twiddle = [] {  // cos and sin of 2π j/168, j = 0…167
    std::array<std::array<double, 2>, week> t{};
    for (std::size_t j = 0; j < week; j++) t[j] = {std::cos(2 * std::numbers::pi * static_cast<double>(j) / week), std::sin(2 * std::numbers::pi * static_cast<double>(j) / week)};
    return t;
  }();
  std::vector<double> p(week / 2 + 1, 0.0);
  for (std::size_t m = 0; m <= week / 2; m++) {
    double re = 0, im = 0;
    for (std::size_t h = 0; h < week; h++) {
      const auto& w = twiddle[(m * h) % week];
      re += x[h] * w[0], im -= x[h] * w[1];
    }
    re /= static_cast<double>(week), im /= static_cast<double>(week);
    p[m] = (m == 0 || m == week / 2 ? 1.0 : 2.0) * (re * re + im * im);
  }
  return p;
}

inline Outcome Problem::evaluate(const Plan& plan, std::span<const double> passers) const {
  validate(plan);
  const std::size_t n = size();
  if (passers.size() != n * week) throw std::invalid_argument("economy::Problem::evaluate: passers for each frontage and hour of the week");
  const auto Dist = distances(plan.links);
  Outcome o;
  o.spending_by_frontage.assign(n, 0.0), o.value_by_frontage.assign(n, 0.0), o.visits_by_frontage.assign(n, 0.0);
  std::vector<double> p(n, 0.0);
  std::vector<char> open(n, 0);
  for (std::size_t h = 0; h < week; h++) {
    for (std::size_t f = 0; f < n; f++) open[f] = plan.open(f, h);
    // The chance a passer stops, with the catalysis of the open neighbours.
    for (std::size_t f = 0; f < n; f++) {
      p[f] = 0;
      if (!open[f] || !trades(plan.kind[f])) continue;
      const auto k = static_cast<std::size_t>(plan.kind[f]);
      double kappa = 0;
      for (const auto& [g, w] : Dist->near[f])
        if (open[g]) kappa += w * P.affinity[static_cast<std::size_t>(plan.kind[g])][k];
      p[f] = P.propensity[k] * P.period[k][band(h)] * std::clamp(1 + P.catalysis * kappa, P.catalysis_floor, P.catalysis_ceiling);
    }
    // Competition on each segment, kind by kind.
    for (std::size_t f = 0; f < n; f++) {
      if (p[f] <= 0) continue;
      double none = 1 - p[f], sum = p[f];
      for (auto g : on_segment[f])
        if (p[g] > 0 && plan.kind[g] == plan.kind[f]) none *= 1 - p[g], sum += p[g];
      const double share = (1 - none) / sum;
      const auto k = static_cast<std::size_t>(plan.kind[f]);
      const double v = passers[f * week + h] * p[f] * share;
      const double s = v * P.spend[k];
      o.visits += v, o.spending += s, o.margin += s * P.margin[k];
      o.spending_by_kind[k] += s;
      o.spending_by_hour[h] += s, o.visits_by_hour[h] += v;
      o.spending_by_frontage[f] += s, o.visits_by_frontage[f] += v;
      o.value_by_frontage[f] += s * P.margin[k];
    }
    for (std::size_t f = 0; f < n; f++) o.footfall_by_hour[h] += passers[f * week + h];
  }
  // Costs.
  for (std::size_t f = 0; f < n; f++) {
    const auto k = static_cast<std::size_t>(plan.kind[f]);
    const auto& x = frontages[f];
    double c = 0;
    if (plan.kind[f] != vacant) {
      const double staffing = P.staff[k] * P.wage * (plan.weekly_hours(f) + plan.days_open(f));
      const double out = P.outgoings[k] * x.trading_area;
      o.staffing += staffing, o.outgoings += out, c += staffing + out;
      o.frontages_by_kind[k] += 1;
      if (plan.weekly_hours(f) > 0) o.open_frontages++;
    } else {
      o.frontages_by_kind[vacant] += 1;
    }
    if (plan.kind[f] != x.kind0) {
      o.changed++;
      if (plan.kind[f] != vacant) {
        const double r = P.refit_weekly(x.trading_area);
        o.refits += r, c += r;
      }
    }
    o.value_by_frontage[f] -= c;
  }
  for (std::size_t k = 0; k < links.size(); k++)
    if (plan.links[k]) o.link_costs += P.link_weekly(links[k].length), o.links++;
  o.value = o.margin - o.staffing - o.outgoings - o.refits - o.link_costs;
  for (double v : o.footfall_by_hour) o.footfall += v;

  // Spread: the week's lines of the spending.
  o.power = week_power(o.spending_by_hour);
  double total = 0, daily = 0;
  for (std::size_t m = 0; m < o.power.size(); m++) {
    total += o.power[m];
    if (m > 0 && m % 7 == 0) daily += o.power[m];
  }
  o.spread = total > 0 ? o.power[0] / total : 0;
  o.effective_hours = week * o.spread;
  o.daily_share = total > 0 ? daily / total : 0;
  o.weekly_share = total > 0 ? 1 - o.spread - o.daily_share : 0;
  {
    double s = 0, s2 = 0;
    for (double v : o.footfall_by_hour) s += v, s2 += v * v;
    o.footfall_spread = s2 > 0 ? (s / week) * (s / week) / (s2 / week) : 0;
  }
  double evening = 0, weekend = 0;
  for (std::size_t h = 0; h < week; h++) {
    if (band(h) == 1) evening += o.spending_by_hour[h];
    if (is_weekend(h)) weekend += o.spending_by_hour[h];
  }
  o.evening_share = o.spending > 0 ? evening / o.spending : 0;
  o.weekend_share = o.spending > 0 ? weekend / o.spending : 0;

  // Mix: H(K | block), the kinds of the frontages (vacant left out) weighted by their length.
  {
    std::vector<std::vector<double>> rows(place_count, std::vector<double>(kinds - 1, 0.0));
    double any = 0;
    for (std::size_t f = 0; f < n; f++)
      if (plan.kind[f] != vacant) {
        const double w = std::max(1.0, frontages[f].length);
        rows[static_cast<std::size_t>(frontages[f].block)][static_cast<std::size_t>(plan.kind[f] - 1)] += w, any += w;
      }
    if (any > 0) {
      std::vector<std::vector<double>> kept;
      for (auto& r : rows)
        if (std::accumulate(r.begin(), r.end(), 0.0) > 0) kept.push_back(r);
      const auto mix = physis::statistical_mechanics::mixing(calx::linalg::Rectangular::from_rows(kept));
      o.mix = mix.conditional, o.mix_entropy = mix.entropy;
    }
  }

  // Access: to the nearest open frontage of each daily need, from where the people are, over the week.
  {
    std::array<std::array<std::vector<double>, 48>, 3> nearest;  // need × (weekday hour, weekend hour) ↦ per frontage
    for (std::size_t q = 0; q < 3; q++)
      for (std::size_t state = 0; state < 48; state++) {
        const std::size_t h = state < 24 ? state : 120 + (state - 24);  // a Monday hour, or a Saturday hour
        auto& d = nearest[q][state];
        d.assign(n, P.access_cap);
        for (std::size_t g = 0; g < n; g++) {
          if (plan.kind[g] != needs[q] || !plan.open(g, h)) continue;
          for (std::size_t f = 0; f < n; f++) d[f] = std::min(d[f], Dist->at(f, g));
        }
      }
    double w = 0;
    std::array<double, 3> a{};
    for (std::size_t h = 0; h < week; h++) {
      const std::size_t state = is_weekend(h) ? 24 + hour_of_day(h) : hour_of_day(h);
      for (std::size_t f = 0; f < n; f++) {
        const double people = passers[f * week + h];
        w += people;
        for (std::size_t q = 0; q < 3; q++) a[q] += people * nearest[q][state][f];
      }
    }
    for (std::size_t q = 0; q < 3; q++) o.access_by_need[q] = w > 0 ? a[q] / w : P.access_cap;
    o.access = (o.access_by_need[0] + o.access_by_need[1] + o.access_by_need[2]) / 3;
  }
  return o;
}

// ---- the canonical relaxation ------------------------------------------------------------------------------------
//
// u ∈ [0, 1]^(F·10 + L): for each frontage, six scores of the kinds (the kind the argmax, the first on a tie) and four
// numbers for its hours (start ⌊24u⌋ ∧ 23, length ⌊25u⌋ ∧ 24, weekday then weekend); then a number for each link (made
// when ≥ ½). encode() gives the centre of the plan's cell, so that decode(encode(p)) = p.

inline Plan Problem::decode(std::span<const double> u) const {
  if (u.size() != dimension()) throw std::invalid_argument("economy::Problem::decode: a point of another dimension");
  Plan p;
  const std::size_t n = size(), w = kinds + hour_genes;
  for (std::size_t f = 0; f < n; f++) {
    const double* x = u.data() + f * w;
    std::size_t k = 0;
    for (std::size_t j = 1; j < kinds; j++)
      if (x[j] > x[k]) k = j;
    p.kind.push_back(static_cast<std::int64_t>(k));
    for (std::size_t j = 0; j < hour_genes; j++) {
      const double v = std::clamp(x[kinds + j], 0.0, 1.0);
      p.hours.push_back(j % 2 == 0 ? std::min<std::int64_t>(23, static_cast<std::int64_t>(std::floor(24 * v))) : std::min<std::int64_t>(24, static_cast<std::int64_t>(std::floor(25 * v))));
    }
  }
  for (std::size_t k = 0; k < links.size(); k++) p.links.push_back(u[n * w + k] >= 0.5 ? 1 : 0);
  return p;
}

inline std::vector<double> Problem::encode(const Plan& p) const {
  validate(p);
  const std::size_t n = size(), w = kinds + hour_genes;
  std::vector<double> u(dimension(), 0.0);
  for (std::size_t f = 0; f < n; f++) {
    double* x = u.data() + f * w;
    for (std::size_t j = 0; j < kinds; j++) x[j] = static_cast<std::int64_t>(j) == p.kind[f] ? 0.75 : 0.25;
    for (std::size_t j = 0; j < hour_genes; j++) {
      const double v = static_cast<double>(p.hours[hour_genes * f + j]);
      x[kinds + j] = j % 2 == 0 ? (v + 0.5) / 24 : (v + 0.5) / 25;
    }
  }
  for (std::size_t k = 0; k < links.size(); k++) u[n * w + k] = p.links[k] ? 0.75 : 0.25;
  return u;
}

// ---- the mixed encoding's kernels (README §3.4), for the studies that search it --------------------------------------

struct Kernels {
  const Problem* problem = nullptr;
  double crossover_probability = 0.9;
  std::int64_t creep = 3;  // hours

  Plan sample(calx::random::Stream& s) const {
    Plan p;
    p.kind = hv::uniform_integers(problem->kind_box(), s);
    p.hours = hv::uniform_integers(problem->hour_box(), s);
    p.links = hv::random_bits(problem->links.size(), s);
    return p;
  }
  // A kind redrawn at one frontage in F; an hour gene crept by at most `creep` at one in F; a link flipped at one in L.
  Plan mutate(const Plan& x, calx::random::Stream& s) const {
    const double n = static_cast<double>(problem->size());
    Plan y;
    y.kind = hv::reset_mutation(x.kind, problem->kind_box(), 1.0 / n, s);
    y.hours = hv::creep_mutation(x.hours, problem->hour_box(), creep, 1.0 / n, s);
    y.links = x.links;
    if (!y.links.empty()) hv::standard_bit_mutation(y.links, 1.0 / static_cast<double>(y.links.size()), s);
    return y;
  }
  // A frontage's kind and hours travel together; frontages are in order along the area, so two cuts exchange a stretch
  // of it (neighbours catalyse and compete with one another). Links cross uniformly.
  std::pair<Plan, Plan> crossover(const Plan& a, const Plan& b, calx::random::Stream& s) const {
    if (s.uniform() >= crossover_probability) return {a, b};
    const std::size_t n = problem->size();
    std::size_t i = s.below(n + 1), j = s.below(n + 1);
    if (i > j) std::swap(i, j);
    Plan x = a, y = b;
    for (std::size_t f = i; f < j; f++) {
      std::swap(x.kind[f], y.kind[f]);
      for (std::size_t q = 0; q < hour_genes; q++) std::swap(x.hours[hour_genes * f + q], y.hours[hour_genes * f + q]);
    }
    for (std::size_t k = 0; k < x.links.size(); k++)
      if (s.uniform() < 0.5) std::swap(x.links[k], y.links[k]);
    return {x, y};
  }
};

// ---- output ----------------------------------------------------------------------------------------------------------

inline Json numbers(std::span<const double> v, int digits = 6) {
  Json a = Json::array();
  const double p = std::pow(10.0, digits);
  for (double x : v) a.push(std::isfinite(x) ? Json(std::round(x * p) / p) : Json());
  return a;
}

inline Json plan_json(const Plan& p) {
  Json k = Json::array(), h = Json::array(), l = Json::array();
  for (auto x : p.kind) k.push(Json(static_cast<double>(x)));
  for (auto x : p.hours) h.push(Json(static_cast<double>(x)));
  for (auto x : p.links) l.push(Json(static_cast<double>(x)));
  return Json::object({{"kind", k}, {"hours", h}, {"links", l}});
}

inline Plan plan_of(const Json& j) {
  Plan p;
  for (const auto& x : j["kind"].as_array()) p.kind.push_back(static_cast<std::int64_t>(x.as_number()));
  for (const auto& x : j["hours"].as_array()) p.hours.push_back(static_cast<std::int64_t>(x.as_number()));
  for (const auto& x : j["links"].as_array()) p.links.push_back(static_cast<std::uint8_t>(x.as_number()));
  return p;
}

// The scalars of an outcome (and, if asked, its profiles and per-frontage values).
inline Json outcome_json(const Outcome& o, bool profiles = false) {
  const auto r = [](double x, double d) { return std::round(x * d) / d; };
  Json j = Json::object({{"value_aud_week", r(o.value, 1)}, {"spending_aud_week", r(o.spending, 1)}, {"margin_aud_week", r(o.margin, 1)}, {"visits_week", r(o.visits, 1)},
                         {"staffing_aud_week", r(o.staffing, 1)}, {"outgoings_aud_week", r(o.outgoings, 1)}, {"refits_aud_week", r(o.refits, 1)}, {"links_aud_week", r(o.link_costs, 1)},
                         {"footfall_week", r(o.footfall, 1)}, {"spread", r(o.spread, 1e6)}, {"effective_hours", r(o.effective_hours, 1e3)}, {"daily_share", r(o.daily_share, 1e6)},
                         {"weekly_share", r(o.weekly_share, 1e6)}, {"footfall_spread", r(o.footfall_spread, 1e6)}, {"evening_share", r(o.evening_share, 1e6)},
                         {"weekend_share", r(o.weekend_share, 1e6)}, {"mix_nats", r(o.mix, 1e6)}, {"kind_entropy_nats", r(o.mix_entropy, 1e6)}, {"access_m", r(o.access, 100)},
                         {"access_by_need_m", numbers(o.access_by_need, 2)}, {"open_frontages", static_cast<double>(o.open_frontages)}, {"links", static_cast<double>(o.links)},
                         {"changed", static_cast<double>(o.changed)}, {"spending_by_kind_aud_week", numbers(o.spending_by_kind, 1)}, {"frontages_by_kind", numbers(o.frontages_by_kind, 0)}});
  if (profiles) {
    j.set("spending_by_hour", numbers(o.spending_by_hour, 2));
    j.set("visits_by_hour", numbers(o.visits_by_hour, 3));
    j.set("footfall_by_hour", numbers(o.footfall_by_hour, 1));
    j.set("spending_by_frontage", numbers(o.spending_by_frontage, 1));
    j.set("value_by_frontage", numbers(o.value_by_frontage, 1));
    j.set("visits_by_frontage", numbers(o.visits_by_frontage, 2));
    j.set("power", numbers(o.power, 4));
  }
  return j;
}

}  // namespace economy
