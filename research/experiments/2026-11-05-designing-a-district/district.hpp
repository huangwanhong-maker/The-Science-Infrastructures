// The redesign of a few blocks of Melbourne's Hoddle Grid as one problem with several kinds of decision and measures
// from several domains (README §2).
//
// The decision. For each parcel i (a building of the city model standing in the blocks between Swanston and Russell
// Streets, Flinders Lane and Bourke Street, its footprint kept as built): its storeys sᵢ ∈ {1, …, bᵢ} and its use
// uᵢ ∈ {dwellings, offices, shops, open space}; for each candidate link k of the walking network, whether it is made,
// pₖ ∈ {0, 1}. A link is a blind lane carried through to the street beyond it (the study's own edit of the network,
// lanes.hpp), or a passage through a block that is a 2-cell (the domain's open_passage). A parcel given to open space
// has no building: its ground becomes public.
//
// The measures, each from the model that owns it:
//   floor area   Σ over built parcels of the footprint's area (architecture::space::massing, after its pooling of
//                corners) times the storeys, less the ground floor a link takes (4 m by its length);
//   winter sun   the sunlit area-hours of the public ground (the streets, the existing open spaces, the parcels given to
//                open space) between 10:00 and 14:00 on 21 June 2026, the Sun from the astronomy domain, the shadows by
//                the study's own model (sunlight.hpp);
//   heating      the heat the buildings need, held at 20 °C through 21 June under the air of an average June day, per
//                square metre of floor: architecture::environment's massing heat (conductances through facades, roofs
//                and the ground, party walls shared) and its sol-air gains of the Sun (massing_gains);
//   walking      the mean walking distance from a home to its nearest shop and to its nearest open space, homes
//                weighted by their floor area: the urban-science domain's district (its walking network as a complex,
//                built again from the network with the lanes carried through, passages opened by open_passage,
//                holdings placed by attach_store), distances by the study's own Dijkstra (walking.hpp);
//   mix          the conditional entropy H(U | X) of the floor's use given its block, the mean local mix of the
//                urban-science domain's Definition 2.6 (theory/street-oceanography.md §2.3), on physis's mixing.
//
// The search sees them as five minimised objectives, each relative to the district as built (so that the district as
// built is the point (−1, −1, 1, 1, −H₀/ln 3)), with four constraints (assumptions of the study): open space at least
// 15% of the land, homes at least 20% of the floor, shops at least 5% of it and at most 20%.

#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <numbers>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <architecture/environment/massing_heat.hpp>
#include <architecture/space/massing.hpp>
#include <calx/linalg/rectangular.hpp>
#include <heuresis/search/evaluation.hpp>
#include <heuresis/variation/binary.hpp>
#include <heuresis/variation/integer.hpp>
#include <physis/computational_geometry/planar.hpp>
#include <physis/statistical_mechanics/mixing.hpp>
#include <science/json.hpp>
#include <urban/landscape/district.hpp>

#include "lanes.hpp"
#include "sunlight.hpp"
#include "walking.hpp"

namespace study {

namespace ul = urban::landscape;
namespace as = architecture::space;
namespace ae = architecture::environment;
namespace hs = heuresis::search;
namespace hv = heuresis::variation;
namespace cg = physis::computational_geometry;
using science::Json;

enum Use : std::int64_t { dwellings = 0, offices = 1, shops = 2, open_space = 3 };
inline const char* use_name(std::int64_t u) {
  static const char* names[] = {"dwellings", "offices", "shops", "open space"};
  return names[std::clamp<std::int64_t>(u, 0, 3)];
}

// The assumptions of the problem (README §2.4).
struct Assumptions {
  double storey = 3.2;            // m: the city model's storey
  double slenderness = 8;         // a building's height at most this many times the square root of its footprint
  double setpoint = 20;           // °C: the buildings held at this
  double passage_width = 4;       // m: the ground floor an arcade through a block takes
  double min_open_share = 0.15;   // of the land (parcels and existing open spaces)
  double min_dwelling_share = 0.20, min_shop_share = 0.05, max_shop_share = 0.20;  // of the floor
  double sample = 4;              // m: the lattice of ground samples
  double margin = 16;             // m: the ground sampled beyond the area's street centrelines (the far half of its streets)
  double min_block = 1000;        // m²: a block large enough for a passage
  double min_link = 15;           // m: a passage, or a blind lane carried through, at least this long
  double reach = 120;             // m: a blind lane carried through at most this far
};

struct Design {
  hv::Integers storeys;  // per parcel
  hv::Integers uses;     // per parcel
  hv::Bits links;        // per candidate link: the lanes carried through, then the passages
};

struct Parcel {
  std::string id, name, clue_property, clue_use, height_source;
  Ring ring;
  double area = 0;
  std::int64_t floors0 = 1, use0 = offices, bound = 1;
  Point centre{};
  double along = 0;   // its position along the area's long axis (the genotype's order)
  int place = -1;     // the block it stands in
};

struct Link {
  std::string kind;      // "lane carried through" or "passage"
  Point a{}, b{};
  double length = 0;
  std::string from, to;  // the lane or street at each end
  std::size_t through = 0;  // the lane's index among the candidates (kind "lane carried through")
};

struct Measures {
  double floor_area = 0, sun = 0, heat = 0, access = 0, mix = 0;  // m², m²·h, kWh/(m² day), m, nats
  double access_shop = 0, access_open = 0, mix_entropy = 0, mix_departure = 0;
  double open_share = 0, dwelling_share = 0, shop_share = 0, plot_ratio = 0;
  std::array<double, 3> floor_by_use{};
  double heat_total = 0;  // kWh over the day
  std::size_t built = 0, open = 0, links = 0;
  bool ok = true;
  std::string reason;
};

// The walking network with a set of links made, with every parcel and open space placed on it.
struct Walk {
  WalkGraph graph;
  std::vector<OnSegment> parcel, space;
};

class Problem {
 public:
  Assumptions A;
  std::vector<Parcel> parcels;
  std::vector<Link> links;                        // the lanes carried through, then the passages
  std::vector<Through> throughs;
  ul::Network plane;                              // the network as the domain planarised it
  std::vector<std::vector<Point>> open_spaces;  // the existing ones inside the area (fixed)
  std::vector<Point> open_space_centre;
  double open_space_area = 0, land = 0;
  ul::District base;                              // the district as the domain builds it, no passage opened
  std::vector<std::string> dropped;               // segments the domain's planarisation dropped
  std::vector<std::uint64_t> places;              // the blocks (2-cells) and then the holes, as places
  std::size_t place_count = 0;
  Shading shade;
  std::vector<Sun> day;                           // the Sun at each hour of the day (heat)
  std::vector<double> air;                        // °C at each hour
  double ground = 15;
  Measures reference;                             // the district as built
  Ring quad;
  // The time spent in each measure (massing, heating, sun, walking, mix), in nanoseconds, summed over the evaluations.
  mutable std::array<std::atomic<long long>, 5> stage_ns{};

  explicit Problem(const Json& in);

  Design existing() const {
    Design d;
    for (const auto& p : parcels) d.storeys.push_back(p.floors0), d.uses.push_back(p.use0);
    d.links.assign(links.size(), 0);
    return d;
  }
  hv::IntegerBox storey_box() const {
    hv::IntegerBox b;
    for (const auto& p : parcels) b.lower.push_back(1), b.upper.push_back(p.bound);
    return b;
  }
  hv::IntegerBox use_box() const { return hv::IntegerBox::cube(parcels.size(), 0, 3); }

  Measures measure(const Design& d) const;
  std::vector<double> heights(const Design& d) const {
    std::vector<double> h(parcels.size(), 0.0);
    for (std::size_t i = 0; i < parcels.size(); i++) h[i] = d.uses[i] == open_space ? 0.0 : A.storey * static_cast<double>(d.storeys[i]);
    return h;
  }
  std::vector<char> opened(const Design& d) const {
    std::vector<char> o(parcels.size());
    for (std::size_t i = 0; i < parcels.size(); i++) o[i] = d.uses[i] == open_space;
    return o;
  }
  // The five objectives relative to the district as built, all minimised.
  std::vector<double> objectives(const Measures& m) const {
    return {-m.floor_area / reference.floor_area, -m.sun / reference.sun, m.heat / reference.heat, m.access / reference.access, -m.mix / std::log(3.0)};
  }
  std::vector<double> constraints(const Measures& m) const {
    return {A.min_open_share - m.open_share, A.min_dwelling_share - m.dwelling_share, A.min_shop_share - m.shop_share, m.shop_share - A.max_shop_share};
  }
  hs::Evaluation evaluate(const Design& d) const {
    const Measures m = measure(d);
    if (!m.ok) return hs::failure(m.reason);
    return hs::evaluation(objectives(m), constraints(m), {m.plot_ratio, m.open_share});
  }
  // Two objectives (floor area and winter sun) under the same constraints.
  hs::Evaluation evaluate_two(const Design& d) const {
    const Measures m = measure(d);
    if (!m.ok) return hs::failure(m.reason);
    const auto f = objectives(m);
    return hs::evaluation({f[0], f[1]}, constraints(m), {m.plot_ratio, m.open_share});
  }
  // For the repertoire: winter sun as the quality, the open share a descriptor rather than a constraint.
  hs::Evaluation evaluate_repertoire(const Design& d) const {
    const Measures m = measure(d);
    if (!m.ok) return hs::failure(m.reason);
    const auto f = objectives(m);
    const auto g = constraints(m);
    return hs::evaluation({f[1]}, {g[1], g[2], g[3]}, {m.plot_ratio, m.open_share});
  }
  std::vector<double> sun_map(const Design& d) const { return sun_hours(shade, heights(d), opened(d)); }

  std::shared_ptr<const Walk> walk(const hv::Bits& open) const;
  // The district with a set of links made: the lanes carried through on the network, the district built from it by
  // the domain, the passages opened in it by the domain.
  ul::District with_links(const hv::Bits& open) const;

 private:
  mutable std::mutex mutex_;
  mutable std::map<std::uint64_t, std::shared_ptr<const Walk>> walks_;
};

// ---- construction ---------------------------------------------------------------------------------------------

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

inline Problem::Problem(const Json& in) {
  quad = ring_of(in["redesign_quad"]);
  A.storey = in["storey_m"].as_number();
  const auto smax = static_cast<std::int64_t>(in["storeys_bound"].as_number());
  // The area's long axis: from the middle of its southern side to the middle of its northern (along Swanston Street).
  const Point s0{(quad[0][0] + quad[1][0]) / 2, (quad[0][1] + quad[1][1]) / 2}, s1{(quad[2][0] + quad[3][0]) / 2, (quad[2][1] + quad[3][1]) / 2};
  const double L = std::hypot(s1[0] - s0[0], s1[1] - s0[1]);
  const Point axis{(s1[0] - s0[0]) / L, (s1[1] - s0[1]) / L};
  for (const auto& j : in["parcels"].as_array()) {
    Parcel p;
    p.id = j["id"].text(), p.name = j["name"].text(), p.clue_property = j["clue_property"].text(), p.clue_use = j["clue_use"].text();
    p.height_source = j["height_source"].text();
    p.ring = ring_of(j["footprint"]);
    p.area = cg::ring_area(p.ring);
    p.floors0 = static_cast<std::int64_t>(j["floors"].as_number());
    p.use0 = static_cast<std::int64_t>(j["use"].as_number());
    p.bound = std::max<std::int64_t>(1, std::min<std::int64_t>(smax, static_cast<std::int64_t>(std::floor(A.slenderness * std::sqrt(p.area) / A.storey))));
    p.centre = centroid(p.ring);
    p.along = (p.centre[0] - s0[0]) * axis[0] + (p.centre[1] - s0[1]) * axis[1];
    parcels.push_back(std::move(p));
  }
  // The genotype's order: along the area's long axis, so that a two-point crossover exchanges a stretch of the area.
  std::stable_sort(parcels.begin(), parcels.end(), [](const Parcel& a, const Parcel& b) { return a.along < b.along; });
  for (const auto& s : in["open_spaces"].as_array()) {
    open_spaces.push_back(ring_of(s["ring"]));
    open_space_centre.push_back(centroid(open_spaces.back()));
    open_space_area += cg::ring_area(open_spaces.back());
  }
  land = open_space_area;
  for (const auto& p : parcels) land += p.area;

  // The district as the urban-science domain builds it from the workbench's network (planarised, as its engine does).
  ul::Network N;
  for (const auto& p : in["network"]["junctions"].as_array()) N.junctions.push_back({p[0].as_number(), p[1].as_number()});
  for (const auto& s : in["network"]["segments"].as_array()) {
    ul::Segment g;
    g.u = static_cast<std::size_t>(s["u"].as_number()), g.v = static_cast<std::size_t>(s["v"].as_number());
    for (const auto& p : s["line"].as_array()) g.line.push_back({p[0].as_number(), p[1].as_number()});
    g.kind = s["kind"].text("street"), g.name = s["name"].text();
    N.segments.push_back(std::move(g));
  }
  plane = ul::planarised(std::move(N), &dropped);
  base = ul::district(plane);

  // Places for the mix: the blocks, then the faces left as holes; each parcel in the one holding its centre.
  for (const auto id : base.of(2)) places.push_back(id);
  place_count = places.size() + base.hole_rings.size();
  for (auto& p : parcels) {
    for (std::size_t k = 0; k < places.size() && p.place < 0; k++)
      if (cg::inside_ring(base.ring(places[k]), p.centre)) p.place = static_cast<int>(k);
    for (std::size_t k = 0; k < base.hole_rings.size() && p.place < 0; k++)
      if (cg::inside_ring(base.hole_rings[k], p.centre)) p.place = static_cast<int>(places.size() + k);
    if (p.place < 0) p.place = static_cast<int>(place_count++);
  }

  // Candidate links: every blind lane of the area carried through (lanes.hpp); then one passage through each block
  // (2-cell) inside the area, across its middle between its two long sides (the domain's plan_passage).
  throughs = lanes_through(plane, quad, A.min_link, A.reach);
  for (std::size_t k = 0; k < throughs.size(); k++) links.push_back({"lane carried through", throughs[k].from, throughs[k].to, throughs[k].length, throughs[k].lane, throughs[k].street, k});
  for (const auto id : base.of(2)) {
    const auto r = base.ring(id);
    const Point c = centroid(r);
    if (!cg::inside_ring(quad, c) || !(base.area(id) >= A.min_block)) continue;
    // The ring's principal axes, from the second moments of points walked along it every metre.
    double sxx = 0, syy = 0, sxy = 0, n = 0;
    for (std::size_t k = 0; k < r.size(); k++) {
      const auto& a = r[k];
      const auto& b = r[(k + 1) % r.size()];
      const double l = std::hypot(b[0] - a[0], b[1] - a[1]);
      for (double t = 0; t < l; t += 1) {
        const double x = a[0] + (b[0] - a[0]) * t / l - c[0], y = a[1] + (b[1] - a[1]) * t / l - c[1];
        sxx += x * x, syy += y * y, sxy += x * y, n += 1;
      }
    }
    const double theta = 0.5 * std::atan2(2 * sxy, sxx - syy);  // the long axis
    const Point across{-std::sin(theta), std::cos(theta)};
    const Point a{c[0] + 500 * across[0], c[1] + 500 * across[1]}, b{c[0] - 500 * across[0], c[1] - 500 * across[1]};
    const auto plan = ul::plan_passage(base, id, a, b);
    if (!plan.refused.empty()) continue;
    const double len = std::hypot(plan.a[0] - plan.b[0], plan.a[1] - plan.b[1]);
    if (len < A.min_link) continue;
    links.push_back({"passage", plan.a, plan.b, len, base.name.at(plan.edge_a), base.name.at(plan.edge_b), 0});
  }

  // The Sun and the air.
  const double deg = std::numbers::pi / 180;
  std::vector<Sun> window;
  for (const auto& s : in["sun"]["window"].as_array()) window.push_back({s["altitude_deg"].as_number() * deg, s["azimuth_deg"].as_number() * deg, 0.5});
  for (const auto& s : in["sun"]["hours"].as_array()) day.push_back({s["altitude_deg"].as_number() * deg, s["azimuth_deg"].as_number() * deg, 1.0});
  for (const auto& t : in["air_by_hour_c"].as_array()) air.push_back(t.as_number());
  ground = in["ground_c"].as_number();

  // The ground sampled on a lattice over the area and the far half of its streets; what stands on it.
  Ring wide;
  {
    const Point c = centroid(quad);
    for (const auto& q : quad) {
      const double dx = q[0] - c[0], dy = q[1] - c[1], r = std::hypot(dx, dy);
      wide.push_back({q[0] + A.margin * std::sqrt(2.0) * dx / r, q[1] + A.margin * std::sqrt(2.0) * dy / r});
    }
  }
  std::vector<Caster> casters;
  std::vector<Ring> context_rings;
  for (const auto& j : in["context"].as_array()) {
    Caster c;
    c.ring = ring_of(j["footprint"]);
    c.height = j["height_m"].as_number();
    context_rings.push_back(c.ring);
    casters.push_back(std::move(c));
  }
  for (std::size_t i = 0; i < parcels.size(); i++) {
    Caster c;
    c.ring = parcels[i].ring, c.parcel = static_cast<int>(i);
    casters.push_back(std::move(c));
  }
  double x0 = 1e300, y0 = 1e300, x1 = -1e300, y1 = -1e300;
  for (const auto& q : wide) x0 = std::min(x0, q[0]), y0 = std::min(y0, q[1]), x1 = std::max(x1, q[0]), y1 = std::max(y1, q[1]);
  std::vector<Point> pts;
  std::vector<int> owner;
  for (double y = std::floor(y0 / A.sample) * A.sample + A.sample / 2; y < y1; y += A.sample)
    for (double x = std::floor(x0 / A.sample) * A.sample + A.sample / 2; x < x1; x += A.sample) {
      const Point p{x, y};
      if (!inside(wide, p)) continue;
      bool built = false;
      for (const auto& r : context_rings)
        if (inside(r, p)) { built = true; break; }
      if (built) continue;
      int own = -1;
      for (std::size_t i = 0; i < parcels.size() && own < 0; i++)
        if (inside(parcels[i].ring, p)) own = static_cast<int>(i);
      pts.push_back(p), owner.push_back(own);
    }
  shade = shading(std::move(pts), std::move(owner), std::move(window), A.sample * A.sample, std::move(casters), 900.0);

  reference = measure(existing());
  if (!reference.ok) throw std::runtime_error("the district as built cannot be measured: " + reference.reason);
}

// ---- the walking network with passages ---------------------------------------------------------------------------

inline ul::District Problem::with_links(const hv::Bits& open) const {
  std::vector<Through> chosen;
  for (std::size_t k = 0; k < links.size(); k++)
    if (open[k] && links[k].kind == "lane carried through") chosen.push_back(throughs[links[k].through]);
  ul::District D = chosen.empty() ? base : ul::district(carry_through(plane, chosen, A.reach));
  for (std::size_t k = 0; k < links.size(); k++) {
    if (!open[k] || links[k].kind != "passage") continue;
    // The block (in this district) that holds the passage's middle.
    const Point mid{(links[k].a[0] + links[k].b[0]) / 2, (links[k].a[1] + links[k].b[1]) / 2};
    std::uint64_t block = 0;
    for (const auto id : D.of(2))
      if (cg::inside_ring(D.ring(id), mid)) { block = id; break; }
    if (!block) throw std::invalid_argument("study::Problem::with_links: no block holds the passage between " + links[k].from + " and " + links[k].to);
    ul::open_passage(D, block, links[k].a, links[k].b);
  }
  return D;
}

inline std::shared_ptr<const Walk> Problem::walk(const hv::Bits& open) const {
  std::uint64_t key = 0;
  for (std::size_t k = 0; k < open.size(); k++)
    if (open[k]) key |= std::uint64_t{1} << k;
  {
    std::lock_guard<std::mutex> lock(mutex_);
    const auto it = walks_.find(key);
    if (it != walks_.end()) return it->second;
  }
  ul::District D = with_links(open);
  auto W = std::make_shared<Walk>();
  std::map<std::size_t, std::size_t> node;
  std::map<std::uint64_t, std::size_t> edge;
  for (std::size_t i = 0; i < D.complex.size(); i++)
    if (D.complex.cells[i].dim == 0) node[i] = W->graph.nodes++;
  for (std::size_t i = 0; i < D.complex.size(); i++) {
    const auto& c = D.complex.cells[i];
    if (c.dim != 1) continue;
    const auto ends = D.complex.ends(i);
    edge[c.id] = W->graph.edges.size();
    W->graph.edges.push_back({node.at(ends[0]), node.at(ends[1]), D.length(c.id)});
  }
  W->graph.finish();
  for (const auto& p : parcels) {
    const auto k = ul::attach_store(D, {p.id, "", p.centre});
    W->parcel.push_back({edge.at(D.stores[k].segment), D.stores[k].s});
  }
  for (const auto& c : open_space_centre) {
    const auto k = ul::attach_store(D, {"open space", "", c});
    W->space.push_back({edge.at(D.stores[k].segment), D.stores[k].s});
  }
  std::lock_guard<std::mutex> lock(mutex_);
  return walks_.emplace(key, std::move(W)).first->second;
}

// ---- the measures -------------------------------------------------------------------------------------------------

inline Measures Problem::measure(const Design& d) const {
  Measures m;
  const std::size_t n = parcels.size();
  auto clock = std::chrono::steady_clock::now();
  const auto tick = [&](int k) {
    const auto now = std::chrono::steady_clock::now();
    stage_ns[static_cast<std::size_t>(k)] += std::chrono::duration_cast<std::chrono::nanoseconds>(now - clock).count();
    clock = now;
  };
  if (d.storeys.size() != n || d.uses.size() != n || d.links.size() != links.size()) throw std::invalid_argument("study::Problem::measure: a design of another district");
  // The buildings as one complex (architecture's massing), and their floor.
  std::vector<as::Footprint> F;
  std::vector<std::size_t> of_building;
  for (std::size_t i = 0; i < n; i++) {
    if (d.uses[i] == open_space) { m.open++; continue; }
    F.push_back({parcels[i].id, parcels[i].name, parcels[i].ring, 0.0, A.storey * static_cast<double>(d.storeys[i])});
    of_building.push_back(i);
  }
  m.built = F.size();
  if (F.empty()) { m.ok = false, m.reason = "no building stands"; return m; }
  const auto M = as::massing(F);
  tick(0);
  if (!M.skipped.empty()) { m.ok = false, m.reason = "the massing skipped " + F[M.skipped[0].first].key + ": " + M.skipped[0].second; return m; }
  std::vector<double> floor(n, 0.0);
  double built_floor = 0;
  for (std::size_t b = 0; b < M.building.size(); b++) {
    const std::size_t i = of_building[M.building[b]];
    floor[i] = M.footprint_area[b] * static_cast<double>(d.storeys[i]);
    m.floor_by_use[static_cast<std::size_t>(d.uses[i])] += floor[i];
    built_floor += floor[i];
  }
  double taken = 0;
  for (std::size_t k = 0; k < links.size(); k++)
    if (d.links[k]) taken += A.passage_width * links[k].length, m.links++;
  m.floor_area = built_floor - taken;
  // Heating: each building held at the setpoint; with every building at one temperature the party walls carry nothing,
  // and what each needs is its loss to the air and the ground less the Sun's gains, never below zero.
  const auto H = ae::massing_heat(M);
  double joules = 0;
  for (std::size_t h = 0; h < day.size(); h++) {
    const auto S = ae::massing_gains(M, day[h].altitude, day[h].azimuth);
    for (std::size_t b = 0; b < H.n; b++)
      joules += 3600.0 * day[h].weight * std::max(0.0, H.to_air[b] * (A.setpoint - air[h]) + H.to_ground[b] * (A.setpoint - ground) - S[b]);
  }
  m.heat_total = joules / 3.6e6;
  m.heat = m.heat_total / built_floor;
  tick(1);
  // Winter sun on the public ground.
  const auto h = heights(d);
  const auto o = opened(d);
  const std::size_t P = shade.points.size(), K = shade.suns.size();
  double sun = 0;
  for (std::size_t p = 0; p < P; p++) {
    if (shade.owner[p] >= 0 && !o[static_cast<std::size_t>(shade.owner[p])]) continue;
    for (std::size_t k = 0; k < K; k++)
      if (sunlit(shade, p, k, h)) sun += shade.suns[k].weight;
  }
  m.sun = sun * shade.cell;
  tick(2);
  // Walking from homes to their nearest shop and open space.
  std::shared_ptr<const Walk> W;
  try {
    W = walk(d.links);
  } catch (const std::exception& e) {
    m.ok = false, m.reason = std::string("the links cannot all be made: ") + e.what();
    return m;
  }
  std::vector<OnSegment> shop_at, open_at = W->space;
  for (std::size_t i = 0; i < n; i++) {
    if (d.uses[i] == shops) shop_at.push_back(W->parcel[i]);
    if (d.uses[i] == open_space) open_at.push_back(W->parcel[i]);
  }
  double weight = 0, to_shop = 0, to_open = 0;
  if (!shop_at.empty() && !open_at.empty()) {
    const auto ds = from_sources(W->graph, shop_at), dopen = from_sources(W->graph, open_at);
    for (std::size_t i = 0; i < n; i++) {
      if (d.uses[i] != dwellings) continue;
      weight += floor[i];
      to_shop += floor[i] * to_target(W->graph, ds, shop_at, W->parcel[i]);
      to_open += floor[i] * to_target(W->graph, dopen, open_at, W->parcel[i]);
    }
  }
  if (!(weight > 0) || shop_at.empty()) { m.ok = false, m.reason = "no homes, or no shop to walk to"; return m; }
  m.access_shop = to_shop / weight, m.access_open = to_open / weight;
  m.access = m.access_shop + m.access_open;
  tick(3);
  if (!std::isfinite(m.access)) { m.ok = false, m.reason = "a home cannot reach a shop or an open space"; return m; }
  // The mix of uses: H(U | X), the floor's use given the block it stands in.
  std::vector<std::array<double, 3>> rows(place_count, {0, 0, 0});
  for (std::size_t i = 0; i < n; i++)
    if (d.uses[i] != open_space) rows[static_cast<std::size_t>(parcels[i].place)][static_cast<std::size_t>(d.uses[i])] += floor[i];
  std::vector<std::vector<double>> kept;
  for (const auto& r : rows)
    if (r[0] + r[1] + r[2] > 0) kept.push_back({r[0], r[1], r[2]});
  const auto mix = physis::statistical_mechanics::mixing(calx::linalg::Rectangular::from_rows(kept));
  m.mix = mix.conditional, m.mix_entropy = mix.entropy, m.mix_departure = mix.departure;
  tick(4);
  // The shares the constraints bound, and the descriptors.
  double open_land = open_space_area;
  for (std::size_t i = 0; i < n; i++)
    if (d.uses[i] == open_space) open_land += parcels[i].area;
  m.open_share = open_land / land;
  m.dwelling_share = m.floor_by_use[0] / built_floor;
  m.shop_share = m.floor_by_use[2] / built_floor;
  m.plot_ratio = m.floor_area / land;
  return m;
}

}  // namespace study
