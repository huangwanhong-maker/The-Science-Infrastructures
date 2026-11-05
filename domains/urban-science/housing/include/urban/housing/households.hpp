// Households as agents among the centre's dwellings: who searches, where they choose to live, and the rents their
// choices form (ADR 0075; the design, development/ideas/2026-11-04-agents-and-decisions-for-urban-planning.md §8).
//
// The city. Buildings, each in an area (a CLUE small area) and a zone (a CLUE block, which lies in one area), with a
// point and a distance to the centre's anchor, and a capacity of dwellings of each of CLUE's three types. A segment is
// a type in an area, and carries one rent. A building's neighbourhood is its zone (the households of its block), or,
// when a graph is given, itself and its neighbours in the graph (as Schelling's cells).
//
// The households. Each has an income class, a tenure (owner or renter) and a size (one person, or two or more), and
// is resident in one dwelling of one building, or searching to enter. Its kind k = (income, tenure, size).
//
// A tick (a quarter):
// 1. `begin`: the dwellings held since the last tick are offered; the segments' utilities at the tick's rents;
// 2. the caller's: the stock (`set_capacity`; a dwelling removed evicts its household, which searches at once) and the
//    arrivals (`arrive`: newcomers of each kind join as searchers);
// 3. `run`: the households, visited in random sequential order (nomos::agents::tick), each by its own stream
//    (seed, "agnt", household, tick), so that a run depends on its seed alone:
//    - a resident leaves the centre with its tenure's hazard of exit (for reasons that are not the centre's: a job
//      elsewhere, studies ended), its dwelling vacated;
//    - a resident searches with its tenure's hazard (a move considered: a lease ending, a household changed), when
//      its neighbourhood's share of its own income class falls below the threshold (Schelling's rule, if asked), or
//      when its rent exceeds the share of its income it can bear (renters);
//    - a searcher chooses among the dwellings on offer, its own (if resident, with a bonus for staying) and the
//      outside option (to leave the centre, or not to enter it), by a logit of the utility of a dwelling in building
//      b of type t,
//        u = δ_seg + γ_{k,t} − α_k ln(r_seg / r̄) − θ d_b + η · like + ε,
//      r̄ the reference rent (the outside option's: the metropolis's, set by the caller each tick), like the share of
//      its own income class among the other households of the neighbourhood (1 when there are none). Dwellings alike are chosen among as one option with the log of their number added; with zones as
//      neighbourhoods, a zone's dwellings of a type differ only by θ d_b, and the option (zone, type) carries
//      ln Σ_b n_b e^{−θ d_b}, its building then drawn in proportion to n_b e^{−θ d_b}: exactly the choice among the
//      dwellings one by one (the logit's nesting with equal scales);
//    - with no taste for neighbours (η = 0) and the logit's scale 1, an option's utility is its segment's for the kind
//      plus the log of its weight, so that the choice is made in two exact steps: a segment (or staying, or the
//      outside option) by the logit of δ_seg + γ − α ln(r/r̄) + ln W_seg, W_seg the segment's weight on offer, then an
//      option within it in proportion to its weight; a searcher then costs the number of segments, not of options;
//    - with an infinite scale, the choice is the best option (θ d_b left out), ties drawn uniformly over dwellings,
//      staying kept on a tie: with no rent, no type and a step utility η [like ≥ f], this is Schelling's model
//      (Schelling 1971, "Dynamic models of segregation", J. Math. Sociol. 1, 143–186) exactly as
//      nomos::agents::Schelling holds it, draw for draw, on a graph;
// 4. `close`: rents formed by excess demand by type and area: each segment's log rent moves by κ ln((D + 1)/(V + 1)),
//    D the tick's notional demand for the segment (the sum over its searchers of their probabilities of taking one
//    of its dwellings among all those posted at the tick's start, as if none had yet been let: the demand that a
//    market of posted vacancies meets, and that exceeds them when searchers outnumber what is posted) and V the
//    dwellings it posted, capped at ± `cap` a tick (Walras's tâtonnement in the logarithm of the
//    price, read on the flow of lettings; Arrow and Hurwicz 1958, "On the stability of the competitive equilibrium,
//    I", Econometrica 26, 522–552, for its stability under gross substitutes).
// A dwelling vacated is offered from the next tick when `letting_lag` is set (a letting takes a quarter): the vacancy
// that results is turnover times the time to let, the natural vacancy of Rosen and Smith (1983, "The price-adjustment
// process for rental housing and the natural vacancy rate", Am. Econ. Rev. 73, 779–786), an outcome of the model
// rather than a parameter.
//
// The record. Lettings are counted by tick and segment, each let at the segment's rent of the tick (Homes Victoria's
// medians are of new lettings); with a log, every move is an event of nomos::agents::Log (rule 0 entered, 1 moved,
// 2 left, 3 evicted, 4 did not enter), touching the dwellings left and taken (b · 3 + type), so that a simulated
// history is a record of events.
//
// Caveats. Rents are one per segment; sitting tenants pay the segment's rent now (no lease is modelled); owners face
// the rent as the user cost of their dwelling and no burden. Preferences are fitted to aggregates (the study), never
// observed: no dwelling's rent or move is in the data.

#pragma once

#include <urban/housing/market.hpp>

#include <nomos/agents/choice.hpp>
#include <nomos/agents/engine.hpp>
#include <nomos/agents/space.hpp>

#include <calx/random/categorical.hpp>
#include <calx/random/stream.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace urban::housing {

enum Tenure : std::uint8_t { owner = 0, renter = 1 };

struct Building {
  std::size_t area = 0, zone = 0;
  double x = 0, y = 0;     // metres, the product's local frame
  double distance = 0;     // km to the anchor
};

struct City {
  std::vector<Building> buildings;
  std::size_t areas = 1, zones = 1;
  nomos::agents::Graph graph;  // empty: neighbourhoods are zones; else a building's own and its graph neighbours'

  std::size_t segments() const { return types * areas; }
  std::size_t segment(std::size_t b, std::size_t t) const { return t * areas + buildings[b].area; }
  bool zoned() const { return graph.size() == 0; }
};

struct Household {
  std::uint8_t income = 0, tenure = renter, size = 1;  // size: 1, or 2 for two or more
  bool active = false;     // resident, or searching to enter
  bool resident = false;
  std::uint32_t building = 0;
  std::uint8_t type = 0;
  std::int32_t since = 0;  // the tick of its letting
  float rent = 0;          // the weekly rent at its letting
};

inline std::size_t kind(const Household& h) { return (static_cast<std::size_t>(h.income) * 2 + h.tenure) * 2 + (h.size >= 2 ? 1 : 0); }

struct Preferences {
  std::vector<double> base;               // δ by segment (type · areas + area)
  std::vector<double> income;             // y by income class, A$ a week
  double alpha = 0, mean_income = 1, income_power = 1;
  std::array<double, types> taste_owner{0, 0, -std::numeric_limits<double>::infinity()};  // γ by type, owners
  std::array<double, types> taste_renter{0, 0, 0};                                         // γ by type, renters
  double family_house = 0;                // added to houses for households of two or more
  bool student_single_only = true;        // student apartments only for one-person households
  double distance_cost = 0;               // θ per km
  double neighbourhood = 0;               // η
  bool step = false;                      // η [like ≥ threshold] in place of η like
  double threshold = 0;                   // f
  double scale = 1;                       // the logit's scale; +∞: the best option
  double stay = 0;                        // added to a resident's own dwelling
  double outside = 0;                     // the outside option's utility; −∞: none

  double rent_coefficient(std::size_t income_class) const {
    return alpha == 0 ? 0.0 : alpha * std::pow(mean_income / income[income_class], income_power);
  }
  double taste(std::size_t k, std::size_t t) const {
    const std::size_t tenure = (k / 2) % 2, size2 = k % 2;
    if (t == student && student_single_only && size2) return -std::numeric_limits<double>::infinity();
    double g = tenure == owner ? taste_owner[t] : taste_renter[t];
    if (t == house && size2) g += family_house;
    return g;
  }
};

struct Dynamics {
  std::array<double, 2> exit{0, 0};       // a tick: owners, renters leave the centre, whatever it offers
  std::array<double, 2> hazard{0, 0};     // a tick: owners, renters consider a move
  bool forced_by_neighbourhood = false;   // below the threshold, search
  double burden = std::numeric_limits<double>::infinity();  // renters above this rent/income search
  double adjustment = 0;                  // κ
  double cap = 0.2;                       // the largest move of a log rent in a tick
  bool letting_lag = false;               // a dwelling vacated is offered from the next tick
};

// Counts of a tick.
struct TickRecord {
  int arrivals = 0, entered = 0, refused = 0, moved = 0, left = 0, evicted = 0, searched = 0;
};

namespace rule {
constexpr std::uint32_t entered = 0, moved = 1, left = 2, evicted = 3, refused = 4;
}

struct State {
  std::vector<Household> households;
  std::vector<std::array<int, types>> capacity, occupied, held;  // by building; held: vacated, offered next tick
  std::vector<std::uint32_t> pending;  // dwellings (b · 3 + t) holding some
  std::vector<int> dwelling_offer;     // by dwelling (b · 3 + t): capacity − occupied − held, as last synchronised
  // Options: (zone, type) with zones as neighbourhoods, else (building, type); each one's dwellings on offer and their
  // weight Σ n_b e^{−θ d_b}.
  std::vector<int> offer;
  std::vector<double> weight;
  std::vector<std::uint32_t> vacant;   // options with a dwelling on offer
  std::vector<std::int32_t> slot;      // each option's place in `vacant`, −1 when it has none
  std::vector<int> members;            // by building and income class, [b C + c]
  std::vector<int> around;             // by neighbourhood (zone, or building with a graph) and income class
  std::vector<int> around_total;
  std::vector<double> rent;            // by segment, A$ a week
  double reference_rent = 1;           // r̄, the outside option's rent
  std::vector<double> demand, offered; // by segment, this tick
  std::vector<std::uint32_t> posted;   // the options on offer at the tick's start, and their weights
  std::vector<double> posted_weight;
  std::vector<double> segment_weight, posted_segment_weight;  // by segment: Σ of its options' weights, now and posted
  std::size_t classes = 1;
  std::int32_t tick = 0;
  // The record.
  std::vector<std::vector<int>> lettings;            // [tick][segment]
  std::vector<std::vector<double>> rents;            // [tick][segment], the rents at which the tick let
  std::vector<std::vector<double>> demands, offers;  // [tick][segment]
  std::vector<TickRecord> ticks;
  // Each tick's utility of each kind for each segment, without the building's terms (kept by `begin`).
  std::vector<double> segment_utility;               // [kind · J + segment]
};

class Simulation {
 public:
  using State = housing::State;

  Simulation(const City& city, Preferences pref, Dynamics dyn) : city_(&city), pref_(std::move(pref)), dyn_(dyn) {
    detail::require(pref_.base.size() == city.segments(), "one base value a segment");
    detail::require(!pref_.income.empty(), "at least one income class");
    const std::size_t B = city.buildings.size();
    exp_w_.resize(B), w_.resize(B);
    for (std::size_t b = 0; b < B; b++) w_[b] = -pref_.distance_cost * city.buildings[b].distance, exp_w_[b] = std::exp(w_[b]);
    if (city.zoned()) {
      zone_buildings_.assign(city.zones, {});
      zone_area_.assign(city.zones, city.areas);
      for (std::size_t b = 0; b < B; b++) {
        const auto& g = city.buildings[b];
        detail::require(g.zone < city.zones, "a building's zone beyond the city's");
        detail::require(zone_area_[g.zone] == city.areas || zone_area_[g.zone] == g.area, "a zone lies in one area");
        zone_area_[g.zone] = g.area, zone_buildings_[g.zone].push_back(b);
      }
    }
  }
  const City& city() const { return *city_; }
  const Preferences& preferences() const { return pref_; }
  const Dynamics& dynamics() const { return dyn_; }

  std::size_t options() const { return (city_->zoned() ? city_->zones : city_->buildings.size()) * types; }
  std::size_t option_of(std::size_t b, std::size_t t) const { return (city_->zoned() ? city_->buildings[b].zone : b) * types + t; }
  std::size_t option_segment(std::size_t o) const {
    return city_->zoned() ? (o % types) * city_->areas + zone_area_[o / types] : city_->segment(o / types, o % types);
  }

  // An empty city at its capacities and rents.
  State empty(const std::vector<std::array<int, types>>& capacity, std::vector<double> rents) const {
    const std::size_t B = city_->buildings.size(), C = pref_.income.size();
    detail::require(capacity.size() == B, "one capacity a building");
    detail::require(rents.size() == city_->segments(), "one rent a segment");
    State s;
    s.classes = C;
    s.capacity = capacity;
    s.occupied.assign(B, {0, 0, 0}), s.held.assign(B, {0, 0, 0});
    s.dwelling_offer.assign(B * types, 0);
    s.offer.assign(options(), 0), s.weight.assign(options(), 0.0), s.slot.assign(options(), -1);
    s.members.assign(B * C, 0);
    const std::size_t N = city_->zoned() ? city_->zones : B;
    s.around.assign(N * C, 0), s.around_total.assign(N, 0);
    s.rent = std::move(rents);
    for (std::size_t b = 0; b < B; b++)
      for (std::size_t t = 0; t < types; t++) sync(s, b, t);
    return s;
  }

  // Places household h (not resident) in building b's type t (for initial states).
  void place(State& s, std::size_t h, std::size_t b, std::size_t t) const { occupy(s, s.households[h], b, t); }

  // Sets the capacity of building b in type t; a dwelling removed that is occupied evicts its latest household, which
  // searches this tick.
  void set_capacity(State& s, std::size_t b, std::size_t t, int n, nomos::agents::Log* log = nullptr) const {
    s.capacity[b][t] = n;
    for (std::size_t h = s.households.size(); h-- > 0 && s.occupied[b][t] > n;) {
      auto& H = s.households[h];
      if (H.resident && H.building == b && H.type == t) {
        vacate(s, H, false);
        if (!s.ticks.empty()) s.ticks.back().evicted++;
        if (log) log->add(s.tick, h, rule::evicted, {b * types + t});
      }
    }
    s.held[b][t] = std::min(s.held[b][t], std::max(0, n - s.occupied[b][t]));
    sync(s, b, t);
  }

  // Newcomers join as searchers.
  void arrive(State& s, std::span<const Household> newcomers) const {
    for (Household h : newcomers) {
      h.active = true, h.resident = false;
      s.households.push_back(h);
    }
    if (!s.ticks.empty()) s.ticks.back().arrivals += static_cast<int>(newcomers.size());
  }

  // Opens a tick: its record, the dwellings held since the last offered, the segments' utilities at its rents.
  void begin(State& s) const {
    s.ticks.push_back({});
    for (const std::uint32_t d : s.pending) {
      s.held[d / types][d % types] = 0;
      sync(s, d / types, d % types);
    }
    s.pending.clear();
    const std::size_t J = city_->segments(), K = pref_.income.size() * 4;
    s.demand.assign(J, 0.0), s.offered.assign(J, 0.0);
    s.lettings.push_back(std::vector<int>(J, 0));
    s.segment_weight.assign(J, 0.0);
    for (const std::uint32_t o : s.vacant) s.segment_weight[option_segment(o)] += s.weight[o];  // afresh: no drift
    s.rents.push_back(s.rent);
    s.segment_utility.assign(K * J, 0.0);
    for (std::size_t k = 0; k < K; k++)
      for (std::size_t j = 0; j < J; j++) {
        const std::size_t t = j / city_->areas;
        const double g = pref_.taste(k, t);
        s.segment_utility[k * J + j] = std::isfinite(g) ? pref_.base[j] + g - pref_.rent_coefficient(k / 4) * std::log(s.rent[j] / s.reference_rent) : -std::numeric_limits<double>::infinity();
      }
  }

  // The households' choices of the tick, in random sequential order (nomos), then the rents.
  void run(State& s, std::uint64_t seed, nomos::agents::Log* log = nullptr) const {
    s.posted = s.vacant, s.posted_weight.clear();
    s.posted_segment_weight = s.segment_weight;
    for (const std::uint32_t o : s.vacant) s.offered[option_segment(o)] += s.offer[o], s.posted_weight.push_back(s.weight[o]);
    nomos::agents::tick(*this, s, nomos::agents::Schedule::random_sequential, seed, static_cast<std::uint32_t>(s.tick), log);
    close(s);
  }

  // The rents' adjustment by the tick's excess demand, and the tick's record.
  void close(State& s) const {
    const std::size_t J = city_->segments();
    std::vector<double> cap(J, 0.0);
    for (std::size_t b = 0; b < city_->buildings.size(); b++)
      for (std::size_t t = 0; t < types; t++) cap[city_->segment(b, t)] += s.capacity[b][t];
    if (dyn_.adjustment != 0)
      for (std::size_t j = 0; j < J; j++)
        if (cap[j] > 0) s.rent[j] *= std::exp(std::clamp(dyn_.adjustment * std::log((s.demand[j] + 1) / (s.offered[j] + 1)), -dyn_.cap, dyn_.cap));
    s.demands.push_back(s.demand), s.offers.push_back(s.offered);
    s.tick++;
  }

  // One whole tick with no change of stock and no arrivals.
  void step(State& s, std::uint64_t seed, nomos::agents::Log* log = nullptr) const {
    begin(s);
    run(s, seed, log);
  }

  // The logit's options of household i at a finite scale: the options on offer it may take, its own dwelling (if
  // resident) and the outside option (if any), in that order, with their scaled utilities (an option's including the
  // log of its weight), so that P = nomos::agents::logit(utility, 1).
  struct Choices {
    std::vector<std::uint32_t> option;
    std::vector<double> utility;
    bool stay = false, outside = false;
  };
  Choices choices(const State& s, std::size_t i) const {
    const Household& H = s.households[i];
    const std::size_t J = city_->segments();
    const double* seg = &s.segment_utility[kind(H) * J];
    Choices out;
    out.utility.reserve(s.vacant.size() + 2), out.option.reserve(s.vacant.size());
    for (const std::uint32_t o : s.vacant) {
      const double base = seg[option_segment(o)];
      if (!std::isfinite(base)) continue;
      double u = base + std::log(s.weight[o]);
      if (pref_.neighbourhood != 0) u += pref_.neighbourhood * like(s, neighbourhood_of_option(o), H.income, H.resident, H.building);
      out.utility.push_back(pref_.scale * u), out.option.push_back(o);
    }
    if (H.resident) out.stay = true, out.utility.push_back(pref_.scale * (own_utility(s, H, seg) + pref_.stay));
    if (std::isfinite(pref_.outside)) out.outside = true, out.utility.push_back(pref_.scale * pref_.outside);
    return out;
  }

  // Adds a searcher's probabilities over the options posted at the tick's start to the segments' demand.
  void notional_demand(State& s, const Household& H) const {
    const std::size_t J = city_->segments();
    const double* seg = &s.segment_utility[kind(H) * J];
    std::vector<double> u;
    std::vector<std::size_t> to;
    u.reserve(s.posted.size() + 2), to.reserve(s.posted.size());
    for (std::size_t k = 0; k < s.posted.size(); k++) {
      const std::uint32_t o = s.posted[k];
      const double base = seg[option_segment(o)];
      if (!std::isfinite(base)) continue;
      double x = base + std::log(s.posted_weight[k]);
      if (pref_.neighbourhood != 0) x += pref_.neighbourhood * like(s, neighbourhood_of_option(o), H.income, H.resident, H.building);
      u.push_back(pref_.scale * x), to.push_back(option_segment(o));
    }
    if (H.resident) u.push_back(pref_.scale * (own_utility(s, H, seg) + pref_.stay));
    if (std::isfinite(pref_.outside)) u.push_back(pref_.scale * pref_.outside);
    if (u.empty()) return;
    const auto P = nomos::agents::logit(u, 1.0);
    for (std::size_t k = 0; k < to.size(); k++) s.demand[to[k]] += P[k];
  }

  // Whether choices are made in two exact steps (no taste for neighbours, the logit's scale 1).
  bool two_step() const { return pref_.neighbourhood == 0 && pref_.scale == 1; }

  // With two_step(): household H's probabilities of each segment (its options now, or those posted at the tick's start),
  // then of staying (if resident) and of the outside option (if any): J + 2 entries, zero where absent.
  std::vector<double> segment_choice(const State& s, const Household& H, bool posted) const {
    const std::size_t J = city_->segments();
    const double* seg = &s.segment_utility[kind(H) * J];
    const auto& W = posted && s.posted_segment_weight.size() == J ? s.posted_segment_weight : s.segment_weight;  // posted: as `run` kept them
    std::vector<double> u(J + 2, -std::numeric_limits<double>::infinity());
    for (std::size_t j = 0; j < J; j++)
      if (W[j] > 0 && std::isfinite(seg[j])) u[j] = seg[j] + std::log(W[j]);
    if (H.resident) u[J] = own_utility(s, H, seg) + pref_.stay;
    if (std::isfinite(pref_.outside)) u[J + 1] = pref_.outside;
    double top = -std::numeric_limits<double>::infinity();
    for (const double x : u) top = std::max(top, x);
    std::vector<double> P(J + 2, 0.0);
    if (!std::isfinite(top)) return P;
    double total = 0;
    for (std::size_t k = 0; k < u.size(); k++) total += (P[k] = std::isfinite(u[k]) ? std::exp(u[k] - top) : 0.0);
    for (double& x : P) x /= total;
    return P;
  }

  // Drops the households no longer active (who left, or did not enter), keeping the others' order: their indices, and
  // so their streams, are renumbered, deterministically. For long runs whose newcomers mostly look elsewhere.
  void compact(State& s) const {
    std::size_t k = 0;
    for (std::size_t i = 0; i < s.households.size(); i++)
      if (s.households[i].active) s.households[k++] = s.households[i];
    s.households.resize(k);
  }

  // ---- nomos's contract ----
  std::size_t agents(const State& s) const { return s.households.size(); }

  void update(std::size_t i, const State&, State& s, calx::random::Stream& rng, nomos::agents::Log* log) const {
    Household& H = s.households[i];
    if (!H.active) return;
    const bool resident = H.resident;
    const std::size_t J = city_->segments();
    const double* seg = &s.segment_utility[kind(H) * J];
    if (resident) {
      const double e = dyn_.exit[H.tenure];
      if (e > 0 && rng.uniform() < e) {
        const std::size_t from = H.building * types + H.type;
        vacate(s, H, dyn_.letting_lag);
        H.active = false;
        s.ticks.back().left++;
        if (log) log->add(s.tick, i, rule::left, {from});
        return;
      }
      bool search = false;
      const double h = dyn_.hazard[H.tenure];
      if (h > 0 && rng.uniform() < h) search = true;
      if (!search && dyn_.forced_by_neighbourhood && like(s, neighbourhood_of_building(H.building), H.income, true, H.building) < pref_.threshold) search = true;
      if (!search && H.tenure == renter && std::isfinite(dyn_.burden) && s.rent[city_->segment(H.building, H.type)] > dyn_.burden * pref_.income[H.income])
        search = true;
      if (!search) return;
    }
    s.ticks.back().searched++;
    std::uint32_t taken = 0;   // the option taken
    long within = -1;          // with the best option: which of its dwellings on offer (in its buildings' order)
    bool take = false, leave = false;
    if (std::isinf(pref_.scale)) {
      // The best option; ties uniformly over dwellings; staying kept on a tie.
      const std::size_t own = resident ? option_of(H.building, H.type) : options();
      const auto best_utility = [&](std::uint32_t o) {
        double u = seg[option_segment(o)];
        if (std::isfinite(u) && pref_.neighbourhood != 0) {
          const double l = like(s, neighbourhood_of_option(o), H.income, resident, H.building);
          u += pref_.neighbourhood * (pref_.step ? (l >= pref_.threshold ? 1.0 : 0.0) : l);
        }
        return u;
      };
      double best = -std::numeric_limits<double>::infinity();
      for (const std::uint32_t o : s.vacant)
        if (!(o == own && !city_->zoned())) best = std::max(best, best_utility(o));
      double current = -std::numeric_limits<double>::infinity();
      if (resident) {
        current = seg[city_->segment(H.building, H.type)];
        if (pref_.neighbourhood != 0) {
          const double l = like(s, neighbourhood_of_building(H.building), H.income, true, H.building);
          current += pref_.neighbourhood * (pref_.step ? (l >= pref_.threshold ? 1.0 : 0.0) : l);
        }
        current += pref_.stay;
      }
      if (!(best > std::max(current, pref_.outside))) {
        leave = !resident || pref_.outside > current;
      } else {
        std::vector<std::uint32_t> tied;
        for (const std::uint32_t o : s.vacant)
          if (!(o == own && !city_->zoned()) && best_utility(o) == best) tied.push_back(o);
        std::sort(tied.begin(), tied.end());
        long total = 0;
        for (const std::uint32_t o : tied) total += s.offer[o];
        long pick = static_cast<long>(rng.below(static_cast<std::uint64_t>(total)));
        for (const std::uint32_t o : tied) {
          if (pick < s.offer[o]) { taken = o, within = pick, take = true; break; }
          pick -= s.offer[o];
        }
        s.demand[option_segment(taken)] += 1;
      }
    } else {
      if (two_step()) {
        // The notional demand, over what was posted; then the choice, a segment and an option within it.
        const auto Pn = segment_choice(s, H, true);
        for (std::size_t j = 0; j < J; j++) s.demand[j] += Pn[j];
        const auto P = segment_choice(s, H, false);
        const std::size_t c = calx::random::categorical(rng, P);
        if (c < J) {
          double u = rng.uniform() * s.segment_weight[c];
          for (const std::uint32_t o : s.vacant) {
            if (option_segment(o) != c) continue;
            taken = o, take = true;
            if ((u -= s.weight[o]) < 0) break;
          }
        } else if (!(resident && c == J)) {
          leave = true;
        }
      } else {
      const Choices C = choices(s, i);
      const std::size_t n = C.option.size();
      if (C.utility.empty()) {
        leave = !resident;
      } else {
        notional_demand(s, H);
        const auto P = nomos::agents::logit(C.utility, 1.0);
        const std::size_t c = calx::random::categorical(rng, P);
        if (c < n) taken = C.option[c], take = true;
        else if (!(C.stay && c == n)) leave = true;
      }
      }
    }
    auto& rec = s.ticks.back();
    if (take) {
      const std::size_t t = taken % types;
      const std::size_t b = city_->zoned() ? building_in(s, taken / types, t, within, rng) : taken / types;
      const std::size_t from = resident ? H.building * types + H.type : 0;
      if (resident) vacate(s, H, dyn_.letting_lag);
      occupy(s, H, b, t);
      if (resident) rec.moved++; else rec.entered++;
      if (log) {
        if (resident) log->add(s.tick, i, rule::moved, {from, b * types + t});
        else log->add(s.tick, i, rule::entered, {b * types + t});
      }
    } else if (leave) {
      if (resident) {
        const std::size_t from = H.building * types + H.type;
        vacate(s, H, dyn_.letting_lag);
        rec.left++;
        if (log) log->add(s.tick, i, rule::left, {from});
      } else {
        rec.refused++;
        if (log) log->add(s.tick, i, rule::refused, {});
      }
      H.active = false;
    }
  }

  // The share of income class c among the other households of neighbourhood n (1 when it holds none), the household
  // itself (in building `self`, if resident) not counted.
  double like(const State& s, std::size_t n, std::size_t c, bool resident, std::size_t self) const {
    const std::size_t C = s.classes;
    int same = s.around[n * C + c], all = s.around_total[n];
    if (city_->zoned()) {
      if (resident && city_->buildings[self].zone == n) same--, all--;
    } else {
      for (std::size_t q = 0; q < C; q++) all += s.members[n * C + q];
      same += s.members[n * C + c];
      if (resident && (self == n || std::binary_search(city_->graph.neighbours(n).begin(), city_->graph.neighbours(n).end(), self))) same--, all--;
    }
    return all <= 0 ? 1.0 : static_cast<double>(same) / all;
  }
  std::size_t neighbourhood_of_building(std::size_t b) const { return city_->zoned() ? city_->buildings[b].zone : b; }
  std::size_t neighbourhood_of_option(std::size_t o) const { return o / types; }

  // A building's dwellings of type t on offer: capacity − occupied − held.
  int on_offer(const State& s, std::size_t b, std::size_t t) const { return s.capacity[b][t] - s.occupied[b][t] - s.held[b][t]; }

 private:
  double own_utility(const State& s, const Household& H, const double* seg) const {
    double u = seg[city_->segment(H.building, H.type)] + w_[H.building];
    if (pref_.neighbourhood != 0) u += pref_.neighbourhood * like(s, neighbourhood_of_building(H.building), H.income, true, H.building);
    return u;
  }

  // The building of a zone's option: the `within`-th dwelling on offer (in the zone's buildings' order) when given,
  // else drawn in proportion to n_b e^{−θ d_b}.
  std::size_t building_in(const State& s, std::size_t z, std::size_t t, long within, calx::random::Stream& rng) const {
    const auto& list = zone_buildings_[z];
    if (within >= 0) {
      for (const std::size_t b : list) {
        const int n = std::max(0, on_offer(s, b, t));
        if (within < n) return b;
        within -= n;
      }
    } else {
      double u = rng.uniform() * s.weight[z * types + t];
      std::size_t last = list.front();
      for (const std::size_t b : list) {
        const int n = std::max(0, on_offer(s, b, t));
        if (n == 0) continue;
        last = b;
        u -= n * exp_w_[b];
        if (u < 0) return b;
      }
      return last;
    }
    detail::require(false, "an option without a dwelling on offer");
    return 0;
  }

  // Brings dwelling (b, t)'s offer into its option.
  void sync(State& s, std::size_t b, std::size_t t) const {
    const std::size_t d = b * types + t, o = option_of(b, t);
    const int now = std::max(0, on_offer(s, b, t)), delta = now - s.dwelling_offer[d];
    if (delta == 0) return;
    s.dwelling_offer[d] = now;
    s.offer[o] += delta;
    const double before = s.weight[o];
    s.weight[o] = s.offer[o] == 0 ? 0.0 : s.weight[o] + delta * exp_w_[b];
    if (!s.segment_weight.empty()) s.segment_weight[option_segment(o)] = std::max(0.0, s.segment_weight[option_segment(o)] + s.weight[o] - before);
    if (s.offer[o] > 0 && s.slot[o] < 0) s.slot[o] = static_cast<std::int32_t>(s.vacant.size()), s.vacant.push_back(static_cast<std::uint32_t>(o));
    if (s.offer[o] == 0 && s.slot[o] >= 0) {
      const std::int32_t at = s.slot[o];
      const std::uint32_t last = s.vacant.back();
      s.vacant[static_cast<std::size_t>(at)] = last, s.slot[last] = at;
      s.vacant.pop_back(), s.slot[o] = -1;
    }
  }

  void move_members(State& s, std::size_t b, std::size_t c, int d) const {
    const std::size_t C = s.classes;
    s.members[b * C + c] += d;
    if (city_->zoned()) {
      const std::size_t z = city_->buildings[b].zone;
      s.around[z * C + c] += d, s.around_total[z] += d;
    } else {
      for (const std::size_t m : city_->graph.neighbours(b)) s.around[m * C + c] += d, s.around_total[m] += d;
    }
  }

  void vacate(State& s, Household& H, bool hold) const {
    const std::size_t b = H.building, t = H.type;
    s.occupied[b][t]--;
    move_members(s, b, H.income, -1);
    H.resident = false;
    if (hold) {
      if (s.held[b][t]++ == 0) s.pending.push_back(static_cast<std::uint32_t>(b * types + t));
    }
    sync(s, b, t);
  }

  void occupy(State& s, Household& H, std::size_t b, std::size_t t) const {
    s.occupied[b][t]++;
    move_members(s, b, H.income, +1);
    H.resident = true, H.active = true, H.building = static_cast<std::uint32_t>(b), H.type = static_cast<std::uint8_t>(t), H.since = s.tick;
    const std::size_t j = city_->segment(b, t);
    H.rent = static_cast<float>(s.rent[j]);
    if (!s.lettings.empty()) s.lettings.back()[j]++;
    sync(s, b, t);
  }

  const City* city_;
  Preferences pref_;
  Dynamics dyn_;
  std::vector<double> w_, exp_w_;                    // −θ d_b and its exponential
  std::vector<std::vector<std::size_t>> zone_buildings_;
  std::vector<std::size_t> zone_area_;
};

// ---- measures of a state ------------------------------------------------------------------------------------------

// The occupied dwellings and the capacity by segment.
struct SegmentCounts {
  std::vector<double> capacity, occupied;
};
inline SegmentCounts segment_counts(const City& city, const State& s) {
  SegmentCounts out{std::vector<double>(city.segments(), 0.0), std::vector<double>(city.segments(), 0.0)};
  for (std::size_t b = 0; b < city.buildings.size(); b++)
    for (std::size_t t = 0; t < types; t++) {
      out.capacity[city.segment(b, t)] += s.capacity[b][t];
      out.occupied[city.segment(b, t)] += s.occupied[b][t];
    }
  return out;
}

}  // namespace urban::housing
