// Owners and occupiers of a district's buildings (the infrastructure's ADR 0071 and ADR 0075; the design,
// development/ideas/2026-11-04-agents-and-decisions-for-urban-planning.md, §10; theory/occupancy.md, §2-3): an
// agent-based model whose baseline law is the nested chain of the district's uses (physis decision 0020, fitted by
// the module `events`), modulated by the neighbours' uses (the module `affinity`'s co-location) and by demand through
// a channel of vacancy.
//
// The agents. Each building is held by an owner and used by an occupier. Its state is its use, one of a
// classification's states (CLUE's space uses in their families: urban::events::prefix_classification). Its rules:
// - the occupier is succeeded at level l, at rate λ_l (l = 0: the building leaves its family; l = 1: it stays in its
//   family, a shop becoming a showroom, an unused building put under renovation): the owner lets it again within its
//   level-l block, choosing the next use by a logit (nomos::agents::logit, McFadden 1974) of the utilities
//       u_i(y) = ln π(y) + β_n A_i(y),
//   π the classification's law, A_i(y) = (1/|N_i|) Σ_{j ∈ N_i} a(u_j, y) the mean affinity of the use y to the
//   neighbours' uses (a(u, y) = ln CLQ_{u→y}, the co-location quotient of y about u; Leslie and Kronenfeld 2011). A
//   redraw that returns the use held is not a change;
// - the occupier of an occupied building leaves with no successor at rate ν e^{−β_d d(t)}: the building falls into the
//   vacant family, unused (V, CLUE's "Unoccupied - Unused") or under works (W: renovation, construction, demolition),
//   by the same logit within that family;
// - the owner lets a vacant building at rate η e^{β_d d(t)}: with probability ρ within the family it last held (the
//   owner's fit-out, which the census does not see: kept in the state as `held`), else among all the occupied uses,
//   each by the same logit;
// d(t) the year's demand. The last two are the channel of vacancy: the record shows buildings falling vacant and let
// about three times as often as the chain's family memory implies, and more than half of its vacancies ending in the
// family that preceded them (theory/occupancy.md, §4.3).
//
// Proposition 1 (the baseline). With ν = η = 0, β_n = 0 and no incentive (ρ then acts on nothing), each building's
// generator is
// Σ_l λ_l (Π_l − I), physis's nested chain on the classification, exactly (demand then acts on nothing); the buildings
// are independent, and the law of a panel is the product of the chain's transitions (theory/occupancy.md,
// Proposition 2.2). The check holds the generator to physis's, entry by entry.
//
// Proposition 2 (three lumps). Let the uses fall in three lumps: occupied, vacant and lettable (V), vacant under works
// (W: the other unoccupied uses). If the vacant uses are exactly one block at level 1 (CLUE's family Unoccupied), then
// for any fixed neighbourhood and demand the model is strongly lumpable onto {O, V, W} (Kemeny and Snell 1960, Theorem
// 6.3.2): a redraw at level 0 lands in a lump with a probability that does not depend on the use held, one at level 1
// never leaves the lump of an occupied use nor the vacant family, the channel of vacancy's rates depend on the lump
// alone, and a letting lands in O whatever the family held. The lumped generator is `lumped_generator`; the planner's
// model is built on it (planning.hpp).
//
// The incentive. An enrolled lettable vacancy (V; not one under works) is let at the extra rate κ η e^{β_d d(t)}, by
// the same law as a letting (κ an assumption declared by the caller: κ = 1 doubles its letting).
//
// The schedule. Asynchronous in continuous time (nomos's Gillespie): the law of the configuration solves its master
// equation, and with the modulation off each building's chain is the nested chain exactly. The demand is constant
// within a year, so a run is a sequence of yearly runs, each from the year's own stream (its seed mixed with the year).
// Every firing appends a Change to the state, aligned with nomos's event log, so that a run's history is exact.
//
// Units: time in years; rates a year. Positions in local metres (x east, y north).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/random/categorical.hpp>
#include <calx/random/stream.hpp>
#include <nomos/agents/choice.hpp>
#include <nomos/agents/engine.hpp>
#include <physis/statistical_mechanics/nested_chain.hpp>
#include <urban/affinity/colocation.hpp>
#include <urban/events/memory.hpp>

namespace urban::occupancy {

namespace sm = physis::statistical_mechanics;
namespace ue = urban::events;

// ---- the uses -------------------------------------------------------------------------------------------------------

enum Lump : std::size_t { occupied = 0, vacant = 1, works = 2 };
constexpr std::size_t lumps = 3;

struct Uses {
  ue::UseClassification U;          // the classification: its levels, law and blocks
  std::vector<Lump> lump;           // by use
  std::vector<std::size_t> family;  // by use: its block at level 1
  std::size_t lettable = std::numeric_limits<std::size_t>::max();  // the vacant, lettable use (V)
  std::size_t states() const { return U.uses.size(); }
  std::size_t levels() const { return U.C.levels(); }
  bool is_vacant(std::size_t y) const { return lump[y] != occupied; }
};

// The uses of a classification with its vacant family (every use of that level-1 block is vacant) and the one vacant
// use that can be let as it stands. Refused when the vacant uses are not one block at level 1 (Proposition 2 then
// fails) or the lettable use is not among them.
inline Uses uses_of(ue::UseClassification U, const std::string& vacant_family = "Unoccupied", const std::string& lettable = "Unoccupied - Unused") {
  Uses S;
  S.U = std::move(U);
  const std::size_t N = S.U.uses.size();
  if (S.U.C.block.size() < 3) throw std::invalid_argument("urban::occupancy::uses_of: a classification with a level between the whole and the uses");
  S.lump.assign(N, occupied);
  S.family.assign(N, 0);
  std::size_t vac_block = std::numeric_limits<std::size_t>::max();
  for (std::size_t y = 0; y < N; y++) {
    S.family[y] = S.U.C.block[1][y];
    if (S.U.block_names[1][y] == vacant_family) vac_block = S.family[y];
  }
  if (vac_block == std::numeric_limits<std::size_t>::max()) throw std::invalid_argument("urban::occupancy::uses_of: no family \"" + vacant_family + "\"");
  for (std::size_t y = 0; y < N; y++)
    if (S.family[y] == vac_block) S.lump[y] = S.U.uses[y] == lettable ? vacant : works;
  const auto it = S.U.index.find(lettable);
  if (it == S.U.index.end() || S.lump[it->second] != vacant) throw std::invalid_argument("urban::occupancy::uses_of: the lettable use \"" + lettable + "\" is not in the vacant family");
  S.lettable = it->second;
  return S;
}

// ---- the district ---------------------------------------------------------------------------------------------------

struct Space {
  std::vector<std::array<double, 2>> at;               // each building's point (m)
  std::vector<std::vector<std::size_t>> neighbours;    // N_i
};

// The k nearest buildings of each (ties by index), within `within` metres (0: no limit).
inline Space nearest(std::vector<std::array<double, 2>> at, std::size_t k, double within = 0) {
  Space S;
  S.at = std::move(at);
  const std::size_t n = S.at.size();
  S.neighbours.resize(n);
  std::vector<std::pair<double, std::size_t>> d;
  for (std::size_t i = 0; i < n; i++) {
    d.clear();
    for (std::size_t j = 0; j < n; j++)
      if (j != i) {
        const double r = std::hypot(S.at[i][0] - S.at[j][0], S.at[i][1] - S.at[j][1]);
        if (within <= 0 || r <= within) d.push_back({r, j});
      }
    const std::size_t m = std::min(k, d.size());
    std::partial_sort(d.begin(), d.begin() + static_cast<std::ptrdiff_t>(m), d.end());
    for (std::size_t q = 0; q < m; q++) S.neighbours[i].push_back(d[q].second);
  }
  return S;
}

// The affinity of uses a(u, y) = ln CLQ_{F(u)→F(y)} of their families, from urban::affinity on a configuration (each
// building's point and use), clamped to [−cap, cap]; a pair the configuration cannot estimate (a family absent) is 0.
// Row-major states × states.
inline std::vector<double> family_affinity(const Uses& S, const Space& X, const std::vector<std::size_t>& use, double radius, double cap = 2.0,
                                           std::size_t permutations = 99, std::uint64_t seed = 1) {
  const std::size_t N = S.states();
  std::size_t F = 0;
  for (auto f : S.family) F = std::max(F, f + 1);
  std::vector<double> x, y;
  std::vector<std::size_t> mark;
  for (std::size_t i = 0; i < use.size(); i++) {
    if (use[i] >= N) continue;
    x.push_back(X.at[i][0]), y.push_back(X.at[i][1]), mark.push_back(S.family[use[i]]);
  }
  std::vector<std::string> names(F);
  for (std::size_t u = 0; u < N; u++) names[S.family[u]] = S.U.block_names[1][u];
  const auto A = urban::affinity::affinity(x, y, mark, names, 1e6, radius, permutations, seed);
  std::vector<double> a(N * N, 0.0);
  for (std::size_t u = 0; u < N; u++)
    for (std::size_t v = 0; v < N; v++) {
      const double l = A.a[S.family[u] * F + S.family[v]];
      a[u * N + v] = std::isfinite(l) ? std::clamp(l, -cap, cap) : 0.0;
    }
  return a;
}

// ---- the model ------------------------------------------------------------------------------------------------------

struct Parameters {
  std::vector<double> lambda;    // λ_l a year, l = 0 … L − 1 (the occupier succeeded at each level)
  double vacate = 0;             // ν a year: an occupier leaving with no successor (at mean demand)
  double let = 0;                // η a year: a vacancy let (at mean demand)
  double returning = 0;          // ρ: a letting's chance of returning to the family last held
  double beta_neighbours = 0;    // β_n
  double beta_demand = 0;        // β_d
  double incentive = 1;          // κ: an enrolled lettable vacancy's extra letting, in units of η e^{β_d d}
};

struct Change {
  std::size_t building = 0, from = 0, to = 0;
  std::uint32_t rule = 0;  // the level redrawn (0 … L − 1); L vacated, L + 1 let, L + 2 let by the incentive
};

struct State {
  std::vector<std::size_t> use;   // by building
  std::vector<char> alive;        // the building stands and is held (between its first and last census)
  std::vector<char> enrolled;     // in the incentive this year
  std::vector<std::size_t> held;  // the family a building last held occupied (none: the classification's states())
  double demand = 0;              // d(t), the year's
  std::vector<Change> changes;    // one per firing, aligned with nomos's event log
};

struct Model {
  const Uses* uses = nullptr;
  const Space* space = nullptr;
  std::vector<double> affinity;  // a(u, y), row-major; empty: no modulation by neighbours
  Parameters theta;

  using State = urban::occupancy::State;
  std::size_t agents(const State& s) const { return s.use.size(); }
  std::size_t rules() const { return uses->levels() + 3; }

  double rate(const State& s, std::size_t i, std::size_t r) const {
    if (!s.alive[i]) return 0;
    const std::size_t L = uses->levels();
    if (r < L) return theta.lambda[r];
    const Lump l = uses->lump[s.use[i]];
    if (r == L) return l == occupied && theta.vacate > 0 ? theta.vacate * std::exp(-theta.beta_demand * s.demand) : 0.0;
    const double let = theta.let > 0 ? theta.let * std::exp(theta.beta_demand * s.demand) : 0.0;
    if (r == L + 1) return l != occupied ? let : 0.0;
    return s.enrolled[i] && l == vacant ? theta.incentive * let : 0.0;
  }

  // A_i(y) for the configuration s.
  double attraction(const State& s, std::size_t i, std::size_t y) const {
    if (affinity.empty() || theta.beta_neighbours == 0) return 0;
    const auto& nb = space->neighbours[i];
    const std::size_t N = uses->states();
    double sum = 0;
    std::size_t k = 0;
    for (const std::size_t j : nb)
      if (s.alive[j]) sum += affinity[s.use[j] * N + y], k++;
    return k ? sum / static_cast<double>(k) : 0.0;
  }

  // The logit of u_i(y) over the uses that `in` admits, zero elsewhere.
  template <class In>
  std::vector<double> logit_over(const State& s, std::size_t i, In&& in) const {
    const auto& C = uses->U.C;
    const std::size_t N = uses->states();
    std::vector<std::size_t> members;
    for (std::size_t y = 0; y < N; y++)
      if (in(y)) members.push_back(y);
    std::vector<double> u(members.size());
    for (std::size_t k = 0; k < members.size(); k++) u[k] = std::log(C.law[members[k]]) + theta.beta_neighbours * attraction(s, i, members[k]);
    std::vector<double> out(N, 0.0);
    if (members.empty()) return out;
    const auto p = nomos::agents::logit(u, 1.0);
    for (std::size_t k = 0; k < members.size(); k++) out[members[k]] = p[k];
    return out;
  }

  // The law over the uses when rule r fires at building i: a redraw at level r < L, within the block; vacated (r = L),
  // within the vacant family; let (r > L), with chance ρ within the family last held, else among the occupied uses.
  std::vector<double> law(const State& s, std::size_t i, std::size_t r) const {
    const auto& C = uses->U.C;
    const std::size_t N = uses->states(), L = uses->levels(), x = s.use[i];
    if (r < L) return logit_over(s, i, [&](std::size_t y) { return C.block[r][y] == C.block[r][x]; });
    if (r == L) return logit_over(s, i, [&](std::size_t y) { return uses->lump[y] != occupied; });
    auto out = logit_over(s, i, [&](std::size_t y) { return uses->lump[y] == occupied; });
    const std::size_t f = s.held.empty() ? N : s.held[i];
    if (theta.returning > 0 && f < N) {
      const auto back = logit_over(s, i, [&](std::size_t y) { return uses->lump[y] == occupied && uses->family[y] == f; });
      double mass = 0;
      for (double v : back) mass += v;
      if (mass > 0)
        for (std::size_t y = 0; y < N; y++) out[y] = (1 - theta.returning) * out[y] + theta.returning * back[y];
    }
    return out;
  }

  void fire(State& s, std::size_t i, std::size_t r, calx::random::Stream& rng, std::vector<std::size_t>& affected) const {
    const auto p = law(s, i, r);
    const std::size_t y = calx::random::categorical(rng, p);
    s.changes.push_back({i, s.use[i], y, static_cast<std::uint32_t>(r)});
    s.use[i] = y;
    if (uses->lump[y] == occupied && !s.held.empty()) s.held[i] = uses->family[y];
    affected.push_back(i);  // rates depend on a building's own use alone; the neighbours' laws are read when they fire
  }
};

// Building i's generator over the uses with the rest of the configuration held: Σ_r rate_r [law_r(x, ·) − δ_x]. Row x
// is the generator's row with i in use x (and enrolled as s says).
inline calx::linalg::Matrix generator(const Model& M, State s, std::size_t i) {
  const std::size_t N = M.uses->states();
  calx::linalg::Matrix G(N);
  s.alive[i] = 1;
  for (std::size_t x = 0; x < N; x++) {
    s.use[i] = x;
    for (std::size_t r = 0; r < M.rules(); r++) {
      const double k = M.rate(s, i, r);
      if (k == 0) continue;
      const auto p = M.law(s, i, r);
      for (std::size_t y = 0; y < N; y++) G(x, y) += k * (p[y] - (y == x ? 1.0 : 0.0));
    }
  }
  return G;
}

// The lumped generator onto {O, V, W} of building i with the rest held (Proposition 2), from a representative use of
// each lump; `enrolled` its enrolment.
inline std::array<std::array<double, lumps>, lumps> lumped_generator(const Model& M, State s, std::size_t i, bool enrolled) {
  const std::size_t N = M.uses->states();
  s.enrolled[i] = enrolled ? 1 : 0;
  const auto G = generator(M, s, i);
  std::array<std::array<double, lumps>, lumps> Q{};
  for (std::size_t a = 0; a < lumps; a++) {
    std::size_t x = N;
    for (std::size_t y = 0; y < N && x == N; y++)
      if (M.uses->lump[y] == a) x = y;
    if (x == N) throw std::invalid_argument("urban::occupancy::lumped_generator: a lump with no use");
    for (std::size_t y = 0; y < N; y++) Q[a][M.uses->lump[y]] += G(x, y);
  }
  return Q;
}

// ---- runs over years ------------------------------------------------------------------------------------------------

// A building's life in the record: its first and last census, and the years it is observed (CLUE's census years).
struct Life {
  int first = 0, last = 0;
  std::size_t initial = 0;     // its use at its first census
  std::vector<int> censuses;   // ascending, first and last among them
};

// The seed of a year's stream: the run's seed mixed with the year (splitmix64's finaliser; Steele, Lea and Flood 2014).
inline std::uint64_t year_seed(std::uint64_t seed, int year) {
  std::uint64_t z = seed + 0x9e3779b97f4a7c15ull * (static_cast<std::uint64_t>(year) + 1);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
  return z ^ (z >> 31);
}

// A timed change.
struct Timed {
  double time = 0;  // years (the census of year y is at y)
  Change change;
};

// A run of the record's window: each building enters at its first census in its observed use, is held until its last,
// and is observed at its censuses; the demand of year y (demand[y − from]) holds over [y, y + 1).
struct History {
  int from = 0, to = 0;
  std::vector<std::vector<std::pair<int, std::size_t>>> seen;  // per building: (census year, use)
  std::vector<Timed> changes;                                  // every redraw, those to the use held included
  State last;                                                  // the configuration at the window's end
};

inline History simulate(const Model& M, const std::vector<Life>& lives, const std::vector<double>& demand, int from, int to, std::uint64_t seed) {
  const std::size_t n = lives.size();
  if (M.space && M.space->at.size() != n) throw std::invalid_argument("urban::occupancy::simulate: one life for each building of the space");
  if (demand.size() < static_cast<std::size_t>(to - from + 1)) throw std::invalid_argument("urban::occupancy::simulate: a demand for every year of the window");
  History H;
  H.from = from, H.to = to;
  H.seen.resize(n);
  State s;
  s.use.assign(n, 0), s.alive.assign(n, 0), s.enrolled.assign(n, 0), s.held.assign(n, M.uses->states());
  const auto enter = [&](std::size_t i) {
    s.use[i] = lives[i].initial;
    s.held[i] = M.uses->lump[s.use[i]] == occupied ? M.uses->family[s.use[i]] : M.uses->states();  // unknown if it enters vacant
  };
  for (std::size_t i = 0; i < n; i++) enter(i);
  std::vector<std::size_t> next(n, 0);  // each building's next census
  for (int y = from; y <= to; y++) {
    for (std::size_t i = 0; i < n; i++) {
      if (lives[i].first == y) s.alive[i] = 1, enter(i);
      if (next[i] < lives[i].censuses.size() && lives[i].censuses[next[i]] == y) {
        H.seen[i].push_back({y, s.use[i]});
        next[i]++;
      }
      if (lives[i].last < y) s.alive[i] = 0;
    }
    if (y == to) break;
    for (std::size_t i = 0; i < n; i++)
      if (lives[i].last <= y) s.alive[i] = 0;  // held until its last census
    s.demand = demand[static_cast<std::size_t>(y - from)];
    s.changes.clear();
    nomos::agents::Log log;
    nomos::agents::gillespie(M, s, static_cast<double>(y), static_cast<double>(y + 1), year_seed(seed, y), &log);
    for (std::size_t k = 0; k < s.changes.size(); k++) H.changes.push_back({log.events[k].time, s.changes[k]});
  }
  s.changes.clear();
  H.last = s;
  return H;
}

// ---- the simulated record ---------------------------------------------------------------------------------------------

// A run written as urban-science's record (ADR 0067): each building a coupling, its censuses observations (type
// observation, kind census, a year's label), each change of use an alteration at its instant (kind "use changed",
// "fell vacant" or "let"), all with the provenance "simulated", so that a simulated history is never read as an
// observed one. The window is [from, to + 1).
inline ue::Record to_record(const Uses& S, const Space& X, const std::vector<std::string>& ids, const History& H, std::uint64_t seed) {
  ue::Record R;
  R.window = {ue::mjd(H.from), ue::mjd(H.to + 1)};
  for (std::size_t i = 0; i < ids.size(); i++) R.couplings.push_back({ids[i], "building", ids[i], X.at[i], 0, 0});
  const std::string source = "urban::occupancy (simulated, seed " + std::to_string(seed) + ")";
  for (std::size_t i = 0; i < H.seen.size(); i++)
    for (const auto& [y, u] : H.seen[i]) {
      ue::Event e;
      e.id = "census " + std::to_string(y) + " " + ids[i], e.type = "observation", e.kind = "census";
      e.when = ue::year_span(y), e.couplings = {i}, e.source = source, e.provenance = "simulated";
      (void)u;
      R.events.push_back(e);
    }
  std::size_t k = 0;
  for (const auto& t : H.changes) {
    if (t.change.from == t.change.to) continue;
    ue::Event e;
    const bool fell = !S.is_vacant(t.change.from) && S.is_vacant(t.change.to), let = S.is_vacant(t.change.from) && !S.is_vacant(t.change.to);
    e.id = "change " + std::to_string(k++) + " " + ids[t.change.building];
    e.type = "alteration", e.kind = fell ? "fell vacant" : let ? "let" : "use changed";
    const int y = static_cast<int>(std::floor(t.time));
    const double day = ue::mjd(y) + (t.time - y) * (ue::mjd(y + 1) - ue::mjd(y));
    e.when = {day, day}, e.couplings = {t.change.building}, e.source = source, e.provenance = "simulated";
    R.events.push_back(e);
  }
  return R;
}

// The model's description by ODD+D (Grimm et al. 2020; Müller et al. 2013), kept beside the model.
inline nomos::agents::Description describe() {
  nomos::agents::Description D;
  D.purpose = "How the uses of a district's buildings change and fall vacant, and what interventions on vacant space do (ADR 0071, the design's §10).";
  D.entities = "Buildings (CLUE's base properties), each held by an owner and used by an occupier; its state its use (CLUE's space uses in their families). The space: each building's k nearest buildings.";
  D.scales = "Continuous time in years, censuses at whole years; the district's buildings as they stand at their censuses.";
  D.process = "At rate λ_l an occupier is succeeded at level l, the owner letting the building within its level-l block; at rate ν e^{-β_d d} an occupier leaves with no successor, the building unused or under works; at rate η e^{β_d d} (κ more when enrolled and unused) the owner lets a vacancy, with chance ρ in the family it last held.";
  D.design_concepts = "Emergence: vacancy's persistence and clustering. Interaction: through the owner's choice, by the neighbours' uses. Stochasticity: every draw from keyed streams.";
  D.decisions = "The owner's choice of the next use: a multinomial logit of ln π(y) + β_n A_i(y) (bounded rationality, no planning). The planner's: a weakly coupled POMDP (planning.hpp).";
  D.initialisation = "Each building in its use at its first census.";
  D.input = "The demand index by year; the classification's law π; the affinity of families from a census.";
  D.submodels = "The nested chain (physis decision 0020) is the baseline (ν = η = β_n = 0); the channel of vacancy with its memory of the family held; the incentive; censuses.";
  D.schedule = "Asynchronous in continuous time (Gillespie's direct method), yearly segments with the year's demand: the master equation is exact and the baseline is physis's nested chain.";
  return D;
}

}  // namespace urban::occupancy
