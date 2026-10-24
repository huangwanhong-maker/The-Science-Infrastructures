// The spacetime complex of a record of transactions (theory note §5).
//
// Agents hold value through time and exchange it at events. The record's encounters (the events
// in time, the agents present, each pair's exchanges) make morphe's worldline complex and its
// refinement Ŵ, whose vertices are the presences (x, v), whose edges are the agents' holdings
// between consecutive events and the exchanges, and whose 2-cells are the lenses (two agents
// between consecutive exchanges) and the triangles of exchanges at one event. Every edge carries a
// transport:
//
//   holding by x from event p to event q:  U = exp(r_x (t_q − t_p)), r_x the agent's rate of holding
//     per unit of time (interest, inflation, depreciation, storage; negative where value decays),
//     one for each agent or one for each agent at each event, held until its next;
//   an exchange from a to b at an event:  U = the value registered by b per unit surrendered by a
//     (for currencies, the rate of exchange).
//
// So Ŵ is a transfer network (economics::transfer) on the presences, with faces. Its curvature,
// the deficit round each 2-cell, has a closed form: round the lens of a pair exchanging a → b at p
// and at q,
//
//   κ = (r_a − r_b)(t_q − t_p) + log U_q − log U_p,
//
// the excess return of holding with a and returning through b, which for currencies is the
// departure from uncovered interest parity over the interval; round a triangle of exchanges at an
// event, the triangular arbitrage of its three rates. The deficit round any loop is the sum of the
// curvatures it encloses (Stokes). Flatness is an intertemporal price system: a potential φ(x, v),
// the log price of x's value at v in a common numeraire, with every transport its difference;
// on a finite record, the absence of arbitrage across agents and dates.
//
// Reach: value surrendered at a presence travels forwards along its holder's worldline and across
// exchanges in their direction. With each exchange's viability yes, no or unknown, a route's
// viability is the least of its steps and a presence's the greatest over the routes reaching it
// (Kleene's strong three-valued logic, no < unknown < yes): the user's regenerative reach (Huang
// 2026, Toward a Pathology of Value Circulation, §4), with its two choices fixed. Exchanges at
// one event may be chained.

#pragma once

#include <economics/transfer/network.hpp>
#include <morphe/topology/worldlines.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace economics {

struct Transactions {
  std::vector<std::string> agents;
  std::vector<double> holding_rate;  // r per unit of time for each agent; empty for 0
  std::vector<std::vector<double>> holding_rates;  // or r for each agent at each event, from that event to its next; empty for holding_rate
  struct Event {
    double time = 0;
    std::vector<std::size_t> present;
    std::string label;
  };
  std::vector<Event> events;  // in time
  struct Exchange {
    std::size_t event = 0, from = 0, to = 0;
    double ratio = 1;     // value registered by `to` per unit surrendered by `from`
    double quantity = 0;  // units surrendered (optional)
  };
  std::vector<Exchange> exchanges;
};

struct Spacetime {
  morphe::topology::Encounters encounters;
  morphe::topology::WorldlineComplex worldlines;
  morphe::topology::RefinedWorldlines refined;
  std::vector<std::vector<std::size_t>> exchange_of;  // pair -> meeting index -> exchange
  TransferNetwork network;                            // on the presences: holdings, exchanges, faces
  Cochain curvature;                                  // the deficit round each 2-cell, negated: d₁ log U
  std::size_t lenses = 0;                             // 2-cells [0, lenses) are lenses, the rest triangles
};

inline Spacetime spacetime(const Transactions& T) {
  const std::size_t A = T.agents.size(), V = T.events.size();
  for (std::size_t v = 1; v < V; v++)
    if (T.events[v].time < T.events[v - 1].time) throw std::invalid_argument("economics::spacetime: events in time");
  Spacetime S;
  S.encounters.participants = A;
  for (const auto& e : T.events) {
    auto here = e.present;
    S.encounters.present.push_back(here);
  }
  auto present = [&](std::size_t x, std::size_t v) {
    return std::find(T.events.at(v).present.begin(), T.events[v].present.end(), x) != T.events[v].present.end();
  };
  // Pairs by unordered agents, each with its exchanges in time, at most one at an event.
  std::map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> by_pair;
  for (std::size_t k = 0; k < T.exchanges.size(); k++) {
    const auto& x = T.exchanges[k];
    if (x.from == x.to || !present(x.from, x.event) || !present(x.to, x.event))
      throw std::invalid_argument("economics::spacetime: an exchange between two agents present at its event");
    if (!(x.ratio > 0)) throw std::invalid_argument("economics::spacetime: an exchange's ratio must be positive");
    by_pair[std::minmax(x.from, x.to)].push_back(k);
  }
  for (auto& [ab, list] : by_pair) {
    std::sort(list.begin(), list.end(), [&](std::size_t i, std::size_t j) { return T.exchanges[i].event < T.exchanges[j].event; });
    morphe::topology::Encounters::Pair p;
    p.a = ab.first;
    p.b = ab.second;
    for (std::size_t k : list) {
      if (!p.meetings.empty() && p.meetings.back() == T.exchanges[k].event)
        throw std::invalid_argument("economics::spacetime: two exchanges of one pair at one event");
      p.meetings.push_back(T.exchanges[k].event);
    }
    S.encounters.pairs.push_back(std::move(p));
    S.exchange_of.push_back(list);
  }
  S.worldlines = morphe::topology::worldline_complex(S.encounters);
  S.refined = morphe::topology::refined_worldlines(S.encounters, S.worldlines);
  const auto& R = S.refined;
  // The network on the presences.
  auto& N = S.network;
  for (const auto& [x, v] : R.presences) N.positions.push_back(T.agents[x] + "@" + (T.events[v].label.empty() ? std::to_string(v) : T.events[v].label));
  for (std::size_t e = 0; e < R.temporal; e++) {
    const auto& s = S.worldlines.segments[e];
    const double r = !T.holding_rates.empty() ? T.holding_rates.at(s.participant).at(s.from) : T.holding_rate.empty() ? 0.0 : T.holding_rate.at(s.participant);
    N.legs.push_back({R.cw.edges[e][0], R.cw.edges[e][1], std::exp(r * (T.events[s.to].time - T.events[s.from].time)), 0.0, "holding"});
  }
  for (std::size_t e = R.temporal; e < R.cw.edges.size(); e++) {
    const auto [t, k] = R.rung_of_edge[e - R.temporal];
    const auto& x = T.exchanges[S.exchange_of[t][k]];
    const bool forwards = x.from == S.encounters.pairs[t].a;  // the rung runs from (a, v) to (b, v)
    N.legs.push_back({R.cw.edges[e][0], R.cw.edges[e][1], forwards ? x.ratio : 1.0 / x.ratio, x.quantity, "exchange"});
  }
  for (const auto& f : R.cw.faces) N.faces.push_back(f.boundary);
  S.lenses = S.worldlines.lenses.size();
  S.curvature = physis::discretisation::d1(R.cw, N.transport());
  return S;
}

// Kleene's strong three-valued logic.
enum class Viability : int { no = 0, unknown = 1, yes = 2 };

// The viability of reaching every presence with value surrendered at (agent, event): forwards
// along worldlines (always viable) and across exchanges in their direction, each with the
// viability given (`viability[k]` for exchange k; empty for all yes).
inline std::vector<Viability> reach(const Spacetime& S, const Transactions& T, std::size_t agent, std::size_t event, const std::vector<Viability>& viability = {}) {
  const auto& R = S.refined;
  std::vector<Viability> best(R.presences.size(), Viability::no);
  const auto src = R.presence_of.find({agent, event});
  if (src == R.presence_of.end()) throw std::invalid_argument("economics::reach: the agent is not present at the event");
  best[src->second] = Viability::yes;
  // Steps out of each presence: (target, viability).
  std::vector<std::vector<std::pair<std::size_t, Viability>>> out(R.presences.size());
  for (std::size_t e = 0; e < R.temporal; e++) out[R.cw.edges[e][0]].push_back({R.cw.edges[e][1], Viability::yes});
  for (std::size_t e = R.temporal; e < R.cw.edges.size(); e++) {
    const auto [t, k] = R.rung_of_edge[e - R.temporal];
    const std::size_t x = S.exchange_of[t][k];
    const auto& ex = T.exchanges[x];
    const Viability w = viability.empty() ? Viability::yes : viability.at(x);
    const std::size_t from = R.presence_of.at({ex.from, ex.event}), to = R.presence_of.at({ex.to, ex.event});
    out[from].push_back({to, w});
  }
  // Presences in the order of their events; within an event, relax to a fixed point (exchanges
  // may chain), then carry forwards.
  std::vector<std::vector<std::size_t>> at(T.events.size());
  for (std::size_t p = 0; p < R.presences.size(); p++) at[R.presences[p][1]].push_back(p);
  auto relax = [&](std::size_t p) {
    bool changed = false;
    for (const auto& [q, w] : out[p]) {
      const Viability through = static_cast<Viability>(std::min(static_cast<int>(best[p]), static_cast<int>(w)));
      if (static_cast<int>(through) > static_cast<int>(best[q])) { best[q] = through; changed = true; }
    }
    return changed;
  };
  for (std::size_t v = 0; v < T.events.size(); v++) {
    for (bool changed = true; changed;) {
      changed = false;
      for (std::size_t p : at[v])
        for (const auto& [q, w] : out[p])
          if (R.presences[q][1] == v) {
            const Viability through = static_cast<Viability>(std::min(static_cast<int>(best[p]), static_cast<int>(w)));
            if (static_cast<int>(through) > static_cast<int>(best[q])) { best[q] = through; changed = true; }
          }
    }
    for (std::size_t p : at[v]) relax(p);
  }
  return best;
}

}  // namespace economics
