// A run of walkers written as a record of events (urban-science's events/record.hpp; ADR 0067; the design of ADR 0071,
// §0.1 claim 6 and §11): a simulated history, analysed as an observed one and never drawn as one.
//
// Couplings are the network's street spaces ("street:<segment>", typed by the segment's kind) and junction spaces
// ("junction:<junction>"); primitives the cross-sections between a street and each of its two junctions
// ("cross:<segment>:<junction>"), as urban::events::city_network names them. Events, each with provenance "simulated":
//   - regime: a street's regime in force over a span (closed or one-way), touching its street spaces and meeting the
//     cross-sections at its ends;
//   - occurrence "crowd": a street space on which, minute after minute for at least a minute, some walkers (ten, by
//     default) stood each in a density at or above Fruin's level E (1.08 m⁻², about them: walkers.hpp), over the span
//     they did; its peak, the most walkers so;
//   - occurrence "queue": walkers waiting at a segment's mouth, a junction space, for at least a minute, at least so many;
//   - observation "count": what a counter would have read in each hour, touching its street space.
// Time is MJD in UTC (events/record.hpp); the run's second t is mjd0 + t / 86400.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

#include <urban/events/record.hpp>

#include <urban/movement/regime.hpp>
#include <urban/movement/walkers.hpp>

namespace urban::movement {

struct RegimeSpan {
  Street street;
  Regime regime = Regime::open;
  double from = 0, to = 0;   // s of the run
};

struct RecordOptions {
  double mjd0 = 0;                 // the MJD (UTC) of the run's second 0
  double crowd = 10;               // walkers at Fruin's level E, on average over a bin
  double least = 60;               // s: an occurrence lasts at least this
  double queue = 20;               // walkers: a queue counts from this many
  bool counts = true;
  std::string source = "urban-science movement: walkers as agents (simulated)";
};

// The record of a run whose tallies kept densities in bins (Options::minute_bins) and the regimes given.
inline urban::events::Record simulated_record(const Simulation& sim, const std::vector<RegimeSpan>& regimes, const RecordOptions& o) {
  namespace ue = urban::events;
  const Walkway& W = sim.walkway();
  const Tallies& T = sim.tallies();
  ue::Record R;
  const double hours = static_cast<double>(T.hours);
  R.window = {o.mjd0, o.mjd0 + hours / 24};
  const std::size_t E = W.edges(), J = W.junctions();
  for (std::size_t e = 0; e < E; e++) {
    ue::Coupling c;
    c.id = "street:" + std::to_string(e), c.type = W.segments[e].kind, c.name = W.segments[e].name;
    c.at = W.point(2 * e, W.segments[e].length / 2);
    R.couplings.push_back(c);
  }
  for (std::size_t j = 0; j < J; j++) {
    ue::Coupling c;
    c.id = "junction:" + std::to_string(j), c.type = "junction", c.at = W.at[j];
    // Named by the streets that meet there ("Bourke Street and Swanston Street"), as a reader names a corner.
    std::vector<std::string> names;
    for (const std::size_t s : W.leaving[j]) {
      const auto& n = W.segments[s / 2].name;
      if (!n.empty() && std::find(names.begin(), names.end(), n) == names.end()) names.push_back(n);
    }
    std::sort(names.begin(), names.end());
    for (std::size_t k = 0; k < names.size(); k++) c.name += (k == 0 ? "" : (k + 1 == names.size() ? " and " : ", ")) + names[k];
    R.couplings.push_back(c);
  }
  std::vector<std::array<std::size_t, 2>> cross(E);
  for (std::size_t e = 0; e < E; e++) {
    std::size_t k = 0;
    for (const std::size_t j : {W.segments[e].u, W.segments[e].v}) {
      ue::Primitive p;
      p.id = "cross:" + std::to_string(e) + ":" + std::to_string(j), p.type = "cross-section";
      p.a = e, p.b = E + j, p.area = W.segments[e].width * 3.0;   // the walkable width times a storey of 3 m, for drawing
      R.primitives.push_back(p);
      cross[e][k++] = R.primitives.size() - 1;
    }
  }
  const auto mjd = [&](double t) { return o.mjd0 + t / 86400; };
  std::size_t serial = 0;
  const auto id = [&](const char* kind) { return std::string("sim-") + kind + "-" + std::to_string(serial++); };
  for (const auto& r : regimes) {
    if (r.regime == Regime::open) continue;
    ue::Event ev;
    ev.id = id("regime"), ev.type = "regime", ev.kind = regime_name(r.regime);
    ev.when = {mjd(r.from), mjd(std::max(r.from, r.to))};
    for (const std::size_t e : r.street.edges) ev.couplings.push_back(e), ev.primitives.push_back(cross[e][0]), ev.primitives.push_back(cross[e][1]);
    ev.source = o.source, ev.provenance = "simulated";
    R.events.push_back(ev);
  }
  // Crowds, from the density bins.
  if (T.bin > 0 && !T.minutes_e.empty()) {
    const std::size_t bins = T.minutes_e.size() / E;
    for (std::size_t e = 0; e < E; e++) {
      std::size_t b = 0;
      while (b < bins) {
        const auto walkers = [&](std::size_t q) { return static_cast<double>(T.minutes_e[q * E + e]) / T.bin; };
        if (walkers(b) < o.crowd) { b++; continue; }
        std::size_t c = b;
        double peak = 0;
        while (c < bins && walkers(c) >= o.crowd) peak = std::max(peak, walkers(c)), c++;
        if (static_cast<double>(c - b) * T.bin >= o.least) {
          ue::Event ev;
          ev.id = id("crowd"), ev.type = "occurrence", ev.kind = "crowd";
          ev.when = {mjd(static_cast<double>(b) * T.bin), mjd(static_cast<double>(c) * T.bin)};
          ev.couplings = {e};
          char buf[96];
          std::snprintf(buf, sizeof buf, "%s; at most %.0f walkers at Fruin's level E", o.source.c_str(), peak);
          ev.source = buf, ev.provenance = "simulated";
          R.events.push_back(ev);
        }
        b = c;
      }
    }
  }
  // Queues, from the hours' longest queue at each mouth.
  if (!T.queue_max.empty()) {
    const std::size_t S = W.states();
    for (std::size_t h = 0; h < T.hours; h++)
      for (std::size_t s = 0; s < S; s++) {
        const double q = T.queue_max[h * S + s];
        if (q < o.queue || T.queue[h * S + s] < o.least * o.queue) continue;
        ue::Event ev;
        ev.id = id("queue"), ev.type = "occurrence", ev.kind = "queue";
        ev.when = {mjd(3600.0 * static_cast<double>(h)), mjd(3600.0 * static_cast<double>(h + 1))};
        ev.couplings = {E + W.tail(s)};
        ev.primitives = {cross[s / 2][W.tail(s) == W.segments[s / 2].u ? 0 : 1]};
        char buf[128];
        std::snprintf(buf, sizeof buf, "%s; at most %.0f waiting, %.0f walker-minutes", o.source.c_str(), q, T.queue[h * S + s] / 60);
        ev.source = buf, ev.provenance = "simulated";
        R.events.push_back(ev);
      }
  }
  if (o.counts)
    for (std::size_t k = 0; k < sim.counters().size(); k++)
      for (std::size_t h = 0; h < T.hours; h++) {
        const auto r = sim.reading(k, h);
        ue::Event ev;
        ev.id = id("count"), ev.type = "observation", ev.kind = "count";
        ev.when = {mjd(3600.0 * static_cast<double>(h)), mjd(3600.0 * static_cast<double>(h + 1))};
        ev.couplings = {sim.counters()[k].edge};
        char buf[160];
        std::snprintf(buf, sizeof buf, "%s; %s reads %.1f and %.1f each way", o.source.c_str(), sim.counters()[k].name.c_str(), r[0], r[1]);
        ev.source = buf, ev.provenance = "simulated";
        R.events.push_back(ev);
      }
  return R;
}

}  // namespace urban::movement
