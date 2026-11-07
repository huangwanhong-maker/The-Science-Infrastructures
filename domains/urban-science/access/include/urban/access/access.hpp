// How far one can go: the ways of a place as a network, and the time to each of their points from a place picked, on foot or by road
// (the method of GOSTnets and of accessibility analysis generally; Hansen 1959; written here from that statement).
//
// A way is a polyline in the local frame (m). Points of different ways are one junction where they coincide to `snap` metres (OSM and
// Overture share vertices where ways join; a bridge over a road shares none and stays apart). Each piece is an arc both ways, its
// time its length over the speed of its way; a way can be closed to a mode by a speed of zero. The start joins the nearest ways at
// walking speed.
//
// Speeds are declared and are the caller's: walking a level 1.34 m/s (Weidmann 1993) by default, a road's by its class. Junctions, signals,
// congestion and gradients are not in it, so that a time is a lower bound on the journey.

#pragma once

#include <morphe/combinatorics/shortest_paths.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace urban::access {

namespace mc = morphe::combinatorics;

struct Way {
  std::vector<double> xy;   // x0, y0, x1, y1, ...
  double speed = 1.34;      // m/s; 0 closes it
};

struct Network {
  std::vector<double> x, y;                 // the junctions and points
  mc::WeightedAdjacency arcs;               // their times, s
  std::vector<std::vector<std::size_t>> points;   // way → its points
};

inline Network network_of(const std::vector<Way>& ways, double snap = 1.0) {
  if (!(snap > 0)) throw std::invalid_argument("access: a snapping distance is positive");
  Network n;
  n.points.resize(ways.size());
  std::unordered_map<std::uint64_t, std::size_t> at;
  const auto key = [&](double x, double y) {
    const auto a = static_cast<std::int64_t>(std::llround(x / snap)), b = static_cast<std::int64_t>(std::llround(y / snap));
    return (static_cast<std::uint64_t>(a + (1LL << 30)) << 32) ^ static_cast<std::uint64_t>(b + (1LL << 30));
  };
  const auto point = [&](double x, double y) {
    const auto k = key(x, y);
    const auto it = at.find(k);
    if (it != at.end()) return it->second;
    const std::size_t id = n.x.size();
    at.emplace(k, id);
    n.x.push_back(x);
    n.y.push_back(y);
    n.arcs.emplace_back();
    return id;
  };
  for (std::size_t w = 0; w < ways.size(); w++) {
    const Way& way = ways[w];
    if (way.xy.size() % 2) throw std::invalid_argument("access: a way's points come in pairs");
    std::size_t last = 0;
    for (std::size_t i = 0; i + 1 < way.xy.size(); i += 2) {
      const std::size_t p = point(way.xy[i], way.xy[i + 1]);
      n.points[w].push_back(p);
      if (i > 0 && way.speed > 0 && p != last) {
        const double t = std::hypot(n.x[p] - n.x[last], n.y[p] - n.y[last]) / way.speed;
        n.arcs[last].push_back({p, t});
        n.arcs[p].push_back({last, t});
      }
      last = p;
    }
  }
  return n;
}

// Joins a way's loose end to the nearest other point within `radius` metres, walking at a speed: a tile cuts a way at its edge without a
// shared point, and the neighbour's piece begins a metre or two away. Returns how many ends were joined. A dead end that happens to
// stop that near another way is joined too, so that the radius is kept small.
inline std::size_t stitch_ends(Network& n, double radius, double speed) {
  if (!(radius > 0)) return 0;
  std::unordered_map<std::uint64_t, std::vector<std::size_t>> cells;
  const auto cell = [&](double x, double y, int di = 0, int dj = 0) {
    const auto a = static_cast<std::int64_t>(std::floor(x / radius)) + di, b = static_cast<std::int64_t>(std::floor(y / radius)) + dj;
    return (static_cast<std::uint64_t>(a + (1LL << 30)) << 32) ^ static_cast<std::uint64_t>(b + (1LL << 30));
  };
  for (std::size_t p = 0; p < n.x.size(); p++) if (!n.arcs[p].empty()) cells[cell(n.x[p], n.y[p])].push_back(p);
  std::size_t joined = 0;
  std::vector<std::pair<std::size_t, std::size_t>> pairs;
  std::vector<double> times;
  for (std::size_t p = 0; p < n.x.size(); p++) {
    if (n.arcs[p].size() != 1) continue;     // one neighbour: a loose end
    const std::size_t self = n.arcs[p][0].to;
    double best = radius;
    std::size_t at = p;
    for (int di = -1; di <= 1; di++) for (int dj = -1; dj <= 1; dj++) {
      const auto it = cells.find(cell(n.x[p], n.y[p], di, dj));
      if (it == cells.end()) continue;
      for (const auto q : it->second) {
        if (q == p || q == self) continue;
        const double d = std::hypot(n.x[q] - n.x[p], n.y[q] - n.y[p]);
        if (d < best) { best = d; at = q; }
      }
    }
    if (at != p) { pairs.push_back({p, at}); times.push_back(best / speed); joined++; }
  }
  for (std::size_t k = 0; k < pairs.size(); k++) {
    n.arcs[pairs[k].first].push_back({pairs[k].second, times[k]});
    n.arcs[pairs[k].second].push_back({pairs[k].first, times[k]});
  }
  return joined;
}

struct Reach {
  std::vector<double> time;       // s to each point; +inf where unreached
  std::size_t start = 0;          // the point it entered at
  double start_gap = 0;           // the metres walked to it
};

// The times from a place, or from several, entering at the nearest point of the network that is joined to anything.
inline Reach reach(const Network& n, const std::vector<std::pair<double, double>>& from, double walk_speed = 1.34, double limit = std::numeric_limits<double>::infinity()) {
  if (n.x.empty() || from.empty()) throw std::invalid_argument("access: there is a network and a place");
  std::vector<std::pair<std::size_t, double>> sources;
  Reach out;
  for (const auto& [x, y] : from) {
    double best = std::numeric_limits<double>::infinity();
    std::size_t at = 0;
    for (std::size_t p = 0; p < n.x.size(); p++) {
      if (n.arcs[p].empty()) continue;
      const double d = std::hypot(n.x[p] - x, n.y[p] - y);
      if (d < best) { best = d; at = p; }
    }
    if (!std::isfinite(best)) throw std::invalid_argument("access: the network joins nothing");
    sources.push_back({at, best / walk_speed});
    out.start = at, out.start_gap = best;
  }
  out.time = mc::dijkstra(n.arcs, sources, limit).distance;
  return out;
}

// The metres of way within each time, cumulative: each piece counted once and in proportion when a limit falls along it.
inline std::vector<double> length_within(const Network& n, const std::vector<double>& time, const std::vector<double>& limits) {
  std::vector<double> out(limits.size(), 0.0);
  for (std::size_t p = 0; p < n.arcs.size(); p++) for (const auto& a : n.arcs[p]) {
    if (a.to < p) continue;
    const double L = std::hypot(n.x[p] - n.x[a.to], n.y[p] - n.y[a.to]);
    const double t0 = std::min(time[p], time[a.to]), t1 = std::max(time[p], time[a.to]);
    for (std::size_t k = 0; k < limits.size(); k++) {
      if (!std::isfinite(t0) || t0 >= limits[k]) continue;
      out[k] += t1 <= limits[k] || !std::isfinite(t1) ? (std::isfinite(t1) ? L : L * std::min(1.0, (limits[k] - t0) / std::max(a.weight, 1e-9))) : L * (limits[k] - t0) / std::max(t1 - t0, 1e-9);
    }
  }
  return out;
}

}  // namespace urban::access
