// Lanes carried through: the study's own provisional edit of the walking network, not a domain's.
//
// The urban-science domain opens a passage through a block only where the block is a 2-cell of its complex
// (urban::landscape::open_passage, on plan_passage). A block holding a dead-end lane is not a disc bounded by a circuit,
// so the domain leaves it as a hole (district.hpp in the domain: "a face with a dead end ... stays a hole"); in the
// blocks this study redesigns, the five largest are such holes, and the domain can open no passage through them. The
// intervention Melbourne's lanes invite there is to carry a blind lane through to the street beyond it. The study writes
// that edit here, on the network the domain is given (not on its complex), and lets the domain build the district from
// the edited network: the faces are traced again, a hole that the new lane closes into circuits becomes 2-cells, and
// the complex is checked regular by the domain as always. It is reported as a gap (README §8): an edit "carry a dead
// end through" belongs beside the domain's cut, passage and join.
//
// The edit. A dead end is a junction of degree one at the end of a laneway. Its lane is carried on in the direction of
// the lane's last piece, as a straight segment, to the first segment of the network the ray meets (not the lane's
// own), within `reach` metres; the segment met is split at that point by a new junction (or the lane meets an existing
// junction there, within a metre). The new segment is a laneway.

#pragma once

#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include <urban/landscape/district.hpp>

#include "sunlight.hpp"

namespace study {

struct Through {
  std::size_t junction = 0;  // the dead end, in the network's numbering (kept by every edit: junctions are only added)
  Point from{}, direction{}, to{};
  double length = 0;
  std::string lane, street;  // the lane carried, and the segment it meets
};

namespace detail {

// The first point at which the ray from p along u meets a segment of N other than `own`: (segment, piece, t, point).
struct Hit { std::size_t segment = 0, piece = 0; double t = 0; Point at{}; };
inline std::optional<Hit> first_hit(const urban::landscape::Network& N, std::size_t own, const Point& p, const Point& u, double reach) {
  std::optional<Hit> best;
  for (std::size_t s = 0; s < N.segments.size(); s++) {
    if (s == own) continue;
    const auto& g = N.segments[s];
    const auto line = g.line.empty() ? std::vector<Point>{N.junctions[g.u], N.junctions[g.v]} : g.line;
    for (std::size_t k = 0; k + 1 < line.size(); k++) {
      const Point& a = line[k];
      const Point& b = line[k + 1];
      const double ex = b[0] - a[0], ey = b[1] - a[1];
      const double den = u[0] * ey - u[1] * ex;
      if (std::fabs(den) < 1e-12) continue;
      const double wx = a[0] - p[0], wy = a[1] - p[1];
      const double t = (wx * ey - wy * ex) / den, w = (wx * u[1] - wy * u[0]) / den;
      if (t > 0.5 && t <= reach && w >= 0 && w <= 1 && (!best || t < best->t)) best = Hit{s, k, t, {p[0] + t * u[0], p[1] + t * u[1]}};
    }
  }
  return best;
}

inline std::size_t segment_of(const urban::landscape::Network& N, std::size_t j) {
  for (std::size_t s = 0; s < N.segments.size(); s++)
    if (N.segments[s].u == j || N.segments[s].v == j) return s;
  return N.segments.size();
}

}  // namespace detail

// The candidates: every dead end of a laneway at least `least` metres long whose end lies inside `area`, carried on to
// the first segment it meets within `reach`, the new lane's middle inside `area`.
inline std::vector<Through> lanes_through(const urban::landscape::Network& N, const Ring& area, double least = 15, double reach = 120) {
  std::vector<std::size_t> degree(N.junctions.size(), 0);
  for (const auto& g : N.segments) degree[g.u]++, degree[g.v]++;
  std::vector<Through> out;
  for (std::size_t j = 0; j < N.junctions.size(); j++) {
    if (degree[j] != 1 || !inside(area, N.junctions[j])) continue;
    const std::size_t s = detail::segment_of(N, j);
    const auto& g = N.segments[s];
    if (g.kind != "laneway") continue;
    auto line = g.line.empty() ? std::vector<Point>{N.junctions[g.u], N.junctions[g.v]} : g.line;
    if (g.v != j) std::reverse(line.begin(), line.end());  // walked toward the dead end
    double L = 0;
    for (std::size_t k = 0; k + 1 < line.size(); k++) L += std::hypot(line[k + 1][0] - line[k][0], line[k + 1][1] - line[k][1]);
    if (L < least) continue;
    const Point& a = line[line.size() - 2];
    const Point& b = line.back();
    const double d = std::hypot(b[0] - a[0], b[1] - a[1]);
    const Point u{(b[0] - a[0]) / d, (b[1] - a[1]) / d};
    const auto hit = detail::first_hit(N, s, b, u, reach);
    if (!hit) continue;
    const Point mid{(b[0] + hit->at[0]) / 2, (b[1] + hit->at[1]) / 2};
    if (!inside(area, mid)) continue;
    out.push_back({j, b, u, hit->at, hit->t, g.name.empty() ? "a lane" : g.name, N.segments[hit->segment].name.empty() ? N.segments[hit->segment].kind : N.segments[hit->segment].name});
  }
  return out;
}

// The network with the chosen lanes carried through, in the order given; each ray cast again on the network as it
// stands, so that a lane may meet one carried before it.
inline urban::landscape::Network carry_through(urban::landscape::Network N, const std::vector<Through>& chosen, double reach = 120) {
  for (const auto& c : chosen) {
    const std::size_t own = detail::segment_of(N, c.junction);
    const auto hit = detail::first_hit(N, own, c.from, c.direction, reach);
    if (!hit) throw std::invalid_argument("study::carry_through: " + c.lane + " meets nothing");
    auto& g = N.segments[hit->segment];
    auto line = g.line.empty() ? std::vector<Point>{N.junctions[g.u], N.junctions[g.v]} : g.line;
    std::size_t meet;
    const auto near = [&](const Point& q) { return std::hypot(q[0] - hit->at[0], q[1] - hit->at[1]) <= 1.0; };
    if (near(N.junctions[g.u])) meet = g.u;
    else if (near(N.junctions[g.v])) meet = g.v;
    else {
      meet = N.junctions.size();
      N.junctions.push_back(hit->at);
      urban::landscape::Segment first = g, second = g;
      first.line.assign(line.begin(), line.begin() + static_cast<std::ptrdiff_t>(hit->piece) + 1);
      first.line.push_back(hit->at), first.v = meet;
      second.line = {hit->at};
      second.line.insert(second.line.end(), line.begin() + static_cast<std::ptrdiff_t>(hit->piece) + 1, line.end());
      second.u = meet;
      g = first;
      N.segments.push_back(second);
    }
    urban::landscape::Segment lane;
    lane.u = c.junction, lane.v = meet, lane.line = {c.from, N.junctions[meet]}, lane.kind = "laneway", lane.name = c.lane + " carried through";
    N.segments.push_back(lane);
  }
  return N;
}

}  // namespace study
