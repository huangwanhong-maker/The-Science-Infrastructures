// The walking network of the walkers (theory/movement.md §2): junctions and segments, each segment walked in two
// directions, the directed segments being the states of a walker's route.
//
// The network is the one the domain and the workbench build (urban::landscape::Network: the city model's ways where
// people walk, cut to a box and split at their shared points; applications/urban-science-workbench/src/district.js).
// To it the walkers add what a route and a step depend on:
//   - a length (the polyline's) and a walkable width, by class of street after Hoddle's grid of 1837, whose streets
//     are 99 ft (30.2 m) wide and whose little streets 33 ft (10.1 m) (`Widths`: assumptions, the same table as the
//     district-economy flow study, research/experiments/2026-11-05-a-district-economy-as-a-flow, README §2.1, so that
//     its logit flow and these walkers stand on one network);
//   - the parent length: the length of the segment as built of which a segment is a piece (itself, unless a caller
//     cuts streets at new junctions), so that a cut segment keeps its weight in the route law;
//   - a rise from tail to head (m), from data on footpaths' steepness where there are any, else zero;
//   - for each junction, the carriageways that meet there (segments of the kind "street"), whether a signal controls
//     its crossings, and its axis (the direction, mod π, of its first carriageway), which splits its crossings into
//     two phases.
//
// Directed segments: state 2e runs from segment e's tail u to its head v, state 2e + 1 from v to u. A crossing at a
// junction x, from state a into state b, counts the carriageways the walker crosses, in expectation: one when x is a
// road intersection (three carriageways or more) and the walker goes straight on (its heading turns by less than 45°),
// one half when it turns there (it crosses one of the two streets, on its side of the street, half the time), none at
// a junction of fewer carriageways or on turning back. Every quantity is in metres and seconds.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <urban/landscape/district.hpp>

namespace urban::movement {

using Point = std::array<double, 2>;
inline constexpr std::size_t none = static_cast<std::size_t>(-1);

// Walkable widths by class of street (assumptions; the flow study's table, its README §2.1).
struct Widths {
  double swanston = 12;     // m: Swanston Street, closed to private cars since 1992, its footpaths widened
  double main = 9;          // the other streets of Hoddle's 99 ft reserve: two footpaths of 4.5 m
  double little = 4;        // the little streets of 33 ft: two footpaths of 2 m
  double street = 4;        // other streets
  double bourke_mall = 12;  // the Bourke Street Mall beside its tramway
  double mall = 5;          // other pedestrian ways and arcades
  double laneway = 3.5;     // the city model's default width of a service way
};

inline bool main_street(const std::string& n) {
  static const std::set<std::string> s{"Collins Street", "Bourke Street", "Elizabeth Street", "Russell Street", "Lonsdale Street", "Flinders Street",
                                       "La Trobe Street", "Queen Street", "William Street", "Exhibition Street", "Spring Street", "King Street"};
  return s.count(n) > 0;
}
inline bool little_street(const std::string& n) {
  static const std::set<std::string> s{"Flinders Lane", "Little Collins Street", "Little Bourke Street", "Little Lonsdale Street"};
  return s.count(n) > 0;
}
inline double width_of(const Widths& W, const std::string& kind, const std::string& name) {
  if (kind == "laneway") return W.laneway;
  if (kind == "mall") return name.rfind("Bourke Street", 0) == 0 ? W.bourke_mall : W.mall;
  if (name == "Swanston Street") return W.swanston;
  if (main_street(name)) return W.main;
  if (little_street(name)) return W.little;
  return W.street;
}

struct Segment {
  std::size_t u = 0, v = 0;   // tail, head
  std::vector<Point> line;    // from u to v
  double length = 0, width = 0, parent = 0, rise = 0;   // m
  std::string kind = "street", name;
};

struct Walkway {
  std::vector<Point> at;                              // junctions, local east and north (m)
  std::vector<Segment> segments;
  std::vector<char> signalled;                        // per junction
  // Derived by index():
  std::vector<std::vector<std::size_t>> leaving, entering;   // per junction, its directed segments out and in
  std::vector<std::size_t> carriageways;              // per junction
  std::vector<double> axis;                           // per junction, radians mod π

  std::size_t junctions() const { return at.size(); }
  std::size_t edges() const { return segments.size(); }
  std::size_t states() const { return 2 * segments.size(); }
  std::size_t tail(std::size_t s) const { return s % 2 ? segments[s / 2].v : segments[s / 2].u; }
  std::size_t head(std::size_t s) const { return s % 2 ? segments[s / 2].u : segments[s / 2].v; }
  static std::size_t reverse(std::size_t s) { return s ^ 1u; }
  double length(std::size_t s) const { return segments[s / 2].length; }
  double width(std::size_t s) const { return segments[s / 2].width; }
  // The grade along the state, rise over run (positive uphill).
  double grade(std::size_t s) const {
    const auto& e = segments[s / 2];
    const double g = e.length > 0 ? e.rise / e.length : 0;
    return s % 2 ? -g : g;
  }
  // The position at arc length `along` from the state's tail.
  Point point(std::size_t s, double along) const {
    const auto& e = segments[s / 2];
    const double d = s % 2 ? e.length - along : along;
    double left = std::clamp(d, 0.0, e.length);
    for (std::size_t k = 0; k + 1 < e.line.size(); k++) {
      const double piece = std::hypot(e.line[k + 1][0] - e.line[k][0], e.line[k + 1][1] - e.line[k][1]);
      if (left <= piece || k + 2 == e.line.size()) {
        const double t = piece > 0 ? std::clamp(left / piece, 0.0, 1.0) : 0.0;
        return {e.line[k][0] + t * (e.line[k + 1][0] - e.line[k][0]), e.line[k][1] + t * (e.line[k + 1][1] - e.line[k][1])};
      }
      left -= piece;
    }
    return e.line.back();
  }
  // Unit headings: leaving the tail, and arriving at the head.
  Point heading_out(std::size_t s) const {
    const auto& l = segments[s / 2].line;
    const Point a = s % 2 ? l[l.size() - 1] : l[0], b = s % 2 ? l[l.size() - 2] : l[1];
    const double n = std::hypot(b[0] - a[0], b[1] - a[1]);
    return n > 0 ? Point{(b[0] - a[0]) / n, (b[1] - a[1]) / n} : Point{1, 0};
  }
  Point heading_in(std::size_t s) const {
    const auto& l = segments[s / 2].line;
    const Point a = s % 2 ? l[1] : l[l.size() - 2], b = s % 2 ? l[0] : l[l.size() - 1];
    const double n = std::hypot(b[0] - a[0], b[1] - a[1]);
    return n > 0 ? Point{(b[0] - a[0]) / n, (b[1] - a[1]) / n} : Point{1, 0};
  }
  bool carriageway(std::size_t e) const { return segments[e].kind == "street"; }

  // Carriageways crossed, in expectation, from state a into state b at their junction (the header).
  double crossings(std::size_t a, std::size_t b) const {
    const std::size_t x = head(a);
    if (b == reverse(a) || carriageways[x] < 3) return 0;
    const Point i = heading_in(a), o = heading_out(b);
    return i[0] * o[0] + i[1] * o[1] > std::cos(std::numbers::pi / 4) ? 1.0 : 0.5;
  }
  // The signal phase that lets a walker leave junction tail(b) along b: 0 along the junction's axis, 1 across it.
  int phase(std::size_t b) const {
    const std::size_t x = tail(b);
    const Point o = heading_out(b);
    const double c = std::fabs(o[0] * std::cos(axis[x]) + o[1] * std::sin(axis[x]));
    return c >= std::cos(std::numbers::pi / 4) ? 0 : 1;
  }

  void index() {
    const std::size_t n = at.size();
    if (signalled.size() != n) signalled.assign(n, 0);
    leaving.assign(n, {}), entering.assign(n, {}), carriageways.assign(n, 0), axis.assign(n, std::nan(""));
    for (std::size_t s = 0; s < states(); s++) leaving[tail(s)].push_back(s), entering[head(s)].push_back(s);
    for (std::size_t e = 0; e < segments.size(); e++) {
      if (!carriageway(e)) continue;
      for (const std::size_t x : {segments[e].u, segments[e].v}) {
        carriageways[x]++;
        if (std::isnan(axis[x])) {
          const Point h = heading_out(x == segments[e].u ? 2 * e : 2 * e + 1);
          double a = std::atan2(h[1], h[0]);
          if (a < 0) a += std::numbers::pi;
          axis[x] = a;
        }
      }
    }
    for (auto& a : axis)
      if (std::isnan(a)) a = 0;
  }
};

// The walkers' network of a district's walking network, with widths by class (and, where given, rises per segment,
// signalled junctions). Segments of no length are refused.
inline Walkway walkway(const urban::landscape::Network& N, const Widths& W = {}) {
  Walkway out;
  for (const auto& p : N.junctions) out.at.push_back({p[0], p[1]});
  for (const auto& s : N.segments) {
    Segment e;
    e.u = s.u, e.v = s.v, e.kind = s.kind, e.name = s.name;
    if (s.line.empty()) e.line = {out.at.at(s.u), out.at.at(s.v)};
    else
      for (const auto& p : s.line) e.line.push_back({p[0], p[1]});
    for (std::size_t k = 0; k + 1 < e.line.size(); k++) e.length += std::hypot(e.line[k + 1][0] - e.line[k][0], e.line[k + 1][1] - e.line[k][1]);
    if (!(e.length > 0)) throw std::invalid_argument("urban::movement::walkway: a segment of no length (" + s.name + ")");
    e.parent = e.length;
    e.width = width_of(W, e.kind, e.name);
    out.segments.push_back(std::move(e));
  }
  out.index();
  return out;
}

// The junctions of the largest connected component (by junctions), as a mask.
inline std::vector<char> largest_component(const Walkway& W) {
  const std::size_t n = W.junctions();
  std::vector<std::size_t> comp(n, none), size;
  for (std::size_t r = 0; r < n; r++) {
    if (comp[r] != none) continue;
    std::vector<std::size_t> st{r};
    comp[r] = size.size();
    std::size_t c = 0;
    while (!st.empty()) {
      const std::size_t j = st.back();
      st.pop_back();
      c++;
      for (const std::size_t s : W.leaving[j])
        if (comp[W.head(s)] == none) comp[W.head(s)] = size.size(), st.push_back(W.head(s));
    }
    size.push_back(c);
  }
  const std::size_t big = size.empty() ? 0 : static_cast<std::size_t>(std::max_element(size.begin(), size.end()) - size.begin());
  std::vector<char> mask(n, 0);
  for (std::size_t j = 0; j < n; j++) mask[j] = comp[j] == big;
  return mask;
}

// The nearest segment to a point (optionally of a named street): its index, arc length from its tail and distance.
struct Nearest {
  std::size_t edge = none;
  double along = 0, distance = std::numeric_limits<double>::infinity();
};
inline Nearest nearest_segment(const Walkway& W, const Point& p, const std::string& street = "") {
  Nearest best;
  for (std::size_t e = 0; e < W.edges(); e++) {
    if (!street.empty() && W.segments[e].name != street) continue;
    const auto& l = W.segments[e].line;
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

// Tobler's hiking function (Tobler 1993, "Three presentations on geographical analysis and modeling", NCGIA Technical
// Report 93-1, p. 2): speed 6 exp(−3.5 |s + 0.05|) km/h on a grade s; here as a factor of the speed on the flat.
inline double tobler(double grade) { return std::exp(-3.5 * (std::fabs(grade + 0.05) - 0.05)); }

}  // namespace urban::movement
