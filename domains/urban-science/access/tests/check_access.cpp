// Checks of reach over a network of ways (access/access.hpp).
//
// - Ways meeting at a shared point are joined; ways that only cross, with no shared point (a bridge), are not.
// - A grid of streets of known spacing at walking speed: every time is the Manhattan distance over the speed.
// - The walk to the network counts, and it enters at the nearest point.
// - A faster way is taken when it is faster, a closed way never.
// - The length within a time: the whole network at a long limit, a known part at a short one, proportionately along a piece.

#include <urban/access/access.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

namespace ua = urban::access;

int main() {
  calx::verification::Suite s;
  {
    const auto n = ua::network_of({{{0, 0, 100, 0}, 1}, {{100, 0.2, 100, 100}, 1}, {{50, -50, 50, 50}, 1}});
    s.require("ways sharing a point within the snapping distance are joined", n.x.size() == 5, std::to_string(n.x.size()));
    const auto r = ua::reach(n, {{0, 0}}, 1);
    s.require("one reaches along them", std::isfinite(r.time[n.points[1][1]]));
    s.require("a way that only crosses them, sharing no point, is not joined", !std::isfinite(r.time[n.points[2][0]]));
  }
  {
    std::vector<ua::Way> ways;
    const int m = 9;
    const double d = 120;
    for (int i = 0; i < m; i++) {
      ua::Way a, b;
      a.speed = b.speed = 1.5;
      for (int j = 0; j < m; j++) { a.xy.push_back(i * d); a.xy.push_back(j * d); b.xy.push_back(j * d); b.xy.push_back(i * d); }
      ways.push_back(a);
      ways.push_back(b);
    }
    const auto n = ua::network_of(ways);
    const auto r = ua::reach(n, {{4 * d, 4 * d}}, 1.5);
    double worst = 0;
    for (std::size_t p = 0; p < n.x.size(); p++) worst = std::max(worst, std::abs(r.time[p] - (std::abs(n.x[p] - 4 * d) + std::abs(n.y[p] - 4 * d)) / 1.5));
    s.check("on a grid of streets every time is the Manhattan distance over the speed", worst, 1e-9);
    const auto within = ua::length_within(n, r.time, {1e6, 120 / 1.5, 60 / 1.5});
    const double total = 2.0 * m * (m - 1) * d;
    s.check("at a long limit the whole network is within reach", within[0] / total - 1, 1e-12);
    s.check("within one block, the four arms from the crossing", within[1] / (4 * d) - 1, 1e-12);
    s.check("within half a block, in proportion along them", within[2] / (2 * d) - 1, 1e-12);
    const auto far = ua::reach(n, {{4 * d - 50, 4 * d}}, 1.0);
    s.check("the walk to the nearest point counts, at its own speed", far.time[far.start] - 50, 1e-9);
    s.check("and it is reported", far.start_gap - 50, 1e-9);
    const auto L = ua::reach(n, {{0, 0}}, 1.5, 200);
    int reached = 0;
    for (const double t : L.time) reached += std::isfinite(t);
    s.require("a limit stops the search", reached == 6, std::to_string(reached));
  }
  {
    // A slow straight way, and a fast way round it.
    const auto n = ua::network_of({{{0, 0, 300, 0}, 1}, {{0, 0, 0, 100, 300, 100, 300, 0}, 10}});
    const auto r = ua::reach(n, {{0, 0}}, 1);
    s.check("the faster way round is taken", r.time[n.points[0][1]] - 500.0 / 10, 1e-9);
    const auto closed = ua::network_of({{{0, 0, 300, 0}, 1}, {{0, 0, 0, 100, 300, 100, 300, 0}, 0}});
    const auto c = ua::reach(closed, {{0, 0}}, 1);
    s.check("a closed way is never taken", c.time[closed.points[0][1]] - 300, 1e-9);
  }
  {
    // Two ways cut at a tile's edge, 3 m apart at their ends.
    auto n = ua::network_of({{{0, 0, 100, 0}, 1}, {{103, 0, 200, 0}, 1}});
    s.require("ways cut apart are not joined at first", !std::isfinite(ua::reach(n, {{0, 0}}, 1).time[n.points[1][1]]));
    s.require("their loose ends are stitched within the radius", ua::stitch_ends(n, 8, 1) >= 1);
    s.check("and one walks across at the stitching's speed", ua::reach(n, {{0, 0}}, 1).time[n.points[1][1]] - 200, 1e-9);
    auto m = ua::network_of({{{0, 0, 100, 0}, 1}, {{120, 0, 200, 0}, 1}});
    s.require("beyond the radius they stay apart", ua::stitch_ends(m, 8, 1) == 0);
  }
  return s.done();
}
