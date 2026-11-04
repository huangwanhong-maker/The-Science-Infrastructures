// Sunlight on the ground of a district: the study's own provisional model, not a domain's.
//
// The architecture domain's `environment` holds the Sun's gains on a building's faces but leaves shadow and daylight
// planned, and its gains ignore shading by neighbours (its README §2; heat.hpp, "a stated limit"). This study needs the
// winter sun that reaches the streets and open spaces, so it writes the simplest model that answers it, here, marked
// as its own, and reports the gap (README §8). It belongs in architecture's `environment` when that module's shadow
// and daylight are written.
//
// The model. Buildings are prisms on flat ground: a footprint (a ring in the local frame, east and north in metres)
// standing from 0 to its height h. The Sun is a direction: azimuth A from north through east and altitude α above the
// horizon (the astronomy domain's, at each instant). A point p of the ground is in the sun when the ray from p toward
// the Sun, p + t (sin A, cos A, tan α) for t ≥ 0 (t the horizontal distance), clears every prism. Because the ray rises
// as it goes, it meets a prism exactly when its height where it first enters the prism's footprint is below the
// prism's top: t_in tan α < h, with t_in the least t at which the horizontal ray crosses the footprint's outer ring.
// So the shading of p at an instant by a building of height h is the test h > z, with the threshold z = t_in tan α
// fixed by the geometry alone: the thresholds are computed once, and a design (its heights) is evaluated by comparisons.
//
// The ground is sampled on a square lattice of spacing `spacing` (each sample standing for spacing² of ground); the
// sunlit area-hours of a set of samples are the samples in the sun, summed over the instants, times spacing² and the
// span each instant stands for. Left out: diffuse light from the sky, light reflected from facades, trees, awnings and
// verandahs, the ground's own relief, and buildings beyond the context box. Exact for prisms on flat ground up to the
// lattice's resolution.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace study {

using Point = std::array<double, 2>;
using Ring = std::vector<Point>;

// The least t ≥ 0 at which the ray p + t u crosses the ring, or +∞.
inline double ray_entry(const Point& p, const Point& u, const Ring& r) {
  double best = std::numeric_limits<double>::infinity();
  for (std::size_t k = 0; k < r.size(); k++) {
    const Point& a = r[k];
    const Point& b = r[(k + 1) % r.size()];
    const double ex = b[0] - a[0], ey = b[1] - a[1];
    const double den = u[0] * ey - u[1] * ex;  // cross(u, e)
    if (std::fabs(den) < 1e-12) continue;
    const double wx = a[0] - p[0], wy = a[1] - p[1];
    const double t = (wx * ey - wy * ex) / den;   // along the ray
    const double s = (wx * u[1] - wy * u[0]) / den;  // along the edge
    if (t >= 0 && s >= 0 && s <= 1 && t < best) best = t;
  }
  return best;
}

inline bool inside(const Ring& r, const Point& p) {
  bool in = false;
  for (std::size_t i = 0, j = r.size() - 1; i < r.size(); j = i++)
    if ((r[i][1] > p[1]) != (r[j][1] > p[1]) && p[0] < (r[j][0] - r[i][0]) * (p[1] - r[i][1]) / (r[j][1] - r[i][1]) + r[i][0]) in = !in;
  return in;
}

struct Sun {
  double altitude = 0, azimuth = 0;  // radians; azimuth from north through east
  double weight = 0;                 // the hours this instant stands for
};

// A building for the precomputation: its outer ring, and either a fixed height (the context) or the index of a
// parcel of the design (whose height the design gives).
struct Caster {
  Ring ring;
  double height = 0;
  int parcel = -1;
  std::array<double, 4> box{};  // x0, y0, x1, y1
};

struct Shading {
  std::vector<Point> points;
  std::vector<int> owner;   // −1: public ground in every design; i ≥ 0: the ground of parcel i (public when i is open)
  std::vector<Sun> suns;
  double cell = 0;          // m² a sample stands for
  // For each sample and instant (index p · suns + k): whether the context alone shades it, and the parcels that could
  // shade it with their thresholds (the parcel shades it when its height exceeds the threshold).
  std::vector<std::uint8_t> context_shaded;
  std::vector<std::vector<std::pair<std::uint32_t, float>>> blockers;
};

// The thresholds of every sample at every instant. `reach` bounds the horizontal distance searched (m).
inline Shading shading(std::vector<Point> points, std::vector<int> owner, std::vector<Sun> suns, double cell, std::vector<Caster> casters, double reach) {
  Shading S;
  S.points = std::move(points), S.owner = std::move(owner), S.suns = std::move(suns), S.cell = cell;
  for (auto& c : casters) {
    c.box = {1e300, 1e300, -1e300, -1e300};
    for (const auto& q : c.ring) c.box = {std::min(c.box[0], q[0]), std::min(c.box[1], q[1]), std::max(c.box[2], q[0]), std::max(c.box[3], q[1])};
  }
  const std::size_t P = S.points.size(), K = S.suns.size();
  S.context_shaded.assign(P * K, 0);
  S.blockers.assign(P * K, {});
  for (std::size_t p = 0; p < P; p++) {
    const Point& x = S.points[p];
    for (std::size_t k = 0; k < K; k++) {
      const auto& sun = S.suns[k];
      auto& list = S.blockers[p * K + k];
      if (!(sun.altitude > 0)) { S.context_shaded[p * K + k] = 1; continue; }
      const Point u{std::sin(sun.azimuth), std::cos(sun.azimuth)};
      const double tana = std::tan(sun.altitude);
      const Point far{x[0] + reach * u[0], x[1] + reach * u[1]};
      const double bx0 = std::min(x[0], far[0]), bx1 = std::max(x[0], far[0]), by0 = std::min(x[1], far[1]), by1 = std::max(x[1], far[1]);
      for (const auto& c : casters) {
        if (c.parcel >= 0 && c.parcel == S.owner[p]) continue;  // a parcel's own ground is open only when it is
        if (c.box[2] < bx0 || c.box[0] > bx1 || c.box[3] < by0 || c.box[1] > by1) continue;
        const double t = ray_entry(x, u, c.ring);
        if (!std::isfinite(t) || t > reach) continue;
        const double z = t * tana;
        if (c.parcel < 0) {
          if (c.height > z) S.context_shaded[p * K + k] = 1;
        } else {
          list.push_back({static_cast<std::uint32_t>(c.parcel), static_cast<float>(z)});
        }
      }
    }
  }
  return S;
}

// Whether sample p is in the sun at instant k under the parcels' heights (0 for a parcel given to open space).
inline bool sunlit(const Shading& S, std::size_t p, std::size_t k, const std::vector<double>& height) {
  const std::size_t K = S.suns.size();
  if (S.context_shaded[p * K + k]) return false;
  for (const auto& [i, z] : S.blockers[p * K + k])
    if (height[i] > z) return false;
  return true;
}

// The hours of sun of each sample over the instants (NaN for ground that is not public in this design).
inline std::vector<double> sun_hours(const Shading& S, const std::vector<double>& height, const std::vector<char>& open) {
  std::vector<double> out(S.points.size(), std::numeric_limits<double>::quiet_NaN());
  for (std::size_t p = 0; p < S.points.size(); p++) {
    if (S.owner[p] >= 0 && !open[static_cast<std::size_t>(S.owner[p])]) continue;
    double h = 0;
    for (std::size_t k = 0; k < S.suns.size(); k++)
      if (sunlit(S, p, k, height)) h += S.suns[k].weight;
    out[p] = h;
  }
  return out;
}

}  // namespace study
