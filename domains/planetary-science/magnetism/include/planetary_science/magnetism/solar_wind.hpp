// The interplanetary magnetic field at a body near the Earth, carried from L1 (ADR 0052;
// theory/magnetism.md §4).
//
// SWPC's real-time magnetometer at the Sun–Earth L1 point measures the interplanetary field minute
// by minute, in GSE axes (x toward the Sun, z the ecliptic's pole). By Taylor's frozen-in hypothesis
// the field is advected with the solar wind unchanged, so a body at distance x∥ sunward of the Earth
// along the Sun–Earth line sees the field measured at L1 a travel time τ = (x_L1 − x∥)/v before;
// a body on the night side has x∥ < 0 and a longer travel. Inside the Earth's magnetotail there is
// no solar wind; the tail is taken as a cylinder of 25 Earth radii about the anti-Sun line.
//
// Caveat. Ballistic carriage ignores the tilt of the field's phase fronts, which moves the arrival by
// minutes (Mailyan, Munteanu and Haaland 2008); the tail's shape is a cylinder only to first order,
// and its lobes' field (some 10 nT) is not given. These belong to a domain of space weather when it
// is founded (README §4.4); until then the carriage of the field to a body is this one's.

#pragma once

#include <array>
#include <cmath>
#include <optional>
#include <vector>

namespace planetary_science::magnetism {

using Vec3 = std::array<double, 3>;

inline constexpr double EARTH_RADIUS_M = 6371.0e3;
inline constexpr double MAGNETOTAIL_RADIUS_M = 25 * EARTH_RADIUS_M;
inline constexpr double L1_DISTANCE_M = 1.5e9;  // from the Earth toward the Sun

// Where a body stands toward the solar wind, and the field's travel from L1 to it.
struct Carriage {
  double along_m;       // x∥: the body's distance sunward of the Earth along the Sun–Earth line
  double off_axis_m;    // its distance from that line
  bool in_magnetotail;  // on the night side, within the tail's radius of the line
  double delay_s;       // (x_L1 − x∥)/v
};

// From the Sun's and the body's positions relative to the Earth, in one frame (m), and the solar
// wind's speed (m/s).
inline Carriage carry_from_l1(const Vec3& sun_from_earth, const Vec3& body_from_earth, double speed_m_s) {
  const double s = std::sqrt(sun_from_earth[0] * sun_from_earth[0] + sun_from_earth[1] * sun_from_earth[1] + sun_from_earth[2] * sun_from_earth[2]);
  const Vec3 x{sun_from_earth[0] / s, sun_from_earth[1] / s, sun_from_earth[2] / s};
  const double along = body_from_earth[0] * x[0] + body_from_earth[1] * x[1] + body_from_earth[2] * x[2];
  const Vec3 lateral{body_from_earth[0] - along * x[0], body_from_earth[1] - along * x[1], body_from_earth[2] - along * x[2]};
  const double off = std::sqrt(lateral[0] * lateral[0] + lateral[1] * lateral[1] + lateral[2] * lateral[2]);
  return {along, off, along < 0 && off < MAGNETOTAIL_RADIUS_M, (L1_DISTANCE_M - along) / speed_m_s};
}

// GSE's axes in the frame of the inputs: x toward the Sun, y = z × x in the ecliptic, z toward the
// ecliptic's pole and orthogonal to x.
struct GseAxes {
  Vec3 x, y, z;
};
inline GseAxes gse_axes(const Vec3& sun_from_earth, const Vec3& ecliptic_pole) {
  const auto unit = [](Vec3 v) {
    const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    return Vec3{v[0] / n, v[1] / n, v[2] / n};
  };
  const auto cross = [](const Vec3& a, const Vec3& b) { return Vec3{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; };
  const Vec3 x = unit(sun_from_earth);
  const Vec3 y = unit(cross(ecliptic_pole, x));
  return {x, y, cross(x, y)};
}
// A vector given in GSE, in the frame of the axes.
inline Vec3 from_gse(const GseAxes& a, const Vec3& v) {
  return {v[0] * a.x[0] + v[1] * a.y[0] + v[2] * a.z[0], v[0] * a.x[1] + v[1] * a.y[1] + v[2] * a.z[1], v[0] * a.x[2] + v[1] * a.y[2] + v[2] * a.z[2]};
}

// The sample [t, Bx, By, Bz] of a series nearest an instant t, if within max_gap_s of it.
inline std::optional<std::array<double, 4>> nearest_sample(const std::vector<std::array<double, 4>>& series, double t, double max_gap_s) {
  if (series.empty()) return std::nullopt;
  const std::array<double, 4>* best = &series.front();
  for (const auto& row : series)
    if (std::fabs(row[0] - t) < std::fabs((*best)[0] - t)) best = &row;
  if (std::fabs((*best)[0] - t) > max_gap_s) return std::nullopt;
  return *best;
}

// North, east and down at a planetocentric latitude and longitude, of a vector on the body's fixed axes.
inline Vec3 north_east_down(const Vec3& b, double latitude, double longitude) {
  const double sl = std::sin(latitude), cl = std::cos(latitude), so = std::sin(longitude), co = std::cos(longitude);
  return {-sl * co * b[0] - sl * so * b[1] + cl * b[2], -so * b[0] + co * b[1], -(cl * co * b[0] + cl * so * b[1] + sl * b[2])};
}

}  // namespace planetary_science::magnetism
