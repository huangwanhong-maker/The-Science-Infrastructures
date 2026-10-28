// Binary planetary constants kernels: the orientation of a body as Euler angles in time, such as
// the Moon's principal-axes frame integrated with the DE440 ephemeris (Park et al. 2021).
//
// A binary PCK is a DAF (daf.hpp) whose type-2 segments give three Euler angles (φ, θ, ψ) of a
// body-fixed frame relative to a base frame (1: the ICRF) as Chebyshev series on records of fixed
// length, the same layout as SPK type 2; the derivatives of the series are the angles' rates
// (NAIF, "PCK Required Reading", N0062, §"Binary PCK Kernel Format"). The change of coordinates
// from the base frame to the body frame is the 3-1-3 sequence of frame rotations
//   M = R_z(ψ) R_x(θ) R_z(φ),   v_body = M v_base,
// with R_z, R_x the frame rotations of SOFA (physis frame_rotation). The time argument is TDB in
// seconds past J2000.0.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <astronomy/ephemerides/daf.hpp>
#include <physis/geometry/kinematics.hpp>
#include <physis/geometry/rotation.hpp>

namespace astronomy::ephemerides {

class Pck {
 public:
  struct Segment {
    int body_frame = 0, base_frame = 0, type = 0;
    std::int64_t begin = 0, end = 0;
    double start_et = 0, end_et = 0, init = 0, interval = 0;
    int record_size = 0, records = 0, coefficients = 0;
  };
  struct Angles {
    std::array<double, 3> angles{};  // φ, θ, ψ, rad
    std::array<double, 3> rates{};   // rad/s
  };

  explicit Pck(Daf daf) : daf_(std::move(daf)) {
    if (daf_.id().rfind("DAF/PCK", 0) != 0) throw std::invalid_argument("PCK: not a binary PCK file (\"" + daf_.id() + "\")");
    if (daf_.nd() != 2 || daf_.ni() != 5) throw std::invalid_argument("PCK: unexpected summary format");
    for (const DafSummary& s : daf_.summaries()) {
      Segment g;
      g.start_et = s.doubles[0], g.end_et = s.doubles[1];
      g.body_frame = s.ints[0], g.base_frame = s.ints[1], g.type = s.ints[2], g.begin = s.ints[3], g.end = s.ints[4];
      if (g.type == 2) {
        g.init = daf_.at(g.end - 3);
        g.interval = daf_.at(g.end - 2);
        g.record_size = int(daf_.at(g.end - 1));
        g.records = int(daf_.at(g.end));
        g.coefficients = (g.record_size - 2) / 3;
      }
      segments_.push_back(g);
    }
  }
  static Pck read(const std::string& path) { return Pck(Daf::read(path)); }
  const std::vector<Segment>& segments() const { return segments_; }

  bool covers(int body_frame, double et) const {
    for (const Segment& s : segments_)
      if (s.body_frame == body_frame && et >= s.start_et && et <= s.end_et) return true;
    return false;
  }
  Angles angles(int body_frame, double et) const {
    for (const Segment& g : segments_) {
      if (g.body_frame != body_frame || et < g.start_et || et > g.end_et) continue;
      if (g.type != 2) throw std::invalid_argument("PCK: segment type " + std::to_string(g.type) + " is not supported");
      if (g.base_frame != 1) throw std::invalid_argument("PCK: only segments relative to the ICRF are read");
      int index = std::clamp(int(std::floor((et - g.init) / g.interval)), 0, g.records - 1);
      const double* r = daf_.doubles(g.begin + std::int64_t(index) * g.record_size);
      const double s = (et - r[0]) / r[1];
      Angles a;
      for (int k = 0; k < 3; k++) {
        const auto [v, d] = chebyshev(r + 2 + k * g.coefficients, g.coefficients, s);
        a.angles[k] = v;
        a.rates[k] = d / r[1];
      }
      return a;
    }
    throw std::domain_error("PCK: no orientation for frame " + std::to_string(body_frame) + " at ET " + std::to_string(et));
  }
  // The change of coordinates from the ICRF to the body frame, with its time derivative.
  physis::geometry::StateTransform icrf_to_body(int body_frame, double et) const {
    using physis::geometry::Axis;
    using physis::geometry::frame_rotation;
    using namespace physis::geometry::operators;
    const Angles a = angles(body_frame, et);
    const double phi = a.angles[0], theta = a.angles[1], psi = a.angles[2];
    const auto z1 = frame_rotation(Axis::z, phi), x = frame_rotation(Axis::x, theta), z2 = frame_rotation(Axis::z, psi);
    // d/dα of a frame rotation about an axis e is −[e]× times it.
    const auto ex = physis::geometry::skew({1, 0, 0}), ez = physis::geometry::skew({0, 0, 1});
    const auto dz1 = (-a.rates[0]) * (ez * z1), dx = (-a.rates[1]) * (ex * x), dz2 = (-a.rates[2]) * (ez * z2);
    physis::geometry::StateTransform t;
    t.r = z2 * x * z1;
    t.r_dot = dz2 * x * z1 + z2 * dx * z1 + z2 * x * dz1;
    return t;
  }

 private:
  Daf daf_;
  std::vector<Segment> segments_;
};

}  // namespace astronomy::ephemerides
