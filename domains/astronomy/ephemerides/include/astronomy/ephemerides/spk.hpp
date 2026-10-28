// The JPL planetary ephemerides: positions and velocities of the Sun, the Moon and the planets from
// SPK kernels such as DE440 (Park et al. 2021, AJ 161, 105).
//
// An SPK is a DAF (daf.hpp) whose segments each give one body (the target) relative to another
// (the centre) over an interval, in a frame (1: the ICRF, "J2000" in SPICE), as Chebyshev series
// of position (type 2) or of position and velocity (type 3) on records of fixed length (NAIF,
// "SPK Required Reading", N0065, §§ "Type 2", "Type 3"). The time argument is TDB in seconds
// past J2000.0, "ephemeris time". A body's barycentric state is the sum along the chain of its
// segments (the Moon: 301 relative to the Earth-Moon barycentre 3, 3 relative to the
// barycentre 0). Positions in km and velocities in km/s, as the kernels hold them, are given here
// in metres and m/s, SI as the infrastructure keeps everything.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <astronomy/ephemerides/daf.hpp>

namespace astronomy::ephemerides {

using Vec3 = std::array<double, 3>;

struct State {
  Vec3 position{0, 0, 0};  // m
  Vec3 velocity{0, 0, 0};  // m/s
};

class Spk {
 public:
  struct Segment {
    int target = 0, center = 0, frame = 0, type = 0;
    std::int64_t begin = 0, end = 0;
    double start_et = 0, end_et = 0;
    double init = 0, interval = 0;
    int record_size = 0, records = 0, coefficients = 0;
    std::string name;
  };

  explicit Spk(Daf daf) : daf_(std::move(daf)) {
    if (daf_.id().rfind("DAF/SPK", 0) != 0) throw std::invalid_argument("SPK: not an SPK file (\"" + daf_.id() + "\")");
    if (daf_.nd() != 2 || daf_.ni() != 6) throw std::invalid_argument("SPK: unexpected summary format");
    for (const DafSummary& s : daf_.summaries()) {
      Segment g;
      g.start_et = s.doubles[0], g.end_et = s.doubles[1];
      g.target = s.ints[0], g.center = s.ints[1], g.frame = s.ints[2], g.type = s.ints[3];
      g.begin = s.ints[4], g.end = s.ints[5];
      g.name = s.name;
      if (g.type == 2 || g.type == 3) {
        g.init = daf_.at(g.end - 3);
        g.interval = daf_.at(g.end - 2);
        g.record_size = int(daf_.at(g.end - 1));
        g.records = int(daf_.at(g.end));
        g.coefficients = (g.record_size - 2) / (g.type == 2 ? 3 : 6);
      }
      segments_.push_back(g);
    }
  }
  static Spk read(const std::string& path) { return Spk(Daf::read(path)); }

  const std::vector<Segment>& segments() const { return segments_; }
  // The interval every body of the kernel covers, TDB seconds past J2000.
  std::pair<double, double> coverage() const {
    double a = -1e300, b = 1e300;
    for (const Segment& s : segments_) a = std::max(a, s.start_et), b = std::min(b, s.end_et);
    return {a, b};
  }
  bool covers(double et) const {
    const auto [a, b] = coverage();
    return et >= a && et <= b;
  }

  // The state of a segment's target relative to its centre, ICRF, at ET.
  State segment_state(const Segment& g, double et) const {
    if (g.type != 2 && g.type != 3) throw std::invalid_argument("SPK: segment type " + std::to_string(g.type) + " is not supported");
    if (g.frame != 1) throw std::invalid_argument("SPK: only segments in the ICRF (frame 1) are read");
    int index = int(std::floor((et - g.init) / g.interval));
    index = std::clamp(index, 0, g.records - 1);
    const double* r = daf_.doubles(g.begin + std::int64_t(index) * g.record_size);
    const double mid = r[0], radius = r[1];
    const double s = (et - mid) / radius;
    const int n = g.coefficients;
    State st;
    for (int k = 0; k < 3; k++) {
      const auto [p, dp] = chebyshev(r + 2 + k * n, n, s);
      st.position[k] = p * 1e3;
      st.velocity[k] = (g.type == 2 ? dp / radius : chebyshev(r + 2 + (3 + k) * n, n, s).first) * 1e3;
    }
    return st;
  }
  const Segment& segment_for(int target, double et) const {
    for (const Segment& s : segments_)
      if (s.target == target && et >= s.start_et && et <= s.end_et) return s;
    throw std::domain_error("SPK: no segment for body " + std::to_string(target) + " at ET " + std::to_string(et));
  }
  // The state of a body relative to the solar-system barycentre (NAIF id 0).
  State barycentric(int target, double et) const {
    State sum;
    int body = target;
    int guard = 0;
    while (body != 0) {
      if (++guard > 16) throw std::logic_error("SPK: the chain of centres does not reach the barycentre");
      const Segment& g = segment_for(body, et);
      const State s = segment_state(g, et);
      for (int k = 0; k < 3; k++) sum.position[k] += s.position[k], sum.velocity[k] += s.velocity[k];
      body = g.center;
    }
    return sum;
  }
  // The geometric state of a target relative to an observer, ICRF axes, at ET. The two chains of
  // segments are summed only up to their nearest common centre, so that the Moon relative to the
  // Earth never passes through the Earth-Moon barycentre's 1.5 × 10¹¹ m from the Sun: the result
  // keeps the precision of its own size (10⁻⁸ m for the Moon, not 10⁻⁵ m).
  State state(int target, int observer, double et) const {
    std::vector<int> up_t{target}, up_o{observer};
    for (int b = target, guard = 0; b != 0 && guard < 16; guard++) up_t.push_back(b = segment_for(b, et).center);
    for (int b = observer, guard = 0; b != 0 && guard < 16; guard++) up_o.push_back(b = segment_for(b, et).center);
    int common = 0;
    for (int b : up_t)
      if (std::find(up_o.begin(), up_o.end(), b) != up_o.end()) {
        common = b;
        break;
      }
    auto sum_to = [&](int body, int stop) {
      State acc;
      while (body != stop) {
        const Segment& g = segment_for(body, et);
        const State s = segment_state(g, et);
        for (int k = 0; k < 3; k++) acc.position[k] += s.position[k], acc.velocity[k] += s.velocity[k];
        body = g.center;
      }
      return acc;
    };
    const State a = sum_to(target, common), b = sum_to(observer, common);
    State d;
    for (int k = 0; k < 3; k++) d.position[k] = a.position[k] - b.position[k], d.velocity[k] = a.velocity[k] - b.velocity[k];
    return d;
  }
  bool has(int body) const {
    for (const Segment& s : segments_)
      if (s.target == body) return true;
    return body == 0;
  }

 private:
  Daf daf_;
  std::vector<Segment> segments_;
};

}  // namespace astronomy::ephemerides
