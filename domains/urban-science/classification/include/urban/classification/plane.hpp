// The classification, the field's first test (theory/street-oceanography.md, Definition 8.1; Huang, Street
// Oceanography, Claim 29.2 and Figure 2): segments placed in the plane of compression ratio C and selectivity
// of response S, and whether they form groups.
//
// The axes are scaled to [0, 1] by their ranges over the segments placed (the classification is invariant to
// a common L_walk, Remark 2.5, and the scaling removes the units' arbitrariness), and Clark and Evans's ratio
// R of the mean nearest-neighbour distance to its expectation under complete spatial randomness on the unit
// square, with Donnelly's correction for its edge (tyche `spatial`), tests the arrangement:
//   - R < 1 significantly: the segments cluster — the vocabulary has sorted its object;
//   - R ≈ 1: they fill the plane without structure — the paper's first condition of failure (§29.3);
//   - R > 1: they spread regularly.
// Caveat 8.2: segments not selected by the design of §26.2 test the vocabulary, not the field's design.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <tyche/spatial/point_patterns.hpp>

namespace urban::classification {

struct Segment {
  std::string id, name;
  double compression = 0;   // C
  double selectivity = 0;   // S
};

struct Plane {
  std::vector<Segment> segments;           // those with both coordinates finite
  std::vector<double> x, y;                // scaled to [0, 1]
  double c_min = 0, c_max = 0, s_min = 0, s_max = 0;
  tyche::spatial::ClarkEvans test;
  std::string reading;                     // clustered | unstructured | regular | too few
};

inline Plane plane(const std::vector<Segment>& all, double alpha = 0.05) {
  Plane p;
  for (const auto& s : all)
    if (std::isfinite(s.compression) && std::isfinite(s.selectivity)) p.segments.push_back(s);
  if (p.segments.size() < 5) { p.reading = "too few"; return p; }
  auto range = [&](auto get, double& lo, double& hi) {
    lo = hi = get(p.segments[0]);
    for (const auto& s : p.segments) lo = std::min(lo, get(s)), hi = std::max(hi, get(s));
    if (hi == lo) hi = lo + 1;
  };
  range([](const Segment& s) { return s.compression; }, p.c_min, p.c_max);
  range([](const Segment& s) { return s.selectivity; }, p.s_min, p.s_max);
  for (const auto& s : p.segments) {
    p.x.push_back((s.compression - p.c_min) / (p.c_max - p.c_min));
    p.y.push_back((s.selectivity - p.s_min) / (p.s_max - p.s_min));
  }
  p.test = tyche::spatial::clark_evans(p.x, p.y, 0, 1, 0, 1);
  p.reading = p.test.p >= alpha ? "unstructured" : (p.test.ratio < 1 ? "clustered" : "regular");
  return p;
}

}  // namespace urban::classification
