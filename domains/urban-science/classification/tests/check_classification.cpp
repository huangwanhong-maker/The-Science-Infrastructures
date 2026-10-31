// Checks of the classification (theory/street-oceanography.md §8): two tight groups of segments read as
// clustered; a regular lattice of segments as regular; segments scattered uniformly at random as unstructured
// in most of 40 draws (the test at α = 0.05); fewer than five segments as too few.

#include <urban/classification/plane.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <string>
#include <vector>

using namespace urban::classification;

int main() {
  calx::verification::Suite s;
  calx::random::Stream rng(6, 6, 6, 6);
  {
    std::vector<Segment> g;
    for (int i = 0; i < 20; i++) g.push_back({std::to_string(i), "", 1 + 0.02 * rng.normal(), 0.2 + 0.01 * rng.normal()});
    for (int i = 0; i < 20; i++) g.push_back({std::to_string(20 + i), "", 3 + 0.02 * rng.normal(), 0.8 + 0.01 * rng.normal()});
    const auto p = plane(g);
    s.require("classification: two tight groups are clustered", p.reading == "clustered", p.reading + ", R = " + std::to_string(p.test.ratio));
  }
  {
    std::vector<Segment> g;
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++) g.push_back({"", "", static_cast<double>(i), static_cast<double>(j)});
    const auto p = plane(g);
    s.require("classification: a lattice of segments is regular", p.reading == "regular", p.reading + ", R = " + std::to_string(p.test.ratio));
  }
  {
    int unstructured = 0;
    for (int k = 0; k < 40; k++) {
      std::vector<Segment> g;
      for (int i = 0; i < 40; i++) g.push_back({"", "", rng.uniform(), rng.uniform()});
      unstructured += plane(g).reading == "unstructured";
    }
    s.require("classification: uniform scatter reads as unstructured in most draws", unstructured >= 30, std::to_string(unstructured) + " of 40");
  }
  {
    const auto p = plane({{"a", "", 1, 1}, {"b", "", 2, 2}});
    s.require("classification: two segments are too few", p.reading == "too few");
  }
  return s.done();
}
