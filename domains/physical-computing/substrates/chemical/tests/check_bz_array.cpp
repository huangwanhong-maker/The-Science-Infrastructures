// Checks of coupled BZ cells (bz_array.hpp).
//
// - Two cells whose f differ by 5%, uncoupled, keep their own periods: the gap between their peaks drifts.
// - Coupled through the activator, they lock: the gap between their peaks settles (its spread over the last
//   peaks below 1% of a period).

#include <physical_computing/chemical/bz_array.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace physical_computing;

int main() {
  calx::verification::Suite s;
  chemical::Oregonator base;
  const auto hopf = chemical::hopf_points_in_f(base);
  const double fm = 0.5 * (hopf[0] + hopf[1]);
  chemical::Oregonator a = base, b = base;
  a.f = fm, b.f = fm * 1.05;
  auto gaps = [&](double coupling) {
    chemical::BZArray arr({a, b}, {{0, 1, coupling}});
    const auto peaks = arr.run(300.0, 2e-3);
    std::vector<double> g;
    if (peaks[0].size() < 8 || peaks[1].size() < 8) return std::pair{g, 0.0};
    const std::size_t n = std::min(peaks[0].size(), peaks[1].size());
    for (std::size_t k = n > 12 ? n - 12 : 0; k < n; k++) g.push_back(peaks[1][k] - peaks[0][k]);
    const double period = (peaks[0].back() - peaks[0][peaks[0].size() - 6]) / 5;
    return std::pair{g, period};
  };
  {
    const auto [g, period] = gaps(0.0);
    s.require("uncoupled: both cells oscillate, and the gap between their peaks drifts", !g.empty() && std::fabs(g.back() - g.front()) > 0.05 * period);
  }
  {
    const auto [g, period] = gaps(2.0);
    s.require("coupled: both cells oscillate", !g.empty());
    if (!g.empty()) {
      const auto [lo, hi] = std::minmax_element(g.begin(), g.end());
      s.check("coupled: the cells lock, the gap's spread over the last peaks (in periods)", (*hi - *lo) / period, 0.01);
    }
  }
  return s.done();
}
