// Checks of the space usage index (usage/usage.hpp).
//
// - The shares of a day sum to one over the places reporting, and do not change when every count of a day is doubled: the index removes what
//   moves the whole record.
// - A place with no record on a day stays missing, and the others' shares are of those that do report.
// - The balanced panel keeps only the places reporting on every day of a period.
// - The baseline: by hand on a short record; Z is zero at the mean and one at a standard deviation above it.
// - A place that never varies has no departure to report, rather than an infinite one.

#include <urban/usage/usage.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <limits>
#include <vector>

namespace uu = urban::usage;

int main() {
  calx::verification::Suite s;
  const double no = std::numeric_limits<double>::quiet_NaN();
  const uu::Matrix counts = {{10, 20, 30, 40}, {30, 60, 30, no}, {60, 120, 40, 60}};
  const auto I = uu::shares(counts);
  s.check("a day's shares sum to one", I[0][0] + I[1][0] + I[2][0] - 1, 1e-12);
  s.check("the share of a place", I[0][2] - 0.3, 1e-12);
  s.check("doubling a whole day changes nothing", I[0][1] - I[0][0], 1e-12);
  s.require("a place with no record stays missing", std::isnan(I[1][3]));
  s.check("and the others' shares are of those that report", I[0][3] - 0.4, 1e-12);
  const auto P = uu::balanced(counts, 0, 4);
  s.require("the balanced panel keeps those reporting on every day", P.size() == 2 && P[0] == 0 && P[1] == 2);
  s.require("over a shorter period, more", uu::balanced(counts, 0, 3).size() == 3);
  const auto B = uu::baseline_of(I, {true, true, true, false});
  const double a = 0.1, b = 0.1, c = 0.3, mean = (a + b + c) / 3;
  const double sd = std::sqrt(((a - mean) * (a - mean) + (b - mean) * (b - mean) + (c - mean) * (c - mean)) / 2);
  s.check("the baseline's mean, by hand", B.mean[0] - mean, 1e-12);
  s.check("its standard deviation, by hand", B.sd[0] - sd, 1e-12);
  const auto Z = uu::departures(I, B);
  s.check("Z on a baseline day, by hand", Z[0][2] - (c - mean) / sd, 1e-12);
  s.check("Z on the day beyond it, by hand", Z[0][3] - (0.4 - mean) / sd, 1e-12);
  const auto flat = uu::departures({{0.2, 0.2, 0.2}}, uu::baseline_of({{0.2, 0.2, 0.2}}, {true, true, true}));
  s.require("a place that never varies reports no departure", std::isnan(flat[0][1]));
  return s.done();
}
