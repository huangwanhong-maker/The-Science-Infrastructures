// Checks of the coarsening foam (foam.hpp).
//
// - A Voronoi start of 300 seeds on a 96 × 96 torus has (nearly) 300 cells with a mean of six sides (within 0.5).
// - Coarsening, the cells of few sides shrink and of many grow: the mean dA/dt of cells of at most 4 sides is
//   negative and of at least 8 positive; a line through the classes with 10 cells or more rises and crosses zero
//   between 5 and 7 sides (von Neumann's law).
// - The mean area grows.

#include <physical_computing/fluids/foam.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

using namespace physical_computing;

int main() {
  calx::verification::Suite s;
  fluids::FoamOptions o;
  o.size = 96, o.cells = 300, o.seed = 4;
  fluids::Foam foam(o);
  const auto c0 = foam.cells();
  s.require("foam: the Voronoi start has nearly all its cells", c0.size() >= 295);
  double sides = 0;
  for (const auto& [id, c] : c0) sides += static_cast<double>(c.sides);
  s.check("foam: a mean of six sides at the start", sides / static_cast<double>(c0.size()) - 6, 0.5);
  const double a0 = foam.mean_area();
  foam.coarsen(30);
  const auto before = foam.cells();
  foam.coarsen(60);
  const auto after = foam.cells();
  const auto rates = fluids::Foam::area_rates(before, after, 60);
  double few = 0, many = 0;
  std::size_t nf = 0, nm = 0;
  std::vector<double> xs, ys, ws;
  for (const auto& [n, r] : rates) {
    if (n <= 4) few += r.first * r.second, nf += r.second;
    if (n >= 8) many += r.first * r.second, nm += r.second;
    if (r.second >= 10) xs.push_back(static_cast<double>(n)), ys.push_back(r.first), ws.push_back(static_cast<double>(r.second));
  }
  s.require("foam: cells of at most four sides shrink", nf > 0 && few / nf < 0);
  s.require("foam: cells of at least eight sides grow", nm > 0 && many / nm > 0);
  double sw = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (std::size_t i = 0; i < xs.size(); i++) sw += ws[i], sx += ws[i] * xs[i], sy += ws[i] * ys[i], sxx += ws[i] * xs[i] * xs[i], sxy += ws[i] * xs[i] * ys[i];
  const double slope = (sw * sxy - sx * sy) / (sw * sxx - sx * sx), intercept = (sy - slope * sx) / sw, zero = -intercept / slope;
  s.require("foam: the rate rises with the number of sides", slope > 0);
  s.check("foam: and crosses zero at six sides (von Neumann)", zero - 6, 1.0);
  s.check("foam: a mean of six sides after coarsening", [&] { double t = 0; for (const auto& [id, c] : after) t += static_cast<double>(c.sides); return t / static_cast<double>(after.size()) - 6; }(), 0.5);
  s.require("foam: the mean area grows", foam.mean_area() > a0);
  return s.done();
}
