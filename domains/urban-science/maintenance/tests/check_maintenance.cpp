// Checks of the maintenance of surfaces (maintenance/deterioration.hpp, planning.hpp) on the Hoddle Grid's road
// surfaces (tests/fixtures/road-condition-hoddle-grid.csv, from the City of Melbourne's open data; make_fixtures.py).
//
// - The deterioration of kerbs, channels, footways and carriageways, fitted to their ratings at their ages, read with
//   the declared error: a family of rates (one, geometric, free) chosen by AIC among those that converge with finite
//   errors; its mean years in each class are reported, with Pearson's χ² of the ratings in bands of five years.
// - The pipeline recovers known rates: classes simulated at the carriageways' real ages from their fitted chain, read
//   with the declared error, and refitted, within three standard errors.
// - Each type's surface model is a POMDP, and every surface's belief now a law.
// - The city's plan for this year, under a budget of 1 % of the stock's renewal: it keeps the budget, and its
//   Lagrangian bound lies above the exact value of doing nothing (a plan within any budget).

#include <urban/maintenance/deterioration.hpp>
#include <urban/maintenance/planning.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace um = urban::maintenance;

namespace {

// Days from 1970-01-01 of an ISO date (proleptic Gregorian; Hinnant's days_from_civil).
long days(const std::string& iso) {
  int y = std::stoi(iso.substr(0, 4)), m = std::stoi(iso.substr(5, 2)), d = std::stoi(iso.substr(8, 2));
  y -= m <= 2;
  const long era = (y >= 0 ? y : y - 399) / 400;
  const long yoe = y - era * 400, doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1, doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
  return era * 146097 + doe - 719468;
}

std::vector<um::Rated> read(const std::string& path, const std::string& now) {
  std::ifstream in(path);
  std::vector<um::Rated> rows;
  std::string line;
  bool header = true;
  const long today = days(now);
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    if (header) { header = false; continue; }
    std::stringstream ss(line);
    std::string type, age, cls, rating, area, evaluated;
    std::getline(ss, type, ','), std::getline(ss, age, ','), std::getline(ss, cls, ','), std::getline(ss, rating, ','), std::getline(ss, area, ','), std::getline(ss, evaluated, ',');
    um::Rated r;
    r.type = type, r.age = age.empty() ? std::nan("") : std::stod(age), r.cls = std::stoi(cls), r.rating = std::stod(rating), r.area = area.empty() ? 0 : std::stod(area);
    r.since = (today - days(evaluated)) / 365.2425;
    rows.push_back(r);
  }
  return rows;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const auto rows = read(std::string(URBAN_FIXTURES) + "/road-condition-hoddle-grid.csv", "2026-11-04");
  s.require("the fixture: the Hoddle Grid's 9630 rated road surfaces", rows.size() == 9630, std::to_string(rows.size()));

  const std::vector<std::string> types = {"Road Kerb", "Road Channel", "Road Footway", "Road Carriageway"};
  std::map<std::string, um::Deterioration> fitted;
  for (const auto& t : types) {
    const auto D = um::fit(rows, t);
    fitted[t] = D;
    char years[260];
    const char* family = D.family == um::Family::common ? "one rate" : D.family == um::Family::geometric ? "rates in geometric progression" : "four free rates";
    std::snprintf(years, sizeof years, "%zu surfaces, %s (AIC of one, geometric, free: %.1f, %.1f, %.1f); mean years in class 5, 4, 3, 2: %.1f, %.1f, %.1f, %.1f", D.surfaces, family,
                  D.tried[0].second, D.tried[1].second, D.tried[2].second, D.residence(5), D.residence(4), D.residence(3), D.residence(2));
    const auto A = um::adequacy(D, rows);
    s.require(t + ": the deterioration fitted (a family converging with finite errors and rates in [10⁻⁴, 10] a year)", D.converged,
              std::string(years) + "; Pearson's χ² " + std::to_string(A.chi2) + " on " + std::to_string(A.df) + " degrees of freedom");
  }

  {
    // Recovery: classes simulated at the carriageways' real ages from their fitted chain, then refitted.
    const auto& D = fitted["Road Carriageway"];
    calx::random::Stream rng(20261104, 0x6d61696eu, 0, 0);  // "main"
    std::vector<um::Rated> sim;
    for (const auto& r : rows) {
      if (r.type != D.type || !std::isfinite(r.age)) continue;
      double t = 0;
      std::size_t k = 0;
      while (k < 4) {
        t += rng.exponential() / D.rate[k];
        if (t > r.age) break;
        k++;
      }
      // Read with the declared error, as the fit assumes.
      const double v = rng.uniform();
      std::size_t j = k;
      if (k > 0 && v < 0.1) j = k - 1;                         // one class better, unless already the best
      else if (k + 1 < 5 && v >= 0.1 && v < 0.2) j = k + 1;    // one class worse, unless already the worst
      um::Rated x = r;
      x.cls = um::class_of(j);
      sim.push_back(x);
    }
    const auto E = um::fit_family(sim, D.type, D.family);
    bool within = E.converged;
    for (std::size_t k = 0; within && k < D.theta.size(); k++) within = std::fabs(E.theta[k] - D.theta[k]) <= 3 * E.error[k];
    std::string got;
    for (std::size_t k = 0; k < D.theta.size(); k++) got += std::to_string(E.theta[k]) + " ± " + std::to_string(E.error[k]) + " (true " + std::to_string(D.theta[k]) + ") ";
    s.require("the fit recovers the carriageways' rates from ratings simulated at their real ages, within three standard errors of its family's parameters", within, got);
  }

  const um::Assumptions A;
  {
    bool valid = true, laws = true;
    for (const auto& [t, D] : fitted) {
      try {
        tyche::decision::validate(um::surface_model(D, A));
      } catch (const std::exception&) {
        valid = false;
      }
    }
    for (const auto& r : rows) {
      if (!fitted.count(r.type)) continue;
      const auto b = um::belief_now(fitted[r.type], r.cls, r.since, A);
      double sum = 0;
      for (const double x : b) sum += x, laws = laws && x >= -1e-15;
      laws = laws && std::fabs(sum - 1) < 1e-12;
    }
    s.require("each type's surface model is a POMDP", valid);
    s.require("every surface's belief now is a law", laws);
  }

  {
    double area = 0;
    for (const auto& r : rows)
      if (fitted.count(r.type)) area += r.area;
    const double budget = 0.01 * A.cost[um::renew] * area;
    const auto P = um::plan_city(rows, fitted, A, budget);
    std::string counts;
    for (const auto& [t, c] : P.counts) counts += t + ": inspect " + std::to_string(c[1]) + ", maintain " + std::to_string(c[2]) + ", renew " + std::to_string(c[3]) + "; ";
    s.require("the city's plan keeps this year's budget of 1 % of the stock's renewal", P.spent <= P.budget + 1e-6,
              std::to_string(P.spent) + " of " + std::to_string(P.budget) + " (" + counts + ")");
    s.require("its Lagrangian bound lies above the exact value of doing nothing", P.upper >= P.nothing_value,
              std::to_string(P.upper) + " at λ = " + std::to_string(P.lambda) + " against " + std::to_string(P.nothing_value));
  }
  return s.done();
}
