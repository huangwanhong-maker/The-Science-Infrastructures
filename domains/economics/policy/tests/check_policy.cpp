// Checks of regimes perceived from growth and of the decisions on them (regimes/regimes.hpp, policy/policy.hpp).
//
// - Regimes: 120 economies' growth simulated for 60 years from three regimes (one series in ten starting late, 5 %
//   missing): the pooled fit recovers the regimes' mean growth within 0.3 and their persistence within 0.05; every
//   reading's filtered and smoothed laws sum to 1; the recession's mean duration is 1 / (1 − A₀₀).
// - Economies of different normals and spreads: standardised, the pooled fit's paths recover each one's regimes.
// - The policymaker's POMDP is a POMDP; its growth bins are laws; support shortens a recession and restraint a boom,
//   as declared; with certainty of recession HSVI's policy supports, with certainty of a boom it restrains, with
//   certainty of normal growth it waits; its interval closes.
// - The fund's model keeps two actions; the cost of support is one a year of support.

#include <economics/policy/policy.hpp>
#include <economics/regimes/regimes.hpp>

#include <tyche/decision/point_based.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

int main() {
  calx::verification::Suite s;
  tyche::processes::GaussianHmm truth;
  truth.states = 3;
  truth.initial = {0.1, 0.7, 0.2};
  truth.transition = {0.45, 0.5, 0.05, 0.08, 0.82, 0.10, 0.05, 0.35, 0.60};
  truth.mean = {-2.5, 3.0, 7.0};
  truth.sd = {2.0, 1.5, 2.0};
  calx::random::Stream rng(20261104, 0x72656769u, 0, 0);  // "regi"
  std::vector<economics::GrowthSeries> series;
  for (int e = 0; e < 120; e++) {
    economics::GrowthSeries g{"E" + std::to_string(e), 1961, {}};
    std::size_t k = 1;
    for (int t = 0; t < 60; t++) {
      const double u = rng.uniform();
      double c = 0;
      for (std::size_t j = 0; j < 3; j++)
        if ((c += truth.a(k, j)) > u) { k = j; break; }
      const bool late = e % 10 == 0 && t < 25;
      g.growth.push_back(late || rng.uniform() < 0.05 ? std::nan("") : rng.normal(truth.mean[k], truth.sd[k]));
    }
    series.push_back(g);
  }
  const auto R = economics::perceive_regimes(series);
  bool ok = R.converged;
  std::string got;
  for (std::size_t k = 0; k < 3; k++) {
    ok = ok && std::fabs(R.model.mean[k] - truth.mean[k]) < 0.3 && std::fabs(R.model.a(k, k) - truth.a(k, k)) < 0.05;
    got += std::to_string(R.model.mean[k]) + " (" + std::to_string(R.model.a(k, k)) + ") ";
  }
  s.require("the pooled fit recovers three regimes' mean growth and persistence, numbered from the worst", ok, got);
  double worst = 0;
  for (const auto& r : R.readings)
    for (std::size_t t = 0; t < r.growth.size(); t++) {
      double f = 0, m = 0;
      for (std::size_t k = 0; k < 3; k++) f += r.filtered[t][k], m += r.smoothed[t][k];
      worst = std::max(worst, std::max(std::fabs(f - 1), std::fabs(m - 1)));
    }
  s.check("every reading's filtered and smoothed laws sum to 1", worst, 1e-12);
  s.check("the recession's mean duration is 1 / (1 − A₀₀)", std::fabs(R.duration(0) - 1 / (1 - R.model.a(0, 0))), 1e-15);
  {
    // Economies of different normals and spreads, each from the same regimes in its own deviations: standardised, the
    // pooled fit's most probable paths recover each economy's regimes.
    std::vector<economics::GrowthSeries> mixed;
    std::vector<std::vector<std::size_t>> truth_path;
    const double zmean[3] = {-1.8, 0.0, 1.3};
    for (int e = 0; e < 120; e++) {
      const double centre = rng.uniform(1, 7), scale = rng.uniform(0.6, 3);
      economics::GrowthSeries g{"M" + std::to_string(e), 1961, {}};
      std::vector<std::size_t> path;
      std::size_t k = 1;
      for (int t = 0; t < 60; t++) {
        const double u = rng.uniform();
        double c = 0;
        for (std::size_t j = 0; j < 3; j++)
          if ((c += truth.a(k, j)) > u) { k = j; break; }
        g.growth.push_back(centre + scale * rng.normal(zmean[k], 0.5));
        path.push_back(k);
      }
      mixed.push_back(g), truth_path.push_back(path);
    }
    const auto S = economics::perceive_regimes(mixed, 3, 10, true);
    double agree = 0, all = 0;
    for (std::size_t e = 0; e < S.readings.size(); e++)
      for (std::size_t t = 0; t < 60; t++) agree += S.readings[e].path[t] == truth_path[e][t], all += 1;
    s.require("standardised, the pooled fit's most probable paths recover each economy's own regimes in 85 % of its years", S.standardised && agree / all > 0.85,
              std::to_string(100 * agree / all) + " %");
  }

  economics::PolicyAssumptions A;
  A.edges = {-2, 0, 2, 4, 6};  // per cent: these regimes are fitted on raw growth
  const auto P = economics::policy_pomdp(R.model, A);
  bool valid = true;
  try {
    tyche::decision::validate(P);
  } catch (const std::exception&) {
    valid = false;
  }
  s.require("the policymaker's POMDP is a POMDP (its growth bins laws)", valid);
  s.require("support shortens a recession, restraint a boom, as declared",
            P.transition(economics::support, 0, 0) < P.transition(economics::wait, 0, 0) && P.transition(economics::restrain, 2, 2) < P.transition(economics::wait, 2, 2));
  const auto H = tyche::decision::hsvi(P, {.epsilon = 0.01});
  const std::vector<double> recession{1, 0, 0}, normal{0, 1, 0}, boom{0, 0, 1};
  s.require("certain of recession it supports, of a boom it restrains, of normal growth it waits",
            H.lower.action(recession) == economics::support && H.lower.action(boom) == economics::restrain && H.lower.action(normal) == economics::wait,
            std::to_string(H.lower.action(recession)) + std::to_string(H.lower.action(normal)) + std::to_string(H.lower.action(boom)));
  s.require("HSVI's interval closes to 0.01", H.converged, "[" + std::to_string(H.lower_value) + ", " + std::to_string(H.upper_value) + "]");
  const auto F = economics::policy_pomdp(R.model, A, 2);
  const auto C = economics::support_cost(P);
  s.require("the fund's model keeps wait and support; support costs one a year", F.actions == 2 && C[economics::support * 3 + 1] == 1 && C[economics::wait * 3 + 1] == 0);
  return s.done();
}
