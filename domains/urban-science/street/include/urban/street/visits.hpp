// The visits of kinds of people (theory/street-oceanography.md, Representation 2.14): which kinds go to which
// holdings, and when. A representation, declared and not fitted (Street Oceanography, §24): it shows how a
// terrain read through unlike weightings sorts arrivals in time as well as in space.
//
// Arrivals of type a are present on the street with a weekly profile N_a(θ) (hourly, 168 values, people per
// hour), enter at s_a, and choose among the holdings open at θ by a multinomial logit (McFadden 1974):
//   V_a(h, θ) = N_a(θ) 1_{O_h}(θ) exp(w_a[u_h] − β_a |s_a − s_h|) / Σ_{h′} 1_{O_{h′}}(θ) exp(w_a[u_{h′}] − β_a |s_a − s_{h′}|).
// An hour with no holding open sends no visits (the arrivals pass through: [SO] §35.1, refusal).
//
// Each holding's weekly visits C_h(θ) = Σ_a V_a(h, θ) have lines Ĉ_h(m) = (1/168) Σ_θ C_h(θ) e^{−2πimθ/168}
// (θ at the hours' centres); the share of their power off the daily harmonics measures how much the
// holding's clientele follows the week rather than the day.

#pragma once

#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include <urban/street/street.hpp>

namespace urban::street {

struct ArrivalType {
  std::string name;
  std::vector<double> presence;   // N_a per hour of the week (168), people per hour
  std::vector<double> affinity;   // w_a per use (index of Street::uses), dimensionless utility
  double distance_cost = 0.005;   // β_a, 1/m
  double entry = 0;               // s_a, m
};

struct Visits {
  std::vector<std::string> types;
  // visits[a][h][θ]: expected visits of type a to holding h in hour θ of the week.
  std::vector<std::vector<std::vector<double>>> visits;
  std::vector<std::vector<double>> by_holding;      // [h][a]: weekly total of type a at h
  std::vector<std::vector<double>> share;           // [h][a]: composition of h's clientele
  std::vector<double> unserved;                     // [a]: arrivals in hours when nothing is open
  std::vector<double> off_daily;                    // [h]: the share of C_h's line power (m ≥ 1) off the daily harmonics
};

inline Visits visits(const Street& st, const std::vector<ArrivalType>& types) {
  const std::size_t H = st.holdings.size(), A = types.size();
  Visits out;
  out.visits.assign(A, std::vector<std::vector<double>>(H, std::vector<double>(168, 0)));
  out.by_holding.assign(H, std::vector<double>(A, 0));
  out.share.assign(H, std::vector<double>(A, 0));
  out.unserved.assign(A, 0);
  for (std::size_t a = 0; a < A; a++) {
    const auto& T = types[a];
    if (T.presence.size() != 168 || T.affinity.size() != st.uses.size())
      throw std::invalid_argument("urban::street::visits: a presence for each hour of the week and an affinity for each use");
    out.types.push_back(T.name);
    for (std::size_t k = 0; k < 168; k++) {
      const double theta = static_cast<double>(k) + 0.5;
      std::vector<double> weight(H, 0);
      double total = 0;
      for (std::size_t h = 0; h < H; h++) {
        const auto& g = st.holdings[h];
        if (!g.declared || !open_at(g, theta)) continue;
        weight[h] = std::exp(T.affinity[g.use] - T.distance_cost * std::fabs(T.entry - g.s));
        total += weight[h];
      }
      if (total <= 0) { out.unserved[a] += T.presence[k]; continue; }
      for (std::size_t h = 0; h < H; h++) {
        const double v = T.presence[k] * weight[h] / total;
        out.visits[a][h][k] = v;
        out.by_holding[h][a] += v;
      }
    }
  }
  out.off_daily.assign(H, 0);
  for (std::size_t h = 0; h < H; h++) {
    double total = 0;
    for (double v : out.by_holding[h]) total += v;
    for (std::size_t a = 0; a < A; a++) out.share[h][a] = total > 0 ? out.by_holding[h][a] / total : 0;
    // Lines of the weekly series, m = 1…84.
    double daily = 0, all = 0;
    for (std::size_t m = 1; m <= 84; m++) {
      std::complex<double> c{0, 0};
      for (std::size_t k = 0; k < 168; k++) {
        double v = 0;
        for (std::size_t a = 0; a < A; a++) v += out.visits[a][h][k];
        c += v * std::exp(std::complex<double>(0, -2 * std::numbers::pi * static_cast<double>(m) * (static_cast<double>(k) + 0.5) / 168.0));
      }
      const double p = std::norm(c / 168.0);
      all += p;
      if (m % 7 == 0) daily += p;
    }
    out.off_daily[h] = all > 0 ? 1 - daily / all : 0;
  }
  return out;
}

}  // namespace urban::street
