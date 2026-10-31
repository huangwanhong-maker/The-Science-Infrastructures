// Currents past a counting line (theory/street-oceanography.md §7.4). A sensor that counts each way is
// read as a current meter's record.
//
// A sensor counts d₁ walkers toward its first direction and d₂ toward its second in each hour. Their mean at
// each civil hour of the week is the record's projection onto the mean and the week's 84 harmonics (exactly so
// for whole weeks of record). From it come:
//   - the transport each way, and the net current u = d₁ − d₂ (people an hour, positive toward the first
//     direction);
//   - the residual current, the mean of u over the week, as the residual of a tidal current is its mean
//     over the tide's cycles (Pugh and Woodworth 2014, Sea-Level Science, ch. 4);
//   - the hours of greatest current each way on the mean weekday: the flood and the ebb;
//   - the irreversibility of the exchange.
//
// The irreversibility. Take the walkers as independent and each a Markov jump process among places. N such
// walkers jump from i to j n_ij = N p_i k_ij times an hour on average. Their entropy production, in units
// of k_B an hour, is Schnakenberg's
//   σ = ½ Σ_{i≠j} (n_ij − n_ji) ln(n_ij / n_ji)
// (Schnakenberg 1976, Rev. Mod. Phys. 48, 571, eq. 7.6, times N; at each instant for a driven process, as the
// total entropy production of the master equation, Seifert 2012, Rep. Prog. Phys. 75, 126001). This is a
// sum over edges, each term (n₊ − n₋) ln(n₊/n₋) ≥ 0, zero exactly when the edge balances. A counting line
// measures one edge's two directed rates, so its term is that edge's share of σ. The terms of the measured
// edges sum to a lower bound on the σ of the whole network of streets, since no term is negative. physis's
// irreversibility computes the same σ from a rate matrix, and the checks hold the two equal on a ring.
// Computed from the hour's mean counts: ln of a small count is biased, so an hour with fewer than one walker
// either way has no term.
//
// The handedness of the walkers. A sensor on one footpath of a street counts along it, and its residual
// current runs one way or the other. Suppose walkers chose their side of the street as left-hand traffic
// chooses its lane. Then the residual on each footpath would run with the street on its right: east on the
// north side, west on the south side, north on the west side, south on the east side. With no handedness,
// each footpath's residual runs either way with probability ½. The number of footpaths agreeing with the
// left-hand rule is then binomial(n, ½) (a sign test; tyche's binomial). The footpaths are not independent:
// two footpaths of one street share its walkers, which the test does not model.
//
// Checks: tests/check_flow.cpp.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <tyche/probability/distributions.hpp>

#include <urban/rhythm/calendar.hpp>
#include <urban/rhythm/record.hpp>

namespace urban::flow {

// One edge's term of the entropy production, (n₊ − n₋) ln(n₊/n₋), for directed rates n₊, n₋ > 0.
inline double edge_entropy_production(double forward, double backward) {
  if (!(forward > 0) || !(backward > 0)) throw std::invalid_argument("urban::flow::edge_entropy_production: both directed rates above zero");
  return (forward - backward) * (std::log(forward) - std::log(backward));
}

struct CurrentMeter {
  std::array<std::vector<double>, 2> way;   // the mean count toward each direction at each civil hour of the week (168)
  std::vector<double> net;                  // u = way₀ − way₁
  std::vector<double> irreversibility;      // (way₀ − way₁) ln(way₀/way₁); NaN below one walker either way
  std::vector<double> weeks;                // values averaged at each hour (where both ways are recorded)
  double residual = 0;                      // the mean of u
  double transport = 0;                     // the mean of way₀ + way₁
  double sigma = 0;                         // the mean of the irreversibility over the hours that have it
  double flood_at = 0, flood = 0;           // the mean weekday's greatest u (toward the first direction) and its civil hour
  double ebb_at = 0, ebb = 0;               // its least u (greatest toward the second direction) and its hour
};

// The current meter of two records on the same hours (the counts toward each direction), over [from, to) UTC.
inline CurrentMeter current_meter(const rhythm::Record& d1, const rhythm::Record& d2, double from = -1e18, double to = 1e18) {
  if (d1.size() != d2.size() || d1.t0 != d2.t0 || d1.step != d2.step) throw std::invalid_argument("urban::flow::current_meter: two records on the same hours");
  CurrentMeter m;
  std::array<std::vector<double>, 2> sum{std::vector<double>(168, 0), std::vector<double>(168, 0)};
  m.weeks.assign(168, 0);
  for (std::size_t i = 0; i < d1.size(); i++) {
    const double t = d1.time(i);
    if (t < from || t >= to || !std::isfinite(d1.y[i]) || !std::isfinite(d2.y[i])) continue;
    const auto k = static_cast<std::size_t>(std::floor(rhythm::week_angle(d1.civil(t)) * 168 + 1e-9)) % 168;
    sum[0][k] += d1.y[i], sum[1][k] += d2.y[i], m.weeks[k] += 1;
  }
  for (std::size_t k = 0; k < 168; k++)
    if (!(m.weeks[k] > 0)) throw std::invalid_argument("urban::flow::current_meter: every hour of the week recorded");
  m.net.resize(168), m.irreversibility.resize(168);
  for (int w = 0; w < 2; w++) {
    m.way[w].resize(168);
    for (std::size_t k = 0; k < 168; k++) m.way[w][k] = sum[w][k] / m.weeks[k];
  }
  double sig = 0, n = 0;
  for (std::size_t k = 0; k < 168; k++) {
    m.net[k] = m.way[0][k] - m.way[1][k];
    m.residual += m.net[k] / 168, m.transport += (m.way[0][k] + m.way[1][k]) / 168;
    m.irreversibility[k] = m.way[0][k] >= 1 && m.way[1][k] >= 1 ? edge_entropy_production(m.way[0][k], m.way[1][k]) : std::nan("");
    if (std::isfinite(m.irreversibility[k])) sig += m.irreversibility[k], n += 1;
  }
  m.sigma = n > 0 ? sig / n : std::nan("");
  // The mean weekday, hour by hour.
  m.flood = -std::numeric_limits<double>::infinity(), m.ebb = std::numeric_limits<double>::infinity();
  for (std::size_t h = 0; h < 24; h++) {
    double u = 0;
    for (std::size_t d = 0; d < 5; d++) u += m.net[24 * d + h] / 5;
    if (u > m.flood) m.flood = u, m.flood_at = static_cast<double>(h) + 0.5;
    if (u < m.ebb) m.ebb = u, m.ebb_at = static_cast<double>(h) + 0.5;
  }
  return m;
}

enum class Compass { north, east, south, west };

struct Handedness {
  std::size_t footpaths = 0;   // footpaths whose side is across the direction they count, with a residual
  std::size_t left = 0;        // those whose residual runs as left-hand traffic would
  double p_left = 1;           // P(X ≥ left), X ~ binomial(footpaths, ½)
  double p_two_sided = 1;      // P(|X − n/2| ≥ |left − n/2|)
};

// The direction left-hand traffic takes on a side of a street: east on the north side, and so on round.
inline Compass left_hand(Compass side) {
  switch (side) {
    case Compass::north: return Compass::east;
    case Compass::east: return Compass::south;
    case Compass::south: return Compass::west;
    default: return Compass::north;
  }
}

// Footpaths by their side of the street, the two directions they count (toward[0] the positive one) and their
// residual currents. A footpath counting along its side (north or south side, east or west counted, or the
// reverse) with a residual not zero enters the test.
inline Handedness handedness(const std::vector<Compass>& side, const std::vector<std::array<Compass, 2>>& toward, const std::vector<double>& residual) {
  if (side.size() != toward.size() || side.size() != residual.size()) throw std::invalid_argument("urban::flow::handedness: a side, two directions and a residual for each footpath");
  auto across = [](Compass a, Compass b) { return (static_cast<int>(a) + static_cast<int>(b)) % 2 == 1; };
  Handedness h;
  for (std::size_t i = 0; i < side.size(); i++) {
    if (!across(side[i], toward[i][0]) || !(residual[i] != 0) || !std::isfinite(residual[i])) continue;
    h.footpaths++;
    h.left += (residual[i] > 0 ? toward[i][0] : toward[i][1]) == left_hand(side[i]);
  }
  if (h.footpaths == 0) return h;
  const tyche::probability::Binomial b{static_cast<std::int64_t>(h.footpaths), 0.5};
  const auto n = static_cast<std::int64_t>(h.footpaths), k = static_cast<std::int64_t>(h.left);
  h.p_left = b.survival(k - 1);
  const std::int64_t far = std::max(k, n - k);
  h.p_two_sided = std::min(1.0, b.survival(far - 1) + b.cdf(n - far));
  if (2 * far == n) h.p_two_sided = 1;
  return h;
}

}  // namespace urban::flow
