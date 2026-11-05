// The market for the centre's dwellings: households of several classes choosing among segments of dwellings (a type
// in an area) by a multinomial logit, the rents that clear it, the utilities that the observed rents and stock reveal,
// and the households' surplus.
//
// Households and segments. A household of class c (an income class, a tenure, a size; its share w_c of the M
// households who consider the centre) values a dwelling of segment j (type t(j) in area a(j)) at
//
//   V_cj = δ_j + γ_{c,t(j)} − α_c ln(r_j / r̄) − θ d_j + ε_cj,   α_c = α (ȳ / y_c)^π,
//
// with δ_j the segment's own value, γ the class's taste for the type (−∞ where the class may not live in it: student
// apartments are for students), r_j the weekly rent, d_j the distance (km) to the centre's anchor, and ε independent
// Gumbel (McFadden 1974, "Conditional logit analysis of qualitative choice behavior", in Zarembka (ed.), Frontiers in
// Econometrics, 105–142). Not choosing the centre (living elsewhere in the metropolis) is the outside option, V = 0, and
// r̄ is the rent of that elsewhere (the reference rent). Rent enters as the logarithm of its ratio to it, so that α is
// the semi-elasticity of a segment's choice with respect to its rent in proportion, free of the general level of
// prices; poorer households are more sensitive (π = 1: in proportion to their income). Then
//
//   P_cj = e^{V_cj} / (1 + Σ_k e^{V_ck}),    D_j = M Σ_c w_c P_cj.
//
// Clearing. The rents clear the market when every segment's demand equals its occupied stock, D_j(r) = (1 − v) S_j at
// a natural vacancy v. The logit's demands are gross substitutes (∂D_j/∂r_j < 0 < ∂D_j/∂r_k), so that the clearing
// rents are unique when they exist (they do when M exceeds the occupied stock of the segments each class may enter);
// they are found by Newton's method on each segment's log demand, z_j ← z_j + (ln D_j − ln O_j) / ᾱ_j, z = ln r and
// ᾱ_j = Σ_c α_c D_cj (1 − P_cj) / D_j its own semi-elasticity (for one class this is Berry's contraction below, read in
// rents). With one class, the clearing rents are explicit: ln(O_j / (M − Σ_k O_k)) = δ_j − α ln(r_j / r̄) − θ d_j.
//
// Inversion (Berry 1994, "Estimating discrete-choice models of product differentiation", RAND J. Econ. 25, 242–262).
// Given the rents and the occupied stock, the δ that make them clear are the fixed point of the contraction
// δ ← δ + ln O − ln D(δ): a segment's own value is what its rent, at its stock, reveals.
//
// Surplus (Small and Rosen 1981, "Applied welfare economics with discrete choice models", Econometrica 49, 105–130).
// A class's expected maximum utility is ln(1 + Σ_j e^{V_cj}); divided by α_c it is in units of ln r: a change of it
// is the proportional change of rent that would leave the class as well off. The surplus in A$ a week multiplies it
// by the class's mean rent paid. By Roy's identity ∂/∂ln r_j [M w_c ln(1 + Σ e^V)/α_c] = −M w_c P_cj.
//
// Caveat. One rent per segment: the dwellings of a segment are alike. A segment's rent is the median of its new
// lettings, as Homes Victoria reports it, not the rent of sitting tenants.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace urban::housing {

namespace detail {
inline void require(bool ok, const std::string& what) {
  if (!ok) throw std::invalid_argument("urban::housing: " + what);
}
}  // namespace detail

// CLUE's types of dwelling (the residential-dwellings dataset's classification).
enum Type : std::size_t { house = 0, apartment = 1, student = 2 };
constexpr std::size_t types = 3;

struct Segment {
  std::size_t type = apartment, area = 0;
  double supply = 0;    // dwellings
  double base = 0;      // δ_j
  double distance = 0;  // km to the anchor
};

struct HouseholdClass {
  double weight = 1;           // w_c, the class's share of the M households
  double income = 1;           // y_c, A$ a week
  std::vector<double> taste;   // γ_{c,t}, one per type; −∞: may not live in the type
};

struct Market {
  std::vector<Segment> segments;
  std::vector<HouseholdClass> classes;
  double households = 0;       // M: households who consider the centre
  double alpha = 1;            // α at the mean income
  double mean_income = 1;      // ȳ
  double income_power = 1;     // π
  double distance_cost = 0;    // θ, per km
  double vacancy = 0.05;       // the natural vacancy v
  double reference_rent = 1;   // r̄, the outside option's rent

  double rent_coefficient(std::size_t c) const { return alpha * std::pow(mean_income / classes[c].income, income_power); }
  double utility(std::size_t c, std::size_t j, double rent) const {
    const auto& s = segments[j];
    const double g = classes[c].taste.empty() ? 0.0 : classes[c].taste[s.type];
    if (!std::isfinite(g)) return -std::numeric_limits<double>::infinity();
    return s.base + g - rent_coefficient(c) * std::log(rent / reference_rent) - distance_cost * s.distance;
  }
  double occupied(std::size_t j) const { return (1 - vacancy) * segments[j].supply; }
};

// P_c over the segments, the outside option last (J + 1 entries), each exp after subtracting the largest utility.
inline std::vector<double> choice(const Market& m, std::size_t c, std::span<const double> rents) {
  const std::size_t J = m.segments.size();
  std::vector<double> v(J + 1), p(J + 1, 0.0);
  double top = 0;  // the outside option's utility
  for (std::size_t j = 0; j < J; j++) top = std::max(top, v[j] = m.utility(c, j, rents[j]));
  v[J] = 0;
  double total = 0;
  for (std::size_t j = 0; j <= J; j++) total += (p[j] = std::isfinite(v[j]) ? std::exp(v[j] - top) : 0.0);
  for (double& x : p) x /= total;
  return p;
}

// D_j = M Σ_c w_c P_cj.
inline std::vector<double> demand(const Market& m, std::span<const double> rents) {
  const std::size_t J = m.segments.size();
  std::vector<double> d(J, 0.0);
  for (std::size_t c = 0; c < m.classes.size(); c++) {
    const auto p = choice(m, c, rents);
    for (std::size_t j = 0; j < J; j++) d[j] += m.households * m.classes[c].weight * p[j];
  }
  return d;
}

struct Clearing {
  std::vector<double> rents;
  int iterations = 0;
  double residual = 0;   // max_j |ln D_j − ln O_j|
  bool converged = false;
};

// The rents at which D_j = (1 − v) S_j for every segment with stock, from `start` (or 1 A$ a week); a segment without
// stock keeps its starting rent. Newton steps on the log rents are capped at one unit.
inline Clearing clear(const Market& m, std::vector<double> start = {}, double tolerance = 1e-12, int max_iterations = 100000) {
  const std::size_t J = m.segments.size();
  detail::require(!m.classes.empty() && J > 0, "a market needs classes and segments");
  double occupied = 0;
  for (std::size_t j = 0; j < J; j++) occupied += m.occupied(j);
  detail::require(m.households > occupied, "the households who consider the centre must outnumber its occupied dwellings");
  Clearing out;
  out.rents = start.empty() ? std::vector<double>(J, 1.0) : std::move(start);
  std::vector<double> z(J);
  for (std::size_t j = 0; j < J; j++) z[j] = std::log(out.rents[j]);
  for (out.iterations = 0; out.iterations < max_iterations; out.iterations++) {
    for (std::size_t j = 0; j < J; j++) out.rents[j] = std::exp(z[j]);
    std::vector<double> d(J, 0.0), slope(J, 0.0);
    for (std::size_t c = 0; c < m.classes.size(); c++) {
      const auto p = choice(m, c, out.rents);
      const double a = m.rent_coefficient(c), w = m.households * m.classes[c].weight;
      for (std::size_t j = 0; j < J; j++) d[j] += w * p[j], slope[j] += a * w * p[j] * (1 - p[j]);
    }
    double worst = 0;
    std::vector<double> step(J, 0.0);
    for (std::size_t j = 0; j < J; j++) {
      const double o = m.occupied(j);
      if (!(o > 0)) continue;
      detail::require(d[j] > 0, "a segment with stock that no class may enter");
      const double gap = std::log(d[j]) - std::log(o);
      worst = std::max(worst, std::fabs(gap));
      step[j] = std::clamp(gap / (slope[j] / d[j]), -1.0, 1.0);
    }
    out.residual = worst;
    if (worst < tolerance) { out.converged = true; break; }
    for (std::size_t j = 0; j < J; j++) z[j] += step[j];
  }
  return out;
}

// The segments' own values δ at which the given rents clear the given occupied stock (Berry's contraction).
struct Inversion {
  std::vector<double> base;
  int iterations = 0;
  double residual = 0;
  bool converged = false;
};

inline Inversion invert(Market m, std::span<const double> rents, std::span<const double> occupied, double tolerance = 1e-12, int max_iterations = 100000) {
  const std::size_t J = m.segments.size();
  double total = 0;
  for (std::size_t j = 0; j < J; j++) total += occupied[j];
  detail::require(m.households > total, "the households who consider the centre must outnumber its occupied dwellings");
  Inversion out;
  for (out.iterations = 0; out.iterations < max_iterations; out.iterations++) {
    const auto d = demand(m, rents);
    double worst = 0;
    for (std::size_t j = 0; j < J; j++) {
      if (!(occupied[j] > 0)) continue;
      detail::require(d[j] > 0, "a segment with stock that no class may enter");
      const double gap = std::log(occupied[j]) - std::log(d[j]);
      worst = std::max(worst, std::fabs(gap));
      m.segments[j].base += gap;
    }
    out.residual = worst;
    if (worst < tolerance) { out.converged = true; break; }
  }
  for (const auto& s : m.segments) out.base.push_back(s.base);
  return out;
}

// Σ_c M w_c ln(1 + Σ_j e^{V_cj}) / α_c: the households' expected maximum utility in units of ln r (the outside
// option's utility zero).
inline double surplus_log_rent(const Market& m, std::span<const double> rents) {
  double s = 0;
  for (std::size_t c = 0; c < m.classes.size(); c++) {
    double top = 0;
    std::vector<double> v(m.segments.size());
    for (std::size_t j = 0; j < v.size(); j++) top = std::max(top, v[j] = m.utility(c, j, rents[j]));
    double total = std::exp(-top);
    for (const double x : v)
      if (std::isfinite(x)) total += std::exp(x - top);
    s += m.households * m.classes[c].weight * (top + std::log(total)) / m.rent_coefficient(c);
  }
  return s;
}

// The surplus in A$ a week: each class's expected maximum utility over α_c, times the mean rent of the segments it
// chooses (the rent at which a proportional change of rent is a change in dollars).
inline double surplus_dollars(const Market& m, std::span<const double> rents) {
  double s = 0;
  for (std::size_t c = 0; c < m.classes.size(); c++) {
    const auto p = choice(m, c, rents);
    double inside = 0, paid = 0;
    for (std::size_t j = 0; j < m.segments.size(); j++) inside += p[j], paid += p[j] * rents[j];
    if (!(inside > 0)) continue;
    double top = 0;
    std::vector<double> v(m.segments.size());
    for (std::size_t j = 0; j < v.size(); j++) top = std::max(top, v[j] = m.utility(c, j, rents[j]));
    double total = std::exp(-top);
    for (const double x : v)
      if (std::isfinite(x)) total += std::exp(x - top);
    s += m.households * m.classes[c].weight * (top + std::log(total)) / m.rent_coefficient(c) * (paid / inside);
  }
  return s;
}

}  // namespace urban::housing
