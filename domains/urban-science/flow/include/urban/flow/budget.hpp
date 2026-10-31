// Person-time, the street's exact budget (theory/street-oceanography.md, Theorem 6.1; Huang, Street
// Oceanography, Claim 27.1): for units present on [a_i, d_i), N(t) = Σ_i 1[a_i ≤ t < d_i] and for any window
//   ∫₀^T N(t) dt = Σ_i |[a_i, d_i) ∩ [0, T)|
// exactly (Fubini). Divided by T it is Little's law L = λW with its boundary terms (Little 1961; Kim and Whitt
// 2013, Oper. Res. 61, 1030, for finite windows): the units present at 0 and at T are counted for the part of
// their stay inside the window.
//
// From cordon counts per interval (in, out) and the number present at the start, N follows by accumulation
// and the person-time by the trapezoid inside each interval only if the crossings are spread evenly in it;
// the identity is exact for the step function of the counts (crossings at the interval's end), which is what
// is computed, and the budget's residual is the miscount of the cordons: N(T) − (N(0) + Σ in − Σ out).

#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace urban::flow {

struct PersonTime {
  double integral = 0;            // ∫₀^T N dt (the left side)
  double sum_of_stays = 0;        // Σ_i |[a_i, d_i) ∩ [0, T)| (the right side)
  double arrival_rate = 0;        // λ: arrivals in the window / T
  double mean_present = 0;        // L = integral / T
  double mean_stay = 0;           // W of the units arriving in the window, their whole stays
  double little_gap = 0;          // L − λW: the boundary terms
  std::size_t units = 0;
};

// Units with arrival and departure times (departure ≥ arrival), over the window [0, T).
inline PersonTime person_time(const std::vector<double>& arrive, const std::vector<double>& depart, double T) {
  if (arrive.size() != depart.size() || !(T > 0)) throw std::invalid_argument("urban::flow::person_time: one departure per arrival, a positive window");
  PersonTime p;
  p.units = arrive.size();
  // The left side: N is a step function changing at the events inside the window.
  std::vector<std::pair<double, int>> events;
  int present = 0;
  for (std::size_t i = 0; i < arrive.size(); i++) {
    if (depart[i] < arrive[i]) throw std::invalid_argument("urban::flow::person_time: a departure before its arrival");
    if (arrive[i] <= 0 && depart[i] > 0) present++;
    if (arrive[i] > 0 && arrive[i] < T) events.emplace_back(arrive[i], +1);
    if (depart[i] > 0 && depart[i] < T) events.emplace_back(depart[i], -1);
    p.sum_of_stays += std::max(0.0, std::min(depart[i], T) - std::max(arrive[i], 0.0));
  }
  std::sort(events.begin(), events.end());
  double t = 0;
  for (const auto& [at, change] : events) {
    p.integral += present * (at - t);
    present += change, t = at;
  }
  p.integral += present * (T - t);
  std::size_t in_window = 0;
  double stays = 0;
  for (std::size_t i = 0; i < arrive.size(); i++)
    if (arrive[i] >= 0 && arrive[i] < T) in_window++, stays += depart[i] - arrive[i];
  p.arrival_rate = static_cast<double>(in_window) / T;
  p.mean_present = p.integral / T;
  p.mean_stay = in_window ? stays / static_cast<double>(in_window) : 0;
  p.little_gap = p.mean_present - p.arrival_rate * p.mean_stay;
  return p;
}

struct CordonBudget {
  std::vector<double> present;    // N at the end of each interval
  double person_time = 0;         // Σ N_k Δt (crossings at the intervals' ends)
  double residual = 0;            // N_end (if given) − accumulated N
};

inline CordonBudget cordon_budget(const std::vector<double>& in, const std::vector<double>& out, double present0, double step, double present_end = -1) {
  if (in.size() != out.size()) throw std::invalid_argument("urban::flow::cordon_budget: counts in and out for each interval");
  CordonBudget b;
  double n = present0;
  for (std::size_t k = 0; k < in.size(); k++) {
    b.person_time += n * step;
    n += in[k] - out[k];
    b.present.push_back(n);
  }
  b.residual = present_end >= 0 ? present_end - n : 0;
  return b;
}

}  // namespace urban::flow
