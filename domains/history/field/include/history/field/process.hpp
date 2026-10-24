// The history in time as a sample path (theory note §8.3–8.4).
//
// The occasions are time-ordered: a history is one term of the Dyson series of a theory whose
// quanta are the beings in the record and whose interactions are the occasions. Its weight
// factorises into the law of the occasions, the law of the strands (which ties begin and end) and
// the law of the signs, the last the spin foam's (transfer.hpp). The first two are estimated here
// as a marked point process, by maximum likelihood, with what a single history can support:
//
//   occasions: an occasion with k beings arriving from before (and any number appearing first)
//     happens at rate λ_k per k-subset of the beings in the record: intensity λ_k C(n(τ), k),
//     n(τ) the beings between their first and last occasion. λ̂_k = n_k / E_k, E_k = ∫ C(n, k) dτ,
//     with Poisson standard error λ̂_k/√n_k. Beings appearing for the first time are created at the
//     vertex, from the vacuum, as a Fock-space vertex creates quanta.
//   strands: at an occasion, two continuing beings with no live tie begin one with probability
//     q_open; a live tie met there ends with probability q_close. Binomial maximum likelihood.
//   binding: a pair of beings in the record meets at rate μ₁ while a tie between them is live and
//     μ₀ while none is; ρ = μ₁/μ₀ says how strongly a tie holds its beings together. A dyadic
//     (composite) likelihood: the pairs of one occasion are not independent, so its interval is
//     indicative.
//
// The checks are the theory's own martingales (Dynkin's formula, the stochastic counterpart of
// the Schwinger–Dyson equations): the count of occasions less its compensator ∫ intensity has
// mean zero over any window, with variance the compensator. Its standardised value on each era of
// the history tests the fitted rates where they were not fitted.

#pragma once

#include <history/construction/complex.hpp>
#include <morphe/arithmetic/integer.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace history {

struct Rate {
  std::size_t k = 0;     // beings arriving from before
  std::size_t count = 0; // occasions of this kind
  double exposure = 0;   // ∫ C(n(τ), k) dτ, in years
  double rate = 0, error = 0, upper = 0;  // λ̂, its standard error, and a 95% upper bound when count = 0
  double mean_new = 0;   // beings appearing for the first time, on average
};

struct Proportion {
  std::size_t trials = 0, successes = 0;
  double p = 0, error = 0;
};

struct EraCheck {
  double from = 0, to = 0;                   // years
  std::size_t occasions = 0;                 // observed
  double expected = 0, z = 0;                // the compensator, and (O − E)/√E
  std::size_t opened = 0;
  double opened_expected = 0, opened_z = 0;  // ties begun, against q_open × opportunities
};

struct Process {
  double from = 0, to = 0;                   // the record's span, years
  std::vector<Rate> rates;                   // by k
  Proportion open, close;
  double mu_tied = 0, mu_untied = 0, binding = 0, binding_low = 0, binding_high = 0;
  std::size_t tied_meetings = 0, untied_meetings = 0;
  double tied_exposure = 0, untied_exposure = 0;
  std::vector<EraCheck> eras;
  double chi2 = 0;                           // Σ z² over the eras, occasions
  std::size_t dof = 0;
  double log_likelihood_occasions = 0, log_likelihood_strands = 0;
};

namespace process_detail {

inline Proportion proportion(std::size_t trials, std::size_t successes) {
  Proportion p;
  p.trials = trials;
  p.successes = successes;
  if (trials) {
    p.p = static_cast<double>(successes) / static_cast<double>(trials);
    p.error = std::sqrt(p.p * (1 - p.p) / static_cast<double>(trials));
  }
  return p;
}

inline double bernoulli_log(double p, std::size_t s, std::size_t n) {
  double l = 0;
  if (s) l += static_cast<double>(s) * std::log(p);
  if (n > s) l += static_cast<double>(n - s) * std::log(1 - p);
  return l;
}

}  // namespace process_detail

inline Process process_of(const HistoryComplex& K, std::size_t eras = 3) {
  using namespace process_detail;
  const auto& r = K.record;
  const std::size_t V = K.order.size(), B = r.beings.size(), T = r.ties.size();
  Process P;
  if (V == 0) return P;
  std::vector<double> year(r.occasions.size());
  for (std::size_t o = 0; o < r.occasions.size(); o++) year[o] = historical_year(r.occasions[o].date.first);
  P.from = year[K.order.front()];
  P.to = year[K.order.back()];
  // Each being's time in the record, and each tie's life, in years.
  std::vector<double> enter(B, 0), leave(B, 0);
  for (std::size_t b = 0; b < B; b++)
    if (!K.worldline[b].empty()) { enter[b] = year[K.worldline[b].front()]; leave[b] = year[K.worldline[b].back()]; }
  std::vector<double> born(T, 0), ends(T, 0);
  for (std::size_t t = 0; t < T; t++)
    if (!K.rungs[t].empty()) { born[t] = year[K.rungs[t].front().occasion]; ends[t] = year[K.rungs[t].back().occasion]; }

  // n(τ) and the tied pairs as step functions: integrate C(n, k) and the pair counts over [a, b).
  auto integrate = [&](double a, double b, std::size_t kmax, std::vector<double>& by_k, double& tied, double& pairs) {
    std::vector<double> cuts{a, b};
    for (std::size_t x = 0; x < B; x++) {
      if (enter[x] > a && enter[x] < b) cuts.push_back(enter[x]);
      if (leave[x] > a && leave[x] < b) cuts.push_back(leave[x]);
    }
    for (std::size_t t = 0; t < T; t++) {
      if (born[t] > a && born[t] < b) cuts.push_back(born[t]);
      if (ends[t] > a && ends[t] < b) cuts.push_back(ends[t]);
    }
    std::sort(cuts.begin(), cuts.end());
    by_k.assign(kmax + 1, 0);
    tied = pairs = 0;
    for (std::size_t i = 0; i + 1 < cuts.size(); i++) {
      const double lo = cuts[i], hi = cuts[i + 1], mid = (lo + hi) / 2, w = hi - lo;
      if (w <= 0) continue;
      long long n = 0;
      for (std::size_t x = 0; x < B; x++) n += enter[x] < mid && mid < leave[x];
      for (std::size_t k = 0; k <= kmax; k++) by_k[k] += w * morphe::arithmetic::binomial(n, static_cast<long long>(k)).to_double();
      std::set<std::pair<std::size_t, std::size_t>> live;
      for (std::size_t t = 0; t < T; t++)
        if (!K.rungs[t].empty() && born[t] < mid && mid < ends[t]) live.insert(std::minmax(r.ties[t].a, r.ties[t].b));
      tied += w * static_cast<double>(live.size());
      pairs += w * morphe::arithmetic::binomial(n, 2).to_double();
    }
  };

  // The occasions: how many arrive from before, how many are new; and the strands at each.
  std::vector<std::size_t> k_in(V, 0), k_new(V, 0);
  std::vector<std::size_t> rank_of(r.occasions.size(), 0);
  for (std::size_t v = 0; v < V; v++) rank_of[K.order[v]] = v;
  std::vector<std::size_t> opp_open(V, 0), opened(V, 0);
  std::size_t close_trials = 0, closed = 0;
  std::size_t kmax = 0;
  for (std::size_t v = 0; v < V; v++) {
    const std::size_t o = K.order[v];
    std::vector<std::size_t> continuing;
    for (const auto& p : r.occasions[o].present) {
      const bool arrives = !K.worldline[p.being].empty() && K.worldline[p.being].front() != o;
      (arrives ? k_in[v] : k_new[v])++;
      if (arrives) continuing.push_back(p.being);
    }
    kmax = std::max(kmax, k_in[v]);
    // Pairs of continuing beings: tied (a tie live just before) or not; new ties begun here.
    std::set<std::pair<std::size_t, std::size_t>> tied_before, begun_here;
    for (std::size_t t = 0; t < T; t++) {
      if (K.rungs[t].empty()) continue;
      const auto key = std::minmax(r.ties[t].a, r.ties[t].b);
      const std::size_t first = rank_of[K.rungs[t].front().occasion], last = rank_of[K.rungs[t].back().occasion];
      if (first < v && last >= v) tied_before.insert(key);
      if (first == v) begun_here.insert(key);
    }
    for (std::size_t i = 0; i < continuing.size(); i++)
      for (std::size_t j = i + 1; j < continuing.size(); j++) {
        const auto key = std::minmax(continuing[i], continuing[j]);
        if (tied_before.count(key)) P.tied_meetings++;
        else {
          P.untied_meetings++;
          opp_open[v]++;
          opened[v] += begun_here.count(key);
        }
      }
  }
  for (std::size_t t = 0; t < T; t++)
    for (std::size_t g = 1; g < K.rungs[t].size(); g++) {
      close_trials++;
      closed += g + 1 == K.rungs[t].size();
    }

  // Rates by k.
  std::vector<double> exposure;
  double tied = 0, pairs = 0;
  integrate(P.from, P.to, kmax, exposure, tied, pairs);
  for (std::size_t k = 0; k <= kmax; k++) {
    Rate R;
    R.k = k;
    R.exposure = k == 0 ? P.to - P.from : exposure[k];
    double news = 0;
    for (std::size_t v = 0; v < V; v++)
      if (k_in[v] == k) { R.count++; news += static_cast<double>(k_new[v]); }
    if (R.count == 0 && R.exposure <= 0) continue;
    R.mean_new = R.count ? news / static_cast<double>(R.count) : 0;
    if (R.exposure > 0) {
      R.rate = static_cast<double>(R.count) / R.exposure;
      R.error = R.count ? R.rate / std::sqrt(static_cast<double>(R.count)) : 0;
      R.upper = R.count ? 0 : 3.0 / R.exposure;
      if (R.count) P.log_likelihood_occasions += static_cast<double>(R.count) * (std::log(R.rate) - 1);
    }
    P.rates.push_back(R);
  }
  // Strands.
  std::size_t opp = 0, ope = 0;
  for (std::size_t v = 0; v < V; v++) { opp += opp_open[v]; ope += opened[v]; }
  P.open = proportion(opp, ope);
  P.close = proportion(close_trials, closed);
  if (P.open.trials && P.open.p > 0 && P.open.p < 1) P.log_likelihood_strands += bernoulli_log(P.open.p, ope, opp);
  if (P.close.trials && P.close.p > 0 && P.close.p < 1) P.log_likelihood_strands += bernoulli_log(P.close.p, closed, close_trials);
  // Binding.
  P.tied_exposure = tied;
  P.untied_exposure = pairs - tied;
  if (P.tied_exposure > 0) P.mu_tied = static_cast<double>(P.tied_meetings) / P.tied_exposure;
  if (P.untied_exposure > 0) P.mu_untied = static_cast<double>(P.untied_meetings) / P.untied_exposure;
  if (P.mu_tied > 0 && P.mu_untied > 0) {
    P.binding = P.mu_tied / P.mu_untied;
    const double se = std::sqrt(1.0 / static_cast<double>(P.tied_meetings) + 1.0 / static_cast<double>(P.untied_meetings));
    P.binding_low = P.binding * std::exp(-1.96 * se);
    P.binding_high = P.binding * std::exp(1.96 * se);
  }
  // The martingale checks by era.
  for (std::size_t q = 0; q < eras && P.to > P.from; q++) {
    EraCheck C;
    C.from = P.from + (P.to - P.from) * static_cast<double>(q) / static_cast<double>(eras);
    C.to = P.from + (P.to - P.from) * static_cast<double>(q + 1) / static_cast<double>(eras);
    std::vector<double> ek;
    double et = 0, ep = 0;
    integrate(C.from, C.to, kmax, ek, et, ep);
    for (const auto& R : P.rates) C.expected += R.rate * (R.k == 0 ? C.to - C.from : ek[R.k]);
    double var_open = 0;
    for (std::size_t v = 0; v < V; v++) {
      const double y = year[K.order[v]];
      const bool inside = (y >= C.from && y < C.to) || (q + 1 == eras && y == C.to);
      if (!inside) continue;
      C.occasions++;
      C.opened += opened[v];
      C.opened_expected += P.open.p * static_cast<double>(opp_open[v]);
      var_open += P.open.p * (1 - P.open.p) * static_cast<double>(opp_open[v]);
    }
    C.z = C.expected > 0 ? (static_cast<double>(C.occasions) - C.expected) / std::sqrt(C.expected) : 0;
    C.opened_z = var_open > 0 ? (static_cast<double>(C.opened) - C.opened_expected) / std::sqrt(var_open) : 0;
    P.chi2 += C.z * C.z;
    P.eras.push_back(C);
  }
  P.dof = P.eras.empty() ? 0 : P.eras.size() - 1;
  return P;
}

}  // namespace history
