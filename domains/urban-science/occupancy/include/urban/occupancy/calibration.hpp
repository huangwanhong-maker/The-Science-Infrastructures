// Calibrating the owners and occupiers by indirect inference (Gouriéroux, Monfort and Renault 1993, "Indirect
// inference", J. Appl. Econometrics 8, S85–S118; Smith 1993, ibid., S63–S84; the design's §2.4 and §10; theory/
// occupancy.md, §4).
//
// The auxiliary model is the one the module `events` fits to an observed record, the nested chain of the uses
// (physis's fit_rates), with the record's event rates beside it: the same statistics ψ are computed on the observed
// panel and on panels simulated from the model at θ, and θ is chosen so that the two agree. ψ has q = L + 8
// components, over one-year pairs of censuses unless said:
//   ψ_l       ln λ̂_l, the chain's rates fitted to every transition of the panel (l = 0 … L − 1);
//   ψ_L       falling vacant: the share of pairs (occupied, then vacant) among those occupied first;
//   ψ_{L+1}   let: the share of pairs (vacant, then occupied) among those vacant first;
//   ψ_{L+2}   falling into works: the share of pairs (occupied, then under works) among those occupied first;
//   ψ_{L+7}   moving between occupied families: the share of pairs (occupied, then occupied in another family) among
//             those occupied first, which the level-0 succession alone drives (it separates λ_0 from ν and η);
//   ψ_{L+3}   returning: of the vacancies seen to begin and end between occupied censuses, the share ending in the
//             family that preceded them;
//   ψ_{L+4}   the neighbours' pull: over the changes of family between consecutive censuses, the mean of the new use's
//             affinity to the neighbours seen at the earlier census, less its mean under π over the uses of the other
//             families (near zero under the chain alone);
//   ψ_{L+5}   the response of falling vacant to demand: the least-squares slope of each year's share falling vacant on
//             the year's demand d(y);
//   ψ_{L+6}   the response of letting to demand, likewise.
// θ = (ln λ_0, …, ln λ_{L−1}, ln ν, ln η, logit ρ, β_n, β_d) (model.hpp); a specification frees some of them and holds
// the rest (the chain alone: ν = η = 0, β = 0).
//
// The estimator. J(θ) = Σ_k (ψ_k^obs − ψ̄_k(θ))² / v_k, ψ̄(θ) the mean over S simulated panels of common random numbers
// (replicate r from the seed `seed + r` at every θ), v_k the variance of ψ_k over the replicates at the start. Its
// minimum in two stages: the binding function ψ̄ inverted by Levenberg and Marquardt's damped Gauss–Newton steps,
// δ = (GᵀWG + μ diag GᵀWG)⁻¹ GᵀW (ψ^obs − ψ̄(θ)), G = ∂ψ̄/∂θ by central differences on the common random numbers
// (steps of 0.08 in a log-rate, 0.15 in logit ρ or a β, far above J's roughness), a step halved until J falls (Marquardt 1963,
// J. SIAM 11, 431–441); then Nelder and Mead (calx) from that point, to polish. A simplex alone, in six dimensions on a
// piecewise constant J, stops short of the minimum (McKinnon 1998). Its covariance (Gouriéroux et al. 1993,
// Proposition 4, with a diagonal weight W):
//     Var θ̂ = (1 + 1/S) (GᵀWG)⁻¹ GᵀWΣWG (GᵀWG)⁻¹,
// G = ∂ψ̄/∂θ by central differences over 2S replicates (common random numbers) and Σ the covariance of ψ over
// S' = max(2S, 10q) replicates at θ̂ (the observed panel's, if the model is right; ten replicates a statistic, since the
// inverse of a covariance estimated from few is unstable). Over-identification (p parameters free): Hansen's statistic
// S'/(1 + S') gᵀ Σ⁻¹ g, g = ψ^obs − ψ̄(θ̂) over the S' replicates, referred to χ² on q − p degrees of freedom; exact
// for the optimal weight, a diagnostic otherwise. Between nested specifications the difference of the statistic at
// their two minima, with the larger model's Σ, is the distance-metric test.
//
// Caveat. The simulated outcomes are discrete, so J is piecewise constant in θ at fixed seeds: Nelder and Mead
// tolerates it, the replicates' number S sets the roughness, and G is taken by differences with steps far above it.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/optimisation/nelder_mead.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <urban/occupancy/model.hpp>

namespace urban::occupancy {

using Seen = std::vector<std::vector<std::pair<int, std::size_t>>>;  // per building: (census year, use), ascending

// A panel's transitions between consecutive observations of each building, aggregated by (x, y, years between): what
// urban::events::transitions_of gives from named observations (the check holds the two equal).
inline std::vector<sm::Transition> transitions(const Seen& seen) {
  std::map<std::tuple<std::size_t, std::size_t, int>, double> agg;
  for (const auto& s : seen)
    for (std::size_t k = 0; k + 1 < s.size(); k++)
      if (s[k + 1].first > s[k].first) agg[{s[k].second, s[k + 1].second, s[k + 1].first - s[k].first}] += 1;
  std::vector<sm::Transition> T;
  for (const auto& [key, n] : agg) T.push_back({std::get<0>(key), std::get<1>(key), static_cast<double>(std::get<2>(key)), n});
  return T;
}

struct Auxiliary {
  std::vector<double> psi;   // L + 8
  sm::NestedFit fit;
  double pairs_occupied = 0, pairs_vacant = 0, changes_of_family = 0, spells_closed = 0;
};

inline std::vector<std::string> auxiliary_names(std::size_t L) {
  std::vector<std::string> n;
  for (std::size_t l = 0; l < L; l++) n.push_back("ln λ̂_" + std::to_string(l));
  n.insert(n.end(), {"falling vacant", "let", "falling into works", "returning to the family", "neighbours' pull", "demand: falling vacant", "demand: let",
                     "moving between occupied families"});
  return n;
}

namespace detail {
inline double slope(const std::vector<double>& x, const std::vector<double>& y, const std::vector<char>& use) {
  double sx = 0, sy = 0, sxx = 0, sxy = 0, m = 0;
  for (std::size_t t = 0; t < x.size(); t++)
    if (use[t]) sx += x[t], sy += y[t], sxx += x[t] * x[t], sxy += x[t] * y[t], m += 1;
  const double den = m * sxx - sx * sx;
  return m > 2 && den > 0 ? (m * sxy - sx * sy) / den : 0.0;
}
}  // namespace detail

// The statistics ψ of a panel. `start`: the chain's rates to start its fit from (positive).
inline Auxiliary auxiliary(const Uses& S, const Space& X, const std::vector<double>& affinity, const Seen& seen, const std::vector<double>& demand, int from,
                           const std::vector<double>& start) {
  const std::size_t L = S.levels(), N = S.states(), n = seen.size();
  Auxiliary A;
  A.fit = sm::fit_rates(S.U.C, transitions(seen), start);
  for (std::size_t l = 0; l < L; l++) A.psi.push_back(std::log(std::max(A.fit.lambda[l], 1e-12)));
  // The panel by year.
  int to = from;
  for (const auto& s : seen)
    for (const auto& o : s) to = std::max(to, o.first);
  const std::size_t Y = static_cast<std::size_t>(to - from + 1);
  std::vector<std::vector<long>> at(n, std::vector<long>(Y, -1));
  for (std::size_t i = 0; i < n; i++)
    for (const auto& [y, u] : seen[i])
      if (y >= from) at[i][static_cast<std::size_t>(y - from)] = static_cast<long>(u);
  const auto vac = [&](long u) { return S.is_vacant(static_cast<std::size_t>(u)); };
  std::vector<double> occ(Y, 0.0), fell(Y, 0.0), vacant_first(Y, 0.0), let_by_year(Y, 0.0);
  double works = 0, moved = 0;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t t = 0; t + 1 < Y; t++) {
      const long a = at[i][t], b = at[i][t + 1];
      if (a < 0 || b < 0) continue;
      if (!vac(a)) {
        occ[t] += 1;
        if (vac(b)) fell[t] += 1;
        if (S.lump[static_cast<std::size_t>(b)] == Lump::works) works += 1;
        if (!vac(b) && S.family[static_cast<std::size_t>(a)] != S.family[static_cast<std::size_t>(b)]) moved += 1;
      } else {
        vacant_first[t] += 1;
        if (!vac(b)) let_by_year[t] += 1;
      }
    }
  double O = 0, F = 0, V = 0, Lt = 0;
  for (std::size_t t = 0; t < Y; t++) O += occ[t], F += fell[t], V += vacant_first[t], Lt += let_by_year[t];
  A.pairs_occupied = O, A.pairs_vacant = V;
  A.psi.push_back(O > 0 ? F / O : 0.0);
  A.psi.push_back(V > 0 ? Lt / V : 0.0);
  A.psi.push_back(O > 0 ? works / O : 0.0);
  // Returning to the family: vacancies between two occupied censuses, consecutive observations.
  double closed = 0, returned = 0;
  for (std::size_t i = 0; i < n; i++) {
    const auto& s = seen[i];
    for (std::size_t k = 0; k + 2 < s.size(); k++) {
      if (vac(static_cast<long>(s[k].second)) || !vac(static_cast<long>(s[k + 1].second))) continue;
      std::size_t e = k + 1;
      while (e < s.size() && vac(static_cast<long>(s[e].second))) e++;
      if (e == s.size()) continue;
      closed += 1, returned += S.family[s[e].second] == S.family[s[k].second];
    }
  }
  A.spells_closed = closed;
  A.psi.push_back(closed > 0 ? returned / closed : 0.0);
  // The neighbours' pull.
  double pull = 0, changes = 0;
  if (!affinity.empty())
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t k = 0; k + 1 < seen[i].size(); k++) {
        const auto [y1, x] = seen[i][k];
        const std::size_t z = seen[i][k + 1].second;
        if (S.family[x] == S.family[z] || y1 < from) continue;
        std::vector<std::size_t> nb;
        for (const std::size_t j : X.neighbours[i]) {
          const long u = at[j][static_cast<std::size_t>(y1 - from)];
          if (u >= 0) nb.push_back(static_cast<std::size_t>(u));
        }
        if (nb.empty()) continue;
        const auto mean_a = [&](std::size_t y) {
          double s = 0;
          for (const std::size_t u : nb) s += affinity[u * N + y];
          return s / static_cast<double>(nb.size());
        };
        double e = 0, w = 0;
        for (std::size_t y = 0; y < N; y++)
          if (S.family[y] != S.family[x]) e += S.U.C.law[y] * mean_a(y), w += S.U.C.law[y];
        pull += mean_a(z) - (w > 0 ? e / w : 0.0);
        changes += 1;
      }
  A.changes_of_family = changes;
  A.psi.push_back(changes > 0 ? pull / changes : 0.0);
  // The responses to demand: the slopes of the yearly shares on d(y).
  std::vector<double> d(Y, 0.0), rf(Y, 0.0), rl(Y, 0.0);
  std::vector<char> uf(Y, 0), ul(Y, 0);
  for (std::size_t t = 0; t + 1 < Y && t < demand.size(); t++) {
    d[t] = demand[t];
    if (occ[t] > 0) rf[t] = fell[t] / occ[t], uf[t] = 1;
    if (vacant_first[t] > 0) rl[t] = let_by_year[t] / vacant_first[t], ul[t] = 1;
  }
  A.psi.push_back(detail::slope(d, rf, uf));
  A.psi.push_back(detail::slope(d, rl, ul));
  A.psi.push_back(O > 0 ? moved / O : 0.0);
  return A;
}

// What a calibration needs of the record.
struct Calibrand {
  const Uses* uses = nullptr;
  const Space* space = nullptr;
  std::vector<double> affinity;
  std::vector<Life> lives;
  std::vector<double> demand;  // by year from `from`
  int from = 0, to = 0;
  double incentive = 1;
};

// θ ↔ parameters: (ln λ_0 … ln λ_{L−1}, ln ν, ln η, logit ρ, β_n, β_d); a rate whose logarithm is −∞ is zero.
inline Parameters parameters_of(const std::vector<double>& theta, std::size_t L, double incentive) {
  Parameters P;
  for (std::size_t l = 0; l < L; l++) P.lambda.push_back(std::exp(theta[l]));
  P.vacate = std::exp(theta[L]), P.let = std::exp(theta[L + 1]);
  P.returning = std::isfinite(theta[L + 2]) ? 1 / (1 + std::exp(-theta[L + 2])) : (theta[L + 2] > 0 ? 1.0 : 0.0);
  P.beta_neighbours = theta[L + 3], P.beta_demand = theta[L + 4], P.incentive = incentive;
  return P;
}

inline std::vector<double> theta_of(const Parameters& P) {
  std::vector<double> th;
  for (double l : P.lambda) th.push_back(std::log(l));
  constexpr double inf = std::numeric_limits<double>::infinity();
  th.push_back(P.vacate > 0 ? std::log(P.vacate) : -inf);
  th.push_back(P.let > 0 ? std::log(P.let) : -inf);
  th.push_back(P.returning <= 0 ? -inf : P.returning >= 1 ? inf : std::log(P.returning / (1 - P.returning)));
  th.push_back(P.beta_neighbours), th.push_back(P.beta_demand);
  return th;
}

// A start for the calibration from a panel's statistics, by the chain's one-year probabilities read backwards:
// - λ_0 from the share moving between occupied families, m ≈ (1 − e^{−λ_0}) c, c = Σ_F w_F (1 − π(F) − π(U)) the
//   chance a level-0 redraw from an occupied family F lands in another occupied family (w_F the occupied families'
//   shares, U the vacant family);
// - ν and η from what the chain leaves of the shares falling vacant and let, (1 − e^{−λ_0}) π(U) and
//   (1 − e^{−λ_0}) (1 − π(U)), floored at 0.002 and 0.01 a year;
// - ρ from the share returning to the family, r = ρ + (1 − ρ) p̄, p̄ = Σ_F w_F² the chance a fresh letting lands in the
//   family held, held to [0.05, 0.95];
// - λ_1 (and deeper levels) the chain's fit; β zero.
inline Parameters starting_parameters(const Uses& S, const Auxiliary& A, double incentive) {
  const std::size_t L = S.levels(), N = S.states();
  std::map<std::size_t, double> mass;
  double pu = 0, occupied_mass = 0;
  for (std::size_t y = 0; y < N; y++) {
    if (S.is_vacant(y)) pu += S.U.C.law[y];
    else mass[S.family[y]] += S.U.C.law[y], occupied_mass += S.U.C.law[y];
  }
  double c = 0, pbar = 0;
  for (const auto& [f, m] : mass) c += m / occupied_mass * (1 - m - pu), pbar += (m / occupied_mass) * (m / occupied_mass);
  const double moved = A.psi[L + 7];
  const double r0 = c > 0 && moved < c ? moved / c : 1 - std::exp(-A.fit.lambda[0]);
  Parameters P;
  P.lambda = A.fit.lambda;
  P.lambda[0] = -std::log(1 - r0);
  P.vacate = std::max(0.002, A.psi[L] - r0 * pu);
  P.let = std::max(0.01, A.psi[L + 1] - r0 * (1 - pu));
  P.returning = std::clamp((A.psi[L + 3] - pbar) / (1 - pbar), 0.05, 0.95);
  P.incentive = incentive;
  return P;
}

// ψ over S replicates at θ: their mean and covariance; replicate r from seed + r. On a pool the replicates run in
// parallel, each from its own seed, so the result does not depend on the threads.
struct Ensemble {
  std::vector<std::vector<double>> psi;  // by replicate
  std::vector<double> mean;
  std::vector<std::vector<double>> covariance;
};

inline Ensemble ensemble(const Calibrand& D, const std::vector<double>& theta, int replicates, std::uint64_t seed, const std::vector<double>& start,
                         calx::scheduler::ThreadPool* pool = nullptr) {
  const std::size_t L = D.uses->levels();
  Model M{D.uses, D.space, D.affinity, parameters_of(theta, L, D.incentive)};
  Ensemble E;
  E.psi.assign(static_cast<std::size_t>(replicates), {});
  const auto one = [&](std::size_t r) {
    const auto H = simulate(M, D.lives, D.demand, D.from, D.to, seed + r);
    E.psi[r] = auxiliary(*D.uses, *D.space, D.affinity, H.seen, D.demand, D.from, start).psi;
  };
  if (pool && replicates > 1)
    calx::scheduler::parallel_for(*pool, static_cast<std::size_t>(replicates), 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t r = b; r < e; r++) one(r);
    });
  else
    for (std::size_t r = 0; r < static_cast<std::size_t>(replicates); r++) one(r);
  const std::size_t q = E.psi[0].size();
  E.mean.assign(q, 0.0);
  for (const auto& p : E.psi)
    for (std::size_t k = 0; k < q; k++) E.mean[k] += p[k] / replicates;
  E.covariance.assign(q, std::vector<double>(q, 0.0));
  for (const auto& p : E.psi)
    for (std::size_t a = 0; a < q; a++)
      for (std::size_t b = 0; b < q; b++) E.covariance[a][b] += (p[a] - E.mean[a]) * (p[b] - E.mean[b]) / std::max(1, replicates - 1);
  return E;
}

struct CalibrationOptions {
  int replicates = 16;           // S
  std::uint64_t seed = 20261105;
  int max_evaluations = 300;     // of J, for Nelder and Mead's polish
  int max_iterations = 8;        // of Levenberg and Marquardt
  bool vacancy = true;           // ν, η and ρ free (else ν = η = 0, and demand acts on nothing)
  bool neighbours = true;        // β_n free (else held at the start's)
  bool demand = true;            // β_d free (else held at the start's)
  double step = 0.1;             // Nelder–Mead's initial step in each coordinate of θ (the polish)
  calx::scheduler::ThreadPool* pool = nullptr;
};

struct Calibration {
  std::vector<double> theta, error;           // θ̂ and its standard errors (0 where held)
  std::vector<bool> free;
  Parameters parameters;
  std::vector<double> psi_observed, psi_mean, psi_sd, weight;
  std::vector<std::vector<double>> covariance;  // of θ̂
  std::vector<std::vector<double>> psi_covariance;  // Σ at θ̂
  double objective = 0;                         // J(θ̂)
  double hansen = 0;                            // S'/(1+S') gᵀ Σ⁻¹ g
  int replicates_sigma = 0;                     // S'
  std::size_t degrees = 0;                      // q − p
  int evaluations = 0;
  bool converged = false;
  std::vector<std::pair<std::vector<double>, double>> trace;  // (θ, J) at every evaluation
};

namespace detail {
inline std::vector<std::vector<double>> inverse(const std::vector<std::vector<double>>& A) {
  const std::size_t n = A.size();
  calx::linalg::Matrix M(n), I(n);
  for (std::size_t i = 0; i < n; i++) {
    I(i, i) = 1;
    for (std::size_t j = 0; j < n; j++) M(i, j) = A[i][j];
  }
  const auto X = calx::linalg::solve(M, I);
  std::vector<std::vector<double>> out(n, std::vector<double>(n));
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) out[i][j] = X(i, j);
  return out;
}

// S'/(1 + S') gᵀ Σ⁻¹ g.
inline double hansen(const std::vector<double>& observed, const std::vector<double>& mean, const std::vector<std::vector<double>>& sigma, int replicates) {
  const auto Si = inverse(sigma);
  double h = 0;
  for (std::size_t a = 0; a < observed.size(); a++)
    for (std::size_t b = 0; b < observed.size(); b++) h += (observed[a] - mean[a]) * Si[a][b] * (observed[b] - mean[b]);
  return static_cast<double>(replicates) / (1.0 + replicates) * h;
}
}  // namespace detail

// Indirect inference from `start`: its rates usually the chain's fit to the observed panel and the channel of
// vacancy's rough rates (positive where free), its β zero.
inline Calibration calibrate(const Calibrand& D, const Seen& observed, const Parameters& start, const CalibrationOptions& o = {}) {
  const std::size_t L = D.uses->levels(), p = L + 5;
  auto theta0 = theta_of(start);
  Calibration C;
  C.free.assign(p, true);
  if (!o.vacancy) {
    C.free[L] = C.free[L + 1] = C.free[L + 2] = C.free[L + 4] = false;
    theta0[L] = theta0[L + 1] = theta0[L + 2] = -std::numeric_limits<double>::infinity(), theta0[L + 4] = 0;  // demand acts through vacancy alone
  }
  if (!o.neighbours) C.free[L + 3] = false;
  if (!o.demand) C.free[L + 4] = false;
  for (std::size_t j = 0; j < p; j++)
    if (C.free[j] && !std::isfinite(theta0[j])) throw std::invalid_argument("urban::occupancy::calibrate: a free parameter starts at a rate of zero");
  const std::vector<double> fit_start = start.lambda;
  C.psi_observed = auxiliary(*D.uses, *D.space, D.affinity, observed, D.demand, D.from, fit_start).psi;
  const std::size_t q = C.psi_observed.size();
  // The weight: the replicates' variances at the start.
  const auto E0 = ensemble(D, theta0, o.replicates, o.seed, fit_start, o.pool);
  C.weight.assign(q, 1.0);
  for (std::size_t k = 0; k < q; k++) C.weight[k] = E0.covariance[k][k] > 0 ? 1 / E0.covariance[k][k] : 1.0;
  std::vector<std::size_t> idx;
  for (std::size_t j = 0; j < p; j++)
    if (C.free[j]) idx.push_back(j);
  const auto full = [&](std::span<const double> x) {
    auto th = theta0;
    for (std::size_t k = 0; k < idx.size(); k++) th[idx[k]] = x[k];
    return th;
  };
  const auto J = [&](std::span<const double> x) {
    const auto th = full(x);
    const auto E = ensemble(D, th, o.replicates, o.seed, fit_start, o.pool);
    double j = 0;
    for (std::size_t k = 0; k < q; k++) j += C.weight[k] * (C.psi_observed[k] - E.mean[k]) * (C.psi_observed[k] - E.mean[k]);
    C.trace.push_back({th, j});
    return j;
  };
  const std::size_t pf = idx.size();
  const auto step_of = [&](std::size_t j) { return j >= L + 2 ? 0.15 : 0.08; };
  // The Jacobian of ψ̄ in the free coordinates at θ, by central differences.
  const auto jacobian = [&](const std::vector<double>& th, int replicates) {
    std::vector<std::vector<double>> G(q, std::vector<double>(pf, 0.0));
    for (std::size_t c = 0; c < pf; c++) {
      const std::size_t j = idx[c];
      const double h = step_of(j);
      auto up = th, dn = th;
      up[j] += h, dn[j] -= h;
      const auto Eu = ensemble(D, up, replicates, o.seed, fit_start, o.pool), Ed = ensemble(D, dn, replicates, o.seed, fit_start, o.pool);
      for (std::size_t k = 0; k < q; k++) G[k][c] = (Eu.mean[k] - Ed.mean[k]) / (2 * h);
    }
    return G;
  };
  std::vector<double> x;
  for (const std::size_t j : idx) x.push_back(theta0[j]);
  double jx = J(x);
  int evaluations = 1;
  bool settled = false;
  for (int it = 0; it < o.max_iterations && !settled; it++) {
    const auto th = full(x);
    const auto E = ensemble(D, th, o.replicates, o.seed, fit_start, o.pool);
    const auto G = jacobian(th, o.replicates);
    evaluations += 1 + 2 * static_cast<int>(pf);
    std::vector<std::vector<double>> A(pf, std::vector<double>(pf, 0.0));
    std::vector<double> b(pf, 0.0);
    for (std::size_t a = 0; a < pf; a++) {
      for (std::size_t k = 0; k < q; k++) b[a] += G[k][a] * C.weight[k] * (C.psi_observed[k] - E.mean[k]);
      for (std::size_t c = 0; c < pf; c++)
        for (std::size_t k = 0; k < q; k++) A[a][c] += G[k][a] * C.weight[k] * G[k][c];
    }
    bool moved = false;
    for (const double mu : {1e-3, 1e-1, 1.0, 10.0}) {
      auto M = A;
      for (std::size_t a = 0; a < pf; a++) M[a][a] *= 1 + mu;
      std::vector<double> delta(pf, 0.0);
      try {
        const auto Mi = detail::inverse(M);
        for (std::size_t a = 0; a < pf; a++)
          for (std::size_t c = 0; c < pf; c++) delta[a] += Mi[a][c] * b[c];
      } catch (const std::exception&) {
        continue;
      }
      double size = 0;
      for (double d : delta) size = std::max(size, std::fabs(d));
      const double cap = size > 1.0 ? 1.0 / size : 1.0;  // at most 1 in any coordinate a step
      for (int half = 0; half < 5 && !moved; half++) {
        std::vector<double> y = x;
        const double f = cap * std::pow(0.5, half);
        for (std::size_t a = 0; a < pf; a++) y[a] += f * delta[a];
        const double jy = J(y);
        evaluations++;
        if (jy < jx) {
          double change = 0;
          for (std::size_t a = 0; a < pf; a++) change = std::max(change, std::fabs(y[a] - x[a]));
          x = y, settled = change < 2e-3 || jx - jy < 1e-3, jx = jy, moved = true;
        }
      }
      if (moved) break;
    }
    if (!moved) settled = true;
  }
  calx::optimisation::NelderMeadOptions nm;
  nm.max_evaluations = o.max_evaluations, nm.initial_step = {o.step}, nm.x_tolerance = 1e-3, nm.f_tolerance = 1e-3;
  const auto R = calx::optimisation::nelder_mead(J, x, nm);
  evaluations += R.evaluations;
  if (R.value < jx) x = R.x, jx = R.value;
  C.theta = full(x), C.objective = jx, C.evaluations = evaluations, C.converged = settled || R.converged;
  C.parameters = parameters_of(C.theta, L, D.incentive);
  // At θ̂: Σ, G, the covariance of θ̂ and Hansen's statistic.
  const int S2 = 2 * o.replicates, Sp = std::max(S2, 10 * static_cast<int>(q));
  const auto Eh = ensemble(D, C.theta, Sp, o.seed, fit_start, o.pool);
  C.psi_mean = Eh.mean, C.psi_covariance = Eh.covariance;
  C.psi_sd.resize(q);
  for (std::size_t k = 0; k < q; k++) C.psi_sd[k] = std::sqrt(std::max(0.0, Eh.covariance[k][k]));
  const auto G = jacobian(C.theta, S2);
  std::vector<std::vector<double>> GWG(pf, std::vector<double>(pf, 0.0)), mid(pf, std::vector<double>(pf, 0.0));
  for (std::size_t a = 0; a < pf; a++)
    for (std::size_t b = 0; b < pf; b++)
      for (std::size_t k = 0; k < q; k++) {
        GWG[a][b] += G[k][a] * C.weight[k] * G[k][b];
        for (std::size_t l = 0; l < q; l++) mid[a][b] += G[k][a] * C.weight[k] * Eh.covariance[k][l] * C.weight[l] * G[l][b];
      }
  C.error.assign(p, 0.0);
  C.covariance.assign(p, std::vector<double>(p, 0.0));
  try {
    const auto Bi = detail::inverse(GWG);
    const double f = 1 + 1.0 / o.replicates;
    for (std::size_t a = 0; a < pf; a++)
      for (std::size_t b = 0; b < pf; b++) {
        double v = 0;
        for (std::size_t c = 0; c < pf; c++)
          for (std::size_t d = 0; d < pf; d++) v += Bi[a][c] * mid[c][d] * Bi[d][b];
        C.covariance[idx[a]][idx[b]] = f * v;
      }
    for (std::size_t j = 0; j < p; j++) C.error[j] = std::sqrt(std::max(0.0, C.covariance[j][j]));
  } catch (const std::exception&) {
    for (const std::size_t j : idx) C.error[j] = std::nan("");
  }
  try {
    C.hansen = detail::hansen(C.psi_observed, Eh.mean, Eh.covariance, Sp);
  } catch (const std::exception&) {
    C.hansen = std::nan("");
  }
  C.degrees = q - pf, C.replicates_sigma = Sp;
  return C;
}

}  // namespace urban::occupancy
