// The spin foam of a history at G = ℤ₂ (theory note §6; the idea note §3): the unified model of
// the social and street spin foams, with the frame group the record measures, as a lattice gauge
// theory of physis on the history's refined complex.
//
// The variables are the rungs of the refined complex K̂ in temporal gauge, u_r = ±1: a tie's sign
// at a meeting. The weight of a configuration is a Gibbs weight, a product over the 2-cells:
//
//   W(u) = Π_lenses e^{β_c u u'} · Π_triangles e^{β_Δ u u' u''},
//
// the first a tie's persistence between consecutive meetings (β_c by the tie's colour), the
// second the balance of three ties at an occasion (β_Δ). W is a class function of holonomies,
// so it is unchanged by switching a being (a gauge transformation).
//
// The character expansion e^{βh} = cosh β (1 + h tanh β) turns the sum over configurations into a
// sum over spin foams (Reisenberger and Rovelli 1997; Baez 1998, §3; Oeckl 2005, ch. 5):
//
//   Z = Π_f cosh β_f · 2^{#rungs} · Σ_{S ∈ Z₂(K̂; ℤ/2)} Π_{f ∈ S} tanh β_f,
//
// where a spin foam S is a set of 2-cells containing every rung an even number of times, which,
// the temporal edges being a forest, is a 2-cycle of K̂ over ℤ/2. Read cell by cell:
//
//   spin j_f ∈ {0, 1}, the irreps of ℤ₂: whether the 2-cell belongs to the foam;
//   face amplitude c₁/c₀ = tanh β_f: how strongly a tie carries its sign across a lens, or a
//     triangle its balance;
//   intertwiners: dim Inv(⊗ V_{j_f}) over the cells at an edge is 1 when the spins there sum to
//     an even number and 0 otherwise: the Gauss law, so the social spin foam's admissibility
//     (its Prop. 6.4) is a theorem;
//   vertex amplitudes: 1 on admissible configurations, ℤ₂ being abelian;
//   the amplitude of a foam: Π_{f ∈ S} tanh β_f.
//
// From a configuration ensemble the spins follow exactly: d log Z / dβ_f gives
//
//   P(f ∈ S) = t_f (⟨h_f⟩ − t_f) / (1 − t_f²),  t_f = tanh β_f,
//
// so a 2-cell in no 2-cycle has ⟨h_f⟩ = t_f and never belongs to a foam.
//
// The prior is W alone: the theory, whose spins and Wilson loops describe the history's
// coherence. The posterior conditions W on the attested signs: the reconstruction ensemble
// (the social spin foam's Def. 8.5), in which each unattested sign is inferred from all the
// others, through persistence and balance at once, rather than carried forward.
//
// The theory is physis's ℤ₂ lattice gauge theory on the 2-complex K̂ (physis gauge), with a
// coupling for each 2-cell: computed exactly by physis's contraction (transfer.hpp), and sampled,
// where the contraction is out of reach, by physis's heat bath (sample(), below). The domain adds
// only what is the record's: which cells, which couplings, which links are fixed by the record.
//
// The configuration side is a probability measure for every real β. The foam side is signed
// where some β < 0 (a colour whose ties flip more often than not, or triangles frustrated more
// often than not): tanh β < 0 makes some foams' amplitudes negative, the social spin foam's
// signed weights (its Rem. 7.5), and P(f ∈ S) is then a signed weight, not a probability.

#pragma once

#include <history/construction/complex.hpp>
#include <history/foam/refined.hpp>
#include <calx/random/stream.hpp>
#include <morphe/topology/cycles.hpp>
#include <physis/field_theory/statistics.hpp>
#include <physis/gauge/field.hpp>
#include <physis/gauge/updates.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace history {

// An estimate of one coupling: n cells with both ends attested, k of them flipped (lenses) or
// frustrated (triangles); p = (k + ½)/(n + 1), the Krichevsky–Trofimov estimate, so that a
// coupling is finite where nothing has flipped; t = 1 − 2p = tanh β.
struct CouplingEstimate {
  std::string of;  // a colour, or "balance"
  std::size_t n = 0, k = 0;
  double p = 0.5, t = 0, beta = 0;
};

inline CouplingEstimate coupling_estimate(std::string of, std::size_t n, std::size_t k) {
  CouplingEstimate c;
  c.of = std::move(of);
  c.n = n;
  c.k = k;
  c.p = (static_cast<double>(k) + 0.5) / (static_cast<double>(n) + 1);
  c.t = 1 - 2 * c.p;
  c.beta = std::atanh(c.t);
  return c;
}

struct Couplings {
  std::map<std::string, double> lens;  // β by colour
  double balance = 0;                  // β_Δ
  std::vector<CouplingEstimate> estimates;
  double of(const std::string& colour) const {
    const auto it = lens.find(colour);
    return it == lens.end() ? lens.at("all") : it->second;
  }
};

// The couplings estimated from the record: each colour's lenses with both meetings attested
// (pooled into "all" where a colour has none), and the triangles at occasions whose three
// signs are all attested. Moment estimates, each as if its cells were independent.
inline Couplings estimate_couplings(const HistoryComplex& K, const RefinedComplex& R) {
  Couplings C;
  std::map<std::string, std::array<std::size_t, 2>> count;
  std::size_t n = 0, k = 0;
  for (std::size_t f = 0; f < K.lenses.size(); f++) {
    const Lens& L = K.lenses[f];
    const auto& g = K.rungs[L.tie];
    if (!g[L.rung_from].attested || !g[L.rung_to].attested) continue;
    auto& c = count[K.record.ties[L.tie].colour];
    c[0]++;
    c[1] += K.flip(f);
    n++;
    k += K.flip(f);
  }
  for (const auto& colour : colours()) {
    if (!count.count(colour)) continue;
    C.estimates.push_back(coupling_estimate(colour, count[colour][0], count[colour][1]));
    C.lens[colour] = C.estimates.back().beta;
  }
  C.estimates.push_back(coupling_estimate("all", n, k));
  C.lens["all"] = C.estimates.back().beta;
  std::size_t tn = 0, tk = 0;
  for (const auto& T : R.triangles) {
    bool attested = true;
    int product = 1;
    for (std::size_t q = 0; q < 3; q++) {
      const Rung& g = K.rungs[T.ties[q]][T.rungs[q]];
      attested = attested && g.attested;
      product *= g.sign;
    }
    if (!attested) continue;
    tn++;
    tk += product < 0;
  }
  C.estimates.push_back(coupling_estimate("balance", tn, tk));
  C.balance = C.estimates.back().beta;
  return C;
}

// The model on the rungs: variables, their attested values, and the 2-cells that couple them.
struct HistoryModel {
  std::vector<std::array<std::size_t, 2>> rung;  // r -> (tie, rung index)
  std::vector<std::vector<std::size_t>> index;   // tie -> rung index -> r
  std::vector<int> attested;                     // r -> ±1 attested, 0 not
  struct Cell { std::vector<std::size_t> rungs; double beta = 0; bool triangle = false; };
  std::vector<Cell> cells;                       // the lenses, in K's order, then the triangles
  std::vector<std::vector<std::size_t>> cells_of;  // r -> cells
  std::vector<std::vector<std::size_t>> presences; // the rungs at each presence (b, o): flipped together by the cluster move

  std::size_t lenses() const { return static_cast<std::size_t>(std::count_if(cells.begin(), cells.end(), [](const Cell& c) { return !c.triangle; })); }
};

inline HistoryModel history_model(const HistoryComplex& K, const RefinedComplex& R, const Couplings& C) {
  HistoryModel M;
  M.index.assign(K.rungs.size(), {});
  for (std::size_t t = 0; t < K.rungs.size(); t++)
    for (std::size_t k = 0; k < K.rungs[t].size(); k++) {
      M.index[t].push_back(M.rung.size());
      M.rung.push_back({t, k});
      M.attested.push_back(K.rungs[t][k].attested ? K.rungs[t][k].sign : 0);
    }
  for (const auto& L : K.lenses)
    M.cells.push_back({{M.index[L.tie][L.rung_from], M.index[L.tie][L.rung_to]}, C.of(K.record.ties[L.tie].colour), false});
  for (const auto& T : R.triangles)
    M.cells.push_back({{M.index[T.ties[0]][T.rungs[0]], M.index[T.ties[1]][T.rungs[1]], M.index[T.ties[2]][T.rungs[2]]}, C.balance, true});
  M.cells_of.assign(M.rung.size(), {});
  for (std::size_t c = 0; c < M.cells.size(); c++)
    for (std::size_t r : M.cells[c].rungs) M.cells_of[r].push_back(c);
  std::map<std::array<std::size_t, 2>, std::vector<std::size_t>> at;  // (being, occasion) -> rungs
  for (std::size_t r = 0; r < M.rung.size(); r++) {
    const auto [t, k] = M.rung[r];
    const std::size_t o = K.rungs[t][k].occasion;
    at[{K.record.ties[t].a, o}].push_back(r);
    at[{K.record.ties[t].b, o}].push_back(r);
  }
  for (auto& [key, rs] : at)
    if (rs.size() >= 2) M.presences.push_back(std::move(rs));
  return M;
}

inline int holonomy(const HistoryModel::Cell& c, const std::vector<int>& u) {
  int h = 1;
  for (std::size_t r : c.rungs) h *= u[r];
  return h;
}

struct Sampling {
  std::size_t burn = 500, sweeps = 4000;
  std::uint64_t seed = 20261024;
};

// A count of frustrated triangles, evaluated on each draw (the network of the moment at a date).
struct Counter {
  std::vector<std::array<std::size_t, 3>> triangles;  // three rungs each
};

// What a sample is asked to measure besides the means of rungs and cells: products of 2-cells'
// holonomies and Wilson loops along closed paths of K̂, which are gauge invariant and so read in
// either ensemble; products of rungs and counts of frustrated triangles, read in the posterior,
// whose temporal gauge gives a rung its sign.
struct Observables {
  std::vector<std::vector<std::size_t>> cells;
  std::vector<morphe::topology::Path> loops;
  std::vector<std::vector<std::size_t>> rungs;
  std::vector<Counter> counters;
};

// What an ensemble yields: the mean of each rung and of each cell's holonomy, the observables asked
// for with their standard errors (binned, physis field_theory), and a histogram of each counter.
struct Ensemble {
  bool conditioned = false;
  std::vector<double> rung;    // ⟨u_r⟩, in temporal gauge (posterior); 0 by symmetry in the prior
  std::vector<double> cell;    // ⟨h_f⟩
  std::vector<double> product, product_error;  // the cells' products asked
  std::vector<double> loop, loop_error;        // the Wilson loops asked
  std::vector<double> rung_product;            // the rungs' products asked
  std::vector<std::vector<std::size_t>> count_hist;
  std::size_t samples = 0;
};

// Samples the prior (conditioned = false) or the posterior by physis's ℤ₂ lattice gauge theory on
// K̂ (physis gauge, field.hpp and updates.hpp): a coupling for each 2-cell, β_f = the cell's, and
// the heat bath of a finite group link by link. The prior leaves every link free, so that the updates of the temporal
// edges move the field along its gauge orbits: flipping a being's edge between two occasions
// changes exactly the lenses along it, the move that keeps a strongly balanced field from freezing,
// and every observable read in it is gauge invariant. The posterior freezes the temporal edges at
// the identity (the temporal gauge, in which a rung is a tie's sign) and every attested rung at its
// sign. Draws are reproducible: a calx Philox stream under a seed.
inline Ensemble sample(const RefinedComplex& R, const HistoryModel& M, bool conditioned, const Sampling& S, const Observables& O = {}) {
  namespace pg = physis::gauge;
  calx::random::Stream z(S.seed, conditioned ? 1u : 0u, 0x4f4d4146u, 0u);
  pg::GaugeField<pg::ZN> U(pg::ZN{2}, R.cw);
  std::vector<double> w(M.cells.size());
  for (std::size_t c = 0; c < M.cells.size(); c++) w[c] = M.cells[c].beta;
  U.set_couplings(std::move(w));
  const std::size_t n = M.rung.size();
  std::vector<std::size_t> edge(n);
  for (std::size_t r = 0; r < n; r++) edge[r] = R.rung_edge[M.rung[r][0]][M.rung[r][1]];
  for (std::size_t e = 0; e < R.cw.edges.size(); e++) U.set(e, z.uniform() < 0.5 ? 0u : 1u);
  if (conditioned) {
    for (std::size_t e = 0; e < R.temporal; e++) U.freeze(e, 0u);
    for (std::size_t r = 0; r < n; r++)
      if (M.attested[r]) U.freeze(edge[r], M.attested[r] > 0 ? 0u : 1u);
  }
  auto sign = [&](std::size_t r) { return U.link(edge[r]) == 0u ? 1 : -1; };
  for (std::size_t k = 0; k < S.burn; k++) pg::heat_bath_sweep(U, 1.0, z);
  Ensemble E;
  E.conditioned = conditioned;
  E.rung.assign(n, 0);
  E.cell.assign(M.cells.size(), 0);
  E.rung_product.assign(O.rungs.size(), 0);
  E.count_hist.assign(O.counters.size(), {});
  for (std::size_t c = 0; c < O.counters.size(); c++) E.count_hist[c].assign(O.counters[c].triangles.size() + 1, 0);
  std::vector<std::vector<double>> products(O.cells.size()), loops(O.loops.size());
  for (std::size_t k = 0; k < S.sweeps; k++) {
    pg::heat_bath_sweep(U, 1.0, z);
    for (std::size_t r = 0; r < n; r++) E.rung[r] += sign(r);
    for (std::size_t c = 0; c < M.cells.size(); c++) E.cell[c] += U.plaquette(c);
    for (std::size_t i = 0; i < O.cells.size(); i++) {
      double h = 1;
      for (std::size_t c : O.cells[i]) h *= U.plaquette(c);
      products[i].push_back(h);
    }
    for (std::size_t i = 0; i < O.loops.size(); i++) loops[i].push_back(U.wilson_loop(O.loops[i]));
    for (std::size_t i = 0; i < O.rungs.size(); i++) {
      int h = 1;
      for (std::size_t r : O.rungs[i]) h *= sign(r);
      E.rung_product[i] += h;
    }
    for (std::size_t c = 0; c < O.counters.size(); c++) {
      std::size_t f = 0;
      for (const auto& t : O.counters[c].triangles) f += sign(t[0]) * sign(t[1]) * sign(t[2]) < 0;
      E.count_hist[c][f]++;
    }
  }
  E.samples = S.sweeps;
  const double N = static_cast<double>(std::max<std::size_t>(1, S.sweeps));
  for (auto& x : E.rung) x /= N;
  for (auto& x : E.cell) x /= N;
  for (auto& x : E.rung_product) x /= N;
  for (const auto& series : products) {
    const auto e = physis::field_theory::binned(series);
    E.product.push_back(e.mean);
    E.product_error.push_back(e.error);
  }
  for (const auto& series : loops) {
    const auto e = physis::field_theory::binned(series);
    E.loop.push_back(e.mean);
    E.loop_error.push_back(e.error);
  }
  return E;
}

// The space of foams: the 2-cycles of K̂ over ℤ/2 (morphe's cycles_mod2), its dimension the number
// of independent foams, b₂(K̂; ℤ/2); a cell outside the support of every 2-cycle is in no foam.
using FoamSpace = morphe::topology::Cycles;
inline FoamSpace foam_space(const RefinedComplex& R) { return morphe::topology::cycles_mod2(R.cw, 2); }

// The spin of each cell in the foam ensemble of the prior: P(f ∈ S) from ⟨h_f⟩ (the exact
// identity), and exactly 0 for a cell outside the support of the foams.
inline std::vector<double> spins(const HistoryModel& M, const std::vector<double>& prior_cell, const FoamSpace* space = nullptr) {
  std::vector<double> j(M.cells.size(), 0);
  for (std::size_t c = 0; c < M.cells.size(); c++) {
    if (space && !space->support[c]) continue;
    const double t = std::tanh(M.cells[c].beta);
    j[c] = std::fabs(t) < 1e-12 ? 0 : t * (prior_cell[c] - t) / (1 - t * t);
  }
  return j;
}

// The divergence of two ensembles' marginals on the given rungs: Σ_r D(Bern(p_r) ‖ Bern(q_r)),
// in nats. It approximates the divergence of the joint distributions (the social spin foam's
// power, Def. 10.3) when the rungs are weakly dependent.
inline double marginal_divergence(const std::vector<double>& p_mean, const std::vector<double>& q_mean, const std::vector<std::size_t>& rungs) {
  double d = 0;
  for (std::size_t r : rungs) {
    const double p = std::clamp((1 + p_mean[r]) / 2, 1e-6, 1 - 1e-6), q = std::clamp((1 + q_mean[r]) / 2, 1e-6, 1 - 1e-6);
    d += p * std::log(p / q) + (1 - p) * std::log((1 - p) / (1 - q));
  }
  return d;
}

}  // namespace history
