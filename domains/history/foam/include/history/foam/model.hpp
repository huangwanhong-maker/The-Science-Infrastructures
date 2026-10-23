// The spin foam of a history at G = ℤ₂ (theory note §6; the idea note §3): the unified model of
// the social and street spin foams, with the frame group the record measures.
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
// Sampling is by single-site heat bath (Glauber dynamics), which leaves W invariant, and after each
// sweep a cluster move at every presence (b, o): all the rungs of b's ties at o flipped at once,
// accepted by Metropolis. Every triangle through the presence holds two of those rungs, so the
// move changes no triangle and only the lenses that begin or end there; where balance is strong
// and the triangles freeze, it carries the chain between configurations the single-site moves
// cannot reach in reasonable time. Its draws are reproducible, a calx Philox stream under a seed. Means of rungs and cells are measured by
// their conditional expectations at a free rung (Rao-Blackwell), unbiased and of far smaller
// variance than the raw values.
//
// The configuration side is a probability measure for every real β. The foam side is signed
// where some β < 0 (a colour whose ties flip more often than not, or triangles frustrated more
// often than not): tanh β < 0 makes some foams' amplitudes negative, the social spin foam's
// signed weights (its Rem. 7.5), and P(f ∈ S) is then a signed weight, not a probability.

#pragma once

#include <history/construction/complex.hpp>
#include <history/foam/refined.hpp>
#include <calx/random/stream.hpp>

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

// What an ensemble yields: the mean of each rung, of each cell's holonomy, and of any products
// of rungs asked for (Wilson loops in temporal gauge), with their standard errors from batches.
struct Ensemble {
  bool conditioned = false;
  std::vector<double> rung;    // ⟨u_r⟩
  std::vector<double> cell;    // ⟨h_f⟩
  std::vector<double> loop;    // ⟨Π_{r ∈ loop} u_r⟩
  std::vector<double> loop_error;
  std::vector<std::vector<std::size_t>> count_hist;  // for each counter, a histogram of its values
  std::size_t samples = 0;
};

// A count of frustrated triangles, evaluated on each draw (the network of the moment at a date).
struct Counter {
  std::vector<std::array<std::size_t, 3>> triangles;  // three rungs each
};

// Samples the prior (conditioned = false) or the posterior (the attested rungs fixed).
inline Ensemble sample(const HistoryModel& M, bool conditioned, const Sampling& S,
                       const std::vector<std::vector<std::size_t>>& loops = {}, const std::vector<Counter>& counters = {}) {
  calx::random::Stream z(S.seed, conditioned ? 1u : 0u, 0x4f4d4146u, 0u);
  const std::size_t n = M.rung.size();
  std::vector<int> u(n, 1);
  std::vector<bool> free(n, true);
  for (std::size_t r = 0; r < n; r++) {
    if (conditioned && M.attested[r]) { u[r] = M.attested[r]; free[r] = false; }
    else u[r] = z.uniform() < 0.5 ? 1 : -1;
  }
  // The local field at a rung: the weight is e^{field · u_r}.
  auto field_at = [&](std::size_t r) {
    double field = 0;
    for (std::size_t c : M.cells_of[r]) {
      const auto& cell = M.cells[c];
      int rest = 1;
      for (std::size_t q : cell.rungs)
        if (q != r) rest *= u[q];
      field += cell.beta * rest;
    }
    return field;
  };
  std::vector<int> flipped(M.cells.size(), 0);
  auto cluster = [&]() {
    for (const auto& g : M.presences) {
      if (std::any_of(g.begin(), g.end(), [&](std::size_t r) { return !free[r]; })) continue;
      std::vector<std::size_t> touched;
      for (std::size_t r : g)
        for (std::size_t c : M.cells_of[r]) {
          if (!flipped[c]) touched.push_back(c);
          flipped[c] ^= 1;
        }
      double delta = 0;  // the change of log W
      for (std::size_t c : touched)
        if (flipped[c]) delta -= 2 * M.cells[c].beta * holonomy(M.cells[c], u);
      for (std::size_t c : touched) flipped[c] = 0;
      if (delta >= 0 || z.uniform() < std::exp(delta))
        for (std::size_t r : g) u[r] = -u[r];
    }
  };
  auto sweep = [&]() {
    for (std::size_t r = 0; r < n; r++)
      if (free[r]) u[r] = z.uniform() < 1 / (1 + std::exp(-2 * field_at(r))) ? 1 : -1;
    cluster();
  };
  // Each cell is measured by its conditional mean at one free rung (Rao-Blackwell):
  // E[h_f | the rest] = (product of its other rungs) · tanh(field at that rung).
  std::vector<std::size_t> pivot(M.cells.size(), n);
  for (std::size_t c = 0; c < M.cells.size(); c++)
    for (std::size_t r : M.cells[c].rungs)
      if (free[r]) { pivot[c] = r; break; }
  for (std::size_t k = 0; k < S.burn; k++) sweep();
  Ensemble E;
  E.conditioned = conditioned;
  E.rung.assign(n, 0);
  E.cell.assign(M.cells.size(), 0);
  E.loop.assign(loops.size(), 0);
  E.loop_error.assign(loops.size(), 0);
  E.count_hist.assign(counters.size(), {});
  for (std::size_t c = 0; c < counters.size(); c++) E.count_hist[c].assign(counters[c].triangles.size() + 1, 0);
  const std::size_t batches = 20, per = std::max<std::size_t>(1, S.sweeps / batches);
  std::vector<std::vector<double>> batch(loops.size(), std::vector<double>(batches, 0));
  for (std::size_t k = 0; k < per * batches; k++) {
    sweep();
    for (std::size_t r = 0; r < n; r++) E.rung[r] += free[r] ? std::tanh(field_at(r)) : u[r];
    for (std::size_t c = 0; c < M.cells.size(); c++) {
      if (pivot[c] == n) { E.cell[c] += holonomy(M.cells[c], u); continue; }
      int rest = 1;
      for (std::size_t q : M.cells[c].rungs)
        if (q != pivot[c]) rest *= u[q];
      E.cell[c] += rest * std::tanh(field_at(pivot[c]));
    }
    for (std::size_t l = 0; l < loops.size(); l++) {
      int w = 1;
      for (std::size_t r : loops[l]) w *= u[r];
      E.loop[l] += w;
      batch[l][k / per] += w;
    }
    for (std::size_t c = 0; c < counters.size(); c++) {
      std::size_t f = 0;
      for (const auto& t : counters[c].triangles) f += u[t[0]] * u[t[1]] * u[t[2]] < 0;
      E.count_hist[c][f]++;
    }
  }
  E.samples = per * batches;
  const double N = static_cast<double>(E.samples);
  for (auto& x : E.rung) x /= N;
  for (auto& x : E.cell) x /= N;
  for (std::size_t l = 0; l < loops.size(); l++) {
    E.loop[l] /= N;
    double v = 0;
    for (double b : batch[l]) v += (b / static_cast<double>(per) - E.loop[l]) * (b / static_cast<double>(per) - E.loop[l]);
    E.loop_error[l] = std::sqrt(v / static_cast<double>(batches * (batches - 1)));
  }
  return E;
}

// The space of foams: the 2-cycles of K̂ over ℤ/2, as a basis found by Gaussian elimination of
// the boundary matrix over GF(2), column by column with the combinations tracked (the kernel of
// ∂₂: a column that reduces to zero gives a cycle). Its dimension is the number of independent
// foams, b₂(K̂; ℤ/2); a cell outside the support of every 2-cycle is in no foam.
struct FoamSpace {
  std::size_t dimension = 0;
  std::vector<std::vector<std::size_t>> basis;  // each a set of cells
  std::vector<bool> support;                    // cell -> in some foam
};

inline FoamSpace foam_space(const RefinedComplex& R) {
  const std::size_t E = R.cw.edges.size(), F = R.cw.faces.size();
  using Bits = std::vector<std::uint64_t>;
  auto set = [](Bits& b, std::size_t i) { b[i / 64] ^= std::uint64_t{1} << (i % 64); };
  auto low = [](const Bits& b) -> long long {
    for (std::size_t w = b.size(); w-- > 0;)
      if (b[w]) return static_cast<long long>(64 * w + 63 - static_cast<std::size_t>(__builtin_clzll(b[w])));
    return -1;
  };
  std::vector<Bits> column(F, Bits((E + 63) / 64, 0)), combination(F, Bits((F + 63) / 64, 0));
  std::vector<long long> owner(E, -1);
  FoamSpace S;
  S.support.assign(F, false);
  for (std::size_t f = 0; f < F; f++) {
    for (const auto& step : R.cw.faces[f].boundary) set(column[f], step.edge);
    set(combination[f], f);
    for (long long l = low(column[f]); l >= 0 && owner[l] >= 0; l = low(column[f])) {
      const auto& o = static_cast<std::size_t>(owner[l]);
      for (std::size_t w = 0; w < column[f].size(); w++) column[f][w] ^= column[o][w];
      for (std::size_t w = 0; w < combination[f].size(); w++) combination[f][w] ^= combination[o][w];
    }
    const long long l = low(column[f]);
    if (l >= 0) { owner[l] = static_cast<long long>(f); continue; }
    std::vector<std::size_t> cycle;
    for (std::size_t g = 0; g < F; g++)
      if ((combination[f][g / 64] >> (g % 64)) & 1) { cycle.push_back(g); S.support[g] = true; }
    S.basis.push_back(std::move(cycle));
  }
  S.dimension = S.basis.size();
  return S;
}

// The spin of each cell in the foam ensemble of the prior: P(f ∈ S) from ⟨h_f⟩ (the exact
// identity), and exactly 0 for a cell outside the support of the foams.
inline std::vector<double> spins(const HistoryModel& M, const Ensemble& prior, const FoamSpace* space = nullptr) {
  std::vector<double> j(M.cells.size(), 0);
  for (std::size_t c = 0; c < M.cells.size(); c++) {
    if (space && !space->support[c]) continue;
    const double t = std::tanh(M.cells[c].beta);
    j[c] = std::fabs(t) < 1e-12 ? 0 : t * (prior.cell[c] - t) / (1 - t * t);
  }
  return j;
}

// Exactly, by enumeration, for small models: ⟨u_r⟩ and ⟨h_f⟩ under the prior or the posterior.
// Refuses more than 24 free rungs.
inline Ensemble exact(const HistoryModel& M, bool conditioned) {
  const std::size_t n = M.rung.size();
  std::vector<std::size_t> vars;
  for (std::size_t r = 0; r < n; r++)
    if (!conditioned || !M.attested[r]) vars.push_back(r);
  if (vars.size() > 24) throw std::invalid_argument("too many free rungs to enumerate");
  std::vector<int> u(n, 1);
  for (std::size_t r = 0; r < n; r++)
    if (conditioned && M.attested[r]) u[r] = M.attested[r];
  Ensemble E;
  E.conditioned = conditioned;
  E.rung.assign(n, 0);
  E.cell.assign(M.cells.size(), 0);
  double Z = 0;
  for (std::uint64_t mask = 0; mask < (std::uint64_t{1} << vars.size()); mask++) {
    for (std::size_t i = 0; i < vars.size(); i++) u[vars[i]] = (mask >> i) & 1 ? -1 : 1;
    double logw = 0;
    for (const auto& c : M.cells) logw += c.beta * holonomy(c, u);
    const double w = std::exp(logw);
    Z += w;
    for (std::size_t r = 0; r < n; r++) E.rung[r] += w * u[r];
    for (std::size_t c = 0; c < M.cells.size(); c++) E.cell[c] += w * holonomy(M.cells[c], u);
  }
  for (auto& x : E.rung) x /= Z;
  for (auto& x : E.cell) x /= Z;
  return E;
}

// The spin foams themselves, enumerated for small models: every set of cells containing every
// rung an even number of times, with its amplitude Π tanh β_f; and from them P(f ∈ S).
// Refuses more than 24 cells.
inline std::vector<double> foam_spins_by_enumeration(const HistoryModel& M, std::size_t* foams = nullptr) {
  if (M.cells.size() > 24) throw std::invalid_argument("too many cells to enumerate foams");
  std::vector<double> in(M.cells.size(), 0);
  double total = 0;
  std::size_t count = 0;
  std::vector<int> parity(M.rung.size());
  for (std::uint64_t mask = 0; mask < (std::uint64_t{1} << M.cells.size()); mask++) {
    std::fill(parity.begin(), parity.end(), 0);
    double a = 1;
    for (std::size_t c = 0; c < M.cells.size(); c++)
      if ((mask >> c) & 1) {
        a *= std::tanh(M.cells[c].beta);
        for (std::size_t r : M.cells[c].rungs) parity[r] ^= 1;
      }
    if (std::any_of(parity.begin(), parity.end(), [](int p) { return p != 0; })) continue;
    count++;
    total += a;
    for (std::size_t c = 0; c < M.cells.size(); c++)
      if ((mask >> c) & 1) in[c] += a;
  }
  for (auto& x : in) x /= total;
  if (foams) *foams = count;
  return in;
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

// SU(2) kinematics: the dimension of the space of invariants of V_{j₁} ⊗ … ⊗ V_{jₙ}, the spins
// given as 2j; by Clebsch–Gordan, multiplicity by multiplicity (exact, in unsigned integers).
inline std::uint64_t invariant_dimension(const std::vector<int>& twice_j) {
  std::map<int, std::uint64_t> m{{0, 1}};  // 2J -> multiplicity
  for (int a : twice_j) {
    std::map<int, std::uint64_t> next;
    for (const auto& [b, k] : m)
      for (int c = std::abs(a - b); c <= a + b; c += 2) next[c] += k;
    m.swap(next);
  }
  const auto it = m.find(0);
  return it == m.end() ? 0 : it->second;
}

}  // namespace history
