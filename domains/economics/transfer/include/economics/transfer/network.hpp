// Transfer networks: value carried between positions, and the invariants of carrying it round
// (theory note §2).
//
// A transfer network has positions and legs; a leg e = (u, v) carries a transport U_e > 0, the
// value registered at v for each unit surrendered at u, and may carry a flow, the units
// surrendered at u per unit of time. Re-describing the value at each position by a factor
// g_v > 0 changes every transport to g_v U_e g_u⁻¹ (Huang 2026, The Geometry of Value, §4): an
// abelian lattice gauge theory with group ℝ_{>0}, which the logarithm makes (ℝ, +) and physis's
// discrete exterior calculus computes (ω = log U, a 1-cochain on the network's complex).
//
//   holonomy of a closed walk γ:  λ(γ) = Π U_e^{±1} = exp ∫_γ ω, invariant under re-description;
//   loop deficit:                 δ(γ) = −log λ(γ), the share of value lost going round
//                                 (δ(γ⁻¹) = −δ(γ); deficits add under concatenation);
//   circuits:                     in the tree gauge of a spanning forest, each leg outside it (a
//                                 chord) carries the holonomy of its fundamental cycle, and these
//                                 are complete invariants: description classes ≅ ℝ_{>0}^{b₁};
//   price potential:              every deficit vanishes iff log U_e = φ(v) − φ(u) for a potential
//                                 φ, unique up to a constant on each component: a consistent
//                                 system of prices, and on a finite network the absence of
//                                 arbitrage round any loop (the gauge theory of exchange: Ilinski
//                                 2001; Young 1999).
//
// Retention (the paper's Lemma 8.4): on a directed cycle with steady flows, the position between
// two legs receives y = U x′ from the leg before and surrenders x on the leg after; the ratios
// r = y / x multiply to λ. Which position "absorbs" a deficit is not an invariant (any split of δ
// whose sum is δ is realised by some flows), and is not computed as if it were.

#pragma once

#include <physis/discretisation/cochains.hpp>
#include <morphe/combinatorics/graph.hpp>
#include <morphe/topology/cw_complex.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace economics {

using morphe::topology::CWComplex;
using morphe::topology::Path;
using physis::discretisation::Cochain;

struct Leg {
  std::size_t from = 0, to = 0;
  double ratio = 1;   // U > 0: value registered at `to` per unit surrendered at `from`
  double flow = 0;    // units surrendered at `from` per unit of time (0 if unknown)
  std::string kind;   // what is carried (money, goods, labour time, …), free text
};

struct TransferNetwork {
  std::vector<std::string> positions;
  std::vector<Leg> legs;
  std::vector<Path> faces;  // 2-cells, if any: closed walks of legs (regions, triangles of places)

  // The network as a CW complex: its positions, its legs (tail → head) and its faces, each based
  // at the tail of its first step.
  CWComplex complex() const {
    CWComplex X;
    X.vertices = positions.size();
    for (const auto& l : legs) {
      if (l.from >= positions.size() || l.to >= positions.size()) throw std::invalid_argument("economics::TransferNetwork: a leg leaves the positions");
      X.edges.push_back({l.from, l.to});
    }
    for (const auto& p : faces) {
      if (p.empty()) throw std::invalid_argument("economics::TransferNetwork: an empty face");
      const auto& first = legs.at(p.front().edge);
      X.faces.push_back({p.front().sign > 0 ? first.from : first.to, p});
    }
    X.validate();
    return X;
  }
  Cochain transport() const {
    Cochain w;
    for (const auto& l : legs) {
      if (!(l.ratio > 0)) throw std::invalid_argument("economics::TransferNetwork: a transport must be positive");
      w.push_back(std::log(l.ratio));
    }
    return w;
  }
  Cochain flows() const {
    Cochain f;
    for (const auto& l : legs) f.push_back(l.flow);
    return f;
  }
};

inline double holonomy(const TransferNetwork& N, const Path& loop) { return std::exp(physis::discretisation::integral(N.transport(), loop)); }
inline double deficit(const TransferNetwork& N, const Path& loop) { return -physis::discretisation::integral(N.transport(), loop); }

// The network re-described by g_v > 0 at each position: U_e ↦ g_head U_e / g_tail.
inline TransferNetwork regauged(TransferNetwork N, const std::vector<double>& g) {
  if (g.size() != N.positions.size()) throw std::invalid_argument("economics::regauged: one factor per position");
  for (auto& l : N.legs) l.ratio = g[l.to] * l.ratio / g[l.from];
  return N;
}

// The circuits of a network: a spanning forest's chords with their fundamental cycles (the chord
// forwards, then the forest back), each cycle's holonomy and deficit, and the potential of the
// forest (log prices, 0 at each root); `defect` is the largest |deficit|, zero iff there is no
// arbitrage.
struct Circuits {
  std::vector<std::size_t> chords;
  std::vector<Path> cycles;
  std::vector<double> holonomy, deficit;
  Cochain log_price;
  double defect = 0;
};

inline Circuits circuits(const TransferNetwork& N) {
  const CWComplex X = N.complex();
  const Cochain w = N.transport();
  const auto P = physis::discretisation::potential(X, w);
  const auto G = X.one_skeleton();
  const auto F = morphe::combinatorics::spanning_forest(G);
  Circuits C;
  C.log_price = P.phi;
  for (std::size_t e = 0; e < N.legs.size(); e++) {
    if (P.in_tree[e]) continue;
    Path cycle{{e, 1}};
    for (const auto& s : F.path(N.legs[e].to, N.legs[e].from, G)) cycle.push_back({s.edge, s.sign});
    const double integral = physis::discretisation::integral(w, cycle);
    C.chords.push_back(e);
    C.cycles.push_back(std::move(cycle));
    C.holonomy.push_back(std::exp(integral));
    C.deficit.push_back(-integral);
    C.defect = std::max(C.defect, std::fabs(integral));
  }
  return C;
}

// Retention along a directed cycle of legs (every step forwards, each leg's head the next leg's
// tail): for the position entered by leg i − 1 and left by leg i, r = U_{i−1} x_{i−1} / x_i.
// Π r = λ(γ) whenever the flows are positive.
inline std::vector<double> retention(const TransferNetwork& N, const Path& cycle) {
  const std::size_t n = cycle.size();
  std::vector<double> r(n);
  for (std::size_t i = 0; i < n; i++) {
    const auto& step = cycle[i];
    const auto& prev = cycle[(i + n - 1) % n];
    if (step.sign < 0 || prev.sign < 0) throw std::invalid_argument("economics::retention: a directed cycle, every leg forwards");
    const Leg& in = N.legs.at(prev.edge);
    const Leg& out = N.legs.at(step.edge);
    if (in.to != out.from) throw std::invalid_argument("economics::retention: the legs do not join");
    if (!(out.flow > 0) || !(in.flow > 0)) throw std::invalid_argument("economics::retention: positive flows on the cycle");
    r[i] = in.ratio * in.flow / out.flow;
  }
  return r;
}

}  // namespace economics
