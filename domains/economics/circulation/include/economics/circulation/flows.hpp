// The circulation of flows on a transfer network (theory note §3).
//
// A flow F, the units surrendered per unit of time on each leg, is a 1-chain. Each leg also has a
// conductance c > 0, its capacity to carry flow (1 unless given; for trade, the gross volume of the
// link). physis's flow decomposition (Ohm and Kirchhoff on the network's complex) splits it into
// three flows, orthogonal in the energy Σ F² / c:
//
//   F = c d₀ α + c δ₁ γ + c h,
//
//   c d₀ α, the potential-driven flow: along each leg c (α(head) − α(tail)), so value runs up the
//     attraction α, from where it is drawn to where it is drawn to, each leg carrying it in
//     proportion to its conductance;
//   c δ₁ γ, the circulation round the network's faces (regions, triangles of places);
//   c h, the circulation round its holes: loops of legs no face fills, a round trip A → B → A
//     among them, of dimension b₁.
//
// The boundary of F, d₀ᵀ F, is the net inflow at each position; only the potential-driven part
// has one (the other two are divergence-free), so a flow is pure circulation exactly when every
// position is balanced. The pairing of a flow with the transports, Σ F log U, needs no metric, and
// with log U decomposed in the same one it splits into the three pairings of like parts: a
// potential-driven flow meets only the price system, and the loop deficits, which live in the
// co-exact and harmonic parts of log U, are met only by circulating flow. Value is extracted round
// loops only by what circulates, whatever the conductances.
//
// Balances by kind: the net inflow at each position of the flows of each kind, the stock–flow
// accounting of the user's typed lifecycle diagnostic (Huang 2026), with stocks changing by
// inflow − outflow per unit of time.

#pragma once

#include <economics/transfer/network.hpp>
#include <physis/discretisation/hodge.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <utility>
#include <string>
#include <vector>

namespace economics {

struct FlowCirculation {
  physis::discretisation::FlowDecomposition parts;  // the driven, local and global flows
  Cochain attraction;     // α, summing to zero on each component
  Cochain net_inflow;     // d₀ᵀ F at each position
  double total = 0;       // the energy Σ F² / c
  double driven = 0, local = 0, global = 0;  // the shares of the energy in the three parts
};

inline FlowCirculation circulation(const TransferNetwork& N, const Cochain& conductance = {}) {
  namespace pd = physis::discretisation;
  const CWComplex X = N.complex();
  const Cochain F = N.flows();
  FlowCirculation C;
  C.parts = pd::flow_decomposition(X, F, conductance);
  C.attraction = C.parts.potential;
  C.net_inflow = pd::d0_transpose(X, F);
  C.total = C.parts.energy;
  if (C.total > 0) {
    C.driven = C.parts.driven_energy / C.total;
    C.local = C.parts.local_energy / C.total;
    C.global = C.parts.global_energy / C.total;
  }
  return C;
}

// Σ F log U and its three terms: potential-driven flow against the price system, circulation
// round faces against the curvature, circulation round holes against the harmonic deficits.
struct Pairing {
  double total = 0, prices = 0, local = 0, global = 0;
};

inline Pairing pairing(const TransferNetwork& N, const Cochain& conductance = {}) {
  namespace pd = physis::discretisation;
  const CWComplex X = N.complex();
  const Cochain F = N.flows(), w = N.transport();
  const auto f = pd::flow_decomposition(X, F, conductance);
  const auto u = pd::hodge_decomposition(X, w, pd::Weights{conductance, {}});
  Pairing P;
  for (std::size_t e = 0; e < F.size(); e++) {
    P.total += F[e] * w[e];
    P.prices += f.driven[e] * u.exact[e];
    P.local += f.local[e] * u.coexact[e];
    P.global += f.global[e] * u.harmonic[e];
  }
  return P;
}

// Reciprocal exchange and net transfer. Flows both ways between two positions are, up to the
// smaller of the two, a round trip: pure circulation, exchange as such. What is left is one net
// flow for each pair of positions, from the one that sends more to the one that sends less, whose
// decomposition says how much is driven and how much circulates among several positions. The
// net network has one leg for each pair that exchanges anything (transport 1: flows only), along
// its net flow (0 where the two directions balance); `gross` holds the volume of each of its legs
// both ways, F_uv + F_vu, the natural conductance of the link; and `reciprocal` the volume
// exchanged both ways, 2 Σ min(F_uv, F_vu).
struct NetTransfer {
  TransferNetwork network;
  Cochain gross;
  double reciprocal = 0, net = 0, total = 0;
};

inline NetTransfer net_transfer(const TransferNetwork& N) {
  std::map<std::pair<std::size_t, std::size_t>, std::array<double, 2>> pair;  // (u < v) -> (u → v, v → u)
  NetTransfer T;
  for (const auto& l : N.legs) {
    if (l.from == l.to) continue;
    auto& f = pair[std::minmax(l.from, l.to)];
    (l.from < l.to ? f[0] : f[1]) += l.flow;
    T.total += l.flow;
  }
  T.network.positions = N.positions;
  T.network.faces.clear();
  for (const auto& [uv, f] : pair) {
    if (!(f[0] + f[1] > 0)) continue;
    T.reciprocal += 2 * std::min(f[0], f[1]);
    const double d = f[0] - f[1];
    if (d >= 0) T.network.legs.push_back({uv.first, uv.second, 1.0, d, "net"});
    else T.network.legs.push_back({uv.second, uv.first, 1.0, -d, "net"});
    T.gross.push_back(f[0] + f[1]);
    T.net += std::fabs(d);
  }
  return T;
}

// The net inflow at each position of the flows of each kind.
inline std::map<std::string, Cochain> balances(const TransferNetwork& N) {
  std::map<std::string, Cochain> B;
  for (const auto& l : N.legs) {
    auto& b = B[l.kind];
    if (b.empty()) b.assign(N.positions.size(), 0.0);
    b[l.to] += l.flow;
    b[l.from] -= l.flow;
  }
  return B;
}

}  // namespace economics
