// The landscape running (theory/landscape-dynamics.md §5 and §8): a district, two stores in it, the walkers on its
// streets, the warmth of its blocks, a permit, a lane and its recognition, as one hybrid process on the strata of
// its configuration — a piecewise-deterministic Markov process (physis dynamics, piecewise_deterministic.hpp) on
// each stratum, the complex changing at the events that rewrite it — with its record, anchored and replayable,
// the causal order of its events, and the balances it keeps.
//
// The augmented state z = (K, x, γ, h). K the district's complex (district.hpp); x the walkers at each junction,
// the temperature of each block, whether each store is active, the permit p and the recognition d of the lane;
// γ the event kinds enabled, which the state decides (the lane can be built only under a permit, removed only
// while not recognised); h the record.
//
// Between events (§5.1), on the stratum of K: walkers hop along the streets toward open stores and leave. A
// walker at junction v crosses segment e to w at the rate (G_e/c) e^{(u_w − u_v)/2}, G_e the segment's
// conductance, c a junction's capacity, u the attraction at the ends of an open store's segment; so
// ṅ = Σ_e (hops in − hops out) + e(t) − n/τ, e(t) the walkers entering at the district's gates in the hour of the
// day. Without attraction this is −Δ₀(n/c), the Hodge Laplacian with stars c on the junctions and G on the
// segments (physis hodge_laplacian.hpp); with it, the hops are a master equation in detailed balance with the
// Boltzmann law n ∝ c e^u (physis statistical_mechanics), and n stays non-negative, as a potential-driven flux
// G(n_v/c − u_v − n_w/c + u_w) would not (it can draw walkers out of an empty junction). Each block warms toward
// the air and the Sun,
// Ṫ_f = (T_air + s · sin⁺ altitude(t) + h n_f − T_f)/τ_f, with a little heat from the walkers about it. The
// continuous part keeps its walkers exactly: d/dt Σ n = Σ e − Σ n/τ, the streets carrying them without loss.
//
// Events (§5.2), each a mechanism with its matches and rates on the state: a store opens at
// α(1 + η·[the other open] + ξ·[the lane there])·e^{−κ (T − T*)⁺} during its hours and closes at β; the permit is
// granted at κ_p; under it the lane is built at κ₊ (three rewrites of K: district.hpp); while unrecognised it is
// removed at κ₋; recognition comes at κ_d times the use of the lane, the walkers at its two ends against a
// reference — an interpretation becoming a state that changes what may happen next (§10.4). At the closing hour
// the stores close: a forced jump, when the clock (in the state) reaches the boundary of the stratum's domain —
// a boundary condition of the generator, not a rate (Davis 1984).
//
// The record (§10.1): the initial state, then each event's time, mechanism, match, what it read and wrote, and
// for a lane its three rewrites; the state between is the flow's, so the record with the flow determines the
// path. Its causal order (morphe traces.hpp) puts each event after those whose writes it read: a store opens after
// the lane it answered to was built; the permit before the lane. Checks: tests/check_landscape.cpp.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <calx/random/stream.hpp>
#include <morphe/formal/traces.hpp>
#include <physis/dynamics/piecewise_deterministic.hpp>

#include "district.hpp"

namespace urban::landscape {

struct LandscapeParams {
  double alpha = 1.6, beta = 0.7, eta = 0.5, xi = 0.8;  // stores, per hour
  double open_from = 7, open_to = 22;                     // the stores' hours
  double heat_threshold = 30, heat_penalty = 0.15;        // °C, per °C above it
  double permit_rate = 0.3, build_rate = 0.5, remove_rate = 0.06, recognise_rate = 0.4;  // per hour
  double recognise_reference = 30;                        // walkers at the lane's ends
  double conductance = 200;                               // persons/h per unit of potential, for a segment of 100 m
  double capacity = 25;                                   // persons at a junction per unit of potential
  double attraction = 1.2;                                // at the ends of an open store's segment: e^u times the walkers
  double stay = 0.6;                                      // hours a walker stays, on average
  std::vector<double> entries;                            // persons/h entering at the gates in each hour (24)
  double air = 16, sun_gain = 9, relax = 1.5, walker_heat = 0.01;  // °C, °C, h, °C a walker
  std::vector<double> altitude;                           // the Sun's altitude in each hour (radians, 24)
  double hours = 24, frame = 0.25, step = 0.005;          // h; the walkers' flow is stiff, of rates up to ~200/h
  std::uint64_t seed = 1;
};

struct LandscapeEvent {
  double time = 0;
  std::string mechanism, match;
  bool forced = false;
  morphe::formal::traces::Footprint footprint;
  std::vector<Step> rewrites;  // for a lane, its rewrites
};

struct Frame {
  double time = 0;
  std::map<std::uint64_t, double> walkers, flow, temperature;  // by junction, segment (net, tail → head), block
  std::map<std::uint64_t, double> forth, back;                 // each segment's hops tail → head and head → tail, per hour
  bool a = false, b = false, permit = false, recognised = false;
  std::uint64_t lane = 0;
  std::size_t stratum = 0;
  std::vector<std::string> enabled;  // the event kinds that can happen now (the registry the state enables)
};

struct Stratum {
  double from = 0, to = 0;
  District district;  // the configuration over the interval (its state at `from`)
  std::uint64_t lane = 0;
  std::vector<std::size_t> betti;
};

struct LandscapeRun {
  std::vector<LandscapeEvent> events;
  std::vector<Frame> frames;
  std::vector<Stratum> strata;
  double entered = 0, left = 0, walkers0 = 0, walkers1 = 0;  // the balance of walkers
  std::string lane_refused;
};

namespace detail {

inline double hourly(const std::vector<double>& v, double t, double fallback) {
  if (v.empty()) return fallback;
  const double h = std::fmod(std::fmod(t, 24.0) + 24.0, 24.0);
  const std::size_t i = static_cast<std::size_t>(h) % v.size(), j = (i + 1) % v.size();
  const double f = h - std::floor(h);
  return v[i] * (1 - f) + v[j] * f;
}

// The stratum's fixed data: junctions and blocks in vector order, segments with their conductances.
struct Layout {
  std::vector<std::uint64_t> junctions, blocks, segments;
  std::map<std::uint64_t, std::size_t> jrow, brow;
  std::vector<std::array<std::size_t, 2>> ends;  // segment → (tail row, head row)
  std::vector<double> G;                         // segment conductances
  std::vector<std::size_t> gates;                // junction rows where walkers enter
  std::vector<std::vector<std::size_t>> block_junctions;  // block row → junction rows on it
  std::vector<double> share;                     // junction row → 1/(number of blocks at it)
};

inline Layout layout(const District& D, const LandscapeParams& P) {
  Layout L;
  for (const auto& c : D.complex.cells) {
    if (c.dim == 0) L.jrow[c.id] = L.junctions.size(), L.junctions.push_back(c.id);
    if (c.dim == 2) L.brow[c.id] = L.blocks.size(), L.blocks.push_back(c.id);
  }
  const auto up = D.complex.cofaces();
  std::vector<std::size_t> blocks_at(L.junctions.size(), 0);
  std::vector<char> gate(L.junctions.size(), 0);
  for (std::size_t i = 0; i < D.complex.size(); i++) {
    const auto& c = D.complex.cells[i];
    if (c.dim != 1) continue;
    const auto e = D.complex.ends(i);
    const std::size_t t = L.jrow.at(D.complex.cells[e[0]].id), h = L.jrow.at(D.complex.cells[e[1]].id);
    L.segments.push_back(c.id);
    L.ends.push_back({t, h});
    L.G.push_back(P.conductance * 100.0 / std::max(30.0, D.length(c.id)));
    if (up[i].size() < 2) gate[t] = gate[h] = 1;  // on the district's edge
  }
  for (std::size_t r = 0; r < L.junctions.size(); r++)
    if (gate[r]) L.gates.push_back(r);
  L.block_junctions.resize(L.blocks.size());
  for (std::size_t b = 0; b < L.blocks.size(); b++) {
    std::set<std::size_t> js;
    for (const auto& [e, s] : D.walk(L.blocks[b])) {
      const auto x = D.complex.ends(D.at(e));
      js.insert(L.jrow.at(D.complex.cells[x[0]].id)), js.insert(L.jrow.at(D.complex.cells[x[1]].id));
    }
    L.block_junctions[b].assign(js.begin(), js.end());
    for (auto j : js) blocks_at[j]++;
  }
  L.share.resize(L.junctions.size());
  for (std::size_t r = 0; r < L.junctions.size(); r++) L.share[r] = 1.0 / std::max<std::size_t>(1, blocks_at[r]);
  return L;
}

}  // namespace detail

// Runs the landscape over P.hours from t = 0 (midnight), with stores 0 and 1 of D the two stores.
inline LandscapeRun run_landscape(District D, const LandscapeParams& P) {
  if (D.stores.size() < 2) throw std::invalid_argument("urban::landscape::run_landscape: two stores");
  LandscapeRun out;
  calx::random::Stream rng(P.seed, 0x1a, 0x2b, 0x3c);
  bool a = false, b = false, permit = false, recognised = false;
  std::uint64_t lane = 0;
  int day = 0;
  // The lane's plan, made again whenever the cells change: closing a lane merges the streets it cut into new
  // segments, and a plan names the cells of the configuration it was made in.
  LanePlan plan = plan_lane(D, 0, 1);
  out.lane_refused = plan.refused;
  for (const auto id : D.of(0)) out.walkers0 += D.value(id, "walkers");
  double t = 0;
  // Footprint items: "a:A", "a:B", "permit", "lane", "recognised", "day", and the rewrites' own.
  const auto fp = [](std::set<std::string> reads, std::set<std::string> writes) { return morphe::formal::traces::Footprint{std::move(reads), std::move(writes)}; };
  std::map<std::uint64_t, double> flow_now;
  while (t < P.hours) {
    const auto L = detail::layout(D, P);
    const std::size_t nj = L.junctions.size(), nb = L.blocks.size();
    // The attraction at each junction from the open stores.
    std::vector<double> u(nj, 0.0);
    const auto attract = [&](const Store& s) {
      const auto x = D.complex.ends(D.at(s.segment));
      for (auto c : x) u[L.jrow.at(D.complex.cells[c].id)] += P.attraction;
    };
    if (a) attract(D.stores[0]);
    if (b) attract(D.stores[1]);
    // The blocks a store's segment bounds.
    const auto store_blocks = [&](const Store& s) {
      std::vector<std::size_t> rows;
      const auto up = D.complex.cofaces()[D.at(s.segment)];
      for (std::size_t c : up) rows.push_back(L.brow.at(D.complex.cells[c].id));
      return rows;
    };
    const auto ba = store_blocks(D.stores[0]), bb = store_blocks(D.stores[1]);
    std::vector<std::size_t> lane_rows;
    if (lane) {
      const auto x = D.complex.ends(D.at(lane));
      lane_rows = {L.jrow.at(D.complex.cells[x[0]].id), L.jrow.at(D.complex.cells[x[1]].id)};
    }
    // x = [t, n (nj), T (nb), entered, left, the next closing]: the closing hour in the state, so that the forced
    // jump at it is a boundary of the domain and its reset (the next day's closing) lands inside again.
    std::vector<double> x0(1 + nj + nb + 3);
    x0[0] = t;
    for (std::size_t r = 0; r < nj; r++) x0[1 + r] = D.value(L.junctions[r], "walkers");
    for (std::size_t q = 0; q < nb; q++) x0[1 + nj + q] = D.value(L.blocks[q], "temperature");
    x0[1 + nj + nb] = out.entered, x0[2 + nj + nb] = out.left, x0[3 + nj + nb] = P.open_to + 24.0 * day;
    const auto flows = [&](const std::vector<double>& x, std::vector<double>& j) {
      j.assign(L.segments.size(), 0.0);
      for (std::size_t e = 0; e < L.segments.size(); e++) {
        const auto [ta, he] = L.ends[e];
        const double k = L.G[e] / P.capacity, tilt = std::exp(0.5 * (u[he] - u[ta]));
        j[e] = k * (tilt * x[1 + ta] - x[1 + he] / tilt);  // net hops tail → head
      }
    };
    // The hops each way, whose difference is the net flow: the walkers a counter on the segment would see.
    const auto each_way = [&](const std::vector<double>& x, std::vector<double>& fw, std::vector<double>& bw) {
      fw.assign(L.segments.size(), 0.0), bw.assign(L.segments.size(), 0.0);
      for (std::size_t e = 0; e < L.segments.size(); e++) {
        const auto [ta, he] = L.ends[e];
        const double k = L.G[e] / P.capacity, tilt = std::exp(0.5 * (u[he] - u[ta]));
        fw[e] = k * tilt * x[1 + ta], bw[e] = k * x[1 + he] / tilt;
      }
    };
    physis::dynamics::PiecewiseDeterministic Pd;
    Pd.flow = [&](const std::vector<double>& x, std::vector<double>& dx) {
      dx.assign(x.size(), 0.0);
      dx[0] = 1;
      std::vector<double> j;
      flows(x, j);
      for (std::size_t e = 0; e < L.segments.size(); e++) dx[1 + L.ends[e][0]] -= j[e], dx[1 + L.ends[e][1]] += j[e];
      const double in = detail::hourly(P.entries, x[0], 600.0), per_gate = L.gates.empty() ? 0 : in / L.gates.size();
      double left = 0;
      for (std::size_t r = 0; r < nj; r++) {
        const double out_r = x[1 + r] / P.stay;
        dx[1 + r] += -out_r;
        left += out_r;
      }
      for (auto g : L.gates) dx[1 + g] += per_gate;
      const double alt = detail::hourly(P.altitude, x[0], 0.0);
      for (std::size_t q = 0; q < nb; q++) {
        double near = 0;
        for (auto r : L.block_junctions[q]) near += x[1 + r] * L.share[r];
        dx[1 + nj + q] = (P.air + P.sun_gain * std::max(0.0, std::sin(alt)) + P.walker_heat * near - x[1 + nj + q]) / P.relax;
      }
      dx[1 + nj + nb] = L.gates.empty() ? 0 : in;
      dx[2 + nj + nb] = left;
    };
    // The occurrences at a state x.
    struct Occ { std::string mechanism, match; double rate; int kind; };
    const auto mean_temp = [&](const std::vector<double>& x, const std::vector<std::size_t>& rows) {
      double s = 0;
      for (auto q : rows) s += x[1 + nj + q];
      return rows.empty() ? P.air : s / rows.size();
    };
    const auto occurrences = [&](const std::vector<double>& x) {
      std::vector<Occ> o;
      const double hour = std::fmod(x[0], 24.0);
      const bool open_hours = hour >= P.open_from && hour < P.open_to;
      const auto open_rate = [&](bool other, const std::vector<std::size_t>& rows) {
        return P.alpha * (1 + P.eta * other + P.xi * (lane != 0)) * std::exp(-P.heat_penalty * std::max(0.0, mean_temp(x, rows) - P.heat_threshold));
      };
      if (!a && open_hours) o.push_back({"a store opens", D.stores[0].name, open_rate(b, ba), 0});
      if (a) o.push_back({"a store closes", D.stores[0].name, P.beta, 1});
      if (!b && open_hours) o.push_back({"a store opens", D.stores[1].name, open_rate(a, bb), 2});
      if (b) o.push_back({"a store closes", D.stores[1].name, P.beta, 3});
      if (!permit) o.push_back({"the permit is granted", "the lane between the stores", P.permit_rate, 4});
      if (permit && !lane && plan.refused.empty()) o.push_back({"the lane is built", "through the block between the stores", P.build_rate, 5});
      if (lane && !recognised) o.push_back({"the lane is removed", "the lane", P.remove_rate, 6});
      if (lane && !recognised) {
        double use = 0;
        for (auto r : lane_rows) use += x[1 + r];
        o.push_back({"the lane is recognised", "the lane", P.recognise_rate * std::min(1.0, use / P.recognise_reference), 7});
      }
      return o;
    };
    Pd.rate = [&](const std::vector<double>& x) {
      double s = 0;
      for (const auto& q : occurrences(x)) s += q.rate;
      return s;
    };
    Pd.jump = [](const std::vector<double>& x, calx::random::Stream&) { return x; };
    const std::size_t next_closing = 3 + nj + nb;
    Pd.inside = [next_closing](const std::vector<double>& x) { return x[next_closing] - x[0]; };
    Pd.reset = [next_closing](const std::vector<double>& x, calx::random::Stream&) { auto y = x; y[next_closing] += 24; return y; };
    // The stratum's samples.
    const auto path = physis::dynamics::simulate(Pd, x0, P.hours - t, P.step, rng, 1, 1);
    Stratum S;
    S.from = t, S.district = D, S.lane = lane;
    S.betti = morphe::algebra::betti(mt::chain_complex(D.complex));
    const std::size_t stratum = out.strata.size();
    // Frames on this stratum: at each multiple of P.frame in [t, end), by linear interpolation of the samples.
    const double end = path.events.empty() ? P.hours : path.events[0].time + t;
    for (double f = std::ceil(t / P.frame - 1e-9) * P.frame; f < end - 1e-9 && f <= P.hours + 1e-9; f += P.frame) {
      std::size_t k = 0;
      while (k + 1 < path.times.size() && path.times[k + 1] + t < f) k++;
      std::vector<double> x = path.states[k];
      if (k + 1 < path.times.size()) {
        const double t0 = path.times[k] + t, t1 = path.times[k + 1] + t, w = t1 > t0 ? (f - t0) / (t1 - t0) : 0;
        for (std::size_t q = 0; q < x.size(); q++) x[q] += w * (path.states[k + 1][q] - x[q]);
      }
      Frame F;
      F.time = f, F.a = a, F.b = b, F.permit = permit, F.recognised = recognised, F.lane = lane, F.stratum = stratum;
      for (std::size_t r = 0; r < nj; r++) F.walkers[L.junctions[r]] = x[1 + r];
      for (std::size_t q = 0; q < nb; q++) F.temperature[L.blocks[q]] = x[1 + nj + q];
      std::vector<double> j;
      flows(x, j);
      for (std::size_t e = 0; e < L.segments.size(); e++) F.flow[L.segments[e]] = j[e];
      std::vector<double> fw, bw;
      each_way(x, fw, bw);
      for (std::size_t e = 0; e < L.segments.size(); e++) F.forth[L.segments[e]] = fw[e], F.back[L.segments[e]] = bw[e];
      std::set<std::string> kinds;
      for (const auto& o : occurrences(x)) kinds.insert(o.mechanism);
      F.enabled.assign(kinds.begin(), kinds.end());
      out.frames.push_back(std::move(F));
    }
    // The state at the end of the stratum, written back to D.
    const std::vector<double>& xe = path.events.empty() ? path.states.back() : path.events[0].before;
    for (std::size_t r = 0; r < nj; r++) D.set(L.junctions[r], "walkers", xe[1 + r]);
    for (std::size_t q = 0; q < nb; q++) D.set(L.blocks[q], "temperature", xe[1 + nj + q]);
    {
      std::vector<double> j;
      flows(xe, j);
      for (std::size_t e = 0; e < L.segments.size(); e++) D.set(L.segments[e], "flow", j[e]);
    }
    out.entered = xe[1 + nj + nb], out.left = xe[2 + nj + nb];
    S.to = path.events.empty() ? P.hours : xe[0];
    out.strata.push_back(std::move(S));
    if (path.events.empty()) { t = P.hours; break; }
    t = xe[0];
    LandscapeEvent E;
    E.time = t;
    if (path.events[0].forced) {
      E.mechanism = "the stores close for the night", E.match = "both stores", E.forced = true;
      E.footprint = fp({"day"}, {"a:A", "a:B", "day"});
      a = b = false;
      day++;
      out.events.push_back(std::move(E));
      continue;
    }
    const auto occ = occurrences(xe);
    double total = 0;
    for (const auto& q : occ) total += q.rate;
    double pick = rng.uniform() * total;
    std::size_t k = 0;
    for (; k + 1 < occ.size() && pick >= occ[k].rate; k++) pick -= occ[k].rate;
    const auto& q = occ[k];
    E.mechanism = q.mechanism, E.match = q.match;
    switch (q.kind) {
      case 0: a = true, E.footprint = fp({"a:A", "a:B", "lane"}, {"a:A"}); break;
      case 1: a = false, E.footprint = fp({"a:A"}, {"a:A"}); break;
      case 2: b = true, E.footprint = fp({"a:B", "a:A", "lane"}, {"a:B"}); break;
      case 3: b = false, E.footprint = fp({"a:B"}, {"a:B"}); break;
      case 4: permit = true, E.footprint = fp({"permit"}, {"permit"}); break;
      case 5: {
        E.rewrites = build_lane(D, plan);
        lane = E.rewrites.back().created[0];
        E.footprint = fp({"permit", "lane"}, {"lane"});
        for (const auto& s : E.rewrites) {
          E.footprint.reads.insert(s.footprint.reads.begin(), s.footprint.reads.end());
          E.footprint.writes.insert(s.footprint.writes.begin(), s.footprint.writes.end());
        }
        break;
      }
      case 6: {
        E.rewrites = close_lane(D, lane);
        lane = 0;
        plan = plan_lane(D, 0, 1);
        E.footprint = fp({"lane", "recognised"}, {"lane"});
        for (const auto& s : E.rewrites) {
          E.footprint.reads.insert(s.footprint.reads.begin(), s.footprint.reads.end());
          E.footprint.writes.insert(s.footprint.writes.begin(), s.footprint.writes.end());
        }
        break;
      }
      case 7: recognised = true, E.footprint = fp({"lane", "recognised"}, {"recognised"}); break;
    }
    out.events.push_back(std::move(E));
  }
  for (const auto id : D.of(0)) out.walkers1 += D.value(id, "walkers");
  return out;
}

}  // namespace urban::landscape
