// The walkers' departures calibrated to the counters, hour by hour and each way (theory/movement.md §5).
//
// A counter on a footpath reads share × (1 ± ρ) of the walkers passing its line each way (walkers.hpp, Counter). The
// expected passes of a stream's trips are the logit flow on the counter's two directed states (route.hpp): linear in
// the stream's rate. So the counts of an hour h are linear in the streams' amplitudes, and grouping the streams into
// patterns (each a fixed law of origins and targets, as the flow study's arrive, leave, through, local and visit) the
// amplitudes a_g(h) ≥ 0 are the non-negative least squares (tyche's NNLS, Lawson and Hanson 1974) of
//
//   y_kw(h) ≈ Σ_g a_g(h) Σ_{s ∈ g} r_s(h) u_skw,     u_skw = share_k (1 + ρ_kw) Σ_pairs p_s(o, A) f_{way_kw}(o → A),
//
// over the counters k and their ways w (a missing count leaves its row out). r_s is the stream's rate as given (its
// shape in time and its scale), u its unit response. A stream that dwells counts its walk back too, in the hour it
// departs: the quasi-static assumption of the flow study, exact for trips shorter than the hour on average; a caller
// that needs the lag (a crowd that stays) composes its own columns from unit_counts().

#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <calx/linalg/rectangular.hpp>
#include <tyche/regression/nnls.hpp>

#include <urban/movement/route.hpp>
#include <urban/movement/walkers.hpp>

namespace urban::movement {

// The expected passes each way at each counter of one trip of the stream (its pairs by their law), and of its walk
// back when it dwells. Unit: passes a trip, before the counter's share.
struct UnitCounts {
  std::vector<std::array<double, 2>> out, back;   // per counter
};

inline UnitCounts unit_counts(const Values& V, const Targets& T, const Stream& s, const std::vector<Counter>& counters) {
  double total = 0;
  for (const double w : s.weight) total += w;
  UnitCounts u;
  u.out.assign(counters.size(), {0.0, 0.0}), u.back.assign(counters.size(), {0.0, 0.0});
  if (!(total > 0)) return u;
  std::vector<Trip> out, back;
  for (std::size_t p = 0; p < s.origin.size(); p++) {
    if (s.weight[p] == 0) continue;
    out.push_back({s.origin[p], s.target[p], s.weight[p] / total});
    if (s.dwell > 0) {
      // The walk back starts where the walk out ended: each junction the target's states enter, by their share.
      const auto& z = V[s.target[p]];
      const auto f = V.routes().visits(s.origin[p], T.states(s.target[p]), z);
      for (const std::size_t c : T.states(s.target[p]))
        if (f[c] > 0) back.push_back({V.routes().walkway().head(c), s.back[p], s.weight[p] / total * f[c]});
    }
  }
  const auto fo = expected_flow(V, T, out);
  const auto fb = back.empty() ? std::vector<double>(fo.size(), 0.0) : expected_flow(V, T, back);
  for (std::size_t k = 0; k < counters.size(); k++)
    for (int w = 0; w < 2; w++) u.out[k][w] = fo[counters[k].way[w]], u.back[k][w] = fb[counters[k].way[w]];
  return u;
}

struct Calibration {
  std::vector<std::vector<double>> amplitude;     // [group][hour]
  std::vector<std::vector<std::array<double, 2>>> fitted;   // [counter][hour], each way, as the counter reads
  double r2 = 0;                                  // over every count used
  std::size_t rows = 0;
};

// observed[k][h][w]: what counter k read each way in hour h (NaN where missing); groups: the streams sharing each
// amplitude (a stream in no group keeps its rate: an offset). Hours: those of the streams' rates.
inline Calibration calibrate(const std::vector<Stream>& streams, const std::vector<UnitCounts>& unit, const std::vector<Counter>& counters,
                             const std::vector<std::vector<std::array<double, 2>>>& observed, const std::vector<std::vector<std::size_t>>& groups) {
  if (unit.size() != streams.size()) throw std::invalid_argument("urban::movement::calibrate: a unit response for each stream");
  std::size_t H = 0;
  for (const auto& s : streams) H = std::max(H, s.rate.size());
  const std::size_t K = counters.size(), G = groups.size();
  std::vector<char> grouped(streams.size(), 0);
  for (const auto& g : groups)
    for (const std::size_t s : g) grouped.at(s) = 1;
  const auto reads = [&](std::size_t s, std::size_t k, int w, std::size_t h) {
    const double r = h < streams[s].rate.size() ? streams[s].rate[h] : 0.0;
    return r * counters[k].share * (1 + counters[k].hand[w]) * (unit[s].out[k][w] + unit[s].back[k][w]);
  };
  Calibration C;
  C.amplitude.assign(G, std::vector<double>(H, 0.0));
  C.fitted.assign(K, std::vector<std::array<double, 2>>(H, {0.0, 0.0}));
  double ss = 0, st = 0, mean = 0, n = 0;
  for (std::size_t k = 0; k < K; k++)
    for (std::size_t h = 0; h < H; h++)
      for (int w = 0; w < 2; w++)
        if (std::isfinite(observed[k][h][w])) mean += observed[k][h][w], n++;
  mean = n > 0 ? mean / n : 0;
  for (std::size_t h = 0; h < H; h++) {
    std::vector<std::vector<double>> rows;
    std::vector<double> b;
    std::vector<std::array<std::size_t, 2>> at;
    for (std::size_t k = 0; k < K; k++)
      for (int w = 0; w < 2; w++) {
        if (!std::isfinite(observed[k][h][w])) continue;
        std::vector<double> a(G, 0.0);
        double offset = 0;
        for (std::size_t g = 0; g < G; g++)
          for (const std::size_t s : groups[g]) a[g] += reads(s, k, w, h);
        for (std::size_t s = 0; s < streams.size(); s++)
          if (!grouped[s]) offset += reads(s, k, w, h);
        rows.push_back(std::move(a)), b.push_back(observed[k][h][w] - offset), at.push_back({k, static_cast<std::size_t>(w)});
      }
    if (rows.empty()) continue;
    const auto x = G > 0 ? tyche::regression::nnls(calx::linalg::Rectangular::from_rows(rows), b).x : std::vector<double>{};
    for (std::size_t g = 0; g < G; g++) C.amplitude[g][h] = x[g];
    for (std::size_t k = 0; k < K; k++)
      for (int w = 0; w < 2; w++) {
        double v = 0;
        for (std::size_t s = 0; s < streams.size(); s++) {
          double a = 1;
          if (grouped[s])
            for (std::size_t g = 0; g < G; g++)
              if (std::find(groups[g].begin(), groups[g].end(), s) != groups[g].end()) a = x[g];
          v += a * reads(s, k, w, h);
        }
        C.fitted[k][h][static_cast<std::size_t>(w)] = v;
      }
    for (std::size_t r = 0; r < rows.size(); r++) {
      const double y = observed[at[r][0]][h][at[r][1]], p = C.fitted[at[r][0]][h][at[r][1]];
      ss += (y - p) * (y - p), st += (y - mean) * (y - mean);
    }
    C.rows += rows.size();
  }
  C.r2 = st > 0 ? 1 - ss / st : 0;
  return C;
}

// The streams with their rates times the fitted amplitudes of their groups.
inline std::vector<Stream> calibrated(std::vector<Stream> streams, const Calibration& C, const std::vector<std::vector<std::size_t>>& groups) {
  for (std::size_t g = 0; g < groups.size(); g++)
    for (const std::size_t s : groups[g])
      for (std::size_t h = 0; h < streams[s].rate.size() && h < C.amplitude[g].size(); h++) streams[s].rate[h] *= C.amplitude[g][h];
  return streams;
}

}  // namespace urban::movement
