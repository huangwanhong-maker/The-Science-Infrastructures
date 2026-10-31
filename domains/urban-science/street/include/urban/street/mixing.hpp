// The departure of a street's arrangement from the fully mixed one (theory/street-oceanography.md §2.3):
// the mutual information M between place and use (Definition 2.6), on physis's mixing, with its exact
// partition among scales (Theorem 2.7).
//
// Places along a street are bins of arc length (`place` metres); the places are grouped into stretches
// (`stretch` metres, a multiple of the place) for the decomposition M = M_between(stretches) +
// Σ_g p_g M_within(g). For points in a plane (a sensor's surroundings), places are the cells of a square
// grid and stretches the cells of a coarser one.
//
// Theorem 2.8 (exchange blind to kind never raises M) is checked on the street's own counts: exchange_step()
// moves a share of each place's holdings to its neighbours, whatever their use.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <calx/linalg/rectangular.hpp>
#include <physis/statistical_mechanics/mixing.hpp>
#include <urban/street/street.hpp>

namespace urban::street {

namespace sm = physis::statistical_mechanics;

struct StreetMixing {
  double place = 0, stretch = 0;                // m
  calx::linalg::Rectangular counts{0, 0};       // places × uses
  std::vector<std::size_t> group;               // the stretch of each place
  sm::Mixing mixing;                            // M, H(U), H(U|X), local terms
  sm::MixingDecomposition scales;               // M = between + within
};

// Counts of uses by place, places given by an index per point; empty places are dropped.
inline StreetMixing mixing_of(const std::vector<std::size_t>& place_of, const std::vector<std::size_t>& use, std::size_t uses,
                              const std::vector<std::size_t>& group_of_place) {
  if (place_of.size() != use.size()) throw std::invalid_argument("urban::street::mixing_of: one place per point");
  std::size_t places = 0;
  for (auto p : place_of) places = std::max(places, p + 1);
  std::vector<std::vector<double>> rows(places, std::vector<double>(uses, 0));
  for (std::size_t i = 0; i < use.size(); i++) rows[place_of[i]][use[i]] += 1;
  StreetMixing out;
  std::vector<std::vector<double>> kept;
  for (std::size_t p = 0; p < places; p++) {
    double total = 0;
    for (double c : rows[p]) total += c;
    if (total <= 0) continue;
    kept.push_back(rows[p]);
    out.group.push_back(p < group_of_place.size() ? group_of_place[p] : 0);
  }
  if (kept.empty()) throw std::invalid_argument("urban::street::mixing_of: no points");
  out.counts = calx::linalg::Rectangular::from_rows(kept);
  out.mixing = sm::mixing(out.counts);
  out.scales = sm::mixing_decomposition(out.counts, out.group);
  return out;
}

// The street's holdings in places of `place` metres grouped into stretches of `stretch` metres.
inline StreetMixing street_mixing(const Street& st, double place, double stretch) {
  if (!(place > 0 && stretch >= place)) throw std::invalid_argument("urban::street::street_mixing: a place shorter than its stretch");
  std::vector<std::size_t> place_of, use;
  for (const auto& h : st.holdings) {
    place_of.push_back(static_cast<std::size_t>(std::max(0.0, std::min(st.length, h.s)) / place));
    use.push_back(h.use);
  }
  std::size_t places = 0;
  for (auto p : place_of) places = std::max(places, p + 1);
  std::vector<std::size_t> group(places);
  for (std::size_t p = 0; p < places; p++) group[p] = static_cast<std::size_t>((static_cast<double>(p) + 0.5) * place / stretch);
  StreetMixing out = mixing_of(place_of, use, st.uses.size(), group);
  out.place = place, out.stretch = stretch;
  return out;
}

// Points in a plane (x, y metres) in cells of `cell` metres, grouped into blocks of `block` metres.
inline StreetMixing plane_mixing(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::size_t>& use, std::size_t uses,
                                 double cell, double block) {
  if (x.size() != y.size() || x.size() != use.size()) throw std::invalid_argument("urban::street::plane_mixing: one position and use per point");
  if (!(cell > 0 && block >= cell)) throw std::invalid_argument("urban::street::plane_mixing: a cell smaller than its block");
  double x0 = 0, y0 = 0;
  for (std::size_t i = 0; i < x.size(); i++) x0 = std::min(x0, x[i]), y0 = std::min(y0, y[i]);
  std::vector<long long> key;
  for (std::size_t i = 0; i < x.size(); i++) key.push_back(static_cast<long long>((x[i] - x0) / cell) * 100000 + static_cast<long long>((y[i] - y0) / cell));
  std::vector<long long> distinct = key;
  std::sort(distinct.begin(), distinct.end());
  distinct.erase(std::unique(distinct.begin(), distinct.end()), distinct.end());
  std::vector<std::size_t> place_of;
  for (auto k : key) place_of.push_back(static_cast<std::size_t>(std::lower_bound(distinct.begin(), distinct.end(), k) - distinct.begin()));
  const long long per = std::max(1LL, static_cast<long long>(block / cell));
  std::vector<long long> block_key;
  for (auto k : distinct) block_key.push_back((k / 100000) / per * 100000 + (k % 100000) / per);
  std::vector<long long> blocks = block_key;
  std::sort(blocks.begin(), blocks.end());
  blocks.erase(std::unique(blocks.begin(), blocks.end()), blocks.end());
  std::vector<std::size_t> group;
  for (auto b : block_key) group.push_back(static_cast<std::size_t>(std::lower_bound(blocks.begin(), blocks.end(), b) - blocks.begin()));
  StreetMixing out = mixing_of(place_of, use, uses, group);
  out.place = cell, out.stretch = block;
  return out;
}

// One step of exchange blind to kind along a line of places: each place keeps 1 − 2ε of its holdings and
// sends ε to each neighbour (reflecting at the ends). Theorem 2.8: M never rises.
inline calx::linalg::Rectangular exchange_step(const calx::linalg::Rectangular& counts, double epsilon) {
  if (!(epsilon >= 0 && epsilon <= 0.5)) throw std::invalid_argument("urban::street::exchange_step: 0 ≤ ε ≤ ½");
  calx::linalg::Rectangular out(counts.rows, counts.cols);
  for (std::size_t p = 0; p < counts.rows; p++)
    for (std::size_t u = 0; u < counts.cols; u++) {
      const double c = counts(p, u);
      const std::size_t left = p == 0 ? p : p - 1, right = p + 1 == counts.rows ? p : p + 1;
      out(p, u) += (1 - 2 * epsilon) * c;
      out(left, u) += epsilon * c;
      out(right, u) += epsilon * c;
    }
  return out;
}

}  // namespace urban::street
