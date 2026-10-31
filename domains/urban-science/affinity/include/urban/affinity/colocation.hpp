// The affinity of uses (theory/street-oceanography.md, Definition 5.1; Huang, Street Oceanography, §11): the
// co-location quotient of every ordered pair of uses, with its significance against two nulls: marks
// permuted freely, and marks permuted within cells of terrain, the null of independent response to a shared
// terrain (Street Oceanography §28.2, Table 5: "variation in composition at fixed terrain").
//
// A_ij = ln CLQ_{i→j} (tyche `spatial`, Leslie and Kronenfeld 2011). Its asymmetry is the paper's requirement
// that A is not symmetric (§11.1). The strata are the cells of a square grid of `cell` metres: uses that merely
// favour the same blocks keep CLQ near 1 under the stratified null, while an attraction within blocks survives it.

#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <tyche/spatial/point_patterns.hpp>

namespace urban::affinity {

struct Affinity {
  std::vector<std::string> uses;
  std::size_t marks = 0;
  tyche::spatial::Colocation free, stratified;   // the two nulls
  std::vector<double> a;                          // ln CLQ, row-major marks × marks (NaN where undefined)
  double cell = 0;                                // m
  double radius = 0;                              // m; 0: nearest neighbours
};

inline Affinity affinity(const std::vector<double>& x, const std::vector<double>& y, const std::vector<std::size_t>& use, const std::vector<std::string>& uses,
                         double cell, double radius = 0, std::size_t permutations = 999, std::uint64_t seed = 1) {
  if (x.size() != y.size() || x.size() != use.size()) throw std::invalid_argument("urban::affinity::affinity: one position and use per point");
  if (!(cell > 0)) throw std::invalid_argument("urban::affinity::affinity: a positive cell");
  double x0 = 0, y0 = 0;
  for (std::size_t i = 0; i < x.size(); i++) x0 = std::min(x0, x[i]), y0 = std::min(y0, y[i]);
  std::vector<std::size_t> stratum;
  for (std::size_t i = 0; i < x.size(); i++)
    stratum.push_back(static_cast<std::size_t>((x[i] - x0) / cell) * 100003u + static_cast<std::size_t>((y[i] - y0) / cell));
  Affinity out;
  out.uses = uses, out.marks = uses.size(), out.cell = cell, out.radius = radius;
  out.free = tyche::spatial::colocation(x, y, use, uses.size(), {}, radius, permutations, seed);
  out.stratified = tyche::spatial::colocation(x, y, use, uses.size(), stratum, radius, permutations, seed + 1);
  for (double c : out.free.clq) out.a.push_back(std::isfinite(c) && c > 0 ? std::log(c) : std::nan(""));
  return out;
}

}  // namespace urban::affinity
