// Evidence combined into breakpoints, LUMPY's way (Layer et al. 2014, Genome Biology 15, R84,
// "Breakpoint clustering" and "Combining evidence").
//
// - Merging. Two pieces of evidence may be of one junction when they have the same variety (the
//   contigs and orientations of their two breakends) and their intervals intersect on both sides.
//   The pieces of a variety are swept in order of their first interval, and those that intersect
//   on both sides are joined (union–find), so that a component is a set of evidence of nearby
//   junctions of one kind.
// - Combining. LUMPY's default combines by the product of the probability vectors over the
//   intersection of the intervals ("trimmed product"); when the intervals have no common point,
//   only the evidence containing the maximum of the summed distribution is kept. Here that rule
//   is applied repeatedly inside a component: the maximum of the summed (normalised)
//   distributions on the first side, then on the second side among the evidence containing it,
//   chooses the evidence of one breakpoint, whose product is its distribution; that evidence is
//   removed and the rest is combined again, while at least `min_support` pieces remain. LUMPY
//   reports the position of maximum probability and the interval holding 95% of the
//   distribution (VCF's CIPOS and CIEND); so does this, after trimming the tails below 0.1%.
// - Precise positions. A split read or an alignment's gap is the junction itself; when a
//   breakpoint has such evidence, its junction is the most frequent of their canonical
//   adjacencies (ties to the first), and the breakpoint is precise; otherwise it is the product's
//   maximum, imprecise, and its interval says how far it may be. An imprecise breakpoint whose
//   intervals (widened by 50 bases) hold a precise breakpoint of its variety is that junction's
//   pairs left over (their intervals missed the exact position), and joins it.
//
// Complexity: O(E log E + Σ |intervals|) for E pieces of evidence, the union–find nearly linear;
// the repeated combination of a component of k pieces is O(k² w) for intervals of width w.

#pragma once

#include <genomics/structural/breakend.hpp>
#include <genomics/structural/evidence.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <tuple>
#include <vector>

namespace genomics::structural {

struct ClusterOptions {
  std::size_t min_support = 4;   // pieces of evidence a breakpoint needs (LUMPY's evidence threshold in its simulations)
  double tail = 1e-3;            // trimmed from each distribution's tails before the interval
  double mass = 0.95;            // the confidence interval's mass
  std::size_t absorb = 50;       // how far outside its interval an imprecise breakpoint joins a precise one
};

struct Breakpoint {
  Adjacency junction;            // canonical when precise
  Interval a, b;                 // the combined distributions (log, unnormalised)
  std::size_t ci_a[2]{0, 0}, ci_b[2]{0, 0};  // the intervals holding `mass` of each side, [lo, hi]
  std::size_t pairs = 0, splits = 0, gaps = 0;
  double weight = 0;             // Σ of the evidence's weights
  bool precise = false;
  Homology homology;
  std::vector<std::uint32_t> evidence;  // indices into the evidence
  std::size_t support() const { return pairs + splits + gaps; }
  Signature signature() const { return structural::signature(junction); }
};

namespace detail {

struct UnionFind {
  std::vector<std::uint32_t> parent;
  explicit UnionFind(std::size_t n) : parent(n) { std::iota(parent.begin(), parent.end(), 0u); }
  std::uint32_t find(std::uint32_t x) {
    while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
    return x;
  }
  void join(std::uint32_t x, std::uint32_t y) {
    x = find(x);
    y = find(y);
    if (x != y) parent[std::max(x, y)] = std::min(x, y);
  }
};

inline bool intersect(const Interval& x, const Interval& y) { return x.begin < y.end() && y.begin < x.end(); }

// log Σ exp of a vector.
inline double log_total(const std::vector<double>& v) {
  double m = -std::numeric_limits<double>::infinity();
  for (double x : v) m = std::max(m, x);
  if (!std::isfinite(m)) return m;
  double s = 0;
  for (double x : v) s += std::exp(x - m);
  return m + std::log(s);
}

// The position of the maximum of Σ_e P_e (normalised) over the given evidence's intervals on one side.
template <class Side>
std::size_t summed_argmax(const std::vector<Evidence>& ev, const std::vector<std::uint32_t>& members, Side side) {
  std::size_t lo = SIZE_MAX, hi = 0;
  for (auto i : members) { lo = std::min(lo, side(ev[i]).begin); hi = std::max(hi, side(ev[i]).end()); }
  std::vector<double> sum(hi - lo, 0.0);
  for (auto i : members) {
    const Interval& v = side(ev[i]);
    const double z = log_total(v.log_p);
    for (std::size_t k = 0; k < v.log_p.size(); k++) sum[v.begin + k - lo] += std::exp(v.log_p[k] - z);
  }
  return lo + static_cast<std::size_t>(std::max_element(sum.begin(), sum.end()) - sum.begin());
}

// The product of the members' distributions on one side over the intersection of their intervals.
template <class Side>
Interval product(const std::vector<Evidence>& ev, const std::vector<std::uint32_t>& members, Side side) {
  std::size_t lo = 0, hi = SIZE_MAX;
  for (auto i : members) { lo = std::max(lo, side(ev[i]).begin); hi = std::min(hi, side(ev[i]).end()); }
  Interval out;
  out.begin = lo;
  if (hi <= lo) return out;
  out.log_p.assign(hi - lo, 0.0);
  for (auto i : members) {
    const Interval& v = side(ev[i]);
    const double z = log_total(v.log_p);
    for (std::size_t x = lo; x < hi; x++) out.log_p[x - lo] += v.log_p[x - v.begin] - z;
  }
  return out;
}

// The argmax and the central interval of `mass` of a distribution after trimming its tails.
inline void summarise(const Interval& v, double tail, double mass, std::size_t& at, std::size_t ci[2]) {
  if (v.log_p.empty()) { at = v.begin; ci[0] = ci[1] = v.begin; return; }
  const double z = log_total(v.log_p);
  std::vector<double> p(v.log_p.size());
  for (std::size_t k = 0; k < p.size(); k++) p[k] = std::exp(v.log_p[k] - z);
  at = v.begin + static_cast<std::size_t>(std::max_element(p.begin(), p.end()) - p.begin());
  // Trim the tails below `tail`, then the central interval of `mass` of what remains.
  std::size_t first = 0, last = p.size();
  double cut = 0;
  while (first < last && cut + p[first] < tail / 2) cut += p[first++];
  cut = 0;
  while (last > first + 1 && cut + p[last - 1] < tail / 2) cut += p[--last];
  double total = 0;
  for (std::size_t k = first; k < last; k++) total += p[k];
  const double side = total * (1 - mass) / 2;
  double acc = 0;
  std::size_t lo = first, hi = last - 1;
  while (lo < hi && acc + p[lo] < side) acc += p[lo++];
  acc = 0;
  while (hi > lo && acc + p[hi] < side) acc += p[hi--];
  ci[0] = v.begin + lo;
  ci[1] = v.begin + hi;
}

}  // namespace detail

// The breakpoints of the evidence, in order of their first breakend.
inline std::vector<Breakpoint> cluster(const std::vector<Codes>& contigs, const std::vector<Evidence>& ev, const ClusterOptions& o = {}) {
  std::vector<Breakpoint> out;
  // Order by variety, then by the first interval.
  const auto variety = [&](const Evidence& e) { return std::make_tuple(e.junction.a.contig, e.junction.b.contig, e.junction.a.plus, e.junction.b.plus); };
  std::vector<std::uint32_t> order(ev.size());
  std::iota(order.begin(), order.end(), 0u);
  std::sort(order.begin(), order.end(), [&](std::uint32_t x, std::uint32_t y) {
    const auto vx = variety(ev[x]), vy = variety(ev[y]);
    if (vx != vy) return vx < vy;
    if (ev[x].a.begin != ev[y].a.begin) return ev[x].a.begin < ev[y].a.begin;
    return x < y;
  });
  detail::UnionFind uf(ev.size());
  for (std::size_t i = 0; i < order.size(); i++) {
    const Evidence& x = ev[order[i]];
    for (std::size_t j = i + 1; j < order.size(); j++) {
      const Evidence& y = ev[order[j]];
      if (variety(y) != variety(x) || y.a.begin >= x.a.end()) break;
      if (detail::intersect(x.b, y.b)) uf.join(order[i], order[j]);
    }
  }
  std::map<std::uint32_t, std::vector<std::uint32_t>> components;
  for (auto i : order) components[uf.find(i)].push_back(i);
  const auto side_a = [](const Evidence& e) -> const Interval& { return e.a; };
  const auto side_b = [](const Evidence& e) -> const Interval& { return e.b; };
  for (auto& [root, all] : components) {
    std::vector<std::uint32_t> left = all;
    while (left.size() >= o.min_support) {
      const std::size_t xa = detail::summed_argmax(ev, left, side_a);
      std::vector<std::uint32_t> with_a;
      for (auto i : left)
        if (ev[i].a.contains(xa)) with_a.push_back(i);
      const std::size_t xb = detail::summed_argmax(ev, with_a, side_b);
      std::vector<std::uint32_t> members, rest;
      for (auto i : left) (ev[i].a.contains(xa) && ev[i].b.contains(xb) ? members : rest).push_back(i);
      if (members.empty()) break;
      if (members.size() >= o.min_support) {
        Breakpoint bp;
        bp.evidence = members;
        bp.a = detail::product(ev, members, side_a);
        bp.b = detail::product(ev, members, side_b);
        std::size_t at_a, at_b;
        detail::summarise(bp.a, o.tail, o.mass, at_a, bp.ci_a);
        detail::summarise(bp.b, o.tail, o.mass, at_b, bp.ci_b);
        const Evidence& first = ev[members.front()];
        bp.junction = {{first.junction.a.contig, at_a, first.junction.a.plus}, {first.junction.b.contig, at_b, first.junction.b.plus}, {}};
        // The most frequent exact junction among splits and gaps.
        std::map<std::tuple<std::size_t, std::size_t, std::string>, std::size_t> exact;
        for (auto i : members) {
          const Evidence& e = ev[i];
          bp.weight += e.weight;
          if (e.source == Source::pair) bp.pairs++;
          else {
            (e.source == Source::split ? bp.splits : bp.gaps)++;
            exact[{e.junction.a.position, e.junction.b.position, e.junction.inserted}]++;
          }
        }
        if (!exact.empty()) {
          auto best = exact.begin();
          for (auto it = exact.begin(); it != exact.end(); ++it)
            if (it->second > best->second) best = it;
          bp.junction.a.position = std::get<0>(best->first);
          bp.junction.b.position = std::get<1>(best->first);
          bp.junction.inserted = std::get<2>(best->first);
          bp.precise = true;
        }
        bp.homology = homology(contigs, bp.junction);
        out.push_back(std::move(bp));
      }
      left = std::move(rest);
    }
  }
  // An imprecise breakpoint whose confidence intervals (widened by `absorb`) hold a precise one of
  // its variety on both sides is the rest of that junction's pairs, whose intervals missed the
  // split reads' exact position: its evidence joins the precise one.
  std::vector<bool> gone(out.size(), false);
  for (std::size_t i = 0; i < out.size(); i++) {
    if (out[i].precise) continue;
    const Breakpoint& x = out[i];
    for (std::size_t j = 0; j < out.size(); j++) {
      Breakpoint& y = out[j];
      if (!y.precise || gone[j] || y.junction.a.contig != x.junction.a.contig || y.junction.b.contig != x.junction.b.contig || y.junction.a.plus != x.junction.a.plus ||
          y.junction.b.plus != x.junction.b.plus)
        continue;
      const auto within = [&](std::size_t p, const std::size_t ci[2]) { return p + o.absorb >= ci[0] && p <= ci[1] + o.absorb; };
      if (!within(y.junction.a.position, x.ci_a) || !within(y.junction.b.position, x.ci_b)) continue;
      y.pairs += x.pairs;
      y.weight += x.weight;
      y.evidence.insert(y.evidence.end(), x.evidence.begin(), x.evidence.end());
      gone[i] = true;
      break;
    }
  }
  std::vector<Breakpoint> kept;
  for (std::size_t i = 0; i < out.size(); i++)
    if (!gone[i]) kept.push_back(std::move(out[i]));
  std::sort(kept.begin(), kept.end(), [](const Breakpoint& x, const Breakpoint& y) {
    if (x.junction.a != y.junction.a) return x.junction.a < y.junction.a;
    return x.junction.b < y.junction.b;
  });
  return kept;
}

}  // namespace genomics::structural
