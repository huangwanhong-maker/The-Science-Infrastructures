// Root-to-tip regression of dated tips: the molecular clock read from a tree of sequences
// sampled at different times (TempEst: Rambaut, Lam, Max Carvalho and Pybus 2016, "Exploring the
// temporal structure of heterochronous sequences using TempEst", Virus Evolution 2, vew007).
//
// Under a strict clock of rate u, a tip sampled at time t_i lies at distance d_i = u (t_i − t_r)
// from the root, sampled at t_r (their eq. 1): the regression of distance on date has slope u
// and meets the date axis at t_r, the date of the root. On an unrooted tree the root is
// unknown: put it on the branch above v at distance x from v, and every tip's distance is
// a_i + s_i x, with a_i its distance from v, s_i = +1 below v and −1 elsewhere. Then the
// regression's sums are polynomials in x: with S_tt = Σ (t_i − t̄)², S_td(x) = c₀ + c₁ x and
// S_dd(x) = v₀ + v₁ x + v₂ x²,
//
//   the residual sum of squares  RSS(x) = S_dd(x) − S_td(x)² / S_tt,  a quadratic in x;
//   the coefficient of determination  R²(x) = S_td(x)² / (S_tt S_dd(x)),
//
// whose derivative vanishes, besides where S_td = 0, at the root of the linear equation
// (2c₁v₀ − c₀v₁) + (c₁v₁ − 2c₀v₂) x = 0. So on each branch the best x is found in closed form
// (the end points, the vertex of the quadratic, or that root), and the best root is the best
// over all branches: the least residual mean square RSS/(n − 2) (TempEst's default, solved
// analytically there too), or the greatest R² or correlation (TempEst's other objectives, which
// it finds numerically). With `positive_rate`, roots that give a negative slope are passed over
// when any gives a positive one. O(n²) for n tips: one traversal per branch.
//
// Dates are numbers on one scale (decimal years, say); a tip without a date (NaN) is left out
// of the regression. The estimates are those of the regression: the points are not
// independent (tips share branches), so no standard errors are given (Rambaut et al. 2016,
// §2.3, use it for exploration, not inference).

#pragma once

#include <evolution/phylogeny/tree.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace evolution::phylogeny {

enum class ClockObjective { residual_mean_square, r_squared, correlation };

struct ClockFit {
  Tree tree;                 // rooted at the best root
  double rate = 0;           // the slope: substitutions per site per unit of time
  double root_date = 0;      // the x-intercept: the date of the root
  double intercept = 0;      // the distance at date 0
  double r_squared = 0, correlation = 0, residual_mean_square = 0;
  int branch = -1;           // the branch of the input tree that holds the root (its lower node)
  double position = 0;       // the root's distance from that node
  std::vector<std::string> tips;   // the dated tips, in the order of the leaves
  std::vector<double> dates, distances, residuals;
};

namespace detail {

struct Regression {
  double rate = 0, intercept = 0, r2 = 0, r = 0, rms = 0;
};
inline Regression regress(const std::vector<double>& t, const std::vector<double>& d) {
  const double n = static_cast<double>(t.size());
  double mt = 0, md = 0;
  for (std::size_t i = 0; i < t.size(); i++) {
    mt += t[i];
    md += d[i];
  }
  mt /= n;
  md /= n;
  double stt = 0, sdd = 0, std_ = 0;
  for (std::size_t i = 0; i < t.size(); i++) {
    stt += (t[i] - mt) * (t[i] - mt);
    sdd += (d[i] - md) * (d[i] - md);
    std_ += (t[i] - mt) * (d[i] - md);
  }
  Regression g;
  g.rate = stt > 0 ? std_ / stt : 0;
  g.intercept = md - g.rate * mt;
  g.r2 = stt > 0 && sdd > 0 ? std_ * std_ / (stt * sdd) : 0;
  g.r = stt > 0 && sdd > 0 ? std_ / std::sqrt(stt * sdd) : 0;
  const double rss = std::max(0.0, sdd - (stt > 0 ? std_ * std_ / stt : 0));
  g.rms = n > 2 ? rss / (n - 2) : 0;
  return g;
}

}  // namespace detail

// The regression on a tree as it is rooted, for dates given by leaf label.
inline ClockFit root_to_tip_at_root(const Tree& t, const std::vector<std::string>& labels, const std::vector<double>& dates) {
  if (labels.size() != dates.size()) throw std::invalid_argument("evolution::phylogeny::root_to_tip: a date for each label");
  ClockFit f;
  f.tree = t;
  const auto depth = depths(t);
  for (int v : t.leaves())
    for (std::size_t k = 0; k < labels.size(); k++)
      if (labels[k] == t[v].label && std::isfinite(dates[k])) {
        f.tips.push_back(labels[k]);
        f.dates.push_back(dates[k]);
        f.distances.push_back(depth[static_cast<std::size_t>(v)]);
      }
  if (f.tips.size() < 3) throw std::invalid_argument("evolution::phylogeny::root_to_tip: at least three dated tips");
  const auto g = detail::regress(f.dates, f.distances);
  f.rate = g.rate;
  f.intercept = g.intercept;
  f.root_date = g.rate != 0 ? -g.intercept / g.rate : 0;
  f.r_squared = g.r2;
  f.correlation = g.r;
  f.residual_mean_square = g.rms;
  for (std::size_t i = 0; i < f.dates.size(); i++) f.residuals.push_back(f.distances[i] - (g.intercept + g.rate * f.dates[i]));
  return f;
}

// The best root of an unrooted (or rooted) tree for the regression of distance on date.
inline ClockFit root_to_tip(const Tree& t, const std::vector<std::string>& labels, const std::vector<double>& dates,
                            ClockObjective objective = ClockObjective::residual_mean_square, bool positive_rate = true) {
  if (labels.size() != dates.size()) throw std::invalid_argument("evolution::phylogeny::root_to_tip: a date for each label");
  // The dated leaves.
  std::vector<int> leaf;
  std::vector<double> date;
  for (int v : t.leaves())
    for (std::size_t k = 0; k < labels.size(); k++)
      if (labels[k] == t[v].label && std::isfinite(dates[k])) {
        leaf.push_back(v);
        date.push_back(dates[k]);
      }
  const std::size_t n = leaf.size();
  if (n < 3) throw std::invalid_argument("evolution::phylogeny::root_to_tip: at least three dated tips");
  double mt = 0;
  for (double x : date) mt += x;
  mt /= static_cast<double>(n);
  double stt = 0;
  for (double x : date) stt += (x - mt) * (x - mt);
  if (!(stt > 0)) throw std::invalid_argument("evolution::phylogeny::root_to_tip: the tips' dates differ");
  struct Best {
    double score = -std::numeric_limits<double>::infinity();
    int branch = -1;
    double x = 0;
    bool positive = false;
  } best;
  std::vector<char> inside(t.size(), 0);
  for (std::size_t vi = 0; vi < t.size(); vi++) {
    const int v = static_cast<int>(vi);
    if (v == t.root) continue;
    const double len = t.length_or(v, 0);
    const auto a_all = distances_from(t, v);
    // Which nodes lie below v.
    std::fill(inside.begin(), inside.end(), 0);
    std::vector<int> stack{v};
    while (!stack.empty()) {
      const int u = stack.back();
      stack.pop_back();
      inside[static_cast<std::size_t>(u)] = 1;
      for (int c : t[u].children) stack.push_back(c);
    }
    double ma = 0, ms = 0;
    std::vector<double> a(n), s(n);
    for (std::size_t i = 0; i < n; i++) {
      a[i] = a_all[static_cast<std::size_t>(leaf[i])];
      s[i] = inside[static_cast<std::size_t>(leaf[i])] ? 1.0 : -1.0;
      ma += a[i];
      ms += s[i];
    }
    ma /= static_cast<double>(n);
    ms /= static_cast<double>(n);
    double c0 = 0, c1 = 0, v0 = 0, v1 = 0, v2 = 0;
    for (std::size_t i = 0; i < n; i++) {
      const double dt = date[i] - mt, da = a[i] - ma, ds = s[i] - ms;
      c0 += dt * da;
      c1 += dt * ds;
      v0 += da * da;
      v1 += 2 * da * ds;
      v2 += ds * ds;
    }
    std::vector<double> xs{0, len};
    if (objective == ClockObjective::residual_mean_square) {
      const double q2 = v2 - c1 * c1 / stt, q1 = v1 - 2 * c0 * c1 / stt;
      if (q2 > 0) xs.push_back(-q1 / (2 * q2));
    } else {
      const double den = c1 * v1 - 2 * c0 * v2;
      if (den != 0) xs.push_back((c0 * v1 - 2 * c1 * v0) / den);
    }
    for (double x : xs) {
      if (!(x >= 0 && x <= len)) continue;
      const double sdt = c0 + c1 * x, sdd = v0 + v1 * x + v2 * x * x;
      const double rss = std::max(0.0, sdd - sdt * sdt / stt);
      double score;
      switch (objective) {
        case ClockObjective::residual_mean_square: score = -rss; break;
        case ClockObjective::r_squared: score = sdd > 0 ? sdt * sdt / (stt * sdd) : 0; break;
        default: score = sdd > 0 ? sdt / std::sqrt(stt * sdd) : 0; break;
      }
      const bool positive = sdt > 0;
      const bool better = positive_rate ? (positive && !best.positive) || (positive == best.positive && score > best.score) : score > best.score;
      if (better) best = {score, v, x, positive};
    }
  }
  if (best.branch < 0) throw std::invalid_argument("evolution::phylogeny::root_to_tip: no branch holds a root");
  // reroot measures x from the branch's lower node, as here.
  ClockFit f = root_to_tip_at_root(reroot(t, best.branch, best.x), labels, dates);
  f.branch = best.branch;
  f.position = best.x;
  return f;
}

}  // namespace evolution::phylogeny
