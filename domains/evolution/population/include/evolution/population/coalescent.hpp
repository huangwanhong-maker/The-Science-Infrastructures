// The coalescent: genealogies of a sample, with mutation, changing population size and
// recombination, simulated backwards in time from a seed.
//
// Time is in units of 2N₀ generations, N₀ the present number of diploid individuals (2N₀
// gene copies), so that two lineages coalesce at rate 1 and k lineages at rate C(k, 2)
// (Kingman 1982, Stoch. Proc. Appl. 13, 235–248; Hudson 1990, Oxford Surveys in Evolutionary
// Biology 7, 1–44). ms (Hudson 2002) measures time in units of 4N₀ generations: divide our
// times by 2 to compare (ms_time). In these units E[T_MRCA] = 2(1 − 1/n) and the expected
// total branch length is 2a_n, a_n = Σ_{i=1}^{n−1} 1/i.
//
// - Population size N(t) = N₀ λ(t) is piecewise exponential back in time (Demography): epochs
//   starting at times t_e, each with a relative size λ_e and a growth rate α_e, λ(t) =
//   λ_e e^{−α_e (t − t_e)}, as ms's -eN and -eG (its α in 4N₀ units is twice ours). Coalescence
//   then runs on the clock Λ(t) = ∫₀ᵗ du/λ(u): the next coalescence among k lineages after t is
//   at Λ⁻¹(Λ(t) + E/C(k, 2)), E exponential (time rescaling; Griffiths and Tavaré 1994, Phil.
//   Trans. R. Soc. B 344, 403–410; for exponential growth Slatkin and Hudson 1991, Genetics 129,
//   555–562). Λ is inverted in closed form on each epoch.
// - Mutation is the infinite-sites model: mutations fall on the branches as a Poisson process
//   of rate θ/2 per unit of branch length and of chromosome (θ = 4N₀μ for the whole chromosome),
//   each at a new position uniform on the branch's interval and at a time uniform on the branch;
//   a mutation on the branch above node v is carried by the samples below v. So
//   E[S] = (θ/2) E[L] = θ a_n (Watterson 1975).
// - Recombination is Hudson's algorithm (Hudson 1983, Theor. Popul. Biol. 23, 183–201; as ms,
//   Hudson 2002, Bioinformatics 18, 337–338) on the continuous chromosome [0, 1) with ρ = 4N₀r
//   for the whole chromosome. Each lineage carries the segments of ancestral material it holds
//   for the sample, each segment tagged with the node it descends to. A lineage whose material
//   spans [L, R) recombines at rate (ρ/2)(R − L), at a breakpoint uniform on (L, R), and splits
//   its material there; recombinations outside the span change nothing the sample can see and
//   are not simulated. When two lineages meet, each interval where both carry material gets a
//   new node (one node per event), with an edge to each child over the interval; an interval
//   where the merged lineage is the last one carrying it has reached its most recent common
//   ancestor and is dropped (the counts of lineages carrying each position decide). The
//   simulation ends when no material is left. The result is the ancestral recombination graph
//   restricted to the sample, as a set of nodes with times and of edges (left, right, parent,
//   child) whose intervals partition the chromosome into marginal trees: the representation of
//   Kelleher, Etheridge and McVean (2016, PLoS Comput. Biol. 12, e1004842), whose Algorithm T
//   sweeps the trees left to right by inserting and removing edges (for_each_tree).
//
// Draws: the genealogy from a stream keyed by (seed, COALESCENT, replicate), the mutations
// from (seed, COALESCENT_MUTATIONS, replicate), so that the genealogy does not depend on θ and
// a replicate is the same on any thread.
//
// Complexity: without recombination n − 1 events of O(1) work, O(n) in all. With
// recombination each event costs O(k + B) for k lineages and B breakpoints, and the number of
// events grows faster than linearly in ρ (about 1100 recombinations for n = 100, ρ = 100). The
// sweep over T trees costs O(E log E + T n) for E edges; mutation O(E + S); the haplotype
// matrix O(S n + T n d), d the depth of a tree.

#pragma once

#include <evolution/population/genotypes.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace evolution::population {

// The purposes of this module's streams.
inline constexpr std::uint32_t COALESCENT = 0x434f414c, COALESCENT_MUTATIONS = 0x4d555453;

inline double ms_time(double t) { return t / 2; }
inline double expected_tmrca(std::size_t n) { return n < 2 ? 0 : 2 * (1 - 1.0 / static_cast<double>(n)); }
inline double expected_tree_length(std::size_t n) {
  double a = 0;
  for (std::size_t i = n; i-- > 1;) a += 1.0 / static_cast<double>(i);
  return 2 * a;
}

// ---- demography --------------------------------------------------------------------------------------

struct Epoch {
  double start = 0;   // back in time, units of 2N₀ generations
  double size = 1;    // λ at the start, relative to N₀
  double growth = 0;  // α: λ(t) = size · e^{−α (t − start)}; α > 0 is growth towards the present
};

class Demography {
 public:
  Demography() : Demography(std::vector<Epoch>{Epoch{}}) {}
  explicit Demography(std::vector<Epoch> epochs) : epochs_(std::move(epochs)) {
    if (epochs_.empty() || epochs_[0].start != 0) throw std::invalid_argument("evolution::population::Demography: the first epoch starts at 0");
    for (std::size_t e = 0; e < epochs_.size(); e++) {
      if (!(epochs_[e].size > 0) || !std::isfinite(epochs_[e].growth)) throw std::invalid_argument("evolution::population::Demography: sizes are positive");
      if (e && !(epochs_[e].start > epochs_[e - 1].start)) throw std::invalid_argument("evolution::population::Demography: epochs in increasing order");
    }
    if (epochs_.back().growth < 0) throw std::invalid_argument("evolution::population::Demography: the oldest epoch's growth is >= 0, so that lineages always meet");
    start_intensity_.assign(epochs_.size(), 0);
    for (std::size_t e = 1; e < epochs_.size(); e++)
      start_intensity_[e] = start_intensity_[e - 1] + within(e - 1, epochs_[e].start - epochs_[e - 1].start);
  }

  static Demography constant() { return {}; }
  // λ(t) = e^{−α t}: a population that has grown exponentially at rate α to its present size.
  static Demography exponential(double growth) { return Demography({Epoch{0, 1, growth}}); }
  // Size `size` (relative) from `start` to `start + duration`, and 1 before and after.
  static Demography bottleneck(double start, double duration, double size) {
    std::vector<Epoch> e;
    if (start > 0) e.push_back({0, 1, 0});
    e.push_back({start, size, 0});
    e.push_back({start + duration, 1, 0});
    return Demography(std::move(e));
  }

  const std::vector<Epoch>& epochs() const { return epochs_; }

  double size(double t) const {
    const Epoch& e = epochs_[epoch_of(t)];
    return e.size * std::exp(-e.growth * (t - e.start));
  }
  // Λ(t) = ∫₀ᵗ du / λ(u).
  double intensity(double t) const {
    const std::size_t e = epoch_of(t);
    return start_intensity_[e] + within(e, t - epochs_[e].start);
  }
  // The s ≥ t with Λ(s) − Λ(t) = x; infinite when Λ stays below.
  double advance(double t, double x) const {
    std::size_t e = epoch_of(t);
    for (;;) {
      const Epoch& ep = epochs_[e];
      const double u = t - ep.start, a = ep.growth, lam = ep.size;
      const double end = e + 1 < epochs_.size() ? epochs_[e + 1].start : std::numeric_limits<double>::infinity();
      // The intensity left in this epoch after t.
      double left = std::numeric_limits<double>::infinity();
      if (std::isfinite(end)) left = a == 0 ? (end - t) / lam : std::exp(a * u) * std::expm1(a * (end - t)) / (a * lam);
      else if (a < 0) left = std::exp(a * u) / (-a * lam);
      if (x < left) {
        if (a == 0) return t + x * lam;
        // e^{a(s − start)} = e^{a u} + a λ x, written to avoid overflow of e^{a u}.
        const double arg = a * lam * x * std::exp(-a * u);
        return t + std::log1p(arg) / a;
      }
      if (!std::isfinite(end)) return std::numeric_limits<double>::infinity();
      x -= left;
      t = end;
      e++;
    }
  }

 private:
  std::size_t epoch_of(double t) const {
    std::size_t e = 0;
    while (e + 1 < epochs_.size() && epochs_[e + 1].start <= t) e++;
    return e;
  }
  // ∫ of 1/λ over the first u of epoch e.
  double within(std::size_t e, double u) const {
    const double a = epochs_[e].growth, lam = epochs_[e].size;
    return a == 0 ? u / lam : std::expm1(a * u) / (a * lam);
  }
  std::vector<Epoch> epochs_;
  std::vector<double> start_intensity_;
};

// ---- the genealogy --------------------------------------------------------------------------------------

struct CoalescentSpec {
  std::size_t samples = 10;
  double theta = 0;  // 4N₀μ for the whole chromosome
  double rho = 0;    // 4N₀r for the whole chromosome
  Demography demography;
};

// Node `parent` is the parent of node `child` on [left, right) of the chromosome.
struct Edge {
  double left = 0, right = 1;
  std::int32_t parent = -1, child = -1;
};

struct Mutation {
  double position = 0, time = 0;
  std::int32_t node = -1;  // the mutation lies on the branch above this node
};

struct Genealogy {
  std::size_t samples = 0;
  std::vector<double> time;          // of each node: the samples 0 … n − 1 at 0, then the ancestors in order of time
  std::vector<Edge> edges;           // sorted by parent, then child, then left; contiguous intervals joined
  std::vector<Mutation> mutations;   // in order of position
  std::size_t recombinations = 0;    // breakpoints within ancestral material
  std::size_t common_ancestors = 0;  // meetings of two lineages, with or without shared material

  std::size_t nodes() const { return time.size(); }
};

namespace detail {

struct Segment {
  double left, right;
  std::int32_t node;
};

// The number of lineages that carry each position: a step function on [0, 1).
class Coverage {
 public:
  explicit Coverage(std::int32_t n) : at_{0.0, 1.0}, count_{n, 0} {}
  // Lowers the count on [l, r) by one, and calls keep(a, b) on each piece still carried by
  // more than one lineage afterwards.
  template <class Keep>
  void merge(double l, double r, Keep&& keep) {
    const std::size_t i = split(l), j = split(r);
    for (std::size_t p = i; p < j; p++)
      if (--count_[p] > 1) keep(at_[p], at_[p + 1]);
  }

 private:
  std::size_t split(double x) {
    const auto it = std::lower_bound(at_.begin(), at_.end(), x);
    const std::size_t i = static_cast<std::size_t>(it - at_.begin());
    if (i < at_.size() && at_[i] == x) return i;
    at_.insert(it, x);
    count_.insert(count_.begin() + static_cast<std::ptrdiff_t>(i), count_[i - 1]);
    return i;
  }
  std::vector<double> at_;
  std::vector<std::int32_t> count_;
};

inline void append(std::vector<Segment>& out, Segment s) {
  if (!out.empty() && out.back().node == s.node && out.back().right == s.left) out.back().right = s.right;
  else out.push_back(s);
}

}  // namespace detail

// The genealogy of a sample (with its mutations), replicate `replicate` of the seed.
inline Genealogy coalescent(const CoalescentSpec& spec, std::uint64_t seed, std::uint32_t replicate = 0) {
  using detail::Segment;
  const std::size_t n = spec.samples;
  if (n == 0 || n > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max() / 4))
    throw std::invalid_argument("evolution::population::coalescent: a sample of at least one");
  if (!(spec.theta >= 0) || !(spec.rho >= 0)) throw std::invalid_argument("evolution::population::coalescent: θ, ρ >= 0");
  Genealogy g;
  g.samples = n;
  g.time.assign(n, 0.0);
  calx::random::Stream r(seed, COALESCENT, replicate, 0);

  std::vector<std::vector<Segment>> lineages(n);
  for (std::size_t i = 0; i < n; i++) lineages[i] = {Segment{0, 1, static_cast<std::int32_t>(i)}};
  detail::Coverage coverage(static_cast<std::int32_t>(n));
  std::vector<Segment> merged;
  double t = 0;
  const double half_rho = spec.rho / 2;

  while (lineages.size() >= 2) {
    const std::size_t k = lineages.size();
    double span = 0;
    if (half_rho > 0)
      for (const auto& l : lineages) span += l.back().right - l.front().left;
    const double rec_rate = half_rho * span;
    const double t_rec = rec_rate > 0 ? t + r.exponential() / rec_rate : std::numeric_limits<double>::infinity();
    const double t_ca = spec.demography.advance(t, r.exponential() / (0.5 * static_cast<double>(k) * static_cast<double>(k - 1)));
    if (!std::isfinite(t_rec) && !std::isfinite(t_ca)) throw std::runtime_error("evolution::population::coalescent: the lineages never meet");

    if (t_rec < t_ca) {
      // A breakpoint uniform on the total span; the lineage that holds it splits there.
      t = t_rec;
      double u = r.uniform() * span;
      std::size_t i = 0;
      for (; i + 1 < k; i++) {
        const double s = lineages[i].back().right - lineages[i].front().left;
        if (u < s) break;
        u -= s;
      }
      std::vector<Segment>& x = lineages[i];
      const double L = x.front().left, R = x.back().right;
      const double at = std::clamp(L + u, std::nextafter(L, R), std::nextafter(R, L));
      std::vector<Segment> right;
      std::vector<Segment> left;
      for (const Segment& s : x) {
        if (s.right <= at) left.push_back(s);
        else if (s.left >= at) right.push_back(s);
        else {
          left.push_back({s.left, at, s.node});
          right.push_back({at, s.right, s.node});
        }
      }
      x = std::move(left);
      lineages.push_back(std::move(right));
      g.recombinations++;
      continue;
    }

    // Two lineages, uniformly, meet.
    t = t_ca;
    g.common_ancestors++;
    std::size_t a = static_cast<std::size_t>(r.below(k)), b = static_cast<std::size_t>(r.below(k - 1));
    if (b >= a) b++;
    const std::vector<Segment>& X = lineages[a];
    const std::vector<Segment>& Y = lineages[b];
    merged.clear();
    std::int32_t parent = -1;
    std::size_t i = 0, j = 0;
    Segment x = X[0], y = Y[0];
    for (;;) {
      if (i == X.size()) {
        for (; j < Y.size(); j++) {
          detail::append(merged, y);
          if (j + 1 < Y.size()) y = Y[j + 1];
        }
        break;
      }
      if (j == Y.size()) {
        for (; i < X.size(); i++) {
          detail::append(merged, x);
          if (i + 1 < X.size()) x = X[i + 1];
        }
        break;
      }
      if (x.right <= y.left) {
        detail::append(merged, x);
        if (++i < X.size()) x = X[i];
      } else if (y.right <= x.left) {
        detail::append(merged, y);
        if (++j < Y.size()) y = Y[j];
      } else if (x.left < y.left) {
        detail::append(merged, {x.left, y.left, x.node});
        x.left = y.left;
      } else if (y.left < x.left) {
        detail::append(merged, {y.left, x.left, y.node});
        y.left = x.left;
      } else {
        // Both carry [l, e): a coalescence there.
        const double l = x.left, e = std::min(x.right, y.right);
        if (parent < 0) {
          parent = static_cast<std::int32_t>(g.time.size());
          g.time.push_back(t);
        }
        g.edges.push_back({l, e, parent, x.node});
        g.edges.push_back({l, e, parent, y.node});
        coverage.merge(l, e, [&](double p, double q) { detail::append(merged, {p, q, parent}); });
        x.left = y.left = e;
        if (x.left == x.right && ++i < X.size()) x = X[i];
        if (y.left == y.right && ++j < Y.size()) y = Y[j];
      }
    }
    // Replace the two by the merged lineage (if it still carries anything), in a fixed order.
    const std::size_t hi = std::max(a, b), lo = std::min(a, b);
    if (hi + 1 != lineages.size()) lineages[hi] = std::move(lineages.back());
    lineages.pop_back();
    if (merged.empty()) {
      if (lo + 1 != lineages.size()) lineages[lo] = std::move(lineages.back());
      lineages.pop_back();
    } else {
      lineages[lo] = merged;
    }
  }

  // Edges sorted and joined where contiguous.
  std::sort(g.edges.begin(), g.edges.end(), [](const Edge& p, const Edge& q) {
    if (p.parent != q.parent) return p.parent < q.parent;
    if (p.child != q.child) return p.child < q.child;
    return p.left < q.left;
  });
  std::vector<Edge> joined;
  joined.reserve(g.edges.size());
  for (const Edge& e : g.edges) {
    if (!joined.empty() && joined.back().parent == e.parent && joined.back().child == e.child && joined.back().right == e.left) joined.back().right = e.right;
    else joined.push_back(e);
  }
  g.edges = std::move(joined);

  // Mutations, edge by edge, then in order of position.
  if (spec.theta > 0) {
    calx::random::Stream m(seed, COALESCENT_MUTATIONS, replicate, 0);
    const double half_theta = spec.theta / 2;
    for (const Edge& e : g.edges) {
      const double t0 = g.time[static_cast<std::size_t>(e.child)], t1 = g.time[static_cast<std::size_t>(e.parent)];
      const std::int64_t k = m.poisson(half_theta * (e.right - e.left) * (t1 - t0));
      for (std::int64_t c = 0; c < k; c++) {
        const double pos = e.left + (e.right - e.left) * m.uniform();
        const double when = t0 + (t1 - t0) * m.uniform();
        g.mutations.push_back({std::min(pos, std::nextafter(e.right, e.left)), when, e.child});
      }
    }
    std::sort(g.mutations.begin(), g.mutations.end(), [](const Mutation& p, const Mutation& q) {
      return p.position != q.position ? p.position < q.position : p.node < q.node;
    });
  }
  return g;
}

// ---- the marginal trees ----------------------------------------------------------------------------------

// One marginal tree during a sweep: parent over all nodes of the genealogy (−1 for nodes not
// in the tree, and for the root).
struct TreeView {
  double left = 0, right = 1;
  std::span<const std::int32_t> parent;
  std::span<const double> time;
  std::size_t samples = 0;
  std::int32_t root = -1;
  double length = 0;  // total branch length
  double tmrca() const { return root < 0 ? 0 : time[static_cast<std::size_t>(root)]; }
};

// Calls visit(view) on each marginal tree from left to right (Kelleher et al. 2016, Algorithm
// T); a visitor returning bool stops the sweep by returning false.
template <class Visit>
void for_each_tree(const Genealogy& g, Visit&& visit) {
  const std::size_t E = g.edges.size(), V = g.nodes();
  std::vector<std::uint32_t> in(E), out(E);
  std::iota(in.begin(), in.end(), 0u);
  std::iota(out.begin(), out.end(), 0u);
  const auto tp = [&](std::uint32_t e) { return g.time[static_cast<std::size_t>(g.edges[e].parent)]; };
  std::sort(in.begin(), in.end(), [&](std::uint32_t p, std::uint32_t q) {
    return g.edges[p].left != g.edges[q].left ? g.edges[p].left < g.edges[q].left : tp(p) < tp(q);
  });
  std::sort(out.begin(), out.end(), [&](std::uint32_t p, std::uint32_t q) {
    return g.edges[p].right != g.edges[q].right ? g.edges[p].right < g.edges[q].right : tp(p) > tp(q);
  });
  std::vector<std::int32_t> parent(V, -1);
  std::vector<std::uint32_t> stamp(V, 0);
  std::uint32_t tree = 0;
  std::size_t j = 0, k = 0;
  double x = 0;
  while (x < 1) {
    while (k < E && g.edges[out[k]].right == x) parent[static_cast<std::size_t>(g.edges[out[k++]].child)] = -1;
    while (j < E && g.edges[in[j]].left == x) {
      const Edge& e = g.edges[in[j++]];
      parent[static_cast<std::size_t>(e.child)] = e.parent;
    }
    double right = 1;
    if (j < E) right = std::min(right, g.edges[in[j]].left);
    if (k < E) right = std::min(right, g.edges[out[k]].right);
    TreeView view;
    view.left = x;
    view.right = right;
    view.parent = parent;
    view.time = g.time;
    view.samples = g.samples;
    // The root, and the length: each branch counted once, walking up from every sample.
    tree++;
    std::int32_t v = 0;
    while (parent[static_cast<std::size_t>(v)] >= 0) v = parent[static_cast<std::size_t>(v)];
    view.root = v;
    double length = 0;
    for (std::size_t s = 0; s < g.samples; s++)
      for (std::int32_t u = static_cast<std::int32_t>(s); parent[static_cast<std::size_t>(u)] >= 0 && stamp[static_cast<std::size_t>(u)] != tree;
           u = parent[static_cast<std::size_t>(u)]) {
        stamp[static_cast<std::size_t>(u)] = tree;
        length += g.time[static_cast<std::size_t>(parent[static_cast<std::size_t>(u)])] - g.time[static_cast<std::size_t>(u)];
      }
    view.length = length;
    if constexpr (std::is_same_v<std::invoke_result_t<Visit&, const TreeView&>, bool>) {
      if (!visit(static_cast<const TreeView&>(view))) return;
    } else {
      visit(static_cast<const TreeView&>(view));
    }
    x = right;
  }
}

inline std::size_t tree_count(const Genealogy& g) {
  std::size_t c = 0;
  for_each_tree(g, [&](const TreeView&) { c++; });
  return c;
}

// The time to the most recent common ancestor at position x of the chromosome.
inline double tmrca_at(const Genealogy& g, double x) {
  double t = 0;
  for_each_tree(g, [&](const TreeView& v) {
    if (x < v.right) {
      t = v.tmrca();
      return false;
    }
    return true;
  });
  return t;
}

// The haplotypes of the sample: one site for each mutation, 1 below it.
inline Haplotypes haplotypes(const Genealogy& g) {
  Haplotypes h;
  h.sequences = g.samples;
  h.sites = g.mutations.size();
  h.alleles.assign(h.sites * h.sequences, 0);
  h.positions.resize(h.sites);
  for (std::size_t m = 0; m < h.sites; m++) h.positions[m] = g.mutations[m].position;
  const std::size_t V = g.nodes();
  std::vector<std::int32_t> head(V, -1), next(h.sites, -1);
  std::vector<std::uint32_t> mark(V, 0);
  std::uint32_t tree = 0;
  std::size_t m = 0;
  for_each_tree(g, [&](const TreeView& v) {
    tree++;
    const std::size_t first = m;
    while (m < h.sites && g.mutations[m].position < v.right) {
      const auto node = static_cast<std::size_t>(g.mutations[m].node);
      if (mark[node] != tree) {
        mark[node] = tree;
        head[node] = -1;
      }
      next[m] = head[node];
      head[node] = static_cast<std::int32_t>(m);
      m++;
    }
    if (m == first) return;
    for (std::size_t s = 0; s < g.samples; s++)
      for (std::int32_t u = static_cast<std::int32_t>(s); u >= 0; u = v.parent[static_cast<std::size_t>(u)])
        if (mark[static_cast<std::size_t>(u)] == tree)
          for (std::int32_t q = head[static_cast<std::size_t>(u)]; q >= 0; q = next[static_cast<std::size_t>(q)])
            h.alleles[static_cast<std::size_t>(q) * g.samples + s] = 1;
  });
  return h;
}

// A marginal tree for drawing: its own nodes, the samples first (node i is sample i), then
// the ancestors by time; children in order of index; the mutations on its branches.
struct TreeMutation {
  double position = 0, time = 0;
  std::int32_t node = -1;  // local
};
struct Tree {
  double left = 0, right = 1;
  std::vector<std::int32_t> node;  // the genealogy's node of each local node
  std::vector<double> time;
  std::vector<std::int32_t> parent;          // local; −1 at the root
  std::vector<std::uint32_t> child_begin;    // children of v: children[child_begin[v] … child_begin[v + 1])
  std::vector<std::int32_t> children;
  std::int32_t root = -1;
  double tmrca = 0, length = 0;
  std::vector<TreeMutation> mutations;
};

namespace detail {

inline Tree make_tree(const Genealogy& g, const TreeView& v) {
  Tree t;
  t.left = v.left;
  t.right = v.right;
  t.tmrca = v.tmrca();
  t.length = v.length;
  std::vector<std::int32_t> local(g.nodes(), -1);
  for (std::size_t s = 0; s < g.samples; s++) {
    local[s] = static_cast<std::int32_t>(s);
    t.node.push_back(static_cast<std::int32_t>(s));
  }
  std::vector<std::int32_t> ancestors;
  for (std::size_t s = 0; s < g.samples; s++)
    for (std::int32_t u = v.parent[s]; u >= 0 && local[static_cast<std::size_t>(u)] == -1; u = v.parent[static_cast<std::size_t>(u)]) {
      local[static_cast<std::size_t>(u)] = -2;
      ancestors.push_back(u);
    }
  std::sort(ancestors.begin(), ancestors.end());  // node ids increase with time
  for (const std::int32_t u : ancestors) {
    local[static_cast<std::size_t>(u)] = static_cast<std::int32_t>(t.node.size());
    t.node.push_back(u);
  }
  const std::size_t k = t.node.size();
  t.time.resize(k);
  t.parent.resize(k);
  std::vector<std::uint32_t> degree(k + 1, 0);
  for (std::size_t i = 0; i < k; i++) {
    const auto u = static_cast<std::size_t>(t.node[i]);
    t.time[i] = g.time[u];
    t.parent[i] = v.parent[u] < 0 ? -1 : local[static_cast<std::size_t>(v.parent[u])];
    if (t.parent[i] < 0) t.root = static_cast<std::int32_t>(i);
    else degree[static_cast<std::size_t>(t.parent[i]) + 1]++;
  }
  t.child_begin.assign(k + 1, 0);
  for (std::size_t i = 0; i < k; i++) t.child_begin[i + 1] = t.child_begin[i] + degree[i + 1];
  t.children.resize(k ? k - 1 : 0);
  std::vector<std::uint32_t> fill(t.child_begin.begin(), t.child_begin.end() - 1);
  for (std::size_t i = 0; i < k; i++)
    if (t.parent[i] >= 0) t.children[fill[static_cast<std::size_t>(t.parent[i])]++] = static_cast<std::int32_t>(i);
  for (const Mutation& m : g.mutations)
    if (m.position >= v.left && m.position < v.right) t.mutations.push_back({m.position, m.time, local[static_cast<std::size_t>(m.node)]});
  return t;
}

}  // namespace detail

// Every marginal tree, left to right (at most `limit` of them).
inline std::vector<Tree> marginal_trees(const Genealogy& g, std::size_t limit = std::numeric_limits<std::size_t>::max()) {
  std::vector<Tree> out;
  for_each_tree(g, [&](const TreeView& v) {
    if (out.size() >= limit) return false;
    out.push_back(detail::make_tree(g, v));
    return true;
  });
  return out;
}

// The marginal tree at position x.
inline Tree tree_at(const Genealogy& g, double x) {
  Tree t;
  for_each_tree(g, [&](const TreeView& v) {
    if (x < v.right) {
      t = detail::make_tree(g, v);
      return false;
    }
    return true;
  });
  return t;
}

// summary(genealogy) for replicates 0 … count − 1 of a seed, in parallel; out[i] is replicate
// i's, the same on any number of threads.
template <class Summary>
auto coalescent_replicates(const CoalescentSpec& spec, std::uint64_t seed, std::size_t count, calx::scheduler::ThreadPool& pool, Summary&& summary) {
  using Result = std::invoke_result_t<Summary&, const Genealogy&>;
  std::vector<Result> out(count);
  calx::scheduler::parallel_for(pool, count, 16, [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t i = begin; i < end; i++) out[i] = summary(static_cast<const Genealogy&>(coalescent(spec, seed, static_cast<std::uint32_t>(i))));
  });
  return out;
}

}  // namespace evolution::population
