// The likelihood of a tree: Felsenstein's pruning algorithm, branch lengths by Newton's
// method, a search of topologies by nearest-neighbour interchanges, and the parameters of the
// model by one-dimensional maximisation; the maximum-likelihood tree.
//
// Pruning (Felsenstein 1981, "Evolutionary trees from DNA sequences: a maximum likelihood
// approach", J. Mol. Evol. 17, 368–376): the probability of a column of states at the leaves is
// a sum over the states of the internal nodes of products along the branches, which the
// distributive law turns into a product of messages passed from the leaves to the root,
//
//   L_v(x) = Π_{c child of v} Σ_y P_xy(t_c) L_c(y),   L = Σ_x π_x L_root(x),
//
// for each site pattern and each rate category (the site's likelihood the weighted sum over
// categories). It is the sum–product algorithm on the tree (theory/phylogeny.md, §5); O(n · s ·
// K · 16) for n nodes, s patterns and K categories. Partial likelihoods are rescaled by 2^256,
// pattern by pattern, whenever their largest entry falls below 2^−256, and the number of
// rescalings is carried, so that no likelihood underflows on any tree.
//
// Reversibility makes the likelihood independent of where the root is (Felsenstein 1981, the
// pulley principle): for a branch v below p, L = Σ_x π_x W_v(x) Σ_y P_xy(t_v) L_v(y), with W_v the
// conditional likelihood at p of everything outside v's subtree. Every branch is then a root,
// and its length is fitted with the rest fixed: with the spectral form of P (models.hpp),
//
//   L_site(t) = Σ_c w_c Σ_k e^{λ_k r_c t} G_ck,  G_ck = (Σ_x π_x W(x) A_xk)(Σ_y B_ky L_v(y)),
//
// so that after O(s K 32) to form G, each Newton step on the log-likelihood with its first and
// second derivatives costs O(s K 4) and no logarithm (the "sum" buffers of Guindon and Gascuel
// 2003, Syst. Biol. 52, 696–704). Newton's method is applied to the derivative, inside a bracket
// of its change of sign, with bisection where a step leaves it (Press et al. 2007, Numerical
// Recipes, §9.4), the length kept in [1e−8, 50]. A sweep visits the branches in preorder,
// keeping the partials from both sides current, so that every branch is fitted on the
// likelihood of the tree as it then is. Only the partials of internal nodes are stored, from
// the leaves (down) and from the rest of the tree (outer); a leaf is read from its states. The
// kernels hold the four states of a pattern and category in calx::cpu::Lanes<double, 4> (a
// message is Σ_y column_y(P) · L(y)); about 5 ns per internal node, pattern and category on
// the development machine (AVX-512, one thread).
//
// Nearest-neighbour interchanges: an internal branch separates four subtrees, (A, B | C, D);
// the two other trees on them, (A, C | B, D) and (A, D | B, C), are scored with the five
// branches around it refitted and the rest fixed, from the conditional likelihoods of the four
// subtrees (Guindon and Gascuel 2003, §"NNI"). The interchanges that gain, on branches far
// enough apart, are made together and the tree refitted; if that does not gain, only the best
// is made, which always gains. The search climbs until no interchange gains more than the
// tolerance: a local maximum, from a neighbour-joining start (Saitou and Nei 1987).
//
// Model parameters (κ, κ₁ and κ₂, the five free rates of GTR, the Γ shape α, the proportion of
// invariant sites) are fitted one at a time, the branch lengths refitted for each value, by
// Brent's method on a logarithmic scale (Brent 1973, ch. 5; calx::optimisation); base
// frequencies are the empirical ones. Limits: nucleotides (four states); a local optimum of the topology, not
// the global one; the topology search assumes a binary tree.

#pragma once

#include <evolution/phylogeny/alignment.hpp>
#include <evolution/phylogeny/distance_trees.hpp>
#include <evolution/phylogeny/models.hpp>
#include <evolution/phylogeny/tree.hpp>

#include <calx/cpu/lanes.hpp>
#include <calx/optimisation/minimise.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

namespace evolution::phylogeny {

namespace detail {

inline constexpr double SCALE_UP = 0x1p256, SCALE_LIMIT = 0x1p-256;
inline constexpr double LOG_SCALE = -256 * std::numbers::ln2;  // the log of one rescaling's factor

// The maximum of f on [a, b] to a width `tol`, by Brent's method (calx::optimisation): golden
// section where the log-likelihood is rough, parabolas where it is smooth near its maximum.
template <class F>
double maximum_on(F&& f, double a, double b, double tol, double& best) {
  const auto m = calx::optimisation::brent_maximise(f, a, b, tol);
  best = m.f;
  return m.x;
}

}  // namespace detail

struct BranchFit {
  double length = 0, log_likelihood = 0;
};

class Likelihood {
 public:
  double min_length = 1e-8, max_length = 50;

  Likelihood(Patterns data, Model model, Rates rates) : data_(std::move(data)), model_(std::move(model)), rates_(std::move(rates)) { size(); }

  const Patterns& data() const { return data_; }
  const Model& model() const { return model_; }
  const Rates& rates() const { return rates_; }
  void set_model(const Model& m) { model_ = m; }
  void set_rates(const Rates& r) { rates_ = r; size(); }
  void set_weights(const std::vector<double>& w) {
    if (w.size() != data_.count) throw std::invalid_argument("evolution::phylogeny::Likelihood::set_weights: a weight for each pattern");
    data_.weight = w;
  }

  // The log-likelihood of the patterns on a tree, rooted anywhere, of any degree.
  double log_likelihood(const Tree& t) {
    prepare(t);
    return root_log_likelihood(t);
  }

  // Fits every branch length, sweeping until a sweep gains less than `tolerance`; returns the
  // log-likelihood. A tree whose root has two children and that has three leaves or more is
  // unrooted first (a reversible model sees only the sum of the two branches), and a tree of two
  // leaves has its first branch held at zero; absent lengths start at 0.1.
  double optimise_branch_lengths(Tree& t, int sweeps = 30, double tolerance = 1e-8) {
    normalise(t);
    prepare(t);
    double L = root_log_likelihood(t);
    for (int s = 0; s < sweeps; s++) {
      const double next = sweep(t);
      const bool done = next - L < tolerance;
      L = std::max(L, next);
      if (done) break;
    }
    return L;
  }

  // Hill-climbing by nearest-neighbour interchanges from t, its branch lengths fitted on the way;
  // returns the log-likelihood of the tree it stops at, which t becomes. `moves` counts the
  // interchanges made.
  double nni_search(Tree& t, int rounds = 100, double tolerance = 1e-6, int* moves = nullptr, int sweeps = 30) {
    double L = optimise_branch_lengths(t, sweeps, tolerance);
    if (!t.is_binary() || t.leaf_count() < 4) return L;
    int made = 0;
    for (int r = 0; r < rounds; r++) {
      refresh(t);
      auto cands = candidates(t, L, tolerance);
      if (cands.empty()) break;
      std::stable_sort(cands.begin(), cands.end(), [](const Candidate& a, const Candidate& b) { return a.gain > b.gain; });
      // Those far enough apart, greedily by gain.
      std::vector<Candidate> chosen;
      std::vector<char> used(t.size(), 0);
      for (const auto& c : cands) {
        bool free = true;
        for (int x : c.nodes)
          if (x >= 0 && used[static_cast<std::size_t>(x)]) free = false;
        if (!free) continue;
        for (int x : c.nodes)
          if (x >= 0) used[static_cast<std::size_t>(x)] = 1;
        chosen.push_back(c);
      }
      Tree trial = t;
      for (const auto& c : chosen) apply(trial, c);
      double Lt = optimise_branch_lengths(trial, sweeps, tolerance);
      int now = static_cast<int>(chosen.size());
      if (chosen.size() > 1 && Lt < L + cands[0].gain * 0.5) {
        trial = t;
        apply(trial, cands[0]);
        Lt = optimise_branch_lengths(trial, sweeps, tolerance);
        now = 1;
      }
      if (!(Lt > L + tolerance)) break;
      t = std::move(trial);
      L = Lt;
      made += now;
    }
    if (moves) *moves += made;
    return L;
  }

 private:
  // A partial as the kernels read it: a leaf's states, or an array of conditional likelihoods
  // (patterns × categories × 4) with its rescalings per pattern.
  using V4 = calx::cpu::Lanes<double, 4>;

  struct Source {
    const Code* tip = nullptr;
    const double* partial = nullptr;
    const int* scale = nullptr;
  };

  Patterns data_;
  Model model_;
  Rates rates_;
  std::size_t P_ = 0, K_ = 0, B_ = 0;
  std::size_t nodes_ = 0, internal_ = 0;
  std::vector<int> slot_;                      // each internal node's place in down_ and outer_, −1 for a leaf
  std::vector<double> down_, outer_;           // internal nodes × B
  std::vector<int> down_scale_, outer_scale_;  // internal nodes × P
  std::vector<Code> tips_;                     // nodes × P: each leaf's states
  std::vector<double> pm_;                     // transition matrices, K × 16 per slot
  std::vector<double> G_, coef_, tipB_;
  std::vector<double> buf_[4];                 // X, Y, W and a leaf's outer partial
  std::vector<int> sbuf_[4];

  struct Candidate {
    int v = -1, swap = 0;
    double gain = 0;
    std::array<double, 5> lengths{};  // q1, q2, q3, q4, the middle
    std::array<int, 6> nodes{};       // the nodes the interchange touches
  };

  void size() {
    P_ = data_.count;
    K_ = rates_.size();
    B_ = P_ * K_ * 4;
    G_.assign(B_, 0);
    coef_.assign(K_ * 12, 0);
    for (auto& b : buf_) b.assign(B_, 0);
    for (auto& x : sbuf_) x.assign(P_, 0);
    nodes_ = 0;
  }

  double* down(int v) { return down_.data() + static_cast<std::size_t>(slot_[static_cast<std::size_t>(v)]) * B_; }
  int* down_scale(int v) { return down_scale_.data() + static_cast<std::size_t>(slot_[static_cast<std::size_t>(v)]) * P_; }
  // The outer partial of v: kept for internal nodes, in scratch for a leaf (used at once).
  double* outer(int v) {
    const int k = slot_[static_cast<std::size_t>(v)];
    return k < 0 ? buf_[3].data() : outer_.data() + static_cast<std::size_t>(k) * B_;
  }
  int* outer_scale(int v) {
    const int k = slot_[static_cast<std::size_t>(v)];
    return k < 0 ? sbuf_[3].data() : outer_scale_.data() + static_cast<std::size_t>(k) * P_;
  }
  Source source(int v) {
    if (slot_[static_cast<std::size_t>(v)] < 0) return {tips_.data() + static_cast<std::size_t>(v) * P_, nullptr, nullptr};
    return {nullptr, down(v), down_scale(v)};
  }

  void normalise(Tree& t) {
    if (t.root >= 0 && t[t.root].children.size() == 2 && t.leaf_count() >= 3) t = unroot(t);
    for (std::size_t v = 0; v < t.size(); v++) {
      if (static_cast<int>(v) == t.root) continue;
      double& l = t.nodes[v].length;
      l = std::isnan(l) ? 0.1 : std::clamp(l, min_length, max_length);
    }
    // A root of degree two that remains (two leaves): its first branch is held at zero and the
    // second carries the whole distance.
    if (t.root >= 0 && t[t.root].children.size() == 2) {
      const int a = t[t.root].children[0], b = t[t.root].children[1];
      t[b].length = std::clamp(t[a].length + t[b].length, min_length, max_length);
      t[a].length = 0;
    }
  }

  // Sizes the buffers for a tree and copies its leaves' states.
  void bind(const Tree& t) {
    std::size_t internal = 0;
    for (const auto& n : t.nodes) internal += !n.children.empty();
    if (t.size() != nodes_ || internal != internal_) {
      nodes_ = t.size();
      internal_ = internal;
      down_.assign(internal_ * B_, 0);
      outer_.assign(internal_ * B_, 0);
      down_scale_.assign(internal_ * P_, 0);
      outer_scale_.assign(internal_ * P_, 0);
      tips_.assign(nodes_ * P_, MISSING);
    }
    slot_.assign(nodes_, -1);
    int k = 0;
    for (std::size_t v = 0; v < nodes_; v++) {
      if (!t.nodes[v].children.empty()) {
        slot_[v] = k++;
        continue;
      }
      const std::size_t r = data_.index(t.nodes[v].label);
      if (r == SIZE_MAX) throw std::invalid_argument("evolution::phylogeny::Likelihood: the leaf " + t.nodes[v].label + " has no sequence");
      Code* x = &tips_[v * P_];
      for (std::size_t p = 0; p < P_; p++) x[p] = data_.states[p * data_.taxa + r];
    }
    // B's columns for a leaf's state; the last row for a missing state, Σ_y B_ky.
    tipB_.assign(20, 0);
    for (int j = 0; j < 4; j++)
      for (int y = 0; y < 4; y++) {
        tipB_[static_cast<std::size_t>(y * 4 + j)] = model_.right[static_cast<std::size_t>(j * 4 + y)];
        tipB_[static_cast<std::size_t>(16 + j)] += model_.right[static_cast<std::size_t>(j * 4 + y)];
      }
  }

  // P(r_c t) for every category, into matrix slot `slot`, by columns: column y (the
  // probabilities of reaching y from each x) is the four values from 4y.
  const double* matrices(double t, std::size_t slot) {
    const std::size_t need = (slot + 1) * K_ * 16;
    if (pm_.size() < need) pm_.resize(need);
    double* M = &pm_[slot * K_ * 16];
    double rows[16];
    for (std::size_t c = 0; c < K_; c++) {
      model_.transition(rates_.rate[c] * t, rows);
      for (int x = 0; x < 4; x++)
        for (int y = 0; y < 4; y++) M[c * 16 + static_cast<std::size_t>(y * 4 + x)] = rows[x * 4 + y];
    }
    return M;
  }

  // dst = Π_i P_i · src_i, the messages of the sources through their matrices, rescaled, with
  // its rescalings: one pass over the patterns, four states to a vector of lanes (a message is
  // Σ_y column_y · src(y)). Matrices must not move during the call.
  void combine(double* dst, int* dscale, const Source* src, const double* const* M, std::size_t n) {
    const std::size_t w = K_ * 4;
    const V4 one = V4::splat(1.0);
    for (std::size_t p = 0; p < P_; p++) {
      double* out = dst + p * w;
      V4 top = V4::zero();
      for (std::size_t c = 0; c < K_; c++) {
        V4 acc = one;
        for (std::size_t i = 0; i < n; i++) {
          const double* m = M[i] + c * 16;
          if (src[i].tip) {
            const Code x = src[i].tip[p];
            if (x < 4) acc = acc * V4::load(m + 4 * x);
          } else {
            const double* d = src[i].partial + p * w + c * 4;
            acc = acc * (V4::load(m) * V4::splat(d[0]) + V4::load(m + 4) * V4::splat(d[1]) + V4::load(m + 8) * V4::splat(d[2]) +
                         V4::load(m + 12) * V4::splat(d[3]));
          }
        }
        acc.store(out + c * 4);
        top = max(top, acc);
      }
      const double largest = std::max(std::max(top[0], top[1]), std::max(top[2], top[3]));
      int sc = 0;
      for (std::size_t i = 0; i < n; i++)
        if (src[i].scale) sc += src[i].scale[p];
      if (largest > 0 && largest < detail::SCALE_LIMIT) {
        for (std::size_t i = 0; i < w; i++) out[i] *= detail::SCALE_UP;
        sc++;
      }
      dscale[p] = sc;
    }
  }

  // The matrices of every slot are computed first, then pointers taken, since pm_ may grow.
  void compute_down(const Tree& t, int v) {
    const auto& ch = t[v].children;
    if (ch.empty()) return;
    src_.resize(ch.size());
    mats_.resize(ch.size());
    for (std::size_t k = 0; k < ch.size(); k++) (void)matrices(t[ch[k]].length, k);
    for (std::size_t k = 0; k < ch.size(); k++) {
      src_[k] = source(ch[k]);
      mats_[k] = &pm_[k * K_ * 16];
    }
    combine(down(v), down_scale(v), src_.data(), mats_.data(), ch.size());
  }
  // W_c: at c's parent v, everything outside c's subtree: v's own outer partial through v's
  // branch (unless v is the root), times the messages of c's siblings.
  void compute_outer(const Tree& t, int c) {
    const int v = t[c].parent;
    src_.clear();
    std::vector<double> lengths;
    if (v != t.root) {
      src_.push_back({nullptr, outer(v), outer_scale(v)});
      lengths.push_back(t[v].length);
    }
    for (int sib : t[v].children) {
      if (sib == c) continue;
      src_.push_back(source(sib));
      lengths.push_back(t[sib].length);
    }
    for (std::size_t k = 0; k < lengths.size(); k++) (void)matrices(lengths[k], k);
    mats_.resize(lengths.size());
    for (std::size_t k = 0; k < lengths.size(); k++) mats_[k] = &pm_[k * K_ * 16];
    combine(outer(c), outer_scale(c), src_.data(), mats_.data(), src_.size());
  }
  std::vector<Source> src_;
  std::vector<const double*> mats_;

  // Partials from the leaves up, for the tree as it is.
  void prepare(const Tree& t) {
    if (t.root < 0 || t.is_leaf(t.root)) throw std::invalid_argument("evolution::phylogeny::Likelihood: a tree of two leaves or more");
    bind(t);
    for (int v : t.postorder()) compute_down(t, v);
  }
  // And from the root down: every internal node's outer partial.
  void refresh(const Tree& t) {
    prepare(t);
    for (int v : t.preorder())
      for (int c : t[v].children)
        if (!t.is_leaf(c)) compute_outer(t, c);
  }

  double root_log_likelihood(const Tree& t) {
    const double* d = down(t.root);
    const int* s = down_scale(t.root);
    double L = 0;
    for (std::size_t p = 0; p < P_; p++) {
      double site = 0;
      for (std::size_t c = 0; c < K_; c++) {
        const double* x = d + (p * K_ + c) * 4;
        site += rates_.weight[c] * (model_.pi[0] * x[0] + model_.pi[1] * x[1] + model_.pi[2] * x[2] + model_.pi[3] * x[3]);
      }
      L += data_.weight[p] * (std::log(std::max(site, std::numeric_limits<double>::min())) + s[p] * detail::LOG_SCALE);
    }
    return L;
  }

  // The best length of the branch between a conditional W (at one end) and the source D (at
  // the other), by Newton's method from t0.
  BranchFit fit_branch(const double* W, const int* sw, const Source& D, double t0) {
    // G_ck = (Σ_x π_x W(x) A_xk) (Σ_y B_ky D(y)), lanes over k: rows of A, columns of B.
    const V4 pi = V4::load(model_.pi.data());
    const V4 A0 = V4::load(&model_.left[0]), A1 = V4::load(&model_.left[4]), A2 = V4::load(&model_.left[8]), A3 = V4::load(&model_.left[12]);
    const V4 B0 = V4::load(&tipB_[0]), B1 = V4::load(&tipB_[4]), B2 = V4::load(&tipB_[8]), B3 = V4::load(&tipB_[12]);
    double constant = 0;
    for (std::size_t p = 0; p < P_; p++) {
      constant += data_.weight[p] * (sw[p] + (D.scale ? D.scale[p] : 0)) * detail::LOG_SCALE;
      const double* tb = D.tip ? &tipB_[static_cast<std::size_t>(std::min<Code>(D.tip[p], 4)) * 4] : nullptr;
      for (std::size_t c = 0; c < K_; c++) {
        const V4 pw = pi * V4::load(W + (p * K_ + c) * 4);
        const V4 left = A0 * V4::splat(pw[0]) + A1 * V4::splat(pw[1]) + A2 * V4::splat(pw[2]) + A3 * V4::splat(pw[3]);
        V4 right;
        if (tb) {
          right = V4::load(tb);
        } else {
          const double* d = D.partial + (p * K_ + c) * 4;
          right = B0 * V4::splat(d[0]) + B1 * V4::splat(d[1]) + B2 * V4::splat(d[2]) + B3 * V4::splat(d[3]);
        }
        (left * right).store(&G_[(p * K_ + c) * 4]);
      }
    }
    // The coefficients w_c e^{λ_k r_c t} and their products with λ_k r_c and its square, by
    // category and eigenvalue, in three rows.
    const std::size_t n = K_ * 4;
    double* e0 = coef_.data();
    double* e1 = e0 + n;
    double* e2 = e1 + n;
    const auto coefficients = [&](double t) {
      for (std::size_t c = 0; c < K_; c++)
        for (int k = 0; k < 4; k++) {
          const std::size_t i = c * 4 + static_cast<std::size_t>(k);
          const double rate = model_.lambda[static_cast<std::size_t>(k)] * rates_.rate[c];
          e0[i] = rates_.weight[c] * std::exp(rate * t);
          e1[i] = e0[i] * rate;
          e2[i] = e1[i] * rate;
        }
    };
    const auto sum = [](V4 x) { return (x[0] + x[1]) + (x[2] + x[3]); };
    // The first two derivatives of the log-likelihood in t (no logarithms needed).
    const auto derivatives = [&](double t, double& f1, double& f2) {
      coefficients(t);
      f1 = 0;
      f2 = 0;
      for (std::size_t p = 0; p < P_; p++) {
        const double* g = &G_[p * n];
        V4 l0 = V4::zero(), l1 = V4::zero(), l2 = V4::zero();
        for (std::size_t i = 0; i < n; i += 4) {
          const V4 gi = V4::load(g + i);
          l0 = l0 + V4::load(e0 + i) * gi;
          l1 = l1 + V4::load(e1 + i) * gi;
          l2 = l2 + V4::load(e2 + i) * gi;
        }
        const double L = std::max(sum(l0), std::numeric_limits<double>::min());
        const double r1 = sum(l1) / L, w = data_.weight[p];
        f1 += w * r1;
        f2 += w * (sum(l2) / L - r1 * r1);
      }
    };
    const auto value = [&](double t) {
      coefficients(t);
      double f = constant;
      for (std::size_t p = 0; p < P_; p++) {
        const double* g = &G_[p * n];
        V4 l0 = V4::zero();
        for (std::size_t i = 0; i < n; i += 4) l0 = l0 + V4::load(e0 + i) * V4::load(g + i);
        f += data_.weight[p] * std::log(std::max(sum(l0), std::numeric_limits<double>::min()));
      }
      return f;
    };
    // Newton's method on the derivative, kept inside a bracket [lo, hi] of its sign change
    // (Press et al. 2007, Numerical Recipes, §9.4, rtsafe): a step that leaves the bracket, or
    // where the curvature is not negative, is replaced by bisection, or by doubling while no
    // upper end is known.
    double t = std::clamp(t0, min_length, max_length), f1 = 0, f2 = 0;
    double lo = min_length, hi = max_length;
    bool hi_known = false;
    derivatives(t, f1, f2);
    for (int it = 0; it < 100; it++) {
      if (f1 > 0) lo = t;
      else {
        hi = t;
        hi_known = true;
      }
      if (hi - lo <= 1e-10 * std::max(0.01, t)) break;
      double tn = f2 < 0 ? t - f1 / f2 : NAN;
      if (!(tn > lo && tn < hi)) tn = hi_known ? (lo + hi) / 2 : std::min(max_length, 2 * t + 0.01);
      if (!hi_known && tn == max_length) hi_known = true;
      const double moved = std::fabs(tn - t);
      t = tn;
      derivatives(t, f1, f2);
      if (moved < 1e-8 * std::max(0.01, t)) break;
    }
    if (f1 > 0 && t >= max_length) t = max_length;
    return {t, value(t)};
  }

  // A branch held fixed: the first of two at a root of degree two.
  static bool fixed(const Tree& t, int v) {
    const int p = t[v].parent;
    return p == t.root && t[p].children.size() == 2 && t[p].children[0] == v;
  }

  // Fits v's branch, then its subtree's, keeping the partials from both sides current.
  void visit(Tree& t, int v) {
    if (!fixed(t, v)) t[v].length = fit_branch(outer(v), outer_scale(v), source(v), t[v].length).length;
    if (t.is_leaf(v)) return;
    for (int c : t[v].children) {
      compute_outer(t, c);
      visit(t, c);
    }
    compute_down(t, v);
  }

  // One sweep over every branch; the partials from the leaves must be current.
  double sweep(Tree& t) {
    for (int c : t[t.root].children) {
      compute_outer(t, c);
      visit(t, c);
    }
    compute_down(t, t.root);
    return root_log_likelihood(t);
  }

  // ---- nearest-neighbour interchanges ----------------------------------------------------------

  struct Part {
    Source source;
    double length = 0;
  };

  // The quartet (q1, q2 | q3, q4) with its five branches refitted: two rounds round the five.
  std::pair<double, std::array<double, 5>> quartet(const std::array<Part, 4>& q, double middle) {
    double* X = buf_[0].data();
    double* Y = buf_[1].data();
    double* W = buf_[2].data();
    int* sX = sbuf_[0].data();
    int* sY = sbuf_[1].data();
    int* sW = sbuf_[2].data();
    std::array<double, 5> len{q[0].length, q[1].length, q[2].length, q[3].length, middle};
    const auto two = [&](double* out, int* so, const Source& a, double la, const Source& b, double lb) {
      (void)matrices(la, 0);
      (void)matrices(lb, 1);
      const Source src[2] = {a, b};
      const double* M[2] = {&pm_[0], &pm_[K_ * 16]};
      combine(out, so, src, M, 2);
    };
    two(X, sX, q[0].source, len[0], q[1].source, len[1]);
    two(Y, sY, q[2].source, len[2], q[3].source, len[3]);
    double f = 0;
    for (int round = 0; round < 2; round++) {
      auto fit = fit_branch(X, sX, Source{nullptr, Y, sY}, len[4]);
      len[4] = fit.length;
      f = fit.log_likelihood;
      for (std::size_t side = 0; side < 2; side++) {
        const Source other = side == 0 ? Source{nullptr, Y, sY} : Source{nullptr, X, sX};
        for (std::size_t j = 0; j < 2; j++) {
          const std::size_t i = side * 2 + j, sib = side * 2 + (1 - j);
          two(W, sW, other, len[4], q[sib].source, len[sib]);
          fit = fit_branch(W, sW, q[i].source, len[i]);
          len[i] = fit.length;
          f = fit.log_likelihood;
        }
        if (side == 0) two(X, sX, q[0].source, len[0], q[1].source, len[1]);
        else two(Y, sY, q[2].source, len[2], q[3].source, len[3]);
      }
    }
    return {f, len};
  }

  // The interchanges that gain, with the outer partials current.
  std::vector<Candidate> candidates(Tree& t, double L, double tolerance) {
    std::vector<Candidate> out;
    for (std::size_t vi = 0; vi < t.size(); vi++) {
      const int v = static_cast<int>(vi);
      if (v == t.root || t.is_leaf(v)) continue;
      const int p = t[v].parent;
      const int a = t[v].children[0], b = t[v].children[1];
      std::vector<int> sibs;
      for (int x : t[p].children)
        if (x != v) sibs.push_back(x);
      const int c = sibs[0];
      const Part pa{source(a), t[a].length}, pb{source(b), t[b].length}, pc{source(c), t[c].length};
      Part px;
      int xnode;
      if (p == t.root) {
        if (sibs.size() != 2) continue;
        xnode = sibs[1];
        px = {source(xnode), t[xnode].length};
      } else {
        xnode = t[p].parent;
        px = {Source{nullptr, outer(p), outer_scale(p)}, t[p].length};
      }
      Candidate best;
      for (int swap = 0; swap < 2; swap++) {
        // swap 0: b ↔ c gives (a, c | b, x); swap 1: a ↔ c gives (c, b | a, x).
        const std::array<Part, 4> q = swap == 0 ? std::array<Part, 4>{pa, pc, pb, px} : std::array<Part, 4>{pc, pb, pa, px};
        const auto [f, len] = quartet(q, t[v].length);
        if (best.v < 0 || f - L > best.gain) {
          best.v = v;
          best.swap = swap;
          best.gain = f - L;
          best.lengths = len;
          best.nodes = {v, p, a, b, c, xnode};
        }
      }
      if (best.gain > tolerance) out.push_back(best);
    }
    return out;
  }

  static void replace(std::vector<int>& xs, int from, int to) {
    for (int& x : xs)
      if (x == from) x = to;
  }
  static void apply(Tree& t, const Candidate& k) {
    const int v = k.v, p = t[v].parent;
    const int a = t[v].children[0], b = t[v].children[1];
    int c = -1, x = -1;
    for (int s : t[p].children)
      if (s != v) {
        if (c < 0) c = s;
        else x = s;
      }
    const int moved = k.swap == 0 ? b : a;  // the child of v exchanged with c
    replace(t[v].children, moved, c);
    replace(t[p].children, c, moved);
    t[c].parent = v;
    t[moved].parent = p;
    if (k.swap == 0) {
      t[a].length = k.lengths[0];
      t[c].length = k.lengths[1];
      t[b].length = k.lengths[2];
    } else {
      t[c].length = k.lengths[0];
      t[b].length = k.lengths[1];
      t[a].length = k.lengths[2];
    }
    if (x >= 0) t[x].length = k.lengths[3];
    else t[p].length = k.lengths[3];
    t[v].length = k.lengths[4];
  }
};

// ---- the maximum-likelihood tree ------------------------------------------------------------

struct SearchOptions {
  bool topology = true;          // search by nearest-neighbour interchanges
  bool model_parameters = true;  // fit κ (κ₁, κ₂, the GTR rates), α and the invariant proportion
  int rounds = 100;              // of interchanges, at most
  int cycles = 6;                // of model parameters and topology, at most
  double tolerance = 1e-6;       // in log-likelihood: the least gain of an interchange or a sweep
  double cycle_tolerance = 1e-3; // the least gain of a cycle of parameters and topology
};

struct MLTree {
  Tree tree;
  double log_likelihood = 0;
  Model model;
  RateSpec rates;
  int interchanges = 0;
};

namespace detail {

// The free exchangeabilities of a model, by their index.
inline std::vector<std::vector<int>> free_rates(ModelKind kind) {
  switch (kind) {
    case ModelKind::k80:
    case ModelKind::hky85: return {{AG, CT}};
    case ModelKind::tn93: return {{AG}, {CT}};
    case ModelKind::gtr: return {{AC}, {AG}, {AT}, {CG}, {CT}};
    default: return {};
  }
}

}  // namespace detail

// The tree, branch lengths and model parameters of greatest likelihood that hill-climbing finds
// from `start`. Models with free frequencies take the data's.
inline MLTree maximum_likelihood(const Patterns& data, Model model, RateSpec rates, const Tree& start, const SearchOptions& o = {}) {
  if (model.kind != ModelKind::jc69 && model.kind != ModelKind::k80) {
    model.pi = base_frequencies(data);
    model.prepare();
  }
  Likelihood lik(data, model, rates.rates());
  MLTree out;
  out.tree = start;
  double L = lik.optimise_branch_lengths(out.tree, 30, o.tolerance);
  // One parameter: f(x) installs the value x and returns the log-likelihood with the branch
  // lengths refitted (two sweeps from the current ones); the value kept is the search's
  // maximum, or the current value if that is better.
  // After the first cycle the search is narrowed to within a factor e of the current value
  // (on the logarithmic scales), and widened back to the bounds if it ends on an edge.
  int cycle = 0;
  const auto fit = [&](auto&& f, double lo, double hi, double current, double reach) {
    const double at_current = f(current);
    double best = -INFINITY;
    double a = lo, b = hi;
    if (cycle > 0) {
      a = std::max(lo, current - reach);
      b = std::min(hi, current + reach);
    }
    double x = detail::maximum_on(f, a, b, 2e-3, best);
    if ((x - a < 4e-3 && a > lo) || (b - x < 4e-3 && b < hi)) x = detail::maximum_on(f, lo, hi, 2e-3, best);
    return best > at_current ? x : current;
  };
  const auto refit = [&](Tree trial) { return lik.optimise_branch_lengths(trial, 2, 1e-6); };
  for (cycle = 0; cycle < o.cycles; cycle++) {
    const double before = L;
    if (o.model_parameters) {
      for (const auto& group : detail::free_rates(model.kind)) {
        const auto f = [&](double x) {
          Model m = model;
          for (int g : group) m.exchange[static_cast<std::size_t>(g)] = std::exp(x);
          m.prepare();
          lik.set_model(m);
          return refit(out.tree);
        };
        const double x = fit(f, std::log(0.01), std::log(200.0), std::log(model.exchange[static_cast<std::size_t>(group[0])]), 1.0);
        for (int g : group) model.exchange[static_cast<std::size_t>(g)] = std::exp(x);
        model.prepare();
        lik.set_model(model);
      }
      if (rates.categories > 1) {
        const auto f = [&](double x) {
          RateSpec r = rates;
          r.alpha = std::exp(x);
          lik.set_rates(r.rates());
          return refit(out.tree);
        };
        rates.alpha = std::exp(fit(f, std::log(0.02), std::log(100.0), std::log(rates.alpha), 1.0));
        lik.set_rates(rates.rates());
      }
      if (rates.estimate_invariant) {
        const auto f = [&](double x) {
          RateSpec r = rates;
          r.invariant = x;
          lik.set_rates(r.rates());
          return refit(out.tree);
        };
        rates.invariant = fit(f, 0.0, 0.95, rates.invariant, 0.1);
        lik.set_rates(rates.rates());
      }
      L = lik.optimise_branch_lengths(out.tree, 30, o.tolerance);
    }
    if (o.topology) L = lik.nni_search(out.tree, o.rounds, o.tolerance, &out.interchanges);
    if (L - before < o.cycle_tolerance) break;
  }
  out.log_likelihood = L;
  out.model = model;
  out.rates = rates;
  return out;
}

// From a neighbour-joining tree of the JC69 distances (Γ ones when the rates are Γ).
inline MLTree maximum_likelihood(const Patterns& data, const Model& model, const RateSpec& rates, const SearchOptions& o = {}) {
  DistanceOptions d;
  d.model = DistanceModel::jc69;
  if (rates.categories > 1) d.alpha = rates.alpha;
  return maximum_likelihood(data, model, rates, neighbour_joining(distances(data, d)), o);
}

}  // namespace evolution::phylogeny
