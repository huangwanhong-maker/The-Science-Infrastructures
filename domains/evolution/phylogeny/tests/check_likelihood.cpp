// Checks of the likelihood: two sequences under JC69 against the closed form, and the maximum-
// likelihood distance equal to the JC69 distance (a theorem); pruning against the sum over every
// assignment of internal states, under GTR + Γ; the same likelihood wherever the root is put
// (the pulley principle); rescaling on a tree of 700 leaves, whose sites would underflow, the
// whole equal to the sum of its sites; branch lengths fitted to at least the likelihood of the
// true ones; maximum likelihood (HKY85 + Γ) and neighbour joining recovering the topology of a
// simulated 12-taxon tree, with κ and α near the truth; and the NNI search from a random start
// climbing to a tree no interchange improves.

#include <evolution/phylogeny/distance.hpp>
#include <evolution/phylogeny/distance_trees.hpp>
#include <evolution/phylogeny/likelihood.hpp>
#include <evolution/phylogeny/newick.hpp>
#include <evolution/phylogeny/simulate.hpp>
#include <evolution/phylogeny/splits.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace evolution::phylogeny;

namespace {

Tree pair_tree(double a, double b) {
  Tree t;
  const int r = t.add(-1);
  t.add(r, a, "x");
  t.add(r, b, "y");
  return t;
}

// The likelihood by brute force: every assignment of states to the internal nodes.
double brute_log_likelihood(const Tree& t, const Patterns& p, const Model& m, const Rates& rates) {
  std::vector<int> internal;
  for (std::size_t v = 0; v < t.size(); v++)
    if (!t.nodes[v].children.empty()) internal.push_back(static_cast<int>(v));
  std::vector<double> P(rates.size() * t.size() * 16);
  for (std::size_t c = 0; c < rates.size(); c++)
    for (std::size_t v = 0; v < t.size(); v++)
      if (static_cast<int>(v) != t.root) m.transition(rates.rate[c] * t.nodes[v].length, &P[(c * t.size() + v) * 16]);
  double L = 0;
  std::vector<int> state(t.size(), 0);
  for (std::size_t q = 0; q < p.count; q++) {
    double site = 0;
    for (std::size_t c = 0; c < rates.size(); c++) {
      double sum = 0;
      const std::size_t total = std::size_t{1} << (2 * internal.size());
      for (std::size_t code = 0; code < total; code++) {
        for (std::size_t k = 0; k < internal.size(); k++) state[static_cast<std::size_t>(internal[k])] = static_cast<int>((code >> (2 * k)) & 3);
        double prod = m.pi[static_cast<std::size_t>(state[static_cast<std::size_t>(t.root)])];
        for (std::size_t v = 0; v < t.size() && prod > 0; v++) {
          if (static_cast<int>(v) == t.root) continue;
          const int x = state[static_cast<std::size_t>(t.nodes[v].parent)];
          if (t.nodes[v].children.empty()) {
            const Code y = p.states[q * p.taxa + p.index(t.nodes[v].label)];
            if (y >= 4) continue;  // missing: sums to 1
            prod *= P[(c * t.size() + v) * 16 + static_cast<std::size_t>(x * 4 + y)];
          } else {
            prod *= P[(c * t.size() + v) * 16 + static_cast<std::size_t>(x * 4 + state[v])];
          }
        }
        sum += prod;
      }
      site += rates.weight[c] * sum;
    }
    L += p.weight[q] * std::log(site);
  }
  return L;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- two sequences ---------------------------------------------------------------------------------
  {
    const Alignment a = simulate_alignment(pair_tree(0.08, 0.12), jc69(), Rates{}, 500, 3);
    const Patterns p = compress(a);
    double same = 0, diff = 0;
    for (std::size_t i = 0; i < a.sites(); i++) (a.rows[0][i] == a.rows[1][i] ? same : diff) += 1;
    Likelihood lik(p, jc69(), Rates{});
    const double t = 0.2, e = std::exp(-4 * t / 3);
    const double closed = same * std::log(0.25 * (0.25 + 0.75 * e)) + diff * std::log(0.25 * (0.25 - 0.25 * e));
    s.check("two sequences under JC69: the log-likelihood is its closed form", lik.log_likelihood(pair_tree(0.05, 0.15)) - closed, 1e-9);
    Tree fit = pair_tree(0.3, 0.3);
    lik.optimise_branch_lengths(fit);
    const double ml = fit.total_length(), jc = jc69_distance(diff / (same + diff)).d;
    s.check("the maximum-likelihood distance of two sequences under JC69 is the JC69 distance", ml - jc, 1e-7);
  }

  // ---- pruning against brute force ----------------------------------------------------------------------
  const Model g = gtr({1.3, 4.1, 0.7, 1.9, 5.2, 1.0}, {0.1, 0.4, 0.3, 0.2});
  const Rates gi = gamma_rates(0.7, 3, 0.15);
  {
    const Tree t = read_newick("((a:0.1,b:0.3):0.05,(c:0.2,(d:0.15,e:0.4):0.1):0.2);").trees[0];
    Alignment al = simulate_alignment(t, g, gi, 60, 4);
    al.rows[2][5] = genomics::sequence::GAP;  // missing data too
    al.rows[4][7] = genomics::sequence::N;
    const Patterns p = compress(al);
    Likelihood lik(p, g, gi);
    s.check("pruning equals the sum over all 4^4 assignments of internal states, under GTR + Γ(3) + I",
            lik.log_likelihood(t) - brute_log_likelihood(t, p, g, gi), 1e-10);
  }

  // ---- the pulley principle ---------------------------------------------------------------------------------
  const Tree truth = random_unrooted_tree(12, 21, 0.08, 0.02);
  {
    const Alignment al = simulate_alignment(truth, g, gi, 1000, 5);
    const Patterns p = compress(al);
    Likelihood lik(p, g, gi);
    const double L0 = lik.log_likelihood(truth);
    double worst = 0;
    for (std::size_t v = 0; v < truth.size(); v++) {
      if (static_cast<int>(v) == truth.root) continue;
      worst = std::max(worst, std::fabs(lik.log_likelihood(reroot(truth, static_cast<int>(v), truth.nodes[v].length * 0.3)) - L0));
      if (!truth.nodes[v].children.empty()) worst = std::max(worst, std::fabs(lik.log_likelihood(root_at(truth, static_cast<int>(v))) - L0));
    }
    s.check("the likelihood is the same with the root on every branch and at every node (pulley principle)", worst, 1e-8);
  }

  // ---- rescaling --------------------------------------------------------------------------------------------
  {
    const Tree big = random_unrooted_tree(700, 8, 0.2, 0.05);
    const Model h = hky85(3, {0.3, 0.2, 0.2, 0.3});
    const Alignment al = simulate_alignment(big, h, gamma_rates(0.5, 4), 12, 9);
    const Patterns p = compress(al);
    Likelihood lik(p, h, gamma_rates(0.5, 4));
    const double L = lik.log_likelihood(big);
    double sum = 0, smallest = 0;
    for (std::size_t site = 0; site < al.sites(); site++) {
      Alignment one;
      one.names = al.names;
      for (const auto& r : al.rows) one.rows.push_back(Codes{r[site]});
      Likelihood l1(compress(one), h, gamma_rates(0.5, 4));
      const double x = l1.log_likelihood(big);
      sum += x;
      smallest = std::min(smallest, x);
    }
    s.require("on 700 leaves a site's likelihood is below the smallest double (its log below −708), and is still computed",
              std::isfinite(L) && smallest < -708, std::to_string(smallest));
    s.check("the log-likelihood of the sites together is the sum of theirs", L - sum, 1e-8 * std::fabs(L));
    s.check("and the same with the root elsewhere", lik.log_likelihood(root_at(big, big[big.leaves()[100]].parent)) - L, 1e-8 * std::fabs(L));
  }

  // ---- branch lengths and maximum likelihood ------------------------------------------------------------------
  {
    const Model h = hky85(4, {0.3, 0.2, 0.25, 0.25});
    const Rates r = gamma_rates(0.8, 4);
    const Alignment al = simulate_alignment(truth, h, r, 3000, 6);
    const Patterns p = compress(al);
    Likelihood lik(p, h, r);
    const double at_truth = lik.log_likelihood(truth);
    Tree fit = truth;
    for (auto& n : fit.nodes) n.length = 0.1;
    const double L = lik.optimise_branch_lengths(fit);
    double err = 0, total = 0;
    for (std::size_t v = 0; v < fit.size(); v++)
      if (static_cast<int>(v) != fit.root) {
        err += std::fabs(fit.nodes[v].length - truth.nodes[v].length);
        total += truth.nodes[v].length;
      }
    s.require("fitted branch lengths have at least the likelihood of the true ones", L >= at_truth - 1e-9);
    s.check("and are near them: the sum of errors over the tree length", err / total, 0.12);

    DistanceOptions d;
    d.model = DistanceModel::tn93;
    d.alpha = 0.8;
    const Tree nj = neighbour_joining(distances(p, d));
    s.require("neighbour joining on TN93 + Γ distances recovers the true topology (RF 0)", robinson_foulds(nj, truth) == 0);
    RateSpec spec;
    spec.categories = 4;
    spec.alpha = 0.5;
    const MLTree ml = maximum_likelihood(p, hky85(2, {0.25, 0.25, 0.25, 0.25}), spec, nj);
    s.require("maximum likelihood (HKY85 + Γ4) recovers the true topology (RF 0)", robinson_foulds(ml.tree, truth) == 0);
    s.check("κ is estimated near its true 4 (relative error)", ml.model.exchange[AG] / 4 - 1, 0.15);
    s.check("α near its true 0.8 (relative error)", ml.rates.alpha / 0.8 - 1, 0.25);
    Likelihood true_model(p, h, r);
    Tree tl = truth;
    const double best_true = true_model.optimise_branch_lengths(tl);
    s.require("the maximum likelihood is at least that of the true tree and model with fitted lengths", ml.log_likelihood >= best_true - 1e-6,
              std::to_string(ml.log_likelihood) + " vs " + std::to_string(best_true));

    // From a random topology, interchanges climb to a tree no interchange improves.
    Tree start = random_unrooted_tree(12, 1234, 0.1);
    Likelihood l2(p, h, r);
    Tree s0 = start;
    const double L_start = l2.optimise_branch_lengths(s0);
    int moves = 0;
    const double L_end = l2.nni_search(start, 200, 1e-6, &moves);
    int more = 0;
    Tree again = start;
    l2.nni_search(again, 200, 1e-6, &more);
    s.require("the NNI search from a random start gains likelihood and makes interchanges", L_end > L_start && moves > 0,
              std::to_string(moves) + " moves, " + std::to_string(L_start) + " → " + std::to_string(L_end));
    s.require("and stops at a local maximum: a second search makes none", more == 0);
    s.require("here the true tree (RF to it, of 18 at most)", robinson_foulds(start, truth) == 0, std::to_string(robinson_foulds(start, truth)));
  }

  return s.done();
}
