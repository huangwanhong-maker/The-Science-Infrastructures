// Checks of the bootstrap: a replicate's weights are L draws of the L sites, a function of the
// seed and the replicate alone, their mean the original counts; neighbour-joining bootstrap
// supports near 1 for the splits of a simulated tree with long internal branches; a split with no
// length behind it poorly supported; and every support and every replicate tree the same on 1
// and 8 threads.

#include <evolution/phylogeny/bootstrap.hpp>
#include <evolution/phylogeny/newick.hpp>
#include <evolution/phylogeny/simulate.hpp>
#include <evolution/phylogeny/splits.hpp>

#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace evolution::phylogeny;

int main() {
  calx::verification::Suite s;

  const Tree truth = random_unrooted_tree(12, 31, 0.06, 0.03);
  const Patterns p = compress(simulate_alignment(truth, k80(3), gamma_rates(1, 4), 1500, 32));
  {
    const auto w = bootstrap_weights(p, 9, 4);
    double total = 0;
    for (double x : w) total += x;
    s.check("a replicate's weights add up to the number of sites", total - p.sites(), 0);
    s.require("and are a function of the seed and the replicate", w == bootstrap_weights(p, 9, 4) && w != bootstrap_weights(p, 9, 5) && w != bootstrap_weights(p, 10, 4));
    std::vector<double> mean(p.count, 0);
    const int B = 400;
    for (int b = 0; b < B; b++) {
      const auto x = bootstrap_weights(p, 9, static_cast<std::uint32_t>(b));
      for (std::size_t q = 0; q < p.count; q++) mean[q] += x[q] / B;
    }
    // Each weight is binomial(L, n_q / L): its mean over B replicates has standard error √(n_q (1 − n_q/L) / B).
    double worst = 0;
    for (std::size_t q = 0; q < p.count; q++) {
      const double sd = std::sqrt(p.weight[q] * (1 - p.weight[q] / p.sites()) / B);
      worst = std::max(worst, std::fabs(mean[q] - p.weight[q]) / sd);
    }
    s.check("their mean is the original counts, every pattern within 5 standard errors", worst, 5);
  }

  DistanceOptions d;
  d.model = DistanceModel::k80;
  calx::scheduler::ThreadPool one(1), eight(8);
  const Bootstrap b1 = bootstrap_neighbour_joining(p, d, 200, 77, one);
  const Bootstrap b8 = bootstrap_neighbour_joining(p, d, 200, 77, eight);
  {
    s.require("the reference is the true tree (RF 0)", robinson_foulds(b1.tree, truth) == 0);
    double lowest = 1;
    std::size_t splits_seen = 0;
    for (double x : b1.support)
      if (!std::isnan(x)) {
        lowest = std::min(lowest, x);
        splits_seen++;
      }
    s.require("every one of its 9 splits has a support", splits_seen == 9);
    s.check("supports of splits behind internal branches of 0.03 or more over 1500 sites are near 1 (1 − the least)", 1 - lowest, 0.05);
    bool same = b1.support.size() == b8.support.size();
    for (std::size_t v = 0; same && v < b1.support.size(); v++)
      same = (std::isnan(b1.support[v]) && std::isnan(b8.support[v])) || b1.support[v] == b8.support[v];
    for (std::size_t r = 0; same && r < b1.replicates.size(); r++) same = write_newick(b1.replicates[r]) == write_newick(b8.replicates[r]);
    s.require("the supports and all 200 replicate trees are the same on 1 and 8 threads", same);
    s.require("the trace records every replicate and the merge", b8.trace.records.size() == 201 && b8.trace.threads == 8);
  }
  {
    // A split with no branch behind it: the internal branch above one cherry set to zero.
    Tree flat = truth;
    int target = -1;
    for (std::size_t v = 0; v < flat.size(); v++)
      if (static_cast<int>(v) != flat.root && !flat.nodes[v].children.empty() && flat.nodes[v].parent != flat.root) {
        target = static_cast<int>(v);
        break;
      }
    flat.nodes[static_cast<std::size_t>(target)].length = 0;
    const Patterns q = compress(simulate_alignment(flat, k80(3), gamma_rates(1, 4), 1500, 33));
    const Bootstrap b = bootstrap_neighbour_joining(q, d, 200, 78, eight);
    // The splits of the reference that are not splits of the flattened tree's resolution are the
    // unsupported one; find the least support.
    double least = 1;
    for (double x : b.support)
      if (!std::isnan(x)) least = std::min(least, x);
    s.require("a split with no branch behind it has a support well below 1 (under 0.9)", least < 0.9, std::to_string(least));
  }
  return s.done();
}
