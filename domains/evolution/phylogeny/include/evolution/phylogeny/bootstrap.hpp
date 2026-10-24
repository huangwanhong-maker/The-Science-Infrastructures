// The nonparametric bootstrap of a tree (Felsenstein 1985, "Confidence limits on phylogenies:
// an approach using the bootstrap", Evolution 39, 783–791): the alignment's L sites resampled
// with replacement, L of them, a tree built from each replicate by the same method, and the
// support of each split of the reference tree the fraction of replicate trees that hold it.
//
// A replicate is the same site patterns with new weights, the number of times each pattern's
// sites were drawn; the draws of replicate b come from the stream (seed, BOOTSTRAP, b), so a
// replicate is a function of the seed and its number alone. The replicates are the tasks of a
// calx::scheduler::TaskGraph, each writing its own tree, and the supports a final task that
// reads them in order, so the result is the same on any number of threads; the trace says where
// and when each ran. Cost: the method's, B times, divided among the threads.
//
// Supports are fractions in [0, 1], written on the reference's internal branches (Node::support).
// A bootstrap support is not the probability that a split is true (Hillis and Bull 1993, Syst.
// Biol. 42, 182–192, find it conservative above 70%, under their conditions).

#pragma once

#include <evolution/phylogeny/alignment.hpp>
#include <evolution/phylogeny/distance.hpp>
#include <evolution/phylogeny/distance_trees.hpp>
#include <evolution/phylogeny/splits.hpp>
#include <evolution/phylogeny/tree.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/task_graph.hpp>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace evolution::phylogeny {

inline constexpr std::uint32_t BOOTSTRAP = 0x424f4f54;

// The weights of replicate b: each of the L sites drawn uniformly, L times.
inline std::vector<double> bootstrap_weights(const Patterns& p, std::uint64_t seed, std::uint32_t replicate) {
  const std::size_t L = p.of_site.size();
  if (L == 0) throw std::invalid_argument("evolution::phylogeny::bootstrap_weights: patterns from an alignment (with of_site)");
  calx::random::Stream r(seed, BOOTSTRAP, replicate, 0);
  std::vector<double> w(p.count, 0);
  for (std::size_t k = 0; k < L; k++) w[p.of_site[static_cast<std::size_t>(r.below(L))]] += 1;
  return w;
}

struct Bootstrap {
  Tree tree;                      // the reference, its branches' supports the bootstrap's
  std::vector<double> support;    // by node of the reference (NaN where no split is made)
  std::vector<Tree> replicates;
  calx::scheduler::Trace trace;
};

// Supports for `reference` from B replicates, each tree built by build(patterns, thread), which
// must be safe to call from several threads at once.
template <class Build>
Bootstrap bootstrap(const Patterns& data, const Tree& reference, std::size_t replicates, std::uint64_t seed, calx::scheduler::ThreadPool& pool,
                    Build&& build) {
  Bootstrap out;
  out.replicates.resize(replicates);
  calx::scheduler::TaskGraph g;
  std::vector<std::size_t> ids;
  ids.reserve(replicates);
  for (std::size_t b = 0; b < replicates; b++)
    ids.push_back(g.add("replicate " + std::to_string(b), [&, b](std::size_t thread) {
      Patterns rep = data;
      rep.weight = bootstrap_weights(data, seed, static_cast<std::uint32_t>(b));
      out.replicates[b] = build(static_cast<const Patterns&>(rep), thread);
    }));
  g.add("supports", [&](std::size_t) {
    out.tree = reference;
    out.support = split_frequencies(reference, out.replicates);
    for (std::size_t v = 0; v < out.tree.size(); v++) out.tree.nodes[v].support = out.support[v];
  }, ids, 0.1);
  out.trace = g.run(pool);
  return out;
}

// The bootstrap of neighbour joining on distances of the given kind: the reference is the tree
// of the data themselves.
inline Bootstrap bootstrap_neighbour_joining(const Patterns& data, const DistanceOptions& d, std::size_t replicates, std::uint64_t seed,
                                             calx::scheduler::ThreadPool& pool) {
  const Tree reference = neighbour_joining(distances(data, d));
  return bootstrap(data, reference, replicates, seed, pool, [&](const Patterns& rep, std::size_t) { return neighbour_joining(distances(rep, d)); });
}

}  // namespace evolution::phylogeny
