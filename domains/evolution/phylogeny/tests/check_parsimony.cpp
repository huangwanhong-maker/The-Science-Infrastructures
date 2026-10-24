// Checks of parsimony: Fitch's score (Hartigan's at nodes of higher degree) against the least
// number of changes over every assignment of states to the internal nodes, on binary trees and
// on a tree with polytomies, with missing data; Sankoff's algorithm with unit costs equal to
// Fitch's, with other costs against the same enumeration, and its traceback attaining its
// score; the exact search by branch and bound against the minimum over all 10,395 unrooted trees
// of 8 taxa; and the heuristic search recovering a simulated tree.

#include <evolution/phylogeny/newick.hpp>
#include <evolution/phylogeny/parsimony.hpp>
#include <evolution/phylogeny/simulate.hpp>
#include <evolution/phylogeny/splits.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <vector>

using namespace evolution::phylogeny;

namespace {

// The least cost of a tree over every assignment of states to its internal nodes (and to leaves
// whose state is missing), site by site.
double brute_force(const Tree& t, const Patterns& p, const CostMatrix& cost) {
  std::vector<int> free;  // internal nodes and missing leaves, per pattern
  double total = 0;
  std::vector<int> state(t.size(), 0);
  for (std::size_t q = 0; q < p.count; q++) {
    free.clear();
    for (std::size_t v = 0; v < t.size(); v++) {
      if (!t.nodes[v].children.empty()) {
        free.push_back(static_cast<int>(v));
        continue;
      }
      const Code c = p.state(q, p.index(t.nodes[v].label));
      if (c < 4) state[v] = c;
      else free.push_back(static_cast<int>(v));
    }
    double best = std::numeric_limits<double>::infinity();
    const std::size_t total_codes = std::size_t{1} << (2 * free.size());
    for (std::size_t code = 0; code < total_codes; code++) {
      for (std::size_t k = 0; k < free.size(); k++) state[static_cast<std::size_t>(free[k])] = static_cast<int>((code >> (2 * k)) & 3);
      double c = 0;
      for (std::size_t v = 0; v < t.size(); v++)
        if (static_cast<int>(v) != t.root) c += cost[static_cast<std::size_t>(state[static_cast<std::size_t>(t.nodes[v].parent)] * 4 + state[v])];
      best = std::min(best, c);
    }
    total += p.weight[q] * best;
  }
  return total;
}

// Every unrooted binary tree on the taxa, by inserting each taxon on every branch in turn.
void every_tree(const std::vector<std::string>& names, const std::function<void(const Tree&)>& f) {
  Tree t;
  const int r = t.add(-1);
  for (int k = 0; k < 3; k++) t.add(r, ABSENT, names[static_cast<std::size_t>(k)]);
  std::function<void(const Tree&, std::size_t)> grow = [&](const Tree& s, std::size_t k) {
    if (k == names.size()) {
      f(s);
      return;
    }
    for (std::size_t u = 0; u < s.size(); u++) {
      if (static_cast<int>(u) == s.root) continue;
      Tree x = s;
      const int p = x[static_cast<int>(u)].parent;
      const int w = static_cast<int>(x.nodes.size());
      Node mid;
      mid.parent = p;
      mid.children = {static_cast<int>(u)};
      x.nodes.push_back(mid);
      for (int& c : x[p].children)
        if (c == static_cast<int>(u)) c = w;
      x[static_cast<int>(u)].parent = w;
      x.add(w, ABSENT, names[k]);
      grow(x, k + 1);
    }
  };
  grow(t, 3);
}

}  // namespace

int main() {
  calx::verification::Suite s;

  const Model k = k80(2);
  // ---- Fitch and Sankoff against enumeration -------------------------------------------------------------
  {
    double worst_fitch = 0, worst_sankoff = 0, worst_unit = 0, worst_trace = 0;
    for (std::uint64_t seed = 1; seed <= 5; seed++) {
      const Tree t = random_unrooted_tree(6, seed, 0.3);
      Alignment al = simulate_alignment(t, k, Rates{}, 40, seed + 100);
      al.rows[1][3] = genomics::sequence::GAP;
      al.rows[4][9] = genomics::sequence::N;
      const Patterns p = compress(al);
      worst_fitch = std::max(worst_fitch, std::fabs(fitch_score(t, p) - brute_force(t, p, unit_costs())));
      const CostMatrix tc = transition_costs(1, 2.5);
      const auto sk = sankoff(t, p, tc);
      worst_sankoff = std::max(worst_sankoff, std::fabs(sk.score - brute_force(t, p, tc)));
      worst_unit = std::max(worst_unit, std::fabs(sankoff(t, p, unit_costs()).score - fitch_score(t, p)));
      // The traceback's assignment costs the score (missing leaves take the state it gives them).
      double c = 0;
      for (std::size_t q = 0; q < p.count; q++)
        for (std::size_t v = 0; v < t.size(); v++)
          if (static_cast<int>(v) != t.root)
            c += p.weight[q] * tc[static_cast<std::size_t>(sk.state(t.nodes[v].parent, q) * 4 + sk.state(static_cast<int>(v), q))];
      worst_trace = std::max(worst_trace, std::fabs(c - sk.score));
    }
    s.check("Fitch's score is the least number of changes over all 4^4 internal assignments, on five 6-taxon trees with missing data", worst_fitch, 0);
    s.check("Sankoff's score with transitions 1 and transversions 2.5 is the least cost by enumeration", worst_sankoff, 1e-12);
    s.check("Sankoff with unit costs is Fitch", worst_unit, 1e-12);
    s.check("Sankoff's traceback gives an assignment that costs its score", worst_trace, 1e-12);
  }
  {
    // A tree with a root of degree four and a node of degree three: Hartigan's rule.
    const Tree t = read_newick("((a,b,c),d,(e,f),g);").trees[0];
    const Alignment al = simulate_alignment(random_unrooted_tree(7, 3, 0.4), k, Rates{}, 60, 7);
    Alignment named = al;
    const std::vector<std::string> letters{"a", "b", "c", "d", "e", "f", "g"};
    named.names = letters;
    const Patterns p = compress(named);
    s.check("Hartigan's rule at polytomies is the least number of changes by enumeration", fitch_score(t, p) - brute_force(t, p, unit_costs()), 0);
  }

  // ---- searches ------------------------------------------------------------------------------------------
  {
    const Tree truth = random_unrooted_tree(8, 41, 0.08, 0.01);
    const Alignment al = simulate_alignment(truth, k, gamma_rates(1, 4), 120, 42);
    const Patterns p = compress(al);
    double best = std::numeric_limits<double>::infinity();
    std::size_t count = 0;
    every_tree(p.names, [&](const Tree& t) {
      count++;
      best = std::min(best, fitch_score(t, p));
    });
    const auto exact = parsimony_exact(p);
    const auto heur = parsimony_search(p);
    s.require("there are 10,395 unrooted binary trees on 8 taxa, (2·8 − 5)!!", count == 10395);
    s.check("branch and bound finds the least score over all of them", exact.score - best, 0);
    s.require("having scored fewer trees than all of them", exact.trees_scored < 10395, std::to_string(exact.trees_scored) + " trees scored");
    s.check("and its tree has that score", fitch_score(exact.tree, p) - exact.score, 0);
    s.require("stepwise addition and interchanges reach a score no lower than the least", heur.score >= best);
  }
  {
    const Tree truth = random_unrooted_tree(14, 5, 0.05, 0.02);
    const Patterns p = compress(simulate_alignment(truth, k, Rates{}, 3000, 6));
    const auto heur = parsimony_search(p);
    s.require("the heuristic search recovers a simulated 14-taxon tree (short branches, 3000 sites; RF 0)", robinson_foulds(heur.tree, truth) == 0,
              std::to_string(robinson_foulds(heur.tree, truth)));
    s.require("whose score is no more than the true tree's", heur.score <= fitch_score(truth, p));
  }
  return s.done();
}
