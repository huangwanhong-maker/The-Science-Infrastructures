// Checks of the coalescent with recombination (Hudson's algorithm): the marginal trees
// partition the chromosome, and each is a Kingman tree (E[T_MRCA] = 2(1 − 1/n) and E[L] = 2a_n
// at a fixed position whatever ρ); the number of distinct marginal trees grows with ρ; for two
// sequences the covariance of the coalescence times at two loci separated by ρ is
// (ρ + 18)/(ρ² + 13ρ + 18) (Griffiths 1981, Theor. Popul. Biol. 19, 169–186; Hudson 1990,
// Oxford Surveys 7; in units of 2N generations, where Var[T] = 1); the sites of
// recombining sequences show four gametes; r² between sites falls with their distance, and the
// drawn trees are consistent with the sweep.

#include <evolution/population/coalescent.hpp>
#include <evolution/population/diversity.hpp>
#include <evolution/population/linkage.hpp>

#include <calx/verification/check.hpp>

#include <chrono>
#include <cmath>
#include <string>
#include <vector>

using namespace evolution::population;

namespace {

struct Moments {
  double n = 0, sum = 0, sum2 = 0;
  void add(double x) {
    n++;
    sum += x;
    sum2 += x * x;
  }
  double mean() const { return sum / n; }
  double se() const { return std::sqrt((sum2 - sum * sum / n) / (n - 1) / n); }
  double z(double expected) const { return (mean() - expected) / se(); }
};

std::string fmt(double x) { return calx::verification::exponential(x); }

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- the trees partition the chromosome, and each is Kingman's ------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 10;
    spec.rho = 20;
    Moments tmrca, length;
    bool partition = true, rooted = true;
    for (std::uint32_t rep = 0; rep < 20000; rep++) {
      const Genealogy g = coalescent(spec, 11, rep);
      double covered = 0, last = 0;
      for_each_tree(g, [&](const TreeView& v) {
        partition = partition && v.left == last && v.right > v.left;
        covered += v.right - v.left;
        last = v.right;
        rooted = rooted && v.root >= static_cast<std::int32_t>(g.samples) && v.parent[static_cast<std::size_t>(v.root)] < 0;
        if (v.left <= 0.3 && 0.3 < v.right) {
          tmrca.add(v.tmrca());
          length.add(v.length);
        }
      });
      partition = partition && last == 1 && std::fabs(covered - 1) < 1e-12;
    }
    s.require("the marginal trees partition [0, 1) into intervals, each with one root above the samples", partition && rooted);
    s.check("at a fixed position E[T_MRCA] = 2(1 − 1/n) whatever ρ (n = 10, ρ = 20; z-score)", tmrca.z(expected_tmrca(10)), 4, "SE");
    s.check("and E[L] = 2a_n (z-score)", length.z(expected_tree_length(10)), 4, "SE");
  }

  // ---- more recombination, more trees -----------------------------------------------------------------
  {
    std::vector<double> mean_trees;
    std::string detail;
    for (const double rho : {0.0, 1.0, 10.0, 100.0}) {
      CoalescentSpec spec;
      spec.samples = 20;
      spec.rho = rho;
      double total = 0;
      for (std::uint32_t rep = 0; rep < 200; rep++) total += static_cast<double>(tree_count(coalescent(spec, 12, rep)));
      mean_trees.push_back(total / 200);
      detail += "ρ = " + fmt(rho) + ": " + fmt(total / 200) + "; ";
    }
    s.require("the number of distinct marginal trees grows with ρ (one tree at ρ = 0)",
              mean_trees[0] == 1 && mean_trees[1] > mean_trees[0] && mean_trees[2] > mean_trees[1] && mean_trees[3] > mean_trees[2], detail);
  }

  // ---- two loci, two sequences: the covariance of their coalescence times ---------------------------------
  for (const double between : {1.0, 5.0}) {
    CoalescentSpec spec;
    spec.samples = 2;
    spec.rho = 2 * between;  // the loci at 0.25 and 0.75 are half the chromosome apart
    const std::size_t R = 40000;
    std::vector<double> ta(R), tb(R);
    for (std::uint32_t rep = 0; rep < R; rep++) {
      const Genealogy g = coalescent(spec, 13 + static_cast<std::uint64_t>(between), rep);
      ta[rep] = tmrca_at(g, 0.25);
      tb[rep] = tmrca_at(g, 0.75);
    }
    // The covariance against its closed form, with the standard error of the mean of the
    // products (T_A − 1)(T_B − 1), the means being 1 exactly.
    Moments product;
    for (std::size_t r = 0; r < R; r++) product.add((ta[r] - 1) * (tb[r] - 1));
    const double expected = (between + 18) / (between * between + 13 * between + 18);
    s.check("Cov[T_A, T_B] = (ρ + 18)/(ρ² + 13ρ + 18) for n = 2 at ρ = " + fmt(between) + " between the loci (z-score)", product.z(expected), 4, "SE");
  }

  // ---- linkage decays with distance -----------------------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 50;
    spec.theta = 100;
    spec.rho = 200;
    std::vector<double> sum(10, 0), pairs(10, 0);
    bool four = false;
    for (std::uint32_t rep = 0; rep < 20; rep++) {
      Haplotypes h = haplotypes(coalescent(spec, 14, rep));
      // Common sites only (minor allele frequency ≥ 0.1), as r² of rare alleles is noise.
      const AlleleCounts c = allele_counts(h);
      Haplotypes common;
      common.sequences = h.sequences;
      for (std::size_t j = 0; j < h.sites; j++) {
        if (c.alt[j] < 5 || c.alt[j] > 45) continue;
        const auto row = h.site(j);
        common.alleles.insert(common.alleles.end(), row.begin(), row.end());
        common.positions.push_back(h.positions[j]);
        common.sites++;
      }
      const auto bins = ld_decay(common, 0.5, 10);
      for (std::size_t b = 0; b < 10; b++) {
        sum[b] += bins[b].mean_r2 * static_cast<double>(bins[b].pairs);
        pairs[b] += static_cast<double>(bins[b].pairs);
      }
      for (std::size_t x = 0; x + 1 < common.sites && !four; x++)
        for (std::size_t y = x + 1; y < common.sites && !four; y++) {
          int seen = 0;
          for (std::size_t k = 0; k < common.sequences; k++) seen |= 1 << (2 * common(k, x) + common(k, y));
          four = seen == 15;
        }
    }
    const double near = sum[0] / pairs[0], mid = sum[2] / pairs[2], far = sum[9] / pairs[9];
    s.require("recombining sequences show four gametes at some pairs of sites", four);
    s.require("mean r² falls with distance: [0, 0.05) > [0.1, 0.15) > [0.45, 0.5) of the chromosome (ρ = 200, n = 50)",
              near > mid && mid > far, fmt(near) + " > " + fmt(mid) + " > " + fmt(far));
    // Ohta and Kimura (1971): E[r²] → σ²_d = (10 + C)/(22 + 13C + C²) with C = ρd, plus about
    // 1/n of sampling noise; at C = 90 this is 0.0099 + 0.02.
    s.require("between distant sites r² is at the level of sampling noise (below 0.05)", far < 0.05, fmt(far));
  }

  // ---- trees for drawing -------------------------------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 8;
    spec.theta = 6;
    spec.rho = 4;
    const Genealogy g = coalescent(spec, 15);
    const std::vector<Tree> trees = marginal_trees(g);
    bool consistent = trees.size() == tree_count(g);
    std::size_t mutations = 0;
    for (const Tree& t : trees) {
      consistent = consistent && t.node.size() == 2 * spec.samples - 1 && t.children.size() == t.node.size() - 1;
      consistent = consistent && t.root >= 0 && t.time[static_cast<std::size_t>(t.root)] == t.tmrca;
      for (std::size_t v = 0; v < t.node.size(); v++) {
        const auto cb = t.child_begin[v], ce = t.child_begin[v + 1];
        consistent = consistent && (v < spec.samples ? ce == cb : ce - cb == 2);
        for (auto c = cb; c < ce; c++)
          consistent = consistent && t.parent[static_cast<std::size_t>(t.children[c])] == static_cast<std::int32_t>(v) &&
                       t.time[static_cast<std::size_t>(t.children[c])] < t.time[v];
      }
      double length = 0;
      for (std::size_t v = 0; v < t.node.size(); v++)
        if (t.parent[v] >= 0) length += t.time[static_cast<std::size_t>(t.parent[v])] - t.time[v];
      consistent = consistent && std::fabs(length - t.length) < 1e-12;
      for (const TreeMutation& m : t.mutations)
        consistent = consistent && m.node >= 0 && t.parent[static_cast<std::size_t>(m.node)] >= 0 && m.time >= t.time[static_cast<std::size_t>(m.node)] &&
                     m.time <= t.time[static_cast<std::size_t>(t.parent[static_cast<std::size_t>(m.node)])];
      mutations += t.mutations.size();
    }
    const Tree mid = tree_at(g, 0.5);
    s.require("the drawn trees are binary, children younger than their parents, lengths and mutations on their branches", consistent && mutations == g.mutations.size(),
              std::to_string(trees.size()) + " trees");
    s.require("the tree at a position is the sweep's", mid.left <= 0.5 && 0.5 < mid.right && mid.tmrca == tmrca_at(g, 0.5));
  }

  // ---- speed --------------------------------------------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 100;
    spec.theta = 100;
    spec.rho = 100;
    const auto t0 = std::chrono::steady_clock::now();
    const Genealogy g = coalescent(spec, 16);
    const Haplotypes h = haplotypes(g);
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    s.require("100 samples with ρ = 100 and θ = 100 (the ARG, its trees and haplotypes) in under 500 ms", ms < 500,
              fmt(ms) + " ms, " + std::to_string(g.recombinations) + " recombinations, " + std::to_string(tree_count(g)) + " trees, " +
                  std::to_string(h.sites) + " sites");
  }
  return s.done();
}
