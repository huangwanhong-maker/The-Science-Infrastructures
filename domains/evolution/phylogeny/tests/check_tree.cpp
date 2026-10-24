// Checks of trees: Newick read and written back unchanged, with quoted labels, comments,
// supports and lengths, and its faults found where they are; rerooting keeps every leaf-to-leaf
// distance; the midpoint root halves the longest path; a binary tree on n leaves has n − 3
// non-trivial splits, which give it back; Robinson–Foulds distances of trees worked by hand;
// majority-rule consensus with its frequencies; layouts; random unrooted topologies uniform
// over the (2n − 5)!! of them (χ²); the mean height of Kingman's coalescent, 2(1 − 1/n), and of
// the Yule process as simulated, Σ_{k=2}^{n} 1/(kλ).

#include <evolution/phylogeny/layout.hpp>
#include <evolution/phylogeny/newick.hpp>
#include <evolution/phylogeny/simulate.hpp>
#include <evolution/phylogeny/splits.hpp>
#include <evolution/phylogeny/tree.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

using namespace evolution::phylogeny;

namespace {

double max_difference(const std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size()) return INFINITY;
  double m = 0;
  for (std::size_t i = 0; i < a.size(); i++) m = std::max(m, std::fabs(a[i] - b[i]));
  return m;
}

Tree one(const std::string& text) {
  auto r = read_newick(text);
  return r.trees.empty() ? Tree{} : r.trees[0];
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- Newick ---------------------------------------------------------------------------------------
  {
    const std::string text = "((A:1,B:2)90:0.5,(C:1,D:1.25)75:0.3,E:2);";
    const auto r = read_newick(text);
    s.require("a tree with lengths and supports is read without problems", r.trees.size() == 1 && r.problems.empty());
    const Tree& t = r.trees[0];
    s.require("it has five leaves, eight nodes, and is unrooted (its root has three children)", t.leaf_count() == 5 && t.size() == 8 && !t.rooted);
    const int ab = t[t.find_leaf("A")].parent, cd = t[t.find_leaf("C")].parent;
    s.require("the internal labels 90 and 75 are the supports of their branches, with lengths 0.5 and 0.3",
              t[ab].support == 90 && t[cd].support == 75 && t[ab].length == 0.5 && t[cd].length == 0.3 && t[ab].label.empty());
    s.require("it is written back exactly as it was read", write_newick(t) == text, write_newick(t));
  }
  {
    const std::string text = "('Homo sapiens':0.1,'it''s (quoted)':0.2,Pan_troglodytes:0.3,'a_b':0.4);";
    const Tree t = one(text);
    s.require("quoted labels keep blanks and doubled quotes; an unquoted underscore is a blank (Felsenstein)",
              t.find_leaf("Homo sapiens") >= 0 && t.find_leaf("it's (quoted)") >= 0 && t.find_leaf("Pan troglodytes") >= 0 && t.find_leaf("a_b") >= 0);
    const Tree back = one(write_newick(t));
    s.require("labels that need quotes are quoted when written, and read back the same", back.leaf_labels() == t.leaf_labels(), write_newick(t));
  }
  {
    const auto r = read_newick("[&R] [a [nested] comment]((A,B)[&&NHX:S=human],(C , D)) ;\n[&U]((A,B),(C,D));");
    s.require("comments, nested or not, and blanks are skipped; [&R] and [&U] say whether a tree is rooted",
              r.problems.empty() && r.trees.size() == 2 && r.trees[0].rooted && !r.trees[1].rooted && r.trees[0].leaf_count() == 4);
  }
  {
    const auto r = read_newick("(A,B,C);\n((A,B),\n  C:x,D);\n(A,B));\n((A,B),C;\n(D,E,F);");
    std::string where;
    for (const auto& p : r.problems) where += std::to_string(p.line) + " " + p.message + "; ";
    s.require("of five trees with three faults, the two good ones are read", r.trees.size() == 2 && r.trees[1].find_leaf("F") >= 0, where);
    s.require("the faults are found on their lines: a length not a number (line 3, column 5), a ')' without its '(' (line 4), a '(' not closed (line 5)",
              r.problems.size() == 3 && r.problems[0].line == 3 && r.problems[0].message.rfind("column 5:", 0) == 0 && r.problems[1].line == 4 &&
                  r.problems[2].line == 5,
              where);
    const auto bad = read_newick("((A,B)[unclosed comment");
    s.require("an unclosed comment is a problem, not an exception", bad.trees.empty() && bad.problems.size() == 1);
  }

  // ---- rerooting -----------------------------------------------------------------------------------
  const Tree t12 = random_unrooted_tree(12, 7, 0.1, 0.01);
  const auto taxa = taxa_of(t12);
  const auto D = leaf_distances(t12, taxa);
  {
    bool same = true;
    double worst = 0;
    for (std::size_t v = 0; v < t12.size(); v++) {
      if (static_cast<int>(v) == t12.root) continue;
      const Tree r = reroot(t12, static_cast<int>(v), t12[static_cast<int>(v)].length * 0.3);
      worst = std::max(worst, max_difference(leaf_distances(r, taxa), D));
      same = same && r.rooted && r[r.root].children.size() == 2 && r.leaf_count() == 12 && robinson_foulds(r, t12) == 0;
      const Tree u = unroot(r);
      same = same && !u.rooted && u[u.root].children.size() == 3 && robinson_foulds(u, t12) == 0 && u.size() == t12.size();
    }
    s.require("rerooting on every branch gives a rooted tree with the same splits, and unrooting it the original shape", same);
    s.check("rerooting keeps every leaf-to-leaf distance", worst, 1e-12);
    const Tree h = root_at(t12, t12[t12.find_leaf("t5")].parent);
    s.check("a tree hung from another internal node keeps them too", max_difference(leaf_distances(h, taxa), D), 1e-12);
  }
  {
    const Tree m = midpoint_root(t12);
    double longest = 0;
    for (double x : D) longest = std::max(longest, x);
    const auto dep = depths(m);
    double deepest = 0;
    std::vector<int> sides;
    for (int v : m.leaves()) deepest = std::max(deepest, dep[static_cast<std::size_t>(v)]);
    // The two deepest leaves lie on different sides of the root.
    const auto side_of = [&](int v) {
      while (m[v].parent != m.root) v = m[v].parent;
      return v;
    };
    for (int v : m.leaves())
      if (std::fabs(dep[static_cast<std::size_t>(v)] - deepest) < 1e-12) sides.push_back(side_of(v));
    std::sort(sides.begin(), sides.end());
    sides.erase(std::unique(sides.begin(), sides.end()), sides.end());
    s.check("the midpoint root is half the longest leaf-to-leaf path from the deepest leaves", deepest - longest / 2, 1e-12);
    s.require("and the deepest leaves lie on both sides of it", sides.size() == 2);
  }

  // ---- splits, Robinson–Foulds, consensus -------------------------------------------------------------
  {
    const auto sp = splits(t12, taxa);
    s.require("a binary unrooted tree on 12 leaves has 9 non-trivial splits", sp.size() == 9);
    bool compat = true;
    for (std::size_t i = 0; i < sp.size(); i++)
      for (std::size_t j = 0; j < sp.size(); j++) compat = compat && compatible(sp[i].side, sp[j].side);
    s.require("the splits of a tree are pairwise compatible (Buneman)", compat);
    const Tree back = tree_from_splits(taxa, splits(t12, taxa, true));
    s.require("the tree of a tree's splits is the tree (RF 0)", robinson_foulds(back, t12) == 0);
    s.check("with the trivial splits' lengths it has the same leaf distances", max_difference(leaf_distances(back, taxa), D), 1e-12);
  }
  {
    const Tree a = one("((A,B),(C,D),E);"), b = one("((A,C),(B,D),E);"), c = one("((B,A),E,(D,C));");
    s.require("RF of ((A,B),(C,D),E) and ((A,C),(B,D),E) is 4, worked by hand: AB|CDE, CD|ABE against AC|BDE, BD|ACE",
              robinson_foulds(a, b) == 4 && normalised_robinson_foulds(a, b) == 1);
    s.require("RF is 0 between two writings of one unrooted tree", robinson_foulds(a, c) == 0);
    const Tree r = one("(((A,B),(C,D)),E);");
    s.require("and between a rooted tree and its unrooted shape", robinson_foulds(a, r) == 0 && r.rooted && splits(r, taxa_of(r)).size() == 2);
    bool threw = false;
    try {
      (void)robinson_foulds(a, one("((A,B),(C,F),E);"));
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    s.require("comparing trees on different leaves is an invalid argument", threw);
  }
  {
    std::vector<Tree> ts{one("((A:1,B:1):2,C:1,(D:1,E:1):1);"), one("((A:1,B:1):4,C:1,(D:1,E:1):1);"), one("((A:1,C:1):1,B:1,(D:1,E:1):1);")};
    const Tree m = majority_consensus(ts);
    const Tree expect = one("((A,B),C,(D,E));");
    // Over the taxa A–E a split is held by its side without A: AB|CDE as {C, D, E}, DE|ABC as {D, E}.
    const auto ms = splits(m, taxa_of(m));
    const Split* ab = nullptr;
    const Split* de = nullptr;
    for (const auto& x : ms) {
      if (x.side[0] == 0b11100) ab = &x;
      if (x.side[0] == 0b11000) de = &x;
    }
    s.require("the majority-rule tree of two ((A,B),C,(D,E)) and one ((A,C),B,(D,E)) is ((A,B),C,(D,E))",
              robinson_foulds(m, expect) == 0 && ab && de && ms.size() == 2);
    s.check("its split AB has frequency 2/3", ab ? ab->support - 2.0 / 3 : 1, 1e-15);
    s.check("its split DE has frequency 1", de ? de->support - 1 : 1, 0);
    s.check("and AB's length is the mean of AB's lengths, 3", ab ? ab->length - 3 : 1, 1e-15);
    const auto freq = split_frequencies(expect, ts);
    s.check("the split frequencies of a reference in the sample agree", freq[static_cast<std::size_t>(expect[expect.find_leaf("A")].parent)] - 2.0 / 3, 1e-15);
  }

  // ---- layout ----------------------------------------------------------------------------------------
  {
    const Tree l = ladderize(t12);
    const Layout lay = layout(l);
    const auto dep = depths(l);
    bool ys = true, between = true;
    for (std::size_t k = 0; k < lay.leaf_order.size(); k++) ys = ys && lay.y[static_cast<std::size_t>(lay.leaf_order[k])] == static_cast<double>(k);
    for (std::size_t v = 0; v < l.size(); v++) {
      const auto& c = l.nodes[v].children;
      if (c.empty()) continue;
      between = between && lay.y[v] >= lay.y[static_cast<std::size_t>(c.front())] && lay.y[v] <= lay.y[static_cast<std::size_t>(c.back())] &&
                lay.x[v] <= lay.x[static_cast<std::size_t>(c.front())];
    }
    s.require("leaves are ranked 0 to n − 1 top to bottom, every internal node between its children and left of them", ys && between);
    s.check("x is the distance from the root", max_difference(lay.x, dep), 1e-12);
    double unrooted = 0;
    for (std::size_t v = 0; v < l.size(); v++)
      if (static_cast<int>(v) != l.root) {
        const int p = l.nodes[v].parent;
        const double dx = lay.ux[v] - lay.ux[static_cast<std::size_t>(p)], dy = lay.uy[v] - lay.uy[static_cast<std::size_t>(p)];
        unrooted = std::max(unrooted, std::fabs(std::hypot(dx, dy) - l.nodes[v].length));
      }
    s.check("in the equal-angle layout every branch is drawn at its length", unrooted, 1e-12);
  }

  // ---- random trees ----------------------------------------------------------------------------------
  {
    // Five leaves: 15 unrooted topologies, each with probability 1/15.
    std::map<std::vector<Bits>, int> count;
    const int N = 15000;
    std::vector<std::string> five{"t1", "t2", "t3", "t4", "t5"};
    for (int k = 0; k < N; k++) {
      const Tree r = random_unrooted_tree(5, 1000 + static_cast<std::uint64_t>(k));
      std::vector<Bits> key;
      for (const auto& sp : splits(r, five)) key.push_back(sp.side);
      count[key]++;
    }
    double chi = 0;
    for (const auto& [key, c] : count) chi += (c - N / 15.0) * (c - N / 15.0) / (N / 15.0);
    s.require("random unrooted trees on 5 leaves take all 15 topologies", count.size() == 15);
    s.check("uniformly: χ² with 14 degrees of freedom below its 0.999 quantile, 36.1", chi, 36.1);
  }
  {
    const int N = 4000, n = 10;
    double h = 0, h2 = 0, y = 0, y2 = 0;
    for (int k = 0; k < N; k++) {
      const Tree c = coalescent_tree(n, 77 + static_cast<std::uint64_t>(k));
      const auto dc = depths(c);
      const double hc = dc[static_cast<std::size_t>(c.leaves()[0])];
      h += hc;
      h2 += hc * hc;
      const Tree yt = yule_tree(n, 5000 + static_cast<std::uint64_t>(k), 2.0);
      const auto dy = depths(yt);
      const double hy = dy[static_cast<std::size_t>(yt.leaves()[0])];
      y += hy;
      y2 += hy * hy;
    }
    h /= N;
    y /= N;
    const double sh = std::sqrt((h2 / N - h * h) / N), sy = std::sqrt((y2 / N - y * y) / N);
    double ey = 0;
    for (int k = 2; k <= n; k++) ey += 1.0 / (2.0 * k);
    s.check("the coalescent's mean height is 2(1 − 1/n) = 1.8 for n = 10, within 4 standard errors", (h - 1.8) / sh, 4);
    s.check("the Yule tree's, Σ_{k=2}^{10} 1/(2k) for λ = 2, within 4 standard errors", (y - ey) / sy, 4);
    const Tree c = coalescent_tree(n, 3);
    const auto dc = depths(c);
    double spread = 0;
    for (int v : c.leaves()) spread = std::max(spread, std::fabs(dc[static_cast<std::size_t>(v)] - dc[static_cast<std::size_t>(c.leaves()[0])]));
    s.check("coalescent trees are ultrametric", spread, 1e-12);
  }

  return s.done();
}
