// Checks of distances and of the trees built from them: Kimura's worked example (1980, p. 116
// and p. 117); the closed forms of Jin and Nei (1990) and of the Γ JC69 distance; TN93 reducing
// to K80 and F81 to JC69 when they should; pairwise deletion; saturation reported, never a NaN;
// JC69, K80, TN93, LogDet and Γ-JC69 distances on sequences simulated at a known length, within
// about four standard errors; neighbour joining on the worked example of Wikipedia's article and
// on Saitou and Nei's (1987) Table 1, and on random additive trees, giving back the tree and its
// lengths exactly; BIONJ likewise; UPGMA on Wikipedia's worked example (heights 8.5, 11, 14,
// 16.5), ultrametric always, and exact on ultrametric matrices.

#include <evolution/phylogeny/distance.hpp>
#include <evolution/phylogeny/distance_trees.hpp>
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

double max_difference(const std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size()) return INFINITY;
  double m = 0;
  for (std::size_t i = 0; i < a.size(); i++) m = std::max(m, std::fabs(a[i] - b[i]));
  return m;
}

// A tree of two leaves at distance d.
Tree pair_tree(double d) {
  Tree t;
  const int r = t.add(-1);
  t.add(r, d / 2, "x");
  t.add(r, d / 2, "y");
  return t;
}

double distance_of(const Alignment& a, DistanceOptions o) { return distances(a, o)(0, 1); }

double ultrametric_spread(const Tree& t) {
  const auto d = depths(t);
  double lo = INFINITY, hi = -INFINITY;
  for (int v : t.leaves()) {
    lo = std::min(lo, d[static_cast<std::size_t>(v)]);
    hi = std::max(hi, d[static_cast<std::size_t>(v)]);
  }
  return hi - lo;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- formulas against their papers -----------------------------------------------------------------
  // Kimura (1980), p. 116: rabbit and chicken β-globin, 438 sites, 58 transitions and 63
  // transversions, K = 0.348; p. 117, the three codon positions of 146 sites each, K1 = 0.300,
  // K2 = 0.195, K3 = 0.635 (published to three places, the third truncated).
  s.check("K80 on Kimura's example (58 transitions, 63 transversions of 438) is his 0.348", k80_distance(58.0 / 438, 63.0 / 438).d - 0.348, 6e-4);
  s.check("K80 at the first codon positions (15, 21 of 146) is his K1 = 0.300", k80_distance(15.0 / 146, 21.0 / 146).d - 0.300, 1e-3);
  s.check("K80 at the second (7, 18 of 146) is his K2 = 0.195", k80_distance(7.0 / 146, 18.0 / 146).d - 0.195, 1e-3);
  s.check("K80 at the third (36, 24 of 146) is his K3 = 0.635", k80_distance(36.0 / 146, 24.0 / 146).d - 0.635, 1e-3);
  // The Γ JC69 distance (3/4) α [(1 − 4p/3)^{−1/α} − 1]: p = 0.3, α = 0.5 gives exactly 2/3.
  s.check("the Γ JC69 distance at p = 0.3, α = 0.5 is 2/3", jc69_distance(0.3, 0.5).d - 2.0 / 3, 1e-15);
  s.check("at α = 1 it is p/(1 − 4p/3) (Jin and Nei 1990, A11)", jc69_distance(0.2, 1).d - 0.2 / (1 - 0.8 / 3), 1e-15);
  {
    const double P = 0.12, Q = 0.07;
    s.check("the Γ K80 distance at α = 1 is (2P + Q)/(2(1 − 2P − Q)) + Q/(2(1 − 2Q)) (Jin and Nei 1990, A8)",
            k80_distance(P, Q, 1).d - ((2 * P + Q) / (2 * (1 - 2 * P - Q)) + Q / (2 * (1 - 2 * Q))), 1e-15);
    s.check("as α grows the Γ K80 distance tends to K80's", k80_distance(P, Q, 1e7).d - k80_distance(P, Q).d, 1e-7);
    const std::array<double, 4> uniform{0.25, 0.25, 0.25, 0.25};
    s.check("TN93 with equal frequencies and P1 = P2 is K80", tn93_distance(P / 2, P / 2, Q, uniform).d - k80_distance(P, Q).d, 1e-15);
    s.check("and so with Γ rates (Tamura and Nei 1993, eq. 15, against Jin and Nei 1990, A4)",
            tn93_distance(P / 2, P / 2, Q, uniform, 0.7).d - k80_distance(P, Q, 0.7).d, 1e-14);
    s.check("F81 with equal frequencies is JC69", f81_distance(0.2, uniform).d - jc69_distance(0.2).d, 1e-15);
    // TN93 with π_G = 0 and no A↔G difference: the purine term vanishes, its limit.
    const auto t0 = tn93_distance(0, 0.05, 0.02, {0.4, 0.3, 0, 0.3});
    s.require("TN93 with a base absent stays defined, finite and positive", !t0.saturated && std::isfinite(t0.d) && t0.d > 0);
  }
  {
    const auto j = jc69_distance(0.8), k = k80_distance(0.3, 0.45);
    s.require("beyond saturation JC69 and K80 report it", j.saturated && k.saturated);
    Alignment a = alignment_from_letters({"a", "b", "c"}, {"ACGTACGTAC", "CATGCATGCA", "----------"});
    DistanceOptions o;
    const auto m = distances(a, o);
    s.require("a saturated pair and a pair with no common site get the cap, flagged, never a NaN",
              m.saturated[0 * 3 + 1] && m.saturated[0 * 3 + 2] && m(0, 1) == o.cap && m(0, 2) == o.cap && m.saturated_pairs() == 3);
  }
  {
    Alignment a = alignment_from_letters({"a", "b"}, {"ACGTACGTAC", "ACGTTCGTA-"});
    DistanceOptions o;
    o.model = DistanceModel::p;
    const auto m = distances(a, o);
    s.check("gaps are deleted pair by pair: 1 difference in the 9 sites both have", m(0, 1) - 1.0 / 9, 1e-15);
    s.require("and the matrix says how many sites it compared", m.sites[1] == 9);
  }

  // ---- distances on simulated pairs -------------------------------------------------------------------
  {
    const double d = 0.3;
    const std::size_t L = 100000;
    DistanceOptions o;
    // JC69: Var(d̂) = p(1 − p) / (L (1 − 4p/3)²), about 0.002² here.
    const Alignment jc = simulate_alignment(pair_tree(d), jc69(), Rates{}, L, 11);
    s.check("JC69 recovers d = 0.3 from 100,000 simulated sites (4 standard errors 0.008)", distance_of(jc, o) - d, 0.008);
    o.model = DistanceModel::logdet;
    s.check("LogDet (paralinear) equals it under JC69", distance_of(jc, o) - d, 0.008);
    o.model = DistanceModel::k80;
    const Alignment k = simulate_alignment(pair_tree(d), k80(4), Rates{}, L, 12);
    s.check("K80 recovers d under K80 with κ = 4", distance_of(k, o) - d, 0.01);
    o.model = DistanceModel::jc69;
    s.require("where JC69 underestimates it", distance_of(k, o) < d);
    o.model = DistanceModel::tn93;
    const Alignment tn = simulate_alignment(pair_tree(d), tn93(6, 3, {0.35, 0.15, 0.2, 0.3}), Rates{}, L, 13);
    s.check("TN93 recovers d under TN93 with κ₁ = 6, κ₂ = 3 and unequal frequencies", distance_of(tn, o) - d, 0.01);
    // Under a stationary model det P(t) = e^{tr(Q) t}, so the paralinear distance is −tr(Q) t / 4,
    // which is t only when the frequencies are equal.
    o.model = DistanceModel::logdet;
    const Model tm = tn93(6, 3, {0.35, 0.15, 0.2, 0.3});
    const double trace = tm.Q[0] + tm.Q[5] + tm.Q[10] + tm.Q[15];
    s.check("LogDet is −tr(Q) d / 4 under TN93 (det P(t) = e^{tr(Q) t})", distance_of(tn, o) + trace * d / 4, 0.012);
    // Continuous Γ rates with α = 0.5: the Γ distance recovers d, the uncorrected one does not.
    SiteRates g;
    g.continuous_alpha = 0.5;
    const Alignment ga = simulate_alignment(pair_tree(d), jc69(), g, 200000, 14);
    o.model = DistanceModel::jc69;
    o.alpha = 0.5;
    s.check("the Γ JC69 distance recovers d under Γ(0.5) rates", distance_of(ga, o) - d, 0.015);
    o.alpha = INFINITY;
    s.require("the JC69 distance without Γ falls short of it by more than 0.03", d - distance_of(ga, o) > 0.03);
  }

  // ---- neighbour joining --------------------------------------------------------------------------------
  {
    // Wikipedia, "Neighbor joining", revision 1373946134 (2026-09-08), the worked example.
    const DistanceMatrix m({"a", "b", "c", "d", "e"}, {0, 5, 9, 9, 8, 5, 0, 10, 10, 9, 9, 10, 0, 8, 7, 9, 10, 8, 0, 3, 8, 9, 7, 3, 0});
    const Tree t = neighbour_joining(m);
    const auto leaf = [&](const char* x) { return t[t.find_leaf(x)].length; };
    std::vector<double> internal;
    for (std::size_t v = 0; v < t.size(); v++)
      if (!t.nodes[v].children.empty() && static_cast<int>(v) != t.root) internal.push_back(t.nodes[v].length);
    std::sort(internal.begin(), internal.end());
    s.require("NJ on Wikipedia's example gives ((a, b), c, (d, e))", robinson_foulds(t, read_newick("((a,b),c,(d,e));").trees[0]) == 0);
    s.check("with its branches a 2, b 3, c 4, d 2, e 1", std::fabs(leaf("a") - 2) + std::fabs(leaf("b") - 3) + std::fabs(leaf("c") - 4) +
                                                           std::fabs(leaf("d") - 2) + std::fabs(leaf("e") - 1), 1e-12);
    s.check("and internal branches 3 (u–v) and 2 (v–w)", internal.size() == 2 ? std::fabs(internal[0] - 2) + std::fabs(internal[1] - 3) : 1, 1e-12);
    s.check("its leaf distances are the matrix, which is additive", max_difference(leaf_distances(t, m.names), m.d), 1e-12);
  }
  {
    // Saitou and Nei (1987), Table 1: the additive matrix of the tree of their Fig. 1 (L_DF = 1, as
    // the matrix requires; the text's 2 on p. 411 is a misprint).
    const std::vector<double> lower{7, 8, 5, 11, 8, 5, 13, 10, 7, 8, 16, 13, 10, 11, 5, 13, 10, 7, 8, 6, 9, 17, 14, 11, 12, 10, 13, 8};
    std::vector<std::string> names{"1", "2", "3", "4", "5", "6", "7", "8"};
    std::vector<double> full(64, 0);
    std::size_t k = 0;
    for (std::size_t i = 1; i < 8; i++)
      for (std::size_t j = 0; j < i; j++) full[i * 8 + j] = full[j * 8 + i] = lower[k++];
    const DistanceMatrix m(names, full);
    const Tree truth = read_newick("('1':5,'2':2,('3':1,('4':3,(('5':1,'6':4):2,('7':2,'8':6):1):2):1):2);").trees[0];
    s.check("the tree of Saitou and Nei's Fig. 1 has the distances of their Table 1", max_difference(leaf_distances(truth, names), full), 1e-12);
    const Tree t = neighbour_joining(m);
    s.require("NJ gives back the tree of Saitou and Nei's Fig. 1 from their Table 1", robinson_foulds(t, truth) == 0);
    s.check("with its lengths: the leaf distances are the table's", max_difference(leaf_distances(t, names), full), 1e-12);
    s.require("BIONJ gives it back too", robinson_foulds(bionj(m), truth) == 0 && max_difference(leaf_distances(bionj(m), names), full) < 1e-12);
  }
  {
    bool topo = true, bio = true;
    double worst = 0, worst_bio = 0;
    for (std::uint64_t seed = 1; seed <= 6; seed++) {
      const Tree truth = random_unrooted_tree(30, seed, 0.1, 0.005);
      const auto taxa = taxa_of(truth);
      const DistanceMatrix m(taxa, leaf_distances(truth, taxa));
      const Tree t = neighbour_joining(m), b = bionj(m);
      topo = topo && robinson_foulds(t, truth) == 0;
      bio = bio && robinson_foulds(b, truth) == 0;
      worst = std::max(worst, max_difference(leaf_distances(t, taxa), m.d));
      worst_bio = std::max(worst_bio, max_difference(leaf_distances(b, taxa), m.d));
    }
    s.require("NJ gives back six random additive trees on 30 leaves exactly (RF 0)", topo);
    s.check("with their branch lengths", worst, 1e-12);
    s.require("so does BIONJ", bio);
    s.check("with their branch lengths", worst_bio, 1e-12);
  }
  {
    // Atteson (1999): NJ is right when every entry is within half the shortest branch. Perturb by 0.45 of it.
    const Tree truth = random_unrooted_tree(25, 99, 0.1, 0.02);
    const auto taxa = taxa_of(truth);
    auto d = leaf_distances(truth, taxa);
    double shortest = INFINITY;
    for (std::size_t v = 0; v < truth.size(); v++)
      if (static_cast<int>(v) != truth.root) shortest = std::min(shortest, truth.nodes[v].length);
    calx::random::Stream r(5, 0x4e4f4953, 0, 0);
    for (std::size_t i = 0; i < 25; i++)
      for (std::size_t j = i + 1; j < 25; j++) d[i * 25 + j] = d[j * 25 + i] = d[i * 25 + j] + (r.uniform() < 0.5 ? -0.45 : 0.45) * shortest;
    s.require("NJ keeps the topology when every distance is off by 0.45 of the shortest branch (within Atteson's radius, 1/2)",
              robinson_foulds(neighbour_joining(DistanceMatrix(taxa, d)), truth) == 0);
  }

  // ---- UPGMA -----------------------------------------------------------------------------------------------
  {
    // Wikipedia, "UPGMA", the worked example (5S rRNA of five bacteria): heights 8.5, 11, 14, 16.5,
    // which follow by hand from the definition.
    const DistanceMatrix m({"a", "b", "c", "d", "e"}, {0, 17, 21, 31, 23, 17, 0, 30, 34, 21, 21, 30, 0, 28, 39, 31, 34, 28, 0, 43, 23, 21, 39, 43, 0});
    const Tree t = upgma(m);
    const Tree expect = read_newick("(((a:8.5,b:8.5):2.5,e:11):5.5,(c:14,d:14):2.5);").trees[0];
    s.require("UPGMA on Wikipedia's example gives (((a, b), e), (c, d))", robinson_foulds(t, expect) == 0 && t.rooted);
    s.check("with the heights 8.5, 11, 14 and 16.5", max_difference(leaf_distances(t, m.names), leaf_distances(expect, m.names)), 1e-12);
    s.check("and its root at 16.5", depths(t)[static_cast<std::size_t>(t.find_leaf("a"))] - 16.5, 1e-12);
    const DistanceMatrix nj({"a", "b", "c", "d", "e"}, {0, 5, 9, 9, 8, 5, 0, 10, 10, 9, 9, 10, 0, 8, 7, 9, 10, 8, 0, 3, 8, 9, 7, 3, 0});
    s.check("UPGMA's tree is ultrametric even on a matrix that is not", ultrametric_spread(upgma(nj)), 1e-12);
    s.check("and so is WPGMA's", ultrametric_spread(upgma(nj, Linkage::weighted)), 1e-12);
  }
  {
    bool topo = true;
    double worst = 0;
    for (std::uint64_t seed = 1; seed <= 6; seed++) {
      const Tree truth = coalescent_tree(25, seed);
      const auto taxa = taxa_of(truth);
      const DistanceMatrix m(taxa, leaf_distances(truth, taxa));
      const Tree t = upgma(m);
      topo = topo && robinson_foulds(t, truth) == 0;
      worst = std::max(worst, max_difference(leaf_distances(t, taxa), m.d));
    }
    s.require("UPGMA gives back six coalescent trees on 25 leaves from their ultrametrics (RF 0)", topo);
    s.check("with their heights", worst, 1e-12);
  }

  return s.done();
}
