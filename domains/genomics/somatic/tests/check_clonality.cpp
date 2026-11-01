// Checks of clonality: the cancer cell fraction of a mutation by hand (McGranahan et al. 2015's
// n_mut, multiplicity and the grid posterior with its interval); clusters of CCFs from simulated
// reads; the trees of clones, (K − 1)! examined, the sum rule and the pigeonhole principle; the
// mutational burden.

#include <genomics/somatic/clonality.hpp>

#include <calx/random/stream.hpp>
#include <calx/special_functions/gamma.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace genomics;

int main() {
  calx::verification::Suite s;

  // ---- one mutation by hand ----------------------------------------------------------------------------
  {
    // ρ = 0.6, a segment 2 + 1, 46 of 100 reads: n_mut = 0.46 (0.6 · 3 + 0.8)/0.6 = 1.9933, m = 2.
    somatic::SomaticCall c;
    c.tumour_depth = 100;
    c.tumour_alt = 46;
    somatic::CnSegment g;
    g.contig = 0;
    g.begin = 0;
    g.end = 1000;
    g.n_a = 2;
    g.n_b = 1;
    const auto m = somatic::cancer_cell_fractions({c}, {g}, 0.6);
    const double n_mut = 0.46 * (0.6 * 3 + 0.8) / 0.6;
    s.require("one mutation, one estimate", m.size() == 1);
    s.check("n_mut = VAF (ρN + 2(1 − ρ))/ρ", m[0].n_mut - n_mut, 1e-12);
    s.require("multiplicity round(n_mut) = 2, at most the major allele's 2", m[0].multiplicity == 2);
    s.check("CCF = n_mut / m", m[0].ccf - n_mut / 2, 1e-12);
    // The grid posterior and its interval by brute force.
    std::vector<double> w(100);
    double total = 0;
    for (int i = 1; i <= 100; i++) {
      const double p = 0.6 * 2 * (i / 100.0) / 2.6;
      w[i - 1] = std::exp(calx::special_functions::log_choose(100, 46) + 46 * std::log(p) + 54 * std::log1p(-p));
      total += w[i - 1];
    }
    double cum = 0, lo = 0, hi = 0;
    for (int i = 1; i <= 100; i++) {
      cum += w[i - 1] / total;
      if (lo == 0 && cum >= 0.025) lo = i / 100.0;
      if (hi == 0 && cum >= 0.975) hi = i / 100.0;
    }
    s.require("the 95% interval of the grid posterior", std::fabs(m[0].ccf_low - lo) < 1e-12 && std::fabs(m[0].ccf_high - hi) < 1e-12,
              std::to_string(m[0].ccf_low) + "–" + std::to_string(m[0].ccf_high));
    s.require("clonal: the interval reaches 1", m[0].clonal == (hi >= 1));
    // A subclonal one: 12 of 100 at 1 + 1, ρ = 0.6: f(φ = 1) = 0.3, so φ̂ = 0.4.
    c.tumour_alt = 12;
    g.n_a = 1;
    const auto q = somatic::cancer_cell_fractions({c}, {g}, 0.6);
    s.check("a subclonal mutation: CCF 0.4", q[0].ccf - 0.4, 1e-12);
    s.require("a subclonal mutation is not clonal", !q[0].clonal, std::to_string(q[0].ccf_high));
  }

  // ---- clusters ----------------------------------------------------------------------------------------
  {
    calx::random::Stream r(9, 0, 0, 0);
    const double rho = 0.7;
    const std::vector<double> ccf{1.0, 0.45};
    std::vector<somatic::SomaticCall> calls;
    std::vector<somatic::CnSegment> segs;
    for (int i = 0; i < 400; i++) {
      somatic::CnSegment g;
      g.begin = static_cast<std::size_t>(i) * 1000;
      g.end = g.begin + 1000;
      g.n_a = r.bernoulli(0.25) ? 2 : 1;
      g.n_b = 1;
      segs.push_back(g);
      const int j = i % 3 == 0 ? 1 : 0;
      somatic::SomaticCall c;
      c.position = g.begin + 500;
      c.tumour_depth = static_cast<std::size_t>(r.poisson(70)) + 10;
      const double f = rho * ccf[static_cast<std::size_t>(j)] / (2 * (1 - rho) + rho * (g.n_a + g.n_b));
      c.tumour_alt = static_cast<std::size_t>(r.binomial(static_cast<std::int64_t>(c.tumour_depth), f));
      calls.push_back(c);
    }
    auto m = somatic::cancer_cell_fractions(calls, segs, rho);
    const auto cl = somatic::cluster_mutations(m, 5);
    std::printf("      clusters:");
    for (const auto& c : cl.clones) std::printf(" CCF %.3f (%zu)", c.ccf, c.mutations);
    std::printf("\n");
    s.require("two clones found", cl.clones.size() == 2, std::to_string(cl.clones.size()));
    s.check("the clonal cluster's CCF", cl.clones.size() == 2 ? cl.clones[0].ccf - 1.0 : 1, 0.03);
    s.check("the subclone's CCF", cl.clones.size() == 2 ? cl.clones[1].ccf - 0.45 : 1, 0.03);
  }

  // ---- trees -----------------------------------------------------------------------------------------
  {
    const auto t = somatic::clone_trees({1.0, 0.6, 0.35, 0.3});
    s.require("four clones: (K − 1)! = 6 trees examined", t.examined == 6, std::to_string(t.examined));
    // The parent vectors (the 0.6 clone under the trunk; the 0.35 under either; the 0.3 under any)
    // satisfying the sum rule within 0.05, counted directly; the star (1.25 under 1.0) is not one.
    std::size_t brute = 0;
    const std::vector<double> c{1.0, 0.6, 0.35, 0.3};
    for (int p2 = 0; p2 < 2; p2++)
      for (int p3 = 0; p3 < 3; p3++) {
        std::vector<double> ch(4, 0);
        ch[0] += c[1];
        ch[static_cast<std::size_t>(p2)] += c[2];
        ch[static_cast<std::size_t>(p3)] += c[3];
        bool ok = true;
        for (int k = 0; k < 4; k++) ok = ok && ch[static_cast<std::size_t>(k)] <= c[static_cast<std::size_t>(k)] + 0.05;
        brute += ok;
      }
    s.require("the trees satisfying the sum rule are those counted by brute force", t.trees.size() == brute, std::to_string(t.trees.size()) + " of " + std::to_string(brute));
    const auto u = somatic::clone_trees({1.0, 0.7, 0.5});
    bool nested = !u.forced.empty();
    for (const auto& tree : u.trees) nested = nested && somatic::descends(tree, 1, 2);
    s.require("pigeonhole: 0.7 + 0.5 > 1 forces the 0.5 clone under the 0.7 one in every tree", nested && u.trees.size() == 1);
    const auto v = somatic::clone_trees({0.8, 0.3});
    s.require("a trunk of CCF 1 is put above clones that all fall short of it", v.clone.size() == 3 && v.clone[0] == -1 && v.ccf[0] == 1);
  }

  // ---- burden ------------------------------------------------------------------------------------------
  {
    somatic::SomaticCalls c;
    c.callable = 2'000'000;
    for (int i = 0; i < 10; i++) {
      somatic::SomaticCall x;
      x.indel = i < 3;
      if (i == 9) x.filters.push_back("strand_bias");
      c.calls.push_back(x);
    }
    const auto b = somatic::mutational_burden(c);
    s.require("burden: 6 substitutions and 3 indels pass, 4.5 per megabase over 2 Mb", b.substitutions == 6 && b.indels == 3 && std::fabs(b.per_megabase - 4.5) < 1e-12);
  }
  return s.done();
}
