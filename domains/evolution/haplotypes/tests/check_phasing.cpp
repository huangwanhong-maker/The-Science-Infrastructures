// Checks of the pair model of a diploid's genotypes: its O(N K²) likelihood and posterior
// against brute-force enumeration of every pair of copying paths, and its phasing (the joint
// maximum over paths and haplotypes) against enumeration of every pair of paths and every
// assignment of the heterozygous alleles; the O(K²) step against the dense O(K⁴) recursion;
// the posterior's symmetry in the two copies; switch errors counted by hand; and phasing on the
// coalescent's diploids, whose switch error rate is measured against random phase.

#include <evolution/haplotypes/phasing.hpp>
#include <evolution/population/coalescent.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace evolution::haplotypes;
using evolution::population::CoalescentSpec;
using evolution::population::coalescent;

namespace {

std::string fmt(double x) { return calx::verification::exponential(x); }

Haplotypes random_panel(std::size_t K, std::size_t N, calx::random::Stream& r) {
  Haplotypes h;
  h.sequences = K;
  h.sites = N;
  h.alleles.resize(K * N);
  for (auto& a : h.alleles) a = static_cast<std::int8_t>(r.bernoulli(0.4));
  double x = 0;
  for (std::size_t j = 0; j < N; j++) {
    h.positions.push_back(x);
    x += r.uniform() * 0.5;
  }
  return h;
}

struct Brute {
  double log_likelihood = NEG_INF, best = NEG_INF;
  std::vector<double> gamma;  // N × K × K
};

// Every pair of paths: the likelihood and marginals from the summed emissions, and the joint
// maximum over the paths and every assignment of alleles to the two copies.
Brute enumerate(const Haplotypes& panel, const Chain& c, const std::vector<std::int8_t>& g) {
  const std::size_t K = c.states, N = c.sites, KK = K * K;
  Brute b;
  b.gamma.assign(N * KK, 0);
  std::vector<std::uint32_t> pairs(N, 0);
  std::vector<std::size_t> het;
  for (std::size_t j = 0; j < N; j++)
    if (g[j] == 1) het.push_back(j);
  Phasing p;
  p.first.resize(N);
  p.second.resize(N);
  p.haplotype_first.resize(N);
  p.haplotype_second.resize(N);
  double total = 0;
  for (;;) {
    for (std::size_t j = 0; j < N; j++) {
      p.first[j] = pairs[j] / static_cast<std::uint32_t>(K);
      p.second[j] = pairs[j] % static_cast<std::uint32_t>(K);
    }
    // The transitions, and the summed emissions.
    double lp = -2 * std::log(static_cast<double>(K));
    for (std::size_t j = 0; j < N; j++) {
      if (j) lp += std::log(c.transition(j, p.first[j - 1], p.first[j])) + std::log(c.transition(j, p.second[j - 1], p.second[j]));
      if (g[j] < 0) continue;
      const double pa = panel(p.first[j], j) ? 1 - c.epsilon : c.epsilon, pb = panel(p.second[j], j) ? 1 - c.epsilon : c.epsilon;
      lp += std::log(g[j] == 0 ? (1 - pa) * (1 - pb) : g[j] == 2 ? pa * pb : pa * (1 - pb) + (1 - pa) * pb);
    }
    const double w = std::exp(lp);
    total += w;
    for (std::size_t j = 0; j < N; j++) b.gamma[j * KK + pairs[j]] += w;
    // Every assignment of the heterozygous alleles; homozygous and missing sites as phase() sets them.
    for (std::uint64_t mask = 0; mask < (1ull << het.size()); mask++) {
      for (std::size_t j = 0; j < N; j++) {
        if (g[j] == 0 || g[j] == 2) p.haplotype_first[j] = p.haplotype_second[j] = static_cast<std::int8_t>(g[j] / 2);
        else if (g[j] < 0) {
          p.haplotype_first[j] = panel(p.first[j], j);
          p.haplotype_second[j] = panel(p.second[j], j);
        }
      }
      for (std::size_t h = 0; h < het.size(); h++) {
        p.haplotype_first[het[h]] = static_cast<std::int8_t>((mask >> h) & 1);
        p.haplotype_second[het[h]] = static_cast<std::int8_t>(1 - p.haplotype_first[het[h]]);
      }
      b.best = std::max(b.best, phased_log_probability(panel, c, g, p));
    }
    std::size_t j = 0;
    while (j < N && ++pairs[j] == KK) pairs[j++] = 0;
    if (j == N) break;
  }
  b.log_likelihood = std::log(total);
  for (double& x : b.gamma) x /= total;
  return b;
}

// The panel of the first K sequences and the diploids of the rest, from one coalescent sample.
struct Split {
  Haplotypes panel;
  std::vector<std::vector<std::int8_t>> haplotypes;  // the targets' haplotypes, two per diploid
};
Split split(const Haplotypes& all, std::size_t K) {
  Split s;
  s.panel.sequences = K;
  s.panel.sites = all.sites;
  s.panel.positions = all.positions;
  s.haplotypes.assign(all.sequences - K, std::vector<std::int8_t>(all.sites));
  for (std::size_t j = 0; j < all.sites; j++) {
    const auto row = all.site(j);
    s.panel.alleles.insert(s.panel.alleles.end(), row.begin(), row.begin() + static_cast<std::ptrdiff_t>(K));
    for (std::size_t i = K; i < all.sequences; i++) s.haplotypes[i - K][j] = row[i];
  }
  return s;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- switch errors by hand ------------------------------------------------------------------------------
  {
    const std::vector<std::int8_t> t = {0, 1, 0, 1, 0}, u = {1, 0, 1, 0, 0}, x = {0, 1, 1, 0, 0}, y = {1, 0, 0, 1, 0};
    const SwitchErrors e = switch_errors(t, u, x, y);
    s.require("switch errors by hand: four heterozygous sites, one switch, rate 1/3", e.heterozygous == 4 && e.switches == 1 && std::fabs(e.rate() - 1.0 / 3) < 1e-15);
  }

  // ---- brute force on small panels -------------------------------------------------------------------------
  {
    calx::random::Stream r(21, 0, 0, 0);
    double e_ll = 0, e_dense = 0, e_post = 0, e_map = 0, e_self = 0;
    bool consistent = true;
    int cases = 0;
    for (std::size_t K : {2u, 3u})
      for (std::size_t N : {3u, 4u})
        for (int rep = 0; rep < 4; rep++) {
          const Haplotypes panel = random_panel(K, N, r);
          CopyingModel m;
          m.rho = 0.5 + 3 * r.uniform();
          m.epsilon = rep == 0 ? -1 : 0.01 + 0.2 * r.uniform();
          const Chain c = copying_chain(panel, m);
          std::vector<std::int8_t> g(N);
          for (auto& a : g) a = static_cast<std::int8_t>(r.below(3));
          if (rep % 2) g[r.below(N)] = -1;
          const Brute b = enumerate(panel, c, g);
          const PairPosterior p = pair_posterior(panel, c, g);
          const Phasing ph = phase(panel, c, g);
          e_ll = std::max(e_ll, calx::verification::relative(pair_log_likelihood(panel, c, g), b.log_likelihood));
          e_ll = std::max(e_ll, calx::verification::relative(p.log_likelihood, b.log_likelihood));
          e_dense = std::max(e_dense, calx::verification::relative(pair_log_likelihood_dense(panel, c, g), b.log_likelihood));
          e_post = std::max(e_post, calx::verification::max_difference(p.gamma, b.gamma));
          e_map = std::max(e_map, calx::verification::relative(ph.log_probability, b.best));
          e_self = std::max(e_self, calx::verification::relative(phased_log_probability(panel, c, g, ph), b.best));
          for (std::size_t j = 0; j < N; j++)
            consistent = consistent && (g[j] < 0 || ph.haplotype_first[j] + ph.haplotype_second[j] == g[j]);
          cases++;
        }
    const std::string n = " (" + std::to_string(cases) + " panels, K = 2, 3, N = 3, 4, against all K^{2N} pairs of paths)";
    s.check("the pair model's log P(g) = log Σ over every pair of paths" + n, e_ll, 1e-12, "relative");
    s.check("and the dense O(K⁴) recursion's", e_dense, 1e-12, "relative");
    s.check("the posterior over pairs = the enumerated marginals", e_post, 1e-12);
    s.check("phasing's log P(A, B, x, y) = the maximum over every pair of paths and assignment of alleles", e_map, 1e-12, "relative");
    s.check("and the paths and haplotypes it returns attain it", e_self, 1e-12, "relative");
    s.require("the phased haplotypes add up to the typed genotypes", consistent);
  }

  // ---- the O(K²) step against the dense O(K⁴) one; symmetry ------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 10;
    spec.theta = 20;
    spec.rho = 20;
    const Split sp = split(evolution::population::haplotypes(coalescent(spec, 22)), 8);
    std::vector<std::int8_t> g(sp.panel.sites);
    for (std::size_t j = 0; j < g.size(); j++) g[j] = static_cast<std::int8_t>(sp.haplotypes[0][j] + sp.haplotypes[1][j]);
    const Chain c = copying_chain(sp.panel, {20, -1, -1});
    s.check("the O(N K²) pair forward against the dense O(N K⁴) one (K = 8, " + std::to_string(g.size()) + " coalescent sites)",
            calx::verification::relative(pair_log_likelihood(sp.panel, c, g), pair_log_likelihood_dense(sp.panel, c, g)), 1e-12, "relative");
    const PairPosterior p = pair_posterior(sp.panel, c, g);
    double asym = 0;
    for (std::size_t j = 0; j < p.sites; j++)
      for (std::size_t a = 0; a < 8; a++)
        for (std::size_t b = 0; b < 8; b++) asym = std::max(asym, std::fabs(p.at(j, a, b) - p.at(j, b, a)));
    s.check("the posterior is symmetric in the two copies, γ(a, b) = γ(b, a)", asym, 1e-12);
  }

  // ---- phasing the coalescent's diploids --------------------------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 140;
    spec.theta = 40;
    spec.rho = 40;
    SwitchErrors total, random;
    std::size_t targets = 0, sites = 0;
    calx::random::Stream coin(23, 0, 0, 0);
    const auto t0 = std::chrono::steady_clock::now();
    for (std::uint32_t rep = 0; rep < 5; rep++) {
      const Split sp = split(evolution::population::haplotypes(coalescent(spec, 23, rep)), 100);
      const Chain c = copying_chain(sp.panel, {40, -1, -1});
      sites += sp.panel.sites;
      for (std::size_t i = 0; i + 1 < sp.haplotypes.size(); i += 2) {
        const auto& t = sp.haplotypes[i];
        const auto& u = sp.haplotypes[i + 1];
        std::vector<std::int8_t> g(t.size()), rx(t.size()), ry(t.size());
        for (std::size_t j = 0; j < g.size(); j++) {
          g[j] = static_cast<std::int8_t>(t[j] + u[j]);
          const bool flip = coin.bernoulli(0.5);
          rx[j] = flip ? u[j] : t[j];
          ry[j] = flip ? t[j] : u[j];
        }
        const Phasing ph = phase(sp.panel, c, g);
        total += switch_errors(t, u, ph.haplotype_first, ph.haplotype_second);
        random += switch_errors(t, u, rx, ry);
        targets++;
      }
    }
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    const std::string n = std::to_string(targets) + " diploids, panels of 100 haplotypes, θ = ρ = 40, " + std::to_string(sites / 5) + " sites on average";
    s.require("phasing the coalescent's diploids: switch error rate below 0.1 (" + n + ")", total.rate() < 0.1,
              fmt(total.rate()) + " (" + std::to_string(total.switches) + " switches among " + std::to_string(total.heterozygous) + " heterozygous sites); random phase " +
                  fmt(random.rate()));
    s.require("random phase switches about half the time", std::fabs(random.rate() - 0.5) < 0.05, fmt(random.rate()));
    s.require("the " + std::to_string(targets) + " diploids phased in under 10 s (O(N K²) each)", ms < 10000, fmt(ms) + " ms");
  }

  // ---- diploids drawn from the model ---------------------------------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 60;
    spec.theta = 40;
    spec.rho = 40;
    const Haplotypes panel = evolution::population::haplotypes(coalescent(spec, 24));
    CopyingModel cm;
    cm.rho = 20;
    cm.epsilon = 0.002;
    const Chain c = copying_chain(panel, cm);
    SwitchErrors total;
    for (std::uint32_t rep = 0; rep < 40; rep++) {
      calx::random::Stream r(24, MOSAIC, rep, 0);
      const Mosaic a = sample_mosaic(panel, c, r), b = sample_mosaic(panel, c, r);
      std::vector<std::int8_t> g(a.alleles.size());
      for (std::size_t j = 0; j < g.size(); j++) g[j] = static_cast<std::int8_t>(a.alleles[j] + b.alleles[j]);
      const Phasing ph = phase(panel, c, g);
      total += switch_errors(a.alleles, b.alleles, ph.haplotype_first, ph.haplotype_second);
    }
    s.require("diploids drawn from the model itself are phased with a switch error rate below 0.05", total.rate() < 0.05,
              fmt(total.rate()) + " over " + std::to_string(total.heterozygous) + " heterozygous sites");
  }
  return s.done();
}
