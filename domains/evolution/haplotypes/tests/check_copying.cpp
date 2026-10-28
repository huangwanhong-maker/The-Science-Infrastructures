// Checks of the Li and Stephens copying model: the O(N K) forward, backward and posterior
// against brute-force enumeration of all K^N copying paths on small panels, the Viterbi path
// against the best of them; the O(N K) forward against the dense O(N K²) recursion on the
// coalescent's haplotypes; the whole model against the genomics domain's HMM on an evenly
// spaced map, where the copying chain is homogeneous; and the posterior's calibration on
// targets drawn from the model itself, by the identity E[γ_j(Z_j)] = E[Σ_k γ_j(k)²] (the true
// state is distributed as the posterior given the data).

#include <evolution/haplotypes/copying.hpp>
#include <evolution/population/coalescent.hpp>

#include <genomics/hmm/model.hpp>

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

// A panel of K haplotypes at N sites with random alleles and a random map.
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
  std::vector<double> gamma;  // N × K
};

// Every one of the K^N paths, its joint probability with x.
Brute enumerate(const Haplotypes& panel, const Chain& c, const std::vector<std::int8_t>& x) {
  const std::size_t K = c.states, N = c.sites;
  Brute b;
  std::vector<std::uint32_t> path(N, 0);
  std::vector<double> mass(N * K, 0);
  double total = 0;
  for (;;) {
    const double lp = path_log_probability(panel, c, x, path);
    const double p = std::exp(lp);
    total += p;
    b.best = std::max(b.best, lp);
    for (std::size_t j = 0; j < N; j++) mass[j * K + path[j]] += p;
    std::size_t j = 0;
    while (j < N && ++path[j] == K) path[j++] = 0;
    if (j == N) break;
  }
  b.log_likelihood = std::log(total);
  b.gamma = mass;
  for (double& g : b.gamma) g /= total;
  return b;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- the parameters ----------------------------------------------------------------------------------
  s.check("θ̃(2) = 1 and θ̃(3) = 2/3 (Σ_{m=1}^{K−1} 1/m inverted)", std::fabs(li_stephens_theta(2) - 1) + std::fabs(li_stephens_theta(3) - 2.0 / 3), 1e-15);
  s.check("ε = θ/(2(K + θ)): 1/6 at K = 2 with θ̃", std::fabs(miscopying(li_stephens_theta(2), 2) - 1.0 / 6), 1e-15);

  // ---- brute force on small panels ----------------------------------------------------------------------
  {
    calx::random::Stream r(1, 0, 0, 0);
    double e_fwd = 0, e_bwd = 0, e_lowmem = 0, e_dense = 0, e_post = 0, e_vit = 0, e_path = 0;
    int cases = 0;
    for (std::size_t K : {2u, 3u, 4u})
      for (std::size_t N : {4u, 5u, 6u})
        for (int rep = 0; rep < 4; rep++) {
          const Haplotypes panel = random_panel(K, N, r);
          CopyingModel m;
          m.rho = 0.5 + 3 * r.uniform();
          m.epsilon = rep == 0 ? -1 : 0.01 + 0.2 * r.uniform();  // the first with θ̃
          const Chain c = copying_chain(panel, m);
          std::vector<std::int8_t> x(N);
          for (auto& a : x) a = static_cast<std::int8_t>(r.bernoulli(0.4));
          if (rep % 2) x[r.below(N)] = -1;  // a missing allele
          const Brute b = enumerate(panel, c, x);
          const Forward f = forward(panel, c, x);
          const Backward bw = backward(panel, c, x);
          const Posterior p = posterior(panel, c, x);
          const Path v = viterbi(panel, c, x);
          e_fwd = std::max(e_fwd, calx::verification::relative(f.log_likelihood, b.log_likelihood));
          e_bwd = std::max(e_bwd, calx::verification::relative(bw.log_likelihood, b.log_likelihood));
          e_lowmem = std::max(e_lowmem, calx::verification::relative(log_likelihood(panel, c, x), b.log_likelihood));
          e_dense = std::max(e_dense, calx::verification::relative(forward_dense(panel, c, x), b.log_likelihood));
          e_post = std::max(e_post, calx::verification::max_difference(p.gamma, b.gamma));
          e_vit = std::max(e_vit, calx::verification::relative(v.log_probability, b.best));
          e_path = std::max(e_path, calx::verification::relative(path_log_probability(panel, c, x, v.states), b.best));
          cases++;
        }
    const std::string n = " (" + std::to_string(cases) + " panels, K = 2 … 4, N = 4 … 6, against all K^N paths)";
    s.check("forward log P(x) = log Σ over every path" + n, e_fwd, 1e-12, "relative");
    s.check("backward computes the same log P(x) on its own", e_bwd, 1e-12, "relative");
    s.check("and so does log_likelihood in O(K) memory", e_lowmem, 1e-12, "relative");
    s.check("and the dense O(K²) recursion", e_dense, 1e-12, "relative");
    s.check("posterior γ_j(k) = the enumerated marginals", e_post, 1e-12);
    s.check("Viterbi's log P(x, Z*) = the best path's", e_vit, 1e-12, "relative");
    s.check("and the path it returns attains it (ties among identical haplotypes allowed)", e_path, 1e-12, "relative");
  }

  // ---- a target no path can emit ---------------------------------------------------------------------------
  {
    Haplotypes panel;
    panel.sequences = 2;
    panel.sites = 3;
    panel.alleles = {0, 0, 1, 0, 0, 1};
    panel.positions = {0, 1, 2};
    CopyingModel m;
    m.epsilon = 0;
    const Chain c = copying_chain(panel, m);
    const std::vector<std::int8_t> x = {1, 0, 0};  // allele 1 at site 0, which no haplotype carries
    const Path v = viterbi(panel, c, x);
    s.require("with ε = 0 an allele no haplotype carries has log P = −∞, an empty path and zero posteriors",
              forward(panel, c, x).log_likelihood == NEG_INF && v.states.empty() && v.log_probability == NEG_INF && posterior(panel, c, x).gamma[0] == 0);
  }

  // ---- O(K) against O(K²) on the coalescent's haplotypes ------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 61;
    spec.theta = 40;
    spec.rho = 40;
    const Haplotypes all = evolution::population::haplotypes(coalescent(spec, 2));
    Haplotypes panel;
    panel.sequences = 60;
    panel.sites = all.sites;
    panel.positions = all.positions;
    std::vector<std::int8_t> x(all.sites);
    calx::random::Stream r(2, 0, 0, 0);
    for (std::size_t j = 0; j < all.sites; j++) {
      const auto row = all.site(j);
      panel.alleles.insert(panel.alleles.end(), row.begin(), row.begin() + 60);
      x[j] = r.bernoulli(0.2) ? -1 : row[60];
    }
    const Chain c = copying_chain(panel, {40, -1, -1});
    const double fast = log_likelihood(panel, c, x), dense = forward_dense(panel, c, x);
    s.check("the O(N K) forward against the dense O(N K²) one on " + std::to_string(all.sites) + " coalescent sites, K = 60, ρ = 40, a fifth missing",
            calx::verification::relative(fast, dense), 1e-12, "relative");
    const Posterior p = posterior(panel, c, x);
    double worst = 0;
    for (std::size_t j = 0; j < p.sites; j++) {
      double t = 0;
      for (const double g : p.site(j)) t += g;
      worst = std::max(worst, std::fabs(t - 1));
    }
    s.check("each site's posterior sums to one", worst, 1e-12);
  }

  // ---- the genomics domain's HMM on an evenly spaced map -------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 31;
    spec.theta = 30;
    spec.rho = 30;
    const Haplotypes all = evolution::population::haplotypes(coalescent(spec, 3));
    const std::size_t K = 30, N = std::min<std::size_t>(all.sites, 120);
    Haplotypes panel;
    panel.sequences = K;
    panel.sites = N;
    std::vector<double> map(N);
    std::vector<std::int8_t> x(N);
    calx::random::Stream r(3, 0, 0, 0);
    for (std::size_t j = 0; j < N; j++) {
      const auto row = all.site(j);
      panel.alleles.insert(panel.alleles.end(), row.begin(), row.begin() + K);
      panel.positions.push_back(all.positions[j]);
      map[j] = static_cast<double>(j) / static_cast<double>(N);
      x[j] = r.bernoulli(0.1) ? -1 : row[K];
    }
    CopyingModel cm;
    cm.rho = 25;
    const Chain c = copying_chain(panel, cm, map);
    // The site folded into the symbol: x_j at site j is the symbol 2j + x_j, emitted by state k
    // with the copying probability of site j; a missing allele is a code beyond the symbols.
    genomics::hmm::Model hm(K, 2 * N);
    for (std::size_t k = 0; k < K; k++) {
      hm.initial[k] = 1.0 / static_cast<double>(K);
      for (std::size_t l = 0; l < K; l++) hm.a(k, l) = c.transition(1, k, l);
      for (std::size_t j = 0; j < N; j++)
        for (int a = 0; a < 2; a++) hm.e(k, 2 * j + static_cast<std::size_t>(a)) = panel(k, j) == a ? 1 - c.epsilon : c.epsilon;
    }
    genomics::hmm::Codes codes(N);
    for (std::size_t j = 0; j < N; j++) codes[j] = x[j] < 0 ? 255 : static_cast<genomics::hmm::Code>(2 * j + static_cast<std::size_t>(x[j]));
    const Posterior p = posterior(panel, c, x);
    const genomics::hmm::Posterior q = genomics::hmm::posterior(hm, codes);
    const Path v = viterbi(panel, c, x);
    const genomics::hmm::Path w = genomics::hmm::viterbi(hm, codes);
    const std::string n = " (K = 30, N = " + std::to_string(N) + ", a tenth missing)";
    s.check("on an evenly spaced map, log P(x) = genomics::hmm's forward" + n, calx::verification::relative(p.log_likelihood, q.log_likelihood), 1e-12, "relative");
    s.check("the posterior = genomics::hmm's", calx::verification::max_difference(p.gamma, q.gamma), 1e-10);
    s.check("Viterbi's log P(x, Z*) = genomics::hmm's", calx::verification::relative(v.log_probability, w.log_probability), 1e-12, "relative");
  }

  // ---- calibration on targets drawn from the model ----------------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 50;
    spec.theta = 30;
    spec.rho = 30;
    const Haplotypes panel = evolution::population::haplotypes(coalescent(spec, 4));
    CopyingModel cm;
    cm.rho = 30;
    cm.epsilon = 0.01;
    const Chain c = copying_chain(panel, cm);
    double n = 0, sum = 0, sum2 = 0, right = 0, sites = 0;
    for (std::uint32_t rep = 0; rep < 400; rep++) {
      calx::random::Stream r(4, MOSAIC, rep, 0);
      Mosaic m = sample_mosaic(panel, c, r);
      for (std::size_t j = 0; j < m.alleles.size(); j++)
        if (r.bernoulli(0.3)) m.alleles[j] = -1;
      const Posterior p = posterior(panel, c, m.alleles);
      const Path v = viterbi(panel, c, m.alleles);
      double d = 0;
      for (std::size_t j = 0; j < p.sites; j++) {
        double sq = 0;
        for (const double g : p.site(j)) sq += g * g;
        d += p.at(j, m.path[j]) - sq;
        right += panel(v.states[j], j) == panel(m.path[j], j);
        sites++;
      }
      d /= static_cast<double>(p.sites);
      n++;
      sum += d;
      sum2 += d * d;
    }
    const double mean = sum / n, se = std::sqrt((sum2 - sum * sum / n) / (n - 1) / n);
    s.check("E[γ_j(Z_j)] = E[Σ_k γ_j(k)²] on 400 mosaics drawn from the model, 30% of alleles hidden (z-score)", mean / se, 4, "SE");
    s.require("the Viterbi path copies the true path's allele at most sites (hidden ones included)", right / sites > 0.95, fmt(right / sites));
  }

  // ---- speed -----------------------------------------------------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 1001;
    spec.theta = 600;
    spec.rho = 600;
    const Haplotypes all = evolution::population::haplotypes(coalescent(spec, 5));
    Haplotypes panel;
    panel.sequences = 1000;
    panel.sites = all.sites;
    panel.positions = all.positions;
    std::vector<std::int8_t> x(all.sites);
    for (std::size_t j = 0; j < all.sites; j++) {
      const auto row = all.site(j);
      panel.alleles.insert(panel.alleles.end(), row.begin(), row.begin() + 1000);
      x[j] = j % 5 ? -1 : row[1000];
    }
    const Chain c = copying_chain(panel, {600, -1, -1});
    const auto t0 = std::chrono::steady_clock::now();
    const Posterior p = posterior(panel, c, x);
    const Path v = viterbi(panel, c, x);
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    s.require("forward–backward and Viterbi over a panel of 1000 haplotypes in under 2 s (O(N K))", ms < 2000 && std::isfinite(p.log_likelihood) && !v.states.empty(),
              fmt(ms) + " ms for " + std::to_string(all.sites) + " sites");
  }
  return s.done();
}
