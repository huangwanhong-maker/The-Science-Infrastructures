// Checks of chains with per-step emissions and transitions (chain.hpp): on a chain whose
// transitions change at every step, Viterbi, the forward and backward likelihoods and the
// posteriors equal the enumeration of every path; the chain read from a symbol model and a
// sequence gives model.hpp's Viterbi path, likelihood and posteriors; scaling a step's emissions
// by a constant leaves the posteriors and the path unchanged and shifts the likelihood by its log.

#include <genomics/hmm/chain.hpp>
#include <genomics/hmm/model.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

using namespace genomics;
using namespace genomics::hmm;

namespace {

// A random distribution of n values, as logs.
std::vector<double> random_log_distribution(calx::random::Stream& r, std::size_t n) {
  std::vector<double> p(n);
  double s = 0;
  for (auto& x : p) s += x = 0.05 + r.uniform();
  for (auto& x : p) x = std::log(x / s);
  return p;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(7, 0x43484b31, 0, 0);

  // ---- a chain of 3 states and 7 steps, every transition different ------------------------------
  {
    const std::size_t S = 3, T = 7;
    Chain c(S, T);
    c.log_initial = random_log_distribution(r, S);
    for (std::size_t t = 0; t < T; t++)
      for (std::size_t k = 0; k < S; k++) c.e(t, k) = std::log(0.01 + r.uniform());  // likelihoods, not distributions
    for (std::size_t t = 0; t + 1 < T; t++)
      for (std::size_t k = 0; k < S; k++) {
        const auto row = random_log_distribution(r, S);
        for (std::size_t l = 0; l < S; l++) c.a(t, k, l) = row[l];
      }
    // Enumeration of all 3⁷ paths.
    double total = NEG_INF, best = NEG_INF;
    std::vector<std::uint32_t> arg;
    std::vector<double> marginal(T * S, NEG_INF);
    std::vector<std::uint32_t> path(T, 0);
    for (std::size_t code = 0; code < 2187; code++) {
      std::size_t x = code;
      for (std::size_t t = 0; t < T; t++) { path[t] = static_cast<std::uint32_t>(x % S); x /= S; }
      const double lp = path_log_probability(c, path);
      total = log_add(total, lp);
      if (lp > best) { best = lp; arg = path; }
      for (std::size_t t = 0; t < T; t++) marginal[t * S + path[t]] = log_add(marginal[t * S + path[t]], lp);
    }
    const Path v = viterbi(c);
    double lf = 0, lb = 0;
    (void)log_forward(c, &lf);
    (void)log_backward(c, &lb);
    const Posterior p = posterior(c);
    double worst = 0;
    for (std::size_t i = 0; i < T * S; i++) worst = std::max(worst, std::fabs(p.gamma[i] - std::exp(marginal[i] - total)));
    s.require("Viterbi finds the best of the 2187 paths of a chain with a transition per step", v.states == arg && std::fabs(v.log_probability - best) < 1e-12);
    s.check("the forward likelihood is the sum over every path", std::fabs(lf - total), 1e-12);
    s.check("and so is the backward one", std::fabs(lb - total), 1e-12);
    s.check("the posteriors are the enumerated marginals", worst, 1e-12);
    // Scaling one step's emissions by a constant.
    Chain d = c;
    for (std::size_t k = 0; k < S; k++) d.e(3, k) += 2.5;
    const Posterior q = posterior(d);
    double moved = 0;
    for (std::size_t i = 0; i < T * S; i++) moved = std::max(moved, std::fabs(q.gamma[i] - p.gamma[i]));
    s.check("scaling a step's emissions leaves the posteriors unchanged", moved, 1e-12);
    s.check("and shifts the likelihood by its logarithm", std::fabs(q.log_likelihood - p.log_likelihood - 2.5), 1e-12);
    s.require("and leaves the Viterbi path", viterbi(d).states == v.states);
  }

  // ---- the chain of a symbol model equals the model ---------------------------------------------
  {
    Model m(4, 4);
    m.initial = {0.4, 0.3, 0.2, 0.1};
    for (std::size_t k = 0; k < 4; k++) {
      double rs = 0, es = 0;
      for (std::size_t l = 0; l < 4; l++) { rs += m.a(k, l) = 0.1 + r.uniform(); es += m.e(k, l) = 0.1 + r.uniform(); }
      for (std::size_t l = 0; l < 4; l++) { m.a(k, l) /= rs; m.e(k, l) /= es; }
    }
    Codes x(3000);
    for (auto& b : x) b = static_cast<Code>(r.below(4));
    x[100] = sequence::N;  // a missing symbol
    const Chain c = chain_of(m, x);
    const Path pm = viterbi(m, x), pc = viterbi(c);
    const Posterior qm = posterior(m, x), qc = posterior(c);
    double worst = 0;
    for (std::size_t i = 0; i < qm.gamma.size(); i++) worst = std::max(worst, std::fabs(qm.gamma[i] - qc.gamma[i]));
    s.require("a model's chain has the model's Viterbi path on 3000 symbols", pm.states == pc.states);
    s.check("its likelihood (log semiring) equals the scaled forward's", std::fabs(qc.log_likelihood - log_likelihood(m, x)) / std::fabs(log_likelihood(m, x)), 1e-12);
    s.check("and its posteriors", worst, 1e-9);
  }
  return s.done();
}
