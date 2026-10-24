// Checks of the hidden Markov model's three problems against their definitions: on tiny
// models every state path is enumerated, and the Viterbi path and score, the forward and
// backward likelihoods and the posterior of every state at every position are held to the
// maximum, the sum and the marginals over the enumeration (with and without an end state, with
// structural zeros, with a missing symbol); the scaled recursions agree with the log semiring on
// a sequence whose probability underflows any floating-point number; Baum–Welch never decreases
// the likelihood (without pseudocounts) or the log posterior (with them), reaches a likelihood
// at least that of the true parameters, and recovers them within five standard errors of
// their sampling, the same on one thread and on four; samples have the frequencies of their
// model.

#include <genomics/hmm/model.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::hmm;

namespace {

// A random distribution of n entries, some forced to zero by `zero_mask` bits.
void random_row(calx::random::Stream& r, double* p, std::size_t n, unsigned zero_mask = 0, double extra_weight = 0, double* extra = nullptr) {
  double s = 0;
  for (std::size_t i = 0; i < n; i++) s += p[i] = (zero_mask >> i) & 1 ? 0 : 0.1 + r.uniform();
  double e = 0;
  if (extra) s += e = extra_weight * (0.1 + r.uniform());
  for (std::size_t i = 0; i < n; i++) p[i] /= s;
  if (extra) *extra = e / s;
}

Model random_model(calx::random::Stream& r, std::size_t S, std::size_t K, bool with_end, bool zeros) {
  Model m(S, K, with_end);
  random_row(r, m.initial.data(), S);
  for (std::size_t k = 0; k < S; k++) {
    const unsigned mask = zeros ? (1u << ((k + 1) % S)) : 0;  // a structural zero in each row
    random_row(r, &m.transition[k * S], S, mask, 0.5, with_end ? &m.end[k] : nullptr);
    random_row(r, &m.emission[k * K], K, zeros && k == 0 ? 1u : 0u);
  }
  return m;
}

struct Enumeration {
  double best = NEG_INF, total = NEG_INF;
  std::vector<std::uint32_t> best_path;
  std::vector<double> marginal;  // T × S, the posterior by enumeration
};

Enumeration enumerate(const Model& m, const Codes& x) {
  const std::size_t S = m.states, T = x.size();
  Enumeration e;
  e.marginal.assign(T * S, 0);
  std::vector<std::uint32_t> path(T, 0);
  std::vector<double> weights;
  std::vector<std::vector<std::uint32_t>> paths;
  for (;;) {
    const double w = path_log_probability(m, x, path);
    if (w > e.best) { e.best = w; e.best_path = path; }
    e.total = log_add(e.total, w);
    weights.push_back(w);
    paths.push_back(path);
    std::size_t t = 0;
    while (t < T && ++path[t] == S) path[t++] = 0;
    if (t == T) break;
  }
  for (std::size_t p = 0; p < paths.size(); p++)
    if (weights[p] > NEG_INF)
      for (std::size_t t = 0; t < T; t++) e.marginal[t * S + paths[p][t]] += std::exp(weights[p] - e.total);
  return e;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x63686b31, 0, 0);

  // ---- tiny models against enumeration of every path ------------------------------------------------
  {
    double worst_v = 0, worst_f = 0, worst_b = 0, worst_g = 0, worst_semiring = 0;
    bool paths_ok = true, decode_ok = true;
    int cases = 0;
    for (int variant = 0; variant < 8; variant++) {
      const bool with_end = variant & 1, zeros = variant & 2, missing = variant & 4;
      const std::size_t S = 3, K = 3, T = 7;
      const Model m = random_model(r, S, K, with_end, zeros);
      m.validate();
      for (int rep = 0; rep < 4; rep++) {
        Codes x(T);
        for (auto& c : x) c = static_cast<Code>(r.below(K));
        if (missing) x[3] = 7;  // a code outside the alphabet: missing
        const Enumeration e = enumerate(m, x);
        const Path v = viterbi(m, x);
        const Posterior p = posterior(m, x);
        const Backward b = backward(m, x);
        cases++;
        worst_v = std::max(worst_v, std::fabs(v.log_probability - e.best));
        worst_f = std::max(worst_f, std::fabs(p.log_likelihood - e.total) / std::fabs(e.total));
        worst_b = std::max(worst_b, std::fabs(b.log_likelihood - e.total) / std::fabs(e.total));
        worst_semiring = std::max({worst_semiring, std::fabs(path_sum<Tropical>(m, x) - e.best), std::fabs(path_sum<LogSum>(m, x) - e.total),
                                   std::fabs(log_likelihood(m, x) - e.total)});
        paths_ok = paths_ok && std::fabs(path_log_probability(m, x, v.states) - v.log_probability) < 1e-12;
        for (std::size_t i = 0; i < T * S; i++) worst_g = std::max(worst_g, std::fabs(p.gamma[i] - e.marginal[i]));
        const auto d = posterior_decode(p);
        for (std::size_t t = 0; t < T; t++) {
          const double* row = &e.marginal[t * S];
          decode_ok = decode_ok && row[d[t]] >= *std::max_element(row, row + S) - 1e-12;
        }
      }
    }
    s.check("Viterbi: the best log probability equals the maximum over every path (3 states, 7 symbols, 32 cases)", worst_v, 1e-12);
    s.require("Viterbi: the path returned has the log probability claimed", paths_ok);
    s.check("forward: log P(x) equals the log of the sum over every path (relative)", worst_f, 1e-12);
    s.check("backward: log P(x) computed from the end equals the same sum (relative)", worst_b, 1e-12);
    s.check("the generic recurrence in the tropical and log semirings, and the O(S)-memory forward, equal the enumeration", worst_semiring, 1e-11);
    s.check("posterior: γ_t(k) equals the marginal over every path, at every position and state", worst_g, 1e-12);
    s.require("posterior decoding takes a state of largest marginal at every position", decode_ok);
  }

  // ---- a long sequence: scaling against the log semiring ----------------------------------------------
  {
    const Model m = random_model(r, 4, 4, true, false);
    Codes x(20000);
    for (auto& c : x) c = static_cast<Code>(r.below(4));
    const double f = log_likelihood(m, x), g = path_sum<LogSum>(m, x), b = backward(m, x).log_likelihood;
    s.check("on 20 000 symbols (P(x) ≈ e^" + std::to_string(static_cast<long>(f)) + ", far below any double), scaled forward equals the log semiring (relative)",
            std::max(std::fabs(f - g), std::fabs(b - g)) / std::fabs(g), 1e-11);
    const Path v = viterbi(m, x);
    s.check("on 20 000 symbols, Viterbi equals the tropical recurrence", std::fabs(v.log_probability - path_sum<Tropical>(m, x)), 1e-8);
    s.require("an impossible sequence has log P = −∞ and no Viterbi path", [&] {
      Model z = m;
      for (std::size_t k = 0; k < 4; k++) {
        z.e(k, 0) = 0;
        double q = z.e(k, 1) + z.e(k, 2) + z.e(k, 3);
        for (std::size_t b2 = 1; b2 < 4; b2++) z.e(k, b2) /= q;
      }
      const Codes y = {1, 0, 2};
      return log_likelihood(z, y) == NEG_INF && viterbi(z, y).states.empty() && backward(z, y).log_likelihood == NEG_INF;
    }());
  }

  // ---- Baum–Welch -------------------------------------------------------------------------------------------
  {
    // A three-state model with distinct emissions (the truth).
    Model truth(3, 4);
    truth.initial = {0.5, 0.3, 0.2};
    truth.transition = {0.90, 0.07, 0.03,  //
                        0.05, 0.90, 0.05,  //
                        0.04, 0.06, 0.90};
    truth.emission = {0.55, 0.15, 0.15, 0.15,  //
                      0.10, 0.40, 0.40, 0.10,  //
                      0.15, 0.10, 0.15, 0.60};
    truth.validate();
    std::vector<Codes> data;
    std::vector<std::vector<std::uint32_t>> paths;
    for (std::uint32_t i = 0; i < 40; i++) {
      calx::random::Stream rs(7, SAMPLE, i, 0);
      Sample smp = sample(truth, 1500, rs);
      data.push_back(smp.symbols);
      paths.push_back(smp.states);
    }
    // A start between the truth and uniform.
    Model start = truth;
    for (auto& p : start.transition) p = 0.5 * p + 0.5 / 3;
    for (auto& p : start.emission) p = 0.5 * p + 0.5 / 4;
    for (auto& p : start.initial) p = 1.0 / 3;
    start.normalise();

    Model plain = start;
    TrainOptions o;
    o.max_iterations = 500;
    o.tolerance = 1e-10;
    const Training t = baum_welch(plain, data, o);
    double worst_drop = 0;
    for (std::size_t i = 1; i < t.log_likelihood.size(); i++) worst_drop = std::max(worst_drop, t.log_likelihood[i - 1] - t.log_likelihood[i]);
    s.check("Baum–Welch without pseudocounts: the likelihood never decreases (" + std::to_string(t.iterations) + " iterations)", std::max(0.0, worst_drop), 1e-7, "nats");
    s.require("Baum–Welch converges", t.converged, std::to_string(t.iterations) + " iterations");

    double at_truth = 0;
    for (const auto& x : data) at_truth += log_likelihood(truth, x);
    s.require("Baum–Welch reaches a likelihood at least that of the true parameters",
              t.log_likelihood.back() >= at_truth - 1e-6, std::to_string(t.log_likelihood.back()) + " against " + std::to_string(at_truth));

    // Standard errors from the expected visits: √(p(1 − p)/n) for each parameter of a state
    // visited n times (the information of a multinomial; the hidden states only add to it).
    std::vector<double> visits(3, 0);
    for (const auto& x : data) {
      const Posterior p = posterior(plain, x);
      for (std::size_t tt = 0; tt + 1 < x.size(); tt++)
        for (std::size_t k = 0; k < 3; k++) visits[k] += p.at(tt, k);
    }
    double worst_z = 0;
    for (std::size_t k = 0; k < 3; k++) {
      for (std::size_t l = 0; l < 3; l++) {
        const double p = truth.a(k, l), se = std::sqrt(p * (1 - p) / visits[k]);
        worst_z = std::max(worst_z, std::fabs(plain.a(k, l) - p) / se);
      }
      for (std::size_t b = 0; b < 4; b++) {
        const double p = truth.e(k, b), se = std::sqrt(p * (1 - p) / visits[k]);
        worst_z = std::max(worst_z, std::fabs(plain.e(k, b) - p) / se);
      }
    }
    s.check("Baum–Welch recovers the transitions and emissions of a simulated model (60 000 symbols)", worst_z, 5, "standard errors");

    Model with_prior = start;
    TrainOptions op = o;
    op.pseudocounts = {1, 1, 1, 0};
    const Training tp = baum_welch(with_prior, data, op);
    double worst_post = 0, worst_lik = 0;
    for (std::size_t i = 1; i < tp.log_posterior.size(); i++) {
      worst_post = std::max(worst_post, tp.log_posterior[i - 1] - tp.log_posterior[i]);
      worst_lik = std::max(worst_lik, tp.log_likelihood[i - 1] - tp.log_likelihood[i]);
    }
    s.check("Baum–Welch with pseudocounts: the log posterior never decreases", std::max(0.0, worst_post), 1e-7, "nats");
    s.require("with pseudocounts the parameters stay positive and the model valid", [&] {
      with_prior.validate(1e-9);
      for (double p : with_prior.emission)
        if (!(p > 0)) return false;
      return true;
    }());

    // Structural zeros are kept.
    Model left_right(3, 4);
    left_right.initial = {1, 0, 0};
    left_right.transition = {0.8, 0.2, 0, 0, 0.7, 0.3, 0, 0, 1};
    left_right.emission = truth.emission;
    TrainOptions ol;
    ol.max_iterations = 20;
    ol.pseudocounts = {0, 1, 1, 0};
    baum_welch(left_right, data, ol);
    s.require("Baum–Welch keeps a model's zeros (a left-to-right model stays one)", left_right.a(1, 0) == 0 && left_right.a(2, 0) == 0 && left_right.a(2, 1) == 0 && left_right.a(0, 2) == 0 && left_right.initial[1] == 0);

    // The same on four threads as on one.
    Model threaded = start;
    calx::scheduler::ThreadPool pool(4);
    TrainOptions ot = o;
    ot.max_iterations = 30;
    Model serial = start;
    const Training a = baum_welch(serial, data, ot);
    const Training b = baum_welch(threaded, data, ot, &pool);
    s.require("Baum–Welch on four threads gives the serial result bit for bit", a.log_likelihood == b.log_likelihood && serial.transition == threaded.transition && serial.emission == threaded.emission);

    // With an end state: lengths and the end probabilities re-estimated.
    Model ending(2, 4, true);
    ending.initial = {0.6, 0.4};
    ending.transition = {0.95, 0.04, 0.03, 0.96};
    ending.end = {0.01, 0.01};
    ending.emission = {0.4, 0.1, 0.1, 0.4, 0.1, 0.4, 0.4, 0.1};
    ending.validate();
    std::vector<Codes> runs;
    double length = 0;
    for (std::uint32_t i = 0; i < 400; i++) {
      calx::random::Stream rs(8, SAMPLE, i, 1);
      runs.push_back(sample_to_end(ending, rs, 100000).symbols);
      length += static_cast<double>(runs.back().size());
    }
    length /= 400;
    s.check("sampling to the end: the mean length is 1/a_k0 = 100 (400 samples, se 5)", (length - 100) / 5, 5, "standard errors");
    Model guess = ending;
    guess.end = {0.05, 0.05};
    guess.normalise();
    guess.transition = {0.9, 0.05, 0.05, 0.9};
    guess.normalise();
    TrainOptions oe;
    oe.max_iterations = 300;
    const Training te = baum_welch(guess, runs, oe);
    double drop = 0;
    for (std::size_t i = 1; i < te.log_likelihood.size(); i++) drop = std::max(drop, te.log_likelihood[i - 1] - te.log_likelihood[i]);
    s.check("Baum–Welch with an end state: the likelihood never decreases", std::max(0.0, drop), 1e-7, "nats");
    s.check("Baum–Welch recovers the end probability (0.01) from 400 sequences", std::fabs(0.5 * (guess.end[0] + guess.end[1]) - 0.01), 0.002);
  }

  // ---- sampling ----------------------------------------------------------------------------------------------
  {
    const Model m = random_model(r, 3, 4, false, false);
    calx::random::Stream rs(9, SAMPLE, 0, 0);
    const Sample smp = sample(m, 300000, rs);
    std::vector<double> n(9, 0), from(3, 0), e(12, 0);
    for (std::size_t t = 0; t < smp.states.size(); t++) {
      e[smp.states[t] * 4 + smp.symbols[t]] += 1;
      if (t + 1 < smp.states.size()) { n[smp.states[t] * 3 + smp.states[t + 1]] += 1; from[smp.states[t]] += 1; }
    }
    double worst = 0;
    for (std::size_t k = 0; k < 3; k++) {
      double visits = 0;
      for (std::size_t b = 0; b < 4; b++) visits += e[k * 4 + b];
      for (std::size_t l = 0; l < 3; l++) {
        const double p = m.a(k, l);
        worst = std::max(worst, std::fabs(n[k * 3 + l] / from[k] - p) / std::sqrt(p * (1 - p) / from[k]));
      }
      for (std::size_t b = 0; b < 4; b++) {
        const double p = m.e(k, b);
        worst = std::max(worst, std::fabs(e[k * 4 + b] / visits - p) / std::sqrt(p * (1 - p) / visits));
      }
    }
    s.check("a sampled path of 300 000 steps has the transition and emission frequencies of its model", worst, 5, "standard errors");
    calx::random::Stream again(9, SAMPLE, 0, 0);
    s.require("a sample is a function of its stream", sample(m, 1000, again).symbols == Codes(smp.symbols.begin(), smp.symbols.begin() + 1000));
  }

  // ---- arguments --------------------------------------------------------------------------------------------
  {
    bool threw = false;
    try {
      Model bad(2, 2);
      bad.initial = {0.5, 0.6};
      bad.validate();
    } catch (const std::invalid_argument&) { threw = true; }
    s.require("validate() rejects initial probabilities that do not sum to one", threw);
  }
  return s.done();
}
