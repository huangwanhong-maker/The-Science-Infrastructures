// Checks of training: the tied M-step is the exact maximiser of the expected complete
// log-likelihood on the tied parameters (its stationarity conditions hold, and no tied parameter
// vector, random or near it, does better); EM never decreases the likelihood, untied or tied, and
// keeps a tied model tied; its E-step is the same on one thread and on four; from structures
// known, the M-step is the complete-data maximum likelihood and recovers the truth; and the
// likelihood-ratio test of the reversal tie on sequences simulated from the grammar: its
// statistic averages near its degrees of freedom when the truth is reversal-symmetric, and
// rejects when it is not.

#include <genomics/rna/scfg.hpp>
#include <genomics/rna/symmetry.hpp>
#include <genomics/rna/training.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <thread>
#include <vector>

using namespace genomics;
using namespace genomics::rna;

namespace {

Codes random_rna(calx::random::Stream& r, std::size_t n) {
  Codes x(n);
  for (auto& c : x) c = static_cast<Code>(r.below(4));
  return x;
}

Parameters random_parameters(const StructureGrammar& sg, calx::random::Stream& r) {
  Parameters t(sg.rules());
  std::vector<double> sum(sg.nonterminals(), 0);
  for (std::size_t k = 0; k < t.size(); k++) sum[sg.info[k].lhs] += t[k] = 0.05 + r.uniform();
  for (std::size_t k = 0; k < t.size(); k++) t[k] /= sum[sg.info[k].lhs];
  return t;
}

// A model of the canonical grammar whose pairs on each left side are in the proportions
// AU : UA : CG : GC : GU : UG = w[0..5] (reversal-symmetric when w pairs them equally).
Parameters truth(const StructureGrammar& sg, const double (&w)[6]) {
  Parameters t = starting_parameters(sg, {0.75, 0.65, 1});
  const Code A = sequence::A, C = sequence::C, G = sequence::G, U = sequence::T;
  const Code pa[6] = {A, U, C, G, G, U}, pb[6] = {U, A, G, C, U, G};
  for (std::size_t side = 1; side < sg.nonterminals(); side++) {
    double mass = 0, total = 0;
    for (std::size_t r = 0; r < sg.rules(); r++)
      if (sg.info[r].lhs == side && sg.info[r].kind == RuleInfo::pair) mass += t[r];
    for (double v : w) total += v;
    for (int k = 0; k < 6; k++) {
      const std::size_t r = side == 1 ? sg.L_pair[pa[k]][pb[k]] : side == 2 ? sg.F_pair[pa[k]][pb[k]] : sg.H_pair[pa[k]][pb[k]];
      t[r] = mass * w[k] / total;
    }
  }
  validate(sg, t);
  return t;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream rng(2026, 103, 0, 0);
  calx::scheduler::ThreadPool pool(std::max(2u, std::thread::hardware_concurrency()));

  // ---- the tied M-step maximises Q on the tied parameters ----------------------------------------------------
  {
    const auto& khr = knudsen_hein_reversible();
    std::vector<Codes> data;
    for (int k = 0; k < 8; k++) data.push_back(random_rna(rng, 20 + rng.below(20)));
    const Parameters at = random_parameters(khr, rng);
    const auto e = expectation(khr, at, data, &pool);
    bool stationary = true, best_random = true, best_near = true, below_untied = true, tied_valid = true;
    std::string detail;
    for (const auto& g : pair_symmetries()) {
      const Tie t = orbit_tie(khr, g);
      const Parameters star = maximisation(khr, t, e.counts, at);
      validate(khr, star);
      tied_valid = tied_valid && is_tied(t, star);
      const double q = expected_complete_log_likelihood(e.counts, star);
      // Stationarity: C_O / (|O| θ_O) equals the left side's total count for every orbit.
      std::vector<double> orbit(t.orbits(), 0), side(khr.nonterminals(), 0);
      for (std::size_t r = 0; r < khr.rules(); r++) { orbit[t.orbit_of[r]] += e.counts[r]; side[khr.info[r].lhs] += e.counts[r]; }
      for (std::size_t r = 0; r < khr.rules(); r++) {
        const std::size_t o = t.orbit_of[r];
        if (orbit[o] > 0) stationary = stationary && std::fabs(orbit[o] / (static_cast<double>(t.size[o]) * star[r]) - side[khr.info[r].lhs]) <= 1e-9 * side[khr.info[r].lhs];
      }
      double gap = 1e300;
      for (int k = 0; k < 500; k++) {
        const Parameters other = project(t, random_parameters(khr, rng));
        gap = std::min(gap, q - expected_complete_log_likelihood(e.counts, other));
        best_random = best_random && expected_complete_log_likelihood(e.counts, other) <= q;
        // Nearby: θ* + ε d, d tied and summing to zero on each left side.
        Parameters d = project(t, random_parameters(khr, rng));
        for (std::size_t r = 0; r < d.size(); r++) d[r] -= star[r];
        const double eps = 1e-3 * rng.uniform();
        Parameters near(star.size());
        for (std::size_t r = 0; r < d.size(); r++) near[r] = star[r] + eps * d[r];
        best_near = best_near && expected_complete_log_likelihood(e.counts, near) <= q + 1e-12 * std::fabs(q);
      }
      below_untied = below_untied && q <= expected_complete_log_likelihood(e.counts, maximisation(khr, untied(khr), e.counts, at)) + 1e-9;
      char line[96];
      std::snprintf(line, sizeof line, "%s%s %.3g", detail.empty() ? "least gap of Q to a random tied θ: " : "; ", g.name.c_str(), gap);
      detail += line;
    }
    s.require("the tied M-step gives a tied probability vector (ρ, κ, κρ, the strand group)", tied_valid);
    s.require("it satisfies the stationarity conditions C_O / (|O| θ_O) = Σ_A C", stationary);
    s.require("no random tied θ has a larger expected complete log-likelihood", best_random, detail);
    s.require("nor any tied θ near it", best_near);
    s.require("the untied maximum is at least the tied one", below_untied);
  }

  // ---- EM on simulated sequences ---------------------------------------------------------------------------------
  const auto& canon = structure_grammar(Shape::reversible, false, canonical_pairs());
  const auto rho = pair_symmetries()[0];
  const double symmetric[6] = {3, 3, 4, 4, 1, 1}, asymmetric[6] = {4, 1, 1, 6, 2, 0.5};
  {
    const Parameters theta = truth(canon, asymmetric);
    std::vector<Codes> data;
    for (int k = 0; k < 30; k++) data.push_back(sample(canon, theta, rng, 25, 60).sequence);
    const auto one = expectation(canon, theta, data, nullptr);
    calx::scheduler::ThreadPool four(4);
    const auto many = expectation(canon, theta, data, &four);
    s.require("the E-step is the same on one thread and on four, bit for bit", one.counts == many.counts && one.log_likelihood == many.log_likelihood);
    TrainOptions o;
    o.max_iterations = 30;
    const auto free = train(canon, untied(canon), starting_parameters(canon), data, o, &pool);
    const Tie t = orbit_tie(canon, rho);
    const auto tied = train(canon, t, starting_parameters(canon), data, o, &pool);
    char line[200];
    std::snprintf(line, sizeof line, "untied %.2f → %.2f in %d iterations; ρ-tied %.2f → %.2f in %d", free.log_likelihood.front(), free.log_likelihood.back(),
                  free.iterations, tied.log_likelihood.front(), tied.log_likelihood.back(), tied.iterations);
    s.require("EM never decreases the likelihood, untied or tied along ρ's orbits", free.monotone && tied.monotone, line);
    s.require("every tied iterate stays tied", is_tied(t, tied.theta, 1e-15));
  }

  // ---- complete data: the M-step is the maximum likelihood of known structures -------------------------------------
  {
    const Parameters theta = truth(canon, asymmetric);
    std::vector<double> counts(canon.rules(), 0);
    double at_truth = 0;
    std::vector<std::pair<Codes, Structure>> known;
    for (int k = 0; k < 2000; k++) {
      const auto smp = sample(canon, theta, rng, 1, 400);
      const auto c = rule_counts(canon, smp.sequence, smp.structure);
      for (std::size_t r = 0; r < c.size(); r++) counts[r] += c[r];
      at_truth += log_joint(canon, theta, smp.sequence, smp.structure);
      known.push_back({smp.sequence, smp.structure});
    }
    const Parameters hat = maximisation(canon, untied(canon), counts, theta);
    double at_hat = 0, worst = 0;
    for (const auto& [x, y] : known) at_hat += log_joint(canon, hat, x, y);
    for (std::size_t r = 0; r < hat.size(); r++) worst = std::max(worst, std::fabs(hat[r] - theta[r]));
    s.require("from known structures the M-step is the complete-data maximum likelihood: at least the truth's", at_hat >= at_truth,
              std::to_string(at_hat) + " against " + std::to_string(at_truth));
    s.check("and recovers the truth from 2000 sampled structures (largest difference of a rule probability)", worst, 0.03);
  }

  // ---- the likelihood-ratio test of the reversal tie ------------------------------------------------------------------
  {
    const Tie t = orbit_tie(canon, rho);
    const std::size_t df = free_parameters(untied(canon)) - free_parameters(t);
    TrainOptions o;
    o.max_iterations = 60;
    const auto fit = [&](const Parameters& theta, double* tied_l, double* free_l) {
      std::vector<Codes> data;
      for (int k = 0; k < 40; k++) data.push_back(sample(canon, theta, rng, 25, 60).sequence);
      const auto a = train(canon, t, theta, data, o, &pool);  // from the truth's projection
      const auto b = train(canon, untied(canon), a.theta, data, o, &pool);
      *tied_l = a.log_likelihood.back();
      *free_l = b.log_likelihood.back();
      return a.monotone && b.monotone;
    };
    double total = 0;
    int replicates = 6;
    bool monotone = true;
    std::string stats;
    for (int k = 0; k < replicates; k++) {
      double a = 0, b = 0;
      monotone = fit(truth(canon, symmetric), &a, &b) && monotone;
      const auto lr = likelihood_ratio(b, free_parameters(untied(canon)), a, free_parameters(t));
      total += lr.statistic;
      char line[48];
      std::snprintf(line, sizeof line, "%s%.1f", stats.empty() ? "" : ", ", lr.statistic);
      stats += line;
    }
    const double mean = total / replicates, se = std::sqrt(2.0 * static_cast<double>(df) / replicates);
    double a = 0, b = 0;
    monotone = fit(truth(canon, asymmetric), &a, &b) && monotone;
    const auto rejected = likelihood_ratio(b, free_parameters(untied(canon)), a, free_parameters(t));
    s.require("the reversal tie removes 9 parameters of the canonical grammar (three pair orbits of two on each of L, H, F)", df == 9, std::to_string(df));
    s.require("EM from the tied fit to the untied one never decreases the likelihood", monotone);
    s.check("with a reversal-symmetric truth Λ averages its degrees of freedom (standard scores of the mean of 6)", (mean - static_cast<double>(df)) / se, 3.0);
    char line[160];
    std::snprintf(line, sizeof line, "Λ = %.1f on %zu degrees of freedom, p = %.2g (symmetric truth: Λ = %s)", rejected.statistic, rejected.df, rejected.p_value, stats.c_str());
    s.require("with an asymmetric truth the test rejects the tie at p < 1e-6", rejected.p_value < 1e-6, line);
  }

  return s.done();
}
