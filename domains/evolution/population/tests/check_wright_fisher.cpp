// Checks of Wright–Fisher populations: without selection an allele fixes with probability equal
// to its frequency, and the mean heterozygosity decays exactly as (1 − 1/(2N))^t; with genic
// selection the fixation probability is Kimura's (1 − e^{−4Nsp})/(1 − e^{−4Ns}) (1962) for
// advantageous and deleterious alleles, and with dominance the diffusion's integral; among
// islands exchanging migrants F_ST settles near 1/(1 + 4Nm) (Wright 1931); mutation holds the
// frequency at its balance ν/(μ + ν); a trajectory is a function of its seed; the replicates'
// summary is the same on 1 and 8 threads. Each expectation is checked within four standard
// errors of the simulation (the z-score is reported).

#include <evolution/population/wright_fisher.hpp>

#include <calx/verification/check.hpp>

#include <chrono>
#include <cmath>
#include <string>
#include <vector>

using namespace evolution::population;

namespace {
std::string fmt(double x) { return calx::verification::exponential(x); }
}  // namespace

int main() {
  calx::verification::Suite s;
  calx::scheduler::ThreadPool pool(8);

  // ---- fixation ----------------------------------------------------------------------------------------
  const double N = 100, p0 = 0.1;
  const std::size_t R = 20000;
  const auto z_of = [&](double observed, double expected) { return (observed - expected) / std::sqrt(expected * (1 - expected) / R); };
  {
    WrightFisherSpec w;
    w.individuals = 100;
    w.p0 = p0;
    w.generations = 100000;
    const WrightFisherSummary neutral = wright_fisher_replicates(w, 1, R, pool, 1000);
    s.require("every replicate fixes or loses the allele", neutral.fixed + neutral.lost == R);
    s.check("without selection P(fixation) = p = 0.1 (N = 100, 20000 replicates; z-score)", z_of(neutral.fixation(), p0), 4, "SE");
    s.check("and the mean time to absorption is −4N(p ln p + q ln q) generations (Kimura and Ohta 1969; relative error)",
            neutral.mean_settling / (-4 * N * (p0 * std::log(p0) + (1 - p0) * std::log(1 - p0))) - 1, 0.03);
    for (const double s_k : {0.0025, 0.005, 0.01, -0.0025}) {
      w.s = 2 * s_k;
      w.h = 0.5;
      const WrightFisherSummary sel = wright_fisher_replicates(w, 2, R, pool, 1000);
      const double u = kimura_fixation(N, s_k, p0);
      s.check("Kimura's u(p) for genic selection, 4Ns = " + fmt(4 * N * s_k) + " (expected " + fmt(u) + "; z-score)", z_of(sel.fixation(), u), 4, "SE");
    }
    w.s = 0.02;
    w.h = 0;
    const WrightFisherSummary recessive = wright_fisher_replicates(w, 3, R, pool, 1000);
    const double u = fixation_probability(N, 0.02, 0, p0);
    s.check("the diffusion's u(p) for a recessive advantage, s = 0.02, h = 0 (expected " + fmt(u) + "; z-score)", z_of(recessive.fixation(), u), 4, "SE");
    s.check("the general integral equals Kimura's closed form at h = ½", fixation_probability(N, 0.01, 0.5, p0) - kimura_fixation(N, 0.005, p0), 1e-10);
  }

  // ---- the decay of heterozygosity ----------------------------------------------------------------------
  {
    WrightFisherSpec w;
    w.individuals = 50;
    w.p0 = 0.5;
    w.generations = 200;
    const WrightFisherSummary sum = wright_fisher_replicates(w, 4, R, pool, 1);
    double worst = 0;
    for (const std::size_t t : {10u, 50u, 100u, 200u}) {
      const double expected = 0.5 * heterozygosity_decay(50, static_cast<double>(t));
      const double var = sum.h_s_sq[t] - sum.h_s[t] * sum.h_s[t];
      const double z = (sum.h_s[t] - expected) / std::sqrt(var / R);
      if (std::fabs(z) > std::fabs(worst)) worst = z;
    }
    s.check("E[2p_t q_t] = 2p₀q₀ (1 − 1/(2N))^t at t = 10, 50, 100, 200 (N = 50; the largest z-score)", worst, 4, "SE");
  }

  // ---- islands -------------------------------------------------------------------------------------------
  {
    WrightFisherSpec w;
    w.individuals = 50;
    w.demes = 100;
    w.m = 0.005;
    w.mu = w.nu = 1e-5;
    w.p0 = 0.5;
    w.generations = 6000;
    const WrightFisherSummary sum = wright_fisher_replicates(w, 5, 16, pool, 1, 2);
    const double fst = sum.fst(1000, 6001);
    s.check("among 100 islands with 4Nm = 1, F_ST = 1 − H_S/H_T settles near 1/(1 + 4Nm) = 0.5 (exact for infinitely many islands " +
                fmt(island_fst_exact(50, 0.005)) + "; absolute error)",
            fst - island_fst(50, 0.005), 0.03);
    w.m = 0.05;
    const double fst_high = wright_fisher_replicates(w, 6, 16, pool, 1, 2).fst(1000, 6001);
    s.check("and near 1/(1 + 4Nm) = 1/11 with 4Nm = 10 (absolute error)", fst_high - island_fst(50, 0.05), 0.02);
  }

  // ---- mutation balance --------------------------------------------------------------------------------
  {
    WrightFisherSpec w;
    w.individuals = 500;
    w.mu = 0.002;
    w.nu = 0.001;
    w.p0 = 0.9;
    w.generations = 4000;
    const WrightFisherSummary sum = wright_fisher_replicates(w, 7, 2000, pool, 100);
    double mean = 0;
    for (std::size_t k = 20; k <= 40; k++) mean += sum.mean_p[k];
    mean /= 21;
    // E[p_t] relaxes to ν/(μ + ν) as (1 − μ − ν)^t exactly; the standard error of the mean over
    // 2000 replicates and 21 generations 100 apart is about 0.0013.
    s.check("mutation holds the mean frequency at ν/(μ + ν) = 1/3 (absolute error, 2000 replicates)", mean - 1.0 / 3, 0.005);
  }

  // ---- determinism -------------------------------------------------------------------------------------
  {
    WrightFisherSpec w;
    w.individuals = 200;
    w.demes = 5;
    w.s = 0.01;
    w.m = 0.01;
    w.mu = w.nu = 1e-4;
    w.generations = 300;
    const Trajectory a = wright_fisher_trajectory(w, 8, 2), b = wright_fisher_trajectory(w, 8, 2), c = wright_fisher_trajectory(w, 8, 3);
    bool bounded = a.generations() == 300;
    for (const double x : a.p) bounded = bounded && x >= 0 && x <= 1;
    s.require("a trajectory is a function of its seed and replicate, frequencies in [0, 1]", a.p == b.p && a.p != c.p && bounded);
    calx::scheduler::ThreadPool one(1);
    const auto t0 = std::chrono::steady_clock::now();
    const WrightFisherSummary x = wright_fisher_replicates(w, 9, 1000, one, 10);
    const double serial = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    const auto t1 = std::chrono::steady_clock::now();
    const WrightFisherSummary y = wright_fisher_replicates(w, 9, 1000, pool, 10);
    const double parallel = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t1).count();
    s.require("the replicates' summary is the same on 1 and 8 threads", x.mean_p == y.mean_p && x.h_s == y.h_s && x.h_t == y.h_t && x.fixed == y.fixed,
              "1000 replicates of 300 generations × 5 demes: " + fmt(serial) + " ms on 1 thread, " + fmt(parallel) + " ms on 8");
    s.require("the task graph ran one task per shard and the merge", y.trace.records.size() == (1000 + 63) / 64 + 1);
  }
  return s.done();
}
