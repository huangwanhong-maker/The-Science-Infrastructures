// Checks of the forward simulator of a selective sweep: the probability that a new beneficial
// allele fixes, against Kimura's diffusion u(1/(2N)) = (1 − e^{−2s})/(1 − e^{−4Ns}) for its
// haploid selection; one generation's law, Binomial(2N, ψ(p)) with ψ(p) = p(1 + s)/(1 + sp),
// by the residuals of the trajectories where conditioning on fixation does not reach; the
// duration of a sweep below the deterministic logistic's 2 ln(2N − 1)/ln(1 + s) generations;
// hitchhiking (Maynard Smith and Haigh 1974), the loss of diversity near the selected site,
// against what drift alone leaves of the same standing variation, (1 − 1/(2N))^G, which far
// from the site is what remains; a partial sweep stopped at its frequency, the selected site
// segregating in the sample at about that frequency; and replicates reproduced from their seeds.

#include <evolution/population/coalescent.hpp>
#include <evolution/population/sweep.hpp>
#include <evolution/population/wright_fisher.hpp>

#include <calx/verification/check.hpp>

#include <chrono>
#include <cmath>
#include <string>
#include <vector>

using namespace evolution::population;

namespace {

std::string fmt(double x) { return calx::verification::exponential(x); }

// π summed over the sites in [lo, hi) of the chromosome, per unit of length.
double pi_in(const Haplotypes& h, double lo, double hi) {
  double pi = 0;
  const double n = static_cast<double>(h.sequences);
  for (std::size_t j = 0; j < h.sites; j++) {
    if (h.positions[j] < lo || h.positions[j] >= hi) continue;
    double k = 0;
    for (const std::int8_t a : h.site(j)) k += a;
    pi += 2 * k * (n - k) / (n * (n - 1));
  }
  return pi / (hi - lo);
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- fixation, duration and one generation's law -------------------------------------------------------
  {
    SweepSpec spec;
    spec.individuals = 250;
    spec.samples = 20;
    spec.theta = 2;
    spec.rho = 10;
    spec.s = 0.05;
    const std::size_t R = 300;
    double attempts = 0, generations = 0, n = 0, sum = 0, sum2 = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (std::uint32_t rep = 0; rep < R; rep++) {
      const Sweep w = sweep(spec, 41, rep);
      attempts += static_cast<double>(w.attempts);
      generations += static_cast<double>(w.generations);
      // Where p ≥ 0.2 the allele is lost with probability about e^{−4Ns p} < 1e-4, so conditioning
      // on fixation leaves the chain's law there unchanged: the residuals
      // (p′ − ψ(p)) / √(ψ(1 − ψ)/(2N)) have mean 0 and variance 1.
      for (std::size_t t = 0; t + 1 < w.trajectory.size(); t++) {
        const double p = w.trajectory[t];
        if (p < 0.2 || p > 0.95) continue;
        const double psi = p * (1 + spec.s) / (1 + spec.s * p), r = (w.trajectory[t + 1] - psi) / std::sqrt(psi * (1 - psi) / 500);
        n++;
        sum += r;
        sum2 += r * r;
      }
    }
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    const double u = kimura_fixation(250, spec.s, 1.0 / 500), estimate = static_cast<double>(R) / attempts;
    // The count of attempts is geometric: Var[attempts/R] = (1 − u)/(u² R), so SE[û] ≈ u √((1 − u)/R).
    const double se = u * std::sqrt((1 - u) / static_cast<double>(R));
    s.check("P(a new beneficial allele fixes) = Kimura's u(1/(2N)) = " + fmt(u) + " (2N = 500, s = 0.05; " + std::to_string(R) + " sweeps, z-score)", (estimate - u) / se, 4, "SE");
    const double logistic = 2 * std::log(499.0) / std::log(1.05), mean = generations / static_cast<double>(R);
    s.require("a sweep that fixes takes fewer generations than the deterministic logistic's 2 ln(2N − 1)/ln(1 + s) from one copy (conditioning on fixation hastens its start)",
              mean < logistic && mean > 0.5 * logistic, fmt(mean) + " against " + fmt(logistic) + "; " + fmt(ms) + " ms in all");
    s.check("the frequency's residual about ψ(p) = p(1 + s)/(1 + sp) has mean 0 (" + fmt(n) + " generations with 0.2 ≤ p ≤ 0.95; z-score)", sum / n * std::sqrt(n), 4, "SE");
    s.check("and variance ψ(1 − ψ)/(2N), Wright–Fisher's binomial (the ratio's departure from 1)", sum2 / n - 1, 0.05);
  }

  // ---- hitchhiking --------------------------------------------------------------------------------------------
  {
    SweepSpec spec;
    spec.individuals = 500;
    spec.samples = 100;
    spec.theta = 60;
    spec.rho = 400;
    spec.s = 0.05;
    // π after the sweep against π of the same standing variation before it (the coalescent's
    // replicate of all 2N chromosomes, as sweep() draws it). Without new mutation, drift alone
    // multiplies the expected heterozygosity of a neutral site by (1 − 1/(2N))^G in G generations
    // (theory/population.md, Proposition 1.3). Near the selected site hitchhiking removes far more;
    // at distances over 0.4 (recombination fractions over 0.08 per generation, against s = 0.05)
    // what is left is drift's share, less the little the sweep still takes there.
    CoalescentSpec standing;
    standing.samples = 1000;
    standing.theta = 60;
    standing.rho = 400;
    double near_after = 0, near_before = 0, far_after = 0, far_before = 0, drift = 0;
    bool monomorphic = true;
    const std::size_t R = 40;
    for (std::uint32_t rep = 0; rep < R; rep++) {
      const Sweep w = sweep(spec, 42, rep);
      monomorphic = monomorphic && w.selected == Sweep::npos && w.frequency == 1;
      const Haplotypes before = haplotypes(coalescent(standing, 42, rep));
      const double d = std::pow(1 - 1.0 / 1000, static_cast<double>(w.generations));
      near_after += pi_in(w.haplotypes, 0.49, 0.51);
      near_before += d * pi_in(before, 0.49, 0.51);
      far_after += pi_in(w.haplotypes, 0, 0.1) + pi_in(w.haplotypes, 0.9, 1);
      far_before += d * (pi_in(before, 0, 0.1) + pi_in(before, 0.9, 1));
      drift += d / R;
    }
    s.require("a complete sweep fixes the selected allele, which then does not segregate in the sample", monomorphic);
    s.require("hitchhiking: π within 0.01 of the selected site falls below 30% of what drift alone leaves (2N = 1000, s = 0.05, θ = 60, ρ = 400; 40 sweeps)",
              near_after < 0.3 * near_before, fmt(near_after / near_before) + " of it; drift's factor (1 − 1/(2N))^G " + fmt(drift));
    s.check("far away (distances over 0.4) π is drift's share of the standing π within 10% (the ratio's departure from 1)", far_after / far_before - 1, 0.1);
  }

  // ---- a partial sweep -------------------------------------------------------------------------------------------
  {
    SweepSpec spec;
    spec.individuals = 500;
    spec.samples = 100;
    spec.theta = 40;
    spec.rho = 80;
    spec.s = 0.05;
    spec.stop = 0.7;
    double worst = 0;
    bool stopped = true, found = true;
    for (std::uint32_t rep = 0; rep < 20; rep++) {
      const Sweep w = sweep(spec, 44, rep);
      const std::size_t T = w.trajectory.size();
      stopped = stopped && w.frequency >= 0.7 && T >= 2 && w.trajectory[T - 2] < 0.7 && w.generations + 1 == T;
      found = found && w.selected != Sweep::npos && w.haplotypes.positions[w.selected] == spec.position;
      if (w.selected == Sweep::npos) continue;
      double k = 0;
      for (const std::int8_t a : w.haplotypes.site(w.selected)) k += a;
      const double f = k / 100, sd = std::sqrt(w.frequency * (1 - w.frequency) / 100);
      worst = std::max(worst, std::fabs(f - w.frequency) / sd);
    }
    s.require("a partial sweep stops at the first generation its allele reaches 0.7", stopped);
    s.require("and its selected site segregates in the sample, at its position", found);
    s.check("the sample frequency of the selected allele is the population's within sampling error (largest |z| of 20)", worst, 4, "SD");
  }

  // ---- reproducibility -------------------------------------------------------------------------------------------
  {
    SweepSpec spec;
    spec.individuals = 200;
    spec.samples = 50;
    spec.theta = 20;
    spec.rho = 40;
    spec.stop = 0.5;
    const Sweep a = sweep(spec, 45, 3), b = sweep(spec, 45, 3), c = sweep(spec, 45, 4);
    s.require("a replicate is reproduced from its seed, and another differs",
              a.haplotypes.alleles == b.haplotypes.alleles && a.haplotypes.positions == b.haplotypes.positions && a.attempts == b.attempts &&
                  a.trajectory == b.trajectory && a.haplotypes.alleles != c.haplotypes.alleles);
  }
  return s.done();
}
