// Checks of the haplotype scans: EHH, iHH and the unstandardised iHS of a small example worked
// by hand, with and without the cut-off; EHH by the partition walk against its definition by
// pairs of carriers, on the coalescent's haplotypes with missing alleles; EHH equal to 1 at the
// core and non-increasing away from it; a scan the same on 1 and 8 threads; iHS standardised to
// mean 0 and variance 1 within each bin of frequency, by construction; and the signal of
// selection: partial sweeps (population/sweep.hpp) against the neutral coalescent of the same
// θ and ρ, their selected sites' iHS strongly negative where neutral sites' |iHS| is rarely large.

#include <evolution/haplotypes/ehh.hpp>
#include <evolution/population/coalescent.hpp>
#include <evolution/population/sweep.hpp>

#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace evolution::haplotypes;
using evolution::population::CoalescentSpec;
using evolution::population::coalescent;
using evolution::population::Sweep;
using evolution::population::SweepSpec;

namespace {

std::string fmt(double x) { return calx::verification::exponential(x); }

// EHH by its definition: the pairs of carriers identical, and typed, at every site from the core to j.
double ehh_by_pairs(const Haplotypes& h, std::size_t core, std::int8_t allele, std::size_t j) {
  std::vector<std::size_t> c;
  for (std::size_t i = 0; i < h.sequences; i++)
    if (h(i, core) == allele) c.push_back(i);
  const std::size_t lo = std::min(core, j), hi = std::max(core, j);
  double same = 0;
  for (std::size_t a = 0; a < c.size(); a++)
    for (std::size_t b = a + 1; b < c.size(); b++) {
      bool equal = true;
      for (std::size_t k = lo; k <= hi && equal; k++) equal = h(c[a], k) >= 0 && h(c[a], k) == h(c[b], k);
      same += equal;
    }
  const double n = static_cast<double>(c.size());
  return same / (n * (n - 1) / 2);
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- a small example by hand ----------------------------------------------------------------------------
  {
    // Six haplotypes at five sites one unit apart; the core is site 2, where h0, h1, h2 carry 1.
    Haplotypes h;
    h.sequences = 6;
    h.sites = 5;
    const std::int8_t by_haplotype[6][5] = {{0, 1, 1, 0, 1}, {0, 1, 1, 0, 0}, {1, 1, 1, 0, 1}, {0, 0, 0, 1, 1}, {0, 0, 0, 1, 1}, {1, 0, 0, 0, 1}};
    h.alleles.resize(30);
    for (std::size_t i = 0; i < 6; i++)
      for (std::size_t j = 0; j < 5; j++) h.alleles[j * 6 + i] = by_haplotype[i][j];
    h.positions = {0, 1, 2, 3, 4};
    // Derived (1): {h0 h1 h2} stays whole to sites 1 and 3, and splits into a pair and one at
    // sites 0 and 4: EHH 1/3, 1, 1, 1, 1/3. Ancestral (0): {h3 h4 h5} splits at site 0 and at
    // site 3: EHH 1/3, 1, 1, 1/3, 1/3.
    const EhhCurve d = ehh(h, 2, 1), a = ehh(h, 2, 0);
    const std::vector<double> ed = {1.0 / 3, 1, 1, 1, 1.0 / 3}, ea = {1.0 / 3, 1, 1, 1.0 / 3, 1.0 / 3};
    s.require("EHH by hand: derived 1/3, 1, 1, 1, 1/3 and ancestral 1/3, 1, 1, 1/3, 1/3 over sites 0 … 4",
              d.ehh.size() == 5 && a.ehh.size() == 5 && calx::verification::max_difference(d.ehh, ed) < 1e-15 && calx::verification::max_difference(a.ehh, ea) < 1e-15);
    // iHH by the trapezoid rule, every walk reaching the ends (cut-off 0.05 never reached):
    // derived (1 + 2/3) + (1 + 2/3) = 10/3, ancestral (1 + 2/3) + (2/3 + 1/3) = 8/3.
    const Ihh id = ihh(h, 2, 1), ia = ihh(h, 2, 0);
    s.require("iHH by hand: 10/3 derived and 8/3 ancestral, both truncated by the ends of the data",
              std::fabs(id.total() - 10.0 / 3) < 1e-15 && std::fabs(ia.total() - 8.0 / 3) < 1e-15 && id.truncated && ia.truncated);
    IhsOptions o;
    o.min_maf = 0;
    o.keep_truncated = true;
    const std::vector<IhsSite> scan = ihs_scan(h, {}, o);
    bool found = false;
    for (const IhsSite& x : scan)
      if (x.site == 2) found = std::fabs(x.unstandardised - std::log(0.8)) < 1e-15 && std::fabs(x.frequency - 0.5) < 1e-15;
    s.require("the unstandardised iHS of the core: ln((8/3)/(10/3)) = ln 0.8, at derived frequency 1/2", found);
    o.keep_truncated = false;
    s.require("truncated cores are left out of the scan unless asked for", ihs_scan(h, {}, o).empty());
    // With the cut-off at 0.5 each walk stops at the first EHH ≤ 0.5, that interval included:
    // derived unchanged (10/3, now not truncated), ancestral (1 + 2/3) + 2/3 = 7/3.
    const Ihh cd = ihh(h, 2, 1, {}, 0.5), ca = ihh(h, 2, 0, {}, 0.5);
    s.require("with the cut-off at 0.5: 10/3 and 7/3, the walks stopped inside the data",
              std::fabs(cd.total() - 10.0 / 3) < 1e-15 && std::fabs(ca.total() - 7.0 / 3) < 1e-15 && !cd.truncated && !ca.truncated);
    // A map stretched twice doubles every integral.
    const std::vector<double> twice = {0, 2, 4, 6, 8};
    s.require("iHH is an integral over the map: twice the distances, twice the iHH", std::fabs(ihh(h, 2, 1, twice).total() - 20.0 / 3) < 1e-14);
  }

  // ---- the walk against the definition by pairs -------------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 40;
    spec.theta = 30;
    spec.rho = 60;
    Haplotypes h = evolution::population::haplotypes(coalescent(spec, 51));
    calx::random::Stream r(51, 0, 0, 0);
    for (auto& x : h.alleles)
      if (r.bernoulli(0.01)) x = -1;
    double worst = 0;
    bool shape = true;
    std::size_t compared = 0;
    for (std::size_t core = 0; core < h.sites; core += 7)
      for (const std::int8_t allele : {std::int8_t{0}, std::int8_t{1}}) {
        if (carriers(h, core, allele) < 2) continue;
        const EhhCurve c = ehh(h, core, allele);
        for (std::size_t k = 0; k < c.sites.size(); k++) {
          worst = std::max(worst, std::fabs(c.ehh[k] - ehh_by_pairs(h, core, allele, c.sites[k])));
          compared++;
          if (c.sites[k] == core) shape = shape && c.ehh[k] == 1;
          // Non-increasing away from the core on each side.
          if (k > 0 && c.sites[k] <= core) shape = shape && c.ehh[k] >= c.ehh[k - 1];
          if (k > 0 && c.sites[k - 1] >= core) shape = shape && c.ehh[k] <= c.ehh[k - 1];
        }
        worst = std::max(worst, std::fabs(ehh_at(h, core, allele, h.sites - 1) - ehh_by_pairs(h, core, allele, h.sites - 1)));
      }
    s.check("EHH by the partition walk = by pairs of carriers identical and typed (" + std::to_string(compared) + " values, 1% of alleles missing)", worst, 1e-15);
    s.require("EHH is 1 at the core and non-increasing away from it", shape);
  }

  // ---- neutral scans: threads and standardisation ------------------------------------------------------------
  CoalescentSpec neutral;
  neutral.samples = 100;
  neutral.theta = 60;
  neutral.rho = 400;
  std::vector<IhsSite> background;
  {
    calx::scheduler::ThreadPool pool(8);
    bool same = true;
    const auto t0 = std::chrono::steady_clock::now();
    for (std::uint32_t rep = 0; rep < 100; rep++) {
      const Haplotypes h = evolution::population::haplotypes(coalescent(neutral, 52, rep));
      const std::vector<IhsSite> scan = ihs_scan(h, {}, {}, &pool);
      if (rep < 3) {
        const std::vector<IhsSite> serial = ihs_scan(h);
        same = same && serial.size() == scan.size();
        for (std::size_t i = 0; same && i < scan.size(); i++)
          same = serial[i].site == scan[i].site && serial[i].unstandardised == scan[i].unstandardised;
      }
      background.insert(background.end(), scan.begin(), scan.end());
    }
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    s.require("a scan is the same on 1 and 8 threads", same);
    s.require("100 neutral scans of 100 haplotypes (θ = 60, ρ = 400) in under 5 s", ms < 5000, fmt(ms) + " ms, " + std::to_string(background.size()) + " core sites");
    std::vector<IhsSite> z = background;
    const IhsStandardisation fit = standardise(z, 20);
    std::vector<double> sum(20, 0), sum2(20, 0), count(20, 0);
    for (const IhsSite& x : z)
      if (std::isfinite(x.ihs)) {
        const std::size_t b = fit.bin(x.frequency);
        sum[b] += x.ihs;
        sum2[b] += x.ihs * x.ihs;
        count[b]++;
      }
    double mean_err = 0, var_err = 0;
    std::size_t bins = 0;
    for (std::size_t b = 0; b < 20; b++)
      if (count[b] >= 2) {
        mean_err = std::max(mean_err, std::fabs(sum[b] / count[b]));
        var_err = std::max(var_err, std::fabs(sum2[b] / count[b] - (sum[b] / count[b]) * (sum[b] / count[b]) - 1));
        bins++;
      }
    s.check("iHS has mean 0 in each of the " + std::to_string(bins) + " bins of derived frequency, by construction (largest |mean|)", mean_err, 1e-12);
    s.check("and variance 1 (largest departure)", var_err, 1e-12);
  }

  // ---- the signal of a partial sweep --------------------------------------------------------------------------
  {
    const IhsStandardisation fit = fit_ihs_standardisation(background, 20);
    std::size_t large = 0, finite = 0;
    for (const IhsSite& x : background) {
      const double v = fit.apply(x.frequency, x.unstandardised);
      if (!std::isfinite(v)) continue;
      finite++;
      large += std::fabs(v) > 2;
    }
    SweepSpec spec;
    spec.individuals = 500;
    spec.samples = 100;
    spec.theta = 60;
    spec.rho = 400;
    spec.s = 0.05;
    spec.stop = 0.7;
    std::size_t sweeps = 0, signal = 0;
    double mean = 0;
    // The same options as the background's: a selected site whose walk reaches the ends is left out.
    for (std::uint32_t rep = 0; rep < 40; rep++) {
      const Sweep w = evolution::population::sweep(spec, 53, rep);
      if (w.selected == Sweep::npos) continue;
      for (const IhsSite& x : ihs_scan(w.haplotypes))
        if (x.site == w.selected) {
          const double v = fit.apply(x.frequency, x.unstandardised);
          if (!std::isfinite(v)) continue;
          sweeps++;
          signal += v < -2;
          mean += v;
        }
    }
    mean /= static_cast<double>(sweeps);
    const double neutral_rate = static_cast<double>(large) / static_cast<double>(finite), sweep_rate = static_cast<double>(signal) / static_cast<double>(sweeps);
    s.require("neutral sites rarely score |iHS| > 2 under the background's standardisation (below 10%)", neutral_rate < 0.1,
              fmt(neutral_rate) + " of " + std::to_string(finite));
    s.require("partial sweeps to 0.7 (2N = 1000, s = 0.05, θ = 60, ρ = 400): most selected sites score iHS < −2", sweep_rate > 0.6,
              fmt(sweep_rate) + " of " + std::to_string(sweeps) + ", mean iHS " + fmt(mean));
  }
  return s.done();
}
