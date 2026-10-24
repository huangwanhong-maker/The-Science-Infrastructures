// Checks of linkage disequilibrium: D, D′ (Lewontin 1964) and r² (Hill and Robertson 1968) on
// haplotypes worked by hand, including complete repulsion and missing calls; bit-packed pairs
// equal the plain count; from unphased genotypes the EM of haplotype frequencies (Excoffier and
// Slatkin 1995; Hill 1974) counts exactly when no individual is a double heterozygote, attains
// the maximum of the likelihood found by a fine grid over x_AB (the likelihood being a function
// of x_AB alone, the allele frequencies fixed), and on genotypes paired from coalescent
// haplotypes recovers the phased r²; the decay of r² with distance is the same from the phased
// genotypes as from their haplotypes, and the same on 1 and 8 threads.

#include <evolution/population/coalescent.hpp>
#include <evolution/population/linkage.hpp>

#include <calx/random/stream.hpp>
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

  // ---- by hand ---------------------------------------------------------------------------------------
  {
    // Ten haplotypes: AB × 4, Ab × 1, aB × 1, ab × 4 (A, B the alternate alleles, coded 1).
    const std::vector<std::int8_t> a = {1, 1, 1, 1, 1, 0, 0, 0, 0, 0}, b = {1, 1, 1, 1, 0, 1, 0, 0, 0, 0};
    const Linkage l = linkage_phased(a, b);
    // p_A = p_B = ½, p_AB = 0.4: D = 0.15, D_max = min(p_A q_B, q_A p_B) = 0.25, r² = 0.15²/0.25².
    s.check("D = p_AB − p_A p_B = 0.15 for AB×4, Ab, aB, ab×4", l.d - 0.15, 1e-15);
    s.check("D′ = D/D_max = 0.6", l.d_prime - 0.6, 1e-15);
    s.check("r² = D²/(p_A q_A p_B q_B) = 0.36", l.r2 - 0.36, 1e-15);
    const Linkage repulsion = linkage_phased(std::vector<std::int8_t>{1, 1, 0, 0}, std::vector<std::int8_t>{0, 0, 1, 1});
    s.require("complete repulsion: D = −¼, D′ = −1, r² = 1", repulsion.d == -0.25 && repulsion.d_prime == -1 && repulsion.r2 == 1);
    const Linkage missing = linkage_phased(std::vector<std::int8_t>{1, -1, 0, 0, 1}, std::vector<std::int8_t>{1, 1, 0, -1, 1});
    s.require("haplotypes missing at either site are left out", missing.n == 3 && std::fabs(missing.r2 - 1) < 1e-12);
    s.require("r² and D′ are undefined at a monomorphic site", std::isnan(linkage_phased(std::vector<std::int8_t>{1, 1}, std::vector<std::int8_t>{0, 1}).r2));
  }

  // ---- packed ------------------------------------------------------------------------------------------
  {
    // 130 haplotypes (three words, the last partly used), 20 sites, 5% missing.
    Haplotypes h;
    h.sequences = 130;
    h.sites = 20;
    calx::random::Stream r(1, 0x4c444c44, 0, 0);
    for (std::size_t j = 0; j < h.sites; j++) {
      const double p = r.uniform(0.1, 0.9);
      for (std::size_t k = 0; k < h.sequences; k++) {
        // Correlated with the previous site, so that pairs are in disequilibrium.
        std::int8_t x = j && r.uniform() < 0.6 && h.alleles[(j - 1) * h.sequences + k] >= 0 ? h.alleles[(j - 1) * h.sequences + k]
                                                                                         : static_cast<std::int8_t>(r.uniform() < p);
        if (r.uniform() < 0.05) x = MISSING;
        h.alleles.push_back(x);
      }
      h.positions.push_back(static_cast<double>(j));
    }
    const PackedHaplotypes packed(h);
    double worst = 0;
    for (std::size_t i = 0; i < h.sites; i++)
      for (std::size_t j = 0; j < h.sites; j++) {
        const Linkage a = packed.linkage(i, j), b = linkage_phased(h.site(i), h.site(j));
        worst = std::max({worst, std::fabs(a.d - b.d), std::fabs(a.r2 - b.r2), std::fabs(a.d_prime - b.d_prime)});
        if (a.n != b.n) worst = 1;
      }
    s.check("bit-packed pairs equal the plain count, across words and with missing calls (largest difference)", worst, 1e-15);
  }

  // ---- the EM of haplotype frequencies ------------------------------------------------------------------
  {
    // No double heterozygotes: the haplotypes are known, and one step counts them.
    const std::vector<std::int8_t> g1 = {2, 2, 1, 0, 0, 1, 2}, g2 = {2, 1, 0, 0, 1, 2, 0};
    const HaplotypeEstimate e = linkage_unphased(g1, g2);
    // The pairs (2, 2), (2, 1), (1, 0), (0, 0), (0, 1), (1, 2), (2, 0) are AB/AB, AB/Ab, Ab/ab,
    // ab/ab, aB/ab, AB/aB, Ab/Ab: AB 4, Ab 4, aB 2, ab 4, fourteen haplotypes in all.
    s.check("without double heterozygotes the EM counts the haplotypes exactly: x_AB = 4/14", e.ld.p_ab - 4.0 / 14, 1e-12);

    // Random tables: the EM's x_AB attains the grid's maximum of the likelihood.
    calx::random::Stream r(2, 0x454d454d, 0, 0);
    double worst_ll = 0, worst_x = 0;
    for (int trial = 0; trial < 50; trial++) {
      GenotypeTable t;
      for (auto& x : t.n) x = r.below(30);
      t.n[4] = 5 + r.below(60);  // double heterozygotes
      const HaplotypeEstimate em = linkage_unphased(t);
      const double pA = em.ld.p_a, pB = em.ld.p_b;
      const double lo = std::max(0.0, pA + pB - 1), hi = std::min(pA, pB);
      double best = -INFINITY, at = lo;
      for (int k = 0; k <= 200000; k++) {
        const double x = lo + (hi - lo) * k / 200000.0;
        const double l = haplotype_log_likelihood(t, {x, pA - x, pB - x, 1 - pA - pB + x});
        if (l > best) {
          best = l;
          at = x;
        }
      }
      worst_ll = std::max(worst_ll, best - em.log_likelihood);
      worst_x = std::max(worst_x, std::fabs(at - em.ld.p_ab));
    }
    s.check("the EM attains the likelihood's maximum over a grid of 200001 values of x_AB (50 random tables; largest shortfall)", worst_ll, 1e-6);
    s.check("at the same x_AB (largest difference)", worst_x, 2e-4);

    // Genotypes from coalescent haplotypes: the EM's r² against the phased r².
    CoalescentSpec spec;
    spec.samples = 400;
    spec.theta = 20;
    spec.rho = 20;
    const Haplotypes h = haplotypes(coalescent(spec, 3));
    const Genotypes g = diploids(h, 1000000);
    double sum = 0, pairs = 0;
    for (std::size_t i = 0; i < g.sites; i++) {
      const auto ci = count_site(g.site(i));
      if (ci.first < 40 || ci.first > 360) continue;
      for (std::size_t j = i + 1; j < std::min(g.sites, i + 30); j++) {
        const auto cj = count_site(g.site(j));
        if (cj.first < 40 || cj.first > 360) continue;
        const double phased = linkage_phased(g.haplotypes_at(i), g.haplotypes_at(j)).r2;
        const double unphased = linkage_unphased(g.site(i), g.site(j)).ld.r2;
        sum += std::fabs(phased - unphased);
        pairs++;
      }
    }
    s.check("on 200 diploids paired from coalescent haplotypes the EM recovers the phased r² (mean absolute difference over " +
                std::to_string(static_cast<int>(pairs)) + " pairs of common sites)",
            sum / pairs, 0.02);
  }

  // ---- decay -----------------------------------------------------------------------------------------
  {
    CoalescentSpec spec;
    spec.samples = 200;
    spec.theta = 200;
    spec.rho = 200;
    const Haplotypes h = haplotypes(coalescent(spec, 4));
    const Genotypes g = diploids(h, 1000000);
    calx::scheduler::ThreadPool one(1), eight(8);
    const auto t0 = std::chrono::steady_clock::now();
    const auto a = ld_decay(h, 0.2, 20, &eight);
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    const auto b = ld_decay(h, 0.2, 20, &one);
    const auto c = ld_decay(g, 200000, 20, &eight);  // the same pairs, in base pairs
    bool same = true, from_genotypes = true;
    std::size_t pairs = 0;
    for (std::size_t k = 0; k < 20; k++) {
      same = same && a[k].pairs == b[k].pairs && (a[k].mean_r2 == b[k].mean_r2 || (std::isnan(a[k].mean_r2) && std::isnan(b[k].mean_r2)));
      from_genotypes = from_genotypes && std::fabs(a[k].mean_r2 - c[k].mean_r2) < 0.02;
      pairs += a[k].pairs;
    }
    s.require("the decay of r² is the same on 1 and 8 threads", same, std::to_string(pairs) + " pairs in " + fmt(ms) + " ms");
    s.require("and from the phased genotypes as from their haplotypes (positions rounded to base pairs)", from_genotypes);
    s.require("r² falls from the nearest bin to the farthest", a[0].mean_r2 > 2 * a[19].mean_r2, fmt(a[0].mean_r2) + " → " + fmt(a[19].mean_r2));
    const std::vector<double> m = ld_matrix(PackedHaplotypes(h), 0, 10);
    bool symmetric = true;
    for (std::size_t i = 0; i < 10; i++)
      for (std::size_t j = 0; j < 10; j++) symmetric = symmetric && (m[i * 10 + j] == m[j * 10 + i] || std::isnan(m[i * 10 + j]));
    s.require("the matrix of r² for drawing is symmetric with 1 on the diagonal", symmetric && std::fabs(m[0] - 1) < 1e-12 && std::fabs(m[11] - 1) < 1e-12);
  }
  return s.done();
}
