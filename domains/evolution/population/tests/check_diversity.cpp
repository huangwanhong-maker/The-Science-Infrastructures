// Checks of diversity and differentiation: S, θ_W, π and θ_H on a sample worked by hand, π
// against the mean of the pairwise differences counted pair by pair; Tajima's D tends to
// Tajima's bounds D_min = (2/n − 1/a₁)/√e₂ and D_max = (n/(2(n − 1)) − 1/a₁)/√e₂ (n even) as S
// grows (Tajima 1989); the projected spectrum is the identity on complete data and equals the
// enumeration of every subsample; folding; polarising by an outgroup; Hudson's F_ST (Bhatia et
// al. 2013) and Weir and Cockerham's variance components on the worked example of scikit-allel's
// documentation (computed there from the same formulas; biallelic sites only); on simulated
// populations diverged by Balding–Nichols' F both estimate F; windows partition the sites, and
// their statistics are the sums over their ranges.

#include <evolution/population/diversity.hpp>
#include <evolution/population/structure.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace evolution::population;

int main() {
  calx::verification::Suite s;

  // ---- by hand -------------------------------------------------------------------------------------
  {
    // Four sequences at four sites; derived counts 1, 2, 3 and 0.
    Haplotypes h;
    h.sequences = 4;
    h.sites = 4;
    h.alleles = {1, 0, 0, 0, /**/ 1, 1, 0, 0, /**/ 0, 1, 1, 1, /**/ 0, 0, 0, 0};
    h.positions = {1, 2, 3, 4};
    const AlleleCounts c = allele_counts(h);
    const double a4 = 1 + 1.0 / 2 + 1.0 / 3;
    s.require("S counts the segregating sites", segregating_sites(c) == 3);
    s.check("θ_W = S/a_n = 3/a₄", theta_watterson(c) - 3 / a4, 1e-15);
    s.check("π = Σ 2k(n − k)/(n(n − 1)) = ½ + ⅔ + ½", nucleotide_diversity(c) - (0.5 + 2.0 / 3 + 0.5), 1e-15);
    double diff = 0;
    for (std::size_t x = 0; x < 4; x++)
      for (std::size_t y = x + 1; y < 4; y++)
        for (std::size_t j = 0; j < 4; j++) diff += h(x, j) != h(y, j);
    s.check("π is the mean number of differences over the six pairs", nucleotide_diversity(c) - diff / 6, 1e-15);
    s.check("θ_H = Σ 2k²/(n(n − 1)) = (2 + 8 + 18)/12", theta_h(c) - 28.0 / 12, 1e-15);
    s.check("Fay and Wu's H = π − θ_H", fay_wu_h(c) - (5.0 / 3 - 28.0 / 12), 1e-15);
    s.check("Tajima's D of the counts is D(n, S, π)", tajima_d(c) - tajima_d(4, 3, 5.0 / 3), 1e-14);
    // Missing calls: θ_W sums 1/a_{n_j} site by site.
    AlleleCounts m = c;
    m.total[1] = 3;
    s.check("with missing calls θ_W sums 1/a_{n_j} over segregating sites", theta_watterson(m) - (2 / a4 + 1 / (1 + 0.5)), 1e-15);
  }

  // ---- Tajima's bounds -----------------------------------------------------------------------------
  {
    const std::size_t n = 10;
    const TajimaConstants t = tajima_constants(n);
    const double S = 1e7;
    const double d_min = (2.0 / n - 1 / t.a1) / std::sqrt(t.e2), d_max = (n / (2.0 * (n - 1)) - 1 / t.a1) / std::sqrt(t.e2);
    // All singletons: π = 2S/n; all at k = n/2: π = S n/(2(n − 1)).
    s.check("D tends to D_min = (2/n − 1/a₁)/√e₂ when every site is a singleton (Tajima 1989)", tajima_d(n, S, 2 * S / n) - d_min, 1e-5);
    s.check("and to D_max = (n/(2(n − 1)) − 1/a₁)/√e₂ when every site is at n/2", tajima_d(n, S, S * n / (2.0 * (n - 1))) - d_max, 1e-5);
    s.require("D is NaN without segregating sites", std::isnan(tajima_d(n, 0, 0)));
  }

  // ---- spectra --------------------------------------------------------------------------------------
  {
    AlleleCounts c;
    c.alt = {1, 3, 5, 0, 2};
    c.total = {6, 6, 6, 6, 6};
    const std::vector<double> xi = site_frequency_spectrum(c, 6);
    s.require("the spectrum of complete data counts the sites by derived count", xi == std::vector<double>{1, 1, 1, 1, 0, 1, 0});
    const std::vector<double> eta = fold(xi);
    s.require("the folded spectrum adds i and n − i", eta == std::vector<double>{1, 2, 1, 1});
    // Projection: a site of 8 sequences with 3 derived, to 5, against every subsample.
    AlleleCounts one;
    one.alt = {3};
    one.total = {8};
    const std::vector<double> projected = site_frequency_spectrum(one, 5);
    std::vector<double> enumerated(6, 0);
    double subsamples = 0;
    for (unsigned mask = 0; mask < 256; mask++) {
      if (__builtin_popcount(mask) != 5) continue;
      enumerated[static_cast<std::size_t>(__builtin_popcount(mask & 0b111))]++;  // the first three sequences derived
      subsamples++;
    }
    double worst = 0;
    for (std::size_t i = 0; i < 6; i++) worst = std::max(worst, std::fabs(projected[i] - enumerated[i] / subsamples));
    s.check("projection to a smaller sample is the hypergeometric law of every subsample (largest difference)", worst, 1e-15);
    AlleleCounts big;
    big.alt = {37, 120, 1};
    big.total = {200, 180, 150};
    double total = 0;
    for (const double x : site_frequency_spectrum(big, 40)) total += x;
    s.check("each projected site contributes one site to the spectrum", total - 3, 1e-12);
    const AlleleCounts d = polarise(c, std::vector<std::int8_t>{0, 1, 1, 0, -1});
    s.require("polarising: an alternate ancestral allele turns k into n − k; an unknown one drops the site",
              d.alt == std::vector<std::uint32_t>{1, 3, 1, 0, 0} && d.total == std::vector<std::uint32_t>{6, 6, 6, 6, 0});
  }

  // ---- F_ST on a worked example --------------------------------------------------------------------
  {
    // scikit-allel (Miles et al.), the docstrings of weir_cockerham_fst and hudson_fst
    // (allel/stats/fst.py, retrieved 2026-10-24): four individuals in two populations {0, 1},
    // {2, 3}; their sites 0, 1, 2 and 4 (site 3 is triallelic):
    //   [0/0 0/0 1/1 1/1], [0/1 0/1 0/1 0/1], [0/0 0/0 0/0 0/0], [0/0 1/1 0/1 ./.]
    // a = 0.5, 0, 0, −0.375; b = 0, −0.25, 0, 0.41666667; c = 0, 0.5, 0, 0.16666667;
    // Hudson's numerators 1, −0.16666667, 0, −0.33333333 and denominators 1, 0.5, 0, 0.5.
    const Genotypes g = genotypes_from_dosages(4, {0, 0, 2, 2, /**/ 1, 1, 1, 1, /**/ 0, 0, 0, 0, /**/ 0, 2, 1, -1}, {0, 1, 2, 4});
    const std::vector<std::uint32_t> pop = {0, 0, 1, 1};
    const double a[4] = {0.5, 0, 0, -0.375}, b[4] = {0, -0.25, 0, 0.41666667}, cc[4] = {0, 0.5, 0, 0.16666667};
    const double num[4] = {1, -0.16666667, 0, -0.33333333}, den[4] = {1, 0.5, 0, 0.5};
    const auto m = members(pop);
    const AlleleCounts c0 = allele_counts(g, m[0]), c1 = allele_counts(g, m[1]);
    double worst_wc = 0, worst_h = 0;
    for (std::size_t j = 0; j < 4; j++) {
      std::vector<std::uint32_t> n(2, 0), alt(2, 0), het(2, 0);
      for (std::size_t i = 0; i < 4; i++)
        if (g(i, j) >= 0) {
          n[pop[i]]++;
          alt[pop[i]] += static_cast<std::uint32_t>(g(i, j));
          het[pop[i]] += g(i, j) == 1;
        }
      const VarianceComponents v = weir_cockerham_site(n, alt, het);
      worst_wc = std::max({worst_wc, std::fabs(v.a - a[j]), std::fabs(v.b - b[j]), std::fabs(v.c - cc[j])});
      const Ratio r = hudson_site(c0.alt[j], c0.total[j], c1.alt[j], c1.total[j]);
      worst_h = std::max({worst_h, std::fabs(r.numerator - num[j]), std::fabs(r.denominator - den[j])});
    }
    s.check("Weir and Cockerham's a, b, c at each site of the worked example (largest difference)", worst_wc, 5e-9);
    s.check("Hudson's numerator and denominator at each site (largest difference)", worst_h, 5e-9);
    s.check("F_ST as a ratio of averages: Σ N / Σ D = 0.5/2", fst_hudson(c0, c1).value() - 0.25, 1e-15);
    s.check("θ̂ = Σ a / Σ (a + b + c) over the four sites", fst_weir_cockerham(g, pop).value() - 0.125 / (0.125 + 0.16666667 + 0.66666667), 1e-8);
  }

  // ---- F_ST of simulated populations ---------------------------------------------------------------
  {
    StructureSpec spec;
    spec.populations = 2;
    spec.sites = 10000;
    spec.fst = {0.1};
    spec.unadmixed = {100, 100};
    const StructuredSample x = simulate_structure(spec, 31);
    const auto m = members(x.population);
    const double hudson = fst_hudson(allele_counts(x.genotypes, m[0]), allele_counts(x.genotypes, m[1])).value();
    const double wc = fst_weir_cockerham(x.genotypes, x.population).value();
    // Two populations each diverged from the ancestor by F: E[(p₁ − p₂)²] = 2F pq and
    // E[p₁q₂ + p₂q₁] = 2pq, so both estimators estimate F.
    s.check("Hudson's F_ST of two populations diverged by Balding–Nichols' F = 0.1 (absolute error)", hudson - 0.1, 0.01);
    s.check("and Weir and Cockerham's θ̂ (absolute error)", wc - 0.1, 0.01);
  }

  // ---- windows -----------------------------------------------------------------------------------------
  {
    StructureSpec spec;
    spec.populations = 2;
    spec.sites = 1000;
    spec.unadmixed = {30, 30};
    spec.spacing = 997;
    const StructuredSample x = simulate_structure(spec, 41);
    const auto w = windows(x.genotypes, 50000, 50000);
    std::size_t covered = 0;
    bool contiguous = true;
    for (std::size_t k = 0; k < w.size(); k++) {
      covered += w[k].last - w[k].first;
      contiguous = contiguous && (k == 0 || w[k].first == w[k - 1].last);
    }
    s.require("adjacent windows partition the sites", covered == spec.sites && contiguous && w.back().end > x.genotypes.positions.back(),
              std::to_string(w.size()) + " windows");
    const auto stats = window_statistics(x.genotypes, 50000, 50000, x.population);
    const AlleleCounts all = allele_counts(x.genotypes);
    double pi = 0;
    for (const auto& st : stats) pi += st.pi * static_cast<double>(st.window.end - st.window.begin);
    s.check("π per base, times the windows' lengths, adds up to π of all sites (relative error)", pi / nucleotide_diversity(all) - 1, 1e-12);
    bool fst = true;
    for (const auto& st : stats) fst = fst && st.sites > 0 && std::isfinite(st.fst_hudson) && std::isfinite(st.fst_weir_cockerham) && std::isfinite(st.tajima_d);
    s.require("each window has its D and both F_ST", fst);
    const auto sliding = windows(x.genotypes, 50000, 10000);
    s.require("sliding windows overlap by size − step", sliding.size() > 4 * w.size() && sliding[1].begin == 10000 && sliding[1].end == 60000);
  }
  return s.done();
}
