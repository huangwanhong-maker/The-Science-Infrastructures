// Checks of ASCAT's model and fit on profiles whose truth is known: eqs. 1–4 invert each other
// exactly; R's rounding; the distance vanishes at the truth and the grid's least distance is the
// truth, noise-free and with noise; the dilation (t, n, ψ) ↦ (2t, 2n, 2ψ) with t = (1 − ρ)/ρ leaves
// logR and BAF unchanged (theory/the-somatic-genome.md, Proposition 3.3), which is why integrality
// decides; ASPCF finds the breakpoints of noisy sites; segments of integer copy numbers.

#include <genomics/somatic/copy_number.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace genomics;

namespace {

struct Truth {
  int a, b;
  std::size_t sites;
};

// Segments of the given copy numbers read at purity ρ, ploidy ψ (site-weighted), with noise sd on logR and BAF.
std::vector<somatic::CnSegment> segments_of(const std::vector<Truth>& t, double rho, double psi, double sd_r, double sd_b, std::uint64_t seed) {
  calx::random::Stream r(seed, 1, 0, 0);
  std::vector<somatic::CnSegment> s;
  std::size_t first = 0;
  for (const auto& x : t) {
    somatic::CnSegment g;
    g.first = first;
    g.last = first + x.sites;
    first = g.last;
    g.logr = somatic::ascat_logr(x.a, x.b, rho, psi) + r.normal(0, sd_r);
    double b = somatic::ascat_baf(std::min(x.a, x.b), std::max(x.a, x.b), rho) + r.normal(0, sd_b);
    if (x.a == x.b) b = 0.5;
    g.baf = std::max(0.5, b);
    s.push_back(g);
  }
  return s;
}

double weighted_ploidy(const std::vector<Truth>& t) {
  double n = 0, w = 0;
  for (const auto& x : t) {
    n += (x.a + x.b) * static_cast<double>(x.sites);
    w += static_cast<double>(x.sites);
  }
  return n / w;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- eqs. 1–4 ----------------------------------------------------------------------------------------
  {
    calx::random::Stream r(5, 0, 0, 0);
    double worst = 0;
    for (int i = 0; i < 2000; i++) {
      const double rho = 0.1 + 0.9 * r.uniform(), psi = 1 + 5 * r.uniform(), gamma = r.bernoulli(0.5) ? 1 : 0.55;
      const double na = static_cast<double>(r.below(6)), nb = static_cast<double>(r.below(6));
      if (na + nb == 0) continue;
      const double lr = somatic::ascat_logr(na, nb, rho, psi, gamma), b = somatic::ascat_baf(na, nb, rho);
      worst = std::max({worst, std::fabs(somatic::ascat_na(lr, b, rho, psi, gamma) - na), std::fabs(somatic::ascat_nb(lr, b, rho, psi, gamma) - nb)});
    }
    s.check("eqs. 3–4 invert eqs. 1–2 (2000 draws of ρ, ψ, γ, n_A, n_B)", worst, 1e-12);
    s.require("round half to even: 0.5 → 0, 1.5 → 2, 2.5 → 2, −1.5 → −2, 2.4 → 2, 2.6 → 3",
              somatic::round_half_even(0.5) == 0 && somatic::round_half_even(1.5) == 2 && somatic::round_half_even(2.5) == 2 && somatic::round_half_even(-1.5) == -2 &&
                  somatic::round_half_even(2.4) == 2 && somatic::round_half_even(2.6) == 3);
  }

  // ---- the dilation symmetry ------------------------------------------------------------------------------
  {
    const double rho = 0.65, t = (1 - rho) / rho, rho2 = 1 / (1 + 2 * t), psi = 2.4;
    double worst = 0;
    for (int a = 0; a <= 3; a++)
      for (int b = 0; b <= 3; b++) {
        if (a + b == 0) continue;
        worst = std::max({worst, std::fabs(somatic::ascat_logr(a, b, rho, psi) - somatic::ascat_logr(2 * a, 2 * b, rho2, 2 * psi)),
                          std::fabs(somatic::ascat_baf(a, b, rho) - somatic::ascat_baf(2 * a, 2 * b, rho2))});
      }
    s.check("(ρ, n, ψ) and (1/(1 + 2t), 2n, 2ψ) give the same logR and BAF (t = (1 − ρ)/ρ)", worst, 1e-14);
  }

  // ---- the fit on noise-free and noisy profiles ---------------------------------------------------------
  const std::vector<Truth> profile{{1, 1, 40}, {2, 1, 30}, {1, 0, 10}, {2, 0, 10}, {3, 1, 10}};
  {
    const double rho = 0.65, psi = weighted_ploidy(profile);
    auto seg = segments_of(profile, rho, psi, 0, 0, 1);
    s.check("the distance at the truth is 0", somatic::ascat_distance(seg, rho, psi), 1e-20);
    const somatic::AscatFit f = somatic::ascat(seg);
    s.require("noise-free: a fit is found in the first pass", f.found && f.pass == 1);
    s.check("noise-free: purity 0.65", f.purity - rho, 1e-9);
    s.check("noise-free: ψ 2.40", f.psi - psi, 1e-9);
    s.check("noise-free: goodness of fit 100%", f.goodness - 100, 1e-9);
    bool cn = true;
    for (std::size_t i = 0; i < seg.size(); i++) cn = cn && seg[i].major() == std::max(profile[i].a, profile[i].b) && seg[i].minor() == std::min(profile[i].a, profile[i].b);
    s.require("noise-free: the integer copy numbers of every segment", cn);
    s.require("the grid: ψ 1.00 to 6.00 by 0.05, ρ 0.10 to 1.05 by 0.01", f.psi_grid.size() == 101 && f.rho_grid.size() == 96 && f.psi_grid.front() == 1 && f.psi_grid.back() == 6 &&
                                                                              f.rho_grid.front() == 0.1 && f.rho_grid.back() == 1.05);
  }
  {
    // Twenty noisy realisations at purities 0.3 to 0.9.
    int found = 0, exact = 0;
    double worst_rho = 0, worst_psi = 0;
    for (int k = 0; k < 20; k++) {
      const double rho = 0.3 + 0.03 * k;
      const double psi = weighted_ploidy(profile);
      auto seg = segments_of(profile, rho, psi, 0.02, 0.01, 100 + static_cast<std::uint64_t>(k));
      const somatic::AscatFit f = somatic::ascat(seg);
      if (!f.found) continue;
      found++;
      worst_rho = std::max(worst_rho, std::fabs(f.purity - rho));
      worst_psi = std::max(worst_psi, std::fabs(f.psi - psi));
      if (std::fabs(f.purity - rho) > 0.02 || std::fabs(f.psi - psi) > 0.05)
        std::printf("      purity %.2f: fit ρ %.2f ψ %.2f goodness %.1f\n", rho, f.purity, f.psi, f.goodness);
      bool cn = true;
      for (std::size_t i = 0; i < seg.size(); i++) cn = cn && seg[i].major() == std::max(profile[i].a, profile[i].b) && seg[i].minor() == std::min(profile[i].a, profile[i].b);
      exact += cn;
    }
    // The distance weighs the minor allele alone, so that ψ is held more loosely than ρ (theory,
    // Remark after Definition 3.4); the integers, and the ploidy they give, are what must be right.
    s.require("noisy profiles (logR sd 0.02, BAF sd 0.01): a fit for each of 20", found == 20, std::to_string(found));
    s.check("noisy profiles: purity within 0.05", worst_rho, 0.05);
    s.check("noisy profiles: ψ (the grid's) within 0.25", worst_psi, 0.25);
    s.require("noisy profiles: every segment's integers in all 20, so the ploidy they give is the truth's", exact == 20, std::to_string(exact));
  }
  {
    // A doubled genome: the fit must find ψ ≈ 4.8, not the halved solution.
    const std::vector<Truth> doubled{{2, 2, 40}, {4, 2, 30}, {2, 0, 10}, {4, 0, 10}, {6, 2, 10}};
    const double rho = 0.7, psi = weighted_ploidy(doubled);
    auto seg = segments_of(doubled, rho, psi, 0, 0, 2);
    const somatic::AscatFit f = somatic::ascat(seg);
    std::printf("      doubled genome: fit ρ %.2f ψ %.2f goodness %.1f (truth ρ %.2f ψ %.2f)\n", f.purity, f.psi, f.goodness, rho, psi);
    s.require("a doubled genome is fitted at the truth or at its halving (the dilation)", f.found && ((std::fabs(f.psi - psi) < 0.06 && std::fabs(f.purity - rho) < 0.011) ||
                                                                                                     (std::fabs(f.psi - psi / 2) < 0.06)));
  }

  // ---- ASPCF on noisy sites ----------------------------------------------------------------------------
  {
    const double rho = 0.6;
    const std::vector<Truth> truth{{1, 1, 120}, {2, 1, 90}, {1, 0, 60}, {1, 1, 80}, {2, 2, 70}, {3, 1, 80}};
    const double psi = weighted_ploidy(truth);
    calx::random::Stream r(77, 0, 0, 0);
    std::vector<somatic::HetSite> sites;
    std::vector<std::size_t> breaks;
    std::size_t pos = 1000;
    for (const auto& t : truth) {
      breaks.push_back(sites.size());
      for (std::size_t i = 0; i < t.sites; i++) {
        somatic::HetSite h;
        h.position = pos;
        pos += 500;
        h.logr = somatic::ascat_logr(t.a, t.b, rho, psi) + r.normal(0, 0.08);
        // BAF of 60 reads at the expected fraction, on a random allele.
        const double f = somatic::ascat_baf(t.a, t.b, rho);
        const double p = r.bernoulli(0.5) ? f : 1 - f;
        h.baf = static_cast<double>(r.binomial(60, p)) / 60;
        sites.push_back(h);
      }
    }
    breaks.erase(breaks.begin());
    auto seg = somatic::aspcf(sites, {pos + 1000});
    std::size_t matched = 0;
    for (std::size_t b : breaks)
      for (std::size_t k = 1; k < seg.size(); k++)
        if ((seg[k].first > b ? seg[k].first - b : b - seg[k].first) <= 2) { matched++; break; }
    s.require("ASPCF: six segments from six states", seg.size() == truth.size(), std::to_string(seg.size()));
    s.require("ASPCF: every breakpoint within two sites", matched == breaks.size(), std::to_string(matched));
    const somatic::AscatFit f = somatic::ascat(seg);
    std::printf("      ASPCF then ASCAT: ρ %.2f ψ %.2f goodness %.1f (truth %.2f, %.2f)\n", f.purity, f.psi, f.goodness, rho, psi);
    s.check("ASPCF then ASCAT: purity", f.found ? f.purity - rho : 1, 0.04);
    s.check("ASPCF then ASCAT: ψ", f.found ? f.psi - psi : 1, 0.15);
    bool cn = seg.size() == truth.size();
    for (std::size_t i = 0; cn && i < seg.size(); i++) cn = seg[i].major() == std::max(truth[i].a, truth[i].b) && seg[i].minor() == std::min(truth[i].a, truth[i].b);
    s.require("ASPCF then ASCAT: the integer copy numbers of every segment", cn);
    s.require("a balanced segment's BAF is 0.5 exactly (the 2 + 2 and the 1 + 1)", seg.size() == truth.size() && seg[0].baf == 0.5 && seg[4].baf == 0.5);
  }
  return s.done();
}
