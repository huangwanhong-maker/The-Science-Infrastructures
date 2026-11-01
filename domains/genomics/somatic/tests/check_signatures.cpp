// Checks of signature analysis on a simulated cohort (simulate::simulate_catalogues, three of the
// synthetic profiles): NNLS refitting exact on noise-free mixtures; de novo extraction by tyche's
// NMF with its stability, the rank chosen, the extracted signatures matched to the truth by cosine;
// refitted exposures against the truth; the bootstrap's intervals and their coverage.

#include <genomics/simulate/tumour.hpp>
#include <genomics/somatic/signatures.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace genomics;

int main() {
  calx::verification::Suite s;
  const auto synthetic = simulate::synthetic_signatures();
  const std::vector<std::vector<double>> truth{synthetic[0].profile, synthetic[1].profile, synthetic[3].profile};

  // ---- refitting -------------------------------------------------------------------------------------
  {
    const std::vector<double> e{300, 120, 80};
    std::vector<double> c(96, 0);
    for (std::size_t j = 0; j < 3; j++)
      for (std::size_t k = 0; k < 96; k++) c[k] += e[j] * truth[j][k];
    const auto f = somatic::refit(c, truth);
    double worst = 0;
    for (std::size_t j = 0; j < 3; j++) worst = std::max(worst, std::fabs(f.exposures[j] - e[j]));
    s.check("NNLS refit of a noise-free mixture: exposures exact", worst, 1e-8);
    s.check("its reconstruction's cosine is 1", f.cosine - 1, 1e-12);
  }

  // ---- a cohort --------------------------------------------------------------------------------------
  calx::random::Stream r(31, 0, 0, 0);
  std::vector<std::vector<double>> exposures;
  for (int i = 0; i < 40; i++) {
    std::vector<double> e(3);
    for (std::size_t j = 0; j < 3; j++) e[j] = r.bernoulli(0.75) ? r.gamma(2, 150) : 0;
    exposures.push_back(e);
  }
  const auto cohort = simulate::simulate_catalogues(truth, exposures, 32);
  tyche::multivariate::StabilityOptions so;
  so.replicates = 20;
  so.nmf.max_iterations = 4000;
  so.nmf.tolerance = 1e-9;
  so.seed = 33;
  const auto ex = somatic::extract_signatures(cohort, 2, 5, so);
  std::printf("      ranks:");
  for (const auto& k : ex.ranks) std::printf(" %zu (stability %.3f, min %.3f, error %.3f)", k.rank, k.average_stability, k.minimum_stability, k.relative_error);
  std::printf("\n");
  s.require("the rank chosen is 3, the largest stable one", ex.stable && ex.solution().rank == 3, std::to_string(ex.solution().rank));
  const auto extracted = somatic::columns_of(ex.solution().signatures);
  const auto m = somatic::match_signatures(extracted, truth);
  std::printf("      cosines to the truth:");
  for (double c : m.cosine) std::printf(" %.4f", c);
  std::printf("\n");
  double least = 1;
  for (double c : m.cosine) least = std::min(least, c);
  s.check("every extracted signature matches its truth (1 − cosine)", 1 - least, 0.03);
  bool distinct = m.reference.size() == 3 && m.reference[0] != m.reference[1] && m.reference[1] != m.reference[2] && m.reference[0] != m.reference[2];
  s.require("the matching is one to one", distinct);

  // ---- exposures of one tumour, refitted and bootstrapped -------------------------------------------
  {
    double rel = 0;
    std::size_t covered = 0, intervals = 0;
    for (std::size_t i = 0; i < 10; i++) {
      const auto& c = cohort[i];
      const auto b = somatic::bootstrap_exposures(c, truth, 200, 40 + i);
      double total = 0, err = 0;
      for (std::size_t j = 0; j < 3; j++) {
        total += exposures[i][j];
        err += std::fabs(b.exposures[j].estimate - exposures[i][j]);
        if (exposures[i][j] > 50) {
          intervals++;
          covered += b.exposures[j].low <= exposures[i][j] && exposures[i][j] <= b.exposures[j].high;
        }
      }
      if (total > 0) rel = std::max(rel, err / total);
    }
    s.check("refitted exposures of ten tumours: Σ|ê − e| / Σe", rel, 0.15);
    std::printf("      bootstrap 95%% intervals covering the exposure: %zu of %zu\n", covered, intervals);
    s.require("bootstrap 95% intervals cover the true exposure at least 80% of the time", intervals > 0 && static_cast<double>(covered) >= 0.8 * static_cast<double>(intervals));
  }
  return s.done();
}
