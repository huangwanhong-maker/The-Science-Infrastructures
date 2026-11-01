// Checks of copy number and runs of homozygosity on observations drawn from their models: Poisson
// counts of mean (c/2) e over bins of varying exposure, with c a step function (2, 1, 2, 3, 2, 0,
// 2, 4, 2), segmented and given integer copy numbers that are the truth's in every bin; the
// integer maximum of the likelihood equal to the brute-force maximum over 0 … 8; the dispersion
// Σ (n − e)²/e within χ² bounds of the number of bins; a contig at copy number 3 everywhere called
// trisomic, one at 2 not. Runs of homozygosity: on sites with clean genotype likelihoods, bcftools
// roh's model enters a run of k homozygous-alternate sites between heterozygous ones exactly when
// k ln(1/f) exceeds −ln P(AZ|HW) − ln P(HW|AZ) (39 sites at the defaults), and a run found begins
// and ends at its first and last site.

#include <genomics/structural/allelic.hpp>
#include <genomics/structural/copy_number.hpp>
#include <genomics/structural/depth.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <array>
#include <cmath>
#include <string>
#include <vector>

using namespace genomics::structural;

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(17, 0x434e4331, 0, 0);

  {
  // ---- a contig of 3000 bins with copy-number steps; a second contig trisomic ---------------------
  const std::vector<std::pair<std::size_t, int>> steps{{0, 2}, {400, 1}, {700, 2}, {1000, 3}, {1300, 2}, {1500, 0}, {1600, 2}, {2000, 4}, {2200, 2}};
  const auto truth = [&](std::size_t b) {
    int c = 2;
    for (const auto& [at, k] : steps) if (b >= at) c = k;
    return c;
  };
  Depth d;
  const std::size_t W = 1000;
  d.first = {0};
  for (std::uint32_t contig = 0; contig < 3; contig++) {
    const std::size_t n = contig == 0 ? 3000 : 800;
    for (std::size_t b = 0; b < n; b++) {
      DepthBin bin;
      bin.contig = contig;
      bin.begin = b * W;
      bin.end = (b + 1) * W;
      bin.exposure = 150 * (0.7 + 0.6 * r.uniform());  // GC and mappability make exposures uneven
      const int c = contig == 0 ? truth(b) : contig == 1 ? 3 : 2;
      bin.count = static_cast<double>(r.poisson((c == 0 ? 0.02 : c / 2.0) * bin.exposure));
      d.bins.push_back(bin);
    }
    d.first.push_back(d.bins.size());
  }
  const auto segments = segment_copy_number(d);
  std::size_t agree = 0, bins = 0;
  for (std::size_t b = 0; b < 3000; b++) {
    bins++;
    agree += copies_at(segments, 0, b * W + W / 2) == truth(b);
  }
  std::size_t on_first = 0;
  for (const auto& g : segments) on_first += g.contig == 0;
  s.require("the steps are found: the integer copy number of every bin is the truth's", agree == bins, std::to_string(agree) + " of " + std::to_string(bins));
  s.require("as nine segments of the first contig", on_first == steps.size(), std::to_string(on_first));
  bool boundaries = true;
  for (std::size_t k = 1; k < steps.size(); k++) {
    bool found = false;
    for (const auto& g : segments) found = found || (g.contig == 0 && g.begin == steps[k].first * W);
    boundaries = boundaries && found;
  }
  s.require("each change of copy number at its bin exactly", boundaries);
  // The integer maximum against brute force.
  bool brute = true;
  for (const auto& g : segments) {
    double best = -INFINITY;
    int arg = -1;
    for (int c = 0; c <= 8; c++) {
      const double mu = c == 0 ? 0.02 : c / 2.0, l = g.count * std::log(mu) - mu * g.exposure;
      if (l > best) { best = l; arg = c; }
    }
    brute = brute && arg == g.copies;
  }
  s.require("each segment's copy number is the likelihood's maximum over 0 … 8 (brute force)", brute);
  // Dispersion at copy number 2 (the third contig): E Σ (n − e)²/e = number of bins under Poisson.
  double disp = 0;
  for (std::size_t b = d.first[2]; b < d.first[3]; b++) disp += (d.bins[b].count - d.bins[b].exposure) * (d.bins[b].count - d.bins[b].exposure) / d.bins[b].exposure;
  const double nb = 800;
  s.check("the dispersion of Poisson counts is the number of bins (within 4 sd of χ²)", std::fabs(disp - nb) / std::sqrt(2 * nb), 4);
  const auto an = find_aneuploidy(d, segments);
  s.require("the contig at three copies everywhere is trisomic, the others are not", an.size() == 1 && an[0].contig == 1 && an[0].copies == 3, std::to_string(an.size()));
  s.require("the first contig's copy number 0 segment is a loss to nothing", copies_at(segments, 0, 1550 * W) == 0);
}

  // ---- runs of homozygosity: the threshold of bcftools roh's model ------------------------------
  {
    const RohOptions o;
    const double threshold = (-std::log(o.hw_to_az) - std::log(o.az_to_hw)) / std::log(1 / o.allele_frequency);
    const std::size_t kmin = static_cast<std::size_t>(std::floor(threshold)) + 1;
    s.require("the closed form gives 39 sites at bcftools' defaults", kmin == 39, std::to_string(kmin));
    // Sites: 200 heterozygous, k homozygous alternate, 200 heterozygous; clean likelihoods (PL 0 and 200).
    const auto run_found = [&](std::size_t k, std::size_t& first, std::size_t& last) {
      std::vector<AlleleSite> sites;
      std::size_t pos = 1000;
      const auto add = [&](int g) {
        AlleleSite x;
        x.position = pos;
        pos += 1000;
        x.genotype = g;
        x.log10_likelihood = g == 1 ? std::array<double, 3>{-20, 0, -20} : std::array<double, 3>{-20, -20, 0};
        x.reference = g == 1 ? 15 : 0;
        x.alternate = g == 1 ? 15 : 30;
        sites.push_back(x);
      };
      for (int i = 0; i < 200; i++) add(1);
      for (std::size_t i = 0; i < k; i++) add(2);
      for (int i = 0; i < 200; i++) add(1);
      RohOptions oo;
      oo.min_sites = 1;
      const auto runs = runs_of_homozygosity(sites, oo);
      if (runs.size() != 1) return false;
      first = runs[0].first;
      last = runs[0].last;
      return true;
    };
    std::size_t f = 0, l = 0;
    const bool below = run_found(kmin - 1, f, l), above = run_found(kmin, f, l);
    s.require("a run of 38 homozygous sites is not entered, one of 39 is", !below && above);
    std::size_t f2 = 0, l2 = 0;
    const bool long_run = run_found(300, f2, l2);
    s.require("a run found begins and ends at its first and last homozygous site", long_run && f2 == 200 && l2 == 499, std::to_string(f2) + "–" + std::to_string(l2));
  }
  return s.done();
}
