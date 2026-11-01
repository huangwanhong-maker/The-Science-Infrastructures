// Checks of the tumour simulator (tumour.hpp): clone trees refused against the sum rule; the cell
// fractions; every mutation's expected allele fraction from its CCF, multiplicity and copy number;
// early mutations at every copy of their allele and none on a lost one; the truth's SBS-96
// catalogue following the signatures' mixture; the depth of tumour reads following copy number;
// the fraction of reads carrying the mutations against the expected one; a shard of reads equal to
// the same reads of the whole sample; cohorts of catalogues with the Poisson means of the model.

#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/mutate.hpp>
#include <genomics/simulate/tumour.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

using namespace genomics;
using sequence::Codes;

namespace {

double cosine(const std::vector<double>& a, const std::vector<double>& b) {
  double ab = 0, aa = 0, bb = 0;
  for (std::size_t i = 0; i < a.size(); i++) { ab += a[i] * b[i]; aa += a[i] * a[i]; bb += b[i] * b[i]; }
  return ab / std::sqrt(aa * bb);
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const Codes ref = simulate::random_sequence(400000, 0.45, 11);
  const std::vector<Codes> contigs{ref};
  simulate::MutationRates rates;
  rates.snv = 1.5e-3;
  rates.indel = 1e-4;
  const simulate::Diploid normal = simulate::mutate_diploid(ref, rates, 12);
  const auto sig = simulate::synthetic_signatures();

  // ---- refusals ----------------------------------------------------------------------------------------
  {
    simulate::TumourSpec bad;
    bad.clones = {{-1, 1.0}, {0, 0.7}, {0, 0.5}};
    bool refused = false;
    try { (void)simulate::make_tumour(contigs, {normal}, bad, 1); } catch (const std::invalid_argument&) { refused = true; }
    s.require("children of CCF 0.7 + 0.5 under a trunk of 1 are refused (the sum rule)", refused);
  }

  simulate::TumourSpec spec;
  spec.purity = 0.6;
  spec.clones = {{-1, 1.0, 3000, 60, {1, 0, 0}}, {0, 0.6, 1500, 30, {0, 1, 0}}, {1, 0.25, 1000, 20, {0, 0, 1}}};
  spec.signatures = {sig[0].profile, sig[1].profile, sig[3].profile};
  spec.copy_number = {{0, 0, 80000, {1, 1}}, {0, 80000, 160000, {2, 1}}, {0, 160000, 220000, {1, 0}}, {0, 220000, 300000, {2, 0}}, {0, 300000, 400000, {3, 1}}};
  spec.early = 0.5;
  const simulate::Tumour t = simulate::make_tumour(contigs, {normal}, spec, 13);

  {
    double total = 0;
    for (double f : t.fraction) total += f;
    s.check("the clones' own fractions sum to 1", total - 1, 1e-12);
    s.check("the fractions: 0.4, 0.35, 0.25", std::fabs(t.fraction[0] - 0.4) + std::fabs(t.fraction[1] - 0.35) + std::fabs(t.fraction[2] - 0.25), 1e-12);
    s.require("lineages: clone 2 descends from 1 and 0", t.lineage[2] == 0b111 && t.lineage[1] == 0b011 && t.lineage[0] == 0b001);
  }
  {
    double worst = 0;
    bool lost = true, early = true;
    std::size_t n_early = 0;
    for (const auto& v : t.truth) {
      const auto& g = t.segment_at(v.contig, v.position);
      const double f = t.purity * v.ccf * v.multiplicity / (2 * (1 - t.purity) + t.purity * g.total());
      worst = std::max(worst, std::fabs(f - v.vaf));
      lost = lost && g.copies[static_cast<std::size_t>(v.haplotype)] >= 1;
      if (v.early) {
        n_early++;
        early = early && v.clone == 0 && v.multiplicity == g.copies[static_cast<std::size_t>(v.haplotype)];
      } else {
        early = early && v.multiplicity == 1;
      }
    }
    s.check("every mutation's VAF is ρ φ m / (2(1 − ρ) + ρ N)", worst, 1e-15);
    s.require("no mutation lies on a lost haplotype", lost);
    s.require("early mutations are the trunk's, at every copy of their allele; the others at one", early && n_early > 0, std::to_string(n_early) + " early");
  }
  {
    // The truth catalogue of each clone against its signature.
    for (int j = 0; j < 3; j++) {
      std::vector<double> cat(96, 0);
      for (const auto& v : t.truth)
        if (v.clone == j && v.sbs96 >= 0) cat[static_cast<std::size_t>(v.sbs96)]++;
      // Each class drawn has the reference trinucleotide of its class.
      s.check("clone " + std::to_string(j) + "'s substitutions follow its signature (1 − cosine)", 1 - cosine(cat, spec.signatures[static_cast<std::size_t>(j)]), 0.02);
    }
    bool context = true;
    for (const auto& v : t.truth)
      if (v.sbs96 >= 0) context = context && sequence::sbs96_at(ref, v.position, sequence::dna_code(v.alt[0])) == v.sbs96 && v.ref[0] == sequence::dna_letter(ref[v.position]);
    s.require("every substitution's class is that of its reference trinucleotide", context);
  }

  // ---- reads -------------------------------------------------------------------------------------------
  simulate::ShortReadModel model;
  const std::size_t count = simulate::reads_for_depth(t, 40, model.length);
  const auto reads = simulate::tumour_reads(t, count, model, 14);
  {
    const auto part = simulate::tumour_reads(t, 300, model, 14, 1000);
    bool same = part.size() == 300;
    for (std::size_t i = 0; i < part.size() && same; i++) same = part[i].bases == reads[1000 + i].bases && part[i].begin == reads[1000 + i].begin && part[i].name == reads[1000 + i].name;
    s.require("reads 1000–1299 alone equal the same reads of the whole sample", same);
  }
  {
    // Depth by segment: read starts per base against 2(1 − ρ) + ρ N, normalised.
    double worst = 0, mean_dna = 0;
    for (const auto& g : t.copy_number) mean_dna += (2 * (1 - t.purity) + t.purity * g.total()) * static_cast<double>(g.end - g.begin);
    mean_dna /= 400000;
    for (const auto& g : t.copy_number) {
      std::size_t n = 0;
      for (const auto& r : reads) n += r.begin >= g.begin && r.begin < g.end;
      const double observed = static_cast<double>(n) / static_cast<double>(g.end - g.begin) / (static_cast<double>(reads.size()) / 400000);
      const double expected = (2 * (1 - t.purity) + t.purity * g.total()) / mean_dna;
      worst = std::max(worst, std::fabs(observed / expected - 1));
    }
    s.check("tumour depth in each segment follows 2(1 − ρ) + ρ N (relative)", worst, 0.03);
  }
  {
    // Allele fractions: reads without errors read straight at substitutions away from indels.
    std::vector<std::size_t> alt(t.truth.size(), 0), depth(t.truth.size(), 0);
    std::vector<char> clean(t.truth.size(), 0);
    const auto near_indel = [&](std::size_t p) {
      for (const auto& v : normal.truth)
        if (!v.snv() && v.position + 200 > p && v.position < p + 200) return true;
      for (const auto& v : t.truth)
        if (!v.snv() && v.position + 200 > p && v.position < p + 200) return true;
      return false;
    };
    for (std::size_t i = 0; i < t.truth.size(); i++) clean[i] = t.truth[i].snv() && !near_indel(t.truth[i].position);
    for (const auto& r : reads) {
      if (r.errors) continue;
      const Codes fwd = r.reverse ? sequence::reverse_complement(r.bases) : r.bases;
      if (r.end - r.begin != fwd.size()) continue;
      auto it = std::lower_bound(t.truth.begin(), t.truth.end(), r.begin, [](const simulate::SomaticVariant& v, std::size_t x) { return v.position < x; });
      for (; it != t.truth.end() && it->position < r.end; ++it) {
        const std::size_t i = static_cast<std::size_t>(it - t.truth.begin());
        if (!clean[i]) continue;
        depth[i]++;
        alt[i] += fwd[it->position - r.begin] == sequence::dna_code(it->alt[0]);
      }
    }
    double obs = 0, expect = 0, var = 0;
    for (std::size_t i = 0; i < t.truth.size(); i++) {
      if (!clean[i] || depth[i] == 0) continue;
      obs += static_cast<double>(alt[i]);
      expect += t.truth[i].vaf * static_cast<double>(depth[i]);
      var += t.truth[i].vaf * (1 - t.truth[i].vaf) * static_cast<double>(depth[i]);
    }
    std::printf("      reads carrying the substitutions: %.0f observed, %.1f expected (sd %.1f)\n", obs, expect, std::sqrt(var));
    s.check("reads carrying the substitutions: observed − expected in standard deviations", (obs - expect) / std::sqrt(var), 4);
  }

  // ---- cohorts of catalogues -----------------------------------------------------------------------
  {
    const std::vector<std::vector<double>> sigs{sig[1].profile, sig[2].profile};
    std::vector<std::vector<double>> e(400, std::vector<double>{100, 50});
    const auto cats = simulate::simulate_catalogues(sigs, e, 15);
    std::vector<double> mean(96, 0);
    double total = 0, sq = 0;
    for (const auto& c : cats) {
      double n = 0;
      for (std::size_t k = 0; k < 96; k++) { mean[k] += c[k] / 400; n += c[k]; }
      total += n / 400;
      sq += n * n / 400;
    }
    std::vector<double> want(96);
    for (std::size_t k = 0; k < 96; k++) want[k] = 100 * sigs[0][k] + 50 * sigs[1][k];
    s.check("cohort: the mean catalogue is Σ e_j s_j (1 − cosine)", 1 - cosine(mean, want), 1e-3);
    s.check("cohort: mutations per tumour have Poisson mean and variance 150", std::fabs(total - 150) / 150 + std::fabs((sq - total * total) - 150) / 150 / 5, 0.1);
  }
  return s.done();
}
