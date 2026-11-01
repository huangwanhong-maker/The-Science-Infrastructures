// Checks of the whole detection on a simulated individual: a 2 Mb genome of two contigs with
// interspersed repeats, carrying deletions, tandem and dispersed duplications, inversions, novel
// insertions, a balanced reciprocal translocation, a large loss and a large gain, a copy-neutral
// loss of heterozygosity and the small variants of a diploid, read as 30× of 2 × 150-base pairs
// from fragments of 450 ± 45 bases with a GC bias, mapped by the mapper and detected. Against the
// truth: precision and recall of every kind by Truvari's matching; the error of the breakpoints of
// precise calls; genotypes; the copy number of the bins and the dispersion of the depth at two
// copies; the loss of heterozygosity found and told from the deletion; the library estimated; the
// same calls on one and eight threads.

#include "scenario.hpp"

#include <calx/verification/check.hpp>

int main() {
  using namespace scenario;
  calx::verification::Suite s;
  calx::scheduler::ThreadPool one(1), eight(8);

  simulate::GenomeSpec spec;
  spec.length = 1'000'000;
  spec.gc = 0.42;
  spec.repeats = {8, 4, 1500, 0.03};
  spec.tandem = {6, 40, 12, 0.02};
  const std::vector<Codes> contigs{simulate::make_genome(spec, 31).sequence, simulate::make_genome(spec, 32).sequence};

  simulate::StructuralSpec ss;
  ss.deletions = 14;
  ss.tandem_duplications = 6;
  ss.inversions = 4;
  ss.insertions = 6;
  ss.dispersed_duplications = 3;
  ss.translocations = 1;
  ss.losses = 1;
  ss.gains = 1;
  ss.losses_of_heterozygosity = 1;
  ss.deletion = {50, 20000};
  ss.insertion = {50, 120};
  ss.loss = {80000, 150000};
  ss.gain = {80000, 150000};
  ss.loss_of_heterozygosity = {150000, 250000};
  const auto svs = simulate::random_structural(contigs, ss, 33);

  simulate::PairModel pm;
  pm.fragment.mean = 450;
  pm.fragment.sd = 45;
  pm.fragment.gc_bias = 0.5;
  pm.fragment.gc_optimum = 0.42;
  pm.fragment.gc_width = 0.08;
  const Run r = run(eight, contigs, svs, 30, pm, 34);
  report(r, "simulated 2 Mb");
  const auto& d = r.detection;
  const auto& cmp = r.comparison;

  s.require("the library's median template is the fragments' (450 ± 5)", std::fabs(d.library.median - 450) <= 5, std::to_string(d.library.median));
  s.check("its MAD estimates the fragments' sd (45) through 1.4826", std::fabs(d.library.sd() - 45), 4, "bases");
  // Limits a margin below the measured values (2026-11-01: recall 0.972, precision 0.946, genotypes
  // 0.914, 31 of 31 precise breakpoints exact, copy number 0.992, dispersion 1.055).
  s.check("recall over every kind (Truvari's matching), at least 0.94", 1 - cmp.overall.recall(), 0.06);
  s.check("precision over every kind, at least 0.90", 1 - cmp.overall.precision(), 0.10);
  s.check("genotype concordance of the matches, at least 0.85", 1 - cmp.overall.genotype_concordance(), 0.15);
  bool kinds = true;
  for (const char* k : {"DEL", "DUP:TANDEM", "INV", "TRA", "DUP:DISPERSED"}) kinds = kinds && cmp.by_kind.count(k) && cmp.by_kind.at(k).recall() == 1;
  s.require("every deletion, tandem duplication, inversion, translocation and dispersed duplication found", kinds);
  s.check("insertions found (recall), at least 0.66", 1 - (cmp.by_kind.count("INS") ? cmp.by_kind.at("INS").recall() : 0), 0.34);
  std::size_t exact = 0;
  for (const auto e : cmp.breakpoint_errors) exact += e == 0;
  s.check("precise breakpoints at their exact junction (microhomology aside), at least 95%", 1 - static_cast<double>(exact) / static_cast<double>(std::max<std::size_t>(1, cmp.breakpoint_errors.size())), 0.05);
  s.check("copy number concordance over bins, at least 0.98", 1 - r.copy.fraction(), 0.02);
  s.check("the dispersion of the depth at two copies, per bin, near its Poisson value 1", std::fabs(dispersion_at_two(r) - 1), 0.3);
  std::vector<std::array<std::size_t, 3>> called, truth;
  loh_intervals(r, called, truth);
  const auto o = structural::base_overlap(called, truth);
  std::printf("# copy-neutral LOH: recall %.3f precision %.3f of %zu bases\n", o.recall(), o.precision(), o.truth);
  s.check("the copy-neutral LOH found (base recall, measured 0.989)", 1 - o.recall(), 0.05);
  s.check("and nothing else called copy-neutral LOH (base precision, measured 0.990)", 1 - o.precision(), 0.05);

  // Determinism: the detection on one thread.
  const auto d1 = structural::detect(one, contigs, r.mates, r.pairs);
  bool same = d1.calls.size() == d.calls.size() && d1.segments.size() == d.segments.size() && d1.small.size() == d.small.size();
  for (std::size_t i = 0; same && i < d.calls.size(); i++)
    same = d1.calls[i].kind == d.calls[i].kind && d1.calls[i].begin == d.calls[i].begin && d1.calls[i].end == d.calls[i].end && d1.calls[i].genotype.genotype == d.calls[i].genotype.genotype &&
           d1.calls[i].qual == d.calls[i].qual;
  for (std::size_t i = 0; same && i < d.depth.bins.size(); i++) same = d1.depth.bins[i].exposure == d.depth.bins[i].exposure;
  s.require("the same calls, segments, exposures and small variants on one thread as on eight", same);
  return s.done();
}
