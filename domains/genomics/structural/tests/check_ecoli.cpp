// Checks of the detection on a real reference: Escherichia coli K-12 MG1655 (NC_000913.3, from the
// data library's cache; skipped when it is absent), whose seven rRNA operons and insertion
// sequences make mapping honest. A simulated diploid of it (labelled so: E. coli is haploid; the
// individual is a construct for the checks) carries deletions, tandem and dispersed duplications,
// inversions, novel insertions, a large loss, a large gain and a copy-neutral loss of
// heterozygosity, placed at random (some fall in repeats, as real ones do), read at 30× of
// 2 × 150 bases and detected. Precision and recall by kind against thresholds measured here;
// copy-number concordance; the loss of heterozygosity found and called copy-neutral.

#include "scenario.hpp"

#include <calx/verification/check.hpp>

#include <filesystem>

int main() {
  using namespace scenario;
  calx::verification::Suite s;
  calx::scheduler::ThreadPool pool(8);
  const auto file = std::filesystem::path(cache()) / "ncbi-genomes" / "NC_000913.3.gb";
  if (!std::filesystem::exists(file)) {
    s.skip("Escherichia coli NC_000913.3 is not in the data cache (data/cache/ncbi-genomes)");
    return s.done();
  }
  const std::vector<Codes> contigs{sequence::encode_dna(sequence::read_genbank(read_file(file.string())).at(0).sequence)};

  simulate::StructuralSpec ss;
  ss.deletions = 16;
  ss.tandem_duplications = 6;
  ss.inversions = 5;
  ss.insertions = 6;
  ss.dispersed_duplications = 3;
  ss.losses = 1;
  ss.gains = 1;
  ss.losses_of_heterozygosity = 1;
  ss.deletion = {50, 20000};
  ss.insertion = {50, 120};
  ss.loss = {100000, 200000};
  ss.gain = {100000, 200000};
  ss.loss_of_heterozygosity = {200000, 400000};
  const auto svs = simulate::random_structural(contigs, ss, 41);
  simulate::PairModel pm;
  pm.fragment.gc_bias = 0.5;
  pm.fragment.gc_optimum = 0.50;
  pm.fragment.gc_width = 0.08;
  const Run r = run(pool, contigs, svs, 30, pm, 42);
  report(r, "E. coli K-12 (simulated diploid)");

  const auto& cmp = r.comparison;
  // Limits a margin below the measured values (2026-11-01: recall 0.947, precision 0.923, copy
  // number 0.981; the false calls are junctions between the rRNA operons' copies).
  s.check("recall over every kind (Truvari's matching), at least 0.90", 1 - cmp.overall.recall(), 0.10);
  s.check("precision over every kind, at least 0.85", 1 - cmp.overall.precision(), 0.15);
  s.check("every deletion found", 1 - (cmp.by_kind.count("DEL") ? cmp.by_kind.at("DEL").recall() : 0), 0.0);
  s.check("copy-number concordance over mappable bins, at least 0.96", 1 - r.copy.fraction(), 0.04);
  std::vector<std::array<std::size_t, 3>> called, truth;
  loh_intervals(r, called, truth);
  const auto o = structural::base_overlap(called, truth);
  std::printf("# copy-neutral LOH: recall %.3f precision %.3f of %zu bases\n", o.recall(), o.precision(), o.truth);
  s.check("the copy-neutral LOH found (base recall, measured 0.985)", 1 - o.recall(), 0.05);
  return s.done();
}
