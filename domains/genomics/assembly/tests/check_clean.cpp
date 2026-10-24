// Checks of assembly/clean.hpp: a tip and a bubble planted in the graph of a genome are found,
// and removed, and nothing else; on reads with errors of a 200 kb genome, the raw graph has
// tips and bubbles, the cleaned one has none, what was removed is almost all erroneous
// k-mers and what remains almost none, the genome's k-mers are kept, and the cleaned graph is
// the graph built from the k-mers that remain; every step is recorded, with its merges; a
// diploid genome's heterozygous SNVs make bubbles, one for each cluster of SNVs closer than k,
// crossed bubbles included, which the cleaning pops, leaving one contig, a mosaic of the
// haplotypes.

#include "truth.hpp"

#include <genomics/assembly/clean.hpp>
#include <genomics/assembly/debruijn.hpp>
#include <genomics/kmers/counter.hpp>
#include <genomics/kmers/spectrum.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/mutate.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <string>
#include <vector>

using namespace genomics;
using assembly::Graph;
using sequence::Codes;
using sequence::Kmer;

namespace {

// Whether a sequence has a k-mer the genome lacks.
bool erroneous(const truth::GenomeIndex& index, const Codes& s) { return index.kmer_fraction(s) < 1; }

// The canonical k-mers of a graph with their counts, sorted: the input of build_graph.
std::vector<kmers::KmerCount> graph_kmers(const Graph& g, const kmers::KmerCounter& counts) {
  std::vector<kmers::KmerCount> v;
  for (const auto& u : g.unitigs) sequence::for_each_kmer(u.sequence, g.k, [&](std::size_t, Kmer f, Kmer r) { v.push_back({std::min(f, r), counts.count(f)}); });
  std::sort(v.begin(), v.end(), [](const auto& x, const auto& y) { return x.kmer < y.kmer; });
  return v;
}

bool same_graph(const Graph& a, const Graph& b) {
  if (a.size() != b.size()) return false;
  for (std::size_t u = 0; u < a.size(); u++)
    if (a.unitigs[u].sequence != b.unitigs[u].sequence || a.unitigs[u].count_sum != b.unitigs[u].count_sum) return false;
  for (std::size_t o = 0; o < a.next.size(); o++) {
    auto x = a.next[o].to, y = b.next[o].to;
    if (a.next[o].size != b.next[o].size) return false;
    std::sort(x.begin(), x.begin() + a.next[o].size);
    std::sort(y.begin(), y.begin() + b.next[o].size);
    if (x != y) return false;
  }
  return true;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::scheduler::ThreadPool pool(8);
  const int k = 31;

  // A planted tip and bubble on a genome's own k-mers (counts 20).
  {
    const Codes genome = simulate::random_sequence(20'000, 0.5, 21);
    std::vector<Codes> pieces = {genome};
    // The tip: 20 bases leaving the genome after position 5000 (k − 1 bases shared), count 2.
    Codes tip(genome.begin() + 5000 - (k - 1), genome.begin() + 5000);
    Codes tail = simulate::random_sequence(20, 0.5, 22);
    tail[0] = static_cast<sequence::Code>((genome[5000] + 1) % 4);  // leaves the genome at its first base
    tip.insert(tip.end(), tail.begin(), tail.end());
    // The bubble: the genome's bases 9970 to 10031 with base 10000 changed, count 3.
    Codes arm(genome.begin() + 10000 - (k - 1), genome.begin() + 10000 + k);
    arm[k - 1] = static_cast<sequence::Code>((arm[k - 1] + 1) % 4);
    kmers::KmerCounter counts;
    for (int i = 0; i < 20; i++) counts.add(pieces, pool);
    for (int i = 0; i < 2; i++) counts.add({tip}, pool);
    for (int i = 0; i < 3; i++) counts.add({arm}, pool);
    const Graph g = assembly::build_graph(counts.sorted(1, pool), k, pool);
    const auto tips = assembly::find_tips(g), bubbles = assembly::find_bubbles(g);
    s.require("clean: the planted tip is the one tip found, 20 k-mers at coverage 2 against 20",
              tips.size() == 1 && tips[0].length == static_cast<std::size_t>(k - 1 + 20) && tips[0].coverage == 2 && tips[0].reference == 20,
              std::to_string(tips.size()) + " tips");
    s.require("clean: the planted bubble is the one bubble found, its arm of k k-mers at coverage 3 against the genome's 20",
              bubbles.size() == 1 && bubbles[0].length == static_cast<std::size_t>(2 * k - 1) && bubbles[0].coverage == 3 && bubbles[0].reference == 20,
              std::to_string(bubbles.size()) + " bubbles");
    const auto c = assembly::clean(g);
    s.require("clean: after the cleaning the graph is the genome, one unitig", c.graph.size() == 1 &&
              (c.graph.unitigs[0].sequence == genome || c.graph.unitigs[0].sequence == sequence::reverse_complement(genome)));
    s.require("clean: two steps recorded, the tip and then the bubble, each with its merges",
              c.steps.size() == 2 && c.steps[0].kind == assembly::StepKind::tips && c.steps[1].kind == assembly::StepKind::bubbles &&
                  c.steps[0].merged_from.size() == c.steps[0].unitigs_after && c.steps[1].unitigs_after == 1);
  }

  // Reads with errors of a 200 kb genome with repeats, at 30×; errors twice the default.
  {
    simulate::GenomeSpec spec;
    spec.length = 200'000;
    spec.repeats = {5, 3, 600, 0.0};
    const auto genome = simulate::make_genome(spec, 30);
    simulate::ShortReadModel model;
    model.error_first = 0.002;
    model.error_last = 0.02;
    const auto reads = truth::reads(simulate::Source::haploid({genome.sequence}), 30 * spec.length / 150, model, 31, pool);
    const auto counts = kmers::count_kmers(reads, pool);
    // The spectrum's threshold (4 here) leaves almost no error to clean at 30×; a threshold of
    // 2 keeps the errors seen twice or more, and the structures they make.
    const auto fit = kmers::fit_spectrum(kmers::spectrum(counts, pool));
    const std::uint32_t threshold = 2;
    const Graph g = assembly::build_graph(counts.sorted(threshold, pool), k, pool);
    const truth::GenomeIndex index(genome.sequence, k);
    const std::size_t tips0 = assembly::find_tips(g).size(), bubbles0 = assembly::find_bubbles(g).size();
    assembly::CleaningOptions o;
    o.keep_graphs = true;
    const auto c = assembly::clean(g, o);
    s.require("clean: the raw graph of reads with errors has tips and bubbles", tips0 > 0 && bubbles0 > 0,
              std::to_string(g.size()) + " unitigs, " + std::to_string(tips0) + " tips, " + std::to_string(bubbles0) + " bubbles, threshold 2 (the spectrum's " + std::to_string(fit.threshold) + ")");
    s.require("clean: the cleaned graph has no tip and no bubble", assembly::find_tips(c.graph).empty() && assembly::find_bubbles(c.graph).empty(),
              std::to_string(c.graph.size()) + " unitigs after " + std::to_string(c.steps.size()) + " steps");
    // What was removed, and what remains.
    std::size_t removed = 0, removed_errors = 0;
    const Graph* before = &g;
    for (const auto& step : c.steps) {
      for (const auto& r : step.removed) {
        removed++;
        removed_errors += erroneous(index, before->unitigs[r.unitig].sequence);
      }
      before = &step.after;
    }
    s.check("clean: removed unitigs that hold only genome k-mers (fraction)", 1 - static_cast<double>(removed_errors) / static_cast<double>(removed), 0.01);
    std::size_t kept = 0, kept_errors = 0, genome_kept = 0;
    for (const auto& u : c.graph.unitigs)
      sequence::for_each_kmer(u.sequence, k, [&](std::size_t, Kmer f, Kmer r) {
        kept++;
        kept_errors += !index.has(std::min(f, r));
      });
    std::size_t genome_kmers = 0;
    const auto in_graph = graph_kmers(c.graph, counts);
    kmers::KmerCounter genome_counts = kmers::count_kmers({genome.sequence}, pool);
    genome_counts.for_each([&](Kmer x, std::uint32_t) {
      genome_kmers++;
      genome_kept += std::binary_search(in_graph.begin(), in_graph.end(), kmers::KmerCount{x, 0}, [](const auto& a, const auto& b2) { return a.kmer < b2.kmer; });
    });
    s.check("clean: erroneous k-mers left in the cleaned graph (fraction)", static_cast<double>(kept_errors) / static_cast<double>(kept), 1e-3);
    s.check("clean: genome k-mers missing from the cleaned graph (fraction)", 1 - static_cast<double>(genome_kept) / static_cast<double>(genome_kmers), 1e-3);
    s.require("clean: the cleaned graph is the graph built from the k-mers it keeps", same_graph(c.graph, assembly::build_graph(in_graph, k, pool)));
    bool merges = true;
    before = &g;
    for (const auto& step : c.steps) {
      std::size_t parts = 0;
      for (const auto& m : step.merged_from) parts += m.size();
      merges = merges && step.merged_from.size() == step.after.size() && parts + step.removed.size() == before->size() && step.unitigs_before == before->size();
      before = &step.after;
    }
    s.require("clean: each step's merges account for every unitig before it, kept or removed", merges);
    s.require("clean: tips, bubbles and low-coverage unitigs were removed", c.removed(assembly::StepKind::tips) > 0 && c.removed(assembly::StepKind::bubbles) > 0,
              std::to_string(c.removed(assembly::StepKind::tips)) + " tips, " + std::to_string(c.removed(assembly::StepKind::bubbles)) + " bubble arms, " +
                  std::to_string(c.removed(assembly::StepKind::low_coverage)) + " low-coverage");
  }

  // A diploid genome: heterozygous SNVs only, error-free reads at 40× over both haplotypes.
  {
    const std::size_t G = 200'000;
    const Codes ref = simulate::random_sequence(G, 0.5, 41);
    simulate::MutationRates rates;
    rates.snv = 1e-3;
    rates.indel = 0;
    rates.heterozygous = 1;
    // The SNVs within 2k of an end are dropped, and the haplotypes rebuilt from the others: an
    // SNV near an end makes a fork of two dead ends, not a bubble.
    simulate::Diploid d = simulate::mutate_diploid(ref, rates, 42);
    std::erase_if(d.truth, [&](const simulate::Variant& v) { return v.position < 2 * static_cast<std::size_t>(k) || v.position + 2 * static_cast<std::size_t>(k) > G; });
    for (int h = 0; h < 2; h++) {
      d.haplotype[h] = ref;
      d.liftover[h] = simulate::Liftover::identity(G);
      for (const auto& v : d.truth)
        if (v.on[h]) d.haplotype[h][v.position] = sequence::dna_code(v.alt[0]);
    }
    const auto reads = truth::reads(simulate::Source::diploid({d}), 40 * G / 150, {150, 0, 0, 0}, 43, pool);
    const auto counts = kmers::count_kmers(reads, pool);
    const Graph g = assembly::build_graph(counts.sorted(2, pool), k, pool);
    // Clusters: SNVs closer than k make one bubble, one arm to pop. Two clusters exactly k apart
    // share no k-mer but a (k − 1)-mer, where the k-mers of their alleles cross: a crossed bubble,
    // two arms to pop, one for each cluster (this seed has one).
    std::size_t clusters = 0, last = 0, crosses = 0;
    bool open = false;
    for (const auto& v : d.truth) {
      if (open && v.position - last == static_cast<std::size_t>(k)) crosses++;
      if (!open || v.position - last >= static_cast<std::size_t>(k)) clusters++;
      open = true;
      last = v.position;
    }
    // Arms of clusters span 2k − 1 bases and the cluster; long enough arms are allowed.
    assembly::CleaningOptions o;
    o.bubble_length = 5;
    const auto bubbles = assembly::find_bubbles(g, o);
    const auto c = assembly::clean(g, o);
    s.require("clean: a diploid genome's heterozygous SNVs make one arm to pop for each cluster of SNVs closer than k, crossed bubbles included", bubbles.size() == clusters,
              std::to_string(d.truth.size()) + " SNVs, " + std::to_string(clusters) + " clusters, " + std::to_string(crosses) + " crosses, " + std::to_string(bubbles.size()) + " arms to pop");
    // The contig is a mosaic of the haplotypes, a few bases short at each end, where fewer than
    // two reads begin.
    const auto haplotypes = kmers::count_kmers({d.haplotype[0], d.haplotype[1]}, pool);
    bool mosaic = c.graph.size() == 1;
    if (mosaic) sequence::for_each_kmer(c.graph.unitigs[0].sequence, k, [&](std::size_t, Kmer f, Kmer) { mosaic = mosaic && haplotypes.count(f) > 0; });
    s.require("clean: the cleaning pops them all, leaving one contig, the genome to within 50 bases at each end, every k-mer from a haplotype",
              c.removed(assembly::StepKind::bubbles) == clusters && mosaic && c.graph.unitigs[0].length() + 100 >= G,
              std::to_string(c.graph.size()) + " unitigs, the longest " + std::to_string(c.graph.size() ? c.graph.unitigs[0].length() : 0));
  }
  return s.done();
}
