// Checks of the gene finder: the Markov chains with phases recover the chain they were trained
// on; on a simulated genome whose genes are drawn from a known codon usage and planted on both
// strands between non-coding sequence of another composition, the finder trained on nothing
// but the genome finds the genes' stops with high sensitivity and specificity and reports how
// many starts it gets right; and, when the data cache has Escherichia coli K-12 MG1655
// (NC_000913.3), its sensitivity and specificity against the CDSs of NCBI's annotation are
// measured and reported.

#include <genomics/annotation/genes.hpp>
#include <genomics/sequence/genbank.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::annotation;

namespace {

std::string fixed(double v, int digits = 3) {
  char b[32];
  std::snprintf(b, sizeof b, "%.*f", digits, v);
  return b;
}

std::string read_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::ostringstream s;
  s << f.rdbuf();
  return s.str();
}

struct TrueGene {
  std::size_t begin = 0, end = 0;
  bool reverse = false;
};

// A genome of planted genes: codons drawn from a codon usage (Dirichlet weights over the 61
// sense codons), starts ATG, GTG and TTG, stops TAA, TGA and TAG, genes of 100 codons and more on
// both strands, separated by non-coding sequence from a first-order chain of its own.
struct Simulated {
  Codes genome;
  std::vector<TrueGene> genes;
};

Simulated simulate_genome(std::size_t count, std::uint64_t seed) {
  calx::random::Stream r(seed, 0x67656e65, 0, 0);
  const auto& code = sequence::genetic_code(11);
  std::vector<int> sense;
  std::vector<double> usage;
  double total = 0;
  for (int i = 0; i < 64; i++)
    if (code.amino_acids[static_cast<std::size_t>(i)] != '*') {
      sense.push_back(i);
      usage.push_back(r.gamma(0.6));
      total += usage.back();
    }
  for (auto& u : usage) u /= total;
  std::array<double, 16> chain{};
  for (int a = 0; a < 4; a++) {
    double t = 0;
    for (int b = 0; b < 4; b++) t += chain[4 * a + b] = r.gamma(2.0);
    for (int b = 0; b < 4; b++) chain[4 * a + b] /= t;
  }
  const auto draw = [&](const double* p, std::size_t n) {
    double u = r.uniform();
    for (std::size_t i = 0; i < n; i++) {
      if (u < p[i]) return i;
      u -= p[i];
    }
    return n - 1;
  };
  Simulated s;
  Code last = sequence::A;
  const auto noncoding = [&](std::size_t length) {
    for (std::size_t i = 0; i < length; i++) s.genome.push_back(last = static_cast<Code>(draw(&chain[4 * last], 4)));
  };
  const double start_p[3] = {0.85, 0.11, 0.04}, stop_p[3] = {0.6, 0.3, 0.1};
  const char* starts[3] = {"ATG", "GTG", "TTG"};
  const char* stops[3] = {"TAA", "TGA", "TAG"};
  noncoding(200);
  for (std::size_t g = 0; g < count; g++) {
    const std::size_t codons = 100 + static_cast<std::size_t>(std::min(1400.0, r.exponential() * 220));
    Codes cds = sequence::encode_dna(starts[draw(start_p, 3)]);
    for (std::size_t c = 1; c + 1 < codons; c++) {
      const auto triple = sequence::codon_of(sense[draw(usage.data(), usage.size())]);
      cds.insert(cds.end(), triple.begin(), triple.end());
    }
    const Codes stop = sequence::encode_dna(stops[draw(stop_p, 3)]);
    cds.insert(cds.end(), stop.begin(), stop.end());
    const bool reverse = r.bernoulli(0.5);
    if (reverse) cds = sequence::reverse_complement(cds);
    TrueGene t;
    t.begin = s.genome.size();
    t.end = t.begin + cds.size();
    t.reverse = reverse;
    s.genome.insert(s.genome.end(), cds.begin(), cds.end());
    s.genes.push_back(t);
    last = s.genome.back();
    noncoding(20 + r.below(300));
  }
  return s;
}

// The stop of a gene: its end on the forward strand, its begin on the reverse.
std::size_t stop_of(std::size_t begin, std::size_t end, bool reverse) { return reverse ? begin : end; }

struct Evaluation {
  double sensitivity = 0, specificity = 0, starts = 0;
  std::size_t matched = 0;
};

Evaluation evaluate(const std::vector<TrueGene>& truth, const std::vector<Gene>& found) {
  std::map<std::pair<std::size_t, bool>, std::size_t> start_of;  // (stop, strand) → start
  for (const auto& t : truth) start_of[{stop_of(t.begin, t.end, t.reverse), t.reverse}] = t.reverse ? t.end : t.begin;
  Evaluation e;
  std::size_t starts = 0;
  for (const auto& g : found) {
    const auto it = start_of.find({stop_of(g.begin, g.end, g.reverse), g.reverse});
    if (it == start_of.end()) continue;
    e.matched++;
    starts += it->second == (g.reverse ? g.end : g.begin);
  }
  e.sensitivity = static_cast<double>(e.matched) / static_cast<double>(truth.size());
  e.specificity = found.empty() ? 0 : static_cast<double>(e.matched) / static_cast<double>(found.size());
  e.starts = e.matched ? static_cast<double>(starts) / static_cast<double>(e.matched) : 0;
  return e;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- the chains ------------------------------------------------------------------------------------------
  {
    // A periodic chain of order 2 learnt from a long sequence drawn from it.
    calx::random::Stream r(3, 0x63686e31, 0, 0);
    std::vector<double> truth(3 * 16 * 4);
    for (std::size_t i = 0; i < truth.size(); i += 4) {
      double t = 0;
      for (int b = 0; b < 4; b++) t += truth[i + static_cast<std::size_t>(b)] = r.gamma(1.0);
      for (int b = 0; b < 4; b++) truth[i + static_cast<std::size_t>(b)] /= t;
    }
    Codes x;
    std::uint32_t ctx = 0;
    for (std::size_t t = 0; t < 600000; t++) {
      const std::size_t phase = t % 3, row = ((phase << 4) + ctx) * 4;
      double u = r.uniform();
      Code b = 3;
      for (Code c = 0; c < 4; c++) {
        if (u < truth[row + c]) { b = c; break; }
        u -= truth[row + c];
      }
      x.push_back(b);
      ctx = ((ctx << 2) | b) & 15;
    }
    PeriodicChain chain(2, 3);
    chain.add(x, 2, x.size(), 2, 0);
    chain.finish(1);
    double worst = 0;
    for (std::size_t i = 0; i < truth.size(); i++) worst = std::max(worst, std::fabs(std::exp2(chain.table()[i]) - truth[i]));
    s.check("a three-periodic chain of order 2 learnt from 600 000 bases drawn from it has its probabilities (within 0.05; the rarest contexts are seen a few thousand times)", worst, 0.05);
  }

  // ---- a simulated genome ------------------------------------------------------------------------------------
  {
    const Simulated sim = simulate_genome(500, 17);
    GeneReport report;
    const auto t0 = std::chrono::steady_clock::now();
    const auto genes = find_genes(sim.genome, {}, &report);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const Evaluation e = evaluate(sim.genes, genes);
    std::string rounds;
    for (auto n : report.genes_per_round) rounds += (rounds.empty() ? "" : ", ") + std::to_string(n);
    s.require("on a simulated genome of " + std::to_string(sim.genome.size() / 1000) + " kb with 500 planted genes, the self-trained finder finds their stops with sensitivity and specificity at least 0.95",
              e.sensitivity >= 0.95 && e.specificity >= 0.95,
              "Sn " + fixed(e.sensitivity) + ", Sp " + fixed(e.specificity) + ", starts right " + fixed(e.starts) + "; " + std::to_string(report.training_orfs) + " long ORFs to begin, genes by round " + rounds + "; " +
                  fixed(seconds, 2) + " s");
    s.require("its starts are right for most of the genes it finds (at least 0.7)", e.starts >= 0.7, fixed(e.starts));
    bool sorted = true, frames = true;
    for (std::size_t i = 0; i < genes.size(); i++) {
      if (i && genes[i].begin < genes[i - 1].begin) sorted = false;
      if ((genes[i].end - genes[i].begin) % 3 || genes[i].score <= 0) frames = false;
    }
    s.require("the genes are sorted, whole codons long, with positive log-likelihood ratios", sorted && frames);
  }

  // ---- Escherichia coli K-12 ----------------------------------------------------------------------------------
  {
    const char* env = std::getenv("SCIENCE_CACHE");
    const std::filesystem::path cache = env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE);
    const std::filesystem::path file = cache / "ncbi-genomes" / "NC_000913.3.gb";
    if (!std::filesystem::exists(file)) {
      s.skip("Escherichia coli K-12 (data/cache/ncbi-genomes/NC_000913.3.gb) is not in the data cache: the finder is not measured on it");
    } else {
      const auto recs = sequence::read_genbank(read_file(file));
      const Codes genome = sequence::encode_dna(recs.at(0).sequence);
      std::vector<TrueGene> cds;
      std::size_t long_enough = 0;
      std::map<std::string, std::size_t> start_codons;
      for (const auto& f : recs[0].features) {
        if (f.key != "CDS" || f.has("pseudo") || f.location.remote || f.location.parts.empty()) continue;
        TrueGene t;
        t.begin = f.location.lower();
        t.end = f.location.upper();
        t.reverse = f.location.reverse();
        cds.push_back(t);
        long_enough += (t.end - t.begin) / 3 >= 80;
        const Codes x = sequence::extract(genome, f.location);
        if (x.size() >= 3) start_codons[sequence::decode_dna(Codes(x.begin(), x.begin() + 3))]++;
      }
      std::string usage;
      for (const char* c : {"ATG", "GTG", "TTG"}) usage += std::string(usage.empty() ? "" : ", ") + c + " " + fixed(static_cast<double>(start_codons[c]) / static_cast<double>(cds.size()));
      GeneReport report;
      const auto t0 = std::chrono::steady_clock::now();
      const auto genes = find_genes(genome, {}, &report);
      const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      const Evaluation e = evaluate(cds, genes);
      std::string rounds;
      for (auto n : report.genes_per_round) rounds += (rounds.empty() ? "" : ", ") + std::to_string(n);
      s.require("E. coli K-12 MG1655: against " + std::to_string(cds.size()) + " annotated CDSs (" + std::to_string(long_enough) + " of 80 codons or more), stops found with sensitivity at least 0.85 and specificity at least 0.8",
                e.sensitivity >= 0.85 && e.specificity >= 0.8,
                "Sn " + fixed(e.sensitivity) + " (" + fixed(static_cast<double>(e.matched) / static_cast<double>(long_enough)) + " of those long enough), Sp " + fixed(e.specificity) + ", " +
                    std::to_string(genes.size()) + " predicted, starts right " + fixed(e.starts) + "; genes by round " + rounds + "; annotated starts " + usage + "; " + fixed(seconds, 2) + " s");
    }
  }
  return s.done();
}
