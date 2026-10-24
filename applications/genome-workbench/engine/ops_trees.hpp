// The operations of the Trees view: whole genomes compared by their sketches (Mash distances) and
// placed on a tree by neighbour joining; a gene of each chosen genome aligned and its tree built
// by distances or by maximum likelihood, with the bootstrap; and the pressure of selection on the
// gene between each pair, dN/dS by Nei and Gojobori on its codon alignment.
//
// Every computation is the domains' (genomics::sketch, genomics::msa, genomics::align,
// evolution::phylogeny, evolution::selection): the engine only gathers the inputs and reports.

#pragma once

#include "genomes.hpp"
#include "shared.hpp"
#include "shared_trees.hpp"

#include <genomics/sketch/distance.hpp>
#include <genomics/sketch/minhash.hpp>
#include <genomics/msa/progressive.hpp>

#include <chrono>
#include <cmath>
#include <string>
#include <vector>

namespace workbench {

namespace detail {

// The genomes a request lists, each loaded (or taken from the worker's cache), with a short name.
struct Listed {
  std::vector<std::shared_ptr<Genome>> genomes;
  std::vector<std::string> names;
};
inline Listed listed_genomes(Worker& worker, const Json& request, std::size_t at_least, std::size_t at_most) {
  const Json& list = request["genomes"];
  if (!list.is_array() || list.size() < at_least) throw Fault("choose at least " + std::to_string(at_least) + " genomes");
  if (list.size() > at_most) throw Fault("at most " + std::to_string(at_most) + " genomes at once");
  Listed out;
  for (const auto& spec : list.as_array()) {
    auto g = worker.genome(Json::object({{"genome", spec}}));
    std::string name = g->organism.empty() ? g->title : g->organism;
    // The same organism twice (two records of one species) keeps both, told apart.
    for (const auto& seen : out.names)
      if (seen == name) name += " (" + g->key.substr(g->key.find_last_of("/\\") + 1) + ")";
    out.names.push_back(name);
    out.genomes.push_back(std::move(g));
  }
  return out;
}

}  // namespace detail

// Every listed genome sketched (the smallest hashes of its canonical k-mers), the Mash distances
// of every pair as a task graph of blocks, and their tree by neighbour joining, midpoint-rooted.
BENCH_OP("trees.sketches") {
  const auto L = detail::listed_genomes(worker, request, 3, 200);
  const int k = static_cast<int>(bounded(request, "k", 21, 9, 31));
  const auto size = static_cast<std::size_t>(bounded(request, "size", 1000, 50, 100000));
  std::vector<std::vector<seq::Codes>> genomes;
  for (const auto& g : L.genomes) {
    std::vector<seq::Codes> contigs;
    for (const auto& c : g->contigs) contigs.push_back(c.codes);
    genomes.push_back(std::move(contigs));
  }
  const auto t0 = std::chrono::steady_clock::now();
  const auto sketches = genomics::sketch::minhash_all(worker.pool, genomes, k, size);
  const auto t1 = std::chrono::steady_clock::now();
  const auto D = genomics::sketch::mash_distances(worker.pool, sketches);
  const auto t2 = std::chrono::steady_clock::now();
  const phy::DistanceMatrix m(L.names, D.values);
  const phy::Tree tree = phy::ladderize(phy::midpoint_root(phy::neighbour_joining(m)));
  Json rows = Json::array();
  for (std::size_t i = 0; i < m.n; i++) {
    Json row = Json::array();
    for (std::size_t j = 0; j < m.n; j++) row.push(detail::finite_or_null(m(i, j)));
    rows.push(row);
  }
  const std::size_t n = L.names.size();
  return Json::object({{"tree", detail::tree_json(tree)}, {"names", L.names}, {"distances", rows}, {"genomes", n}, {"pairs", n * (n - 1) / 2},
                       {"k", k}, {"size", size}, {"tasks", D.tasks},
                       {"sketch_seconds", std::chrono::duration<double>(t1 - t0).count()}, {"distance_seconds", std::chrono::duration<double>(t2 - t1).count()}});
}

namespace detail {

// The coding sequence of a gene in each listed genome, by its name.
inline std::vector<seq::Codes> gene_of_each(const Listed& L, const std::string& gene) {
  std::vector<seq::Codes> out;
  for (std::size_t i = 0; i < L.genomes.size(); i++) {
    const auto& g = *L.genomes[i];
    const auto at = find_gene(g, gene);
    if (!at) throw Fault("no gene " + gene + " in " + L.names[i]);
    const auto& c = g.contigs[at->contig];
    out.push_back(seq::extract(c.codes, c.features[at->feature].location));
  }
  return out;
}

}  // namespace detail

// A gene of each listed genome aligned (progressive alignment, MUSCLE's three stages), and its
// tree: by UPGMA or neighbour joining on distances corrected by a model, or by maximum
// likelihood under the model (with four Γ rates when asked), its branches supported by the
// bootstrap, each replicate a task on the worker's threads.
BENCH_OP("trees.gene") {
  const auto L = detail::listed_genomes(worker, request, 3, 40);
  const std::string gene = str(request, "gene");
  if (flag(request, "protein", false)) throw Fault("the models here are of nucleotides; trees of proteins need models of amino acids, not yet written");
  const auto o = detail::tree_options(request);
  using clock = std::chrono::steady_clock;

  // The alignment.
  const auto t = clock::now();
  const auto seqs = detail::gene_of_each(L, gene);
  const auto aligned = genomics::msa::progressive(seqs, L.names, genomics::msa::Alphabet::dna, genomics::msa::default_options(genomics::msa::Alphabet::dna), worker.pool);
  const double align_seconds = std::chrono::duration<double>(clock::now() - t).count();
  std::vector<std::string> rows;
  for (const auto& r : aligned.msa.rows) rows.push_back(seq::decode_dna(r));

  // The tree, its supports and what it is compared with.
  auto built = detail::tree_of_alignment(L.names, rows, o, worker.pool);
  Json out = std::move(built.json);
  Json alignment = detail::alignment_columns(rows);
  alignment.set("names", L.names);
  alignment.set("rows", rows);
  alignment.set("method", "progressive (MUSCLE's k-mer guide tree, re-estimated tree and refinement); sum of pairs " + std::to_string(aligned.sp));
  alignment.set("seconds", align_seconds);
  out.set("alignment", alignment);
  out.set("protein", false);
  out.set("seconds", Json::object({{"align", align_seconds}, {"tree", built.seconds}}));
  return out;
}

// dN/dS of a gene between every pair of the listed genomes: the coding sequences aligned through
// their proteins (whole codons), then synonymous and non-synonymous sites and differences by Nei
// and Gojobori (1986), Jukes–Cantor corrected, under the genomes' genetic code.
BENCH_OP("trees.dnds") {
  const auto L = detail::listed_genomes(worker, request, 2, 30);
  const std::string gene = str(request, "gene");
  const auto cds = detail::gene_of_each(L, gene);
  int table = 0;
  for (std::size_t i = 0; i < L.genomes.size(); i++) {
    const auto& g = *L.genomes[i];
    const auto at = detail::find_gene(g, gene);
    if (g.contigs[at->contig].features[at->feature].key != "CDS") throw Fault(gene + " of " + L.names[i] + " is not a coding sequence");
    if (table && table != g.table) throw Fault("these genomes use different genetic codes (" + std::to_string(table) + " and " + std::to_string(g.table) + ")");
    table = g.table;
  }
  const auto& code = seq::genetic_code(table);
  const evolution::selection::CodonTable codons(code);
  Json pairs = Json::array();
  for (std::size_t i = 0; i < cds.size(); i++)
    for (std::size_t j = i + 1; j < cds.size(); j++) {
      Json pair = Json::object({{"a", i}, {"b", j}});
      for (auto& [k, v] : detail::dnds_pair(cds[i], cds[j], code, codons).as_object()) pair.set(k, std::move(v));
      pairs.push(std::move(pair));
    }
  Json lengths = Json::array();
  for (const auto& x : cds) lengths.push(x.size());
  return Json::object({{"names", L.names}, {"gene", gene}, {"table", table}, {"code", code.name}, {"lengths", lengths}, {"pairs", pairs}});
}

}  // namespace workbench
