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

#include <genomics/align/codon.hpp>
#include <genomics/sketch/distance.hpp>
#include <genomics/sketch/minhash.hpp>
#include <genomics/msa/progressive.hpp>

#include <evolution/phylogeny/distance.hpp>
#include <evolution/phylogeny/alignment.hpp>
#include <evolution/phylogeny/bootstrap.hpp>
#include <evolution/phylogeny/distance_trees.hpp>
#include <evolution/phylogeny/likelihood.hpp>
#include <evolution/phylogeny/models.hpp>
#include <evolution/phylogeny/splits.hpp>
#include <evolution/phylogeny/tree.hpp>
#include <evolution/selection/dnds.hpp>

#include <chrono>
#include <functional>
#include <cmath>
#include <string>
#include <vector>

namespace atlas {

namespace phy = evolution::phylogeny;

namespace detail {

inline Json finite_or_null(double x, double scale = 1e6) { return std::isfinite(x) ? Json(std::round(x * scale) / scale) : Json(); }

// A tree as the page draws it: each node's parent (−1 at the root), branch length, label and
// support (null where absent).
inline Json tree_json(const phy::Tree& t) {
  Json nodes = Json::array();
  double total = 0;
  for (std::size_t v = 0; v < t.nodes.size(); v++) {
    const auto& n = t.nodes[v];
    nodes.push(Json::object({{"parent", n.parent}, {"length", finite_or_null(n.length)}, {"label", n.label}, {"support", finite_or_null(n.support, 1e4)}}));
    if (n.parent >= 0 && std::isfinite(n.length)) total += n.length;
  }
  return Json::object({{"nodes", nodes}, {"length", total}, {"rooted", t.rooted}});
}

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
ATLAS_OP("trees.sketches") {
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

// A substitution model by name, its frequencies the alignment's and its ratios to be fitted.
inline phy::Model model_named(const std::string& name, const std::array<double, 4>& pi) {
  if (name == "JC69") return phy::jc69();
  if (name == "K80") return phy::k80(2);
  if (name == "F81") return phy::f81(pi);
  if (name == "HKY85") return phy::hky85(2, pi);
  if (name == "TN93") return phy::tn93(2, 2, pi);
  if (name == "GTR") return phy::gtr({1, 1, 1, 1, 1, 1}, pi);
  throw Fault("a model is JC69, K80, F81, HKY85, TN93 or GTR");
}

// The distance that corrects for the model's kind of substitution (no distance of its own for
// HKY85 and GTR: TN93's, which holds HKY85 and approximates GTR).
inline phy::DistanceModel distance_named(const std::string& name) {
  if (name == "JC69") return phy::DistanceModel::jc69;
  if (name == "K80") return phy::DistanceModel::k80;
  if (name == "F81") return phy::DistanceModel::f81;
  return phy::DistanceModel::tn93;
}

}  // namespace detail

// A gene of each listed genome aligned (progressive alignment, MUSCLE's three stages), and its
// tree: by UPGMA or neighbour joining on distances corrected by a model, or by maximum
// likelihood under the model (with four Γ rates when asked), its branches supported by the
// bootstrap, each replicate a task on the worker's threads.
ATLAS_OP("trees.gene") {
  const auto L = detail::listed_genomes(worker, request, 3, 40);
  const std::string gene = str(request, "gene");
  if (flag(request, "protein", false)) throw Fault("the models here are of nucleotides; trees of proteins need models of amino acids, not yet written");
  const std::string method = str(request, "method", "nj"), model_name = str(request, "model", "HKY85");
  if (method != "upgma" && method != "nj" && method != "ml") throw Fault("a tree is built by upgma, nj or ml");
  const bool gamma = flag(request, "gamma", true);
  const auto B = static_cast<std::size_t>(bounded(request, "bootstrap", 100, 0, 1000));
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
  using clock = std::chrono::steady_clock;
  const auto since = [](clock::time_point t) { return std::chrono::duration<double>(clock::now() - t).count(); };

  // The alignment.
  auto t = clock::now();
  const auto seqs = detail::gene_of_each(L, gene);
  const auto aligned = genomics::msa::progressive(seqs, L.names, genomics::msa::Alphabet::dna, genomics::msa::default_options(genomics::msa::Alphabet::dna), worker.pool);
  const double align_seconds = since(t);
  std::vector<std::string> rows;
  for (const auto& r : aligned.msa.rows) rows.push_back(seq::decode_dna(r));
  const std::size_t columns = rows.empty() ? 0 : rows[0].size();
  std::size_t variable = 0, informative = 0;
  std::vector<double> variable_columns(columns, 0);
  for (std::size_t j = 0; j < columns; j++) {
    std::array<int, 4> n{};
    for (const auto& r : rows) {
      const auto k = std::string_view("ACGT").find(r[j]);
      if (k != std::string_view::npos) n[k]++;
    }
    int kinds = 0, twice = 0;
    for (int c : n) { kinds += c > 0; twice += c >= 2; }
    if (kinds > 1) { variable++; variable_columns[j] = 1; }
    if (twice >= 2) informative++;
  }
  const phy::Patterns data = phy::compress(phy::alignment_from_letters(L.names, rows));

  // The distances and the tree.
  t = clock::now();
  phy::DistanceOptions dopt;
  dopt.model = detail::distance_named(model_name);
  const phy::DistanceMatrix D = phy::distances(data, dopt, &worker.pool);
  const phy::Tree nj = phy::neighbour_joining(D), up = phy::upgma(D);
  Json out = Json::object();
  phy::Tree tree;
  std::function<phy::Tree(const phy::Patterns&, std::size_t)> build;
  if (method == "upgma") {
    tree = up;
    build = [dopt](const phy::Patterns& rep, std::size_t) { return phy::upgma(phy::distances(rep, dopt)); };
  } else if (method == "nj") {
    tree = nj;
    build = [dopt](const phy::Patterns& rep, std::size_t) { return phy::neighbour_joining(phy::distances(rep, dopt)); };
  } else {
    phy::RateSpec rates;
    rates.categories = gamma ? 4 : 1;
    rates.alpha = 0.5;
    const auto fit = phy::maximum_likelihood(data, detail::model_named(model_name, phy::base_frequencies(data)), rates, nj);
    tree = fit.tree;
    out.set("log_likelihood", fit.log_likelihood);
    Json par = Json::object();
    const auto& m = fit.model;
    if (fit.model.kind == phy::ModelKind::k80 || fit.model.kind == phy::ModelKind::hky85) par.set("κ", m.exchange[phy::AG]);
    if (fit.model.kind == phy::ModelKind::tn93) { par.set("κ purines", m.exchange[phy::AG]); par.set("κ pyrimidines", m.exchange[phy::CT]); }
    if (fit.model.kind == phy::ModelKind::gtr)
      for (const auto& [k, name] : std::vector<std::pair<int, const char*>>{{phy::AC, "AC"}, {phy::AG, "AG"}, {phy::AT, "AT"}, {phy::CG, "CG"}, {phy::CT, "CT"}}) par.set(name, m.exchange[static_cast<std::size_t>(k)] / m.exchange[phy::GT]);
    if (gamma) par.set("α", fit.rates.alpha);
    out.set("parameters", par);
    out.set("interchanges", fit.interchanges);
    // Each replicate's tree by likelihood under the model fitted here, its topology searched
    // from its own neighbour-joining start.
    phy::SearchOptions so;
    so.model_parameters = false;
    so.cycles = 2;
    const phy::Model fitted = fit.model;
    const phy::RateSpec fitted_rates = fit.rates;
    build = [fitted, fitted_rates, so, dopt](const phy::Patterns& rep, std::size_t) {
      return phy::maximum_likelihood(rep, fitted, fitted_rates, phy::neighbour_joining(phy::distances(rep, dopt)), so).tree;
    };
  }
  const double tree_seconds = since(t);
  if (B > 0) {
    t = clock::now();
    const auto bs = phy::bootstrap(data, tree, B, seed, worker.pool, build);
    tree = bs.tree;
    out.set("bootstrap", Json::object({{"replicates", B}, {"seconds", since(t)}, {"threads", worker.pool.size()}, {"trace", detail::trace_json(bs.trace)}}));
  }
  if (method != "upgma") tree = phy::midpoint_root(tree);
  tree = phy::ladderize(tree);
  out.set("tree", detail::tree_json(tree));
  Json compare = Json::object();
  if (method != "nj") compare.set("neighbour joining", phy::robinson_foulds(tree, nj));
  if (method != "upgma") compare.set("UPGMA", phy::robinson_foulds(tree, up));
  out.set("compare", compare);
  Json drows = Json::array();
  for (std::size_t i = 0; i < D.n; i++) {
    Json row = Json::array();
    for (std::size_t j = 0; j < D.n; j++) row.push(detail::finite_or_null(D(i, j)));
    drows.push(row);
  }
  out.set("distances", drows);
  out.set("distance_model", model_name == "JC69" || model_name == "K80" || model_name == "F81" ? model_name : std::string("TN93"));
  out.set("alignment", Json::object({{"names", L.names}, {"rows", rows}, {"columns", columns}, {"variable", variable}, {"informative", informative},
                                     {"variable_columns", columns <= 20000 ? Json(variable_columns) : Json()},
                                     {"method", "progressive (MUSCLE's k-mer guide tree, re-estimated tree and refinement); sum of pairs " + std::to_string(aligned.sp)},
                                     {"seconds", align_seconds}}));
  out.set("protein", false);
  out.set("seconds", Json::object({{"align", align_seconds}, {"tree", tree_seconds}}));
  return out;
}

// dN/dS of a gene between every pair of the listed genomes: the coding sequences aligned through
// their proteins (whole codons), then synonymous and non-synonymous sites and differences by Nei
// and Gojobori (1986), Jukes–Cantor corrected, under the genomes' genetic code.
ATLAS_OP("trees.dnds") {
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
      const auto a = genomics::align::codon_alignment(cds[i], cds[j], code);
      const auto r = evolution::selection::nei_gojobori(a.x, a.y, codons);
      pairs.push(Json::object({{"a", i}, {"b", j}, {"S", detail::finite_or_null(r.S, 1e3)}, {"N", detail::finite_or_null(r.N, 1e3)}, {"Sd", detail::finite_or_null(r.Sd, 1e3)},
                               {"Nd", detail::finite_or_null(r.Nd, 1e3)}, {"dS", detail::finite_or_null(r.dS)}, {"dN", detail::finite_or_null(r.dN)},
                               {"omega", detail::finite_or_null(r.omega)}, {"saturated_S", r.saturated_S}, {"saturated_N", r.saturated_N},
                               {"codons", r.codons}, {"skipped", r.skipped}, {"gap_codons", a.x.size() / 3 - a.pairs}}));
    }
  Json lengths = Json::array();
  for (const auto& x : cds) lengths.push(x.size());
  return Json::object({{"names", L.names}, {"gene", gene}, {"table", table}, {"code", code.name}, {"lengths", lengths}, {"pairs", pairs}});
}

}  // namespace atlas
