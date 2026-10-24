// What the Trees view (trees.gene, trees.dnds) and the toolbox (tool.tree.build, tool.align.dnds)
// compute alike: the tree of an alignment of DNA, by UPGMA or neighbour joining on distances
// corrected by a model or by maximum likelihood under it, with the bootstrap; and dN/dS of two
// coding sequences. They differ only in where their sequences come from (a gene of chosen genomes,
// or the user's sequences), so the cores are here, in a header both units include, and not in
// either ops header.
//
// Every computation is the domains' (evolution::phylogeny, evolution::selection, genomics::align):
// these only choose and report.

#pragma once

#include "genomes.hpp"
#include "shared.hpp"

#include <genomics/align/codon.hpp>

#include <evolution/phylogeny/alignment.hpp>
#include <evolution/phylogeny/bootstrap.hpp>
#include <evolution/phylogeny/distance.hpp>
#include <evolution/phylogeny/distance_trees.hpp>
#include <evolution/phylogeny/likelihood.hpp>
#include <evolution/phylogeny/models.hpp>
#include <evolution/phylogeny/splits.hpp>
#include <evolution/phylogeny/tree.hpp>
#include <evolution/selection/dnds.hpp>

#include <array>
#include <chrono>
#include <cmath>
#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace workbench {

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

// How a tree is to be built, as a request asks.
struct TreeOptions {
  std::string method = "nj", model = "HKY85";
  bool gamma = true;
  std::size_t bootstrap = 100;
  std::uint64_t seed = 1;
};
inline TreeOptions tree_options(const Json& request) {
  TreeOptions o;
  o.method = str(request, "method", "nj");
  o.model = str(request, "model", "HKY85");
  if (o.method != "upgma" && o.method != "nj" && o.method != "ml") throw Fault("a tree is built by upgma, nj or ml");
  model_named(o.model, {0.25, 0.25, 0.25, 0.25});  // a fault for a model it does not know
  o.gamma = flag(request, "gamma", true);
  o.bootstrap = static_cast<std::size_t>(bounded(request, "bootstrap", 100, 0, 1000));
  o.seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
  return o;
}

// What the columns of aligned rows of DNA letters hold: the variable ones (two bases or more) and
// the parsimony-informative ones (two bases or more, each in two rows or more).
inline Json alignment_columns(const std::vector<std::string>& rows) {
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
  return Json::object({{"columns", columns}, {"variable", variable}, {"informative", informative}, {"variable_columns", columns <= 20000 ? Json(variable_columns) : Json()}});
}

// The tree of aligned rows of DNA letters ('-' a gap), its supports from the bootstrap (each
// replicate a task on the pool), midpoint-rooted unless by UPGMA, and ladderized; with the
// distances, the fitted model by likelihood, and the Robinson–Foulds distances to the trees of the
// other methods.
struct TreeResult {
  Json json;        // tree, compare, distances, distance_model, and for ml log_likelihood, parameters, interchanges; bootstrap
  phy::Tree tree;
  double seconds = 0;  // the tree's, before the bootstrap
};
inline TreeResult tree_of_alignment(const std::vector<std::string>& names, const std::vector<std::string>& rows, const TreeOptions& o, calx::scheduler::ThreadPool& pool) {
  using clock = std::chrono::steady_clock;
  const auto since = [](clock::time_point t) { return std::chrono::duration<double>(clock::now() - t).count(); };
  const phy::Patterns data = phy::compress(phy::alignment_from_letters(names, rows));
  auto t = clock::now();
  phy::DistanceOptions dopt;
  dopt.model = distance_named(o.model);
  const phy::DistanceMatrix D = phy::distances(data, dopt, &pool);
  const phy::Tree nj = phy::neighbour_joining(D), up = phy::upgma(D);
  TreeResult res;
  Json& out = res.json;
  out = Json::object();
  phy::Tree tree;
  std::function<phy::Tree(const phy::Patterns&, std::size_t)> build;
  if (o.method == "upgma") {
    tree = up;
    build = [dopt](const phy::Patterns& rep, std::size_t) { return phy::upgma(phy::distances(rep, dopt)); };
  } else if (o.method == "nj") {
    tree = nj;
    build = [dopt](const phy::Patterns& rep, std::size_t) { return phy::neighbour_joining(phy::distances(rep, dopt)); };
  } else {
    phy::RateSpec rates;
    rates.categories = o.gamma ? 4 : 1;
    rates.alpha = 0.5;
    const auto fit = phy::maximum_likelihood(data, model_named(o.model, phy::base_frequencies(data)), rates, nj);
    tree = fit.tree;
    out.set("log_likelihood", fit.log_likelihood);
    Json par = Json::object();
    const auto& m = fit.model;
    if (fit.model.kind == phy::ModelKind::k80 || fit.model.kind == phy::ModelKind::hky85) par.set("κ", m.exchange[phy::AG]);
    if (fit.model.kind == phy::ModelKind::tn93) { par.set("κ purines", m.exchange[phy::AG]); par.set("κ pyrimidines", m.exchange[phy::CT]); }
    if (fit.model.kind == phy::ModelKind::gtr)
      for (const auto& [k, name] : std::vector<std::pair<int, const char*>>{{phy::AC, "AC"}, {phy::AG, "AG"}, {phy::AT, "AT"}, {phy::CG, "CG"}, {phy::CT, "CT"}}) par.set(name, m.exchange[static_cast<std::size_t>(k)] / m.exchange[phy::GT]);
    if (o.gamma) par.set("α", fit.rates.alpha);
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
  res.seconds = since(t);
  if (o.bootstrap > 0) {
    t = clock::now();
    const auto bs = phy::bootstrap(data, tree, o.bootstrap, o.seed, pool, build);
    tree = bs.tree;
    out.set("bootstrap", Json::object({{"replicates", o.bootstrap}, {"seconds", since(t)}, {"threads", pool.size()}, {"trace", trace_json(bs.trace)}}));
  }
  if (o.method != "upgma") tree = phy::midpoint_root(tree);
  tree = phy::ladderize(tree);
  out.set("tree", tree_json(tree));
  Json compare = Json::object();
  if (o.method != "nj") compare.set("neighbour joining", phy::robinson_foulds(tree, nj));
  if (o.method != "upgma") compare.set("UPGMA", phy::robinson_foulds(tree, up));
  out.set("compare", compare);
  Json drows = Json::array();
  for (std::size_t i = 0; i < D.n; i++) {
    Json row = Json::array();
    for (std::size_t j = 0; j < D.n; j++) row.push(finite_or_null(D(i, j)));
    drows.push(row);
  }
  out.set("distances", drows);
  out.set("distance_model", o.model == "JC69" || o.model == "K80" || o.model == "F81" ? o.model : std::string("TN93"));
  res.tree = std::move(tree);
  return res;
}

// dN/dS of two coding sequences under a genetic code: aligned through their proteins (whole
// codons, a final stop and a trailing partial codon left out), then synonymous and non-synonymous
// sites and differences by Nei and Gojobori (1986), Jukes–Cantor corrected. With `proteins`, the
// two proteins as their codons were aligned.
inline Json dnds_pair(const seq::Codes& a, const seq::Codes& b, const seq::GeneticCode& code, const evolution::selection::CodonTable& codons, bool proteins = false) {
  const auto al = genomics::align::codon_alignment(a, b, code);
  const auto r = evolution::selection::nei_gojobori(al.x, al.y, codons);
  Json out = Json::object({{"S", finite_or_null(r.S, 1e3)}, {"N", finite_or_null(r.N, 1e3)}, {"Sd", finite_or_null(r.Sd, 1e3)}, {"Nd", finite_or_null(r.Nd, 1e3)},
                           {"pS", finite_or_null(r.pS)}, {"pN", finite_or_null(r.pN)}, {"dS", finite_or_null(r.dS)}, {"dN", finite_or_null(r.dN)},
                           {"omega", finite_or_null(r.omega)}, {"saturated_S", r.saturated_S}, {"saturated_N", r.saturated_N},
                           {"codons", r.codons}, {"skipped", r.skipped}, {"gap_codons", al.x.size() / 3 - al.pairs}});
  if (proteins) {
    out.set("protein_a", al.protein_x);
    out.set("protein_b", al.protein_y);
    out.set("partial", Json::array({al.partial_x, al.partial_y}));
    out.set("stops", Json::array({al.stop_x, al.stop_y}));
  }
  return out;
}

}  // namespace detail

}  // namespace workbench
