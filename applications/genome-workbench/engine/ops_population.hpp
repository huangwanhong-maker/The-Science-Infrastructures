// The operations of the Populations view, on the evolution domain's population module: a
// coalescent sample with its genealogy (the marginal tree at a point of the chromosome, its
// mutations), its haplotypes, its site frequency spectrum and the estimators of θ; Wright–Fisher
// replicates with a few trajectories to draw, against Kimura's fixation probability; and a
// structured, admixed sample with its principal components and the ancestry the EM recovers,
// beside the truth.

#pragma once

#include "protocol.hpp"
#include "genomes.hpp"
#include "shared.hpp"

#include <evolution/population/coalescent.hpp>
#include <evolution/population/diversity.hpp>
#include <evolution/population/genotypes.hpp>
#include <evolution/population/structure.hpp>
#include <evolution/population/wright_fisher.hpp>

#include <string>
#include <vector>

namespace workbench {

namespace pop = evolution::population;

namespace detail {

inline pop::Demography demography_of(const Json& d) {
  const std::string kind = str(d, "kind", "constant");
  if (kind == "growth") return pop::Demography::exponential(num(d, "alpha", 5));
  if (kind == "bottleneck") return pop::Demography::bottleneck(num(d, "start", 0.05), num(d, "duration", 0.05), num(d, "size", 0.05));
  return pop::Demography::constant();
}

inline Json tree_json(const pop::Tree& t) {
  Json nodes = Json::array();
  for (std::size_t v = 0; v < t.time.size(); v++) nodes.push(Json::array({std::round(t.time[v] * 1e6) / 1e6, t.parent[v]}));
  Json muts = Json::array();
  for (const auto& m : t.mutations) muts.push(Json::array({m.node, std::round(m.time * 1e6) / 1e6, std::round(m.position * 1e6) / 1e6}));
  return Json::object({{"left", t.left}, {"right", t.right}, {"root", t.root}, {"tmrca", t.tmrca}, {"length", t.length}, {"nodes", nodes}, {"mutations", muts}});
}

}  // namespace detail

BENCH_OP("pop.coalescent") {
  pop::CoalescentSpec spec;
  spec.samples = static_cast<std::size_t>(bounded(request, "samples", 20, 2, 500));
  spec.theta = num(request, "theta", 10);
  spec.rho = num(request, "rho", 0);
  if (request["demography"].is_object()) spec.demography = detail::demography_of(request["demography"]);
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
  const auto g = pop::coalescent(spec, seed);
  // The tree at the point asked (or the middle), and the run of marginal trees with their T_MRCA.
  const double x = std::clamp(num(request, "at", 0.5), 0.0, 1.0 - 1e-12);
  Json out = Json::object({{"samples", spec.samples}, {"theta", spec.theta}, {"rho", spec.rho}, {"recombinations", g.recombinations}, {"trees", pop::tree_count(g)}});
  out.set("tree", detail::tree_json(pop::tree_at(g, x)));
  Json run = Json::array();
  for (const auto& t : pop::marginal_trees(g, 600)) run.push(Json::array({t.left, t.right, std::round(t.tmrca * 1e5) / 1e5}));
  out.set("tmrca_along", run);
  // Haplotypes, the spectrum and the estimators.
  const auto h = pop::haplotypes(g);
  const auto counts = pop::allele_counts(h);
  const auto xi = pop::site_frequency_spectrum(counts, spec.samples);
  const auto s = pop::summarise(xi);
  const double S = s.S, a = pop::watterson_a(spec.samples);
  out.set("spectrum", numbers(xi, 3));
  out.set("segregating", S);
  out.set("theta_w", a > 0 ? S / a : 0);
  out.set("pi", s.pi);
  out.set("theta_h", s.theta_h);
  out.set("tajima_d", S > 0 ? Json(pop::tajima_d(spec.samples, S, s.pi)) : Json());
  out.set("fay_wu_h", s.pi - s.theta_h);
  out.set("expected_tmrca", pop::expected_tmrca(spec.samples));
  if (h.sites * h.sequences <= 400'000) {
    Json rows = Json::array();
    for (std::size_t q = 0; q < h.sequences; q++) {
      std::string row(h.sites, '0');
      for (std::size_t j = 0; j < h.sites; j++) row[j] = h(q, j) == 1 ? '1' : '0';
      rows.push(row);
    }
    out.set("haplotypes", rows);
    out.set("positions", numbers(h.positions, 5));
  }
  return out;
}

BENCH_OP("pop.wright_fisher") {
  pop::WrightFisherSpec w;
  w.individuals = static_cast<std::size_t>(bounded(request, "individuals", 100, 2, 1'000'000));
  w.demes = static_cast<std::size_t>(bounded(request, "demes", 1, 1, 50));
  w.s = num(request, "s", 0);
  w.h = num(request, "h", 0.5);
  w.mu = num(request, "mu", 0);
  w.nu = num(request, "nu", 0);
  w.m = num(request, "m", 0);
  w.p0 = num(request, "p0", 0.1);
  w.generations = static_cast<std::size_t>(bounded(request, "generations", 400, 1, 200000));
  const auto replicates = static_cast<std::size_t>(bounded(request, "replicates", 400, 1, 200000));
  const auto draw = static_cast<std::size_t>(bounded(request, "draw", 40, 0, 400));
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
  const std::size_t every = std::max<std::size_t>(1, w.generations / 400);
  const auto sum = pop::wright_fisher_replicates(w, seed, replicates, worker.pool, every);
  Json paths = Json::array();
  for (std::size_t r = 0; r < std::min(draw, replicates); r++) {
    const auto t = pop::wright_fisher_trajectory(w, seed, static_cast<std::uint32_t>(r));
    std::vector<double> mean;
    for (std::size_t gen = 0; gen <= t.generations(); gen += every) {
      double p = 0;
      for (std::size_t d = 0; d < t.demes; d++) p += t.at(gen, d);
      mean.push_back(p / static_cast<double>(t.demes));
    }
    paths.push(numbers(mean, 4));
  }
  return Json::object({{"individuals", w.individuals}, {"demes", w.demes}, {"s", w.s}, {"h", w.h}, {"p0", w.p0}, {"generations", w.generations}, {"every", every},
                       {"replicates", sum.replicates}, {"fixed", sum.fixed}, {"lost", sum.lost}, {"fixation", sum.fixation()},
                       {"kimura", pop::fixation_probability(static_cast<double>(w.individuals * w.demes), w.s, w.h, w.p0)},
                       {"mean_settling", sum.mean_settling}, {"mean_p", numbers(sum.mean_p, 5)}, {"h_s", numbers(sum.h_s, 5)}, {"h_t", numbers(sum.h_t, 5)},
                       {"fst", w.demes > 1 ? Json(sum.fst(sum.h_t.size() / 2, sum.h_t.size())) : Json()},
                       {"island_fst", w.demes > 1 && w.m > 0 ? Json(pop::island_fst(static_cast<double>(w.individuals), w.m)) : Json()},
                       {"paths", paths}, {"trace", detail::trace_json(sum.trace)}});
}

BENCH_OP("pop.structure") {
  pop::StructureSpec spec;
  spec.populations = static_cast<std::size_t>(bounded(request, "populations", 3, 2, 8));
  spec.sites = static_cast<std::size_t>(bounded(request, "sites", 4000, 100, 100000));
  spec.fst = {num(request, "fst", 0.1)};
  spec.unadmixed.assign(spec.populations, static_cast<std::size_t>(bounded(request, "each", 25, 2, 500)));
  spec.admixed = static_cast<std::size_t>(bounded(request, "admixed", 30, 0, 2000));
  spec.alpha = num(request, "alpha", 1);
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
  const auto sample = pop::simulate_structure(spec, seed);
  const auto pc = pop::pca(sample.genotypes, 4, &worker.pool);
  pop::AdmixtureOptions o;
  o.populations = static_cast<std::size_t>(bounded(request, "k", static_cast<long long>(spec.populations), 1, 10));
  o.seed = seed;
  o.restarts = 1;
  const auto mix = pop::admixture(sample.genotypes, o, &worker.pool);
  const std::size_t n = sample.genotypes.samples;
  Json coords = Json::array();
  for (std::size_t i = 0; i < n; i++) {
    std::vector<double> c(pc.vectors.begin() + static_cast<std::ptrdiff_t>(i * pc.components), pc.vectors.begin() + static_cast<std::ptrdiff_t>((i + 1) * pc.components));
    coords.push(numbers(c, 5));
  }
  Json q = Json::array(), truth = Json::array();
  for (std::size_t i = 0; i < n; i++) {
    q.push(numbers(std::vector<double>(mix.q.begin() + static_cast<std::ptrdiff_t>(i * mix.populations), mix.q.begin() + static_cast<std::ptrdiff_t>((i + 1) * mix.populations)), 4));
    truth.push(numbers(std::vector<double>(sample.ancestry.begin() + static_cast<std::ptrdiff_t>(i * spec.populations), sample.ancestry.begin() + static_cast<std::ptrdiff_t>((i + 1) * spec.populations)), 4));
  }
  return Json::object({{"samples", n}, {"sites", spec.sites}, {"populations", spec.populations}, {"k", mix.populations},
                       {"pca", Json::object({{"coordinates", coords}, {"values", numbers(pc.values, 4)}, {"tracy_widom", numbers(pc.tracy_widom, 3)}})},
                       {"admixture", Json::object({{"q", q}, {"log_likelihood", mix.log_likelihood}, {"iterations", mix.iterations}, {"converged", mix.converged}})},
                       {"truth", truth}, {"population", Json(std::vector<double>(sample.population.begin(), sample.population.end()))},
                       {"admixed", Json(std::vector<double>(sample.admixed.begin(), sample.admixed.end()))}});
}

}  // namespace workbench
