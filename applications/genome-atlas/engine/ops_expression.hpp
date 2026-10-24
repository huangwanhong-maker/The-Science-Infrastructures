// The operations of the Expression view: an RNA-seq sample of a simulated transcriptome
// pseudoaligned and quantified, against the truth; and a simulated study of two conditions tested
// gene by gene.
//
// Every computation is the domain's (genomics::expression): the engine draws the inputs a study
// would have (base means, dispersions, which genes change), times the steps and reports.

#pragma once

#include "genomes.hpp"

#include <genomics/expression/differential.hpp>
#include <genomics/expression/pseudoalign.hpp>
#include <genomics/expression/quantify.hpp>
#include <genomics/expression/transcriptome.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/parallel_for.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <vector>

namespace atlas {

namespace gex = genomics::expression;

namespace detail {

inline double rounded(double x, double scale = 1e4) { return std::isfinite(x) ? std::round(x * scale) / scale : 0; }

// Pearson's correlation of log(1 + x) and log(1 + y).
inline double log_correlation(const std::vector<double>& x, const std::vector<double>& y) {
  const std::size_t n = x.size();
  double mx = 0, my = 0;
  for (std::size_t i = 0; i < n; i++) { mx += std::log1p(x[i]); my += std::log1p(y[i]); }
  mx /= static_cast<double>(n);
  my /= static_cast<double>(n);
  double sxy = 0, sxx = 0, syy = 0;
  for (std::size_t i = 0; i < n; i++) {
    const double a = std::log1p(x[i]) - mx, b = std::log1p(y[i]) - my;
    sxy += a * b; sxx += a * a; syy += b * b;
  }
  return sxy / std::sqrt(sxx * syy);
}

// The median of |estimate − truth| / truth over transcripts with a true count of at least 10.
inline double median_relative_error(const std::vector<double>& est, const std::vector<double>& truth) {
  std::vector<double> e;
  for (std::size_t i = 0; i < est.size(); i++)
    if (truth[i] >= 10) e.push_back(std::fabs(est[i] - truth[i]) / truth[i]);
  if (e.empty()) return 0;
  std::nth_element(e.begin(), e.begin() + static_cast<std::ptrdiff_t>(e.size() / 2), e.end());
  return e[e.size() / 2];
}

}  // namespace detail

// A transcriptome (genes of exons, isoforms sharing them, paralogs), its abundances and an RNA-seq
// sample of it; the sample pseudoaligned, quantified by EM and by an equal split of each class,
// and bootstrapped; every estimate against the truth.
ATLAS_OP("expression.quantify") {
  using clock = std::chrono::steady_clock;
  const auto since = [](clock::time_point t) { return std::chrono::duration<double>(clock::now() - t).count(); };
  gex::TranscriptomeSpec spec;
  spec.genes = static_cast<std::size_t>(bounded(request, "genes", 400, 10, 5000));
  spec.paralogs = num(request, "paralogs", 0.15);
  if (!(spec.paralogs >= 0 && spec.paralogs <= 0.5)) throw Fault("a fraction of paralogs between 0 and 0.5");
  const auto N = static_cast<std::size_t>(bounded(request, "fragments", 400000, 1000, 5'000'000));
  const auto B = static_cast<std::size_t>(bounded(request, "bootstrap", 30, 0, 200));
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
  gex::ReadModel model;
  model.length = static_cast<std::size_t>(bounded(request, "length", 75, 31, 300));
  model.paired = flag(request, "paired", false);
  model.error = num(request, "error", 0.003);
  if (!(model.error >= 0 && model.error <= 0.05)) throw Fault("an error rate between 0 and 0.05");

  auto t = clock::now();
  const auto sim = gex::simulate_transcriptome(spec, seed);
  const auto& tx = sim.transcriptome;
  const std::size_t T = tx.size();
  const auto rho = gex::simulate_abundances(T, num(request, "sigma", 1.5), num(request, "zero", 0.1), seed + 1);
  const auto d = gex::FragmentLengths::normal(200, 25);
  // The fragments in chunks on the threads: fragment i draws from its own stream.
  constexpr std::size_t CHUNK = 32768;
  std::vector<std::vector<gex::Fragment>> parts((N + CHUNK - 1) / CHUNK);
  calx::scheduler::parallel_for(worker.pool, parts.size(), 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t c = b; c < e; c++) parts[c] = gex::simulate_fragments(tx, rho, d, model, std::min(CHUNK, N - c * CHUNK), seed + 2, c * CHUNK);
  });
  std::vector<gex::Fragment> fragments;
  fragments.reserve(N);
  for (auto& p : parts) std::move(p.begin(), p.end(), std::back_inserter(fragments));
  const double simulate_seconds = since(t);

  t = clock::now();
  const gex::Index index(tx, 31);
  const double index_seconds = since(t);
  t = clock::now();
  const auto pa = gex::pseudoalign(index, fragments, worker.pool);
  const double pseudo_seconds = since(t);
  const auto eff = gex::effective_lengths(tx.lengths(), d);
  t = clock::now();
  const auto ab = gex::quantify(pa, eff, gex::EmOptions{}, worker.pool);
  const double em_seconds = since(t);
  const auto eq = gex::equal_split(pa.classes, pa.counts, T);
  t = clock::now();
  const auto bs = B ? gex::bootstrap(pa.classes, pa.counts, eff, B, seed + 3, gex::EmOptions{}, worker.pool) : gex::Bootstrap{};
  const double bootstrap_seconds = since(t);

  // The truth: each transcript's fragments that pseudoaligned, and the share of them in classes
  // of their own transcript alone (what makes a transcript identifiable from classes).
  std::vector<double> truth(T, 0), own(T, 0);
  for (std::size_t i = 0; i < fragments.size(); i++) {
    const auto c = pa.fragment_class[i];
    if (c == gex::ClassTable::NONE) continue;
    const auto tr = fragments[i].transcript;
    truth[tr] += 1;
    if (pa.classes.members(c).size() == 1) own[tr] += 1;
  }
  const auto true_tpm = gex::tpm(truth, eff);
  // Class sizes: how many classes hold 1, 2, 3, ... transcripts, and the fragments in them.
  std::vector<double> class_sizes(9, 0), class_fragments(9, 0);
  for (std::uint32_t c = 0; c < pa.classes.size(); c++) {
    const auto k = std::min<std::size_t>(pa.classes.members(c).size(), 8);
    class_sizes[k]++;
    if (c < pa.counts.size()) class_fragments[k] += pa.counts[c];
  }
  Json rows = Json::array();
  for (std::size_t i = 0; i < T; i++)
    rows.push(Json::array({tx.transcripts[i].name, tx.gene_of[i], tx.transcripts[i].sequence.size(), detail::rounded(truth[i], 1), detail::rounded(ab.counts[i], 10),
                           detail::rounded(eq[i], 10), detail::rounded(true_tpm[i], 100), detail::rounded(ab.tpm[i], 100), B ? detail::rounded(bs.sd[i], 100) : 0,
                           truth[i] > 0 ? detail::rounded(own[i] / truth[i], 1000) : 0, Json(std::vector<double>(sim.exons[i].begin(), sim.exons[i].end()))}));
  return Json::object({{"transcripts", rows}, {"genes", tx.genes.size()}, {"fragments", N}, {"aligned", pa.aligned}, {"classes", pa.classes.size()},
                       {"kmers", index.kmers()}, {"index_bytes", index.memory_bytes()}, {"class_sizes", Json(class_sizes)}, {"class_fragments", Json(class_fragments)},
                       {"rounds", ab.rounds}, {"converged", ab.converged}, {"bootstrap", B},
                       {"em", Json::object({{"correlation", detail::rounded(detail::log_correlation(ab.counts, truth))}, {"median_error", detail::rounded(detail::median_relative_error(ab.counts, truth))}})},
                       {"split", Json::object({{"correlation", detail::rounded(detail::log_correlation(eq, truth))}, {"median_error", detail::rounded(detail::median_relative_error(eq, truth))}})},
                       {"seconds", Json::object({{"simulate", simulate_seconds}, {"index", index_seconds}, {"pseudoalign", pseudo_seconds}, {"em", em_seconds}, {"bootstrap", bootstrap_seconds}})},
                       {"threads", worker.pool.size()}});
}

// A study of two conditions: negative binomial counts of genes in n samples each (base means
// spread over orders of magnitude, dispersions a₀ + a₁/μ, a fraction of genes changed by a fold),
// tested by the domain's DESeq2-like analysis; every call against which genes changed.
ATLAS_OP("expression.differential") {
  const auto G = static_cast<std::size_t>(bounded(request, "genes", 4000, 100, 30000));
  const auto n = static_cast<std::size_t>(bounded(request, "samples", 4, 2, 20));
  const double changed = num(request, "changed", 0.1), fold = num(request, "fold", 1.0);
  const double a0 = num(request, "a0", 0.05), a1 = num(request, "a1", 3.0);
  if (!(changed >= 0 && changed <= 0.5) || !(fold > 0 && fold <= 5)) throw Fault("a fraction changed up to 0.5 and a log2 fold up to 5");
  const auto seed = static_cast<std::uint64_t>(integer(request, "seed", 1));
  // The study's inputs, drawn here: what the experiment would be, not how it is analysed.
  std::vector<double> base(G), dispersion(G), log_fold(G, 0.0), size(2 * n);
  std::vector<int> group(2 * n);
  for (std::size_t g = 0; g < G; g++) {
    calx::random::Stream r(seed, 0x53545544, static_cast<std::uint32_t>(g), 0);  // "STUD"
    base[g] = std::exp(r.normal(4.5, 1.6));
    dispersion[g] = std::exp(std::log(a0 + a1 / base[g]) + r.normal(0, 0.4));
    if (r.uniform() < changed) log_fold[g] = (r.uniform() < 0.5 ? -1 : 1) * fold * std::log(2.0);
  }
  for (std::size_t j = 0; j < 2 * n; j++) {
    calx::random::Stream r(seed, 0x53495a45, static_cast<std::uint32_t>(j), 0);  // "SIZE"
    group[j] = j < n ? 0 : 1;
    size[j] = std::exp(r.normal(0, 0.25));
  }
  const auto counts = gex::simulate_counts(base, dispersion, log_fold, group, size, seed + 1);
  const auto t0 = std::chrono::steady_clock::now();
  const auto res = gex::differential(counts, gex::Design::two_groups(group), gex::DifferentialOptions{}, worker.pool);
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::size_t discoveries = 0, true_discoveries = 0, de = 0;
  Json genes = Json::array();
  for (std::size_t g = 0; g < G; g++) {
    const auto& r = res.genes[g];
    const bool is_de = log_fold[g] != 0, called = std::isfinite(r.padj) && r.padj < 0.1;
    de += is_de;
    discoveries += called;
    true_discoveries += called && is_de;
    genes.push(Json::array({detail::rounded(r.base_mean, 100), detail::rounded(r.log2_fold_change), std::isfinite(r.p_value) ? Json(r.p_value) : Json(),
                            std::isfinite(r.padj) ? Json(r.padj) : Json(), detail::rounded(r.dispersion_gene, 1e5), detail::rounded(r.dispersion_trend, 1e5),
                            detail::rounded(r.dispersion, 1e5), r.dispersion_outlier, detail::rounded(log_fold[g] / std::log(2.0), 1e3), detail::rounded(dispersion[g], 1e5)}));
  }
  Json sf = Json::array();
  for (std::size_t j = 0; j < 2 * n; j++) sf.push(Json::array({detail::rounded(size[j]), detail::rounded(res.size_factors[j])}));
  return Json::object({{"genes", genes}, {"samples", 2 * n}, {"de", de}, {"discoveries", discoveries}, {"true_discoveries", true_discoveries},
                       {"fdr", discoveries ? static_cast<double>(discoveries - true_discoveries) / static_cast<double>(discoveries) : 0},
                       {"power", de ? static_cast<double>(true_discoveries) / static_cast<double>(de) : 0},
                       {"trend", Json::object({{"a0", res.trend_a0}, {"a1", res.trend_a1}, {"parametric", res.trend_parametric}, {"true_a0", a0}, {"true_a1", a1}})},
                       {"prior_variance", res.prior_variance}, {"size_factors", sf}, {"seconds", seconds}, {"threads", worker.pool.size()}});
}

}  // namespace atlas
