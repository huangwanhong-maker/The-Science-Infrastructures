// Checks against biology, on Escherichia coli's tRNA genes (the data cache's NC_000913.3, read
// with genomics/sequence's GenBank reader; skipped without it). The distinct genes are split
// alternately into a training and a held-out half. The grammar (reversible shape, canonical pairs,
// no lonely pairs) is trained by EM from a start that weighs every canonical pair alike, and its
// structures of the held-out tRNAs are held to what is measured of them without any typed
// reference structure: the acceptor stem, whose seven pairs join the first seven bases to the
// seven before the 3′ discriminator base and CCA, is compared with the same positions in
// shuffles of the same composition, as are the grammar's likelihoods. Then the equivariant EM:
// the grammar tied along the orbits of reversal ρ, of the complement κ and of the reverse
// complement κρ (on the pair set each needs), each compared with the untied model by held-out
// likelihood and by the likelihood-ratio test on the training half; the verdict is printed, not
// assumed. Last, Knudsen and Hein's grammar as published (every pair, lonely pairs allowed),
// trained the same way, for comparison.

#include <genomics/rna/scfg.hpp>
#include <genomics/rna/symmetry.hpp>
#include <genomics/rna/training.hpp>
#include <genomics/sequence/genbank.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace genomics;
using namespace genomics::rna;
namespace gs = genomics::sequence;

namespace {

std::filesystem::path cache() {
  const char* env = std::getenv("SCIENCE_CACHE");
  return env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE);
}

double seconds_since(std::chrono::steady_clock::time_point t) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count(); }

// The acceptor stem's pairs (k, n − 5 − k), k = 0..6, present in a structure.
int acceptor_pairs(const Structure& y) {
  int k = 0;
  for (std::size_t i = 0; i < 7; i++) k += y.has(i, y.size() - 5 - i);
  return k;
}
double acceptor_probability(const PairProbabilities& bp) {
  double p = 0;
  for (std::size_t i = 0; i < 7; i++) p += bp.at(i, bp.n - 5 - i);
  return p / 7;
}

Codes shuffled(const Codes& x, calx::random::Stream& r) {
  Codes s = x;
  for (std::size_t i = s.size(); i > 1; i--) std::swap(s[i - 1], s[r.below(i)]);
  return s;
}

struct Folds {
  std::size_t viterbi6 = 0, mea6 = 0, viterbi7 = 0, count = 0;
  double probability = 0;
};

void fold(const StructureGrammar& g, const Parameters& theta, const Codes& x, Folds& f) {
  const int v = acceptor_pairs(viterbi(g, theta, x).structure);
  PairProbabilities bp;
  const int m = acceptor_pairs(mea_fold(g, theta, x, 1.0, &bp).structure);
  f.viterbi6 += v >= 6;
  f.viterbi7 += v == 7;
  f.mea6 += m >= 6;
  f.probability += acceptor_probability(bp);
  f.count++;
}

std::string pair_table(const StructureGrammar& g, const Parameters& theta) {
  std::string out;
  for (std::size_t side = 1; side < g.nonterminals(); side++) {
    double mass = 0;
    for (std::size_t r = 0; r < g.rules(); r++)
      if (g.info[r].lhs == side && g.info[r].kind == RuleInfo::pair) mass += theta[r];
    out += (out.empty() ? "" : "; ") + g.grammar.nonterminals[side] + ":";
    for (std::size_t r = 0; r < g.rules(); r++)
      if (g.info[r].lhs == side && g.info[r].kind == RuleInfo::pair) {
        char cell[24];
        std::snprintf(cell, sizeof cell, " %c%c %.3f", "ACGU"[g.info[r].a], "ACGU"[g.info[r].b], theta[r] / mass);
        out += cell;
      }
  }
  return out;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const auto path = cache() / "ncbi-genomes" / "NC_000913.3.gb";
  if (!std::filesystem::exists(path)) {
    s.skip("Escherichia coli is not in the data cache: its tRNAs are not folded and no grammar is trained on them");
    return s.done();
  }
  std::ifstream in(path, std::ios::binary);
  std::stringstream text;
  text << in.rdbuf();
  const auto records = gs::read_genbank(text.str());
  const auto& rec = records.at(0);
  const auto genome = gs::encode_dna(rec.sequence);
  std::size_t genes = 0, with_cca = 0;
  std::vector<Codes> distinct;
  std::set<Codes> seen;
  for (const auto& f : rec.features) {
    if (f.key != "tRNA") continue;
    const auto x = gs::extract(genome, f.location);  // T reads as U: the same code
    if (std::any_of(x.begin(), x.end(), [](Code c) { return c > 3; })) continue;
    genes++;
    if (seen.insert(x).second) {
      distinct.push_back(x);
      with_cca += x.size() > 12 && x[x.size() - 3] == gs::C && x[x.size() - 2] == gs::C && x[x.size() - 1] == gs::A;
    }
  }
  std::vector<Codes> train_set, test_set;
  for (std::size_t i = 0; i < distinct.size(); i++) (i % 2 ? test_set : train_set).push_back(distinct[i]);
  std::size_t train_bases = 0, test_bases = 0;
  for (const auto& x : train_set) train_bases += x.size();
  for (const auto& x : test_set) test_bases += x.size();
  char line[400];
  std::snprintf(line, sizeof line, "%zu genes, %zu distinct sequences, %zu ending in CCA; %zu to train on (%zu bases), %zu held out (%zu bases)", genes,
                distinct.size(), with_cca, train_set.size(), train_bases, test_set.size(), test_bases);
  s.require("Escherichia coli's tRNA genes: at least 80, at least 40 distinct, every distinct one ending in CCA", genes >= 80 && distinct.size() >= 40 &&
                                                                                                                       with_cca == distinct.size(),
            line);

  calx::scheduler::ThreadPool pool(std::max(2u, std::thread::hardware_concurrency()));
  TrainOptions o;
  o.max_iterations = 300;
  const auto& G = structure_grammar(Shape::reversible, false, canonical_pairs());

  // ---- the trained grammar and its structures ----------------------------------------------------------------
  auto clock = std::chrono::steady_clock::now();
  const auto fit = train(G, untied(G), starting_parameters(G), train_set, o, &pool);
  const double held = log_likelihood(G, fit.theta, test_set, &pool);
  std::snprintf(line, sizeof line, "%s: log-likelihood %.1f → %.1f in %d iterations (%.1f s); held out %.1f nats, %.3f bits a base", G.name.c_str(),
                fit.log_likelihood.front(), fit.log_likelihood.back(), fit.iterations, seconds_since(clock), held,
                -held / static_cast<double>(test_bases) / std::log(2.0));
  s.require("EM on the training half never decreases the likelihood and converges", fit.monotone && fit.converged, line);
  // Base composition alone, from the training half, for comparison.
  double f[4] = {0, 0, 0, 0}, composition = 0;
  for (const auto& x : train_set)
    for (Code c : x) f[c] += 1.0 / static_cast<double>(train_bases);
  for (const auto& x : test_set)
    for (Code c : x) composition += std::log(f[c]);
  calx::random::Stream rng(2026, 104, 0, 0);
  Folds real, shuffles;
  std::size_t preferred = 0, comparisons = 0;
  double log_odds = 0;
  for (const auto& x : test_set) {
    fold(G, fit.theta, x, real);
    const double lx = log_likelihood(G, fit.theta, x);
    for (int k = 0; k < 3; k++) {
      const Codes y = shuffled(x, rng);
      fold(G, fit.theta, y, shuffles);
      const double ly = log_likelihood(G, fit.theta, y);
      preferred += lx > ly;
      log_odds += lx - ly;
      comparisons++;
    }
  }
  const auto frac = [](std::size_t a, std::size_t b) { return static_cast<double>(a) / static_cast<double>(b); };
  std::snprintf(line, sizeof line, "%zu of %zu, mean log-odds %.1f nats; base composition alone gives the held-out half %.1f nats against the grammar's %.1f",
                preferred, comparisons, log_odds / static_cast<double>(comparisons), composition, held);
  s.require("the trained grammar finds each held-out tRNA more probable than its shuffles (in at least 90% of them)", frac(preferred, comparisons) >= 0.9, line);
  std::snprintf(line, sizeof line, "Viterbi: %zu of %zu tRNAs with at least 6 of the 7 pairs (%zu with all 7), shuffles %zu of %zu; MEA (γ = 1): %zu and %zu",
                real.viterbi6, real.count, real.viterbi7, shuffles.viterbi6, shuffles.count, real.mea6, shuffles.mea6);
  s.require("held-out tRNAs fold with their acceptor stem far more often than shuffles (at least half of them, at most 5% of the shuffles)",
            frac(real.viterbi6, real.count) >= 0.5 && frac(shuffles.viterbi6, shuffles.count) <= 0.05 && frac(real.mea6, real.count) >= 0.5 &&
                frac(shuffles.mea6, shuffles.count) <= 0.05,
            line);
  std::snprintf(line, sizeof line, "%.3f for the tRNAs, %.3f for the shuffles", real.probability / static_cast<double>(real.count),
                shuffles.probability / static_cast<double>(shuffles.count));
  s.require("the mean posterior probability of an acceptor-stem pair is at least 0.4 in the tRNAs and at most 0.1 in the shuffles",
            real.probability / static_cast<double>(real.count) >= 0.4 && shuffles.probability / static_cast<double>(shuffles.count) <= 0.1, line);

  // ---- the equivariant EM and the data's verdict ---------------------------------------------------------------------
  const auto symmetries = pair_symmetries();
  const std::size_t expected_df[3] = {9, 14, 8};
  for (std::size_t k = 0; k < 3; k++) {
    const auto& sym = symmetries[k];
    clock = std::chrono::steady_clock::now();
    const auto& G2 = structure_grammar(Shape::reversible, false, closure(sym, canonical_pairs()));
    const Tie t = orbit_tie(G2, sym);
    const auto tied = train(G2, t, starting_parameters(G2), train_set, o, &pool);
    const auto free = train(G2, untied(G2), tied.theta, train_set, o, &pool);
    double best_free = free.log_likelihood.back();
    if (G2.pairs == G.pairs) best_free = std::max(best_free, fit.log_likelihood.back());
    const auto lr = likelihood_ratio(best_free, free_parameters(untied(G2)), tied.log_likelihood.back(), free_parameters(t));
    const double held_tied = log_likelihood(G2, tied.theta, test_set, &pool), held_free = log_likelihood(G2, free.theta, test_set, &pool);
    std::snprintf(line, sizeof line, "pairs %s; %zu orbits; %zu free parameters against %zu", pair_names(G2.pairs).c_str(), t.orbits(), free_parameters(t),
                  free_parameters(untied(G2)));
    s.require("tie " + sym.name + ": EM tied and then untied from it never decreases the likelihood, the tied fit stays tied, and the tie removes " +
                  std::to_string(expected_df[k]) + " parameters",
              tied.monotone && free.monotone && is_tied(t, tied.theta, 1e-15) && lr.df == expected_df[k], line);
    std::snprintf(line, sizeof line, "training: tied %.1f, untied %.1f, Λ = %.1f on %zu df, p = %.2g; held out: tied %.1f, untied %.1f (%s by %.1f nats); %.1f s",
                  tied.log_likelihood.back(), best_free, lr.statistic, lr.df, lr.p_value, held_tied, held_free, held_free > held_tied ? "untied" : "tied",
                  std::fabs(held_free - held_tied), seconds_since(clock));
    std::string verdict = line;
    if (k == 0) verdict += "; the pairs of each left side, untied: " + pair_table(G2, free.theta) + "; tied: " + pair_table(G2, tied.theta);
    s.require("the data's verdict on " + sym.name, std::isfinite(lr.statistic) && lr.p_value >= 0 && lr.p_value <= 1, verdict);
  }

  // ---- Knudsen and Hein's grammar as published -------------------------------------------------------------------------
  {
    clock = std::chrono::steady_clock::now();
    const auto& kh = knudsen_hein();
    const auto kfit = train(kh, untied(kh), starting_parameters(kh), train_set, o, &pool);
    Folds kreal;
    for (const auto& x : test_set) fold(kh, kfit.theta, x, kreal);
    std::snprintf(line, sizeof line, "log-likelihood %.1f → %.1f in %d iterations; held out %.1f; acceptor stem (≥ 6 pairs) in %zu of %zu Viterbi structures, mean pair probability %.3f; %.1f s",
                  kfit.log_likelihood.front(), kfit.log_likelihood.back(), kfit.iterations, log_likelihood(kh, kfit.theta, test_set, &pool), kreal.viterbi6,
                  kreal.count, kreal.probability / static_cast<double>(kreal.count), seconds_since(clock));
    s.require("Knudsen and Hein's grammar as published, from the same uniform start: EM never decreases the likelihood", kfit.monotone, line);
  }
  return s.done();
}
