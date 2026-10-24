// Checks of the multiple aligners against the truth of simulated families: profile alignment of
// two single sequences equal to Gotoh's global alignment (score and path); families of DNA
// evolved by simulate::evolve_family at substitution rates 0.01 to 0.3 and indel rates 0.003 to
// 0.05 per branch, of 8 to 40 members, and families of protein evolved on BLOSUM62's target
// frequencies, aligned by the progressive method (MUSCLE's three stages, and its first alone),
// the centre star and partial orders; their Q and TC against the true alignment and their SP
// scores against the truth's, reported with the truth's pairwise identities; every row without
// its gaps its sequence; the progressive alignment accurate at low divergence and at least as
// accurate as the centre star wherever the families are not saturated; its SP score the sum of
// its stages' gains; the same alignments on 1 and 8 threads; Gusfield's bound (2 − 2/k) on the
// centre star under edit distance, held against the lower bound Σ D(i, j) of the optimum.

#include <genomics/align/gotoh.hpp>
#include <genomics/align/scoring.hpp>
#include <genomics/msa/centre_star.hpp>
#include <genomics/msa/partial_order.hpp>
#include <genomics/msa/progressive.hpp>
#include <genomics/msa/scores.hpp>
#include <genomics/simulate/family.hpp>
#include <genomics/simulate/genome.hpp>

#include "families.hpp"

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace genomics::msa;
namespace align = genomics::align;
namespace seq = genomics::sequence;
namespace sim = genomics::simulate;

namespace {

struct Case {
  std::size_t members;
  double substitution, indel;
  std::size_t length;
};

struct Scored {
  Accuracy acc;
  std::int64_t sp = 0;
  bool rows = false;
};

Scored score(const Msa& m, const sim::Family& f, const align::Scoring& s, calx::scheduler::ThreadPool& pool) {
  Scored out;
  out.rows = m.size() == f.sequences.size() && m.rectangular();
  for (std::size_t r = 0; out.rows && r < m.size(); r++) out.rows = ungapped(m.rows[r]) == f.sequences[r];
  if (!out.rows) return out;
  out.acc = accuracy(m, test_families::truth(f));
  out.sp = sp_score(m, s, false, pool);
  return out;
}

// The mean and least identity of the pairs of rows of the true alignment.
std::pair<double, double> identities(const Msa& truth, int letters) {
  double mean = 0, least = 1;
  std::size_t pairs = 0;
  for (std::size_t i = 0; i < truth.size(); i++)
    for (std::size_t j = i + 1; j < truth.size(); j++, pairs++) {
      const double id = 1 - p_distance(truth.rows[i], truth.rows[j], letters);
      mean += id;
      least = std::min(least, id);
    }
  return {pairs ? mean / static_cast<double>(pairs) : 1.0, least};
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x50524f47, 0, 0);
  calx::scheduler::ThreadPool one(1), eight(8);

  // ---- two sequences: the profile alignment is Gotoh's ------------------------------------------------------------------------------------
  {
    bool same = true, path = true;
    for (int t = 0; t < 300; t++) {
      const bool protein = t % 3 == 0;
      const align::Scoring sc = protein ? align::blastp() : align::scoring(align::Matrix::dna(2, -3, -1), static_cast<int>(r.below(6)), 1 + static_cast<int>(r.below(3)));
      Codes x(r.below(60)), y(r.below(60));
      for (auto& c : x) c = static_cast<Code>(r.below(protein ? 20 : 4));
      for (auto& c : y) c = static_cast<Code>(r.below(protein ? 20 : 4));
      const ProfileAlignment pa = align_profiles(make_profile(single(x), {}, sc.matrix.size), make_profile(single(y), {}, sc.matrix.size), sc, true);
      const std::int64_t g = align::gotoh_score(x, y, sc, align::Mode::global);
      same = same && std::fabs(pa.score - static_cast<double>(g)) < 1e-9;
      const Msa m = merge(single(x), single(y), pa.path);
      path = path && ungapped(m.rows[0]) == x && ungapped(m.rows[1]) == y && pair_score(m.rows[0], m.rows[1], sc, true) == g;
    }
    s.require("the profile alignment of two single sequences scores as Gotoh's global optimum (300 DNA and protein pairs, gap opens 0 to 5)", same);
    s.require("and its path, merged, is an alignment of the two at that score", path);
  }

  // ---- DNA families ------------------------------------------------------------------------------------
  // The rates are per branch of a balanced tree, so the leaves of n members are up to 2 log2 n
  // branches apart: the table gives the true alignment's mean and least pairwise identity.
  {
    const std::vector<Case> cases = {{8, 0.01, 0.003, 300}, {8, 0.03, 0.01, 300},  {8, 0.05, 0.01, 300},  {16, 0.03, 0.005, 300}, {32, 0.01, 0.003, 300},
                                     {32, 0.03, 0.01, 300}, {32, 0.06, 0.003, 300}, {40, 0.06, 0.01, 250}, {8, 0.1, 0.02, 300},   {16, 0.2, 0.03, 250},
                                     {8, 0.3, 0.05, 250}};
    const ProgressiveOptions full = default_options(Alphabet::dna);
    ProgressiveOptions first = full;
    first.second_stage = false;
    first.refine_rounds = 0;
    std::printf("DNA families (+2/-1, open 4, extend 1): Q / TC against the truth; SP minus the truth's SP, per pair of rows\n");
    std::printf("  n   sub   indel  identity    | progressive    first stage    centre star    partial order | dSP prog     cs     poa | ms\n");
    bool rows = true, not_worse = true, bookkeeping = true, stages = true;
    double q_low = 1, tc_low = 1;
    for (std::size_t c = 0; c < cases.size(); c++) {
      const Case& k = cases[c];
      const std::uint64_t seed = 1000 + c;
      const sim::Family f = sim::evolve_family(sim::random_sequence(k.length, 0.5, seed), k.members, k.substitution, k.indel, seed);
      const Msa truth = test_families::truth(f);
      const auto [mean, least] = identities(truth, 4);
      const auto t0 = std::chrono::steady_clock::now();
      const ProgressiveResult p = progressive(f.sequences, f.names, Alphabet::dna, full, eight);
      const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
      const ProgressiveResult p1 = progressive(f.sequences, f.names, Alphabet::dna, first, eight);
      const CentreStar cs = centre_star(f.sequences, f.names, full.scoring, eight);
      const Msa poa = partial_order_alignment(f.sequences, f.names, full.scoring, p.tree.leaves_under(p.tree.root()));
      const Scored a = score(p.msa, f, full.scoring, eight), a1 = score(p1.msa, f, full.scoring, eight), b = score(cs.msa, f, full.scoring, eight),
                   d = score(poa, f, full.scoring, eight);
      const double truth_sp = static_cast<double>(sp_score(truth, full.scoring, false, eight));
      const double np = static_cast<double>(truth.size() * (truth.size() - 1) / 2);
      std::printf("%3zu  %.2f  %.3f  %.3f %.3f | %.3f / %.3f  %.3f / %.3f  %.3f / %.3f  %.3f / %.3f | %+7.1f %+7.1f %+7.1f | %.0f\n", k.members, k.substitution, k.indel, mean,
                  least, a.acc.q, a.acc.tc, a1.acc.q, a1.acc.tc, b.acc.q, b.acc.tc, d.acc.q, d.acc.tc, (static_cast<double>(a.sp) - truth_sp) / np,
                  (static_cast<double>(b.sp) - truth_sp) / np, (static_cast<double>(d.sp) - truth_sp) / np, ms);
      rows = rows && a.rows && a1.rows && b.rows && d.rows;
      if (mean >= 0.5) not_worse = not_worse && a.acc.q >= b.acc.q - 0.005;
      bookkeeping = bookkeeping && p.sp == sp_score(p.msa, full.scoring, false, one) && p.sp >= p.sp_first;
      stages = stages && a.sp >= a1.sp;
      if (k.members <= 16 && k.substitution <= 0.05 && k.indel <= 0.01) {
        q_low = std::min(q_low, a.acc.q);
        tc_low = std::min(tc_low, a.acc.tc);
      }
    }
    s.require("every row of every alignment (progressive, first stage, centre star, partial order) without its gaps is its sequence", rows);
    s.check("on the four families of at most 16 members at substitution ≤ 0.05 and indels ≤ 0.01 per branch, the progressive alignment's least Q: 1 − Q", 1 - q_low, 0.05);
    s.check("and its least TC: 1 − TC", 1 - tc_low, 0.25);
    s.require("the progressive alignment's Q is never below the centre star's (less 0.005) on the families of mean identity 0.5 or more", not_worse);
    s.require("the SP score the refinement reports is the alignment's, and never below the first stage's", bookkeeping);
    s.require("the full progressive alignment's SP score is at least the first stage's on every family", stages);
  }

  // ---- protein families ------------------------------------------------------------------------------------
  {
    const test_families::ProteinModel model;
    const std::vector<Case> cases = {{8, 0.03, 0.01, 200},  {8, 0.08, 0.01, 200},   {8, 0.15, 0.01, 200}, {32, 0.03, 0.01, 200},
                                     {32, 0.08, 0.003, 200}, {32, 0.15, 0.003, 200}, {24, 0.3, 0.02, 200}};
    const ProgressiveOptions full = default_options(Alphabet::protein);
    std::printf("protein families (BLOSUM62, open 6, extend 1): Q / TC\n");
    std::printf("  n   sub   indel  identity    | progressive    centre star    partial order | dSP prog     cs     poa\n");
    bool rows = true, not_worse = true;
    double q_low = 1;
    for (std::size_t c = 0; c < cases.size(); c++) {
      const Case& k = cases[c];
      const std::uint64_t seed = 2000 + c;
      const sim::Family f = test_families::evolve_protein_family(model, test_families::random_protein(model, k.length, seed), k.members, k.substitution, k.indel, seed);
      const Msa truth = test_families::truth(f);
      const auto [mean, least] = identities(truth, 20);
      const ProgressiveResult p = progressive(f.sequences, f.names, Alphabet::protein, full, eight);
      const CentreStar cs = centre_star(f.sequences, f.names, full.scoring, eight);
      const Msa poa = partial_order_alignment(f.sequences, f.names, full.scoring, p.tree.leaves_under(p.tree.root()));
      const Scored a = score(p.msa, f, full.scoring, eight), b = score(cs.msa, f, full.scoring, eight), d = score(poa, f, full.scoring, eight);
      const double truth_sp = static_cast<double>(sp_score(truth, full.scoring, false, eight));
      const double np = static_cast<double>(truth.size() * (truth.size() - 1) / 2);
      std::printf("%3zu  %.2f  %.3f  %.3f %.3f | %.3f / %.3f  %.3f / %.3f  %.3f / %.3f | %+7.1f %+7.1f %+7.1f\n", k.members, k.substitution, k.indel, mean, least, a.acc.q,
                  a.acc.tc, b.acc.q, b.acc.tc, d.acc.q, d.acc.tc, (static_cast<double>(a.sp) - truth_sp) / np, (static_cast<double>(b.sp) - truth_sp) / np,
                  (static_cast<double>(d.sp) - truth_sp) / np);
      rows = rows && a.rows && b.rows && d.rows;
      if (mean >= 0.25) not_worse = not_worse && a.acc.q >= b.acc.q - 0.005;
      if (k.substitution <= 0.08) q_low = std::min(q_low, a.acc.q);
    }
    s.require("every row of every protein alignment without its gaps is its sequence", rows);
    s.check("at substitution ≤ 0.08 per branch the progressive alignment of protein recovers the true pairs: 1 − least Q", 1 - q_low, 0.05);
    s.require("the progressive alignment's Q is never below the centre star's (less 0.005) on the protein families of mean identity 0.25 or more", not_worse);
  }

  // ---- threads ------------------------------------------------------------------------------------
  {
    const std::uint64_t seed = 77;
    const sim::Family f = sim::evolve_family(sim::random_sequence(220, 0.5, seed), 33, 0.03, 0.01, seed);
    const ProgressiveOptions o = default_options(Alphabet::dna);
    const ProgressiveResult a = progressive(f.sequences, f.names, Alphabet::dna, o, one), b = progressive(f.sequences, f.names, Alphabet::dna, o, eight);
    s.require("the progressive alignment of 33 sequences is the same on 1 and 8 threads, row for row", a.msa.rows == b.msa.rows && a.sp == b.sp &&
                                                                                                         a.refinement.accepted == b.refinement.accepted);
    const CentreStar c1 = centre_star(f.sequences, f.names, o.scoring, one), c8 = centre_star(f.sequences, f.names, o.scoring, eight);
    s.require("the centre star is the same on 1 and 8 threads", c1.msa.rows == c8.msa.rows && c1.centre == c8.centre);
  }

  // ---- Gusfield's bound ------------------------------------------------------------------------------------
  {
    // Edit distance as a score: matches 0, mismatches −1, a gap of length L costing L.
    const align::Scoring edit = align::scoring(align::Matrix::dna(0, -1, -1), 0, 1);
    bool bound = true, exact = true;
    double worst = 0;
    for (int t = 0; t < 30; t++) {
      const std::size_t k = 3 + r.below(10);
      const sim::Family f = sim::evolve_family(sim::random_sequence(40 + r.below(60), 0.5, 500 + t), k, 0.1 + 0.3 * r.uniform(), 0.05 * r.uniform(), 500 + t);
      const CentreStar cs = centre_star(f.sequences, f.names, edit, eight);
      std::int64_t lower = 0;  // Σ_{i<j} D(i, j) ≤ the optimal SP cost
      for (std::size_t i = 0; i < k; i++)
        for (std::size_t j = i + 1; j < k; j++) lower += -cs.scores[i * k + j];
      const std::int64_t cost = -sp_score(cs.msa, edit, true, one);
      bound = bound && static_cast<double>(cost) <= (2.0 - 2.0 / static_cast<double>(k)) * static_cast<double>(lower) + 1e-9;
      if (lower > 0) worst = std::max(worst, static_cast<double>(cost) / static_cast<double>(lower));
      for (std::size_t j = 0; j < k; j++)
        if (j != cs.centre) exact = exact && pair_score(cs.msa.rows[cs.centre], cs.msa.rows[j], edit, true) == cs.scores[cs.centre * k + j];
    }
    s.require("the centre star induces each optimal alignment to the centre exactly (30 families)", exact);
    s.require("its SP cost under edit distance is within (2 − 2/k) of the lower bound Σ D(i, j) on the optimum", bound, "largest ratio " + std::to_string(worst));
  }
  return s.done();
}
