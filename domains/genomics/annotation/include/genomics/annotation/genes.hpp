// A self-trained finder of protein-coding genes in prokaryotic genomes: the statistics of
// codons, learnt from the genome itself.
//
// Coding sequence has a period of three: the base at each position of a codon follows its own
// distribution, conditioned on the bases before it (Fickett 1982, Nucleic Acids Res. 10,
// 5303–5318; Borodovsky and McIninch 1993, Computers Chem. 17, 123–133). Non-coding sequence
// has none. So a three-periodic inhomogeneous Markov chain of order k (the phase of a base is
// its position in its codon) against a homogeneous chain of the same order tells a reading
// frame that codes from one that does not, and the log-likelihood ratio of an ORF,
// Σ_t log₂ P_c(x_t | x_{t−k..t−1}, phase) − log₂ P_n(x_t | x_{t−k..t−1}), is the evidence that it
// is a gene.
//
// - Training without annotation (the idea of GeneMarkS, Besemer, Lomsadze and Borodovsky 2001,
//   Nucleic Acids Res. 29, 2607–2618, and of Glimmer's long ORFs, Delcher, Harmon, Kasif, White
//   and Salzberg 1999, Nucleic Acids Res. 27, 4636–4641): the first coding model is learnt from
//   the long ORFs (≥ 300 codons by default, from sequence::open_reading_frames, those overlapping
//   a longer one dropped), which are nearly all genes because a random frame meets a stop every
//   21 codons or so; the non-coding model from the rest of the genome on both strands. Then the
//   genome is annotated, and both models are re-learnt from the prediction, for a few rounds or
//   until the prediction no longer changes.
// - The chains have orders up to k (5 by default, the order of GeneMark), estimated with a
//   Dirichlet prior centred on the chain of order one less (P_j(b | c) = (n(c, b) + α
//   P_{j−1}(b | c′))/(n(c) + α), c′ the context without its most distant base; MacKay and Peto
//   1995, Nat. Lang. Eng. 1, 289–307, the hierarchical Dirichlet model), so that a context seen
//   rarely falls back on a shorter one, the purpose Glimmer's interpolation serves (Salzberg,
//   Delcher, Kasif and White 1998, Nucleic Acids Res. 26, 544–548).
// - Candidates are the ORFs of both strands of at least `min_codons` codons, with every in-frame
//   ATG, GTG or TTG (those the genetic code allows) as a possible start. The start of a gene is
//   the one maximising its log-likelihood ratio plus the log prior of its start codon; upstream
//   of the true start the sequence is non-coding, and lengthening there lowers the ratio.
// - Genes are the candidates whose ratio exceeds a threshold, of which a set of largest total
//   ratio with no two sharing more than `max_overlap` bases is chosen, exactly, by weighted
//   interval scheduling on the intervals shrunk by max_overlap/2 at each end (two shrunk
//   intervals are disjoint exactly when the originals overlap by at most max_overlap).
//
// *Caveat.* This is a simplification of GeneMark.hmm and GeneMarkS, which decode a hidden
// semi-Markov model with gene length distributions and learn a ribosome-binding-site motif to
// choose starts, and of Glimmer3; it has neither, so its stops are found far better than its
// starts (the checks report both).
//
// Complexity: per round O(n (k + 1)) to train and O(n) to score (every base of the six frames
// is visited once per candidate it lies in, about twice in all), O(G log G) for the selection
// of G candidates. Memory O(n) for the reverse complement, O(3 · 4^{k+1}) for the tables. On
// E. coli K-12 (4.6 Mb), six rounds take about a second on one core of an i7-11370H, and find
// the stops of 95% of the annotated CDSs of 80 codons or more, 97.5% of the predictions being
// annotated CDSs (the checks).
//
// Limits: genes crossing the origin of a circular genome are not found (the genome is read as
// linear); frameshifts, selenocysteine and other recodings are invisible to it.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/genetic_code.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace genomics::annotation {

using sequence::Code;
using sequence::Codes;

// ---- Markov chains with phases -----------------------------------------------------------------------
// A chain of order k with P phases; the context of a base is the k bases before it, packed two
// bits each with the nearest in the lowest bits, so that dropping the most distant base is a mask.
class PeriodicChain {
 public:
  PeriodicChain(int order = 5, int phases = 3) : k_(order), P_(phases) {
    if (order < 0 || order > 8 || phases < 1) throw std::invalid_argument("genomics::annotation::PeriodicChain: order 0 to 8 and at least one phase");
    counts_.resize(static_cast<std::size_t>(k_) + 1);
    log_.resize(static_cast<std::size_t>(k_) + 1);
    for (int j = 0; j <= k_; j++) counts_[static_cast<std::size_t>(j)].assign(static_cast<std::size_t>(P_) << (2 * (j + 1)), 0);
  }
  int order() const { return k_; }
  int phases() const { return P_; }

  void clear() {
    for (auto& c : counts_) std::fill(c.begin(), c.end(), 0);
  }
  // Counts the bases of x[from, to): the base at t is in phase (t − from + first_phase) mod P,
  // with the context of the bases before it from `context_from` on (an N restarts the context).
  void add(const Codes& x, std::size_t from, std::size_t to, int first_phase = 0, std::size_t context_from = SIZE_MAX) {
    to = std::min(to, x.size());
    if (context_from == SIZE_MAX) context_from = from;
    std::uint32_t ctx = 0;
    int have = 0;
    for (std::size_t t = context_from; t < from; t++) push(ctx, have, x[t]);
    int phase = first_phase % P_;
    for (std::size_t t = from; t < to; t++) {
      const Code b = x[t];
      if (b < 4)
        for (int j = 0; j <= std::min(k_, have); j++) counts_[static_cast<std::size_t>(j)][index(j, phase, ctx & mask(j), b)] += 1;
      push(ctx, have, b);
      if (++phase == P_) phase = 0;
    }
  }
  // The estimates and their log₂ tables, with prior strength alpha.
  void finish(double alpha = 4) {
    std::vector<double> lower;  // P_{j−1}, for each (phase, context of j − 1), 4 values
    for (int j = 0; j <= k_; j++) {
      const auto& n = counts_[static_cast<std::size_t>(j)];
      std::vector<double> p(n.size());
      const std::uint32_t contexts = 1u << (2 * j);
      for (int ph = 0; ph < P_; ph++)
        for (std::uint32_t c = 0; c < contexts; c++) {
          const std::size_t at = index(j, ph, c, 0);
          double total = 0;
          for (int b = 0; b < 4; b++) total += n[at + static_cast<std::size_t>(b)];
          const std::size_t below = j ? index(j - 1, ph, c & ((1u << (2 * (j - 1))) - 1), 0) : 0;
          for (int b = 0; b < 4; b++) {
            const double prior = j ? lower[below + static_cast<std::size_t>(b)] : 0.25;
            p[at + static_cast<std::size_t>(b)] = (n[at + static_cast<std::size_t>(b)] + alpha * prior) / (total + alpha);
          }
        }
      auto& l = log_[static_cast<std::size_t>(j)];
      l.resize(p.size());
      for (std::size_t i = 0; i < p.size(); i++) l[i] = std::log2(p[i]);
      lower = std::move(p);
    }
  }
  // log₂ P(b | context, phase) with a context of `have` known bases (the order used is min(k, have)).
  double log2p(int phase, std::uint32_t ctx, int have, Code b) const {
    const int j = std::min(k_, have);
    return log_[static_cast<std::size_t>(j)][index(j, phase, ctx & mask(j), b)];
  }
  // The table of the full order, P × 4^k × 4, for scanning.
  const std::vector<double>& table() const { return log_[static_cast<std::size_t>(k_)]; }

  std::uint32_t mask(int j) const { return j >= 16 ? ~0u : (1u << (2 * j)) - 1; }
  void push(std::uint32_t& ctx, int& have, Code b) const {
    if (b > 3) { ctx = 0; have = 0; return; }
    ctx = ((ctx << 2) | b) & mask(k_);
    have = std::min(have + 1, k_);
  }
  std::size_t index(int j, int phase, std::uint32_t ctx, Code b) const {
    return ((static_cast<std::size_t>(phase) << (2 * j)) + ctx) * 4 + b;
  }

 private:
  int k_, P_;
  std::vector<std::vector<double>> counts_;  // per order j: P × 4^j × 4
  std::vector<std::vector<double>> log_;
};

// ---- the gene finder -------------------------------------------------------------------------------------
struct GeneOptions {
  int table = 11;                     // NCBI's genetic code
  int order = 5;
  double alpha = 4;                   // the prior strength of the back-off
  std::size_t min_codons = 80;        // candidates: ORFs of at least this many codons, the stop included
  std::size_t train_min_codons = 300; // the long ORFs of the first round
  int rounds = 6;                     // at most; stops earlier when the prediction is unchanged
  double threshold = 0;               // bits
  std::size_t max_overlap = 60;       // bases
  // ATG, GTG, TTG: the frequencies of the starts of E. coli K-12's annotated CDSs (0.901, 0.079,
  // 0.019 in NC_000913.3, measured by the checks).
  std::array<double, 3> start_prior = {0.90, 0.08, 0.02};
};

struct Gene {
  std::size_t begin = 0, end = 0;          // forward-strand coordinates, [begin, end), the stop codon included
  bool reverse = false;
  double score = 0;                        // the log-likelihood ratio, bits
  std::size_t orf_begin = 0, orf_end = 0;  // the longest ORF, from the first start after the previous stop
  int start = 0;                           // the start codon: 0 ATG, 1 GTG, 2 TTG
  std::size_t codons() const { return (end - begin) / 3; }
};

struct GeneReport {
  std::size_t candidates = 0, training_orfs = 0;
  std::vector<std::size_t> genes_per_round;
  int rounds = 0;
  bool converged = false;
};

namespace detail {

struct Candidate {
  bool reverse = false;
  std::size_t stop_end = 0;            // strand coordinates: the end of the stop codon
  std::uint32_t first = 0, count = 0;  // its starts in the flat array, in increasing order
};

// The candidates of one strand y: ORFs with every allowed start, at least min_codons long.
inline void candidates_of(const Codes& y, bool reverse, const sequence::GeneticCode& code, const std::array<bool, 3>& allowed, std::size_t min_codons,
                          std::vector<Candidate>& out, std::vector<std::uint32_t>& starts, std::vector<std::uint8_t>& kinds) {
  const std::size_t n = y.size();
  std::vector<std::uint32_t> open;
  std::vector<std::uint8_t> open_kind;
  for (int frame = 0; frame < 3; frame++) {
    open.clear();
    open_kind.clear();
    for (std::size_t i = static_cast<std::size_t>(frame); i + 3 <= n; i += 3) {
      const Code a = y[i], b = y[i + 1], c = y[i + 2];
      if (sequence::is_stop(code, a, b, c)) {
        const std::size_t end = i + 3;
        Candidate cand;
        cand.reverse = reverse;
        cand.stop_end = end;
        cand.first = static_cast<std::uint32_t>(starts.size());
        for (std::size_t s = 0; s < open.size(); s++)
          if ((end - open[s]) / 3 >= min_codons) {
            starts.push_back(open[s]);
            kinds.push_back(open_kind[s]);
          }
        cand.count = static_cast<std::uint32_t>(starts.size()) - cand.first;
        if (cand.count) out.push_back(cand);
        open.clear();
        open_kind.clear();
        continue;
      }
      if (b == sequence::T && c == sequence::G) {
        const int kind = a == sequence::A ? 0 : a == sequence::G ? 1 : a == sequence::T ? 2 : -1;
        if (kind >= 0 && allowed[static_cast<std::size_t>(kind)]) {
          open.push_back(static_cast<std::uint32_t>(i));
          open_kind.push_back(static_cast<std::uint8_t>(kind));
        }
      }
    }
  }
}

// The log-ratio tables of the full order: coding (3 phases) minus non-coding, bits.
struct Ratio {
  int k = 5;
  std::vector<double> full;  // 3 × 4^k × 4
};

}  // namespace detail

inline std::vector<Gene> find_genes(const Codes& genome, const GeneOptions& o = {}, GeneReport* report = nullptr) {
  const sequence::GeneticCode& code = sequence::genetic_code(o.table);
  const std::size_t n = genome.size();
  const Codes rc = sequence::reverse_complement(genome);
  std::array<bool, 3> allowed{};
  const Code first[3] = {sequence::A, sequence::G, sequence::T};
  for (int s = 0; s < 3; s++) allowed[static_cast<std::size_t>(s)] = s == 0 || sequence::is_start(code, first[s], sequence::T, sequence::G);

  // The candidates of both strands.
  std::vector<detail::Candidate> cands;
  std::vector<std::uint32_t> starts;
  std::vector<std::uint8_t> kinds;
  detail::candidates_of(genome, false, code, allowed, o.min_codons, cands, starts, kinds);
  detail::candidates_of(rc, true, code, allowed, o.min_codons, cands, starts, kinds);
  if (report) report->candidates = cands.size();

  const auto to_forward = [&](bool reverse, std::size_t b, std::size_t e) { return reverse ? std::pair{n - e, n - b} : std::pair{b, e}; };

  // The genes of the first round: the long ORFs, the longest first, none overlapping a kept one by more than max_overlap.
  std::vector<Gene> genes;
  {
    auto orfs = sequence::open_reading_frames(genome, code, o.train_min_codons, false);
    std::sort(orfs.begin(), orfs.end(), [](const sequence::Orf& a, const sequence::Orf& b) { return a.end - a.begin != b.end - b.begin ? a.end - a.begin > b.end - b.begin : a.begin < b.begin; });
    std::vector<std::pair<std::size_t, std::size_t>> kept;
    for (const auto& orf : orfs) {
      bool ok = true;
      for (const auto& [b, e] : kept)
        if (std::min(e, orf.end) > std::max(b, orf.begin) && std::min(e, orf.end) - std::max(b, orf.begin) > o.max_overlap) { ok = false; break; }
      if (!ok) continue;
      kept.push_back({orf.begin, orf.end});
      Gene g;
      g.begin = orf.begin;
      g.end = orf.end;
      g.reverse = orf.reverse;
      genes.push_back(g);
    }
    if (report) report->training_orfs = genes.size();
  }

  PeriodicChain coding(o.order, 3), noncoding(o.order, 1);
  std::vector<double> r;  // per-base ratios of one candidate
  std::vector<Gene> previous;
  int round = 0;
  for (; round < o.rounds; round++) {
    // Train on the current genes and the rest of the genome.
    coding.clear();
    noncoding.clear();
    std::vector<std::uint8_t> covered(n, 0);
    for (const Gene& g : genes) {
      const Codes& y = g.reverse ? rc : genome;
      const auto [b, e] = g.reverse ? std::pair{n - g.end, n - g.begin} : std::pair{g.begin, g.end};
      coding.add(y, b, e - 3, 0, b >= static_cast<std::size_t>(o.order) ? b - static_cast<std::size_t>(o.order) : 0);
      std::fill(covered.begin() + static_cast<std::ptrdiff_t>(g.begin), covered.begin() + static_cast<std::ptrdiff_t>(g.end), 1);
    }
    for (std::size_t i = 0; i < n;) {
      if (covered[i]) { i++; continue; }
      std::size_t j = i;
      while (j < n && !covered[j]) j++;
      noncoding.add(genome, i, j, 0);
      noncoding.add(rc, n - j, n - i, 0);
      i = j;
    }
    coding.finish(o.alpha);
    noncoding.finish(o.alpha);

    // Score every candidate at every start; keep its best start.
    std::vector<Gene> scored;
    scored.reserve(cands.size());
    for (const auto& c : cands) {
      const Codes& y = c.reverse ? rc : genome;
      const std::size_t s0 = starts[c.first], coding_end = c.stop_end - 3;
      r.assign(coding_end - s0, 0);
      std::uint32_t ctx = 0;
      int have = 0;
      const std::size_t from = s0 >= static_cast<std::size_t>(o.order) ? s0 - static_cast<std::size_t>(o.order) : 0;
      for (std::size_t t = from; t < s0; t++) coding.push(ctx, have, y[t]);
      int phase = 0;
      for (std::size_t t = s0; t < coding_end; t++) {
        const Code b = y[t];
        r[t - s0] = b < 4 ? coding.log2p(phase, ctx, have, b) - noncoding.log2p(0, ctx, have, b) : 0;
        coding.push(ctx, have, b);
        if (++phase == 3) phase = 0;
      }
      // Suffix sums from the stop back to each start.
      double best = -INFINITY, best_llr = 0;
      std::size_t best_start = s0;
      int best_kind = 0;
      double acc = 0;
      std::size_t t = coding_end;
      for (std::uint32_t k = c.count; k-- > 0;) {
        const std::size_t s = starts[c.first + k];
        while (t > s) acc += r[--t - s0];
        const int kind = kinds[c.first + k];
        const double v = acc + std::log2(o.start_prior[static_cast<std::size_t>(kind)]);
        if (v > best) { best = v; best_llr = acc; best_start = s; best_kind = kind; }
      }
      if (!(best_llr > o.threshold)) continue;
      Gene g;
      g.reverse = c.reverse;
      std::tie(g.begin, g.end) = to_forward(c.reverse, best_start, c.stop_end);
      std::tie(g.orf_begin, g.orf_end) = to_forward(c.reverse, s0, c.stop_end);
      g.score = best_llr;
      g.start = best_kind;
      scored.push_back(g);
    }

    // The largest total score with overlaps of at most max_overlap: weighted interval scheduling.
    const std::size_t h = o.max_overlap / 2;
    std::sort(scored.begin(), scored.end(), [](const Gene& a, const Gene& b) { return a.end != b.end ? a.end < b.end : a.begin < b.begin; });
    const std::size_t G = scored.size();
    std::vector<std::size_t> lo(G), hi(G);
    for (std::size_t i = 0; i < G; i++) {
      lo[i] = scored[i].begin + h;
      hi[i] = scored[i].end - h;
    }
    std::vector<std::size_t> order(G);
    for (std::size_t i = 0; i < G; i++) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return hi[a] != hi[b] ? hi[a] < hi[b] : a < b; });
    std::vector<std::size_t> ends(G);
    for (std::size_t i = 0; i < G; i++) ends[i] = hi[order[i]];
    std::vector<double> best(G + 1, 0);
    std::vector<std::size_t> pred(G);
    for (std::size_t i = 0; i < G; i++) {
      const std::size_t g = order[i];
      // The number of intervals (in the order) ending at or before lo[g].
      pred[i] = static_cast<std::size_t>(std::upper_bound(ends.begin(), ends.begin() + static_cast<std::ptrdiff_t>(i), lo[g]) - ends.begin());
      const double take = scored[g].score + best[pred[i]];
      best[i + 1] = take > best[i] ? take : best[i];
    }
    std::vector<Gene> chosen;
    for (std::size_t i = G; i > 0;) {
      const std::size_t g = order[i - 1];
      if (best[i] != best[i - 1] && scored[g].score + best[pred[i - 1]] == best[i]) {
        chosen.push_back(scored[g]);
        i = pred[i - 1];
      } else {
        i--;
      }
    }
    std::sort(chosen.begin(), chosen.end(), [](const Gene& a, const Gene& b) { return a.begin != b.begin ? a.begin < b.begin : a.end < b.end; });
    if (report) report->genes_per_round.push_back(chosen.size());
    const bool same = chosen.size() == previous.size() && std::equal(chosen.begin(), chosen.end(), previous.begin(), [](const Gene& a, const Gene& b) {
                        return a.begin == b.begin && a.end == b.end && a.reverse == b.reverse;
                      });
    previous = chosen;
    genes = std::move(chosen);
    if (same) {
      if (report) report->converged = true;
      round++;
      break;
    }
  }
  if (report) report->rounds = round;
  return genes;
}

}  // namespace genomics::annotation
