// Checks of the scores of a multiple alignment: the score of two rows against hand-worked
// projections (two runs of gaps in different rows, a terminal run with and without its opening,
// a row of gaps only); the pair score of random rows equal to align::rescore of their
// projection as a pairwise alignment; the SP score the sum over pairs, the same on 1 and 8
// threads, and its part across a split; the Q, TC and precision of a hand-worked test against
// its reference (Q = 6/8, TC = 3/4, precision 1, BAliBASE's definitions in Thompson, Plewniak
// and Poch 1999), of an alignment against itself (all 1), and the duality Q(A, B) = precision(B, A)
// on random pairs of alignments of the same sequences.

#include <genomics/align/scoring.hpp>
#include <genomics/align/types.hpp>
#include <genomics/msa/msa.hpp>
#include <genomics/msa/scores.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace genomics::msa;
namespace align = genomics::align;
namespace seq = genomics::sequence;

namespace {

Codes row(std::string_view letters) {
  Codes out;
  for (char c : letters) out.push_back(code(c, Alphabet::dna));
  return out;
}
Msa from_letters(const std::vector<std::string>& rows) {
  Msa m;
  for (std::size_t r = 0; r < rows.size(); r++) {
    m.names.push_back("r" + std::to_string(r + 1));
    m.rows.push_back(row(rows[r]));
  }
  return m;
}

// The projection of two rows as a pairwise alignment of their sequences: a the query, b the target.
align::Alignment projection(const Codes& a, const Codes& b) {
  align::Alignment x;
  for (std::size_t c = 0; c < a.size(); c++) {
    const bool ga = a[c] == GAP, gb = b[c] == GAP;
    if (ga && gb) continue;
    align::push(x.cigar, ga ? 'D' : gb ? 'I' : 'M');
    if (!ga) x.query_end++;
    if (!gb) x.target_end++;
  }
  return x;
}

// A random alignment of random sequences: residues placed in order among gap columns.
Msa random_alignment(calx::random::Stream& r, const std::vector<Codes>& seqs, std::size_t width) {
  Msa m;
  for (std::size_t i = 0; i < seqs.size(); i++) {
    m.names.push_back("s" + std::to_string(i));
    // Choose which of `width` columns hold the residues.
    Codes out(width, GAP);
    std::size_t need = seqs[i].size(), left = width, k = 0;
    for (std::size_t c = 0; c < width; c++, left--)
      if (need && r.below(left) < need) { out[c] = seqs[i][k++]; need--; }
    m.rows.push_back(out);
  }
  return m;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x53434f52, 0, 0);
  const align::Scoring sc = align::scoring(align::Matrix::dna(2, -4, -1), 4, 2);

  // ---- the pair score by hand ------------------------------------------------------------------------------------
  {
    // AC--GT / A--TGT projects to A/A C/- -/T G/G T/T: 2 + 2 + 2 − 2(4 + 2) = −6.
    s.require("two rows with a run of gaps in each score 3 pairs minus two openings: −6", pair_score(row("AC--GT"), row("A--TGT"), sc) == -6,
              std::to_string(pair_score(row("AC--GT"), row("A--TGT"), sc)));
    // --ACGT / TTACGT: four pairs (8) and a terminal run of 2: −(4 + 4) with its opening, −4 without.
    s.require("a terminal run pays its opening only when terminal gaps are charged: 0 and 4",
              pair_score(row("--ACGT"), row("TTACGT"), sc, true) == 0 && pair_score(row("--ACGT"), row("TTACGT"), sc, false) == 4);
    s.require("a row of gaps only against two bases: one run of 2, −8, or −4 when its opening is free",
              pair_score(row("AC"), row("--"), sc, true) == -8 && pair_score(row("AC"), row("--"), sc, false) == -4);
    s.require("columns of gaps in both rows are not seen: A-C/A-C scores as AC/AC", pair_score(row("A-C"), row("A-C"), sc) == 4);
    // A run at each end, both refunded: --AC-- / GGACTT, 4 − 2·(4 + 4) = −12, or −4 without openings.
    s.require("runs at both ends are each refunded their opening", pair_score(row("--AC--"), row("GGACTT"), sc, false) == -4 &&
                                                                         pair_score(row("--AC--"), row("GGACTT"), sc, true) == -12);
  }

  // ---- the pair score is the score of the projection ------------------------------------------------------------------------------------
  {
    bool same = true;
    for (int t = 0; t < 500 && same; t++) {
      const std::size_t w = 1 + r.below(60);
      Codes a(w), b(w);
      for (std::size_t c = 0; c < w; c++) {
        a[c] = r.uniform() < 0.3 ? GAP : static_cast<Code>(r.below(5));
        b[c] = r.uniform() < 0.3 ? GAP : static_cast<Code>(r.below(5));
      }
      const auto x = projection(a, b);
      same = pair_score(a, b, sc, true) == align::rescore(x, ungapped(a), ungapped(b), sc);
    }
    s.require("the pair score of 500 random pairs of rows is align::rescore of their projection", same);
  }

  // ---- the SP score ------------------------------------------------------------------------------------
  {
    std::vector<Codes> seqs;
    for (int i = 0; i < 17; i++) {
      Codes x(40 + r.below(30));
      for (auto& c : x) c = static_cast<Code>(r.below(4));
      seqs.push_back(x);
    }
    const Msa m = random_alignment(r, seqs, 90);
    std::int64_t by_hand = 0;
    for (std::size_t i = 0; i < m.size(); i++)
      for (std::size_t j = i + 1; j < m.size(); j++) by_hand += pair_score(m.rows[i], m.rows[j], sc, false);
    calx::scheduler::ThreadPool one(1), eight(8);
    s.require("the SP score is the sum of the pair scores, the same on 1 and 8 threads",
              sp_score(m, sc, false, one) == by_hand && sp_score(m, sc, false, eight) == by_hand);
    const std::vector<std::size_t> left = {0, 3, 4, 9, 16}, right = {1, 2, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15};
    const std::int64_t across = sp_across(m, left, right, sc, false, eight);
    const std::int64_t inside = sp_score(sub_alignment(m, left), sc, false) + sp_score(sub_alignment(m, right), sc, false);
    s.require("the SP score is the part across a split plus the SP scores of the two sub-alignments", across + inside == by_hand);
  }

  // ---- accuracy against a reference ------------------------------------------------------------------------------------
  {
    const Msa ref = from_letters({"ACGT", "AC-T", "A-GT"});
    const Msa test = from_letters({"ACGT-", "AC--T", "A-GT-"});
    const Accuracy a = accuracy(test, ref);
    s.require("a hand-worked test recovers 6 of the reference's 8 aligned pairs (Q = 0.75) and 3 of its 4 columns (TC = 0.75), with no wrong pair",
              a.pairs_reference == 8 && a.pairs_correct == 6 && a.pairs_test == 6 && a.columns_reference == 4 && a.columns_correct == 3 && a.q == 0.75 &&
                  a.tc == 0.75 && a.precision == 1.0);
    const Accuracy b = accuracy(ref, test);
    s.require("against the test as reference, the reference recovers every pair (Q = 1) at precision 6/8", b.q == 1.0 && b.precision == 0.75);
    const Accuracy self = accuracy(ref, ref);
    s.require("an alignment against itself has Q = TC = precision = 1", self.q == 1 && self.tc == 1 && self.precision == 1);

    bool dual = true;
    for (int t = 0; t < 100 && dual; t++) {
      std::vector<Codes> seqs;
      const std::size_t n = 2 + r.below(10);
      for (std::size_t i = 0; i < n; i++) {
        Codes x(1 + r.below(25));
        for (auto& c : x) c = static_cast<Code>(r.below(4));
        seqs.push_back(x);
      }
      const Msa x = random_alignment(r, seqs, 30), y = random_alignment(r, seqs, 35);
      const Accuracy p = accuracy(x, y), q = accuracy(y, x);
      dual = p.pairs_correct == q.pairs_correct && std::fabs(p.q - q.precision) < 1e-15 && std::fabs(p.precision - q.q) < 1e-15;
    }
    s.require("Q(A against B) = precision(B against A) on 100 random pairs of alignments of the same sequences", dual);
    bool threw = false;
    try {
      accuracy(from_letters({"ACGT", "AC-T"}), from_letters({"ACGT", "AC-A"}));
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    s.require("alignments of different sequences are refused", threw);
  }
  return s.done();
}
