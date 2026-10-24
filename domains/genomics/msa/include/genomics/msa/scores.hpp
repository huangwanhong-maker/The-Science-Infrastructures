// Scores of a multiple alignment: its sum-of-pairs score under a scoring, the objective of
// progressive alignment and its refinement; and its accuracy against a reference alignment,
// by the SP (Q) and TC scores of BAliBASE.
//
// - The sum-of-pairs score (Carrillo and Lipman 1988, SIAM J. Appl. Math. 48, 1073–1082;
//   Altschul 1989, J. Theor. Biol. 138, 297–309, for affine gaps) is Σ_{i<j} S(A_i, A_j), where
//   the score of two rows is that of their projection: the pairwise alignment left when the
//   columns in which both rows have a gap are removed, its pairs scored by the matrix and each
//   maximal run of gaps in one row by −(open + L · extend) (the "natural" gap cost; Altschul
//   1989). With terminal_open false, a run that touches either end of the projection pays only
//   its extension (ClustalW's and MUSCLE's leniency to end gaps). O(N² L) for N rows of L
//   columns, in parallel by rows above 2^20 column visits; the scores are integers, so the sum
//   is exact in any order.
// - Accuracy against a reference (Thompson, Plewniak and Poch 1999, Nucleic Acids Research 27,
//   2682–2690, §Scoring): two residues are aligned when they share a column. The SP score of
//   BAliBASE (the Q score of Edgar 2004) is the fraction of the pairs of residues aligned in the
//   reference that the test aligns too; its precision ("modeler", Sauder, Arthur and Dunbrack
//   2000, Proteins 40, 6–22) the fraction of the test's aligned pairs that the reference
//   aligns. The TC score is the fraction of the reference's columns reproduced exactly: a
//   column whose residues all share one test column, and which holds no others. Columns of one
//   residue align nothing and are not counted (as in Edgar's qscore). The whole alignment is
//   scored (BAliBASE scores its core blocks only; a simulated truth has no such blocks). Each
//   residue is placed in both alignments once, O(N L log N).

#pragma once

#include <genomics/align/types.hpp>
#include <genomics/msa/msa.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace genomics::msa {

// The score of the projection of two aligned rows of equal length.
inline std::int64_t pair_score(const Codes& a, const Codes& b, const align::Scoring& s, bool terminal_open = true) {
  if (a.size() != b.size()) throw std::invalid_argument("genomics::msa::pair_score: rows of different lengths");
  std::int64_t score = 0;
  int state = 0;               // 0 a pair, 1 a run of gaps in a, 2 in b
  bool started = false, first_is_gap = false;
  std::size_t runs = 0;
  const std::int64_t open = s.gap_open, extend = s.gap_extend;
  for (std::size_t c = 0; c < a.size(); c++) {
    const bool ga = a[c] == GAP, gb = b[c] == GAP;
    if (ga && gb) continue;
    if (!ga && !gb) {
      score += s.matrix(a[c], b[c]);
      state = 0;
    } else {
      const int now = ga ? 1 : 2;
      if (now != state) { score -= open; runs++; }
      score -= extend;
      if (!started) first_is_gap = true;
      state = now;
    }
    started = true;
  }
  if (!terminal_open && runs) {
    // Refund the openings of the runs that touch an end (one run, if it spans the whole projection).
    const bool last_is_gap = state != 0;
    std::int64_t refund = (first_is_gap ? 1 : 0) + (last_is_gap ? 1 : 0);
    if (runs == 1 && first_is_gap && last_is_gap) refund = 1;
    score += refund * open;
  }
  return score;
}

// Σ_{i<j} of the pair scores.
namespace detail {

// Below about 2^20 column visits a parallel loop costs more in wake-ups than it saves.
inline calx::scheduler::ThreadPool& runner_for(std::size_t work, calx::scheduler::ThreadPool& pool, calx::scheduler::ThreadPool& serial) {
  return work >= (std::size_t{1} << 20) ? pool : serial;
}

}  // namespace detail

inline std::int64_t sp_score(const Msa& m, const align::Scoring& s, bool terminal_open, calx::scheduler::ThreadPool& pool) {
  require_rectangular(m, "genomics::msa::sp_score");
  const std::size_t n = m.size();
  std::vector<std::int64_t> row(n, 0);
  calx::scheduler::ThreadPool serial(1);
  calx::scheduler::parallel_for(detail::runner_for(n * n / 2 * m.columns(), pool, serial), n, 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t i = b; i < e; i++)
      for (std::size_t j = i + 1; j < n; j++) row[i] += pair_score(m.rows[i], m.rows[j], s, terminal_open);
  });
  std::int64_t total = 0;
  for (std::int64_t x : row) total += x;
  return total;
}
inline std::int64_t sp_score(const Msa& m, const align::Scoring& s, bool terminal_open = true) {
  calx::scheduler::ThreadPool one(1);
  return sp_score(m, s, terminal_open, one);
}

// Σ over the pairs of a row of `left` and a row of `right`: the part of the SP score that
// re-aligning the two sides of a split can change.
inline std::int64_t sp_across(const Msa& m, const std::vector<std::size_t>& left, const std::vector<std::size_t>& right, const align::Scoring& s, bool terminal_open,
                              calx::scheduler::ThreadPool& pool) {
  std::vector<std::int64_t> row(left.size(), 0);
  calx::scheduler::ThreadPool serial(1);
  calx::scheduler::parallel_for(detail::runner_for(left.size() * right.size() * m.columns(), pool, serial), left.size(), 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t i = b; i < e; i++)
      for (std::size_t j : right) row[i] += pair_score(m.rows[left[i]], m.rows[j], s, terminal_open);
  });
  std::int64_t total = 0;
  for (std::int64_t x : row) total += x;
  return total;
}

// ---- accuracy against a reference ------------------------------------------------------------------------------------
struct Accuracy {
  double q = 0;           // SP score of BAliBASE: reference pairs recovered over reference pairs
  double tc = 0;          // reference columns reproduced over reference columns (of two residues or more)
  double precision = 0;   // test pairs the reference aligns over test pairs
  std::size_t pairs_correct = 0, pairs_reference = 0, pairs_test = 0;
  std::size_t columns_correct = 0, columns_reference = 0;
};

namespace detail {

// For each row, the column of each of its residues.
inline std::vector<std::vector<std::uint32_t>> residue_columns(const Msa& m) {
  std::vector<std::vector<std::uint32_t>> out(m.size());
  for (std::size_t r = 0; r < m.size(); r++)
    for (std::size_t c = 0; c < m.rows[r].size(); c++)
      if (m.rows[r][c] != GAP) out[r].push_back(static_cast<std::uint32_t>(c));
  return out;
}

inline std::size_t pairs_of(std::size_t k) { return k ? k * (k - 1) / 2 : 0; }

}  // namespace detail

// The test against the reference, row by row: row r of each must be the same sequence.
inline Accuracy accuracy(const Msa& test, const Msa& reference) {
  require_rectangular(test, "genomics::msa::accuracy");
  require_rectangular(reference, "genomics::msa::accuracy");
  if (test.size() != reference.size()) throw std::invalid_argument("genomics::msa::accuracy: the alignments have different numbers of rows");
  for (std::size_t r = 0; r < test.size(); r++)
    if (ungapped(test.rows[r]) != ungapped(reference.rows[r]))
      throw std::invalid_argument("genomics::msa::accuracy: row " + std::to_string(r) + " is not the same sequence in both alignments");
  const std::size_t N = test.size();
  const auto ref_col = detail::residue_columns(reference), test_col = detail::residue_columns(test);
  Accuracy a;
  // Pairs: the residues of each test column, grouped by their reference column.
  std::vector<std::size_t> next(N, 0);
  std::vector<std::uint32_t> group;
  std::vector<std::size_t> test_count(test.columns(), 0);
  for (std::size_t t = 0; t < test.columns(); t++) {
    group.clear();
    for (std::size_t r = 0; r < N; r++)
      if (test.rows[r][t] != GAP) group.push_back(ref_col[r][next[r]++]);
    test_count[t] = group.size();
    a.pairs_test += detail::pairs_of(group.size());
    std::sort(group.begin(), group.end());
    for (std::size_t i = 0; i < group.size();) {
      std::size_t j = i;
      while (j < group.size() && group[j] == group[i]) j++;
      a.pairs_correct += detail::pairs_of(j - i);
      i = j;
    }
  }
  // Columns: a reference column is reproduced when its residues share one test column of
  // exactly as many residues.
  std::fill(next.begin(), next.end(), 0);
  for (std::size_t c = 0; c < reference.columns(); c++) {
    std::size_t count = 0;
    std::uint32_t where = 0;
    bool together = true;
    for (std::size_t r = 0; r < N; r++) {
      if (reference.rows[r][c] == GAP) continue;
      const std::uint32_t t = test_col[r][next[r]++];
      if (count == 0) where = t;
      else if (t != where) together = false;
      count++;
    }
    a.pairs_reference += detail::pairs_of(count);
    if (count < 2) continue;
    a.columns_reference++;
    if (together && test_count[where] == count) a.columns_correct++;
  }
  a.q = a.pairs_reference ? static_cast<double>(a.pairs_correct) / static_cast<double>(a.pairs_reference) : 1.0;
  a.precision = a.pairs_test ? static_cast<double>(a.pairs_correct) / static_cast<double>(a.pairs_test) : 1.0;
  a.tc = a.columns_reference ? static_cast<double>(a.columns_correct) / static_cast<double>(a.columns_reference) : 1.0;
  return a;
}

}  // namespace genomics::msa
