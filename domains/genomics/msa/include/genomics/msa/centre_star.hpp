// The centre-star alignment (Gusfield 1993, Bulletin of Mathematical Biology 55, 141–154;
// Gusfield 1997, Algorithms on Strings, Trees, and Sequences, §14.6.2): the sequence whose
// optimal global alignments to all the others score best in sum is the centre; every other
// sequence is aligned to it alone, and the pairwise alignments are merged by "once a gap,
// always a gap", so that the multiple alignment induces each of them exactly.
//
// The bound (Gusfield 1993; Gusfield 1997, §14.6.2; the theory note, Proposition 3.3). For a cost D on
// pairwise alignments that satisfies the triangle inequality and vanishes on a sequence with
// itself (edit distance; linear gaps and a metric substitution cost), the sum-of-pairs cost of
// the centre-star alignment of k sequences is at most (2 − 2/k) times the optimum. With affine
// gaps or similarity scores the triangle inequality fails and the bound is not claimed; the
// method is still what it is, the star of optimal pairwise alignments.
//
// The merge. Let c be the centre. For each position p = 0 .. |c| of c, the gap before c_p (p =
// |c| the end) receives in the alignment of x to c some residues of x against gaps, I_x(p) of
// them; the multiple alignment gives that gap max_x I_x(p) columns, and places x's residues
// there left-justified. Every pairwise alignment to c is then the projection of the result
// (columns of gaps in both rows removed), so the SP score of the result counts each pair (c, x)
// at its optimum.
//
// Complexity: the n(n − 1)/2 pairwise scores O(n² L²) time and O(L) memory each (align::gotoh_end),
// in parallel; then n − 1 alignments with traceback, O(n L²); the merge O(n L).

#pragma once

#include <genomics/align/gotoh.hpp>
#include <genomics/align/types.hpp>
#include <genomics/msa/msa.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace genomics::msa {

struct CentreStar {
  Msa msa;
  std::size_t centre = 0;
  std::vector<std::int64_t> scores;  // n × n: the optimal global scores of the pairs (0 on the diagonal)
};

inline CentreStar centre_star(const std::vector<Codes>& seqs, const std::vector<std::string>& names_in, const align::Scoring& s,
                              calx::scheduler::ThreadPool& pool) {
  const std::size_t n = seqs.size();
  if (!names_in.empty() && names_in.size() != n) throw std::invalid_argument("genomics::msa::centre_star: one name per sequence");
  std::vector<std::string> names = names_in;
  if (names.empty())
    for (std::size_t i = 0; i < n; i++) names.push_back("s" + std::to_string(i + 1));
  CentreStar out;
  if (n == 0) return out;
  out.scores.assign(n * n, 0);
  // Every pair, a row of the upper triangle at a time.
  calx::scheduler::parallel_for(pool, n, 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t i = b; i < e; i++)
      for (std::size_t j = i + 1; j < n; j++) {
        const std::int64_t x = align::gotoh_score(seqs[i], seqs[j], s, align::Mode::global);
        out.scores[i * n + j] = out.scores[j * n + i] = x;
      }
  });
  std::int64_t best = 0;
  for (std::size_t c = 0; c < n; c++) {
    std::int64_t sum = 0;
    for (std::size_t j = 0; j < n; j++) sum += out.scores[c * n + j];
    if (c == 0 || sum > best) { best = sum; out.centre = c; }
  }
  const std::size_t c = out.centre;
  const Codes& centre = seqs[c];
  const std::size_t L = centre.size();
  // Each sequence aligned to the centre: the query x, the target c; I a residue of x against a
  // gap in c, D a residue of c against a gap in x.
  std::vector<align::Alignment> to_centre(n);
  calx::scheduler::parallel_for(pool, n, 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t x = b; x < e; x++)
      if (x != c) to_centre[x] = align::gotoh(seqs[x], centre, s, align::Mode::global);
  });
  // The width of each gap of the centre.
  std::vector<std::size_t> width(L + 1, 0), here(L + 1);
  for (std::size_t x = 0; x < n; x++) {
    if (x == c) continue;
    std::fill(here.begin(), here.end(), 0);
    std::size_t p = 0;
    for (const auto& op : to_centre[x].cigar) {
      if (op.op == 'I') here[p] += op.length;
      else p += op.length;  // M or D: centre residues
    }
    for (std::size_t q = 0; q <= L; q++) width[q] = std::max(width[q], here[q]);
  }
  // The column of each centre residue, and where each gap of the centre begins.
  std::vector<std::size_t> gap_at(L + 1), residue_at(L);
  std::size_t col = 0;
  for (std::size_t q = 0; q <= L; q++) {
    gap_at[q] = col;
    col += width[q];
    if (q < L) residue_at[q] = col++;
  }
  out.msa.names = names;
  out.msa.rows.assign(n, Codes(col, GAP));
  for (std::size_t q = 0; q < L; q++) out.msa.rows[c][residue_at[q]] = centre[q];
  for (std::size_t x = 0; x < n; x++) {
    if (x == c) continue;
    Codes& row = out.msa.rows[x];
    std::size_t p = 0, i = 0, used = 0;  // the centre position, x's position, residues already in gap p
    for (const auto& op : to_centre[x].cigar)
      for (std::uint32_t k = 0; k < op.length; k++) {
        if (op.op == 'I') { row[gap_at[p] + used++] = seqs[x][i++]; continue; }
        if (op.op == 'M') row[residue_at[p]] = seqs[x][i++];
        p++;
        used = 0;
      }
  }
  return out;
}

}  // namespace genomics::msa
