// Extension of an alignment from a seed with the X-drop rule, ungapped and gapped, and minimap2's
// Z-drop.
//
// Ungapped (Altschul et al. 1990, J. Mol. Biol. 215, 403–410; Altschul et al. 1997, Nucleic Acids
// Res. 25, 3389–3402, §2): from a seed q[qpos, qpos + length) ~ t[tpos, tpos + length) the
// diagonal is extended to the right and to the left, a running sum kept in each direction and the
// extension stopped when the sum falls more than X below the best seen; each side keeps its best
// prefix. O(L) for an extension of length L.
//
// Gapped (Zhang, Schwartz, Wagner and Miller 2000, A greedy algorithm for aligning DNA sequences,
// J. Comput. Biol. 7, 203–214, §2, the X-drop dynamic programme; Altschul et al. 1997, §3, the
// gapped extension of BLAST 2): Gotoh's recurrences from the anchor vertex (qpos, tpos) to the
// right, row by row, with every cell whose score is more than X below the best seen so far (in
// row-major order) set to −∞, and each row computed only over the columns that the live cells of
// the row above can reach (one column beyond them on the diagonal, and as far right as a gap from
// a live cell stays alive). The extension stops when a row has no live cell, and returns the
// path to the best cell. To the left the same runs on the reversed prefixes. With X = ∞ the
// extension is the best alignment of q[qpos, ·) against t[tpos, ·) anchored at its start, with
// its end free; a finite X explores only a band that follows the alignment, at the cost of
// missing an extension that would recover after a drop of more than X.
//
// Z-drop (Li 2018, Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics 34,
// 3094–3100, §2.3.3): with (i′, j′) the best cell so far, of score S′, the extension also stops
// after a row whose best cell (i, j), of score S, has S′ − S > Z + e |(i − i′) − (j − j′)|, so that
// a long gap is not crossed where a run of mismatches would have stopped it, which separates an
// alignment from a neighbouring one across an inversion or a large indel. minimap2 tests this on
// each antidiagonal; here it is tested on each row, the order in which this schedule computes.
//
// Complexity: O(W) per row for a band of W live columns, O(L W) in all, with the directions of the
// visited cells kept for traceback. Deterministic: the tie-breaking of gotoh.hpp (pair, then I,
// then D; opening before extending), the best cell the first in row-major order.

#pragma once

#include <genomics/align/gotoh.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace genomics::align {

// The drop thresholds of an extension: X ≥ 0 always; Z < 0 disables the Z-drop.
struct Drop {
  std::int64_t x = 50;
  std::int64_t z = -1;
};

// ---- ungapped ---------------------------------------------------------------------------------------
// The seed of `length` pairs at (qpos, tpos) extended both ways along its diagonal.
inline Alignment ungapped_xdrop(const Codes& query, const Codes& target, std::size_t qpos, std::size_t tpos, std::size_t length, const Matrix& matrix, std::int64_t x) {
  if (qpos + length > query.size() || tpos + length > target.size()) throw std::invalid_argument("genomics::align::ungapped_xdrop: the seed leaves a sequence");
  if (x < 0) throw std::invalid_argument("genomics::align::ungapped_xdrop: X must be nonnegative");
  std::int64_t seed = 0;
  for (std::size_t k = 0; k < length; k++) seed += matrix(query[qpos + k], target[tpos + k]);
  // Right of the seed.
  std::int64_t run = 0, best_right = 0;
  std::size_t right = 0;
  for (std::size_t k = 0; qpos + length + k < query.size() && tpos + length + k < target.size(); k++) {
    run += matrix(query[qpos + length + k], target[tpos + length + k]);
    if (run > best_right) { best_right = run; right = k + 1; }
    if (run < best_right - x) break;
  }
  // Left of it.
  run = 0;
  std::int64_t best_left = 0;
  std::size_t left = 0;
  for (std::size_t k = 1; k <= qpos && k <= tpos; k++) {
    run += matrix(query[qpos - k], target[tpos - k]);
    if (run > best_left) { best_left = run; left = k; }
    if (run < best_left - x) break;
  }
  Alignment a;
  a.score = seed + best_left + best_right;
  a.query_begin = qpos - left;
  a.target_begin = tpos - left;
  a.query_end = qpos + length + right;
  a.target_end = tpos + length + right;
  push(a.cigar, 'M', static_cast<std::uint32_t>(left + length + right));
  return a;
}

// ---- gapped -------------------------------------------------------------------------------------------
namespace detail {

// The X-drop extension of q[qpos, ·) against t[tpos, ·) anchored at (qpos, tpos), in coordinates
// relative to the anchor; letters through qa(i) and tb(j) so that the left extension runs on the
// reversed prefixes without copying.
template <class QA, class TB>
Alignment extend(std::size_t m, std::size_t n, QA qa, TB tb, const Scoring& s, Drop drop) {
  const std::int64_t e = s.gap_extend, oe = static_cast<std::int64_t>(s.gap_open) + s.gap_extend;
  const int sigma = s.matrix.size;
  const int* score = s.matrix.score.data();
  std::vector<std::int64_t> H(n + 1, NEG), I(n + 1, NEG);
  Trace tr;
  std::int64_t best = 0;
  std::size_t bi = 0, bj = 0;
  const auto alive = [&](std::int64_t v) { return v >= best - drop.x && v > NEG / 2; };
  // Row 0: the anchor and the gaps along the target that stay alive.
  tr.lo.push_back(0);
  tr.offset.push_back(0);
  tr.dir.push_back(FROM_START);
  H[0] = 0;
  std::size_t lo = 0, hi = 0;
  {
    std::int64_t hl = 0, dl = NEG;
    for (std::size_t j = 1; j <= n; j++) {
      const std::int64_t dv_open = hl - oe, dv_ext = dl - e;
      const bool dext = dv_ext > dv_open;
      const std::int64_t dv = dext ? dv_ext : dv_open;
      if (!alive(dv)) break;
      H[j] = dv;
      hl = dl = dv;
      hi = j;
      tr.dir.push_back(static_cast<std::uint8_t>(FROM_D | (dext ? D_EXTENDS : 0)));
    }
    tr.hi.push_back(hi);
  }
  for (std::size_t i = 1; i <= m; i++) {
    const int* prof = score + static_cast<std::ptrdiff_t>(qa(i - 1)) * sigma;
    std::int64_t diag = lo > 0 ? H[lo - 1] : NEG;  // NEG: left of the window above
    std::int64_t hl = NEG, dl = NEG, row_best = NEG;
    std::size_t row_lo = std::numeric_limits<std::size_t>::max(), row_hi = 0, row_bj = 0;
    tr.lo.push_back(lo);
    tr.offset.push_back(tr.dir.size());
    std::size_t j = lo;
    for (; j <= n; j++) {
      const std::int64_t dv_open = hl - oe, dv_ext = dl - e;
      const bool dext = dv_ext > dv_open;
      const std::int64_t dv = dext ? dv_ext : dv_open;
      // Beyond the reach of the row above only a gap along the row can live.
      if (j > hi + 1 && !alive(dv)) break;
      const std::int64_t up_h = H[j], up_i = I[j];
      const std::int64_t iv_open = up_h - oe, iv_ext = up_i - e;
      const bool iext = iv_ext > iv_open;
      std::int64_t iv = iext ? iv_ext : iv_open;
      const std::int64_t mv = j > 0 ? diag + prof[tb(j - 1)] : NEG;
      diag = up_h;
      std::int64_t h = mv;
      std::uint8_t src = FROM_PAIR;
      if (iv > h) { h = iv; src = FROM_I; }
      if (dv > h) { h = dv; src = FROM_D; }
      std::int64_t dkeep = dv;
      if (!alive(h)) {
        h = iv = dkeep = NEG;
      } else {
        row_lo = std::min(row_lo, j);
        row_hi = j;
        if (h > row_best) { row_best = h; row_bj = j; }
        if (h > best) { best = h; bi = i; bj = j; }
      }
      H[j] = h;
      I[j] = iv;
      hl = h;
      dl = dkeep;
      tr.dir.push_back(static_cast<std::uint8_t>(src | (iext ? I_EXTENDS : 0) | (dext ? D_EXTENDS : 0)));
    }
    const std::size_t computed_hi = j - 1;
    tr.hi.push_back(computed_hi);
    // Cells of the window above that this row did not reach, and dead cells, are −∞ for the next.
    for (std::size_t c = computed_hi + 1; c <= hi && c <= n; c++) H[c] = I[c] = NEG;
    if (row_lo == std::numeric_limits<std::size_t>::max()) break;
    for (std::size_t c = lo; c < row_lo; c++) H[c] = I[c] = NEG;
    for (std::size_t c = row_hi + 1; c <= computed_hi; c++) H[c] = I[c] = NEG;
    lo = row_lo;
    hi = row_hi;
    if (drop.z >= 0) {
      const std::int64_t di = static_cast<std::int64_t>(i) - static_cast<std::int64_t>(bi), dj = static_cast<std::int64_t>(row_bj) - static_cast<std::int64_t>(bj);
      if (best - row_best > drop.z + e * (di > dj ? di - dj : dj - di)) break;
    }
  }
  // The path to the best cell.
  Alignment a;
  a.score = best;
  a.query_end = bi;
  a.target_end = bj;
  std::size_t i = bi, j = bj;
  enum { IN_H, IN_I, IN_D } state = IN_H;
  Cigar rev;
  for (;;) {
    const std::uint8_t d = tr.dir[tr.offset[i] + (j - tr.lo[i])];
    if (state == IN_H) {
      const std::uint8_t src = d & 3;
      if (src == FROM_START) break;
      if (src == FROM_PAIR) { push(rev, 'M'); i--; j--; continue; }
      state = src == FROM_I ? IN_I : IN_D;
      continue;
    }
    if (state == IN_I) {
      push(rev, 'I');
      state = (d & I_EXTENDS) ? IN_I : IN_H;
      i--;
    } else {
      push(rev, 'D');
      state = (d & D_EXTENDS) ? IN_D : IN_H;
      j--;
    }
  }
  a.cigar = reversed(std::move(rev));
  return a;
}

inline void check_extension(const Codes& query, const Codes& target, std::size_t qpos, std::size_t tpos, const Scoring& s, Drop drop, const char* who) {
  validate(query, target, s, who);
  if (qpos > query.size() || tpos > target.size()) throw std::invalid_argument(std::string(who) + ": the anchor lies outside the sequences");
  if (drop.x < 0) throw std::invalid_argument(std::string(who) + ": X must be nonnegative");
}

}  // namespace detail

// The extension to the right of the anchor vertex (qpos, tpos): an alignment from it to the best
// cell the X-drop (and Z-drop) rule reaches.
inline Alignment extend_right(const Codes& query, const Codes& target, std::size_t qpos, std::size_t tpos, const Scoring& s, Drop drop = {}) {
  detail::check_extension(query, target, qpos, tpos, s, drop, "genomics::align::extend_right");
  Alignment a = detail::extend(query.size() - qpos, target.size() - tpos, [&](std::size_t i) { return query[qpos + i]; }, [&](std::size_t j) { return target[tpos + j]; }, s, drop);
  a.query_begin = qpos;
  a.target_begin = tpos;
  a.query_end += qpos;
  a.target_end += tpos;
  return a;
}

// The extension to the left: an alignment ending at the anchor vertex (qpos, tpos).
inline Alignment extend_left(const Codes& query, const Codes& target, std::size_t qpos, std::size_t tpos, const Scoring& s, Drop drop = {}) {
  detail::check_extension(query, target, qpos, tpos, s, drop, "genomics::align::extend_left");
  Alignment r = detail::extend(qpos, tpos, [&](std::size_t i) { return query[qpos - 1 - i]; }, [&](std::size_t j) { return target[tpos - 1 - j]; }, s, drop);
  Alignment a;
  a.score = r.score;
  a.query_begin = qpos - r.query_end;
  a.target_begin = tpos - r.target_end;
  a.query_end = qpos;
  a.target_end = tpos;
  a.cigar = reversed(std::move(r.cigar));
  return a;
}

// The two extensions joined at the anchor. The score is the joined path's (rescored, which may
// exceed the sum of the two when a gap runs through the anchor and is opened only once).
inline Alignment xdrop(const Codes& query, const Codes& target, std::size_t qpos, std::size_t tpos, const Scoring& s, Drop drop = {}) {
  const Alignment l = extend_left(query, target, qpos, tpos, s, drop), r = extend_right(query, target, qpos, tpos, s, drop);
  Alignment a;
  a.query_begin = l.query_begin;
  a.target_begin = l.target_begin;
  a.query_end = r.query_end;
  a.target_end = r.target_end;
  a.cigar = l.cigar;
  for (const auto& o : r.cigar) push(a.cigar, o.op, o.length);
  a.score = rescore(a, query, target, s);
  return a;
}

}  // namespace genomics::align
