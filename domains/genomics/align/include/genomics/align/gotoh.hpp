// Optimal alignment with affine gaps by Gotoh's recurrences, with full traceback, in the four
// modes; the optimal score and the end of an optimal path in linear memory; and the whole field
// of the recurrences, for drawing the edit graph (gotoh_field).
//
// The answer is defined by morphe's generic sum over the paths of the edit graph in the
// tropical semiring (morphe/formal/alignment.hpp, with three states); this header is one
// schedule of that sum, row by row, for one semiring and integer scores (theory note
// alignment.md §4). With the query as rows (i = 0 .. m) and the target as columns
// (j = 0 .. n), a gap of length L scoring −(o + L e) (types.hpp), Gotoh (1982, J. Mol. Biol.
// 162, 705–708, eqs. 3–5, in the corrected form of Altschul and Erickson 1986, Bull. Math.
// Biol. 48, 603–616) keeps for every vertex (i, j):
//
//   I(i, j) = max(H(i − 1, j) − o − e, I(i − 1, j) − e)        a path ending with a query base
//                                                              against a gap (CIGAR I)
//   D(i, j) = max(H(i, j − 1) − o − e, D(i, j − 1) − e)        ending with a target base
//                                                              against a gap (CIGAR D)
//   H(i, j) = max(start(i, j), H(i − 1, j − 1) + s(q_i, t_j), I(i, j), D(i, j))
//
// where start(i, j) = 0 at the vertices where the mode lets a path begin, and −∞ elsewhere:
// global (0, 0); local every vertex (Smith and Waterman 1981: the semiring's one at every
// cell); semiglobal the first row (the target's start free); overlap the first row and the
// first column. The optimum is the maximum of H over the vertices where the mode lets a path
// end: global (m, n); local every vertex; semiglobal the last row; overlap the last row and
// the last column. With o ≥ 0 a gap is never split, so I and D are exact (Altschul and
// Erickson 1986); a negative open is refused.
//
// Deterministic tie-breaking, so that a pair of sequences has one answer on every machine:
//   - the end: the first optimal vertex in row-major order (smallest query end, then smallest
//     target end);
//   - at H: start, then the pair (M), then I, then D, the first that attains the maximum: the
//     traceback stops as soon as it may, and gaps are placed as far left as they can go
//     (left-aligned indels, the normal form of VCF);
//   - in a gap: opening (from H) before extending, so that runs are as short as the score
//     allows.
// A local alignment of score 0 is the empty alignment at (0, 0).
//
// Complexity: O(mn) time; the traceback keeps one byte of directions per vertex
// ((m + 1)(n + 1) bytes, 100 MB for 10 kb × 10 kb), the scores two rows of O(n). gotoh_score
// and gotoh_end use O(n) memory. Scores are 64-bit; there is no overflow below 10^15 cells.
// Codes must be below the matrix's size (sequence::GAP is not a code).

#pragma once

#include <genomics/align/types.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace genomics::align {

// The score of an optimal path and the vertex where it ends (the query and target ends of the
// alignment).
struct Endpoint {
  std::int64_t score = 0;
  std::size_t query_end = 0, target_end = 0;
};

// The whole field of Gotoh's recurrences, for drawing the edit graph: (|query| + 1) × (|target| + 1)
// values each, row-major, `unreachable` where no path of the mode reaches the vertex in that
// state.
struct ScoreField {
  static constexpr std::int64_t unreachable = std::numeric_limits<std::int64_t>::min();
  std::size_t rows = 0, cols = 0;
  std::vector<std::int64_t> best;      // H(i, j): the best path from where the mode starts to (i, j)
  std::vector<std::int64_t> ending_i;  // I(i, j): the best such path ending with a query base against a gap
  std::vector<std::int64_t> ending_d;  // D(i, j): the best ending with a target base against a gap
  std::int64_t at(std::size_t i, std::size_t j) const { return best[i * cols + j]; }
};

namespace detail {

inline constexpr std::int64_t NEG = std::numeric_limits<std::int64_t>::min() / 4;

// Directions of one vertex: bits 0–1 the source of H (0 start, 1 pair, 2 I, 3 D); bit 2 set when
// I extends I(i − 1, j); bit 3 when D extends D(i, j − 1).
inline constexpr std::uint8_t FROM_START = 0, FROM_PAIR = 1, FROM_I = 2, FROM_D = 3, I_EXTENDS = 4, D_EXTENDS = 8;

// The directions of the vertices of a band, row by row: row i holds the columns
// [lo[i], hi[i]] at offset[i].
struct Trace {
  std::vector<std::uint8_t> dir;
  std::vector<std::size_t> lo, hi, offset;
  std::uint8_t at(std::size_t i, std::size_t j) const { return dir[offset[i] + (j - lo[i])]; }
};

struct Best {
  std::int64_t score = NEG;
  std::size_t i = 0, j = 0;
  bool found = false;
};

inline void validate(const Codes& query, const Codes& target, const Scoring& s, const char* who) {
  if (s.matrix.size <= 0 || s.matrix.score.size() != static_cast<std::size_t>(s.matrix.size) * static_cast<std::size_t>(s.matrix.size))
    throw std::invalid_argument(std::string(who) + ": the substitution matrix is not square");
  if (s.gap_open < 0) throw std::invalid_argument(std::string(who) + ": a negative gap open would split gaps");
  for (Code c : query)
    if (c >= s.matrix.size) throw std::invalid_argument(std::string(who) + ": a query code outside the matrix");
  for (Code c : target)
    if (c >= s.matrix.size) throw std::invalid_argument(std::string(who) + ": a target code outside the matrix");
}

// The columns of row i inside the diagonals [low, high] (j − i), clipped to [0, n]; lo > hi
// when the row has none.
inline void row_range(std::size_t i, std::size_t n, std::int64_t low, std::int64_t high, std::size_t& lo, std::size_t& hi) {
  const std::int64_t a = static_cast<std::int64_t>(i) + low, b = static_cast<std::int64_t>(i) + high;
  const std::int64_t nn = static_cast<std::int64_t>(n);
  if (b < 0 || a > nn) { lo = 1; hi = 0; return; }
  lo = static_cast<std::size_t>(a < 0 ? 0 : a);
  hi = static_cast<std::size_t>(b > nn ? nn : b);
}

// Gotoh's recurrences over the vertices of a band, row by row, keeping two rows of scores;
// the directions too when `trace` is given. Paths begin where the mode `from` lets them and end
// where the mode `to` does (the two differ only in the reverse pass of hirschberg.hpp, which
// fixes the start). Returns the optimal end by the rule above. With Field, every vertex's H, I
// and D are kept in `field` too (gotoh_field).
template <bool Keep, bool Field = false>
Best fill(const Codes& query, const Codes& target, const Scoring& s, Mode from, Mode to, std::int64_t low, std::int64_t high, Trace* trace, ScoreField* field = nullptr) {
  const std::size_t m = query.size(), n = target.size();
  const std::int64_t e = s.gap_extend, oe = static_cast<std::int64_t>(s.gap_open) + s.gap_extend;
  const int sigma = s.matrix.size;
  const int* score = s.matrix.score.data();
  std::vector<std::int64_t> H(n + 1, NEG), I(n + 1, NEG);
  Best best;
  const bool all_sinks = to == Mode::local;
  if constexpr (Keep) {
    trace->lo.assign(m + 1, 1);
    trace->hi.assign(m + 1, 0);
    trace->offset.assign(m + 1, 0);
    std::size_t cells = 0;
    for (std::size_t i = 0; i <= m; i++) {
      row_range(i, n, low, high, trace->lo[i], trace->hi[i]);
      trace->offset[i] = cells;
      if (trace->lo[i] <= trace->hi[i]) cells += trace->hi[i] - trace->lo[i] + 1;
    }
    trace->dir.assign(cells, 0);
  }
  if constexpr (Field) {
    field->rows = m + 1;
    field->cols = n + 1;
    field->best.assign((m + 1) * (n + 1), ScoreField::unreachable);
    field->ending_i.assign((m + 1) * (n + 1), ScoreField::unreachable);
    field->ending_d.assign((m + 1) * (n + 1), ScoreField::unreachable);
  }
  const auto record = [&](std::size_t i, std::size_t j, std::int64_t h, std::int64_t iv, std::int64_t dv) {
    if constexpr (Field) {
      const auto keep = [](std::int64_t v) { return v > NEG / 2 ? v : ScoreField::unreachable; };
      const std::size_t c = i * (n + 1) + j;
      field->best[c] = keep(h);
      field->ending_i[c] = keep(iv);
      field->ending_d[c] = keep(dv);
    } else {
      (void)i, (void)j, (void)h, (void)iv, (void)dv;
    }
  };
  // A vertex no path reaches holds −∞ plus a little, never more than NEG / 2.
  const auto consider = [&](std::int64_t h, std::size_t i, std::size_t j) {
    if (h > best.score && h > NEG / 2) { best.score = h; best.i = i; best.j = j; best.found = true; }
  };
  // The sinks other than local's (which are every vertex, met in the loops), in row-major
  // order: (i, n) of the overlap's last column, the last row for semiglobal and overlap,
  // (m, n) for global.
  const auto row_sinks = [&](std::size_t i, std::size_t lo, std::size_t hi) {
    if (all_sinks) return;
    if (i < m) {
      if (to == Mode::overlap && hi == n) consider(H[n], i, n);
    } else if (to == Mode::global) {
      if (hi == n) consider(H[n], m, n);
    } else {
      for (std::size_t c = lo; c <= hi; c++) consider(H[c], m, c);
    }
  };
  for (std::size_t i = 0; i <= m; i++) {
    std::size_t lo, hi;
    row_range(i, n, low, high, lo, hi);
    if (lo > hi) continue;
    std::uint8_t* dir = nullptr;
    if constexpr (Keep) dir = trace->dir.data() + trace->offset[i];
    std::int64_t hl = NEG, dl = NEG;  // H and D of the left neighbour, (i, j − 1)
    std::size_t j = lo;
    if (i == 0) {
      // The first row: paths of D steps only, and starts where the mode allows them.
      const std::int64_t s0 = (from == Mode::global) ? NEG : 0;
      for (; j <= hi; j++) {
        const std::int64_t start = (j == 0) ? 0 : s0;
        const std::int64_t dv_open = hl - oe, dv_ext = dl - e;
        const bool dext = dv_ext > dv_open;
        const std::int64_t dv = j > lo ? (dext ? dv_ext : dv_open) : NEG;
        std::int64_t h = start;
        std::uint8_t src = FROM_START;
        if (dv > h) { h = dv; src = FROM_D; }
        H[j] = h;
        I[j] = NEG;
        hl = h;
        dl = dv;
        if constexpr (Keep) dir[j - lo] = static_cast<std::uint8_t>(src | (dext && j > lo ? D_EXTENDS : 0));
        record(i, j, h, NEG, dv);
        if (all_sinks) consider(h, i, j);
      }
      row_sinks(i, lo, hi);
      continue;
    }
    const std::int64_t s0 = from == Mode::local ? 0 : NEG;  // start at an inner vertex
    const int* prof = score + static_cast<std::ptrdiff_t>(query[i - 1]) * sigma;
    std::int64_t diag = lo > 0 ? H[lo - 1] : NEG;  // H(i − 1, j − 1)
    if (j == 0) {
      // The first column: I steps only, and a start in local and overlap modes.
      const std::int64_t start = (from == Mode::local || from == Mode::overlap) ? 0 : NEG;
      const std::int64_t iv_open = H[0] - oe, iv_ext = I[0] - e;
      const bool iext = iv_ext > iv_open;
      const std::int64_t iv = iext ? iv_ext : iv_open;
      std::int64_t h = start;
      std::uint8_t src = FROM_START;
      if (iv > h) { h = iv; src = FROM_I; }
      diag = H[0];
      H[0] = h;
      I[0] = iv;
      hl = h;
      dl = NEG;
      if constexpr (Keep) dir[0] = static_cast<std::uint8_t>(src | (iext ? I_EXTENDS : 0));
      record(i, 0, h, iv, NEG);
      if (all_sinks) consider(h, i, 0);
      j = 1;
    }
    for (; j <= hi; j++) {
      const std::int64_t up_h = H[j], up_i = I[j];
      const std::int64_t iv_open = up_h - oe, iv_ext = up_i - e;
      const bool iext = iv_ext > iv_open;
      const std::int64_t iv = iext ? iv_ext : iv_open;
      const std::int64_t dv_open = hl - oe, dv_ext = dl - e;
      const bool dext = dv_ext > dv_open;
      const std::int64_t dv = dext ? dv_ext : dv_open;
      const std::int64_t mv = diag + prof[target[j - 1]];
      diag = up_h;
      // The first term attaining the maximum, in the order start, pair, I, D.
      std::int64_t h = s0;
      std::uint8_t src = FROM_START;
      if (mv > h) { h = mv; src = FROM_PAIR; }
      if (iv > h) { h = iv; src = FROM_I; }
      if (dv > h) { h = dv; src = FROM_D; }
      H[j] = h;
      I[j] = iv;
      hl = h;
      dl = dv;
      if constexpr (Keep) dir[j - lo] = static_cast<std::uint8_t>(src | (iext ? I_EXTENDS : 0) | (dext ? D_EXTENDS : 0));
      record(i, j, h, iv, dv);
      if (all_sinks) consider(h, i, j);
    }
    row_sinks(i, lo, hi);
  }
  return best;
}

// The path into the end found by fill, read back through the directions.
inline Alignment traceback(const Trace& trace, const Best& best) {
  Alignment a;
  a.score = best.score;
  a.query_end = best.i;
  a.target_end = best.j;
  std::size_t i = best.i, j = best.j;
  enum { IN_H, IN_I, IN_D } state = IN_H;
  Cigar rev;
  for (;;) {
    const std::uint8_t d = trace.at(i, j);
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
  a.query_begin = i;
  a.target_begin = j;
  a.cigar = reversed(std::move(rev));
  return a;
}

}  // namespace detail

// An optimal alignment of the query to the target in the mode, with its path.
inline Alignment gotoh(const Codes& query, const Codes& target, const Scoring& s, Mode m) {
  detail::validate(query, target, s, "genomics::align::gotoh");
  detail::Trace trace;
  const auto best = detail::fill<true>(query, target, s, m, m, -static_cast<std::int64_t>(query.size()), static_cast<std::int64_t>(target.size()), &trace);
  return detail::traceback(trace, best);
}

// The optimal score and the end of the optimal path gotoh() returns, in O(n) memory.
inline Endpoint gotoh_end(const Codes& query, const Codes& target, const Scoring& s, Mode m) {
  detail::validate(query, target, s, "genomics::align::gotoh_end");
  const auto best = detail::fill<false>(query, target, s, m, m, -static_cast<std::int64_t>(query.size()), static_cast<std::int64_t>(target.size()), nullptr);
  return {best.score, best.i, best.j};
}

// The whole field H, I, D of the recurrences in the mode, for drawing (O(mn) memory: meant for
// sequences of a few hundred letters). Its value at the end of gotoh()'s alignment is the score.
inline ScoreField gotoh_field(const Codes& query, const Codes& target, const Scoring& s, Mode m) {
  detail::validate(query, target, s, "genomics::align::gotoh_field");
  ScoreField field;
  detail::fill<false, true>(query, target, s, m, m, -static_cast<std::int64_t>(query.size()), static_cast<std::int64_t>(target.size()), nullptr, &field);
  return field;
}

// The optimal score, in O(n) memory.
inline std::int64_t gotoh_score(const Codes& query, const Codes& target, const Scoring& s, Mode m) { return gotoh_end(query, target, s, m).score; }

// The end of an optimal local alignment and its score, in O(n) memory: the columns of
// Smith and Waterman's matrix kept one row at a time.
inline Endpoint local_end(const Codes& query, const Codes& target, const Scoring& s) { return gotoh_end(query, target, s, Mode::local); }

}  // namespace genomics::align
