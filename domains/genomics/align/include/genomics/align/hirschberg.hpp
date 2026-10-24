// Optimal alignment with affine gaps in linear space: Myers and Miller's divide and conquer
// (Myers and Miller 1988, Optimal alignments in linear space, CABIOS 4, 11–17), which carries
// Hirschberg's (1975, CACM 18, 341–343) to Gotoh's three states; and the other modes by locating
// the ends first.
//
// Global. Split the query at its middle row i* = ⌊(a + b)/2⌋. A forward pass of Gotoh's
// recurrences over the rows a .. i* gives, for every column j, the best score F(j) of a path from
// the start to (i*, j) and FI(j) of one that ends with an I step; a backward pass over the
// reversed rows b .. i* gives B(j) and BI(j) for paths from (i*, j) to the end, and those that
// begin with an I step. An optimal path crosses row i* either at a vertex, the last of row i* it
// visits (type 1), or inside a vertical gap that runs through it (type 2), which F and B would
// each have opened once:
//
//   optimum = max_j max(F(j) + B(j), FI(j) + BI(j) + o).
//
// Type 1 recurses on (a .. i*, c .. j) and (i* .. b, j .. d). Type 2 takes out the two I steps
// into and out of (i*, j), which carry the gap's open, and recurses on (a .. i* − 1, c .. j) and
// (i* + 1 .. b, j .. d) with the open of a vertical gap at the shared boundary set to 0: each
// sub-problem has a top and a bottom open cost (o, or 0 when the gap there continues one already
// paid), Myers and Miller's tb and te. A sub-problem of one row, or of no row or no column, is
// solved directly (Myers and Miller 1988, §3). Any split of a horizontal gap at row i* costs an
// extra open, so the maximum is attained where no gap is split, and the halves' optima concatenate
// to an optimal path.
//
// Other modes. The end of an optimal path is found by gotoh_end (O(n) memory), its start by the
// same recurrences run backwards from that end with the start free as the mode allows, and the
// path between them is the global alignment of the two substrings, whose optimum equals the
// mode's: every path between the two vertices is a path of the mode, and the optimal one is
// among them.
//
// Complexity: O(mn) time (about twice Gotoh's: each level of the recursion re-reads half the
// area of the level above), O(m + n) memory, recursion depth ⌈log₂ m⌉. The path is optimal; it
// may differ from gotoh()'s among co-optimal paths.

#pragma once

#include <genomics/align/gotoh.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace genomics::align {

namespace detail {

class MyersMiller {
 public:
  MyersMiller(const Codes& q, const Codes& t, const Scoring& s)
      : q_(q), t_(t), s_(s), o_(s.gap_open), e_(s.gap_extend), F_(t.size() + 1), FI_(t.size() + 1), B_(t.size() + 1), BI_(t.size() + 1), H_(t.size() + 1), I_(t.size() + 1) {}

  // The optimal path of q[a, b) against t[c, d), appended to `out`.
  void solve(std::size_t a, std::size_t b, std::size_t c, std::size_t d, std::int64_t top_open, std::int64_t bottom_open, Cigar& out) {
    const std::size_t M = b - a, N = d - c;
    if (N == 0) {
      push(out, 'I', static_cast<std::uint32_t>(M));
      return;
    }
    if (M == 0) {
      push(out, 'D', static_cast<std::uint32_t>(N));
      return;
    }
    if (M == 1) {
      one_row(a, c, d, top_open, bottom_open, out);
      return;
    }
    const std::size_t mid = (a + b) / 2;
    // Forward over q[a, mid) and t[c, d); backward over q[mid, b) and t[c, d), reversed.
    pass(mid - a, N, [&](std::size_t k) { return q_[a + k]; }, [&](std::size_t l) { return t_[c + l]; }, top_open);
    for (std::size_t l = 0; l <= N; l++) { F_[l] = H_[l]; FI_[l] = I_[l]; }
    pass(b - mid, N, [&](std::size_t k) { return q_[b - 1 - k]; }, [&](std::size_t l) { return t_[d - 1 - l]; }, bottom_open);
    for (std::size_t l = 0; l <= N; l++) { B_[N - l] = H_[l]; BI_[N - l] = I_[l]; }
    std::int64_t best = NEG;
    std::size_t at = 0;
    bool gap = false;
    for (std::size_t l = 0; l <= N; l++) {
      const std::int64_t one = F_[l] + B_[l], two = FI_[l] + BI_[l] + o_;
      if (one > best) { best = one; at = l; gap = false; }
      if (two > best) { best = two; at = l; gap = true; }
    }
    const std::size_t j = c + at;
    if (!gap) {
      solve(a, mid, c, j, top_open, o_, out);
      solve(mid, b, j, d, o_, bottom_open, out);
    } else {
      solve(a, mid - 1, c, j, top_open, 0, out);
      push(out, 'I', 2);
      solve(mid + 1, b, j, d, 0, bottom_open, out);
    }
  }

 private:
  // Gotoh's recurrences over rows k = 0 .. rows and columns l = 0 .. cols of the letters qa(k),
  // tb(l), from the vertex (0, 0) only, a vertical gap there costing top_open instead of o; the
  // last row's H and I left in H_ and I_.
  template <class QA, class TB>
  void pass(std::size_t rows, std::size_t cols, QA qa, TB tb, std::int64_t top_open) {
    const std::int64_t oe = o_ + e_;
    const int sigma = s_.matrix.size;
    const int* score = s_.matrix.score.data();
    H_[0] = 0;
    I_[0] = -top_open;  // the virtual gap at the start, so that the first I step costs top_open + e
    for (std::size_t l = 1; l <= cols; l++) {
      H_[l] = -(o_ + static_cast<std::int64_t>(l) * e_);
      I_[l] = NEG;
    }
    for (std::size_t k = 1; k <= rows; k++) {
      const int* prof = score + static_cast<std::ptrdiff_t>(qa(k - 1)) * sigma;
      std::int64_t diag = H_[0];
      const std::int64_t i0 = std::max(H_[0] - oe, I_[0] - e_);
      H_[0] = i0;
      I_[0] = i0;
      std::int64_t hl = i0, dl = NEG;
      for (std::size_t l = 1; l <= cols; l++) {
        const std::int64_t iv = std::max(H_[l] - oe, I_[l] - e_);
        const std::int64_t dv = std::max(hl - oe, dl - e_);
        const std::int64_t mv = diag + prof[tb(l - 1)];
        diag = H_[l];
        const std::int64_t h = std::max(mv, std::max(iv, dv));
        H_[l] = h;
        I_[l] = iv;
        hl = h;
        dl = dv;
      }
    }
    if (rows == 0) I_[0] = NEG;  // the virtual gap is not a path
  }

  // One query letter q[a] against t[c, d), d > c.
  void one_row(std::size_t a, std::size_t c, std::size_t d, std::int64_t top_open, std::int64_t bottom_open, Cigar& out) {
    const auto gap = [&](std::size_t L) { return L ? -(o_ + static_cast<std::int64_t>(L) * e_) : std::int64_t(0); };
    const std::size_t N = d - c;
    std::int64_t best = NEG;
    int choice = 0;  // 0: q[a] aligned with t[k]; 1: I then D^N; 2: D^N then I
    std::size_t k_best = c;
    for (std::size_t k = c; k < d; k++) {
      const std::int64_t v = gap(k - c) + s_.matrix(q_[a], t_[k]) + gap(d - 1 - k);
      if (v > best) { best = v; k_best = k; choice = 0; }
    }
    const std::int64_t first = -(top_open + e_) + gap(N), last = gap(N) - (bottom_open + e_);
    if (first > best) { best = first; choice = 1; }
    if (last > best) { best = last; choice = 2; }
    if (choice == 0) {
      push(out, 'D', static_cast<std::uint32_t>(k_best - c));
      push(out, 'M');
      push(out, 'D', static_cast<std::uint32_t>(d - 1 - k_best));
    } else if (choice == 1) {
      push(out, 'I');
      push(out, 'D', static_cast<std::uint32_t>(N));
    } else {
      push(out, 'D', static_cast<std::uint32_t>(N));
      push(out, 'I');
    }
  }

  const Codes& q_;
  const Codes& t_;
  const Scoring& s_;
  std::int64_t o_, e_;
  std::vector<std::int64_t> F_, FI_, B_, BI_, H_, I_;
};

}  // namespace detail

// An optimal alignment in O(m + n) memory, with its path.
inline Alignment hirschberg(const Codes& query, const Codes& target, const Scoring& s, Mode m = Mode::global) {
  detail::validate(query, target, s, "genomics::align::hirschberg");
  std::size_t qb = 0, qe = query.size(), tb = 0, te = target.size();
  if (m != Mode::global) {
    const Endpoint end = gotoh_end(query, target, s, m);
    qe = end.query_end;
    te = end.target_end;
    // The start: the same recurrences on the reversed prefixes, from the end, free to stop where
    // the mode lets a path begin. Reversal carries the first row (and column) of the original to
    // the last row (and column) of the reversed problem, so the reversed sinks are the mode's own.
    const Codes rq(query.rbegin() + static_cast<std::ptrdiff_t>(query.size() - qe), query.rend());
    const Codes rt(target.rbegin() + static_cast<std::ptrdiff_t>(target.size() - te), target.rend());
    const auto start = detail::fill<false>(rq, rt, s, Mode::global, m, -static_cast<std::int64_t>(rq.size()), static_cast<std::int64_t>(rt.size()), nullptr);
    qb = qe - start.i;
    tb = te - start.j;
  }
  Alignment a;
  a.query_begin = qb;
  a.query_end = qe;
  a.target_begin = tb;
  a.target_end = te;
  detail::MyersMiller mm(query, target, s);
  mm.solve(qb, qe, tb, te, s.gap_open, s.gap_open, a.cigar);
  a.score = rescore(a, query, target, s);
  return a;
}

}  // namespace genomics::align
