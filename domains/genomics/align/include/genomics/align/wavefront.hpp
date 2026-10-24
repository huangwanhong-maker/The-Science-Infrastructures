// The wavefront alignment algorithm (WFA) with gap-affine penalties, its bidirectional form in
// O(s) memory (BiWFA), the adaptive heuristic, and the exchange of penalties and scores.
//
// Penalties (types.hpp) are minimised: 0 for a match, x for a mismatch, o + L e for a gap of
// length L. A match is two equal codes below 4; N matches nothing (as align::penalty counts).
// With the query (length n) as rows v and the target (length m) as columns h, a diagonal is
// k = h − v and a point on it is named by its offset h. The wavefronts of a penalty s are the
// furthest offsets reached on each diagonal by a path of penalty exactly s, ending with a pair
// (M), a query base against a gap (I: v + 1, k − 1) or a target base against a gap (D: h + 1,
// k + 1). Marco-Sola, Moure, Moreto and Espinosa (2021, Bioinformatics 37, 456–463):
//
//   Ĩ(s, k) = max(M̃(s − o − e, k + 1), Ĩ(s − e, k + 1))
//   D̃(s, k) = max(M̃(s − o − e, k − 1), D̃(s − e, k − 1)) + 1
//   M̃(s, k) = extend(max(M̃(s − x, k) + 1, Ĩ(s, k), D̃(s, k)))
//
// where extend follows the diagonal while the bases match (here 8 bytes at a time: the xor of
// two words and its trailing zeros), and every candidate is dropped when it leaves the matrix.
// The first s at which M̃(s, m − n) = m is the optimal penalty of the global alignment; the
// furthest point dominates every nearer one on its diagonal because the cost to the end never
// increases along a diagonal (theory note fast-alignment.md §5). This is the (min, +) sum over
// the edit graph computed in the order of its values: the cells of penalty s are visited only
// after all cells of smaller penalty, and only the frontier is kept.
//
// - Modes: global; semiglobal (ends-free on the target: the query end to end, the target's
//   ends free, WFA2-lib's text_begin_free = text_end_free = m), whose wavefront 0 holds every
//   diagonal 0 .. m and which stops at the first s with a point on the last row (the smallest
//   target end among them).
// - Traceback: every wavefront is kept, and the path is read back from the end, preferring at
//   M the start, then a mismatch, then I, then D, and in a gap its opening before its extension
//   (the order of gotoh.hpp). Memory O(s²/e) offsets in the worst case.
// - BiWFA (Marco-Sola, Eizenga, Guarracino, Paten, Garrison and Moreto 2023, Bioinformatics 39,
//   btad074): a forward and a reverse wavefront, on the sequences and on their reverses, are
//   advanced in turn, keeping only the last max(x, o + e) + 1 of each, until they overlap on a
//   diagonal (offsets summing to m) in M, or in I or D with the gap's opening counted twice;
//   the least such sum is the optimal penalty, and the point is a breakpoint of an optimal path.
//   The two halves are aligned recursively: at a gap breakpoint the first half ends with a step
//   of the gap, which is peeled off, so that what precedes it ends inside the open gap (its last
//   run refunded its opening) and the second half begins inside it (its first run paying none).
//   Every sub-problem then ends freely or in an open gap, the ends toward which the furthest
//   point of a diagonal dominates the nearer ones (theory note §5.4); a constraint that forbade
//   an ending (a last step that must be a gap) would break that on the last row. Halves of
//   small penalty (≤ 256) by the WFA above. Memory O(s), time O((n + m) s) at worst.
// - WF-adaptive (Marco-Sola et al. 2021, its heuristic; WFA2-lib's wfadaptive with its defaults 10, 50,
//   1): off by default. After each wavefront of at least 10 diagonals, the diagonals whose
//   remaining distance max(n − v, m − h) exceeds the least by more than 50 are dropped. A
//   heuristic: the alignment is valid and its penalty an upper bound, not always the optimum.
// - Penalties and scores (Eizenga and Paten 2022, bioRxiv 2022.01.12.476087; also
//   Marco-Sola et al. 2021, Suppl.): a global alignment of lengths n and m with a matches, b
//   mismatches, g gaps and L gapped bases has 2a + 2b + L = n + m, so a score with match α,
//   mismatch −β and gaps −(q + L r) is S = (α(n + m) − P)/2 with the penalties x = 2(α + β),
//   o = 2q, e = 2r + α, P = 2(α + β) b + 2q g + (2r + α) L; the penalties divided by their
//   gcd. The two orders of global alignments are the same, so the WFA computes Gotoh's global
//   optimum (penalty_form, wavefront_global).
//
// Limits: sequences shorter than 2³⁰; penalties x, e > 0, o ≥ 0. Every function also takes an
// optional WavefrontStats (memory, penalties, offsets), and wavefront() a WavefrontRecord.

#pragma once

#include <genomics/align/types.hpp>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <deque>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace genomics::align {

struct WavefrontOptions {
  bool adaptive = false;            // WF-adaptive: a heuristic, off by default
  int min_wavefront_length = 10;    // wavefronts narrower than this are never cut
  int max_distance_threshold = 50;  // diagonals further behind the best are dropped
  int steps_between_cutoffs = 1;
  std::int64_t max_penalty = -1;    // give up beyond this penalty (−1: never)
};

// What a run cost: the most bytes of wavefronts held at once, the penalties visited, the
// offsets computed.
struct WavefrontStats {
  std::size_t peak_bytes = 0, current_bytes = 0;
  std::int64_t scores = 0;
  std::size_t cells = 0;
  void add(std::ptrdiff_t bytes) {
    current_bytes = static_cast<std::size_t>(static_cast<std::ptrdiff_t>(current_bytes) + bytes);
    peak_bytes = std::max(peak_bytes, current_bytes);
  }
};

// The wavefronts of a run, for drawing the schedule on the edit graph: for each penalty with a
// wavefront, its M̃ offsets on the diagonals lo .. hi (−1 where none), in the order computed;
// and the work done, against the (n + 1)(m + 1) cells of the full matrix.
struct WavefrontRecord {
  struct Front {
    std::int64_t score = 0, lo = 0, hi = 0;
    std::vector<std::int64_t> offsets;
  };
  std::vector<Front> fronts;
  std::size_t cells_extended = 0;    // matching bases walked by the extension
  std::size_t offsets_computed = 0;  // diagonals of every wavefront computed
};

// A scoring as penalties: score = (match (n + m) − scale · penalty) / 2 for every global
// alignment of a query of length n and a target of length m (without N, unless the matrix
// scores N as a mismatch).
struct PenaltyForm {
  Penalties penalties;
  std::int64_t match = 0, scale = 1;
};

namespace detail {

using WfOffset = std::int32_t;
inline constexpr WfOffset WF_NULL = std::numeric_limits<std::int32_t>::min() / 2;
enum WfComp : int { WF_M = 0, WF_I = 1, WF_D = 2 };
// The state at an end of a sub-problem: plain (M), or a gap of one kind open across it. At the
// beginning an open gap lets the path's first run of that kind pay no opening (it may also not
// continue it); at the end it refunds the opening of the path's last run of that kind (the
// mirror image). An empty path inside an open gap is in that gap at both ends.
enum class WfBegin { M, I_open, D_open };

struct WfSlot {
  std::int64_t score = -1;
  bool present = false;
  std::int32_t lo = 1, hi = 0, base = 0, width = 0;
  std::vector<WfOffset> buf;  // component c at buf[c · width + (k − base)]
  WfOffset at(int c, std::int32_t k) const { return present && k >= lo && k <= hi ? buf[static_cast<std::size_t>(c * width + (k - base))] : WF_NULL; }
  WfOffset* comp(int c) { return buf.data() + static_cast<std::size_t>(c * width); }
  const WfOffset* comp(int c) const { return buf.data() + static_cast<std::size_t>(c * width); }
};

// The wavefronts by penalty: all of them, or the last `window` in a ring.
class WfStore {
 public:
  WfStore(std::size_t window, WavefrontStats* stats) : window_(window), stats_(stats) {
    if (window_) ring_.resize(window_);
  }
  ~WfStore() {
    if (stats_) stats_->add(-static_cast<std::ptrdiff_t>(bytes_));
  }
  WfStore(const WfStore&) = delete;
  WfStore& operator=(const WfStore&) = delete;

  WfSlot& make(std::int64_t s, std::int32_t lo, std::int32_t hi) {
    WfSlot& w = slot_for(s);
    const std::size_t before = w.buf.capacity();
    w.score = s;
    w.present = true;
    w.lo = w.base = lo;
    w.hi = hi;
    w.width = hi - lo + 1;
    w.buf.assign(static_cast<std::size_t>(3 * w.width), WF_NULL);
    account(before, w.buf.capacity());
    return w;
  }
  void make_null(std::int64_t s) {
    WfSlot& w = slot_for(s);
    w.score = s;
    w.present = false;
  }
  const WfSlot* get(std::int64_t s) const {
    if (s < 0) return nullptr;
    if (window_) {
      const WfSlot& w = ring_[static_cast<std::size_t>(s) % window_];
      return w.score == s && w.present ? &w : nullptr;
    }
    if (static_cast<std::size_t>(s) >= all_.size()) return nullptr;
    const WfSlot& w = all_[static_cast<std::size_t>(s)];
    return w.present ? &w : nullptr;
  }
  WfSlot* get_mutable(std::int64_t s) { return const_cast<WfSlot*>(get(s)); }

 private:
  WfSlot& slot_for(std::int64_t s) {
    if (window_) return ring_[static_cast<std::size_t>(s) % window_];
    while (all_.size() <= static_cast<std::size_t>(s)) all_.emplace_back();
    return all_[static_cast<std::size_t>(s)];
  }
  void account(std::size_t before, std::size_t after) {
    const std::ptrdiff_t d = static_cast<std::ptrdiff_t>((after - before) * sizeof(WfOffset));
    bytes_ += static_cast<std::size_t>(d);
    if (stats_) stats_->add(d);
  }
  std::size_t window_;
  WavefrontStats* stats_;
  std::size_t bytes_ = 0;
  std::deque<WfSlot> all_;
  std::vector<WfSlot> ring_;
};

// The sequences as the extension reads them: codes ≥ 4 replaced by a byte no other sequence
// holds, and 8 bytes of a sentinel after the end, so that a word compared past the end stops.
inline std::vector<std::uint8_t> wf_padded(const Code* x, std::size_t n, bool reverse, std::uint8_t other, std::uint8_t sentinel) {
  std::vector<std::uint8_t> p(n + 16, sentinel);
  for (std::size_t i = 0; i < n; i++) {
    const Code c = reverse ? x[n - 1 - i] : x[i];
    p[i] = c < 4 ? c : other;
  }
  return p;
}

// The wavefronts of one direction, advanced one penalty at a time.
class WfEngine {
 public:
  WfEngine(const std::vector<std::uint8_t>& q, std::int32_t n, const std::vector<std::uint8_t>& t, std::int32_t m, const Penalties& p, WfBegin begin,
           bool ends_free, std::size_t window, const WavefrontOptions& options, WavefrontStats* stats, WavefrontRecord* record = nullptr)
      : q_(q.data()), t_(t.data()), n_(n), m_(m), x_(p.mismatch), oe_(p.gap_open + p.gap_extend), e_(p.gap_extend), begin_(begin), ends_free_(ends_free),
        options_(options), stats_(stats), record_(record), store_(window, stats) {}

  std::int64_t score() const { return s_; }
  const WfSlot* at(std::int64_t s) const { return store_.get(s); }
  std::int32_t n() const { return n_; }
  std::int32_t m() const { return m_; }

  // Computes the wavefront of the next penalty (0 first); null when it has none.
  const WfSlot* next() {
    s_++;
    const std::int64_t s = s_;
    if (stats_) stats_->scores++;
    const WfSlot* A = store_.get(s - x_);
    const WfSlot* B = store_.get(s - oe_);
    const WfSlot* C = store_.get(s - e_);
    // The diagonals the sources reach, and the seeds of the beginning.
    std::int32_t lo = std::numeric_limits<std::int32_t>::max(), hi = std::numeric_limits<std::int32_t>::min();
    const auto cover = [&](std::int32_t a, std::int32_t b) { lo = std::min(lo, a); hi = std::max(hi, b); };
    if (A) cover(A->lo, A->hi);
    if (B) cover(B->lo - 1, B->hi + 1);
    if (C) cover(C->lo - 1, C->hi + 1);
    const bool seed0 = s == 0;
    if (seed0) cover(0, ends_free_ ? m_ : 0);
    lo = std::max(lo, -n_);
    hi = std::min(hi, m_);
    if (lo > hi) {
      store_.make_null(s);
      return nullptr;
    }
    // make() moves no source: the deque keeps its elements, and a ring slot is another penalty's.
    WfSlot& w = store_.make(s, lo, hi);
    WfOffset* M = w.comp(WF_M);
    WfOffset* I = w.comp(WF_I);
    WfOffset* D = w.comp(WF_D);
    const std::int32_t n = n_, m = m_;
    if (B) {
      const WfOffset* bm = B->comp(WF_M);
      const std::int32_t bb = B->base;
      for (std::int32_t k = std::max(lo, B->lo - 1), end = std::min(hi, B->hi - 1); k <= end; k++) {
        const WfOffset h = bm[k + 1 - bb];
        I[k - lo] = (h >= 0 && h - k <= n) ? h : WF_NULL;
      }
      for (std::int32_t k = std::max(lo, B->lo + 1), end = std::min(hi, B->hi + 1); k <= end; k++) {
        const WfOffset h = bm[k - 1 - bb] + 1;
        D[k - lo] = (h > 0 && h <= m) ? h : WF_NULL;
      }
    }
    if (C) {
      const WfOffset* ci = C->comp(WF_I);
      const WfOffset* cd = C->comp(WF_D);
      const std::int32_t cb = C->base;
      for (std::int32_t k = std::max(lo, C->lo - 1), end = std::min(hi, C->hi - 1); k <= end; k++) {
        const WfOffset h = ci[k + 1 - cb];
        if (h >= 0 && h - k <= n) I[k - lo] = std::max(I[k - lo], h);
      }
      for (std::int32_t k = std::max(lo, C->lo + 1), end = std::min(hi, C->hi + 1); k <= end; k++) {
        const WfOffset h = cd[k - 1 - cb] + 1;
        if (h > 0 && h <= m) D[k - lo] = std::max(D[k - lo], h);
      }
    }
    if (A) {
      const WfOffset* am = A->comp(WF_M);
      const std::int32_t ab = A->base;
      for (std::int32_t k = std::max(lo, A->lo), end = std::min(hi, A->hi); k <= end; k++) {
        const WfOffset h = am[k - ab] + 1;
        if (h > 0 && h <= m && h - k <= n) M[k - lo] = h;
      }
    }
    const std::int32_t width = hi - lo + 1;
    for (std::int32_t i = 0; i < width; i++) M[i] = std::max(M[i], std::max(I[i], D[i]));
    if (seed0) {
      if (ends_free_)
        for (std::int32_t k = 0; k <= m; k++) M[k - lo] = std::max(M[k - lo], k);
      else M[-lo] = std::max(M[-lo], 0);
      if (begin_ == WfBegin::I_open) I[-lo] = std::max(I[-lo], 0);
      if (begin_ == WfBegin::D_open) D[-lo] = std::max(D[-lo], 0);
    }
    extend(M, lo, hi);
    if (stats_) stats_->cells += static_cast<std::size_t>(width);
    if (options_.adaptive) cut(w);
    if (record_) {
      record_->offsets_computed += static_cast<std::size_t>(width);
      WavefrontRecord::Front f;
      f.score = s;
      f.lo = w.lo;
      f.hi = w.hi;
      f.offsets.reserve(static_cast<std::size_t>(w.hi - w.lo + 1));
      for (std::int32_t k = w.lo; k <= w.hi; k++) {
        const WfOffset h = M[k - lo];
        f.offsets.push_back(h >= 0 ? h : -1);
      }
      record_->fronts.push_back(std::move(f));
    }
    return &w;
  }

 private:
  // Follows each diagonal while the bases match, a word of 8 at a time: the first differing
  // byte is the trailing zeros of the xor over 8. The sentinels stop it at the ends.
  void extend(WfOffset* M, std::int32_t lo, std::int32_t hi) const {
    for (std::int32_t k = lo; k <= hi; k++) {
      WfOffset h = M[k - lo];
      if (h < 0) continue;
      const WfOffset h0 = h;
      std::int64_t v = static_cast<std::int64_t>(h) - k;
      for (;;) {
        std::uint64_t a, b;
        std::memcpy(&a, q_ + v, 8);
        std::memcpy(&b, t_ + h, 8);
        const std::uint64_t d = a ^ b;
        if (d == 0) { h += 8; v += 8; continue; }
        h += static_cast<WfOffset>(std::countr_zero(d) >> 3);
        break;
      }
      M[k - lo] = h;
      if (record_) record_->cells_extended += static_cast<std::size_t>(h - h0);
    }
  }
  // WF-adaptive: drop the outer diagonals that lag the leader by more than the threshold.
  void cut(WfSlot& w) const {
    if (w.hi - w.lo + 1 < options_.min_wavefront_length || (options_.steps_between_cutoffs > 1 && s_ % options_.steps_between_cutoffs != 0)) return;
    const WfOffset* M = w.comp(WF_M);
    const auto distance = [&](std::int32_t k) -> std::int64_t {
      const WfOffset h = M[k - w.base];
      if (h < 0) return std::numeric_limits<std::int64_t>::max() / 4;
      const std::int64_t left_v = n_ - (static_cast<std::int64_t>(h) - k), left_h = m_ - static_cast<std::int64_t>(h);
      return ends_free_ ? left_v : std::max(left_v, left_h);
    };
    std::int64_t least = std::numeric_limits<std::int64_t>::max();
    for (std::int32_t k = w.lo; k <= w.hi; k++) least = std::min(least, distance(k));
    while (w.lo < w.hi && distance(w.lo) - least > options_.max_distance_threshold) w.lo++;
    while (w.hi > w.lo && distance(w.hi) - least > options_.max_distance_threshold) w.hi--;
  }

  const std::uint8_t* q_;
  const std::uint8_t* t_;
  std::int32_t n_, m_;
  std::int64_t x_, oe_, e_;
  WfBegin begin_;
  bool ends_free_;
  WavefrontOptions options_;
  WavefrontStats* stats_;
  WavefrontRecord* record_;
  WfStore store_;
  std::int64_t s_ = -1;
};

inline void wf_validate(const Codes& q, const Codes& t, const Penalties& p, const char* who) {
  if (p.mismatch <= 0 || p.gap_extend <= 0 || p.gap_open < 0) throw std::invalid_argument(std::string(who) + ": the penalties need x > 0, e > 0 and o >= 0");
  if (q.size() >= (std::size_t{1} << 30) || t.size() >= (std::size_t{1} << 30)) throw std::invalid_argument(std::string(who) + ": sequences of 2^30 or more");
}

// Where a run ended: whether it did, its penalty (less a refunded opening, so possibly
// negative), the penalty of the wavefront that holds the end, and the diagonal, offset and
// component of the end.
struct WfEnd {
  bool found = false;
  std::int64_t score = 0, raw = 0;
  std::int32_t k = 0;
  WfOffset h = 0;
  int comp = WF_M;
};

inline int wf_comp_of(WfBegin b) { return b == WfBegin::I_open ? WF_I : b == WfBegin::D_open ? WF_D : WF_M; }
inline WfBegin wf_open(int c) { return c == WF_I ? WfBegin::I_open : c == WF_D ? WfBegin::D_open : WfBegin::M; }

// The largest penalty any path of the problem can have: every base a gap of its own.
inline std::int64_t wf_ceiling(std::size_t n, std::size_t m, const Penalties& p) {
  return static_cast<std::int64_t>(n + m + 2) * std::max<std::int64_t>(p.mismatch, p.gap_open + p.gap_extend);
}

// Runs a direction to its end. Global: the end (n, m), reached in any component, or, when
// `end` is an open gap, also in that gap with its opening refunded. Ends-free: the last row, at
// the smallest target end. Gives up past max_penalty (not found).
inline WfEnd wf_run(WfEngine& e, bool ends_free, WfBegin end, std::int64_t gap_open, std::int64_t max_penalty, std::int64_t ceiling) {
  const std::int32_t n = e.n(), m = e.m(), kend = m - n;
  const int gap = wf_comp_of(end);
  const std::int64_t refund = gap == WF_M ? 0 : gap_open;
  WfEnd best;
  for (;;) {
    // A later penalty s gives ends of at least s − refund.
    if (best.found && e.score() + 1 - refund >= best.score) return best;
    if (max_penalty >= 0 && e.score() + 1 - refund > max_penalty) return best;
    if (e.score() > ceiling) throw std::logic_error("genomics::align::wavefront: no path reached the end");
    const WfSlot* w = e.next();
    if (!w) continue;
    const std::int64_t s = e.score();
    if (ends_free) {
      const WfOffset* M = w->comp(WF_M);
      for (std::int32_t k = w->lo; k <= w->hi; k++) {
        const WfOffset h = M[k - w->base];
        if (h >= 0 && h - k == n) return {true, s, s, k, h, WF_M};
      }
      continue;
    }
    if (gap != WF_M && w->at(gap, kend) == m && (!best.found || s - refund < best.score)) best = {true, s - refund, s, kend, m, gap};
    if (w->at(WF_M, kend) == m && (!best.found || s < best.score)) best = {true, s, s, kend, m, WF_M};
  }
}

// The path into a run's end, read back through the kept wavefronts; sets the start offset.
inline Cigar wf_backtrace(const WfEngine& e, const Penalties& p, WfBegin begin, bool ends_free, const WfEnd& end, std::size_t& target_begin) {
  const std::int64_t x = p.mismatch, oe = p.gap_open + p.gap_extend, ge = p.gap_extend;
  const std::int32_t n = e.n(), m = e.m();
  Cigar rev;
  std::int64_t s = end.raw;
  std::int32_t k = end.k;
  WfOffset h = end.h;
  int comp = end.comp;
  const auto get = [&](std::int64_t score, int c, std::int32_t kk) -> WfOffset {
    const WfSlot* w = score >= 0 ? e.at(score) : nullptr;
    return w ? w->at(c, kk) : WF_NULL;
  };
  for (;;) {
    if (comp == WF_M) {
      // The beginning: (0, 0), or the first row in ends-free mode.
      WfOffset init = WF_NULL;
      if (s == 0) init = ends_free ? (k >= 0 && k <= m ? k : WF_NULL) : (k == 0 ? 0 : WF_NULL);
      WfOffset mm = WF_NULL;
      const WfOffset a = get(s - x, WF_M, k);
      if (a >= 0 && a + 1 <= m && a + 1 - k <= n) mm = a + 1;
      const WfOffset ins = get(s, WF_I, k), del = get(s, WF_D, k);
      const WfOffset h0 = std::max(std::max(init, mm), std::max(ins, del));
      if (h0 < 0 || h0 > h) throw std::logic_error("genomics::align::wavefront: the backtrace lost its path");
      push(rev, 'M', static_cast<std::uint32_t>(h - h0));
      h = h0;
      if (init == h0) break;
      if (mm == h0) { push(rev, 'M'); s -= x; h = h0 - 1; continue; }
      comp = ins == h0 ? WF_I : WF_D;
      continue;
    }
    if (comp == WF_I) {
      if (begin == WfBegin::I_open && s == 0 && k == 0 && h == 0) break;  // the open gap of the beginning
      push(rev, 'I');
      if (get(s - oe, WF_M, k + 1) == h) { comp = WF_M; s -= oe; k++; continue; }
      if (get(s - ge, WF_I, k + 1) == h) { s -= ge; k++; continue; }
      throw std::logic_error("genomics::align::wavefront: the backtrace lost its gap (I)");
    }
    if (begin == WfBegin::D_open && s == 0 && k == 0 && h == 0) break;
    push(rev, 'D');
    if (get(s - oe, WF_M, k - 1) == h - 1) { comp = WF_M; s -= oe; k--; h--; continue; }
    if (get(s - ge, WF_D, k - 1) == h - 1) { s -= ge; k--; h--; continue; }
    throw std::logic_error("genomics::align::wavefront: the backtrace lost its gap (D)");
  }
  target_begin = static_cast<std::size_t>(h);
  return reversed(std::move(rev));
}

// The penalty and a path of the sub-problem q[0, n) × t[0, m) between the states `begin` and
// `end` (plain or an open gap at each), by the WFA with its traceback.
inline Cigar wf_solve(const Code* q, std::size_t n, const Code* t, std::size_t m, const Penalties& p, WfBegin begin, WfBegin end, WavefrontStats* stats,
                      std::int64_t* penalty) {
  const auto qp = wf_padded(q, n, false, 0xF0, 0xF1), tp = wf_padded(t, m, false, 0xF2, 0xF3);
  WfEngine e(qp, static_cast<std::int32_t>(n), tp, static_cast<std::int32_t>(m), p, begin, false, 0, WavefrontOptions{}, stats);
  const WfEnd stop = wf_run(e, false, end, p.gap_open, -1, wf_ceiling(n, m, p));
  if (penalty) *penalty = stop.score;
  std::size_t tb = 0;
  return wf_backtrace(e, p, begin, false, stop, tb);
}

// A breakpoint of an optimal path of a sub-problem: the optimum, the penalties of the two
// halves (the first including the gap step into the point, when the point is in a gap), the
// point and the component in which the path passes it.
struct WfBreakpoint {
  bool found = false;
  std::int64_t score = 0, left = 0, right = 0;
  std::size_t v = 0, h = 0;
  int comp = WF_M;
};

// BiWFA's search (Marco-Sola et al. 2023): the forward wavefronts from the beginning, the
// reverse ones on the reversed sequences from the end (whose open gap is their beginning), in
// turn, each new wavefront compared with the other direction's last max(x, o + e) + 1.
inline WfBreakpoint wf_breakpoint(const Code* q, std::size_t n, const Code* t, std::size_t m, const Penalties& p, WfBegin begin, WfBegin end,
                                  WavefrontStats* stats) {
  const std::int64_t o = p.gap_open, maxp = std::max<std::int64_t>(p.mismatch, p.gap_open + p.gap_extend);
  const std::size_t window = static_cast<std::size_t>(maxp) + 1;
  const auto qf = wf_padded(q, n, false, 0xF0, 0xF1), tf = wf_padded(t, m, false, 0xF2, 0xF3);
  const auto qr = wf_padded(q, n, true, 0xF0, 0xF1), tr = wf_padded(t, m, true, 0xF2, 0xF3);
  const std::int32_t nn = static_cast<std::int32_t>(n), mm = static_cast<std::int32_t>(m);
  WfEngine F(qf, nn, tf, mm, p, begin, false, window, WavefrontOptions{}, stats);
  WfEngine R(qr, nn, tr, mm, p, end, false, window, WavefrontOptions{}, stats);
  WfBreakpoint best;
  const std::int32_t kk = mm - nn;  // the forward diagonal k is the reverse diagonal kk − k
  const auto overlap = [&](std::int64_t a, std::int64_t b) {
    const WfSlot* f = F.at(a);
    const WfSlot* r = R.at(b);
    if (!f || !r) return;
    const std::int32_t lo = std::max(f->lo, kk - r->hi), hi = std::min(f->hi, kk - r->lo);
    for (int c = WF_M; c <= WF_D; c++) {
      // In a gap both directions paid its opening (or the open gap of their beginning paid none).
      const std::int64_t sum = a + b - (c == WF_M ? 0 : o);
      if (best.found && sum >= best.score) continue;
      const WfOffset* fc = f->comp(c);
      const WfOffset* rc = r->comp(c);
      for (std::int32_t k = lo; k <= hi; k++) {
        const WfOffset hf = fc[k - f->base], hr = rc[kk - k - r->base];
        if (hf >= 0 && hr >= 0 && hf + hr >= mm) {
          best = {true, sum, a, b - (c == WF_M ? 0 : o), static_cast<std::size_t>(hf - k), static_cast<std::size_t>(hf), c};
          break;
        }
      }
    }
  };
  F.next();
  R.next();
  overlap(0, 0);
  const std::int64_t ceiling = wf_ceiling(n, m, p);
  for (;;) {
    if (best.found && F.score() + R.score() >= best.score + o + maxp) break;
    if (F.score() > ceiling) throw std::logic_error("genomics::align::biwfa: the wavefronts never met");
    F.next();
    for (std::int64_t b = std::max<std::int64_t>(0, R.score() - maxp); b <= R.score(); b++) overlap(F.score(), b);
    if (best.found && F.score() + R.score() >= best.score + o + maxp) break;
    R.next();
    for (std::int64_t a = std::max<std::int64_t>(0, F.score() - maxp); a <= F.score(); a++) overlap(a, R.score());
  }
  return best;
}

// BiWFA's recursion: a breakpoint, then each half; halves known to have a penalty of at most
// `small`, or of at most `tiny` bases, by the WFA with its traceback. At a gap breakpoint the
// first half ends with a step of the gap: that step is peeled off, and what precedes it ends
// inside the open gap (its last run refunded, the step paying the opening: o + e in all, or e
// when it extends that run), so that every sub-problem ends freely or in an open gap, the two
// ends whose cost to reach never rises along a diagonal (theory note fast-alignment.md §5.4).
inline void wf_bidirectional(const Code* q, std::size_t n, const Code* t, std::size_t m, const Penalties& p, WfBegin begin, WfBegin end, bool known,
                             std::int64_t hint, Cigar& out, WavefrontStats* stats, int depth, std::int64_t small, std::size_t tiny) {
  const auto solve = [&] {
    const Cigar c = wf_solve(q, n, t, m, p, begin, end, stats, nullptr);
    for (const auto& op : c) push(out, op.op, op.length);
  };
  if ((known && hint <= small) || n + m <= tiny || n == 0 || m == 0 || depth > 200) { solve(); return; }
  const WfBreakpoint bp = wf_breakpoint(q, n, t, m, p, begin, end, stats);
  const WfBegin rbegin = wf_open(bp.comp);
  // A split that leaves the whole problem, with its ends, on one side would recurse forever.
  if (!bp.found || (bp.comp == WF_M && bp.v == n && bp.h == m && end == WfBegin::M) || (bp.v == 0 && bp.h == 0 && rbegin == begin)) { solve(); return; }
  const std::int64_t oe = p.gap_open + p.gap_extend;
  if (bp.comp == WF_M) {
    wf_bidirectional(q, bp.v, t, bp.h, p, begin, WfBegin::M, true, bp.left, out, stats, depth + 1, small, tiny);
  } else if (bp.comp == WF_I) {
    wf_bidirectional(q, bp.v - 1, t, bp.h, p, begin, WfBegin::I_open, true, bp.left - oe, out, stats, depth + 1, small, tiny);
    push(out, 'I');
  } else {
    wf_bidirectional(q, bp.v, t, bp.h - 1, p, begin, WfBegin::D_open, true, bp.left - oe, out, stats, depth + 1, small, tiny);
    push(out, 'D');
  }
  wf_bidirectional(q + bp.v, n - bp.v, t + bp.h, m - bp.h, p, rbegin, end, true, bp.right, out, stats, depth + 1, small, tiny);
}

// Halves of at most this penalty are aligned by the WFA with its traceback: about
// (s/e)² offsets, a few hundred kilobytes.
inline std::int64_t wf_small(const Penalties& p) {
  return std::max<std::int64_t>(256, 4 * (std::max<std::int64_t>(p.mismatch, p.gap_open + p.gap_extend) + p.gap_open) + 8);
}

}  // namespace detail

// An optimal alignment under the penalties, global or semiglobal (ends-free on the target),
// with its path; score = −penalty. Empty when options.max_penalty is passed. With `record`,
// every wavefront's M̃ offsets and the work done are kept, for drawing.
inline std::optional<Alignment> wavefront(const Codes& query, const Codes& target, const Penalties& p, Mode mode = Mode::global, const WavefrontOptions& options = {},
                                          WavefrontStats* stats = nullptr, WavefrontRecord* record = nullptr) {
  detail::wf_validate(query, target, p, "genomics::align::wavefront");
  if (mode != Mode::global && mode != Mode::semiglobal) throw std::invalid_argument("genomics::align::wavefront: global or semiglobal");
  const bool ends_free = mode == Mode::semiglobal;
  const auto qp = detail::wf_padded(query.data(), query.size(), false, 0xF0, 0xF1), tp = detail::wf_padded(target.data(), target.size(), false, 0xF2, 0xF3);
  detail::WfEngine e(qp, static_cast<std::int32_t>(query.size()), tp, static_cast<std::int32_t>(target.size()), p, detail::WfBegin::M, ends_free, 0, options, stats,
                     record);
  const detail::WfEnd end =
      detail::wf_run(e, ends_free, detail::WfBegin::M, p.gap_open, options.max_penalty, detail::wf_ceiling(query.size(), target.size(), p));
  if (!end.found) return std::nullopt;
  Alignment a;
  std::size_t tb = 0;
  a.cigar = detail::wf_backtrace(e, p, detail::WfBegin::M, ends_free, end, tb);
  a.score = -end.score;
  a.query_begin = 0;
  a.query_end = query.size();
  a.target_begin = tb;
  a.target_end = static_cast<std::size_t>(end.h);
  return a;
}

// The optimal penalty alone, keeping only the last max(x, o + e) + 1 wavefronts: O(s) memory.
// −1 when options.max_penalty is passed.
inline std::int64_t wavefront_penalty(const Codes& query, const Codes& target, const Penalties& p, Mode mode = Mode::global, const WavefrontOptions& options = {},
                                      WavefrontStats* stats = nullptr) {
  detail::wf_validate(query, target, p, "genomics::align::wavefront_penalty");
  if (mode != Mode::global && mode != Mode::semiglobal) throw std::invalid_argument("genomics::align::wavefront_penalty: global or semiglobal");
  const bool ends_free = mode == Mode::semiglobal;
  const auto qp = detail::wf_padded(query.data(), query.size(), false, 0xF0, 0xF1), tp = detail::wf_padded(target.data(), target.size(), false, 0xF2, 0xF3);
  const std::size_t window = static_cast<std::size_t>(std::max(p.mismatch, p.gap_open + p.gap_extend)) + 1;
  detail::WfEngine e(qp, static_cast<std::int32_t>(query.size()), tp, static_cast<std::int32_t>(target.size()), p, detail::WfBegin::M, ends_free, window, options,
                     stats);
  const detail::WfEnd end = detail::wf_run(e, ends_free, detail::WfBegin::M, p.gap_open, options.max_penalty, detail::wf_ceiling(query.size(), target.size(), p));
  return end.found ? end.score : -1;
}

// An optimal global alignment by BiWFA, in O(s) memory; score = −penalty.
inline Alignment biwfa(const Codes& query, const Codes& target, const Penalties& p, WavefrontStats* stats = nullptr) {
  detail::wf_validate(query, target, p, "genomics::align::biwfa");
  Alignment a;
  detail::wf_bidirectional(query.data(), query.size(), target.data(), target.size(), p, detail::WfBegin::M, detail::WfBegin::M, false, 0, a.cigar, stats, 0,
                           detail::wf_small(p), 64);
  a.query_end = query.size();
  a.target_end = target.size();
  a.score = -penalty(a, query, target, p);
  return a;
}

// ---- penalties and scores -------------------------------------------------------------------------
// The penalties of a scoring whose matrix has one match score α and one mismatch score −β on A,
// C, G, T (α > −β, 2r + α > 0), reduced by their gcd (Eizenga and Paten 2022).
inline PenaltyForm penalty_form(const Scoring& s) {
  if (s.matrix.size < 4) throw std::invalid_argument("genomics::align::penalty_form: a DNA matrix");
  const int alpha = s.matrix(0, 0), mis = s.matrix(0, 1);
  for (Code a = 0; a < 4; a++)
    for (Code b = 0; b < 4; b++)
      if (s.matrix(a, b) != (a == b ? alpha : mis)) throw std::invalid_argument("genomics::align::penalty_form: the matrix is not one match and one mismatch score on ACGT");
  const std::int64_t beta = -static_cast<std::int64_t>(mis);
  std::int64_t x = 2 * (alpha + beta), o = 2 * static_cast<std::int64_t>(s.gap_open), e = 2 * static_cast<std::int64_t>(s.gap_extend) + alpha;
  if (x <= 0 || e <= 0 || o < 0) throw std::invalid_argument("genomics::align::penalty_form: the scoring has no penalties (match <= mismatch, or 2 extend + match <= 0)");
  const std::int64_t g = std::gcd(std::gcd(x, e), o);
  PenaltyForm f;
  f.penalties = {static_cast<int>(x / g), static_cast<int>(o / g), static_cast<int>(e / g)};
  f.match = alpha;
  f.scale = g;
  return f;
}
// The score of a global alignment from its penalty under the form.
inline std::int64_t score_from_penalty(const PenaltyForm& f, std::int64_t penalty, std::size_t query_length, std::size_t target_length) {
  return (f.match * static_cast<std::int64_t>(query_length + target_length) - f.scale * penalty) / 2;
}
inline std::int64_t penalty_from_score(const PenaltyForm& f, std::int64_t score, std::size_t query_length, std::size_t target_length) {
  return (f.match * static_cast<std::int64_t>(query_length + target_length) - 2 * score) / f.scale;
}

// An optimal global alignment under a scoring, by the WFA (or BiWFA) on its penalties; the
// score is the scoring's. Refuses sequences with N when the matrix scores N otherwise than a
// mismatch, for which the two orders differ.
inline Alignment wavefront_global(const Codes& query, const Codes& target, const Scoring& s, bool bidirectional = false, WavefrontStats* stats = nullptr) {
  const PenaltyForm f = penalty_form(s);
  if (sequence::count_ambiguous(query) + sequence::count_ambiguous(target) > 0)
    for (Code a = 0; a < s.matrix.size; a++)
      for (Code b = 0; b < s.matrix.size; b++)
        if ((a >= 4 || b >= 4) && s.matrix(a, b) != s.matrix(0, 1))
          throw std::invalid_argument("genomics::align::wavefront_global: N is scored otherwise than a mismatch; the penalty form would not be exact");
  Alignment a = bidirectional ? biwfa(query, target, f.penalties, stats) : *wavefront(query, target, f.penalties, Mode::global, {}, stats);
  a.score = score_from_penalty(f, -a.score, query.size(), target.size());
  return a;
}

}  // namespace genomics::align
