// Checks of the wavefront algorithm: WFA2-lib's published example; the optimal penalty, global
// and ends-free, equal to Gotoh's recurrences (gotoh.hpp, on the matrix of the penalties) on
// pairs of lengths 0 to 3000 at divergences from 0 to unrelated, under several sets of
// penalties, with every CIGAR holding the penalty it claims; with unit penalties the edit
// distance of Myers' algorithm; the score-only run in O(s) memory equal to the full one;
// BiWFA equal to the WFA everywhere, on long similar sequences too, in a fraction of its
// memory; the penalty form of a scoring an identity on every global path, and the WFA on it
// equal to Gotoh's optimum under the scoring; the adaptive heuristic valid, never better than
// the optimum, and equal to it on similar sequences.

#include <genomics/align/bitparallel.hpp>
#include <genomics/align/gotoh.hpp>
#include <genomics/align/wavefront.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::align;

namespace {

constexpr std::uint32_t NOISE = 0x57464e53;

Codes noisy(const Codes& x, double rate, std::uint64_t seed, std::uint32_t object, double gap_share = 2.0 / 3) {
  calx::random::Stream r(seed, NOISE, object, 0);
  Codes y;
  y.reserve(x.size() + x.size() / 8 + 4);
  for (std::size_t i = 0; i < x.size(); i++) {
    const double u = r.uniform();
    if (u < rate * (1 - gap_share)) { y.push_back(simulate::substitute(r, x[i])); continue; }
    if (u < rate * (1 - gap_share / 2)) {  // a deletion of 1 to 4 bases
      i += r.below(4);
      continue;
    }
    if (u < rate) {
      const std::size_t len = 1 + r.below(4);
      for (std::size_t k = 0; k < len; k++) y.push_back(static_cast<Code>(r.below(4)));
    }
    y.push_back(x[i]);
  }
  return y;
}

// The scoring whose optimum is minus the penalty: 0 for a match of A, C, G, T, −x otherwise.
Scoring scoring_of(const Penalties& p) {
  Scoring s;
  s.matrix = Matrix::dna(0, -p.mismatch, -p.mismatch);
  s.gap_open = p.gap_open;
  s.gap_extend = p.gap_extend;
  return s;
}

// The least penalty of a path of q × t between the states asked, by Gotoh's three-state
// recurrences with the boundary set accordingly: X(i, j) the least penalty of a path into
// (i, j) whose last step is X; the empty path at (0, 0); an open gap at the beginning as the
// state I or D at (0, 0) with penalty 0; an open gap at the end as the refund of an opening to
// a path that ends in it.
std::int64_t constrained(const Codes& q, const Codes& t, const Penalties& p, detail::WfBegin begin, detail::WfBegin end) {
  const std::int64_t INF = std::int64_t{1} << 50, oe = p.gap_open + p.gap_extend, e = p.gap_extend;
  const std::size_t n = q.size(), m = t.size(), W = m + 1;
  std::vector<std::int64_t> M((n + 1) * W, INF), I((n + 1) * W, INF), D((n + 1) * W, INF), H((n + 1) * W, INF);
  I[0] = begin == detail::WfBegin::I_open ? 0 : INF;
  D[0] = begin == detail::WfBegin::D_open ? 0 : INF;
  H[0] = 0;
  for (std::size_t i = 0; i <= n; i++)
    for (std::size_t j = 0; j <= m; j++) {
      if (i == 0 && j == 0) continue;
      const std::size_t c = i * W + j;
      if (i && j) M[c] = H[c - W - 1] + ((q[i - 1] == t[j - 1] && q[i - 1] < 4) ? 0 : p.mismatch);
      if (i) I[c] = std::min(H[c - W] + oe, I[c - W] + e);
      if (j) D[c] = std::min(H[c - 1] + oe, D[c - 1] + e);
      H[c] = std::min({M[c], I[c], D[c]});
    }
  const std::size_t last = n * W + m;
  if (end == detail::WfBegin::I_open) return std::min(H[last], I[last] - p.gap_open);
  if (end == detail::WfBegin::D_open) return std::min(H[last], D[last] - p.gap_open);
  return H[last];
}

// The penalty of a CIGAR of a sub-problem: its first run free when it continues the open gap
// of the beginning, its last run refunded when the end's gap is open after it (an empty path
// being in the beginning's state); −2^40 when it does not span q × t.
std::int64_t constrained_penalty(const Cigar& c, const Codes& q, const Codes& t, const Penalties& p, detail::WfBegin begin, detail::WfBegin end) {
  if (query_length(c) != q.size() || target_length(c) != t.size()) return -(std::int64_t{1} << 40);
  Alignment a;
  a.cigar = c;
  a.query_end = q.size();
  a.target_end = t.size();
  std::int64_t pen = penalty(a, q, t, p);
  const char opened = begin == detail::WfBegin::I_open ? 'I' : begin == detail::WfBegin::D_open ? 'D' : 0;
  const char open_after = end == detail::WfBegin::I_open ? 'I' : end == detail::WfBegin::D_open ? 'D' : 0;
  if (!c.empty() && opened && c.front().op == opened) pen -= p.gap_open;
  const char last = c.empty() ? opened : c.back().op;
  if (open_after && last == open_after) pen -= p.gap_open;
  return pen;
}


bool spans(const Alignment& a, const Codes& q, const Codes& t) {
  return a.query_begin == 0 && a.query_end == q.size() && a.query_begin + query_length(a.cigar) == a.query_end && a.target_begin + target_length(a.cigar) == a.target_end &&
         a.target_end <= t.size();
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const Penalties standard{4, 6, 2};

  // ---- WFA2-lib's example ----------------------------------------------------------------------------
  {
    // The example pair and penalties of WFA2-lib's README (Marco-Sola et al., github.com/smarco/
    // WFA2-lib, §2.1, retrieved 2026-10-24: mismatch 4, gap opening 6, extension 2). The README
    // prints no score; its optimum is 24 (two mismatches and two gaps of one base), which the
    // three algorithms must agree on.
    const Codes q = sequence::encode_dna("TCTTTACTCGCGCGTTGGAGAAATACAATAGT"), t = sequence::encode_dna("TCTATACTGCGCGTTTGGAGAAATAAAATAGT");
    const auto a = wavefront(q, t, standard);
    const auto b = biwfa(q, t, standard);
    s.require("WFA2-lib's example aligns with penalty 24, by the WFA, by BiWFA and by Gotoh's recurrences",
              a && a->score == -24 && penalty(*a, q, t, standard) == 24 && b.score == -24 && gotoh_score(q, t, scoring_of(standard), Mode::global) == -24,
              a ? cigar_string(explicit_matches(a->cigar, q, t)) : "none");
  }

  // ---- against Gotoh's recurrences ---------------------------------------------------------------------
  const std::vector<Penalties> sets = {standard, {1, 0, 1}, {2, 3, 1}, {6, 5, 3}, {5, 10, 1}, {3, 1, 4}, {9, 2, 2}};
  const std::vector<std::size_t> lengths = {0, 1, 2, 7, 33, 100, 257, 600, 1500, 3000};
  const std::vector<double> rates = {0.0, 0.01, 0.05, 0.1, 0.2, 1.0};
  int global_ok = 0, semi_ok = 0, cigar_ok = 0, score_ok = 0, bi_ok = 0, unit_ok = 0, total = 0, semi_total = 0, unit_total = 0;
  std::string bad;
  std::uint32_t object = 0;
  for (const auto& p : sets) {
    const Scoring sc = scoring_of(p);
    for (std::size_t m : lengths)
      for (double rate : rates) {
        object++;
        if (m > 1500 && (p.mismatch != 4 || rate > 0.1)) continue;  // the O(mn) reference is the cost
        Codes q = simulate::random_sequence(m, 0.5, 41, object);
        Codes t = rate >= 1.0 ? simulate::random_sequence(m + m / 5, 0.5, 42, object) : noisy(q, rate, 43, object);
        if (object % 5 == 0 && !q.empty()) q[q.size() / 2] = sequence::N;  // N matches nothing
        total++;
        const auto a = wavefront(q, t, p);
        const std::int64_t want = -gotoh_score(q, t, sc, Mode::global);
        if (a && -a->score == want) global_ok++;
        else if (bad.empty()) bad = "global m=" + std::to_string(m) + " rate=" + std::to_string(rate) + ": " + std::to_string(a ? -a->score : -1) + " vs " + std::to_string(want);
        if (a && spans(*a, q, t) && a->target_begin == 0 && a->target_end == t.size() && penalty(*a, q, t, p) == -a->score) cigar_ok++;
        else if (bad.empty()) bad = "global CIGAR m=" + std::to_string(m);
        if (wavefront_penalty(q, t, p) == want) score_ok++;
        const Alignment b = biwfa(q, t, p);
        if (-b.score == want && spans(b, q, t) && b.target_end == t.size()) bi_ok++;
        else if (bad.empty()) bad = "biwfa m=" + std::to_string(m) + " rate=" + std::to_string(rate) + ": " + std::to_string(-b.score) + " vs " + std::to_string(want);
        if (p.mismatch == 1 && p.gap_open == 0 && p.gap_extend == 1) {
          unit_total++;
          if (want == myers_distance(q, t)) unit_ok++;
        }
        // ends-free: the query placed in a longer target
        if (m <= 1500) {
          Codes text = simulate::random_sequence(m / 2 + 9, 0.5, 44, object);
          text.insert(text.end(), t.begin(), t.end());
          const Codes tail = simulate::random_sequence(m / 3 + 5, 0.5, 45, object);
          text.insert(text.end(), tail.begin(), tail.end());
          semi_total++;
          const auto e = wavefront(q, text, p, Mode::semiglobal);
          const Endpoint ge = gotoh_end(q, text, sc, Mode::semiglobal);
          if (e && -e->score == -ge.score && e->target_end == ge.target_end && spans(*e, q, text) && penalty(*e, q, text, p) == -e->score &&
              wavefront_penalty(q, text, p, Mode::semiglobal) == -ge.score)
            semi_ok++;
          else if (bad.empty()) bad = "semiglobal m=" + std::to_string(m) + " rate=" + std::to_string(rate);
        }
      }
  }
  s.require("the WFA's global penalty equals Gotoh's optimum on " + std::to_string(total) + " pairs under 7 sets of penalties", global_ok == total, bad);
  s.require("every global CIGAR spans both sequences and holds the penalty it claims", cigar_ok == total, bad);
  s.require("the ends-free penalty and the smallest optimal target end equal Gotoh's semiglobal ones (" + std::to_string(semi_total) + " pairs), their CIGARs exact",
            semi_ok == semi_total, bad);
  s.require("the score-only WFA, in O(s) memory, gives the same penalties", score_ok == total);
  s.require("BiWFA's penalty equals the optimum on every pair", bi_ok == total, bad);
  s.require("with mismatch 1, opening 0 and extension 1, the penalty is Myers' edit distance (" + std::to_string(unit_total) + " pairs)", unit_ok == unit_total);

  // ---- the sub-problems of BiWFA: every beginning and end, and recursion to single bases ------------
  {
    using detail::WfBegin;
    const WfBegin states[] = {WfBegin::M, WfBegin::I_open, WfBegin::D_open};
    const std::int64_t never = -(std::int64_t{1} << 40);
    calx::random::Stream r(81, NOISE, 1, 0);
    int solve_ok = 0, split_ok = 0, deep_ok = 0, cases = 0;
    std::string why;
    const auto name = [](std::size_t n, std::size_t m, WfBegin b, WfBegin e) {
      return "n=" + std::to_string(n) + " m=" + std::to_string(m) + " begin=" + std::to_string(int(b)) + " end=" + std::to_string(int(e));
    };
    for (int trial = 0; trial < 900; trial++) {
      const std::size_t n = r.below(trial < 600 ? 9 : 60), m = r.below(trial < 600 ? 9 : 60);
      Codes q(n), t(m);
      for (auto& c : q) c = static_cast<Code>(r.below(trial % 3 == 0 ? 2 : 4));
      for (auto& c : t) c = static_cast<Code>(r.below(trial % 3 == 0 ? 2 : 4));
      const Penalties& p = sets[static_cast<std::size_t>(trial) % sets.size()];
      const std::int64_t oe = p.gap_open + p.gap_extend;
      for (WfBegin b : states)
        for (WfBegin e : states) {
          cases++;
          const std::int64_t want = constrained(q, t, p, b, e);
          std::int64_t got = never;
          const Cigar c = detail::wf_solve(q.data(), n, t.data(), m, p, b, e, nullptr, &got);
          if (got == want && constrained_penalty(c, q, t, p, b, e) == want) solve_ok++;
          else if (why.empty()) why = "solve " + name(n, m, b, e) + ": " + std::to_string(got) + " vs " + std::to_string(want);
          if (n + m >= 1) {
            const auto bp = detail::wf_breakpoint(q.data(), n, t.data(), m, p, b, e, nullptr);
            // The halves, solved alone with the states the breakpoint gives them, cost what it says;
            // at a gap the first half is its peeled last step and what precedes it, open after.
            const WfBegin rbegin = detail::wf_open(bp.comp);
            const Codes qr(q.begin() + static_cast<std::ptrdiff_t>(bp.v), q.end()), tr(t.begin() + static_cast<std::ptrdiff_t>(bp.h), t.end());
            bool left_ok;
            if (bp.v == 0 && bp.h == 0 && rbegin == b) left_ok = bp.left == 0;  // the beginning's own open gap
            else {
              const std::size_t v = bp.v - (bp.comp == detail::WF_I ? 1 : 0), h = bp.h - (bp.comp == detail::WF_D ? 1 : 0);
              const Codes ql(q.begin(), q.begin() + static_cast<std::ptrdiff_t>(v)), tl(t.begin(), t.begin() + static_cast<std::ptrdiff_t>(h));
              left_ok = bp.comp == detail::WF_M ? constrained(ql, tl, p, b, WfBegin::M) == bp.left : constrained(ql, tl, p, b, rbegin) + oe == bp.left;
            }
            if (bp.found && bp.score == want && left_ok && bp.left + bp.right == want && constrained(qr, tr, p, rbegin, e) == bp.right) split_ok++;
            else if (why.empty()) why = "split " + name(n, m, b, e);
          } else split_ok++;
          Cigar deep;
          detail::wf_bidirectional(q.data(), n, t.data(), m, p, b, e, false, 0, deep, nullptr, 0, -1, 1);
          if (constrained_penalty(deep, q, t, p, b, e) == want) deep_ok++;
          else if (why.empty()) why = "deep " + name(n, m, b, e);
        }
    }
    s.require("the WFA solves every sub-problem, each end plain or inside an open gap, at the three-state recurrence's optimum (" + std::to_string(cases) + ")",
              solve_ok == cases, why);
    s.require("every breakpoint has the optimal penalty, and its halves alone cost what it says they do", split_ok == cases, why);
    s.require("BiWFA recursed down to single bases gives an optimal path of every sub-problem", deep_ok == cases, why);
  }

  // ---- the record of the wavefronts, for drawing -----------------------------------------------------------
  {
    const Codes q = simulate::random_sequence(800, 0.5, 91), t = noisy(q, 0.05, 92, 0);
    WavefrontRecord rec;
    const auto plain = wavefront(q, t, standard), drawn = wavefront(q, t, standard, Mode::global, {}, nullptr, &rec);
    bool increasing = true, consistent = true;
    for (std::size_t i = 0; i < rec.fronts.size(); i++) {
      const auto& f = rec.fronts[i];
      if (i && f.score <= rec.fronts[i - 1].score) increasing = false;
      if (f.offsets.size() != static_cast<std::size_t>(f.hi - f.lo + 1)) consistent = false;
    }
    const auto& last = rec.fronts.back();
    const std::int64_t kend = static_cast<std::int64_t>(t.size()) - static_cast<std::int64_t>(q.size());
    s.require("a recorded run gives the same alignment; its fronts rise by penalty and the last reaches (n, m)",
              plain && drawn && plain->cigar == drawn->cigar && increasing && consistent && last.score == -drawn->score && kend >= last.lo && kend <= last.hi &&
                  last.offsets[static_cast<std::size_t>(kend - last.lo)] == static_cast<std::int64_t>(t.size()));
    s.require("and it touched a small part of the matrix: offsets and extended cells under 5% of (n + 1)(m + 1)",
              (rec.offsets_computed + rec.cells_extended) * 20 < (q.size() + 1) * (t.size() + 1),
              std::to_string(rec.offsets_computed) + " offsets, " + std::to_string(rec.cells_extended) + " extended");
  }

  // ---- long similar sequences: BiWFA against the WFA, and its memory ---------------------------------------
  {
    const Codes q = simulate::random_sequence(30'000, 0.5, 51);
    const Codes t = noisy(q, 0.01, 52, 0);
    WavefrontStats full, bi;
    const auto a = wavefront(q, t, standard, Mode::global, {}, &full);
    const Alignment b = biwfa(q, t, standard, &bi);
    s.require("on two 30 kb sequences at 1% divergence BiWFA's alignment has the WFA's optimal penalty",
              a && b.score == a->score && penalty(b, q, t, standard) == -b.score && spans(b, q, t), std::to_string(-b.score) + " vs " + std::to_string(a ? -a->score : -1));
    s.require("and holds at most a tenth of the WFA's wavefront memory", bi.peak_bytes * 10 <= full.peak_bytes,
              std::to_string(bi.peak_bytes) + " B vs " + std::to_string(full.peak_bytes) + " B");
    // Long and more divergent: the halves after a gap breakpoint start inside the gap.
    const Codes q2 = simulate::random_sequence(12'000, 0.5, 53);
    const Codes t2 = noisy(q2, 0.06, 54, 0, 0.9);
    const std::int64_t want = wavefront_penalty(q2, t2, standard);
    const Alignment b2 = biwfa(q2, t2, standard);
    s.require("on 12 kb at 6% divergence, mostly indels, BiWFA's penalty is the optimum", -b2.score == want && penalty(b2, q2, t2, standard) == want && spans(b2, q2, t2),
              std::to_string(-b2.score) + " vs " + std::to_string(want));
  }

  // ---- penalties and scores ------------------------------------------------------------------------------
  {
    struct Case { int match, mismatch, open, extend; };
    const std::vector<Case> cases = {{2, -4, 4, 2}, {1, -1, 0, 1}, {2, -3, 5, 2}, {5, -4, 10, 1}, {1, -4, 6, 1}, {3, -2, 1, 3}};
    int identity_ok = 0, identity_total = 0, opt_ok = 0, opt_total = 0;
    std::string why;
    std::uint32_t obj = 1000;
    for (const auto& c : cases) {
      Scoring sc;
      sc.matrix = Matrix::dna(c.match, c.mismatch, c.mismatch);
      sc.gap_open = c.open;
      sc.gap_extend = c.extend;
      const PenaltyForm f = penalty_form(sc);
      for (std::size_t m : {1, 10, 80, 300, 900})
        for (double rate : {0.02, 0.1, 0.3}) {
          obj++;
          const Codes q = simulate::random_sequence(m, 0.5, 61, obj), t = noisy(q, rate, 62, obj);
          // The identity on a path that is not optimal for either order: Gotoh's under another scoring.
          Scoring other = sc;
          other.gap_open = c.open + 7;
          for (const Scoring* by : {&sc, &other}) {
            const Alignment g = gotoh(q, t, *by, Mode::global);
            identity_total++;
            if (rescore(g, q, t, sc) == score_from_penalty(f, penalty(g, q, t, f.penalties), q.size(), t.size())) identity_ok++;
          }
          opt_total++;
          const Alignment w = wavefront_global(q, t, sc), wb = wavefront_global(q, t, sc, true);
          const std::int64_t want = gotoh_score(q, t, sc, Mode::global);
          if (w.score == want && rescore(w, q, t, sc) == want && wb.score == want && rescore(wb, q, t, sc) == want) opt_ok++;
          else if (why.empty()) why = "match " + std::to_string(c.match) + " m=" + std::to_string(m) + ": " + std::to_string(w.score) + " vs " + std::to_string(want);
        }
    }
    s.require("S = (α(n + m) − scale · P)/2 holds on every global path of " + std::to_string(identity_total) + " (Eizenga and Paten 2022)", identity_ok == identity_total);
    s.require("the WFA and BiWFA on the penalty form find Gotoh's global optimum under 6 scorings (" + std::to_string(opt_total) + " pairs)", opt_ok == opt_total, why);
    const PenaltyForm f = penalty_form(Scoring{});
    s.require("the default scoring (2, −4, 4, 2) is the penalties (6, 4, 3) at scale 2", f.penalties.mismatch == 6 && f.penalties.gap_open == 4 && f.penalties.gap_extend == 3 && f.scale == 2);
  }

  // ---- the adaptive heuristic ------------------------------------------------------------------------------
  {
    WavefrontOptions adaptive;
    adaptive.adaptive = true;
    int valid = 0, equal = 0, similar = 0, n = 0;
    for (std::uint32_t i = 0; i < 40; i++) {
      const double rate = i < 20 ? 0.01 : 0.15;
      const Codes q = simulate::random_sequence(2000 + 50 * i, 0.5, 71, i), t = noisy(q, rate, 72, i);
      const std::int64_t opt = wavefront_penalty(q, t, standard);
      const auto a = wavefront(q, t, standard, Mode::global, adaptive);
      n++;
      if (a && penalty(*a, q, t, standard) == -a->score && -a->score >= opt && spans(*a, q, t)) valid++;
      if (i < 20) { similar++; equal += a && -a->score == opt; }
    }
    s.require("WF-adaptive's alignments are valid and never better than the optimum (40 pairs at 1% and 15%)", valid == n);
    s.require("and on the 20 pairs at 1% divergence it finds the optimum", equal == similar, std::to_string(equal) + " of " + std::to_string(similar));
  }

  return s.done();
}
