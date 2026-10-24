// Checks of xdrop.hpp: with X = ∞ the gapped extension to the right (left) is the best alignment
// anchored at its start (end) with the other end free, by brute force over every prefix (suffix)
// pair; for every X the extension never beats that, begins (ends) at the anchor and rescores to its
// claim; the joined extension spans the anchor; the ungapped extension is the best prefix sum on
// each side with X = ∞ and grows with X; a planted homologous core between random flanks is
// recovered from a seed in its middle to within a few bases of its ends; and the Z-drop splits an
// alignment at a divergent block where the X-drop alone crosses it.

#include <genomics/align/gotoh.hpp>
#include <genomics/align/scoring.hpp>
#include <genomics/align/xdrop.hpp>
#include <genomics/simulate/family.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::align;

namespace {

Codes random_codes(calx::random::Stream& r, std::size_t n) {
  Codes x(n);
  for (auto& c : x) c = static_cast<Code>(r.below(4));
  return x;
}

Codes slice(const Codes& x, std::size_t b, std::size_t e) { return Codes(x.begin() + static_cast<std::ptrdiff_t>(b), x.begin() + static_cast<std::ptrdiff_t>(e)); }

// The best global score of q[qpos, qpos + i) against t[tpos, tpos + j) over all i, j.
std::int64_t anchored_right(const Codes& q, const Codes& t, std::size_t qpos, std::size_t tpos, const Scoring& s) {
  std::int64_t best = 0;
  for (std::size_t i = qpos; i <= q.size(); i++)
    for (std::size_t j = tpos; j <= t.size(); j++) best = std::max(best, gotoh_score(slice(q, qpos, i), slice(t, tpos, j), s, Mode::global));
  return best;
}
std::int64_t anchored_left(const Codes& q, const Codes& t, std::size_t qpos, std::size_t tpos, const Scoring& s) {
  std::int64_t best = 0;
  for (std::size_t i = 0; i <= qpos; i++)
    for (std::size_t j = 0; j <= tpos; j++) best = std::max(best, gotoh_score(slice(q, i, qpos), slice(t, j, tpos), s, Mode::global));
  return best;
}

constexpr std::int64_t INF = 1'000'000'000;

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x78647270, 0, 0);
  const Scoring sc = scoring(Matrix::dna(2, -3, -1), 5, 2);

  {
    bool right = true, left = true, bounded = true, anchored = true, joined = true;
    for (int rep = 0; rep < 120; rep++) {
      // Related pairs, so that extensions are long: a sequence and a mutated copy.
      const Codes root = random_codes(r, 6 + r.below(18));
      const auto f = simulate::evolve_family(root, 2, 0.15, 0.08, 500 + static_cast<std::uint64_t>(rep));
      const Codes& q = f.sequences[0];
      const Codes& t = f.sequences[1];
      const std::size_t qpos = r.below(q.size() + 1), tpos = r.below(t.size() + 1);
      const Alignment er = extend_right(q, t, qpos, tpos, sc, {INF, -1}), el = extend_left(q, t, qpos, tpos, sc, {INF, -1});
      right = right && er.score == anchored_right(q, t, qpos, tpos, sc);
      left = left && el.score == anchored_left(q, t, qpos, tpos, sc);
      for (std::int64_t x : {0, 3, 8, 20}) {
        const Alignment a = extend_right(q, t, qpos, tpos, sc, {x, -1}), b = extend_left(q, t, qpos, tpos, sc, {x, -1});
        bounded = bounded && a.score <= er.score && b.score <= el.score && a.score >= 0 && b.score >= 0;
        anchored = anchored && a.query_begin == qpos && a.target_begin == tpos && b.query_end == qpos && b.target_end == tpos && rescore(a, q, t, sc) == a.score &&
                   rescore(b, q, t, sc) == b.score;
        const Alignment j = xdrop(q, t, qpos, tpos, sc, {x, -1});
        joined = joined && j.query_begin <= qpos && j.query_end >= qpos && j.score >= a.score + b.score && rescore(j, q, t, sc) == j.score;
      }
    }
    s.require("gapped X-drop, X = ∞: the right extension is the best alignment anchored at its start, by brute force over prefixes", right, "120 related pairs of length ≤ 30, random anchors");
    s.require("gapped X-drop, X = ∞: the left extension is the best alignment anchored at its end, by brute force over suffixes", left);
    s.require("gapped X-drop: for X = 0, 3, 8, 20 the extension never beats X = ∞ and is never negative", bounded);
    s.require("gapped X-drop: each extension begins (ends) at the anchor and rescores to its claim", anchored);
    s.require("gapped X-drop: the joined extension spans the anchor and scores at least the sum of its halves", joined);
  }
  {
    bool best = true, monotone = true;
    for (int rep = 0; rep < 200; rep++) {
      const Codes root = random_codes(r, 60);
      const auto f = simulate::evolve_family(root, 2, 0.2, 0.0, 900 + static_cast<std::uint64_t>(rep));
      const Codes& q = f.sequences[0];
      const Codes& t = f.sequences[1];
      const std::size_t len = 4, qpos = 20 + r.below(30), tpos = qpos;
      const Matrix& M = sc.matrix;
      // The best prefix sums of the diagonal on each side, by brute force.
      std::int64_t seed = 0, run = 0, br = 0, bl = 0;
      for (std::size_t k = 0; k < len; k++) seed += M(q[qpos + k], t[tpos + k]);
      for (std::size_t k = qpos + len; k < q.size() && k < t.size(); k++) br = std::max(br, run += M(q[k], t[k]));
      run = 0;
      for (std::size_t k = qpos; k-- > 0;) bl = std::max(bl, run += M(q[k], t[k]));
      best = best && ungapped_xdrop(q, t, qpos, tpos, len, M, INF).score == seed + br + bl;
      std::int64_t last = -INF;
      for (std::int64_t x : {0, 1, 2, 4, 8, 16, 32, 64}) {
        const Alignment a = ungapped_xdrop(q, t, qpos, tpos, len, M, x);
        monotone = monotone && a.score >= last && rescore(a, q, t, sc) == a.score;
        last = a.score;
      }
    }
    s.require("ungapped X-drop, X = ∞: the seed plus the best prefix sum on each side of the diagonal", best, "200 pairs of length 60 with 20% substitutions");
    s.require("ungapped X-drop: the score grows with X (0, 1, 2, …, 64) and rescores to its claim", monotone);
  }
  {
    // A homologous core of 600 bases (two branches of 2% substitutions and 0.5% indels each) between
    // unrelated flanks of 400. Near an end, a last mismatch or indel is kept only when the matches
    // after it repay it, so the ends are found to within a few bases, not exactly.
    bool found = true;
    std::string worst;
    std::size_t slack = 0;
    for (std::uint64_t seed = 1; seed <= 10; seed++) {
      const Codes core = simulate::random_sequence(600, 0.5, seed);
      const auto f = simulate::evolve_family(core, 2, 0.02, 0.005, seed);
      Codes q = simulate::random_sequence(400, 0.5, seed + 1000), t = simulate::random_sequence(400, 0.5, seed + 2000);
      q.insert(q.end(), f.sequences[0].begin(), f.sequences[0].end());
      t.insert(t.end(), f.sequences[1].begin(), f.sequences[1].end());
      const Codes qf = simulate::random_sequence(400, 0.5, seed + 3000), tf = simulate::random_sequence(400, 0.5, seed + 4000);
      q.insert(q.end(), qf.begin(), qf.end());
      t.insert(t.end(), tf.begin(), tf.end());
      // The anchor: the middle of the core, on the true alignment (a column of the family where
      // both members have a base).
      std::size_t col = f.alignment[0].size() / 2, qi = 0, ti = 0;
      while (f.alignment[0][col] == sequence::GAP || f.alignment[1][col] == sequence::GAP) col++;
      for (std::size_t c = 0; c < col; c++) {
        qi += f.alignment[0][c] != sequence::GAP;
        ti += f.alignment[1][c] != sequence::GAP;
      }
      const Alignment a = xdrop(q, t, 400 + qi, 400 + ti, bwa_mem(), {40, -1});
      const std::size_t qend = 400 + f.sequences[0].size(), tend = 400 + f.sequences[1].size();
      const auto off = [](std::size_t u, std::size_t v) { return u > v ? u - v : v - u; };
      const std::size_t err = std::max({off(a.query_begin, 400), off(a.target_begin, 400), off(a.query_end, qend), off(a.target_end, tend)});
      slack = std::max(slack, err);
      found = found && err <= 12;
    }
    s.require("gapped X-drop: a planted core of 600 bases is recovered from its middle to within 12 bases of each end", found, "largest error " + std::to_string(slack) + " bases (BWA-MEM's scoring, X = 40)");
  }
  {
    // The Z-drop (expensive gaps, so the drop is driven by mismatches): A (300) · divergent block
    // (100, unrelated in the two) · B (800). The X-drop alone crosses the block, since B repays it;
    // with Z = 100 the extension stops at the end of A.
    const Scoring z = scoring(Matrix::dna(2, -4, -1), 40, 2);
    const Codes A = simulate::random_sequence(300, 0.5, 31), B = simulate::random_sequence(800, 0.5, 32);
    const Codes R1 = simulate::random_sequence(100, 0.5, 33), R2 = simulate::random_sequence(100, 0.5, 34);
    Codes q = A, t = A;
    q.insert(q.end(), R1.begin(), R1.end());
    t.insert(t.end(), R2.begin(), R2.end());
    q.insert(q.end(), B.begin(), B.end());
    t.insert(t.end(), B.begin(), B.end());
    const Alignment xonly = extend_right(q, t, 0, 0, z, {INF, -1}), zdrop = extend_right(q, t, 0, 0, z, {INF, 100}), zhuge = extend_right(q, t, 0, 0, z, {INF, INF});
    s.require("Z-drop: without it the extension crosses the divergent block to the end of B", xonly.query_end == q.size() && xonly.target_end == t.size(),
              "ends at " + std::to_string(xonly.query_end) + " of " + std::to_string(q.size()));
    s.require("Z-drop: with Z = 100 the extension stops at the end of A, before the block", zdrop.query_end <= 305 && zdrop.query_end >= 295, "ends at " + std::to_string(zdrop.query_end));
    s.require("Z-drop: an infinite Z changes nothing", zhuge.score == xonly.score && zhuge.cigar == xonly.cigar);
  }
  return s.done();
}
