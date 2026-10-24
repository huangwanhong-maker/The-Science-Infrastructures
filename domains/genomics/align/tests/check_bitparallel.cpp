// Checks of Myers' bit-vector algorithm: the published distances of the textbook examples; the
// global distance and every semiglobal end with its distance equal to the O(mn) recurrence on
// random and on noisy pairs of lengths 1 to 5000, the queries on both sides of every 64-bit
// block boundary; the band never changes an answer within its threshold; the threshold refuses
// what lies beyond it; every CIGAR consumes both sequences as its mode says and costs what it
// claims (align::edit_distance); N matches nothing.

#include <genomics/align/bitparallel.hpp>
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

constexpr std::uint32_t NOISE = 0x4e4f4953;  // the purpose of the noise streams

// C(m, j) for j = 0 .. n, by the O(mn) recurrence of the edit distance (Wagner and Fischer
// 1974; Sellers 1980 for the free start in the target), codes below sigma matching themselves.
std::vector<std::int64_t> last_row(const Codes& q, const Codes& t, bool global, int sigma = 4) {
  const std::size_t m = q.size(), n = t.size();
  std::vector<std::int64_t> col(m + 1), prev(m + 1), out(n + 1);
  for (std::size_t i = 0; i <= m; i++) col[i] = static_cast<std::int64_t>(i);
  out[0] = static_cast<std::int64_t>(m);
  for (std::size_t j = 1; j <= n; j++) {
    std::swap(col, prev);
    col[0] = global ? static_cast<std::int64_t>(j) : 0;
    for (std::size_t i = 1; i <= m; i++) {
      const bool same = q[i - 1] == t[j - 1] && q[i - 1] < sigma;
      col[i] = std::min({prev[i - 1] + (same ? 0 : 1), prev[i] + 1, col[i - 1] + 1});
    }
    out[j] = col[m];
  }
  return out;
}

// A copy of x with substitutions, insertions and deletions, each at rate/3 per base.
Codes noisy(const Codes& x, double rate, std::uint64_t seed, std::uint32_t object) {
  calx::random::Stream r(seed, NOISE, object, 0);
  Codes y;
  y.reserve(x.size() + x.size() / 8 + 4);
  for (std::size_t i = 0; i < x.size(); i++) {
    const double u = r.uniform();
    if (u < rate / 3) { y.push_back(simulate::substitute(r, x[i])); continue; }
    if (u < 2 * rate / 3) continue;
    if (u < rate) y.push_back(static_cast<Code>(r.below(4)));
    y.push_back(x[i]);
  }
  return y;
}

Codes letters(const std::string& s) {
  Codes c;
  for (char ch : s) c.push_back(static_cast<Code>(ch - 'a'));
  return c;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- textbook examples -------------------------------------------------------------------------
  {
    // Levenshtein's kitten → sitting: three edits (two substitutions and an insertion).
    const auto d1 = myers_distance(letters("kitten"), letters("sitting"), Mode::global, -1, 26);
    // Gusfield 1997, Algorithms on Strings, Trees and Sequences, §11.2: vintner becomes writers
    // by the transcript RIMDMDMMI, of cost 5, the edit distance.
    const auto d2 = myers_distance(letters("vintner"), letters("writers"), Mode::global, -1, 26);
    s.require("the edit distance of kitten and sitting is 3, and of vintner and writers 5 (Gusfield 1997, §11.2)", d1 == 3 && d2 == 5,
              std::to_string(d1) + ", " + std::to_string(d2));
    // Approximate matching, the running example of Navarro (2001, ACM Comput. Surv. 33, 31–88):
    // "survey" in "surgery" with at most 2 edits; "surgery"'s prefix "surge" at 2 edits.
    const auto hits = myers_search(letters("survey"), letters("surgery"), 2, 26);
    const auto ref = last_row(letters("survey"), letters("surgery"), false, 26);
    std::vector<EditHit> want;
    for (std::size_t j = 0; j < ref.size(); j++)
      if (ref[j] <= 2) want.push_back({j, ref[j]});
    s.require("the matches of survey in surgery within 2 edits are the recurrence's, among them surge at 2", hits == want && !hits.empty() &&
              std::find(hits.begin(), hits.end(), EditHit{5, 2}) != hits.end(), std::to_string(hits.size()) + " hits");
    const Codes n4(4, sequence::N);
    s.require("N matches nothing: NNNN against NNNN is at distance 4", myers_distance(n4, n4) == 4 && myers_distance(sequence::encode_dna("ACGN"), sequence::encode_dna("ACGT")) == 1);
  }

  // ---- distances against the recurrence ---------------------------------------------------------
  const std::vector<std::size_t> lengths = {1, 2, 5, 31, 63, 64, 65, 100, 127, 128, 129, 191, 192, 193, 255, 256, 257, 500, 1000, 2000, 5000};
  const std::vector<double> rates = {0.0, 0.01, 0.05, 0.2, 1.0};  // 1.0: an unrelated target
  int global_ok = 0, global_total = 0, semi_ok = 0, semi_total = 0, cigar_ok = 0, cigar_total = 0, limit_ok = 0, limit_total = 0;
  std::string first_bad;
  std::uint32_t object = 0;
  for (std::size_t m : lengths)
    for (double rate : rates) {
      object++;
      const Codes q = simulate::random_sequence(m, 0.5, 7, object);
      Codes t = rate >= 1.0 ? simulate::random_sequence(m + m / 7 + 3, 0.5, 8, object) : noisy(q, rate, 9, object);
      if (t.empty()) t.push_back(sequence::A);
      // global
      const auto rg = last_row(q, t, true);
      const std::int64_t dg = myers_distance(q, t, Mode::global);
      global_total++;
      if (dg == rg.back()) global_ok++;
      else if (first_bad.empty()) first_bad = "global m=" + std::to_string(m) + ": " + std::to_string(dg) + " vs " + std::to_string(rg.back());
      const auto ag = myers_align(q, t, Mode::global);
      cigar_total++;
      if (ag && query_length(ag->cigar) == q.size() && target_length(ag->cigar) == t.size() && ag->target_begin == 0 &&
          static_cast<std::int64_t>(edit_distance(*ag, q, t)) == dg && ag->score == -dg)
        cigar_ok++;
      else if (first_bad.empty()) first_bad = "global CIGAR m=" + std::to_string(m);
      // the threshold: refused just below the distance, found at it
      limit_total++;
      if (myers_distance(q, t, Mode::global, dg) == dg && (dg == 0 || myers_distance(q, t, Mode::global, dg - 1) == -1) && (dg == 0 || !myers_align(q, t, Mode::global, dg - 1)))
        limit_ok++;
      // semiglobal: every end within k, with its distance; the query inside a longer target
      Codes text = simulate::random_sequence(m / 2 + 17, 0.5, 10, object);
      text.insert(text.end(), t.begin(), t.end());
      const Codes tail = simulate::random_sequence(m / 3 + 11, 0.5, 11, object);
      text.insert(text.end(), tail.begin(), tail.end());
      const auto rs = last_row(q, text, false);
      const std::int64_t best = *std::min_element(rs.begin(), rs.end());
      for (std::int64_t k : {best, best + 3, best + 40, static_cast<std::int64_t>(m) / 2}) {
        if (k < 0) continue;
        std::vector<EditHit> want;
        for (std::size_t j = 0; j < rs.size(); j++)
          if (rs[j] <= k) want.push_back({j, rs[j]});
        semi_total++;
        if (myers_search(q, text, k) == want) semi_ok++;
        else if (first_bad.empty()) first_bad = "semiglobal m=" + std::to_string(m) + " k=" + std::to_string(k);
      }
      const auto as = myers_align(q, text, Mode::semiglobal);
      const std::size_t leftmost = static_cast<std::size_t>(std::find(rs.begin(), rs.end(), best) - rs.begin());
      cigar_total++;
      if (as && myers_distance(q, text, Mode::semiglobal) == best && as->score == -best && as->target_end == leftmost && query_length(as->cigar) == q.size() &&
          as->target_begin + target_length(as->cigar) == as->target_end && static_cast<std::int64_t>(edit_distance(*as, q, text)) == best)
        cigar_ok++;
      else if (first_bad.empty()) first_bad = "semiglobal CIGAR m=" + std::to_string(m);
    }
  s.require("the global distance equals the O(mn) recurrence on " + std::to_string(global_total) + " pairs of lengths 1 to 5000, across every block boundary",
            global_ok == global_total, first_bad);
  s.require("every semiglobal end within k, and its distance, equals the recurrence's (" + std::to_string(semi_total) + " searches)", semi_ok == semi_total, first_bad);
  s.require("every CIGAR spans its mode's sequences and has the edit distance it claims; the semiglobal end is the leftmost best (" + std::to_string(cigar_total) + ")",
            cigar_ok == cigar_total, first_bad);
  s.require("a threshold below the distance refuses the pair, and at the distance finds it", limit_ok == limit_total);

  // ---- many small pairs over two letters, where ties abound, with thresholds of every size ------------
  {
    calx::random::Stream r(31, NOISE, 1, 0);
    int ok = 0, total = 0;
    std::string bad;
    for (int trial = 0; trial < 3000; trial++) {
      const std::size_t m = 1 + r.below(trial < 1500 ? 70 : 300), n = 1 + r.below(trial < 1500 ? 70 : 300);
      Codes q(m), t(n);
      for (auto& c : q) c = static_cast<Code>(r.below(2));
      for (auto& c : t) c = static_cast<Code>(r.below(2));
      if (trial % 7 == 0) q[r.below(m)] = sequence::N;
      const std::int64_t k = static_cast<std::int64_t>(r.below(static_cast<std::uint64_t>(std::max(m, n)) + 2));
      const auto rg = last_row(q, t, true), rs = last_row(q, t, false);
      std::vector<EditHit> want;
      for (std::size_t j = 0; j < rs.size(); j++)
        if (rs[j] <= k) want.push_back({j, rs[j]});
      const std::int64_t dg = rg.back() <= k ? rg.back() : -1;
      const auto a = myers_align(q, t, Mode::global, k);
      total++;
      if (myers_search(q, t, k) == want && myers_distance(q, t, Mode::global, k) == dg && (dg < 0 ? !a : (a && static_cast<std::int64_t>(edit_distance(*a, q, t)) == dg)))
        ok++;
      else if (bad.empty()) bad = "m=" + std::to_string(m) + " n=" + std::to_string(n) + " k=" + std::to_string(k);
    }
    s.require("on 3000 small pairs over two letters, every threshold gives the recurrence's ends, distances and paths", ok == total, bad);
  }

  // ---- the band changes nothing within its threshold ---------------------------------------------
  {
    const Codes q = simulate::random_sequence(30'000, 0.45, 21);
    const Codes t = noisy(q, 0.02, 22, 0);
    const std::int64_t banded = myers_distance(q, t);
    // The whole matrix: a threshold no cell exceeds, so the band is every block.
    std::int64_t whole = -1;
    detail::myers_pass(MyersPattern(q), t.data(), t.size(), true, static_cast<std::int64_t>(q.size() + t.size()), nullptr, [&](std::size_t j, std::int64_t d) {
      if (j == t.size()) whole = d;
    });
    MyersStats work;
    const auto a = myers_align(q, t, Mode::global, -1, 4, &work);
    s.require("on 30 kb at 2% divergence, the banded distance equals the unbanded one, and its CIGAR has it", banded == whole && a && static_cast<std::int64_t>(edit_distance(*a, q, t)) == banded,
              std::to_string(banded) + " vs " + std::to_string(whole));
    const double share = static_cast<double>(work.cells()) / (static_cast<double>(q.size() + 1) * static_cast<double>(t.size() + 1));
    s.require("and its passes, every threshold of the doubling and the traceback's, computed under a fifth of the matrix", share < 0.2,
              std::to_string(work.passes) + " passes, " + std::to_string(share) + " of the cells");
    // A read placed in a genome: the best end is where it came from.
    const Codes genome = simulate::random_sequence(200'000, 0.5, 23);
    const Codes read(genome.begin() + 123'456, genome.begin() + 123'456 + 1500);
    const Codes noisy_read = noisy(read, 0.05, 24, 0);
    const auto hit = myers_align(noisy_read, genome, Mode::semiglobal);
    s.require("a 1.5 kb read with 5% noise is found where it came from in a 200 kb genome",
              hit && hit->target_end >= 123'456 + 1490 && hit->target_end <= 123'456 + 1510 && hit->target_begin >= 123'446 && hit->target_begin <= 123'466 &&
                  static_cast<std::int64_t>(edit_distance(*hit, noisy_read, genome)) == -hit->score,
              hit ? std::to_string(hit->target_begin) + ".." + std::to_string(hit->target_end) : "none");
  }

  return s.done();
}
