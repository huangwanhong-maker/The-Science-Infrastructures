// Checks of the striped Smith–Waterman kernel: its best local score and its end (the first
// cell of that score in row-major order) equal Gotoh's recurrences (gotoh_end in local mode) on
// DNA and protein pairs of every length across the stripe boundaries, under the scorings of
// BLAST, BWA-MEM and minimap2 and a linear one; scores past 255 drive the 16-bit pass, past
// 32767 the 32-bit one, with the same answers; the full alignment from the reversed prefixes and
// the banded region holds its score (align::rescore) and lies inside the sequences; the lazy-F
// loop is exercised by long gaps.

#include <genomics/align/gotoh.hpp>
#include <genomics/align/scoring.hpp>
#include <genomics/align/striped.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cstdint>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::align;

namespace {

constexpr std::uint32_t NOISE = 0x53545250;

// A copy of x with substitutions and indels of 1 to `longest` bases, at `rate` per base.
Codes noisy(const Codes& x, double rate, std::uint64_t seed, std::uint32_t object, int alphabet = 4, std::size_t longest = 3) {
  calx::random::Stream r(seed, NOISE, object, 0);
  Codes y;
  for (std::size_t i = 0; i < x.size(); i++) {
    const double u = r.uniform();
    if (u < rate / 3) { y.push_back(static_cast<Code>((x[i] + 1 + r.below(static_cast<std::uint64_t>(alphabet - 1))) % static_cast<std::uint64_t>(alphabet))); continue; }
    if (u < 2 * rate / 3) { i += r.below(longest); continue; }
    if (u < rate)
      for (std::size_t k = 0, len = 1 + r.below(longest); k < len; k++) y.push_back(static_cast<Code>(r.below(static_cast<std::uint64_t>(alphabet))));
    y.push_back(x[i]);
  }
  return y;
}

Codes random_codes(std::size_t n, int alphabet, std::uint64_t seed, std::uint32_t object) {
  calx::random::Stream r(seed, NOISE + 1, object, 0);
  Codes x(n);
  for (auto& c : x) c = static_cast<Code>(r.below(static_cast<std::uint64_t>(alphabet)));
  return x;
}

// A target holding a diverged copy of a piece of the query, between random flanks.
Codes related(const Codes& q, double rate, int alphabet, std::uint64_t seed, std::uint32_t object) {
  calx::random::Stream r(seed, NOISE + 2, object, 0);
  const std::size_t a = q.empty() ? 0 : r.below(q.size()), b = a + (q.size() - a) / 2 + 1;
  Codes t = random_codes(r.below(200) + 1, alphabet, seed + 1, object);
  const Codes piece(q.begin() + static_cast<std::ptrdiff_t>(a), q.begin() + static_cast<std::ptrdiff_t>(std::min(b, q.size())));
  const Codes copy = noisy(piece, rate, seed + 2, object, alphabet);
  t.insert(t.end(), copy.begin(), copy.end());
  const Codes tail = random_codes(r.below(200) + 1, alphabet, seed + 3, object);
  t.insert(t.end(), tail.begin(), tail.end());
  return t;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  s.skip(std::string("lanes: ") + calx::cpu::ISA + ", " + std::to_string(calx::cpu::LANE_BYTES) + " bytes (" + std::to_string(calx::cpu::LANE_BYTES) + " lanes of 8 bits)");

  struct Named {
    std::string name;
    Scoring scoring;
    int alphabet;
  };
  Scoring linear;
  linear.matrix = Matrix::dna(1, -1, -1);
  linear.gap_open = 0;
  linear.gap_extend = 1;
  const std::vector<Named> scorings = {{"minimap2 map-ont", minimap2_map_ont(), 4}, {"blastn", blastn(), 5},   {"bwa-mem", bwa_mem(), 4},
                                       {"linear (1, -1, 0, 1)", linear, 4},         {"blastp", blastp(), 20}, {"BLOSUM80, 10, 1", scoring(matrix("BLOSUM80"), 10, 1), 20}};
  const std::vector<std::size_t> lengths = {1, 2, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 127, 128, 129, 300, 1000};

  // ---- scores and ends against Gotoh's recurrences ---------------------------------------------------
  int ok = 0, total = 0, full_ok = 0;
  std::string bad;
  std::uint32_t object = 0;
  for (const auto& sc : scorings)
    for (std::size_t m : lengths)
      for (int kind = 0; kind < 3; kind++) {
        object++;
        const Codes q = random_codes(m, sc.alphabet, 1, object);
        const Codes t = kind == 0 ? random_codes(m + 40, sc.alphabet, 2, object) : related(q, kind == 1 ? 0.05 : 0.25, sc.alphabet, 3, object);
        const Endpoint want = gotoh_end(q, t, sc.scoring, Mode::local);
        StripedStats st;
        const Endpoint got = striped_end(StripedProfile(q, sc.scoring), t, &st);
        total++;
        if (got.score == want.score && (want.score == 0 || (got.query_end == want.query_end && got.target_end == want.target_end))) ok++;
        else if (bad.empty())
          bad = sc.name + " m=" + std::to_string(m) + ": " + std::to_string(got.score) + " at (" + std::to_string(got.query_end) + ", " + std::to_string(got.target_end) +
                ") vs " + std::to_string(want.score) + " at (" + std::to_string(want.query_end) + ", " + std::to_string(want.target_end) + ")";
        const Alignment a = striped_local(q, t, sc.scoring);
        if (a.score == want.score &&
            (a.score == 0 || (a.query_end == want.query_end && a.target_end == want.target_end && rescore(a, q, t, sc.scoring) == a.score && a.query_begin < a.query_end &&
                              a.target_begin < a.target_end)))
          full_ok++;
        else if (bad.empty()) bad = "alignment " + sc.name + " m=" + std::to_string(m);
      }
  s.require("the striped score and end equal Gotoh's local optimum and its first cell on " + std::to_string(total) + " DNA and protein pairs under 6 scorings",
            ok == total, bad);
  s.require("every full alignment, its start from the reversed prefixes and its path from the banded region, holds its score", full_ok == total, bad);

  // ---- the widths: 8 bits, then 16, then 32 ---------------------------------------------------------------
  {
    const Scoring sc = minimap2_map_ont();  // match 2
    const Codes g = simulate::random_sequence(40'000, 0.5, 11);
    struct Case { std::size_t length; int bits; };
    for (const Case c : {Case{100, 8}, Case{1'000, 16}, Case{20'000, 32}}) {
      const Codes q(g.begin(), g.begin() + static_cast<std::ptrdiff_t>(c.length));
      Codes t = simulate::random_sequence(500, 0.5, 12);
      const Codes copy = noisy(q, 0.01, 13, static_cast<std::uint32_t>(c.length));
      t.insert(t.end(), copy.begin(), copy.end());
      StripedStats st;
      const Endpoint got = striped_end(StripedProfile(q, sc), t, &st);
      const Endpoint want = gotoh_end(q, t, sc, Mode::local);
      s.require("a local score of " + std::to_string(want.score) + " is found in " + std::to_string(c.bits) + "-bit lanes (" + std::to_string(st.passes) + " passes), as Gotoh's",
                got.score == want.score && got.query_end == want.query_end && got.target_end == want.target_end && st.bits == c.bits,
                std::to_string(got.score) + " in " + std::to_string(st.bits) + " bits");
    }
  }

  // ---- long gaps exercise the lazy-F loop -----------------------------------------------------------------
  {
    Scoring sc = minimap2_map_ont();
    sc.gap_open = 2;
    sc.gap_extend = 1;
    int lazy_ok = 0, lazy_total = 0;
    std::uint64_t lazy = 0;
    for (std::uint32_t i = 0; i < 20; i++) {
      // A query with a long insertion relative to the target: a vertical gap that crosses lanes.
      const Codes left = simulate::random_sequence(300, 0.5, 21, i), mid = simulate::random_sequence(150 + 10 * i, 0.5, 22, i), right = simulate::random_sequence(300, 0.5, 23, i);
      Codes q = left;
      q.insert(q.end(), mid.begin(), mid.end());
      q.insert(q.end(), right.begin(), right.end());
      Codes t = left;
      t.insert(t.end(), right.begin(), right.end());
      StripedStats st;
      const Alignment a = striped_local(StripedProfile(q, sc), t, &st);
      lazy += st.lazy_f;
      lazy_total++;
      if (a.score == gotoh_score(q, t, sc, Mode::local) && rescore(a, q, t, sc) == a.score) lazy_ok++;
    }
    s.require("queries with insertions of 150 to 340 bases align at Gotoh's optimum, through " + std::to_string(lazy) + " lazy-F corrections", lazy_ok == lazy_total && lazy > 0);
  }

  return s.done();
}
