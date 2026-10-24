// Checks of hirschberg.hpp: Myers and Miller's linear-space alignment with affine gaps gives
// Gotoh's optimal score, with a path that rescores to it and spans what its mode allows, on random
// DNA and protein pairs of every length up to 60 (one-row, empty and unbalanced sub-problems all
// reached), with linear gaps, with large opens that make vertical gaps cross the middle row, and
// with free extension; in all four modes; and on a homologous pair of 4 kb.

#include <genomics/align/gotoh.hpp>
#include <genomics/align/hirschberg.hpp>
#include <genomics/align/scoring.hpp>
#include <genomics/simulate/family.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::align;

namespace {

Codes random_codes(calx::random::Stream& r, std::size_t n, std::size_t sigma) {
  Codes x(n);
  for (auto& c : x) c = static_cast<Code>(r.below(sigma));
  return x;
}

bool spans(const Alignment& a, std::size_t m, std::size_t n, Mode mode) {
  if (query_length(a.cigar) != a.query_end - a.query_begin || target_length(a.cigar) != a.target_end - a.target_begin) return false;
  switch (mode) {
    case Mode::global: return a.query_begin == 0 && a.target_begin == 0 && a.query_end == m && a.target_end == n;
    case Mode::local: return a.query_end <= m && a.target_end <= n;
    case Mode::semiglobal: return a.query_begin == 0 && a.query_end == m && a.target_end <= n;
    case Mode::overlap: return (a.query_begin == 0 || a.target_begin == 0) && (a.query_end == m || a.target_end == n);
  }
  return false;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x6869726b, 0, 0);
  const std::vector<std::pair<std::string, Scoring>> scorings = {
      {"affine 2/−3/5/2", scoring(Matrix::dna(2, -3, -1), 5, 2)}, {"linear 1/−1/0/2", scoring(Matrix::dna(1, -1, -1), 0, 2)},
      {"open 20, extend 1", scoring(Matrix::dna(3, -2, -1), 20, 1)}, {"free extension 2/−3/6/0", scoring(Matrix::dna(2, -3, -1), 6, 0)}, {"BLOSUM62 11/1", blastp()}};
  for (Mode mode : {Mode::global, Mode::local, Mode::semiglobal, Mode::overlap}) {
    bool score = true, rescored = true, span = true;
    std::size_t pairs = 0;
    for (const auto& [name, sc] : scorings) {
      const std::size_t sigma = sc.matrix.size == 24 ? 20 : 4;
      for (int rep = 0; rep < 160; rep++) {
        const std::size_t m = r.below(61), n = rep % 5 == 0 ? r.below(4) : r.below(61);
        const Codes x = random_codes(r, m, sigma), y = random_codes(r, n, sigma);
        const Alignment h = hirschberg(x, y, sc, mode);
        score = score && h.score == gotoh_score(x, y, sc, mode);
        rescored = rescored && rescore(h, x, y, sc) == h.score;
        span = span && spans(h, m, n, mode);
        pairs++;
      }
    }
    const std::string tag = std::string("hirschberg, ") + mode_name(mode) + ": ";
    s.require(tag + "Myers–Miller's score is Gotoh's", score, std::to_string(pairs) + " random pairs of length ≤ 60 under five scorings");
    s.require(tag + "its path rescores to that score", rescored);
    s.require(tag + "its path spans what the mode allows", span);
  }
  {
    // Repeats force long vertical gaps through the middle row: a sequence against itself with a
    // block deleted.
    bool ok = true;
    for (std::uint64_t seed = 1; seed <= 40; seed++) {
      const Codes x = simulate::random_sequence(80, 0.5, seed);
      Codes y(x.begin(), x.begin() + 20);
      y.insert(y.end(), x.begin() + 20 + static_cast<std::ptrdiff_t>(seed % 37), x.end());
      for (const auto& [name, sc] : scorings) {
        if (sc.matrix.size != 5) continue;
        const Alignment h = hirschberg(x, y, sc, Mode::global);
        ok = ok && h.score == gotoh_score(x, y, sc, Mode::global) && rescore(h, x, y, sc) == h.score;
      }
    }
    s.require("hirschberg: a block deleted from a sequence (a gap crossing the middle row) is aligned optimally", ok);
  }
  {
    const Codes root = simulate::random_sequence(4000, 0.45, 9);
    const auto f = simulate::evolve_family(root, 2, 0.08, 0.02, 9);
    bool ok = true;
    for (Mode mode : {Mode::global, Mode::local}) {
      const Alignment h = hirschberg(f.sequences[0], f.sequences[1], minimap2_map_ont(), mode);
      ok = ok && h.score == gotoh_score(f.sequences[0], f.sequences[1], minimap2_map_ont(), mode) && rescore(h, f.sequences[0], f.sequences[1], minimap2_map_ont()) == h.score;
    }
    s.require("hirschberg: a homologous pair of 4 kb (minimap2's map-ont scoring), global and local, equals Gotoh", ok);
  }
  return s.done();
}
