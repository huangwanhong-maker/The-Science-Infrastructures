// Checks of banded.hpp: a band holding the whole edit graph is Gotoh; a band holding the diagonals
// of one of Gotoh's optimal paths gives Gotoh's score, in every mode; a narrower band never scores
// more than Gotoh, and its path stays inside the band and rescores to its claim; the certificate
// banded_is_optimal is never wrong when it says yes; banded_global (doubling bands) equals Gotoh's
// global score on homologous pairs with long indels; bands that hold no path are refused.

#include <genomics/align/banded.hpp>
#include <genomics/align/gotoh.hpp>
#include <genomics/align/scoring.hpp>
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

Codes random_codes(calx::random::Stream& r, std::size_t n, std::size_t sigma) {
  Codes x(n);
  for (auto& c : x) c = static_cast<Code>(r.below(sigma));
  return x;
}

// The least and greatest diagonal j − i a path visits.
std::pair<std::int64_t, std::int64_t> diagonals(const Alignment& a) {
  std::int64_t i = static_cast<std::int64_t>(a.query_begin), j = static_cast<std::int64_t>(a.target_begin);
  std::int64_t lo = j - i, hi = j - i;
  for (const auto& o : a.cigar)
    for (std::uint32_t k = 0; k < o.length; k++) {
      if (consumes_query(o.op)) i++;
      if (consumes_target(o.op)) j++;
      lo = std::min(lo, j - i);
      hi = std::max(hi, j - i);
    }
  return {lo, hi};
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x62616e64, 0, 0);
  const Scoring sc = scoring(Matrix::dna(2, -3, -1), 5, 2);
  for (Mode mode : {Mode::global, Mode::local, Mode::semiglobal, Mode::overlap}) {
    bool whole = true, holds = true, never_more = true, inside = true;
    for (int rep = 0; rep < 300; rep++) {
      const std::size_t m = 1 + r.below(50), n = 1 + r.below(50);
      const Codes x = random_codes(r, m, 4), y = random_codes(r, n, 4);
      const Alignment g = gotoh(x, y, sc, mode);
      const auto full = banded(x, y, sc, mode, -static_cast<std::int64_t>(m), static_cast<std::int64_t>(n));
      whole = whole && full.score == g.score && full.cigar == g.cigar;
      const auto [lo, hi] = diagonals(g);
      holds = holds && banded(x, y, sc, mode, lo, hi).score == g.score;
      // A random band containing the corners a global path needs.
      const std::int64_t d = static_cast<std::int64_t>(n) - static_cast<std::int64_t>(m);
      const std::int64_t blo = std::min<std::int64_t>(0, d) - static_cast<std::int64_t>(r.below(4)), bhi = std::max<std::int64_t>(0, d) + static_cast<std::int64_t>(r.below(4));
      const Alignment b = banded(x, y, sc, mode, blo, bhi);
      const auto [plo, phi] = diagonals(b);
      never_more = never_more && b.score <= g.score;
      inside = inside && plo >= blo && phi <= bhi && rescore(b, x, y, sc) == b.score;
    }
    const std::string tag = std::string("banded, ") + mode_name(mode) + ": ";
    s.require(tag + "the band of all diagonals is Gotoh, path for path", whole, "300 random pairs of length ≤ 50");
    s.require(tag + "the band of an optimal path's diagonals reaches the optimum", holds);
    s.require(tag + "a narrower band never scores more than the whole graph", never_more);
    s.require(tag + "its path stays in the band and rescores to its claim", inside);
  }
  {
    bool sound = true;
    int certified = 0;
    for (int rep = 0; rep < 400; rep++) {
      const Codes root = random_codes(r, 30 + r.below(60), 4);
      const auto f = simulate::evolve_family(root, 2, 0.1, 0.03, 1000 + static_cast<std::uint64_t>(rep));
      const Codes& x = f.sequences[0];
      const Codes& y = f.sequences[1];
      const std::int64_t d = static_cast<std::int64_t>(y.size()) - static_cast<std::int64_t>(x.size());
      const std::int64_t w = static_cast<std::int64_t>(r.below(12));
      const std::int64_t lo = std::min<std::int64_t>(0, d) - w, hi = std::max<std::int64_t>(0, d) + w;
      const Alignment b = banded(x, y, sc, Mode::global, lo, hi);
      if (banded_is_optimal(x.size(), y.size(), sc, lo, hi, b.score)) {
        certified++;
        sound = sound && b.score == gotoh_score(x, y, sc, Mode::global);
      }
    }
    s.require("certificate: whenever banded_is_optimal says yes, the banded score is Gotoh's", sound && certified > 0, std::to_string(certified) + " of 400 bands certified");
  }
  {
    bool ok = true;
    for (std::uint64_t seed = 1; seed <= 10; seed++) {
      const Codes x = simulate::random_sequence(600, 0.5, seed);
      Codes y = x;
      // A deletion of 40 + seed bases and an insertion of 25: the optimum leaves narrow bands.
      y.erase(y.begin() + 100, y.begin() + 140 + static_cast<std::ptrdiff_t>(seed));
      const Codes ins = simulate::random_sequence(25, 0.5, seed + 100);
      y.insert(y.begin() + 400, ins.begin(), ins.end());
      ok = ok && banded_global(x, y, sc, 4).score == gotoh_score(x, y, sc, Mode::global);
    }
    s.require("banded_global: doubling bands from half-width 4 reach Gotoh's optimum across indels of 41–50 bases", ok);
  }
  {
    const Codes x = random_codes(r, 20, 4), y = random_codes(r, 30, 4);
    bool refused = false;
    try {
      banded(x, y, sc, Mode::global, -5, 5);  // n − m = 10 lies outside
    } catch (const std::invalid_argument&) { refused = true; }
    s.require("banded: a global band that misses the corner (m, n) is refused", refused);
  }
  return s.done();
}
