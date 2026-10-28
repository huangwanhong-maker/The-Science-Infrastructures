// Checks of maximal unique matches: on 400 small random pairs of genomes sharing segments, some
// reverse-complemented, some repeated (so that they are not unique), with bases N, the MUMs of
// both strands equal those found by brute force (every pair of positions extended to its maximal
// exact match, kept when left-maximal, long enough, and occurring once in each genome); the
// matches are exact on their strand; no MUM nests inside another, in A or in B (theory,
// Proposition 4.4); and a genome against itself has one MUM, the whole genome, when it holds no
// repeat.

#include <genomics/comparative/mums.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <set>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::comparative;
using sequence::Codes;

namespace {

Codes random_dna(calx::random::Stream& r, std::size_t n) {
  Codes x(n);
  for (auto& c : x) c = static_cast<sequence::Code>(r.below(4));
  return x;
}

std::size_t occurrences(const Codes& x, const Codes& w) {
  std::size_t n = 0;
  for (std::size_t i = 0; i + w.size() <= x.size(); i++)
    if (std::equal(w.begin(), w.end(), x.begin() + static_cast<std::ptrdiff_t>(i))) n++;
  return n;
}

// MUMs by brute force: maximal exact matches between A and the strand Bs, left-maximal by
// construction, of at least m bases, whose word occurs once in A and once in Bs.
std::vector<Match> brute(const Codes& A, const Codes& B, std::size_t m, bool reverse) {
  const Codes Bs = reverse ? sequence::reverse_complement(B) : B;
  std::vector<Match> out;
  for (std::size_t p = 0; p < A.size(); p++)
    for (std::size_t q = 0; q < Bs.size(); q++) {
      if (A[p] > 3 || A[p] != Bs[q]) continue;
      if (p > 0 && q > 0 && A[p - 1] <= 3 && A[p - 1] == Bs[q - 1]) continue;
      std::size_t l = 0;
      while (p + l < A.size() && q + l < Bs.size() && A[p + l] <= 3 && A[p + l] == Bs[q + l]) l++;
      if (l < m) continue;
      const Codes w(A.begin() + static_cast<std::ptrdiff_t>(p), A.begin() + static_cast<std::ptrdiff_t>(p + l));
      if (occurrences(A, w) != 1 || occurrences(Bs, w) != 1) continue;
      out.push_back({p, reverse ? B.size() - q - l : q, l, reverse});
    }
  std::sort(out.begin(), out.end());
  return out;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(1999, 0x4d554d53, 0, 0);

  bool equal = true, exact = true, nested = false;
  std::size_t total = 0;
  for (int trial = 0; trial < 400; trial++) {
    // Shared segments, placed in both genomes, some on the other strand, some twice in one genome.
    std::vector<Codes> segments;
    for (int i = 0; i < 5; i++) segments.push_back(random_dna(r, 5 + r.below(30)));
    const auto build = [&] {
      Codes g = random_dna(r, r.below(30));
      for (int j = 0; j < 6; j++) {
        Codes seg = segments[r.below(segments.size())];
        if (r.below(3) == 0) seg = sequence::reverse_complement(seg);
        if (r.below(4) == 0 && !seg.empty()) seg[r.below(seg.size())] = static_cast<sequence::Code>(r.below(4));  // a substitution
        g.insert(g.end(), seg.begin(), seg.end());
        const Codes gap = random_dna(r, r.below(25));
        g.insert(g.end(), gap.begin(), gap.end());
      }
      for (auto& c : g)
        if (r.below(80) == 0) c = sequence::N;
      return g;
    };
    const Codes A = build(), B = build();
    const std::size_t m = 3 + r.below(10);
    for (bool reverse : {false, true}) {
      const auto got = maximal_unique_matches(A, B, m, reverse);
      equal = equal && got == brute(A, B, m, reverse);
      total += got.size();
      const Codes rc = sequence::reverse_complement(B);
      for (const auto& x : got) {
        for (std::size_t t = 0; t < x.length; t++) {
          const sequence::Code b = reverse ? sequence::complement(B[x.b + x.length - 1 - t]) : B[x.b + t];
          exact = exact && A[x.a + t] <= 3 && A[x.a + t] == b;
        }
        for (const auto& y : got)
          if (!(x == y) && ((y.a <= x.a && x.a_end() <= y.a_end()) || (y.b <= x.b && x.b_end() <= y.b_end()))) nested = true;
      }
    }
  }
  s.require("on 400 random pairs sharing segments on both strands, with repeats and bases N, the MUMs of each strand equal those found by brute force (" + std::to_string(total) + " in all)", equal);
  s.require("every MUM is an exact match on its strand, free of N", exact);
  s.require("no MUM's interval contains another's, in A or in B (Proposition 4.4)", !nested);

  // A genome against itself.
  const Codes x = random_dna(r, 5000);
  const auto self = maximal_unique_matches(x, x, 20);
  s.require("a random genome of 5 kb against itself has one forward MUM, itself", self.size() == 1 && self[0].a == 0 && self[0].b == 0 && self[0].length == x.size());
  s.require("and no reverse MUM of 20 bases or more (none expected by chance: 25 × 10⁶ / 4²⁰ ≈ 2 × 10⁻⁵)", maximal_unique_matches(x, x, 20, true).empty());
  return s.done();
}
