// Checks of k-mers and their samples: k-mers from 1 to 32 encoded, decoded and
// reverse-complemented against the codes; every k-mer free of N visited once with its
// reverse complement; hash64 a bijection of [0, 4^k); ntHash rolled against ntHash computed
// directly, across runs of N, for k from 1 to 100, and its canonical value the same on the
// reverse strand; minimizers against their definition by brute force (every window's
// minimal positions, ties kept, palindromes skipped), the same on the reverse complement, and
// their density 2/(w + 1) on a random megabase; closed and open syncmers against their
// definition and their densities 2/(k − s + 1) and 1/(k − s + 1).

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/kmer.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace genomics::sequence;

namespace {

Codes random_dna(calx::random::Stream& r, std::size_t n, double n_rate = 0) {
  Codes x(n);
  for (auto& c : x) c = r.uniform() < n_rate ? N : static_cast<Code>(r.below(4));
  return x;
}

bool free_of_n(const Codes& x, std::size_t i, std::size_t k) {
  for (std::size_t j = i; j < i + k; j++)
    if (x[j] > 3) return false;
  return true;
}

// Minimizers by their definition: every window of w consecutive k-mers within a run free of N
// (a shorter run being one window), its minimal hashes, ties kept, palindromic k-mers skipped.
std::set<std::pair<std::uint32_t, std::uint64_t>> brute_minimizers(const Codes& x, int k, int w) {
  std::set<std::pair<std::uint32_t, std::uint64_t>> out;
  const Kmer mask = kmer_mask(k);
  // The k-mers of each run.
  std::size_t i = 0;
  while (i < x.size()) {
    if (x[i] > 3) { i++; continue; }
    std::size_t j = i;
    while (j < x.size() && x[j] <= 3) j++;
    // Run [i, j): k-mers at i .. j - k.
    if (j - i >= static_cast<std::size_t>(k)) {
      const std::size_t count = j - i - static_cast<std::size_t>(k) + 1;
      std::vector<std::pair<std::uint64_t, bool>> h(count);  // hash, valid
      for (std::size_t p = 0; p < count; p++) {
        const Kmer f = encode_kmer(&x[i + p], k), rv = reverse_complement(f, k);
        h[p] = {hash64(std::min(f, rv), mask), f != rv};
      }
      const std::size_t windows = count >= static_cast<std::size_t>(w) ? count - static_cast<std::size_t>(w) + 1 : 1;
      const std::size_t width = std::min<std::size_t>(count, static_cast<std::size_t>(w));
      for (std::size_t a = 0; a < windows; a++) {
        std::uint64_t m = UINT64_MAX;
        bool any = false;
        for (std::size_t p = a; p < a + width; p++)
          if (h[p].second) { m = std::min(m, h[p].first); any = true; }
        if (!any) continue;
        for (std::size_t p = a; p < a + width; p++)
          if (h[p].second && h[p].first == m) out.insert({static_cast<std::uint32_t>(i + p), m});
      }
    }
    i = j;
  }
  return out;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x6b6d7221, 0, 0);

  // ---- k-mers ----------------------------------------------------------------------------------
  {
    bool code = true, rc = true, canon = true;
    for (int k = 1; k <= 32; k++)
      for (int trial = 0; trial < 20; trial++) {
        const Codes x = random_dna(r, static_cast<std::size_t>(k));
        const Kmer v = encode_kmer(x.data(), k);
        code = code && kmer_string(v, k) == decode_dna(x) && encode_kmer(decode_dna(x)) == v;
        const Codes y = reverse_complement(x);
        rc = rc && reverse_complement(v, k) == encode_kmer(y.data(), k) && reverse_complement(reverse_complement(v, k), k) == v;
        canon = canon && canonical(v, k) == canonical(encode_kmer(y.data(), k), k) && canonical(v, k) <= v;
      }
    s.require("k-mers of k = 1 to 32: encoded, decoded, the first base most significant", code);
    s.require("the reverse complement of a k-mer agrees with that of its codes, and is an involution", rc);
    s.require("a k-mer and its reverse complement have one canonical k-mer", canon);
    s.require("the order of k-mers is the lexicographic order of their letters", encode_kmer("ACGT") < encode_kmer("AGAA") && encode_kmer("TTTT") == 255);

    const Codes x = random_dna(r, 5000, 0.01);
    std::size_t visited = 0, expected = 0;
    bool same = true;
    for (int k : {1, 7, 21, 32}) {
      std::vector<std::size_t> positions;
      for_each_kmer(x, k, [&](std::size_t p, Kmer f, Kmer rv) {
        positions.push_back(p);
        same = same && f == encode_kmer(&x[p], k) && rv == reverse_complement(f, k);
      });
      std::vector<std::size_t> want;
      for (std::size_t p = 0; p + static_cast<std::size_t>(k) <= x.size(); p++)
        if (free_of_n(x, p, static_cast<std::size_t>(k))) want.push_back(p);
      visited += positions.size();
      expected += want.size();
      same = same && positions == want;
    }
    s.require("every k-mer free of N is visited once, with its reverse complement", same && visited == expected, std::to_string(visited) + " k-mers");
  }

  // ---- hash64 ----------------------------------------------------------------------------------
  {
    bool bijective = true;
    for (int k = 1; k <= 8; k++) {
      const Kmer mask = kmer_mask(k);
      std::vector<bool> hit(mask + 1, false);
      for (Kmer v = 0; v <= mask; v++) {
        const auto h = hash64(v, mask);
        if (h > mask || hit[h]) bijective = false;
        else hit[h] = true;
      }
    }
    s.require("hash64 is a bijection of [0, 4^k) for k = 1 to 8", bijective);
  }

  // ---- ntHash ----------------------------------------------------------------------------------
  {
    const Codes x = random_dna(r, 3000, 0.005);
    bool rolled = true, strands = true;
    std::size_t count = 0;
    for (unsigned k : {1u, 5u, 31u, 32u, 33u, 64u, 100u}) {
      std::vector<std::pair<std::size_t, std::uint64_t>> fwd;
      for_each_nthash(x, k, [&](std::size_t p, std::uint64_t f, std::uint64_t rv, std::uint64_t c) {
        const auto [df, dr] = nthash(&x[p], k);
        rolled = rolled && f == df && rv == dr && c == std::min(df, dr) && free_of_n(x, p, k);
        fwd.push_back({p, c});
        count++;
      });
      // On the reverse strand the k-mer at p is the one at n − k − p, with the same canonical hash.
      const Codes y = reverse_complement(x);
      std::vector<std::pair<std::size_t, std::uint64_t>> rev;
      for_each_nthash(y, k, [&](std::size_t p, std::uint64_t, std::uint64_t, std::uint64_t c) { rev.push_back({x.size() - k - p, c}); });
      std::sort(rev.begin(), rev.end());
      strands = strands && rev == fwd;
    }
    s.require("ntHash rolled equals ntHash computed directly, k from 1 to 100, across runs of N", rolled && count > 0, std::to_string(count) + " k-mers");
    s.require("ntHash's canonical value is the same on the reverse strand", strands);
  }

  // ---- minimizers ----------------------------------------------------------------------------------
  {
    bool definition = true, mirror = true;
    for (int trial = 0; trial < 30; trial++) {
      const int k = 3 + static_cast<int>(r.below(20)), w = 1 + static_cast<int>(r.below(12));
      // Small alphabets of k make ties and palindromes common.
      const Codes x = random_dna(r, 400 + r.below(600), trial % 3 == 0 ? 0.01 : 0);
      const auto mins = minimizers(x, k, w);
      std::set<std::pair<std::uint32_t, std::uint64_t>> got;
      bool ordered = true;
      for (std::size_t i = 0; i < mins.size(); i++) {
        got.insert({mins[i].position, mins[i].hash});
        if (i && mins[i].position <= mins[i - 1].position) ordered = false;
      }
      definition = definition && ordered && got == brute_minimizers(x, k, w);
      // On the reverse complement: the same hashes at mirrored positions.
      const auto back = minimizers(reverse_complement(x), k, w);
      std::set<std::pair<std::uint32_t, std::uint64_t>> mirrored;
      for (const auto& m : back) mirrored.insert({static_cast<std::uint32_t>(x.size() - static_cast<std::size_t>(k) - m.position), m.hash});
      mirror = mirror && mirrored == got;
    }
    s.require("minimizers agree with their definition by brute force (ties kept, palindromes skipped, N respected)", definition);
    s.require("the minimizers of the reverse strand are those of the forward strand, mirrored", mirror);
    const Codes mb = random_dna(r, 1'000'000);
    for (auto [w, k] : {std::pair{10, 21}, std::pair{5, 21}, std::pair{19, 19}}) {
      const double density = static_cast<double>(minimizers(mb, k, w).size()) / static_cast<double>(mb.size() - static_cast<std::size_t>(k) + 1);
      s.check("minimizer density on a random megabase, (w, k) = (" + std::to_string(w) + ", " + std::to_string(k) + "), against 2/(w + 1)",
              density / (2.0 / (w + 1)) - 1, 0.01);
    }
    // At k = 15 hash64's last step vanishes modulo 2^30 (kmer.hpp, Remark): the excess is bounded.
    const double d15 = static_cast<double>(minimizers(mb, 15, 10).size()) / static_cast<double>(mb.size() - 14);
    s.check("minimizer density at (w, k) = (10, 15), where hash64 is weaker: within 4 per cent above 2/(w + 1)", d15 / (2.0 / 11) - 1.02, 0.02);
  }

  // ---- syncmers ----------------------------------------------------------------------------------
  {
    bool definition = true;
    for (int trial = 0; trial < 20; trial++) {
      const int k = 8 + static_cast<int>(r.below(10)), sm = 3 + static_cast<int>(r.below(5));
      const int t = trial % 2 ? -1 : static_cast<int>(r.below(static_cast<std::uint64_t>(k - sm + 1)));
      const Codes x = random_dna(r, 600, 0.01);
      const auto got = syncmers(x, k, sm, t);
      std::vector<std::uint32_t> want;
      const Kmer mask = kmer_mask(sm);
      for (std::size_t p = 0; p + static_cast<std::size_t>(k) <= x.size(); p++) {
        if (!free_of_n(x, p, static_cast<std::size_t>(k))) continue;
        std::uint64_t best = UINT64_MAX;
        std::size_t at = 0;
        for (std::size_t j = 0; j + static_cast<std::size_t>(sm) <= static_cast<std::size_t>(k); j++) {
          const Kmer f = encode_kmer(&x[p + j], sm), rv = reverse_complement(f, sm);
          const auto h = hash64(std::min(f, rv), mask);
          if (h < best) { best = h; at = j; }
        }
        const std::size_t last = static_cast<std::size_t>(k - sm);
        if (t < 0 ? (at == 0 || at == last) : at == static_cast<std::size_t>(t)) want.push_back(static_cast<std::uint32_t>(p));
      }
      definition = definition && got == want;
    }
    s.require("closed and open syncmers agree with their definition", definition);
    const Codes mb = random_dna(r, 1'000'000);
    const double kmers = static_cast<double>(mb.size() - 14);
    const double k21 = static_cast<double>(mb.size() - 20);
    s.check("closed syncmer density (k = 21, s = 17) against 2/(k − s + 1)", static_cast<double>(syncmers(mb, 21, 17).size()) / k21 / (2.0 / 5) - 1, 0.01);
    s.check("open syncmer density (k = 21, s = 17, t = 2) against 1/(k − s + 1)", static_cast<double>(syncmers(mb, 21, 17, 2).size()) / k21 / (1.0 / 5) - 1, 0.01);
    (void)kmers;
  }

  return s.done();
}
