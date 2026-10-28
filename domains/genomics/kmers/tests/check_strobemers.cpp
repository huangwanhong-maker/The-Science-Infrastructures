// Checks of the seeds that span more than they read: minstrobes and randstrobes of orders 2 to 4
// against their definitions by brute force, on random sequences with bases N; a constructed case
// in which a strobemer survives a deletion between its strobes while the k-mer of the same span
// does not, and the rates of survival measured over many random indels; and the density of the
// closed syncmers of sequence/kmer.hpp computed exactly, by enumerating every k-mer, for random
// sequence of equal base frequencies, against a random sequence of two megabases, and beside
// Edgar's 2/(k − s + 1), the density for a random order of distinct s-mers, from which a fixed
// hash and repeated s-mers make it depart a little (measured 2026-10-24: 0.33408 against 1/3 at
// (k, s) = (9, 4), 0.40376 and 0.39952 against 2/5 at (11, 7) and (12, 8)).

#include <genomics/kmers/strobemers.hpp>
#include <genomics/sequence/kmer.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::kmers;
using sequence::Codes;
using sequence::Kmer;

namespace {

std::string fixed(double v, int digits = 3) {
  char b[32];
  std::snprintf(b, sizeof b, "%.*f", digits, v);
  return b;
}

Codes random_dna(calx::random::Stream& r, std::size_t n, double n_rate = 0) {
  Codes x(n);
  for (auto& c : x) c = r.uniform() < n_rate ? sequence::N : static_cast<sequence::Code>(r.below(4));
  return x;
}

bool free_of_n(const Codes& x, std::size_t p, std::size_t len) {
  if (p + len > x.size()) return false;
  for (std::size_t j = p; j < p + len; j++)
    if (x[j] > 3) return false;
  return true;
}

std::uint64_t lmer_hash(const Codes& x, std::size_t p, int len) { return sequence::hash64(sequence::encode_kmer(&x[p], len), sequence::kmer_mask(len)); }

// Strobemers by their definition, every position, every window scanned.
std::vector<Strobemer> brute(const Codes& x, const StrobemerShape& s, bool random) {
  std::vector<Strobemer> out;
  const auto L = static_cast<std::size_t>(s.length);
  for (std::size_t i = 0; i + static_cast<std::size_t>((s.order - 1) * s.w_max) + L <= x.size(); i++) {
    if (!free_of_n(x, i, L)) continue;
    Strobemer m;
    m.position = static_cast<std::uint32_t>(i);
    m.order = s.order;
    m.strobes[0] = static_cast<std::uint32_t>(i);
    std::uint64_t H = lmer_hash(x, i, s.length);
    bool ok = true;
    for (int j = 2; j <= s.order && ok; j++) {
      const std::size_t a = i + static_cast<std::size_t>(s.w_min + (j - 2) * s.w_max), b = i + static_cast<std::size_t>((j - 1) * s.w_max);
      bool any = false;
      std::uint64_t best = 0, next = 0;
      std::size_t at = 0;
      for (std::size_t p = a; p <= b; p++) {
        if (!free_of_n(x, p, L)) continue;
        const std::uint64_t h = lmer_hash(x, p, s.length);
        const std::uint64_t key = random ? mix64(H ^ h) : h;
        if (!any || key < best) { best = key; at = p; next = mix64(H ^ h); any = true; }
      }
      ok = any;
      m.strobes[static_cast<std::size_t>(j - 1)] = static_cast<std::uint32_t>(at);
      H = next;
    }
    if (!ok) continue;
    m.hash = H;
    out.push_back(m);
  }
  return out;
}

bool same(const std::vector<Strobemer>& a, const std::vector<Strobemer>& b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); i++)
    if (a[i].position != b[i].position || a[i].hash != b[i].hash || a[i].order != b[i].order || a[i].strobes != b[i].strobes) return false;
  return true;
}

// The strobes' bases of a strobemer, end to end.
Codes strobe_bases(const Codes& x, const Strobemer& m, int length) {
  Codes out;
  for (int j = 0; j < m.order; j++) out.insert(out.end(), x.begin() + m.strobes[static_cast<std::size_t>(j)], x.begin() + m.strobes[static_cast<std::size_t>(j)] + length);
  return out;
}

bool occurs(const Codes& x, const Codes& w) { return std::search(x.begin(), x.end(), w.begin(), w.end()) != x.end(); }

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(2021, 0x5354524f, 0, 0);

  // ---- the definitions --------------------------------------------------------------------------
  {
    bool mins = true, rands = true;
    std::size_t total = 0;
    for (int trial = 0; trial < 60; trial++) {
      StrobemerShape sh;
      sh.order = 2 + trial % 3;
      sh.length = 5 + static_cast<int>(r.below(16));
      sh.w_min = 1 + static_cast<int>(r.below(30));
      sh.w_max = sh.w_min + static_cast<int>(r.below(30));
      const Codes x = random_dna(r, 300 + r.below(1500), trial % 4 ? 0.003 : 0.03);
      const auto a = minstrobes(x, sh), b = randstrobes(x, sh);
      mins = mins && same(a, brute(x, sh, false));
      rands = rands && same(b, brute(x, sh, true));
      total += a.size() + b.size();
    }
    s.require("minstrobes of orders 2 to 4 equal their definition (each window's least ℓ-mer hash, leftmost) on 60 random sequences with bases N", mins);
    s.require("randstrobes equal their definition (each window's ℓ-mer minimising mix64(H ⊕ h), H the strobes before), " + std::to_string(total) + " seeds in all", rands);
  }

  // ---- surviving an indel between strobes ------------------------------------------------------------
  {
    // The constructed case: order 2, ℓ = 15, window [i + 40, i + 70]; a deletion in the gap between
    // the first strobe and the window, at a strobemer whose second strobe is neither the window's
    // first ℓ-mer nor beaten by the ℓ-mer after the window's end.
    const StrobemerShape sh{2, 15, 40, 70};
    bool found = false, survives = false, kmer_lost = false;
    std::string detail;
    for (int trial = 0; trial < 50 && !found; trial++) {
      const Codes x = random_dna(r, 400);
      for (const auto& m : minstrobes(x, sh)) {
        const std::size_t i = m.position, p2 = m.strobes[1], first = i + 40, after = i + 71;
        if (after + 15 > x.size() || p2 == first || lmer_hash(x, after, 15) < lmer_hash(x, p2, 15)) continue;
        // Delete base i + 20, between the first strobe (i .. i + 14) and the window.
        Codes y = x;
        y.erase(y.begin() + static_cast<std::ptrdiff_t>(i + 20));
        const Codes span(x.begin() + static_cast<std::ptrdiff_t>(i), x.begin() + static_cast<std::ptrdiff_t>(p2 + 15));
        const auto ys = minstrobes(y, sh);
        const auto it = std::find_if(ys.begin(), ys.end(), [&](const Strobemer& n) { return n.position == i; });
        survives = it != ys.end() && it->hash == m.hash && strobe_bases(y, *it, 15) == strobe_bases(x, m, 15) && it->strobes[1] + 1 == p2;
        kmer_lost = !occurs(y, span);
        detail = "strobemer at " + std::to_string(i) + " with strobes at " + std::to_string(i) + " and " + std::to_string(p2) + ", the k-mer of its span " + std::to_string(span.size()) + " bases";
        found = true;
        break;
      }
    }
    s.require("constructed: a minstrobe survives a deletion between its strobes (the same strobes and hash, its second strobe one base earlier)", found && survives, detail);
    s.require("constructed: the k-mer of the same span, which covers the deletion, occurs nowhere in the mutated sequence", found && kmer_lost, detail);

    // Rates over random indels of one base in the gap before the window, for both kinds.
    std::size_t trials = 0, min_kept = 0, rand_kept = 0, kmers_kept = 0;
    for (int trial = 0; trial < 3000; trial++) {
      const Codes x = random_dna(r, 160);
      const std::size_t i = r.below(10);
      const std::size_t d = i + 15 + r.below(25);  // in [i + 15, i + 40): after strobe 1, before the window
      Codes y = x;
      if (trial % 2) y.erase(y.begin() + static_cast<std::ptrdiff_t>(d));
      else y.insert(y.begin() + static_cast<std::ptrdiff_t>(d), static_cast<sequence::Code>(r.below(4)));
      const auto at = [&](const std::vector<Strobemer>& v, std::size_t pos) { return std::find_if(v.begin(), v.end(), [&](const Strobemer& n) { return n.position == pos; }); };
      const auto xm = minstrobes(x, sh), ym = minstrobes(y, sh), xr = randstrobes(x, sh), yr = randstrobes(y, sh);
      const auto a = at(xm, i), b = at(ym, i), c = at(xr, i), e = at(yr, i);
      if (a == xm.end() || c == xr.end()) continue;
      trials++;
      min_kept += b != ym.end() && b->hash == a->hash;
      rand_kept += e != yr.end() && e->hash == c->hash;
      const Codes span(x.begin() + static_cast<std::ptrdiff_t>(i), x.begin() + static_cast<std::ptrdiff_t>(a->strobes[1] + 15));
      kmers_kept += occurs(y, span);
    }
    const double pm = static_cast<double>(min_kept) / trials, pr = static_cast<double>(rand_kept) / trials, pk = static_cast<double>(kmers_kept) / trials;
    // Measured 2026-10-24: 0.942 and 0.937 kept, against 1 − 2/(w + 1) = 0.938 for w = 31 (theory, Proposition 3.4); k-mers 0.
    s.require("over " + std::to_string(trials) + " random one-base indels between the first strobe and the window (ℓ = 15, w = 31): minstrobes and randstrobes kept at least 0.9 of the time, k-mers of the same span never",
              pm >= 0.9 && pr >= 0.9 && kmers_kept == 0, "minstrobes " + fixed(pm) + ", randstrobes " + fixed(pr) + ", k-mers " + fixed(pk) + "; 1 − 2/(w + 1) = " + fixed(1 - 2.0 / 32));
  }

  // ---- the exact density of closed syncmers --------------------------------------------------------
  for (auto [k, sm] : {std::pair{9, 4}, std::pair{11, 7}, std::pair{12, 8}}) {
    // Every k-mer at once: it is a closed syncmer when its first minimal canonical s-mer (by hash64)
    // is its first or its last. The density on a random sequence of equal frequencies is the
    // fraction of k-mers that are, exactly, by linearity of expectation.
    const int m = k - sm + 1;
    const Kmer smask = sequence::kmer_mask(sm);
    std::uint64_t closed = 0;
    const std::uint64_t all = 1ULL << (2 * k);
    for (std::uint64_t x = 0; x < all; x++) {
      std::uint64_t best = UINT64_MAX;
      int at = 0;
      for (int j = 0; j < m; j++) {
        const Kmer f = (x >> (2 * (k - sm - j))) & smask;
        const Kmer rc = sequence::reverse_complement(f, sm);
        const std::uint64_t h = sequence::hash64(f < rc ? f : rc, smask);
        if (h < best) { best = h; at = j; }
      }
      closed += at == 0 || at == m - 1;
    }
    const double exact = static_cast<double>(closed) / static_cast<double>(all), edgar = 2.0 / m;
    const Codes mb = random_dna(r, 2'000'000);
    const double positions = static_cast<double>(mb.size() - static_cast<std::size_t>(k) + 1);
    const double measured = static_cast<double>(sequence::syncmers(mb, k, sm).size()) / positions;
    // The positions are (k − 1)-dependent, so the variance of the count is at most (2k − 1) n p (1 − p).
    const double sigma = std::sqrt((2.0 * k - 1) * exact * (1 - exact) / positions);
    s.check("closed syncmers (k, s) = (" + std::to_string(k) + ", " + std::to_string(sm) + "): density on a random 2 Mb against its exact value by enumeration of 4^k k-mers, " + fixed(exact, 5) +
                " (Edgar's 2/(k − s + 1) = " + fixed(edgar, 5) + "), in σ",
            (measured - exact) / sigma, 5.0);
  }
  return s.done();
}
