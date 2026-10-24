// Checks of kmers/filters.hpp: a Bloom filter has no false negatives, and its false-positive
// rate on keys never inserted lies within the binomial spread of (1 − e^{−hn/m})^h, for h
// from 1 to 12 hashes and at the optimum (m/n) ln 2, and within the spread of ρ^h given its
// fill ρ; the fill against 1 − e^{−hn/m}; a filter sized for a rate gives it; shards merge to
// the whole. A count-min sketch never under-counts a stream of k-mers drawn from a skewed law,
// over-counts by more than ε‖a‖₁ for at most a fraction δ of the keys, has a mean excess
// within ‖a‖₁/w, and sketches of shards sum to the sketch of the whole.

#include <genomics/kmers/filters.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

using namespace genomics;

int main() {
  calx::verification::Suite s;

  // Distinct random 62-bit keys (canonical 31-mers are below 2^62).
  calx::random::Stream r(1970, 0x4b455953, 0, 0);
  const std::size_t n = 200'000, queries = 2'000'000;
  std::vector<std::uint64_t> keys(n), others(queries);
  for (auto& x : keys) x = r.next64() >> 2;
  for (auto& x : others) x = (r.next64() >> 2) | (1ULL << 62);  // never among the keys

  // False-positive rates for h = 1 to 12 at m = 10 n bits.
  {
    const std::size_t m = 10 * n;
    bool negatives = true;
    double worst = 0, worst_fill = 0;
    for (unsigned h = 1; h <= 12; h++) {
      kmers::BloomFilter f(m, h, 3);
      for (auto x : keys) f.insert(x);
      for (auto x : keys) negatives = negatives && f.contains(x);
      std::size_t fp = 0;
      for (auto x : others) fp += f.contains(x);
      const double p = kmers::bloom_false_positive_rate(static_cast<double>(n), static_cast<double>(m), h);
      const double sd = std::sqrt(static_cast<double>(queries) * p * (1 - p));
      worst = std::max(worst, std::fabs(static_cast<double>(fp) - static_cast<double>(queries) * p) / sd);
      // Given the fill, each query is a false positive with probability ρ^h.
      const double rho = f.fill();
      const double fill_sd = std::sqrt(std::exp(-static_cast<double>(h * n) / static_cast<double>(m)) * (1 - std::exp(-static_cast<double>(h * n) / static_cast<double>(m))) / static_cast<double>(m));
      worst_fill = std::max(worst_fill, std::fabs(rho - (1 - std::exp(-static_cast<double>(h * n) / static_cast<double>(m)))) / fill_sd);
    }
    s.require("Bloom: no false negatives among 200000 keys, for h = 1 to 12", negatives);
    s.check("Bloom: the false-positive rate on 2·10^6 absent keys against (1 − e^{−hn/m})^h, h = 1 to 12, m = 10n (largest deviation in binomial σ)", worst, 4);
    s.check("Bloom: the fill against 1 − e^{−hn/m} (largest deviation in σ)", worst_fill, 5);
  }

  // The rate given the fill, at the optimal number of hashes.
  {
    const std::size_t m = 8 * n;
    const unsigned h = kmers::bloom_optimal_hashes(static_cast<double>(n), static_cast<double>(m));
    kmers::BloomFilter f(m, h, 4);
    for (auto x : keys) f.insert(x);
    std::size_t fp = 0;
    for (auto x : others) fp += f.contains(x);
    const double q = std::pow(f.fill(), h);
    const double sd = std::sqrt(static_cast<double>(queries) * q * (1 - q));
    s.require("Bloom: the optimal number of hashes for m = 8n is round(8 ln 2) = 6", h == 6);
    s.check("Bloom: the false-positive rate against ρ^h given the fill ρ (deviation in binomial σ)", (static_cast<double>(fp) - static_cast<double>(queries) * q) / sd, 4);
  }

  // Sized for a rate.
  {
    for (double p : {0.1, 0.01, 0.001}) {
      auto f = kmers::BloomFilter::for_capacity(n, p, 5);
      for (auto x : keys) f.insert(x);
      std::size_t fp = 0;
      for (auto x : others) fp += f.contains(x);
      const double want = kmers::bloom_false_positive_rate(static_cast<double>(n), static_cast<double>(f.bits()), f.hashes());
      const double sd = std::sqrt(static_cast<double>(queries) * want * (1 - want));
      s.check("Bloom: a filter for rate " + calx::verification::shortest(p) + " has it (relative to p)", calx::verification::relative(static_cast<double>(fp) / static_cast<double>(queries), p),
              0.05 + 4 * sd / static_cast<double>(queries) / p);
    }
  }

  // Shards merge to the whole.
  {
    kmers::BloomFilter whole(1 << 20, 5, 6), a(1 << 20, 5, 6), b(1 << 20, 5, 6);
    for (std::size_t i = 0; i < n; i++) {
      whole.insert(keys[i]);
      (i % 2 ? a : b).insert(keys[i]);
    }
    a.merge(b);
    bool same = true;
    for (std::size_t i = 0; i < 100000; i++) same = same && a.contains(others[i]) == whole.contains(others[i]);
    s.require("Bloom: the union of two shards' filters answers as the filter of all the keys", same && a.set_bits() == whole.set_bits());
  }

  // Count-min: a stream of 10^6 occurrences of 10^5 keys from a Zipf law of exponent 1.1.
  {
    const std::size_t distinct = 100'000, total = 1'000'000;
    std::vector<double> cdf(distinct);
    double z = 0;
    for (std::size_t i = 0; i < distinct; i++) cdf[i] = z += std::pow(static_cast<double>(i + 1), -1.1);
    for (auto& c : cdf) c /= z;
    calx::random::Stream q(2005, 0x5a495046, 0, 0);
    std::vector<std::uint64_t> stream(total);
    for (auto& x : stream) {
      const std::size_t i = static_cast<std::size_t>(std::lower_bound(cdf.begin(), cdf.end(), q.uniform()) - cdf.begin());
      x = keys[std::min(i, distinct - 1)];
    }
    const double epsilon = 1e-4, delta = 0.01;
    auto cms = kmers::CountMinSketch::for_error(epsilon, delta, 7);
    std::unordered_map<std::uint64_t, std::uint32_t> truth;
    for (auto x : stream) {
      cms.add(x);
      truth[x]++;
    }
    std::size_t under = 0, over = 0;
    double excess = 0;
    const double bound = cms.epsilon() * static_cast<double>(cms.total());
    for (const auto& [x, c] : truth) {
      const std::uint32_t e = cms.estimate(x);
      under += e < c;
      over += static_cast<double>(e - c) > bound;
      excess += static_cast<double>(e) - c;
    }
    s.require("count-min: w = ⌈e/ε⌉ = 27183 and d = ⌈ln(1/δ)⌉ = 5 for ε = 10^-4, δ = 0.01", cms.width() == 27183 && cms.depth() == 5);
    s.require("count-min: never under-counts, over " + std::to_string(truth.size()) + " keys of a skewed stream", under == 0);
    const double frac = static_cast<double>(over) / static_cast<double>(truth.size());
    s.check("count-min: the fraction of keys over-counted by more than ε‖a‖₁, at most δ = 0.01", frac, delta);
    s.check("count-min: the mean excess of the estimate, at most ‖a‖₁/w for one row", excess / static_cast<double>(truth.size()), static_cast<double>(total) / static_cast<double>(cms.width()));
    // Absent keys: their estimates are excesses alone.
    std::size_t absent_over = 0;
    for (std::size_t i = 0; i < 100000; i++) absent_over += cms.estimate(others[i]) > bound;
    s.check("count-min: absent keys estimated above ε‖a‖₁ (fraction), at most δ", static_cast<double>(absent_over) / 100000.0, delta);
    // One row of width w: a uniform hash gives a key the mean excess Σ_{y ≠ x} a(y)/w.
    {
      kmers::CountMinSketch row(1000, 1, 8);
      for (auto x : stream) row.add(x);
      double mean = 0, want = 0;
      for (const auto& [x, c] : truth) {
        mean += static_cast<double>(row.estimate(x)) - c;
        want += static_cast<double>(total - c) / 1000.0;
      }
      s.check("count-min: one row's mean excess against Σ_{y ≠ x} a(y)/w, the pairwise independent hash (relative)", calx::verification::relative(mean, want), 0.05);
    }
    // Shards.
    auto a = kmers::CountMinSketch::for_error(epsilon, delta, 7), b = kmers::CountMinSketch::for_error(epsilon, delta, 7);
    for (std::size_t i = 0; i < total; i++) (i % 3 ? a : b).add(stream[i]);
    a.merge(b);
    bool same = a.total() == cms.total();
    for (std::size_t i = 0; i < 1000; i++) same = same && a.estimate(keys[i]) == cms.estimate(keys[i]);
    s.require("count-min: the sum of two shards' sketches estimates as the sketch of the whole stream", same);
  }
  return s.done();
}
