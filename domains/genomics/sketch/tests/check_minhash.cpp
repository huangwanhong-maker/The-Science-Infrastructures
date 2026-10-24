// Checks of MinHash and FracMinHash sketches on simulated genomes of a megabase: the bottom-s sketch
// is the s smallest hashes of the exact k-mer set, the same on the reverse strand; the Jaccard
// estimate within its binomial error of the exact Jaccard of the k-mer sets, for pairs diverged
// by substitutions at known rates; the Mash distance within Ondov et al.'s Table 1 bounds of the
// distance of the exact Jaccard; the Mash distance tracking −ln(1 − p) under independent
// substitutions and p under Mash's Poisson model, for p ≤ 0.2; the p-value's binomial tail against
// a closed form and a direct sum, and its null model of chance matches among canonical k-mers;
// FracMinHash's sketch size, the containment of a subset (1
// exactly), containment within Hera et al.'s Theorem 3 variance, and ANI = C^{1/k} recovering
// 1 − p.

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/kmer.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/sketch/minhash.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace genomics;
using sequence::Codes;

namespace {

constexpr std::uint32_t CHECK = 0x4d484153;

// The canonical k-mers of a sequence, sorted and distinct.
std::vector<sequence::Kmer> kmer_set(const Codes& s, int k) {
  std::vector<sequence::Kmer> v;
  sequence::for_each_kmer(s, k, [&](std::size_t, sequence::Kmer f, sequence::Kmer r) { v.push_back(std::min(f, r)); });
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
  return v;
}
std::size_t common(const std::vector<sequence::Kmer>& a, const std::vector<sequence::Kmer>& b) {
  std::size_t i = 0, j = 0, n = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i] < b[j]) i++;
    else if (b[j] < a[i]) j++;
    else { n++; i++; j++; }
  }
  return n;
}

// Each base substituted with probability q (simulate::substitute: a different base, uniformly).
Codes mutate(const Codes& x, double q, std::uint64_t seed, std::uint32_t object) {
  calx::random::Stream r(seed, CHECK, object, 1);
  Codes y = x;
  for (auto& c : y)
    if (r.uniform() < q) c = simulate::substitute(r, c);
  return y;
}

// Ondov et al. 2016, Table 1 (k = 21, s = 1000): the 0.99 bound of the Mash distance's error at
// the tabulated distance at or above D.
double table1_bound(double D) {
  if (D <= 0.05) return 0.0068;
  if (D <= 0.10) return 0.0158;
  if (D <= 0.15) return 0.0323;
  return 0.0630;
}

std::string fixed(double x, int digits = 4) {
  char b[48];
  std::snprintf(b, sizeof b, "%.*f", digits, x);
  return b;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const int k = 21;
  const std::size_t L = 1'000'000, size = 1000;
  const Codes a = simulate::random_sequence(L, 0.5, 20261024);
  const auto set_a = kmer_set(a, k);
  const auto sketch_a = sketch::minhash(a, k, size);

  // ---- the sketch ----------------------------------------------------------------------------------
  {
    std::vector<std::uint64_t> all;
    for (auto x : set_a) all.push_back(sketch::kmer_hash(x, 42));
    std::sort(all.begin(), all.end());
    s.require("the bottom-s sketch is the s smallest hashes of the exact set of canonical k-mers",
              sketch_a.hashes == std::vector<std::uint64_t>(all.begin(), all.begin() + static_cast<std::ptrdiff_t>(size)) && sketch_a.kmers == L - k + 1);
    s.require("the sketch of the reverse strand is the same", sketch::minhash(sequence::reverse_complement(a), k, size).hashes == sketch_a.hashes);
    s.check("the distinct k-mers estimated from the sketch, s·2^64/v, against the exact count (relative, 4/√s)",
            (sketch_a.cardinality() - static_cast<double>(set_a.size())) / static_cast<double>(set_a.size()), 4 / std::sqrt(static_cast<double>(size)));
    const auto small = sketch::minhash(Codes(a.begin(), a.begin() + 500), k, size);
    s.require("a sketch that is not full counts its k-mers exactly", !small.full() && small.cardinality() == static_cast<double>(kmer_set(Codes(a.begin(), a.begin() + 500), k).size()));
  }

  // ---- Jaccard and the Mash distance at known divergence ----------------------------------------------------------------------
  {
    std::uint32_t object = 0;
    for (bool poisson : {false, true})
      for (double p : {0.001, 0.01, 0.05, 0.10, 0.15, 0.20}) {
        // The simple model substitutes each base with probability p; the Poisson model alters a
        // base when a Poisson process of rate p has an event there, with probability 1 − e^{−p}.
        const double q = poisson ? -std::expm1(-p) : p;
        const Codes b = mutate(a, q, 20261024, object++);
        const auto set_b = kmer_set(b, k);
        const std::size_t both = common(set_a, set_b);
        const double J = static_cast<double>(both) / static_cast<double>(set_a.size() + set_b.size() - both);
        const auto c = sketch::compare(sketch_a, sketch::minhash(b, k, size));
        const double sd = std::sqrt(J * (1 - J) / static_cast<double>(size));
        const std::string model = poisson ? "Poisson model, d = " : "substitutions, p = ";
        s.check("Jaccard from the sketches against the exact Jaccard (" + model + fixed(p, 3) + ", J = " + fixed(J) + "), in binomial standard errors",
                (c.jaccard - J) / std::max(sd, 1e-9), 3.5, "σ");
        const double D = sketch::mash_distance(J, k);
        // Table 1 bounds the error up to D = 0.20 for s = 1000 (beyond, a sketch may share nothing).
        if (D <= 0.20)
          s.check("the Mash distance of the sketches against that of the exact Jaccard (" + model + fixed(p, 3) + "), within Table 1 of Ondov et al.",
                  c.distance - D, table1_bound(D));
        // The model's distance and the noise of the mutation process: the conserved fraction w is
        // a sum of about L w k-mers in runs of about k, so D carries a relative error of about
        // √(k / (L w)) / k.
        const double w = std::pow(1 - q, k), want = poisson ? p : -std::log1p(-p);
        const double noise = 4 * std::sqrt(static_cast<double>(k) / (static_cast<double>(L) * w)) / k + 2e-4;
        s.check(poisson ? "under the Poisson model the Mash distance of the exact sets tracks d (d = " + fixed(p, 3) + ", D = " + fixed(D) + ")"
                        : "under substitutions the Mash distance of the exact sets tracks −ln(1 − p) > p (p = " + fixed(p, 3) + ", D = " + fixed(D) + ")",
                D - want, noise);
      }
  }

  // ---- the p-value ----------------------------------------------------------------------------------
  {
    s.check("the binomial tail P(X ≥ 5) for X ~ Bin(10, 1/2) = 638/1024", sketch::binomial_upper_tail(5, 10, 0.5) - 638.0 / 1024, 1e-14);
    double direct = 0;
    for (std::size_t i = 40; i <= 1000; i++)
      direct += std::exp(std::lgamma(1001.0) - std::lgamma(i + 1.0) - std::lgamma(1001.0 - i) + static_cast<double>(i) * std::log(0.03) + static_cast<double>(1000 - i) * std::log1p(-0.03));
    s.check("the binomial tail P(X ≥ 40) for Bin(1000, 0.03) against the direct sum (relative)", sketch::binomial_upper_tail(40, 1000, 0.03) / direct - 1, 1e-10);
    // Unrelated genomes at k = 12, where chance matches are common: r ≈ 0.03.
    const int k12 = 12;
    const Codes x = simulate::random_sequence(L, 0.5, 1), y = simulate::random_sequence(L, 0.5, 2);
    const auto sx = sketch::minhash(x, k12, size), sy = sketch::minhash(y, k12, size);
    const auto c = sketch::compare(sx, sy);
    // The expected Jaccard index r of eq. 5, with the space of canonical 12-mers, and with Mash's 4^k.
    const auto expected = [&](double space) {
      const double px = sketch::random_kmer_probability(sx.cardinality(), space), py = sketch::random_kmer_probability(sy.cardinality(), space);
      return px * py / (px + py - px * py) * static_cast<double>(c.considered);
    };
    const double canonical = expected(sketch::canonical_kmer_space(k12)), mash = expected(std::pow(4.0, k12));
    s.require("two unrelated genomes share about r·s hashes by chance, r from the canonical 12-mers, and are not significant (k = 12)",
              std::fabs(static_cast<double>(c.shared) - canonical) < 4 * std::sqrt(canonical) && c.p_value > 1e-3,
              std::to_string(c.shared) + " shared, r·s = " + fixed(canonical, 1) + " (with Mash's 4^k, " + fixed(mash, 1) + "), p = " + fixed(c.p_value, 3));
    s.require("with Mash's space of 4^k k-mers the same count would be significant (the Caveat of minhash.hpp)",
              sketch::mash_p_value(c.shared, c.considered, sx.cardinality(), sy.cardinality(), std::pow(4.0, k12)) < 1e-3);
    const auto related = sketch::compare(sketch_a, sketch::minhash(mutate(a, 0.05, 7, 99), k, size));
    s.require("a genome and its copy diverged by 5 per cent are significant (k = 21)", related.p_value < 1e-100, std::to_string(related.shared) + " shared");
  }

  // ---- FracMinHash ----------------------------------------------------------------------------------
  {
    const std::uint64_t scaled = 200;
    const double frac = 1.0 / static_cast<double>(scaled);
    const auto fa = sketch::frac_minhash(a, k, scaled);
    bool definition = true;
    std::size_t expected = 0;
    for (auto x : set_a) expected += sketch::kmer_hash(x, 42) <= fa.max_hash();
    definition = fa.hashes.size() == expected && std::all_of(fa.hashes.begin(), fa.hashes.end(), [&](std::uint64_t h) { return h <= fa.max_hash(); });
    s.require("FracMinHash keeps every hash at most H/scaled of the k-mer set", definition, std::to_string(fa.hashes.size()) + " hashes");
    const double n = static_cast<double>(set_a.size());
    s.check("its size against |A|·s, in binomial standard errors", (static_cast<double>(fa.hashes.size()) - n * frac) / std::sqrt(n * frac * (1 - frac)), 4, "σ");

    // A subset: the first 400 kb of A is contained in A, exactly.
    const Codes part(a.begin(), a.begin() + 400000);
    const auto fp = sketch::frac_minhash(part, k, scaled);
    s.require("the containment of a subset in its superset is 1 exactly", sketch::containment(fp, fa) == 1);
    const auto set_p = kmer_set(part, k);
    const double C = static_cast<double>(set_p.size()) / n;  // C(A, part) = |part| / |A|
    const double var = sketch::containment_variance(static_cast<double>(set_p.size()), n - static_cast<double>(set_p.size()), frac);
    s.check("the containment of A in its 40 per cent prefix against |A ∩ B|/|A|, in standard errors of Theorem 3",
            (sketch::containment(fa, fp) - C) / std::sqrt(var), 4, "σ");
    s.check("the debiased containment (eq. 3) is the same for a set this large, to 10⁻¹²",
            sketch::containment_debiased(fa, fp, n) - sketch::containment(fa, fp), 1e-12);

    for (double p : {0.001, 0.01, 0.05, 0.1}) {
      const Codes b = mutate(a, p, 11, static_cast<std::uint32_t>(p * 1e4));
      const auto set_b = kmer_set(b, k);
      const auto fb = sketch::frac_minhash(b, k, scaled);
      const double shared = static_cast<double>(common(set_a, set_b));
      const double c_exact = shared / n;
      const double sd = std::sqrt(sketch::containment_variance(shared, n - shared, frac));
      const double c_hat = sketch::containment(fa, fb);
      s.check("FracMinHash containment against the exact containment (p = " + fixed(p, 3) + "), in standard errors of Theorem 3", (c_hat - c_exact) / sd, 4, "σ");
      // ANI = C^{1/k}: the error of Ĉ carried through the k-th root, and the mutation process's.
      const double ani = sketch::ani(c_hat, k);
      const double tol = 4 * sd * ani / (k * c_hat) + 4 * std::sqrt(static_cast<double>(k) / (static_cast<double>(L) * c_exact)) / k + 1e-4;
      s.check("ANI = Ĉ^{1/k} recovers 1 − p (p = " + fixed(p, 3) + ", ANI = " + fixed(ani, 5) + ")", ani - (1 - p), tol);
      const double J = shared / (n + static_cast<double>(set_b.size()) - shared);
      const double u = static_cast<double>(fa.hashes.size() + fb.hashes.size()) - static_cast<double>(sketch::intersection(fa, fb));
      s.check("FracMinHash Jaccard against the exact Jaccard (p = " + fixed(p, 3) + "), in binomial standard errors",
              (sketch::frac_jaccard(fa, fb) - J) / std::max(std::sqrt(J * (1 - J) / u), 1e-9), 4, "σ");
    }
  }

  return s.done();
}
