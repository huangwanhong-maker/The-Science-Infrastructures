// Checks of HyperLogLog: over 300 sets of distinct elements from 10 to 3·10⁵, the relative error of
// Ertl's improved estimate has a root mean square within 1.04/√m (the bound of Flajolet et al.)
// and no bias beyond its sampling error, at p = 10 and p = 14; the classical estimate the same away
// from its transition; inserting again changes nothing; the sketch of a union is the register-wise
// maximum; Jaccard by inclusion and exclusion against the exact Jaccard; and the distinct k-mers of
// a simulated genome, the same on its reverse strand.

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/kmer.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/sketch/hyperloglog.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace genomics;

namespace {

constexpr std::uint32_t CHECK = 0x484c4c43;

std::string fixed(double x, int digits = 4) {
  char b[48];
  std::snprintf(b, sizeof b, "%.*f", digits, x);
  return b;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, CHECK, 0, 0);

  // ---- the error over many sets ----------------------------------------------------------------------------------
  for (int p : {10, 14}) {
    const double m = std::ldexp(1.0, p), bound = 1.04 / std::sqrt(m);
    double sum2 = 0, sum = 0, sum2_classic = 0;
    std::size_t sets = 0, classic_sets = 0;
    std::uint64_t base = 0;
    for (int t = 0; t < 300; t++) {
      // Cardinalities spread on a log scale; elements are distinct integers, hashed (a bijection).
      const auto n = static_cast<std::uint64_t>(std::exp(std::log(10.0) + r.uniform() * (std::log(3e5) - std::log(10.0))));
      sketch::HyperLogLog h(p, 1000 + static_cast<std::uint64_t>(t));
      for (std::uint64_t i = 0; i < n; i++) h.insert(sketch::kmer_hash(base + i, h.seed()));
      base += n;
      const double e = (h.estimate() - static_cast<double>(n)) / static_cast<double>(n);
      sum += e;
      sum2 += e * e;
      sets++;
      if (static_cast<double>(n) > 6 * m) {
        const double c = (h.estimate_classic() - static_cast<double>(n)) / static_cast<double>(n);
        sum2_classic += c * c;
        classic_sets++;
      }
    }
    const double rms = std::sqrt(sum2 / static_cast<double>(sets)), mean = sum / static_cast<double>(sets);
    s.check("p = " + std::to_string(p) + ": the root mean square relative error of Ertl's estimate over 300 sets, against 1.04/√m = " + fixed(bound),
            std::max(0.0, rms / bound - 1), 0.15);
    s.check("p = " + std::to_string(p) + ": its mean relative error, in standard errors of the mean", mean / (bound / std::sqrt(static_cast<double>(sets))), 4, "σ");
    if (classic_sets > 20)
      s.check("p = " + std::to_string(p) + ": the classical estimate above 6m, root mean square against 1.04/√m (" + std::to_string(classic_sets) + " sets)",
              std::max(0.0, std::sqrt(sum2_classic / static_cast<double>(classic_sets)) / bound - 1), 0.2);
  }

  // ---- idempotence, union, Jaccard ----------------------------------------------------------------------------------
  {
    const int p = 14;
    sketch::HyperLogLog a(p);
    const std::uint64_t n = 100000;
    for (std::uint64_t i = 0; i < n; i++) a.insert(sketch::kmer_hash(i, 42));
    const auto before = a.registers();
    for (std::uint64_t i = 0; i < n; i++) a.insert(sketch::kmer_hash(i, 42));
    s.require("inserting the same elements again changes no register", a.registers() == before);
    bool jaccard_ok = true;
    std::string detail;
    for (double J : {0.05, 0.2, 0.5, 0.8, 1.0}) {
      // |A| = |B| = n with |A ∩ B| = 2nJ/(1 + J).
      const auto shared = static_cast<std::uint64_t>(std::llround(2 * static_cast<double>(n) * J / (1 + J)));
      sketch::HyperLogLog x(p), y(p), u(p);
      for (std::uint64_t i = 0; i < n; i++) x.insert(sketch::kmer_hash(i, 42));
      for (std::uint64_t i = n - shared; i < 2 * n - shared; i++) y.insert(sketch::kmer_hash(i, 42));
      for (std::uint64_t i = 0; i < 2 * n - shared; i++) u.insert(sketch::kmer_hash(i, 42));
      jaccard_ok = jaccard_ok && sketch::union_of(x, y).registers() == u.registers();
      const double exact = static_cast<double>(shared) / static_cast<double>(2 * n - shared);
      const double est = sketch::jaccard(x, y);
      // Three estimates of relative error ε each: |ΔJ| ≲ (1 + J)·ε·√3.
      const double eps = 1.04 / std::sqrt(std::ldexp(1.0, p));
      s.check("Jaccard by inclusion and exclusion against the exact J = " + fixed(exact, 3) + " (p = 14)", est - exact, 4 * (1 + exact) * eps * std::sqrt(3.0));
      detail += (detail.empty() ? "" : ", ") + fixed(est, 3);
    }
    s.require("the sketch of a union is the register-wise maximum of the sketches, exactly", jaccard_ok, "estimates " + detail);
  }

  // ---- k-mers of a genome ----------------------------------------------------------------------------------
  {
    const sequence::Codes g = simulate::random_sequence(2'000'000, 0.45, 5);
    std::vector<sequence::Kmer> v;
    sequence::for_each_kmer(g, 21, [&](std::size_t, sequence::Kmer f, sequence::Kmer rv) { v.push_back(std::min(f, rv)); });
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    sketch::HyperLogLog h(12), rc(12);
    h.add(g, 21);
    rc.add(sequence::reverse_complement(g), 21);
    s.check("the distinct 21-mers of a 2 Mb genome, in relative standard errors 1.04/√m (p = 12)",
            (h.estimate() - static_cast<double>(v.size())) / static_cast<double>(v.size()) / h.standard_error(), 4, "σ");
    s.require("the sketch of the reverse strand is the same", h.registers() == rc.registers());
    s.require("an empty sketch estimates 0", sketch::HyperLogLog(10).estimate() == 0 && sketch::HyperLogLog(10).estimate_classic() == 0);
  }

  return s.done();
}
