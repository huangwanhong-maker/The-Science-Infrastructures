// The Wedderburn channels of a genome's statistics: the isotypic decomposition of its k-mer
// frequency vector under a group acting on k-mers (theory/group-grammars.md §7).
//
// A group G acting on the set X of k-mers acts linearly on functions f: X → ℝ by
// (g·f)(x) = f(g⁻¹·x), a permutation representation, and the central idempotents of ℂ[G] project
// it onto its isotypic components,
//
//   P_χ f = (χ(1)/|G|) Σ_g χ(g)‾ (g·f),     f = Σ_χ P_χ f,   ‖f‖² = Σ_χ ‖P_χ f‖²
//
// (Serre 1977, Linear Representations of Finite Groups, §2.6, Theorem 8; the components are
// orthogonal for the inner product of functions on X, which the permutation action preserves).
// This is where the draft's "parallel channels, one per Wedderburn block" holds exactly: a linear
// object, the frequency vector, splits into |Ĝ| orthogonal channels, and the energy in each says
// how much of the statistics transforms by that character. Under the group of order 8 of the
// bases and strands (strands.hpp), which is abelian with ±1 characters, a character with
// χ(κρ) = −1 collects the part of f that is odd under the reverse complement: Chargaff's second
// rule says it is small, and the channels measure how small.
//
// The characters come exactly from morphe's character table (Dixon–Schneider); the frequencies
// are real numbers, so the projections are computed in double precision, to rounding.

#pragma once

#include <morphe/algebra/characters.hpp>

#include <cmath>
#include <complex>
#include <cstdint>
#include <string>
#include <vector>

#include "strands.hpp"

namespace genomics::grammar {

struct Channel {
  std::size_t character = 0;
  std::vector<int> signs;  // χ(g) for each element, when the character is ±1-valued (else empty)
  double energy = 0;       // ‖P_χ f‖²
  double fraction = 0;     // of ‖f‖²
};

// The isotypic decomposition of the k-mer frequencies of x under a group acting on words.
inline std::vector<Channel> kmer_channels(const sequence::Codes& x, std::size_t k, const WordGroup& w, double* total = nullptr) {
  const auto f = kmer_frequencies(x, k);
  const auto table = ma::character_table(w.group);
  const std::size_t n = w.group.order();
  // The permutation of the k-mers by each element.
  std::vector<std::vector<std::size_t>> perm(n, std::vector<std::size_t>(f.size()));
  for (std::uint32_t g = 0; g < n; g++)
    for (std::size_t u = 0; u < f.size(); u++) perm[g][u] = act_on_kmer(w.words, g, u, k);
  double norm = 0;
  for (double v : f) norm += v * v;
  if (total) *total = norm;
  std::vector<Channel> out;
  for (std::size_t chi = 0; chi < table.size(); chi++) {
    std::vector<std::complex<double>> p(f.size(), 0.0);
    Channel c;
    c.character = chi;
    bool signs = true;
    for (std::uint32_t g = 0; g < n; g++) {
      const std::complex<double> value = table.value(chi, g).to_complex();
      const std::complex<double> weight = static_cast<double>(table.degree(chi)) / static_cast<double>(n) * std::conj(value);
      // (g·f)(x) = f(g⁻¹·x): the value at g·u is f(u).
      for (std::size_t u = 0; u < f.size(); u++) p[perm[g][u]] += weight * f[u];
      if (std::abs(value.imag()) > 1e-12 || std::abs(std::abs(value.real()) - 1) > 1e-12) signs = false;
      c.signs.push_back(value.real() > 0 ? 1 : -1);
    }
    if (!signs) c.signs.clear();
    for (const auto& v : p) c.energy += std::norm(v);
    c.fraction = norm > 0 ? c.energy / norm : 0;
    out.push_back(c);
  }
  return out;
}

}  // namespace genomics::grammar
