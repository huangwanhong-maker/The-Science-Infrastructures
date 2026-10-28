// The genetic code as a function on a group: codons are a torsor of V₄³ = (C₂)⁶ (nucleotides.hpp),
// a code assigns each codon an amino acid or a stop, and its symmetries are the elements m with
// aa(c ⊕ m) = aa(c) for every codon c (theory/group-grammars.md §6).
//
// The invariance profile I(m) = #{c : aa(c ⊕ m) = aa(c)}/64 measures how far each element is a
// symmetry; the elements with I(m) = 1 form a subgroup, the stabiliser of the code. The profile
// is the autocorrelation of the indicator functions of the amino acids, I(m) = (1/64) Σ_a Σ_c
// f_a(c) f_a(c ⊕ m), so its Walsh–Hadamard transform is their power spectrum, Σ_a |f̂_a|²
// (the Wiener–Khinchin theorem on the group (C₂)⁶; Terras 1999, Fourier Analysis on Finite Groups
// and Applications, Cambridge, ch. 2): the spectrum says how much of the code's structure lies in
// characters that ignore the third base.
//
// Exact facts, checked against NCBI's tables (sequence/genetic_code.hpp):
// - the standard code (table 1) is invariant under no element but 1, and the third-position
//   transition τ₃ keeps the amino acid of 60 codons of 64, failing only at AUA/AUG (Ile/Met) and
//   UGA/UGG (stop/Trp);
// - the vertebrate mitochondrial code (table 2), which reads AUA as Met and UGA as Trp, is
//   invariant under τ₃: its stabiliser is {1, τ₃}.
// Both follow from the tables; they are not typed in.

#pragma once

#include <genomics/sequence/genetic_code.hpp>

#include <morphe/algebra/abelian_fourier.hpp>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "nucleotides.hpp"

namespace genomics::grammar {

// The amino acid (or '*') of each codon in the six-bit numbering of nucleotides.hpp.
inline std::array<char, 64> code_on_group(const sequence::GeneticCode& code) {
  std::array<char, 64> out{};
  for (std::uint32_t c = 0; c < 64; c++) {
    const auto b = codon_bases(c);
    out[c] = sequence::amino_acid(code, b[0], b[1], b[2]);
  }
  return out;
}

// I(m) as the number of codons (of 64) whose amino acid m keeps.
inline std::array<int, 64> invariance_profile(const std::array<char, 64>& aa) {
  std::array<int, 64> out{};
  for (std::uint32_t m = 0; m < 64; m++)
    for (std::uint32_t c = 0; c < 64; c++) out[m] += aa[c] == aa[c ^ m];
  return out;
}

// The stabiliser of the code: the elements with I(m) = 64.
inline std::vector<std::uint32_t> stabiliser(const std::array<char, 64>& aa) {
  const auto p = invariance_profile(aa);
  std::vector<std::uint32_t> out;
  for (std::uint32_t m = 0; m < 64; m++)
    if (p[m] == 64) out.push_back(m);
  return out;
}

// The power spectrum of the code: P(u) = Σ_a f̂_a(u)², f_a the indicator of amino acid a (or the
// stop) on (C₂)⁶ and f̂ its Walsh–Hadamard transform (morphe, exact over ℤ). By the
// Wiener–Khinchin identity on (C₂)⁶, P is the Walsh–Hadamard transform of the invariance profile:
// P(u) = Σ_m (−1)^{u·m} I(m), exactly in integers.
inline std::array<long long, 64> code_power_spectrum(const std::array<char, 64>& aa) {
  std::array<long long, 64> p{};
  std::string letters;
  for (char a : aa)
    if (letters.find(a) == std::string::npos) letters += a;
  for (char a : letters) {
    std::vector<long long> f(64, 0);
    for (std::uint32_t c = 0; c < 64; c++) f[c] = aa[c] == a;
    const auto fh = ma::walsh_hadamard(f);
    for (std::size_t u = 0; u < 64; u++) p[u] += fh[u] * fh[u];
  }
  return p;
}

// The characters u that ignore the bases at the given positions (their two bits of u are zero).
inline bool ignores_positions(std::uint32_t u, const std::vector<int>& positions) {
  for (int p : positions)
    if ((u >> (2 * (2 - p))) & 3u) return false;
  return true;
}

// The codons where an element fails, as pairs (c, c ⊕ m) with c < c ⊕ m, in letters.
inline std::vector<std::string> failures(const std::array<char, 64>& aa, std::uint32_t m) {
  std::vector<std::string> out;
  for (std::uint32_t c = 0; c < 64; c++) {
    const std::uint32_t d = c ^ m;
    if (c < d && aa[c] != aa[d]) {
      std::string s;
      for (auto b : codon_bases(c)) s += sequence::dna_letter(b);
      s += '/';
      for (auto b : codon_bases(d)) s += sequence::dna_letter(b);
      out.push_back(s);
    }
  }
  return out;
}

}  // namespace genomics::grammar
