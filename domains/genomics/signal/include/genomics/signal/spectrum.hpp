// The spectra of DNA: the indicator sequences of the four bases, their discrete Fourier
// transforms, the spectrum S(k) = Σ_b |U_b(k)|², and the period-3 signal of protein-coding DNA,
// over a whole sequence or in windows along a genome.
//
// - Indicator sequences (Voss 1992, "Evolution of long-range fractal correlations and 1/f noise
//   in DNA base sequences", Phys. Rev. Lett. 68, 3805–3808): u_b(n) = [x_n = b] for the four bases
//   b; a base N is 0 in all four. Their transforms U_b(k) = Σ_n u_b(n) e^{−2πi nk/N}, k ∈ ℤ_N, are
//   computed by calx's FFT (radix 2, or Bluestein's chirp for other N), and S(k) is the sum of
//   the four power spectra. S(N − k) = S(k), and S is the same for a sequence and its reverse
//   complement (theory §1.2).
// - Parseval: Σ_k S(k) = N ν, ν the number of bases other than N, so the mean S̄ of the spectrum
//   over the N frequencies is ν exactly.
// - The period-3 signal (Tiwari, Ramachandran, Bhattacharya, Bhattacharya and Ramaswamy 1997,
//   "Prediction of probable genes by Fourier analysis of genomic sequences", CABIOS 13, 263–270;
//   Anastassiou 2001, "Genomic signal processing", IEEE Signal Processing Magazine 18(4), 8–20):
//   in protein-coding DNA the bases are used unequally at the three positions of codons, and S
//   has a peak at k = N/3. When 3 | N the character n ↦ e^{−2πi n/3} of ℤ_N factors through
//   ℤ_N → ℤ_3, so that U_b(N/3) = c_{b,0} + c_{b,1} ζ + c_{b,2} ζ², ζ = e^{−2πi/3}, with c_{b,r} the
//   number of b at the positions ≡ r (mod 3). Its squared modulus is the norm of an Eisenstein
//   integer, c₀² + c₁² + c₂² − c₀c₁ − c₁c₂ − c₂c₀ = ½ Σ_{r<s} (c_r − c_s)², so S(N/3) is an integer
//   computed exactly from twelve counts, with no transform (theory §1.3). A window sliding by one
//   base updates two counts, so the profile along a genome of n bases costs O(n).
// - Its relative strength is P = S(N/3) / S̄ = S(N/3) / ν, the peak over the mean of the spectrum
//   (the measure of Tiwari et al. 1997). For independent bases of probabilities p_b,
//   E S(N/3) = N (1 − Σ_b p_b²) (theory, Proposition 1.7): P is near 3/4 in random sequence of
//   equal frequencies, and larger in coding sequence.
//
// The transforms are in double precision, S(k) with an error of O(ε log N) relative to N², as
// calx's FFT; they hold four transforms of N complex numbers in memory, O(N log N) time. The
// period-3 values of windows are exact integers.

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <calx/transforms/fft.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace genomics::signal {

using sequence::Code;
using sequence::Codes;
using calx::transforms::Complex;

// The indicator sequence of base b (0 to 3): 1 where x is b, 0 elsewhere.
inline std::vector<double> indicator(const Codes& x, Code b) {
  if (b > 3) throw std::invalid_argument("genomics::signal::indicator: a base, 0 to 3");
  std::vector<double> u(x.size());
  for (std::size_t n = 0; n < x.size(); n++) u[n] = x[n] == b ? 1.0 : 0.0;
  return u;
}

// U_b(k), k = 0 … N − 1, for the four bases.
inline std::array<std::vector<Complex>, 4> indicator_transforms(const Codes& x) {
  std::array<std::vector<Complex>, 4> U;
  for (Code b = 0; b < 4; b++) U[b] = calx::transforms::fft(indicator(x, b));
  return U;
}

// S(k) = Σ_b |U_b(k)|², k = 0 … N − 1.
inline std::vector<double> spectrum(const Codes& x) {
  const auto U = indicator_transforms(x);
  std::vector<double> S(x.size(), 0.0);
  for (Code b = 0; b < 4; b++)
    for (std::size_t k = 0; k < x.size(); k++) S[k] += std::norm(U[b][k]);
  return S;
}

// ---- the period-3 signal ----------------------------------------------------------------------------

// The counts c_{b,r} of a stretch: base b at the positions ≡ r (mod 3), positions counted from
// any fixed origin (S(N/3) does not depend on which, theory §1.3).
struct CodonCounts {
  std::array<std::array<std::int64_t, 3>, 4> c{};
  std::int64_t bases() const {
    std::int64_t n = 0;
    for (const auto& b : c) n += b[0] + b[1] + b[2];
    return n;
  }
  // S(N/3) = Σ_b N(c_{b,0} + c_{b,1} ζ + c_{b,2} ζ²), N(·) the norm of ℤ[ζ].
  std::int64_t power() const {
    std::int64_t s = 0;
    for (const auto& b : c) s += b[0] * b[0] + b[1] * b[1] + b[2] * b[2] - b[0] * b[1] - b[1] * b[2] - b[2] * b[0];
    return s;
  }
};

inline CodonCounts codon_counts(const Codes& x, std::size_t begin, std::size_t length) {
  if (begin + length > x.size()) throw std::invalid_argument("genomics::signal::codon_counts: a stretch of the sequence");
  CodonCounts k;
  for (std::size_t n = begin; n < begin + length; n++)
    if (x[n] < 4) k.c[x[n]][n % 3]++;
  return k;
}

// S(N/3) of the stretch [begin, begin + length), exactly; length a multiple of 3.
inline std::int64_t period3_power(const Codes& x, std::size_t begin, std::size_t length) {
  if (length % 3) throw std::invalid_argument("genomics::signal::period3_power: a length divisible by 3");
  return codon_counts(x, begin, length).power();
}

// One window of the profile: its first base, S(W/3), and its bases other than N (= S̄).
struct Period3 {
  std::size_t begin = 0;
  std::int64_t power = 0;
  std::int64_t bases = 0;
  double strength() const { return bases > 0 ? static_cast<double>(power) / static_cast<double>(bases) : 0.0; }
};

// The windows [i, i + window) for i = 0, step, 2 step, … while they fit, window a multiple of 3.
// O(n) whatever the step: the counts slide by one base at a time.
inline std::vector<Period3> period3_profile(const Codes& x, std::size_t window, std::size_t step = 1) {
  if (window == 0 || window % 3) throw std::invalid_argument("genomics::signal::period3_profile: a window divisible by 3");
  if (step == 0) throw std::invalid_argument("genomics::signal::period3_profile: step >= 1");
  std::vector<Period3> out;
  if (x.size() < window) return out;
  out.reserve((x.size() - window) / step + 1);
  CodonCounts k = codon_counts(x, 0, window);
  for (std::size_t i = 0;; i++) {
    if (i % step == 0) out.push_back({i, k.power(), k.bases()});
    if (i + window >= x.size()) break;
    if (x[i] < 4) k.c[x[i]][i % 3]--;
    if (x[i + window] < 4) k.c[x[i + window]][(i + window) % 3]++;
  }
  return out;
}

}  // namespace genomics::signal
