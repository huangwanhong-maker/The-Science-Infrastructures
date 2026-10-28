// Checks of the spectra of DNA: the transforms of the indicator sequences by calx's FFT against
// the direct sum over ℤ_N, for a length that is a power of two and one that is not, with bases
// N; Parseval's identity, the mean of the spectrum the number of bases; S(N − k) = S(k) and the
// spectrum of the reverse complement the same; a sequence of period p supported on the multiples
// of N/p, and (ACG)^m peaking exactly at N/3 with S(N/3) = N²/3; the exact period-3 power from
// twelve counts against the transform at N/3, and the sliding profile against counts made afresh
// for every window; and the mean strength of random windows against E P = 1 − Σ_b p_b²
// (theory, Proposition 1.7) at three GC contents.

#include <genomics/signal/spectrum.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <string>
#include <vector>

using namespace genomics;
using signal::Complex;

namespace {

sequence::Codes random_dna(calx::random::Stream& r, std::size_t n, double n_rate) {
  sequence::Codes x(n);
  for (auto& c : x) c = r.uniform() < n_rate ? sequence::N : static_cast<sequence::Code>(r.below(4));
  return x;
}

// U_b(k) = Σ_n u_b(n) e^{−2πi nk/N} directly, the phase nk reduced mod N in integers.
std::vector<Complex> direct(const sequence::Codes& x, sequence::Code b) {
  const std::size_t N = x.size();
  std::vector<Complex> U(N);
  for (std::size_t k = 0; k < N; k++) {
    Complex s = 0;
    for (std::size_t n = 0; n < N; n++)
      if (x[n] == b) {
        const double t = -2 * std::numbers::pi * static_cast<double>(n * k % N) / static_cast<double>(N);
        s += Complex(std::cos(t), std::sin(t));
      }
    U[k] = s;
  }
  return U;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(1992, 0x5350454b, 0, 0);

  // ---- the transforms against the direct sum ---------------------------------------------------------
  for (std::size_t N : {1024, 999}) {
    const sequence::Codes x = random_dna(r, N, 0.02);
    const auto U = signal::indicator_transforms(x);
    const auto S = signal::spectrum(x);
    double eu = 0, es = 0;
    std::vector<double> Sd(N, 0.0);
    for (sequence::Code b = 0; b < 4; b++) {
      const auto D = direct(x, b);
      for (std::size_t k = 0; k < N; k++) {
        eu = std::max(eu, std::abs(U[b][k] - D[k]));
        Sd[k] += std::norm(D[k]);
      }
    }
    for (std::size_t k = 0; k < N; k++) es = std::max(es, std::abs(S[k] - Sd[k]));
    const std::string n = std::to_string(N) + (N == 1024 ? " (radix 2)" : " (Bluestein)");
    s.check("N = " + n + ": U_b(k) by calx's FFT against the direct sum over Z_N, the largest difference over N", eu / static_cast<double>(N), 1e-12);
    s.check("N = " + n + ": S(k) against the direct sum, the largest difference over N²", es / static_cast<double>(N * N), 1e-12);
    // Parseval: Σ_k S(k) = N ν.
    double total = 0;
    for (double v : S) total += v;
    const double nu = static_cast<double>(N - sequence::count_ambiguous(x));
    s.check("N = " + n + ": Parseval, the mean of S over the frequencies equals the bases other than N (relative error)", total / static_cast<double>(N) / nu - 1, 1e-12);
    // Symmetry and the reverse complement.
    double sym = 0, rc = 0;
    const auto Src = signal::spectrum(sequence::reverse_complement(x));
    for (std::size_t k = 1; k < N; k++) sym = std::max(sym, std::abs(S[k] - S[N - k]));
    for (std::size_t k = 0; k < N; k++) rc = std::max(rc, std::abs(S[k] - Src[k]));
    s.check("N = " + n + ": S(N − k) = S(k), the largest difference over N²", sym / static_cast<double>(N * N), 1e-12);
    s.check("N = " + n + ": the reverse complement has the same spectrum, the largest difference over N²", rc / static_cast<double>(N * N), 1e-12);
  }

  // ---- periodic sequences ----------------------------------------------------------------------------
  {
    // (ACG)^m: every indicator has period 3, S(N/3) = 3 m² = N²/3, and S vanishes off 0, N/3, 2N/3.
    const std::size_t m = 400, N = 3 * m;
    sequence::Codes x;
    for (std::size_t i = 0; i < m; i++) x.insert(x.end(), {sequence::A, sequence::C, sequence::G});
    const auto S = signal::spectrum(x);
    std::size_t peak = 1;
    double off = 0;
    for (std::size_t k = 1; k <= N / 2; k++) {
      if (S[k] > S[peak]) peak = k;
      if (k != N / 3) off = std::max(off, S[k]);
    }
    s.require("(ACG)^400: the spectrum peaks at k = N/3 exactly", peak == N / 3, "k = " + std::to_string(peak));
    s.check("(ACG)^400: S(N/3) = N²/3 (relative error)", S[N / 3] / (static_cast<double>(N * N) / 3) - 1, 1e-12);
    s.check("(ACG)^400: S vanishes at every other frequency of (0, N/2], relative to N²", off / static_cast<double>(N * N), 1e-12);
    s.require("(ACG)^400: the exact period-3 power from the counts is N²/3", signal::period3_power(x, 0, N) == static_cast<std::int64_t>(N * N / 3));
    // A random word of period p repeated m times: S(k) = 0 unless m | k.
    bool support = true;
    double worst = 0;
    for (std::size_t p : {5, 7, 10, 21}) {
      const std::size_t reps = 60, L = p * reps;
      const sequence::Codes unit = random_dna(r, p, 0);
      sequence::Codes y;
      for (std::size_t i = 0; i < reps; i++) y.insert(y.end(), unit.begin(), unit.end());
      const auto T = signal::spectrum(y);
      std::size_t top = 1;
      for (std::size_t k = 1; k <= L / 2; k++) {
        if (k % reps) worst = std::max(worst, T[k] / static_cast<double>(L * L));
        if (T[k] > T[top]) top = k;
      }
      support = support && top % reps == 0;
    }
    s.check("words of period p = 5, 7, 10, 21 repeated 60 times: S vanishes off the multiples of N/p, relative to N²", worst, 1e-12);
    s.require("and their highest peak in (0, N/2] lies at a multiple of N/p", support);
  }

  // ---- the period-3 power, exactly ----------------------------------------------------------------
  {
    double worst = 0;
    for (int trial = 0; trial < 50; trial++) {
      const std::size_t L = 3 * (20 + r.below(200));
      const sequence::Codes x = random_dna(r, L + 50, 0.01);
      const std::size_t at = r.below(50);
      const sequence::Codes w(x.begin() + static_cast<std::ptrdiff_t>(at), x.begin() + static_cast<std::ptrdiff_t>(at + L));
      const auto S = signal::spectrum(w);
      const std::int64_t p = signal::period3_power(x, at, L);
      worst = std::max(worst, std::abs(S[L / 3] - static_cast<double>(p)) / static_cast<double>(L * L));
    }
    s.check("on 50 random stretches, the exact period-3 power from the counts (Eisenstein norms) equals S(N/3) by the FFT, relative to N²", worst, 1e-12);
    // The sliding profile against counts made afresh.
    const sequence::Codes x = random_dna(r, 20000, 0.005);
    bool same = true;
    for (std::size_t W : {3, 60, 120, 351})
      for (std::size_t step : {1, 3, 7}) {
        const auto prof = signal::period3_profile(x, W, step);
        same = same && prof.size() == (x.size() - W) / step + 1;
        for (const auto& w : prof) {
          const auto c = signal::codon_counts(x, w.begin, W);
          same = same && w.begin % step == 0 && w.power == c.power() && w.bases == c.bases();
        }
      }
    s.require("the profile sliding by one base equals the counts made afresh for every window (W = 3 to 351, steps 1, 3, 7)", same);
  }

  // ---- random sequence: E P = 1 − Σ p_b² ------------------------------------------------------------
  for (double gc : {0.3, 0.5, 0.7}) {
    const std::size_t W = 120;
    const sequence::Codes x = simulate::random_sequence(3'000'000, gc, 2001, static_cast<std::uint32_t>(gc * 10));
    const auto prof = signal::period3_profile(x, W, W);  // disjoint windows, independent
    double mean = 0, sq = 0;
    for (const auto& w : prof) mean += w.strength();
    mean /= static_cast<double>(prof.size());
    for (const auto& w : prof) sq += (w.strength() - mean) * (w.strength() - mean);
    const double se = std::sqrt(sq / static_cast<double>(prof.size() - 1) / static_cast<double>(prof.size()));
    const double pa = (1 - gc) / 2, pg = gc / 2;
    const double expected = 1 - 2 * pa * pa - 2 * pg * pg;
    s.check("random sequence at GC " + std::to_string(gc).substr(0, 3) + ": the mean strength of " + std::to_string(prof.size()) + " windows of 120 against 1 − Σ p_b² = " + std::to_string(expected).substr(0, 5) + ", in standard errors",
            (mean - expected) / se, 4.0);
  }
  return s.done();
}
