// The effective terrain in time (theory/street-oceanography.md §2.4; Huang, Street Oceanography, Claim 22.1
// and §22.4): what the declared hours of a street's holdings say about the street over the week.
//
// The opening forcing F(θ) = Σ_h 1_{O_h}(θ) counts the holdings open at hour θ of the week. Its Fourier
// lines F̂(m) = (1/168) ∫_W F(θ) e^{−2πimθ/168} dθ are exact for a sum of intervals:
// ∫_a^b e^{−iνθ} dθ = (e^{−iνa} − e^{−iνb}) / (iν), ν = 2πm/168 (a wrapped interval integrates as it
// stands, the exponential having period 168). Lines with m a multiple of 7 are the day's harmonics; the
// others carry the week's modulation of the days (Definition 2.12).
//
// The open fraction ω(s, θ) is the kernel-weighted share of holdings near s open at θ, with a Gaussian
// kernel of declared width.
//
// Ascent and descent (Definition 2.13). For a weekday, the opening and closing times of the holdings open
// that day: their dispersions σ_o and σ_c, the asymmetry A = (σ_c − σ_o)/(σ_c + σ_o), and Kendall's τ_b
// between the two orders (Kendall 1945, Biometrika 33, 239, with ties), with the permutation null that
// exchanges the closing times among the holdings (Monte Carlo, calx's seeded stream; p two-sided,
// (1 + #{|τ*| ≥ |τ|}) / (1 + permutations)).

#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <vector>

#include <calx/random/stream.hpp>
#include <urban/street/street.hpp>

namespace urban::street {

// F(θ) at the centres of n equal steps of the week (θ_k = (k + ½)·168/n), over the holdings whose use is in
// `uses` (empty: every use). Holdings without declared hours are left out.
inline std::vector<double> opening_forcing(const Street& st, std::size_t n = 168, const std::vector<std::size_t>& uses = {}) {
  std::vector<double> f(n, 0);
  for (const auto& h : st.holdings) {
    if (!h.declared) continue;
    if (!uses.empty() && std::find(uses.begin(), uses.end(), h.use) == uses.end()) continue;
    for (std::size_t k = 0; k < n; k++)
      if (open_at(h, (static_cast<double>(k) + 0.5) * week_hours / static_cast<double>(n))) f[k] += 1;
  }
  return f;
}

struct WeekLine {
  std::size_t m = 0;              // cycles a week
  double frequency = 0;           // cycles per day, m/7
  std::complex<double> value;     // F̂(m)
  double amplitude = 0;           // 2|F̂(m)| for m > 0 (the cosine's amplitude), |F̂(0)| for the mean
  bool daily = false;             // m a multiple of 7: a harmonic of the civil day
};

// The exact lines m = 0…m_max of the opening forcing.
inline std::vector<WeekLine> week_lines(const Street& st, std::size_t m_max = 84, const std::vector<std::size_t>& uses = {}) {
  std::vector<WeekLine> out;
  for (std::size_t m = 0; m <= m_max; m++) {
    std::complex<double> sum{0, 0};
    const double nu = 2 * std::numbers::pi * static_cast<double>(m) / week_hours;
    for (const auto& h : st.holdings) {
      if (!h.declared) continue;
      if (!uses.empty() && std::find(uses.begin(), uses.end(), h.use) == uses.end()) continue;
      for (const auto& [a, b] : h.hours) {
        if (m == 0) sum += b - a;
        else sum += (std::exp(std::complex<double>(0, -nu * a)) - std::exp(std::complex<double>(0, -nu * b))) / std::complex<double>(0, nu);
      }
    }
    WeekLine l;
    l.m = m, l.frequency = static_cast<double>(m) / 7.0, l.value = sum / week_hours;
    l.amplitude = (m == 0 ? 1.0 : 2.0) * std::abs(l.value);
    l.daily = m > 0 && m % 7 == 0;
    out.push_back(l);
  }
  return out;
}

// ω(s, θ) on a grid of positions (rows) and hours (columns); NaN where no holding has weight.
inline std::vector<std::vector<double>> open_fraction(const Street& st, const std::vector<double>& s, const std::vector<double>& theta, double width) {
  std::vector<std::vector<double>> out(s.size(), std::vector<double>(theta.size(), std::nan("")));
  for (std::size_t i = 0; i < s.size(); i++) {
    std::vector<double> weight;
    std::vector<const Holding*> near;
    for (const auto& h : st.holdings) {
      if (!h.declared) continue;
      const double z = (s[i] - h.s) / width;
      if (std::fabs(z) > 4) continue;
      weight.push_back(std::exp(-0.5 * z * z));
      near.push_back(&h);
    }
    double total = 0;
    for (double w : weight) total += w;
    if (total <= 0) continue;
    for (std::size_t k = 0; k < theta.size(); k++) {
      double open = 0;
      for (std::size_t j = 0; j < near.size(); j++)
        if (open_at(*near[j], theta[k])) open += weight[j];
      out[i][k] = open / total;
    }
  }
  return out;
}

// ---- ascent and descent ---------------------------------------------------------------------------------

// Kendall's τ_b between two lists (ties in either excluded pairwise and corrected in the denominator).
inline double kendall_tau(const std::vector<double>& a, const std::vector<double>& b) {
  long long concordant = 0, discordant = 0, ties_a = 0, ties_b = 0;
  for (std::size_t i = 0; i < a.size(); i++)
    for (std::size_t j = i + 1; j < a.size(); j++) {
      const double da = a[i] - a[j], db = b[i] - b[j];
      if (da == 0 && db == 0) continue;
      if (da == 0) { ties_a++; continue; }
      if (db == 0) { ties_b++; continue; }
      (da > 0) == (db > 0) ? concordant++ : discordant++;
    }
  const double n1 = static_cast<double>(concordant + discordant + ties_a), n2 = static_cast<double>(concordant + discordant + ties_b);
  if (n1 <= 0 || n2 <= 0) return 0;
  return static_cast<double>(concordant - discordant) / std::sqrt(n1 * n2);
}

struct AscentDescent {
  int day = 0;                     // 0 = Monday
  std::size_t holdings = 0;        // open that day with one opening and one closing
  std::vector<double> open, close; // hours of the day (close may exceed 24)
  double mean_open = 0, mean_close = 0, sd_open = 0, sd_close = 0;
  double asymmetry = 0;            // (σ_c − σ_o)/(σ_c + σ_o)
  double tau = 0;                  // Kendall's τ_b of the opening and closing orders
  double p = 1;                    // permutation p, two-sided
  std::size_t permutations = 0;
};

// Each holding's day on a weekday: its first opening and its last closing among the intervals that start that
// day (a restaurant open 11:30–14:30 and 17:30–22:00 opens at 11:30 and closes at 22:00).
inline AscentDescent ascent_descent(const Street& st, int day, std::size_t permutations = 999, std::uint64_t seed = 1) {
  AscentDescent out;
  out.day = day;
  for (const auto& h : st.holdings) {
    if (!h.declared) continue;
    double first = 1e9, close = 0;
    for (const auto& [a, b] : h.hours)
      if (a >= 24.0 * day && a < 24.0 * (day + 1)) first = std::min(first, a), close = std::max(close, b);
    if (first > 1e8) continue;
    out.open.push_back(first - 24.0 * day);
    out.close.push_back(close - 24.0 * day);
  }
  out.holdings = out.open.size();
  if (out.holdings < 3) return out;
  auto moments = [](const std::vector<double>& v, double& mean, double& sd) {
    mean = 0;
    for (double x : v) mean += x;
    mean /= static_cast<double>(v.size());
    double ss = 0;
    for (double x : v) ss += (x - mean) * (x - mean);
    sd = std::sqrt(ss / static_cast<double>(v.size() - 1));
  };
  moments(out.open, out.mean_open, out.sd_open);
  moments(out.close, out.mean_close, out.sd_close);
  out.asymmetry = out.sd_open + out.sd_close > 0 ? (out.sd_close - out.sd_open) / (out.sd_close + out.sd_open) : 0;
  out.tau = kendall_tau(out.open, out.close);
  calx::random::Stream rng(seed, 0x0ad, static_cast<std::uint32_t>(day), 0);
  std::vector<double> shuffled = out.close;
  std::size_t extreme = 0;
  for (std::size_t p = 0; p < permutations; p++) {
    for (std::size_t i = shuffled.size() - 1; i > 0; i--) std::swap(shuffled[i], shuffled[rng.below(i + 1)]);
    if (std::fabs(kendall_tau(out.open, shuffled)) >= std::fabs(out.tau) - 1e-12) extreme++;
  }
  out.permutations = permutations;
  out.p = (1.0 + static_cast<double>(extreme)) / (1.0 + static_cast<double>(permutations));
  return out;
}

}  // namespace urban::street
