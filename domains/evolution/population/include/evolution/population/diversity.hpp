// Diversity within populations and differentiation between them, from allele counts: the
// segregating sites, the estimators of θ = 4Nμ and the tests built on their differences, the
// site frequency spectra, and F_ST; in windows along a chromosome.
//
// For a sample of n sequences and a site where k carry the derived (or alternate) allele:
//
// - S, the number of segregating sites (0 < k < n), and Watterson's θ_W = S / a_n with
//   a_n = Σ_{i=1}^{n−1} 1/i (Watterson 1975, Theor. Popul. Biol. 7, 256–276): under
//   the neutral coalescent E[S] = θ a_n and Var[S] = θ a_n + θ² b_n, b_n = Σ_{i=1}^{n−1} 1/i².
// - π, the mean number of differences between two sequences (Nei and Li 1979, PNAS 76,
//   5269–5273; Tajima 1983, Genetics 105, 437–460): Σ_sites 2k(n − k) / (n(n − 1)); E[π] = θ.
// - θ_H = Σ_sites 2k² / (n(n − 1)) (Fay and Wu 2000, Genetics 155, 1405–1413), which weighs
//   high-frequency derived alleles and needs the ancestral state (an outgroup); E[θ_H] = θ,
//   and Fay and Wu's H = π − θ_H, unnormalised as they define it.
// - Tajima's D = (π − S/a₁) / √(e₁ S + e₂ S(S − 1)) (Tajima 1989, Genetics 123, 585–595, with
//   the constants as printed there):
//     a₁ = Σ_{i=1}^{n−1} 1/i,  a₂ = Σ 1/i²,  b₁ = (n + 1)/(3(n − 1)),  b₂ = 2(n² + n + 3)/(9n(n − 1)),
//     c₁ = b₁ − 1/a₁,  c₂ = b₂ − (n + 2)/(a₁ n) + a₂/a₁²,  e₁ = c₁/a₁,  e₂ = c₂/(a₁² + a₂);
//   Var[π − S/a₁] = c₁θ + c₂θ², of which e₁S + e₂S(S − 1) is the unbiased estimate.
// - The unfolded spectrum ξ_i (sites with i derived copies, i = 0 … n) and the folded one
//   η_i = ξ_i + ξ_{n−i} (i ≤ n/2, without double counting i = n/2); under the neutral
//   coalescent E[ξ_i] = θ/i (Fu 1995, Theor. Popul. Biol. 48, 172–197).
// - Missing calls make n vary from site to site. θ_W and π then sum each site's own terms
//   (1/a_{n_j} and 2k(n_j − k)/(n_j(n_j − 1))); the spectra, and the D and H read from them, are
//   projected to a common n by hypergeometric sampling without replacement (Marth et al. 2004,
//   Genetics 166, 351–372; Gutenkunst et al. 2009, PLoS Genetics 5, e1000695), which is the
//   identity where nothing is missing.
// - F_ST between two populations by Hudson, Slatkin and Maddison (1992, Genetics 132, 583–589)
//   as Bhatia, Patterson, Sankararaman and Price (2013, Genome Research 23, 1514–1521) write it
//   for samples: at each site
//     N = (p₁ − p₂)² − p₁q₁/(n₁ − 1) − p₂q₂/(n₂ − 1),   D = p₁q₂ + p₂q₁,
//   and F_ST = Σ N / Σ D over sites, a ratio of averages (their recommendation), never an
//   average of ratios.
// - θ of Weir and Cockerham (1984, Evolution 38, 1358–1370) for r ≥ 2
//   populations of n_i individuals with allele frequencies p̃_i and heterozygote frequencies h̃_i:
//     n̄ = Σ n_i / r,  n_c = (r n̄ − Σ n_i²/(r n̄))/(r − 1),  p̄ = Σ n_i p̃_i/(r n̄),
//     s² = Σ n_i (p̃_i − p̄)²/((r − 1) n̄),  h̄ = Σ n_i h̃_i/(r n̄),
//     a = n̄/n_c [s² − (p̄(1 − p̄) − (r − 1)s²/r − h̄/4)/(n̄ − 1)],
//     b = n̄/(n̄ − 1) [p̄(1 − p̄) − (r − 1)s²/r − (2n̄ − 1) h̄/(4n̄)],   c = h̄/2,
//   and θ̂ = Σ a / Σ (a + b + c) over sites. It needs genotypes (h̃), not only allele counts.
//
// Complexity: one pass over the sites, O(m); a projection costs O(n) per site whose sample
// size differs from the target. Estimates are sums over sites, so windows are cheap: a window
// is a range of sites, found by binary search on the positions.

#pragma once

#include <evolution/population/genotypes.hpp>

#include <calx/special_functions/gamma.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace evolution::population {

inline constexpr double NOT_A_NUMBER = std::numeric_limits<double>::quiet_NaN();

// a_n = Σ_{i=1}^{n−1} 1/i and b_n = Σ_{i=1}^{n−1} 1/i².
inline double watterson_a(std::size_t n) {
  double s = 0;
  for (std::size_t i = n; i-- > 1;) s += 1.0 / static_cast<double>(i);
  return s;
}
inline double watterson_b(std::size_t n) {
  double s = 0;
  for (std::size_t i = n; i-- > 1;) s += 1.0 / (static_cast<double>(i) * static_cast<double>(i));
  return s;
}

inline std::size_t segregating_sites(CountsView c) {
  std::size_t s = 0;
  for (std::size_t j = 0; j < c.sites(); j++) s += c.alt[j] > 0 && c.alt[j] < c.total[j];
  return s;
}

// θ_W = Σ over segregating sites of 1/a_{n_j}: S/a_n when no call is missing.
inline double theta_watterson(CountsView c) {
  double s = 0;
  std::uint32_t last_n = 0;
  double last_a = 0;
  for (std::size_t j = 0; j < c.sites(); j++) {
    if (c.alt[j] == 0 || c.alt[j] >= c.total[j]) continue;
    if (c.total[j] != last_n) {
      last_n = c.total[j];
      last_a = watterson_a(last_n);
    }
    s += 1 / last_a;
  }
  return s;
}

// π: the expected number of differences between two sequences drawn without replacement,
// summed over sites (divide by the length for π per base).
inline double nucleotide_diversity(CountsView c) {
  double s = 0;
  for (std::size_t j = 0; j < c.sites(); j++) {
    const double n = c.total[j], k = c.alt[j];
    if (n >= 2) s += 2 * k * (n - k) / (n * (n - 1));
  }
  return s;
}

// θ_H, from derived-allele counts.
inline double theta_h(CountsView derived) {
  double s = 0;
  for (std::size_t j = 0; j < derived.sites(); j++) {
    const double n = derived.total[j], k = derived.alt[j];
    if (n >= 2 && k < n) s += 2 * k * k / (n * (n - 1));
  }
  return s;
}

struct TajimaConstants {
  double a1 = 0, a2 = 0, b1 = 0, b2 = 0, c1 = 0, c2 = 0, e1 = 0, e2 = 0;
};

inline TajimaConstants tajima_constants(std::size_t n) {
  if (n < 2) throw std::invalid_argument("evolution::population::tajima_constants: n >= 2");
  TajimaConstants t;
  const double x = static_cast<double>(n);
  t.a1 = watterson_a(n);
  t.a2 = watterson_b(n);
  t.b1 = (x + 1) / (3 * (x - 1));
  t.b2 = 2 * (x * x + x + 3) / (9 * x * (x - 1));
  t.c1 = t.b1 - 1 / t.a1;
  t.c2 = t.b2 - (x + 2) / (t.a1 * x) + t.a2 / (t.a1 * t.a1);
  t.e1 = t.c1 / t.a1;
  t.e2 = t.c2 / (t.a1 * t.a1 + t.a2);
  return t;
}

// Tajima's D from n, S and π (NaN when S = 0).
inline double tajima_d(std::size_t n, double S, double pi) {
  if (!(S > 0)) return NOT_A_NUMBER;
  const TajimaConstants t = tajima_constants(n);
  return (pi - S / t.a1) / std::sqrt(t.e1 * S + t.e2 * S * (S - 1));
}

// ---- spectra --------------------------------------------------------------------------------------

namespace detail {

// Adds to xi the law of the derived count in a subsample of n drawn without replacement from N
// sequences of which k are derived: the hypergeometric law, computed from its mode by the
// ratio P(i + 1)/P(i) = (k − i)(n − i) / ((i + 1)(N − k − n + i + 1)).
inline void add_hypergeometric(std::vector<double>& xi, std::uint32_t N, std::uint32_t k, std::uint32_t n) {
  if (N == n) {
    xi[k] += 1;
    return;
  }
  const std::int64_t lo = std::max<std::int64_t>(0, static_cast<std::int64_t>(n) - (static_cast<std::int64_t>(N) - k));
  const std::int64_t hi = std::min<std::int64_t>(k, n);
  std::int64_t mode = static_cast<std::int64_t>((static_cast<double>(n) + 1) * (static_cast<double>(k) + 1) / (static_cast<double>(N) + 2));
  mode = std::clamp(mode, lo, hi);
  using calx::special_functions::log_choose;
  const double Nd = N, kd = k, nd = n;
  const double at_mode =
      std::exp(log_choose(kd, static_cast<double>(mode)) + log_choose(Nd - kd, nd - static_cast<double>(mode)) - log_choose(Nd, nd));
  xi[static_cast<std::size_t>(mode)] += at_mode;
  double p = at_mode;
  for (std::int64_t i = mode; i < hi && p > 0; i++) {
    const double id = static_cast<double>(i);
    p *= (kd - id) * (nd - id) / ((id + 1) * (Nd - kd - nd + id + 1));
    xi[static_cast<std::size_t>(i + 1)] += p;
  }
  p = at_mode;
  for (std::int64_t i = mode; i > lo && p > 0; i--) {
    const double id = static_cast<double>(i);
    p *= id * (Nd - kd - nd + id) / ((kd - id + 1) * (nd - id + 1));
    xi[static_cast<std::size_t>(i - 1)] += p;
  }
}

}  // namespace detail

// The unfolded spectrum ξ_0 … ξ_n of derived counts, every site projected to n sequences
// (sites with fewer than n typed are left out). n = 0 takes the smallest sample size of a site.
inline std::vector<double> site_frequency_spectrum(CountsView derived, std::size_t n = 0) {
  if (n == 0) {
    std::uint32_t least = std::numeric_limits<std::uint32_t>::max();
    for (std::size_t j = 0; j < derived.sites(); j++) least = std::min(least, derived.total[j]);
    n = derived.sites() ? least : 0;
  }
  std::vector<double> xi(n + 1, 0.0);
  for (std::size_t j = 0; j < derived.sites(); j++) {
    if (derived.total[j] < n || derived.alt[j] > derived.total[j]) continue;
    detail::add_hypergeometric(xi, derived.total[j], derived.alt[j], static_cast<std::uint32_t>(n));
  }
  return xi;
}

// The folded spectrum η_0 … η_{⌊n/2⌋} of an unfolded one.
inline std::vector<double> fold(std::span<const double> xi) {
  if (xi.empty()) return {};
  const std::size_t n = xi.size() - 1;
  std::vector<double> eta(n / 2 + 1, 0.0);
  for (std::size_t i = 0; i <= n; i++) eta[std::min(i, n - i)] += xi[i];
  return eta;
}

// S, π and θ_H of an unfolded spectrum of n = xi.size() − 1.
struct SpectrumSummary {
  double S = 0, pi = 0, theta_h = 0;
};
inline SpectrumSummary summarise(std::span<const double> xi) {
  SpectrumSummary s;
  if (xi.size() < 3) return s;
  const double n = static_cast<double>(xi.size() - 1);
  for (std::size_t i = 1; i + 1 < xi.size(); i++) {
    const double k = static_cast<double>(i);
    s.S += xi[i];
    s.pi += 2 * k * (n - k) / (n * (n - 1)) * xi[i];
    s.theta_h += 2 * k * k / (n * (n - 1)) * xi[i];
  }
  return s;
}

// Tajima's D of a set of sites, from its spectrum projected to the smallest sample size.
inline double tajima_d(CountsView c) {
  const std::vector<double> xi = site_frequency_spectrum(c);
  if (xi.size() < 3) return NOT_A_NUMBER;
  const SpectrumSummary s = summarise(xi);
  return tajima_d(xi.size() - 1, s.S, s.pi);
}

// Fay and Wu's H = π − θ_H of derived counts, from the projected spectrum.
inline double fay_wu_h(CountsView derived) {
  const std::vector<double> xi = site_frequency_spectrum(derived);
  if (xi.size() < 3) return NOT_A_NUMBER;
  const SpectrumSummary s = summarise(xi);
  return s.pi - s.theta_h;
}

// Derived-allele counts from alternate-allele counts and the outgroup's allele at each site
// (0: the reference is ancestral; 1: the alternate is; −1: unknown, and the site is dropped by
// giving it no observations).
inline AlleleCounts polarise(CountsView c, std::span<const std::int8_t> outgroup) {
  if (outgroup.size() != c.sites()) throw std::invalid_argument("evolution::population::polarise: one outgroup allele for each site");
  AlleleCounts d;
  d.alt.resize(c.sites());
  d.total.resize(c.sites());
  for (std::size_t j = 0; j < c.sites(); j++) {
    if (outgroup[j] == 0) {
      d.alt[j] = c.alt[j];
      d.total[j] = c.total[j];
    } else if (outgroup[j] == 1) {
      d.alt[j] = c.total[j] - c.alt[j];
      d.total[j] = c.total[j];
    } else {
      d.alt[j] = d.total[j] = 0;
    }
  }
  return d;
}

// ---- F_ST ---------------------------------------------------------------------------------------------

// A ratio of sums over sites.
struct Ratio {
  double numerator = 0, denominator = 0;
  double value() const { return denominator > 0 ? numerator / denominator : NOT_A_NUMBER; }
  Ratio& operator+=(const Ratio& o) {
    numerator += o.numerator;
    denominator += o.denominator;
    return *this;
  }
};

// Hudson's N and D at one site, from the alternate counts and sample sizes of two populations;
// zero where either has fewer than two sequences.
inline Ratio hudson_site(std::uint32_t alt1, std::uint32_t n1, std::uint32_t alt2, std::uint32_t n2) {
  if (n1 < 2 || n2 < 2) return {};
  const double p1 = static_cast<double>(alt1) / n1, p2 = static_cast<double>(alt2) / n2;
  const double num = (p1 - p2) * (p1 - p2) - p1 * (1 - p1) / (n1 - 1.0) - p2 * (1 - p2) / (n2 - 1.0);
  return {num, p1 * (1 - p2) + p2 * (1 - p1)};
}

inline Ratio fst_hudson(CountsView a, CountsView b) {
  if (a.sites() != b.sites()) throw std::invalid_argument("evolution::population::fst_hudson: the same sites in both populations");
  Ratio r;
  for (std::size_t j = 0; j < a.sites(); j++) r += hudson_site(a.alt[j], a.total[j], b.alt[j], b.total[j]);
  return r;
}

// Weir and Cockerham's variance components at one site: between populations (a), between
// individuals within them (b), and between the gametes of an individual (c).
struct VarianceComponents {
  double a = 0, b = 0, c = 0;
};

// From each population's individuals typed, copies of the alternate allele, and heterozygotes.
// Populations with nobody typed are left out; with fewer than two left, or n̄ ≤ 1, all are 0.
inline VarianceComponents weir_cockerham_site(std::span<const std::uint32_t> individuals, std::span<const std::uint32_t> alt,
                                              std::span<const std::uint32_t> het) {
  double sum_n = 0, sum_n2 = 0, r = 0;
  for (std::size_t i = 0; i < individuals.size(); i++)
    if (individuals[i]) {
      r += 1;
      sum_n += individuals[i];
      sum_n2 += static_cast<double>(individuals[i]) * individuals[i];
    }
  if (r < 2) return {};
  const double n_bar = sum_n / r;
  if (!(n_bar > 1)) return {};
  const double n_c = (r * n_bar - sum_n2 / (r * n_bar)) / (r - 1);
  double p_bar = 0, h_bar = 0;
  for (std::size_t i = 0; i < individuals.size(); i++)
    if (individuals[i]) {
      p_bar += alt[i] / 2.0;  // n_i p̃_i
      h_bar += het[i];        // n_i h̃_i
    }
  p_bar /= r * n_bar;
  h_bar /= r * n_bar;
  double s2 = 0;
  for (std::size_t i = 0; i < individuals.size(); i++)
    if (individuals[i]) {
      const double p = alt[i] / (2.0 * individuals[i]);
      s2 += individuals[i] * (p - p_bar) * (p - p_bar);
    }
  s2 /= (r - 1) * n_bar;
  const double pq = p_bar * (1 - p_bar);
  VarianceComponents v;
  v.a = n_bar / n_c * (s2 - (pq - (r - 1) / r * s2 - h_bar / 4) / (n_bar - 1));
  v.b = n_bar / (n_bar - 1) * (pq - (r - 1) / r * s2 - (2 * n_bar - 1) / (4 * n_bar) * h_bar);
  v.c = h_bar / 2;
  return v;
}

// θ̂ of Weir and Cockerham over the sites [first, last) (all by default), the samples labelled
// by population: Σ a / Σ (a + b + c).
inline Ratio fst_weir_cockerham(const Genotypes& g, std::span<const std::uint32_t> population, std::size_t first = 0,
                                std::size_t last = static_cast<std::size_t>(-1)) {
  if (population.size() != g.samples) throw std::invalid_argument("evolution::population::fst_weir_cockerham: one population for each sample");
  last = std::min(last, g.sites);
  std::uint32_t r = 0;
  for (const auto p : population) r = std::max(r, p + 1);
  std::vector<std::uint32_t> n(r), alt(r), het(r);
  Ratio total;
  for (std::size_t j = first; j < last; j++) {
    std::fill(n.begin(), n.end(), 0);
    std::fill(alt.begin(), alt.end(), 0);
    std::fill(het.begin(), het.end(), 0);
    const auto row = g.site(j);
    for (std::size_t i = 0; i < g.samples; i++) {
      const std::int8_t x = row[i];
      if (x < 0) continue;
      n[population[i]]++;
      alt[population[i]] += static_cast<std::uint32_t>(x);
      het[population[i]] += x == 1;
    }
    const VarianceComponents v = weir_cockerham_site(n, alt, het);
    total += Ratio{v.a, v.a + v.b + v.c};
  }
  return total;
}

// ---- windows --------------------------------------------------------------------------------------

// A window [begin, end) of a chromosome, holding the sites [first, last).
struct Window {
  std::uint32_t chromosome = 0;
  std::uint64_t begin = 0, end = 0;
  std::size_t first = 0, last = 0;
};

// Windows of `size` bases every `step` bases along each chromosome, from 0 to its last site.
inline std::vector<Window> windows(const Genotypes& g, std::uint64_t size, std::uint64_t step) {
  if (size == 0 || step == 0) throw std::invalid_argument("evolution::population::windows: a positive size and step");
  std::vector<Window> out;
  for (std::uint32_t c = 0; c < g.chromosomes.size(); c++) {
    const auto [lo, hi] = g.chromosome_sites(c);
    if (lo == hi) continue;
    const std::uint64_t end_of = g.positions[hi - 1] + 1;
    const auto first = g.positions.begin() + static_cast<std::ptrdiff_t>(lo), last = g.positions.begin() + static_cast<std::ptrdiff_t>(hi);
    for (std::uint64_t b = 0; b < end_of; b += step) {
      Window w;
      w.chromosome = c;
      w.begin = b;
      w.end = b + size;
      w.first = static_cast<std::size_t>(std::lower_bound(first, last, w.begin) - g.positions.begin());
      w.last = static_cast<std::size_t>(std::lower_bound(first, last, w.end) - g.positions.begin());
      out.push_back(w);
      if (b + size >= end_of) break;
    }
  }
  return out;
}

struct WindowStatistics {
  Window window;
  std::size_t sites = 0, segregating = 0;
  double pi = 0, theta_w = 0;  // per base
  double tajima_d = NOT_A_NUMBER;
  double fst_hudson = NOT_A_NUMBER, fst_weir_cockerham = NOT_A_NUMBER;
};

// Diversity in windows over all samples; with population labels, F_ST too (Hudson's between
// populations 0 and 1, Weir and Cockerham's among all).
inline std::vector<WindowStatistics> window_statistics(const Genotypes& g, std::uint64_t size, std::uint64_t step,
                                                       std::span<const std::uint32_t> population = {}) {
  const AlleleCounts all = allele_counts(g);
  AlleleCounts first_pop, second_pop;
  if (!population.empty()) {
    const auto m = members(population);
    if (m.size() >= 2) {
      first_pop = allele_counts(g, m[0]);
      second_pop = allele_counts(g, m[1]);
    }
  }
  std::vector<WindowStatistics> out;
  for (const Window& w : windows(g, size, step)) {
    WindowStatistics s;
    s.window = w;
    const CountsView c = all.slice(w.first, w.last);
    const double length = static_cast<double>(w.end - w.begin);
    s.sites = w.last - w.first;
    s.segregating = segregating_sites(c);
    s.pi = nucleotide_diversity(c) / length;
    s.theta_w = theta_watterson(c) / length;
    s.tajima_d = tajima_d(c);
    if (!first_pop.alt.empty()) {
      s.fst_hudson = fst_hudson(first_pop.slice(w.first, w.last), second_pop.slice(w.first, w.last)).value();
      s.fst_weir_cockerham = fst_weir_cockerham(g, population, w.first, w.last).value();
    }
    out.push_back(s);
  }
  return out;
}

}  // namespace evolution::population
