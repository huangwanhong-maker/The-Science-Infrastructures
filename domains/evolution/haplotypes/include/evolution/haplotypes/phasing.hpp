// Unphased diploid genotypes as a pair of copying processes over a reference panel: the
// likelihood, the posterior of the pair of haplotypes copied, and phasing as the jointly most
// probable pair of copying paths and haplotypes, each in O(N K²) for N sites and K panel
// haplotypes.
//
// The model. A diploid's two haplotypes x and y each copy the panel as in copying.hpp,
// independently and with the same chain: the hidden state at site j is the ordered pair
// (A_j, B_j) ∈ {1, …, K}², x copying h_{A_j} and y copying h_{B_j}. The genotype g_j = x_j + y_j
// is observed, so the emission sums over the alleles that give it:
//
//   P(g_j | a, b) = Σ_{x + y = g_j} P(x | h_a(j)) P(y | h_b(j)),   P(x | h) = 1 − ε if x = h, else ε,
//
// and a missing genotype (−1) is emitted with probability 1. A hidden Markov model over pairs
// of panel haplotypes is how the first version of IMPUTE imputed unphased genotypes (Marchini,
// Howie, Myers, McVean and Donnelly 2007, Nature Genetics 39, 906–913); here each copy follows
// Li and Stephens's chain exactly as in copying.hpp.
//
// The step. With u = 1 − r_j and q = r_j/K, the pair's transition is a product,
// P((a, b) → (c, d)) = (u[a = c] + q)(u[b = d] + q), so the sum over the K² predecessors is
//
//   Σ_{a,b} α(a, b) P((a, b) → (c, d)) = u² α(c, d) + u q (Σ_b α(c, b) + Σ_a α(a, d)) + q² Σ_{a,b} α(a, b),
//
// four terms read from the row sums, the column sums and the total: O(K²) per site, not O(K⁴);
// the backward pass is the same identity transposed.
//
// Phasing. The phase of x and y is read from the joint maximum a posteriori of the copying
// paths and the haplotypes, argmax over (A, B, x, y) with x + y = g of P(A, B, x, y): the
// Viterbi recursion over pairs with each emission replaced by its largest term,
// max_{x + y = g} P(x | h_a) P(y | h_b), and the alleles chosen with it. The maximum over the
// predecessors is again four terms, with s = u + q the probability of staying:
//
//   max{ s² δ(c, d), s q max_b δ(c, b), s q max_a δ(a, d), q² max_{a,b} δ(a, b) },
//
// exact because q ≤ s: a term that includes a predecessor with a smaller weight than its own
// is never larger than that predecessor's own term. Ties go to the first of the four, and a
// heterozygous site whose two copied haplotypes carry the same allele gives allele 1 to x. At a
// missing genotype the alleles are those copied (x = h_a, y = h_b), which multiplies every pair's
// emission by the same (1 − ε)² and so leaves the path unchanged. The model is symmetric in x
// and y, so the marginal posterior that x carries allele 1 at a heterozygous site is ½ whatever
// the data: phase is a property of a path, which is why phasing takes the joint maximum.
//
// Switch errors, the usual measure of phasing: among the sites heterozygous in both the truth
// and the inference, a switch is a pair of consecutive such sites whose relative phase differs;
// the rate is switches / (sites − 1). Random phase gives ½.
//
// Complexity: O(N K²) time for the likelihood, posterior and phasing; O(K²) memory for the
// likelihood, O(N K²) doubles for the stored forward pass of the posterior and O(N K²) bytes for
// phasing's back-pointers. pair_log_likelihood_dense costs O(N K⁴) and is the reference. The
// pair model is for small panels (K up to a few hundred); larger ones are phased by iterating
// haploid copying (SHAPEIT, Beagle), not written here.

#pragma once

#include <evolution/haplotypes/copying.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace evolution::haplotypes {

namespace detail {

inline void check_genotypes(const Chain& c, std::span<const std::int8_t> g, const char* where) {
  if (g.size() != c.sites) fail(where, "one genotype for each site of the panel");
  for (const std::int8_t a : g)
    if (a < -1 || a > 2) fail(where, "a genotype is 0, 1, 2 or −1 (missing)");
}

// P(allele 1 | copying h_k(j)) for every k.
inline void allele_one(const Haplotypes& panel, std::size_t j, double eps, double* out) {
  const std::int8_t* h = panel.alleles.data() + j * panel.sequences;
  for (std::size_t k = 0; k < panel.sequences; k++) out[k] = h[k] ? 1 - eps : eps;
}

// P(g | a, b) from p = P(1 | a) and r = P(1 | b); its largest term when `largest`.
inline double genotype_emission(std::int8_t g, double p, double r, bool largest = false) {
  switch (g) {
    case 0: return (1 - p) * (1 - r);
    case 2: return p * r;
    case 1: return largest ? std::max(p * (1 - r), (1 - p) * r) : p * (1 - r) + (1 - p) * r;
    default: return 1;
  }
}

// The K × K emissions of site j into e (row a, column b).
inline void pair_emissions(const Haplotypes& panel, std::size_t j, std::int8_t g, double eps, std::vector<double>& p1, double* e) {
  const std::size_t K = panel.sequences;
  if (g < 0) {
    std::fill(e, e + K * K, 1.0);
    return;
  }
  allele_one(panel, j, eps, p1.data());
  for (std::size_t a = 0; a < K; a++)
    for (std::size_t b = 0; b < K; b++) e[a * K + b] = genotype_emission(g, p1[a], p1[b]);
}

// One forward pass over pairs, calling row(j, α̂_j) on each normalised row; returns log P(g).
template <class Row>
double pair_forward_pass(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> g, Row&& row) {
  const std::size_t K = c.states, N = c.sites, KK = K * K;
  std::vector<double> a(KK), e(KK), p1(K), rows(K), cols(K);
  double logp = 0;
  for (std::size_t j = 0; j < N; j++) {
    pair_emissions(panel, j, g[j], c.epsilon, p1, e.data());
    double total = 0;
    if (j == 0) {
      const double u = 1 / static_cast<double>(KK);
      for (std::size_t i = 0; i < KK; i++) total += a[i] = u * e[i];
    } else {
      std::fill(rows.begin(), rows.end(), 0.0);
      std::fill(cols.begin(), cols.end(), 0.0);
      for (std::size_t x = 0; x < K; x++)
        for (std::size_t y = 0; y < K; y++) {
          rows[x] += a[x * K + y];
          cols[y] += a[x * K + y];
        }
      const double u = 1 - c.jump[j], q = c.move(j), uu = u * u, uq = u * q, qq = q * q;
      for (std::size_t x = 0; x < K; x++)
        for (std::size_t y = 0; y < K; y++) {
          const std::size_t i = x * K + y;
          total += a[i] = e[i] * (uu * a[i] + uq * (rows[x] + cols[y]) + qq);
        }
    }
    if (!(total > 0)) return NEG_INF;
    const double inv = 1 / total;
    for (double& v : a) v *= inv;
    logp += std::log(total);
    row(j, static_cast<const std::vector<double>&>(a));
  }
  return logp;
}

}  // namespace detail

// log P(g) under the pair model, in O(K²) memory.
inline double pair_log_likelihood(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> g) {
  detail::check_genotypes(c, g, "pair_log_likelihood");
  return detail::pair_forward_pass(panel, c, g, [](std::size_t, const std::vector<double>&) {});
}

// The same by the general recursion over the K² × K² transitions: O(N K⁴), the reference.
inline double pair_log_likelihood_dense(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> g) {
  detail::check_genotypes(c, g, "pair_log_likelihood_dense");
  const std::size_t K = c.states, KK = K * K;
  std::vector<double> a(KK), next(KK), e(KK), p1(K);
  double logp = 0;
  for (std::size_t j = 0; j < c.sites; j++) {
    detail::pair_emissions(panel, j, g[j], c.epsilon, p1, e.data());
    double total = 0;
    for (std::size_t cd = 0; cd < KK; cd++) {
      double s = 0;
      if (j == 0) s = 1 / static_cast<double>(KK);
      else
        for (std::size_t ab = 0; ab < KK; ab++) s += a[ab] * c.transition(j, ab / K, cd / K) * c.transition(j, ab % K, cd % K);
      total += next[cd] = s * e[cd];
    }
    if (!(total > 0)) return NEG_INF;
    for (std::size_t i = 0; i < KK; i++) a[i] = next[i] / total;
    logp += std::log(total);
  }
  return logp;
}

// Calls visit(j, γ_j) for j = N − 1 down to 0, γ_j(a, b) = P(A_j = a, B_j = b | g) as K × K
// (row a, column b); returns log P(g), −∞ (and no visits) when g cannot be emitted.
template <class Visit>
double for_each_pair_posterior(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> g, Visit&& visit) {
  detail::check_genotypes(c, g, "for_each_pair_posterior");
  const std::size_t K = c.states, N = c.sites, KK = K * K;
  std::vector<double> alpha(N * KK);
  const double logp = detail::pair_forward_pass(panel, c, g, [&](std::size_t j, const std::vector<double>& a) {
    std::copy(a.begin(), a.end(), alpha.begin() + static_cast<std::ptrdiff_t>(j * KK));
  });
  if (logp == NEG_INF) return logp;
  std::vector<double> beta(KK, 1 / static_cast<double>(KK)), w(KK), e(KK), p1(K), gamma(KK), rows(K), cols(K);
  for (std::size_t j = N; j-- > 0;) {
    const double* a = &alpha[j * KK];
    double s = 0;
    for (std::size_t i = 0; i < KK; i++) s += gamma[i] = a[i] * beta[i];
    for (double& v : gamma) v /= s;
    visit(j, std::span<const double>(gamma));
    if (j == 0) break;
    // β̂_{j−1} from β̂_j.
    detail::pair_emissions(panel, j, g[j], c.epsilon, p1, e.data());
    std::fill(rows.begin(), rows.end(), 0.0);
    std::fill(cols.begin(), cols.end(), 0.0);
    double W = 0;
    for (std::size_t x = 0; x < K; x++)
      for (std::size_t y = 0; y < K; y++) {
        const std::size_t i = x * K + y;
        w[i] = e[i] * beta[i];
        rows[x] += w[i];
        cols[y] += w[i];
        W += w[i];
      }
    const double u = 1 - c.jump[j], q = c.move(j), uu = u * u, uq = u * q, qq = q * q;
    double d = 0;
    for (std::size_t x = 0; x < K; x++)
      for (std::size_t y = 0; y < K; y++) {
        const std::size_t i = x * K + y;
        d += beta[i] = uu * w[i] + uq * (rows[x] + cols[y]) + qq * W;
      }
    for (double& v : beta) v /= d;
  }
  return logp;
}

struct PairPosterior {
  std::size_t states = 0, sites = 0;
  std::vector<double> gamma;  // sites × K × K
  double log_likelihood = NEG_INF;
  double at(std::size_t j, std::size_t a, std::size_t b) const { return gamma[(j * states + a) * states + b]; }
};

// Every site's posterior over pairs (N K² doubles: for small panels and checks).
inline PairPosterior pair_posterior(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> g) {
  PairPosterior p;
  p.states = c.states;
  p.sites = c.sites;
  const std::size_t KK = c.states * c.states;
  p.gamma.assign(c.sites * KK, 0);
  p.log_likelihood = for_each_pair_posterior(panel, c, g, [&](std::size_t j, std::span<const double> gamma) {
    std::copy(gamma.begin(), gamma.end(), p.gamma.begin() + static_cast<std::ptrdiff_t>(j * KK));
  });
  return p;
}

// ---- phasing ---------------------------------------------------------------------------------------------

struct Phasing {
  std::vector<std::uint32_t> first, second;  // the copying paths A and B
  std::vector<std::int8_t> haplotype_first, haplotype_second;  // x and y, x_j + y_j = g_j where g_j is typed
  double log_probability = NEG_INF;  // log P(A, B, x, y) over the typed sites
};

// The jointly most probable copying paths and haplotypes of a diploid's genotypes.
inline Phasing phase(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> g) {
  detail::check_genotypes(c, g, "phase");
  const std::size_t K = c.states, N = c.sites, KK = K * K;
  std::vector<std::uint8_t> code(N * KK, 0);  // which of the four terms won at (j, c, d)
  std::vector<std::uint32_t> row_arg(N * K), col_arg(N * K), top_arg(N);
  std::vector<double> v(KK), e(KK), row_max(K), col_max(K);
  // log max_{x + y = g} P(x | h_a) P(y | h_b) takes four values at a site, one for each pair of
  // copied alleles (h_a, h_b): computed once, then read for every pair of haplotypes.
  const auto log_emissions = [&](std::size_t j) {
    double table[2][2];
    for (int ha = 0; ha < 2; ha++)
      for (int hb = 0; hb < 2; hb++)
        table[ha][hb] = safe_log(detail::genotype_emission(g[j], ha ? 1 - c.epsilon : c.epsilon, hb ? 1 - c.epsilon : c.epsilon, true));
    const std::int8_t* h = panel.alleles.data() + j * K;
    for (std::size_t x = 0; x < K; x++)
      for (std::size_t y = 0; y < K; y++) e[x * K + y] = table[h[x]][h[y]];
  };
  // The argmaxes of a row of δ over its rows, columns and all, for the step after site j.
  const auto summarise = [&](std::size_t j) {
    std::fill(row_max.begin(), row_max.end(), NEG_INF);
    std::fill(col_max.begin(), col_max.end(), NEG_INF);
    double top = NEG_INF;
    std::uint32_t* ra = &row_arg[j * K];
    std::uint32_t* ca = &col_arg[j * K];
    std::fill(ra, ra + K, 0u);
    std::fill(ca, ca + K, 0u);
    top_arg[j] = 0;
    for (std::size_t x = 0; x < K; x++)
      for (std::size_t y = 0; y < K; y++) {
        const double d = v[x * K + y];
        if (d > row_max[x]) {
          row_max[x] = d;
          ra[x] = static_cast<std::uint32_t>(y);
        }
        if (d > col_max[y]) {
          col_max[y] = d;
          ca[y] = static_cast<std::uint32_t>(x);
        }
        if (d > top) {
          top = d;
          top_arg[j] = static_cast<std::uint32_t>(x * K + y);
        }
      }
    return top;
  };
  log_emissions(0);
  const double l0 = -2 * std::log(static_cast<double>(K));
  for (std::size_t i = 0; i < KK; i++) v[i] = l0 + e[i];
  double top = summarise(0);
  for (std::size_t j = 1; j < N; j++) {
    log_emissions(j);
    const double ls = safe_log(c.stay(j)), lq = safe_log(c.move(j));
    const double w0 = 2 * ls, w1 = ls + lq, w3 = 2 * lq + top;
    std::uint8_t* cj = &code[j * KK];
    for (std::size_t x = 0; x < K; x++)
      for (std::size_t y = 0; y < K; y++) {
        const std::size_t i = x * K + y;
        double best = w0 + v[i];
        std::uint8_t which = 0;
        const double t1 = w1 + row_max[x], t2 = w1 + col_max[y];
        if (t1 > best) {
          best = t1;
          which = 1;
        }
        if (t2 > best) {
          best = t2;
          which = 2;
        }
        if (w3 > best) {
          best = w3;
          which = 3;
        }
        cj[i] = which;
        v[i] = best + e[i];
      }
    top = summarise(j);
  }
  Phasing out;
  out.log_probability = top;
  if (top == NEG_INF) return out;
  out.first.resize(N);
  out.second.resize(N);
  std::uint32_t at = top_arg[N - 1];
  for (std::size_t j = N; j-- > 0;) {
    const std::uint32_t x = at / static_cast<std::uint32_t>(K), y = at % static_cast<std::uint32_t>(K);
    out.first[j] = x;
    out.second[j] = y;
    if (j == 0) break;
    switch (code[j * KK + at]) {
      case 0: break;
      case 1: at = x * static_cast<std::uint32_t>(K) + row_arg[(j - 1) * K + x]; break;
      case 2: at = col_arg[(j - 1) * K + y] * static_cast<std::uint32_t>(K) + y; break;
      default: at = top_arg[j - 1]; break;
    }
  }
  out.haplotype_first.resize(N);
  out.haplotype_second.resize(N);
  for (std::size_t j = 0; j < N; j++) {
    const std::int8_t ha = panel.alleles[j * K + out.first[j]], hb = panel.alleles[j * K + out.second[j]];
    std::int8_t x = ha, y = hb;
    if (g[j] == 0) x = y = 0;
    else if (g[j] == 2) x = y = 1;
    else if (g[j] == 1) {
      // The larger of P(1 | h_a) P(0 | h_b) and P(0 | h_a) P(1 | h_b); allele 1 to x on a tie.
      const double p = ha ? 1 - c.epsilon : c.epsilon, r = hb ? 1 - c.epsilon : c.epsilon;
      x = p * (1 - r) >= (1 - p) * r;
      y = static_cast<std::int8_t>(1 - x);
    }
    out.haplotype_first[j] = x;
    out.haplotype_second[j] = y;
  }
  return out;
}

// log P(A, B, x, y) of given paths and haplotypes, over the typed sites (x_j + y_j = g_j there).
inline double phased_log_probability(const Haplotypes& panel, const Chain& c, std::span<const std::int8_t> g, const Phasing& p) {
  detail::check_genotypes(c, g, "phased_log_probability");
  const std::size_t K = c.states, N = c.sites;
  if (p.first.size() != N || p.second.size() != N || p.haplotype_first.size() != N || p.haplotype_second.size() != N)
    detail::fail("phased_log_probability", "paths and haplotypes of every site");
  const auto copy = [&](std::size_t j, std::uint32_t k, std::int8_t a) {
    return safe_log(panel.alleles[j * K + k] == a ? 1 - c.epsilon : c.epsilon);
  };
  double s = -2 * std::log(static_cast<double>(K));
  for (std::size_t j = 0; j < N; j++) {
    if (j) s += safe_log(c.transition(j, p.first[j - 1], p.first[j])) + safe_log(c.transition(j, p.second[j - 1], p.second[j]));
    if (g[j] < 0) continue;
    if (p.haplotype_first[j] + p.haplotype_second[j] != g[j]) return NEG_INF;
    s += copy(j, p.first[j], p.haplotype_first[j]) + copy(j, p.second[j], p.haplotype_second[j]);
  }
  return s;
}

// ---- switch errors ---------------------------------------------------------------------------------------

struct SwitchErrors {
  std::size_t heterozygous = 0;  // sites heterozygous in both the truth and the inference
  std::size_t switches = 0;
  double rate() const { return heterozygous > 1 ? static_cast<double>(switches) / static_cast<double>(heterozygous - 1) : 0; }
  SwitchErrors& operator+=(const SwitchErrors& o) {
    heterozygous += o.heterozygous;
    switches += o.switches;
    return *this;
  }
};

// The switches between the inferred haplotypes (x, y) and the true ones (t, u).
inline SwitchErrors switch_errors(std::span<const std::int8_t> t, std::span<const std::int8_t> u, std::span<const std::int8_t> x,
                                  std::span<const std::int8_t> y) {
  if (u.size() != t.size() || x.size() != t.size() || y.size() != t.size()) detail::fail("switch_errors", "four haplotypes of one length");
  SwitchErrors s;
  int last = -1;
  for (std::size_t j = 0; j < t.size(); j++) {
    if (t[j] < 0 || u[j] < 0 || x[j] < 0 || y[j] < 0 || t[j] == u[j] || x[j] == y[j]) continue;
    const int phase = x[j] == t[j];
    if (last >= 0 && phase != last) s.switches++;
    last = phase;
    s.heterozygous++;
  }
  return s;
}

}  // namespace evolution::haplotypes
