// Indices of a gene's adaptation to the codon usage of a reference set of genes: the codon
// adaptation index, the frequency of optimal codons, MILC, and Karlin and Mrázek's codon bias and
// predicted expression.
//
// - CAI (Sharp and Li 1987, Nucleic Acids Research 15, 1281–1295). Over the reference set, the
//   relative synonymous codon usage RSCU_ij = X_ij / ((1/n_i) Σ_j X_ij) (eq. 1) and the relative
//   adaptiveness w_ij = RSCU_ij / RSCU_imax = X_ij / X_imax (eq. 2), X_ij the count of codon j of
//   amino acid i; a codon absent from the reference set is given X = 0.5 (their p. 1285, so that
//   w is never 0). CAI = (Π_k w_k)^{1/L} over the gene's L codons (eq. 3), computed as exp of the
//   mean log (eq. 4), the codons of amino acids with one codon (Met and Trp in the standard code)
//   and stops left out, w being 1 for them by definition. Under a code where Met or Trp has two
//   codons (the vertebrate mitochondrial code, 2) their codons count.
// - Fop (Ikemura 1981, Journal of Molecular Biology 151, 389–409): the fraction of a gene's codons
//   that are optimal, over the codons of the amino acids that have an optimal codon. The optimal
//   codons are given, or found from a reference set as the codons of w = 1 (the most used codon of
//   each amino acid), or, as CodonW finds them (Peden 1999, §"Fop"), as the codons significantly
//   more frequent in highly than in lowly expressed genes: for codon c of amino acid a, the 2 × 2
//   table (c, a − c) × (high, low) by Pearson's χ² test on one degree of freedom, at level α.
// - MILC (Supek and Vlahoviček 2005, BMC Bioinformatics 6, 182): the contribution of each amino acid
//   M_a = 2 Σ_c O_c ln(O_c/E_c), the G statistic of the gene's codons against those expected from the
//   reference frequencies g_c within the amino acid, E_c = g_c Σ_{c∈a} O_c; MILC = Σ_a M_a / L − C
//   with C = Σ_a (r_a − 1)/L − 0.5 over the amino acids present (r_a its number of codons), the
//   expectation of the G statistic under the null correcting the bias of short genes. *Remark.*
//   The paper prints M_a without the factor 2 (its eq. 1) and C with + 0.5 (its eq. 4); the G
//   statistic, whose null expectation is r_a − 1, has the 2, and the authors' own implementation
//   (coRdon, Elek, Kuzman and Vlahoviček 2019, Bioconductor; its MILC method) uses 2 and − 0.5,
//   which we follow. The reference frequencies carry a pseudocount (0.5, as CAI's absent codons),
//   so that a codon absent from the reference set does not make MILC infinite. MELP = MILC(genome)
//   / MILC(reference) predicts expression (their eq. 5).
// - Codon bias of one set against another (Karlin and Mrázek 2000, Journal of Bacteriology 182,
//   5238–5250, eq. 1): B(F|G) = Σ_a p_a(F) Σ_{codons of a} |f(xyz) − g(xyz)|, f and g the codon
//   frequencies of F and G normalised within each amino acid and p_a(F) F's amino acid frequencies.
//   Predicted expression against a standard S, E_S(g) = B(g|C)/B(g|S) (eq. 2), C all genes; E(g) =
//   B(g|C) / (½ B(g|RP) + ¼ B(g|CH) + ¼ B(g|TF)) (eq. 3), RP the ribosomal proteins, CH the
//   chaperones and TF the translation and transcription factors. A gene is predicted highly
//   expressed (PHX, Definition I) when at least two of E_RP, E_CH, E_TF exceed 1.05 and E ≥ 1.00.
//   The paper applies it to genes of at least 100 codons. An amino acid absent from G is taken as
//   used uniformly there.
//
// Complexity: O(64) per gene; the contrast of optimal codons O(64).

#pragma once

#include <genomics/codon/usage.hpp>

#include <calx/special_functions/gamma.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace genomics::codon {

// ---- the codon adaptation index ----------------------------------------------------------------------
struct Adaptiveness {
  std::array<double, 64> w{};          // w_ij; 1 for single-codon amino acids, 0 for stops
  std::array<bool, 64> counted{};      // the codons CAI averages over
};

inline Adaptiveness relative_adaptiveness(const CodonCounts& reference, const GeneticCode& code, double absent = 0.5) {
  const Families fam = synonymous_families(code, false);
  Adaptiveness a;
  for (std::size_t f = 0; f < fam.size(); f++) {
    const auto& cs = fam.codons[f];
    double top = 0;
    for (int c : cs) top = std::max(top, reference[static_cast<std::size_t>(c)] > 0 ? reference[static_cast<std::size_t>(c)] : absent);
    for (int c : cs) {
      const double x = reference[static_cast<std::size_t>(c)] > 0 ? reference[static_cast<std::size_t>(c)] : absent;
      a.w[static_cast<std::size_t>(c)] = top > 0 ? x / top : 0;
      a.counted[static_cast<std::size_t>(c)] = cs.size() > 1;
    }
  }
  return a;
}

// CAI of a gene (its counts); 0 when it has no counted codon.
inline double cai(const CodonCounts& n, const Adaptiveness& a) {
  double s = 0, L = 0;
  for (std::size_t c = 0; c < 64; c++)
    if (a.counted[c] && n[c] > 0) {
      s += n[c] * std::log(a.w[c]);
      L += n[c];
    }
  return L > 0 ? std::exp(s / L) : 0;
}

// The w of each codon along a gene (0 for codons not counted), to paint it codon by codon.
inline std::vector<double> adaptiveness_profile(const Codes& cds, const Adaptiveness& a) {
  std::vector<double> out;
  out.reserve(cds.size() / 3);
  for (std::size_t i = 0; i + 3 <= cds.size(); i += 3) {
    const int c = sequence::codon_index(cds[i], cds[i + 1], cds[i + 2]);
    out.push_back(c >= 0 && a.counted[static_cast<std::size_t>(c)] ? a.w[static_cast<std::size_t>(c)] : 0);
  }
  return out;
}

// ---- optimal codons and Fop ---------------------------------------------------------------------------
using CodonSet = std::array<bool, 64>;

// The codons of w ≥ threshold among amino acids of more than one codon.
inline CodonSet optimal_by_adaptiveness(const Adaptiveness& a, double threshold = 1.0) {
  CodonSet s{};
  for (std::size_t c = 0; c < 64; c++) s[c] = a.counted[c] && a.w[c] >= threshold;
  return s;
}

struct OptimalContrast {
  CodonSet optimal{};
  std::array<double, 64> chi2{}, p{};          // Pearson's χ² of the 2 × 2 table and its p-value
  std::array<double, 64> high{}, low{};        // the codon's share of its amino acid in each set
};

inline OptimalContrast optimal_by_contrast(const CodonCounts& high, const CodonCounts& low, const GeneticCode& code, double alpha = 0.01) {
  const Families fam = synonymous_families(code, false);
  OptimalContrast out;
  for (std::size_t f = 0; f < fam.size(); f++) {
    const auto& cs = fam.codons[f];
    if (cs.size() < 2) continue;
    double H = 0, Lo = 0;
    for (int c : cs) { H += high[static_cast<std::size_t>(c)]; Lo += low[static_cast<std::size_t>(c)]; }
    if (!(H > 0 && Lo > 0)) continue;
    for (int c : cs) {
      const std::size_t k = static_cast<std::size_t>(c);
      const double a = high[k], b = H - a, cc = low[k], d = Lo - cc, n = H + Lo;
      out.high[k] = a / H;
      out.low[k] = cc / Lo;
      const double r1 = a + b, r2 = cc + d, c1 = a + cc, c2 = b + d;
      if (!(c1 > 0 && c2 > 0)) continue;
      const double chi2 = n * (a * d - b * cc) * (a * d - b * cc) / (r1 * r2 * c1 * c2);
      out.chi2[k] = chi2;
      out.p[k] = calx::special_functions::chi_square_q(chi2, 1);
      out.optimal[k] = out.high[k] > out.low[k] && out.p[k] < alpha;
    }
  }
  return out;
}

// Fop: optimal codons over the codons of amino acids that have an optimal codon.
inline double fop(const CodonCounts& n, const CodonSet& optimal, const GeneticCode& code) {
  const Families fam = synonymous_families(code, false);
  double opt = 0, all = 0;
  for (std::size_t f = 0; f < fam.size(); f++) {
    bool has = false;
    for (int c : fam.codons[f]) has = has || optimal[static_cast<std::size_t>(c)];
    if (!has) continue;
    for (int c : fam.codons[f]) {
      all += n[static_cast<std::size_t>(c)];
      if (optimal[static_cast<std::size_t>(c)]) opt += n[static_cast<std::size_t>(c)];
    }
  }
  return all > 0 ? opt / all : 0;
}

// ---- MILC --------------------------------------------------------------------------------------------
// The reference frequencies g_c within each amino acid, with a pseudocount per codon.
inline std::array<double, 64> within_family_frequencies(const CodonCounts& reference, const Families& fam, double pseudocount) {
  std::array<double, 64> g{};
  for (const auto& cs : fam.codons) {
    double total = 0;
    for (int c : cs) total += reference[static_cast<std::size_t>(c)] + pseudocount;
    for (int c : cs) g[static_cast<std::size_t>(c)] = total > 0 ? (reference[static_cast<std::size_t>(c)] + pseudocount) / total : 1.0 / static_cast<double>(cs.size());
  }
  return g;
}

inline double milc(const CodonCounts& n, const CodonCounts& reference, const GeneticCode& code, double pseudocount = 0.5) {
  const Families fam = synonymous_families(code, false);
  const auto g = within_family_frequencies(reference, fam, pseudocount);
  double M = 0, L = 0, correction = 0;
  for (const auto& cs : fam.codons) {
    double Oa = 0;
    for (int c : cs) Oa += n[static_cast<std::size_t>(c)];
    if (!(Oa > 0)) continue;
    L += Oa;
    correction += static_cast<double>(cs.size()) - 1;
    for (int c : cs) {
      const double O = n[static_cast<std::size_t>(c)];
      if (O > 0) M += 2 * O * std::log(O / (Oa * g[static_cast<std::size_t>(c)]));
    }
  }
  return L > 0 ? M / L - (correction / L - 0.5) : 0;
}

// ---- Karlin and Mrázek ------------------------------------------------------------------------------
inline double codon_bias(const CodonCounts& F, const CodonCounts& G, const GeneticCode& code) {
  const Families fam = synonymous_families(code, false);
  double total = 0;
  for (const auto& cs : fam.codons) {
    for (int c : cs) total += F[static_cast<std::size_t>(c)];
  }
  if (!(total > 0)) return 0;
  double B = 0;
  for (const auto& cs : fam.codons) {
    double nf = 0, ng = 0;
    for (int c : cs) { nf += F[static_cast<std::size_t>(c)]; ng += G[static_cast<std::size_t>(c)]; }
    if (!(nf > 0)) continue;
    double d = 0;
    for (int c : cs) {
      const double f = F[static_cast<std::size_t>(c)] / nf;
      const double g = ng > 0 ? G[static_cast<std::size_t>(c)] / ng : 1.0 / static_cast<double>(cs.size());
      d += std::fabs(f - g);
    }
    B += nf / total * d;
  }
  return B;
}

struct PredictedExpression {
  double b_all = 0, b_rp = 0, b_ch = 0, b_tf = 0;  // B(g|C), B(g|RP), B(g|CH), B(g|TF)
  double e_rp = 0, e_ch = 0, e_tf = 0, e = 0;      // eqs. 2 and 3
  bool phx = false;                                // Definition I
};

struct ExpressionStandards {
  CodonCounts all{}, ribosomal{}, chaperones{}, factors{};
};

inline PredictedExpression predicted_expression(const CodonCounts& gene, const ExpressionStandards& s, const GeneticCode& code) {
  PredictedExpression p;
  p.b_all = codon_bias(gene, s.all, code);
  p.b_rp = codon_bias(gene, s.ribosomal, code);
  p.b_ch = codon_bias(gene, s.chaperones, code);
  p.b_tf = codon_bias(gene, s.factors, code);
  const auto ratio = [](double a, double b) { return b > 0 ? a / b : 0; };
  p.e_rp = ratio(p.b_all, p.b_rp);
  p.e_ch = ratio(p.b_all, p.b_ch);
  p.e_tf = ratio(p.b_all, p.b_tf);
  p.e = ratio(p.b_all, 0.5 * p.b_rp + 0.25 * p.b_ch + 0.25 * p.b_tf);
  const int above = (p.e_rp > 1.05) + (p.e_ch > 1.05) + (p.e_tf > 1.05);
  p.phx = above >= 2 && p.e >= 1.0;
  return p;
}

}  // namespace genomics::codon
