// The usage of synonymous codons by a gene or a set of genes, under any of NCBI's genetic codes:
// the synonymous families of a code, the composition of third positions, the effective number of
// codons in Wright's form and in Sun, Yang and Xia's, and the codon bias scores of Roymondal et al.
// and of Sabi and Tuller, which need no reference set.
//
// - Synonymous families. The sense codons of a code partitioned by the amino acid they encode;
//   with split_by_box, an amino acid whose codons lie in several boxes (codons sharing their first
//   two bases) is split into one family per box: in the standard code leucine into TTR and CTN,
//   serine into TCN and AGY, arginine into CGN and AGR (Sun, Yang and Xia 2013, Molecular Biology
//   and Evolution 30, 191–196, who call these compound families and give the reason: a two-fold and
//   a four-fold family may be under quite different selection). Stops belong to no family.
// - Third positions. GC3, the fraction of G or C at the third position of sense codons; GC3s, the
//   same over the codons of amino acids with more than one codon (CodonW's "silent" third
//   positions: Met, Trp and stops excluded in the standard code; Peden 1999, CodonW documentation).
// - Wright's effective number of codons (Wright 1990, Gene 87, 23–29): for each amino acid with
//   n ≥ 2 observed codons, F̂ = (n Σ p² − 1)/(n − 1) (eq. 1); F̄_k the mean of F̂ over the amino acids
//   of degeneracy k; N_c = Σ_k K_k / F̄_k (eq. 2 for the standard code: 2 + 9/F̄₂ + 1/F̄₃ + 5/F̄₄ +
//   3/F̄₆). This is sequence::effective_number_of_codons, reused; uncapped, as dos Reis, Savva and
//   Wernisch (2004) needed it, N_c may exceed the number of sense codons in short genes.
// - Sun, Yang and Xia's N_c (2013, §"Materials and Methods"): families split by box; in each family
//   of m codons with counts n_i (n = Σ n_i), F_CF = Σ_i ((n_i + 1)/(n + m))², a pseudocount of one
//   for each codon; F̄_m the average of F_CF over the families of m codons weighted by their n;
//   N_c = Σ_m K_m / F̄_m. Its range is exactly [number of families, number of sense codons]: a
//   gene using its codons equally gives F_CF = 1/m for any n, and N_c the number of sense codons; a
//   gene using one codon a family gives F_CF → 1 as n grows, and N_c the number of families. A class
//   with no observed codon takes F̄_m = 1/m (no evidence of bias).
// - RCBS (Roymondal, Das and Sahoo 2009, DNA Research 16, 13–30, eqs. 2–3, as Sabi and Tuller 2014
//   write them): d_xyz = (f(xyz) − f₁(x) f₂(y) f₃(z)) / (f₁(x) f₂(y) f₃(z)), the frequency of a codon
//   against the product of its bases' frequencies at their positions in the same gene, and RCBS =
//   (Π_i (1 + d_i))^{1/L} − 1 over the gene's L codons. DCBS (Sabi and Tuller 2014, DNA Research 21,
//   511–526, eqs. 4–5): d_xyz = max(f/(f₁f₂f₃), f₁f₂f₃/f), so that a rare codon and a frequent one
//   both raise it, and DCBS = (1/L) Σ_i d_i. Both are computed over sense codons, stops excluded.
//
// Complexity: O(64) per gene once its counts are known.

#pragma once

#include <genomics/codon/genes.hpp>
#include <genomics/sequence/composition.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <vector>

namespace genomics::codon {

struct Families {
  int table = 1;
  bool split = false;
  std::array<int, 64> family{};            // the family of each codon; −1 for a stop
  std::vector<char> amino_acid;            // of each family
  std::vector<std::vector<int>> codons;    // of each family, in NCBI's order
  std::size_t size() const { return codons.size(); }
  std::size_t sense() const {
    std::size_t n = 0;
    for (int f : family) n += f >= 0;
    return n;
  }
};

inline Families synonymous_families(const GeneticCode& code, bool split_by_box = false) {
  Families out;
  out.table = code.id;
  out.split = split_by_box;
  out.family.fill(-1);
  std::map<std::pair<char, int>, int> index;  // (amino acid, box or −1) → family
  for (int c = 0; c < 64; c++) {
    const char aa = code.amino_acids[static_cast<std::size_t>(c)];
    if (aa == '*') continue;
    const int box = c >> 2;  // the first two bases
    bool several = false;
    if (split_by_box)
      for (int d = 0; d < 64; d++)
        if (code.amino_acids[static_cast<std::size_t>(d)] == aa && (d >> 2) != box) several = true;
    const std::pair<char, int> key{aa, several ? box : -1};
    auto it = index.find(key);
    if (it == index.end()) {
      it = index.emplace(key, static_cast<int>(out.codons.size())).first;
      out.codons.emplace_back();
      out.amino_acid.push_back(aa);
    }
    out.family[static_cast<std::size_t>(c)] = it->second;
    out.codons[static_cast<std::size_t>(it->second)].push_back(c);
  }
  return out;
}

// ---- third positions -------------------------------------------------------------------------------
struct ThirdPositions {
  double gc3 = 0, gc3s = 0;   // fractions
  double codons = 0, silent = 0;
};

inline ThirdPositions third_positions(const CodonCounts& n, const GeneticCode& code) {
  const Families fam = synonymous_families(code, false);
  ThirdPositions t;
  double gc = 0, gcs = 0;
  for (int c = 0; c < 64; c++) {
    const int f = fam.family[static_cast<std::size_t>(c)];
    if (f < 0) continue;
    const Code third = sequence::codon_of(c)[2];
    const double k = n[static_cast<std::size_t>(c)];
    const bool strong = third == sequence::C || third == sequence::G;
    t.codons += k;
    gc += strong ? k : 0;
    if (fam.codons[static_cast<std::size_t>(f)].size() > 1) {
      t.silent += k;
      gcs += strong ? k : 0;
    }
  }
  t.gc3 = t.codons > 0 ? gc / t.codons : 0;
  t.gc3s = t.silent > 0 ? gcs / t.silent : 0;
  return t;
}

// ---- effective numbers of codons ---------------------------------------------------------------------
// Wright's N_c, sequence::effective_number_of_codons on rounded counts.
inline double enc_wright(const CodonCounts& n, const GeneticCode& code, bool cap = true) {
  std::array<std::size_t, 64> k{};
  for (std::size_t i = 0; i < 64; i++) k[i] = static_cast<std::size_t>(std::llround(std::max(0.0, n[i])));
  return sequence::effective_number_of_codons(k, code, cap);
}

// Sun, Yang and Xia's N_c.
inline double enc_sun(const CodonCounts& n, const GeneticCode& code, bool split_by_box = true, double pseudocount = 1) {
  const Families fam = synonymous_families(code, split_by_box);
  std::map<std::size_t, std::pair<double, double>> by_m;  // m → (Σ n F, Σ n)
  std::map<std::size_t, int> K;                           // m → families
  for (std::size_t f = 0; f < fam.size(); f++) {
    const auto& cs = fam.codons[f];
    const std::size_t m = cs.size();
    K[m]++;
    if (m == 1) continue;
    double total = 0;
    for (int c : cs) total += n[static_cast<std::size_t>(c)];
    if (!(total > 0)) continue;
    double F = 0;
    for (int c : cs) {
      const double p = (n[static_cast<std::size_t>(c)] + pseudocount) / (total + pseudocount * static_cast<double>(m));
      F += p * p;
    }
    by_m[m].first += total * F;
    by_m[m].second += total;
  }
  double nc = 0;
  for (const auto& [m, families] : K) {
    if (m == 1) { nc += families; continue; }
    const auto it = by_m.find(m);
    const double Fbar = it != by_m.end() && it->second.second > 0 ? it->second.first / it->second.second : 1.0 / static_cast<double>(m);
    nc += families / Fbar;
  }
  return nc;
}

// ---- codon bias scores without a reference set --------------------------------------------------------
namespace detail {

// The frequency of each sense codon in a gene and the expected f₁(x) f₂(y) f₃(z) from its bases at
// their positions (sense codons only).
struct PositionalExpectation {
  std::array<double, 64> observed{}, expected{};
  double codons = 0;
};
inline PositionalExpectation positional(const CodonCounts& n, const GeneticCode& code) {
  PositionalExpectation e;
  std::array<std::array<double, 4>, 3> base{};
  for (int c = 0; c < 64; c++) {
    if (code.amino_acids[static_cast<std::size_t>(c)] == '*') continue;
    const double k = n[static_cast<std::size_t>(c)];
    const auto b = sequence::codon_of(c);
    for (int p = 0; p < 3; p++) base[static_cast<std::size_t>(p)][b[static_cast<std::size_t>(p)]] += k;
    e.codons += k;
  }
  if (!(e.codons > 0)) return e;
  for (int c = 0; c < 64; c++) {
    if (code.amino_acids[static_cast<std::size_t>(c)] == '*') continue;
    const auto b = sequence::codon_of(c);
    e.observed[static_cast<std::size_t>(c)] = n[static_cast<std::size_t>(c)] / e.codons;
    e.expected[static_cast<std::size_t>(c)] = base[0][b[0]] / e.codons * base[1][b[1]] / e.codons * base[2][b[2]] / e.codons;
  }
  return e;
}

}  // namespace detail

inline double rcbs(const CodonCounts& n, const GeneticCode& code) {
  const auto e = detail::positional(n, code);
  if (!(e.codons > 0)) return 0;
  double s = 0;
  for (std::size_t c = 0; c < 64; c++)
    if (n[c] > 0 && e.expected[c] > 0) s += n[c] * std::log(e.observed[c] / e.expected[c]);  // log(1 + d) = log(f/e)
  return std::exp(s / e.codons) - 1;
}

inline double dcbs(const CodonCounts& n, const GeneticCode& code) {
  const auto e = detail::positional(n, code);
  if (!(e.codons > 0)) return 0;
  double s = 0;
  for (std::size_t c = 0; c < 64; c++)
    if (n[c] > 0 && e.expected[c] > 0) s += n[c] * std::max(e.observed[c] / e.expected[c], e.expected[c] / e.observed[c]);
  return s / e.codons;
}

// RSCU of each codon (sequence::rscu on rounded counts; stops form a family of their own there).
inline std::array<double, 64> rscu(const CodonCounts& n, const GeneticCode& code) {
  std::array<std::size_t, 64> k{};
  for (std::size_t i = 0; i < 64; i++) k[i] = static_cast<std::size_t>(std::llround(std::max(0.0, n[i])));
  return sequence::rscu(k, code);
}

}  // namespace genomics::codon
