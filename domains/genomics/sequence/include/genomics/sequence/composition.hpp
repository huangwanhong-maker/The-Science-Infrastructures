// The composition of a sequence, along it and over its genes: base counts and GC content;
// the GC skew and its cumulative sum, which locate a bacterial chromosome's origin and
// terminus of replication; CpG observed against expected; the entropy of k-mers; and the
// usage of codons, as relative synonymous codon usage and the effective number of codons.
//
// - GC skew (G − C)/(G + C) in windows (Lobry 1996, Molecular Biology and Evolution 13,
//   660–665). Its cumulative sum Σ (G − C) (Grigoriev 1998, Nucleic Acids Research 26,
//   2286–2290) falls along the lagging strand's template and rises along the leading one's, so
//   that in most bacteria its minimum marks the origin of replication and its maximum the
//   terminus (Compeau and Pevzner 2015, Bioinformatics Algorithms, ch. 1).
// - CpG observed/expected = n_CG · L / (n_C · n_G) (Gardiner-Garden and Frommer 1987,
//   Journal of Molecular Biology 196, 261–282), whose islands have o/e above 0.6 and GC above
//   0.5 over 200 bases.
// - The Shannon entropy of the k-mer distribution of a window, in bits, at most 2k.
// - RSCU, the relative synonymous codon usage: a codon's count over the mean count of the
//   codons of its amino acid (Sharp, Tuohy and Mosurski 1986, Nucleic Acids Research 14,
//   5125–5143).
// - The effective number of codons N_c (Wright 1990, Gene 87, 23–29): for each amino acid
//   with n ≥ 2 observed codons, F = (n Σ p_i² − 1)/(n − 1); the mean F of each class of
//   degeneracy; N_c = Σ_classes (amino acids in the class)/F̄, from 20 (one codon each) to 61
//   (no preference). Wright's handling of isoleucine when F̄₃ is missing (the mean of F̄₂ and
//   F̄₄) is followed; N_c is capped at 61.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/genetic_code.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <vector>

namespace genomics::sequence {

// The counts of A, C, G, T and N.
inline std::array<std::size_t, 5> base_counts(const Codes& x, std::size_t begin = 0, std::size_t end = SIZE_MAX) {
  std::array<std::size_t, 5> c{};
  end = std::min(end, x.size());
  for (std::size_t i = begin; i < end; i++) c[x[i] < 4 ? x[i] : 4]++;
  return c;
}

// G + C over A + C + G + T (N not counted); 0 for a sequence of Ns.
inline double gc_fraction(const Codes& x, std::size_t begin = 0, std::size_t end = SIZE_MAX) {
  const auto c = base_counts(x, begin, end);
  const double acgt = static_cast<double>(c[A] + c[C] + c[G] + c[T]);
  return acgt > 0 ? static_cast<double>(c[C] + c[G]) / acgt : 0;
}

// A value of a window [begin, end).
struct Window {
  std::size_t begin = 0, end = 0;
  double value = 0;
};

namespace detail {
template <class F>
std::vector<Window> windows(const Codes& x, std::size_t size, std::size_t step, F&& f) {
  std::vector<Window> out;
  if (size == 0 || step == 0) return out;
  for (std::size_t b = 0; b < x.size(); b += step) {
    const std::size_t e = std::min(x.size(), b + size);
    out.push_back({b, e, f(b, e)});
    if (e == x.size()) break;
  }
  return out;
}
}  // namespace detail

inline std::vector<Window> gc_windows(const Codes& x, std::size_t size, std::size_t step) {
  return detail::windows(x, size, step, [&](std::size_t b, std::size_t e) { return gc_fraction(x, b, e); });
}

inline std::vector<Window> gc_skew_windows(const Codes& x, std::size_t size, std::size_t step) {
  return detail::windows(x, size, step, [&](std::size_t b, std::size_t e) {
    const auto c = base_counts(x, b, e);
    const double s = static_cast<double>(c[G] + c[C]);
    return s > 0 ? (static_cast<double>(c[G]) - static_cast<double>(c[C])) / s : 0;
  });
}

// The cumulative skew Σ_{j < i} (G − C) at every `step`-th position i, from 0 to the length.
struct CumulativeSkew {
  std::vector<std::size_t> position;
  std::vector<double> value;
  std::size_t minimum = 0, maximum = 0;  // the positions of its extremes: the origin and terminus
};
inline CumulativeSkew cumulative_gc_skew(const Codes& x, std::size_t step = 1) {
  CumulativeSkew out;
  step = std::max<std::size_t>(1, step);
  double s = 0, lo = 0, hi = 0;
  out.position.push_back(0);
  out.value.push_back(0);
  for (std::size_t i = 0; i < x.size(); i++) {
    s += (x[i] == G) - (x[i] == C);
    if (s < lo) { lo = s; out.minimum = i + 1; }
    if (s > hi) { hi = s; out.maximum = i + 1; }
    if ((i + 1) % step == 0 || i + 1 == x.size()) {
      out.position.push_back(i + 1);
      out.value.push_back(s);
    }
  }
  return out;
}

// n_CG L / (n_C n_G) on [begin, end); 0 where C or G is absent.
inline double cpg_observed_expected(const Codes& x, std::size_t begin = 0, std::size_t end = SIZE_MAX) {
  end = std::min(end, x.size());
  std::size_t cg = 0, c = 0, g = 0, n = 0;
  for (std::size_t i = begin; i < end; i++) {
    if (x[i] > 3) continue;
    n++;
    c += x[i] == C;
    g += x[i] == G;
    if (i + 1 < end && x[i] == C && x[i + 1] == G) cg++;
  }
  return c && g ? static_cast<double>(cg) * static_cast<double>(n) / (static_cast<double>(c) * static_cast<double>(g)) : 0;
}

inline std::vector<Window> cpg_windows(const Codes& x, std::size_t size, std::size_t step) {
  return detail::windows(x, size, step, [&](std::size_t b, std::size_t e) { return cpg_observed_expected(x, b, e); });
}

// The Shannon entropy, in bits, of the k-mers (k ≤ 12) free of N in [begin, end).
inline double kmer_entropy(const Codes& x, int k, std::size_t begin = 0, std::size_t end = SIZE_MAX) {
  if (k < 1 || k > 12) throw std::invalid_argument("genomics::sequence::kmer_entropy: 1 <= k <= 12");
  end = std::min(end, x.size());
  std::vector<std::uint32_t> count(std::size_t{1} << (2 * k), 0);
  std::size_t total = 0;
  std::uint32_t code = 0;
  const std::uint32_t mask = (1u << (2 * k)) - 1;
  int valid = 0;
  for (std::size_t i = begin; i < end; i++) {
    if (x[i] > 3) { valid = 0; continue; }
    code = ((code << 2) | x[i]) & mask;
    if (++valid >= k) { count[code]++; total++; }
  }
  double h = 0;
  for (std::uint32_t c : count)
    if (c) {
      const double p = static_cast<double>(c) / static_cast<double>(total);
      h -= p * std::log2(p);
    }
  return h;
}

// ---- codon usage ------------------------------------------------------------------------------------
// The counts of the 64 codons, in NCBI's order, over coding sequences (each read from frame 0
// in whole codons; codons with N are skipped).
inline std::array<std::size_t, 64> codon_counts(const std::vector<Codes>& cds) {
  std::array<std::size_t, 64> n{};
  for (const auto& x : cds)
    for (std::size_t i = 0; i + 3 <= x.size(); i += 3) {
      const int k = codon_index(x[i], x[i + 1], x[i + 2]);
      if (k >= 0) n[static_cast<std::size_t>(k)]++;
    }
  return n;
}

// RSCU of each codon: its count over the mean count of its amino acid's codons (0 where the
// amino acid is never used; stops are treated as an amino acid of their own).
inline std::array<double, 64> rscu(const std::array<std::size_t, 64>& n, const GeneticCode& code) {
  std::array<double, 64> out{};
  std::map<char, std::pair<double, int>> by_aa;  // total count and number of codons
  for (int i = 0; i < 64; i++) {
    auto& t = by_aa[code.amino_acids[static_cast<std::size_t>(i)]];
    t.first += static_cast<double>(n[static_cast<std::size_t>(i)]);
    t.second++;
  }
  for (int i = 0; i < 64; i++) {
    const auto& t = by_aa[code.amino_acids[static_cast<std::size_t>(i)]];
    out[static_cast<std::size_t>(i)] = t.first > 0 ? static_cast<double>(n[static_cast<std::size_t>(i)]) * t.second / t.first : 0;
  }
  return out;
}

// Wright's effective number of codons.
inline double effective_number_of_codons(const std::array<std::size_t, 64>& n, const GeneticCode& code) {
  // Amino acids (stops excluded) by their degeneracy.
  std::map<char, std::vector<int>> codons;
  for (int i = 0; i < 64; i++)
    if (code.amino_acids[static_cast<std::size_t>(i)] != '*') codons[code.amino_acids[static_cast<std::size_t>(i)]].push_back(i);
  std::map<std::size_t, std::pair<double, int>> F;  // degeneracy -> sum of F, number of amino acids with F
  std::map<std::size_t, int> members;                // degeneracy -> number of amino acids
  for (const auto& [aa, cs] : codons) {
    const std::size_t d = cs.size();
    members[d]++;
    if (d < 2) continue;
    double total = 0, sum_sq = 0;
    for (int c : cs) total += static_cast<double>(n[static_cast<std::size_t>(c)]);
    if (total < 2) continue;
    for (int c : cs) {
      const double p = static_cast<double>(n[static_cast<std::size_t>(c)]) / total;
      sum_sq += p * p;
    }
    auto& f = F[d];
    f.first += (total * sum_sq - 1) / (total - 1);
    f.second++;
  }
  double nc = 0;
  for (const auto& [d, count] : members) {
    if (d == 1) { nc += count; continue; }
    double mean;
    if (F.count(d) && F[d].second > 0) mean = F[d].first / F[d].second;
    else if (d == 3 && F.count(2) && F.count(4) && F[2].second && F[4].second) mean = 0.5 * (F[2].first / F[2].second + F[4].first / F[4].second);
    else mean = 1.0 / static_cast<double>(d);  // no evidence: as if unbiased
    nc += count / std::max(mean, 1e-12);
  }
  return std::min(nc, 61.0);
}

}  // namespace genomics::sequence
