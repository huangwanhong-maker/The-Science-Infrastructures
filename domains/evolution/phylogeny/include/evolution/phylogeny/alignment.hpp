// Aligned nucleotide sequences as phylogenetics reads them: rows of codes with their names,
// and the alignment compressed to its distinct columns (site patterns) with their weights.
//
// The rows are genomics::sequence codes (A = 0, C = 1, G = 2, T = 3), with N for any other
// letter and sequence::GAP for a gap; phylogenetics treats N and gaps alike as missing (the
// state unknown), as PHYLIP, PAUP* and RAxML do by default. A pattern is a column up to its
// sites: likelihoods, parsimony scores and pairwise distances are sums over sites of a
// function of the column, so they are sums over patterns weighted by their counts, and the
// nonparametric bootstrap is the same patterns with resampled weights (Felsenstein 2004,
// Inferring Phylogenies, on likelihood and on the bootstrap).
//
// compress() sorts the columns (O(L log L · n) for L sites and n rows) and numbers the
// patterns in the order of their first site, so that compression is deterministic.

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace evolution::phylogeny {

using genomics::sequence::Code;
using genomics::sequence::Codes;

inline constexpr Code MISSING = 4;  // N or a gap, in a pattern

struct Alignment {
  std::vector<std::string> names;
  std::vector<Codes> rows;

  std::size_t taxa() const { return rows.size(); }
  std::size_t sites() const { return rows.empty() ? 0 : rows[0].size(); }
};

// An alignment from letters: A, C, G, T (U), '-' and '.' as gaps, every other letter N.
inline Alignment alignment_from_letters(std::vector<std::string> names, const std::vector<std::string>& rows) {
  if (names.size() != rows.size()) throw std::invalid_argument("evolution::phylogeny::alignment_from_letters: a name for each row");
  Alignment a;
  a.names = std::move(names);
  for (const auto& r : rows) {
    Codes x = genomics::sequence::encode_dna(r);
    for (std::size_t i = 0; i < r.size(); i++)
      if (r[i] == '-' || r[i] == '.') x[i] = genomics::sequence::GAP;
    a.rows.push_back(std::move(x));
  }
  return a;
}

inline void validate(const Alignment& a) {
  if (a.names.size() != a.rows.size()) throw std::invalid_argument("evolution::phylogeny: an alignment has a name for each row");
  for (const auto& r : a.rows)
    if (r.size() != a.sites()) throw std::invalid_argument("evolution::phylogeny: the rows of an alignment have one length");
}

struct Patterns {
  std::vector<std::string> names;
  std::size_t taxa = 0, count = 0;
  std::vector<Code> states;           // count × taxa, pattern by pattern: 0–3, or MISSING
  std::vector<double> weight;         // count
  std::vector<std::uint32_t> of_site; // the pattern of each site of the alignment

  Code state(std::size_t pattern, std::size_t taxon) const { return states[pattern * taxa + taxon]; }
  double sites() const {
    double s = 0;
    for (double w : weight) s += w;
    return s;
  }
  // The row of a name, or SIZE_MAX.
  std::size_t index(std::string_view name) const {
    for (std::size_t i = 0; i < names.size(); i++)
      if (names[i] == name) return i;
    return SIZE_MAX;
  }
};

inline Patterns compress(const Alignment& a) {
  validate(a);
  Patterns p;
  p.names = a.names;
  p.taxa = a.taxa();
  const std::size_t L = a.sites(), n = p.taxa;
  std::vector<Code> column(L * n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t s = 0; s < L; s++) column[s * n + i] = a.rows[i][s] < 4 ? a.rows[i][s] : MISSING;
  std::vector<std::uint32_t> order(L);
  for (std::size_t s = 0; s < L; s++) order[s] = static_cast<std::uint32_t>(s);
  const auto less = [&](std::uint32_t x, std::uint32_t y) {
    const int c = n ? std::memcmp(&column[x * n], &column[y * n], n) : 0;
    return c != 0 ? c < 0 : x < y;
  };
  std::sort(order.begin(), order.end(), less);
  // Classes of equal columns, then numbered by first site.
  std::vector<std::uint32_t> cls(L), first;
  std::uint32_t k = 0;
  for (std::size_t r = 0; r < L; r++) {
    if (r > 0 && (n == 0 || std::memcmp(&column[order[r] * n], &column[order[r - 1] * n], n) != 0)) k++;
    if (r == 0 || k == first.size()) first.push_back(order[r]);
    cls[order[r]] = k;
  }
  std::vector<std::uint32_t> number(first.size());
  {
    std::vector<std::uint32_t> by(first.size());
    for (std::uint32_t c = 0; c < first.size(); c++) by[c] = c;
    std::sort(by.begin(), by.end(), [&](std::uint32_t x, std::uint32_t y) { return first[x] < first[y]; });
    for (std::uint32_t r = 0; r < by.size(); r++) number[by[r]] = r;
  }
  p.count = first.size();
  p.states.assign(p.count * n, 0);
  p.weight.assign(p.count, 0);
  p.of_site.resize(L);
  for (std::size_t c = 0; c < first.size(); c++)
    std::memcpy(&p.states[number[c] * n], &column[first[c] * n], n);
  for (std::size_t s = 0; s < L; s++) {
    p.of_site[s] = number[cls[s]];
    p.weight[number[cls[s]]] += 1;
  }
  return p;
}

// The frequencies of A, C, G, T over every known state, weighted by the patterns' weights; a
// base never seen gets a frequency of 1e-6 (renormalised), so that models stay defined.
inline std::array<double, 4> base_frequencies(const Patterns& p) {
  std::array<double, 4> f{0, 0, 0, 0};
  for (std::size_t q = 0; q < p.count; q++)
    for (std::size_t i = 0; i < p.taxa; i++) {
      const Code c = p.state(q, i);
      if (c < 4) f[c] += p.weight[q];
    }
  double s = 0;
  for (double x : f) s += x;
  if (s == 0) return {0.25, 0.25, 0.25, 0.25};
  double t = 0;
  for (double& x : f) {
    x = std::max(x / s, 1e-6);
    t += x;
  }
  for (double& x : f) x /= t;
  return f;
}

}  // namespace evolution::phylogeny
