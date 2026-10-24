// Synonymous and non-synonymous divergence of two coding sequences: the method of Nei and
// Gojobori (1986, "Simple methods for estimating the numbers of synonymous and nonsynonymous
// nucleotide substitutions", Mol. Biol. Evol. 3, 418–426), under any NCBI genetic code.
//
// Sites. Of the three changes possible at each position of a codon, let f_i be the fraction that
// are synonymous; the codon has s = Σ f_i synonymous sites and n = 3 − s non-synonymous ones
// (NG86, p. 419: TTA has s = 2/3, n = 7/3). A change to a stop codon counts as non-synonymous,
// as in the paper. A sequence of r codons has S = Σ s, N = 3r − S; for a pair, the mean of the
// two sequences'.
//
// Differences. Two codons that differ at one position give one synonymous or one
// non-synonymous difference; at k = 2 or 3 positions, the k! pathways of single changes between
// them are weighed equally and s_d, n_d are the means of their counts (NG86: TTT and GTA, by
// TTT → GTT → GTA and TTT → TTA → GTA, give s_d = 0.5, n_d = 1.5). Pathways through a stop codon
// are left out, as in PAML (Yang 2007, Mol. Biol. Evol. 24, 1586–1591; tools.c, difcodonNG) and
// MATLAB's dnds; where every pathway passes a stop, PAML's rule is followed (two differences:
// s_d = 0, n_d = 2; three: s_d = 1, n_d = 2). The paper does not say how it treated them.
//
// Distances. p_S = S_d / S and p_N = N_d / N (eqs. 1 and 2), corrected for multiple hits by
// Jukes and Cantor, d = −(3/4) ln(1 − 4p/3) (eq. 3); ω = d_N / d_S. A codon with an ambiguous
// base or a gap, or that is a stop, in either sequence is left out (pairwise deletion of codons).
//
// Limits: the method assumes equal rates of transition and transversion and equal codon
// frequencies, and so is biased when they are not (Yang and Nielsen 2000, Mol. Biol. Evol. 17,
// 32–43); the Jukes–Cantor correction treats synonymous and non-synonymous sites as two
// independent JC69 processes. Saturation (p ≥ 3/4) is reported, the distance infinite.
// Complexity: a table of 64 codons and 64² pairs per genetic code, then O(codons).

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/genetic_code.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace evolution::selection {

using genomics::sequence::Code;
using genomics::sequence::Codes;
using genomics::sequence::GeneticCode;

// Sites and pathway differences of every codon of a genetic code, indexed in NCBI's order.
class CodonTable {
 public:
  explicit CodonTable(const GeneticCode& code) : code_(code) {
    for (int i = 0; i < 64; i++) aa_[static_cast<std::size_t>(i)] = code.amino_acids[static_cast<std::size_t>(i)];
    for (int i = 0; i < 64; i++) {
      const auto c = genomics::sequence::codon_of(i);
      double s = 0;
      for (int pos = 0; pos < 3; pos++)
        for (Code b = 0; b < 4; b++) {
          if (b == c[static_cast<std::size_t>(pos)]) continue;
          auto m = c;
          m[static_cast<std::size_t>(pos)] = b;
          s += amino(index(m)) == amino(i) ? 1.0 / 3 : 0.0;
        }
      syn_sites_[static_cast<std::size_t>(i)] = s;
    }
    for (int i = 0; i < 64; i++)
      for (int j = 0; j < 64; j++) pathways(i, j);
  }

  const GeneticCode& code() const { return code_; }
  char amino(int i) const { return aa_[static_cast<std::size_t>(i)]; }
  bool stop(int i) const { return amino(i) == '*'; }
  double synonymous_sites(int i) const { return syn_sites_[static_cast<std::size_t>(i)]; }
  double nonsynonymous_sites(int i) const { return 3 - syn_sites_[static_cast<std::size_t>(i)]; }
  double synonymous_differences(int i, int j) const { return sd_[static_cast<std::size_t>(i * 64 + j)]; }
  double nonsynonymous_differences(int i, int j) const { return nd_[static_cast<std::size_t>(i * 64 + j)]; }

  static int index(const std::array<Code, 3>& c) { return genomics::sequence::codon_index(c[0], c[1], c[2]); }

 private:
  GeneticCode code_;
  std::array<char, 64> aa_{};
  std::array<double, 64> syn_sites_{};
  std::array<double, 64 * 64> sd_{}, nd_{};

  // The mean counts over the pathways from i to j that avoid stop codons.
  void pathways(int i, int j) {
    const auto a = genomics::sequence::codon_of(i), b = genomics::sequence::codon_of(j);
    std::array<int, 3> diff{};
    int k = 0;
    for (int pos = 0; pos < 3; pos++)
      if (a[static_cast<std::size_t>(pos)] != b[static_cast<std::size_t>(pos)]) diff[static_cast<std::size_t>(k++)] = pos;
    double s = 0, n = 0;
    int valid = 0;
    if (k > 0) {
      std::array<int, 3> order{0, 1, 2};  // the order in which the k differing positions change
      do {
        auto cur = a;
        double ps = 0, pn = 0;
        bool through_stop = false;
        for (int step = 0; step < k; step++) {
          const int pos = diff[static_cast<std::size_t>(order[static_cast<std::size_t>(step)])];
          auto next = cur;
          next[static_cast<std::size_t>(pos)] = b[static_cast<std::size_t>(pos)];
          const int ni = index(next);
          if (step + 1 < k && stop(ni)) through_stop = true;
          (amino(index(cur)) == amino(ni) ? ps : pn) += 1;
          cur = next;
        }
        if (!through_stop) {
          s += ps;
          n += pn;
          valid++;
        }
      } while (std::next_permutation(order.begin(), order.begin() + k));
    }
    if (valid > 0) {
      s /= valid;
      n /= valid;
    } else if (k == 2) {
      s = 0;
      n = 2;
    } else if (k == 3) {
      s = 1;
      n = 2;
    }
    sd_[static_cast<std::size_t>(i * 64 + j)] = s;
    nd_[static_cast<std::size_t>(i * 64 + j)] = n;
  }
};

struct NeiGojobori {
  double S = 0, N = 0;     // synonymous and non-synonymous sites
  double Sd = 0, Nd = 0;   // synonymous and non-synonymous differences
  double pS = 0, pN = 0;   // the proportions of difference
  double dS = 0, dN = 0;   // with the Jukes–Cantor correction (infinite when saturated)
  double omega = std::numeric_limits<double>::quiet_NaN();  // dN / dS, NaN when not defined
  bool saturated_S = false, saturated_N = false;
  std::size_t codons = 0, skipped = 0;
};

namespace detail {

inline double jukes_cantor(double p, bool& saturated) {
  const double x = 1 - 4 * p / 3;
  saturated = !(x > 0);
  return saturated ? std::numeric_limits<double>::infinity() : -0.75 * std::log(x);
}

}  // namespace detail

// Two aligned coding sequences, read in frame from their first base (a trailing partial codon
// is ignored).
inline NeiGojobori nei_gojobori(const Codes& x, const Codes& y, const CodonTable& table) {
  if (x.size() != y.size()) throw std::invalid_argument("evolution::selection::nei_gojobori: aligned sequences of one length");
  NeiGojobori r;
  double sx = 0, sy = 0;
  for (std::size_t i = 0; i + 3 <= x.size(); i += 3) {
    const int a = genomics::sequence::codon_index(x[i], x[i + 1], x[i + 2]);
    const int b = genomics::sequence::codon_index(y[i], y[i + 1], y[i + 2]);
    if (a < 0 || b < 0 || table.stop(a) || table.stop(b)) {
      r.skipped++;
      continue;
    }
    r.codons++;
    sx += table.synonymous_sites(a);
    sy += table.synonymous_sites(b);
    r.Sd += table.synonymous_differences(a, b);
    r.Nd += table.nonsynonymous_differences(a, b);
  }
  r.S = (sx + sy) / 2;
  r.N = 3.0 * static_cast<double>(r.codons) - r.S;
  r.pS = r.S > 0 ? r.Sd / r.S : 0;
  r.pN = r.N > 0 ? r.Nd / r.N : 0;
  r.dS = detail::jukes_cantor(r.pS, r.saturated_S);
  r.dN = detail::jukes_cantor(r.pN, r.saturated_N);
  if (!r.saturated_S && !r.saturated_N && r.dS > 0) r.omega = r.dN / r.dS;
  else if (!r.saturated_N && r.dS == 0 && r.dN > 0) r.omega = std::numeric_limits<double>::infinity();
  return r;
}

inline NeiGojobori nei_gojobori(const Codes& x, const Codes& y, const GeneticCode& code = genomics::sequence::genetic_code(1)) {
  return nei_gojobori(x, y, CodonTable(code));
}

}  // namespace evolution::selection
