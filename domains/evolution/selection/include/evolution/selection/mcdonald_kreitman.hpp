// The McDonald–Kreitman test (McDonald and Kreitman 1991, "Adaptive protein evolution at the Adh
// locus in Drosophila", Nature 351, 652–654): under neutrality the ratio of replacement to
// synonymous change is the same for differences fixed between species as for polymorphisms
// within them, since both are the same mutations, sampled at different times; adaptive
// substitution raises the fixed ratio, slightly deleterious polymorphism the polymorphic one.
//
// The table. From aligned coding sequences of two species (several of one, one or more of the
// other), each nucleotide site of each codon is
//
//   polymorphic   when it varies within either species (one count per site);
//   fixed         when each species is monomorphic there and they differ;
//
// and synonymous when, in every codon observed in that column, changing the site between its
// observed bases leaves the amino acid unchanged, a replacement otherwise (a change to a stop
// codon is a replacement). For a column with one variable site this is McDonald and Kreitman's
// classification; a column with several (counted as `complex_codons`) is classified site by site
// in every context, so that a change is called a replacement if any context makes it one.
// Columns with an ambiguous base, a gap or a stop codon in any sequence are left out.
//
// The statistics. NI = (P_n/P_s)/(D_n/D_s), the neutrality index (Rand and Kann 1996, Mol. Biol.
// Evol. 13, 735–748); α = 1 − NI, the proportion of replacement substitutions fixed by
// selection (Smith and Eyre-Walker 2002, Nature 415, 1022–1024); DoS = D_n/(D_n + D_s) −
// P_n/(P_n + P_s), the direction of selection (Stoletzki and Eyre-Walker 2011, Mol. Biol. Evol. 28,
// 63–69). Fisher's exact test of the 2 × 2 table, two-sided: the sum of the hypergeometric
// probabilities of the tables with the observed margins that are no more probable than the
// observed one (with a relative tolerance of 10⁻⁷, as R's fisher.test); and the G-test of
// independence with Williams' correction, q = 1 + (n Σ 1/R_i − 1)(n Σ 1/C_j − 1)/(6n) for a
// 2 × 2 table (Sokal and Rohlf 1995, Biometry, 3rd ed., ch. 17; Williams 1976, Biometrika 63,
// 33–37), referred to χ² with one degree of freedom. McDonald and Kreitman used the G-test, and
// found P = 0.006 for Adh (McDonald, Handbook of Biological Statistics, 3rd ed., 2014, "Fisher's
// exact test").
//
// The hypergeometric probabilities are formed from calx's log Γ (calx has no hypergeometric
// distribution); exact for counts in the tens of thousands. O(codons × sequences).

#pragma once

#include <evolution/selection/dnds.hpp>

#include <calx/special_functions/gamma.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace evolution::selection {

struct MKTable {
  std::size_t Dn = 0, Ds = 0;  // fixed: replacement, synonymous
  std::size_t Pn = 0, Ps = 0;  // polymorphic: replacement, synonymous
  std::size_t codons = 0;          // codon columns compared
  std::size_t complex_codons = 0;  // of them, with more than one variable site
  std::size_t skipped = 0;         // columns left out
};

struct MKTest {
  MKTable table;
  double neutrality_index = 0, alpha = 0, direction_of_selection = 0;
  double fisher_p = 1;               // two-sided
  double g = 0, g_p = 1;             // the G-test without correction
  double g_williams = 0, g_williams_p = 1;
};

// The table of fixed and polymorphic, replacement and synonymous sites between two species.
inline MKTable mk_table(const std::vector<Codes>& species1, const std::vector<Codes>& species2, const CodonTable& table) {
  if (species1.empty() || species2.empty()) throw std::invalid_argument("evolution::selection::mk_table: sequences of both species");
  const std::size_t L = species1[0].size();
  for (const auto* g : {&species1, &species2})
    for (const auto& r : *g)
      if (r.size() != L) throw std::invalid_argument("evolution::selection::mk_table: aligned sequences of one length");
  MKTable t;
  std::vector<int> codons1, codons2;
  for (std::size_t i = 0; i + 3 <= L; i += 3) {
    codons1.clear();
    codons2.clear();
    bool ok = true;
    const auto read = [&](const std::vector<Codes>& g, std::vector<int>& out) {
      for (const auto& r : g) {
        const int c = genomics::sequence::codon_index(r[i], r[i + 1], r[i + 2]);
        if (c < 0 || table.stop(c)) ok = false;
        else out.push_back(c);
      }
    };
    read(species1, codons1);
    read(species2, codons2);
    if (!ok) {
      t.skipped++;
      continue;
    }
    t.codons++;
    int variable = 0;
    for (int pos = 0; pos < 3; pos++) {
      // The bases at this position in each species, as 4-bit sets.
      unsigned a1 = 0, a2 = 0;
      for (int c : codons1) a1 |= 1u << genomics::sequence::codon_of(c)[static_cast<std::size_t>(pos)];
      for (int c : codons2) a2 |= 1u << genomics::sequence::codon_of(c)[static_cast<std::size_t>(pos)];
      const bool poly = std::popcount(a1) > 1 || std::popcount(a2) > 1;
      const bool fixed = !poly && a1 != a2;
      if (!poly && !fixed) continue;
      variable++;
      // Synonymous when no observed codon changes its amino acid with this site among its bases.
      const unsigned bases = a1 | a2;
      bool synonymous = true;
      for (const auto* cs : {&codons1, &codons2})
        for (int c : *cs)
          for (Code b = 0; b < 4 && synonymous; b++) {
            if (!((bases >> b) & 1)) continue;
            auto m = genomics::sequence::codon_of(c);
            m[static_cast<std::size_t>(pos)] = b;
            if (table.amino(CodonTable::index(m)) != table.amino(c)) synonymous = false;
          }
      if (poly) (synonymous ? t.Ps : t.Pn)++;
      else (synonymous ? t.Ds : t.Dn)++;
    }
    t.complex_codons += variable > 1;
  }
  return t;
}

inline MKTable mk_table(const std::vector<Codes>& species1, const std::vector<Codes>& species2,
                        const GeneticCode& code = genomics::sequence::genetic_code(1)) {
  return mk_table(species1, species2, CodonTable(code));
}

// The two-sided p-value of Fisher's exact test of the table [[a, b], [c, d]].
inline double fisher_exact(std::size_t a, std::size_t b, std::size_t c, std::size_t d) {
  namespace sf = calx::special_functions;
  const double r1 = static_cast<double>(a + b), r2 = static_cast<double>(c + d), c1 = static_cast<double>(a + c);
  const double n = r1 + r2;
  if (n == 0) return 1;
  const auto logp = [&](double x) { return sf::log_choose(r1, x) + sf::log_choose(r2, c1 - x) - sf::log_choose(n, c1); };
  const double lo = std::max(0.0, c1 - r2), hi = std::min(r1, c1);
  const double observed = logp(static_cast<double>(a));
  double p = 0;
  for (double x = lo; x <= hi; x += 1) {
    const double l = logp(x);
    if (l <= observed + 1e-7) p += std::exp(l);
  }
  return std::min(1.0, p);
}

inline MKTest mcdonald_kreitman(const MKTable& t) {
  MKTest r;
  r.table = t;
  const double Dn = static_cast<double>(t.Dn), Ds = static_cast<double>(t.Ds), Pn = static_cast<double>(t.Pn), Ps = static_cast<double>(t.Ps);
  r.neutrality_index = Ps > 0 && Dn > 0 ? (Pn / Ps) / (Dn / Ds) : std::numeric_limits<double>::quiet_NaN();
  r.alpha = 1 - r.neutrality_index;
  r.direction_of_selection = (Dn + Ds > 0 ? Dn / (Dn + Ds) : 0) - (Pn + Ps > 0 ? Pn / (Pn + Ps) : 0);
  r.fisher_p = fisher_exact(t.Dn, t.Ds, t.Pn, t.Ps);
  // G = 2 Σ O ln(O / E), E = row × column / n.
  const double n = Dn + Ds + Pn + Ps;
  const std::array<double, 4> O{Dn, Ds, Pn, Ps};
  const std::array<double, 2> R{Dn + Ds, Pn + Ps}, C{Dn + Pn, Ds + Ps};
  if (n > 0 && R[0] > 0 && R[1] > 0 && C[0] > 0 && C[1] > 0) {
    double g = 0;
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++) {
        const double o = O[static_cast<std::size_t>(i * 2 + j)], e = R[static_cast<std::size_t>(i)] * C[static_cast<std::size_t>(j)] / n;
        if (o > 0) g += o * std::log(o / e);
      }
    r.g = 2 * g;
    const double q = 1 + (n * (1 / R[0] + 1 / R[1]) - 1) * (n * (1 / C[0] + 1 / C[1]) - 1) / (6 * n);
    r.g_williams = r.g / q;
    r.g_p = calx::special_functions::chi_square_q(r.g, 1);
    r.g_williams_p = calx::special_functions::chi_square_q(r.g_williams, 1);
  }
  return r;
}

}  // namespace evolution::selection
