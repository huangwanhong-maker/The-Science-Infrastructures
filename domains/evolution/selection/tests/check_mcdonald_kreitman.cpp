// Checks of the McDonald–Kreitman test: the Adh table of McDonald and Kreitman (1991), D_n = 7,
// D_s = 17, P_n = 2, P_s = 42, gives their G-test's P = 0.006 (with Williams' correction,
// G = 7.4327), Fisher's two-sided P = 0.007327 (as R's fisher.test gives), NI = 0.1156 and
// α = 0.8844; Fisher's tea-tasting table [[3, 1], [1, 3]] gives P = 0.4857 (Fisher 1935); and a
// table read from aligned sequences built by hand, with a column of two variable sites and
// columns left out.

#include <evolution/selection/mcdonald_kreitman.hpp>

#include <genomics/sequence/alphabet.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace evolution::selection;
using genomics::sequence::encode_dna;

int main() {
  calx::verification::Suite s;

  {
    MKTable adh;
    adh.Dn = 7;
    adh.Ds = 17;
    adh.Pn = 2;
    adh.Ps = 42;
    const MKTest t = mcdonald_kreitman(adh);
    // McDonald, Handbook of Biological Statistics (3rd ed., 2014), "Fisher's exact test": "we used
    // a G–test of independence in the original McDonald and Kreitman [1991] paper ... P=0.006".
    s.check("Adh: the G-test with Williams' correction gives McDonald and Kreitman's P = 0.006", t.g_williams_p - 0.006, 5e-4);
    s.check("its statistic is 7.4327 (q = 1.0638)", t.g_williams - 7.432694, 5e-6);
    s.check("uncorrected, G = 7.9071 and P = 0.004924", std::fabs(t.g - 7.907148) + std::fabs(t.g_p - 0.004924), 5e-6);
    s.check("Fisher's exact test, two-sided: P = 0.007327 (R's fisher.test)", t.fisher_p - 0.0073272, 5e-7);
    s.check("the neutrality index (2/42)/(7/17) = 0.115646", t.neutrality_index - (2.0 / 42) / (7.0 / 17), 1e-15);
    s.check("α = 1 − NI = 0.884354", t.alpha - 0.884354, 5e-7);
    s.check("the direction of selection 7/24 − 2/44", t.direction_of_selection - (7.0 / 24 - 2.0 / 44), 1e-15);
  }
  s.check("Fisher's tea-tasting table [[3, 1], [1, 3]]: P = 0.4857 (Fisher 1935, The Design of Experiments, §2)", fisher_exact(3, 1, 1, 3) - 34.0 / 70, 1e-14);
  s.check("a table without association has P = 1", fisher_exact(5, 5, 5, 5) - 1, 1e-14);

  {
    // Columns, by codon: 1 GCT | GCC fixed synonymous; 2 Lys AAA/AGA polymorphic in species 1, a
    // replacement; 3 Leu CTA/TTA polymorphic, synonymous; 4 GAT | GAA, Asp | Glu fixed, a
    // replacement; 5 ATG unvaried; 6 TTT/TTC polymorphic synonymous and T | G fixed at the first
    // site, a replacement (Phe | Val): two variable sites; 7 an N, left out; 8 CAT/CAC polymorphic
    // in species 2, synonymous; 9 a stop in one sequence, left out.
    const std::vector<Codes> one{encode_dna("GCTAAACTAGATATGTTTCATCATTGG"), encode_dna("GCTAAATTAGATATGTTCCNTCATTGG"),
                                 encode_dna("GCTAGACTAGATATGTTTCATCATTGA")};
    const std::vector<Codes> two{encode_dna("GCCAAACTAGAAATGGTTCATCATTGG"), encode_dna("GCCAAACTAGAAATGGTTCATCACTGG")};
    const MKTable t = mk_table(one, two);
    s.require("the table read from sequences: D_n = 2, D_s = 1, P_n = 1, P_s = 3",
              t.Dn == 2 && t.Ds == 1 && t.Pn == 1 && t.Ps == 3,
              std::to_string(t.Dn) + " " + std::to_string(t.Ds) + " " + std::to_string(t.Pn) + " " + std::to_string(t.Ps));
    s.require("with one column of two variable sites, seven columns compared and two left out", t.complex_codons == 1 && t.codons == 7 && t.skipped == 2);
  }
  return s.done();
}
