// Checks of Nei and Gojobori's method: the sites of TTA (2/3 and 7/3, NG86 p. 419), of a
// fourfold-degenerate codon, of ATG and TGG, and TGG under the vertebrate mitochondrial code;
// s + n = 3 for every sense codon; the differences of GTT and GTA (1, 0) and of TTT and GTA (0.5,
// 1.5), NG86's examples; TCA and TGG with the pathway through TGA left out (1, 1); symmetry; a
// pair of twelve-codon sequences worked by hand, S = 67/6, N = 149/6, S_d = 3.5, N_d = 2.5 and
// its distances in closed form; synonymous change alone giving ω = 0; codons with gaps and stops
// left out; saturation reported.

#include <evolution/selection/dnds.hpp>

#include <genomics/sequence/alphabet.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <string>

using namespace evolution::selection;
using genomics::sequence::encode_dna;

namespace {

int codon(const char* s) { return CodonTable::index({genomics::sequence::dna_code(s[0]), genomics::sequence::dna_code(s[1]), genomics::sequence::dna_code(s[2])}); }

}  // namespace

int main() {
  calx::verification::Suite s;
  const CodonTable standard(genomics::sequence::genetic_code(1));

  s.check("TTA (Leu) has 2/3 synonymous sites (Nei and Gojobori 1986, p. 419)", standard.synonymous_sites(codon("TTA")) - 2.0 / 3, 1e-15);
  s.check("and 7/3 non-synonymous ones", standard.nonsynonymous_sites(codon("TTA")) - 7.0 / 3, 1e-15);
  s.check("GCT (Ala, fourfold degenerate) has exactly 1 synonymous site", standard.synonymous_sites(codon("GCT")) - 1, 1e-15);
  s.check("ATG (Met) and TGG (Trp) have none in the standard code",
          standard.synonymous_sites(codon("ATG")) + standard.synonymous_sites(codon("TGG")), 0);
  const CodonTable mito(genomics::sequence::genetic_code(2));
  s.check("TGG has 1/3 under the vertebrate mitochondrial code, where TGA is Trp too", mito.synonymous_sites(codon("TGG")) - 1.0 / 3, 1e-15);
  {
    bool three = true, symmetric = true;
    for (int i = 0; i < 64; i++) {
      if (standard.stop(i)) continue;
      three = three && std::fabs(standard.synonymous_sites(i) + standard.nonsynonymous_sites(i) - 3) < 1e-15;
      for (int j = 0; j < 64; j++)
        if (!standard.stop(j))
          symmetric = symmetric && standard.synonymous_differences(i, j) == standard.synonymous_differences(j, i) &&
                      standard.nonsynonymous_differences(i, j) == standard.nonsynonymous_differences(j, i);
    }
    s.require("s + n = 3 for every sense codon", three);
    s.require("the differences of two codons are symmetric", symmetric);
  }
  s.require("GTT and GTA: one synonymous difference (NG86)",
            standard.synonymous_differences(codon("GTT"), codon("GTA")) == 1 && standard.nonsynonymous_differences(codon("GTT"), codon("GTA")) == 0);
  s.require("TTT and GTA: 0.5 synonymous and 1.5 non-synonymous, by the two pathways (NG86)",
            standard.synonymous_differences(codon("TTT"), codon("GTA")) == 0.5 && standard.nonsynonymous_differences(codon("TTT"), codon("GTA")) == 1.5);
  s.require("TCA and TGG: the pathway through the stop TGA left out, one of each",
            standard.synonymous_differences(codon("TCA"), codon("TGG")) == 1 && standard.nonsynonymous_differences(codon("TCA"), codon("TGG")) == 1);
  {
    // Worked by hand. x: GCT×8 TTA AAA TTT GTT; y: GCC GCT×7 CTA AGA GTA GTA.
    // Sites: x 8 + 2/3 + 1/3 + 1/3 + 1 = 31/3, y 1 + 7 + 4/3 + 2/3 + 1 + 1 = 12, S = 67/6, N = 36 − S = 149/6.
    // Differences: GCT/GCC s 1, TTA/CTA s 1, AAA/AGA n 1, TTT/GTA s 0.5 n 1.5, GTT/GTA s 1.
    const auto x = encode_dna("GCTGCTGCTGCTGCTGCTGCTGCTTTAAAATTTGTT");
    const auto y = encode_dna("GCCGCTGCTGCTGCTGCTGCTGCTCTAAGAGTAGTA");
    const NeiGojobori r = nei_gojobori(x, y, standard);
    s.check("S = 67/6 for the pair worked by hand", r.S - 67.0 / 6, 1e-13);
    s.check("N = 149/6", r.N - 149.0 / 6, 1e-13);
    s.check("S_d = 3.5", r.Sd - 3.5, 1e-15);
    s.check("N_d = 2.5", r.Nd - 2.5, 1e-15);
    s.check("d_S = −(3/4) ln(39/67)", r.dS + 0.75 * std::log(39.0 / 67), 1e-14);
    s.check("d_N = −(3/4) ln(129/149)", r.dN + 0.75 * std::log(129.0 / 149), 1e-14);
    s.check("ω = d_N / d_S", r.omega - std::log(129.0 / 149) / std::log(39.0 / 67), 1e-14);
    s.require("twelve codons compared", r.codons == 12 && r.skipped == 0);
  }
  {
    const auto x = encode_dna("GCTGGTCTGAAAATTCGTGCTGGTCTG");
    const auto y = encode_dna("GCCGGCCTAAAGATCCGCGCTGGTCTG");
    const NeiGojobori r = nei_gojobori(x, y, standard);
    s.require("synonymous change alone: N_d = 0 and ω = 0", r.Nd == 0 && r.Sd == 6 && r.omega == 0);
    Codes g = x, h = y;
    g[3] = genomics::sequence::GAP;
    h[6] = genomics::sequence::T;  // CTG → TTG in y: still Leu
    h[7] = genomics::sequence::A;  // TTG → TAG: a stop
    const NeiGojobori q = nei_gojobori(g, h, standard);
    s.require("codons with a gap, or a stop, in either sequence are left out", q.codons == 7 && q.skipped == 2);
  }
  {
    const auto x = encode_dna("TTTTTTTTTTTT"), y = encode_dna("TTCTTCTTCTTC");
    const NeiGojobori r = nei_gojobori(x, y, standard);
    s.require("p_S ≥ 3/4 is reported as saturation, d_S infinite, never a NaN", r.saturated_S && std::isinf(r.dS) && !std::isnan(r.dN));
  }
  return s.done();
}
