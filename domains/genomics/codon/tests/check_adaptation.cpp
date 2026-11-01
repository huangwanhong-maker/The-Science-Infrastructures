// Checks of the indices of adaptation to a reference set: CAI against Biopython 1.88 on Escherichia
// coli K-12 and the human mitochondrion (tests/fixtures/codon-*.tsv, make_codon_fixtures.py); Fop,
// the χ² contrast of optimal codons, MILC and Karlin and Mrázek's B and E against closed forms; the
// predicted highly expressed genes of E. coli.

#include "common.hpp"

#include <genomics/codon/adaptation.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <map>
#include <string>

using namespace genomics;
using namespace genomics::codon;

namespace {
std::size_t at(const char* c) { return static_cast<std::size_t>(sequence::codon_index(sequence::dna_code(c[0]), sequence::dna_code(c[1]), sequence::dna_code(c[2]))); }
}  // namespace

int main() {
  calx::verification::Suite s;
  const auto& standard = sequence::genetic_code(1);
  const auto& bacterial = sequence::genetic_code(11);

  // ---- CAI: closed forms ---------------------------------------------------------------------------------
  {
    CodonCounts ref{};
    ref[at("TTT")] = 30;
    ref[at("TTC")] = 10;  // Phe: w(TTT) = 1, w(TTC) = 1/3
    ref[at("GCT")] = 8;   // Ala: GCC, GCA, GCG absent: 0.5/8
    const auto a = relative_adaptiveness(ref, standard);
    s.check("w = X/X_max (Sharp and Li, eq. 2): w(TTC) = 1/3", a.w[at("TTC")] - 1.0 / 3, 1e-15);
    s.check("a codon absent from the reference set has X = 0.5: w(GCC) = 1/16", a.w[at("GCC")] - 1.0 / 16, 1e-15);
    s.require("Met and Trp, one codon each, and the stops are not counted", !a.counted[at("ATG")] && !a.counted[at("TGG")] && !a.counted[at("TAA")] && a.counted[at("TTC")]);
    CodonCounts gene{};
    gene[at("TTT")] = 2;
    gene[at("TTC")] = 1;
    gene[at("GCC")] = 1;
    gene[at("ATG")] = 5;
    s.check("CAI = (Π w)^{1/L} over counted codons (eq. 3): (1 · 1 · 1/3 · 1/16)^{1/4}", cai(gene, a) - std::pow(1.0 / 48, 0.25), 1e-15);
    const auto mt = relative_adaptiveness(ref, sequence::genetic_code(2));
    s.require("under the vertebrate mitochondrial code Met (ATA, ATG) and Trp (TGA, TGG) have two codons and are counted", mt.counted[at("ATG")] && mt.counted[at("TGA")]);
  }

  // ---- Fop and the contrast of optimal codons -------------------------------------------------------------
  {
    CodonSet opt{};
    opt[at("TTC")] = true;
    CodonCounts gene{};
    gene[at("TTT")] = 1;
    gene[at("TTC")] = 3;
    gene[at("GCT")] = 7;  // Ala has no optimal codon here: not counted
    s.check("Fop: optimal over the codons of amino acids that have one, 3/4", fop(gene, opt, standard) - 0.75, 1e-15);
    CodonCounts high{}, low{};
    high[at("TTC")] = 90;
    high[at("TTT")] = 10;
    low[at("TTC")] = 40;
    low[at("TTT")] = 60;
    const auto c = optimal_by_contrast(high, low, standard);
    // The 2 × 2 table ((90, 10), (40, 60)): χ² = 200 (90·60 − 10·40)² / (100 · 100 · 130 · 70) = 54.945…
    s.check("the contrast's Pearson χ² of the 2 × 2 table, by hand", c.chi2[at("TTC")] - 200.0 * 5000 * 5000 / (100.0 * 100 * 130 * 70), 1e-9);
    s.require("TTC, more used in the high set (0.9 against 0.4, p < 0.01), is optimal; TTT is not", c.optimal[at("TTC")] && !c.optimal[at("TTT")] && c.p[at("TTC")] < 1e-10);
  }

  // ---- MILC ------------------------------------------------------------------------------------------
  {
    CodonCounts gene{};
    gene[at("TTT")] = 6;
    gene[at("TTC")] = 2;
    gene[at("GCT")] = 4;
    gene[at("GCA")] = 4;
    // Against itself (pseudocount 0): every M_a = 0, MILC = −C = −(Σ (r_a − 1)/L − 0.5) = −((1 + 3)/16 − 0.5).
    s.check("MILC of a gene against its own usage is −C = 0.25", milc(gene, gene, standard, 0) - 0.25, 1e-15);
    CodonCounts ref{};
    ref[at("TTT")] = 1;
    ref[at("TTC")] = 1;
    ref[at("GCT")] = 1;
    ref[at("GCC")] = 1;
    ref[at("GCA")] = 1;
    ref[at("GCG")] = 1;
    // M_Phe = 2 (6 ln(6/4) + 2 ln(2/4)), M_Ala = 2 (4 ln(4/2) + 4 ln(4/2)).
    const double M = 2 * (6 * std::log(1.5) + 2 * std::log(0.5)) + 2 * (8 * std::log(2.0));
    s.check("MILC against uniform usage, by hand (the G statistic over L, less C)", milc(gene, ref, standard, 0) - (M / 16 - (4.0 / 16 - 0.5)), 1e-14);
  }

  // ---- Karlin and Mrázek's codon bias ------------------------------------------------------------------
  {
    CodonCounts F{}, G{};
    F[at("TTT")] = 3;
    F[at("TTC")] = 1;
    F[at("GCT")] = 4;
    G[at("TTT")] = 1;
    G[at("TTC")] = 1;
    G[at("GCT")] = 1;
    G[at("GCC")] = 1;
    // p_Phe = p_Ala = 1/2; Phe: |3/4 − 1/2| + |1/4 − 1/2| = 1/2; Ala: |1 − 1/2| + |0 − 1/2| = 1.
    s.check("B(F|G) by hand (eq. 1): ½ · ½ + ½ · 1", codon_bias(F, G, standard) - 0.75, 1e-15);
    s.check("B(F|F) = 0", codon_bias(F, F, standard), 1e-15);
  }

  // ---- E. coli against Biopython ----------------------------------------------------------------------
  const auto ecoli = check::genome("NC_000913.3");
  if (!ecoli) {
    s.skip("the NCBI genomes are not in the data cache (data/cache/ncbi-genomes): CAI is not compared with Biopython");
    return s.done();
  }
  const auto set = coding_genes(*ecoli);
  const auto rp = genes_by_product(set.genes, RIBOSOMAL_PROTEINS);
  s.require("the reference set: E. coli's 56 ribosomal protein genes by their /product", rp.size() == 56, std::to_string(rp.size()));
  const auto w = relative_adaptiveness(pooled_counts(set.genes, rp), bacterial);
  {
    double worst = 0;
    std::size_t n = 0;
    double worst_mt = 0;
    std::size_t n_mt = 0;
    const auto human = check::genome("NC_012920.1");
    Adaptiveness wmt;
    if (human) {
      const auto mt = coding_genes(*human, {2, 0});
      wmt = relative_adaptiveness(pooled_counts(mt.genes), sequence::genetic_code(2));
    }
    for (const auto& r : check::fixture("codon-weights.tsv", false)) {
      const auto& c = r[1];
      if (c.size() != 3) continue;
      const std::size_t k = at(c.c_str());
      if (r[0] == "cai-ecoli-ribosomal" && bacterial.amino_acids[k] != '*') {
        worst = std::max(worst, std::fabs(w.w[k] - std::stod(r[2])));
        n++;
      }
      if (r[0] == "cai-human-mitochondrion" && human && sequence::genetic_code(2).amino_acids[k] != '*') {
        worst_mt = std::max(worst_mt, std::fabs(wmt.w[k] - std::stod(r[2])));
        n_mt++;
      }
    }
    s.check("w of E. coli's " + std::to_string(n) + " sense codons from its ribosomal proteins against Biopython", worst, 1e-15);
    if (human) s.check("w of the human mitochondrion's " + std::to_string(n_mt) + " sense codons from its 13 genes (code 2) against Biopython", worst_mt, 1e-15);
  }
  {
    std::map<std::string, double> ref;
    for (const auto& r : check::fixture("codon-ecoli-genes.tsv")) ref[r[0]] = std::stod(r[2]);
    double worst = 0;
    std::size_t n = 0;
    for (const auto& g : set.genes) {
      const auto it = ref.find(check::key(g));
      if (it == ref.end()) continue;
      worst = std::max(worst, std::fabs(cai(codon_counts(g), w) - it->second));
      n++;
    }
    s.check("CAI of " + std::to_string(n) + " E. coli genes against Biopython's CodonAdaptationIndex", worst, 1e-12);
  }
  // Karlin and Mrázek's predicted highly expressed genes, for genes of at least 100 codons.
  {
    ExpressionStandards st;
    st.all = pooled_counts(set.genes);
    st.ribosomal = pooled_counts(set.genes, rp);
    st.chaperones = pooled_counts(set.genes, genes_by_name(set.genes, ecoli_chaperones()));
    st.factors = pooled_counts(set.genes, genes_by_name(set.genes, ecoli_translation_factors()));
    std::size_t long_genes = 0, phx = 0, rp_long = 0, rp_phx = 0;
    double top = 0;
    for (std::size_t i = 0; i < set.genes.size(); i++) {
      if (set.genes[i].codons() < 100) continue;
      long_genes++;
      const auto p = predicted_expression(codon_counts(set.genes[i]), st, bacterial);
      phx += p.phx;
      top = std::max(top, p.e);
      if (std::find(rp.begin(), rp.end(), i) != rp.end()) { rp_long++; rp_phx += p.phx; }
    }
    s.require("Karlin and Mrázek on E. coli: " + std::to_string(phx) + " of " + std::to_string(long_genes) + " genes of ≥ 100 codons predicted highly expressed (" +
                  std::to_string(100.0 * phx / long_genes).substr(0, 4) + "%; the paper, on the 1997 annotation: 306 of 3898, 8%), largest E " + std::to_string(top).substr(0, 4) +
                  " (2.66)",
              phx > 150 && phx < 500 && top > 1.8 && top < 3.5);
    s.require("most ribosomal proteins of ≥ 100 codons are predicted highly expressed", rp_phx * 10 >= rp_long * 8, std::to_string(rp_phx) + " of " + std::to_string(rp_long));
  }
  return s.done();
}
