// Checks of codon usage: the coding genes of real genomes, the synonymous families of the genetic
// codes, third positions, Wright's and Sun, Yang and Xia's effective numbers of codons, RCBS and
// DCBS; against closed forms, and against codon-bias 0.5.0 on Escherichia coli K-12
// (tests/fixtures/codon-ecoli-genes.tsv, codon-weights.tsv; make_codon_fixtures.py).

#include "common.hpp"

#include <genomics/codon/usage.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <map>
#include <string>

using namespace genomics;
using namespace genomics::codon;

int main() {
  calx::verification::Suite s;
  const auto& standard = sequence::genetic_code(1);
  const auto& bacterial = sequence::genetic_code(11);
  const auto& mito = sequence::genetic_code(2);

  // ---- families ------------------------------------------------------------------------------------
  {
    const auto f = synonymous_families(standard);
    const auto fs = synonymous_families(standard, true);
    const auto f2 = synonymous_families(mito);
    const auto f2s = synonymous_families(mito, true);
    std::map<std::size_t, int> sizes;
    for (const auto& c : f.codons) sizes[c.size()]++;
    s.require("standard code: 61 sense codons in 20 families, of sizes 1 (×2), 2 (×9), 3 (×1), 4 (×5), 6 (×3), Wright's classes",
              f.sense() == 61 && f.size() == 20 && sizes[1] == 2 && sizes[2] == 9 && sizes[3] == 1 && sizes[4] == 5 && sizes[6] == 3);
    std::map<std::size_t, int> split;
    for (const auto& c : fs.codons) split[c.size()]++;
    s.require("standard code split by boxes: Leu, Ser and Arg each into a 2-fold and a 4-fold family, 23 families (Sun, Yang and Xia 2013)",
              fs.size() == 23 && split[2] == 12 && split[4] == 8 && split[6] == 0, std::to_string(fs.size()));
    std::map<std::size_t, int> m2;
    for (const auto& c : f2.codons) m2[c.size()]++;
    s.require("vertebrate mitochondrial code: 60 sense codons (AGA, AGG stops), Met and Trp two-fold, Ile two-fold, 20 families",
              f2.sense() == 60 && f2.size() == 20 && m2[1] == 0 && f2.codons[static_cast<std::size_t>(f2.family[static_cast<std::size_t>(sequence::codon_index(sequence::A, sequence::T, sequence::A))])].size() == 2);
    s.require("vertebrate mitochondrial code split: Ser into TCN and AGY, Leu into TTR and CTN, 22 families", f2s.size() == 22, std::to_string(f2s.size()));
  }

  // ---- third positions -------------------------------------------------------------------------------
  {
    CodonCounts n{};
    const auto at = [](const char* c) { return static_cast<std::size_t>(sequence::codon_index(sequence::dna_code(c[0]), sequence::dna_code(c[1]), sequence::dna_code(c[2]))); };
    n[at("ATG")] = 4;  // Met: counts in GC3 (G), not in GC3s
    n[at("TGG")] = 2;  // Trp: likewise
    n[at("GCC")] = 3;  // Ala, C
    n[at("GCA")] = 1;  // Ala, A
    n[at("TAA")] = 5;  // a stop: neither
    const auto t = third_positions(n, standard);
    s.check("GC3 of (ATG ×4, TGG ×2, GCC ×3, GCA) is 9/10", t.gc3 - 0.9, 1e-15);
    s.check("GC3s leaves out Met, Trp and stops: 3/4", t.gc3s - 0.75, 1e-15);
  }

  // ---- effective numbers: closed forms --------------------------------------------------------------
  {
    const auto f = synonymous_families(standard);
    CodonCounts uniform{}, one{};
    for (const auto& cs : f.codons) {
      for (int c : cs) uniform[static_cast<std::size_t>(c)] = 12;  // 12 of every codon: p = 1/k exactly
      one[static_cast<std::size_t>(cs.front())] = 50;
    }
    // Wright: F̂ = (n Σp² − 1)/(n − 1) with Σp² = 1/k gives (n/k − 1)/(n − 1), a little under 1/k.
    s.check("Wright's N_c of a gene using one codon per amino acid is 20", enc_wright(one, standard) - 20, 1e-12);
    s.require("Wright's N_c of a gene using its codons equally is capped at 61, and above it uncapped (F̂ < 1/k)",
              enc_wright(uniform, standard) == 61 && enc_wright(uniform, standard, false) > 61, std::to_string(enc_wright(uniform, standard, false)));
    s.check("Sun's N_c of a gene using its codons equally is the number of sense codons, 61, for any count", enc_sun(uniform, standard) - 61, 1e-12);
    CodonCounts few{};
    for (const auto& cs : f.codons)
      for (int c : cs) few[static_cast<std::size_t>(c)] = 1;
    s.check("... and with one of each codon", enc_sun(few, standard) - 61, 1e-12);
    const double big = enc_sun([&] { CodonCounts x{}; for (const auto& cs : synonymous_families(standard, true).codons) x[static_cast<std::size_t>(cs.front())] = 1e7; return x; }(), standard);
    s.check("Sun's N_c of one codon per family tends to the number of families, 23, as counts grow (1e7 each)", big - 23, 1e-4);
    const auto f2 = synonymous_families(mito, true);
    CodonCounts u2{};
    for (const auto& cs : f2.codons)
      for (int c : cs) u2[static_cast<std::size_t>(c)] = 7;
    s.check("Sun's N_c under the vertebrate mitochondrial code: 60 sense codons", enc_sun(u2, mito) - 60, 1e-12);
    // Wright's eq. by hand for a gene of two amino acids.
    CodonCounts small{};
    const auto at = [](const char* c) { return static_cast<std::size_t>(sequence::codon_index(sequence::dna_code(c[0]), sequence::dna_code(c[1]), sequence::dna_code(c[2]))); };
    small[at("TTT")] = 3;
    small[at("TTC")] = 1;  // Phe: n = 4, Σp² = 10/16, F = (4·10/16 − 1)/3 = 1/2
    small[at("GCT")] = 2;
    small[at("GCC")] = 2;  // Ala: n = 4, Σp² = 1/2, F = (2 − 1)/3 = 1/3
    // N_c = 2 (Met, Trp) + 9/F₂ + 1/F₃ + 5/F₄ + 3/F₆ with F₂ = 1/2, F₄ = 1/3, F₃ = (F₂ + F₄)/2 = 5/12, F₆ = 1/6.
    s.check("Wright's N_c by hand: 2 + 18 + 12/5 + 15 + 18 = 55.4", enc_wright(small, standard) - 55.4, 1e-12);
  }

  // ---- RCBS and DCBS: closed forms ------------------------------------------------------------------
  {
    // A gene whose codons are exactly the products of its positional base frequencies: d = 0 for RCBS, 1 for DCBS.
    CodonCounts n{};
    const auto at = [](const char* c) { return static_cast<std::size_t>(sequence::codon_index(sequence::dna_code(c[0]), sequence::dna_code(c[1]), sequence::dna_code(c[2]))); };
    for (const char* c : {"GCT", "GCC", "GTT", "GTC", "ACT", "ACC", "ATT", "ATC"}) n[at(c)] = 5;  // {G,A} × {C,T} × {T,C}
    s.check("RCBS of a gene whose codons are independent at their positions is 0", rcbs(n, standard), 1e-12);
    s.check("DCBS of it is 1", dcbs(n, standard) - 1, 1e-12);
    CodonCounts m{};
    m[at("GCT")] = 3;
    m[at("ACC")] = 1;  // f(GCT) = 3/4, f₁(G) f₂(C) f₃(T) = 3/4 · 1 · 3/4; f(ACC) = 1/4, expected 1/4 · 1 · 1/4
    const double d1 = (0.75) / (0.75 * 0.75), d2 = 0.25 / (0.25 * 0.25);
    s.check("DCBS by hand: (3 · max(4/3, 3/4) + max(4, 1/4))/4", dcbs(m, standard) - (3 * d1 + d2) / 4, 1e-12);
    s.check("RCBS by hand: (4/3)^{3/4} · 4^{1/4} − 1", rcbs(m, standard) - (std::pow(d1, 0.75) * std::pow(d2, 0.25) - 1), 1e-12);
  }

  // ---- real genomes ------------------------------------------------------------------------------------
  const auto ecoli = check::genome("NC_000913.3");
  if (!ecoli) {
    s.skip("the NCBI genomes are not in the data cache (data/cache/ncbi-genomes): E. coli's genes are not compared with codon-bias");
    return s.done();
  }
  const auto set = coding_genes(*ecoli);
  s.require("E. coli K-12: 4294 coding genes counted, the same as Biopython's parse with the same exclusions (fixture)",
            set.genes.size() == 4294 && set.problems.empty(), std::to_string(set.genes.size()) + " genes; left out " + std::to_string(set.skipped.total()) +
                                                                  " (pseudo " + std::to_string(set.skipped.pseudo) + ", partial " + std::to_string(set.skipped.partial) +
                                                                  ", exceptions " + std::to_string(set.skipped.exception) + ", recoded " + std::to_string(set.skipped.transl_except) +
                                                                  ", internal stops " + std::to_string(set.skipped.internal_stop) + ")");
  {
    const auto rows = check::fixture("codon-ecoli-genes.tsv");
    std::map<std::string, std::vector<std::string>> by;
    for (const auto& r : rows) by[r[0]] = r;
    std::size_t found = 0, compared = 0;
    double worst = 0;
    // Where a degeneracy class has no amino acid seen twice, the two fill it differently (codon-bias
    // takes a missing F̄₃ from F̄₂ and F̄₄ after imputing them as 1/k; Wright's rule, kept here, only
    // from observed ones): the genes compared have every class of 2, 4 and 6 codons observed.
    const auto fam = synonymous_families(bacterial);
    const auto classes_observed = [&](const CodonCounts& n) {
      std::map<std::size_t, bool> seen = {{2, false}, {4, false}, {6, false}};
      for (const auto& cs : fam.codons) {
        double t = 0;
        for (int c : cs) t += n[static_cast<std::size_t>(c)];
        if (t >= 2 && seen.count(cs.size())) seen[cs.size()] = true;
      }
      return seen[2] && seen[4] && seen[6];
    };
    for (const auto& g : set.genes) {
      const auto it = by.find(check::key(g));
      if (it == by.end()) continue;
      found += std::stoul(it->second[1]) == g.codons();
      if (it->second[4] != "1" || !classes_observed(codon_counts(g))) continue;
      compared++;
      worst = std::max(worst, std::fabs(enc_wright(codon_counts(g), bacterial) - std::stod(it->second[3])));
    }
    s.require("every gene of the fixture found by its key, with its length in codons", found == rows.size() && found == set.genes.size(), std::to_string(found));
    s.check("Wright's N_c of " + std::to_string(compared) + " genes against codon-bias (genes whose classes are all observed and that its extra rule leaves alone)", worst, 1e-9);
    // Sun, Yang and Xia's N_c of the pooled genome, families unsplit, against codon-bias, which weights
    // each family by n + m (pseudocounts included) where we weight by n: equal to O(m/n).
    double pooled_sun = 0;
    for (const auto& r : check::fixture("codon-weights.tsv", false))
      if (r[0] == "nc-sun-ecoli-pooled") pooled_sun = std::stod(r[2]);
    const double mine = enc_sun(pooled_counts(set.genes), bacterial, false);
    s.check("Sun's N_c of E. coli's pooled genes (unsplit) against codon-bias, relative", (mine - pooled_sun) / pooled_sun, 1e-4);
    s.require("Sun's N_c of the pooled genome, split as the paper splits it, " + std::to_string(enc_sun(pooled_counts(set.genes), bacterial)) + ", lies in [23, 61]",
              enc_sun(pooled_counts(set.genes), bacterial) > 23 && enc_sun(pooled_counts(set.genes), bacterial) < 61);
  }
  {
    const auto human = check::genome("NC_012920.1");
    if (human) {
      const auto mt = coding_genes(*human, {2, 0});
      std::size_t codons = 0;
      for (const auto& g : mt.genes) codons += g.codons();
      s.require("the human mitochondrion: its 13 protein genes under code 2, the incomplete stops completed by polyadenylation dropped", mt.genes.size() == 13,
                std::to_string(mt.genes.size()) + " genes, " + std::to_string(codons) + " codons");
    }
  }
  return s.done();
}
