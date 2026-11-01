// Checks of the tRNA adaptation index and the tests of translational selection: the tRNA genes of
// Escherichia coli K-12 read from RefSeq, tAI's weights against codon-bias 0.5.0
// (tests/fixtures/codon-weights.tsv), the wobble rules by hand, Wright's expected N_c and dos Reis
// et al.'s refit recovered by simulating genes without selection, the S value with its permutation
// test, and stAI's hill climbing.

#include "common.hpp"

#include <genomics/codon/trna.hpp>

#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <chrono>
#include <cmath>
#include <map>
#include <string>

using namespace genomics;
using namespace genomics::codon;

namespace {
std::size_t at(const char* c) { return static_cast<std::size_t>(sequence::codon_index(sequence::dna_code(c[0]), sequence::dna_code(c[1]), sequence::dna_code(c[2]))); }
Isoacceptor iso(char aa, const char* anticodon, double copies) {
  Isoacceptor t;
  t.amino_acid = aa;
  t.anticodon = {sequence::dna_code(anticodon[0]), sequence::dna_code(anticodon[1]), sequence::dna_code(anticodon[2])};
  t.copies = copies;
  t.lysidine = aa == 'I' && std::string(anticodon) == "CAT";
  return t;
}
}  // namespace

int main() {
  calx::verification::Suite s;
  const auto& bacterial = sequence::genetic_code(11);

  // ---- the wobble rules, by hand ----------------------------------------------------------------------
  {
    // Arg ACG (I₃₄) ×4 reads CGU, CGC (s_I:C), CGA (s_I:A); Arg CCG ×1 reads CGG; Phe GAA ×2 reads UUC
    // and UUU (s_G:U); Ile CAU ×2 (lysidine) reads AUA only; Met CAU ×2 reads AUG only.
    const std::vector<Isoacceptor> pool = {iso('R', "ACG", 4), iso('R', "CCG", 1), iso('F', "GAA", 2), iso('I', "CAT", 2), iso('M', "CAT", 2)};
    const Penalties p;
    const auto a = trna_adaptiveness(pool, bacterial, p);
    s.check("W(CGT) = 4 (I:U)", a.W[at("CGT")] - 4, 1e-15);
    s.check("W(CGC) = 4 (1 − 0.28) (I:C)", a.W[at("CGC")] - 4 * 0.72, 1e-15);
    s.check("W(CGA) = 4 (1 − 0.9999) (I:A)", a.W[at("CGA")] - 4 * (1 - 0.9999), 1e-15);
    s.check("W(CGG) = 1 (C:G)", a.W[at("CGG")] - 1, 1e-15);
    s.check("W(TTT) = 2 (1 − 0.41) (G:U)", a.W[at("TTT")] - 2 * 0.59, 1e-15);
    s.check("W(ATA) = 2 (1 − 0.89): lysidine reads A", a.W[at("ATA")] - 2 * 0.11, 1e-15);
    s.require("lysidine does not read AUG, and tRNA-Met does not read AUA (cognate reads)", a.W[at("ATG")] == 2 && !a.counted[at("ATG")]);
    s.check("w = W/W_max (eq. 2): w(CGC) = 0.72", a.w[at("CGC")] - 0.72, 1e-15);
    // Codons no tRNA reads take the geometric mean of the others' w.
    double logs = 0;
    int k = 0;
    for (std::size_t c = 0; c < 64; c++)
      if (a.counted[c] && a.W[c] > 0) { logs += std::log(a.w[c]); k++; }
    s.check("unread codons take the geometric mean of the read ones' w (eq. 2)", a.w[at("GGG")] - std::exp(logs / k), 1e-15);
    CodonCounts g{};
    g[at("CGT")] = 1;
    g[at("CGC")] = 1;
    g[at("ATG")] = 3;  // methionine: not counted (eq. 3)
    s.check("tAI = (Π w)^{1/l} without methionine: √(1 · 0.72)", tai(g, a) - std::sqrt(0.72), 1e-15);
  }

  // ---- Wright's curve and dos Reis et al.'s ------------------------------------------------------------
  {
    s.check("Wright's f(0.5) = 2 + 0.5 + 29/0.5 = 60.5 (dos Reis et al., eq. 5)", wright_expected_nc(0.5) - 60.5, 1e-12);
    const NcCurve suggested;
    s.check("their f₁ with Table 3's suggested constants at x = 0.5: −6 + 0.5 + 34/(0.25 + 0.525²)", suggested(0.5) - (-5.5 + 34 / (0.25 + 0.525 * 0.525)), 1e-12);
  }

  const auto ecoli = check::genome("NC_000913.3");
  if (!ecoli) {
    s.skip("the NCBI genomes are not in the data cache (data/cache/ncbi-genomes): E. coli's tRNAs and genes are not read");
    return s.done();
  }
  const auto trnas = trna_genes(*ecoli);
  std::size_t initiators = 0, sec = 0;
  for (const auto& t : trnas.genes) { initiators += t.initiator; sec += t.amino_acid == "Sec"; }
  s.require("E. coli K-12's 86 tRNA genes read from RefSeq's notes: 4 initiators, 1 for selenocysteine", trnas.genes.size() == 86 && trnas.unread == 0 && initiators == 4 && sec == 1,
            std::to_string(trnas.genes.size()) + " read, " + std::to_string(trnas.unread) + " unread, " + std::to_string(initiators) + " initiators");
  const auto pool = trna_pool(trnas);
  double copies = 0;
  for (const auto& t : pool) copies += t.copies;
  s.require("its elongator pool: 81 genes in " + std::to_string(pool.size()) + " isoacceptors", copies == 81, std::to_string(copies));
  const auto a = trna_adaptiveness(pool, bacterial);
  {
    // codon-bias (with its eukaryotic rules: its prokaryotic row lets every C₃₄ read NNA as lysidine)
    // reads each anticodon for every codon it pairs with; away from AUA, which only lysidine reads,
    // and AUG, which tAI leaves out, the two agree up to the normalisation, so w ∝ W there.
    std::map<std::size_t, double> cb;
    for (const auto& r : check::fixture("codon-weights.tsv", false))
      if (r[0] == "tai-ecoli") cb[at(r[1].c_str())] = std::stod(r[2]);
    const std::size_t ref = at("CTG");
    double worst = 0;
    std::size_t n = 0;
    for (std::size_t c = 0; c < 64; c++) {
      if (!a.counted[c] || !(a.W[c] > 0) || c == at("ATA") || !cb.count(c)) continue;
      worst = std::max(worst, std::fabs(a.W[c] / a.W[ref] - cb[c] / cb[ref]));
      n++;
    }
    s.check("tAI's W of " + std::to_string(n) + " codons, relative to CUG's, against codon-bias's", worst, 1e-12);
    s.check("W(AUA) from E. coli's two tRNA-Ile(CAU) genes, lysidine: 2 (1 − 0.89)", a.W[at("ATA")] - 0.22, 1e-12);
    s.require("every counted codon is read by some tRNA in E. coli", a.unread == 0, std::to_string(a.unread));
  }

  const auto set = coding_genes(*ecoli);
  std::vector<CodonCounts> counts;
  for (const auto& g : set.genes) counts.push_back(codon_counts(g));
  calx::scheduler::ThreadPool threads(4);

  // ---- the expected N_c under no selection (dos Reis et al., eqs. 8–9) ------------------------------------
  {
    // The fit recovers a planted curve from noisy points.
    calx::random::Stream r(3, 0x504c4e54, 0, 0);
    const NcCurve planted{-5.0, 30.0, 1.04};
    std::vector<double> xs, ys;
    for (int i = 0; i < 4000; i++) {
      const double x = 0.05 + 0.9 * r.uniform();
      xs.push_back(x);
      ys.push_back(planted(x) + r.normal(0, 2));
    }
    const auto back = fit_nc_curve(xs, ys);
    double gap = 0;
    for (double x = 0.1; x <= 0.9001; x += 0.05) gap = std::max(gap, std::fabs(back.curve(x) - planted(x)));
    s.check("eq. 9's least squares recover a planted f₁ (−5, 30, 1.04) from 4000 points with noise of sd 2, largest gap over [0.1, 0.9]", gap, 0.3);
    // E. coli's genes simulated without selection.
    const auto fit = fit_nc_curve(counts, bacterial, 3, 7);
    const NcCurve table3{-6.459, 34.01, 1.023};  // E. coli K-12, their Table 3
    double worst = 0;
    for (double x = 0.2; x <= 0.8001; x += 0.05) worst = std::max(worst, std::fabs(fit.curve(x) - table3(x)));
    s.require("E. coli's genes (≥ 100 codons) simulated without selection, " + std::to_string(fit.nc.size()) + " of them: N_c at GC3s ½, where codons are used uniformly, is " +
                  std::to_string(fit.curve(0.5)).substr(0, 5) + ", near 61 (Table 3's E. coli curve gives " + std::to_string(table3(0.5)).substr(0, 5) + "; largest gap over [0.2, 0.8] " +
                  std::to_string(worst).substr(0, 4) + ")",
              std::fabs(fit.curve(0.5) - 61) < 2.5);
    double wright_gap = -1e9;
    for (double x = 0.2; x <= 0.8001; x += 0.05) wright_gap = std::max(wright_gap, std::fabs(wright_expected_nc(x) - fit.curve(x)));
    s.require("the refit and Wright's f differ by at most " + std::to_string(wright_gap).substr(0, 4) + " codons over [0.2, 0.8]", wright_gap < 8);
  }

  // ---- the S value -----------------------------------------------------------------------------------
  {
    const auto t0 = std::chrono::steady_clock::now();
    const auto sv = s_value(counts, a, bacterial, NcCurve{}, 999, 11, &threads);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    s.require("E. coli's S value, " + std::to_string(sv.s).substr(0, 5) + ", is above 0.5 and beyond all 999 permutations of w (p = 0.001; E. coli among the genomes of strongest translational selection, their Figure 3), " +
                  std::to_string(sec).substr(0, 4) + " s",
              sv.s > 0.5 && sv.as_extreme == 0 && std::fabs(sv.p_value - 0.001) < 1e-12);
    const auto again = s_value(counts, a, bacterial, NcCurve{}, 999, 11, nullptr);
    s.require("the permutations are the same on one thread and four", again.as_extreme == sv.as_extreme && again.s == sv.s);
  }

  // ---- stAI ------------------------------------------------------------------------------------------
  {
    std::vector<double> target;
    for (const auto& c : counts) target.push_back(dcbs(c, bacterial));
    const auto t0 = std::chrono::steady_clock::now();
    const auto r = optimise_penalties(counts, target, pool, bacterial, {}, &threads);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const auto w = r.penalties.wobble();
    bool inside = true;
    for (double v : w) inside = inside && v >= 0 && v <= 1;
    s.require("stAI on E. coli: ρ(tAI, DCBS) raised from " + std::to_string(r.start_correlation).substr(0, 5) + " at dos Reis et al.'s penalties to " +
                  std::to_string(r.correlation).substr(0, 5) + " (" + std::to_string(r.starts) + " starts, " + std::to_string(r.evaluations) + " evaluations, " + std::to_string(sec).substr(0, 4) +
                  " s), every penalty in [0, 1]",
              r.correlation >= r.start_correlation && inside && r.starts == 35);
    const auto again = optimise_penalties(counts, target, pool, bacterial, {}, nullptr);
    s.require("the hill climbing is the same on one thread and four", again.correlation == r.correlation && again.penalties.wobble() == w);
  }
  return s.done();
}
