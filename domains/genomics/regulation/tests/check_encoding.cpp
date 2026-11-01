// Checks of the encodings of sequence and of the thermodynamic model: one-hot by hand; k-mer and
// gapped k-mer counts against brute force, with the number of canonical classes in closed form;
// the occupancy of promoters against Bintu et al.'s closed forms (2005, eqs. 4, 5, 8, 9, 10, and
// the dual repression of their Table 1).

#include <genomics/regulation/encoding.hpp>
#include <genomics/regulation/thermodynamic.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <map>
#include <string>

using namespace genomics;
using namespace genomics::regulation;

int main() {
  calx::verification::Suite s;

  // ---- one-hot ------------------------------------------------------------------------------------------
  {
    const auto x = one_hot<double>(sequence::encode_dna("ACGTN"));
    const std::vector<double> want = {1, 0, 0, 0, .25, 0, 1, 0, 0, .25, 0, 0, 1, 0, .25, 0, 0, 0, 1, .25};
    s.require("one-hot of ACGTN: rows A, C, G, T over the positions, N a quarter in each", x == want);
  }

  // ---- k-mers ---------------------------------------------------------------------------------------------
  {
    calx::random::Stream r(3, 0x4b4d4552, 0, 0);
    std::string text;
    for (int i = 0; i < 300; i++) text += "ACGT"[r.below(4)];
    text[100] = 'N';
    const Codes x = sequence::encode_dna(text);
    for (const bool canonical : {false, true}) {
      const KmerFeatures f(2, 4, canonical);
      std::vector<double> row(f.columns(), 0);
      f.count(x, row.data());
      // Brute force: every window without N, its word or the lesser of it and its reverse complement.
      std::map<std::string, double> brute;
      for (std::size_t k = 2; k <= 4; k++)
        for (std::size_t i = 0; i + k <= text.size(); i++) {
          std::string w = text.substr(i, k);
          if (w.find('N') != std::string::npos) continue;
          if (canonical) w = std::min(w, sequence::reverse_complement(std::string_view(w)));
          brute[w] += 1;
        }
      bool same = true;
      for (std::size_t c = 0; c < f.columns(); c++) same = same && row[c] == (brute.count(f.names()[c]) ? brute[f.names()[c]] : 0);
      // Canonical classes: (4^k + 4^{k/2})/2 for even k, 4^k/2 for odd.
      const std::size_t classes = canonical ? (16 + 4) / 2 + 64 / 2 + (256 + 16) / 2 : 16 + 64 + 256;
      s.require(std::string("k-mer counts (k = 2 … 4") + (canonical ? ", canonical" : "") + ") against brute force, with " + std::to_string(f.columns()) + " columns", same && f.columns() == classes);
    }
    for (const bool canonical : {false, true}) {
      const GappedKmerFeatures g(6, 4, canonical);
      std::vector<double> row(g.columns(), 0);
      g.count(x, row.data());
      // Brute force: windows of 6, masks with both ends informative and two of the four inner positions.
      std::map<std::string, double> brute;
      for (std::size_t i = 0; i + 6 <= text.size(); i++)
        for (int m = 0; m < 64; m++) {
          if (!(m & 1) || !(m & 32) || __builtin_popcount(static_cast<unsigned>(m)) != 4) continue;
          std::string w(6, 'N');
          bool ok = true;
          for (int j = 0; j < 6; j++)
            if (m >> j & 1) {
              w[static_cast<std::size_t>(j)] = text[i + static_cast<std::size_t>(j)];
              ok = ok && w[static_cast<std::size_t>(j)] != 'N';
            }
          if (!ok) continue;
          if (canonical) w = std::min(w, sequence::reverse_complement(std::string_view(w)));
          brute[w] += 1;
        }
      bool same = true;
      double total = 0;
      for (std::size_t c = 0; c < g.columns(); c++) {
        // A class's name is one of its two members; the brute force keys it by the lesser string.
        std::string key = g.names()[c];
        if (canonical) key = std::min(key, sequence::reverse_complement(std::string_view(key)));
        same = same && row[c] == (brute.count(key) ? brute[key] : 0);
        total += row[c];
      }
      s.require(std::string("gapped 4-mers in windows of 6 (C(4, 2) = 6 masks") + (canonical ? ", canonical" : "") + ") against brute force, " + std::to_string(g.columns()) + " columns, " +
                    std::to_string(static_cast<long>(total)) + " counted",
                same && g.masks() == 6 && (canonical || g.columns() == 6 * 256));
    }
  }

  // ---- the thermodynamic model ---------------------------------------------------------------------------
  {
    const double N = 5e6, P = 2000, e_pd = -6.0;
    Architecture basal;
    basal.nonspecific = N;
    basal.species = {{"RNAP", P, true}};
    basal.sites = {{0, e_pd, -35, 0, "promoter"}};
    s.check("RNA polymerase alone: p_bound = 1/(1 + (N_NS/P) e^{Δε}) (Bintu et al., eq. 4)", occupancy(basal).p_bound - p_bound_basal(P, N, e_pd), 1e-15);

    // An activator recruiting polymerase: F_reg = (1 + a f)/(1 + a) (eq. 8), p_bound by eq. 5.
    const double A = 50, e_ad = -13.0, e_ap = -4.0;
    Architecture act = basal;
    act.species.push_back({"activator", A, false});
    act.sites.push_back({1, e_ad, -70, -50, "CRP site"});
    act.interactions.push_back({0, 1, e_ap});
    const double a = A / N * std::exp(-e_ad), f = std::exp(-e_ap);
    s.check("activation: p_bound = eq. 5 with F_reg = (1 + a f)/(1 + a) (eq. 8)", occupancy(act).p_bound - p_bound_regulated(P, N, e_pd, regulation_activation(a, f)), 1e-15);
    // Weak promoter: the fold change is F_reg.
    Architecture weak = act;
    weak.sites[0].energy = 8.0;  // p_bound about 10⁻⁷ alone, 6 × 10⁻⁶ activated
    s.check("in the weak-promoter limit the fold change is F_reg, relative (p_bound F_reg ≈ 6 × 10⁻⁶)", fold_change(weak) / regulation_activation(a, f) - 1, 1e-4);

    // A repressor whose operator overlaps the promoter: F_reg = 1/(1 + r) (eq. 9).
    const double R = 10, e_rd = -15.0;
    Architecture rep = basal;
    rep.species.push_back({"repressor", R, false});
    rep.sites.push_back({1, e_rd, -10, 10, "operator"});
    const double r = R / N * std::exp(-e_rd);
    s.check("repression: eq. 5 with F_reg = 1/(1 + r) (eq. 9)", occupancy(rep).p_bound - p_bound_regulated(P, N, e_pd, regulation_repression(r)), 1e-15);

    // Two repressors, both overlapping the promoter but not each other, interacting with ω.
    Architecture dual = rep;
    dual.species.push_back({"repressor 2", 30, false});
    dual.sites.push_back({2, -14.0, 0, 20, "operator 2"});
    dual.sites[1].end = 0;  // the first operator −10 … 0, the second 0 … 20: both overlap the promoter, not each other
    dual.sites[0].end = 5;
    dual.interactions.push_back({1, 2, -2.0});
    const double r2 = 30 / N * std::exp(14.0), w = std::exp(2.0);
    s.check("two repressors interacting: F_reg = 1/(1 + r₁ + r₂ + ω r₁ r₂) (Bintu et al., Table 1)", occupancy(dual).p_bound - p_bound_regulated(P, N, e_pd, 1 / (1 + r + r2 + w * r * r2)), 1e-15);
    dual.interactions.clear();
    s.check("... and without interaction 1/((1 + r₁)(1 + r₂))", occupancy(dual).p_bound - p_bound_regulated(P, N, e_pd, 1 / ((1 + r) * (1 + r2))), 1e-15);

    // A helper recruiting the activator: its effective number A' = A (1 + h ω)/(1 + h) (eq. 10).
    Architecture help = act;
    help.species.push_back({"helper", 20, false});
    help.sites.push_back({2, -12.0, -100, -80, "helper site"});
    help.interactions.push_back({1, 2, -3.0});
    const double h = 20 / N * std::exp(12.0), wh = std::exp(3.0);
    const double a_eff = a * (1 + h * wh) / (1 + h);
    // With the helper, polymerase's weight is unchanged and the activator's a becomes a_eff: F_reg = (1 + a' f)/(1 + a').
    s.check("a helper: eq. 8 with the effective activator of eq. 10", occupancy(help).p_bound - p_bound_regulated(P, N, e_pd, regulation_activation(a_eff, f)), 1e-12);
    // Occupancies are probabilities that sum, over a site's states, to at most one.
    const auto o = occupancy(help);
    bool ok = true;
    for (double p : o.site) ok = ok && p >= 0 && p <= 1;
    s.require("every site's occupancy lies in [0, 1], " + std::to_string(o.states.size()) + " states of " + std::to_string(1u << help.sites.size()) + " allowed", ok && o.states.size() == 8);
    s.check("site energies from PWM scores: −Δε = λ ln 2 · bits + μ", site_energy(10, 0.5, 1) + (0.5 * std::log(2.0) * 10 + 1), 1e-15);
  }
  return s.done();
}
