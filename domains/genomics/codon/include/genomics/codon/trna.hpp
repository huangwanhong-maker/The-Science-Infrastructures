// The adaptation of genes to the tRNA pool of their genome: the tRNA adaptation index, the S value
// of translational selection, the expected effective number of codons under no selection, and the
// optimisation of the wobble penalties from the genome alone.
//
// - The tRNA pool. The tRNA genes of the genome (genes.hpp), grouped by amino acid and anticodon:
//   isoacceptors with their gene copy numbers tGCN. Initiator tRNAs (fMet) read only the start codon
//   and selenocysteine's tRNA only a recoded UGA; neither takes part in elongation, and both are
//   left out.
// - Absolute and relative adaptiveness (dos Reis, Savva and Wernisch 2004, Nucleic Acids Research
//   32, 5036–5044, eqs. 1–2): W_i = Σ_j (1 − s_ij) tGCN_ij over the isoacceptors j that read codon i,
//   s_ij the penalty of the codon–anticodon pairing; w_i = W_i/W_max, and for codons with W_i = 0
//   the geometric mean of the other w (eq. 2). Who reads what is their Table 1, Crick's wobble
//   rules (Crick 1966, J. Mol. Biol. 19, 548–555) with inosine: an anticodon N₃₄N₃₅N₃₆ reads the
//   codons comp(N₃₆) comp(N₃₅) X with
//     I₃₄ (a genomic A, deaminated):  X = U (s_I:U), C (s_I:C), A (s_I:A)
//     G₃₄:                            X = C (s_G:C), U (s_G:U)
//     U₃₄:                            X = A (s_U:A), G (s_U:G)
//     C₃₄:                            X = G (s_C:G)
//     L₃₄ (lysidine, the C₃₄ of a bacterial tRNA-Ile):  X = A (s_L:A) only
//   and, here, only codons of the tRNA's own amino acid under the genome's code (a cognate read; it
//   keeps tRNA-Met from AUA and tRNA-Ile2 from AUG). The penalties default to their Table 2's
//   optimised values: s_G:U = 0.41, s_I:C = 0.28, s_I:A = 0.9999, s_U:G = 0.68, s_L:A = 0.89, the
//   Watson–Crick ones 0. *Remark.* Their R package (tAI 0.21, get.ws) sets W(AUA) = 1 − s_L:A in a
//   bacterium, as for one lysidine tRNA whatever the genome holds; here the genome's tRNA-Ile(CAU)
//   genes are counted (two in E. coli K-12, ileX and ileY).
// - tAI_g = (Π_k w_{i_kg})^{1/l_g} over the gene's codons, stops and methionine's codons left out
//   (eq. 3; the package's get.tai).
// - Wright's expected N_c under no selection, f(x) = 2 + x + 29/(x² + (1 − x)²) with x the gene's
//   GC3s (dos Reis et al.'s eq. 5, after Wright 1990), and their refit f₁(x) = a + x + b/(x² + (c −
//   x)²) (eq. 8), whose constants they fitted by least squares (eq. 9) to genes simulated with their
//   amino acids kept and their codons drawn by a silent GC content chosen uniformly; their Table 3
//   gives a = −6.459, b = 34.01, c = 1.023 for E. coli and suggests −6.0, 34.0, 1.025. The
//   simulation is fit_nc_curve below: each synonymous codon drawn with probability proportional to
//   x or 1 − x as its third base is strong or weak, our reading of "according to a silent GC
//   content"; only genes of at least min_codons codons whose N_c is defined (every degeneracy class
//   with a positive mean F̂) enter the fit. *Caveat.* Our simulation does not reproduce Table 3:
//   under it a gene of GC3s ½ uses its codons uniformly and its N_c averages about 61, where
//   Table 3's curve gives 59; the paper does not say enough of its simulation (nor of CodonW's
//   treatment of an F̄ ≤ 0) to find why. The S value therefore uses Table 3's suggested constants.
// - The S value (eqs. 6–7): ψ_g = f₁(GC3s_g) − N_c,g, the reduction of a gene's N_c below its
//   expectation, and S = cor(tAI, ψ), Pearson's correlation over the genes; its significance by
//   permuting the w among the codons (their §"The relationship of tAI to Nc"), the p-value (b + 1)/(B
//   + 1) (Phipson and Smyth 2010, Stat. Appl. Genet. Mol. Biol. 9, art. 39, eq. 1). N_c is Wright's,
//   capped at the number of sense codons as Wright and CodonW cap it. *Remark.* dos Reis et al.
//   lifted the cap; but F̂ = (n Σp² − 1)/(n − 1) is ≤ 0 whenever an amino acid seen n times uses n
//   different codons, and for isoleucine, alone in its class, that makes 1/F̄₃ and N_c undefined
//   (in E. coli K-12 the uncapped values reach 10¹², and S falls from 0.65 to 0.15): the cap is the
//   definition that keeps N_c finite.
// - stAI (Sabi and Tuller 2014, DNA Research 21, 511–526, §2.3): the wobble penalties chosen, for a
//   genome without measured expression, to maximise Spearman's correlation between tAI and DCBS
//   (usage.hpp) over the genes, by hill climbing from starting points: the original values, all
//   zeros, all ones, and one point drawn from each combination of the halves [0, ½), [½, 1] of the
//   penalties' range (2⁵ for the bacterial five, 2⁴ for the four of eukaryotes); at each step size,
//   from 0.3 down to 0.001, every penalty and direction is tried in an order drawn at random, the
//   first that raises the correlation taken, and the step divided by 1.35 when none does.
//
// Limits. The wobble rules are those of bacterial and eukaryotic nuclear decoding; mitochondria
// read four-fold boxes with a single tRNA whose U₃₄ reads all four codons (superwobble), which
// Table 1 does not describe, so the tAI of a mitochondrial genome is not meaningful (the engine
// refuses it). Complexity: the weights O(64 · isoacceptors); tAI O(64) per gene; an S value O(G)
// per permutation; stAI O(starts · steps · G log G).

#pragma once

#include <genomics/codon/usage.hpp>

#include <tyche/descriptive/rank.hpp>
#include <tyche/descriptive/summary.hpp>

#include <calx/optimisation/nelder_mead.hpp>
#include <calx/random/stream.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <stdexcept>
#include <vector>

namespace genomics::codon {

inline constexpr std::uint32_t PERMUTE_W = 0x5045524d;  // "PERM": the S value's permutations
inline constexpr std::uint32_t NC_NULL = 0x4e434e55;    // "NCNU": genes simulated without selection
inline constexpr std::uint32_t STAI = 0x53544149;       // "STAI": the hill climbing's starts and orders

struct Isoacceptor {
  char amino_acid = 0;              // one letter, as the code tables write it
  std::array<Code, 3> anticodon{};  // 5′→3′: N₃₄ N₃₅ N₃₆
  double copies = 0;                // tGCN
  bool lysidine = false;            // a C₃₄ read as lysidine (tRNA-Ile with anticodon CAU)
};

inline std::vector<Isoacceptor> trna_pool(const TrnaGenes& trnas) {
  std::map<std::pair<char, int>, Isoacceptor> by;
  for (const auto& t : trnas.genes) {
    if (t.initiator) continue;
    const char aa = amino_acid_letter(t.amino_acid);
    if (aa == 0 || aa == 'U' || aa == 'O') continue;  // selenocysteine, pyrrolysine: recoded stops
    const int key = 16 * t.anticodon[0] + 4 * t.anticodon[1] + t.anticodon[2];
    auto& iso = by[{aa, key}];
    iso.amino_acid = aa;
    iso.anticodon = {t.anticodon[0], t.anticodon[1], t.anticodon[2]};
    iso.copies += 1;
    iso.lysidine = aa == 'I' && t.anticodon[0] == sequence::C && t.anticodon[1] == sequence::A && t.anticodon[2] == sequence::T;
  }
  std::vector<Isoacceptor> out;
  for (auto& [k, v] : by) out.push_back(v);
  return out;
}

// The penalties s of codon–anticodon pairings (dos Reis et al. 2004, Table 2).
struct Penalties {
  double IU = 0, GC = 0, UA = 0, CG = 0;                    // Watson–Crick, I:U
  double GU = 0.41, IC = 0.28, IA = 0.9999, UG = 0.68, LA = 0.89;
  std::array<double, 5> wobble() const { return {GU, IC, IA, UG, LA}; }
  void set_wobble(const std::array<double, 5>& s) { GU = s[0]; IC = s[1]; IA = s[2]; UG = s[3]; LA = s[4]; }
};

struct TrnaAdaptiveness {
  std::array<double, 64> W{}, w{};
  std::array<bool, 64> counted{};   // sense codons other than methionine's
  double w_mean = 0;                // the geometric mean that replaces w where W = 0
  std::size_t unread = 0;           // counted codons no tRNA reads (W = 0)
};

// The codons an isoacceptor reads, with their penalties.
inline std::vector<std::pair<int, double>> reads(const Isoacceptor& t, const GeneticCode& code, const Penalties& s) {
  using namespace sequence;
  const Code b1 = complement(t.anticodon[2]), b2 = complement(t.anticodon[1]);
  std::vector<std::pair<Code, double>> third;
  switch (t.anticodon[0]) {
    case A: third = {{T, s.IU}, {C, s.IC}, {A, s.IA}}; break;  // inosine
    case G: third = {{C, s.GC}, {T, s.GU}}; break;
    case T: third = {{A, s.UA}, {G, s.UG}}; break;
    case C:
      if (t.lysidine) third = {{A, s.LA}};
      else third = {{G, s.CG}};
      break;
    default: break;
  }
  std::vector<std::pair<int, double>> out;
  for (const auto& [x, penalty] : third) {
    const int c = codon_index(b1, b2, x);
    if (c >= 0 && code.amino_acids[static_cast<std::size_t>(c)] == t.amino_acid) out.push_back({c, penalty});
  }
  return out;
}

inline TrnaAdaptiveness trna_adaptiveness(const std::vector<Isoacceptor>& pool, const GeneticCode& code, const Penalties& s = {}) {
  TrnaAdaptiveness a;
  for (std::size_t c = 0; c < 64; c++) {
    const char aa = code.amino_acids[c];
    a.counted[c] = aa != '*' && aa != 'M';
  }
  for (const auto& t : pool)
    for (const auto& [c, penalty] : reads(t, code, s)) a.W[static_cast<std::size_t>(c)] += (1 - penalty) * t.copies;
  double top = 0;
  for (std::size_t c = 0; c < 64; c++)
    if (a.counted[c]) top = std::max(top, a.W[c]);
  if (!(top > 0)) throw std::invalid_argument("genomics::codon::trna_adaptiveness: no tRNA reads any counted codon");
  double logs = 0;
  std::size_t nonzero = 0;
  for (std::size_t c = 0; c < 64; c++)
    if (a.counted[c] && a.W[c] > 0) {
      a.w[c] = a.W[c] / top;
      logs += std::log(a.w[c]);
      nonzero++;
    }
  a.w_mean = std::exp(logs / static_cast<double>(nonzero));
  for (std::size_t c = 0; c < 64; c++)
    if (a.counted[c] && !(a.W[c] > 0)) { a.w[c] = a.w_mean; a.unread++; }
  return a;
}

inline double tai(const CodonCounts& n, const TrnaAdaptiveness& a) {
  double s = 0, L = 0;
  for (std::size_t c = 0; c < 64; c++)
    if (a.counted[c] && n[c] > 0) {
      s += n[c] * std::log(a.w[c]);
      L += n[c];
    }
  return L > 0 ? std::exp(s / L) : 0;
}

// The tAI of many genes from their log w: one pass over their counts.
namespace detail {
inline void tai_many(const std::vector<CodonCounts>& genes, const std::array<double, 64>& logw, const std::array<bool, 64>& counted, std::vector<double>& out) {
  out.resize(genes.size());
  for (std::size_t g = 0; g < genes.size(); g++) {
    double s = 0, L = 0;
    for (std::size_t c = 0; c < 64; c++)
      if (counted[c]) { s += genes[g][c] * logw[c]; L += genes[g][c]; }
    out[g] = L > 0 ? std::exp(s / L) : 0;
  }
}
}  // namespace detail

// ---- the expected effective number of codons ----------------------------------------------------------
inline double wright_expected_nc(double gc3s) { return 2 + gc3s + 29 / (gc3s * gc3s + (1 - gc3s) * (1 - gc3s)); }

struct NcCurve {
  double a = -6.0, b = 34.0, c = 1.025;  // dos Reis et al. (2004), Table 3, suggested
  double operator()(double x) const { return a + x + b / (x * x + (c - x) * (c - x)); }
};

struct NcCurveFit {
  NcCurve curve;
  double rss = 0;                       // Σ (f₁(x) − N_c)², eq. 9
  std::vector<double> gc3s, nc;         // the simulated genes
  int evaluations = 0;
  bool converged = false;
};

// Fits eq. 8 to points (x, N_c) by Nelder and Mead from Table 3's suggested constants (eq. 9).
inline NcCurveFit fit_nc_curve(std::vector<double> gc3s, std::vector<double> nc) {
  NcCurveFit fit;
  fit.gc3s = std::move(gc3s);
  fit.nc = std::move(nc);
  const auto rss = [&](std::span<const double> p) {
    const NcCurve f{p[0], p[1], p[2]};
    double s = 0;
    for (std::size_t i = 0; i < fit.nc.size(); i++) {
      const double e = f(fit.gc3s[i]) - fit.nc[i];
      s += e * e;
    }
    return s;
  };
  calx::optimisation::NelderMeadOptions o;
  o.max_evaluations = 20000;
  const auto r = calx::optimisation::nelder_mead(rss, {-6.0, 34.0, 1.025}, o);
  fit.curve = {r.x[0], r.x[1], r.x[2]};
  fit.rss = r.value;
  fit.evaluations = r.evaluations;
  fit.converged = r.converged;
  return fit;
}

// Wright's N_c when it is defined, with no cap; NaN when a class's mean F̂ is not positive (the
// sequence module's form then returns a value above 10¹¹, from its floor of 10⁻¹² on F̄).
inline double enc_wright_defined(const CodonCounts& n, const GeneticCode& code) {
  const double v = enc_wright(n, code, false);
  return v < 1e6 ? v : std::numeric_limits<double>::quiet_NaN();
}

// Simulates each gene of at least `min_codons` codons `replicates` times without selection (its
// amino acids kept, its codons drawn by a silent GC content x ~ U(0, 1)) and fits eq. 8 to those
// whose N_c is defined.
inline NcCurveFit fit_nc_curve(const std::vector<CodonCounts>& genes, const GeneticCode& code, std::size_t replicates = 3, std::uint64_t seed = 1,
                               std::size_t min_codons = 100) {
  const Families fam = synonymous_families(code, false);
  std::vector<double> xs, ys;
  for (std::size_t g = 0; g < genes.size(); g++) {
    double length = 0;
    for (std::size_t c = 0; c < 64; c++) length += fam.family[c] >= 0 ? genes[g][c] : 0;
    if (length < static_cast<double>(min_codons)) continue;
    for (std::size_t r = 0; r < replicates; r++) {
      calx::random::Stream rng(seed, NC_NULL, static_cast<std::uint32_t>(g), static_cast<std::uint32_t>(r));
      const double x = rng.uniform();
      CodonCounts sim{};
      for (const auto& cs : fam.codons) {
        double count = 0;
        for (int c : cs) count += genes[g][static_cast<std::size_t>(c)];
        if (!(count > 0)) continue;
        std::vector<double> p(cs.size());
        double total = 0;
        for (std::size_t k = 0; k < cs.size(); k++) {
          const Code third = sequence::codon_of(cs[k])[2];
          total += p[k] = (third == sequence::C || third == sequence::G) ? x : 1 - x;
        }
        for (long long draw = 0; draw < std::llround(count); draw++) {
          double u = rng.uniform() * total;
          std::size_t k = 0;
          while (k + 1 < cs.size() && u >= p[k]) u -= p[k++];
          sim[static_cast<std::size_t>(cs[k])] += 1;
        }
      }
      const double nc = enc_wright_defined(sim, code);
      if (!std::isfinite(nc)) continue;
      xs.push_back(third_positions(sim, code).gc3s);
      ys.push_back(nc);
    }
  }
  return fit_nc_curve(std::move(xs), std::move(ys));
}

// ---- the S value ------------------------------------------------------------------------------------
struct SValue {
  double s = 0;                   // cor(tAI, ψ)
  std::vector<double> tai, psi, nc, gc3s;
  double p_value = 1;             // by permutations of w over the counted codons
  std::size_t permutations = 0, as_extreme = 0;
};

inline SValue s_value(const std::vector<CodonCounts>& genes, const TrnaAdaptiveness& a, const GeneticCode& code, const NcCurve& curve = {}, std::size_t permutations = 0,
                      std::uint64_t seed = 1, calx::scheduler::ThreadPool* pool = nullptr) {
  SValue out;
  const std::size_t G = genes.size();
  out.tai.resize(G);
  out.psi.resize(G);
  out.nc.resize(G);
  out.gc3s.resize(G);
  for (std::size_t g = 0; g < G; g++) {
    out.tai[g] = tai(genes[g], a);
    out.nc[g] = enc_wright(genes[g], code, true);
    out.gc3s[g] = third_positions(genes[g], code).gc3s;
    out.psi[g] = curve(out.gc3s[g]) - out.nc[g];
  }
  out.s = tyche::descriptive::correlation(out.tai, out.psi);
  if (permutations == 0) return out;
  std::vector<std::size_t> idx;
  for (std::size_t c = 0; c < 64; c++)
    if (a.counted[c]) idx.push_back(c);
  std::vector<double> sp(permutations, 0);
  const auto run = [&](std::size_t b, std::size_t e) {
    std::vector<double> t;
    for (std::size_t k = b; k < e; k++) {
      calx::random::Stream rng(seed, PERMUTE_W, static_cast<std::uint32_t>(k), 0);
      std::vector<double> w;
      for (std::size_t c : idx) w.push_back(a.w[c]);
      for (std::size_t i = w.size(); i > 1; i--) std::swap(w[i - 1], w[rng.below(i)]);  // Fisher–Yates
      std::array<double, 64> logw{};
      for (std::size_t i = 0; i < idx.size(); i++) logw[idx[i]] = std::log(w[i]);
      detail::tai_many(genes, logw, a.counted, t);
      sp[k] = tyche::descriptive::correlation(t, out.psi);
    }
  };
  if (pool) calx::scheduler::parallel_for(*pool, permutations, 8, [&](std::size_t b, std::size_t e, std::size_t) { run(b, e); });
  else run(0, permutations);
  for (double v : sp) out.as_extreme += v >= out.s;
  out.permutations = permutations;
  out.p_value = (static_cast<double>(out.as_extreme) + 1) / (static_cast<double>(permutations) + 1);
  return out;
}

// ---- stAI ------------------------------------------------------------------------------------------
struct StaiOptions {
  bool prokaryote = true;          // five wobble penalties (with lysidine's) or four
  double first_step = 0.3, last_step = 0.001, shrink = 1.35;
  std::uint64_t seed = 1;
};

struct StaiResult {
  Penalties penalties;
  double correlation = 0;          // Spearman's ρ(tAI, target) at the optimum
  double start_correlation = 0;    // at dos Reis et al.'s penalties
  std::size_t starts = 0, evaluations = 0;
};

// Maximises Spearman's ρ between the genes' tAI and a target (DCBS, or measured expression).
inline StaiResult optimise_penalties(const std::vector<CodonCounts>& genes, const std::vector<double>& target, const std::vector<Isoacceptor>& pool,
                                     const GeneticCode& code, const StaiOptions& o = {}, calx::scheduler::ThreadPool* threads = nullptr) {
  if (genes.size() != target.size() || genes.size() < 3) throw std::invalid_argument("genomics::codon::optimise_penalties: a target for each gene, at least 3");
  const std::size_t d = o.prokaryote ? 5 : 4;
  const auto objective = [&](const std::array<double, 5>& s, std::vector<double>& t) {
    Penalties p;
    p.set_wobble(s);
    const auto a = trna_adaptiveness(pool, code, p);
    std::array<double, 64> logw{};
    for (std::size_t c = 0; c < 64; c++) logw[c] = a.counted[c] ? std::log(a.w[c]) : 0;
    detail::tai_many(genes, logw, a.counted, t);
    return tyche::descriptive::spearman(t, target);
  };
  // The starting points (Sabi and Tuller 2014, §2.3).
  std::vector<std::array<double, 5>> starts;
  starts.push_back(Penalties{}.wobble());
  starts.push_back({0, 0, 0, 0, 0});
  starts.push_back({1, 1, 1, 1, 1});
  calx::random::Stream draw(o.seed, STAI, 0xFFFFFFFFu, 0);
  for (std::size_t m = 0; m < (std::size_t{1} << d); m++) {
    std::array<double, 5> s = Penalties{}.wobble();
    for (std::size_t k = 0; k < d; k++) s[k] = ((m >> k) & 1 ? 0.5 : 0.0) + 0.5 * draw.uniform();
    starts.push_back(s);
  }
  if (!o.prokaryote)
    for (auto& s : starts) s[4] = Penalties{}.LA;
  struct Climb {
    std::array<double, 5> s{};
    double value = -2;
    std::size_t evaluations = 0;
  };
  std::vector<Climb> climbs(starts.size());
  const auto climb = [&](std::size_t k) {
    calx::random::Stream rng(o.seed, STAI, static_cast<std::uint32_t>(k), 1);
    std::vector<double> t;
    Climb c;
    c.s = starts[k];
    c.value = objective(c.s, t);
    c.evaluations = 1;
    std::vector<std::size_t> moves(2 * d);
    for (double h = o.first_step; h >= o.last_step;) {
      for (std::size_t i = 0; i < moves.size(); i++) moves[i] = i;
      for (std::size_t i = moves.size(); i > 1; i--) std::swap(moves[i - 1], moves[rng.below(i)]);
      bool improved = false;
      for (std::size_t m : moves) {
        auto s = c.s;
        const std::size_t which = m / 2;
        s[which] = std::clamp(s[which] + (m % 2 ? -h : h), 0.0, 1.0);
        if (s[which] == c.s[which]) continue;
        const double v = objective(s, t);
        c.evaluations++;
        if (v > c.value) {
          c.s = s;
          c.value = v;
          improved = true;
          break;
        }
      }
      if (!improved) h /= o.shrink;
    }
    climbs[k] = c;
  };
  if (threads) calx::scheduler::parallel_for(*threads, climbs.size(), 1, [&](std::size_t b, std::size_t e, std::size_t) { for (std::size_t k = b; k < e; k++) climb(k); });
  else
    for (std::size_t k = 0; k < climbs.size(); k++) climb(k);
  StaiResult r;
  std::size_t best = 0;
  for (std::size_t k = 0; k < climbs.size(); k++) {
    r.evaluations += climbs[k].evaluations;
    if (climbs[k].value > climbs[best].value) best = k;
  }
  r.penalties.set_wobble(climbs[best].s);
  r.correlation = climbs[best].value;
  std::vector<double> t;
  r.start_correlation = objective(Penalties{}.wobble(), t);
  r.starts = climbs.size();
  return r;
}

}  // namespace genomics::codon
