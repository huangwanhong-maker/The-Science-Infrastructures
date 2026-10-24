// Checks of the CpG-island model: the chains of Durbin et al. (1998, p. 50) are distributions,
// their log-odds is the book's table β (p. 51) to its rounding, the eight-state model is a valid
// HMM whose labels have the sojourns asked for, a chain is re-estimated from its own samples, and
// islands planted in a long sequence drawn from the model are found by Viterbi and by posterior
// decoding, with their sensitivity and specificity at the level of bases and of islands.

#include <genomics/hmm/cpg.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::hmm;

namespace {

std::string fixed(double v, int digits = 3) {
  char b[32];
  std::snprintf(b, sizeof b, "%.*f", digits, v);
  return b;
}

struct Score {
  double sensitivity = 0, specificity = 0, found = 0, called_true = 0;
};

// Base level: Sn = TP/(TP + FN), Sp = TP/(TP + FP) (Burset and Guigó 1996, Genomics 34, 353–367).
// Island level: the fraction of planted islands overlapped by a call, and of calls overlapping a
// planted island.
Score score(const std::vector<bool>& truth, const std::vector<Island>& calls) {
  std::vector<bool> called(truth.size(), false);
  for (const auto& c : calls)
    for (std::size_t i = c.begin; i < c.end; i++) called[i] = true;
  double tp = 0, fp = 0, fn = 0;
  for (std::size_t i = 0; i < truth.size(); i++) {
    tp += truth[i] && called[i];
    fp += !truth[i] && called[i];
    fn += truth[i] && !called[i];
  }
  Score s;
  s.sensitivity = tp / (tp + fn);
  s.specificity = tp / (tp + fp);
  // Planted islands as runs of truth.
  double planted = 0, hit = 0;
  for (std::size_t i = 0; i < truth.size();) {
    if (!truth[i]) { i++; continue; }
    std::size_t j = i;
    bool any = false;
    for (; j < truth.size() && truth[j]; j++) any = any || called[j];
    planted++;
    hit += any;
    i = j;
  }
  double good = 0;
  for (const auto& c : calls) {
    bool any = false;
    for (std::size_t i = c.begin; i < c.end && !any; i++) any = truth[i];
    good += any;
  }
  s.found = hit / planted;
  s.called_true = calls.empty() ? 0 : good / static_cast<double>(calls.size());
  return s;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- the chains and their log-odds ------------------------------------------------------------------
  {
    double worst = 0;
    for (const auto* c : {&CPG_PLUS, &CPG_MINUS})
      for (int row = 0; row < 4; row++) {
        double sum = 0;
        for (int t = 0; t < 4; t++) sum += c->a[4 * row + t];
        worst = std::max(worst, std::fabs(sum - 1));
      }
    s.check("the + and − chains of Durbin et al. (p. 50): every row sums to one within its rounding", worst, 0.0015);
    // The book's table of β = log₂(a⁺/a⁻) in bits, Durbin et al. 1998, p. 51 (computed by the
    // authors from unrounded counts; transcribed by the Stanford Stats 366 notes,
    // web.stanford.edu/class/stats366/exs/HMM1.html, retrieved 2026-10-24, whose G→A reads
    // −0.623 where other transcriptions read −0.624).
    const double book[16] = {-0.740, 0.419, 0.580, -0.803,  //
                             -0.913, 0.302, 1.812, -0.685,  //
                             -0.623, 0.461, 0.331, -0.730,  //
                             -1.169, 0.573, 0.393, -0.679};
    const auto beta = log_odds_table();
    double d = 0;
    for (int i = 0; i < 16; i++) d = std::max(d, std::fabs(beta[i] - book[i]));
    s.check("the log-odds of the chains is the book's table β (p. 51), to the rounding of the chains' three decimals", d, 0.006, "bits");
    s.require("CG is the most island-like pair (β_CG = 1.812 bits) and TA the least", std::max_element(beta.begin(), beta.end()) - beta.begin() == 4 * 1 + 2 &&
                                                                                    std::min_element(beta.begin(), beta.end()) - beta.begin() == 4 * 3 + 0);
    // A sequence's log-odds is the sum of its pairs': CGCG has CG, GC, CG.
    const Codes x = sequence::encode_dna("CGCGNTA");
    s.check("the log-odds of CGCGNTA is β_CG + β_GC + β_CG + β_TA (pairs with N skipped)", std::fabs(log_odds(x) - (2 * beta[6] + beta[9] + beta[12])), 1e-12, "bits");
  }

  // ---- the eight-state model --------------------------------------------------------------------------
  {
    const Model m = cpg_model(1000, 100000);
    bool valid = true;
    try { m.validate(1e-12); } catch (...) { valid = false; }
    s.require("the eight-state model is a valid HMM, each state emitting its own base", valid && m.e(1, 1) == 1 && m.e(5, 1) == 1 && m.e(5, 2) == 0);
    double stay_plus = 0, stay_minus = 0;
    for (int t = 0; t < 4; t++) { stay_plus += m.a(0, t); stay_minus += m.a(4, 4 + t); }
    s.check("leaving an island has probability 1/1000 and entering one 1/100 000", std::fabs(1 - stay_plus - 1e-3) + std::fabs(1 - stay_minus - 1e-5), 1e-15);
  }

  // ---- re-estimating a chain from its samples ----------------------------------------------------------------
  {
    Model chain(4, 4);
    const auto st = CPG_PLUS.normalised().stationary();
    for (int k = 0; k < 4; k++) {
      chain.initial[k] = st[k];
      chain.e(k, k) = 1;
      for (int t = 0; t < 4; t++) chain.a(k, t) = CPG_PLUS.normalised()(k, t);
    }
    calx::random::Stream r(31, SAMPLE, 0, 0);
    const Codes x = sample(chain, 200000, r).symbols;
    const MarkovChain e = estimate_chain({x});
    double worst = 0;
    std::array<double, 4> from{};
    for (std::size_t i = 0; i + 1 < x.size(); i++) from[x[i]] += 1;
    for (int k = 0; k < 4; k++)
      for (int t = 0; t < 4; t++) {
        const double p = CPG_PLUS.normalised()(k, t);
        worst = std::max(worst, std::fabs(e(k, t) - p) / std::sqrt(p * (1 - p) / from[k]));
      }
    s.check("the maximum-likelihood chain of 200 000 bases drawn from the + chain is the + chain", worst, 5, "standard errors");
  }

  // ---- planted islands ------------------------------------------------------------------------------------
  {
    // A sequence drawn from the model itself, islands of mean 1 kb every 20 kb or so: the truth is
    // the labels of its path.
    const Model m = cpg_model(1000, 20000);
    calx::random::Stream r(41, SAMPLE, 0, 0);
    const Sample smp = sample(m, 600'000, r);
    std::vector<bool> truth(smp.states.size());
    for (std::size_t i = 0; i < truth.size(); i++) truth[i] = smp.states[i] < 4;
    IslandOptions o;
    o.mean_island = 1000;
    o.mean_background = 20000;
    const auto viterbi_calls = cpg_islands(smp.symbols, o);
    o.posterior = true;
    const auto posterior_calls = cpg_islands(smp.symbols, o);
    const Score v = score(truth, viterbi_calls), p = score(truth, posterior_calls);
    std::size_t planted = 0;
    for (std::size_t i = 0; i < truth.size(); i++) planted += truth[i] && (i == 0 || !truth[i - 1]);
    const std::string vd = "Sn " + fixed(v.sensitivity) + ", Sp " + fixed(v.specificity) + "; islands found " + fixed(v.found) + ", calls true " + fixed(v.called_true) +
                           " (" + std::to_string(viterbi_calls.size()) + " calls, " + std::to_string(planted) + " planted)";
    const std::string pd = "Sn " + fixed(p.sensitivity) + ", Sp " + fixed(p.specificity) + "; islands found " + fixed(p.found) + ", calls true " + fixed(p.called_true) +
                           " (" + std::to_string(posterior_calls.size()) + " calls)";
    s.require("Viterbi decoding finds planted islands in 600 kb, base-level sensitivity and specificity at least 0.85", v.sensitivity >= 0.85 && v.specificity >= 0.85, vd);
    s.require("posterior decoding finds them too, at least 0.85 each", p.sensitivity >= 0.85 && p.specificity >= 0.85, pd);
    double inside = 0, outside = 0, ni = 0, no = 0;
    for (const auto& w : log_odds_windows(smp.symbols, 200, 200)) {
      std::size_t t = 0;
      for (std::size_t i = w.begin; i < w.end; i++) t += truth[i];
      if (t == w.end - w.begin) { inside += w.value; ni++; }
      if (t == 0) { outside += w.value; no++; }
    }
    s.require("windows of 200 bases inside islands score positive log-odds per base on average, and outside negative",
              inside / ni > 0 && outside / no < 0, fixed(inside / ni) + " and " + fixed(outside / no) + " bits per base");
    bool stats = !viterbi_calls.empty();
    for (const auto& c : viterbi_calls) stats = stats && c.end > c.begin && c.gc >= 0 && c.gc <= 1;
    s.require("each island carries its interval, log-odds, GC and CpG o/e", stats);
  }
  return s.done();
}
