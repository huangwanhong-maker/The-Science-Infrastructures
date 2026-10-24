// Checks of motifs: the information of columns whose answer is known (a single base 2 bits,
// two bases 1 bit, all four 0); Schneider's small-sample correction exact where it can be done
// by hand (e(1) = 2, e(2) = 5/4) and agreeing with Miller's first-order formula at n = 150; the
// distribution of PWM scores by the lattice equal, point for point, to the enumeration of every
// word of short PWMs, and the unrounded tails within the bounds the rounding allows; scanning
// finds a site on the reverse strand; and the Gibbs sampler and MEME's EM (OOPS, one strand and
// both) recover a motif planted in random sequences at the right positions in most of them, the
// sampler the same on one thread and on four; EM never lowers its log posterior, and its ZOOPS fit
// finds every planted site at a likelihood at least that of the planted truth.

#include <genomics/annotation/motifs.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::annotation;

namespace {

std::string fixed(double v, int digits = 3) {
  char b[32];
  std::snprintf(b, sizeof b, "%.*f", digits, v);
  return b;
}

struct Planted {
  std::vector<Codes> sequences;
  std::vector<std::size_t> position;  // SIZE_MAX where none was planted
  std::vector<bool> reverse;
};

// Random sequences with a motif planted: each site keeps the consensus base with probability
// `keep` and takes another base uniformly otherwise.
Planted plant(std::size_t count, std::size_t length, const Codes& consensus, double keep, double fraction, bool both, std::uint64_t seed) {
  Planted p;
  calx::random::Stream r(seed, 0x706c6e74, 0, 0);
  const std::size_t W = consensus.size();
  for (std::size_t s = 0; s < count; s++) {
    Codes x = simulate::random_sequence(length, 0.5, seed, static_cast<std::uint32_t>(s));
    const bool here = r.uniform() < fraction;
    std::size_t at = SIZE_MAX;
    bool rev = false;
    if (here) {
      at = r.below(length - W + 1);
      rev = both && r.bernoulli(0.5);
      Codes site = consensus;
      for (auto& c : site)
        if (r.uniform() >= keep) c = simulate::substitute(r, c);
      if (rev) site = sequence::reverse_complement(site);
      std::copy(site.begin(), site.end(), x.begin() + static_cast<std::ptrdiff_t>(at));
    }
    p.sequences.push_back(std::move(x));
    p.position.push_back(at);
    p.reverse.push_back(rev);
  }
  return p;
}

// The fraction of planted sites found at their place (and strand), and of reported sites that are planted ones.
std::pair<double, double> recovered(const Planted& p, const Motif& m) {
  double planted = 0, right = 0;
  for (std::size_t s = 0; s < p.sequences.size(); s++) planted += p.position[s] != SIZE_MAX;
  for (const auto& site : m.sites)
    right += p.position[site.sequence] == site.position && p.reverse[site.sequence] == site.reverse;
  return {right / planted, m.sites.empty() ? 0 : right / static_cast<double>(m.sites.size())};
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x6d6f7466, 0, 0);

  // ---- information ------------------------------------------------------------------------------------------
  {
    // Columns: all A; A and C equally; all four equally (counts 4 each), without pseudocounts.
    const Pwm m = pwm_from_counts({8, 0, 0, 0, 4, 4, 0, 0, 2, 2, 2, 2}, 0);
    const auto ic = information_content(m);
    s.check("a column of a single base carries 2 bits, two equal bases 1 bit, four 0 bits (uniform background)", std::fabs(ic[0] - 2) + std::fabs(ic[1] - 1) + std::fabs(ic[2]), 1e-15, "bits");
    s.require("a single-base column's weight is 2 bits for its base and −∞ for the others", m.w(0, sequence::A) == 2 && m.w(0, sequence::C) == MOTIF_NEG_INF);
    s.check("Schneider's e(1) = 2: one base never has entropy", std::fabs(small_sample_correction(1) - 2), 1e-12, "bits");
    s.check("e(2) = 5/4: two bases differ with probability 3/4, and then have 1 bit", std::fabs(small_sample_correction(2) - 1.25), 1e-12, "bits");
    const double exact = small_sample_correction(150), miller = 3.0 / (2 * 150 * std::numbers::ln2);
    s.check("e(150) by the sum over compositions agrees with Miller's (s − 1)/(2 n ln 2) (relative)", std::fabs(exact - miller) / miller, 0.01);
    const auto lg = logo(m, false);
    s.require("the logo of the single-base column is one letter 2 bits tall; of the uniform one, nothing", std::fabs(lg[0].height[0] - 2) < 1e-15 && lg[0].height[1] == 0 && lg[2].information == 0);
    const auto corrected = logo(m, true);
    s.check("with the correction, eight sites of one base carry 2 − e(8) bits", std::fabs(corrected[0].information - (2 - small_sample_correction(8))), 1e-12, "bits");
  }

  // ---- the distribution of scores ------------------------------------------------------------------------------
  {
    double worst_mass = 0, worst_total = 0;
    bool inside = true, thresholds = true;
    for (int rep = 0; rep < 4; rep++) {
      const std::size_t W = rep < 2 ? 6 : 8;
      std::vector<double> counts(4 * W);
      for (auto& c : counts) c = static_cast<double>(r.below(12));
      const Background q = rep % 2 ? Background{0.3, 0.2, 0.2, 0.3} : UNIFORM;
      const Pwm m = pwm_from_counts(counts, 1, q);
      const ScoreDistribution d = score_distribution(m, rep < 2 ? 0.01 : 0.002);
      // Enumerate every word.
      std::vector<double> mass(d.mass.size(), 0), scores, probs;
      const std::size_t words = std::size_t{1} << (2 * W);
      Codes w(W);
      for (std::size_t k = 0; k < words; k++) {
        double p = 1;
        for (std::size_t i = 0; i < W; i++) {
          w[i] = static_cast<Code>((k >> (2 * i)) & 3);
          p *= q[w[i]];
        }
        const std::int64_t lat = d.rounded(w, 0);
        mass[static_cast<std::size_t>(lat - d.lowest)] += p;
        scores.push_back(m.score(w, 0));
        probs.push_back(p);
      }
      double total = 0;
      for (std::size_t j = 0; j < mass.size(); j++) {
        worst_mass = std::max(worst_mass, std::fabs(mass[j] - d.mass[j]));
        total += d.mass[j];
      }
      worst_total = std::max(worst_total, std::fabs(total - 1));
      // The unrounded tail P(S ≥ t) at thresholds across the range, within the rounding's bounds.
      for (int k = 0; k <= 20; k++) {
        const double t = m.min_score() + (m.max_score() - m.min_score()) * k / 20.0;
        double tail = 0;
        for (std::size_t i = 0; i < scores.size(); i++)
          if (scores[i] >= t) tail += probs[i];
        const auto [lo, hi] = d.bounds(t);
        inside = inside && lo <= tail + 1e-14 && tail <= hi + 1e-14;
      }
      for (double p : {1e-1, 1e-2, 1e-3}) {
        const double t = d.threshold(p);
        thresholds = thresholds && d.tail(t) <= p && d.tail(t - d.granularity) > p;
      }
    }
    s.check("the lattice distribution equals the enumeration of every word, point by point (widths 6 and 8, two backgrounds)", worst_mass, 1e-15);
    s.check("its mass sums to one", worst_total, 1e-13);
    s.require("every unrounded tail P(S ≥ t) lies within the bounds P(S' ≥ t ± Wε/2)", inside);
    s.require("the threshold of a p-value is the least lattice score whose tail is at most it", thresholds);
  }

  // ---- scanning ------------------------------------------------------------------------------------------------------
  {
    const Codes consensus = sequence::encode_dna("TTGACAGCTAGC");
    std::vector<Codes> sites(20, consensus);
    const Pwm m = pwm_from_sites(sites, 1);
    Codes x = simulate::random_sequence(500, 0.5, 3);
    const Codes rc = sequence::reverse_complement(consensus);
    std::copy(rc.begin(), rc.end(), x.begin() + 300);
    const ScoreDistribution d = score_distribution(m);
    const auto hits = scan(m, x, m.max_score() - 1e-9, true, &d);
    s.require("scanning both strands finds the consensus planted reversed at 300, with the largest score and its p-value 4^−12 under the rounded matrix",
              hits.size() == 1 && hits[0].position == 300 && hits[0].reverse && std::fabs(hits[0].score - m.max_score()) < 1e-9 && std::fabs(hits[0].p_value - std::pow(0.25, 12)) < 1e-15,
              hits.empty() ? "no hit" : std::to_string(hits[0].position) + (hits[0].reverse ? " reverse" : " forward") + ", p " + std::to_string(hits[0].p_value));
  }

  // ---- discovery --------------------------------------------------------------------------------------------------------
  {
    const Codes consensus = sequence::encode_dna("GATTACATGCCGTA");
    const Planted oops = plant(25, 200, consensus, 0.9, 1.0, false, 11);
    GibbsOptions go;
    go.width = consensus.size();
    go.seed = 5;
    const Motif g = gibbs_motif(oops.sequences, go);
    const auto [gs, gp] = recovered(oops, g);
    s.require("the Gibbs sampler recovers a planted motif at its exact place in at least 80% of 25 sequences", gs >= 0.8, fixed(gs) + " (F = " + fixed(g.score, 1) + " bits)");
    calx::scheduler::ThreadPool pool(4);
    const Motif g4 = gibbs_motif(oops.sequences, go, &pool);
    bool same = g4.score == g.score && g4.sites.size() == g.sites.size();
    for (std::size_t k = 0; same && k < g.sites.size(); k++) same = g4.sites[k].position == g.sites[k].position;
    s.require("the Gibbs sampler gives the same motif on four threads as on one", same);

    MemeOptions mo;
    mo.width = consensus.size();
    const Motif e = meme_motif(oops.sequences, mo, &pool);
    const auto [es, ep] = recovered(oops, e);
    s.require("MEME's EM (OOPS) recovers it at its exact place in at least 80% of the sequences", es >= 0.8, fixed(es) + " after " + std::to_string(e.iterations) + " iterations");

    double drop = 0;
    for (std::size_t i = 1; i < e.trace.size(); i++) drop = std::max(drop, e.trace[i - 1] - e.trace[i]);
    s.check("MEME's EM never lowers the log posterior (OOPS, " + std::to_string(e.trace.size()) + " iterations)", drop, 1e-9, "bits");

    const Planted zoops = plant(30, 200, consensus, 0.9, 0.6, false, 12);
    std::size_t planted = 0;
    std::vector<Codes> true_sites;
    for (std::size_t k = 0; k < zoops.sequences.size(); k++)
      if (zoops.position[k] != SIZE_MAX) {
        planted++;
        true_sites.emplace_back(zoops.sequences[k].begin() + static_cast<std::ptrdiff_t>(zoops.position[k]), zoops.sequences[k].begin() + static_cast<std::ptrdiff_t>(zoops.position[k] + consensus.size()));
      }
    mo.model = MotifModel::zoops;
    const Motif z = meme_motif(zoops.sequences, mo, &pool);
    const auto [zs, zp] = recovered(zoops, z);
    drop = 0;
    for (std::size_t i = 1; i < z.trace.size(); i++) drop = std::max(drop, z.trace[i - 1] - z.trace[i]);
    s.check("MEME's EM never lowers the log posterior (ZOOPS)", drop, 1e-9, "bits");
    // The planted truth: the matrix of the true sites, and the true fraction.
    const double frac = static_cast<double>(planted) / 30.0;
    const double truth = motif_log_likelihood(zoops.sequences, pwm_from_sites(true_sites, 1, base_frequencies(zoops.sequences)), MotifModel::zoops, frac);
    s.require("MEME's EM (ZOOPS) on sequences of which some lack the motif finds every planted site, at a likelihood at least that of the planted motif", zs >= 0.95 && z.score >= truth,
              "Sn " + fixed(zs) + "; " + fixed(z.score, 2) + " bits against the truth's " + fixed(truth, 2) + "; precision " + fixed(zp) + " (" + std::to_string(z.sites.size()) +
                  " sites for " + std::to_string(planted) + " planted), γ " + fixed(z.gamma) + " for " + fixed(frac) + ": the maximum-likelihood fit takes weak sites where none was planted");

    const Planted both = plant(25, 200, consensus, 0.9, 1.0, true, 13);
    go.both_strands = true;
    const Motif gb = gibbs_motif(both.sequences, go, &pool);
    mo.model = MotifModel::oops;
    mo.both_strands = true;
    const Motif eb = meme_motif(both.sequences, mo, &pool);
    // A motif found on both strands may be reported as the reverse complement of the planted one,
    // with every strand flipped: either reading counts.
    const auto either = [&](const Motif& m) {
      double same_strand = 0, flipped = 0;
      for (const auto& site : m.sites) {
        same_strand += both.position[site.sequence] == site.position && both.reverse[site.sequence] == site.reverse;
        flipped += both.position[site.sequence] == site.position && both.reverse[site.sequence] != site.reverse;
      }
      return std::max(same_strand, flipped) / 25.0;
    };
    s.require("on both strands, the Gibbs sampler and MEME place at least 80% of the sites, strand included", either(gb) >= 0.8 && either(eb) >= 0.8,
              "Gibbs " + fixed(either(gb)) + ", MEME " + fixed(either(eb)));
  }
  return s.done();
}
