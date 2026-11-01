// The genotype of a structural variant from the reads that support each allele: SVTyper's
// Bayesian model (Chiang et al. 2015, "SpeedSeq: ultra-fast personal genome analysis and
// interpretation", Nature Methods 12, 966–968; its implementation svtyper at 3c5218a,
// singlesample.py and statistics.py, bayes_gt).
//
// - Support. The alternate allele is supported by the junction's evidence (cluster.hpp): its split
//   reads and gaps, each weighted by the probability that its parts are placed right (the product
//   of 1 − 10^(−Q/10) over them), and its discordant pairs, weighted by the product over both
//   mates. The reference allele at a breakend is supported by the mates whose alignment crosses
//   the reference junction there with at least `min_aligned` (20, SVTyper's) bases on each side,
//   and by the proper pairs that straddle it (one mate wholly on each side), each weighted alike and,
//   across a deletion or an insertion, by SVTyper's p_conc (pairs_across).
//   *Change.* SVTyper adds the reference reads of the two breakends of an event; a read at one
//   breakend and a read at the other are then two observations of the reference against one of
//   the single alternate junction, and the reference is counted twice. Here the reference support
//   of a junction is the mean over its two breakends, which is the depth of the reference allele
//   at the junction.
// - Likelihood (bayes_gt). With R and A the reference and alternate support, the genotypes 0/0,
//   0/1, 1/1 have log₁₀ L(g) = log₁₀ C(R + A, A) + A log₁₀ p_g + R log₁₀(1 − p_g), p_g the
//   fraction of alternate observations each genotype expects: SVTyper's 10⁻³, 1/2 and 0.9 for an
//   event that leaves one copy at its junction. *Change.* Here the fractions follow from the copies
//   at the junction: an adjacency made on h of the two haplotypes, where the reference adjacency
//   survives on r copies, expects h/(h + r); so a deletion, an inversion's or a translocation's
//   junction expects 1/2 (heterozygous) and 1 − 10⁻³ (homozygous), and a tandem duplication of k
//   extra copies, whose reference junctions survive on both haplotypes, expects k/(k + 2) and
//   k/(k + 1) (SVTyper takes 0.2 and 1/3 for any duplication). Real-valued supports are allowed
//   (C by log Γ).
// - Qualities as SVTyper's: GQ = −10 (log₁₀ L₂ − log₁₀ L₁) for the best and second genotypes,
//   held to 200; SQ = −10 log₁₀ P(0/0 | data) under a flat prior, the quality of the site.
//
// Complexity: O(log n + k) a breakend for n reads and k of them near it.

#pragma once

#include <genomics/align/types.hpp>
#include <genomics/structural/breakend.hpp>
#include <genomics/structural/pairs.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace genomics::structural {

struct SvGenotype {
  double reference = 0, alternate = 0;  // the weighted supports
  std::array<double, 3> log10_likelihood{};
  int genotype = 0;    // alternate copies: 0, 1, 2
  int gq = 0;
  double qual = 0;     // SQ
  double allele_fraction() const { return reference + alternate > 0 ? alternate / (reference + alternate) : 0; }
};

// SVTyper's bayes_gt with the expected fractions p = (ε, p_het, p_hom).
inline SvGenotype genotype_of(double reference, double alternate, double p_het, double p_hom, double error = 1e-3) {
  SvGenotype g;
  g.reference = reference;
  g.alternate = alternate;
  const double n = reference + alternate;
  const double combo = (std::lgamma(n + 1) - std::lgamma(alternate + 1) - std::lgamma(reference + 1)) / std::log(10.0);
  const std::array<double, 3> p{error, p_het, std::min(p_hom, 1 - error)};
  for (int k = 0; k < 3; k++) g.log10_likelihood[k] = combo + alternate * std::log10(p[k]) + reference * std::log10(1 - p[k]);
  std::array<int, 3> order{0, 1, 2};
  std::stable_sort(order.begin(), order.end(), [&](int x, int y) { return g.log10_likelihood[x] > g.log10_likelihood[y]; });  // ties as Python's sort
  g.genotype = order[0];
  g.gq = static_cast<int>(std::min(200.0, -10 * (g.log10_likelihood[order[1]] - g.log10_likelihood[order[0]])));
  const double top = g.log10_likelihood[order[0]];
  double s = 0;
  for (double l : g.log10_likelihood) s += std::pow(10.0, l - top);
  const double total = top + std::log10(s);
  g.qual = std::fabs(-10 * (g.log10_likelihood[0] - total));
  return g;
}

// The mapped mates and the proper pairs, by position, for counting the reference's support.
class ReadIndex {
 public:
  ReadIndex(const MateMappings& mates, const std::vector<PairRecord>& pairs, const Library& lib, std::size_t contigs, int min_mapq = 20, std::size_t min_gap = 10)
      : reads_(contigs), spans_(contigs), median_(lib.median), sd_(std::max(1.0, lib.sd())) {
    // Each mate's aligned blocks between its gaps of at least `min_gap` bases: a mate carrying a
    // deletion or an insertion inside its alignment does not cross the reference there.
    for (const auto& ms : mates)
      for (const auto& m : ms) {
        if (!m.primary || m.mapq < min_mapq) continue;
        const double w = placement_probability(m.mapq);
        std::size_t j = m.alignment.target_begin, start = j;
        for (const auto& op : m.alignment.cigar) {
          if ((op.op == 'D' || op.op == 'I') && op.length >= min_gap) {
            if (j > start) reads_[m.contig].push_back({start, j, false, false, w});
            start = j + (op.op == 'D' ? op.length : 0);
          }
          if (align::consumes_target(op.op)) j += op.length;
        }
        if (m.alignment.target_end > start) reads_[m.contig].push_back({start, m.alignment.target_end, false, false, w});
      }
    for (const auto& p : pairs) {
      if (p.type != PairClass::proper || !p.normal) continue;
      const Mapping& f = p.mate[0]->reverse ? *p.mate[1] : *p.mate[0];
      const Mapping& r = p.mate[0]->reverse ? *p.mate[0] : *p.mate[1];
      spans_[f.contig].push_back({f.target_end, r.target_begin, false, false, placement_probability(f.mapq) * placement_probability(r.mapq), static_cast<double>(p.length)});
    }
    for (auto& v : reads_) std::sort(v.begin(), v.end(), [](const Span& x, const Span& y) { return x.begin < y.begin; });
    for (auto& v : spans_) std::sort(v.begin(), v.end(), [](const Span& x, const Span& y) { return x.begin < y.begin; });
    for (const auto& v : reads_)
      for (const auto& s : v) longest_ = std::max(longest_, s.end - s.begin);
    for (const auto& v : spans_)
      for (const auto& s : v) longest_span_ = std::max(longest_span_, s.end > s.begin ? s.end - s.begin : 0);
  }

  // The boundary between bases J − 1 and J that a breakend cuts.
  static std::size_t boundary(const Breakend& e) { return e.plus ? e.position + 1 : e.position; }

  // Σ of the placement probabilities of the mates crossing the boundary with `min_aligned` bases each side.
  double reads_across(const Breakend& e, std::size_t min_aligned) const {
    const std::size_t J = boundary(e);
    double s = 0;
    for_each(reads_[e.contig], longest_, J, [&](const Span& x) {
      if (x.begin + min_aligned <= J && x.end >= J + min_aligned) s += x.weight;
    });
    return s;
  }
  // Σ over the proper pairs with one mate wholly on each side of the boundary, each weighted by
  // SVTyper's p_conc when the alternate allele changes a template by `delta` bases (+L across a
  // deletion of L, −L across an insertion): P(ref | t) = f(t)/(f(t) + f(t − δ)) for the template t
  // seen on the reference, f the library's law (normal, of its median and 1.4826 MAD).
  double pairs_across(const Breakend& e, double delta = 0) const {
    const std::size_t J = boundary(e);
    double s = 0;
    for_each(spans_[e.contig], longest_span_, J, [&](const Span& x) {
      if (x.begin > J || x.end < J) return;
      if (delta == 0) { s += x.weight; return; }
      const auto f = [&](double t) { const double z = (t - median_) / sd_; return std::exp(-0.5 * z * z); };
      const double fr = f(x.length), fa = f(x.length - delta);
      s += fr + fa > 0 ? x.weight * fr / (fr + fa) : 0;
    });
    return s;
  }
  // The weighted depth of mates over [b, e) of a contig (their aligned bases).
  double depth(std::uint32_t contig, std::size_t b, std::size_t e) const {
    double s = 0;
    for_each(reads_[contig], longest_, b, [&](const Span& x) { if (x.end > b && x.begin < e) s += x.weight; }, e);
    return s;
  }

 private:
  struct Span {
    std::size_t begin = 0, end = 0;
    bool clip_left = false, clip_right = false;
    double weight = 1;
    double length = 0;  // a pair's template
  };
  // f(x) for every span that may contain [J, to): begin in [J − longest, to].
  template <class F>
  static void for_each(const std::vector<Span>& v, std::size_t longest, std::size_t J, F&& f, std::size_t to = SIZE_MAX) {
    const std::size_t from = J > longest ? J - longest : 0;
    const std::size_t last = to == SIZE_MAX ? J : to;
    auto it = std::lower_bound(v.begin(), v.end(), from, [](const Span& s, std::size_t x) { return s.begin < x; });
    for (; it != v.end() && it->begin <= last; ++it) f(*it);
  }
  std::vector<std::vector<Span>> reads_, spans_;
  std::size_t longest_ = 0, longest_span_ = 0;
  double median_ = 0, sd_ = 1;
};

struct GenotypeOptions {
  std::size_t min_aligned = 20;  // SVTyper's
  double error = 1e-3;           // p of an alternate observation at 0/0, SVTyper's
};

// The reference support of a junction: the mean over its two breakends of the mates and proper
// pairs crossing the reference there (for a novel insertion, whose breakends are one site, that site).
inline double reference_support(const ReadIndex& index, const Adjacency& j, const GenotypeOptions& o = {}, double delta = 0) {
  const double a = index.reads_across(j.a, o.min_aligned) + index.pairs_across(j.a, delta);
  const double b = index.reads_across(j.b, o.min_aligned) + index.pairs_across(j.b, delta);
  const bool one_site = j.a.contig == j.b.contig && ReadIndex::boundary(j.a) == ReadIndex::boundary(j.b);
  return one_site ? a : (a + b) / 2;
}

}  // namespace genomics::structural
