// The evidence of junctions in mapped pairs: discordant pairs, split reads, long gaps inside one
// alignment, and soft clips; each piece of evidence two intervals of positions with a probability
// on each, LUMPY's breakpoint (Layer et al. 2014, Genome Biology 15, R84, "Breakpoint
// representation": b = ⟨E, l, r, v⟩ with probability vectors over intervals l and r and a variety
// v, here the adjacency's contigs and orientations).
//
// - A discordant pair (pairs.hpp) points from each mate to a breakend: a mate mapped forward at
//   unclipped start s to a + breakend at some i ≥ its last aligned base, a mate mapped reverse
//   with unclipped end e to a − breakend at some i ≤ its first aligned base. Were the junction at
//   (i, j), the fragment would have length (i − s + 1) + (e′ − j) for the other mate's unclipped
//   end e′ (or the same terms mirrored), so the probability of a placement is the probability of
//   that length. LUMPY's paired-end module gives each side the marginal P(i) ∝ P(F > i − s)
//   ("Paired-end module": the proportion of fragments longer than the distance); here the other
//   mate's whole length ℓ′ is added, P(i) ∝ S(i − s + 1 + ℓ′) with S(x) = P(F ≥ x) the library's
//   empirical survival (*Change*: the other mate occupies at least its own length of the
//   fragment), over the positions where S is above its tail (10⁻³), and back into the mate by
//   LUMPY's back distance (20 bases, its simulations' setting) at the probability of the mate's
//   end, since a junction's microhomology is aligned on both of its sides.
// - A split read (a mate with a supplementary part) gives the junction exactly: its parts, in the
//   read's order, leave the reference at the end of one and enter it at the start of the next; an
//   overlap of the two parts on the read is microhomology (taken by the first part), a gap between
//   them novel bases. The adjacency is written canonically (breakend.hpp), and its probability is
//   LUMPY's split-read shape: a peak at the junction falling exponentially (here e^(−|d|/2) over
//   ±`split_slop` bases).
// - A gap of at least `min_gap` (50) bases inside one alignment, a deletion D or an insertion I of
//   the CIGAR, is a junction of a deletion or a novel insertion seen whole (as Sniffles reads long
//   reads' alignments; Smolka et al. 2024, Nature Biotechnology 42, 1571–1580), exact and with its
//   bases.
// - A soft clip of at least `min_clip` (25, DELLY's -c) bases at an end of a primary alignment not
//   explained by a supplementary part is a breakend without a partner, exact: the clipped bases
//   are the other side of the junction (read in the derivative's direction), which is how
//   insertions longer than a read, and junctions whose other side maps nowhere, are seen.
//
// Each piece carries the weight of its placement: SVTyper's prob_mapq, 1 − 10^(−Q/10), multiplied
// over the mapped parts (Chiang et al. 2015, Nature Methods 12, 966–968; svtyper/singlesample.py).
// Complexity O(N + Σ ℓ) for N pairs and the CIGARs' lengths.

#pragma once

#include <genomics/align/types.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/structural/breakend.hpp>
#include <genomics/structural/pairs.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace genomics::structural {

enum class Source : std::uint8_t { pair, split, gap };

inline const char* source_name(Source s) { return s == Source::pair ? "pair" : s == Source::split ? "split" : "gap"; }

// An interval of positions with relative log probabilities.
struct Interval {
  std::size_t begin = 0;
  std::vector<double> log_p;  // position begin + k has log probability log_p[k] (unnormalised)
  std::size_t end() const { return begin + log_p.size(); }
  bool contains(std::size_t x) const { return x >= begin && x < end(); }
};

struct Evidence {
  Source source = Source::pair;
  Adjacency junction;      // the most likely junction (exact for splits and gaps)
  Interval a, b;           // over the positions of junction.a and junction.b
  std::uint32_t read = 0;  // the pair (pair evidence) or the mate (split, gap)
  double weight = 1;       // the probability that its placements are right
};

// A soft clip: an exact breakend and the bases on the junction's other side.
struct Clip {
  Breakend at;
  std::string bases;   // in the derivative's direction away from `at` (for −: the bases before it, ending at the junction)
  std::uint32_t read = 0;
  int mapq = 0;
};

struct EvidenceOptions {
  int min_mapq = 20;
  std::size_t min_gap = 50;     // the shortest D or I of a CIGAR taken as a junction
  std::size_t min_clip = 25;    // DELLY's minimum clipping length
  std::size_t min_part = 20;    // the shortest part of a split read
  std::size_t split_slop = 10;  // half-width of a split read's interval
};

namespace detail {

// The interval of a mate's breakend: the side it points to, from its mapping and the other mate's length.
inline Interval mate_interval(const Mapping& m, std::int64_t other_length, const Library& lib, std::size_t contig_length, Breakend& e, std::int64_t back = 20) {
  Interval v;
  e.contig = m.contig;
  e.plus = !m.reverse;
  const std::int64_t n = static_cast<std::int64_t>(contig_length);
  if (e.plus) {
    const std::int64_t s = unclipped_begin(m), first = static_cast<std::int64_t>(m.target_end) - 1;
    // LUMPY's back distance: the junction may lie inside the mate's last bases (microhomology is
    // aligned on both of a junction's sides), at the probability of the mate's end.
    const std::int64_t lo = std::max<std::int64_t>(0, first - back);
    v.begin = static_cast<std::size_t>(lo);
    const double at_end = lib.survival_at(first - s + 1 + other_length);
    if (at_end <= 0) return v;
    for (std::int64_t i = lo; i < first; i++) v.log_p.push_back(std::log(at_end));
    std::int64_t i = std::max<std::int64_t>(first, 0);
    for (; i < n; i++) {
      const double p = lib.survival_at(i - s + 1 + other_length);
      if (p <= 0) break;
      v.log_p.push_back(std::log(p));
    }
  } else {
    const std::int64_t e2 = unclipped_end(m), first = static_cast<std::int64_t>(m.target_begin);
    const double at_end = lib.survival_at(e2 - first + other_length);
    if (at_end <= 0) return v;
    std::vector<double> rev;
    const std::int64_t hi = std::min(first + back, n - 1);
    for (std::int64_t i = hi; i > first; i--) rev.push_back(std::log(at_end));
    for (std::int64_t i = std::min(first, n - 1); i >= 0; i--) {
      const double p = lib.survival_at(e2 - i + other_length);
      if (p <= 0) break;
      rev.push_back(std::log(p));
    }
    v.begin = static_cast<std::size_t>(hi - static_cast<std::int64_t>(rev.size()) + 1);
    v.log_p.assign(rev.rbegin(), rev.rend());
  }
  // The most likely position: the mate's own end.
  e.position = e.plus ? static_cast<std::size_t>(std::max<std::int64_t>(0, static_cast<std::int64_t>(m.target_end) - 1)) : m.target_begin;
  return v;
}

inline Interval peak(std::size_t at, std::size_t slop, std::size_t contig_length) {
  Interval v;
  v.begin = at > slop ? at - slop : 0;
  const std::size_t e = std::min(contig_length, at + slop + 1);
  for (std::size_t x = v.begin; x < e; x++) v.log_p.push_back(-std::fabs(static_cast<double>(x) - static_cast<double>(at)) / 2.0);
  return v;
}

// The read interval of a part in the read's own coordinates.
inline void read_span(const Mapping& m, std::size_t& b, std::size_t& e) {
  const std::size_t L = m.query_length;
  b = m.reverse ? L - m.alignment.query_end : m.alignment.query_begin;
  e = m.reverse ? L - m.alignment.query_begin : m.alignment.query_end;
}

inline Breakend exit_of(const Mapping& m) { return m.reverse ? Breakend{m.contig, m.target_begin, false} : Breakend{m.contig, m.target_end - 1, true}; }
inline Breakend entry_of(const Mapping& m) { return m.reverse ? Breakend{m.contig, m.target_end - 1, true} : Breakend{m.contig, m.target_begin, false}; }

}  // namespace detail

// The evidence of discordant pairs, split reads and long gaps, and the soft clips; `bases` gives
// each mate's bases as sequenced (by read index 2i + m), for the novel bases of splits and clips.
template <class Bases>
void gather_evidence(const std::vector<Codes>& contigs, const MateMappings& mates, const std::vector<PairRecord>& pairs, const Library& lib, Bases&& bases,
                     std::vector<Evidence>& evidence, std::vector<Clip>& clips, const EvidenceOptions& o = {}) {
  // Discordant pairs.
  for (const auto& p : pairs) {
    if (!p.discordant()) continue;
    const Mapping &x = *p.mate[0], &y = *p.mate[1];
    Breakend ex, ey;
    const Interval ix = detail::mate_interval(x, unclipped_end(y) - unclipped_begin(y), lib, contigs[x.contig].size(), ex);
    const Interval iy = detail::mate_interval(y, unclipped_end(x) - unclipped_begin(x), lib, contigs[y.contig].size(), ey);
    if (ix.log_p.empty() || iy.log_p.empty()) continue;
    Evidence e;
    e.source = Source::pair;
    e.read = p.pair;
    e.weight = placement_probability(x.mapq) * placement_probability(y.mapq);
    if (ey < ex) { e.junction = {ey, ex, {}}; e.a = iy; e.b = ix; }
    else { e.junction = {ex, ey, {}}; e.a = ix; e.b = iy; }
    evidence.push_back(std::move(e));
  }
  // Split reads, gaps and clips, mate by mate.
  for (std::size_t r = 0; r < mates.size(); r++) {
    const auto& ms = mates[r];
    std::vector<const Mapping*> parts;
    for (const auto& m : ms)
      if (m.primary && m.mapq >= o.min_mapq) parts.push_back(&m);
    if (parts.empty()) continue;
    const auto& read = bases(r);  // the mate's bases as sequenced
    const std::size_t L = read.size();
    // Gaps inside an alignment.
    for (const Mapping* m : parts) {
      std::size_t i = m->alignment.query_begin, j = m->alignment.target_begin;
      for (std::size_t k = 0; k < m->alignment.cigar.size(); k++) {
        const auto& op = m->alignment.cigar[k];
        if ((op.op == 'D' || op.op == 'I') && op.length >= o.min_gap && j > 0) {
          Evidence e;
          e.source = Source::gap;
          e.read = static_cast<std::uint32_t>(r);
          e.weight = placement_probability(m->mapq);
          if (op.op == 'D') e.junction = {{m->contig, j - 1, true}, {m->contig, j + op.length, false}, {}};
          else {
            e.junction = {{m->contig, j - 1, true}, {m->contig, j, false}, {}};
            for (std::uint32_t t = 0; t < op.length; t++) {
              const std::size_t q = i + t;  // on the strand mapped
              e.junction.inserted += sequence::dna_letter(m->reverse ? sequence::complement(read[L - 1 - q]) : read[q]);
            }
          }
          e.junction = canonical(contigs, e.junction);
          e.a = detail::peak(e.junction.a.position, o.split_slop, contigs[e.junction.a.contig].size());
          e.b = detail::peak(e.junction.b.position, o.split_slop, contigs[e.junction.b.contig].size());
          evidence.push_back(std::move(e));
        }
        if (align::consumes_query(op.op)) i += op.length;
        if (align::consumes_target(op.op)) j += op.length;
      }
    }
    // Split reads: consecutive parts in the read's order.
    std::vector<std::pair<std::size_t, const Mapping*>> order;
    for (const Mapping* m : parts) {
      std::size_t b, e;
      detail::read_span(*m, b, e);
      if (e - b >= o.min_part) order.push_back({b, m});
    }
    std::sort(order.begin(), order.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
    for (std::size_t k = 0; k + 1 < order.size(); k++) {
      const Mapping &P = *order[k].second, &Q = *order[k + 1].second;
      std::size_t s1, e1, s2, e2;
      detail::read_span(P, s1, e1);
      detail::read_span(Q, s2, e2);
      Adjacency j{detail::exit_of(P), detail::entry_of(Q), {}};
      if (e1 > s2) j = Adjacency{j.a, detail::moved(j.b, static_cast<std::int64_t>(e1 - s2), false), {}};
      else
        for (std::size_t t = e1; t < s2; t++) j.inserted += sequence::dna_letter(read[t]);
      Evidence e;
      e.source = Source::split;
      e.read = static_cast<std::uint32_t>(r);
      e.weight = placement_probability(P.mapq) * placement_probability(Q.mapq);
      e.junction = canonical(contigs, j);
      e.a = detail::peak(e.junction.a.position, o.split_slop, contigs[e.junction.a.contig].size());
      e.b = detail::peak(e.junction.b.position, o.split_slop, contigs[e.junction.b.contig].size());
      evidence.push_back(std::move(e));
    }
    // Soft clips of the read's ends not explained by another part.
    if (order.empty()) continue;
    const auto clip_of = [&](const Mapping& m, bool first_in_read) {
      // The read's own first bases clipped (first_in_read) or its last.
      std::size_t b, e;
      detail::read_span(m, b, e);
      const std::size_t n = first_in_read ? b : L - e;
      if (n < o.min_clip) return;
      Clip c;
      c.read = static_cast<std::uint32_t>(r);
      c.mapq = m.mapq;
      // Read's start clipped: the junction is entered at the part's entry; else left at its exit.
      c.at = first_in_read ? detail::entry_of(m) : detail::exit_of(m);
      std::string s;
      if (first_in_read) for (std::size_t t = 0; t < b; t++) s += sequence::dna_letter(read[t]);
      else for (std::size_t t = e; t < L; t++) s += sequence::dna_letter(read[t]);
      // In the derivative read through the part forward (leaving at a + exit, entering at a − entry),
      // the read's orientation is the derivative's when the part is forward.
      c.bases = m.reverse ? sequence::reverse_complement(s) : s;
      clips.push_back(std::move(c));
    };
    clip_of(*order.front().second, true);
    clip_of(*order.back().second, false);
  }
}

}  // namespace genomics::structural
