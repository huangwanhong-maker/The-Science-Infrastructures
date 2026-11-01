// From junctions to rearrangements: breakpoints (cluster.hpp) and soft-clip piles interpreted as
// the events that make them, each with its size, its support, its genotype (genotype.hpp) and,
// where the depth measures it, its copy number (theory/structural-variants.md, §2–3).
//
// - One junction is an event by its signature (breakend.hpp): on one contig, a (+, −) adjacency
//   (p, +) ~ (q, −) is a deletion of [p + 1, q), a (−, +) adjacency (p, −) ~ (q, +) a tandem
//   duplication of [p, q + 1); novel bases between two breakends one base apart are an insertion.
// - Two junctions are one event when their four breakends pair up into two sides, each side one +
//   and one − breakend of one contig, and each side either adjacent (the − one base after the +,
//   within a tolerance: a cut of the reference) or spanning (the − at a and the + at b − 1 > a: a
//   segment [a, b) whose two ends are joined elsewhere):
//     adjacent, adjacent on one contig, (+,+) and (−,−)   an inversion of the segment between the cuts
//     adjacent, adjacent on two contigs                   a balanced reciprocal translocation, t(A;B)
//     adjacent, spanning                                  a dispersed duplication: the spanning segment
//                                                         copied into the adjacent cut, reversed when
//                                                         the junctions are (+,+) and (−,−)
//   (the patterns of the rewrites of simulate/structural.hpp; GRIDSS2 and LINX read such chains of
//   breakends more generally, Cameron et al. 2021, Genome Biology 22, 202). Inversion and
//   translocation junctions left unpaired are an inversion seen at one end and a breakend (BND).
// - Insertions from soft clips. A pile of reads clipped after base x (+) and a pile clipped before
//   base y (−) within `insertion_window` bases (y = x + 1 for a plain insertion, y ≤ x for a target
//   site duplication), neither explained by a breakpoint, is an insertion of sequence absent from the
//   reference (or a mobile element's copy the mapper did not place). Its bases: the + pile's clips
//   read after the junction and the − pile's read before it, each a column-wise majority consensus;
//   the insertion is resolved when the + consensus runs into the reference after y, or the −
//   consensus into the reference before x, or the − consensus ends with what the + consensus begins
//   with (a short insertion seen whole from both sides), or the two overlap by ≥ 12 bases;
//   otherwise its length is only bounded below, and it is imprecise. One-end-anchored pairs pointing
//   to it add support (the mates inside a long insertion are unmapped). This is a two-sided
//   consensus of the junction (DELLY's split-read refinement aligns such consensuses to the
//   reference); Manta and GRIDSS assemble every breakend region (Chen et al. 2016; Cameron et al.
//   2017), which is not written here.
// - Copy number. For a deletion or duplication spanning whole mappable bins, its copy number is
//   2 Σ n / Σ e over them (depth.hpp), and a duplication's extra copies k come from it.
// - Genotype (genotype.hpp): the alternate support is the event's evidence (for two junctions the
//   mean over them), the reference support the mean over its junctions of the reads and proper pairs
//   crossing the reference there; the expected fractions 1/2 and 1 − ε (k/(k + 2) and k/(k + 1) for
//   a duplication of k extra copies). An event of QUAL below `min_qual` is filtered (LowQual).
//
// Complexity: O(B² ) for B breakpoints of one variety pair (they are few), O(C log C) for C clips.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/structural/breakend.hpp>
#include <genomics/structural/cluster.hpp>
#include <genomics/structural/depth.hpp>
#include <genomics/structural/evidence.hpp>
#include <genomics/structural/genotype.hpp>
#include <genomics/structural/pairs.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace genomics::structural {

enum class SvKind : std::uint8_t { deletion, duplication, inversion, insertion, translocation, dispersed_duplication, breakend, copy_loss, copy_gain };

// The names VCF's SVTYPE gives (VCFv4.2 §3; DUP:TANDEM and DUP:DISPERSED are its subtypes, CNV the
// copy-number one).
inline const char* sv_kind_name(SvKind k) {
  switch (k) {
    case SvKind::deletion: return "DEL";
    case SvKind::duplication: return "DUP:TANDEM";
    case SvKind::inversion: return "INV";
    case SvKind::insertion: return "INS";
    case SvKind::translocation: return "TRA";
    case SvKind::dispersed_duplication: return "DUP:DISPERSED";
    case SvKind::breakend: return "BND";
    case SvKind::copy_loss: return "CNV:LOSS";
    case SvKind::copy_gain: return "CNV:GAIN";
  }
  return "?";
}

struct SvCall {
  SvKind kind = SvKind::deletion;
  std::uint32_t contig = 0;
  std::size_t begin = 0, end = 0;           // the reference bases rewritten [begin, end); an insertion's site begin = end
  std::uint32_t partner = 0;                // a translocation's other contig, a dispersed copy's destination contig
  std::size_t partner_position = 0;         // its breakpoint (the first base after the cut)
  bool inverted = false;                    // a dispersed copy (or translocation) joined reversed
  std::string sequence;                     // an insertion's bases, as far as resolved
  bool resolved = false;                    // the insertion's sequence is whole
  std::size_t length = 0;                   // SVLEN's magnitude (for an unresolved insertion, a lower bound)
  bool precise = false;
  std::array<std::size_t, 2> ci_begin{0, 0}, ci_end{0, 0};  // confidence intervals of begin and end
  std::vector<Adjacency> junctions;         // the junctions it makes (one or two)
  std::vector<std::uint32_t> breakpoints;   // their indices
  std::size_t pairs = 0, splits = 0, clips = 0, anchored = 0;  // support
  double copies = -1;                       // the copy number the depth measures over it (−1: none)
  SvGenotype genotype;
  double qual = 0;
  std::vector<std::string> filters;
  bool passed() const { return filters.empty(); }
};

struct CallOptions {
  std::size_t min_size = 50;           // GIAB's and Truvari's smallest structural variant
  std::size_t precise_tolerance = 20;  // how far two precise breakends of one cut may lie
  double imprecise_sds = 3;            // ... and imprecise ones: this many library sds
  std::size_t min_pile = 2;            // clipped reads a side of an insertion needs
  std::size_t min_insertion_support = 5;
  std::size_t insertion_window = 20;
  std::size_t explained = 30;          // clips this near a breakpoint's breakend are its own
  double min_qual = 20;
  std::size_t depth_check = 5000;                           // events this long are held to the depth
  double deletion_copies = 1.6, duplication_copies = 2.4;   // the copy numbers they must stay below, above
  GenotypeOptions genotype;
};

namespace detail {

// The column-wise majority of left-aligned (or right-aligned) strings, kept while ≥ 2 agree on ≥ 60%.
inline std::string consensus(std::vector<std::string> s, bool right_aligned) {
  if (right_aligned)
    for (auto& x : s) std::reverse(x.begin(), x.end());
  std::string out;
  for (std::size_t col = 0;; col++) {
    std::array<std::size_t, 5> n{};
    std::size_t total = 0;
    for (const auto& x : s)
      if (col < x.size()) { n[sequence::dna_code(x[col])]++; total++; }
    if (total < 2) break;
    const auto best = static_cast<std::size_t>(std::max_element(n.begin(), n.begin() + 4) - n.begin());
    if (static_cast<double>(n[best]) < 0.6 * static_cast<double>(total)) break;
    out += sequence::DNA_LETTERS[best];
  }
  if (right_aligned) std::reverse(out.begin(), out.end());
  return out;
}

inline std::string letters(const Codes& x, std::size_t b, std::size_t e) {
  std::string s;
  for (std::size_t i = b; i < std::min(e, x.size()); i++) s += sequence::dna_letter(x[i]);
  return s;
}

// The insertion between bases x and y from the + consensus p (bases after x) and the − consensus m (before y).
inline void resolve_insertion(const Codes& ref, std::size_t x, std::size_t y, const std::string& p, const std::string& m, SvCall& call) {
  const std::size_t k = 12;
  const std::string right = letters(ref, y, y + k), left = y >= 1 && x + 1 >= k ? letters(ref, x + 1 - k, x + 1) : std::string();
  // The + consensus running into the reference after the insertion (at any site the two piles
  // allow), or the − consensus starting in the reference before it.
  (void)right;
  (void)left;
  const std::size_t w = 20;
  std::size_t best_t = std::string::npos;
  for (std::size_t q = y > w ? y - w : 0; q <= y + w; q++) {
    const std::string r = letters(ref, q, q + k);
    if (r.size() != k) continue;
    const std::size_t t = p.find(r);
    if (t != std::string::npos && (best_t == std::string::npos || t < best_t)) best_t = t;
  }
  if (best_t != std::string::npos) { call.sequence = p.substr(0, best_t); call.resolved = true; return; }
  std::size_t best_u = std::string::npos;
  for (std::size_t q = x + 1 > w + k ? x + 1 - w - k : 0; q <= x + 1 + w && q + k <= ref.size(); q++) {
    const std::string l = letters(ref, q, q + k);
    const std::size_t t = m.rfind(l);
    if (t != std::string::npos && (best_u == std::string::npos || t + k > best_u)) best_u = t + k;
  }
  if (best_u != std::string::npos) { call.sequence = m.substr(best_u); call.resolved = true; return; }
  // A short insertion seen whole from both sides: the − consensus ends with it and the + consensus
  // begins with it (m = …flank · u, p = u · flank…).
  for (std::size_t o = std::min(p.size(), m.size()); o >= k; o--)
    if (m.compare(m.size() - o, o, p, 0, o) == 0) { call.sequence = p.substr(0, o); call.resolved = true; return; }
  // A long one seen in two parts that overlap: p = u[0, |p|), m = u[|u| − |m|, |u|).
  for (std::size_t o = std::min(p.size(), m.size()); o >= k; o--)
    if (p.compare(p.size() - o, o, m, 0, o) == 0) { call.sequence = p + m.substr(o); call.resolved = true; return; }
  call.sequence = p + m;
  call.resolved = false;
}

struct Side {
  std::uint32_t contig;
  Breakend plus, minus;
};

// The two sides of a pair of junctions: their first breakends together, their second together.
inline bool sides_of(const Adjacency& u, const Adjacency& v, Side& s1, Side& s2) {
  const auto side = [](const Breakend& x, const Breakend& y, Side& s) {
    if (x.contig != y.contig || x.plus == y.plus) return false;
    s.contig = x.contig;
    s.plus = x.plus ? x : y;
    s.minus = x.plus ? y : x;
    return true;
  };
  return side(u.a, v.a, s1) && side(u.b, v.b, s2);
}

inline bool adjacent(const Side& s, double tol) { return std::fabs(static_cast<double>(s.minus.position) - static_cast<double>(s.plus.position) - 1) <= tol; }
inline bool spanning(const Side& s, std::size_t min_size) { return s.plus.position + 1 >= s.minus.position + min_size; }

}  // namespace detail

// The copy number the depth measures over [b, e) of a contig: 2 Σ n / Σ e over its whole mappable bins, or −1.
inline double copies_over(const Depth& d, std::uint32_t contig, std::size_t b, std::size_t e) {
  if (contig + 1 >= d.first.size()) return -1;
  double n = 0, x = 0;
  for (std::size_t i = d.first[contig]; i < d.first[contig + 1]; i++) {
    const auto& bin = d.bins[i];
    if (bin.begin >= b && bin.end <= e && bin.mappable) { n += bin.count; x += bin.exposure; }
  }
  return x > 0 ? 2 * n / x : -1;
}

// The events of the breakpoints and clips. `depth` may be empty (no bins): no copy numbers then.
inline std::vector<SvCall> interpret(const std::vector<Codes>& contigs, const std::vector<Evidence>& evidence, const std::vector<Breakpoint>& bps, const std::vector<Clip>& clips,
                                     const std::vector<PairRecord>& pairs, const Library& lib, const ReadIndex& index, const Depth& depth, const CallOptions& o = {}) {
  std::vector<SvCall> out;
  std::vector<bool> used(bps.size(), false);
  const auto tol = [&](const Breakpoint& u, const Breakpoint& v) {
    return u.precise && v.precise ? static_cast<double>(o.precise_tolerance) : o.imprecise_sds * lib.sd() + static_cast<double>(o.precise_tolerance);
  };
  const auto add_support = [&](SvCall& c, std::size_t i) {
    c.breakpoints.push_back(static_cast<std::uint32_t>(i));
    c.junctions.push_back(bps[i].junction);
    c.pairs += bps[i].pairs;
    c.splits += bps[i].splits + bps[i].gaps;
    c.precise = c.precise || bps[i].precise;
  };
  // ---- pairs of junctions -----------------------------------------------------------------------
  for (std::size_t i = 0; i < bps.size(); i++) {
    if (used[i]) continue;
    for (std::size_t j = i + 1; j < bps.size() && !used[i]; j++) {
      if (used[j]) continue;
      const Adjacency &u = bps[i].junction, &v = bps[j].junction;
      detail::Side s1, s2;
      if (!detail::sides_of(u, v, s1, s2)) continue;
      const double t = tol(bps[i], bps[j]);
      const bool a1 = detail::adjacent(s1, t), a2 = detail::adjacent(s2, t);
      const bool sp1 = detail::spanning(s1, o.min_size), sp2 = detail::spanning(s2, o.min_size);
      const bool same_orientation = u.a.plus == u.b.plus;  // (+,+) and (−,−) junctions
      SvCall c;
      if (a1 && a2) {
        if (s1.contig == s2.contig && same_orientation) {
          c.kind = SvKind::inversion;
          c.contig = s1.contig;
          c.begin = s1.minus.position;
          c.end = s2.minus.position;
        } else if (s1.contig != s2.contig) {
          c.kind = SvKind::translocation;
          c.contig = s1.contig;
          c.begin = s1.minus.position;
          c.partner = s2.contig;
          c.partner_position = s2.minus.position;
          c.inverted = same_orientation;
        } else {
          continue;
        }
      } else if ((a1 && sp2) || (a2 && sp1)) {
        const detail::Side& cut = a1 ? s1 : s2;
        const detail::Side& seg = a1 ? s2 : s1;
        // The copy must not contain its destination.
        if (cut.contig == seg.contig && cut.minus.position >= seg.minus.position && cut.minus.position <= seg.plus.position) continue;
        c.kind = SvKind::dispersed_duplication;
        c.contig = seg.contig;
        c.begin = seg.minus.position;
        c.end = seg.plus.position + 1;
        c.partner = cut.contig;
        c.partner_position = cut.minus.position;
        c.inverted = same_orientation;
      } else {
        continue;
      }
      add_support(c, i);
      add_support(c, j);
      c.length = c.kind == SvKind::translocation ? 0 : c.end - c.begin;
      used[i] = used[j] = true;
      out.push_back(std::move(c));
    }
  }
  // ---- single junctions -----------------------------------------------------------------------
  for (std::size_t i = 0; i < bps.size(); i++) {
    if (used[i]) continue;
    const Adjacency& j = bps[i].junction;
    SvCall c;
    switch (bps[i].signature()) {
      case Signature::deletion: {
        const std::size_t gap = j.b.position > j.a.position ? j.b.position - j.a.position - 1 : 0;
        if (gap < o.min_size && j.inserted.size() >= o.min_size) {
          c.kind = SvKind::insertion;
          c.begin = c.end = j.a.position + 1;
          c.sequence = j.inserted;
          c.resolved = true;
          c.length = j.inserted.size();
        } else if (gap >= o.min_size) {
          c.kind = SvKind::deletion;
          c.begin = j.a.position + 1;
          c.end = j.b.position;
          c.length = gap;
          c.sequence = j.inserted;
        } else {
          continue;
        }
        c.contig = j.a.contig;
        break;
      }
      case Signature::duplication:
        if (j.b.position + 1 < j.a.position + o.min_size) continue;
        c.kind = SvKind::duplication;
        c.contig = j.a.contig;
        c.begin = j.a.position;
        c.end = j.b.position + 1;
        c.length = c.end - c.begin;
        break;
      case Signature::inversion_plus:
      case Signature::inversion_minus: {
        const bool plus = j.a.plus;
        c.kind = SvKind::inversion;
        c.contig = j.a.contig;
        c.begin = plus ? j.a.position + 1 : j.a.position;
        c.end = plus ? j.b.position + 1 : j.b.position;
        if (c.end < c.begin + o.min_size) continue;
        c.length = c.end - c.begin;
        break;
      }
      case Signature::translocation:
        c.kind = SvKind::breakend;
        c.contig = j.a.contig;
        c.begin = c.end = j.a.position;
        c.partner = j.b.contig;
        c.partner_position = j.b.position;
        c.inverted = j.a.plus == j.b.plus;
        break;
    }
    add_support(c, i);
    used[i] = true;
    out.push_back(std::move(c));
  }
  // ---- insertions from clip piles ----------------------------------------------------------------
  {
    // Clips by exact breakend.
    std::map<std::tuple<std::uint32_t, std::size_t, bool>, std::vector<std::size_t>> piles;
    for (std::size_t i = 0; i < clips.size(); i++) piles[{clips[i].at.contig, clips[i].at.position, clips[i].at.plus}].push_back(i);
    const auto explained = [&](std::uint32_t contig, std::size_t p) {
      for (const auto& b : bps)
        for (const Breakend* e : {&b.junction.a, &b.junction.b})
          if (e->contig == contig && (e->position > p ? e->position - p : p - e->position) <= o.explained) return true;
      return false;
    };
    std::map<std::tuple<std::uint32_t, std::size_t, bool>, bool> done;
    for (const auto& [key, plus_pile] : piles) {
      const auto [contig, x, plus] = key;
      if (!plus || plus_pile.size() < o.min_pile || done[key] || explained(contig, x)) continue;
      // The best − pile near x + 1.
      const std::vector<std::size_t>* minus_pile = nullptr;
      std::size_t y = 0;
      const std::size_t from = x + 1 > o.insertion_window ? x + 1 - o.insertion_window : 0;
      for (auto it = piles.lower_bound({contig, from, false}); it != piles.end() && std::get<0>(it->first) == contig && std::get<1>(it->first) <= x + 1 + o.insertion_window; ++it) {
        if (std::get<2>(it->first) || it->second.size() < o.min_pile || done[it->first]) continue;
        if (!minus_pile || it->second.size() > minus_pile->size()) { minus_pile = &it->second; y = std::get<1>(it->first); }
      }
      if (!minus_pile) continue;
      done[key] = true;
      done[{contig, y, false}] = true;
      SvCall c;
      c.kind = SvKind::insertion;
      c.contig = contig;
      c.begin = c.end = x + 1;
      c.precise = true;
      std::vector<std::string> ps, ms;
      for (auto i : plus_pile) ps.push_back(clips[i].bases);
      for (auto i : *minus_pile) ms.push_back(clips[i].bases);
      c.clips = ps.size() + ms.size();
      detail::resolve_insertion(contigs[contig], x, y, detail::consensus(ps, false), detail::consensus(ms, true), c);
      c.length = c.sequence.size();
      // One-end-anchored pairs pointing to the site.
      for (const auto& p : pairs) {
        if (p.type != PairClass::one_end) continue;
        const Mapping& m = *(p.mate[0] ? p.mate[0] : p.mate[1]);
        if (m.contig != contig) continue;
        const auto near = static_cast<std::int64_t>(lib.max_normal);
        if (!m.reverse && unclipped_begin(m) <= static_cast<std::int64_t>(x) && static_cast<std::int64_t>(x) - unclipped_begin(m) <= near) c.anchored++;
        if (m.reverse && unclipped_end(m) > static_cast<std::int64_t>(y) && unclipped_end(m) - static_cast<std::int64_t>(y) <= near) c.anchored++;
      }
      if (c.clips + c.anchored < o.min_insertion_support) continue;
      c.junctions.push_back({{contig, x, true}, {contig, y, false}, c.sequence});
      out.push_back(std::move(c));
    }
  }
  // ---- copy numbers, genotypes, qualities ------------------------------------------------------------
  for (auto& c : out) {
    if ((c.kind == SvKind::deletion || c.kind == SvKind::duplication) && !depth.bins.empty()) c.copies = copies_over(depth, c.contig, c.begin, c.end);
    if (c.kind == SvKind::dispersed_duplication && !depth.bins.empty()) c.copies = copies_over(depth, c.contig, c.begin, c.end);
    // Support.
    double alt = 0, ref = 0;
    if (c.kind == SvKind::insertion && c.breakpoints.empty()) {
      for (const auto& cl : clips)
        if (cl.at.contig == c.contig && ((cl.at.plus && cl.at.position + 1 == c.begin) || (!cl.at.plus && cl.at.position + o.insertion_window >= c.begin && cl.at.position <= c.begin + o.insertion_window)))
          alt += placement_probability(cl.mapq);
      alt += static_cast<double>(c.anchored);
      ref = reference_support(index, c.junctions.front(), o.genotype, -static_cast<double>(c.length));
    } else {
      const double delta = c.kind == SvKind::deletion ? static_cast<double>(c.length) : c.kind == SvKind::insertion ? -static_cast<double>(c.length) : 0;
      for (auto b : c.breakpoints) {
        for (auto e : bps[b].evidence) alt += evidence[e].weight;
        ref += reference_support(index, bps[b].junction, o.genotype, delta);
      }
      alt /= static_cast<double>(c.breakpoints.size());
      ref /= static_cast<double>(c.breakpoints.size());
    }
    double p_het = 0.5, p_hom = 1;
    if (c.kind == SvKind::duplication) {
      const double k = c.copies > 0 ? std::max(1.0, std::round(c.copies - 2)) : 1.0;
      const double kh = c.copies > 0 ? std::max(1.0, std::round((c.copies - 2) / 2)) : 1.0;
      p_het = k / (k + 2);
      p_hom = kh / (kh + 1);
    }
    c.genotype = genotype_of(ref, alt, p_het, p_hom, o.genotype.error);
    // A deletion spanning bins whose depth says how many copies remain.
    if (c.kind == SvKind::deletion && c.copies >= 0) {
      if (c.copies < 0.5 && c.genotype.genotype == 1) c.genotype.genotype = 2;
      if (c.copies > 0.5 && c.copies < 1.5 && c.genotype.genotype == 2) c.genotype.genotype = 1;
    }
    c.qual = c.genotype.qual;
    if (c.qual < o.min_qual || c.genotype.genotype == 0) c.filters.push_back("LowQual");
    // The depth against the event, as duphold holds calls to the fold change of depth (Pedersen and
    // Quinlan 2019, GigaScience 8, giz040): a deletion of whole bins that keeps two copies, or a
    // duplication that gains none, is a junction of another event read alone.
    if (c.length >= o.depth_check && c.copies >= 0 &&
        ((c.kind == SvKind::deletion && c.copies > o.deletion_copies) || (c.kind == SvKind::duplication && c.copies < o.duplication_copies)))
      c.filters.push_back("DepthConflict");
    // Confidence intervals of the event's ends from its first junction.
    if (!c.breakpoints.empty()) {
      const Breakpoint& b = bps[c.breakpoints.front()];
      c.ci_begin = {b.ci_a[0], b.ci_a[1]};
      c.ci_end = {b.ci_b[0], b.ci_b[1]};
    } else {
      c.ci_begin = {c.begin, c.begin};
      c.ci_end = {c.end, c.end};
    }
  }
  std::sort(out.begin(), out.end(), [](const SvCall& x, const SvCall& y) {
    if (x.contig != y.contig) return x.contig < y.contig;
    if (x.begin != y.begin) return x.begin < y.begin;
    return x.kind < y.kind;
  });
  return out;
}

// Copy-number segments other than the baseline that no deletion or duplication call explains (half
// of the segment inside the call), and not on an aneuploid contig: copy-number calls of their own.
template <class Segments, class Aneuploidies>
std::vector<SvCall> copy_number_calls(const Segments& segments, const Aneuploidies& aneuploid, const std::vector<SvCall>& calls, std::size_t min_length = 10000, double overlap = 0.5) {
  std::vector<SvCall> out;
  for (const auto& s : segments) {
    if (s.copies == 2 || s.end - s.begin < min_length) continue;
    bool whole = false;
    for (const auto& a : aneuploid) whole = whole || a.contig == s.contig;
    if (whole) continue;
    bool explained = false;
    for (const auto& c : calls) {
      if (c.contig != s.contig || (c.kind != SvKind::deletion && c.kind != SvKind::duplication)) continue;
      const std::size_t lo = std::max(c.begin, s.begin), hi = std::min(c.end, s.end);
      if (hi <= lo) continue;
      // The segment mostly inside the call (a large event's depth may be cut by a repeat inside it).
      if (static_cast<double>(hi - lo) >= overlap * static_cast<double>(s.end - s.begin)) explained = true;
    }
    if (explained) continue;
    SvCall c;
    c.kind = s.copies < 2 ? SvKind::copy_loss : SvKind::copy_gain;
    c.contig = s.contig;
    c.begin = s.begin;
    c.end = s.end;
    c.length = s.end - s.begin;
    c.copies = s.estimate();
    c.qual = s.quality;
    c.genotype.genotype = s.copies == 0 ? 2 : 1;
    c.ci_begin = {s.begin, s.begin};
    c.ci_end = {s.end, s.end};
    out.push_back(c);
  }
  return out;
}

}  // namespace genomics::structural
