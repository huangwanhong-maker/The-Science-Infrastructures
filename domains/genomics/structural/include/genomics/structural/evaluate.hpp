// Calls against a truth set: Truvari's matching of structural variants (English et al. 2022,
// "Truvari: refined structural variant comparison preserves allelic diversity", Genome Biology 23,
// 271; its bench, defaults of version 4), precision and recall by type and size as GIAB's benchmark
// reports them (Zook et al. 2020, Nature Biotechnology 38, 1347–1355), the accuracy of breakpoints,
// and the base-level concordance of copy number and of loss of heterozygosity.
//
// - Matching (truvari bench). A call and a truth variant of one type are compared when they lie on
//   one contig with |Δstart| ≤ refdist and |Δend| ≤ refdist (500), and match when their size
//   similarity min/max ≥ pctsize (0.7) and their reciprocal overlap, the shared bases over the
//   longer span, ≥ pctovl (0); sequence similarity (pctseq) needs both alleles' bases and is not
//   used (Truvari's --pctseq 0). Matches are chosen one to one ("pick single") in decreasing
//   TruScore, the mean of the size similarity and the overlap. Variants shorter than sizemin (50)
//   are not compared. Precision = TP_call/(TP_call + FP), recall = TP_base/(TP_base + FN), F₁
//   their harmonic mean, genotype concordance over the matches.
// - Kinds Truvari does not compare, matched here by the same distance: an insertion by its site
//   (its size only when the call's sequence is resolved; an unresolved call is a lower bound that
//   must not exceed the truth); a translocation by both breakpoints; a dispersed duplication by
//   its copied interval and its destination. Copy-number calls (from depth alone, whose ends are
//   known to a bin) match a truth deletion or duplication by reciprocal overlap ≥ 0.5, the
//   convention of copy-number benchmarks (*Change*, so that a bin's imprecision is not a miss).
// - Breakpoints: for every match of a precise call, the distance of each of its junctions from the
//   nearest junction of the truth, both written canonically (breakend.hpp), so that microhomology
//   costs nothing.
// - Copy number: over the mappable bins, the fraction whose called integer copy number equals the
//   truth's at the bin's middle, and the confusion of the two.
// - LOH: over bases, the called copy-neutral LOH against the truth's (recall, precision).
//
// Complexity: O(C T) for C calls and T truth variants, which are few.

#pragma once

#include <genomics/structural/allelic.hpp>
#include <genomics/structural/breakend.hpp>
#include <genomics/structural/call.hpp>
#include <genomics/structural/copy_number.hpp>
#include <genomics/structural/depth.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace genomics::structural {

// A truth variant in the calls' language.
struct TruthEvent {
  SvKind kind = SvKind::deletion;
  std::uint32_t contig = 0;
  std::size_t begin = 0, end = 0;
  std::uint32_t partner = 0;
  std::size_t partner_position = 0;
  bool inverted = false;
  std::size_t length = 0;
  std::string sequence;
  int genotype = 1;
  int copies = 1;     // extra copies of a tandem duplication
  std::size_t index = 0;  // in the truth it came from
};

// The structural variants of a simulated sample (simulate::StructuralVariant, or any type with its
// members) as truth events; aneuploidy and LOH are not junctions and are compared apart.
template <class V>
std::vector<TruthEvent> truth_events(const std::vector<V>& variants) {
  std::vector<TruthEvent> out;
  for (std::size_t i = 0; i < variants.size(); i++) {
    const auto& v = variants[i];
    TruthEvent t;
    const std::string type = sv_type_name(v.type);
    if (type == "DEL") t.kind = SvKind::deletion;
    else if (type == "DUP:TANDEM") t.kind = SvKind::duplication;
    else if (type == "DUP:DISPERSED") t.kind = SvKind::dispersed_duplication;
    else if (type == "INV") t.kind = SvKind::inversion;
    else if (type == "INS") t.kind = SvKind::insertion;
    else if (type == "TRA") t.kind = SvKind::translocation;
    else continue;
    t.contig = v.contig;
    t.begin = v.begin;
    t.end = v.end;
    t.partner = v.partner;
    t.partner_position = v.partner_position;
    t.inverted = v.inverted;
    t.length = v.length();
    t.sequence = v.sequence;
    t.genotype = v.on[0] && v.on[1] ? 2 : 1;
    t.copies = v.copies;
    t.index = i;
    out.push_back(t);
  }
  return out;
}

// The junctions a truth event makes, in the calls' convention (one copy of each).
inline std::vector<Adjacency> junctions_of(const TruthEvent& t) {
  const std::uint32_t c = t.contig;
  switch (t.kind) {
    case SvKind::deletion:
    case SvKind::copy_loss: return {{{c, t.begin - 1, true}, {c, t.end, false}, {}}};
    case SvKind::duplication:
    case SvKind::copy_gain: return {{{c, t.begin, false}, {c, t.end - 1, true}, {}}};
    case SvKind::inversion: return {{{c, t.begin - 1, true}, {c, t.end - 1, true}, {}}, {{c, t.begin, false}, {c, t.end, false}, {}}};
    case SvKind::insertion: return {{{c, t.begin - 1, true}, {c, t.begin, false}, t.sequence}};
    case SvKind::translocation: {
      const std::uint32_t d = t.partner;
      const std::size_t p = t.begin, q = t.partner_position;
      return {{{c, p - 1, true}, {d, q, false}, {}}, {{d, q - 1, true}, {c, p, false}, {}}};
    }
    case SvKind::dispersed_duplication: {
      const std::uint32_t d = t.partner;
      const std::size_t p = t.partner_position;
      if (!t.inverted) return {{{d, p - 1, true}, {c, t.begin, false}, {}}, {{c, t.end - 1, true}, {d, p, false}, {}}};
      return {{{d, p - 1, true}, {c, t.end - 1, true}, {}}, {{c, t.begin, false}, {d, p, false}, {}}};
    }
    case SvKind::breakend: return {};
  }
  return {};
}

struct MatchOptions {
  std::size_t refdist = 500;
  double pctsize = 0.7;
  double pctovl = 0.0;
  std::size_t sizemin = 50;
  bool pass_only = true;
  double cnv_overlap = 0.5;
};

struct SvMatch {
  std::size_t call = 0, truth = 0;
  double size_similarity = 1, overlap = 0, score = 0;
  bool genotype = false;
  std::size_t breakpoint_error = 0;  // the largest distance of a junction from the truth's (precise calls)
  bool precise = false;
};

struct Tally {
  std::size_t tp_base = 0, tp_call = 0, fp = 0, fn = 0, genotype = 0;
  double precision() const { return tp_call + fp ? static_cast<double>(tp_call) / static_cast<double>(tp_call + fp) : 1; }
  double recall() const { return tp_base + fn ? static_cast<double>(tp_base) / static_cast<double>(tp_base + fn) : 1; }
  double f1() const { const double p = precision(), r = recall(); return p + r > 0 ? 2 * p * r / (p + r) : 0; }
  double genotype_concordance() const { return tp_call ? static_cast<double>(genotype) / static_cast<double>(tp_call) : 0; }
};

// GIAB's size strata, extended to copy-number sizes.
inline constexpr std::array<std::size_t, 6> SIZE_EDGES{50, 100, 1000, 10000, 100000, SIZE_MAX};
inline const char* size_bin_name(std::size_t k) {
  static const char* names[] = {"50–99", "100–999", "1k–9.9k", "10k–99k", "≥100k"};
  return k < 5 ? names[k] : "?";
}
inline std::size_t size_bin(std::size_t length) {
  for (std::size_t k = 0; k + 1 < SIZE_EDGES.size(); k++)
    if (length < SIZE_EDGES[k + 1]) return k;
  return 4;
}

struct SvComparison {
  std::vector<SvMatch> matches;
  std::vector<std::size_t> false_positives, false_negatives;  // call and truth indices
  Tally overall;
  std::map<std::string, Tally> by_kind;    // by the truth's kind (calls' kind for FPs)
  std::array<Tally, 5> by_size{};
  std::vector<std::size_t> breakpoint_errors;  // of precise matches
};

namespace detail {

inline const char* kind_class(SvKind k) {
  switch (k) {
    case SvKind::copy_loss: return "DEL";
    case SvKind::copy_gain: return "DUP:TANDEM";
    default: return sv_kind_name(k);
  }
}

inline std::size_t distance(std::size_t a, std::size_t b) { return a > b ? a - b : b - a; }

inline std::size_t junction_error(const std::vector<Codes>& contigs, const std::vector<Adjacency>& calls, const std::vector<Adjacency>& truth) {
  std::size_t worst = 0;
  for (const auto& c : calls) {
    const Adjacency x = canonical(contigs, c);
    std::size_t best = SIZE_MAX;
    for (const auto& t : truth) {
      const Adjacency y = canonical(contigs, t);
      if (x.a.contig != y.a.contig || x.b.contig != y.b.contig || x.a.plus != y.a.plus || x.b.plus != y.b.plus) continue;
      best = std::min(best, std::max(distance(x.a.position, y.a.position), distance(x.b.position, y.b.position)));
    }
    worst = std::max(worst, best);
  }
  return worst;
}

// The distance of an insertion's site from the sites of the truth's equivalent placements: u
// inserted before base p makes the haplotype of its rotation inserted k bases later while
// ref[p + j] = u[j mod |u|], or earlier while ref[p − 1 − j] = u[|u| − 1 − (j mod |u|)].
inline std::size_t insertion_error(const Codes& ref, const TruthEvent& t, std::size_t site) {
  const std::string& u = t.sequence;
  if (u.empty()) return distance(site, t.begin);
  std::size_t lo = t.begin, hi = t.begin;
  while (hi < ref.size() && sequence::dna_letter(ref[hi]) == u[(hi - t.begin) % u.size()]) hi++;
  while (lo > 0 && sequence::dna_letter(ref[lo - 1]) == u[u.size() - 1 - ((t.begin - lo) % u.size())]) lo--;
  return site < lo ? lo - site : site > hi ? site - hi : 0;
}

}  // namespace detail

// Whether a call and a truth event may match, and how well (score < 0: not).
inline double match_score(const SvCall& c, const TruthEvent& t, const MatchOptions& o, double& size_similarity, double& overlap) {
  size_similarity = 1;
  overlap = 0;
  const bool cnv = c.kind == SvKind::copy_loss || c.kind == SvKind::copy_gain;
  if (std::string(detail::kind_class(c.kind)) != detail::kind_class(t.kind)) return -1;
  const auto span_overlap = [&](std::size_t b1, std::size_t e1, std::size_t b2, std::size_t e2) {
    const std::size_t lo = std::max(b1, b2), hi = std::min(e1, e2);
    const double shared = hi > lo ? static_cast<double>(hi - lo) : 0;
    return shared / static_cast<double>(std::max<std::size_t>(1, std::max(e1 - b1, e2 - b2)));
  };
  switch (t.kind) {
    case SvKind::deletion:
    case SvKind::duplication:
    case SvKind::inversion: {
      if (c.contig != t.contig) return -1;
      overlap = span_overlap(c.begin, c.end, t.begin, t.end);
      const double a = static_cast<double>(c.end - c.begin), b = static_cast<double>(t.end - t.begin);
      size_similarity = std::min(a, b) / std::max(a, b);
      if (cnv) {
        // Reciprocal overlap: the shared bases over each of the two.
        const std::size_t lo = std::max(c.begin, t.begin), hi = std::min(c.end, t.end);
        const double shared = hi > lo ? static_cast<double>(hi - lo) : 0;
        if (shared < o.cnv_overlap * a || shared < o.cnv_overlap * b) return -1;
      } else if (detail::distance(c.begin, t.begin) > o.refdist || detail::distance(c.end, t.end) > o.refdist || size_similarity < o.pctsize || overlap < o.pctovl) {
        return -1;
      }
      return (size_similarity + overlap) / 2;
    }
    case SvKind::insertion: {
      if (c.contig != t.contig || detail::distance(c.begin, t.begin) > o.refdist) return -1;
      const double a = static_cast<double>(std::max<std::size_t>(1, c.length)), b = static_cast<double>(std::max<std::size_t>(1, t.length));
      if (c.resolved) {
        size_similarity = std::min(a, b) / std::max(a, b);
        if (size_similarity < o.pctsize) return -1;
      } else if (a > b / o.pctsize) {
        return -1;
      }
      return (size_similarity + 1 - static_cast<double>(detail::distance(c.begin, t.begin)) / static_cast<double>(o.refdist + 1)) / 2;
    }
    case SvKind::translocation: {
      const bool same = c.contig == t.contig && c.partner == t.partner && detail::distance(c.begin, t.begin) <= o.refdist && detail::distance(c.partner_position, t.partner_position) <= o.refdist;
      const bool swapped = c.contig == t.partner && c.partner == t.contig && detail::distance(c.begin, t.partner_position) <= o.refdist && detail::distance(c.partner_position, t.begin) <= o.refdist;
      if (!same && !swapped) return -1;
      const std::size_t d = same ? std::max(detail::distance(c.begin, t.begin), detail::distance(c.partner_position, t.partner_position))
                                 : std::max(detail::distance(c.begin, t.partner_position), detail::distance(c.partner_position, t.begin));
      return 1 - static_cast<double>(d) / static_cast<double>(o.refdist + 1);
    }
    case SvKind::dispersed_duplication: {
      if (c.contig != t.contig || c.partner != t.partner || detail::distance(c.begin, t.begin) > o.refdist || detail::distance(c.end, t.end) > o.refdist ||
          detail::distance(c.partner_position, t.partner_position) > o.refdist)
        return -1;
      const double a = static_cast<double>(c.end - c.begin), b = static_cast<double>(t.end - t.begin);
      size_similarity = std::min(a, b) / std::max(a, b);
      overlap = span_overlap(c.begin, c.end, t.begin, t.end);
      if (size_similarity < o.pctsize) return -1;
      return (size_similarity + overlap) / 2;
    }
    default: return -1;
  }
}

inline SvComparison compare_sv(const std::vector<Codes>& contigs, const std::vector<SvCall>& calls, const std::vector<TruthEvent>& truth, const MatchOptions& o = {}) {
  SvComparison r;
  const auto counted_call = [&](const SvCall& c) {
    if (o.pass_only && !c.passed()) return false;
    return c.kind == SvKind::translocation || c.kind == SvKind::breakend || c.length >= o.sizemin;
  };
  const auto counted_truth = [&](const TruthEvent& t) { return t.kind == SvKind::translocation || t.length >= o.sizemin; };
  struct Candidate { double score, size, overlap; std::size_t c, t; };
  std::vector<Candidate> cand;
  for (std::size_t i = 0; i < calls.size(); i++) {
    if (!counted_call(calls[i])) continue;
    for (std::size_t j = 0; j < truth.size(); j++) {
      if (!counted_truth(truth[j])) continue;
      double size = 1, ovl = 0;
      const double s = match_score(calls[i], truth[j], o, size, ovl);
      if (s >= 0) cand.push_back({s, size, ovl, i, j});
    }
  }
  std::sort(cand.begin(), cand.end(), [](const Candidate& x, const Candidate& y) {
    if (x.score != y.score) return x.score > y.score;
    return x.c != y.c ? x.c < y.c : x.t < y.t;
  });
  std::vector<bool> call_used(calls.size(), false), truth_used(truth.size(), false);
  for (const auto& k : cand) {
    if (call_used[k.c] || truth_used[k.t]) continue;
    call_used[k.c] = truth_used[k.t] = true;
    SvMatch m;
    m.call = k.c;
    m.truth = k.t;
    m.score = k.score;
    m.size_similarity = k.size;
    m.overlap = k.overlap;
    m.genotype = calls[k.c].genotype.genotype == truth[k.t].genotype;
    m.precise = calls[k.c].precise && !calls[k.c].junctions.empty();
    if (m.precise && truth[k.t].kind == SvKind::insertion) {
      m.breakpoint_error = detail::insertion_error(contigs[truth[k.t].contig], truth[k.t], calls[k.c].begin);
      r.breakpoint_errors.push_back(m.breakpoint_error);
    } else if (m.precise) {
      m.breakpoint_error = detail::junction_error(contigs, calls[k.c].junctions, junctions_of(truth[k.t]));
      r.breakpoint_errors.push_back(m.breakpoint_error);
    }
    r.matches.push_back(m);
  }
  for (const auto& m : r.matches) {
    const std::string kind = detail::kind_class(truth[m.truth].kind);
    for (Tally* t : {&r.overall, &r.by_kind[kind], &r.by_size[size_bin(truth[m.truth].length)]}) {
      t->tp_base++;
      t->tp_call++;
      t->genotype += m.genotype;
    }
  }
  for (std::size_t i = 0; i < calls.size(); i++)
    if (counted_call(calls[i]) && !call_used[i]) {
      r.false_positives.push_back(i);
      for (Tally* t : {&r.overall, &r.by_kind[detail::kind_class(calls[i].kind)], &r.by_size[size_bin(calls[i].length)]}) t->fp++;
    }
  for (std::size_t j = 0; j < truth.size(); j++)
    if (counted_truth(truth[j]) && !truth_used[j]) {
      r.false_negatives.push_back(j);
      for (Tally* t : {&r.overall, &r.by_kind[detail::kind_class(truth[j].kind)], &r.by_size[size_bin(truth[j].length)]}) t->fn++;
    }
  return r;
}

// ---- copy number and LOH, base by base -------------------------------------------------------------
struct CopyConcordance {
  std::size_t bins = 0, agree = 0;
  std::map<std::pair<int, int>, std::size_t> confusion;  // (truth, called) → bins
  double fraction() const { return bins ? static_cast<double>(agree) / static_cast<double>(bins) : 0; }
};

// `truth(contig, position)` → the true copy number.
template <class Truth>
CopyConcordance copy_concordance(const Depth& d, const std::vector<CopySegmentCall>& segments, Truth&& truth) {
  CopyConcordance r;
  for (const auto& b : d.bins) {
    if (!b.mappable) continue;
    const std::size_t mid = (b.begin + b.end) / 2;
    const int t = truth(b.contig, mid), c = copies_at(segments, b.contig, mid);
    r.bins++;
    r.agree += t == c;
    r.confusion[{t, c}]++;
  }
  return r;
}

struct OverlapTally {
  std::size_t truth = 0, called = 0, shared = 0;
  double recall() const { return truth ? static_cast<double>(shared) / static_cast<double>(truth) : 1; }
  double precision() const { return called ? static_cast<double>(shared) / static_cast<double>(called) : 1; }
};

// Bases of called intervals against true ones (each list on one or more contigs, disjoint within itself).
inline OverlapTally base_overlap(const std::vector<std::array<std::size_t, 3>>& called, const std::vector<std::array<std::size_t, 3>>& truth) {
  OverlapTally r;
  for (const auto& c : called) r.called += c[2] - c[1];
  for (const auto& t : truth) r.truth += t[2] - t[1];
  for (const auto& c : called)
    for (const auto& t : truth)
      if (c[0] == t[0]) {
        const std::size_t lo = std::max(c[1], t[1]), hi = std::min(c[2], t[2]);
        if (hi > lo) r.shared += hi - lo;
      }
  return r;
}

}  // namespace genomics::structural
