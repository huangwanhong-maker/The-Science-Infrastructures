// A minimap2-like mapper of reads to a genome: minimizer anchors, chains, base-level alignment
// of the chain, and the read's ends extended and soft-clipped (Li 2018, "Minimap2: pairwise
// alignment for nucleotide sequences", Bioinformatics 34, 3094–3100, §2.1–2.2).
//
// For a read of length L:
//
// 1. Seeds. The read's (w, k)-minimizers are looked up in the genome's minimizer index
//    (index/minimizer_index.hpp, repetitive minimizers ignored as minimap2's -f), and each
//    anchor becomes a forward seed of k bases on the strand mapped (chain.hpp). Or, seeding by
//    SMEMs (BWA-MEM's idea; Li 2013): the read's super-maximal exact matches of at least
//    19 bases (-k) in the genome's FMD-index (index/fmd.hpp), each occurrence of one occurring at
//    most 500 times (-c) a seed of the match's length, without BWA-MEM's re-seeding of long
//    SMEMs. Seeds of either kind are chained alike: α counts the bases a seed adds whatever its
//    span, and w̄ is the seeds' mean span.
// 2. Chains. The seeds are chained by minimap2's DP with its gap cost, backtracked, and divided
//    into primary and secondary chains by their overlap on the read (chain.hpp).
// 3. Base-level alignment of each kept chain, on the read's strand (the reverse complement for a
//    chain on the reverse strand). Each seed is an exact match of its span. They are taken in
//    order: a seed that starts at or after the current position of the alignment is joined to
//    it by a global alignment of the gap between them, then matched; a seed overlapping the
//    current position on the same diagonal extends the run of matches; any other seed is
//    passed over (minimap2 likewise skips anchors that would make the path go back). A gap whose
//    two sides have one length n is matched ungapped when that path scores at least
//    a (n − 1) − 2 (o + e) (a the best pair score): a path with gaps has an insertion and a
//    deletion at least, so it cannot score more, and the ungapped path is optimal (most gaps of
//    short reads, a mismatch between two seeds, end here). Other gaps are
//    filled by the gap-affine wavefront algorithm (align/wavefront.hpp, Marco-Sola et al. 2021,
//    on the scoring's penalty form; Eizenga and Paten 2022), or by doubling bands certified
//    optimal (align/banded.hpp) when it has an N, when its two sides differ by more than
//    `gap_wfa_limit` bases (a long indel, which the band's O(m · |d|) handles in less than the
//    wavefront's O(s (m + n))), or when the scoring has no penalty form.
// 4. Ends. Before the first seed and after the last, the read is extended by the X-drop rule with
//    minimap2's Z-drop (align/xdrop.hpp; Li 2018), over a window of the
//    target a quarter longer than the read's end plus 32 bases; what the extension does not
//    reach is soft-clipped. With a clipping penalty p ≥ 0 (BWA-MEM's -L, 5 for short reads;
//    Li 2013, arXiv:1303.3997), the end is also aligned to the read's last base (a
//    semiglobal alignment of the read's end, with the seed's k bases on the inner side forcing
//    the path through the seed, in a band of 16 + c/8 diagonals for an end of c bases), and kept
//    when its score exceeds the X-drop's best less p: a mismatch near the end is then aligned
//    rather than clipped.
// 5. Mapping quality. The chain's (chain.hpp) is refined by the alignments: among a primary and
//    its kept secondaries, the one of the best alignment score becomes primary (minimap2 sorts
//    its hits by DP score after alignment likewise), and with s₁ ≥ s₂ the best two alignment
//    scores and x_c the chains' ratio f₂/f₁ (f₂ at least the minimum chain score),
//    mapQ = 40 · (1 − x_c · s₂/s₁) · min{1, m/10} · ln f₁, held to [0, 60]: the product x_c ·
//    s₂/s₁ is how minimap2's implementation (mm_set_mapq) combines the two ratios; without an
//    aligned secondary s₂/s₁ is taken as 1, which is Li's formula. Reads in an exact repeat,
//    whose copies give equal chains and equal alignments, get 0. Seeded by SMEMs, a chain has
//    one to a few seeds and min{1, m/10} would cap every read's quality; the quality is then
//    BWA-MEM's approximation from the alignment scores (bwamem.c, mem_approx_mapq_se, without
//    its length coefficient and its count of sub-optimal hits): with s₂ at least 19 times the
//    match score a, ℓ the longer side of the alignment and ι = 1 − (ℓa − s₁)/((a + b)ℓ) its
//    identity (b the mismatch penalty), mapQ = 6.02 (s₁ − s₂)/a · ι², held to [0, 60]. The 6.02
//    is 10 log₁₀ 4: Proposition 3.2 of the theory note with a match four times likelier on the
//    right placement than on a random one.
//
// Reads are mapped in parallel by calx::scheduler::parallel_for, each into its own slot, so that
// the output is in read order and the same on any number of threads. The mapping of one read
// is a pure function of the read, the genome and the options.
//
// Presets: short_reads() — k = 21, w = 11, minimap2's sr scoring (2, −8, 12, 2), chains of two
// seeds and score 25 within G = 200 and a band of 100 (chosen here for reads of 100–300 bases),
// X-drop and Z-drop 100 (BWA-MEM's -d), clipping penalty 5. short_reads_smem() — the same with
// SMEM seeds, chains of one seed and score 25. long_reads() — k = 15, w = 10 as
// minimap2's map-ont, its scoring (2, −4, 4, 2), chains of three seeds and score 40 within
// G = 5000 and a band of 500 (minimap2's -n, -m, -g, -r), X-drop and Z-drop 400 (-z), no
// clipping penalty.
//
// Limits. One affine gap, the first piece of minimap2's two-piece cost; no Z-drop inside the
// gaps between seeds, so an inversion or a large structural variant inside one chain is aligned
// through rather than split; no rescue of mates (reads are single-end); no re-seeding of SMEMs.
// Complexity per read: O(L) minimizers (or O(L) bidirectional extensions for the SMEMs, and
// O(r) steps of LF to locate each occurrence at the sample rate r), O(N h) chaining of N seeds,
// and the alignment of the gaps between seeds, O(Σ s_g (m_g + n_g)) for gaps of penalty s_g,
// which is O(L e) for a read with an error rate e.

#pragma once

#include <genomics/align/banded.hpp>
#include <genomics/align/scoring.hpp>
#include <genomics/align/types.hpp>
#include <genomics/align/wavefront.hpp>
#include <genomics/align/xdrop.hpp>
#include <genomics/index/fmd.hpp>
#include <genomics/index/minimizer_index.hpp>
#include <genomics/mapping/chain.hpp>
#include <genomics/sequence/alphabet.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace genomics::mapping {

using sequence::Code;
using sequence::Codes;

enum class Seeding { minimizers, smems };

struct MapperOptions {
  Seeding seeding = Seeding::minimizers;
  index::MinimizerOptions index;
  std::size_t smem_min_length = 19;   // SMEMs shorter than this give no seeds (BWA-MEM's -k)
  std::uint64_t smem_max_occ = 500;   // nor do those occurring more often (-c)
  ChainOptions chain;
  align::Scoring scoring = align::minimap2_map_ont();
  align::Drop drop{400, 400};
  int clip_penalty = -1;            // < 0: the ends by the X-drop alone
  bool secondary = true;            // report secondary mappings
  std::int64_t gap_wfa_limit = 64;  // gaps whose sides differ by more go to the band

  static MapperOptions short_reads() {
    MapperOptions o;
    o.index.k = 21;
    o.index.w = 11;
    o.chain.max_gap = 200;
    o.chain.bandwidth = 100;
    o.chain.min_anchors = 2;
    o.chain.min_score = 25;
    o.scoring = align::minimap2_sr();
    o.drop = {100, 100};
    o.clip_penalty = 5;
    return o;
  }
  static MapperOptions short_reads_smem() {
    MapperOptions o = short_reads();
    o.seeding = Seeding::smems;
    o.chain.min_anchors = 1;
    return o;
  }
  static MapperOptions long_reads() {
    MapperOptions o;
    o.index.k = 15;
    o.index.w = 10;
    return o;
  }
};

// One alignment of a read: where, on which strand, how well, and the path.
struct Mapping {
  std::uint32_t read = 0;     // the read's index in its batch
  std::uint32_t contig = 0;
  std::size_t target_begin = 0, target_end = 0;
  bool reverse = false;
  int mapq = 0;
  bool primary = true;        // false: a secondary alignment
  bool supplementary = false; // a primary of another part of the read
  std::size_t query_length = 0;
  align::Alignment alignment;  // of the read on the strand mapped (its reverse complement when reverse), without clipping
  std::size_t nm = 0;          // mismatches and gap columns
  std::int64_t chain_score = 0;
  std::uint32_t seeds = 0;     // seeds of the chain

  std::size_t clip_left() const { return alignment.query_begin; }
  std::size_t clip_right() const { return query_length - alignment.query_end; }
  // The CIGAR with the soft clips, as SAM writes it.
  align::Cigar cigar() const {
    align::Cigar c;
    align::push(c, 'S', static_cast<std::uint32_t>(clip_left()));
    for (const auto& o : alignment.cigar) align::push(c, o.op, o.length);
    align::push(c, 'S', static_cast<std::uint32_t>(clip_right()));
    return c;
  }
};

class Mapper {
 public:
  Mapper(calx::scheduler::ThreadPool& pool, std::vector<Codes> contigs, MapperOptions options = MapperOptions::long_reads())
      : options_(std::move(options)),
        contigs_(std::move(contigs)),
        index_(options_.seeding == Seeding::minimizers ? index::MinimizerIndex(pool, contigs_, options_.index) : index::MinimizerIndex()),
        fmd_(options_.seeding == Seeding::smems ? index::FmdIndex(contigs_) : index::FmdIndex()) {
    try {
      (void)align::penalty_form(options_.scoring);
      wfa_ = true;
    } catch (const std::invalid_argument&) {
      wfa_ = false;
    }
  }

  const MapperOptions& options() const { return options_; }
  const std::vector<Codes>& contigs() const { return contigs_; }
  const index::MinimizerIndex& minimizer_index() const { return index_; }
  const index::FmdIndex& fmd_index() const { return fmd_; }

  // The seeds of a read, on the strand mapped, sorted for chaining.
  std::vector<Seed> seeds(const Codes& read) const {
    if (options_.seeding == Seeding::minimizers) return seeds_of(index_.anchors(read), read.size(), options_.index.k);
    std::vector<Seed> out;
    const auto L = static_cast<std::int64_t>(read.size());
    for (const auto& m : fmd_.smems(read, options_.smem_min_length)) {
      if (m.interval.s > options_.smem_max_occ) continue;
      for (const auto& h : fmd_.locate(m.interval, m.length())) {
        if (fmd_.overlaps_ambiguous(h, m.length())) continue;
        const std::int64_t q = h.reverse ? L - static_cast<std::int64_t>(m.end) : static_cast<std::int64_t>(m.begin);
        out.push_back({h.contig, h.reverse, static_cast<std::int64_t>(h.position), q, static_cast<std::int32_t>(m.length())});
      }
    }
    std::sort(out.begin(), out.end());
    return out;
  }

  // The mappings of one read: the primary first, then supplementary parts, then secondaries in
  // decreasing alignment score; empty when the read is unmapped.
  std::vector<Mapping> map(const Codes& read, std::uint32_t id = 0) const {
    std::vector<Mapping> out;
    const std::size_t shortest = options_.seeding == Seeding::minimizers ? static_cast<std::size_t>(options_.index.k) : options_.smem_min_length;
    if (read.size() < shortest) return out;
    std::vector<Chain> chains = chain(seeds(read), options_.chain);
    if (chains.empty()) return out;
    select_chains(chains, options_.chain);
    Codes rc;
    for (const auto& c : chains)
      if (c.reverse) { rc = sequence::reverse_complement(read); break; }
    std::vector<Mapping> aligned(chains.size());
    for (std::size_t c = 0; c < chains.size(); c++) {
      aligned[c] = align_chain(chains[c], chains[c].reverse ? rc : read);
      aligned[c].read = id;
    }
    // Groups: a primary chain and its secondaries.
    std::vector<Mapping> secondaries;
    for (std::size_t p = 0; p < chains.size(); p++) {
      if (!chains[p].primary) continue;
      std::vector<std::size_t> group{p};
      for (std::size_t c = p + 1; c < chains.size(); c++)
        if (!chains[c].primary && chains[c].parent == static_cast<int>(p)) group.push_back(c);
      // The best alignment heads the group; ties stay with the best chain.
      std::size_t head = p;
      for (const std::size_t c : group)
        if (aligned[c].alignment.score > aligned[head].alignment.score) head = c;
      // Two chains whose alignments land on the same bases are one placement.
      const auto same = [&](std::size_t c) {
        const Mapping &x = aligned[c], &y = aligned[head];
        return c == head || (x.contig == y.contig && x.reverse == y.reverse && x.target_begin == y.target_begin && x.target_end == y.target_end);
      };
      std::int64_t f2 = std::max(chains[p].sub_score, options_.chain.min_score), s2 = 0;
      bool aligned_secondary = false;
      for (const std::size_t c : group) {
        if (same(c)) continue;
        f2 = std::max(f2, chains[c].score);
        s2 = std::max(s2, aligned[c].alignment.score);
        aligned_secondary = true;
      }
      const double f1 = static_cast<double>(chains[head].score), s1 = static_cast<double>(aligned[head].alignment.score);
      const double xc = std::clamp(static_cast<double>(f2) / f1, 0.0, 1.0);
      const double xa = aligned_secondary && s1 > 0 ? std::clamp(static_cast<double>(s2) / s1, 0.0, 1.0) : 1.0;
      Mapping m = aligned[head];
      m.primary = true;
      m.supplementary = !out.empty();
      m.mapq = options_.seeding == Seeding::smems ? bwa_quality(m, aligned_secondary ? s2 : 0) : mapping_quality(f1, xc * xa * f1, chains[head].seeds.size());
      out.push_back(std::move(m));
      for (const std::size_t c : group)
        if (!same(c) && options_.secondary) {
          Mapping x = aligned[c];
          x.primary = false;
          x.supplementary = false;
          x.mapq = 0;
          secondaries.push_back(std::move(x));
        }
    }
    std::stable_sort(secondaries.begin(), secondaries.end(), [](const Mapping& a, const Mapping& b) { return a.alignment.score > b.alignment.score; });
    for (auto& x : secondaries) out.push_back(std::move(x));
    return out;
  }

  // The mappings of many reads, in read order, the same on any number of threads. A read is a
  // Codes or anything with a member `bases` (simulate::Read).
  template <class Read>
  std::vector<std::vector<Mapping>> map_all(calx::scheduler::ThreadPool& pool, const std::vector<Read>& reads, std::size_t grain = 16) const {
    std::vector<std::vector<Mapping>> out(reads.size());
    calx::scheduler::parallel_for(pool, reads.size(), grain, [&](std::size_t lo, std::size_t hi, std::size_t) {
      for (std::size_t i = lo; i < hi; i++) {
        if constexpr (std::is_same_v<Read, Codes>) out[i] = map(reads[i], static_cast<std::uint32_t>(i));
        else out[i] = map(reads[i].bases, static_cast<std::uint32_t>(i));
      }
    });
    return out;
  }

  // The alignment of one chain, on the strand mapped (`q` is the read or its reverse complement).
  Mapping align_chain(const Chain& c, const Codes& q) const {
    const Codes& t = contigs_.at(c.contig);
    align::Cigar mid;
    std::int64_t cq = c.seeds.front().query, ct = c.seeds.front().target;
    const std::int64_t q0 = cq, t0 = ct;
    const std::int64_t span0 = c.seeds.front().span;
    std::int64_t last_span = span0;
    for (const auto& s : c.seeds) {
      if (s.query >= cq && s.target >= ct) {
        fill(mid, q, cq, s.query, t, ct, s.target);
        align::push(mid, 'M', static_cast<std::uint32_t>(s.span));
        cq = s.query + s.span;
        ct = s.target + s.span;
        last_span = s.span;
      } else if (s.query - s.target == cq - ct && s.query + s.span > cq) {
        align::push(mid, 'M', static_cast<std::uint32_t>(s.query + s.span - cq));
        ct += s.query + s.span - cq;
        cq = s.query + s.span;
        last_span = s.span;
      }
    }
    const End left = left_end(q, static_cast<std::size_t>(q0), t, static_cast<std::size_t>(t0), static_cast<std::size_t>(span0));
    const End right = right_end(q, static_cast<std::size_t>(cq), t, static_cast<std::size_t>(ct), static_cast<std::size_t>(last_span));
    Mapping m;
    m.contig = c.contig;
    m.reverse = c.reverse;
    m.query_length = q.size();
    m.chain_score = c.score;
    m.seeds = static_cast<std::uint32_t>(c.seeds.size());
    align::Alignment& a = m.alignment;
    a.query_begin = left.query;
    a.target_begin = left.target;
    a.query_end = right.query;
    a.target_end = right.target;
    a.cigar = left.cigar;
    for (const auto& o : mid) align::push(a.cigar, o.op, o.length);
    for (const auto& o : right.cigar) align::push(a.cigar, o.op, o.length);
    a.score = align::rescore(a, q, t, options_.scoring);
    m.nm = align::edit_distance(a, q, t);
    m.target_begin = a.target_begin;
    m.target_end = a.target_end;
    return m;
  }

 private:
  // An end of the alignment: its path and the outer vertex it reaches.
  struct End {
    align::Cigar cigar;
    std::size_t query = 0, target = 0;
  };

  // BWA-MEM's approximate mapping quality of a primary alignment against its best rival's score.
  int bwa_quality(const Mapping& m, std::int64_t s2) const {
    const double a = options_.scoring.matrix(0, 0), b = -static_cast<double>(options_.scoring.matrix(0, 1));
    const double s1 = static_cast<double>(m.alignment.score);
    const double sub = std::max(static_cast<double>(s2), static_cast<double>(options_.smem_min_length) * a);
    if (!(a > 0) || sub >= s1) return 0;
    const double l = static_cast<double>(std::max(m.alignment.query_end - m.alignment.query_begin, m.alignment.target_end - m.alignment.target_begin));
    const double identity = 1 - (l * a - s1) / ((a + b) * l);
    return std::clamp(static_cast<int>(std::floor(6.02 * (s1 - sub) / a * identity * identity + 0.499)), 0, 60);
  }

  static Codes slice(const Codes& x, std::size_t b, std::size_t e) {
    return Codes(x.begin() + static_cast<std::ptrdiff_t>(b), x.begin() + static_cast<std::ptrdiff_t>(e));
  }

  // The gap q[qa, qb) × t[ta, tb) between two seeds, aligned globally.
  void fill(align::Cigar& out, const Codes& q, std::int64_t qa, std::int64_t qb, const Codes& t, std::int64_t ta, std::int64_t tb) const {
    const std::int64_t lq = qb - qa, lt = tb - ta;
    if (lq == 0 && lt == 0) return;
    if (lq == 0) { align::push(out, 'D', static_cast<std::uint32_t>(lt)); return; }
    if (lt == 0) { align::push(out, 'I', static_cast<std::uint32_t>(lq)); return; }
    if (lq == lt) {
      // Sides of one length: a path with gaps has an I and a D at least, so it scores at most
      // a (n − 1) − 2 (o + e) for the best pair score a; the ungapped path, when it reaches that,
      // is optimal.
      const auto& s = options_.scoring;
      std::int64_t ungapped = 0;
      for (std::int64_t k = 0; k < lq; k++) ungapped += s.matrix(q[static_cast<std::size_t>(qa + k)], t[static_cast<std::size_t>(ta + k)]);
      const std::int64_t bound = static_cast<std::int64_t>(s.matrix.max_score()) * (lq - 1) - 2 * (static_cast<std::int64_t>(s.gap_open) + s.gap_extend);
      if (ungapped >= bound) { align::push(out, 'M', static_cast<std::uint32_t>(lq)); return; }
    }
    const Codes qs = slice(q, static_cast<std::size_t>(qa), static_cast<std::size_t>(qb));
    const Codes ts = slice(t, static_cast<std::size_t>(ta), static_cast<std::size_t>(tb));
    const bool clean = sequence::count_ambiguous(qs) + sequence::count_ambiguous(ts) == 0;
    const std::int64_t d = lq > lt ? lq - lt : lt - lq;
    const align::Alignment g = (wfa_ && clean && d <= options_.gap_wfa_limit) ? align::wavefront_global(qs, ts, options_.scoring)
                                                                              : align::banded_global(qs, ts, options_.scoring, 16);
    for (const auto& o : g.cigar) align::push(out, o.op, o.length);
  }

  // The left end: the read's bases q[0, qs) before the first seed, which starts at (qs, ts) and
  // spans `span` bases.
  End left_end(const Codes& q, std::size_t qs, const Codes& t, std::size_t ts, std::size_t span) const {
    if (qs == 0) return {{}, 0, ts};
    const std::size_t c = qs, slack = c / 4 + 32;
    const std::size_t wb = ts - std::min(ts, c + slack);
    const Codes window = slice(t, wb, ts + span);
    const align::Alignment loc = align::extend_left(q, window, qs, ts - wb, options_.scoring, options_.drop);
    End best{loc.cigar, loc.query_begin, wb + loc.target_begin};
    if (options_.clip_penalty >= 0 && loc.query_begin > 0) {
      const Codes Q = slice(q, 0, qs + span);
      const std::int64_t n = static_cast<std::int64_t>(window.size());
      const std::int64_t d_end = n - static_cast<std::int64_t>(qs + span), b = 16 + static_cast<std::int64_t>(c / 8);
      try {
        align::Alignment e = align::banded(Q, window, options_.scoring, align::Mode::semiglobal, d_end - b, d_end + b);
        if (e.target_end == window.size() && !e.cigar.empty() && e.cigar.back().op == 'M' && e.cigar.back().length >= span) {
          e.cigar.back().length -= static_cast<std::uint32_t>(span);
          if (e.cigar.back().length == 0) e.cigar.pop_back();
          e.query_end -= span;
          e.target_end -= span;
          if (align::rescore(e, Q, window, options_.scoring) > loc.score - options_.clip_penalty) best = {e.cigar, e.query_begin, wb + e.target_begin};
        }
      } catch (const std::invalid_argument&) {
        // The band missed the edit graph (the read overhangs the contig's start): keep the X-drop's.
      }
    }
    return best;
  }

  // The right end: the read's bases q[qe, L) after the last seed, which ends at (qe, te) and
  // spans `span` bases.
  End right_end(const Codes& q, std::size_t qe, const Codes& t, std::size_t te, std::size_t span) const {
    if (qe == q.size()) return {{}, qe, te};
    const std::size_t c = q.size() - qe, slack = c / 4 + 32;
    const std::size_t we = std::min(t.size(), te + c + slack), wb = te - span;
    const Codes window = slice(t, wb, we);
    const align::Alignment loc = align::extend_right(q, window, qe, span, options_.scoring, options_.drop);
    End best{loc.cigar, loc.query_end, wb + loc.target_end};
    if (options_.clip_penalty >= 0 && loc.query_end < q.size()) {
      const Codes Q = slice(q, qe - span, q.size());
      const std::int64_t b = 16 + static_cast<std::int64_t>(c / 8);
      try {
        align::Alignment e = align::banded(Q, window, options_.scoring, align::Mode::semiglobal, -b, b);
        if (e.target_begin == 0 && !e.cigar.empty() && e.cigar.front().op == 'M' && e.cigar.front().length >= span) {
          e.cigar.front().length -= static_cast<std::uint32_t>(span);
          if (e.cigar.front().length == 0) e.cigar.erase(e.cigar.begin());
          e.query_begin += span;
          e.target_begin += span;
          if (align::rescore(e, Q, window, options_.scoring) > loc.score - options_.clip_penalty) best = {e.cigar, q.size(), wb + e.target_end};
        }
      } catch (const std::invalid_argument&) {
        // No path of the band reaches the read's end inside the contig: keep the X-drop's.
      }
    }
    return best;
  }

  MapperOptions options_;
  std::vector<Codes> contigs_;
  index::MinimizerIndex index_;
  index::FmdIndex fmd_;
  bool wfa_ = false;
};

}  // namespace genomics::mapping
