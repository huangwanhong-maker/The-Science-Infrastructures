// The FMD-index of a genome: an FM-index of its contigs and of their reverse complements together,
// which, because that set of words is closed under reverse complementation, can extend a match
// in both directions (Li 2012, "Exploring single-sample SNP and INDEL calling with whole-genome de
// novo assembly", Bioinformatics 28, 1838–1844, §§2.3–2.4 and 2.6, Algorithms 2, 3 and 5). It is the seed index of
// BWA-MEM: the super-maximal exact matches of a read.
//
// - The text. T = c₀ $₀ c̄₀ $₁ c₁ $₂ c̄₁ $₃ …, the contigs c_i and their reverse complements c̄_i,
//   each closed by a sentinel; the sentinels are distinct and smaller than every base, ordered by
//   position, as Li's (§2.3.1). A base N is replaced, before the reverse complement is taken, by a
//   pseudo-random base drawn from a stream keyed by the seed, FMD_N and the contig, as BWA's
//   index replaces its holes; the runs replaced are recorded, and a hit can be tested against
//   them. The suffix array is morphe's SA-IS (morphe::formal::suffix_array), O(|T|), over the
//   alphabet of the 2m sentinels and the four bases; T is held as bytes when 2m + 4 ≤ 256.
// - The BWT B over {$, A, C, G, T} is kept as BWA keeps it (bwt.h, bwt_occ): in blocks of 128
//   symbols, each block four 64-bit counts of A, C, G, T before it and four 64-bit words of 2-bit
//   symbols, 64 bytes, one cache line, 4 bits a symbol. A $ is stored as A and its rows (2m of
//   them) are kept sorted apart; Occ(A, i) subtracts those before i. Occ(c, i) is a count from the
//   block and popcounts of the masks of equal 2-bit groups: O(1), and Occ of the four bases at
//   once costs four popcounts a word.
// - Backward search: the rows of cP are [C(c) + Occ(c, k), C(c) + Occ(c, k + s)) when those of
//   P are [k, k + s) (Ferragina and Manzini 2000). With C(c) the number of symbols of B smaller
//   than c, and the 2m sentinel rows first.
// - Bidirectional search. A bi-interval (k, l, s) holds the rows [k, k + s) of P and
//   [l, l + s) of its reverse complement P̄, which occur equally often since T is closed under
//   reverse complementation. Extending P backwards by a gives the rows of aP from Occ as above,
//   and those of P̄ā are the sub-interval of P̄'s rows followed by ā; the rows of P̄ split by
//   their next symbol in the order $ < A < C < G < T, which is the order T, G, C, A, $ of the
//   bases a that precede P (Li 2012, Algorithm 2):
//     l_$ = l,  l_T = l_$ + s_$,  l_G = l_T + s_T,  l_C = l_G + s_G,  l_A = l_C + s_C,
//   with s_$ = s − Σ_b s_b the occurrences of P at the start of a sequence. Extending P forwards by
//   a is extending P̄ backwards by ā, with k and l exchanged (Algorithm 3).
// - Super-maximal exact matches (SMEMs). A maximal exact match of a query Q is a factor Q[b, e)
//   occurring in T that cannot be extended to either side; an SMEM is one contained in no other
//   (Li 2012, §2.6). Algorithm 5 finds every SMEM containing a position x: extend Q[x] forwards,
//   keeping each bi-interval at which the count drops; then extend all of them backwards in
//   step, each surviving while it occurs, and output a match when it can be extended no further
//   and no longer one survived. Starting again from the end of the longest match from x finds
//   every SMEM of Q, each once. O(|Q|) extensions per x in the worst case; in practice O(|Q|)
//   for the whole query when matches are long. This is BWA-MEM's bwt_smem1 without re-seeding.
// - Locate. Suffix-array values are sampled at the text positions divisible by the sample rate,
//   and at the start of every sequence (the rows whose symbol is $), marked in a bit vector with
//   rank (morphe::formal::BitVector); LF never leaves a row of $, and a position is recovered in
//   at most rate − 1 steps. A text position becomes a hit (contig, position, reverse): on c̄_i at
//   offset o, a match of length m is the reverse strand of c_i at |c_i| − o − m.
//
// Space: |T|/2 bytes for B and its counts, |T|/8 · 1.25 bytes of marks, 4 bytes a sample. For
// Escherichia coli's 4.64 Mb, |T| = 9.28 M: 7.3 MB at rate 32. Construction holds the 32-bit
// suffix array (4|T| bytes) and SA-IS's work arrays. The text must be shorter than 2^32 − 1
// (a genome below 2.1 Gb); larger genomes need 64-bit indices and an external-memory
// construction, not written. Exact: no probabilistic step, apart from the documented bases that
// stand for N.

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <morphe/formal/succinct.hpp>
#include <morphe/formal/suffix_array.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace genomics::index {

using sequence::Code;
using sequence::Codes;

// The purpose of the stream that draws the bases standing for N.
inline constexpr std::uint32_t FMD_N = 0x464d444e;

// Rows [k, k + s) of a pattern P and [l, l + s) of its reverse complement.
struct BiInterval {
  std::uint64_t k = 0, l = 0, s = 0;
  friend bool operator==(const BiInterval&, const BiInterval&) = default;
};

// An occurrence: on the forward strand of a contig at [position, position + length), or, when
// reverse, the reverse complement of the pattern there.
struct Hit {
  std::uint32_t contig = 0;
  std::uint64_t position = 0;
  bool reverse = false;
  friend bool operator==(const Hit&, const Hit&) = default;
  friend bool operator<(const Hit& a, const Hit& b) {
    return a.contig != b.contig ? a.contig < b.contig : a.position != b.position ? a.position < b.position : a.reverse < b.reverse;
  }
};

// A super-maximal exact match Q[begin, end), with its bi-interval (s occurrences in T).
struct Smem {
  std::uint32_t begin = 0, end = 0;
  BiInterval interval;
  std::uint32_t length() const { return end - begin; }
  friend bool operator==(const Smem&, const Smem&) = default;
};

// A run of N in a contig, [begin, end), replaced in the index by drawn bases.
struct AmbiguousRun {
  std::uint32_t contig = 0;
  std::uint64_t begin = 0, end = 0;
  friend bool operator==(const AmbiguousRun&, const AmbiguousRun&) = default;
};

struct FmdOptions {
  std::uint32_t sample = 32;   // the suffix-array sample rate
  std::uint64_t seed = 11;     // the stream of bases for N
};

class FmdIndex {
 public:
  static constexpr Code DOLLAR = 4;

  FmdIndex() = default;

  explicit FmdIndex(const std::vector<Codes>& contigs, FmdOptions options = {}) : options_(options) {
    if (options.sample == 0) throw std::invalid_argument("genomics::index::FmdIndex: a sample rate of at least 1");
    if (contigs.empty()) throw std::invalid_argument("genomics::index::FmdIndex: at least one contig");
    // The contigs with N replaced, as the index sees them.
    std::uint64_t total = 0;
    for (const auto& c : contigs) total += 2 * (c.size() + 1);
    if (total >= std::numeric_limits<std::uint32_t>::max()) throw std::invalid_argument("genomics::index::FmdIndex: the text exceeds 32-bit indices");
    lengths_.reserve(contigs.size());
    for (std::uint32_t i = 0; i < contigs.size(); i++) lengths_.push_back(contigs[i].size());
    const std::size_t m = contigs.size(), seqs = 2 * m;
    starts_.reserve(seqs + 1);
    if (seqs + 4 <= 256) build<std::uint8_t>(contigs, total);
    else build<std::uint32_t>(contigs, total);
  }

  // ---- the genome -------------------------------------------------------------------------------------------
  std::size_t contigs() const { return lengths_.size(); }
  std::uint64_t contig_length(std::size_t i) const { return lengths_.at(i); }
  // |T|, the rows of the BWT.
  std::uint64_t size() const { return n_; }
  const std::vector<AmbiguousRun>& ambiguous() const { return ambiguous_; }
  const FmdOptions& options() const { return options_; }
  // Where each sequence (c₀, c̄₀, c₁, …) begins in T; the last entry is |T|.
  const std::vector<std::uint64_t>& sequence_starts() const { return starts_; }

  // Whether [position, position + length) of a hit's contig meets a run of N.
  bool overlaps_ambiguous(const Hit& h, std::uint64_t length) const {
    auto it = std::lower_bound(ambiguous_.begin(), ambiguous_.end(), h.contig, [](const AmbiguousRun& r, std::uint32_t c) { return r.contig < c; });
    for (; it != ambiguous_.end() && it->contig == h.contig && it->begin < h.position + length; ++it)
      if (it->end > h.position) return true;
    return false;
  }

  // ---- the BWT and its occurrences --------------------------------------------------------------------------------
  // The symbol of B at a row: 0..3 for a base, DOLLAR for a sentinel.
  Code bwt(std::uint64_t i) const {
    if (is_dollar(i)) return DOLLAR;
    const std::uint64_t* b = &occ_[(i >> 7) * 8];
    return static_cast<Code>((b[4 + ((i & 127) >> 5)] >> (2 * (i & 31))) & 3);
  }
  // C(c): the symbols of B smaller than c (c = 0..4, C(4) = |T|).
  std::uint64_t C(Code c) const { return c_[c]; }

  // Occ(c, i): the occurrences of base c in B[0, i), 0 ≤ i ≤ |T|.
  std::uint64_t occ(Code c, std::uint64_t i) const {
    const std::uint64_t* b = &occ_[(i >> 7) * 8];
    std::uint64_t n = b[c];
    const unsigned r = static_cast<unsigned>(i & 127);
    const std::uint64_t pattern = static_cast<std::uint64_t>(c) * 0x5555555555555555ULL;
    for (unsigned w = 0; w < (r >> 5); w++) n += equal_groups(b[4 + w] ^ pattern, ODD);
    if (r & 31) n += equal_groups(b[4 + (r >> 5)] ^ pattern, ODD & ((1ULL << (2 * (r & 31))) - 1));
    if (c == 0) n -= dollars_before(i);
    return n;
  }
  // Occ of the four bases at once.
  std::array<std::uint64_t, 4> occ4(std::uint64_t i) const {
    const std::uint64_t* b = &occ_[(i >> 7) * 8];
    std::array<std::uint64_t, 4> n{b[0], b[1], b[2], b[3]};
    const unsigned r = static_cast<unsigned>(i & 127);
    const auto add = [&](std::uint64_t x, std::uint64_t mask) {
      const std::uint64_t lo = x & ODD, hi = (x >> 1) & ODD;
      n[0] += static_cast<std::uint64_t>(std::popcount(~hi & ~lo & mask));
      n[1] += static_cast<std::uint64_t>(std::popcount(~hi & lo & mask));
      n[2] += static_cast<std::uint64_t>(std::popcount(hi & ~lo & mask));
      n[3] += static_cast<std::uint64_t>(std::popcount(hi & lo & mask));
    };
    for (unsigned w = 0; w < (r >> 5); w++) add(b[4 + w], ODD);
    if (r & 31) add(b[4 + (r >> 5)], ODD & ((1ULL << (2 * (r & 31))) - 1));
    n[0] -= dollars_before(i);
    return n;
  }

  // LF(i) for a row whose symbol is a base: the row of the suffix one position earlier.
  std::uint64_t lf(std::uint64_t i) const {
    const Code c = bwt(i);
    return c_[c] + occ(c, i);
  }

  // ---- bidirectional search --------------------------------------------------------------------------------
  // The bi-interval of one base.
  BiInterval init(Code a) const {
    if (a > 3) return {};
    return {c_[a], c_[3 - a], c_[a + 1] - c_[a]};
  }

  // The bi-interval of aP from that of P (Li 2012, Algorithm 2).
  BiInterval backward(const BiInterval& x, Code a) const {
    if (a > 3 || x.s == 0) return {};
    const auto lo = occ4(x.k), hi = occ4(x.k + x.s);
    std::array<std::uint64_t, 4> k, s;
    for (int b = 0; b < 4; b++) {
      k[b] = c_[b] + lo[b];
      s[b] = hi[b] - lo[b];
    }
    const std::uint64_t dollar = x.s - (s[0] + s[1] + s[2] + s[3]);
    std::array<std::uint64_t, 4> l;
    l[3] = x.l + dollar;
    l[2] = l[3] + s[3];
    l[1] = l[2] + s[2];
    l[0] = l[1] + s[1];
    return {k[a], l[a], s[a]};
  }

  // The bi-interval of Pa (Li 2012, Algorithm 3): P̄ extended backwards by ā.
  BiInterval forward(const BiInterval& x, Code a) const {
    if (a > 3) return {};
    const BiInterval y = backward({x.l, x.k, x.s}, static_cast<Code>(3 - a));
    return {y.l, y.k, y.s};
  }

  // The bi-interval of a pattern by backward search; s = 0 when it does not occur or holds N.
  BiInterval search(const Codes& pattern) const {
    if (pattern.empty()) return {0, 0, n_};
    BiInterval x = init(pattern.back());
    for (std::size_t i = pattern.size() - 1; i-- > 0 && x.s > 0;) x = backward(x, pattern[i]);
    return x;
  }
  // The occurrences on both strands.
  std::uint64_t count(const Codes& pattern) const { return pattern.empty() ? 0 : search(pattern).s; }

  // ---- locate --------------------------------------------------------------------------------
  // SA[row]: the text position of the suffix in a row.
  std::uint64_t suffix(std::uint64_t row) const {
    std::uint64_t steps = 0;
    while (!marked_[row]) {
      row = lf(row);
      steps++;
    }
    return samples_[marked_.rank1(row)] + steps;
  }

  // The hit of a match of the given length at a text position.
  Hit hit(std::uint64_t t, std::uint64_t length) const {
    const auto it = std::upper_bound(starts_.begin(), starts_.end(), t);
    const std::size_t q = static_cast<std::size_t>(it - starts_.begin()) - 1;
    const std::uint32_t contig = static_cast<std::uint32_t>(q / 2);
    const std::uint64_t offset = t - starts_[q];
    if (q % 2 == 0) return {contig, offset, false};
    return {contig, lengths_[contig] - offset - length, true};
  }

  // Every occurrence of the rows [k, k + s), sorted; at most `max` of them.
  std::vector<Hit> locate(const BiInterval& x, std::uint64_t length, std::size_t max = SIZE_MAX) const {
    std::vector<Hit> out;
    const std::uint64_t n = std::min<std::uint64_t>(x.s, max);
    out.reserve(n);
    for (std::uint64_t r = x.k; r < x.k + n; r++) out.push_back(hit(suffix(r), length));
    std::sort(out.begin(), out.end());
    return out;
  }
  std::vector<Hit> locate(const Codes& pattern, std::size_t max = SIZE_MAX) const {
    if (pattern.empty()) return {};
    return locate(search(pattern), pattern.size(), max);
  }

  // ---- super-maximal exact matches --------------------------------------------------------------------------------
  // Every SMEM of the query of at least min_length bases occurring at least min_occ times,
  // sorted by begin.
  std::vector<Smem> smems(const Codes& query, std::size_t min_length = 1, std::uint64_t min_occ = 1) const {
    std::vector<Smem> all, found;
    std::vector<Smem> a, b;
    for (std::size_t x = 0; x < query.size();) {
      found.clear();
      x = smems_at(query, x, min_occ, found, a, b);
      for (const auto& m : found)
        if (m.length() >= min_length) all.push_back(m);
    }
    return all;
  }

  // The SMEMs containing x (Li 2012, Algorithm 5), appended to out in order of begin; returns the
  // end of the longest match starting at x (x + 1 when Q[x] is N). prev and curr are scratch.
  std::size_t smems_at(const Codes& q, std::size_t x, std::uint64_t min_occ, std::vector<Smem>& out, std::vector<Smem>& prev, std::vector<Smem>& curr) const {
    const std::size_t n = q.size();
    if (q[x] > 3) return x + 1;
    if (min_occ < 1) min_occ = 1;
    prev.clear();
    Smem ik{static_cast<std::uint32_t>(x), static_cast<std::uint32_t>(x + 1), init(q[x])};
    std::size_t i = x + 1;
    // Forward: keep the bi-interval each time the count is about to drop.
    for (; i < n; i++) {
      if (q[i] > 3) { prev.push_back(ik); break; }
      const BiInterval ok = forward(ik.interval, q[i]);
      if (ok.s != ik.interval.s) {
        prev.push_back(ik);
        if (ok.s < min_occ) break;
      }
      ik.interval = ok;
      ik.end = static_cast<std::uint32_t>(i + 1);
    }
    if (i == n) prev.push_back(ik);
    std::reverse(prev.begin(), prev.end());  // the longest match first
    const std::size_t ret = prev.front().end;
    const std::size_t first = out.size();
    // Backward: extend every candidate in step; one that cannot go further is an SMEM when no
    // longer candidate went further at this step.
    for (std::ptrdiff_t j = static_cast<std::ptrdiff_t>(x) - 1; j >= -1; j--) {
      const int c = j < 0 ? -1 : q[static_cast<std::size_t>(j)] < 4 ? q[static_cast<std::size_t>(j)] : -1;
      curr.clear();
      for (const Smem& p : prev) {
        BiInterval ok{};
        if (c >= 0) ok = backward(p.interval, static_cast<Code>(c));
        if (c < 0 || ok.s < min_occ) {
          if (curr.empty() && (out.size() == first || static_cast<std::size_t>(j + 1) < out.back().begin))
            out.push_back({static_cast<std::uint32_t>(j + 1), p.end, p.interval});
        } else if (curr.empty() || ok.s != curr.back().interval.s) {
          curr.push_back({static_cast<std::uint32_t>(j), p.end, ok});
        }
      }
      if (curr.empty()) break;
      std::swap(prev, curr);
    }
    std::reverse(out.begin() + static_cast<std::ptrdiff_t>(first), out.end());
    return ret;
  }

  // The SMEMs of many queries on the pool, each query's in its own slot.
  std::vector<std::vector<Smem>> smems(calx::scheduler::ThreadPool& pool, const std::vector<Codes>& queries, std::size_t min_length = 1, std::uint64_t min_occ = 1) const {
    std::vector<std::vector<Smem>> out(queries.size());
    calx::scheduler::parallel_for(pool, queries.size(), 16, [&](std::size_t begin, std::size_t end, std::size_t) {
      for (std::size_t i = begin; i < end; i++) out[i] = smems(queries[i], min_length, min_occ);
    });
    return out;
  }

  std::size_t bytes() const {
    return occ_.size() * 8 + marked_.bytes() + samples_.size() * 4 + dollars_.size() * 8 + starts_.size() * 8 + ambiguous_.size() * sizeof(AmbiguousRun);
  }

 private:
  static constexpr std::uint64_t ODD = 0x5555555555555555ULL;

  // The 2-bit groups of x that are zero, among those the mask selects.
  static std::uint64_t equal_groups(std::uint64_t x, std::uint64_t mask) {
    return static_cast<std::uint64_t>(std::popcount(~(x | (x >> 1)) & mask));
  }
  std::uint64_t dollars_before(std::uint64_t i) const {
    return static_cast<std::uint64_t>(std::lower_bound(dollars_.begin(), dollars_.end(), i) - dollars_.begin());
  }
  bool is_dollar(std::uint64_t i) const { return std::binary_search(dollars_.begin(), dollars_.end(), i); }

  template <class Letter>
  void build(const std::vector<Codes>& contigs, std::uint64_t total) {
    const std::size_t m = contigs.size(), seqs = 2 * m;
    // T: sentinel j is the letter j, base b the letter 2m + b.
    std::vector<Letter> text;
    text.reserve(total);
    Codes forward;
    for (std::uint32_t i = 0; i < m; i++) {
      forward = contigs[i];
      calx::random::Stream r(options_.seed, FMD_N, i, 0);
      for (std::uint64_t p = 0; p < forward.size(); p++) {
        if (forward[p] <= 3) continue;
        if (!ambiguous_.empty() && ambiguous_.back().contig == i && ambiguous_.back().end == p) ambiguous_.back().end = p + 1;
        else ambiguous_.push_back({i, p, p + 1});
        forward[p] = static_cast<Code>(r.below(4));
      }
      starts_.push_back(text.size());
      for (Code b : forward) text.push_back(static_cast<Letter>(seqs + b));
      text.push_back(static_cast<Letter>(2 * i));
      starts_.push_back(text.size());
      for (std::size_t p = forward.size(); p-- > 0;) text.push_back(static_cast<Letter>(seqs + 3 - forward[p]));
      text.push_back(static_cast<Letter>(2 * i + 1));
    }
    starts_.push_back(text.size());
    n_ = text.size();
    const auto sa = morphe::formal::suffix_array<std::uint32_t>(text, seqs + 4);

    // B, packed, with the counts before each block; the rows of $ apart.
    const std::uint64_t blocks = n_ / 128 + 1;
    occ_.assign(blocks * 8, 0);
    std::array<std::uint64_t, 5> count{};
    marked_ = morphe::formal::BitVector(n_);
    for (std::uint64_t i = 0; i < n_; i++) {
      if ((i & 127) == 0)
        for (int c = 0; c < 4; c++) occ_[(i >> 7) * 8 + static_cast<std::uint64_t>(c)] = count[static_cast<std::size_t>(c)] + (c == 0 ? count[4] : 0);
      const std::uint64_t p = sa[i];
      const Letter prior = text[p == 0 ? n_ - 1 : p - 1];
      Code sym;
      if (prior < seqs) {
        sym = DOLLAR;
        dollars_.push_back(i);
        count[4]++;
      } else {
        sym = static_cast<Code>(prior - seqs);
        count[sym]++;
        occ_[(i >> 7) * 8 + 4 + ((i & 127) >> 5)] |= static_cast<std::uint64_t>(sym) << (2 * (i & 31));
      }
      if (p % options_.sample == 0 || sym == DOLLAR) marked_.set(i);
    }
    if ((n_ & 127) == 0)
      for (int c = 0; c < 4; c++) occ_[(n_ >> 7) * 8 + static_cast<std::uint64_t>(c)] = count[static_cast<std::size_t>(c)] + (c == 0 ? count[4] : 0);
    marked_.build();
    samples_.reserve(marked_.ones());
    for (std::uint64_t i = 0; i < n_; i++)
      if (marked_[i]) samples_.push_back(sa[i]);
    // C over $ < A < C < G < T.
    c_[0] = count[4];
    for (int c = 0; c < 4; c++) c_[static_cast<std::size_t>(c) + 1] = c_[static_cast<std::size_t>(c)] + count[static_cast<std::size_t>(c)];
  }

  FmdOptions options_;
  std::uint64_t n_ = 0;
  std::vector<std::uint64_t> lengths_, starts_;
  std::vector<AmbiguousRun> ambiguous_;
  std::vector<std::uint64_t> occ_;       // blocks of 128: counts of A, C, G, T before, then 4 words of symbols
  std::vector<std::uint64_t> dollars_;   // the rows whose symbol is $, increasing
  std::array<std::uint64_t, 5> c_{};     // C(A), C(C), C(G), C(T), |T|
  morphe::formal::BitVector marked_;
  std::vector<std::uint32_t> samples_;
};

}  // namespace genomics::index
