// Whole-genome alignment of two genomes by their maximal unique matches: the chain of MUMs that
// lie in the same order in both genomes, the gaps between its anchors, the collinear blocks, and
// the points of a dot plot (Delcher et al. 1999, Nucleic Acids Research 27, 2369–2376: MUMmer's
// alignment by MUMs, their longest increasing subsequence, and the gaps between them).
//
// - The chain. Sorted by their starts in A, the MUMs of one strand whose starts in B also increase
//   form a chain: an alignment that keeps both orders. The chain of greatest weight, the sum of the
//   MUMs' lengths, is the heaviest increasing subsequence of the MUMs' starts in B, found by
//   morphe's heaviest_chain in O(m log m) for m MUMs (morphe/combinatorics/increasing.hpp,
//   Jacobson and Vo 1992). On the reverse strand B is read on its reverse complement, so that an
//   inverted segment is a chain like any other. MUMs never nest (mums.hpp), so starts increasing in
//   both genomes means ends increasing too; consecutive anchors may still overlap, in a repeat or
//   beside an indel in a homopolymer, and the chain's weight then counts the shared bases twice.
//   *Rationale.* The mapper's chaining (mapping/chain.hpp) does not fit: it is minimap2's DP over
//   the h nearest predecessors with a cut-off G on gaps and a bandwidth on diagonals, made for the
//   seeds of one read; between two genomes the MUMs are few, the gaps between them may be as long
//   as a rearrangement, and the chain wanted is the global one, exactly.
// - The gaps. Between consecutive anchors of a chain lie length_a bases of A and length_b of B
//   (negative when the anchors overlap). Equal lengths of 1 are a single-base substitution (MUMmer's
//   SNP); equal lengths above 1 a stretch of substitutions; unequal lengths an insertion or deletion
//   of |length_b − length_a| bases, with or without substitutions beside it. A gap is told by its
//   lengths alone; its bases are not aligned here (align/ does that, when asked).
// - Blocks. A chain is cut where a gap exceeds max_gap in either genome; each piece is a collinear
//   block, forward or reversed. One chain a strand recovers one collinear backbone a strand: a
//   genome against a copy of itself with an inversion gives two forward blocks and one reversed
//   between them; several rearrangements of one strand need local chains (theory §4.3).
// - The dot plot: each MUM as a segment from (a, b) to (a + length, b + length), or from
//   (a, b + length) to (a + length, b) on the reverse strand, marked when it is chained.
//
// Cost: the MUMs of both strands, O(|A| + |B|); the chains, O(m log m).

#pragma once

#include <genomics/comparative/mums.hpp>

#include <morphe/combinatorics/increasing.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace genomics::comparative {

enum class GapKind { snp, substitutions, indel };

struct Gap {
  std::int64_t a = 0, b = 0;                // where it begins in A, and in B on the chain's strand
  std::int64_t length_a = 0, length_b = 0;  // bases between the anchors; negative where they overlap
  GapKind kind() const { return length_a != length_b ? GapKind::indel : length_a == 1 ? GapKind::snp : GapKind::substitutions; }
  std::int64_t indel() const { return length_b - length_a; }  // > 0: bases in B that A lacks
};

struct Chain {
  bool reverse = false;
  std::vector<Match> anchors;  // increasing a (and b on the chain's strand)
  std::vector<Gap> gaps;       // gaps[i] lies between anchors i and i + 1
  std::int64_t weight = 0;     // the sum of the anchors' lengths
};

// The start of a match in B on its own strand: b, or |B| − b − length on the reverse strand.
inline std::int64_t strand_b(const Match& m, std::size_t length_b) {
  return m.reverse ? static_cast<std::int64_t>(length_b - m.b - m.length) : static_cast<std::int64_t>(m.b);
}

// The heaviest chain of matches of one strand.
inline Chain chain_matches(const std::vector<Match>& matches, std::size_t length_b) {
  Chain c;
  if (matches.empty()) return c;
  c.reverse = matches.front().reverse;
  std::vector<std::int64_t> x, y, w;
  x.reserve(matches.size());
  for (const auto& m : matches) {
    if (m.reverse != c.reverse) throw std::invalid_argument("genomics::comparative::chain_matches: matches of one strand");
    x.push_back(static_cast<std::int64_t>(m.a));
    y.push_back(strand_b(m, length_b));
    w.push_back(static_cast<std::int64_t>(m.length));
  }
  for (const std::size_t i : morphe::combinatorics::heaviest_chain(x, y, w)) {
    c.anchors.push_back(matches[i]);
    c.weight += w[i];
  }
  for (std::size_t i = 1; i < c.anchors.size(); i++) {
    const Match& p = c.anchors[i - 1];
    const Match& q = c.anchors[i];
    Gap g;
    g.a = static_cast<std::int64_t>(p.a_end());
    g.b = strand_b(p, length_b) + static_cast<std::int64_t>(p.length);
    g.length_a = static_cast<std::int64_t>(q.a) - g.a;
    g.length_b = strand_b(q, length_b) - g.b;
    c.gaps.push_back(g);
  }
  return c;
}

struct Block {
  bool reverse = false;
  std::size_t a_begin = 0, a_end = 0, b_begin = 0, b_end = 0;  // the span of its anchors; B forward
  std::size_t anchors = 0, matched = 0;                          // its anchors and their bases
};

// The collinear blocks of a chain: cut where a gap exceeds max_gap in A or in B.
inline std::vector<Block> blocks(const Chain& c, std::int64_t max_gap) {
  std::vector<Block> out;
  for (std::size_t i = 0; i < c.anchors.size(); i++) {
    const Match& m = c.anchors[i];
    const bool cut = i == 0 || c.gaps[i - 1].length_a > max_gap || c.gaps[i - 1].length_b > max_gap;
    if (cut) out.push_back({c.reverse, m.a, m.a_end(), m.b, m.b_end(), 0, 0});
    Block& b = out.back();
    b.a_begin = std::min(b.a_begin, m.a);
    b.a_end = std::max(b.a_end, m.a_end());
    b.b_begin = std::min(b.b_begin, m.b);
    b.b_end = std::max(b.b_end, m.b_end());
    b.anchors++;
    b.matched += m.length;
  }
  return out;
}

struct GenomeOptions {
  std::size_t min_length = 20;   // MUMmer's default minimum length of a match
  std::int64_t max_gap = 1000;   // a longer gap cuts a block
  bool both_strands = true;
};

struct GenomeAlignment {
  std::size_t length_a = 0, length_b = 0;
  std::vector<Match> forward, reverse;  // every MUM of each strand
  Chain forward_chain, reverse_chain;
  std::vector<Block> blocks;            // of both chains, in increasing a_begin

  // The fraction of a genome's bases inside the anchors of the chains (each base once).
  double coverage_a() const { return covered(true) / static_cast<double>(length_a ? length_a : 1); }
  double coverage_b() const { return covered(false) / static_cast<double>(length_b ? length_b : 1); }

 private:
  double covered(bool in_a) const {
    std::vector<std::pair<std::size_t, std::size_t>> iv;
    for (const Chain* c : {&forward_chain, &reverse_chain})
      for (const auto& m : c->anchors) iv.push_back(in_a ? std::make_pair(m.a, m.a_end()) : std::make_pair(m.b, m.b_end()));
    std::sort(iv.begin(), iv.end());
    std::size_t total = 0, end = 0;
    for (const auto& [b, e] : iv) {
      if (e <= end) continue;
      total += e - std::max(b, end);
      end = e;
    }
    return static_cast<double>(total);
  }
};

inline GenomeAlignment align_genomes(const Codes& A, const Codes& B, const GenomeOptions& o = {}) {
  GenomeAlignment g;
  g.length_a = A.size();
  g.length_b = B.size();
  g.forward = maximal_unique_matches(A, B, o.min_length, false);
  if (o.both_strands) g.reverse = maximal_unique_matches(A, B, o.min_length, true);
  g.forward_chain = chain_matches(g.forward, B.size());
  g.reverse_chain = chain_matches(g.reverse, B.size());
  g.reverse_chain.reverse = true;
  for (const Chain* c : {&g.forward_chain, &g.reverse_chain})
    for (const auto& b : blocks(*c, o.max_gap)) g.blocks.push_back(b);
  std::sort(g.blocks.begin(), g.blocks.end(), [](const Block& x, const Block& y) { return x.a_begin < y.a_begin; });
  return g;
}

struct DotSegment {
  std::size_t x0 = 0, y0 = 0, x1 = 0, y1 = 0;  // A on x, B (forward) on y
  bool reverse = false, chained = false;
};

inline std::vector<DotSegment> dot_plot(const GenomeAlignment& g) {
  std::vector<DotSegment> out;
  for (const Chain* c : {&g.forward_chain, &g.reverse_chain}) {
    const std::vector<Match>& all = c == &g.forward_chain ? g.forward : g.reverse;
    for (const auto& m : all) {
      const bool chained = std::binary_search(c->anchors.begin(), c->anchors.end(), m, [](const Match& x, const Match& y) { return x.a < y.a; });
      if (m.reverse) out.push_back({m.a, m.b_end(), m.a_end(), m.b, true, chained});
      else out.push_back({m.a, m.b, m.a_end(), m.b_end(), false, chained});
    }
  }
  return out;
}

}  // namespace genomics::comparative
