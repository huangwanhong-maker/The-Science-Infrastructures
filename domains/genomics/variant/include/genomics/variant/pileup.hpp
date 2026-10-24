// The pileup of aligned reads over a region of a reference: at each position, the bases the
// reads put there with their qualities and strands; the insertions and deletions that start
// there; and the reference interval each read covers.
//
// - Which reads. Mappings on the region's contig that overlap it, of mapping quality at least
//   `min_mapping_quality` (20), primary or supplementary (secondary alignments are another
//   placement of a read counted once already); bases of quality below `min_base_quality` (13,
//   samtools' mpileup -Q) are left out, as bcftools mpileup leaves them (Li 2011,
//   Bioinformatics 27, 2987–2993, §2.1 assumes the qualities are calibrated).
// - Strands. A read mapped to the reverse strand is aligned as its reverse complement; its
//   bases and qualities are read on that strand, and every observation keeps the strand.
// - Indels. Each read's gaps are first moved as far left as they go without changing the
//   alignment's score (normalize.hpp's left_align_gaps), so that every read carrying one indel
//   in a repeat reports it at the same place. An indel is recorded at its anchor, the reference
//   base before it, as VCF anchors it: an insertion of bases S between reference positions p
//   and p + 1 at p (length +|S|, bases S); a deletion of reference bases [p + 1, p + L] at p
//   (length −L, bases the deleted ones). Gaps at a read's ends, not flanked by pairs on both
//   sides, are not recorded: such a gap says little and is usually an artefact of the end.
// - Spans. The reference interval of each read used, so that the reads spanning an indel's
//   context (those that could show it) can be counted.
//
// Layout: the bases of position p are the slice [offset[p − begin], offset[p − begin + 1]) of
// one flat array, filled by two passes (count, then place) with no allocation per base.
// Complexity O(Σ ℓ_r + n) for reads of aligned lengths ℓ_r over a region of n positions, plus the
// left-alignment's O(ℓ s).

#pragma once

#include <genomics/align/types.hpp>
#include <genomics/mapping/mapper.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/fasta.hpp>
#include <genomics/variant/normalize.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::variant {

struct PileupOptions {
  int min_base_quality = 13;
  int min_mapping_quality = 20;
  int default_quality = 20;   // the quality of a base of a read given without qualities
  bool left_align = true;
};

struct Observation {
  sequence::Code base = 0;
  std::uint8_t quality = 0, mapq = 0;
  bool reverse = false;
  std::uint32_t read = 0;
};

struct IndelObservation {
  std::size_t position = 0;  // the anchor, the reference base before the event
  std::int32_t length = 0;   // > 0 inserted, < 0 deleted bases
  std::string bases;         // the inserted bases, or the deleted reference bases, in letters
  bool reverse = false;
  std::uint8_t mapq = 0;
  std::uint32_t read = 0;
};

struct ReadSpan {
  std::size_t begin = 0, end = 0;
  bool reverse = false;
  std::uint32_t read = 0;
};

// A read as the pileup sees it: its bases and qualities as sequenced (Phred + 33, or empty).
struct ReadView {
  std::span<const sequence::Code> bases;
  std::string_view quality;
};

template <class R>
std::vector<ReadView> read_views(const std::vector<R>& reads) {
  std::vector<ReadView> v;
  v.reserve(reads.size());
  for (const auto& r : reads) v.push_back({std::span<const sequence::Code>(r.bases), std::string_view(r.quality)});
  return v;
}

struct Pileup {
  std::uint32_t contig = 0;
  std::size_t begin = 0, end = 0;
  std::vector<std::uint32_t> offset;     // end − begin + 1 entries
  std::vector<Observation> bases;
  std::vector<IndelObservation> indels;  // anchors in [begin, end), by (position, length, bases, read)
  std::vector<ReadSpan> spans;           // by begin
  std::size_t longest = 0;               // the longest span

  std::span<const Observation> at(std::size_t p) const {
    if (p < begin || p >= end) return {};
    return {bases.data() + offset[p - begin], offset[p - begin + 1] - offset[p - begin]};
  }
  std::size_t depth(std::size_t p) const { return at(p).size(); }
  // f(span) for each read whose aligned interval contains [b, e), in the order of the spans.
  template <class F>
  void for_each_spanning(std::size_t b, std::size_t e, F&& f) const {
    const std::size_t from = b > longest ? b - longest : 0;
    auto it = std::lower_bound(spans.begin(), spans.end(), from, [](const ReadSpan& s, std::size_t x) { return s.begin < x; });
    for (; it != spans.end() && it->begin <= b; ++it)
      if (it->end >= e) f(*it);
  }
  // The reads whose aligned interval contains [b, e), and how many of them are on each strand.
  std::size_t spanning(std::size_t b, std::size_t e, std::size_t* forward = nullptr, std::size_t* reverse = nullptr) const {
    std::size_t f = 0, r = 0;
    for_each_spanning(b, e, [&](const ReadSpan& s) { (s.reverse ? r : f)++; });
    if (forward) *forward = f;
    if (reverse) *reverse = r;
    return f + r;
  }
};

// The pileup of [begin, end) of a contig from mappings (each naming its read by Mapping::read in
// `reads`) of the read's strand to the reference.
inline Pileup pileup(const Codes& reference, std::uint32_t contig, std::size_t begin, std::size_t end, std::span<const mapping::Mapping> mappings,
                     std::span<const ReadView> reads, const PileupOptions& o = {}) {
  if (begin > end || end > reference.size()) throw std::invalid_argument("genomics::variant::pileup: a region of the reference");
  Pileup p;
  p.contig = contig;
  p.begin = begin;
  p.end = end;
  const std::size_t n = end - begin;
  // The mappings used, with their left-aligned paths.
  struct Used {
    const mapping::Mapping* m;
    align::Cigar cigar;
  };
  std::vector<Used> used;
  for (const auto& m : mappings) {
    if (m.contig != contig || !m.primary || m.mapq < o.min_mapping_quality || m.target_end <= begin || m.target_begin >= end) continue;
    if (m.read >= reads.size() || reads[m.read].bases.size() != m.query_length)
      throw std::invalid_argument("genomics::variant::pileup: a mapping whose read is not among the reads");
    Used u{&m, {}};
    if (o.left_align) {
      const auto& b = reads[m.read].bases;
      Codes q(b.begin(), b.end());
      if (m.reverse) q = sequence::reverse_complement(q);
      u.cigar = left_align_gaps(m.alignment, q, reference);
    } else {
      for (const auto& op : m.alignment.cigar) align::push(u.cigar, (op.op == '=' || op.op == 'X') ? 'M' : op.op, op.length);
    }
    used.push_back(std::move(u));
  }
  // A base of read r at position i of the strand mapped.
  const auto base_at = [&](const mapping::Mapping& m, std::size_t i) -> sequence::Code {
    const auto& b = reads[m.read].bases;
    return m.reverse ? sequence::complement(b[b.size() - 1 - i]) : b[i];
  };
  const auto quality_at = [&](const mapping::Mapping& m, std::size_t i) -> int {
    const auto& q = reads[m.read].quality;
    if (q.size() != m.query_length) return o.default_quality;
    return sequence::phred_of(m.reverse ? q[q.size() - 1 - i] : q[i]);
  };
  // Pass 1: counts; pass 2: placement.
  std::vector<std::uint32_t> count(n + 1, 0);
  for (int pass = 0; pass < 2; pass++) {
    if (pass == 1) {
      p.offset.assign(n + 1, 0);
      for (std::size_t x = 0; x < n; x++) p.offset[x + 1] = p.offset[x] + count[x];
      p.bases.resize(p.offset[n]);
      std::fill(count.begin(), count.end(), 0);
    }
    for (const auto& u : used) {
      const mapping::Mapping& m = *u.m;
      std::size_t i = m.alignment.query_begin, j = m.alignment.target_begin;
      for (std::size_t k = 0; k < u.cigar.size(); k++) {
        const auto& op = u.cigar[k];
        if (op.op == 'M') {
          for (std::uint32_t t = 0; t < op.length; t++, i++, j++) {
            if (j < begin || j >= end) continue;
            const int q = quality_at(m, i);
            if (q < o.min_base_quality) continue;
            const std::size_t x = j - begin;
            if (pass == 0) count[x]++;
            else p.bases[p.offset[x] + count[x]++] = {base_at(m, i), static_cast<std::uint8_t>(std::clamp(q, 0, 255)), static_cast<std::uint8_t>(std::clamp(m.mapq, 0, 255)), m.reverse, m.read};
          }
          continue;
        }
        const bool flanked = k > 0 && k + 1 < u.cigar.size() && u.cigar[k - 1].op == 'M' && u.cigar[k + 1].op == 'M';
        if (pass == 1 && flanked && j >= begin + 1 && j - 1 < end) {
          IndelObservation e;
          e.position = j - 1;
          e.reverse = m.reverse;
          e.mapq = static_cast<std::uint8_t>(std::clamp(m.mapq, 0, 255));
          e.read = m.read;
          if (op.op == 'I') {
            e.length = static_cast<std::int32_t>(op.length);
            for (std::uint32_t t = 0; t < op.length; t++) e.bases += sequence::dna_letter(base_at(m, i + t));
          } else if (op.op == 'D') {
            e.length = -static_cast<std::int32_t>(op.length);
            for (std::uint32_t t = 0; t < op.length; t++) e.bases += sequence::dna_letter(reference[j + t]);
          }
          if (e.length != 0) p.indels.push_back(std::move(e));
        }
        if (align::consumes_query(op.op)) i += op.length;
        if (align::consumes_target(op.op)) j += op.length;
      }
    }
  }
  std::sort(p.indels.begin(), p.indels.end(), [](const IndelObservation& a, const IndelObservation& b) {
    if (a.position != b.position) return a.position < b.position;
    if (a.length != b.length) return a.length < b.length;
    if (a.bases != b.bases) return a.bases < b.bases;
    return a.read < b.read;
  });
  p.spans.reserve(used.size());
  for (const auto& u : used) {
    p.spans.push_back({u.m->target_begin, u.m->target_end, u.m->reverse, u.m->read});
    p.longest = std::max(p.longest, u.m->target_end - u.m->target_begin);
  }
  std::sort(p.spans.begin(), p.spans.end(), [](const ReadSpan& a, const ReadSpan& b) { return a.begin != b.begin ? a.begin < b.begin : a.read < b.read; });
  return p;
}

}  // namespace genomics::variant
