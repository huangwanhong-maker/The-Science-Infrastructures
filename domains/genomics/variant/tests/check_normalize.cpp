// Checks of normalisation: a CA deleted from a CA repeat, written in ten ways, every
// representation of it normalising to one (the problem of Tan, Abecasis and Kang's Fig. 1);
// trimming of substitutions and padding at the first base; random indels in random tandem
// repeats, every representation of the haplotype
// normalising to one form, the parsimonious representation of least position (by brute force
// over all representations); and the left-alignment of the gaps of alignments, which keeps the
// score, is idempotent and leaves no gap that could move further left.

#include <genomics/align/gotoh.hpp>
#include <genomics/align/scoring.hpp>
#include <genomics/align/types.hpp>
#include <genomics/align/wavefront.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/variant/normalize.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <string>
#include <tuple>
#include <vector>

using namespace genomics;
using sequence::Codes;
using variant::Normalized;

namespace {

constexpr std::uint32_t CHECK = 0x4e4f524d;

struct Rep {
  std::size_t p;
  std::string ref, alt;
};

// Every (p, REF, ALT) with non-empty alleles, REF ≠ ALT and |REF| ≤ max_ref, that makes the haplotype h of the reference.
std::vector<Rep> representations(const Codes& reference, const Codes& h, std::size_t max_ref) {
  std::vector<Rep> out;
  const std::size_t n = reference.size();
  for (std::size_t p = 0; p < n; p++)
    for (std::size_t r = 1; r <= max_ref && p + r <= n; r++) {
      const std::size_t suffix = n - p - r;
      if (h.size() < p + suffix) continue;
      const std::size_t a = h.size() - p - suffix;
      if (a == 0) continue;
      if (!std::equal(reference.begin(), reference.begin() + static_cast<std::ptrdiff_t>(p), h.begin())) continue;
      if (!std::equal(reference.end() - static_cast<std::ptrdiff_t>(suffix), reference.end(), h.end() - static_cast<std::ptrdiff_t>(suffix))) continue;
      Rep x{p, sequence::decode_dna(Codes(reference.begin() + static_cast<std::ptrdiff_t>(p), reference.begin() + static_cast<std::ptrdiff_t>(p + r))),
            sequence::decode_dna(Codes(h.begin() + static_cast<std::ptrdiff_t>(p), h.begin() + static_cast<std::ptrdiff_t>(p + a)))};
      if (x.ref != x.alt) out.push_back(std::move(x));
    }
  return out;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, CHECK, 0, 0);

  // ---- a CA deleted from a CA repeat -----------------------------------------------------------------------
  {
    // The problem Tan, Abecasis and Kang (2015) open with (their Fig. 1 draws one indel in a
    // repeat written in several ways), on a CA deleted from the repeat CACACA: reference
    // GGGCACACAGGG, haplotype GGGCACAGGG. Written at four places, and with more context, it is
    // one variant; its normal form is G CA deleted after position 2: (2, GCA, G).
    const Codes ref = sequence::encode_dna("GGGCACACAGGG");
    const std::vector<Rep> written = {{2, "GCA", "G"},   {3, "CAC", "C"},   {4, "ACA", "A"},   {5, "CAC", "C"},     {6, "ACA", "A"},
                                      {3, "CACA", "CA"}, {2, "GCAC", "GC"}, {5, "CACAG", "CAG"}, {1, "GGCA", "GG"}, {0, "GGGCACACAGGG", "GGGCACAGGG"}};
    bool one = true, same = true;
    const Codes h = variant::apply(ref, 2, "GCA", "G");
    for (const auto& w : written) {
      const Normalized n = variant::normalize(w.p, w.ref, w.alt, ref);
      one = one && n.position == 2 && n.ref == "GCA" && n.alts.size() == 1 && n.alts[0] == "G" && n.ref_matches;
      same = same && variant::apply(ref, w.p, w.ref, w.alt) == h;
    }
    s.require("ten representations of a CA deleted from CACACA are one haplotype", same);
    s.require("every one normalises to (2, GCA, G): left-aligned, parsimonious", one);
    const Normalized ins = variant::normalize(8, "A", "ACA", ref);
    s.require("an inserted CA at the repeat's end normalises to its start, (2, G, GCA)", ins.position == 2 && ins.ref == "G" && ins.alts[0] == "GCA");
    const Normalized mnp = variant::normalize(3, "CAC", "CTC", ref);
    s.require("a substitution written with context is trimmed to the base, (4, A, T)", mnp.position == 4 && mnp.ref == "A" && mnp.alts[0] == "T");
    const Normalized multi = variant::normalize(3, "CACA", std::vector<std::string>{"CA", "CACACA"}, ref);
    s.require("a multi-allelic indel normalises its alleles together, (2, GCA, [G, GCACA])",
              multi.position == 2 && multi.ref == "GCA" && multi.alts == std::vector<std::string>{"G", "GCACA"}, multi.ref + " " + std::to_string(multi.position));
    const Codes aaa = sequence::encode_dna("AAACGT");
    const Normalized first = variant::normalize(1, "AA", "A", aaa);
    s.require("an A deleted from the leading AAA, where nothing precedes it, is padded with the base after as VCF pads position 1: (0, AA, A)",
              first.position == 0 && first.ref == "AA" && first.alts[0] == "A" && variant::apply(aaa, 0, first.ref, first.alts[0]) == variant::apply(aaa, 1, "AA", "A"),
              std::to_string(first.position) + " " + first.ref + " " + first.alts[0]);
    s.require("a REF that is not the reference's is reported", !variant::normalize(3, "GG", "G", ref).ref_matches);
  }

  // ---- random indels in random repeats ----------------------------------------------------------------------
  {
    bool unique = true, minimal = true, equivalent = true;
    std::size_t reps = 0;
    for (int trial = 0; trial < 400; trial++) {
      // A random flank, a tandem repeat of a unit of 1–4 bases, another flank.
      Codes ref;
      for (int i = 0; i < 6; i++) ref.push_back(static_cast<sequence::Code>(r.below(4)));
      Codes unit;
      for (std::uint64_t i = 0, u = 1 + r.below(4); i < u; i++) unit.push_back(static_cast<sequence::Code>(r.below(4)));
      const std::size_t copies = 2 + r.below(5);
      for (std::size_t c = 0; c < copies; c++) ref.insert(ref.end(), unit.begin(), unit.end());
      for (int i = 0; i < 6; i++) ref.push_back(static_cast<sequence::Code>(r.below(4)));
      // An indel inside the repeat: a deletion or insertion of k units, or of a random string.
      const std::size_t start = 6 + r.below(copies * unit.size());
      Codes h = ref;
      const bool deletion = r.bernoulli(0.5);
      const std::size_t len = r.bernoulli(0.7) ? unit.size() * (1 + r.below(2)) : 1 + r.below(3);
      if (deletion) {
        if (start + len > ref.size()) continue;
        h.erase(h.begin() + static_cast<std::ptrdiff_t>(start), h.begin() + static_cast<std::ptrdiff_t>(start + len));
      } else {
        Codes ins;
        for (std::size_t i = 0; i < len; i++) ins.push_back(r.bernoulli(0.7) ? unit[i % unit.size()] : static_cast<sequence::Code>(r.below(4)));
        h.insert(h.begin() + static_cast<std::ptrdiff_t>(start), ins.begin(), ins.end());
      }
      if (h == ref) continue;
      const auto all = representations(ref, h, 14);
      if (all.empty()) continue;
      reps += all.size();
      const Normalized n0 = variant::normalize(all[0].p, all[0].ref, all[0].alt, ref);
      std::size_t least = SIZE_MAX, least_p = SIZE_MAX;
      for (const auto& x : all) {
        const Normalized n = variant::normalize(x.p, x.ref, x.alt, ref);
        unique = unique && n.position == n0.position && n.ref == n0.ref && n.alts == n0.alts;
        const std::size_t size = x.ref.size() + x.alt.size();
        if (size < least || (size == least && x.p < least_p)) { least = size; least_p = x.p; }
      }
      equivalent = equivalent && variant::apply(ref, n0.position, n0.ref, n0.alts[0]) == h;
      minimal = minimal && n0.ref.size() + n0.alts[0].size() == least && n0.position == least_p;
    }
    s.require("every representation of a random indel in a random repeat normalises to one form (400 trials)", unique, std::to_string(reps) + " representations");
    s.require("the normal form makes the same haplotype", equivalent);
    s.require("the normal form is the shortest representation of least position (brute force)", minimal);
  }

  // ---- the gaps of alignments ---------------------------------------------------------------------------------
  {
    // A read with a CA deleted, its gap written at the repeat's right end.
    const Codes target = sequence::encode_dna("TTGGGCACACAGGGTT");
    const Codes query = sequence::encode_dna("TTGGGCACAGGGTT");
    align::Alignment a;
    a.query_end = query.size();
    a.target_end = target.size();
    a.cigar = align::parse_cigar("9M2D5M");
    const auto sc = align::minimap2_sr();
    const auto left = variant::left_align_gaps(a, query, target);
    align::Alignment b = a;
    b.cigar = left;
    s.require("a deletion written at a repeat's right end moves to its left end: 9M2D5M → 5M2D9M", align::cigar_string(left) == "5M2D9M", align::cigar_string(left));
    s.require("and keeps its score", align::rescore(a, query, target, sc) == align::rescore(b, query, target, sc));

    bool kept = true, idempotent = true, stuck = true;
    for (int trial = 0; trial < 300; trial++) {
      Codes t;
      for (int i = 0; i < 40; i++) t.push_back(static_cast<sequence::Code>(r.below(4)));
      // Repeats where gaps can move.
      for (int i = 0; i < 6; i++) t.insert(t.begin() + 20, {sequence::C, sequence::A});
      Codes q = t;
      for (int e = 0; e < 4; e++) {
        const std::size_t p = r.below(q.size());
        const auto kind = r.below(3);
        if (kind == 0) q[p] = static_cast<sequence::Code>(r.below(4));
        else if (kind == 1) q.erase(q.begin() + static_cast<std::ptrdiff_t>(p));
        else q.insert(q.begin() + static_cast<std::ptrdiff_t>(p), static_cast<sequence::Code>(r.below(4)));
      }
      const auto wf = align::wavefront_global(q, t, sc);
      align::Alignment x = wf;
      x.cigar = variant::left_align_gaps(wf, q, t);
      kept = kept && align::rescore(x, q, t, sc) >= align::rescore(wf, q, t, sc) && align::query_length(x.cigar) == q.size() && align::target_length(x.cigar) == t.size();
      idempotent = idempotent && variant::left_align_gaps(x, q, t) == x.cigar;
      // No gap preceded by a pair satisfies its condition.
      std::size_t i = 0, j = 0;
      for (std::size_t k = 0; k < x.cigar.size(); k++) {
        const auto& o = x.cigar[k];
        if ((o.op == 'D' || o.op == 'I') && k > 0 && x.cigar[k - 1].op == 'M') {
          const bool can = o.op == 'D' ? t[j - 1] == t[j + o.length - 1] : q[i - 1] == q[i + o.length - 1];
          stuck = stuck && !can;
        }
        if (align::consumes_query(o.op)) i += o.length;
        if (align::consumes_target(o.op)) j += o.length;
      }
    }
    s.require("left-aligned gaps never lower the score and keep the path's lengths (300 wavefront alignments)", kept);
    s.require("left-alignment is idempotent", idempotent);
    s.require("no left-aligned gap can move further left", stuck);
  }

  return s.done();
}
