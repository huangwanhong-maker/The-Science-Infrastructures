// Checks of whole-genome alignment by MUMs. The heaviest chain against the quadratic recurrence
// over the MUMs of random related pairs. A simulated genome of 200 kb against a copy of itself with
// a 30 kb inversion and then SNVs and indels (genomics/simulate, with their truth): every chained
// anchor on the truth (forward anchors where the copy's bases came from, reversed anchors where the
// inversion sent them), but for bases near its ends that match by chance or place an indel in a
// repeat elsewhere; the inversion recovered as the one reversed block, its breakpoints within
// 10 bases, and the flanks as two forward blocks; every isolated SNV a single-base gap at its
// place and every single-base gap an SNV; every isolated indel a gap of its length. The dot plot's
// segments, one a MUM, the chained ones marked. And, when the data cache has them, the
// mitochondria of human (NC_012920.1) and chimpanzee (NC_001643.1), measured on 2026-10-24: one
// forward chain over 0.964 of each (MUMs of 20 bases or more), no reversed anchor, and every forward
// MUM on one diagonal modulo the chimpanzee's length, as the two circles begin 582 bases apart; with
// the chimpanzee's circle rotated by that offset, the chain spans 0.999 of both.

#include <genomics/comparative/whole_genome.hpp>
#include <genomics/sequence/genbank.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/mutate.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::comparative;
using sequence::Codes;

namespace {

std::string fixed(double v, int digits = 3) {
  char b[32];
  std::snprintf(b, sizeof b, "%.*f", digits, v);
  return b;
}

std::string read_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::ostringstream s;
  s << f.rdbuf();
  return s.str();
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(1999, 0x57474131, 0, 0);

  // ---- the chain against the quadratic recurrence -----------------------------------------------------
  {
    bool heaviest = true;
    for (int trial = 0; trial < 40; trial++) {
      const Codes A = simulate::random_sequence(3000 + r.below(3000), 0.5, 7000 + static_cast<std::uint64_t>(trial));
      simulate::MutationRates rates;
      rates.snv = 0.02;
      rates.indel = 0.005;
      rates.heterozygous = 0;
      Codes B = simulate::mutate_diploid(A, rates, 8000 + static_cast<std::uint64_t>(trial)).haplotype[0];
      // A few segments moved, so that the MUMs are not all collinear.
      for (int j = 0; j < 3; j++) {
        const std::size_t len = 100 + r.below(300), from = r.below(B.size() - len), to = r.below(B.size() - len);
        Codes seg(B.begin() + static_cast<std::ptrdiff_t>(from), B.begin() + static_cast<std::ptrdiff_t>(from + len));
        B.erase(B.begin() + static_cast<std::ptrdiff_t>(from), B.begin() + static_cast<std::ptrdiff_t>(from + len));
        B.insert(B.begin() + static_cast<std::ptrdiff_t>(std::min(to, B.size())), seg.begin(), seg.end());
      }
      const auto mums = maximal_unique_matches(A, B, 12);
      const Chain c = chain_matches(mums, B.size());
      std::vector<std::int64_t> f(mums.size());
      std::int64_t best = 0;
      for (std::size_t i = 0; i < mums.size(); i++) {
        std::int64_t m = 0;
        for (std::size_t j = 0; j < i; j++)
          if (mums[j].a < mums[i].a && mums[j].b < mums[i].b) m = std::max(m, f[j]);
        f[i] = static_cast<std::int64_t>(mums[i].length) + m;
        best = std::max(best, f[i]);
      }
      bool increasing = true;
      for (std::size_t i = 1; i < c.anchors.size(); i++) increasing = increasing && c.anchors[i - 1].a < c.anchors[i].a && c.anchors[i - 1].b < c.anchors[i].b;
      heaviest = heaviest && increasing && c.weight == best && c.gaps.size() + 1 == c.anchors.size();
    }
    s.require("on 40 related pairs with moved segments, the chain increases in both genomes and weighs what the quadratic recurrence over its MUMs gives", heaviest);
  }

  // ---- an inversion and small variants, with their truth -------------------------------------------
  {
    const std::size_t n = 200'000, lo = 60'000, hi = 90'000, m = 20;
    const Codes A = simulate::random_sequence(n, 0.5, 1999);
    Codes R = A;  // A with [lo, hi) inverted
    const Codes inv = sequence::reverse_complement(Codes(A.begin() + lo, A.begin() + hi));
    std::copy(inv.begin(), inv.end(), R.begin() + lo);
    simulate::MutationRates rates;
    rates.snv = 1e-3;
    rates.indel = 2e-4;
    rates.heterozygous = 0;
    const simulate::Diploid d = simulate::mutate_diploid(R, rates, 2000);
    const Codes& B = d.haplotype[0];
    const simulate::Liftover& lift = d.liftover[0];
    // A's position of a base of R (the inversion undone).
    const auto to_a = [&](std::size_t q) { return q >= lo && q < hi ? lo + hi - 1 - q : q; };
    // Whether base h of B is a base of R (not inserted), and which.
    const auto block_of = [&](std::size_t h, std::size_t& q) {
      for (const auto& b : lift.blocks)
        if (h >= b.hap && h < b.hap + b.length) { q = b.ref + (h - b.hap); return true; }
      return false;
    };
    const GenomeAlignment g = align_genomes(A, B, {m, 1000, true});
    // Every base of every chained anchor where the truth puts it, but near the anchor's ends, where
    // a base inserted or across a breakpoint may match by chance, and where an indel in a repeat
    // has several placements (the anchor then ends where another placement puts it).
    std::size_t bases = 0, off = 0, off_inside = 0;
    for (const Chain* c : {&g.forward_chain, &g.reverse_chain})
      for (const auto& x : c->anchors)
        for (std::size_t t = 0; t < x.length; t++) {
          const std::size_t h = x.reverse ? x.b + x.length - 1 - t : x.b + t;
          std::size_t q = 0;
          bases++;
          if (block_of(h, q) && to_a(q) == x.a + t && ((q >= lo && q < hi) == x.reverse)) continue;
          off++;
          off_inside += std::min(t, x.length - 1 - t) >= 8;
        }
    s.require("200 kb with a 30 kb inversion, then " + std::to_string(d.truth.size()) + " SNVs and indels: every base of the chained anchors (" + std::to_string(g.forward_chain.anchors.size()) + " forward, " +
                  std::to_string(g.reverse_chain.anchors.size()) + " reversed) lies where the truth puts it, but for a few within 8 bases of an anchor's end",
              off_inside == 0, std::to_string(off) + " of " + std::to_string(bases) + " bases off the truth, " + std::to_string(off_inside) + " of them further inside");
    // The blocks: forward [0, lo), reversed [lo, hi), forward [hi, n), within a few bases.
    std::string blocks;
    for (const auto& b : g.blocks) blocks += std::string(blocks.empty() ? "" : ", ") + (b.reverse ? "−" : "+") + "[" + std::to_string(b.a_begin) + ", " + std::to_string(b.a_end) + ")";
    const auto near = [](std::size_t x, std::size_t y) { return (x > y ? x - y : y - x) <= 10; };  // measured 2026-10-24: within 2
    const bool three = g.blocks.size() == 3 && !g.blocks[0].reverse && g.blocks[1].reverse && !g.blocks[2].reverse;
    s.require("the alignment has three blocks: forward, reversed, forward", three, blocks);
    s.require("the reversed block is the inversion, its breakpoints in A within 10 bases of 60000 and 90000, and the forward blocks its flanks",
              three && near(g.blocks[1].a_begin, lo) && near(g.blocks[1].a_end, hi) && near(g.blocks[0].a_begin, 0) && near(g.blocks[0].a_end, lo) && near(g.blocks[2].a_begin, hi) && near(g.blocks[2].a_end, n),
              blocks);
    // The gaps against the variants, in A's coordinates.
    std::map<std::size_t, const simulate::Variant*> snvs;
    std::vector<const simulate::Variant*> isolated;
    for (std::size_t i = 0; i < d.truth.size(); i++) {
      const auto& v = d.truth[i];
      const std::size_t p = v.position;
      bool alone = p >= 2 * m && p + 2 * m < n && (p < lo || p >= lo + 2 * m) && (p + 2 * m < lo || p >= hi + 2 * m) && (p + 2 * m < hi || p >= hi);
      if (i > 0) alone = alone && d.truth[i - 1].position + d.truth[i - 1].ref.size() + 2 * m <= p;
      if (i + 1 < d.truth.size()) alone = alone && p + v.ref.size() + 2 * m <= d.truth[i + 1].position;
      if (alone) isolated.push_back(&v);
      if (v.snv()) snvs[to_a(p)] = &v;
    }
    std::size_t snv_found = 0, snv_isolated = 0, indel_found = 0, indel_isolated = 0, single = 0, single_true = 0;
    for (const auto* v : isolated) {
      const std::size_t pa = to_a(v->position);
      const Chain& c = v->position >= lo && v->position < hi ? g.reverse_chain : g.forward_chain;
      if (v->snv()) {
        snv_isolated++;
        for (const auto& gap : c.gaps) snv_found += gap.kind() == GapKind::snp && static_cast<std::size_t>(gap.a) == pa;
      } else {
        indel_isolated++;
        const std::int64_t len = static_cast<std::int64_t>(v->alt.size()) - static_cast<std::int64_t>(v->ref.size());
        bool found = false;
        for (const auto& gap : c.gaps) found = found || (gap.indel() == len && std::abs(gap.a - static_cast<std::int64_t>(pa)) <= static_cast<std::int64_t>(2 * m));
        indel_found += found;
      }
    }
    for (const Chain* c : {&g.forward_chain, &g.reverse_chain})
      for (const auto& gap : c->gaps)
        if (gap.kind() == GapKind::snp) {
          single++;
          single_true += snvs.count(static_cast<std::size_t>(gap.a)) > 0;
        }
    s.require("every isolated SNV (" + std::to_string(snv_isolated) + ") is a single-base gap at its place in A", snv_found == snv_isolated, std::to_string(snv_found) + " found");
    s.require("every single-base gap of the chains (" + std::to_string(single) + ") lies at a true SNV", single_true == single, std::to_string(single_true) + " true");
    s.require("every isolated indel (" + std::to_string(indel_isolated) + ") is a gap of its length near its place", indel_found == indel_isolated, std::to_string(indel_found) + " found");
    // The dot plot.
    const auto dots = dot_plot(g);
    std::size_t chained = 0;
    bool shapes = dots.size() == g.forward.size() + g.reverse.size();
    for (const auto& p : dots) {
      chained += p.chained;
      shapes = shapes && p.x1 > p.x0 && (p.reverse ? p.y0 > p.y1 : p.y1 > p.y0) && p.x1 - p.x0 == (p.reverse ? p.y0 - p.y1 : p.y1 - p.y0);
    }
    s.require("the dot plot has one segment a MUM (" + std::to_string(dots.size()) + "), forward ones rising and reversed ones falling at 45°, the chained ones marked",
              shapes && chained == g.forward_chain.anchors.size() + g.reverse_chain.anchors.size());
    s.require("the anchors cover nearly all of both genomes (at least 0.99)", g.coverage_a() >= 0.99 && g.coverage_b() >= 0.99, "A " + fixed(g.coverage_a()) + ", B " + fixed(g.coverage_b()));
  }

  // ---- human and chimpanzee mitochondria ---------------------------------------------------------------
  {
    const char* env = std::getenv("SCIENCE_CACHE");
    const std::filesystem::path dir = (env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE)) / "ncbi-genomes";
    if (!std::filesystem::exists(dir / "NC_012920.1.gb") || !std::filesystem::exists(dir / "NC_001643.1.gb")) {
      s.skip("the human and chimpanzee mitochondria (NC_012920.1, NC_001643.1) are not in the data cache: they are not aligned");
      return s.done();
    }
    const Codes human = sequence::encode_dna(sequence::read_genbank(read_file(dir / "NC_012920.1.gb")).at(0).sequence);
    const Codes chimp = sequence::encode_dna(sequence::read_genbank(read_file(dir / "NC_001643.1.gb")).at(0).sequence);
    // The two circles begin at different places: the chain of the linear sequences leaves out the
    // stretch across the origin, whose MUMs lie on the same diagonal modulo the length of B.
    // Rotating B by the chain's offset closes the circle.
    const auto measure = [&](const Codes& B, std::size_t m, const std::string& what, double need, std::int64_t& offset) {
      const GenomeAlignment g = align_genomes(human, B, {m, 1000, true});
      const Chain& c = g.forward_chain;
      const std::size_t span_a = c.anchors.empty() ? 0 : c.anchors.back().a_end() - c.anchors.front().a;
      const std::size_t span_b = c.anchors.empty() ? 0 : c.anchors.back().b_end() - c.anchors.front().b;
      std::vector<std::int64_t> offs;
      for (const auto& x : c.anchors) offs.push_back(static_cast<std::int64_t>(x.b) - static_cast<std::int64_t>(x.a));
      std::sort(offs.begin(), offs.end());
      offset = offs.empty() ? 0 : offs[offs.size() / 2];
      // Every forward MUM within 20 of the chain's diagonal, modulo |B|.
      const auto nB = static_cast<std::int64_t>(B.size());
      std::size_t on = 0;
      for (const auto& x : g.forward) {
        std::int64_t d = ((static_cast<std::int64_t>(x.b) - static_cast<std::int64_t>(x.a) - offset) % nB + nB) % nB;
        on += std::min(d, nB - d) <= 20;
      }
      std::size_t rev_bases = 0, fwd_bases = 0, snp = 0, subst = 0, indel = 0;
      for (const auto& x : g.reverse_chain.anchors) rev_bases += x.length;
      for (const auto& x : g.forward) fwd_bases += x.length;
      for (const auto& gap : c.gaps) (gap.kind() == GapKind::snp ? snp : gap.kind() == GapKind::substitutions ? subst : indel)++;
      const std::string detail = std::to_string(g.forward.size()) + " forward MUMs (" + std::to_string(fwd_bases) + " bases), " + std::to_string(on) + " on the diagonal modulo |B|, " + std::to_string(c.anchors.size()) +
                                 " chained over " + std::to_string(span_a) + " of " + std::to_string(human.size()) + " human and " + std::to_string(span_b) + " of " + std::to_string(B.size()) +
                                 " chimpanzee bases; coverage " + fixed(g.coverage_a()) + " and " + fixed(g.coverage_b()) + "; offsets b − a from " + std::to_string(offs.empty() ? 0 : offs.front()) + " to " +
                                 std::to_string(offs.empty() ? 0 : offs.back()) + "; gaps: " + std::to_string(snp) + " SNPs, " + std::to_string(subst) + " of substitutions, " + std::to_string(indel) + " indels; " +
                                 std::to_string(g.reverse_chain.anchors.size()) + " reversed anchors (" + std::to_string(rev_bases) + " bases); " + std::to_string(g.blocks.size()) + " blocks";
      s.require("human and chimpanzee mitochondria" + what + ", MUMs of " + std::to_string(m) + " or more: one forward chain spanning at least " + fixed(need, 2) + " of both genomes, every forward MUM on its diagonal modulo |B|, none reversed",
                span_a >= need * static_cast<double>(human.size()) && span_b >= need * static_cast<double>(B.size()) && on == g.forward.size() && g.reverse.empty(), detail);
    };
    std::int64_t offset = 0, again = 0;
    measure(chimp, 15, "", 0.95, again);
    measure(chimp, 20, "", 0.95, offset);
    // Rotate the chimpanzee's circle so that its base b + offset comes first... B'[i] = B[(i + offset) mod |B|].
    const auto nB = static_cast<std::int64_t>(chimp.size());
    const std::size_t shift = static_cast<std::size_t>(((offset % nB) + nB) % nB);
    Codes rotated(chimp.begin() + static_cast<std::ptrdiff_t>(shift), chimp.end());
    rotated.insert(rotated.end(), chimp.begin(), chimp.begin() + static_cast<std::ptrdiff_t>(shift));
    measure(rotated, 20, " (the chimpanzee's circle rotated by the chain's offset, " + std::to_string(offset) + ")", 0.99, again);
  }
  return s.done();
}
