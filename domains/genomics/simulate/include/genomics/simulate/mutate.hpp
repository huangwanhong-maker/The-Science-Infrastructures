// A diploid individual mutated from a reference, with the truth of every change: single
// nucleotide variants and short insertions and deletions, each heterozygous (on one
// haplotype, chosen at random) or homozygous (on both), never overlapping; the two haplotypes
// that result; and for each haplotype the map from its positions back to the reference's,
// which says where a simulated read truly comes from.
//
// Variants are written as VCF writes them (Danecek et al. 2011, Bioinformatics 27,
// 2156–2158): at a 0-based position of the reference, a REF and an ALT allele, an indel
// anchored on the base before it (a deletion of L bases after p has REF the L + 1 bases from
// p and ALT the base at p; an insertion has REF the base at p and ALT that base followed by
// the inserted ones). An indel in a repeat has several such descriptions; it is written where
// it was drawn, and comparison with calls is done after both are left-aligned (the variant
// module's normalisation).
//
// Draws: a stream keyed by the seed and MUTATIONS. Along the reference each position starts
// an SNV with probability `snv`, else an indel with probability `indel`, half insertions and
// half deletions, of a length geometric on 1, 2, ... with mean `indel_mean`, at most
// `indel_max`; then the zygosity; after a variant the next `spacing` bases stay unmutated.

#pragma once

#include <genomics/simulate/genome.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace genomics::simulate {

struct Variant {
  std::size_t position = 0;     // 0-based on the reference, of REF's first base
  std::string ref, alt;         // letters
  std::array<bool, 2> on{};     // which haplotypes carry it
  bool heterozygous() const { return on[0] != on[1]; }
  bool snv() const { return ref.size() == 1 && alt.size() == 1; }
  bool insertion() const { return alt.size() > ref.size(); }
  bool deletion() const { return ref.size() > alt.size(); }
};

struct MutationRates {
  double snv = 1e-3;
  double indel = 1e-4;
  double heterozygous = 0.6;  // the fraction of variants on one haplotype only
  double indel_mean = 2;
  std::size_t indel_max = 20;
  std::size_t spacing = 2;
};

// Maps positions of a haplotype to the reference: blocks [hap, hap + length) aligned to
// [ref, ref + length); inserted bases lie between blocks.
struct Liftover {
  struct Block { std::size_t hap = 0, ref = 0, length = 0; };
  std::vector<Block> blocks;
  std::size_t hap_length = 0, ref_length = 0;

  // The reference position of a haplotype position; an inserted base maps to the reference
  // base after the insertion.
  std::size_t to_reference(std::size_t h) const {
    auto it = std::upper_bound(blocks.begin(), blocks.end(), h, [](std::size_t x, const Block& b) { return x < b.hap; });
    if (it == blocks.begin()) return 0;
    const Block& b = *std::prev(it);
    return h < b.hap + b.length ? b.ref + (h - b.hap) : b.ref + b.length;
  }
  static Liftover identity(std::size_t n) {
    Liftover l;
    l.blocks.push_back({0, 0, n});
    l.hap_length = l.ref_length = n;
    return l;
  }
};

struct Diploid {
  std::array<Codes, 2> haplotype;
  std::array<Liftover, 2> liftover;
  std::vector<Variant> truth;  // in order of position
};

inline Diploid mutate_diploid(const Codes& reference, const MutationRates& rates, std::uint64_t seed, std::uint32_t contig = 0) {
  calx::random::Stream r(seed, MUTATIONS, contig, 0);
  Diploid d;
  const std::size_t n = reference.size();
  for (std::size_t p = 1; p + 1 < n; p++) {
    if (reference[p] > 3) continue;
    const double u = r.uniform();
    Variant v;
    v.position = p;
    if (u < rates.snv) {
      v.ref = std::string(1, sequence::dna_letter(reference[p]));
      v.alt = std::string(1, sequence::dna_letter(substitute(r, reference[p])));
    } else if (u < rates.snv + rates.indel) {
      const double q = 1 / std::max(1.0, rates.indel_mean);
      const std::size_t len = std::min<std::size_t>(rates.indel_max, 1 + r.geometric(q));
      if (r.bernoulli(0.5)) {
        if (p + len + 1 >= n) continue;
        bool clean = true;
        for (std::size_t i = p; i <= p + len; i++) clean = clean && reference[i] <= 3;
        if (!clean) continue;
        for (std::size_t i = p; i <= p + len; i++) v.ref += sequence::dna_letter(reference[i]);
        v.alt = v.ref.substr(0, 1);
      } else {
        v.ref = std::string(1, sequence::dna_letter(reference[p]));
        v.alt = v.ref;
        for (std::size_t i = 0; i < len; i++) v.alt += sequence::DNA_LETTERS[r.below(4)];
      }
    } else {
      continue;
    }
    if (r.uniform() < rates.heterozygous) {
      const bool which = r.bernoulli(0.5);
      v.on = {!which, which};
    } else {
      v.on = {true, true};
    }
    p += v.ref.size() - 1 + rates.spacing;
    d.truth.push_back(std::move(v));
  }
  // Apply each haplotype's variants, recording its blocks.
  for (int h = 0; h < 2; h++) {
    Codes& out = d.haplotype[h];
    Liftover& l = d.liftover[h];
    out.reserve(n + n / 100);
    std::size_t ref = 0;
    const auto copy_to = [&](std::size_t end) {
      if (end > ref) {
        l.blocks.push_back({out.size(), ref, end - ref});
        out.insert(out.end(), reference.begin() + static_cast<std::ptrdiff_t>(ref), reference.begin() + static_cast<std::ptrdiff_t>(end));
        ref = end;
      }
    };
    for (const auto& v : d.truth) {
      if (!v.on[h]) continue;
      if (v.snv()) {
        copy_to(v.position);
        l.blocks.push_back({out.size(), ref, 1});
        out.push_back(sequence::dna_code(v.alt[0]));
        ref++;
      } else if (v.insertion()) {
        copy_to(v.position + 1);
        for (std::size_t i = 1; i < v.alt.size(); i++) out.push_back(sequence::dna_code(v.alt[i]));
      } else {
        copy_to(v.position + 1);
        ref += v.ref.size() - 1;  // the deleted bases
      }
    }
    copy_to(n);
    l.hap_length = out.size();
    l.ref_length = n;
  }
  return d;
}

}  // namespace genomics::simulate
