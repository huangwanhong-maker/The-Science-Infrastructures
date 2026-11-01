// A diploid individual whose genome departs from its reference by structural variants, with the
// truth of every change: deletions, tandem and dispersed duplications, inversions, novel
// insertions, balanced reciprocal translocations between contigs, large gains and losses of
// copies, whole-contig aneuploidy and copy-neutral loss of heterozygosity, each heterozygous or
// homozygous, on top of the small variants of mutate.hpp.
//
// - The model. A rearranged chromosome is a word over oriented pieces of the haplotypes: each
//   piece an interval [begin, end) of the reference contig `contig`, taken from haplotype h of
//   that contig (mutate_diploid's, with its SNVs and short indels) forward or reversed and
//   complemented, or a novel sequence. The unrearranged individual is two chromosomes per contig,
//   each one piece; every structural variant rewrites the words, and the sequence and the
//   liftover of a chromosome are those of its word. This is the presentation of a rearranged
//   genome as a set of adjacencies of oriented breakends (theory/structural-variants.md, §1): two
//   consecutive pieces whose reference ends are not consecutive bases of one strand make a
//   junction, the pair of breakends that a read crossing it joins.
// - Breakends. A breakend (c, p, +) is the junction after base p of contig c, the reference to its
//   left kept (reads mapped forward upstream of it point to it); (c, p, −) the junction before
//   base p, the reference to its right kept. The rewrites, on reference intervals [a, b):
//     deletion            … [·, a)[b, ·) …                      junction (a − 1, +) ~ (b, −)
//     tandem duplication  … [·, b)[a, b)^k[b, ·) …               junction (b − 1, +) ~ (a, −), k times
//     inversion           … [·, a) rev[a, b) [b, ·) …            (a − 1, +) ~ (b − 1, +) and (a, −) ~ (b, −)
//     novel insertion     … [·, a) novel [a, ·) …                (a − 1, +) ~ (a, −) through the novel bases
//     dispersed duplic.   … [·, p) [a, b) or rev[a, b) [p, ·) …  two junctions at p, the copy's ends
//     reciprocal transl.  der(A) = A[0, pA) B[pB, |B|), der(B) = B[0, pB) A[pA, |A|)
//     gain, loss          a tandem duplication of k copies, a deletion, of a large interval
//     aneuploidy          a chromosome of the contig copied (+k) or removed (−1)
//     copy-neutral LOH    … [·, a) [a, b) of the other haplotype [b, ·) …   no junction at all
//   A heterozygous variant rewrites one haplotype (`on`), a homozygous one both; a translocation
//   joins haplotype h of contig A with haplotype h of contig B for each h it is on.
// - The truth: the variants as given (with the bases of novel insertions drawn), each
//   chromosome's pieces, sequence, liftover (simulate::Liftover with oriented blocks of any
//   contig) and junctions, and the copy number of every reference base as a step function by
//   haplotype of origin (`CopySegment`): from it the expected depth of reads and the expected
//   fraction of the B allele at each germline heterozygous site follow.
//
// Random variants (`random_structural`) are drawn from a stream keyed by the seed and STRUCTURAL,
// placed apart from one another, from the contigs' ends and from runs of N by `margin` bases; novel
// bases from a stream keyed by NOVEL and the variant's index. Every variant is checked on entry:
// intervals inside their contig and clear of N, pairwise disjoint (with the insertion sites of
// dispersed copies and the breakpoints of translocations) on each contig, at most one
// translocation per chromosome and no aneuploidy of a translocated contig; std::invalid_argument
// otherwise.
//
// Complexity: O(n + v log v) for a genome of n bases and v variants. Limits: variants never
// nest (a duplication inside an inversion, a chain of translocations); complex rearrangements
// (chromothripsis, templated insertions) are words of the same kind but are not drawn here.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/mutate.hpp>

#include <calx/random/stream.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace genomics::simulate {

inline constexpr std::uint32_t STRUCTURAL = 0x53545256, NOVEL = 0x4e4f564c;  // "STRV", "NOVL"

enum class SvType : std::uint8_t {
  deletion,
  tandem_duplication,
  dispersed_duplication,
  inversion,
  insertion,
  translocation,
  aneuploidy,
  loss_of_heterozygosity,
};

inline const char* sv_type_name(SvType t) {
  switch (t) {
    case SvType::deletion: return "DEL";
    case SvType::tandem_duplication: return "DUP:TANDEM";
    case SvType::dispersed_duplication: return "DUP:DISPERSED";
    case SvType::inversion: return "INV";
    case SvType::insertion: return "INS";
    case SvType::translocation: return "TRA";
    case SvType::aneuploidy: return "ANEUPLOIDY";
    case SvType::loss_of_heterozygosity: return "LOH";
  }
  return "?";
}

struct StructuralVariant {
  SvType type = SvType::deletion;
  std::uint32_t contig = 0;
  std::size_t begin = 0, end = 0;     // the reference interval rewritten; an insertion's site is begin (= end)
  std::uint32_t partner = 0;          // a dispersed copy's destination contig; a translocation's other contig
  std::size_t partner_position = 0;   // the copy's insertion site; the other contig's breakpoint
  bool inverted = false;              // a dispersed copy inserted reverse-complemented
  int copies = 1;                     // extra tandem copies; for aneuploidy the change in copies (+k or −1)
  std::size_t insert_length = 0;      // the length of a novel insertion to draw, when `sequence` is empty
  std::string sequence;               // a novel insertion's bases
  std::array<bool, 2> on{true, false};  // the haplotypes carrying it; for LOH the haplotype overwritten
  bool heterozygous() const { return on[0] != on[1]; }
  // The bases it changes: deleted, duplicated, inverted, copied, inserted or made homozygous.
  std::size_t length() const {
    switch (type) {
      case SvType::insertion: return sequence.empty() ? insert_length : sequence.size();
      case SvType::translocation: return 0;
      default: return end - begin;
    }
  }
};

// A piece of a chromosome: reference bases [begin, end) of `contig`, from its haplotype
// `haplotype`, forward or reversed; or (novel) bases [begin, end) of the sample's novel sequence
// `contig`.
struct Piece {
  bool novel = false;
  std::uint32_t contig = 0;
  int haplotype = 0;
  std::size_t begin = 0, end = 0;
  bool reverse = false;
};

// The adjacency at a junction of a chromosome: the breakend before it and the one after, with the
// novel bases between them (an insertion's).
struct Junction {
  std::uint32_t chromosome = 0;
  std::size_t at = 0;                    // the chromosome position of the first base after the junction (after the novel bases)
  std::uint32_t left_contig = 0, right_contig = 0;
  std::size_t left = 0, right = 0;       // the breakends' reference positions
  bool left_plus = true, right_plus = false;  // their orientations (+: the reference to the left of the position kept)
  std::string inserted;                  // novel bases between the two breakends
  std::uint32_t variant = 0;             // the truth variant it belongs to
};

struct Chromosome {
  std::string name;
  std::uint32_t contig = 0;  // the contig of its first piece (its centromere's, as der(A) is A's)
  int haplotype = 0;
  std::vector<Piece> pieces;
  Codes sequence;
  Liftover map;
  std::vector<std::size_t> junctions;  // positions where a junction lies (Junction::at, and before novel bases)
};

// The copies of reference bases [begin, end) of a contig in the sample, by the haplotype they come from.
struct CopySegment {
  std::uint32_t contig = 0;
  std::size_t begin = 0, end = 0;
  std::array<int, 2> copies{1, 1};
  int total() const { return copies[0] + copies[1]; }
};

struct StructuralSample {
  std::vector<Diploid> diploids;           // the small variants, one per contig, as mutate_diploid made them
  std::vector<StructuralVariant> truth;    // in the order given, novel bases drawn
  std::vector<Chromosome> chromosomes;
  std::vector<Junction> junctions;
  std::vector<CopySegment> copy_number;    // by contig and position, covering every contig
  std::vector<std::string> novel;          // novel sequences, by piece index
  std::size_t length() const {
    std::size_t n = 0;
    for (const auto& c : chromosomes) n += c.sequence.size();
    return n;
  }
  // The copies of reference base p of a contig, by haplotype.
  std::array<int, 2> copies_at(std::uint32_t contig, std::size_t p) const {
    auto it = std::upper_bound(copy_number.begin(), copy_number.end(), std::pair{contig, p}, [](const std::pair<std::uint32_t, std::size_t>& x, const CopySegment& s) {
      return x.first != s.contig ? x.first < s.contig : x.second < s.begin;
    });
    if (it == copy_number.begin()) return {0, 0};
    --it;
    return it->contig == contig && p < it->end ? it->copies : std::array<int, 2>{0, 0};
  }
};

namespace detail {

// The haplotype position of reference position x (a deleted base maps to the haplotype base after the deletion).
inline std::size_t to_haplotype(const Liftover& l, std::size_t x) {
  if (x >= l.ref_length) return l.hap_length;
  auto it = std::upper_bound(l.blocks.begin(), l.blocks.end(), x, [](std::size_t v, const Liftover::Block& b) { return v < b.ref; });
  if (it == l.blocks.begin()) return l.blocks.empty() ? x : l.blocks.front().hap;
  const auto& b = *std::prev(it);
  if (x < b.ref + b.length) return b.hap + (x - b.ref);
  return it == l.blocks.end() ? l.hap_length : it->hap;
}

inline bool clean(const Codes& c, std::size_t a, std::size_t b) {
  for (std::size_t i = a; i < b; i++)
    if (c[i] > 3) return false;
  return true;
}

// The breakend leaving a piece at its right end, and the one entering a piece at its left end.
inline void exit_breakend(const Piece& p, std::uint32_t& contig, std::size_t& position, bool& plus) {
  contig = p.contig;
  if (!p.reverse) { position = p.end - 1; plus = true; }
  else { position = p.begin; plus = false; }
}
inline void entry_breakend(const Piece& p, std::uint32_t& contig, std::size_t& position, bool& plus) {
  contig = p.contig;
  if (!p.reverse) { position = p.begin; plus = false; }
  else { position = p.end - 1; plus = true; }
}

inline bool contiguous(const Piece& a, const Piece& b) {
  if (a.novel || b.novel || a.contig != b.contig || a.reverse != b.reverse) return false;
  return a.reverse ? b.end == a.begin : a.end == b.begin;
}

}  // namespace detail

// Throws std::invalid_argument unless the variants can be applied together (see the header).
inline void validate_structural(const std::vector<Codes>& contigs, const std::vector<StructuralVariant>& svs) {
  const auto fail = [](std::size_t i, const std::string& m) { throw std::invalid_argument("genomics::simulate structural variant " + std::to_string(i) + ": " + m); };
  // Taken intervals per contig, to keep the variants apart.
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> taken(contigs.size());
  std::vector<int> translocated(contigs.size() * 2, 0);
  std::vector<bool> aneuploid(contigs.size(), false);
  const auto take = [&](std::size_t i, std::uint32_t c, std::size_t a, std::size_t b) {
    for (const auto& [x, y] : taken[c])
      if (a < y + 1 && x < b + 1) fail(i, "overlaps or touches another variant");
    taken[c].push_back({a, b});
  };
  for (std::size_t i = 0; i < svs.size(); i++) {
    const auto& v = svs[i];
    if (v.contig >= contigs.size()) fail(i, "no such contig");
    const Codes& c = contigs[v.contig];
    if (!v.on[0] && !v.on[1]) fail(i, "on no haplotype");
    switch (v.type) {
      case SvType::deletion:
      case SvType::tandem_duplication:
      case SvType::inversion:
      case SvType::loss_of_heterozygosity:
        if (!(v.begin > 0 && v.begin < v.end && v.end < c.size())) fail(i, "an interval strictly inside its contig");
        if (!detail::clean(c, v.begin - 1, v.end + 1)) fail(i, "an interval clear of N");
        if (v.type == SvType::tandem_duplication && v.copies < 1) fail(i, "at least one extra copy");
        if (v.type == SvType::loss_of_heterozygosity && !v.heterozygous()) fail(i, "a loss of heterozygosity overwrites one haplotype");
        take(i, v.contig, v.begin, v.end);
        break;
      case SvType::insertion:
        if (!(v.begin > 0 && v.begin < c.size()) || v.end != v.begin) fail(i, "an insertion site strictly inside its contig, begin = end");
        if (v.sequence.empty() && v.insert_length == 0) fail(i, "an insertion of at least one base");
        if (!detail::clean(c, v.begin - 1, v.begin + 1)) fail(i, "a site clear of N");
        take(i, v.contig, v.begin, v.begin);
        break;
      case SvType::dispersed_duplication: {
        if (!(v.begin > 0 && v.begin < v.end && v.end < c.size())) fail(i, "a copied interval strictly inside its contig");
        if (v.partner >= contigs.size()) fail(i, "no destination contig");
        const Codes& d = contigs[v.partner];
        if (!(v.partner_position > 0 && v.partner_position < d.size())) fail(i, "a destination strictly inside its contig");
        if (!detail::clean(c, v.begin - 1, v.end + 1) || !detail::clean(d, v.partner_position - 1, v.partner_position + 1)) fail(i, "clear of N");
        take(i, v.contig, v.begin, v.end);
        take(i, v.partner, v.partner_position, v.partner_position);
        break;
      }
      case SvType::translocation: {
        if (v.partner >= contigs.size() || v.partner == v.contig) fail(i, "a translocation joins two different contigs");
        const Codes& d = contigs[v.partner];
        if (!(v.begin > 0 && v.begin < c.size()) || !(v.partner_position > 0 && v.partner_position < d.size())) fail(i, "breakpoints strictly inside their contigs");
        if (!detail::clean(c, v.begin - 1, v.begin + 1) || !detail::clean(d, v.partner_position - 1, v.partner_position + 1)) fail(i, "breakpoints clear of N");
        take(i, v.contig, v.begin, v.begin);
        take(i, v.partner, v.partner_position, v.partner_position);
        for (int h = 0; h < 2; h++)
          if (v.on[h]) {
            if (translocated[v.contig * 2 + h]++ || translocated[v.partner * 2 + h]++) fail(i, "a chromosome in two translocations");
          }
        break;
      }
      case SvType::aneuploidy:
        if (v.copies == 0 || v.copies < -1) fail(i, "a change of +k or −1 copies");
        if (v.heterozygous() == false) fail(i, "aneuploidy names the one haplotype gained or lost");
        if (aneuploid[v.contig]) fail(i, "one aneuploidy per contig");
        aneuploid[v.contig] = true;
        break;
    }
  }
  for (const auto& v : svs)
    if (v.type == SvType::aneuploidy && (translocated[v.contig * 2] || translocated[v.contig * 2 + 1]))
      throw std::invalid_argument("genomics::simulate structural variants: an aneuploid contig is not also translocated");
}

// The individual: small variants by mutate_diploid (seed, contig), then the structural variants.
inline StructuralSample make_structural_sample(const std::vector<Codes>& contigs, std::vector<StructuralVariant> svs, const MutationRates& rates, std::uint64_t seed,
                                               const std::vector<std::string>& names = {}) {
  validate_structural(contigs, svs);
  StructuralSample s;
  for (std::uint32_t c = 0; c < contigs.size(); c++) s.diploids.push_back(mutate_diploid(contigs[c], rates, seed, c));
  // Novel bases.
  for (std::size_t i = 0; i < svs.size(); i++) {
    auto& v = svs[i];
    if (v.type != SvType::insertion) continue;
    if (v.sequence.empty()) {
      calx::random::Stream r(seed, NOVEL, static_cast<std::uint32_t>(i), 0);
      for (std::size_t k = 0; k < v.insert_length; k++) v.sequence += sequence::DNA_LETTERS[r.below(4)];
    }
    v.insert_length = v.sequence.size();
  }
  // The events of each (contig, haplotype), in reference order: what the sweep emits at a position.
  struct Event {
    std::size_t at = 0;    // the reference position where it starts
    std::size_t index = 0; // the variant
    bool destination = false;  // a dispersed copy's insertion site, or a translocation's partner breakpoint
  };
  const std::size_t C = contigs.size();
  std::vector<std::vector<Event>> events(C * 2);
  for (std::size_t i = 0; i < svs.size(); i++) {
    const auto& v = svs[i];
    for (int h = 0; h < 2; h++) {
      if (!v.on[h]) continue;
      switch (v.type) {
        case SvType::aneuploidy: break;
        case SvType::dispersed_duplication:
        case SvType::translocation:
          if (v.type == SvType::translocation) events[v.contig * 2 + h].push_back({v.begin, i, false});
          events[v.partner * 2 + h].push_back({v.partner_position, i, true});
          break;
        default: events[v.contig * 2 + h].push_back({v.begin, i, false});
      }
    }
  }
  // Word of pieces of each (contig, haplotype), and the cut at a translocation's breakpoint.
  std::vector<std::vector<Piece>> words(C * 2);
  std::vector<std::size_t> cut(C * 2, SIZE_MAX);   // the number of pieces before the breakpoint
  std::vector<std::size_t> cut_variant(C * 2, SIZE_MAX);
  for (std::uint32_t c = 0; c < C; c++)
    for (int h = 0; h < 2; h++) {
      auto& ev = events[c * 2 + h];
      std::sort(ev.begin(), ev.end(), [](const Event& a, const Event& b) { return a.at < b.at; });
      auto& w = words[c * 2 + h];
      std::size_t p = 0;
      const auto emit = [&](std::size_t a, std::size_t b) {
        if (b > a) w.push_back({false, c, h, a, b, false});
      };
      for (const auto& e : ev) {
        const auto& v = svs[e.index];
        if (e.destination) {
          emit(p, e.at);
          p = e.at;
          if (v.type == SvType::dispersed_duplication) w.push_back({false, v.contig, h, v.begin, v.end, v.inverted});
          else { cut[c * 2 + h] = w.size(); cut_variant[c * 2 + h] = e.index; }
          continue;
        }
        switch (v.type) {
          case SvType::deletion: emit(p, v.begin); p = v.end; break;
          case SvType::tandem_duplication:
            emit(p, v.end);
            for (int k = 0; k < v.copies; k++) w.push_back({false, c, h, v.begin, v.end, false});
            p = v.end;
            break;
          case SvType::inversion:
            emit(p, v.begin);
            w.push_back({false, c, h, v.begin, v.end, true});
            p = v.end;
            break;
          case SvType::insertion:
            emit(p, v.begin);
            w.push_back({true, static_cast<std::uint32_t>(s.novel.size()), h, 0, v.sequence.size(), false});
            s.novel.push_back(v.sequence);
            p = v.begin;
            break;
          case SvType::loss_of_heterozygosity:
            emit(p, v.begin);
            w.push_back({false, c, 1 - h, v.begin, v.end, false});
            p = v.end;
            break;
          case SvType::translocation:
            emit(p, v.begin);
            p = v.begin;
            cut[c * 2 + h] = w.size();
            cut_variant[c * 2 + h] = e.index;
            break;
          default: break;
        }
      }
      emit(p, contigs[c].size());
    }
  // The chromosomes: translocations exchange the arms of two words.
  struct Word { std::uint32_t contig; int haplotype; std::vector<Piece> pieces; };
  std::vector<Word> chromosomes;
  for (std::uint32_t c = 0; c < C; c++)
    for (int h = 0; h < 2; h++) {
      const std::size_t k = c * 2 + h;
      if (cut[k] == SIZE_MAX) { chromosomes.push_back({c, h, words[k]}); continue; }
      const auto& v = svs[cut_variant[k]];
      const std::uint32_t other = v.contig == c ? v.partner : v.contig;
      const std::size_t j = other * 2 + h;
      const auto& a = words[k];
      const auto& b = words[j];
      Word der{c, h, {}};
      der.pieces.assign(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(cut[k]));
      der.pieces.insert(der.pieces.end(), b.begin() + static_cast<std::ptrdiff_t>(cut[j]), b.end());
      chromosomes.push_back(std::move(der));
    }
  // Aneuploidy: copies of the chromosome of (contig, haplotype) added, or it removed.
  for (const auto& v : svs) {
    if (v.type != SvType::aneuploidy) continue;
    const int h = v.on[0] ? 0 : 1;
    const auto it = std::find_if(chromosomes.begin(), chromosomes.end(), [&](const Word& w) { return w.contig == v.contig && w.haplotype == h; });
    const std::ptrdiff_t at = it - chromosomes.begin();
    if (v.copies < 0) chromosomes.erase(it);
    else {
      const Word w = *it;
      for (int k = 0; k < v.copies; k++) chromosomes.insert(chromosomes.begin() + at + 1, w);
    }
  }
  // Sequences, liftovers and junctions.
  const auto variant_of_junction = [&](const Piece& a, const Piece& b) -> std::uint32_t {
    // The variant whose rewrite made the adjacency a·b: the nearest by position on either piece's contig.
    std::uint32_t best = 0;
    std::size_t best_d = SIZE_MAX;
    std::uint32_t ca, cb;
    std::size_t pa, pb;
    bool oa, ob;
    detail::exit_breakend(a, ca, pa, oa);
    detail::entry_breakend(b, cb, pb, ob);
    for (std::size_t i = 0; i < svs.size(); i++) {
      const auto& v = svs[i];
      if (v.type == SvType::aneuploidy || v.type == SvType::loss_of_heterozygosity) continue;
      const auto dist = [&](std::uint32_t c, std::size_t p) {
        std::size_t d = SIZE_MAX;
        const auto near = [&](std::uint32_t vc, std::size_t lo, std::size_t hi) {
          if (vc != c) return;
          const std::size_t x = p + 1 < lo ? lo - p - 1 : p > hi ? p - hi : 0;
          d = std::min(d, x);
        };
        near(v.contig, v.begin, v.end);
        if (v.type == SvType::dispersed_duplication || v.type == SvType::translocation) near(v.partner, v.partner_position, v.partner_position);
        return d;
      };
      const std::size_t d = std::min(dist(a.novel ? cb : ca, a.novel ? pb : pa), dist(b.novel ? ca : cb, b.novel ? pa : pb));
      if (d < best_d) { best_d = d; best = static_cast<std::uint32_t>(i); }
    }
    return best;
  };
  for (std::size_t w = 0; w < chromosomes.size(); w++) {
    const Word& word = chromosomes[w];
    Chromosome out;
    out.contig = word.contig;
    out.haplotype = word.haplotype;
    out.pieces = word.pieces;
    out.name = (names.size() > word.contig ? names[word.contig] : "contig" + std::to_string(word.contig)) + "#" + std::to_string(word.haplotype + 1);
    if (word.pieces.size() > 1 && std::any_of(word.pieces.begin(), word.pieces.end(), [&](const Piece& p) { return !p.novel && p.contig != word.contig; })) out.name = "der(" + out.name + ")";
    for (std::size_t k = 0; k < word.pieces.size(); k++) {
      const Piece& p = word.pieces[k];
      const std::size_t offset = out.sequence.size();
      if (p.novel) {
        out.sequence.reserve(out.sequence.size() + p.end - p.begin);
        for (std::size_t i = p.begin; i < p.end; i++) out.sequence.push_back(sequence::dna_code(s.novel[p.contig][i]));
      } else {
        const Diploid& d = s.diploids[p.contig];
        const Codes& hap = d.haplotype[p.haplotype];
        const Liftover& lift = d.liftover[p.haplotype];
        const std::size_t ha = detail::to_haplotype(lift, p.begin), hb = detail::to_haplotype(lift, p.end);
        Codes part(hap.begin() + static_cast<std::ptrdiff_t>(ha), hap.begin() + static_cast<std::ptrdiff_t>(hb));
        if (p.reverse) part = sequence::reverse_complement(part);
        out.sequence.insert(out.sequence.end(), part.begin(), part.end());
        std::vector<Liftover::Block> blocks;
        for (const auto& b : lift.blocks) {
          const std::size_t y0 = std::max(b.hap, ha), y1 = std::min(b.hap + b.length, hb);
          if (y0 >= y1) continue;
          const std::size_t ref = b.ref + (y0 - b.hap);
          if (!p.reverse) blocks.push_back({offset + (y0 - ha), ref, y1 - y0, p.contig, false});
          else blocks.push_back({offset + (hb - y1), ref, y1 - y0, p.contig, true});
        }
        if (p.reverse) std::reverse(blocks.begin(), blocks.end());
        out.map.blocks.insert(out.map.blocks.end(), blocks.begin(), blocks.end());
      }
      // A junction before this piece: from the last reference piece before it.
      if (k > 0 && !p.novel) {
        std::size_t prev = k - 1;
        std::string inserted;
        while (prev > 0 && word.pieces[prev].novel) prev--;
        for (std::size_t q = prev + 1; q < k; q++) inserted += s.novel[word.pieces[q].contig].substr(word.pieces[q].begin, word.pieces[q].end - word.pieces[q].begin);
        const Piece& a = word.pieces[prev];
        if (!a.novel && (!inserted.empty() || !detail::contiguous(a, p))) {
          Junction j;
          j.chromosome = static_cast<std::uint32_t>(w);
          j.at = offset;
          detail::exit_breakend(a, j.left_contig, j.left, j.left_plus);
          detail::entry_breakend(p, j.right_contig, j.right, j.right_plus);
          j.inserted = inserted;
          j.variant = variant_of_junction(a, p);
          s.junctions.push_back(j);
          if (!inserted.empty()) out.junctions.push_back(offset - inserted.size());
          out.junctions.push_back(offset);
        }
      }
    }
    out.map.hap_length = out.sequence.size();
    out.map.ref_length = 0;
    s.chromosomes.push_back(std::move(out));
  }
  // Copy number by haplotype of origin: a sweep over every piece.
  for (std::uint32_t c = 0; c < C; c++) {
    std::vector<std::pair<std::size_t, std::array<int, 2>>> delta;  // position, change
    for (const auto& ch : s.chromosomes)
      for (const auto& p : ch.pieces) {
        if (p.novel || p.contig != c) continue;
        std::array<int, 2> up{0, 0}, down{0, 0};
        up[p.haplotype] = 1;
        down[p.haplotype] = -1;
        delta.push_back({p.begin, up});
        delta.push_back({p.end, down});
      }
    std::sort(delta.begin(), delta.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    std::array<int, 2> now{0, 0};
    std::size_t at = 0, i = 0;
    const std::size_t n = contigs[c].size();
    while (at < n) {
      while (i < delta.size() && delta[i].first == at) { now[0] += delta[i].second[0]; now[1] += delta[i].second[1]; i++; }
      const std::size_t next = i < delta.size() ? std::min(n, delta[i].first) : n;
      if (!s.copy_number.empty() && s.copy_number.back().contig == c && s.copy_number.back().copies == now) s.copy_number.back().end = next;
      else s.copy_number.push_back({c, at, next, now});
      at = next;
    }
    for (auto& ch : s.chromosomes)
      if (ch.contig == c) ch.map.ref_length = n;
  }
  s.truth = std::move(svs);
  return s;
}

// ---- random structural variants ---------------------------------------------------------------
struct SizeRange {
  std::size_t low = 0, high = 0;  // drawn log-uniformly in [low, high]
};

struct StructuralSpec {
  std::size_t deletions = 0, tandem_duplications = 0, dispersed_duplications = 0, inversions = 0, insertions = 0;
  std::size_t translocations = 0, gains = 0, losses = 0, losses_of_heterozygosity = 0;
  SizeRange deletion{50, 20000}, tandem_duplication{300, 20000}, dispersed_duplication{500, 5000}, inversion{1000, 20000}, insertion{50, 200};
  SizeRange gain{100000, 400000}, loss{100000, 400000}, loss_of_heterozygosity{200000, 800000};
  double homozygous = 0.3;            // the fraction of variants on both haplotypes (never for LOH)
  std::size_t margin = 2000;          // bases kept free around every variant, contig end and run of N
  std::vector<int> aneuploidy;        // per contig: 0, or the change of copies (+1 trisomy, −1 monosomy)
};

// Variants drawn at random on the contigs, apart from one another and from N (see the header).
inline std::vector<StructuralVariant> random_structural(const std::vector<Codes>& contigs, const StructuralSpec& spec, std::uint64_t seed) {
  calx::random::Stream r(seed, STRUCTURAL, 0, 0);
  std::vector<StructuralVariant> out;
  // Forbidden intervals: N runs, ends, and every variant's span, each widened by the margin.
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> taken(contigs.size());
  std::size_t total = 0;
  for (std::uint32_t c = 0; c < contigs.size(); c++) {
    const Codes& x = contigs[c];
    total += x.size();
    taken[c].push_back({0, std::min(x.size(), spec.margin)});
    taken[c].push_back({x.size() > spec.margin ? x.size() - spec.margin : 0, x.size()});
    for (std::size_t i = 0; i < x.size();) {
      if (x[i] <= 3) { i++; continue; }
      std::size_t j = i;
      while (j < x.size() && x[j] > 3) j++;
      taken[c].push_back({i > spec.margin ? i - spec.margin : 0, std::min(x.size(), j + spec.margin)});
      i = j;
    }
  }
  if (total == 0) return out;
  const auto free = [&](std::uint32_t c, std::size_t a, std::size_t b) {
    for (const auto& [x, y] : taken[c])
      if (a < y && x < b) return false;
    return true;
  };
  const auto hold = [&](std::uint32_t c, std::size_t a, std::size_t b) {
    taken[c].push_back({a > spec.margin ? a - spec.margin : 0, std::min(contigs[c].size(), b + spec.margin)});
  };
  const auto size_of = [&](const SizeRange& s) {
    if (s.high <= s.low) return s.low;
    const double u = r.uniform();
    return static_cast<std::size_t>(std::llround(std::exp(std::log(static_cast<double>(s.low)) + u * (std::log(static_cast<double>(s.high)) - std::log(static_cast<double>(s.low))))));
  };
  // A place for an interval of `len` bases: a contig in proportion to its length, then a start.
  const auto place = [&](std::size_t len, std::uint32_t& contig, std::size_t& begin, std::uint32_t avoid = UINT32_MAX) {
    for (int attempt = 0; attempt < 2000; attempt++) {
      std::size_t u = r.below(total);
      std::uint32_t c = 0;
      while (u >= contigs[c].size()) { u -= contigs[c].size(); c++; }
      if (c == avoid || contigs[c].size() < len + 2) continue;
      const std::size_t b = r.below(contigs[c].size() - len);
      if (!free(c, b, b + len + 1)) continue;
      contig = c;
      begin = b;
      return true;
    }
    return false;
  };
  const auto zygosity = [&](StructuralVariant& v) {
    if (r.uniform() < spec.homozygous) v.on = {true, true};
    else { const bool h = r.bernoulli(0.5); v.on = {!h, h}; }
  };
  const auto interval = [&](SvType t, const SizeRange& range, std::size_t count, bool het_only) {
    for (std::size_t k = 0; k < count; k++) {
      StructuralVariant v;
      v.type = t;
      const std::size_t len = std::max<std::size_t>(1, size_of(range));
      if (!place(len, v.contig, v.begin)) continue;
      v.end = v.begin + len;
      if (het_only) { const bool h = r.bernoulli(0.5); v.on = {!h, h}; }
      else zygosity(v);
      hold(v.contig, v.begin, v.end);
      out.push_back(v);
    }
  };
  // The largest first, while there is room for them.
  interval(SvType::loss_of_heterozygosity, spec.loss_of_heterozygosity, spec.losses_of_heterozygosity, true);
  interval(SvType::deletion, spec.loss, spec.losses, true);
  interval(SvType::tandem_duplication, spec.gain, spec.gains, true);
  interval(SvType::deletion, spec.deletion, spec.deletions, false);
  interval(SvType::tandem_duplication, spec.tandem_duplication, spec.tandem_duplications, false);
  interval(SvType::inversion, spec.inversion, spec.inversions, false);
  for (std::size_t k = 0; k < spec.insertions; k++) {
    StructuralVariant v;
    v.type = SvType::insertion;
    if (!place(1, v.contig, v.begin)) continue;
    v.end = v.begin;
    v.insert_length = std::max<std::size_t>(1, size_of(spec.insertion));
    zygosity(v);
    hold(v.contig, v.begin, v.begin);
    out.push_back(v);
  }
  for (std::size_t k = 0; k < spec.dispersed_duplications; k++) {
    StructuralVariant v;
    v.type = SvType::dispersed_duplication;
    const std::size_t len = std::max<std::size_t>(1, size_of(spec.dispersed_duplication));
    if (!place(len, v.contig, v.begin)) continue;
    v.end = v.begin + len;
    hold(v.contig, v.begin, v.end);
    if (!place(1, v.partner, v.partner_position)) { taken[v.contig].pop_back(); continue; }
    hold(v.partner, v.partner_position, v.partner_position);
    v.inverted = r.bernoulli(0.5);
    zygosity(v);
    out.push_back(v);
  }
  // Translocations: between two different contigs, heterozygous, at most one per chromosome.
  std::vector<bool> translocated(contigs.size(), false);
  for (std::size_t k = 0; k < spec.translocations && contigs.size() > 1; k++) {
    StructuralVariant v;
    v.type = SvType::translocation;
    bool ok = false;
    for (int attempt = 0; attempt < 50 && !ok; attempt++) {
      if (!place(1, v.contig, v.begin) || translocated[v.contig]) continue;
      if (!place(1, v.partner, v.partner_position, v.contig) || translocated[v.partner]) continue;
      ok = true;
    }
    if (!ok) continue;
    v.end = v.begin;
    translocated[v.contig] = translocated[v.partner] = true;
    hold(v.contig, v.begin, v.begin);
    hold(v.partner, v.partner_position, v.partner_position);
    const bool h = r.bernoulli(0.5);
    v.on = {!h, h};
    out.push_back(v);
  }
  for (std::uint32_t c = 0; c < spec.aneuploidy.size() && c < contigs.size(); c++) {
    if (spec.aneuploidy[c] == 0 || translocated[c]) continue;
    StructuralVariant v;
    v.type = SvType::aneuploidy;
    v.contig = c;
    v.begin = 0;
    v.end = contigs[c].size();
    v.copies = spec.aneuploidy[c];
    const bool h = r.bernoulli(0.5);
    v.on = {!h, h};
    out.push_back(v);
  }
  return out;
}

}  // namespace genomics::simulate
