// Checks of the simulation of rearranged genomes and of paired reads: every base of every
// rearranged chromosome lifted to the reference is the reference's base there (complemented on a
// reversed block), and every inserted base is novel; the lengths of the chromosomes change by
// what the variants add and remove; the copy number of each reference base, counted by brute
// force over the lifted bases, is the step function the sample reports; each kind of variant
// makes its junctions, with the breakends the theory note gives; a copy-neutral loss of
// heterozygosity makes none and leaves the copy number at two, from one haplotype; random
// variants are as many as asked, valid together and apart; paired reads have the fragment law
// asked for, FR orientation, mates that are their chromosome's bases at their truth up to their
// errors, junction flags where junctions are, and a GC bias that moves the fragments' GC toward
// the optimum; a shard of pairs is the same pairs of the whole sample.

#include <genomics/sequence/alphabet.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/pairs.hpp>
#include <genomics/simulate/structural.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <map>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::simulate;
using sequence::Codes;

int main() {
  calx::verification::Suite s;

  // Two contigs of a random genome, a run of N in the second.
  GenomeSpec spec;
  spec.length = 300'000;
  spec.gc = 0.45;
  std::vector<Codes> contigs{make_genome(spec, 3).sequence, make_genome(spec, 4).sequence};
  for (std::size_t i = 150'000; i < 160'000; i++) contigs[1][i] = sequence::N;

  // One of every kind, by hand.
  std::vector<StructuralVariant> svs;
  const auto add = [&](SvType t, std::uint32_t c, std::size_t a, std::size_t b, std::array<bool, 2> on) {
    StructuralVariant v;
    v.type = t;
    v.contig = c;
    v.begin = a;
    v.end = b;
    v.on = on;
    svs.push_back(v);
    return &svs.back();
  };
  add(SvType::deletion, 0, 10'000, 12'000, {true, false});
  add(SvType::deletion, 0, 20'000, 20'300, {true, true});
  add(SvType::tandem_duplication, 0, 30'000, 33'000, {false, true})->copies = 2;
  add(SvType::inversion, 0, 40'000, 45'000, {true, false});
  add(SvType::insertion, 0, 50'000, 50'000, {false, true})->insert_length = 120;
  {
    auto* v = add(SvType::dispersed_duplication, 0, 60'000, 61'500, {true, false});
    v->partner = 1;
    v->partner_position = 20'000;
    v->inverted = true;
  }
  {
    auto* v = add(SvType::translocation, 0, 200'000, 200'000, {false, true});
    v->partner = 1;
    v->partner_position = 250'000;
  }
  add(SvType::loss_of_heterozygosity, 0, 100'000, 150'000, {false, true});
  add(SvType::deletion, 1, 60'000, 120'000, {false, true});            // a large loss
  add(SvType::tandem_duplication, 1, 170'000, 230'000, {true, false});  // a large gain

  MutationRates none;
  none.snv = 0;
  none.indel = 0;
  const StructuralSample x = make_structural_sample(contigs, svs, none, 5, {"chrA", "chrB"});

  // ---- the liftover: every lifted base is the reference's ------------------------------------------
  {
    std::size_t lifted = 0, wrong = 0, inserted = 0;
    std::vector<std::vector<int>> count(2);
    for (std::size_t c = 0; c < 2; c++) count[c].assign(contigs[c].size(), 0);
    for (const auto& ch : x.chromosomes)
      for (std::size_t i = 0; i < ch.sequence.size(); i++) {
        const Locus l = ch.map.locate(i);
        if (l.inserted) { inserted++; continue; }
        lifted++;
        const auto ref = contigs[l.contig][l.position];
        wrong += (l.reverse ? sequence::complement(ref) : ref) != ch.sequence[i];
        count[l.contig][l.position]++;
      }
    s.require("every base of every rearranged chromosome lifts to the reference base it was taken from", wrong == 0 && lifted > 0,
              std::to_string(wrong) + " of " + std::to_string(lifted) + " wrong");
    s.require("the inserted bases are exactly the novel insertion's", inserted == 120, std::to_string(inserted));
    std::size_t disagree = 0;
    for (std::uint32_t c = 0; c < 2; c++)
      for (std::size_t p = 0; p < contigs[c].size(); p++) {
        const auto k = x.copies_at(c, p);
        disagree += k[0] + k[1] != count[c][p];
      }
    s.require("the copy number of every reference base is the number of lifted bases there (brute force)", disagree == 0, std::to_string(disagree));
    const auto at = [&](std::uint32_t c, std::size_t p) { const auto k = x.copies_at(c, p); return std::to_string(k[0]) + "+" + std::to_string(k[1]); };
    s.require("copy numbers by haplotype: 1+0 in a heterozygous deletion, 0+0 in a homozygous one, 1+3 in a triplication, 2+0 after loss of heterozygosity of haplotype 2",
              at(0, 11'000) == "0+1" && at(0, 20'100) == "0+0" && at(0, 31'000) == "1+3" && at(0, 120'000) == "2+0" && at(0, 99'999) == "1+1",
              at(0, 11'000) + " " + at(0, 20'100) + " " + at(0, 31'000) + " " + at(0, 120'000));
    s.require("the dispersed copy adds a copy of its source; the large loss and gain change one haplotype",
              at(0, 60'500) == "2+1" && at(1, 90'000) == "1+0" && at(1, 200'000) == "2+1" && at(1, 140'000) == "1+1",
              at(0, 60'500) + " " + at(1, 90'000) + " " + at(1, 200'000));
  }

  // ---- lengths and junctions ---------------------------------------------------------------------
  {
    std::size_t total = 0;
    for (const auto& ch : x.chromosomes) total += ch.sequence.size();
    const long expected = 2L * 600'000 - 2'000 - 2 * 300 + 2 * 3'000 + 120 + 1'500 - 60'000 + 60'000;
    s.require("the chromosomes' lengths change by what the variants add and remove", static_cast<long>(total) == expected, std::to_string(total) + " vs " + std::to_string(expected));
    std::map<std::uint32_t, int> per;
    for (const auto& j : x.junctions) per[j.variant]++;
    s.require("junctions: 1 per deletion haplotype, 2 for a triplication, 2 for an inversion, 1 for an insertion, 2 for a dispersed copy and a translocation, none for loss of heterozygosity",
              per[0] == 1 && per[1] == 2 && per[2] == 2 && per[3] == 2 && per[4] == 1 && per[5] == 2 && per[6] == 2 && per[7] == 0 && per[8] == 1 && per[9] == 1,
              [&] { std::string t; for (auto& [k, n] : per) t += std::to_string(k) + ":" + std::to_string(n) + " "; return t; }());
    bool breakends = true;
    for (const auto& j : x.junctions) {
      const auto& v = x.truth[j.variant];
      switch (v.type) {
        case SvType::deletion: breakends = breakends && j.left == v.begin - 1 && j.left_plus && j.right == v.end && !j.right_plus; break;
        case SvType::tandem_duplication: breakends = breakends && j.left == v.end - 1 && j.left_plus && j.right == v.begin && !j.right_plus; break;
        case SvType::inversion:
          breakends = breakends && ((j.left == v.begin - 1 && j.left_plus && j.right == v.end - 1 && j.right_plus) || (j.left == v.begin && !j.left_plus && j.right == v.end && !j.right_plus));
          break;
        case SvType::insertion: breakends = breakends && j.left == v.begin - 1 && j.right == v.begin && j.inserted == v.sequence; break;
        case SvType::translocation:
          breakends = breakends && j.left_contig != j.right_contig &&
                      ((j.left_contig == 0 && j.left == v.begin - 1 && j.right == v.partner_position) || (j.left_contig == 1 && j.left == v.partner_position - 1 && j.right == v.begin));
          break;
        default: break;
      }
    }
    s.require("every junction's breakends are those of its kind (deletion (a−1,+)~(b,−), tandem (b−1,+)~(a,−), inversion (a−1,+)~(b−1,+) and (a,−)~(b,−), ...)", breakends);
    bool der = false;
    for (const auto& ch : x.chromosomes) der = der || ch.name.rfind("der(", 0) == 0;
    s.require("the translocation makes two derivative chromosomes", der && x.chromosomes.size() == 4);
  }

  // ---- aneuploidy and small variants -------------------------------------------------------------
  {
    std::vector<StructuralVariant> tri(1);
    tri[0].type = SvType::aneuploidy;
    tri[0].contig = 1;
    tri[0].copies = 1;
    tri[0].on = {true, false};
    const auto y = make_structural_sample(contigs, tri, MutationRates{}, 6);
    s.require("a trisomy adds a third chromosome of its contig, copy number 2+1 everywhere", y.chromosomes.size() == 5 && y.copies_at(1, 5) == std::array<int, 2>{2, 1} && y.copies_at(0, 5) == std::array<int, 2>{1, 1});
    tri[0].copies = -1;
    const auto z = make_structural_sample(contigs, tri, MutationRates{}, 6);
    s.require("a monosomy removes one, copy number 0+1", z.chromosomes.size() == 3 && z.copies_at(1, 5) == std::array<int, 2>{0, 1});
    // With small variants, a lifted base differs from the reference at about the SNV rate.
    std::size_t lifted = 0, wrong = 0;
    for (const auto& ch : y.chromosomes)
      for (std::size_t i = 0; i < ch.sequence.size(); i++) {
        const Locus l = ch.map.locate(i);
        if (l.inserted) continue;
        lifted++;
        wrong += contigs[l.contig][l.position] != ch.sequence[i];
      }
    // An SNV is on a given haplotype with probability 0.4 + 0.6/2 (homozygous, or heterozygous on it).
    s.check("with SNVs at 10⁻³ a lifted base differs from the reference at 0.7 · 10⁻³", std::fabs(static_cast<double>(wrong) / static_cast<double>(lifted) - 0.7e-3) / 0.7e-3, 0.15);
  }

  // ---- random variants ----------------------------------------------------------------------------
  {
    StructuralSpec rs;
    rs.deletions = 8;
    rs.tandem_duplications = 4;
    rs.inversions = 3;
    rs.insertions = 4;
    rs.dispersed_duplications = 2;
    rs.translocations = 1;
    rs.losses = 1;
    rs.gains = 1;
    rs.losses_of_heterozygosity = 1;
    rs.gain = {30'000, 50'000};
    rs.loss = {30'000, 50'000};
    rs.loss_of_heterozygosity = {30'000, 60'000};
    const auto r = random_structural(contigs, rs, 9);
    bool valid = true;
    try { validate_structural(contigs, r); } catch (const std::invalid_argument&) { valid = false; }
    s.require("random variants: as many as asked (25), valid together", r.size() == 25 && valid, std::to_string(r.size()));
    s.require("random variants are a function of the seed", random_structural(contigs, rs, 9).size() == r.size() && random_structural(contigs, rs, 9)[3].begin == r[3].begin);
    bool clear = true;
    for (const auto& v : r)
      if (v.contig == 1 && v.begin < 160'000 + rs.margin && v.end + rs.margin > 150'000) clear = false;
    s.require("random variants keep their margin from the run of N", clear);
    const auto y = make_structural_sample(contigs, r, MutationRates{}, 10);
    s.require("a sample of random variants is made", y.chromosomes.size() >= 4);
  }

  // ---- paired reads --------------------------------------------------------------------------------
  {
    PairModel pm;
    pm.fragment.mean = 450;
    pm.fragment.sd = 45;
    pm.read.indel = 0;  // substitutions only: a mate is its truth base for base
    const std::size_t n = pairs_for_coverage(x.chromosomes, 10, pm.read.length);
    const auto pairs = read_pairs(x.chromosomes, n, pm, 11);
    double m = 0, m2 = 0;
    std::size_t fr = 0, junction = 0, mismatches = 0, errors = 0, bases = 0;
    bool truth = true;
    for (const auto& p : pairs) {
      const double f = static_cast<double>(p.fragment_length);
      m += f;
      m2 += f * f;
      fr += p.mates[0].reverse != p.mates[1].reverse;
      const auto& ch = x.chromosomes[p.chromosome];
      for (const auto& mate : p.mates) {
        junction += mate.junction;
        // The mate against its chromosome bases (no indel errors in the model below → base for base).
        Codes t(ch.sequence.begin() + static_cast<std::ptrdiff_t>(mate.begin), ch.sequence.begin() + static_cast<std::ptrdiff_t>(mate.end));
        if (mate.reverse) t = sequence::reverse_complement(t);
        if (t.size() == mate.bases.size())
          for (std::size_t i = 0; i < t.size(); i++) mismatches += t[i] != mate.bases[i];
        errors += mate.errors;
        bases += mate.bases.size();
        truth = truth && (mate.first.position == ch.map.locate(mate.begin).position) && mate.end - mate.begin >= mate.bases.size() - mate.errors;
      }
    }
    const double N = static_cast<double>(pairs.size());
    m /= N;
    const double sd = std::sqrt(m2 / N - m * m);
    s.require("pairs: as many as asked for 10×", pairs.size() == n, std::to_string(pairs.size()));
    s.check("the fragments' mean is 450", std::fabs(m - 450), 1.5, "bases");
    s.check("and their standard deviation 45", std::fabs(sd - 45), 1.0, "bases");
    s.require("every pair's mates are on opposite strands (FR)", fr == pairs.size());
    s.require("every mate is its chromosome's bases at its truth but for the errors made", mismatches == errors && bases > 0, std::to_string(mismatches) + " vs " + std::to_string(errors));
    s.require("each mate's truth locus is its chromosome's liftover at its interval", truth);
    // Expected mates across a junction: a read of 150 crosses a junction when it starts in the 149 bases before it.
    std::size_t junctions = 0;
    for (const auto& ch : x.chromosomes) junctions += ch.junctions.size();
    const double expected = static_cast<double>(junctions) * 149.0 * 2.0 * N / static_cast<double>(x.length());
    s.check("mates cross junctions as often as junctions predict", std::fabs(static_cast<double>(junction) - expected) / expected, 0.3);
    s.require("a shard of pairs is the same pairs of the whole sample", [&] {
      const auto part = read_pairs(x.chromosomes, 100, pm, 11, 500);
      return part.size() == 100 && part[7].mates[1].bases == pairs[507].mates[1].bases && part[7].fragment_begin == pairs[507].fragment_begin;
    }());
    // GC bias: the fragments' mean GC moves toward the optimum from the genome's 0.45.
    PairModel biased = pm;
    biased.fragment.gc_bias = 1;
    biased.fragment.gc_optimum = 0.40;
    biased.fragment.gc_width = 0.03;
    double g0 = 0, g1 = 0;
    const auto a = read_pairs(x.chromosomes, 4000, pm, 12), b = read_pairs(x.chromosomes, 4000, biased, 12);
    for (const auto& p : a) g0 += gc_or_n(x.chromosomes[p.chromosome].sequence, p.fragment_begin, p.fragment_begin + p.fragment_length);
    for (const auto& p : b) g1 += gc_or_n(x.chromosomes[p.chromosome].sequence, p.fragment_begin, p.fragment_begin + p.fragment_length);
    g0 /= static_cast<double>(a.size());
    g1 /= static_cast<double>(b.size());
    s.require("a GC-biased library draws fragments of GC nearer its optimum (0.40) than an unbiased one", g1 < g0 - 0.01, std::to_string(g0) + " → " + std::to_string(g1));
    // No fragment from the run of N.
    bool none_from_n = true;
    for (const auto& p : a) none_from_n = none_from_n && gc_or_n(x.chromosomes[p.chromosome].sequence, p.fragment_begin, p.fragment_begin + p.fragment_length) >= 0;
    s.require("no fragment touches an N", none_from_n);
  }
  return s.done();
}
