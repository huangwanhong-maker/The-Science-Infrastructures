// Checks of the simulator: a genome is a function of its seed, of the length and GC asked for,
// with its repeat copies where it says, diverged at the rate asked; a diploid individual's
// haplotypes are the reference with exactly its variants applied, its liftover maps unmutated
// bases home, and its rates and zygosity come out as asked; short reads are the length asked,
// a shard of a sample equals the same reads of the whole sample, each read matches its source
// at its true place up to its errors, the error rate rises along the read as the model says,
// and qualities track errors; long reads have the mean length and error rate asked; a family's
// true alignment, with its gaps removed, gives back each member, and its columns are in an
// order every member respects.

#include <genomics/simulate/family.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/mutate.hpp>
#include <genomics/simulate/reads.hpp>
#include <genomics/sequence/composition.hpp>

#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::simulate;

int main() {
  calx::verification::Suite s;

  // ---- genomes ----------------------------------------------------------------------------------
  GenomeSpec spec;
  spec.length = 200'000;
  spec.gc = 0.62;
  spec.repeats = {4, 5, 3000, 0.01};
  spec.tandem = {2, 37, 40, 0.02};
  const Genome g = make_genome(spec, 11);
  {
    s.require("a genome is a function of its seed", make_genome(spec, 11).sequence == g.sequence && make_genome(spec, 12).sequence != g.sequence);
    s.check("its GC content is the one asked for (0.62)", sequence::count_ambiguous(g.sequence) == 0 ? std::fabs(sequence::gc_fraction(g.sequence) - 0.62) : 1, 0.005);
    s.require("it has its 20 interspersed copies and 2 tandem arrays, never overlapping", g.repeats.size() == 22 && [&] {
      for (std::size_t i = 1; i < g.repeats.size(); i++)
        if (g.repeats[i].begin < g.repeats[i - 1].end) return false;
      return true;
    }());
    // Copies of one element differ by about twice the divergence (each diverged from the element).
    std::map<std::size_t, std::vector<RepeatCopy>> by;
    for (const auto& c : g.repeats)
      if (!c.tandem) by[c.element].push_back(c);
    double diff = 0, sites = 0;
    for (auto& [e, copies] : by) {
      const auto seq_of = [&](const RepeatCopy& c) {
        sequence::Codes x(g.sequence.begin() + static_cast<std::ptrdiff_t>(c.begin), g.sequence.begin() + static_cast<std::ptrdiff_t>(c.end));
        return c.reverse ? sequence::reverse_complement(x) : x;
      };
      const auto a = seq_of(copies[0]);
      for (std::size_t k = 1; k < copies.size(); k++) {
        const auto b = seq_of(copies[k]);
        for (std::size_t i = 0; i < a.size(); i++) diff += a[i] != b[i];
        sites += static_cast<double>(a.size());
      }
    }
    s.check("copies of an element differ at 2p − 4p²/3 for divergence p = 0.01", diff / sites - (0.02 - 4 * 0.0001 / 3), 0.003);
  }

  // ---- a diploid individual ----------------------------------------------------------------------------------
  const sequence::Codes& ref = g.sequence;
  MutationRates rates;
  rates.snv = 2e-3;
  rates.indel = 4e-4;
  const Diploid d = mutate_diploid(ref, rates, 5);
  {
    // Apply the truth, as a VCF would be applied, and compare.
    bool applied = true;
    for (int h = 0; h < 2; h++) {
      std::string hap;
      std::size_t at = 0;
      const std::string text = sequence::decode_dna(ref);
      for (const auto& v : d.truth) {
        if (!v.on[h]) continue;
        applied = applied && text.compare(v.position, v.ref.size(), v.ref) == 0;
        hap += text.substr(at, v.position - at) + v.alt;
        at = v.position + v.ref.size();
      }
      hap += text.substr(at);
      applied = applied && hap == sequence::decode_dna(d.haplotype[h]);
    }
    s.require("the haplotypes are the reference with exactly their variants applied, each REF matching the reference", applied);
    std::size_t snv = 0, indel = 0, het = 0;
    for (const auto& v : d.truth) {
      snv += v.snv();
      indel += !v.snv();
      het += v.heterozygous();
    }
    s.check("SNVs at the rate asked (per base, 2e-3)", static_cast<double>(snv) / static_cast<double>(ref.size()) / 2e-3 - 1, 0.1);
    // Indels are rarer: their rate is measured on 2 Mb (about 800 expected, a Poisson spread of 3.5%).
    const auto big = mutate_diploid(random_sequence(2'000'000, 0.5, 8), rates, 6);
    std::size_t big_indels = 0;
    for (const auto& v : big.truth) big_indels += !v.snv();
    s.check("indels at the rate asked (4e-4, on 2 Mb)", static_cast<double>(big_indels) / 2e6 / 4e-4 - 1, 0.15);
    s.require("this genome has indels and SNVs", indel > 0 && snv > 0);
    s.check("heterozygous at the fraction asked (0.6)", static_cast<double>(het) / static_cast<double>(d.truth.size()) - 0.6, 0.05);
    // The liftover sends every base of a block home.
    bool home = true;
    for (int h = 0; h < 2; h++)
      for (const auto& b : d.liftover[h].blocks)
        for (std::size_t k = 0; k < b.length; k += 97) {
          const bool snv_here = std::any_of(d.truth.begin(), d.truth.end(), [&](const Variant& v) { return v.on[h] && v.snv() && v.position == b.ref + k; });
          home = home && d.liftover[h].to_reference(b.hap + k) == b.ref + k && (snv_here || d.haplotype[h][b.hap + k] == ref[b.ref + k]);
        }
    s.require("the liftover maps the bases of each haplotype to their places on the reference", home);
  }

  // ---- short reads ----------------------------------------------------------------------------------
  {
    const Source src = Source::diploid({d});
    ShortReadModel m;
    m.length = 150;
    m.error_first = 0.002;
    m.error_last = 0.02;
    const auto reads = short_reads(src, 4000, m, 9);
    const auto shard = short_reads(src, 500, m, 9, 1200);
    bool same = shard.size() == 500;
    for (std::size_t i = 0; i < shard.size() && same; i++) same = shard[i].bases == reads[1200 + i].bases && shard[i].begin == reads[1200 + i].begin;
    s.require("a shard of a sample (reads 1200 to 1699) equals the same reads of the whole sample", same);
    bool lengths = true;
    for (const auto& r : reads) lengths = lengths && r.bases.size() == 150 && r.quality.size() == 150;
    s.require("every short read has 150 bases and 150 qualities", lengths);
    // Errors along the read, by comparison with the haplotype it came from (reads without indels).
    std::vector<double> err(150, 0), n(150, 0);
    std::size_t total_errors = 0, counted = 0;
    double q_err = 0, q_ok = 0, n_err = 0, n_ok = 0;
    bool placed = true;
    for (const auto& r : reads) {
      if (r.errors == 0) { /* exact */ }
      // Reads free of indels align without gaps to the haplotype.
      const auto& hap = d.haplotype[r.haplotype];
      // Find the haplotype interval: the read covers [b, b + 150) of the haplotype, mapped from the reference.
      std::size_t hb = SIZE_MAX;
      for (const auto& blk : d.liftover[r.haplotype].blocks)
        if (r.begin >= blk.ref && r.begin < blk.ref + blk.length) { hb = blk.hap + (r.begin - blk.ref); break; }
      if (hb == SIZE_MAX || hb + 150 > hap.size()) continue;
      sequence::Codes t(hap.begin() + static_cast<std::ptrdiff_t>(hb), hap.begin() + static_cast<std::ptrdiff_t>(hb + 150));
      if (r.reverse) t = sequence::reverse_complement(t);
      std::size_t mism = 0;
      for (std::size_t i = 0; i < 150; i++) mism += t[i] != r.bases[i];
      if (mism != r.errors) continue;  // a read with an indel, or a liftover across one: skipped
      counted++;
      total_errors += mism;
      for (std::size_t i = 0; i < 150; i++) {
        n[i]++;
        const bool e = t[i] != r.bases[i];
        err[i] += e;
        (e ? q_err : q_ok) += sequence::phred_of(r.quality[i]);
        (e ? n_err : n_ok) += 1;
      }
      placed = placed && r.end > r.begin && r.end - r.begin <= 170;
    }
    double first = 0, last = 0, nf = 0, nl = 0;
    for (std::size_t i = 0; i < 20; i++) { first += err[i]; nf += n[i]; last += err[149 - i]; nl += n[149 - i]; }
    s.require("most reads match their source exactly up to their substitutions", counted > 3000, std::to_string(counted) + " of 4000");
    s.check("the error rate near the first bases (0.002 to 0.0044)", first / nf / 0.0032 - 1, 0.25);
    s.check("the error rate near the last bases (0.0176 to 0.02)", last / nl / 0.0188 - 1, 0.15);
    s.require("qualities track errors: erroneous bases have lower mean quality", q_err / n_err < q_ok / n_ok, std::to_string(q_err / n_err) + " < " + std::to_string(q_ok / n_ok));
    s.require("each read's reference interval is about its length", placed);
  }

  // ---- long reads ----------------------------------------------------------------------------------
  {
    const Source src = Source::haploid({g.sequence});
    LongReadModel m;
    m.mean_length = 6000;
    m.error = 0.1;
    const auto reads = long_reads(src, 600, m, 3);
    double span = 0, errors = 0, bases = 0;
    for (const auto& r : reads) {
      span += static_cast<double>(r.end - r.begin);
      errors += static_cast<double>(r.errors);
      bases += static_cast<double>(r.end - r.begin);
    }
    s.check("long reads: the mean span is the mean length asked (6000)", span / static_cast<double>(reads.size()) / 6000 - 1, 0.1);
    s.check("long reads: errors at the rate asked (0.1 per base of template)", errors / bases / (0.1 / (1 - 0.1 * 0.35)) - 1, 0.05);
  }

  // ---- a family ----------------------------------------------------------------------------------
  {
    const auto root = random_sequence(400, 0.5, 21);
    const Family f = evolve_family(root, 7, 0.08, 0.02, 4);
    bool back = f.sequences.size() == 7 && f.alignment.size() == 7;
    std::size_t width = f.alignment.empty() ? 0 : f.alignment[0].size();
    for (std::size_t m = 0; m < f.sequences.size(); m++) {
      sequence::Codes ungapped;
      for (auto c : f.alignment[m]) if (c != sequence::GAP) ungapped.push_back(c);
      back = back && ungapped == f.sequences[m] && f.alignment[m].size() == width;
    }
    s.require("a family's true alignment, its gaps removed, gives back each member", back);
    bool no_empty = true;
    for (std::size_t c = 0; c < width; c++) {
      bool any = false;
      for (const auto& row : f.alignment) any = any || row[c] != sequence::GAP;
      no_empty = no_empty && any;
    }
    s.require("every column of the true alignment holds a base of some member", no_empty);
    s.require("the family has indels: its alignment is wider than any member", width > std::max_element(f.sequences.begin(), f.sequences.end(), [](const auto& a, const auto& b) { return a.size() < b.size(); })->size());
    s.require("the family's tree is written in Newick with its seven leaves", std::count(f.newick.begin(), f.newick.end(), 'm') == 7 && f.newick.back() == ';', f.newick);
  }

  return s.done();
}
