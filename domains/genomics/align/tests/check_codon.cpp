// Checks of codon.hpp: a coding sequence against itself; codon insertions and deletions of a
// planted size found at their place, as whole codons; every codon of each sequence once and in
// order, every gap three gap codes on a codon boundary, on random pairs of related genes; final
// stops and the partial codons of mitochondrial genes left out and counted; synonymous changes
// leaving the protein alignment without gaps; and the vertebrate mitochondrial code (AGA a stop,
// TGA tryptophan).

#include <genomics/align/codon.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::align;

namespace {

const sequence::GeneticCode& standard() { return sequence::genetic_code(1); }

// A random gene of n sense codons (no stop) under the code, and a final stop.
Codes random_gene(calx::random::Stream& r, std::size_t n, const sequence::GeneticCode& code) {
  Codes x;
  while (x.size() < 3 * n) {
    const Code a = static_cast<Code>(r.below(4)), b = static_cast<Code>(r.below(4)), c = static_cast<Code>(r.below(4));
    if (sequence::is_stop(code, a, b, c)) continue;
    x.insert(x.end(), {a, b, c});
  }
  x.insert(x.end(), {sequence::T, sequence::A, sequence::A});
  return x;
}

// The rows without their gaps.
Codes ungapped(const Codes& row) {
  Codes out;
  for (Code c : row)
    if (c != sequence::GAP) out.push_back(c);
  return out;
}

// Every gap three gap codes starting on a codon boundary.
bool codon_gaps(const Codes& row) {
  if (row.size() % 3) return false;
  for (std::size_t k = 0; k < row.size(); k += 3) {
    const int gaps = (row[k] == sequence::GAP) + (row[k + 1] == sequence::GAP) + (row[k + 2] == sequence::GAP);
    if (gaps != 0 && gaps != 3) return false;
  }
  return true;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x636f646e, 0, 0);

  {
    const Codes g = random_gene(r, 120, standard());
    const CodonAlignment a = codon_alignment(g, g, standard());
    const Codes body(g.begin(), g.end() - 3);
    s.require("a gene against itself: its codons aligned one to one, the final stop left out", a.x == body && a.y == body && a.pairs == 120 && a.stop_x && a.stop_y && a.protein_x == a.protein_y);
  }
  {
    // A deletion of codons 40..44 and an insertion of three codons after codon 90 of the copy.
    const Codes g = random_gene(r, 150, standard());
    Codes h(g.begin(), g.begin() + 3 * 40);
    h.insert(h.end(), g.begin() + 3 * 45, g.begin() + 3 * 90);
    const Codes extra = {sequence::G, sequence::C, sequence::A, sequence::T, sequence::G, sequence::G, sequence::A, sequence::A, sequence::C};
    h.insert(h.end(), extra.begin(), extra.end());
    h.insert(h.end(), g.begin() + 3 * 90, g.end());
    const CodonAlignment a = codon_alignment(g, h, standard());
    std::size_t gx = 0, gy = 0;
    for (std::size_t k = 0; k < a.x.size(); k += 3) {
      gx += a.x[k] == sequence::GAP;
      gy += a.y[k] == sequence::GAP;
    }
    // The deleted codons stand against gaps in the copy, at codons 40..44 of the original.
    bool placed = true;
    for (std::size_t c = 40; c < 45; c++) placed = placed && a.y[3 * c] == sequence::GAP;
    s.require("planted indels of codons: five codons deleted and three inserted, found as whole codons at their place", gx == 3 && gy == 5 && placed && codon_gaps(a.x) && codon_gaps(a.y),
              std::to_string(gx) + " and " + std::to_string(gy) + " gap codons");
  }
  {
    bool order = true, frames = true, lengths = true;
    for (int rep = 0; rep < 200; rep++) {
      const std::size_t n = 10 + r.below(200);
      const Codes g = random_gene(r, n, standard());
      // A relative: substitutions at 10%, codon indels at 3%.
      Codes h;
      for (std::size_t c = 0; c + 3 < g.size(); c += 3) {
        if (r.uniform() < 0.03) continue;
        Codes cod(g.begin() + static_cast<std::ptrdiff_t>(c), g.begin() + static_cast<std::ptrdiff_t>(c + 3));
        for (auto& b : cod)
          if (r.uniform() < 0.1) b = static_cast<Code>(r.below(4));
        if (sequence::is_stop(standard(), cod[0], cod[1], cod[2])) cod = {sequence::G, sequence::G, sequence::C};
        h.insert(h.end(), cod.begin(), cod.end());
        if (r.uniform() < 0.03) h.insert(h.end(), {sequence::C, sequence::T, sequence::G});
      }
      if (h.empty()) continue;
      const CodonAlignment a = codon_alignment(g, h, standard());
      order = order && ungapped(a.x) == Codes(g.begin(), g.end() - 3) && ungapped(a.y) == h;
      frames = frames && codon_gaps(a.x) && codon_gaps(a.y);
      lengths = lengths && a.x.size() == a.y.size() && a.x.size() == 3 * a.protein_x.size() && a.protein_x.size() == a.protein_y.size();
    }
    s.require("every codon of each gene once and in order, on 200 related pairs with codon indels", order);
    s.require("every gap three gap codes on a codon boundary", frames);
    s.require("the rows of equal length, three bases to each column of the protein alignment", lengths);
  }
  {
    // Synonymous changes only: the third positions of fourfold-degenerate codons.
    const Codes g = random_gene(r, 80, standard());
    Codes h = g;
    std::size_t changed = 0;
    for (std::size_t c = 0; c + 3 < h.size(); c += 3) {
      const char aa = sequence::amino_acid(standard(), h[c], h[c + 1], h[c + 2]);
      const Code other = static_cast<Code>((h[c + 2] + 1) % 4);
      if (sequence::amino_acid(standard(), h[c], h[c + 1], other) == aa) { h[c + 2] = other; changed++; }
    }
    const CodonAlignment a = codon_alignment(g, h, standard());
    s.require("synonymous changes only: the same protein, no gap, codon against codon", changed > 10 && a.protein_x == a.protein_y && a.pairs == 80 && ungapped(a.x) == a.x,
              std::to_string(changed) + " codons changed");
  }
  {
    // A mitochondrial gene ending in an incomplete stop, T--, and one ending in TAA; table 2.
    const auto& mito = sequence::genetic_code(2);
    Codes g = {sequence::A, sequence::T, sequence::G, sequence::T, sequence::G, sequence::A, sequence::A, sequence::G, sequence::C, sequence::T};
    Codes h = {sequence::A, sequence::T, sequence::G, sequence::T, sequence::G, sequence::A, sequence::A, sequence::G, sequence::C, sequence::T, sequence::A, sequence::A};
    const CodonAlignment a = codon_alignment(g, h, mito);
    s.require("the vertebrate mitochondrial code: TGA is tryptophan, and the partial codon T and the stop TAA are left out", a.protein_x == "MWS" && a.protein_y == "MWS" &&
                  a.partial_x == 1 && a.partial_y == 0 && !a.stop_x && a.stop_y && a.x == a.y && a.x.size() == 9,
              a.protein_x + " / " + a.protein_y);
    const CodonAlignment b = codon_alignment(g, h, standard());
    s.require("the same genes under the standard code read TGA as a stop", b.protein_x == "M*S", b.protein_x);
    const Codes aga = {sequence::A, sequence::T, sequence::G, sequence::C, sequence::C, sequence::C, sequence::A, sequence::G, sequence::A};
    s.require("AGA ends a vertebrate mitochondrial gene (left out as its stop)", codon_alignment(aga, aga, mito).stop_x && codon_alignment(aga, aga, mito).protein_x == "MP");
  }
  return s.done();
}
