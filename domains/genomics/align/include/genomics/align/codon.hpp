// Codon alignment: two coding sequences aligned through their proteins, the gaps of the protein
// alignment carried back to the codons (the method of PAL2NAL; Suyama, Torrents and Bork 2006,
// Nucleic Acids Research 34, W609–W612). Coding sequences change by codons: aligned as bases, a
// gap may cut a codon, and the reading frame, with it the classification of a difference as
// synonymous or not (evolution::selection), is lost. Aligned as proteins, they keep both.
//
// - Codons are read from the first base, whole codons only: a trailing partial codon (the
//   incomplete stops of mitochondrial genes, completed by polyadenylation) is left out and
//   counted, and a final stop codon is left out of both.
// - The proteins, translated under the given genetic code (a stop inside a gene is *, which
//   matches nothing), are aligned globally by Gotoh's algorithm under a protein scoring.
// - A pair of aligned residues becomes a pair of aligned codons; a residue against a gap becomes
//   its codon against three gap codes (sequence::GAP).
//
// Every codon of each sequence appears once and in order, and every gap is three gap codes on a
// codon boundary.

#pragma once

#include <genomics/align/gotoh.hpp>
#include <genomics/align/scoring.hpp>
#include <genomics/align/types.hpp>
#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/genetic_code.hpp>

#include <cstddef>
#include <stdexcept>
#include <string>

namespace genomics::align {

struct CodonAlignment {
  Codes x, y;                          // the aligned codons, sequence::GAP in gaps; equal lengths, multiples of 3
  Alignment protein;                   // the alignment of the proteins
  std::string protein_x, protein_y;    // the proteins aligned, '-' in gaps
  std::size_t codons_x = 0, codons_y = 0;   // the codons aligned of each (final stop left out)
  std::size_t pairs = 0;                    // columns with a codon of each
  std::size_t partial_x = 0, partial_y = 0; // trailing bases left out (0, 1 or 2)
  bool stop_x = false, stop_y = false;      // a final stop left out
};

namespace detail {

// The whole codons of a coding sequence, without a final stop, and its translation.
struct Coding {
  Codes codons;
  std::string protein;
  std::size_t partial = 0;
  bool stop = false;
};

inline Coding coding_part(const Codes& cds, const sequence::GeneticCode& code) {
  Coding c;
  c.partial = cds.size() % 3;
  c.codons.assign(cds.begin(), cds.end() - static_cast<std::ptrdiff_t>(c.partial));
  c.protein = sequence::translate(c.codons, code);
  if (!c.protein.empty() && c.protein.back() == '*') {
    c.protein.pop_back();
    c.codons.resize(c.codons.size() - 3);
    c.stop = true;
  }
  return c;
}

}  // namespace detail

inline CodonAlignment codon_alignment(const Codes& cds_x, const Codes& cds_y, const sequence::GeneticCode& code, const Scoring& protein = blastp()) {
  const detail::Coding a = detail::coding_part(cds_x, code), b = detail::coding_part(cds_y, code);
  if (a.protein.empty() || b.protein.empty()) throw std::invalid_argument("genomics::align::codon_alignment: a coding sequence of at least one codon");
  CodonAlignment out;
  out.partial_x = a.partial;
  out.partial_y = b.partial;
  out.stop_x = a.stop;
  out.stop_y = b.stop;
  out.codons_x = a.protein.size();
  out.codons_y = b.protein.size();
  out.protein = gotoh(sequence::encode_protein(a.protein), sequence::encode_protein(b.protein), protein, Mode::global);
  const Codes gap(3, sequence::GAP);
  std::size_t i = 0, j = 0;  // the next codon of each
  const auto codon = [](const Codes& s, std::size_t k) { return Codes(s.begin() + static_cast<std::ptrdiff_t>(3 * k), s.begin() + static_cast<std::ptrdiff_t>(3 * k + 3)); };
  for (const CigarOp& op : out.protein.cigar)
    for (std::uint32_t n = 0; n < op.length; n++) {
      const bool takes_x = consumes_query(op.op), takes_y = consumes_target(op.op);
      const Codes cx = takes_x ? codon(a.codons, i) : gap, cy = takes_y ? codon(b.codons, j) : gap;
      out.x.insert(out.x.end(), cx.begin(), cx.end());
      out.y.insert(out.y.end(), cy.begin(), cy.end());
      out.protein_x += takes_x ? a.protein[i] : '-';
      out.protein_y += takes_y ? b.protein[j] : '-';
      if (takes_x && takes_y) out.pairs++;
      i += takes_x;
      j += takes_y;
    }
  if (i != out.codons_x || j != out.codons_y) throw std::logic_error("genomics::align::codon_alignment: the protein alignment does not cover both proteins");
  return out;
}

}  // namespace genomics::align
