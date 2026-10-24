// The genetic codes: the twenty-seven translation tables of NCBI's gc.prt (version 4.6,
// retrieved 2026-10-24 from https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt), each as its
// amino acids and its start codons over the 64 codons in NCBI's order, TCAG on each of the
// three bases (TTT, TTC, TTA, TTG, TCT, ..., GGG); and translation, in any frame, of codes or
// of letters, with the ORFs of both strands.
//
// A codon with an ambiguous base translates to the amino acid that every reading of it
// gives, if there is one (GCN is alanine in every code), and to X otherwise. A codon is a
// start when the table marks it M in the start row; translate() writes the first codon of a
// coding sequence as M when it is a start, as GenBank's translations do (the INSDC Feature
// Table, §7.4.5, /translation), and drops the final stop.
//
// An open reading frame here is a stretch of whole codons from a start to the first stop in
// frame, on either strand: the first start after the previous stop in the frame, so that each
// stop ends at most one ORF (the convention of NCBI's ORFfinder). Coordinates are those of the
// forward strand, [begin, end), the stop codon included.

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::sequence {

struct GeneticCode {
  int id = 1;
  std::string_view name;
  std::string_view amino_acids;  // 64 letters, NCBI's codon order
  std::string_view starts;       // 64 characters: M for a start codon, * for a stop, - otherwise
};

inline const std::vector<GeneticCode>& genetic_codes() {
  static const std::vector<GeneticCode> codes = {
    {1, "Standard", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M---------------M----------------------------"},
    {2, "Vertebrate Mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG", "----------**--------------------MMMM----------**---M------------"},
    {3, "Yeast Mitochondrial", "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**----------------------MM---------------M------------"},
    {4, "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--MM------**-------M------------MMMM---------------M------------"},
    {5, "Invertebrate Mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", "---M------**--------------------MMMM---------------M------------"},
    {6, "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear", "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"},
    {9, "Echinoderm Mitochondrial; Flatworm Mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------"},
    {10, "Euplotid Nuclear", "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**-----------------------M----------------------------"},
    {11, "Bacterial, Archaeal and Plant Plastid", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M------------MMMM---------------M------------"},
    {12, "Alternative Yeast Nuclear", "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*----M---------------M----------------------------"},
    {13, "Ascidian Mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", "---M------**----------------------MM---------------M------------"},
    {14, "Alternative Flatworm Mitochondrial", "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "-----------*-----------------------M----------------------------"},
    {15, "Blepharisma Macronuclear", "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------*---*--------------------M----------------------------"},
    {16, "Chlorophycean Mitochondrial", "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------*---*--------------------M----------------------------"},
    {21, "Trematode Mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------"},
    {22, "Scenedesmus obliquus Mitochondrial", "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "------*---*---*--------------------M----------------------------"},
    {23, "Thraustochytrium Mitochondrial", "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--*-------**--*-----------------M--M---------------M------------"},
    {24, "Rhabdopleuridae Mitochondrial", "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "---M------**-------M---------------M---------------M------------"},
    {25, "Candidate Division SR1 and Gracilibacteria", "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**-----------------------M---------------M------------"},
    {26, "Pachysolen tannophilus Nuclear", "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*----M---------------M----------------------------"},
    {27, "Karyorelict Nuclear", "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"},
    {28, "Condylostoma Nuclear", "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*--------------------M----------------------------"},
    {29, "Mesodinium Nuclear", "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"},
    {30, "Peritrich Nuclear", "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"},
    {31, "Blastocrithidia Nuclear", "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**-----------------------M----------------------------"},
    {32, "Balanophoraceae Plastid", "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------*---*----M------------MMMM---------------M------------"},
    {33, "Cephalodiscidae Mitochondrial", "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "---M-------*-------M---------------M---------------M------------"},
  };
  return codes;
}

inline const GeneticCode& genetic_code(int id) {
  for (const auto& g : genetic_codes())
    if (g.id == id) return g;
  throw std::invalid_argument("genomics::sequence::genetic_code: NCBI has no translation table " + std::to_string(id));
}

// The position of a codon in NCBI's TCAG order, or -1 if a base is not A, C, G or T.
constexpr int codon_index(Code b1, Code b2, Code b3) {
  if (b1 > 3 || b2 > 3 || b3 > 3) return -1;
  constexpr int rank[4] = {2, 1, 3, 0};  // A, C, G, T in the order T, C, A, G
  return 16 * rank[b1] + 4 * rank[b2] + rank[b3];
}
// The codon at an index of NCBI's order, as three codes.
constexpr std::array<Code, 3> codon_of(int index) {
  constexpr Code base[4] = {T, C, A, G};
  return {base[(index >> 4) & 3], base[(index >> 2) & 3], base[index & 3]};
}

// The amino acid of a codon (a letter), resolving ambiguous bases when every reading agrees.
inline char amino_acid(const GeneticCode& code, Code b1, Code b2, Code b3) {
  const int i = codon_index(b1, b2, b3);
  if (i >= 0) return code.amino_acids[i];
  char found = 0;
  for (Code x = 0; x < 4; x++) {
    if (b1 < 4 && x != b1) continue;
    for (Code y = 0; y < 4; y++) {
      if (b2 < 4 && y != b2) continue;
      for (Code z = 0; z < 4; z++) {
        if (b3 < 4 && z != b3) continue;
        const char a = code.amino_acids[codon_index(x, y, z)];
        if (found && a != found) return 'X';
        found = a;
      }
    }
  }
  return found ? found : 'X';
}
inline bool is_start(const GeneticCode& code, Code b1, Code b2, Code b3) {
  const int i = codon_index(b1, b2, b3);
  return i >= 0 && code.starts[i] == 'M';
}
inline bool is_stop(const GeneticCode& code, Code b1, Code b2, Code b3) {
  const int i = codon_index(b1, b2, b3);
  return i >= 0 && code.amino_acids[i] == '*';
}

// The translation of x from `frame` (0, 1 or 2) to its end, in letters. With `coding`, the first
// codon is written M when it is a start codon, and a final stop is dropped.
inline std::string translate(const Codes& x, const GeneticCode& code, std::size_t frame = 0, bool coding = false) {
  std::string out;
  if (frame >= x.size()) return out;
  out.reserve((x.size() - frame) / 3);
  for (std::size_t i = frame; i + 3 <= x.size(); i += 3) out += amino_acid(code, x[i], x[i + 1], x[i + 2]);
  if (coding && !out.empty()) {
    if (is_start(code, x[frame], x[frame + 1], x[frame + 2])) out[0] = 'M';
    if (out.back() == '*') out.pop_back();
  }
  return out;
}

struct Orf {
  std::size_t begin = 0, end = 0;  // forward-strand coordinates, the stop codon included
  bool reverse = false;
  int frame = 0;                   // the frame of the start on its own strand, 0 to 2
  std::size_t codons() const { return (end - begin) / 3; }
};

// The ORFs of both strands with at least `min_codons` codons (the stop included). With
// `any_start`, every start codon of the table opens a frame; otherwise only ATG.
inline std::vector<Orf> open_reading_frames(const Codes& x, const GeneticCode& code, std::size_t min_codons = 100, bool any_start = false) {
  std::vector<Orf> out;
  const std::size_t n = x.size();
  const auto scan = [&](const Codes& s, bool reverse) {
    for (int frame = 0; frame < 3; frame++) {
      std::size_t open = SIZE_MAX;
      for (std::size_t i = static_cast<std::size_t>(frame); i + 3 <= n; i += 3) {
        const Code a = s[i], b = s[i + 1], c = s[i + 2];
        if (open == SIZE_MAX && (any_start ? is_start(code, a, b, c) : (a == A && b == T && c == G))) open = i;
        if (is_stop(code, a, b, c)) {
          if (open != SIZE_MAX && (i + 3 - open) / 3 >= min_codons) {
            Orf o;
            o.reverse = reverse;
            o.frame = frame;
            o.begin = reverse ? n - (i + 3) : open;
            o.end = reverse ? n - open : i + 3;
            out.push_back(o);
          }
          open = SIZE_MAX;
        }
      }
    }
  };
  scan(x, false);
  scan(reverse_complement(x), true);
  std::sort(out.begin(), out.end(), [](const Orf& p, const Orf& q) { return p.begin != q.begin ? p.begin < q.begin : p.end < q.end; });
  return out;
}

}  // namespace genomics::sequence
