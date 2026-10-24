// Checks of the alphabets, packed DNA, the genetic codes and composition: letters and codes
// round-trip, U reads as T and the complement is an involution; the fifteen IUPAC sets and
// their complements; packed sequences with runs of N of every length against their codes,
// base by base, count by count, k-mer by k-mer and reverse complement by reverse complement;
// the twenty-seven NCBI tables well formed, the standard code's stops and starts, the
// vertebrate mitochondrial code's AGA, AGG, TGA and ATA, the alternative starts of table 11,
// ambiguous codons resolved only when every reading agrees; a peptide back-translated and
// translated again; ORFs that start, stop and are maximal on both strands; composition: GC,
// CpG o/e near 1 in i.i.d. sequence and 0 without CpG, entropy of k-mers, the cumulative
// skew's minimum at a planted origin, RSCU 1 and N_c = 61 for even usage, N_c = 20 for one
// codon per amino acid.

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/composition.hpp>
#include <genomics/sequence/genetic_code.hpp>
#include <genomics/sequence/kmer.hpp>
#include <genomics/sequence/packed.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace genomics::sequence;

namespace {

Codes random_dna(calx::random::Stream& r, std::size_t n, double n_rate = 0) {
  Codes x(n);
  for (auto& c : x) c = r.uniform() < n_rate ? N : static_cast<Code>(r.below(4));
  return x;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x73657121, 0, 0);

  // ---- alphabets ----------------------------------------------------------------------------------
  {
    s.require("DNA letters round-trip through codes, in either case", decode_dna(encode_dna("ACGTNacgtn")) == "ACGTNACGTN");
    s.require("U reads as T, and other letters as N", decode_dna(encode_dna("UuRYX-")) == "TTNNNN");
    bool inv = true;
    for (Code c = 0; c < 5; c++) inv = inv && complement(complement(c)) == c;
    s.require("the complement is an involution, A↔T, C↔G, N fixed", inv && complement(A) == T && complement(C) == G && complement(N) == N);
    const std::string text = "ACGGTTNACRYKM";
    s.require("the reverse complement of letters keeps IUPAC sets: ACGGTTNACRYKM → KMRYGTNAACCGT", reverse_complement(text) == "KMRYGTNAACCGT", reverse_complement(text));
    s.require("the reverse complement of codes agrees with that of letters", decode_dna(reverse_complement(encode_dna("ACGGTTNAC"))) == reverse_complement(std::string("ACGGTTNAC")));
    bool iupac = true;
    const std::map<char, std::string> sets = {{'A', "A"}, {'C', "C"}, {'G', "G"}, {'T', "T"}, {'R', "AG"}, {'Y', "CT"}, {'S', "CG"}, {'W', "AT"}, {'K', "GT"}, {'M', "AC"},
                                              {'B', "CGT"}, {'D', "AGT"}, {'H', "ACT"}, {'V', "ACG"}, {'N', "ACGT"}};
    for (const auto& [letter, bases] : sets) {
      std::uint8_t m = 0;
      for (char b : bases) m |= static_cast<std::uint8_t>(1u << dna_code(b));
      iupac = iupac && iupac_mask(letter) == m && iupac_letter(m) == letter && iupac_mask(static_cast<char>(letter + 32)) == m;
      for (Code b = 0; b < 4; b++) iupac = iupac && iupac_matches(letter, b) == (bases.find(DNA_LETTERS[b]) != std::string::npos);
    }
    s.require("the fifteen IUPAC letters are the sets of the NC-IUB recommendations", iupac);
    bool comp = true;
    for (std::uint8_t m = 1; m < 16; m++) comp = comp && iupac_complement(iupac_complement(m)) == m;
    comp = comp && iupac_letter(iupac_complement(iupac_mask('R'))) == 'Y' && iupac_letter(iupac_complement(iupac_mask('K'))) == 'M' &&
           iupac_letter(iupac_complement(iupac_mask('B'))) == 'V' && iupac_letter(iupac_complement(iupac_mask('D'))) == 'H' &&
           iupac_letter(iupac_complement(iupac_mask('S'))) == 'S' && iupac_letter(iupac_complement(iupac_mask('W'))) == 'W';
    s.require("IUPAC complements: R↔Y, K↔M, B↔V, D↔H, S and W fixed, an involution", comp);
    s.require("protein letters round-trip in the order of the NCBI matrices", decode_protein(encode_protein("ARNDCQEGHILKMFPSTWYVBZX*")) == "ARNDCQEGHILKMFPSTWYVBZX*" && protein_code('W') == 17);
    s.require("unknown amino acids (U, O, J) read as X", decode_protein(encode_protein("UOJ")) == "XXX");
  }

  // ---- packed DNA ----------------------------------------------------------------------------------
  {
    bool round = true, index = true, counts = true, rc = true, kmers = true;
    for (std::size_t n : {0, 1, 2, 31, 32, 33, 63, 64, 65, 100, 1000, 4099}) {
      Codes x = random_dna(r, n);
      // Runs of N of lengths 1, 31, 32, 33 where they fit.
      std::size_t at = 3;
      for (std::size_t len : {1, 31, 32, 33})
        if (at + len < n) { std::fill(x.begin() + static_cast<std::ptrdiff_t>(at), x.begin() + static_cast<std::ptrdiff_t>(at + len), N); at += len + 50; }
      const PackedDna p(x);
      round = round && p.codes() == x && p.size() == n;
      for (std::size_t i = 0; i < n; i++) index = index && p[i] == x[i];
      const auto c = p.counts();
      const auto d = base_counts(x);
      counts = counts && std::equal(c.begin(), c.end(), d.begin());
      rc = rc && p.reverse_complement().codes() == reverse_complement(x) && p.reverse_complement().reverse_complement() == p;
      if (n >= 32)
        for (int trial = 0; trial < 50; trial++) {
          const int k = 1 + static_cast<int>(r.below(32));
          const std::size_t i = r.below(n - static_cast<std::size_t>(k) + 1);
          Codes w(x.begin() + static_cast<std::ptrdiff_t>(i), x.begin() + static_cast<std::ptrdiff_t>(i) + k);
          for (auto& b : w) b &= 3;  // N reads as A inside the words
          kmers = kmers && p.kmer(i, k) == encode_kmer(w.data(), k);
        }
    }
    s.require("packed DNA: the codes come back, N runs of 1, 31, 32 and 33 included", round);
    s.require("packed DNA: every base reads back", index);
    s.require("packed DNA: counts by popcount agree with counting", counts);
    s.require("packed DNA: the word-wise reverse complement agrees with the base-wise one, and is an involution", rc);
    s.require("packed DNA: a k-mer at any position, k from 1 to 32, across word boundaries", kmers);
    const PackedDna big(random_dna(r, 1'000'000));
    s.require("packed DNA: a megabase takes a quarter of a byte a base", big.bytes() <= 250'008, std::to_string(big.bytes()) + " bytes");
  }

  // ---- genetic codes ----------------------------------------------------------------------------------
  {
    const auto& all = genetic_codes();
    bool well = all.size() == 27;
    std::set<int> ids;
    for (const auto& g : all) {
      well = well && g.amino_acids.size() == 64 && g.starts.size() == 64 && ids.insert(g.id).second;
      for (char c : g.starts) well = well && (c == 'M' || c == '*' || c == '-');
    }
    s.require("the twenty-seven NCBI translation tables are well formed, each id once", well, std::to_string(all.size()) + " tables");
    const auto& std1 = genetic_code(1);
    const auto aa = [&](const GeneticCode& g, const char* codon) { const Codes c = encode_dna(codon); return amino_acid(g, c[0], c[1], c[2]); };
    const auto start = [&](const GeneticCode& g, const char* codon) { const Codes c = encode_dna(codon); return is_start(g, c[0], c[1], c[2]); };
    s.require("the standard code: TAA, TAG and TGA stop; ATG is M; TGG is W; AGA is R",
              aa(std1, "TAA") == '*' && aa(std1, "TAG") == '*' && aa(std1, "TGA") == '*' && aa(std1, "ATG") == 'M' && aa(std1, "TGG") == 'W' && aa(std1, "AGA") == 'R');
    s.require("the standard code starts at ATG, CTG and TTG", start(std1, "ATG") && start(std1, "CTG") && start(std1, "TTG") && !start(std1, "GTG"));
    const auto& mito = genetic_code(2);
    s.require("the vertebrate mitochondrial code: AGA and AGG stop, TGA is W, ATA is M",
              aa(mito, "AGA") == '*' && aa(mito, "AGG") == '*' && aa(mito, "TGA") == 'W' && aa(mito, "ATA") == 'M');
    const auto& bact = genetic_code(11);
    s.require("the bacterial code starts at ATG, GTG, TTG, CTG, ATT, ATC and ATA",
              start(bact, "ATG") && start(bact, "GTG") && start(bact, "TTG") && start(bact, "CTG") && start(bact, "ATT") && start(bact, "ATC") && start(bact, "ATA"));
    s.require("ambiguous codons: GCN is A and GGN is G in every reading; TTN is X (F or L); NNN is X",
              aa(std1, "GCN") == 'A' && aa(std1, "GGN") == 'G' && aa(std1, "TTN") == 'X' && aa(std1, "NNN") == 'X' && aa(std1, "TTY") == 'X');
    s.require("codon_of inverts codon_index over the 64 codons", [] {
      for (int i = 0; i < 64; i++) { const auto c = codon_of(i); if (codon_index(c[0], c[1], c[2]) != i) return false; }
      return true;
    }());
    // Back-translate a peptide with the first codon of each amino acid, and translate it again.
    const std::string peptide = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWAVKTWKVKGYPDLAI*";
    Codes cds;
    for (char a : peptide) {
      int found = -1;
      for (int i = 0; i < 64 && found < 0; i++)
        if (std1.amino_acids[static_cast<std::size_t>(i)] == a) found = i;
      const auto c = codon_of(found);
      cds.insert(cds.end(), c.begin(), c.end());
    }
    const std::string got = translate(cds, std1, 0, true);
    s.require("a peptide back-translated and translated again (the final stop dropped, as a coding sequence)", got + "*" == peptide, got.substr(0, 20));
    // Frames.
    const Codes f = encode_dna("CATGGCCTAA");
    s.require("translation in frames 0, 1 and 2", translate(f, std1, 0) == "HGL" && translate(f, std1, 1) == "MA*" && translate(f, std1, 2) == "WP");
  }

  // ---- open reading frames ----------------------------------------------------------------------------------
  {
    const auto& g = genetic_code(11);
    bool starts = true, stops = true, inside = true, maximal = true, any = false, both = false;
    std::set<bool> strands;
    for (int trial = 0; trial < 20; trial++) {
      const Codes x = random_dna(r, 3000);
      const auto orfs = open_reading_frames(x, g, 30);
      for (const auto& o : orfs) {
        any = true;
        strands.insert(o.reverse);
        Codes y(x.begin() + static_cast<std::ptrdiff_t>(o.begin), x.begin() + static_cast<std::ptrdiff_t>(o.end));
        if (o.reverse) y = reverse_complement(y);
        starts = starts && y[0] == A && y[1] == T && y[2] == G;
        stops = stops && is_stop(g, y[y.size() - 3], y[y.size() - 2], y[y.size() - 1]) && o.codons() >= 30 && y.size() % 3 == 0;
        for (std::size_t i = 0; i + 3 < y.size(); i += 3) inside = inside && !is_stop(g, y[i], y[i + 1], y[i + 2]);
        // Maximal: no ATG in frame between the previous stop (or the strand's start) and this start.
        const Codes& strand = o.reverse ? reverse_complement(x) : x;
        const std::size_t at = o.reverse ? x.size() - o.end : o.begin;
        for (std::size_t i = at; i >= 3; i -= 3) {
          const Code a = strand[i - 3], b = strand[i - 2], c = strand[i - 1];
          if (is_stop(g, a, b, c)) break;
          if (a == A && b == T && c == G) { maximal = false; break; }
        }
      }
    }
    both = strands.size() == 2;
    s.require("ORFs were found in random sequence, on both strands", any && both);
    s.require("every ORF starts with ATG on its own strand", starts);
    s.require("every ORF ends at its first stop in frame, with at least 30 codons", stops && inside);
    s.require("every ORF begins at the first start after the previous stop", maximal);
  }

  // ---- composition ----------------------------------------------------------------------------------
  {
    const Codes x = encode_dna("GGGCCCAATTNNGC");
    const auto c = base_counts(x);
    s.require("base counts: A 2, C 4, G 4, T 2, N 2", c[A] == 2 && c[C] == 4 && c[G] == 4 && c[T] == 2 && c[N] == 2);
    s.check("GC fraction 8/12, N not counted", gc_fraction(x) - 8.0 / 12, 1e-15);
    const Codes iid = random_dna(r, 400'000);
    s.check("CpG o/e of i.i.d. sequence is near 1", cpg_observed_expected(iid) - 1, 0.03);
    Codes depleted = iid;
    for (std::size_t i = 1; i < depleted.size(); i++)
      if (depleted[i - 1] == C && depleted[i] == G) depleted[i] = A;
    s.check("CpG o/e without a CpG is 0", cpg_observed_expected(depleted), 0);
    s.check("entropy of 2-mers of i.i.d. sequence: 4 bits", kmer_entropy(iid, 2) - 4, 0.001);
    s.check("entropy of a homopolymer: 0 bits", kmer_entropy(Codes(500, A), 3), 0);
    const auto w = gc_windows(iid, 1000, 500);
    s.require("GC windows of 1000 every 500 tile the sequence", w.size() == 799 && w.front().begin == 0 && w.back().end == iid.size());
    // A chromosome whose leading strands are G-rich: G > C after the origin, C > G before it.
    const std::size_t L = 200'000, ori = 60'000, ter = 160'000;
    Codes chr(L);
    for (std::size_t i = 0; i < L; i++) {
      const bool leading = i >= ori && i < ter;
      const double u = r.uniform();
      // G and C at 30% and 20% on one replichore, swapped on the other.
      chr[i] = u < 0.25 ? A : u < 0.5 ? T : u < (leading ? 0.8 : 0.7) ? G : C;
    }
    const auto skew = cumulative_gc_skew(chr, 1000);
    s.require("the cumulative skew's minimum marks the origin and its maximum the terminus (within 2 kb)",
              std::abs(static_cast<long>(skew.minimum) - static_cast<long>(ori)) < 2000 && std::abs(static_cast<long>(skew.maximum) - static_cast<long>(ter)) < 2000,
              std::to_string(skew.minimum) + ", " + std::to_string(skew.maximum));
    const auto& g = genetic_code(1);
    std::array<std::size_t, 64> even{};
    for (auto& v : even) v = 1000;
    const auto rs = rscu(even, g);
    s.require("RSCU of even usage is 1 for every codon", std::all_of(rs.begin(), rs.end(), [](double v) { return std::fabs(v - 1) < 1e-12; }));
    s.check("N_c of even usage is 61", effective_number_of_codons(even, g) - 61, 0.5);
    std::array<std::size_t, 64> one{};
    std::set<char> used;
    for (int i = 0; i < 64; i++)
      if (used.insert(g.amino_acids[static_cast<std::size_t>(i)]).second) one[static_cast<std::size_t>(i)] = 500;
    s.check("N_c with one codon for each amino acid is 20", effective_number_of_codons(one, g) - 20, 1e-9);
  }

  return s.done();
}
