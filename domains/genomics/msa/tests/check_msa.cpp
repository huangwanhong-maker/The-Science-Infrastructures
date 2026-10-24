// Checks of the multiple alignment type and its formats: aligned FASTA and Clustal written and
// read back exactly, for DNA and protein, with rows of unequal length padded and reported at
// their lines; a file in the form Clustal Omega writes read whole; ClustalW's conservation
// marks on columns built to fall in its strong and weak groups; columns, gap fractions, slices
// and sub-alignments; the consensus; the entropy of columns against hand-worked values (a
// uniform column of four bases 2 bits, two bases 1 bit) and the information as its complement.

#include <genomics/msa/msa.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace genomics::msa;
namespace seq = genomics::sequence;

namespace {

Msa from_letters(const std::vector<std::string>& rows, Alphabet a) {
  Msa m;
  for (std::size_t r = 0; r < rows.size(); r++) {
    m.names.push_back("r" + std::to_string(r + 1));
    Codes row;
    for (char c : rows[r]) row.push_back(code(c, a));
    m.rows.push_back(row);
  }
  return m;
}

Msa random_msa(calx::random::Stream& r, std::size_t rows, std::size_t cols, Alphabet a, double gaps) {
  Msa m;
  for (std::size_t i = 0; i < rows; i++) {
    m.names.push_back("seq_" + std::to_string(i));
    Codes row(cols);
    for (auto& c : row) c = r.uniform() < gaps ? GAP : static_cast<Code>(r.below(a == Alphabet::dna ? 5 : 24));
    m.rows.push_back(row);
  }
  return m;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x4d534131, 0, 0);

  // ---- formats ------------------------------------------------------------------------------------
  {
    bool fasta = true, clustal = true;
    for (int t = 0; t < 40; t++) {
      const Alphabet a = t % 2 ? Alphabet::protein : Alphabet::dna;
      const Msa m = random_msa(r, 1 + r.below(9), r.below(200), a, 0.3);
      std::vector<Problem> p1, p2;
      const Msa f = read_aligned_fasta(write_aligned_fasta(m, a, 1 + r.below(80)), a, &p1);
      const Msa c = read_clustal(write_clustal(m, a, 1 + r.below(80)), a, &p2);
      fasta = fasta && f.names == m.names && f.rows == m.rows && p1.empty();
      // A Clustal row with no columns writes no residues, and reads back as nothing.
      if (m.columns() > 0) clustal = clustal && c.names == m.names && c.rows == m.rows && p2.empty();
    }
    s.require("aligned FASTA written and read back gives the alignment exactly (40 random DNA and protein alignments)", fasta);
    s.require("Clustal written and read back gives the alignment exactly (40 random DNA and protein alignments)", clustal);

    std::vector<Problem> p;
    const Msa ragged = read_aligned_fasta(">a\nAC-GT\n>b\nACG\n>c\nA-\nCGTT\n", Alphabet::dna, &p);
    s.require("a short row of aligned FASTA is padded with gaps to the longest and reported at its header's line",
              ragged.rectangular() && ragged.columns() == 6 && ragged.rows[1] == Codes{seq::A, seq::C, seq::G, GAP, GAP, GAP} && p.size() == 2 && p[0].line == 1 &&
                  p[1].line == 3,
              std::to_string(p.size()) + " problems");

    // Human α and β globin's first residues (UniProt P69905, P68871), aligned by hand and written in
    // the form Clustal Omega writes: a header, blocks with running counts, a line of marks below each.
    const std::string clustalo =
        "CLUSTAL O(1.2.4) multiple sequence alignment\n"
        "\n"
        "\n"
        "sp|P69905|HBA   MV-LSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF 47\n"
        "sp|P68871|HBB   MVHLTPEEKSAVTALWGKVN--VDEVGGEALGRLLVVYPWTQRFFESF 46\n"
        "                **:*:*.:*. * * ****.     * *.*** *:: .:* *:.:* .*\n"
        "\n"
        "sp|P69905|HBA   DLSH 51\n"
        "sp|P68871|HBB   GDLS 50\n"
        "                 ** \n";
    std::vector<Problem> pc;
    const Msa hb = read_clustal(clustalo, Alphabet::protein, &pc);
    s.require("a file in Clustal Omega's form is read whole: two rows of 52 columns, names as written, gaps where it has them",
              hb.size() == 2 && hb.columns() == 52 && hb.names[0] == "sp|P69905|HBA" && hb.rows[0][2] == GAP && hb.rows[1][20] == GAP && hb.rows[1][21] == GAP &&
                  ungapped(hb.rows[0]).size() == 51 && pc.empty());
    std::vector<Problem> pn;
    const Msa bare = read_clustal("x ACGT\ny AC-T\n", Alphabet::dna, &pn);
    s.require("a Clustal file without its header is read whole and reported at line 1", bare.size() == 2 && bare.columns() == 4 && pn.size() == 1 && pn[0].line == 1);
  }

  // ---- the conservation line ------------------------------------------------------------------------------------
  {
    // Columns: identical; S,T,A (strong STA); C,S,A (weak CSA only); W,L (no group); a gap; M,I,L,V (strong); H,F,Y (weak HFY).
    const Msa m = from_letters({"ASCWAMH", "ATSL-IF", "AAAWAVY", "ASCLAL" "H"}, Alphabet::protein);
    const std::string marks = clustal_conservation(m, Alphabet::protein);
    s.require("Clustal's marks: '*' for one residue, ':' for a strong group, '.' for a weak one, ' ' for none or a gap", marks == "*:.  :.", "'" + marks + "'");
    const Msa d = from_letters({"ACGT-", "ACGA-", "ACGTT"}, Alphabet::dna);
    s.require("for DNA only '*' is written", clustal_conservation(d, Alphabet::dna) == "***  ");
  }

  // ---- columns, slices, sub-alignments, consensus ------------------------------------------------------------------------------------
  {
    const Msa m = from_letters({"AC-GT-", "AC-GA-", "-T-GAC", "ACAG--"}, Alphabet::dna);
    s.require("a column read down the rows", column(m, 1) == Codes{seq::C, seq::C, seq::T, seq::C});
    s.check("the gap fraction of a column (2 of 4 rows) is one half", gap_fraction(m, 5) - 0.75, 1e-15);
    const Msa w = without_gap_columns(slice(m, 0, 3));
    s.require("a slice keeps its columns, and removing columns of gaps only keeps the rest", w.columns() == 3 && slice(m, 2, 3).columns() == 1);
    const Msa sub = sub_alignment(m, {0, 1});
    s.require("a sub-alignment of rows 1 and 2 drops the columns that are gaps in both", sub.columns() == 4 && row_string(sub.rows[0], Alphabet::dna) == "ACGT" &&
                                                                                           row_string(sub.rows[1], Alphabet::dna) == "ACGA");
    s.require("the sequences of a sub-alignment are those of its rows", sequences(sub) == std::vector<Codes>{ungapped(m.rows[0]), ungapped(m.rows[1])});
    const Codes cons = consensus(m, 4, 0.5);
    s.require("the consensus takes each column's most frequent base (ties to the smaller code) and leaves columns of more than half gaps empty",
              row_string(cons, Alphabet::dna) == "AC-GA-", row_string(cons, Alphabet::dna));
  }

  // ---- entropy ------------------------------------------------------------------------------------
  {
    const Msa m = from_letters({"AAAA", "ACCA", "AGA-", "ATC-"}, Alphabet::dna);
    const auto h = column_entropy(m, 4);
    const auto hg = column_entropy(m, 4, true);
    const auto info = column_information(m, 4);
    s.check("a conserved column has entropy 0 bits", h[0], 1e-15, "bits");
    s.check("a column of A, C, G, T has entropy 2 bits (log2 4)", h[1] - 2.0, 1e-15, "bits");
    s.check("a column of two bases, twice each, has entropy 1 bit", h[2] - 1.0, 1e-15, "bits");
    s.check("a column of A, A and two gaps has entropy 0 bits without gaps and 1 bit with the gap a symbol", std::fabs(h[3]) + std::fabs(hg[3] - 1.0), 1e-15, "bits");
    s.check("the information of a column is log2 4 minus its entropy", std::fabs(info[0] - 2.0) + std::fabs(info[1]), 1e-15, "bits");
    // Protein: twenty residues once each.
    Msa p;
    for (int a = 0; a < 20; a++) {
      p.names.push_back("p" + std::to_string(a));
      p.rows.push_back(Codes{static_cast<Code>(a)});
    }
    s.check("a protein column of the twenty residues once each has entropy log2 20", column_entropy(p, 20)[0] - std::log2(20.0), 1e-13, "bits");
  }
  return s.done();
}
