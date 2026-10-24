// The multiple alignment of sequences: its rows and columns, what a column holds (residues,
// gaps, the consensus, the entropy of its residues), sub-alignments, and the two text formats
// it travels in, aligned FASTA and Clustal.
//
// - An Msa is a matrix of codes, one row per sequence and one column per aligned position, with
//   sequence::GAP where a row has no residue; every row has the same length, and a row without
//   its gaps is its sequence. Codes are those of sequence/alphabet.hpp (DNA A, C, G, T, N = 0 to
//   4; protein in the NCBI order ARNDCQEGHILKMFPSTWYVBZX*), so a column indexes a substitution
//   matrix directly.
// - Aligned FASTA (Pearson and Lipman 1988, with '-' or '.' for a gap; the format of MUSCLE,
//   MAFFT and Clustal Omega's --outfmt=fa). Letters are read case-insensitively.
// - Clustal (Thompson, Higgins and Gibson 1994, Nucleic Acids Research 22, 4673–4680; the
//   format ClustalW 1.83 and Clustal Omega write): a first line beginning with CLUSTAL, then
//   blocks of lines "name  residues [count]", the blocks separated by blank lines, each with a
//   line of conservation marks below it. The marks are ClustalW's (alnout.c of ClustalW 1.83;
//   the Clustal Omega documentation): '*' a column of one residue and no gap; ':' a column
//   whose residues all lie in one of the strong groups STA NEQK NHQK NDEQ QHRK MILV MILF HY FYW;
//   '.' one of the weak groups CSA ATV SAG STNK STPA SGND SNDEQK NDEQHK NEQHRK FVLIM HFY;
//   ' ' otherwise. The groups are those of residues whose Gonnet PAM 250 scores exceed 0.5 and
//   0 respectively (the ClustalX help); for DNA only '*' is written.
// - Entropy: the Shannon entropy of the residues of a column, H = −Σ p_a log₂ p_a over the
//   letters below `letters` (ambiguous codes are not counted), at most log₂ letters; its
//   complement log₂ letters − H is the information of the column, the height of a sequence
//   logo's column without the small-sample correction (Schneider and Stephens 1990, Nucleic
//   Acids Research 18, 6097–6100).
//
// Readers never throw on content: rows of different lengths are padded with gaps to the
// longest, and said so, with their lines. Complexity: linear in the size of the matrix.

#pragma once

#include <genomics/sequence/alphabet.hpp>
#include <genomics/sequence/fasta.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::msa {

using sequence::Code;
using sequence::Codes;
using sequence::GAP;
using sequence::Problem;

enum class Alphabet { dna, protein };

// The letters a residue is read as: 4 for DNA (A, C, G, T; N is ambiguous), 20 for protein.
constexpr int letters(Alphabet a) { return a == Alphabet::dna ? 4 : 20; }

struct Msa {
  std::vector<std::string> names;
  std::vector<Codes> rows;  // equal lengths, GAP in the gaps

  std::size_t size() const { return rows.size(); }
  std::size_t columns() const { return rows.empty() ? 0 : rows[0].size(); }
  bool rectangular() const {
    for (const auto& r : rows)
      if (r.size() != columns()) return false;
    return names.size() == rows.size();
  }
  Code at(std::size_t row, std::size_t column) const { return rows[row][column]; }
};

inline void require_rectangular(const Msa& m, const char* who) {
  if (!m.rectangular()) throw std::invalid_argument(std::string(who) + ": the rows are not of one length, or the names do not match the rows");
}

// A row without its gaps: the sequence.
inline Codes ungapped(const Codes& row) {
  Codes out;
  out.reserve(row.size());
  for (Code c : row)
    if (c != GAP) out.push_back(c);
  return out;
}
inline std::vector<Codes> sequences(const Msa& m) {
  std::vector<Codes> out;
  out.reserve(m.size());
  for (const auto& r : m.rows) out.push_back(ungapped(r));
  return out;
}

// An alignment of one sequence: a single row.
inline Msa single(const Codes& sequence, std::string name = "") {
  Msa m;
  m.names.push_back(std::move(name));
  m.rows.push_back(sequence);
  return m;
}

// ---- columns ------------------------------------------------------------------------------------
inline Codes column(const Msa& m, std::size_t c) {
  Codes out(m.size());
  for (std::size_t r = 0; r < m.size(); r++) out[r] = m.rows[r][c];
  return out;
}
inline std::size_t residues_in_column(const Msa& m, std::size_t c) {
  std::size_t n = 0;
  for (const auto& r : m.rows) n += r[c] != GAP;
  return n;
}
inline double gap_fraction(const Msa& m, std::size_t c) {
  return m.size() ? 1.0 - static_cast<double>(residues_in_column(m, c)) / static_cast<double>(m.size()) : 0;
}
inline std::vector<double> gap_fractions(const Msa& m) {
  std::vector<double> out(m.columns());
  for (std::size_t c = 0; c < m.columns(); c++) out[c] = gap_fraction(m, c);
  return out;
}

// The columns [begin, end).
inline Msa slice(const Msa& m, std::size_t begin, std::size_t end) {
  require_rectangular(m, "genomics::msa::slice");
  end = std::min(end, m.columns());
  begin = std::min(begin, end);
  Msa out;
  out.names = m.names;
  for (const auto& r : m.rows) out.rows.emplace_back(r.begin() + static_cast<std::ptrdiff_t>(begin), r.begin() + static_cast<std::ptrdiff_t>(end));
  return out;
}

// The alignment without its columns of gaps only.
inline Msa without_gap_columns(const Msa& m) {
  require_rectangular(m, "genomics::msa::without_gap_columns");
  std::vector<std::size_t> keep;
  for (std::size_t c = 0; c < m.columns(); c++)
    if (residues_in_column(m, c)) keep.push_back(c);
  Msa out;
  out.names = m.names;
  for (const auto& r : m.rows) {
    Codes row(keep.size());
    for (std::size_t k = 0; k < keep.size(); k++) row[k] = r[keep[k]];
    out.rows.push_back(std::move(row));
  }
  return out;
}

// The rows chosen, in the order given, without the columns that are gaps in all of them: the
// alignment the full one induces on a subset of its sequences.
inline Msa sub_alignment(const Msa& m, const std::vector<std::size_t>& chosen) {
  require_rectangular(m, "genomics::msa::sub_alignment");
  Msa out;
  for (std::size_t r : chosen) {
    if (r >= m.size()) throw std::invalid_argument("genomics::msa::sub_alignment: a row outside the alignment");
    out.names.push_back(m.names[r]);
    out.rows.push_back(m.rows[r]);
  }
  return without_gap_columns(out);
}

// The consensus: in each column the most frequent residue below `letters` (ties to the smaller
// code), or GAP where more than max_gap of the rows are gaps or no residue below `letters` is
// present. The aligned row; ungapped() gives the consensus sequence.
inline Codes consensus(const Msa& m, int letters_, double max_gap = 0.5) {
  require_rectangular(m, "genomics::msa::consensus");
  Codes out(m.columns(), GAP);
  std::vector<std::size_t> count(static_cast<std::size_t>(std::max(letters_, 1)), 0);
  for (std::size_t c = 0; c < m.columns(); c++) {
    if (gap_fraction(m, c) > max_gap) continue;
    std::fill(count.begin(), count.end(), 0);
    for (const auto& r : m.rows)
      if (r[c] < letters_) count[r[c]]++;
    std::size_t best = 0;
    for (std::size_t a = 1; a < count.size(); a++)
      if (count[a] > count[best]) best = a;
    if (count[best] > 0) out[c] = static_cast<Code>(best);
  }
  return out;
}

// The Shannon entropy in bits of the residues below `letters` in each column (0 for a column
// with none); with count_gaps, the gap is one more symbol and the bound is log₂(letters + 1).
inline std::vector<double> column_entropy(const Msa& m, int letters_, bool count_gaps = false) {
  require_rectangular(m, "genomics::msa::column_entropy");
  std::vector<double> out(m.columns(), 0);
  std::vector<std::size_t> count(static_cast<std::size_t>(letters_) + 1, 0);
  for (std::size_t c = 0; c < m.columns(); c++) {
    std::fill(count.begin(), count.end(), 0);
    std::size_t total = 0;
    for (const auto& r : m.rows) {
      const Code x = r[c];
      if (x < letters_) { count[x]++; total++; }
      else if (x == GAP && count_gaps) { count[static_cast<std::size_t>(letters_)]++; total++; }
    }
    double h = 0;
    for (std::size_t n : count)
      if (n) {
        const double p = static_cast<double>(n) / static_cast<double>(total);
        h -= p * std::log2(p);
      }
    out[c] = h;
  }
  return out;
}

// log₂ letters − H: the information of each column, 0 for a uniform column and log₂ letters
// for a conserved one.
inline std::vector<double> column_information(const Msa& m, int letters_) {
  auto h = column_entropy(m, letters_);
  const double top = std::log2(static_cast<double>(letters_));
  for (double& x : h) x = top - x;
  return h;
}

// ---- the conservation line of Clustal ------------------------------------------------------------
namespace detail {

inline constexpr std::array<std::string_view, 9> STRONG = {"STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW"};
inline constexpr std::array<std::string_view, 11> WEAK = {"CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK", "NEQHRK", "FVLIM", "HFY"};

// Whether every letter of `present` (a set of letters as a bit mask over 'A'..'Z') lies in one
// of the groups.
template <std::size_t K>
bool within_a_group(std::uint32_t present, const std::array<std::string_view, K>& groups) {
  for (const auto& g : groups) {
    std::uint32_t mask = 0;
    for (char c : g) mask |= 1u << (c - 'A');
    if ((present & ~mask) == 0) return true;
  }
  return false;
}

}  // namespace detail

// One character per column: '*', ':', '.' or ' ' (Clustal's marks; for DNA only '*').
inline std::string clustal_conservation(const Msa& m, Alphabet alphabet) {
  require_rectangular(m, "genomics::msa::clustal_conservation");
  std::string out(m.columns(), ' ');
  for (std::size_t c = 0; c < m.columns(); c++) {
    std::uint32_t present = 0;
    bool gap = false, ambiguous = false;
    for (const auto& r : m.rows) {
      const Code x = r[c];
      if (x == GAP) { gap = true; break; }
      if (x >= letters(alphabet)) { ambiguous = true; continue; }
      const char l = alphabet == Alphabet::dna ? sequence::dna_letter(x) : sequence::protein_letter(x);
      present |= 1u << (l - 'A');
    }
    if (gap || ambiguous || present == 0 || m.size() == 0) continue;
    if ((present & (present - 1)) == 0) out[c] = '*';
    else if (alphabet == Alphabet::protein && detail::within_a_group(present, detail::STRONG)) out[c] = ':';
    else if (alphabet == Alphabet::protein && detail::within_a_group(present, detail::WEAK)) out[c] = '.';
  }
  return out;
}

// ---- letters ------------------------------------------------------------------------------------
inline char letter(Code x, Alphabet a) {
  if (x == GAP) return '-';
  return a == Alphabet::dna ? sequence::dna_letter(x) : sequence::protein_letter(x);
}
inline Code code(char c, Alphabet a) {
  if (c == '-' || c == '.') return GAP;
  return a == Alphabet::dna ? sequence::dna_code(c) : sequence::protein_code(c);
}
inline std::string row_string(const Codes& row, Alphabet a) {
  std::string s(row.size(), '-');
  for (std::size_t i = 0; i < row.size(); i++) s[i] = letter(row[i], a);
  return s;
}

namespace detail {

template <class F>
void lines(std::string_view text, F&& f) {
  std::size_t number = 0, i = 0;
  while (i < text.size()) {
    std::size_t j = text.find('\n', i);
    if (j == std::string_view::npos) j = text.size();
    std::string_view line = text.substr(i, j - i);
    if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
    f(line, ++number);
    i = j + 1;
  }
}

// Rows of unequal length padded with gaps to the longest, each reported.
inline void square(Msa& m, const std::vector<std::size_t>& line_of, std::vector<Problem>* problems) {
  std::size_t width = 0;
  for (const auto& r : m.rows) width = std::max(width, r.size());
  for (std::size_t k = 0; k < m.rows.size(); k++)
    if (m.rows[k].size() < width) {
      if (problems)
        problems->push_back({line_of[k], "row " + m.names[k] + " has " + std::to_string(m.rows[k].size()) + " columns where the longest has " +
                                             std::to_string(width) + ": padded with gaps"});
      m.rows[k].resize(width, GAP);
    }
}

}  // namespace detail

// ---- aligned FASTA ------------------------------------------------------------------------------------
inline Msa read_aligned_fasta(std::string_view text, Alphabet alphabet, std::vector<Problem>* problems = nullptr) {
  Msa m;
  std::vector<std::size_t> line_of;
  bool warned = false;
  detail::lines(text, [&](std::string_view line, std::size_t number) {
    if (line.empty() || line.front() == ';') return;
    if (line.front() == '>') {
      std::string_view h = line.substr(1);
      while (!h.empty() && (h.front() == ' ' || h.front() == '\t')) h.remove_prefix(1);
      const std::size_t sp = h.find_first_of(" \t");
      m.names.emplace_back(h.substr(0, sp));
      m.rows.emplace_back();
      line_of.push_back(number);
      if (m.names.back().empty() && problems) problems->push_back({number, "a record without a name"});
      return;
    }
    if (m.rows.empty()) {
      if (!warned && problems) problems->push_back({number, "residues before the first '>' header: skipped"});
      warned = true;
      return;
    }
    for (char c : line)
      if (c != ' ' && c != '\t') m.rows.back().push_back(code(c, alphabet));
  });
  detail::square(m, line_of, problems);
  return m;
}

inline std::string write_aligned_fasta(const Msa& m, Alphabet alphabet, std::size_t width = 60) {
  require_rectangular(m, "genomics::msa::write_aligned_fasta");
  if (width == 0) width = SIZE_MAX;
  std::string out;
  for (std::size_t r = 0; r < m.size(); r++) {
    out += '>' + m.names[r] + '\n';
    const std::string s = row_string(m.rows[r], alphabet);
    for (std::size_t i = 0; i < s.size(); i += width) out += s.substr(i, width) + '\n';
    if (s.empty()) out += '\n';
  }
  return out;
}

// ---- Clustal ------------------------------------------------------------------------------------
// Reads the blocks of a Clustal file: the first line must begin with CLUSTAL (or MUSCLE, which
// writes the same format), lines beginning with white space are conservation lines and are
// skipped, and a trailing residue count on a line is ignored. Rows keep the order in which
// their names first appear.
inline Msa read_clustal(std::string_view text, Alphabet alphabet, std::vector<Problem>* problems = nullptr) {
  Msa m;
  std::vector<std::size_t> line_of;
  bool header = false;
  detail::lines(text, [&](std::string_view line, std::size_t number) {
    if (!header) {
      if (line.empty()) return;
      header = true;
      if (line.substr(0, 7) != "CLUSTAL" && line.substr(0, 6) != "MUSCLE" && problems)
        problems->push_back({number, "the first line does not begin with CLUSTAL"});
      if (line.substr(0, 7) == "CLUSTAL" || line.substr(0, 6) == "MUSCLE") return;
    }
    if (line.empty() || line.front() == ' ' || line.front() == '\t') return;
    const std::size_t sp = line.find_first_of(" \t");
    if (sp == std::string_view::npos) {
      if (problems) problems->push_back({number, "a line with a name and no residues"});
      return;
    }
    const std::string name(line.substr(0, sp));
    std::string_view rest = line.substr(sp);
    while (!rest.empty() && (rest.front() == ' ' || rest.front() == '\t')) rest.remove_prefix(1);
    const std::size_t end = rest.find_first_of(" \t");
    std::string_view residues = rest.substr(0, end);
    std::size_t k = 0;
    while (k < m.names.size() && m.names[k] != name) k++;
    if (k == m.names.size()) {
      m.names.push_back(name);
      m.rows.emplace_back();
      line_of.push_back(number);
    }
    for (char c : residues) m.rows[k].push_back(code(c, alphabet));
  });
  if (!header && problems) problems->push_back({0, "an empty file"});
  detail::square(m, line_of, problems);
  return m;
}

inline std::string write_clustal(const Msa& m, Alphabet alphabet, std::size_t width = 60) {
  require_rectangular(m, "genomics::msa::write_clustal");
  if (width == 0) width = 60;
  std::size_t name_width = 10;
  for (const auto& n : m.names) name_width = std::max(name_width, n.size() + 1);
  name_width += 5;
  std::string out = "CLUSTAL multiple sequence alignment (genomics::msa)\n\n";
  const std::string marks = clustal_conservation(m, alphabet);
  std::vector<std::size_t> count(m.size(), 0);
  for (std::size_t b = 0; b < m.columns(); b += width) {
    const std::size_t e = std::min(m.columns(), b + width);
    out += '\n';
    for (std::size_t r = 0; r < m.size(); r++) {
      std::string line = m.names[r];
      line.resize(name_width, ' ');
      for (std::size_t c = b; c < e; c++) {
        line += letter(m.rows[r][c], alphabet);
        count[r] += m.rows[r][c] != GAP;
      }
      out += line + ' ' + std::to_string(count[r]) + '\n';
    }
    out += std::string(name_width, ' ') + marks.substr(b, e - b) + '\n';
  }
  return out;
}

}  // namespace genomics::msa
