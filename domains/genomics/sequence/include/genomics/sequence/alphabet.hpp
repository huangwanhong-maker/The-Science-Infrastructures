// The alphabets of the genome and its products, as small integer codes.
//
// Every module of the domain computes on codes, never on letters: a sequence is a vector of
// bytes, read from letters at the boundary (a file, a request) and written back to letters
// there. The codes are fixed once, here:
//
//   DNA (and RNA)   A = 0, C = 1, G = 2, T = 3 (U reads as T), N = 4 for every other letter;
//                   the complement of x < 4 is 3 − x, and that of N is N. With this order the
//                   two bits of a base pack four to a byte, the complement of a packed word is
//                   its bitwise negation, and A < C < G < T is the lexicographic order of the
//                   letters, so that sorting codes sorts sequences as text.
//   IUPAC           a set of bases as a 4-bit mask, A = 1, C = 2, G = 4, T = 8: R = A|G,
//                   Y = C|T, S = C|G, W = A|T, K = G|T, M = A|C, B = C|G|T, D = A|G|T,
//                   H = A|C|T, V = A|C|G, N = A|C|G|T (Cornish-Bowden 1985, Nucleic Acids
//                   Research 13, 3021–3030; the NC-IUB recommendations).
//   protein         the order of the NCBI substitution matrices, ARNDCQEGHILKMFPSTWYVBZX*
//                   (0 to 23): the twenty amino acids, B (D or N), Z (E or Q), X (any) and the
//                   stop *; U (selenocysteine), O (pyrrolysine) and J (I or L) read as X.
//   aligned rows    GAP = 0xFF stands for a gap in either alphabet.
//
// Reading letters never fails: a letter outside the alphabet reads as N (DNA) or X
// (protein), and a caller that must know counts them with is_dna_letter / is_protein_letter.

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace genomics::sequence {

using Code = std::uint8_t;
using Codes = std::vector<Code>;

inline constexpr Code A = 0, C = 1, G = 2, T = 3, N = 4;
inline constexpr Code GAP = 0xFF;
inline constexpr std::string_view DNA_LETTERS = "ACGTN";
inline constexpr std::string_view PROTEIN_LETTERS = "ARNDCQEGHILKMFPSTWYVBZX*";
inline constexpr Code PROTEIN_X = 22, PROTEIN_STOP = 23;

namespace detail {

constexpr std::array<Code, 256> dna_table() {
  std::array<Code, 256> t{};
  for (auto& x : t) x = N;
  t['A'] = t['a'] = A;
  t['C'] = t['c'] = C;
  t['G'] = t['g'] = G;
  t['T'] = t['t'] = T;
  t['U'] = t['u'] = T;
  return t;
}
constexpr std::array<Code, 256> protein_table() {
  std::array<Code, 256> t{};
  for (auto& x : t) x = PROTEIN_X;
  for (std::size_t i = 0; i < PROTEIN_LETTERS.size(); i++) {
    const char c = PROTEIN_LETTERS[i];
    t[static_cast<unsigned char>(c)] = static_cast<Code>(i);
    if (c >= 'A' && c <= 'Z') t[static_cast<unsigned char>(c - 'A' + 'a')] = static_cast<Code>(i);
  }
  return t;
}
constexpr std::array<std::uint8_t, 256> iupac_table() {
  std::array<std::uint8_t, 256> t{};
  const char* letters = "ACMGRSVTWYHKDBN";  // the masks 1 to 15, in order
  for (int m = 1; m <= 15; m++) {
    const char c = letters[m - 1];
    t[static_cast<unsigned char>(c)] = static_cast<std::uint8_t>(m);
    t[static_cast<unsigned char>(c - 'A' + 'a')] = static_cast<std::uint8_t>(m);
  }
  t['U'] = t['u'] = 8;
  return t;
}
inline constexpr auto DNA_TABLE = dna_table();
inline constexpr auto PROTEIN_TABLE = protein_table();
inline constexpr auto IUPAC_TABLE = iupac_table();

}  // namespace detail

// ---- DNA ------------------------------------------------------------------------------------
constexpr Code dna_code(char c) { return detail::DNA_TABLE[static_cast<unsigned char>(c)]; }
constexpr char dna_letter(Code x) { return x < 5 ? DNA_LETTERS[x] : x == GAP ? '-' : 'N'; }
constexpr bool is_dna_letter(char c) { return detail::DNA_TABLE[static_cast<unsigned char>(c)] != N || c == 'N' || c == 'n'; }
constexpr Code complement(Code x) { return x < 4 ? static_cast<Code>(3 - x) : x; }

inline Codes encode_dna(std::string_view s) {
  Codes out(s.size());
  for (std::size_t i = 0; i < s.size(); i++) out[i] = dna_code(s[i]);
  return out;
}
inline std::string decode_dna(const Codes& x) {
  std::string s(x.size(), 'N');
  for (std::size_t i = 0; i < x.size(); i++) s[i] = dna_letter(x[i]);
  return s;
}
inline Codes reverse_complement(const Codes& x) {
  Codes out(x.size());
  for (std::size_t i = 0; i < x.size(); i++) out[x.size() - 1 - i] = complement(x[i]);
  return out;
}
inline std::string reverse_complement(std::string_view s) {
  static constexpr auto table = [] {
    std::array<char, 256> t{};
    for (int i = 0; i < 256; i++) t[i] = static_cast<char>(i);
    const char* from = "ACGTUNacgtunRYKMBVDHrykmbvdhSWswX";
    const char* to = "TGCAANtgcaanYRMKVBHDyrmkvbhdSWswX";
    for (int i = 0; from[i]; i++) t[static_cast<unsigned char>(from[i])] = to[i];
    return t;
  }();
  std::string out(s.size(), 'N');
  for (std::size_t i = 0; i < s.size(); i++) out[s.size() - 1 - i] = table[static_cast<unsigned char>(s[i])];
  return out;
}
// The number of positions that are not A, C, G or T.
inline std::size_t count_ambiguous(const Codes& x) {
  std::size_t n = 0;
  for (Code c : x) n += c > 3;
  return n;
}

// ---- IUPAC ------------------------------------------------------------------------------------
// The set of bases a letter stands for (0 for a letter that is not IUPAC).
constexpr std::uint8_t iupac_mask(char c) { return detail::IUPAC_TABLE[static_cast<unsigned char>(c)]; }
constexpr char iupac_letter(std::uint8_t mask) { return mask >= 1 && mask <= 15 ? "ACMGRSVTWYHKDBN"[mask - 1] : '-'; }
// The complement of a set: A↔T and C↔G exchanged, bit by bit.
constexpr std::uint8_t iupac_complement(std::uint8_t mask) {
  return static_cast<std::uint8_t>(((mask & 1) << 3) | ((mask & 2) << 1) | ((mask & 4) >> 1) | ((mask & 8) >> 3));
}
// Whether a base (a code below 4) belongs to the set a letter stands for.
constexpr bool iupac_matches(char letter, Code base) { return base < 4 && (iupac_mask(letter) >> base) & 1; }

// ---- protein ------------------------------------------------------------------------------------
constexpr Code protein_code(char c) { return detail::PROTEIN_TABLE[static_cast<unsigned char>(c)]; }
constexpr char protein_letter(Code x) { return x < PROTEIN_LETTERS.size() ? PROTEIN_LETTERS[x] : x == GAP ? '-' : 'X'; }
constexpr bool is_protein_letter(char c) { return protein_code(c) != PROTEIN_X || c == 'X' || c == 'x'; }

inline Codes encode_protein(std::string_view s) {
  Codes out(s.size());
  for (std::size_t i = 0; i < s.size(); i++) out[i] = protein_code(s[i]);
  return out;
}
inline std::string decode_protein(const Codes& x) {
  std::string s(x.size(), 'X');
  for (std::size_t i = 0; i < x.size(); i++) s[i] = protein_letter(x[i]);
  return s;
}

}  // namespace genomics::sequence
