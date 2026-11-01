// The classes of substitutions in their sequence context, as mutational catalogues count them.
//
// A substitution of one base of a double-stranded genome is the same event read on either
// strand: C>T at 5′-ACG-3′ on one strand is G>A at 5′-CGT-3′ on the other. The classification of
// Alexandrov et al. (2013, "Signatures of mutational processes in human cancer", Nature 500,
// 415–421, Methods; COSMIC's SBS-96) names each event once, on the strand where the mutated base
// is a pyrimidine (C or T): six substitution types C>A, C>G, C>T, T>A, T>C, T>G, each with the
// base 5′ and the base 3′ of the mutated one, 6 · 4 · 4 = 96 classes. In the language of
// grammar.hpp, the strand group C₂ acts on the 192 directed trinucleotide substitutions freely
// (no trinucleotide is its own reverse complement, the middle base being exchanged with its
// complement), and the pyrimidine-centred class is the orbit's representative.
//
// - Index: k = 16 s + 4 l + r, s the substitution type in the order above, l and r the 5′ and 3′
//   bases in the order A, C, G, T: the order in which the profiles are drawn, by type. Names as
//   COSMIC and SigProfilerMatrixGenerator write them, "A[C>A]A"; SigProfilerMatrixGenerator's
//   matrices list the classes in the lexicographic order of their names (5′ base first), which
//   sbs96_lexicographic gives. Tables are read by name, so either order is accepted.
// - Doublets (DBS-78, Alexandrov et al. 2020, "The repertoire of mutational signatures in human
//   cancer", Nature 578, 94–101, Methods): two adjacent bases both substituted, 16 · 9 = 144
//   directed events identified with their reverse complements. The four dinucleotides equal to
//   their reverse complement (AT, CG, GC, TA) have orbits of their nine alternates of sizes 1 and
//   2 (six orbits each); the other twelve pair off into six classes of nine; 4 · 6 + 6 · 9 = 78.
//   Which member names an orbit is a convention: the list of SigProfilerMatrixGenerator 1.4.0
//   (its mutation_types and mutation_types_non_tsb), in the lexicographic order of its DBS78
//   matrices (pandas' sorted groupby), held to that source in the checks.
//
// Everything here is a pure function of codes; N and other codes give no class (−1).

#pragma once

#include <genomics/sequence/alphabet.hpp>

#include <array>
#include <cstddef>
#include <string>
#include <string_view>

namespace genomics::sequence {

// ---- SBS-96 ----------------------------------------------------------------------------------------

inline constexpr std::size_t SBS96 = 96;

// The six substitution types, the reference base a pyrimidine: (from, to) codes.
inline constexpr std::array<std::array<Code, 2>, 6> SBS_TYPES{{{C, A}, {C, G}, {C, T}, {T, A}, {T, C}, {T, G}}};

struct SubstitutionContext {
  Code left = N, ref = N, alt = N, right = N;  // 5′ base, reference base, alternate, 3′ base, on the pyrimidine strand
};

// The class of the substitution ref>alt between the bases left (5′) and right (3′) on the strand
// given, or −1 when a code is not a base or alt = ref.
constexpr int sbs96_index(Code left, Code ref, Code alt, Code right) {
  if (left > 3 || ref > 3 || alt > 3 || right > 3 || ref == alt) return -1;
  if (ref == A || ref == G) {  // a purine: read the other strand
    const Code l = complement(right), r = complement(left);
    left = l;
    right = r;
    ref = complement(ref);
    alt = complement(alt);
  }
  int s = -1;
  for (int t = 0; t < 6; t++)
    if (SBS_TYPES[static_cast<std::size_t>(t)][0] == ref && SBS_TYPES[static_cast<std::size_t>(t)][1] == alt) s = t;
  return 16 * s + 4 * left + right;
}

constexpr SubstitutionContext sbs96_context(std::size_t k) {
  if (k >= SBS96) return {};
  const auto& t = SBS_TYPES[k / 16];
  return {static_cast<Code>((k % 16) / 4), t[0], t[1], static_cast<Code>(k % 4)};
}

// "A[C>A]A".
inline std::string sbs96_name(std::size_t k) {
  const SubstitutionContext c = sbs96_context(k);
  std::string s;
  s += dna_letter(c.left);
  s += '[';
  s += dna_letter(c.ref);
  s += '>';
  s += dna_letter(c.alt);
  s += ']';
  s += dna_letter(c.right);
  return s;
}

// The class of a substitution at position p of x (0 < p < |x| − 1), or −1.
inline int sbs96_at(const Codes& x, std::size_t p, Code alt) {
  if (p == 0 || p + 1 >= x.size()) return -1;
  return sbs96_index(x[p - 1], x[p], alt, x[p + 1]);
}

// "C>A".
inline std::string sbs_type_name(std::size_t s) {
  if (s >= 6) return "?";
  return std::string{dna_letter(SBS_TYPES[s][0]), '>', dna_letter(SBS_TYPES[s][1])};
}

// The class at position i of the lexicographic order of the names (SigProfilerMatrixGenerator's):
// 5′ base, then the substitution, then the 3′ base.
constexpr std::size_t sbs96_lexicographic(std::size_t i) {
  if (i >= SBS96) return SBS96;
  const std::size_t l = i / 24, s = (i % 24) / 4, r = i % 4;
  return 16 * s + 4 * l + r;
}

// ---- DBS-78 ----------------------------------------------------------------------------------------

inline constexpr std::size_t DBS78 = 78;

inline constexpr std::array<std::string_view, DBS78> DBS78_NAMES{
    "AC>CA", "AC>CG", "AC>CT", "AC>GA", "AC>GG", "AC>GT", "AC>TA", "AC>TG", "AC>TT", "AT>CA", "AT>CC", "AT>CG", "AT>GA", "AT>GC", "AT>TA", "CC>AA",
    "CC>AG", "CC>AT", "CC>GA", "CC>GG", "CC>GT", "CC>TA", "CC>TG", "CC>TT", "CG>AT", "CG>GC", "CG>GT", "CG>TA", "CG>TC", "CG>TT", "CT>AA", "CT>AC",
    "CT>AG", "CT>GA", "CT>GC", "CT>GG", "CT>TA", "CT>TC", "CT>TG", "GC>AA", "GC>AG", "GC>AT", "GC>CA", "GC>CG", "GC>TA", "TA>AT", "TA>CG", "TA>CT",
    "TA>GC", "TA>GG", "TA>GT", "TC>AA", "TC>AG", "TC>AT", "TC>CA", "TC>CG", "TC>CT", "TC>GA", "TC>GG", "TC>GT", "TG>AA", "TG>AC", "TG>AT", "TG>CA",
    "TG>CC", "TG>CT", "TG>GA", "TG>GC", "TG>GT", "TT>AA", "TT>AC", "TT>AG", "TT>CA", "TT>CC", "TT>CG", "TT>GA", "TT>GC", "TT>GG"};

namespace detail {

// The class of each directed doublet (r₀r₁ > a₀a₁ as 64 r₀ + 16 r₁ + 4 a₀ + a₁) named in the list, else −1.
constexpr std::array<int, 256> dbs78_direct() {
  std::array<int, 256> t{};
  for (auto& x : t) x = -1;
  for (std::size_t k = 0; k < DBS78; k++) {
    const auto n = DBS78_NAMES[k];
    const int key = 64 * dna_code(n[0]) + 16 * dna_code(n[1]) + 4 * dna_code(n[3]) + dna_code(n[4]);
    t[static_cast<std::size_t>(key)] = static_cast<int>(k);
  }
  return t;
}
inline constexpr auto DBS78_DIRECT = dbs78_direct();

}  // namespace detail

// The class of the doublet r₀r₁ > a₀a₁ (both bases changed) read on either strand, or −1.
constexpr int dbs78_index(Code r0, Code r1, Code a0, Code a1) {
  if (r0 > 3 || r1 > 3 || a0 > 3 || a1 > 3 || r0 == a0 || r1 == a1) return -1;
  const int k = detail::DBS78_DIRECT[static_cast<std::size_t>(64 * r0 + 16 * r1 + 4 * a0 + a1)];
  if (k >= 0) return k;
  return detail::DBS78_DIRECT[static_cast<std::size_t>(64 * complement(r1) + 16 * complement(r0) + 4 * complement(a1) + complement(a0))];
}

inline std::string dbs78_name(std::size_t k) { return k < DBS78 ? std::string(DBS78_NAMES[k]) : std::string("?"); }

}  // namespace genomics::sequence
