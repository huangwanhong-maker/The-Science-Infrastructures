// SPDI, NCBI's data model of a variant (Holmes, Moyer, Phan, Maglott and Kattman 2020, "SPDI: data
// model for variants and applications at NCBI", Bioinformatics 36, 1902–1907): a Sequence, a
// Position, a Deletion and an Insertion, "start at the boundary before the first nucleotide in the
// sequence S, advance P nucleotides, delete D nucleotides, then insert the nucleotides in the string
// I" (§"Variant Model"). P is a 0-based interbase coordinate; D is a length or, redundantly, the
// deleted sequence; I may be empty. `NC_000001.11:12345:A:G`, `Seq1:4:1:G`, `Seq1:5::C`.
//
// - VOCA, the Variant Overprecision Correction Algorithm (§"Correcting Overprecision"), written here
//   as the paper gives it: shrinking trims the deleted and inserted sequences' common 3′ ends, then
//   their 5′ ends (moving P); if both remain non-empty the variant is final. Otherwise the non-empty
//   one is the bud; the blossom repeats it from P in both directions, blossom(i) = bud[(i − P) mod |bud|];
//   the largest interval [a, b) of the reference that matches the blossom and contains P becomes the
//   deleted sequence, and the shrunken variant applied to it the inserted. The paper's Table 2: in
//   ATGAGACT a deleted GA becomes 2:GAGA:GA; GA inserted into ATGACT becomes 2:GA:GAGA.
// - Theorem 2 of theory/nomenclature.md: VOCA's interval is the hull of the change's placements, so
//   voca() equals describe.hpp's spdi_form(); the checks hold the two equal by brute force.
// - Conversions: to a Change (describe.hpp), to VCF fields (POS 1-based; a padding base before an
//   empty allele, VCF 4.3 §1.6.1), and to the HGVS description.
//
// Reading never throws on content: parse_spdi returns what is wrong.

#pragma once

#include <genomics/nomenclature/describe.hpp>

#include <charconv>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>

namespace genomics::nomenclature {

struct SpdiText {
  bool ok = false;
  Spdi spdi;
  bool deletion_as_length = false;  // D was written as a number
  std::size_t deleted_length = 0;
  std::string error;
};

inline SpdiText parse_spdi(std::string_view text) {
  SpdiText out;
  std::string_view f[4];
  std::size_t k = 0, i = 0;
  for (; k < 3; k++) {
    const std::size_t j = text.find(':', i);
    if (j == std::string_view::npos) { out.error = "an SPDI has four fields separated by ':' (sequence:position:deletion:insertion)"; return out; }
    f[k] = text.substr(i, j - i);
    i = j + 1;
  }
  f[3] = text.substr(i);
  if (f[3].find(':') != std::string_view::npos) { out.error = "an SPDI has four fields, not more"; return out; }
  if (f[0].empty()) { out.error = "the sequence is empty"; return out; }
  out.spdi.sequence = std::string(f[0]);
  std::size_t p = 0;
  if (f[1].empty() || std::from_chars(f[1].data(), f[1].data() + f[1].size(), p).ptr != f[1].data() + f[1].size()) {
    out.error = "the position is not a non-negative integer";
    return out;
  }
  out.spdi.position = p;
  const auto bases = [](std::string_view s) {
    for (char c : s)
      if (!sequence::is_dna_letter(c)) return false;
    return true;
  };
  if (!f[2].empty() && f[2].find_first_not_of("0123456789") == std::string_view::npos) {
    out.deletion_as_length = true;
    std::from_chars(f[2].data(), f[2].data() + f[2].size(), out.deleted_length);
  } else if (!bases(f[2])) {
    out.error = "the deletion is neither a length nor a sequence of bases";
    return out;
  } else {
    out.spdi.deleted = detail::upper(f[2]);
    out.deleted_length = f[2].size();
  }
  if (!bases(f[3])) { out.error = "the insertion is not a sequence of bases"; return out; }
  out.spdi.inserted = detail::upper(f[3]);
  out.ok = true;
  return out;
}

// The deleted sequence of an SPDI written with a length, read from the reference.
inline Spdi resolve(const SpdiText& t, const Codes& reference) {
  Spdi s = t.spdi;
  if (t.deletion_as_length) {
    if (s.position + t.deleted_length > reference.size()) throw std::invalid_argument("genomics::nomenclature::resolve: the deletion runs past the sequence");
    s.deleted = detail::slice(reference, s.position, s.position + t.deleted_length);
  }
  return s;
}

inline Change change_of(const Spdi& s) { return Change{s.position, s.deleted, s.inserted}; }

// VOCA as Holmes et al. state it: shrink, then grow by the blossom of the bud.
inline Spdi voca(const Codes& reference, Spdi s) {
  if (!deleted_agrees(reference, change_of(s))) throw std::invalid_argument("genomics::nomenclature::voca: the deleted bases are not the reference's");
  s.deleted = detail::upper(s.deleted);
  s.inserted = detail::upper(s.inserted);
  // Shrinking: the 3′ ends, then the 5′ ends.
  while (!s.deleted.empty() && !s.inserted.empty() && s.deleted.back() == s.inserted.back()) { s.deleted.pop_back(); s.inserted.pop_back(); }
  std::size_t k = 0;
  while (k < s.deleted.size() && k < s.inserted.size() && s.deleted[k] == s.inserted[k]) k++;
  s.deleted.erase(0, k);
  s.inserted.erase(0, k);
  s.position += k;
  if (!s.deleted.empty() && !s.inserted.empty()) return s;
  if (s.deleted.empty() && s.inserted.empty()) return s;  // an identity: nothing to grow
  // Growing.
  const bool deletion = !s.deleted.empty();
  const std::string bud = deletion ? s.deleted : s.inserted;
  const std::size_t m = bud.size();
  const std::size_t from = s.position;
  const auto blossom = [&](std::size_t i) {
    const std::ptrdiff_t d = static_cast<std::ptrdiff_t>(i) - static_cast<std::ptrdiff_t>(from);
    const std::ptrdiff_t r = ((d % static_cast<std::ptrdiff_t>(m)) + static_cast<std::ptrdiff_t>(m)) % static_cast<std::ptrdiff_t>(m);
    return bud[static_cast<std::size_t>(r)];
  };
  std::size_t a = from, b = from;
  while (b < reference.size() && sequence::dna_letter(reference[b]) == blossom(b)) b++;
  while (a > 0 && sequence::dna_letter(reference[a - 1]) == blossom(a - 1)) a--;
  Spdi out;
  out.sequence = s.sequence;
  out.position = a;
  out.deleted = detail::slice(reference, a, b);
  // The shrunken variant applied to the new deleted sequence.
  const std::size_t at = from - a;
  out.inserted = deletion ? out.deleted.substr(0, at) + out.deleted.substr(at + m) : out.deleted.substr(0, at) + bud + out.deleted.substr(at);
  return out;
}

inline std::string spdi_text(const Spdi& s) { return s.text(); }

struct VcfFields {
  std::size_t pos = 0;  // 1-based
  std::string ref, alt;
  friend bool operator==(const VcfFields&, const VcfFields&) = default;
};

// VCF fields of an SPDI as written (both alleles non-empty, or padded with the base before; at the
// sequence's start, the base after, VCF 4.3 §1.6.1).
inline VcfFields vcf_fields(const Spdi& s, const Codes& reference) {
  VcfFields v{s.position + 1, s.deleted, s.inserted};
  if (v.ref.empty() || v.alt.empty()) {
    if (s.position > 0) {
      const char b = detail::letter(reference, s.position - 1);
      v.ref.insert(v.ref.begin(), b);
      v.alt.insert(v.alt.begin(), b);
      v.pos--;
    } else if (s.deleted.size() < reference.size()) {
      const char b = detail::letter(reference, s.deleted.size());
      v.ref += b;
      v.alt += b;
    }
  }
  return v;
}

}  // namespace genomics::nomenclature
