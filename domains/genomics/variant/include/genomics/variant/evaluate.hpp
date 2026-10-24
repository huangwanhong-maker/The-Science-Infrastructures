// The evaluation of calls against a truth set: true and false positives and false negatives,
// precision, recall and F₁, by type of variant (SNV, insertion, deletion, other), and the
// concordance of genotypes among the true positives.
//
// - Matching. Both sides are normalised (normalize.hpp: left-aligned and parsimonious) against
//   the reference, so that two representations of one indel in a repeat meet; a call matches a
//   truth variant when the two have the same contig, position, REF and ALT. A multi-allelic call
//   is split into one call per ALT, each with the number of copies of that ALT its genotype has.
//   *Remark.* This is allele matching after normalisation, the first level of the comparisons of
//   the GA4GH benchmarking team (Krusche et al. 2019, Nature Biotechnology 37, 555–560); it
//   does not compare haplotypes, so that a complex variant written as an SNV and an indel next
//   to each other in one set and as one block in the other does not match.
// - Precision TP / (TP + FP), recall TP / (TP + FN), F₁ their harmonic mean; genotype
//   concordance, of the true positives whose dosage (1 heterozygous, 2 homozygous) is the
//   truth's. By default only calls whose FILTER is PASS count.
//
// Complexity: O((c + t) log t) for c calls and t truth variants, plus their normalisation. The
// truth is any type with position, ref, alt and on[2] (simulate::Variant), so that the caller
// does not depend on the simulator.

#pragma once

#include <genomics/variant/caller.hpp>
#include <genomics/variant/normalize.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

namespace genomics::variant {

enum class VariantType { snv, insertion, deletion, other };

inline const char* type_name(VariantType t) {
  switch (t) {
    case VariantType::snv: return "SNV";
    case VariantType::insertion: return "insertion";
    case VariantType::deletion: return "deletion";
    case VariantType::other: return "other";
  }
  return "?";
}

// The type of a normalised biallelic variant.
inline VariantType type_of(std::string_view ref, std::string_view alt) {
  if (ref.size() == 1 && alt.size() == 1) return VariantType::snv;
  if (ref.size() == 1 && alt.size() > 1 && alt[0] == ref[0]) return VariantType::insertion;
  if (alt.size() == 1 && ref.size() > 1 && ref[0] == alt[0]) return VariantType::deletion;
  return VariantType::other;
}

struct TruthVariant {
  std::uint32_t contig = 0;
  std::size_t position = 0;
  std::string ref, alt;
  int dosage = 1;  // copies of ALT: 1 heterozygous, 2 homozygous
};

template <class V>
TruthVariant truth_variant(const V& v, std::uint32_t contig = 0) {
  return {contig, static_cast<std::size_t>(v.position), std::string(v.ref), std::string(v.alt), static_cast<int>(v.on[0]) + static_cast<int>(v.on[1])};
}
template <class V>
std::vector<TruthVariant> truth_set(const std::vector<V>& variants, std::uint32_t contig = 0) {
  std::vector<TruthVariant> t;
  t.reserve(variants.size());
  for (const auto& v : variants) t.push_back(truth_variant(v, contig));
  return t;
}

struct Counts {
  std::size_t tp = 0, fp = 0, fn = 0, genotype_match = 0;
  double precision() const { return tp + fp ? static_cast<double>(tp) / static_cast<double>(tp + fp) : 1; }
  double recall() const { return tp + fn ? static_cast<double>(tp) / static_cast<double>(tp + fn) : 1; }
  double f1() const {
    const double p = precision(), r = recall();
    return p + r > 0 ? 2 * p * r / (p + r) : 0;
  }
  double genotype_concordance() const { return tp ? static_cast<double>(genotype_match) / static_cast<double>(tp) : 1; }
};

struct Comparison {
  Counts snv, insertion, deletion, other, all;
  std::vector<TruthVariant> false_positives, false_negatives;  // normalised
  Counts& of(VariantType t) {
    switch (t) {
      case VariantType::snv: return snv;
      case VariantType::insertion: return insertion;
      case VariantType::deletion: return deletion;
      default: return other;
    }
  }
  const Counts& of(VariantType t) const { return const_cast<Comparison*>(this)->of(t); }
};

inline Comparison compare(const std::vector<Call>& calls, const std::vector<TruthVariant>& truth, const std::vector<Codes>& contigs, bool pass_only = true) {
  const auto norm = [&](TruthVariant v) {
    if (v.contig >= contigs.size()) throw std::invalid_argument("genomics::variant::compare: a variant on a contig not given");
    Normalized n = normalize(v.position, v.ref, v.alt, contigs[v.contig]);
    v.position = n.position;
    v.ref = std::move(n.ref);
    v.alt = std::move(n.alts[0]);
    return v;
  };
  const auto key = [](const TruthVariant& v) { return std::tie(v.contig, v.position, v.ref, v.alt); };
  std::vector<TruthVariant> t;
  t.reserve(truth.size());
  for (const auto& v : truth) t.push_back(norm(v));
  std::sort(t.begin(), t.end(), [&](const TruthVariant& a, const TruthVariant& b) { return key(a) < key(b); });
  std::vector<char> found(t.size(), 0);
  Comparison out;
  for (const auto& c : calls) {
    if (pass_only && !c.passed()) continue;
    for (std::size_t a = 0; a < c.alts.size(); a++) {
      const int dosage = (c.genotype[0] == static_cast<int>(a + 1)) + (c.genotype[1] == static_cast<int>(a + 1));
      if (dosage == 0) continue;
      const TruthVariant v = norm({c.contig, c.position, c.ref, c.alts[a], dosage});
      const auto it = std::lower_bound(t.begin(), t.end(), v, [&](const TruthVariant& x, const TruthVariant& y) { return key(x) < key(y); });
      if (it != t.end() && key(*it) == key(v) && !found[static_cast<std::size_t>(it - t.begin())]) {
        found[static_cast<std::size_t>(it - t.begin())] = 1;
        Counts& k = out.of(type_of(it->ref, it->alt));
        k.tp++;
        out.all.tp++;
        if (it->dosage == v.dosage) { k.genotype_match++; out.all.genotype_match++; }
      } else {
        out.of(type_of(v.ref, v.alt)).fp++;
        out.all.fp++;
        out.false_positives.push_back(v);
      }
    }
  }
  for (std::size_t i = 0; i < t.size(); i++)
    if (!found[i]) {
      out.of(type_of(t[i].ref, t[i].alt)).fn++;
      out.all.fn++;
      out.false_negatives.push_back(t[i]);
    }
  return out;
}

}  // namespace genomics::variant
