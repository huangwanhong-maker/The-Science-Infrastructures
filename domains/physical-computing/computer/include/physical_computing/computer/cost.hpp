// The cost of a verified result (theory note §1.5; ADR 0076 §2.3).
//
// cost = (E_fabrication + E_materials) / N_results over the lifetime + E_operation + E_conversion + E_verification,
//
// in joules, each term with how it is known: measured, derived from a model, or declared as an assumption,
// and its source. The physical floors come from physis: Landauer's bound per bit erased and the thermodynamic
// uncertainty relation's bound per precision of a measured current. Conversion to the digital world is
// estimated by an analogue-to-digital converter's figure of merit (Walden's: the energy per conversion step,
// E = FoM · 2^ENOB per sample), whose value is declared with its source (Murmann's survey of converters). A
// term not known is not set to zero: it is absent, and the report says so.
#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <physis/statistical_mechanics/fluctuations.hpp>

namespace physical_computing::computer {

enum class Basis { measured, derived, declared };

enum class Stage { fabrication, materials, operation, conversion, verification };

struct CostTerm {
  std::string name;
  Stage stage = Stage::operation;
  double joules = 0;  // per result for operation, conversion and verification; in total for fabrication and materials
  Basis basis = Basis::declared;
  std::string source;
};

class CostAccount {
 public:
  void add(CostTerm term) {
    if (!(term.joules >= 0) || !std::isfinite(term.joules)) throw std::invalid_argument("physical_computing::computer::CostAccount: a finite, non-negative energy");
    terms_.push_back(std::move(term));
  }
  const std::vector<CostTerm>& terms() const { return terms_; }

  // The energy of one verified result, the fixed costs spread over a lifetime of results.
  double per_verified_result(double lifetime_results) const {
    if (!(lifetime_results > 0)) throw std::invalid_argument("physical_computing::computer::CostAccount: a positive number of results over the lifetime");
    double fixed = 0, each = 0;
    for (const auto& t : terms_) (t.stage == Stage::fabrication || t.stage == Stage::materials ? fixed : each) += t.joules;
    return fixed / lifetime_results + each;
  }
  // Whether every stage has at least one term: an account with a stage missing is incomplete, not cheap.
  std::vector<Stage> missing_stages() const {
    std::vector<Stage> out;
    for (Stage s : {Stage::fabrication, Stage::materials, Stage::operation, Stage::conversion, Stage::verification}) {
      bool found = false;
      for (const auto& t : terms_) found = found || t.stage == s;
      if (!found) out.push_back(s);
    }
    return out;
  }
  // The share of declared (not measured or derived) energy in a result's cost.
  double declared_share(double lifetime_results) const {
    CostAccount d;
    for (const auto& t : terms_)
      if (t.basis == Basis::declared) d.terms_.push_back(t);
    const double total = per_verified_result(lifetime_results);
    return total > 0 ? d.per_verified_result(lifetime_results) / total : 0;
  }

 private:
  std::vector<CostTerm> terms_;
};

// The floors, from physis.
inline double landauer_floor(double temperature, double bits_erased) {
  return bits_erased * physis::statistical_mechanics::landauer_bound(temperature);
}
inline double precision_floor(double temperature, double relative_precision) {
  return physis::statistical_mechanics::precision_cost_bound(temperature, relative_precision);
}

// The energy of one conversion by an ADC of `effective_bits`, at Walden's figure of merit (J per step).
inline double conversion_energy(double effective_bits, double figure_of_merit) { return figure_of_merit * std::pow(2.0, effective_bits); }

}  // namespace physical_computing::computer
