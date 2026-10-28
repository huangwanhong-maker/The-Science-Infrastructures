// Schmidt semi-normalised associated Legendre functions of sin φ′, and their derivatives in φ′: the
// functions of the WMM2025 report's eq. 5 (Chulliat et al. 2025, §1.2), with Pₙ,ₘ the associated
// Legendre functions without the Condon–Shortley phase.
//
// They are calx's (calx/special_functions/legendre.hpp, where their recurrences and checks are),
// since they name no object in nature; this header gives them the names the module uses. For a
// geocentric latitude φ′, schmidt_functions(n, φ′) holds P̆ₙᵐ(sin φ′), dP̆ₙᵐ/dφ′ and P̆ₙᵐ/cos φ′
// (finite at the poles, which the east component of the field needs; the report's eq. 11 and §1.4).

#pragma once

#include <calx/special_functions/legendre.hpp>

namespace geophysics::geomagnetism {

using SchmidtFunctions = calx::special_functions::SchmidtLegendre;

inline SchmidtFunctions schmidt_functions(int degree, double latitude_rad) {
  return calx::special_functions::schmidt_legendre(degree, latitude_rad);
}

}  // namespace geophysics::geomagnetism
