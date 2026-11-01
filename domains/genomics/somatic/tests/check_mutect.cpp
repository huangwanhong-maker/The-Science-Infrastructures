// Checks of MuTect's statistics by hand: the tumour and normal LODs of small pileups computed from
// the model written out term by term; the thresholds from their priors (6.3, 2.2, 5.5, Cibulskis et
// al. 2013); the sensitivity against a sum over every count of reads; the depths of the callable
// territory (14 and 8) as the least that reach their powers; the indel form on counts; the median
// and median absolute deviation; the read offset of a reference base through a CIGAR.

#include <genomics/somatic/mutect.hpp>

#include <calx/special_functions/gamma.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace genomics;
using sequence::A;
using sequence::C;
using sequence::G;
using sequence::T;

int main() {
  calx::verification::Suite s;

  // ---- LOD_T and LOD_N of a pileup ---------------------------------------------------------------
  // Reference A, candidate C: seven A at Q30, three C at Q20, one G at Q25.
  std::vector<variant::BaseCall> b;
  for (int i = 0; i < 7; i++) b.push_back({A, 30});
  for (int i = 0; i < 3; i++) b.push_back({C, 20});
  b.push_back({G, 25});
  {
    const double e30 = 1e-3, e20 = 1e-2, e25 = std::pow(10.0, -2.5), f = 3.0 / 11;
    const double alt1 = std::log10(f * (1 - e20) + (1 - f) * e20 / 3), alt0 = std::log10(e20 / 3);
    const double ref1 = std::log10(f * e30 / 3 + (1 - f) * (1 - e30)), ref0 = std::log10(1 - e30);
    // The other base is e/3 under both models and cancels.
    const double hand = 3 * (alt1 - alt0) + 7 * (ref1 - ref0);
    s.check("LOD_T of 7 A, 3 C and a G against the model term by term", somatic::tumour_lod(b, A, C) - hand, 1e-12);
    const double other = std::log10(e25 / 3);
    s.check("log10 L(M_f) sums the other base's e/3", somatic::log10_likelihood(b, A, C, f) - (3 * alt1 + 7 * ref1 + other), 1e-12);
    const double half_alt = std::log10(0.5 * (1 - e20) + 0.5 * e20 / 3), half_ref = std::log10(0.5 * e30 / 3 + 0.5 * (1 - e30));
    const double lodn = 3 * (alt0 - half_alt) + 7 * (ref0 - half_ref);
    s.check("LOD_N against the model at f = 1/2 term by term", somatic::normal_lod(b, A, C) - lodn, 1e-12);
    s.require("a normal with three reads of C of eleven is not somatic (LOD_N < 2.2)", somatic::normal_lod(b, A, C) < 2.2);
  }

  // ---- the thresholds from their priors ------------------------------------------------------------
  {
    const double theta_t = std::log10(2.0) - std::log10(3e-6 / 3);
    s.check("θ_T = log10 δ_T − log10 P(m, f) = 6.30 (δ_T = 2, P(m) = 3e-6/3)", theta_t - 6.3, 0.005);
    const double theta_n = std::log10(10.0) + std::log10(5e-5 / 3e-6), theta_known = std::log10(10.0) + std::log10(0.095 / 3e-6);
    s.check("θ_N = log10 δ_N + log10 P(germline)/P(somatic) = 2.2 off dbSNP", theta_n - 2.2, 0.03);
    s.check("θ_N = 5.5 at a dbSNP site (P(germline) = 0.095)", theta_known - 5.5, 0.005);
    const somatic::MutectOptions o;
    s.require("the options carry 6.3, 2.2 and 5.5", o.theta_tumour == 6.3 && o.theta_normal == 2.2 && o.theta_normal_known == 5.5);
    // One read of m can never pass: its LOD is below log10(3/e) at Q41.
    std::vector<variant::BaseCall> one{{C, 41}};
    s.require("one read of m at Q41 stays below θ_T", somatic::tumour_lod(one, A, C) < 6.3, std::to_string(somatic::tumour_lod(one, A, C)));
  }

  // ---- sensitivity -------------------------------------------------------------------------------------
  {
    double worst = 0;
    for (std::size_t n : {5u, 14u, 30u, 80u})
      for (double f : {0.05, 0.1, 0.3, 0.5})
        for (double e : {1e-3, 1e-2}) {
          double brute = 0;
          const double pm = f * (1 - e) + (1 - f) * e / 3;
          for (std::size_t k = 0; k <= n; k++)
            if (somatic::lod_at(k, n, e) >= 6.3)
              brute += std::exp(calx::special_functions::log_choose(static_cast<double>(n), static_cast<double>(k)) + static_cast<double>(k) * std::log(pm) +
                                static_cast<double>(n - k) * std::log1p(-pm));
          worst = std::max(worst, std::fabs(somatic::power(n, f, e, 6.3) - brute));
        }
    s.check("power: P(LOD ≥ 6.3) against the sum over every count", worst, 1e-10);
    s.require("14 tumour reads detect f = 0.3 at Q30 with power ≥ 0.8, 13 do not",
              somatic::power(14, 0.3, 1e-3, 6.3) >= 0.8 && somatic::power(13, 0.3, 1e-3, 6.3) < 0.8,
              std::to_string(somatic::power(13, 0.3, 1e-3, 6.3)) + " " + std::to_string(somatic::power(14, 0.3, 1e-3, 6.3)));
    std::vector<variant::BaseCall> n7(7, {A, 30}), n8(8, {A, 30});
    s.require("8 reference reads at Q30 reach LOD_N ≥ 2.2, 7 do not", somatic::normal_lod(n8, A, C) >= 2.2 && somatic::normal_lod(n7, A, C) < 2.2,
              std::to_string(somatic::normal_lod(n7, A, C)) + " " + std::to_string(somatic::normal_lod(n8, A, C)));
  }

  // ---- the counting form -------------------------------------------------------------------------------
  {
    const double eps = 0.005, f = 4.0 / 20;
    const double hand = 4 * (std::log10(f * (1 - eps) + (1 - f) * eps) - std::log10(eps)) + 16 * (std::log10(1 - (f * (1 - eps) + (1 - f) * eps)) - std::log10(1 - eps));
    s.check("indel LOD_T of 4 carriers in 20 against the counting model", somatic::tumour_lod_counts(4, 20, eps) - hand, 1e-12);
    const double hn = 20 * (std::log10(1 - eps) - std::log10(0.5));
    s.check("indel LOD_N of 20 reads without the indel", somatic::normal_lod_counts(0, 20, eps) - hn, 1e-12);
  }

  // ---- median and MAD; offsets in reads ---------------------------------------------------------------
  {
    const auto [m, d] = somatic::median_mad({1, 2, 3, 4, 100});
    s.require("median and MAD of {1, 2, 3, 4, 100}: 3 and 1", m == 3 && d == 1);
    const auto [m2, d2] = somatic::median_mad({4, 1, 3, 2});
    s.require("median and MAD of {1, 2, 3, 4}: 2.5 and 1", m2 == 2.5 && d2 == 1);
    mapping::Mapping x;
    x.alignment.query_begin = 2;
    x.alignment.target_begin = 100;
    x.alignment.cigar = {{'M', 10}, {'I', 2}, {'M', 5}, {'D', 3}, {'M', 10}};
    x.target_begin = 100;
    x.target_end = 128;
    bool ok = somatic::query_offset(x, 100) == 2 && somatic::query_offset(x, 109) == 11 && somatic::query_offset(x, 110) == 14 && somatic::query_offset(x, 114) == 18 &&
              somatic::query_offset(x, 115) == -1 && somatic::query_offset(x, 117) == -1 && somatic::query_offset(x, 118) == 19 && somatic::query_offset(x, 127) == 28 &&
              somatic::query_offset(x, 128) == -1 && somatic::query_offset(x, 99) == -1;
    s.require("query offsets through 10M2I5M3D10M after two clipped bases", ok);
  }
  return s.done();
}
