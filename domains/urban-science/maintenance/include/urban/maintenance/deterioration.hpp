// The deterioration of a city's surfaces through classes of condition, fitted to one rating a surface at its age.
//
// Classes run from 5 (good) to 1 (failed), the City of Melbourne's scale for its road surfaces. A surface is laid in
// class 5 and wears through the classes in order, as a sequential Markov chain in continuous time: out of class k at
// the rate μ_k a year (k = 5, 4, 3, 2), class 1 absorbing until renewed (Madanat and Ben-Akiva 1994, "Optimal
// inspection and repair policies for infrastructure facilities", Transp. Sci. 28, 55–62, for Markov deterioration of
// latent condition). The chain's state index is 5 − class.
//
// Each surface is rated once, at its age when rated: current-status data, whose likelihood is that of tyche's panel
// of a sequential chain observed from class 5 after the age (tyche processes, Kalbfleisch and Lawless 1985). The rates
// are fitted by type of surface (kerbs, channels, footways, carriageways), their errors by Laplace, and the mean years
// in each class are 1/μ_k. The fit's adequacy is read against the observed classes in bands of age (Pearson's χ²):
// a Markov chain with constant rates is a model, and the data may reject it.
//
// A rating reads the condition with error: one class better or worse with probability `confusion` each side (the
// planning model's inspection, declared there), so that the likelihood is tyche's panel with misclassification
// (Jackson et al. 2003). Without it, a surface rated in class 4 on the day it was laid would be impossible, and the
// Hoddle Grid's kerbs hold such surfaces.
//
// Caveats. A renewal that resets the installation date makes a surface young again: the chain then describes surfaces
// since their last renewal, as a renewal should. The ratings are a real number rounded to their class; the noise of a
// rating is not identified from one rating a surface, and is declared.

#pragma once

#include <tyche/processes/panel_markov.hpp>

#include <calx/linalg/dense.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace urban::maintenance {

constexpr std::size_t classes = 5;  // 5 … 1

inline std::size_t state_of(int cls) { return static_cast<std::size_t>(5 - cls); }
inline int class_of(std::size_t state) { return 5 - static_cast<int>(state); }

struct Rated {
  std::string type;
  double age = std::nan("");  // years since installed, when known
  int cls = 0;                // the class of the rating, 1 … 5
  double rating = 0, area = 0;
  double since = 0;           // years from the rating to the time of planning
};

// The rates' family: one rate for every stage (1 parameter), rates in geometric progression μ_k = μ ρ^k (2), or four
// free rates (4). One rating a surface identifies the later stages poorly where few surfaces have reached them; the
// family is chosen by Akaike's criterion among the fits that converge with finite errors and rates between 10⁻⁴ and 10
// a year.
enum class Family { common = 1, geometric = 2, free = 4 };

inline std::vector<double> log_rates(Family f, const std::vector<double>& theta) {
  std::vector<double> r(4);
  for (std::size_t k = 0; k < 4; k++) r[k] = f == Family::common ? theta[0] : f == Family::geometric ? theta[0] + k * theta[1] : theta[k];
  return r;
}

struct Deterioration {
  std::string type;
  Family family = Family::free;
  std::array<double, 4> rate{};          // μ out of classes 5, 4, 3, 2, a year
  std::vector<double> theta, error;      // the family's parameters (log μ, and log ρ) and Laplace's errors
  double log_likelihood = 0, aic = 0;
  std::size_t surfaces = 0;
  bool converged = false;
  std::vector<std::pair<Family, double>> tried;  // each family's AIC (NaN if refused)

  calx::linalg::Matrix generator() const {
    std::vector<double> theta(4);
    for (std::size_t k = 0; k < 4; k++) theta[k] = std::log(rate[k]);
    return tyche::processes::sequential_generator(theta);
  }
  // The law of the class (states 5 − class) after t years from class `from`.
  std::vector<double> law(double t, int from = 5) const { return tyche::processes::transition_row(generator(), state_of(from), t); }
  // Mean years in class k (5, 4, 3, 2).
  double residence(int k) const { return 1 / rate[state_of(k)]; }
};

// The misclassification of a rating: one class better or worse with probability `confusion` each side.
inline tyche::processes::Misclassification misreading(double confusion) {
  tyche::processes::Misclassification E(5, std::vector<double>(5, 0.0));
  for (std::size_t k = 0; k < 5; k++) {
    const double up = k > 0 ? confusion : 0, down = k + 1 < 5 ? confusion : 0;
    E[k][k] = 1 - up - down;
    if (k > 0) E[k][k - 1] = up;
    if (k + 1 < 5) E[k][k + 1] = down;
  }
  return E;
}

// The fit within one family; converged only with finite errors and rates in [10⁻⁴, 10] a year (its AIC NaN else).
// The families are nested, and each starts from the fit of the one inside it (one rate, then (that rate, ρ = 1), then
// the geometric fit's four rates), so that BFGS begins near its maximum: from a distant start its first step, scaled
// by a gradient of thousands of observations, can leave the region where the likelihood is finite.
inline Deterioration fit_family(const std::vector<Rated>& rows, const std::string& type, Family f, double confusion = 0.1, std::vector<double> start = {}) {
  std::vector<tyche::processes::PanelObservation> data;
  for (const auto& r : rows)
    if (r.type == type && std::isfinite(r.age) && r.cls >= 1 && r.cls <= 5) data.push_back({0, state_of(r.cls), r.age, 1});
  Deterioration D;
  D.type = type, D.family = f, D.surfaces = data.size(), D.aic = std::nan("");
  if (data.empty()) return D;
  const std::size_t p = static_cast<std::size_t>(f);
  const tyche::processes::Generator Q = [f](const std::vector<double>& th) { return tyche::processes::sequential_generator(log_rates(f, th)); };
  if (start.empty()) {
    if (f == Family::common) start = {std::log(0.05)};
    else {
      const auto inner = fit_family(rows, type, f == Family::geometric ? Family::common : Family::geometric, confusion);
      if (f == Family::geometric) start = {inner.converged ? inner.theta[0] : std::log(0.05), 0.0};
      else {
        start.assign(4, std::log(0.05));
        if (inner.converged)
          for (std::size_t k = 0; k < 4; k++) start[k] = std::log(inner.rate[k]);
      }
    }
  }
  try {
    const auto F = tyche::processes::fit_panel(Q, data, start, 500, misreading(confusion));
    const auto lr = log_rates(f, F.theta);
    bool ok = F.converged;
    for (const double e : F.error) ok = ok && std::isfinite(e);
    for (const double x : lr) ok = ok && x >= std::log(1e-4) && x <= std::log(10.0);
    for (std::size_t k = 0; k < 4; k++) D.rate[k] = std::exp(lr[k]);
    D.theta = F.theta, D.error = F.error, D.log_likelihood = F.log_likelihood, D.converged = ok;
    if (ok) D.aic = 2 * static_cast<double>(p) - 2 * F.log_likelihood;
  } catch (const std::invalid_argument&) {
  }
  return D;
}

// The family chosen by AIC among those that converge.
inline Deterioration fit(const std::vector<Rated>& rows, const std::string& type, double confusion = 0.1) {
  Deterioration best;
  best.type = type;
  std::vector<std::pair<Family, double>> tried;
  bool any = false;
  for (const Family f : {Family::common, Family::geometric, Family::free}) {
    auto D = fit_family(rows, type, f, confusion);
    tried.push_back({f, D.aic});
    best.surfaces = D.surfaces;
    if (std::isfinite(D.aic) && (!any || D.aic < best.aic)) best = std::move(D), any = true;
  }
  best.tried = tried;
  if (!any) best.converged = false;
  return best;
}

// The fit against the classes observed in bands of age: Pearson's χ² over the cells with at least five expected, and
// its degrees of freedom (cells less one a band, less the four rates).
struct Adequacy {
  double chi2 = 0;
  int df = 0, cells = 0;
  std::vector<std::array<double, 5>> observed, expected;  // per band, by state (class 5 … 1)
  std::vector<double> band_start;
};

// The ratings expected under the fitted chain and the misreading, against those observed.
inline Adequacy adequacy(const Deterioration& D, const std::vector<Rated>& rows, double band = 5, double confusion = 0.1) {
  Adequacy A;
  double top = 0;
  for (const auto& r : rows)
    if (r.type == D.type && std::isfinite(r.age)) top = std::max(top, r.age);
  const std::size_t bands = static_cast<std::size_t>(std::floor(top / band)) + 1;
  A.observed.assign(bands, {}), A.expected.assign(bands, {});
  for (std::size_t b = 0; b < bands; b++) A.band_start.push_back(b * band);
  for (const auto& r : rows) {
    if (r.type != D.type || !std::isfinite(r.age)) continue;
    const std::size_t b = static_cast<std::size_t>(std::floor(r.age / band));
    A.observed[b][state_of(r.cls)] += 1;
    const auto p = D.law(r.age);
    const auto E = misreading(confusion);
    for (std::size_t s = 0; s < 5; s++)
      for (std::size_t j = 0; j < 5; j++) A.expected[b][j] += p[s] * E[s][j];
  }
  int used_bands = 0;
  for (std::size_t b = 0; b < bands; b++) {
    int cells = 0;
    for (std::size_t s = 0; s < 5; s++)
      if (A.expected[b][s] >= 5) A.chi2 += (A.observed[b][s] - A.expected[b][s]) * (A.observed[b][s] - A.expected[b][s]) / A.expected[b][s], cells++;
    if (cells > 0) A.cells += cells, A.df += cells - 1, used_bands++;
  }
  A.df = std::max(1, A.df - 4);
  return A;
}

}  // namespace urban::maintenance
