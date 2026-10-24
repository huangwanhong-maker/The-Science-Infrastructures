// Checks of the substitution models and rates across sites: for JC69, K80, F81, HKY85, TN93 and
// GTR, the rate at stationarity is 1, P(0) = I, every row of P(t) sums to 1, P(s + t) =
// P(s)P(t), detailed balance π_i P_ij = π_j P_ji, P(t) → Π, P(t) against calx::linalg::expm(Qt),
// and dP/dt = QP, d²P/dt² = Q²P; JC69's and K80's closed forms; Yang's (1994) discrete Γ rates
// for α = 0.5 and K = 4, their exact values at α = 1, the median variant, and the invariant class.

#include <evolution/phylogeny/models.hpp>

#include <calx/linalg/expm.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

using namespace evolution::phylogeny;

namespace {

using M4 = std::array<double, 16>;

M4 product(const M4& a, const M4& b) {
  M4 c{};
  for (int i = 0; i < 4; i++)
    for (int k = 0; k < 4; k++)
      for (int j = 0; j < 4; j++) c[i * 4 + j] += a[i * 4 + k] * b[k * 4 + j];
  return c;
}

double max_diff(const M4& a, const M4& b) {
  double m = 0;
  for (int i = 0; i < 16; i++) m = std::max(m, std::fabs(a[i] - b[i]));
  return m;
}

M4 P(const Model& m, double t) {
  M4 x{};
  m.transition(t, x.data());
  return x;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  const std::array<double, 4> pi{0.35, 0.15, 0.2, 0.3};
  const std::vector<std::pair<std::string, Model>> models{
      {"JC69", jc69()},
      {"K80 (κ = 3)", k80(3)},
      {"F81", f81(pi)},
      {"HKY85 (κ = 5)", hky85(5, pi)},
      {"TN93 (κ₁ = 8, κ₂ = 2)", tn93(8, 2, pi)},
      {"GTR", gtr({1.3, 4.1, 0.7, 1.9, 5.2, 1.0}, {0.1, 0.4, 0.3, 0.2})},
  };
  const std::vector<double> times{0.001, 0.05, 0.3, 1.0, 4.0};
  for (const auto& [name, m] : models) {
    double rate = 0;
    for (int i = 0; i < 4; i++) rate -= m.pi[static_cast<std::size_t>(i)] * m.Q[static_cast<std::size_t>(i * 4 + i)];
    s.check(name + ": the rate at stationarity is 1", rate - 1, 1e-14);
    M4 I{};
    for (int i = 0; i < 4; i++) I[static_cast<std::size_t>(i * 5)] = 1;
    s.check(name + ": P(0) = I", max_diff(P(m, 0), I), 1e-14);
    double rows = 0, balance = 0, semigroup = 0, vs_expm = 0, d1 = 0, d2 = 0;
    for (double t : times) {
      const M4 p = P(m, t);
      for (int i = 0; i < 4; i++) {
        double r = 0;
        for (int j = 0; j < 4; j++) {
          r += p[static_cast<std::size_t>(i * 4 + j)];
          balance = std::max(balance, std::fabs(m.pi[static_cast<std::size_t>(i)] * p[static_cast<std::size_t>(i * 4 + j)] -
                                                m.pi[static_cast<std::size_t>(j)] * p[static_cast<std::size_t>(j * 4 + i)]));
        }
        rows = std::max(rows, std::fabs(r - 1));
      }
      semigroup = std::max(semigroup, max_diff(P(m, t + 0.37), product(p, P(m, 0.37))));
      calx::linalg::Matrix Qt(4);
      for (int i = 0; i < 16; i++) Qt.a[static_cast<std::size_t>(i)] = m.Q[static_cast<std::size_t>(i)] * t;
      const auto E = calx::linalg::expm(Qt);
      M4 e{};
      for (int i = 0; i < 16; i++) e[static_cast<std::size_t>(i)] = E.a[static_cast<std::size_t>(i)];
      vs_expm = std::max(vs_expm, max_diff(p, e));
      M4 q{}, dp{}, d2p{}, pp{};
      m.transition(t, pp.data(), dp.data(), d2p.data());
      for (int i = 0; i < 16; i++) q[static_cast<std::size_t>(i)] = m.Q[static_cast<std::size_t>(i)];
      d1 = std::max(d1, max_diff(dp, product(q, p)));
      d2 = std::max(d2, max_diff(d2p, product(product(q, q), p)));
    }
    s.check(name + ": every row of P(t) sums to 1", rows, 1e-14);
    s.check(name + ": detailed balance, π_i P_ij(t) = π_j P_ji(t)", balance, 1e-15);
    s.check(name + ": P(s + t) = P(s) P(t)", semigroup, 1e-14);
    s.check(name + ": P(t) = expm(Qt) by Padé (calx::linalg::expm)", vs_expm, 1e-13);
    s.check(name + ": dP/dt = Q P(t)", d1, 1e-13);
    s.check(name + ": d²P/dt² = Q² P(t)", d2, 1e-13);
    double limit = 0;
    const M4 far = P(m, 400);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) limit = std::max(limit, std::fabs(far[static_cast<std::size_t>(i * 4 + j)] - m.pi[static_cast<std::size_t>(j)]));
    s.check(name + ": P(t) → Π, every row π", limit, 1e-14);
  }

  // Closed forms. JC69: P_ii(t) = 1/4 + 3/4 e^{−4t/3}. K80 with κ, normalised: α = κ/(κ + 2) to
  // the transition, β = 1/(κ + 2) to each transversion, P_ts = 1/4 + 1/4 e^{−4βt} − 1/2 e^{−2(α+β)t},
  // P_tv = 1/4 − 1/4 e^{−4βt} (Kimura 1980; Yang 2014, ch. 1).
  {
    double jc = 0, kts = 0, ktv = 0;
    const Model k = k80(3);
    for (double t : times) {
      const M4 p = P(jc69(), t);
      jc = std::max({jc, std::fabs(p[0] - (0.25 + 0.75 * std::exp(-4 * t / 3))), std::fabs(p[1] - (0.25 - 0.25 * std::exp(-4 * t / 3)))});
      const double a = 3.0 / 5, b = 1.0 / 5;
      const M4 q = P(k, t);
      kts = std::max(kts, std::fabs(q[0 * 4 + 2] - (0.25 + 0.25 * std::exp(-4 * b * t) - 0.5 * std::exp(-2 * (a + b) * t))));
      ktv = std::max(ktv, std::fabs(q[0 * 4 + 1] - (0.25 - 0.25 * std::exp(-4 * b * t))));
    }
    s.check("JC69's P(t) is its closed form", jc, 1e-15);
    s.check("K80's transition probability is its closed form", kts, 1e-15);
    s.check("K80's transversion probability is its closed form", ktv, 1e-15);
  }

  // Yang (1994): with α = 0.5 and four categories, the mean rates 0.0334, 0.2519, 0.8203, 2.8944
  // (published to four places; to six, 0.033388, 0.251916, 0.820268, 2.894428, from the same
  // formulas evaluated by SciPy).
  {
    const auto r = discrete_gamma(0.5, 4);
    const double pub[4] = {0.0334, 0.2519, 0.8203, 2.8944}, six[4] = {0.033388, 0.251916, 0.820268, 2.894428};
    double e4 = 0, e6 = 0, mean = 0;
    for (int k = 0; k < 4; k++) {
      e4 = std::max(e4, std::fabs(r[static_cast<std::size_t>(k)] - pub[k]));
      e6 = std::max(e6, std::fabs(r[static_cast<std::size_t>(k)] - six[k]));
      mean += r[static_cast<std::size_t>(k)] / 4;
    }
    s.check("Yang's discrete Γ rates for α = 0.5, K = 4, as published (four places)", e4, 5e-5);
    s.check("and to six places", e6, 5e-7);
    s.check("their mean is 1", mean - 1, 1e-13);
    // α = 1: the exponential distribution, cut at ln(4/3), ln 2, ln 4; the means are
    // 1 − 3 ln(4/3), 1 + 3 ln(4/3) − 2 ln 2, 1 and 1 + ln 4.
    const auto e = discrete_gamma(1, 4);
    const double l43 = std::log(4.0 / 3), l2 = std::log(2.0);
    const double exact[4] = {1 - 3 * l43, 1 + 3 * l43 - 2 * l2, 1, 1 + 2 * l2};
    double ex = 0;
    for (int k = 0; k < 4; k++) ex = std::max(ex, std::fabs(e[static_cast<std::size_t>(k)] - exact[k]));
    s.check("at α = 1 they are exactly 1 − 3 ln(4/3), 1 + 3 ln(4/3) − 2 ln 2, 1, 1 + ln 4", ex, 1e-12);
    // The median variant, rescaled to mean 1 (as PAML computes it): 0.029078, 0.280715, 0.924773, 2.765435.
    const auto m = discrete_gamma(0.5, 4, true);
    const double med[4] = {0.029078, 0.280715, 0.924773, 2.765435};
    double em = 0;
    for (int k = 0; k < 4; k++) em = std::max(em, std::fabs(m[static_cast<std::size_t>(k)] - med[k]));
    s.check("the median rates for α = 0.5, K = 4, rescaled to mean 1", em, 5e-7);
    const auto big = discrete_gamma(200, 8);
    double spread = 0;
    for (double x : big) spread = std::max(spread, std::fabs(x - 1));
    s.require("as α grows the rates close in on 1", spread < 0.25 && spread > 0);
  }
  {
    const Rates r = gamma_rates(0.5, 4, 0.2);
    double w = 0, mean = 0;
    for (std::size_t k = 0; k < r.size(); k++) {
      w += r.weight[k];
      mean += r.weight[k] * r.rate[k];
    }
    s.require("Γ + I: five classes, the last of rate 0 and weight 0.2", r.size() == 5 && r.rate[4] == 0 && r.weight[4] == 0.2);
    s.check("the weights sum to 1", w - 1, 1e-15);
    s.check("and the mean rate is 1", mean - 1, 1e-13);
  }
  {
    bool threw = false;
    try {
      (void)f81({0.5, 0.5, 0, 0});
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    s.require("a frequency of zero is an invalid argument", threw);
  }
  return s.done();
}
