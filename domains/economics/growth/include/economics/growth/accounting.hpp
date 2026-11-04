// Growth accounting (Solow 1957) and the perpetual inventory of capital.
//
// With Y = A K^α L^(1−α) (or shares varying in time), the growth of output splits into the
// contributions of capital, of labour and of total factor productivity, the Solow residual:
//
//   Δ ln Y_t = ᾱ_t Δ ln K_t + (1 − ᾱ_t) Δ ln L_t + Δ ln A_t,                     (Solow 1957, eq. 2a)
//
// ᾱ_t the capital share, constant or the average of the two periods' shares (the Törnqvist index,
// exact for a translog technology: Diewert 1976). The decomposition adds up exactly by construction:
// the residual is what the factors leave. Per worker, Δ ln(Y/L) = ᾱ Δ ln(K/L) + Δ ln A.
//
// The capital stock, where only investment is recorded, by the perpetual inventory method
// (OECD 2009, Measuring Capital, ch. 6): K_{t+1} = (1 − δ) K_t + I_t, started at the steady stock
// of the first investment, K_0 = I_0 / (g + δ), g the growth of investment over the first years
// (Harberger 1978).

#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

namespace economics::growth {

struct GrowthAccounts {
  std::vector<double> output, capital, labour, tfp;  // contributions to Δ ln Y, one for each step t−1 → t
  std::vector<double> share;                         // ᾱ_t used at each step
  double mean_output = 0, mean_capital = 0, mean_labour = 0, mean_tfp = 0;  // a period, on average
  double tfp_share = 0;                              // the residual's share of output growth over the span
};

// Y, K and L over the same periods (positive); capital shares one for all, or one for each period.
inline GrowthAccounts growth_accounts(const std::vector<double>& Y, const std::vector<double>& K, const std::vector<double>& L,
                                      const std::vector<double>& capital_share) {
  const std::size_t n = Y.size();
  if (n < 2 || K.size() != n || L.size() != n) throw std::invalid_argument("economics::growth::growth_accounts: Y, K and L over the same two or more periods");
  if (capital_share.size() != 1 && capital_share.size() != n) throw std::invalid_argument("economics::growth::growth_accounts: one capital share, or one for each period");
  for (std::size_t t = 0; t < n; t++)
    if (!(Y[t] > 0 && K[t] > 0 && L[t] > 0)) throw std::invalid_argument("economics::growth::growth_accounts: output, capital and labour are positive");
  auto share = [&](std::size_t t) { return capital_share.size() == 1 ? capital_share[0] : capital_share[t]; };
  GrowthAccounts G;
  for (std::size_t t = 1; t < n; t++) {
    const double a = 0.5 * (share(t - 1) + share(t));
    const double gy = std::log(Y[t] / Y[t - 1]), gk = std::log(K[t] / K[t - 1]), gl = std::log(L[t] / L[t - 1]);
    G.share.push_back(a);
    G.output.push_back(gy);
    G.capital.push_back(a * gk);
    G.labour.push_back((1 - a) * gl);
    G.tfp.push_back(gy - a * gk - (1 - a) * gl);
  }
  const double m = static_cast<double>(n - 1);
  for (std::size_t i = 0; i + 1 < n; i++) {
    G.mean_output += G.output[i] / m;
    G.mean_capital += G.capital[i] / m;
    G.mean_labour += G.labour[i] / m;
    G.mean_tfp += G.tfp[i] / m;
  }
  G.tfp_share = G.mean_output != 0 ? G.mean_tfp / G.mean_output : 0;
  return G;
}

// The capital stock from investment by the perpetual inventory method; `growth` the growth of
// investment used for the first stock (by default the mean log growth of the first five periods).
inline std::vector<double> perpetual_inventory(const std::vector<double>& investment, double depreciation, double growth = NAN) {
  const std::size_t n = investment.size();
  if (n == 0 || !(depreciation > 0 && depreciation < 1)) throw std::invalid_argument("economics::growth::perpetual_inventory: investment, and a depreciation in (0, 1)");
  if (std::isnan(growth)) {
    const std::size_t m = std::min<std::size_t>(5, n - 1);
    growth = m && investment[0] > 0 && investment[m] > 0 ? std::log(investment[m] / investment[0]) / static_cast<double>(m) : 0;
  }
  if (!(growth + depreciation > 0)) throw std::invalid_argument("economics::growth::perpetual_inventory: g + δ is positive");
  std::vector<double> K(n);
  K[0] = investment[0] / (growth + depreciation);
  for (std::size_t t = 1; t < n; t++) K[t] = (1 - depreciation) * K[t - 1] + investment[t - 1];
  return K;
}

}  // namespace economics::growth
