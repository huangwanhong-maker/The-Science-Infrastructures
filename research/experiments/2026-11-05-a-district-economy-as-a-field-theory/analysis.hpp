// What the action's Gaussian measure says of a plan's value (README §4): the posterior uncertainty of the passers at
// the frontages, given the counters (universal kriging's covariance, line by line); the fluctuations of a real week
// about the mean week (a common deviation of the whole district, and a local field with the action's propagator as
// its spatial correlation); and the free energies of the value under them.
//
// For a plan held fixed the value is linear in the passers, V = Σ_{f,h} A_f(h) P_f(h) − costs (economy.hpp), so under
// a Gaussian law of P it is Gaussian, and
//
//   F_T = −T ln E[e^{−V/T}] = E V − Var V / (2T)
//
// exactly: the free energy of the value as a source of the people's field at the temperature T, the certainty
// equivalent of a planner with constant absolute risk aversion 1/T. A free-energy difference ranks two plans by their
// expected values less half their variances (and covariances) over T.

#pragma once

#include "model.hpp"

namespace field {

// The posterior covariance of the frontages' passers at each line (F × F, complex Hermitian), given the counters:
// C_m = σ²(m) [K_FF − K_FO Σ₁⁻¹ K_OF + Q N⁻¹ Q^H],  Q = X_F − K_FO Σ₁⁻¹ X,  N = X^H Σ₁⁻¹ X.
struct Posterior {
  std::size_t F = 0;
  std::vector<std::vector<Complex>> C;  // lines × (F·F)
  std::vector<double> static_correlation;  // F × F: the prior's spatial correlation at the weekly mean (m = 0)
};

inline Posterior posterior(const Kriging& K, const std::vector<std::vector<Complex>>& UA_unused, const std::vector<Observation>& front, const std::vector<double>& sigma2) {
  (void)UA_unused;
  const Graph& G = *K.G;
  Posterior out;
  const std::size_t F = front.size(), n = K.obs.size(), Gr = K.groups, V = G.V;
  out.F = F;
  out.C.assign(lines, std::vector<Complex>(F * F, Complex{0, 0}));
  const auto BF = read_modes(G, front);
  std::vector<std::size_t> all(n);
  std::iota(all.begin(), all.end(), std::size_t{0});
  for (std::size_t m = 0; m < lines; m++) {
    const auto& t = K.t[m];
    std::vector<double> t2(V);
    for (std::size_t k = 0; k < V; k++) t2[k] = std::norm(t[k]);
    std::vector<double> KFF(F * F, 0.0), KFO(F * n, 0.0);
    for (std::size_t f = 0; f < F; f++) {
      for (std::size_t g = f; g < F; g++) {
        double s = 0;
        for (std::size_t k = 0; k < V; k++) s += BF[f][k] * t2[k] * BF[g][k];
        KFF[f * F + g] = KFF[g * F + f] = s;
      }
      for (std::size_t c = 0; c < n; c++) {
        double s = 0;
        for (std::size_t k = 0; k < V; k++) s += BF[f][k] * t2[k] * K.B[c][k];
        KFO[f * n + c] = s;
      }
    }
    if (m == 0) {
      out.static_correlation.assign(F * F, 0.0);
      for (std::size_t f = 0; f < F; f++)
        for (std::size_t g = 0; g < F; g++) out.static_correlation[f * F + g] = KFF[f * F + g] / std::sqrt(KFF[f * F + f] * KFF[g * F + g]);
    }
    std::vector<std::vector<Complex>> XF(F, std::vector<Complex>(Gr, Complex{0, 0}));
    for (std::size_t f = 0; f < F; f++)
      for (std::size_t g = 0; g < Gr; g++) {
        Complex s{0, 0};
        for (std::size_t k = 0; k < V; k++) s += BF[f][k] * t[k] * K.UA[m][g][k];
        XF[f][g] = s;
      }
    const RealSpd S(K.sigma1(m, all));
    // Z = Σ₁⁻¹ K_OF (n × F), W = Σ₁⁻¹ X (n × G).
    std::vector<std::vector<Complex>> Z(F), W(Gr);
    for (std::size_t f = 0; f < F; f++) {
      std::vector<Complex> col(n);
      for (std::size_t c = 0; c < n; c++) col[c] = KFO[f * n + c];
      Z[f] = S.solve(col);
    }
    for (std::size_t g = 0; g < Gr; g++) {
      std::vector<Complex> col(n);
      for (std::size_t c = 0; c < n; c++) col[c] = K.X[m][c][g];
      W[g] = S.solve(col);
    }
    std::vector<std::vector<Complex>> N(Gr, std::vector<Complex>(Gr, Complex{0, 0}));
    for (std::size_t a = 0; a < Gr; a++)
      for (std::size_t b = 0; b < Gr; b++)
        for (std::size_t c = 0; c < n; c++) N[a][b] += std::conj(K.X[m][c][a]) * W[b][c];
    double tr = 0;
    for (std::size_t g = 0; g < Gr; g++) tr += N[g][g].real();
    for (std::size_t g = 0; g < Gr; g++) N[g][g] += 1e-10 * tr / static_cast<double>(Gr);
    // Q = X_F − K_FO W (F × G); then N⁻¹ Q^H column by column.
    std::vector<std::vector<Complex>> Q(F, std::vector<Complex>(Gr, Complex{0, 0}));
    for (std::size_t f = 0; f < F; f++)
      for (std::size_t g = 0; g < Gr; g++) {
        Complex s = XF[f][g];
        for (std::size_t c = 0; c < n; c++) s -= KFO[f * n + c] * W[g][c];
        Q[f][g] = s;
      }
    std::vector<std::vector<Complex>> NQ(F);  // per frontage g: N⁻¹ (row g of Q)^H
    for (std::size_t g = 0; g < F; g++) {
      std::vector<Complex> col(Gr);
      for (std::size_t a = 0; a < Gr; a++) col[a] = std::conj(Q[g][a]);
      NQ[g] = hermitian_solve(N, col);
    }
    for (std::size_t f = 0; f < F; f++)
      for (std::size_t g = 0; g < F; g++) {
        double kz = 0;
        for (std::size_t c = 0; c < n; c++) kz += KFO[f * n + c] * Z[g][c].real();
        Complex q{0, 0};
        for (std::size_t a = 0; a < Gr; a++) q += Q[f][a] * NQ[g][a];
        out.C[m][f * F + g] = sigma2[m] * (KFF[f * F + g] - kz + q);
      }
  }
  return out;
}

// The posterior variance of Σ_{f,h} a_f(h) P⁰_f(h): 168² [a₀ᵀC₀a₀ + a₈₄ᵀC₈₄a₈₄ + 2 Σ_{m=1}^{83} a_m^H C_m a_m].
inline double posterior_variance(const Posterior& Q, const std::vector<double>& a) {
  const std::size_t F = Q.F;
  const auto& W = the_week();
  std::vector<std::vector<Complex>> lines_of(F);
  for (std::size_t f = 0; f < F; f++) lines_of[f] = W.forward(std::span<const double>(a.data() + f * week, week));
  double v = 0;
  for (std::size_t m = 0; m < lines; m++) {
    Complex s{0, 0};
    for (std::size_t f = 0; f < F; f++) {
      if (lines_of[f][m] == Complex{0, 0}) continue;
      Complex row{0, 0};
      for (std::size_t g = 0; g < F; g++) row += Q.C[m][f * F + g] * lines_of[g][m];
      s += std::conj(lines_of[f][m]) * row;
    }
    v += (m == 0 || m == lines - 1 ? 1.0 : 2.0) * s.real();
  }
  return static_cast<double>(week) * static_cast<double>(week) * v;
}

// The fluctuations of a week about the mean week (README §4.2): the district's common relative deviation ε(h), with
// its covariance between the hours of a week over the year's 52 weeks; and a local relative field ζ_f(h) with
// variance T_ζ, the counters' pooled autocorrelation in time, and the action's static correlation in space.
struct Fluctuations {
  std::vector<double> common;   // 168 × 168
  std::vector<double> rho;      // lags 0 … 167
  double T_zeta = 0;
  std::vector<double> relative_variance;  // per counter of the district
};

inline Fluctuations fluctuations(const economy::Json& counters, const std::vector<Counter>& inside) {
  Fluctuations F;
  F.common.assign(week * week, 0.0);
  const auto& Cw = counters["common"]["weekly_covariance"].as_array();
  for (std::size_t a = 0; a < week; a++)
    for (std::size_t b = 0; b < week; b++) {
      const auto& x = Cw[a][b];
      F.common[a * week + b] = x.is_null() ? 0.0 : x.as_number();
    }
  for (const auto& x : counters["local"]["autocorrelation"].as_array()) F.rho.push_back(x.is_null() ? 0.0 : x.as_number());
  F.rho.resize(week, 0.0);
  const auto& var = counters["local"]["variance"].as_array();
  const auto& ids = counters["district_counters"].as_array();
  std::vector<double> rel;
  for (std::size_t i = 0; i < ids.size(); i++) {
    const long long id = static_cast<long long>(ids[i].as_number());
    for (const auto& c : inside)
      if (c.id == id) {
        double m2 = 0;
        for (double y : c.week) m2 += y * y / week;
        rel.push_back(var[i].as_number() / m2);
      }
  }
  F.relative_variance = rel;
  auto s = rel;
  std::sort(s.begin(), s.end());
  F.T_zeta = s.empty() ? 0.0 : s[s.size() / 2];
  return F;
}

// The week-to-week variance of Σ_{f,h} q_f(h) (1 + ε(h) + ζ_f(h)) about its mean: the common part gᵀ C_ε g with
// g(h) = Σ_f q_f(h), and the local part T_ζ Σ_{f,g} ρ_s(f, g) Σ_{h,h'} q_f(h) q_g(h') ρ(|h − h'|).
struct WeekVariance {
  double common = 0, local = 0;
  double total() const { return common + local; }
};

inline WeekVariance week_variance(const Fluctuations& Fl, const Posterior& Q, const std::vector<double>& q) {
  WeekVariance W;
  const std::size_t F = Q.F;
  std::vector<double> g(week, 0.0);
  for (std::size_t f = 0; f < F; f++)
    for (std::size_t h = 0; h < week; h++) g[h] += q[f * week + h];
  for (std::size_t a = 0; a < week; a++)
    for (std::size_t b = 0; b < week; b++) W.common += g[a] * g[b] * Fl.common[a * week + b];
  std::vector<double> conv(F * week, 0.0);
  for (std::size_t f = 0; f < F; f++)
    for (std::size_t h = 0; h < week; h++) {
      const double x = q[f * week + h];
      if (x == 0) continue;
      for (std::size_t k = 0; k < week; k++) conv[f * week + k] += x * Fl.rho[h > k ? h - k : k - h];
    }
  for (std::size_t f = 0; f < F; f++)
    for (std::size_t gg = 0; gg < F; gg++) {
      const double c = Q.static_correlation[f * F + gg];
      if (std::fabs(c) < 1e-6) continue;
      double s = 0;
      for (std::size_t k = 0; k < week; k++) s += conv[f * week + k] * q[gg * week + k];
      W.local += c * s;
    }
  W.local *= Fl.T_zeta;
  return W;
}

}  // namespace field
