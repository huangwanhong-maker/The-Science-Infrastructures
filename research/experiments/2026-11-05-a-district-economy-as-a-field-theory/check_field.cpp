// Checks of the study's field theory (README §7), against closed forms, direct solves, physis's own propagator and the
// shared problem's evaluate():
//
//   1. the week's modes: the inverse transform undoes the forward one;
//   2. the stationary path in the modes solves the conservation law: (iω⋆₀ + L + μ⋆₀) φ̂ = ŝ, to rounding;
//   3. with unit stars on the domain's cells, the mode sum is physis::field_theory::free_propagator (Δ + m²)⁻¹ on the
//      same complex; the subdivided lattice's chain complex is valid and has the domain's homology;
//   4. Noether: the number of people, Σ ⋆₀ φ̂, obeys (iω + μ) N̂ = Σ ŝ at every line;
//   5. Hodge: Δ₀ has one zero mode per component, with a constant eigenvector; Δ₁'s harmonic dimension is
//      b₁ = E − V + b₀ − F (the district's 2-cells filled), and Δ₂ has none;
//   6. kriging with no nugget interpolates the counters exactly, and the posterior mean of the sources, carried by the
//      stationary path, gives back the kriging predictor at every counter;
//   7. the attraction's response R(h) equals a direct solve of (L + μ⋆₀) δφ = b L_φ̄ W δp;
//   8. a new edge's exact effect by Sherman–Morrison equals a direct solve with the edge added; and the model's links,
//      by Woodbury's identity, equal a direct solve with two links' Whitney edges added;
//   9. the plan field's value equals the problem's evaluate() on the field's passers, for the plan as built and random
//      plans; after every heat-bath update (at T > 0 and T = 0) it still equals a rebuild and evaluate(); and an update
//      at T = 0 never lowers the value.
//
//   check_field <problem-inputs.json> <problem.json> <sensors.json> <counters.json>

#include <cstdio>

#include <morphe/topology/regular_complex.hpp>
#include <physis/field_theory/scalar.hpp>

#include "context.hpp"
#include "plan.hpp"

namespace {

using namespace field;
int failures = 0, checks = 0;
void check(bool ok, const std::string& what, double detail = 0) {
  checks++;
  std::printf("%s  %s (%.3g)\n", ok ? "ok  " : "FAIL", what.c_str(), detail);
  if (!ok) failures++;
}

// (iω⋆₀ + D L₁ + μ⋆₀) φ, L₁ = d₀ᵀ diag(c/ℓ) d₀.
std::vector<Complex> apply(const Graph& G, const Physics& ph, std::size_t m, const std::vector<Complex>& phi) {
  std::vector<Complex> out(G.V, Complex{0, 0});
  const double D = ph.diffusivity();
  for (std::size_t v = 0; v < G.V; v++) out[v] = Complex(ph.leave, omega(m)) * G.star0[v] * phi[v];
  for (std::size_t e = 0; e < G.E; e++) {
    const auto a = G.ends[e][0], b = G.ends[e][1];
    const Complex j = D * G.conductance[e] * (phi[a] - phi[b]);
    out[a] += j, out[b] -= j;
  }
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 5) {
    std::fprintf(stderr, "usage: check_field <problem-inputs.json> <problem.json> <sensors.json> <counters.json>\n");
    return 2;
  }
  Context C = load(argv[1], argv[2], argv[3], argv[4]);
  const Problem& P = *C.problem;
  Physics ph;
  ph.length = 45, ph.little = 1, ph.lane = 1, ph.nugget = 0.3, ph.relative = true;
  Model M(P, ph);
  const Graph& G = M.base;
  calx::random::Stream z(20261105, 3, 0, 0);

  // 1. The week's modes.
  {
    std::vector<double> x(week);
    for (auto& v : x) v = z.normal();
    const auto y = the_week().inverse(the_week().forward(x));
    double e = 0;
    for (std::size_t h = 0; h < week; h++) e = std::max(e, std::fabs(y[h] - x[h]));
    check(e < 1e-12, "the week's inverse transform undoes the forward one", e);
  }
  // 2. The stationary path solves the conservation law; 4. Noether.
  {
    double worst = 0, noether = 0;
    for (std::size_t m : {0, 1, 7, 14, 42, 84}) {
      std::vector<Complex> s(G.V);
      for (auto& v : s) v = Complex(z.normal(), m ? z.normal() : 0.0);
      const auto phi = respond(G, transfer(G, ph, m), s);
      const auto back = apply(G, ph, m, phi);
      double e = 0, n = 0;
      for (std::size_t v = 0; v < G.V; v++) e = std::max(e, std::abs(back[v] - s[v])), n = std::max(n, std::abs(s[v]));
      worst = std::max(worst, e / n);
      Complex N{0, 0}, S{0, 0};
      for (std::size_t v = 0; v < G.V; v++) N += G.star0[v] * phi[v], S += s[v];
      noether = std::max(noether, std::abs(Complex(ph.leave, omega(m)) * N - S) / std::abs(S));
    }
    check(worst < 1e-9, "the mode sum solves (iω⋆₀ + L + μ⋆₀) φ = s", worst);
    check(noether < 1e-10, "Noether: (iω + μ) Σ⋆₀φ̂ = Σŝ at every line", noether);
  }
  // 3. Unit stars: physis's free propagator.
  {
    const double m2 = 0.37;
    const auto X = morphe::topology::cw(P.base.complex);
    const auto Gf = physis::field_theory::free_propagator(X, m2);
    const Graph G0 = graph(P.base, 1, 1, 0);  // the domain's cells, not subdivided
    const auto H = physis::discretisation::hodge_laplacian(G0.chain, 0, {});
    double e = 0;
    for (std::size_t a = 0; a < G0.V; a += 7)
      for (std::size_t b = 0; b < G0.V; b++) {
        double s = 0;
        for (std::size_t k = 0; k < G0.V; k++) s += H.vectors(a, k) * H.vectors(b, k) / (H.values[k] + m2);
        e = std::max(e, std::fabs(s - Gf(a, b)));
      }
    check(e < 1e-9, "with unit stars the mode sum is physis::field_theory::free_propagator", e);
    std::printf("      the lattice: %zu vertices, %zu pieces (the domain's %zu vertices, %zu segments), %zu blocks\n", G.V, G.E, G0.V, G0.E, G.chain.n[2]);
  }
  // 5. Hodge: the components by a search of the graph, then the kernels.
  std::size_t b0 = 0;
  {
    std::vector<int> comp(G.V, -1);
    for (std::size_t s = 0; s < G.V; s++) {
      if (comp[s] >= 0) continue;
      std::vector<std::size_t> stack{s};
      comp[s] = static_cast<int>(b0);
      while (!stack.empty()) {
        const auto v = stack.back();
        stack.pop_back();
        for (const auto& [e, u] : G.at[v])
          if (comp[u] < 0) comp[u] = static_cast<int>(b0), stack.push_back(u);
      }
      b0++;
    }
    const auto zero = static_cast<std::size_t>(std::count_if(G.lambda.begin(), G.lambda.end(), [&](double x) { return x < 1e-9 * G.lambda.back(); }));
    // Each zero mode is constant on every component (a combination of their indicators).
    double spread = 0;
    for (std::size_t k = 0; k < zero; k++)
      for (std::size_t c = 0; c < b0; c++) {
        double lo = 1e300, hi = -1e300;
        for (std::size_t v = 0; v < G.V; v++)
          if (comp[v] == static_cast<int>(c)) lo = std::min(lo, G.U(v, k)), hi = std::max(hi, G.U(v, k));
        spread = std::max(spread, hi - lo);
      }
    check(zero == b0 && spread < 1e-8, "Δ₀ has one zero mode per component (" + std::to_string(b0) + "), constant on each: the numbers of people", spread);
    std::size_t F = G.chain.n.size() > 2 ? G.chain.n[2] : 0;
    std::vector<double> w2(F, 1.0);
    const auto H1 = physis::discretisation::hodge_laplacian(G.chain, 1, {G.star0, G.conductance, w2});
    const auto H2 = physis::discretisation::hodge_laplacian(G.chain, 2, {G.star0, G.conductance, w2});
    const long long b1 = static_cast<long long>(G.E) - static_cast<long long>(G.V) + static_cast<long long>(b0) - static_cast<long long>(F);
    check(static_cast<long long>(H1.harmonic(1e-9)) == b1 && H2.harmonic(1e-9) == 0, "Δ₁'s harmonic dimension is E − V + b₀ − F, Δ₂ has none", static_cast<double>(b1));
  }
  // 6. Kriging: exact interpolation with no nugget, and the sources' posterior mean carried back.
  {
    Physics p0 = ph;
    p0.nugget = 1e-9, p0.relative = false;
    std::vector<Observation> obs;
    std::vector<std::vector<Complex>> y;
    for (const auto& c : C.inside) obs.push_back(observe(G, counter_foot(G, c))), y.push_back(the_week().forward(c.week));
    // Counters at one site would make the system singular without a nugget: keep one counter a site.
    const auto sites = sites_of(C.inside);
    std::vector<Observation> o1;
    std::vector<std::vector<Complex>> y1;
    for (const auto& s : sites) o1.push_back(obs[s[0]]), y1.push_back(y[s[0]]);
    const auto gr = source_groups(P, G, C.inputs, C.outside, false, false);
    const Kriging K(G, p0, o1, y1, gr);
    std::vector<std::size_t> all(o1.size());
    std::iota(all.begin(), all.end(), std::size_t{0});
    double interp = 0, carried = 0;
    for (std::size_t m : {0, 7, 14, 21}) {
      const auto F = K.fit(m, all);
      const auto s = K.source(m, F, gr);
      const auto phi = respond(G, transfer(G, p0, m), s);
      for (std::size_t q = 0; q < o1.size(); q++) {
        const Complex pr = K.predict(m, F, q);
        const Complex at = o1[q].wa * phi[o1[q].a] + o1[q].wb * phi[o1[q].b];
        interp = std::max(interp, std::abs(pr - y1[q][m]) / std::abs(y1[q][0]));
        carried = std::max(carried, std::abs(at - pr) / std::abs(y1[q][0]));
      }
    }
    check(interp < 1e-5, "kriging with no nugget interpolates the counters", interp);
    check(carried < 1e-8, "the sources' posterior mean, carried by the stationary path, is the kriging predictor", carried);
  }
  // A calibrated model at the chosen parameters, for the rest.
  {
    std::vector<Observation> obs;
    std::vector<std::vector<Complex>> y;
    for (const auto& c : C.inside) obs.push_back(observe(G, counter_foot(G, c))), y.push_back(the_week().forward(c.week));
    const auto gr = source_groups(P, G, C.inputs, C.outside, false, false);
    const Kriging K(G, ph, obs, y, gr);
    std::vector<std::size_t> all(obs.size());
    std::iota(all.begin(), all.end(), std::size_t{0});
    M.source.assign(lines, {}), M.sigma2.assign(lines, 0.0);
    for (std::size_t m = 0; m < lines; m++) {
      const auto F = K.fit(m, all);
      M.source[m] = K.source(m, F, gr), M.sigma2[m] = F.sigma2;
    }
    M.p0 = M.chances(P.as_built(), *P.distances(economy::hv::Bits(P.links.size(), 0)));
    M.prepare();
  }
  // 7. The attraction's response against a direct solve.
  {
    const auto L = M.build(economy::hv::Bits(P.links.size(), 0), true);
    const double D = ph.diffusivity();
    double worst = 0;
    for (std::size_t trial = 0; trial < 4; trial++) {
      const std::size_t h = 8 + 37 * trial % 150, g = z.below(P.size());
      // b L_φ̄ W_g, then (L + μ⋆₀) δφ = it, by calx's dense solve.
      std::vector<double> y(G.V, 0.0);
      const auto& o = L.front[g];
      for (const auto& [v, w] : std::vector<std::pair<std::size_t, double>>{{o.a, o.wa}, {o.b, o.wb}})
        for (const auto& [e, u] : G.at[v]) {
          const double c = std::max(0.0, D * G.conductance[e] * 0.5 * (L.phi[G.ends[e][0] * week + h] + L.phi[G.ends[e][1] * week + h]));
          y[v] += ph.attraction * c * w, y[u] -= ph.attraction * c * w;
        }
      calx::linalg::Matrix A(G.V), Bm(G.V);
      for (std::size_t v = 0; v < G.V; v++) A(v, v) += ph.leave * G.star0[v], Bm(v, 0) = y[v];
      for (std::size_t e = 0; e < G.E; e++) {
        const auto a = G.ends[e][0], b = G.ends[e][1];
        const double c = D * G.conductance[e];
        A(a, a) += c, A(b, b) += c, A(a, b) -= c, A(b, a) -= c;
      }
      const auto X = calx::linalg::solve(A, Bm);
      for (std::size_t f = 0; f < P.size(); f++) {
        const auto& of = L.front[f];
        const double direct = of.wa * X(of.a, 0) + of.wb * X(of.b, 0);
        const double viaR = M.R[(h * P.size() + f) * P.size() + g];
        worst = std::max(worst, std::fabs(direct - viaR) / (std::fabs(M.R[(h * P.size() + g) * P.size() + g]) + 1e-300));
      }
    }
    check(worst < 1e-8, "the attraction's response R(h) equals a direct solve", worst);
  }
  // 8. Sherman–Morrison for a new edge, at the weekly mean (m = 0).
  {
    const double D = ph.diffusivity();
    const std::size_t a = z.below(G.V), b = (a + 101) % G.V;
    const double kappa = D * 1.0 / 40.0;
    std::vector<Complex> s(G.V);
    for (auto& v : s) v = z.normal();
    const auto t = transfer(G, ph, 0);
    const auto phi = respond(G, t, s);
    std::vector<Complex> bv(G.V, Complex{0, 0});
    bv[a] = 1, bv[b] = -1;
    const auto Gb = respond(G, t, bv);
    const Complex bphi = phi[a] - phi[b], bGb = Gb[a] - Gb[b];
    calx::linalg::Matrix A(G.V), Bm(G.V);
    for (std::size_t v = 0; v < G.V; v++) A(v, v) += ph.leave * G.star0[v], Bm(v, 0) = s[v].real();
    for (std::size_t e = 0; e < G.E; e++) {
      const auto u = G.ends[e][0], w = G.ends[e][1];
      const double c = D * G.conductance[e];
      A(u, u) += c, A(w, w) += c, A(u, w) -= c, A(w, u) -= c;
    }
    A(a, a) += kappa, A(b, b) += kappa, A(a, b) -= kappa, A(b, a) -= kappa;
    const auto X = calx::linalg::solve(A, Bm);
    double e = 0, n = 0;
    for (std::size_t v = 0; v < G.V; v++) {
      const Complex sm = phi[v] - kappa * Gb[v] * bphi / (1.0 + kappa * bGb);
      e = std::max(e, std::fabs(sm.real() - X(v, 0))), n = std::max(n, std::fabs(X(v, 0)));
    }
    check(e / n < 1e-9, "Sherman–Morrison: a new edge's exact effect equals a direct solve", e / n);
  }
  // 8b. The model's links, by Woodbury, against a direct solve with two links' Whitney edges added (line m = 0).
  {
    const std::vector<std::size_t> S{2, 8};
    calx::linalg::Matrix A(G.V), Bm(G.V);
    const double D = ph.diffusivity();
    for (std::size_t v = 0; v < G.V; v++) A(v, v) += ph.leave * G.star0[v], Bm(v, 0) = M.source[0][v].real();
    for (std::size_t e = 0; e < G.E; e++) {
      const auto u = G.ends[e][0], w = G.ends[e][1];
      const double c = D * G.conductance[e];
      A(u, u) += c, A(w, w) += c, A(u, w) -= c, A(w, u) -= c;
    }
    for (auto j : S) {
      std::vector<std::pair<std::size_t, double>> b{{M.link_a[j].a, M.link_a[j].wa}, {M.link_a[j].b, M.link_a[j].wb}, {M.link_b[j].a, -M.link_b[j].wa}, {M.link_b[j].b, -M.link_b[j].wb}};
      for (const auto& [x, wx] : b)
        for (const auto& [y, wy] : b) A(x, y) += M.link_kappa[j] * wx * wy;
    }
    const auto X = calx::linalg::solve(A, Bm);
    const std::size_t L = P.links.size();
    std::vector<std::vector<Complex>> H2(2, std::vector<Complex>(2));
    std::vector<Complex> b2(2);
    for (std::size_t i = 0; i < 2; i++) {
      for (std::size_t j = 0; j < 2; j++) H2[i][j] = M.H[0][S[i] * L + S[j]] + (i == j ? 1.0 / M.link_kappa[S[i]] : 0.0);
      b2[i] = M.Bphi[0][S[i]];
    }
    const auto c = complex_solve(H2, b2);
    double e = 0, n = 0;
    for (std::size_t f = 0; f < P.size(); f++) {
      Complex x = M.P0hat[0][f];
      for (std::size_t i = 0; i < 2; i++) x -= M.TF[0][f * L + S[i]] * c[i];
      const auto& o = M.front_base[f];
      const double direct = o.wa * X(o.a, 0) + o.wb * X(o.b, 0);
      e = std::max(e, std::fabs(x.real() - direct)), n = std::max(n, std::fabs(direct));
    }
    check(e / n < 1e-9, "Woodbury: two links made equal a direct solve with their Whitney edges", e / n);
  }
  // 9. The plan field against the problem's evaluate().
  {
    economy::Kernels K{&P};
    double worst = 0;
    std::vector<Plan> plans{P.as_built()};
    for (int i = 0; i < 12; i++) plans.push_back(K.sample(z));
    for (const auto& p : plans) {
      const PlanState S(M, p);
      const double v = M.evaluate(p).value;
      worst = std::max(worst, std::fabs(S.value - v) / std::max(1.0, std::fabs(v)));
    }
    check(worst < 1e-9, "the plan field's value equals evaluate() on the field's passers", worst);
    double drift = 0, against = 0;
    bool never_lower = true;
    LocalWork W;
    for (int trial = 0; trial < 6; trial++) {
      PlanState S(M, trial == 0 ? P.as_built() : K.sample(z));
      for (int u = 0; u < 25; u++) {
        const std::size_t f = z.below(P.size());
        const double T = u % 2 ? 0.0 : 500.0;
        const double before = S.value;
        S.update(f, T, z, W);
        if (T == 0 && S.value < before - 1e-6) never_lower = false;
      }
      const double kept = S.value;
      std::size_t clamped = 0;
      M.passers(S.plan, &clamped);
      std::printf("      trial %d: %zu of %zu passers clamped at zero\n", trial, clamped, P.size() * week);
      PlanState R(M, S.plan);
      const double ev = M.evaluate(S.plan).value;
      drift = std::max(drift, std::fabs(kept - R.value) / std::max(1.0, std::fabs(R.value)));
      against = std::max(against, std::fabs(kept - ev) / std::max(1.0, std::fabs(ev)));
    }
    check(drift < 1e-8, "after heat-bath updates the incremental value equals a rebuild", drift);
    check(against < 1e-6, "after heat-bath updates the incremental value equals evaluate()", against);
    check(never_lower, "an update at T = 0 never lowers the value");
  }
  std::printf("%d checks, %d failed\n", checks, failures);
  return failures ? 1 : 0;
}
