// Checks of study (2)'s model of flow (flow.hpp) and its structured search (structured.hpp), against direct solves,
// conservation, closed forms and the shared problem's own value:
//
//   1. walks: Woodbury's update for a set of links equals the direct inverse (I − M)⁻¹; M's spectral radius is below
//      one with every link made, so that the walk sums converge and Z is positive;
//   2. a trip's flow is conserved: at every junction but its ends what flows in flows out, one leaves its origin and one
//      reaches its destination;
//   3. the absorbing destination's closed form equals the flow computed from (I − M⁽ᵈ⁾)⁻¹ with d's row removed;
//   4. the families of trips (arrive, leave, local, visits) equal the sums of their trips' flows, written apart;
//   5. Thomson's principle: a flow driven by Kirchhoff's law through the conductances has no circulation, and the driven
//      part of any flow is the Kirchhoff flow of its divergence;
//   6. the passers' week lines solved as phasors from the sources' lines equal the lines of the hourly passers;
//   7. the value read hour by hour (structured::Engine) equals economy::Problem::evaluate's, as built and for random plans;
//   8. a frontage's best response is exact: the value of the plan changes by exactly what the response predicts.
//
//   check_flow <problem-inputs.json> <problem.json> <sensors.json>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>

#include <calx/random/stream.hpp>

#include "economy.hpp"
#include "flow.hpp"
#include "structured.hpp"

namespace {

int failures = 0, checks = 0;
void check(bool ok, const std::string& what) {
  checks++;
  if (!ok) failures++;
  std::printf("%s  %s\n", ok ? "ok  " : "FAIL", what.c_str());
}
std::string sci(double x) {
  char b[32];
  std::snprintf(b, sizeof b, "%.2e", x);
  return b;
}
std::string read_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read " + path);
  std::stringstream b;
  b << in.rdbuf();
  return b.str();
}

}  // namespace

int main(int argc, char** argv) {
  using namespace economy;
  using flow::none;
  namespace la = calx::linalg;
  constexpr std::size_t H = economy::week;
  if (argc < 4) {
    std::cerr << "usage: check_flow <problem-inputs.json> <problem.json> <sensors.json>\n";
    return 2;
  }
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  const Json in = Json::parse(read_file(argv[1]));
  const Json pj = Json::parse(read_file(argv[2]));
  const Json sj = Json::parse(read_file(argv[3]));
  const Problem P(in, &pj);
  flow::Options opt;
  flow::Model M(P, in, sj, opt);
  const auto& N = M.net;
  const std::size_t n = N.nodes(), L = P.links.size(), F = P.size();
  std::printf("network: %zu junctions, %zu segments (%zu of links), %zu entrances, %zu counters inside, theta %.3f/m\n", n, N.edges.size(),
              static_cast<std::size_t>(std::count_if(N.edges.begin(), N.edges.end(), [](const auto& e) { return e.link >= 0; })), M.entrances.size(), M.counters.size(), M.o.theta);
  calx::random::Stream s(20261105, 2, 0, 0);

  // 1. Walks.
  const std::size_t S = N.states();
  {
    check(M.radius_all < 1, "spectral radius of the transitions with every link made below one (" + sci(M.radius_all) + ")");
    double worst = 0;
    for (int t = 0; t < 6; t++) {
      hv::Bits made(L, 0);
      for (auto& b : made) b = t == 0 ? 1 : (s.uniform() < 0.5);
      const auto A = M.walks(made);
      const auto B = flow::walks_direct(N, M.omega, made);
      double scale = 0;
      for (double z : B.z) scale = std::max(scale, std::fabs(z));
      for (std::size_t i = 0; i < S * S; i++) worst = std::max(worst, std::fabs(A.z[i] - B.z[i]) / scale);
    }
    check(worst < 1e-10, "Woodbury's update equals the direct inverse (largest difference relative to the largest entry " + sci(worst) + ")");
    bool positive = true;
    std::size_t zeros = 0;
    const auto A = M.walks(hv::Bits(L, 1));
    for (double z : A.z) positive = positive && z > -1e-12, zeros += z == 0;
    check(positive, "no walk sum negative with every link made (" + std::to_string(zeros) + " exactly zero: segments no walk joins)");
  }
  // The junctions of the largest component.
  std::vector<std::size_t> main_nodes;
  for (std::size_t j = 0; j < n; j++)
    if (M.component[j] == 0) main_nodes.push_back(j);
  // 2 and 3. Conservation, and the absorbing destination.
  {
    const hv::Bits made(L, 1);
    const auto T = M.trips(made);
    double worst = 0, worst_abs = 0;
    for (int t = 0; t < 8; t++) {
      const std::size_t o = main_nodes[s.below(main_nodes.size())], d = main_nodes[s.below(main_nodes.size())];
      if (o == d) continue;
      std::vector<double> net(n, 0.0);  // outflow − inflow at each junction
      std::vector<double> f(S, 0.0);
      for (std::size_t c = 0; c < S; c++) {
        if (!N.present(c / 2, made)) continue;
        f[c] = T.visits(o, d, c);
        net[N.tail(c)] += f[c], net[N.head(c)] -= f[c];
      }
      for (std::size_t j = 0; j < n; j++) worst = std::max(worst, std::fabs(net[j] - (j == o ? 1.0 : j == d ? -1.0 : 0.0)));
      // d absorbing, directly: (I − M^(d))⁻¹ with the rows of the segments entering d removed.
      la::Matrix A = flow::transitions(N, M.omega, made);
      for (auto a : T.entering[d])
        for (std::size_t j = 0; j < S; j++) A(a, j) = 0;
      for (std::size_t i = 0; i < S; i++)
        for (std::size_t j = 0; j < S; j++) A(i, j) = (i == j ? 1.0 : 0.0) - A(i, j);
      const auto Zd = la::solve(A, la::Matrix::identity(S));
      std::vector<double> pi(S, 0.0), sigma(S, 0.0);
      for (auto b : T.leaving[o])
        for (std::size_t c = 0; c < S; c++) pi[c] += M.omega[b / 2] * Zd(b, c);
      for (std::size_t c = 0; c < S; c++)
        for (auto a : T.entering[d]) sigma[c] += Zd(c, a);
      double W = 0;
      for (auto a : T.entering[d]) W += pi[a];
      for (std::size_t c = 0; c < S; c++)
        if (N.present(c / 2, made)) worst_abs = std::max(worst_abs, std::fabs(pi[c] * sigma[c] / W - f[c]));
    }
    check(worst < 1e-9, "a trip is conserved at every junction (largest imbalance " + sci(worst) + ")");
    check(worst_abs < 1e-9, "the absorbing destination by Woodbury equals the direct solve (largest difference " + sci(worst_abs) + ")");
  }
  // 4. The families against their trips.
  {
    const hv::Bits made(L, 0);
    const auto T = M.trips(made);
    std::vector<std::size_t> edges;
    for (int t = 0; t < 6; t++) edges.push_back(s.below(N.edges.size()));
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    edges.erase(std::remove_if(edges.begin(), edges.end(), [&](std::size_t e) { return !N.present(e, made); }), edges.end());
    const auto R = M.terms(T, edges);
    double sw = 0;
    for (double w : M.work) sw += w;
    double worst = 0;
    for (std::size_t k = 0; k < edges.size(); k++) {
      const std::size_t c = 2 * edges[k];
      // arrive (work) from the first entrance, and leave (work) to it.
      const std::size_t o = M.entrances[0].node;
      double a = 0, l = 0;
      for (std::size_t b = 0; b < n; b++) {
        if (M.work[b] <= 0 || b == o) continue;
        a += M.work[b] / sw * T.visits(o, b, c);
        l += M.work[b] / sw * T.visits(b, o, c);
      }
      worst = std::max({worst, std::fabs(a - R.arrive_leave[0][0][k][0]), std::fabs(l - R.arrive_leave[2][0][k][0])});
      // visits to the first frontage and back.
      double v = 0;
      for (const auto& [j, wgt] : M.frontage_nodes[0])
        for (std::size_t oo = 0; oo < n; oo++)
          if (M.origin[oo] > 0 && oo != j) v += wgt * M.origin[oo] * (T.visits(oo, j, c) + T.visits(j, oo, c));
      worst = std::max(worst, std::fabs(v - R.frontage[0][k][0]));
    }
    check(worst < 1e-9, "arrive, leave and visits equal the sums of their trips (largest difference " + sci(worst) + ")");
  }
  // 5. Thomson's principle, on the network as built.
  {
    const hv::Bits made(L, 0);
    std::vector<std::size_t> present;
    for (std::size_t e = 0; e < N.edges.size(); e++)
      if (N.present(e, made)) present.push_back(e);
    // Kirchhoff's flow for a source at one entrance and a sink at another, through the conductances w/ℓ.
    la::Matrix Lp(n);
    for (auto e : present) {
      const auto& x = N.edges[e];
      const double c = x.width / x.length;
      Lp(x.u, x.u) += c, Lp(x.v, x.v) += c, Lp(x.u, x.v) -= c, Lp(x.v, x.u) -= c;
    }
    const std::size_t a = M.entrances[0].node, b = M.entrances[1].node;
    // Ground b, and one junction of every other component of the network (it has several: lanes on their own).
    std::vector<std::size_t> comp(n, none);
    std::vector<std::vector<std::size_t>> adj(n);
    for (auto e : present) adj[N.edges[e].u].push_back(N.edges[e].v), adj[N.edges[e].v].push_back(N.edges[e].u);
    std::size_t comps = 0;
    std::vector<std::size_t> roots{b};
    for (std::size_t j = 0; j < n; j++) roots.push_back(j);
    for (const std::size_t root : roots) {
      if (comp[root] != none) continue;
      std::vector<std::size_t> stack{root};
      comp[root] = comps;
      while (!stack.empty()) {
        const auto j = stack.back();
        stack.pop_back();
        for (auto k : adj[j])
          if (comp[k] == none) comp[k] = comps, stack.push_back(k);
      }
      for (std::size_t j = 0; j < n; j++) Lp(root, j) = 0;
      Lp(root, root) = 1;
      comps++;
    }
    std::printf("      the network as built has %zu components\n", comps);
    la::Matrix rhs(n);
    rhs(a, 0) = 1;
    const auto phi = la::solve(Lp, rhs);
    std::vector<std::array<double, 2>> fl(N.edges.size(), {0.0, 0.0});
    for (auto e : present) {
      const auto& x = N.edges[e];
      const double q = x.width / x.length * (phi(x.u, 0) - phi(x.v, 0));
      fl[e] = {std::max(q, 0.0), std::max(-q, 0.0)};
    }
    const auto C = flow::circulation(N, made, fl);
    check(C.parts.local + C.parts.global < 1e-10, "a Kirchhoff flow has no circulation (local and global shares " + sci(C.parts.local + C.parts.global) + ")");
    // A walkers' flow: its driven part is the Kirchhoff flow of its divergence.
    const auto Tr = M.trips(made);
    std::vector<std::array<double, 2>> wf(N.edges.size(), {0.0, 0.0});
    for (auto e : present) wf[e] = {Tr.visits(a, b, 2 * e), Tr.visits(a, b, 2 * e + 1)};
    const auto W = flow::circulation(N, made, wf);
    double worst = 0, scale = 0;
    for (std::size_t k = 0; k < W.edges.size(); k++) {
      const auto& x = N.edges[W.edges[k]];
      const double q = x.width / x.length * (phi(x.u, 0) - phi(x.v, 0));
      worst = std::max(worst, std::fabs(W.parts.parts.driven[k] - q)), scale = std::max(scale, std::fabs(q));
    }
    check(worst / scale < 1e-8, "the driven part of the walkers' flow is the Kirchhoff flow of its divergence (relative " + sci(worst / scale) + ")");
    std::printf("      one trip between two entrances: driven %.3f, round the blocks %.3f, round the holes %.3f of the energy; %zu blocks, %zu holes\n", W.parts.driven,
                W.parts.local, W.parts.global, W.faces.size(), W.holes);
  }
  // Calibrate (the next checks need passers).
  {
    std::vector<std::size_t> edges;
    for (const auto& c : M.counters) edges.push_back(c.edge);
    for (auto e : M.fronted_edges()) edges.push_back(e);
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    const auto T = M.terms(M.trips(hv::Bits(L, 0)), edges);
    const auto fit = M.calibrate(T);
    std::printf("      calibration at theta %.3f: R2 %.3f (leave one counter out %.3f)\n", M.o.theta, fit.r2, fit.loo_r2);
  }
  // 6. Phasors.
  {
    const Plan built = P.as_built();
    const auto R = M.response(hv::Bits(L, 0));
    const auto pass = M.passers(built, R);
    // The sources' lines: P_f(h) = B_f(h) + Σ_g c(h) draw_g(h) R_fg; B is itself linear in its amplitudes' series.
    double worst = 0;
    const auto draw = M.frontage_draw(built);
    for (std::size_t m : {1ul, 7ul, 14ul}) {
      std::vector<std::complex<double>> src(F);
      for (std::size_t g = 0; g < F; g++) {
        std::vector<double> sg(H);
        for (std::size_t h = 0; h < H; h++) sg[h] = M.x[flow::visit][h] / M.W_ref * draw[g * H + h];
        src[g] = flow::line(sg, m);
      }
      for (std::size_t f = 0; f < F; f++) {
        std::complex<double> ph = flow::line(std::span<const double>(R.B.data() + f * H, H), m);
        for (std::size_t g = 0; g < F; g++) ph += src[g] * R.R[f * F + g];
        const auto direct = flow::line(std::span<const double>(pass.data() + f * H, H), m);
        worst = std::max(worst, std::abs(ph - direct) / std::max(1.0, std::abs(direct)));
      }
    }
    check(worst < 1e-10, "the passers' lines solved as phasors equal the lines of the hourly passers (relative " + sci(worst) + ")");
  }
  // 7 and 8. The value hour by hour, and the best response's exactness.
  {
    const structured::Engine E(P, M);
    const economy::Kernels K{&P};
    double worst = 0, worst_br = 0;
    for (int t = 0; t < 20; t++) {
      Plan p = t == 0 ? P.as_built() : K.sample(s);
      if (t % 3 == 1) p.links.assign(L, 0);
      const auto R = M.response(p.links);
      const auto o = P.evaluate(p, M.passers(p, R));
      const double v = E.value(p, R);
      worst = std::max(worst, std::fabs(o.value - v) / std::max(1.0, std::fabs(o.value)));
      if (t < 8) {
        const std::size_t f = s.below(F);
        const auto G = E.gains(p, R, f);
        const auto [best, per] = E.choose(p, f, G);
        (void)per;
        const double predicted = best.value - E.present(p, f, G);
        Plan q = p;
        q.kind[f] = best.kind;
        for (std::size_t j = 0; j < 4; j++) q.hours[4 * f + j] = best.hours[j];
        const double actual = P.evaluate(q, M.passers(q, M.response(q.links))).value - o.value;
        worst_br = std::max(worst_br, std::fabs(predicted - actual) / std::max(1.0, std::fabs(actual)));
      }
    }
    check(worst < 1e-10, "the value hour by hour equals the problem's evaluate (relative " + sci(worst) + ")");
    check(worst_br < 1e-9, "a best response changes the value by exactly what it predicts (relative " + sci(worst_br) + ")");
  }
  std::printf("%d checks, %d failed\n", checks, failures);
  return failures ? 1 : 0;
}
