// The people's field composed with the shared problem (README §2.4, §4): the calibrated sources carried to the
// network of any set of links, the passers at every frontage and hour, and their response to the plan's attraction.
//
// For a plan with links L, the passers are
//
//   P_f(h) = P⁰_f(h; L) + Σ_g R_fg(h) (p_g(h) − p⁰_g(h)),
//
// P⁰ the stationary path of the calibrated sources on the lattice with the links L made (README §2.4): each link a
// new 1-cell of a lane's conductance κ_j = D c_lane / ℓ_j between its ends' Whitney points, b_j = W(a_j) − W(b_j), so
// that the operator gains B K Bᵀ and, by Woodbury's identity, line by line,
//
//   φ̂' = φ̂ − G B (K⁻¹ + Bᵀ G B)⁻¹ Bᵀ φ̂,   G = U diag(t) Uᵀ,
//
// exact on the lattice and a k × k solve for k links; p_g(h) the plan's chance to
// stop at g (the problem's, with its catalysis), p⁰ the plan as built's, and R the static linear response of the
// field to the attraction −b p_g of each open frontage (the drift current −⋆₁ φ̄ d₀U, README §2.2):
//
//   R(h) = b O_F G₀ L_{φ̄(h)} W,   G₀ = U diag(1/(Dλ + μ)) Uᵀ,   L_φ̄ = d₀ᵀ diag(⋆₁ φ̄_e) d₀,
//
// with W the frontages' Whitney weights and φ̄_e the edge's mean of the field at hour h, computed once on the network
// as built: a link changes the response only near itself, and the response is itself a few per cent of the passers, so
// the error is of second order (README §2.4). The response is quasi-static:
// the district relaxes at a rate of at least μ = 6/h, so an hour's attraction acts within that hour (the error is of
// order (ω/μ)², 0.2% at the daily line). The passers are kept non-negative.

#pragma once

#include <chrono>
#include <memory>
#include <mutex>

#include "field.hpp"

namespace field {

using economy::Plan;
using economy::Problem;

struct Counter {
  std::size_t row = 0;  // in sensors.json
  long long id = 0;
  std::string name;
  Point at{}, foot{};
  std::vector<double> week;  // the mean week, people an hour
};

// The base vertices' patterns of the source groups (README §3.1): the area's edges (each fed, if asked, by the mean
// week of the nearest counter outside the area, so that the city around the district is its boundary condition),
// the buildings' volume (each parcel's and neighbour's footprint × height on its nearest segment), and the streets'
// and the lanes' length (⋆₀ by kind of way). Static patterns are the same at every line.
inline Groups source_groups(const Problem& P, const Graph& G, const economy::Json& inputs, const std::vector<Counter>& outside, bool edge_from_counters, bool trams) {
  Groups gr;
  const auto& box = inputs["district"]["district_box"];
  const double x0 = box["x0"].as_number(), x1 = box["x1"].as_number(), y0 = box["y0"].as_number(), y1 = box["y1"].as_number();
  std::vector<double> build(G.V, 0.0), main(G.V, 0.0), little(G.V, 0.0), lane(G.V, 0.0), tram(G.V, 0.0);
  // The streets of the district that carry trams (Public Transport Victoria's network map): an assumption read from
  // the public map, not from the data library.
  static const std::set<std::string> tram_streets{"Swanston Street", "Elizabeth Street", "Collins Street", "Bourke Street", "Flinders Street"};
  std::vector<double> edge_weight(G.V, 0.0);
  std::vector<std::size_t> boundary;
  for (std::size_t v = 0; v < G.V; v++) {
    const auto& p = G.pos[v];
    const double d = std::min({std::fabs(p[0] - x0), std::fabs(p[0] - x1), std::fabs(p[1] - y0), std::fabs(p[1] - y1)});
    if (G.at[v].size() == 1 && d < 3) boundary.push_back(v), edge_weight[v] = G.conductance[G.at[v][0].first] * G.length[G.at[v][0].first];
  }
  for (std::size_t e = 0; e < G.E; e++) {
    auto& t = G.cls[e] == main_street ? main : G.cls[e] == little_street ? little : lane;
    t[G.ends[e][0]] += 0.5 * G.length[e], t[G.ends[e][1]] += 0.5 * G.length[e];
    if (G.cls[e] == main_street && tram_streets.count(G.name[e])) tram[G.ends[e][0]] += 0.5 * G.length[e], tram[G.ends[e][1]] += 0.5 * G.length[e];
  }
  const auto add_building = [&](const economy::Ring& ring, double height) {
    const double area = std::fabs(cg::ring_area(ring));
    if (!(area > 0) || !(height > 0)) return;
    double cx = 0, cy = 0;
    for (const auto& q : ring) cx += q[0], cy += q[1];
    const Point c{cx / static_cast<double>(ring.size()), cy / static_cast<double>(ring.size())};
    const auto o = observe(G, c);
    if (o.offset > 120) return;  // a building beyond the network's reach
    build[o.a] += o.wa * area * height, build[o.b] += o.wb * area * height;
  };
  for (const auto& f : P.frontages) add_building(f.footprint, f.height);
  for (const auto& j : inputs["district"]["context"].as_array()) {
    economy::Ring r;
    for (const auto& q : j["footprint"].as_array()) r.push_back({q[0].as_number(), q[1].as_number()});
    add_building(r, j["height_m"].as_number());
  }
  const auto still = [&](std::vector<double> v) {
    double s = 0;
    for (double x : v) s += x;
    std::vector<std::vector<Complex>> out(lines, std::vector<Complex>(G.V, Complex{0, 0}));
    for (std::size_t m = 0; m < lines; m++)
      for (std::size_t i = 0; i < G.V; i++) out[m][i] = s > 0 ? v[i] / s : 0.0;
    return out;
  };
  std::vector<std::vector<Complex>> edge(lines, std::vector<Complex>(G.V, Complex{0, 0}));
  if (edge_from_counters && !outside.empty()) {
    const auto& W = the_week();
    std::vector<std::vector<Complex>> L;
    for (const auto& c : outside) L.push_back(W.forward(c.week));
    double total = 0;
    std::vector<std::size_t> pick;
    for (auto v : boundary) {
      std::size_t best = 0;
      double bd = 1e300;
      for (std::size_t i = 0; i < outside.size(); i++) {
        const double d = std::hypot(outside[i].at[0] - G.pos[v][0], outside[i].at[1] - G.pos[v][1]);
        if (d < bd) bd = d, best = i;
      }
      pick.push_back(best), total += L[best][0].real();
    }
    for (std::size_t m = 0; m < lines; m++)
      for (std::size_t i = 0; i < boundary.size(); i++) edge[m][boundary[i]] = L[pick[i]][m] / total;
  } else {
    double total = 0;
    for (auto v : boundary) total += edge_weight[v];
    for (std::size_t m = 0; m < lines; m++)
      for (auto v : boundary) edge[m][v] = edge_weight[v] / total;
  }
  gr.names = {edge_from_counters ? "the area's edges, fed by the counters outside" : "the area's edges, by the width of their ways", "the buildings' volume", "the main streets' length",
              "the little streets' length", "the lanes' length"};
  gr.pattern = {edge, still(build), still(main), still(little), still(lane)};
  if (trams) gr.names.push_back("the tram streets' length"), gr.pattern.push_back(still(tram));
  return gr;
}

struct Flows {
  Graph G;                          // kept only when asked (keep), as are the next three
  std::vector<double> phi;          // V × 168, the field at the vertices
  std::vector<double> R;            // 168 × F × F
  std::vector<std::vector<Complex>> lines_at_vertices;  // lines × V
  std::vector<Observation> front;   // the frontages on G
  std::vector<double> P0;           // F × 168
  std::size_t vertices = 0, mapped = 0;  // G's vertices; base vertices found on G
  double build_seconds = 0;
};

class Model {
 public:
  const Problem* problem = nullptr;
  Physics phys;
  Graph base;
  Groups groups;
  std::vector<Point> front_foot;                // each frontage's foot on its main segment's centreline
  std::vector<std::vector<Complex>> source;     // lines × base V: the calibrated sources
  std::vector<double> sigma2;                   // lines: the residual's scale
  std::vector<double> p0;                       // F × 168: the plan as built's chances to stop
  std::vector<double> R;                        // 168 × F × F: the response to the attraction, on the network as built
  // Woodbury's links: their ends on the lattice, their conductances, and per line the frontages' lines P̂⁰, the
  // links' G-products TF = O_F G B, H = Bᵀ G B and Bᵀ φ̂.
  std::vector<Observation> link_a, link_b, front_base;
  std::vector<double> link_kappa;
  std::vector<std::vector<Complex>> TF, H, Bphi, P0hat;

  Model(const Problem& P, const Physics& ph) : problem(&P), phys(ph) {
    base = graph(P.base, phys.little, phys.lane);
    for (const auto& f : P.frontages) {
      const auto& l = P.base.line.at(f.segment);
      const auto [along, d] = cg::nearest_on_polyline(l, f.anchor);
      (void)d;
      front_foot.push_back(point_at(l, along));
    }
  }

  // The chances to stop p_f(h) of a plan, as the problem's evaluate() computes them (economy.hpp).
  std::vector<double> chances(const Plan& plan, const economy::Distances& D) const {
    const std::size_t n = problem->size();
    const auto& Q = problem->P;
    std::vector<double> p(n * week, 0.0);
    std::vector<char> open(n);
    for (std::size_t h = 0; h < week; h++) {
      for (std::size_t f = 0; f < n; f++) open[f] = plan.open(f, h);
      for (std::size_t f = 0; f < n; f++) {
        if (!open[f] || !economy::trades(plan.kind[f])) continue;
        const auto k = static_cast<std::size_t>(plan.kind[f]);
        double kappa = 0;
        for (const auto& [g, w] : D.near[f])
          if (open[g]) kappa += w * Q.affinity[static_cast<std::size_t>(plan.kind[g])][k];
        p[f * week + h] = Q.propensity[k] * Q.period[k][economy::band(h)] * std::clamp(1 + Q.catalysis * kappa, Q.catalysis_floor, Q.catalysis_ceiling);
      }
    }
    return p;
  }

  std::uint64_t key(const economy::hv::Bits& made) const {
    std::uint64_t k = 0;
    for (std::size_t i = 0; i < made.size(); i++)
      if (made[i]) k |= std::uint64_t{1} << i;
    return k;
  }

  std::shared_ptr<const Flows> flows(const economy::hv::Bits& made) const {
    const auto k = key(made);
    {
      std::lock_guard<std::mutex> lock(mutex_);
      const auto it = cache_.find(k);
      if (it != cache_.end()) return it->second;
    }
    Flows L;
    bool any = false;
    for (auto b : made) any = any || b;
    if (any) L.P0 = woodbury(made), L.front = front_base, L.vertices = base.V, L.mapped = base.V;
    else L = build(made, false);
    auto F = std::make_shared<Flows>(std::move(L));
    std::lock_guard<std::mutex> lock(mutex_);
    return cache_.emplace(k, std::move(F)).first->second;
  }

  // The response R on the network as built, and Woodbury's products for the links, once the sources are calibrated.
  void prepare() {
    R = build(economy::hv::Bits(problem->links.size(), 0), true).R;
    const Graph& G = base;
    const std::size_t L = problem->links.size(), n = problem->size(), V = G.V;
    link_a.clear(), link_b.clear(), link_kappa.clear(), front_base.clear();
    for (const auto& lk : problem->links) {
      link_a.push_back(observe(G, lk.a)), link_b.push_back(observe(G, lk.b));
      link_kappa.push_back(phys.diffusivity() * way_conductance(lane_way, phys.little, phys.lane) / std::max(1.0, lk.length));
    }
    for (std::size_t f = 0; f < n; f++) front_base.push_back(observe(G, front_foot[f]));
    std::vector<double> UB(V * L, 0.0);  // (Uᵀ B)_kj
    for (std::size_t j = 0; j < L; j++)
      for (const auto& [v, w] : std::vector<std::pair<std::size_t, double>>{{link_a[j].a, link_a[j].wa}, {link_a[j].b, link_a[j].wb}, {link_b[j].a, -link_b[j].wa}, {link_b[j].b, -link_b[j].wb}})
        for (std::size_t k = 0; k < V; k++) UB[k * L + j] += w * G.U(v, k);
    const auto BF = read_modes(G, front_base);
    TF.assign(lines, {}), H.assign(lines, {}), Bphi.assign(lines, {}), P0hat.assign(lines, {});
    std::vector<Complex> us(V);
    for (std::size_t m = 0; m < lines; m++) {
      const auto t = transfer(G, phys, m);
      for (std::size_t k = 0; k < V; k++) {
        Complex s{0, 0};
        for (std::size_t v = 0; v < V; v++) s += G.U(v, k) * source[m][v];
        us[k] = t[k] * s;
      }
      TF[m].assign(n * L, Complex{0, 0}), H[m].assign(L * L, Complex{0, 0}), Bphi[m].assign(L, Complex{0, 0}), P0hat[m].assign(n, Complex{0, 0});
      for (std::size_t k = 0; k < V; k++) {
        for (std::size_t i = 0; i < L; i++) {
          Bphi[m][i] += UB[k * L + i] * us[k];
          for (std::size_t j = 0; j < L; j++) H[m][i * L + j] += UB[k * L + i] * t[k] * UB[k * L + j];
        }
        for (std::size_t f = 0; f < n; f++) {
          P0hat[m][f] += BF[f][k] * us[k];
          for (std::size_t j = 0; j < L; j++) TF[m][f * L + j] += BF[f][k] * t[k] * UB[k * L + j];
        }
      }
    }
  }

  // The frontages' passers P⁰ with a set of links made, by Woodbury's identity on the lattice.
  std::vector<double> woodbury(const economy::hv::Bits& made) const {
    const std::size_t L = problem->links.size(), n = problem->size();
    std::vector<std::size_t> S;
    for (std::size_t j = 0; j < L; j++)
      if (made[j]) S.push_back(j);
    const std::size_t k = S.size();
    std::vector<std::vector<Complex>> lines_of(n, std::vector<Complex>(lines));
    for (std::size_t m = 0; m < lines; m++) {
      std::vector<std::vector<Complex>> A(k, std::vector<Complex>(k));
      std::vector<Complex> b(k);
      for (std::size_t i = 0; i < k; i++) {
        for (std::size_t j = 0; j < k; j++) A[i][j] = H[m][S[i] * L + S[j]] + (i == j ? 1.0 / link_kappa[S[i]] : 0.0);
        b[i] = Bphi[m][S[i]];
      }
      const auto c = complex_solve(A, b);
      for (std::size_t f = 0; f < n; f++) {
        Complex x = P0hat[m][f];
        for (std::size_t i = 0; i < k; i++) x -= TF[m][f * L + S[i]] * c[i];
        lines_of[f][m] = x;
      }
    }
    std::vector<double> P0(n * week);
    for (std::size_t f = 0; f < n; f++) {
      const auto x = the_week().inverse(lines_of[f]);
      for (std::size_t h = 0; h < week; h++) P0[f * week + h] = std::max(0.0, x[h]);
    }
    return P0;
  }

  Flows build(const economy::hv::Bits& made, bool keep) const {
    const auto T0 = std::chrono::steady_clock::now();
    Flows L;
    bool any = false;
    for (auto b : made) any = any || b;
    L.G = any ? graph(problem->with_links(made), phys.little, phys.lane) : base;
    const Graph& G = L.G;
    const std::size_t n = problem->size(), V = G.V;
    for (std::size_t f = 0; f < n; f++) L.front.push_back(observe(G, front_foot[f]));
    // The calibrated sources carried to G: at a vertex of the same position, or else (a subdivision point of a segment
    // the links split) shared between the two ends of G's piece through that point by Whitney's weights; a link's own
    // vertices have none.
    std::vector<std::array<std::pair<std::size_t, double>, 2>> to(base.V);
    for (std::size_t v = 0; v < base.V; v++) {
      const auto u = G.vertex_at(base.pos[v]);
      if (u < V) to[v] = {{{u, 1.0}, {u, 0.0}}}, L.mapped++;
      else {
        const auto o = observe(G, base.pos[v]);
        to[v] = {{{o.a, o.wa}, {o.b, o.wb}}};
      }
    }
    L.lines_at_vertices.assign(lines, std::vector<Complex>(V, Complex{0, 0}));
    for (std::size_t m = 0; m < lines; m++) {
      std::vector<Complex> s(V, Complex{0, 0});
      for (std::size_t v = 0; v < base.V; v++)
        for (const auto& [u, w] : to[v]) s[u] += w * source[m][v];
      L.lines_at_vertices[m] = respond(G, transfer(G, phys, m), s);
    }
    L.phi.assign(V * week, 0.0);
    const auto& W = the_week();
    for (std::size_t v = 0; v < V; v++) {
      std::vector<Complex> X(lines);
      for (std::size_t m = 0; m < lines; m++) X[m] = L.lines_at_vertices[m][v];
      const auto x = W.inverse(X);
      for (std::size_t h = 0; h < week; h++) L.phi[v * week + h] = x[h];
    }
    L.P0.assign(n * week, 0.0);
    for (std::size_t f = 0; f < n; f++)
      for (std::size_t h = 0; h < week; h++) {
        const auto& o = L.front[f];
        L.P0[f * week + h] = std::max(0.0, o.wa * L.phi[o.a * week + h] + o.wb * L.phi[o.b * week + h]);
      }
    L.vertices = V;
    if (!keep) {
      L.G = Graph{}, L.phi.clear(), L.lines_at_vertices.clear();
      L.build_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count();
      return L;
    }
    // The static response to the attraction: M = O_F G₀ (F × V), then R(h) = b M L_φ̄(h) W.
    const double D = phys.diffusivity();
    std::vector<double> M(n * V, 0.0);
    {
      std::vector<double> c(V);
      for (std::size_t f = 0; f < n; f++) {
        const auto& o = L.front[f];
        for (std::size_t k = 0; k < V; k++) c[k] = (o.wa * G.U(o.a, k) + o.wb * G.U(o.b, k)) / (D * G.lambda[k] + phys.leave);
        for (std::size_t v = 0; v < V; v++) {
          double s = 0;
          for (std::size_t k = 0; k < V; k++) s += G.U(v, k) * c[k];
          M[f * V + v] = s;
        }
      }
    }
    L.R.assign(week * n * n, 0.0);
    for (std::size_t h = 0; h < week; h++)
      for (std::size_t g = 0; g < n; g++) {
        // y = L_φ̄ W_g, on the vertices of g's segment and their neighbours.
        std::map<std::size_t, double> y;
        const auto& o = L.front[g];
        const auto bump = [&](std::size_t v, double w) {
          if (w == 0) return;
          for (const auto& [e, u] : G.at[v]) {
            const double ce = D * G.conductance[e] * 0.5 * (L.phi[G.ends[e][0] * week + h] + L.phi[G.ends[e][1] * week + h]);
            const double c = std::max(0.0, ce);
            y[v] += c * w, y[u] -= c * w;
          }
        };
        bump(o.a, o.wa), bump(o.b, o.wb);
        for (std::size_t f = 0; f < n; f++) {
          double s = 0;
          for (const auto& [v, val] : y) s += M[f * V + v] * val;
          L.R[(h * n + f) * n + g] = phys.attraction * s;
        }
      }
    L.build_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count();
    return L;
  }

  // The passers of a plan, F × 168, for the problem's evaluate().
  std::vector<double> passers(const Plan& plan, std::size_t* clamped = nullptr) const {
    const auto L = flows(plan.links);
    const auto D = problem->distances(plan.links);
    const auto p = chances(plan, *D);
    const std::size_t n = problem->size();
    std::vector<double> P(L->P0);
    std::vector<double> dp(n);
    std::size_t c = 0;
    for (std::size_t h = 0; h < week; h++) {
      bool any = false;
      for (std::size_t g = 0; g < n; g++) dp[g] = p[g * week + h] - p0[g * week + h], any = any || dp[g] != 0;
      if (!any) continue;
      for (std::size_t f = 0; f < n; f++) {
        const double* r = &R[(h * n + f) * n];
        double s = 0;
        for (std::size_t g = 0; g < n; g++) s += r[g] * dp[g];
        P[f * week + h] += s;
      }
    }
    for (auto& x : P)
      if (x < 0) x = 0, c++;
    if (clamped) *clamped = c;
    return P;
  }

  economy::Outcome evaluate(const Plan& plan) const { return problem->evaluate(plan, passers(plan)); }

  std::size_t cached() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return cache_.size();
  }

 private:
  mutable std::mutex mutex_;
  mutable std::map<std::uint64_t, std::shared_ptr<const Flows>> cache_;
};

}  // namespace field
