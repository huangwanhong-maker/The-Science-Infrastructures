// The people of the district as a field on its space–time complex (README §2): the base, its Hodge stars and modes,
// the stationary path of the action in the modes, its calibration to the counters by Gaussian conditioning line by
// line, and the response of the field to the plan's attraction.
//
// The base. X is the 1-skeleton of the urban-science domain's district (urban::landscape::district, with the plan's
// links made as the shared problem makes them), crossed with the cycle Z₁₆₈ of the week's hours. On X the Hodge
// stars are ⋆₀(v) = ½ Σ_{e∋v} ℓ_e (the dual length of a vertex, m) and ⋆₁(e) = D c_e / ℓ_e (a conductance; c_e = 1 on
// streets and malls, c_lane on laneways), and Δ₀ = ⋆₀⁻¹ d₀ᵀ ⋆₁ d₀ is physis's Hodge Laplacian of degree 0
// (physis::discretisation::hodge_laplacian), self-adjoint for ⟨a, b⟩₀ = Σ ⋆₀ a b, with eigenpairs Δ₀ u_k = λ_k u_k,
// u_kᵀ ⋆₀ u_l = δ_kl. On Z₁₆₈ the modes are e^{iω_m h}, ω_m = 2πm/168.
//
// The field. φ_v(t), people passing a point an hour on one footpath (a counter's reading; a frontage's passers),
// proportional to the line density of walkers (φ = v₀ ρ / 2 at walking speed v₀, both footpaths alike). Its law is
// the conservation ⋆₀ ∂_t φ = d₀ᵀ J + s − μ ⋆₀ φ with Fick's current J = −⋆₁ d₀ φ (and, from the plan, the drift
// −⋆₁ φ̄ d₀ U): the stationary path of the action S of README §2.2. In the modes it is diagonal:
//
//   φ̂(m) = Σ_k u_k t_k(m) u_kᵀ ŝ(m),   t_k(m) = 1 / (iω_m + D λ_k + μ),
//
// exact for sources band-limited to the week's 84 lines (the continuous-time solution read at the hours).
//
// Calibration (README §3). The sources are not observed. They are a sum of source groups (the area's edges, the
// buildings' volume, the streets' and the lanes' length), each with its own line ŝ_γ(m), and a residual source,
// white in space per unit length with variance σ²(m) ⋆₀ at each line: a Gaussian field whose covariance is the
// action's. The counters observe φ at their points (Whitney's interpolation on their segment), with a nugget
// τ² of the prior variance. Line by line this is universal kriging: generalised least squares for the groups,
// Gaussian conditioning for the residual, and the posterior mean of the sources follows in closed form. The lines
// decouple because the prior is stationary in time and every counter reports every hour of its mean week.
//
// Provisional code of the study's own (README §8): the complex Gaussian conditioning line by line (tyche's spatial
// conditioning is real and one observation at a time), and the mode sums of the weighted propagator (physis's
// field_theory has the unweighted graph Laplacian only).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <functional>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/linalg/cholesky.hpp>
#include <calx/linalg/dense.hpp>
#include <morphe/topology/regular_complex.hpp>
#include <physis/computational_geometry/planar.hpp>
#include <physis/discretisation/hodge_laplacian.hpp>

#include "economy.hpp"

namespace field {

using Complex = std::complex<double>;
using economy::Point;
namespace ul = urban::landscape;
namespace cg = physis::computational_geometry;

inline constexpr std::size_t week = 168;
inline constexpr std::size_t lines = week / 2 + 1;  // m = 0 … 84
inline double omega(std::size_t m) { return 2 * std::numbers::pi * static_cast<double>(m) / week; }

// ---- the parameters (README §3.3: each fitted, sourced or declared an assumption) --------------------------------------

struct Physics {
  double leave = 6;        // μ, 1/h: a walker's stay on the district's streets, ten minutes (assumption)
  double length = 150;     // ξ = √(D/μ), m: the field's screening length (fitted by leave-one-site-out)
  double little = 1;       // c_little: a little street's conductance relative to a main street's (fitted)
  double lane = 1;         // c_lane: a lane's conductance relative to a main street's (fitted)
  double nugget = 0.1;     // τ²: the counters' nugget, relative to the prior variance there (fitted)
  bool relative = false;   // the nugget in proportion to each counter's squared level too, so that the fit weighs errors relatively (fitted)
  double attraction = 5;   // b: an open frontage's potential, −b p_f(h) in units of the walkers' temperature (assumption)
  double speed = 1.34 * 3600;  // v₀, m/h: a walker's mean speed (Weidmann 1993, 1.34 m/s)
  double diffusivity() const { return leave * length * length; }  // D, m²/h
};

// ---- the base: the district's 1-skeleton with its stars and modes ---------------------------------------------------------

// The Hoddle Grid's hierarchy of ways (Hoddle's survey of 1837): its main streets, 99 feet wide; its little streets,
// 33 feet; and the lanes, malls and arcades between them (and every link a plan makes).
enum Way : int { main_street = 0, little_street = 1, lane_way = 2 };
inline int way_class(const std::string& name, const std::string& type) {
  static const std::set<std::string> mains{"Flinders Street", "Collins Street", "Bourke Street", "Lonsdale Street", "La Trobe Street", "Swanston Street",
                                           "Elizabeth Street", "Russell Street", "Queen Street", "William Street", "Exhibition Street", "Spring Street",
                                           "King Street", "Spencer Street", "Bourke Street Mall"};
  static const std::set<std::string> littles{"Flinders Lane", "Little Collins Street", "Little Bourke Street", "Little Lonsdale Street"};
  if (type != "laneway" && mains.count(name)) return main_street;
  if (type != "laneway" && littles.count(name)) return little_street;
  return lane_way;
}
inline double way_conductance(int cls, double little, double lane) { return cls == main_street ? 1.0 : cls == little_street ? little : lane; }

inline bool made_by_link(const std::string& nm) {
  const std::string tail = " carried through";
  return nm == "a passage" || (nm.size() >= tail.size() && nm.compare(nm.size() - tail.size(), tail.size(), tail) == 0);
}

struct Graph {
  std::size_t V = 0, E = 0;
  std::vector<Point> pos;
  std::vector<std::array<std::size_t, 2>> ends;  // tail, head
  std::vector<double> length;
  std::vector<std::string> type, name;
  std::vector<char> link;  // made by a link of the plan
  std::vector<int> cls;    // the way's class (Way)
  std::vector<std::vector<Point>> line;
  std::vector<double> star0;        // ⋆₀, m
  std::vector<double> conductance;  // c_e / ℓ_e, 1/m: ⋆₁ per unit diffusivity
  std::vector<std::vector<std::pair<std::size_t, std::size_t>>> at;  // vertex ↦ (edge, the other end)
  std::vector<double> lambda;  // Δ₀'s eigenvalues per unit diffusivity, 1/m², ascending
  calx::linalg::Matrix U;      // U(v, k) = u_k(v), ⋆₀-orthonormal
  morphe::algebra::ChainComplex chain;
  std::map<std::pair<long long, long long>, std::size_t> by_position;  // vertices by their position to the millimetre

  static std::pair<long long, long long> key(Point p) { return {std::llround(p[0] * 1000), std::llround(p[1] * 1000)}; }
  std::size_t vertex_at(Point p) const {
    const auto it = by_position.find(key(p));
    return it == by_position.end() ? V : it->second;
  }
};

inline Point point_at(const std::vector<Point>& l, double along) {
  double run = 0;
  for (std::size_t i = 0; i + 1 < l.size(); i++) {
    const double d = std::hypot(l[i + 1][0] - l[i][0], l[i + 1][1] - l[i][1]);
    if (run + d >= along || i + 2 == l.size()) {
      const double t = d > 0 ? std::clamp((along - run) / d, 0.0, 1.0) : 0.0;
      return {l[i][0] + t * (l[i + 1][0] - l[i][0]), l[i][1] + t * (l[i + 1][1] - l[i][1])};
    }
    run += d;
  }
  return l.front();
}

// The lattice: the domain's district with every segment cut into equal pieces no longer than h (a subdivision of the
// complex, which changes no homology), so that the lattice resolves the field's screening length ξ = 45 m; h ≤ 0
// keeps the domain's cells. Each block's boundary runs through its segments' pieces with the segments' signs, so the
// chain complex (vertices, pieces, blocks) is the subdivided one and d₁d₀ = 0 still holds.
inline constexpr double piece_length = 15;  // m (README §2.1)

inline Graph graph(const ul::District& D, double little, double lane, double h = piece_length) {
  Graph G;
  const auto local = D.complex.local_index();
  std::vector<std::size_t> vertex_of(D.complex.size(), 0);
  for (std::size_t i = 0; i < D.complex.size(); i++) {
    const auto& c = D.complex.cells[i];
    if (c.dim != 0) continue;
    if (local[i] != G.pos.size()) throw std::logic_error("field::graph: the 0-cells' local order");
    vertex_of[i] = G.pos.size();
    G.pos.push_back(D.point.at(c.id));
  }
  // The segments, each cut into k pieces along its polyline.
  std::vector<std::vector<std::size_t>> pieces_of;  // per segment of the domain (local order), its pieces in order tail → head
  for (std::size_t i = 0; i < D.complex.size(); i++) {
    const auto& c = D.complex.cells[i];
    if (c.dim != 1) continue;
    if (local[i] != pieces_of.size()) throw std::logic_error("field::graph: the 1-cells' local order");
    const auto e = D.complex.ends(i);
    const auto& poly = D.line.at(c.id);
    const double L = std::max(0.5, cg::polyline_length(poly));
    const std::size_t k = h > 0 ? std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(L / h - 1e-9))) : 1;
    const std::string nm = D.name.count(c.id) ? D.name.at(c.id) : "";
    std::vector<std::size_t> chain_of_vertices{vertex_of[e[0]]};
    std::vector<std::vector<Point>> lines(k);
    {
      // Split the polyline at arc lengths L j / k.
      std::vector<double> cum{0};
      for (std::size_t q = 0; q + 1 < poly.size(); q++) cum.push_back(cum.back() + std::hypot(poly[q + 1][0] - poly[q][0], poly[q + 1][1] - poly[q][1]));
      const double total = cum.back();
      for (std::size_t j = 0; j < k; j++) {
        const double s0 = total * static_cast<double>(j) / static_cast<double>(k), s1 = total * static_cast<double>(j + 1) / static_cast<double>(k);
        auto& out = lines[j];
        out.push_back(j == 0 ? poly.front() : point_at(poly, s0));
        for (std::size_t q = 1; q + 1 < poly.size(); q++)
          if (cum[q] > s0 && cum[q] < s1) out.push_back(poly[q]);
        out.push_back(j + 1 == k ? poly.back() : point_at(poly, s1));
      }
    }
    for (std::size_t j = 1; j < k; j++) {
      chain_of_vertices.push_back(G.pos.size());
      G.pos.push_back(lines[j].front());
    }
    chain_of_vertices.push_back(vertex_of[e[1]]);
    std::vector<std::size_t> mine;
    for (std::size_t j = 0; j < k; j++) {
      mine.push_back(G.ends.size());
      G.ends.push_back({chain_of_vertices[j], chain_of_vertices[j + 1]});
      G.length.push_back(std::max(0.5, L / static_cast<double>(k)));
      G.type.push_back(c.type);
      G.name.push_back(nm);
      G.link.push_back(made_by_link(nm));
      G.cls.push_back(way_class(nm, c.type));
      G.line.push_back(lines[j]);
    }
    pieces_of.push_back(mine);
  }
  G.V = G.pos.size();
  G.E = G.ends.size();
  G.star0.assign(G.V, 0.0);
  G.at.assign(G.V, {});
  for (std::size_t e = 0; e < G.E; e++) {
    G.conductance.push_back(way_conductance(G.cls[e], little, lane) / G.length[e]);
    for (int q = 0; q < 2; q++) G.star0[G.ends[e][q]] += 0.5 * G.length[e];
    G.at[G.ends[e][0]].push_back({e, G.ends[e][1]});
    G.at[G.ends[e][1]].push_back({e, G.ends[e][0]});
  }
  for (std::size_t v = 0; v < G.V; v++) {
    if (!(G.star0[v] > 0)) G.star0[v] = 0.5;  // an isolated vertex keeps a small star
    G.by_position[Graph::key(G.pos[v])] = v;
  }
  // The subdivided chain complex.
  const auto C = morphe::topology::chain_complex(D.complex);
  G.chain.n = {G.V, G.E, C.n.size() > 2 ? C.n[2] : 0};
  G.chain.boundary.resize(3);
  G.chain.boundary[1] = morphe::algebra::SparseMatrix(G.V, G.E);
  for (std::size_t e = 0; e < G.E; e++) {
    G.chain.boundary[1].add(G.ends[e][0], e, morphe::algebra::Integer(-1));
    G.chain.boundary[1].add(G.ends[e][1], e, morphe::algebra::Integer(1));
  }
  G.chain.boundary[2] = morphe::algebra::SparseMatrix(G.E, G.chain.n[2]);
  if (C.n.size() > 2)
    for (std::size_t f = 0; f < C.n[2]; f++)
      for (const auto& en : C.boundary[2].columns[f])
        for (auto piece : pieces_of[en.index]) G.chain.boundary[2].add(piece, f, en.value);
  G.chain.validate();
  const auto H = physis::discretisation::hodge_laplacian(G.chain, 0, {G.star0, G.conductance});
  G.lambda = H.values;
  for (auto& l : G.lambda) l = std::max(0.0, l);
  G.U = H.vectors;
  return G;
}

// ---- observing the field at a point -------------------------------------------------------------------------------------

struct Observation {
  std::size_t edge = 0, a = 0, b = 0;  // the segment and its tail and head
  double wa = 1, wb = 0;               // Whitney's weights: φ(point) = wa φ_a + wb φ_b
  double offset = 0;                   // the point's distance from the segment's centreline, m
  Point on{};                          // the point's foot on the centreline
};

// The nearest segment to p among those `allowed` accepts, and Whitney's weights of p's foot on it.
inline Observation observe(const Graph& G, Point p, const std::function<bool(std::size_t)>& allowed) {
  Observation o;
  double best = std::numeric_limits<double>::infinity(), along = 0;
  for (std::size_t e = 0; e < G.E; e++) {
    if (!allowed(e)) continue;
    const auto [s, d] = cg::nearest_on_polyline(G.line[e], p);
    if (d < best) best = d, o.edge = e, along = s;
  }
  if (!std::isfinite(best)) throw std::invalid_argument("field::observe: no segment allowed");
  const double L = std::max(1e-9, cg::polyline_length(G.line[o.edge]));
  const double s = std::clamp(along / L, 0.0, 1.0);
  o.a = G.ends[o.edge][0], o.b = G.ends[o.edge][1], o.wa = 1 - s, o.wb = s, o.offset = best;
  o.on = point_at(G.line[o.edge], along);
  return o;
}
inline Observation observe(const Graph& G, Point p) {
  return observe(G, p, [&](std::size_t e) { return !G.link[e]; });
}

// ---- the modes on the week ------------------------------------------------------------------------------------------------

// x̂(m) = (1/168) Σ_h x_h e^{−iω_m h}, m = 0 … 84; and back: x_h = Re Σ_m c_m x̂(m) e^{iω_m h}, c_0 = c_84 = 1, else 2.
struct Week {
  std::vector<Complex> e;  // e[m * 168 + h] = e^{iω_m h}
  Week() : e(lines * week) {
    for (std::size_t m = 0; m < lines; m++)
      for (std::size_t h = 0; h < week; h++) e[m * week + h] = std::polar(1.0, 2 * std::numbers::pi * static_cast<double>((m * h) % week) / week);
  }
  std::vector<Complex> forward(std::span<const double> x) const {
    std::vector<Complex> X(lines);
    for (std::size_t m = 0; m < lines; m++) {
      Complex s{0, 0};
      for (std::size_t h = 0; h < week; h++) s += x[h] * std::conj(e[m * week + h]);
      X[m] = s / static_cast<double>(week);
    }
    return X;
  }
  std::vector<double> inverse(std::span<const Complex> X) const {
    std::vector<double> x(week, 0.0);
    for (std::size_t m = 0; m < lines; m++) {
      const double c = (m == 0 || m == lines - 1) ? 1.0 : 2.0;
      for (std::size_t h = 0; h < week; h++) x[h] += c * (X[m] * e[m * week + h]).real();
    }
    return x;
  }
};

inline const Week& the_week() {
  static const Week w;
  return w;
}

// t_k(m) = 1 / (iω_m + D λ_k + μ).
inline std::vector<Complex> transfer(const Graph& G, const Physics& phys, std::size_t m) {
  std::vector<Complex> t(G.V);
  const double D = phys.diffusivity();
  for (std::size_t k = 0; k < G.V; k++) t[k] = 1.0 / Complex(D * G.lambda[k] + phys.leave, omega(m));
  return t;
}

// The field at the vertices from a source at line m: φ̂ = U diag(t) Uᵀ ŝ.
inline std::vector<Complex> respond(const Graph& G, const std::vector<Complex>& t, std::span<const Complex> s) {
  std::vector<Complex> c(G.V, Complex{0, 0}), phi(G.V, Complex{0, 0});
  for (std::size_t v = 0; v < G.V; v++)
    if (s[v] != Complex{0, 0})
      for (std::size_t k = 0; k < G.V; k++) c[k] += G.U(v, k) * s[v];
  for (std::size_t k = 0; k < G.V; k++) c[k] *= t[k];
  for (std::size_t v = 0; v < G.V; v++) {
    Complex x{0, 0};
    for (std::size_t k = 0; k < G.V; k++) x += G.U(v, k) * c[k];
    phi[v] = x;
  }
  return phi;
}

// B = O U: the modes read at the observations (n × V).
inline std::vector<std::vector<double>> read_modes(const Graph& G, const std::vector<Observation>& O) {
  std::vector<std::vector<double>> B(O.size(), std::vector<double>(G.V, 0.0));
  for (std::size_t i = 0; i < O.size(); i++)
    for (std::size_t k = 0; k < G.V; k++) B[i][k] = O[i].wa * G.U(O[i].a, k) + O[i].wb * G.U(O[i].b, k);
  return B;
}

// ---- small complex linear algebra (the study's own: README §8) --------------------------------------------------------

// A real symmetric positive-definite system with complex right-hand sides, by calx's Cholesky on the real and
// imaginary parts.
struct RealSpd {
  calx::linalg::Cholesky f;
  explicit RealSpd(const calx::linalg::Matrix& A) : f(calx::linalg::cholesky(A)) {}
  std::vector<Complex> solve(const std::vector<Complex>& b) const {
    std::vector<double> re(b.size()), im(b.size());
    for (std::size_t i = 0; i < b.size(); i++) re[i] = b[i].real(), im[i] = b[i].imag();
    const auto x = calx::linalg::cholesky_solve(f, re), y = calx::linalg::cholesky_solve(f, im);
    std::vector<Complex> out(b.size());
    for (std::size_t i = 0; i < b.size(); i++) out[i] = {x[i], y[i]};
    return out;
  }
};

// A small complex system A x = b by Gaussian elimination with partial pivoting.
inline std::vector<Complex> complex_solve(std::vector<std::vector<Complex>> A, std::vector<Complex> b) {
  const std::size_t n = b.size();
  for (std::size_t c = 0; c < n; c++) {
    std::size_t p = c;
    for (std::size_t r = c + 1; r < n; r++)
      if (std::abs(A[r][c]) > std::abs(A[p][c])) p = r;
    if (std::abs(A[p][c]) == 0) throw std::runtime_error("field::complex_solve: singular");
    std::swap(A[p], A[c]), std::swap(b[p], b[c]);
    for (std::size_t r = c + 1; r < n; r++) {
      const Complex f = A[r][c] / A[c][c];
      for (std::size_t j = c; j < n; j++) A[r][j] -= f * A[c][j];
      b[r] -= f * b[c];
    }
  }
  std::vector<Complex> x(n);
  for (std::size_t r = n; r-- > 0;) {
    Complex s = b[r];
    for (std::size_t j = r + 1; j < n; j++) s -= A[r][j] * x[j];
    x[r] = s / A[r][r];
  }
  return x;
}

// A complex Hermitian positive-definite system H x = b, as the real symmetric [[Re, −Im], [Im, Re]].
inline std::vector<Complex> hermitian_solve(const std::vector<std::vector<Complex>>& H, const std::vector<Complex>& b) {
  const std::size_t n = b.size();
  calx::linalg::Matrix A(2 * n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = 0; j < n; j++) {
      A(i, j) = A(n + i, n + j) = H[i][j].real();
      A(i, n + j) = -H[i][j].imag();
      A(n + i, j) = H[i][j].imag();
    }
  for (std::size_t i = 0; i < 2 * n; i++)
    for (std::size_t j = i + 1; j < 2 * n; j++) A(i, j) = A(j, i) = 0.5 * (A(i, j) + A(j, i));
  std::vector<double> r(2 * n);
  for (std::size_t i = 0; i < n; i++) r[i] = b[i].real(), r[n + i] = b[i].imag();
  const auto x = calx::linalg::cholesky_solve(calx::linalg::cholesky(A), r);
  std::vector<Complex> out(n);
  for (std::size_t i = 0; i < n; i++) out[i] = {x[i], x[n + i]};
  return out;
}

// ---- calibration: universal kriging line by line ---------------------------------------------------------------------------

struct Groups {
  std::vector<std::string> names;
  std::vector<std::vector<std::vector<Complex>>> pattern;  // [group][line][vertex]: a pattern on the base graph's vertices
};

// One line's fit on a set of counters: the groups' lines, the residual's weights w = Σ₁⁻¹ r, the scale σ²(m).
struct LineFit {
  std::vector<Complex> beta;   // G
  std::vector<Complex> w;      // n (in the counters' order of the set)
  std::vector<std::vector<Complex>> normal;  // X^H Σ₁⁻¹ X, for the posterior's group term
  double sigma2 = 0;
  std::vector<std::size_t> set;
};

struct Kriging {
  const Graph* G = nullptr;
  Physics phys;
  std::vector<Observation> obs;                  // the counters
  std::vector<std::vector<double>> B;            // n × V
  std::vector<std::vector<std::vector<Complex>>> UA;  // lines × G × V: Uᵀ a_γ(m)
  std::vector<std::vector<Complex>> y;           // n × lines
  std::vector<std::vector<Complex>> t;           // lines × V
  std::vector<std::vector<std::vector<Complex>>> X;  // lines × n × G
  std::vector<std::vector<std::vector<double>>> K;   // lines × n × n
  std::vector<double> level;                          // per counter, its squared weekly mean over the counters' median (for a relative nugget)
  std::size_t groups = 0;

  Kriging(const Graph& graph, const Physics& p, std::vector<Observation> o, const std::vector<std::vector<Complex>>& lines_of, const Groups& gr)
      : G(&graph), phys(p), obs(std::move(o)), y(lines_of) {
    B = read_modes(*G, obs);
    groups = gr.pattern.size();
    UA.assign(lines, std::vector<std::vector<Complex>>(groups, std::vector<Complex>(G->V, Complex{0, 0})));
    for (std::size_t m = 0; m < lines; m++)
      for (std::size_t g = 0; g < groups; g++)
        for (std::size_t v = 0; v < G->V; v++)
          if (gr.pattern[g][m][v] != Complex{0, 0})
            for (std::size_t k = 0; k < G->V; k++) UA[m][g][k] += G->U(v, k) * gr.pattern[g][m][v];
    const std::size_t n = obs.size();
    {
      std::vector<double> means;
      for (std::size_t i = 0; i < n; i++) means.push_back(std::abs(y[i][0]));
      auto sorted = means;
      std::sort(sorted.begin(), sorted.end());
      const double med = sorted[n / 2];
      for (double x : means) level.push_back(x * x / (med * med));
    }
    t.resize(lines), X.resize(lines), K.resize(lines);
    for (std::size_t m = 0; m < lines; m++) {
      t[m] = transfer(*G, phys, m);
      X[m].assign(n, std::vector<Complex>(groups, Complex{0, 0}));
      K[m].assign(n, std::vector<double>(n, 0.0));
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t g = 0; g < groups; g++) {
          Complex s{0, 0};
          for (std::size_t k = 0; k < G->V; k++) s += B[i][k] * t[m][k] * UA[m][g][k];
          X[m][i][g] = s;
        }
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = i; j < n; j++) {
          double s = 0;
          for (std::size_t k = 0; k < G->V; k++) s += B[i][k] * std::norm(t[m][k]) * B[j][k];
          K[m][i][j] = K[m][j][i] = s;
        }
    }
  }

  // Σ₁ = K + τ² diag(K) on a set.
  calx::linalg::Matrix sigma1(std::size_t m, const std::vector<std::size_t>& set) const {
    calx::linalg::Matrix S(set.size());
    for (std::size_t a = 0; a < set.size(); a++)
      for (std::size_t b = 0; b < set.size(); b++)
        S(a, b) = K[m][set[a]][set[b]] + (a == b ? phys.nugget * K[m][set[a]][set[a]] * (phys.relative && !level.empty() ? level[set[a]] : 1.0) : 0.0);
    return S;
  }

  LineFit fit(std::size_t m, const std::vector<std::size_t>& set) const {
    LineFit F;
    F.set = set;
    const std::size_t n = set.size();
    const RealSpd S(sigma1(m, set));
    std::vector<std::vector<Complex>> SX(groups);  // Σ₁⁻¹ X, column by column
    for (std::size_t g = 0; g < groups; g++) {
      std::vector<Complex> col(n);
      for (std::size_t a = 0; a < n; a++) col[a] = X[m][set[a]][g];
      SX[g] = S.solve(col);
    }
    std::vector<Complex> ys(n);
    for (std::size_t a = 0; a < n; a++) ys[a] = y[set[a]][m];
    const auto Sy = S.solve(ys);
    F.normal.assign(groups, std::vector<Complex>(groups, Complex{0, 0}));
    std::vector<Complex> rhs(groups, Complex{0, 0});
    for (std::size_t g = 0; g < groups; g++) {
      for (std::size_t h = 0; h < groups; h++)
        for (std::size_t a = 0; a < n; a++) F.normal[g][h] += std::conj(X[m][set[a]][g]) * SX[h][a];
      for (std::size_t a = 0; a < n; a++) rhs[g] += std::conj(X[m][set[a]][g]) * Sy[a];
    }
    // A ridge of 10⁻¹⁰ of the trace keeps the normal equations definite when a group's line is not seen.
    double tr = 0;
    for (std::size_t g = 0; g < groups; g++) tr += F.normal[g][g].real();
    for (std::size_t g = 0; g < groups; g++) F.normal[g][g] += 1e-10 * tr / static_cast<double>(groups);
    F.beta = hermitian_solve(F.normal, rhs);
    std::vector<Complex> r(n);
    for (std::size_t a = 0; a < n; a++) {
      Complex s = ys[a];
      for (std::size_t g = 0; g < groups; g++) s -= X[m][set[a]][g] * F.beta[g];
      r[a] = s;
    }
    F.w = S.solve(r);
    double q = 0;
    for (std::size_t a = 0; a < n; a++) q += (std::conj(r[a]) * F.w[a]).real();
    F.sigma2 = n > groups ? q / static_cast<double>(n - groups) : 0.0;
    return F;
  }

  // The prediction at counter q (in or out of the fit's set): X_q β + K_{q,set} w.
  Complex predict(std::size_t m, const LineFit& F, std::size_t q) const {
    Complex s{0, 0};
    for (std::size_t g = 0; g < groups; g++) s += X[m][q][g] * F.beta[g];
    for (std::size_t a = 0; a < F.set.size(); a++) s += K[m][q][F.set[a]] * F.w[a];
    return s;
  }

  // The posterior mean of the sources on the graph's vertices: Σ_γ a_γ β_γ + ⋆₀ ⊙ U (t̄ ⊙ Bᵀ w).
  std::vector<Complex> source(std::size_t m, const LineFit& F, const Groups& gr) const {
    std::vector<Complex> s(G->V, Complex{0, 0});
    for (std::size_t g = 0; g < groups; g++)
      for (std::size_t v = 0; v < G->V; v++) s[v] += gr.pattern[g][m][v] * F.beta[g];
    std::vector<Complex> c(G->V, Complex{0, 0});
    for (std::size_t a = 0; a < F.set.size(); a++)
      for (std::size_t k = 0; k < G->V; k++) c[k] += B[F.set[a]][k] * F.w[a];
    for (std::size_t k = 0; k < G->V; k++) c[k] *= std::conj(t[m][k]);
    for (std::size_t v = 0; v < G->V; v++) {
      Complex x{0, 0};
      for (std::size_t k = 0; k < G->V; k++) x += G->U(v, k) * c[k];
      s[v] += G->star0[v] * x;
    }
    return s;
  }
};

}  // namespace field
