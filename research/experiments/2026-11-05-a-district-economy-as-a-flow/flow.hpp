// People as a fluid on the district's streets: study (2)'s own provisional model of flow (README §2), on the shared
// problem "A district's economy" (../2026-11-05-a-district-economy-the-problem).
//
// THE NETWORK. The district's walking network as the urban-science domain builds it (economy::Problem::with_links), with
// every candidate link made: the superset network. A link not made is that network with the link's segment removed;
// the junctions the link would make on the streets it meets stay as points along them. Each segment has a walkable
// width (Widths: assumptions by class of street, after Hoddle's grid of 1837, whose streets are 99 ft wide and whose
// little streets 33 ft) and a length.
//
// THE FLUID. People are many streams, each a trip from an origin junction to a destination junction (a class of walkers
// in the sense of Hughes 2002, "A continuum theory for the flow of pedestrians", Transp. Res. B 36, 507–535: walkers who
// share a destination share a potential, here the log of the sum of the weights of the ways left to it). A stream is
// conserved: at every junction but its two ends what flows in flows out. Its walkers choose their way by a logit over
// every walk (Bell 1995, "Alternatives to Dial's logit assignment algorithm", Transp. Res. B 29, 287–295; Akamatsu 1996,
// "Cyclic flows, Markov process and stochastic traffic assignment", Transp. Res. B 30, 369–386), the walk's states being
// the directed segments, so that a walker does not turn back except at a dead end (the link-based form of Fosgerau,
// Frejinger and Karlström 2013, "A link based network route choice model with unrestricted choice set", Transp. Res. B
// 56, 70–80). Entering a directed segment e has the weight
//
//   ω_e = exp(−θ ℓ_e) (w_e / w_ref)^(ℓ_e / ℓ_parent),
//
// so that a walk's log-weight is −θ (its length) plus the log of the widths it takes: between two parallel segments of
// the network as built the walkers divide in proportion to width, as a fluid's flux divides by conductance (ℓ_parent is
// the length of the segment as built the piece belongs to, so that a segment cut by a link's junction keeps its weight).
// With M the transitions between directed segments (M_ab = ω_b when b leaves the junction a enters, b not a reversed
// unless that junction is a dead end) and Z = (I − M)⁻¹ = Σ_k M^k, the sum over every walk, convergent when M's
// spectral radius is below one (checked), a trip from junction o that ends on first entering junction d visits the
// directed segment c, on average,
//
//   f_c(o → d) = π_od(c) σ_d(c) / W_od,
//   π_od = s_oᵀ Z⁽ᵈ⁾,  σ_d = Z⁽ᵈ⁾ 1_{A_d},  W_od = s_oᵀ Z⁽ᵈ⁾ 1_{A_d},
//
// s_o the weights of the segments leaving o, A_d the segments entering d, and Z⁽ᵈ⁾ the walks with A_d absorbing: by
// Woodbury, Z⁽ᵈ⁾ = Z − Z[:, A_d] Z[A_d, A_d]⁻¹ (Z − I)[A_d, :], so that σ_d = Z[:, A_d] Z[A_d, A_d]⁻¹ 1 (the walks through
// c, split there, over all walks). The flow is linear in the trips: the streams superpose.
//
// THE SOURCES AND SINKS (each hour h of the mean week), on the network's largest component:
//   entrances: the network's dead ends on the edge of its box, where the streets leave it; each with a gross flow
//     g_i(h) (people an hour, both ways) from the nearest of the City's counters outside the box on the same street, else
//     (for a street or mall) from one near it scaled by width, else the median of the counters on streets scaled by width
//     (assumptions); counters at crossings, which count people crossing a street, are left out;
//   buildings: the parcels and the context standing in the box, their floor (footprint × storeys) as mass on the two
//     ends of their nearest segment, offices and shops (work) apart from dwellings (home);
//   attractors: CLUE's establishments of the trading kinds in the box outside the 61 parcels, open by their kind's
//     typical hours, each of weight w_k (assumption), and the plan's frontages, open by the plan's hours.
// The trips of an hour are a sum of patterns, each a fixed spatial law of origins and destinations times an amplitude
// x(h) ≥ 0 fitted to the counters inside the box (non-negative least squares, tyche):
//   arrive (work, home): entrance i → buildings, g_i(h) × the buildings' shares of floor;
//   leave (work, home): buildings → entrance i, likewise;
//   through: entrance i → entrance j, g_i(h) g_j(h) / Σ g;
//   local: building a → building b in proportion to m_a m_b exp(−d_ab / λ);
//   visit: from the entrances and the buildings (half and half) to each open attractor and back.
// The plan enters through the visits: an open frontage of kind k draws w_k x_V(h) / W visits an hour, W the weight of
// every attractor as built. So the passers at a frontage are B_f(h; links) + c(h) Σ_g w_{k_g} [g open at h] R_fg(links):
// affine in the plan's openings, with the links entering through Z, which a made link changes by a low-rank update
// (Woodbury: the transitions into its two directed segments, and the turn back at a dead end it carries on), so that
// each of the 2,048 sets of links is a small solve, not a new inverse.
//
// WHAT A FRONTAGE AND A COUNTER SEE. The share σ of a segment's walkers that pass one footpath (and so one frontage or
// one counter): ½ on streets and malls, which have two sides, 1 on laneways (an assumption).
//
// THE CIRCULATION. The net flow F = f(u → v) − f(v → u) of all the streams on the network's complex (the junctions, the
// segments, and the blocks physis traces as 2-cells) is decomposed by the economics domain's circulation (physis's flow
// with conductances c_e = w_e / ℓ_e): a driven flow up an attraction α, which is the flow a Darcy fluid with the same
// sources would carry (Thomson's principle: of all flows with given divergence the driven part dissipates least), the
// circulation round the blocks (eddies), and the circulation round the holes. In a Darcy fluid the last two vanish; in
// the walkers they measure how their choice of way departs from it.
//
// Everything here is the study's own and provisional: a model the domains lack (urban-science's `movement` is planned).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numbers>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/linalg/rectangular.hpp>
#include <economics/circulation/flows.hpp>
#include <physis/computational_geometry/planar.hpp>
#include <tyche/regression/nnls.hpp>
#include <urban/flow/currents.hpp>
#include <urban/rhythm/record.hpp>

#include "economy.hpp"

namespace flow {

using economy::Json;
using economy::Point;
namespace ul = urban::landscape;
namespace cg = physis::computational_geometry;
namespace la = calx::linalg;
namespace hv = heuresis::variation;
inline constexpr std::size_t H = economy::week;
inline constexpr std::size_t none = static_cast<std::size_t>(-1);

// A file whose path the inputs name (the data library's cache), read whole.
inline std::string read_text(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("flow: cannot read " + path);
  std::stringstream b;
  b << in.rdbuf();
  return b.str();
}

// ---- walkable widths (assumptions) ------------------------------------------------------------------------------------

struct Widths {
  double swanston = 12;     // m: Swanston Street, closed to private cars since 1992, its footpaths widened
  double main = 9;          // the other streets of Hoddle's 99 ft (30.2 m) reserve: two footpaths of 4.5 m
  double little = 4;        // the little streets of 33 ft (10.1 m): two footpaths of 2 m
  double street = 4;        // other streets (Royal Lane, Russell Place, unnamed)
  double bourke_mall = 12;  // the Bourke Street Mall beside its tramway
  double mall = 5;          // other pedestrian ways and arcades
  double laneway = 3.5;     // the city model's default width of a service way
  double lane_link = 3.5;   // a lane carried through
  double passage = 5;       // a passage through a block (an arcade)
  double reference = 12;    // w_ref
};

inline bool main_street(const std::string& n) {
  static const std::set<std::string> s{"Collins Street", "Bourke Street", "Elizabeth Street", "Russell Street", "Lonsdale Street", "Flinders Street",
                                       "La Trobe Street", "Queen Street", "William Street", "Exhibition Street", "Spring Street", "King Street"};
  return s.count(n) > 0;
}
inline bool little_street(const std::string& n) {
  static const std::set<std::string> s{"Flinders Lane", "Little Collins Street", "Little Bourke Street", "Little Lonsdale Street"};
  return s.count(n) > 0;
}
inline double width_of(const Widths& W, const std::string& kind, const std::string& name) {
  if (kind == "laneway") return W.laneway;
  if (kind == "mall") return name.rfind("Bourke Street", 0) == 0 ? W.bourke_mall : W.mall;
  if (name == "Swanston Street") return W.swanston;
  if (main_street(name)) return W.main;
  if (little_street(name)) return W.little;
  return W.street;
}
// The share of a segment's walkers that pass one footpath, one frontage or one counter (assumption).
inline double side_share(const std::string& kind) { return kind == "laneway" ? 1.0 : 0.5; }

// ---- the superset network ------------------------------------------------------------------------------------------------

struct Edge {
  std::size_t u = 0, v = 0;  // tail, head: the polyline runs from u to v
  std::vector<Point> line;
  double length = 0, width = 0, parent = 0;  // m; parent: the length of the segment as built it is a piece of
  std::string kind, name;
  int link = -1;             // the candidate link it makes, or −1
};

struct Network {
  std::vector<Point> at;
  std::vector<Edge> edges;
  std::vector<std::vector<std::size_t>> link_edges;  // per candidate link, its edges
  std::size_t nodes() const { return at.size(); }
  bool present(std::size_t e, const hv::Bits& made) const { return edges[e].link < 0 || made[static_cast<std::size_t>(edges[e].link)]; }
  // Directed segments: state 2e runs u → v, 2e + 1 runs v → u.
  std::size_t states() const { return 2 * edges.size(); }
  std::size_t tail(std::size_t s) const { return s % 2 ? edges[s / 2].v : edges[s / 2].u; }
  std::size_t head(std::size_t s) const { return s % 2 ? edges[s / 2].u : edges[s / 2].v; }
};

inline double point_segment_distance(const Point& p, const Point& a, const Point& b) {
  const double ex = b[0] - a[0], ey = b[1] - a[1], L2 = ex * ex + ey * ey;
  const double t = L2 > 0 ? std::clamp(((p[0] - a[0]) * ex + (p[1] - a[1]) * ey) / L2, 0.0, 1.0) : 0.0;
  return std::hypot(p[0] - a[0] - t * ex, p[1] - a[1] - t * ey);
}
inline Point point_along(const std::vector<Point>& l, double s) {
  for (std::size_t k = 0; k + 1 < l.size(); k++) {
    const double d = std::hypot(l[k + 1][0] - l[k][0], l[k + 1][1] - l[k][1]);
    if (s <= d || k + 2 == l.size()) {
      const double t = d > 0 ? std::clamp(s / d, 0.0, 1.0) : 0.0;
      return {l[k][0] + t * (l[k + 1][0] - l[k][0]), l[k][1] + t * (l[k + 1][1] - l[k][1])};
    }
    s -= d;
  }
  return l.back();
}

inline Network superset(const economy::Problem& P, const Widths& W) {
  const ul::District D = P.with_links(hv::Bits(P.links.size(), 1));
  Network N;
  std::map<std::size_t, std::size_t> node;
  for (std::size_t i = 0; i < D.complex.size(); i++)
    if (D.complex.cells[i].dim == 0) node[i] = N.at.size(), N.at.push_back(D.point.at(D.complex.cells[i].id));
  N.link_edges.assign(P.links.size(), {});
  for (std::size_t i = 0; i < D.complex.size(); i++) {
    const auto& c = D.complex.cells[i];
    if (c.dim != 1) continue;
    const auto ends = D.complex.ends(i);
    Edge e;
    e.u = node.at(ends[0]), e.v = node.at(ends[1]);
    e.line = D.line.at(c.id);
    const auto d = [](const Point& a, const Point& b) { return std::hypot(a[0] - b[0], a[1] - b[1]); };
    if (d(e.line.front(), N.at[e.u]) > d(e.line.front(), N.at[e.v])) std::reverse(e.line.begin(), e.line.end());
    e.length = cg::polyline_length(e.line);
    e.kind = c.type;
    e.name = D.name.count(c.id) ? D.name.at(c.id) : "";
    const std::string tail = " carried through";
    const bool is_link = e.name == "a passage" || (e.name.size() >= tail.size() && e.name.compare(e.name.size() - tail.size(), tail.size(), tail) == 0);
    if (is_link) {
      double best = std::numeric_limits<double>::infinity();
      for (std::size_t k = 0; k < P.links.size(); k++) {
        const double dk = std::max(point_segment_distance(e.line.front(), P.links[k].a, P.links[k].b), point_segment_distance(e.line.back(), P.links[k].a, P.links[k].b));
        if (dk < best) best = dk, e.link = static_cast<int>(k);
      }
      if (!(best < 1.5)) throw std::logic_error("flow::superset: a link's segment on no candidate link");
      e.width = P.links[static_cast<std::size_t>(e.link)].kind == "passage" ? W.passage : W.lane_link;
    } else {
      e.width = width_of(W, e.kind, e.name);
    }
    N.edges.push_back(std::move(e));
  }
  for (std::size_t k = 0; k < N.edges.size(); k++) {
    auto& e = N.edges[k];
    if (e.link >= 0) {
      N.link_edges[static_cast<std::size_t>(e.link)].push_back(k);
      continue;
    }
    const Point mid = point_along(e.line, e.length / 2);
    double best = std::numeric_limits<double>::infinity(), parent = e.length;
    for (const auto& [id, l] : P.base.line) {
      const auto [s, dist] = cg::nearest_on_polyline(l, mid);
      (void)s;
      if (dist < best) best = dist, parent = cg::polyline_length(l);
    }
    e.parent = best < 0.5 ? std::max(parent, e.length) : e.length;
  }
  for (std::size_t k = 0; k < P.links.size(); k++) {
    double total = 0;
    for (auto e : N.link_edges[k]) total += N.edges[e].length;
    for (auto e : N.link_edges[k]) N.edges[e].parent = total;
    if (N.link_edges[k].empty()) throw std::logic_error("flow::superset: a candidate link with no segment");
  }
  return N;
}

// ---- walks among the directed segments -------------------------------------------------------------------------------------

// ω for each edge at θ (both directions alike).
inline std::vector<double> weights(const Network& N, double theta, const Widths& W) {
  std::vector<double> w(N.edges.size());
  for (std::size_t e = 0; e < N.edges.size(); e++) {
    const auto& x = N.edges[e];
    w[e] = std::exp(-theta * x.length) * std::pow(x.width / W.reference, x.length / x.parent);
  }
  return w;
}

// The degree of each junction with the links made.
inline std::vector<std::size_t> degrees(const Network& N, const hv::Bits& made) {
  std::vector<std::size_t> d(N.nodes(), 0);
  for (std::size_t e = 0; e < N.edges.size(); e++)
    if (N.present(e, made)) d[N.edges[e].u]++, d[N.edges[e].v]++;
  return d;
}

// The transitions M_ab = ω_b between directed segments, with the links made: b leaves the junction a enters, and is not
// a reversed unless that junction is a dead end. A link not made keeps the transitions out of its segments (onto the
// streets as built) but none into them, so that nobody reaches it: then a link made adds only the transitions into it,
// and its rows need no update.
inline la::Matrix transitions(const Network& N, const std::vector<double>& w, const hv::Bits& made) {
  const std::size_t S = N.states();
  const auto deg = degrees(N, made);
  std::vector<std::vector<std::size_t>> leaving(N.nodes()), leaving_built(N.nodes());
  for (std::size_t s = 0; s < S; s++) {
    if (N.present(s / 2, made)) leaving[N.tail(s)].push_back(s);
    if (N.edges[s / 2].link < 0) leaving_built[N.tail(s)].push_back(s);
  }
  la::Matrix M(S);
  for (std::size_t a = 0; a < S; a++) {
    const std::size_t x = N.head(a);
    if (N.present(a / 2, made)) {
      for (auto b : leaving[x])
        if (b != (a ^ 1) || deg[x] == 1) M(a, b) = w[b / 2];
    } else {
      for (auto b : leaving_built[x]) M(a, b) = w[b / 2];
    }
  }
  return M;
}

struct Walks {
  std::size_t n = 0;
  std::vector<double> z;  // row-major, over the directed segments
  double operator()(std::size_t a, std::size_t b) const { return z[a * n + b]; }
};

inline Walks walks_direct(const Network& N, const std::vector<double>& w, const hv::Bits& made) {
  const std::size_t S = N.states();
  la::Matrix A = transitions(N, w, made);
  for (std::size_t i = 0; i < S; i++)
    for (std::size_t j = 0; j < S; j++) A(i, j) = (i == j ? 1.0 : 0.0) - A(i, j);
  const la::Matrix Z = la::solve(A, la::Matrix::identity(S));
  return {S, Z.a};
}

// An estimate of the spectral radius of a non-negative matrix: ‖M^k 1‖^(1/k) for large k.
inline double spectral_radius(const la::Matrix& M, std::size_t iterations = 4000) {
  const std::size_t n = M.n;
  std::vector<double> x(n, 1.0), y(n);
  double log_growth = 0;
  for (std::size_t it = 0; it < iterations; it++) {
    for (std::size_t i = 0; i < n; i++) {
      double s = 0;
      for (std::size_t j = 0; j < n; j++) s += M(i, j) * x[j];
      y[i] = s;
    }
    double norm = 0;
    for (double v : y) norm = std::max(norm, v);
    if (norm == 0) return 0;
    if (it >= iterations / 2) log_growth += std::log(norm);
    for (std::size_t i = 0; i < n; i++) x[i] = y[i] / norm;
  }
  return std::exp(log_growth / static_cast<double>(iterations - iterations / 2));
}

// Z for the links made, from Z0 (none made) by Woodbury. A made link adds, for each of its directed segments s, the
// column of transitions into s (rank one: u = M[:, s], v = e_s), and removes the turning back at a junction it ends
// that was a dead end (rank one). With M = M0 + U Vᵀ, (I − M)⁻¹ = Z0 + Z0 U (I − Vᵀ Z0 U)⁻¹ Vᵀ Z0.
inline Walks walks_updated(const Walks& Z0, const Network& N, const std::vector<double>& w, const hv::Bits& made) {
  bool any = false;
  for (auto b : made) any = any || b;
  if (!any) return Z0;
  const std::size_t S = N.states();
  const hv::Bits none_made(made.size(), 0);
  const auto deg0 = degrees(N, none_made), deg = degrees(N, made);
  std::vector<std::vector<std::size_t>> entering(N.nodes());
  for (std::size_t s = 0; s < S; s++)
    if (N.present(s / 2, made)) entering[N.head(s)].push_back(s);
  // Each rank-one term: u sparse (pairs), v = e_col.
  struct Term { std::vector<std::pair<std::size_t, double>> u; std::size_t col; };
  std::vector<Term> terms;
  for (std::size_t k = 0; k < N.link_edges.size(); k++) {
    if (!made[k]) continue;
    for (auto e : N.link_edges[k])
      for (std::size_t s : {2 * e, 2 * e + 1}) {
        Term t;
        t.col = s;
        const std::size_t x = N.tail(s);
        for (auto a : entering[x])
          if (a != (s ^ 1) || deg[x] == 1) t.u.push_back({a, w[e]});
        terms.push_back(std::move(t));
      }
  }
  // Turning back at a junction that was a dead end and is no longer one.
  for (std::size_t x = 0; x < N.nodes(); x++) {
    if (!(deg0[x] == 1 && deg[x] > 1)) continue;
    for (std::size_t s = 0; s < S; s++)
      if (N.edges[s / 2].link < 0 && N.head(s) == x) terms.push_back({{{s, -w[s / 2]}}, s ^ 1});
  }
  const std::size_t r = terms.size(), n = Z0.n;
  // Z0 U (n × r) and K = I − Vᵀ Z0 U.
  std::vector<double> ZU(n * r, 0.0);
  for (std::size_t q = 0; q < r; q++)
    for (const auto& [a, val] : terms[q].u)
      for (std::size_t i = 0; i < n; i++) ZU[i * r + q] += Z0(i, a) * val;
  la::Matrix K(r);
  for (std::size_t p = 0; p < r; p++)
    for (std::size_t q = 0; q < r; q++) K(p, q) = (p == q ? 1.0 : 0.0) - ZU[terms[p].col * r + q];
  const la::Matrix Kinv = la::solve(K, la::Matrix::identity(r));
  // W = K⁻¹ Vᵀ Z0 (r × n): rows col_q of Z0 mixed by K⁻¹.
  std::vector<double> Wm(r * n, 0.0);
  for (std::size_t p = 0; p < r; p++)
    for (std::size_t q = 0; q < r; q++) {
      const double k = Kinv(p, q);
      if (k == 0) continue;
      const double* row = Z0.z.data() + terms[q].col * n;
      double* out = Wm.data() + p * n;
      for (std::size_t b = 0; b < n; b++) out[b] += k * row[b];
    }
  Walks Z = Z0;
  for (std::size_t a = 0; a < n; a++) {
    double* out = Z.z.data() + a * n;
    for (std::size_t p = 0; p < r; p++) {
      const double c = ZU[a * r + p];
      if (c == 0) continue;
      const double* row = Wm.data() + p * n;
      for (std::size_t b = 0; b < n; b++) out[b] += c * row[b];
    }
  }
  return Z;
}

// The trips of one set of links: Z, and what the formula of a trip's flow needs at its ends (README §2.2).
class Trips {
 public:
  const Network& N;
  hv::Bits made;
  Walks Z;
  std::vector<double> w;
  std::vector<std::vector<std::size_t>> leaving, entering;  // per junction, present directed segments

  Trips(const Network& net, const std::vector<double>& weights, const hv::Bits& m, Walks z) : N(net), made(m), Z(std::move(z)), w(weights) {
    leaving.assign(N.nodes(), {}), entering.assign(N.nodes(), {});
    for (std::size_t s = 0; s < N.states(); s++)
      if (N.present(s / 2, made)) leaving[N.tail(s)].push_back(s), entering[N.head(s)].push_back(s);
  }

  // Per destination d: Q = Z[A, A]⁻¹ and its row sums q = Q 1.
  struct Destination {
    std::vector<std::size_t> A;
    std::vector<double> Q, q;
  };
  const Destination& destination(std::size_t d) const {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = dest_.find(d);
    if (it != dest_.end()) return it->second;
    Destination D;
    D.A = entering[d];
    const std::size_t r = D.A.size();
    if (r) {
      la::Matrix B(r);
      for (std::size_t i = 0; i < r; i++)
        for (std::size_t j = 0; j < r; j++) B(i, j) = Z(D.A[i], D.A[j]);
      const auto Qm = la::solve(B, la::Matrix::identity(r));
      D.Q = Qm.a;
      D.q.assign(r, 0.0);
      for (std::size_t i = 0; i < r; i++)
        for (std::size_t j = 0; j < r; j++) D.q[i] += Qm(i, j);
    }
    return dest_.emplace(d, std::move(D)).first->second;
  }
  // Y[o, c] = Σ_{b leaving o} ω_b Z[b, c].
  double Y(std::size_t o, std::size_t c) const {
    double s = 0;
    for (auto b : leaving[o]) s += w[b / 2] * Z(b, c);
    return s;
  }
  // σ_d(c) = Z[c, A_d] q_d: the weight of the walks from c to d.
  double sigma(const Destination& D, std::size_t c) const {
    double s = 0;
    for (std::size_t i = 0; i < D.A.size(); i++) s += Z(c, D.A[i]) * D.q[i];
    return s;
  }
  // (Y[o, A] Q): the absorbed weights, and W_od their sum.
  std::vector<double> yq(std::size_t o, const Destination& D) const {
    const std::size_t r = D.A.size();
    std::vector<double> y(r), out(r, 0.0);
    for (std::size_t i = 0; i < r; i++) y[i] = Y(o, D.A[i]);
    for (std::size_t j = 0; j < r; j++)
      for (std::size_t i = 0; i < r; i++) out[j] += y[i] * D.Q[i * r + j];
    return out;
  }
  // π_od(c) = Y[o, c] − (Y[o, A] Q)·(Z[A, c] − δ_{A c}).
  double pi(const Destination& D, const std::vector<double>& yqv, std::size_t c, double y_oc) const {
    double s = y_oc;
    for (std::size_t i = 0; i < D.A.size(); i++) s -= yqv[i] * (Z(D.A[i], c) - (D.A[i] == c ? 1.0 : 0.0));
    return s;
  }
  // The mean visits of one trip o → d to the directed segment c.
  double visits(std::size_t o, std::size_t d, std::size_t c) const {
    if (o == d) return 0;
    const auto& D = destination(d);
    const auto q = yq(o, D);
    double W = 0;
    for (double v : q) W += v;
    if (!(W > 0)) return 0;
    return pi(D, q, c, Y(o, c)) * sigma(D, c) / W;
  }

 private:
  mutable std::mutex mutex_;
  mutable std::map<std::size_t, Destination> dest_;
};

// ---- the model ----------------------------------------------------------------------------------------------------------

struct Options {
  double theta = 0.02;          // 1/m: the route choice's sensitivity to length (calibrated over a grid)
  Widths widths;
  double entrance_edge = 10;    // m: a dead end this near the box's edge is an entrance
  double same_street = 400;     // m: an entrance's counter on its own street within this
  double near_counter = 120;    // m: else, for a street or mall, any counter within this
  double local_length = 200;    // m: λ of the local trips (assumption)
  std::array<double, economy::kinds> attraction{0, 1, 1, 0.5, 0, 1};  // w_k: visits drawn by an open establishment (assumption)
  double visit_share = 0.1;     // ν: the visits' share of the hour's trips where the counters cannot tell it (assumption)
  bool handedness = true;       // fit each footpath's share of each way (ρ)
};

struct Entrance {
  std::size_t node = 0;
  std::string street, kind;
  double width = 0;
  std::array<double, H> gross{};   // people an hour, both ways
  std::string source;              // "same street", "near", "default"
  long long counter = -1;
  double distance = 0;
};

struct Counter {
  long long id = 0;
  std::string name, street;
  Point at{};
  std::size_t edge = 0;
  double share = 0.5;
  std::array<double, H> week{};                 // the mean week, both ways (the problem's sensors.json)
  bool directed = false;                        // each way read and placed on the segment's directions
  std::array<std::string, 2> toward;            // the City's two directions
  std::array<std::size_t, 2> state{};           // the directed segment (2e or 2e + 1) each way runs along
  std::array<std::array<double, H>, 2> way{};   // the mean count each way at each hour of the week (urban::flow::current_meter)
  double cosine = 0;                            // of the first way's compass direction with the segment's u → v
  double hand = 0;                              // ρ: the footpath carries (1 + ρ)/2 of the first way's walkers and (1 − ρ)/2 of the second's
  std::array<double, H> gross() const {
    std::array<double, H> g{};
    for (std::size_t h = 0; h < H; h++) g[h] = directed ? way[0][h] + way[1][h] : week[h];
    return g;
  }
};

// The patterns of trips (README §2.2) and their names.
enum Mode : std::size_t { arrive_work = 0, arrive_home, leave_work, leave_home, through, local, visit, modes };
inline const char* mode_name(std::size_t m) {
  static const char* n[] = {"arrive (work)", "arrive (home)", "leave (work)", "leave (home)", "through", "local", "visit"};
  return n[m];
}

// The unit responses of every family of trips on a list of edges: the flow (u → v, v → u) on each.
struct Terms {
  std::vector<std::size_t> edges;
  using Flows = std::vector<std::array<double, 2>>;
  std::array<std::vector<Flows>, 4> arrive_leave;  // [arrive_work, arrive_home, leave_work, leave_home][entrance]
  std::vector<Flows> through;                       // [i * E + j]
  Flows local;
  std::array<Flows, economy::kinds> context;        // per kind: the visits to (and from) the context's open attractors, w_k included
  std::vector<Flows> frontage;                      // per frontage: the visits to it and back, per unit weight
};

// Per set of links: what the search needs (README §3).
struct Response {
  std::vector<double> B;  // n × H: passers from everything but the plan's own frontages
  std::vector<double> R;  // n × n: σ_f × the gross flow on f's segment of one visit an hour to frontage g and back
  std::shared_ptr<const economy::Distances> distances;
};

class Model {
 public:
  using Amplitudes = std::array<std::array<double, H>, modes>;
  const economy::Problem& P;
  Options o;
  Network net;
  std::vector<double> omega;
  Walks Z0;
  double radius_none = 0, radius_all = 0;
  std::vector<std::size_t> component;     // per junction, its component of the network as built (0: the largest)
  std::vector<Entrance> entrances;
  std::vector<Counter> counters;          // inside the box: the calibration's
  std::vector<double> work, home;         // floor per node, m² (offices and shops; dwellings)
  std::array<std::vector<double>, economy::kinds> attract;  // the context's establishments per node, by kind (counts)
  std::array<std::array<char, H>, economy::kinds> open_kind{};  // the kinds' typical hours as built
  std::vector<double> origin;             // visits' origins π, summing to 1
  std::vector<std::size_t> frontage_edge;
  std::vector<std::array<std::pair<std::size_t, double>, 2>> frontage_nodes;
  std::vector<double> frontage_share;
  double W_ref = 0;                        // the weight of every attractor as built
  Amplitudes x{};                          // the fitted amplitudes
  std::array<double, H> gsum{};
  std::vector<std::size_t> building_nodes, attractor_nodes, origin_nodes;

  Model(const economy::Problem& problem, const Json& inputs, const Json& sensors, Options options) : P(problem), o(options) {
    net = superset(P, o.widths);
    components();
    read_terminals(inputs, sensors);
    set_theta(o.theta);
  }

  void set_theta(double theta) {
    o.theta = theta;
    omega = weights(net, theta, o.widths);
    Z0 = walks_direct(net, omega, hv::Bits(P.links.size(), 0));
    radius_none = spectral_radius(transitions(net, omega, hv::Bits(P.links.size(), 0)), 600);
    radius_all = spectral_radius(transitions(net, omega, hv::Bits(P.links.size(), 1)), 600);
  }

  Walks walks(const hv::Bits& made) const { return walks_updated(Z0, net, omega, made); }
  Trips trips(const hv::Bits& made) const { return Trips(net, omega, made, walks(made)); }

  void components() {
    const std::size_t n = net.nodes();
    std::vector<std::vector<std::size_t>> adj(n);
    for (const auto& e : net.edges)
      if (e.link < 0) adj[e.u].push_back(e.v), adj[e.v].push_back(e.u);
    std::vector<std::size_t> comp(n, none), size;
    for (std::size_t r = 0; r < n; r++) {
      if (comp[r] != none) continue;
      std::vector<std::size_t> st{r};
      comp[r] = size.size();
      std::size_t c = 0;
      while (!st.empty()) {
        const auto j = st.back();
        st.pop_back();
        c++;
        for (auto k : adj[j])
          if (comp[k] == none) comp[k] = size.size(), st.push_back(k);
      }
      size.push_back(c);
    }
    const std::size_t big = static_cast<std::size_t>(std::max_element(size.begin(), size.end()) - size.begin());
    component.assign(n, 1);
    for (std::size_t j = 0; j < n; j++) component[j] = comp[j] == big ? 0 : 1;
  }
  bool on_main(std::size_t e) const { return component[net.edges[e].u] == 0; }

  // The nearest edge of the network as built (no link's), on its largest component, to a point; its arc length and
  // distance. With a street's name, only that street's edges.
  std::tuple<std::size_t, double, double> nearest_edge(const Point& p, const std::string& street = "") const {
    std::size_t best = none;
    double bs = 0, bd = std::numeric_limits<double>::infinity();
    for (std::size_t e = 0; e < net.edges.size(); e++) {
      if (net.edges[e].link >= 0 || !on_main(e)) continue;
      if (!street.empty() && net.edges[e].name != street) continue;
      const auto [s, d] = cg::nearest_on_polyline(net.edges[e].line, p);
      if (d < bd) bd = d, bs = s, best = e;
    }
    return {best, bs, bd};
  }

  void read_terminals(const Json& in, const Json& sj) {
    const std::size_t n = net.nodes(), F = P.size();
    const Json& box = in["district"]["district_box"];
    const double x0 = box["x0"].as_number(), x1 = box["x1"].as_number(), y0 = box["y0"].as_number(), y1 = box["y1"].as_number();
    const auto inside = [&](const Point& p) { return p[0] >= x0 && p[0] <= x1 && p[1] >= y0 && p[1] <= y1; };
    struct Raw { long long id; std::string name, street, kind; Point at; std::array<double, H> week; bool crossing; };
    std::vector<Raw> raw;
    for (const auto& s : sj["sensors"].as_array()) {
      Raw r;
      r.id = static_cast<long long>(s["id"].as_number()), r.name = s["name"].text(), r.street = s["way"]["name"].text(), r.kind = s["way"]["kind"].text();
      r.at = {s["x"].as_number(), s["y"].as_number()};
      const auto& w = s["mean_week"].as_array();
      for (std::size_t h = 0; h < H; h++) r.week[h] = std::max(0.0, w[h].as_number());
      r.crossing = r.name.find("Crossing") != std::string::npos;
      raw.push_back(std::move(r));
    }
    // Each way of the counters inside the box: the City's directional counts over the reference year, public holidays
    // left out as the problem's tide leaves them, read as a current meter by the urban-science domain.
    std::map<long long, const Json*> inputs_of;
    for (const auto& s : in["sensors"].as_array()) inputs_of[static_cast<long long>(s["id"].as_number())] = &s;
    std::set<long long> holidays;
    {
      const Json Hj = Json::parse(read_text(in["files"]["holidays"].text()));
      for (const auto& h : Hj["holidays"].as_array()) holidays.insert(static_cast<long long>(std::floor(urban::rhythm::parse_utc(h["date"].text()) / urban::rhythm::day_s)));
    }
    const double ref0 = urban::rhythm::parse_utc(in["reference"]["start_utc"].text());
    const double ref1 = ref0 + in["reference"]["weeks"].as_number() * 168 * 3600.0;
    const auto compass = [](const std::string& c) -> Point {
      if (c == "North") return {0, 1};
      if (c == "South") return {0, -1};
      if (c == "East") return {1, 0};
      return {-1, 0};
    };
    for (const auto& r : raw) {
      if (!inside(r.at) || r.crossing) continue;
      Counter c;
      c.id = r.id, c.name = r.name, c.street = r.street, c.at = r.at, c.week = r.week;
      auto [e, s, d] = nearest_edge(r.at, r.street);
      if (e == none || d > 40) std::tie(e, s, d) = nearest_edge(r.at);
      c.edge = e, c.share = side_share(net.edges[e].kind);
      const Json* si = inputs_of.count(r.id) ? inputs_of.at(r.id) : nullptr;
      if (si && (*si)["directions"].is_array()) {
        const Json series = Json::parse(read_text((*si)["series"].text()));
        if (series.has("directions")) {
          const Json& dj = series["directions"];
          urban::rhythm::Record d1, d2;
          d1.id = d2.id = std::to_string(r.id);
          d1.t0 = d2.t0 = urban::rhythm::parse_utc(dj["t0_utc"].text());
          d1.step = d2.step = dj["step_s"].as_number();
          for (const auto& o : series["offsets"].as_array()) d1.offsets.push_back({urban::rhythm::parse_utc(o[0].text()), o[1].as_number()});
          d2.offsets = d1.offsets;
          for (const auto& v : dj["d1"].as_array()) d1.y.push_back(v.is_null() ? std::nan("") : v.as_number());
          for (const auto& v : dj["d2"].as_array()) d2.y.push_back(v.is_null() ? std::nan("") : v.as_number());
          for (std::size_t i = 0; i < d1.size(); i++)
            if (holidays.count(d1.civil_day(d1.time(i)))) d1.y[i] = d2.y[i] = std::nan("");
          const auto cm = urban::flow::current_meter(d1, d2, ref0, ref1);
          for (std::size_t h = 0; h < H; h++) c.way[0][h] = cm.way[0][h], c.way[1][h] = cm.way[1][h];
          c.toward = {(*si)["directions"][0].text(), (*si)["directions"][1].text()};
          // The segment's direction where the counter stands, against the first way's compass direction.
          const auto& l = net.edges[e].line;
          const Point a = point_along(l, std::max(0.0, s - 3)), b = point_along(l, std::min(net.edges[e].length, s + 3));
          const double tx = b[0] - a[0], ty = b[1] - a[1], tn = std::hypot(tx, ty);
          const Point cv = compass(c.toward[0]);
          c.cosine = tn > 0 ? (tx * cv[0] + ty * cv[1]) / tn : 0;
          if (std::fabs(c.cosine) >= 0.5) {
            c.directed = true;
            c.state = c.cosine > 0 ? std::array<std::size_t, 2>{2 * e, 2 * e + 1} : std::array<std::size_t, 2>{2 * e + 1, 2 * e};
          }
        }
      }
      counters.push_back(std::move(c));
    }
    // Entrances: dead ends of the network as built on the edge of its box, on its largest component.
    std::vector<std::size_t> degree(n, 0);
    for (const auto& e : net.edges)
      if (e.link < 0) degree[e.u]++, degree[e.v]++;
    std::array<double, H> median{};
    {
      std::vector<std::array<double, H>> ext;
      for (const auto& r : raw)
        if (!inside(r.at) && !r.crossing && r.kind == "street") {
          std::array<double, H> g{};
          for (std::size_t h = 0; h < H; h++) g[h] = r.week[h] / side_share("street");
          ext.push_back(g);
        }
      for (std::size_t h = 0; h < H; h++) {
        std::vector<double> v;
        for (const auto& g : ext) v.push_back(g[h]);
        std::sort(v.begin(), v.end());
        median[h] = v.empty() ? 0 : (v.size() % 2 ? v[v.size() / 2] : 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]));
      }
    }
    for (std::size_t j = 0; j < n; j++) {
      if (degree[j] != 1 || component[j] != 0) continue;
      const Point& p = net.at[j];
      const double dd = std::min({std::fabs(p[0] - x0), std::fabs(p[0] - x1), std::fabs(p[1] - y0), std::fabs(p[1] - y1)});
      if (dd > o.entrance_edge) continue;
      Entrance E;
      E.node = j;
      for (const auto& e : net.edges)
        if (e.link < 0 && (e.u == j || e.v == j)) E.street = e.name, E.kind = e.kind, E.width = e.width;
      const Raw* best = nullptr;
      double bd = std::numeric_limits<double>::infinity(), scale = 1;
      for (const auto& r : raw)
        if (!inside(r.at) && !r.crossing && !E.street.empty() && r.street == E.street) {
          const double d = std::hypot(r.at[0] - p[0], r.at[1] - p[1]);
          if (d < bd && d <= o.same_street) bd = d, best = &r;
        }
      if (best) E.source = "same street";
      else if (E.kind != "laneway")
        for (const auto& r : raw)
          if (!inside(r.at) && !r.crossing) {
            const double d = std::hypot(r.at[0] - p[0], r.at[1] - p[1]);
            if (d < bd && d <= o.near_counter) bd = d, best = &r, E.source = "near";
          }
      if (best) {
        const std::string ck = best->kind == "laneway" ? "laneway" : (best->kind == "mall" ? "mall" : "street");
        if (E.source == "near") scale = E.width / width_of(o.widths, ck, best->street);
        for (std::size_t h = 0; h < H; h++) E.gross[h] = scale * best->week[h] / side_share(ck);
        E.counter = best->id, E.distance = bd;
      } else {
        E.source = "default";
        for (std::size_t h = 0; h < H; h++) E.gross[h] = median[h] * E.width / o.widths.main;
      }
      entrances.push_back(std::move(E));
    }
    for (std::size_t h = 0; h < H; h++) {
      gsum[h] = 0;
      for (const auto& E : entrances) gsum[h] += E.gross[h];
    }
    // Buildings: floor on the two ends of the nearest segment.
    work.assign(n, 0.0), home.assign(n, 0.0);
    const auto attach = [&](const Point& p, double mass, std::vector<double>& to) {
      const auto [e, s, d] = nearest_edge(p);
      (void)d;
      const auto& x = net.edges[e];
      const double t = x.length > 0 ? std::clamp(s / x.length, 0.0, 1.0) : 0.5;
      to[x.u] += mass * (1 - t), to[x.v] += mass * t;
    };
    const Json& D = in["district"];
    for (const char* group : {"parcels", "context"})
      for (const auto& b : D[group].as_array()) {
        if (std::string(group) == "context" && !b["in_district_box"].as_bool()) continue;
        const double floors = b["floors"].is_null() ? std::max(1.0, std::round(b["height_m"].as_number() / 3.2)) : b["floors"].as_number();
        const double floor = b["area_m2"].as_number() * floors;
        const Point c = economy::detail::centroid(economy::detail::ring_of(b["footprint"]));
        attach(c, floor, static_cast<int>(b["use"].as_number()) == 0 ? home : work);
      }
    // The kinds' typical hours as built (the problem's), for the context's establishments.
    for (std::size_t k = 1; k < economy::kinds; k++) {
      const Json& t = in["typical_hours"][economy::kind_name(static_cast<std::int64_t>(k))];
      for (std::size_t h = 0; h < H; h++) {
        const Json& wd = t[economy::is_weekend(h) ? "weekend" : "weekday"];
        const double hd = static_cast<double>(economy::hour_of_day(h));
        open_kind[k][h] = !wd.is_null() && hd >= wd[0].as_number() && hd < wd[1].as_number();
      }
    }
    // Attractors: CLUE's establishments of the trading kinds in the box, outside the 61 parcels.
    for (auto& a : attract) a.assign(n, 0.0);
    for (const auto& r : in["establishments"]["rows"].as_array()) {
      const Point p{r[0].as_number(), r[1].as_number()};
      const auto k = static_cast<std::size_t>(r[2].as_number());
      if (!economy::trades(static_cast<std::int64_t>(k)) || !inside(p)) continue;
      bool in_parcel = false;
      for (const auto& f : P.frontages)
        if (cg::inside_ring(f.footprint, p)) { in_parcel = true; break; }
      if (in_parcel) continue;
      attach(p, 1.0, attract[k]);
    }
    // The frontages: on the nearest segment of the streets as built, as the problem places them.
    for (std::size_t f = 0; f < F; f++) {
      const auto [e, s, d] = nearest_edge(P.frontages[f].anchor);
      (void)d;
      const auto& x = net.edges[e];
      const double t = x.length > 0 ? std::clamp(s / x.length, 0.0, 1.0) : 0.5;
      frontage_edge.push_back(e);
      frontage_nodes.push_back({std::pair{x.u, 1 - t}, std::pair{x.v, t}});
      frontage_share.push_back(side_share(x.kind));
    }
    // The visits' origins: half the entrances (by their mean flow), half the buildings (by floor).
    origin.assign(n, 0.0);
    double gs = 0, ms = 0;
    for (const auto& E : entrances)
      for (double g : E.gross) gs += g;
    for (std::size_t j = 0; j < n; j++) ms += work[j] + home[j];
    for (const auto& E : entrances) {
      double g = 0;
      for (double v : E.gross) g += v;
      origin[E.node] += 0.5 * g / gs;
    }
    for (std::size_t j = 0; j < n; j++) origin[j] += 0.5 * (work[j] + home[j]) / ms;
    W_ref = 0;
    for (std::size_t k = 0; k < economy::kinds; k++)
      for (std::size_t j = 0; j < n; j++) W_ref += o.attraction[k] * attract[k][j];
    for (const auto& f : P.frontages) W_ref += o.attraction[static_cast<std::size_t>(f.kind0)];
    for (std::size_t j = 0; j < n; j++) {
      if (work[j] + home[j] > 0) building_nodes.push_back(j);
      double a = 0;
      for (std::size_t k = 0; k < economy::kinds; k++) a += attract[k][j];
      if (a > 0) attractor_nodes.push_back(j);
      if (origin[j] > 0) origin_nodes.push_back(j);
    }
  }

  // ---- the unit responses ------------------------------------------------------------------------------------------

  // The flows of every family of trips on the given edges, with the links of `T`.
  Terms terms(const Trips& T, const std::vector<std::size_t>& edges) const {
    const std::size_t n = net.nodes(), E = entrances.size(), Ne = edges.size(), F = P.size();
    Terms R;
    R.edges = edges;
    const auto zero = [&] { return Terms::Flows(Ne, {0.0, 0.0}); };
    // The directed segments asked for, and every terminal's Y[o, c] on them.
    std::vector<std::size_t> cs;
    for (auto e : edges) cs.push_back(2 * e), cs.push_back(2 * e + 1);
    const std::size_t C = cs.size();
    std::set<std::size_t> nodes_set(building_nodes.begin(), building_nodes.end());
    for (auto j : attractor_nodes) nodes_set.insert(j);
    for (auto j : origin_nodes) nodes_set.insert(j);
    for (const auto& en : entrances) nodes_set.insert(en.node);
    for (const auto& fn : frontage_nodes)
      for (const auto& [j, wg] : fn)
        if (wg > 0) nodes_set.insert(j);
    const std::vector<std::size_t> nodes(nodes_set.begin(), nodes_set.end());
    std::vector<std::size_t> index(n, none);
    for (std::size_t i = 0; i < nodes.size(); i++) index[nodes[i]] = i;
    const std::size_t Nn = nodes.size();
    std::vector<double> Yc(Nn * C);
    for (std::size_t i = 0; i < Nn; i++)
      for (std::size_t k = 0; k < C; k++) Yc[i * C + k] = T.Y(nodes[i], cs[k]);
    // Per destination: σ_d(c) and Z[A_d, c] − δ.
    struct Dest { const Trips::Destination* D; std::vector<double> sigma, G; };
    std::vector<Dest> dest(Nn);
    for (std::size_t i = 0; i < Nn; i++) {
      auto& d = dest[i];
      d.D = &T.destination(nodes[i]);
      const std::size_t r = d.D->A.size();
      d.sigma.assign(C, 0.0), d.G.assign(C * r, 0.0);
      for (std::size_t k = 0; k < C; k++) {
        d.sigma[k] = T.sigma(*d.D, cs[k]);
        for (std::size_t a = 0; a < r; a++) d.G[k * r + a] = T.Z(d.D->A[a], cs[k]) - (d.D->A[a] == cs[k] ? 1.0 : 0.0);
      }
    }
    // Y[o, A_d] Q_d and W_od for every terminal pair, computed on demand and kept.
    std::vector<std::vector<double>> YQ(Nn * Nn);
    std::vector<double> Wod(Nn * Nn, -1.0);
    const auto pair = [&](std::size_t oi, std::size_t di) -> const std::vector<double>& {
      auto& v = YQ[oi * Nn + di];
      if (Wod[oi * Nn + di] < 0) {
        v = T.yq(nodes[oi], *dest[di].D);
        double W = 0;
        for (double q : v) W += q;
        Wod[oi * Nn + di] = W;
      }
      return v;
    };
    // The flow on segment k of a trip oi → di, given its YQ and W.
    const auto trip = [&](std::size_t oi, std::size_t di, std::size_t k) {
      if (oi == di) return 0.0;
      const auto& q = pair(oi, di);
      const double W = Wod[oi * Nn + di];
      if (!(W > 0)) return 0.0;
      const auto& d = dest[di];
      const std::size_t r = d.D->A.size();
      double p = Yc[oi * C + k];
      for (std::size_t a = 0; a < r; a++) p -= q[a] * d.G[k * r + a];
      return p * d.sigma[k] / W;
    };
    // A pattern: Σ T_od × trip(o, d) over (origin, destination, weight) triples.
    const auto pattern = [&](const std::vector<std::tuple<std::size_t, std::size_t, double>>& od) {
      Terms::Flows fl = zero();
      for (const auto& [oi, di, t] : od) {
        if (t == 0) continue;
        for (std::size_t k = 0; k < C; k++) fl[k / 2][k % 2] += t * trip(oi, di, k);
      }
      return fl;
    };
    // Normalised masses.
    double sw = 0, sh = 0;
    for (std::size_t j = 0; j < n; j++) sw += work[j], sh += home[j];
    for (std::size_t q = 0; q < 4; q++) R.arrive_leave[q].assign(E, zero());
    for (std::size_t i = 0; i < E; i++) {
      const std::size_t oi = index[entrances[i].node];
      std::vector<std::tuple<std::size_t, std::size_t, double>> aw, ah, lw, lh;
      for (auto b : building_nodes) {
        const std::size_t bi = index[b];
        if (work[b] > 0) aw.push_back({oi, bi, work[b] / sw}), lw.push_back({bi, oi, work[b] / sw});
        if (home[b] > 0) ah.push_back({oi, bi, home[b] / sh}), lh.push_back({bi, oi, home[b] / sh});
      }
      R.arrive_leave[0][i] = pattern(aw), R.arrive_leave[1][i] = pattern(ah);
      R.arrive_leave[2][i] = pattern(lw), R.arrive_leave[3][i] = pattern(lh);
    }
    R.through.assign(E * E, zero());
    for (std::size_t i = 0; i < E; i++)
      for (std::size_t j = 0; j < E; j++)
        if (i != j) R.through[i * E + j] = pattern({{index[entrances[i].node], index[entrances[j].node], 1.0}});
    {
      std::vector<std::tuple<std::size_t, std::size_t, double>> od;
      double norm = 0;
      for (auto a : building_nodes)
        for (auto b : building_nodes) {
          if (a == b) continue;
          const double t = (work[a] + home[a]) * (work[b] + home[b]) * std::exp(-std::hypot(net.at[a][0] - net.at[b][0], net.at[a][1] - net.at[b][1]) / o.local_length);
          od.push_back({index[a], index[b], t}), norm += t;
        }
      for (auto& [a, b, t] : od) t /= norm;
      R.local = pattern(od);
    }
    // Visits: to a node j from the origins π and back, per unit.
    std::set<std::size_t> visit_nodes(attractor_nodes.begin(), attractor_nodes.end());
    for (const auto& fn : frontage_nodes)
      for (const auto& [j, wg] : fn)
        if (wg > 0) visit_nodes.insert(j);
    std::map<std::size_t, Terms::Flows> per_node;
    for (auto j : visit_nodes) {
      std::vector<std::tuple<std::size_t, std::size_t, double>> od;
      for (auto oo : origin_nodes)
        if (oo != j) od.push_back({index[oo], index[j], origin[oo]}), od.push_back({index[j], index[oo], origin[oo]});
      per_node.emplace(j, pattern(od));
    }
    for (std::size_t q = 0; q < economy::kinds; q++) {
      R.context[q] = zero();
      if (o.attraction[q] <= 0) continue;
      for (auto j : attractor_nodes) {
        const double a = attract[q][j] * o.attraction[q];
        if (a <= 0) continue;
        const auto& fl = per_node.at(j);
        for (std::size_t k = 0; k < Ne; k++) R.context[q][k][0] += a * fl[k][0], R.context[q][k][1] += a * fl[k][1];
      }
    }
    R.frontage.assign(F, zero());
    for (std::size_t f = 0; f < F; f++)
      for (const auto& [j, wg] : frontage_nodes[f]) {
        if (wg <= 0) continue;
        const auto& fl = per_node.at(j);
        for (std::size_t k = 0; k < Ne; k++) R.frontage[f][k][0] += wg * fl[k][0], R.frontage[f][k][1] += wg * fl[k][1];
      }
    return R;
  }

  // ---- the hour's directed flows on the terms' edges ----------------------------------------------------------------

  // draw[g * H + h]: the weight with which frontage g draws visits at h (its kind's w_k when open, else 0).
  std::vector<double> frontage_draw(const economy::Plan& plan) const {
    const std::size_t F = P.size();
    std::vector<double> d(F * H, 0.0);
    for (std::size_t g = 0; g < F; g++)
      for (std::size_t h = 0; h < H; h++)
        if (plan.open(g, h)) d[g * H + h] = o.attraction[static_cast<std::size_t>(plan.kind[g])];
    return d;
  }

  // The directed flows (u → v, v → u) on the terms' edges at hour h, by mode (the plan's own visits within `visit`).
  std::array<Terms::Flows, modes> mode_flows(const Terms& T, const std::vector<double>& draw, std::size_t h) const { return mode_flows(T, draw, h, x); }
  std::array<Terms::Flows, modes> mode_flows(const Terms& T, const std::vector<double>& draw, std::size_t h, const Amplitudes& a) const {
    const std::size_t Ne = T.edges.size(), E = entrances.size(), F = P.size();
    std::array<Terms::Flows, modes> out;
    for (auto& f : out) f.assign(Ne, {0.0, 0.0});
    const auto add = [&](Terms::Flows& to, const Terms::Flows& from, double s) {
      if (s == 0) return;
      for (std::size_t k = 0; k < Ne; k++) to[k][0] += s * from[k][0], to[k][1] += s * from[k][1];
    };
    for (std::size_t i = 0; i < E; i++) {
      const double g = entrances[i].gross[h];
      add(out[arrive_work], T.arrive_leave[0][i], a[arrive_work][h] * g);
      add(out[arrive_home], T.arrive_leave[1][i], a[arrive_home][h] * g);
      add(out[leave_work], T.arrive_leave[2][i], a[leave_work][h] * g);
      add(out[leave_home], T.arrive_leave[3][i], a[leave_home][h] * g);
    }
    if (gsum[h] > 0 && a[through][h] != 0)
      for (std::size_t i = 0; i < E; i++)
        for (std::size_t j = 0; j < E; j++)
          if (i != j) add(out[through], T.through[i * E + j], a[through][h] * entrances[i].gross[h] * entrances[j].gross[h] / gsum[h]);
    add(out[local], T.local, a[local][h]);
    const double c = a[visit][h] / W_ref;
    for (std::size_t q = 0; q < economy::kinds; q++)
      if (open_kind[q][h]) add(out[visit], T.context[q], c);
    for (std::size_t g = 0; g < F; g++) add(out[visit], T.frontage[g], c * draw[g * H + h]);
    return out;
  }

  // ---- calibration -----------------------------------------------------------------------------------------------------

  // The columns of the hour's least squares at the counters, per unit amplitude: σ_c × the gross flow, and σ_c × the
  // flow along each of the counter's two ways.
  std::vector<std::array<std::array<double, modes>, 3>> columns(const Terms& T, const std::vector<double>& draw, std::size_t h) const {
    Amplitudes unit{};
    for (std::size_t m = 0; m < modes; m++) unit[m][h] = 1;
    const auto fl = mode_flows(T, draw, h, unit);
    std::vector<std::array<std::array<double, modes>, 3>> cols(counters.size());
    for (std::size_t c = 0; c < counters.size(); c++) {
      const auto& k = counters[c];
      const auto j = static_cast<std::size_t>(std::find(T.edges.begin(), T.edges.end(), k.edge) - T.edges.begin());
      if (j == T.edges.size()) throw std::logic_error("flow::Model::columns: a counter's edge is not among the terms'");
      for (std::size_t m = 0; m < modes; m++) {
        cols[c][0][m] = k.share * (fl[m][j][0] + fl[m][j][1]);
        for (std::size_t w = 0; w < 2; w++) cols[c][1 + w][m] = k.share * fl[m][j][k.state[w] % 2];
      }
    }
    return cols;
  }

  // The trips an hour of each mode per unit of its amplitude, as built.
  std::array<double, modes> trips_per_unit(std::size_t h, const std::vector<double>& draw) const {
    std::array<double, modes> t{};
    double g = 0, gg = 0;
    for (const auto& E : entrances) g += E.gross[h];
    for (std::size_t i = 0; i < entrances.size(); i++)
      for (std::size_t j = 0; j < entrances.size(); j++)
        if (i != j && gsum[h] > 0) gg += entrances[i].gross[h] * entrances[j].gross[h] / gsum[h];
    t[arrive_work] = t[arrive_home] = t[leave_work] = t[leave_home] = g;
    t[through] = gg, t[local] = 1;
    double open = 0;
    for (std::size_t q = 0; q < economy::kinds; q++)
      if (open_kind[q][h])
        for (std::size_t j = 0; j < net.nodes(); j++) open += o.attraction[q] * attract[q][j];
    for (std::size_t f = 0; f < P.size(); f++) open += draw[f * H + h];
    t[visit] = 2 * open / W_ref;  // to an attractor and back
    return t;
  }

  struct Fit {
    Amplitudes x{};
    std::vector<std::array<double, H>> fitted, loo;             // per counter: gross
    std::vector<std::array<std::array<double, H>, 2>> fitted_way, loo_way;  // per counter: each way
    double r2 = 0, loo_r2 = 0, r2_log = 0, loo_r2_log = 0;      // of the gross, every counter
    double net_r2 = 0, loo_net_r2 = 0;                          // of the net current u = way₀ − way₁, the counters read each way
    std::vector<double> counter_r2, counter_loo_r2, hand;
    double visit_share = 0;                                     // the visits' mean share of the trips, as fitted or set
  };

  // Fits the amplitudes hour by hour by NNLS (tyche) to each way of the counters read each way and to the gross of the
  // others, with every counter, and leaving each out in turn. `use` masks the modes allowed (all by default).
  //   The visits: the counters cannot tell visits from the other trips (their amplitude is fitted to zero), so with
  //     ν = visit_share > 0 the visits' amplitude is set each hour to make them ν of the hour's trips, and the other
  //     modes are fitted to what is left (a fixed point, four rounds).
  //   Handedness: a footpath counter sees (1 + ρ_c)/2 of one way's walkers and (1 − ρ_c)/2 of the other's (walkers keep
  //     to a side: Street Oceanography's handedness, urban::flow::handedness); ρ_c of each counter on a street is fitted
  //     by least squares over the week, alternating with the amplitudes (three rounds). Left out, a counter's gross is
  //     predicted with ρ = 0.
  Fit calibrate(const Terms& T, bool apply = true, std::array<char, modes> use = {1, 1, 1, 1, 1, 1, 1}) {
    const std::size_t C = counters.size();
    const std::vector<double> draw = frontage_draw(P.as_built());
    const bool fixed_visits = o.visit_share > 0 && use[visit];
    Fit fit;
    fit.fitted.assign(C, {}), fit.loo.assign(C, {}), fit.fitted_way.assign(C, {}), fit.loo_way.assign(C, {});
    std::vector<std::size_t> cols_used;
    for (std::size_t m = 0; m < modes; m++)
      if (use[m] && !(fixed_visits && m == visit)) cols_used.push_back(m);
    std::vector<std::vector<std::array<std::array<double, modes>, 3>>> cols(H);
    std::vector<std::array<double, modes>> per_unit(H);
    for (std::size_t h = 0; h < H; h++) cols[h] = columns(T, draw, h), per_unit[h] = trips_per_unit(h, draw);
    std::vector<double> rho(C, 0.0);
    const auto solve = [&](std::size_t h, std::size_t leave_out) {
      std::array<double, modes> full{};
      for (int round = 0; round < (fixed_visits ? 4 : 1); round++) {
        std::vector<std::vector<double>> rows;
        std::vector<double> b;
        for (std::size_t c = 0; c < C; c++) {
          if (c == leave_out) continue;
          const auto& k = counters[c];
          const auto row = [&](std::size_t which, double factor, double y) {
            std::vector<double> a(cols_used.size());
            for (std::size_t m = 0; m < cols_used.size(); m++) a[m] = factor * cols[h][c][which][cols_used[m]];
            rows.push_back(std::move(a)), b.push_back(y - factor * cols[h][c][which][visit] * (fixed_visits ? full[visit] : 0.0));
          };
          if (k.directed) row(1, 1 + rho[c], k.way[0][h]), row(2, 1 - rho[c], k.way[1][h]);
          else row(0, 1, k.week[h]);
        }
        const auto xs = tyche::regression::nnls(la::Rectangular::from_rows(rows), b).x;
        const double xv = full[visit];
        full = {};
        for (std::size_t m = 0; m < cols_used.size(); m++) full[cols_used[m]] = xs[m];
        if (fixed_visits) {
          double others = 0;
          for (std::size_t m = 0; m < modes; m++)
            if (m != visit) others += full[m] * per_unit[h][m];
          full[visit] = per_unit[h][visit] > 0 ? o.visit_share / (1 - o.visit_share) * others / per_unit[h][visit] : 0.0;
          (void)xv;
        }
      }
      return full;
    };
    const auto predict = [&](std::size_t h, std::size_t c, std::size_t which, const std::array<double, modes>& xs) {
      double v = 0;
      for (std::size_t m = 0; m < modes; m++) v += cols[h][c][which][m] * xs[m];
      return v;
    };
    std::vector<std::array<double, modes>> xh(H);
    for (int round = 0; round < (o.handedness ? 3 : 1); round++) {
      for (std::size_t h = 0; h < H; h++) xh[h] = solve(h, none);
      if (!o.handedness) break;
      for (std::size_t c = 0; c < C; c++) {
        const auto& k = counters[c];
        if (!k.directed || net.edges[k.edge].kind == "laneway") continue;
        double num = 0, den = 0;
        for (std::size_t h = 0; h < H; h++) {
          const double f1 = predict(h, c, 1, xh[h]), f2 = predict(h, c, 2, xh[h]);
          num += f1 * (k.way[0][h] - f1) - f2 * (k.way[1][h] - f2), den += f1 * f1 + f2 * f2;
        }
        rho[c] = den > 0 ? std::clamp(num / den, -0.9, 0.9) : 0.0;
      }
    }
    for (std::size_t h = 0; h < H; h++) xh[h] = solve(h, none);
    double vt = 0, at = 0;
    for (std::size_t h = 0; h < H; h++) {
      for (std::size_t m = 0; m < modes; m++) fit.x[m][h] = xh[h][m], at += xh[h][m] * per_unit[h][m];
      vt += xh[h][visit] * per_unit[h][visit];
      for (std::size_t c = 0; c < C; c++) {
        const auto xl = solve(h, c);
        const double r1 = 1 + rho[c], r2 = 1 - rho[c];
        const double w1 = predict(h, c, 1, xh[h]), w2 = predict(h, c, 2, xh[h]);
        fit.fitted_way[c][0][h] = r1 * w1, fit.fitted_way[c][1][h] = r2 * w2;
        fit.fitted[c][h] = counters[c].directed ? r1 * w1 + r2 * w2 : predict(h, c, 0, xh[h]);
        fit.loo_way[c][0][h] = predict(h, c, 1, xl), fit.loo_way[c][1][h] = predict(h, c, 2, xl);
        fit.loo[c][h] = predict(h, c, 0, xl);
      }
    }
    fit.hand = rho;
    fit.visit_share = at > 0 ? vt / at : 0;
    const auto r2_of = [&](const std::vector<std::array<double, H>>& pred, bool logs, std::vector<double>* per) {
      double ss = 0, st = 0, mean = 0;
      std::size_t count = 0;
      const auto tr = [&](double v) { return logs ? std::log1p(std::max(0.0, v)) : v; };
      for (std::size_t c = 0; c < C; c++) {
        const auto g = counters[c].gross();
        for (std::size_t h = 0; h < H; h++) mean += tr(g[h]), count++;
      }
      mean /= static_cast<double>(count);
      for (std::size_t c = 0; c < C; c++) {
        const auto g = counters[c].gross();
        double cs = 0, ct = 0, cm = 0;
        for (std::size_t h = 0; h < H; h++) cm += g[h] / H;
        for (std::size_t h = 0; h < H; h++) {
          const double y = tr(g[h]), p = tr(pred[c][h]);
          ss += (y - p) * (y - p), st += (y - mean) * (y - mean);
          cs += (g[h] - pred[c][h]) * (g[h] - pred[c][h]), ct += (g[h] - cm) * (g[h] - cm);
        }
        if (per) per->push_back(1 - cs / ct);
      }
      return 1 - ss / st;
    };
    const auto net_r2_of = [&](const std::vector<std::array<std::array<double, H>, 2>>& pred) {
      double ss = 0, st = 0, mean = 0, count = 0;
      for (std::size_t c = 0; c < C; c++)
        if (counters[c].directed)
          for (std::size_t h = 0; h < H; h++) mean += counters[c].way[0][h] - counters[c].way[1][h], count++;
      if (count == 0) return std::nan("");
      mean /= count;
      for (std::size_t c = 0; c < C; c++)
        if (counters[c].directed)
          for (std::size_t h = 0; h < H; h++) {
            const double y = counters[c].way[0][h] - counters[c].way[1][h], p = pred[c][0][h] - pred[c][1][h];
            ss += (y - p) * (y - p), st += (y - mean) * (y - mean);
          }
      return 1 - ss / st;
    };
    fit.r2 = r2_of(fit.fitted, false, &fit.counter_r2);
    fit.loo_r2 = r2_of(fit.loo, false, &fit.counter_loo_r2);
    fit.r2_log = r2_of(fit.fitted, true, nullptr);
    fit.loo_r2_log = r2_of(fit.loo, true, nullptr);
    fit.net_r2 = net_r2_of(fit.fitted_way), fit.loo_net_r2 = net_r2_of(fit.loo_way);
    if (apply) {
      x = fit.x;
      for (std::size_t c = 0; c < C; c++) counters[c].hand = rho[c];
    }
    return fit;
  }

  // ---- the response for a set of links --------------------------------------------------------------------------

  std::vector<std::size_t> fronted_edges() const {
    std::set<std::size_t> s(frontage_edge.begin(), frontage_edge.end());
    return {s.begin(), s.end()};
  }

  Response response(const hv::Bits& made) const {
    const std::size_t F = P.size();
    const Trips T = trips(made);
    const auto edges = fronted_edges();
    const Terms R0 = terms(T, edges);
    std::vector<std::size_t> at(F);
    for (std::size_t f = 0; f < F; f++) at[f] = static_cast<std::size_t>(std::find(edges.begin(), edges.end(), frontage_edge[f]) - edges.begin());
    Response R;
    R.B.assign(F * H, 0.0), R.R.assign(F * F, 0.0);
    const std::vector<double> no_draw(F * H, 0.0);
    for (std::size_t h = 0; h < H; h++) {
      const auto fl = mode_flows(R0, no_draw, h);
      for (std::size_t f = 0; f < F; f++) {
        double g = 0;
        for (std::size_t m = 0; m < modes; m++) g += fl[m][at[f]][0] + fl[m][at[f]][1];
        R.B[f * H + h] = frontage_share[f] * g;
      }
    }
    for (std::size_t f = 0; f < F; f++)
      for (std::size_t g = 0; g < F; g++) R.R[f * F + g] = frontage_share[f] * (R0.frontage[g][at[f]][0] + R0.frontage[g][at[f]][1]);
    R.distances = P.distances(made);
    return R;
  }

  // The passers of a plan: B + c(h) Σ_g w_{k_g} [g open at h] R_fg.
  std::vector<double> passers(const economy::Plan& plan, const Response& R) const {
    const std::size_t F = P.size();
    std::vector<double> out(R.B);
    const auto draw = frontage_draw(plan);
    for (std::size_t h = 0; h < H; h++) {
      const double c = x[visit][h] / W_ref;
      if (c == 0) continue;
      for (std::size_t g = 0; g < F; g++) {
        const double d = draw[g * H + h];
        if (d == 0) continue;
        for (std::size_t f = 0; f < F; f++) out[f * H + h] += c * d * R.R[f * F + g];
      }
    }
    return out;
  }
};

// A cache of responses by the set of links made, as the problem caches its distances.
class Responses {
 public:
  explicit Responses(const Model& m) : m_(m) {}
  std::shared_ptr<const Response> operator()(const hv::Bits& made) const {
    std::uint64_t key = 0;
    for (std::size_t k = 0; k < made.size(); k++)
      if (made[k]) key |= std::uint64_t{1} << k;
    {
      std::lock_guard<std::mutex> lock(mutex_);
      const auto it = cache_.find(key);
      if (it != cache_.end()) return it->second;
    }
    auto r = std::make_shared<const Response>(m_.response(made));
    std::lock_guard<std::mutex> lock(mutex_);
    return cache_.emplace(key, std::move(r)).first->second;
  }
  std::size_t size() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return cache_.size();
  }

 private:
  const Model& m_;
  mutable std::mutex mutex_;
  mutable std::map<std::uint64_t, std::shared_ptr<const Response>> cache_;
};

// ---- the circulation of the net flow --------------------------------------------------------------------------------------

struct Circulation {
  std::vector<std::size_t> edges;                 // the network's edges present (made), in the complex's order
  std::vector<std::vector<std::pair<std::size_t, int>>> faces;  // the blocks: (index into edges, sign)
  std::vector<double> face_area;
  std::size_t holes = 0;
  economics::FlowCirculation parts;
  std::vector<double> net, gross;
};

// The decomposition of the net flow (u → v minus v → u) of directed flows given on every edge of the network.
inline Circulation circulation(const Network& N, const hv::Bits& made, const std::vector<std::array<double, 2>>& flows) {
  Circulation C;
  std::vector<std::size_t> node(N.nodes(), none);
  std::vector<Point> pts;
  for (std::size_t e = 0; e < N.edges.size(); e++) {
    if (!N.present(e, made)) continue;
    C.edges.push_back(e);
    for (auto j : {N.edges[e].u, N.edges[e].v})
      if (node[j] == none) node[j] = pts.size(), pts.push_back(N.at[j]);
  }
  std::vector<cg::PlaneEdge> pe;
  for (auto e : C.edges) pe.push_back({node[N.edges[e].u], node[N.edges[e].v], N.edges[e].line});
  const auto faces = cg::plane_faces(pts, pe);
  economics::TransferNetwork T;
  for (std::size_t j = 0; j < pts.size(); j++) T.positions.push_back(std::to_string(j));
  economics::Cochain conductance;
  for (std::size_t k = 0; k < C.edges.size(); k++) {
    const auto& x = N.edges[C.edges[k]];
    const double net = flows[C.edges[k]][0] - flows[C.edges[k]][1];
    T.legs.push_back({node[x.u], node[x.v], 1.0, net, "people"});
    conductance.push_back(x.width / x.length);
    C.net.push_back(net), C.gross.push_back(flows[C.edges[k]][0] + flows[C.edges[k]][1]);
  }
  for (const auto& f : faces) {
    if (!(f.area > 0)) continue;
    if (!f.simple) { C.holes++; continue; }
    morphe::topology::Path p;
    std::vector<std::pair<std::size_t, int>> w;
    for (const auto& [e, s] : f.edges) p.push_back({e, s}), w.push_back({e, s});
    T.faces.push_back(std::move(p));
    C.faces.push_back(std::move(w));
    C.face_area.push_back(f.area);
  }
  C.parts = economics::circulation(T, conductance);
  return C;
}

// ---- the week's lines ---------------------------------------------------------------------------------------------------

// X̂(m) = (1/168) Σ_h x_h e^{−2πimh/168}, the convention of economy::week_power.
inline std::complex<double> line(std::span<const double> x, std::size_t m) {
  std::complex<double> s{0, 0};
  for (std::size_t h = 0; h < H; h++) s += x[h] * std::polar(1.0, -2 * std::numbers::pi * static_cast<double>(m * h) / H);
  return s / static_cast<double>(H);
}
// A line's amplitude (2|X̂|, the constituent's) and the hour of its crest within its period: x ≈ A cos(2π m h / 168 + φ),
// crest where the argument is 0 (hours counted from Monday 00:00; the series' values are at the start of each hour).
inline double amplitude(std::complex<double> c) { return 2 * std::abs(c); }
inline double crest_hour(std::complex<double> c, std::size_t m) {
  const double period = static_cast<double>(H) / static_cast<double>(m);
  double t = -std::arg(c) / (2 * std::numbers::pi) * period;
  while (t < 0) t += period;
  while (t >= period) t -= period;
  return t;
}

}  // namespace flow
