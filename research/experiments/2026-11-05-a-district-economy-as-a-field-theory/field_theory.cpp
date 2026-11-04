// A district's economy as a field theory (README.md): the people of the district as a field on its space–time
// complex, calibrated to the counters; the plan's economy evaluated by the shared problem on the field's passers; the
// action's modes, fluctuations and free energies; and searches that use the action's structure against black-box
// searches at equal budgets.
//
//   field_theory <problem-inputs.json> <problem.json> <sensors.json> <counters.json> <results directory> [stage]
//
// From the repository's root. Stages: calibrate (the field alone), all (the default).

#include <chrono>
#include <cstdio>
#include <thread>
#include <iostream>
#include <numeric>

#include <calx/scheduler/thread_pool.hpp>
#include <heuresis/annealing/annealing.hpp>
#include <heuresis/distributions/distributions.hpp>
#include <heuresis/genetic/genetic.hpp>
#include <heuresis/search/search.hpp>

#include <calx/scheduler/parallel_for.hpp>

#include "analysis.hpp"
#include "context.hpp"
#include "plan.hpp"
#include "searches.hpp"

namespace {

using namespace field;
using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }
double r(double x, double d) { return std::isfinite(x) ? std::round(x * d) / d : 0.0; }
Json num(double x, double d = 1e4) { return std::isfinite(x) ? Json(r(x, d)) : Json(); }

Json skill_json(const Skill& S, const std::vector<Counter>& cs, bool series) {
  Json per = Json::array();
  for (std::size_t c = 0; c < cs.size(); c++) {
    Json o = Json::object({{"id", static_cast<double>(cs[c].id)}, {"name", cs[c].name}, {"r2", num(S.r2[c])}, {"total_ratio", num(S.ratio[c])}});
    if (series) o.set("predicted", economy::numbers(S.predicted[c], 1));
    per.push(o);
  }
  return Json::object({{"normalised_sse", num(S.nsse, 1e5)}, {"median_r2", num(S.median_r2)}, {"log_total_rms", num(S.log_total_rms)}, {"log_total_correlation", num(S.total_corr)}, {"counters", per}});
}

// The calibration (README §3): the grid of ξ, c_lane, τ² and the edge group, each scored by leave-one-site-out.
struct Calibration {
  Fitted best;
  Json grid = Json::array();
  Skill nearest;  // the baseline: each counter's mean week taken from its nearest other counter of another site
  Groups groups;
  std::vector<Observation> obs;
  std::shared_ptr<Kriging> kriging;  // the final fit's, on the base graph as calibrated
};

Calibration calibrate(Model& M, const Context& C) {
  Calibration out;
  const auto& cs = C.inside;
  double best = 1e300;
  const std::vector<double> lengths{30, 45, 65, 95, 140, 200, 290, 420};
  const std::vector<std::pair<double, double>> ways{{1.0, 1.0}, {1.0 / 3, 1.0 / 3}, {1.0 / 3, 1.0 / 6}, {1.0 / 3, 1.0 / 10}};
  const std::vector<double> nuggets{0.1, 0.3, 1.0};
  for (const auto& [little, lane] : ways) {
    const Graph G = graph(M.problem->base, little, lane);
    std::vector<Observation> obs;
    for (const auto& c : cs) obs.push_back(observe(G, counter_foot(G, c)));
    std::vector<std::vector<Complex>> y;
    for (const auto& c : cs) y.push_back(the_week().forward(c.week));
    for (int variant = 0; variant < 4; variant++) {
      const bool edge = variant % 2 == 1, trams = variant >= 2;
      const auto gr = source_groups(*M.problem, G, C.inputs, C.outside, edge, trams);
      for (double xi : lengths)
        for (double tau : nuggets)
          for (int rel = 0; rel < 2; rel++) {
          Physics ph = M.phys;
          ph.length = xi, ph.little = little, ph.lane = lane, ph.nugget = tau, ph.relative = rel == 1;
          const Kriging K(G, ph, obs, y, gr);
          const auto [loo, reg] = leave_out(K, cs);
          out.grid.push(Json::object({{"length_m", xi}, {"little", r(little, 1e4)}, {"lane", r(lane, 1e4)}, {"nugget", tau}, {"relative", rel == 1}, {"edge_from_counters", edge}, {"trams", trams},
                                      {"loo_nsse", num(loo.nsse, 1e5)}, {"loo_median_r2", num(loo.median_r2)}, {"regression_nsse", num(reg.nsse, 1e5)},
                                      {"loo_log_total_rms", num(loo.log_total_rms)}, {"loo_log_total_correlation", num(loo.total_corr)}}));
          if (loo.nsse < best) {
            best = loo.nsse;
            out.best.phys = ph, out.best.edge_from_counters = edge, out.best.trams = trams, out.best.loo = loo, out.best.regression_loo = reg;
            out.groups = gr, out.obs = obs;
          }
        }
    }
  }
  // The final fit, on every counter, at the chosen parameters.
  M.phys = out.best.phys;
  M.base = graph(M.problem->base, M.phys.little, M.phys.lane);
  std::vector<std::vector<Complex>> y;
  for (const auto& c : cs) y.push_back(the_week().forward(c.week));
  out.kriging = std::make_shared<Kriging>(M.base, M.phys, out.obs, y, out.groups);
  const Kriging& K = *out.kriging;
  std::vector<std::size_t> all(cs.size());
  std::iota(all.begin(), all.end(), std::size_t{0});
  M.source.assign(lines, {});
  M.sigma2.assign(lines, 0.0);
  std::vector<std::vector<Complex>> pin(cs.size(), std::vector<Complex>(lines));
  for (std::size_t m = 0; m < lines; m++) {
    const auto F = K.fit(m, all);
    out.best.fits.push_back(F);
    M.source[m] = K.source(m, F, out.groups);
    M.sigma2[m] = F.sigma2;
    for (std::size_t q = 0; q < cs.size(); q++) pin[q][m] = K.predict(m, F, q);
  }
  std::vector<std::vector<double>> a;
  for (const auto& v : pin) a.push_back(the_week().inverse(v));
  out.best.in_sample = skill_of(cs, a);
  // The baseline: the nearest counter of another site.
  const auto sites = sites_of(cs);
  std::vector<int> site(cs.size());
  for (std::size_t s = 0; s < sites.size(); s++)
    for (auto q : sites[s]) site[q] = static_cast<int>(s);
  std::vector<std::vector<double>> nb;
  for (std::size_t q = 0; q < cs.size(); q++) {
    std::size_t best_j = q;
    double bd = 1e300;
    for (std::size_t j = 0; j < cs.size(); j++) {
      if (site[j] == site[q]) continue;
      const double d = std::hypot(cs[j].at[0] - cs[q].at[0], cs[j].at[1] - cs[q].at[1]);
      if (d < bd) bd = d, best_j = j;
    }
    nb.push_back(cs[best_j].week);
  }
  out.nearest = skill_of(cs, nb);
  return out;
}

// ---- the people's field: its modes, its conservation law, its currents ------------------------------------------------

Json series(const std::vector<double>& v, int digits = 2) { return economy::numbers(v, digits); }

Point compass(const std::string& name) {
  if (name == "North") return {0, 1};
  if (name == "South") return {0, -1};
  if (name == "East") return {1, 0};
  if (name == "West") return {-1, 0};
  return {0, 0};
}

Json people_field(const Model& M, const Flows& L, const Context& C, const Calibration& cal, Json& geom) {
  const Graph& G = L.G;
  const Problem& P = *M.problem;
  const std::size_t V = G.V, E = G.E, n = P.size();
  const double D = M.phys.diffusivity(), mu = M.phys.leave, v0 = M.phys.speed;
  Json out = Json::object();
  // Geometry, for drawing (ODbL extracts; field-geometry.json is not committed).
  {
    Json vx = Json::array(), ed = Json::array(), ff = Json::array(), cf = Json::array(), lk = Json::array();
    for (const auto& p : G.pos) vx.push(Json::array({r(p[0], 100), r(p[1], 100)}));
    for (std::size_t e = 0; e < E; e++) {
      Json line = Json::array();
      for (const auto& q : G.line[e]) line.push(Json::array({r(q[0], 100), r(q[1], 100)}));
      ed.push(Json::object({{"ends", Json::array({static_cast<double>(G.ends[e][0]), static_cast<double>(G.ends[e][1])})}, {"class", static_cast<double>(G.cls[e])},
                            {"name", G.name[e]}, {"type", G.type[e]}, {"line", line}}));
    }
    for (const auto& o : L.front) ff.push(Json::array({r(o.on[0], 100), r(o.on[1], 100)}));
    for (const auto& c : C.inside) {
      const auto foot = counter_foot(G, c);
      cf.push(Json::object({{"id", static_cast<double>(c.id)}, {"at", Json::array({r(c.at[0], 100), r(c.at[1], 100)})}, {"foot", Json::array({r(foot[0], 100), r(foot[1], 100)})}}));
    }
    for (const auto& l : P.links) lk.push(Json::object({{"a", Json::array({r(l.a[0], 100), r(l.a[1], 100)})}, {"b", Json::array({r(l.b[0], 100), r(l.b[1], 100)})}, {"kind", l.kind}}));
    geom.set("vertices", vx), geom.set("edges", ed), geom.set("frontage_feet", ff), geom.set("counters", cf), geom.set("links", lk);
    geom.set("credit", "(c) OpenStreetMap contributors, ODbL 1.0");
  }
  // The field at the vertices: its weekly mean and five hours.
  std::vector<double> mean(V, 0.0);
  for (std::size_t v = 0; v < V; v++)
    for (std::size_t h = 0; h < week; h++) mean[v] += L.phi[v * week + h] / week;
  Json at = Json::object();
  for (const auto& [name, h] : std::vector<std::pair<std::string, std::size_t>>{{"Monday 03:00", 3}, {"Wednesday 08:00", 56}, {"Wednesday 13:00", 61}, {"Friday 21:00", 117}, {"Saturday 13:00", 133}}) {
    std::vector<double> x(V);
    for (std::size_t v = 0; v < V; v++) x[v] = L.phi[v * week + h];
    at.set(name, series(x, 1));
  }
  out.set("vertex_weekly_mean", series(mean, 2)), out.set("vertex_hours", at), out.set("vertex_star_m", series(G.star0, 2));
  Json fp = Json::array();
  for (std::size_t f = 0; f < n; f++) fp.push(series(std::vector<double>(L.P0.begin() + f * week, L.P0.begin() + (f + 1) * week), 1));
  out.set("frontage_passers", fp);
  // The modes of Δ₀: their rates D λ_k + μ, and the first ten.
  std::vector<double> rate(V);
  for (std::size_t k = 0; k < V; k++) rate[k] = D * G.lambda[k];
  out.set("mode_rates_per_hour", series(rate, 6));
  Json modes = Json::array();
  for (std::size_t k = 0; k < std::min<std::size_t>(V, 24); k++) {
    std::vector<double> u(V);
    for (std::size_t v = 0; v < V; v++) u[v] = G.U(v, k);
    modes.push(series(u, 6));
  }
  out.set("modes", modes);
  // The field's power by mode and line: c_k(m) = Σ_v u_k(v) ⋆₀(v) φ̂_v(m).
  // The kernel first (the numbers of people in each of the network's components), then octaves of the modes above it.
  const std::size_t b0 = static_cast<std::size_t>(std::count_if(G.lambda.begin(), G.lambda.end(), [&](double x) { return x < 1e-9 * G.lambda.back(); }));
  std::vector<std::size_t> edges{0, b0};
  for (std::size_t w = 2; edges.back() + w < V; w *= 2) edges.push_back(edges.back() + w);
  edges.push_back(V);
  Json bins = Json::array();
  std::vector<std::vector<double>> power(edges.size() - 1, std::vector<double>(lines, 0.0));
  std::vector<double> by_mode(V, 0.0), by_line(lines, 0.0);
  for (std::size_t m = 0; m < lines; m++) {
    const double w = (m == 0 || m == lines - 1) ? 1.0 : 2.0;
    for (std::size_t k = 0; k < V; k++) {
      Complex c{0, 0};
      for (std::size_t v = 0; v < V; v++) c += G.U(v, k) * G.star0[v] * L.lines_at_vertices[m][v];
      const double p = w * std::norm(c);
      by_mode[k] += p, by_line[m] += p;
      for (std::size_t b = 0; b + 1 < edges.size(); b++)
        if (k >= edges[b] && k < edges[b + 1]) power[b][m] += p;
    }
  }
  for (std::size_t b = 0; b + 1 < edges.size(); b++)
    bins.push(Json::object({{"modes", Json::array({static_cast<double>(edges[b]), static_cast<double>(std::min(edges[b + 1], V) - 1)})},
                            {"rate_per_hour", Json::array({r(rate[edges[b]], 1e4), r(rate[std::min(edges[b + 1], V) - 1], 1e4)})}, {"power", series(power[b], 3)}}));
  out.set("spectrum", Json::object({{"bins", bins}, {"by_mode", series(by_mode, 3)}, {"by_line", series(by_line, 3)}}));
  // The conservation law (README §2.3): the number of people N = (2/v₀) Σ ⋆₀ φ obeys dN/dt = (2/v₀) Σ s − μ N, at every line.
  std::vector<Complex> Nh(lines), Sh(lines);
  double worst = 0, scale = 0;
  for (std::size_t m = 0; m < lines; m++) {
    Complex a{0, 0}, s{0, 0};
    for (std::size_t v = 0; v < V; v++) a += G.star0[v] * L.lines_at_vertices[m][v];
    for (std::size_t v = 0; v < M.base.V; v++) s += M.source[m][v];
    Nh[m] = a * 2.0 / v0, Sh[m] = s * 2.0 / v0;
    worst = std::max(worst, std::abs(Complex(mu, omega(m)) * Nh[m] - Sh[m]));
    scale = std::max(scale, std::abs(Sh[m]));
  }
  const auto Nt = the_week().inverse(Nh), St = the_week().inverse(Sh);
  out.set("noether", Json::object({{"people_on_the_streets", series(Nt, 1)}, {"inflow_per_hour", series(St, 1)}, {"largest_relative_residual", worst / scale},
                                   {"statement", "dN/dt = (2/v0) sum_v s_v - mu N, N = (2/v0) sum_v star0(v) phi_v: the k = 0 mode of the stationary path, at every line m"}}));
  // Fick's current on the edges, people an hour on both footpaths, tail to head.
  const auto current = [&](std::size_t e, std::size_t h) {
    return -D * G.conductance[e] * (L.phi[G.ends[e][1] * week + h] - L.phi[G.ends[e][0] * week + h]) * 2.0 / v0;
  };
  std::vector<double> jmean(E, 0.0), jsat(E, 0.0);
  for (std::size_t e = 0; e < E; e++) {
    for (std::size_t h = 0; h < week; h++) jmean[e] += std::fabs(current(e, h)) / week;
    jsat[e] = current(e, 133);
  }
  out.set("current_weekly_mean_abs", series(jmean, 2)), out.set("current_saturday_13", series(jsat, 2));
  // Against the counters' directions: the net flow on the counter's footpath along the segment, tail to head.
  Json dirs = Json::array();
  const auto& D2 = C.counters["directions"];
  double sum_corr = 0, cnt = 0;
  for (const auto& c : C.inside) {
    const std::string key = std::to_string(c.id);
    if (!D2.has(key)) continue;
    const auto& d = D2[key];
    const auto o = observe(G, counter_foot(G, c));
    const Point t{G.pos[o.b][0] - G.pos[o.a][0], G.pos[o.b][1] - G.pos[o.a][1]};
    const Point k = compass(d["names"][0].text());
    const double sign = (t[0] * k[0] + t[1] * k[1]) >= 0 ? 1.0 : -1.0;
    std::vector<double> obs(week), mod(week);
    double so = 0, sm = 0, sxy = 0, sxx = 0, syy = 0, mo = 0, mm = 0;
    for (std::size_t h = 0; h < week; h++) {
      const auto& a = d["d1"][h];
      const auto& b = d["d2"][h];
      obs[h] = (a.is_null() || b.is_null()) ? 0.0 : sign * (a.as_number() - b.as_number());
      mod[h] = 0.5 * current(o.edge, h);
      mo += obs[h] / week, mm += mod[h] / week;
    }
    for (std::size_t h = 0; h < week; h++) sxy += (obs[h] - mo) * (mod[h] - mm), sxx += (obs[h] - mo) * (obs[h] - mo), syy += (mod[h] - mm) * (mod[h] - mm), so += obs[h] * obs[h], sm += mod[h] * mod[h];
    const double corr = sxx > 0 && syy > 0 ? sxy / std::sqrt(sxx * syy) : 0.0;
    sum_corr += corr, cnt++;
    dirs.push(Json::object({{"id", static_cast<double>(c.id)}, {"name", c.name}, {"along", d["names"][0].text()}, {"street", G.name[o.edge]}, {"observed", series(obs, 1)}, {"model", series(mod, 2)},
                            {"correlation", r(corr, 1e4)}, {"rms_observed", r(std::sqrt(so / week), 10)}, {"rms_model", r(std::sqrt(sm / week), 100)}}));
  }
  out.set("directions", Json::object({{"counters", dirs}, {"mean_correlation", r(cnt > 0 ? sum_corr / cnt : 0.0, 1e4)}}));
  // The Hodge Laplacians of every degree (stars ⋆₀, ⋆₁ per unit D, and 1/area on the blocks): the density's modes,
  // the currents' (gradient, curl round the blocks, harmonic round the holes) and the blocks'.
  {
    std::vector<double> area;
    for (std::size_t i = 0; i < P.base.complex.size(); i++)
      if (P.base.complex.cells[i].dim == 2) area.push_back(std::max(1.0, P.base.area(P.base.complex.cells[i].id)));
    std::vector<double> w2(area.size());
    for (std::size_t i = 0; i < area.size(); i++) w2[i] = 1.0 / area[i];
    const auto H1 = physis::discretisation::hodge_laplacian(G.chain, 1, {G.star0, G.conductance, w2});
    const auto H2 = physis::discretisation::hodge_laplacian(G.chain, 2, {G.star0, G.conductance, w2});
    std::vector<double> l1 = H1.values, l2 = H2.values;
    out.set("hodge", Json::object({{"delta0", series(G.lambda, 9)}, {"delta1", series(l1, 9)}, {"delta2", series(l2, 9)},
                                   {"betti", Json::array({static_cast<double>(std::count_if(G.lambda.begin(), G.lambda.end(), [&](double x) { return x < 1e-9 * G.lambda.back(); })),
                                                          static_cast<double>(H1.harmonic(1e-9)), static_cast<double>(H2.harmonic(1e-9))})},
                                   {"cells", Json::array({static_cast<double>(V), static_cast<double>(E), static_cast<double>(area.size())})}}));
  }
  (void)cal;
  return out;
}

// ---- the fluctuations: the model against the counters' residuals ----------------------------------------------------------

Json fluctuation_json(const Fluctuations& Fl, const Kriging& K, const Context& C) {
  const auto& ids = C.counters["district_counters"].as_array();
  const auto& obs = C.counters["local"]["correlation"].as_array();
  const std::size_t n = ids.size();
  std::vector<std::size_t> row(n);
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t c = 0; c < C.inside.size(); c++)
      if (C.inside[c].id == static_cast<long long>(ids[i].as_number())) row[i] = c;
  Json pairs = Json::array();
  double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0, m = 0;
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t j = i + 1; j < n; j++) {
      if (obs[i][j].is_null()) continue;
      const std::size_t a = row[i], b = row[j];
      const double pred = K.K[0][a][b] / std::sqrt(K.K[0][a][a] * K.K[0][b][b]);
      const double o = obs[i][j].as_number();
      const double d = std::hypot(C.inside[a].at[0] - C.inside[b].at[0], C.inside[a].at[1] - C.inside[b].at[1]);
      pairs.push(Json::array({static_cast<double>(C.inside[a].id), static_cast<double>(C.inside[b].id), r(d, 10), r(pred, 1e4), r(o, 1e4)}));
      sx += pred, sy += o, sxx += pred * pred, syy += o * o, sxy += pred * o, m++;
    }
  const double corr = (sxy - sx * sy / m) / std::sqrt((sxx - sx * sx / m) * (syy - sy * sy / m));
  std::vector<double> sd(week);
  for (std::size_t h = 0; h < week; h++) sd[h] = std::sqrt(std::max(0.0, Fl.common[h * week + h]));
  return Json::object({{"T_zeta", r(Fl.T_zeta, 1e6)}, {"relative_variance_by_counter", series(Fl.relative_variance, 5)}, {"local_autocorrelation", series(Fl.rho, 4)},
                       {"common_sd_by_hour", series(sd, 4)}, {"pairs", pairs}, {"pairs_are", "counter a, counter b, distance m, predicted correlation (the action's, m = 0), observed correlation of the local residual"},
                       {"correlation_agreement", r(corr, 1e4)}});
}

// ---- the links: the linear response of the field (FDT) against the network rebuilt -----------------------------------

double normal_cdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

Json link_response(const Model& M, const Flows& L0, const Plan& plan0, const Plan& best) {
  const Problem& P = *M.problem;
  const Graph& G = L0.G;
  const std::size_t n = P.size(), V = G.V;
  const double D = M.phys.diffusivity();
  const PlanState S0(M, plan0);
  const double V0 = M.evaluate(plan0).value, Vb = M.evaluate(best).value;
  const PlanState Sb(M, best);
  Json out = Json::array();
  std::vector<std::vector<Complex>> t(lines);
  for (std::size_t m = 0; m < lines; m++) t[m] = transfer(G, M.phys, m);
  for (std::size_t j = 0; j < P.links.size(); j++) {
    economy::hv::Bits made(P.links.size(), 0);
    made[j] = 1;
    const auto Lj = M.flows(made);                 // the lattice with the link: Woodbury
    const Flows Rj = M.build(made, false);          // the domain's complex rebuilt with the link, subdivided alike
    std::vector<double> exact(n * week), rebuilt(n * week);
    for (std::size_t i = 0; i < n * week; i++) exact[i] = Lj->P0[i] - L0.P0[i], rebuilt[i] = Rj.P0[i] - L0.P0[i];
    // The link as one new edge between its ends' Whitney points, of a lane's conductance D c_lane / ℓ.
    const auto& lk = P.links[j];
    const auto oa = observe(G, lk.a), ob = observe(G, lk.b);
    std::vector<Complex> b(V, Complex{0, 0});
    b[oa.a] += oa.wa, b[oa.b] += oa.wb, b[ob.a] -= ob.wa, b[ob.b] -= ob.wb;
    const double kappa = D * way_conductance(lane_way, M.phys.little, M.phys.lane) / std::max(1.0, lk.length);
    std::vector<std::vector<Complex>> first(n, std::vector<Complex>(lines)), resummed(n, std::vector<Complex>(lines));
    for (std::size_t m = 0; m < lines; m++) {
      const auto Gb = respond(G, t[m], b);
      Complex bphi{0, 0}, bGb{0, 0};
      for (std::size_t v = 0; v < V; v++) bphi += b[v] * L0.lines_at_vertices[m][v], bGb += b[v] * Gb[v];
      for (std::size_t f = 0; f < n; f++) {
        const auto& o = L0.front[f];
        const Complex g = o.wa * Gb[o.a] + o.wb * Gb[o.b];
        first[f][m] = -kappa * g * bphi;
        resummed[f][m] = first[f][m] / (1.0 + kappa * bGb);
      }
    }
    std::vector<double> lin(n * week), sm(n * week);
    for (std::size_t f = 0; f < n; f++) {
      const auto a = the_week().inverse(first[f]), c = the_week().inverse(resummed[f]);
      for (std::size_t h = 0; h < week; h++) lin[f * week + h] = a[h], sm[f * week + h] = c[h];
    }
    const auto weekly = [&](const std::vector<double>& x) {
      std::vector<double> w(n, 0.0);
      for (std::size_t f = 0; f < n; f++)
        for (std::size_t h = 0; h < week; h++) w[f] += x[f * week + h];
      return w;
    };
    const auto we = weekly(exact), wl = weekly(lin), ws = weekly(sm), wr = weekly(rebuilt);
    const auto corr = [&](const std::vector<double>& x, const std::vector<double>& y) {
      double mx = 0, my = 0;
      for (std::size_t i = 0; i < x.size(); i++) mx += x[i] / static_cast<double>(x.size()), my += y[i] / static_cast<double>(x.size());
      double a = 0, bb = 0, c = 0;
      for (std::size_t i = 0; i < x.size(); i++) a += (x[i] - mx) * (y[i] - my), bb += (x[i] - mx) * (x[i] - mx), c += (y[i] - my) * (y[i] - my);
      return bb > 0 && c > 0 ? a / std::sqrt(bb * c) : 0.0;
    };
    const auto dv = [&](const PlanState& S, const std::vector<double>& x) {
      double s = 0;
      for (std::size_t i = 0; i < n * week; i++) s += S.A[i] * x[i];
      return s;
    };
    double se = 0, sl = 0, ss = 0, sr = 0;
    for (std::size_t f = 0; f < n; f++) se += we[f], sl += wl[f], ss += ws[f], sr += wr[f];
    Plan with = plan0;
    with.links[j] = 1;
    Plan bt = best;
    bt.links[j] ^= 1;
    out.push(Json::object({{"link", static_cast<double>(j)}, {"joins", Json::array({lk.from, lk.to})}, {"kind", lk.kind}, {"length_m", r(lk.length, 10)}, {"conductance", r(kappa, 1e3)},
                           {"weekly_passers_change", Json::object({{"exact", r(se, 1)}, {"linear_response", r(sl, 1)}, {"resummed", r(ss, 1)}, {"domain_rebuilt", r(sr, 1)}})},
                           {"frontage_correlation", Json::object({{"linear_response", r(corr(we, wl), 1e4)}, {"resummed", r(corr(we, ws), 1e4)}, {"domain_rebuilt", r(corr(we, wr), 1e4)}})},
                           {"frontage_weekly_change", Json::object({{"exact", series(we, 1)}, {"linear_response", series(wl, 1)}, {"resummed", series(ws, 1)}, {"domain_rebuilt", series(wr, 1)}})},
                           {"value_change_as_built", Json::object({{"flow_exact", r(dv(S0, exact), 1)}, {"flow_linear_response", r(dv(S0, lin), 1)}, {"flow_resummed", r(dv(S0, sm), 1)},
                                                                    {"flow_domain_rebuilt", r(dv(S0, rebuilt), 1)}, {"all_exact", r(M.evaluate(with).value - V0, 1)}, {"cost", r(P.P.link_weekly(lk.length), 1)}})},
                           {"value_change_at_best", Json::object({{"toggled_exact", r(M.evaluate(bt).value - Vb, 1)}, {"made_in_best", static_cast<double>(best.links[j])},
                                                                  {"flow_linear_response", r((best.links[j] ? -1.0 : 1.0) * dv(Sb, lin), 1)}})}}));
  }
  return out;
}

// ---- the plans: their values, uncertainties and free energies ------------------------------------------------------------

Json plans_json(const Model& M, const Posterior& Q, const Fluctuations& Fl, const std::vector<std::pair<std::string, Plan>>& named) {
  const Problem& P = *M.problem;
  const std::size_t n = P.size();
  Json out = Json::array();
  std::vector<PlanState> states;
  for (const auto& [nm, pl] : named) states.emplace_back(M, pl);
  const auto& S0 = states[0];
  const std::vector<double> temps{1e4, 3e4, 1e5};
  for (std::size_t i = 0; i < named.size(); i++) {
    const auto& S = states[i];
    const auto o = M.evaluate(named[i].second);
    std::vector<double> q(n * week), dq(n * week), da(n * week);
    for (std::size_t k = 0; k < n * week; k++) q[k] = S.A[k] * S.passers(k), dq[k] = q[k] - S0.A[k] * S0.passers(k), da[k] = S.A[k] - S0.A[k];
    const double vpost = posterior_variance(Q, S.A);
    const auto wv = week_variance(Fl, Q, q);
    const double dpost = posterior_variance(Q, da);
    const auto dwv = week_variance(Fl, Q, dq);
    Json fe = Json::array();
    for (double T : temps)
      fe.push(Json::object({{"T", T}, {"free_energy_week", r(o.value - wv.total() / (2 * T), 1)}, {"free_energy_week_and_posterior", r(o.value - (wv.total() + vpost) / (2 * T), 1)},
                            {"difference_from_as_built", r((o.value - M.evaluate(named[0].second).value) - (dwv.total() + dpost) / (2 * T), 1)}}));
    const double dE = o.value - M.evaluate(named[0].second).value;
    Json pass = Json::array();
    for (std::size_t f = 0; f < n; f++) {
      std::vector<double> x(week);
      for (std::size_t h = 0; h < week; h++) x[h] = S.passers(f * week + h);
      pass.push(series(x, 1));
    }
    out.push(Json::object({{"name", named[i].first}, {"plan", economy::plan_json(named[i].second)}, {"outcome", economy::outcome_json(o, true)},
                           {"posterior_sd", r(std::sqrt(std::max(0.0, vpost)), 1)}, {"week_sd", r(std::sqrt(std::max(0.0, wv.total())), 1)},
                           {"week_sd_common", r(std::sqrt(std::max(0.0, wv.common)), 1)}, {"week_sd_local", r(std::sqrt(std::max(0.0, wv.local)), 1)},
                           {"difference_from_as_built", Json::object({{"expected", r(dE, 1)}, {"posterior_sd", r(std::sqrt(std::max(0.0, dpost)), 1)},
                                                                      {"week_sd", r(std::sqrt(std::max(0.0, dwv.total())), 1)},
                                                                      {"probability_better", i == 0 ? Json() : Json(r(normal_cdf(dE / std::sqrt(std::max(1e-12, dpost))), 1e5))}})},
                           {"free_energies", fe}, {"passers", pass}}));
  }
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 6) {
    std::cerr << "usage: field_theory <problem-inputs.json> <problem.json> <sensors.json> <counters.json> <results directory> [stage]\n";
    return 2;
  }
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  const auto T0 = Clock::now();
  const std::string outdir = argv[5];
  const std::string stage = argc > 6 ? argv[6] : "all";
  Context C = load(argv[1], argv[2], argv[3], argv[4]);
  const Problem& P = *C.problem;
  Physics phys;
  Model M(P, phys);
  std::printf("district: %zu frontages, %zu links; base graph %zu vertices, %zu edges; %zu counters on the network, %zu around it (%.1f s)\n", P.size(), P.links.size(), M.base.V,
              M.base.E, C.inside.size(), C.outside.size(), seconds_since(T0));

  // ---- 1. calibration ----
  auto T1 = Clock::now();
  const auto cal = calibrate(M, C);
  const auto& B = cal.best;
  std::printf("calibration: xi = %.0f m, little %.2f, lane %.2f, nugget %.2f%s, trams %d, edges %s; LOO nsse %.4f, median R2 %.3f, log-total rms %.3f (regression alone nsse %.4f; nearest counter nsse %.4f, median R2 %.3f); in-sample median R2 %.3f (%.1f s)\n",
              B.phys.length, B.phys.little, B.phys.lane, B.phys.nugget, B.phys.relative ? " (relative)" : "", B.trams ? 1 : 0, B.edge_from_counters ? "from counters" : "uniform", B.loo.nsse, B.loo.median_r2, B.loo.log_total_rms, B.regression_loo.nsse,
              cal.nearest.nsse, cal.nearest.median_r2, B.in_sample.median_r2, seconds_since(T1));
  M.p0 = M.chances(P.as_built(), *P.distances(economy::hv::Bits(P.links.size(), 0)));
  M.prepare();
  Json calj = Json::object({{"chosen", Json::object({{"length_m", B.phys.length}, {"little", r(B.phys.little, 1e4)}, {"lane", r(B.phys.lane, 1e4)}, {"nugget", B.phys.nugget},
                                                     {"relative_nugget", B.phys.relative}, {"edge_from_counters", B.edge_from_counters}, {"trams", B.trams},
                                                     {"groups", [&] { Json g = Json::array(); for (const auto& x : cal.groups.names) g.push(x); return g; }()},
                                                     {"leave_per_hour", B.phys.leave}, {"diffusivity_m2_per_h", r(B.phys.diffusivity(), 1)}, {"attraction", B.phys.attraction},
                                                     {"speed_m_per_h", B.phys.speed}})},
                            {"loo", skill_json(B.loo, C.inside, true)}, {"regression_loo", skill_json(B.regression_loo, C.inside, false)},
                            {"in_sample", skill_json(B.in_sample, C.inside, true)}, {"nearest_counter", skill_json(cal.nearest, C.inside, true)}, {"grid", cal.grid},
                            {"sigma2_by_line", economy::numbers(M.sigma2, 3)}});
  write_file(outdir + "/calibration.json", calj.dump());
  if (stage == "calibrate") return 0;

  // ---- 2. the people's field ----
  T1 = Clock::now();
  const economy::hv::Bits none(P.links.size(), 0);
  const Flows Lk = M.build(none, true);
  const Graph& G = Lk.G;
  Json geom = Json::object();
  const Json fieldj = people_field(M, Lk, C, cal, geom);
  write_file(outdir + "/field.json", fieldj.dump());
  std::printf("the people's field: %zu modes, Noether residual %.2e; (%.1f s)\n", G.V, fieldj["noether"]["largest_relative_residual"].as_number(), seconds_since(T1));
  if (stage == "field") return 0;
  if (stage == "links") {
    // The links' comparison alone, at the best plan of a finished run (results/search.json), written back into it.
    Json S = Json::parse(read_file(outdir + "/search.json"));
    Plan best;
    for (const auto& p : S["plans"].as_array())
      if (p["name"].text() == "best heat bath") best = economy::plan_of(p["plan"]);
    const Json linksj = link_response(M, Lk, P.as_built(), best);
    S.set("links", linksj);
    write_file(outdir + "/search.json", S.dump());
    std::printf("links rewritten (%.1f s)\n", seconds_since(T0));
    return 0;
  }

  // ---- 3. the posterior and the fluctuations ----
  T1 = Clock::now();
  const Posterior Q = posterior(*cal.kriging, {}, Lk.front, M.sigma2);
  const Fluctuations Fl = fluctuations(C.counters, C.inside);
  Json flj = fluctuation_json(Fl, *cal.kriging, C);
  {
    // Each frontage's weekly passers as built and their posterior standard deviation given the counters (the line m = 0).
    std::vector<double> mean_total(P.size(), 0.0), sd_total(P.size(), 0.0);
    for (std::size_t f = 0; f < P.size(); f++) {
      for (std::size_t h = 0; h < week; h++) mean_total[f] += Lk.P0[f * week + h];
      sd_total[f] = static_cast<double>(week) * std::sqrt(std::max(0.0, Q.C[0][f * P.size() + f].real()));
    }
    flj.set("frontage_weekly_total", economy::numbers(mean_total, 1));
    flj.set("frontage_weekly_total_posterior_sd", economy::numbers(sd_total, 1));
  }
  std::printf("posterior and fluctuations: T_zeta %.4f (relative variance of the local field); spatial correlation of residuals, predicted against observed r = %.3f (%.1f s)\n",
              Fl.T_zeta, flj["correlation_agreement"].as_number(), seconds_since(T1));

  // ---- 4. the plan as built and the searches ----
  const Plan plan0 = P.as_built();
  {
    PlanState S0(M, plan0);
    const auto o = M.evaluate(plan0);
    std::printf("as built: V = %.1f A$/week (the plan field's own sum %.1f); footfall %.0f, spread %.3f, mix %.3f, access %.1f m\n", o.value, S0.value, o.footfall, o.spread, o.mix, o.access);
  }
  // Warm the caches of every link set (the problem's distances and the field's passers), on four threads.
  T1 = Clock::now();
  calx::scheduler::ThreadPool pool(4);
  {
    const std::size_t sets = stage == "quick" ? 1 : std::size_t{1} << P.links.size();
    calx::scheduler::parallel_for(pool, sets, 8, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t k = b; k < e; k++) {
        economy::hv::Bits made(P.links.size(), 0);
        for (std::size_t j = 0; j < P.links.size(); j++) made[j] = (k >> j) & 1;
        P.distances(made);
        M.flows(made);
      }
    });
  }
  std::printf("warmed %zu link sets in %.1f s\n", M.cached(), seconds_since(T1));
  economy::Kernels K{&P};
  double eval_seconds = 0;
  {
    calx::random::Stream s(4242, 1, 0, 0);
    std::vector<Plan> sample;
    for (int i = 0; i < 400; i++) sample.push_back(K.sample(s));
    const auto t = Clock::now();
    double sink = 0;
    for (const auto& p : sample) sink += M.evaluate(p).value;
    eval_seconds = seconds_since(t) / static_cast<double>(sample.size());
    std::printf("one evaluation (passers and the problem's evaluate): %.3f ms (mean value of random plans %.0f)\n", 1000 * eval_seconds, sink / static_cast<double>(sample.size()));
  }
  // The stage "quick" (a test of the pipeline end to end) runs one seed at a small budget.
  const double budget = stage == "quick" ? 300 : 20000;
  const std::vector<std::uint64_t> seeds = stage == "quick" ? std::vector<std::uint64_t>{1} : std::vector<std::uint64_t>{1, 2, 3, 4, 5};
  std::vector<Run> runs;
  // Separable CMA-ES spends most of its time on its own updates in 621 dimensions: its seeds run side by side, each on
  // one thread.
  std::vector<Run> cma_runs(seeds.size());
  {
    const auto t = Clock::now();
    std::vector<std::thread> th;
    for (std::size_t i = 0; i < seeds.size(); i++) th.emplace_back([&, i] { cma_runs[i] = run_cma(M, seeds[i], budget, nullptr); });
    for (auto& x : th) x.join();
    std::printf("CMA-ES, %zu seeds side by side: %.1f s\n", seeds.size(), seconds_since(t));
  }
  for (std::size_t si = 0; si < seeds.size(); si++) {
    const auto seed = seeds[si];
    T1 = Clock::now();
    runs.push_back(run_ga(M, K, seed, budget, pool, false));
    runs.push_back(run_ga(M, K, seed, budget, pool, true));
    runs.push_back(cma_runs[si]);
    runs.push_back(run_sa(M, K, seed, budget, 2000));
    runs.push_back(run_umda(M, seed, budget, pool));
    runs.push_back(run_random(M, K, seed, budget, pool));
    HeatBathOptions hb;
    runs.push_back(run_heat_bath(M, plan0, seed, hb, eval_seconds, "heat bath (from as built)"));
    calx::random::Stream s(seed, 0x7374, 0x6172, 0x74);
    runs.push_back(run_heat_bath(M, K.sample(s), seed, hb, eval_seconds, "heat bath (from a random plan)"));
    std::printf("seed %llu:", static_cast<unsigned long long>(seed));
    for (std::size_t i = runs.size() - 8; i < runs.size(); i++) std::printf("  %s %.0f (%.0f ev, %.1f s)", runs[i].algorithm.c_str(), runs[i].value, runs[i].evaluations, runs[i].seconds);
    std::printf("  [%.1f s]\n", seconds_since(T1));
  }
  // The best of each family, checked against the problem's own evaluate().
  const auto best_of = [&](const std::string& prefix) {
    const Run* b = nullptr;
    for (const auto& x : runs)
      if (x.algorithm.rfind(prefix, 0) == 0 && (!b || x.value > b->value)) b = &x;
    return b;
  };
  const Run* best_bb = nullptr;
  for (const auto& x : runs)
    if (x.algorithm.rfind("heat bath", 0) != 0 && (!best_bb || x.value > best_bb->value)) best_bb = &x;
  const Run* best_hb = best_of("heat bath");
  {
    const auto a = M.evaluate(best_hb->plan), b = M.evaluate(best_bb->plan);
    std::printf("best heat bath %.1f (evaluate %.1f); best black box %s %.1f (evaluate %.1f)\n", best_hb->value, a.value, best_bb->algorithm.c_str(), best_bb->value, b.value);
  }

  // ---- 5. the links: linear response against the network rebuilt ----
  T1 = Clock::now();
  const Json linksj = link_response(M, Lk, plan0, best_hb->plan);
  std::printf("links: linear response against exact (%.1f s)\n", seconds_since(T1));

  // ---- 6. free energies, and the attraction's sensitivity ----
  T1 = Clock::now();
  std::vector<std::pair<std::string, Plan>> named{{"as built", plan0}, {"best heat bath", best_hb->plan}, {"best black box (" + best_bb->algorithm + ")", best_bb->plan}};
  for (const char* fam : {"GA", "CMA-ES", "SA", "UMDA", "random"})
    if (const Run* x = best_of(fam)) named.push_back({std::string("best ") + fam, x->plan});
  const Json plansj = plans_json(M, Q, Fl, named);
  Json sens = Json::array();
  for (double b : {0.0, 20.0}) {
    Physics ph = M.phys;
    ph.attraction = b;
    Model Mb(P, ph);
    Mb.source = M.source, Mb.sigma2 = M.sigma2, Mb.p0 = M.p0;
    Mb.prepare();
    HeatBathOptions hb;
    const Run x = run_heat_bath(Mb, plan0, 1, hb, eval_seconds, "heat bath");
    Json vals = Json::object();
    for (const auto& [nm, pl] : named) vals.set(nm, r(Mb.evaluate(pl).value, 1));
    sens.push(Json::object({{"attraction", b}, {"heat_bath_from_as_built", r(x.value, 1)}, {"plan", economy::plan_json(x.plan)}, {"values_of_the_plans", vals},
                            {"kinds_changed_from_best", [&] { double d = 0; for (std::size_t f = 0; f < P.size(); f++) d += x.plan.kind[f] != best_hb->plan.kind[f]; return d; }()}}));
  }
  std::printf("free energies and sensitivity (%.1f s)\n", seconds_since(T1));

  Json runsj = Json::array();
  for (const auto& x : runs) {
    Json curve = Json::array();
    for (const auto& c : x.curve) curve.push(Json::array({r(c[0], 10), r(c[1], 1)}));
    runsj.push(Json::object({{"algorithm", x.algorithm}, {"seed", static_cast<double>(x.seed)}, {"value", r(x.value, 1)}, {"evaluations", r(x.evaluations, 1)},
                             {"seconds", r(x.seconds, 100)}, {"curve", curve}, {"plan", economy::plan_json(x.plan)}}));
  }
  write_file(outdir + "/search.json",
             Json::object({{"budget_evaluations", budget}, {"evaluation_ms", r(1000 * eval_seconds, 1e4)}, {"seeds", Json::array({1.0, 2.0, 3.0, 4.0, 5.0})}, {"runs", runsj},
                           {"plans", plansj}, {"links", linksj}, {"attraction_sensitivity", sens}, {"fluctuations", flj}})
                 .dump());
  // Geometry for drawing: extracts of OpenStreetMap (ODbL), not committed.
  write_file(outdir + "/field-geometry.json", geom.dump());
  std::printf("done in %.1f s\n", seconds_since(T0));
  return 0;
}
