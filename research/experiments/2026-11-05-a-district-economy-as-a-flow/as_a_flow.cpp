// A district's economy as a flow: study (2) of the shared problem (README.md; the problem in
// ../2026-11-05-a-district-economy-the-problem). People are a conserved flow on the district's streets (flow.hpp); the
// search uses the model's structure (structured.hpp).
//
//   as_a_flow <problem-inputs.json> <problem.json> <sensors.json> <results directory> [--threads 4] [--seeds 5] [--scale 1]
//
// 1. Calibration. The amplitudes of the trips' patterns, hour by hour, by NNLS at the City's counters inside the box,
//    over a grid of the route choice's θ; leave-one-counter-out prediction measures the fit; the boundary's consistency
//    and two ablations are reported.
// 2. Every set of links. The response of the passers to the plan (B, R) for each of the 2,048 sets of links, by
//    Woodbury's update of one inverse.
// 3. Searches at equal budgets over seeds 1…5, as study (1) runs them (400 and 20,000 evaluations): the genetic algorithm
//    on the mixed encoding, separable CMA-ES on the canonical relaxation, random search (the black box); and three that
//    use the structure: block-coordinate ascent of exact best responses with a descent over the links (BR), the same
//    with a (1+1) EA over the links (BR + EA), and the GA with kernels biased by each frontage's marginal value
//    (GA, marginal kernels). Work is counted in evaluations of the value (168 margins of an hour each).
// 4. The flows of the plans: the people on every segment, their week's lines (as phasors), the circulation of their net
//    flow (driven, round the blocks, round the holes), what each link does, the spending as flux through the frontages,
//    and the arrangement the plans use.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <heuresis/genetic/genetic.hpp>
#include <heuresis/search/search.hpp>
#include <heuresis/strategies/strategies.hpp>

#include "economy.hpp"
#include "flow.hpp"
#include "structured.hpp"

namespace {

using economy::Json;
using economy::numbers;
using economy::Plan;
namespace hs = heuresis::search;
namespace hv = heuresis::variation;
constexpr std::size_t H = economy::week;
using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

std::string read_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read " + path);
  std::stringstream b;
  b << in.rdbuf();
  return b.str();
}
void write_file(const std::string& path, const std::string& text) {
  std::ofstream f(path, std::ios::binary);
  f << text;
  if (!f) throw std::runtime_error("cannot write " + path);
}
double r(double x, double d) { return std::isfinite(x) ? std::round(x * d) / d : std::nan(""); }
Json jr(double x, double d) { return std::isfinite(x) ? Json(r(x, d)) : Json(); }
hv::Bits bits_of(std::uint64_t key, std::size_t L) {
  hv::Bits b(L, 0);
  for (std::size_t k = 0; k < L; k++) b[k] = (key >> k) & 1u;
  return b;
}

// A black-box run's tape: every evaluation's value by ticket, and its best plan.
struct Tape {
  std::mutex m;
  std::map<hs::Ticket, double> values;
  double best = -std::numeric_limits<double>::infinity();
  Plan best_plan;
  hs::Ticket best_ticket = 0;
  void add(hs::Ticket t, double v, const Plan& p) {
    std::lock_guard<std::mutex> lock(m);
    values[t] = v;
    if (v > best || (v == best && t < best_ticket)) best = v, best_plan = p, best_ticket = t;
  }
};

struct Run {
  Run() = default;
  Run(std::string a, std::size_t b, std::uint64_t s) : algorithm(std::move(a)), budget(b), seed(s) {}
  std::string algorithm;
  std::size_t budget = 0;
  std::uint64_t seed = 0;
  double seconds = 0, evaluations = 0;
  double best = -std::numeric_limits<double>::infinity();
  Plan best_plan;
  std::vector<std::pair<double, double>> anytime;  // (evaluations, best so far)
  std::size_t inner = 0;                           // structured: link sets visited
};

// The anytime curve of a tape, offset by work spent before it (in evaluations).
std::vector<std::pair<double, double>> anytime_of(const Tape& t, double offset, std::size_t step) {
  std::vector<std::pair<double, double>> out;
  double best = -std::numeric_limits<double>::infinity();
  std::size_t i = 0;
  for (const auto& [ticket, v] : t.values) {
    (void)ticket;
    best = std::max(best, v), i++;
    if (i % step == 0 || i == t.values.size()) out.push_back({offset + static_cast<double>(i), best});
  }
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 5) {
    std::cerr << "usage: as_a_flow <problem-inputs.json> <problem.json> <sensors.json> <results directory> [--threads N] [--seeds k] [--scale s]\n";
    return 2;
  }
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  std::size_t threads = 4, seeds_n = 5;
  double scale = 1;
  bool calibrate_only = false, extra_only = false;
  for (int i = 5; i < argc; i++) {
    if (std::strcmp(argv[i], "--calibrate") == 0) calibrate_only = true;
    if (std::strcmp(argv[i], "--extra") == 0) extra_only = true;
    if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) threads = static_cast<std::size_t>(std::stoul(argv[++i]));
    else if (std::strcmp(argv[i], "--seeds") == 0 && i + 1 < argc) seeds_n = static_cast<std::size_t>(std::stoul(argv[++i]));
    else if (std::strcmp(argv[i], "--scale") == 0 && i + 1 < argc) scale = std::stod(argv[++i]);
  }
  const std::string out = argv[4];
  const auto T0 = Clock::now();
  const Json in = Json::parse(read_file(argv[1]));
  const Json pj = Json::parse(read_file(argv[2]));
  const Json sj = Json::parse(read_file(argv[3]));
  const economy::Problem P(in, &pj);
  const std::size_t F = P.size(), L = P.links.size();
  flow::Options opt;
  flow::Model M(P, in, sj, opt);
  const auto& N = M.net;
  const hv::Bits none_made(L, 0), all_made(L, 1);
  std::size_t main_edges = 0;
  for (std::size_t e = 0; e < N.edges.size(); e++) main_edges += M.on_main(e) && N.edges[e].link < 0;
  double floor_work = 0, floor_home = 0;
  for (double w : M.work) floor_work += w;
  for (double h : M.home) floor_home += h;
  std::array<double, economy::kinds> attractors{};
  for (std::size_t k = 0; k < economy::kinds; k++)
    for (double a : M.attract[k]) attractors[k] += a;
  std::printf("network: %zu junctions, %zu segments (%zu of the 11 links; %zu on the largest component as built), %zu entrances, %zu counters inside the box\n", N.nodes(),
              N.edges.size(), static_cast<std::size_t>(std::count_if(N.edges.begin(), N.edges.end(), [](const auto& e) { return e.link >= 0; })), main_edges, M.entrances.size(),
              M.counters.size());
  std::printf("floor: %.0f m2 of work, %.0f m2 of homes on %zu junctions; attractors (CLUE 2024, outside the 61 parcels): retail %.0f, food %.0f, services %.0f, culture %.0f; W = %.1f\n",
              floor_work, floor_home, M.building_nodes.size(), attractors[1], attractors[2], attractors[3], attractors[5], M.W_ref);

  // ---- 1. calibration ------------------------------------------------------------------------------------------------
  std::vector<std::size_t> cal_edges;
  for (const auto& c : M.counters) cal_edges.push_back(c.edge);
  for (auto e : M.fronted_edges()) cal_edges.push_back(e);
  std::vector<std::size_t> entrance_edge;
  for (const auto& E : M.entrances)
    for (std::size_t e = 0; e < N.edges.size(); e++)
      if (N.edges[e].link < 0 && (N.edges[e].u == E.node || N.edges[e].v == E.node)) { entrance_edge.push_back(e), cal_edges.push_back(e); break; }
  std::sort(cal_edges.begin(), cal_edges.end());
  cal_edges.erase(std::unique(cal_edges.begin(), cal_edges.end()), cal_edges.end());
  Json grid = Json::array();
  double best_theta = 0, best_loo = -1e9;
  for (double th : {0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.08}) {
    M.set_theta(th);
    if (!(M.radius_all < 1)) {
      std::printf("theta %.3f/m: spectral radius %.3f with every link: the walk sums diverge; skipped\n", th, M.radius_all);
      grid.push(Json::object({{"theta_per_m", th}, {"radius_all_links", jr(M.radius_all, 1e4)}, {"diverges", true}}));
      continue;
    }
    const auto T = M.trips(none_made);
    const auto R = M.terms(T, cal_edges);
    const auto fit = M.calibrate(R, false);
    std::printf("theta %.3f/m: radius %.3f (all links %.3f); gross R2 %.3f, leave-one-counter-out %.3f; in ln(1 + count) %.3f, %.3f; net current R2 %.3f, %.3f\n", th, M.radius_none, M.radius_all,
                fit.r2, fit.loo_r2, fit.r2_log, fit.loo_r2_log, fit.net_r2, fit.loo_net_r2);
    grid.push(Json::object({{"theta_per_m", th}, {"radius_none", jr(M.radius_none, 1e4)}, {"radius_all_links", jr(M.radius_all, 1e4)}, {"r2", jr(fit.r2, 1e4)},
                            {"loo_r2", jr(fit.loo_r2, 1e4)}, {"r2_log", jr(fit.r2_log, 1e4)}, {"loo_r2_log", jr(fit.loo_r2_log, 1e4)}, {"net_r2", jr(fit.net_r2, 1e4)},
                            {"loo_net_r2", jr(fit.loo_net_r2, 1e4)}}));
    if (fit.loo_r2 > best_loo) best_loo = fit.loo_r2, best_theta = th;
  }
  M.set_theta(best_theta);
  const auto Tc = M.trips(none_made);
  const auto Rc = M.terms(Tc, cal_edges);
  // The modes: all, or without the homes' (which may only overfit), chosen by the leave-one-counter-out fit.
  std::array<char, flow::modes> use{1, 1, 1, 1, 1, 1, 1};
  {
    const auto all = M.calibrate(Rc, false), no_homes = M.calibrate(Rc, false, {1, 0, 1, 0, 1, 1, 1});
    if (no_homes.loo_r2 > all.loo_r2) use = {1, 0, 1, 0, 1, 1, 1};
    std::printf("modes: with the homes' leave-one-counter-out %.3f, without %.3f: %s\n", all.loo_r2, no_homes.loo_r2, use[1] ? "kept" : "left out");
  }
  const auto fit = M.calibrate(Rc, true, use);
  std::printf("chosen theta %.3f/m: gross R2 %.3f, leave-one-counter-out %.3f; net current R2 %.3f, leave-one-counter-out %.3f; %zu of %zu counters read each way\n", best_theta, fit.r2, fit.loo_r2,
              fit.net_r2, fit.loo_net_r2, static_cast<std::size_t>(std::count_if(M.counters.begin(), M.counters.end(), [](const auto& c) { return c.directed; })), M.counters.size());
  std::printf("  visits: %.3f of the trips (set: %.2f); footpaths' handedness:", fit.visit_share, M.o.visit_share);
  for (std::size_t c = 0; c < M.counters.size(); c++) std::printf(" %lld %+.2f", M.counters[c].id, fit.hand[c]);
  std::printf("\n");
  for (std::size_t m = 0; m < flow::modes; m++) {
    double sm = 0, mx = 0;
    std::size_t on = 0;
    for (double v : M.x[m]) sm += v, mx = std::max(mx, v), on += v > 0;
    std::printf("  amplitude of %-14s mean %.4g, largest %.4g, positive in %zu hours\n", flow::mode_name(m), sm / H, mx, on);
  }
  // Ablations: modes left out.
  Json ablations = Json::array();
  {
    const auto ab = [&](const char* name, std::array<char, flow::modes> mask) {
      for (std::size_t m = 0; m < flow::modes; m++) mask[m] = mask[m] && use[m];
      const auto f = M.calibrate(Rc, false, mask);
      std::printf("  without %s: R2 %.3f, leave-one-counter-out %.3f\n", name, f.r2, f.loo_r2);
      ablations.push(Json::object({{"without", name}, {"r2", jr(f.r2, 1e4)}, {"loo_r2", jr(f.loo_r2, 1e4)}}));
    };
    ab("visits", {1, 1, 1, 1, 1, 1, 0});
    ab("local trips", {1, 1, 1, 1, 1, 0, 1});
    ab("through trips", {1, 1, 1, 1, 0, 1, 1});
    ab("homes", {1, 0, 1, 0, 1, 1, 1});
  }
  // Widths all alike (the conductance's dependence on width removed).
  double uniform_r2 = 0, uniform_loo = 0;
  {
    flow::Options u = opt;
    u.theta = best_theta;
    for (double* w : {&u.widths.swanston, &u.widths.main, &u.widths.little, &u.widths.street, &u.widths.bourke_mall, &u.widths.mall, &u.widths.laneway, &u.widths.lane_link, &u.widths.passage})
      *w = 6;
    u.widths.reference = 6;
    flow::Model Mu(P, in, sj, u);
    std::vector<std::size_t> ce;
    for (const auto& c : Mu.counters) ce.push_back(c.edge);
    for (auto e : Mu.fronted_edges()) ce.push_back(e);
    std::sort(ce.begin(), ce.end());
    ce.erase(std::unique(ce.begin(), ce.end()), ce.end());
    const auto fu = Mu.calibrate(Mu.terms(Mu.trips(none_made), ce), false, use);
    uniform_r2 = fu.r2, uniform_loo = fu.loo_r2;
    std::printf("  widths all alike: R2 %.3f, leave-one-counter-out %.3f\n", fu.r2, fu.loo_r2);
  }
  // The boundary: the modelled gross flow on each entrance's segment against its counter's.
  Json boundary = Json::array();
  std::vector<double> ratios;
  {
    const auto draw = M.frontage_draw(P.as_built());
    std::vector<double> model(M.entrances.size(), 0.0), given(M.entrances.size(), 0.0);
    for (std::size_t h = 0; h < H; h++) {
      const auto fl = M.mode_flows(Rc, draw, h);
      for (std::size_t i = 0; i < M.entrances.size(); i++) {
        const auto k = static_cast<std::size_t>(std::find(Rc.edges.begin(), Rc.edges.end(), entrance_edge[i]) - Rc.edges.begin());
        for (std::size_t m = 0; m < flow::modes; m++) model[i] += fl[m][k][0] + fl[m][k][1];
        given[i] += M.entrances[i].gross[h];
      }
    }
    for (std::size_t i = 0; i < M.entrances.size(); i++) {
      const auto& E = M.entrances[i];
      ratios.push_back(model[i] / given[i]);
      boundary.push(Json::object({{"street", E.street}, {"kind", E.kind}, {"source", E.source}, {"counter", E.counter >= 0 ? Json(static_cast<double>(E.counter)) : Json()},
                                  {"distance_m", jr(E.distance, 10)}, {"width_m", E.width}, {"given_week", jr(given[i], 1)}, {"model_week", jr(model[i], 1)}}));
    }
    std::vector<double> s = ratios;
    std::sort(s.begin(), s.end());
    std::printf("boundary: modelled over given weekly gross at the %zu entrances: median %.2f, quartiles %.2f to %.2f\n", s.size(), s[s.size() / 2], s[s.size() / 4], s[3 * s.size() / 4]);
  }
  if (calibrate_only) return 0;

  // ---- 2. every set of links -------------------------------------------------------------------------------------
  calx::scheduler::ThreadPool pool(threads);
  flow::Responses resp(M);
  {
    const auto t = Clock::now();
    const std::size_t masks = std::size_t{1} << L;
    calx::scheduler::parallel_for(pool, masks, 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t k = b; k < e; k++) resp(bits_of(k, L));
    });
    std::printf("responses of the passers for all %zu sets of links: %.1f s on %zu threads\n", masks, seconds_since(t), threads);
  }
  const Plan built = P.as_built();
  const auto R0 = resp(none_made);
  const auto pass0 = M.passers(built, *R0);
  const auto O0 = P.evaluate(built, pass0);
  std::printf("as built: V = %.0f A$/week, spending %.0f, footfall %.0f people-passings/week at the frontages, spread %.3f (%.1f h), mix %.3f, access %.1f m\n", O0.value, O0.spending,
              O0.footfall, O0.spread, O0.effective_hours, O0.mix, O0.access);
  const auto measure = [&](const Plan& p) { return P.evaluate(p, M.passers(p, *resp(p.links))); };

  // ---- 3. searches -------------------------------------------------------------------------------------------------
  const economy::Kernels ker{&P};
  const std::size_t large = static_cast<std::size_t>(std::max(1000.0, std::round(20000 * scale))), small = 400;
  std::vector<std::uint64_t> seeds;
  for (std::size_t s = 1; s <= seeds_n; s++) seeds.push_back(s);
  std::vector<Run> runs;
  const auto finish = [&](Run& run) {
    std::printf("  %-22s budget %5zu seed %llu: best %.0f A$/week, %.0f evaluations, %.1f s\n", run.algorithm.c_str(), run.budget, static_cast<unsigned long long>(run.seed), run.best, run.evaluations,
                run.seconds);
    runs.push_back(run);
  };
  // The marginal value of each frontage's decision at the plan as built (no link made): the table of the biased kernels.
  std::vector<std::array<structured::Engine::Best, economy::kinds>> table(F);
  double table_cost = 0;
  {
    const structured::Engine E(P, M);
    for (std::size_t f = 0; f < F; f++) table[f] = E.best_response(built, *R0, f).second;
    table_cost = E.evaluations();
  }
  const double tau = 1000;  // A$/week: the temperature of the kernels' choice of kind (assumption)
  const auto draw_kind = [&](std::size_t f, calx::random::Stream& s) {
    double mx = -std::numeric_limits<double>::infinity();
    for (const auto& b : table[f]) mx = std::max(mx, b.value);
    std::array<double, economy::kinds> w{};
    double sum = 0;
    for (std::size_t k = 0; k < economy::kinds; k++) sum += w[k] = std::exp((table[f][k].value - mx) / tau);
    double u = s.uniform() * sum;
    for (std::size_t k = 0; k < economy::kinds; k++) {
      if (u < w[k]) return k;
      u -= w[k];
    }
    return economy::kinds - 1;
  };
  const auto set_from_table = [&](Plan& p, std::size_t f, std::size_t k) {
    p.kind[f] = static_cast<std::int64_t>(k);
    for (std::size_t j = 0; j < 4; j++) p.hours[4 * f + j] = table[f][k].hours[j];
  };

  // --extra: study (1)'s GA started from the plan as built (its first generation the plan as built mutated five times
  // over), added to study (1) after this study's main runs; run alone, into runs-extra.csv and anytime-extra.csv.
  if (extra_only) {
    std::ostringstream rc, ac;
    rc << "algorithm,budget,seed,evaluations,seconds,best_value,links_made,link_sets_visited\n";
    ac << "algorithm,budget,seed,evaluations,best_value\n";
    for (const std::size_t budget : {std::size_t{400}, static_cast<std::size_t>(std::max(1000.0, std::round(20000 * scale)))}) {
      const bool is_small = budget == 400;
      for (std::uint64_t seed = 1; seed <= seeds_n; seed++) {
        Tape tape;
        heuresis::genetic::GaOptions o;
        o.mu = is_small ? 20 : 100, o.lambda = o.mu, o.elites = 2, o.crossover_rate = 0.9, o.tournament = 2;
        heuresis::genetic::GeneticAlgorithm<Plan> ga({[&](calx::random::Stream& s) {
                                                        Plan p = built;
                                                        for (int k = 0; k < 5; k++) p = ker.mutate(p, s);
                                                        return p;
                                                      },
                                                      [&](const Plan& x, calx::random::Stream& s) { return ker.mutate(x, s); },
                                                      [&](const Plan& a, const Plan& b, calx::random::Stream& s) { return ker.crossover(a, b, s); }},
                                                     o, hs::Draws(seed));
        const auto t0 = Clock::now();
        hs::run(ga, [&](const Plan& p, hs::Ticket t) { const double v = measure(p).value; tape.add(t, v, p); return hs::scalar(-v); }, hs::Budget{static_cast<double>(budget)}, hs::Pool(pool));
        const double secs = seconds_since(t0);
        std::size_t made = 0;
        for (auto b : tape.best_plan.links) made += b;
        rc << "\"GA from as built\"," << budget << "," << seed << "," << tape.values.size() << "," << r(secs, 100) << "," << r(tape.best, 10) << "," << made << ",0\n";
        for (const auto& [e, v] : anytime_of(tape, 0, is_small ? 10 : 100)) ac << "\"GA from as built\"," << budget << "," << seed << "," << r(e, 10) << "," << r(v, 10) << "\n";
        std::printf("  GA from as built       budget %5zu seed %llu: best %.0f A$/week, %zu evaluations, %.1f s\n", budget, static_cast<unsigned long long>(seed), tape.best, tape.values.size(), secs);
      }
    }
    write_file(out + "/runs-extra.csv", rc.str());
    write_file(out + "/anytime-extra.csv", ac.str());
    return 0;
  }

  // Sweeps of best responses until none changes or the budget is spent.
  const auto settle = [&](const structured::Engine& E, Plan& p, const flow::Response& R, calx::random::Stream& s, double limit, Run& run, Plan& best_plan) {
    double v = -std::numeric_limits<double>::infinity();
    for (;;) {
      std::vector<std::size_t> order(F);
      std::iota(order.begin(), order.end(), std::size_t{0});
      for (std::size_t i = F; i-- > 1;) std::swap(order[i], order[s.below(i + 1)]);
      std::size_t changed = 0;
      bool out_of_budget = false;
      for (auto f : order) {
        if (E.evaluations() + 7 > limit) { out_of_budget = true; break; }
        const auto G = E.gains(p, R, f);
        const auto [best, per] = E.choose(p, f, G);
        (void)per;
        if (best.value > E.present(p, f, G) + 1e-6) {
          p.kind[f] = best.kind;
          for (std::size_t j = 0; j < 4; j++) p.hours[4 * f + j] = best.hours[j];
          changed++;
        }
      }
      v = E.value(p, R);
      if (v > run.best) run.best = v, best_plan = p, run.anytime.push_back({E.evaluations(), v});
      if (changed == 0 || out_of_budget || E.evaluations() + 7 > limit) return v;
    }
  };

  for (const std::size_t budget : {small, large}) {
    const bool is_small = budget == small;
    const std::size_t step = is_small ? 10 : 100;
    for (auto seed : seeds) {
      // The black box, as study (1) runs it. GA on the mixed encoding.
      {
        Tape tape;
        Run run{"GA", budget, seed};
        heuresis::genetic::GaOptions o;
        o.mu = is_small ? 20 : 100, o.lambda = o.mu, o.elites = 2, o.crossover_rate = 0.9, o.tournament = 2;
        heuresis::genetic::GeneticAlgorithm<Plan> ga({[&](calx::random::Stream& s) { return ker.sample(s); }, [&](const Plan& x, calx::random::Stream& s) { return ker.mutate(x, s); },
                                                      [&](const Plan& a, const Plan& b, calx::random::Stream& s) { return ker.crossover(a, b, s); }},
                                                     o, hs::Draws(seed));
        const auto t0 = Clock::now();
        hs::run(ga, [&](const Plan& p, hs::Ticket t) { const double v = measure(p).value; tape.add(t, v, p); return hs::scalar(-v); }, hs::Budget{static_cast<double>(budget)}, hs::Pool(pool));
        run.seconds = seconds_since(t0), run.evaluations = static_cast<double>(tape.values.size()), run.best = tape.best, run.best_plan = tape.best_plan, run.anytime = anytime_of(tape, 0, step);
        finish(run);
      }
      // Separable CMA-ES on the canonical relaxation, restarted from a new random mean when it stops.
      {
        Tape tape;
        Run run{"CMA-ES", budget, seed};
        const std::size_t D = P.dimension();
        const auto t0 = Clock::now();
        std::uint64_t serial = 0, restart = 0;
        double spent = 0;
        while (spent < static_cast<double>(budget)) {
          heuresis::strategies::CmaOptions o;
          o.separable = true, o.bounds = hv::Box::cube(D, 0, 1);
          auto ms = hs::Draws(seed)(hs::tag("study.cma.mean"), restart, 0);
          hv::Real mean(D);
          for (auto& x : mean) x = 0.1 + 0.8 * ms.uniform();
          heuresis::strategies::Cma cma(mean, 0.3, o, hs::Draws(seed * 1000 + restart), serial);
          const auto oc = hs::run(cma, [&](const hv::Real& u, hs::Ticket t) { const Plan p = P.decode(u); const double v = measure(p).value; tape.add(t, v, p); return hs::scalar(-v); },
                                  hs::Budget{static_cast<double>(budget) - spent}, hs::Pool(pool));
          spent += oc.evaluations, serial = cma.next_serial(), restart++;
          if (oc.evaluations == 0) break;
        }
        run.seconds = seconds_since(t0), run.evaluations = static_cast<double>(tape.values.size()), run.best = tape.best, run.best_plan = tape.best_plan, run.anytime = anytime_of(tape, 0, step);
        finish(run);
      }
      // Random search.
      {
        Tape tape;
        Run run{"random", budget, seed};
        hs::RandomSearch<Plan> rs([&](calx::random::Stream& s) { return ker.sample(s); }, hs::Draws(seed), 32);
        const auto t0 = Clock::now();
        hs::run(rs, [&](const Plan& p, hs::Ticket t) { const double v = measure(p).value; tape.add(t, v, p); return hs::scalar(-v); }, hs::Budget{static_cast<double>(budget)}, hs::Pool(pool), 32);
        run.seconds = seconds_since(t0), run.evaluations = static_cast<double>(tape.values.size()), run.best = tape.best, run.best_plan = tape.best_plan, run.anytime = anytime_of(tape, 0, step);
        finish(run);
      }
      // The structure. GA with kernels biased by the marginal value of each frontage's decision.
      {
        Tape tape;
        Run run{"GA, marginal kernels", budget, seed};
        heuresis::genetic::GaOptions o;
        o.mu = is_small ? 20 : 100, o.lambda = o.mu, o.elites = 2, o.crossover_rate = 0.9, o.tournament = 2;
        const auto sample = [&](calx::random::Stream& s) {
          Plan p = ker.sample(s);
          for (std::size_t f = 0; f < F; f++)
            if (s.uniform() < 0.5) set_from_table(p, f, draw_kind(f, s));
          return p;
        };
        const auto mutate = [&](const Plan& x, calx::random::Stream& s) {
          Plan y = ker.mutate(x, s);
          for (std::size_t f = 0; f < F; f++)
            if (s.uniform() < 1.0 / static_cast<double>(F)) set_from_table(y, f, draw_kind(f, s));
          return y;
        };
        heuresis::genetic::GeneticAlgorithm<Plan> ga({sample, mutate, [&](const Plan& a, const Plan& b, calx::random::Stream& s) { return ker.crossover(a, b, s); }}, o, hs::Draws(seed));
        const auto t0 = Clock::now();
        const double left = std::max(0.0, static_cast<double>(budget) - std::ceil(table_cost));
        hs::run(ga, [&](const Plan& p, hs::Ticket t) { const double v = measure(p).value; tape.add(t, v, p); return hs::scalar(-v); }, hs::Budget{left}, hs::Pool(pool));
        run.seconds = seconds_since(t0), run.evaluations = std::ceil(table_cost) + static_cast<double>(tape.values.size()), run.best = tape.best, run.best_plan = tape.best_plan;
        run.anytime = anytime_of(tape, std::ceil(table_cost), step);
        finish(run);
      }
    }
    // Best responses, with a descent over the links, and with a (1+1) EA over the links: one seed a thread.
    std::vector<Run> br(seeds.size()), brea(seeds.size());
    calx::scheduler::parallel_for(pool, seeds.size(), 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t i = b; i < e; i++) {
        const std::uint64_t seed = seeds[i];
        // BR with a descent over the links.
        {
          Run run{"BR", budget, seed};
          const structured::Engine E(P, M);
          auto s = hs::Draws(seed)(hs::tag("study.br"), 0, 0);
          const auto t0 = Clock::now();
          Plan p = built, best_plan = built;
          const double limit = static_cast<double>(budget);
          double v = settle(E, p, *resp(p.links), s, limit, run, best_plan);
          run.inner = 1;
          bool improved = true;
          while (improved && E.evaluations() + 7 <= limit) {
            improved = false;
            std::vector<std::size_t> order(L);
            std::iota(order.begin(), order.end(), std::size_t{0});
            for (std::size_t k = L; k-- > 1;) std::swap(order[k], order[s.below(k + 1)]);
            for (auto k : order) {
              if (E.evaluations() + 7 > limit) break;
              Plan q = p;
              q.links[k] ^= 1;
              const double vq = settle(E, q, *resp(q.links), s, limit, run, best_plan);
              run.inner++;
              if (vq > v + 1e-6) p = q, v = vq, improved = true;
            }
          }
          run.seconds = seconds_since(t0), run.evaluations = E.evaluations(), run.best_plan = best_plan;
          run.anytime.push_back({run.evaluations, run.best});
          br[i] = std::move(run);
        }
        // BR with a (1+1) EA over the links.
        {
          Run run{"BR + (1+1) EA", budget, seed};
          const structured::Engine E(P, M);
          auto s = hs::Draws(seed)(hs::tag("study.br.inner"), 0, 0);
          const auto t0 = Clock::now();
          Plan incumbent = built, best_plan = built;
          double incumbent_value = -std::numeric_limits<double>::infinity();
          const double limit = static_cast<double>(budget);
          heuresis::genetic::OnePlusLambda<hv::Bits> ea([&](calx::random::Stream&) { return hv::Bits(L, 0); },
                                                        [&](const hv::Bits& x, calx::random::Stream& st) {
                                                          hv::Bits y = x;
                                                          hv::standard_bit_mutation_positive(y, 1.0 / static_cast<double>(L), st);
                                                          return y;
                                                        },
                                                        1, hs::Draws(seed));
          hs::run(ea, [&](const hv::Bits& links) {
                    const double before = E.evaluations();
                    Plan q = incumbent;
                    q.links = links;
                    const double v = settle(E, q, *resp(links), s, limit, run, best_plan);
                    run.inner++;
                    if (v >= incumbent_value) incumbent = q, incumbent_value = v;
                    return hs::scalar(-v, E.evaluations() - before);
                  }, hs::Budget{std::numeric_limits<double>::infinity(), limit - 7});
          run.seconds = seconds_since(t0), run.evaluations = E.evaluations(), run.best_plan = best_plan;
          run.anytime.push_back({run.evaluations, run.best});
          brea[i] = std::move(run);
        }
      }
    });
    for (auto& run : br) {
      // The best plan's value read again by the problem's own evaluate.
      run.best = measure(run.best_plan).value;
      finish(run);
    }
    for (auto& run : brea) {
      run.best = measure(run.best_plan).value;
      finish(run);
    }
  }

  // ---- reading the runs ------------------------------------------------------------------------------------------------
  std::ostringstream rcsv, acsv;
  rcsv << "algorithm,budget,seed,evaluations,seconds,best_value,links_made,link_sets_visited\n";
  acsv << "algorithm,budget,seed,evaluations,best_value\n";
  std::map<std::string, std::vector<double>> finals;
  const Run* best_run = nullptr;
  const Run* best_black = nullptr;
  for (const auto& run : runs) {
    std::size_t made = 0;
    for (auto b : run.best_plan.links) made += b;
    rcsv << '"' << run.algorithm << "\"," << run.budget << "," << run.seed << "," << r(run.evaluations, 10) << "," << r(run.seconds, 100) << "," << r(run.best, 10) << "," << made << "," << run.inner << "\n";
    for (const auto& [e, v] : run.anytime) acsv << '"' << run.algorithm << "\"," << run.budget << "," << run.seed << "," << r(e, 10) << "," << r(v, 10) << "\n";
    finals[run.algorithm + "|" + std::to_string(run.budget)].push_back(run.best);
    if (run.budget == large) {
      if (!best_run || run.best > best_run->best) best_run = &run;
      if ((run.algorithm == "GA" || run.algorithm == "CMA-ES" || run.algorithm == "random") && (!best_black || run.best > best_black->best)) best_black = &run;
    }
  }
  write_file(out + "/runs.csv", rcsv.str());
  write_file(out + "/anytime.csv", acsv.str());
  Json table_json = Json::array();
  for (const auto& [key, v] : finals) {
    double mean = 0, sd = 0;
    for (double x : v) mean += x / static_cast<double>(v.size());
    for (double x : v) sd += (x - mean) * (x - mean) / static_cast<double>(std::max<std::size_t>(1, v.size() - 1));
    const auto bar = key.find('|');
    table_json.push(Json::object({{"algorithm", key.substr(0, bar)}, {"budget", std::stod(key.substr(bar + 1))}, {"mean", jr(mean, 1)}, {"sd", jr(std::sqrt(sd), 1)},
                                  {"best", jr(*std::max_element(v.begin(), v.end()), 1)}, {"worst", jr(*std::min_element(v.begin(), v.end()), 1)}}));
    std::printf("%-22s budget %5s: mean %.0f ± %.0f A$/week (best %.0f)\n", key.substr(0, bar).c_str(), key.substr(bar + 1).c_str(), mean, std::sqrt(sd), *std::max_element(v.begin(), v.end()));
  }

  // ---- 4. the flows of the plans --------------------------------------------------------------------------------------
  const Plan best_plan = best_run->best_plan;
  const Plan black_plan = best_black->best_plan;
  std::printf("best plan: %s seed %llu, V = %.0f A$/week; best black-box plan: %s seed %llu, V = %.0f\n", best_run->algorithm.c_str(), static_cast<unsigned long long>(best_run->seed), best_run->best,
              best_black->algorithm.c_str(), static_cast<unsigned long long>(best_black->seed), best_black->best);
  std::vector<std::size_t> all_edges(N.edges.size());
  std::iota(all_edges.begin(), all_edges.end(), std::size_t{0});
  struct PlanFlows {
    std::string name;
    Plan plan;
    economy::Outcome outcome;
    std::vector<double> passers;
    std::vector<std::array<std::array<double, 2>, H>> by_edge;  // per edge, per hour, (u → v, v → u)
    std::array<std::vector<double>, flow::modes> mode_gross;    // per mode, per edge: the week's gross
  };
  const auto flows_of = [&](const std::string& name, const Plan& p) {
    PlanFlows pf;
    pf.name = name, pf.plan = p;
    pf.passers = M.passers(p, *resp(p.links));
    pf.outcome = P.evaluate(p, pf.passers);
    std::vector<std::size_t> present;
    for (auto e : all_edges)
      if (N.present(e, p.links)) present.push_back(e);
    const auto T = M.trips(p.links);
    const auto R = M.terms(T, present);
    const auto draw = M.frontage_draw(p);
    pf.by_edge.assign(N.edges.size(), {});
    for (auto& m : pf.mode_gross) m.assign(N.edges.size(), 0.0);
    for (std::size_t h = 0; h < H; h++) {
      const auto fl = M.mode_flows(R, draw, h);
      for (std::size_t k = 0; k < present.size(); k++)
        for (std::size_t m = 0; m < flow::modes; m++) {
          pf.by_edge[present[k]][h][0] += fl[m][k][0], pf.by_edge[present[k]][h][1] += fl[m][k][1];
          pf.mode_gross[m][present[k]] += fl[m][k][0] + fl[m][k][1];
        }
    }
    return pf;
  };
  std::vector<PlanFlows> plans;
  plans.push_back(flows_of("as built", built));
  plans.push_back(flows_of("best", best_plan));
  plans.push_back(flows_of("best black box", black_plan));
  std::printf("flows on every segment for %zu plans: %.0f s since the start\n", plans.size(), seconds_since(T0));

  // The circulation of the net flow, hour by hour and at chosen hours.
  const std::vector<std::pair<std::string, std::size_t>> moments{{"Wednesday 08:00", 2 * 24 + 8}, {"Wednesday 12:00", 2 * 24 + 12}, {"Wednesday 17:00", 2 * 24 + 17}, {"Saturday 14:00", 5 * 24 + 14}};
  Json plan_json = Json::array(), maps = Json::array();
  for (const auto& pf : plans) {
    Json shares = Json::array();
    std::vector<std::array<double, 4>> share_h(H);  // driven, local, global, energy
    for (std::size_t h = 0; h < H; h++) {
      std::vector<std::array<double, 2>> fl(N.edges.size());
      for (std::size_t e = 0; e < N.edges.size(); e++) fl[e] = pf.by_edge[e][h];
      const auto C = flow::circulation(N, pf.plan.links, fl);
      share_h[h] = {C.parts.driven, C.parts.local, C.parts.global, C.parts.total};
      shares.push(Json::array({r(C.parts.driven, 1e4), r(C.parts.local, 1e4), r(C.parts.global, 1e4), r(C.parts.total, 1)}));
    }
    // Per edge: the week's gross, the lines of the gross, and at the moments the net flow and its parts.
    Json edges = Json::array();
    std::map<std::string, flow::Circulation> at;
    for (const auto& [label, h] : moments) {
      std::vector<std::array<double, 2>> fl(N.edges.size());
      for (std::size_t e = 0; e < N.edges.size(); e++) fl[e] = pf.by_edge[e][h];
      at.emplace(label, flow::circulation(N, pf.plan.links, fl));
    }
    {
      // The week's mean flow too.
      std::vector<std::array<double, 2>> fl(N.edges.size(), {0.0, 0.0});
      for (std::size_t e = 0; e < N.edges.size(); e++)
        for (std::size_t h = 0; h < H; h++) fl[e][0] += pf.by_edge[e][h][0] / H, fl[e][1] += pf.by_edge[e][h][1] / H;
      at.emplace("the week's mean", flow::circulation(N, pf.plan.links, fl));
    }
    Json moments_json = Json::object();
    for (auto& [label, C] : at) {
      std::vector<double> netv(N.edges.size(), 0.0), drv(N.edges.size(), 0.0), loc(N.edges.size(), 0.0), glo(N.edges.size(), 0.0), pot(N.nodes(), 0.0);
      for (std::size_t k = 0; k < C.edges.size(); k++)
        netv[C.edges[k]] = C.net[k], drv[C.edges[k]] = C.parts.parts.driven[k], loc[C.edges[k]] = C.parts.parts.local[k], glo[C.edges[k]] = C.parts.parts.global[k];
      Json faces = Json::array();
      for (std::size_t q = 0; q < C.faces.size(); q++) {
        Json ring = Json::array();
        for (const auto& [k, sg] : C.faces[q]) ring.push(Json::array({static_cast<double>(C.edges[k]), static_cast<double>(sg)}));
        faces.push(Json::object({{"edges", ring}, {"circulation", jr(C.parts.parts.circulation[q], 1000)}, {"area_m2", jr(C.face_area[q], 1)}}));
      }
      moments_json.set(label, Json::object({{"net", numbers(netv, 2)}, {"driven", numbers(drv, 2)}, {"local", numbers(loc, 2)}, {"global", numbers(glo, 2)},
                                            {"shares", Json::array({r(C.parts.driven, 1e4), r(C.parts.local, 1e4), r(C.parts.global, 1e4)})}, {"energy", jr(C.parts.total, 1)},
                                            {"faces", faces}, {"holes", static_cast<double>(C.holes)}}));
    }
    std::vector<double> week_gross(N.edges.size(), 0.0);
    std::vector<double> d1a(N.edges.size(), 0.0), d1c(N.edges.size(), 0.0), w1a(N.edges.size(), 0.0), w1c(N.edges.size(), 0.0), d2a(N.edges.size(), 0.0), mean(N.edges.size(), 0.0);
    std::vector<double> noon(N.edges.size(), 0.0), evening(N.edges.size(), 0.0), saturday(N.edges.size(), 0.0), morning(N.edges.size(), 0.0);
    for (std::size_t e = 0; e < N.edges.size(); e++) {
      std::vector<double> g(H);
      for (std::size_t h = 0; h < H; h++) g[h] = pf.by_edge[e][h][0] + pf.by_edge[e][h][1], week_gross[e] += g[h];
      const auto l7 = flow::line(g, 7), l1 = flow::line(g, 1), l14 = flow::line(g, 14);
      mean[e] = week_gross[e] / H, d1a[e] = flow::amplitude(l7), d1c[e] = flow::crest_hour(l7, 7), w1a[e] = flow::amplitude(l1), w1c[e] = flow::crest_hour(l1, 1), d2a[e] = flow::amplitude(l14);
      morning[e] = g[2 * 24 + 8], noon[e] = g[2 * 24 + 12], evening[e] = g[2 * 24 + 17], saturday[e] = g[5 * 24 + 14];
    }
    Json modes_j = Json::object();
    for (std::size_t m = 0; m < flow::modes; m++) modes_j.set(flow::mode_name(m), numbers(pf.mode_gross[m], 1));
    maps.push(Json::object({{"plan", pf.name}, {"links", economy::plan_json(pf.plan)["links"]}, {"week_gross", numbers(week_gross, 1)}, {"mean_gross", numbers(mean, 2)},
                            {"d1_amplitude", numbers(d1a, 2)}, {"d1_crest_hour", numbers(d1c, 3)}, {"w1_amplitude", numbers(w1a, 2)}, {"w1_crest_hour", numbers(w1c, 3)},
                            {"d2_amplitude", numbers(d2a, 2)}, {"gross_wed_08", numbers(morning, 2)}, {"gross_wed_12", numbers(noon, 2)}, {"gross_wed_17", numbers(evening, 2)},
                            {"gross_sat_14", numbers(saturday, 2)}, {"mode_week_gross", modes_j}, {"moments", moments_json}}));
    // The plan's outcome, its frontages' passers and lines, and its shares of circulation by hour.
    Json fr = Json::array();
    for (std::size_t f = 0; f < F; f++) {
      const std::span<const double> pf_h(pf.passers.data() + f * H, H);
      double week = 0;
      for (double v : pf_h) week += v;
      const auto l7 = flow::line(pf_h, 7), l1 = flow::line(pf_h, 1);
      fr.push(Json::object({{"passers_week", jr(week, 1)}, {"d1_amplitude", jr(flow::amplitude(l7), 100)}, {"d1_crest_hour", jr(flow::crest_hour(l7, 7), 1000)},
                            {"w1_amplitude", jr(flow::amplitude(l1), 100)}, {"spending_week", jr(pf.outcome.spending_by_frontage[f], 1)}, {"value_week", jr(pf.outcome.value_by_frontage[f], 1)},
                            {"visits_week", jr(pf.outcome.visits_by_frontage[f], 10)}}));
    }
    std::vector<double> totals(H, 0.0);
    for (std::size_t f = 0; f < F; f++)
      for (std::size_t h = 0; h < H; h++) totals[h] += pf.passers[f * H + h];
    std::vector<double> open_by_hour(H, 0.0);
    for (std::size_t h = 0; h < H; h++)
      for (std::size_t f = 0; f < F; f++) open_by_hour[h] += pf.plan.open(f, h) && economy::trades(pf.plan.kind[f]);
    Json lines = Json::object();
    for (std::size_t m : {1ul, 2ul, 6ul, 7ul, 8ul, 14ul, 21ul}) {
      const auto cpass = flow::line(totals, m);
      const auto cspend = flow::line(std::span<const double>(pf.outcome.spending_by_hour.data(), H), m);
      lines.set(std::to_string(m), Json::object({{"passers_amplitude", jr(flow::amplitude(cpass), 1)}, {"passers_crest_hour", jr(flow::crest_hour(cpass, m), 1000)},
                                                 {"spending_amplitude", jr(flow::amplitude(cspend), 10)}, {"spending_crest_hour", jr(flow::crest_hour(cspend, m), 1000)}}));
    }
    plan_json.push(Json::object({{"plan", pf.name}, {"genotype", economy::plan_json(pf.plan)}, {"outcome", economy::outcome_json(pf.outcome, true)}, {"frontages", fr},
                                 {"passers_by_hour", numbers(totals, 1)}, {"open_trading_by_hour", numbers(open_by_hour, 0)}, {"lines", lines}, {"circulation_by_hour", shares}}));
    double dm = 0, lm = 0, gm = 0;
    for (const auto& s : share_h) dm += s[0] / H, lm += s[1] / H, gm += s[2] / H;
    const auto& wk = at.at("the week's mean");
    std::printf("%s: V = %.0f; circulation of the net flow, mean over the hours: driven %.3f, round the blocks %.3f, round the holes %.3f; of the week's mean flow %.3f, %.3f, %.3f\n",
                pf.name.c_str(), pf.outcome.value, dm, lm, gm, wk.parts.driven, wk.parts.local, wk.parts.global);
  }

  // What each link does: alone, with the plan as built (the flows rerouted, the plan not changed), and from the best
  // plan, flipped and the frontages re-optimised by best responses.
  Json link_json = Json::array();
  {
    const structured::Engine E(P, M);
    const double v_best = E.value(best_plan, *resp(best_plan.links));
    std::vector<double> base_pass(F, 0.0);
    for (std::size_t f = 0; f < F; f++)
      for (std::size_t h = 0; h < H; h++) base_pass[f] += pass0[f * H + h];
    for (std::size_t k = 0; k < L; k++) {
      Plan q = built;
      q.links[k] = 1;
      const auto Rq = resp(q.links);
      const auto pq = M.passers(q, *Rq);
      const auto Oq = P.evaluate(q, pq);
      double gained = 0, lost = 0;
      for (std::size_t f = 0; f < F; f++) {
        double w = 0;
        for (std::size_t h = 0; h < H; h++) w += pq[f * H + h];
        const double d = w - base_pass[f];
        (d > 0 ? gained : lost) += d;
      }
      // The flow on the link itself (Wednesday 17:00) and the circulation's shares then.
      std::vector<std::size_t> present;
      for (std::size_t e = 0; e < N.edges.size(); e++)
        if (N.present(e, q.links)) present.push_back(e);
      const auto Tq = M.trips(q.links);
      const auto Rt = M.terms(Tq, present);
      const auto draw = M.frontage_draw(q);
      double on_link_week = 0;
      std::vector<std::array<double, 2>> fl17(N.edges.size(), {0.0, 0.0});
      for (std::size_t h = 0; h < H; h++) {
        const auto fl = M.mode_flows(Rt, draw, h);
        for (std::size_t kk = 0; kk < present.size(); kk++) {
          double fw = 0, bw = 0;
          for (std::size_t m = 0; m < flow::modes; m++) fw += fl[m][kk][0], bw += fl[m][kk][1];
          if (N.edges[present[kk]].link == static_cast<int>(k)) on_link_week += fw + bw;
          if (h == 2 * 24 + 17) fl17[present[kk]] = {fw, bw};
        }
      }
      const auto C17 = flow::circulation(N, q.links, fl17);
      // From the best plan: flip the link and settle.
      Plan b = best_plan;
      b.links[k] ^= 1;
      auto s = hs::Draws(77)(hs::tag("study.links"), k, 0);
      Run scratch;
      Plan bp = b;
      const double vb = settle(E, b, *resp(b.links), s, E.evaluations() + 20000, scratch, bp);
      link_json.push(Json::object({{"link", static_cast<double>(k)}, {"kind", P.links[k].kind}, {"joins", Json::array({P.links[k].from, P.links[k].to})}, {"length_m", jr(P.links[k].length, 10)},
                                   {"weekly_cost_aud", jr(P.P.link_weekly(P.links[k].length), 10)}, {"on_link_week", jr(on_link_week, 1)}, {"frontage_passers_gained_week", jr(gained, 1)},
                                   {"frontage_passers_lost_week", jr(lost, 1)}, {"value_as_built_plus_link", jr(Oq.value, 1)}, {"value_change_as_built", jr(Oq.value - O0.value, 1)},
                                   {"in_best_plan", static_cast<double>(best_plan.links[k])}, {"value_change_from_best_flipped_and_settled", jr(vb - v_best, 1)},
                                   {"circulation_wed_17", Json::array({r(C17.parts.driven, 1e4), r(C17.parts.local, 1e4), r(C17.parts.global, 1e4)})}}));
      std::printf("link %2zu (%s, %s to %s): %.0f walkers a week on it; frontages gain %.0f and lose %.0f passings; V as built + link %+.0f; from the best plan flipped %+.0f\n", k, P.links[k].kind.c_str(),
                  P.links[k].from.c_str(), P.links[k].to.c_str(), on_link_week, gained, -lost, Oq.value - O0.value, vb - v_best);
    }
  }

  // Catalysis's part in each plan's value (the same flows, the affinities set to zero), and the arrangement used.
  Json arrangement = Json::array();
  {
    const economy::Problem P0(in);  // no affinities: catalysis off
    for (const auto& pf : plans) {
      const auto o0 = P0.evaluate(pf.plan, pf.passers);
      std::size_t same_segment = 0, catalysing = 0, party_same = 0, party_trading = 0;
      const auto D = P.distances(pf.plan.links);
      for (std::size_t f = 0; f < F; f++) {
        if (!economy::trades(pf.plan.kind[f])) continue;
        for (auto g : P.on_segment[f])
          if (g > f && pf.plan.kind[g] == pf.plan.kind[f]) same_segment++;
        for (const auto& [g, w] : D->near[f])
          if (g > f && economy::trades(pf.plan.kind[g]) && P.P.affinity[static_cast<std::size_t>(pf.plan.kind[g])][static_cast<std::size_t>(pf.plan.kind[f])] > 0) catalysing++;
      }
      for (const auto& pw : P.party) {
        const auto a = pf.plan.kind[pw[0]], b = pf.plan.kind[pw[1]];
        if (economy::trades(a) && economy::trades(b)) party_trading++, party_same += a == b;
      }
      arrangement.push(Json::object({{"plan", pf.name}, {"value", jr(pf.outcome.value, 1)}, {"value_without_catalysis", jr(o0.value, 1)}, {"same_kind_same_segment_pairs", static_cast<double>(same_segment)},
                                     {"positively_catalysing_pairs_in_reach", static_cast<double>(catalysing)}, {"party_walls_both_trading", static_cast<double>(party_trading)},
                                     {"party_walls_same_kind", static_cast<double>(party_same)}}));
      std::printf("%s: catalysis adds %.0f A$/week; %zu same-kind pairs on one segment, %zu positively catalysing pairs within reach, %zu of %zu party walls between trading frontages (%zu of one kind)\n",
                  pf.name.c_str(), pf.outcome.value - o0.value, same_segment, catalysing, party_trading, P.party.size(), party_same);
    }
  }

  // ---- the record -------------------------------------------------------------------------------------------------------
  Json counters = Json::array();
  for (std::size_t c = 0; c < M.counters.size(); c++) {
    const auto& k = M.counters[c];
    Json lines = Json::object();
    for (std::size_t m : {1ul, 7ul, 14ul}) {
      const auto g = k.gross();
      const auto lo = flow::line(g, m), lf = flow::line(fit.fitted[c], m);
      lines.set(std::to_string(m), Json::object({{"observed_amplitude", jr(flow::amplitude(lo), 10)}, {"observed_crest_hour", jr(flow::crest_hour(lo, m), 1000)},
                                                 {"model_amplitude", jr(flow::amplitude(lf), 10)}, {"model_crest_hour", jr(flow::crest_hour(lf, m), 1000)}}));
    }
    counters.push(Json::object({{"id", static_cast<double>(k.id)}, {"name", k.name}, {"street", N.edges[k.edge].name}, {"edge", static_cast<double>(k.edge)}, {"share", k.share},
                                {"observed", numbers(k.gross(), 2)}, {"mean_week", numbers(k.week, 2)}, {"fitted", numbers(fit.fitted[c], 2)}, {"loo", numbers(fit.loo[c], 2)},
                                {"r2", jr(fit.counter_r2[c], 1e4)}, {"loo_r2", jr(fit.counter_loo_r2[c], 1e4)}, {"lines", lines}, {"directed", k.directed},
                                {"toward", Json::array({k.toward[0], k.toward[1]})}, {"cosine", jr(k.cosine, 1000)}, {"way_observed", Json::array({numbers(k.way[0], 2), numbers(k.way[1], 2)})},
                                {"way_fitted", Json::array({numbers(fit.fitted_way[c][0], 2), numbers(fit.fitted_way[c][1], 2)})}, {"hand", jr(fit.hand[c], 1e4)}}));
  }
  Json amps = Json::object();
  for (std::size_t m = 0; m < flow::modes; m++) amps.set(flow::mode_name(m), numbers(M.x[m], 6));
  write_file(out + "/calibration.json",
             Json::object({{"theta_grid", grid}, {"theta_per_m", best_theta}, {"r2", jr(fit.r2, 1e4)}, {"loo_r2", jr(fit.loo_r2, 1e4)}, {"r2_log", jr(fit.r2_log, 1e4)},
                           {"loo_r2_log", jr(fit.loo_r2_log, 1e4)}, {"net_r2", jr(fit.net_r2, 1e4)}, {"loo_net_r2", jr(fit.loo_net_r2, 1e4)},
                           {"visit_share", jr(fit.visit_share, 1e4)}, {"visit_share_set", M.o.visit_share}, {"modes_used", Json::array({static_cast<double>(use[0]), static_cast<double>(use[1]), static_cast<double>(use[2]), static_cast<double>(use[3]), static_cast<double>(use[4]),
                                                       static_cast<double>(use[5]), static_cast<double>(use[6])})}, {"ablations", ablations}, {"uniform_widths", Json::object({{"r2", jr(uniform_r2, 1e4)}, {"loo_r2", jr(uniform_loo, 1e4)}})},
                           {"counters", counters}, {"amplitudes", amps}, {"boundary", boundary}, {"W_ref", jr(M.W_ref, 100)},
                           {"floor_m2", Json::object({{"work", jr(floor_work, 1)}, {"home", jr(floor_home, 1)}})},
                           {"attractors", Json::object({{"retail", attractors[1]}, {"food and drink", attractors[2]}, {"services", attractors[3]}, {"culture", attractors[5]}})},
                           {"radius", Json::object({{"none", jr(M.radius_none, 1e4)}, {"all_links", jr(M.radius_all, 1e4)}})}})
                 .dump());
  write_file(out + "/plans.json", Json::object({{"plans", plan_json}, {"links", link_json}, {"arrangement", arrangement}, {"table_tau_aud_week", tau}}).dump());
  // The maps hold the segments' geometry (OpenStreetMap, ODbL): not committed.
  {
    Json edges = Json::array();
    for (const auto& e : N.edges) {
      Json line = Json::array();
      for (const auto& q : e.line) line.push(Json::array({r(q[0], 100), r(q[1], 100)}));
      edges.push(Json::object({{"u", static_cast<double>(e.u)}, {"v", static_cast<double>(e.v)}, {"line", line}, {"kind", e.kind}, {"name", e.name}, {"width_m", e.width},
                               {"length_m", jr(e.length, 100)}, {"link", static_cast<double>(e.link)}, {"main", static_cast<double>(M.on_main(&e - N.edges.data()))}}));
    }
    Json nodes = Json::array();
    for (std::size_t j = 0; j < N.nodes(); j++) nodes.push(Json::array({r(N.at[j][0], 100), r(N.at[j][1], 100), r(M.work[j], 1), r(M.home[j], 1)}));
    Json ent = Json::array();
    for (const auto& E : M.entrances) ent.push(Json::object({{"node", static_cast<double>(E.node)}, {"street", E.street}, {"source", E.source}, {"gross", numbers(E.gross, 1)}}));
    Json att = Json::array();
    for (auto j : M.attractor_nodes) {
      Json a = Json::array();
      for (std::size_t k = 0; k < economy::kinds; k++) a.push(Json(M.attract[k][j]));
      att.push(Json::object({{"node", static_cast<double>(j)}, {"count", a}}));
    }
    Json fe = Json::array();
    for (std::size_t f = 0; f < F; f++) fe.push(Json::object({{"edge", static_cast<double>(M.frontage_edge[f])}, {"anchor", numbers(std::vector<double>{P.frontages[f].anchor[0], P.frontages[f].anchor[1]}, 2)}}));
    Json cs = Json::array();
    for (const auto& c : M.counters) cs.push(Json::object({{"id", static_cast<double>(c.id)}, {"at", numbers(std::vector<double>{c.at[0], c.at[1]}, 2)}, {"edge", static_cast<double>(c.edge)}}));
    write_file(out + "/flow-maps.json", Json::object({{"edges", edges}, {"nodes", nodes}, {"entrances", ent}, {"attractors", att}, {"frontages", fe}, {"counters", cs}, {"plans", maps},
                                                      {"credit", "(c) OpenStreetMap contributors, ODbL 1.0; City of Melbourne (CC BY 4.0)"}})
                                            .dump());
  }
  Json summary = Json::object({{"as_built", economy::outcome_json(O0)}, {"best", economy::outcome_json(plans[1].outcome)}, {"best_black_box", economy::outcome_json(plans[2].outcome)},
                               {"best_algorithm", best_run->algorithm}, {"best_seed", static_cast<double>(best_run->seed)}, {"best_black_box_algorithm", best_black->algorithm},
                               {"budgets", Json::array({static_cast<double>(small), static_cast<double>(large)})}, {"seeds", static_cast<double>(seeds_n)}, {"algorithms", table_json},
                               {"table_cost_evaluations", jr(table_cost, 100)}, {"theta_per_m", best_theta}, {"calibration_r2", jr(fit.r2, 1e4)}, {"calibration_loo_r2", jr(fit.loo_r2, 1e4)},
                               {"seconds", jr(seconds_since(T0), 10)}, {"threads", static_cast<double>(threads)}});
  write_file(out + "/summary.json", summary.dump());
  std::printf("done in %.0f s\n", seconds_since(T0));
  return 0;
}
