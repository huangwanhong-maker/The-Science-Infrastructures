// A body learns to move: the rabbit of the robotics domain, hopping on the Moon, its controller's gains tuned by
// heuresis's searches (ADR 0058's body; heuresis decision 0003's protocol).
//
// Modes (each writes to the results directory):
//   check    the experiment's loop (hopper.hpp) against the domain's Simulation with the default gains, bit for bit;
//            the cost of a trial
//   explore  the default gains over a grid of commands: what the build does untuned
//   tune     CMA-ES, DE, PSO, the (1+1)-ES and random search on ten gains, 10 seeds each, one task
//   front    NSGA-II and SMS-EMOA on the gains and the speed asked: speed against the cost of transport
//   motions  the bodies' motions for chosen points of the front, and for the default and the tuned gains
//
// Usage: rabbit <mode> <results directory>

#include "hopper.hpp"

#include <heuresis/multiobjective/multiobjective.hpp>
#include <heuresis/search/search.hpp>
#include <heuresis/strategies/strategies.hpp>
#include <heuresis/swarms/swarms.hpp>
#include <robotics/legged/simulation.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numbers>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace rl = robotics::legged;
namespace hs = heuresis::search;
namespace hv = heuresis::variation;
namespace hst = heuresis::strategies;
namespace hsw = heuresis::swarms;
namespace hm = heuresis::multiobjective;
using experiment::Hopper;

constexpr double MOON = 1.62;  // m/s², the robotics domain's default site (legged/simulation.hpp)

// ---- the gains searched, on the unit box --------------------------------------------------------------------------

struct Parameter {
  const char* name;
  const char* unit;
  double lower, upper;
  double rl::Gains::*field;
};

// Ten of the controller's gains (legged/controller.hpp, Gains), each over a range about its default.
const std::vector<Parameter> PARAMETERS = {
    {"leg_stiffness", "N/m", 200, 800, &rl::Gains::leg_stiffness},
    {"leg_damping", "N s/m", 0, 10, &rl::Gains::leg_damping},
    {"leg_rest", "m", 0.17, 0.215, &rl::Gains::leg_rest},
    {"hop_ankle", "rad", 55 * rl::DEG, 85 * rl::DEG, &rl::Gains::hop_ankle},
    {"attitude_kp", "N m/rad", 0.2, 3.0, &rl::Gains::attitude_kp},
    {"momentum_rate", "1/s", 5, 80, &rl::Gains::momentum_rate},
    {"placement", "s", 0.02, 0.20, &rl::Gains::placement},
    {"speed_integral", "", 0.0, 1.0, &rl::Gains::speed_integral},
    {"height_gain", "N/m", 5, 60, &rl::Gains::height_gain},
    {"thrust_max", "N", 4, 20, &rl::Gains::thrust_max},
};
constexpr std::size_t D = 10;

rl::Gains gains_of(const std::vector<double>& x) {
  rl::Gains g;
  for (std::size_t i = 0; i < D; i++) g.*(PARAMETERS[i].field) = PARAMETERS[i].lower + std::clamp(x[i], 0.0, 1.0) * (PARAMETERS[i].upper - PARAMETERS[i].lower);
  return g;
}

std::vector<double> unit_of(const rl::Gains& g) {
  std::vector<double> x(D);
  for (std::size_t i = 0; i < D; i++) x[i] = (g.*(PARAMETERS[i].field) - PARAMETERS[i].lower) / (PARAMETERS[i].upper - PARAMETERS[i].lower);
  return x;
}

// ---- a trial: stand, then hop as asked; measured over the last seconds ----------------------------------------------

struct Task {
  double speed = 0.4, height = 0.12;  // asked: m/s, m
  double stand = 0.5, hop = 8.0, window = 4.0;  // s
};

struct Trial {
  bool failed = false;
  std::string reason;
  double speed = 0, cot = 0, power = 0, copper = 0;  // over the window: m/s, 1, W, W
  double apex = 0, apex_sd = 0, pitch = 0;           // m, m, rad (largest over the hop)
  int landings = 0;                                  // in the window
  double simulated = 0;                              // s
};

Trial trial(const rl::Rabbit& r, const rl::Gains& g, const Task& task, std::vector<experiment::Pose>* poses = nullptr, bool feel = false) {
  Trial out;
  Hopper h(r, g, MOON, feel);
  if (poses) h.keep_poses(1.0 / 30);
  h.advance(task.stand);
  h.command({rl::Behaviour::hop, task.speed, task.height});
  if (!h.controller().refusal().empty()) {
    out.failed = true, out.reason = "refused";
    return out;
  }
  h.reset_max_pitch();
  const double t_window = task.stand + task.hop - task.window, t_end = task.stand + task.hop;
  double x0 = 0, e0 = 0, c0 = 0;
  std::size_t a0 = 0;
  bool in_window = false;
  while (h.time() < t_end - 1e-9) {
    const double chunk = std::min(0.1, (in_window ? t_end : t_window) - h.time());
    h.advance(chunk);
    if (!h.finite()) return out.failed = true, out.reason = "diverged", out.simulated = h.time(), out;
    if (h.fallen()) return out.failed = true, out.reason = "fell", out.simulated = h.time(), out;
    if (h.max_pitch() > 60 * rl::DEG) return out.failed = true, out.reason = "tipped over", out.simulated = h.time(), out;
    if (!in_window && h.time() >= t_window - 1e-9) {
      in_window = true, x0 = h.q()[rl::X], e0 = h.electrical_energy(), c0 = h.copper_energy(), a0 = h.apexes().size();
    }
  }
  out.simulated = h.time();
  const double dx = h.q()[rl::X] - x0, energy = h.electrical_energy() - e0;
  out.speed = dx / task.window;
  out.power = energy / task.window;
  out.copper = (h.copper_energy() - c0) / task.window;
  out.cot = std::fabs(dx) > 0.05 ? energy / (h.mass() * MOON * std::fabs(dx)) : std::numeric_limits<double>::infinity();
  out.pitch = h.max_pitch();
  const auto& ap = h.apexes();
  out.landings = int(ap.size() - a0);
  if (out.landings < 3) return out.failed = true, out.reason = "stopped hopping", out;
  double s = 0, s2 = 0;
  for (std::size_t k = a0; k < ap.size(); k++) s += ap[k].second, s2 += ap[k].second * ap[k].second;
  out.apex = s / out.landings;
  out.apex_sd = std::sqrt(std::max(s2 / out.landings - out.apex * out.apex, 0.0));
  if (poses) *poses = h.poses();
  return out;
}

// ---- the objectives -------------------------------------------------------------------------------------------------

// One task: hop forward at the speed and height asked, at the least cost of transport. Constraints (satisfied at ≤ 0):
// the speed within 3 cm/s of the one asked, the apex within 2 cm, the pitch within 20°. A fall, a refusal, a
// divergence or a stop is a failure with its reason.
hs::Evaluation single(const rl::Rabbit& r, const std::vector<double>& x, const Task& task) {
  const Trial t = trial(r, gains_of(x), task);
  if (t.failed) return hs::failure(t.reason, t.simulated);
  return hs::evaluation({t.cot}, {std::fabs(t.speed - task.speed) / 0.03 - 1, std::fabs(t.apex - task.height) / 0.02 - 1, t.pitch / (20 * rl::DEG) - 1}, {}, t.simulated);
}

// Two objectives: the speed reached (maximised, so its negative minimised) and the mean electrical power drawn (the
// battery's drain, which sets the endurance). The cost of transport is not the second objective: copper loss is most
// of the power and nearly constant, so the cost per metre falls with speed up to the gait's limit, and a front of speed
// against it collapses to a point. The speed asked is the eleventh coordinate, from 0.1 to 1.0 m/s; the apex asked
// stays 0.12 m. Constraints: the apex within 2 cm of the one asked, the pitch within 20°.
constexpr double SPEED_LO = 0.1, SPEED_HI = 1.0;
double asked_speed(const std::vector<double>& x) { return SPEED_LO + std::clamp(x[D], 0.0, 1.0) * (SPEED_HI - SPEED_LO); }
hs::Evaluation two(const rl::Rabbit& r, const std::vector<double>& x) {
  Task task;
  task.speed = asked_speed(x);
  const Trial t = trial(r, gains_of(std::vector<double>(x.begin(), x.begin() + D)), task);
  if (t.failed) return hs::failure(t.reason, t.simulated);
  return hs::evaluation({-t.speed, t.power}, {std::fabs(t.apex - task.height) / 0.02 - 1, t.pitch / (20 * rl::DEG) - 1}, {}, t.simulated);
}

// ---- output -------------------------------------------------------------------------------------------------------

std::string num(double v) {
  if (!std::isfinite(v)) return "null";
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.10g", v);
  return buf;
}

std::string json_string(const std::string& s) {
  std::string o = "\"";
  for (char c : s) {
    if (c == '"' || c == '\\') o += '\\';
    o += c;
  }
  return o + "\"";
}

std::string vec(const std::vector<double>& x) {
  std::string s = "[";
  for (std::size_t i = 0; i < x.size(); i++) s += (i ? ", " : "") + num(x[i]);
  return s + "]";
}

std::string trial_json(const Trial& t) {
  return "{\"failed\": " + std::string(t.failed ? "true" : "false") + ", \"reason\": " + json_string(t.reason) + ", \"speed\": " + num(t.speed) + ", \"cot\": " + num(t.cot) +
         ", \"power\": " + num(t.power) + ", \"copper\": " + num(t.copper) + ", \"apex\": " + num(t.apex) + ", \"apex_sd\": " + num(t.apex_sd) + ", \"pitch_deg\": " +
         num(t.pitch / rl::DEG) + ", \"landings\": " + std::to_string(t.landings) + "}";
}

// Runs jobs on the machine's threads.
void parallel(std::size_t n, const std::function<void(std::size_t)>& job, unsigned threads = 8) {
  std::atomic<std::size_t> next{0};
  std::vector<std::thread> pool;
  for (unsigned w = 0; w < std::max(1u, std::min(threads, std::thread::hardware_concurrency())); w++)
    pool.emplace_back([&] {
      for (std::size_t i; (i = next++) < n;) job(i);
    });
  for (auto& th : pool) th.join();
}

double seconds_since(std::chrono::steady_clock::time_point t0) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }

// ---- modes ----------------------------------------------------------------------------------------------------------

int check(const std::filesystem::path& dir) {
  const rl::Rabbit r = rl::rabbit();
  rl::SimulationConfig c;
  c.method = rl::Method::semi_implicit;
  rl::Simulation sim(r, c);
  Hopper h(r, rl::Gains{}, MOON);
  sim.advance(1.0), h.advance(1.0);
  sim.command({rl::Behaviour::hop, 0.3, 0.15}), h.command({rl::Behaviour::hop, 0.3, 0.15});
  sim.advance(9.0), h.advance(9.0);
  double worst = 0;
  for (std::size_t k = 0; k < rl::DOF; k++) worst = std::max({worst, std::fabs(sim.q()[k] - h.q()[k]), std::fabs(sim.v()[k] - h.v()[k])});
  const bool same = std::memcmp(sim.q().data(), h.q().data(), sizeof(double) * rl::DOF) == 0 && std::memcmp(sim.v().data(), h.v().data(), sizeof(double) * rl::DOF) == 0;
  const double e_sim = sim.electrical_energy(), e_loop = h.electrical_energy();
  std::fprintf(stderr, "check: after 10 s the loop's state %s the Simulation's (largest difference %.3e); electrical energy %.9f J against %.9f J; %d hops\n",
               same ? "equals, bit for bit," : "differs from", worst, e_loop, e_sim, sim.controller().hops());
  // Without the feet's readings: the same motion while the friction it asks (0.5) stays below 3/4 of what it believes.
  Hopper n(r, rl::Gains{}, MOON, false);
  n.advance(1.0), n.command({rl::Behaviour::hop, 0.3, 0.15}), n.advance(9.0);
  const bool same_unfelt = std::memcmp(sim.q().data(), n.q().data(), sizeof(double) * rl::DOF) == 0;
  double worst_unfelt = 0;
  for (std::size_t k = 0; k < rl::DOF; k++) worst_unfelt = std::max(worst_unfelt, std::fabs(sim.q()[k] - n.q()[k]));
  std::fprintf(stderr, "check: without the feet's readings the state %s (largest difference %.3e)\n", same_unfelt ? "is the same, bit for bit" : "differs", worst_unfelt);
  // The cost of a trial, with and without the feet's readings.
  auto t0 = std::chrono::steady_clock::now();
  const Trial t = trial(r, rl::Gains{}, Task{}, nullptr, true);
  const double cost = seconds_since(t0);
  t0 = std::chrono::steady_clock::now();
  const Trial u = trial(r, rl::Gains{}, Task{}, nullptr, false);
  const double cost_unfelt = seconds_since(t0);
  std::fprintf(stderr, "a trial of %.1f s simulated costs %.3f s here (%.4f s a second); without the feet's readings %.3f s (%.4f s a second)\n", t.simulated, cost, cost / t.simulated,
               cost_unfelt, cost_unfelt / u.simulated);
  std::ofstream f(dir / "check.json");
  f << "{\"bitwise_equal\": " << (same ? "true" : "false") << ", \"largest_difference\": " << num(worst) << ", \"electrical_energy_loop\": " << num(e_loop)
    << ", \"electrical_energy_simulation\": " << num(e_sim) << ", \"hops\": " << sim.controller().hops() << ", \"unfelt_bitwise_equal\": " << (same_unfelt ? "true" : "false")
    << ", \"unfelt_largest_difference\": " << num(worst_unfelt) << ", \"trial_seconds\": " << num(cost) << ", \"trial_seconds_unfelt\": " << num(cost_unfelt)
    << ", \"trial_simulated\": " << num(t.simulated) << "}\n";
  return same ? 0 : 1;
}

int explore(const std::filesystem::path& dir) {
  const rl::Rabbit r = rl::rabbit();
  std::vector<std::pair<double, double>> grid;
  for (double h : {0.08, 0.12, 0.16})
    for (double v : {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}) grid.push_back({v, h});
  std::vector<Trial> out(grid.size());
  const auto t0 = std::chrono::steady_clock::now();
  parallel(grid.size(), [&](std::size_t i) {
    Task task;
    task.speed = grid[i].first, task.height = grid[i].second;
    out[i] = trial(r, rl::Gains{}, task);
  });
  std::fprintf(stderr, "explore: %zu trials in %.1f s\n", grid.size(), seconds_since(t0));
  std::ofstream f(dir / "explore.json");
  f << "[\n";
  for (std::size_t i = 0; i < grid.size(); i++) {
    const auto& t = out[i];
    std::fprintf(stderr, "asked %.2f m/s, %.2f m: %s%s speed %.3f  CoT %.2f  power %.2f W (copper %.2f W)  apex %.3f ± %.3f  pitch %.1f°  landings %d\n", grid[i].first, grid[i].second,
                 t.failed ? "FAILED " : "", t.reason.c_str(), t.speed, t.cot, t.power, t.copper, t.apex, t.apex_sd, t.pitch / rl::DEG, t.landings);
    f << "  {\"asked_speed\": " << num(grid[i].first) << ", \"asked_height\": " << num(grid[i].second) << ", \"trial\": " << trial_json(t) << "}" << (i + 1 < grid.size() ? "," : "") << "\n";
  }
  f << "]\n";
  return 0;
}

// ---- tune: five searches on one task ----------------------------------------------------------------------------------

// The task: hop forward at 0.6 m/s with an apex of 0.12 m, where the default gains pitch to 55°.
Task tune_task() {
  Task t;
  t.speed = 0.6, t.height = 0.12;
  return t;
}

hv::Box unit_box(std::size_t n) { return hv::Box{std::vector<double>(n, 0.0), std::vector<double>(n, 1.0)}; }

struct Checkpoint {
  std::size_t evaluations;
  double best_feasible, best_violation;  // +∞ until met
  std::size_t failures;
};

struct RunRecord {
  std::string algorithm;
  std::uint64_t seed = 0;
  std::vector<Checkpoint> checkpoints;
  std::map<std::string, std::size_t> reasons;
  std::size_t evaluations = 0, feasible = 0;
  std::vector<double> best_x;
  hs::Evaluation best;
  std::string ending;
  double seconds = 0;
};

template <class A>
RunRecord drive(A& algorithm, const std::string& name, std::uint64_t seed, double budget, const rl::Rabbit& r, const Task& task) {
  RunRecord rec;
  rec.algorithm = name, rec.seed = seed;
  double best_f = std::numeric_limits<double>::infinity(), best_v = best_f;
  std::size_t failures = 0;
  auto f = [&](const std::vector<double>& x) {
    auto e = single(r, x, task);
    rec.evaluations++;
    if (e.failed()) failures++, rec.reasons[e.reason]++;
    else {
      best_v = std::min(best_v, e.violation());
      if (e.violation() == 0) rec.feasible++, best_f = std::min(best_f, e.objectives[0]);
    }
    if (rec.evaluations % 10 == 0) rec.checkpoints.push_back({rec.evaluations, best_f, best_v, failures});
    return e;
  };
  const auto t0 = std::chrono::steady_clock::now();
  const auto out = hs::run(algorithm, f, hs::Budget{budget});
  rec.seconds = seconds_since(t0);
  if (out.best) rec.best_x = out.best->genotype, rec.best = out.best->evaluation;
  rec.ending = out.ending == hs::Ending::budget ? "budget" : out.ending == hs::Ending::stopped ? "stopped" : out.ending == hs::Ending::stalled ? "stalled" : "target";
  return rec;
}

const std::vector<std::string> ALGORITHMS = {"random", "(1+1)-ES", "CMA-ES", "DE", "PSO"};

RunRecord tune_once(const std::string& name, std::uint64_t seed, double budget, const rl::Rabbit& r, const Task& task) {
  const hv::Box box = unit_box(D);
  const hs::Draws draws(seed);
  // A start for the strategies that need one: a uniform point of the box, the seed's own.
  calx::random::Stream s = draws(hs::tag("experiment.start"), 0, 0);
  const std::vector<double> x0 = hv::uniform_point(box, s);
  if (name == "random") {
    hs::RandomSearch<hv::Real> a([box](calx::random::Stream& st) { return hv::uniform_point(box, st); }, draws);
    return drive(a, name, seed, budget, r, task);
  }
  if (name == "(1+1)-ES") {
    hst::OnePlusOneOptions o;
    o.bounds = box;
    hst::OnePlusOneEs a(x0, 0.25, o, draws);
    return drive(a, name, seed, budget, r, task);
  }
  if (name == "CMA-ES") {
    hst::CmaOptions o;
    o.bounds = box;
    hst::Cma a(x0, 0.25, o, draws);
    return drive(a, name, seed, budget, r, task);
  }
  if (name == "DE") {
    hst::DeOptions o;
    o.population = 20;
    hst::DifferentialEvolution a(box, o, draws);
    return drive(a, name, seed, budget, r, task);
  }
  hsw::PsoOptions o;
  o.swarm = 20, o.rule = hsw::PsoRule::inertia, o.topology = hsw::PsoTopology::ring;
  hsw::ParticleSwarm a(box, o, draws);
  return drive(a, name, seed, budget, r, task);
}

std::string poses_json(const std::vector<experiment::Pose>& ps) {
  std::ostringstream o;
  o << "[";
  for (std::size_t i = 0; i < ps.size(); i++) {
    o << (i ? ",\n" : "\n") << "{\"t\": " << num(ps[i].t) << ", \"phase\": " << json_string(rl::phase_name(ps[i].phase)) << ", \"p\": [";
    for (std::size_t k = 0; k < ps[i].p.size(); k++) o << (k ? ", " : "") << "[" << num(ps[i].p[k][0]) << ", " << num(ps[i].p[k][1]) << "]";
    o << "]}";
  }
  o << "]";
  return o.str();
}

int tune(const std::filesystem::path& dir, int seeds, double budget) {
  const rl::Rabbit r = rl::rabbit();
  const Task task = tune_task();
  std::vector<std::pair<std::string, std::uint64_t>> jobs;
  for (int s = 0; s < seeds; s++)
    for (const auto& a : ALGORITHMS) jobs.push_back({a, std::uint64_t(2000 + s)});
  std::vector<RunRecord> out(jobs.size());
  const auto t0 = std::chrono::steady_clock::now();
  std::mutex m;
  parallel(jobs.size(), [&](std::size_t i) {
    out[i] = tune_once(jobs[i].first, jobs[i].second, budget, r, task);
    std::lock_guard<std::mutex> lock(m);
    std::fprintf(stderr, "%-9s seed %llu: best %s %.4f (violation %.3f), %zu feasible of %zu, %s, %.0f s\n", out[i].algorithm.c_str(), (unsigned long long)out[i].seed,
                 out[i].best.failed() ? "failed" : "", out[i].best.value(), out[i].best.failed() ? 0.0 : out[i].best.violation(), out[i].feasible, out[i].evaluations,
                 out[i].ending.c_str(), out[i].seconds);
  });
  std::fprintf(stderr, "tune: %zu runs in %.1f s\n", jobs.size(), seconds_since(t0));
  // The default gains on the task, and the best gains found, measured again, with their motions.
  std::vector<experiment::Pose> p_default, p_best;
  const Trial t_default = trial(r, rl::Gains{}, task, &p_default);
  std::size_t best = 0;
  for (std::size_t i = 1; i < out.size(); i++)
    if (hs::better(out[i].best, out[best].best)) best = i;
  const Trial t_best = trial(r, gains_of(out[best].best_x), task, &p_best);
  std::ofstream f(dir / "tune.json");
  f << "{\n  \"task\": {\"speed\": " << num(task.speed) << ", \"height\": " << num(task.height) << ", \"stand\": " << num(task.stand) << ", \"hop\": " << num(task.hop)
    << ", \"window\": " << num(task.window) << "},\n  \"budget\": " << num(budget) << ",\n  \"seeds\": " << seeds << ",\n  \"parameters\": [";
  for (std::size_t i = 0; i < D; i++)
    f << (i ? ", " : "") << "{\"name\": " << json_string(PARAMETERS[i].name) << ", \"unit\": " << json_string(PARAMETERS[i].unit) << ", \"lower\": " << num(PARAMETERS[i].lower)
      << ", \"upper\": " << num(PARAMETERS[i].upper) << "}";
  f << "],\n  \"default_x\": " << vec(unit_of(rl::Gains{})) << ",\n  \"default\": " << trial_json(t_default) << ",\n  \"best_run\": " << best << ",\n  \"best\": " << trial_json(t_best)
    << ",\n  \"runs\": [\n";
  for (std::size_t i = 0; i < out.size(); i++) {
    const auto& rec = out[i];
    f << "    {\"algorithm\": " << json_string(rec.algorithm) << ", \"seed\": " << rec.seed << ", \"evaluations\": " << rec.evaluations << ", \"feasible\": " << rec.feasible
      << ", \"ending\": " << json_string(rec.ending) << ", \"seconds\": " << num(rec.seconds) << ", \"best_value\": " << num(rec.best.value())
      << ", \"best_failed\": " << (rec.best.failed() ? "true" : "false") << ", \"best_violation\": " << num(rec.best.failed() ? NAN : rec.best.violation()) << ", \"best_x\": " << vec(rec.best_x)
      << ", \"reasons\": {";
    bool first = true;
    for (const auto& [k, v] : rec.reasons) f << (first ? "" : ", ") << json_string(k) << ": " << v, first = false;
    f << "}, \"checkpoints\": [";
    for (std::size_t k = 0; k < rec.checkpoints.size(); k++) {
      const auto& c = rec.checkpoints[k];
      f << (k ? ", " : "") << "[" << c.evaluations << ", " << num(c.best_feasible) << ", " << num(c.best_violation) << ", " << c.failures << "]";
    }
    f << "]}" << (i + 1 < out.size() ? "," : "") << "\n";
  }
  f << "  ]\n}\n";
  std::ofstream g(dir / "tune_motions.json");
  g << "{\"default\": " << poses_json(p_default) << ",\n\"best\": " << poses_json(p_best) << "}\n";
  std::fprintf(stderr, "default gains on the task: %s %s speed %.3f CoT %.3f pitch %.1f°; best (%s, seed %llu): speed %.3f CoT %.3f power %.2f W apex %.3f pitch %.1f°\n",
               t_default.failed ? "FAILED" : "", t_default.reason.c_str(), t_default.speed, t_default.cot, t_default.pitch / rl::DEG, out[best].algorithm.c_str(),
               (unsigned long long)out[best].seed, t_best.speed, t_best.cot, t_best.power, t_best.apex, t_best.pitch / rl::DEG);
  return 0;
}

// ---- front: two objectives ---------------------------------------------------------------------------------------------

struct FrontRecord {
  std::string algorithm;
  std::uint64_t seed = 0;
  std::vector<std::pair<std::vector<double>, hs::Evaluation>> feasible;  // every feasible point met, with its genotype
  std::vector<std::pair<std::size_t, double>> hypervolume;               // after n evaluations, of the feasible points met
  std::map<std::string, std::size_t> reasons;
  std::size_t evaluations = 0;
  double seconds = 0;
};

// The hypervolume's reference point: speed 0, power 6 W.
const hm::Point REFERENCE = {0.0, 6.0};

template <class A>
FrontRecord drive_front(A& algorithm, const std::string& name, std::uint64_t seed, double budget, const rl::Rabbit& r) {
  FrontRecord rec;
  rec.algorithm = name, rec.seed = seed;
  std::vector<hm::Point> pts;
  auto f = [&](const std::vector<double>& x) {
    auto e = two(r, x);
    rec.evaluations++;
    if (e.failed()) rec.reasons[e.reason]++;
    else if (e.violation() == 0) rec.feasible.push_back({x, e}), pts.push_back(e.objectives);
    if (rec.evaluations % 50 == 0) rec.hypervolume.push_back({rec.evaluations, pts.empty() ? 0.0 : hm::hypervolume(pts, REFERENCE)});
    return e;
  };
  const auto t0 = std::chrono::steady_clock::now();
  hs::run(algorithm, f, hs::Budget{budget});
  rec.seconds = seconds_since(t0);
  return rec;
}

int front(const std::filesystem::path& dir, int seeds, double budget) {
  const rl::Rabbit r = rl::rabbit();
  std::vector<std::pair<std::string, std::uint64_t>> jobs;
  for (int s = 0; s < seeds; s++)
    for (const char* a : {"NSGA-II", "SMS-EMOA"}) jobs.push_back({a, std::uint64_t(3000 + s)});
  std::vector<FrontRecord> out(jobs.size());
  const auto t0 = std::chrono::steady_clock::now();
  std::mutex m;
  parallel(jobs.size(), [&](std::size_t i) {
    const hs::Draws draws(jobs[i].second);
    const auto kernels = hm::real_kernels(unit_box(D + 1), 15, 20);
    if (jobs[i].first == "NSGA-II") {
      hm::Nsga2<hv::Real> a(kernels, 40, draws);
      out[i] = drive_front(a, jobs[i].first, jobs[i].second, budget, r);
    } else {
      hm::SmsEmoa<hv::Real> a(kernels, 40, draws);
      out[i] = drive_front(a, jobs[i].first, jobs[i].second, budget, r);
    }
    std::lock_guard<std::mutex> lock(m);
    std::fprintf(stderr, "%-8s seed %llu: %zu feasible of %zu, hypervolume %.4f, %.0f s\n", out[i].algorithm.c_str(), (unsigned long long)out[i].seed, out[i].feasible.size(),
                 out[i].evaluations, out[i].hypervolume.empty() ? 0.0 : out[i].hypervolume.back().second, out[i].seconds);
  });
  std::fprintf(stderr, "front: %zu runs in %.1f s\n", jobs.size(), seconds_since(t0));

  // The union of every run's feasible points, and its non-dominated set.
  std::vector<std::pair<std::vector<double>, hs::Evaluation>> all;
  std::vector<std::string> from;
  for (const auto& rec : out)
    for (const auto& p : rec.feasible) all.push_back(p), from.push_back(rec.algorithm);
  std::vector<hm::Point> pts;
  for (const auto& p : all) pts.push_back(p.second.objectives);
  auto nd = hm::non_dominated(pts);
  std::sort(nd.begin(), nd.end(), [&](std::size_t a, std::size_t b) { return pts[a][0] > pts[b][0]; });  // slowest first
  // Each run's own non-dominated set.
  std::ofstream f(dir / "front.json");
  f << "{\n  \"budget\": " << num(budget) << ", \"seeds\": " << seeds << ", \"population\": 40, \"reference\": " << vec(REFERENCE) << ", \"speed_asked\": [" << num(SPEED_LO) << ", "
    << num(SPEED_HI) << "],\n  \"union_front\": [\n";
  for (std::size_t k = 0; k < nd.size(); k++) {
    const auto& [x, e] = all[nd[k]];
    f << "    {\"speed\": " << num(-e.objectives[0]) << ", \"power\": " << num(e.objectives[1]) << ", \"asked\": " << num(asked_speed(x)) << ", \"from\": " << json_string(from[nd[k]])
      << ", \"x\": " << vec(x) << "}" << (k + 1 < nd.size() ? "," : "") << "\n";
  }
  f << "  ],\n  \"runs\": [\n";
  for (std::size_t i = 0; i < out.size(); i++) {
    const auto& rec = out[i];
    std::vector<hm::Point> own;
    for (const auto& p : rec.feasible) own.push_back(p.second.objectives);
    f << "    {\"algorithm\": " << json_string(rec.algorithm) << ", \"seed\": " << rec.seed << ", \"evaluations\": " << rec.evaluations << ", \"feasible\": " << rec.feasible.size()
      << ", \"seconds\": " << num(rec.seconds) << ", \"reasons\": {";
    bool first = true;
    for (const auto& [k, v] : rec.reasons) f << (first ? "" : ", ") << json_string(k) << ": " << v, first = false;
    f << "}, \"hypervolume\": [";
    for (std::size_t k = 0; k < rec.hypervolume.size(); k++) f << (k ? ", " : "") << "[" << rec.hypervolume[k].first << ", " << num(rec.hypervolume[k].second) << "]";
    f << "], \"front\": [";
    bool first_point = true;
    for (std::size_t k : hm::non_dominated(own)) f << (first_point ? "" : ", ") << "[" << num(-own[k][0]) << ", " << num(own[k][1]) << "]", first_point = false;
    f << "]}" << (i + 1 < out.size() ? "," : "") << "\n";
  }
  f << "  ]\n}\n";

  // The motions of chosen points of the union front: the slowest, the fastest and two between, by speed.
  std::vector<std::size_t> chosen;
  if (!nd.empty()) {
    const double lo = -pts[nd.front()][0], hi = -pts[nd.back()][0];
    for (double q : {0.0, 1.0 / 3, 2.0 / 3, 1.0}) {
      const double want = lo + q * (hi - lo);
      std::size_t best = nd.front();
      for (std::size_t k : nd)
        if (std::fabs(-pts[k][0] - want) < std::fabs(-pts[best][0] - want)) best = k;
      if (std::find(chosen.begin(), chosen.end(), best) == chosen.end()) chosen.push_back(best);
    }
  }
  std::ofstream g(dir / "front_motions.json");
  g << "[\n";
  for (std::size_t c = 0; c < chosen.size(); c++) {
    const auto& x = all[chosen[c]].first;
    Task task;
    task.speed = asked_speed(x);
    std::vector<experiment::Pose> poses;
    const Trial t = trial(r, gains_of(std::vector<double>(x.begin(), x.begin() + D)), task, &poses);
    g << "{\"asked\": " << num(task.speed) << ", \"trial\": " << trial_json(t) << ", \"x\": " << vec(x) << ", \"poses\": " << poses_json(poses) << "}" << (c + 1 < chosen.size() ? ",\n" : "\n");
    std::fprintf(stderr, "front point %zu: asked %.3f m/s, reached %.3f m/s at %.2f W (CoT %.2f), apex %.3f, pitch %.1f°\n", c, task.speed, t.speed, t.power, t.cot, t.apex, t.pitch / rl::DEG);
  }
  g << "]\n";
  return 0;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "usage: rabbit <check|explore|tune|front> <results directory> [seeds] [budget]\n";
    return 2;
  }
  const std::string mode = argv[1];
  const std::filesystem::path dir = argv[2];
  std::filesystem::create_directories(dir);
  if (mode == "check") return check(dir);
  if (mode == "explore") return explore(dir);
  if (mode == "tune") return tune(dir, argc > 3 ? std::atoi(argv[3]) : 10, argc > 4 ? std::atof(argv[4]) : 1500);
  if (mode == "front") return front(dir, argc > 3 ? std::atoi(argv[3]) : 3, argc > 4 ? std::atof(argv[4]) : 3000);
  std::cerr << "unknown mode " << mode << "\n";
  return 2;
}
