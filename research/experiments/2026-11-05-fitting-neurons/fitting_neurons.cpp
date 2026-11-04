// Fitting neurons, and the repertoire of their firing (README.md).
//
// The neuroscience domain's adaptive exponential integrate-and-fire neuron (neuroscience::snn::models::Adex;
// Brette and Gerstner 2005), simulated by the domain's network engine on its 0.1 ms grid, is
//
//   (a) fitted to a behaviour it produced itself at hidden parameters: four parameters of its adaptation and
//       reset, (a, τ_w, b, V_r), from the initial and steady firing rates at six currents (the f–I curves and
//       their adaptation), and from the spike times themselves, by CMA-ES, differential evolution, a particle
//       swarm and random search under one budget; with the error along two planes through the truth;
//   (b) illuminated by MAP-Elites: six parameters of the cell varied under one current, every cell of a grid
//       of (log firing rate, coefficient of variation of the intervals) holding the cell nearest a reference
//       neuron that fires so, against random search with the same budget; each elite named by its firing class.
//
// Every number of the record is written to <results>/ by this program.
//
//   fitting_neurons <results directory> [time|fit|landscape|repertoire|exemplars|refine|all]

#include <heuresis/diversity/diversity.hpp>
#include <heuresis/search/search.hpp>
#include <heuresis/strategies/strategies.hpp>
#include <heuresis/swarms/swarms.hpp>
#include <heuresis/variation/variation.hpp>
#include <neuroscience/snn/snn.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace sn = neuroscience::snn;
namespace u = physis::units::unit;
namespace hs = heuresis::search;
namespace hv = heuresis::variation;
using Real = hv::Real;
using Adex = sn::models::Adex;

namespace {

constexpr double H = 0.1e-3;  // the engine's grid, s
std::filesystem::path out_dir;
// The integrator: the model's own (Dormand and Prince's adaptive method), at a tolerance relative to each
// component's typical size (the engine's default is 1e-6).
namespace pd = physis::dynamics;
pd::Stepping stepping{pd::Method::exact, 1, {1e-6, 1e-6}};

double seconds_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// The spike times (s) of a population of AdEx neurons of parameters p, neuron k driven by the constant current
// currents[k] (A) for T seconds, on the domain's engine.
std::vector<std::vector<double>> respond(const Adex::Params& p, const std::vector<double>& currents, double T) {
  sn::Network net({H, 1, 1, stepping});
  auto& pop = net.add<Adex>("n", currents.size(), p);
  net.set_bias(pop, currents);
  auto& rec = net.record_spikes(pop);
  net.run(T);
  std::vector<std::vector<double>> out(currents.size());
  for (const auto& s : rec.spikes) out[s.gid - pop.first].push_back(static_cast<double>(s.time) * H);
  for (auto& t : out) std::sort(t.begin(), t.end());
  return out;
}

// ---- (a) fitting ---------------------------------------------------------------------------------------------

// The four fitted parameters, searched in the unit cube: a = 8 nS·x₀, τ_w = 20 ms·25^{x₁} (logarithmic, 20 to
// 500 ms), b = 200 pA·x₂, V_r = −70 mV + 25 mV·x₃.
struct Fitted {
  double a, tau_w, b, V_r;
};
Fitted decode(const Real& x) { return {8 * u::nS * x[0], 20 * u::ms * std::pow(25.0, x[1]), 200 * u::pA * x[2], -70 * u::mV + 25 * u::mV * x[3]}; }
Real encode(const Fitted& f) { return {f.a / (8 * u::nS), std::log(f.tau_w / (20 * u::ms)) / std::log(25.0), f.b / (200 * u::pA), (f.V_r + 70 * u::mV) / (25 * u::mV)}; }

// The cell: Naud et al.'s (2008, Table 1) adapting neuron, of which (a, τ_w, b, V_r) are fitted.
Adex::Params fitted_cell(const Fitted& f) {
  Adex::Params p;
  p.C_m = 200 * u::pF;
  p.g_L = 12 * u::nS;
  p.E_L = -70 * u::mV;
  p.V_T = -50 * u::mV;
  p.Delta_T = 2 * u::mV;
  p.V_peak = 0;
  p.t_ref = 0;
  p.a = f.a;
  p.tau_w = f.tau_w;
  p.b = f.b;
  p.V_reset = f.V_r;
  return p;
}
const Fitted truth{2 * u::nS, 300 * u::ms, 60 * u::pA, -58 * u::mV};  // Naud et al. 2008, "adaptation"
const std::vector<double> currents{300e-12, 400e-12, 500e-12, 600e-12, 700e-12, 800e-12};
constexpr double T_fit = 1.0;  // s

// The initial rate 1/ISI₁ and the steady rate (the count in [0.5, 1) s over 0.5 s) of one train.
struct Rates {
  double initial = 0, steady = 0;
};
Rates rates(const std::vector<double>& t) {
  Rates r;
  if (t.size() >= 2) r.initial = 1 / (t[1] - t[0]);
  r.steady = static_cast<double>(std::count_if(t.begin(), t.end(), [](double x) { return x >= 0.5 && x < 1.0; })) / 0.5;
  return r;
}

// Σ_ij exp(−|x_i − y_j|/τ): van Rossum's (2001) inner product of two trains filtered by a causal exponential
// of time constant τ, up to the factor τ/2.
double vr_product(const std::vector<double>& x, const std::vector<double>& y, double tau) {
  double s = 0;
  for (double a : x)
    for (double b : y) s += std::exp(-std::fabs(a - b) / tau);
  return s;
}

struct Target {
  std::vector<std::vector<double>> spikes;
  std::vector<Rates> rates;
  std::vector<double> norms;  // each target train's squared van Rossum norm
};
constexpr double tau_vr = 20e-3;  // s

Target make_target() {
  Target t;
  t.spikes = respond(fitted_cell(truth), currents, T_fit);
  for (const auto& s : t.spikes) {
    t.rates.push_back(rates(s));
    t.norms.push_back(vr_product(s, s, tau_vr));
  }
  return t;
}

// The two errors of a response against the target:
//   features: the mean over the 2K rates of log²((f + 1 Hz)/(f* + 1 Hz)), the f–I curves and their adaptation;
//   spikes: the mean over the K currents of D²(s, s*)/‖s*‖², van Rossum's distance with τ = 20 ms, relative.
struct Errors {
  double features = 0, spikes = 0;
};
Errors errors(const std::vector<std::vector<double>>& s, const Target& t) {
  Errors e;
  for (std::size_t k = 0; k < s.size(); k++) {
    const Rates r = rates(s[k]);
    const double li = std::log((r.initial + 1) / (t.rates[k].initial + 1)), ls = std::log((r.steady + 1) / (t.rates[k].steady + 1));
    e.features += 0.5 * (li * li + ls * ls);
    const double d2 = vr_product(s[k], s[k], tau_vr) + t.norms[k] - 2 * vr_product(s[k], t.spikes[k], tau_vr);
    e.spikes += std::max(0.0, d2) / t.norms[k];
  }
  e.features /= static_cast<double>(s.size());
  e.spikes /= static_cast<double>(s.size());
  return e;
}

const hv::Box unit4 = hv::Box::cube(4, 0, 1);

struct RunResult {
  std::string algorithm, objective;
  std::uint64_t seed = 0;
  std::vector<hs::Progress> progress;
  Real best;
  double best_value = 0, evaluations = 0, wall = 0;
};

template <class A, class F>
RunResult drive(A& algorithm, F&& f, double budget, std::size_t batch) {
  RunResult r;
  const auto t0 = std::chrono::steady_clock::now();
  const auto out = hs::run(algorithm, f, hs::Budget{budget}, hs::Serial{}, batch);
  r.wall = seconds_since(t0);
  r.progress = out.progress;
  r.best = out.best->genotype;
  r.best_value = out.best->evaluation.value();
  r.evaluations = out.evaluations;
  return r;
}

void fit(calx::scheduler::ThreadPool& pool, const Target& target, std::size_t seeds, double budget) {
  const std::vector<std::string> algorithms{"CMA-ES", "DE", "PSO", "random"};
  const std::vector<std::string> objectives{"features", "spikes"};
  struct Job {
    std::size_t algorithm, objective;
    std::uint64_t seed;
  };
  std::vector<Job> jobs;
  for (std::size_t o = 0; o < objectives.size(); o++)
    for (std::size_t a = 0; a < algorithms.size(); a++)
      for (std::uint64_t s = 1; s <= seeds; s++) jobs.push_back({a, o, s});
  std::vector<RunResult> results(jobs.size());
  const auto t0 = std::chrono::steady_clock::now();
  calx::scheduler::parallel_for(pool, jobs.size(), 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t j = b; j < e; j++) {
      const Job job = jobs[j];
      const bool spikes = job.objective == 1;
      auto f = [&target, spikes](const Real& x) -> hs::Evaluation {
        try {
          const Errors err = errors(respond(fitted_cell(decode(x)), currents, T_fit), target);
          return hs::scalar(spikes ? err.spikes : err.features);
        } catch (const std::exception& ex) {
          return hs::failure(ex.what());
        }
      };
      // Each algorithm's draws keyed by its seed and the objective, so that runs are independent.
      const hs::Draws draws(1000 * job.seed + 17 * job.objective + job.algorithm);
      RunResult r;
      if (job.algorithm == 0) {
        heuresis::strategies::CmaOptions o;
        o.bounds = unit4;
        heuresis::strategies::Cma alg(Real(4, 0.5), 0.3, o, draws);
        r = drive(alg, f, budget, hs::unlimited);
      } else if (job.algorithm == 1) {
        heuresis::strategies::DeOptions o;
        o.population = 20;
        heuresis::strategies::DifferentialEvolution alg(unit4, o, draws);
        r = drive(alg, f, budget, hs::unlimited);
      } else if (job.algorithm == 2) {
        heuresis::swarms::PsoOptions o;
        o.swarm = 20;
        o.rule = heuresis::swarms::PsoRule::inertia;
        o.topology = heuresis::swarms::PsoTopology::ring;
        heuresis::swarms::ParticleSwarm alg(unit4, o, draws);
        r = drive(alg, f, budget, hs::unlimited);
      } else {
        hs::RandomSearch<Real> alg([](calx::random::Stream& s) { return hv::uniform_point(unit4, s); }, draws);
        r = drive(alg, f, budget, 1);
      }
      r.algorithm = algorithms[job.algorithm];
      r.objective = objectives[job.objective];
      r.seed = job.seed;
      results[j] = std::move(r);
    }
  });
  const double wall = seconds_since(t0);

  FILE* fc = std::fopen((out_dir / "convergence.csv").string().c_str(), "w");
  std::fprintf(fc, "objective,algorithm,seed,evaluations,best\n");
  for (const auto& r : results)
    for (const auto& p : r.progress) std::fprintf(fc, "%s,%s,%llu,%.0f,%.6e\n", r.objective.c_str(), r.algorithm.c_str(), static_cast<unsigned long long>(r.seed), p.evaluations, p.value);
  std::fclose(fc);

  // The fits: each run's best point in physical units, and both errors there (recomputed).
  FILE* ff = std::fopen((out_dir / "fits.csv").string().c_str(), "w");
  std::fprintf(ff, "objective,algorithm,seed,evaluations,best,error_features,error_spikes,a_nS,tau_w_ms,b_pA,V_r_mV,wall_s\n");
  for (const auto& r : results) {
    const Fitted p = decode(r.best);
    const Errors e = errors(respond(fitted_cell(p), currents, T_fit), target);
    std::fprintf(ff, "%s,%s,%llu,%.0f,%.6e,%.6e,%.6e,%.5f,%.4f,%.4f,%.5f,%.2f\n", r.objective.c_str(), r.algorithm.c_str(), static_cast<unsigned long long>(r.seed), r.evaluations,
                 r.best_value, e.features, e.spikes, p.a / u::nS, p.tau_w / u::ms, p.b / u::pA, p.V_r / u::mV, r.wall);
  }
  std::fclose(ff);

  // The f–I curves of the target and of every run's best fit, for drawing.
  FILE* fi = std::fopen((out_dir / "fi_curves.csv").string().c_str(), "w");
  std::fprintf(fi, "who,objective,seed,current_pA,initial_Hz,steady_Hz\n");
  const std::vector<double> dense = [] {
    std::vector<double> I;
    for (int k = 0; k <= 30; k++) I.push_back((250 + 20 * k) * 1e-12);
    return I;
  }();
  {
    const auto s = respond(fitted_cell(truth), dense, T_fit);
    for (std::size_t k = 0; k < dense.size(); k++) {
      const Rates q = rates(s[k]);
      std::fprintf(fi, "truth,,0,%.0f,%.4f,%.4f\n", dense[k] / u::pA, q.initial, q.steady);
    }
  }
  for (const auto& r : results) {
    const auto s = respond(fitted_cell(decode(r.best)), dense, T_fit);
    for (std::size_t k = 0; k < dense.size(); k++) {
      const Rates q = rates(s[k]);
      std::fprintf(fi, "%s,%s,%llu,%.0f,%.4f,%.4f\n", r.algorithm.c_str(), r.objective.c_str(), static_cast<unsigned long long>(r.seed), dense[k] / u::pA, q.initial, q.steady);
    }
  }
  std::fclose(fi);
  std::printf("fit: %zu runs of %.0f evaluations in %.1f s\n", results.size(), budget, wall);
}

// ---- two stages: the rates first, then the spike times from there ------------------------------------------

// The fits of fits.csv: (objective, algorithm, seed) → the fitted parameters.
struct FitRow {
  std::string objective, algorithm;
  std::uint64_t seed;
  Fitted p;
};
std::vector<FitRow> read_fits() {
  std::vector<FitRow> rows;
  FILE* f = std::fopen((out_dir / "fits.csv").string().c_str(), "r");
  if (!f) throw std::runtime_error("refine needs fits.csv: run the fit first");
  char line[1024];
  if (!std::fgets(line, sizeof line, f)) throw std::runtime_error("fits.csv is empty");
  while (std::fgets(line, sizeof line, f)) {
    std::vector<std::string> c;
    std::string cell;
    for (const char* q = line; *q && *q != '\n'; q++) {
      if (*q == ',') c.push_back(cell), cell.clear();
      else cell += *q;
    }
    c.push_back(cell);
    if (c.size() < 11) continue;
    rows.push_back({c[0], c[1], std::stoull(c[2]), {std::stod(c[7]) * u::nS, std::stod(c[8]) * u::ms, std::stod(c[9]) * u::pA, std::stod(c[10]) * u::mV}});
  }
  std::fclose(f);
  return rows;
}

// CMA-ES on the spike times, started at each CMA-ES fit of the rates with a small step (σ = 0.05 of the box).
void refine(calx::scheduler::ThreadPool& pool, const Target& target, double budget) {
  const auto rows = read_fits();
  std::vector<FitRow> starts;
  for (const auto& r : rows)
    if (r.objective == "features" && r.algorithm == "CMA-ES") starts.push_back(r);
  std::vector<RunResult> results(starts.size());
  calx::scheduler::parallel_for(pool, starts.size(), 1, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t j = b; j < e; j++) {
      auto f = [&target](const Real& x) -> hs::Evaluation {
        try {
          return hs::scalar(errors(respond(fitted_cell(decode(x)), currents, T_fit), target).spikes);
        } catch (const std::exception& ex) {
          return hs::failure(ex.what());
        }
      };
      heuresis::strategies::CmaOptions o;
      o.bounds = unit4;
      heuresis::strategies::Cma alg(encode(starts[j].p), 0.05, o, hs::Draws(900 + starts[j].seed));
      results[j] = drive(alg, f, budget, hs::unlimited);
      results[j].seed = starts[j].seed;
    }
  });
  FILE* fr = std::fopen((out_dir / "refine.csv").string().c_str(), "w");
  std::fprintf(fr, "seed,stage,evaluations,error_features,error_spikes,a_nS,tau_w_ms,b_pA,V_r_mV\n");
  FILE* fc = std::fopen((out_dir / "refine_convergence.csv").string().c_str(), "w");
  std::fprintf(fc, "seed,evaluations,best\n");
  for (std::size_t j = 0; j < starts.size(); j++) {
    const Fitted p0 = starts[j].p, p1 = decode(results[j].best);
    const Errors e0 = errors(respond(fitted_cell(p0), currents, T_fit), target), e1 = errors(respond(fitted_cell(p1), currents, T_fit), target);
    std::fprintf(fr, "%llu,rates,0,%.6e,%.6e,%.5f,%.4f,%.4f,%.5f\n", static_cast<unsigned long long>(starts[j].seed), e0.features, e0.spikes, p0.a / u::nS, p0.tau_w / u::ms, p0.b / u::pA, p0.V_r / u::mV);
    std::fprintf(fr, "%llu,spike times,%.0f,%.6e,%.6e,%.5f,%.4f,%.4f,%.5f\n", static_cast<unsigned long long>(starts[j].seed), results[j].evaluations, e1.features, e1.spikes, p1.a / u::nS, p1.tau_w / u::ms,
                 p1.b / u::pA, p1.V_r / u::mV);
    for (const auto& q : results[j].progress) std::fprintf(fc, "%llu,%.0f,%.6e\n", static_cast<unsigned long long>(starts[j].seed), q.evaluations, q.value);
  }
  std::fclose(fr);
  std::fclose(fc);

  // The spike trains at the six currents: the truth; the CMA-ES fit of the rates of seed 3 (its a furthest from
  // the truth among the exact fits); that fit refined on the spike times; CMA-ES on the spike times alone, seed 1.
  FILE* ft = std::fopen((out_dir / "trains.csv").string().c_str(), "w");
  std::fprintf(ft, "who,current_pA,spike_s\n");
  auto write = [&](const char* who, const Fitted& p) {
    const auto s = respond(fitted_cell(p), currents, T_fit);
    for (std::size_t k = 0; k < currents.size(); k++)
      for (double t : s[k]) std::fprintf(ft, "%s,%.0f,%.4f\n", who, currents[k] / u::pA, t);
  };
  write("truth", truth);
  for (std::size_t j = 0; j < starts.size(); j++)
    if (starts[j].seed == 3) {
      write("rates fit", starts[j].p);
      write("two-stage fit", decode(results[j].best));
    }
  for (const auto& r : rows)
    if (r.objective == "spikes" && r.algorithm == "CMA-ES" && r.seed == 1) write("spike-time fit", r.p);
  std::fclose(ft);
  std::printf("refine: %zu runs\n", starts.size());
}

// ---- the error on two planes through the truth -------------------------------------------------------------

void landscape(calx::scheduler::ThreadPool& pool, const Target& target, std::size_t n) {
  struct Point {
    int plane;
    double x, y;  // unit coordinates of the plane's two parameters
    Errors e;
  };
  std::vector<Point> pts;
  for (int plane = 0; plane < 2; plane++)
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < n; j++) pts.push_back({plane, double(i) / double(n - 1), double(j) / double(n - 1), {}});
  const Real t = encode(truth);
  const auto t0 = std::chrono::steady_clock::now();
  calx::scheduler::parallel_for(pool, pts.size(), 8, [&](std::size_t b, std::size_t e, std::size_t) {
    for (std::size_t k = b; k < e; k++) {
      Real x = t;
      // plane 0: (a, b) with τ_w and V_r at the truth; plane 1: (τ_w, b) with a and V_r at the truth.
      if (pts[k].plane == 0) x[0] = pts[k].x, x[2] = pts[k].y;
      else x[1] = pts[k].x, x[2] = pts[k].y;
      pts[k].e = errors(respond(fitted_cell(decode(x)), currents, T_fit), target);
    }
  });
  FILE* f = std::fopen((out_dir / "landscape.csv").string().c_str(), "w");
  std::fprintf(f, "plane,x_unit,y_unit,a_nS,tau_w_ms,b_pA,error_features,error_spikes\n");
  for (const auto& p : pts) {
    Real x = t;
    if (p.plane == 0) x[0] = p.x, x[2] = p.y;
    else x[1] = p.x, x[2] = p.y;
    const Fitted q = decode(x);
    std::fprintf(f, "%s,%.5f,%.5f,%.5f,%.4f,%.4f,%.6e,%.6e\n", p.plane == 0 ? "a-b" : "tau_w-b", p.x, p.y, q.a / u::nS, q.tau_w / u::ms, q.b / u::pA, p.e.features, p.e.spikes);
  }
  std::fclose(f);
  std::printf("landscape: %zu points in %.1f s\n", pts.size(), seconds_since(t0));
}

// ---- (b) the repertoire of firing ---------------------------------------------------------------------------

// Six parameters in the unit cube: C_m = 100 + 200x₀ pF, g_L = 5 + 15x₁ nS, a = 10x₂ nS, τ_w = 10·50^{x₃} ms,
// b = 250x₄ pA, V_r = −70 + 30x₅ mV; E_L = −70 mV, V_T = −50 mV, Δ_T = 2 mV, V_peak = 0, under 500 pA for 1 s.
Adex::Params repertoire_cell(const Real& x) {
  Adex::Params p;
  p.C_m = (100 + 200 * x[0]) * u::pF;
  p.g_L = (5 + 15 * x[1]) * u::nS;
  p.a = 10 * x[2] * u::nS;
  p.tau_w = 10 * std::pow(50.0, x[3]) * u::ms;
  p.b = 250 * x[4] * u::pA;
  p.V_reset = (-70 + 30 * x[5]) * u::mV;
  p.E_L = -70 * u::mV;
  p.V_T = -50 * u::mV;
  p.Delta_T = 2 * u::mV;
  p.V_peak = 0;
  p.t_ref = 0;
  return p;
}
// The reference: Naud et al.'s tonic neuron (C 200 pF, g_L 10 nS, a 2 nS, τ_w 30 ms, b 0, V_r −58 mV).
const Real reference{0.5, 1.0 / 3.0, 0.2, std::log(3.0) / std::log(50.0), 0.0, 0.4};
constexpr double I_rep = 500e-12, T_rep = 1.0;
const hv::Box unit6 = hv::Box::cube(6, 0, 1);
// The descriptor box: log₁₀ of the rate in [0.3, 2.5] (2 to 316 Hz), the intervals' CV in [0, 1.5].
const std::vector<double> d_lo{0.3, 0.0}, d_hi{2.5, 1.5};
const std::vector<std::size_t> d_res{44, 30};

double cv_of(const std::vector<double>& d) {
  if (d.size() < 2) return 0;
  const double m = std::accumulate(d.begin(), d.end(), 0.0) / double(d.size());
  double v = 0;
  for (double x : d) v += (x - m) * (x - m);
  return std::sqrt(v / double(d.size() - 1)) / m;
}

// The firing class of a train under a constant current, from its intervals d₁ … d_m (m ≥ 2):
//   bursting:      intervals in two groups, short within bursts and long between them (long: > half the
//                  longest), the longest more than three times the shortest, at least two long intervals and
//                  no fewer short than long, repeatedly;
//   initial burst: the largest ratio of successive intervals, above 2.5, among the first three, and every
//                  later interval long;
//   adapting:      the last interval more than 1.5 times the first;
//   fast:          otherwise, at 100 Hz or more;
//   tonic:         otherwise.
// These follow the patterns of Naud et al. (2008, Fig. 4) and the domain's checks of them (snn/tests/check_models.cpp).
std::string classify(const std::vector<double>& t) {
  std::vector<double> d;
  for (std::size_t i = 1; i < t.size(); i++) d.push_back(t[i] - t[i - 1]);
  const std::size_t m = d.size();
  const double mx = *std::max_element(d.begin(), d.end()), mn = *std::min_element(d.begin(), d.end());
  std::size_t longs = 0;
  for (double x : d) longs += x > 0.5 * mx;
  const std::size_t shorts = m - longs;
  if (m >= 4 && mx > 3 * mn && longs >= 2 && shorts >= longs) return "bursting";
  if (m >= 4) {
    std::size_t at = 0;
    double jump = 0;
    for (std::size_t i = 0; i + 1 < m; i++)
      if (d[i + 1] / d[i] > jump) jump = d[i + 1] / d[i], at = i;
    bool later_long = true;
    for (std::size_t i = at + 1; i < m; i++) later_long = later_long && d[i] > 0.5 * mx;
    if (jump > 2.5 && at < 3 && later_long) return "initial burst";
  }
  if (d.back() > 1.5 * d.front()) return "adapting";
  if (double(t.size()) / T_rep >= 100) return "fast";
  return "tonic";
}

struct Behaviour {
  bool ok = false;
  std::string reason, kind;
  double log_rate = 0, cv = 0;
  std::size_t spikes = 0;
};
Behaviour behaviour(const std::vector<double>& t) {
  Behaviour b;
  b.spikes = t.size();
  if (t.size() < 3) { b.reason = "fewer than three spikes"; return b; }
  std::vector<double> d;
  for (std::size_t i = 1; i < t.size(); i++) d.push_back(t[i] - t[i - 1]);
  b.log_rate = std::log10(double(t.size()) / T_rep);
  b.cv = cv_of(d);
  if (b.log_rate < d_lo[0] || b.log_rate > d_hi[0] || b.cv > d_hi[1]) { b.reason = "outside the descriptor box"; return b; }
  b.kind = classify(t);
  b.ok = true;
  return b;
}

hs::Evaluation repertoire_evaluation(const Real& x) {
  try {
    const auto s = respond(repertoire_cell(x), {I_rep}, T_rep)[0];
    const Behaviour b = behaviour(s);
    if (!b.ok) return hs::failure(b.reason);
    return hs::evaluation({hv::euclidean(x, reference)}, {}, {b.log_rate, b.cv});
  } catch (const std::exception& ex) {
    return hs::failure(ex.what());
  }
}

void repertoire(calx::scheduler::ThreadPool& pool, std::size_t budget) {
  using Grid = heuresis::diversity::Grid;
  heuresis::diversity::Kernels<Real> k;
  k.initialise = [](calx::random::Stream& s) { return hv::uniform_point(unit6, s); };
  k.mutate = [](const Real& x, calx::random::Stream& s) {
    Real y = hv::gaussian_mutation(x, 0.05, s);
    hv::repair(y, unit6, hv::Repair::reflect, s);
    return y;
  };
  k.vary = [](const Real& x, const Real& y, calx::random::Stream& s) {
    Real c = hv::iso_line(x, y, 0.01, 0.2, s);
    hv::repair(c, unit6, hv::Repair::reflect, s);
    return c;
  };
  typename heuresis::diversity::MapElites<Real, Grid>::Options o;
  o.initial = 500;
  o.variation_probability = 0.5;
  o.breadth = 64;
  heuresis::diversity::MapElites<Real, Grid> me(k, Grid(d_lo, d_hi, d_res), o, hs::Draws(7));
  hs::RandomSearch<Real> rs([](calx::random::Stream& s) { return hv::uniform_point(unit6, s); }, hs::Draws(8));
  heuresis::diversity::Archive<Real, Grid> random_archive(Grid(d_lo, d_hi, d_res));
  hs::Pool executor(pool);
  const double offset = 2.5;  // above √6, the largest distance in the unit cube

  FILE* fc = std::fopen((out_dir / "coverage.csv").string().c_str(), "w");
  std::fprintf(fc, "algorithm,evaluations,coverage,qd_score,filled\n");
  std::map<std::string, std::size_t> failures_me, failures_rs;
  const auto t0 = std::chrono::steady_clock::now();
  std::size_t done = 0;
  while (done < budget) {
    auto cs = me.ask(std::min<std::size_t>(64 * 8, budget - done));
    auto es = executor(std::span<const hs::Candidate<Real>>(cs), repertoire_evaluation);
    std::vector<hs::Told> told;
    for (std::size_t i = 0; i < cs.size(); i++) {
      if (es[i].failed()) failures_me[es[i].reason]++;
      told.push_back({cs[i].ticket, es[i]});
    }
    me.tell(told);
    done += cs.size();
    std::fprintf(fc, "MAP-Elites,%zu,%.5f,%.4f,%zu\n", done, me.archive().coverage(), me.archive().qd_score(offset), me.archive().size());
  }
  const double wall_me = seconds_since(t0);
  done = 0;
  while (done < budget) {
    auto cs = rs.ask(std::min<std::size_t>(512, budget - done));
    auto es = executor(std::span<const hs::Candidate<Real>>(cs), repertoire_evaluation);
    std::vector<hs::Told> told;
    for (std::size_t i = 0; i < cs.size(); i++) {
      if (es[i].failed()) failures_rs[es[i].reason]++;
      told.push_back({cs[i].ticket, es[i]});
      random_archive.offer({cs[i].ticket, cs[i].genotype, es[i], cs[i].origin});
    }
    rs.tell(told);
    done += cs.size();
    std::fprintf(fc, "random,%zu,%.5f,%.4f,%zu\n", done, random_archive.coverage(), random_archive.qd_score(offset), random_archive.size());
  }
  std::fclose(fc);

  // The elites, each with its class.
  std::map<std::string, std::size_t> classes_me, classes_rs;
  FILE* fe = std::fopen((out_dir / "elites.csv").string().c_str(), "w");
  std::fprintf(fe, "algorithm,cell,log_rate,cv,distance,class,spikes,C_pF,g_L_nS,a_nS,tau_w_ms,b_pA,V_r_mV\n");
  auto dump = [&](const char* name, const auto& archive, std::map<std::string, std::size_t>& counts) {
    for (auto c : archive.filled()) {
      const auto& e = *archive.elite(c);
      const auto s = respond(repertoire_cell(e.genotype), {I_rep}, T_rep)[0];
      const Behaviour b = behaviour(s);
      counts[b.kind]++;
      const auto p = repertoire_cell(e.genotype);
      std::fprintf(fe, "%s,%zu,%.5f,%.5f,%.5f,%s,%zu,%.3f,%.4f,%.4f,%.3f,%.3f,%.4f\n", name, c, e.evaluation.descriptors[0], e.evaluation.descriptors[1], e.evaluation.value(), b.kind.c_str(),
                   b.spikes, p.C_m / u::pF, p.g_L / u::nS, p.a / u::nS, p.tau_w / u::ms, p.b / u::pA, p.V_reset / u::mV);
    }
  };
  dump("MAP-Elites", me.archive(), classes_me);
  dump("random", random_archive, classes_rs);
  std::fclose(fe);

  FILE* fj = std::fopen((out_dir / "repertoire.json").string().c_str(), "w");
  std::fprintf(fj, "{\n  \"budget\": %zu,\n  \"cells\": %zu,\n  \"wall_map_elites_s\": %.1f,\n", budget, d_res[0] * d_res[1], wall_me);
  std::fprintf(fj, "  \"map_elites\": {\"filled\": %zu, \"coverage\": %.5f, \"qd_score\": %.3f},\n", me.archive().size(), me.archive().coverage(), me.archive().qd_score(offset));
  std::fprintf(fj, "  \"random\": {\"filled\": %zu, \"coverage\": %.5f, \"qd_score\": %.3f},\n", random_archive.size(), random_archive.coverage(), random_archive.qd_score(offset));
  auto counts = [&](const char* key, const std::map<std::string, std::size_t>& m, bool last) {
    std::fprintf(fj, "  \"%s\": {", key);
    bool first = true;
    for (const auto& [k2, v] : m) std::fprintf(fj, "%s\"%s\": %zu", first ? "" : ", ", k2.c_str(), v), first = false;
    std::fprintf(fj, "}%s\n", last ? "" : ",");
  };
  counts("classes_map_elites", classes_me, false);
  counts("classes_random", classes_rs, false);
  counts("failures_map_elites", failures_me, false);
  counts("failures_random", failures_rs, true);
  std::fprintf(fj, "}\n");
  std::fclose(fj);
  std::printf("repertoire: MAP-Elites %zu cells (%.1f%%), random %zu cells (%.1f%%), %.1f s\n", me.archive().size(), 100 * me.archive().coverage(), random_archive.size(),
              100 * random_archive.coverage(), wall_me);
}

// ---- exemplars of the classes -------------------------------------------------------------------------------

// For each firing class, the MAP-Elites elite (elites.csv) whose descriptors lie nearest the class's median
// descriptors, each scaled by its box's width; its potential over the first 600 ms, every 0.3 ms, and its spikes.
void exemplars() {
  struct Row {
    std::string kind;
    double log_rate, cv;
    Adex::Params p;
  };
  std::vector<Row> rows;
  FILE* f = std::fopen((out_dir / "elites.csv").string().c_str(), "r");
  if (!f) throw std::runtime_error("exemplars needs elites.csv: run the repertoire first");
  char line[1024];
  if (!std::fgets(line, sizeof line, f)) throw std::runtime_error("elites.csv is empty");
  while (std::fgets(line, sizeof line, f)) {
    std::vector<std::string> c;
    std::string cell;
    for (const char* q = line; *q && *q != '\n'; q++) {
      if (*q == ',') c.push_back(cell), cell.clear();
      else cell += *q;
    }
    c.push_back(cell);
    if (c.size() < 13 || c[0] != "MAP-Elites") continue;
    Adex::Params p = repertoire_cell(Real(6, 0.0));
    p.C_m = std::stod(c[7]) * u::pF;
    p.g_L = std::stod(c[8]) * u::nS;
    p.a = std::stod(c[9]) * u::nS;
    p.tau_w = std::stod(c[10]) * u::ms;
    p.b = std::stod(c[11]) * u::pA;
    p.V_reset = std::stod(c[12]) * u::mV;
    rows.push_back({c[5], std::stod(c[2]), std::stod(c[3]), p});
  }
  std::fclose(f);
  std::map<std::string, std::vector<const Row*>> by;
  for (const auto& r : rows) by[r.kind].push_back(&r);
  FILE* ft = std::fopen((out_dir / "traces.csv").string().c_str(), "w");
  FILE* fs = std::fopen((out_dir / "trace_spikes.csv").string().c_str(), "w");
  FILE* fx = std::fopen((out_dir / "exemplars.csv").string().c_str(), "w");
  std::fprintf(ft, "class,t_ms,V_mV\n");
  std::fprintf(fs, "class,t_ms\n");
  std::fprintf(fx, "class,elites,median_log_rate,median_cv,log_rate,cv,C_pF,g_L_nS,a_nS,tau_w_ms,b_pA,V_r_mV,class_again,log_rate_again,cv_again\n");
  auto median = [](std::vector<double> v) {
    std::sort(v.begin(), v.end());
    return v.size() % 2 ? v[v.size() / 2] : 0.5 * (v[v.size() / 2 - 1] + v[v.size() / 2]);
  };
  for (const auto& [kind, members] : by) {
    std::vector<double> lr, cv;
    for (const auto* r : members) lr.push_back(r->log_rate), cv.push_back(r->cv);
    const double ml = median(lr), mc = median(cv);
    const Row* best = members.front();
    double bd = std::numeric_limits<double>::infinity();
    for (const auto* r : members) {
      const double d = std::hypot((r->log_rate - ml) / (d_hi[0] - d_lo[0]), (r->cv - mc) / (d_hi[1] - d_lo[1]));
      if (d < bd) bd = d, best = r;
    }
    sn::Network net({H, 1, 1, stepping});
    auto& pop = net.add<Adex>("n", 1, best->p);
    net.set_bias(pop, {I_rep});
    auto& rec = net.record_state(pop, {0}, {0}, 3);
    auto& sp = net.record_spikes(pop);
    net.run(0.6);
    for (std::size_t i = 0; i < rec.times.size(); i++) std::fprintf(ft, "%s,%.2f,%.3f\n", kind.c_str(), double(rec.times[i]) * H * 1e3, rec.values[0][i] / u::mV);
    for (const auto& s2 : sp.spikes) std::fprintf(fs, "%s,%.2f\n", kind.c_str(), double(s2.time) * H * 1e3);
    // The parameters as elites.csv rounds them, simulated again for 1 s: the class and descriptors they give.
    const auto& q = best->p;
    const Behaviour again = behaviour(respond(q, {I_rep}, T_rep)[0]);
    std::fprintf(fx, "%s,%zu,%.4f,%.4f,%.4f,%.4f,%.3f,%.4f,%.4f,%.3f,%.3f,%.4f,%s,%.4f,%.4f\n", kind.c_str(), members.size(), ml, mc, best->log_rate, best->cv, q.C_m / u::pF, q.g_L / u::nS,
                 q.a / u::nS, q.tau_w / u::ms, q.b / u::pA, q.V_reset / u::mV, again.kind.c_str(), again.log_rate, again.cv);
  }
  std::fclose(ft);
  std::fclose(fs);
  std::fclose(fx);
  std::printf("exemplars: %zu classes\n", by.size());
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: fitting_neurons <results directory> [time|fit|landscape|repertoire|exemplars|refine|all]\n");
    return 2;
  }
  out_dir = argv[1];
  std::filesystem::create_directories(out_dir);
  const std::string mode = argc > 2 ? argv[2] : "all";
  const std::size_t threads = std::max(1u, std::thread::hardware_concurrency());
  calx::scheduler::ThreadPool pool(threads);
  const Target target = make_target();

  if (mode == "time") {
    for (const pd::Stepping how : {pd::Stepping{pd::Method::exact, 1, {1e-6, 1e-6}}, pd::Stepping{pd::Method::exact, 1, {1e-4, 1e-4}}, pd::Stepping{pd::Method::exact, 1, {1e-3, 1e-3}},
                                   pd::Stepping{pd::Method::rk4, 1, {}}, pd::Stepping{pd::Method::rk4, 4, {}}}) {
      stepping = pd::Stepping{pd::Method::exact, 1, {1e-6, 1e-6}};
      const auto ref = respond(fitted_cell(truth), currents, T_fit);
      stepping = how;
      const auto t2 = std::chrono::steady_clock::now();
      std::vector<std::vector<double>> s2;
      for (int i = 0; i < 5; i++) s2 = respond(fitted_cell(truth), currents, T_fit);
      const double ms = 1e3 * seconds_since(t2) / 5;
      std::size_t differ = 0, total = 0;
      double worst = 0;
      for (std::size_t k = 0; k < ref.size(); k++) {
        total += ref[k].size();
        if (s2[k].size() != ref[k].size()) { differ += 1000; continue; }
        for (std::size_t i = 0; i < ref[k].size(); i++) worst = std::max(worst, std::fabs(ref[k][i] - s2[k][i]));
      }
      std::printf("stepping method %d substeps %d tol %.0e: %.2f ms; count mismatch %zu of %zu spikes; worst shift %.2f ms\n", int(how.method), how.substeps, how.tolerance.relative, ms, differ, total, worst * 1e3);
    }
    stepping = pd::Stepping{pd::Method::exact, 1, {1e-6, 1e-6}};
    const auto t0 = std::chrono::steady_clock::now();
    const int n = 20;
    for (int i = 0; i < n; i++) (void)errors(respond(fitted_cell(truth), currents, T_fit), target);
    std::printf("one fitting evaluation (6 currents x 1 s): %.2f ms\n", 1e3 * seconds_since(t0) / n);
    const auto t1 = std::chrono::steady_clock::now();
    for (int i = 0; i < n; i++) (void)repertoire_evaluation(reference);
    std::printf("one repertoire evaluation (1 current x 1 s): %.2f ms\n", 1e3 * seconds_since(t1) / n);
    for (std::size_t k = 0; k < currents.size(); k++)
      std::printf("  target at %.0f pA: %zu spikes, initial %.1f Hz, steady %.1f Hz\n", currents[k] / u::pA, target.spikes[k].size(), target.rates[k].initial, target.rates[k].steady);
    const auto s = respond(repertoire_cell(reference), {I_rep}, T_rep)[0];
    std::printf("  reference under 500 pA: %zu spikes, class %s\n", s.size(), behaviour(s).kind.c_str());
    return 0;
  }
  if (mode == "fit" || mode == "all") fit(pool, target, 5, 1000);
  if (mode == "landscape" || mode == "all") landscape(pool, target, 41);
  if (mode == "repertoire" || mode == "all") repertoire(pool, 20000);
  if (mode == "exemplars" || mode == "all") exemplars();
  if (mode == "refine" || mode == "all") refine(pool, target, 600);
  return 0;
}
