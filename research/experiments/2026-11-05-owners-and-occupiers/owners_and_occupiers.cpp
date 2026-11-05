// Owners and occupiers of the Hoddle Grid's buildings, 2002-2024: the urban-science domain's `occupancy` on the City of
// Melbourne's CLUE (README.md). Four tasks:
//   (1) the uses of the buildings, simulated against the censuses after calibration by indirect inference;
//   (2) vacancy: where and why it persists, simulated histories beside the observed one;
//   (3) interventions on vacant space under a budget: the POMDP's ranked policy against simple rules, with the
//       Lagrangian bound;
//   (4) what the neighbours' modulation adds over the chain alone.
// Every number of the record is written here to results/, from which figures.py draws.
//
// Usage (from the repository's root):
//   owners_and_occupiers <occupancy-hoddle-grid.json> <city-events-melbourne.json> <results directory> [--quick]

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <numbers>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <calx/scheduler/thread_pool.hpp>
#include <science/json.hpp>
#include <urban/events/memory.hpp>
#include <urban/occupancy/calibration.hpp>
#include <urban/occupancy/model.hpp>
#include <urban/occupancy/planning.hpp>

namespace uo = urban::occupancy;
namespace ue = urban::events;
namespace sm = physis::statistical_mechanics;
namespace td = tyche::decision;
using science::Json;

namespace {

constexpr int FROM = 2002, TO = 2024;
constexpr std::size_t YEARS = TO - FROM + 1;
constexpr std::size_t K_NEIGHBOURS = 6;       // assumption: a building's neighbours are its six nearest
constexpr double AFFINITY_RADIUS = 60;        // m: the co-location's radius (about half a Hoddle block)
constexpr double INCENTIVE = 1.0;             // κ, assumption: an enrolled lettable vacancy is offered at twice the rate

double seconds_since(std::chrono::steady_clock::time_point t0) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); }

Json read_json(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read " + path);
  std::stringstream ss;
  ss << in.rdbuf();
  return Json::parse(ss.str());
}

void write_json(const std::string& path, const Json& j) {
  std::ofstream out(path, std::ios::binary);
  out << j.dump() << "\n";
}

Json numbers(const std::vector<double>& v) {
  Json a = Json::array();
  for (double x : v) a.push(std::isfinite(x) ? Json(std::round(x * 1e6) / 1e6) : Json());
  return a;
}

Json num(double x) { return std::isfinite(x) ? Json(x) : Json(); }

double quantile(std::vector<double> v, double q) {
  if (v.empty()) return std::nan("");
  std::sort(v.begin(), v.end());
  const double h = (v.size() - 1) * q;
  const std::size_t lo = static_cast<std::size_t>(std::floor(h)), hi = std::min(v.size() - 1, lo + 1);
  return v[lo] + (h - lo) * (v[hi] - v[lo]);
}

// ---- the record ------------------------------------------------------------------------------------------------------

struct Data {
  std::vector<std::string> ids, building, zone;
  std::vector<std::array<double, 2>> at;
  std::vector<long> block;
  std::vector<std::map<int, std::string>> use;  // by building: year ↦ use
  std::vector<std::map<int, long>> floors;
  std::vector<double> demand_raw;                // establishments outside the box, occupied, by year
  std::vector<std::pair<std::size_t, int>> permits;  // (building, first year)
  std::set<std::string> district;                // Event Spacetime's properties
};

Data read_data(const Json& P, const Json& E) {
  Data D;
  std::map<std::string, std::size_t> index;
  for (const auto& r : P["properties"]["rows"].as_array()) {
    index[r[0].text()] = D.ids.size();
    D.ids.push_back(r[0].text());
    D.block.push_back(r[2].is_number() ? static_cast<long>(r[2].as_number()) : -1);
    D.at.push_back({r[3][0].as_number(), r[3][1].as_number()});
    D.building.push_back(r[4].text());
    D.zone.push_back(r[6].text());
  }
  D.use.resize(D.ids.size()), D.floors.resize(D.ids.size());
  for (const auto& r : P["censuses"]["rows"].as_array()) {
    const auto it = index.find(r[1].text());
    if (it == index.end()) continue;
    const int y = static_cast<int>(r[0].as_number());
    if (y < FROM || y > TO) continue;
    D.use[it->second][y] = r[2].text();
    if (r[4].is_number()) D.floors[it->second][y] = static_cast<long>(r[4].as_number());
  }
  std::map<int, double> outside;
  for (const auto& r : P["demand"]["rows"].as_array()) outside[static_cast<int>(r[0].as_number())] = r[1].as_number();
  for (int y = FROM; y <= TO; y++) D.demand_raw.push_back(outside.at(y));
  for (const auto& r : P["permits"]["rows"].as_array()) {
    const auto it = index.find(r[1].text());
    if (it == index.end()) continue;
    D.permits.push_back({it->second, std::stoi(r[2].text().substr(0, 4))});
  }
  for (const auto& r : E["properties"]["rows"].as_array()) D.district.insert(r[0].text());
  return D;
}

// ---- what a panel shows: the statistics every figure compares --------------------------------------------------------

struct Shows {
  std::vector<std::vector<double>> count;    // [year][use]
  std::vector<double> vacant_years, lettable_years;  // by building: censuses vacant (V or W), lettable (V)
  std::vector<double> at_risk, exits;        // vacancy spells (whose start is seen) by their age in years, 1 … 22
  std::vector<double> nb_vac_risk, nb_vac_exit, nb_occ_risk, nb_occ_fall;  // by the neighbours' vacant share (4 bins)
  std::vector<double> clustering;            // by year: vacant neighbours' share about the vacant less about the occupied
  std::vector<double> turnover, family_turnover;  // by year (pairs y → y + 1): the share changing use, family
  std::vector<double> period_vacant, period_lettable;  // the pooled periodogram over complete series, k = 1 … 11
  double complete = 0;                       // buildings with a complete series
};

std::size_t bin_of(double share) { return share <= 0 ? 0 : share <= 0.2 ? 1 : share <= 0.4 ? 2 : 3; }

Shows show(const uo::Uses& S, const uo::Space& X, const uo::Seen& seen) {
  const std::size_t n = seen.size(), N = S.states();
  Shows W;
  std::vector<std::vector<long>> at(n, std::vector<long>(YEARS, -1));
  for (std::size_t i = 0; i < n; i++)
    for (const auto& [y, u] : seen[i])
      if (y >= FROM && y <= TO) at[i][static_cast<std::size_t>(y - FROM)] = static_cast<long>(u);
  W.count.assign(YEARS, std::vector<double>(N, 0.0));
  W.vacant_years.assign(n, 0.0), W.lettable_years.assign(n, 0.0);
  W.at_risk.assign(YEARS, 0.0), W.exits.assign(YEARS, 0.0);
  W.nb_vac_risk.assign(4, 0.0), W.nb_vac_exit.assign(4, 0.0), W.nb_occ_risk.assign(4, 0.0), W.nb_occ_fall.assign(4, 0.0);
  W.clustering.assign(YEARS, std::nan("")), W.turnover.assign(YEARS - 1, 0.0), W.family_turnover.assign(YEARS - 1, 0.0);
  const auto vac = [&](long u) { return u >= 0 && S.is_vacant(static_cast<std::size_t>(u)); };
  const auto share = [&](std::size_t i, std::size_t t, bool& any) {
    double v = 0, k = 0;
    for (const std::size_t j : X.neighbours[i])
      if (at[j][t] >= 0) k += 1, v += vac(at[j][t]);
    any = k > 0;
    return any ? v / k : 0.0;
  };
  std::vector<double> pairs(YEARS - 1, 0.0);
  for (std::size_t i = 0; i < n; i++) {
    long age = 0;  // the current vacancy spell's age in years, 0 if not in one, −1 if its start is not seen
    for (std::size_t t = 0; t < YEARS; t++) {
      const long u = at[i][t];
      if (u < 0) { age = 0; continue; }
      W.count[t][static_cast<std::size_t>(u)] += 1;
      W.vacant_years[i] += vac(u), W.lettable_years[i] += S.lump[static_cast<std::size_t>(u)] == uo::vacant;
      if (vac(u)) {
        if (age == 0) age = (t > 0 && at[i][t - 1] >= 0 && !vac(at[i][t - 1])) ? 1 : -1;
        else if (age > 0) age++;
      } else age = 0;
      if (t + 1 < YEARS && at[i][t + 1] >= 0) {
        const long w = at[i][t + 1];
        pairs[t] += 1;
        W.turnover[t] += u != w, W.family_turnover[t] += S.family[static_cast<std::size_t>(u)] != S.family[static_cast<std::size_t>(w)];
        if (vac(u) && age > 0) W.at_risk[static_cast<std::size_t>(age)] += 1, W.exits[static_cast<std::size_t>(age)] += !vac(w);
        bool any = false;
        const double sh = share(i, t, any);
        if (any) {
          const std::size_t b = bin_of(sh);
          if (vac(u)) W.nb_vac_risk[b] += 1, W.nb_vac_exit[b] += !vac(w);
          else W.nb_occ_risk[b] += 1, W.nb_occ_fall[b] += vac(w);
        }
      }
    }
  }
  for (std::size_t t = 0; t + 1 < YEARS; t++)
    if (pairs[t] > 0) W.turnover[t] /= pairs[t], W.family_turnover[t] /= pairs[t];
  for (std::size_t t = 0; t < YEARS; t++) {
    double sv = 0, nv = 0, so = 0, no = 0;
    for (std::size_t i = 0; i < n; i++) {
      if (at[i][t] < 0) continue;
      bool any = false;
      const double sh = share(i, t, any);
      if (!any) continue;
      if (vac(at[i][t])) sv += sh, nv += 1;
      else so += sh, no += 1;
    }
    if (nv > 0 && no > 0) W.clustering[t] = sv / nv - so / no;
  }
  // The pooled periodogram of the centred indicators over the buildings seen every year.
  const std::size_t K = YEARS / 2;
  W.period_vacant.assign(K, 0.0), W.period_lettable.assign(K, 0.0);
  for (std::size_t i = 0; i < n; i++) {
    bool full = true;
    for (std::size_t t = 0; t < YEARS && full; t++) full = at[i][t] >= 0;
    if (!full) continue;
    W.complete += 1;
    for (int which = 0; which < 2; which++) {
      std::vector<double> x(YEARS);
      double mean = 0;
      for (std::size_t t = 0; t < YEARS; t++) {
        const auto u = static_cast<std::size_t>(at[i][t]);
        x[t] = which == 0 ? S.is_vacant(u) : S.lump[u] == uo::vacant;
        mean += x[t] / YEARS;
      }
      for (std::size_t k = 1; k <= K; k++) {
        std::complex<double> z = 0;
        for (std::size_t t = 0; t < YEARS; t++) z += (x[t] - mean) * std::polar(1.0, -2 * std::numbers::pi * static_cast<double>(k * t) / YEARS);
        (which == 0 ? W.period_vacant : W.period_lettable)[k - 1] += std::norm(z) / YEARS;
      }
    }
  }
  if (W.complete > 0)
    for (std::size_t k = 0; k < K; k++) W.period_vacant[k] /= W.complete, W.period_lettable[k] /= W.complete;
  return W;
}

// The chain's spectrum of an indicator f at annual sampling: Σ_{m≥1} c_m (1 − r²)/(1 − 2 r cos ω + r²), r = e^{−μ_m},
// c_m = ⟨f, Q_m f⟩_π (theory/occupancy.md, Proposition 3.1).
std::vector<double> chain_spectrum(const sm::NestedChain& X, const std::vector<double>& f, std::vector<double>* weights = nullptr) {
  const std::size_t N = X.C.states(), L = X.levels();
  std::vector<double> c(L + 1, 0.0);
  for (std::size_t m = 1; m <= L; m++)
    for (std::size_t x = 0; x < N; x++)
      for (std::size_t y = 0; y < N; y++) c[m] += X.C.law[x] * f[x] * X.mode(m, x, y) * f[y];
  if (weights) *weights = c;
  std::vector<double> out;
  for (std::size_t k = 1; k <= YEARS / 2; k++) {
    const double w = 2 * std::numbers::pi * k / YEARS;
    double s = 0;
    for (std::size_t m = 1; m <= L; m++) {
      const double r = std::exp(-X.mu[m]);
      s += c[m] * (1 - r * r) / (1 - 2 * r * std::cos(w) + r * r);
    }
    out.push_back(s);
  }
  return out;
}

// Quantiles (5, 50, 95 %) across replicates of a vector quantity, and the mean.
Json bands(const std::vector<std::vector<double>>& reps) {
  const std::size_t m = reps.empty() ? 0 : reps[0].size();
  std::vector<double> lo(m), mid(m), hi(m), mean(m, 0.0);
  for (std::size_t k = 0; k < m; k++) {
    std::vector<double> v;
    for (const auto& r : reps)
      if (std::isfinite(r[k])) v.push_back(r[k]), mean[k] += r[k];
    mean[k] = v.empty() ? std::nan("") : mean[k] / v.size();
    lo[k] = quantile(v, 0.05), mid[k] = quantile(v, 0.5), hi[k] = quantile(v, 0.95);
  }
  return Json::object({{"q05", numbers(lo)}, {"q50", numbers(mid)}, {"q95", numbers(hi)}, {"mean", numbers(mean)}});
}

Json calibration_json(const uo::Calibration& C, std::size_t L) {
  Json trace = Json::array();
  for (const auto& [th, j] : C.trace) trace.push(Json::object({{"theta", numbers(th)}, {"J", j}}));
  Json cov = Json::array();
  for (const auto& r : C.covariance) cov.push(numbers(r));
  Json names = Json::array();
  for (const auto& s : uo::auxiliary_names(L)) names.push(s);
  return Json::object({{"theta", numbers(C.theta)}, {"error", numbers(C.error)}, {"covariance", cov}, {"lambda", numbers(C.parameters.lambda)},
                       {"vacate", C.parameters.vacate}, {"let", C.parameters.let}, {"returning", C.parameters.returning},
                       {"beta_neighbours", C.parameters.beta_neighbours}, {"beta_demand", C.parameters.beta_demand},
                       {"psi_names", names}, {"psi_observed", numbers(C.psi_observed)}, {"psi_mean", numbers(C.psi_mean)}, {"psi_sd", numbers(C.psi_sd)},
                       {"weight", numbers(C.weight)}, {"J", num(C.objective)}, {"hansen", num(C.hansen)}, {"degrees", static_cast<double>(C.degrees)},
                       {"replicates_sigma", C.replicates_sigma},
                       {"evaluations", C.evaluations}, {"converged", C.converged}, {"trace", trace}});
}

Json memory_json(const ue::UseMemory& M) {
  return Json::object({{"lambda", numbers(M.fit.lambda)}, {"lambda_error", numbers(M.fit.lambda_error)}, {"mu", numbers(M.fit.mu)},
                       {"half_life", numbers(M.half_life)}, {"half_life_error", numbers(M.half_life_error)}, {"transitions", M.transitions},
                       {"units", static_cast<double>(M.units)}, {"uses", static_cast<double>(M.uses)}, {"ratio", num(M.ratio)},
                       {"degrees", static_cast<double>(M.degrees)}, {"log_likelihood", M.fit.log_likelihood}, {"agree", numbers(M.agree)}});
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "usage: owners_and_occupiers <occupancy-hoddle-grid.json> <city-events-melbourne.json> <results directory> [--quick]\n";
    return 2;
  }
  const bool quick = argc > 4 && std::string(argv[4]) == "--quick";
  const std::string out = argv[3];
  const auto t_start = std::chrono::steady_clock::now();
  const auto P = read_json(argv[1]), E = read_json(argv[2]);
  const Data D = read_data(P, E);
  const std::size_t n = D.ids.size();
  calx::scheduler::ThreadPool pool(4);

  // ---- the uses, the space, the affinity, the lives -----------------------------------------------------------------
  std::vector<ue::Observation> O;
  for (std::size_t i = 0; i < n; i++)
    for (const auto& [y, u] : D.use[i]) O.push_back({D.ids[i], u, y});
  const uo::Uses S = uo::uses_of(ue::prefix_classification(ue::use_counts(O)));
  const std::size_t N = S.states(), L = S.levels();
  const uo::Space X = uo::nearest(D.at, K_NEIGHBOURS);
  uo::Seen seen(n);
  std::vector<uo::Life> lives(n);
  for (std::size_t i = 0; i < n; i++) {
    for (const auto& [y, u] : D.use[i]) seen[i].push_back({y, S.U.of(u)});
    if (seen[i].empty()) continue;
    lives[i].first = seen[i].front().first, lives[i].last = seen[i].back().first, lives[i].initial = seen[i].front().second;
    for (const auto& o : seen[i]) lives[i].censuses.push_back(o.first);
  }
  std::vector<std::size_t> use2002(n, N);
  for (std::size_t i = 0; i < n; i++)
    if (D.use[i].count(FROM)) use2002[i] = S.U.of(D.use[i].at(FROM));
  const auto affinity = uo::family_affinity(S, X, use2002, AFFINITY_RADIUS, 2.0, 99, 20261105);
  // The demand index: the occupied establishments outside the box, ln, standardised over the window.
  std::vector<double> demand(YEARS);
  {
    double m = 0, v = 0;
    for (double x : D.demand_raw) m += std::log(x) / YEARS;
    for (double x : D.demand_raw) v += (std::log(x) - m) * (std::log(x) - m) / (YEARS - 1);
    for (std::size_t t = 0; t < YEARS; t++) demand[t] = (std::log(D.demand_raw[t]) - m) / std::sqrt(v);
  }
  {
    Json uses = Json::array();
    for (std::size_t y = 0; y < N; y++)
      uses.push(Json::object({{"use", S.U.uses[y]}, {"family", S.U.block_names[1][y]}, {"lump", static_cast<double>(S.lump[y])}, {"law", S.U.C.law[y]}}));
    std::vector<double> a(affinity.begin(), affinity.end());
    write_json(out + "/uses.json", Json::object({{"uses", uses}, {"affinity", numbers(a)}, {"demand", numbers(demand)}, {"demand_raw", numbers(D.demand_raw)},
                                                 {"years", Json(std::vector<double>{FROM, TO})}, {"neighbours", static_cast<double>(K_NEIGHBOURS)},
                                                 {"affinity_radius", AFFINITY_RADIUS}, {"buildings", static_cast<double>(n)}}));
  }
  std::cout << n << " buildings, " << O.size() << " observations, " << N << " uses in " << S.U.C.block[1].size() << " states\n";

  // ---- (1) the memory of the uses, and the calibration ----------------------------------------------------------------
  const auto T_all = ue::transitions_of(S.U, O);
  const auto Mem = ue::remember(S.U, T_all, n);
  std::vector<ue::Observation> Od;
  std::set<std::string> units_d;
  for (const auto& o : O)
    if (D.district.count(o.unit)) Od.push_back(o), units_d.insert(o.unit);
  const auto Ud = ue::prefix_classification(ue::use_counts(Od));
  const auto Mem_d = ue::remember(Ud, ue::transitions_of(Ud, Od), units_d.size());
  std::printf("memory: λ %.4f %.4f, half-lives %.2f ± %.2f, %.2f ± %.2f years; district %.2f, %.2f\n", Mem.fit.lambda[0], Mem.fit.lambda[1], Mem.half_life[0], Mem.half_life_error[0],
              Mem.half_life[1], Mem.half_life_error[1], Mem_d.half_life[0], Mem_d.half_life[1]);
  write_json(out + "/memory.json", Json::object({{"hoddle_grid", memory_json(Mem)}, {"district", memory_json(Mem_d)}}));

  uo::Calibrand Cd{&S, &X, affinity, lives, demand, FROM, TO, INCENTIVE};
  const auto A_obs = uo::auxiliary(S, X, affinity, seen, demand, FROM, Mem.fit.lambda);
  const auto start = uo::starting_parameters(S, A_obs, INCENTIVE);
  std::printf("start: λ %.4f %.4f ν %.4f η %.4f ρ %.3f\n", start.lambda[0], start.lambda[1], start.vacate, start.let, start.returning);
  uo::CalibrationOptions co;
  co.replicates = quick ? 6 : 16, co.seed = 20261105, co.max_evaluations = quick ? 40 : 300, co.max_iterations = quick ? 2 : 8, co.pool = &pool;
  const auto report = [&](const char* name, const uo::Calibration& C, double secs) {
    std::printf("calibrated (%s) in %.0f s: λ %.4f %.4f ν %.4f η %.4f ρ %.3f βn %.3f ± %.3f βd %.3f ± %.3f; J %.2f, Hansen %.2f on %zu\n", name, secs,
                C.parameters.lambda[0], C.parameters.lambda[1], C.parameters.vacate, C.parameters.let, C.parameters.returning, C.theta[L + 3], C.error[L + 3],
                C.theta[L + 4], C.error[L + 4], C.objective, C.hansen, C.degrees);
    std::fflush(stdout);
  };
  auto t0 = std::chrono::steady_clock::now();
  const auto Cm = uo::calibrate(Cd, seen, start, co);
  report("the full model", Cm, seconds_since(t0));
  co.neighbours = false;
  t0 = std::chrono::steady_clock::now();
  const auto Cv = uo::calibrate(Cd, seen, start, co);
  report("without the neighbours", Cv, seconds_since(t0));
  co.vacancy = false, co.demand = false;
  t0 = std::chrono::steady_clock::now();
  const auto C0 = uo::calibrate(Cd, seen, start, co);
  report("the chain alone", C0, seconds_since(t0));
  // The distance-metric tests at the full model's Σ: without the neighbours, and the chain alone, against it.
  double distance_neighbours = 0, distance_chain = 0;
  {
    const int Sp = Cm.replicates_sigma;
    const auto Eh = uo::ensemble(Cd, Cm.theta, Sp, co.seed, start.lambda, &pool);
    const auto Ev = uo::ensemble(Cd, Cv.theta, Sp, co.seed, start.lambda, &pool);
    const auto E0 = uo::ensemble(Cd, C0.theta, Sp, co.seed, start.lambda, &pool);
    const auto q = [&](const std::vector<double>& mean) { return uo::detail::hansen(Cm.psi_observed, mean, Cm.psi_covariance, Sp); };
    distance_neighbours = q(Ev.mean) - q(Eh.mean), distance_chain = q(E0.mean) - q(Eh.mean);
  }
  std::printf("distance metric: without the neighbours %.1f (1 df), the chain alone %.1f (5 df)\n", distance_neighbours, distance_chain);
  write_json(out + "/calibration.json", Json::object({{"modulated", calibration_json(Cm, L)}, {"without_neighbours", calibration_json(Cv, L)},
                                                      {"chain_alone", calibration_json(C0, L)}, {"distance_neighbours", num(distance_neighbours)},
                                                      {"distance_chain", num(distance_chain)}, {"replicates", co.replicates}, {"incentive", INCENTIVE},
                                                      {"start", numbers(uo::theta_of(start))}}));

  // ---- ensembles at the two calibrations: tasks (1), (2), (4) and the spectra ------------------------------------------
  const int R = quick ? 20 : 200;
  const auto observed = show(S, X, seen);
  struct Run {
    std::vector<Shows> shows;
    std::vector<double> mean_vacant_years, mean_lettable_years;  // by building
  };
  const auto run_ensemble = [&](const uo::Parameters& theta, std::uint64_t seed) {
    Run out;
    out.shows.resize(static_cast<std::size_t>(R));
    uo::Model M{&S, &X, affinity, theta};
    calx::scheduler::parallel_for(pool, static_cast<std::size_t>(R), 1, [&](std::size_t b, std::size_t e, std::size_t) {
      for (std::size_t r = b; r < e; r++) out.shows[r] = show(S, X, uo::simulate(M, lives, demand, FROM, TO, seed + r).seen);
    });
    out.mean_vacant_years.assign(n, 0.0), out.mean_lettable_years.assign(n, 0.0);
    for (const auto& w : out.shows)
      for (std::size_t i = 0; i < n; i++) out.mean_vacant_years[i] += w.vacant_years[i] / R, out.mean_lettable_years[i] += w.lettable_years[i] / R;
    return out;
  };
  t0 = std::chrono::steady_clock::now();
  const Run Rm = run_ensemble(Cm.parameters, 90000), Rv = run_ensemble(Cv.parameters, 90000), R0 = run_ensemble(C0.parameters, 90000);
  std::printf("ensembles of %d in %.0f s\n", R, seconds_since(t0));
  const auto ensemble_json = [&](const Run& run) {
    std::vector<std::vector<double>> vac_share, let_share, clustering, turnover, fam_turnover, pv, pl, haz, nbv, nbo;
    std::vector<std::vector<double>> per_use(N);
    for (const auto& w : run.shows) {
      std::vector<double> v(YEARS), l(YEARS);
      for (std::size_t t = 0; t < YEARS; t++) {
        double tot = 0, vv = 0, ll = 0;
        for (std::size_t y = 0; y < N; y++) tot += w.count[t][y], vv += S.is_vacant(y) ? w.count[t][y] : 0, ll += S.lump[y] == uo::vacant ? w.count[t][y] : 0;
        v[t] = vv / tot, l[t] = ll / tot;
      }
      vac_share.push_back(v), let_share.push_back(l), clustering.push_back(w.clustering), turnover.push_back(w.turnover), fam_turnover.push_back(w.family_turnover);
      pv.push_back(w.period_vacant), pl.push_back(w.period_lettable);
      std::vector<double> h(YEARS, std::nan("")), b1(4, std::nan("")), b2(4, std::nan(""));
      for (std::size_t k = 1; k < YEARS; k++)
        if (w.at_risk[k] >= 5) h[k] = w.exits[k] / w.at_risk[k];
      for (std::size_t b = 0; b < 4; b++) {
        if (w.nb_vac_risk[b] > 0) b1[b] = w.nb_vac_exit[b] / w.nb_vac_risk[b];
        if (w.nb_occ_risk[b] > 0) b2[b] = w.nb_occ_fall[b] / w.nb_occ_risk[b];
      }
      haz.push_back(h), nbv.push_back(b1), nbo.push_back(b2);
      for (std::size_t y = 0; y < N; y++) {
        std::vector<double> c(YEARS);
        for (std::size_t t = 0; t < YEARS; t++) c[t] = w.count[t][y];
        per_use[y].insert(per_use[y].end(), c.begin(), c.end());
      }
    }
    // Per use and year, the mean count across replicates.
    Json uses = Json::array();
    for (std::size_t y = 0; y < N; y++) {
      std::vector<double> mean(YEARS, 0.0);
      for (std::size_t r = 0; r < run.shows.size(); r++)
        for (std::size_t t = 0; t < YEARS; t++) mean[t] += run.shows[r].count[t][y] / run.shows.size();
      uses.push(numbers(mean));
    }
    return Json::object({{"vacant_share", bands(vac_share)}, {"lettable_share", bands(let_share)}, {"clustering", bands(clustering)}, {"turnover", bands(turnover)},
                         {"family_turnover", bands(fam_turnover)}, {"periodogram_vacant", bands(pv)}, {"periodogram_lettable", bands(pl)},
                         {"exit_by_age", bands(haz)}, {"exit_by_neighbours", bands(nbv)}, {"fall_by_neighbours", bands(nbo)}, {"use_counts", uses}});
  };
  {
    std::vector<double> v(YEARS), l(YEARS);
    for (std::size_t t = 0; t < YEARS; t++) {
      double tot = 0, vv = 0, ll = 0;
      for (std::size_t y = 0; y < N; y++) tot += observed.count[t][y], vv += S.is_vacant(y) ? observed.count[t][y] : 0, ll += S.lump[y] == uo::vacant ? observed.count[t][y] : 0;
      v[t] = vv / tot, l[t] = ll / tot;
    }
    Json counts = Json::array();
    for (std::size_t t = 0; t < YEARS; t++) counts.push(numbers(observed.count[t]));
    // The chain's closed-form spectra (the observed fit), with each mode's weight.
    const sm::NestedChain chain(S.U.C, Mem.fit.lambda);
    std::vector<double> fv(N), fl(N), wv, wl;
    for (std::size_t y = 0; y < N; y++) fv[y] = S.is_vacant(y), fl[y] = S.lump[y] == uo::vacant;
    const auto sv = chain_spectrum(chain, fv, &wv), sl = chain_spectrum(chain, fl, &wl);
    write_json(out + "/ensembles.json",
               Json::object({{"observed", Json::object({{"vacant_share", numbers(v)}, {"lettable_share", numbers(l)}, {"clustering", numbers(observed.clustering)},
                                                       {"turnover", numbers(observed.turnover)}, {"family_turnover", numbers(observed.family_turnover)},
                                                       {"periodogram_vacant", numbers(observed.period_vacant)}, {"periodogram_lettable", numbers(observed.period_lettable)},
                                                       {"at_risk", numbers(observed.at_risk)}, {"exits", numbers(observed.exits)},
                                                       {"nb_vac_risk", numbers(observed.nb_vac_risk)}, {"nb_vac_exit", numbers(observed.nb_vac_exit)},
                                                       {"nb_occ_risk", numbers(observed.nb_occ_risk)}, {"nb_occ_fall", numbers(observed.nb_occ_fall)},
                                                       {"counts", counts}, {"complete", observed.complete}})},
                             {"modulated", ensemble_json(Rm)}, {"without_neighbours", ensemble_json(Rv)}, {"chain_alone", ensemble_json(R0)}, {"replicates", R},
                             {"spectrum", Json::object({{"vacant", numbers(sv)}, {"lettable", numbers(sl)}, {"weights_vacant", numbers(wv)}, {"weights_lettable", numbers(wl)},
                                                        {"mu", numbers(chain.mu)}})}}));
  }

  // ---- the permits for changes of use against the censuses ----------------------------------------------------------
  {
    // A matter in year y at a building seen at y − 1 and y + 1: did its use change between them? Against every such
    // pair of censuses two years apart.
    double with = 0, changed = 0, base_pairs = 0, base_changed = 0, vacant_before = 0, base_vacant = 0;
    for (const auto& [i, y] : D.permits) {
      const auto a = D.use[i].find(y - 1), b = D.use[i].find(y + 1);
      if (a == D.use[i].end() || b == D.use[i].end()) continue;
      with += 1, changed += a->second != b->second, vacant_before += S.is_vacant(S.U.of(a->second));
    }
    for (std::size_t i = 0; i < n; i++)
      for (const auto& [y, u] : D.use[i]) {
        const auto b = D.use[i].find(y + 2);
        if (b == D.use[i].end()) continue;
        base_pairs += 1, base_changed += u != b->second, base_vacant += S.is_vacant(S.U.of(u));
      }
    std::printf("permits: %zu matters; %.0f with censuses about them, %.3f changed use (base %.3f); vacant before %.3f (base %.3f)\n", D.permits.size(), with,
                changed / with, base_changed / base_pairs, vacant_before / with, base_vacant / base_pairs);
    write_json(out + "/permits.json", Json::object({{"matters", static_cast<double>(D.permits.size())}, {"with_censuses", with}, {"changed", changed},
                                                    {"vacant_before", vacant_before}, {"base_pairs", base_pairs}, {"base_changed", base_changed}, {"base_vacant", base_vacant}}));
  }

  // ---- (3) the planner, from the census of 2024 ----------------------------------------------------------------------
  const int START = TO + 1;
  uo::State s24;
  s24.use.assign(n, 0), s24.alive.assign(n, 0), s24.enrolled.assign(n, 0), s24.held.assign(n, N), s24.demand = demand.back();
  std::vector<char> planned(n, 0);
  for (std::size_t i = 0; i < n; i++)
    if (!seen[i].empty()) {
      s24.use[i] = seen[i].back().second;
      for (std::size_t k = seen[i].size(); k-- > 0;)  // the family it last held occupied, as the owner knows it
        if (S.lump[seen[i][k].second] == uo::occupied) { s24.held[i] = S.family[seen[i][k].second]; break; }
      if (seen[i].back().first == TO) s24.alive[i] = 1, planned[i] = 1;
    }
  const uo::PlanAssumptions A;
  const auto plan_for = [&](const uo::Parameters& theta, const std::vector<double>& budgets, const std::string& label) {
    uo::Model M{&S, &X, affinity, theta};
    std::vector<td::Belief> beliefs(n, td::Belief(uo::lumps, 1.0 / uo::lumps));
    uo::Knowledge K;
    K.lump.assign(n, -1), K.year.assign(n, TO), K.since.assign(n, TO), K.belief.assign(n, td::Belief(uo::lumps, 1.0 / uo::lumps));
    for (std::size_t i = 0; i < n; i++) {
      if (!planned[i]) continue;
      const std::size_t lump = S.lump[s24.use[i]];
      K.lump[i] = static_cast<int>(lump);
      if (lump == uo::vacant) {
        int since = TO;
        for (std::size_t k = seen[i].size(); k-- > 0;) {
          if (S.lump[seen[i][k].second] != uo::vacant) break;
          since = seen[i][k].first;
        }
        K.since[i] = since;
      }
      beliefs[i] = K.belief[i] = uo::belief_after(uo::building_model(M, s24, i, A), lump, 1);
    }
    Json rows = Json::array();
    const std::vector<uo::Rule> rules = {uo::Rule::nothing, uo::Rule::every_vacancy, uo::Rule::longest_vacant, uo::Rule::by_neighbours, uo::Rule::pomdp};
    Json first_actions = Json::object();
    for (const double B : budgets) {
      auto tp = std::chrono::steady_clock::now();
      const auto C = uo::plan_city(M, s24, planned, beliefs, B, A);
      const double t_plan = seconds_since(tp);
      Json per_rule = Json::array();
      std::vector<std::vector<double>> abm_paired(rules.size()), dec_paired(rules.size());
      const int episodes_abm = quick ? 4 : 160, episodes_dec = quick ? 4 : 80, horizon_dec = quick ? 10 : 90, horizon_abm = 10;
      for (std::size_t r = 0; r < rules.size(); r++) {
        tp = std::chrono::steady_clock::now();
        // Both evaluations episode by episode on common seeds, so that each rule is paired against doing nothing; the
        // episodes run in parallel, each from its own seed.
        std::vector<double> dec(static_cast<std::size_t>(episodes_dec)), years(static_cast<std::size_t>(episodes_abm));
        double spent = 0;
        std::vector<double> spent_e(static_cast<std::size_t>(episodes_dec));
        calx::scheduler::parallel_for(pool, static_cast<std::size_t>(episodes_dec), 1, [&](std::size_t b, std::size_t e, std::size_t) {
          for (std::size_t k = b; k < e; k++) {
            const auto Ed = uo::evaluate_decoupled(C, rules[r], X, K, START, horizon_dec, 1, 4242 + k);
            dec[k] = Ed.mean, spent_e[k] = Ed.spent;
          }
        });
        calx::scheduler::parallel_for(pool, static_cast<std::size_t>(episodes_abm), 1, [&](std::size_t b, std::size_t e, std::size_t) {
          for (std::size_t k = b; k < e; k++) years[k] = uo::evaluate_abm(M, C, rules[r], s24, K, START, horizon_abm, 1, 7000 + k).years;
        });
        for (double x : spent_e) spent += x / episodes_dec;
        uo::AbmRun shares;
        uo::evaluate_abm(M, C, rules[r], s24, K, START, horizon_abm, quick ? 2 : 10, 7000, &shares);
        abm_paired[r] = years, dec_paired[r] = dec;
        double my = 0, md = 0, sd = 0;
        for (double y : years) my += y / years.size();
        for (double y : dec) md += y / dec.size();
        for (double y : dec) sd += (y - md) * (y - md) / (dec.size() - 1);
        per_rule.push(Json::object({{"rule", uo::rule_name(rules[r])}, {"decoupled", md}, {"decoupled_error", std::sqrt(sd / dec.size())},
                                    {"spent", spent}, {"abm_years", my}, {"abm_episodes", numbers(years)}, {"abm_shares", numbers(shares.vacant_by_year)},
                                    {"seconds", seconds_since(tp)}}));
        if (B == budgets[1]) {
          const auto act = uo::decide(rules[r], C, K, X, START, 4242);
          std::vector<double> a(n, -1.0);
          for (std::size_t k = 0; k < C.planned.size(); k++) a[C.planned[k]] = static_cast<double>(act[k]);
          first_actions.set(uo::rule_name(rules[r]), numbers(a));
        }
      }
      // Averted years, paired against doing nothing: in the agent-based model (undiscounted, 2025-2034), and in the
      // planner's decoupled model (discounted; the policy for 90 years, then nothing: γ^90 < 0.01 of its value left out).
      const auto paired = [&](const std::vector<std::vector<double>>& x, bool gain) {
        Json out = Json::array();
        for (std::size_t r = 0; r < rules.size(); r++) {
          double m = 0, q = 0;
          const std::size_t ne = x[r].size();
          for (std::size_t e = 0; e < ne; e++) {
            const double d = gain ? x[r][e] - x[0][e] : x[0][e] - x[r][e];
            m += d, q += d * d;
          }
          m /= ne;
          const double se = ne > 1 ? std::sqrt(std::max(0.0, (q - ne * m * m) / (ne - 1) / ne)) : 0;
          out.push(Json::object({{"rule", uo::rule_name(rules[r])}, {"years", m}, {"error", se}}));
        }
        return out;
      };
      rows.push(Json::object({{"budget", B}, {"bound", C.solution.upper}, {"lambda", C.solution.lambda}, {"nothing", C.nothing_value},
                              {"models", static_cast<double>(C.models.size())}, {"planned", static_cast<double>(C.planned.size())}, {"plan_seconds", t_plan},
                              {"rules", per_rule}, {"abm_averted", paired(abm_paired, false)}, {"decoupled_averted", paired(dec_paired, true)},
                              {"episodes_decoupled", episodes_dec}, {"episodes_abm", episodes_abm}, {"horizon_decoupled", horizon_dec}}));
      std::printf("%s B = %.0f: bound %.1f, nothing %.1f (%zu models, %.0f s)\n", label.c_str(), B, C.solution.upper, C.nothing_value, C.models.size(), t_plan);
      for (const auto& x : per_rule.as_array())
        std::printf("   %-28s decoupled %.1f ± %.1f   abm %.2f years over 10\n", x["rule"].text().c_str(), x["decoupled"].as_number(), x["decoupled_error"].as_number(),
                    x["abm_years"].as_number());
      std::fflush(stdout);
    }
    std::vector<double> bv(n, std::nan(""));
    for (std::size_t i = 0; i < n; i++)
      if (planned[i]) bv[i] = K.belief[i][uo::vacant];
    return Json::object({{"budgets", rows}, {"first_actions", first_actions}, {"belief_vacant", numbers(bv)}});
  };
  t0 = std::chrono::steady_clock::now();
  const std::vector<double> budgets = quick ? std::vector<double>{5, 10} : std::vector<double>{5, 10, 20, 40};
  const auto Pm = plan_for(Cm.parameters, budgets, "the full model");
  const auto Pv = plan_for(Cv.parameters, budgets, "without the neighbours");
  std::printf("planning in %.0f s\n", seconds_since(t0));
  write_json(out + "/planning.json", Json::object({{"modulated", Pm}, {"without_neighbours", Pv}, {"start", START},
                                                   {"assumptions", Json::object({{"survey", A.survey}, {"cost", numbers({A.cost[0], A.cost[1], A.cost[2]})},
                                                                                 {"discount", A.discount}, {"incentive", INCENTIVE}})}}));

  // ---- the buildings ---------------------------------------------------------------------------------------------------
  {
    std::ofstream f(out + "/buildings.csv");
    f << "id,x,y,block,zone,building,first,last,censuses,use_first,use_last,observed_vacant,observed_lettable,expected_vacant_modulated,"
         "expected_vacant_without_neighbours,expected_vacant_chain,expected_lettable_modulated,permits,district\n";
    std::vector<int> permits(n, 0);
    for (const auto& [i, y] : D.permits) permits[i]++;
    for (std::size_t i = 0; i < n; i++) {
      if (seen[i].empty()) continue;
      char line[512];
      std::snprintf(line, sizeof line, "%s,%.1f,%.1f,%ld,%s,%s,%d,%d,%zu,\"%s\",\"%s\",%.0f,%.0f,%.3f,%.3f,%.3f,%.3f,%d,%d\n", D.ids[i].c_str(), D.at[i][0], D.at[i][1],
                    D.block[i], D.zone[i].c_str(), D.building[i].c_str(), lives[i].first, lives[i].last, seen[i].size(), S.U.uses[seen[i].front().second].c_str(),
                    S.U.uses[seen[i].back().second].c_str(), observed.vacant_years[i], observed.lettable_years[i], Rm.mean_vacant_years[i], Rv.mean_vacant_years[i],
                    R0.mean_vacant_years[i], Rm.mean_lettable_years[i], permits[i], D.district.count(D.ids[i]) ? 1 : 0);
      f << line;
    }
  }

  // ---- a simulated record of Event Spacetime's district, for drawing beside the observed one --------------------------
  {
    uo::Model M{&S, &X, affinity, Cm.parameters};
    const std::uint64_t seed = 20261105;
    const auto H = uo::simulate(M, lives, demand, FROM, TO, seed);
    const auto Rec = uo::to_record(S, X, D.ids, H, seed);
    const auto bad = ue::problems(Rec);
    Json props = Json::array(), cens = Json::array(), events = Json::array();
    std::size_t kept = 0;
    for (std::size_t i = 0; i < n; i++) {
      if (!D.district.count(D.ids[i]) || H.seen[i].empty()) continue;
      kept++;
      props.push(Json::array({D.ids[i], numbers({D.at[i][0], D.at[i][1]})}));
      for (const auto& [y, u] : H.seen[i]) cens.push(Json::array({y, D.ids[i], S.U.uses[u]}));
    }
    std::size_t k = 0;
    for (const auto& t : H.changes) {
      const auto& c = t.change;
      if (c.from == c.to || !D.district.count(D.ids[c.building])) continue;
      const int y = static_cast<int>(std::floor(t.time));
      const double day = (t.time - y) * (ue::mjd(y + 1) - ue::mjd(y));
      const auto date = [&](double mjd) {
        // The civil date of an MJD (Fliegel and Van Flandern's inverse).
        long l = static_cast<long>(std::floor(mjd)) + 2400001 + 68569;
        const long nn = 4 * l / 146097;
        l = l - (146097 * nn + 3) / 4;
        const long i = 4000 * (l + 1) / 1461001;
        l = l - 1461 * i / 4 + 31;
        const long j = 80 * l / 2447, d = l - 2447 * j / 80;
        l = j / 11;
        const long m = j + 2 - 12 * l, yy = 100 * (nn - 49) + i + l;
        char s[48];
        std::snprintf(s, sizeof s, "%04ld-%02ld-%02ld", yy, m, d);
        return std::string(s);
      };
      const std::string when = date(ue::mjd(y) + day);
      const bool fell = !S.is_vacant(c.from) && S.is_vacant(c.to), let = S.is_vacant(c.from) && !S.is_vacant(c.to);
      events.push(Json::array({"simulated-" + std::to_string(k++), "alteration", fell ? "fell vacant" : let ? "let" : "use changed", when, when,
                               Json::object({{"property", D.ids[c.building]}}), "urban::occupancy", "simulated",
                               Json::object({{"from", S.U.uses[c.from]}, {"to", S.U.uses[c.to]}})}));
    }
    write_json(out + "/simulated-record-district.json",
               Json::object({{"product", "city-events"}, {"id", "melbourne-simulated"}, {"provenance", "simulated"},
                             {"note", "A history simulated by urban-science's occupancy at the calibrated parameters (research/experiments/2026-11-05-owners-and-occupiers): "
                                      "never an observed record. Its censuses and events carry the provenance simulated."},
                             {"model", Json::object({{"lambda", numbers(Cm.parameters.lambda)}, {"vacate", Cm.parameters.vacate}, {"let", Cm.parameters.let},
                                                     {"returning", Cm.parameters.returning},
                                                     {"beta_neighbours", Cm.parameters.beta_neighbours},
                                                     {"beta_demand", Cm.parameters.beta_demand}, {"seed", static_cast<double>(seed)}})},
                             {"window", Json::object({{"from", "2002-01-01"}, {"to", "2025-01-01"}})},
                             {"properties", Json::object({{"columns", Json::array({"id", "at"})}, {"rows", props}})},
                             {"censuses", Json::object({{"columns", Json::array({"year", "property", "predominant_use"})}, {"rows", cens}})},
                             {"events", Json::object({{"columns", Json::array({"id", "type", "kind", "from", "to", "match", "source", "provenance", "detail"})}, {"rows", events}})},
                             {"checked", Json::object({{"problems", static_cast<double>(bad.size())}, {"record_events", static_cast<double>(Rec.events.size())}})}}));
    std::printf("simulated record: %zu properties of the district, %zu changes, the whole record's problems: %zu\n", kept, k, bad.size());
  }
  std::printf("done in %.0f s\n", seconds_since(t_start));
  return 0;
}
