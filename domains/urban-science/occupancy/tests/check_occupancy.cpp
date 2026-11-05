// Checks of the owners and occupiers (occupancy/model.hpp, calibration.hpp, planning.hpp; theory/occupancy.md) on a
// district made for the purpose: 300 buildings on a square lattice, seven of CLUE's space uses in four families.
//
// - The baseline: with the channel of vacancy shut and the modulation off, each building's generator is physis's
//   nested chain, entry by entry, and
//   the agent-based model's one-year law from each use is the chain's kernel (Monte Carlo, four standard errors).
// - Three lumps: with the modulation on, the chain is strongly lumpable onto occupied, lettable vacancy and works.
// - The simulated record: sound by the record's axioms, every event of provenance "simulated", its complex built.
// - Reproducibility: one seed, one history; two seeds, two; an ensemble on four threads equal to one on none.
// - The panel's transitions equal the events module's.
// - Indirect inference recovers known parameters (the channel of vacancy open, the modulation on) from four simulated
//   records: every estimate within four standard errors, each mean within three; the chain alone fits worse.
// - The planner: each building's model a POMDP; the Lagrangian bound above the product's value on three buildings
//   (HSVI on the product: brute force), the ranked policy within its budget and below the bound.

#include <urban/occupancy/calibration.hpp>
#include <urban/occupancy/model.hpp>
#include <urban/occupancy/planning.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>
#include <tyche/decision/point_based.hpp>
#include <urban/events/complex.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

namespace uo = urban::occupancy;
namespace ue = urban::events;
namespace sm = physis::statistical_mechanics;
namespace td = tyche::decision;

namespace {

std::string fmt(const char* f, double a, double b = 0, double c = 0, double d = 0) {
  char s[256];
  std::snprintf(s, sizeof s, f, a, b, c, d);
  return s;
}

uo::Uses the_uses() {
  const std::map<std::string, double> counts = {{"Office", 30},
                                                {"Retail - Shop", 15},
                                                {"Retail - Showroom", 5},
                                                {"Residential Apartment", 25},
                                                {"Unoccupied - Unused", 15},
                                                {"Unoccupied - Under Renovation", 6},
                                                {"Unoccupied - Under Construction", 4}};
  return uo::uses_of(ue::prefix_classification(counts));
}

// Like attracts like: +1 between uses of one family, −0.4 otherwise (a stand-in for ln CLQ).
std::vector<double> like_affinity(const uo::Uses& S) {
  const std::size_t N = S.states();
  std::vector<double> a(N * N);
  for (std::size_t u = 0; u < N; u++)
    for (std::size_t y = 0; y < N; y++) a[u * N + y] = S.family[u] == S.family[y] ? 1.0 : -0.4;
  return a;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const auto S = the_uses();
  const std::size_t N = S.states(), L = S.levels();
  s.require("the uses: seven in four families, two levels of memory, the vacant family lumped into lettable and works", N == 7 && L == 2 &&
                std::count(S.lump.begin(), S.lump.end(), uo::vacant) == 1 && std::count(S.lump.begin(), S.lump.end(), uo::works) == 2);

  // The district: 300 buildings on a 20 × 15 lattice 30 m apart, their six nearest neighbours.
  std::vector<std::array<double, 2>> at;
  for (int r = 0; r < 15; r++)
    for (int c = 0; c < 20; c++) at.push_back({30.0 * c, 30.0 * r});
  const auto X = uo::nearest(at, 6);
  const std::size_t n = at.size();
  calx::random::Stream init(20261105, 0x696e6974u, 0, 0);  // "init"
  std::vector<std::size_t> start_use(n);
  for (auto& u : start_use) u = calx::random::categorical(init, S.U.C.law);
  const auto affinity = like_affinity(S);

  uo::State s0;
  s0.use = start_use, s0.alive.assign(n, 1), s0.enrolled.assign(n, 0), s0.demand = 0.7;
  s0.held.assign(n, 0);
  for (std::size_t i = 0; i < n; i++) s0.held[i] = S.lump[start_use[i]] == uo::occupied ? S.family[start_use[i]] : N;

  // ---- the baseline ---------------------------------------------------------------------------------------------------
  const std::vector<double> lambda = {0.10, 0.15};
  {
    uo::Model M{&S, &X, affinity, {lambda, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};
    const auto G = sm::NestedChain(S.U.C, lambda).generator();
    double worst = 0;
    for (std::size_t i = 0; i < n; i += 37) {
      const auto Gi = uo::generator(M, s0, i);
      for (std::size_t x = 0; x < N; x++)
        for (std::size_t y = 0; y < N; y++) worst = std::max(worst, std::fabs(Gi(x, y) - G(x, y)));
    }
    s.check("the modulation off (ν = η = β = 0): each building's generator is physis's nested chain (max |G − Σ λ_l (Π_l − I)|)", worst, 1e-15);

    // One building's one-year law from each use, by Gillespie, against the chain's kernel.
    const sm::NestedChain chain(S.U.C, lambda);
    uo::Space lone;
    lone.at = {{0.0, 0.0}}, lone.neighbours = {std::vector<std::size_t>{}};
    uo::Model M1{&S, &lone, {}, {lambda, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};
    double z = 0;
    const int runs = 20000;
    for (std::size_t x = 0; x < N; x++) {
      std::vector<double> count(N, 0.0);
      for (int k = 0; k < runs; k++) {
        uo::State one;
        one.use = {x}, one.alive = {1}, one.enrolled = {0}, one.held = {N};
        nomos::agents::gillespie(M1, one, 0.0, 1.0, uo::year_seed(1000 + x, k));
        count[one.use[0]] += 1;
      }
      for (std::size_t y = 0; y < N; y++) {
        const double p = chain.kernel(x, y, 1.0), se = std::sqrt(p * (1 - p) / runs);
        if (se > 0) z = std::max(z, std::fabs(count[y] / runs - p) / se);
      }
    }
    s.require("... and the agent-based model's one-year law from every use is the chain's kernel (20,000 runs each; max |z| < 4)", z < 4, fmt("max |z| %.2f", z));
  }

  // ---- three lumps ----------------------------------------------------------------------------------------------------
  {
    uo::Model M{&S, &X, affinity, {lambda, 0.03, 0.2, 0.5, 1.2, 0.8, 1.0}};
    double worst = 0;
    for (std::size_t i = 0; i < n; i += 23)
      for (const bool enrolled : {false, true}) {
        auto st = s0;
        st.enrolled[i] = enrolled;
        const auto G = uo::generator(M, st, i);
        for (std::size_t A = 0; A < uo::lumps; A++)
          for (std::size_t B = 0; B < uo::lumps; B++) {
            double first = std::nan("");
            for (std::size_t x = 0; x < N; x++) {
              if (S.lump[x] != A) continue;
              double sum = 0;
              for (std::size_t y = 0; y < N; y++)
                if (S.lump[y] == B) sum += G(x, y);
              if (std::isnan(first)) first = sum;
              else worst = std::max(worst, std::fabs(sum - first));
            }
          }
      }
    s.check("the modulation on, the channel of vacancy open: the model is strongly lumpable onto occupied, lettable vacancy and works (max spread of a lump's rates)", worst, 1e-13);
  }

  // ---- runs, records and reproducibility -------------------------------------------------------------------------------
  const int from = 2002, to = 2022;
  std::vector<uo::Life> lives(n);
  for (std::size_t i = 0; i < n; i++) {
    lives[i].first = from + (i % 10 == 0 ? 5 : 0), lives[i].last = to - (i % 17 == 0 ? 3 : 0), lives[i].initial = start_use[i];
    for (int y = lives[i].first; y <= lives[i].last; y++)
      if (!(i % 7 == 3 && y % 2 == 1) || y == lives[i].first || y == lives[i].last) lives[i].censuses.push_back(y);
  }
  std::vector<double> demand;
  for (int y = from; y <= to; y++) demand.push_back(std::sin(0.6 * (y - from)) + 0.03 * (y - from) - 0.3);
  const uo::Parameters truth{lambda, 0.03, 0.2, 0.5, 1.0, 0.5, 1.0};
  uo::Model Mt{&S, &X, affinity, truth};
  const auto H = uo::simulate(Mt, lives, demand, from, to, 777);
  {
    std::vector<std::string> ids;
    for (std::size_t i = 0; i < n; i++) ids.push_back("b" + std::to_string(i));
    const auto R = uo::to_record(S, X, ids, H, 777);
    const auto bad = ue::problems(R);
    std::size_t observations = 0, alterations = 0, changes = 0, censuses = 0, simulated = 0;
    for (const auto& e : R.events) {
      observations += e.type == "observation", alterations += e.type == "alteration", simulated += e.provenance == "simulated";
    }
    for (const auto& c : H.changes) changes += c.change.from != c.change.to;
    for (const auto& l : lives) censuses += l.censuses.size();
    s.require("the simulated record is sound by the record's axioms", bad.empty(), bad.empty() ? "" : bad.front());
    s.require("... every event of provenance \"simulated\", one observation a census and one alteration a change of use",
              simulated == R.events.size() && observations == censuses && alterations == changes,
              std::to_string(observations) + " observations, " + std::to_string(alterations) + " alterations");
    bool built = true;
    std::size_t vertices = 0;
    try {
      const auto K = ue::construct(R);
      vertices = K.K.cw.vertices;
    } catch (const std::exception&) {
      built = false;
    }
    s.require("... and its event complex is built (K's vertices: the events and the window's two ends)", built && vertices == R.events.size() + 2, std::to_string(vertices));
  }
  {
    const auto H2 = uo::simulate(Mt, lives, demand, from, to, 777), H3 = uo::simulate(Mt, lives, demand, from, to, 778);
    bool same = H2.changes.size() == H.changes.size(), differ = H3.changes.size() != H.changes.size();
    for (std::size_t k = 0; same && k < H.changes.size(); k++)
      same = H.changes[k].time == H2.changes[k].time && H.changes[k].change.to == H2.changes[k].change.to && H.changes[k].change.building == H2.changes[k].change.building;
    for (std::size_t k = 0; !differ && k < H.changes.size(); k++) differ = H.changes[k].time != H3.changes[k].time;
    s.require("reproducibility: one seed, one history to the bit; another seed, another", same && differ, std::to_string(H.changes.size()) + " redraws");
    uo::Calibrand D{&S, &X, affinity, lives, demand, from, to, 1.0};
    calx::scheduler::ThreadPool pool(4);
    const auto th = uo::theta_of(truth);
    const auto E1 = uo::ensemble(D, th, 6, 31, lambda), E4 = uo::ensemble(D, th, 6, 31, lambda, &pool);
    bool equal = true;
    for (std::size_t r = 0; r < 6; r++) equal = equal && E1.psi[r] == E4.psi[r];
    s.require("... and an ensemble's statistics on four threads equal those on one", equal);
  }

  // ---- the panel ------------------------------------------------------------------------------------------------------
  {
    std::vector<ue::Observation> O;
    for (std::size_t i = 0; i < n; i++)
      for (const auto& [y, u] : H.seen[i]) O.push_back({"b" + std::to_string(i), S.U.uses[u], y});
    auto a = uo::transitions(H.seen), b = ue::transitions_of(S.U, O);
    const auto key = [](const sm::Transition& t) { return std::make_tuple(t.x, t.y, t.tau, t.count); };
    std::sort(a.begin(), a.end(), [&](auto& p, auto& q) { return key(p) < key(q); });
    std::sort(b.begin(), b.end(), [&](auto& p, auto& q) { return key(p) < key(q); });
    bool equal = a.size() == b.size();
    for (std::size_t k = 0; equal && k < a.size(); k++) equal = key(a[k]) == key(b[k]);
    s.require("the panel's transitions equal the events module's (transitions_of)", equal, std::to_string(a.size()) + " (x, y, τ)");
  }

  // ---- indirect inference ---------------------------------------------------------------------------------------------
  // Four records simulated at the truth (seeds 777 … 780), each calibrated from its own statistics on its own seeds:
  // every estimate within four of its standard errors of the truth, and each parameter's mean over the four within three
  // standard errors of that mean (their spread over √4) of the truth; the chain alone fits the first record far worse.
  {
    const auto t0 = std::chrono::steady_clock::now();
    uo::Calibrand D{&S, &X, affinity, lives, demand, from, to, 1.0};
    calx::scheduler::ThreadPool pool(4);
    const auto want = uo::theta_of(truth);
    const std::size_t p = want.size();
    std::vector<std::vector<double>> est;
    double zmax = 0, j_first = 0, j_chain = 0;
    std::string got;
    bool finite = true;
    for (int k = 0; k < 4; k++) {
      const auto Hk = uo::simulate(Mt, lives, demand, from, to, 777 + static_cast<std::uint64_t>(k));
      const auto A0 = uo::auxiliary(S, X, affinity, Hk.seen, demand, from, lambda);
      const auto start = uo::starting_parameters(S, A0, 1.0);
      uo::CalibrationOptions o;
      o.replicates = 12, o.seed = 5000 + 1000 * static_cast<std::uint64_t>(k), o.pool = &pool;
      const auto C = uo::calibrate(D, Hk.seen, start, o);
      est.push_back(C.theta);
      for (std::size_t j = 0; j < p; j++) {
        finite = finite && std::isfinite(C.error[j]) && C.error[j] > 0;
        zmax = std::max(zmax, std::fabs(C.theta[j] - want[j]) / C.error[j]);
      }
      if (k == 0) {
        j_first = C.objective;
        o.vacancy = false, o.neighbours = false, o.demand = false;
        j_chain = uo::calibrate(D, Hk.seen, start, o).objective;
      }
    }
    bool centred = true;
    for (std::size_t j = 0; j < p; j++) {
      double m = 0, v = 0;
      for (const auto& e : est) m += e[j] / 4;
      for (const auto& e : est) v += (e[j] - m) * (e[j] - m) / 3;
      const double se = std::sqrt(v / 4);
      centred = centred && std::fabs(m - want[j]) <= 3 * se;
      got += fmt("%.3f ± %.3f (true %.3f); ", m, se, want[j]);
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    s.require("indirect inference recovers (ln λ_0, ln λ_1, ln ν, ln η, logit ρ, β_n, β_d) from four simulated records: every estimate within four standard errors",
              finite && zmax <= 4, fmt("max |z| %.2f", zmax));
    s.require("... and each parameter's mean over the four within three standard errors of the truth", centred, got + fmt("%.0f s", secs));
    s.require("... and the chain alone (ν = η = β = 0) fits the first record far worse: its J far above the model's", j_chain > j_first + 9,
              fmt("J %.2f against %.2f", j_chain, j_first));
  }

  // ---- the planner ----------------------------------------------------------------------------------------------------
  {
    uo::Model M{&S, &X, affinity, {lambda, 0.03, 0.2, 0.5, 1.2, 0.8, 1.0}};
    const uo::PlanAssumptions A;
    bool valid = true;
    double rewards_ok = 1;
    for (std::size_t i = 0; i < n; i += 41) {
      const auto P = uo::building_model(M, s0, i, A);
      try {
        td::validate(P);
      } catch (const std::exception&) {
        valid = false;
      }
      for (std::size_t a = 0; a < 3; a++) rewards_ok = std::min(rewards_ok, P.reward(uo::occupied, a) > P.reward(uo::vacant, a) && P.reward(uo::vacant, a) >= -1 ? 1.0 : 0.0);
      rewards_ok = std::min(rewards_ok, P.reward(uo::vacant, uo::incentive) > P.reward(uo::vacant, uo::nothing) ? 1.0 : 0.0);
    }
    s.require("each building's model is a POMDP; a year lettable-vacant costs at most 1, less when enrolled", valid && rewards_ok == 1);

    // Three buildings of different neighbourhoods (one in a vacant street), planned by the bound and by brute force.
    auto st = s0;
    for (const std::size_t j : X.neighbours[0]) st.use[j] = S.lettable;
    std::vector<char> planned(n, 0);
    const std::size_t pick[3] = {0, 150, 299};
    std::vector<td::Belief> beliefs(n, td::Belief(uo::lumps, 1.0 / uo::lumps));
    const std::vector<td::Belief> b3 = {{0.1, 0.8, 0.1}, {0.5, 0.4, 0.1}, {0.3, 0.3, 0.4}};
    for (int k = 0; k < 3; k++) planned[pick[k]] = 1, beliefs[pick[k]] = b3[static_cast<std::size_t>(k)];
    const double budget = 1.0;
    const auto C = uo::plan_city(M, st, planned, beliefs, budget, A, {.epsilon = 0.02, .max_trials = 400});
    const auto J = td::joint_pomdp(C.components, budget);
    const auto joint = td::hsvi(J, {.epsilon = 0.2});
    s.require("three buildings: the product (27 states, its joint actions within the budget) solved by HSVI", joint.converged,
              fmt("[%.4f, %.4f], %.0f joint actions", joint.lower_value, joint.upper_value, static_cast<double>(J.actions)));
    s.require("the Lagrangian bound lies above V* (above the product's lower bound), and doing nothing below it",
              C.solution.upper >= joint.lower_value - 1e-9 && C.nothing_value <= joint.upper_value + 1e-9,
              fmt("bound %.4f at λ = %.4f; product [%.4f, %.4f]", C.solution.upper, C.solution.lambda, joint.lower_value, joint.upper_value) +
                  fmt("; nothing %.4f", C.nothing_value));
    uo::Knowledge K;
    K.lump.assign(n, -1), K.year.assign(n, 2024), K.since.assign(n, 2024), K.belief = beliefs;
    const auto E = uo::evaluate_decoupled(C, uo::Rule::pomdp, X, K, 2025, 150, 400, 99);
    s.require("the ranked policy keeps its budget and its value lies below the bound, within 10 % of the product's lower bound",
              E.spent <= budget + 1e-12 && E.mean <= C.solution.upper + 4 * E.error && E.mean >= joint.lower_value - 0.1 * std::fabs(joint.lower_value) - 4 * E.error - E.truncation,
              fmt("%.4f ± %.4f (spent %.2f a year; truncation %.1e)", E.mean, E.error, E.spent, E.truncation));
  }
  return s.done();
}
