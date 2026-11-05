// Checks of the walkers (movement/*.hpp) on a small grid district and two corridors, against the flow study's law,
// exact expectations, closed forms and an exactly solved POMDP.
//
// - The network: its directed segments, crossings and phases.
// - The route law is the district-economy flow study's: its transitions, and its walk sums and visits computed as that
//   study computes them (Z = (I − M)⁻¹, the absorbing destination by Woodbury), against the module's Gauss–Seidel.
// - The walkers' law is a Markov chain whose mean visits are the logit flow exactly (its own fundamental matrix).
// - The walkers' mean against the logit flow: entries by departure hour over replicate runs, in the exact standard
//   errors of a compound Poisson count (the second moments of the visits from the chain's fundamental matrix).
// - Conservation: walkers born = arrived + failed + walking; at every junction the walks in and out balance.
// - Reproducibility: one thread and four give the same run, bit for bit; another seed another run.
// - Calibration recovers the amplitudes of exact counts.
// - A narrow segment passes no more than its capacity and a queue forms at its mouth.
// - A signal's mean wait for walkers arriving at random: (1 − g)² C / 2.
// - Regimes: closed to walkers passing through, one-way, and the street's visitors still reaching it.
// - The record of a run: sound, every event simulated.
// - The corridor: its agents have the law of its POMDP; POMCP on the agents earns tyche's exact interval.

#include <urban/movement/movement.hpp>

#include <calx/linalg/dense.hpp>
#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>
#include <tyche/decision/alpha.hpp>
#include <tyche/decision/online.hpp>
#include <tyche/decision/point_based.hpp>
#include <urban/events/record.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace um = urban::movement;
namespace la = calx::linalg;
namespace td = tyche::decision;

namespace {

std::string num(double x, int digits = 3) {
  char b[64];
  std::snprintf(b, sizeof b, "%.*g", digits, x);
  return b;
}

// A grid of nx × ny blocks of 80 m × 60 m, its streets named as the Hoddle Grid's (so that their widths are the
// table's), and a dead-end lane into one block.
urban::landscape::Network grid() {
  urban::landscape::Network N;
  const std::size_t nx = 3, ny = 2;
  const double bx = 80, by = 60;
  const char* rows[] = {"Flinders Lane", "Collins Street", "Little Collins Street"};
  const char* cols[] = {"Swanston Street", "Russell Place", "Russell Street", "Exhibition Street"};
  for (std::size_t j = 0; j <= ny; j++)
    for (std::size_t i = 0; i <= nx; i++) N.junctions.push_back({bx * static_cast<double>(i), by * static_cast<double>(j)});
  const auto v = [&](std::size_t i, std::size_t j) { return j * (nx + 1) + i; };
  for (std::size_t j = 0; j <= ny; j++)
    for (std::size_t i = 0; i < nx; i++) N.segments.push_back({v(i, j), v(i + 1, j), {}, "street", rows[j]});
  for (std::size_t i = 0; i <= nx; i++)
    for (std::size_t j = 0; j < ny; j++) N.segments.push_back({v(i, j), v(i, j + 1), {}, "street", cols[i]});
  N.junctions.push_back({100, 25});
  N.segments.push_back({v(1, 0), N.junctions.size() - 1, {}, "laneway", "a lane"});
  return N;
}

// The flow study's transitions (its flow.hpp, transitions()), written again here: M_ab = ω_b for b leaving the
// junction a enters, b not a reversed unless that junction is a dead end.
la::Matrix study_transitions(const um::Walkway& W, double theta, double wref) {
  const std::size_t S = W.states();
  std::vector<std::size_t> degree(W.junctions(), 0);
  for (const auto& e : W.segments) degree[e.u]++, degree[e.v]++;
  la::Matrix M(S);
  for (std::size_t a = 0; a < S; a++)
    for (const std::size_t b : W.leaving[W.head(a)]) {
      if (b == (a ^ 1u) && degree[W.head(a)] != 1) continue;
      const auto& e = W.segments[b / 2];
      M(a, b) = std::exp(-theta * e.length) * std::pow(e.width / wref, e.length / e.parent);
    }
  return M;
}

// The walkers' chain toward a target: P_ab = M_ab z(b)/z(a), and its start from an origin; the mean visits and their
// second moments from its fundamental matrix N = (I − P_TT)⁻¹ (an absorbing chain: E V_c² = Σ_i p_i N_ic (2 N_cc − 1)).
struct ChainMoments {
  std::vector<double> mean, second;
};
ChainMoments chain_moments(const um::Routes& R, std::size_t origin, const std::vector<std::size_t>& target, const std::vector<double>& z) {
  const std::size_t S = R.states();
  std::vector<char> in(S, 0);
  for (auto s : target) in[s] = 1;
  std::vector<std::size_t> T;
  std::vector<std::size_t> pos(S, um::none);
  for (std::size_t s = 0; s < S; s++)
    if (!in[s] && z[s] > 0) pos[s] = T.size(), T.push_back(s);
  const std::size_t n = T.size();
  la::Matrix A(n);
  for (std::size_t i = 0; i < n; i++) A(i, i) = 1;
  for (std::size_t i = 0; i < n; i++)
    for (const auto& law : R.next_law(T[i], z))
      if (pos[law.first] != um::none) A(i, pos[law.first]) -= law.second;
  const la::Matrix N = la::solve(A, la::Matrix::identity(n));
  std::vector<double> p0(S, 0.0);
  for (const auto& [b, p] : R.first_law(origin, z)) p0[b] = p;
  ChainMoments m;
  m.mean.assign(S, 0.0), m.second.assign(S, 0.0);
  for (std::size_t j = 0; j < n; j++) {
    double v = 0;
    for (std::size_t i = 0; i < n; i++) v += p0[T[i]] * N(i, j);
    m.mean[T[j]] = v, m.second[T[j]] = v * (2 * N(j, j) - 1);
  }
  // Absorbing states: entered once at most.
  for (std::size_t c = 0; c < S; c++)
    if (in[c]) {
      double v = p0[c];
      for (std::size_t i = 0; i < n; i++)
        for (const auto& law : R.next_law(T[i], z))
          if (law.first == c) v += m.mean[T[i]] * law.second;
      m.mean[c] = v, m.second[c] = v;
    }
  return m;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  auto W = std::make_shared<um::Walkway>(um::walkway(grid()));
  const std::size_t S = W->states();
  s.require("the grid district: 13 junctions, 18 segments, 36 directed segments", W->junctions() == 13 && W->edges() == 18 && S == 36,
            std::to_string(W->junctions()) + ", " + std::to_string(W->edges()) + ", " + std::to_string(S));
  {
    // Straight on along Collins Street through Russell Place: one carriageway crossed; turning there: one half; the
    // lane's mouth, a junction of three carriageways and a lane, is crossed straight on.
    const std::size_t a = 2 * 3, b = 2 * 4;   // (0,1)→(1,1), (1,1)→(2,1)
    std::size_t up = um::none;
    for (auto t : W->leaving[W->head(a)])
      if (W->segments[t / 2].name == "Russell Place" && W->head(t) != W->tail(a) && W->at[W->head(t)][1] > 60) up = t;
    s.require("crossings: straight on 1, a turn ½, at a junction of four carriageways; the phase of a turn is the other",
              W->crossings(a, b) == 1 && up != um::none && W->crossings(a, up) == 0.5 && W->phase(b) != W->phase(up),
              num(W->crossings(a, b)) + ", " + num(W->crossings(a, up)));
  }

  // ---- the law is the flow study's -----------------------------------------------------------------------------------
  const um::RouteLaw law;   // θ = 0.02 per metre, w_ref = 12 m
  auto R = std::make_shared<const um::Routes>(W, law);
  um::Targets T;
  {
    const la::Matrix M = study_transitions(*W, law.theta, law.reference_width);
    double dm = 0;
    for (std::size_t a = 0; a < S; a++) {
      std::vector<double> row(S, 0.0);
      for (const auto& arc : R->out(a)) row[arc.to] = arc.w;
      for (std::size_t b = 0; b < S; b++) dm = std::max(dm, std::fabs(row[b] - M(a, b)));
    }
    s.check("the transitions are the flow study's (max |ΔM|)", dm, 1e-15);
    const double radius = R->spectral_radius();
    s.require("the spectral radius of the transitions is below one: the walk sums converge", radius < 1, num(radius, 4));
    la::Matrix A(S);
    for (std::size_t i = 0; i < S; i++)
      for (std::size_t j = 0; j < S; j++) A(i, j) = (i == j ? 1.0 : 0.0) - M(i, j);
    const la::Matrix Z = la::solve(A, la::Matrix::identity(S));
    double dz = 0, df = 0, top = 0;
    for (std::size_t d = 0; d < 13; d++) {
      const std::size_t t = T.junction(*W, d);
      const auto& Ad = T.states(t);
      const std::size_t r = Ad.size();
      la::Matrix B(r);
      for (std::size_t i = 0; i < r; i++)
        for (std::size_t j = 0; j < r; j++) B(i, j) = Z(Ad[i], Ad[j]);
      const la::Matrix Q = la::solve(B, la::Matrix::identity(r));
      std::vector<double> q(r, 0.0), sigma(S, 0.0);
      for (std::size_t i = 0; i < r; i++)
        for (std::size_t j = 0; j < r; j++) q[i] += Q(i, j);
      for (std::size_t c = 0; c < S; c++)
        for (std::size_t i = 0; i < r; i++) sigma[c] += Z(c, Ad[i]) * q[i];
      const auto z = R->value(Ad);
      for (std::size_t c = 0; c < S; c++) dz = std::max(dz, std::fabs(z[c] - sigma[c]) / std::max(1.0, sigma[c]));
      for (std::size_t o = 0; o < 13; o++) {
        if (o == d) continue;
        // π_od(c) = Y[o, c] − (Y[o, A] Q)·(Z[A, c] − δ), Y[o, c] = Σ_{b out of o} ω_b Z[b, c]; f = π σ / W.
        const auto Y = [&](std::size_t c) { double v = 0; for (auto b : W->leaving[o]) v += R->weight(b) * Z(b, c); return v; };
        std::vector<double> yq(r, 0.0);
        for (std::size_t j = 0; j < r; j++)
          for (std::size_t i = 0; i < r; i++) yq[j] += Y(Ad[i]) * Q(i, j);
        double Wod = 0;
        for (double v : yq) Wod += v;
        const auto f = R->visits(o, Ad, z);
        for (std::size_t c = 0; c < S; c++) {
          double pi = Y(c);
          for (std::size_t i = 0; i < r; i++) pi -= yq[i] * (Z(Ad[i], c) - (Ad[i] == c ? 1.0 : 0.0));
          const double fs = pi * sigma[c] / Wod;
          df = std::max(df, std::fabs(f[c] - fs)), top = std::max(top, fs);
        }
      }
    }
    s.check("the walk sums z by Gauss–Seidel against the flow study's σ_d = Z[:, A] Z[A, A]⁻¹ 1, 13 destinations, the lane's dead end among them (relative)", dz, 1e-12);
    s.check("a trip's visits against the flow study's π_od σ_d / W_od, 156 trips (largest error over largest visit)", df / top, 1e-12);
  }

  // ---- a long ladder: walk sums of 10⁻¹⁴, held state by state ----------------------------------------------------------
  {
    urban::landscape::Network N;
    const std::size_t n = 41;
    for (std::size_t row = 0; row < 2; row++)
      for (std::size_t i = 0; i < n; i++) N.junctions.push_back({40.0 * static_cast<double>(i), 10.0 * static_cast<double>(row)});
    for (std::size_t i = 0; i + 1 < n; i++) N.segments.push_back({i, i + 1, {}, "street", "Collins Street"}), N.segments.push_back({n + i, n + i + 1, {}, "street", "Flinders Lane"});
    for (std::size_t i = 0; i < n; i++) N.segments.push_back({i, n + i, {}, "laneway", "a lane"});
    auto Lw = std::make_shared<um::Walkway>(um::walkway(N));
    const um::Routes RL(Lw, law);
    const std::size_t SL = Lw->states(), far = 2 * n - 1;
    const auto& A = Lw->entering[far];
    const auto z = RL.value(A);
    // The flow study's formulas on a dense inverse.
    la::Matrix I_M(SL);
    for (std::size_t a = 0; a < SL; a++) {
      I_M(a, a) += 1;
      for (const auto& arc : RL.out(a)) I_M(a, arc.to) -= arc.w;
    }
    const la::Matrix Z = la::solve(I_M, la::Matrix::identity(SL));
    const std::size_t r = A.size();
    la::Matrix B(r);
    for (std::size_t i = 0; i < r; i++)
      for (std::size_t j = 0; j < r; j++) B(i, j) = Z(A[i], A[j]);
    const la::Matrix Q = la::solve(B, la::Matrix::identity(r));
    std::vector<double> q(r, 0.0);
    for (std::size_t i = 0; i < r; i++)
      for (std::size_t j = 0; j < r; j++) q[i] += Q(i, j);
    double dz = 0;
    for (std::size_t c = 0; c < SL; c++) {
      double sg = 0;
      for (std::size_t i = 0; i < r; i++) sg += Z(c, A[i]) * q[i];
      if (sg > 0) dz = std::max(dz, std::fabs(z[c] - sg) / sg);
    }
    const auto f = RL.visits(0, A, z);
    double into = 0;
    for (const auto c : A) into += f[c];
    s.check("a ladder 1.6 km long (W = " + num(RL.total(0, z), 3) + "): the walk sums against the dense formulas, each state relative", dz, 1e-10);
    s.check("… and a trip along it enters its target exactly once (Σ of its visits over the target − 1)", into - 1, 1e-10);
  }

  // ---- the walkers' chain: its mean visits are the logit flow ------------------------------------------------------------
  {
    double worst = 0;
    for (std::size_t d : {0u, 5u, 11u, 12u}) {
      const std::size_t t = T.junction(*W, d);
      const auto z = R->value(T.states(t));
      for (std::size_t o : {3u, 6u, 12u}) {
        const auto m = chain_moments(*R, o, T.states(t), z);
        const auto f = R->visits(o, T.states(t), z);
        for (std::size_t c = 0; c < S; c++) worst = std::max(worst, std::fabs(m.mean[c] - f[c]));
      }
    }
    s.check("the walkers' Markov chain (P_ab = M_ab z_b / z_a) has the logit flow's mean visits (its fundamental matrix)", worst, 1e-12);
  }

  // ---- the walkers' mean against the logit flow ---------------------------------------------------------------------------
  const std::size_t lane_end = 12;
  std::vector<um::Stream> streams(3);
  streams[0] = {"across", 0, 0, {0, 3, 8}, {T.junction(*W, 11), T.junction(*W, 8), T.junction(*W, 0)}, {}, {1, 2, 1}, {700, 900, 500, 0, 0}, 0};
  streams[1] = {"to the lane", 1, 0, {0, 11}, {T.junction(*W, lane_end), T.junction(*W, lane_end)}, {}, {1, 1}, {300, 300, 300, 0, 0}, 0};
  streams[2] = {"from the lane", 2, 0, {lane_end}, {T.junction(*W, 7)}, {}, {1}, {400, 0, 400, 0, 0}, 0};
  auto V = std::make_shared<const um::Values>(R, T);
  auto Tp = std::make_shared<const um::Targets>(T);
  um::Options opt;
  opt.hours = 6;
  opt.kinematics.capacity = false;
  opt.kinematics.signals = false;
  std::vector<um::Counter> counters{{"Collins Street east of Russell Place", 4, 40, 0.5, {0.1, -0.1}, {8, 9}},
                                    {"Russell Place south", 11, 30, 0.5, {0, 0}, {22, 23}}};
  {
    // The expectation and the exact variance of each state's entries by departure hour over the run.
    std::vector<double> mu(S, 0.0), var(S, 0.0);
    for (const auto& st : streams) {
      double rate = 0, wsum = 0;
      for (double r : st.rate) rate += r;
      for (double w : st.weight) wsum += w;
      for (std::size_t p = 0; p < st.origin.size(); p++) {
        const auto m = chain_moments(*R, st.origin[p], T.states(st.target[p]), (*V)[st.target[p]]);
        const double lam = rate * st.weight[p] / wsum;
        for (std::size_t c = 0; c < S; c++) mu[c] += lam * m.mean[c], var[c] += lam * m.second[c];
      }
    }
    const int runs = 6;
    std::vector<double> sum(S, 0.0);
    double chi = 0, worst = 0;
    std::size_t used = 0;
    for (int k = 0; k < runs; k++) {
      um::Options o = opt;
      o.seed = 1000 + static_cast<std::uint64_t>(k);
      um::Simulation sim(Tp, {V}, streams, counters, o);
      sim.run(6 * 3600);
      const auto& t = sim.tallies();
      for (std::size_t h = 0; h < t.hours; h++)
        for (std::size_t c = 0; c < S; c++) sum[c] += t.entries_by_departure[h * S + c];
      if (k == 0) {
        s.require("every walker arrived: born = arrived + failed + walking, none failed, none left", sim.totals().departed == sim.totals().arrived + sim.totals().failed + static_cast<double>(sim.active()) &&
                  sim.totals().failed == 0 && sim.active() == 0, num(sim.totals().departed, 8) + " departed, " + num(sim.totals().arrived, 8) + " arrived");
        // At every junction the walks in and out balance.
        double worst_balance = 0;
        for (std::size_t x = 0; x < W->junctions(); x++) {
          double in = 0, out = 0;
          for (std::size_t h = 0; h < t.hours; h++) {
            for (auto c : W->entering[x]) in += t.entries[h * S + c];
            for (auto c : W->leaving[x]) out += t.entries[h * S + c];
          }
          worst_balance = std::max(worst_balance, std::fabs(in - t.ends[x] + t.starts[x] - out));
        }
        s.check("conservation at every junction: entries in − walks ended + walks begun = entries out", worst_balance, 0);
      }
    }
    for (std::size_t c = 0; c < S; c++) {
      if (mu[c] < 20) continue;
      const double z = (sum[c] / runs - mu[c]) / std::sqrt(var[c] / runs);
      chi += z * z, worst = std::max(worst, std::fabs(z)), used++;
    }
    s.check("the walkers' mean entries against the logit flow: the largest of " + std::to_string(used) + " standardised deviations (exact compound-Poisson errors)", worst, 4.5);
    s.check("… and their χ² over its degrees of freedom, as a normal deviate", (chi - static_cast<double>(used)) / std::sqrt(2.0 * static_cast<double>(used)), 3.5);
  }

  // ---- reproducibility on any number of threads -----------------------------------------------------------------------
  {
    um::Options o = opt;
    o.kinematics.capacity = true, o.kinematics.signals = true;
    for (std::size_t x = 0; x < W->junctions(); x++) W->signalled[x] = W->carriageways[x] >= 3;
    um::Simulation one(Tp, {V}, streams, counters, o), four = one, other(Tp, {V}, streams, counters, [&] { auto p = o; p.seed++; return p; }());
    calx::scheduler::ThreadPool pool(4);
    one.run(3 * 3600), four.run(3 * 3600, &pool), other.run(3 * 3600);
    const auto& a = one.tallies(), &b = four.tallies(), &c = other.tallies();
    s.require("one thread and four give the same run, bit for bit (entries, counts, occupancy, totals)",
              a.entries == b.entries && a.counts == b.counts && a.occupancy == b.occupancy && one.totals().present == four.totals().present && one.totals().delay == four.totals().delay,
              num(one.totals().arrived, 8) + " trips");
    s.require("another seed, another run", a.entries != c.entries);
  }

  // ---- calibration ----------------------------------------------------------------------------------------------------
  {
    std::vector<um::Stream> base = streams;
    for (auto& st : base) st.rate = {100, 100, 100, 100, 100};
    std::vector<um::UnitCounts> unit;
    for (const auto& st : base) unit.push_back(um::unit_counts(*V, T, st, counters));
    const std::vector<std::vector<double>> truth{{2.0, 7.5, 3.1, 0.0, 1.2}, {0.4, 0.9, 2.2, 1.0, 0.0}, {1.0, 0.0, 0.5, 3.0, 2.0}};
    std::vector<std::vector<std::array<double, 2>>> obs(counters.size(), std::vector<std::array<double, 2>>(5, {0.0, 0.0}));
    for (std::size_t k = 0; k < counters.size(); k++)
      for (std::size_t h = 0; h < 5; h++)
        for (int w = 0; w < 2; w++)
          for (std::size_t g = 0; g < 3; g++)
            obs[k][h][w] += truth[g][h] * 100 * counters[k].share * (1 + counters[k].hand[w]) * (unit[g].out[k][w] + unit[g].back[k][w]);
    const auto C = um::calibrate(base, unit, counters, obs, {{0}, {1}, {2}});
    double err = 0;
    for (std::size_t g = 0; g < 3; g++)
      for (std::size_t h = 0; h < 5; h++) err = std::max(err, std::fabs(C.amplitude[g][h] - truth[g][h]));
    s.check("calibration by NNLS recovers the amplitudes of exact counts, three streams, four rows an hour", err, 1e-8);
  }

  // ---- capacity and queues ----------------------------------------------------------------------------------------------
  {
    urban::landscape::Network N;
    N.junctions = {{0, 0}, {60, 0}, {100, 0}, {160, 0}};
    N.segments = {{0, 1, {}, "street", "Collins Street"}, {1, 2, {}, "laneway", "a narrow lane"}, {2, 3, {}, "street", "Collins Street"}};
    auto L = std::make_shared<um::Walkway>(um::walkway(N));
    auto RL = std::make_shared<const um::Routes>(L, law);
    um::Targets TL;
    const std::size_t to = TL.junction(*L, 3);
    auto VL = std::make_shared<const um::Values>(RL, TL);
    std::vector<um::Stream> st{{"through the lane", 0, 0, {0}, {to}, {}, {1}, {9000, 0, 0}, 0}};
    um::Options o;
    o.hours = 3, o.minute_bins = 60;
    um::Simulation sim(std::make_shared<const um::Targets>(TL), {VL}, st, {}, o);
    sim.set_width(1, 1.0);
    sim.run(3 * 3600);
    const auto& t = sim.tallies();
    double through = 0, upstream = 0;
    for (std::size_t h = 0; h < 3; h++) through = std::max(through, t.entries[h * L->states() + 2]);
    for (std::size_t b = 0; b * L->edges() < t.minutes.size(); b++) upstream = std::max(upstream, static_cast<double>(t.minutes[b * L->edges()]) / 60 / (60 * 9.0));
    const double cap = 1.34 * 5.4 * 1.0 / 4 * 3600;
    s.require("a lane 1 m wide, its demand above capacity, passes its capacity v ρ_j w / 4 in an hour within 5% and no more; the rest queue upstream above the critical density",
              through <= cap + 1 && through >= 0.95 * cap && upstream > 5.4 / 2,
              num(through, 6) + " in its busiest hour against " + num(cap, 6) + " (demand 9000); upstream, a minute's density up to " + num(upstream, 3) + " m⁻²");
    const auto rec = um::simulated_record(sim, {}, {.mjd0 = 61000});
    bool simulated = !rec.events.empty();
    std::size_t crowds = 0, queues = 0;
    for (const auto& e : rec.events) simulated = simulated && e.provenance == "simulated", crowds += e.kind == "crowd", queues += e.kind == "queue";
    const auto problems = urban::events::problems(rec);
    s.require("the run's record is sound, every event simulated, with crowds and queues among them", problems.empty() && simulated && crowds > 0 && queues > 0,
              std::to_string(rec.events.size()) + " events, " + std::to_string(crowds) + " crowds, " + std::to_string(queues) + " queues" + (problems.empty() ? "" : "; " + problems[0]));
  }

  // ---- a signal's mean wait ---------------------------------------------------------------------------------------------
  {
    urban::landscape::Network N;
    N.junctions = {{0, 0}, {100, 0}, {200, 0}, {100, 100}, {100, -100}};
    N.segments = {{0, 1, {}, "street", "Collins Street"}, {1, 2, {}, "street", "Collins Street"}, {1, 3, {}, "street", "Russell Street"}, {1, 4, {}, "street", "Russell Street"}};
    auto X = std::make_shared<um::Walkway>(um::walkway(N));
    X->signalled[1] = 1;
    auto RX = std::make_shared<const um::Routes>(X, um::RouteLaw{.theta = 0.1});   // no walker strays up a side street and back
    um::Targets TX;
    const std::size_t to = TX.junction(*X, 2);
    auto VX = std::make_shared<const um::Values>(RX, TX);
    std::vector<um::Stream> st{{"across Russell Street", 0, 0, {0}, {to}, {}, {1}, std::vector<double>(8, 500), 0}};
    um::Options o;
    o.hours = 9, o.kinematics.capacity = false;
    um::Simulation sim(std::make_shared<const um::Targets>(TX), {VX}, st, {}, o);
    for (std::size_t e = 0; e < 4; e++) sim.set_width(e, 100);   // so wide that the crowd does not slow them
    sim.run(9 * 3600);
    const double n = sim.totals().arrived, wait = sim.totals().delay / n;
    const double sd = std::sqrt(0.5 * 3600.0 / 3 - 15.0 * 15.0);   // wait 0 half the time, else uniform on (0, 60 s)
    const double se = sd / std::sqrt(n);
    s.require("a signal's mean wait for walkers arriving at random is (1 − g)² C / 2 = 15 s, within four standard errors and 0.3 s of rounding",
              std::fabs(wait - 15) <= 4 * se + 0.3, num(wait, 4) + " s ± " + num(se, 2) + " over " + num(n, 5) + " walkers");
  }

  // ---- regimes ---------------------------------------------------------------------------------------------------------
  {
    um::Street street{"Collins Street, Russell Place to Russell Street", {4}, 5};
    um::Targets TR = T;
    const std::size_t on_street = TR.segments(street.edges);
    std::vector<um::Stream> st = streams;
    st[0].rate = {1500, 0, 0}, st[1].rate = {600, 0, 0}, st[2].rate = {0, 0, 0};
    st.push_back({"to the street", 3, 1, {0, 11}, {on_street, on_street}, {T.junction(*W, 0), T.junction(*W, 11)}, {1, 1}, {600, 0, 0}, 300});
    auto TRp = std::make_shared<const um::Targets>(TR);
    auto RR = std::make_shared<const um::RegimeRoutes>(um::regime_routes(W, law, street, TR));
    um::Options o = opt;
    o.hours = 3;
    std::array<std::vector<double>, 3> e;
    std::array<double, 3> visitors{};
    for (int r = 0; r < 3; r++) {
      auto only_through = st;
      only_through[3].rate = {0, 0, 0};
      um::Simulation a(TRp, {RR->values[static_cast<std::size_t>(r)][0], RR->values[static_cast<std::size_t>(r)][1]}, only_through, {}, o);
      a.run(3 * 3600);
      e[static_cast<std::size_t>(r)] = {0, 0};
      for (std::size_t h = 0; h < 3; h++) e[static_cast<std::size_t>(r)][0] += a.tallies().entries[h * S + 8], e[static_cast<std::size_t>(r)][1] += a.tallies().entries[h * S + 9];
      um::Simulation b(TRp, {RR->values[static_cast<std::size_t>(r)][0], RR->values[static_cast<std::size_t>(r)][1]}, st, {}, o);
      b.run(3 * 3600);
      visitors[static_cast<std::size_t>(r)] = b.tallies().trips[3];
      visitors[static_cast<std::size_t>(r)] = b.totals().failed == 0 ? visitors[static_cast<std::size_t>(r)] : -1;
    }
    s.require("closed to walkers passing through, the street takes none; one-way, none against its way; open, both ways",
              e[2][0] == 0 && e[2][1] == 0 && e[1][1] == 0 && e[1][0] > 0 && e[0][0] > 0 && e[0][1] > 0,
              "open " + num(e[0][0], 5) + " / " + num(e[0][1], 5) + ", one-way " + num(e[1][0], 5) + " / " + num(e[1][1], 5) + ", closed " + num(e[2][0], 5) + " / " + num(e[2][1], 5));
    s.require("the street's visitors reach it and walk back under every regime, none lost", visitors[0] > 0 && visitors[1] > 0 && visitors[2] > 0,
              num(visitors[0], 4) + ", " + num(visitors[1], 4) + ", " + num(visitors[2], 4) + " in the first hour");
    // The regime decision on the walkers: POMCP plans on the ABM and is advanced on the world's observation.
    um::RegimeModel::Setup setup;
    setup.routes = RR, setup.epoch = 900, setup.epochs = 4, setup.sensed = {0}, setup.bins = {50, 150}, setup.scaled = {3};
    setup.scales = {0.5, 1.5}, setup.prior = {0.5, 0.5};
    um::Options light = o;
    light.light = true;
    um::Simulation start(TRp, {RR->values[0][0], RR->values[0][1]}, st, counters, light);
    um::RegimeModel model(start, setup);
    td::Pomcp<um::RegimeModel> planner(model, {.simulations = 60, .reward_span = 50, .depth = 4, .particles = 40, .refill_attempts = 200, .seed = 7});
    um::Simulation world = start;
    world.set_scale(3, 1.5);
    bool ok = true;
    std::string acts;
    for (std::size_t k = 0; k < 2 && ok; k++) {
      const std::size_t a = planner.plan();
      const auto [obs, reward] = model.act(world, a);
      acts += std::string(acts.empty() ? "" : ", ") + um::regime_name(static_cast<um::Regime>(a));
      ok = a < 3 && std::isfinite(reward) && planner.advance(a, obs);
    }
    s.require("POMCP plans the street's regime on the walkers and follows the world's counts", ok, acts);
  }

  // ---- the corridor: exact ---------------------------------------------------------------------------------------------
  {
    const um::Corridor C;
    const td::Pomdp P = um::corridor_pomdp(C);
    td::validate(P);
    const um::CorridorAgents A{C};
    // The agents' law of one step against the table's, every state and action.
    double worst = 0;
    const int n = 20000;
    for (std::size_t a = 0; a < 3; a++)
      for (std::size_t k = 0; k < 2; k++)
        for (std::size_t g = 0; g <= C.capacity; g++) {
          um::CorridorAgents::State x;
          x.level = k;
          for (std::uint32_t q = 0; q < g; q++) x.groups.push_back(q);
          x.next = static_cast<std::uint32_t>(g);
          std::vector<double> freq(P.states * P.observations, 0.0);
          calx::random::Stream rng(99, static_cast<std::uint32_t>(a), static_cast<std::uint32_t>(k), static_cast<std::uint32_t>(g));
          for (int i = 0; i < n; i++) {
            const auto st = A.step(x, a, rng);
            freq[A.index(st.next) * P.observations + st.observation] += 1.0 / n;
          }
          const std::size_t s0 = C.state(k, g);
          for (std::size_t s2 = 0; s2 < P.states; s2++)
            for (std::size_t o = 0; o < P.observations; o++) {
              const double p = P.transition(a, s0, s2) * P.observe(a, s2, o);
              if (p < 1e-4) continue;
              worst = std::max(worst, std::fabs(freq[s2 * P.observations + o] - p) / std::sqrt(p * (1 - p) / n));
            }
        }
    s.check("the corridor's agents have the law of its POMDP (largest standardised deviation over every step's (s′, o))", worst, 5);
    const auto hs = td::hsvi(P, {.epsilon = 1e-3});
    const auto blind = td::blind_lower_bound(P);
    const int episodes = 40, steps = 60;
    double sum = 0, sq = 0;
    for (int e = 0; e < episodes; e++) {
      calx::random::Stream world(20261105, 0x776f726cu, static_cast<std::uint32_t>(e), 0);
      td::Pomcp<um::CorridorAgents> planner(A, {.simulations = 3000, .reward_span = P.reward_span(), .particles = 1000, .seed = 500u + static_cast<std::uint64_t>(e)});
      planner.set_leaf([&](const um::CorridorAgents::State& x) {
        double v = -1e300;
        for (const auto& b : blind.vectors) v = std::max(v, b.values[A.index(x)]);
        return v;
      });
      auto state = A.initial(world);
      double g = 0, disc = 1;
      for (int t = 0; t < steps; t++) {
        const std::size_t act = planner.plan();
        const auto st = A.step(state, act, world);
        g += disc * st.reward, disc *= C.discount;
        state = st.next;
        if (!planner.advance(act, st.observation)) break;
      }
      sum += g, sq += g * g;
    }
    const double mean = sum / episodes, se = std::sqrt((sq - sum * sum / episodes) / (episodes - 1) / episodes);
    const double truncation = std::pow(C.discount, steps) * P.reward_span() / (1 - C.discount);
    s.require("POMCP on the corridor's agents earns the exact value (HSVI's interval) within four standard errors and the truncation",
              mean >= hs.lower_value - 4 * se - truncation && mean <= hs.upper_value + 4 * se,
              num(mean, 5) + " ± " + num(se, 2) + " against [" + num(hs.lower_value, 6) + ", " + num(hs.upper_value, 6) + "]");
  }
  return s.done();
}
