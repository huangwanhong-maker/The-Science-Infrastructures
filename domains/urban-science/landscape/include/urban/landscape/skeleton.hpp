// The event skeleton of a landscape (theory/landscape-dynamics.md §6): the discrete part of a district's state —
// two stores open or shut, a rule in force (a permit for a lane), the rule carried out (the lane built), a meaning
// adopted (the lane recognised) — composed from five open systems at their interfaces, and what can be said of it
// exactly: its stationary law and its irreversibility, which coarser descriptions forget something exactly and
// which do not (with the error of the best coarse chain when they do not, and the memory they keep), how
// interventions interact on an outcome, and what a record of it says about the stores' coupling.
//
// The parts (physis composition.hpp). A and B, the stores (shut 0, open 1); P the permit (0, 1); W the lane (0, 1);
// R its recognition (0, 1). The interfaces are relations between the parts' states: the lane can stand only under a
// permit (w ≤ p) and be recognised only while it stands (r ≤ w), so a permit cannot be withdrawn while its lane
// stands nor a recognised lane removed — those local moves are refused at the interface, not forbidden by hand.
// The mechanisms, each declared once on the parts it reads:
//   a store opens        on {A, B, W}: α(1 + η·[the other open] + ξ·[the lane there])·π (π a promotion, 1 for none)
//   a store shuts        on itself, at β
//   the permit           granted at κ_g, withdrawn at κ_r
//   the lane             built on {P, W} under a permit at κ_b, removed at κ_d
//   its recognition      gained on {A, B, W, R} while it stands at κ_u·u·(1 + ν·(a + b)), lapsing at κ_l: u the
//                        lane's use with both stores shut (its walkers against a reference, the continuous part of
//                        the landscape held at a reference value), ν what each open store adds to it
// The two orders of composition — ((A ⊙ B) ⊙ P) ⊙ W) ⊙ R and A ⊙ (B ⊙ (P ⊙ (W ⊙ R))) — give the same rates for
// every mechanism between every two joint states (each joint mechanism added at the first composition holding all
// its parts): 2 · 2 · 4 = 16 joint states.
//
// What follows from the rates. The institutions' chain (P, W, R) is a path of four states, a tree, so always in
// detailed balance; the stores answer the lane (ξ) and its recognition answers the stores (ν), so the whole is in
// detailed balance exactly when ξ_A = ξ_B = ν = 0 and η_A = η_B, its stationary law then the product of the
// stores' Gibbs law e^{−S}/Z and the institutions' own. A meaning adopted through use is a driven process. A coarser description that forgets a part lumps exactly when no remaining
// mechanism reads what it forgets: the stores can forget the institutions when ξ = 0; the permit cannot be
// forgotten (building the lane reads it), nor the recognition (removing the lane is blocked by it) — an adopted
// meaning changes what may happen next. When a description does not lump, its aggregation by the stationary law
// (physis lumping.hpp) keeps the stationary law of what it retains, and over a horizon T its whole path law is
// within 1 − e^{−δT} of the fine one's projection, δ the outward discrepancy: within 10% for T ≤ −ln 0.9/δ. The anchored Möbius interactions of
// interventions (morphe mobius.hpp) say what promoting both stores, or granting permits faster, does beyond each
// alone, on an outcome chosen. Checks: tests/check_landscape.cpp.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <calx/linalg/dense.hpp>
#include <calx/random/stream.hpp>
#include <morphe/combinatorics/mobius.hpp>
#include <physis/statistical_mechanics/composition.hpp>
#include <physis/statistical_mechanics/jump_process.hpp>
#include <physis/statistical_mechanics/lumping.hpp>
#include <physis/statistical_mechanics/master_equation.hpp>

namespace urban::landscape {

namespace sm = physis::statistical_mechanics;
using calx::linalg::Matrix;

struct SkeletonParams {
  double alpha_a = 1.6, alpha_b = 1.6, beta_a = 0.7, beta_b = 0.7;  // per hour
  double eta_a = 0.5, eta_b = 0.5, xi_a = 0.8, xi_b = 0.8;
  double grant = 0.3, revoke = 0.05, build = 0.5, remove = 0.06, recognise = 0.4, lapse = 0.02, use = 1.0, draw = 0.5;
  double promote_a = 1, promote_b = 1, promote_grant = 1;
};

inline const std::vector<std::string>& skeleton_parts() {
  static const std::vector<std::string> p{"A", "B", "P", "W", "R"};
  return p;
}

// The skeleton composed in one of two orders, "left" or "right".
inline sm::OpenSystem skeleton(const SkeletonParams& q, const std::string& order = "left") {
  const auto A = sm::atomic("A", 2, {}, {{"shuts", 1, 0, q.beta_a}});
  const auto B = sm::atomic("B", 2, {}, {{"shuts", 1, 0, q.beta_b}});
  const auto P = sm::atomic("P", 2, {{"p→w", {0, 1}}}, {{"granted", 0, 1, q.grant * q.promote_grant}, {"withdrawn", 1, 0, q.revoke}});
  const auto W = sm::atomic("W", 2, {{"w←p", {0, 1}}, {"w→r", {0, 1}}}, {{"removed", 1, 0, q.remove}});
  const auto R = sm::atomic("R", 2, {{"r←w", {0, 1}}}, {{"lapses", 1, 0, q.lapse}});
  const auto under = [](int above, int below) { return below <= above; };  // the second part's state at most the first's
  const sm::JointMechanism a_opens{"A opens", {"A", "B", "W"}, [q](const sm::ComponentState& z) {
    std::vector<std::pair<sm::ComponentState, double>> out;
    if (z.at("A") == 0) out.push_back({{{"A", 1}}, q.alpha_a * (1 + q.eta_a * z.at("B") + q.xi_a * z.at("W")) * q.promote_a});
    return out;
  }};
  const sm::JointMechanism b_opens{"B opens", {"A", "B", "W"}, [q](const sm::ComponentState& z) {
    std::vector<std::pair<sm::ComponentState, double>> out;
    if (z.at("B") == 0) out.push_back({{{"B", 1}}, q.alpha_b * (1 + q.eta_b * z.at("A") + q.xi_b * z.at("W")) * q.promote_b});
    return out;
  }};
  const sm::JointMechanism built{"W built", {"P", "W"}, [q](const sm::ComponentState& z) {
    std::vector<std::pair<sm::ComponentState, double>> out;
    if (z.at("P") == 1 && z.at("W") == 0) out.push_back({{{"W", 1}}, q.build});
    return out;
  }};
  const sm::JointMechanism gained{"R gained", {"A", "B", "W", "R"}, [q](const sm::ComponentState& z) {
    std::vector<std::pair<sm::ComponentState, double>> out;
    if (z.at("W") == 1 && z.at("R") == 0) out.push_back({{{"R", 1}}, q.recognise * q.use * (1 + q.draw * (z.at("A") + z.at("B")))});
    return out;
  }};
  const std::vector<sm::JointMechanism> J{a_opens, b_opens, built, gained};
  const auto any = [](int, int) { return true; };
  if (order == "left") {
    auto x = sm::compose(A, B, "", "", any, J).system;
    x = sm::compose(x, P, "", "", any, J).system;
    x = sm::compose(x, W, "p→w", "w←p", under, J).system;
    return sm::compose(x, R, "w→r", "r←w", under, J).system;
  }
  const auto wr = sm::compose(W, R, "w→r", "r←w", under, J).system;
  const auto pwr = sm::compose(P, wr, "p→w", "w←p", under, J).system;
  const auto bpwr = sm::compose(B, pwr, "", "", any, J).system;
  return sm::compose(A, bpwr, "", "", any, J).system;
}

// A coarser description: the parts it keeps, and for the stores whether it keeps which is open or only how many.
struct View {
  std::string name;
  std::set<std::string> keep;
  bool count_stores = false;
};
inline std::vector<View> skeleton_views() {
  return {{"the stores alone", {"A", "B"}, false},
          {"without the lane", {"A", "B", "P", "R"}, false},
          {"without the permit", {"A", "B", "W", "R"}, false},
          {"without the recognition", {"A", "B", "P", "W"}, false},
          {"how many stores are open, with the institutions", {"P", "W", "R"}, true}};
}
// The block of each joint state under a view.
inline std::vector<std::size_t> view_blocks(const sm::OpenSystem& S, const View& v) {
  std::map<std::vector<int>, std::size_t> index;
  std::vector<std::size_t> block;
  for (const auto& z : S.states) {
    std::vector<int> key;
    if (v.count_stores) key.push_back(static_cast<int>(z.at("A") + z.at("B")));
    for (const auto& part : skeleton_parts())
      if (v.keep.count(part)) key.push_back(static_cast<int>(z.at(part)));
    block.push_back(index.emplace(key, index.size()).first->second);
  }
  return block;
}

struct ViewVerdict {
  std::string name;
  std::size_t blocks = 0;
  bool eventwise = false, whole = false;
  double discrepancy = 0;   // δ of the stationary aggregation, per hour
  double horizon = 0;       // hours over which its path law stays within 10% of the fine one's (−ln 0.9/δ; ∞ if exact)
  // Where it fails worst, when it does: a mechanism that leaves two states the description cannot tell apart, `state`
  // and `other`, for one coarse state at rates differing by `by` per hour — what it reads, the description forgot.
  std::string mechanism, reason;
  std::size_t state = 0, other = 0;
  double by = 0;
};

struct SkeletonAnalysis {
  sm::OpenSystem system;
  std::vector<double> pi;
  bool reversible_rates = false;   // every rate has its reverse
  double entropy_production = std::numeric_limits<double>::quiet_NaN();
  bool detailed_balance = false;
  std::vector<sm::Cycle> cycles;
  bool coherent = false;           // the two orders of composition agree
  std::vector<ViewVerdict> views;
};

inline SkeletonAnalysis analyse_skeleton(const SkeletonParams& q) {
  SkeletonAnalysis out;
  out.system = skeleton(q, "left");
  out.coherent = out.system.rates_by_components() == skeleton(q, "right").rates_by_components();
  const auto k = out.system.rates();
  out.pi = sm::stationary(k);
  out.reversible_rates = true;
  for (std::size_t i = 0; i < k.n; i++)
    for (std::size_t j = 0; j < k.n; j++)
      if (i != j && (k(i, j) > 0) != (k(j, i) > 0)) out.reversible_rates = false;
  if (out.reversible_rates) {
    const auto irr = sm::irreversibility(k);
    out.entropy_production = irr.entropy_production, out.detailed_balance = irr.detailed_balance, out.cycles = irr.cycles;
  }
  std::set<std::string> mechanisms;
  for (const auto& m : out.system.moves) mechanisms.insert(m.mechanism);
  for (const auto& v : skeleton_views()) {
    ViewVerdict r;
    r.name = v.name;
    const auto block = view_blocks(out.system, v);
    r.blocks = sm::block_count(block, k.n);
    std::vector<Matrix> mech;
    double worst = 0;
    for (const auto& name : mechanisms) {
      const auto km = out.system.rates(name);
      const auto L = sm::strongly_lumpable(km, block);
      mech.push_back(km);
      if (!L.lumpable && L.discrepancy > worst) {
        worst = L.discrepancy;
        r.mechanism = name, r.state = L.from, r.other = L.other, r.by = L.discrepancy;
        const auto text = [&](std::size_t i) {
          std::string t;
          for (const auto& part : skeleton_parts()) t += (t.empty() ? "" : " ") + part + "=" + std::to_string(out.system.states[i].at(part));
          return t;
        };
        r.reason = "\"" + name + "\" from " + text(L.from) + " and from " + text(L.other) + ", one state to this description, differ by " +
                   std::to_string(L.discrepancy) + " per hour";
      }
    }
    r.eventwise = sm::eventwise_lumpable(mech, block).lumpable;
    r.whole = sm::strongly_lumpable(k, block).lumpable;
    const auto co = sm::aggregated(k, block, out.pi);
    r.discrepancy = sm::outward_discrepancy(k, block, co);
    double top = 0;
    for (double x : k.a) top = std::max(top, x);
    if (r.discrepancy <= 1e-12 * top) r.discrepancy = 0;  // rounding of an exact aggregation
    r.horizon = r.discrepancy > 0 ? -std::log(0.9) / r.discrepancy : std::numeric_limits<double>::infinity();
    out.views.push_back(std::move(r));
  }
  return out;
}

// An outcome at stationarity: "both" P(both open), "open" E[number open], "recognised" P(r = 1), "lane" P(w = 1).
inline double skeleton_outcome(const SkeletonParams& q, const std::string& outcome) {
  const auto S = skeleton(q);
  const auto pi = sm::stationary(S.rates());
  double y = 0;
  for (std::size_t i = 0; i < S.states.size(); i++) {
    const auto& z = S.states[i];
    const double a = static_cast<double>(z.at("A")), b = static_cast<double>(z.at("B"));
    y += pi[i] * (outcome == "both" ? a * b : outcome == "open" ? a + b : outcome == "recognised" ? z.at("R") : z.at("W"));
  }
  return y;
}

// The anchored interactions of three interventions — promoting A, promoting B (each opening rate × (1 + gain)),
// and granting permits faster (× (1 + gain)) — on an outcome: indexed by subset mask (bit 0 A, bit 1 B, bit 2 the
// permits), g[0] the outcome with none, g[7] the three-way interaction.
inline std::vector<double> intervention_interactions(const SkeletonParams& q, double gain, const std::string& outcome) {
  return morphe::combinatorics::anchored_interactions<double>(3, [&](std::uint32_t S) {
    SkeletonParams r = q;
    r.promote_a = (S & 1) ? 1 + gain : 1;
    r.promote_b = (S & 2) ? 1 + gain : 1;
    r.promote_grant = (S & 4) ? 1 + gain : 1;
    return skeleton_outcome(r, outcome);
  });
}

// The time a coarse state lasts under a view, entered with the stationary law on its fine states: phase-type.
inline sm::PhaseType skeleton_dwell(const SkeletonAnalysis& a, const View& v, std::size_t which_block) {
  const auto block = view_blocks(a.system, v);
  std::vector<double> alpha;
  double z = 0;
  for (std::size_t i = 0; i < block.size(); i++)
    if (block[i] == which_block) alpha.push_back(a.pi[i]), z += a.pi[i];
  for (auto& x : alpha) x /= z;
  return sm::dwell_law(a.system.rates(), block, which_block, alpha);
}

// A record of the skeleton over `hours` (every change of every part observed), simulated from a start.
struct SkeletonRecord {
  std::size_t start = 0;
  std::vector<std::pair<double, std::size_t>> changes;
  double end = 0;
};
inline SkeletonRecord simulate_skeleton(const SkeletonParams& q, double hours, std::uint64_t seed) {
  const auto S = skeleton(q);
  const auto k = S.rates();
  calx::random::Stream rng(seed, 5, 7, 11);
  const auto occ = [&](const std::size_t& z) {
    std::vector<sm::Occurrence<std::size_t>> o;
    for (std::size_t j = 0; j < k.n; j++)
      if (j != z && k(z, j) > 0) o.push_back({"change", "", k(z, j), [j](const std::size_t&, calx::random::Stream&) { return j; }});
    return o;
  };
  const auto R = sm::simulate<std::size_t>(0, occ, hours, rng);
  SkeletonRecord rec;
  rec.end = hours;
  for (const auto& e : R.events) rec.changes.push_back({e.time, e.successor});
  return rec;
}
// The log-density of a record under parameters q.
inline double skeleton_log_density(const SkeletonParams& q, const SkeletonRecord& R) {
  const auto k = skeleton(q).rates();
  std::vector<std::size_t> id(k.n);
  for (std::size_t s = 0; s < k.n; s++) id[s] = s;
  std::vector<double> p0(k.n, 0.0);
  p0[R.start] = 1;
  sm::CoarsePath path;
  double t = 0;
  std::size_t here = R.start;
  for (const auto& [time, next] : R.changes) path.blocks.push_back(here), path.dwell.push_back(time - t), t = time, here = next;
  path.blocks.push_back(here), path.dwell.push_back(R.end - t);
  return sm::coarse_log_density(k, id, p0, path);
}

// The posterior over a common coupling η (η_A = η_B) of the worlds given, from a record and a prior.
struct Posterior {
  std::vector<double> worlds, probability;
  template <class Phi>
  double of(Phi phi) const {
    double s = 0;
    for (std::size_t k = 0; k < worlds.size(); k++)
      if (phi(worlds[k])) s += probability[k];
    return s;
  }
};
inline Posterior coupling_posterior(SkeletonParams q, const SkeletonRecord& R, const std::vector<double>& etas, const std::vector<double>& prior) {
  Posterior out;
  out.worlds = etas;
  std::vector<double> lp(etas.size());
  double top = -std::numeric_limits<double>::infinity();
  for (std::size_t k = 0; k < etas.size(); k++) {
    q.eta_a = q.eta_b = etas[k];
    lp[k] = std::log(prior[k]) + skeleton_log_density(q, R);
    top = std::max(top, lp[k]);
  }
  double z = 0;
  for (double& x : lp) z += x = std::exp(x - top);
  for (double x : lp) out.probability.push_back(x / z);
  return out;
}

}  // namespace urban::landscape
