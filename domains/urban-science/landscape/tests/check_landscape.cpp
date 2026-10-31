// Checks of the landscape module (district.hpp, skeleton.hpp, process.hpp, value.hpp). The building that was here is
// the architecture domain's (space, environment).
//
// The district: a 4 × 3 grid of junctions 100 m by 80 m apart (one street curved), with a dead end in one block.
// - Its complex: 13 junctions, 18 segments, 5 blocks, the block with the dead end left a hole (b₁ = 1), Euler
//   characteristic 0, regular; a network whose segments cross without a junction refused.
// - A lane between two stores facing one block: three rewrites, two junctions, three segments and a block more,
//   the homology kept, the two new blocks' areas and heat capacities summing to the old, the temperature copied;
//   the subdivisions concurrent and the split after both (two records in the trace); closed again, the district
//   back up to isomorphism, its capacities, temperatures and walkers kept, the stores on the merged streets.
//   Two stores facing no common block, or one segment, refused.
// The landscape's event skeleton (two stores, a permit, a lane, its recognition, composed at their interfaces):
// - 16 admissible joint states, the two orders of composition agreeing on every rate;
// - driven while the stores answer the lane; in detailed balance when they do not and answer each other equally, its
//   law then the stores' Gibbs law times the institutions' own;
// - the stores alone forget the institutions exactly when they do not answer the lane; the permit and the
//   recognition cannot be forgotten; an exact description's best coarse chain has no discrepancy and stays within
//   10% for ever, an inexact one for a finite horizon, and remembers (both shut lasts a non-exponential time);
// - interventions: stores that answer nothing show no interaction on the number open, one between the promotions
//   on both open, and where they answer the lane faster permits interact with the promotions;
// - a long record with η = 1 puts most of the posterior there, coherently over overlapping propositions.
// The landscape running for a day: walkers kept (Σ n changes by what entered less what left, the transfers at
// rewrites included); every stratum regular with b₁ the holes; the stores closed at 22:00 by a forced jump; the
// lane after the permit in the causal order; frames every quarter hour.
// Value: transport, conversion, copying and a ledger's claims, and the lattice of what they conserve.
// Walkers as a reaction network: the Doi–Peliti operator is the master equation, and its first moment the
// mean-field flow.

#include <cmath>
#include <numbers>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <physis/chemistry/doi_peliti.hpp>
#include <physis/statistical_mechanics/jump_process.hpp>
#include <urban/landscape/district.hpp>
#include <urban/landscape/process.hpp>
#include <urban/landscape/skeleton.hpp>
#include <urban/landscape/value.hpp>

using calx::verification::Suite;
namespace ul = urban::landscape;
namespace mt = morphe::topology;
namespace tr = morphe::formal::traces;
namespace sm = physis::statistical_mechanics;

namespace {

std::string text(const std::vector<std::size_t>& v) {
  std::string s;
  for (std::size_t x : v) s += (s.empty() ? "" : " ") + std::to_string(x);
  return s;
}

ul::Network grid_network() {
  ul::Network N;
  for (int y = 0; y < 3; y++)
    for (int x = 0; x < 4; x++) N.junctions.push_back({100.0 * x, 80.0 * y});
  const auto at = [](int x, int y) { return static_cast<std::size_t>(y * 4 + x); };
  for (int y = 0; y < 3; y++)
    for (int x = 0; x < 3; x++) {
      ul::Segment s{at(x, y), at(x + 1, y), {}, "street", "Street " + std::to_string(y)};
      if (y == 0 && x == 0) s.line = {{0, 0}, {50, -6}, {100, 0}};  // a curved street
      N.segments.push_back(s);
    }
  for (int y = 0; y < 2; y++)
    for (int x = 0; x < 4; x++) N.segments.push_back({at(x, y), at(x, y + 1), {}, "street", "Avenue " + std::to_string(x)});
  N.junctions.push_back({250, 120});
  N.segments.push_back({at(3, 2), N.junctions.size() - 1, {}, "footway", "a dead end"});
  return N;
}

double total(const ul::District& D, const std::string& type, const std::string& f) {
  double s = 0;
  for (const auto& c : D.complex.cells)
    if (c.type == type || (type == "junction" && c.type == "portal")) s += D.value(c.id, f);
  return s;
}

}  // namespace

int main() {
  Suite s;

  // ---- the district ----
  auto D = ul::district(grid_network());
  s.require("the district: 13 junctions, 18 segments, 5 blocks", text(D.complex.counts()) == "13 18 5", text(D.complex.counts()));
  s.require("the block with a dead end is a hole: b₁ = 1, Euler characteristic 0", D.holes == 1 && text(morphe::algebra::betti(mt::chain_complex(D.complex))) == "1 1 0" && D.complex.euler_characteristic() == 0);
  {
    ul::Network X;
    X.junctions = {{0, 0}, {10, 10}, {0, 10}, {10, 0}};
    X.segments = {{0, 1, {}, "street", ""}, {2, 3, {}, "street", ""}};
    bool refused = false;
    try { ul::district(X); } catch (const std::invalid_argument&) { refused = true; }
    s.require("segments that cross without a junction are refused", refused);
    std::vector<std::string> dropped;
    const auto Y = ul::district(ul::planarised(X, &dropped));
    s.require("or, planarised, the shorter of two crossing segments dropped and named", dropped.size() == 1 && Y.complex.counts().size() == 2 && Y.complex.counts()[1] == 1);
  }
  const std::size_t A = ul::attach_store(D, {"the bakery", "food", {150, -5}});
  const std::size_t B = ul::attach_store(D, {"the bookshop", "retail", {140, 85}});
  s.require("each store faces its nearest segment", std::fabs(D.stores[A].s - 50) < 1e-9 && std::fabs(D.stores[B].s - 40) < 1e-9 && D.name.at(D.stores[A].segment) == "Street 0");

  // ---- the lane ----
  {
    for (const auto id : D.of(0)) D.set(id, "walkers", 10.0 + id % 7);
    const auto plan = ul::plan_lane(D, A, B);
    s.require("the lane runs through the block both stores face", plan.refused.empty() && plan.block != 0, plan.refused);
    ul::District L = D;
    const double C0 = total(L, "block", "capacity"), W0 = total(L, "junction", "walkers"), area0 = L.area(plan.block), T0 = L.value(plan.block, "temperature");
    const auto steps = ul::build_lane(L, plan);
    s.require("built by three rewrites: two junctions, three segments and a block more", steps.size() == 3 && text(L.complex.counts()) == "15 21 6", text(L.complex.counts()));
    s.require("regular, its homology kept", mt::admissible(L.complex).regular && text(morphe::algebra::betti(mt::chain_complex(L.complex))) == "1 1 0");
    const std::uint64_t lane = steps.back().created[0], f1 = steps.back().created[1], f2 = steps.back().created[2];
    s.check("the two new blocks' areas sum to the old", L.area(f1) + L.area(f2) - area0, 1e-9, "m²");
    s.check("their capacities too", total(L, "block", "capacity") - C0, 1e-6, "J/K");
    s.require("the temperature copied to both", L.value(f1, "temperature") == T0 && L.value(f2, "temperature") == T0);
    s.check("the walkers kept (a new junction has none)", total(L, "junction", "walkers") - W0, 1e-12);
    const auto ends = L.complex.ends(L.at(lane));
    s.check("the lane's ends are the two points facing the stores", std::hypot(L.line.at(lane).front()[0] - plan.a[0], L.line.at(lane).front()[1] - plan.a[1]) + std::hypot(L.line.at(lane).back()[0] - plan.b[0], L.line.at(lane).back()[1] - plan.b[1]), 1e-9, "m");
    s.require("whose junctions are portals", L.complex.cells[ends[0]].type == "portal" && L.complex.cells[ends[1]].type == "portal");
    s.require("the two subdivisions refine the district (the same streets, cut finer), the split intervenes", steps[0].refinement && steps[1].refinement && !steps[2].refinement);
    const auto P = tr::causal_order({steps[0].footprint, steps[1].footprint, steps[2].footprint});
    s.require("the two subdivisions are concurrent, the split after both: two records in the trace", P.concurrent(0, 1) && P.before[0][2] && P.before[1][2] && tr::linearisations(P).str() == "2");
    s.require("the stores face the new pieces of their streets", L.stores[A].segment != D.stores[A].segment && L.line.count(L.stores[A].segment));
    auto back = L;
    const auto undo = ul::close_lane(back, lane);
    s.require("closed by three rewrites, the district back up to isomorphism", undo.size() == 3 && mt::isomorphic(back.complex, D.complex), std::to_string(undo.size()));
    s.check("its capacities kept", total(back, "block", "capacity") - C0, 1e-6, "J/K");
    {
      auto stale = back;
      bool refused = false;
      try { ul::build_lane(stale, plan); } catch (const std::invalid_argument&) { refused = true; }
      auto again = back;
      const auto replan = ul::plan_lane(again, 0, 1);
      s.require("closed, the old plan is refused (its streets are merged anew), and a new plan builds the same lane again",
                refused && replan.refused.empty() && ul::build_lane(again, replan).size() == 3 && mt::isomorphic(again.complex, L.complex), replan.refused);
      s.check("from the same two points", std::hypot(replan.a[0] - plan.a[0], replan.a[1] - plan.a[1]) + std::hypot(replan.b[0] - plan.b[0], replan.b[1] - plan.b[1]), 1e-9, "m");
    }
    s.check("and its walkers", total(back, "junction", "walkers") - W0, 1e-12);
    s.require("the stores face the merged streets", back.line.count(back.stores[A].segment) && std::fabs(back.stores[A].s - 50) < 1e-6, std::to_string(back.stores[A].s));
    ul::District F = D;
    const std::size_t far = ul::attach_store(F, {"far", "retail", {290, 165}});
    s.require("two stores that face no common block are refused", !ul::plan_lane(F, A, far).refused.empty());
    const std::size_t near = ul::attach_store(F, {"next door", "retail", {160, -5}});
    s.require("and two on one segment need no lane", !ul::plan_lane(F, A, near).refused.empty());
  }

  // ---- the landscape's event skeleton ----
  {
    ul::SkeletonParams q;
    const auto a = ul::analyse_skeleton(q);
    s.require("the skeleton: 16 joint states of the stores, the permit, the lane and its recognition", a.system.states.size() == 16);
    s.require("composed in either order, the same rate for every mechanism between every two states", a.coherent);
    bool admissible = true;
    for (const auto& z : a.system.states) admissible = admissible && z.at("W") <= z.at("P") && z.at("R") <= z.at("W");
    s.require("every joint state admissible: a lane only under a permit, recognised only while it stands", admissible);
    s.require("the stores answer the lane, it does not answer them: driven (entropy production > 0)", a.reversible_rates && !a.detailed_balance && a.entropy_production > 1e-6, std::to_string(a.entropy_production));
    ul::SkeletonParams r = q;
    r.xi_a = r.xi_b = 0, r.eta_b = r.eta_a;
    {
      const auto used = ul::analyse_skeleton(r);
      s.require("a meaning adopted through use is driven: with ξ = 0 and η_A = η_B, recognition drawn by the open stores breaks detailed balance",
                !used.detailed_balance && used.entropy_production > 1e-9, std::to_string(used.entropy_production));
    }
    r.draw = 0;
    const auto b = ul::analyse_skeleton(r);
    s.require("with ξ = ν = 0 and η_A = η_B: detailed balance", b.detailed_balance);
    // Its law: the stores' Gibbs law times the institutions' path chain.
    const double ka = r.alpha_a / r.beta_a, kb = r.alpha_b / r.beta_b;
    std::map<std::vector<int>, double> inst;
    {
      const double w1 = r.grant / r.revoke, w2 = w1 * r.build / r.remove, w3 = w2 * r.recognise * r.use / r.lapse;
      const double z = 1 + w1 + w2 + w3;
      inst[{0, 0, 0}] = 1 / z, inst[{1, 0, 0}] = w1 / z, inst[{1, 1, 0}] = w2 / z, inst[{1, 1, 1}] = w3 / z;
    }
    double zs = 0, worst = 0;
    for (int x = 0; x < 2; x++)
      for (int y = 0; y < 2; y++) zs += std::pow(ka, x) * std::pow(kb, y) * std::pow(1 + r.eta_a, x * y);
    for (std::size_t i = 0; i < b.system.states.size(); i++) {
      const auto& z = b.system.states[i];
      const int x = static_cast<int>(z.at("A")), y = static_cast<int>(z.at("B"));
      const double want = std::pow(ka, x) * std::pow(kb, y) * std::pow(1 + r.eta_a, x * y) / zs * inst.at({static_cast<int>(z.at("P")), static_cast<int>(z.at("W")), static_cast<int>(z.at("R"))});
      worst = std::max(worst, std::fabs(b.pi[i] - want));
    }
    s.check("and its law the stores' Gibbs law times the institutions' own", worst, 1e-12);
    const auto view = [](const ul::SkeletonAnalysis& x, const std::string& n) {
      for (const auto& v : x.views)
        if (v.name == n) return v;
      return ul::ViewVerdict{};
    };
    {
      const auto& v = view(a, "the stores alone");
      const auto &x = a.system.states[v.state], &y = a.system.states[v.other];
      s.require("where the stores alone fail, the two states differ only in what they forget, and the mechanism reads it",
                x.at("A") == y.at("A") && x.at("B") == y.at("B") && x.at("W") != y.at("W") && v.mechanism.find("opens") != std::string::npos, v.reason);
    }
    s.require("the stores alone forget the institutions exactly when the stores do not answer the lane",
              view(b, "the stores alone").eventwise && !view(a, "the stores alone").eventwise && !view(a, "the stores alone").whole);
    s.require("the permit cannot be forgotten: building the lane reads it", !view(b, "without the permit").whole);
    s.require("nor the recognition: removing the lane is refused while it holds", !view(b, "without the recognition").whole);
    s.check("an exact description's best coarse chain has no discrepancy", view(b, "the stores alone").discrepancy, 1e-12);
    s.require("an exact one stays within 10% for ever, an inexact one for a finite horizon", std::isinf(view(b, "the stores alone").horizon) && view(a, "the stores alone").discrepancy > 0 && std::isfinite(view(a, "the stores alone").horizon) && view(a, "the stores alone").horizon > 0);
    const auto alone = ul::skeleton_views()[0];
    const auto cv2 = [](const sm::PhaseType& P) { const double m1 = P.moment(1), m2 = P.moment(2); return m2 / (m1 * m1) - 1; };
    const auto block_shut = [&](const ul::SkeletonAnalysis& x) {
      const auto blk = ul::view_blocks(x.system, alone);
      for (std::size_t i = 0; i < blk.size(); i++)
        if (x.system.states[i].at("A") == 0 && x.system.states[i].at("B") == 0) return blk[i];
      return std::size_t{0};
    };
    s.check("both stores shut lasts an exponential time when the description is exact", cv2(ul::skeleton_dwell(b, alone, block_shut(b))) - 1, 1e-10);
    s.require("and not when it is not: the stores alone remember the lane", std::fabs(cv2(ul::skeleton_dwell(a, alone, block_shut(a))) - 1) > 1e-4);
    ul::SkeletonParams ind = q;
    ind.eta_a = ind.eta_b = ind.xi_a = ind.xi_b = 0;
    const auto open = ul::intervention_interactions(ind, 0.5, "open"), both = ul::intervention_interactions(ind, 0.5, "both");
    double none = 0;
    for (std::uint32_t S = 0; S < 8; S++)
      if (__builtin_popcount(S) >= 2) none = std::max(none, std::fabs(open[S]));
    s.check("stores that answer nothing: the number open has no interaction of any two or three interventions", none, 1e-12);
    s.require("the chance both are open has one between the two promotions", std::fabs(both[3]) > 1e-4);
    const auto three = ul::intervention_interactions(q, 0.5, "both");
    s.require("and where the stores answer the lane, permits granted faster interact with the promotions", std::fabs(three[5]) > 1e-6 && std::fabs(three[7]) > 1e-9, std::to_string(three[7]));
    ul::SkeletonParams truth = q;
    truth.eta_a = truth.eta_b = 1;
    const auto rec = ul::simulate_skeleton(truth, 3000, 13);
    const auto post = ul::coupling_posterior(q, rec, {0.0, 0.5, 1.0, 2.0}, {0.25, 0.25, 0.25, 0.25});
    s.require("a long record with η = 1 puts most of the posterior there", post.probability[2] > 0.8, std::to_string(post.probability[2]));
    s.require("P(η > 0) ≥ P(η ≥ 1): the second's worlds are among the first's", post.of([](double e) { return e > 0; }) >= post.of([](double e) { return e >= 1; }));
  }

  // ---- the landscape running ----
  {
    ul::LandscapeParams P;
    P.seed = 3;
    P.entries.assign(24, 0.0);
    P.altitude.assign(24, 0.0);
    for (int h = 0; h < 24; h++) {
      P.entries[h] = 200 + 1400 * std::exp(-0.5 * std::pow((h - 13) / 4.0, 2));
      P.altitude[h] = std::max(-0.3, std::sin(std::numbers::pi * (h - 6) / 14.0)) * 0.9;
    }
    const auto run = ul::run_landscape(D, P);
    s.check("walkers kept: Σ n changes by what entered less what left (relative)", (run.walkers1 - run.walkers0 - (run.entered - run.left)) / std::max(1.0, run.entered), 1e-9);
    bool regular = true;
    for (const auto& S : run.strata) regular = regular && mt::admissible(S.district.complex).regular && S.betti.size() >= 2 && S.betti[1] == S.district.holes;
    s.require("every stratum regular, its b₁ the holes", regular && !run.strata.empty());
    bool closed = false, quiet = true, lane = false;
    std::size_t permit = run.events.size(), built = run.events.size();
    for (std::size_t k = 0; k < run.events.size(); k++) {
      const auto& e = run.events[k];
      if (e.forced) closed = closed || std::fabs(e.time - 22) < 1e-6;
      if (e.mechanism == "a store opens" && e.time >= 22) quiet = false;
      if (e.mechanism == "the permit is granted") permit = k;
      if (e.mechanism == "the lane is built") built = k, lane = true;
    }
    s.require("the stores close at 22:00 by a forced jump, and open no more that night", closed && quiet);
    s.require("this day the lane is built", lane, std::to_string(run.events.size()) + " events");
    std::vector<tr::Footprint> fps;
    for (const auto& e : run.events) fps.push_back(e.footprint);
    const auto order = tr::causal_order(fps);
    s.require("the lane after the permit in the causal order", lane && order.before[permit][built]);
    bool after = true;
    for (std::size_t k = built + 1; k < run.events.size(); k++)
      if (run.events[k].mechanism == "a store opens") after = after && order.before[built][k];
    s.require("and every opening after it reads it", after);
    bool frames = run.frames.size() >= 96, positive = true;
    for (const auto& F : run.frames)
      for (const auto& [id, n] : F.walkers) positive = positive && n > -1e-6;
    s.require("a frame every quarter hour, walkers never negative", frames && positive, std::to_string(run.frames.size()));
    double each = 0;
    for (const auto& F : run.frames)
      for (const auto& [id, f] : F.flow) each = std::max(each, std::fabs(F.forth.at(id) - F.back.at(id) - f) / std::max(1.0, F.forth.at(id) + F.back.at(id)));
    s.check("each way's hops differ by the net flow, segment by segment (relative)", each, 1e-12);
    bool registry = true;
    for (const auto& F : run.frames)
      for (const auto& k : F.enabled) {
        if (k == "the lane is built") registry = registry && F.permit && F.lane == 0;
        if (k == "the lane is removed" || k == "the lane is recognised") registry = registry && F.lane != 0 && !F.recognised;
      }
    s.require("each frame's registry of events: the lane buildable only under a permit, removable only while unrecognised", registry);
  }

  // ---- value ----
  {
    // Places: store A, the street, store B. Types: goods, money, knowledge, claim, debt.
    ul::Resources R;
    R.places = 3, R.types = 5;
    R.type_names = {"goods", "money", "knowledge", "claim", "debt"};
    R.links = {{0, 1}, {1, 2}};
    R.transported = {true, true, false, false, false};
    const auto idx = [&](std::size_t p, std::size_t t) { return R.index(p, t); };
    std::vector<long long> make(15, 0), copy(15, 0), sale(15, 0);
    make[idx(0, 0)] = -1, make[idx(0, 1)] = 2;     // one good made into two of money's units (a conversion)
    copy[idx(2, 2)] = 1;                            // A's knowledge copied to B: B gains, A keeps
    sale[idx(0, 3)] = 1, sale[idx(2, 4)] = 1;       // a claim at A, a debt at B
    R.channels = {make, copy, sale};
    const auto W = R.conserved();
    const auto M = R.moves();
    bool kernel = true, no_knowledge = true;
    for (const auto& w : W) {
      for (std::size_t j = 0; j < M.cols; j++) {
        long long dot = 0;
        for (std::size_t i = 0; i < M.rows; i++) dot += w[i] * M(i, j).to_int64();
        kernel = kernel && dot == 0;
      }
      for (std::size_t p = 0; p < 3; p++) no_knowledge = no_knowledge && (p == 2 ? w[idx(p, 2)] == 0 : true);
    }
    s.require("every conserved weighting is orthogonal to every move", kernel && !W.empty());
    const auto in_lattice = [&](const std::vector<long long>& w) {
      for (std::size_t j = 0; j < M.cols; j++) {
        long long dot = 0;
        for (std::size_t i = 0; i < M.rows; i++) dot += w[i] * M(i, j).to_int64();
        if (dot) return false;
      }
      return true;
    };
    std::vector<long long> two_goods_money(15, 0), claims(15, 0), knowledge(15, 0);
    for (std::size_t p = 0; p < 3; p++) two_goods_money[idx(p, 0)] = 2, two_goods_money[idx(p, 1)] = 1, knowledge[idx(p, 2)] = 1;
    claims[idx(0, 3)] = 1, claims[idx(2, 4)] = -1;
    s.require("conserved: 2 goods + money everywhere (one good makes two of money)", in_lattice(two_goods_money));
    s.require("and claims less debts, not their number", in_lattice(claims) && !in_lattice([&] { auto c = claims; c[idx(2, 4)] = 1; return c; }()));
    s.require("knowledge copied is conserved in no count of holders", !in_lattice(knowledge) && no_knowledge);
    calx::random::Stream rng(9, 9, 9, 9);
    std::vector<double> j(4), r(3), u(15, 0.0);
    for (auto& x : j) x = rng.normal();
    for (auto& x : r) x = std::fabs(rng.normal());
    u[idx(1, 0)] = 0.7;
    const auto dm = R.rate(j, r, u);
    double dw = 0;
    for (std::size_t i = 0; i < 15; i++) dw += two_goods_money[i] * dm[i];
    s.check("d(wᵀm)/dt = wᵀu for a conserved w", dw - 2 * 0.7, 1e-12);
  }

  // ---- walkers as a reaction network ----
  {
    namespace pc = physis::chemistry;
    pc::Network net;
    net.species = {"X0", "X1", "X2"};
    net.add("0→1", {{"X0", 1}}, {{"X1", 1}}, pc::Rational(3, 2));
    net.add("1→0", {{"X1", 1}}, {{"X0", 1}}, pc::Rational(1, 2));
    net.add("1→2", {{"X1", 1}}, {{"X2", 1}}, pc::Rational(1, 1));
    net.add("2→1", {{"X2", 1}}, {{"X1", 1}}, pc::Rational(2, 1));
    const pc::Fock F{3, 4};
    const auto Lop = pc::network_operator(F, net);
    const auto S = pc::reachable(net, {4, 0, 0});
    const auto Q = pc::generator<double>(net, S);
    double worst = 0, drift = 0;
    for (std::size_t x = 0; x < S.states.size(); x++) {
      const std::size_t kx = F.index({S.states[x][0], S.states[x][1], S.states[x][2]});
      for (std::size_t y = 0; y < S.states.size(); y++) {
        if (y == x) continue;
        const auto it = Q[x].find(y);
        worst = std::max(worst, std::fabs(Lop(F.index({S.states[y][0], S.states[y][1], S.states[y][2]}), kx) - (it == Q[x].end() ? 0.0 : it->second)));
      }
      // The first moment: ⟨1| a₁† a₁ L |n⟩ against the mean-field flow into junction 1.
      const auto N1 = pc::creation(F, 1) * pc::annihilation(F, 1) * Lop;
      double lhs = 0;
      for (std::size_t m = 0; m < F.size(); m++) lhs += N1(m, kx);
      const double n0 = S.states[x][0], n1 = S.states[x][1], n2 = S.states[x][2];
      drift = std::max(drift, std::fabs(lhs - (1.5 * n0 - 0.5 * n1 - 1.0 * n1 + 2.0 * n2)));
    }
    s.require("four walkers on three junctions: 15 states", S.states.size() == 15);
    s.check("walking as a reaction network: the Doi–Peliti operator is the master equation", worst, 1e-12);
    s.check("and its first moment is the mean-field flow, exactly (first order)", drift, 1e-12);
  }

  return s.done();
}
