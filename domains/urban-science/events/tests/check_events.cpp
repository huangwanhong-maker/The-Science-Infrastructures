// Checks of the events module (record.hpp, complex.hpp, network.hpp): the record of a city's events, its event complex
// and refinement on morphe's worldline complexes, the cut, and a district's volumes and interfaces.
//
// Known by hand.
// - The hinge: two buildings A and B sharing a wall, both fronting a street S, the three interfaces alive through the
//   window and no event between: K has the two boundary events, three segments and three lenses, the theta-sphere
//   (Betti 1 0 1); K̂ six presences, nine edges, three squares and the two triangles at the ends, a prism's surface
//   (Betti 1 0 1).
// - An alteration of A cuts A's worldline and no lens: one vertex and one edge more, the same homology (a
//   refinement). A regime at the frontage A–S cuts its lens: one vertex, two edges and one face more, the same
//   homology. An observation of A, attended by A alone, is a refinement too.
// - A building C born at an event with its wall to B and its frontage to S: C's worldline begins there; the cut
//   before it has three couplings and three interfaces, after it four and five.
// - Faulty records refused: an interface outliving its building, an event touching a building after its end, an
//   event outside the window, an interface whose two sides are one.
// - Global slicing (the design's Theorem 2.8, its plaquettes aside): ticks attended by every coupling and meeting no
//   interface give K̂ the couplings times the path of the instants, the interfaces' squares only subdivided along
//   their sides, and the homology of the record with no tick (unsubdividing is a refinement).
// - The volumes and interfaces of a small district: two buildings standing against each other on one side of a
//   street, a third across it; the party wall, three frontages and the street's two cross-sections, each primitive
//   joining the right two couplings, the shared wall's normal from the first building toward the second.

#include <array>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <urban/events/complex.hpp>
#include <urban/events/memory.hpp>
#include <urban/events/network.hpp>
#include <urban/events/record.hpp>

using calx::verification::Suite;
namespace ue = urban::events;
namespace sm = physis::statistical_mechanics;

namespace {

std::string text(const ue::Topology& T) {
  std::string s = std::to_string(T.cells[0]) + " " + std::to_string(T.cells[1]) + " " + std::to_string(T.cells[2]) + " | " + std::to_string(T.b0);
  s += " " + (T.b1 ? std::to_string(*T.b1) : std::string("?")) + " " + (T.b2 ? std::to_string(*T.b2) : std::string("?"));
  return s;
}

// Two buildings sharing a wall, both fronting one street, over 2002–2024.
ue::Record hinge() {
  ue::Record R;
  R.window = {ue::mjd(2002), ue::mjd(2025)};
  R.couplings = {{"A", "building", "", {0, 0}}, {"B", "building", "", {10, 0}}, {"S", "street", "", {5, -12}}};
  R.primitives = {{"AB", "party-wall", 0, 1}, {"AS", "frontage", 0, 2}, {"BS", "frontage", 1, 2}};
  return R;
}

bool valid(const morphe::topology::CWComplex& C) {
  try {
    C.validate();
  } catch (const std::exception&) {
    return false;
  }
  return true;
}

bool refused(const ue::Record& R) {
  try {
    (void)ue::construct(R);
  } catch (const std::invalid_argument&) {
    return true;
  }
  return false;
}

}  // namespace

int main() {
  Suite s;
  // The calendar.
  s.require("MJD of 2000-01-01 is 51544, of 1858-11-17 is 0", ue::mjd(2000) == 51544 && ue::mjd(1858, 11, 17) == 0);
  // The hinge.
  {
    const auto R = hinge();
    s.require("the hinge is a sound record", ue::problems(R).empty());
    const auto X = ue::construct(R);
    s.require("its event complex is the theta-sphere: 2 events, 3 segments, 3 lenses, Betti 1 0 1", text(ue::topology(X.K.cw)) == "2 3 3 | 1 0 1", text(ue::topology(X.K.cw)));
    s.require("its refinement a prism's surface: 6 presences, 9 edges, 3 squares and 2 triangles, Betti 1 0 1",
              text(ue::topology(X.refined.cw)) == "6 9 5 | 1 0 1" && X.refined.triangles.size() == 2, text(ue::topology(X.refined.cw)));
    s.require("both complexes valid (closed words, ∂∂ = 0)", valid(X.K.cw) && valid(X.refined.cw));
  }
  // Events that refine and events that cut.
  {
    auto R = hinge();
    R.events.push_back({"alter A", "alteration", "refurbished", ue::year_span(2010), {0}});
    const auto X = ue::construct(R);
    s.require("an alteration of A: 3 vertices, 4 segments, the same 3 lenses, Betti 1 0 1", text(ue::topology(X.K.cw)) == "3 4 3 | 1 0 1", text(ue::topology(X.K.cw)));
    auto Q = hinge();
    Q.events.push_back({"close AS", "regime", "closed", {ue::mjd(2015, 6, 1), ue::mjd(2015, 6, 2)}, {}, {1}});
    const auto Y = ue::construct(Q);
    s.require("a regime at the frontage A–S: 3 vertices, 5 segments, 4 lenses, Betti 1 0 1", text(ue::topology(Y.K.cw)) == "3 5 4 | 1 0 1", text(ue::topology(Y.K.cw)));
    s.require("... and A and S attend it, B does not", Y.encounters.present[1] == std::vector<std::size_t>{0, 2});
    auto O = hinge();
    for (int y = 2003; y <= 2023; y += 2) O.events.push_back({"census " + std::to_string(y), "observation", "census", ue::year_span(y), {0}});
    const auto Z = ue::construct(O);
    s.require("eleven observations of A refine its worldline and keep the homology (Betti 1 0 1)", ue::topology(Z.K.cw).b2 == 1 && ue::topology(Z.K.cw).b1 == 0 && Z.K.segments_of[0].size() == 12);
  }
  // A birth, and the cut.
  {
    auto R = hinge();
    R.couplings.push_back({"C", "building", "", {20, 0}});
    R.primitives.push_back({"BC", "party-wall", 1, 3});
    R.primitives.push_back({"CS", "frontage", 3, 2});
    ue::Event born{"C built", "birth", "constructed", ue::year_span(2012), {1}};
    born.begins_couplings = {3}, born.begins_primitives = {3, 4};
    R.events.push_back(born);
    s.require("a building born with its wall and frontage: a sound record", ue::problems(R).empty(), ue::problems(R).empty() ? "" : ue::problems(R).front());
    const auto X = ue::construct(R);
    s.require("C's worldline begins at its birth", X.K.worldline[3].front() == X.vertex_of[0] && X.life_c[3][0] == 1);
    const auto before = ue::cut(X, ue::mjd(2008)), after = ue::cut(X, ue::mjd(2016));
    s.require("the cut before it: 3 couplings, 3 interfaces; after it: 4 and 5",
              before.couplings.size() == 3 && before.primitives.size() == 3 && after.couplings.size() == 4 && after.primitives.size() == 5);
    s.require("B and S attend C's birth (their interfaces begin there), A does not", X.encounters.present[1] == std::vector<std::size_t>{1, 2, 3});
    const auto T = ue::topology(X.K.cw);
    s.require("the complex with C: valid, connected, its Euler characteristic V − E + F", valid(X.K.cw) && T.b0 == 1 && T.euler == 3 - 6 + 5, text(T));
  }
  // Faulty records.
  {
    auto R = hinge();
    ue::Event gone{"A demolished", "death", "demolished", ue::year_span(2010), {0}};
    gone.ends_couplings = {0};
    R.events.push_back(gone);
    s.require("an interface outliving its building is refused", refused(R));
    auto Q = hinge();
    ue::Event g2{"A demolished", "death", "demolished", ue::year_span(2010), {0}};
    g2.ends_couplings = {0}, g2.ends_primitives = {0, 1};
    Q.events.push_back(g2);
    s.require("... and accepted when the interfaces end with it", ue::problems(Q).empty());
    Q.events.push_back({"A altered", "alteration", "refurbished", ue::year_span(2014), {0}});
    s.require("an event touching a building after its end is refused", refused(Q));
    auto O = hinge();
    O.events.push_back({"late", "occurrence", "", ue::year_span(2030), {0}});
    s.require("an event outside the window is refused", refused(O));
    auto L = hinge();
    L.primitives.push_back({"AA", "party-wall", 0, 0});
    s.require("an interface whose two sides are one is refused", refused(L));
  }
  // Global slicing (Theorem 2.8 without plaquettes).
  {
    auto R = hinge();
    const int n = 7;
    for (int k = 0; k < n; k++) R.events.push_back({"tick " + std::to_string(k), "tick", "clock", ue::year_span(2003 + 3 * k), {0, 1, 2}});
    const auto X = ue::construct(R);
    const auto T = ue::topology(X.refined.cw);
    // Each coupling at each of n + 2 instants; n + 1 periods of each worldline and the rungs at the two ends only
    // (the ticks meet no interface); the squares of the interfaces are not cut, their boundaries subdivided.
    const bool counts = T.cells[0] == 3u * (n + 2) && T.cells[1] == 3u * (n + 1) + 3u * 2 && T.cells[2] == 3u + 2u;
    s.require("seven global ticks: K̂'s cells those of the hinge's network times the path of its instants", counts, text(T));
    s.require("... and its homology that of the hinge with no event (Betti 1 0 1)", T.b0 == 1 && T.b1 == 0 && T.b2 == 1, text(T));
  }
  // The volumes and interfaces of a small district.
  {
    namespace as = architecture::space;
    namespace ul = urban::landscape;
    // A street along y = 0 from x = −40 to 40, with junctions at its ends; two buildings side by side north of it,
    // a third south of it.
    ul::Network net;
    net.junctions = {{-40, 0}, {40, 0}};
    net.segments = {{0, 1, {}, "street", "Main Street"}};
    const auto D = ul::district(net);
    std::vector<as::Footprint> F = {
        {"a", "A", {{-20, 5}, {0, 5}, {0, 25}, {-20, 25}}, 0, 12},
        {"b", "B", {{0, 5}, {20, 5}, {20, 25}, {0, 25}}, 0, 18},
        {"c", "C", {{-10, -25}, {10, -25}, {10, -6}, {-10, -6}}, 0, 9},
    };
    const auto M = as::massing(F);
    const auto N = ue::city_network(M, F, D);
    s.require("three buildings, one street space and two junction spaces", N.buildings == 3 && N.streets == 1 && N.junctions == 2);
    s.require("one party wall, three frontages, two cross-sections", N.party_walls == 1 && N.frontages == 3 && N.cross_sections == 2,
              std::to_string(N.party_walls) + " " + std::to_string(N.frontages) + " " + std::to_string(N.cross_sections));
    bool ends = true;
    std::set<std::string> fronted;
    for (const auto& p : N.primitives) {
      const auto& A = N.couplings[p.a];
      const auto& B = N.couplings[p.b];
      if (p.type == "party-wall") {
        ends = ends && A.type == "building" && B.type == "building";
        // The wall at x = 0 seen from the first building toward the second: along ±x, its area 20 m × 12 m (the lower).
        const double dir = (B.at[0] - A.at[0]) * p.normal[0];
        ends = ends && dir > 0 && std::fabs(p.area - 240) < 1e-9 && std::fabs(std::fabs(p.normal[0]) - 1) < 1e-9;
      } else if (p.type == "frontage") {
        ends = ends && A.type == "building" && B.type == "street";
        fronted.insert(A.name);
        // Each building fronts the street with the one facade facing it, 20 m wide.
        const double height = A.name == "A" ? 12 : A.name == "B" ? 18 : 9;
        ends = ends && std::fabs(p.area - 20 * height) < 1e-9;
      } else if (p.type == "cross-section") {
        ends = ends && A.type == "street" && B.type == "junction";
      }
    }
    s.require("each interface joins the right two volumes, with its area; the party wall's normal from the first building to the second", ends);
    s.require("every building fronts the street", fronted == std::set<std::string>{"A", "B", "C"});
    // The district's record with no event: a sound record, its K connected.
    ue::Record R;
    R.window = {ue::mjd(2002), ue::mjd(2025)};
    R.couplings = N.couplings, R.primitives = N.primitives;
    const auto X = ue::construct(R);
    const auto T = ue::topology(X.K.cw);
    s.require("the district's record with no event: its event complex connected and valid", T.b0 == 1 && valid(X.K.cw), text(T));
  }
  // What a city remembers of its uses.
  {
    const std::map<std::string, double> counts = {{"Retail - Shop", 30}, {"Retail - Showroom", 10}, {"Office", 40}, {"Unoccupied - Unused", 20}};
    const auto U = ue::prefix_classification(counts);
    s.require("CLUE's uses by their own prefixes: three families (Office, Retail, Unoccupied) over four uses, the law their frequencies",
              U.uses.size() == 4 && U.C.levels() == 2 && std::fabs(U.C.law[U.of("Office")] - 0.4) < 1e-12
                && U.C.block[1][U.of("Retail - Shop")] == U.C.block[1][U.of("Retail - Showroom")] && U.C.block[1][U.of("Office")] != U.C.block[1][U.of("Retail - Shop")]);
    const auto A = ue::anzsic_classification({{"4511", "H", "45", "451"}, {"4512", "H", "45", "451"}, {"4520", "H", "45", "452"}, {"6932", "M", "69", "693"}, {"0000", "", "00", "000"}},
                                             {{"4511", 5}, {"4512", 3}, {"4520", 2}, {"6932", 6}, {"0000", 4}});
    s.require("ANZSIC: division, subdivision, group, class, vacant space a division of its own",
              A.C.levels() == 4 && A.C.block[1][A.of("0000")] != A.C.block[1][A.of("4511")] && A.C.block[3][A.of("4511")] == A.C.block[3][A.of("4512")]
                && A.C.block[3][A.of("4520")] != A.C.block[3][A.of("4511")] && A.C.block[2][A.of("4520")] == A.C.block[2][A.of("4511")]);
    // A panel simulated from known rates: 400 units over 22 censuses.
    const std::vector<double> truth = {0.04, 0.5};
    const sm::NestedChain X(U.C, truth);
    calx::random::Stream rng(20261104, 21, 0, 0);
    const auto draw = [&](const std::vector<double>& w) { double u = rng.uniform(), t = 0; for (std::size_t i = 0; i < w.size(); i++) if ((t += w[i]) >= u) return i; return w.size() - 1; };
    std::vector<ue::Observation> O;
    for (int k = 0; k < 400; k++) {
      std::size_t x = draw(U.C.law);
      for (int y = 2002; y <= 2023; y++) {
        O.push_back({"unit " + std::to_string(k), U.uses[x], y});
        std::vector<double> w(U.uses.size());
        for (std::size_t z = 0; z < w.size(); z++) w[z] = X.kernel(x, z, 1.0);
        x = draw(w);
      }
    }
    const auto T = ue::transitions_of(U, O);
    double n = 0;
    for (const auto& t : T) n += t.count;
    const auto M = ue::remember(U, T, 400);
    s.require("a simulated panel: 400 units, 21 transitions each, aggregated by pair", std::fabs(n - 400 * 21) < 1e-9);
    s.require("... its rates recovered within three standard errors, the family kept longer than the use",
              std::fabs(M.fit.lambda[0] - truth[0]) < 3 * M.fit.lambda_error[0] && std::fabs(M.fit.lambda[1] - truth[1]) < 3 * M.fit.lambda_error[1] && M.half_life[0] > M.half_life[1],
              std::to_string(M.fit.lambda[0]) + " " + std::to_string(M.fit.lambda[1]));
    s.require("... and the flat chain rejected (the levels are remembered apart)", M.ratio > 20 && M.degrees == 1, std::to_string(M.ratio));
    // The chains on the complex: the hinge's two buildings observed at four censuses, one refurbished between them.
    auto R = hinge();
    const std::vector<std::pair<int, std::array<const char*, 2>>> census = {{2004, {"Office", "Retail - Shop"}}, {2008, {"Office", "Retail - Showroom"}}, {2012, {"Unoccupied - Unused", "Retail - Showroom"}}, {2020, {"Office", "Office"}}};
    for (const auto& [y, uses] : census)
      for (std::size_t c = 0; c < 2; c++) R.events.push_back({"census " + std::to_string(y) + std::to_string(c), "observation", "census", ue::year_span(y), {c}});
    R.events.push_back({"refurbished", "alteration", "refurbished", ue::year_span(2010), {0}});
    const auto XC = ue::construct(R);
    std::map<std::size_t, std::map<std::size_t, std::size_t>> seen;
    for (std::size_t e = 0; e < R.events.size(); e++) {
      if (R.events[e].type != "observation") continue;
      const std::size_t c = R.events[e].couplings[0];
      const int y = std::stoi(R.events[e].id.substr(7, 4));
      for (const auto& [yy, uses] : census)
        if (yy == y) seen[c][XC.vertex_of[e]] = U.of(uses[c]);
    }
    const auto chains = ue::chains_on_complex(XC, X, seen);
    const double z = sm::contract(chains).log_z;
    double closed = 0;
    for (std::size_t c = 0; c < 2; c++) {
      closed += std::log(U.C.law[U.of(census[0].second[c])]);
      for (std::size_t k = 0; k + 1 < census.size(); k++)
        closed += std::log(X.kernel(U.of(census[k].second[c]), U.of(census[k + 1].second[c]), (ue::mjd(census[k + 1].first) - ue::mjd(census[k].first)) / 365.2425));
    }
    s.require("the buildings' uses as chains on the event complex: Z by the contraction equals the product of their transitions (the refurbishment a refinement)",
              std::fabs(z - closed) < 1e-10, std::to_string(z) + " vs " + std::to_string(closed));
    // The rates of the record's kinds of events.
    const auto rates = ue::event_rates(R, XC);
    const auto ref = std::find_if(rates.begin(), rates.end(), [](const ue::EventRate& r) { return r.kind == "refurbished"; });
    s.require("one refurbishment over two buildings' 23 years: 1/46 a building-year", ref != rates.end() && std::fabs(ref->rate - 1.0 / 46) < 1e-3 && ref->of == "building",
              ref == rates.end() ? "none" : std::to_string(ref->rate));
  }
  return s.done();
}
