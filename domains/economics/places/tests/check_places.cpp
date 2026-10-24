// Checks of the faces a map gives a network of places (map.hpp).
//
// Six places at the vertices of an octahedron, with legs only between neighbours: the map fills
// all eight triangles, the network becomes a sphere, and any flow on it is driven or circulates
// locally, never globally; a leg between antipodes adds a hole, round which flow circulates
// globally. The curvature of the transports on the triangles adds up, by Stokes, to the deficit
// round any loop of neighbours, and the neutral legs added between neighbours change no holonomy.
// Places on one side of the world (a region) give no faces on the other side.

#include <economics/circulation/flows.hpp>
#include <economics/places/map.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

using namespace economics;

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(20261024, 0x706c6163u, 0, 0);  // "plac"
  constexpr double pi = std::numbers::pi;
  // The octahedron: ±x, ±y, ±z as places.
  const std::vector<Place> places{{"px", 0, 0}, {"nx", 0, pi}, {"py", 0, pi / 2}, {"ny", 0, -pi / 2}, {"pz", pi / 2, 0}, {"nz", -pi / 2, 0}};
  TransferNetwork N;
  for (const auto& p : places) N.positions.push_back(p.name);
  // Legs between neighbours only (not between antipodes), random transports and flows.
  for (std::size_t u = 0; u < 6; u++)
    for (std::size_t v = u + 1; v < 6; v++)
      if (v != (u ^ 1)) N.legs.push_back({u, v, std::exp(z.uniform() - 0.5), z.uniform() * 10, "goods"});
  MapFaces info;
  const auto M = with_map_faces(N, places, &info);
  s.require("the octahedron's places: eight triangles, no neutral leg, covering the sphere", info.triangles.size() == 8 && info.neutral_legs == 0 && info.covers_sphere);
  const auto C = circulation(M);
  s.check("... a flow on the filled sphere has no global circulation", C.global, 1e-12);
  s.check("... its driven and local shares add to 1", C.driven + C.local - 1, 1e-12);
  // Stokes: the curvatures of the four triangles round the north pole add up to the deficit round
  // the equator, px → py → nx → ny → px.
  const Cochain w = M.transport();
  const Cochain curvature = physis::discretisation::d1(M.complex(), w);
  auto leg = [&](std::size_t u, std::size_t v) {
    for (std::size_t e = 0; e < M.legs.size(); e++) {
      if (M.legs[e].from == u && M.legs[e].to == v) return Path::value_type{e, 1};
      if (M.legs[e].from == v && M.legs[e].to == u) return Path::value_type{e, -1};
    }
    throw std::logic_error("no leg");
  };
  const Path equator{leg(0, 2), leg(2, 1), leg(1, 3), leg(3, 0)};
  double north = 0;
  for (std::size_t f = 0; f < info.triangles.size(); f++) {
    const auto& t = info.triangles[f];
    if (t[0] == 4 || t[1] == 4 || t[2] == 4) north += curvature[f];
  }
  s.check("Stokes: the curvatures round the north pole add up to the integral round the equator",
          north - physis::discretisation::integral(w, equator), 1e-12);
  // A long-distance leg between antipodes adds a hole.
  TransferNetwork L = N;
  L.legs.push_back({0, 1, 1.0, 50, "goods"});
  const auto CL = circulation(with_map_faces(L, places));
  s.require("a leg between antipodes opens a hole, round which flow circulates globally", CL.global > 1e-3, std::to_string(CL.global));
  // Neutral legs: remove one neighbour leg and the map adds it back, changing no holonomy.
  TransferNetwork G = N;
  G.legs.erase(G.legs.begin());  // px – py
  MapFaces gi;
  const auto MG = with_map_faces(G, places, &gi);
  s.require("a missing neighbour leg is added as a neutral leg", gi.neutral_legs == 1 && MG.legs.back().ratio == 1.0 && MG.legs.back().flow == 0.0);
  // A region: places in the northern hemisphere only.
  std::vector<Place> region;
  TransferNetwork R;
  for (int i = 0; i < 12; i++) {
    region.push_back({"r" + std::to_string(i), 0.2 + 1.0 * z.uniform(), 2 * pi * z.uniform()});
    R.positions.push_back(region.back().name);
  }
  MapFaces ri;
  (void)with_map_faces(R, region, &ri);
  s.require("places on one side of the world: fewer than 2V − 4 faces, none across the empty side", !ri.covers_sphere && ri.triangles.size() < 2 * 12 - 4,
            std::to_string(ri.triangles.size()) + " faces");
  return s.done();
}
