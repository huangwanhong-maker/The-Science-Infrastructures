// Places, and the faces a map gives a network of them (theory note §3.2).
//
// A transfer network among places has loops but, as a graph, no faces: its holonomies cannot be
// located. Geography supplies faces. The spherical Delaunay triangulation of the places (physis
// computational_geometry) tiles the map with triangles of neighbouring places, each a 2-cell
// bounded by legs between its corners; where two neighbours have no leg, a neutral one is added
// (transport 1, no flow), which changes no holonomy and no flow. On a network so filled:
//
//   the curvature of the transports on a triangle locates the deficits of the loops round it,
//     and the deficit round any loop of neighbours is the sum of the curvatures it encloses;
//   the circulation of flows splits into local eddies round the triangles and global
//     circulation round the holes the map leaves: loops through legs between places that are not
//     neighbours (long-distance links), which no triangle fills.
//
// When the places do not surround the centre of the sphere (a region, not the world), the hull's
// triangles whose planes leave the centre outside cover the empty side of the sphere, not the
// region, and are not faces of the map.

#pragma once

#include <economics/transfer/network.hpp>
#include <physis/computational_geometry/sphere.hpp>

#include <array>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace economics {

struct Place {
  std::string name;
  double latitude = 0, longitude = 0;  // geocentric, radians
};

struct MapFaces {
  std::vector<std::array<std::size_t, 3>> triangles;  // places, anticlockwise seen from outside
  std::size_t neutral_legs = 0;                       // legs added between neighbours with none
  bool covers_sphere = false;
};

// The network with the triangles of its places' spherical Delaunay triangulation as faces (the
// positions of N are the places, in order).
inline TransferNetwork with_map_faces(TransferNetwork N, const std::vector<Place>& places, MapFaces* info = nullptr) {
  namespace cg = physis::computational_geometry;
  if (places.size() != N.positions.size()) throw std::invalid_argument("economics::with_map_faces: one place per position");
  std::vector<cg::Vec3> points;
  for (const auto& p : places) points.push_back(cg::unit_vector(p.latitude, p.longitude));
  const auto T = cg::spherical_delaunay(points);
  for (std::size_t i = 0; i < places.size(); i++)
    if (T.point_of[i] != i) throw std::invalid_argument("economics::with_map_faces: two places at one point");
  // A leg for each side: an existing one either way, else a neutral one.
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> leg_of;
  for (std::size_t e = 0; e < N.legs.size(); e++) {
    const auto key = std::minmax(N.legs[e].from, N.legs[e].to);
    if (N.legs[e].from != N.legs[e].to && !leg_of.count(key)) leg_of[key] = e;
  }
  MapFaces M;
  M.covers_sphere = T.covers_sphere;
  for (const auto& t : T.triangles) {
    // A face of the map only if its plane has the centre inside.
    const cg::Vec3 &a = T.points[t[0]], &b = T.points[t[1]], &c = T.points[t[2]];
    const cg::Vec3 normal = cg::cross({b[0] - a[0], b[1] - a[1], b[2] - a[2]}, {c[0] - a[0], c[1] - a[1], c[2] - a[2]});
    if (cg::dot(normal, a) <= 0) continue;
    Path face;
    for (int k = 0; k < 3; k++) {
      const std::size_t u = t[k], v = t[(k + 1) % 3];
      const auto key = std::minmax(u, v);
      auto it = leg_of.find(key);
      if (it == leg_of.end()) {
        it = leg_of.emplace(key, N.legs.size()).first;
        N.legs.push_back({u, v, 1.0, 0.0, "neighbours"});
        M.neutral_legs++;
      }
      face.push_back({it->second, N.legs[it->second].from == u ? 1 : -1});
    }
    N.faces.push_back(std::move(face));
    M.triangles.push_back(t);
  }
  if (info) *info = M;
  return N;
}

}  // namespace economics
