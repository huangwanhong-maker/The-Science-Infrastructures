// Sections of the body: each part cut by a cardinal plane, in the coordinates of the section's display, and
// the vertebral level of a transverse plane (ADR 0066; theory/anatomy.md §5).
//
// A section is cut by physis's plane_section with the display's normal w (frames.hpp), so that the loops
// round each piece run counter-clockwise in the image's (u, v) and the holes clockwise; a piece's area is the
// sum of its loops' shoelace areas (Green's theorem), which is the part's cross-sectional area there. A part
// whose mesh is not closed is cut into open chains, drawn as lines, with no area. Only parts whose bounds
// straddle the plane are cut.
//
// Cross-sectional anatomy names a transverse plane by the vertebra at its height, as atlases and imaging do
// (the sternal angle at T4/T5, the transpyloric plane at L1; Standring 2020, ch. 52 and 59). The level of a
// height is the vertebra whose body it passes through, read here from the body itself: the centres of the
// vertebrae's meshes, each vertebra the interval between the midpoints to its neighbours. Checks:
// tests/check_sections.cpp.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include <anatomy/frames/frames.hpp>
#include <physis/computational_geometry/surfaces.hpp>

namespace anatomy::sections {

namespace cg = physis::computational_geometry;

struct Piece {
  std::size_t part = 0;                         // the index of the surface cut
  std::vector<std::vector<cg::Point2>> loops;   // closed, counter-clockwise round the solid
  std::vector<std::vector<cg::Point2>> chains;  // open, where the mesh is not closed
  double area = 0;                              // mm², of the closed loops
};

// The plane of a view at a level: the coordinate along the plane's axis in LPS (x for sagittal, y for
// coronal, z for transverse), mm.
inline cg::Plane3 plane(const frames::View& V, double level) {
  cg::Plane3 P;
  P.origin = {0, 0, 0};
  P.origin[frames::axis(V.plane)] = level;
  P.normal = V.w;
  return P;
}

// The sections of the surfaces whose bounds straddle the plane; `cut[k]` false leaves surface k out.
inline std::vector<Piece> section(const std::vector<const cg::TriangleSurface*>& surfaces, const std::vector<cg::Box3>& bounds,
                                  const frames::View& V, double level, const std::vector<char>& cut = {}) {
  const cg::Plane3 P = plane(V, level);
  const int a = frames::axis(V.plane);
  std::vector<Piece> out;
  for (std::size_t k = 0; k < surfaces.size(); k++) {
    if (!surfaces[k] || (!cut.empty() && !cut[k])) continue;
    if (!(bounds[k].lower[a] <= level && level <= bounds[k].upper[a])) continue;
    const auto chains = cg::plane_section(*surfaces[k], P);
    if (chains.empty()) continue;
    Piece piece;
    piece.part = k;
    for (const auto& c : chains) {
      auto q = cg::in_plane(c.points, {0, 0, 0}, V.u, V.v);
      if (c.closed) {
        piece.area += cg::ring_area(q);
        piece.loops.push_back(std::move(q));
      } else {
        piece.chains.push_back(std::move(q));
      }
    }
    out.push_back(std::move(piece));
  }
  return out;
}

// ---- vertebral levels ------------------------------------------------------------------------------------------

struct Vertebra {
  std::string_view label, fma, name;
};

inline constexpr std::array<Vertebra, 25> VERTEBRAE{{
    {"C1", "FMA12519", "atlas"}, {"C2", "FMA12520", "axis"}, {"C3", "FMA12521", "third cervical vertebra"},
    {"C4", "FMA12522", "fourth cervical vertebra"}, {"C5", "FMA12523", "fifth cervical vertebra"},
    {"C6", "FMA12524", "sixth cervical vertebra"}, {"C7", "FMA12525", "seventh cervical vertebra"},
    {"T1", "FMA9165", "first thoracic vertebra"}, {"T2", "FMA9187", "second thoracic vertebra"},
    {"T3", "FMA9209", "third thoracic vertebra"}, {"T4", "FMA9248", "fourth thoracic vertebra"},
    {"T5", "FMA9922", "fifth thoracic vertebra"}, {"T6", "FMA9945", "sixth thoracic vertebra"},
    {"T7", "FMA9968", "seventh thoracic vertebra"}, {"T8", "FMA9991", "eighth thoracic vertebra"},
    {"T9", "FMA10014", "ninth thoracic vertebra"}, {"T10", "FMA10037", "tenth thoracic vertebra"},
    {"T11", "FMA10059", "eleventh thoracic vertebra"}, {"T12", "FMA10081", "twelfth thoracic vertebra"},
    {"L1", "FMA13072", "first lumbar vertebra"}, {"L2", "FMA13073", "second lumbar vertebra"},
    {"L3", "FMA13074", "third lumbar vertebra"}, {"L4", "FMA13075", "fourth lumbar vertebra"},
    {"L5", "FMA13076", "fifth lumbar vertebra"}, {"S", "FMA16202", "sacrum"},
}};

struct Level {
  std::string label;
  double height = 0;  // the centre's z, mm
};

// The level at a height: the vertebra whose interval holds it, the intervals bounded by the midpoints
// between neighbouring centres and, at the ends, by as far again as half the gap to the neighbour. Empty
// above or below the column.
inline std::string level_at(std::vector<Level> levels, double z) {
  if (levels.empty()) return {};
  std::sort(levels.begin(), levels.end(), [](const Level& a, const Level& b) { return a.height > b.height; });
  for (std::size_t k = 0; k < levels.size(); k++) {
    const double up = k == 0 ? levels[0].height + (levels.size() > 1 ? (levels[0].height - levels[1].height) / 2 : 0)
                             : (levels[k - 1].height + levels[k].height) / 2;
    const double down = k + 1 == levels.size() ? levels[k].height - (k > 0 ? (levels[k - 1].height - levels[k].height) / 2 : 0)
                                               : (levels[k].height + levels[k + 1].height) / 2;
    if (z <= up && z >= down) return levels[k].label;
  }
  return {};
}

}  // namespace anatomy::sections
