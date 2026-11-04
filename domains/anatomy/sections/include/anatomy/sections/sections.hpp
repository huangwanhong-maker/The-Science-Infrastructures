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
#include <anatomy/parts/parts.hpp>
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

// The Human Reference Atlas's woman (body-hra-f): her vertebrae by the model's nodes. She has six lumbar
// vertebrae, L1 to L6, as the Visible Human Female's column was modelled (a lumbar sixth, a variant found in a
// few per cent of people), so her levels below L5 are her own.
inline constexpr std::array<Vertebra, 26> REFERENCE_VERTEBRAE{{
    {"C1", "VH_F_cervical_vertebra_1", "cervical vertebra 1"},
    {"C2", "VH_F_cervical_vertebra_2", "cervical vertebra 2"},
    {"C3", "VH_F_cervical_vertebra_3", "cervical vertebra 3"},
    {"C4", "VH_F_cervical_vertebra_4", "cervical vertebra 4"},
    {"C5", "VH_F_cervical_vertebra_5", "cervical vertebra 5"},
    {"C6", "VH_F_cervical_vertebra_6", "cervical vertebra 6"},
    {"C7", "VH_F_cervical_vertebra_7", "cervical vertebra 7"},
    {"T1", "VH_F_thoracic_vertebra_1", "thoracic vertebra 1"},
    {"T2", "VH_F_thoracic_vertebra_2", "thoracic vertebra 2"},
    {"T3", "VH_F_thoracic_vertebra_3", "thoracic vertebra 3"},
    {"T4", "VH_F_thoracic_vertebra_4", "thoracic vertebra 4"},
    {"T5", "VH_F_thoracic_vertebra_5", "thoracic vertebra 5"},
    {"T6", "VH_F_thoracic_vertebra_6", "thoracic vertebra 6"},
    {"T7", "VH_F_thoracic_vertebra_7", "thoracic vertebra 7"},
    {"T8", "VH_F_thoracic_vertebra_8", "thoracic vertebra 8"},
    {"T9", "VH_F_thoracic_vertebra_9", "thoracic vertebra 9"},
    {"T10", "VH_F_thoracic_vertebra_10", "thoracic vertebra 10"},
    {"T11", "VH_F_thoracic_vertebra_11", "thoracic vertebra 11"},
    {"T12", "VH_F_thoracic_vertebra_12", "thoracic vertebra 12"},
    {"L1", "VH_F_lumbar_vertebra_1", "lumbar vertebra 1"},
    {"L2", "VH_F_lumbar_vertebra_2", "lumbar vertebra 2"},
    {"L3", "VH_F_lumbar_vertebra_3", "lumbar vertebra 3"},
    {"L4", "VH_F_lumbar_vertebra_4", "lumbar vertebra 4"},
    {"L5", "VH_F_lumbar_vertebra_5", "lumbar vertebra 5"},
    {"L6", "VH_F_lumbar_vertebra_6", "lumbar vertebra 6"},
    {"S", "VH_F_sacrum", "sacrum"},
}};

// The vertebrae an atlas names: FMA's for BodyParts3D's man, the model's nodes for the Atlas's woman.
inline std::vector<Vertebra> vertebrae_for(const parts::Atlas& A) {
  if (parts::is_reference_atlas(A)) return {REFERENCE_VERTEBRAE.begin(), REFERENCE_VERTEBRAE.end()};
  return {VERTEBRAE.begin(), VERTEBRAE.end()};
}

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
