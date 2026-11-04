// The frames of the body: the anatomical directions, the patient frames LPS and RAS, the cardinal planes and
// how a section through each is displayed (ADR 0066; theory/anatomy.md §2).
//
// The anatomical position. The body stands upright, facing forward, arms at the sides with the palms
// forward. Its directions are named from the body itself, never from an observer: left and right are the
// subject's own; anterior (ventral) is the side of the face, posterior (dorsal) the back; superior (cranial)
// toward the head, inferior (caudal) toward the feet (Terminologia Anatomica, 2nd ed., FIPAT 2019, §1.1;
// Standring 2020, Gray's Anatomy, 42nd ed., ch. 1).
//
// Patient frames. Two right-handed Cartesian frames are in use, named by the directions of +x, +y, +z:
// - LPS: +x left, +y posterior, +z superior. DICOM's patient coordinate system (PS3.3 2024, C.7.6.2.1.1),
//   used by CT and MRI scanners and ITK, and the frame of BodyParts3D's meshes (the database's
//   coordinate_system.png: left x > 0, posterior y > 0, superior z greater).
// - RAS: +x right, +y anterior, +z superior. NIfTI's scanner and template spaces, MNI and Talairach, used by
//   the Brain Atlas.
// Both are right-handed (L × P = S and R × A = S), and the one is the other turned by π about the superior
// axis: (x, y, z)_RAS = (−x, −y, z)_LPS, an involution of determinant +1. Conversion is explicit, here, and
// never implicit (CLAUDE.md §4.1).
//
// The cardinal planes. A sagittal plane is normal to left–right (the median plane is the sagittal plane
// through the midline), a coronal (frontal) plane normal to anterior–posterior, a transverse (axial, or
// horizontal in the anatomical position) plane normal to superior–inferior.
//
// Displaying a section. An image of a section is drawn with u to the right and v upward, and the viewer
// looks at it from the side w = u × v. Two conventions are in use (Hounsfield's scanners set the first):
// - radiological: transverse sections seen from the feet (the subject's left on the image's right,
//   anterior up), coronal sections seen from the front (the subject's left on the image's right);
// - neurological: transverse sections seen from above, coronal from behind (the subject's right on the
//   image's right).
// Sagittal sections are seen from the subject's left in both, anterior on the image's left. Each display is
// a basis (u, v, w) of the patient frame with u × v = w; a section cut by physis's plane_section with the
// normal w has its loops counter-clockwise in (u, v) round the solid, so that their shoelace areas are the
// areas of the pieces cut. Checks: tests/check_frames.cpp.

#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace anatomy::frames {

using Vec3 = std::array<double, 3>;

enum class Direction { left, right, anterior, posterior, superior, inferior };
enum class Frame { LPS, RAS };
enum class Plane { sagittal, coronal, transverse };
enum class Convention { radiological, neurological };

inline std::string_view name(Direction d) {
  switch (d) {
    case Direction::left: return "left";
    case Direction::right: return "right";
    case Direction::anterior: return "anterior";
    case Direction::posterior: return "posterior";
    case Direction::superior: return "superior";
    case Direction::inferior: return "inferior";
  }
  return "";
}
inline std::string_view name(Plane p) {
  switch (p) {
    case Plane::sagittal: return "sagittal";
    case Plane::coronal: return "coronal";
    case Plane::transverse: return "transverse";
  }
  return "";
}
inline Plane plane_named(std::string_view s) {
  if (s == "sagittal") return Plane::sagittal;
  if (s == "coronal" || s == "frontal") return Plane::coronal;
  if (s == "transverse" || s == "axial" || s == "horizontal") return Plane::transverse;
  throw std::invalid_argument("anatomy::frames::plane_named: sagittal, coronal or transverse, not " + std::string(s));
}

inline Direction opposite(Direction d) {
  switch (d) {
    case Direction::left: return Direction::right;
    case Direction::right: return Direction::left;
    case Direction::anterior: return Direction::posterior;
    case Direction::posterior: return Direction::anterior;
    case Direction::superior: return Direction::inferior;
    case Direction::inferior: return Direction::superior;
  }
  return d;
}

// The unit vector of a direction, in a frame.
inline Vec3 unit(Direction d, Frame f) {
  Vec3 lps{};
  switch (d) {
    case Direction::left: lps = {1, 0, 0}; break;
    case Direction::right: lps = {-1, 0, 0}; break;
    case Direction::posterior: lps = {0, 1, 0}; break;
    case Direction::anterior: lps = {0, -1, 0}; break;
    case Direction::superior: lps = {0, 0, 1}; break;
    case Direction::inferior: lps = {0, 0, -1}; break;
  }
  return f == Frame::LPS ? lps : Vec3{-lps[0], -lps[1], lps[2]};
}

// A point or vector carried from one patient frame to the other.
inline Vec3 convert(const Vec3& p, Frame from, Frame to) {
  if (from == to) return p;
  return {-p[0], -p[1], p[2]};
}

// The axis normal to a cardinal plane, as the index of its coordinate (the same in LPS and RAS).
inline int axis(Plane p) { return p == Plane::sagittal ? 0 : p == Plane::coronal ? 1 : 2; }

// The display of a section: u to the right, v upward, w = u × v toward the viewer; LPS components.
struct View {
  Plane plane = Plane::transverse;
  Convention convention = Convention::radiological;
  Vec3 u{}, v{}, w{};
  Direction right{}, up{};  // the directions the image's right and top show
};

inline View view(Plane p, Convention c = Convention::radiological) {
  View V;
  V.plane = p;
  V.convention = c;
  const bool radiological = c == Convention::radiological;
  switch (p) {
    case Plane::transverse:
      V.right = radiological ? Direction::left : Direction::right;
      V.up = Direction::anterior;
      break;
    case Plane::coronal:
      V.right = radiological ? Direction::left : Direction::right;
      V.up = Direction::superior;
      break;
    case Plane::sagittal:
      V.right = Direction::posterior;
      V.up = Direction::superior;
      break;
  }
  V.u = unit(V.right, Frame::LPS);
  V.v = unit(V.up, Frame::LPS);
  V.w = {V.u[1] * V.v[2] - V.u[2] * V.v[1], V.u[2] * V.v[0] - V.u[0] * V.v[2], V.u[0] * V.v[1] - V.u[1] * V.v[0]};
  return V;
}

// Where the viewer stands: the direction w points to.
inline Direction seen_from(const View& V) {
  for (Direction d : {Direction::left, Direction::right, Direction::anterior, Direction::posterior, Direction::superior, Direction::inferior}) {
    const Vec3 e = unit(d, Frame::LPS);
    if (e[0] * V.w[0] + e[1] * V.w[1] + e[2] * V.w[2] > 0.5) return d;
  }
  throw std::logic_error("anatomy::frames::seen_from: not a cardinal view");
}

// Where a point lies from a reference, on each axis whose difference exceeds a tolerance (mm): at most one
// of left/right, anterior/posterior and superior/inferior, in that order. Both points in LPS.
inline std::vector<Direction> relative(const Vec3& p, const Vec3& reference, double tolerance) {
  std::vector<Direction> out;
  const Vec3 d{p[0] - reference[0], p[1] - reference[1], p[2] - reference[2]};
  if (std::fabs(d[0]) > tolerance) out.push_back(d[0] > 0 ? Direction::left : Direction::right);
  if (std::fabs(d[1]) > tolerance) out.push_back(d[1] > 0 ? Direction::posterior : Direction::anterior);
  if (std::fabs(d[2]) > tolerance) out.push_back(d[2] > 0 ? Direction::superior : Direction::inferior);
  return out;
}

// Medial or lateral: nearer to or farther from the median plane x = midline than a reference, by more
// than a tolerance (mm); 0 when neither.
inline int mediolateral(const Vec3& p, const Vec3& reference, double midline, double tolerance) {
  const double a = std::fabs(p[0] - midline), b = std::fabs(reference[0] - midline);
  return a < b - tolerance ? -1 : a > b + tolerance ? +1 : 0;  // −1 medial, +1 lateral
}

}  // namespace anatomy::frames
