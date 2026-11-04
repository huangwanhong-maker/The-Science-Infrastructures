// Checks of the frames of the body (frames.hpp).
//
// - LPS and RAS are right-handed (L × P = S, R × A = S), and the conversion between them is an involution
//   of determinant +1 that keeps superior and swaps left with right, anterior with posterior.
// - Every display is a right-handed basis, u × v = w, of unit vectors.
// - The radiological transverse section is seen from the feet, with the subject's left on the image's right
//   and anterior up; the neurological one from above, the subject's right on the right; the radiological
//   coronal section from the front, the neurological from behind; the sagittal section from the subject's
//   left, anterior on the image's left.
// - A point's place relative to another, and medial or lateral to it.

#include <cmath>
#include <vector>

#include <anatomy/frames/frames.hpp>
#include <calx/verification/check.hpp>

using calx::verification::Suite;
namespace af = anatomy::frames;
using af::Direction;
using af::Vec3;

namespace {
Vec3 cross(const Vec3& a, const Vec3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
double distance(const Vec3& a, const Vec3& b) { return std::hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]); }
}  // namespace

int main() {
  Suite s;
  {
    using F = af::Frame;
    s.check("LPS is right-handed: L × P = S", distance(cross(af::unit(Direction::left, F::LPS), af::unit(Direction::posterior, F::LPS)), af::unit(Direction::superior, F::LPS)), 0);
    s.check("RAS is right-handed: R × A = S", distance(cross(af::unit(Direction::right, F::RAS), af::unit(Direction::anterior, F::RAS)), af::unit(Direction::superior, F::RAS)), 0);
    bool same = true;
    for (Direction d : {Direction::left, Direction::right, Direction::anterior, Direction::posterior, Direction::superior, Direction::inferior})
      same = same && distance(af::convert(af::unit(d, F::LPS), F::LPS, F::RAS), af::unit(d, F::RAS)) == 0 &&
             distance(af::unit(af::opposite(d), F::LPS), Vec3{-af::unit(d, F::LPS)[0], -af::unit(d, F::LPS)[1], -af::unit(d, F::LPS)[2]}) == 0;
    s.require("each direction converts from LPS to RAS as itself, and its opposite is its negative", same);
    const Vec3 p{12.5, -40, 1300};
    s.check("the conversion is an involution", distance(af::convert(af::convert(p, F::LPS, F::RAS), F::RAS, F::LPS), p), 0);
    const Vec3 x = af::convert({1, 0, 0}, F::LPS, F::RAS), y = af::convert({0, 1, 0}, F::LPS, F::RAS), z = af::convert({0, 0, 1}, F::LPS, F::RAS);
    s.check("of determinant +1, a turn by π about the superior axis", x[0] * (y[1] * z[2] - y[2] * z[1]) - x[1] * (y[0] * z[2] - y[2] * z[0]) + x[2] * (y[0] * z[1] - y[1] * z[0]) - 1, 0);
  }
  {
    double worst = 0;
    for (auto p : {af::Plane::sagittal, af::Plane::coronal, af::Plane::transverse})
      for (auto c : {af::Convention::radiological, af::Convention::neurological}) {
        const auto V = af::view(p, c);
        worst = std::max(worst, distance(cross(V.u, V.v), V.w));
        worst = std::max(worst, std::fabs(std::hypot(V.u[0], V.u[1], V.u[2]) - 1) + std::fabs(std::hypot(V.v[0], V.v[1], V.v[2]) - 1));
        const int a = af::axis(p);
        worst = std::max(worst, 1 - std::fabs(V.w[a]));
      }
    s.check("every display is right-handed, u × v = w, its w along the plane's axis", worst, 0);
    const auto rt = af::view(af::Plane::transverse, af::Convention::radiological);
    const auto nt = af::view(af::Plane::transverse, af::Convention::neurological);
    s.require("the radiological transverse section is seen from the feet, the subject's left on the right, anterior up",
              af::seen_from(rt) == Direction::inferior && rt.right == Direction::left && rt.up == Direction::anterior);
    s.require("the neurological transverse section is seen from above, the subject's right on the right",
              af::seen_from(nt) == Direction::superior && nt.right == Direction::right && nt.up == Direction::anterior);
    s.require("the radiological coronal section is seen from the front, the neurological from behind",
              af::seen_from(af::view(af::Plane::coronal, af::Convention::radiological)) == Direction::anterior &&
                  af::seen_from(af::view(af::Plane::coronal, af::Convention::neurological)) == Direction::posterior);
    const auto sg = af::view(af::Plane::sagittal);
    s.require("the sagittal section is seen from the subject's left, anterior on the image's left, superior up",
              af::seen_from(sg) == Direction::left && sg.right == Direction::posterior && sg.up == Direction::superior);
    s.require("planes are read by their names, axial and frontal among them",
              af::plane_named("axial") == af::Plane::transverse && af::plane_named("frontal") == af::Plane::coronal && af::plane_named("sagittal") == af::Plane::sagittal);
  }
  {
    const auto r = af::relative({30, -20, 1000}, {0, 0, 1100}, 5);
    s.require("a point 30 mm to the left, 20 mm anterior and 100 mm inferior of another", r == std::vector<Direction>{Direction::left, Direction::anterior, Direction::inferior});
    s.require("differences within the tolerance are not named", af::relative({2, 1, 1000}, {0, 0, 1000}, 5).empty());
    s.require("medial and lateral are read from the median plane, on either side",
              af::mediolateral({-20, 0, 0}, {60, 0, 0}, 0, 1) == -1 && af::mediolateral({-80, 0, 0}, {60, 0, 0}, 0, 1) == +1 && af::mediolateral({-60, 0, 0}, {60, 0, 0}, 0, 1) == 0);
  }
  return s.done();
}
