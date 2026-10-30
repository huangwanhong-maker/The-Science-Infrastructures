// The Earth's ground as a foot meets it: its class presumed from what CRUST1.0 holds at the place, each
// class's contact law presumed, and a realisation of it, the ground a simulated world stands on (ADR 0059;
// theory/ground.md).
//
// The class. CRUST1.0 describes each 1° cell by its layers (structure/crust.hpp): where it has water the
// place is sea, not ground; where it has ice, the ground is ice (Greenland, Antarctica, the ice caps); where it
// has sediments, the ground is sediment (soil, sand, gravel); else rock. A cell of a degree is the whole of
// what is known: a rock outcrop in a sedimentary basin, or the snow on rock, is not seen.
//
// The laws: physis's compliant contact, f_n = k δⁿ (1 + c δ̇) and −μ f_n tanh(v_t/v_s), its three parameters
// presumed for a rubber foot, each log-normal with a median and a spread of its logarithm:
//
//   class     friction μ          load at 1 mm, k (1 mm)ⁿ    damping c
//   ice       0.15, spread 0.45   30 N, spread 0.4          1.5 s/m, spread 0.3
//   sediment  0.60, spread 0.25   15 N, spread 0.6          2.0 s/m, spread 0.3
//   rock      0.80, spread 0.15   40 N, spread 0.4          1.5 s/m, spread 0.3
//
// All are assumptions, stated as such in the fields' sources: rubber on ice slides at a few hundredths near
// the melting point and at several tenths in deep cold, so its spread is wide; rubber on dry soil and on
// rock grips at the order of 0.5–1. On ice and rock the contact's stiffness and damping are those of the
// foot's rubber pad, in series with a ground far stiffer than it (a stiff, undamped contact on rock made the
// feet ring); on sediment the ground's own, softer. A measured law is the next step.
//
// The realisation. Each parameter varies as its class's median times exp(s z(x)), z a stationary Gaussian
// field of unit variance with a correlation length of 1 m (physis analysis, StationaryField) read at the
// point's place in the Earth's fixed frame, its own stream of the world's seed (synthetic).

#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <geodesy/figure/ellipsoid.hpp>
#include <geophysics/structure/crust.hpp>
#include <physis/analysis/stationary_field.hpp>
#include <physis/fields/field.hpp>
#include <physis/fields/sampling.hpp>

namespace geophysics::ground {

enum class GroundClass { ice, sediment, rock, sea };
inline const char* class_name(GroundClass c) { return c == GroundClass::ice ? "ice" : c == GroundClass::sediment ? "sediment" : c == GroundClass::rock ? "rock" : "sea"; }

struct Presumed {
  double median = 1, log_sd = 0;
  double at(double z) const { return median * std::exp(log_sd * z); }
};

struct GroundLaw {
  Presumed friction, load, damping;  // load: N at 1 mm
  double exponent = 1.5;
  std::string basis;
};

inline GroundLaw law_of(GroundClass c) {
  switch (c) {
    case GroundClass::ice: return {{0.15, 0.45}, {30, 0.4}, {1.5, 0.3}, 1.5, "assumed: rubber on ice, a few hundredths near the melting point to several tenths in deep cold"};
    case GroundClass::sediment: return {{0.60, 0.25}, {15, 0.6}, {2.0, 0.3}, 1.5, "assumed: a rubber foot on dry soil, sand or gravel"};
    case GroundClass::rock: return {{0.80, 0.15}, {40, 0.4}, {1.5, 0.3}, 1.5, "assumed: a rubber foot on dry rock"};
    case GroundClass::sea: break;
  }
  throw std::invalid_argument("geophysics::ground: the sea has no ground's law");
}

inline GroundClass class_of(const structure::CrustColumn& c) {
  if (c.water_m > 0) return GroundClass::sea;
  if (c.ice_m > 0) return GroundClass::ice;
  if (c.sediments_m > 0) return GroundClass::sediment;
  return GroundClass::rock;
}

// The parameters of the contact at a place.
struct GroundAt {
  GroundClass cls = GroundClass::rock;
  double friction = 0, load_1mm_n = 0, damping_s_m = 0, exponent = 1.5;
  double stiffness() const { return load_1mm_n / std::pow(1e-3, exponent); }
};

class EarthGround {
 public:
  using Point = physis::analysis::StationaryField<3>::Point;
  static constexpr double CORRELATION_M = 1.0;
  EarthGround(std::shared_ptr<const structure::Crust1> crust, std::uint64_t seed)
      : crust_(std::move(crust)), z_{field(1, seed), field(2, seed), field(3, seed)} {
    if (!crust_) throw std::invalid_argument("geophysics::ground::EarthGround: CRUST1.0");
  }
  GroundClass class_at(double lat_deg, double lon_deg) const { return class_of(crust_->column(lat_deg, lon_deg)); }
  // What is presumed of the ground at a place: its class's medians.
  GroundAt presumed(double lat_deg, double lon_deg) const {
    const auto c = class_at(lat_deg, lon_deg);
    if (c == GroundClass::sea) return {c, 0, 0, 0, 1.5};
    const auto l = law_of(c);
    return {c, l.friction.median, l.load.median, l.damping.median, l.exponent};
  }
  // The realisation at a point of the Earth's fixed frame (m).
  GroundAt at(const Point& x) const {
    static const geodesy::figure::Ellipsoid wgs84 = geodesy::figure::Ellipsoid::wgs84();
    const auto g = wgs84.geodetic({x[0], x[1], x[2]});
    const double d = 180 / 3.141592653589793;
    const auto c = class_at(g.latitude * d, g.longitude * d);
    if (c == GroundClass::sea) return {c, 0, 0, 0, 1.5};
    const auto l = law_of(c);
    return {c, l.friction.at(z_[0](x)), l.load.at(z_[1](x)), l.damping.at(z_[2](x)), l.exponent};
  }

 private:
  static physis::analysis::StationaryField<3> field(std::uint32_t stream, std::uint64_t seed) { return {1.0, CORRELATION_M, seed, 0x6EA0u + stream}; }
  std::shared_ptr<const structure::Crust1> crust_;
  std::array<physis::analysis::StationaryField<3>, 3> z_;
};

// The ground as fields of the world: a parameter presumed (model) or realised (synthetic), on the land's
// ground in the Earth's fixed frame; at sea refused.
class GroundField {
 public:
  using value_type = double;
  using point_type = physis::fields::Vec3;
  static constexpr std::uint32_t extra_capabilities = 0;
  enum class Quantity { friction, load, damping };
  GroundField(std::shared_ptr<const EarthGround> ground, Quantity q, bool presumed) : ground_(std::move(ground)), q_(q), presumed_(presumed) {
    namespace pf = physis::fields;
    static constexpr const char* names[] = {"friction", "load", "damping"};
    static constexpr const char* titles[] = {"The friction of a foot on the Earth's ground", "The load that presses a foot 1 mm into the Earth's ground",
                                             "The damping of a foot's contact with the Earth's ground"};
    d_.name = std::string("earth.ground.") + names[int(q)] + (presumed ? ".presumed" : "");
    d_.title = std::string(titles[int(q)]) + (presumed ? ", presumed" : "");
    d_.frame = "ITRS";
    d_.time = "none";
    d_.base = {pf::BaseKind::surface, 2, "the land's ground: points within 20 km of the WGS 84 ellipsoid where CRUST1.0 has no water, read at their place"};
    const physis::units::Dimension dims[] = {physis::units::dim::one, physis::units::dim::force, physis::units::dim::time / physis::units::dim::length};
    d_.fibre = {pf::FibreKind::scalar, dims[int(q)], 0, "", {}};
    d_.provenance = presumed ? pf::Provenance::model : pf::Provenance::synthetic;
    d_.source = presumed ? "presumed by the ground's class from CRUST1.0 (ice, sediment or rock): each class's median, an assumption"
                         : "a realisation of the presumed law of the place's class, its median × exp(s z), z stationary over 1 m";
    d_.representation = presumed ? pf::Representation{pf::RepresentationKind::uniform, "the class's median, by CRUST1.0's 1° cells"}
                                 : pf::Representation{pf::RepresentationKind::analytic, "sums of 256 random waves, the world's seed"};
  }
  const physis::fields::Descriptor& descriptor() const { return d_; }
  physis::fields::Evaluation<double> at(const physis::fields::Query<point_type>& q) const {
    namespace pf = physis::fields;
    static const geodesy::figure::Ellipsoid wgs84 = geodesy::figure::Ellipsoid::wgs84();
    const auto g = wgs84.geodetic(q.x);
    if (!(std::fabs(g.height) <= 2e4)) return pf::refused<double>(pf::Quality::outside, "not on the Earth's ground (more than 20 km from the ellipsoid)");
    const double d = 180 / 3.141592653589793;
    const auto a = presumed_ ? ground_->presumed(g.latitude * d, g.longitude * d) : ground_->at({q.x[0], q.x[1], q.x[2]});
    if (a.cls == GroundClass::sea) return pf::refused<double>(pf::Quality::outside, "not ground: CRUST1.0 has water at this place, the sea");
    const double v = q_ == Quantity::friction ? a.friction : q_ == Quantity::load ? a.load_1mm_n : a.damping_s_m;
    static constexpr const char* notes[] = {"ice, by CRUST1.0", "sediment, by CRUST1.0", "rock, by CRUST1.0"};
    return pf::valid(v, notes[int(a.cls)]);
  }
  template <physis::fields::Sampling G>
  physis::fields::Samples<double> sample(const G& g, double t) const {
    const auto q = physis::fields::queries_of(g, t);
    std::vector<physis::fields::Evaluation<double>> a(q.size());
    for (std::size_t i = 0; i < q.size(); i++) a[i] = at(q[i]);
    return physis::fields::samples_of(std::span<const physis::fields::Evaluation<double>>(a));
  }

 private:
  std::shared_ptr<const EarthGround> ground_;
  Quantity q_;
  bool presumed_;
  physis::fields::Descriptor d_;
};

}  // namespace geophysics::ground
