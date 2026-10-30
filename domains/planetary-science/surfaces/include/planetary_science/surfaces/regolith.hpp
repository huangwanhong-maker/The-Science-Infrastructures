// The Moon's ground as a robot's foot meets it: the regolith's contact law, presumed where nothing has
// measured it, and a realisation of that presumption, the ground a simulated world stands on (ADR 0059;
// theory/surfaces.md).
//
// The law. A foot presses into the ground by δ and slides on it at v_t; physis's compliant contact gives the
// normal force f_n = k δⁿ (1 + c δ̇) and the friction −μ f_n tanh(v_t/v_s). Three parameters name the ground:
//  - μ, the friction of a foot on regolith. Presumed tan φ, with φ the regolith's angle of internal friction,
//    30° to 50° in the returned samples and in situ (Carrier, Olhoeft and Mendell 1991, Lunar Sourcebook,
//    ch. 9): a log-normal law whose 95 % interval is tan 30° to tan 50°, its median their geometric mean,
//    √(tan 30° tan 50°) = 0.830, and the spread of its logarithm ln(tan 50°/tan 30°)/(2 × 1.96) = 0.185;
//  - the load that presses a foot 1 mm in, k (1 mm)ⁿ, with Hertz's exponent n = 3/2. Assumed: 9.5 N, the
//    rabbit's contact as built, with a spread of its logarithm of 0.5 (a factor of 2.7 at one deviation),
//    since no measurement of a small foot on lunar regolith is adopted here;
//  - c, the damping of the contact. Assumed: 1.5 s/m, a spread of 0.3.
// The load is given at 1 mm, not as k, whose unit N/m^{3/2} has a half-integer exponent of length.
//
// Its variation. Each parameter varies over the ground as its median times exp(s z(x)), s its spread and z a
// stationary Gaussian field of unit variance and a correlation length of 1 m (assumed: a robot's step), read
// at the point's place in the Moon's fixed frame (physis analysis, StationaryField). Independent fields for
// the three, each its own stream of the world's seed. The presumption is the median (provenance model); a
// realisation is the ground of a simulated world (synthetic), which a robot's feet read and the presumption
// does not know.
//
// Caveat. Regolith is not elastic: it compacts under a foot and keeps the print (the Apollo boot prints).
// The compliant law is its elastic part; its plastic sinkage is open (theory §4).

#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <physis/analysis/stationary_field.hpp>
#include <physis/fields/field.hpp>
#include <physis/fields/sampling.hpp>

namespace planetary_science::surfaces {

inline constexpr double MOON_RADIUS_M = 1737.4e3;
inline constexpr double GROUND_LAYER_M = 2.0e4;  // the ground lies within 20 km of the sphere (its height is not checked)

// A parameter presumed where nothing has measured it: its median, the standard deviation of its logarithm,
// and the basis of both.
struct Presumed {
  double median = 1, log_sd = 0;
  std::string basis;
  double at(double z) const { return median * std::exp(log_sd * z); }  // its value at a standard normal z
};

// The parameters of a contact law at a place.
struct GroundContact {
  double friction = 0;   // μ
  double load_1mm_n = 0; // k (1 mm)ⁿ, N
  double damping_s_m = 0;
  double exponent = 1.5;
  double stiffness() const { return load_1mm_n / std::pow(1e-3, exponent); }  // k, N/m^n
};

struct RegolithLaw {
  Presumed friction{0.829493, 0.184884, "tan φ, φ the regolith's angle of internal friction, 30°–50° (Lunar Sourcebook, ch. 9): 95 % within tan 30°–tan 50°, the median their geometric mean"};
  Presumed load{9.5, 0.5, "assumed: the rabbit's contact as built, Hertz's exponent; a factor of e^{±1} for one deviation"};
  Presumed damping{1.5, 0.3, "assumed: a restitution of about 0.9 at 0.1 m/s"};
  double exponent = 1.5;
  double correlation_m = 1.0;  // assumed: a robot's step
  GroundContact median() const { return {friction.median, load.median, damping.median, exponent}; }
};

// The Moon's ground: the presumed law, and its realisation for a world's seed.
class RegolithGround {
 public:
  using Point = physis::analysis::StationaryField<3>::Point;
  explicit RegolithGround(RegolithLaw law = {}, std::uint64_t seed = 0)
      : law_(std::move(law)), z_{field(1, seed), field(2, seed), field(3, seed)} {}

  const RegolithLaw& law() const { return law_; }
  GroundContact presumed() const { return law_.median(); }
  // The ground at a point of the Moon's fixed frame (m): the realisation.
  GroundContact at(const Point& x) const {
    return {law_.friction.at(z_[0](x)), law_.load.at(z_[1](x)), law_.damping.at(z_[2](x)), law_.exponent};
  }
  // The standard normal z of each parameter at x: (friction, load, damping).
  std::array<double, 3> z(const Point& x) const { return {z_[0](x), z_[1](x), z_[2](x)}; }

 private:
  physis::analysis::StationaryField<3> field(std::uint32_t stream, std::uint64_t seed) const {
    return physis::analysis::StationaryField<3>(1.0, law_.correlation_m, seed, 0x6E00u + stream);
  }
  RegolithLaw law_;
  std::array<physis::analysis::StationaryField<3>, 3> z_;
};

// The ground as fields of the world: a parameter presumed (model) or its realisation (synthetic), a scalar
// field on the Moon's ground in its fixed frame.
class RegolithField {
 public:
  using value_type = double;
  using point_type = physis::fields::Vec3;
  static constexpr std::uint32_t extra_capabilities = 0;
  enum class Quantity { friction, load, damping };
  enum class Kind { presumed, realised };

  RegolithField(std::shared_ptr<const RegolithGround> ground, Quantity q, Kind kind) : ground_(std::move(ground)), q_(q), kind_(kind) {
    if (!ground_) throw std::invalid_argument("planetary_science::surfaces::RegolithField: a ground");
    namespace pf = physis::fields;
    static constexpr const char* names[] = {"friction", "load", "damping"};
    static constexpr const char* titles[] = {"The friction of a foot on the Moon's regolith", "The load that presses a foot 1 mm into the Moon's regolith",
                                             "The damping of a foot's contact with the Moon's regolith"};
    d_.name = std::string("moon.ground.") + names[int(q)] + (kind == Kind::presumed ? ".presumed" : "");
    d_.title = std::string(titles[int(q)]) + (kind == Kind::presumed ? ", presumed" : "");
    d_.frame = "MOON_ME";
    d_.time = "none";
    d_.base = {pf::BaseKind::surface, 2, "the Moon's ground: points within 20 km of its sphere of 1737.4 km, read at their place (their height is not checked)"};
    const physis::units::Dimension dims[] = {physis::units::dim::one, physis::units::dim::force, physis::units::dim::time / physis::units::dim::length};
    d_.fibre = {pf::FibreKind::scalar, dims[int(q)], 0, "", {}};
    const Presumed& p = parameter();
    if (kind == Kind::presumed) {
      d_.provenance = pf::Provenance::model;
      d_.source = "presumed: median " + std::to_string(p.median).substr(0, 6) + ", " + p.basis;
      d_.representation = {pf::RepresentationKind::uniform, "the presumption's median, the same everywhere"};
    } else {
      d_.provenance = pf::Provenance::synthetic;
      d_.source = "a realisation of the presumed law: median × exp(" + std::to_string(p.log_sd).substr(0, 5) + " z), z stationary with a correlation length of " +
                  std::to_string(ground_->law().correlation_m).substr(0, 4) + " m (physis StationaryField)";
      d_.representation = {pf::RepresentationKind::analytic, "a sum of 256 random waves for each parameter, the world's seed"};
    }
  }

  const physis::fields::Descriptor& descriptor() const { return d_; }
  const Presumed& parameter() const {
    const auto& l = ground_->law();
    return q_ == Quantity::friction ? l.friction : q_ == Quantity::load ? l.load : l.damping;
  }

  physis::fields::Evaluation<double> at(const physis::fields::Query<point_type>& q) const {
    namespace pf = physis::fields;
    const double r = std::hypot(q.x[0], q.x[1], q.x[2]);
    if (!(std::fabs(r - MOON_RADIUS_M) <= GROUND_LAYER_M)) return pf::refused<double>(pf::Quality::outside, "not on the Moon's ground (more than 20 km from its sphere)");
    if (kind_ == Kind::presumed) return pf::valid(parameter().median, "presumed: the same wherever nothing has measured it");
    const auto c = ground_->at({q.x[0], q.x[1], q.x[2]});
    return pf::valid(q_ == Quantity::friction ? c.friction : q_ == Quantity::load ? c.load_1mm_n : c.damping_s_m);
  }

  template <physis::fields::Sampling G>
  physis::fields::Samples<double> sample(const G& g, double t) const {
    const auto q = physis::fields::queries_of(g, t);
    std::vector<physis::fields::Evaluation<double>> a(q.size());
    for (std::size_t i = 0; i < q.size(); i++) a[i] = at(q[i]);
    return physis::fields::samples_of(std::span<const physis::fields::Evaluation<double>>(a));
  }

 private:
  std::shared_ptr<const RegolithGround> ground_;
  Quantity q_;
  Kind kind_;
  physis::fields::Descriptor d_;
};

}  // namespace planetary_science::surfaces
