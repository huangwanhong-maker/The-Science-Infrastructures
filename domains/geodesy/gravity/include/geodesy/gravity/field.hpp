// Normal gravity as a field (physis fields; the infrastructure's ADR 0053).
//
// NormalGravityField is γ = ∇U of a level ellipsoid on its Earth-fixed axes (normal.hpp): the gravity of
// the rotating Earth's reference figure, its attraction and the centrifugal acceleration together, a
// vector field over the ellipsoid's exterior and, with a note, over its analytic continuation down to
// 11 km below it, the deepest ocean floor (normal.hpp, Validity). It is static on the rotating axes.
// Being a gradient it has no curl, and outside the masses its divergence is that of the centrifugal
// part alone: ∇·γ = ∇²U = ∇²(½ ω² (x² + y²)) = 2ω², which the checks hold.

#pragma once

#include <utility>

#include <geodesy/gravity/normal.hpp>
#include <physis/fields/field.hpp>

namespace geodesy::gravity {

// The depth below the ellipsoid to which the field is offered: the deepest ocean floor (the Challenger
// Deep, some 10.9 km below sea level).
inline constexpr double NORMAL_GRAVITY_LOWEST_M = -11e3;

class NormalGravityField {
 public:
  using value_type = physis::fields::Vec3;
  using point_type = physis::fields::Vec3;

  explicit NormalGravityField(LevelEllipsoid e = LevelEllipsoid::wgs84()) : e_(std::move(e)) {
    namespace pf = physis::fields;
    d_.name = "earth.gravity.normal";
    d_.title = "Normal gravity of the level ellipsoid of " + e_.name();
    d_.frame = "ITRS";
    d_.time = "none";
    d_.base = {pf::BaseKind::region, 3, "the exterior of the level ellipsoid, and its continuation to 11 km below it"};
    d_.fibre = {pf::FibreKind::vector, physis::units::dim::acceleration, 0, "frame", {}};
    d_.representation = {pf::RepresentationKind::analytic, "closed form in ellipsoidal-harmonic coordinates (NGA eqs. 4-4 to 4-13)"};
    d_.provenance = pf::Provenance::model;
    d_.source = "NGA.STND.0036_1.0.0_WGS84 (2014), §4; Moritz 2000, Geodetic Reference System 1980";
  }

  const physis::fields::Descriptor& descriptor() const { return d_; }
  const LevelEllipsoid& ellipsoid() const { return e_; }

  physis::fields::Evaluation<value_type> at(const physis::fields::Query<point_type>& q) const {
    namespace pf = physis::fields;
    const Geodetic p = e_.figure().geodetic(q.x);
    if (!(p.height >= NORMAL_GRAVITY_LOWEST_M)) return pf::refused<value_type>(pf::Quality::outside, "more than 11 km below the ellipsoid, beneath the field's continuation");
    return pf::valid(e_.at(p).ecef_m_s2, p.height < 0 ? "below the ellipsoid: its analytic continuation, as free-air reduction takes it" : "");
  }

 private:
  LevelEllipsoid e_;
  physis::fields::Descriptor d_;
};

}  // namespace geodesy::gravity
