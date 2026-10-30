// The ocean's water as a body in it meets it: its density and its current near the surface, presumed where
// nothing has measured them, and a realisation of that presumption, the water a simulated world holds
// (ADR 0059; theory/water.md).
//
// The density. Seawater near the surface weighs 1020 to 1030 kg/m³ over the world's oceans, by its
// temperature and salinity (Talley, Pickard, Emery and Swift 2011, Descriptive Physical Oceanography, ch. 3).
// Presumed log-normal with a median of 1025 kg/m³ and a spread of its logarithm of 0.0025, so that 95 % lies
// within 1020–1030 kg/m³. Its equation of state (TEOS-10) from the water's temperature and salinity is the
// domain's planned `seawater`, which will replace the presumption where those are known.
//
// The current. Each horizontal component near the surface is presumed normal with mean 0 and a deviation of
// 0.1 m/s, varying over a kilometre. Assumed: the order of the open ocean's slow currents away from its
// boundary currents; HYCOM's analysis of the place, which the data library fetches on demand, is the next
// source.
//
// The realisation. Each quantity varies as its median times exp(s z(x)), or as s z(x) for the current's
// components, z a stationary Gaussian field of unit variance (physis analysis, StationaryField) read at the
// point's place in the Earth's fixed frame, its own stream of the world's seed (synthetic).
//
// The base: the upper ocean, within 200 m of mean sea level. Where the sea is, and where the land, is not
// this module's to tell (the geophysics domain's CRUST1.0 knows the water's depth): it answers wherever it is
// asked within that layer.

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

namespace oceanography::water {

inline constexpr double EARTH_MEAN_RADIUS_M = 6371.0e3;
inline constexpr double UPPER_OCEAN_M = 200.0;

struct WaterLaw {
  double density_median = 1025, density_log_sd = 0.0025;  // kg/m³
  double density_correlation_m = 1.0e4;                   // assumed: water masses vary over kilometres
  double current_sd = 0.1;                                // m/s, each horizontal component
  double current_correlation_m = 1.0e3;                   // assumed
  std::string density_basis = "surface seawater, 1020–1030 kg/m³ over the world's oceans (Talley et al. 2011, ch. 3): 95 % within that";
  std::string current_basis = "assumed: each horizontal component normal, mean 0 and 0.1 m/s, varying over 1 km";
};

// The water near a point: its density, and its current on the local east and north.
struct WaterAt {
  double density = 0;
  double east = 0, north = 0;  // m/s
};

class UpperOcean {
 public:
  using Point = physis::analysis::StationaryField<3>::Point;
  explicit UpperOcean(WaterLaw law = {}, std::uint64_t seed = 0)
      : law_(std::move(law)), rho_(1.0, law_.density_correlation_m, seed, 0x0CE1u), ue_(1.0, law_.current_correlation_m, seed, 0x0CE2u),
        un_(1.0, law_.current_correlation_m, seed, 0x0CE3u) {}
  const WaterLaw& law() const { return law_; }
  WaterAt presumed() const { return {law_.density_median, 0.0, 0.0}; }
  // The realisation at a point of the Earth's fixed frame (m).
  WaterAt at(const Point& x) const {
    return {law_.density_median * std::exp(law_.density_log_sd * rho_(x)), law_.current_sd * ue_(x), law_.current_sd * un_(x)};
  }
  // The current as a vector of the Earth's fixed frame at x: east and north turned by the local frame of the
  // geocentric direction (the difference from the geodetic normal is 0.2° at most, nothing for a current).
  std::array<double, 3> current_fixed(const Point& x, const WaterAt& w) const {
    const double r = std::hypot(x[0], x[1], x[2]), lon = std::atan2(x[1], x[0]), lat = std::asin(x[2] / r);
    const std::array<double, 3> E{-std::sin(lon), std::cos(lon), 0.0}, N{-std::sin(lat) * std::cos(lon), -std::sin(lat) * std::sin(lon), std::cos(lat)};
    return {w.east * E[0] + w.north * N[0], w.east * E[1] + w.north * N[1], w.east * E[2] + w.north * N[2]};
  }

 private:
  WaterLaw law_;
  physis::analysis::StationaryField<3> rho_, ue_, un_;
};

namespace detail {
inline bool in_upper_ocean(const physis::fields::Vec3& x) {
  // The ellipsoid's radius at the point's latitude, to a metre's error: within 200 m of it.
  const double a = 6378137.0, b = 6356752.314245, r = std::hypot(x[0], x[1], x[2]);
  const double s = x[2] / r, c = std::sqrt(1 - s * s);
  const double R = a * b / std::sqrt(b * b * c * c + a * a * s * s);
  return std::fabs(r - R) <= UPPER_OCEAN_M;
}
}  // namespace detail

// The water's density as a field of the world: presumed (model) or realised (synthetic).
class DensityField {
 public:
  using value_type = double;
  using point_type = physis::fields::Vec3;
  static constexpr std::uint32_t extra_capabilities = 0;
  DensityField(std::shared_ptr<const UpperOcean> ocean, bool presumed) : ocean_(std::move(ocean)), presumed_(presumed) {
    namespace pf = physis::fields;
    d_.name = presumed ? "earth.ocean.density.presumed" : "earth.ocean.density";
    d_.title = presumed ? "The density of the upper ocean's water, presumed" : "The density of the upper ocean's water";
    d_.frame = "ITRS";
    d_.time = "none";
    d_.base = {pf::BaseKind::region, 3, "the upper ocean: within 200 m of the WGS 84 ellipsoid (where the sea is, is not checked)"};
    d_.fibre = {pf::FibreKind::scalar, physis::units::dim::mass / physis::units::dim::volume, 0, "", {}};
    d_.provenance = presumed ? pf::Provenance::model : pf::Provenance::synthetic;
    d_.source = presumed ? "presumed: median 1025 kg/m³, " + ocean_->law().density_basis : "a realisation of the presumed law: median × exp(0.0025 z), z stationary over 10 km";
    d_.representation = presumed ? pf::Representation{pf::RepresentationKind::uniform, "the presumption's median"} : pf::Representation{pf::RepresentationKind::analytic, "a sum of 256 random waves, the world's seed"};
  }
  const physis::fields::Descriptor& descriptor() const { return d_; }
  physis::fields::Evaluation<double> at(const physis::fields::Query<point_type>& q) const {
    namespace pf = physis::fields;
    if (!detail::in_upper_ocean(q.x)) return pf::refused<double>(pf::Quality::outside, "not in the upper ocean (more than 200 m from mean sea level)");
    if (presumed_) return pf::valid(ocean_->law().density_median, "presumed: the same wherever nothing has measured it");
    return pf::valid(ocean_->at({q.x[0], q.x[1], q.x[2]}).density);
  }
  template <physis::fields::Sampling G>
  physis::fields::Samples<double> sample(const G& g, double t) const {
    const auto q = physis::fields::queries_of(g, t);
    std::vector<physis::fields::Evaluation<double>> a(q.size());
    for (std::size_t i = 0; i < q.size(); i++) a[i] = at(q[i]);
    return physis::fields::samples_of(std::span<const physis::fields::Evaluation<double>>(a));
  }

 private:
  std::shared_ptr<const UpperOcean> ocean_;
  bool presumed_;
  physis::fields::Descriptor d_;
};

// The water's current as a field of the world, a vector of the Earth's fixed frame: presumed (none) or realised.
class CurrentField {
 public:
  using value_type = physis::fields::Vec3;
  using point_type = physis::fields::Vec3;
  static constexpr std::uint32_t extra_capabilities = 0;
  CurrentField(std::shared_ptr<const UpperOcean> ocean, bool presumed) : ocean_(std::move(ocean)), presumed_(presumed) {
    namespace pf = physis::fields;
    d_.name = presumed ? "earth.ocean.current.presumed" : "earth.ocean.current";
    d_.title = presumed ? "The upper ocean's current, presumed" : "The upper ocean's current";
    d_.frame = "ITRS";
    d_.time = "none";
    d_.base = {pf::BaseKind::region, 3, "the upper ocean: within 200 m of the WGS 84 ellipsoid (where the sea is, is not checked)"};
    d_.fibre = {pf::FibreKind::vector, physis::units::dim::velocity, 0, "frame", {}};
    d_.provenance = presumed ? pf::Provenance::model : pf::Provenance::synthetic;
    d_.source = presumed ? "presumed: none, " + ocean_->law().current_basis : "a realisation of the presumed law, each horizontal component 0.1 z m/s, z stationary over 1 km";
    d_.representation = presumed ? pf::Representation{pf::RepresentationKind::uniform, "zero"} : pf::Representation{pf::RepresentationKind::analytic, "sums of 256 random waves, the world's seed"};
  }
  const physis::fields::Descriptor& descriptor() const { return d_; }
  physis::fields::Evaluation<value_type> at(const physis::fields::Query<point_type>& q) const {
    namespace pf = physis::fields;
    if (!detail::in_upper_ocean(q.x)) return pf::refused<value_type>(pf::Quality::outside, "not in the upper ocean (more than 200 m from mean sea level)");
    if (presumed_) return pf::valid(value_type{0, 0, 0}, "presumed: none wherever nothing has measured it");
    const UpperOcean::Point x{q.x[0], q.x[1], q.x[2]};
    const auto c = ocean_->current_fixed(x, ocean_->at(x));
    return pf::valid(value_type{c[0], c[1], c[2]}, "horizontal: on the local east and north");
  }
  template <physis::fields::Sampling G>
  physis::fields::Samples<value_type> sample(const G& g, double t) const {
    const auto q = physis::fields::queries_of(g, t);
    std::vector<physis::fields::Evaluation<value_type>> a(q.size());
    for (std::size_t i = 0; i < q.size(); i++) a[i] = at(q[i]);
    return physis::fields::samples_of(std::span<const physis::fields::Evaluation<value_type>>(a));
  }

 private:
  std::shared_ptr<const UpperOcean> ocean_;
  bool presumed_;
  physis::fields::Descriptor d_;
};

}  // namespace oceanography::water
