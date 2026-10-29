// The World Magnetic Model as a field (physis fields; the infrastructure's ADR 0053).
//
// WmmField is B of WMM2025 on the Earth-fixed axes (the ITRS, as WGS 84 realises them), a vector field
// whose base is the region the model is made for: from 1 km below the WGS 84 ellipsoid to 850 km above
// it, during its lifetime of five years from its epoch (Chulliat et al. 2025, §1.1). Outside either a
// query is refused, with the reason. The model's time is a decimal year of UTC; the world's time
// coordinate is the caller's, who gives the conversion (in Cosmopoiesis, TT seconds from J2000.0 to UTC
// through the leap seconds, metrology's).
//
// The model is a potential field of degree 12, so ∇·B = 0 and ∇×B = 0 wherever it holds; the checks
// hold both by central differences.
//
// Its extent. The model is made for its base; its potential, the core's and the crust's field of long
// wavelength, holds outward as far as no other sources intervene. Continued outward (Extent::continued),
// the field is the core's field beyond 850 km, without the currents of the ionosphere and the
// magnetosphere: near the Earth what field lines follow, and only a sketch beyond some four Earth radii,
// where the solar wind shapes the field. Its answers above 850 km carry a note that says so; it is
// another field, named apart.

#pragma once

#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <geophysics/geomagnetism/wmm.hpp>
#include <physis/fields/field.hpp>

namespace geophysics::geomagnetism {

class WmmField {
 public:
  using value_type = physis::fields::Vec3;
  using point_type = physis::fields::Vec3;

  enum class Extent { model, continued };

  WmmField(std::shared_ptr<const Coefficients> c, std::function<double(double t)> decimal_year_of, Extent extent = Extent::model)
      : c_(std::move(c)), year_(std::move(decimal_year_of)), extent_(extent) {
    if (!c_) throw std::invalid_argument("geophysics::geomagnetism::WmmField: coefficients");
    namespace pf = physis::fields;
    const bool model = extent == Extent::model;
    d_.name = model ? "earth.magnetic" : "earth.magnetic.core";
    d_.title = model ? "The geomagnetic field (" + c_->model + ")" : "The core's field (" + c_->model + "), continued outward";
    d_.frame = "ITRS";
    d_.time = "the caller's, converted to a decimal year of UTC";
    d_.base = {pf::BaseKind::region, 3,
               model ? "the Earth from 1 km below the WGS 84 ellipsoid to 850 km above it, during the model's five years"
                     : "from 1 km below the WGS 84 ellipsoid outward, during the model's five years; beyond 850 km the core's field alone"};
    d_.fibre = {pf::FibreKind::vector, physis::units::dim::magnetic_field, 0, "frame", {}};
    d_.representation = {pf::RepresentationKind::spectral, "Gauss coefficients to degree 12 with their secular variation, Schmidt semi-normalised"};
    d_.provenance = pf::Provenance::model;
    d_.source = "Chulliat et al. 2025, The US/UK World Magnetic Model for 2025–2030, NOAA NCEI and BGS";
  }

  const physis::fields::Descriptor& descriptor() const { return d_; }
  const Coefficients& coefficients() const { return *c_; }
  Extent extent() const { return extent_; }

  physis::fields::Evaluation<value_type> at(const physis::fields::Query<point_type>& q) const {
    namespace pf = physis::fields;
    const Geodetic place = wgs84().geodetic(q.x);
    if (!(place.height >= LOWEST_HEIGHT_M)) return pf::refused<value_type>(pf::Quality::outside, "more than 1 km below the WGS 84 ellipsoid, beneath WMM's base");
    const bool above = !(place.height <= HIGHEST_HEIGHT_M);
    if (above && extent_ == Extent::model) return pf::refused<value_type>(pf::Quality::outside, "more than 850 km above the WGS 84 ellipsoid, above WMM's base");
    const double year = year_(q.t);
    if (!(year >= c_->valid_from_year() && year <= c_->valid_until_year())) return pf::refused<value_type>(pf::Quality::outside, "outside the model's five years");
    return pf::valid(field(*c_, place, year).itrs_T, above ? "above 850 km: the core's field continued outward, without the magnetosphere's currents" : "");
  }

 private:
  std::shared_ptr<const Coefficients> c_;
  std::function<double(double)> year_;
  Extent extent_;
  physis::fields::Descriptor d_;
};

}  // namespace geophysics::geomagnetism
