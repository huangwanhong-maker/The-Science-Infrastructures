// The U.S. Standard Atmosphere 1976 as fields (physis fields; the infrastructure's ADR 0053).
//
// The standard atmosphere is a function of height alone: its base is a line, geometric height above
// mean sea level, from −5 km′ of geopotential height to 86 km (NOAA, NASA and USAF 1976, p. 12 and
// Table 4). StandardAtmosphereField is one of its quantities on that line, a scalar field with its unit;
// beyond the line a query is refused. A world carries it into space as a pullback along each point's
// height above its geoid (physis fields' pullback), which needs a figure and a geoid that are the
// world's, not the atmosphere's.

#pragma once

#include <array>
#include <string>
#include <string_view>

#include <meteorology/atmosphere/standard.hpp>
#include <physis/fields/field.hpp>

namespace meteorology::atmosphere {

enum class Quantity { temperature, pressure, density, speed_of_sound };

class StandardAtmosphereField {
 public:
  using value_type = double;
  using point_type = double;

  explicit StandardAtmosphereField(Quantity q) : q_(q) {
    namespace pf = physis::fields;
    namespace dim = physis::units::dim;
    struct Named {
      std::string_view name, title;
      physis::units::Dimension dimension;
    };
    constexpr std::array<Named, 4> named{{{"temperature", "Kinetic temperature", dim::temperature},
                                          {"pressure", "Pressure", dim::pressure},
                                          {"density", "Density", dim::mass / dim::volume},
                                          {"speed_of_sound", "Speed of sound", dim::velocity}}};
    const Named& n = named[static_cast<std::size_t>(q)];
    d_.name = "atmosphere." + std::string(n.name);
    d_.title = std::string(n.title) + " of the U.S. Standard Atmosphere 1976";
    d_.frame = "geometric height above mean sea level";
    d_.time = "none";
    d_.base = {pf::BaseKind::curve, 1, "geometric height from −5 km′ geopotential to 86 km above mean sea level"};
    d_.fibre = {pf::FibreKind::scalar, n.dimension, 0, "", {}};
    d_.representation = {pf::RepresentationKind::analytic, "the 1976 standard's seven layers of linear molecular temperature in geopotential height"};
    d_.provenance = pf::Provenance::model;
    d_.source = "NOAA, NASA and USAF 1976, U.S. Standard Atmosphere, 1976, NOAA-S/T 76-1562";
  }

  const physis::fields::Descriptor& descriptor() const { return d_; }

  physis::fields::Evaluation<double> at(const physis::fields::Query<double>& q) const {
    namespace pf = physis::fields;
    const auto s = at_geometric_height(q.x);
    if (!s) return pf::refused<double>(pf::Quality::outside, "outside the 1976 standard's heights, −5 km′ to 86 km");
    switch (q_) {
      case Quantity::temperature: return pf::valid(s->temperature_K);
      case Quantity::pressure: return pf::valid(s->pressure_Pa);
      case Quantity::density: return pf::valid(s->density_kg_m3);
      case Quantity::speed_of_sound: return pf::valid(s->speed_of_sound_m_s);
    }
    return pf::refused<double>(pf::Quality::invalid, "an unknown quantity");
  }

 private:
  Quantity q_;
  physis::fields::Descriptor d_;
};

}  // namespace meteorology::atmosphere
