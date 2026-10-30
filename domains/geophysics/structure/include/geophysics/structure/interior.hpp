// The Earth's interior beneath a place: CRUST1.0 from the top of its column to its Moho, PREM below,
// as one column, as a profile from the surface to the centre, and as fields (physis fields; the
// infrastructure's ADR 0053) of density and of the speeds of P and S waves.
//
// The composition is the one seismology commonly makes of a crustal model over a radial one: depth is
// measured below mean sea level, along the normal of the ellipsoid, and PREM's spherical depth is
// taken as that depth (no correction for the Earth's ellipticity, an approximation of the order of its
// flattening, at most 21 km near the centre). Where the local Moho is shallower than PREM's (24.4 km,
// under the oceans), the mantle beneath it is PREM's lid continued upward, its polynomial evaluated
// at the depth; PREM's own crust and ocean, global averages, are not used. Above the column's top (in
// the air) a query is refused, as it is beyond the centre.

#pragma once

#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
#include <string>
#include <string_view>
#include <vector>

#include <geophysics/structure/crust.hpp>
#include <geophysics/structure/prem.hpp>
#include <physis/fields/field.hpp>

namespace geophysics::structure {

// A place: geodetic latitude and longitude (rad, WGS 84) and its height above mean sea level (m),
// negative below.
struct InteriorPlace {
  double latitude = 0, longitude = 0, height_m = 0;
};

struct InteriorValues {
  double density_kg_m3 = 0, vp_m_s = 0, vs_m_s = 0;
  std::string_view layer;  // CRUST1.0's layer, or PREM's shell
  bool crust1 = false;     // CRUST1.0's (true) or PREM's (false)
};

// PREM's lid, index 9 among PREM_SHELLS: the mantle continued upward to a Moho shallower than PREM's.
inline constexpr int PREM_LID = 9;

// PREM at a depth below sea level beneath a Moho: its own shell (at a boundary the one below, or with
// `above` the one above), or above 24.4 km its lid continued upward.
inline InteriorValues prem_beneath(double depth_m, bool above = false) {
  if (depth_m < 24.4e3) {
    const auto& s = PREM_SHELLS[PREM_LID];
    const double x = (PREM_RADIUS_M - depth_m) / PREM_RADIUS_M;
    return {prem_detail::poly(s.rho, x) * 1e3, prem_detail::poly(s.vp, x) * 1e3, prem_detail::poly(s.vs, x) * 1e3, "lid (PREM, continued upward)", false};
  }
  const PremValues p = prem_at_depth(std::min(depth_m, PREM_RADIUS_M), above);
  return {p.density_kg_m3, p.vp_m_s, p.vs_m_s, p.shell, false};
}

// The composite at a height above sea level beneath a column; nothing (layer empty) above its top or
// beyond the centre.
inline InteriorValues interior_at(const CrustColumn& c, double height_m) {
  if (height_m > c.surface_m || -height_m > PREM_RADIUS_M) return {};
  const int k = c.layer_at(height_m);
  if (k >= 0 && k < int(c.layers.size())) {
    const auto& l = c.layers[std::size_t(k)];
    return {l.density_kg_m3, l.vp_m_s, l.vs_m_s, l.name, true};
  }
  return prem_beneath(-height_m);
}

// One sample of a profile: a depth below sea level (m) and the values there. At a boundary there are
// two samples of the same depth, the one above first.
struct ProfileSample {
  double depth_m = 0;
  double density_kg_m3 = 0, vp_m_s = 0, vs_m_s = 0;
  std::string_view layer;
  bool crust1 = false;
};

// The profile beneath a column, from its top to the centre: CRUST1.0's layers, each constant, then
// PREM's shells, each sampled at its ends and at most every `step_m` between.
inline std::vector<ProfileSample> interior_profile(const CrustColumn& c, double step_m = 50e3) {
  std::vector<ProfileSample> out;
  for (const auto& l : c.layers) {
    out.push_back({-l.top_m, l.density_kg_m3, l.vp_m_s, l.vs_m_s, l.name, true});
    out.push_back({-l.bottom_m, l.density_kg_m3, l.vp_m_s, l.vs_m_s, l.name, true});
  }
  // PREM from the Moho down: the boundaries of its mantle and core below the Moho (not of its crust,
  // which CRUST1.0 replaces), each shell sampled between them.
  std::vector<double> cuts{c.moho_depth_m};
  for (int k = 1; k <= PREM_LID; k++) {
    const double d = PREM_RADIUS_M - PREM_SHELLS[std::size_t(k)].inner_km * 1e3;
    if (d > c.moho_depth_m + 1) cuts.push_back(d);
  }
  cuts.push_back(PREM_RADIUS_M);
  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
  for (std::size_t i = 0; i + 1 < cuts.size(); i++) {
    const double top = cuts[i], bottom = cuts[i + 1];
    const int n = std::max(1, int(std::ceil((bottom - top) / step_m)));
    for (int j = 0; j <= n; j++) {
      const double d = top + (bottom - top) * j / n;
      // Within the shell: its top taken from below the boundary, its bottom from above the next.
      const InteriorValues v = prem_beneath(d, j == n);
      out.push_back({d, v.density_kg_m3, v.vp_m_s, v.vs_m_s, v.layer, false});
    }
  }
  return out;
}

enum class InteriorQuantity { density, vp, vs };

// One quantity of the interior as a scalar field on the places beneath the surface.
class InteriorField {
 public:
  using value_type = double;
  using point_type = InteriorPlace;

  InteriorField(std::shared_ptr<const Crust1> crust, InteriorQuantity q) : crust_(std::move(crust)), q_(q) {
    namespace pf = physis::fields;
    namespace dim = physis::units::dim;
    static constexpr std::string_view names[] = {"density", "vp", "vs"};
    static constexpr std::string_view titles[] = {"Density", "Speed of P waves", "Speed of S waves"};
    const auto i = static_cast<std::size_t>(q);
    d_.name = "interior." + std::string(names[i]);
    d_.title = std::string(titles[i]) + " of the Earth's interior (CRUST1.0 over PREM)";
    d_.frame = "geodetic latitude and longitude (WGS 84), height above mean sea level";
    d_.time = "none";
    d_.base = {pf::BaseKind::region, 3, "the Earth beneath its ground and its seas, from the top of CRUST1.0's column to the centre"};
    d_.fibre = {pf::FibreKind::scalar, q == InteriorQuantity::density ? dim::mass / dim::volume : dim::velocity, 0, "", {}};
    d_.representation = {pf::RepresentationKind::composite, "CRUST1.0's cells of 1° (constant in each layer) to the Moho; PREM's polynomials in the radius below"};
    d_.provenance = pf::Provenance::model;
    d_.source = "CRUST1.0 (Laske et al. 2013) over PREM (Dziewonski and Anderson 1981, isotropic, 1 s)";
  }

  const physis::fields::Descriptor& descriptor() const { return d_; }

  physis::fields::Evaluation<double> at(const physis::fields::Query<InteriorPlace>& q) const {
    namespace pf = physis::fields;
    constexpr double deg = 180 / std::numbers::pi;
    if (-q.x.height_m > PREM_RADIUS_M) return pf::refused<double>(pf::Quality::outside, "beyond the centre of the Earth");
    const CrustColumn c = crust_->column(q.x.latitude * deg, q.x.longitude * deg);
    if (q.x.height_m > c.surface_m) return pf::refused<double>(pf::Quality::outside, "above the ground or the sea: in the air");
    const InteriorValues v = interior_at(c, q.x.height_m);
    const double x = q_ == InteriorQuantity::density ? v.density_kg_m3 : q_ == InteriorQuantity::vp ? v.vp_m_s : v.vs_m_s;
    return pf::valid(x, v.crust1 ? std::string_view{} : std::string_view{"PREM, below CRUST1.0's Moho"});
  }

 private:
  std::shared_ptr<const Crust1> crust_;
  InteriorQuantity q_;
  physis::fields::Descriptor d_;
};

}  // namespace geophysics::structure
