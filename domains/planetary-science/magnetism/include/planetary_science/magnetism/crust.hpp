// The crustal magnetic field of a body without a dynamo (ADR 0052; theory/magnetism.md §1–§3).
//
// The Moon's field at its surface is the remanence of its crust. Its potential is known to a degree
// from an assimilation of orbital magnetometry: T2015_449, the expansion by Wieczorek (2018) of
// Tsunakawa et al. (2015)'s surface vector mapping of Kaguya's and Lunar Prospector's data, to
// degree 449 on a sphere of 1737.4 km, in tesla, in the Moon's mean-Earth frame (the data library's
// moon-magnetic-t2015). A CrustalField holds such a model with the law of its spectrum, fitted on a
// band of degrees, and the realisations of what it leaves unresolved, each made once for its degree
// and seed and kept. The methods are physis's (analysis): this module names the body, the model,
// its band and its lowest ground.

#pragma once

#include <algorithm>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <utility>

#include <physis/analysis/potential_field.hpp>
#include <physis/analysis/random_field.hpp>

namespace planetary_science::magnetism {

namespace pa = physis::analysis;

// The band on which the Moon's law is fitted: above the largest scales, whose power is not the
// crust's statistics alone, and below the degrees that T2015's regularisation smooths.
inline constexpr int LUNAR_LAW_FROM = 30, LUNAR_LAW_TO = 250;
// The Moon's lowest ground below its reference sphere: the floor of South Pole–Aitken, some 9 km
// below 1737.4 km (Smith et al. 2010, LOLA); a magnetometer reads from there up.
inline constexpr double LUNAR_LOWEST_GROUND_M = 10e3;
// The key of the stream of the Moon's unresolved part ("MOON").
inline constexpr std::uint32_t LUNAR_STREAM = 0x4d4f4f4e;

class CrustalField {
 public:
  CrustalField(pa::HarmonicPotential model, int law_from, int law_to, std::string source, std::uint32_t stream)
      : model_(std::move(model)), law_(pa::fit_spectral_law(pa::lowes_spectrum(model_), law_from, law_to)),
        source_(std::move(source)), stream_(stream) {}

  const pa::HarmonicPotential& model() const { return model_; }
  const pa::SpectralLaw& law() const { return law_; }
  // ρ < 1: the law's power falls at high degree, so its unresolved part is bounded.
  bool law_bounded() const { return law_.ratio < 1; }
  const std::string& source() const { return source_; }
  // The effective depth of the sources (and of the model's regularisation together), m.
  double depth_m() const { return law_.depth(model_.radius); }

  // The part the model leaves unresolved to `degree`, from `seed`: the law's power less the model's,
  // and the law's own above the model's degree (physis unresolved_part). Synthetic.
  std::shared_ptr<const pa::HarmonicPotential> unresolved(int degree, std::uint64_t seed) const {
    std::lock_guard<std::mutex> lock(mutex_);
    const auto key = std::make_pair(degree, seed);
    auto it = drawn_.find(key);
    if (it == drawn_.end())
      it = drawn_.emplace(key, std::make_shared<const pa::HarmonicPotential>(pa::unresolved_part(model_, law_, degree, seed, stream_))).first;
    return it->second;
  }

 private:
  pa::HarmonicPotential model_;
  pa::SpectralLaw law_;
  std::string source_;
  std::uint32_t stream_;
  mutable std::mutex mutex_;
  mutable std::map<std::pair<int, std::uint64_t>, std::shared_ptr<const pa::HarmonicPotential>> drawn_;
};

// T2015_449 from the data library's file (SHTOOLS's format, its header in metres, its coefficients
// in tesla), with its law on degrees 30–250.
inline std::shared_ptr<CrustalField> read_lunar_t2015(const std::string& path) {
  return std::make_shared<CrustalField>(pa::read_shtools(path, 1.0), LUNAR_LAW_FROM, LUNAR_LAW_TO,
                                        "T2015_449 (Wieczorek 2018, from Tsunakawa et al. 2015)", LUNAR_STREAM);
}

// The crust's field at a mount, from its planetocentric radius, latitude and longitude.
struct CrustReading {
  pa::FieldValue field;   // north, east, down in the body's local horizon
  double radius_read_m;   // where the expansion was evaluated
  double below_sphere_m;  // how far the mount is below the reference sphere (0 above it)
};

// An expansion of the field outside its sources is not continued downward: there degree n grows as
// (a/r)ⁿ⁺², and the scales the model does not resolve with it. A mount on low ground, below the
// reference sphere, reads the model at the sphere above it.
inline CrustReading crust_at(const pa::HarmonicPotential& p, double r_m, double latitude, double longitude) {
  const double r = std::max(r_m, p.radius);
  return {pa::field_at(p, r, latitude, longitude), r, std::max(0.0, p.radius - r_m)};
}

}  // namespace planetary_science::magnetism
