// The Moon's crustal field as fields (physis fields; the infrastructure's ADR 0053; theory/magnetism.md).
//
// LunarCrustalField is B of the Moon's crust on the axes of its mean-Earth frame, from its lowest ground
// to 10 000 km above its reference sphere, in one of two parts: the model (T2015_449, assimilated) or the
// part it leaves unresolved (drawn from its spectral law to a degree, with a seed: synthetic). Both are
// physis's potential fields, sharing their expansions with the CrustalField that holds them; on low
// ground below the sphere each is read at the sphere above the point (theory §3), with a note. Their sum
// is a composite (physis fields' sum), of the more derived provenance of the two: synthetic.

#pragma once

#include <cmath>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <physis/fields/potential.hpp>
#include <physis/fields/sampling.hpp>
#include <planetary_science/magnetism/crust.hpp>

namespace planetary_science::magnetism {

// The height above the reference sphere to which the field is offered: well inside the Moon's Hill
// sphere, and where the crust's field has long fallen below the solar wind's.
inline constexpr double LUNAR_CRUST_FIELD_TOP_M = 1e7;

class LunarCrustalField {
 public:
  using value_type = physis::fields::Vec3;
  using point_type = physis::fields::Vec3;
  static constexpr std::uint32_t extra_capabilities = static_cast<std::uint32_t>(physis::fields::Capability::sphere_grid);

  enum class Part { model, unresolved };

  explicit LunarCrustalField(const std::shared_ptr<const CrustalField>& f, Part part = Part::model, int degree = 720, std::uint64_t seed = 0)
      : part_(part), inner_(described(f, part, degree, seed), expansion(f, part, degree, seed)) {}

  const physis::fields::Descriptor& descriptor() const { return inner_.descriptor(); }
  Part part() const { return part_; }
  const physis::analysis::HarmonicPotential& expansion() const { return inner_.expansion(); }
  double radius() const { return inner_.expansion().radius; }

  physis::fields::Evaluation<value_type> at(const physis::fields::Query<point_type>& q) const {
    namespace pf = physis::fields;
    const double r = std::hypot(q.x[0], q.x[1], q.x[2]), a = radius();
    if (!(r >= a - LUNAR_LOWEST_GROUND_M)) return pf::refused<value_type>(pf::Quality::outside, "below the Moon's lowest ground");
    if (!(r <= a + LUNAR_CRUST_FIELD_TOP_M)) return pf::refused<value_type>(pf::Quality::outside, "more than 10 000 km above the Moon, beyond the field's base");
    if (r >= a) return inner_.at(q);
    // On the sphere itself, at the point's latitude and longitude.
    const double psi = std::atan2(q.x[2], std::hypot(q.x[0], q.x[1])), lambda = std::atan2(q.x[1], q.x[0]);
    return pf::valid(pf::PotentialField::cartesian(physis::analysis::field_at(expansion(), a, psi, lambda, inner_.degree()), psi, lambda),
                     "below the model's sphere: read at the sphere above, not continued downward");
  }

  // A sphere about the Moon's centre within the field's base above its reference sphere, row by row.
  physis::fields::Samples<value_type> sample(const physis::fields::SphereGrid& g, double t) const {
    if (g.centre == value_type{} && g.radius >= radius() && g.radius <= radius() + LUNAR_CRUST_FIELD_TOP_M) return inner_.sample(g, t);
    return generic(g, t);
  }
  template <physis::fields::Sampling G>
  physis::fields::Samples<value_type> sample(const G& g, double t) const { return generic(g, t); }

 private:
  static std::shared_ptr<const physis::analysis::HarmonicPotential> expansion(const std::shared_ptr<const CrustalField>& f, Part part, int degree,
                                                                              std::uint64_t seed) {
    if (!f) throw std::invalid_argument("planetary_science::magnetism::LunarCrustalField: a crustal field");
    // The model is held by f: an aliasing pointer keeps f alive without copying its coefficients.
    if (part == Part::model) return {f, &f->model()};
    return f->unresolved(degree, seed);
  }
  static physis::fields::Descriptor described(const std::shared_ptr<const CrustalField>& f, Part part, int degree, std::uint64_t seed) {
    namespace pf = physis::fields;
    if (!f) throw std::invalid_argument("planetary_science::magnetism::LunarCrustalField: a crustal field");
    pf::Descriptor d;
    d.frame = "MOON_ME";
    d.time = "none";
    d.base = {pf::BaseKind::region, 3, "the Moon from its lowest ground, 10 km below its sphere of 1737.4 km, to 10 000 km above it"};
    d.fibre = {pf::FibreKind::vector, physis::units::dim::magnetic_field, 0, "frame", {}};
    if (part == Part::model) {
      d.name = "moon.magnetic.crust";
      d.title = "The Moon's crustal field (T2015_449)";
      d.provenance = pf::Provenance::assimilated;
      d.source = f->source();
      d.representation.about = "spherical harmonics to degree " + std::to_string(f->model().degree) + " on a sphere of 1737.4 km";
    } else {
      d.name = "moon.magnetic.unresolved";
      d.title = "The scales of the Moon's crust its model does not resolve, drawn from their law";
      d.provenance = pf::Provenance::synthetic;
      d.source = "the law A nᵝ ρ²ⁿ fitted to T2015_449 on degrees 30–250 (planetary-science magnetism, theory §2)";
      d.representation.about = "spherical harmonics to degree " + std::to_string(degree) + ", a realisation of seed " + std::to_string(seed);
    }
    return d;
  }
  template <physis::fields::Sampling G>
  physis::fields::Samples<value_type> generic(const G& g, double t) const {
    const auto q = physis::fields::queries_of(g, t);
    std::vector<physis::fields::Evaluation<value_type>> a(q.size());
    for (std::size_t i = 0; i < q.size(); i++) a[i] = at(q[i]);
    return physis::fields::samples_of(std::span<const physis::fields::Evaluation<value_type>>(a));
  }

  Part part_;
  physis::fields::PotentialField inner_;
};

}  // namespace planetary_science::magnetism
