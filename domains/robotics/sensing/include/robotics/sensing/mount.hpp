// Sensors on a vehicle: where a mount puts a sensor and which way it faces, and the noise of its readings
// (theory/vehicles.md §6).
//
// A mount is the sensor's place on the vehicle (its lever arm, in the vehicle's axes from the vehicle's
// reference point, the centre of its footprint on the ground) and its orientation (the rotation from the
// sensor's axes to the vehicle's). With the vehicle's pose A and place q, the sensor is at q + A ℓ and its
// axes are A O in the sphere's frame; a vector v of the world reads as (A O)ᵀ v in the sensor's axes. A mast
// is a lever arm along the vehicle's up.
//
// White noise is drawn by its index: sample k of a channel comes from calx's counter-based stream at
// (seed, channel, k), so that it does not depend on what else was drawn or in which order, and a run is the
// same on every machine (calx random, Salmon et al. 2011).

#pragma once

#include <cstdint>

#include <calx/random/stream.hpp>
#include <robotics/vehicles/surface.hpp>

namespace robotics::sensing {

using robotics::vehicles::Rotation;
using robotics::vehicles::Sphere;
using robotics::vehicles::SurfacePose;
using robotics::vehicles::Vec3;

struct Mount {
  Vec3 lever_arm{};      // m, in the vehicle's axes (forward, left, up)
  Rotation orientation;  // the sensor's axes → the vehicle's
};

// A sensor's place and axes in the sphere's frame.
struct SensorPose {
  Vec3 place{};
  Rotation axes;  // the sensor's axes → the sphere's frame
  Vec3 to_sensor(const Vec3& v) const { return axes.inverse().apply(v); }
  Vec3 from_sensor(const Vec3& v) const { return axes.apply(v); }
};

inline SensorPose sensor_pose(const Vec3& place, const SurfacePose& pose, const Mount& m) {
  using namespace physis::geometry::operators;
  return {place + pose.attitude().apply(m.lever_arm), pose.attitude() * m.orientation};
}
inline SensorPose sensor_pose(const SurfacePose& pose, const Sphere& s, const Mount& m) { return sensor_pose(pose.place(s), pose, m); }

// Independent normal errors of a standard deviation, by index.
struct WhiteNoise {
  double sd = 0;
  std::uint64_t seed = 0;
  std::uint32_t channel = 0;

  double sample(std::uint64_t k) const {
    calx::random::Stream z(seed, channel, static_cast<std::uint32_t>(k), static_cast<std::uint32_t>(k >> 32));
    return sd * z.normal();
  }
  // Three independent samples, for a vector sensor's axes.
  Vec3 vector(std::uint64_t k) const {
    calx::random::Stream z(seed, channel, static_cast<std::uint32_t>(k), static_cast<std::uint32_t>(k >> 32));
    const double a = z.normal(), b = z.normal(), c = z.normal();
    return {sd * a, sd * b, sd * c};
  }
};

}  // namespace robotics::sensing
