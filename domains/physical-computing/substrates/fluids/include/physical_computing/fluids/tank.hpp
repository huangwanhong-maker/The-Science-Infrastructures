// A tank of water as a reservoir: Fernando and Sojakka's "pattern recognition in a bucket" (2003)
// (substrates/fluids/README.md, simulator 1).
//
// The surface of a shallow tank obeys the damped wave equation of physis's Basin with the long-wave speed
// c = √(g h). Motors at the edges (or within) are forces spread over small Gaussian footprints; the input of
// a motor is its forcing (m/s²), held over each step of the schedule. A camera sees the surface as a coarse
// grid of pixels, each the mean height under it, or thresholded (Fernando and Sojakka used thresholded images
// of the waves as the reservoir's state, read by a linear perceptron). The tank forgets its start at the
// damping rate γ: the distance between two runs under the same inputs falls as e^{−γt}, the echo state property
// of a linear substrate, for every departure of the surface but one. The mean level is the volume of water, which
// the tank keeps: the uniform mode is neutral, the one thing the tank never forgets. The forcing is therefore made
// volume-conserving, and the property holds for surfaces of equal volume.
//
// Limits (declared): long waves only, without dispersion (√(gh) at every wavelength, true for kh ≪ 1); no
// capillarity, no nonlinearity of the surface; the damping is one declared rate, the walls' and floor's boundary
// layers and the surface film folded into it. Checks: tests/check_tank.cpp.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

#include <physis/continuum/waves.hpp>

#include <physical_computing/fluids/fluid.hpp>
#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::fluids {

struct Motor {
  double x = 0, y = 0;  // m
  double radius = 0.01; // m, the footprint's standard deviation
};

struct WaveTankOptions {
  double length = 0.4, width = 0.3;  // m
  double depth = 0.005;              // m: c = √(g h) ≈ 0.22 m/s
  std::size_t nx = 64, ny = 48;
  double damping = 0.8;              // γ, 1/s, declared
  std::vector<Motor> motors{{0.08, 0.15, 0.012}, {0.32, 0.15, 0.012}};
  std::size_t camera_x = 16, camera_y = 12;
  double threshold = 0;              // m, for binary images
  Fluid fluid = Fluid::water();
};

class WaveTank final : public runtime::Backend {
 public:
  explicit WaveTank(const WaveTankOptions& o) : o_(o), basin_(basin_options(o)) {
    if (o.motors.empty() || o.camera_x == 0 || o.camera_y == 0 || o.camera_x > o.nx || o.camera_y > o.ny)
      throw std::invalid_argument("physical_computing::fluids::WaveTank: at least one motor, and a camera no finer than the grid");
    for (const auto& m : o.motors) {
      std::vector<double> w(o.nx * o.ny);
      double total = 0;
      for (std::size_t j = 0; j < o.ny; j++)
        for (std::size_t i = 0; i < o.nx; i++) {
          const double x = (i + 0.5) * basin_.dx() - m.x, y = (j + 0.5) * basin_.dy() - m.y;
          total += (w[j * o.nx + i] = std::exp(-(x * x + y * y) / (2 * m.radius * m.radius)));
        }
      // Normalised so that a unit input forces the footprint's mean at 1 m/s² over its area π r², then made
      // volume-conserving: a plunger displaces water and adds none, so the forcing's mean over the tank is zero.
      const double cell = basin_.dx() * basin_.dy(), norm = 3.14159265358979 * m.radius * m.radius / (total * cell);
      double mean = 0;
      for (double& v : w) v *= norm, mean += v;
      mean /= static_cast<double>(w.size());
      for (double& v : w) v -= mean;
      footprints_.push_back(std::move(w));
    }
  }

  static physis::continuum::Basin::Options basin_options(const WaveTankOptions& o) {
    physis::continuum::Basin::Options b;
    b.nx = o.nx, b.ny = o.ny, b.lx = o.length, b.ly = o.width;
    b.speed = std::sqrt(standard_gravity * o.depth), b.damping = o.damping;
    return b;
  }

  double wave_speed() const { return basin_.options().speed; }
  const physis::continuum::Basin& basin() const { return basin_; }
  physis::continuum::Basin& basin() { return basin_; }
  const WaveTankOptions& options() const { return o_; }

  runtime::TargetDescription describe() const override {
    runtime::TargetDescription d;
    d.family = "fluids", d.substrate = "wave tank";
    d.accepts = {tasks::Kind::filter, tasks::Kind::phenomenon};
    for (std::size_t k = 0; k < o_.motors.size(); k++) d.inputs.push_back({"motor " + std::to_string(k), "m/s²"});
    for (std::size_t q = 0; q < o_.camera_x * o_.camera_y; q++) d.outputs.push_back({"pixel " + std::to_string(q), "m"});
    d.time_scale = o_.length / wave_speed();
    d.power = 0;  // not declared: the motors' and camera's power depend on the build
    d.materials.materials = {{"water", "the medium"}, {"acrylic (PMMA)", "the tank"}, {"copper", "the motors' windings"}};
    d.materials.reuse = computer::Reuse::reusable;
    d.materials.end_of_life = computer::EndOfLife::recyclable;
    d.model = "the damped wave equation of a shallow surface, c = √(g h) (physis::continuum::Basin); Fernando and Sojakka 2003";
    d.limits = {"long waves without dispersion", "no capillarity or surface nonlinearity", "one declared damping rate"};
    return d;
  }
  void reset(std::uint64_t seed) override {
    seed_ = seed;
    basin_.set(std::vector<double>(o_.nx * o_.ny, 0.0));
  }
  // Holds the motors' forcing for `duration` s: whole steps of the basin, at least one.
  void apply(std::span<const double> inputs, double duration) override {
    if (inputs.size() != o_.motors.size()) throw std::invalid_argument("physical_computing::fluids::WaveTank::apply: one input per motor");
    std::vector<double> f(o_.nx * o_.ny, 0.0);
    for (std::size_t k = 0; k < inputs.size(); k++)
      for (std::size_t c = 0; c < f.size(); c++) f[c] += inputs[k] * footprints_[k][c];
    const auto steps = std::max<std::size_t>(1, static_cast<std::size_t>(std::llround(duration / basin_.step_size())));
    for (std::size_t s = 0; s < steps; s++) basin_.step(&f);
  }
  // The camera: the mean height under each pixel.
  std::vector<double> read() const override {
    std::vector<double> px(o_.camera_x * o_.camera_y, 0.0);
    std::vector<int> count(px.size(), 0);
    for (std::size_t j = 0; j < o_.ny; j++)
      for (std::size_t i = 0; i < o_.nx; i++) {
        const std::size_t q = (j * o_.camera_y / o_.ny) * o_.camera_x + (i * o_.camera_x / o_.nx);
        px[q] += basin_.at(i, j), count[q]++;
      }
    for (std::size_t q = 0; q < px.size(); q++) px[q] /= count[q];
    return px;
  }
  // The camera thresholded: 1 where the surface stands above the threshold.
  std::vector<double> binary_image() const {
    auto px = read();
    for (double& v : px) v = v > o_.threshold ? 1.0 : 0.0;
    return px;
  }
  double time() const override { return basin_.time(); }
  runtime::Provenance provenance() const override {
    auto p = simulator_provenance("wave tank", "physis::continuum::Basin, c = √(g h)", seed_);
    p.conditions = {{"depth (m)", o_.depth}, {"damping (1/s)", o_.damping}};
    return p;
  }

 private:
  WaveTankOptions o_;
  physis::continuum::Basin basin_;
  std::vector<std::vector<double>> footprints_;
  std::uint64_t seed_ = 0;
};

}  // namespace physical_computing::fluids
