// The ground under a robot in its plane: the truth its feet press on, and its belief of it from what its
// feet feel (ADR 0059; theory/legged.md §8).
//
// The truth. Along the robot's line the ground's contact law varies: a function x ↦ (μ, k, n, c), the world's
// (in Cosmopoiesis, a realisation of the presumed law of the place's ground). Evaluating a field of the world
// at every foot, step and stage is costly; `TabulatedGround` reads it once a centimetre, as the robot goes,
// and interpolates linearly. The law varies on a metre's scale, so the table's error is a part in 10⁴.
//
// The belief. The robot is given a presumption of its ground (a virtual instrument reading the world's model
// of the place: a median and a spread for each parameter, and a correlation length), and learns from its
// feet. Its feet are its sensors:
//  - a foot that slides (|v_t| above five slip speeds) feels the friction itself: Coulomb's law gives
//    |f_t| = μ f_n while sliding, so |f_t|/f_n observes μ, through the foot's force sensor (2 % noise);
//  - a foot at its deepest (δ̇ changing sign) feels the stiffness: the damping term vanishes there, and
//    f_n (1 mm/δ)ⁿ observes the load at 1 mm. Its sinkage δ is known only from the leg's kinematics against
//    the ground's expected height, to 10 %, so the load to 15 %; and the deepest instant is caught at the
//    controller's samples, 2 ms apart: 20 % in all.
// The whitened innovations' mean square is the belief's own test of its presumption: near 1 on a ground
// drawn from the presumed law, far above it where the ground is not (a patch of ice on regolith).
// A foot that holds says only that μ exceeds |f_t|/f_n; that bound is not used. Each parameter's logarithm
// is a Gaussian field along x with the presumption's median, spread and correlation length; the readings
// condition it (tyche spatial, one observation at a time). Readings within 1 cm of the last taken of the
// same parameter are skipped, so that the observations do not repeat one another.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <span>
#include <stdexcept>
#include <vector>

#include <physis/mechanics/contact.hpp>
#include <tyche/spatial/conditioning.hpp>

namespace robotics::legged {

namespace pm = physis::mechanics;

// The truth along the line, read once a centimetre.
class TabulatedGround {
 public:
  using Law = std::function<pm::CompliantContact(double x)>;
  explicit TabulatedGround(Law law, double spacing = 0.01) : law_(std::move(law)), h_(spacing) {
    if (!law_ || !(spacing > 0)) throw std::invalid_argument("robotics::legged::TabulatedGround: a law and a positive spacing");
  }
  pm::CompliantContact at(double x) {
    const double u = std::floor(x / h_);
    const long i = long(u);
    const auto& a = node(i);
    const auto& b = node(i + 1);
    const double w = x / h_ - u;
    pm::CompliantContact c = a;
    c.stiffness = a.stiffness + w * (b.stiffness - a.stiffness);
    c.damping = a.damping + w * (b.damping - a.damping);
    c.friction = a.friction + w * (b.friction - a.friction);
    return c;
  }

 private:
  const pm::CompliantContact& node(long i) {
    if (nodes_.empty()) first_ = i, nodes_.push_back(law_(double(i) * h_));
    while (i < first_) first_--, nodes_.insert(nodes_.begin(), law_(double(first_) * h_));
    while (i >= first_ + long(nodes_.size())) nodes_.push_back(law_(double(first_ + long(nodes_.size())) * h_));
    return nodes_[std::size_t(i - first_)];
  }
  Law law_;
  double h_;
  long first_ = 0;
  std::vector<pm::CompliantContact> nodes_;
};

// A parameter presumed: its median, the spread of its logarithm, and the correlation length of its variation.
struct PresumedParameter {
  double median = 1, log_sd = 0.5;
};
struct GroundPresumption {
  PresumedParameter friction{0.8295, 0.185}, load{9.5, 0.5};  // load: N at 1 mm
  double exponent = 1.5;
  double correlation_m = 1.0;
};

// What the robot believes of one parameter along its line: a Gaussian field of its logarithm.
class ParameterBelief {
 public:
  ParameterBelief(PresumedParameter p, double correlation) : p_(p), l_(correlation) {}
  // A reading of the parameter at x, with its relative noise.
  bool observe(double x, double value, double relative_noise) {
    if (!(value > 0) || !std::isfinite(x)) return false;
    if (!xs_.empty() && std::fabs(x - xs_.back()) < 0.01) return false;
    if (xs_.size() >= 400) return false;
    std::vector<double> c(xs_.size());
    for (std::size_t i = 0; i < xs_.size(); i++) c[i] = cov(xs_[i], x);
    const double r = std::log1p(relative_noise);
    k_.observe(c, p_.log_sd * p_.log_sd, r * r, std::log(value / p_.median));
    xs_.push_back(x);
    return true;
  }
  // The belief at x: the median it implies and the spread of its logarithm.
  struct At { double median, log_sd; };
  At at(double x) const {
    std::vector<double> c(xs_.size());
    for (std::size_t i = 0; i < xs_.size(); i++) c[i] = cov(xs_[i], x);
    const auto post = k_.predict(c, 0.0, p_.log_sd * p_.log_sd);
    return {p_.median * std::exp(post.mean), std::sqrt(std::max(post.variance, 0.0))};
  }
  std::size_t readings() const { return xs_.size(); }
  const std::vector<double>& places() const { return xs_; }
  // The whitened innovations' mean square: 1 when the presumption is calibrated.
  double whitened_mean_square() const {
    const auto& z = k_.whitened();
    if (z.empty()) return 0;
    double s = 0;
    for (double v : z) s += v * v;
    return s / double(z.size());
  }
  const PresumedParameter& presumed() const { return p_; }

 private:
  double cov(double a, double b) const { return p_.log_sd * p_.log_sd * std::exp(-(a - b) * (a - b) / (2 * l_ * l_)); }
  PresumedParameter p_;
  double l_;
  tyche::spatial::GaussianConditioning k_;
  std::vector<double> xs_;
};

// The robot's belief of its ground: friction and load.
struct GroundBelief {
  explicit GroundBelief(const GroundPresumption& g = {}) : presumption(g), friction(g.friction, g.correlation_m), load(g.load, g.correlation_m) {}
  GroundPresumption presumption;
  ParameterBelief friction, load;
};

}  // namespace robotics::legged
