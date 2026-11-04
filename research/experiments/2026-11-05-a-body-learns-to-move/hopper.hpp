// The rabbit of the robotics domain (robotics/legged) moved by its controller with gains given, by semi-implicit
// Euler, as the domain's Simulation moves it.
//
// Why this loop exists. robotics::legged::Simulation builds its Controller with the default Gains: its
// SimulationConfig carries the method, the step, the ground and the presumption, but no gains, so a search cannot
// tune the controller through it. This loop is the Simulation's own for the semi-implicit method, step for step,
// from the domain's public parts (Rabbit, RabbitForces, potential_force, Controller, GroundBelief) and physis's
// semi_implicit_euler: the controller sampled every 2 ms with its torques held, the feet read as sensors and the
// belief of the ground conditioned on them (the same noise stream), the motors' electrical energy counted over each
// held period. With the default gains it reproduces the Simulation's state bit for bit (the experiment's `check`
// mode). What it adds is only what the experiment measures: falls, apexes, a window of speed and energy, frames.

#pragma once

#include <robotics/legged/controller.hpp>
#include <robotics/legged/ground.hpp>
#include <robotics/legged/rabbit.hpp>

#include <calx/random/stream.hpp>
#include <physis/mechanics/contact.hpp>
#include <physis/mechanics/multibody.hpp>

#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace experiment {

namespace rl = robotics::legged;
namespace pm = physis::mechanics;

// A frame for drawing: the points of the body in its plane (x, z).
struct Pose {
  double t = 0;
  // torso corners (4), hip, knee, ankle, toe, heel, shoulder, elbow, paw, centre of mass
  std::array<std::array<double, 2>, 13> p{};
  rl::Phase phase = rl::Phase::standing;
};

class Hopper {
 public:
  Hopper(const rl::Rabbit& rabbit, const rl::Gains& gains, double gravity = 1.62, bool feel = true, double step = 5e-4, double period = 2e-3)
      : r_(rabbit), g_(gravity), feel_(feel), h_(step), period_(period), tree_(pm::Base::floating, r_.links, gravity), forces_{&tree_, &r_, pm::GroundLine{}, &torque_, nullptr},
        ctl_(r_, gravity, gains), belief_(rl::GroundPresumption{}), noise_(1, 0xF007u, 0u, 0u) {
    r_.ground_law.validate();
    ctl_.attach(tree_);
    q_ = ctl_.standing_state(0.0, 0.0);
    v_.fill(0);
  }

  void command(const rl::Command& c) { ctl_.command(c); }
  const rl::Controller& controller() const { return ctl_; }
  const rl::Coordinates& q() const { return q_; }
  const rl::Coordinates& v() const { return v_; }
  double time() const { return t_; }
  double electrical_energy() const { return electrical_; }
  double copper_energy() const { return copper_; }
  double mass() const { return r_.mass(); }
  double gravity() const { return g_; }
  bool fallen() const { return fallen_; }
  double fall_time() const { return fall_t_; }
  double max_pitch() const { return max_pitch_; }
  void reset_max_pitch() { max_pitch_ = 0; }
  const std::vector<std::pair<double, double>>& apexes() const { return apexes_; }  // (time, apex) at each landing
  int landings() const { return int(apexes_.size()); }
  void keep_poses(double every) { pose_every_ = every, next_pose_ = t_; }
  const std::vector<Pose>& poses() const { return poses_; }

  void advance(double seconds) {
    const double end = t_ + seconds;
    while (t_ < end - 0.5 * h_) {
      if (t_ >= next_control_ - 1e-12) control();
      step();
      if (pose_every_ > 0 && t_ >= next_pose_ - 1e-12) pose(), next_pose_ += pose_every_;
    }
  }

  bool finite() const {
    for (std::size_t k = 0; k < rl::DOF; k++)
      if (!std::isfinite(q_[k]) || !std::isfinite(v_[k])) return false;
    return true;
  }

 private:
  void control() {
    rl::Observation o;
    o.t = t_, o.q = q_, o.v = v_;
    const auto cs = forces_.contacts(q_, v_);
    for (std::size_t i = 0; i < 3; i++) o.foot_force[i] = cs[i].normal;
    if (feel_) {
      sense(cs);
      const auto mu = belief_.friction.at(cs[0].place[0]);
      o.friction = mu.median * std::exp(-mu.log_sd);
    } else {
      // Without its feet's readings, what it was told: the presumption, less one deviation.
      o.friction = belief_.presumption.friction.median * std::exp(-belief_.presumption.friction.log_sd);
    }
    const rl::Phase before = ctl_.phase();
    torque_ = ctl_.step(o);
    next_control_ += period_;
    for (std::size_t j = 0; j < rl::JOINTS; j++) {
      const auto& m = r_.actuators[j];
      const double p = m.electrical_power(torque_[j], v_[rl::HIP + j]);
      electrical_ += std::max(p, 0.0) * period_;
      const double i = torque_[j] / m.torque_constant;
      copper_ += i * i * m.resistance * period_;
    }
    // What the experiment watches: the torso touching the ground (its four corners and its head are feet 3–6), the
    // pitch, and each landing's apex.
    for (std::size_t i = 3; i < cs.size(); i++)
      if (cs[i].touching && !fallen_) fallen_ = true, fall_t_ = t_;
    max_pitch_ = std::max(max_pitch_, std::fabs(q_[rl::PITCH]));
    if (before == rl::Phase::flight && ctl_.phase() == rl::Phase::compression) apexes_.push_back({t_, ctl_.last_apex()});
  }

  // The feet as sensors, as the domain's Simulation reads them (legged/simulation.hpp, sense).
  void sense(const std::vector<rl::FootContact>& cs) {
    for (std::size_t i = 0; i < 3; i++) {
      const auto& f = cs[i];
      const auto law = forces_.law_at(f.place[0]);
      if (f.touching && f.normal > 0.2 && std::fabs(f.slip) > 5 * law.slip_speed) {
        const double mu = std::fabs(f.tangential) / f.normal * (1 + 0.02 * noise_.normal());
        belief_.friction.observe(f.place[0], mu, 0.02);
      }
      if (f.touching && pen_[i][1] > pen_[i][0] && pen_[i][1] > f.penetration && fn_[i] > 0.5) {
        const double load = fn_[i] * std::pow(1e-3 / (pen_[i][1] * (1 + 0.10 * noise_.normal())), law.exponent);
        belief_.load.observe(x_[i], load, 0.20);
      }
      pen_[i][0] = pen_[i][1], pen_[i][1] = f.touching ? f.penetration : 0.0, fn_[i] = f.normal, x_[i] = f.place[0];
    }
  }

  void step() {
    const auto Q = [&](double t, const rl::Coordinates& q, const rl::Coordinates& v) {
      auto f = forces_(t, q, v);
      const auto p = rl::potential_force(r_, q);
      for (std::size_t k = 0; k < rl::DOF; k++) f[k] += p[k];
      return f;
    };
    pm::semi_implicit_euler(tree_, Q, t_, q_, v_, h_);
  }

  void pose() {
    Pose p;
    p.t = t_, p.phase = ctl_.phase();
    const auto st = tree_.kinematics(q_);
    const std::array<std::array<double, 2>, 4> corners = {{{-0.15, -0.05}, {0.15, -0.05}, {0.15, 0.05}, {-0.15, 0.05}}};
    for (int k = 0; k < 4; k++) p.p[k] = tree_.point(st, rl::TORSO, corners[k]);
    p.p[4] = tree_.point(st, rl::TORSO, r_.hip);
    p.p[5] = tree_.point(st, rl::THIGH, {0.10, 0.0});
    p.p[6] = tree_.point(st, rl::SHANK, {0.12, 0.0});
    p.p[7] = tree_.point(st, rl::FOOT, {0.10, 0.0});
    p.p[8] = tree_.point(st, rl::FOOT, {-0.01, 0.0});
    p.p[9] = tree_.point(st, rl::TORSO, r_.shoulder);
    p.p[10] = tree_.point(st, rl::ARM, {0.07, 0.0});
    p.p[11] = tree_.point(st, rl::FOREARM, {0.08, 0.0});
    p.p[12] = tree_.centre_of_mass(q_);
    poses_.push_back(p);
  }

  rl::Rabbit r_;
  double g_;
  bool feel_;
  double h_, period_;
  pm::PlanarTree<rl::DOF> tree_;
  std::array<double, rl::JOINTS> torque_{};
  rl::RabbitForces forces_;
  rl::Controller ctl_;
  rl::GroundBelief belief_;
  calx::random::Stream noise_;
  double t_ = 0, next_control_ = 0;
  rl::Coordinates q_{}, v_{};
  std::array<std::array<double, 2>, 3> pen_{};
  std::array<double, 3> fn_{}, x_{};
  double electrical_ = 0, copper_ = 0;
  bool fallen_ = false;
  double fall_t_ = 0, max_pitch_ = 0;
  std::vector<std::pair<double, double>> apexes_;
  double pose_every_ = 0, next_pose_ = 0;
  std::vector<Pose> poses_;
};

}  // namespace experiment
