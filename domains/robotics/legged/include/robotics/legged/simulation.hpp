// The rabbit simulated: its body, its controller and the ground, advanced by one of three methods of one
// motion (ADR 0058; theory/legged.md §4).
//
// The methods. One Lagrangian (physis multibody.hpp: the tree's kinetic energy and gravity, with the tendon
// and the limits as RabbitPotential) and one set of forces (RabbitForces: the motors' held torques, the
// gears' friction, the limits' damping, the ground), advanced by
//  - the variational integrator of the discrete action (physis variational.hpp, the midpoint rule, the
//    forces by Lagrange–d'Alembert): symplectic, of order 2, a momentum map kept exactly when its symmetry
//    is; implicit, each step solved by Newton's method on jets: the accurate and costly end;
//  - semi-implicit Euler on Newton–Euler's equations: of order 1, cheap, the real-time end;
//  - the classical Runge–Kutta method on them: of order 4, four evaluations a step.
// They are the same motion: the Euler–Lagrange equations of the Lagrangian are Newton–Euler's, to rounding
// (physis's checks). What differs is how the step is taken, its cost and what it keeps.
//
// The controller runs at its own period (2 ms), its torques held over the steps between (zero-order hold),
// as a robot's are; its observation is the state at its sample.
//
// The account of energy. The robot's energy E = T + V (gravity, tendon, limits) changes by the work of the
// motors, the gears' friction, the limits' damping and the ground: E(t) − E(0) = ∫ Σ Q_i v_i dt, the power of
// RabbitForces' generalized forces. Its residual, measured with the power by the trapezoidal rule over each
// step, is the method's own error of energy, shown with the motion.
//
// The angular momentum about the centre of mass, in flight, where only gravity acts: constant in the motion,
// and measured by how much each flight changes it.
//
// The ground (ground.hpp). Where the configuration gives the ground's law along the line (a world's), each
// foot presses on it where it touches; else on the rabbit's own law, the same everywhere. The feet are
// sensors, sampled with the controller: a sliding foot reads the friction, a foot at its deepest the load at
// 1 mm, each with its sensor's noise (a stream of its own), and the robot's belief of its ground, which
// starts from the presumption it was given, is conditioned on them. The controller asks of the ground no
// more friction than it believes is there.

#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <physis/mechanics/contact.hpp>
#include <physis/mechanics/multibody.hpp>
#include <physis/mechanics/variational.hpp>
#include <robotics/legged/controller.hpp>
#include <robotics/legged/rabbit.hpp>

namespace robotics::legged {

enum class Method { variational, semi_implicit, runge_kutta };
inline const char* method_name(Method m) { return m == Method::variational ? "variational" : m == Method::semi_implicit ? "semi_implicit" : "runge_kutta"; }
inline Method method_of(const std::string& s) {
  if (s == "variational") return Method::variational;
  if (s == "runge_kutta" || s == "rk4") return Method::runge_kutta;
  if (s == "semi_implicit" || s == "euler") return Method::semi_implicit;
  throw std::invalid_argument("robotics::legged: the method is variational, semi_implicit or runge_kutta");
}

struct SimulationConfig {
  Method method = Method::semi_implicit;
  double gravity = 1.62;          // m/s², the body's at the site
  double step = 5e-4;             // s
  double control_period = 2e-3;   // s
  double frame_period = 1.0 / 60; // s, between frames kept for drawing
  pm::GroundLine ground{};
  // The ground's law along the line (a world's), or none: the rabbit's own law everywhere.
  std::function<pm::CompliantContact(double x)> ground_at;
  // What the robot is told of its ground before it touches it (a virtual instrument's reading of the
  // world's model), and the seed of its feet's sensor noise.
  GroundPresumption presumption{};
  std::uint64_t sensor_seed = 1;
};

// A reading of a foot: where, which parameter (0 friction, 1 load) and its value.
struct FootReading {
  double t, x;
  int parameter;
  double value;
};

// A frame of the motion, for drawing and for the readings.
struct Frame {
  double t = 0;
  Coordinates q{};
  std::array<double, 3> normal{};
  std::array<double, JOINTS> torque{};
  std::array<double, 2> com{};
  Phase phase = Phase::standing;
};

class Simulation {
  using L = pm::TreeLagrangian<DOF, RabbitPotential>;
  using VI = pm::VariationalIntegrator<DOF, L, RabbitForces>;

 public:
  Simulation(const Rabbit& rabbit, SimulationConfig c)
      : r_(rabbit), c_(c), tree_(pm::Base::floating, r_.links, c.gravity), forces_{&tree_, &r_, c.ground, &torque_}, ctl_(r_, c.gravity),
        belief_(c.presumption), noise_(c.sensor_seed, 0xF007u, 0u, 0u) {
    if (!(c.step > 0) || !(c.control_period >= c.step)) throw std::invalid_argument("robotics::legged::Simulation: a positive step within the controller's period");
    r_.ground_law.validate();
    if (c_.ground_at) {
      table_ = std::make_unique<TabulatedGround>(c_.ground_at);
      forces_.table = table_.get();
    }
    ctl_.attach(tree_);
    q_ = ctl_.standing_state(0.0, c.ground.z0);
    // Settle the start onto the ground: the standing height puts the feet on it, pressed in by their share
    // of the weight.
    v_.fill(0);
    if (c_.method == Method::variational) start_variational();
    e0_ = energy();
    record(true);
  }

  void command(const Command& cmd) { ctl_.command(cmd); }
  const Controller& controller() const { return ctl_; }
  const pm::PlanarTree<DOF>& tree() const { return tree_; }
  const Rabbit& rabbit() const { return r_; }
  const SimulationConfig& config() const { return c_; }
  double time() const { return t_; }
  const Coordinates& q() const { return q_; }
  const Coordinates& v() const { return v_; }
  const std::vector<Frame>& frames() const { return frames_; }
  std::size_t frames_dropped() const { return dropped_; }
  // The variational steps that fell back to semi-implicit Euler (their implicit equation unsolved).
  std::size_t fallbacks() const { return fallbacks_; }

  // Advance by `seconds` of the robot's time.
  void advance(double seconds) {
    const auto wall0 = std::chrono::steady_clock::now();
    const double end = t_ + seconds;
    while (t_ < end - 0.5 * c_.step) {
      if (t_ >= next_control_ - 1e-12) control();
      const double P0 = power();
      step();
      const double P1 = power();
      work_ += 0.5 * c_.step * (P0 + P1);
      account();
      record(false);
    }
    wall_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
    simulated_ += seconds;
  }

  // The energy, its account and what the motors drew.
  double energy() const { return tree_.energy(q_, v_) + RabbitPotential{&r_}(q_); }
  double energy_residual() const { return energy() - e0_ - work_; }
  double motor_work() const { return motor_work_; }
  double electrical_energy() const { return electrical_; }
  double peak_power() const { return peak_power_; }
  double distance() const { return q_[X] - x_start_; }
  // Electrical energy over m g d: the cost of transport.
  double cost_of_transport() const { return std::fabs(distance()) > 0.05 ? electrical_ / (r_.mass() * c_.gravity * std::fabs(distance())) : 0.0; }
  double flight_momentum_change() const { return flight_dL_; }
  double max_penetration() const { return max_pen_; }
  // The machine's seconds per second of the robot's time: the method's cost here.
  double cost() const { return simulated_ > 0 ? wall_ / simulated_ : 0.0; }
  const std::array<double, JOINTS>& torque() const { return torque_; }
  std::vector<FootContact> contacts() const { return forces_.contacts(q_, v_); }
  // The ground: its law at x (the truth the feet press on), the robot's belief of it, its feet's readings.
  pm::CompliantContact ground_at(double x) const { return forces_.law_at(x); }
  const GroundBelief& belief() const { return belief_; }
  const std::vector<FootReading>& readings() const { return readings_; }
  std::size_t readings_dropped() const { return readings_dropped_; }
  void trim_readings(std::size_t keep_from) {
    const std::size_t k = std::min(keep_from, readings_.size());
    readings_.erase(readings_.begin(), readings_.begin() + std::ptrdiff_t(k));
    readings_dropped_ += k;
  }

  // Keep only the frames after index i (the page has read them).
  void trim_frames(std::size_t keep_from) {
    if (keep_from == 0 || frames_.empty()) return;
    const std::size_t k = std::min(keep_from, frames_.size());
    frames_.erase(frames_.begin(), frames_.begin() + std::ptrdiff_t(k));
    dropped_ += k;
  }

 private:
  void start_variational() {
    vi_ = std::make_unique<VI>(pm::DiscreteLagrangian<DOF, L>(L{&tree_, RabbitPotential{&r_}}, c_.step), forces_, 1e-11, 60);
    vs_ = vi_->start(t_, q_, v_);
  }

  void control() {
    Observation o;
    o.t = t_, o.q = q_, o.v = v_;
    const auto cs = forces_.contacts(q_, v_);
    for (std::size_t i = 0; i < 3; i++) o.foot_force[i] = cs[i].normal;
    sense(cs);
    // The friction the controller may count on at its toe: its belief, less one deviation.
    const auto mu = belief_.friction.at(cs[0].place[0]);
    o.friction = mu.median * std::exp(-mu.log_sd);
    torque_ = ctl_.step(o);
    next_control_ += c_.control_period;
    // The electrical energy of the motors, over the held period.
    for (std::size_t j = 0; j < JOINTS; j++) {
      const double p = r_.actuators[j].electrical_power(torque_[j], v_[HIP + j]);
      electrical_ += std::max(p, 0.0) * c_.control_period;
      peak_power_ = std::max(peak_power_, p);
    }
  }

  // The feet as sensors: a sliding foot reads the friction (|f_t|/f_n, 2 % noise); a foot at its deepest
  // reads the load at 1 mm (f_n (1 mm/δ)ⁿ, its sinkage known to 10 %, so the load to 15 %, and taken at the
  // sample nearest the deepest instant: 20 % in all, as the belief is told).
  void sense(const std::vector<FootContact>& cs) {
    for (std::size_t i = 0; i < 3; i++) {
      const auto& f = cs[i];
      const auto law = forces_.law_at(f.place[0]);
      if (f.touching && f.normal > 0.2 && std::fabs(f.slip) > 5 * law.slip_speed) {
        const double mu = std::fabs(f.tangential) / f.normal * (1 + 0.02 * noise_.normal());
        if (belief_.friction.observe(f.place[0], mu, 0.02)) record_reading(f.place[0], 0, mu);
      }
      // the deepest: the sinkage rose to the last sample and falls now
      if (f.touching && pen_[i][1] > pen_[i][0] && pen_[i][1] > f.penetration && fn_[i] > 0.5) {
        const double load = fn_[i] * std::pow(1e-3 / (pen_[i][1] * (1 + 0.10 * noise_.normal())), law.exponent);
        if (belief_.load.observe(x_[i], load, 0.20)) record_reading(x_[i], 1, load);
      }
      pen_[i][0] = pen_[i][1], pen_[i][1] = f.touching ? f.penetration : 0.0, fn_[i] = f.normal, x_[i] = f.place[0];
    }
  }
  void record_reading(double x, int parameter, double value) {
    readings_.push_back({t_, x, parameter, value});
    if (readings_.size() > 4000) trim_readings(2000);
  }

  void step() {
    const auto Q = [&](double t, const Coordinates& q, const Coordinates& v) {
      auto f = forces_(t, q, v);
      const auto p = potential_force(r_, q);
      for (std::size_t k = 0; k < DOF; k++) f[k] += p[k];
      return f;
    };
    if (c_.method == Method::variational) {
      // The step's implicit equation is nonsmooth where a foot touches or a joint meets its limit, and Newton's
      // method may not converge there; that step is then taken by semi-implicit Euler, counted, and the
      // variational state restarted from it.
      try {
        vi_->step(vs_);
        // The velocity from the discrete momentum by the Legendre transform, p = M(q) v (the Lagrangian is
        // quadratic in v): the motion's own velocity at the step, not the difference quotient of the step.
        t_ = vs_.t, q_ = vs_.q, v_ = calx::linalg::solve<DOF>(tree_.mass_matrix(q_), vs_.p);
      } catch (const std::runtime_error&) {
        fallbacks_++;
        pm::semi_implicit_euler(tree_, Q, t_, q_, v_, c_.step);
        vs_ = vi_->start(t_, q_, v_);
      }
      return;
    }
    if (c_.method == Method::semi_implicit) pm::semi_implicit_euler(tree_, Q, t_, q_, v_, c_.step);
    else pm::runge_kutta4(tree_, Q, t_, q_, v_, c_.step);
  }

  // The power of the generalized forces, Σ Q_i v_i, and the motors' share.
  double power() {
    const auto Q = forces_(t_, q_, v_);
    double P = 0;
    for (std::size_t k = 0; k < DOF; k++) P += Q[k] * v_[k];
    return P;
  }

  void account() {
    double pm = 0;
    for (std::size_t j = 0; j < JOINTS; j++) pm += torque_[j] * v_[HIP + j];
    motor_work_ += pm * c_.step;
    const auto cs = forces_.contacts(q_, v_);
    bool any = false;
    for (const auto& f : cs) any = any || f.touching, max_pen_ = std::max(max_pen_, f.penetration);
    // The angular momentum about the centre of mass while nothing touches.
    if (!any) {
      const double L = tree_.angular_momentum(q_, v_, tree_.centre_of_mass(q_));
      if (!in_flight_) in_flight_ = true, flight_L0_ = L, flight_scale_ = std::max(std::fabs(L), 1e-3);
      else flight_dL_ = std::max(flight_dL_, std::fabs(L - flight_L0_) / flight_scale_);
    } else {
      in_flight_ = false;
    }
  }

  void record(bool force) {
    if (!force && t_ < next_frame_ - 1e-12) return;
    next_frame_ = t_ + c_.frame_period;
    Frame f;
    f.t = t_, f.q = q_, f.torque = torque_, f.com = tree_.centre_of_mass(q_), f.phase = ctl_.phase();
    const auto cs = forces_.contacts(q_, v_);
    for (std::size_t i = 0; i < 3; i++) f.normal[i] = cs[i].normal;
    frames_.push_back(f);
    if (frames_.size() > 20000) frames_.erase(frames_.begin(), frames_.begin() + 10000), dropped_ += 10000;
  }

  Rabbit r_;
  SimulationConfig c_;
  pm::PlanarTree<DOF> tree_;
  std::array<double, JOINTS> torque_{};
  RabbitForces forces_;
  Controller ctl_;
  std::unique_ptr<VI> vi_;
  pm::PhaseState<DOF> vs_;
  double t_ = 0, next_control_ = 0, next_frame_ = 0;
  Coordinates q_{}, v_{};
  double e0_ = 0, work_ = 0, motor_work_ = 0, electrical_ = 0, peak_power_ = 0, x_start_ = 0;
  double flight_L0_ = 0, flight_scale_ = 1, flight_dL_ = 0, max_pen_ = 0;
  bool in_flight_ = false;
  double wall_ = 0, simulated_ = 0;
  std::vector<Frame> frames_;
  std::size_t dropped_ = 0, fallbacks_ = 0;
  std::unique_ptr<TabulatedGround> table_;
  GroundBelief belief_;
  calx::random::Stream noise_;
  std::array<std::array<double, 2>, 3> pen_{};
  std::array<double, 3> fn_{}, x_{};
  std::vector<FootReading> readings_;
  std::size_t readings_dropped_ = 0;
};

}  // namespace robotics::legged
