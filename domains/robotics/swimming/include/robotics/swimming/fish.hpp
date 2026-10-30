// A fish: a swimming robot of 1.85 kg in its horizontal plane, its body, its controller and its simulation by
// three methods of one motion (ADR 0059; theory/swimming.md).
//
// The body. A rigid head (the base, holding the battery and the electronics) and a chain of four links behind
// it: two body segments, the peduncle and the caudal fin, joined by four driven joints. Its configuration
// space is SE(2) × T⁴: q = (x, y, ψ; q₁ … q₄), the head's centre and heading in the plane of swimming and the
// joints' angles. It is neutrally buoyant, so gravity does not act in that plane. The tree's plane (x, z) is
// the horizontal plane here: x forward along the heading at rest, the second axis to the fish's left.
//
// The water (physis fluid.hpp): each link a slender part of its length and its height (the body's depth;
// the fin's span), the water's kinetic energy (added mass ρπh²/4 a unit length, across each strip) in the
// Lagrangian, the drag of the flow across each strip and the friction along it as forces, in the water's
// density and current where it swims. The drag coefficients are assumed: C_n = 1.2 for a plate across the
// flow, C_f = 0.03 along the body for its skin friction and its form together.
//
// The motors (actuation/motor.hpp): geared, with a torque–speed line, armature, friction and copper loss.
//
// The controller: a travelling wave (a central pattern generator, Ijspeert 2008, Neural Networks 21,
// 642–653): each joint follows q*_i = A_i sin(φ − iλ) + b_i, its amplitude growing toward the tail
// (carangiform swimming), the phase φ advancing at the frequency f, λ the lag between joints; the speed is
// servoed by f, the heading by a bias b on the front joints. The wave runs at 200 Hz and gives each joint's
// servo its angle and rate, held between; the servo, as a hobby servo's driver does, closes its own loop at
// every instant: τ = τ_p tanh((k_p (q* − q) + k_d (q̇* − q̇))/τ_p), saturating smoothly at its peak, its
// gains set by a bandwidth of 120 rad/s and the joint's reduced inertia with the water's added mass,
// 1/(M⁻¹)_jj. (A loop closed at 200 Hz with its torque held was too weak against the water to be stable and
// stiff at once.) The servo's torque–speed line is not applied: the joints turn at a fraction of their
// no-load speeds.
//
// The simulation: the variational integrator of the discrete action of the body and the water (the drag
// and the motors by Lagrange–d'Alembert), or semi-implicit Euler or Runge–Kutta on the Euler–Lagrange
// equations solved by jets (physis lagrangian.hpp). The energy's account: E = T + T_f, the body's and the
// water's kinetic energy, changes by the power of the drag, the motors and the gears' friction; its residual
// is the method's own error. The wake (what the fish and the water make together): the force each strip of
// the peduncle and the fin gives the water, where it gives it, kept with each frame.

#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include <physis/mechanics/fluid.hpp>
#include <physis/mechanics/lagrangian.hpp>
#include <physis/mechanics/multibody.hpp>
#include <physis/mechanics/variational.hpp>
#include <robotics/actuation/motor.hpp>

namespace robotics::swimming {

namespace pm = physis::mechanics;

inline constexpr std::size_t DOF = 7, JOINTS = 4;
using Coordinates = std::array<double, DOF>;
enum Coordinate : std::size_t { X = 0, Y = 1, HEADING = 2, J1 = 3, J2 = 4, J3 = 5, J4 = 6 };
inline constexpr std::array<const char*, JOINTS> JOINT_NAMES = {"neck", "waist", "peduncle", "tail"};

struct Fish {
  std::vector<pm::PlanarLink> links;
  std::vector<pm::SlenderPart> parts;  // one a link
  std::array<actuation::GearedMotor, JOINTS> motors;
  std::array<double, JOINTS> limit{};  // |q| at most (rad)
  double coulomb_speed = 0.05;
  double mass() const {
    double m = 0;
    for (const auto& l : links) m += l.mass;
    return m;
  }
  double length() const {
    double s = 0;
    for (const auto& p : parts) s += p.length;
    return s;
  }
};

inline Fish fish() {
  Fish f;
  const auto rod = [](double m, double l) { return m * l * l / 12; };
  f.links = {
      {"head", -1, {0, 0}, 0.0, 1.20, {0.0, 0.0}, 1.20 * (0.25 * 0.25 + 0.10 * 0.10) / 12, 0},
      {"body", 0, {-0.125, 0.0}, std::numbers::pi, 0.35, {0.07, 0.0}, rod(0.35, 0.14), 3e-4},
      {"tail body", 1, {0.14, 0.0}, 0.0, 0.20, {0.06, 0.0}, rod(0.20, 0.12), 3e-4},
      {"peduncle", 2, {0.12, 0.0}, 0.0, 0.06, {0.035, 0.0}, rod(0.06, 0.07), 2e-4},
      {"caudal fin", 3, {0.07, 0.0}, 0.0, 0.04, {0.05, 0.0}, rod(0.04, 0.10), 1e-4},
  };
  f.parts = {{0, {-0.125, 0.0}, 0.25, 0.12}, {1, {0, 0}, 0.14, 0.11}, {2, {0, 0}, 0.12, 0.09}, {3, {0, 0}, 0.07, 0.04}, {4, {0, 0}, 0.10, 0.18}};
  f.motors = {{{3.0, 8, 3e-4, 0.01, 0.02, 1.0, 1.0}, {3.0, 8, 3e-4, 0.01, 0.02, 1.0, 1.0}, {2.0, 10, 2e-4, 0.005, 0.01, 0.6, 1.2}, {1.5, 12, 1e-4, 0.003, 0.01, 0.4, 1.5}}};
  f.limit = {0.6, 0.8, 1.0, 1.1};
  return f;
}

// The water where the fish swims: its density and current (the world's), and the drag the fish's shape has.
inline pm::Fluid water(double density = 1025, std::array<double, 2> current = {0, 0}) {
  pm::Fluid w;
  w.density = density, w.normal_drag = 1.2, w.friction_drag = 0.03, w.current = current, w.strips = 4;
  return w;
}

// ---- the controller ----------------------------------------------------------------------------------------

enum class Behaviour { hold, swim, stop };
inline const char* behaviour_name(Behaviour b) { return b == Behaviour::hold ? "hold" : b == Behaviour::swim ? "swim" : "stop"; }

struct Command {
  Behaviour behaviour = Behaviour::hold;
  double speed = 0.3;  // m/s
  double turn = 0;     // −1 … 1: right … left
};

struct Gains {
  std::array<double, JOINTS> amplitude = {0.05, 0.10, 0.20, 0.32};  // rad
  double lag = 1.1;                   // rad between joints: a wave of about one body length
  double frequency = 1.2;             // Hz, to start with
  double frequency_min = 0.4, frequency_max = 3.0;
  double speed_integral = 0.8;        // Hz per (m/s) per second
  double turn_bias = 0.35;            // rad on the neck and the waist at a full turn
  double bandwidth = 120, damping_ratio = 0.9;  // rad/s, the servos'
  double ramp = 1.0;                  // s, for the wave's amplitude to grow
};

class Controller {
 public:
  // inertia: each joint's, with the water's added mass (the Lagrangian's Hessian in v at rest).
  Controller(const Fish& f, std::array<double, JOINTS> inertia, Gains g = {}) : f_(&f), inertia_(inertia), k_(g), freq_(g.frequency) {}
  void command(const Command& c) {
    if (c.behaviour == Behaviour::swim && cmd_.behaviour != Behaviour::swim) started_ = -1;
    cmd_ = c;
  }
  const Command& current() const { return cmd_; }
  double frequency() const { return cmd_.behaviour == Behaviour::swim ? freq_ : 0.0; }
  double phase() const { return phase_; }

  // The servos' targets at time t: each joint's angle and rate; `forward` the head's speed along its heading.
  struct Targets {
    std::array<double, JOINTS> angle{}, rate{};
  };
  Targets step(double t, double forward, double dt) {
    std::array<double, JOINTS> target{}, rate{};
    if (cmd_.behaviour == Behaviour::swim) {
      if (started_ < 0) started_ = t;
      freq_ = std::clamp(freq_ + k_.speed_integral * (cmd_.speed - forward) * dt, k_.frequency_min, k_.frequency_max);
      phase_ += 2 * std::numbers::pi * freq_ * dt;
      const double grow = std::clamp((t - started_) / k_.ramp, 0.0, 1.0);
      for (std::size_t j = 0; j < JOINTS; j++) {
        const double a = grow * k_.amplitude[j], ph = phase_ - double(j) * k_.lag;
        // a turn to the left (turn > 0) biases the front joints clockwise: so biased, the beating body turns to
        // the left (the checks see 60° in 8 s); the sign is found, not derived
        target[j] = a * std::sin(ph) - (j < 2 ? k_.turn_bias * std::clamp(cmd_.turn, -1.0, 1.0) : 0.0);
        rate[j] = a * 2 * std::numbers::pi * freq_ * std::cos(ph);
      }
    }
    return {target, rate};
  }
  // The servos' gains: k_p = I ω², k_d = 2ζ I ω.
  std::array<double, JOINTS> kp() const {
    std::array<double, JOINTS> k;
    for (std::size_t j = 0; j < JOINTS; j++) k[j] = inertia_[j] * k_.bandwidth * k_.bandwidth;
    return k;
  }
  std::array<double, JOINTS> kd() const {
    std::array<double, JOINTS> k;
    for (std::size_t j = 0; j < JOINTS; j++) k[j] = 2 * k_.damping_ratio * inertia_[j] * k_.bandwidth;
    return k;
  }

 private:
  const Fish* f_;
  std::array<double, JOINTS> inertia_;
  Gains k_;
  Command cmd_{};
  double freq_, phase_ = 0, started_ = -1;
};

// ---- the simulation ----------------------------------------------------------------------------------------

enum class Method { variational, semi_implicit, runge_kutta };
inline const char* method_name(Method m) { return m == Method::variational ? "variational" : m == Method::semi_implicit ? "semi_implicit" : "runge_kutta"; }
inline Method method_of(const std::string& s) {
  if (s == "variational") return Method::variational;
  if (s == "runge_kutta" || s == "rk4") return Method::runge_kutta;
  if (s == "semi_implicit" || s == "euler") return Method::semi_implicit;
  throw std::invalid_argument("robotics::swimming: the method is variational, semi_implicit or runge_kutta");
}

struct SimulationConfig {
  Method method = Method::runge_kutta;
  double step = 1e-3, control_period = 5e-3, frame_period = 1.0 / 30;
  double density = 1025;               // kg/m³, the water's where it swims
  std::array<double, 2> current{0, 0}; // m/s, in the plane's axes
};

// The force a strip gives the water, where: a mark of the wake.
struct WakeMark {
  double x, y, fx, fy;
};

struct Frame {
  double t = 0;
  Coordinates q{};
  std::array<double, JOINTS> torque{};
  std::vector<WakeMark> wake;
};

// The forces on the fish but its inertia and the water's kinetic energy: the drag, the servos' torques
// toward their held targets, the gears' friction; over any number type.
struct FishForces {
  const pm::PlanarTree<DOF>* tree;
  const Fish* fish;
  const pm::Fluid* water;
  const Controller::Targets* target;
  std::array<double, JOINTS> kp{}, kd{};
  template <class S>
  S servo(std::size_t j, const S& q, const S& v) const {
    using std::tanh;
    const double peak = fish->motors[j].peak;
    return tanh((kp[j] * (target->angle[j] - q) + kd[j] * (target->rate[j] - v)) / peak) * peak;
  }
  template <class S>
  std::array<S, DOF> operator()(double, const std::array<S, DOF>& q, const std::array<S, DOF>& v) const {
    using std::tanh;
    auto Q = pm::fluid_drag<DOF>(*tree, tree->kinematics(q, v), fish->parts, *water);
    for (std::size_t j = 0; j < JOINTS; j++) {
      const auto& m = fish->motors[j];
      const S& w = v[J1 + j];
      Q[J1 + j] += servo(j, q[J1 + j], w) - w * m.viscous - tanh(w / fish->coulomb_speed) * m.coulomb;
      // the joint's stop: a stiff spring beyond its range
      const S& x = q[J1 + j];
      if (x > fish->limit[j]) Q[J1 + j] -= (x - fish->limit[j]) * 20.0;
      if (x < -fish->limit[j]) Q[J1 + j] -= (x + fish->limit[j]) * 20.0;
    }
    return Q;
  }
};

class Simulation {
  using L = pm::TreeInFluid<DOF>;
  using VI = pm::VariationalIntegrator<DOF, L, FishForces>;

 public:
  Simulation(const Fish& f, SimulationConfig c)
      : f_(f), c_(c), tree_(pm::Base::floating, f_.links, 0.0), water_(water(c.density, c.current)), forces_{&tree_, &f_, &water_, &targets_},
        ctl_(f_, joint_inertia()) {
    forces_.kp = ctl_.kp(), forces_.kd = ctl_.kd();
    if (!(c.step > 0) || !(c.control_period >= c.step)) throw std::invalid_argument("robotics::swimming::Simulation: a positive step within the controller's period");
    q_.fill(0), v_.fill(0);
    if (c_.method == Method::variational) vi_ = std::make_unique<VI>(pm::DiscreteLagrangian<DOF, L>(lagrangian(), c_.step), forces_, 1e-11, 60), vs_ = vi_->start(t_, q_, v_);
    e0_ = energy();
    record(true);
  }

  void command(const Command& cmd) { ctl_.command(cmd); }
  const Controller& controller() const { return ctl_; }
  const pm::PlanarTree<DOF>& tree() const { return tree_; }
  const Fish& fish() const { return f_; }
  const SimulationConfig& config() const { return c_; }
  double time() const { return t_; }
  const Coordinates& q() const { return q_; }
  const Coordinates& v() const { return v_; }
  const std::vector<Frame>& frames() const { return frames_; }
  std::size_t frames_dropped() const { return dropped_; }
  void trim_frames(std::size_t keep_from) {
    const std::size_t k = std::min(keep_from, frames_.size());
    frames_.erase(frames_.begin(), frames_.begin() + std::ptrdiff_t(k));
    dropped_ += k;
  }
  std::size_t fallbacks() const { return fallbacks_; }

  void advance(double seconds) {
    const auto wall0 = std::chrono::steady_clock::now();
    const double end = t_ + seconds;
    while (t_ < end - 0.5 * c_.step) {
      if (t_ >= next_control_ - 1e-12) control();
      const double P0 = power();
      step();
      work_ += 0.5 * c_.step * (P0 + power());
      odometer_ += c_.step * std::hypot(v_[X], v_[Y]);
      // the servos' torques now, and the energy their motors draw
      for (std::size_t j = 0; j < JOINTS; j++) {
        torque_[j] = forces_.servo(j, q_[J1 + j], v_[J1 + j]);
        const double p = f_.motors[j].electrical_power(torque_[j], v_[J1 + j]);
        electrical_ += std::max(p, 0.0) * c_.step;
        peak_power_ = std::max(peak_power_, p);
      }
      record(false);
    }
    wall_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
    simulated_ += seconds;
  }

  // The energy of the body and the water, and its account.
  double energy() const { return lagrangian()(0.0, q_, v_); }  // L = T + T_f, quadratic in v: E = L
  double energy_residual() const { return energy() - e0_ - work_; }
  double electrical_energy() const { return electrical_; }
  double peak_power() const { return peak_power_; }
  double odometer() const { return odometer_; }
  double cost() const { return simulated_ > 0 ? wall_ / simulated_ : 0.0; }
  const std::array<double, JOINTS>& torque() const { return torque_; }
  // The head's speed along its heading, relative to the ground (the plane), and its heading.
  double forward() const { return v_[X] * std::cos(q_[HEADING]) + v_[Y] * std::sin(q_[HEADING]); }
  double cost_of_transport() const { return odometer_ > 0.05 ? electrical_ / (f_.mass() * 9.81 * odometer_) : 0.0; }
  // The tail's tip, in the plane.
  std::array<double, 2> tail_tip() const {
    const auto st = tree_.kinematics(q_);
    return tree_.point(st, 4, {0.10, 0.0});
  }
  // Each joint's inertia with the water's added mass, at rest straight, as the joint feels it on a free body:
  // 1/(M⁻¹)_jj, M = ∂²L/∂v∂v, the reduced inertia (the head turns back as the tail turns), not M_jj, which
  // holds the head still and overstates it many times; gains from M_jj made the held torque unstable.
  std::array<double, JOINTS> joint_inertia() const {
    std::array<double, 2 * DOF> z{};
    const auto x = calx::differentiation::variables<2 * DOF, 2>(z);
    std::array<calx::differentiation::Taylor<2 * DOF, 2>, DOF> qj, vj;
    for (std::size_t k = 0; k < DOF; k++) qj[k] = x[k], vj[k] = x[DOF + k];
    const auto l = lagrangian()(0.0, qj, vj);
    std::array<std::array<double, DOF>, DOF> M;
    for (std::size_t i = 0; i < DOF; i++)
      for (std::size_t j = 0; j < DOF; j++) M[i][j] = l.h[DOF + i][DOF + j];
    std::array<double, JOINTS> I;
    for (std::size_t j = 0; j < JOINTS; j++) {
      Coordinates e{};
      e[J1 + j] = 1;
      I[j] = 1 / calx::linalg::solve<DOF>(M, e)[J1 + j];
    }
    return I;
  }

 private:
  L lagrangian() const { return L{&tree_, &f_.parts, &water_}; }

  void control() {
    targets_ = ctl_.step(t_, forward(), c_.control_period);
    next_control_ += c_.control_period;
  }

  void step() {
    const auto Q = [&](double t, const Coordinates& q, const Coordinates& v) { return forces_(t, q, v); };
    if (c_.method == Method::variational) {
      try {
        vi_->step(vs_);
        t_ = vs_.t, q_ = vs_.q;
        v_ = solve_velocity(q_, vs_.p);
      } catch (const std::runtime_error&) {
        fallbacks_++;
        pm::semi_implicit_euler_lagrangian<DOF>(lagrangian(), Q, t_, q_, v_, c_.step);
        vs_ = vi_->start(t_, q_, v_);
      }
      return;
    }
    if (c_.method == Method::semi_implicit) pm::semi_implicit_euler_lagrangian<DOF>(lagrangian(), Q, t_, q_, v_, c_.step);
    else pm::runge_kutta4_lagrangian<DOF>(lagrangian(), Q, t_, q_, v_, c_.step);
  }

  // The velocity from the discrete momentum p = ∂L/∂v (L quadratic in v: p = M(q) v, M the Hessian in v).
  Coordinates solve_velocity(const Coordinates& q, const Coordinates& p) const {
    std::array<double, 2 * DOF> z{};
    for (std::size_t k = 0; k < DOF; k++) z[k] = q[k];
    const auto x = calx::differentiation::variables<2 * DOF, 2>(z);
    std::array<calx::differentiation::Taylor<2 * DOF, 2>, DOF> qj, vj;
    for (std::size_t k = 0; k < DOF; k++) qj[k] = x[k], vj[k] = x[DOF + k];
    const auto l = lagrangian()(0.0, qj, vj);
    std::array<std::array<double, DOF>, DOF> M;
    for (std::size_t i = 0; i < DOF; i++)
      for (std::size_t j = 0; j < DOF; j++) M[i][j] = l.h[DOF + i][DOF + j];
    return calx::linalg::solve<DOF>(M, p);
  }

  double power() const {
    const auto Q = forces_(t_, q_, v_);
    double P = 0;
    for (std::size_t k = 0; k < DOF; k++) P += Q[k] * v_[k];
    return P;
  }

  void record(bool force) {
    if (!force && t_ < next_frame_ - 1e-12) return;
    next_frame_ = t_ + c_.frame_period;
    Frame fr;
    fr.t = t_, fr.q = q_, fr.torque = torque_;
    // The wake: the force each strip of the peduncle and the fin gives the water (minus the water's on it).
    const auto st = tree_.kinematics(q_, v_);
    for (std::size_t i = 3; i < f_.parts.size(); i++) {
      const auto& b = f_.parts[i];
      const double ds = b.length / water_.strips;
      for (int k = 0; k < water_.strips; k++) {
        const std::array<double, 2> Xs{b.from[0] + (k + 0.5) * ds, b.from[1]};
        const auto p = tree_.point(st, b.link, Xs);
        std::vector<pm::SlenderPart> one = {{b.link, {Xs[0] - ds / 2, Xs[1]}, ds, b.height}};
        pm::Fluid w = water_;
        w.strips = 1;
        // the strip's force from the drag's generalized force on the base's translation: its x and y
        const auto g = pm::fluid_drag<DOF>(tree_, st, one, w);
        fr.wake.push_back({p[0], p[1], -g[X], -g[Y]});
      }
    }
    frames_.push_back(std::move(fr));
    if (frames_.size() > 12000) trim_frames(6000);
  }

  Fish f_;
  SimulationConfig c_;
  pm::PlanarTree<DOF> tree_;
  pm::Fluid water_;
  std::array<double, JOINTS> torque_{};
  Controller::Targets targets_{};
  FishForces forces_;
  Controller ctl_;
  std::unique_ptr<VI> vi_;
  pm::PhaseState<DOF> vs_;
  double t_ = 0, next_control_ = 0, next_frame_ = 0;
  Coordinates q_{}, v_{};
  double e0_ = 0, work_ = 0, electrical_ = 0, peak_power_ = 0, odometer_ = 0, wall_ = 0, simulated_ = 0;
  std::vector<Frame> frames_;
  std::size_t dropped_ = 0, fallbacks_ = 0;
};

}  // namespace robotics::swimming
