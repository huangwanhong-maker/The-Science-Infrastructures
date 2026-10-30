// A dog: a four-legged robot of 10.2 kg in its sagittal plane, its body, its trotting controller and its
// simulation by three methods of one motion, on a ground it learns (ADR 0059; theory/quadruped.md).
//
// The body. A torso (the base) and four legs, each a thigh and a shank with a point foot: the left and right
// fore legs at the shoulders, the left and right hind legs at the hips, all in the torso's plane, the left ones
// drawn behind the right. Its configuration space is SE(2) × T⁸. The two legs of a pair are separate links, so
// that it can trot: its diagonal pairs (left fore with right hind, right fore with left hind) stand in turn, and
// at every instant a fore foot and a hind foot are down, spanning its body in the plane.
//
// The motors (actuation/motor.hpp), geared; the feet on physis's compliant contact, on the ground along the
// line (the world's: legged/ground.hpp's tabulated truth), the dog's belief of it from its feet as the rabbit's.
//
// The controller, at 1 kHz:
//  - standing: each joint to its posture by a PD law, the weight carried by feed-forward, −Jᵀ f for each
//    foot's share f;
//  - trotting: a gait clock of period T; the diagonal pair in stance gives the body the force and the moment
//    asked, the weight, the height and the pitch held, the speed servoed, distributed between its fore foot
//    and its hind foot (each the same share of the forward force, the vertical ones solving the balance of the
//    moment), within the friction it believes is there; the swing pair follows an arc of lift to Raibert's
//    placement ahead of the hip, x_f = v T_s/2 + k_v (v − v*), by each leg's inverse kinematics and a PD law;
//  - stopping: the speed asked is zero, and when slow it stands.
//
// The simulation: semi-implicit Euler, Runge–Kutta (on Newton–Euler's equations) or the variational
// integrator of the discrete action (the contacts and the motors by Lagrange–d'Alembert), with the energy's
// account, as the rabbit's (legged/simulation.hpp).

#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <functional>
#include <memory>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <calx/linalg/small.hpp>
#include <calx/random/stream.hpp>
#include <physis/mechanics/contact.hpp>
#include <physis/mechanics/multibody.hpp>
#include <physis/mechanics/variational.hpp>
#include <robotics/actuation/motor.hpp>
#include <robotics/legged/ground.hpp>

namespace robotics::quadruped {

namespace pm = physis::mechanics;
using legged::GroundBelief;
using legged::GroundPresumption;
using legged::TabulatedGround;

inline constexpr std::size_t DOF = 11, LEGS = 4, JOINTS = 8;
inline constexpr double DEG = std::numbers::pi / 180;
using Coordinates = std::array<double, DOF>;
enum Coordinate : std::size_t { X = 0, Z = 1, PITCH = 2 };
// the legs, and the coordinates of each leg's hip and knee: 3 + 2 leg, 4 + 2 leg
enum Leg : std::size_t { RF = 0, LF = 1, RH = 2, LH = 3 };
inline constexpr std::array<const char*, LEGS> LEG_NAMES = {"right fore", "left fore", "right hind", "left hind"};
inline constexpr std::size_t hip(std::size_t leg) { return 3 + 2 * leg; }
inline constexpr std::size_t knee(std::size_t leg) { return 4 + 2 * leg; }
inline constexpr bool fore(std::size_t leg) { return leg < 2; }

struct Dog {
  std::vector<pm::PlanarLink> links;  // torso, then for each leg its thigh and shank
  std::array<actuation::GearedMotor, JOINTS> motors;
  double thigh = 0.20, shank = 0.20, shoulder_x = 0.18, hip_x = -0.18;
  double width = 0.20;                // between the left and right legs, for drawing
  double coulomb_speed = 0.05;
  double mass() const {
    double m = 0;
    for (const auto& l : links) m += l.mass;
    return m;
  }
  std::array<double, 2> mount(std::size_t leg) const { return {fore(leg) ? shoulder_x : hip_x, 0.0}; }
};

inline Dog dog() {
  Dog d;
  const auto rod = [](double m, double l) { return m * l * l / 12; };
  d.links = {{"torso", -1, {0, 0}, 0.0, 7.0, {0.0, 0.0}, 7.0 * (0.40 * 0.40 + 0.10 * 0.10) / 12, 0}};
  for (std::size_t l = 0; l < LEGS; l++) {
    const int parent = int(d.links.size());
    d.links.push_back({std::string(LEG_NAMES[l]) + " thigh", 0, d.mount(l), 0.0, 0.60, {0.09, 0.0}, rod(0.60, 0.20), 4e-3});
    d.links.push_back({std::string(LEG_NAMES[l]) + " shank", parent, {0.20, 0.0}, 0.0, 0.20, {0.08, 0.0}, rod(0.20, 0.20), 4e-3});
  }
  for (std::size_t j = 0; j < JOINTS; j++)
    d.motors[j] = j % 2 == 0 ? actuation::GearedMotor{20.0, 25, 4e-3, 0.02, 0.1, 1.2, 0.3} : actuation::GearedMotor{30.0, 20, 4e-3, 0.02, 0.1, 1.8, 0.3};
  return d;
}

// A two-link leg's inverse kinematics in the torso's frame, the knee behind (as a robot dog's): the angles of
// the thigh and the shank from the torso's x that put the foot at T from the hip.
inline std::optional<std::array<double, 2>> leg_ik(double l1, double l2, std::array<double, 2> T) {
  double d = std::hypot(T[0], T[1]);
  if (d > (l1 + l2) * 0.999) T = {T[0] * (l1 + l2) * 0.999 / d, T[1] * (l1 + l2) * 0.999 / d}, d = (l1 + l2) * 0.999;
  if (d < 0.05) return std::nullopt;
  const double base = std::atan2(T[1], T[0]), c = std::clamp((l1 * l1 + d * d - l2 * l2) / (2 * l1 * d), -1.0, 1.0);
  const double a1 = base - std::acos(c);  // the knee on the side behind the hip–foot line
  const double a2 = std::atan2(T[1] - l1 * std::sin(a1), T[0] - l1 * std::cos(a1));
  return std::array<double, 2>{a1, a2};
}

// ---- the controller -------------------------------------------------------------------------------------

enum class Behaviour { stand, trot, stop };
inline const char* behaviour_name(Behaviour b) { return b == Behaviour::stand ? "stand" : b == Behaviour::trot ? "trot" : "stop"; }

struct Command {
  Behaviour behaviour = Behaviour::stand;
  double speed = 0;  // m/s, forward
};

struct Gains {
  double height = 0.30;         // m, the hips above the ground
  double period = 0.36;         // s, a stride: each pair stands half of it
  double lift = 0.06;           // m, a swing foot's rise
  double kz = 900, bz = 60;     // N/m, N s/m: the height
  double kth = 120, bth = 12;   // N m/rad, N m s/rad: the pitch
  double kv = 4;                // 1/s: the speed's error to a forward acceleration
  double ki = 4;                // 1/s²: its integral's (Raibert's placement and the force leave a steady error)
  double placement = 0.06;      // s, k_v
  double friction = 0.6;        // the friction it asks of the ground, at most
  std::array<double, 2> inertia = {0.03, 0.008};  // kg m², about the hip and the knee with the leg's links
  double stand_bandwidth = 30, swing_bandwidth = 40, damping_ratio = 0.9;
};

struct Observation {
  double t = 0;
  Coordinates q{}, v{};
  std::array<double, LEGS> foot_force{};
  double friction = 0;  // what it believes it may count on, or 0
};

class Controller {
 public:
  Controller(const Dog& d, double g, Gains k = {}) : d_(&d), g_(g), k_(k) {}
  void attach(const pm::PlanarTree<DOF>& tree) { tree_ = &tree; }
  void command(const Command& c) {
    if (c.behaviour == Behaviour::trot && cmd_.behaviour != Behaviour::trot) clock0_ = -1;
    cmd_ = c;
  }
  const Command& current() const { return cmd_; }
  bool trotting() const { return cmd_.behaviour != Behaviour::stand; }
  // Whether a leg stands now (by the gait's clock when trotting).
  bool in_stance(std::size_t leg) const { return stance_[leg]; }

  // The standing posture: each foot beneath its mount, the hips `height` above it.
  std::array<double, JOINTS> stand_posture() const {
    std::array<double, JOINTS> a{};
    for (std::size_t l = 0; l < LEGS; l++)
      if (const auto ik = leg_ik(d_->thigh, d_->shank, {0.0, -k_.height})) a[2 * l] = (*ik)[0], a[2 * l + 1] = (*ik)[1] - (*ik)[0];
    return a;
  }
  Coordinates standing_state(double x0, double ground) const {
    Coordinates q{};
    q[X] = x0, q[Z] = ground + k_.height - 0.002, q[PITCH] = 0;
    const auto a = stand_posture();
    for (std::size_t j = 0; j < JOINTS; j++) q[3 + j] = a[j];
    return q;
  }

  std::array<double, JOINTS> step(const Observation& o) {
    const auto& tree = *tree_;
    const auto st = tree.kinematics(o.q, o.v);
    const auto com = tree.centre_of_mass(o.q);
    const auto P = tree.linear_momentum(o.q, o.v);
    const double m = d_->mass(), vx = P[0] / m, vz = P[1] / m;
    const double th = o.q[PITCH], om = o.v[PITCH];
    const double mu = o.friction > 0 ? std::min(k_.friction, 0.75 * o.friction) : k_.friction;
    double ground = 0;  // the feet's level: the lowest foot's height
    std::array<std::array<double, 2>, LEGS> foot{};
    for (std::size_t l = 0; l < LEGS; l++) foot[l] = tree.point(st, 2 + 2 * l, {d_->shank, 0.0});
    ground = std::min(std::min(foot[0][1], foot[1][1]), std::min(foot[2][1], foot[3][1]));
    // the gait's clock
    double phase = 0;
    if (cmd_.behaviour == Behaviour::trot || cmd_.behaviour == Behaviour::stop) {
      if (clock0_ < 0) clock0_ = o.t;
      phase = std::fmod((o.t - clock0_) / k_.period, 1.0);
      const bool a = phase < 0.5;  // pair A: right fore and left hind; pair B: left fore and right hind
      stance_ = {a, !a, !a, a};
      if (cmd_.behaviour == Behaviour::stop && std::fabs(vx) < 0.05 && std::fabs(phase - 0.5) < 0.02) cmd_.behaviour = Behaviour::stand;
    }
    if (cmd_.behaviour == Behaviour::stand) stance_ = {true, true, true, true};
    // The force and moment asked of the feet down, about the centre of mass.
    const double hz = o.q[Z] - ground - k_.height;  // the torso's height above the feet, against the one asked
    double Fz = m * g_ + k_.kz * (0.0 - hz) - k_.bz * vz;
    const double v_target = cmd_.behaviour == Behaviour::trot ? cmd_.speed : 0.0;
    const double dt = last_t_ < 0 ? 0.0 : std::clamp(o.t - last_t_, 0.0, 0.01);
    last_t_ = o.t;
    if (cmd_.behaviour == Behaviour::stand) speed_integral_ = 0;
    else speed_integral_ = std::clamp(speed_integral_ + (v_target - vx) * dt, -0.2, 0.2);
    double Fx = m * (k_.kv * (v_target - vx) + k_.ki * speed_integral_);
    const double My = -k_.kth * th - k_.bth * om;
    Fz = std::max(Fz, 0.0);
    Fx = std::clamp(Fx, -mu * Fz * 0.8, mu * Fz * 0.8);
    std::array<std::array<double, 2>, LEGS> f{};
    int n = 0;
    for (std::size_t l = 0; l < LEGS; l++) n += stance_[l];
    // Distribute: each stance foot an equal share of the forward force; the vertical forces of the fore and
    // the hind feet so that their moment about the centre is My (a foot's share split between the pair's
    // two legs when all four stand).
    {
      double xf = 0, xh = 0, zf = 0, zh = 0;
      int nf = 0, nh = 0;
      for (std::size_t l = 0; l < LEGS; l++)
        if (stance_[l]) {
          if (fore(l)) xf += foot[l][0] - com[0], zf += foot[l][1] - com[1], nf++;
          else xh += foot[l][0] - com[0], zh += foot[l][1] - com[1], nh++;
        }
      if (nf && nh) {
        xf /= nf, xh /= nh, zf /= nf, zh /= nh;
        const double fx = Fx / 2;  // the fore and the hind each half
        // fz_f + fz_h = Fz; xf fz_f − zf fx + xh fz_h − zh fx = My
        const double fzf = std::clamp((My + (zf + zh) * fx - xh * Fz) / (xf - xh), 0.0, Fz), fzh = Fz - fzf;
        for (std::size_t l = 0; l < LEGS; l++)
          if (stance_[l]) {
            const double share = 1.0 / (fore(l) ? nf : nh);
            const double fz = (fore(l) ? fzf : fzh) * share, fxl = std::clamp(fx * share, -mu * fz, mu * fz);
            f[l] = {fxl, fz};
          }
      }
    }
    std::array<double, JOINTS> tau{};
    const auto pd = [&](std::size_t j, double target, double rate, double w) {
      const double I = k_.inertia[j % 2];
      return I * w * w * (target - o.q[3 + j]) + 2 * k_.damping_ratio * I * w * (rate - o.v[3 + j]);
    };
    const auto stand = stand_posture();
    for (std::size_t l = 0; l < LEGS; l++) {
      const std::size_t h = 2 * l, k = 2 * l + 1;
      if (stance_[l]) {
        // −Jᵀ f at the foot, and a light posture servo for the slack of the model
        const auto J = tree.jacobian(st, 2 + 2 * l, foot[l]);
        tau[h] = -(J[0][hip(l)] * f[l][0] + J[1][hip(l)] * f[l][1]);
        tau[k] = -(J[0][knee(l)] * f[l][0] + J[1][knee(l)] * f[l][1]);
        if (cmd_.behaviour == Behaviour::stand) tau[h] += pd(h, stand[h], 0, k_.stand_bandwidth), tau[k] += pd(k, stand[k], 0, k_.stand_bandwidth);
        lift_off_[l] = {foot[l][0] - o.q[X], foot[l][1] - o.q[Z]};
      } else {
        // the swing: an arc from where it left the ground to the placement ahead of its mount, down to the
        // feet's level, in world axes from the torso's origin, turned into the torso's frame for the leg's
        // inverse kinematics
        const double s = std::clamp((phase < 0.5 ? phase : phase - 0.5) / 0.5, 0.0, 1.0);
        const double c = std::cos(th), sn = std::sin(th);
        const auto mnt = d_->mount(l);
        const double mx = c * mnt[0] - sn * mnt[1], mz = sn * mnt[0] + c * mnt[1];  // the mount from the torso's origin
        const double xf = mx + std::clamp(vx * k_.period / 4 + k_.placement * (vx - v_target), -0.10, 0.10);
        const double x = lift_off_[l][0] + (xf - lift_off_[l][0]) * (0.5 - 0.5 * std::cos(std::numbers::pi * s));
        const double z = (ground - o.q[Z]) + k_.lift * std::sin(std::numbers::pi * s);
        const std::array<double, 2> rel = {x - mx, z - mz};
        const std::array<double, 2> T = {c * rel[0] + sn * rel[1], -sn * rel[0] + c * rel[1]};
        if (const auto ik = leg_ik(d_->thigh, d_->shank, T)) {
          tau[h] = pd(h, (*ik)[0], 0, k_.swing_bandwidth);
          tau[k] = pd(k, (*ik)[1] - (*ik)[0], 0, k_.swing_bandwidth);
        }
      }
    }
    for (std::size_t j = 0; j < JOINTS; j++) tau[j] = d_->motors[j].available(tau[j], o.v[3 + j]);
    return tau;
  }

 private:
  const Dog* d_;
  const pm::PlanarTree<DOF>* tree_ = nullptr;
  double g_;
  Gains k_;
  Command cmd_{};
  double clock0_ = -1, last_t_ = -1, speed_integral_ = 0;
  std::array<bool, LEGS> stance_{true, true, true, true};
  std::array<std::array<double, 2>, LEGS> lift_off_{};
};

// ---- the simulation -----------------------------------------------------------------------------------------

enum class Method { variational, semi_implicit, runge_kutta };
inline const char* method_name(Method m) { return m == Method::variational ? "variational" : m == Method::semi_implicit ? "semi_implicit" : "runge_kutta"; }
inline Method method_of(const std::string& s) {
  if (s == "variational") return Method::variational;
  if (s == "runge_kutta" || s == "rk4") return Method::runge_kutta;
  if (s == "semi_implicit" || s == "euler") return Method::semi_implicit;
  throw std::invalid_argument("robotics::quadruped: the method is variational, semi_implicit or runge_kutta");
}

struct SimulationConfig {
  Method method = Method::semi_implicit;
  double gravity = 9.81, step = 2.5e-4, control_period = 1e-3, frame_period = 1.0 / 60;
  pm::GroundLine ground{};
  std::function<pm::CompliantContact(double x)> ground_at;  // the world's ground along the line, or a rock's
  GroundPresumption presumption{};
  std::uint64_t sensor_seed = 1;
};

struct Frame {
  double t = 0;
  Coordinates q{};
  std::array<double, LEGS> normal{};
  std::array<double, JOINTS> torque{};
  std::array<double, 2> com{};
  std::array<bool, LEGS> stance{};
};

struct FootReading {
  double t, x;
  int parameter;
  double value;
};

// The contact law of a rubber foot on rock, where no ground is given (geophysics ground's median for rock:
// the pad's stiffness and damping).
inline pm::CompliantContact rock() {
  pm::CompliantContact c;
  c.stiffness = 40 / std::pow(1e-3, 1.5), c.exponent = 1.5, c.damping = 1.5, c.friction = 0.8, c.slip_speed = 0.03;
  return c;
}

struct DogForces {
  const pm::PlanarTree<DOF>* tree;
  const Dog* dog;
  pm::GroundLine ground;
  const std::array<double, JOINTS>* torque;
  TabulatedGround* table = nullptr;
  // The ground's law at x, its slip speed at least 0.03 m/s: below it Coulomb's regularised friction is a damper
  // of μ f_n/v_s on the foot, 800 N s/m under a dog's 30 N on a foot of 0.2 kg; an explicit step integrates it only
  // while h μ f_n/(v_s m) < 2 (0.8 here at 0.25 ms; at 0.01 m/s and 0.5 ms, 5: the feet chattered).
  pm::CompliantContact law_at(double x) const {
    auto c = table ? table->at(x) : rock();
    c.slip_speed = std::max(c.slip_speed, 0.03);
    return c;
  }
  template <class S>
  std::array<S, DOF> operator()(double, const std::array<S, DOF>& q, const std::array<S, DOF>& v) const {
    using std::tanh;
    std::array<S, DOF> Q;
    Q.fill(S(0.0));
    for (std::size_t j = 0; j < JOINTS; j++) {
      const auto& m = dog->motors[j];
      const S& w = v[3 + j];
      Q[3 + j] += (*torque)[j] - w * m.viscous - tanh(w / dog->coulomb_speed) * m.coulomb;
    }
    const auto st = tree->kinematics(q, v);
    for (std::size_t l = 0; l < LEGS; l++) {
      const std::size_t link = 2 + 2 * l;
      const auto p = tree->point(st, link, {dog->shank, 0.0});
      const auto c = pm::contact_force(ground, law_at(calx::differentiation::value(p[0])), p, tree->point_velocity(st, link, {dog->shank, 0.0}));
      if (!c.touching) continue;
      const auto g = tree->generalized_force(st, link, p, c.force);
      for (std::size_t k = 0; k < DOF; k++) Q[k] += g[k];
    }
    return Q;
  }
  struct Contact {
    std::array<double, 2> place{}, force{};
    double normal = 0, tangential = 0, slip = 0, penetration = 0;
    bool touching = false;
  };
  std::array<Contact, LEGS> contacts(const Coordinates& q, const Coordinates& v) const {
    const auto st = tree->kinematics(q, v);
    const auto t = ground.tangent();
    std::array<Contact, LEGS> out;
    for (std::size_t l = 0; l < LEGS; l++) {
      const std::size_t link = 2 + 2 * l;
      auto& o = out[l];
      o.place = tree->point(st, link, {dog->shank, 0.0});
      const auto pv = tree->point_velocity(st, link, {dog->shank, 0.0});
      const auto c = pm::contact_force(ground, law_at(o.place[0]), o.place, pv);
      o.force = c.force, o.normal = c.normal, o.penetration = c.penetration, o.touching = c.touching;
      o.slip = pv[0] * t[0] + pv[1] * t[1], o.tangential = c.force[0] * t[0] + c.force[1] * t[1];
    }
    return out;
  }
};

class Simulation {
  using L = pm::TreeLagrangian<DOF, pm::NoPotential>;
  using VI = pm::VariationalIntegrator<DOF, L, DogForces>;

 public:
  Simulation(const Dog& d, SimulationConfig c)
      : d_(d), c_(c), tree_(pm::Base::floating, d_.links, c.gravity), forces_{&tree_, &d_, c.ground, &torque_}, ctl_(d_, c.gravity),
        belief_(c.presumption), noise_(c.sensor_seed, 0xD06u, 0u, 0u) {
    if (!(c.step > 0) || !(c.control_period >= c.step)) throw std::invalid_argument("robotics::quadruped::Simulation: a positive step within the controller's period");
    if (c_.ground_at) table_ = std::make_unique<TabulatedGround>(c_.ground_at), forces_.table = table_.get();
    ctl_.attach(tree_);
    q_ = ctl_.standing_state(0.0, c.ground.z0);
    v_.fill(0);
    if (c_.method == Method::variational) vi_ = std::make_unique<VI>(pm::DiscreteLagrangian<DOF, L>(L{&tree_, {}}, c_.step), forces_, 1e-11, 60), vs_ = vi_->start(t_, q_, v_);
    e0_ = energy();
    record(true);
  }

  void command(const Command& cmd) { ctl_.command(cmd); }
  const Controller& controller() const { return ctl_; }
  const pm::PlanarTree<DOF>& tree() const { return tree_; }
  const Dog& dog() const { return d_; }
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
  const GroundBelief& belief() const { return belief_; }
  const std::vector<FootReading>& readings() const { return readings_; }
  std::size_t readings_dropped() const { return readings_dropped_; }
  void trim_readings(std::size_t keep_from) {
    const std::size_t k = std::min(keep_from, readings_.size());
    readings_.erase(readings_.begin(), readings_.begin() + std::ptrdiff_t(k));
    readings_dropped_ += k;
  }
  pm::CompliantContact ground_at(double x) const { return forces_.law_at(x); }
  std::array<DogForces::Contact, LEGS> contacts() const { return forces_.contacts(q_, v_); }

  void advance(double seconds) {
    const auto wall0 = std::chrono::steady_clock::now();
    const double end = t_ + seconds;
    while (t_ < end - 0.5 * c_.step) {
      if (t_ >= next_control_ - 1e-12) control();
      const double P0 = power();
      step();
      work_ += 0.5 * c_.step * (P0 + power());
      record(false);
    }
    wall_ += std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
    simulated_ += seconds;
  }

  double energy() const { return tree_.energy(q_, v_); }
  double energy_residual() const { return energy() - e0_ - work_; }
  double electrical_energy() const { return electrical_; }
  double peak_power() const { return peak_power_; }
  double distance() const { return q_[X]; }
  double cost_of_transport() const { return std::fabs(distance()) > 0.2 ? electrical_ / (d_.mass() * c_.gravity * std::fabs(distance())) : 0.0; }
  double cost() const { return simulated_ > 0 ? wall_ / simulated_ : 0.0; }
  double max_penetration() const { return max_pen_; }
  const std::array<double, JOINTS>& torque() const { return torque_; }

 private:
  void control() {
    Observation o;
    o.t = t_, o.q = q_, o.v = v_;
    const auto cs = forces_.contacts(q_, v_);
    for (std::size_t l = 0; l < LEGS; l++) o.foot_force[l] = cs[l].normal;
    sense(cs);
    // the friction it may count on: its belief at its fore feet, less one deviation
    const auto mu = belief_.friction.at(0.5 * (cs[RF].place[0] + cs[LF].place[0]));
    o.friction = mu.median * std::exp(-mu.log_sd);
    torque_ = ctl_.step(o);
    next_control_ += c_.control_period;
    for (std::size_t j = 0; j < JOINTS; j++) {
      const double p = d_.motors[j].electrical_power(torque_[j], v_[3 + j]);
      electrical_ += std::max(p, 0.0) * c_.control_period;
      peak_power_ = std::max(peak_power_, p);
    }
  }
  // The feet as the rabbit's are sensors (legged/ground.hpp): a sliding foot reads the friction, a foot at its
  // deepest the load at 1 mm.
  void sense(const std::array<DogForces::Contact, LEGS>& cs) {
    for (std::size_t l = 0; l < LEGS; l++) {
      const auto& f = cs[l];
      const auto law = forces_.law_at(f.place[0]);
      if (f.touching && f.normal > 2.0 && std::fabs(f.slip) > 5 * law.slip_speed) {
        const double mu = std::fabs(f.tangential) / f.normal * (1 + 0.02 * noise_.normal());
        if (belief_.friction.observe(f.place[0], mu, 0.02)) record_reading(f.place[0], 0, mu);
      }
      if (f.touching && pen_[l][1] > pen_[l][0] && pen_[l][1] > f.penetration && fn_[l] > 5.0) {
        const double load = fn_[l] * std::pow(1e-3 / (pen_[l][1] * (1 + 0.10 * noise_.normal())), law.exponent);
        if (belief_.load.observe(x_[l], load, 0.20)) record_reading(x_[l], 1, load);
      }
      pen_[l][0] = pen_[l][1], pen_[l][1] = f.touching ? f.penetration : 0.0, fn_[l] = f.normal, x_[l] = f.place[0];
      max_pen_ = std::max(max_pen_, f.penetration);
    }
  }
  void record_reading(double x, int parameter, double value) {
    readings_.push_back({t_, x, parameter, value});
    if (readings_.size() > 4000) trim_readings(2000);
  }
  void step() {
    const auto Q = [&](double t, const Coordinates& q, const Coordinates& v) { return forces_(t, q, v); };
    if (c_.method == Method::variational) {
      try {
        vi_->step(vs_);
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
  double power() const {
    const auto Q = forces_(t_, q_, v_);
    double P = 0;
    for (std::size_t k = 0; k < DOF; k++) P += Q[k] * v_[k];
    return P;
  }
  void record(bool force) {
    if (!force && t_ < next_frame_ - 1e-12) return;
    next_frame_ = t_ + c_.frame_period;
    Frame f;
    f.t = t_, f.q = q_, f.torque = torque_, f.com = tree_.centre_of_mass(q_);
    const auto cs = forces_.contacts(q_, v_);
    for (std::size_t l = 0; l < LEGS; l++) f.normal[l] = cs[l].normal, f.stance[l] = ctl_.in_stance(l);
    frames_.push_back(f);
    if (frames_.size() > 20000) trim_frames(10000);
  }

  Dog d_;
  SimulationConfig c_;
  pm::PlanarTree<DOF> tree_;
  std::array<double, JOINTS> torque_{};
  DogForces forces_;
  Controller ctl_;
  std::unique_ptr<VI> vi_;
  pm::PhaseState<DOF> vs_;
  double t_ = 0, next_control_ = 0, next_frame_ = 0;
  Coordinates q_{}, v_{};
  double e0_ = 0, work_ = 0, electrical_ = 0, peak_power_ = 0, max_pen_ = 0, wall_ = 0, simulated_ = 0;
  std::vector<Frame> frames_;
  std::size_t dropped_ = 0, fallbacks_ = 0;
  std::unique_ptr<TabulatedGround> table_;
  GroundBelief belief_;
  calx::random::Stream noise_;
  std::array<std::array<double, 2>, LEGS> pen_{};
  std::array<double, LEGS> fn_{}, x_{};
  std::vector<FootReading> readings_;
  std::size_t readings_dropped_ = 0;
};

}  // namespace robotics::quadruped
