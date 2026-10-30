// The rabbit's controller: standing, hopping in place or forward at a speed, and stopping (ADR 0058;
// theory/legged.md §3).
//
// What it reads, as a robot does: its joints' angles and rates (encoders), its torso's pitch, pitch rate
// and velocity (an inertial unit and a state estimator), and whether each foot presses the ground (a foot
// switch, its normal force above a threshold). What it gives: a torque for each motor, which the motor
// clips to its torque–speed line. It runs at its own rate (500 Hz), its torques held between samples.
//
// Standing. The feet's places are chosen (the hind foot flat under the body, the fore paw under the
// shoulder) and the joints' angles found by the legs' inverse kinematics; each joint follows its angle by a
// proportional–derivative law whose gains are set by a bandwidth and the joint's inertia (Gains), and the
// weight is carried by feed-forward, −Jᵀ f for each foot's share f by the lever rule about the centre.
//
// Hopping: Raibert's three parts (Raibert 1986, Legged Robots That Balance, ch. 2–3), on a virtual leg from
// the toe to the centre of mass, so that the ground's force asked along it passes through the centre and
// turns nothing; the ankle held at its hopping angle and the forelegs lifted:
//  - height: in stance the leg is a spring k (L₀ − L) − b L̇ along it, and while it extends a thrust is
//    added, raised or lowered after each flight by the error of the apex reached against the one asked;
//    the force f asked at the toe becomes joint torques by τ = −Jᵀ f, J the toe's Jacobian in the hip, knee
//    and ankle (virtual work: the leg pushes the ground with −f, the ground the body with f);
//  - speed: in flight the toe is placed ahead of the centre by x_f = v T_s/2 + k_v (v − v* − b), the neutral
//    point of the last stance's duration T_s and a correction of the speed's error, b an integral of that
//    error hop by hop (the placement alone leaves a steady error);
//  - attitude: in stance a force across the leg, within the friction asked, gives the moment about the
//    centre M = −k_θ θ − k_L L_c, the pitch back to level and the angular momentum about the centre to zero
//    (a flight keeps what the stance leaves it, and the Moon's flights are long).
// A hop starts from a crouch (the ready phase: the body lowered over the flat hind foot, the forelegs
// lifted, eased in over ready_time), and the first push is a launch: the force m g (1 + h/ΔL), which does
// the work m g (ΔL + h) over the leg's extension ΔL, asked at the point of the sole beneath the centre while
// the heel is down. The command is refused, with its reason, when the ankle cannot hold that force over the
// foot's length (on the Earth, for this build).
// Phases: standing → ready → thrust (the launch) → flight → compression (the toe touches) → thrust (the leg
// starts to lengthen) → flight. Stopping asks for zero speed and no thrust, and stands when the rabbit is
// down and slow.
//
// The gains are set for the robot's 1.69 kg at the Moon's 1.62 m/s², scaled with g where the weight sets
// them.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <string>

#include <robotics/legged/rabbit.hpp>

namespace robotics::legged {

enum class Behaviour { stand, hop, stop };
enum class Phase { standing, ready, flight, compression, thrust };

inline const char* behaviour_name(Behaviour b) { return b == Behaviour::stand ? "stand" : b == Behaviour::hop ? "hop" : "stop"; }
inline const char* phase_name(Phase p) {
  return p == Phase::standing ? "standing" : p == Phase::ready ? "ready" : p == Phase::flight ? "flight" : p == Phase::compression ? "compression" : "thrust";
}

struct Observation {
  double t = 0;
  Coordinates q{}, v{};
  std::array<double, 3> foot_force{};  // toe, heel, paw: normal force (N)
  double friction = 0;                  // the friction it may count on at its toe (its belief), or 0: unknown
};

struct Command {
  Behaviour behaviour = Behaviour::stand;
  double speed = 0;     // m/s, forward
  double height = 0.15; // m, the apex of the torso's centre above its standing height
};

// A two-link leg's inverse kinematics in the torso's frame: the angles (φ₁ of the first link, φ₂ of the
// second, both from the torso's x) that put the far end at T from the first joint, the knee on the side
// `forward` (+1: the knee ahead of the line, as a rabbit's hind leg; −1: behind, as its elbow).
inline std::optional<std::array<double, 2>> two_link(double l1, double l2, std::array<double, 2> T, double forward) {
  double d = std::hypot(T[0], T[1]);
  const double reach = l1 + l2;
  if (d > reach * 0.999) {
    const double s = reach * 0.999 / d;
    T = {T[0] * s, T[1] * s}, d = reach * 0.999;
  }
  if (d < std::fabs(l1 - l2) * 1.001) return std::nullopt;
  const double base = std::atan2(T[1], T[0]);
  const double c = std::clamp((l1 * l1 + d * d - l2 * l2) / (2 * l1 * d), -1.0, 1.0);
  const double phi1 = base + forward * std::acos(c);
  const double phi2 = std::atan2(T[1] - l1 * std::sin(phi1), T[0] - l1 * std::cos(phi1));
  return std::array<double, 2>{phi1, phi2};
}

// The joints' gains are set by a bandwidth: kp = I ω², kd = 2ζ I ω, I the joint's inertia (its armature and
// the links it moves, about it), so that a light joint is not driven harder than its 2 ms of held torque
// allows (the product 2ζω T stays well below 1).
struct Gains {
  std::array<double, JOINTS> inertia = {3.0e-3, 1.6e-3, 4.5e-4, 5.0e-4, 3.0e-4};  // kg m², about each joint
  double stand_bandwidth = 50, swing_bandwidth = 60, arm_bandwidth = 15, damping_ratio = 0.9;  // rad/s
  double leg_stiffness = 420;      // N/m, the virtual leg at the Moon's g (scaled with g)
  double leg_damping = 3;          // N s/m (scaled with √g)
  double leg_rest = 0.20;          // m, from the toe to the centre of mass
  double hop_ankle = 70 * DEG;     // the ankle's angle in a hop, the foot stretched
  double attitude_kp = 1.0, attitude_kd = 0.5;    // N m/rad, N m s/rad (scaled with g)
  double momentum_rate = 30;       // 1/s: in stance the angular momentum about the centre is driven to zero at this rate
  double placement = 0.09;         // s, k_v
  double speed_integral = 0.4;     // the speed bias's change per hop per m/s of error
  double height_gain = 20;         // N/m: the thrust's change per metre of apex error (scaled with g)
  double thrust_max = 12;          // N (scaled with g)
  double friction = 0.5;           // the friction the controller asks of the ground, at most
  double touch = 0.3;              // N, the foot switch's threshold (scaled with g)
  double ready_time = 0.5;         // s, from standing to the crouch a hop starts from
  double ready_leg = 0.15;         // m, the virtual leg's length in that crouch
  double kp(std::size_t j, double w) const { return inertia[j] * w * w; }
  double kd(std::size_t j, double w) const { return 2 * damping_ratio * inertia[j] * w; }
};

class Controller {
 public:
  Controller(const Rabbit& r, double gravity, Gains g = {}) : r_(&r), g_(gravity), k_(g) {
    standing_height_ = 0.155;
    thrust_ = 0.15 * r.mass() * gravity;  // a first thrust; each flight's apex corrects it
    stance_time_ = 0.2;
  }

  void attach(const pm::PlanarTree<DOF>& tree) { tree_ = &tree; }
  // A command, unless the body cannot do it: a hop needs the launch's force, m g (1 + h/ΔL), at the toe,
  // which the ankle must hold over the foot's length; on the Earth this build's ankle cannot, and the
  // command is refused with the reason (refusal()).
  void command(Command c) {
    refusal_.clear();
    if (c.behaviour == Behaviour::hop) {
      const double F = r_->mass() * g_ * (1 + c.height / std::max(k_.leg_rest - k_.ready_leg, 0.02));
      const double need = F * 0.10, have = r_->actuators[2].peak;
      if (need > have) {
        refusal_ = "a hop of " + std::to_string(int(std::lround(c.height * 100))) + " cm needs " + std::to_string(F).substr(0, 4) +
                   " N at the toe, " + std::to_string(need).substr(0, 4) + " N m at the ankle; its motor gives " + std::to_string(have).substr(0, 4) + " N m";
        return;
      }
    }
    cmd_ = c;
  }
  const std::string& refusal() const { return refusal_; }
  const Command& current() const { return cmd_; }
  Phase phase() const { return phase_; }
  double thrust() const { return thrust_; }
  double last_apex() const { return last_apex_; }
  double stance_time() const { return stance_time_; }
  int hops() const { return hops_; }
  double standing_height() const { return standing_height_; }
  const Gains& gains() const { return k_; }

  // The standing posture's joint angles, the torso level: sitting as a rabbit sits, the hind foot flat under
  // the body with its toe 0.01 m ahead of the centre (so that the foot alone could carry the weight), the
  // fore paw 0.01 m ahead of the shoulder, both on the ground standing_height below the centre.
  std::array<double, JOINTS> stand_posture() const {
    std::array<double, JOINTS> a{};
    const double ground = -standing_height_;
    const std::array<double, 2> ankle = {0.01 - 0.10 - r_->hip[0], ground - r_->hip[1]};
    if (const auto ik = two_link(0.10, 0.12, ankle, +1)) a[0] = (*ik)[0], a[1] = (*ik)[1] - (*ik)[0], a[2] = -(*ik)[1];
    const std::array<double, 2> paw = {0.01, ground - r_->shoulder[1]};
    if (const auto ik = two_link(0.07, 0.08, paw, -1)) a[3] = (*ik)[0], a[4] = (*ik)[1] - (*ik)[0];
    return a;
  }
  // The crouch a hop starts from: the hind foot flat under the body as in sitting, the body lowered so that
  // the toe is ready_leg from the centre, the forelegs lifted. The sole, from heel to toe, holds the centre
  // above it, so the crouch needs no balancing on a point.
  std::array<double, JOINTS> ready_posture() const {
    std::array<double, JOINTS> a{};
    const double ground = -std::sqrt(k_.ready_leg * k_.ready_leg - 0.01 * 0.01);
    const std::array<double, 2> ankle = {0.01 - 0.10 - r_->hip[0], ground - r_->hip[1]};
    if (const auto ik = two_link(0.10, 0.12, ankle, +1)) a[0] = (*ik)[0], a[1] = (*ik)[1] - (*ik)[0], a[2] = -(*ik)[1];
    a[3] = -105 * DEG, a[4] = 95 * DEG;  // the forelegs lifted
    return a;
  }

  Coordinates standing_state(double x0, double ground) const {
    Coordinates q{};
    q[X] = x0, q[Z] = ground + standing_height_ - 0.0009, q[PITCH] = 0;
    const auto a = stand_posture();
    for (std::size_t j = 0; j < JOINTS; j++) q[HIP + j] = a[j];
    return q;
  }

  // The torques for an observation.
  std::array<double, JOINTS> step(const Observation& o) {
    const double scale = g_ / 1.62;
    const auto& tree = *tree_;
    const auto st = tree.kinematics(o.q, o.v);
    const auto com = tree.centre_of_mass(o.q);
    const auto mom = tree.linear_momentum(o.q, o.v);
    const double m = r_->mass();
    const std::array<double, 2> comv = {mom[0] / m, mom[1] / m};
    const auto toe = tree.point(st, FOOT, r_->feet[0].point), toev = tree.point_velocity(st, FOOT, r_->feet[0].point);
    const std::array<double, 2> rl = {com[0] - toe[0], com[1] - toe[1]};  // the virtual leg, toe to centre
    const double L = std::hypot(rl[0], rl[1]);
    const std::array<double, 2> u = {rl[0] / L, rl[1] / L};
    const double Ldot = u[0] * (comv[0] - toev[0]) + u[1] * (comv[1] - toev[1]);
    const bool touching = o.foot_force[0] > k_.touch * scale;
    const double theta = o.q[PITCH], omega = o.v[PITCH];

    // The phases.
    if (cmd_.behaviour == Behaviour::hop && phase_ == Phase::standing) {
      phase_ = Phase::ready, ready_start_ = o.t, ground_ = toe[1];
      for (std::size_t j = 0; j < JOINTS; j++) ready_from_[j] = o.q[HIP + j];
    }
    if (phase_ == Phase::ready && o.t - ready_start_ >= k_.ready_time) phase_ = Phase::thrust, stance_start_ = o.t, launching_ = true;
    if (cmd_.behaviour == Behaviour::stand) phase_ = Phase::standing;
    if (phase_ == Phase::flight) {
      apex_ = std::max(apex_, com[1]);
      if (touching && o.t - liftoff_ > 0.04) {
        phase_ = Phase::compression, stance_start_ = o.t;
        last_apex_ = apex_ - ground_ - standing_height_;
        if (cmd_.behaviour == Behaviour::hop) {
          thrust_ = std::clamp(thrust_ + k_.height_gain * scale * (cmd_.height - last_apex_), 0.0, k_.thrust_max * scale);
          // Raibert's placement leaves a steady error of speed; an integral of it, hop by hop, removes it.
          speed_bias_ = std::clamp(speed_bias_ + k_.speed_integral * (cmd_.speed - comv[0]), -0.3, 0.3);
        }
      }
    } else if (phase_ == Phase::compression) {
      if (Ldot >= 0 && o.t - stance_start_ > 0.01) phase_ = Phase::thrust;
    } else if (phase_ == Phase::thrust) {
      if (!touching && o.t - stance_start_ > 0.03) {
        phase_ = Phase::flight, liftoff_ = o.t, apex_ = com[1], hops_++, launching_ = false;
        stance_time_ = std::clamp(o.t - stance_start_, 0.05, 0.6);
      }
    }
    if (cmd_.behaviour == Behaviour::stop && phase_ == Phase::compression && std::hypot(comv[0], comv[1]) < 0.25 * std::sqrt(scale))
      phase_ = Phase::standing, cmd_.behaviour = Behaviour::stand;

    std::array<double, JOINTS> tau{};
    const auto pd = [&](std::size_t j, double target, double w) { return k_.kp(j, w) * (target - o.q[HIP + j]) - k_.kd(j, w) * o.v[HIP + j]; };

    if (phase_ == Phase::standing || phase_ == Phase::ready) {
      // Each joint to its posture (in the ready phase eased from where it was to the crouch), and the weight
      // carried by feed-forward: the hind foot and the paw share it by the lever rule about the centre, each
      // leg's joints give −Jᵀ f for its share f, and the tendon's and the limits' torques are cancelled.
      std::array<double, JOINTS> a = stand_posture();
      bool on_toe = false;
      if (phase_ == Phase::ready) {
        const double x = std::clamp((o.t - ready_start_) / k_.ready_time, 0.0, 1.0), e = x * x * (3 - 2 * x);
        const auto b = ready_posture();
        for (std::size_t j = 0; j < JOINTS; j++) a[j] = ready_from_[j] + e * (b[j] - ready_from_[j]);
        on_toe = e > 0.5;  // the forelegs lifted: the hind foot carries the weight
      }
      for (std::size_t j = 0; j < JOINTS; j++) tau[j] = pd(j, a[j], k_.stand_bandwidth);
      const auto passive = potential_force(*r_, o.q);
      for (std::size_t j = 0; j < JOINTS; j++) tau[j] -= passive[HIP + j];
      const auto heel = tree.point(st, FOOT, r_->feet[1].point), paw = tree.point(st, FOREARM, r_->feet[2].point);
      const std::array<double, 2> mid = {(toe[0] + heel[0]) / 2, (toe[1] + heel[1]) / 2};
      const double span = paw[0] - mid[0];
      const double share = on_toe ? 1.0 : span > 0.02 ? std::clamp((paw[0] - com[0]) / span, 0.0, 1.0) : 0.5;
      const double W = m * g_;
      const auto Jh = tree.jacobian(st, FOOT, mid), Jf = tree.jacobian(st, FOREARM, paw);
      for (std::size_t j = 0; j < 3; j++) tau[j] -= Jh[1][HIP + j] * W * share;
      // While the forelegs lift, the hip holds the torso level against the planted foot.
      if (phase_ == Phase::ready) tau[0] += k_.attitude_kp * scale * theta + k_.attitude_kd * std::sqrt(scale) * omega;
      for (std::size_t j = 3; j < 5; j++) tau[j] -= Jf[1][HIP + j] * W * (1 - share);
      ground_ = toe[1];
      return limit(tau, o);
    }

    // The forelegs lifted a little from standing, slowly: a fast swing would turn the torso.
    tau[3] = pd(3, -105 * DEG, k_.arm_bandwidth);
    tau[4] = pd(4, 95 * DEG, k_.arm_bandwidth);

    if (phase_ == Phase::flight) {
      // Place the toe for the next landing: ahead of the centre by the neutral point and the speed's error, the
      // leg at its rest length; the joints by the leg's inverse kinematics from the hip, the ankle stretched.
      const double v_target = cmd_.behaviour == Behaviour::hop ? cmd_.speed + speed_bias_ : 0.0;
      const double xf = std::clamp(comv[0] * stance_time_ / 2 + k_.placement * (comv[0] - v_target), -0.09, 0.09);
      const double zf = -std::sqrt(k_.leg_rest * k_.leg_rest - xf * xf);
      const double c = std::cos(theta), s = std::sin(theta);
      const std::array<double, 2> hipw = {c * r_->hip[0] - s * r_->hip[1], s * r_->hip[0] + c * r_->hip[1]};
      const std::array<double, 2> rel = {xf - hipw[0], zf - hipw[1]};                            // toe from hip, world axes
      const std::array<double, 2> target = {c * rel[0] + s * rel[1], -s * rel[0] + c * rel[1]};  // torso axes
      const double lx = 0.12 + 0.10 * std::cos(k_.hop_ankle), lz = 0.10 * std::sin(k_.hop_ankle);
      const double l2 = std::hypot(lx, lz), beta = std::atan2(lz, lx);
      if (const auto ik = two_link(0.10, l2, target, +1)) {
        tau[0] = pd(0, (*ik)[0], k_.swing_bandwidth);
        tau[1] = pd(1, (*ik)[1] - beta - (*ik)[0], k_.swing_bandwidth);
      }
      tau[2] = pd(2, k_.hop_ankle, k_.swing_bandwidth);
      return limit(tau, o);
    }

    // Stance: the ground's force asked along the virtual leg, through the centre of mass (no moment about it),
    // the leg a spring with a thrust while it extends; and across it a force whose moment about the centre
    // servoes the pitch, within the friction asked. The joints give −Jᵀ f at the toe; the ankle adds a
    // servo to its hopping angle, the leg's third joint.
    double F = k_.leg_stiffness * scale * (k_.leg_rest - L) - k_.leg_damping * std::sqrt(scale) * Ldot;
    if (phase_ == Phase::thrust && cmd_.behaviour == Behaviour::hop) F += thrust_;
    // The launch from the crouch: the force that does m g (ΔL + h) of work over the extension ΔL to the rest
    // length, the hop's height h asked; the leg's spring takes over from the first landing.
    if (launching_) F = L < k_.leg_rest ? m * g_ * (1 + cmd_.height / std::max(k_.leg_rest - k_.ready_leg, 0.02)) : 0.0;
    F = std::max(F, 0.0);
    // The moment about the centre asked: the pitch back to level, and the angular momentum about the centre
    // to zero, since a flight keeps what the stance leaves it (and the Moon's flights are long).
    const double Lcom = tree.angular_momentum(o.q, o.v, com);
    const double Md = -k_.attitude_kp * scale * theta - k_.momentum_rate * Lcom;
    // The point the force is asked at: the toe, or in the launch, while the heel is down, the point of the sole
    // beneath the centre (the centre of pressure that puts the force through it), and the leg from there.
    std::array<double, 2> pc = toe;
    if (launching_ && o.foot_force[1] > k_.touch * scale) {
      const auto heel = tree.point(st, FOOT, r_->feet[1].point);
      const double w = std::clamp((com[0] - heel[0]) / std::max(toe[0] - heel[0], 1e-3), 0.0, 1.0);
      pc = {heel[0] + w * (toe[0] - heel[0]), heel[1] + w * (toe[1] - heel[1])};
    }
    const std::array<double, 2> rc = {com[0] - pc[0], com[1] - pc[1]};
    const double Lc = std::hypot(rc[0], rc[1]);
    const std::array<double, 2> uc = {rc[0] / Lc, rc[1] / Lc};
    const std::array<double, 2> r = {-rc[0], -rc[1]};
    double ft = Md / (Lc * Lc);  // f_t = M (−r_z, r_x)/|r|²: its moment r × f_t is M
    // No more friction than it believes is there (three quarters of it), nor than its gains ask.
    const double mu = o.friction > 0 ? std::min(k_.friction, 0.75 * o.friction) : k_.friction;
    ft = std::clamp(ft, -mu * F / Lc, mu * F / Lc);
    const std::array<double, 2> f = {F * uc[0] - ft * r[1], F * uc[1] + ft * r[0]};
    const auto J = tree.jacobian(st, FOOT, pc);
    for (std::size_t j = 0; j < 3; j++) tau[j] = -(J[0][HIP + j] * f[0] + J[1][HIP + j] * f[1]);
    if (!launching_) tau[2] += pd(2, k_.hop_ankle, k_.swing_bandwidth);
    ground_ = toe[1];
    return limit(tau, o);
  }

 private:
  std::array<double, JOINTS> limit(std::array<double, JOINTS> tau, const Observation& o) const {
    for (std::size_t j = 0; j < JOINTS; j++) tau[j] = r_->actuators[j].available(tau[j], o.v[HIP + j]);
    return tau;
  }

  const Rabbit* r_;
  const pm::PlanarTree<DOF>* tree_ = nullptr;
  double g_;
  Gains k_;
  Command cmd_{};
  Phase phase_ = Phase::standing;
  double standing_height_, thrust_, stance_time_;
  double stance_start_ = 0, liftoff_ = -1, apex_ = 0, last_apex_ = 0, ground_ = 0, ready_start_ = 0;
  std::array<double, JOINTS> ready_from_{};
  bool launching_ = false;
  double speed_bias_ = 0;
  std::string refusal_;
  int hops_ = 0;
};

}  // namespace robotics::legged
