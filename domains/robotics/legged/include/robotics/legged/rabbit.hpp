// A rabbit: a small hopping robot's body, its motors, tendon and feet, in its own vertical plane (ADR 0058;
// theory/legged.md §1–2).
//
// The body. A rabbit-sized robot of 1.69 kg, its proportions a European rabbit's (Oryctolagus cuniculus, of
// 1.5–2.5 kg: a long hind leg of thigh, shank and a long foot, a short foreleg). It moves in its sagittal
// plane: the two hind legs move together and are one leg of their summed mass, as are the forelegs. The
// links and their masses are a build's, not an animal's: the torso carries the battery, the electronics
// and the hip and knee motors (the knee driven from the hip by a belt, so that the leg is light), and each
// leg link's mass is its structure's and its distal motor's. The battery is two packs at the torso's ends
// (0.4 kg each, 0.12 m from the centre), which raises its inertia from a uniform box's 1.08 × 10⁻² to
// 1.6 × 10⁻² kg m², so that the legs' swing in flight turns it less.
//
//   link     parent  joint (parent frame)  length  mass     centre    inertia (kg m²)
//   torso    —       —                      0.30 m  1.30 kg  (0, 0)    1.6 × 10⁻²    box 0.30 × 0.10 m, the head in it
//   thigh    torso   hip (−0.11, −0.03)     0.10    0.16     0.045     1.33 × 10⁻⁴
//   shank    thigh   knee (0.10, 0)         0.12    0.08     0.050     9.6 × 10⁻⁵
//   foot     shank   ankle (0.12, 0)        0.10    0.05     0.040     4.2 × 10⁻⁵    toe (0.10, 0), heel (−0.01, 0)
//   arm      torso   shoulder (0.10, −0.035) 0.07   0.06     0.030     2.45 × 10⁻⁵
//   forearm  arm     elbow (0.07, 0)        0.08    0.04     0.040     2.1 × 10⁻⁵    paw (0.08, 0)
//
// A link's frame points from its joint to its far end; its angle is its parent's plus its joint's
// (physis multibody.hpp). The coordinates are q = (x, z, θ, q_hip, q_knee, q_ankle, q_shoulder,
// q_elbow): the torso's centre of mass, its pitch (nose up positive), and the joints.
//
// The motors (one actuator stands for the pair of a joint's two sides). Each is a brushless motor and a
// gear: at the joint, a peak torque the driver's current limit allows, a no-load speed the supply's voltage
// allows, and the torque–speed line between them: the torque at speed ω lies in
// [max(−τ_peak, −τ_stall (1 + ω/ω₀)), min(τ_peak, τ_stall (1 − ω/ω₀))], τ_stall = 2.5 τ_peak. Its rotor adds
// an armature G² J_r to its coordinate; its gear loses a viscous b ω and a Coulomb τ_c tanh(ω/ω_s). Its
// electrical power is τ ω + (τ/k_τ)² R, the mechanical and the copper loss (k_τ the torque constant
// through the gear, R the winding's resistance); a negative power is braking, which the battery is not
// assumed to absorb (the supply's diode).
//
//   joint     τ_peak    ω₀         armature        b                τ_c     k_τ          R
//   hip       5.0 N m   30 rad/s   8 × 10⁻⁴ kg m²  2 × 10⁻³ N m s   0.02    0.45 N m/A   0.35 Ω
//   knee      5.0       30         8 × 10⁻⁴        2 × 10⁻³         0.02    0.45         0.35
//   ankle     2.5       40         3 × 10⁻⁴        1 × 10⁻³         0.01    0.25         0.45
//   shoulder  2.0       40         2 × 10⁻⁴        1 × 10⁻³         0.01    0.20         0.45
//   elbow     2.0       40         2 × 10⁻⁴        1 × 10⁻³         0.01    0.20         0.45
//
// The tendon. A spring across the ankle, as the Achilles tendon is, of 1.5 N m/rad about the foot stretched
// as in a hop (75°): slack in flight, so the motor need not hold the foot there; stretched as the ankle bends
// on landing, taking the landing's energy and returning it in the push-off; and carrying part of the weight
// when the rabbit sits, its foot flat. It is part of the potential, so the action sees it.
//
// The joints' limits. Each joint stops softly at its range: a stiff spring (20 N m/rad) beyond it, in the
// potential, and a damper (0.05 N m s/rad) while it is beyond.
//
// The feet. Three points touch the ground: the hind toe, the hind heel and the fore paw (and, if it falls, the
// torso's four corners and its head), each by physis's
// compliant contact (contact.hpp) with a footpad on compacted regolith or soil: k = 3 × 10⁵ N/m^{3/2},
// c = 1.5 s/m, μ = 0.7 (regolith's internal friction is 30–50°, tan φ ≈ 0.6–1.2), v_s = 0.01 m/s.
//
// What it is not: a rabbit's spine, which bends in a bound; its two sides apart (no roll or yaw, no turning);
// a soil that sinks; a motor's electrical dynamics (its current follows the command within the step).

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

#include <calx/differentiation/jet.hpp>
#include <physis/mechanics/contact.hpp>
#include <physis/mechanics/multibody.hpp>

#include <robotics/actuation/motor.hpp>
#include <robotics/legged/ground.hpp>

namespace robotics::legged {

namespace pm = physis::mechanics;

inline constexpr std::size_t DOF = 8;
using Coordinates = pm::Coordinates<DOF>;
inline constexpr double DEG = std::numbers::pi / 180;

// The coordinates.
enum Coordinate : std::size_t { X = 0, Z = 1, PITCH = 2, HIP = 3, KNEE = 4, ANKLE = 5, SHOULDER = 6, ELBOW = 7 };
enum Link : std::size_t { TORSO = 0, THIGH = 1, SHANK = 2, FOOT = 3, ARM = 4, FOREARM = 5 };
inline constexpr std::size_t JOINTS = 5;  // HIP … ELBOW
inline constexpr std::array<const char*, JOINTS> JOINT_NAMES = {"hip", "knee", "ankle", "shoulder", "elbow"};

// A joint's geared motor: the domain's (actuation/motor.hpp).
using Actuator = robotics::actuation::GearedMotor;

struct JointLimit {
  double lower = 0, upper = 0;  // rad
};

// A point that touches the ground, on a link.
struct Foot {
  const char* name;
  Link link;
  std::array<double, 2> point;
};

// How a link is drawn: a box (a length along its frame's x from `from`, a height across) or a capsule from
// its joint along x; and how wide it is across the plane, for the pair of sides it stands for.
struct Shape {
  const char* kind;   // "box" | "capsule" | "sphere"
  std::array<double, 2> from;
  double length, thickness, width;
  bool paired;        // drawn twice, at ±width apart (a left and a right leg)
};

struct Rabbit {
  std::vector<pm::PlanarLink> links;
  std::array<Actuator, JOINTS> actuators;
  std::array<JointLimit, JOINTS> limits;
  double tendon_stiffness = 1.5;     // N m/rad
  double tendon_rest = 0;            // rad, the ankle's angle at which the tendon is slack
  double limit_stiffness = 20;       // N m/rad
  double limit_damping = 0.05;       // N m s/rad
  double coulomb_speed = 0.05;       // rad/s, the gears' friction regularised below it
  std::vector<Foot> feet;
  pm::CompliantContact ground_law;
  std::vector<Shape> shapes;         // one a link
  std::array<double, 2> hip{}, shoulder{};

  double mass() const {
    double m = 0;
    for (const auto& l : links) m += l.mass;
    return m;
  }
};

inline Rabbit rabbit() {
  Rabbit r;
  r.hip = {-0.11, -0.03};
  r.shoulder = {0.10, -0.035};
  const auto rod = [](double m, double l) { return m * l * l / 12; };
  r.links = {
      {"torso", -1, {0, 0}, 0.0, 1.30, {0.0, 0.0}, 1.6e-2, 0},
      {"thigh", 0, r.hip, 0.0, 0.16, {0.045, 0.0}, rod(0.16, 0.10), 8e-4},
      {"shank", 1, {0.10, 0.0}, 0.0, 0.08, {0.050, 0.0}, rod(0.08, 0.12), 8e-4},
      {"foot", 2, {0.12, 0.0}, 0.0, 0.05, {0.040, 0.0}, rod(0.05, 0.10), 3e-4},
      {"arm", 0, r.shoulder, 0.0, 0.06, {0.030, 0.0}, rod(0.06, 0.07), 2e-4},
      {"forearm", 4, {0.07, 0.0}, 0.0, 0.04, {0.040, 0.0}, rod(0.04, 0.08), 2e-4},
  };
  r.actuators = {{{5.0, 30, 8e-4, 2e-3, 0.02, 0.45, 0.35},
                  {5.0, 30, 8e-4, 2e-3, 0.02, 0.45, 0.35},
                  {2.5, 40, 3e-4, 1e-3, 0.01, 0.25, 0.45},
                  {2.0, 40, 2e-4, 1e-3, 0.01, 0.20, 0.45},
                  {2.0, 40, 2e-4, 1e-3, 0.01, 0.20, 0.45}}};
  // The joints' ranges: the hip from the thigh raised forward to stretched back; the knee folded (the shank
  // turned back under the thigh); the ankle from the foot folded to stretched; the shoulder and elbow.
  r.limits = {{{-150 * DEG, 30 * DEG}, {-165 * DEG, -15 * DEG}, {20 * DEG, 175 * DEG}, {-175 * DEG, -30 * DEG}, {-10 * DEG, 150 * DEG}}};
  r.tendon_rest = 75 * DEG;
  r.feet = {{"toe", FOOT, {0.10, 0.0}}, {"heel", FOOT, {-0.01, 0.0}}, {"paw", FOREARM, {0.08, 0.0}},
            {"rump", TORSO, {-0.15, -0.05}}, {"chest", TORSO, {0.15, -0.05}}, {"back", TORSO, {-0.15, 0.05}}, {"head", TORSO, {0.17, 0.06}}};
  r.ground_law.stiffness = 3e5, r.ground_law.exponent = 1.5, r.ground_law.damping = 1.5, r.ground_law.friction = 0.7, r.ground_law.slip_speed = 0.01;
  r.shapes = {{"box", {-0.15, 0.0}, 0.30, 0.10, 0.12, false},
              {"capsule", {0, 0}, 0.10, 0.030, 0.05, true},
              {"capsule", {0, 0}, 0.12, 0.022, 0.05, true},
              {"capsule", {-0.01, 0}, 0.11, 0.018, 0.05, true},
              {"capsule", {0, 0}, 0.07, 0.020, 0.07, true},
              {"capsule", {0, 0}, 0.08, 0.016, 0.07, true}};
  return r;
}

// The potential beyond gravity: the tendon and the joints' limits, over any number type.
struct RabbitPotential {
  const Rabbit* r;
  template <class S, std::size_t N>
  S operator()(const std::array<S, N>& q) const {
    S V = S(0.0);
    const S a = q[ANKLE] - r->tendon_rest;
    V += a * a * (0.5 * r->tendon_stiffness);
    for (std::size_t j = 0; j < JOINTS; j++) {
      const S& x = q[HIP + j];
      if (x > r->limits[j].upper) {
        const S d = x - r->limits[j].upper;
        V += d * d * (0.5 * r->limit_stiffness);
      } else if (x < r->limits[j].lower) {
        const S d = r->limits[j].lower - x;
        V += d * d * (0.5 * r->limit_stiffness);
      }
    }
    return V;
  }
};

// The generalized force of that potential, −∂V/∂q, for the explicit methods, which take every force but
// inertia and gravity as a force (the variational integrator has it in its Lagrangian).
inline Coordinates potential_force(const Rabbit& r, const Coordinates& q) {
  const auto z = calx::differentiation::variables<DOF, 1>(q);
  const auto V = RabbitPotential{&r}(z);
  Coordinates f;
  for (std::size_t k = 0; k < DOF; k++) f[k] = -V.d[k];
  return f;
}

// What the ground does at each foot, on doubles: for the readings and the energy's account.
struct FootContact {
  std::array<double, 2> place{}, force{};
  double normal = 0, penetration = 0, slip = 0, tangential = 0;
  bool touching = false;
};

// The generalized forces of everything but inertia, gravity and the potential: the motors' torques held
// over the step, the gears' friction, the limits' damping and the ground at each foot.
struct RabbitForces {
  const pm::PlanarTree<DOF>* tree;
  const Rabbit* r;
  pm::GroundLine ground;
  const std::array<double, JOINTS>* torque = nullptr;  // the motors' torques, held between the controller's samples
  TabulatedGround* table = nullptr;  // the ground's law along the line, where it varies; else the rabbit's own

  pm::CompliantContact law_at(double x) const { return table ? table->at(x) : r->ground_law; }

  template <class S>
  std::array<S, DOF> operator()(double, const std::array<S, DOF>& q, const std::array<S, DOF>& v) const {
    using std::tanh;
    std::array<S, DOF> Q;
    Q.fill(S(0.0));
    for (std::size_t j = 0; j < JOINTS; j++) {
      const auto& m = r->actuators[j];
      const S& w = v[HIP + j];
      Q[HIP + j] += (torque ? (*torque)[j] : 0.0) - w * m.viscous - tanh(w / r->coulomb_speed) * m.coulomb;
      const S& x = q[HIP + j];
      if (x > r->limits[j].upper || x < r->limits[j].lower) Q[HIP + j] -= w * r->limit_damping;
    }
    const auto st = tree->kinematics(q, v);
    for (const auto& f : r->feet) {
      const auto p = tree->point(st, f.link, f.point);
      const auto c = pm::contact_force(ground, law_at(calx::differentiation::value(p[0])), p, tree->point_velocity(st, f.link, f.point));
      if (!c.touching) continue;
      const auto g = tree->generalized_force(st, f.link, p, c.force);
      for (std::size_t k = 0; k < DOF; k++) Q[k] += g[k];
    }
    return Q;
  }

  std::vector<FootContact> contacts(const Coordinates& q, const Coordinates& v) const {
    const auto st = tree->kinematics(q, v);
    std::vector<FootContact> out;
    const auto t = ground.tangent();
    for (const auto& f : r->feet) {
      FootContact fc;
      fc.place = tree->point(st, f.link, f.point);
      const auto pv = tree->point_velocity(st, f.link, f.point);
      const auto c = pm::contact_force(ground, law_at(fc.place[0]), fc.place, pv);
      fc.force = c.force, fc.normal = c.normal, fc.penetration = c.penetration, fc.touching = c.touching;
      fc.slip = pv[0] * t[0] + pv[1] * t[1];
      fc.tangential = c.force[0] * t[0] + c.force[1] * t[1];
      out.push_back(fc);
    }
    return out;
  }
};

}  // namespace robotics::legged
