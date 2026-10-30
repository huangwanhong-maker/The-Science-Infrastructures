// A joint's geared motor, as a robot's joints are driven (ADR 0058; theory/legged.md §2).
//
// The standard model of a DC motor through a gear (Siciliano, Sciavicco, Villani and Oriolo 2009, Robotics,
// §5.2): at joint speed ω it gives at most τ_s (1 − ω/ω₀) forward and −τ_s (1 + ω/ω₀) backward, τ_s the
// stall torque and ω₀ the no-load speed, both capped at the driver's peak τ_p (τ_s = 2.5 τ_p here); its
// rotor's inertia through the gear, G² J_r, is an armature on the joint's coordinate; its gears lose a
// viscous b ω and a Coulomb τ_c tanh(ω/ω_c); it draws the electrical power τω + (τ/k_τ)² R, the mechanical
// power and the copper loss. A command is clipped to what is available, so near the no-load speed a motor
// cannot drive, and against its motion it brakes at its full peak. The rabbit's joints and the fish's are
// driven so.

#pragma once

#include <algorithm>

namespace robotics::actuation {

struct GearedMotor {
  double peak = 0;         // N m, the driver's current limit at the joint
  double no_load = 0;      // rad/s
  double armature = 0;     // kg m²
  double viscous = 0;      // N m s
  double coulomb = 0;      // N m
  double torque_constant = 0;  // N m/A at the joint
  double resistance = 0;   // Ω
  double stall() const { return 2.5 * peak; }
  // The torque the motor can give at speed ω, the command clipped to its torque–speed line.
  double available(double command, double omega) const {
    const double hi = std::min(peak, stall() * (1 - omega / no_load)), lo = std::max(-peak, -stall() * (1 + omega / no_load));
    return std::clamp(command, std::min(lo, hi), std::max(lo, hi));
  }
  double electrical_power(double torque, double omega) const {
    const double i = torque / torque_constant;
    return torque * omega + i * i * resistance;
  }
};

}  // namespace robotics::actuation
