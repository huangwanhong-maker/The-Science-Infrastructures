# robotics

*Status: draft. Field, height 0.* Systems that sense and act: their bodies on the ground, the controllers that drive them, the plans they follow and the sensors they carry ([ADR 0054](../../development/decisions/0054-bodies-dynamics-and-knowledge.md)). It was founded at the user's request ("you may add a domain robotics, if needed"), earlier than its place in tier 4 of the order (README §4), for the embodied agents of [Cosmopoiesis](../../applications/cosmopoiesis-lab/): first a rover surveying the Moon's magnetic field. It serves use and learning.

## 1. Concepts

- **A robot is a body with a controller.** Physically a vehicle is a body like any other, and its mechanics is the foundations': constrained motion (calx's RATTLE), rigid bodies and discrete actions (physis `mechanics`). What makes it a robot is the loop:
  - its sensors read the world through mounts on its body;
  - its controller turns what it believes into commands;
  - its motors turn the commands into forces, within what the ground allows.
- **A pose on a sphere is a rotation.** $SO(3)$ is the unit tangent bundle of $S^2$, so a place and a heading are one rotation. Driving and turning is $\dot A = A\,(0, v/R, \omega)^\wedge$ ([theory §1](theory/vehicles.md)).
- **One vehicle at two fidelities.** The two models take the same commands:
  - the kinematic vehicle rolls without slipping and keeps the speed it is given, its step exact for a constant speed ([theory §2](theory/vehicles.md));
  - the dynamic vehicle is a mass held to the ground by its weight. Its wheels drive it through a speed loop against the ground's Coulomb friction, each kick solved exactly, about a RATTLE step ([theory §3](theory/vehicles.md)).

  The kinematic vehicle is the dynamic one's limit. Each state translates into the other exactly, and the dynamic vehicle lags by its loop's $v^*\tau$, less a step.
- **Precision far from the centre.** On the Moon a place is $1.7 \times 10^6$ m from the sphere's centre. Both vehicles keep it as a local origin and a displacement, so a long drive gathers no rounding of $R$.
- **A controller sees beliefs, not the world.** The waypoint follower is given a pose and returns a command, whatever gave the pose ([theory §4](theory/vehicles.md)).
- **A survey is planned by a pattern or by what is unknown.** A lawnmower carried onto the sphere by the exponential map, or the greedy choice of the candidate of greatest value less cost. With the value a posterior variance, this is uncertainty sampling ([theory §5](theory/vehicles.md)).

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `vehicles` | A vehicle on a sphere: its pose as a rotation (`surface.hpp`), the kinematic vehicle (`kinematic.hpp`), the dynamic vehicle with its wheels, speed loop and the ground's friction (`dynamic.hpp`), the translation between them ([`vehicles/`](vehicles/include/robotics/vehicles/)) | draft |
| `control` | Following waypoints on a sphere: the heading error, a proportional heading controller, arrival ([`control/`](control/include/robotics/control/)) | draft |
| `planning` | Surveys: the lawnmower on the sphere, rings of candidates, the best candidate by value less cost ([`planning/`](planning/include/robotics/planning/)) | draft |
| `sensing` | Mounts (a sensor's place and axes on its vehicle) and white noise drawn by index ([`sensing/`](sensing/include/robotics/sensing/)) | draft |
| `estimation` | The vehicle's belief of its own pose: odometry, sun and star sensors, filters | planned |
| `terrain` | Vehicles on terrain heights, with slopes; wheels on soil (Bekker–Wong) | planned |
| `manipulation` | Arms and their kinematic chains, on physis's rigid bodies | planned |

## 3. Conventions

- **Axes.** A vehicle's axes are x forward, y left, z up (ROS's REP 103), in the frame of the sphere it moves on (for the Moon, its mean-Earth frame). Yaw is positive to the left, about up. An azimuth is measured from north, clockwise towards east.
- **Units.** SI (CLAUDE.md §4.1): metres, seconds, kilograms, newtons, watts, radians.
- **The sphere.** A radius and a GM, the caller's. A terrain's heights are not yet carried.

## 4. Dependencies

The foundations alone (`module.json`): physis (`geometry`: rotations) and calx (`integrators`: RATTLE; `random`: counter-based streams). It names no body of the Solar System: the Moon's radius and GM come from those who use it.

## 5. Checks

`node domains/robotics/tests/run-tests.mjs` builds the C++ checks with CMake and runs them on calx's verification. On the Moon's sphere:

- `vehicles`:
  - the kinematic vehicle driving 2 km straight on its great circle;
  - its turning circle at the angle $\arctan(v/(R\omega))$ to its fixed axis, closing after $2\pi/|w|$;
  - its speed slewing, and the odometer both ways;
  - the dynamic vehicle on the sphere through 3 km and two changes of origin, its momentum tangent and its normal force $m(g - v^2/R)$;
  - coasting and sideslip stopping in $v_0^2/(2\mu g)$;
  - its speed loop exactly the discrete loop, and the traction limit $(\mu - C_{rr})g$;
  - the translation exact, and the lag $v^*(\tau - h)$ to $10^{-10}$ m.
- `control`: the heading error's sign and size; both vehicles reaching a waypoint behind them after turning in place; both driving a lawnmower to its end.
- `planning`: the lawnmower's lines, their spacing and directions; a ring's radius and order; the best candidate.
- `sensing`: a mast's height; a turned sensor's axes; noise drawn by index, its mean, deviation and independence.

## 6. References

See [theory/vehicles.md](theory/vehicles.md). Chief among them:

- Hairer, E., Lubich, C. and Wanner, G. (2006), *Geometric Numerical Integration*, 2nd ed., Springer;
- Wong, J. Y. (2008), *Theory of Ground Vehicles*, 4th ed., Wiley;
- Coulter, R. C. (1992), *Implementation of the Pure Pursuit Path Tracking Algorithm*, CMU-RI-TR-92-01;
- Shewry, M. C. and Wynn, H. P. (1987), Maximum entropy sampling, *J. Appl. Statist.* 14, 165–170.

## 7. Open Questions

- **Terrain and soil.** Heights from the data library's tiles, slopes the dynamic vehicle feels, and terramechanics in place of Coulomb friction ([theory §7](theory/vehicles.md)).
- **The chassis as a rigid body**, with suspension and contact, on physis's Lie group integrator.
- **Estimation.** The controller is given the true pose. An estimator would give it a belief, with tyche's planned `assimilation`.
- **Where the agent loop lives.** The loop of sensing, believing, deciding and acting is written in Cosmopoiesis's engine for its first rover. It moves here when a second work needs it.
