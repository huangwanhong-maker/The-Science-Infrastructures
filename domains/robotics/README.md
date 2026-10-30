# robotics

*Status: draft. Field, height 0.* Systems that sense and act: their bodies on the ground, the controllers that drive them, the plans they follow and the sensors they carry ([ADR 0054](../../development/decisions/0054-bodies-dynamics-and-knowledge.md)). It was founded at the user's request ("you may add a domain robotics, if needed"), earlier than its place in tier 4 of the order (README §4), for the embodied agents of [Cosmopoiesis](../../applications/cosmopoiesis-lab/): first a rover surveying the Moon's magnetic field, then a hopping rabbit with joints ([ADR 0058](../../development/decisions/0058-a-body-with-joints.md)). It serves use and learning.

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
- **A body with joints is a kinematic tree.** The rabbit is six rigid links on a floating base, its configuration space $SE(2)\times T^5$ in its plane. Its Lagrangian and its Newton–Euler equations are physis's, one motion computed two ways. What is the domain's: its links as built, its motors with their torque–speed lines, its tendon, its feet, and its controller ([theory: legged](theory/legged.md)).
- **One motion, three methods.** The variational integrator of the discrete action, semi-implicit Euler and Runge–Kutta move the same body by the same forces. They agree on the motion to millimetres and differ in what they keep (the energy, a flight's angular momentum) and what they cost ([theory: legged §5](theory/legged.md)).
- **A body learns its ground.** The robot is told what the world presumes of its ground (a median and a spread for each parameter of the contact). Its feet read the truth: a sliding foot feels the friction, a foot at its deepest the stiffness. Its belief is the presumption conditioned on those readings, and the controller acts on it ([theory: legged §7](theory/legged.md), ADR 0059).
- **A fish swims by its strips.** The water acts on each strip of its body: the kinetic energy of the water it moves is in its Lagrangian, its drag a force. It swims and turns by a travelling wave. It swims by drag, at a beat several times a real fish's, since the model leaves out the momentum the tail sheds into its wake ([theory: swimming](theory/swimming.md)).
- **A dog trots on diagonal pairs.** In its plane its four legs are separate, so a fore foot and a hind foot are always down. It distributes the force its body needs between them, within the friction it believes its ground gives ([theory: quadruped](theory/quadruped.md)).
- **A controller that refuses.** A command the body cannot do (a hop the ankle cannot launch, on the Earth) is refused with its reason, not attempted ([theory: legged §4.5](theory/legged.md)).
- **A survey is planned by a pattern or by what is unknown.** A lawnmower carried onto the sphere by the exponential map, or the greedy choice of the candidate of greatest value less cost. With the value a posterior variance, this is uncertainty sampling ([theory §5](theory/vehicles.md)).

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `vehicles` | A vehicle on a sphere: its pose as a rotation (`surface.hpp`), the kinematic vehicle (`kinematic.hpp`), the dynamic vehicle with its wheels, speed loop and the ground's friction (`dynamic.hpp`), the translation between them ([`vehicles/`](vehicles/include/robotics/vehicles/)) | draft |
| `control` | Following waypoints on a sphere: the heading error, a proportional heading controller, arrival ([`control/`](control/include/robotics/control/)) | draft |
| `planning` | Surveys: the lawnmower on the sphere, rings of candidates, the best candidate by value less cost ([`planning/`](planning/include/robotics/planning/)) | draft |
| `sensing` | Mounts (a sensor's place and axes on its vehicle) and white noise drawn by index ([`sensing/`](sensing/include/robotics/sensing/)) | draft |
| `legged` | A rabbit: a hopping robot of 1.69 kg in its plane, its links, motors, tendon, limits and feet (`rabbit.hpp`); its controller, standing, hopping in place and forward by Raibert's parts on a virtual leg from the toe to the centre, stopping, and refusing what it cannot do (`controller.hpp`); its simulation by three methods of one motion, with the energy's account (`simulation.hpp`); the ground along its line, its feet as sensors and its belief of friction and load conditioned on their readings (`ground.hpp`, ADR 0059) ([`legged/`](legged/include/robotics/legged/), [theory](theory/legged.md)) | draft |
| `estimation` | The vehicle's belief of its own pose: odometry, sun and star sensors, filters | planned |
| `terrain` | Vehicles on terrain heights, with slopes; wheels on soil (Bekker–Wong) | planned |
| `manipulation` | Arms and their kinematic chains, on physis's rigid bodies | planned |
| `swimming` | A fish of 1.85 kg in its horizontal plane, its head and four driven links in water (physis's strip theory), its servos and its travelling wave (swimming, turning, stopping), its simulation by three methods (`fish.hpp`, ADR 0059) ([`swimming/`](swimming/include/robotics/swimming/), [theory](theory/swimming.md)) | draft |
| `quadruped` | A dog of 10.2 kg in its sagittal plane, its four legs separate so that it trots on diagonal pairs; standing, trotting in place and forward, stopping; its stance pair's force distributed by the moment's balance within the friction it believes, its swing pair to Raibert's placement; three methods; its feet learning its ground (`dog.hpp`, ADR 0059) ([`quadruped/`](quadruped/include/robotics/quadruped/), [theory](theory/quadruped.md)) | draft |
| `actuation` | A joint's geared motor: its torque–speed line, armature, friction and copper loss, the rabbit's and the fish's (`motor.hpp`) ([`actuation/`](actuation/include/robotics/actuation/)) | draft |

## 3. Conventions

- **Axes.** A vehicle's axes are x forward, y left, z up (ROS's REP 103), in the frame of the sphere it moves on (for the Moon, its mean-Earth frame). Yaw is positive to the left, about up. An azimuth is measured from north, clockwise towards east.
- **A planar robot** moves in its plane $(x, z)$, the vehicle's axes without $y$. Seen from its right side, a positive angle turns $x$ toward $z$; a joint's positive torque turns its link so.
- **Units.** SI (CLAUDE.md §4.1): metres, seconds, kilograms, newtons, watts, radians.
- **The sphere.** A radius and a GM, the caller's. A terrain's heights are not yet carried.

## 4. Dependencies

The foundations alone (`module.json`): tyche (`spatial`: Gaussian conditioning), physis (`geometry`: rotations; `mechanics`: planar trees and compliant contact; `analysis`: stationary random fields, in the checks) and calx (`integrators`: RATTLE and Runge–Kutta; `random`: counter-based streams; `differentiation`: jets; `linalg`: small dense systems). It names no body of the Solar System: the Moon's radius and GM come from those who use it.

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
- `legged` (36), on the Moon's gravity:
  - the body: its mass, the legs' inverse kinematics, the standing posture on its feet, a motor's torque–speed line and its electrical power, the potential's force against differences;
  - standing still and level, the feet carrying the weight, the energy accounted for;
  - hopping in place (the apex brought to the one asked) and forward at the speed asked, the pitch held within 15°, stopping to sit;
  - the three methods agreeing on the first apex, the variational integrator's energy best kept, Runge–Kutta's flight momentum to $10^{-8}$;
  - on the Earth the hop refused with its reason;
  - the ground it learns: calibrated on a ground drawn from its presumption; across a slippery patch it keeps hopping, reads the friction there, and its belief says the presumption was wrong.
- `swimming` (15): the fish's body; neutrally buoyant, still; swimming at the speed asked on its course, its wake pushing the water backward; turning; gliding to rest; the three methods' distances; carried by a current.
- `actuation` (6): a motor's torque–speed line, braking and electrical power.
- `quadruped` (18): the dog's body and legs; standing on rock; trotting in place and at the speed asked, level and at its height; stopping to stand; the three methods' distances; on ice, gathering speed more gently.

## 6. References

See [theory/vehicles.md](theory/vehicles.md) and [theory/legged.md](theory/legged.md). Chief among them:

- Hairer, E., Lubich, C. and Wanner, G. (2006), *Geometric Numerical Integration*, 2nd ed., Springer;
- Raibert, M. H. (1986), *Legged Robots That Balance*, MIT Press;
- Featherstone, R. (2008), *Rigid Body Dynamics Algorithms*, Springer;
- Wong, J. Y. (2008), *Theory of Ground Vehicles*, 4th ed., Wiley;
- Coulter, R. C. (1992), *Implementation of the Pure Pursuit Path Tracking Algorithm*, CMU-RI-TR-92-01;
- Shewry, M. C. and Wynn, H. P. (1987), Maximum entropy sampling, *J. Appl. Statist.* 14, 165–170.

## 7. Open Questions

- **Terrain and soil.** Heights from the data library's tiles, slopes the dynamic vehicle feels, and terramechanics in place of Coulomb friction ([theory §7](theory/vehicles.md)).
- **The chassis as a rigid body**, with suspension and contact, on physis's Lie group integrator.
- **Estimation.** The controller is given the true pose. An estimator would give it a belief, with tyche's planned `assimilation`.
- **A reduced model and a typed transition.** The rabbit reduced to a spring-loaded inverted pendulum, and the translation between the two with what it keeps ([theory: legged §8](theory/legged.md)).
- **Where the agent loop lives.** The loop of sensing, believing, deciding and acting is written in Cosmopoiesis's engine for its first rover. It moves here when a second work needs it.
