# A Body with Joints: the Rabbit, Its Controller, and Three Methods of One Motion

*Status: draft (2026-10-30). Serves use, learning and research.* The theory of the module `legged` ([ADR 0058](../../../development/decisions/0058-a-body-with-joints.md)): `rabbit.hpp`, `controller.hpp` and `simulation.hpp`. Its mechanics is physis's planar trees and compliant contact ([physis decision 0014](https://github.com/huangwanhong-maker/physis/blob/main/development/decisions/0014-articulated-bodies-in-the-plane-and-compliant-contact.md), its theory note §6–7). What is here is what names a robot: its body as built, its motors, its feet, its controller, and the comparison of the methods that move it.

## 1. The Body as Built

**1.1 A kinematic tree.** The rabbit is a planar tree of six rigid links on a floating base. Its configuration space is $SE(2)\times T^5$, with coordinates
$$q = (x, z, \theta;\ q_{\text{hip}}, q_{\text{knee}}, q_{\text{ankle}}, q_{\text{shoulder}}, q_{\text{elbow}}).$$
The torso is the base. The thigh, shank and foot hang from the hip, and the arm and forearm from the shoulder:

| Link | Parent | Joint in the parent's frame (m) | Length (m) | Mass (kg) | Inertia about its centre (kg m²) |
|---|---|---|---|---|---|
| torso | — | — | 0.30 × 0.10 × 0.12 (a box) | 1.30 | $1.6\times10^{-2}$ |
| thigh | torso | (−0.11, −0.03) | 0.10 | 0.16 | $ml^2/12$ |
| shank | thigh | (0.10, 0) | 0.12 | 0.08 | $ml^2/12$ |
| foot | shank | (0.12, 0) | 0.10 | 0.05 | $ml^2/12$ |
| arm | torso | (0.10, −0.035) | 0.07 | 0.06 | $ml^2/12$ |
| forearm | arm | (0.07, 0) | 0.08 | 0.04 | $ml^2/12$ |

The robot weighs 1.69 kg in all. The torso's inertia is that of a box with its battery packs at the ends ($I > m(a^2 + b^2)/12 = 1.1\times10^{-2}$). This slows its pitch for a given moment, which a hopping body needs. The two legs of each pair move together: the model is planar, and a pair is one link of twice the width.

**1.2 Conventions.** The plane is $(x, z)$, $x$ forward and $z$ up. Seen from the robot's right side, as the page draws it, a positive angle turns $x$ toward $z$ (counter-clockwise). These are the domain's axes (README §3) restricted to the plane: $y$ is to the left.

**1.3 The potential beyond gravity.** Two elastic parts are written into the Lagrangian (`RabbitPotential`), and so are integrated with its structure:
- an *Achilles tendon* at the ankle, $V = \tfrac12 k_t(q_{\text{ankle}} - q_0)^2$, with $k_t = 1.5$ N m/rad and $q_0 = 75°$. A tendon stores a landing's energy and returns it, as an animal's does (Alexander 2003, ch. 2);
- the joints' *limits*, one-sided springs of 20 N m/rad beyond each range, with 0.05 N m s/rad of damping as a force.

The explicit methods take $-\partial V/\partial q$ as a force, computed by calx's jets from the same function.

## 2. The Motors

Each joint is driven by a geared motor (`Actuator`), in the standard model of a DC motor through a gear (Siciliano, Sciavicco, Villani and Oriolo 2009, §5.2):
- **The torque–speed line.** At joint speed $\omega$ the motor can give at most $\tau_s(1 - \omega/\omega_0)$ forward and $-\tau_s(1 + \omega/\omega_0)$ backward, with $\tau_s$ the stall torque and $\omega_0$ the no-load speed. A driver's current limit caps both at the peak $\tau_p$. The model takes $\tau_s = 2.5\,\tau_p$. A command is clipped to what is available, so near the no-load speed a motor cannot drive, and against its motion it brakes at its full peak.
- **Armature.** $G^2J_r$ on its coordinate, in the tree's kinetic energy (physis's theory §6, Remark).
- **Friction in the gears.** Viscous, $b\,\omega$, and Coulomb, $\tau_c\tanh(\omega/\omega_c)$.
- **Electrical power.** $P = \tau\omega + (\tau/k_\tau)^2R$, the mechanical power and the copper loss, with $k_\tau$ the torque constant at the joint and $R$ the winding's resistance. Energy drawn is counted where $P > 0$; nothing is recovered when braking.

| Joint | Peak (N m) | No-load (rad/s) | Armature (kg m²) | $k_\tau$ (N m/A) | $R$ (Ω) | Range |
|---|---|---|---|---|---|---|
| hip, knee | 5.0 | 30 | $8\times10^{-4}$ | 0.45 | 0.35 | −150° to 30°; −165° to −15° |
| ankle | 2.5 | 40 | $3\times10^{-4}$ | 0.25 | 0.45 | 20° to 175° |
| shoulder, elbow | 2.0 | 40 | $2\times10^{-4}$ | 0.20 | 0.45 | −175° to −30°; −10° to 150° |

*Caveat.* These are the values of a plausible build with small geared brushless motors, not a catalogue's. A built robot would replace them with its datasheets' and measure its armature and friction.

**Sampling.** The controller runs at 500 Hz, and its torques are held between its samples (a zero-order hold), as a robot's are.

## 3. The Feet and the Ground

Seven points touch the ground: the toe, the heel, the fore paw, and four corners and the head of the torso, which matter only in a fall. Each feels physis's compliant contact (its theory §7). The parameters:
- stiffness $k = 3\times10^5$ N/m$^{3/2}$, with Hertz's exponent $3/2$;
- damping $c = 1.5$ s/m, a restitution of about $0.9$ at 0.1 m/s (physis's Proposition 7.2);
- friction $\mu = 0.7$, regularised at $v_s = 0.01$ m/s.

Standing, each foot is pressed in by about 0.9 mm. The ground is a level line at the terrain's height under the start.

## 4. The Controller

The controller reads what a robot reads: its joints' angles and rates (encoders), its torso's pitch, rate and velocity, and its angular momentum about its centre (an inertial unit and a state estimator, here given the true state), and whether each foot presses the ground (a foot switch at 0.3 N). It returns a torque for each motor, which the motor clips (§2).

**4.1 Gains by bandwidth.** A joint follows an angle by $\tau = k_p(q^* - q) - k_d\dot q$, with
$$k_p = I\omega_b^2, \qquad k_d = 2\zeta I\omega_b ,$$
$I$ the inertia the joint moves (its armature and the links beyond it, about it), $\omega_b$ the bandwidth asked and $\zeta = 0.9$. A loop with the torque held over $T = 2$ ms loses stability as $\omega_bT$ grows. Gains equal for all joints drove the light ones (the elbow's $I = 3\times10^{-4}$) far past it, and the standing robot chattered and gained energy. Gains by inertia keep $2\zeta\omega_bT \le 0.22$ for every joint.

**4.2 Standing.** The feet are placed (the hind foot flat under the body, its toe 1 cm ahead of the centre; the paw under the shoulder). The joints' angles follow from the legs' inverse kinematics: a two-link leg reaching a point, its knee on a side. The weight is carried by feed-forward: the hind foot and the paw share $mg$ by the lever rule about the centre, each leg's joints give $-J^{\mathsf T}f$ for its share $f$, and the tendon's and the limits' torques are cancelled.

**4.3 Hopping.** Raibert's three parts (Raibert 1986, ch. 2–3), adapted to a body whose hip is far from its centre:
- **The virtual leg runs from the toe to the centre of mass**, not from the hip. A force along it passes through the centre and turns nothing. A first version had the leg from the hip, 11 cm behind the centre, and every push pitched the rabbit over.
- **Height.** In stance the leg is a spring, $F = k(L_0 - L) - b\dot L$. While it lengthens a thrust is added, raised or lowered after each flight by the error of the apex reached against the apex asked.
- **Speed.** In flight the toe is placed ahead of the centre by the neutral point and a correction,
$$x_f = \frac{\dot x\,T_s}2 + k_v(\dot x - \dot x^* - b),$$
with $T_s$ the last stance's duration and $b$ a bias raised hop by hop by $0.4(\dot x^* - \dot x)$. The placement alone left a steady error of the speed.
- **Attitude.** In stance a force across the leg, within the friction asked ($0.5\,F$), gives the moment about the centre
$$M = -k_\theta\,\theta - k_L\,L_c ,$$
the pitch back to level and the angular momentum about the centre $L_c$ to zero. In flight nothing acts on $L_c$, and the Moon's flights are long (0.9 s for an apex of 15 cm), so a stance that left a little would turn the rabbit far.

The force $f$ asked at the toe becomes joint torques $\tau = -J^{\mathsf T}f$ in the hip, knee and ankle.

**4.4 The launch.** A hop starts from a crouch (the *ready* phase, eased in over 0.5 s): the body lowered over the flat hind foot until the virtual leg is $L_r = 0.15$ m, the forelegs lifted. The first push extends the leg by $\Delta L = L_0 - L_r = 0.05$ m. To leave the ground with the energy for an apex $h$ above the rest height, the leg does the work $mg(\Delta L + h)$ over $\Delta L$, so
$$F = mg\Big(1 + \frac h{\Delta L}\Big).$$
It is asked at the point of the sole beneath the centre while the heel is down, the centre of pressure that puts the force through the centre. A point toe would tip back.

**4.5 Refusal.** The launch's force must be held by the ankle over the foot's length, $\tau_{\text{ankle}} = F\ell_{\text{foot}}$. A command the body cannot do is refused, with its reason:
- on the Moon, $F = 1.69 \times 1.624 \times (1 + 0.15/0.05) = 11.0$ N, and the ankle needs 1.10 N m of its 2.5;
- on the Earth a hop of 5 cm needs 33.1 N at the toe and 3.31 N m at the ankle. The rabbit stands and says so.

This build is a Moon robot. An Earth robot would need an ankle twice as strong, or a lighter body.

**4.6 Stopping.** Asked to stop, the speed asked is zero and the thrust is withdrawn. When the rabbit lands slowly ($|v| < 0.25\sqrt{g/1.62}$ m/s) it stands.

## 5. Three Methods of One Motion

`Simulation` advances the robot by one of three methods on one Lagrangian and one set of forces:

| Method | What it is | Order | Cost (s/s) |
|---|---|---|---|
| variational | the midpoint discrete action's stationary path (physis §1–3), the forces by Lagrange–d'Alembert, each step by Newton's method on jets | 2, symplectic | 0.74 |
| semi-implicit Euler | Newton–Euler's equations, the velocity first | 1 | 0.006 |
| Runge–Kutta | Newton–Euler's equations by the classical method | 4 | 0.015 |

The step is 0.5 ms for all three. The costs are the machine's seconds per second of the robot's, on one core of the development machine. Where Newton's method does not meet the variational step's equation (a contact's first touch), that step is taken by semi-implicit Euler and counted, and the variational state restarts from it.

**5.1 The energy's account.** Let $E = T + V$, with $V$ gravity's, the tendon's and the limits'. Then
$$E(t) - E(0) = \int_0^t P\,dt, \qquad P = \sum_i Q_i\dot q_i ,$$
where $Q$ is the generalized force of everything else: the motors, the gears' friction, the limits' damping and the ground. The *residual* $E - E_0 - \int P\,dt$, with $P$ integrated by the trapezoidal rule over each step, is the method's own error of energy.

**5.2 What the methods keep.** On the first 2.6 s of hopping in place, after a second of standing (the crouch, the launch, flights and landings):

| | variational | Runge–Kutta | semi-implicit Euler |
|---|---|---|---|
| the first apex, against Runge–Kutta's | 3.5 mm | — | 2.0 mm |
| the residual of energy | 0.021 J | 0.115 J | 0.216 J |
| the angular momentum about the centre in flight, its largest relative change | | $3.6\times10^{-9}$ | |

The three agree on the motion to millimetres, and they differ in what they keep. The energy's residual does not rank as the orders do. Contact is stiff (the foot's 50 g on the ground, pressed in by a landing's 2 mm, rings near 100 Hz) and nonsmooth (a first touch, a slip speed), and an explicit method's order is not reached across it. The variational integrator, symplectic, does not drift between contacts. Runge–Kutta keeps the flight's angular momentum to its truncation, $3.6\times10^{-9}$.

This is the multi-fidelity of ADR 0054 for a body with joints: one model, three integrators, and their costs and their errors shown beside the motion. The page lets one choose.

## 6. What the Checks Hold

`check_legged` (26 checks):
- **The body:** its mass; the legs' inverse kinematics against their forward kinematics, either knee, to $10^{-15}$ m; the standing posture's toe, heel and paw on the ground, pressed in by their share of the weight; a motor's torque–speed line, its braking and its electrical power; the potential's force against differences.
- **Standing on the Moon:** still (the centre moves by $3\times10^{-5}$ m in a second), level ($2\times10^{-5}$ rad), its feet carrying its weight ($10^{-5}$), its energy accounted for (1 mJ).
- **Hopping in place:** at least 6 hops in 8 s (9), the pitch within 15° (11°), the apex brought to the 0.15 m asked (within 8 mm), the speed near zero.
- **Hopping forward:** the mean speed over the last 5 s of 10 at the 0.3 m/s asked (within 1 cm/s), the pitch within 15°.
- **Stopping:** from 0.2 m/s it lands and sits, still and level.
- **The three methods:** the first apex's agreement, the ranking of their residuals, Runge–Kutta's flight momentum, their costs stated.
- **The Earth:** the hop refused with its reason, and the rabbit standing, its feet carrying its weight.

## 7. The Ground It Learns

*([ADR 0059](../../../development/decisions/0059-robots-in-their-environments.md); `ground.hpp`.)*

**7.1 The truth and what the robot is told.** Along the robot's line the ground's contact law varies, $x \mapsto (\mu, k, n, c)$. In Cosmopoiesis it is the world's: the Moon's regolith presumed by the planetary-science domain and realised with the world's seed. The robot is told the presumption: for each parameter a median $p_{50}$, a spread $s$ of its logarithm and a correlation length $\ell$. Evaluating the world's ground at every foot and stage is costly, so `TabulatedGround` reads it once a centimetre as the robot goes and interpolates linearly. Its error is a part in $10^4$ of a law that varies over a metre.

**7.2 The feet as sensors.** Sampled with the controller, at 500 Hz:
- **Sliding.** A foot sliding faster than five slip speeds has $\tanh(v_t/v_s) > 0.9999$, so Coulomb's law gives $|f_t|/f_n = \mu$. The reading carries the foot's force sensor's noise, 2 %.
- **At its deepest.** A foot at its deepest has $\dot\delta = 0$, so $f_n = k\,\delta^n$ and $f_n\,(1\text{ mm}/\delta)^n$ is the load at 1 mm. Its sinkage is known only from the leg's kinematics against the ground's expected height, to 10 %, so the load to 15 %. The deepest instant is caught at the samples, 2 ms apart. The reading is taken to carry 20 % in all.

A foot that holds says only that $\mu > |f_t|/f_n$; that bound is not used.

**7.3 The belief.** $\ln(p/p_{50})$ is a Gaussian field along $x$ with variance $s^2$ and covariance $s^2 e^{-(x-y)^2/2\ell^2}$. A reading is an observation of it with noise variance $\ln(1 + r)^2$. tyche's conditioning adds the readings one at a time; a reading within 1 cm of the last of its parameter is skipped. The controller counts on the friction $\mu_{50}e^{-s_\mu}$ its belief gives at the toe: one deviation below its median. It asks the ground for at most three quarters of that.

**7.4 Calibration.** The whitened innovations $z = L^{-1}d$ are independent standard normals when the presumption and the noise are the truth. So their mean square is near 1 on a ground drawn from the presumption, and far above where the ground is not. The checks:
- **A ground drawn from the presumption.** The mean squares are 0.64 for the friction and 0.97 for the load. Where its feet slid, the belief knows $\mu$ to 1.2 % (the rms of $\log$ belief/truth).
- **A slippery patch.** On $\mu = 0.3$ where 0.83 was presumed, a jump of 5.5 presumed deviations, the rabbit keeps hopping across, its pitch within 7°. Its feet read 0.3 within 2.6 %, and its mean square of 152 says the presumption was wrong there.
- **At Reiner Gamma.** On the world's regolith, through the engine, the belief at the toe is 0.937 against a truth of 0.937 after ten hops.

## 8. Open Questions

- **A reduced model and a typed transition.** The spring-loaded inverted pendulum (SLIP) is the hopping body reduced to its centre and a massless leg. A transition from the tree to SLIP and back, with what it keeps and loses stated, is the multi-formalism the user described. It belongs in the domain once two models of one body exist.
- **Rigid contact.** Complementarity or nonsmooth variational collision integrators in place of compliance (physis §9).
- **Terrain.** The ground is a level line. A slope along the path is physis's `GroundLine` already; heights from the data library's tiles are not.
- **Three dimensions.** The spatial tree, $SE(3)\times T^k$. The fish, and a rabbit that turns.
- **The ankle while sitting: a flaw of this build.** Sitting, the ankle is folded far from the tendon's rest. Its motor holds the stretched tendon with 1.31 N m continuously: 5.3 A, and 12.4 W of copper loss, nearly all of the 12.5 W the rabbit draws while it sits. Hopping forward it draws 3.3 W on average, a cost of transport of 5.0 at 0.24 m/s. A small motor would overheat holding that current. A built rabbit would need one of three changes: a tendon that goes slack when sitting (one-sided, its rest at the sitting posture), a brake at the ankle, or a sitting posture near the tendon's rest. The model answers this question; the choice is the design's.
- **State estimation.** The controller is given the true state. An estimator on the inertial unit and the encoders would give it a belief.

## References

- Alexander, R. McN. (2003), *Principles of Animal Locomotion*, Princeton University Press.
- Featherstone, R. (2008), *Rigid Body Dynamics Algorithms*, Springer.
- Hunt, K. H. and Crossley, F. R. E. (1975), Coefficient of restitution interpreted as damping in vibroimpact, *J. Appl. Mech.* 42, 440–445.
- Marsden, J. E. and West, M. (2001), Discrete mechanics and variational integrators, *Acta Numerica* 10, 357–514.
- Raibert, M. H. (1986), *Legged Robots That Balance*, MIT Press.
- Siciliano, B., Sciavicco, L., Villani, L. and Oriolo, G. (2009), *Robotics: Modelling, Planning and Control*, Springer.
