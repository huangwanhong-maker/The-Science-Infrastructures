# A Fish: a Body with Joints in Water

*Status: draft (2026-10-31). Serves use, learning and research.* The theory of the module `swimming` ([ADR 0059](../../../development/decisions/0059-robots-in-their-environments.md); `fish.hpp`). Its mechanics is physis's planar trees (decision 0014) and bodies in a fluid (decision 0017, its theory note §8). What is here names the fish: its body as built, its servos, its controller, and what its swimming shows of the model.

## 1. The Body as Built

A rigid head of 1.20 kg (the battery and the electronics) and a chain of four links behind it, joined by four driven joints:

| Link | Length (m) | Height (m) | Mass (kg) |
|---|---|---|---|
| head | 0.25 | 0.12 | 1.20 |
| body | 0.14 | 0.11 | 0.35 |
| tail body | 0.12 | 0.09 | 0.20 |
| peduncle | 0.07 | 0.04 | 0.06 |
| caudal fin | 0.10 | 0.18 (its span) | 0.04 |

The fish is 0.68 m long and weighs 1.85 kg. Its configuration space is $SE(2)\times T^4$ in the horizontal plane: the head's place and heading, and the joints' angles. It is neutrally buoyant, so gravity does not act in that plane. The water's density and current are the world's where it swims; the drag coefficients are assumed: $C_n = 1.2$ across, $C_f = 0.03$ along, for skin friction and form together.

## 2. The Servos

Each joint is a geared motor (the domain's `actuation`) driven as a hobby servo is. The controller gives it an angle and a rate 200 times a second. Its driver closes its own loop at every instant:
$$\tau = \tau_p\tanh\!\Big(\frac{k_p(q^* - q) + k_d(\dot q^* - \dot q)}{\tau_p}\Big),$$
saturating smoothly at its peak. The gains are $k_p = I\omega^2$ and $k_d = 2\zeta I\omega$, with $\omega = 120$ rad/s, $\zeta = 0.9$, and $I$ the joint's inertia as it feels it:

**Remark (the reduced inertia).** On a free body a joint turns its links one way and the head the other. The inertia it feels is $1/(M^{-1})_{jj}$, the reduced inertia, where $M$ is the mass matrix of the body and the water's added mass. It is not the diagonal $M_{jj}$, which holds the head still: at the neck the diagonal is 0.45 kg m² and the reduced inertia 0.005, ninety times less. Gains set from the diagonal, with the torque held 5 ms, were unstable at once. Gains set from the reduced inertia with the torque held were stable but too weak against the water: the tail barely moved. A servo that closes its loop continuously is both stable and stiff.

## 3. The Controller

A travelling wave (a central pattern generator, Ijspeert 2008). Joint $i$ follows
$$q_i^* = A_i\sin(\varphi - i\lambda) + b_i, \qquad \dot\varphi = 2\pi f,$$
with amplitudes 0.05, 0.10, 0.20 and 0.32 rad from the neck to the tail, as in carangiform swimming, and a lag of $\lambda = 1.1$ rad between joints, a wave of about one body length.
- **Speed:** the frequency is integrated toward the speed asked, $\dot f = 0.8(v^* - v)$, between 0.4 and 3 Hz.
- **Turning:** a bias on the neck and the waist. Its sign is found, not derived: biased clockwise, the beating body turns left.
- **Stopping:** the wave stops and the body straightens; it glides to rest.

## 4. Three Methods of One Motion

The variational integrator moves the discrete action of the body and the water, the drag and the servos by Lagrange–d'Alembert. Semi-implicit Euler and Runge–Kutta move the Euler–Lagrange equations solved by jets (physis's `lagrangian.hpp`). Over 6 s of swimming, from the checks:

| | variational | Runge–Kutta | semi-implicit Euler |
|---|---|---|---|
| distance swum, against Runge–Kutta's | $1.6\times10^{-5}$ | — | $-1.9\,\%$ |
| cost (s/s) | 1.1 | 1.0 | 0.26 |

Runge–Kutta's account of the energy of the body and the water, against the work of the drag, the servos and the gears, closes to 4 mJ.

## 5. What Its Swimming Shows

The checks:
- **Speed.** Asked 0.15 m/s, it swims at 0.125 m/s on average (the mean over the last 8 s of 20).
- **Course.** Its heading stays within 5° of its course and its line within 2 mm.
- **Turning and stopping.** It turns 60° to the left in 8 s, and glides to 0.034 m/s within 10 s of stopping.
- **Current.** Holding still in a current of 0.2 m/s, it is carried at 0.194 m/s.
- **Wake.** The force its tail gives the water points backward on average, $-0.017$ N. That is the momentum it sheds into its wake, which its drag balances.

**It swims by drag, not as a fish does.** To swim at 0.18 m/s it beats at 2.5 Hz with its tail 0.18 m from side to side. That is a Strouhal number $fA/U$ near 2.5, where real fish swim at 0.2 to 0.4 (Triantafyllou, Triantafyllou and Grosenbaugh 1993). The reason is the model's: strip theory gives each strip its reactive and resistive forces but not the momentum the trailing edge sheds into the wake, which is most of a real fish's thrust (Lighthill 1971; physis's theory, §8). The fish is a faithful model of a body pushed by drag, and the measure of what strip theory leaves out.

## 6. Open Questions

- **Lighthill's reactive thrust.** Large-amplitude elongated-body theory adds the trailing edge's momentum flux. With it the beat should fall toward a fish's.
- **The water it learns.** Its belief of the current, from its own motion against its inertial unit, as the rabbit learns its ground.
- **Three dimensions.** Pitch, depth and buoyancy: the spatial tree.

## References

- Ijspeert, A. J. (2008), Central pattern generators for locomotion control in animals and robots: a review, *Neural Networks* 21, 642–653.
- Lighthill, M. J. (1971), Large-amplitude elongated-body theory of fish locomotion, *Proc. R. Soc. B* 179, 125–138.
- Triantafyllou, G. S., Triantafyllou, M. S. and Grosenbaugh, M. A. (1993), Optimal thrust development in oscillating foils with application to fish propulsion, *J. Fluids Struct.* 7, 205–224.
