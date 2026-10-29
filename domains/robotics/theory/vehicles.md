# Vehicles on a Sphere: Their Poses, Their Motion at Two Fidelities, and What Guides Them

*Status: draft (2026-10-29). Serves learning and use.* The theory beneath the robotics domain's first modules: `vehicles`, `control`, `planning` and `sensing` ([ADR 0054](../../../development/decisions/0054-bodies-dynamics-and-knowledge.md)). The mechanics it uses is physis's (constrained motion by RATTLE in calx, rigid bodies and discrete actions in physis `mechanics`). What this note adds is what names a vehicle: its pose on a surface, its wheels and the ground, its controller and its sensors.

## 1. The Pose of a Vehicle on a Sphere

**Definition 1.1.** A *pose* on the sphere $S_R = \{q \in \mathbb R^3 : |q| = R\}$ is a pair $(q, e)$ of a place $q \in S_R$ and a heading $e$, a unit vector tangent there ($q \cdot e = 0$). Its *axes* are forward $e$, left $n = \hat u \times e$ and up $\hat u = q/R$ (x forward, y left, z up: ROS's REP 103).

**Proposition 1.2.** The map $(q, e) \mapsto A = [\,e \;\; n \;\; \hat u\,]$ is a bijection from the poses on $S_R$ onto $SO(3)$.

*Proof.* The columns of $A$ are orthonormal and $e \times n = e \times (\hat u \times e) = \hat u$, so $A \in SO(3)$. Conversely a rotation's third column is a unit vector $\hat u$, which gives the place $R\hat u$, and its first column is a unit vector orthogonal to it, the heading. The two maps are inverse to each other. $\square$

*Remark.* This is the classical identification of $SO(3)$ with the unit tangent bundle of $S^2$. A pose is therefore one rotation, with no constraint to keep: a vehicle's configuration space on a sphere is a Lie group.

**Proposition 1.3 (motion).** A vehicle driving forward at speed $v$ and turning left at the rate $\omega$ about its up axis has $\dot A = A\, w^\wedge$ with the body angular velocity $w = (0, v/R, \omega)$. For $v$ and $\omega$ constant, $A(t) = A(0) \exp(t\, w^\wedge)$. Its place moves on a circle of the sphere, at the constant angle $\arctan(v/(R\omega))$ from the fixed axis $A(0) w/|w|$, and returns to its start, heading included, after $2\pi/|w|$.

*Proof.* With $\dot A = A w^\wedge$, $\dot{\hat u} = A (w \times e_z) = (v/R)\,A e_x = (v/R)\, e$, the motion forward at speed $v$, and $\dot e = A(w \times e_x) = \omega\, n - (v/R)\, \hat u$, a turn to the left at $\omega$ and the dip that keeps $e$ tangent. For constant $w$ the solution is the exponential. The spatial angular velocity $A(t) w = A(0) \exp(t w^\wedge) w = A(0) w$ is constant, so the motion is a rotation about a fixed axis, whose angle with $\hat u$ has cosine $w_z/|w| = \omega/|w|$, that is tangent $(v/R)/\omega$. The rotation is the identity again after the angle $2\pi$, at $t = 2\pi/|w|$. $\square$

*Remark.* The period is $2\pi/|w| = 2\pi/\sqrt{\omega^2 + (v/R)^2}$, not $2\pi/\omega$: the turn also carries the vehicle round the sphere. On the Moon, at $v = 0.2$ m/s and $\omega = 0.01$ rad/s, the two differ by $6.6 \times 10^{-11}$ of a turn.

**Lemma 1.4 (precision).** Composing rotations step after step, a pose gathers about $\varepsilon$ of angle a step, which is $R\varepsilon \approx 4 \times 10^{-10}$ m of place on the Moon. A place kept as a local origin $o \in S_R$ and a displacement $x$ (so $q = o + x$) gathers only the rounding of the displacements, about $\varepsilon\,|\Delta x|$ a step.

Both vehicles keep their place so (§2, §3). Their displacement in a step is computed without cancellation: $S e_z - e_z = 2 s_0 (s \times e_z) + 2 s \times (s \times e_z)$ for the step's rotation $S$ with quaternion $(s_0, s)$. The origin moves to the vehicle when $|x|$ exceeds a kilometre. Its rounding is then that of one sum, $\mathrm{ulp}(R) \approx 2.3 \times 10^{-10}$ m.

## 2. The Kinematic Vehicle

**Definition 2.1.** The *kinematic vehicle* has the state (pose, forward speed) and takes a command $(v^*, \omega^*)$. Over a step $h$:
- its speed moves toward $v^*$ by at most $a_{\max} h$;
- its yaw rate is $\omega^*$, clamped to $\omega_{\max}$;
- its pose turns by $S = \exp(h\,(0, \bar v/R, \omega)^\wedge)$, where $\bar v$ is the step's mean speed.

By Proposition 1.3 the step is exact when the speed is constant over it. While the speed slews it is of order 2: the curvature $\omega/v$ then varies within the step.

It is the model of a vehicle whose wheels roll without slipping and whose motors keep the speed they are given. It is fast, and it is the limit of §3's vehicle (Proposition 3.4).

## 3. The Dynamic Vehicle

**Definition 3.1.** The *dynamic vehicle* is a point mass $m$ at $q$, held to $S_R$ by the holonomic constraint $g(q) = 0$. These act on it:
- gravity, $-GMm\, q/|q|^3$;
- the ground's normal force, the constraint's multiplier;
- the wheels' drive force $F$ along the heading;
- the ground's Coulomb friction: rolling resistance $R_r = C_{rr} N$ against the forward motion, and lateral friction $\mu_{\text{lat}} N$ against sideslip.

**The scheme.** A Strang splitting.
1. A half kick by the wheels and the ground.
2. A RATTLE step (Andersen 1983; Hairer, Lubich and Wanner 2006, §VII.1.4) for gravity and the constraint. It is symplectic on the constraint manifold and of order 2.
3. A second half kick.

Each kick is solved exactly for the forces held over the step.

**Proposition 3.2 (Coulomb kicks).** The solution of $\dot v = (F - R_r \operatorname{sgn} v)/m$ with stiction has three cases:
- it moves at the constant rate $(F - R_r \operatorname{sgn} v)/m$ while $v \ne 0$;
- reaching $0$, it stays there while $|F| \le R_r$;
- from $0$ it moves off at $(F - R_r \operatorname{sgn} F)/m$ when $|F| > R_r$.

The sideways kick is the case $F = 0$.

*Proof.* On each interval where $\operatorname{sgn} v$ is constant the right-hand side is constant. At $v = 0$ static friction balances any drive within $R_r$, which is the definition of Coulomb's law with stiction. $\square$

**Corollary 3.3.** A drive that changes the speed at a constant rate is integrated by the trapezoidal rule: the drift moves the vehicle at the speed after the first half kick, which is the mean of the step's two speeds. Coasting from $v_0$, the vehicle stops in $v_0^2/(2 C_{rr} g)$. Sideslip from $v_0$ dies in $v_0^2/(2\mu_{\text{lat}} g)$.

**The constraint about a local origin.** With $q = o + x$, $|o| = R$:

$$g(x) = \frac{|x|^2 + 2\, o \cdot x}{2R}, \qquad \nabla g = \frac{o + x}{R}.$$

Here $g = |o + x| - R$ to first order, a height in metres, and no large numbers cancel. The multipliers of RATTLE's two halves are $\lambda$ and $\mu$, and the ground's normal force over the step is $N = -(\lambda + \mu)/2$, since $\nabla g$ is the outward normal. On the sphere at constant speed, $N = m(g - v^2/R)$.

**The wheels.** A motor's speed loop asks for

$$F = m\,\frac{v^* - v}{\tau} + C_{rr} N \operatorname{sgn} v^*,$$

the error over the loop's time constant plus the rolling resistance it must overcome. $|F|$ is limited by the ground's traction $\mu N$, the motors' force, and their power over the speed. Under this loop, away from the limits, the speed follows the discrete loop exactly: $v_k = v^*\,(1 - (1 - h/\tau)^k)$. When asked too much, the vehicle accelerates at $(\mu - C_{rr})\, g$.

**The heading.** The heading turns at the commanded yaw rate (a skid-steered vehicle's yaw is much faster than its translation). It is carried from place to place by the rotation that takes one up to the next, the discrete parallel transport along the step's great circle, and made tangent again.

**Proposition 3.4 (the fidelities and their translation).** The two vehicles share the state's form: origin, place from it, pose, speed. The translation between them is the following.
- *project*: forward speed $=(p\cdot e)/m$.
- *lift*: $p = m v e$, with $N = m(g - v^2/R)$.

Lifting then projecting is the identity. Driving straight from rest at $v^*$, with the kinematic vehicle reaching its speed within a step, the dynamic vehicle lags the kinematic one by

$$v^*(\tau - h)$$

once the loop has settled: its own lag $v^*\tau$, less a step.

*Proof.* The dynamic vehicle's speed deficit sums to $h\sum_k (v^* - v_k) = v^* h \sum_k (1 - h/\tau)^k = v^*\tau$. Its trapezoidal drift returns half a step of it, $v^* h/2$. The kinematic vehicle's first step has the mean speed $v^*/2$, a deficit of $v^* h/2$. The difference is $v^*\tau - v^* h/2 - v^* h/2$. $\square$

So the kinematic vehicle is the dynamic one's limit $\tau \to 0$ on firm ground. The engine chooses between them by what a run needs: the kinematic vehicle for planning and real time, the dynamic one where slip, slopes or loads matter.

## 4. Following Waypoints

A controller sees what its vehicle believes, not the world. Given a pose it returns a command, the same command for both vehicles. Toward a waypoint $w$ from $(q, e)$:
- the bearing is $b = \widehat{w - (w\cdot\hat u)\hat u}$;
- the heading error is $\psi = \operatorname{atan2}((e \times b)\cdot \hat u,\ e \cdot b)$;
- the yaw rate is $k_\psi \psi$;
- the speed is $\min(v_c \max(0, \cos\psi),\ k_d\, d)$, where $d$ is the distance left.

The heading error decays at the rate $k_\psi$ while the bearing turns more slowly than that, so a vehicle facing away turns in place before it goes. Near the end the distance left decays at the rate $k_d$, and the waypoint is taken within its arrival radius in a finite time. It is a proportional heading controller, a relative of pure pursuit (Coulter 1992), on the sphere's geometry. Computed from unit vectors, a bearing between places $d$ apart carries the rounding $R\varepsilon/d$.

## 5. Planning a Survey

**A lawnmower** covers a rectangle with parallel lines driven alternately. It is laid out in the tangent plane at its centre and carried to the sphere by the exponential map $\exp_c(d) = R(\cos(|d|/R)\,\hat c + \sin(|d|/R)\,\hat d)$, which keeps distances and directions from the centre exactly. A line's midpoint lies at $\tan\rho/\rho \approx 1 + \rho^2/3$ of its offset, $\rho$ the half-line's angle. Over a rectangle a hundredth of the radius across, the spacing changes by less than $10^{-4}$ of itself.

**An adaptive survey** chooses each next place from candidates by a value less $\lambda$ times a cost. When the value is the posterior variance of the field to be mapped, this is greedy maximum-entropy sampling, which for a Gaussian field is uncertainty sampling (Shewry and Wynn 1987). The value is the caller's: planning names no field.

## 6. Sensors on a Vehicle

A mount is a lever arm $\ell$ in the vehicle's axes and an orientation $O$, the rotation from the sensor's axes to the vehicle's. With the vehicle's pose $A$ and place $q$, the sensor is at $q + A\ell$ and its axes are $AO$. A vector $v$ of the world reads $(AO)^{\mathsf T} v$.

White noise is drawn by index from calx's counter-based stream at (seed, channel, $k$): sample $k$ is the same whenever it is drawn, and a run is the same on every machine (Salmon et al. 2011).

## 7. Open Questions

1. **Terrain.** The surface is a sphere. Heights from the data library's tiles would make the constraint $|q| - R(\hat q) = 0$, with slopes that the dynamic vehicle feels and the kinematic one does not.
2. **Wheels and soil.** Bekker–Wong terramechanics (slip, sinkage, the drawbar pull of regolith) in place of Coulomb friction.
3. **The chassis as a rigid body** (physis `mechanics`'s Lie group integrator) with suspension and contact, for rough ground.
4. **Estimation.** The controller is given the true pose. An estimator (odometry, a sun sensor, star tracking, the Kalman filters of tyche's planned `assimilation`) would give it its belief.

## References

- Andersen, H. C. (1983), Rattle: a "velocity" version of the Shake algorithm for molecular dynamics calculations, *J. Comput. Phys.* 52, 24–34.
- Coulter, R. C. (1992), *Implementation of the Pure Pursuit Path Tracking Algorithm*, CMU-RI-TR-92-01, Carnegie Mellon University.
- Hairer, E., Lubich, C. and Wanner, G. (2006), *Geometric Numerical Integration*, 2nd ed., Springer.
- Open Robotics, *REP 103: Standard Units of Measure and Coordinate Conventions*, https://www.ros.org/reps/rep-0103.html.
- Salmon, J. K., Moraes, M. A., Dror, R. O. and Shaw, D. E. (2011), Parallel random numbers: as easy as 1, 2, 3, *SC '11*.
- Shewry, M. C. and Wynn, H. P. (1987), Maximum entropy sampling, *J. Appl. Statist.* 14, 165–170.
- Wong, J. Y. (2008), *Theory of Ground Vehicles*, 4th ed., Wiley (rolling resistance, traction, terramechanics).
