# 0054. Bodies, Their Dynamics, and What Is Known: Discrete Mechanics, the Robotics Domain, and a Mission That Learns

- **Status:** Accepted, provisionally. One vehicle on a sphere at two fidelities, one instrument, a Gaussian knowledge. Terrain, contact, the estimation of a body's own pose and several bodies sharing a knowledge are open.
- **Date:** 2026-10-29
- **Decided by:** the user. After the unified fields (ADR 0053), they asked for "physical agent(s) building and rendering, simulation". Agents are not distinguished from other physical objects. Each has its own local coordinates, dynamics and fields, coupled in one world whose fields are "gaussed/assimilated" while "robot can provide us new info": "We like a god!", a "triple-relation". They asked whether dynamics should come from an action, "the field theory one", and whether a "path integral [is] even needed". They asked for "different dynamics system, and even [to] translate between them", since some uses need accuracy and some real time. Then: "you may add a domain robotics, if needed".

  A discussion with another assistant was read critically. Kept: its embodied agents, its triple relation, its experiment of a rover discovering a hidden synthetic field, and its fidelities. Not kept: its symbolic variation of actions, its general solver framework, and its reordering of the work's slices. Claude designed and wrote the work and this record.

  The foundations' records are physis's [decision 0013](https://github.com/huangwanhong-maker/physis/blob/main/development/decisions/0013-discrete-mechanics-and-the-law-of-a-field-known-in-part.md) and calx's [decision 0012](https://github.com/huangwanhong-maker/calx/blob/main/development/decisions/0012-symplectic-integrators-and-small-systems.md). The work's record is Cosmopoiesis's [decision 0008](https://github.com/huangwanhong-maker/Cosmopoiesis-Lab/blob/main/development/decisions/0008-a-body-that-senses-and-a-mission-that-learns.md).

## 1. Context

Cosmopoiesis's formalism held agents as bodies with ports, sensors, actuators and controllers (its Definition 8.1), and a field known in part (its §8.7). Nothing yet moved through a world and learnt from it. The foundations lacked four things:
- a way to state a body's dynamics once and integrate it without losing what it conserves;
- rigid bodies;
- constrained motion with its forces;
- the law of a field between two points, with the statistics that condition it on readings.

No domain held a vehicle, its controller, its plan or its sensors' mounts.

## 2. Decision

### 2.1 What Went Where

| Where | What |
|---|---|
| calx `integrators`, `linalg` | Störmer–Verlet and RATTLE, the constraint's forces returned (the normal force a wheel's friction depends on); small dense systems without allocation (calx decision 0012) |
| physis `mechanics` | motion from a discrete action: midpoint and trapezoidal discrete Lagrangians, the discrete Euler–Lagrange equations by Newton's method on calx's jets, forces by Lagrange–d'Alembert, the Legendre transform, discrete Noether, the Euclidean action; rigid bodies by the Lie group variational integrator, mass properties composed from parts, wrenches from a potential (physis decision 0013) |
| physis `analysis` | the covariance of a random potential field at any two points, of its potential and its field, by the addition theorem and the chain rule; tabulated at one radius (physis decision 0013) |
| tyche `spatial` (new) | Gaussian conditioning (simple kriging; the analysis step of optimal interpolation), extended one observation at a time by a row of the Cholesky factor; whitened innovations, the likelihood. A separable module of the inference foundation, whose ADR 0048 is in progress in another session |
| the robotics domain (new) | vehicles on a sphere at two fidelities and the translation between them; a waypoint follower; survey plans; sensors' mounts and noise drawn by index (§2.3) |
| Cosmopoiesis | the triple relation (its formalism §8.8); a mission: a rover surveying the Moon's crust, a hidden truth, the mission's Gaussian knowledge, the map, calibration over many truths; the page ROVER (its decision 0008) |

### 2.2 Dynamics from an Action; the Path Integral the Same Action

A body's dynamics is a Lagrangian, written once over its number type. Its motion is the stationary path of a discrete action, so symmetries' momenta are conserved exactly and coupling is a term of the action: two bodies and their interaction, a charge and its vector potential.

The user's question is answered in physis's theory note (§5). The Euclidean path integral is the same discrete action with imaginary time, `field_theory`'s free lattice field on the cycle of the time steps for a quadratic Lagrangian. It is not needed for classical motion, and it is there for a body's quantum statistics.

The check that the trapezoidal rule's integrator is Störmer–Verlet joins physis's discrete mechanics to calx's method.

### 2.3 The Robotics Domain

`domains/robotics` is founded at the user's word. It is a field at height 0 that stands on the foundations alone (calx, physis); the README's §4.4 had it in tier 4, on mechanical and electrical engineering, which are not yet founded. Its first modules:
- `vehicles`:
  - a pose on a sphere as one rotation ($SO(3)$ is the unit tangent bundle of $S^2$);
  - the kinematic unicycle, exact for a constant speed;
  - the dynamic vehicle, a mass held to the sphere by RATTLE about a local origin, its speed loop and the ground's Coulomb friction solved exactly;
  - the translation between them, with the dynamic one lagging by $v^*(\tau - h)$;
- `control`: following waypoints on a sphere;
- `planning`: the lawnmower carried by the exponential map, rings of candidates, the best candidate by value less cost (uncertainty sampling when the value is a posterior variance);
- `sensing`: mounts and white noise drawn by index.

A robot is a body with a controller. Its mechanics is physis's, and what names a vehicle, its wheels and its ground is the domain's. It names no body of the Solar System: the Moon's radius and GM come from those who use it.

### 2.4 The Triple Relation, and Calibration by Construction

A world's truth, a body in it and a knowledge of the world are kept apart:
- the body reads the truth only through its instruments;
- the knowledge changes only by delivered readings;
- a god's view sets the truth beside the knowledge without the knowledge reading it.

When the truth is a draw from the knowledge's prior, the posterior is calibrated by construction. Its errors over its standard deviations are standard normal, and so are its whitened innovations (Cosmopoiesis's formalism, Proposition 8.14). This is how an inference engine is tested without real data. The first mission does it at Reiner Gamma, with the truth T2015_449 and its unresolved part to degree 1440.

## 3. Alternatives Considered

- **A general framework of fidelities, or a multi-fidelity optimiser.** Deferred: the first case, the vehicle at two fidelities, is written in its domain. The pattern is abstracted when it appears twice (CLAUDE.md §2.2).
- **The agent loop in a foundation or in robotics now.** Deferred: sensing, believing, deciding and acting are written once, in Cosmopoiesis's mission. They move to robotics when a second work needs them.
- **Robotics on mechanical and electrical engineering**, as the README's map has it. Not possible yet; the domain declares them when they are founded.
- **MuJoCo, Chrono or another engine behind the interfaces** (CLAUDE.md §2.2, borrow implementations). Deferred to contact and articulated bodies, which are not yet needed.
- **A rover drawn on the three-dimensional Moon.** Not yet. It drives on the reference sphere, which lies above the maria's terrain, and would float there.

## 4. Consequences

- New checks:
  - calx: 14 (symplectic 10, small 4);
  - physis: 49 (discrete mechanics 25, rigid bodies 18, covariance 6);
  - tyche `spatial`: 13;
  - robotics: 48 (vehicles 21, control 9, planning 9, sensing 9);
  - Cosmopoiesis: 26 (the page's arithmetic, a survey at both fidelities and by both plans, the truth hidden, the posterior's z over 24 truths).
- Findings worth keeping:
  - The free rigid body's discrete energy is exact: Moser and Veselov's integrals, with $J_d^2 = \alpha I + \beta J^{-1}$ in three dimensions.
  - The textbook pitch libration is a small body's limit: a box at $3\times10^{-2}$ of its orbit's radius is shifted by $2\times10^{-3}$, through the orbit's epicycle near resonance.
  - The law of T2015's unresolved part puts its sources at 6.6 km and its power at degree 450. A rover's kilometre sees the field's value and gradients there, not finer structure: at Reiner Gamma the mission's error in $\lvert B\rvert$ is 0.3 nT, against T2015's 4 nT.
- The domains are sixteen: CLAUDE.md, the README and `domains/README.md` gain robotics at height 0.
- tyche is published at the user's word (`https://github.com/huangwanhong-maker/tyche`), with the modules of its founding session and `spatial`, and attached here as the submodule `foundations/tyche`; its ADR 0048 remains that session's to write.
- Open:
  - terrain following;
  - contact and the chassis as a rigid body;
  - estimation of a body's own pose (tyche `assimilation`);
  - a non-stationary prior;
  - several bodies sharing one knowledge;
  - the solar wind in the readings.
