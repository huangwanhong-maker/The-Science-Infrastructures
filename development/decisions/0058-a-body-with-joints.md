# 0058. A Body with Joints: a Rabbit on the Moon, and Three Methods of One Motion

- **Status:** Accepted, provisionally. One robot in its plane, compliant contact on a level ground, three methods of integration. The spatial tree, rigid contact, terrain, a reduced model with a typed transition, and the fish are open.
- **Date:** 2026-10-30
- **Decided by:** the user. They asked whether coupled dynamics of objects was implemented, then: "make 1 or 2 physical robot as examples … and can visualize on the UI", "a rabbit robot, and a fish robot", "I love rabbit for moon." On the representation they asked whether a 2-complex, triangulated, would be insufficient and raise problems of correspondence. A kinematic graph of links and joints was proposed: the state in $SE(2)\times T^k$ first, the geometry attached to the links. Then: "so theoretical physics also adopt here. But also consider practical, and simulation reality. As we may really implement such robot in the real world … (1) make robot example models and their physical property, and has a theoretical physics base, (2) visualize the robot, (3) implement multiple the dynamics coupling, as we said we may has many methods, like action S based, (4) the robot controller, to do some example behaviors." Claude designed and wrote the work and this record while the user slept.

  The foundation's record is physis's [decision 0014](https://github.com/huangwanhong-maker/physis/blob/main/development/decisions/0014-articulated-bodies-in-the-plane-and-compliant-contact.md). The work's record is Cosmopoiesis's [decision 0011](https://github.com/huangwanhong-maker/Cosmopoiesis-Lab/blob/main/development/decisions/0011-a-rabbit-on-the-moon.md).

## 1. Context

ADR 0054 founded the robotics domain with vehicles: a body that rolls on a sphere, a point with a heading. A robot that could be built needs more:
- a body of several rigid parts, joined, with the motors that drive the joints;
- feet that touch the ground and leave it;
- a controller that reads what a robot reads and does something worth watching;
- the dynamics computed more than one way, with the costs and errors of each shown: the "many methods" of the user's request, the action among them.

physis's `mechanics` held the discrete action and the free rigid body (its decision 0013). It held no articulated body and no contact.

## 2. Decision

### 2.1 What Went Where

| Where | What |
|---|---|
| physis `mechanics` | planar trees of rigid links on a fixed or floating base: their kinematics and Jacobians, their Lagrangian over any number type, Newton–Euler (the recursive algorithm, the mass matrix), forward dynamics, momenta and energy, semi-implicit Euler and Runge–Kutta on them; compliant contact with a ground line, Hunt and Crossley's normal force and a regularised Coulomb friction (physis decision 0014) |
| the robotics domain, `legged` (new) | the rabbit as built (links, motors with their torque–speed lines, armature, gears' friction, copper loss, a tendon, the joints' limits, seven points of contact, shapes for drawing); its controller at 500 Hz (standing; hopping in place and forward; stopping; refusing); its simulation by three methods with the energy's account (§2.2–2.4) |
| Cosmopoiesis | the robot as a session of its engine (`robot.models`, `.start`, `.command`, `.step`, `.state`, `.stop`, `.run`); the page ROBOT: the side view, the body's numbers, its motors, feet and energy; the rabbit in the world at its site, at true scale, lit by the Sun (its decision 0011) |
| the visualization layer | the globe's controls' clearance above the ground a property, so a camera may come within a metre of a robot |

### 2.2 The Representation: a Kinematic Tree with Attached Geometry

A robot is a tree of links and joints: a graph, its links the vertices and its joints the edges. Its state is a point of $SE(2)\times T^k$ in its plane, and later of $SE(3)\times T^k$. Each link carries its mass properties and, apart from them, its shapes for drawing and its points of contact.

A triangulated 2-complex was the alternative raised. It describes a surface, not a mechanism: its vertices would have to be tied to the links again, and a joint is not a cell of it. The tree is the configuration space's own structure. The geometry rides on it, so the drawing, the contacts and the dynamics never need to correspond by search.

### 2.3 The Rabbit as Built

A hopping robot of 1.69 kg ([theory](../../domains/robotics/theory/legged.md) §1–3):
- a torso of 1.30 kg with its battery at the ends;
- two-link hind legs and feet, two-link forelegs;
- five geared motors with a torque–speed line, the driver's peak current, armature, friction and the copper loss;
- an Achilles tendon at the ankle and the joints' limits, both elastic and in the Lagrangian;
- compliant feet on a regolith-like ground.

Its parameters are those of a plausible build, not a catalogue's; a built robot would put its datasheets' in their place.

On the Moon it stands, hops in place to the apex asked, hops forward at the speed asked, and stops to sit. On the Earth this build's ankle cannot launch a hop, and the controller refuses with the arithmetic: a hop of 5 cm needs 3.31 N m at the ankle, and its motor gives 2.50. It is a Moon robot.

### 2.4 The Controller

Raibert's three parts (height, speed, attitude), with what a body whose hip is far from its centre needs (theory §4):
- a virtual leg from the toe to the centre of mass, so the ground's force turns nothing;
- the pitch and the angular momentum about the centre servoed in stance, since the Moon's flights are long;
- an integral of the speed's error;
- a crouch and a launch whose force does the work the apex asks;
- gains set by each joint's inertia and a bandwidth, so the 2 ms of held torque stays stable.

It reads what a robot reads: encoders, an inertial unit (here the true state), foot switches.

### 2.5 Three Methods of One Motion

One Lagrangian and one set of forces, three integrators (theory §5):

| | variational | semi-implicit Euler | Runge–Kutta |
|---|---|---|---|
| order | 2, symplectic | 1 | 4 |
| energy's residual (2.6 s of hopping) | 0.021 J | 0.216 J | 0.115 J |
| cost, s per s | 0.74 | 0.006 | 0.015 |

They agree on the first apex within 3.5 mm. The residual of the energy's account, $E - E_0 - \int P\,dt$, is shown live as each method's own error. This is the answer to the user's "many methods": the action is one of them, the structured end. The explicit methods are the cheap end, which a controller is designed and run on.

## 3. Alternatives Considered

- **A 2-complex (a triangulated surface) as the representation.** Rejected (§2.2).
- **Three dimensions first.** Deferred: the plane shows every part of the problem with checks in closed form. The spatial tree comes with the fish.
- **A physics engine behind the interface** (MuJoCo, Bullet, Drake). Not adopted for this: the rabbit is how articulated dynamics, contact and control are learnt, and its three methods are compared on the infrastructure's own action (CLAUDE.md §7, our own C++ where it is how a method is learnt). An engine may sit behind the interface later, held to these checks.
- **Rigid contact.** Deferred (physis decision 0014).
- **A fish first.** The user chose the rabbit for the Moon. The fish needs the spatial tree and a fluid's forces.

## 4. Consequences

- physis `mechanics` gains `check_physis_multibody` (16). The robotics domain gains `legged` and `check_legged` (26). Cosmopoiesis's checks gain the robot's arithmetic (13) and its engine's answers (12).
- Cosmopoiesis's engine answers seven more operations.
- Open: the spatial tree; rigid contact; terrain under the feet; the rabbit reduced to a spring-loaded inverted pendulum with a typed transition (the multi-formalism the user described); state estimation; the fish.
