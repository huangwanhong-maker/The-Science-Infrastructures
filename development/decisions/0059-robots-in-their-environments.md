# 0059. Robots in Their Environments: What a World Presumes of a Place, What a Body Feels There, and the Fields at a Point

- **Status:** Accepted, provisionally. Written: every field at a point; the Moon's ground presumed and realised, and a rabbit that learns it; the ocean's water presumed and realised, and a fish that swims in it; the Earth's ground presumed by CRUST1.0's classes and realised, and a dog that trots on it and learns it (2026-10-31).
- **Date:** 2026-10-30
- **Decided by:** the user. After the rabbit (ADR 0058) failed to hop on Greenland's ice, they observed that a simulation needs the physical properties of each place (the friction of ice, of ground), which no data give. They proposed that a simulated world "create virtual sensor instance, to 'guess' some physical quantities of a spacetime point … presume some of their physical property, as we have no real data"; that "model produce data, and data as a virtual sensor … When real sensor attach, we can utilize real data"; and they recalled that "a global coordination [is] associated with a integrated field (like robot + env co-produce field), and has dynamics model(s)". They asked for three robots: "rabbit for moon, fish for ocean, a dog for the Earth ground", the rabbit to stay a Moon robot, and "a panel [that] can show all dimensions of the physical fields' values of a selected point". Then: "please begin design and implementations." Claude designed the work and wrote this record.

  The foundation's records are physis's [decision 0015](https://github.com/huangwanhong-maker/physis/blob/main/development/decisions/0015-stationary-random-fields.md) and [decision 0017](https://github.com/huangwanhong-maker/physis/blob/main/development/decisions/0017-bodies-in-a-fluid-and-lagrangian-dynamics.md). The work's record is Cosmopoiesis's [decision 0012](https://github.com/huangwanhong-maker/Cosmopoiesis-Lab/blob/main/development/decisions/0012-the-ground-a-robot-learns-and-the-fields-at-a-point.md).

## 1. Context

A world holds fields: the magnetic, gravitational and atmospheric ones, the interior, the Moon's crust (ADR 0053). A robot in it meets properties no field held: the friction of the ground under a foot, the stiffness of regolith, the density of the water a fish swims in. Nothing has measured them at most places. The rabbit stood on one contact law everywhere, the same on regolith and on ice.

Cosmopoiesis already treats a model as an instrument. A virtual instrument reads a domain's model and delivers readings with provenance `model`, in the same envelope as a physical one (its decision 0004). And a world already keeps what is known of a field apart from its truth: the Moon's crust is T2015_449 plus a hidden realisation of what it leaves unresolved (ADR 0052), and the rover's mission learns it from readings (ADR 0054).

## 2. Decision

### 2.1 What a World Presumes of a Place

Where nothing has measured a property, the world presumes it: a law for the property at a place, from what the world knows of the place. This is a model with its uncertainty, kept as fields of the world with provenance `model`:
- a material or class at the place (regolith on the Moon; on the Earth, ice, sea, sediment or rock, from what the world holds there);
- the class's law for each parameter: a median and a spread, published values where they exist and labelled assumptions where they do not;
- how the parameter varies over the place: a correlation length.

The law belongs to the domain that names the material: regolith to planetary science, the Earth's ground and ice to geophysics, which holds CRUST1.0 (§2.8), seawater to oceanography.

### 2.2 The Truth of a Simulated World, and of the Real One

A simulated world needs a truth for the body to press on, and the presumption is not one. If the physics used the presumption, a robot could never slip or be surprised. The truth is a realisation of the presumed law, drawn with the world's seed (provenance `synthetic`), as the Moon's unresolved crust is. On a real robot, the truth is the real ground, and nothing else changes.

### 2.3 What a Body Feels There

A robot is told the presumption: a virtual instrument's reading of the world's model. Its own sensors read the truth:
- a sliding foot feels the friction, $|f_t|/f_n = \mu$ by Coulomb's law;
- a foot at its deepest feels the stiffness.

Its belief of the place is the presumption conditioned on its readings (tyche's Gaussian conditioning, ADR 0054). The whitened innovations test the presumption: near 1 on a ground drawn from it, far above where the ground is not (ice on regolith). The controller acts on the belief.

### 2.4 The Fields at a Point, and What a Body and Its Place Make Together

At each point of a world, its fields together are one integrated field (Cosmopoiesis's formalism §3). The page shows every field of the world at the point chosen, every component with its unit, frame and provenance, or why it has none there. A robot and its place make fields together: the forces of its feet on the ground and their readings, kept where they were made; later a footprint's compaction and a fish's wake. They are drawn where they were made.

### 2.5 Three Robots

| Robot | Environment | Its body | Its environment's law | Status |
|---|---|---|---|---|
| a rabbit (ADR 0058) | the Moon's regolith | a hopping tree in its sagittal plane | friction, stiffness and damping of contact (planetary science `surfaces`) | written |
| a fish | the Earth's ocean | a head and four driven links in its horizontal plane | the upper ocean's density and current (oceanography `water`); the water's forces on slender bodies by strip theory (physis) | written |
| a dog | the Earth's ground | a quadruped in its sagittal plane, its four legs separate | ice, sediment or rock under its feet, by CRUST1.0 (geophysics `ground`) | written |

### 2.6 What Went Where (the First Part)

| Where | What |
|---|---|
| physis `analysis` | a stationary random field of a squared-exponential covariance, read at any point by the spectral representation (physis decision 0015) |
| planetary science `surfaces` (new) | the Moon's regolith as a foot meets it: the friction presumed from the angle of internal friction (Lunar Sourcebook), stiffness and damping assumed, each varying over a metre; its realisation; the ground as fields |
| robotics `legged` | the ground along the robot's line (tabulated), the feet as sensors, the belief of friction and load conditioned on their readings, the controller's friction from the belief; checks of calibration on a ground drawn from the presumption and of the misfit on a slippery patch |
| Cosmopoiesis | `field.at` and the POINT page's fields at a point; the world's fields `moon.ground.*`; the robot on the world's regolith, its ground told, felt, believed and (revealed) true, on the ROBOT page and as marks in the world |

### 2.7 What Went Where (the Fish)

| Where | What |
|---|---|
| physis `mechanics` | the Euler–Lagrange equations of any Lagrangian solved for the acceleration by jets, with semi-implicit Euler and Runge–Kutta on them; a tree in a fluid by strip theory, the water's kinetic energy (added mass) in the Lagrangian and its drag a force (physis decision 0017) |
| oceanography `water` (new, the domain's first C++) | the upper ocean's density presumed from the range of surface seawater, its near-surface current assumed, each realised for a world's seed, as fields |
| robotics `actuation` (new), `swimming` (new) | the geared motor, moved out of `legged` now that two robots use it; the fish: its body, its servos closing their own loops on the reduced inertia, its travelling wave, three methods, its wake |
| Cosmopoiesis | the world's fields `earth.ocean.*`; the fish as a session of the engine, refused away from the sea (CRUST1.0's water), its water the world's at the site; the ROBOT page's choice of robot, the fish seen from above with its wake and the current, and the fish in the world below the sea's surface |

### 2.8 What Went Where (the Dog)

| Where | What |
|---|---|
| geophysics `ground` (new) | the land's ground: its class by CRUST1.0's layers (ice, sediment, rock; the sea not ground), each class's presumed contact law (assumed, labelled so), realised for a world's seed, as fields |
| robotics `quadruped` (new) | the dog: its body, its trot on diagonal pairs (the stance pair's force by the moment's balance within the friction it believes; the swing pair to Raibert's placement; an integral of the speed's error), three methods, its feet learning its ground with the rabbit's belief |
| Cosmopoiesis | the world's fields `earth.ground.*`; the dog as a session, refused at sea; the ROBOT page's dog (its sites on Greenland's ice and Paris's sediment, its side view of four legs) and its drawing in the world |

The ground is set apart from the geography domain's planned work: geophysics holds CRUST1.0, from which the class is read. A finer classification (land cover from imagery, snow from a forecast) would come from geography and the data library.

**What the fish shows.** It swims by drag. At its speeds it beats with a Strouhal number near 2.5 where a real fish's is 0.2 to 0.4, because strip theory leaves out the momentum a tail sheds into its wake, most of a real fish's thrust. That is a finding about the model, recorded as the fish's first open question rather than hidden by tuning.

## 3. Alternatives Considered

- **One property per place, the same everywhere** (as the rabbit had). Rejected: a robot could not meet ice, and nothing could learn.
- **The presumption as the truth.** Rejected (§2.2).
- **A virtual sensor that guesses at run time, without a law.** Rejected: a guess without its spread cannot be conditioned or tested; the presumption is a law, and the virtual instrument reads it.
- **A grid of the ground's properties** over each site. Rejected for the spectral representation: it is read anywhere, needs no grid, and is the same on every machine.
- **Learning the ground by a neural network** from the feet's signals. Not now: the Gaussian belief is calibrated by construction, and its test says when the presumption fails.

## 4. Consequences

- The checks gain the stationary field's law (physis, 9), the regolith's presumption and realisation (planetary science, 11), the ground the rabbit learns (robotics, 10 more), and `field.at` and the robot's ground through the engine (Cosmopoiesis); then bodies in a fluid (physis, 9), the ocean's water (oceanography, 9), the motor (robotics, 6) and the fish (robotics, 15).
- Then the Earth's ground (geophysics, 13) and the dog (robotics, 18).
- **Two lessons of the contacts.**
  - A foot's contact on ice or rock is the rubber pad's, not the ground's; a stiff, undamped one made the dog's feet ring.
  - Coulomb's regularised friction is a damper $\mu f_n/v_s$ that an explicit step integrates only while $h\mu f_n/(v_sm) < 2$; the dog's feet slip at $v_s \ge 0.03$ m/s at a step of 0.25 ms.
- Open: the Earth's ground classified more finely (land cover, snow, the climate's classes, temperature for ice); plastic sinkage (a footprint that stays); temperature in ice friction; Lighthill's reactive thrust for the fish, and its wake as a field of the water.
