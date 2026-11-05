# 0076. Physical Computing: Computation by the Dynamics of Matter, its Pipeline and its Substrates

- **Status:** Accepted, provisionally; being built (2026-11-06: the formalism and the simulators of the fluids, chemical and molecular families, on physis's decision 0022)
- **Date:** 2026-11-06
- **Decided by:**
  - **The request.** The user proposed that "emerging computation tech" deserves a domain, after passing on two assessments of candidate substrates. Their premise was that "any dynamical system, or a system that can be changed under some constraints, can potentially be utilized as a computation system". They gave three criteria: complex systems (chaos, nonlinearity, changing attractors, critical transitions, so that AI may discover knowledge in a virtual world), "high generativity capacity (so that we can simulate various phenomena)", and kindness to the environment ("important for climate, Earth environment, ecosystem").
  - **The discussion.** Claude proposed the domain (the idea of 2026-11-05) with a precise definition of a physical computer, and assessed the candidates of both lists and its own. The user agreed ("Ok, utilize the full list, also add" the seven candidates Claude added). Seeing that the tokens of one round cannot build them all, the user asked for "folders and planning for each". For the software they asked for a simulator of each substrate first, "like TVM pipeline (like including intermediate language/compiler/scheduler/backend etc …), and simulator is one of backend. And we not implement intermediate language, make it optional and future implement."
  - **What Claude did.** Claude founded the domain as planned: its pipeline, its fourteen families of substrates, and the plan of each family's simulators, with the checks each simulator will be held to. No code is written.

## 1. Context

The infrastructure computes on digital machines through calx, whose retargetable stack separates what is computed from how and where (calx's decision 0003, after TVM, Halide and MLIR). Matter also computes in its own way. A tank of water, a reaction, a memristor, a network of oscillators and a culture of cells each have dynamics that, encoded and read in the right way, carry out a computation. Two of the three criteria cannot be met on digital machines alone. A world for AI to discover in should show genuine critical transitions and chaos, and the energy and materials of computing should fall. Nothing in the infrastructure holds such substrates as computers. Neuroscience holds the brain, robotics bodies, physis the laws, and none of them holds the encoding, the readout, the tolerance and the cost that make a substrate a computer.

## 2. Decision

### 2.1 The Domain

A domain, `physical-computing`, holds computation by the dynamics of matter: the substrates as computers, how tasks are carried onto them, and what they cost. It is a field at height 0, on calx, morphe, physis and tyche, with heuresis when a substrate is configured by search. It is C++ from its first code, listed under engineering in the README's §4.4. Its status is `planned`. It rises to height 1 when its organoids declare neuroscience and its swarms declare robotics (§2.5).

### 2.2 The Definition of a Physical Computer

A physical computer is a substrate with its evolution under a control, an abstract model, an encoding $E$, a readout $R$ and a tolerance $\varepsilon$. It computes the model on a set $K$ when the square commutes within $\varepsilon$: $d(R(\Phi_u(E(m))), C(m)) \le \varepsilon$ for every $m \in K$. It is non-trivial when $E$ and $R$ are of a class fixed in advance and cheaper than the computation (Horsman, Stepney, Wagner and Kendon 2014). The theory note ([`domains/physical-computing/theory/physical-computation.md`](../../domains/physical-computing/theory/physical-computation.md)) states the rest of the theory:

- the four modes: analogue by conjugacy, reservoir by training, programmed by design, and learning by the substrate's own rules;
- why simulating and producing need opposite properties;
- what transfers from a substrate to its target (its normal forms and universality classes);
- the cost of a verified result, in energy and in materials;
- generativity as breadth (the class of laws a substrate can be made to obey, greatest for the analogue computer, reaction networks and DNA, where it is a theorem) and as depth (whether the substrate *is* an instance of its law).

### 2.3 The Pipeline

The domain's software follows calx's stack and TVM: what is computed is kept apart from how it is carried onto a substrate, when, and on what. It has six stages, each a directory of the domain:

| Stage | Directory | Role | First work |
|---|---|---|---|
| Task | `tasks/` | What is computed, substrate-neutral, with its tolerance and its check. There are six kinds: filters, flows, energy problems, functions, phenomena and learning | planned |
| Intermediate representation | `ir/` | A versioned text form of tasks, programs and schedules | **optional, future**: not implemented; until it exists, tasks are C++ values |
| Compiler | `compiler/` | Passes that lower a task onto a target: matching capabilities, scaling, encoding, readout, fitting and search. The target-specific lowering lives with each family of substrates | planned |
| Schedule | `schedule/` | When: injection, washout, virtual nodes, sampling, annealing, training and testing, repetitions; kept apart from the program, as Halide keeps it | planned |
| Runtime | `runtime/` | Runs a compiled program under its schedule on a backend, and returns readings with their provenance | planned |
| Backends | `substrates/<family>/` | One family of substrates per directory, each with its backends: a **simulator first**, later a device and a learned surrogate | planned; simulators first |

The contract and the accounts, in `computer/`, stand beside the pipeline. They hold the square and its tolerance, the invariants by which chaos is compared (statistics, not trajectories), the cost of a verified result, and the materials' five marks. Every run reports them.

The stack keeps six invariants:

1. **One directory per family.** A new substrate is a new directory under `substrates/`, or a new simulator in its family, and nothing else changes.
2. **The simulator is the reference.** A device is held to its family's simulator, as calx's CPU interpreter defines correct results.
3. **A declared tolerance.** Every task declares its tolerance, which is never silently loosened, as calx never silently reduces precision.
4. **The cost reported.** Every result reports its cost per verified result and its materials.
5. **Reproducible draws.** Noise is drawn from calx's streams, keyed by what it is for, so a run is the same on any number of threads.
6. **Laws from the foundations.** A simulator takes its laws from the foundations, and its directory holds only what names the substrate.

### 2.4 The Substrates

Every candidate is placed in one of fourteen families. The candidates are those of the two assessments the user passed on and those Claude added, among them the seven the user asked to add. Each family's README plans its simulators: the model, its source, its noise and non-idealities, what it stands on, the checks against closed forms, and the tasks it accepts.

| Family | Substrates |
|---|---|
| `analogue-electronics` | analogue computers (Shannon's general-purpose analogue computer), field-programmable analogue arrays, CMOS oscillators and oscillator-based Ising machines |
| `spintronics` | spin-torque oscillators, alone and coupled; magnetic networks and artificial spin ices; stochastic magnetic tunnel junctions |
| `memristive` | memristors and crossbars; phase-transition switches and oscillators (VO₂, NbO₂); phase-change memories (GeSbTe) |
| `probabilistic` | p-bits and probabilistic computers; thermodynamic sampling and linear algebra |
| `photonics` | interferometer meshes; delay-loop reservoirs; coherent Ising machines |
| `chemical` | the Belousov–Zhabotinsky medium and its programmable arrays; chemical reaction networks as programs; electrochemical oscillators |
| `molecular` | DNA strand displacement; enzymatic DNA networks (the PEN toolbox) |
| `fluids` | water in a tank; convection cells; rotating, stratified tanks; Faraday waves and walking droplets; capillary and foam networks; bubble logic; Phillips's hydraulic machine |
| `soft-matter` | ionic membranes and hydrogels, self-oscillating gels; liquid crystals and active nematics; active colloids |
| `mechanical` | mechanical metamaterials and hysterons; mechanical oscillators (metronomes); soft bodies; granular media and experimental landscapes |
| `learning-networks` | physical learning networks of resistors, flows and springs |
| `living` | mycelium; *Physarum*; ecological microcosms (chemostats, plankton, microbial food webs); synthetic gene circuits; biofilms; neural cultures and organoids |
| `quantum` | quantum reservoirs |
| `swarms` | robot swarms |

Beside them stands no substrate of ours: the Earth itself, observed through the data library, is the greenest substrate and the one of last appeal (theory note §3.3).

### 2.5 Substrates Whose Material Another Domain Studies

A substrate is this domain's object as a device: a memristor is a computer of titanium oxide, a chemostat a computer of algae and rotifers. The science of its material belongs to that material's domain. Where that domain is built, the simulator uses it and the domain declares it: the organoids will use neuroscience's `snn`, and the swarms robotics's bodies, in a module `swarm` that robotics will gain. Where that domain is only planned (ecology, cell-biology, microbiology, condensed-matter), the least model the simulator needs is written in the family, marked provisional. It moves to that domain when the domain is founded, leaving a pointer, as the Chemistry Workbench keeps its chemistry until a chemistry domain exists.

This amends the proposal (theory note §2.2), which would have composed living substrates only in works. That rule would have left a living substrate with no simulator until its domain was founded, and it would have put ecology at height 5 beneath any work that used a microcosm.

### 2.6 What goes to the Foundations

The general methods are written in the foundations first, each by its own decision, when a simulator first needs them:

| Where | What |
|---|---|
| physis `dynamics` (planned there) | discrete-time and delay dynamics; Lyapunov spectra and Kaplan–Yorke dimension; invariant measures; bifurcations, their normal forms and continuation; Floquet multipliers; test systems with known invariants (Lorenz 63 and 96, Rössler, Hénon, the logistic map, Mackey–Glass, Kuramoto–Sivashinsky) |
| physis `control` (planned there) | driven systems; the echo state property; conditional Lyapunov exponents |
| physis `statistical_mechanics` | Landauer's bound and the thermodynamic uncertainty relation, on the master equations; Ising and Potts models on graphs; sandpiles (Dhar's abelian sandpile, the Oslo model); Vicsek's self-propelled particles |
| physis `chemistry` | polynomial ODEs compiled to mass-action networks by dual rails |
| physis `discretisation` | Kirchhoff networks on graphs, whose conductances adapt by local rules (*Physarum*, coupled learning) |
| physis `continuum` | linear waves in a basin; the linear stability of Rayleigh–Bénard convection; two-layer quasi-geostrophic flow |
| physis `quantum` | registers of qubits, Pauli Hamiltonians, density matrices and their evolution |
| physis `symmetry` | the factorisation of a unitary into a mesh of two-mode rotations (Reck; Clements) |
| tyche `learning` | reservoirs and readouts; memory and information-processing capacity; next-generation reservoir computing, the digital baseline; NARMA |
| tyche `time_series` and `descriptive` | delay embeddings, false nearest neighbours, the correlation dimension, early warning signals; the energy distance |

### 2.7 The Order of Work

1. **Simulators for every family**, in the order the families' READMEs give: first the most rooted (analogue electronics, chemical, probabilistic, learning networks), each simulator with its checks, and the foundations' methods beneath them as each needs them. Then `computer/`, `tasks/` and a minimal compiler, schedule and runtime, enough to run one task end to end on each simulator.
2. **Compilers** for each family: the target-specific lowering of each kind of task it accepts.
3. **The intermediate representation**, when two or more families' compilers show what it must carry. Whether it is calx's IR extended (continuous time, a tolerance in place of a precision) is calx's decision.
4. **Devices**, each an instrument behind `observation/`, held to its family's simulator.

The case study for the first criterion follows the first step. A substrate is tuned towards a fold or a Hopf point, and an AI is asked to find its law and to warn of its transition from its readings alone.

## 3. Alternatives Considered

- **Substrates as calx backends.** Rejected. A calx backend executes calx's IR at a declared precision, and a substrate computes by its own physics, with an error that is part of that physics. The domain's stages follow calx's so that a common IR remains possible (§2.7, step 3).
- **A foundation of physical computing.** Rejected. Its general methods do go to the foundations (§2.6), but the substrates name objects in nature.
- **A first slice of three substrates** (the proposal's §4). Superseded at the user's request: a plan for every family now, and then a simulator for each.
- **Living substrates only in works** (the proposal's §2.2). Amended (§2.5).
- **One directory for each of the forty-odd substrates.** Rejected for families of shared physics. The substrates of a family share their simulators' laws and their lowering, and each family's README plans each substrate in it.

## 4. Consequences

- The domain `domains/physical-computing/` has its README and `module.json`. It has six directories of the pipeline and a theory note, each with its plan, and fourteen directories under `substrates/`, each family's plan.
- No code is written yet. The first code of each family is its simulator, and the foundations' methods it needs go to their libraries, each by a decision of its own.
- The root README (§4.4, §4.6), `domains/README.md`, CLAUDE.md (§1, §3.2, §5, §6.3, §9) and the decision log record the domain.
- Open: the C++ namespace (provisionally `physical_computing::<module>`, CLAUDE.md §8.1); whether living neural tissue and engineered organisms are in scope beyond simulation, and under what rule of ethics and containment; which work draws the domain (a workbench of substrates, or Cosmopoiesis).
