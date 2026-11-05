# physical-computing: Computation by the Dynamics of Matter

*Status: planned ([ADR 0076](../../development/decisions/0076-physical-computing.md)); no code yet.* Purpose: learning, use and research.

## 1. Purpose

Matter computes. A tank of water, a chemical medium, a memristor, a network of oscillators and a culture of cells each carry out a computation when their inputs are encoded and their states read in the right way. This domain holds such substrates as computers. For each it holds how a task is carried onto it, how the result is read and checked, and what a verified result costs in energy and in materials. It serves three aims, the user's:

1. **Worlds for discovery.** Substrates that simulate, or really produce, chaos, nonlinearity, changing attractors and critical transitions, so that AI may learn and discover in them.
2. **Generativity.** Substrates that can be made to obey many laws (breadth) or that are instances of a law (depth).
3. **Kindness to the environment.** Computation whose energy and materials are counted, per verified result.

The theory is in [`theory/physical-computation.md`](theory/physical-computation.md). It defines a physical computer as a square that commutes within a tolerance, with an encoding and a readout of a class fixed in advance. It then sets out the four modes (analogue, reservoir, programmed, learning), why simulating and producing ask for opposite properties, what transfers from a substrate to its target, the cost of a verified result, and generativity as breadth and depth.

## 2. The Pipeline

The software follows calx's retargetable stack and TVM, keeping apart *what* is computed, *how* it is carried onto a substrate, *when*, and *on what*:

```
task ──▶ [ir: optional, future] ──▶ compiler ──▶ schedule ──▶ runtime ──▶ backend
                                       │                                  │
                          target-specific lowering,          simulator (first) · device (later,
                          in each substrates/<family>/        through observation/) · surrogate
                                       └──────────── computer/: the square, its tolerance, the cost ─┘
```

| Stage | Directory | Role | Status |
|---|---|---|---|
| Task | [`tasks/`](tasks/) | What is computed, substrate-neutral, with its tolerance and check: filters, flows, energy problems, functions, phenomena, learning | planned |
| Intermediate representation | [`ir/`](ir/) | A versioned text form of tasks, programs and schedules | optional, future |
| Compiler | [`compiler/`](compiler/) | Passes lowering a task onto a target: capabilities, scaling, encoding, readout, fitting and search | planned |
| Schedule | [`schedule/`](schedule/) | When: injection, washout, virtual nodes, sampling, annealing, training and testing, repetitions | planned |
| Runtime | [`runtime/`](runtime/) | Runs a compiled program under its schedule on a backend; readings with their provenance | planned |
| Backends | [`substrates/`](substrates/) | One family of substrates per directory, each with its simulator first | planned |
| The contract | [`computer/`](computer/) | The square and its tolerance, invariants, the cost of a verified result, the materials' marks | planned |

The pipeline keeps six invariants:

1. A new substrate is a new directory or a new simulator in its family, and nothing else changes.
2. A family's simulator is the reference its devices are held to.
3. Every task declares its tolerance, which is never silently loosened.
4. Every result reports its cost and its materials.
5. Noise comes from calx's keyed streams, so that runs are reproducible.
6. Laws come from the foundations; a family holds only what names its substrates.

## 3. The Substrates

Fourteen families hold every candidate of the assessments the user passed on, and those Claude added ([`substrates/`](substrates/)). Each family's README plans its simulators: the model and its source, its noise and non-idealities, what it stands on, its checks against closed forms, and the tasks it accepts.

| Family | Substrates | Breadth, depth | Status |
|---|---|---|---|
| [`analogue-electronics`](substrates/analogue-electronics/) | analogue computers, field-programmable analogue arrays, CMOS oscillators and oscillator-based Ising machines | broadest: polynomial ODEs | planned |
| [`spintronics`](substrates/spintronics/) | spin-torque oscillators; magnetic networks and artificial spin ices; stochastic tunnel junctions | networks of one unit | planned |
| [`memristive`](substrates/memristive/) | memristors and crossbars; phase-transition switches and oscillators; phase-change memories | networks of one unit | planned |
| [`probabilistic`](substrates/probabilistic/) | p-bits and probabilistic computers; thermodynamic sampling and linear algebra | any Boltzmann law on spins; Gaussians | planned |
| [`photonics`](substrates/photonics/) | interferometer meshes; delay-loop reservoirs; coherent Ising machines | any unitary; delay dynamics | planned |
| [`chemical`](substrates/chemical/) | the BZ medium and its programmable arrays; reaction networks as programs; electrochemical oscillators | broadest: polynomial ODEs | planned |
| [`molecular`](substrates/molecular/) | DNA strand displacement; enzymatic DNA networks | broadest: any reaction network | planned |
| [`fluids`](substrates/fluids/) | water in a tank; convection cells; rotating, stratified tanks; Faraday waves and walking droplets; capillary and foam networks; bubble logic; Phillips's hydraulic machine | deepest for the atmosphere and ocean | planned |
| [`soft-matter`](substrates/soft-matter/) | ionic membranes and hydrogels; liquid crystals and active nematics; active colloids | deep | planned |
| [`mechanical`](substrates/mechanical/) | metamaterials and hysterons; metronomes; soft bodies; granular media and landscapes | deep; programmable graphs | planned |
| [`learning-networks`](substrates/learning-networks/) | networks of resistors, flows and springs that learn by local rules | the filters they learn | planned |
| [`living`](substrates/living/) | mycelium; *Physarum*; ecological microcosms; gene circuits; biofilms; neural cultures and organoids | deep for populations and cells | planned |
| [`quantum`](substrates/quantum/) | quantum reservoirs | open quantum dynamics | planned |
| [`swarms`](substrates/swarms/) | robot swarms | network dynamics of a programmable unit | planned |

The greenest substrate is none of these: the Earth itself, read through the data library, which a built substrate must outdo in control, repetition or speed to earn its place (theory note §3.3).

## 4. The Order of Work

1. **A simulator for every family**, each with its checks, in this order:
   - analogue electronics, chemical, probabilistic and learning networks, the most rooted;
   - then memristive, spintronics, photonics, mechanical and molecular;
   - then fluids, soft matter, living, quantum and swarms, whose laws the foundations must first gain (two-layer flow, Potts models, qubit registers; robotics's `swarm`).

   Beside them come `computer/`, `tasks/` and a minimal compiler, schedule and runtime, enough to run one task end to end on each simulator.
2. **Compilers**: each family's lowering of each kind of task it accepts.
3. **The intermediate representation**, once two families' compilers show what it must carry.
4. **Devices** through `observation/`, each held to its family's simulator.

## 5. Dependencies

The foundations calx, morphe, physis and tyche, and heuresis when a substrate is configured by search. The general methods the simulators need go to the foundations first (ADR 0076 §2.6):

- in physis: chaos and bifurcations, driven systems, the thermodynamic uncertainty relation, spin models and sandpiles, dual-rail compilation, adaptive Kirchhoff networks, waves and two-layer flow, qubit registers, unitary meshes;
- in tyche: reservoirs and capacities, embeddings, dimensions and early warnings.

No domain is declared yet. The organoids will declare neuroscience (its `snn`), and the swarms robotics (a module `swarm` it will gain), raising the domain to height 1. A substrate whose material's domain is only planned holds the least model its simulator needs, provisionally, until that domain is founded (ADR 0076 §2.5).

## 6. Open Questions

1. The C++ namespace: provisionally `physical_computing::<module>`.
2. Whether the IR is calx's, extended with continuous time and a tolerance in place of a precision (calx's decision).
3. Whether living neural tissue and engineered organisms are in scope beyond their simulation, and under what rule of ethics and containment.
4. Which work draws the domain: a workbench of substrates, or Cosmopoiesis, which constructs worlds.

## 7. References

The references are collected in the theory note, and each family's README adds its own. All were written from memory, and are to be checked against the papers before a module cites them.
