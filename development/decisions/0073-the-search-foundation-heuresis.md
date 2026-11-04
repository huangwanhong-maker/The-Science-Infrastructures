# 0073. The Search Foundation heuresis: Evolutionary Computation and Search over Spaces Too Large to Enumerate

- **Status:** Accepted, provisionally
- **Date:** 2026-11-05
- **Decided by:**
  - **The request.** The user proposed "a new module, the evolutionary computing model", for the stage of discovery that searches for solutions: "even in chemistry domain, it can be so challenging, we even want integrated with other domains". The user's framing was that the infrastructure serves AI as well as people, because "AI more comfortable on simulated, virtual world". They first proposed a domain named optimisation, "and evolutionary methods as one".
  - **The discussion.** Claude proposed a foundation instead: search names no object in nature, and a foundation cannot use a domain. The user agreed ("it should be a foundation") and chose the name *heuresis*. They asked that the first slice be "comprehensive as possible" and that the design be "generalizable, extendable" for "many future evolutionary algorithms, emerging computation devices, like even distributed evolutionary computing, or quantum evolutionary algorithms". They will create the repository on GitHub.
  - **What Claude did.** Claude designed the contracts and seams (heuresis's decision 0003) and chose the first eighteen modules (its decision 0004). It wrote the core and the benchmarks, and the other modules in parallel sessions on those contracts, all checked.

## 1. Context

Discovery has stages: a question, the search for a solution, an experiment, an evaluation.
- The infrastructure holds the experiments: the domains' models and simulations.
- It holds the evaluation: tyche.
- It did not hold the search. calx's `optimisation` has local, deterministic methods on ℝⁿ, and morphe has exact combinatorial optimisation.

The spaces the sciences need searched are too large to enumerate and have no gradient: designs, plans, controllers, programs, sequences, molecules. Such spaces are searched by populations, search distributions, archives of diverse solutions and surrogate models.

These methods name no object in nature, and they serve every domain. By the infrastructure's principle (§2.1; ADR 0042), they are written once, in a foundation.

The No Free Lunch theorems say that no search does better than random search on average over all problems. A search's power is the structure of its problem that it exploits. That is the argument for an integrative infrastructure: the domains supply the structure, as representations, constraints and models.

## 2. Decision

### 2.1 A Seventh Foundation

heuresis (Greek *εὕρεσις*, "finding, discovery") is a foundation at `foundations/heuresis`, a repository of its own (`https://github.com/huangwanhong-maker/heuresis`), attached as a submodule. It depends on calx, morphe, physis and tyche, and nothing in the foundations depends on it.

Its boundaries (its decision 0002):
- calx keeps local minimisers and linear programming;
- morphe keeps exact combinatorial optimisation, grammars and graph rewriting;
- tyche keeps statistical models and tests;
- physis keeps physical law and quantum states;
- the domains state their problems;
- the edges and works hold evaluators on devices, transports between machines, and proposers that call language models.

### 2.2 The Search Process and Its Seams

A search algorithm is a pure state machine with `ask` and `tell`, and it never calls its objective (heuresis's decision 0003).
- **Generations.** Generational algorithms update once a whole generation is told, in ticket order. Their states do not depend on the order of evaluations, the executor or the machines.
- **Keyed draws.** Every draw is keyed by what it is for, on calx's counter-based streams.
- **Provenance.** Every individual records its parents and the operation that made it, so a run is a genealogy.

Five seams keep the design open:
1. **Executors**: where evaluations run (serial, a thread pool, simulated asynchronous workers in the core; clusters and devices at the edges).
2. **Representations and kernels.**
3. **Search distributions and their samplers.** A quantum device, an annealer or a probabilistic-bit machine is a sampler of a law, and the core's exact classical emulation of that law is its reference, as calx's reference interpreter is for its backends.
4. **Migration between islands**, coordinated in the core and transported at the edges.
5. **External proposers**: a person, a language model or an agent through the API (ADR 0030), whose proposals enter as injected individuals with their origin. This is the shape of FunSearch and AlphaEvolve, with the infrastructure's simulations as their evaluators.

### 2.3 The First Modules

`search`, `variation`, `selection`, `benchmarks`, `genetic`, `local`, `strategies`, `swarms`, `distributions`, `quantum`, `annealing`, `multiobjective`, `diversity`, `programs`, `populations`, `models`, `landscapes`, `analysis` (heuresis's decision 0004, and its README for what each holds).

Each is checked against closed forms where the theory gives them (runtimes, invariances, stationary laws, takeover times, hypervolumes), against brute force on small spaces, or against reference implementations by their deterministic parts: 492 checks, all passing. CMA-ES's update is held to pycma 4.5.0 by telling both the same samples (to 5 × 10⁻¹⁶).

Planned: `graphs` (morphe's rewriting), `coevolution`, `neuroevolution` (tyche's `learning`), `noise`, `dynamics`, `configuration`.

### 2.4 What the Domains Gain

A domain that searches states its problem (a space, objectives, constraints, descriptors) and composes heuresis.
- robotics: controllers and gaits under its simulators;
- genomics: sequence design;
- urban science: plans under budgets, with several objectives;
- architecture: layouts;
- economics: designs of mechanisms;
- a future molecular science: molecules as graphs rewritten under valence rules.

Programs typed by physis's dimensions let symbolic regression discover laws in their physical form from a domain's data.

### 2.5 Requested from the Siblings

Each is to be written when its working tree is free of other sessions' changes:
- calx: the fast Walsh–Hadamard transform and low-discrepancy sequences;
- tyche: rank tests and effect sizes, k-means, Latin hypercube designs;
- physis: registers of qubits, gates and measurement;
- morphe: the algebraic simplification of terms.

## 3. Alternatives Considered

- **A domain named optimisation.** The user's first proposal. Rejected by the user on reflection: a domain holds objects in nature, a foundation cannot use a domain, and every domain and some foundations need search.
- **A module of calx.** Rejected: calx works on floating-point arrays, while evolutionary computation searches bit strings, permutations, programs and graphs, and needs morphe, physis and tyche.
- **Borrowing a framework** (pagmo, DEAP, pymoo, Nevergrad, EvoX). Rejected for the core:
  - their algorithms mostly call their objectives, so evaluations cannot be spread, delayed or made by something that is not a function;
  - their interfaces are not ours;
  - several are copyleft.

  They remain references; pycma and pymoo may give fixtures.

## 4. Consequences

- The foundations are seven: calx, morphe, physis, tyche, nomos, eikon and heuresis. The commit order is calx and morphe, then physis, tyche, nomos, eikon and heuresis, then the neuroscience domain, Cosmopoiesis and this repository.
- `docs/mathematics.md` §9 says where optimisation lives:
  - calx's when continuous and local;
  - morphe's when combinatorial and exact;
  - heuresis's when global, stochastic or black-box, over any representation;
  - physis's when it is the control of a dynamical system;
  - tyche's when it is the fitting of a model.
- The dependency checks know heuresis (`check-dependencies.mjs`, `check-alone.mjs`). The root build adds it after eikon.
