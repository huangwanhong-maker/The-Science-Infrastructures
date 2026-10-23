# 0036. Networks of Spiking Neurons, and the C++ Among the Checks

- **Status:** Accepted
- **Date:** 2026-10-23
- **Amends:** [0028](0028-neuroscience-and-the-brain-atlas.md) (the neuroscience domain gains C++)
- **Decided by:** the user asked for spiking neural networks next, with many models including the latest, abstract dynamics and states implemented by models, the high-performance foundation as the main effort, in the neuroscience domain but not the Brain Atlas, in C++ or Python rather than JavaScript, with original papers cited in the code and recent models surveyed; and, along the way, for quantum spiking models (Ising-like and open-system, GKSL) to be considered. Claude designed the module and wrote this record. The domain's own record is its [decision 0002](https://github.com/huangwanhong-maker/neuroscience/blob/main/development/decisions/0002-spiking-networks-in-cpp.md).

## 1. Context

The roadmap's Phase 4 named spiking networks first among the next domains' work, as a module of neuroscience shown in the Brain Atlas. The user moved its presentation out of the Brain Atlas and asked for its engine in a compiled language on foundations designed for the whole of mathematics (ADR 0035). The infrastructure's runner had built C++ only for the API's clients.

## 2. Decision

1. **`domains/neuroscience/snn/`**, C++20 on calx and physis, as the domain's decision 0002 describes: neurons as physis hybrid systems, exact propagators, a parallel engine deterministic across threads, plasticity evaluated lazily, JSON specifications with units, a command-line simulator, a benchmark, and checks against closed forms, published firing patterns and mean-field theory.
2. **C++ among the checks.** `development/tools/checks.mjs` builds a CMake project (here, or with WSL's GCC on Windows), runs its CTest programs and counts their lines as a suite's own (`cppChecks`); `run-tests.mjs` runs each foundation's C++ as a suite (`run-cpp.mjs`); a domain's suite runs its own C++; `check-alone.mjs` copies each domain and foundation whole, so that a domain's C++ is built alone with only the foundations beneath it, in a build directory removed afterwards.
3. **The next models** come from a verified survey kept as the domain's idea of 2026-10-23: GIF, exact AdLIF, linear multi-compartment and complex resonate-and-fire neurons, dendritic burst neurons, conductance-based interneurons, neural masses, neuromorphic neurons; e-prop, three-factor, voltage-based, inhibitory and calcium-based plasticity; and quantum neurons with exact classical routes (a qubit's Bloch flow with quantum-jump spikes, open quantum Hopfield networks, dissipative reservoirs), for which physis plans a `quantum` module (docs/mathematics.md).
4. **An application of its own** to explore networks, apart from the Brain Atlas, comes after the library.

## 3. Alternatives Considered

- **Showing the networks in the Brain Atlas**, as the roadmap first said: the user asked for an application apart.
- **Skipping the domain's C++ when it runs alone:** the alone run is what shows that a domain stands on what it declares; its C++ must be part of it.

## 4. Consequences

- The whole run builds three CMake projects (calx, physis, neuroscience) and one more for the domain alone; the first build takes some minutes, later ones reuse `~/.cache/science-build/`.
- The roadmap's Phase 4 records the networks as begun and their application as apart from the Brain Atlas.
