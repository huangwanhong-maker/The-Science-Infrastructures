# 0035. The Mathematical Foundations: a Map, an Order, and the First C++

- **Status:** Accepted, provisionally (the C++ conventions in §2.4 and §2.5 settle CLAUDE.md §8's open decisions 1 and 2 until Phase 0 confirms them)
- **Date:** 2026-10-23
- **Decided by:** while the spiking networks were being designed, the user asked that mathematics common to several sciences be separated into the foundations and implemented now where the networks need it; that the mathematical package structure be designed first, since the mathematics is broad (topology, geometry, category and type theory, differential equations, linear algebra, statistics, numerical computation and simulation) and must be well structured and comprehensive, though only what the networks need is written now; and that mature packages such as NumPy be reused where they already do the work well, with reimplementation in C++ where it allows better optimisation and, later, heterogeneous computing. Claude designed the map and the conventions and wrote this record.

## 1. Context

The three foundations were designed in October 2026 (calx's, morphe's and physis's ADRs 0002 and 0003), with planned modules and no code. The spiking networks of the neuroscience domain, the first C++ of the infrastructure, need random streams, matrix exponentials, integrators, sparse storage, parallel execution, special functions, quadrature, root finding and units. Written inside the domain, these would be found again in every domain that simulates; written in the foundations without a plan, they would be placed by convenience. Three areas the user named had no home in the plans: logic and type theory, combinatorics and graphs, and dynamical systems as objects in their own right.

## 2. Decision

1. **A map of the mathematics onto the libraries,** in [docs/mathematics.md](../../docs/mathematics.md): morphe for exact, symbolic and discrete structures; calx for numerical computation and the computation stack; physis for the continuum (units, geometry, analysis, dynamics, control, discretisation) and the physical law written in it; *inference* (named provisionally) for probability, statistics and learning; and the methods of the social sciences, planned above it. Every area has one owning module; an area with an exact and a numerical side lives on each side with its mode of computation (linear algebra exact in morphe, numerical in calx).
2. **Comprehensive in plan:** every class of the Mathematics Subject Classification 2020 has a home in the map (docs/mathematics.md §8), so that whatever a science asks for has a place before it is written.
3. **New modules in the plans:** `morphe/logic` (logic and type theory) and `morphe/combinatorics` (graphs and hypergraphs as C-sets); `calx/random`, `calx/signal`, `calx/optimisation` and `calx/verification`; `physis/analysis`, `physis/dynamics`, `physis/control`, `physis/computational_geometry`, `physis/quantum` and `physis/statistical_mechanics`; and the sixteen planned modules of inference, from `probability` to `learning`. Each foundation records its part (calx ADR 0004, morphe ADR 0003, physis ADR 0004).
4. **The order:** morphe and calx independent; physis on both; inference on physis, calx and morphe; the social methods on inference; the domains on all. Within each library, the module order of the map.
5. **The C++ of the libraries:** a module is `<library>/[<area>/]<module>/` with `include/<library>/<module>/`, a CMake target `<library>::<module>`, and its checks in `tests/`; namespaces follow the headers (`calx::random`, `physis::dynamics`), and a domain's C++ is `<domain>::<module>` (`neuroscience::snn`).
6. **The build:** CMake 3.20 or later and C++20, the infrastructure's root `CMakeLists.txt` composing the libraries in their order; on the development machine GCC in WSL, on the server GCC.
7. **Own and borrowed implementations:** interfaces are owned; behind them, mature licence-compatible libraries may be used (Eigen, OpenBLAS, pocketfft, SUNDIALS), and Python work uses NumPy and SciPy directly. Our own implementations are written where they serve several devices through calx's IR, exploit structure a general library cannot, give control a library does not (reproducible streams, declared precision), or are how a method is learnt; each is checked against a reference.
8. **Written now, only what the spiking networks need** (status `draft`): `calx/random`, `calx/linalg` (small dense matrices, LU solve, the exponential, compressed sparse rows), `calx/integrators` (Euler, Heun, RK4, Dormand–Prince 5(4)), `calx/special_functions` (erfcx), `calx/quadrature` (adaptive Gauss–Kronrod), `calx/solvers` (Brent), `calx/verification`, `calx/compute/scheduler` (a pool of host threads); `physis/units` and `physis/dynamics`. The rest of the map stays `planned`.

## 3. Alternatives Considered

- **One mathematics library.** Exact and floating-point computation differ in what an algorithm must state (termination against accuracy) and in what they depend on; the three libraries already separate them, and a fourth (inference) separates the stochastic.
- **The networks' mathematics inside the domain,** to be moved later. It would be moved by every domain that simulates; the user asked for it to be placed now.
- **Dynamical systems in calx.** calx knows arrays and functions, not systems with events and invariants; a hybrid system's guard and jump, and the structure of a flow, are mathematics of the continuum, and physis is where structures meet numbers.
- **Eigen for the small matrices.** The networks' matrices are of order below twenty, exponentiated once per population; a hundred lines with a checked Padé approximant suffice, and leave Eigen for when large dense problems arrive.
- **A borrowed test framework** (doctest, Catch2). The checks must report as the JavaScript harness does, with values against limits, so that one runner adds them up; `calx/verification` is that report, and small.

## 4. Consequences

- The foundations have code for the first time; calx, physis and the neuroscience domain gain CMake projects, and the infrastructure a root `CMakeLists.txt`.
- The runner builds and runs the C++ checks where a compiler is found (locally, or WSL's GCC on Windows), and skips them, saying so, where none is.
- morphe gains no code yet; its new modules are planned.
- CLAUDE.md §8's open decisions 1 and 2 are settled provisionally; Phase 0 confirms or replaces them.
