# The Mathematics of the Infrastructure

*Status: draft, 2026-10-23.* This document maps the mathematics the sciences need onto the infrastructure's foundation libraries, sets the order among them and among their modules, checks the map against the whole of mathematics, and fixes how the libraries' C++ is laid out, named, built and checked. The decision is [ADR 0035](../development/decisions/0035-the-mathematical-foundations.md); each foundation records its part in its own decision log. Only the modules the first simulations need are written (§11), and since [ADR 0039](../development/decisions/0039-complexes-symmetry-and-fields.md) those of complexes, symmetry and fields; the rest of the map is the plan, `planned` until a domain needs it.

## 1. Principles of the Division

The sciences share their mathematics, so it is written once, below the domains, in libraries that name no object in nature. Six principles divide it among them.

1. **By the kind of object and the mode of computation.** Exact, symbolic and discrete structures are one library, *morphe*, where an algorithm states when it terminates. Approximate computation on floating-point numbers and arrays, and the machinery that runs it on devices, is another, *calx*, where an algorithm states its accuracy and cost. The continuum, where structures meet numbers (geometry, analysis, dynamics) and the laws of physics written in it, is a third, *physis*. Chance (probability, statistics, inference, learning) is a fourth, *inference*, named provisionally.
2. **One home for each area.** Every area has exactly one module that owns it. An area with an exact and a numerical side lives on each side with its mode of computation, and each module names its counterpart: linear algebra over the integers in morphe and in floating point in calx; a finite group in morphe and a Lie group in physis.
3. **An order, and only downward dependence.** The libraries and their modules form a directed acyclic graph; a module uses only modules beneath it, through their public headers. morphe and calx are independent of each other.
4. **Generic, not particular.** A foundation knows no body, dataset or organism. A neuron is a hybrid dynamical system to physis, and becomes a neuron only in the neuroscience domain.
5. **Comprehensive in plan, written by need.** The plan covers the whole of mathematics as the Mathematics Subject Classification divides it (§8), so that whatever a science asks for has a place before it is written. A module is written when a domain or a sibling needs it; an abstraction is introduced once its pattern has appeared twice.
6. **Own the interfaces, borrow where it serves** (§10.3).

## 2. The Libraries and Their Order

```
   the domains and the works        may use every library beneath them
        │
        ▼
   social methods (planned)         measurement · econometrics · games · social choice
        │
        ▼
   inference (planned)              chance: probability · statistics · learning
        │
        ▼
   physis                           the continuum: units · geometry · analysis · dynamics · physical law
        │
        ▼
   morphe            calx           structures, exact  ·  computation, numerical   (independent of each other)
```

Each level may use every level beneath it, and none above.

| Library | Namespace | Kind of mathematics | Depends on | State |
|---|---|---|---|---|
| morphe | `morphe` | Structures and their equivalences: exact, symbolic, discrete | — | designed; first slice written (arithmetic, algebra, combinatorics, topology) |
| calx | `calx` | Numerical methods on floating-point numbers and arrays; the computation stack | — | designed; first modules written |
| physis | `physis` | Units, geometry, analysis, dynamics, control, discretisation, and physical law | calx, morphe | designed; units, dynamics, quantum, symmetry, geometry (in part), spacetime, field theory, gauge, quantum gravity written |
| inference (name open) | `inference` | Probability, stochastic processes, statistics, inference, learning | calx, morphe, physis | designed here (§6); no repository yet |
| social methods (name open) | — | Measurement, surveys, econometrics, games, social choice, agent-based models | inference, morphe, calx | named in README §4.3 (§7) |

*Remark.* inference depends on physis because data assimilation and the uncertainty of models need dynamical systems; physis never depends on inference. A stochastic differential equation is therefore a system of physis (`dynamics`), integrated by calx (`integrators`), whose statistics are inference's (`processes`).

*Remark.* physis has a mathematical half (geometry, analysis, dynamics, control, discretisation) and a physical half (mechanics to radiation). They share one library because the physics is written in that mathematics and they grow together; if the mathematical half grows large enough to stand alone, it can be extracted with its history under the graduation rule of ADR 0026.

## 3. Structures (morphe)

| Module | Content | Status |
|---|---|---|
| `logic` *(new)* | Propositional, first-order and modal logic; set theory as a foundation; computability; the simply and dependently typed λ-calculi, proof terms and normalisation (Curry–Howard–Lambek); satisfiability | planned |
| `categories` | Categories, functors, natural transformations; limits and colimits; adjunctions; monoidal and enriched categories; algebras and coalgebras of functors; toposes | planned |
| `equivalence` | Canonical forms; partition refinement; bisimulation; morphisms of coalgebras | planned |
| `arithmetic` | ℤ, ℚ, finite fields, algebraic numbers; polynomial rings; computational number theory (primes, factorisation, lattices of ℤⁿ) | **draft**: integers without bound, rationals, residues |
| `algebra` | Groups (finite, permutation, presented), rings, fields, modules, associative and Lie algebras over exact rings, representations; orders and lattices; universal algebra; homological algebra; **exact linear algebra** over ℤ, ℚ and finite fields (elimination, Hermite and Smith normal forms, kernels, determinants without rounding) | **draft**: the Smith normal form, homology with generators, presented groups, finite groups and homomorphisms into them |
| `combinatorics` *(new)* | Enumeration, generating functions, designs; **graphs** and hypergraphs as C-sets (paths, cuts, flows, matchings, colourings, isomorphism by refinement); matroids; exact combinatorial optimisation | **draft**: multigraphs, spanning forests |
| `topology` | General topology of finite spaces; simplicial, cubical and CW complexes; homology and cohomology; persistence; fundamental groups; homotopy equivalence; combinatorial manifolds and knots; dual complexes; coloured graphs (crystallisations) | **draft**: CW complexes with attaching words, simplicial complexes, constructions, homology, π₁, surfaces, manifolds, collapses, dual 2-complexes, coloured graphs, persistence over ℤ/2 of integer filtrations |
| `algebraic_geometry` | Polynomial ideals and Gröbner bases (Buchberger, F4); elimination; dimension, Hilbert series, free resolutions; varieties over ℚ and finite fields; toric varieties; real algebraic geometry; tropical geometry; schemes as far as computation reaches | planned (its order of work: morphe decision 0004, §2.3) |
| `formal` | Words and trees; automata, grammars and transducers; rewriting systems; transition systems and process calculi; bimorphisms | planned |

Order: `logic`, `categories` and `arithmetic` stand alone; `equivalence` uses `categories`; `algebra` uses `arithmetic`, `categories`, `equivalence`; `combinatorics` uses `categories`, `equivalence`; `topology` uses `algebra`, `combinatorics`, `equivalence`; `algebraic_geometry` uses `arithmetic`, `algebra`; `formal` uses `categories`, `equivalence`, `logic`.

## 4. Computation (calx)

| Module | Content | Status |
|---|---|---|
| `compute/*` | The representation (IR), schedules, compiler, runtime, task scheduler and backends of the computation stack; automatic differentiation as a transformation of the IR (calx's open decision 1) | `scheduler` **draft**: a pool of host threads for bulk-synchronous phases; the rest planned |
| `number_systems` | IEEE formats, double-double, interval arithmetic, compensated summation | planned |
| `linalg` | **Numerical linear algebra.** *Dense:* the kernels of BLAS, LU, Cholesky, QR, SVD, eigenproblems, least squares. *Sparse:* compressed rows and columns, sparse products, iterative solvers (conjugate gradients, GMRES) and preconditioners. *Functions of matrices:* exponential, logarithm, square root. *Structured:* banded, tridiagonal and tree systems; tensors | **draft**: small dense matrices, LU solve, the exponential, compressed sparse rows |
| `random` *(new)* | Counter-based generators; streams keyed by purpose and entity; sampling of the standard distributions; quasi-random sequences | **draft**: Philox4x32-10, uniform, normal, exponential, geometric, Poisson, bounded integers |
| `approximation` | Polynomial, Chebyshev, rational and spline approximation; interpolation in one and many dimensions; radial basis functions | planned |
| `special_functions` | Error and gamma functions and their relatives; orthogonal polynomials; Legendre, Bessel, hypergeometric and elliptic functions | **draft**: the scaled complementary error function |
| `transforms` | Fast Fourier transforms (real, complex, many-dimensional); discrete cosine and wavelet transforms; spherical-harmonic transforms | planned |
| `signal` *(new)* | Digital filters (FIR, IIR, zero-phase), windows, convolution and correlation, resampling, spectral estimation by periodograms | planned |
| `quadrature` | Integration on intervals (Gauss, Gauss–Kronrod, Clenshaw–Curtis), on the sphere, in many dimensions (cubature, Monte Carlo and quasi-Monte Carlo) | **draft**: adaptive Gauss–Kronrod (7, 15) |
| `differentiation` *(new)* | Forward differentiation by Taylor arithmetic (jets of order one and two) of functions written as templates; automatic differentiation as a transformation of the IR stays calx's open decision | **draft**: `Jet<N>`, `Dual<N>`, gradients and Hessians |
| `integrators` | Ordinary differential equations: explicit, embedded adaptive, implicit and IMEX; differential-algebraic and delay equations; stochastic differential equations (Euler–Maruyama, Milstein) | **draft**: Euler, Heun, RK4, Dormand–Prince 5(4) |
| `solvers` | Roots of functions; nonlinear systems (Newton, Broyden); fixed points; continuation | **draft**: bracketed roots (Brent) |
| `optimisation` *(new)* | Unconstrained (gradient, Newton, quasi-Newton, trust region); constrained (interior point, SQP); linear, quadratic and conic programming; stochastic and derivative-free methods; nonlinear least squares | planned |
| `verification` *(new)* | Comparison with references within declared limits; convergence orders; the report every C++ check prints | **draft** |

Order: every numerics module may use `compute/` for bulk kernels; `integrators`, `solvers` and `optimisation` may use `linalg`; `quadrature` may use `special_functions`; `signal` uses `transforms`; nothing in `compute/` uses `numerics/`.

## 5. The Continuum (physis)

| Module | Content | Status |
|---|---|---|
| `units` | Dimensions as the free abelian group ℤ⁷ over the SI base units; quantities read and written with their units; CODATA constants, versioned | **draft**: dimensions and quantities |
| `symmetry` *(new)* | Lie groups and algebras as physics uses them (U(1), SU(N), SO(N), Sp(N), Lorentz, Poincaré, SL(2, ℂ)); exponential maps, Haar measures, characters, the Weyl integration formula; representations, Clebsch–Gordan series and recoupling (3nj) symbols and their quantum deformations; symmetry breaking, Noether currents | **draft**: U(1), SU(2), SU(N), the Lorentz group and SL(2, ℂ), exact 3j, 6j, 9j, quantum 6j, Weyl integration |
| `geometry` | Euclidean, affine and projective geometry; spherical and ellipsoidal geometry; rotations (SO(3), SE(3), quaternions, from `symmetry`); manifolds, metrics, connections, geodesics, curvature, differential forms; information geometry; map projections | **draft**: metrics, curvature by forward differentiation, geodesics |
| `computational_geometry` *(new)* | Convex hulls, triangulations, Voronoi diagrams, meshes and their generation; spatial indices (k-d trees, R-trees); robust geometric predicates | planned |
| `analysis` *(new)* | Real and complex analysis as computation needs it (conformal maps, asymptotics); function spaces and operators, spectral theory; Fourier and harmonic analysis on groups and spheres (spherical harmonics, Wigner D-functions); integral transforms and integral equations; calculus of variations | planned |
| `dynamics` *(new)* | State spaces, vector fields and flows; linear systems and their exact propagators; hybrid systems (flow, guard, jump); discrete-time maps and difference equations; delay and stochastic dynamics; equilibria, stability, bifurcation, Lyapunov exponents, invariant measures | **draft**: linear propagators, hybrid systems |
| `control` *(new)* | Systems with inputs and outputs: controllability and observability, feedback, optimal control (LQR, Pontryagin), model-predictive control, system identification | planned |
| `discretisation` | Field equations in space and time: finite differences, volumes and elements; spectral methods; boundary elements; discrete exterior calculus on morphe's complexes; geometric and Lie-group integrators. A discretised field equation is a dynamical system on a finite state (the method of lines). | planned; a JavaScript stand-in |
| `quantum` *(new)* | Hilbert spaces, states and operators (pure and mixed); spins and qubits; unitary dynamics; open systems by the GKSL (Lindblad) master equation, vectorised as Liouvillians or in the affine coordinates of a generalised Bloch vector, integrated exactly as linear systems; quantum trajectories (Monte Carlo wave function: jumps at a hazard, as escape events); measurement; mean-field (product-state) closures | draft: the open two-level system and its trajectories |
| `statistical_mechanics` *(new)* | Ensembles and partition functions; lattice and spin models (Ising, Potts, Heisenberg; with a transverse field, quantum); Glauber and Metropolis dynamics as master equations; mean-field theory and phase transitions | planned |
| `spacetime` *(new)* | Spacetimes as objects: Minkowski space and causal structure; exact solutions in charts; horizons; discrete spacetimes (Regge geometry, causal sets, causal dynamical triangulations) | **draft**: Minkowski space, a catalogue of exact solutions, Regge geometry, causal sets, spectral and Hausdorff dimensions |
| `field_theory` *(new)* | Classical fields and Noether's theorem; lattice field theory on complexes; Monte Carlo and its statistics; perturbation theory (Feynman graphs, symmetry factors, integrals); free fields; renormalisation group; conformal and topological field theories; fields on curved spacetimes | **draft**: scalar fields on complexes, Markov-chain statistics, Feynman graphs, zero-dimensional field theory; invariant theories of every rank (scalar, vector, matrix, tensor), their exact series in N and large-N solutions |
| `gauge` *(new)* | Connections and curvature; lattice gauge theory on any 2-complex; Monte Carlo of gauge fields; Wilson and Polyakov loops; character and strong-coupling expansions; exact two-dimensional Yang–Mills; BF theory, flat connections; Chern–Simons | **draft**: gauge fields of finite groups, U(1), SU(2), SU(N) on CW complexes, updates, observables, Migdal–Witten |
| `quantum_gravity` *(new)* | Spin networks; spin foams (Ponzano–Regge, Turaev–Viro, Barrett–Crane, EPRL); group field theory and tensor models; matrix models and random surfaces; causal dynamical triangulations; causal-set dynamics | **draft**: Turaev–Viro state sums, BF amplitudes, coloured tensor models, matrix-model gluings; group field theories, spin foams with their quantum tetrahedra, random quadrangulations and melonic graphs, condensate cosmology |
| `electromagnetism` *(new)* | Maxwell's equations, potentials and gauges, radiation from sources, fields in media; on complexes by discrete exterior calculus | planned |
| `mechanics` · `gravitation` · `thermodynamics` · `continuum` · `chemistry` · `radiation` | The laws of physics, as physis's architecture plans them | planned |

Order: `units` stands alone; `symmetry` uses calx and morphe `arithmetic`; `geometry` uses `units` and `symmetry`; `spacetime` uses `geometry` and morphe `topology`; `field_theory` uses `spacetime`; `gauge` uses `field_theory` and `symmetry`; `quantum_gravity` uses `gauge`; `electromagnetism` uses `geometry`, `discretisation` and `gauge` (physis decision 0005, §2.2); `analysis` and `dynamics` use `units`; `computational_geometry` uses `geometry`; `control` uses `dynamics`; `discretisation` uses `geometry`, `analysis` and `dynamics`; `quantum` uses `geometry` (its groups) and `dynamics` (an open system's flow is linear on its density matrix); `statistical_mechanics` uses `thermodynamics`, `dynamics` and `quantum`; the physics modules as physis's architecture lists them, `mechanics` using `dynamics` too.

## 6. Chance (inference, planned)

The library of probability, statistics and learning, with its modules in the order they build on each other. It reads what calx samples and what physis simulates.

| Module | Content |
|---|---|
| `probability` | Probability spaces and measures; distributions as objects (densities, distribution functions, quantiles, moments, characteristic functions); transformations, mixtures and copulas; limit theorems and concentration inequalities |
| `information` | Entropy, mutual information, divergences; coding and channel capacity; the information geometry of statistical models (with physis `geometry`) |
| `processes` | Markov chains in discrete and continuous time; point processes (Poisson, renewal, Hawkes); Gaussian processes and random fields; diffusions and the statistics of stochastic differential equations; random graphs |
| `descriptive` | Summaries, robust statistics, empirical distributions, kernel density estimation, visual diagnostics' numbers |
| `estimation` | Likelihood, M-estimators, the method of moments; the bootstrap and other resampling; confidence sets |
| `testing` | Hypothesis tests, parametric and permutation; multiple comparisons (family-wise error, false discovery rate); power and sample size |
| `regression` | Linear, generalised linear and mixed models; nonparametric and semiparametric regression; regularisation (ridge, lasso) |
| `multivariate` | Principal components, factor analysis, canonical correlation, independent components; clustering and mixture models; dimension reduction |
| `time_series` | Autoregressive and moving-average models; spectral analysis; state-space models; change points; nonstationarity |
| `spatial` | Geostatistics and kriging; spatial point patterns; areal models; spatio-temporal statistics; statistics on the sphere and on manifolds |
| `bayes` | Bayesian inference: conjugate models, Markov chain Monte Carlo (Metropolis–Hastings, Gibbs, Hamiltonian), sequential Monte Carlo, variational inference, model comparison |
| `assimilation` | Data assimilation: Kalman filters (extended, ensemble), smoothers, variational assimilation (3D-Var, 4D-Var), particle filters |
| `uncertainty` | Propagation of uncertainty through models; sensitivity analysis (Sobol); surrogate models and emulators; calibration of models |
| `causal` | Causal graphs and the do-calculus; potential outcomes; identification and estimation of effects |
| `design` | Design of experiments; survey sampling and weighting; optimal design |
| `learning` | Statistical learning: kernels, trees and ensembles, neural networks and their training, reinforcement learning; generalisation and validation |

Order: `probability` stands alone; `information` and `processes` use it; `descriptive`, `estimation` and `testing` use `probability`; `regression`, `multivariate` and `time_series` use `estimation`; `spatial` uses `regression` and `processes`; `bayes` uses `estimation` and `processes`; `assimilation` and `uncertainty` use `bayes` and physis `dynamics`; `causal` uses `regression` and morphe `combinatorics`; `design` uses `estimation`; `learning` uses `regression`, `multivariate` and calx `optimisation`.

*Remark.* Its name and first module are an open decision (CLAUDE.md §8, open decision 3). The first module will be the one a domain needs first; for the neuroscience domain that is `processes` (spike trains as point processes) and `testing` (statistics of maps).

## 7. The Methods of the Social Sciences (planned)

Measurement (psychometrics, item response theory, indices), surveys, econometrics (panels, instruments, discrete choice), game theory and mechanism design, social choice and voting, social networks, and agent-based models. They rest on inference (estimation, causal inference, design), morphe (games and choice as combinatorial and order structures, networks as graphs) and calx (optimisation, simulation). Whether they form a library of their own or modules of inference is open (README §4.3).

## 8. Coverage: the Mathematics Subject Classification

Every class of the Mathematics Subject Classification 2020 (zbMATH Open and Mathematical Reviews) and its home. A class may be split between homes by its exact and numerical sides; a class of applications lives in the domains, which use the foundations.

| MSC 2020 | Class | Home |
|---|---|---|
| 00, 01, 97 | General topics; history and biography; mathematics education | the notes and theory of every module; not a library |
| 03 | Mathematical logic and foundations | morphe `logic` |
| 05 | Combinatorics | morphe `combinatorics` |
| 06 | Order, lattices, ordered algebraic structures | morphe `algebra` |
| 08 | General algebraic systems | morphe `algebra` (universal algebra), `categories` |
| 11 | Number theory | morphe `arithmetic` |
| 12, 13 | Field theory and polynomials; commutative algebra | morphe `algebra`, `algebraic_geometry` |
| 14 | Algebraic geometry | morphe `algebraic_geometry` |
| 15 | Linear and multilinear algebra; matrix theory | morphe `algebra` (exact); calx `linalg` (numerical) |
| 16, 17 | Associative rings and algebras; nonassociative rings and algebras | morphe `algebra`; Lie algebras of physical symmetries in physis `symmetry` |
| 18 | Category theory; homological algebra | morphe `categories`, `algebra` |
| 19 | K-theory | morphe `algebra`, `topology` (research only) |
| 20 | Group theory and generalisations | morphe `algebra` |
| 22 | Topological groups, Lie groups | physis `symmetry` |
| 26 | Real functions | physis `analysis`; calx `approximation` |
| 28 | Measure and integration | inference `probability`; calx `quadrature` |
| 30, 32 | Functions of one and several complex variables | physis `analysis`; calx `special_functions` |
| 31 | Potential theory | physis `analysis`, `gravitation` |
| 33 | Special functions | calx `special_functions` |
| 34 | Ordinary differential equations | physis `dynamics` (systems); calx `integrators` (methods) |
| 35 | Partial differential equations | physis `discretisation`, `continuum` |
| 37 | Dynamical systems and ergodic theory | physis `dynamics` |
| 39 | Difference and functional equations | physis `dynamics` |
| 40 | Sequences, series, summability | calx `approximation`, `number_systems` (acceleration, compensated sums) |
| 41 | Approximations and expansions | calx `approximation` |
| 42, 43 | Harmonic analysis on Euclidean spaces; abstract harmonic analysis | physis `analysis`; calx `transforms` |
| 44 | Integral transforms, operational calculus | physis `analysis`; calx `transforms` |
| 45 | Integral equations | physis `analysis`, `discretisation` (boundary elements) |
| 46, 47 | Functional analysis; operator theory | physis `analysis` |
| 49 | Calculus of variations and optimal control; optimisation | physis `analysis` (variations), `control`; calx `optimisation` |
| 51 | Geometry | physis `geometry` |
| 52 | Convex and discrete geometry | physis `computational_geometry`; morphe `combinatorics` (discrete) |
| 53 | Differential geometry | physis `geometry` |
| 54 | General topology | morphe `topology` |
| 55 | Algebraic topology | morphe `topology` |
| 57 | Manifolds and cell complexes | morphe `topology` (combinatorial); physis `geometry` (smooth) |
| 58 | Global analysis, analysis on manifolds | physis `geometry`, `analysis` |
| 60 | Probability theory and stochastic processes | inference `probability`, `processes`; calx `random` (sampling) |
| 62 | Statistics | inference (all its modules but `learning`) |
| 65 | Numerical analysis | calx `numerics` (all modules); physis `discretisation` (65M, 65N) |
| 68 | Computer science | calx `compute` (algorithms, parallel computation); morphe `formal`, `logic` (languages, semantics, complexity); inference `learning` (68T) |
| 70 | Mechanics of particles and systems | physis `mechanics` |
| 74 | Mechanics of deformable solids | physis `continuum` |
| 76 | Fluid mechanics | physis `continuum` |
| 78 | Optics, electromagnetic theory | physis `radiation`, `electromagnetism` (planned) |
| 80 | Classical thermodynamics, heat transfer | physis `thermodynamics` |
| 81 | Quantum theory | physis `quantum`; `field_theory`, `gauge` (81T, quantum field theory) |
| 82 | Statistical mechanics, structure of matter | physis `statistical_mechanics`, `thermodynamics` |
| 83 | Relativity and gravitational theory | physis `spacetime`, `gravitation`; `quantum_gravity` (83C45) |
| 85, 86 | Astronomy and astrophysics; geophysics | the domains (astronomy; geophysics, oceanography, meteorology) |
| 90 | Operations research, mathematical programming | calx `optimisation`; morphe `combinatorics`; inference `processes` (queues) |
| 91 | Game theory, economics, social and behavioural sciences | the social methods (§7); the social-science domains |
| 92 | Biology and other natural sciences | the domains (neuroscience, the life sciences, chemistry) |
| 93 | Systems theory; control | physis `control`, `dynamics`; inference `assimilation` (filtering) |
| 94 | Information and communication theory, circuits | inference `information`; calx `signal` |

## 9. How the Areas Meet

- **Differential equations** are owned by what is asked of them. A method for x' = f(x) on vectors, knowing nothing of f, is `calx/integrators`. A system as an object (its state space, its flow, its events, its invariants) is `physis/dynamics`, which integrates by calx. A field equation in space and time is discretised by `physis/discretisation` into a dynamical system. A stochastic one is a system of `physis/dynamics` whose statistics are `inference/processes`.
- **Linear algebra** is exact in morphe (over ℤ, ℚ and finite fields: the Smith normal form gives homology) and numerical in calx (in floating point, dense and sparse, on any device). A Lie algebra, linear algebra with a bracket acting on space, is physis's geometry. A linear system x' = A x + B u + c is integrated exactly by the exponential of an augmented matrix (Van Loan 1978): calx computes the exponential, physis builds the propagator.
- **Topology** is exact in morphe (homology over ℤ) and meets numbers in persistence, whose filtrations are real; geometry (metric, curvature) is physis's.
- **Groups** are morphe's when finite or presented, physis's when they are Lie groups acting on space.
- **Graphs** are morphe's as structures (`combinatorics`), calx's as sparse matrices (`linalg`), and inference's as random models (`processes`). A network of neurons is stored by rows, simulated, and analysed as a graph.
- **Randomness** is sampled by calx (`random`), described and inferred by inference.
- **Optimisation** is calx's when continuous, morphe's when combinatorial and exact, physis's when it is the control of a dynamical system, inference's when it is the fitting of a model.
- **Signals** are filtered and transformed by calx (`signal`, `transforms`) and estimated by inference (`time_series`).
- **Quantum neurons** meet the classical machinery exactly: an open quantum system's density matrix (or a qubit's Bloch vector) flows linearly under the GKSL equation, so `physis/quantum` hands `physis/dynamics` a linear system whose propagator `calx/linalg` exponentiates; a quantum-jump unravelling fires at a hazard, as an escape-noise neuron does; Ising-like networks are `statistical_mechanics`' spin models with Glauber or Lindblad dynamics. The neuroscience domain's spiking networks have the first, a two-level emitter along exact quantum trajectories, and plan more (its survey of 2026-10-23).
- **Complexes and fields** meet in physis: a scalar field is a 0-cochain on a morphe complex; a lattice gauge field a function from its 1-cells to a group of `symmetry`, whose holonomy around a 2-cell is the product along that cell's attaching word; the Feynman graphs of a group field theory are morphe's coloured graphs, whose bicoloured cycles bound the 2-cells of the dual 2-skeleton on which their amplitudes are BF or spin-foam state sums; a discrete spacetime is a simplicial complex with edge lengths (Regge) or a partial order (causal sets). The exact side (homology, fundamental groups, counting homomorphisms) is morphe's; the numerical side (Monte Carlo, recoupling in floating point, curvature) physis's (ADR 0039).
- **Spacetime from field theories** ([ADR 0040](../development/decisions/0040-the-generation-of-spacetime.md), physis decision 0006): an action's Feynman graphs are morphe's coloured graphs, glued from bubbles; their 2-complexes and pseudomanifolds are the discrete spacetimes, spin foams when the theory is a group field theory; their large-N or condensed phase (`field_theory/large_n`, `quantum_gravity/condensate`) and its continuum (critical exponents, `spacetime/dimensions`) are physis's. Scalar, vector, matrix and tensor theories and group field theories run through the same five stages.
- **Simulation** is the composition: a domain states its systems through `physis/dynamics`, integrates them with `calx/integrators` and `calx/linalg`, draws its noise from `calx/random`, runs in parallel on `calx/compute`, and is compared with data through inference.

## 10. The C++ of the Libraries

### 10.1 Layout and Names

A module of a library is a directory with its README, public headers, implementation and checks (architecture §6.2):

```
<library>/<area>/<module>/            calx's areas are compute/ and numerics/; the others have none
├── README.md
├── CMakeLists.txt                    the target <library>::<module>
├── include/<library>/<module>/       #include <calx/random/philox.hpp>
├── src/                              when the module is compiled
├── tests/                            checks on calx/verification
├── bindings/python/                  later
└── benchmarks/                       when speed is a claim
```

- **Namespaces** follow the headers: `calx::random`, `calx::linalg`, `physis::units`, `physis::dynamics`. A domain's C++ is `<domain>::<module>`, so the spiking networks are `neuroscience::snn` (provisionally; CLAUDE.md §8, open decision 1).
- **Header-only first.** Templates and small inline functions are header-only `INTERFACE` targets; a module gains `src/` when it has code worth compiling once.
- A module's directory is created when its first content is written; until then it is a `planned` row in its library's README.

### 10.2 Build

- **CMake 3.20 or later, C++20,** every library a project of its own, every module a target. Inside the infrastructure, the root `CMakeLists.txt` adds the libraries in their order, and each later library uses the targets of the earlier ones that already exist; built alone, a library would fetch its siblings at pinned versions (not yet written).
- **Compilers.** GCC 11 or later, Clang 14 or later; MSVC when the Windows toolchain is chosen. On the development machine the compiler is GCC in WSL (CLAUDE.md §8, open decision 2, decided provisionally).

### 10.3 Own and Borrowed Implementations

The libraries own their interfaces; behind an interface an implementation may be borrowed (CLAUDE.md §2.2).

- **Borrow** a mature, licence-compatible library where it already does the work well and the infrastructure would gain nothing by redoing it: Eigen (MPL 2.0) or a BLAS and LAPACK (OpenBLAS, BSD) for large dense linear algebra on the CPU; pocketfft (BSD) for transforms; SUNDIALS (BSD) for stiff and differential-algebraic systems. In Python, the orchestration, notebooks and data library use NumPy and SciPy directly and do not wrap them; statistical work in Python may use SciPy and statsmodels until inference has its own.
- **Write our own** where the infrastructure gains by it:
  - kernels that must run on more than one device, which calx expresses in its IR so that one definition serves the CPU and the accelerators to come (calx architecture §7);
  - structure that a general library cannot exploit: small matrices of fixed order, event-driven sparse updates, a propagator computed once and applied to millions of states;
  - control a library does not give: reproducible random streams independent of the number of threads, declared precision never silently reduced;
  - and, this being a learning project, where writing the method is how it is understood; such an implementation is checked against a borrowed one or a published reference.
- **Copyleft** libraries (FFTW, SHTns, FLINT, GMP) need an explicit decision and are never linked by default.

For the spiking networks every need is of the second or third kind, and the first implementations are calx's own.

### 10.4 Checks

- Each module's checks are a C++ program on `calx/verification`, printing what the JavaScript harness prints: `ok` or `FAIL`, the name, the value and its limit, then the count. The infrastructure's runner builds every library with CMake, locally or in WSL, runs the programs, and adds them to its total.
- **Reference values** live with the checks that use them, each with its source.
- **Determinism.** A computation with random input is reproducible from its seed for any number of threads (`calx/random`'s streams are keyed by what they are for, not by who draws them).

## 11. What Is Built First

The spiking networks of the neuroscience domain are the first C++ of the infrastructure. They need, from the foundations:

| Need | Module |
|---|---|
| Random input and connectivity, reproducible across threads | `calx/random` |
| Exact propagators of linear models | `calx/linalg` (the exponential), `physis/dynamics` |
| Nonlinear models (adaptive exponential, Hodgkin–Huxley) | `calx/integrators`, `physis/dynamics` |
| Synapses stored by rows | `calx/linalg` (compressed sparse rows) |
| Parallel update and delivery | `calx/compute/scheduler` |
| Mean-field theory (Siegert's formula, fixed points) | `calx/special_functions`, `calx/quadrature`, `calx/solvers` |
| Parameters with units | `physis/units` |
| Checks | `calx/verification` |

The second slice ([ADR 0039](../development/decisions/0039-complexes-symmetry-and-fields.md)) is the research on complexes, symmetry and fields: morphe's `arithmetic`, `algebra`, `combinatorics` and `topology`; calx's `differentiation`; physis's `symmetry`, `geometry` (metrics and curvature), `spacetime`, `field_theory`, `gauge` and `quantum_gravity`; and the Complex Atlas, which shows them.

Everything else in this map is planned. The roadmap continues with morphe's arithmetic, algebra and topology toward discrete exterior calculus, and physis's shallow-water model (README §4.6); inference begins when a domain first needs its statistics.

## 12. Open Questions

1. **The names of inference and of the social methods,** and the boundary between them (CLAUDE.md §8).
2. **Whether the mathematical half of physis** (geometry, analysis, dynamics, control, discretisation) becomes a library of its own as it grows (§2).
3. **Typed units.** Dimensions are checked at run time now; compile-time types, and calx's typed IR annotations, may follow (physis's open decision 1).
4. **Automatic differentiation** in calx's IR, which `optimisation`, `bayes` and `learning` would use (calx's open decision 1).
