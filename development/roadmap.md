# Roadmap

*Status: draft, 2026-10-22.* The phases are ordered intentions, not commitments. Revise them as the work teaches us.

## Progress

- **2026-10-22. The first applications, ahead of Phase 1.** At the user's request, the Celestial Atlas and then the Earth Atlas were built before Phase 0 ([ADR 0011](decisions/0011-first-application-celestial-atlas.md), [ADR 0014](decisions/0014-earth-atlas-separated-data-and-typed-layers.md)). Their provisional JavaScript already covers Phase 1 steps 1–3:
  - time scales;
  - reference frames;
  - the evaluator for the JPL DE format and the registry entry for DE440.

  It is validated against JPL Horizons. When the C++ is built, it takes over this work and must pass the same checks.
- **2026-10-22. The umbrella.** The Science Infrastructures was founded, with the first laboratory as a submodule ([ADR 0025](decisions/0025-an-umbrella-of-laboratories.md)).
- **2026-10-22. Domains as libraries** ([ADR 0026](decisions/0026-domains-as-libraries.md)):
  - The first laboratory's history was joined to this repository, and the libraries were raised to `foundations/`.
  - Its knowledge became the domains: one library per science, each declaring its dependencies, each passing its checks alone.
  - Geophysics is a domain.
  - The word "laboratory" is retired.
- **2026-10-23. The Earth Atlas as a forecasting system, and an API** ([ADR 0030](decisions/0030-the-api-edge.md), [ADR 0031](decisions/0031-observations-scores-and-hazards.md)): observations (METAR, buoys, the user's sensors), forecasts kept and scored on a schedule, ECMWF's open forecasts, hazard models (fire weather, heat, heavy rain, tsunami travel times, aftershocks, volcanic ash), people and regions with exposure, drainage, Sentinel-2, assimilation of real observations; and `api/`, through which programs in Python, C and C++ and JavaScript reach the data, the domains' functions, the sensors in real time and the plug-ins. The research arc of prediction is planned, and deferred.
- **2026-10-23. The order of the next work,** set by the user: first the foundations, morphe (the algebraic foundation: category theory, algebra, topology, geometry) and physis (theoretical physics); then five domains (Phase 4). The README's §4 maps every science, the social sciences and the humanities included, and plans a foundation for the quantitative social sciences.
- **2026-10-22. Neuroscience, the first domain beyond the Earth and the sky** ([ADR 0028](decisions/0028-neuroscience-and-the-brain-atlas.md)), in a repository of its own as a submodule: a brain atlas in the sense of the Human Brain Project's (reference spaces, parcellations, maps, for the human, the mouse and the rat), and the head with the 10–20 electrode system; the Brain Atlas application draws it and maps data onto it. Phase 4's choice of the next domains has begun.

## Phase 0. Foundations

- **Repositories.** Done (2026-10-22): this repository and the three foundations are published on GitHub.
- **Names.**
  - Settled: The Science Infrastructures; calx, morphe and physis.
  - Still open:
    - the C++ namespaces of the domains;
    - the name of the inference foundation.
- **C++ toolchain.** As of 2026-10-22, the Windows side has no C++ compiler and no CMake. The options are the MSVC Build Tools, LLVM/clang, or WSL.
- **Build.**
  - A CMake superbuild across `foundations/` and `domains/`, in which each domain is a target that links only its down-set.
  - A C++ test framework.
  - Python bindings.
  - pytest.
- **Continuous integration** in every repository. Here it runs `node development/tools/run-tests.mjs`: the dependency rules, every suite, and each domain alone.
- **Conventions.** A glossary, with the orientation and normalisation conventions; a `references.bib`. Each foundation keeps its own.

## Phase 1. First Vertical Slice: Sun, Moon, Earth and the Tides

This slice touches every tier at small scale. It reaches astronomy, geodesy, oceanography and the tides bridge, and atmospheric tides are a natural extension.

| Step | Where | Work |
|---|---|---|
| 1 | physis `units`; `domains/metrology` | Units; the time scales UTC, TAI, TT, TDB and UT1; frames for the Sun and the Moon |
| 2 | calx `numerics`; physis `geometry` | Chebyshev evaluation and associated Legendre functions; rotations |
| 3 | `domains/astronomy`, `data/` | An evaluator for the JPL DE format; a registry entry for DE440 |
| 4 | physis `gravitation`; `domains/tides` | Spherical-harmonic potential theory; the tide-generating potential; the equilibrium tide; Doodson arguments |
| 5 | calx `compute` | A minimal IR and CPU reference interpreter. The first IR kernel is the equilibrium tide on a global grid. |
| 6 | the inference foundation; `data/` | Least-squares harmonic analysis of a real NOAA tide-gauge record |
| 7 | `visualization/` | Scene model v0, a static renderer and a web renderer, showing the Earth–Moon system with the tidal bulge |
| 8 | the domains' `tests/` | The ephemeris checked against JPL Horizons; the tides checked against NOAA predictions |
| 9 | the domains' and foundations' notes | A note for each derivation used |

## Phase 2. Breadth and a Second Device

- **calx: a second backend.** Either CUDA (on the RTX 3050 Ti) or SYCL with Level Zero (on the Iris Xe). This tests whether the stack really is retargetable.
- **morphe: two slices.**
  - Arithmetic, then algebra (Smith normal form), then topology (homology of cell complexes). This is the foundation for discrete exterior calculus in physis.
  - A research slice: partition refinement and bisimulation.
- **physis: a shallow-water model on the rotating sphere**, in C++, replacing the stand-in. It serves both ocean tides and atmospheric waves, and runs on both backends.
- **The domains.**
  - An energy-balance model driven by a Milankovitch insolation calculation (climate, astronomy).
  - New applications: an orbital sandbox, a digital tide-predicting machine, and a sight-reduction trainer.

## Phase 3. Research Capacity

- **Ocean dynamics beyond tides:** circulation, waves, seawater thermodynamics and biogeochemistry.
- **Meteorology:** a dry dynamical core first ([the idea of a primitive-equation model](ideas/2026-10-22-primitive-equation-model.md)), then moist processes and radiation.
- **Geophysics:** seismology, gravity and magnetism; slabs where they can be fetched.
- **Geography:** terrain and drainage analysis, remote sensing, and a land-use atlas.
- **Coupling:** the atmosphere and ocean coupled through `air-sea`, leading towards a simple coupled climate model.
- **Inference:** data assimilation and uncertainty quantification, in the inference foundation; the steps are planned in [the research arc of prediction](ideas/2026-10-23-the-research-arc.md).
- **Machine learning:** surrogate models.
- **Observation and field work.**
- **Visualization:** star and galaxy catalogues at full scale in the renderers.

## Phase 4. The Next Domains

- **After the foundations,** morphe and physis (README §4.6), which the user set first. Then **in the order the user set**, each reusing what the one before made, chosen where their structures meet those already here (see [the idea on shared structures](ideas/2026-10-22-shared-structures-across-the-sciences.md)). Neuroscience came first (ADR 0028, ADR 0032).
  1. **Spiking neural networks,** a module of neuroscience in C++ on the foundations (`snn`, ADR 0036): *begun 2026-10-23* with fourteen neuron models, plasticity, a parallel engine and its mean-field checks. Next, from the [survey](https://github.com/huangwanhong-maker/neuroscience/blob/main/development/ideas/2026-10-23-a-survey-of-spiking-models.md): recent neuron models (twenty-two written by the end of the day, the first quantum neuron among them), online learning rules (e-prop, three-factor), more quantum neurons (open quantum Hopfield networks), Python bindings, and an application of their own, apart from the Brain Atlas. Later, and only planned: networks compiled through calx to GPUs, neuromorphic chips and quantum processors ([idea](https://github.com/huangwanhong-maker/neuroscience/blob/main/development/ideas/2026-10-23-spiking-networks-on-devices.md)).
  2. **Computational linguistics:** automata, context-free grammars, parsing, grammar inference, on morphe's formal structures.
  3. **Genomics:** sequence alignment, k-mers, assembly, variant analysis, reusing the sequence, grammar, graph and inference algorithms.
  4. **Molecular science:** molecular graphs, conformers, reaction networks, on physis's chemistry and thermodynamics.
  5. **Drug discovery,** a bridge over genomics and molecular science: target–compound evidence, molecular properties, pharmacological models.
- **A foundation for the quantitative social sciences** (README §4.3): measurement, surveys, causal inference, panels, networks, games and social choice, agent-based models; its name and its boundary with the inference foundation to be decided.
- **Then tiers 2 and 3** (README §4.6, proposed): the inference foundation and the bridges tides, air–sea and navigation, as Phases 1 and 3 plan them; then chemistry, evolution, cell biology, the cryosphere, geology, planetary science, ecology, the quantitative social methods, demography and epidemiology.
- **The rest of the map** (README §4.4, tier 4), the social sciences and the humanities included, as the work reaches each.
- **Place each in the order** by what it may depend on, with its README, `module.json`, code, checks and notes.
- **A cross-index of structures:** for each shared structure, the places it appears in each domain.
