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
- **Inference:** data assimilation and uncertainty quantification, in the inference foundation.
- **Machine learning:** surrogate models.
- **Observation and field work.**
- **Visualization:** star and galaxy catalogues at full scale in the renderers.

## Phase 4. The Next Domains

- **Choose them,** with the user, by where their structures meet those already here (see [the idea on shared structures](ideas/2026-10-22-shared-structures-across-the-sciences.md)). Candidates:
  - **chemistry**, where physis's chemistry and thermodynamics would become a domain's foundation;
  - **the life sciences:** genomics and sequence analysis would meet morphe's formal structures, and ecology the climate;
  - **the cryosphere**, between climate, oceanography and geophysics.
- **Place each in the order** by what it may depend on, with its README, `module.json`, code, checks and notes.
- **A cross-index of structures:** for each shared structure, the places it appears in each domain.
