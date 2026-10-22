# Architecture

*Status: draft, 2026-10-21.* This document describes how the laboratory and its libraries are structured: the repositories, the layers, the rules that connect them, and the representation stack. The decisions behind it are recorded in [development/decisions/](../development/decisions/).

## 1. Principles

### 1.1 Pure Core, Effectful Edges

The libraries in `lib/` and the laboratory's `core/` perform pure computation: given inputs, they return results. They have no network access, no instrument I/O, no user interface and no hidden global state. They read only files whose paths the caller passes in.

Everything that touches the outside world lives at the *edges* and calls the core. The edges are `data/`, `observation/`, `visualization/` and `applications/`.

*Rationale.* A pure core can be tested exhaustively and reproduced exactly. It can be deployed to any device, including devices without an operating system or a file system.

### 1.2 Layered, Acyclic Dependencies

Each layer depends only on the layers beneath it (§3). Within a layer or a library, a module may depend on another when two conditions hold:

- the dependency is declared;
- the module graph stays acyclic.

### 1.3 Separation of Content from Realisation

Two parts of the project face a set of technologies that is open-ended and keeps changing:

- the hardware that executes computations;
- the media that present results.

In both, a stable, device-neutral intermediate representation separates *what* from *how*, and each technology is a backend that can be added without modifying anything else. For computation this is the stack in calx (§4). For presentation it is the representation layer (§5).

### 1.4 Own the Interfaces, Borrow the Implementations

The project defines its own small, stable interfaces. Mature external libraries may sit behind them and be replaced later. Every external dependency must have a licence compatible with MIT, or be adopted by an explicit decision.

### 1.5 Independent Libraries

General knowledge lives in independent libraries, each in its own repository ([ADR 0008](../development/decisions/0008-independent-libraries.md)):

- numerics and computation;
- mathematical structures;
- geometry and physics.

The laboratory composes these libraries and adds what is specific to its subject.

### 1.6 Growth by Vertical Slices

Interfaces are extracted from working code that runs end to end. They are not designed completely in advance. Each slice passes through every layer at small scale ([roadmap](../development/roadmap.md)).

## 2. Repositories and Layer Map

```
Edges        applications/    labs and tools, by field and bridge
             visualization/   scene model · renderers · styles
             observation/     instruments · field work
             data/            registry · acquisition · formats · provenance
                    │ call into
                    ▼
core/        inference        estimation · filtering · assimilation · uncertainty · learning
             domains          fields: astronomy · geography · oceanography · meteorology · geodesy
                              bridges: tides · navigation · air_sea · climate
             metrology        time scales · reference frames · Earth orientation
                    │ use
                    ▼
lib/         physis           units · geometry · discretisation · mechanics · gravitation ·
                              thermodynamics · continuum · chemistry · radiation
             calx             compute stack · numerical methods
             morphe           arithmetic · categories · equivalence · algebra · topology ·
                              algebraic geometry · formal structures

Knowledge    theory/ · notebooks/ · research/
Process      development/   (ideas · decisions · roadmap)
Assurance    validation/
```

The libraries depend on each other as follows:

```
laboratory ──► physis ──► calx
     │            └─────► morphe
     └──► calx, morphe        (direct use is allowed)

calx and morphe do not depend on each other.
```

## 3. Dependency Rules

1. **Libraries.**
   - calx and morphe depend on no project library.
   - physis depends on calx and morphe.
   - No library depends on the laboratory.
2. **Inside `core/`.** The layers are ordered metrology < domains < inference. Each layer depends only on earlier layers and on `lib/`.
3. **Within a layer or a library.** Modules may depend on each other when the dependency is declared and the graph stays acyclic.
4. **The core stands alone.** Nothing in `core/` depends on an edge, knowledge, process or assurance directory.
5. **Edges.**
   - `data/` and `observation/` depend only on `core/` and `lib/`.
   - `visualization/` depends on those, and may also depend on `data/`.
   - `applications/` may depend on everything.
6. **Who may depend on what.**
   - Only `applications/`, `notebooks/` and `research/` may use `visualization/`.
   - Nothing depends on `applications/`.
7. **Leaf directories.** `notebooks/`, `research/` and `validation/` may depend on anything, and nothing depends on them.

A violation of these rules is a design error. Fix it in one of two ways:

- move the code down;
- or define an interface in the lower layer that the upper layer implements (dependency inversion).

Never fix it with an upward import. When the right place for the code is a library, make the change in that library's repository.

*Example.* Tidal computations need lunar and solar ephemerides, but the core never fetches files.

- `core/domains/astronomy` defines an ephemeris interface and an evaluator for the JPL DE format.
- `data/` locates, verifies and records the DE440 file, and passes its path to the evaluator.

## 4. Computation Stack (calx)

calx keeps the parts of its computation stack separate:

- representation (`ir/`);
- schedule (`schedule/`);
- compiler (`compiler/`);
- runtime (`runtime/`);
- task scheduler (`scheduler/`);
- backends (`backends/`).

The stack follows four rules:

- A new device means a new backend directory, and nothing else changes.
- The CPU reference interpreter defines correct results.
- Every computation declares its precision.
- Precision is never silently reduced.

The full design is in calx's own [architecture document](https://github.com/huangwanhong-maker/calx/blob/main/docs/architecture.md); locally it is at `lib/calx/docs/architecture.md`.

The laboratory writes bulk kernels through calx, and never addresses a device directly.

## 5. Representation Stack (`visualization/`)

The representation layer applies the same separation of content from realisation to output. Only `applications/`, `notebooks/` and `research/` may use it. Nothing in `core/`, `lib/`, `data/` or `observation/` depends on it.

| Component | Responsibility |
|---|---|
| `scene/` | A device-neutral, serialisable description of what is shown: bodies, fields, trajectories, reference frames, annotations, cameras, time. It carries units, time scales and frames from the core, so renderers never guess them. |
| `renderers/` | One per technique. Examples: static figures (Matplotlib, SVG, PDF); web (WebGPU, three.js, CesiumJS); native (Vulkan); immersive (WebXR); physical (3D printing); non-visual (sonification). |
| `styles/` | Palette, typography and marks, shared by all renderers |

*Remark.* Astronomical scenes span roughly 10⁻³ m to 10²⁶ m, far beyond the range of single-precision graphics coordinates. The scene model therefore uses a hierarchy of local frames with floating origins, not one global coordinate system, as OpenSpace does.

## 6. Module Layout

A C++ module, in any of the four repositories, has this layout:

```
<module>/
├── README.md                     purpose, status, dependencies, references, open questions
├── include/<namespace>/<module>/ public headers (calx, morphe, physis; the laboratory's is not yet named)
├── src/                          implementation
├── tests/                        unit tests
├── bindings/python/              Python bindings (optional)
└── benchmarks/                   performance measurements (optional)
```

Until a module has content, it exists only as a `planned` row in its parent README.

## 7. Placement Guide

| Content | Location |
|---|---|
| A kernel that must run on several devices | IR through calx (`lib/calx/compute/`) |
| A new device or accelerator | `lib/calx/compute/backends/<target>/` |
| Linear algebra, FFT, spherical-harmonic transform, ODE integrator | `lib/calx/numerics/` |
| Exact arithmetic, groups, representations, lattices, presentations | morphe: `arithmetic`, `algebra` |
| Categories, C-sets, bisimulation, automata, rewriting | morphe: `categories`, `equivalence`, `formal` |
| Cell complexes, homology, Gröbner bases | morphe: `topology`, `algebraic_geometry` |
| Units and physical constants | physis: `units` |
| Rotations, geodesics, map projections, differential geometry | physis: `geometry` |
| PDE discretisation, discrete exterior calculus, geometric integrators | physis: `discretisation` |
| n-body dynamics, potential theory, geophysical fluid dynamics, radiative transfer, reaction kinetics | physis |
| Astronomical time scales, celestial reference frames, Earth orientation | `core/metrology/` |
| Terrain analysis, spatial data models, land-cover classification, population and land use | `core/domains/geography/` |
| Seawater equation of state (TEOS-10) | `core/domains/oceanography/` |
| Moist thermodynamics, cloud microphysics, atmospheric chemistry | `core/domains/meteorology/` |
| Tide-generating potential, tidal prediction, atmospheric tides | `core/domains/tides/` |
| Bulk formulae for air–sea fluxes, wind waves | `core/domains/air_sea/` |
| Energy-balance models, orbital forcing of insolation | `core/domains/climate/` |
| Least squares, Kalman filters, 4D-Var, surrogate models | `core/inference/` |
| Downloading and registering DE440 or a tide-gauge record | `data/` |
| Telescope or sensor control | `observation/` |
| A 3D globe or planetarium renderer | `visualization/renderers/` |
| An orbital sandbox, a digital tide-predicting machine | `applications/<field or bridge>/`, e.g. `applications/astronomy/`, `applications/tides/` |
| Reference values from JPL Horizons or NOAA | `validation/` |
| A derivation about a field or bridge | `theory/<subject>/` |
| A derivation about mathematics, numerics or general physics | the library that implements it |
| A half-formed idea | `development/ideas/`, or the library's own if it concerns that library only |

## 8. Working Across Repositories

- **Submodules.** `lib/` holds the three libraries as git submodules, each pinned to an exact commit.
- **Editing a library.** Edit it in place under `lib/<name>`, commit there, then commit the updated pointer in the laboratory.
- **Commit order.** Commit in dependency order: calx and morphe first, then physis, then the laboratory.
- **Local rules.** Each library's CLAUDE.md governs work inside that library.
- **Links between repositories.** Use absolute GitHub URLs.
- **Provenance.** Every result records the commit of each library it used.
- **Build.** The laboratory's build adds the three libraries as subprojects, and physis uses the calx and morphe targets that already exist. Built on its own, physis fetches them. The build system is set up in Phase 0.
