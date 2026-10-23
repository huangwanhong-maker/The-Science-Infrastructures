# Architecture

*Status: draft, 2026-10-22.* This document describes how The Science Infrastructures is structured: its tiers, the order of its domains, the rules that connect them, and the computation and representation stacks. The decisions behind it are recorded in [development/decisions/](../development/decisions/); the present form is [ADR 0026](../development/decisions/0026-domains-as-libraries.md).

## 1. Principles

### 1.1 Knowledge by Subject, Work by Purpose

Two kinds of content are kept apart.

- **Knowledge is organised by subject.** It holds the objects of one science and their models: time scales, the ocean column, the motion of the plates. It is pure and reusable. It lives in the **domains**, on the **foundations** that hold what no single science owns.
- **Work is organised by purpose.** An atlas, a forecast experiment or a study is integrative: it draws on several domains and touches the world. It lives in the **applications** and in **research**.

*Rationale.* A unit that holds both must partition the sciences among its instances, and the sciences do not partition: geodesy serves astronomy and geophysics alike. Works may overlap freely in the domains they use; they *cover* the domains, and no domain is owned by a work.

### 1.2 Pure Core, Effectful Edges

The foundations and the domains perform pure computation: given inputs, they return results. They have no network access, no instrument I/O, no user interface and no hidden global state, and they read only files whose paths the caller passes in.

Everything that touches the outside world lives at the *edges* and calls the core: `data/`, `observation/`, `visualization/`, and the works.

*Rationale.* A pure core can be tested exhaustively and reproduced exactly. It can be deployed to any device, including devices without an operating system or a file system.

### 1.3 An Order of Domains

Each domain declares the domains it may depend on. The declared relation must be acyclic, so that its reflexive and transitive closure ≤ is a partial order on the domains.

- **A domain's library is its down-set** ↓d = { e : e ≤ d }. It builds and passes its checks with ↓d and its declared foundations alone.
- **Its height** is the length of the longest chain beneath it, the measure of how integrative it is.
- **A bridge** is a domain with several fields beneath it.
- **A work** that uses the domains S needs the order ideal ↓S = ⋃_{d∈S} ↓d.

### 1.4 Separation of Content from Realisation

Two parts of the infrastructure face a set of technologies that is open-ended and keeps changing:

- the hardware that executes computations;
- the media that present results.

In both, a stable, device-neutral intermediate representation separates *what* from *how*, and each technology is a backend that can be added without modifying anything else. For computation this is the stack in calx (§4). For presentation it is the representation layer (§5).

### 1.5 Own the Interfaces, Borrow the Implementations

The infrastructure defines its own small, stable interfaces. Mature external libraries may sit behind them and be replaced later. Every external dependency must have a licence compatible with MIT, or be adopted by an explicit decision.

### 1.6 Independent Foundations

General knowledge lives in foundation libraries, each in its own repository ([ADR 0008](../development/decisions/0008-independent-libraries.md), as amended by 0026):

- numerics and computation (calx);
- mathematical structures (morphe);
- geometry and physics (physis);
- inference under uncertainty (planned).

The test of placement: content that names an object in nature belongs to a domain; content that does not, and that more than one domain needs, belongs to a foundation.

### 1.7 Growth by Vertical Slices

Interfaces are extracted from working code that runs end to end. They are not designed completely in advance. Each slice passes through every tier at small scale ([roadmap](../development/roadmap.md)).

## 2. Tiers and Repositories

```
Works        applications/    one folder per application; each composes domains
             research/        questions · experiments · papers
             notebooks/       exploration and teaching
                    │ call into
Edges        data/            registry · acquisition · formats · products · requests on demand
             visualization/   scene model · renderers · styles   (used by works only)
             observation/     instruments · field work
                    │ call into
Domains      height 4   climate
(domains/)   height 3   tides · navigation · air-sea
             height 2   astronomy · geophysics · geography · oceanography · meteorology
             height 1   geodesy
             height 0   metrology · neuroscience (a submodule)
                    │ use
Foundations  physis           units · symmetry · geometry · dynamics · quantum · discretisation ·
(foundations/,                spacetime · field theory · gauge · quantum gravity · mechanics ·
 submodules)                  gravitation · thermodynamics · continuum · chemistry · radiation
             calx             compute stack · numerical methods
             morphe           arithmetic · categories · equivalence · algebra · topology ·
                              algebraic geometry · formal structures
             (inference)      estimation · filtering · assimilation · uncertainty · learning

Process      development/     decisions · ideas · roadmap · tools
Documents    docs/            architecture · gallery
```

The foundations depend on each other as follows; calx and morphe do not depend on each other.

```
physis ──► calx
   └─────► morphe
```

This repository holds the domains, the edges, the works and the documents. The foundations are repositories of their own, attached as git submodules, because they are general enough to be used outside the programme.

*Remark.* A domain is independent by the order and the checks of §3, not by a repository boundary. A domain or application that grows large, or gains collaborators of its own, may be extracted with its history into its own repository and attached as a submodule (ADR 0026, the graduation rule). Neuroscience is such a domain from its birth ([ADR 0028](../development/decisions/0028-neuroscience-and-the-brain-atlas.md)): `domains/neuroscience` is the submodule of `https://github.com/huangwanhong-maker/neuroscience`, with its own CLAUDE.md and logs, and the checks run it like any other domain.

## 3. Dependency Rules

1. **Foundations.** calx and morphe depend on no project library; physis depends on calx and morphe. No foundation depends on a domain, an edge, a work or this repository.
2. **Domains.** A domain depends only on the domains and foundations its `module.json` declares, and the declared order is acyclic.
3. **Through interfaces.** A unit reaches another only through its public interface: today `js/index.js`, later the C++ headers.
4. **The core stands alone.** Nothing in `foundations/` or `domains/` depends on an edge, a work, the process or the documents.
5. **Edges.**
   - `data/` and `observation/` depend only on the domains and foundations.
   - `visualization/` depends on those it declares.
6. **Works.**
   - Applications may depend on everything below them, and nothing depends on an application. Two applications may link to each other's views, but neither imports the other.
   - Only applications, `notebooks/` and `research/` may use `visualization/`.
7. **Tests.** A domain's tests import only what the domain may import, the harness of the checks and Node's built-in modules. A check that joins domains which do not stand on each other belongs to a work.

The rules for the JavaScript are checked by `node development/tools/check-dependencies.mjs`; `node development/tools/check-alone.mjs` runs each domain's checks with nothing but its down-set.

A violation of these rules is a design error. Fix it in one of two ways:

- move the code down;
- or define an interface in the lower tier that the upper tier implements (dependency inversion).

Never fix it with an upward import. When the right place for the code is a foundation, make the change in that library's repository.

*Example.* The Celestial Atlas needs the planetary ephemeris, but the domains never fetch files.

- `domains/astronomy` reads the JPL SPK format from bytes it is given.
- `data/` fetches DE440s, records its provenance and puts it in the cache.
- The atlas fetches the file from its server and passes the bytes to astronomy.

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

The full design is in calx's own [architecture document](https://github.com/huangwanhong-maker/calx/blob/main/docs/architecture.md); locally it is at `foundations/calx/docs/architecture.md`. The domains write bulk kernels through calx, and never address a device directly.

## 5. Representation Stack (`visualization/`)

The representation layer applies the same separation of content from realisation to output. Only applications, `notebooks/` and `research/` may use it.

| Component | Responsibility |
|---|---|
| `scene/` | A device-neutral, serialisable description of what is shown: bodies, fields, trajectories, reference frames, annotations, cameras, time. It carries units, time scales and frames from the domains, so renderers never guess them. |
| `renderers/` | One per technique. Examples: static figures (Matplotlib, SVG, PDF); web (WebGPU, three.js, CesiumJS); native (Vulkan); immersive (WebXR); physical (3D printing); non-visual (sonification). |
| `styles/` | Palette, typography and marks, shared by all renderers |
| `js/` | Provisional: the drawing helpers the atlases share (contours, particles on the sphere, display units, map tiles), until the scene model exists |

*Remark.* Astronomical scenes span roughly 10⁻³ m to 10²⁶ m, far beyond the range of single-precision graphics coordinates. The scene model therefore uses a hierarchy of local frames with floating origins, not one global coordinate system, as OpenSpace does.

## 6. Module Layout

### 6.1 A Domain

```
domains/<name>/
├── README.md        purpose, status, scope by module, dependencies, references, open questions
├── module.json      name · kind (field | bridge) · status · depends · foundations
├── js/              provisional JavaScript; js/index.js is the public interface
├── tests/           run-tests.mjs, and fixtures/ with each reference value's source and date
├── theory/          notes, in Markdown with LaTeX
└── (later)          the C++ modules, laid out as in §6.2
```

### 6.2 A C++ Module

A C++ module, in a foundation or a domain, has this layout:

```
<module>/
├── README.md                     purpose, status, dependencies, references, open questions
├── CMakeLists.txt                the target <library>::<module>, and its checks
├── include/<library>/<module>/   public headers: #include <calx/random/stream.hpp>, <neuroscience/snn/snn.hpp>
├── src/                          implementation, when the module is compiled
├── tests/                        checks: programs on calx's verification, registered with CTest
├── bindings/python/              Python bindings (optional)
└── benchmarks/                   performance measurements (optional)
```

- **Names** (provisional, ADR 0035): the namespace follows the headers, `calx::random`, `physis::dynamics`, and a domain's C++ is `<domain>::<module>` (`neuroscience::snn`).
- **Build:** CMake 3.20 or later and C++20; every library is a project, and the root `CMakeLists.txt` composes them in their order. A domain's `CMakeLists.txt` adds the foundations beneath it from `../../foundations` when no larger project has.
- **Where each area of mathematics goes,** and when to borrow a library rather than write one, is [docs/mathematics.md](mathematics.md).

Until a module has content, it exists only as a `planned` row in its parent README.

### 6.3 Provisional JavaScript in the Browser

The applications run in the browser, as ES modules without a build step. They import the domains, physis's stand-in and the visualization layer by relative paths, which resolve to the same files in Node and in the browser. Each application's local server serves its page at `/` and mounts the repository's `domains/`, `foundations/` and `visualization/` under the same names, beside the data cache at `/data/`.

## 7. Placement Guide

| Content | Location |
|---|---|
| A kernel that must run on several devices | IR through calx (`foundations/calx/compute/`) |
| A new device or accelerator | `foundations/calx/compute/backends/<target>/` |
| Linear algebra, FFT, ODE integrator | calx: `numerics/` |
| Exact arithmetic, groups, representations, lattices, presentations | morphe: `arithmetic`, `algebra` |
| Categories, C-sets, bisimulation, automata, rewriting | morphe: `categories`, `equivalence`, `formal` |
| Cell complexes, homology, fundamental groups, coloured graphs, Gröbner bases | morphe: `topology`, `algebraic_geometry` |
| Units and physical constants | physis: `units` |
| Lie groups, their representations, recoupling (3nj) symbols | physis: `symmetry` |
| Rotations, geodesics, curvature, map projections, differential geometry | physis: `geometry` |
| Exact solutions of general relativity, causal structure, Regge geometry, causal sets | physis: `spacetime` |
| Lattice fields, gauge fields on complexes, spin foams, group field theory | physis: `field_theory`, `gauge`, `quantum_gravity` |
| The spherical-harmonic transform, PDE discretisation, discrete exterior calculus, geometric integrators | physis: `discretisation` |
| n-body dynamics, potential theory, geophysical fluid dynamics, radiative transfer, reaction kinetics | physis |
| Least squares, Kalman filters, 4D-Var, surrogate models | the inference foundation (planned) |
| Astronomical time scales, celestial reference frames, Earth orientation | `domains/metrology/` |
| The figure of the Earth, its gravity field, reference systems, places and directions | `domains/geodesy/` |
| Ephemerides, apparent places, the rotation of the bodies, galaxies | `domains/astronomy/` |
| Plate motion, earthquake mechanisms, the crust, seismology | `domains/geophysics/` |
| Terrain analysis, spatial data models, land-cover classification, population and land use | `domains/geography/` |
| Seawater equation of state (TEOS-10), the ocean column, waves | `domains/oceanography/` |
| Atmospheric dynamics, moist thermodynamics, cyclone tracking | `domains/meteorology/` |
| Tide-generating potential, tidal prediction, atmospheric tides | `domains/tides/` |
| Bulk formulae for air–sea fluxes, wind waves | `domains/air-sea/` |
| Energy-balance models, orbital forcing of insolation, climate classification | `domains/climate/` |
| Reference spaces of the brain, parcellations and their maps, the head and the 10–20 system, readers of NIfTI and NRRD | `domains/neuroscience/` |
| A new science, such as chemistry or ecology | a new domain in `domains/`, placed in the order |
| Downloading and registering DE440 or a tide-gauge record | `data/` |
| Telescope or sensor control | `observation/` |
| A 3D globe or planetarium renderer | `visualization/renderers/` |
| A surface extracted from a volume, for any application to draw | `visualization/js/isosurface.js` |
| An orbital sandbox, a digital tide-predicting machine | `applications/<name>/` |
| Reference values from JPL Horizons or NOAA | the tests of the domain they verify, in `tests/fixtures/` |
| A derivation about a domain | `domains/<name>/theory/` |
| A derivation about mathematics, numerics or general physics | the foundation that implements it |
| A half-formed idea | `development/ideas/`, or a foundation's own if it concerns that library only |

## 8. Working with the Foundations

- **Submodules.** `foundations/` holds the libraries as git submodules, each pinned to an exact commit; `domains/neuroscience` is a submodule too.
- **Editing a foundation.** Edit it in place under `foundations/<name>`, commit there, then commit the updated pointer here. The same for a domain that is a submodule.
- **Commit order.** Commit in dependency order: calx and morphe first, then physis, then the domain submodules, then this repository.
- **Local rules.** Each foundation's CLAUDE.md governs work inside it.
- **Links between repositories.** Use absolute GitHub URLs.
- **Provenance.** Every result records the commit of each foundation it used.
- **Build.** The build adds the foundations as subprojects, and physis uses the calx and morphe targets that already exist. Built on its own, physis fetches them. The build system is set up in Phase 0.
