# CLAUDE.md

This file guides Claude Code in the repository of The Science Infrastructures ("the infrastructure"). It records the programme's purpose, design philosophy, architecture and working conventions.

This file can be revised: when a decision recorded here changes, update this file in the same change.

- The detailed design is in [docs/architecture.md](docs/architecture.md).
- Decisions and the reasons for them are in [development/decisions/](development/decisions/). The present structure is [ADR 0026](development/decisions/0026-domains-as-libraries.md).
- Each foundation library under `foundations/` has its own CLAUDE.md, which governs work inside it.

## 1. Purpose and Scope

This is an integrative programme to learn, use and research all domains of science, in theory and in practice, as one connected body of knowledge.

- **Three tiers** ([ADR 0026](development/decisions/0026-domains-as-libraries.md)):
  - **Foundations:** calx (computation and numerical methods), morphe (mathematical structures and their equivalences), physis (units, geometry, analysis, dynamics and theoretical physics), and two planned libraries: one for inference, one for the methods of the quantitative social sciences (measurement, surveys, causal inference, networks, games and social choice, agent-based models). They hold laws and methods that name no object in nature. Their mathematics is mapped, area by area against the whole Mathematics Subject Classification, in [docs/mathematics.md](docs/mathematics.md) ([ADR 0035](development/decisions/0035-the-mathematical-foundations.md)); a module is written when a domain needs it.
  - **Domains:** one library per science, holding its objects and their models. A domain is a field or a bridge between fields.
  - **Works:** applications, research and notebooks, organised by purpose, each composing the domains it needs.
- **Shared edges** touch the world for every work: data, visualization, observation, and the API through which other programs reach the infrastructure (ADR 0030).
- **The first domains** are those of the Earth and the sky:
  - **Fields:** metrology, geodesy, astronomy, geophysics, geography, oceanography and meteorology.
  - **Bridges:** tides, navigation, air–sea exchange and climate.
- **Scope of the first fields:**
  - **Geography** covers all its branches: physical geography, cartography and GIS, remote sensing, human geography, and others as the work requires.
  - **Meteorology** is atmospheric science in the broad sense; ocean and atmosphere are one coupled system ([ADR 0007](development/decisions/0007-atmosphere-and-ocean-as-coupled-system.md)).
  - **Geodesy** supports the other Earth sciences.
  - **Geophysics** is the solid Earth: plate motion, earthquakes and the crust ([ADR 0024](development/decisions/0024-the-solid-earth.md)).
- **Neuroscience** is the first domain beyond the Earth and the sky ([ADR 0028](development/decisions/0028-neuroscience-and-the-brain-atlas.md)): the brain's reference spaces, parcellations and maps, and the head with its electrodes, drawn by the Brain Atlas. It is a repository of its own, attached as a submodule.
- **History** is the first domain of the humanities ([ADR 0041](development/decisions/0041-the-history-domain.md)), built at the user's request for the History Complex. It covers records of histories coded by a codebook, their chronology across calendars, the 2-complex a record makes, and the spin foam of its signed ties, computed exactly on physis; and the history read as a field theory, one Feynman graph of a group field theory and, in time, one sample path ([ADR 0042](development/decisions/0042-the-history-as-a-field-theory.md)).
- **Genomics and evolution** are the first domains of the life sciences ([ADR 0043](development/decisions/0043-the-genome-infrastructure.md)), built at the user's request for genome analysis and drawn by the Genome Atlas. Genomics (height 0) holds genomes and their sequences: the alignment of sequences as a sum over an edit graph in a semiring with its fast schedules, indexes and sketches, k-mers and assembly, mapping and variants, HMMs and annotation, multiple alignment, expression, and simulation with its truth; evolution (height 1, on genomics) holds phylogeny, population genetics and selection. Their general mathematics is in the foundations: words, their indexes and their alignment in morphe's `formal`; vector lanes, task graphs and the statistics they need in calx.
- **Economics** is the first domain of the social sciences beside history ([ADR 0044](development/decisions/0044-value-in-space-and-time.md)), built at the user's request for the circulation of value in space and in time and drawn by the Economy Atlas. It holds value carried between positions: transports and their loop deficits (prices exactly when every deficit vanishes), the circulation of flows on a network whose faces geography supplies, production as a circular flow, and the spacetime complex of a record of transactions. Its general mathematics is in the foundations: discrete exterior calculus and points on the sphere in physis, conjugate gradients in calx, the worldline complex in morphe.
- **Further domains** join the same order as the work reaches them. Every science is meant to join, the social sciences and the humanities included; the README's §4 lists them all once, as foundations, domains or modules, with a priority for each. The foundations come next: morphe, the algebraic foundation (category theory, algebra, topology, geometry), and physis, the foundation of theoretical physics, whose first slices of complexes, symmetry and fields are written ([ADR 0039](development/decisions/0039-complexes-symmetry-and-fields.md)). Then the next five domains, in order: spiking neural networks (a module of neuroscience), computational linguistics, genomics, molecular science, and drug discovery (a bridge over genomics and molecular science).
- **No laboratories.** The word is retired as the name of what we are. A themed grouping of works, if ever needed, is a *programme*: a README and a roadmap, owning no code. Do not use "laboratory" for any part of the infrastructure; proper names (the Nevada Geodetic Laboratory) and the decision records keep it.

Every artefact serves one of three activities, and says which:

- **Learning:** derivations, worked examples and expositions, written to be understood.
- **Use:** tested tools whose numbers can be relied on, such as ephemerides, frame transformations and tidal predictions.
- **Research:** open questions, experiments and new results, built on the first two.

*Rationale.* The sciences share their structures: symmetry and its representations, conservation laws and the equations that express them, networks, sequences and grammars, and inference under uncertainty. The first domains share more still: spherical geometry and coordinate systems, the rotation group, potential theory and harmonic analysis, and the dynamics of rotating, stratified fluids. They meet in tides, navigation, the figure and mapping of the Earth, the exchange between ocean and atmosphere, and the astronomical forcing of climate. One infrastructure makes these connections explicit, and lets one tested implementation of each method serve every science.

## 2. Design Philosophy

The infrastructure will grow over a long time and will change direction as new ideas emerge. Its governing principle is **modular, extensible, revisable**.

### 2.1 Modularity

- **Knowledge by subject, work by purpose.** Domains hold knowledge; works compose it. Works cover the domains; no work owns a domain.
- **An order of domains.** Each domain declares in `module.json` the domains and foundations it may use. The declared relation is acyclic, so that its closure is a partial order, and a domain's library is its down-set.
- **Public interfaces.** Units use each other only through their public interfaces: today `js/index.js`, later the C++ headers.
- **Reuse before writing.** A domain uses what the foundations already hold and never writes it again; a general method a domain needs (one that names no object in nature) is written in the foundation where it belongs, and the domain then uses it. A reusable, integrative foundation is the purpose of the repository (the user's principle, [ADR 0042](development/decisions/0042-the-history-as-a-field-theory.md)).
- **Independence is checked.** The rules are verified by `development/tools/check-dependencies.mjs`, and each domain's checks run alone with its down-set by `development/tools/check-alone.mjs`.

### 2.2 Extensibility

- A new science is a new domain; a new method, device or presentation technique is a new module. Adding one does not mean editing existing ones.
- Computing devices and presentation media keep changing. In both places, a stable, device-neutral representation separates *what* from *how*, and each technology is a pluggable backend (§3.5, §3.6).
- Own the interfaces and borrow the implementations: mature external libraries may sit behind the infrastructure's interfaces.
- Prefer small, composable units. Introduce an abstraction only after its pattern has appeared twice.

### 2.3 Revisability

- Every module, note and model carries a status:
  - `planned`: named in the architecture, with no content yet;
  - `seed`: an idea, possibly wrong;
  - `draft`: being worked out;
  - `stable`: tested and relied upon;
  - `deprecated`: superseded, and kept with a pointer to its replacement.
- Significant decisions go into the decision log (§6.3).
- Prefer reversible choices, and mark provisional ones where they are made.
- Capture ideas early (§6.2) and promote them once they mature.

### 2.4 Unity of Theory and Computation

- Every implemented formula cites its derivation or source, with equation numbers.
- Every theoretical result that can be checked numerically should be checked by code.
- Notation in the notes corresponds to names in the code, and the glossary maps one to the other.
- A domain's notes live beside its code, in its `theory/`.

### 2.5 Verifiability

- Numerical code is validated against authoritative references, such as IAU SOFA/ERFA, JPL Horizons, the IERS Conventions and NOAA tide predictions.
- Reference values live in the tests of the domain they verify (`tests/fixtures/`), each with its source and retrieval date.
- Docstrings state limits of accuracy: model error, range of validity, and approximations made.

## 3. Architecture

This section summarises [docs/architecture.md](docs/architecture.md). Read that document before adding a domain or a module.

### 3.1 Repositories

| Repository | Path | Scope | Depends on |
|---|---|---|---|
| The Science Infrastructures | `.` | Domains, shared edges, works, documents, process | the foundations |
| calx | `foundations/calx` | Retargetable computation stack; numerical methods on floating-point numbers and arrays. C++ written: random, linalg (with the symmetric eigenproblem and determinants), integrators, special functions (with the gamma family), quadrature, solvers, optimisation (Brent's minimiser), differentiation, verification, the scheduler's thread pool, task graphs and dynamic loops, the CPU backend's target and vector lanes, conjugate gradients on sparse rows | — |
| morphe | `foundations/morphe` | Mathematical structures and their equivalences: logic, arithmetic, algebra, combinatorics, categories, topology, algebraic geometry, formal structures. C++ written: arithmetic (integers without bound, rationals), algebra (the Smith normal form, homology with generators, presented and finite groups), combinatorics (graphs, canonical forms, Eulerian circuits and the BEST theorem), topology (CW and simplicial complexes, π₁, surfaces, collapses, dual 2-complexes, coloured graphs, persistent homology and cycles over ℤ/2, Q-analysis, worldline complexes), formal (words, suffix arrays, the BWT, FM-indexes, Aho–Corasick, the alignment of words over semirings, the polytope algebra) | — |
| physis | `foundations/physis` | Units, geometry, analysis, dynamics, control, discretisation, and physical law: mechanics, gravitation, thermodynamics, continuum, chemistry, radiation. C++ written: units, dynamics, quantum (the open two-level system), symmetry (Lie groups, recoupling), geometry (metrics, curvature), spacetime (exact solutions, Regge geometry, causal sets), field theory (lattice fields, invariant theories of every rank and their large-N solutions), gauge (lattice gauge theory on 2-complexes, the exact contraction of ℤ₂ theories), quantum gravity (Turaev–Viro, group field theories, spin foams, random geometries, condensates, tensor and matrix models), discretisation (cochains, the Hodge decomposition, flows with conductances; a provisional JavaScript stand-in for spectral transforms), computational geometry (points on the sphere, their spherical Delaunay triangulation) | calx, morphe |
| neuroscience | `domains/neuroscience` | A domain in its own repository (ADR 0028): the brain's reference spaces, parcellations, maps and the head, its signals, and networks of spiking neurons in C++ (`snn/`) | calx, physis (declared; no domain) |

The foundations and the neuroscience domain are git submodules. A domain or application that grows large, or gains collaborators of its own, may be extracted with its history into its own repository (ADR 0026, the graduation rule); neuroscience was created so at the user's request.

### 3.2 The Order of the Domains

```
height 4   climate
height 3   tides · navigation · air-sea
height 2   astronomy · geophysics · geography · oceanography · meteorology
height 1   geodesy · evolution
height 0   metrology · neuroscience · history · genomics · economics
```

A bridge is a domain with several fields beneath it. A new domain is placed in the order by what it may depend on.

### 3.3 Pure Core, Effectful Edges

The foundations and the domains perform pure computation:

- no network access;
- no instrument I/O;
- no user interface;
- no hidden global state;
- only files whose paths are passed in are read.

The edges (`data/`, `observation/`, `visualization/`, `api/`) and the works (`applications/`, `research/`, `notebooks/`) deal with the outside world and call the core.

### 3.4 Dependency Rules

- calx and morphe are independent of each other; physis depends on both. No foundation depends on a domain, an edge or a work.
- A domain depends only on what its `module.json` declares, through public interfaces.
- Nothing in `foundations/` or `domains/` depends on an edge or a work.
- `applications/` may depend on everything below it, and nothing depends on an application.
- Only `applications/`, `notebooks/` and `research/` may use `visualization/`.
- A domain's tests import only what the domain may import, the harness (`development/tools/checks.mjs`) and Node's built-in modules. A check that joins domains which do not stand on each other belongs to a work.
- An upward dependency is a design error. Fix it by moving code down, or by dependency inversion.

### 3.5 Computation Stack (calx)

- calx separates representation, schedule, compiler, runtime, task scheduler and backends, following TVM, Halide and MLIR.
- A new device is a new backend directory, with no other change.
- The CPU reference interpreter defines correct results.
- Precision is declared for each computation and never silently reduced.
- The domains never address a device directly.

### 3.6 Representation Layer

`visualization/` separates three things:

- a device-neutral scene model (`scene/`);
- one renderer per technique (`renderers/`);
- shared styles (`styles/`).

Units, time scales and frames travel with the data into the scene. Until the scene model exists, drawing helpers that any application could use live in `visualization/js/`.

### 3.7 Placement

- Content that names an object in nature belongs to a domain. General content that more than one domain needs belongs to a foundation; make the change in that library's repository.
- Content that serves one purpose, and touches the world, belongs to a work or an edge.
- Use [architecture §7](docs/architecture.md) for the details.

## 4. Scientific and Mathematical Conventions

### 4.1 Units, Time Scales and Reference Frames

These are the most common source of silent error in this subject, so they are always explicit.

- **Units.** Use SI units internally. Use radians for angles internally, and convert only at input and output boundaries. When a quantity's unit is not obvious, carry it in the name or the type.
- **Time.** A time value always states its scale (UTC, TAI, TT, TDB or UT1) and its representation (JD, MJD, two-part JD or ISO 8601). Epochs are named, for example J2000.0 = JD 2451545.0 TT.
- **Reference frames.**
  - Every coordinate states its frame, and where relevant its epoch and equinox.
  - Examples of frames: ICRS, GCRS, CIRS, ITRS, ecliptic of date, topocentric, local horizon.
  - Transformations between frames are explicit functions and never happen implicitly.
- **Figure of the Earth.** State the ellipsoid or geoid in use, such as WGS 84 or GRS 80.
- **Orientation conventions.** Define each of these once in the glossary and follow it everywhere:
  - handedness;
  - the origin of azimuth;
  - the sign of longitude;
  - whether rotations are active or passive.

Domains of other sciences add their own conventions in their READMEs, such as the nomenclature of chemistry.

### 4.2 Notation, Rigour and Glossary

- Write notes in LaTeX, with the precision of a mathematics text:
  - state the hypotheses;
  - distinguish definitions from results;
  - give proofs or cite them.
- Where it is natural, show algebraic structure rather than hiding it: groups, representations, lattices, characters and complexes.
- `docs/glossary.md` lists every symbol with its meaning, units, convention and code name. It is created in Phase 0.

### 4.3 References

Cite primary sources and standard references, and collect them in `docs/references.bib`. Examples:

- the *Explanatory Supplement to the Astronomical Almanac*;
- the IERS Conventions and IAU resolutions;
- Meeus, *Astronomical Algorithms*;
- Pugh and Woodworth, *Sea-Level Science*;
- Gill, *Atmosphere–Ocean Dynamics*;
- Vallis, *Atmospheric and Oceanic Fluid Dynamics*;
- Holton and Hakim, *An Introduction to Dynamic Meteorology*;
- Aki and Richards, *Quantitative Seismology*.

### 4.4 Writing Style

These conventions are shared with the foundation libraries.

- Headings are declarative, numbered and academic, for example "3.2 Precession of the Equinoxes". Never use question forms.
- Present concepts and philosophy before structure and usage.
- Label asides with *Rationale.*, *Remark.* or *Caveat.*
- Code identifiers, docstrings and comments are in English.

## 5. Repository Layout

```
.
├── foundations/       calx · morphe · physis (git submodules); inference planned
├── domains/           one library per science (ADR 0026)
│   ├── metrology/     time scales · reference frames · Earth orientation
│   ├── geodesy/       figure of the Earth · places · directions · reference systems
│   ├── astronomy/     ephemerides · apparent places · rotation of the bodies · galaxies
│   ├── geophysics/    plates · earthquake mechanisms · crust · sections
│   ├── geography/     cartography · GIS · physical and human geography · remote sensing
│   ├── oceanography/  the ocean column · waves · seawater · sea floor
│   ├── meteorology/   dynamics · our models · cyclones · wind
│   ├── tides/ · navigation/ · air-sea/     bridges (planned)
│   ├── climate/       bridge: climate classes · energy balance · orbital forcing
│   ├── neuroscience/  reference spaces · parcellations · maps · the head (git submodule, ADR 0028)
│   ├── history/       the record and its chronology · the history complex · its readings and spin foam · its field theory · canons (ADR 0041, 0042)
│   ├── genomics/      sequences · simulation · alignment · indexes · sketches · k-mers · assembly · mapping · variants · HMMs · annotation · multiple alignment · expression (ADR 0043)
│   ├── evolution/     phylogeny · population genetics · selection (ADR 0043)
│   └── economics/     transfer networks · circulation · production · the spacetime of transactions · places (ADR 0044)
├── data/              dataset registry · acquisition · formats · products · provenance
├── visualization/     scene model · renderers · styles (ADR 0005); provisional js/
├── observation/       instruments · field work · sensors
├── api/               the service and its clients: Python · C and C++ · JavaScript (ADR 0030)
├── applications/      one folder per application: celestial-atlas · earth-atlas · brain-atlas · spiking-simulator · complex-atlas · history-complex · genome-atlas · economy-atlas
├── research/          questions · experiments · papers
├── notebooks/         exploration and teaching
├── docs/              architecture · mathematics · gallery (later: glossary · references)
├── development/       decisions · ideas · release-notes · roadmap (ADR 0009) · tools
├── assets/            images · fonts · style tokens
└── CMakeLists.txt     the C++ of the foundations and the domains, in their order (ADR 0035)
```

- Every directory has a README giving its purpose, its planned modules with their status, and its dependencies.
- A subdirectory is created when its first content is written. Until then it is a `planned` row in its parent README.
- The exceptions are:
  - the domain directories of `domains/`, each with its README and `module.json` from the start;
  - the component directories of calx's `compute/`.
- Paths are lowercase and hyphenated, with no spaces (`air-sea`).
- When the layout changes, update this section, the README's §6 and the affected parent READMEs in the same change.

## 6. Working Procedures

### 6.1 The Contract of a Domain and a Module

A domain directory contains:

- `README.md`: purpose, status (§2.3), scope by module, dependencies, references and open questions;
- `module.json`: `name`, `kind` (`field` or `bridge`), `status`, `depends` and `foundations`;
- `js/`: provisional JavaScript, behind `js/index.js`, until the C++ exists;
- `tests/`: its checks (`run-tests.mjs`, on the shared harness), with reference values in `tests/fixtures/`;
- `theory/`: its notes;
- later, its C++ modules, which follow [architecture §6](docs/architecture.md).

### 6.2 Idea Capture

- Write each new idea in `development/ideas/` as `YYYY-MM-DD-slug.md`, with status `seed`, and add it to the index.
- An idea that concerns only one foundation goes in that library's own `development/ideas/`.
- Promote an idea by moving its substance and leaving a link behind.
- Never delete an abandoned idea. Mark it `abandoned` and give one sentence explaining why.

### 6.3 Decision Log

- Record decisions that shape structure, conventions or tooling in `development/decisions/NNNN-slug.md`, using the template in the index.
- Never rewrite a record. A new record supersedes or amends it, and the old record's status is updated to point to the new one. Paths in old records are as they were written; ADR 0026 maps them to the present.
- Records accepted so far: 0001–0044. Records 0001–0024 were written in the first laboratory, which the infrastructure absorbed; 0025 founded this repository.

### 6.4 Release Notes

- At each commit, add to `development/release-notes/YYYY-MM-DD.md` of each repository committed (one file per day, one short section per commit, newest first) a clear, concise summary in plain prose of what is new, what one can now do and what works better, and index the day in that folder's README ([ADR 0037](development/decisions/0037-release-notes.md)).
- Not a list of changed files or functions: the git log keeps those.
- Not between commits.

### 6.5 Change Discipline

- Before restructuring directories, adding or renaming a domain, moving a foundation, or changing a convention in this file, propose the change and wait for the user's agreement.
- Keep each change to one concept.
- When adding a domain or module, update:
  - its README and `module.json`;
  - its parent README's table;
  - the layout (§5), if it adds a directory here;
  - the glossary, if it introduces notation.
- If this file no longer matches the repository, say so and propose a fix.

### 6.6 Checks

- `node development/tools/run-tests.mjs` checks the dependency rules, runs every suite (the foundations' C++, built with CMake here or in WSL; the foundations' JavaScript; the domains in their order, each with its C++; visualization, the data library's Python checks in `data/tests/`, the API, applications), then runs each domain alone with its down-set. Run it before every commit that touches code.
- C++ checks are programs on calx's `verification`, registered with CTest; `development/tools/checks.mjs` builds a CMake project and counts them (`cppChecks`), keeping its builds in `~/.cache/science-build/`. One library's C++ alone: `node development/tools/run-cpp.mjs foundations/calx foundations/calx`.
- A domain's own suite: `node domains/<name>/tests/run-tests.mjs`.
- Check the applications visually by driving headless Chrome over the DevTools protocol, with an absolute `--user-data-dir`.

### 6.7 Working with the Foundations

- **Submodule commits.** Edit a foundation in place under `foundations/<name>` and commit there first. Then commit the updated submodule pointer here. Commit only the pointer, by path, when other changes are staged. The same holds for the neuroscience domain under `domains/neuroscience`.
- **Commit order.** calx and morphe first, then physis, then neuroscience, then this repository.
- **Each foundation's own CLAUDE.md** governs work inside it.
- **Links between repositories.** Use absolute GitHub URLs.
- **Shared conventions.** The writing style, status vocabulary, and decision, idea and release-note conventions are shared by every repository. Change them in all of them together.
- **Pushing.** Push only with the user's agreement.

### 6.8 Collaboration Style

- The principal investigator is fluent in abstract algebra and formal language theory. Prefer precise mathematical formulations.
- This is also a learning project. When introducing physics or conventions specific to a domain, give the reasoning behind them and cite their source.
- Visualizations should look warm and humane, not like a generic dashboard.

## 7. Technology

- **Decided** ([ADR 0003](development/decisions/0003-cpp-computational-core.md)): C++ for computational code, in the foundations and the domains, wherever possible. Python or another suitable language for data acquisition, orchestration, notebooks, machine learning and visualization.
- **Provisional** (in use since [ADR 0035](development/decisions/0035-the-mathematical-foundations.md), to be confirmed in roadmap Phase 0): C++20 without extensions; CMake 3.20 or later, each library a project and each module a target `<library>::<module>`, composed by the root `CMakeLists.txt`; GCC 11 or later (in WSL on the development machine, GCC 13 on the server); checks on calx's `verification` rather than a borrowed test framework; nanobind for Python bindings (not yet); pytest; Python 3.12 or later.
- **Own and borrowed implementations** ([docs/mathematics.md §10.3](docs/mathematics.md)): the libraries own their interfaces; mature, licence-compatible libraries may sit behind them (Eigen, OpenBLAS, pocketfft, SUNDIALS), and Python work uses NumPy and SciPy directly; our own C++ is written where it serves several devices, exploits structure, gives control (reproducible streams, declared precision) or is how a method is learnt.
- **Provisional JavaScript.** Until the C++ exists, the domains' code is ES modules in their `js/`, as is physis's stand-in for discretisation (physis ADR 0003); the C++ must pass the same checks.
- **Web applications:**
  - Use ES modules without a build step.
  - Load three.js from jsDelivr, pinned to an exact version (0.170.0).
  - Serve locally with Python's standard library. Each application's server mounts `domains/`, `foundations/` and `visualization/` beside its page, and the data cache at `/data/`; it reports the application's data recipe (`recipe.json`) at `/setup` and runs it on the page's request.
  - Test with Node (v22 is installed).
  - Each application has its own local server: the Celestial Atlas on port 8770, the Earth Atlas on 8790, the Brain Atlas on 8810, the Spiking Simulator on 8830, the Complex Atlas on 8850, the History Complex on 8870, the Genome Atlas on 8890, the Economy Atlas on 8910.
  - The Spiking Simulator's server runs the neuroscience domain's C++ engine, built by its recipe with g++ (natively, or in WSL) into `data/cache/engines/snn/` ([ADR 0038](development/decisions/0038-the-spiking-simulator.md)). The Complex Atlas's server runs an engine of its own on morphe and physis, built the same way into `data/cache/engines/complex/` ([ADR 0039](development/decisions/0039-complexes-symmetry-and-fields.md)). The History Complex's runs one on the history domain, morphe and physis, built into `data/cache/engines/history/` ([ADR 0041](development/decisions/0041-the-history-domain.md), [ADR 0042](development/decisions/0042-the-history-as-a-field-theory.md)). The Genome Atlas's engine, on the genomics and evolution domains, is built by its recipe into `data/cache/engines/genome/`, a unit for each group of operations compiled one at a time (one unit holding both domains needed more memory than WSL gave without failing), and kept alive as workers: local ones, and others started by a command the user writes (ssh), among which jobs are cut into shards ([ADR 0043](development/decisions/0043-the-genome-infrastructure.md)). The Economy Atlas's engine, on the economics domain, physis and the API's JSON, is built the same way into `data/cache/engines/economy/`; its staleness is judged by the headers it includes, followed transitively ([ADR 0044](development/decisions/0044-value-in-space-and-time.md)).
  - A public instance runs on the user's server, `eu_region_server` (173.212.213.207): the Celestial Atlas on port 9200, the Earth Atlas on 9201, in the tmux session `science`, served with `--host 0.0.0.0 --sibling`; the Brain Atlas on 9202, in the window `brain-atlas`, served with `--host 0.0.0.0` after its recipe. The READMEs link to it; keep the links current if it moves.
- **Visual style:**
  - The Celestial Atlas uses Dusty Orchid (#CA8DBE) as its main colour, chosen by the user: plum nights, orchid accents and orchid-tinted paper for charts.
  - The Earth Atlas shares that palette for now (ADR 0014), with data in perceptually uniform colour maps (cmocean, and our own orchid map); the user may give it a colour of its own.
  - The Complex Atlas uses Old Rose (#C08081), provisionally, until the user chooses its colour.
  - The History Complex uses a faded mauve (#9C6B85) on parchment, provisionally; its foam is drawn on a plum night.
  - The Economy Atlas uses a warm ochre (#A8743A) on ivory, provisionally: value drawn in terracotta, value sent out slate.
  - The Genome Atlas uses mulberry (#A8577E) on blush paper, provisionally, with the four bases in softened conventional colours (A sage, C dusty blue, G ochre, T rose).
  - All interfaces use serif type (Source Serif 4).
  - Each application may carry its own main colour within the warm, gentle family the user prefers.
- **Version control.**
  - This repository's remote is `https://github.com/huangwanhong-maker/The-Science-Infrastructures`; the foundations' are `https://github.com/huangwanhong-maker/{calx,morphe,physis}`; the neuroscience domain's is `https://github.com/huangwanhong-maker/neuroscience`.
  - Each repository fetches over HTTPS and pushes over the SSH alias `serendip_wanhong` (its `pushurl`).
  - The first laboratory's repository, `https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory`, is superseded by this one (ADR 0026); archiving it is the user's decision.
- **Licences.** Every repository is MIT-licensed.
  - Check every external library for licence compatibility before adopting it. Copyleft examples: FFTW, SHTns, FLINT, GMP, and parts of CGAL.
  - Datasets keep their own licences, and are fetched into `data/cache/`, never committed.
- **Model environment.**
  - AI models run in `.venv` at the repository root (Python 3.12, not committed), with the packages pinned in `data/products/requirements-models.txt` (ADR 0022).
  - The rest of the data library needs the standard library, numpy and Pillow (`data/requirements.txt`).
  - Model weights under non-commercial licences (Pangu-Weather: CC BY-NC-SA 4.0) are fetched into the cache and never committed.
- **Development machine** (as of 2026-10-22):
  - Windows 11; Intel i7-11370H; 32 GB of memory.
  - NVIDIA RTX 3050 Ti Laptop GPU (4 GB, CUDA 12.9) and Intel Iris Xe.
  - Python 3.12, Node 22, Julia and WSL are installed. There is no C++ compiler or CMake on the Windows side; WSL has g++ 11.5, make, Python 3.12 with numpy and Node 18.19, and its networking is mirrored. The C and C++ checks (the API's clients, the C++ plug-ins) build with WSL's g++ when no native compiler is found.
  - The repository's path, `The Science`, contains a space, so quote paths.

## 8. Open Decisions

1. **The C++ namespaces of the domains**, and of the infrastructure as a whole. *Provisionally* `<library>::<module>` and `<domain>::<module>` (ADR 0035).
2. **C++ toolchain** on Windows (MSVC, clang/LLVM or WSL), together with the build and packaging tools. *Provisionally* GCC in WSL with CMake (ADR 0035).
3. **The inference foundation:** its name and its first module; and the foundation for the quantitative social sciences: its name, and its boundary with the inference foundation.
4. **Automatic differentiation** in calx, which would be built as an IR transformation.
5. **Precision policy details** (calx).
6. **Notes format.** Markdown with LaTeX (the current default), Jupyter or Quarto; LaTeX for papers.
7. **Languages of the documents.** English, with abstracts in Chinese and Japanese as in the README, or more.
8. **Licence for prose.** MIT covers software; notes may suit CC BY 4.0.
9. **First renderer technology** for the scene model.
10. **The next domains:** the order is set (README §4.6): first the foundations morphe and physis, then spiking neural networks, computational linguistics, genomics, molecular science, drug discovery (tier 1, the user's); genomics has begun, with evolution beside it (ADR 0043); tiers 2 and 3 are proposed; each domain is placed in the order, with its repository (here or its own), when its work begins. Which modules of theoretical physics physis adds beyond those it plans is physis's decision to record; its decision 0005 added symmetry, spacetime, field theory, gauge, quantum gravity and electromagnetism (ADR 0039).
11. **The first laboratory's repository on GitHub:** archive it with a pointer here, or leave it.

## 9. Current State

As of 2026-10-24:

- **The structure of ADR 0026 is in place.** The first laboratory's history is joined to this repository, and the foundations are submodules under `foundations/`.
- **The domains.** Thirteen domains have code: metrology, geodesy, astronomy, geophysics, geography, oceanography, meteorology, climate, neuroscience, history, genomics, evolution and economics. Tides, navigation and air–sea are planned. The domains of the Earth and the sky are provisional JavaScript, with 422 checks against published references and closed-form answers; neuroscience's networks of spiking neurons and the whole of the history domain are C++, and their suites run 274 and 375 checks; economics is C++ too, with 68. With the foundations (their C++ among them), the edges and the applications, the whole run is 4689 checks in twenty-eight suites. Every domain passes its checks alone.
- **The Brain Atlas** (`applications/brain-atlas`, [ADR 0028](development/decisions/0028-neuroscience-and-the-brain-atlas.md), [ADR 0032](development/decisions/0032-the-anatomy-of-the-brain-atlas.md)) shows the human brain in the round from an anatomy product (`data/products/brain_anatomy.py`): the cortex by gyrus, the deep nuclei, ventricles and fibre bundles, every part with a card from the domain's knowledge base (152 structures, with references). It cuts the brain and shows sections, colours the cortex by Julich-Brain or Schaefer, maps data by structure and by coordinate, places the 10–20 electrodes on the scalp and carries their values onto the cortex, scales the brain to a subject's head, and shows the Allen mouse brain; five guided tours. Its Signals mode ([ADR 0034](development/decisions/0034-biosignals-on-the-brain.md)) carries EEG onto the cortex through a three-shell forward model and sLORETA (band power, the potential at an instant, contrasts of conditions, with the traces playing beneath), shows one current's EEG potential and MEG field side by side, and reads fMRI maps through the cortex, smoothed and with their peaks named.
- **The Spiking Simulator** (`applications/spiking-simulator`, [ADR 0038](development/decisions/0038-the-spiking-simulator.md)) composes networks of spiking neurons from the engine's catalogue of twenty-two models, classical and quantum, edits them as the engine's own JSON specification, and runs them in the C++ engine on the server: eight examples (an adapting neuron, Brunel's four states, the classes of excitability in one run, quantum emitters, STDP). Its visualizer and profiler draw a run's raster and rates, its statistics (rates, CV and LV, Fano factors, synchrony χ, spectra, f–I curves), states, weights and network, and the engine's profile phase by phase and thread by thread; the statistics and the profile are the engine's own (`snn/analysis.hpp`).
- **The Complex Atlas** (`applications/complex-atlas`, [ADR 0039](development/decisions/0039-complexes-symmetry-and-fields.md)) shows cell complexes, their invariants and the fields on them: twenty-five families (surfaces glued by words, projective and lens spaces, presentation complexes, lattices, triangulations, dual 2-complexes, coloured graphs) or the user's own words, facets or JSON; three views (a realisation in three dimensions, the net of polygons, the diagram of cells); homology with generators drawn on the complex, π₁ raw and simplified, the classification of surfaces, collapses, manifolds and Gurau degrees; lattice gauge fields of ℤ_N, U(1), SU(2) and SU(N) by Monte Carlo against Migdal–Witten's exact values; and Turaev–Viro state sums. Its engine is C++ on morphe and physis, run by its server.
- **The History Complex** (`applications/history-complex`, [ADR 0041](development/decisions/0041-the-history-domain.md)) builds histories into 2-complexes from the history domain's canons: occasions are vertices, the lives of beings edges, and each tie between two of its meetings a lens, a 2-cell. Seven views: the foam in three dimensions (lives rising through time, cut at any date, opened into the refined complex K̂, coloured by discipline, spin, flips or inference, with guided tours of each canon), the loom, the moment (the spin network at a date, its exact reconstruction, and the spin networks along which the foams cross it), the complex (loops, shells, the other complexes, persistent homology in time), the field (the history as one Feynman graph of a group field theory, its bubbles, action and exact amplitude, and in time its rates and martingale checks), the record, and the readings (persistence and balance fitted by their exact likelihood, the labels of the spin foam, Wilson loops, power, SU(2) kinematics). Six canons, each with a guided tour: Athens 470–322 BC, Florence and Rome 1401–1564, the Scientific Revolution 1543–1727, Europe 1785–1828, the European alliances 1871–1914 and Vienna 1890–1938, and all of them joined. Its engine is C++ on the history domain, morphe, physis and calx, run by its server.
- **The history domain** ([ADR 0041](development/decisions/0041-the-history-domain.md)), the first of the humanities, is C++ from the start: records coded by a codebook (v1), chronology across the Gregorian, Julian and Old Style calendars (a reader's dates in the historians' convention, Julian before 1582), the join of records (an id identifies across canons only where they agree who it is, by Wikidata item or name), the history complex and its readings (persistence, flips, frustrated triangles, closed shells, moments), and the spin foam of the signed ties at ℤ₂ (the theory note, §6): the Gibbs weight on K̂ as physis's lattice gauge theory, its character expansion over the ℤ/2 2-cycles, spins, the Gauss law, Wilson loops, the posterior and power, all computed exactly by physis's contraction along the order of time on every single canon (the joined canons by physis's heat bath), with the couplings the maximum a posteriori of the exact likelihood ([ADR 0042](development/decisions/0042-the-history-as-a-field-theory.md)); the history as a field theory (§8): one Feynman graph of a group field theory over ℤ₂, whose amplitude is 2^(F + E) Z exactly, and in time one sample path with its rates and martingale checks; and other complexes of a record (shells filled, Dowker's complex with Atkin's Q-analysis, the time filtration). The domain writes no general method of its own: morphe gained persistent homology, cycles over ℤ/2, Q-analysis and canonical forms for it, and physis the exact contraction of ℤ₂ gauge theories, the heat bath of finite groups and the invariants of SU(2) (physis decision 0007). The canons are seeds drafted by a language model, unchecked; `data/products/history_check.py` checks them against Wikidata.
- **The Genome Atlas** (`applications/genome-atlas`, [ADR 0043](development/decisions/0043-the-genome-infrastructure.md)) draws genomes and computes with them, in nine views. Genome: the circular map and a zoomable strip, the origin of replication from the cumulative skew, and genes, CpG islands and repeats found from the sequence alone. Align: two sequences as a path through their edit graph, with the footprint of each schedule (the full matrix, the band, the WFA's wavefronts), every schedule held to Gotoh's optimum and morphe's generic sum, other semirings, the parametric polygon and its fan, and partial orders. Index: suffix arrays, the BWT and the FMD-index. Reads: a diploid read, mapped and called against its truth, a pileup browser, mapping accuracy from a job across the workers, SAM and VCF. Assembly: the k-mer spectrum, the de Bruijn graph cleaned step by step, contigs, the BEST theorem. Expression: pseudoalignment and the EM against the truth, and a two-condition study. Trees: gene trees by distances or likelihood with the bootstrap, dN/dS, whole genomes by their sketches. Populations: the coalescent, Wright–Fisher, PCA and admixture. Run: the machine, the workers, jobs on a timeline, threads, and a WebGPU kernel of Smith–Waterman checked against the CPU's three schedules. Fourteen RefSeq genomes come through the data library (`data/registry/ncbi-genomes.json`).
- **The genomics and evolution domains** ([ADR 0043](development/decisions/0043-the-genome-infrastructure.md)) are C++ from the start, header-only, each module checked against published examples and brute force: genomics' thirteen modules (sequence, simulate, kmers, align, index, sketch, mapping, assembly, variant, hmm, msa, expression, annotation) and evolution's three (phylogeny, population, selection), whose suites run 916 and 394 checks. Among their results: every CDS of fourteen RefSeq genomes translated as RefSeq translates it; Farrar's striped and inter-sequence Smith–Waterman at 26 and 31 GCUPS on eight threads; BiWFA in memory proportional to the score; a minimap2-like mapper placing *E. coli*'s reads 98.7% correctly; a caller with SNV precision 1.000 and recall 0.985 on a simulated diploid at 30×; a de Bruijn assembler giving *E. coli* at N50 20.6 kb; maximum likelihood under GTR + Γ recovering simulated trees. Their general mathematics went to the foundations: morphe's `formal` (decision 0005), and calx's lanes, task graphs, eigenproblem, gamma family and optimisation (decision 0006).
- **The Economy Atlas** (`applications/economy-atlas`, [ADR 0044](development/decisions/0044-value-in-space-and-time.md)) draws value carried between places and through time, in four views. Map: the world's trade in goods (OECD BIMTS, 2010–2024) split into reciprocal exchange and net flows, and the net flows, each link conducting its gross trade, into a driven part (up an attraction, the potential from which value is drawn to where it is drawn to), local eddies round the spherical Delaunay triangles of the economies, and global circulation through distant links. Spacetime: twenty currencies held at their policy rates (BIS; the euro at the ECB's deposit facility rate) and exchanged with the euro (ECB) month by month, every lens's curvature the departure from uncovered interest parity, measured against the curvature of the exchange rates alone. Production: an editable input–output table and its Leontief model, with the circulation of deliveries. Circuits: an editable transfer network with the examples of the user's paper on loop deficits.
- **The economics domain** ([ADR 0044](development/decisions/0044-value-in-space-and-time.md)) is C++ from the start, header-only, in five modules: `transfer` (holonomies, loop deficits, re-description, circuits, prices, retention), `circulation` (the decomposition of flows with conductances, exchange and net transfer, the pairing of flows and transports, stocks), `production` (Leontief), `spacetime` (morphe's worldline complex with holding and exchange transports, curvature, three-valued reach) and `places` (the faces of a map). Its theory note proves that prices exist exactly when every loop deficit vanishes and that only circulating flow meets a deficit, and says which claims of the user's papers are not yet precise. Its general mathematics went to the foundations: physis's discretisation and computational geometry (decision 0008), calx's conjugate gradients (decision 0007) and morphe's worldline complex, on which the history domain's construction now also stands.
- **The Celestial Atlas** (`applications/celestial-atlas`, [ADR 0011](development/decisions/0011-first-application-celestial-atlas.md)) is a web application with five views:
  - the Solar System;
  - the stellar neighbourhood, from HYG or Gaia;
  - galaxies and their relations;
  - the sky chart;
  - orbits.

  It shows any instant from 1849 to 2150, or live. The bodies carry maps turned by the IAU rotation models ([ADR 0016](development/decisions/0016-maps-of-bodies-and-iau-rotation.md)), with finer detail at close range ([ADR 0020](development/decisions/0020-detail-of-bodies-at-close-range.md)). The observer is chosen by clicking the Earth, and the two atlases link to each other's view of the same place and instant ([ADR 0021](development/decisions/0021-the-sky-from-a-place.md)).
- **The Earth Atlas** (`applications/earth-atlas`, [ADR 0014](development/decisions/0014-earth-atlas-separated-data-and-typed-layers.md)) shows the atmosphere, ocean, land and people on a globe or on equal-area and equirectangular maps, at one instant.
  - **Data and layers.** Its data are separate products, joined only in the application, and its layers are typed by the shape of their data.
  - **Forecasts.** Forecasts from GFS, AIGFS and AIGEFS can be compared (ADR 0015), with Pangu-Weather run here on the CPU ([ADR 0022](development/decisions/0022-ai-forecast-cyclones-and-disasters.md)) and our two models in the page: the barotropic vorticity and shallow-water equations ([ADR 0019](development/decisions/0019-shallow-water-model.md)).
  - **Cyclones and disasters.** Tropical cyclones are tracked in every forecast, beside GDACS's observed tracks; GDACS's alerts are a layer.
  - **The ocean beneath a place** comes from HYCOM, on demand ([ADR 0023](development/decisions/0023-ocean-column-and-data-on-demand.md)).
  - **The solid Earth** ([ADR 0024](development/decisions/0024-the-solid-earth.md)): plates and their motion, GPS velocities and daily records, Global CMT mechanisms, CRUST1.0 and sections.
  - **Plug-ins** ([ADR 0027](development/decisions/0027-earth-atlas-plugins.md), [ADR 0029](development/decisions/0029-plugins-in-python-and-cpp.md)): scripts loaded on request that derive from `Diagnostic`, `ForecastModel`, `DataAssimilation` or `Trajectories`, run on the library and are shown through the layers; in JavaScript in the page, or in Python or C++ in the server, under one contract; all four kinds exist in each language ([ADR 0033](development/decisions/0033-cpp-plugin-classes-completed.md)), the classes in `api/python/science/plugins.py` and `api/cpp/plugin.hpp`. The examples in `applications/earth-atlas/plugins/` are the guide.
  - **Imagery and relief.** Satellite imagery or relief forms the base ([ADR 0018](development/decisions/0018-fine-relief-at-close-range.md)).
  - **Climate:** Köppen–Geiger classes from 1901 to 2099, normals, climate diagrams, and our own classification ([ADR 0017](development/decisions/0017-climate-classes-and-normals.md)).
- **`data/` is a general data library** ([ADR 0012](development/decisions/0012-selectable-parameterised-datasets.md), ADR 0014, ADR 0023):
  - **Manifests** choose which data, and how much, from fixed files, TAP queries, ERDDAP, THREDDS (NCSS) and OPeNDAP subsets, or web requests.
  - **Products** are of twenty-one kinds, listed in an index that applications read.
  - **Data requested on demand** go through each application's server.
- **The Earth Atlas as a forecasting system** ([ADR 0031](development/decisions/0031-observations-scores-and-hazards.md)): observations (METAR, NDBC buoys, the user's sensors through `observation/sensors`), forecasts kept and scored on a schedule (`data/acquisition/schedule.py`, `applications/earth-atlas/schedule.json`), ECMWF's open forecasts (IFS, AIFS, ENS probabilities), hazard plug-ins (fire weather, heat, heavy rain, tsunami travel times, aftershocks, volcanic ash), people and regions with exposure, drainage, Sentinel-2 NDVI, and assimilation of real observations. The research arc (inference foundation, C++ toolchain, a primitive-equation model with a real assimilation cycle) is planned in `development/ideas/2026-10-23-the-research-arc.md` and deferred.
- **The generation of spacetime by field theories** ([ADR 0040](development/decisions/0040-the-generation-of-spacetime.md), physis decision 0006): one framework for every theory whose expansion generates discrete spacetimes, in five stages (the action; Feynman graphs as coloured graphs with exact series in N; their 2-complexes, pseudomanifolds and spin foams; the large-N or condensed phase; the continuum), for scalar, vector, matrix and tensor field theories and for group field theories (Boulatov, Ooguri, quartic melonic). The matrix model recovers 2D quantum gravity's γ = −1/2, the tensor models the branched polymers' γ = 1/2; spin-foam amplitudes sum to Turaev–Viro; GFT condensates give Friedmann's equation with a bounce at ρ_c; random quadrangulations and melonic graphs carry spectral and Hausdorff dimensions. The Complex Atlas walks the five stages in its Generation mode and shows spin foams vertex by vertex (boundary spin network, quantum tetrahedron, 6j against Ponzano–Regge).
- **Complexes, symmetry and fields** ([ADR 0039](development/decisions/0039-complexes-symmetry-and-fields.md)): morphe's first slice (85 checks): exact integers, the Smith normal form and homology with generators, presented and finite groups, and the topology of CW and simplicial complexes, π₁, surfaces, manifolds, collapses, dual 2-complexes and coloured graphs; calx's forward differentiation; physis's symmetry (Lie groups, exact and quantum recoupling), geometry (curvature of any metric), spacetime (exact solutions, Regge geometry, causal sets), field theory (lattice fields, Feynman graphs), gauge (lattice gauge theory on any 2-complex, exact two-dimensional Yang–Mills) and quantum gravity (Turaev–Viro, BF amplitudes, coloured tensor models, matrix models), 202 checks; algebraic geometry is planned (morphe decision 0004).
- **The first C++** ([ADR 0035](development/decisions/0035-the-mathematical-foundations.md), [ADR 0036](development/decisions/0036-spiking-networks.md)): calx's first numerical modules (61 checks), physis's units, dynamics and quantum theory of the open qubit (20), and the neuroscience domain's networks of spiking neurons, `snn` (112): twenty-two neuron models as hybrid systems, among them the classical models of excitability, each checked against the bifurcation of its class, and a quantum emitter followed along exact quantum trajectories, plasticity, a parallel engine giving the same spikes on any number of threads, JSON specifications and a command-line simulator, checked against closed forms, firing patterns and mean-field theory; Brunel's network of 12,500 neurons runs a second in 2.6 s on eight threads. The mathematics of the foundations is mapped in [docs/mathematics.md](docs/mathematics.md); further spiking models, quantum ones among them, are planned from a verified survey.
