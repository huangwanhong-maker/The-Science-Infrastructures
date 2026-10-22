# CLAUDE.md

This file guides Claude Code in the repository of The Science Infrastructures ("the infrastructure"). It records the programme's purpose, design philosophy, architecture and working conventions.

This file can be revised: when a decision recorded here changes, update this file in the same change.

- The detailed design is in [docs/architecture.md](docs/architecture.md).
- Decisions and the reasons for them are in [development/decisions/](development/decisions/). The present structure is [ADR 0026](development/decisions/0026-domains-as-libraries.md).
- Each foundation library under `foundations/` has its own CLAUDE.md, which governs work inside it.

## 1. Purpose and Scope

This is an integrative programme to learn, use and research all domains of science, in theory and in practice, as one connected body of knowledge.

- **Three tiers** ([ADR 0026](development/decisions/0026-domains-as-libraries.md)):
  - **Foundations:** calx (computation and numerical methods), morphe (mathematical structures and their equivalences), physis (geometry and theoretical physics), and a planned library for inference. They hold laws and methods that name no object in nature.
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
- **Further domains** join the same order as the work reaches them: chemistry, the life sciences, health, the social sciences, engineering.
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
| calx | `foundations/calx` | Retargetable computation stack; numerical methods on floating-point numbers and arrays | — |
| morphe | `foundations/morphe` | Mathematical structures and their equivalences: arithmetic, algebra, categories, topology, algebraic geometry, formal structures | — |
| physis | `foundations/physis` | Geometry and theoretical physics: units, geometry, discretisation, mechanics, gravitation, thermodynamics, continuum, chemistry, radiation; a provisional JavaScript stand-in for discretisation | calx, morphe |
| neuroscience | `domains/neuroscience` | A domain in its own repository (ADR 0028): the brain's reference spaces, parcellations, maps and the head | calx, physis (declared; no domain) |

The foundations and the neuroscience domain are git submodules. A domain or application that grows large, or gains collaborators of its own, may be extracted with its history into its own repository (ADR 0026, the graduation rule); neuroscience was created so at the user's request.

### 3.2 The Order of the Domains

```
height 4   climate
height 3   tides · navigation · air-sea
height 2   astronomy · geophysics · geography · oceanography · meteorology
height 1   geodesy
height 0   metrology · neuroscience
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
│   └── neuroscience/  reference spaces · parcellations · maps · the head (git submodule, ADR 0028)
├── data/              dataset registry · acquisition · formats · products · provenance
├── visualization/     scene model · renderers · styles (ADR 0005); provisional js/
├── observation/       instruments · field work · sensors
├── api/               the service and its clients: Python · C and C++ · JavaScript (ADR 0030)
├── applications/      one folder per application: celestial-atlas · earth-atlas · brain-atlas
├── research/          questions · experiments · papers
├── notebooks/         exploration and teaching
├── docs/              architecture · gallery (later: glossary · references)
├── development/       decisions · ideas · roadmap (ADR 0009) · tools
└── assets/            images · fonts · style tokens
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
- Records accepted so far: 0001–0033. Records 0001–0024 were written in the first laboratory, which the infrastructure absorbed; 0025 founded this repository.

### 6.4 Change Discipline

- Before restructuring directories, adding or renaming a domain, moving a foundation, or changing a convention in this file, propose the change and wait for the user's agreement.
- Keep each change to one concept.
- When adding a domain or module, update:
  - its README and `module.json`;
  - its parent README's table;
  - the layout (§5), if it adds a directory here;
  - the glossary, if it introduces notation.
- If this file no longer matches the repository, say so and propose a fix.

### 6.5 Checks

- `node development/tools/run-tests.mjs` checks the dependency rules, runs every suite (foundations, domains in their order, visualization, the data library's Python checks in `data/tests/`, the API, applications), then runs each domain alone with its down-set. Run it before every commit that touches code.
- A domain's own suite: `node domains/<name>/tests/run-tests.mjs`.
- Check the applications visually by driving headless Chrome over the DevTools protocol, with an absolute `--user-data-dir`.

### 6.6 Working with the Foundations

- **Submodule commits.** Edit a foundation in place under `foundations/<name>` and commit there first. Then commit the updated submodule pointer here. Commit only the pointer, by path, when other changes are staged. The same holds for the neuroscience domain under `domains/neuroscience`.
- **Commit order.** calx and morphe first, then physis, then neuroscience, then this repository.
- **Each foundation's own CLAUDE.md** governs work inside it.
- **Links between repositories.** Use absolute GitHub URLs.
- **Shared conventions.** The writing style, status vocabulary, and decision and idea conventions are shared by every repository. Change them in all of them together.
- **Pushing.** Push only with the user's agreement.

### 6.7 Collaboration Style

- The principal investigator is fluent in abstract algebra and formal language theory. Prefer precise mathematical formulations.
- This is also a learning project. When introducing physics or conventions specific to a domain, give the reasoning behind them and cite their source.
- Visualizations should look warm and humane, not like a generic dashboard.

## 7. Technology

- **Decided** ([ADR 0003](development/decisions/0003-cpp-computational-core.md)): C++ for computational code, in the foundations and the domains, wherever possible. Python or another suitable language for data acquisition, orchestration, notebooks, machine learning and visualization.
- **Provisional** (to be confirmed in roadmap Phase 0): C++20, CMake, nanobind for Python bindings, a C++ test framework, pytest, and Python 3.12 or later.
- **Provisional JavaScript.** Until the C++ exists, the domains' code is ES modules in their `js/`, as is physis's stand-in for discretisation (physis ADR 0003); the C++ must pass the same checks.
- **Web applications:**
  - Use ES modules without a build step.
  - Load three.js from jsDelivr, pinned to an exact version (0.170.0).
  - Serve locally with Python's standard library. Each application's server mounts `domains/`, `foundations/` and `visualization/` beside its page, and the data cache at `/data/`; it reports the application's data recipe (`recipe.json`) at `/setup` and runs it on the page's request.
  - Test with Node (v22 is installed).
  - Each application has its own local server: the Celestial Atlas on port 8770, the Earth Atlas on 8790, the Brain Atlas on 8810.
  - A public instance runs on the user's server, `eu_region_server` (173.212.213.207): the Celestial Atlas on port 9200, the Earth Atlas on 9201, in the tmux session `science`, served with `--host 0.0.0.0 --sibling`. The READMEs link to it; keep the links current if it moves.
- **Visual style:**
  - The Celestial Atlas uses Dusty Orchid (#CA8DBE) as its main colour, chosen by the user: plum nights, orchid accents and orchid-tinted paper for charts.
  - The Earth Atlas shares that palette for now (ADR 0014), with data in perceptually uniform colour maps (cmocean, and our own orchid map); the user may give it a colour of its own.
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

1. **The C++ namespaces of the domains**, and of the infrastructure as a whole.
2. **C++ toolchain** on Windows (MSVC, clang/LLVM or WSL), together with the build and packaging tools.
3. **The inference foundation:** its name and its first module.
4. **Automatic differentiation** in calx, which would be built as an IR transformation.
5. **Precision policy details** (calx).
6. **Notes format.** Markdown with LaTeX (the current default), Jupyter or Quarto; LaTeX for papers.
7. **Languages of the documents.** English, with abstracts in Chinese and Japanese as in the README, or more.
8. **Licence for prose.** MIT covers software; notes may suit CC BY 4.0.
9. **First renderer technology** for the scene model.
10. **The next domains,** to be discussed now that the structure is settled.
11. **The first laboratory's repository on GitHub:** archive it with a pointer here, or leave it.

## 9. Current State

As of 2026-10-22:

- **The structure of ADR 0026 is in place.** The first laboratory's history is joined to this repository, and the foundations are submodules under `foundations/`.
- **The domains.** Nine domains have code: metrology, geodesy, astronomy, geophysics, geography, oceanography, meteorology, climate and neuroscience. Tides, navigation and air–sea are planned. All the code is provisional JavaScript, with 380 checks against published references and closed-form answers. With the foundations, the edges and the applications, the whole run is 1051 checks in sixteen suites. Every domain passes its checks alone.
- **The Brain Atlas** (`applications/brain-atlas`, [ADR 0028](development/decisions/0028-neuroscience-and-the-brain-atlas.md), [ADR 0032](development/decisions/0032-the-anatomy-of-the-brain-atlas.md)) shows the human brain in the round from an anatomy product (`data/products/brain_anatomy.py`): the cortex by gyrus, the deep nuclei, ventricles and fibre bundles, every part with a card from the domain's knowledge base (152 structures, with references). It cuts the brain and shows sections, colours the cortex by Julich-Brain or Schaefer, maps data by structure and by coordinate, places the 10–20 electrodes on the scalp and carries their values onto the cortex, scales the brain to a subject's head, and shows the Allen mouse brain; five guided tours.
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
  - **Products** are of fourteen kinds, listed in an index that applications read.
  - **Data requested on demand** go through each application's server.
- **The Earth Atlas as a forecasting system** ([ADR 0031](development/decisions/0031-observations-scores-and-hazards.md)): observations (METAR, NDBC buoys, the user's sensors through `observation/sensors`), forecasts kept and scored on a schedule (`data/acquisition/schedule.py`, `applications/earth-atlas/schedule.json`), ECMWF's open forecasts (IFS, AIFS, ENS probabilities), hazard plug-ins (fire weather, heat, heavy rain, tsunami travel times, aftershocks, volcanic ash), people and regions with exposure, drainage, Sentinel-2 NDVI, and assimilation of real observations. The research arc (inference foundation, C++ toolchain, a primitive-equation model with a real assimilation cycle) is planned in `development/ideas/2026-10-23-the-research-arc.md` and deferred.
- **There is no C++ code yet.** The next step is roadmap Phase 0.
