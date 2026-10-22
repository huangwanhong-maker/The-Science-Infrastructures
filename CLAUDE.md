# CLAUDE.md

This file guides Claude Code in the repository of The Science Infrastructures ("the infrastructure"). It records the programme's purpose, design philosophy, architecture and working conventions.

This file can be revised: when a decision recorded here changes, update this file in the same change.

- Decisions and the reasons for them are in [development/decisions/](development/decisions/).
- Each laboratory and each library has its own CLAUDE.md, which governs work inside it. Read it before working there.

## 1. Purpose and Scope

This is an integrative programme to learn, use and research all domains of science, in theory and in practice, as one connected body of knowledge.

- **Laboratories.** The work is organised as laboratories. Each covers one domain or a group of closely connected domains, lives in its own repository, and is attached here as a git submodule.
- **Shared foundations.** Mathematics, computation and physics are written once, in shared libraries, and used by every laboratory:
  - calx: computation and numerical methods;
  - morphe: mathematical structures and their equivalences;
  - physis: geometry and theoretical physics.
- **Bridges.** Topics that join laboratories, such as biogeochemistry, or climate and health, are bridges (planned).
- **This repository** holds what concerns all of them:
  - the map of the sciences and of the laboratories;
  - the shared conventions;
  - the decisions that cross laboratories;
  - the roadmap.

Every artefact serves one of three activities, and says which:

- **Learning:** derivations, worked examples and expositions, written to be understood.
- **Use:** tested tools whose numbers can be relied on.
- **Research:** open questions, experiments and new results, built on the first two.

*Rationale.* The sciences share their structures: symmetry and its representations, conservation laws and the partial differential equations that express them, networks, sequences and grammars, and inference under uncertainty. One infrastructure makes these connections explicit. It also lets one tested implementation of each method serve every field.

## 2. Design Philosophy

The governing principle is that of the first laboratory: **modular, extensible, revisable**.

- **Modularity.**
  - A laboratory is a module: it states its scope, declares its dependencies, and exposes its products (data products, tools, notes).
  - Laboratories use each other only through the shared libraries and these published products.
- **Extensibility.**
  - A new domain is a new laboratory, and a new general method is added to a shared library. Neither requires editing what exists.
  - Introduce an abstraction only after its pattern has appeared in two places.
- **Revisability.**
  - Every module, note and model carries a status: `planned`, `seed`, `draft`, `stable` or `deprecated`.
  - Decisions are recorded with their reasons and never rewritten.
  - Prefer reversible choices, and mark provisional ones.
- **Verifiability.**
  - Numerical work is checked against authoritative references.
  - Results that can be checked by computation are.

## 3. Architecture

### 3.1 Repositories

| Repository | Path here | Scope | Depends on |
|---|---|---|---|
| The Science Infrastructures | `.` | The map, shared conventions, cross-laboratory decisions, roadmap | every laboratory (as submodules) |
| Integrative Astronomy, Geography, Oceanography Laboratory | `astronomy and geometry/` | Astronomy, geography, oceanography, meteorology, geodesy (geophysics proposed) | calx, morphe, physis |
| calx, morphe, physis | today `astronomy and geometry/lib/` | Computation; mathematical structures; geometry and physics | calx and morphe: none; physis: both |

### 3.2 Dependency Rules

- Shared libraries depend on nothing outside themselves (physis on calx and morphe).
- A laboratory depends on the shared libraries, and on another laboratory only through its published products.
- No laboratory depends on this repository.
- An upward dependency is a design error. Fix it by moving code down, into a shared library, or by inverting the dependency.

### 3.3 Placement

- **A shared library** takes content that is general: mathematics, computation, physics that more than one laboratory needs.
- **A laboratory** takes content particular to its domains.
- **This repository** takes only what concerns all laboratories: their map, the conventions, cross-laboratory decisions and the roadmap. It holds no domain content.

## 4. Conventions

The conventions of the first laboratory's CLAUDE.md apply to every repository of the infrastructure, and change in all of them together:

- **Scientific conventions:** SI units, explicit time scales and reference frames.
- **Writing style:**
  - headings declarative, numbered and academic, never questions;
  - concepts before structure and usage;
  - asides labelled *Rationale.*, *Remark.* or *Caveat.*;
  - code identifiers, docstrings and comments in English.
- **Notes:** LaTeX mathematics, with hypotheses stated and proofs given or cited.
- **Status vocabulary:** the five statuses of §2.
- **Decision records and idea capture:** as in §5.

Laboratories may add conventions particular to their domains, such as the orientation conventions of geodesy or the nomenclature of chemistry.

## 5. Working Procedures

- **Decision log.** Record decisions that shape the infrastructure's structure, conventions or tooling in `development/decisions/NNNN-slug.md`, using the template in its index. Never rewrite a record; a new record supersedes or amends it. A decision that concerns one laboratory or library is recorded in that repository's own log. Records accepted so far: 0001.
- **Idea capture.** Write each new idea in `development/ideas/` as `YYYY-MM-DD-slug.md`, with status `seed`, and add it to the index. An idea that concerns one laboratory goes in that laboratory's own `development/ideas/`.
- **Change discipline.**
  - Before restructuring directories, adding or renaming a laboratory, moving a library, or changing a convention, propose the change and wait for the user's agreement.
  - Keep each change to one concept.
  - When adding a laboratory, update the README's table of laboratories, §3.1 of this file, and the roadmap.
- **Submodules.**
  - Work on a laboratory or library inside its own directory, and commit there first.
  - Then commit the updated submodule pointer here.
  - Commit in dependency order: libraries, then laboratories, then this repository.
- **Links between repositories.** Use absolute GitHub URLs: a relative link into a submodule does not resolve on GitHub.
- **Pushing.** Push only with the user's agreement.
- **Collaboration style.**
  - The principal investigator is fluent in abstract algebra and formal language theory; prefer precise mathematical formulations.
  - This is also a learning project. Give the reasoning behind conventions specific to a domain, and cite their sources.
  - Visualizations should look warm and humane, not like a generic dashboard.

## 6. Repository Layout

```
.
├── astronomy and geometry/   the Integrative Astronomy, Geography, Oceanography Laboratory (submodule)
├── development/              decisions · ideas · roadmap
├── CLAUDE.md
├── LICENSE
└── README.md
```

- Every directory has a README giving its purpose and contents.
- A directory is created when its first content is written. Until then it is a `planned` row in the README.
- The path `astronomy and geometry` contains spaces, so quote it.

## 7. Technology

- **Languages** (decided in the first laboratory, ADR 0003 there):
  - C++ for computational cores and the shared libraries;
  - Python for data acquisition, orchestration, notebooks and machine learning;
  - web applications as ES modules without a build step.
- **Licences.**
  - Every repository is MIT-licensed.
  - Check external libraries for licence compatibility before adopting them; copyleft libraries such as FFTW, FLINT and GMP need care.
  - Datasets keep their own licences, and are fetched into caches rather than committed.
- **Version control.**
  - This repository's remote is `https://github.com/huangwanhong-maker/The-Science-Infrastructures`.
  - The first laboratory's is `https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory`.
  - The libraries' are `https://github.com/huangwanhong-maker/{calx,morphe,physis}`.
  - Each repository fetches over HTTPS and pushes over the SSH alias `serendip_wanhong` (its `pushurl`).
- **Development machine** (as of 2026-10-22):
  - Windows 11, Intel i7-11370H, 32 GB of memory.
  - NVIDIA RTX 3050 Ti Laptop GPU (4 GB) and Intel Iris Xe.
  - Python 3.12, Node 22, Julia and WSL. No C++ compiler on the Windows side.

## 8. Open Decisions

1. **Where the shared libraries live.** calx, morphe and physis are general, but today they are submodules of the first laboratory. Raising them to this repository's `lib/` would let every laboratory share them. It changes the first laboratory's structure (its ADR 0008), so it awaits the user's agreement (ADR 0001 here).
2. **The second laboratory,** and the map of laboratories in the README.
3. **Geophysics.** The first laboratory proposes it as a field (its ADR 0024). It might instead begin an Earth sciences laboratory of its own.
4. **The first laboratory's path.** `astronomy and geometry` contains spaces and names only part of its scope. A slug such as `astronomy-geography-oceanography` would be safer.
5. **Shared edges.** The first laboratory's data library and visualization layer may serve every laboratory, and could become shared repositories.
6. **Languages of the documents.** English, with abstracts in Chinese and Japanese as in the READMEs, or more.

## 9. Current State

As of 2026-10-22:

- This repository is founded, with the first laboratory as a submodule ([ADR 0001](development/decisions/0025-an-umbrella-of-laboratories.md)).
- The first laboratory has decision records 0001–0024. It has two web applications, the Celestial Atlas (66 checks) and the Earth Atlas (151 checks), and a general data library. It has no C++ code yet.
- The shared libraries are scaffolded, inside the first laboratory.
