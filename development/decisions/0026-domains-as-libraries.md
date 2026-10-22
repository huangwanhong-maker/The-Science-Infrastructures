# 0026. Domains as Libraries, and Works That Compose Them

- **Status:** Accepted
- **Date:** 2026-10-22
- **Decided by:** the user, in discussion. They agreed to raise the shared libraries and to rename the first laboratory's folder. They then asked that each domain's core be separated, so that it could even run on its own as a library, with geophysics a domain in its own right and applications, experiments and studies integrative. Of the proposal that followed, they chose:
  - one repository for the domains, the shared edges and the works;
  - the names `foundations/`, `domains/`, `applications/` and `research/`;
  - a single decision log;
  - to drop the word "laboratory".

  Claude proposed the design and wrote this record.

## 1. Context

[0025](0025-an-umbrella-of-laboratories.md) made this repository an umbrella of laboratories:

- each laboratory a submodule grouping several domains;
- the first grouping astronomy, geography, oceanography, meteorology and geodesy.

[0024](0024-the-solid-earth.md) then proposed geophysics as a field of that laboratory, and the question of where geophysics belonged exposed a flaw. A laboratory *partitions* the domains, each domain belonging to exactly one laboratory. The sciences do not partition: geodesy serves astronomy and geophysics alike, and time scales serve every science that timestamps its data.

Two kinds of content were mixed in the laboratory:

- **Knowledge, organised by subject:** time scales, plate kinematics, the ocean column. It is pure, reusable, and belongs to one science.
- **Work, organised by purpose:** the Earth Atlas, a forecast experiment, a study. It is integrative, drawing on several sciences, and it touches the world.

The knowledge itself was scattered:

- notes were in `theory/`;
- planned C++ modules were in `core/`;
- the only working code, provisional JavaScript, was inside the two applications' `src/core/`.

Some of it was written twice: the Sun's position, and places on the Earth.

## 2. Decision

1. **Three tiers, and shared edges.**
   - **Foundations** (`foundations/`): laws and methods that name no object in nature. These are calx, morphe and physis, each its own repository attached as a submodule, with an inference library planned.
   - **Domains** (`domains/`): one library per science. Each holds a science's objects and their models, and is either a field or a bridge.
   - **Works**, organised by purpose and composing domains: `applications/`, `research/` and `notebooks/`.
   - **Shared edges**, which touch the world for every work: `data/`, `visualization/` and `observation/`.
   - **The test of placement:** content that names a thing in nature belongs to a domain.
2. **The order of the domains.**
   - Each domain declares in `module.json` the domains and foundations it may use. The relation "may depend on" must be acyclic, so that its closure is a partial order.
   - A domain's library is its down-set.
   - A bridge is a domain with several fields beneath it.
   - Geophysics is a field, above geodesy, accepting the proposal of 0024.
3. **The contract of a domain:**
   - a README;
   - `module.json`;
   - its provisional JavaScript in `js/`, behind the public interface `js/index.js`;
   - its tests in `tests/`, with their reference values;
   - its notes in `theory/`;
   - later, its C++.

   A domain is pure: no network access, no instrument I/O, no user interface, no hidden global state.
4. **Independence is checked, not hoped for.** Three tools in `development/tools/` enforce it:
   - `check-dependencies.mjs` verifies that the order is acyclic, that every import reaches another unit only through its public interface and only if declared, and that domains and foundations are pure.
   - `check-alone.mjs` copies each domain with nothing but its down-set into an empty directory, and runs its checks there.
   - `run-tests.mjs` runs everything.
5. **One repository for the domains, the edges and the works; the foundations stay apart.**
   - *Remark.* Independence is a property of the dependency graph and the build, not of the repository boundary.
   - **Graduation rule.** A domain or application that grows large, or gains collaborators of its own, may be extracted with its history into a repository of its own and attached as a submodule. The move is reversible.
6. **The first laboratory joins this repository.**
   - Its history is merged at the root, so that every file keeps its past through the moves.
   - Its own repository on GitHub is superseded; whether to archive it is the user's decision.
7. **One decision log.** The laboratory's records keep their numbers, 0001–0024, because code and documents cite them. This repository's 0001 becomes 0025, and this record is 0026.
8. **No laboratories.**
   - The word is retired as the name of what we are.
   - A themed grouping of works, if one is ever needed, is a *programme*: a README and a roadmap, owning no code.
   - Proper names keep the word (the Nevada Geodetic Laboratory), and so do the records, which are history.
9. **Where the provisional JavaScript goes.**
   - General numerics go down to physis, as its JavaScript stand-in ([physis ADR 0003](https://github.com/huangwanhong-maker/physis/blob/main/development/decisions/0003-provisional-javascript-stand-in.md)):
     - the spectral transform on the sphere;
     - longitude-latitude grids.
   - Drawing helpers that any application could use go to `visualization/js/`: contours, particles, display units and map tiles.
   - What only one application needs stays with it: the readers of its products, and the fetching of its data.
10. **Names.** Paths are lowercase and hyphenated, without spaces (`air_sea` becomes `air-sea`).

## 3. Alternatives Considered

- **Laboratories as the unit (0025).** This is the partition the sciences refuse, and every shared domain would be owned by one laboratory and borrowed by the others.
- **A repository per domain.**
  - It would give strong boundaries and versions of their own.
  - But it would mean eleven repositories today, most of them empty.
  - Every change across domains would be ordered commits in several repositories, followed by pointer updates. That friction was already felt with four.
  - The checks of §2.4 give the independence without it.
- **The foundations in the same repository.** They are general enough to be used outside this programme, and have their own conventions and logs.
- **Leaving the JavaScript in the applications until the C++ exists.** The domains would have stayed names without content. The duplication would also have stayed invisible.

## 4. Consequences

- **Earlier records.** This record changes the following:

  | Record | Change |
  |---|---|
  | [0002](0002-layered-architecture-with-pure-core.md) | Amended: the layers of `core/` become the foundations and the order of the domains. The pure core and the effectful edges remain. |
  | [0006](0006-theory-notes-in-subject-tree.md) | Superseded: notes live in each domain's `theory/`, beside its code. |
  | [0008](0008-independent-libraries.md) | Amended: the libraries move from `lib/` to `foundations/`; what was the laboratory's is the domains'. |
  | [0010](0010-domain-organised-fields-and-applications.md) | Superseded in part: the fields and bridges are the domains; applications have one directory each, not one per field. |
  | [0024](0024-the-solid-earth.md) | Its proposal accepted: geophysics is a domain. |
  | [0025](0025-an-umbrella-of-laboratories.md) | Superseded in part: there are no laboratories; the libraries are raised as it proposed. |

- **Paths in earlier records** are as they were when written. The map:

  | Before | After |
  |---|---|
  | `lib/calx`, `lib/morphe`, `lib/physis` | `foundations/calx`, `foundations/morphe`, `foundations/physis` |
  | `core/metrology/` | `domains/metrology/` |
  | `core/domains/<field or bridge>/` | `domains/<name>/` (`air_sea` as `air-sea`) |
  | `core/inference/` | the planned inference foundation |
  | `theory/<subject>/` | `domains/<subject>/theory/`; plate kinematics in geophysics, particles in `visualization/theory/` |
  | `validation/` | each domain's `tests/` and `tests/fixtures/` |
  | `applications/<field>/<application>/` | `applications/<application>/` |
  | `applications/*/src/core/` | the domains' `js/`, physis's `js/`, `visualization/js/`, and `applications/earth-atlas/src/data/fields.js` |
  | `applications/earth-atlas/src/layers/labmodel.js` | `applications/earth-atlas/src/layers/ourmodels.js` |

- **The checks.**
  - The 217 checks of the two atlases now run in twelve suites, with one more check added: physis refuses an operation without the radius of its sphere.
  - Every domain passes its checks alone with its down-set.
- **The local servers** mount `domains/`, `foundations/` and `visualization/` beside each application's page, so that the browser loads the same modules as Node does.
- *Caveat.*
  - The domains' code is still provisional JavaScript.
  - Their independence holds for that code. The C++ must keep it, through the same rules and the build.
  - The Sun is still computed twice: once from DE440, and once by a low-precision theory. This is an open question in astronomy's README.
