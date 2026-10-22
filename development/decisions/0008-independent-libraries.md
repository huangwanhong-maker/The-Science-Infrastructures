# 0008. Independent Libraries in Separate Repositories

- **Status:** Accepted
- **Date:** 2026-10-21
- **Decided by:**
  - The user proposed separate repositories: first one for numerical computing and one for theoretical physics, then a third for algebraic structures.
  - Claude drew the boundaries.
  - The library names (calx, morphe, physis) are Claude's proposals, adopted provisionally.

## 1. Context

After [0002](0002-layered-architecture-with-pure-core.md), `core/` held seven layers: compute, numerics, mathematics, metrology, physics, domains and inference.

The lower layers know nothing of astronomy, the ocean or the atmosphere. They are useful beyond this project, for example in group-theoretic work on formal languages and sequences. They also change on a different rhythm from the laboratory's domain code. Inside one repository, their independence rested only on convention.

## 2. Decision

1. **Three libraries.** Three libraries, each in its own repository, are included here as git submodules under `lib/`:

| Library | Scope | Depends on |
|---|---|---|
| [calx](https://github.com/huangwanhong-maker/calx) | The retargetable computation stack; numerical methods on floating-point numbers and arrays | — |
| [morphe](https://github.com/huangwanhong-maker/morphe) | Exact arithmetic, algebra, category theory, topology, algebraic geometry, formal structures, equivalences | — |
| [physis](https://github.com/huangwanhong-maker/physis) | Units and constants; geometry as physics uses it; structure-preserving discretisation; mechanics, gravitation, thermodynamics, continuum mechanics, chemistry, radiation | calx, morphe |

2. **One rule draws the boundaries.**
   - Numbers and arrays go to calx.
   - Structures go to morphe.
   - Geometry and physical laws go to physis.
   - A method that needs both structures and numbers goes to physis.
3. **Independence.** calx and morphe are independent of each other. No library depends on the laboratory.
4. **What the laboratory keeps.** The laboratory keeps what is specific to its subject:
   - astronomical time scales, celestial reference frames and Earth orientation (`core/metrology`);
   - the domains and inference;
   - the edges;
   - validation, theory and research.
5. **Each library stands alone.** Each has its own README, CLAUDE.md, MIT licence, decision log, ideas, tests and CI, together with the theory notes for what it implements.
6. **No nested submodules.**
   - The laboratory provides all three libraries.
   - When physis is built on its own, it fetches calx and morphe through the build system.
7. **Mapping from the previous layout:**

| Before | After |
|---|---|
| `core/compute/` | `lib/calx/compute/` |
| `core/numerics/` (floating point) | `lib/calx/numerics/` |
| `core/numerics/`: exact integers and rationals | morphe, `arithmetic` |
| `core/numerics/`: discretisation | physis, `discretisation` |
| `core/mathematics/`: algebra, topology | morphe |
| `core/mathematics/`: geometry | physis, `geometry` |
| `core/physics/` | physis |
| `core/metrology/`: units, physical constants | physis, `units` |

8. **Earlier records amended.**
   - [0002](0002-layered-architecture-with-pure-core.md): the layers inside `core/`.
   - [0004](0004-retargetable-computation-stack.md): the stack now lives in calx, where its decision 0003 carries it forward.
   - [0006](0006-theory-notes-in-subject-tree.md): notes on the theory a library implements now travel with that library.
   - [0007](0007-atmosphere-and-ocean-as-coupled-system.md): the physics shared by ocean and atmosphere now lives in physis.

## 3. Alternatives Considered

- **One repository (the status quo).** Rejected. The boundaries rest only on convention, and the libraries cannot be reused alone.
- **Two libraries, with mathematics inside physics.** Rejected. Structures are broader than physics, and they are objects of research in their own right.
- **Splitting geometry from physics.** Rejected as premature.
- **Nested submodules** (physis containing calx and morphe). Rejected, because the laboratory would then hold two copies of each.

## 4. Consequences

- **Cross-repository changes.** A change that spans repositories is committed in dependency order: calx and morphe first, then physis, then the laboratory's submodule pointers.
- **Versioning.** Library versions stay at 0.x until their interfaces settle.
- **Provenance.** Each result records the commit of each library it used.
- **Shared conventions.** The writing style, status vocabulary, and the conventions for decisions and ideas are duplicated in four repositories. They must be changed in all four together.
- **Remote repositories.** All three must be created on GitHub before the laboratory can be cloned together with its submodules.
