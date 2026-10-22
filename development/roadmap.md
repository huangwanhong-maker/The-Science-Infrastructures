# Roadmap

*Status: draft, 2026-10-22.* The phases are ordered intentions, not commitments. Revise them as the work teaches us. Each laboratory keeps its own roadmap; this one covers the infrastructure as a whole.

## Phase 0. Foundations

- **This repository.** Done (2026-10-22): founded, with the first laboratory as a submodule ([ADR 0001](decisions/0001-an-umbrella-of-laboratories.md)).
- **The shared libraries.** Decide whether calx, morphe and physis move from the first laboratory to `lib/` here. Proposed; awaits the user.
- **The first laboratory's path.** Decide whether `astronomy and geometry` becomes a slug without spaces that names its scope.
- **Continuous integration.** Each repository runs its own tests; this repository checks that its submodules' recorded commits exist and pass.

## Phase 1. Shared Layers

- **The libraries** as the foundation of every laboratory: the first laboratory's roadmap Phase 0 (a C++ toolchain, a CMake superbuild, tests and bindings) carried out for them.
- **Shared edges.** Decide whether the first laboratory's data library (manifests, fetching, provenance, products, requests on demand) and its visualization layer become shared repositories. Every laboratory needs data and pictures.
- **A glossary and references** shared across laboratories, with each laboratory's own notation mapped to them.

## Phase 2. The Second Laboratory

- **Choose it.** Candidates:
  - **Life sciences.** Genomics and sequence analysis would meet morphe's formal structures, and ecology the first laboratory's climate.
  - **Physics and chemistry,** where physis would become a laboratory's foundation in full.
  - **Earth sciences, solid Earth,** beginning from the first laboratory's geophysics.
- **Found it** with its own repository, CLAUDE.md, decision log and first application, on the shared libraries.

## Phase 3. Bridges

- **The first bridges** between laboratories, chosen by where their structures and phenomena meet (see [the idea on shared structures](ideas/2026-10-22-shared-structures-across-the-sciences.md)). Examples: biogeochemistry, climate and health, the physics of the atmosphere and ocean.
