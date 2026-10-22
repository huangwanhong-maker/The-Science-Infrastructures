# 0006. Theory Notes in a Subject Tree

- **Status:** Accepted; amended by [0008](0008-independent-libraries.md)
- **Date:** 2026-10-21
- **Decided by:** designed by Claude at the user's request; settles the open question of where theory notes live

## 1. Context

The initial skeleton had `theoretical astronomy/` alongside `core/`. Every subject needs theory notes, though, and many notes span several modules. For example, the equilibrium tide draws on potential theory, spherical harmonics and lunar ephemerides.

## 2. Decision

- **One tree, organised by subject.** All theory notes live under a top-level `theory/` directory, with subject directories `mathematics/`, `physics/`, `astronomy/`, `oceanography/`, `geodesy/` and `computation/`.
- **Two-way links between notes and code.** A note links to the code that implements it through its `implemented-by` front matter, and code links back to the note through source citations. The links do not depend on matching directory paths.

## 3. Alternatives Considered

- **Notes beside the code, in each module.** Rejected: notes that span modules would have no home. Notes also have a different audience and medium from code.
- **A `theory/` tree that mirrors `core/` exactly.** Rejected: knowledge is organised by subject, while code is organised by layer. For example, metrology is a layer of the code, not a subject.

## 4. Consequences

- `theory/README.md` defines the front-matter convention.
- Tooling can later check that the links resolve in both directions.
