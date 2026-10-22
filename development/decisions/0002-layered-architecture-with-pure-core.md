# 0002. Layered Architecture with a Pure Core

- **Status:** Accepted; amended by [0008](0008-independent-libraries.md) and [0010](0010-domain-organised-fields-and-applications.md)
- **Date:** 2026-10-21
- **Decided by:** designed by Claude at the user's request, following the discussion of 2026-10-21

## 1. Context

The initial skeleton grouped code by field. It had `core/astronomy`, `core/geometry` and `core/oceangraphy`, alongside `theoretical astronomy/`, `visualizations/`, `notebooks/` and `assets/`.

The discussion of 2026-10-21 widened the scope to include:

- a retargetable computation stack;
- mathematics and theoretical-physics libraries, including algebra, algebraic topology and geometry;
- visualization with real data;
- laboratories and applications;
- units, time and frames;
- data infrastructure;
- inference and machine learning;
- validation;
- a knowledge layer;
- observation;
- research management.

A field-by-field layout cannot hold all of this. The fields share their foundations, and the most interesting topics, such as tides and navigation, lie between fields.

## 2. Decision

1. **Layers.** Organise the project as layers with acyclic, downward dependencies, as described in [architecture §2–3](../../docs/architecture.md).
2. **A pure core.** `core/` does no network access, no instrument control, no user interface and keeps no hidden global state. It reads only files whose paths are passed in. Effects live at the edges: `data/`, `observation/`, `visualization/` and `applications/`.
3. **Layers inside the core.** In order: compute, foundations (numerics and mathematics), metrology, physics, domains and inference.
4. **Fields and bridges.** Both are siblings under `core/domains/`. A bridge, such as `tides`, declares the fields it joins.
5. **Directories are created on demand.** A directory is created when its first content is written. Until then, a module is a `planned` row in its parent README. The exception is `core/compute/`, whose component READMEs are its interface specification.
6. **Mapping from the initial skeleton.** All of the old directories were empty.

| Before | After |
|---|---|
| `core/astronomy/` | `core/domains/astronomy` (planned) |
| `core/geometry/` | `core/mathematics/geometry` (planned) |
| `core/oceangraphy/` | `core/domains/oceanography` (planned; spelling corrected) |
| `theoretical astronomy/` | `theory/` (see [0006](0006-theory-notes-in-subject-tree.md)) |
| `visualizations/` | `visualization/` (see [0005](0005-visualization-as-representation-layer.md)) |

## 3. Alternatives Considered

- **Fields at the top level** (`astronomy/`, `geometry/`, `oceanography/`, each with its own stack). Rejected: it duplicates the shared foundations and hides the bridges.
- **Languages at the top level** (`cpp/`, `python/`). Rejected: it splits one concept across two trees.
- **Everything in `core/`, including data access and rendering.** Rejected: the core would then carry network, file-system and graphics dependencies, which would stop it from being deployed to arbitrary devices.

## 4. Consequences

- The location of new code follows the placement guide ([architecture §7](../../docs/architecture.md)).
- Some content needs dependency inversion. For example, the ephemeris interface lives in the core, while file acquisition lives in `data/`.
- Paths inside `core/` are deeper (for example `core/domains/tides/`) than in a flat layout.
