# Decision Log

Each record captures one decision that shapes the structure, conventions or tooling of The Science Infrastructures.

- Records are numbered and are never rewritten.
- To change a decision, write a new record that supersedes or amends the old one, and update the old record's status to point to the new one.

Decisions that concern only one library are recorded in that library's own `development/decisions/`.

*Remark.* The log joins two logs. Records 0001–0024 were written in the Integrative Astronomy, Geography, Oceanography Laboratory, where "the laboratory" means that repository. Record 0025 founded the umbrella repository as its own 0001, and was renumbered when the two histories were joined. The numbers 0001–0024 are kept because code and documents cite them.

## 1. Index

| No. | Decision | Status |
|---|---|---|
| [0001](0001-record-architecture-decisions.md) | Record architecture decisions | Accepted; amended by 0009 |
| [0002](0002-layered-architecture-with-pure-core.md) | Layered architecture with a pure core | Accepted; amended by 0008, 0010 |
| [0003](0003-cpp-computational-core.md) | C++ for the computational core | Accepted |
| [0004](0004-retargetable-computation-stack.md) | Retargetable computation stack | Accepted; amended by 0008 |
| [0005](0005-visualization-as-representation-layer.md) | Visualization as a representation layer | Accepted |
| [0006](0006-theory-notes-in-subject-tree.md) | Theory notes in a subject tree | Accepted; amended by 0008 |
| [0007](0007-atmosphere-and-ocean-as-coupled-system.md) | Atmosphere and ocean as one coupled system | Accepted; amended by 0008, 0010 |
| [0008](0008-independent-libraries.md) | Independent libraries in separate repositories | Accepted |
| [0009](0009-development-folder.md) | A development folder for process material | Accepted |
| [0010](0010-domain-organised-fields-and-applications.md) | Domain-organised fields and applications | Accepted |
| [0011](0011-first-application-celestial-atlas.md) | First application: Celestial Atlas, with provisional computation in JavaScript | Accepted |
| [0012](0012-selectable-parameterised-datasets.md) | Selectable, parameterised datasets | Accepted |
| [0013](0013-galaxies-and-their-relations.md) | Galaxies and their relations in the Celestial Atlas | Accepted |
| [0014](0014-earth-atlas-separated-data-and-typed-layers.md) | The Earth Atlas: separated data and typed layers | Accepted, provisionally; amended by 0015 |
| [0015](0015-earth-atlas-models-imagery-detail.md) | Forecast models, a laboratory model, imagery and detail in the Earth Atlas | Accepted, provisionally; amended by 0018 |
| [0016](0016-maps-of-bodies-and-iau-rotation.md) | Maps of the bodies, IAU rotation, and imagery shared across applications | Accepted, provisionally; amended by 0020 |
| [0017](0017-climate-classes-and-normals.md) | Climate classes, climate normals, and two new ways of acquiring data | Accepted, provisionally |
| [0018](0018-fine-relief-at-close-range.md) | Fine relief at close range, and a finer base map | Accepted, provisionally |
| [0019](0019-shallow-water-model.md) | A second laboratory model: the shallow-water equations | Accepted, provisionally |
| [0020](0020-detail-of-bodies-at-close-range.md) | Detail of the Moon, Mars and Mercury at close range | Accepted, provisionally |
| [0021](0021-the-sky-from-a-place.md) | The sky from a place: choosing where to stand, and links between the atlases | Accepted, provisionally |
| [0022](0022-ai-forecast-cyclones-and-disasters.md) | An AI weather model run in the laboratory, tropical cyclones tracked, and disaster alerts | Accepted, provisionally |
| [0023](0023-ocean-column-and-data-on-demand.md) | The ocean beneath a place, and data requested on demand | Accepted, provisionally |
| [0024](0024-the-solid-earth.md) | The solid Earth: plates, GPS, earthquake mechanisms, the crust and sections | Accepted, provisionally; geophysics as a field proposed |
| [0025](0025-an-umbrella-of-laboratories.md) | An umbrella of laboratories for all the sciences (the umbrella's 0001) | Accepted; raising the shared libraries proposed |

## 2. Template

```markdown
# NNNN. Title

- **Status:** Proposed | Accepted | Amended by NNNN | Superseded by NNNN
- **Date:** YYYY-MM-DD
- **Decided by:** …

## 1. Context
## 2. Decision
## 3. Alternatives Considered
## 4. Consequences
```
