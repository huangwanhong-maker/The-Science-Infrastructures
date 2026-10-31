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
| [0002](0002-layered-architecture-with-pure-core.md) | Layered architecture with a pure core | Accepted; amended by 0008, 0010, 0026 |
| [0003](0003-cpp-computational-core.md) | C++ for the computational core | Accepted |
| [0004](0004-retargetable-computation-stack.md) | Retargetable computation stack | Accepted; amended by 0008 |
| [0005](0005-visualization-as-representation-layer.md) | Visualization as a representation layer | Accepted |
| [0006](0006-theory-notes-in-subject-tree.md) | Theory notes in a subject tree | Superseded by 0026 |
| [0007](0007-atmosphere-and-ocean-as-coupled-system.md) | Atmosphere and ocean as one coupled system | Accepted; amended by 0008, 0010 |
| [0008](0008-independent-libraries.md) | Independent libraries in separate repositories | Accepted; amended by 0026 |
| [0009](0009-development-folder.md) | A development folder for process material | Accepted |
| [0010](0010-domain-organised-fields-and-applications.md) | Domain-organised fields and applications | Superseded in part by 0026 |
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
| [0024](0024-the-solid-earth.md) | The solid Earth: plates, GPS, earthquake mechanisms, the crust and sections | Accepted, provisionally; geophysics a domain by 0026; its crust and sections in C++, with PREM beneath, by 0055 |
| [0025](0025-an-umbrella-of-laboratories.md) | An umbrella of laboratories for all the sciences (the umbrella's 0001) | Superseded in part by 0026 |
| [0026](0026-domains-as-libraries.md) | Domains as libraries, and works that compose them | Accepted |
| [0027](0027-earth-atlas-plugins.md) | Plug-ins of the Earth Atlas: algorithms loaded on request, shown through the layers | Accepted, provisionally |
| [0028](0028-neuroscience-and-the-brain-atlas.md) | Neuroscience as a domain, in a repository of its own, and the Brain Atlas | Accepted; amended by 0032, 0034 and 0036 |
| [0029](0029-plugins-in-python-and-cpp.md) | Plug-ins in Python and C++: the same contract, run by the server | Accepted, provisionally; amends 0027; amended by 0030 and 0033 |
| [0030](0030-the-api-edge.md) | The API: an edge through which other programs reach the infrastructure | Accepted, provisionally; amends 0029 |
| [0031](0031-observations-scores-and-hazards.md) | Observations, scores and hazards: the Earth Atlas as a forecasting system | Accepted, provisionally |
| [0032](0032-the-anatomy-of-the-brain-atlas.md) | The anatomy of the Brain Atlas: surfaces, a tree of structures and what is known of each | Accepted; amends 0028; amended by 0034 |
| [0033](0033-cpp-plugin-classes-completed.md) | The C++ plug-in classes completed | Accepted, provisionally; amends 0029 |
| [0034](0034-biosignals-on-the-brain.md) | Biosignals on the brain: EEG, MEG and fMRI carried onto the template | Accepted; amends 0028 and 0032 |
| [0035](0035-the-mathematical-foundations.md) | The mathematical foundations: a map, an order, and the first C++ | Accepted, provisionally |
| [0036](0036-spiking-networks.md) | Networks of spiking neurons, and the C++ among the checks | Accepted; amends 0028 |
| [0037](0037-release-notes.md) | Release notes at each commit | Accepted |
| [0038](0038-the-spiking-simulator.md) | The Spiking Simulator | Accepted |
| [0039](0039-complexes-symmetry-and-fields.md) | Complexes, symmetry and fields in the foundations; the Complex Atlas | Accepted, provisionally |
| [0040](0040-the-generation-of-spacetime.md) | The generation of spacetime by field theories; spin foams in the Complex Atlas | Accepted, provisionally |
| [0041](0041-the-history-domain.md) | The history domain and the History Complex | Accepted; amended by 0042 |
| [0042](0042-the-history-as-a-field-theory.md) | The history as a field theory, on the foundations | Accepted, provisionally; amends 0041 |
| [0043](0043-the-genome-infrastructure.md) | The genome infrastructure: the genomics and evolution domains and the Genome Atlas | Accepted; amended by 0045 and 0047 |
| [0044](0044-value-in-space-and-time.md) | Value in space and time: the economics domain and the Economy Atlas | Accepted, provisionally; amended by 0046 |
| [0045](0045-the-genome-workbench.md) | The Genome Workbench: the Genome Atlas renamed, its Toolbox, the user's files and a shell built from registries | Accepted, provisionally; amends 0043 |
| [0046](0046-the-economy-atlas-data-and-plugins.md) | The Economy Atlas: one's own data, plug-ins over networks of places, and a ledger in gold | Accepted, provisionally; amends 0044 |
| [0047](0047-group-grammars-and-advanced-genome-analysis.md) | Group grammars for genomes (the user's GE-CSG draft made exact), and a survey of advanced genome analysis | Accepted, provisionally; amends 0043 |
| [0049](0049-the-spacetime-of-the-earth-and-the-moon.md) | The spacetime of the Earth and the Moon: the domains' C++ of time, frames and figures, the navigation bridge, and the globe | Accepted, provisionally |
| [0050](0050-instruments-and-their-acquisition.md) | Instruments and their acquisition: addresses of events, the models instruments read, and scheduling with deadlines | Accepted, provisionally |
| [0051](0051-the-chemistry-workbench.md) | The Chemistry Workbench, and the manuscript on chemical structure generation made exact | Accepted, provisionally |
| [0052](0052-fields-known-in-part.md) | Fields known in part: potential fields, the Moon's magnetism, the planetary-science domain, the kernel asked, and instruments declared | Accepted, provisionally |
| [0053](0053-fields-and-their-visual-expression.md) | Fields and their visual expression: one model of fields, their laws apart, and the scene model | Accepted, provisionally |
| [0054](0054-bodies-dynamics-and-knowledge.md) | Bodies, their dynamics, and what is known: discrete mechanics, the robotics domain, and a mission that learns | Accepted, provisionally |
| [0055](0055-the-earth-beneath.md) | The Earth beneath: CRUST1.0 over PREM, sections along geodesics, and the world below the ground and the sea | Accepted, provisionally |
| [0057](0057-a-project-of-serendip-commons-society.md) | A project of Serendip Commons Society: stewardship without ownership, and the record of contributors | Accepted |
| [0058](0058-a-body-with-joints.md) | A body with joints: a rabbit on the Moon, and three methods of one motion | Accepted, provisionally |
| [0059](0059-robots-in-their-environments.md) | Robots in their environments: what a world presumes of a place, what a body feels there, and the fields at a point | Accepted, provisionally (the rabbit, the fish and the dog written) |
| [0060](0060-streets-as-open-systems.md) | Streets as open systems: street oceanography made exact, the urban-science domain and the Urban Building Workbench | Accepted, provisionally |

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
