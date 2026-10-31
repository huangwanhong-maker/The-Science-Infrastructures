# architecture

*Status: planned. Field, height 0.* Buildings as objects of a science: their form, the spaces they enclose and join, the frontage at which a building meets its street, and the rules that bound them. Founded as a domain of its own at the user's word (2026-10-31, [ADR 0060](../../development/decisions/0060-streets-as-open-systems.md) §2.8).

## 1. Purpose

A building is a thing in the world with models of its own. Its massing stands in the city's light and casts its shadows. Its plan is a complex of rooms, courts and passages, which a person reaches by paths and sees across. Its frontage is the interface through which a holding meets the street. Its form answers to rules: heights, setbacks, plot ratios. No other domain holds these.

Urban science, the domain of streets and cities as open systems, will stand on this one when a module of its own needs buildings: frontages in the terrain of a street, heights in the shadows of a district. The Urban Building Workbench's plug-in *Architecture* draws it.

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `form` | Buildings as solids and voids: footprints, heights and their bases, massing, setbacks, the envelope; the solids of a city model and what can be measured of them (volume, floor area, surface, compactness) | planned |
| `space` | Plans and districts as cell complexes (morphe's topology): rooms, courts and blocks as cells, their adjacency, the paths through them, visibility and reach (Hillier and Hanson 1984, *The Social Logic of Space*), computed on the complex | planned |
| `frontage` | The frontage as an interface between a holding and the street: openings, depths and thresholds measured at the holding, its physical component kept apart from its symbolic one | planned |
| `rules` | The rules that bound form: heights, setbacks, plot ratios, the envelopes they allow, and whether a building lies within them | planned |
| `environment` | A building in its environment: sun and shadow on its faces and the street, daylight, and later its heat and air (physis's continuum) | planned |

## 3. Dependencies

| Source | What architecture uses from it |
|---|---|
| morphe | `topology`: cell complexes, their homology and paths |
| physis | `geometry` (solids, frames), `continuum` (heat, later air) |
| calx | linear algebra, quadrature |

The Sun's position for `environment` will come from astronomy, declared when that module begins; civil engineering, when it is founded, will hold structures, and architecture will meet it there.

## 4. Used By

The Urban Building Workbench (its plug-in *Architecture*, planned), and urban science when it declares architecture.

## 5. Open Questions

- Where the line between architecture and civil engineering falls: form and space here, structures and their loads there.
- Which records of buildings come first: the city model's footprints and heights (OpenStreetMap; measured heights from the City of Melbourne's outlines), or plans of single buildings.
