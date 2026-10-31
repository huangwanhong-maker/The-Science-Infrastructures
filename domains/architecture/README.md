# architecture

*Status: draft. Field, height 0.* Buildings as objects of a science: their form, the spaces they enclose and join, the frontage at which a building meets its street, and the rules that bound them. Founded as a domain of its own at the user's word (2026-10-31, [ADR 0060](../../development/decisions/0060-streets-as-open-systems.md) §2.8); its first modules, `space` and `environment`, written for the formalism of landscapes ([ADR 0062](../../development/decisions/0062-landscapes-as-cellular-configurations.md)), C++ from the start, 17 checks.

## 1. Purpose

A building is a thing in the world with models of its own. Its massing stands in the city's light and casts its shadows. Its plan is a complex of rooms, courts and passages, which a person reaches by paths and sees across. Its frontage is the interface through which a holding meets the street. Its form answers to rules: heights, setbacks, plot ratios. No other domain holds these.

Urban science, the domain of streets and cities as open systems, will stand on this one when a module of its own needs buildings: frontages in the terrain of a street, heights in the shadows of a district. The Urban Building Workbench's plug-in *Architecture* draws it.

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `form` | Buildings as solids and voids: footprints, heights and their bases, massing, setbacks, the envelope; the solids of a city model and what can be measured of them (volume, floor area, surface, compactness) | planned |
| `space` | Plans and districts as cell complexes (morphe's topology): rooms, courts and blocks as cells, their adjacency, the paths through them, visibility and reach (Hillier and Hanson 1984, *The Social Logic of Space*), computed on the complex | draft: a building as a box of rooms fitted to its footprint, a regular complex (a ball), what each room carries; doors and stairs opened (the same cells read otherwise), rooms merged and parted (new cells, the heat kept); the rooms reached from the entrance (`rooms.hpp`) |
| `frontage` | The frontage as an interface between a holding and the street: openings, depths and thresholds measured at the holding, its physical component kept apart from its symbolic one | planned |
| `rules` | The rules that bound form: heights, setbacks, plot ratios, the envelopes they allow, and whether a building lies within them | planned |
| `environment` | A building in its environment: sun and shadow on its faces and the street, daylight, and later its heat and air (physis's continuum) | draft: the heat of its rooms as the Hodge Laplacian Δ₃ with the fabric in its stars, the Sun's sol-air gains through glass and opaque fabric, the building's thermal modes, each hour propagated exactly (`heat.hpp`); shadow and daylight planned |

## 3. Dependencies

| Source | What architecture uses from it |
|---|---|
| morphe | `topology`: regular complexes, cubical grids, their rewriting and homology; `formal`: footprints of rewrites |
| physis | `discretisation` (Hodge Laplacians with stars); later `geometry` (solids, frames) and `continuum` (air) |
| calx | linear algebra, matrix exponentials |

The Sun's altitude and azimuth are given to `environment` by its caller (the workbench takes them from astronomy), so the domain stays at height 0; civil engineering, when it is founded, will hold structures, and architecture will meet it there.

## 4. Used By

The Urban Building Workbench: its plug-in *Landscape Dynamics* fills a building picked on the map with rooms and follows its heat through the day (the analysis "A building's rooms"); its plug-in *Architecture* is planned. Urban science, when it declares architecture.

## 5. Open Questions

- Where the line between architecture and civil engineering falls: form and space here, structures and their loads there.
- Which records of buildings come first: the city model's footprints and heights (OpenStreetMap; measured heights from the City of Melbourne's outlines), or plans of single buildings.
