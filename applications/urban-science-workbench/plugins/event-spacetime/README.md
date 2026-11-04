# Event Spacetime

*A district's real events as a 2-complex, drawn above a city that can be seen through.* A plug-in of the [Urban Science Workbench](../../README.md). Status: **draft** (2026-11-04). Activity: research. [ADR 0067](../../../../development/decisions/0067-event-spacetime.md); the design and its formalism, [development/ideas/2026-11-04-event-based-gis-on-a-2-complex.md](../../../../development/ideas/2026-11-04-event-based-gis-on-a-2-complex.md), for the planned paper *Event-Based 3D GIS on a 2-Complex: A Spin-Foam Approach*.

## 1. Purpose

A geographic state is the realization of a history of events. The plug-in holds a district's record of real events as the complex its design defines (the design, §1–2):
- **Vertices are events.** The two ends of the record's window are vertices too.
- **Edges are couplings.** These are volumes: buildings, the space of a street between two junctions, the space of a junction. Each lives between two events it attends.
- **Faces are primitives.** These are interfaces between two volumes: party walls, frontages, the cross-sections where a street meets a junction. Each is a lens between two events that meet it, attached up one side's worldline and down the other's.

It draws that complex in spacetime above the district, with time rising, and the network of the moment at the clock's date. The city beneath can be made to fade, so that the complex shows through its towers.

## 2. The Record

Melbourne's Hoddle Grid, the middle of Landscape Dynamics' district (Elizabeth–Russell, Flinders–Lonsdale), 2002–2024. Every fact attaches to the cell it names by identifier. The data library compiles them (`data/products/city_events.py`, `city-events-melbourne.json`):

| Source | What it gives the record |
|---|---|
| CLUE's buildings, by `base_property_id` | each building's construction (a birth, when in the window), refurbishments (alterations) and its state at every census (an observation, a local tick) |
| The Development Activity Monitor | developments completed: a birth, or an alteration of a building that stood |
| The Victorian Heritage Register | registrations (an institution touching the building) and their amendments |
| Event permits, 2014–2018 | the streets' events, a state of the street spaces they name from their start to their end |

Births end nothing in this record: the buildings are those standing now (OpenStreetMap's), so a demolition before 2025 has no geometry to stand on. Melbourne's permits do not say whether a street was closed to vehicles, so they touch the street spaces and meet none of their interfaces. A record that says so (New York's closures, to come) meets the cross-sections, whose permeability it changes.

## 3. Layers and Analyses

**The district's events in spacetime.** Time rises 4, 10 or 20 m a year, from the ground or from above the tallest building. It shows:
- worldlines for the volumes, their vertices marked;
- sheets for the interfaces' lenses;
- beads for the events, coloured by kind;
- a plane at the clock's date with the network of what stood then, nodes and links.

Two ways to see one complex: refined (K̂, each event split into its volumes, joined by rungs), or as events (K, each event one vertex, the worldlines bending into it). Pick a bead, a worldline or a sheet for what it is.

**The record and its complex.** Six tabs:
- **The record:** what the record holds and what it left out, with its reasons and its sources.
- **The complex:** the cells and homology of K and K̂.
- **The moment:** the network at the clock's date.
- **What it remembers:** the buildings' uses and the premises' trades, fitted as chains on their classifications: each level's rate, its memory and how long it is half kept, against one level alone; and the buildings' chains on the event complex, the record's probability by the contraction against the product of their transitions.
- **Rates:** each kind of event per year of exposure.
- **The drawing:** how it is drawn, which cells, the city's opacity, and the buildings' worldlines coloured by their unseen changes of use.

## 4. The Engine

C++ on urban-science's `events` (the record, its complexes, the cut), its `landscape` (the district) and architecture's `space` (the massing), on morphe's worldline complexes (`engine/main.cpp`, `engine/ops.hpp`):
- **`record`.** Takes the district (its network and buildings, sent by the page) and the area's record. It returns the couplings, interfaces and events for drawing, the topology, and what was left out.
- **`cut`.** Returns the network of the moment at an instant.

The workbench's checks give the engine a district of three buildings and a street, with a record of their own, and hold every count to its value by hand.

## 5. What the District Remembers

The design's stages E3 and E4: physis's chains on a classified set and its contraction of a record (decision 0020), fitted to the censuses of 2002–2024.
- **The buildings' uses.** 404 properties, 8,455 transitions between censuses, 22 of CLUE's space uses in their families. A building half keeps its family (Retail, Unoccupied, Parking, …) for 7.3 ± 0.3 years and, within it, its exact use for 2.9 ± 0.3 years. One level alone is rejected, at twice the log-likelihood ratio 128 on one degree of freedom.
- **The premises' trades.** 213 premises that one establishment holds, 3,154 transitions, 68 ANZSIC classes. A premises that changes trade mostly draws a new one from the whole district's mix, at 0.093 ± 0.006 a year, so its division, subdivision and group are half kept for 7.4 ± 0.5 years. Within its group it changes class at a further 0.020 ± 0.007 a year, so its exact class is half kept for 6.1 ± 0.5 years. The division and the subdivision have no rate of their own, at the bound: a premises that changes hands keeps its group or loses everything. One level alone is rejected, at twice the log-likelihood ratio 17.9 on three degrees of freedom (p = 5 × 10⁻⁴).
- **On the complex.** The buildings' chains on the worldlines of K, contracted through every event they attend (8,272 cells), give log Z = −4191.535138, equal to the product of their transitions.

## 6. What Comes Next

- the selection rule at splits of properties (only six in the record);
- the posterior of unrecorded events drawn as vertices;
- New York's Financial District, with its closures and Open Streets as regimes of the interfaces;
- the gauge models of relational states.
