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

**The record and its complex.** Four tabs:
- **The record:** what the record holds and what it left out, with its reasons and its sources.
- **The complex:** the cells and homology of K and K̂.
- **The moment:** the network at the clock's date.
- **The drawing:** how it is drawn, which cells, and the city's opacity.

## 4. The Engine

C++ on urban-science's `events` (the record, its complexes, the cut), its `landscape` (the district) and architecture's `space` (the massing), on morphe's worldline complexes (`engine/main.cpp`, `engine/ops.hpp`):
- **`record`.** Takes the district (its network and buildings, sent by the page) and the area's record. It returns the couplings, interfaces and events for drawing, the topology, and what was left out.
- **`cut`.** Returns the network of the moment at an instant.

The workbench's checks give the engine a district of three buildings and a street, with a record of their own, and hold every count to its value by hand.

## 5. What Comes Next

The amplitudes (the design, §4), in the order of the design's stages:
- **E3:** the event-coupled Markov model's contraction (physis);
- **E4:** the memories of each level of the classification of uses, fitted to CLUE's censuses; the selection rule at splits; the posterior of unrecorded changes;
- then the foams drawn on the complex.

After them: New York's Financial District, with its closures and Open Streets as regimes of the interfaces.
