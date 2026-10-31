# 0063. The District as a Complex to Work On: Its Cells, Buildings Among Them, the Freedom on Top and the Dynamics

- **Status:** Accepted, provisionally. Amends [ADR 0062](0062-landscapes-as-cellular-configurations.md).
- **Date:** 2026-10-31
- **Decided by:**
  - **The request.** After ADR 0062 was pushed, the user wrote: "I think the utilization of our model still wierd. We should can visualize the complex, and freedom on top of it.. And dynamics." And then: "I means even architecture can fully convert into complex, They are space..."
  - **What Claude did.** It read the first message as asking for three things: the complex shown as a complex, the degrees of freedom on its cells set by hand together with the complex rewritten at will, and the dynamics run on what results. It read the second as asking that the buildings be cells of that complex too. It wrote this record.

## 1. Context

The plug-in of ADR 0062 used the formalism through one scenario: two stores picked for the user, one lane between them, a day run from midnight. The complex was drawn as a map, not as a complex. Nothing on it could be set or changed except through that scenario. Its buildings were not cells: one building at a time could be filled with rooms, apart from the district.

## 2. Decision

### 2.1 The Complex Shown as a Complex

The layer "the district as a complex" draws the cells by dimension, in two views:
- **On the ground:**
  - each block (2-cell) drawn inset from its streets;
  - each segment (1-cell) darker the more blocks it bounds, with an arrow for its orientation;
  - junctions (0-cells) as posts, and holes tinted.
- **By dimension:** the junctions on the ground, the segments, blocks and buildings at rising heights. Each cell is joined to its faces by threads coloured by the sign of the incidence, the face poset of the complex standing over the map.

A cell picked shows its closure and its star.

### 2.2 Buildings Are Cells

The architecture domain's `space` gains `massing.hpp`: the buildings of a district as one regular 3-complex.
- Each building is a 3-cell bounded by its ground face, its roof and its walls.
- Two buildings standing against each other share the wall between them, one 2-cell with two cofaces, so that buildings join into terraces as streets join blocks.
- Corners are pooled, and T-junctions are inserted, so that shared stretches of wall are shared exactly.
- Walls are cut at every height where something begins or ends on their edge.
- Every face is walked so that each building is oriented outward.

`environment` gains `massing_heat.hpp`: the buildings' heat as the Hodge Laplacian on these 3-cells, through party walls between neighbours and to the air and the ground, warmed by the Sun's sol-air gains on facades and roofs, each step Crank–Nicolson (calx's θ-method). On the Hoddle Grid's middle: 440 buildings, 565 shared walls, 98 terraces, a regular complex of Euler characteristic 77.

### 2.3 The Freedom on Top

**The state on the cells.** What the cells carry at the start can be set: walkers at a junction, a block's warmth, a building's warmth, the stores open or shut, a rule in force.

**The cells themselves can be rewritten.** The urban-science domain gains free edits on the district, each a rewrite checked for regularity and refused with its reason:
- `cut_street` cuts a street at a point (a refinement);
- `open_passage` opens a passage through a block between two points of its boundary (an intervention);
- `join_blocks` joins the two blocks on either side of a street (an intervention).

The edits are kept as a list that every request carries, so the engine rebuilds the same complex with the same cell ids. The engine gives their order (a trace).

### 2.4 The Dynamics

The landscape runs on the edited complex from the state set, from the clock's hour (`LandscapeParams::start` and the stores' and permit's state at the start). The buildings' heat runs beside it over the same hours. The map draws both at the clock's hour: the walkers, and the buildings coloured against the air.

### 2.5 The Workspace

The analysis "The landscape and its events" is replaced by "The district: complex, freedom, dynamics". It shows:
- the complex, with its counts, homology and holes, its views, and the cell picked;
- the freedom: the state on the cells, and the edits with their history;
- the dynamics: a run, its record, and what could happen at each hour and what did.

### 2.6 Corrections Made Along the Way

- **The shell's layer race.** A layer switched off while it loads no longer leaves its drawing behind.
- **Frames by time.** A run may start at any hour, so frames are looked up by time, not by index from midnight.

## 3. Consequences

- **Checks:**
  - urban science: 73 landscape checks, among them the edits and a run from noon;
  - architecture: 33 checks, among them the massing's regularity, its shared walls and T-junctions, ∂∂ = 0 and the party wall's exchange of heat;
  - the workbench: 36.
- **Open:**
  - the ground beneath the buildings is not yet subdivided by their footprints, so the buildings and the blocks are one picture of space but two complexes, joined by where the buildings stand;
  - a courtyard is not carved, since a ground face with a hole is not a cell.
