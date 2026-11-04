# 0067. Event Spacetime: a City's Events as a 2-Complex, Urban Science as the Bridge over Architecture, and the Urban Science Workbench

- **Status:** Accepted, provisionally.
  - Amends [ADR 0060](0060-streets-as-open-systems.md): the workbench's name; urban science's place in the order of the domains (§2.8 there).
  - Builds on [ADR 0062](0062-landscapes-as-cellular-configurations.md) and [ADR 0063](0063-the-district-as-a-complex-to-work-on.md): the district, its events and its spacetime.
  - Builds on [ADR 0041](0041-the-history-domain.md) and [ADR 0042](0042-the-history-as-a-field-theory.md): the carrier and its spin foam.
- **Date:** 2026-11-04
- **Decided by:**
  - **The request.** On 2026-11-04 the user wrote:
    - "the name is not good. 'Urban Science Workbench' is better";
    - asked whether the spacetime and heterogeneous events can be visualised "utilize 2-Complex Spacetime structure (when visualize, the city transparancy can adjust to make theoretical physics models like spin foam more clear)";
    - "Try utilize real data as possible";
    - "I think we have dual representation … vertex are events … Faces are some primitive elements, Edge are coupling of such elements";
    - "Serious design first, and record the designed formalism somewhere", for a planned paper, *Event-Based 3D GIS On 2-Complex: A Spin-Foam Approach*, in which "we will present the full formalism we want, and what we really implement in engineering";
    - and pointed to their discussion of the paper.
  - **The design.** Claude surveyed the foundations' and the domains' complexes and spin foams and the open data of both cities. It then wrote the design, [the formalism and its engineering](../ideas/2026-11-04-event-based-gis-on-a-2-complex.md), and listed five choices.
  - **The user's choices:**
    - "Event Spacetime is good";
    - "volumes and interfaces";
    - "Melbourne";
    - the discussion moved beside the design and renamed;
    - on urban science and architecture: "B. In urban science."

## 1. Context

**The workbench already draws a district's day as spacetime.** Landscape Dynamics (ADR 0062, ADR 0063) does so:
- streets as slices at each hour;
- stores as worldlines, a lane as a sheet;
- events as beads, joined by their causal order.

It is a drawing of one model's run. It is not yet a structure in which any record of a city's events, of every kind, can be held, computed with and drawn.

**The infrastructure already has the carrier.** morphe's worldline complex (`topology/worldlines.hpp`) has:
- events as vertices;
- participants' segments as edges;
- the lenses of pairs as 2-cells.

The history complex (ADR 0041) and the economics domain's spacetime of transactions (ADR 0044) stand on it. The history domain also reads a record's signs as a $\mathbb Z_2$ spin foam, by the character expansion of a gauge-invariant weight (ADR 0042).

**The user's discussion asks open questions:**
- which complex;
- whether event categories need spins;
- what an amplitude means;
- whether events are more fundamental than the primitives they change.

The user noted that the last of these, taken literally, would make the engine generate complexes.

## 2. Decision

### 2.1 The Name

The Urban Building Workbench is the **Urban Science Workbench**, at `applications/urban-science-workbench`. Its port stays 8970, and its engines' cache stays `data/cache/engines/urban/`.

Current documents, manifests and the data library's registry use the new name. Records and release notes before this one keep the old name and paths, as they were written (§6.3 of CLAUDE.md).

### 2.2 The Formalism

The design note is the formalism, and the paper's core. In brief:

1. **The cells.**
   - Vertices are events.
   - Edges are couplings: persistent things that hold primitives together.
   - Faces are primitives: relations or interfaces between exactly two couplings, over the interval between two events that meet them.

   This is the spin foam's own reading: a face is the history of a link, an edge the history of a node. It is a CW complex, the cyclic order of a face's boundary being time.
   - A property of one thing is a state of its coupling.
   - A relation among three or more is held by a coupling of its own.
2. **The complexes.** The refined complex $\hat K$ and the event complex $K$ are morphe's refined and plain worldline complexes of the record. The spin network of the moment is a cut. Up to refinement, the complex is the spatial dual of the city's cellular spacetime (the note's Theorem 2.8).
3. **Spins are modes, not categories.**
   - A category is a rule: its footprint, its rewrite, its kernel, its rate.
   - A spin labels a mode of a cell's state under the symmetry of its dynamics, found from the rates.
4. **The amplitude is a likelihood.** It is the likelihood of the record under the *event-coupled Markov model*, in which cells interact only at events. Factorised over $K$ (Theorem 4.2) and expanded in modes, it is a sum over foams (Theorem 4.3), with:
   - face amplitude: a mode's survival over the face's duration;
   - vertex amplitude: the event's kernel between modes;
   - selection rules where the dynamics is symmetric.

   It is invariant under refinement (Theorem 4.5). Non-trivial foams end on observations (Proposition 4.6). Gibbs weights on relational states (the history domain's kind) are the second family, for later.
5. **Two levels for events.**
   - In the record, events are primary, and the geographic state at a date is its realization.
   - In the theory, they are the interactions of a field theory whose Feynman graphs are event complexes.
   - The engine reads records and never generates complexes to explain them. Unrecorded changes of a single cell are summed by its semigroup; latent interactions are order zero.

### 2.3 Volumes and Interfaces; Urban Science the Bridge over Architecture

**The network of the moment.** A city's spin network of the moment is the dual of its decomposition of space:
- couplings are volumes: buildings, street spaces, junction spaces;
- primitives are interfaces: party walls, frontages, the cross-sections of streets, doors.

**Urban science declares architecture.** Its `module.json` reads `"kind": "bridge"` and `"depends": ["architecture"]`, so:
- it stands at height 1;
- its new module `events` joins architecture's buildings (`space/massing.hpp`) to its own streets as one network, under its own checks;
- the root build adds architecture before urban science.

Architecture still depends on nothing but the foundations, and stays the domain of its own that ADR 0060 founded.

### 2.4 The Plug-in *Event Spacetime*

The workbench gains a plug-in, `plugins/event-spacetime/`, on the module `events`. It draws a district's record as its event complex in spacetime:
- **couplings** as worldlines rising from their volumes;
- **primitives** as ribbons between them, cut into lenses;
- **events** as beads at their place and time, coloured by category and marked by provenance;
- **the cut** at the clock's date: a plane carrying the spin network of the moment;
- **the city as realised** at that date;
- **a switch** between $\hat K$ (spacetime points) and $K$ (events);
- **the link of an event**, and the foams and the unrecorded changes read from the model.

The shell gains a control for every plug-in: the city's opacity, continuous from solid to a ghost, so that a model drawn above the city can be seen through it.

### 2.5 The General Methods

- **The contraction (physis decision 0020).** The contraction of an event-coupled Markov model on a worldline complex names no object in nature, so it goes to physis's `statistical_mechanics`:
  - forward and backward messages along time;
  - biorthogonal modes and their truncation;
  - the laws at cuts, and gradients.
- **The order of time moves down when it has a second user.** The history domain's order of a worldline complex's cells in time (`foam/transfer.hpp`, `time_order`) moves down beside the contraction when the gauge models of relational states (the design's §4.7) need it. The event-coupled Markov model needs no such order: K's vertices are its events in the order of time.
- **Isotypic projectors.** If morphe's `algebra` lacks the isotypic projectors of a permutation module, they are added there.

- **Triangles found through incidences.** morphe's refined worldline complex found its triangles by trying every triple of the pairs meeting at an event, O(m³). At a city's boundary event about two thousand interfaces meet, so it now finds them through the participants' incidences, the same triangles in the same order. Its check holds the new enumeration to the old on random records with parallel pairs.

### 2.6 Real Data, Melbourne First

The first district is Landscape Dynamics' middle of the Hoddle Grid, 2002–2024. The data library compiles its record (`data/products/city_events.py`). Couplings, states and observations:
- CLUE's buildings, by `base_property_id`, with use and floors at every census;
- CLUE's establishments, with their ANZSIC classes and their occupancies, re-fetched with `base_property_id`.

Events:
- the Development Activity Monitor's completions;
- building permits;
- road works;
- registrations on the Victorian Heritage Register;
- event permits, for the day's scale.

Observations by the hour: the pedestrian counts.

The Financial District follows, with buildings by BIN, streets by LION, closures for events, Open Streets, DOB filings, LPC designations and the storefront reports. Every dataset was checked for its identifiers and its span on 2026-11-04 (the note, §9).

## 3. Alternatives Considered

- **Keeping the name.** Rejected at the user's word: the workbench holds a science of the city, not only its buildings.
- **The walking network as the dictionary** (junctions as couplings, segments as primitives). Rejected for the first implementation. It holds streets but not buildings, whereas volumes and interfaces hold both in one network, as loop quantum gravity's spin networks hold volume and area. The walking network remains the 1-skeleton dual to it.
- **The join in the plug-in, the domains independent.** It would leave the network of volumes and interfaces tested only by the workbench and unusable by other works. Rejected by the user ("B").
- **Faces as facts of any arity** (the social spin foam). Rejected, because the cyclic order of a face's boundary would then mean nothing, as the history domain's idea note found.
- **Spins encoding event categories.** Rejected, because it would give names the meaning of representations.
- **Stipulated weights on faces, edges and vertices.** Replaced by a factorisation proved for a declared model.
- **Complex amplitudes with a quantum reading.** Kept out of scope. Complex coefficients arise from cyclic states (the hour of a daily crest), and $Z$ stays a probability.
- **An engine that generates complexes from interactions.** Rejected: the record is read, and the generative level is estimated from its rates.

## 4. Consequences

- **The order of the domains.** Urban science stands at height 1 in the order of the domains (CLAUDE.md §3.2, the README's and `domains/README.md`'s diagrams). Its checks run alone with architecture in their down-set (`check-alone.mjs`).
- **New in the domain.** The domain gains the module `events`; physis gains the contraction of decision 0020.
- **The workbench.** It gains *Event Spacetime* and the city's opacity.
- **The paper.** The paper's formalism is the design note. Its implementation is what the note's §8 marks as built, and its numbers come only from the engine.
- **Proposed to the history domain, not made.** A lens's coupling should depend on its duration, for invariance under refinement (the note, §4.4).
- **The order of the work.**
  1. the data;
  2. the module and its checks;
  3. the contraction;
  4. the fit;
  5. the plug-in.

  Each stage is checked as the note's §13 says.
