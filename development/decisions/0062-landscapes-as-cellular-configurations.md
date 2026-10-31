# 0062. Landscapes as Cellular Configurations: Their Events, Their Composition and Their Spacetime

- **Status:** Accepted, provisionally. The rates of the running landscape are declared, and the plug-in is a draft.
- **Date:** 2026-10-31
- **Decided by:**
  - **The request.** The user gave Claude their manuscript on the dynamics of urban landscapes in a unified framework (v0.9, then v2.3), and asked: "Try implementing the idea, and visualization. I think some module may resuable, like CW-complex, QFT, etc. Be carefully, high-quality, comprehensive." Then: "Be high-quality, and make efforts. UI also need well design."
  - **The user's correction, during the work.** "And not need reproduce the paper examples, or mention the paper's details. What we need is implement the formalism in the system, also visualizations, UI, etc.. For the visualization, the spacetime representation and emergence philosophy in this paper can also adopt." The work was redirected from the manuscript's worked examples to the formalism in general, on Melbourne's real district, with the district drawn as a spacetime.
  - **What Claude did.** It designed the formalism's implementation and wrote this record. It had reviewed an earlier draft (v0.5) for the user before this work began.

## 1. Context

The urban-science domain (ADR 0060) holds streets as open systems: their holdings, their rhythm, their currents. It had no account of the place itself changing: a passage cut through a block, a rule coming into force and being carried out, a reading of a place becoming part of what the place is. Nor did it have one of a building as a complex of rooms.

The foundations had most of what such an account needs:
- morphe: CW complexes and graph rewriting by double pushouts (its decision 0007);
- physis: master equations, entropy production and phase-type passage (its decision 0016), and the Hodge decomposition of flows;
- calx: the θ-method and matrix exponentials.

They lacked six things:
- regular complexes decided by their face posets, and rewriting that keeps a complex regular;
- the causal order of a record of events;
- Möbius inversion over a poset;
- the composition of open stochastic systems at interfaces;
- lumpability, with its error when inexact;
- piecewise-deterministic processes, exact linear elimination (memory), Doi–Peliti operators and weighted Hodge Laplacians of any degree.

## 2. Decision

### 2.1 The Standing of the Formalism

The formalism is implemented as mathematics in the foundations and as a module of the domain. The domain's theory note, `domains/urban-science/theory/landscape-dynamics.md`, states it self-contained: definitions, results with proofs or sources, and the line between what holds of any district and what is a model with declared rates. The manuscript's examples are not reproduced, and its apparatus does not appear in the page.

### 2.2 The General Methods, in the Foundations

- **morphe** (its decision 0008):
  - regular complexes decided by their face posets, exactly up to dimension 3, with orientation, chain complexes and cubical grids (`topology/regular_complex.hpp`);
  - their rewriting on incidence graphs, with rules built at a site and a regularity test after each rewrite (`topology/complex_rewriting.hpp`);
  - traces: footprints, the causal order of a record, its Foata steps and its linearisations (`formal/traces.hpp`);
  - Möbius inversion and anchored interactions (`combinatorics/mobius.hpp`).
- **physis** (its decision 0018):
  - the composition of open systems by fibre products at relations, joint mechanisms added once (`statistical_mechanics/composition.hpp`);
  - strong and eventwise lumpability, aggregation by the stationary law, its outward discrepancy and the coupling bound (`lumping.hpp`);
  - anchored jump processes (`jump_process.hpp`);
  - piecewise-deterministic Markov processes with forced jumps (`dynamics/piecewise_deterministic.hpp`);
  - exact linear elimination and its memory kernel (`dynamics/elimination.hpp`);
  - Doi–Peliti operators (`chemistry/doi_peliti.hpp`);
  - Hodge Laplacians of any degree with diagonal stars (`discretisation/hodge_laplacian.hpp`);
  - the faces of a plane graph (`computational_geometry/planar.hpp`).

The domain writes none of these again.

### 2.3 The Domain's Module

The module `landscape` of urban science holds what names a place:
- **District** (`district.hpp`): the district as a typed cellular configuration, its schema and the transfer of state, and the lane as three rewrites.
- **Running landscape** (`process.hpp`): the walkers, the blocks' warmth, the stores and the institutions as one piecewise-deterministic process on the strata of the complex, with its record, its registry of enabled events and each street's hops each way.
- **Event skeleton** (`skeleton.hpp`): five parts composed, with:
  - its law and irreversibility;
  - what coarser descriptions forget;
  - interventions and their interactions;
  - evidence.
- **Resources** (`value.hpp`): resources carried and converted, with the lattice of conserved weightings.

Its checks number 62.

### 2.4 Buildings, in Architecture

A building as a complex of rooms belongs to the architecture domain, which the user made a domain of its own (ADR 0060 §2.8: "Only architecture independent in domain"). Its first two modules are written here, on the foundations alone, at height 0:
- `space` (`rooms.hpp`): a box of rooms fitted to a building, its rooms the 3-cells, walls, floors and facades the 2-cells, what each room carries, and its events (doors and stairs opened, rooms merged and parted), with the rooms reached from the entrance;
- `environment` (`heat.hpp`): the heat of its rooms as the Hodge Laplacian Δ₃ with the fabric in its stars, the Sun's sol-air gains, its thermal modes, each hour propagated exactly.

Their checks number 17. The architecture domain's status moves from planned to draft; `form`, `frontage` and `rules` remain planned.

### 2.5 Meaning as State

The formalism distinguishes three things about the lane: a rule in force (the permit), the rule carried out (the lane built) and a meaning adopted (the lane recognised). The last is a state that changes the set of possible events, not only their rates. Two results make it precise:
- no coarser description can forget it;
- when it is drawn by use, the whole process is driven, with no detailed balance.

### 2.6 The Plug-in and Its Spacetime

The Urban Building Workbench gains its second plug-in, *Landscape Dynamics*, with an engine in C++ on urban science, architecture and the foundations.
- **Its layers:**
  - the district as cells;
  - the landscape running, with walkers walking each way at their counts;
  - the district's spacetime;
  - a building's rooms.
- **Its analyses:**
  - the landscape and its events;
  - the landscape's events, exactly;
  - a building's rooms;
  - scales and memory.

The spacetime layer draws the emergence the user named. Time rises above the district. Each hour the streets are a slice coloured by their traffic, so the day's tide is a pattern of the whole (dynamical). The stores' worldlines thicken while open, and the lane's sheet begins in mid-air with the junctions it made (structural). The lane deepens once recognised (interpretive). The events are beads joined by the threads of their causal order, and a plane marks the clock's hour.

### 2.7 What the Work Corrected

- **The lane planned again.** Closing a lane merges the streets it cut into new segments. A lane built again later in the same day, from a plan made at the start, asked for a cell that no longer existed. The process now plans again after every rewrite, `build_lane` refuses a plan made for another configuration, and closing restores each street's polyline.
- **Net flow and traffic.** The net flow of a street, which nearly cancels when walkers go both ways, had been drawn as if it were its traffic. The frames now carry each way's hops.
- **Sunlight.** Sunlight absorbed by opaque fabric was all put into the rooms. It is now the sol-air gain, and the rooms stay within a few degrees of the air.
- **Memory.** Memory on the real district was computed with dense exponentials of about 700 variables. It is now the θ-method on the detailed-balance form, in under a second.

## 3. Consequences

- The domain can say exactly what a rewrite conserves, what a description forgets, and how far a coarse chain drifts, for any district.
- The spacetime of a district is a view any later work on change in a place can draw on.
- The rates of the running landscape are declared, not learnt. Learning them from records with unseen changes is open (the theory note, §11).
- CLAUDE.md §1, §5 and §9, the READMEs of urban science, architecture and the workbench record the modules and the plug-in.
