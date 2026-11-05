# 0075. Agents in the City: Walkers, Households, Owners and Occupiers, Simulated and Drawn

- **Status:** Accepted, provisionally (2026-11-06). Amends ADR 0071, whose planned `movement`, `housing` and `occupancy` it builds.
- **Date:** 2026-11-05 to 2026-11-06
- **Decided by:**
  - **The request.** The user asked whether agent-based models are "also simulated and visualized", and then: "I think we should try ABM in urban science with several tasks, and visualize it."
  - **The proposal.** Claude proposed building the urban-science modules ADR 0071 planned on nomos's `agents` and tyche's `decision`, following the design in [development/ideas/2026-11-04-agents-and-decisions-for-urban-planning.md](../ideas/2026-11-04-agents-and-decisions-for-urban-planning.md) §§8–11. Each module would carry several tasks in a case study, and the walkers would be drawn in the Urban Science Workbench's planned *Human flow* plug-in.
  - **The work.** Parallel sessions wrote the three modules, their case studies and the plug-in on that design. Claude integrated them.

## 1. Context

ADR 0071 gave the infrastructure agent-based models (nomos `agents`) and decisions under partial observation (tyche `decision`), and built the first urban model on them (`maintenance`). It planned three more:
- walkers on the pedestrian network (`movement`);
- households among the centre's dwellings (`housing`);
- owners and occupiers of the district's buildings (`occupancy`).

The workbench listed *Human flow* as a planned plug-in. The district-economy case studies of ADR 0073 had calibrated a logit over walks to the City's counters. Agents choosing their way by the same law have that flow as their mean, so the macroscopic and microscopic models can be held to each other exactly.

## 2. Decision

### 2.1 The Domain Declares nomos

The urban-science domain adds nomos to its foundations (`module.json`) for its agents. The three modules are targets `urban::movement`, `urban::housing` and `urban::occupancy`.

### 2.2 `movement`: Walkers

**The model.**
- Walkers are nomos agents with an origin, a target and a way chosen junction by junction by the logit over walks. The law is read as a Markov chain P_ab = M_ab z_b / z_a, so the walkers' mean is the logit flow exactly.
- Walking speeds follow Weidmann, slowed by the density ahead.
- Admission to a segment is limited by Godunov's supply, and queues stand behind a segment's end.
- Signals run in two phases, and Fruin's levels of service are tallied.
- Arrivals are calibrated to the counters by hour and direction (tyche's NNLS).

**The street's regime.**
- A street can be open, one-way, or closed to walkers passing through.
- The choice is a POMDP solved online by tyche's POMCP, with the walkers themselves as the generative model.
- A run is written as urban-science's record with provenance "simulated".

24 checks. Among them:
- the walkers against the logit flow;
- conservation;
- runs identical on any number of threads;
- a lane's capacity and a signal's mean wait;
- POMCP against HSVI on an exactly solvable corridor.

### 2.3 `housing`: Households

**The model.**
- Classes of households (income, tenure, size) choose segments (CLUE type by small area) by a logit.
- Rents clear the market exactly. Berry's inversion recovers the segments' values, and Roy's identity gives the surplus.
- Agents move, arrive and leave, quarter by quarter, as rents follow excess demand.
- Schelling's sorting is an exact special case.
- The planner's permitting decision is an MOMDP: stock seen, demand hidden.

**Calibration.** CLUE's dwellings, Homes Victoria's rents of new lettings, the 2021 Census at SA1, and the Development Activity Monitor.

25 checks.

### 2.4 `occupancy`: Owners and Occupiers

**The model.**
- Buildings are agents on nomos's Gillespie engine, following physis's nested chain of uses (decision 0020).
- They are modulated by their neighbours' uses (affinity's co-location quotients) and by demand.
- The record required a vacancy channel the design did not name: falling vacant, letting, and returning to the family last held. With that channel and the neighbours shut off, the model is the nested chain exactly.
- It is lumpable onto occupied, unused and under works.
- It is calibrated by indirect inference.

**The planner's decision.** A weakly coupled POMDP over buildings under a budget, with tyche's Lagrangian bound.

17 checks.

### 2.5 The Case Studies

Each module has a record in `research/experiments/`, with several tasks reported in space, time and frequency:

| Record | Tasks |
|---|---|
| *Walkers of the Hoddle Grid* | a mean week (3.16 million trips; R² 0.720 against the counts, 0.998 against the flow); Little Bourke Street closed; a festival's crowd; the regime decided online |
| *Households of the Centre* | 2002–2025 calibrated (rents to log error 0.105, the Census's sorting); supply shocks; Schelling's limit; the planner's permits |
| *Owners and Occupiers* | the uses of 1,535 buildings over 2002–2024; where and why vacancy persists; interventions under a budget (97–98 % of the Lagrangian bound against 74–81 % for simple rules); what the neighbours add |

### 2.6 Their Visualisation

*Human flow* (status draft) draws every walker on the workbench's 3-D map, moving with the clock, with each street's density and queues.
- **Counters.** Their simulated readings stand beside the observed counts.
- **Tools.**
  - A street closed, one-way or open for chosen hours, run against the same hours open with the same draws.
  - An event's crowd from the City's permits.
  - POMCP's advice on the street's regime.
- **Engine.** C++ on `urban::movement`, built by the workbench's recipe.
- **Gallery.** 39 pictures in the gallery's §14: 29 figures from the records, and 10 screenshots of the workbench.

### 2.7 What Is Simulated Is Marked

Every simulated history, crowd and queue is written with provenance "simulated", and drawn apart from what was observed (`events/record.hpp`).

## 3. Alternatives Considered

- **Agents only inside the case studies.** Rejected: the models name a city's walkers, households and buildings, so they are the domain's, and the workbench and later studies reuse them.
- **A social-force model of walkers.** Deferred: the design keeps social force as a later, finer scale (§9). The network walkers are calibrated to what the data hold, the counters.
- **One plug-in for all three models.** Deferred: *Human flow* draws the walkers first. Households and buildings may enter as layers once their records are drawn in Event Spacetime.

## 4. Consequences

- The urban-science domain declares nomos (`module.json`), and its suite builds the three modules.
- The data library gains eight manifests (CLUE dwellings, blocks and floor space; building permits; Vicmap Planning's zones; Homes Victoria's rents; the 2021 Census at SA1 and its boundaries) and two products (`housing`, `occupancy`). All are CC BY 4.0, and their caches are not committed.
- **What the records found against their own models**, each stated in its record:
  - **Walkers.**
    - The logit over walks gives one daily tide for the whole district.
    - The model calls a closure a "saving": closing a street removes the wandering walks through its lanes.
    - The street's regime decision is degenerate, since closure wins whatever the crowd's size.
  - **Households.** The planner's plan turns on its assumed cost of a dwelling.
  - **Occupancy.** Vacancy's jump after 2020 and its clustering are not reproduced.
- **What the work asks of other units** is recorded in [development/ideas/2026-11-05-what-the-search-case-studies-ask-of-the-libraries.md](../ideas/2026-11-05-what-the-search-case-studies-ask-of-the-libraries.md) §4:
  - the logit over walks belongs in a foundation;
  - POMCP needs beliefs handed in and common random numbers;
  - a district builder belongs in `movement`;
  - the one-way regime strands walkers on one street;
  - the workbench's staleness check stalls requests.
