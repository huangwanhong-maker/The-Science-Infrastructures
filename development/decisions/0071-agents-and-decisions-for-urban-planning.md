# 0071. Agents and Decisions for Urban Planning

- **Status:** Accepted, provisionally. Built: tyche's `decision` and `processes`, nomos's `agents`, calx's categorical draws and urban-science's `maintenance`. Housing, movement and occupancy are planned. Amended by [ADR 0075](0075-agents-in-the-city.md), which builds housing, movement and occupancy.
- **Date:** 2026-11-04
- **Decided by:**
  - **The request.** The user asked to "further implement urban science lib ... two important module for modern urban planning. One is ABM, another one is POMDP (many variant). Both of two can be reuse. So I think it would be seperate from the applications."
  - **The user's choices.** Asked, they chose:
    - agent-based models in nomos's planned `agents`, and POMDPs in a new tyche module `decision`;
    - all four urban problems offered (maintenance and heritage, housing and residential choice, pedestrians and streets, land use and vacancy): "All you think is needed. You may also search review papers.";
    - all four families of POMDP (exact and point-based, online Monte Carlo, factored with mixed observability, multi-agent and constrained).
  - **What Claude did.** It wrote the design ([development/ideas/2026-11-04-agents-and-decisions-for-urban-planning.md](../ideas/2026-11-04-agents-and-decisions-for-urban-planning.md)), after surveys of the literature and of Melbourne's open data. It then wrote the foundations' modules (tyche decisions 0009 and 0010, nomos decision 0004), the first urban model, and this record.

## 1. Context

Planning a city means deciding under uncertainty about things one sees only in part:
- the condition of its assets between inspections;
- the uses of its buildings between censuses;
- the crowd on streets without sensors;
- the demand for its dwellings.

It also means understanding populations whose individual decisions compose: households, owners, walkers. The infrastructure held neither POMDPs nor agent-based models. Both name no city, so by the rule that general methods go to the foundations (CLAUDE.md §2.1, the user's principle) they belong below the domain.

## 2. Decision

### 2.1 Placement

| Unit | Holds |
|---|---|
| tyche `decision` (tyche decision 0009) | finite POMDPs and their bounds; exact value iteration; PBVI, Perseus and HSVI; POMCP; MOMDPs; components joined by a budget with the Lagrangian bound, grouped by shared models and weighted; constrained POMDPs by duality, with finite-state controllers; Dec-POMDPs and JESP |
| tyche `processes` (tyche decision 0010) | continuous-time Markov chains observed at discrete times: the panel likelihood by uniformisation, its maximum and errors; sequential chains |
| nomos `agents` (nomos decision 0004) | schedules (synchronous on a pool, random sequential, Gillespie), draws keyed by agent and tick, event logs, ensembles, spaces, logit choice, Schelling's model, the ODD+D description |
| calx `random` | `categorical.hpp`: categorical draws by inversion and by Walker's alias table |
| urban-science `maintenance` | the deterioration of the city's surfaces fitted to their ratings, each surface a POMDP, the city planned under a budget |
| urban-science `housing`, `movement`, `occupancy` | planned (the design's §§8–10) |

### 2.2 The First Urban Model: the Maintenance of Surfaces

**The data.** The City of Melbourne's road surfaces in the Hoddle Grid: 9,630 rated, each with its latest condition rating, its evaluation date and its installation date. This is a new manifest, `melbourne-road-surfaces`, and a new product, `road_condition.py`.

- **The deterioration.** A sequential chain through the five classes, fitted by type to the ratings at their ages: current-status data, on tyche's panel likelihood. Installation dates that are placeholders (1 January or 25 December) are excluded, and counted.
- **Each surface.** A POMDP on its fitted chain, with four actions (nothing, inspect, maintain, renew). Its belief comes from its last rating, carried forward by the years since.
- **The city.** The surfaces' weighted sum, joined by the year's budget. The Lagrangian bound is an upper bound on any plan within the budget, and the ranked plan keeps it.
- **The assumptions.** Service values, costs and an inspection's noise are declared, not fitted.

## 3. Alternatives Considered

- **Both methods in urban-science** (the user's first thought: kept apart from the applications). Taken further: kept apart from the domain too, since robotics, missions, economics and history need them.
- **POMDPs in nomos or in physis's `control`.** Rejected (tyche decision 0009, §3).
- **Borrowed frameworks** (APPL/SARSOP, AI-Toolbox and the MADP toolbox are GPL; Mesa and NetLogo are Python and Java). Rejected for the engines. The GPL tools may serve as external oracles in fixture scripts.
- **A deterioration fitted only by the design life.** Rejected. The data give each surface's rating at its age, and the design life is the City's assumption, not an observation.

## 4. Consequences

- **Checks.**
  - tyche `decision`: 57, among them the published tiger value at discount 0.75 and Dec-Tiger's optima at horizons 2 and 3.
  - tyche `processes`: 5.
  - nomos `agents`: 14, among them Gillespie against the exact master equation and Kurtz's mean field.
  - calx `random`: 7.
  - urban-science `maintenance`: on the real surfaces.
- **The order of the domains** is unchanged. urban-science will declare nomos in its `module.json` when its first agent model is written.
- **The data library** gains a manifest and a product. Its survey also found that the pedestrian manifests' licence should now read CC BY 4.0. That is noted for the data library, not changed here.
- **Next.**
  - `occupancy`: the district's uses as a POMDP over the record, the fitted nested chain its transition, with agents calibrated by indirect inference.
  - `movement`: walkers on the pedestrian network, calibrated to the counters, with street regimes as decisions.
  - `housing`: households on CLUE's dwellings, calibrated to Homes Victoria's rents.
  - The heritage overlay as a constraint of maintenance.
