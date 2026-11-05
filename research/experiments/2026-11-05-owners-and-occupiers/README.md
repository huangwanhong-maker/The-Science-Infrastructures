# Owners and Occupiers of the Hoddle Grid's Buildings

*Status: draft (2026-11-05). Research.* A case study of the urban-science domain's `occupancy` ([theory/occupancy.md](../../../domains/urban-science/theory/occupancy.md); [ADR 0071](../../../development/decisions/0071-agents-and-decisions-for-urban-planning.md), [ADR 0075](../../../development/decisions/0075-agents-in-the-city.md)): an agent-based model of the uses of the 1,535 buildings of Melbourne's Hoddle Grid, 2002–2024, calibrated to the City of Melbourne's censuses, and a planner's interventions on vacant space as a weakly coupled POMDP.

This is a study of methods on a real place, not a planning proposal. Its values, costs and the effect of an incentive are assumptions (§4), and nothing in it says what any building should hold or which owner should be offered anything.

## 1. The Question

The user asked for agent-based models "in urban science with several tasks", visualised. This record takes the third of ADR 0075's models, owners and occupiers, through four tasks:

1. **The uses.** Can agents whose baseline is the record's own memory (physis's nested chain, decision 0020) reproduce the uses of the Hoddle Grid's buildings at CLUE's censuses 2002–2024, once calibrated to them?
2. **Vacancy.** Where does it persist, and why, judged by simulated histories beside the observed one?
3. **Interventions.** Which buildings should a budget for vacant space go to? The POMDP's ranked policy is set against simple rules (every vacancy, the longest vacant, by neighbours), under the Lagrangian bound.
4. **The neighbours.** What does the neighbours' modulation add over the chain alone?

Each is reported in space (maps of the buildings), in time (the censuses) and, where it means something, in frequency (the spectrum of a building's vacancy).

## 2. Method

1. **The record** (`data/products/occupancy_record.py`, the data library's `occupancy-hoddle-grid` product). It holds:
   - CLUE's base properties in the Hoddle Grid's box, each with its point, block, planning zone and footprint;
   - each property's use at every census: the predominant space use of its smallest property id with a use, the rule of Event Spacetime's engine;
   - the vacant floor space of the blocks;
   - the occupied establishments outside the box by year, the demand index;
   - the building permits naming a change of use, matched to properties by address.
2. **The model** (`urban::occupancy`, model.hpp). Each building's use changes in continuous time by four rules:
   - succession within the nested chain's levels, at λ0 and λ1;
   - an occupier leaving with no successor, at ν e^(−βd d);
   - a vacancy let, at η e^(βd d), back into the family it last held with probability ρ;
   - an incentive, κ η e^(βd d) more for an enrolled lettable vacancy.

   Each new use is chosen by a logit pulled by the six nearest neighbours' affinity: ln CLQ of the families at the 2002 census, within 60 m, clamped to ±2. With ν = η = βn = 0 the model is the nested chain exactly (checked).
3. **Calibration** (calibration.hpp). Indirect inference on ten statistics: the nested chain's two rates and the record's event rates (falling vacant, let, falling into works, returning to the family, the neighbours' pull, two responses to demand, moving between occupied families). There are 16 replicates of common random numbers per value of J, and the minimum is found by Levenberg and Marquardt, then Nelder and Mead. Three specifications are fitted: the chain alone; the model without the neighbours (βn = 0); the full model.
4. **Ensembles.** Two hundred histories of each specification, each building entering at its first census in its observed use and observed at its censuses, as the record is.
5. **Planning** (planning.hpp). From the 2024 census, each building is a POMDP on its lumped chain (occupied, unused, under works), its neighbours and the demand held at 2024. The year's budget joins them: tyche's Lagrangian bound and its ranked policy, against the simple rules. Every policy is evaluated in two worlds, paired against doing nothing on common seeds:
   - the planner's decoupled model: 80 episodes, the policy for 90 years and then nothing, valued exactly;
   - the agent-based model: 160 episodes over 2025–2034.

## 3. Code Version

The infrastructure at `c4adfdf` with this study's and the module's working files (not committed); calx `90aa710`, morphe `14b3fb8`, physis `e2b98d8`, tyche `b4034b8`, nomos `5382beb`.

```
cmake -S research/experiments/2026-11-05-owners-and-occupiers -B ~/.cache/science-build/exp-owners-occupiers -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-owners-occupiers -j4
~/.cache/science-build/exp-owners-occupiers/owners_and_occupiers data/cache/products/occupancy-hoddle-grid.json \
    data/cache/products/city-events-melbourne.json research/experiments/2026-11-05-owners-and-occupiers/results
python research/experiments/2026-11-05-owners-and-occupiers/figures.py
```

The run takes 4.7 minutes on four threads in WSL.

## 4. Parameters

| Parameter | Value | Source |
|---|---|---|
| Neighbours | the six nearest buildings | assumption |
| Affinity | ln CLQ of families, 2002 census, 60 m, clamped to ±2 | assumption; `urban::affinity` |
| Demand d(y) | ln of the occupied establishments outside the box, standardised over 2002–2024 | CLUE; the choice of index an assumption |
| Replicates | 16 a value of J; Σ from 100 | declared |
| Survey | a census surveys a building with chance 0.5 a year | assumption, from CLUE's "at least once every two years" |
| Costs | nothing 0, inspect 0.1, incentive 1 (incentive-years) | assumption |
| κ | 1: an incentive doubles an enrolled vacancy's letting | assumption |
| Discount | 0.95 a year | assumption |
| Budgets | 5, 10, 20 and 40 incentive-years a year | declared |

## 5. Data

| Dataset | Licence | Retrieved |
|---|---|---|
| CLUE building information (`melbourne-clue-buildings`, the Hoddle Grid's box) | CC BY 4.0 | 2026-11-04 |
| CLUE business establishments (`melbourne-clue-establishments`) | CC BY 4.0 | 2026-11-04 |
| CLUE floor space by use by block (`melbourne-clue-floor-space`, new) | CC BY 4.0 | 2026-11-05 |
| CLUE blocks (`melbourne-clue-blocks`, new, shared with the housing study) | CC BY 4.0 | 2026-11-05 |
| Building permits naming a change of use (`melbourne-building-permits`, new) | CC BY 4.0 | 2026-11-05 |
| Vicmap Planning zones (`vicmap-planning-zones`, new) | CC BY 4.0 (State of Victoria) | 2026-11-05 |
| The city model's buildings (`osm-city`) | ODbL 1.0, © OpenStreetMap contributors | 2026-10-31 |

The record holds 1,535 properties and 33,083 census rows; 1,283 footprints are joined (1,455 properties inside one, 67 by the nearest within 15 m, 13 none). Zones: CCZ1 1,091 properties, CCZ2 368, MUZ 52, others 24. The Capital City Zone admits almost every use, so zoning does not constrain the model within the grid. Of the 4,129 permit documents, 1,709 match a property in the box: 1,384 by an exact range, 300 within a range, 25 overlapping one. These make 699 matters.

## 6. Results

### 6.1 The Uses (Task 1)

**Memory.** The nested chain fitted to the record's 31,548 transitions half keeps a building's family for 7.68 ± 0.16 years and its use for 2.94 ± 0.14. On Event Spacetime's district alone (573 properties) the same fit gives 7.26 and 2.90, the plug-in's numbers. One level alone is rejected (twice the log-likelihood ratio 674).

**Calibration.**

| | J | Hansen (df) |
|---|---|---|
| the chain alone | 11,808.6 | 4,290.9 (8) |
| without the neighbours | 126.2 | 251.2 (4) |
| the full model | 8.35 | 62.1 (3) |

The full model's estimates:

| λ0 | λ1 | ν | η | ρ | βn | βd |
|---|---|---|---|---|---|---|
| 0.028 ± 0.001 | 0.153 ± 0.013 | 0.033 ± 0.001 | 0.230 ± 0.010 | 0.51 ± 0.03 | 2.08 ± 0.18 | −0.058 ± 0.028 |

All rates are a year.

**What the record demanded** (theory/occupancy.md §4.3):
- **Vacancy as transit.** The chain alone gives 0.011 falling vacant and 0.079 let in a year; the record has 0.033 and 0.220.
- **Memory through vacancy.** 54% of vacancies end in the family before them, against 16% under the chain alone.
- **The neighbours.** The new use's pull towards its neighbours is 0.070 in the record, 0.024 without the neighbours' term and 0.072 with it.

**The full model fits each statistic within 2.8 simulated standard deviations, but Hansen's test rejects it.** Most of the misfit is in falling into works (0.0074 modelled against 0.0060).

**The uses in time** (agents-23):
- **Shares.** Office fell from 26% of the buildings in 2002 to 18% in 2024 and residential rose from 15% to 20%; the model follows both trends but flattens them.
- **The jump.** Vacancy (unused or under works) held near 10% from 2006 to 2019, then jumped to 17% in 2020 and 21% in 2022–2023. Every model drifts smoothly between the two levels and holds neither: the demand index carries none of the jump (βd ≈ 0), and the pooled rates average the two regimes.

**A simulated history** (agents-21). Event Spacetime's district shows 998 changes of use between consecutive censuses in the record and 970 in one simulated history; 33% and 30% of its buildings never change. The record concentrates its changes in fewer buildings and in certain census years, while the model, with the same rates for every building, scatters them.

### 6.2 Vacancy: Where and Why (Task 2)

- **Where** (agents-24). 56% of the buildings were vacant at some census, and 329 at five or more. Observed and expected years vacant correlate at 0.44 building by building: the model's expectation, which carries each building's first use and neighbours, misses most of where vacancy lasts.
- **How long** (agents-25). A vacancy is let within a year with probability 0.06 after its first year and 0.48 after its second, then 0.16, 0.31, 0.20. Every model gives about 0.21 at every age. The saw-tooth is CLUE's survey cadence: a property is surveyed at least every two years and its use carried forward between surveys. It is a property of the record that the model does not represent.
- **Why: neighbours** (agents-26). The chance of falling vacant rises with the share of the six nearest buildings vacant: 0.029 with none vacant, 0.035 up to a fifth, 0.040 up to two fifths and 0.059 beyond. Every model is flat at about 0.032. The chance of being let does not depend on the neighbours (0.21 to 0.23). The record's vacancy stays clustered, a vacancy's neighbours 0.02 to 0.07 more vacant than an occupied building's; in every model the clustering of 2002 dissolves within four years. Vacancy spreads through falling vacant near vacancy, which the neighbours' term, acting on the owner's choice of use, cannot reproduce.
- **Block floor space** (agents-25). The vacant floor space of the grid's blocks (a lower bound) shows the same jump: 5.0% in 2019, 7.8% in 2020, 9.0% in 2024.
- **Permits.** Of the 559 permit matters with censuses either side, 28.8% are followed by a change of use within the two years, against 15.5% of all two-year pairs; 15.0% were vacant before, against 11.6%. The permits mark changes of use, independently of the censuses.
- **Frequency** (agents-27). Over the 1,308 buildings seen every year, the record's pooled periodogram of vacancy falls more steeply than the models'. It has less power at the two-year period (0.0115 against 0.0184), as carried-forward censuses would smooth it, and less at the longest period (0.32 against 0.41). The chain's closed form puts the family's indicator in one mode (memory 0.090 a year) and the unused use's in two (0.090 and 0.236).

### 6.3 Interventions under a Budget (Task 3)

In the planner's model (agents-28, left), the ceiling is the Lagrangian bound less doing nothing; the other columns are discounted unused building-years averted, paired, ± standard error.

| budget | ceiling | the POMDP | every vacancy | by neighbours | the longest vacant |
|---|---|---|---|---|---|
| 5 | 49.4 | 48.5 ± 1.0 (98%) | 39.0 ± 1.0 | 38.8 ± 0.9 | 29.6 ± 0.7 |
| 10 | 98.8 | 96.2 ± 1.5 (97%) | 79.4 ± 1.4 | 72.8 ± 1.2 | 60.2 ± 1.1 |
| 20 | 197.5 | 193.2 ± 2.3 (98%) | 154.2 ± 1.7 | 147.1 ± 1.7 | 118.7 ± 1.4 |
| 40 | 395.0 | 382.2 ± 3.5 (97%) | 302.0 ± 2.2 | 290.7 ± 2.3 | 241.1 ± 2.0 |

The ranked policy reaches 97–98% of the bound, so the best plan's value is known to within 2–3%. The simple rules reach 74–81% (every vacancy, by neighbours) and about 61% (the longest vacant).

**Why the longest vacant does worst.** In a Markov chain a long vacancy is no more responsive than a new one, and the longest-vacant buildings were often seen longest ago, so they are the least likely to be still unused.

In the agent-based model, 2025–2034 (agents-28, right), the table gives building-years averted, paired against doing nothing, and the POMDP's margin over every vacancy, paired:

| budget | the POMDP | every vacancy | by neighbours | the longest vacant | the POMDP − every vacancy |
|---|---|---|---|---|---|
| 5 | 25.1 ± 7.2 | 13.3 ± 7.1 | 16.1 ± 7.6 | 3.6 ± 7.4 | 11.8 ± 6.5 |
| 10 | 42.3 ± 7.5 | 23.0 ± 7.8 | 27.3 ± 7.3 | 18.1 ± 7.6 | 19.4 ± 8.6 |
| 20 | 71.1 ± 7.1 | 59.6 ± 7.8 | 55.4 ± 7.3 | 44.2 ± 6.7 | 11.5 ± 8.0 |
| 40 | 144.9 ± 7.6 | 115.8 ± 7.1 | 117.3 ± 7.4 | 83.2 ± 7.6 | 29.1 ± 7.3 |

The ranking holds in the agent-based model, where the neighbours' feedback the planner ignored is kept. The POMDP's margin over every vacancy is 1.4 to 4.0 standard errors.

**Where the margin comes from** (agents-29). In 2025 every candidate was seen unused in 2024 (belief 0.75), so the ranked policy picks among them by their own models. Its margin comes later, as the census skips buildings and the policy weighs what it believes against what it last saw.

### 6.4 What the Neighbours Add (Task 4)

- **The fit.** Without the neighbours J is 126.2, with them 8.35; the distance metric is 129.4 on one degree of freedom. The neighbours' term is what makes the model reproduce the new uses' pull towards their neighbours.
- **What it does not reproduce.** It does not reproduce vacancy's contagion or its clustering (§6.2): both models are flat in falling vacant by the neighbours' vacancy, and both lose the clustering.
- **Planning.** In the planner's decoupled model the two specifications give the same values (POMDP 96.2 and 96.1 at budget 10), since the neighbours are held. In the agent-based model the POMDP averts more under the full model at budgets 5 and 10 (25.1 against 14.2; 42.3 against 28.4), and about the same or less at 20 and 40 (71.1 against 72.2; 144.9 against 162.8). Targeting neighbourhoods ("by neighbours") gains nothing over every vacancy in either: from −5.6 to +5.7 building-years, all within one standard error. That is consistent with a model in which vacancy does not spread.

The next modulation the record asks for is the neighbours' vacancy in the rate of vacating, not only in the choice of use.

## 7. Figures

`assets/gallery/agents-20-the-uses-of-the-grid.jpg` to `agents-29-what-the-neighbours-add.jpg`, drawn by `figures.py` from `results/` and the occupancy product's footprints; captions in [gallery.md](gallery.md).

## 8. Caveats

- **Values and effects.** The values, costs and an incentive's effect (κ) are assumptions; no rent, lease or owner is observed.
- **Rejected by Hansen's test.** The full model fits each statistic but misses a combination of them, mainly falling into works.
- **The survey cadence.** CLUE's biennial survey carries uses forward. The model reads every census as a fresh observation, so its one-year rates average over carried and surveyed years.
- **2020–2024.** The jump in vacancy from 2020 is not in the model; its pooled rates average two regimes, and the demand index chosen does not carry the change.
- **Building heterogeneity.** The model gives every building the same rates; the record's churn concentrates in fewer buildings.
- **β_d's estimator.** On simulated records β_d is estimated about 0.1 high on average (the module's check), a finite-sample bias of the demand statistics. Here β_d is near zero either way.
- **The planner's model.** It holds the neighbours at 2024, so its bound is for the decoupled problem; spillovers appear only in the agent-based evaluation.
- **The simulated record.** `results/simulated-record-district.json` is one history at the calibrated parameters, written in the city-events record's tables with every event of provenance *simulated*. Event Spacetime does not yet read a simulated record beside the observed one.
