# The Households of the Centre

*Status: draft (2026-11-05).* Research. A case study of urban-science's `housing` ([theory](../../../domains/urban-science/theory/housing.md), [ADR 0075](../../../development/decisions/0075-agents-in-the-city.md)) on the City of Melbourne. **This is a study of methods on a real place, not a housing policy proposal**: its preferences are fitted to aggregates, several of its parameters are assumptions, and its scenarios are what the model does, not forecasts.

## 1. The Question

How do households sort among the centre's dwellings by rent, type and place, and what do changes of supply do to rents and to who lives where? Four tasks:

1. **The households of the centre, calibrated**: an agent-based model of households among the City of Melbourne's 11,036 residential buildings, run quarter by quarter from 2002 to 2025, against Homes Victoria's rents of new lettings and the 2021 Census.
2. **A supply shock**: the Development Activity Monitor's pipeline of towers, or of student apartments, completed 2026–2030.
3. **Schelling's limit**: the model reduced exactly to Schelling's on a lattice; in the city, how far the Census's sorting is from the sorting that neighbours alone would make.
4. **The planner's permits under uncertain demand**: a MOMDP whose seen part is the stock and hidden part demand, solved by HSVI and set against simple policies, with the cycles of rents each makes.

The infrastructure had agent-based models and POMDPs in its foundations (ADR 0071) and no model of households. The study asks what the two can say together on a real city, and where the data stop them.

## 2. The Method

**The market** (`urban/housing/market.hpp`). Households of class c choose a segment (a CLUE dwelling type in a CLUE small area) by a logit in V = δ + γ − α_c ln(r / r̄) − θ d, r̄ the metropolis's rent (the outside option), α_c = α (ȳ / y_c)^π.

**What the rents reveal** (`calibration.hpp`). A segment's log rent regressed on its log stock with quarter and segment effects, on Homes Victoria's flats (and houses where enough are let) in the six suburb groups inside the municipality, 2003Q1–2025Q3. The slope is −1/α. The quarter effects, less the log of the metropolis's rent and times α, are the log of the demand index: the households who consider the centre.

**The agents** (`households.hpp`, on nomos's random sequential schedule and keyed draws).
- Each quarter, residents leave the centre at an exit hazard, or consider a move.
- Newcomers arrive in proportion to the demand index.
- Searchers choose among the dwellings on offer (a block's dwellings of a type are one exact option), their own and the outside option.
- Each segment's rent moves by κ ln((D + 1)/(V + 1)), D the searchers' notional demand for what was posted and V what was posted.
- A dwelling vacated is let from the next quarter, so vacancy is turnover's.

**Calibration.**
- α from the regression.
- Classes, tenure, size, exits and moves from the 2021 Census.
- The segments' values δ and the owners' taste for houses by a burn-in at 2002's rents: Berry's contraction on the flow of lettings, and the Census's owners among houses.
- κ by the rents' log error over a grid.
- π by the Census's dissimilarity of low from high incomes across SA1s.

**The supply shock.** The 2025Q3 state of each seed carried to 2030Q4 under four scenarios, with common seeds: demand and the metropolis's rent grow at their 2010–2019 means, and the pipeline's dwellings are phased in at their buildings.

**Schelling's limit.**
- On a 40 × 40 torus: the households' model in its limit (no rent, type or distance, a step utility at f, the best option) against `nomos::agents::Schelling`, tick by tick.
- In the city of 2021Q3: the same limit by block with three income classes, and the calibrated model with a taste η for like neighbours.

**The planner** (`planning.hpp`, tyche's MOMDP and HSVI).
- The state: the stock in units of 10 % of 2024's apartments (seen), demand in the same units (hidden).
- An action: 0, 1 or 2 units permitted, completed the next year.
- An observation: the band of the clearing rent.
- The reward: the market's welfare (households' surplus and rents) less the cost of the units permitted.
- The demand's yearly steps are counted from the revealed demand index.

**Frequencies.** Thomson's multitaper spectra (calx `signal`) of the quarterly growth of the centre's rents, reported and simulated. Periodograms of the clearing rent under each of the planner's policies.

## 3. The Code Version

The infrastructure at c4adfdf, with the uncommitted `housing` module, the data library's new manifests and product, and this record. The foundations: calx 90aa710, morphe 14b3fb8, physis e2b98d8, tyche b4034b8, nomos 5382beb.

```
python data/acquisition/fetch.py melbourne-clue-dwellings melbourne-clue-blocks homes-victoria-rents abs-census-2021-sa1-vic abs-asgs-2021-sa1
python data/acquisition/fetch.py melbourne-development-activity --set south=-37.86 --set west=144.88 --set north=-37.77 --set east=145.0
python data/products/housing_records.py
cmake -S research/experiments/2026-11-05-households-of-the-centre -B ~/.cache/science-build/exp-households -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-households -j2
~/.cache/science-build/exp-households/households data/cache/products/housing-melbourne.json research/experiments/2026-11-05-households-of-the-centre/results all
.venv/Scripts/python.exe research/experiments/2026-11-05-households-of-the-centre/figures.py
```

`households ... shock,planner reuse` reruns chosen tasks with the κ and π an earlier run chose (results/trials.json). The whole run took 21 minutes on eight threads, the η sweep most of it.

## 4. The Parameters

| Parameter | Value | Source |
|---|---|---|
| α, the semi-elasticity of choice in rent | 7.14 ± 1.49 (slope −0.140, clustered s.e. 0.029; HC1 0.010) | fitted: 728 segment-quarters, 8 segments |
| π, how rent's weight falls with income | 0.35 (grid 0, 0.25, 0.35, 0.5, 0.75) | fitted: the Census's dissimilarity 0.401 |
| κ, the rents' adjustment | 1.6 a quarter, capped at ±0.15 in log (grid 0.1–1.6; the log error falls from 0.133 to 0.109 across it, so κ lies at its edge) | fitted: rents of new lettings |
| income classes | low < A$650 a week (18.7 %), middle A$650–1,999 (39.6 %), high ≥ A$2,000 (41.7 %); class medians A$316, 1,307, 3,153 | Census G33, the municipality's 353 SA1s; the open top band taken as A$4,000–6,000 (assumption) |
| tenure | owners 30.6 % (52.6 % of houses, 27.0 % of flats); newcomers' owners 9.9 %, so that exits keep the residents' share | Census G37; the newcomers' share derived |
| size | one person 43.1 % | Census G35 |
| exits and moves | persons at a different address a year before 42.7 %: from another SA2, overseas or not stated 35.7 % (exits), within the SA2 6.9 % (moves); a quarter, renters 0.136 and 0.023, owners a quarter of those | Census G44; the owners' ratio an assumption |
| newcomers | 46,966 a quarter in 2002, twenty times the exits, then in proportion to the demand index | assumption: one in twenty of those who look settle, as P₀ ≈ 1 requires |
| θ, distance to the Town Hall | 0.3 a km | assumption (not identified by the aggregates) |
| initial occupancy, 2002 | 96 % | assumption |
| student apartments | open to one-person households; their 2002 rent the group's one-bedroom flats | assumption: the report holds no student apartments |
| the planner | a unit 9,562 apartments; A$410,000 a dwelling (A$3.92 billion a unit); discount 0.93 a year; rent bands within ±1 % of 2025Q3's A$655; a reading's error 0.15 each side; one in twenty of the market's households in the centre | assumptions; the cost chosen so that a unit pays at some demands and not at others (at A$380,000 the plan permits at every demand, at A$430,000 almost never) |
| seeds | 11–14 (the history and the shocks), 1 (the calibration trials), 2026 (the lattices) | — |

## 5. The Data

| Dataset | Manifest | Licence | Retrieved |
|---|---|---|---|
| CLUE residential dwellings by building and type, 2002–2024 (219,680 rows; 11,036 buildings; 2,869 rows without a point, placed by their property's other years or their block) | `melbourne-clue-dwellings` | CC BY 4.0 | 2026-11-05 |
| CLUE blocks (603) | `melbourne-clue-blocks` | CC BY 4.0 | 2026-11-05 |
| Homes Victoria, moving annual median rents of new lettings by suburb and type, 2000Q1–2025Q3, the September 2025 workbook | `homes-victoria-rents` | CC BY 4.0 | 2026-11-05, from the Wayback Machine's copy (the publisher's server refuses this library's requests) |
| ABS 2021 Census, General Community Profile, Victoria's SA1s: G02, G33, G35, G36, G37, G40, G41, G44 | `abs-census-2021-sa1-vic` | CC BY 4.0 | 2026-11-05, members of the DataPack by range requests |
| ABS ASGS 2021 SA1 boundaries about the municipality | `abs-asgs-2021-sa1` | CC BY 4.0 | 2026-11-05 |
| The Development Activity Monitor, the municipality (1,451 developments; under construction 8,196 dwellings and 3,678 student apartments, approved 21,337 and 5,254) | `melbourne-development-activity` (a box set to the municipality) | CC BY 4.0 | 2026-11-05 |

The product `data/products/housing_records.py` joins them in the city model's local frame (`housing-melbourne.json`). ABS data used with permission from the Australian Bureau of Statistics.

**What the data cannot give.**
- No rent or move of a single dwelling.
- No vacancy between censuses.
- No rents of student apartments.
- No incomes crossed with tenure or type at SA1.
- One Census, taken on 10 August 2021 in a lockdown, when 26,131 of 96,878 private dwellings (27 %) were unoccupied.
- Preferences are fitted to aggregates, never observed.

## 6. The Results

### 6.1 The Households of the Centre, Calibrated

- **What rents reveal.** Across the eight segments, a segment whose stock grows 1 % against the others' lets for 0.140 % less (clustered s.e. 0.029), so α = 7.14 ± 1.49.
  - The revealed demand, the metropolis's rents taken out, was nearly flat from 2003 to 2019 (1.45 in 2019Q3, 2003Q1 = 1) while the stock tripled and the centre's rents fell against the metropolis's.
  - It fell to 0.32 in 2021Q3, when the borders were closed to international students, and stood at 0.68 in 2025Q3.
- **Rents.** The model, calibrated only at 2002's rents and driven by CLUE's stock and the revealed demand, follows the six groups' flats to 2025.
  - The centre's index has a log error of 0.105 (seeds 0.098–0.110), the groups' 0.139.
  - It overshoots the CBD, Southbank and Docklands in 2013–2019 by up to A$90 a week, and undershoots North and West Melbourne throughout.
  - It falls a year late and too far in the closure: A$425–459 in 2021Q3 against 377, A$329–343 in 2022Q3 against 447.
  - In 2025Q3 it lets at A$625–657 against 655.
  - κ sits at the edge of its grid, and the model's rents move in steps: the variance of their quarterly growth is 2.6–3.0 times the report's.
- **The 2021 Census against the model's 2021Q3** (four seeds).
  - The dissimilarity of low from high incomes across SA1s is 0.395 ± 0.018 against the Census's 0.401, π fitted to it.
  - Renters are 73.7 % against 69.4 %.
  - Households that moved in within the year are 45.0 %, against 42.7 % of persons at a different address.
  - The model keeps too few low incomes, 10.9 % against 18.7 %, and too many high, 53.5 % against 41.7 %. The Census's low incomes in the centre are largely students, whose rents are not paid from the incomes the Census counts. A rent coefficient by income cannot hold them at the centre's rents.
  - Across SA1s with more than 20 renters, the model's renters' median rent correlates with the Census's at 0.29, and its median is A$360 against A$401: the model's rents were in their trough of the closure in 2021Q3.

### 6.2 A Supply Shock

Mean of 2027Q1–2030Q4 against the baseline, common seeds (mean ± s.d. over four seeds):

| Scenario | Added | The centre's flats, A$ a week | Student apartments | Low-income households, 2030Q4 |
|---|---|---|---|---|
| baseline | — | 714.9 | 554.7 | 14,937 |
| towers under construction | 8,196 apartments | −11.9 ± 12.5 (−1.7 %) | −5.8 ± 9.1 | +1,089 |
| every approved tower | 29,533 apartments | −50.3 ± 11.5 (−7.0 %) | −30.9 ± 9.3 | +5,920 (+40 %) |
| student apartments | 8,932 student apartments | −19.6 ± 10.0 (−2.7 %) | −58.3 ± 3.9 (−10.5 %) | +2,065 |

- **Rents fall while the new dwellings fill, then return.** Every approved tower holds the centre's rents near A$650 to 2029, while the baseline's rise to A$747. They regain it by 2030Q4, as the flow of newcomers fills the stock.
- **The lasting effect is who lives there.** Low-income households are 40 % more numerous with every approved tower: they stay and enter where the baseline prices them out.
- **Who takes the new dwellings.** Of the 25,774 households in the new towers in 2030Q4, 14 % have low, 36 % middle and 50 % high incomes. The student apartments take 1,996 low-income households of 8,014.
- **The student apartments' own rent** falls from A$512 to A$379 in the half year after 2027Q1, when the pipeline's student apartments begin to arrive, then recovers to A$574 by 2030Q3 (the baseline A$595): the large κ makes the model's rents ring.

### 6.3 Schelling's Limit

- **On the lattice the reduction is exact.** At f = 1/3 and 1/2, none of 40 ticks differs between the households' model and nomos's Schelling: the like-neighbour share runs from 0.497 to 0.720 and to 0.867 in both, draw for draw.
- **In the city, sorting by price alone gives the Census's sorting.** The calibrated model with no taste for neighbours (η = 0) has D = 0.406; the Census has 0.401.
- **The neighbourhood taste over-sorts quickly.** With η = 2, 4 and 8 the model's D is 0.457, 0.752 and 0.918, and its rents fit worse (log error 0.123, 0.131, 0.151 against 0.105).
- **Sorting by neighbours alone**, from 2021's households, prices, types and distance off, settles in three ticks: D 0.578 at f = 1/3 (content 0.795 to 0.944), 0.751 at f = 1/2 (0.426 to 0.866). The city sits at the sorting that prices give: its D is 69 % of the limit's at f = 1/3 and 53 % of it at f = 1/2.

### 6.4 The Planner's Permits

- **The market behind the plan.** One segment, the municipality's 95,625 apartments of 2024, and the Census's three classes at α and π as calibrated. Its values are set so that 2025Q3's stock clears at the report's A$655.
  - Over seven levels of stock (+0 to +57,375 apartments) and eleven of demand (−30 % to +70 %), the clearing rent runs from A$575 to A$714.
  - The year's welfare runs from A$3.41 to A$5.84 billion.
- **Demand's steps.** The revealed demand's 22 yearly changes, in units: −2 (9 %), −1 (27 %), 0 (32 %), +1 (9 %), +2 (23 %). The changes of 2020–2023 are as large as −15 and +7 units, clipped at ±2.
- **The plan.** With A$410,000 a dwelling and a discount of 0.93, a unit pays only where demand runs ahead of the stock.
  - Where demand is known, the plan permits one unit at the 2024 stock when demand stands 10 % above today's, and two units from 20 %.
  - At +19,125 apartments it permits from 30 %; at the largest stock, never.
- **The policies** (3,000 simulated episodes of 150 years, two standard errors):
  - HSVI's interval for the optimum from today's belief is [52,734.7, 52,739.6] A$ million after 16 trials.
  - HSVI's own policy earns 52,746.7 ± 62.8, QMDP's 52,745.9 ± 62.2 and the certainty equivalent's 52,746.4 ± 62.6. They are indistinguishable from each other and from the interval, and permit in 3–5 % of years.
  - Never permitting earns 52,434.0 ± 51.8: planning under uncertainty adds A$312 million, 0.6 %.
  - Following the rents (permitting when the expected band is high) earns 46,728.8 ± 198.4. It permits in a quarter of the years, whenever rents read high, and loses A$5.7 billion against never permitting.
  - One unit a year earns 16,247.5 ± 73.0.
- **What the hidden part is worth here.** Little: the rent bands within ±1 % reveal demand well enough that the certainty equivalent does as well as the plan. The planner's problem is the timing of irreversible supply, and the policy that reacts to rents is the one that loses.

### 6.5 Rents and Their Cycles

- **The report's rents.** The multitaper spectrum of the centre's quarterly rent growth (89 quarters, 2003Q3–2025Q3) peaks at 0.225 cycles a year (4.4 years), the closure and its rebound dominating. 92 % of its power lies below half a cycle a year.
- **The model's rents** peak at the same frequency in every seed, with 72–83 % below half a cycle a year and 2.6–3.0 times the power: the steps of a large κ.

- **The planner's rents.** The clearing rent, simulated a year at a time over 400 paths of 120 years and read after 20, has a red spectrum under every policy, rising to the longest period with no peak: no cobweb, since supply cannot be withdrawn and the lag is a year.
  - Planning lowers the rent's variance from 0.00144 (never permitting) to 0.00122 (HSVI's policy, 15 % less).
  - Following the rents lowers it to 0.00125.

## 7. The Figures

| Figure | Content |
|---|---|
| [agents-10](../../../assets/gallery/agents-10-the-centres-dwellings.jpg) | The centre's dwellings by building in 2024, the stock by type 2002–2024, and by small area |
| [agents-11](../../../assets/gallery/agents-11-what-rents-reveal.jpg) | The inverse demand: within-segment log rent against log stock; the revealed demand and the metropolis's rent |
| [agents-12](../../../assets/gallery/agents-12-rents-reported-and-simulated.jpg) | Six groups' rents of new lettings, reported and simulated, 2002–2025 |
| [agents-13](../../../assets/gallery/agents-13-the-census-and-the-model.jpg) | The low incomes' share by SA1, Census and model; small areas' income classes and renters |
| [agents-14](../../../assets/gallery/agents-14-a-supply-shock-in-time.jpg) | The supply scenarios in time: rents, student rents, low incomes resident |
| [agents-15](../../../assets/gallery/agents-15-a-supply-shock-in-space.jpg) | The supply scenarios in space: low incomes by block against the baseline; who moves into the new buildings |
| [agents-16](../../../assets/gallery/agents-16-schellings-limit-exactly.jpg) | Schelling's limit on a lattice: the two programs' trajectories and lattices |
| [agents-17](../../../assets/gallery/agents-17-how-far-from-schellings-limit.jpg) | The city's low incomes by block in the model and in Schelling's limit; D against η |
| [agents-18](../../../assets/gallery/agents-18-the-planners-permits.jpg) | The planner's plan over supply and demand, and the policies' values |
| [agents-19](../../../assets/gallery/agents-19-rents-and-their-cycles.jpg) | Spectra of rents: reported and simulated; under the planner's policies |

## 8. The Caveats

- **Demand is read from rents.** The demand index comes from the same rents the model is held to: the centre's aggregate path is fitted in part by construction. The tests are the groups' paths, the Census's cross-section and the scenarios' differences.
- **The identifying assumption** (stock predetermined) is a model. Developers build where demand grew, which would bias α.
- **Newcomers settle one in twenty**, an assumption that P₀ ≈ 1 requires. With one in two, the model could not absorb the closure's fall of demand without a collapse of rents.
- **κ lies at its grid's edge.** The rents ring after large shocks (the student apartments), and their quarterly growth is three times as variable as the report's.
- **The model cannot hold the Census's low incomes**, who are largely students.
- **Rents are one per segment.** Sitting tenants pay the segment's current rent; no lease is modelled.
- **The 2021 Census was taken in a lockdown.**
- **The planner's costs, discount, bands and reading error are assumptions**, and its welfare counts rents as surplus to their owners.
- **The scenarios are what the model does under its assumptions**, not forecasts or advice.
