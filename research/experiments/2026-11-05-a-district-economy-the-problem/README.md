# A District's Economy: the Problem

*Status: draft (2026-11-05). Research.* The shared problem of three case studies of planning a district's economy, on the search foundation heuresis ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)) and the urban-science, architecture and economics of a real place: the ground floors of a few blocks of Melbourne's Hoddle Grid.

This is a study of methods on a real place, not a planning proposal. Its economic parameters are assumptions (§4.3), its model of a shop is a sketch, and nothing in it says what any frontage should hold.

## 1. The Question and the Three Studies

The user asked for three plans of a district's economy, each considering space, time and frequency, the people's flow, and the arrangement of the buildings and their relations:

1. **without a structured model**: the flows come from the data alone ([study 1](../2026-11-05-a-district-economy-without-a-model/README.md));
2. **with a fluid model**: the flows are a fluid on the network, conserved, with their circulation, and the search uses the model's structure (study 2);
3. **with a field theory**: the flows are a field whose law is the stationary point of an action S (study 3).

This record defines the one problem the three share, so that they differ only in how much structure their model of the city has. A study supplies one thing: the people passing each frontage at each hour of the mean week, P_f(h), for any plan. Everything else is here and in [`economy.hpp`](economy.hpp): the district, the decisions and their encoding, the value and the other measures, the calibration of the people from the counts, and what every study reports.

*Rationale.* A comparison of models of flow is fair only if the models are asked the same question. A study that changed the value, the decisions or the district along with its model would confound the two. So the economics is fixed once, with its assumptions in the open, and only the flows are left to the models.

## 2. The District

### 2.1 The Area

The area is that of the case study [Designing a District](../2026-11-05-designing-a-district/README.md): the blocks between Swanston and Russell Streets, Flinders Lane and Bourke Street.

- **The walking network.** It is Landscape Dynamics' district of the Urban Science Workbench, prepared by the designing study's `prepare.py`, whose output (`district-inputs.json`) this problem reads rather than makes again. It has 356 junctions and 399 segments, and the urban-science domain's complex has 356 0-cells, 399 1-cells and 37 2-cells (`urban::landscape::planarised`, `district`), with no segment dropped.
- **The parcels.** The 61 buildings of the city model whose footprint's centre lies in the area, in the designing study's order along the area's long axis (south to north).
- **The context.** The 524 other outlines of the box stand as built.

### 2.2 The Frontages

Each parcel has one ground-floor frontage, the unit of decision.

- **From the domain.** The parcels and their context are massed as built by architecture's `massing`, and the urban-science domain's city network (`urban::events::city_network`) reads from the massing which outer walls face which street segment (within 25 m, facing it).
- **The main facade.** A frontage is the parcel's facade on the street or mall it fronts most, if it fronts one by at least 3 m; else the laneway it fronts most. A building's back on a service lane is not its shopfront. Its length is the facade's area over the building's height, and its anchor the middle of its wall pieces' base lines.
- **The result.** 60 of the 61 parcels front a street; one (on Russell Street, enclosed by its neighbours) fronts none and stands on its nearest segment with no length. The frontages face 29 segments: Collins Street (14 frontages), Swanston Street (10), Little Collins Street (9), Russell Street (7), Bourke Street (7), Russell Place (6), Flinders Lane (4), Royal Lane (3) and Rainbow Alley (1). Their main facades are 1,277 m long.
- **The trading floor.** A frontage's trading floor is the lesser of its footprint and its length times 20 m (an assumption): 23,252 m² in all.

### 2.3 The Relations of the Buildings

The arrangement of the buildings is held as three relations, which the value reads (§4) and every study draws (§6).

| Relation | What it is | Count |
|---|---|---|
| party walls | pairs of parcels that share a wall: 2-cells of the massing with two cofaces | 54 pairs; none with the context, since the area is whole blocks |
| one segment | frontages that face the same segment of the network: they compete (§4.1) | 29 segments fronted |
| reach | pairs of frontages within 3ξ = 127.6 m of each other along the network: they catalyse each other (§4.2) | 335 pairs with no link made |
| blocks | the place each parcel stands in (the domain's 2-cells, then the faces it leaves as holes): the mix is read within them (§4.4) | the designing study's places |

The walking distances between frontages are computed on the domain's district with the plan's links made, each frontage on its nearest segment of the streets that were there before the links (a link is a way past them, not their frontage). The distances are cached by the set of links made.

## 3. The Decisions

### 3.1 The Kinds

A frontage holds one of six kinds. CLUE's establishments are mapped to them by their ANZSIC 2006 class (`prepare.py`, `kind_of`).

| Kind | ANZSIC |
|---|---|
| vacant | CLUE's 0000, Vacant Space; or no establishment |
| retail | division G, Retail Trade (39–43) |
| food and drink | subdivision 45, Food and Beverage Services |
| services | 62 banking, 641–642 auxiliary finance, 672 real estate agents, 722 travel agents, 85 medical and allied health, 94–95 repair and personal care, 51 postal |
| office lobby | every other establishment: offices, education, public administration, accommodation; a ground floor that is a building's entrance |
| culture | division R, Arts and Recreation (89–92), and 551 motion picture exhibition |

### 3.2 The Plan as Built

The ground floor's kind as built is read from CLUE's 2024 census.

- **The rule (an assumption).** A frontage takes the street-facing kind (food and drink, retail, services, culture) with most establishments in its parcel, ties broken in that order. Failing one, it is an office lobby if the parcel has any other establishment, and vacant otherwise. CLUE does not say which floor an establishment is on.
- **The parcel's establishments.** They are those of the CLUE base property the data library's city-events product joins to the building. Where the join gives none, they are the establishments whose CLUE position lies inside the footprint. Where neither gives any, the frontage is vacant.
- **The result.** 18 frontages of food and drink, 14 of retail, 13 of services, 4 office lobbies and 12 vacant. Of the vacant, 11 are parcels CLUE does not record (among them the City Square's entrance, the Council Chamber, the Visitor Services kiosk and several small buildings), and one is recorded as vacant space. Being unrecorded is read as vacant, which is an assumption.

### 3.3 The Hours

A frontage's hours are two pairs of integers: on a weekday, a start s ∈ {0, …, 23} and a length ℓ ∈ {0, …, 24}; and the same at the weekend. It is open over [s, min(24, s + ℓ)) civil hours and closed when ℓ = 0. Every weekday is alike, and both weekend days are alike. An end hour past midnight is not offered.

The hours as built are each kind's median opening and closing on a Wednesday and on a Saturday, from the opening hours declared in OpenStreetMap for the holdings of five Melbourne streets (the data library's `street-holdings` product: Swanston Street, the Bourke Street Mall, Lygon, Errol and Macaulay Streets).

| Kind | Weekday | Weekend | Holdings that declare hours |
|---|---|---|---|
| retail | 9–18 | 10–18 | 44 |
| food and drink | 11–22 | 11–22 | 60 |
| services | 9–18 | 9–15 | 22 |
| office lobby | 9–17 | 9–13 | 6 |
| culture | 10–22 | 10–22 | 2: an assumption |

### 3.4 The Interventions

The interventions are the designing study's eleven candidate links, each an arcade at ground level, made or not. Pedestrianising an existing service lane is not among them.

| Link | Kind | Length (m) | Weekly cost (A$) |
|---|---|---|---|
| Ramsden Place to Collins Street | lane carried through | 91.1 | 2,272 |
| a lane to Regent Place | lane carried through | 28.9 | 721 |
| Rainbow Alley to Bourke Street | lane carried through | 66.5 | 1,657 |
| Athenaeum Place to Collins Street | lane carried through | 65.6 | 1,637 |
| Baptist Place to Collins Street | lane carried through | 36.1 | 901 |
| a lane to Athenaeum Place | lane carried through | 19.5 | 487 |
| Portland Lane to Russell Place | lane carried through | 22.5 | 562 |
| Donaldson Lane to Russell Place | lane carried through | 21.3 | 531 |
| Rainbow Alley to Little Collins Street | passage | 30.8 | 768 |
| Russell Street to Portland Lane | passage | 36.3 | 906 |
| Russell Place to Royal Lane | passage | 39.1 | 975 |

### 3.5 The Encoding

- **The mixed encoding.** A plan has 61 categorical kinds, 244 integer hours and 11 bits (`economy::Plan`). Every value of it is a valid plan.
- **The canonical relaxation.** For continuous searches, u ∈ [0, 1]^621 (`Problem::decode`). For each frontage it has six kind scores, decoded by argmax (the first on a tie), and four hour numbers: a start ⌊24u⌋ capped at 23 and a length ⌊25u⌋ capped at 24. For each link it has one number, made when u ≥ ½. `encode` gives the centre of a plan's cell, so that decode(encode(p)) = p.
- **The kernels.** `economy::Kernels` holds the mixed encoding's kernels, which a study may use:
  - **sampling:** uniform;
  - **mutation:** a kind redrawn at one frontage in 61 (heuresis `reset_mutation`), an hour crept by at most 3 at one gene in 61 (`creep_mutation`), and a link flipped at one in 11 (`standard_bit_mutation`);
  - **crossover:** a frontage's kind and hours travel together, and two cuts along the area exchange a stretch of it, since neighbours catalyse and compete with one another; links cross uniformly.

## 4. The Value

### 4.1 The Definition

The value is the week's gross operating surplus of the district's frontages, in A$ per week of the mean week. It is the margin on the spending they capture, less the costs of opening them and of the interventions:

V = Σ_f Σ_{h open} P_f(h) · p_f(h) · c_f(h) · s_k · m_k − Σ_f (staffing_f + outgoings_f + refit_f) − Σ_{links made} cost

- **The people passing**, P_f(h): people an hour passing frontage f in hour h of the mean week (Monday 00:00 civil first). This is the study's model of flow.
- **The chance to stop**, p_f(h) = π_k · φ_k(band of h) · clamp(1 + γ Σ_g w_fg a(k_g → k_f), ½, 2): the kind's propensity π_k by day, its multiplier φ_k for the day (07–18), the evening (18–24) or the night (00–07), and the catalysis of the open frontages g ≠ f within 3ξ (§4.2).
- **Competition**, c_f(h) = [1 − Π_{g∈G} (1 − p_g)] / Σ_{g∈G} p_g, where G is the set of open frontages of f's kind on f's segment. Each passer has an independent chance at each of them, and the stops are shared in proportion to the chances. When f is alone, c = 1. Two equal frontages each take p(1 − p/2) of the passers.
- **Spending and margin.** s_k is the spend of a visit and m_k its gross margin, the share not paid for goods.
- **Staffing**: staff_k × wage × (open hours + one hour for each day opened, to set up and close).
- **Outgoings**: per m² of trading floor, for any kind but vacant.
- **Refit**: a change of kind to any but vacant, priced per m² of trading floor and amortised. A change to vacant costs nothing.
- **Links**: their cost per metre, amortised.

Rent is a transfer within the district and is left out.

### 4.2 Catalysis, from CLUE's Co-location

The affinity a(i → j) of kind j to kind i is ln CLQ_{i→j}, the co-location quotient of j among i's neighbours (Leslie and Kronenfeld 2011), clamped to [−1, 1]. It is computed by the urban-science domain's `affinity` (on tyche's `colocation`) from CLUE's 10,534 establishments of 2024 in the Hoddle Grid's box. Neighbours are taken within 50 m, with 199 permutations under a free null and a null within 200 m cells.

| CLQ (i → j) | vacant | retail | food and drink | services | office lobby | culture |
|---|---|---|---|---|---|---|
| vacant | 1.15 | 0.77 | 0.92 | 0.98 | 1.02 | 0.73 |
| retail | 0.77 | **2.67** | **1.22** | 0.85 | 0.56 | 1.03 |
| food and drink | 0.90 | **1.22** | **1.54** | 0.85 | 0.81 | **1.11** |
| services | 0.99 | 0.84 | 0.86 | **1.26** | 1.00 | 0.78 |
| office lobby | 1.02 | 0.56 | 0.83 | 0.99 | **1.26** | 0.81 |
| culture | 0.78 | **1.07** | **1.17** | 0.80 | 0.83 | **7.71** |

Bold marks a quotient above one at p ≤ 0.025 under the free null (one-sided). The same-kind quotients and retail with food and drink survive the stratified null too (p = 0.005); food and drink among culture's neighbours does not (p = 0.175).

The weights are w_fg = exp(−d_fg/ξ) over walking distances d_fg ≤ 3ξ. ξ = 42.5 m is the coupling length of the domain's own catalysis representation at its defaults (`urban::affinity::walk_coupling_length`: wandering diffusivity 1 m²/s, a leave rate of one per half hour, frontages 10 m apart). It is taken as that representation states it, and it too is an assumption.

*Caveat.* A co-location quotient measures where kinds are found together, not why. Reading it as catalysis assumes that kinds found together help each other, rather than that they answer to the same terrain or rent.

### 4.3 The Parameters

Every economic parameter is an assumption, chosen to be of the right order for central Melbourne in 2026 and stated so that a study can vary it. None is fitted to data.

| Parameter | retail | food and drink | services | culture | office lobby | Source |
|---|---|---|---|---|---|---|
| propensity π_k (a passer's chance to stop, by day) | 0.006 | 0.010 | 0.0025 | 0.0015 | 0 | assumption |
| multiplier φ_k: day, evening, night | 1, 0.7, 0.3 | 1, 1.3, 0.6 | 1, 0.4, 0.1 | 1, 2.0, 0.5 | — | assumption |
| spend a visit s_k (A$) | 45 | 18 | 55 | 30 | — | assumption |
| gross margin m_k | 0.40 | 0.65 | 0.70 | 0.75 | — | assumption |
| staff while open | 2 | 3 | 2 | 2 | 0 | assumption |
| outgoings (A$/m² of trading floor a week) | 2 | 2 | 2 | 2 | 0 | assumption |

| Parameter | Value | Source |
|---|---|---|
| wage, with on-costs | A$35 an hour | assumption, near the adult rates of the retail and hospitality awards with superannuation |
| refit for a change of kind | A$1,500/m², over 10 years at 5%: A$3.72/m² a week | assumption |
| link (an arcade through a building) | A$20,000/m, over 30 years at 5%: A$24.9/m a week | assumption |
| trading floor | at most the frontage's length × 20 m, and at most the footprint | assumption |
| catalysis γ, and the clamp of 1 + κ | 0.5; [½, 2] | assumption |
| coupling length ξ | 42.5 m | the urban-science domain's catalysis at its defaults (§4.2) |
| affinities a(i → j) | the table of §4.2 | CLUE 2024 (CC BY 4.0) |
| access cap | 400 m | assumption: the walk out of the area when no frontage of a need is open |

### 4.4 The Other Measures

Each is computed in `economy.hpp` from the plan and the flows, the same for every study.

| Measure | Definition | Sense |
|---|---|---|
| spread | The spending S(h) over the week's 168 hours, read on the week's circle: its lines X̂(m) = (1/168) Σ_h S(h) e^{−2πimh/168}, m = 0…84. By Parseval, Σ_m (folded) \|X̂(m)\|² is the mean of S². The spread is the mean line's share of the power, \|X̂(0)\|²/Σ = 1/(1 + CV²), and 168 × spread is the effective number of hours the spending fills. The shares of the daily lines (m = 7j) and of the week's other lines (the modulation of the days, weekday against weekend) are reported beside it | maximised |
| evening share, weekend share | the spending's share in 18–24 h; on Saturday and Sunday | descriptors |
| mix | H(K \| block): the conditional entropy of the kinds of the non-vacant frontages, weighted by their length, given the place they stand in (physis `mixing`; the urban-science domain's Definition 2.6) | maximised |
| access | For each daily need (retail, food and drink, services) and each hour, the walk from each frontage to the nearest open frontage of that need (capped at 400 m). It is weighted by the people passing there and then, and averaged over the week's hours and the three needs | minimised, in m |

A study searches value alone, or value against spread and mix, or value over the evening and weekend shares. Access is reported for every plan.

## 5. The People

### 5.1 The Counters and the Reference Year

The people are calibrated from the City of Melbourne's pedestrian counters, through the data library's `street-sensors` and `street-counts` products.

- **The counters.** There are 60 counters in the Hoddle Grid's box, each placed in the city model's frame by the data library's exact east–north on WGS 84, with the kind of the nearest way (street, laneway, mall).
- **The reference year.** It is the 52 weeks from Monday 27 October 2025 to Sunday 25 October 2026, civil time (AEDT and AEST, daylight saving as the record keeps it). It is the most recent year, and 55 counters cover at least 80% of its hours; the pandemic's years and the City's gap of 2022–2024 are outside it.

`calibrate.cpp` reads each of the 55 with the urban-science domain:

- **The tide** (`urban::rhythm::tide`). The street constituents are fitted by least squares to the counts, people an hour, over the reference year with Victoria's public holidays left out: 46 constituents a counter, the week's lines to six cycles a day and the year's. The fit explains a median 88% of the variance.
- **The mean week.** These are the week's constituents alone (k_M = k_Y = 0), read at the middle of each of the 168 hours. This is the counter's typical week, from which a study's flows are made.
- **The readings.** These are high and low water on the mean weekday and weekend day, and Courtier's form. By it 44 counters are diurnal, 10 mixed and mainly diurnal, and 1 mixed and mainly semidiurnal: the centre has one crest a day, at midday or in the evening, with the commuters' two crests left only at its edges. On Swanston Street at the Town Hall (West), the weekday's high water is 3,400 people an hour at 13:15 and the weekend's 3,841.
- **The spectrum** (`urban::rhythm::spectrum`). Thomson's multitaper estimate (calx `signal`; NW = 4, seven tapers) of ln(1 + counts) is taken over the year's longest stretch with no gap above six hours, with tyche's red-noise level. Its peaks are read on the calendars' lattice of the week and the year. At all 17 counters of the district's box the strongest line is the day's D1 and every one of the four strongest is forced. They are harmonics of the day (D2 to D5) at 15 counters, and a weekly sideband of D2 (D2−W1) at two.
- **The response** (`urban::response::respond`). ln(1 + counts) is regressed on the hours of the week, the year's first two harmonics, rain (the hour and the hour before), the air's temperature and the public holidays, with errors clustered by civil day (Open-Meteo's ERA5 series; Nager.Date's holidays, as the data library keeps them for Victoria).
  - **Rain.** Its median susceptibility is −17% a mm, from −3% to −43%. It is significant (|z| > 2) at 54 of 55 counters.
  - **Holidays.** Their median is +4%, from −56% at Collins Place, where offices empty, to +55% at Hamer Hall, where people come for leisure. They are significant at 33 counters.
  - **Temperature.** Its median is +0.2% a degree.

### 5.2 What a Study Supplies

A study supplies P_f(h; plan) for the mean week, from the counters' mean weeks (`results/sensors.json`) and its own model. Three things are its own, and it states them:

- how the counts become passers at the frontages, which do not hold counters;
- how a plan changes them: opening hours that draw people, links that reroute them;
- what happens on weeks that are not typical. The responses above let it give a wet week or a week with a holiday.

The value and the measures are then this problem's.

## 6. What Every Study Reports

1. **Spatial:** maps of the district, with the frontages by kind and by the value they capture, the flows on the network, and the links made.
2. **Temporal:** hourly and weekly profiles of the people, the spending and the frontages open.
3. **Frequency:** the spectra and harmonic constituents of the flows (from §5.1) and the week's lines of the spending (§4.4), as built and as planned.
4. **People flow:** what the study's model says the people do, and how the plan changes it.
5. **Arrangement:** the buildings and their relations (party walls, one segment, reach) as the plan uses them: where kinds cluster, and which catalysis and competition the plan exploits.
6. **The search:** the algorithms compared at equal budgets over seeds, and the plans found, each with its value, spread, mix and access, beside the plan as built.

## 7. The Machine-Readable Form

| File | What it holds | Committed |
|---|---|---|
| [`prepare.py`](prepare.py) | gathers the inputs from the data library's cache | yes |
| `results/problem-inputs.json` | the district (the designing study's), the counters' places, CLUE's establishments by kind, the parcels' kinds as built, the typical hours | no: extracts of OpenStreetMap (ODbL), made again by `prepare.py` |
| [`calibrate.cpp`](calibrate.cpp) | builds the district's frontages and relations, calibrates the counters, computes the affinities, and writes the problem | yes |
| [`results/problem.json`](results/problem.json) | the problem: kinds and parameters with their sources, the encoding, the value's definition, the objectives, the 61 frontages (keyed by OpenStreetMap's way ids and CLUE's base properties) with their kinds and hours as built, the relations, the links, the affinities, the plan as built | yes |
| [`results/sensors.json`](results/sensors.json) | the 55 counters' calibration: mean weeks, tides and their readings, the 16 largest constituents, spectral peaks, responses; and the spectra of the 17 counters of the district's box | yes (derived from the City's counts, CC BY 4.0) |
| `results/district-geometry.json` | footprints, anchors, segments, blocks, links, for drawing | no: ODbL |
| [`economy.hpp`](economy.hpp) | the shared header: `Problem` (the district, its relations and distances), `Plan`, `Parameters`, `evaluate` (the value and the measures given the flows), `decode` and `encode` (the relaxation), `Kernels` (the mixed encoding's), and JSON helpers | yes |
| [`check_economy.cpp`](check_economy.cpp) | 13 checks of the header | yes |

**Using it.** A study's CMakeLists.txt adds the foundations and domains as the designing study does, includes this directory, and links what `district_economy` links. It constructs `economy::Problem(inputs, &problem)` from the two JSON files, and calls `evaluate(plan, passers)` with its flows.

**The checks** (`check_economy`, all 13 passing):

- **The encoding.** decode(encode(p)) = p for the plan as built and 200 random plans, and every random point of the cube decodes to a valid plan.
- **The distances.** They are zero on the diagonal, symmetric, satisfy the triangle inequality on every triple, and none grows with every link made.
- **Parseval.** The week's lines sum to the mean square, the three shares sum to one, and the spread equals 1/(1 + CV²).
- **The value.** It equals a recomputation written apart for 30 plans under random flows, to a relative 2.6 × 10⁻¹⁵.
- **Closed forms.** Every frontage vacant gives V = 0 and access at the cap. One frontage alone captures P π φ s m, and two of one kind on one segment share 1 − (1 − p)² of the passers.
- **Linearity.** Twice the passers give twice the margin.

```
.venv/Scripts/python.exe research/experiments/2026-11-05-a-district-economy-the-problem/prepare.py
cmake -S research/experiments/2026-11-05-a-district-economy-the-problem -B ~/.cache/science-build/exp-district-economy-problem -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-district-economy-problem -j4
~/.cache/science-build/exp-district-economy-problem/calibrate research/experiments/2026-11-05-a-district-economy-the-problem/results/problem-inputs.json research/experiments/2026-11-05-a-district-economy-the-problem/results
~/.cache/science-build/exp-district-economy-problem/check_economy research/experiments/2026-11-05-a-district-economy-the-problem/results/problem-inputs.json research/experiments/2026-11-05-a-district-economy-the-problem/results/problem.json
.venv/Scripts/python.exe research/experiments/2026-11-05-a-district-economy-the-problem/figures.py
```

(The programs are run from the repository's root, where the inputs' paths to the cache resolve. Calibration takes about 10 s.)

## 8. Code Version

| Library | Commit |
|---|---|
| The Science Infrastructures (root) | `512477b` |
| heuresis | `5730e70` |
| calx | `90aa710` |
| morphe | `14b3fb8` |
| physis | `e2b98d8` |
| tyche | `b4034b8` |
| domains/urban-science, domains/architecture | unchanged since `485f82c` |
| the designing study (its `district-inputs.json`, `lanes.hpp`, `walking.hpp`) | as committed at `4153cb9` |

Built with g++ 11.5 and CMake 3.28 in WSL, in Release, without warnings from this problem's code.

## 9. Data and Provenance

| Data | Source | Licence | Retrieved |
|---|---|---|---|
| buildings, ways, the walking network | the data library's city model of Melbourne, from OpenStreetMap, through the designing study | ODbL 1.0, © OpenStreetMap contributors | 2026-10-31 |
| hourly pedestrian counts, the counters' places | City of Melbourne, Pedestrian Counting System (`melbourne-pedestrian-counts`, `melbourne-pedestrian-archive`; the data library's `street-sensors` and `street-counts`, built 2026-10-31) | the counts' metadata state no licence; the City's open data are CC BY 4.0, and they are used on that understanding, attributed to the City | 2026-10-30 |
| establishments, floors | City of Melbourne, Census of Land Use and Employment (`melbourne-clue-establishments`; joined to the city model by the data library's `city-events-melbourne`, built 2026-11-04) | CC BY 4.0 | 2026-11-04 |
| declared opening hours | OpenStreetMap, the data library's `street-holdings` (built 2026-10-31) | ODbL 1.0 | 2026-10-30 |
| weather | Open-Meteo, from Copernicus ERA5 (`open-meteo-hourly`, the data library's `street-weather-melbourne`) | CC BY 4.0 | 2026-10-30 |
| public holidays | Nager.Date (`nager-holidays`, the data library's `street-holidays-victoria`) | Nager.Date's API terms: private and non-profit use | 2026-10-30 |

## 10. Figures

`figures.py` draws two figures from `results/` (and the uncommitted geometry), in the house style.

| Figure | What it shows |
|---|---|
| `search-60-district-economy-the-problem` | the district's frontages by kind as built, the party walls, the reach of catalysis, the candidate links and the counters, and the affinities of the kinds |
| `search-61-the-people-of-the-district` | the counters' mean weeks, the week's power by line, their spectra on the calendars' lattice, and their response to rain and holidays |

## 11. Caveats

- **The economics is a sketch.** Shoppers are passers with fixed chances to stop. They have no purpose, no budget across stops, and no choice of district. A real frontage's takings depend on its brand, its price and its rent, none of which is here.
- **Competition and catalysis are local.** Competition acts only among frontages of one kind on one segment, and catalysis only within 128 m. A shop across the street or round the corner competes with nothing.
- **The as-built plan is read, not observed.** Kinds come from establishments with no floor, and hours from five other streets' medians. Eleven parcels that CLUE does not record are read as vacant.
- **One frontage a parcel.** Corner buildings have one shopfront; a building's arcade or second street front is not a frontage.
- **The mean week.** The value is of a typical week. Holidays, rain and the seasons enter only through the responses a study chooses to apply (§5.2).
- **Links cost little.** At A$500–2,300 a week, a link is cheap beside the A$ tens of thousands a busy frontage turns over. Whether a link is worth making depends almost entirely on what the study's model of flow says it does to the people.

### 11.1 What the Domains Gave and Where the Problem Worked Round Them

- **`urban::events::city_network` assigns every outer wall to the nearest segment it faces.** A building's back on a service lane is then one of its frontages, and for the Town Hall and others the largest one. The problem takes a street's or mall's facade first (§2.2). A frontage of the domain's own, with its length and its main street, would serve every study of ground floors.
- **`urban::landscape::attach_store` attaches to whichever segment is nearest, a new link included.** With a passage made, a frontage's anchor could move onto the passage and its distances grow. The problem attaches the frontages only to the segments that were there before the links (`Problem::distances`), and the check that no distance grows with every link made holds.
- **Eleven of the 61 parcels have no CLUE base property** joined by the data library's city-events product, among them the Council Chamber and Federal Pacific House. They are read as vacant (§3.2).
