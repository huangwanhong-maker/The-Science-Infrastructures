# Designing a District

*Status: draft (2026-11-05). Research.* A case study of the search foundation heuresis ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)) on a real place: a few blocks of Melbourne's Hoddle Grid, redesigned as one problem with decisions of several kinds and measures from several sciences.

This is a study of methods on a real place, not a planning proposal. Its models are the domains' (and three provisional ones of its own, §2.6), its constraints are assumptions, and nothing in it says what should be built there.

## 1. The Question

The design of a district is a standard example of a heterogeneous, many-objective problem. Its decisions are of different kinds:

- how tall each building stands (an integer);
- what it is used for (a category);
- which links the street network gains (a yes or no).

Its measures come from different sciences, and they disagree:

- floor area is the developer's;
- the winter sun on the streets is the walker's;
- the heat a building needs is the engineer's;
- the walk to a shop or a park is the resident's;
- the mix of uses is the planner's.

No design is best in all of them, so the answer is a front of designs, each better than every other in something.

The study asks four questions.

1. What front do heuresis's many-objective algorithms find when the decisions are a real district's heights, uses and links, and every measure is the model of the domain that owns it?
2. How do NSGA-III, MOEA/D and random search compare on it, by hypervolume?
3. What does a two-objective search (NSGA-II) find that the five-objective one misses?
4. What does a repertoire (MAP-Elites) show of the designs possible at each density and share of open space?

*Rationale.* It matters because the infrastructure's domains were built to be composed. A search that reads five models through one ask–tell interface tests that composition at once: of architecture and urban science, of the astronomy domain's Sun, and of the data library's city model and census.

## 2. Method

### 2.1 The District

- **The district.** It is Landscape Dynamics' district of the Urban Science Workbench: the Hoddle Grid's middle, Elizabeth to Russell Street and Flinders Street to Lonsdale Street.
- **The walking network.** It is built as the workbench builds it, and `prepare.py` follows `applications/urban-science-workbench/src/district.js` line for line:
  - the box is carried into the city model's local frame by the data library's exact east–north on WGS 84 (`street_records.enu`);
  - the city model's walkable ways are cut to the box and split at every shared point.
  - The port was checked once, by running the workbench's own `networkOf` on the same city model in Node: the two networks are identical, 356 junctions and 399 segments.
  - The study does not import the workbench, because nothing may depend on an application (CLAUDE.md §3.4).
- **The domain's complex.** The urban-science domain builds the district from that network (`urban::landscape::planarised`, then `district`). It drops no segment. Its complex has 356 0-cells, 399 1-cells and 37 2-cells, and it is regular.
- **The redesign area.** The study redesigns the blocks between Swanston and Russell Streets, Flinders Lane and Bourke Street. They are bounded by the junctions where those streets meet; where a street is two carriageways, the inner junction is used.
- **The parcels.** The 61 buildings of the city model whose footprint's centre lies inside the area are the parcels; the footprint of each is kept as built.
  - The city model's 42 `building:part` entries in the box and its margin are left out. A part lies inside its outline, and a redesign gives a parcel one height.
  - 50 of the parcels are joined to a CLUE base property by the data library's city-events product. Their floors above ground and predominant use come from CLUE's latest census, 2023 or 2024.
  - The other 11 take the city model's height, divided by its 3.2 m storey.
- **The context.** The other 524 outlines in the box, widened by the product's 60 m margin, stand as built. 431 of them take CLUE's floors × 3.2 m, the rest the city model's height. They cast shadows, and they are drawn.
- **Existing open spaces.** The city model's open spaces inside the area are kept as public open space in every design: City Square's plaza and six small parks, 3,621 m² in all.
- **The land.** The parcels and these spaces are 53,561 m².

### 2.2 The Decisions

A design is one genotype of three parts, composed by hand as heuresis's `variation` README (§6) says mixed spaces are.

- **Storeys.** For each parcel i, its storeys sᵢ ∈ {1, …, bᵢ}, an integer. The bound bᵢ is the lesser of two figures:
  - 48: the 90th percentile of the floors of the 248 developments that the City's Development Activity Monitor records as completed in the Hoddle Grid's box;
  - 8 √(footprint) / 3.2 m: a slenderness bound, an assumption.
- **Use.** For each parcel, its use uᵢ ∈ {dwellings, offices and other workplaces, shops (shops, food and entertainment), open space}. CLUE's uses map to these as `prepare.py` states. A parcel given to open space has no building: its ground becomes public.
- **Links.** Eleven candidate links of the walking network, each made or not (pₖ ∈ {0, 1}).
  - **Eight blind lanes carried through** to the street beyond them, 19–91 m long. Among them: Ramsden Place to Collins Street, Rainbow Alley to Bourke Street, Athenaeum Place and Baptist Place to Collins Street, and Portland Lane and Donaldson Lane to Russell Place.
  - **Three passages through blocks that are 2-cells,** 31–39 m long. They are placed across each block's middle by the domain's `plan_passage` and opened by its `open_passage`.
  - Every link is an arcade at ground level that takes 4 m of ground floor along its length.

Altogether a design has 61 integer, 61 categorical and 11 binary variables.

The kernels are these:
- **Sampling.** Each part is drawn uniformly.
- **Mutation.** heuresis's `creep_mutation` moves the storeys by at most 6, at 2 parcels in 61 on average. Its `reset_mutation` redraws a use at 1 parcel in 61. Its `standard_bit_mutation` flips a link at 1 in 11.
- **Crossover** (probability 0.9).
  - A parcel's storeys and use travel together.
  - The parcels are ordered along the area's long axis, so `two_point_crossover` exchanges a stretch of the area between two cuts across it. This assumes, as k-point crossover does, that neighbouring positions belong together: neighbouring parcels shade one another and share walls.
  - The links cross uniformly.

### 2.3 The Measures, Each from Its Model

| Measure | Definition | Model |
|---|---|---|
| floor area (maximised) | Σ over built parcels of the footprint's area × storeys, less 4 m × length for each link made | `architecture::space::massing`: its footprint areas after pooling corners and inserting T-junctions |
| winter sun (maximised) | the sunlit area-hours of the public ground (all unbuilt ground sampled every 4 m over the area and the far half of its streets, the existing open spaces, and the parcels given to open space) at the eight half-hours between 10:00 and 14:00 AEST on 21 June 2026 | the Sun: the astronomy domain's subsolar point and altitude (`domains/astronomy/js/sun.js`), as the workbench's clock takes it; the shadows: the study's own model (`sunlight.hpp`, §2.6) |
| heating (minimised) | the heat the buildings need on 21 June, held at 20 °C, per m² of floor: each hour, Σ max(0, G_air (20 − T_air) + G_ground (20 − T_ground) − S), with the Sun's gains S | `architecture::environment`: `massing_heat` (conductances by face: facades 1.2, roofs 0.6, ground 0.4 W/m²K; party walls shared and, at one temperature, carrying nothing) and `massing_gains` (sol-air gains through glass and opaque fabric, Meinel's clear sky) at the Sun of each hour |
| walking (minimised) | the mean walk from a home to its nearest shop, plus the mean walk to its nearest open space, homes weighted by their floor | `urban::landscape`: the district built from the network with the lanes carried through, passages opened by `open_passage`, every parcel and open space placed on its nearest segment by `attach_store`; distances by the study's own Dijkstra (`walking.hpp`, §2.6) |
| mix (maximised) | H(U \| X), the conditional entropy of the floor's use given the block it stands in: the mean local mix | the urban-science domain's Definition 2.6 (theory/street-oceanography.md §2.3), computed by `physis::statistical_mechanics::mixing` on floor areas (the domain's `mixing_of` counts holdings one each, §8) |

The air is the mean of each civil hour over the Junes of 2009–2025 at Melbourne, from 7.8 to 13.4 °C, taken from Open-Meteo's ERA5 series in the data library. The ground is at the record's mean, 14.9 °C.

The search minimises five objectives, each relative to the district as built: −floor/floor₀, −sun/sun₀, heat/heat₀, walk/walk₀ and −H/ln 3. The district as built is then the point (−1, −1, 1, 1, −0.41).

### 2.4 The Constraints and Assumptions

The constraints are the study's assumptions, not the planning scheme's. Each is satisfied at ≤ 0, and the searches handle them by Deb's constrained dominance (heuresis `multiobjective`).
- open space at least 15% of the land;
- homes at least 20% of the floor;
- shops at least 5% of the floor and at most 20%.

The district as built meets only the last, with 15.2% shops. Its open space is 6.8% of the land and its homes 4.1% of the floor.

Other assumptions:
- a storey of 3.2 m (the city model's);
- buildings as prisms on flat ground, their footprints kept;
- uses that do not change a building's heating;
- links that are arcades;
- destinations, shops and open space, counted inside the area only: the area answers for itself.

### 2.5 The Searches

Every search runs on heuresis's ask–tell protocol. Evaluations run on a pool of four threads, and every draw is keyed by its purpose (`search::Draws`), so a run is the same on any number of threads.

| Search | Objectives | Population | Generations | Evaluations a run | Seeds |
|---|---|---|---|---|---|
| NSGA-III | five | 210, on Das and Dennis's 210 directions (m = 5, p = 6) | 210 | 44,100 | 1, 2, 3 |
| MOEA/D, Tchebycheff | five | the same 210 directions; 20 neighbours, at most 2 replaced | 210 | 44,100 | 1, 2, 3 |
| random search | five | — | — | 44,100 | 1, 2, 3 |
| NSGA-II | floor area and winter sun | 100 | 210 | 21,000 | 101, 102, 103 |
| MAP-Elites | winter sun its quality | a grid of 16 × 12 cells: plot ratio 0–32, open share 0.05–0.65 | — | 84,000 (1,000 random first; half bred by crossover; 64 per ask) | 7 |

- **The repertoire.** MAP-Elites keeps the homes and shops constraints; the open share is its descriptor, not a constraint.
- **Hypervolume.** The five-objective hypervolumes are exact (WFG). Each objective is normalised by the least and greatest values on the union of the final feasible fronts, with the reference point 1.1 in each.
- **The fronts compared.** The front of NSGA-III and of MOEA/D is the feasible non-dominated set of its final population. The front of random search is its whole feasible non-dominated archive, which is generous to the baseline.
- **Two objectives.** In two objectives, NSGA-II's fronts and NSGA-III's fronts projected onto the same plane are normalised together, with the reference 1.1.
- **The chosen designs.** Six designs are taken from NSGA-III's front (seed 1): the best in each measure, and a balanced design, whose largest normalised objective is the least.

### 2.6 What the Study Writes Itself

Three provisional pieces fill gaps in the domains and foundations. Each is marked in its header and reported in §8.
- **`sunlight.hpp`: shadows on the ground.**
  - A sample is in the sun when the ray toward the Sun clears every prism.
  - A prism of height h blocks the ray exactly when t_in tan α < h, where t_in is the ray's first entry into the footprint.
  - So the thresholds are computed once, and a design is evaluated by comparisons.
- **`walking.hpp`: shortest walks.** Dijkstra's algorithm, from sources placed along segments.
- **`lanes.hpp`: blind lanes carried through.** An edit of the network the domain is given. The domain then builds the district again from it, and checks it regular.

## 3. Code Version

| Library | Commit |
|---|---|
| The Science Infrastructures (root; heuresis a submodule) | `4153cb9` |
| heuresis | `6b3cc88` |
| calx | `90aa710` |
| morphe | `14b3fb8` |
| physis | `e2b98d8` |
| tyche | `218cf0b` (heuresis's CMake adds it; nothing of it is compiled into the study) |
| domains/urban-science, domains/architecture | unchanged since `485f82c` |

The C++ was built with g++ 11.5 and CMake 3.28 in WSL, in Release, without warnings from the study's own code; morphe's `modular.hpp` warns of `__int128` only under `-Wpedantic`, which the study leaves off. Node 22.17 ran `sun.mjs`, and the repository's Python 3.12 ran `prepare.py` and `figures.py`.

```
.venv/Scripts/python.exe research/experiments/2026-11-05-designing-a-district/prepare.py
cmake -S research/experiments/2026-11-05-designing-a-district -B ~/.cache/science-build/exp-designing-a-district -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-designing-a-district -j4
~/.cache/science-build/exp-designing-a-district/designing_a_district results/district-inputs.json results --threads 4
.venv/Scripts/python.exe research/experiments/2026-11-05-designing-a-district/figures.py
```

(The study's binary is run from its directory. `--scale` multiplies every budget, and `--profile` times 400 random designs and stops.)

## 4. Parameters

All parameters are in §2.2–2.5 and in `district.hpp` (`Assumptions`). The ground is sampled every 4 m, giving 5,831 samples, and the shadows reach at most 900 m.

The final run took 711 s. One evaluation cost about 3.5 ms on four threads, with walking (while its cache of networks filled), heating and massing the largest parts (`results/run.log`).

## 5. Data and Provenance

| Data | Source | Licence | Retrieved |
|---|---|---|---|
| buildings, ways, open spaces | the data library's city model of Melbourne (`city-model-melbourne.json`, built 2026-10-31) from OpenStreetMap | ODbL 1.0, © OpenStreetMap contributors | 2026-10-31 |
| floors and uses | City of Melbourne, Census of Land Use and Employment (CLUE), joined to the city model by the data library's `city-events-melbourne.json` (built 2026-11-04) | CC BY 4.0 | 2026-11-04 |
| the bound on storeys | City of Melbourne, Development Activity Monitor (`melbourne-development-activity`, profile `hoddle-grid`) | CC BY 4.0 | 2026-11-04 |
| air temperature | Open-Meteo, from Copernicus ERA5 (`street-weather-melbourne.json`) | CC BY 4.0 | 2026-10-30 |
| the Sun | the astronomy domain's low-precision theory (Astronomical Almanac, §C5), 0.01° | — | — |

`results/district-inputs.json` holds the prepared extract, credited as above: 61 parcels, 524 context outlines and the network. It and `results/district-cells.json` (the district's blocks, holes and segments) are extracts of OpenStreetMap under ODbL, so they are not committed (`.gitignore`): `prepare.py` makes the first from the data library's cache, and the study's program writes the second. The other results hold only the study's own decisions and measures, keyed by OpenStreetMap's identifiers.

## 6. Results

### 6.1 The District as Built

| Measure | Value |
|---|---|
| floor area | 440,046 m² (plot ratio 8.2) |
| winter sun | 50,512 m²·h |
| heating | 72.7 Wh per m² of floor on 21 June |
| walk | 234.2 m (37.8 m to a shop, 196.4 m to open space) |
| mix H(U \| X) | 0.453 nats |

### 6.2 The Five-Objective Front

**NSGA-III's front.** NSGA-III's final front (seed 1) holds 209 designs. 192 of them are better than the district as built in all five measures, and all meet the constraints that the district as built fails.

**The six chosen designs.**

| Design | Floor (m²) | Plot ratio | Winter sun (m²·h) | Heating (Wh/m²) | Walk (m) | Mix (nats) | Open | Homes | Shops | Links |
|---|---|---|---|---|---|---|---|---|---|---|
| as built | 440,046 | 8.2 | 50,512 | 72.7 | 234.2 | 0.453 | 6.8% | 4.1% | 15.2% | 0 |
| the most floor | 1,573,871 | 29.4 | 50,400 | 52.0 | 97.9 | 0.804 | 25.3% | 36.1% | 18.2% | 3 |
| the most winter sun | 744,850 | 13.9 | 196,400 | 46.7 | 174.6 | 0.688 | 63.9% | 64.6% | 9.3% | 8 |
| the least heating | 768,488 | 14.3 | 179,824 | 41.5 | 34.9 | 0.594 | 64.6% | 32.4% | 13.1% | 6 |
| the shortest walks | 840,811 | 15.7 | 154,168 | 48.9 | 22.7 | 0.656 | 56.5% | 22.4% | 11.2% | 4 |
| the most mixed | 774,372 | 14.5 | 126,440 | 57.0 | 119.1 | 1.033 | 56.9% | 41.3% | 18.2% | 4 |
| a balanced design | 1,068,546 | 19.9 | 164,856 | 47.7 | 156.8 | 0.870 | 50.0% | 52.6% | 14.2% | 5 |

**Where the front puts height.** The front gives the north half-block, Little Collins to Bourke Street, to the public, and builds tall to the south.
- Over the 209 designs, the median share of that half-block's land given to open space is 94%. Its median built height is 4 storeys, and it holds a median 0.5% of the floor.
- The half-block Collins to Little Collins gives 34% of its land and builds at a median of 42 storeys.
- The half-block Flinders Lane to Collins gives 19% and builds at 29 storeys.

*Remark.* In the Southern Hemisphere's winter the Sun stands low in the north: 28.5° at noon on 21 June. Ground on the northern edge of the area is shaded by nothing of the area, and towers to its south cast their shadows over each other and over streets already in shadow. In the domains' models, a winter-sunlit district is one that keeps its northern edge low.

**The links.** Each link is made in between 28% and 58% of the front's designs, so no link is decisive in these measures. The walks depend far more on where the shops and open spaces are than on the links.

**The extremes.** Across the three seeds, MOEA/D's fronts reach further than NSGA-III's: up to 2,088,522 m² of floor against 1,723,466 m², 273,280 m²·h of winter sun against 211,048, and 35.2 Wh/m² of heating against 41.5. NSGA-III's designs spread through the middle of the front (§6.3).

### 6.3 The Algorithms Compared

Normalised hypervolume of the final feasible fronts, five objectives:

| Search | Seed 1 | Seed 2 | Seed 3 | Mean | Front size |
|---|---|---|---|---|---|
| MOEA/D | 0.685 | 0.703 | 0.686 | 0.691 | 111, 134, 127 |
| NSGA-III | 0.567 | 0.567 | 0.592 | 0.575 | 209, 209, 209 |
| random search | 0.143 | 0.139 | 0.136 | 0.139 | 334, 325, 388 (archive) |

- **Random search.** Both evolutionary searches far exceed it: four to five times its hypervolume at the same 44,100 evaluations. Random search's archive barely improves after the first few thousand evaluations.
- **MOEA/D.** It leads at the end. Its hypervolume rises fastest and is within a few hundredths of its final value by about 30,000 evaluations. Its 210 Tchebycheff directions include the corners of the simplex, and its front gathers at the extremes of the five measures, which the hypervolume rewards.
- **NSGA-III.** It is still rising at 44,100 evaluations, from 0.41–0.43 at 21,000 to a mean of 0.575. It spreads its 209 designs through the middle of the front. Which leads at a larger budget is not settled here.

**Two objectives.** Asked only for floor area and winter sun, NSGA-II finds much better trade-offs in that plane than NSGA-III's five-objective fronts projected onto it.
- **Hypervolume in that plane.** NSGA-II scores 0.755, 0.710 and 0.711 (mean 0.725 ± 0.021) on 21,000 evaluations. NSGA-III's projection scores 0.490, 0.502 and 0.536 (mean 0.509 ± 0.019) on 44,100.
- **The seeds agree.** NSGA-II's three seeds draw nearly one curve, from 652k m² of floor with 228k m²·h of winter sun to 2,063k m² with 49.5k m²·h.
- **The end of the curve.** At that end the area holds 4.7 times its floor and still has the winter sun it has as built. That sun includes the sun on the 15% of the land the constraint requires to be open.
- **The cost.** That focus is paid for in what the plane does not show. NSGA-II's designs walk 114–412 m and reach at most 0.91 nats of mix, where NSGA-III's front reaches walks of 4.6 m and mixes of 1.04 nats.

### 6.4 The Repertoire

MAP-Elites filled 156 of its 192 cells, 81%, in 84,000 evaluations, every elite feasible. Its QD-score is 426.8: the sum over the filled cells of each elite's winter sun relative to the district as built. Coverage passed 75% at about 30,000 evaluations and grew slowly after.

The elites' winter sun ranges from 0.81 times the district as built's, in the densest cell with the least open space, to 5.1 times.
- **At the density as built.** Near the district as built's plot ratio, the best design found keeps 1.47 times its winter sun: plot ratio 9.2 and open share 9.7%, the cell holding the district as built (6.8%).
- **With more open space.** With 15.0% open space it keeps 1.88 times (plot ratio 9.7). With 19.9% it keeps 2.18 times (plot ratio 8.6).
- **The rest of the repertoire.** The repertoire also holds designs denser or sparser than any on the five-objective front. Fourteen of NSGA-III's 209 designs lie beyond its grid, with more than 65% open space.

## 7. Figures

The figures are drawn by `figures.py` from `results/` only. `gallery.md` holds them with their captions.

| Figure | What it shows |
|---|---|
| `search-50-district-as-built` | the district as built in three dimensions, by CLUE's uses, with its shadows at noon on 21 June, and its measures |
| `search-51-decisions-and-sun` | the parcels' storey bounds, the domain's walking network, blocks and holes, the eleven candidate links; the hours of winter sun on the ground as built |
| `search-52-front-parallel` | NSGA-III's front in parallel coordinates, the chosen designs and the district as built |
| `search-53-front-pairs` | the ten pairs of measures for NSGA-III's and MOEA/D's fronts |
| `search-54-six-designs` | the six chosen designs in three dimensions |
| `search-55-two-objectives` | NSGA-II's fronts in floor area and winter sun beside NSGA-III's projection, with three designs along them |
| `search-56-hypervolume` | hypervolume against evaluations for NSGA-III, MOEA/D and random search |
| `search-57-repertoire-atlas` | MAP-Elites' repertoire as an atlas of districts |
| `search-58-repertoire-map` | the repertoire's winter sun over plot ratio and open share; coverage and QD-score |
| `search-59-where-height-goes` | where the front puts its height, and what it gives to open space |

`results/designs.json` holds the district as built and the six chosen designs for the Urban Science Workbench's GIS. Each is keyed by the city model's building ids (OpenStreetMap ways) and CLUE's base properties, with storeys, heights, uses, the links made (with their endpoints in the local frame) and the measures.

## 8. Caveats

### 8.1 What the Models Leave Out

- **The shadow model.** Shadows are cast by prisms on flat ground. The study's shadow model leaves out:
  - diffuse light from the sky, and light reflected from facades;
  - trees, verandahs and awnings;
  - buildings beyond the context box.

  "Public ground" is all ground not covered by a footprint. It includes narrow gaps between footprints that do not quite touch in the city model, which are nearly always in shade.
- **The heating model.** The domain's heating leaves out:
  - ventilation and infiltration, a large part of real heating;
  - internal gains of people, lighting and equipment, which differ by use;
  - shading of facades by neighbours (the domain's stated limit);
  - diffuse and reflected sunlight.

  It is a measure of envelope losses less direct solar gains, at one fabric for every building. Taller and terraced buildings lose less per m² of floor in it, as they do in fact, but the absolute figures are not a building's energy bill.
- **The walking model.** Walking counts only destinations inside the area, along the network's centrelines. It has:
  - no crossings, gradients or waiting at lights;
  - no shop or park beyond Swanston, Russell, Flinders Lane or Bourke Street.
- **The mix.** It counts floor area by three uses in 49 places (the domain's blocks and holes). It knows nothing of what the floor holds.
- **The design.** Footprints are kept as built, a use fills a whole building, and links are arcades costing 4 m of ground floor. Heritage, structure, cost, ownership, light inside buildings, wind and the planning scheme are all absent. The front would move with any of them.
- **The budgets.** The budgets are modest: 44,100 evaluations for a problem of 133 variables. NSGA-III had not converged (§6.3).
- **The district as built.** Its heights are CLUE's floors × 3.2 m: 50 of 61 parcels and 431 of 524 context buildings. The rest are the city model's heights, many of them defaults of their kind.

### 8.2 Gaps Found in the Domains and Foundations

- **architecture `environment`.** It has no model of shadow and daylight: it lists them as planned. Its `massing_gains` ignores shading by neighbours. The study wrote its own ground-shadow model, `sunlight.hpp`. Shadows on the ground and on facades belong in `environment`, and its gains would then see them.
- **morphe `combinatorics`.** It gives distances in edges only, by breadth-first search. There is no single-source shortest distance on a weighted multigraph, although morphe's `MinPlus` semiring gives its algebra. The study wrote Dijkstra's algorithm, `walking.hpp`.
- **urban-science `landscape`.**
  - Passages open only through 2-cells, and a block holding a dead-end lane is left as a hole. In this area the five largest blocks, 7,845–24,613 m², are holes, so the domain could open passages only through three small blocks.
  - There is no edit "carry a blind lane through". The study wrote one on the network, `lanes.hpp`, and let the domain rebuild the district.
- **urban-science `street`.** `mixing_of` counts holdings one each; a weighted form would serve floor areas. The study called physis's `mixing` directly.
- **architecture `form` and `rules`.** Both are planned. Floor area by storeys, height and slenderness bounds and plot ratios are the study's. They belong there.
- **Two observations, not gaps.**
  - `architecture::environment::massing_gains` recomputes the cofaces of the whole complex at every call, so a day of gains rebuilds them at each hour.
  - heuresis's MOEA/D scalarises the objectives as they come, without normalising them. The study gave it objectives relative to the district as built, so that their scales were comparable.

### 8.3 A Bug Found in the Workbench

**`footprintsIn` sends parts with their outlines.** `footprintsIn` in `applications/urban-science-workbench/src/district.js` passes the city model's `building:part` entries (kind `part:…`) to the engine's massing together with the outlines that hold them. In Landscape Dynamics' district, 30 of the 440 footprints sent are parts. One is `way/140838470`, a `part:office` 169 m tall, inside its outline `way/46398815`. The massing then holds overlapping 3-cells.

Reproduction: run the workbench's `networkOf` and `footprintsIn` in Node on `data/cache/products/city-model-melbourne.json` with the Melbourne area's district box, then count the returned keys whose city-model `kind` starts with `part:`. A filter on `kind` (or on `part_of`) would fix it. The study does not edit the application.
