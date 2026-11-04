# A District's Economy without a Model

*Status: draft (2026-11-05). Research.* Study (1) of the three on [the shared problem](../2026-11-05-a-district-economy-the-problem/README.md), the ground floors of the blocks between Swanston and Russell Streets, Flinders Lane and Bourke Street: the people from the data alone, and the value searched as a black box on heuresis ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)).

This is a study of methods on a real place, not a planning proposal. The economics are the shared problem's assumptions, the flows are interpolations and one stated correlation, and nothing in it says what any frontage should hold.

## 1. The Question

What does a plan of a district's economy look like when the people are taken from the counts and nothing else? This study fixes one answer to that and asks four questions of it.

- **The flows.** No conservation law, no network of flow and no circulation. The passers at each frontage are interpolated from the counters, and their change under a plan is an empirical rule stated as such.
- **The searches.** The value is treated as a black box.

1. How do heuresis's searches of different kinds compare on it, at equal budgets over seeds? They are the genetic algorithm on the mixed encoding, CMA-ES on the canonical relaxation, Bayesian optimisation, MAP-Elites, NSGA-II and random search.
2. What plan do they find, in space, in time and in frequency?
3. What does the plan give up in spread over the week and in mix, and what can the week hold in the evening and at the weekend?
4. Where is a plan without a model blind?

*Rationale.* Studies (2) and (3) give the same problem more structure: a fluid with its circulation, then a field with its action. This one is their baseline. What it finds without a model, and what it cannot see, is what the models must earn their keep against.

## 2. Method

### 2.1 The People from the Data

The study's own provisional code ([`flows.hpp`](flows.hpp)) reads the people from the data in two steps.

**Interpolation (kriging).**
- **The field.** At each hour h of the mean week, y_i(h) = ln(1 + c_i(h)) at the 55 calibrated counters is a Gaussian random field in the plane. Its mean μ(h) is the hour's mean over the counters. Its deviation has covariance s²[(1 − η) exp(−d/ℓ) + η δ]: an exponential (Matérn ν = ½) with a nugget η, which absorbs the difference between a street's two footpaths and each counter's own scatter.
- **The conditioning.** The deviation at a frontage's anchor is its posterior mean given the counters, by tyche's Gaussian conditioning (`tyche::spatial::GaussianConditioning`; simple kriging, Cressie 1993 §3.4.1). One factorisation serves all 168 hours, each hour's innovations whitened by it.
- **The passers** are the posterior median of the count, exp(μ̂) − 1.
- **The hyperparameters.** ℓ and η are chosen on a grid (ℓ from 20 to 900 m, η from 0 to 0.6) by the marginal likelihood, with s² profiled out. The 168 hours are taken as independent replicates of the field, which is an approximation.
- **The check.** Each counter's week is predicted from the other 54, leaving one out.

**How a plan changes them (an empirical rule).**
- **The regression.** Across the counters, ln(1 + c_i(h)) is regressed on ln(1 + A_i(h)). A_i(h) is the number of CLUE's 2024 trading establishments (retail, food and drink, services, culture) within 50 m of the counter that are open at h by their kind's typical hours (the shared problem's hours as built). The 168 hours' means are absorbed, and errors are clustered by counter (tyche's least squares). The slope ε is a cross-sectional elasticity.
- **Its use.** A plan multiplies a frontage's passers by ((1 + A_f(h; plan))/(1 + A_f(h; as built)))^ε, clamped to [½, 2]. A_f(h; plan) is CLUE's count about the frontage with the area's ground floors replaced by the plan's: the trading ground floors within 50 m of f open at h as built are taken out, and the plan's put in.
- **Links.** A link made changes nobody's passers. It enters the value only through its cost, and the shared measures only through access.

### 2.2 The Value

The value and every other measure are the shared problem's (`economy.hpp`, `Problem::evaluate`), under the flows of §2.1: the week's gross operating surplus, the spending's spread over the week, the mix H(K | block), access, and the evening and weekend shares.

### 2.3 The Searches

Every search runs on heuresis's ask–tell protocol, on a pool of four threads, with every draw keyed by its purpose (`search::Draws`). Every evaluation of every run is kept by its ticket: its value, spread, mix, access, and evening and weekend shares. Each run is then read by the same measures: its best value along the way, the front of everything it evaluated, and the repertoire of everything it evaluated.

| Search | Encoding | Settings |
|---|---|---|
| GA (heuresis `genetic`) | the mixed encoding (61 kinds, 244 hours, 11 links) | μ = λ = 100 (20 at the small budget); tournaments of two; crossover 0.9; 2 elites; the shared problem's kernels: a kind redrawn at 1 frontage in 61, an hour crept by ≤ 3 at 1 gene in 61, a link flipped at 1 in 11, two-point crossover along the area with uniform crossover of the links |
| CMA-ES (heuresis `strategies`) | the canonical relaxation, [0, 1]^621 | separable (Ros and Hansen 2008), λ = 23, σ₀ = 0.3, the box by Hansen's transformation; a random mean in [0.1, 0.9]^621; restarted from a new mean when it stops |
| BO (heuresis `models`) | a random embedding of the relaxation (REMBO; Wang et al. 2016) | z ∈ [−√10, √10]^10, u = (clip(Az, −1, 1) + 1)/2 with A a 621 × 10 Gaussian matrix keyed by the seed; Matérn 5/2 with a length scale per coordinate; expected improvement; 22 initial points, then 4 a batch by the kriging believer |
| MAP-Elites (heuresis `diversity`) | the mixed encoding | a 12 × 12 grid of the evening share (18–24 h) × the weekend share of the spending over [0, 0.6]²; quality the value; 500 random plans first (100 at the small budget), half then bred by crossover, 32 a batch (8) |
| NSGA-II (heuresis `multiobjective`) | the mixed encoding | value, spread and mix, all maximised; μ = 100 (20) |
| random search | the mixed encoding | uniform plans, 32 a batch |
| *GA from the plan as built* | the mixed encoding | the GA above, its first generation the plan as built mutated five times over: outside the comparison of equal knowledge, as a planner would start |

The budgets are equal within each comparison, over seeds 1–5:
- **large:** 20,000 evaluations for the GA, CMA-ES, MAP-Elites, NSGA-II and random search;
- **small:** 400 evaluations for all six. Bayesian optimisation refits its surrogate at each batch, and cannot afford the large budget.

Every search starts from uniform random plans (CMA-ES from a random mean), not from the plan as built. *Remark.* A planner would start from the plan as built. Here the searches are compared as black boxes, from the same knowledge of the space.

### 2.4 What Is Read from the Runs

- **Value.** The best value found against evaluations: the median and range over seeds.
- **Fronts.** The non-dominated plans in (value, spread, mix) among all that a run evaluated, and their exact hypervolume (WFG). Each objective is normalised over the union of the large budget's fronts. The value's axis runs from the best found to the plan as built, so that a plan worse than as built adds nothing; the reference point is 1.1.
- **Repertoires.** Every evaluation of a run is binned into MAP-Elites' grid: the cells reached, and a QD-score, Σ over cells of max(0, the best value in the cell), in A$ thousand a week.
- **Blindness.**
  - **Links.** The people's shortest walks between the frontages are found by Dijkstra on the urban-science domain's district, for every pair weighted by the product of their passers as built, with each link made alone. Two things are read: the share of the walks that pass through the link, and how many frontages' own segment would gain or lose a tenth or more of the walks past it. Beside these stands the value's change with the link made.
  - **Conservation.** The best plan's people against those as built.
  - **Lanes.** The kriged passers on the area's lanes against the laneway counters'.
  - **Atypical weeks.** The value under a Monday holiday and a wet week, by the median responses of the district's counters (from the shared calibration).

## 3. Code Version

| Library | Commit |
|---|---|
| The Science Infrastructures (root) | `512477b` |
| heuresis | `5730e70` (with its decision 0005) |
| calx | `90aa710` |
| morphe | `14b3fb8` |
| physis | `e2b98d8` |
| tyche | `b4034b8` |
| domains/urban-science, domains/architecture | unchanged since `485f82c` |
| the shared problem | `../2026-11-05-a-district-economy-the-problem` (this change; its `economy.hpp`, `results/problem.json`, `results/sensors.json`) |

The C++ was built with g++ 11.5 and CMake 3.28 in WSL, in Release, without warnings from the study's own code. The repository's Python 3.12 ran `figures.py`.

```
cmake -S research/experiments/2026-11-05-a-district-economy-without-a-model -B ~/.cache/science-build/exp-district-economy-without -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-district-economy-without -j4
P=research/experiments/2026-11-05-a-district-economy-the-problem/results
~/.cache/science-build/exp-district-economy-without/without_a_model $P/problem-inputs.json $P/problem.json $P/sensors.json research/experiments/2026-11-05-a-district-economy-without-a-model/results --threads 4
.venv/Scripts/python.exe research/experiments/2026-11-05-a-district-economy-without-a-model/figures.py
```

The shared problem's `prepare.py` and `calibrate` must have been run first (its README §7); its `problem-inputs.json` and `district-geometry.json` are extracts of OpenStreetMap and are not committed. The programs are run from the repository's root. `--profile` times 200 random plans and stops.

## 4. Parameters

- **The value.** Every economic parameter is the shared problem's ([its §4.3](../2026-11-05-a-district-economy-the-problem/README.md)).
- **The study's own.** The kriging grid (§2.1); the radius of the empirical rule, 50 m, with 25 m and 100 m fitted beside it; the rule's clamp, [½, 2]; and the searches' settings (§2.3).
- **The runs.** 25 runs at the large budget (five searches × five seeds) and 30 at the small (six × five); the GA from the plan as built adds five to each. The whole took about half an hour on four threads, most of it CMA-ES (§7.3).

## 5. Data and Provenance

The study reads only the shared problem's prepared data ([its §9](../2026-11-05-a-district-economy-the-problem/README.md)).

| Data | Source | Licence | Retrieved |
|---|---|---|---|
| the counters' mean weeks and responses | City of Melbourne, Pedestrian Counting System, calibrated by the shared problem (`results/sensors.json`) | the City's open data, CC BY 4.0, on the understanding the data library records | 2026-10-30 |
| the trading establishments and their kinds | City of Melbourne, CLUE 2024 | CC BY 4.0 | 2026-11-04 |
| the district, its frontages and its network | OpenStreetMap, through the data library's city model and the designing study | ODbL 1.0 | 2026-10-31 |

The committed results hold the study's measures and plans, keyed by the shared problem's frontage indices (its `results/problem.json` keys them to OpenStreetMap's way ids and CLUE's base properties). They also hold the kriged passers at the 61 frontages and at the middles of the district's segments, keyed by segment id and with no geometry. The maps are drawn from the shared problem's uncommitted geometry.

## 6. Results

### 6.1 The People from the Data

**The kriging.**
- **The fit.** The likelihood chooses ℓ = 300 m and a nugget of half the variance (η = 0.5, s² = 1.19 in ln(1 + count)), on a grid that reaches 900 m and 0.6.
- **Its smoothness.** The field is smooth over the area and rough at each counter. At 300 m the whole area is one neighbourhood, and half of each counter's deviation from its hour's mean is its own.
- **Leave one out.**
  - The error is 0.94 in ln(1 + count) over all hours, with R² 0.73.
  - The median ratio of a counter's predicted week to its observed one is 0.99, with a quartile range of 0.57 to 1.48 and extremes of 0.27 and 10.5.
  - The predictions are drawn towards the mean (figure 62): the quietest counters are over-predicted and the busiest under-predicted.
- **The frontages.** Their passers run from 40,700 to 113,800 people a week, with a median of 72,200.

**The empirical rule.**
- **The elasticity.** It is 0.145 within 50 m (standard error 0.092, t = 1.58, errors clustered over the 55 counters), 0.129 within 25 m (t = 0.89) and 0.077 within 100 m (t = 1.05). None is significant at 5%, and the counters' hour-demeaned counts are barely explained by the open establishments about them (within-R² 0.015).
- **Its use.** It is used as estimated, 0.145, because the study's question is what a plan without a model does with the data it has. Under it a plan moves the district's people by about one per cent.

### 6.2 The Plan as Built

Under these flows and the shared value, the plan as built earns A$ 9,385 a week. A margin of 329,857 on spending of 599,661 is nearly spent on staffing (283,640) and outgoings (36,833).

| Measure | As built |
|---|---|
| people passing the frontages | 4.73 million passings a week |
| spending | A$ 599,661 a week: food and drink 268,716, retail 254,899, services 76,046 |
| spread | 0.417, 70.0 effective hours; the day's lines 0.572, the week's other lines 0.012 |
| evening share, weekend share | 0.147, 0.280 |
| mix H(K \| block) | 0.883 nats |
| access | 151.0 m: retail 160, food and drink 101, services 192 |

### 6.3 The Searches Compared

The best value found, A$ a week, five seeds each.

| Search | 400 evaluations: median (range) | 20,000 evaluations: median (range) | Time a run at 20,000 |
|---|---|---|---|
| GA | −80,542 (−87,390 to −79,234) | 143,120 (120,464 to 159,959) | 4 s |
| CMA-ES (separable, relaxation) | −94,277 (−95,974 to −90,205) | **160,115** (151,873 to 172,473) | 253 s |
| MAP-Elites | −100,963 (−102,904 to −75,955) | 63,563 (59,672 to 76,221) | 4 s |
| NSGA-II | −97,114 (−123,517 to −76,091) | 53,033 (36,758 to 77,788) | 4 s |
| random search | −121,531 (−134,218 to −118,277) | −109,083 (−110,098 to −99,313) | 5 s |
| BO (REMBO) | **−44,743** (−59,045 to −39,712) | — | 95 s for 400 |
| *GA from the plan as built* | *40,525 (32,311 to 45,305)* | *187,726 (180,039 to 192,275)* | 5 s |

- **Random search never reaches the plan as built.** Its best in 20,000 evaluations is a loss of A$ 99,313 a week. A random plan refits most frontages and opens them at random hours, and the plans worth more than the district's own are a small part of the space.
- **At 400 evaluations every cold search is below the plan as built.** BO on its ten-dimensional embedding is the least bad, at −44,743. The embedding's clipping sends most coordinates to 0 or 1: frontages closed or open all day, kinds at the corners. That is a coarse but cheap region of plans.
- **At 20,000 the vector searches lead.** CMA-ES on the relaxation leads (median 160,115), and the GA on the mixed encoding follows (143,120). The medians pass the plan as built at 3,400 evaluations (GA) and 3,800 (CMA-ES), MAP-Elites at 5,900 and NSGA-II at 6,800. CMA-ES was still climbing at the end and spread less over seeds (standard deviation 7,044 against 15,620).
- **The searches of other objectives keep less value.** MAP-Elites (63,563) and NSGA-II (53,033) spend their evaluations on the repertoire and the front, which is what they are for (§6.6, §6.7).
- **A GA started from the plan as built does best.** It passes the plan as built within its first generation and ends at 187,726. Knowing the district's own plan is worth more than any choice among the cold searches.

### 6.4 The Best Plan

The best plan found, by the GA from the plan as built (seed 4), earns A$ 192,275 a week (figure 64).

- **What it does.** It turns the district to retail: 39 frontages of retail (14 as built), 2 of food and drink (18), 1 of services (13), 2 office lobbies and 17 vacant. It changes 37 frontages and makes two links.
- **The value's parts.**
  - **Spending** is A$ 1,212,950 a week, 1,197,099 of it retail, and the margin on it 489,341.
  - **The costs** are staffing 232,015, outgoings 28,755, refits 35,044 and links 1,252.
  - **Where it trades.** Retail gathers at the north end, along Swanston Street from Little Collins to Bourke Street and along Bourke Street. That is where the kriged people are most numerous.
- **Why retail.** The shared value makes it so. A retail visit keeps A$ 0.108 of margin per passer against food and drink's 0.117, but on two staff instead of three, so it pays from about 650 passers an hour against 900. Retail's catalysis by retail is the strongest of the kinds' (CLQ 2.67, a = 0.98), so retail beside retail doubles its chance to stop, up to the clamp. Competition along a segment takes back less than catalysis gives. The cold searches find the same monoculture: CMA-ES's best (seed 2) has 28 retail and no food.
- **What it gives up.**
  - **The mix** falls from 0.883 to 0.302 nats.
  - **Access to food and drink** rises from 101 m to 300 m, and to services from 192 m to 366 m. Retail, at 44 m, is everywhere.
  - **The spread** is 0.409 (68.8 effective hours), against 0.417 as built.

The plan is the economics' assumptions taken to their end, not a finding about Melbourne. A study that wanted mixed streets would have to put the mix into the value or the constraints, as NSGA-II does (§6.7).

### 6.5 Time and Frequency

- **The hours** (figure 65). The best plan opens its retail through the day, mostly from 7–9 until 17–20, and at the weekend from the late morning. It shuts what the evening's thinner crowds do not pay for: retail's evening multiplier (0.7) costs more staff than it earns.
- **The spending.** It doubles at every hour of the day. Its profile follows the people's, with one crest at midday and a second at the evening commute on weekdays.
- **The week's lines** (figure 66). The best plan keeps the week's shape as built:

| Plan | Spread (mean line) | The day's lines | The week's other lines | Effective hours | Evening share | Weekend share |
|---|---|---|---|---|---|---|
| as built | 0.417 | 0.572 | 0.012 | 70.0 | 0.147 | 0.280 |
| the best value | 0.409 | 0.573 | 0.017 | 68.8 | 0.155 | 0.262 |
| the most spread of NSGA-II's front worth as much as built | 0.565 | 0.410 | 0.025 | 94.9 | 0.250 | 0.264 |

The most spread plan earns A$ 14,494 a week. It moves a sixth of the spending's power from the day's harmonics to the mean by opening into the evening: a quarter of its spending comes then, with four cultural frontages and its food and retail open late. The power left outside the day's lines (the weekday against the weekend) stays at 1–2.5% in every plan. The people's own week has little of it (the shared problem, §5.1), and a plan without a model cannot add people at the weekend who are not counted there.

### 6.6 The Repertoire of Evenings and Weekends

- **MAP-Elites** (figure 67, seed 1) fills 112 of its 144 cells, 101 of them with a plan worth more than nothing.
  - **The best cells** lie at an evening share of 0.25–0.40 and a weekend share of 0.20–0.35, at up to A$ 59,672 a week.
  - **The evening.** 50 cells with an evening share of 0.3 or more keep positive value, at most A$ 53,776.
  - **The weekend.** 36 cells with a weekend share of 0.4 or more keep at most 42,171.
  - **The extremes.** Spending almost wholly at the weekend, or none at all then, costs money.
- **Illumination.** Every evaluation of every run is binned into the same grid.

| Search | QD-score, A$k (median, range) | Cells reached |
|---|---|---|
| MAP-Elites | **3,245** (2,492 to 3,997) | 0.78 |
| CMA-ES | 1,393 (1,347 to 1,424) | 0.33 |
| GA | 943 (906 to 1,149) | 0.24 |
| GA from as built | 786 (744 to 886) | 0.03 |
| NSGA-II | 389 (143 to 588) | 0.29 |
| random search | 0 | 0.35 |

  - MAP-Elites illuminates the week's balances two to eight times better than the others, as it is meant to.
  - The GA from the plan as built never leaves the plan's own balance of the week: it reaches 5 of 144 cells.

### 6.7 Value against Spread and Mix

- **The fronts** (figure 68). These are the non-dominated plans in (value, spread, mix) among every evaluation of each run.
  - **Their hypervolume**, normalised over the union with the value's axis from the best found to the plan as built, ranks CMA-ES first: 0.498 (0.449 to 0.517). Then come the GA (0.312), NSGA-II (0.227), the GA from the plan as built (0.221) and MAP-Elites (0.161); random search scores 0.
  - **The two kinds of front.** The GA from the plan as built climbs in value while keeping the plan as built's week (spreads of 0.40–0.45). The cold searches reach spreads above 0.5.
  - **NSGA-II ranks below the GA and CMA-ES on its own three objectives.** Its fronts reach far in spread and mix, to spreads of 0.67 and mixes of 1.4. But of the 100 plans of its final population (seed 1) only 41 are worth the plan as built, and they lie below the value the single-objective searches reach. The value's range dominates the normalised volume. At 20,000 evaluations NSGA-II has not climbed the value's axis far enough for its spread to count.
- **What mix costs.** The most mixed of NSGA-II's final front worth as much as built has a mix of 1.261 nats: 14 retail, 5 food and drink, 4 services, 9 office lobbies, 4 culture and 25 vacant. It earns A$ 22,036 a week, a ninth of the best plan's. Its access, 167 m, is near the plan as built's.

### 6.8 Where the Plan Without a Model Is Blind

- **Links reroute people, and the value does not see it** (figure 69).
  - **On the network**, the people's shortest walks between frontages, weighted by their passers, would take a link in proportion to what it shortens. Athenaeum Place carried through to Collins Street would carry 14.0% of them, and Baptist Place to Collins Street 11.6%. A lane to Athenaeum Place would carry 7.5%, Rainbow Alley to Little Collins Street 5.7%, and Portland Lane to Russell Place 5.5%.
  - **The frontages affected.** Making the first would change the walks past 34 of the 61 frontages by a tenth or more, and the second past 44.
  - **In the value**, every link only costs: A$ 324 to 2,272 a week less on the best plan. The one exception is that a link brings catalysing frontages nearer, which never pays for its arcade.
  - **The best plan's two links** are made for that catalysis in the search's plan, not for the people.
  - **What a model would see.** A model with a network would reroute the walkers through the links, and their frontages' passers with them.
- **Nobody is conserved.** Under the best plan 1.5% more people pass the frontages than as built (4.80 against 4.73 million passings a week). The rule draws them from nowhere and takes them from no other street.
- **Lanes are as busy as streets.**
  - **The kriging.** It gives the one frontage on a lane (Rainbow Alley) 103,097 passers a week, as many as Bourke Street's (100,338). Nothing in it knows a lane from a street: the lanes take their neighbours' streets' numbers (figure 62).
  - **The counters.** The City's laneway counters count from 4,081 (Westwood Place, a lane of the centre) to 113,198 (the Flinders Street underpass), with a median of 55,700. None of them stands in the area.
- **Weeks that are not typical.** By the median responses of the 17 counters of the district's box:
  - **A public holiday on the Monday** (+6.0% people that day) raises the best plan's value to A$ 196,132 and the plan as built's to 11,915.
  - **A wet week** (1 mm of rain in each daylight hour, −15.2% people an hour) lowers them to 123,839 and −35,073.
  - **The plan as built** loses money in the rain under these assumptions. The best plan keeps two thirds of its value.
  - **No search was asked** for a plan that holds up in the rain.

## 7. Figures

The figures are drawn by `figures.py` from `results/` and from the shared problem's results (its geometry, uncommitted, for the maps). `gallery.md` holds them with their captions.

| Figure | What it shows |
|---|---|
| `search-62-people-without-a-model` | the kriged passers on the district's segments at four hours, the counters' measured values, the leave-one-out check and the profile likelihood |
| `search-63-searches-compared` | the best value against evaluations for every search at both budgets, five seeds |
| `search-64-best-plan-map` | the plan as built and the best plan on the map, each frontage sized by its spending; the value's parts; frontages by kind |
| `search-65-week-of-the-best-plan` | each frontage's hours through the week, as built and planned; the district's spending and people, hour by hour |
| `search-66-spending-lines` | the week's lines of the spending (Parseval), their shares, and the mean days, for three plans |
| `search-67-repertoire-evenings-weekends` | MAP-Elites' repertoire over the evening and weekend shares; the QD-score and cells reached of every search |
| `search-68-value-spread-mix` | the fronts in value, spread and mix of every search, coloured by mix; their hypervolumes |
| `search-69-where-it-is-blind` | the share of the people's shortest walks each link would carry, the walks past each frontage that one link would move, and what the value says of each link |

### 7.1 Results Files

| File | What it holds |
|---|---|
| `people.json` | the kriging's grid, choice and leave-one-out check; the elasticities; the passers as built at each frontage and hour; kriged passers at the segments' middles at four hours |
| `anytime.csv`, `runs.csv` | each run's best value along the way, and its best plan's measures and time |
| `fronts.csv`, `hypervolume.csv` | each large run's front of all it evaluated, and its hypervolume |
| `repertoire.json`, `repertoire-scores.csv` | MAP-Elites' archive (seed 1), and every run's coverage and QD-score in the same grid |
| `nsga-front.csv` | NSGA-II's final population's front (seed 1) |
| `plans.json` | the plan as built, the best plan, the most spread and the most mixed of NSGA-II's front, and each search's best (seed 1): their plans, outcomes with hourly profiles and week's lines, and per-frontage passers, spending and value |
| `blindness.json` | each link's share of the walks and effect on the frontages and the value; the people added; the lanes; the atypical weeks |
| `summary.json`, `run.log` | the measures of the plan as built and the best plan, every run's best and time, the normalisation |

### 7.2 Reproducibility

Every draw is keyed by its purpose (heuresis `search::Draws`). The study was run twice. The second run added the GA from the plan as built and corrected the selection of NSGA-II's plans. Every run common to both gave the same best values to the dollar.

### 7.3 Time

The whole took 1,944 s on four threads.
- **Evaluation.** A plan costs 0.18 ms on four threads once the distances of its links are cached, and about 2.7 ms when they are not (a district rebuilt and 61 shortest-path trees). Each of the 2,048 sets of links is built once.
- **CMA-ES** took 253 s a large run against 4–5 s for the others (§8.3).
- **Bayesian optimisation** took 95 s for 400 evaluations: its surrogate is refitted at every batch.

## 8. Caveats

### 8.1 What the Flows Leave Out

- **Kriging in the plane** knows nothing of the network. A lane, a street and a mall at one distance from the counters get the same people. The area holds only one frontage on a lane, but the passages and carried lanes would be lanes.
- **The hours are replicates.** The 168 hours are treated as independent replicates of one spatial field, which overstates the likelihood's information. The 168 fields are interpolated with one covariance, so their correlation in time is the counters'.
- **The passers are medians.** They are the posterior median exp(μ̂) − 1. The mean would be larger by a factor up to e^{s²/2} where the counters are far, about 1.8 here.
- **The rule is a correlation.** It is a cross-sectional association, not significant, read as a response. It conserves nobody and routes nobody.

### 8.2 What the Searches Leave Out

- **Cold starts are not how a planner works.** Every search but one starts from random plans, which is fair to the comparison and far from practice. The GA from the plan as built shows what that knowledge is worth.
- **The budgets are modest** for 316 decisions: 20,000 evaluations. CMA-ES was still climbing.
- **NSGA-II was not run longer.** Its hypervolume in value, spread and mix would rise if it were.
- **BO searched a ten-dimensional embedding,** not the space. In 621 dimensions a Gaussian process's lengths could not be fitted with its budget.
- **The value's assumptions decide the plan.** The retail monoculture is the shared problem's propensities, spends, margins, staffing and catalysis taken to their optimum (§6.4).

### 8.3 Gaps and Bugs Found

- **heuresis `strategies`: the separable CMA-ES is not O(n).** Its header says the separable variant costs O(n) a sample. But `CmaState::decompose` builds C^{−1/2} as a full n × n matrix by a triple loop, O(n³), even when B is the identity. And `whiten` and `colour` multiply full n × n matrices, O(n²) a sample. At n = 621 a large run took 253 s, 60 times the GA's. In the separable case B = I and C^{−1/2} = diag(1/D), so all three are O(n). The study did not edit the foundation.
- **The urban-science domain has no interpolation of the counts.** It reads each counter's tide at the counter. A frontage, a segment or a junction between counters has no flow of its own until a model or a kriging gives one. The study's kriging on tyche's Gaussian conditioning is its own (`flows.hpp`). It belongs with the domain's `flow` or `rhythm`, or with the planned `movement`, as the passers of a place.
- **tyche's `GaussianConditioning` takes one innovation vector.** A field observed at many times with one covariance (168 hours here) is conditioned once and each time's innovations whitened by `whiten`. That works, but `mean` and `log_likelihood` assume the stored innovations. A form with one factor and many right-hand sides would serve fields in space and time.
- **The shared problem's own gaps** are in its README (§11 and below).
