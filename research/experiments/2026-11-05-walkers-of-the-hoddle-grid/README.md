# The Walkers of the Hoddle Grid

*Status: draft (2026-11-05). Research.* The first case study of urban science's walkers as agents (`domains/urban-science/movement`; [ADR 0071](../../../development/decisions/0071-agents-and-decisions-for-urban-planning.md), [ADR 0075](../../../development/decisions/0075-agents-in-the-city.md)), on the district of the shared problem [A District's Economy](../2026-11-05-a-district-economy-the-problem/README.md): the walking network of a few blocks of Melbourne's Hoddle Grid, between Elizabeth and Russell Streets and between Flinders Street and Lonsdale Street.

This is a study of a model on a real place, not a planning proposal. The walkers are a representation in Street Oceanography's sense: a mechanism offered for exclusion, its assumptions in the open (§4). Nothing here says how Little Bourke Street should be run.

## 1. The Question and Why

The user asked to "try ABM in urban science with several tasks, and visualize it". The district-economy studies read the people of the district as a flow: the logit over walks of [the flow study](../2026-11-05-a-district-economy-as-a-flow/README.md), calibrated hour by hour to the City's counters each way. A flow has no individuals. It holds no queue, no crowd slowing on a full street, no walker who waits at a red light, and no decision taken while the crowd is only partly seen. The questions:

1. **The week.** Do walkers as agents, each choosing its way by the flow study's law, reproduce that flow? And how do they compare with the counters, in space, in the hours of the week, and in frequency?
2. **A street closed.** Little Bourke Street was closed for the Lunar New Year of 2017 (the City's event permits, 27–29 January). Closed to walkers passing through, how do they redistribute, and at what cost?
3. **An event's crowd.** On Sunday 29 January 2017 the counters about Chinatown recorded a crowd. Fitted to that excess, how does the crowd arrive, dwell and leave, and what densities and queues does it make?
4. **A decision under partial observation.** If the street's regime (open, one-way, closed to walkers passing through) is decided every half hour from the counts at two counters alone, does deciding online on the walkers do better than holding one regime all day?

## 2. The Method

The model is the module `movement` (its theory note, [domains/urban-science/theory/movement.md](../../../domains/urban-science/theory/movement.md)). This study composes it with the flow study's calibrated trips; `walkers.cpp` is the program.

### 2.1 The Network

The network is the flow study's, as built: 369 junctions and 412 segments, 824 directed. Its widths come by class of street, and its parent lengths are kept, so that the route law is the flow study's exactly (`Bridge` in `walkers.cpp`). Its signals are assumed at the 47 junctions where three carriageways or more meet, since the data library holds no record of them (§8).

### 2.2 The Trips

The trips are the flow study's patterns, its own calibration run again by its own code: θ = 0.02 per metre, the homes' modes left out, and the visits set at ν = 0.1 of the trips. Its fit is reproduced exactly: gross R² 0.7214 (0.5293 with each counter left out), net R² 0.5494. Its amplitudes agree with its `results/calibration.json` to 4.3 × 10⁻⁴, which is that file's rounding to four figures.

Each pattern becomes a stream of the module: a fixed law of (origin, target) pairs times an hourly rate. The week has 666 streams:
- arrivals at work and departures, from each entrance;
- through trips, for each pair of entrances;
- local trips;
- visits to each kind of establishment;
- visits to each of the 61 frontages, as built.

Departures are Poisson within each minute.

### 2.3 The Walkers

Each walker is an agent of nomos's engine on a synchronous schedule, its draws keyed by its identity:
- it chooses its way junction by junction by the logit over walks, as a Markov chain on directed segments;
- it walks at its own speed (Weidmann's 1.34 ± 0.26 m/s), slowed by the density in the 10 m ahead of it (Greenshields, jam at 5.4 m⁻²);
- it waits for its phase at signals (a cycle of 120 s, half and half);
- it is admitted to a segment by Godunov's supply, so it queues at full ones;
- it may dwell, then walk back.

The step is 2 s. A run tallies:
- entries into each directed segment, by the hour of entry and by the hour of departure;
- each counter's passes each way;
- the person-seconds on each segment, all and at Fruin's levels E and F (1.08 and 2.0 m⁻²);
- the queues;
- the trips' times and delays.

### 2.4 Micro against Macro

Counted by the hour of departure, the walkers' entries have exactly the flow's mean in every hour (the theory note, Corollary 3.5). The program computes the flow twice:
- by the flow study's own `mode_flows`, on every segment;
- by the module's `expected_flow`.

It then compares the walkers with it state-hour by state-hour.

### 2.5 The Street, the Festival and the Regimes

- **The street.** Little Bourke Street from Swanston Street past Russell Street to the edge of the district: 15 segments, 251 m. The permit of 27–29 January 2017 names the street.
- **The street's own walkers.** On the Sunday, a pair of the trips walks as the street's own class (1) when it starts or ends on the street, enters the district at the street's east end, or has no walk to its target with the street closed (its lanes are reached only through it). Closed to walkers passing through, the street still admits them.
- **The festival's visitors.** Two patterns: from the district's entrances (by their flow) and from its buildings (by their floor). Each visitor is bound for a junction along the street, chosen by the length about it, dwells there, and walks back.
- **Fitting the crowd.** Their arrivals by hour (08:00 to 22:00, each pattern) and their mean dwell are fitted by NNLS to the excess of Sunday 29 January 2017 over the six Sundays about it, at the counters whose counts of those days are whole. A visitor's walk back falls in the hours after its arrival by the law of an exponential dwell.
- **The festival day.** The mean Sunday's walkers and the fitted visitors together, the street given to walkers at 8 m (its 10.1 m reserve less 2 m of stalls, an assumption).

**The regime decision.**
- **Policies.** Every half hour from 09:00 to 17:00 (16 decisions), the street's regime is chosen:
  - by POMCP planning on the walkers themselves;
  - or held all day open, one-way from Swanston Street, or closed to walkers passing through.
- **Hidden.** The crowd's scale, 0.6, 1 or 1.6 times the fitted (prior ¼, ½, ¼), and the configuration of every walker on the links no counter sees.
- **Observed.** The half hour's readings at two counters: Little Bourke Street at Swanston Street (19) and Russell Street (21). Each is binned at the terciles of an open festival day's.
- **Reward.** −(person-hours present) − 10 × (person-hours at Fruin's level F); the weight is an assumption.
- **The planner.** It runs lightened copies of the walkers at 4 s steps. Each decision gets 12 simulations of depth 2 with a leaf value of zero, and the belief is 6 particles refilled by rejection on the observed bins (at most 12 attempts). If the belief is lost, the planner starts again from the world's configuration, and each such restart is counted.
- **Comparison.** Each scale and each of 2 days is run under all four policies with common random numbers.

## 3. Code Version and Running

| Library | Commit |
|---|---|
| The Science Infrastructures (root) | `c4adfdf`, with urban-science's `movement` written in this work (not yet committed) |
| calx | `90aa710` |
| morphe | `14b3fb8` |
| physis | `e2b98d8` |
| tyche | `b4034b8` |
| nomos | `5382beb` |
| heuresis | `e878c3c` |
| the shared problem and the flow study | as at `c4adfdf` |

Built with g++ 11.5 and CMake 3.28 in WSL, in Release, without warnings. From the repository's root:

```
cmake -S research/experiments/2026-11-05-walkers-of-the-hoddle-grid -B ~/.cache/science-build/exp-walkers -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-walkers -j2
P=research/experiments/2026-11-05-a-district-economy-the-problem/results
~/.cache/science-build/exp-walkers/walkers $P/problem-inputs.json $P/problem.json $P/sensors.json research/experiments/2026-11-05-walkers-of-the-hoddle-grid/results --threads 4 --task all
.venv/Scripts/python.exe research/experiments/2026-11-05-walkers-of-the-hoddle-grid/figures.py
```

`--task week` runs the week alone; `--task sunday` the closure, the festival and the regimes; `--task regime` the regimes alone; `--task diag` compares single trips with the flow study's `Trips::visits`.

The regime decision is resumable. Each case is appended to `results/regime-cases.jsonl` when it ends, and a case found there is not run again; `--regime-cases N` runs at most N of the rest. A case depends only on its scale, its day and its policy (its index seeds POMCP), so a resumed run gives what one whole run would. The machine's WSL restarted twice during this work. The first 17 cases were then recovered from the lost run's log (`results/run-sunday.log`, by `recover_journal.py`), with their values, totals and regimes but not their rewards epoch by epoch. Those records carry `"from_log": true`.

## 4. Parameters

| Parameter | Value | Source |
|---|---|---|
| θ, the route's sensitivity to length; widths; w_ref | 0.02 per m; by class; 12 m | the flow study (calibrated; assumptions) |
| trips and their amplitudes | the flow study's patterns, fitted to the counters each way | the flow study's calibration, rerun |
| free speeds | normal, 1.34 ± 0.26 m/s, clipped to 0.5–2.5 | Weidmann 1993 |
| jam density; the speed's law | 5.4 m⁻²; Greenshields, 1 − ρ/ρ_j | Weidmann 1993; as urban::flow |
| the window ahead; a crowd's creep | 10 m; 0.05 m/s | assumptions |
| a standing queue's density | 2.5 m⁻² | assumption, after Fruin's levels for queues |
| signals | at junctions of 3+ carriageways; 120 s, half and half; fixed offsets | assumptions (no data) |
| Fruin's levels E and F | 1.08, 2.0 m⁻² | Fruin 1971 |
| the step | 2 s (the planner's copies 4 s) | choice |
| the festival street's walkable width | 8 m | assumption |
| the ordinary Sundays of 2017 | 8, 15, 22 January; 5, 12, 19 February | choice: the six about the festival |
| the festival's dwell | chosen from 15 min to 3 h by the fit | fitted |
| the crowd's scales and prior | 0.6, 1, 1.6; ¼, ½, ¼ | assumption |
| the reward's weight of a person-hour at level F | 10 | assumption |
| POMCP | 12 simulations a decision, depth 2, leaf value 0, 6 particles, 12 refills, c = 300 | choice, bounded by time on a shared machine |

## 5. Data and Provenance

| Data | Source | Licence | Retrieved |
|---|---|---|---|
| the walking network, buildings | OpenStreetMap through the data library's city model, via the shared problem | ODbL 1.0, © OpenStreetMap contributors | 2026-10-31 |
| hourly pedestrian counts each way, 2009–2026; the counters' places | City of Melbourne, Pedestrian Counting System (the data library's `street-counts`, `street-sensors`) | CC BY 4.0 | 2026-10-30 |
| establishments; floors | City of Melbourne, CLUE 2024 (the shared problem's inputs) | CC BY 4.0 | 2026-11-04 |
| event permits, 2014–2018 | City of Melbourne (the data library's `melbourne-event-permits`) | CC BY 4.0 | 2026-11-04 |
| public holidays | Nager.Date (the data library's `street-holidays-victoria`) | Nager.Date's API terms | 2026-10-30 |

No data were added to the library. `results/maps.json` holds the segments' geometry and the walkers' places on it, an extract of OpenStreetMap. It is not committed, and `walkers` makes it again. Everything else in `results/` is committed.

## 6. Results

### 6.1 The Week (task 1)

**Micro against macro.**
- **The two references.** The module's expected flow equals the flow study's own `mode_flows` to 9.3 × 10⁻¹⁴ of the largest flow (4,977 walker-segments an hour), over every directed segment and hour of the week: 53,580,571 walker-segments.
- **The walkers.** One simulated week holds 3,160,120 trips, run in 203 s on two threads (the same totals as on four): 1,554,585 local, 1,190,395 through, 279,382 visits, 126,932 arrivals at work and 8,826 departures. Every walker arrived and none failed. At most 7,295 walked at once (Saturday afternoon).
- **Counted by departure.** Their entries into segments total 53,610,965. Over the 118,722 segment-hours whose expectation is at least 5, the standardised deviations (n − μ)/√μ have mean 0.005 and variance 1.027, with 0.29% beyond ±3 (the normal's is 0.27%).
- **The variance.** It exceeds 1 by what Corollary 3.5 predicts: a walker can enter a segment more than once, so a segment's entries are compound Poisson.
- **Counted by the hour of entry.** The walking takes time. The entries then lag the quasi-static flow by up to 28,000 an hour (16% of that hour's), at Wednesday morning's rise. Over the week the two agree with R² 0.9988.

**Against the counters.**
- **Overall.** Over the 16 counters × 168 hours, the walkers' readings have R² 0.7202 against the counts: the flow's own fit (0.7214), less the noise of one simulated week. Against the flow they have R² 0.9983.
- **Counter by counter.** Elizabeth Street at 197 (0.98) and 124 (0.93) and the Town Hall on Swanston Street (0.95) are walked as counted. Flinders Street at ACMI (0.03), Royal Lane at 231 Bourke Street (0.03) and Flinders Lane at Degraves Street's south side (0.05) are not, as in the flow, where its patterns of trips do not reach.

**Time.**
- **The trips.** A trip takes 549.5 s on average. Of that, 84.6 s is delay over walking the same way alone at its free speed: the signals and the crowd.
- **By the hour.** The delay ranges from 54 s to 108 s over the hours with more than a thousand trips. In the busiest hour (Monday 13:00, 49,433 trips) it is 83.5 s.
- **Waiting.** 96,031 of the week's 482,385 person-hours are spent waiting, at red or for room.
- **Crowding.** 540 person-hours are at Fruin's level E or worse and 24 at F, counting walkers walking (a queue's density is assumed, not read). Of them, 441 are on Little Collins Street, 4 m wide, mostly on two of its segments, and 63 on Bourke Street.

**Frequency.**
- **The spectrum.** The spectrum of the walkers on the streets, minute by minute over the week, holds the counted week's structure: the day D1 strongest, its harmonics D2 to D6, the week W1, and the day's sidebands D1 ± W1 (the weekend's difference).
- **The daily line's amplitude.** At the counters it follows the counts where the flow does, and departs where the flow does.
- **The crest.** The walkers' crests all fall between 14:50 and 15:34, where the counted crests range from 13:32 to 17:29. Like the flow, the walkers have one tide for the whole district. The flow study found this of its own model (its §6.4), and walking does not change it: it comes from the trips' patterns, whose spatial laws are fixed.

### 6.2 A Street Closed (task 2)

The mean Sunday has 418,801 trips. In the logit flow, the street (Little Bourke Street from Swanston Street past Russell Street to the district's edge, 15 segments, 251 m) carries 123,009 walker-passes.

- **What moves.** Closed to walkers passing through, the exact flow moves only 1,752 passes onto other streets. Most of the street's walkers are its own, in three ways:
  - they start or end on it, at its buildings and establishments;
  - they enter the district at its east end, where it leaves the box;
  - they can reach their lanes only through it.

  The walkers' change, segment by segment, lies on the flow's.
- **What it costs.** The day costs the walkers 189 person-hours *fewer* (63,351 against 63,540), and the mean trip shortens from 546.2 s to 544.6 s.
- **Why it costs nothing.** This is a property of the logit over walks, not a saving a closure would make. The law weights every walk, and some of them wander through the street's lanes before going on. Closing the street removes those walks, so the expected walk to a target is shorter. Its largest drops are on the street itself (−21,654 passes on one segment) and on Russell Street beside it (−8,186), and they are not made up elsewhere.
- **The reading.** A street that mostly serves its own walkers can be closed to the rest at no cost to them, by this model. The model's walkers do not wander for pleasure, but its law lets some walks do so.

### 6.3 The Festival's Crowd (task 3)

**What the counters saw.** On Sunday 29 January 2017, five counters inside the district have whole counts that day and on the six Sundays about it: the Bourke Street Mall (2), Chinatown at Swanston Street (19), Chinatown on Little Bourke Street east of Russell Street (20), Russell Street (21) and Collins Street (53). Their excess from 08:00 to midnight was:

| Counter | Excess counted (people) | The walkers' excess |
|---|---|---|
| Chinatown, Swanston Street (19) | 6,604 | 9,452 |
| Chinatown, Little Bourke Street (20) | 4,790 | 1,470 |
| Russell Street (21) | 12,395 | 4,082 |
| Bourke Street Mall (2) | 1,700 | 4,792 |
| Collins Street (53) | 433 | 1,320 |

**The fit.**
- **The visitors.** 16,341 of them: 8,245 from the entrances and 8,096 from the buildings. Their arrivals peak at 2,601 between 11:00 and 12:00.
- **The dwell.** The fit's R² against the excess is 0.295 at a mean dwell of 15 minutes, the grid's shortest. It falls only to 0.286 at three hours: the counters cannot tell the dwell, since a visitor's walk back passes them as its walk out did.
- **What the fit misplaces.** The model cannot put the crowd on Russell Street and east of Little Bourke Street while keeping it off the Mall and Swanston Street. Its visitors' spatial laws (the flow study's origins, the logit's ways) bring them from the west. The counts say the crowd came along Russell Street and from Chinatown's eastern blocks, beyond the district's edge.

**The day.**
- **The walkers.** An ordinary Sunday's walkers and the fitted visitors together make 435,332 trips. They spend 73,505 person-hours on the streets, against 63,540 on the ordinary Sunday. The mean trip is 574 s against 546 s.
- **The street.** At most 799 walkers are on it at once, at 13:32. Its busiest segment holds at most 0.63 m⁻² over a minute, Fruin's level C.
- **Crowding.** Across the district, 85 person-hours are spent at level E or worse (45 on the ordinary Sunday) and 2.2 at F (0.7). The crowds of level E stand on Little Collins Street, 4 m wide, not on the festival's street, which is given 8 m.
- **Queues.** The longest queue at the festival street's mouths is 31 walkers; the longest at the signals of Bourke Street is 99. At its fitted size, the festival is a busy street, not a crush.

**The record.** The day, written in urban science's record, holds 1,658 events: 80 crowds, 1,178 queues and 400 hourly counts. All are simulated, and the record is sound by `events::problems`.

### 6.4 The Regime Decided Online (task 4)

**What was compared.** Three sizes of the crowd, two days each and four policies make 24 cases, with common random numbers within a day. A case's value is its reward summed over the 16 half hours from 09:00 to 17:00. A half hour's reading at each sensed counter is binned at 788 and 973 walkers. Holding the street open is worth −42,507 and −42,663 on the two days at 0.6×, −45,747 and −45,910 at 1×, and −50,638 and −50,764 at 1.6×. Against it:

| The crowd | One-way all day | Closed all day | POMCP |
|---|---|---|---|
| 0.6× its fitted size | −446, −436 | +113, +107 | −165, −132 |
| 1× | −711, −722 | +97, +106 | −136, −218 |
| 1.6× | −1,165, −1,160 | +99, +87 | −762, −462 (8 and 1 restarts) |

- **Closed is best at every size.** It gains 87 to 113 a day, and its gain does not grow with the crowd. It is the saving of §6.2: the walks that wandered through the street's lanes are removed. It is not relief of a crush. The person-hours at level F are about the same closed as open: 1.49 and 1.72 against 1.46 and 1.81 at 0.6×; 4.5 and 4.8 against 4.35 and 4.5 at 1.6×.
- **One-way is worst, and worse as the crowd grows.** It loses 436 to 446 at 0.6× and 1,160 to 1,165 at 1.6×. Every walker bound west goes round by the streets beside it, the visitors walking back among them. Its person-hours at level F are two to four times open's (17.6 against 4.35 at 1.6×).
- **The decision does not depend on what is hidden.** In these runs closed is best at every size of the crowd, on both days, so no reading of the counters can change the best regime. The value of information is zero, and the best policy is to close the street at 09:00 and keep it closed. For this street, the POMDP is well posed but degenerate.
- **POMCP does worse than every fixed regime but one-way.** It loses 132 to 218 against open at 0.6× and 1×, and 462 and 762 at 1.6×.
  - **Its choices are near uniform.** Of its 96 decisions, 36 are open, 35 one-way and 25 closed.
  - **Why.** It has 12 simulations a decision, about four for each regime. Each draws its half hour of walkers afresh, from a particle whose crowd may be any size the belief still holds.
    - **The noise.** Between two days of the same crowd (1.6×, closed), a half hour's return spreads by 23 person-hours (one standard deviation). Between sizes of the crowd it differs by 203 (0.6× to 1×) and 306 (1× to 1.6×) a half hour.
    - **The differences.** Under common random numbers (1.6×, day 1), closed gains only −2 to 11 a half hour over open (5.5 on average), and one-way loses up to 121.
    - **The result.** The planner cannot tell the regimes apart through that noise, and it pays for the half hours of one-way it chooses.
- **The belief was lost at the largest crowd.** At 1.6×, the particles failed to explain the readings 8 times on day 0 and once on day 1. Each time, the planner started again from the world's configuration (§2.5). At 0.6× and 1× its belief was never lost.
- **Cost.** A day of POMCP took 375 to 1,724 s on one thread, 7,291 s for the six. A fixed regime's day took 35 to 141 s.

**The reading.** For this street the decision is not hard: the closure is best whatever the crowd, by a margin the route law sets. A planner that compared the actions at a node through the same draws (common random numbers, §8.2) would see the differences that 12 independent simulations cannot. Where the best regime does depend on what is hidden, in the module's exactly solvable corridor, POMCP on agents finds HSVI's value: −14.32 ± 0.42 against −14.26.

## 7. Figures

`figures.py` draws nine figures from `results/` (and `results/maps.json`, which `walkers` makes and which is not committed). `gallery.md` holds them with their captions. The place is drawn in the Urban Science Workbench's colours: a gentle sea blue for the walkers and sand for the shore. The series' colours are deeper steps of both, with plum and green, checked for colour vision deficiencies by the dataviz validator.

| Figure | What it shows |
|---|---|
| `agents-01-the-walkers-of-a-week` | the district's segments by the week's walkers; every walker at Wednesday 08:30 and Saturday 13:30, walking or waiting |
| `agents-02-the-walkers-and-the-flow` | micro against macro: each segment's week, the standardised deviations of every segment-hour, the lag of walking hour by hour |
| `agents-03-the-counters-week` | six counters' weeks: counted, the logit flow, the walkers |
| `agents-04-the-lines-of-the-week` | the spectrum of the walkers on the streets minute by minute; the daily line's amplitude and crest at the counters |
| `agents-05-a-street-closed` | where the walkers go when Little Bourke Street is closed to walkers passing through; the walkers' change against the flow's; the cost hour by hour |
| `agents-06-the-festival-crowd` | the counters' excess on Sunday 29 January 2017 against the walkers'; the visitors' arrivals and the dwell fitted |
| `agents-07-densities-and-queues` | the crowd's density along the street minute by minute; walkers on it and at Fruin's level E; the longest queues |
| `agents-08-the-regime-decided` | each policy's value against holding the street open, for each scale of the crowd; the regimes POMCP chose |
| `agents-09-the-simulated-record` | the festival day as a record of simulated events: crowds and queues on the street and junction spaces |

## 8. Caveats

### 8.1 What the Walkers Inherit, and What They Assume

- **The trips are the flow study's.** The walkers reproduce its fit to the counters, and its misfits: Russell Street, Flinders Street at ACMI, Royal Lane. They also keep its one tide for the whole district. A walker's choice of way is the flow's, so walking adds time, queues and crowds, not new destinations.
- **The route law reads no density.** A walker chooses its way by the logit over walks whatever the crowd ahead. It is slowed and held by the crowd but does not avoid it. This is what makes the walkers' mean exactly the flow's (§2.4). A walker that rerouted round a jam would be another model, and would break that identity.
- **The kinematics are a one-dimensional local reading of the kinematic wave.** There are no lanes, no overtaking, no friction of counterflow beyond density, and no social force. The window ahead (10 m), the creep (0.05 m/s) and a standing queue's density (2.5 m⁻²) are assumptions.
- **The signals are assumed.** They stand at every junction of three carriageways or more, with a cycle of 120 s and fixed offsets. The data library holds no record of them (§8.2).
- **The festival's visitors are a sketch.** Two patterns, entrances and buildings, each with one spatial law, are fitted to the excess at five counters, the only ones whose counts of those days are whole. The festival of 2017 is read through the trips of 2025–26. The street's 8 m of walkable width is an assumption.
- **The regime decision is small.** POMCP has 12 simulations a decision on lightened copies of the walkers at 4 s steps, so its values are estimates, not bounds. Each policy has two days for each scale of the crowd, so differences between policies are read against the days' spread, with common random numbers. The reward's weight of a person-hour at level F (10) is an assumption. The horizon is short: depth 2 with a leaf value of zero sees an hour ahead. Seventeen of the 24 cases are kept from the log of a run lost when WSL restarted (§3): their values, totals and regimes, without their rewards epoch by epoch.

### 8.2 Gaps and Bugs Found

**In the module, found by this study and fixed.** Each was found by a check that failed, and each now has a check.
- **The walk sums' tolerance.** The first stopping rule compared each sweep's change with the largest value. On a trip whose total weight of walks is 2.7 × 10⁻¹³, it left the visits converged to 5 × 10⁻⁴, and the expected flow differed from the flow study's by 1.8 × 10⁻⁴. Each state is now held to its own relative tolerance. Single trips agree with the flow study's `Trips::visits` to 1.4 × 10⁻¹³, the week's flow to 9.3 × 10⁻¹⁴.
- **The density of a segment taken whole.** Fed at capacity, a narrow lane filled past the critical density and passed 2,591 an hour against its capacity of 6,512. Each walker now reads the density ahead of it, and admission is Godunov's supply. The lane passes 6,409.
- **A queue at one point.** Walkers waiting at red stood all at the segment's end and read as an infinite density. The district locked by mid-afternoon of the first week run, with 115,214 walkers held. Waiting walkers now stand as a queue behind the end.
- **Dwellers in the stream.** A festival's visitors dwelling counted as walkers ahead and as room taken at the mouths, and they locked the festival day. They now stand aside: they slow the stream through their density but do not stop it.
- **Queues counted at Fruin's walkway levels.** A queue stands at its assumed density, so every walker in one read as level F. Only walkers walking or dwelling are now tallied.
- **A counter at a segment's very end.** Such a counter missed one way. Two of the district's counters lost a way entirely. Counting lines are now kept half a metre inside.

**For other units of the infrastructure.**
- **A general method in a domain.** The logit over walks as a Markov chain (`route.hpp`: the walk sums, the walkers' chain, the expected flow) names no city. It belongs in a foundation, beside physis's `first_passage` or nomos's `choice`. It is written in the domain because this work may not edit the foundations.
- **nomos `agents`.** `tick` keys each agent's draws by its index, not its identity. Walkers are born and die, so the module keys its draws by identity inside its rule (nomos's open question 3).
- **tyche `decision`.** POMCP cannot be handed a belief from outside. A planner whose particles are lost must be rebuilt, losing its tree. Its refill by rejection runs the generative model once per attempt, which for an agent-based model is the costly part. A weighted particle filter (as POMCPOW's) would serve such models better. And each simulation draws its walkers afresh, so the actions at a node are compared through independent noise. Here that noise (23 person-hours a half hour between draws of one crowd, 203 to 306 between the crowds a belief holds) is larger than the differences between the regimes (−2 to 11 a half hour for closed against open, under common random numbers; §6.4). Common random numbers across the actions at a node, the scenarios of PEGASUS (Ng and Jordan 2000) and of DESPOT, would make such differences visible to a planner with a small budget.
- **urban-science `events`.** `record.hpp`'s comment lists provenance as recorded, inferred, modelled or synthetic. ADR 0071 and this module write *simulated*.
- **The data library.**
  - It holds no traffic signals, no footpath steepness and no pedestrian network of the City's own. Manifests for Victoria's traffic signals and the City's `footpath-steepness` and `pedestrian-network` would replace the assumptions.
  - Of the 16 counters inside the district, only 5 have whole counts on the festival's Sunday and the six Sundays about it.
- **The study's own bug.** Building the streams, the pair whose origin is its own target was dropped but the others' law renormalised, inflating the arrivals by up to the dropped share. The module's expected flow then differed from the flow study's by 0.3% until it was found.
