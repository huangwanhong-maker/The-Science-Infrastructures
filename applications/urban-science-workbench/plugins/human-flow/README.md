# Human Flow

*The district's walkers as agents, on the 3D map, moving with the clock.* A plug-in of the [Urban Science Workbench](../../README.md). Status: **draft** (2026-11-05). Activity: research. [ADR 0075](../../../../development/decisions/0075-agents-in-the-city.md); the model is the urban-science domain's `movement` module.

## 1. Purpose

The plug-in draws every walker of Melbourne's Hoddle Grid district as an agent on the map. Each walker has an origin and a target, chooses its way as it goes, walks at its own speed, slows in a crowd, waits at red signals and at full streets. The day runs in the engine; the map shows the walkers of the clock's hour, moving.

On that base it answers three questions:
- **Who is on the streets?** At the clock's hour and through the day, by pattern: passing through, arriving at work, leaving, walking between buildings, visiting shops.
- **What if a street changed?** A street, or one block of it, closed to walkers passing through or made one-way for a few hours, and how the walkers redistribute.
- **What does an event's crowd do?** A crowd brought to a street on the day of one of the City's event permits; the crowds and queues it makes, written as a record of simulated events; and the street's regime decided online as the crowd arrives, by POMCP on the walkers themselves.

Every result comes from a model with declared parameters. The walkers' departures are set by what the City's counters count each way, so the model's walkers pass each counter as often as the calibration expects; nothing the walkers do is drawn or written as observed.

## 2. How to Use It

1. **Open Melbourne** and, under *Human flow* in the rail, switch on **walkers on the streets**. The first time, the engine prepares the district (about ten seconds) and starts running the day; the legend says how far it has come. The walkers of the clock's hour appear once the run reaches it.
2. **Move the clock.** Each quarter hour has its own walkers; while the clock stands, they walk on for a minute and the minute plays round. Play the day to watch the morning's arrivals, the lunchtime crowd and the evening's departures.
3. **Colour them** by pattern, or by what they do: walking, waiting at a red signal, waiting to enter a full street, standing at an event. The street beneath each walker is coloured by how many walk on it per square metre.
4. **Pick a street** on the map for its walkers, its density and its level of service, and for *close it or make it one-way*.
5. **Compare with the counters:** switch on **counters: counted and simulated**. Sand columns are what each counter counts at that hour of the week; sea-blue columns are what the walkers make it read. Pick a counter for its day.
6. **Analyses:**
   - **Walkers of the district:** the walkers now and through the day; every counter's day, counted, expected and simulated; how the departures were set.
   - **A street closed or one-way:** choose the street (or the block you picked), its regime and its hours, and *Run the walkers*. The map then draws the walkers with the change, and each street by how many more (brick rose) or fewer (sea) walk it than as built. *all three regimes* puts open, one-way and closed side by side.
   - **An event's crowd:** choose one of the City's permits on the district's streets (2014–2018). Beside each is what the nearest counter counted that day against its usual day. Set the crowd's size, its hours and how long it stays, and the street's regime meanwhile, then *Bring the crowd*. The crowd walks in rose. *The record* tab lists the crowds and queues the walkers made, all simulated and drawn on the map in pale violet glass. *Its regime, online* asks POMCP for the street's regime every ten minutes while the crowd's true size is hidden from it.

A permit gives a date and a place, not hours or a crowd. You set those.

## 3. What It Stands On

| Part | What it is | Where |
|---|---|---|
| The walkers | Agents on nomos's engine, run synchronously (the same on any number of threads); their way chosen by the logit over every walk (θ = 0.02 per metre, widths by Hoddle's classes of street); Weidmann's speeds, Greenshields's slowing, capacity at a street's mouth, two-phase signals | `domains/urban-science/movement` |
| Their departures | The district-economy flow study's patterns (through, arrive, leave, local, visits), calibrated hour by hour of the week to the counters each way by tyche's NNLS; visits set at a tenth of the trips; each footpath's handedness fitted | `engine/district.hpp`; the study, `research/experiments/2026-11-05-a-district-economy-as-a-flow` |
| The counters | The City's counters in the district that count each way, over the 52 weeks from 27 October 2025, Victoria's public holidays left out | the data library's `street-sensors`, `street-counts` |
| The attractors | CLUE's 2024 establishments of the trading kinds, open by their kind's typical hours | `street-terrain-melbourne` |
| The permits | The City's event permits on the district's streets, from the record of the Hoddle Grid's events | `city-events-melbourne` |
| The regime decision | The movement module's `RegimeModel`, solved online by tyche's POMCP | `domains/urban-science/movement/regime.hpp` |

**Assumptions, stated where they enter.** Signals stand where two of Hoddle's main streets cross, with fixed, arbitrary offsets. A street's stubs where it leaves the district stay open under a closure, since walkers enter the district there. A visit and its walk back are two trips, as in the flow study. A crowd standing at an event counts in its street's density.

## 4. The Engine

`engine/main.cpp`, C++ on the movement module, built by the workbench's recipe. Its operations:

| Operation | What it answers |
|---|---|
| `prepare` | the walkway, signals, counters, entrances, streets, streams and their calibration |
| `run` | a run of a day's window started, or found: how far it has come |
| `frame` | a quarter hour's walkers every ten seconds for a minute, the density on each segment and the queues |
| `result` | a done run's totals, each hour's walkers and trips, each counter's reading beside its count and expectation |
| `compare` | two runs of one window, segment by segment |
| `permits` | the City's permits on the district's streets, and their days' counts |
| `advice` | POMCP's regime for an event's street, epoch by epoch, beside each regime held throughout |
| `record` | a run as a record of events, every one simulated |

Runs and advice are jobs on the engine's worker thread. The page asks again until they are done, and the map draws a run while it is still running.

```
python applications/urban-science-workbench/engine.py build human-flow
```

## 5. Checks

`tests/run-tests.mjs` checks, on Melbourne's district as the page sends it:
- the walkway and its counters;
- the departures calibrated each way;
- the walkers' counters reading what the calibrated flow expects;
- a quarter's frames inside the district;
- a street closed to walkers passing through;
- an event's crowd and its sound, simulated record;
- the permits;
- one decision of POMCP.

It also checks that the walkers' colours stay apart in every vision. The model itself is checked in the urban-science domain.
