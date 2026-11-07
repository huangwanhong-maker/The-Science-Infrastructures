# Urban Science Workbench

*A three-dimensional GIS from the whole Earth to a doorstep, any city standing on it, with the analysis of streets and buildings on top, built of plug-ins.* Status: **draft** (2026-11-07; [ADR 0078](../../development/decisions/0078-the-whole-earth-in-the-urban-workbench.md)).

**Since ADR 0078.** One view runs from the globe to a street. Any place opens as one comes down to it, from the world's archives: fly to a name in the gazetteer (click the place's title) and its buildings rise on the real terrain. The page puts a question first, by five lenses (People, Movement, Fabric, Risk, Change), each layer showing the scales at which it means something. It is pink and purple, on a plum night or lilac paper. Three plug-ins join: *The territory* (the world's people, poverty, lights, settlements and hazards on one grid of hexagons), *Flood risk* (Freetown's assessment, and hazards by name anywhere) and *Reach and use* (a wave of minutes down the streets from any place, and how a city's walkers use its places). What was taken from the list of tools and data the user gave, what was left and why, is in the tool *Sources and what stands behind them*.
 Activities: research, learning. [ADR 0060](../../development/decisions/0060-streets-as-open-systems.md), [ADR 0065](../../development/decisions/0065-new-york-in-the-workbench.md), [ADR 0067](../../development/decisions/0067-event-spacetime.md), [ADR 0075](../../development/decisions/0075-agents-in-the-city.md). Named the Urban Building Workbench until the user renamed it (ADR 0067).

## 1. Purpose

The workbench is a GIS first. A city stands in three dimensions — its buildings at their mapped heights, its streets at their widths, its parks, water and squares — in an exact local frame on WGS 84. A clock runs in the city's own time zone and lights the city by the real Sun of the instant. On that base, plug-ins put layers on the map and analyses beside it. Two cities open, each where its records are:
- **Melbourne's centre,** from Southbank to Carlton, where the City's hundred pedestrian sensors count every hour;
- **Lower Manhattan,** from the Battery to Canal Street and from the Hudson to the Brooklyn Bridge ([ADR 0065](../../development/decisions/0065-new-york-in-the-workbench.md)): its buildings at the roof heights the city's records give, its trees and its shore from Valdrada's compiled city of New York, its streets from OpenStreetMap, and the Brooklyn Bridge's promenade, where the city counted its walkers each way every hour from 2017 to 2019.

The first plug-in, **Street Oceanography**, reads streets as open systems, and draws them on the map as a sea:
- the tide of a street's walkers, rising and falling with the clock, and its surge on a date;
- its currents, walkers each way at their true numbers, and the commuters' flood and ebb;
- the waves of a crowd at the crossings of a street (a model);
- how a street answers the weather and the calendar;
- how its uses mix and work together, a store's visitors spreading toward another on the map.

It reads real records of Melbourne, and of New York's bridge. Every result says whether it was measured from records, proved for any street, or comes from a model with declared parameters, and gives how it was obtained and its limits.

The second plug-in, **Landscape Dynamics** ([ADR 0062](../../development/decisions/0062-landscapes-as-cellular-configurations.md), [ADR 0063](../../development/decisions/0063-the-district-as-a-complex-to-work-on.md)), takes a district (the Hoddle Grid's middle; New York's Financial District) as cells — junctions, streets, blocks, and its buildings as 3-cells sharing their walls — to be seen as a complex, set and rewritten by hand, and run as one process in which structure, use and meaning change one another:
- walkers each way along the streets and the warmth of the blocks, between events;
- two stores of the district opening and shutting; a rule coming into force, carried out as a lane that rewrites the cells, and the lane taking on a meaning as it is used;
- the day drawn as the district's spacetime, time rising above the map: the streets' traffic slice by slice, the stores' worldlines, the lane's sheet beginning in mid-air and deepening once recognised, the events joined by their causal order;
- the same events solved exactly: what coarser descriptions forget, and how interventions interact;
- a building picked on the map filled with rooms, its heat through the day under the Sun;
- the district's own frequencies, its flows round the holes, and what the streets remember of the shops.

The third plug-in, **Event Spacetime** ([ADR 0067](../../development/decisions/0067-event-spacetime.md); its design, [development/ideas/2026-11-04-event-based-gis-on-a-2-complex.md](../../development/ideas/2026-11-04-event-based-gis-on-a-2-complex.md)), draws a district's record of real events as a 2-complex in spacetime. The district is the Hoddle Grid's middle, and the record runs from 2002 to 2024.
- **The complex.**
  - Its volumes rise as worldlines: buildings, street spaces and junction spaces.
  - The interfaces between them stand as sheets between their worldlines, cut into lenses at the events that meet them: party walls, frontages, and the ends of streets.
  - The events are vertices at their place and date: buildings built, refurbished or completed, places registered as heritage, street events under permits, and the censuses as ticks.
- **Its two forms.** Seen as its refinement, each event is split into its volumes. Seen as the event complex, each event is one vertex, its worldlines bending into it.
- **The moment.** A plane at the clock's date carries the network of what stood then.
- **The record and its complex.** The analysis gives what the record holds and what it left out, and the cells and homology of both complexes.

The fourth plug-in, **Human flow** ([ADR 0075](../../development/decisions/0075-agents-in-the-city.md); [its README](plugins/human-flow/README.md)), draws the Hoddle Grid district's walkers as agents, the urban-science domain's movement model, on the map:
- **Every walker, moving with the clock.** Each has an origin and a target, chooses its way as it goes, slows in a crowd and waits at red signals and at full streets. They are coloured by pattern (passing through, arriving at work, leaving, between buildings, visiting shops) or by what they are doing, over the streets coloured by their density.
- **Set by the counters.** Their departures are the district-economy flow study's patterns, calibrated hour by hour to what the City's counters count each way. Each counter's reading by the walkers stands beside its count.
- **A street changed.** A street or a block closed to walkers passing through, or made one-way, for a few hours: the walkers run again, and the map draws how they redistribute.
- **An event's crowd.** A crowd is brought to a street on the day of one of the City's event permits, and the crowds and queues it makes are written as a record of events, every one simulated. The street's regime can be decided online by POMCP on the walkers, the crowd's size hidden from it.

The city can be seen through, for this and every plug-in: the rail's slider lowers its opacity.

## 2. The Map

| Part | What it draws | Source |
|---|---|---|
| Buildings | Extruded footprints with their holes, at their mapped height (from `height`, else levels × 3.2 m); paler where no height is mapped and a default stands. They can stand as built, lowered to a fifth of their height, or of glass, so that the streets show between towers; and the whole city (buildings, ways, spaces, trees) can be seen through, its opacity lowered by the rail's slider, so that what a plug-in draws above or within it shows | Melbourne: OpenStreetMap; New York: Valdrada's compiled city, the roof heights of NYC's Building Footprints (data library: `city-model.json`, a file an area) |
| Streets and ways | Ribbons at their mapped or estimated widths; pedestrian ways tinted | OpenStreetMap |
| Parks, water, squares | Flat areas, with their holes (piers in the river) | OpenStreetMap; New York's water traced on Valdrada's ground |
| Trees | A crown on its trunk, its breadth by the trunk's diameter (New York's) | Valdrada: NYC's street-tree census and NYC Parks' tree points |
| The clock | A date and an hour in the city's time zone (daylight saving its own), playing through the day; the Sun's altitude and azimuth there, lighting the city and casting shadows | the browser's time-zone database; the astronomy domain's subsolar point |

Drag to pan, right-drag to turn and tilt, and wheel to zoom. Picking a building, a way or a feature of a layer shows it, with actions that open the analyses about it.

## 3. Plug-ins

| Plug-in | Content | Status |
|---|---|---|
| [Street Oceanography](plugins/street-oceanography/README.md) | Six layers and eight analyses (§4) | draft |
| [Landscape Dynamics](plugins/landscape-dynamics/README.md) | Four layers and four analyses: the district as a complex (its buildings among its cells), the landscape running, its spacetime, a building's rooms; the district's complex, freedom and dynamics; its events, exactly; scales and memory | draft |
| [Event Spacetime](plugins/event-spacetime/README.md) | One layer and one analysis: the district's real events as a 2-complex in spacetime (volumes as worldlines, interfaces as sheets, events as vertices, the network of the moment at the clock's date); the record and its complex. The amplitudes and their modes are being written (the design, §8) | draft |
| [Human flow](plugins/human-flow/README.md) | Two layers and three analyses: the district's walkers as agents on the map, moving with the clock (by pattern or by what they do, over the streets' density, with queues at full streets); the counters, counted beside the walkers' reading; the walkers of the district; a street closed or one-way; an event's crowd, its simulated record and its street's regime decided online by POMCP. Melbourne first | draft |
| [The territory](plugins/territory/README.md) | Seventeen layers and three analyses: the world's people, poverty, lights, settlements and hazards, on one grid of hexagons and on the ground, anywhere; any indicator; one's own composite; what stands behind them | draft |
| [Flood risk](plugins/flood-risk/README.md) | Four layers and two analyses: Freetown's flood maps, buildings, wards' losses and the Council's records, its loss tables recomputed; hazards by name anywhere | draft |
| [Reach and use](plugins/access/README.md) | Two layers: a wave of minutes down the streets from any place, and the space usage index on Melbourne's counters | draft |
| Architecture | Form, space and order; the frontage as an interface; plans and districts as cell complexes | planned |
| Cultural analysis | Display and description; the record of manners (a survey at fixed positions, after Kon Wajirō) | planned |

A plug-in is a directory `plugins/<name>/`, and nothing in the shell names one. A layer or a tool may name the areas it serves (`areas`), and is offered there alone. In New York, Street Oceanography offers the tide, the currents and the sensors' weeks, its rhythm and its weather and holidays, all on the Brooklyn Bridge's counter; Landscape Dynamics is offered whole.
- **The manifest.** `plugin.json` names the plug-in's layers and tools, or, if it is planned, what it will hold.
- **A layer** (`layers/*.js`) is given the map, the frame, the clock and the engine. It draws on the map, follows the clock, and answers for what is picked.
- **A tool** (`tools/*.js`) opens in the drawer beside the map. It can select, fly the map, and take what is picked.
- **An engine.** A plug-in that computes has one in C++, `engine/main.cpp`.

## 4. Street Oceanography

**Layers.** The street as a sea:
- **The tide on the streets.** Each sensor's block of street is a channel of water standing to the level of its walkers, a metre for every hundred an hour. The level is the tide the calendars predict, fitted to the counts since October 2024; it rises and falls as the clock runs. Three more readings: the count observed on a date in the record, coloured by its surge against the tide, with a tide mark at the prediction; the hour of high water; and the kind of tide, two crests a day or one.
- **Currents: walkers each way.** At the 92 sensors that count each way, two lanes of walkers, as many as walk there on average at the clock's hour (N = λL/v), moving in real time, keeping left. Coral walk toward the centre of the counts, sea blue away; an arrow is the net current. On a date in the record, that day's counts.
- **Waves on Swanston Street (a model).** The crowd walking north from Flinders Street as a kinematic wave, second by second. Its blocks and crossings come from the map, with a signal at each street. Walkers enter at the counted rate and queue at red. The carriageway can be opened, and the crowd doubled or quadrupled.

The records:
- **Pedestrians at the sensors.** Melbourne's sensors as columns, the mean count at the clock's hour of the week (2024–2026).
- **Shops and services.** Holdings along the named streets, coloured by use. Those open at the clock's hour by their declared hours stand tall.
- **Businesses by use.** Melbourne's census of land use about the sensors.

**Analyses:**
0. **Tides and currents of the city.** When the city's high water comes; which streets have two crests a day; the commuters' flood toward the centre and its ebb, with what the counting lines take in; the arrow of time on the streets, the walkers' entropy production bounded below by the measured edges; which way walkers keep (29 of 41 footpaths keep left); and the surge on the clock's date.
1. **Rhythm of a street.**
   - The record and its average week.
   - The spectrum against red noise, its lines read on the calendars (the week with the day as its harmonics, the lunar month, the year).
   - The harmonic constants of a standard list of constituents, and how selectively the street answers the calendars.
   - Whose clock it keeps when daylight saving begins and ends, measured at every change since 2009.
   - The coupling left in what the calendars do not explain.
2. **Weather and holidays.** How much rain, warmth and cold, and holidays change the counts, with the hour of the week and the season taken out first; every street side by side.
3. **Mix of uses.** How far a street departs from an even mix, and at what scale. Also: the distance over which its uses change, the census years, and why no single score ranks streets.
4. **Opening hours.** Where and when a street is open; the rhythm of its hours (a shop closing at five and a bar at two are different frequencies); who opens and closes first; and which kinds of people its open shops draw at each hour (a model).
5. **Stores that work together.** Pick two stores on the map, one that raises a need (a bridal store) and one that meets it (a jeweller). The pathway between them is solved exactly at their real separation and declared hours: how the chance falls with distance and with legibility, and what one creates for the other. On the map, the first store's visitors spread along the street minute by minute, rose where the need has been raised; an arc joins the stores, as thick as their coupling. Also: which uses keep company beyond the blocks they share.
6. **Crowds and circulation.** Swanston Street on an ordinary day and with its carriageway opened (kinematic waves, the crowd conserved to rounding), and circulation in a district of blocks.
7. **Kinds of streets.** Melbourne's streets placed by how finely their uses change and how selectively they answer the calendars, and whether they fall into groups. A point picked flies the map there.

## 5. Running It

```
python applications/urban-science-workbench/serve.py --open
```

Open http://127.0.0.1:8970/. The data lamp runs the recipe (`recipe.json`):
- it compiles the engine with g++ (natively, or in WSL);
- it fetches the records: Melbourne's pedestrian counts (2009–2022 and the live two years), sensors and census of land use, Open-Meteo's hourly weather, Victoria's public holidays, and OpenStreetMap's holdings and city model; for New York, Valdrada's tiles of Lower Manhattan (about 26 MB, at a pinned revision), the Brooklyn Bridge's counts, New York's weather and holidays;
- it builds the products.

The records take several minutes; they go to `data/cache/`, and none is committed. The server binds 127.0.0.1 unless `--host` says otherwise. When it starts, it asks the engine for every street's tide and current meter (the requests a manifest names under `warm`), so that the map finds them ready.

The engine can be asked directly:

```
python applications/urban-science-workbench/engine.py ask street-oceanography '{"op": "profiles"}'
```

An operation that lists sensors takes an area, Melbourne's when none is named: `'{"op": "tides", "area": "new-york"}'`.

## 6. Design

- **Colour.** A gentle sea blue (#5B8FB9), chosen by the user, with sand (#C98B6B), the shore beside it, as its minor colour; misty paper by day and a deep-sea night.
- **The city.** Pale stone and mist by day, deep sea at night, lit by the Sun with filmic tone mapping.
- **Data.** In the visualization layer's validated maps: `sea` for the water of the tide and the crowd, `sea-brick` for a surge, `warm` (cyclic) for the hour of high water, `gold` for counts, `categorical-apart` for uses. Walkers toward the centre are coral and away sea blue, apart in every vision.
- **Type and lines.** Source Serif 4, and hairlines.

## 7. Checks

`node applications/urban-science-workbench/tests/run-tests.mjs` checks:
- the manifests, the areas, the recipe and every page module's syntax;
- when the engines are built, their answers (Landscape Dynamics on a grid district with a hole: cells, rewrites, a day's run, the skeleton, the scales, a building; Human flow on Melbourne's district: the calibration each way, the walkers' counters reading what the calibrated flow expects, a quarter's frames, a closure, an event's crowd and its simulated record, the permits, one decision of POMCP);
- Human flow's colours, apart in every vision;
- when the records are built, a sensor's week, rhythm and clock, a street's hours and mixing, every street's tide and a day of the sea (the counts within 15% of the tide on an ordinary Wednesday, each way summing to the total), the current meters and their handedness, a street's crowd second by second, and a visitor's plume.

The GIS's arithmetic (the frame against the WGS 84 meridian arc, the round trip, rings, ribbons) is checked in the visualization layer's suite. The mathematics is checked in the foundations and in the urban-science domain.

## 8. Dependencies

- the urban-science and architecture domains and the foundations calx, morphe, physis, tyche and nomos (Human flow's walkers stand on nomos's agents);
- the astronomy domain (the Sun);
- the visualization layer (`gis`, the styles);
- three.js 0.170.0 from jsDelivr;
- the API's C++ JSON;
- the data library's manifests and products.
