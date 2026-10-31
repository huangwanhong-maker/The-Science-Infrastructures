# Urban Building Workbench

*A three-dimensional GIS of cities, with the analysis of streets and buildings on top, built of plug-ins.* Status: **draft** (2026-10-31). Activities: research, learning. [ADR 0060](../../development/decisions/0060-streets-as-open-systems.md).

## 1. Purpose

The workbench is a GIS first. A city stands in three dimensions — its buildings at their mapped heights, its streets at their widths, its parks, water and squares — in an exact local frame on WGS 84. A clock runs in the city's own time zone and lights the city by the real Sun of the instant. On that base, plug-ins put layers on the map and analyses beside it. The city is Melbourne's centre, where the records are.

The first plug-in, **Street Oceanography**, reads streets as open systems, and draws them on the map as a sea:
- the tide of a street's walkers, rising and falling with the clock, and its surge on a date;
- its currents, walkers each way at their true numbers, and the commuters' flood and ebb;
- the waves of a crowd at the crossings of a street (a model);
- how a street answers the weather and the calendar;
- how its uses mix and work together, a store's visitors spreading toward another on the map.

It reads real records of Melbourne. Every result says whether it was measured from records, proved for any street, or comes from a model with declared parameters, and gives how it was obtained and its limits.

## 2. The Map

| Part | What it draws | Source |
|---|---|---|
| Buildings | Extruded footprints with their holes, at their mapped height (from `height`, else levels × 3.2 m); paler where no height is mapped and a default stands. They can stand as built, lowered to a fifth of their height, or of glass, so that the streets show between towers | OpenStreetMap (data library: `city-model.json`) |
| Streets and ways | Ribbons at their mapped or estimated widths; pedestrian ways tinted | OpenStreetMap |
| Parks, water, squares | Flat areas | OpenStreetMap |
| The clock | A date and an hour in the city's time zone (daylight saving its own), playing through the day; the Sun's altitude and azimuth there, lighting the city and casting shadows | the browser's time-zone database; the astronomy domain's subsolar point |

Drag to pan, right-drag to turn and tilt, and wheel to zoom. Picking a building, a way or a feature of a layer shows it, with actions that open the analyses about it.

## 3. Plug-ins

| Plug-in | Content | Status |
|---|---|---|
| [Street Oceanography](plugins/street-oceanography/README.md) | Six layers and eight analyses (§4) | draft |
| Human flow | Trajectories and origin–destination; what people carry at entry and departure; crowds beyond the kinematic wave | planned |
| Architecture | Form, space and order; the frontage as an interface; plans and districts as cell complexes | planned |
| Cultural analysis | Display and description; the record of manners (a survey at fixed positions, after Kon Wajirō) | planned |

A plug-in is a directory `plugins/<name>/`, and nothing in the shell names one.
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
python applications/urban-building-workbench/serve.py --open
```

Open http://127.0.0.1:8970/. The data lamp runs the recipe (`recipe.json`):
- it compiles the engine with g++ (natively, or in WSL);
- it fetches the records: Melbourne's pedestrian counts (2009–2022 and the live two years), sensors and census of land use, Open-Meteo's hourly weather, Victoria's public holidays, and OpenStreetMap's holdings and city model;
- it builds the products.

The records take several minutes; they go to `data/cache/`, and none is committed. The server binds 127.0.0.1 unless `--host` says otherwise. When it starts, it asks the engine for every street's tide and current meter (the requests a manifest names under `warm`), so that the map finds them ready.

The engine can be asked directly:

```
python applications/urban-building-workbench/engine.py ask street-oceanography '{"op": "profiles"}'
```

## 6. Design

- **Colour.** A gentle sea blue (#5B8FB9), chosen by the user, with sand (#C98B6B), the shore beside it, as its minor colour; misty paper by day and a deep-sea night.
- **The city.** Pale stone and mist by day, deep sea at night, lit by the Sun with filmic tone mapping.
- **Data.** In the visualization layer's validated maps: `sea` for the water of the tide and the crowd, `sea-brick` for a surge, `warm` (cyclic) for the hour of high water, `gold` for counts, `categorical-apart` for uses. Walkers toward the centre are coral and away sea blue, apart in every vision.
- **Type and lines.** Source Serif 4, and hairlines.

## 7. Checks

`node applications/urban-building-workbench/tests/run-tests.mjs` checks:
- the manifests, the areas, the recipe and every page module's syntax;
- when the engine is built, its answers;
- when the records are built, a sensor's week, rhythm and clock, a street's hours and mixing, every street's tide and a day of the sea (the counts within 15% of the tide on an ordinary Wednesday, each way summing to the total), the current meters and their handedness, a street's crowd second by second, and a visitor's plume.

The GIS's arithmetic (the frame against the WGS 84 meridian arc, the round trip, rings, ribbons) is checked in the visualization layer's suite. The mathematics is checked in the foundations and in the urban-science domain.

## 8. Dependencies

- the urban-science domain and the foundations calx, morphe, physis and tyche;
- the astronomy domain (the Sun);
- the visualization layer (`gis`, the styles);
- three.js 0.170.0 from jsDelivr;
- the API's C++ JSON;
- the data library's manifests and products.
