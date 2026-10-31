# Street Oceanography

*Streets as open systems, on the map.* A plug-in of the [Urban Building Workbench](../../README.md). Status: **draft** (2026-10-31).

The analyses come from the user's working papers on street oceanography. The urban-science domain states them exactly, with their proofs ([theory/street-oceanography.md](../../../../domains/urban-science/theory/street-oceanography.md)). This plug-in draws them on the workbench's 3D map and beside it, in plain terms.

## 1. Layers and Tools

| Layer | Records |
|---|---|
| `layers/tide.js`: the tide on the streets, its surge, high water and kind | City of Melbourne, hourly counts (2024–2026); Open-Meteo; Victoria's holidays |
| `layers/currents.js`: walkers each way, the net current | City of Melbourne, hourly counts each way (2024–2026) |
| `layers/waves.js`: waves on Swanston Street (a model) | the city model's blocks and crossings; the counts at the street's start |
| `layers/sensors.js`: pedestrians at Melbourne's sensors, by the clock's hour | City of Melbourne, hourly counts (2024–2026) |
| `layers/holdings.js`: shops and services, open or closed at the clock's hour | OpenStreetMap, holdings and `opening_hours` |
| `layers/establishments.js`: businesses by use | City of Melbourne, census of land use (CLUE) |

| Tool | Kind of result | Records |
|---|---|---|
| Tides and currents of the city | measured, proved | counts, each way |
| Rhythm of a street | measured | counts |
| Weather and holidays | measured | counts; Open-Meteo; Victoria's holidays |
| Mix of uses | measured, proved | OSM; CLUE |
| Opening hours | measured; a model of visits | OSM `opening_hours` |
| Stores that work together | a model (the pathway); measured (co-location) | OSM; CLUE |
| Crowds and circulation | a model; proved | the city model's blocks; the counts at Swanston Street's start |
| Kinds of streets | measured | counts and CLUE |

## 2. The Engine

`engine/main.cpp` is C++20 on the urban-science domain, the foundations and the API's JSON. It reads only the data library's products, under the directory the server passes (`--data`). Its operations:
- the catalogue: `status`, `sensors`, `streets`;
- the map: `profiles`, `establishments`;
- terrain and hours: `street`, `visits`, `terrain`, `dominance`;
- affinity: `affinity`, `spacing`;
- records: `rhythm`, `clock`, `bicoherence`, `constructed`, `response`, `summary`, `plane`;
- the sea: `tides` (every street's tide and its readings), `sea` (a date: each street's predicted level, observed counts and counts each way, the day's holiday and weather), `currents` (the current meters and the handedness test);
- models: `pathway`, `plume`, `streetday`, `streetwave`, `district`, `persontime`.

`main.cpp` documents each operation's request.

## 3. Limits

- **What a count holds.** Melbourne's sensors count crossings. What the people carry (purpose, language, acquaintance with the place) is not measured anywhere here, and appears only in models with declared parameters.
- **Sampling of streets.** The sensors were placed for counting pedestrians, not chosen to span kinds of streets.
- **Labels of use.** The bridal use is read from CLUE's trading names; OpenStreetMap's tags are as their mappers left them.
- **Weather.** It is a reanalysis at the city centre, not a gauge on each street.
- **The tide's era.** The tides are fitted to the live record (October 2024 to October 2026), and a date outside it shows the tide alone.
- **Cross-sections, not a cordon.** The counting lines are a sample: people also arrive by train beneath the streets, so what the lines take in is not the city's whole prism.
- **The crowd model.** Nobody turns into a cross street, so its flow further up Swanston Street overstates the counts there.
