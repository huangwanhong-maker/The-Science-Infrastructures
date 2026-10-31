# 0065. New York in the Urban Building Workbench: Lower Manhattan from Valdrada's Compiled City, and the Brooklyn Bridge's Counts

- **Status:** Accepted, provisionally. Amends [ADR 0060](0060-streets-as-open-systems.md) (the workbench's one area) and [ADR 0063](0063-the-district-as-a-complex-to-work-on.md) (its district).
- **Date:** 2026-11-01
- **Decided by:**
  - **The request.** On 2026-10-31 the user asked whether using Valdrada's digital twin of New York ("Boundless NYC") would be sensitive. Claude answered that it is not, in the sense of privacy, and that its licences are open on conditions (below), and proposed New York as a second area of the workbench. On 2026-11-01 the user wrote: "do the New York one".
  - **What Claude did.** It surveyed Valdrada's tiles and New York's open records of pedestrians, built the area and wrote this record.

## 1. Context

The workbench opened one city, Melbourne; Ginza had been removed at the user's word (ADR 0060). Its city model was one file of OpenStreetMap's buildings, ways and spaces. Street Oceanography read Melbourne's counting system, its weather and its holidays, and assumed them throughout. Landscape Dynamics took one district, the Hoddle Grid's middle.

Valdrada (M. K. Turkcan, https://github.com/mkturkcan/valdrada) compiles New York City from its public records into binary tiles, about 13.6 GB for the five boroughs. Its code is MIT-licensed. Its compiled city is a Derivative Database under ODbL 1.0, because it takes building colours, missing footprints and other parts from OpenStreetMap. Anything drawn from it must carry "© OpenStreetMap contributors", and a database derived from it stays ODbL. Its pedestrian models contain MetaHuman components, whose licence forbids building databases or training and testing models with them. Its renderer is under Epic's licence.

New York has nothing like Melbourne's counting system. Business districts publish their counts as PDF reports only, and the city's twice-yearly counts cover three windows a day. The one open record of pedestrians counted each way, every hour, in Manhattan is DOT's counter on the Brooklyn Bridge's promenade at its Manhattan approach, October 2017 to December 2019, now historic.

## 2. Decision

### 2.1 The Area

**New York: Lower Manhattan,** from the Battery to Canal Street and from the Hudson to the Brooklyn Bridge's Manhattan tower: the box 40.6995–40.7190 N, 74.0195–73.9930 W. The bridge's counter stands in it. Landscape Dynamics takes the Financial District's old streets about Wall Street, Stone Street and Hanover Square as its district: irregular colonial blocks, pedestrian lanes and towers rather than a grid.

### 2.2 Each Source What It Holds Best

- **From Valdrada's tiles, fetched at a pinned revision with their checksums (`data/registry/valdrada-nyc.json`):**
  - **buildings**, at the roof heights NYC's Building Footprints record (from aerial survey), with their light courts and the compiler's style;
  - **trees**, from the street-tree census and NYC Parks' tree points;
  - **water**, traced by marching squares on the tiles' ground where it is connected to the harbour, with its piers as holes.
- **From OpenStreetMap:** the streets and the open spaces, by Melbourne's rules, without tunnels.

  *Rationale.* The tiles' street centrelines were made for the renderer. Their header lists only some of the junctions, and most centrelines stop short of the junction they meet, so they do not form a network. OpenStreetMap's ways share their junctions' nodes. OpenStreetMap's Manhattan is mostly without heights, which is where Valdrada's are needed.

- **The reader:** `data/products/valdrada_city.py`, which documents the tiles' layout (the compiler's, v0.6.2) and inverts Valdrada's equirectangular frame exactly to WGS 84 before projecting as Melbourne is projected.
- **Not used:** Valdrada's models, pedestrians, renderers and textures; PLUTO itself (whose owner names would be personal data).

### 2.3 The City Model, an Area to a File

The city model becomes an index (`city-model.json`) and a file for each area (`city-model-<id>.json`). Each area's file carries:
- its sources and its credit, which the map's legend shows;
- the conventions by which it was read;
- its trees, where it has them.

An area whose source is not fetched is left out, with the command that fetches it. The page reads an area's file when the area is opened.

### 2.4 The Records of New York's Streets

- **The products.** The bridge's counts become a sensor of its own area (`street-sensors-new-york.json`). The wall-clock record is placed in UTC by the United States' rule of daylight saving since 2007:
  - the spring change's missing hour is dropped;
  - the autumn change's repeated hour is read twice where it appears twice;
  - days the counter read zero for twelve hours or more are null.

  Beside it stand New York's hourly weather (Open-Meteo) and its public holidays (Nager.Date). An index, `street-areas.json`, names each area's sensors, weather, holidays and time zone.
- **The engine.**
  - A sensor belongs to an area, and its weather and holidays are its area's.
  - The operations that list sensors take an area, Melbourne's when none is named, so every earlier caller is unchanged.
  - A record that ended before Melbourne's live record began is read whole, for the tide and for the recent era.
- **The bridge as a tide.** One crest a day: a form number of 3.1, with weekday high water at 15:30. The calendars explain 79% of its variance. It is a promenade's tide, not a commuter street's two.

### 2.5 The Workbench by Area

- **Areas declared.** A tool, like a layer, may name the areas it serves, and the rail offers it there alone. In New York:
  - Street Oceanography offers its tide, its currents and the sensors' weeks, and its rhythm and weather-and-holidays analyses;
  - its parts that rest on Melbourne's land-use census, its holdings along named streets or its hundred sensors remain Melbourne's;
  - Landscape Dynamics is offered whole.
- **The district.** Each area names its district in `areas.json`. The plug-in's state belongs to the area: opening another forgets the district, its edits, its runs and its buildings.
- **The view.** An area's first view may name its camera's range, as Lower Manhattan's towers ask. Among towers taller than a fifth of the district's width, the workspace looks down more steeply and from higher.

## 3. Consequences

- **What one can do.** Open New York beside Melbourne, both drawn from the same kind of product:
  - see Lower Manhattan's buildings at their surveyed heights, its trees and its shore;
  - work on the Financial District as a complex, its 168 buildings 3-cells sharing 201 walls in 62 terraces: set what its cells carry, rewrite them and run the landscape on them;
  - read the Brooklyn Bridge's promenade as a tide and its two ways as currents, with New York's weather and holidays.
- **Limits.**
  - One sensor, historic, with a five-month gap in 2018. Its neighbourhood has no land-use census or declared hours in the workbench, so the district's stores are declared and the walkers entering it follow a declared day.
  - The weather is the reanalysis grid point 4.4 km away, across the Hudson.
  - A Valdrada height of 4–8 m may be the compiler's default rather than a record; the tiles do not say which.
- **Licence.**
  - The map credits "Valdrada (M. K. Turkcan) from NYC Open Data · © OpenStreetMap contributors".
  - `city-model-new-york.json` is a database derived from ODbL data, and stays ODbL if it is ever published. Like every dataset, it lives in the cache and is not committed.
- **Revisability.**
  - Valdrada's format is at version 0.6.2 and may change. The revision is pinned and its tiles are checked by their checksums.
  - Should the format change, the buildings can be read from NYC's Building Footprints directly.
  - Further areas join by an entry in `areas.json` and a reader of their source.
