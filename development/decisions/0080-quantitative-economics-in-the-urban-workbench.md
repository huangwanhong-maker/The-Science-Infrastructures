# 0080. Quantitative Economics in the Urban Science Workbench: One Toolbox, Its Answers on the Map

- **Status:** Accepted, provisionally (2026-11-07). Amends ADR 0056 (the toolbox's home) and ADR 0079.
- **Date:** 2026-11-07
- **Decided by:**
  - **The request.** The user: "Next, please integrate the (quantitative) economics toolbox into it", the Urban Science Workbench.
  - **The design.** Claude, rooted in what the infrastructure already held: the Economy Atlas's 27 tools (ADR 0056), the workbench's line of time (ADR 0079), and the rule that nothing depends on an application.

## 1. Context

The toolbox stood inside the Economy Atlas: its engine headers in `engine/tools/`, its catalogue, charts and tables in the atlas's `src/` and `tables.py`. The workbench could not use them without depending on another application (CLAUDE.md §3.4), and copying them would have made two toolboxes to keep.

A toolbox in a GIS also owes the map more than a drawer of numbers. Its rows are often places, so its answers should stand where those places are. A regression on places reports standard errors that assume the places are independent, and neighbours seldom are.

## 2. Decision

### 2.1 The Toolbox Moves Down

| Part | From | To |
|---|---|---|
| The engine's tools | `applications/economy-atlas/engine/tools*.hpp` | `api/cpp/science/toolbox.hpp`, `toolbox/` (namespace `science::toolbox`) |
| The catalogue, form, answer and charts | the atlas's `src/tools.js`, `src/toolbox.js`, `src/charts.js` | `visualization/js/toolbox/` |
| The library's tables | the atlas's `tables.py` | `data/products/economic_tables.py`, kept in `data/cache/products/tables/` |

They were moved by `git mv`, so that their history follows. The Economy Atlas now uses them, and keeps only its rail, its own tables and its map. Its suite is unchanged and passes. The charts take their colours from the page's stylesheet (`--tb-series-1` and so on), so that each application keeps its palette.

### 2.2 What a Map Asks of the Toolbox

- **Values by row.** A regression asked for them (`by_row`) gives its residuals and fitted values, and a panel its units' effects, for the rows it kept.
- **For each place.** Any tool on a table can be run for each part of it (`"each": [columns]`): its rows are split, the tool is run on each part, and the numbers of its answers come back as measures by key. A part it cannot answer is left out and counted by reason, never given a zero. Every tool gains this, and none need know of it.
- **Neighbours in space.** `moran`, on tyche's new spatial autocorrelation (tyche's decision 0013, held to PySAL's esda), gives Moran's I with hot and cold spots, neighbours being found on the sphere. A regression whose rows are distinct places also tests its residuals, and says when its standard errors are too small.
- **Location quotients.** `specialisation`, on nomos's measurement: where each kind of activity gathers, and how far each kind lies apart from the rest.

### 2.3 The Plug-in *Quantitative Economics*

It stands under a seventh lens, *Economy*, as nine ways in: questions of the place, and the toolbox's groups.

- **Questions in plain words**, each a tool's form filled in, for Melbourne, Freetown, anywhere, and the world. Each says in a sentence what it shows and what it does not.
- **Any tool, on any table here:**
  - Melbourne's blocks by census year (floor space by use, businesses by kind, dwellings);
  - its counters by day with the weather, and by month since 2009;
  - Freetown's wards;
  - the Space2Stats cells about the point looked at, and those cells year by year;
  - the world's economies and currencies;
  - what is shown in the studio;
  - a table pasted.
- **Answers on the map.** What a tool finds of places becomes a record of the line of time (ADR 0079): cells for polygons, small hexagons for points, the world's mosaic for economies (in which an area in dispute has no economy). Its time steps run on the line of time, and it carries a card, a hover and a list of the highest and lowest.
- **Drawn with care.** A city's blocks lie as a carpet beneath the towers, and counters stand as pillars above them. By night, a value near zero sinks into the dark instead of glaring in a diverging map's pale middle. A record of one moment opens no line of time.

### 2.4 The City's Tables

`data/products/city_tables.py` makes them from what the data library already holds:

- **What the census withheld.** A cell withheld for confidentiality stays null, never zero.
- **A civil day** is one in Melbourne, counted only if it has 23 hours or more. A day on which a counter counted nobody is a counter out of order, and is left out.
- **The gap.** The City published nothing from November 2022 to October 2024, and nothing stands in those months.

## 3. Results

Every question answers on the real tables. What they find:

- **The return after the lockdowns.** Melbourne's counters came back unevenly from 2019 to 2025: of 45 counters, the log ratio runs from −1.28 (Flinders Underpass) to +1.03.
- **Rain.** A millimetre of rain takes 1.4% of a day's walkers (counters' fixed effects, 54,869 counter-days, s.e. 0.07%).
- **Businesses and floor space.** A block's businesses rise with its floor space at an elasticity of 0.78 (R² 0.43, 519 blocks in 2024). Its residuals are spatially autocorrelated (Moran's I 0.44, p 0.002), so that the reported standard errors are too small, and the page says so.
- **Food and drink** cluster (I 0.54): 67 blocks are hot spots, against 118 cold.
- **The spread of businesses over blocks** grew more unequal, its Gini coefficient rising from 0.70 in 2002 to 0.74 in 2024.
- **Freetown's loss per km²** does not cluster (I −0.02, p 0.48), and the page says that too.
- **The world:**
  - poorer economies grew faster from 2000 to 2019 (β −0.0079 a year, p 8 × 10⁻¹⁴);
  - gravity's distance elasticity is −0.94.

- **A found fault: questions waited behind tiles.**
  - **The cause.** A browser gives a host six connections. The globe asks for up to twelve terrain tiles at once, and the city's tiles come from the same server, so a question asked while the camera moved could wait minutes.
  - **The fix.** Tiles are now asked at low priority and the engines at high, so that a question takes the next connection to free.
  - **Measured** in the same headless page: a question after a camera move answered in 11 s, against more than 150 s before. The first question after the page opens still waits about a minute while the server builds the terrain it has not built before.

Pictures are in the [gallery's §18](../../docs/gallery.md).

## 4. Consequences

- `science::toolbox` is the API's, and any application's engine may answer it: the Economy Atlas's and the workbench's do.
- The Economy Atlas gains `moran` and `specialisation` in its catalogue, and the economies table its places' longitude and latitude.
- **Checks:**
  - tyche's 44;
  - `check_economics.py` in the workbench's suite (every question, each held to its parts, residuals by row, location quotients and Moran's I by hand, the tables against their sources);
  - the visualization layer's catalogue.

## 5. Open

- A tool for each place runs its parts one after another. A table of thousands of places with a slow tool (ARIMA) takes a while.
- The answers are not kept between sessions.
- Local tests are not corrected for their number, and spatial regression models (lag, error) are not yet estimated as such.
- New York has no table of its own yet: the cells about the point looked at serve it.
- The page's tiles and its engines share one server's six connections; a second origin for tiles would part them entirely.
