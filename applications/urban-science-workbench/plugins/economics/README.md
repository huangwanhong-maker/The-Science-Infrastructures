# Quantitative Economics

*The toolbox of quantitative economics on the place open, its answers on the map.* Status: **draft**. Activity: use. [ADR 0080](../../../../development/decisions/0080-quantitative-economics-in-the-urban-workbench.md).

## 1. Purpose

A question put to a place's records in plain words, or any tool of the toolbox on any table here, answered by the foundations and drawn where the places are.

## 2. Ways In

Under the lens *Economy*:

- **Questions of this place.** Each is a tool's form filled in ([`questions.js`](questions.js)), saying in a sentence what it shows and what it does not.
- **The toolbox's groups:**
  - regression and panels;
  - inequality and concentration;
  - neighbours in space;
  - time series;
  - summary and correlation;
  - trade;
  - growth and the short run;
  - games and choice.

## 3. Tables

[`tables.js`](tables.js) reads the tables from where one stands:

| Table | Where | Its rows |
|---|---|---|
| `melbourne-blocks` | Melbourne | each CLUE block in each census year, 2002–2024 |
| `melbourne-counters-day` | Melbourne | each counter on each day since October 2024, with the weather |
| `melbourne-counters-month` | Melbourne | each counter in each month since May 2009 |
| `freetown-wards` | Freetown | each ward, with the assessment's annual loss |
| `cells-here`, `cells-years` | anywhere | Space2Stats's cells of the 4° square about the point looked at |
| `economies`, `currencies`, `currency-rates` | the world | the data library's tables |
| shown | anywhere | each record shown in the studio |
| pasted | anywhere | a table pasted, its rows places if it has columns that say where |

The library's tables are made by [`data/products/city_tables.py`](../../../../data/products/city_tables.py) and [`economic_tables.py`](../../../../data/products/economic_tables.py).

## 4. Answers on the Map

[`answers.js`](answers.js) makes what a tool finds of places into a record of the line of time (ADR 0079):

- residuals and fitted values by row;
- a panel's effects;
- a tool's numbers for each place;
- an economy's value on its cells of the world's mosaic.

Any tool on a table with places can be run **for each place**, each moment, or each place at each moment.

## 5. Engine

[`engine/main.cpp`](engine/main.cpp) answers the API's toolbox (`api/cpp/science/toolbox.hpp`). It reads no file but the world's trade, for the tools that need it.

## 6. Limits

- A tool for each place runs its parts one after another.
- Local tests of Moran's I are not corrected for their number.
- The answers are not kept between sessions.
