# The Economy Atlas

*Status: draft.* An application of the economics domain ([ADR 0044](../../development/decisions/0044-value-in-space-and-time.md), [ADR 0046](../../development/decisions/0046-the-economy-atlas-data-and-plugins.md), [ADR 0056](../../development/decisions/0056-nomos-and-a-toolbox-of-quantitative-economics.md), [ADR 0068](../../development/decisions/0068-the-economy-atlas-as-a-scholarly-instrument.md)). It draws value carried between places and through time:

- the world's trade on a map, or flows of one's own between any places, split into exchange, drive and circulation;
- currencies held and exchanged as a spacetime complex, whose curvature is the departure from a price system;
- production as a circular flow;
- the loop deficits of transfer networks;
- a toolbox of quantitative economics, on tables of the economies, the currencies, or one's own.

Everything shown is computed by the domain's C++ on the foundations. Plug-ins in Python and C++ analyse the flows the map shows and draw their results on it. Serves research, learning and use.

## 1. What It Shows

- **The map.** Flows between places in a year on an Equal Earth map, the geography domain's projection: by default the world's trade in goods (the OECD's balanced bilateral statistics, 2010–2024), or data of one's own (§2.2).
  - What two places send each other both ways is reciprocal exchange. What is left, one net flow for each pair, is decomposed (the domain's theory note, §3) on a complex whose faces are the spherical Delaunay triangles of the places.
  - Each link conducts in proportion to its gross flow. The *driven* part between two places is their gross flow times their difference in attraction α: the potential runs from the places that send value to those that draw it, a share of gross flow without unit.
  - Places are coloured by α; by the flows' other measures (net inflow, exports, imports, the share circulating); by any of the World Bank's indicators in the toolbox's table of economies (GDP per capita, population, life expectancy, CO₂ per person, the Gini index and others), in the map's year or the latest of the five before it, on a log scale where they span decades; or by the values a plug-in or a tool returned.
  - *Local circulation* turns round the triangles of neighbours, drawn as rose eddies. *Global circulation* runs round the loops of distant partners no triangle fills.
  - Arcs are the largest net flows, their driven or circulating parts, or a plug-in's arcs. Each wears the colour of its largest part: indigo driven, rose round the triangles, gold round the holes.
  - The map zooms (the wheel, the buttons) and pans (dragging), and fits a regional dataset. Places without a country of their own are circles. The pointer brings a card for a place or an arc. A place is found by its name or code, or clicked: its card gives its flows, attraction and largest net flows, each with the share of it driven and circulating, and its trade over the years beside the indicator the map is coloured by.
  - A strip of years beneath the map gives every year's shares of the net flows' energy; a click, or ← and →, chooses the year, and **play** goes through them, each year's decomposition fetched while the one before is shown.
  - The dock and the inspector fold away (‹ controls, inspector ›), and the map takes the room. Its instruments name the projection, the longitude and latitude under the pointer and the zoom; each legend ends with its source, year, unit and projection.
  - The year's places and net flows are saved as CSV, and the map as a figure composed for publication: a title and a line of what is drawn, the map as shown with its arcs solid and headed, the legend, and a line of source, projection, date and link; as SVG with Source Serif 4 written into it, or as PNG at 3600 pixels wide.
- **The toolbox.** Quantitative economics on tables (§2.4): the library's (economies by year with their trade and the World Bank's indicators; currencies by month), one's own, or one pasted. Each tool is a form built from its declaration, run by the engine on the foundations' methods, and answered in figures, tables, charts and, for values by economy, a layer of the map.
- **Spacetime.** Currencies held month by month at their policy rates and exchanged with a base at the month's rate, as the spacetime complex of the theory note, §5. The library's record holds twenty currencies against the euro: the ECB's rates, the BIS's policy rates, and for the euro the ECB's deposit facility rate. A record of one's own uses its own base.
  - A sheet shows every lens's curvature, the base's excess return over a currency in a month, the departure from uncovered interest parity.
  - Lines show its sum in time.
  - The panel measures the curvature against that of the exchange rates alone: the share of it the rates of holding remove. This is a gauge invariant; the share of the transports a price system would carry is not, and is not shown.
- **Production.** An input–output table, edited in place (Miller and Blair's example, and a three-sector circular economy), solved by the Leontief model:
  - technical coefficients, the inverse, multipliers, prices and linkages;
  - cost-push and demand shocks;
  - the deliveries split into exchange, drive and circulation.
- **Circuits.** A transfer network, edited in place, with the examples of the user's paper on loop deficits (its triangle, its five positions, a fee round a circuit, a market with prices, a round trip):
  - every loop's holonomy and deficit in the tree gauge;
  - the price potential where the deficits vanish;
  - retention round directed cycles;
  - the circulation of flows and their pairing with the transports.

The page is a scholarly instrument in porcelain and gold (ADR 0068): porcelain surfaces and graphite ink, a champagne gold only for what is chosen or in focus, rose red as its minor colour; graphite by night, when chosen (the button at the top right). Its data's colours are kept apart from the interface's and validated for colour vision deficiencies on both themes (ADR 0046, §2.1; ADR 0056, §2.5). Source Serif 4 throughout, with tabular figures; 3 px corners and hairlines.

## 2. Running

From the repository root, make the engine and the data once (g++ 11 or later, natively or in WSL; about 60 MB from the OECD, the ECB, the BIS and the World Bank), and serve the page:

```
python data/acquisition/recipe.py applications/economy-atlas/recipe.json
python applications/economy-atlas/serve.py --open
```

The page is at http://127.0.0.1:8910/. It offers to run the recipe itself when something is missing.

### 2.1 Links

The address holds the view and what it shows, and on the map follows every change, so that a view, or a figure, can be linked and drawn again:

- `#map&year=2010&arcs=driven&economy=CHN`:
  - `year`;
  - `arcs`: `all`, `driven`, `circulating` or `none`;
  - `shown`: how many arcs;
  - `faces=0`, `labels=0`, `flowing=0`: hide the triangles, the names or the movement of the arcs;
  - `economy`: a place's id, the ISO alpha-3 code for the library's;
  - `colour`: `attraction`, `balance`, `exports`, `imports`, `circulating`, or `wb:` and a column of the economies table (`wb:gdp_per_capita`);
  - `dataset`: a dataset of one's own, by its id.
- `#toolbox&tool=ols`: a tool, by its id, run on its example.
- `#spacetime&from=2015-01&to=2026-06&currency=MXN`, with `dataset` as well.
- `#production&table=three` (or `miller-blair`).
- `#circuits&network=five` (or `triangle`, `fee`, `prices`, `round-trip`).

### 2.2 One's Own Data

**+ your data**, on the map or in spacetime's list of records, takes CSV or JSON files ([datasets.py](datasets.py)). They are kept as products of the library's own schemas under `data/cache/economy-atlas/datasets/` and drawn, and given to the plug-ins, as the library's are.

- **Flows between places.**
  - A file of flows: `from,to,value[,year]`.
  - For places that are not the library's economies, a file of places: `id,name,latitude,longitude[,group]`, in degrees on WGS 84.
  - Ids that are ISO alpha-3 codes or the names of the library's economies need no place, and every spelling of one economy (CHN, chn, China) is one place. Flows without a year belong to one undated year.
  - Or one JSON: `{"title", "unit", "places": [...], "flows": [...]}`.
- **Currencies held and exchanged.**
  - `month,currency,rate,policy`: the month as YYYY-MM; the rate as units of the currency for one unit of the base; the policy rate in per cent a year.
  - The base's rows give its own policy rate, with a rate of 1 or none, and the base is found from them, or named.

A file that cannot be read is refused with what to change: an unknown id, a negative flow, missing columns, a word for a number, two places at one point. [examples/](examples/) holds invented data to try the formats. From the command line: `python applications/economy-atlas/datasets.py add flows "<title>" places.csv flows.csv`.

### 2.3 Plug-ins

A plug-in is a Python script or a C++ program in [plugins/](plugins/), found by the page (**look again**) and described by itself: its title, what it does and its parameters, from which its form is built. **Run on this map** gives it the network the map shows:

- its places;
- its gross flows in the year;
- the engine's decomposition of them: the attraction of each place, the driven, local and global parts of each net leg, and the circulation round each triangle.

What it returns becomes layers:

- `values`, a number for each place, colour the places (diverging or sequential);
- `arcs` are drawn between places;
- `sketch`es (lines, markers, labels in longitude and latitude) are drawn over the map;
- `table`s and `note`s go to the panel.

The contract is the API's ([api/README.md §2.7](../../api/README.md), `science.networks` and `api/cpp/network_plugin.hpp`), and [plugins/README.md](plugins/README.md) shows a skeleton in each language. A Python plug-in runs as a process of its own, so it cannot take the server down and every run can be stopped. A C++ one is compiled at its first use with the economics domain's and the foundations' include directories, so it may use the domain itself.

Three examples:

- `gravity.py`, Tinbergen's gravity equation;
- `concentration.py`, the concentration of each economy's exports;
- `stationary.cpp`, where value settles if each economy passed its stock on in its export shares, on `economics::stationary`.

Only files in `plugins/` are run. On a public address, plug-ins and datasets need the server's token, as the recipe does.

### 2.4 The Toolbox

Twenty-seven tools of quantitative economics ([ADR 0056](../../development/decisions/0056-nomos-and-a-toolbox-of-quantitative-economics.md)), in seven groups, each wearing the name of the library that computes it:

| Group | Tools | On |
|---|---|---|
| Statistics | summary and correlation (Pearson and Spearman) | tyche |
| Econometrics | linear regression (classical, HC0–HC3 or clustered errors, weights); logit, probit and Poisson (PPML) with average marginal effects; panels (fixed effects one- or two-way, random effects, Hausman's test); instrumental variables (2SLS, the first stage's F, Sargan and Hansen, least squares beside it); differences in differences (2 × 2 and two-way fixed effects); the convergence of incomes (β and σ) | tyche, nomos |
| Time series | autocorrelation (ACF, PACF, Ljung–Box); unit roots (ADF, KPSS); ARIMA and forecasts; trend and cycle (Hodrick–Prescott); vector autoregressions (order, stability, impulse responses); cointegration (Engle and Granger) | tyche |
| Inequality and indices | inequality (Gini, Theil, Atkinson, Palma, the Lorenz curve, Theil's split by group); poverty (FGT, Watts, the incidence curve); concentration (HHI, concentration ratios); index numbers (Laspeyres, Paasche, Fisher, Törnqvist, direct or chained) | nomos |
| Trade | the gravity of trade (PPML with exporter and importer effects, PPML with sizes, least squares on logarithms); the concentration of each economy's partners; trade intensity with every partner | economics, on nomos and tyche |
| Growth and macroeconomics | the Solow model; IS–LM with a policy's shift; growth accounting | economics |
| Games and choice | two-player games; zero-sum games; stable matching; voting rules | nomos |

A tool's form is built from its declaration ([`src/tools.js`](src/tools.js)): a table and its rows (`where year = 2019`), the columns its roles take, and its choices. Every tool opens on an example and runs at once. Where the library's data cannot carry a valid design, the example is simulated and says so: markets whose cost of supply is a valid instrument for their price, regions treated from 2020.

- **The tables.** The library's ([`tables.py`](tables.py)): the economies by year, 2000–2024, with their trade (OECD BIMTS) beside eleven of the World Bank's indicators and its regions and income groups; the currencies by month, long and wide. One's own: **+ your data** takes a CSV or JSON table (a column is numbers when every filled cell is one). Or a table pasted into the form.
- **The answers** are figures, tables (copied, or saved as CSV), charts (lines with their intervals, bars, scatters with their fit, matrices) and notes. Values for the economies (growth of income, a tool's effects, a partner's intensity) go to the map as a layer, by **Show on the map**.

`#toolbox&tool=gravity` opens a tool by its id.

## 3. The Protocol

The engine is `engine/economy.cpp`, header-only C++20 on:

- the economics domain (`domains/economics/`);
- physis: discrete exterior calculus, flows with conductances, the spherical Delaunay triangulation;
- morphe: complexes, the worldline construction;
- calx: linear algebra, integrators and solvers;
- tyche: statistics, regression and time series (the toolbox);
- nomos: the methods of the social sciences (the toolbox);
- the API's JSON.

The server runs it once per request: one JSON request on its standard input, one JSON answer on its standard output, an error as `{"error": "..."}`. Products are passed by their paths; the engine reads nothing else.

### 3.1 Requests

| Request | Answer |
|---|---|
| `{"op": "map", "flows": "<trade-flows.json>", "year": y}` | the year's places, net flows and the map's triangles (§3.2) |
| `{"op": "timeline", "flows": "<trade-flows.json>"}` | `{"years": [{"year", "total", "reciprocal", "net", "driven", "local", "global"}], "unit"}`: every year's totals and shares |
| `{"op": "spacetime", "record": "<currency-spacetime.json>", "first": "YYYY-MM", "last": "YYYY-MM", "currencies"?: [...]}` | the lenses, statistics and price system of the currencies with complete records over the span, the record's base first (§3.3) |
| `{"op": "production", "table": {"sectors", "z", "final_demand", "value_added"}, "value_added_coefficients"?, "final_demand"?}` | the Leontief model and the deliveries' circulation (§3.4) |
| `{"op": "circuits", "network": {"positions", "legs": [{"from", "to", "ratio", "flow", "kind"}], "faces"?}}` | the network's cycles, prices and circulation (§3.5) |
| `{"op": "tool", "tool", "columns"?: {name: [...]}, "roles": {role: column(s)}, "params", "flows"?, "sizes"?}` | a tool's answer: `{"title", "blocks": [{"kind": "figures" \| "table" \| "lines" \| "bars" \| "scatter" \| "matrix" \| "map" \| "note", …}]}`, as `engine/tools.hpp` describes |

### 3.2 The Answer to `map`

```
{"year", "years", "unit",
 "economies": [{"iso3", "name", "continent", "latitude", "longitude", "exports", "imports",
                "balance", "attraction", "circulating"}],                // iso3: a place's id
 "legs": [[from, to, net flow, driven, local, global]],              // net legs, in the product's unit
 "triangles": [{"places": [i, j, k], "curl"}],                         // the eddy round each
 "totals": {"total", "reciprocal", "net", "driven", "local", "global"},  // the last three: shares of the energy Σ F²/c
 "faces": {"triangles", "neutral_legs", "covers_sphere", "neutral_conductance"}, "method", "source"}
// attraction: a share of gross flow, without unit; circulating: of the driven and circulating parts on a place's links, the share that circulates
```

Fewer than four places make no map's faces: their flows are decomposed on their links alone.

### 3.3 The Answer to `spacetime`

```
{"base", "currencies", "names", "months", "times",
 "lenses": [[currency, month from, month to, curvature]],
 "stats": [{"currency", "mean", "volatility", "removed", "cumulative": [...]}],  // a year; the base's excess return over the currency
 "presences": [[currency, month, log price]],                         // the nearest price system's potential
 "policy", "rates",
 "curvature": {"rms", "rms_exchange", "removed", "distance"},         // against the exchange rates' alone; removed = 1 − Σκ²/Σκ₀²
 "cells": {"presences", "edges", "lenses", "triangles"}, "method", "iterations", "source"}
```

### 3.4 The Answer to `production`

```
{"sectors", "A", "L", "output", "multipliers", "value_added_coefficients", "prices", "backward", "forward",
 "productive", "imbalance",
 "deliveries": {"legs": [[from, to, net, driven, local, global]], "attraction", "reciprocal", "net", "driven", "local", "global"},
 "shocked_prices"?, "shocked_output"?}
```

### 3.5 The Answer to `circuits`

```
{"positions", "cycles": [{"chord", "path": [[leg, ±1]], "holonomy", "deficit", "retention"?}],
 "log_price", "defect", "curvature",
 "flows": {"attraction", "net_inflow", "driven", "local", "global", "parts": {"driven", "local", "global"}},
 "pairing": {"total", "prices", "local", "global"}}
```

### 3.6 The Server

| Path | |
|---|---|
| `GET /api/engine` | whether the engine is built and current, and the limits of a request |
| `GET /api/data` | the products built, and the datasets of one's own |
| `POST /api/map` `{"year", "dataset"?}` | the engine's `map` |
| `POST /api/timeline` `{"dataset"?}` | the engine's `timeline` |
| `POST /api/spacetime` `{"first", "last", "currencies"?, "dataset"?}` | the engine's `spacetime` |
| `POST /api/production` `{"table", …}` | the engine's `production` |
| `POST /api/circuits` `{"network"}` | the engine's `circuits` |
| `POST /api/cancel` | stops the computations in progress |
| `GET /api/datasets`, `POST /api/datasets` `{"kind", "title", "files": [{name, text}], "base"?}`, `DELETE /api/datasets/<id>` | the datasets |
| `GET /api/plugins`, `GET /api/plugins/describe?file=` | the plug-ins and their descriptions |
| `POST /api/plugins/run` `{"file", "params", "dataset"?, "year"?}` | starts a run; `{"job"}` |
| `GET /api/plugins/job/<id>`, `POST /api/plugins/job/<id>/cancel` | a run's state, progress and results; stopping it |
| `GET /api/tables`, `GET /api/tables/<id>?rows=n` | the toolbox's tables, without or with their data |
| `POST /api/tool` `{"tool", "table"? \| "columns"?, "where"?, "roles", "params"}` | the engine's `tool`, on the columns the roles name, in the rows `where` keeps; trade's tools on the flows, gravity with the economies' GDP |
| `GET /setup`, `POST /setup/run` | the recipe: the engine and the data |

Answers for a year or a span are kept in `data/cache/economy-atlas/answers/`, keyed by the engine's build and the product.

## 4. Files

| Path | |
|---|---|
| `index.html`, `styles.css`, `src/` | The page: ES modules, no build step; its colours as tokens of day and night |
| `serve.py`, `engine.py`, `recipe.json` | The server, the engine's builder, the recipe |
| `engine/economy.cpp` | The engine |
| `engine/tools.hpp`, `engine/tools/` | The toolbox's engine: reading a request, writing an answer, and a header for each group of tools |
| `datasets.py` | One's own data, read into the library's schemas and kept |
| `tables.py` | The toolbox's tables, made from the library's products and the World Bank's indicators |
| `src/toolbox.js`, `src/tools.js`, `src/charts.js` | The toolbox's view, its catalogue of tools, and its charts |
| `plugins/` | The plug-ins, with their README and three examples |
| `examples/` | Invented data to try the formats |
| `src/map.js` | The map: places, arcs, eddies and the plug-ins' layers on Equal Earth; zoom, the card under the pointer, the strip of years; its instruments, the address it keeps, and its figures |
| `src/plugins.js`, `src/data.js`, `src/state.js` | The plug-ins' forms and runs, the dialog for one's own data, what the views share |
| `src/spacetime.js` | The currencies' sheet of lenses and their cumulative curvature |
| `src/production.js`, `src/circuits.js`, `src/presets.js` | The editable table and network, and their examples |
| `src/palette.js` | The data's colours, read from the theme |
| `tests/` | The application's checks, among them those of one's own data (`check_datasets.py`) and of the example plug-ins (`check_plugins.py`) |
