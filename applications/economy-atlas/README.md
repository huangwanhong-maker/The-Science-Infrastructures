# The Economy Atlas

*Status: draft.* An application of the economics domain ([ADR 0044](../../development/decisions/0044-value-in-space-and-time.md), [ADR 0046](../../development/decisions/0046-the-economy-atlas-data-and-plugins.md)). It draws value carried between places and through time:

- the world's trade on a map, or flows of one's own between any places, split into exchange, drive and circulation;
- currencies held and exchanged as a spacetime complex, whose curvature is the departure from a price system;
- production as a circular flow;
- the loop deficits of transfer networks.

Everything shown is computed by the domain's C++ on the foundations. Plug-ins in Python and C++ analyse the flows the map shows and draw their results on it. Serves research, learning and use.

## 1. What It Shows

- **The map.** Flows between places in a year on an Equal Earth map, the geography domain's projection: by default the world's trade in goods (the OECD's balanced bilateral statistics, 2010–2024), or data of one's own (§2.2).
  - What two places send each other both ways is reciprocal exchange. What is left, one net flow for each pair, is decomposed (the domain's theory note, §3) on a complex whose faces are the spherical Delaunay triangles of the places.
  - Each link conducts in proportion to its gross flow. The *driven* part between two places is their gross flow times their difference in attraction α: the potential runs from the places that send value to those that draw it, a share of gross flow without unit. Places are coloured by α, or by the values a plug-in returned.
  - *Local circulation* turns round the triangles of neighbours, drawn as rose eddies. *Global circulation* runs round the loops of distant partners no triangle fills.
  - Arcs are the largest net flows, their driven or circulating parts, or a plug-in's arcs. Each wears the colour of its largest part: indigo driven, rose round the triangles, gold round the holes.
  - The map zooms (the wheel, the buttons) and pans (dragging), and fits a regional dataset. Places without a country of their own are circles. The pointer brings a card for a place or an arc. Clicking a place opens its card: its flows, attraction and largest net flows, each with the share of it driven and circulating.
  - A strip of years beneath the map gives every year's shares of the net flows' energy; a click, or ← and →, chooses the year.
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

The page is a ledger in gold: weak light gold on ivory by day and on deep ink by night, as the system prefers or as chosen (the button at the top right), with rose red as its minor colour. Its data's colours are validated for colour vision deficiencies on both (ADR 0046, §2.1).

## 2. Running

From the repository root, make the engine and the data once (g++ 11 or later, natively or in WSL; about 60 MB from the OECD, the ECB and the BIS), and serve the page:

```
python data/acquisition/recipe.py applications/economy-atlas/recipe.json
python applications/economy-atlas/serve.py --open
```

The page is at http://127.0.0.1:8910/. It offers to run the recipe itself when something is missing.

### 2.1 Links

The address holds the view and what it shows, so that a picture can be linked and taken again:

- `#map&year=2010&arcs=driven&economy=CHN`:
  - `year`;
  - `arcs`: `all`, `driven`, `circulating` or `none`;
  - `shown`: how many arcs;
  - `faces=0`, `labels=0`, `flowing=0`: hide the triangles, the names or the movement of the arcs;
  - `economy`: a place's id, the ISO alpha-3 code for the library's;
  - `dataset`: a dataset of one's own, by its id.
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

## 3. The Protocol

The engine is `engine/economy.cpp`, header-only C++20 on:

- the economics domain (`domains/economics/`);
- physis: discrete exterior calculus, flows with conductances, the spherical Delaunay triangulation;
- morphe: complexes, the worldline construction;
- calx: linear algebra;
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
| `GET /setup`, `POST /setup/run` | the recipe: the engine and the data |

Answers for a year or a span are kept in `data/cache/economy-atlas/answers/`, keyed by the engine's build and the product.

## 4. Files

| Path | |
|---|---|
| `index.html`, `styles.css`, `src/` | The page: ES modules, no build step; its colours as tokens of day and night |
| `serve.py`, `engine.py`, `recipe.json` | The server, the engine's builder, the recipe |
| `engine/economy.cpp` | The engine |
| `datasets.py` | One's own data, read into the library's schemas and kept |
| `plugins/` | The plug-ins, with their README and three examples |
| `examples/` | Invented data to try the formats |
| `src/map.js` | The map: places, arcs, eddies and the plug-ins' layers on Equal Earth; zoom, the card under the pointer, the strip of years |
| `src/plugins.js`, `src/data.js`, `src/state.js` | The plug-ins' forms and runs, the dialog for one's own data, what the views share |
| `src/spacetime.js` | The currencies' sheet of lenses and their cumulative curvature |
| `src/production.js`, `src/circuits.js`, `src/presets.js` | The editable table and network, and their examples |
| `src/palette.js` | The data's colours, read from the theme |
| `tests/` | The application's checks, among them those of one's own data (`check_datasets.py`) and of the example plug-ins (`check_plugins.py`) |
