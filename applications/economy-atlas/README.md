# The Economy Atlas

*Status: draft.* An application of the economics domain ([ADR 0044](../../development/decisions/0044-value-in-space-and-time.md)). It draws value carried between places and through time:

- the world's trade on a map, split into exchange, drive and circulation;
- currencies held and exchanged as a spacetime complex, whose curvature is the departure from a price system;
- production as a circular flow;
- the loop deficits of transfer networks.

Everything shown is computed by the domain's C++ on the foundations. Serves research and learning.

## 1. What It Shows

- **The map.** The world's trade in goods in a year (the OECD's balanced bilateral statistics, 2010–2024) on an Equal Earth map, the geography domain's projection:
  - what two economies send each other both ways is reciprocal exchange; what is left, one net flow for each pair, is decomposed (the domain's theory note, §3) on a complex whose faces are the spherical Delaunay triangles of the economies' places;
  - each link conducts in proportion to its gross trade, so the *driven* part between two economies is their gross trade times their difference in attraction α, the potential from the economies that send value to those that draw it (a share of gross trade, without unit): countries are coloured by α;
  - *local circulation* turns round the triangles of neighbouring economies, drawn as eddies;
  - *global circulation* runs round the loops of distant partners no triangle fills;
  - arcs are the largest net flows, or their driven or circulating parts;
  - a country's card gives its exports, imports, attraction, the share of the parts on its links that circulates, and its largest net flows.
- **Spacetime.** Twenty currencies held month by month at their central banks' policy rates (BIS; the euro at the ECB's deposit facility rate, which euro money has earned since 2015) and exchanged with the euro at the month's average rate (ECB), as the spacetime complex of the theory note, §5:
  - a sheet shows every lens's curvature, the euro's excess return over a currency in a month, the departure from uncovered interest parity;
  - lines show its sum in time;
  - the panel measures the curvature against that of the exchange rates alone: the share of it the rates of holding remove (a gauge invariant; the share of the transports a price system would carry is not, and is not shown).
- **Production.** An input–output table, edited in place (Miller and Blair's example, and a three-sector circular economy), solved by the Leontief model:
  - technical coefficients, the inverse, multipliers, prices and linkages;
  - cost-push and demand shocks;
  - the deliveries split into exchange, drive and circulation.
- **Circuits.** A transfer network, edited in place, with the examples of the user's paper on loop deficits (its triangle, its five positions, a fee round a circuit, a market with prices, a round trip):
  - every loop's holonomy and deficit in the tree gauge;
  - the price potential where the deficits vanish;
  - retention round directed cycles;
  - the circulation of flows and their pairing with the transports.

## 2. Running

From the repository root, make the engine and the data once (g++ 11 or later, natively or in WSL; about 60 MB from the OECD, the ECB and the BIS), and serve the page:

```
python data/acquisition/recipe.py applications/economy-atlas/recipe.json
python applications/economy-atlas/serve.py --open
```

The page is at http://127.0.0.1:8910/. It offers to run the recipe itself when something is missing.

## 3. The Protocol

The engine is `engine/economy.cpp`, header-only C++20 on the economics domain (`domains/economics/`), physis (discrete exterior calculus, the spherical Delaunay triangulation), morphe (complexes, the worldline construction), calx (linear algebra) and the API's JSON. The server runs it once per request: one JSON request on its standard input, one JSON answer on its standard output, an error as `{"error": "..."}`. Products are passed by their paths; the engine reads nothing else.

### 3.1 Requests

| Request | Answer |
|---|---|
| `{"op": "map", "flows": "<trade-flows.json>", "year": y}` | the year's economies, net flows and the map's triangles (§3.2) |
| `{"op": "spacetime", "record": "<currency-spacetime.json>", "first": "YYYY-MM", "last": "YYYY-MM", "currencies"?: [...]}` | the lenses, statistics and price system of the currencies with complete records over the span (§3.3) |
| `{"op": "production", "table": {"sectors", "z", "final_demand", "value_added"}, "value_added_coefficients"?, "final_demand"?}` | the Leontief model and the deliveries' circulation (§3.4) |
| `{"op": "circuits", "network": {"positions", "legs": [{"from", "to", "ratio", "flow", "kind"}], "faces"?}}` | the network's cycles, prices and circulation (§3.5) |

### 3.2 The Answer to `map`

```
{"year", "years", "unit",
 "economies": [{"iso3", "name", "continent", "latitude", "longitude", "exports", "imports",
                "balance", "attraction", "circulating"}],
 "legs": [[from, to, net flow, driven, local, global]],              // net legs, millions of US dollars
 "triangles": [{"places": [i, j, k], "curl"}],
 "totals": {"total", "reciprocal", "net", "driven", "local", "global"},  // the last three: shares of ‖net‖²
 "faces": {"triangles", "neutral_legs", "covers_sphere"}, "method", "source"}
```

### 3.3 The Answer to `spacetime`

```
{"currencies", "names", "months", "times",
 "lenses": [[currency, month from, month to, curvature]],
 "stats": [{"currency", "mean", "volatility", "cumulative": [...]}],  // a year; the euro's excess return over the currency
 "presences": [[currency, month, log price]],                         // the exact part's potential
 "policy", "rates", "shares": {"exact", "coexact", "harmonic"},
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
 "flows": {"attraction", "net_inflow", "driven", "local", "global", "exact", "coexact", "harmonic"},
 "pairing": {"total", "prices", "local", "global"}}
```

### 3.6 The Server

| Path | |
|---|---|
| `GET /api/engine` | whether the engine is built and current, and the limits of a request |
| `GET /api/data` | the products built: the trade flows' years and economies, the currencies' months |
| `POST /api/map` `{"year"}` | the engine's `map` |
| `POST /api/spacetime` `{"first", "last", "currencies"?}` | the engine's `spacetime` |
| `POST /api/production` `{"table", …}` | the engine's `production` |
| `POST /api/circuits` `{"network"}` | the engine's `circuits` |
| `POST /api/cancel` | stops the computations in progress |
| `GET /setup`, `POST /setup/run` | the recipe: the engine and the data |

Answers for a year or a span are kept in `data/cache/economy-atlas/answers/`, keyed by the engine's build and the product.

## 4. Files

| Path | |
|---|---|
| `index.html`, `styles.css`, `src/` | The page: ES modules, no build step |
| `serve.py`, `engine.py`, `recipe.json` | The server, the engine's builder, the recipe |
| `engine/economy.cpp` | The engine |
| `src/map.js` | The map: countries, arcs and eddies on Equal Earth |
| `src/spacetime.js` | The currencies' sheet of lenses and their cumulative curvature |
| `src/production.js`, `src/circuits.js`, `src/presets.js` | The editable table and network, and their examples |
| `tests/` | The application's checks |

Its colour is a warm ochre (#A8743A), provisionally, until the user chooses one.
