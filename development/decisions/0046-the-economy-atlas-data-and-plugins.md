# 0046. The Economy Atlas: One's Own Data, Plug-ins over Networks of Places, and a Ledger in Gold

- **Status:** Accepted, provisionally (the colours are the user's to confirm)
- **Date:** 2026-10-24
- **Amends:** [ADR 0044](0044-value-in-space-and-time.md)
- **Decided by:** the user asked to revise the Economy Atlas's interface ("we love weak light gold main colour, rose red as minor colour … please revise"), to let it use new data of one's own, and to give it plug-ins, "custom written scripts (support Python and C++ as plug-in) to analyse data and render in the economics GIS", and to make it more elegant and creative. Claude designed the palette with the dataviz method's validator, the formats of one's own data, the plug-in contract and the new map, and wrote this record. A subagent wrote the API's module and header and the example plug-ins to the contract stated here.

## 1. Context

The Economy Atlas drew only the data library's two products, in a provisional ochre. The Earth Atlas has plug-ins ([ADR 0027](0027-earth-atlas-plugins.md), [ADR 0029](0029-plugins-in-python-and-cpp.md)), but their contract is one of fields on grids: their needs are products and levels, their C++ file contract carries float32 grids, and their results are field products. A network of places is another kind of input: places with positions, flows between them, and the domain's decomposition of the flows.

## 2. Decision

### 2.1 A Ledger in Gold

The interface is weak light gold on ivory by day and on deep ink by night (the system's choice unless the reader chooses), with rose red as its minor colour. The data's colours were chosen by the job they do and validated for colour vision deficiencies on both surfaces (OKLab ΔE ≥ 8 between neighbours under protanopia and deuteranopia; text ≥ 4.5:1):

| Job | Day | Night |
|---|---|---|
| value sent out (diverging pole) | steel blue `#3A6EA5` | `#5E92D4` |
| value drawn in (diverging pole) | rose red `#C4506E` | `#CF5F7A` |
| driven flow | indigo `#4B4F8F` | `#7178C6` |
| circulation round the triangles, the eddies | rose `#C4506E` | `#CF5F7A` |
| circulation round the holes | gold `#A98A2C` | `#AE8D36` |
| a magnitude (sequential) | the neutral to deep gold | the neutral to bright gold |

The cool pole stays blue: rose against gold, both warm, failed for protanopia as the poles of a diverging map.

### 2.2 One's Own Data

`applications/economy-atlas/datasets.py` reads data given as CSV or JSON and keeps them as products of the library's own schemas, so that the engine reads them unchanged and every view and plug-in works on them:

- **Flows between places** (`trade-flows/1`): a file of flows, `from,to,value[,year]`, and for places that are not the library's economies a file of places, `id,name,latitude,longitude[,group]`. Ids that are ISO alpha-3 codes or names of the library's economies need no place, and every spelling of one economy becomes one place.
- **Currencies held and exchanged** (`currency-spacetime/1`): `month,currency,rate,policy`, the rate against a base the record names (or finds: the currency whose rate is always 1 or absent). The engine's spacetime takes the base from the record instead of the euro.

Every refusal says what to change: an unknown id, a negative flow, missing columns, a word for a number, two places at one point. The data are kept under `data/cache/economy-atlas/datasets/`, never committed.

### 2.3 Plug-ins over Networks of Places

A new module of the API, `api/python/science/networks.py`, and a header, `api/cpp/network_plugin.hpp`, beside the Earth Atlas's plug-in classes and reusing their compiler and results:

- **The input**: the network of a dataset in a year. It has places (id, name, latitude, longitude, group), the gross flows, and the engine's decomposition of them: the attraction of each place, the driven, local and global parts of each net leg, and the circulation round each triangle.
- **The results**:
  - `values`, a number for each place, drawn as its colour (diverging or sequential);
  - `arcs` between places;
  - the Earth Atlas's `sketch`, `table` and `note`.
- **One file contract for both languages**: `input.json` in, progress lines on standard output, `output.json` out. Python plug-ins run as processes too (`python -m science.networks run`), so a script cannot take the server down and every run can be stopped or timed out. A C++ plug-in is compiled with the economics domain's and the foundations' include directories, so it can use the domain itself.
- **The host** (`networks.Host`) runs only files in the application's `plugins/` directory. On a public address, running them needs the server's token, as the recipe does.

Three examples: Tinbergen's gravity equation (Python, least squares), the concentration of each economy's exports (Python), and where value settles if each economy passed its stock on in its export shares (C++, on `economics::stationary`).

### 2.4 The Map

- **The dock** holds the data (the library's or one's own, with **+ your data**), the layers (what colours the places, which arcs, how many, the triangles, the eddies, the names, the arcs flowing), and the plug-ins, each with a form built from its parameters, a run with progress and a stop.
- **The map** zooms and pans. It fits a regional dataset, draws places without a country of their own as circles, names the largest without overlap, and answers the pointer with a card for a place or an arc. Each arc wears the colour of its largest part.
- **A strip of years** beneath the map gives every year's shares, from a new engine operation, `timeline`, that decomposes every year of a product in one run.
- **Links**: the address holds the view and what it shows, so that every picture can be linked and taken again.

## 3. Alternatives Considered

- **The Earth Atlas's field contract** (`hosting.py`, `plugin.hpp`). Its needs, files and results are those of grids; a network fitted into it would be a grid of one row. The new module reuses its compiler and its result kinds instead.
- **Python plug-ins imported into the server**, as the Earth Atlas's are. Rejected for isolation, time limits and one contract for both languages, at the cost of 0.4 s of process start.
- **Plug-in code uploaded through the page.** Not done: plug-ins are files in `plugins/`, as the Earth Atlas's are, which keeps running code a matter of the local disk.
- **Gold against rose as the diverging poles.** Failed the validator for protanopia.

## 4. Consequences

- **The API** gains `science.networks` and `network_plugin.hpp`, with their checks (`api/tests/check_networks.py`), and §2.7 of its README.
- **The Economy Atlas** gains `datasets.py`, `plugins/`, `examples/` (invented data, to try the formats) and its checks (`tests/check_datasets.py`, `tests/check_plugins.py`). Its engine gains `timeline` and spacetime on any base. Its README tells the formats and the plug-in contract.
- **CLAUDE.md §7** records the Economy Atlas's colours; the gallery's pictures of the atlas are taken again in them.
