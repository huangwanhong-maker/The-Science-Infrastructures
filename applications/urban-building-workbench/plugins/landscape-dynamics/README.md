# Landscape Dynamics

*A district as cells, and what happens in it.* A plug-in of the [Urban Building Workbench](../../README.md). Status: **draft** (2026-10-31). [ADR 0062](../../../../development/decisions/0062-landscapes-as-cellular-configurations.md).

A district of Melbourne's Hoddle Grid (Elizabeth to Russell Street, Flinders to Lonsdale Street) is taken as a cell complex: its junctions, the segments of its streets and lanes, and the blocks they enclose. Its structure, its use and how it is read change one another in one process. The urban-science domain states the formalism exactly ([theory/landscape-dynamics.md](../../../../domains/urban-science/theory/landscape-dynamics.md)); this plug-in draws it on the map and beside it, in plain terms.

## 1. Layers and Tools

| Layer | What it shows |
|---|---|
| `layers/cells.js`: junctions, streets and blocks | the district as a regular complex; blocks with a dead end left as holes; the two stores and where a lane between them would run |
| `layers/landscape.js`: the landscape running (a model) | the day at the clock's hour: walkers each way along the streets, an arrow where more walk one way, the blocks' warmth, the stores open or shut, the lane when it stands |
| `layers/spacetime.js`: the district's spacetime (a model) | the day with time as height: each hour's streets coloured by their traffic, the stores' worldlines thick while open, the lane's sheet from the instant it is built, deep once recognised, the blocks outlined where the cells change, the events as beads joined by their causal order, a plane at the clock's hour |
| `layers/rooms.js`: a building's rooms (a model) | a building picked on the map as rooms on its floors, coloured by their temperature at the clock's hour; doors, stairs and merged rooms as they open |

| Tool | Kind of result |
|---|---|
| The landscape and its events | proved (the cells, the lane's three rewrites and their order); a model (the day's run, its record, what could happen at each hour and what did) |
| The landscape's events, exactly | proved (composition, the law and irreversibility, what coarser descriptions forget and for how long, interventions' interactions); a model (evidence from a record) |
| A building's rooms | proved (the complex, its modes, the heat each event keeps); a model (the day's heat under the Sun and a declared air) |
| Scales and memory | proved (the district's modes, its flows round the holes, particles and flow one description); a model (what the streets remember of the shops) |

Every rate of the models is declared, not measured. They explore how structure, use and meaning can change one another; they do not predict a particular district.

## 2. The Records

- the walking network, the buildings and their levels: OpenStreetMap, through the data library's city model;
- the holdings and their uses: Street Oceanography's holdings (OpenStreetMap), the two stores the nearest pair of holdings a passer-by walks into that a lane of 15 m or more could join, or two picked on the map;
- the walkers entering the district each hour: half the City of Melbourne's counts at the sensors inside the box, on the clock's weekday;
- the Sun: the astronomy domain, at the clock's date.

## 3. The Engine

`engine/main.cpp`, C++20 on the urban-science domain's `landscape`, the architecture domain's `space` and `environment`, the foundations and the API's JSON, built by the workbench's `engine.py` into `data/cache/engines/urban/landscape-dynamics/`. Its operations: `district`, `pairs`, `lane`, `run`, `skeleton`, `building`, `scales` (`main.cpp` documents each request).

## 4. Limits

- The district is cut by a box: some of its pieces are cut off from the rest, and the footways beside streets are left out so that no street counts twice.
- Walkers hop junction to junction at walking pace; nobody plans a route.
- A building's rooms are a regular grid fitted to its footprint's smallest rectangle, not its plan; its fabric is typical, not measured.
- Shops' stays and the rate of entering them are declared by use.
