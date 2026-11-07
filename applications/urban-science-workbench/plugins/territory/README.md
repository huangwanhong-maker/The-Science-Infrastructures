# The territory

*What is known of the world's places, on one grid of hexagons and on the ground.* Status: **draft**. Activity: research. [ADR 0078](../../../../development/decisions/0078-the-whole-earth-in-the-urban-workbench.md).

## 1. Purpose

Each scale asks its own question, and this plug-in answers with what means something at that scale: the world's people and economies from orbit, a region's cells from a few hundred kilometres, a city's urbanisation, settlement and floods on the ground, a city's built volume by hexagon. All of it lies on physis's H3 where a grid is wanted, so that a population and a building can be compared where they stand.

## 2. Layers

| Lens | Layers | Source |
|---|---|---|
| People | where people live; poverty at $3 a day; inequality; people, 2015–2030; urban centres and suburbs; degree of urbanisation | GHS-POP and GHS-SMOD (JRC), the Poverty and Inequality Platform, Space2Stats |
| Change | settlement spreading, 1985–2015; night lights, 2012–2024; built-up land, 1975–2030; the Earth by night | WSF Evolution (DLR), Space2Stats, Black Marble (NASA) |
| Risk | river floods by return period; people in flood zones; cyclones; landslides; drought; wildfires | GloFAS (CEMS), Space2Stats |
| Fabric | built volume by hexagon | whatever buildings are in view |

## 3. Tools

- **Any indicator on the globe:** Data360 searched by meaning, an indicator raised over the economies.
- **A composite of place:** one's own index of a region, by nomos's composite. It is a description by one's weights, and not a measure.
- **Sources and what stands behind them:** what was taken from each tool and dataset the workbench was asked to draw on, what was left and why.

## 4. Engine

`engine/main.cpp`, on physis's H3 and nomos's composite: `earth_people`, `bin_points`, `composite`. Regional cells come from the data library's `services.py`, a 4 degree tile at a time, kept once asked.

## 5. Limits

- An economy's value says nothing of the cells within it: the hexagons only show where its people live.
- A cell holding any part of an area in dispute (Natural Earth's list) takes no economy's value or name, whichever holds it, under the naming policy of [ADR 0069](../../../../development/decisions/0069-names-of-places-and-a-2d-gis.md). Its people are still counted where people are drawn.
- A survey's year differs from economy to economy, and income and consumption are not strictly alike.
- River floods leave out the sea and small streams.
