# Reach and use

*How far one can go from any place, and how a city's walkers use its places.* Status: **draft**. Activity: research. [ADR 0078](../../../../development/decisions/0078-the-whole-earth-in-the-urban-workbench.md).

## 1. Layers

- **How far in so many minutes** (anywhere). Pick a place on the map: the minutes to every way about it run down the streets as a wave, on foot or by road, with the length of way reached within each time. On the ways of the tiles in view, so that it ends where they end.
- **Where the walkers go, against their habit** (Melbourne). The Urban Space Usage Index on H3 cells: each cell's share of the day's walkers at the counters, and its departure from the same weekday over the twelve weeks before. It follows the clock's date. A date in a gap of the record takes the nearest before it.

## 2. Engine

`engine/main.cpp`, on urban-science's `access` and `usage`, morphe's shortest paths and physis's H3: `reach`, `usage`.

## 3. Limits

- Speeds are declared in `layers/reach.js`, and are assumptions: walking at 1.34 m/s (Weidmann 1993), a road's free flow by its class. No signals, gradients, one-way streets or traffic: a time is a lower bound on the journey.
- Ways that end within 8 m of another point are joined, since a tile cuts a way at its edge. A dead end that stops that near another way is joined too.
- The counters are a sample of the centre's footpaths and not of the city, and the shares are of those reporting.
