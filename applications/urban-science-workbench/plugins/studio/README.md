# The studio

*Your own data in space and time, and your own code over it.* Status: **draft**. Activity: research. [ADR 0079](../../../../development/decisions/0079-data-in-space-and-time.md).

## 1. Purpose

Anything with a place and a time goes on the workbench's map, and runs through time above the clock. A script is given what is shown, and answers what to draw.

## 2. Data

- **A file**, dropped anywhere on the page, which opens the studio to take it: a table with a place and a time (CSV, TSV, JSON lines), GeoJSON, or a record of the contract. Its columns are told by their names and contents, and asked for only if they cannot be. What cannot be placed is left out, and counted by reason on its card.
- **A record of the library:** earthquakes since 1976 (Global CMT), earthquakes this month (USGS), tropical cyclones (IBTrACS), and in Melbourne its walkers through a mean week.

## 3. Code

One contract in three languages, set out in [scripts/README.md](scripts/README.md), with an example in each:

| Example | Language | What it does |
|---|---|---|
| [`density.py`](scripts/density.py) | Python | a Gaussian kernel density of events, through time, on the ground |
| [`hex_bins.cpp`](scripts/hex_bins.cpp) | C++ | events counted on physis's H3 cells, through time, as columns |
| [`track_speed.js`](scripts/track_speed.js) | JavaScript | the speed along every track, with the fastest and a distribution |

## 4. For a Plug-in

A plug-in of the workbench reaches the same line of time from its layers by `ctx.st`:

- `ctx.st.show(record, style)` draws a record, and answers a handle with `remove()` and `setStyle()`;
- `ctx.st.records()`, `ctx.st.now()` and `ctx.st.window()` read what is shown and the moment;
- `ctx.st.timeline` plays, pauses and seeks;
- `ctx.st.acceptFiles(take)` receives the files dropped on the page, and answers what to call when it stops.

`style` may name `mark`, `value`, `map`, `scale` (`linear`, `sqrt`, `log` or `diverging`), `opacity`, `window_ms`, `credit`, `top` (the height of the tallest column, m), `zone` (the IANA zone its moments are read in) and `weekly` (a nominal week, read by its days).

## 5. Limits

- A record is sent to a server script as JSON, which is slow beyond a few hundred thousand rows.
- A record shown is not kept between sessions, though a script kept is.
- Only one raster can lie on the ground at a time.
