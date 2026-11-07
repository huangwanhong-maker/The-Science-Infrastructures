# scripts: Your Code over Data in Space and Time

*Status: draft.* Scripts that the workbench's studio runs on what is shown ([ADR 0079](../../../../../development/decisions/0079-data-in-space-and-time.md)). The examples here are the guide: [`density.py`](density.py), [`hex_bins.cpp`](hex_bins.cpp) and [`track_speed.js`](track_speed.js). A script kept from the studio is written here.

## 1. The Contract

One contract in three languages: Python ([`api/python/science/spacetime.py`](../../../../../api/python/science/spacetime.py)), C++ ([`api/cpp/spacetime_plugin.hpp`](../../../../../api/cpp/spacetime_plugin.hpp)) and JavaScript ([`visualization/js/gis/spacetime.js`](../../../../../visualization/js/gis/spacetime.js)).

A script is **given** the records shown, its parameters defaulted and typed, and the context of the view: where it looks, the moment and its window, the place open.

It **answers** a list of results: records, which go on the map and on the line of time; and tables, curves and notes, which go in the panel.

Space is degrees on WGS 84, longitude east and latitude north. Time is milliseconds since 1970-01-01 in UTC.

| Record | What it is | Its columns | Drawn as |
|---|---|---|---|
| `events` | points in space and time | `time, lat, lon, values, labels?` | blooms, or columns |
| `tracks` | moving things | `tracks: [{id, name?, time, lat, lon, values?}]` | threads |
| `cells` | places with a series | `time, cells: [{id, ring, name?}], values: {name: [[each cell] at each time]}` | columns rising and falling |
| `grid` | a raster through time | `time, box, width, height, values: {name: [[width × height, rows from the south] at each time]}` | on the ground |
| `flows` | what passes between places | `places: [{id, lat, lon}], flows: [{from, to, time, until?, value?}]` | threads rising between them |

A record may carry `units` and `note`, and `show`, how it asks to be drawn: `{mark, value, map, scale, window_ms}`, where `scale` is `linear`, `sqrt`, `log` or `diverging`, chosen by the data if not said.

## 2. How Each Language Runs

- **JavaScript** runs in a worker of the page, with a time limit. It needs nothing of the server.
- **Python** runs as a process of the server, with numpy, with a time limit.
- **C++** is compiled on first use (natively, or by g++ in WSL), against `api/cpp`, the foundations and the domains, so that it may use any of them. It then runs as a process with a time limit.

## 3. Who May Run What

A Python or C++ script runs with the authority of whoever started the server. So the server runs only files in this directory, and writes here or runs them only when asked from its own page, on its own address: a request from any other page is refused. On a public address it refuses unless started with `--token`.

## 4. A Plug-in of the Workbench

A plug-in directory (`plugins/<name>/`) reaches the same line of time from its layers, in a few lines:

```js
export async function layer(ctx) {
  const record = { kind: 'events', title: 'Mine', time, lat, lon, values: { depth } };
  const shown = ctx.st.show(record, { value: 'depth', map: 'plum' });
  return { dispose() { shown.remove(); } };
}
```
