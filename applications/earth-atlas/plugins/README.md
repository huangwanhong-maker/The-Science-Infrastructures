# plugins: Algorithms Loaded on Request

*Status: draft (2026-10-22).* Decision record: [ADR 0027](../../../development/decisions/0027-earth-atlas-plugins.md). The classes are in [src/plugins/plugin.js](../src/plugins/plugin.js); the host that loads and runs plug-ins is [src/plugins/host.js](../src/plugins/host.js).

A plug-in is a script that the Earth Atlas loads when asked, runs on the data it holds, and shows on the Earth. It derives from one of the atlas's classes, says which fields it needs, and returns data, which the atlas draws like anything else. Nothing is loaded without a click, and a loaded script runs with the page's authority: read a script before loading it, as any code.

## 1. The Contract

A plug-in is an ES module with one export:

```js
export function create(api) {
  class Mine extends api.Diagnostic { /* ... */ }
  return new Mine();
}
```

`api` is [src/plugins/index.js](../src/plugins/index.js): the classes, the result makers and a few numerical helpers. A script written against it needs no import of the atlas; it may import the domains and foundations by URL for the science it uses, as the examples here do with relative paths.

## 2. The Classes

Each class is named by what the algorithm is, and writes `run()` from smaller hooks that a subclass overrides. All of them share:

| Member | Meaning |
|---|---|
| `id`, `title`, `description`, `domain` | Identity. `id` names the product (`plugin-<id>`); a second run replaces the first. |
| `parameters` | `{ name: { label, default, min, max, step, unit } }` or `{ name: { label, default, choices: [[value, label]] } }`. The panel offers them and passes their values to `run()` as `params`. |
| `needs()` | The fields needed: `[{ name, label, standard_name` or `standard_names: [u, v], level?, times?, distinct? }]`. The host finds every product and level that carries them; the user chooses the source; `run()` receives `inputs[name] = { product, names, variables, meta, model, level }`. `distinct` asks for a source other than the needs before it; `times` for at least that many times. |
| `run(inputs, params, context)` | The algorithm, returning results. The subclasses below write it from their hooks. |
| `breathe(context, fraction, text)` | A resting point in a long loop: progress for the panel, and a turn for the page. |

| Class | What it is | Hooks to override |
|---|---|---|
| `Diagnostic` | A field computed from fields, at every time of the first need, on its grid | `outputs()`: the variables returned; `compute(inputs, params, index, context)`: `{ name: Float32Array }` for one time; `vectors()`, `show()` |
| `ForecastModel` | A state advanced in time from an analysis, started at the input time at or before the instant shown, recorded every few hours, verified against the input's later times beside persistence (20–80° N) | `grid()`, `outputs()`, `initial(inputs, index, params, context)`, `advance(state, seconds, params)`, `record(state, params)`, `verify()`; the fields `dt`, `every` |
| `DataAssimilation` | A background and observations to an analysis and its increment | `observations(inputs, index, params, context)`: `[{ lon, lat, value }]`; `analyse(background, grid, observations, params, context)`: the analysis; `markers()` |
| `Trajectories` | Parcels carried through a flow by a Runge–Kutta step on the sphere | `seeds(inputs, params, context)`: `[{ lon, lat, colour? }]`; `velocity(inputs, params)`: `(lon, lat, ms) -> [u, v]` or `null` |

`Plugin` itself may be derived from for anything else, with `run()` written in full.

## 3. Results

`run()` returns a list of results made by the makers in `api`. They are data; the atlas does the drawing.

| Maker | Result | Shown as |
|---|---|---|
| `field(product, show)` | A field product (`makeProduct`), with `show = { surface, flow, lines }` naming what to select by variable or vector name | A model among the others of its quantity: drawn, contoured, differenced with "minus", carried by particles, read in the place card; listed in the verification table when it carries one |
| `sketch({ title, polylines, markers, labels, moving })` | Polylines `{ coords: [lon, lat, ...], colour, width, alpha, dashed }`, markers `{ lon, lat, colour, size }`, labels `{ lon, lat, text }`, moving points `{ times: [ms], coords, colour, size }` | Drawn as they are over the globe, in every projection; a moving point is drawn where it is at the time shown |
| `table(title, columns, rows, note)` | Headings and rows | In the panel, under the plug-in |
| `note(text)` | A sentence | In the panel |

`makeProduct({ name, title, model, domain, grid, times, variables, vectors, source, verification })` builds a product from `variables: [{ name, long_name, standard_name, units, level?, slices: [Float32Array per time] }]`; `grid` is `{ lon0, dlon, nlon, lat0, dlat, nlat, periodic }` with rows north to south, values row-major, `NaN` where missing. The catalogue recognises standard names and levels, so a product with `geopotential_height` at `500 hPa` joins the height maps of the forecast models. A variable that is a difference of two quantities of its units, such as an analysis increment in kelvin, carries `difference: true`: it is then drawn on a symmetric diverging scale, and the offset of its units (kelvin to degrees Celsius) is not applied.

## 4. The Context

`context` is what a plug-in sees of the atlas while it runs: `time`, the instant shown in ms UTC; `pinned`, the chosen place as `[lon, lat]` or `null`; `fields()`, the metadata of every field product; `field(name)`; `geography()`; and `progress(fraction, text)`. The same view is made from files when a plug-in is checked without a page.

## 5. Helpers in `api`

`sampleOnGrid(product, name, index, grid)` samples a variable onto another grid; `indexAtOrBefore(times, ms)`; `lattice(n)` spreads `n` points evenly over the sphere; `parcelStep(flow, lon, lat, ms, dt)` is one Runge–Kutta step on the sphere; `iso(ms)`.

## 6. Loading

The panel's Plug-ins section lists the built-in plug-ins of [index.json](index.json) to load with a click, and takes the URL of any other script. A loaded plug-in shows its sources (one choice per need), its parameters, a Run button, its progress, its tables and the titles of its sketches; unload removes its sketches. From the console, `atlas.plugins.load('plugins/vorticity.js')` and `atlas.plugins.run('vorticity')` do the same.

## 7. Checking Without a Page

The classes touch neither the page nor three.js. A plug-in runs in Node from synthetic fields:

```js
const api = await import('../src/plugins/index.js');
const plugin = (await import('../plugins/vorticity.js')).create(api);
const product = api.makeProduct({ name: 'zonal', title: 'Zonal flow', grid, times, variables: [...] });
const [result] = await plugin.run({ wind: { product, names: ['u', 'v'], variables: product.meta.variables, meta: product.meta, model: 'test' } }, {}, context);
```

[tests/run-tests.mjs](../tests/run-tests.mjs) checks each example against an answer known in closed form.

## 8. The Examples

| File | Class | Shows | Checked by |
|---|---|---|---|
| [vorticity.js](vorticity.js) | `Diagnostic` | Relative vorticity and divergence of the wind at the level chosen, by centred differences | A solid-body rotation: ζ = 2U sin φ / a, δ = 0 |
| [trajectories.js](trajectories.js) | `Trajectories` | Parcels from a ring around the chosen place, or a lattice, through the forecast wind; the paths, and the parcels moving with the time shown | A parcel in a solid-body rotation keeps its latitude and turns by U t / a |
| [optimal-interpolation.js](optimal-interpolation.js) | `DataAssimilation` | A background analysed toward a second model's values at a lattice of points; the increment, and the analysis among the models of the quantity | Exact observations are reproduced where they are, and the background kept far from them |
| [barotropic.js](barotropic.js) | `ForecastModel` | The barotropic vorticity equation from the 500 hPa wind, the height in balance, verified against the source | A solid-body rotation is steady |

## 9. Limits and Next Steps

- A plug-in runs in the page's thread; `breathe()` keeps the page alive but a long model is slow there. A plug-in may start its own worker, as the atlas's models do.
- The two models the atlas runs itself (ADR 0015, 0019) are the pattern's other appearance; migrating them to `ForecastModel` plug-ins would delete their hand-written wiring.
- The classes could be mirrored in Python, with the server as the context, for algorithms too heavy for the page.
