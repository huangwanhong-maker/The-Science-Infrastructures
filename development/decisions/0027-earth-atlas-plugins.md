# 0027. Plug-ins of the Earth Atlas: Algorithms Loaded on Request, Shown Through the Layers

- **Status:** Accepted, provisionally: designed and built at the user's request; open to revision on review. Extends [ADR 0014](0014-earth-atlas-separated-data-and-typed-layers.md) and [ADR 0015](0015-earth-atlas-models-imagery-detail.md) §2.3.
- **Date:** 2026-10-22
- **Decided by:** the user asked for a plug-in design for weather prediction, data assimilation and other algorithms, with scripts that expose functions of well-designed classes, loaded dynamically and visualised on the Earth. Claude chose the design below.

## 1. Context

The Earth Atlas draws data typed by shape through layers typed the same way (ADR 0014), and it runs two models of its own in the page (ADR 0015, [ADR 0019](0019-shallow-water-model.md)). ADR 0015 set the rule for anything computed in the page: *a computed product has the same interface as a library product, and is drawn, contoured, compared and verified like any other.* Those two models, however, are wired by hand: the panel names them, `app.js` registers their output and selects it, and a third model would mean editing both.

The user wants algorithms to be written as scripts that the atlas loads when asked, runs on the data it holds, and shows on the Earth: weather prediction, data assimilation, diagnostics and visualisation tasks. The scripts should derive from well-designed classes and expose specific functions, and examples should show how.

## 2. Decision

1. **A plug-in is an ES module that exports `create(api)`.** The function receives the atlas's plug-in API and returns an instance of one of its classes. The host loads the module by dynamic `import()`, from the atlas's own `plugins/` directory (listed in `plugins/index.json`) or from any URL the user types, and never without a click. A script written this way needs no import of the atlas; it may import the domains and foundations by URL for the science it uses.
2. **Classes named by what the algorithm is.** `Plugin` holds identity, parameters, the fields the plug-in needs, and `run()`. Four subclasses write `run()` from smaller hooks that a script overrides:
   - `Diagnostic`: a field computed from fields, at every time of the input, on its grid (`outputs()`, `compute()`);
   - `ForecastModel`: a state advanced in time from an analysis, recorded every few hours and verified against the input's later times beside persistence (`grid()`, `initial()`, `advance()`, `record()`, `verify()`);
   - `DataAssimilation`: a background and observations to an analysis and its increment (`observations()`, `analyse()`);
   - `Trajectories`: parcels carried through a flow by a Runge–Kutta step on the sphere (`seeds()`, `velocity()`).
3. **Results are data, drawn by the atlas.** A plug-in returns a list of results and never draws:
   - a *field product* (`ComputedField`) joins the library and the catalogue, so it is a model beside the others under the same quantity: drawn, contoured, differenced with "minus", carried by particles, read in the place card, and listed in the verification table when it carries one;
   - a *sketch*, polylines, markers, labels and moving points in longitude and latitude, is drawn as it is by a new `SketchLayer`, the one place that turns such data into geometry;
   - *tables* and *notes* appear in the panel.
4. **Needs, resolved by the host.** A plug-in declares the fields it needs by standard name and level. The host finds every product and level that carries them, the panel lets the user choose the source model for each need, and `run()` receives the chosen products. The context a plug-in sees is narrow: the instant shown, the chosen place, the library's field metadata, and a progress hook.
5. **Where things live.**
   - The framework in `applications/earth-atlas/src/plugins/` (`plugin.js`, the classes; `context.js`; `host.js`; `index.js`, the API) and `src/layers/sketch.js`.
   - The built-in plug-ins, which are also the examples, in `applications/earth-atlas/plugins/`, with a README that is the guide to writing one. They obey the dependency rules; the checker scans `plugins/` as it scans `src/`.
6. **Checked without a page.** The classes touch neither the DOM nor three.js, so a plug-in runs in Node from synthetic fields. Each example has a check with an analytic answer: the vorticity of a solid-body rotation, a parcel in a uniform zonal wind, an analysis that reproduces exact observations and leaves the far field alone, a steady zonal flow that the barotropic model keeps.

## 3. Alternatives Considered

- **Plug-ins as layers, with three.js and the stage in hand.** The most powerful, and the most coupling: every script would depend on the rendering internals, and the rule that data are typed by shape and drawn by type would be lost. The sketch covers the visualisation that products do not.
- **A declarative recipe without code.** Fine for choosing and styling data; too weak for an assimilation or a model.
- **Python plug-ins run by the server.** The data library already runs jobs for the page (requests on demand, recipes), and heavy algorithms belong there eventually. The user asked for algorithms in the page first; the same classes could be mirrored in Python later, with the server as a `Context`.
- **Isolation in an iframe or a worker with a message protocol.** Safer, and much heavier; a plug-in would lose direct access to the products. The decision is that a plug-in runs with the page's authority, and that nothing loads without the user's click. Scripts from elsewhere are to be read before they are loaded, as any code is.

## 4. Consequences

- **Four examples** in `plugins/`: relative vorticity and divergence of the wind (`Diagnostic`), parcels through the forecast wind (`Trajectories`), optimal interpolation of a background toward a network of observations (`DataAssimilation`), and the barotropic vorticity model run through the class (`ForecastModel`).
- **The two models of ADR 0015 and 0019** are now the pattern's second appearance; migrating them to `ForecastModel` plug-ins would delete the hand-written wiring. It is the next step, not taken here so that this change stays one concept.
- **The panel** gains a Plug-ins section: the loaded plug-ins with their sources, parameters, progress, tables and sketches; the built-in ones to load; a field for a URL.
- **Security.** A loaded script runs as the page does. The page holds no secrets, and the local server exposes only reads and the data recipe, but a script can do what the page can. The README says so.
- **`window.atlas.plugins`** exposes the host for inspection and automated checks.
