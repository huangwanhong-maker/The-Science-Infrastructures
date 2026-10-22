# visualization: Representation Layer

*Status: draft.* Edge, shared by every application. For the design, see [architecture §5](../docs/architecture.md) and [ADR 0005](../development/decisions/0005-visualization-as-representation-layer.md).

## 1. Purpose

This layer presents results to people, using:

- scientific figures;
- 3D scenes of the planets and the universe, built from real data;
- interactive views;
- techniques that do not exist yet.

Only `applications/`, `notebooks/` and `research/` may use this layer. Nothing in `foundations/`, `domains/`, `data/` or `observation/` depends on it. It may itself use the domains its `module.json` declares (today geodesy), through their public interfaces, since units, time scales and frames travel with the data into the scene.

## 2. Structure

| Directory | Role |
|---|---|
| `scene/` | A device-neutral, serialisable description of what is shown. It carries units, time scales and reference frames. |
| `renderers/` | One renderer per technique |
| `styles/` | Palette, typography and marks |

## 3. Provisional JavaScript

Until the scene model and the renderers exist, the parts of the two atlases' drawing that any application could use live here, in `js/`, behind the public interface `js/index.js` ([ADR 0026](../development/decisions/0026-domains-as-libraries.md)):

| File | Content |
|---|---|
| `contour.js` | Isolines of a gridded field by marching squares, and its highs and lows |
| `flow.js` | Particles carried by a vector field on the sphere, to make a flow visible ([note](theory/particles-on-the-sphere.md)) |
| `units.js` | The units in which quantities are shown, and numbers written for reading |
| `tiles.js` | The tiles of an equirectangular map service that cover a region at a resolution |

Their checks run with `node visualization/tests/run-tests.mjs`.

## 4. Planned Renderers

| Renderer | Technique |
|---|---|
| `static` | Matplotlib, SVG, PDF |
| `web` | WebGPU, three.js, CesiumJS |
| `native` | Vulkan |
| `immersive` | WebXR |
| `sonification` | Sound |

## 5. Style

The visual style is warm and humane, not that of a generic dashboard. The palette is defined in `styles/`, together with the first renderer.

## 6. Prior Art

Study these before building: OpenSpace, NASA Eyes, CesiumJS, Stellarium and Celestia.
