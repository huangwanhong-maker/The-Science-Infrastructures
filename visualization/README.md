# visualization: Representation Layer

*Status: planned.* Edge. For the design, see [architecture §5](../docs/architecture.md) and [ADR 0005](../development/decisions/0005-visualization-as-representation-layer.md).

## 1. Purpose

This layer presents results to people, using:

- scientific figures;
- 3D scenes of the planets and the universe, built from real data;
- interactive views;
- techniques that do not exist yet.

Only `applications/`, `notebooks/` and `research/` may use this layer. Nothing in `core/`, `data/` or `observation/` depends on it.

## 2. Structure

| Directory | Role |
|---|---|
| `scene/` | A device-neutral, serialisable description of what is shown. It carries units, time scales and reference frames. |
| `renderers/` | One renderer per technique |
| `styles/` | Palette, typography and marks |

## 3. Planned Renderers

| Renderer | Technique |
|---|---|
| `static` | Matplotlib, SVG, PDF |
| `web` | WebGPU, three.js, CesiumJS |
| `native` | Vulkan |
| `immersive` | WebXR |
| `sonification` | Sound |

## 4. Style

The visual style is warm and humane, not that of a generic dashboard. The palette is defined in `styles/`, together with the first renderer.

## 5. Prior Art

Study these before building: OpenSpace, NASA Eyes, CesiumJS, Stellarium and Celestia.
