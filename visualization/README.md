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
| `isosurface.js` | A surface at a level of a volume (a brain, a region), as a triangle mesh, by surface nets |
| `complex.js` | Cell complexes drawn: a CW complex realised in space by stress majorisation, its 2-cells unfolded into a net, its cells as a Hasse diagram; surface words and facets read into Complex JSON; coloured graphs (the Feynman graphs of tensor models and group field theories) drawn in the plane and turned into their 2-complexes; a tetrahedron set in space from its six edge lengths (the quantum tetrahedra of spin foams); warm cyclic, sequential and diverging colour maps, and a colour for each index of a coloured graph |
| `stars.js` | The colour of a star from its colour index (B−V or Gaia's BP−RP) through its temperature (Ballesteros 2012; Mucciarelli and Bellazzini 2020) and a blackbody's sRGB (Helland 2012); first written for the Celestial Atlas |
| `globe/` | Virtual globes, exported as the namespace `globe` (ADR 0049): the geographic tile pyramid (`tiling.js`, the data library's); the figures of the Earth and the Moon with geodetic coordinates, ray hits and Cesium's horizon test and horizon culling point (`figure.js`); terrain tiles decoded and built into meshes with skirts at their edges and under the lines between their quarters, of the visible surface or the sea floor, their triangles ordered by quarter (`terrain.js`); the level of detail (`lod.js`: screen-space error with a falloff toward the horizon, frustum and horizon culling, a walk through tiles not yet loaded so that a close view asks for all its levels at once, parents drawn in the quarters their children cannot yet cover, imagery from ancestors by subregion, eviction); map services cut to the pyramid (`imagery.js`: EOX's Sentinel-2 cloudless 2016, NASA GIBS's Blue Marble, VIIRS, Black Marble, the sea's temperature and its anomaly, rain, coastlines and names, NASA Moon Trek); the camera of a virtual globe and its gestures (`controls.js`); a pool of workers that fetch and build terrain (`workers.js`, `worker.js`); rotations and paths sampled by an engine, interpolated by slerp with the Earth's spin factored out and by cubic Hermite polynomials (`motion.js`); a three.js renderer of a globe (`renderer-three.js`: camera-relative tiles in double precision, a logarithmic depth buffer, the Sun's light, lunar-Lambert shading for the Moon, the sea's glint, night lights, a field and a reference laid over, the provenance of the terrain) and of the sky (`sky-three.js`: the stars of a catalogue, the Sun's disc and the Milky Way in galactic coordinates); a body's map as the pyramid's two root tiles, finer from a WMTS (`imagery.js`, `mapLayer`); the Sun's own lighting with limb darkening; flights of the camera between bodies (`transfer.js`: the distance in its logarithm, raised to see both; [ADR 0050](../development/decisions/0050-instruments-and-their-acquisition.md)) |

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
