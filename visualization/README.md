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
| [`scene/`](scene/) | A device-neutral, serialisable description of what is shown. It carries units, time scales and reference frames. |
| [`renderers/`](renderers/) | One renderer per technique, beside the arithmetic they share |
| [`styles/`](styles/) | Colour maps and their validation; later typography and marks |

Code outside the layer reaches them through its public interface, `js/index.js`, which exports them as the namespaces `scene`, `styles` and `renderers`, and the three.js marks as `fieldMarks` and `canvasOf` (ADR 0026's rule of public interfaces, checked by `development/tools/check-dependencies.mjs`).

### 2.1 The Scene Model

A drawing of a field is a *visual expression* $\mathsf V=\mathsf M\circ\rho\circ r$ (Cosmopoiesis's formalism, Definition 3.10; [ADR 0053](../development/decisions/0053-fields-and-their-visual-expression.md)):

- a restriction $r$ of the field to a sub-base it is drawn on;
- a reduction $\rho$ of its fibre to what is shown;
- a mark $\mathsf M$, whose geometry the sub-base decides and whose channels the reduced fibre decides.

The field is given by its descriptor, as an engine serving physis's fields gives it (base, fibre with its unit, representation, provenance, parts and source). An expression is plain data, checked by the formalism's rules before anything is drawn: a vector's colour needs a reduction, a component is never drawn without its frame, a magnitude has no sign to diverge about, a label takes categorical hues.

| File | Content |
|---|---|
| `scene/grid.js` | The samplings of physis's fields (sphere, plane, volume, points): the places of their points and their orders, the same formulas as physis's `sampling.hpp` |
| `scene/samples.js` | A field's answers decoded: values (float32 in base64, or JSON), each point's quality, its reason or note; integral curves with |v| along them |
| `scene/scale.js` | Scales: linear, logarithmic with a floor, diverging about a centre, diverging and logarithmic in the size of the departure; nice ticks and decades |
| `scene/expression.js` | `express(field, restriction, options)` with the formalism's defaults; `validate(expression)` with its rules; `offers(field)`, the expressions a field's base and fibre allow |
| `scene/legend.js` | A legend as data: what is shown, its unit with the prefix of its values, its frame, its scale's ticks, the map's gradient, the provenance and parts, where it could not be read and why |
| `scene/volume.js` | The arithmetic of a field in a volume: the cells shown (all, the shell, three slices), the linear fit and its gradient, the value at the centre, the data of a level surface |

### 2.2 Styles

Every colour map is made in OKLab (Ottosson 2020) and validated in normal vision and under protanopia, deuteranopia and tritanopia (Machado, Oliveira and Fernandes 2009):

- **sequential maps** are monotone in lightness, so that their order reads without hue;
- **diverging maps** have two monotone arms, lightest in the middle, and ends apart;
- **categorical hues** are ordered so that neighbours stay apart.

The simulation is held to colorspacious's (`tests/fixtures/cvd.json`, made by `tests/make_cvd_fixtures.py`).

| File | Content |
|---|---|
| `styles/colour.js` | sRGB and OKLab, ΔE, lightness; ramps through stops, tables for textures, gradients for legends |
| `styles/palettes.js` | The warm family's maps: `plum`, `gold`, `orchid` (sequential); `lavender-rust`, `orchid-gold` (diverging); `warm` (cyclic); `categorical`, the Complex Atlas's order, and `categorical-apart`, the same hues ordered so that neighbours stay apart in every vision (at least 0.158 in normal vision and 0.071 under each deficiency) |
| `styles/cvd.js` | The three deficiencies simulated; a map's lightness profile in every vision; the least separation of a list of colours |

### 2.3 Renderers

| File | Content |
|---|---|
| `renderers/raster.js` | The arithmetic of a renderer, without a device: a value's colour under an expression and a scale (with opacity rising from a clear low end or middle); a drape on a sphere or a plane as an image; curves as segments coloured at their ends |
| `renderers/three/index.js` | The marks in three.js, which is passed in: lines coloured along their length, cells as shaded boxes, a level surface, a drape on a sphere in the field's frame; a drape's image on a canvas |

## 3. Provisional JavaScript

Until the scene model and the renderers exist, the parts of the two atlases' drawing that any application could use live here, in `js/`, behind the public interface `js/index.js` ([ADR 0026](../development/decisions/0026-domains-as-libraries.md)):

| File | Content |
|---|---|
| `contour.js` | Isolines of a gridded field by marching squares, and its highs and lows |
| `flow.js` | Particles carried by a vector field on the sphere, to make a flow visible ([note](theory/particles-on-the-sphere.md)) |
| `units.js` | The units in which quantities are shown, and numbers written for reading |
| `tiles.js` | The tiles of an equirectangular map service that cover a region at a resolution |
| `isosurface.js` | A surface at a level of a volume (a brain, a region), as a triangle mesh, by surface nets |
| `complex.js` | Cell complexes drawn: a CW complex realised in space by stress majorisation, its 2-cells unfolded into a net, its cells as a Hasse diagram; surface words and facets read into Complex JSON; coloured graphs (the Feynman graphs of tensor models and group field theories) drawn in the plane and turned into their 2-complexes; a tetrahedron set in space from its six edge lengths (the quantum tetrahedra of spin foams); a colour for each index of a coloured graph, from the colour maps of `styles/`, which it re-exports under their old names |
| `stars.js` | The colour of a star from its colour index (B−V or Gaia's BP−RP) through its temperature (Ballesteros 2012; Mucciarelli and Bellazzini 2020) and a blackbody's sRGB (Helland 2012); first written for the Celestial Atlas |
| `globe/` | Virtual globes, exported as the namespace `globe` (ADR 0049): the geographic tile pyramid (`tiling.js`, the data library's); the figures of the Earth and the Moon with geodetic coordinates, ray hits and Cesium's horizon test and horizon culling point (`figure.js`); terrain tiles decoded and built into meshes with skirts at their edges and under the lines between their quarters, of the visible surface or the sea floor, their triangles ordered by quarter (`terrain.js`); the level of detail (`lod.js`: screen-space error with a falloff toward the horizon, frustum and horizon culling, a walk through tiles not yet loaded so that a close view asks for all its levels at once, parents drawn in the quarters their children cannot yet cover, imagery from ancestors by subregion, eviction); map services cut to the pyramid (`imagery.js`: EOX's Sentinel-2 cloudless 2016, NASA GIBS's Blue Marble, VIIRS, Black Marble, the sea's temperature and its anomaly, rain, coastlines and names, NASA Moon Trek); the camera of a virtual globe and its gestures (`controls.js`); a pool of workers that fetch and build terrain (`workers.js`, `worker.js`); rotations and paths sampled by an engine, interpolated by slerp with the Earth's spin factored out and by cubic Hermite polynomials (`motion.js`); a three.js renderer of a globe (`renderer-three.js`: camera-relative tiles in double precision, a logarithmic depth buffer, the Sun's light, lunar-Lambert shading for the Moon, the sea's glint, night lights, a field and a reference laid over, the provenance of the terrain) and of the sky (`sky-three.js`: the stars of a catalogue, the Sun's disc and the Milky Way in galactic coordinates); a body's map as the pyramid's two root tiles, finer from a WMTS (`imagery.js`, `mapLayer`); the Sun's own lighting with limb darkening; flights of the camera between bodies (`transfer.js`: the distance in its logarithm, raised to see both; [ADR 0050](../development/decisions/0050-instruments-and-their-acquisition.md)) |

Their checks, and those of the scene model, the styles and the renderers' arithmetic, run with `node visualization/tests/run-tests.mjs`.

## 4. Planned Renderers

| Renderer | Technique |
|---|---|
| `static` | Matplotlib, SVG, PDF |
| `web` | three.js (`renderers/three/`, written); WebGPU, CesiumJS |
| `native` | Vulkan |
| `immersive` | WebXR |
| `sonification` | Sound |

## 5. Style

The visual style is warm and humane, not that of a generic dashboard. The colour maps are defined in `styles/` and validated for colour-vision deficiencies by the checks; each application keeps its own main colour and its provenance hues in its own style sheet.

## 6. Prior Art

Study these before building: OpenSpace, NASA Eyes, CesiumJS, Stellarium and Celestia.
