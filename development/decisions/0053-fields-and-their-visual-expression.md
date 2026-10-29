# 0053. Fields and Their Visual Expression: One Model of Fields, Their Laws Apart, and the Scene Model

- **Status:** Accepted, provisionally. Boxes are the first geometric complexes and heat conduction the first balance law; the scene model has one renderer, three.js.
- **Date:** 2026-10-29
- **Decided by:** the user saw that "our model basically a field theory based model", and asked as the next step for "a unified way, model various fields (tensor, scalar, even manifold or complexes), and their formalism. Then the visualizations for them". They added that "the evolution should also contains structured way, like PDE. Also other ways. But it is seperated". They gave a reviewer's design to be read critically, and asked for both ambitions to be implemented elegantly and without mistakes. Claude designed it and wrote it, with this record.

  The foundations' records are physis's [decision 0012](https://github.com/huangwanhong-maker/physis/blob/main/development/decisions/0012-fields-as-sections.md) and calx's [decision 0011](https://github.com/huangwanhong-maker/calx/blob/main/development/decisions/0011-the-theta-method.md). The work's record is Cosmopoiesis's [decision 0007](https://github.com/huangwanhong-maker/Cosmopoiesis-Lab/blob/main/development/decisions/0007-the-worlds-fields-and-their-drawing.md).

## 1. Context

The domains held their models as functions of their own:

- WMM2025 at a geodetic place;
- the standard atmosphere at a height;
- the Moon's crust at a radius and two angles;
- normal gravity at a place.

Each work that read or drew one wrote its own grid, loop, refusal and colouring. Cosmopoiesis traced field lines by a loop of its own, mapped the Moon by another and filled its cube by a third. Its formalism had named what a field and a drawing are (§3.5, §3.6: a phenomenon's kind is base, fibre and process; a drawing is mark ∘ reduction ∘ restriction), and the code had no such type. ADR 0005's scene model, styles and renderers existed only as a plan.

The reviewer's design had these strengths, kept here:

- the base kept apart from the representation;
- capabilities rather than one interface;
- the conserved quantity as the state;
- provenance and versions at the level of the field.

It had these faults, which were not kept:

- a vector-calculus language with a general `solve()`;
- numerics in physis;
- nine kinds of governing model;
- capabilities as virtual interfaces.

## 2. Decision

### 2.1 What Went Where

| Where | What |
|---|---|
| calx `integrators` | the θ-method for $M\dot y=-Ay+c$ with a diagonal mass and a sparse stiffness, solved for the increment by conjugate gradients (calx decision 0011) |
| physis `fields` (new) | a field as a section over its base, with its fibre, representation and provenance; evaluations that refuse with a reason (outside, missing, invalid); capabilities as concepts; formulas with exact derivatives by dual numbers; potential fields read by rows on spheres; uniform and sampled fields; reductions with their invariance; sums, pullbacks, changes of frame and restrictions; samplings; integral curves; versions in time; type erasure for registries (physis decision 0012, its theory note `docs/theory/fields.md`) |
| physis `discretisation` | complexes realised in space (coboundaries in every degree, measures, facets and the orientation identity), boxes, cochains on boxes as fields, and their transfer to a refinement commuting with $d$ |
| physis `continuum` (new) | heat conduction as a balance law: the two-point flux, four kinds of boundary, the θ-step, the energy updated by the fluxes, conserved to rounding whatever the solve's tolerance |
| geodesy `gravity` | `NormalGravityField`: WGS 84's normal gravity on the ITRS axes, its continuation to 11 km below the ellipsoid noted; its divergence is $2\omega^2$, the check's |
| geophysics `geomagnetism` | `WmmField`: WMM2025 on its base, and continued outward (`earth.magnetic.core`), the core's field beyond 850 km, noted |
| meteorology `atmosphere` | `StandardAtmosphereField`: each quantity of the 1976 standard on its line of heights; meteorology now uses physis |
| planetary-science `magnetism` | `LunarCrustalField`: the Moon's crust (T2015_449) or its unresolved part, on low ground read at the sphere above with a note, sharing the expansion it reads |
| visualization `scene/` (new) | the device-neutral scene model: the samplings' places and orders; a field's answers decoded; scales (linear, logarithmic, diverging, diverging with a floor) and their ticks; visual expressions $\mathsf V=\mathsf M\circ\rho\circ r$, their defaults by base and fibre, and the formalism's rules; legends; the arithmetic of a volume |
| visualization `styles/` (new) | colour in OKLab; the warm family's colour maps (plum, gold, orchid; lavender–rust, orchid–gold; warm; categorical); colour-vision deficiencies simulated after Machado et al. (2009), against colorspacious's own simulation |
| visualization `renderers/` (new) | the arithmetic of a renderer without a device (a value's colour, drapes as images, curves as coloured segments), and `three/`, its three.js marks (lines coloured along their length, cells, a level surface, a drape on a sphere) |
| Cosmopoiesis | the world's fields, by name, with the operations `fields`, `field.sample` and `field.curves`; `field_lines`, `field_map` and a field's cube (`probe_volume` with `field`) on them; the cube on physis's volume grid; the page drawing through the scene model, its LAYERS listing the world's fields and drawing a field on a sphere about the Earth, its FIELD cube read by an instrument or a field (its decision 0007) |

### 2.2 One Model of Fields

A field is a section $\varphi\in\Gamma(E)$ over a base, its four axes kept apart:

| Axis | Values |
|---|---|
| base | a region, a surface, a curve, a worldline, time alone, phase space, a complex |
| fibre | a scalar, a vector on named axes, a tensor, a $k$-form, a label, a measure, a connection, each with its unit |
| representation | a formula, a spectral expansion, samples with a reconstruction, a cochain, a composite, a uniform value |
| provenance | physical, archive, model, assimilated, synthetic |

The unification asked for is in the base and the fibre, not in a class per kind of field. A Regge geometry (edge lengths, a metric) and a circulation (a 1-cochain) share a base and differ in fibre and representation. A field has no value outside its base, and says why. What a representation can do is a set of capabilities, concepts read off its type, so that numerical kernels keep their static types. A registry erases them once per batch.

### 2.3 The Evolution Apart

A law that advances a field is separate from the field and from the time integrator:

- **the state** is the conserved quantity on the cells;
- **the law** gives the fluxes and the rate (physis `continuum`);
- **the integrator** is calx's;
- **the result** is published as committed versions of a field, which a reader holds while the next is committed (formalism Proposition 8.6).

Conservation is the law's, not the solver's: the conserved quantity is updated by the fluxes of the solved state, so a loose solve costs accuracy, never conservation. Other ways of evolving a field (a jump, a reaction network, an agent's state) are other laws over the same fields. The first is heat conduction; the next are the ones the works reach.

### 2.4 The Scene Model

A drawing is a visual expression $\mathsf V=\mathsf M\circ\rho\circ r$ (Cosmopoiesis's formalism, Definition 3.10):

- a **restriction** $r$ to a sub-base: a sphere, a plane, a volume, points, or curves through the region;
- a **reduction** $\rho$ of the fibre: none, magnitude, a component, the radial or horizontal part, the components on given axes;
- a **mark** $\mathsf M$: a drape, cells, slices, a shell, a level surface, arrows, curves, glyphs;
- a **colour**: a map and a scale.

The defaults follow the formalism's tables. The rules are checked before anything is drawn:

- a vector's colour needs a reduction;
- integral curves are of a vector field;
- a component is never drawn without its frame;
- a magnitude does not diverge about zero;
- an absolute temperature has no natural zero;
- a label takes categorical hues;
- a diverging map needs a diverging scale;
- a signed value has no logarithmic scale.

Every colour map is made in OKLab and validated in normal vision and under protanopia, deuteranopia and tritanopia:

- sequential maps are monotone in lightness;
- diverging maps have monotone arms and ends apart;
- the categorical hues are ordered so that neighbours stay apart.

The simulation is held to colorspacious (`visualization/tests/fixtures/cvd.json`, made by `make_cvd_fixtures.py`). The renderers' arithmetic is device-neutral and checked in Node, and the three.js renderer takes three.js as an argument, as the globe's does.

## 3. Alternatives Considered

- **A class for each kind of field.** Rejected: the kind is two axes, and a hierarchy would multiply by the representations.
- **One general solver behind a vector-calculus language.** Rejected: a solver depends on the law's structure, which a general interface hides. Each law states its discretisation and uses calx's integrators.
- **Capabilities as virtual interfaces.** Rejected: a virtual call in every evaluation. Concepts keep the kernels' static types, and `AnyField` erases them once per batch.
- **Drawing in each work.** Rejected by ADR 0005, and now replaced: Cosmopoiesis's drawings of fields go through the scene model, the styles and the renderers.
- **Changing the categorical palette's order in place.** Rejected: the Complex Atlas draws its indices with it. The validated order is a palette of its own, `categorical-apart`, the default for new expressions.

## 4. Consequences

- **New checks.**
  - physis: three programs, `check_physis_fields` (49), `check_physis_geometric` (18) and `check_physis_heat` (17), for 31 in all.
  - calx: `check_calx_theta` (8), for 27.
  - The domains' field checks (geodesy 5, geophysics 11, meteorology 5, planetary science 9). Each holds the field to the domain's references (NOAA's test values, pyshtools, the 1976 standard's tables) and to the field equations: no divergence or curl for WMM2025, divergence $2\omega^2$ and no curl for normal gravity.
  - The visualization layer's 60 new checks, for 169: colour-vision simulation against colorspacious, every map in every vision, samplings, answers, scales and ticks, expressions and their rules, legends, volumes, drapes and curves.
  - Cosmopoiesis's 19 new checks of its field operations, for 240, each against a truth outside the operation.
- **Corrections found while writing.**
  - Duplicate entries of the conductance Laplacian were stored apart, and the explicit bound read a part of the diagonal.
  - The meteorology header said $Z(-5\text{ km}')=-5003.9$ m. The standard's geopotential height gives $r_0H/(r_0-H)=-4996.1$ m.
  - A point of the Moon's low ground scaled up to the reference sphere could fall below it by rounding and be refused. It is now read on the sphere itself, at the point's latitude and longitude.
  - physis called a composite's provenance "as certain as its least certain part". The order is of derivation, not of accuracy, and is renamed `most_derived`.
- **The categorical palette.** The ten warm hues of the Complex Atlas, in their order, are adjacent within 0.035 in OKLab under deuteranopia. In `categorical-apart`, the best of the 10! orders, neighbours are at least 0.158 apart in normal vision and 0.071 under each deficiency. Some pairs of the ten stay within 0.003 under deuteranopia, so identity is carried by labels as well as by hue.
- **Field lines end at the base.** The Earth's field lines stop where they first leave WMM2025's base, 1 km below the ellipsoid, or the shell from its reference sphere of 6371.2 km. Near the poles that is the sphere, as before. At low latitudes the base lies above the sphere, by up to some 6 km at the equator, and the lines end there. Each end now says why it ends.
