# Complex Atlas

*Status: draft (2026-10-23). Serves: learning and research. Draws on: the foundations morphe (arithmetic, algebra, topology) and physis (symmetry, field theory, gauge, quantum gravity, spacetime), through its C++ engine.* Decision records: the infrastructure's [ADR 0039](../../development/decisions/0039-complexes-symmetry-and-fields.md) and ADR 0040; morphe's decision 0004; physis's decisions 0005 and 0006.

The Complex Atlas shows cell complexes, their invariants and the fields on them. A complex is built from a catalogue (spheres and surfaces, presentation complexes, lattices, projective spaces, triangulations, the coloured graphs of tensor models), or written by its user as surface words, as the facets of a triangulation, or as JSON. It is drawn three ways: realised in space, unfolded into a net of polygons, and as a Hasse diagram of its cells. The engine, the foundations' C++, computes its invariants exactly (homology and cohomology with generators, Betti numbers over finite fields, the fundamental group, the surface it is, whether it is a manifold, whether it collapses), samples gauge fields on it by Monte Carlo, evaluates state sums on triangulated 3-manifolds, lays spin foams on them, and measures the dimensions a complex has as a space.

Its second mode, **Generation**, walks the generation of spacetime by a field theory (physis decision 0006) in five stations: an action S[φ] with its interactions drawn as bubbles; its Feynman graphs at an order, each a coloured graph weighted by N to the number of its faces, with the series of Z and log Z; a graph's discrete spacetime (its 2-complex, its pseudomanifold, its spin foam), opened in the complex views; the large-N or condensed phase (the eigenvalue density and planar free energy of matrices, the melonic resummation of tensors, the gap equation of vectors, the condensate of a group field theory and its bounce); and the continuum, a random geometry of that phase with its spectral and Hausdorff dimensions. The same stations serve quantum field theories of every rank (scalar, vector, matrix, tensor) and group field theories (Boulatov's, Ooguri's, a melonic one).

## 1. Concepts

- **A complex is combinatorial.** Vertices are counted; each edge has a tail and a head; each 2-cell has a base vertex and an attaching word, a closed path of edges (lowercase forward, uppercase backward, as group theorists write); a higher cell has its cellular boundary. This is what the fundamental group needs (the words), what homology needs (their abelianisation), and what a gauge field needs (the holonomy of a 2-cell is the product of its links along its word).
- **A complex has no geometry of its own, so the page makes one** (`visualization/js/complex.js`). The *realisation* places the vertices, the midpoints of the edges and the centres of the faces by stress majorisation of their distances in the graph (Gansner, Koren and North 2004), so that loops and parallel edges stay apart; each face is a soap film spanned on its boundary. A complex whose cells are embedded (no loops, no parallel edges, no face meeting a vertex twice) is laid out by its vertices alone, with straight edges; a complex that brings positions for its vertices is drawn where they are. The *net* draws each 2-cell as a regular polygon with unit sides and glues neighbours along shared sides wherever they do not overlap, as a surface is cut open to be classified: the torus becomes its square a b a⁻¹ b⁻¹, a 4 × 4 lattice its periodic fundamental domain. Corners are dotted in the colour of the vertex they are, so that the identifications of corners can be read. The *Hasse diagram* sets the cells in layers by dimension, each joined to its boundary, with incidences of multiplicity other than one written on their links.
- **Invariants are exact.** Homology over ℤ comes from the Smith normal form of the boundary matrices, after the unit pivots of a sparse reduction; each generator is a chain that the page can draw. The fundamental group comes from a spanning tree of the 1-skeleton: a generator for each edge off the tree, a relator for each 2-cell, then Tietze moves. Its generators are loops at the base vertex, drawn on the complex.
- **Fields live on the cells.** A lattice gauge field puts a group element (ℤ/N, U(1), SU(2), SU(N)) on each edge; the Wilson action weighs each face by β (1 − Re tr U_f / N). The engine thermalises and samples it, and returns the last configuration, which colours the edges by the class angle of their link and the faces by their plaquette, the history of the mean plaquette against the exact value of two-dimensional Yang–Mills on a closed surface (Migdal 1975; Witten 1991), and Wilson loops along the generators of the fundamental group.
- **State sums.** On a triangulated closed 3-manifold, the Turaev–Viro sum over colourings of the edges by spins, with quantum 6j symbols at q = e^{iπ/r}, is a topological invariant (Turaev and Viro 1992): the Ponzano–Regge model of three-dimensional quantum gravity, made finite.
- **Spin foams.** On a closed triangulated 3-manifold, or a 4-coloured graph (a Feynman graph of Boulatov's group field theory), a spin foam sits on the dual 2-complex: a spin j on each 2-cell (an edge of the triangulation, of length j + ½), and at each vertex (a tetrahedron) the 6j symbol of its six spins, quantum at q = exp(iπ/r). The Spin foam panel samples spins (Metropolis with weight |A|, its sign kept) or sets them uniform; the realisation colours each face by its spin and sizes each vertex by its amplitude (gold positive, lilac negative); a vertex opens an inset with its boundary spin network (the tetrahedral graph, its nodes the four triangles, its links the six edges) and its quantum tetrahedron in space, with its volume, its dihedral angles, and the 6j symbol beside the Ponzano–Regge asymptotic cos(Σ (j + ½) θ + π/4)/√(12πV), the bridge to the Regge action of discrete gravity.
- **Emergent dimensions.** For any complex the Invariants panel measures the spectral dimension (from the return probability of a lazy random walk) and the Hausdorff dimension (from the growth of balls), on its 1-skeleton or its dual graph: the lattice's dimension on a lattice; 4 and 2 on a random planar quadrangulation (the Brownian map of two-dimensional quantum gravity); 2 and 4/3 on branched polymers, the melonic phase of tensor models.
- **Your own complexes.** Every view takes what the user writes: words (`abAB cdCD`, `a1b1A1B1a2b2A2B2`), facets (`0 1 2; 0 2 3`), or Complex JSON with positions; sizes of vertices, edges and arrows, the opacity and smoothing of faces, names and the layout's seed are set in *Display*. A complex can be saved as JSON and a realisation as a picture.

## 2. Running

From the repository root, build the engine once (g++ 11 or later, natively or in WSL; about a minute) and serve the page:

```
python data/acquisition/recipe.py applications/complex-atlas/recipe.json
python applications/complex-atlas/serve.py --open
```

The page is at http://127.0.0.1:8850/. Without the engine the views still draw the examples and the user's complexes, and the page offers to build it. With `?fixtures` in the address (http://127.0.0.1:8850/?fixtures) the page answers from `tests/fixtures/`, the engine's own answers for a few complexes (torus, surface of genus 2, projective plane, lattice torus with a U(1) run, boundary of the 4-simplex, a melonic graph, the last two with their Turaev–Viro sums), to see every panel without the engine; `python applications/complex-atlas/tests/make-fixtures.py` writes them again from the built engine.

The address may name the first complex: `#example=torus-7`, `#family=surface&genus=3`, `#words=abAB`, `#facets=0 1 2; 0 2 3`; `&foam=1` lays a spin foam on it (`#family=melon&D=3&insertions=2&seed=1&foam=1`). `#mode=generation&theory=matrix&order=2` opens the Generation mode on a theory, expanded to an order, and `&stage=large-n` (or `action`, `graphs`, `spacetime`, `continuum`) at a station. Its passages open it at good points: the BIPZ matrix model at order 2, coloured graphs and their degree, Boulatov's GFT and its spin foams, a condensate with one spin and its bounce, a random planar quadrangulation. With `?fixtures` the Generation mode answers from fixtures too (`expand-<theory>-<order>`, `large-n-<theory>`, `condensate`, `random-geometry-<kind>`, `dimensions-<name>`, `spin-foam-<name>`).

- **In the realisation:** drag to turn, scroll to approach, right-drag to move; the cell under the pointer is named below.
- **In the net and the Hasse diagram:** drag to move, scroll to zoom.
- **Generators:** click a homology generator, a generator of the fundamental group or a Wilson loop to draw it on the complex; click it again to clear it.

## 3. The Engine and Its Protocol

The engine is `engine/complex.cpp`, header-only C++20 on morphe, calx and physis, with the API's JSON (`api/cpp/science/json.hpp`). `engine.py` compiles it into `data/cache/engines/complex/` with the hash of its sources (any change makes it stale) and writes the catalogue there; `serve.py` runs it once per request, a JSON request on its standard input and a JSON answer on its standard output, an error as `{"error": "..."}`.

**A complex:**

```
{"name": "torus", "title": "Torus",
 "vertices": 1,
 "edges": [[0, 0], [0, 0]],                          [tail, head]
 "faces": [{"base": 0, "boundary": [1, 2, -1, -2]}], signed 1-based edges: k is edge k-1 forward, -k backward
 "cells": [[{"boundary": [[face, coefficient], ...]}, ...], ...],   cells[k-3]: the k-cells
 "positions": [[x, y, z], ...],                      optional, one per vertex
 "labels": {"edges": ["a", "b"], "faces": [...]},    optional
 "simplices": [[0, 1, 2], ...],                      optional: the facets, when it is simplicial
 "coloured_graph": {"colours": 4, "edges": [...]}}   optional: per colour, the white neighbour of each black vertex
```

**Requests** (the server's routes in brackets):

| Request | Answer |
|---|---|
| `{"op": "catalogue"}` (GET `/api/catalogue`, cached at build) | `{"families": [{"id", "title", "group", "description", "params": [{"name", "type": "int" \| "number" \| "text", "min", "max", "default"}]}]}` |
| `{"op": "build", "family", "params"}` (POST `/api/build`) | `{"complex"}` |
| `{"op": "analyse", "complex"}` (POST `/api/analyse`) | `{"analysis": {"dimension", "counts", "euler", "components", "homology": [{"dimension", "rank", "torsion", "generators": [{"order": 0 if free, else the torsion order, "chain": [[cell, coefficient]]}]}], "cohomology": [{"dimension", "rank", "torsion"}], "betti": {"2": [...], "3": [...]}, "fundamental_group": {"generators", "relators", "simplified": {"generators", "relators"}, "loops", "abelianisation": {"rank", "torsion"}} or null, "surface": {"closed", "orientable", "genus", "crosscaps", "boundary_components", "name", "normal_form"} or null, "manifold": {"pseudomanifold", "manifold", "dimension"} or null, "collapse": {"collapsible", "remaining", "steps"}, "coloured_graph": {"colours", "vertices", "faces", "bubbles", "jackets", "degree", "melonic"} or null}}` |
| `{"op": "gauge", "complex", "group": {"kind": "Z", "n"} \| {"kind": "U1"} \| {"kind": "SU2"} \| {"kind": "SU", "n"}, "beta", "sweeps", "thermalise", "seed", "update": "heatbath" \| "metropolis", "loops": [words]}` (POST `/api/gauge`) | `{"group", "link_angle": [per edge, in [0, 2π)], "links": [per edge, raw], "plaquettes": [Re tr U_f / N per face], "plaquette_angle", "history": {"plaquette"}, "mean": {"plaquette", "error", "tau"}, "exact": {"plaquette"} or null, "loops": [{"loop", "mean", "error"}], "charge" or null, "acceptance"}` |
| `{"op": "turaev_viro", "complex", "r"}` (POST `/api/turaev-viro`) | `{"value", "colourings", "r"}` |
| `{"op": "theories"}` (GET `/api/theories`) | `{"theories": [{"id", "title", "family": "QFT" \| "GFT", "rank", "description", "action", "params", "bubbles": [{"name", "coloured_graph": {"colours", "p", "sigma"}}]}]}`; ids `scalar`, `vector`, `matrix`, `tensor-melonic`, `tensor-wheel`, `tensor-necklace`, `coloured-simplicial`, `gft-boulatov`, `gft-ooguri`, `gft-melonic`; the coloured models have no bubbles (every colour is a propagator, each vertex a simplex) |
| `{"op": "expand", "theory", "params", "order"}` (POST `/api/expand`) | `{"graphs": [{"id", "coloured_graph", "counts", "connected", "n_power", "degree" or "genus" or null, "automorphisms", "weight": "p/q", "amplitude" or null}], "total", "truncated", "series": {"Z": [{"order", "terms": [{"power", "coeff": "p/q"}]}], "F": [...]}, "leading": {"power", "coeff"}}`; `n_power` counts the faces of colours (0, i), colour 0 the propagators, or all faces in a model without bubbles, and `n_power_scaled` subtracts the couplings' large-N scaling; `total` counts the Wick contractions, `graphs` their classes, the leading first; the series run over orders 1..n as exact coefficients of λⁿ (of (λλ̄)ⁿ, `variable`, in the coloured models, whose Z is the exponential of F) |
| `{"op": "graph_complex", "coloured_graph", "view": "2-complex" \| "pseudomanifold"}` (POST `/api/graph-complex`) | `{"complex"}`, carrying its `coloured_graph` |
| `{"op": "large_n", "theory", "params"}` (POST `/api/large-n`) | matrices: `{"kind": "matrix", "g", "a2", "density": {"x", "rho"}, "free_energy": {"g", "E0"}, "critical": {"g_c", "gamma", "gamma_estimate"}, "series"}` (the Hermitian quartic model V = M²/2 + gM⁴, g_c = −1/48); tensors: `{"kind": "melonic", "m", "lambda", "G", "critical": {"lambda_c", "gamma", "gamma_estimate"}, "series", "equation"}` (G = 1 + x G^m); vectors: `{"kind": "vector", "coupling", "gap"}`; group field theories: the melonic answer with `"kind": "condensate"`. `gamma` is the exact exponent, `gamma_estimate` the one read from 400 exact coefficients; `series` are exact rationals as text |
| `{"op": "condensate", "spins": [{"two_j", "Q", "E"}], "G", "phi": [min, max], "points"}` (POST `/api/condensate`) | `{"phi", "V", "H" (V′/3V), "rho" (matter density), "classical": {"V"}, "bounce": {"phi", "V", "rho", "rho_c"} or null, "friedmann": {"rho", "H2", "classical_H2"}}` |
| `{"op": "random_geometry", "kind": "quadrangulation" \| "melonic", "size", "D", "seed"}` (POST `/api/random-geometry`) | `{"complex"}` |
| `{"op": "dimensions", "complex", "graph": "skeleton" \| "dual", "steps", "starts", "seed"}` (POST `/api/dimensions`) | `{"spectral": {"t", "P", "ds"}, "hausdorff": {"r", "ball", "dh"}}`, each curve only where its estimator is defined (its centred differences need room at the ends); `"graph": "coloured"` measures a coloured graph itself |
| `{"op": "spin_foam", "complex", "r", "mode": "sample" \| "uniform" \| "given", "two_j", "spins", "sweeps", "seed"}` (POST `/api/spin-foam`) | `{"foam" (its 2-complex), "face_spins": [two_j per face], "face_edge", "amplitude", "sign", "vertices": [{"faces": [6, in the order 01, 02, 03, 12, 13, 23], "spins", "amplitude", "asymptotic", "tetrahedron": {"lengths", "volume", "euclidean", "dihedral"}}], "statistics": {"mean_two_j", "sign_mean", "acceptance", "samples"}, "tv" or null, "regge": {"deficits", "action"} or null}` |

The server also answers `GET /api/engine` (whether the engine is built from the present sources, and the limits of a request), `POST /api/cancel`, and the recipe at `/setup`. On a public address (`--host 0.0.0.0`) building needs `--token`, and complexes, sweeps, the order of a state sum or an expansion and the size of a random geometry are held within smaller limits.

## 4. Files

| File | Content |
|---|---|
| `index.html`, `styles.css` | The page: plum nights, Old Rose (#C08081, provisional) for this atlas, Source Serif 4 |
| `src/app.js` | The state of the page and its wiring |
| `src/views/realisation.js`, `net.js`, `hasse.js` | The three views (three.js 0.170.0; SVG) |
| `src/invariants.js`, `fields.js`, `statesums.js`, `spinfoam.js`, `dimensions.js`, `catalogue.js` | The panels: invariants (with the emergent dimensions), gauge fields, state sums, spin foams (with the inset of a vertex), the catalogue |
| `src/generation.js` | The Generation mode: theories, the five stations, their charts and drawings |
| `src/graphs.js`, `plot.js`, `views/tetrahedron.js` | Drawings of bubbles, Feynman graphs and spin networks (SVG); charts of curves and histograms; the quantum tetrahedron (three.js) |
| `src/input.js`, `examples.js`, `format.js`, `selection.js`, `chart.js`, `api.js`, `dom.js` | Reading complexes; the examples; writing groups and words; selections; the history chart; the server's interface |
| `engine.py`, `serve.py`, `recipe.json` | Building and running the engine; the server; the recipe |
| `engine/complex.cpp`, `complexes.hpp`, `fields.hpp` | The engine: its requests; complexes to and from JSON, the catalogue and the analysis (morphe); gauge fields and state sums (physis) |
| `tests/run-tests.mjs`, `tests/fixtures/`, `tests/make-fixtures.py` | The atlas's own checks; the engine's answers for a few complexes, for the page without the engine, and the script that writes them |

## 5. Checks

`node applications/complex-atlas/tests/run-tests.mjs`: how groups, words and presentations are written; what is read as a complex and what is refused; that every example is whole and has the Euler characteristic of its space; that every fixture is whole and obeys the Euler–Poincaré formula; the fixtures of the generation of spacetime (the nine theories with their bubbles; every Feynman graph weighted by N to its own faces, drawn, its 2-complex whole; the matrix model at first order, 1 − 2N³λ, and the vector model at second, N(N+1)(N+2)(N+3)/2; the eigenvalue density of mass 1; the condensate classical at low density and still at its bounce; random geometries whole, the quadrangulation a sphere of squares; spin foams whose vertices carry their faces' spins and edges of length j + ½); that the page's modules and the import map are in place; and, when the engine is built, every family of the catalogue built and analysed by the engine, each a whole complex obeying Euler–Poincaré, and, where the engine has them, its theories, its expansion of the matrix model and its spin foam on ∂Δ⁴. The drawing is checked with the visualization layer (`node visualization/tests/run-tests.mjs`), the mathematics in morphe and physis.

## 6. References

- Gansner, E. R., Koren, Y. and North, S. (2004), Graph drawing by stress majorization, *Graph Drawing 2004*, LNCS 3383, 239–250.
- Brandes, U. and Pich, C. (2006), Eigensolver methods for progressive multidimensional scaling of large data, *Graph Drawing 2006*, LNCS 4372, 42–53.
- Massey, W. S. (1967), *Algebraic Topology: An Introduction*, Harcourt, Brace & World, ch. 1 (surfaces cut into polygons).
- Hatcher, A. (2002), *Algebraic Topology*, Cambridge University Press, ch. 1–2 (CW complexes, cellular homology, the fundamental group of a 2-complex).
- Migdal, A. A. (1975), Recursion equations in gauge field theories, *Sov. Phys. JETP* 42, 413–418; Witten, E. (1991), On quantum gauge theories in two dimensions, *Commun. Math. Phys.* 141, 153–209.
- Turaev, V. G. and Viro, O. Y. (1992), State sum invariants of 3-manifolds and quantum 6j-symbols, *Topology* 31, 865–902.

## 7. Open Questions

1. Views of spacetimes beside complexes: causal sets as Hasse diagrams of their order, geodesics in exact solutions, Regge geometries coloured by deficit angle.
2. A realisation that respects the topology where it can: an embedding of a surface by its genus, or a harmonic map onto a model surface, in place of stress majorisation.
3. Fields of other kinds on the same cells: scalar fields on the vertices, cochains from the cohomology.
4. Spin foams in four dimensions: the 15j vertex of BF theory and the EPRL vertex on SL(2, ℂ), their 4-simplices drawn as the tetrahedra are now.
