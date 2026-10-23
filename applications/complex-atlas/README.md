# Complex Atlas

*Status: draft (2026-10-23). Serves: learning and research. Draws on: the foundations morphe (arithmetic, algebra, topology) and physis (symmetry, gauge, quantum gravity), through its C++ engine.* Decision records: the infrastructure's [ADR 0039](../../development/decisions/0039-complexes-symmetry-and-fields.md); morphe's decision 0004; physis's decision 0005.

The Complex Atlas shows cell complexes, their invariants and the fields on them. A complex is built from a catalogue (spheres and surfaces, presentation complexes, lattices, projective spaces, triangulations, the coloured graphs of tensor models), or written by its user as surface words, as the facets of a triangulation, or as JSON. It is drawn three ways: realised in space, unfolded into a net of polygons, and as a Hasse diagram of its cells. The engine, the foundations' C++, computes its invariants exactly (homology and cohomology with generators, Betti numbers over finite fields, the fundamental group, the surface it is, whether it is a manifold, whether it collapses), samples gauge fields on it by Monte Carlo, and evaluates state sums on triangulated 3-manifolds.

## 1. Concepts

- **A complex is combinatorial.** Vertices are counted; each edge has a tail and a head; each 2-cell has a base vertex and an attaching word, a closed path of edges (lowercase forward, uppercase backward, as group theorists write); a higher cell has its cellular boundary. This is what the fundamental group needs (the words), what homology needs (their abelianisation), and what a gauge field needs (the holonomy of a 2-cell is the product of its links along its word).
- **A complex has no geometry of its own, so the page makes one** (`visualization/js/complex.js`). The *realisation* places the vertices, the midpoints of the edges and the centres of the faces by stress majorisation of their distances in the graph (Gansner, Koren and North 2004), so that loops and parallel edges stay apart; each face is a soap film spanned on its boundary. A complex whose cells are embedded (no loops, no parallel edges, no face meeting a vertex twice) is laid out by its vertices alone, with straight edges; a complex that brings positions for its vertices is drawn where they are. The *net* draws each 2-cell as a regular polygon with unit sides and glues neighbours along shared sides wherever they do not overlap, as a surface is cut open to be classified: the torus becomes its square a b a⁻¹ b⁻¹, a 4 × 4 lattice its periodic fundamental domain. Corners are dotted in the colour of the vertex they are, so that the identifications of corners can be read. The *Hasse diagram* sets the cells in layers by dimension, each joined to its boundary, with incidences of multiplicity other than one written on their links.
- **Invariants are exact.** Homology over ℤ comes from the Smith normal form of the boundary matrices, after the unit pivots of a sparse reduction; each generator is a chain that the page can draw. The fundamental group comes from a spanning tree of the 1-skeleton: a generator for each edge off the tree, a relator for each 2-cell, then Tietze moves. Its generators are loops at the base vertex, drawn on the complex.
- **Fields live on the cells.** A lattice gauge field puts a group element (ℤ/N, U(1), SU(2), SU(N)) on each edge; the Wilson action weighs each face by β (1 − Re tr U_f / N). The engine thermalises and samples it, and returns the last configuration, which colours the edges by the class angle of their link and the faces by their plaquette, the history of the mean plaquette against the exact value of two-dimensional Yang–Mills on a closed surface (Migdal 1975; Witten 1991), and Wilson loops along the generators of the fundamental group.
- **State sums.** On a triangulated closed 3-manifold, the Turaev–Viro sum over colourings of the edges by spins, with quantum 6j symbols at q = e^{iπ/r}, is a topological invariant (Turaev and Viro 1992): the Ponzano–Regge model of three-dimensional quantum gravity, made finite.
- **Your own complexes.** Every view takes what the user writes: words (`abAB cdCD`, `a1b1A1B1a2b2A2B2`), facets (`0 1 2; 0 2 3`), or Complex JSON with positions; sizes of vertices, edges and arrows, the opacity and smoothing of faces, names and the layout's seed are set in *Display*. A complex can be saved as JSON and a realisation as a picture.

## 2. Running

From the repository root, build the engine once (g++ 11 or later, natively or in WSL; about a minute) and serve the page:

```
python data/acquisition/recipe.py applications/complex-atlas/recipe.json
python applications/complex-atlas/serve.py --open
```

The page is at http://127.0.0.1:8850/. Without the engine the views still draw the examples and the user's complexes, and the page offers to build it. With `?fixtures` in the address (http://127.0.0.1:8850/?fixtures) the page answers from `tests/fixtures/`, the engine's own answers for a few complexes (torus, surface of genus 2, projective plane, lattice torus with a U(1) run, boundary of the 4-simplex, a melonic graph, the last two with their Turaev–Viro sums), to see every panel without the engine; `python applications/complex-atlas/tests/make-fixtures.py` writes them again from the built engine.

The address may name the first complex: `#example=torus-7`, `#family=surface&genus=3`, `#words=abAB`, `#facets=0 1 2; 0 2 3`.

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

The server also answers `GET /api/engine` (whether the engine is built from the present sources, and the limits of a request), `POST /api/cancel`, and the recipe at `/setup`. On a public address (`--host 0.0.0.0`) building needs `--token`, and complexes, sweeps and the order of a state sum are held within smaller limits.

## 4. Files

| File | Content |
|---|---|
| `index.html`, `styles.css` | The page: plum nights, Old Rose (#C08081, provisional) for this atlas, Source Serif 4 |
| `src/app.js` | The state of the page and its wiring |
| `src/views/realisation.js`, `net.js`, `hasse.js` | The three views (three.js 0.170.0; SVG) |
| `src/invariants.js`, `fields.js`, `statesums.js`, `catalogue.js` | The panels |
| `src/input.js`, `examples.js`, `format.js`, `selection.js`, `chart.js`, `api.js`, `dom.js` | Reading complexes; the examples; writing groups and words; selections; the history chart; the server's interface |
| `engine.py`, `serve.py`, `recipe.json` | Building and running the engine; the server; the recipe |
| `engine/complex.cpp`, `complexes.hpp`, `fields.hpp` | The engine: its requests; complexes to and from JSON, the catalogue and the analysis (morphe); gauge fields and state sums (physis) |
| `tests/run-tests.mjs`, `tests/fixtures/`, `tests/make-fixtures.py` | The atlas's own checks; the engine's answers for a few complexes, for the page without the engine, and the script that writes them |

## 5. Checks

`node applications/complex-atlas/tests/run-tests.mjs`: how groups, words and presentations are written; what is read as a complex and what is refused; that every example is whole and has the Euler characteristic of its space; that every fixture is whole and obeys the Euler–Poincaré formula; that the page's modules and the import map are in place; and, when the engine is built, every family of the catalogue built and analysed by the engine, each a whole complex obeying Euler–Poincaré. The drawing is checked with the visualization layer (`node visualization/tests/run-tests.mjs`), the mathematics in morphe and physis.

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
3. Fields of other kinds on the same cells: scalar fields on the vertices, cochains from the cohomology, spin labels of spin foams on the faces.
