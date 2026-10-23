# 0039. Complexes, Symmetry and Fields in the Foundations; the Complex Atlas

- **Status:** Accepted, provisionally (the plans of morphe's topology and algebraic geometry and of physis's six new modules are first designs; each foundation records its part)
- **Date:** 2026-10-23
- **Decided by:** the user asked for the theoretical-physics, algebraic-geometry and topology libraries to be implemented in parallel, with their structure planned and only what is needed now written: first, the analysis and visualization of 2-complexes and CW complexes, spacetime models, the infrastructure of field theory (quantum field theory and group field theory especially), and tools for research on symmetry and gauge. Claude designed the structure across the foundations, wrote this record, and divided the work among parallel agents, each on one module.

## 1. Context

ADR 0035 mapped the mathematics onto the foundations and wrote the first C++ for the spiking networks. morphe had no code; physis had units, dynamics and the open qubit. The user's next research, quantum gravity and gauge theory in the manner of spin foams and group field theory, lives on 2-complexes: their fundamental groups carry gauge holonomies, their 2-cells curvature; the Feynman graphs of group field theory are coloured graphs whose bicoloured cycles bound the 2-cells of the dual 2-skeleton of a triangulation; and a discrete spacetime is a complex with a geometry. The infrastructure's open decision 10 left physis's modules of theoretical physics to physis.

## 2. Decision

1. **morphe's first slice** ([morphe decision 0004](https://github.com/huangwanhong-maker/morphe/blob/main/development/decisions/0004-the-first-slice-cell-complexes.md)): `arithmetic` (integers without bound, rationals, residues), `algebra` (the Smith normal form, homology with generators, presented and finite groups, homomorphisms into finite groups), `combinatorics` (multigraphs, spanning forests) and `topology` (CW complexes with attaching words, simplicial complexes, a catalogue of constructions, homology, cohomology, fundamental groups, surfaces, manifolds, collapses, dual 2-complexes, coloured graphs). The plans of topology and of algebraic geometry are recorded there; algebraic geometry stays planned.
2. **physis's theoretical physics** ([physis decision 0005](https://github.com/huangwanhong-maker/physis/blob/main/development/decisions/0005-symmetry-spacetime-fields-and-gauge.md)): six new modules, in the order `symmetry` → `geometry` → `spacetime` → `field_theory` → `gauge` → `quantum_gravity`, with `electromagnetism` planned. Lie groups move from `geometry` to `symmetry`. Written now: `symmetry` (U(1), SU(2), SU(N), the Lorentz group and SL(2, ℂ), exact 3j, 6j and 9j symbols, quantum 6j symbols, the Weyl integration formula), `geometry` in part (metrics, curvature, geodesics), `spacetime` (exact solutions, Minkowski space, Regge geometry, causal sets), `field_theory` (scalar fields on complexes, Monte Carlo statistics, Feynman graphs, zero-dimensional field theory), `gauge` (lattice gauge theory on any 2-complex, exact two-dimensional Yang–Mills, flat connections) and `quantum_gravity` (Turaev–Viro state sums, BF amplitudes, coloured tensor models, matrix-model gluings).
3. **calx's forward differentiation** ([calx decision 0005](https://github.com/huangwanhong-maker/calx/blob/main/development/decisions/0005-forward-differentiation.md)): `numerics/differentiation`, Taylor arithmetic of order two, so that a metric written once gives its curvature to rounding.
4. **The order of the foundations is unchanged:** morphe and calx independent; physis on both, now using morphe (its `symmetry` the exact integers, its `gauge` and `quantum_gravity` the complexes). The root `CMakeLists.txt` adds morphe before physis, and the runner of the checks builds morphe as it builds the other foundations.
5. **The Complex Atlas** (`applications/complex-atlas/`, port 8850): a work for learning and research that shows a complex (realised in three dimensions, unfolded into a net of polygons, or as the diagram of its cells), its invariants with generators drawn on it, and the fields on it: gauge fields by Monte Carlo, coloured by their holonomies and compared with exact values, and Turaev–Viro state sums. Its engine is a C++ program on the server, built by its recipe as the Spiking Simulator's is (ADR 0038), composing morphe and physis; the drawing of complexes that any application could use (realisation by stress majorisation, nets, Hasse diagrams, parsing of words and facets) lives in `visualization/js/complex.js`. Its main colour is Old Rose (#C08081), within the warm family, until the user chooses another.

## 3. Alternatives Considered

- **A new library for topology or for theoretical physics.** morphe and physis were designed for exactly this (ADR 0026, 0035); a third would split one subject in two.
- **Group field theory inside `field_theory`.** Its amplitudes are gauge theories on 2-complexes, so it stands above `gauge`, beside the spin foams it generates (physis decision 0005, §3).
- **A domain for quantum gravity.** The theories name no object in nature; they are laws, and belong to a foundation. A domain may come for the observations that test them (cosmology, gravitational waves).
- **Drawing complexes in C++ or in each application.** Drawing is presentation, the visualization layer's; computation is the foundations'.

## 4. Consequences

- The foundations have four new morphe modules, six new physis modules (one only in part) and one calx module, header-only, each with its checks; the whole run of the infrastructure counts them.
- The map of the mathematics (docs/mathematics.md) gives the new modules and moves MSC 22 (Lie groups) to physis `symmetry`; MSC 57 (cell complexes, manifolds), 70, 81 and 83 name their new homes.
- CLAUDE.md §3.1, §5, §7 and §9 and the README name the Complex Atlas and what the foundations hold.
- The applications are five.
