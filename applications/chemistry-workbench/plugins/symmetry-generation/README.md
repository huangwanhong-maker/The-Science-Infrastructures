# Symmetry and Generation

*Status: draft (2026-10-29). Activity: research.* The first plug-in of the Chemistry Workbench ([ADR 0051](../../../../development/decisions/0051-the-chemistry-workbench.md)).

## 1. Purpose

The user's manuscript *An Inquiry into Chemical Structure Generation through Symmetry Breaking: Two Formalisms Drawing on Group-Equivariant Grammar and Group Field Theory* (W. Huang, v1.4, 2026-10-28) proposes two symmetry-first formalisms:
- a group-equivariant chemical graph grammar with its stochastic dynamics;
- a group-field-inspired generator of 2-complexes read chemically through a readout.

This plug-in executes both on the foundations and verifies the manuscript's formalism claim by claim. The manuscript's source is kept beside the workbench as background (`../../Chemistry_Symmetry_Manuscript_v1.4_…zip`); it is not part of the code.

## 2. What It Found

[`theory/formal-verification.md`](theory/formal-verification.md) proves what holds and lists what the revision should change (§7). The engine's `checks` operation runs the computations behind each claim, and the workbench's tests require every one to pass.

- **Four errors.**
  - Definition 4.2 cannot change a kept atom's charge; it needs relabelling.
  - For abelian G, the symmetry §5.3 breaks is the gauge redundancy §5.2 sets aside.
  - §5.3's η lies along a flat direction of the interaction; (−1)^b gives S = 4μη² + 64λη⁴.
  - The two-vertex complex comes from the Wick expansion, not from a field configuration.
- **Seven amendments.** They concern:
  - η's incompleteness;
  - the necessary regulator;
  - the expansion over stranded graphs;
  - which of the 24 gluings is meant;
  - τ relative to an order, X₊ = (R);
  - the measure;
  - Conjecture 7.1's known instances; notation and two sentences.
- **Three stronger statements.**
  - Stabilizers are constant along flows.
  - Frank's excess keeps its sign and tends to ±1.
  - The interaction's symmetry is G × G.

## 3. Views

| View | Shows |
|---|---|
| Grammar | rules L ⇒ R with what they keep, and derivation spaces with counted derivations. Examples: ethanol, glycine backbones (orbits of S₃), the stereocentre (closure and equivariance), the acid–base edit (relabelling, atom ledger, matches against occurrences) |
| Stereocentre | X₊ and X₋ in 3D as mirror images, τ and its CIP descriptor, congruences; (Z)- and (E)-1,2-dichloroethene, one graph, one dipole allowed |
| Point groups | 15 idealized molecules with their symmetry elements, point groups, character tables, and the fixed polar and axial vectors |
| Selection | Frank's autocatalysis: mass action and the excess m with its equation; the Markov chain's exactly symmetric law and asymmetric paths; the pitchfork |
| Complexes | the six-cycle, its filled complex and the two-vertex complex, with boundary matrices and homology |
| Field | the tetrahedral interaction over Z₂ on a field one sets, its G × G symmetry, the characters and order parameters with and without the gauge projection, the action along (−1)^a and (−1)^b, and the Wick expansion to order 3 with every class's 2-complex, homology and amplitudes |

## 4. Engine

`engine/main.cpp` answers JSON lines. Its operations:

| Operation | What it answers |
|---|---|
| `catalogue` | the molecules with their point groups |
| `grammar {example}` | a grammar's derivation space |
| `complexes` | the three complexes of Case Study II |
| `field {group, psi}` | the interaction over a finite group |
| `wick {n, group}` | the Feynman classes at order n |
| `frank {k, ki, k0, l0, r0, t, precursor}` | Frank's network, deterministic and stochastic |
| `pitchfork {b}` | the cubic normal form |
| `landau {mu, lambda}` | the action along the characters |
| `checks` | the 55 checks of the manuscript's claims |

Its headers:

| Header | What it holds |
|---|---|
| `molecules.hpp` | the provisional chemistry (elements, geometries, τ and CIP) and point groups as JSON |
| `grammars.hpp` | Formalism I on morphe |
| `fields.hpp` | Formalism II on physis and morphe |
| `selection.hpp` | kinetics, the normal form, Landau |
| `checks.hpp` | the checks of the claims |

## 5. References

The manuscript's own; and those of the note (Boulatov 1992; Habel and Plump 2002; Frank 1953; Gillespie 1977; Kondepudi et al. 1990; Soai et al. 1995; Nye 1957; Cotton 1990; IUPAC 2013).
