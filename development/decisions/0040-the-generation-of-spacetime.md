# 0040. The Generation of Spacetime by Field Theories; Spin Foams in the Complex Atlas

- **Status:** Accepted, provisionally
- **Date:** 2026-10-23
- **Decided by:** the user set out a formalism of emergent spacetime (an action S[φ], its perturbative expansion into 2-complexes, spin foams and discrete spacetimes, a condensed phase, a continuous spacetime), asked whether the libraries were designed to investigate it and to show it, spin foams above all, and asked that group field theory and quantum field theory alike be supported. Claude assessed the libraries, designed the framework ([physis decision 0006](https://github.com/huangwanhong-maker/physis/blob/main/development/decisions/0006-the-generation-of-spacetime-by-field-theories.md)) and this record, and divided the work among parallel agents.

## 1. Context

ADR 0039 gave the foundations complexes, symmetry, fields, gauge theory and state sums. The assessment found the layering sound (structures in morphe, physics in physis, presentation in the application) and the passage between the stages missing: no action as an object from which Feynman rules follow, no spin foam as an object with spins and amplitudes, nothing of the condensed phase, and no measure of the continuum that results.

## 2. Decision

1. **One framework for every theory** whose expansion generates discrete spacetimes, in five stages (physis decision 0006, §2.1): the action; the perturbative expansion into coloured Feynman graphs with exact series; the discrete spacetime (2-complex, pseudomanifold, spin foam); the large-N or condensed phase; the continuum (critical exponents, effective equations, emergent dimensions). Theories by the rank of their field: scalar, vector, matrix and tensor field theories, and group field theories, whose amplitudes are spin foams.
2. **In the foundations:** bubbles and their gluing in morphe `topology`; invariant theories and their large-N solutions in physis `field_theory`; group field theories, spin foams, random geometries and condensates in physis `quantum_gravity`; emergent dimensions in physis `spacetime`.
3. **In the Complex Atlas:** a Generation mode that walks the five stages for a chosen theory, from the action and its bubbles through the gallery of Feynman graphs and their series to the 2-complex of a chosen graph, the large-N density or the condensate's bounce, and a random geometry of the continuum with its dimensions; and a Spin foam panel that colours a foam's faces by spin and opens each vertex onto its boundary spin network and its quantum tetrahedron.

## 3. Alternatives Considered

- **Group field theory alone.** Matrix and tensor models are where the continuum limit is known exactly; with them the machinery is checked before group field theory relies on it, and the user asked for quantum field theory as well.
- **A separate application for spin foams.** The foams are 2-complexes; the atlas already draws, analyses and runs fields on them, and the passage from a Feynman graph to its foam is one click within one application.

## 4. Consequences

- The foundations' new files carry their checks against exact results; the whole run counts them.
- The atlas's engine gains operations for theories, expansions, large-N solutions, condensates, random geometries, dimensions and spin foams; its README gives the protocol.
