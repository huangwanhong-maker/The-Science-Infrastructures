# computer: The Contract, its Tolerance and its Cost

*Status: draft (2026-11-06): the square (`square.hpp`), the cost account with physis's floors (`cost.hpp`), the materials' marks with the list of critical raw materials (`materials.hpp`); `tests/check_computer.cpp`. The invariants of chaos use physis's Lyapunov spectra; the energy distance waits for tyche.* What makes a substrate a computer, and what a result costs. Every run on every backend reports it.

## 1. The Square

The tuple $(\mathcal{S}, \Phi_u, \mathcal{M}, C, E, R, \varepsilon)$ and its check: $d(R(\Phi_u(E(m))), C(m)) \le \varepsilon$ on the task's set $K$, with $E$ and $R$ of the class the task allows (theory note §1.1). For chaotic flows the check is by invariants, not trajectories (theory note §1.4):

- the Lyapunov spectrum and the Kaplan–Yorke dimension (physis `dynamics`);
- the energy distance between invariant measures (tyche `descriptive`);
- power spectra (tyche `time_series`);
- the places and kinds of bifurcations.

## 2. The Account of Energy

The cost of a verified result: fabrication and materials spread over the results of a lifetime, then operation, conversion and verification. Each term carries its source, or is declared an assumption (theory note §1.5). Two floors and one estimate:

- Landauer's bound, from physis `statistical_mechanics`;
- the thermodynamic uncertainty relation's bound on a precision, from the same module;
- conversion by an ADC's figure of merit, the energy per conversion step from Murmann's survey.

## 3. The Marks of Materials

| Mark | Data |
|---|---|
| Supply | the critical raw materials of Regulation (EU) 2024/1252, Annex II, as a table of the data library; its caveat, that it measures supply risk rather than harm, kept with it |
| Toxicity and persistence | each material's classification, with its source |
| Reuse | whether a run consumes the substrate |
| Conditions | temperature, vacuum |
| End of life | biodegradable, recyclable, or neither |

## 4. The Report

Every run returns:

- the square's error against its tolerance;
- the invariants used;
- the energy by term, with the floors;
- the materials' marks;
- the readings' provenance (see [`runtime/`](../runtime/)).
