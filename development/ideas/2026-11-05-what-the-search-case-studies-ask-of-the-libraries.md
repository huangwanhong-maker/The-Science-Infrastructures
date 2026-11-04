---
title: What the Search Case Studies Ask of the Libraries and the Domains
status: seed
date: 2026-11-05
related: [heuresis (its decision 0005), calx, morphe, physis, tyche, eikon, the urban-science, architecture and neuroscience domains, the Urban Science Workbench, research/experiments/2026-11-05-*]
---

# What the Search Case Studies Ask of the Libraries and the Domains

*Fourteen case studies of heuresis (ADR 0073) found what is missing beneath them. What concerns heuresis itself is fixed in its decision 0005. What concerns the other libraries and the domains is listed here: each item is a general method that belongs in the place named, not in a study, by the rule that general methods are written once, in the foundation where they belong.*

## 1. The Foundations

| Where | What | Asked by |
|---|---|---|
| calx | the fast Walsh–Hadamard transform, for the Walsh spectra of landscapes | heuresis `landscapes` |
| calx | low-discrepancy sequences (Sobol, Halton), for space-filling designs and quasi-random sampling of behaviour spaces | heuresis `models`, `diversity` |
| tyche | rank tests and effect sizes (Mann–Whitney, Wilcoxon signed-rank, Vargha–Delaney Â₁₂), and bootstrapped intervals of expected running times | heuresis `analysis`; the benchmark of continuous search |
| tyche | k-means (Lloyd's, with k-means++ seeding), for the centroids of Voronoi cells | CVT-MAP-Elites; the repertoires of an arm |
| tyche | Latin hypercube designs | heuresis `models` |
| morphe | weighted shortest paths (Dijkstra), beside the breadth-first distances of `combinatorics` | the district's walking distances, written provisionally in the study (`walking.hpp`) |
| morphe | the algebraic simplification of terms | programs found by symbolic regression ((a+a)π√(a/GM) for 2π√(a³/GM)) |
| physis | registers of qubits, gates and measurement | a simulated device behind heuresis's Grover sampler, to be held to its exact emulation |
| eikon | a cap on the threads of `grid::parallel_for` | a pipeline evaluated on heuresis's thread pool oversubscribes the machine |

## 2. The Domains

- **architecture `environment`:** shadows and daylight are planned, and the buildings' solar gains ignore their neighbours' shade. The district's study wrote a ground-shadow model for prisms (`sunlight.hpp`), marked as its own.
- **urban-science `landscape`:** a passage can be opened only through a block that is a 2-cell. A block holding a dead-end lane is a hole, and the largest blocks of the Hoddle Grid are holes. The study wrote "carry a blind lane through" as an edit of the network (`lanes.hpp`).
- **urban-science `flow`:** `mixing_of` counts each holding once; mixing weighted by floor area needed physis directly.
- **architecture `form` and `rules`** are planned; floor areas, height bounds and plot ratios were the study's own.
- **neuroscience `snn`:** the response of one neuron to a constant current exists only in the tests, and `analysis` has no f–I curve.
- **robotics:** the simulations' controller gains are now set through their configurations (done, 2026-11-05).

## 3. A Bug in the Urban Science Workbench

`footprintsIn` (`applications/urban-science-workbench/src/district.js`) sends OpenStreetMap's `building:part` outlines to the massing together with the outlines that contain them, so a district's 3-complex holds overlapping cells. In Landscape Dynamics' district, 30 of the 440 footprints are parts; one is `way/140838470`, an office part inside `way/46398815`. Fixing it changes the district's published counts (440 buildings, 565 walls, 98 terraces), so it is left for the workbench's own session and the user's decision.
