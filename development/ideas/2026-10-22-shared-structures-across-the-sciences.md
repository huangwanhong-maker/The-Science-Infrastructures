---
title: Shared Structures Across the Sciences
status: seed
date: 2026-10-22
related: [morphe, calx, physis]
---

# Shared Structures Across the Sciences

## 1. Idea

If the sciences are to be one body of knowledge, the shared libraries should be organised by the structures the sciences have in common, not by the sciences themselves. Five structures recur.

- **Symmetry.**
  - Groups act on crystals and molecules (point and space groups), on the sphere (SO(3): spherical harmonics in geodesy, astronomy and quantum mechanics), on spacetime, and on sequences.
  - Their representations give spectra, selection rules and normal modes.
  - Home: morphe (groups and representations), physis (their action on fields).
- **Conservation and transport.**
  - The same balance laws, ∂ₜu + ∇·F(u) = S, move heat in the ocean, solutes in a cell, carbon through the biosphere and individuals across a landscape (reaction–diffusion).
  - The numerics (finite volumes, spectral methods) and their invariants are shared.
  - Home: physis (continuum), calx (discretisation and solvers).
- **Networks and complexes.**
  - Food webs, metabolic and gene-regulatory networks, neural circuits, trade and transport are graphs with dynamics on them.
  - Meshes and simplicial complexes generalise them; homology measures their holes.
  - Home: morphe (graphs, complexes, topology).
- **Sequences and grammars.**
  - Genomes, proteins, languages and programs are strings with structure.
  - Formal language theory classifies them: regular, context-free and context-sensitive grammars, with symmetry groups acting on them.
  - Home: morphe (formal structures).
- **Inference under uncertainty.**
  - Every science estimates hidden states from noisy, incomplete observations: data assimilation in weather, phylogenetics in biology, econometrics, Bayesian cosmology.
  - Filtering, sampling and variational methods are shared.
  - Home: calx (numerics), with a shared inference layer still to be placed.

## 2. Consequences to Explore

- **A cross-index.** For each structure, the places it appears in each domain, kept in this repository as the map of connections. Where the bridges join domains through shared phenomena, the cross-index would join them through shared structure.
- **Tests that cross fields.** One spectral solver, verified once on a standard problem (Williamson's shallow-water tests, say), then used by an ocean model and a reaction–diffusion model of a population.
- **Teaching by structure.** A learning path that meets the rotation group in geodesy, in the hydrogen atom and in molecular symmetry, in that order.

## 3. The Place of Inference

The question this idea raised, where the shared inference layer (estimation, filtering, assimilation, uncertainty, learning) belongs, is settled by [ADR 0026](../decisions/0026-domains-as-libraries.md): it becomes a fourth foundation, `planned` until its first code, and named then.
