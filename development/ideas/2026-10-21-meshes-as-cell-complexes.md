---
title: Meshes as Cell Complexes
status: seed
date: 2026-10-21
related: [lib/morphe/topology, lib/physis/discretisation, lib/physis/continuum]
---

# Meshes as Cell Complexes

## 1. Idea

Use one topology library both for pure research (homology, CW complexes) and for simulation.

- A simulation mesh is a cell complex, and the fields on it are cochains.
- Discrete exterior calculus (DEC) and finite-element exterior calculus (FEEC) build PDE schemes on such complexes. In these schemes the discrete exterior derivative satisfies d∘d = 0 exactly.
- As a result, conservation laws, and some topological invariants, hold exactly in the discrete system.

## 2. Relevance

- **Ocean models on the sphere.** Shallow-water and ocean models conserve mass and vorticity on cell complexes of the sphere, such as icosahedral and cubed-sphere grids.
- **Islands and circulation.** Coastlines make the ocean a surface that is not simply connected.
  - Its first Betti number counts the islands, and equals the dimension of the space of harmonic flows in the Hodge decomposition.
  - These harmonic flows carry the circulation around each island, which is the subject of the island rule.

## 3. Questions

- **First cell structures.** Which should we support first: simplicial (triangles) or polygonal (Voronoi) complexes?
- **One data structure.** Can mesh refinement and algebraic operations on complexes share a single data structure?
