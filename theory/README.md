# theory: Notes and Derivations

*Status: draft.* Knowledge. See [ADR 0006](../development/decisions/0006-theory-notes-in-subject-tree.md), as amended by [ADR 0008](../development/decisions/0008-independent-libraries.md).

## 1. Purpose

This directory holds expositions, derivations and proofs about the laboratory's subjects. They are written to be understood, and are organised by subject.

## 2. Subjects

- **Fields:** `astronomy/`, `geography/`, `oceanography/`, `meteorology/` and `geodesy/`.
- **Proposed:** `geophysics/`, for the solid Earth's motion and structure (earthquakes, the crust). It has its first note; making geophysics a field awaits the user's agreement (ADR 0024).
- **Bridging topics:** `tides/`, `navigation/`, `air_sea/` and `climate/`.

A subject directory is created together with its first note.

Notes on mathematics, numerics and general physics live in the library that implements them:
- morphe, for mathematics;
- calx, for numerics;
- physis, for general physics.

## 3. Note Format

Notes are written in Markdown with LaTeX mathematics. This format is provisional (CLAUDE.md §8).

Each note begins with front matter:

```yaml
---
title: The Equilibrium Tide
status: draft              # planned | seed | draft | stable | deprecated
implemented-by:            # code that implements results from this note
  - core/domains/tides/
references: [Pugh2014]     # keys in docs/references.bib
---
```

When writing a note:
- state the hypotheses;
- separate definitions from results;
- number equations, so that code can cite them;
- follow the glossary.
