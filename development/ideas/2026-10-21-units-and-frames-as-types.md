---
title: Units and Frames as Types in the IR
status: seed
date: 2026-10-21
related: [lib/calx/compute/ir, lib/morphe/algebra, lib/physis/units, core/metrology]
---

# Units and Frames as Types in the IR

## 1. Idea

Make dimensional correctness and frame correctness properties that the compiler checks, rather than conventions that people follow.

## 2. Sketch

- **Dimensions.** Over the seven SI base units, physical dimensions form the free abelian group ℤ⁷, or ℚ⁷ if fractional powers are allowed.
  - A quantity's type is an element of this group.
  - Multiplying two quantities adds their exponents.
  - Two quantities can be added only if their dimensions are equal.
  - Buckingham's π theorem becomes a computation: the kernel of the dimension matrix over ℚ.
- **Frames.** Reference frames form a groupoid. Its objects are frames, and its morphisms are transformations between them, which may depend on time. A coordinate whose type names its frame can be combined only with coordinates in the same frame, or after an explicit morphism.
- **Keeping calx domain-agnostic.** The IR in calx would provide only generic mechanisms:
  - type annotations valued in an abelian group supplied by the user;
  - tags whose combination is governed by a groupoid.

  physis would then instantiate these with SI dimensions, and the laboratory's `core/metrology` with astronomical frames. The groups and groupoids themselves are studied and implemented in morphe.

## 3. Questions

- **Type inference.** Unification in abelian groups is decidable (Kennedy's dimension types). How should inference work for generic kernels?
- **Time scales.** These could form a further groupoid (UTC, TAI, TT, TDB). Some of its morphisms depend on data, such as leap seconds and ΔT. How should that be represented?
- **Formal language design.** Could we design a small modelling language whose grammar and type system are developed together?
