# core: The Laboratory's Scientific Library

## 1. Purpose

`core/` holds the laboratory's own pure computation: the parts specific to its subject. General numerics, mathematics and physics live in the libraries under `lib/` ([ADR 0008](../development/decisions/0008-independent-libraries.md)).

Like those libraries, `core/` is pure:
- It has no network access, no instrument I/O, no user interface and no hidden global state.
- It reads only the files whose paths the caller passes in.

Computational code is written in C++ ([ADR 0003](../development/decisions/0003-cpp-computational-core.md)).

## 2. Layers

| Order | Layer | Scope | Status |
|---|---|---|---|
| 0 | `lib/` (calx, morphe, physis) | Numerics and compute; structures; geometry and physics | planned |
| 1 | [`metrology/`](metrology/) | Astronomical time scales, celestial reference frames, Earth orientation | planned |
| 2 | [`domains/`](domains/) | Fields: astronomy, geography, oceanography, meteorology, geodesy. Bridges: tides, navigation, air–sea, climate | planned |
| 3 | [`inference/`](inference/) | Estimation, filtering, assimilation, uncertainty, learning | planned |

## 3. Dependency Rule

- A layer may depend on any layer with a lower order, including any library in `lib/`.
- Within a layer, one module may depend on another if the dependency is declared and the dependency graph stays acyclic.
- Nothing here depends on anything outside `core/` and `lib/`.

The full rules are in [architecture §3](../docs/architecture.md).
