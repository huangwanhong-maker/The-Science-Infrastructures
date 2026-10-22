# 0005. Visualization as a Representation Layer

- **Status:** Accepted
- **Date:** 2026-10-21
- **Decided by:** the user ("visualization can be various … a representation layer … remain extensible to support future techniques")

## 1. Context

Visualization techniques change quickly. They range from static figures to GPU rendering on the web, native engines, immersive displays, and techniques still to come. The project needs scientific figures, and also 3D views of the planets and the universe built from real data.

## 2. Decision

`visualization/` is a representation layer at the edge of the architecture. Only `applications/`, `notebooks/` and `research/` may use it. Nothing in `core/`, `data/` or `observation/` depends on it.

It has three parts ([architecture §5](../../docs/architecture.md)):

- **Scene model** (`scene/`): device-neutral.
- **Renderers** (`renderers/`): pluggable, one per technique.
- **Styles** (`styles/`): shared by all renderers.

Units, time scales and reference frames travel with the data into the scene, so renderers never have to guess them.

## 3. Alternatives Considered

- **Committing to one renderer**, for example three.js. Rejected: it would lock the project to one technology.
- **Rendering inside each application.** Rejected: it duplicates work and lets the visual style drift.

## 4. Consequences

- **The scene model is one more interface to design.** It starts minimal, in Phase 1.
- **A new technique means a new renderer**, with no change elsewhere.
- **Visual style is defined once**, in `styles/`.
