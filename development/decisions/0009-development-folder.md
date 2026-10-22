# 0009. A Development Folder for Process Material

- **Status:** Accepted
- **Date:** 2026-10-21
- **Decided by:** the user ("ideas should not be in the root folder, but can be in something like development/ideas")

## 1. Context

Ideas sat at the root of the repository, and the decision log and roadmap sat in `docs/`. Material about *how the project develops* was therefore mixed with the product and with its documentation.

## 2. Decision

- `development/` holds the process material:
  - `ideas/`;
  - `decisions/`;
  - `roadmap.md`.
- `docs/` holds documentation of the product. For now that is the architecture; later it will also hold the glossary and references.
- calx, morphe and physis follow the same convention.

## 3. Alternatives Considered

- **Keep ideas at the root.** Rejected: it clutters the root with process material.
- **Put everything in `docs/`.** Rejected: it mixes documentation of the product with records of its development.

## 4. Consequences

- Three paths change:
  - `ideas/` became `development/ideas/`;
  - `docs/decisions/` became `development/decisions/`;
  - `docs/roadmap.md` became `development/roadmap.md`.
- This record amends [0001](0001-record-architecture-decisions.md), which gave the location of the decision log.
