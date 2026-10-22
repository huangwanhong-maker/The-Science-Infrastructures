# 0011. First Application: Celestial Atlas, with Provisional Computation in JavaScript

- **Status:** Accepted
- **Date:** 2026-10-22
- **Decided by:** the user asked for an application first: 3D views of planets and stars from real data, with 2D diagrams at a chosen time or in real time. Claude chose the form.

## 1. Context

The roadmap begins with Phase 0 (toolchain and build), followed by a vertical slice through the C++ core. The user asked for an application first. No C++ toolchain is installed on this machine, and the core does not exist yet.

## 2. Decision

1. **Form.** Build `applications/astronomy/celestial-atlas` as a static web application. It uses ES modules and three.js 0.170 from a CDN, and is served locally by a Python server that uses only the standard library.
2. **Data.** Data are registered in `data/registry/`. `data/acquisition/fetch.py` fetches them into `data/cache/` and records their checksums. Derived products are built by `data/products/`.
3. **Provisional computation.** The computation is written in JavaScript, in the application's `src/core/`. It covers time scales, reference frames, the DE440 evaluator, apparent places and osculating orbits. These are provisional stand-ins for:
   - `core/metrology`;
   - `core/domains/astronomy`;
   - physis.

   Each file names the module it anticipates.
4. **Testing.** The stand-ins are tested against JPL Horizons and against Meeus's worked examples. The same reference values will test the C++ core. The C++ core may later reach the browser as WebAssembly.
5. **Views.** The views are written directly with three.js and Canvas. The scene model for `visualization/` will be extracted from them ([architecture §1.6](../../docs/architecture.md)).

## 3. Alternatives Considered

- **Wait for the C++ core.** Rejected: the user asked for an application now, and Phase 0 has not begun.
- **Compute positions in a Python server.** Rejected: real-time animation needs the computation in the browser, and Python is not the core language either.
- **Precomputed position tables.** Rejected: they cannot serve arbitrary instants or real time.

## 4. Consequences

- **Duplicate code.** For a while, some computation exists twice: in JavaScript now and in C++ later. The shared tests keep the two consistent.
- **Download size.** The browser loads the 32 MB DE440 kernel. This is acceptable locally, but must be revisited for a public deployment, for example by serving a subset of years.
- **Data layer.** The data layer (registry, acquisition and products) now exists in working form.
