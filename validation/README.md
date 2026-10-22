# validation: Verification Against References

*Status: planned.* Assurance.

## 1. Purpose

This directory shows that the results are correct, not merely that the code runs.

- **Unit tests** live with each module.
- **Each library** verifies its own methods in its own repository.
- **This directory** validates the laboratory's results against external references, and holds checks that span several modules.

## 2. Planned Content

| Directory | Content |
|---|---|
| `references/` | Authoritative values from sources such as SOFA/ERFA, JPL Horizons and NOAA predictions. Each value records its source and retrieval date. |
| `problems/` | Standard test problems. Examples: Williamson's shallow-water tests on the sphere; the Held–Suarez and Jablonowski–Williamson tests for atmospheric dynamical cores; energy drift in n-body runs. |
| `benchmarks/` | Performance for each backend and each problem size |
| `suites/` | Verification suites that span several modules |
