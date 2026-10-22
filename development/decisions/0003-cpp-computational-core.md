# 0003. C++ for the Computational Core

- **Status:** Accepted
- **Date:** 2026-10-21
- **Decided by:** the user

## 1. Context

The core must be fast, must run on CPUs and other devices, and must last for many years. Four options were discussed:

- a C++ core with Python orchestration;
- Julia throughout;
- Rust;
- Python with JAX.

## 2. Decision

- Write the core in C++ wherever possible, especially numerical and computational code.
- Other parts may use Python, or any language suited to the task. These parts include data acquisition, orchestration, notebooks, machine learning and visualization.

## 3. Alternatives Considered

- **Julia throughout.** It offers a strong scientific ecosystem (for example Oceananigans.jl and KernelAbstractions.jl) and one language for everything. Not chosen: C++ reaches more devices, toolchains and HPC centres, and is the natural host for a compiler stack.
- **Rust.** Safe and modern, but its HPC and GPU ecosystem is younger.
- **Python with JAX.** Excellent for differentiable computing. Not chosen: it ties the core to one compiler (XLA) and one runtime.

## 4. Consequences

- **A C++ toolchain is needed.** As of 2026-10-21 the Windows side has none (roadmap Phase 0).
- **Python bindings are needed**, so that notebooks and the edges can call the core.
- **Provisional details**, to be confirmed in Phase 0:
  - C++20;
  - CMake;
  - nanobind for the bindings.
